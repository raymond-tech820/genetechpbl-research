"""Resumable GSE176206 download, h5ad adaptation, and external validation."""

from __future__ import annotations

import torch  # Import torch before h5py/anndata to avoid a Windows CUDA DLL load-order conflict.

import gzip
import h5py
import logging
import re
import shutil
import time
import zlib
from pathlib import Path
from typing import Any, Sequence

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import requests
import scipy.sparse as sp
from scipy.stats import spearmanr
from tqdm import tqdm
from anndata._core.sparse_dataset import sparse_dataset
from anndata._io.specs import read_elem

from .config import load_yaml
from .data import load_prepared_bundle
from .features import log_normalize_counts
from .scoring import YouthScoreEnsemble, align_named_matrix
from .utils import ensure_directory, read_json, sha256_file, utc_now, write_json


LOGGER = logging.getLogger(__name__)


def download_resumable(
    url: str,
    destination: Path,
    chunk_size: int = 8 * 1024 * 1024,
    max_retries: int = 100,
) -> Path:
    """Download a large file with byte-range resume and automatic reconnects."""

    if destination.exists() and destination.stat().st_size > 0:
        return destination
    ensure_directory(destination.parent)
    partial = destination.with_suffix(destination.suffix + ".part")
    attempt = 0
    expected_total: int | None = None
    while attempt <= max_retries:
        existing = partial.stat().st_size if partial.exists() else 0
        headers = {"Range": f"bytes={existing}-"} if existing else {}
        try:
            with requests.get(url, headers=headers, stream=True, timeout=(30, 300)) as response:
                if response.status_code == 416:
                    match = re.search(r"\*/(\d+)", response.headers.get("Content-Range", ""))
                    if match and existing == int(match.group(1)):
                        partial.replace(destination)
                        return destination
                    response.raise_for_status()
                if response.status_code not in {200, 206}:
                    response.raise_for_status()

                if response.status_code == 206:
                    match = re.fullmatch(
                        r"bytes (\d+)-(\d+)/(\d+|\*)",
                        response.headers.get("Content-Range", ""),
                    )
                    if not match or int(match.group(1)) != existing:
                        raise RuntimeError(
                            f"Invalid resume response for {destination.name}: "
                            f"{response.headers.get('Content-Range')!r}"
                        )
                    if match.group(3) != "*":
                        expected_total = int(match.group(3))
                    mode = "ab"
                else:
                    if existing:
                        raise requests.exceptions.ConnectionError(
                            f"Server ignored the Range request for {destination.name}; "
                            "preserving the partial file and retrying"
                        )
                    mode = "wb"
                    length = int(response.headers.get("Content-Length", 0))
                    expected_total = length or expected_total

                with partial.open(mode) as handle:
                    with tqdm(
                        total=expected_total,
                        initial=existing,
                        unit="B",
                        unit_scale=True,
                        desc=destination.name,
                    ) as progress:
                        for chunk in response.iter_content(chunk_size=chunk_size):
                            if chunk:
                                handle.write(chunk)
                                progress.update(len(chunk))

            downloaded = partial.stat().st_size
            if expected_total is None or downloaded == expected_total:
                partial.replace(destination)
                return destination
            if downloaded > expected_total:
                raise RuntimeError(
                    f"Downloaded file exceeds advertised size: {downloaded} > {expected_total}"
                )
            raise requests.exceptions.ChunkedEncodingError(
                f"Connection ended at {downloaded}/{expected_total} bytes"
            )
        except (requests.RequestException, OSError) as error:
            attempt += 1
            if attempt > max_retries:
                raise RuntimeError(
                    f"Unable to download {destination.name} after {max_retries} reconnects"
                ) from error
            delay = min(2 ** min(attempt, 5), 30)
            LOGGER.warning(
                "Download interrupted for %s at %d bytes; reconnect %d/%d in %ds: %s",
                destination.name,
                partial.stat().st_size if partial.exists() else 0,
                attempt,
                max_retries,
                delay,
                error,
            )
            time.sleep(delay)

    raise AssertionError("Download retry loop terminated unexpectedly")


def decompress_gzip(source: Path, destination: Path, chunk_size: int = 16 * 1024 * 1024) -> Path:
    if destination.exists() and destination.stat().st_size > 0:
        return destination
    partial = destination.with_suffix(destination.suffix + ".part")
    ensure_directory(destination.parent)
    try:
        with gzip.open(source, "rb") as compressed, partial.open("wb") as output:
            shutil.copyfileobj(compressed, output, length=chunk_size)
    except (OSError, EOFError, zlib.error):
        partial.unlink(missing_ok=True)
        raise
    partial.replace(destination)
    return destination


def download_gse176206(config_path: str | Path, keys: Sequence[str] | None = None) -> Path:
    config = load_yaml(config_path)
    root = Path(config_path).resolve().parents[1]
    external_dir = ensure_directory(root / "data" / "external" / config["accession"])
    selected = set(keys or config["files"].keys())
    manifest: dict[str, Any] = {"accession": config["accession"], "created_at": utc_now(), "files": {}}
    for key, item in config["files"].items():
        if key not in selected:
            continue
        compressed = external_dir / item["filename"]
        h5ad_path = compressed.with_suffix("")
        for integrity_attempt in range(2):
            if not compressed.exists():
                download_resumable(item["url"], compressed)
            try:
                decompress_gzip(compressed, h5ad_path)
                break
            except (OSError, EOFError, zlib.error) as error:
                if integrity_attempt == 1:
                    raise RuntimeError(
                        f"Downloaded {compressed.name} failed gzip integrity twice"
                    ) from error
                LOGGER.warning(
                    "Gzip integrity failed for %s; deleting the corrupt copy and downloading again: %s",
                    compressed.name,
                    error,
                )
                compressed.unlink(missing_ok=True)
                compressed.with_suffix(compressed.suffix + ".part").unlink(missing_ok=True)
        manifest["files"][key] = {
            "url": item["url"],
            "compressed": str(compressed),
            "compressed_bytes": compressed.stat().st_size,
            "compressed_sha256": sha256_file(compressed),
            "h5ad": str(h5ad_path),
            "h5ad_bytes": h5ad_path.stat().st_size,
            "h5ad_sha256": sha256_file(h5ad_path),
            "model_dataset_id": item["model_dataset_id"],
        }
    write_json(external_dir / "manifest.json", manifest)
    return external_dir


def _unique_strings(series: pd.Series, limit: int = 500) -> list[str]:
    values = series.dropna().astype(str).unique()
    return [str(value) for value in values[:limit]]


def _find_semantic_column(obs: pd.DataFrame, semantic: str) -> str | None:
    keywords = {
        "age": ("age",),
        "treatment": ("treatment", "condition", "stim", "perturb", "protocol"),
        "replicate": ("donor", "mouse", "animal", "replicate", "sample", "batch"),
    }[semantic]
    value_tokens = {
        "age": ("young", "aged", "old"),
        "treatment": ("sokm", "oskm", "control", "ctrl", "tg+", "tg-"),
        "replicate": (),
    }[semantic]
    candidates = []
    for column in obs.columns:
        column_name = str(column).lower()
        if semantic == "replicate" and any(
            token in column_name
            for token in ("treatment", "condition", "age", "cell_type", "celltype")
        ):
            continue
        name_score = sum(token in column_name for token in keywords)
        if semantic == "age" and column_name == "age":
            name_score += 10
        if semantic == "treatment" and column_name in {"treatment", "sample_treatment"}:
            name_score += 10
        if semantic == "replicate" and any(
            token in column_name for token in ("donor", "mouse", "animal")
        ):
            name_score += 10
        values = _unique_strings(obs[column])
        value_score = sum(any(token in value.lower() for token in value_tokens) for value in values)
        unique_count = len(values)
        if semantic == "replicate":
            valid = name_score > 0 and 2 <= unique_count <= 1_000
        else:
            valid = value_score > 0
        if valid:
            candidates.append((name_score, value_score, -unique_count, str(column)))
    return max(candidates)[-1] if candidates else None


def infer_external_schema(obs: pd.DataFrame) -> dict[str, str | None]:
    return {
        "age_column": _find_semantic_column(obs, "age"),
        "treatment_column": _find_semantic_column(obs, "treatment"),
        "replicate_column": _find_semantic_column(obs, "replicate"),
    }


def _map_age(value: Any) -> str:
    text = str(value).lower()
    if "young" in text:
        return "young"
    if "aged" in text or "old" in text:
        return "aged"
    return "unknown"


def _map_treatment(value: Any) -> str:
    text = str(value).lower()
    if (
        "control" in text
        or "ctrl" in text
        or "negctrl" in text
        or "dox-" in text
        or "tg-" in text
    ):
        return "control"
    if "sokm" in text or "oskm" in text or "tg+" in text:
        return "SOKM"
    return "unknown"


def _matrix_node(node):
    return sparse_dataset(node) if isinstance(node, h5py.Group) else node


def _open_external_h5ad(path: Path):
    """Open only obs, var, and one count matrix without materializing h5ad layers."""

    handle = h5py.File(path, "r")
    try:
        obs = read_elem(handle["obs"])
        for layer in ("counts", "raw_counts", "count"):
            node_path = f"layers/{layer}"
            if node_path in handle:
                return handle, _matrix_node(handle[node_path]), obs, read_elem(handle["var"]), node_path
        if "raw/X" in handle and "raw/var" in handle:
            return handle, _matrix_node(handle["raw/X"]), obs, read_elem(handle["raw/var"]), "raw/X"
        if "X" in handle:
            return handle, _matrix_node(handle["X"]), obs, read_elem(handle["var"]), "X"
        raise ValueError(f"No expression matrix was found in {path}")
    except Exception:
        handle.close()
        raise


def _resolve_gene_names(var: pd.DataFrame, var_names: Sequence[str], model_genes: set[str]) -> tuple[list[str], str]:
    candidates: list[tuple[str, Sequence[Any]]] = [("var_names", var_names)]
    for column in var.columns:
        if any(token in str(column).lower() for token in ("gene", "symbol", "feature")):
            candidates.append((str(column), var[column].astype(str).tolist()))
    scored = []
    for name, values in candidates:
        overlap = len(set(map(str, values)) & model_genes)
        scored.append((overlap, name, list(map(str, values))))
    overlap, name, values = max(scored, key=lambda item: item[0])
    if overlap == 0:
        return list(map(str, var_names)), "var_names"
    return values, name


def _as_csr(value) -> sp.csr_matrix:
    if sp.issparse(value):
        return value.tocsr()
    return sp.csr_matrix(np.asarray(value))


def _sample_backed_matrix(matrix, rows: int, maximum_cells: int = 4_096) -> sp.csr_matrix:
    if rows <= maximum_cells:
        return _as_csr(matrix[:rows])
    block_size = maximum_cells // 4
    starts = np.linspace(0, rows - block_size, num=4, dtype=int)
    return sp.vstack(
        [_as_csr(matrix[start : start + block_size]) for start in starts],
        format="csr",
    )


def _identity_evidence(obs: pd.DataFrame, h5ad_path: Path) -> dict[str, Any]:
    identity_columns = [
        str(column)
        for column in obs.columns
        if any(
            token in str(column).lower()
            for token in ("cell_type", "celltype", "annotation", "identity", "tissue")
        )
    ]
    observed = {
        column: _unique_strings(obs[column], limit=50)
        for column in identity_columns
    }
    filename = h5ad_path.name.lower()
    expected_terms = ("adipo", "adipogenic") if "adipo" in filename else ("msc", "mesenchymal")
    metadata_text = " ".join(
        str(value).lower()
        for values in observed.values()
        for value in values
    )
    if any(term in metadata_text for term in expected_terms):
        status = "metadata_match"
    elif any(term in filename for term in expected_terms):
        status = "source_file_declared"
    else:
        status = "unresolved"
    return {
        "status": status,
        "expected_terms": list(expected_terms),
        "identity_columns": identity_columns,
        "observed_identity_values": observed,
    }


def _external_preflight(
    matrix,
    gene_names: Sequence[str],
    obs: pd.DataFrame,
    h5ad_path: Path,
    scorer: YouthScoreEnsemble,
) -> dict[str, Any]:
    input_genes = set(map(str, gene_names))
    fold_overlaps = {
        str(fold.fold): len(input_genes & set(fold.state.gene_names)) / len(fold.state.gene_names)
        for fold in scorer.folds
    }
    model_genes = sorted({gene for fold in scorer.folds for gene in fold.state.gene_names})
    project_root = Path(scorer.selection["config"]["dataset"]["project_root"])
    dataset_id = str(scorer.selection["dataset_id"])
    reference = load_prepared_bundle(
        project_root / "data" / "processed" / "youth_score" / dataset_id
    )
    reference_gene_names = reference.genes["gene_name"].astype(str).tolist()
    shared = sorted(input_genes & set(reference_gene_names) & set(model_genes))
    external_sample = _sample_backed_matrix(matrix, len(obs))
    correlation = float("nan")
    if len(shared) >= 3:
        reference_aligned, _ = align_named_matrix(reference.counts, reference_gene_names, shared)
        external_aligned, _ = align_named_matrix(external_sample, gene_names, shared)
        reference_profile = np.asarray(log_normalize_counts(reference_aligned).mean(axis=0)).ravel()
        external_profile = np.asarray(log_normalize_counts(external_aligned).mean(axis=0)).ravel()
        valid = np.isfinite(reference_profile) & np.isfinite(external_profile)
        if valid.sum() >= 3:
            correlation = float(spearmanr(reference_profile[valid], external_profile[valid]).statistic)
    return {
        "fold_gene_overlap": fold_overlaps,
        "minimum_gene_overlap": float(min(fold_overlaps.values())),
        "shared_reference_genes": len(shared),
        "reference_expression_spearman": correlation,
        "reference_expression_status": "positive"
        if np.isfinite(correlation) and correlation > 0
        else "warning_nonpositive_or_unresolved",
        "reference_dataset_id": dataset_id,
        "sampled_external_cells": int(external_sample.shape[0]),
        "identity": _identity_evidence(obs, h5ad_path),
    }


def _external_summary(scores_path: Path) -> dict[str, Any]:
    scores = pd.read_parquet(
        scores_path,
        columns=["external_age_group", "external_treatment", "external_replicate", "youth_score"],
    )
    valid = scores[
        scores["external_age_group"].isin(["young", "aged"])
        & scores["external_treatment"].isin(["control", "SOKM"])
        & ~scores["external_replicate"].astype(str).str.lower().isin(
            ["", "unknown", "nan", "none", "negative", "unassigned"]
        )
    ].copy()
    replicate = (
        valid.groupby(["external_age_group", "external_treatment", "external_replicate"], as_index=False)
        .agg(youth_score=("youth_score", "mean"), cells=("youth_score", "size"))
    )
    condition = (
        replicate.groupby(["external_age_group", "external_treatment"], as_index=False)
        .agg(youth_score=("youth_score", "mean"), replicates=("external_replicate", "nunique"), cells=("cells", "sum"))
    )
    rng = np.random.default_rng(20260717)

    def bootstrap_mean_interval(
        values: np.ndarray,
        draws: int = 2_000,
    ) -> tuple[float | None, float | None]:
        values = np.asarray(values, dtype=float)
        if len(values) < 2:
            return None, None
        means = rng.choice(values, size=(draws, len(values)), replace=True).mean(axis=1)
        low, high = np.quantile(means, [0.025, 0.975])
        return float(low), float(high)

    interval_rows = []
    for _, row in condition.iterrows():
        values = replicate[
            replicate["external_age_group"].eq(row["external_age_group"])
            & replicate["external_treatment"].eq(row["external_treatment"])
        ]["youth_score"].to_numpy(dtype=float)
        low, high = bootstrap_mean_interval(values)
        interval_rows.append((low, high))
    condition["ci_low"] = [value[0] for value in interval_rows]
    condition["ci_high"] = [value[1] for value in interval_rows]

    def mean_for(age: str, treatment: str) -> float | None:
        value = condition[(condition["external_age_group"] == age) & (condition["external_treatment"] == treatment)]
        return None if value.empty else float(value["youth_score"].iloc[0])

    young_control = mean_for("young", "control")
    aged_control = mean_for("aged", "control")
    aged_sokm = mean_for("aged", "SOKM")
    contrasts = {
        "young_control_minus_aged_control": None
        if young_control is None or aged_control is None
        else young_control - aged_control,
        "aged_sokm_minus_aged_control": None
        if aged_sokm is None or aged_control is None
        else aged_sokm - aged_control,
    }

    def values_for(age: str, treatment: str) -> np.ndarray:
        return replicate[
            replicate["external_age_group"].eq(age)
            & replicate["external_treatment"].eq(treatment)
        ]["youth_score"].to_numpy(dtype=float)

    def bootstrap_difference_interval(
        first: np.ndarray,
        second: np.ndarray,
        draws: int = 2_000,
    ) -> tuple[float, float] | None:
        if len(first) < 2 or len(second) < 2:
            return None
        first_means = rng.choice(first, size=(draws, len(first)), replace=True).mean(axis=1)
        second_means = rng.choice(second, size=(draws, len(second)), replace=True).mean(axis=1)
        low, high = np.quantile(first_means - second_means, [0.025, 0.975])
        return float(low), float(high)

    contrast_intervals = {
        "young_control_minus_aged_control": bootstrap_difference_interval(
            values_for("young", "control"), values_for("aged", "control")
        ),
        "aged_sokm_minus_aged_control": bootstrap_difference_interval(
            values_for("aged", "SOKM"), values_for("aged", "control")
        ),
    }
    replicate.to_csv(scores_path.with_name("replicate_scores.csv"), index=False)
    condition.to_csv(scores_path.with_name("condition_scores.csv"), index=False)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    plot_data = replicate.copy()
    plot_data["condition"] = (
        plot_data["external_age_group"].str.title()
        + " / "
        + plot_data["external_treatment"]
    )
    order = ["Young / control", "Aged / control", "Aged / SOKM", "Young / SOKM"]
    order = [value for value in order if value in set(plot_data["condition"])]
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(8, 5))
    sns.boxplot(data=plot_data, x="condition", y="youth_score", order=order, color="#9ecae1")
    sns.stripplot(data=plot_data, x="condition", y="youth_score", order=order, color="#08306b", size=7)
    plt.xlabel("External condition")
    plt.ylabel("Donor/replicate Youth Score")
    plt.title("GSE176206 external validation")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    plt.savefig(scores_path.with_name("condition_scores.png"), dpi=180, bbox_inches="tight")
    plt.close()
    minimum_replicates = int(condition["replicates"].min()) if not condition.empty else 0
    condition_records = (
        condition.astype(object)
        .where(pd.notna(condition), None)
        .to_dict(orient="records")
    )
    return {
        "conditions": condition_records,
        "contrasts": contrasts,
        "contrast_intervals_95": contrast_intervals,
        "replication_status": "replicated" if minimum_replicates >= 2 else "replication_limited",
        "minimum_replicates_per_condition": minimum_replicates,
        "young_control_sanity_pass": contrasts["young_control_minus_aged_control"] is not None
        and contrasts["young_control_minus_aged_control"] > 0,
        "aged_sokm_rejuvenation_direction": contrasts["aged_sokm_minus_aged_control"] is not None
        and contrasts["aged_sokm_minus_aged_control"] > 0,
    }


def validate_external_file(
    h5ad_path: str | Path,
    model_output_directory: str | Path,
    output_directory: str | Path,
    minimum_gene_overlap: float = 0.70,
    chunk_size: int = 2_048,
    device: str = "auto",
    force: bool = False,
) -> Path:
    h5ad_path = Path(h5ad_path)
    output_directory = ensure_directory(output_directory)
    scores_path = output_directory / "cell_scores.parquet"
    if (
        not force
        and scores_path.exists()
        and (output_directory / "validation_summary.json").exists()
    ):
        return output_directory
    scorer = YouthScoreEnsemble(model_output_directory, device=device)
    model_genes = {gene for fold in scorer.folds for gene in fold.state.gene_names}
    h5ad_handle, matrix, obs, var, matrix_source = _open_external_h5ad(h5ad_path)
    var_names = list(map(str, var.index))
    gene_names, gene_name_source = _resolve_gene_names(var, var_names, model_genes)
    schema = infer_external_schema(obs)
    write_json(output_directory / "metadata_schema.json", {**schema, "matrix_source": matrix_source, "gene_name_source": gene_name_source})
    preflight = _external_preflight(matrix, gene_names, obs, h5ad_path, scorer)
    write_json(output_directory / "preflight.json", preflight)
    if preflight["minimum_gene_overlap"] < minimum_gene_overlap:
        raise ValueError(
            f"Model gene overlap {preflight['minimum_gene_overlap']:.3f} is below the "
            f"required {minimum_gene_overlap:.3f}."
        )

    missing_semantics = [
        name for name, value in schema.items() if value is None
    ]
    if missing_semantics:
        raise ValueError(
            "External validation requires explicit age, treatment, and biological-replicate metadata; "
            f"missing {missing_semantics}. Available columns: {list(map(str, obs.columns))}"
        )

    age_source = obs[str(schema["age_column"])].map(_map_age)
    treatment_source = obs[str(schema["treatment_column"])].map(_map_treatment)
    replicate_source = obs[str(schema["replicate_column"])].astype(str)
    if float(age_source.eq("unknown").mean()) > 0.05:
        raise ValueError("More than 5% of external cells have an unrecognized age label.")
    if float(treatment_source.eq("unknown").mean()) > 0.05:
        raise ValueError("More than 5% of external cells have an unrecognized treatment label.")
    observed_conditions = set(zip(age_source, treatment_source, strict=True))
    expected_conditions = {
        ("young", "control"),
        ("young", "SOKM"),
        ("aged", "control"),
        ("aged", "SOKM"),
    }
    if not expected_conditions.issubset(observed_conditions):
        raise ValueError(
            f"External file is missing required conditions: {sorted(expected_conditions - observed_conditions)}"
        )
    writer: pq.ParquetWriter | None = None
    minimum_observed_overlap = 1.0
    try:
        for start in tqdm(range(0, matrix.shape[0], chunk_size), desc=f"Scoring {h5ad_path.name}"):
            stop = min(start + chunk_size, matrix.shape[0])
            chunk = _as_csr(matrix[start:stop])
            cell_ids = obs.index[start:stop].astype(str).tolist()
            scored = scorer.score_counts(chunk, gene_names, cell_ids)
            minimum_observed_overlap = min(minimum_observed_overlap, float(scored["gene_overlap"].min()))
            if minimum_observed_overlap < minimum_gene_overlap:
                raise ValueError(
                    f"Model gene overlap {minimum_observed_overlap:.3f} is below the required {minimum_gene_overlap:.3f}."
                )
            scored["external_age_group"] = age_source.iloc[start:stop].to_numpy()
            scored["external_treatment"] = treatment_source.iloc[start:stop].to_numpy()
            scored["external_replicate"] = replicate_source.iloc[start:stop].astype(str).to_numpy()
            table = pa.Table.from_pandas(scored, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(scores_path, table.schema, compression="zstd")
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()
        h5ad_handle.close()
    summary = _external_summary(scores_path)
    report_data = {
        "created_at": utc_now(),
        "source_h5ad": str(h5ad_path),
        "model_output_directory": str(model_output_directory),
        "official_model": scorer.official_model,
        "minimum_gene_overlap": minimum_observed_overlap,
        "preflight": preflight,
        "metadata_schema": schema,
        **summary,
    }
    write_json(output_directory / "validation_summary.json", report_data)
    condition_table = pd.DataFrame(summary["conditions"]).round(4).to_markdown(index=False)
    report = f"""# GSE176206 External Validation Report

Generated: {report_data['created_at']}

- Source: `{h5ad_path}`
- Model: `{scorer.official_model}` from `{model_output_directory}`
- Minimum model-gene overlap: `{minimum_observed_overlap:.3f}`
- TMS/external reference-expression Spearman: `{preflight['reference_expression_spearman']:.3f}`
- Cell-identity evidence: `{preflight['identity']['status']}`
- External replication status: `{summary['replication_status']}`
- Young-control sanity direction: `{summary['young_control_sanity_pass']}`
- Aged-SOKM rejuvenation direction: `{summary['aged_sokm_rejuvenation_direction']}`

## Donor/replicate-aggregated conditions

{condition_table}

The external files were not used for feature selection, fitting, early stopping, or calibration. Confidence intervals use 2,000 biological-replicate bootstrap draws only when at least two identifiable replicates are available in each contrasted condition. A `replication_limited` result is descriptive because cell-level variation is not a substitute for biological replication. These results measure expression similarity to the TMS young reference and do not establish safety, preserved cell identity, causality, or treatment efficacy.
"""
    (output_directory / "validation_report.md").write_text(report, encoding="utf-8")
    return output_directory


def validate_gse176206(
    config_path: str | Path,
    device: str = "auto",
    force: bool = False,
) -> Path:
    config_path = Path(config_path).resolve()
    config = load_yaml(config_path)
    root = config_path.parents[1]
    external_dir = root / "data" / "external" / config["accession"]
    output_root = ensure_directory(root / "outputs" / "youth_score" / "external" / config["accession"])
    summaries = {}
    for key, item in config["files"].items():
        compressed = external_dir / item["filename"]
        h5ad_path = compressed.with_suffix("")
        if not h5ad_path.exists():
            raise FileNotFoundError(f"External h5ad file is missing: {h5ad_path}")
        model_output = root / "outputs" / "youth_score" / item["model_dataset_id"]
        target = validate_external_file(
            h5ad_path,
            model_output,
            output_root / key,
            minimum_gene_overlap=float(config["minimum_gene_overlap"]),
            chunk_size=int(config["chunk_size"]),
            device=device,
            force=force,
        )
        summaries[key] = read_json(target / "validation_summary.json")
    write_json(output_root / "combined_summary.json", summaries)
    return output_root

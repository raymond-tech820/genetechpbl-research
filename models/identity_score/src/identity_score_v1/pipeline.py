"""End-to-end Identity Score v1 pipeline."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tomllib

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from .aggregation import aggregate_identity_scores, write_aggregated_outputs
from .h5ad_io import open_h5ad
from .metadata import normalize_metadata
from .resources import load_frozen_resources
from .scoring import (
    PRIMARY_COLUMN,
    RANK_COLUMN,
    log_normalize_counts,
    prepare_identity_plan,
    score_chunk,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _resolve(base: Path, value: str | Path) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def _load_config(path: str | Path) -> tuple[dict, Path]:
    config_path = Path(path).resolve()
    with config_path.open("rb") as handle:
        config = tomllib.load(handle)
    project_root = _resolve(
        config_path.parent, config.get("project_root", ".")
    )
    return config, project_root


def _reference_log_means(
    backed,
    reference_mask: np.ndarray,
    chunk_size: int,
    target_sum: float,
) -> tuple[np.ndarray, int]:
    totals = np.zeros(backed.shape[1], dtype=np.float64)
    cells = 0
    for start, stop, counts in backed.iter_chunks(chunk_size):
        local = reference_mask[start:stop]
        if not local.any():
            continue
        normalized = log_normalize_counts(counts[local], target_sum=target_sum)
        totals += np.asarray(normalized.sum(axis=0)).ravel()
        cells += int(local.sum())
    if cells == 0:
        raise ValueError("No known control cells were available for background matching.")
    return totals / cells, cells


def _score_cells(
    input_data,
    metadata: pd.DataFrame,
    plan,
    output_path: Path,
    chunk_size: int,
    target_sum: float,
) -> None:
    temporary = output_path.with_suffix(".parquet.part")
    temporary.unlink(missing_ok=True)
    writer: pq.ParquetWriter | None = None
    try:
        for start, stop, counts in input_data.counts.iter_chunks(chunk_size):
            cell_ids = metadata.iloc[start:stop]["cell_id"].tolist()
            scores = score_chunk(
                counts,
                plan,
                cell_ids=cell_ids,
                target_sum=target_sum,
            )
            frame = metadata.iloc[start:stop].reset_index(drop=True).merge(
                scores,
                on="cell_id",
                validate="one_to_one",
            )
            table = pa.Table.from_pandas(frame, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(
                    temporary,
                    table.schema,
                    compression="zstd",
                )
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()
    if writer is None:
        raise RuntimeError("No cell scores were written.")
    temporary.replace(output_path)


def run_pipeline(
    config_path: str | Path,
    *,
    force: bool = False,
    input_h5ad: str | Path | None = None,
    output_directory: str | Path | None = None,
) -> Path:
    config, project_root = _load_config(config_path)
    package = config["package"]
    input_config = config["input"]
    scoring = config["scoring"]
    metadata_config = config["metadata"]
    validation = config.get("validation", {})

    h5ad_path = (
        Path(input_h5ad).resolve()
        if input_h5ad is not None
        else _resolve(project_root, input_config["h5ad"])
    )
    output = (
        Path(output_directory).resolve()
        if output_directory is not None
        else _resolve(project_root, config["output"]["directory"])
    )
    output.mkdir(parents=True, exist_ok=True)
    raw_cell_path = output / "identity_cell_scores_raw.parquet"

    obs_columns = list(
        dict.fromkeys(
            [
                metadata_config.get("age_column", "age"),
                metadata_config.get("treatment_column", "treatment"),
                metadata_config.get("animal_column", "animal"),
                metadata_config.get("state_column", "state"),
                metadata_config.get("sample_column", "sample"),
                metadata_config.get("batch_column", "batch"),
            ]
        )
    )
    identity_genes, all_module_genes = load_frozen_resources()
    chunk_size = int(scoring.get("chunk_size", 1024))
    target_sum = float(scoring.get("target_sum", 10_000.0))

    with open_h5ad(
        h5ad_path,
        counts_layer=input_config.get("counts_layer", "counts"),
        obs_columns=obs_columns,
        gene_name_column=input_config.get("gene_name_column", "gene_name"),
    ) as input_data:
        metadata = normalize_metadata(input_data.obs, metadata_config)
        reference_mask = (
            metadata["known_animal"]
            & metadata["analysis_role"].eq("control")
        ).to_numpy()
        reference_means, reference_cells = _reference_log_means(
            input_data.counts,
            reference_mask,
            chunk_size=chunk_size,
            target_sum=target_sum,
        )
        plan, qc = prepare_identity_plan(
            input_data.gene_names,
            identity_genes,
            all_module_genes,
            reference_means,
            minimum_gene_coverage=float(
                scoring.get("minimum_gene_coverage", 0.70)
            ),
            n_control_sets=int(scoring.get("n_control_sets", 100)),
            control_bins=int(scoring.get("control_bins", 24)),
            seed=int(scoring.get("seed", 20260729)),
        )
        if plan.status != "pass":
            raise ValueError("Identity Score failed the frozen gene-coverage threshold.")
        if force or not raw_cell_path.exists():
            _score_cells(
                input_data,
                metadata,
                plan,
                raw_cell_path,
                chunk_size=chunk_size,
                target_sum=target_sum,
            )
        matrix_source = input_data.matrix_source
        gene_name_source = input_data.gene_name_source
        input_shape = input_data.counts.shape

    cells = pd.read_parquet(raw_cell_path)
    expected_cells = validation.get("expected_cells")
    if expected_cells is not None and len(cells) != int(expected_cells):
        raise AssertionError(
            f"Expected {expected_cells} cells, observed {len(cells)}."
        )
    tables = aggregate_identity_scores(
        cells,
        dataset_id=package["dataset_id"],
        species=package["species"],
        cell_population=package["cell_population"],
        unit_prefix=package["unit_prefix"],
        gene_coverage=float(qc.iloc[0]["gene_coverage"]),
        expected_biological_units=(
            int(validation["expected_biological_units"])
            if "expected_biological_units" in validation
            else None
        ),
    )
    write_aggregated_outputs(tables, output)
    qc.to_csv(output / "identity_gene_set_qc.csv", index=False)

    run_manifest = {
        "package": "Identity_Score_v1",
        "package_version": "1.0.0",
        "frozen_date": "2026-08-06",
        "config": str(Path(config_path).resolve()),
        "input_h5ad": str(h5ad_path),
        "input_bytes": h5ad_path.stat().st_size,
        "input_shape": list(input_shape),
        "matrix_source": matrix_source,
        "gene_name_source": gene_name_source,
        "known_control_reference_cells": reference_cells,
        "normalization": f"total-count normalize to {target_sum:g}, then log1p",
        "primary_method": (
            "mean identity-gene expression minus 100 expression-matched "
            "background module scores"
        ),
        "sensitivity_method": "within-cell expression-rank score",
        "primary_cell_field": PRIMARY_COLUMN,
        "sensitivity_cell_field": RANK_COLUMN,
        "known_biological_units": len(tables["donors"]),
        "cross_arm_pairing_status": "unpaired",
    }
    (output / "run_manifest.json").write_text(
        json.dumps(run_manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    validation_summary = {
        "status": "pass",
        "all_cell_scores": len(cells),
        "known_animal_cells": int(cells["known_animal"].sum()),
        "unknown_animal_cells": int((~cells["known_animal"]).sum()),
        "known_biological_units": len(tables["donors"]),
        "gene_coverage": float(qc.iloc[0]["gene_coverage"]),
        "missing_identity_scores": int(
            cells[[PRIMARY_COLUMN, RANK_COLUMN]].isna().sum().sum()
        ),
        "cross_arm_pairing_status": "unpaired",
    }
    (output / "validation_summary.json").write_text(
        json.dumps(validation_summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return output

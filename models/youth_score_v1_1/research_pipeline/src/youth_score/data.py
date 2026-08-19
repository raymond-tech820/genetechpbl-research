"""TMS extraction, validation, and training-bundle preparation."""

from __future__ import annotations

import re
import shutil
import sys
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
import scipy.sparse as sp
import zstandard as zstd
from sklearn.model_selection import StratifiedKFold

from .config import DatasetConfig
from .constants import ARCHIVE_NAMES, EXPECTED_MATRIX_SHAPES, MATRIX_DIR_NAMES, OBS_NAMES
from .utils import ensure_directory, runtime_info, sha256_file, utc_now, write_json


KNOWN_Y_GENES = {
    "Ddx3y",
    "Eif2s3y",
    "Gm29650",
    "Kdm5d",
    "Rbmy",
    "Rbmy1a1",
    "Sry",
    "Tmsb4y",
    "Usp9y",
    "Uty",
    "Zfy1",
    "Zfy2",
}


@dataclass(frozen=True)
class PreparedBundle:
    counts: sp.csr_matrix
    cells: pd.DataFrame
    genes: pd.DataFrame
    folds: pd.DataFrame
    directory: Path


def _safe_members(tar: tarfile.TarFile, destination: Path) -> Iterator[tarfile.TarInfo]:
    destination = destination.resolve()
    for member in tar:
        if member.issym() or member.islnk():
            raise ValueError(f"Archive links are not allowed: {member.name}")
        target = (destination / member.name).resolve()
        if destination != target and destination not in target.parents:
            raise ValueError(f"Unsafe archive member path: {member.name}")
        yield member


def extract_bpcells_archive(archive: Path, destination: Path, matrix_dir_name: str) -> Path:
    """Stream a tar.zst archive to disk with path traversal protection."""

    matrix_dir = destination / matrix_dir_name
    required = {"shape", "row_names", "col_names", "version", "val_data", "index_data"}
    if matrix_dir.is_dir() and required.issubset({p.name for p in matrix_dir.iterdir()}):
        return matrix_dir

    ensure_directory(destination)
    temporary = destination / f".{matrix_dir_name}.extracting"
    if temporary.exists():
        shutil.rmtree(temporary)
    ensure_directory(temporary)

    with archive.open("rb") as compressed:
        with zstd.ZstdDecompressor().stream_reader(compressed) as stream:
            with tarfile.open(fileobj=stream, mode="r|") as tar:
                tar.extractall(path=temporary, members=_safe_members(tar, temporary))

    extracted = temporary / matrix_dir_name
    if not extracted.is_dir():
        raise RuntimeError(f"Expected matrix directory was not found after extraction: {extracted}")
    if matrix_dir.exists():
        shutil.rmtree(matrix_dir)
    shutil.move(str(extracted), str(matrix_dir))
    shutil.rmtree(temporary, ignore_errors=True)
    return matrix_dir


def read_name_file(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8") as handle:
        return [line.rstrip("\r\n") for line in handle]


def open_bpcells_matrix(path: Path):
    try:
        from bpcells.experimental import DirMatrix
    except ImportError:
        from bpcells import DirMatrix  # type: ignore[attr-defined,no-redef]
    return DirMatrix(str(path))


def slice_bpcells_columns(matrix, positions: np.ndarray, chunk_size: int = 1_024) -> sp.csr_matrix:
    """Read selected cells in conservative chunks.

    The experimental BPCells Windows binding can access-violate when a large,
    multithreaded random-column request is issued. Single-threaded chunking is
    stable and still keeps peak memory small.
    """

    if hasattr(matrix, "threads"):
        matrix.threads = 1 if sys.platform == "win32" else min(8, max(1, len(positions) // chunk_size))
    chunks = []
    for start in range(0, len(positions), chunk_size):
        chunk_positions = positions[start : start + chunk_size]
        chunks.append(matrix[:, chunk_positions].T.tocsr().astype(np.uint32))
    if not chunks:
        return sp.csr_matrix((0, int(matrix.shape[0])), dtype=np.uint32)
    return sp.vstack(chunks, format="csr")


def validate_matrix(matrix_dir: Path, modality: str, metadata: pd.DataFrame) -> tuple[list[str], list[str]]:
    genes = read_name_file(matrix_dir / "row_names")
    cells = read_name_file(matrix_dir / "col_names")
    expected_shape = EXPECTED_MATRIX_SHAPES[modality]
    if (len(genes), len(cells)) != expected_shape:
        raise ValueError(
            f"Unexpected {modality} matrix names shape: {(len(genes), len(cells))}; "
            f"expected {expected_shape}."
        )
    if metadata.shape[0] != expected_shape[1]:
        raise ValueError(
            f"Unexpected {modality} metadata rows: {metadata.shape[0]}; expected {expected_shape[1]}."
        )
    metadata_ids = metadata["index"].astype(str).tolist()
    if cells != metadata_ids:
        first_mismatch = next((i for i, (a, b) in enumerate(zip(cells, metadata_ids)) if a != b), None)
        raise ValueError(f"Matrix and metadata cell IDs are not aligned; first mismatch={first_mismatch}.")
    return genes, cells


def parse_age_months(age: str) -> float:
    match = re.fullmatch(r"(\d+(?:\.\d+)?)m", str(age).strip())
    if not match:
        raise ValueError(f"Unsupported TMS age label: {age}")
    return float(match.group(1))


def gene_exclusion_reason(gene: str) -> str:
    lower = gene.lower()
    if lower.startswith("mt-"):
        return "mitochondrial"
    if re.match(r"^(Rpl|Rps)\d", gene):
        return "ribosomal"
    if re.match(r"^Hb[ab]-", gene) or re.match(r"^Hb[ab]\d", gene):
        return "hemoglobin"
    if gene in KNOWN_Y_GENES or lower.endswith("y") and gene in KNOWN_Y_GENES:
        return "y_chromosome"
    if gene in {"Xist", "Tsix"}:
        return "sex_associated"
    return ""


def build_gene_table(genes: list[str]) -> pd.DataFrame:
    table = pd.DataFrame({"gene_index": np.arange(len(genes), dtype=np.int32), "gene_name": genes})
    table["exclusion_reason"] = table["gene_name"].map(gene_exclusion_reason)
    table["eligible_primary"] = table["exclusion_reason"].eq("")
    return table


def assign_donor_folds(cells: pd.DataFrame, n_splits: int, seed: int) -> pd.DataFrame:
    labeled = cells.loc[cells["label"].notna(), ["mouse_id", "label", "age_group"]].drop_duplicates()
    if labeled["mouse_id"].duplicated().any():
        raise ValueError("Each mouse must have exactly one age label.")
    class_counts = labeled.groupby("label")["mouse_id"].nunique()
    if len(class_counts) != 2 or class_counts.min() < n_splits:
        raise ValueError(f"At least {n_splits} donors per class are required; observed {class_counts.to_dict()}.")
    splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    donor_table = labeled.sort_values("mouse_id").reset_index(drop=True)
    donor_table["outer_fold"] = -1
    dummy = np.zeros((len(donor_table), 1), dtype=np.float32)
    for fold, (_, test_idx) in enumerate(splitter.split(dummy, donor_table["label"].astype(int))):
        donor_table.loc[test_idx, "outer_fold"] = fold
    if (donor_table["outer_fold"] < 0).any():
        raise RuntimeError("Incomplete donor fold assignment.")
    fold_classes = donor_table.groupby("outer_fold")["label"].nunique()
    if not fold_classes.eq(2).all():
        raise RuntimeError("Every outer fold must contain young and old donors.")
    return donor_table.sort_values(["outer_fold", "label", "mouse_id"]).reset_index(drop=True)


def _selection_mask(metadata: pd.DataFrame, config: DatasetConfig) -> pd.Series:
    base = metadata["tissue"].eq(config.tissue) & metadata["cell_ontology_class"].eq(config.cell_type)
    if config.excluded_subtissues:
        if "subtissue" not in metadata:
            raise ValueError("The requested subtissue exclusion cannot be applied because metadata lacks 'subtissue'.")
        base &= ~metadata["subtissue"].fillna("").isin(config.excluded_subtissues)
    if config.role == "sensitivity":
        return base & metadata["age"].isin(config.sensitivity_ages)
    return base & metadata["age"].isin((*config.young_ages, *config.old_ages))


def _validate_expected_selection(selected: pd.DataFrame, config: DatasetConfig) -> None:
    young = selected["age"].isin(config.young_ages)
    old = selected["age"].isin(config.old_ages)
    observed = {
        "cells": len(selected),
        "young_cells": int(young.sum()),
        "old_cells": int(old.sum()),
        "young_mice": int(selected.loc[young, "mouse.id"].nunique()),
        "old_mice": int(selected.loc[old, "mouse.id"].nunique()),
    }
    expected = {
        "cells": config.expected_cells,
        "young_cells": config.expected_young_cells,
        "old_cells": config.expected_old_cells,
        "young_mice": config.expected_young_mice,
        "old_mice": config.expected_old_mice,
    }
    if observed != expected:
        raise ValueError(f"Dataset selection does not match the declared counts. observed={observed}, expected={expected}")


def prepare_tms_bundle(config: DatasetConfig, force: bool = False) -> Path:
    """Create a compact, validated sparse training bundle for one TMS target."""

    output_dir = config.processed_dir
    manifest_path = output_dir / "manifest.json"
    if manifest_path.exists() and not force:
        return output_dir

    source_dir = config.project_root / "data" / "tms_compact"
    archive = source_dir / ARCHIVE_NAMES[config.modality]
    obs_path = source_dir / OBS_NAMES[config.modality]
    if not archive.exists() or not obs_path.exists():
        raise FileNotFoundError(f"Missing TMS source files for {config.modality} in {source_dir}")

    extracted_root = config.project_root / "data" / "extracted"
    matrix_dir = extract_bpcells_archive(archive, extracted_root, MATRIX_DIR_NAMES[config.modality])
    metadata = pd.read_csv(obs_path, low_memory=False)
    genes, matrix_cells = validate_matrix(matrix_dir, config.modality, metadata)

    selection = _selection_mask(metadata, config)
    selected = metadata.loc[selection].copy()
    _validate_expected_selection(selected, config)
    selected_positions = np.flatnonzero(selection.to_numpy())

    matrix = open_bpcells_matrix(matrix_dir)
    if tuple(int(v) for v in matrix.shape) != EXPECTED_MATRIX_SHAPES[config.modality]:
        raise ValueError(f"BPCells matrix reports an unexpected shape: {matrix.shape}")
    counts = slice_bpcells_columns(matrix, selected_positions)
    counts.eliminate_zeros()

    total_counts = np.asarray(counts.sum(axis=1)).ravel().astype(np.int64)
    detected_genes = np.diff(counts.indptr).astype(np.int32)
    qc_pass = (total_counts > 0) & (detected_genes >= config.minimum_detected_genes)
    counts = counts[qc_pass].tocsr()
    selected = selected.iloc[np.flatnonzero(qc_pass)].copy().reset_index(drop=True)
    total_counts = total_counts[qc_pass]
    detected_genes = detected_genes[qc_pass]

    selected["cell_id"] = selected["index"].astype(str)
    selected["mouse_id"] = selected["mouse.id"].astype(str)
    selected["dataset_id"] = config.dataset_id
    selected["modality"] = config.modality
    selected["age_months"] = selected["age"].map(parse_age_months).astype(np.float32)
    selected["age_group"] = np.select(
        [selected["age"].isin(config.young_ages), selected["age"].isin(config.old_ages)],
        ["young", "old"],
        default="sensitivity",
    )
    selected["label"] = selected["age_group"].map({"young": 1.0, "old": 0.0})
    selected["total_counts_matrix"] = total_counts
    selected["detected_genes_matrix"] = detected_genes
    selected["qc_pass"] = True
    selected["source_matrix_column"] = selected_positions[qc_pass]
    selected["sex"] = selected["sex"].fillna("unknown").astype(str)

    genes_table = build_gene_table(genes)
    if config.role == "sensitivity":
        folds = pd.DataFrame(columns=["mouse_id", "label", "age_group", "outer_fold"])
        selected["outer_fold"] = -1
    else:
        folds = assign_donor_folds(selected, config.folds, config.seeds[0])
        selected = selected.merge(folds[["mouse_id", "outer_fold"]], on="mouse_id", how="left", validate="many_to_one")
        if selected["outer_fold"].isna().any():
            raise RuntimeError("Labeled cells are missing an outer fold assignment.")
        selected["outer_fold"] = selected["outer_fold"].astype(np.int8)

    ensure_directory(output_dir)
    sp.save_npz(output_dir / "counts.npz", counts, compressed=True)
    selected.to_parquet(output_dir / "cells.parquet", index=False)
    genes_table.to_parquet(output_dir / "genes.parquet", index=False)
    folds.to_csv(output_dir / "donor_folds.csv", index=False)

    manifest = {
        "created_at": utc_now(),
        "dataset_id": config.dataset_id,
        "role": config.role,
        "source": {
            "archive": str(archive),
            "archive_sha256": sha256_file(archive),
            "metadata": str(obs_path),
            "metadata_sha256": sha256_file(obs_path),
            "matrix_dir": str(matrix_dir),
            "matrix_shape": list(EXPECTED_MATRIX_SHAPES[config.modality]),
            "matrix_first_cell": matrix_cells[0],
            "matrix_last_cell": matrix_cells[-1],
        },
        "selection": {
            "tissue": config.tissue,
            "cell_type": config.cell_type,
            "excluded_subtissues": list(config.excluded_subtissues),
            "young_ages": list(config.young_ages),
            "old_ages": list(config.old_ages),
            "sensitivity_ages": list(config.sensitivity_ages),
            "cells_before_qc": int(config.expected_cells),
            "cells_after_qc": int(counts.shape[0]),
            "genes": int(counts.shape[1]),
            "mice": int(selected["mouse_id"].nunique()),
        },
        "outputs": {
            "counts_sha256": sha256_file(output_dir / "counts.npz"),
            "cells_sha256": sha256_file(output_dir / "cells.parquet"),
            "genes_sha256": sha256_file(output_dir / "genes.parquet"),
        },
        "runtime": runtime_info(),
    }
    write_json(manifest_path, manifest)
    return output_dir


def load_prepared_bundle(path: str | Path) -> PreparedBundle:
    directory = Path(path)
    required = ["counts.npz", "cells.parquet", "genes.parquet", "donor_folds.csv", "manifest.json"]
    missing = [name for name in required if not (directory / name).exists()]
    if missing:
        raise FileNotFoundError(f"Prepared bundle is incomplete: {missing}")
    counts = sp.load_npz(directory / "counts.npz").tocsr()
    cells = pd.read_parquet(directory / "cells.parquet")
    genes = pd.read_parquet(directory / "genes.parquet")
    folds = pd.read_csv(directory / "donor_folds.csv")
    if counts.shape != (len(cells), len(genes)):
        raise ValueError(f"Bundle shape mismatch: counts={counts.shape}, cells={len(cells)}, genes={len(genes)}")
    return PreparedBundle(counts=counts, cells=cells, genes=genes, folds=folds, directory=directory)

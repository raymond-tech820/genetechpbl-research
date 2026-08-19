"""Backed CSR h5ad reader for large single-cell count matrices."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


def _decode_strings(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind not in {"O", "S", "U"}:
        return values
    return np.asarray(
        [
            value.decode("utf-8")
            if isinstance(value, (bytes, np.bytes_))
            else str(value)
            for value in values
        ],
        dtype=object,
    )


def _read_series(root: h5py.File, dataset: h5py.Dataset) -> np.ndarray:
    values = dataset[:]
    categories_reference = dataset.attrs.get("categories")
    if categories_reference is not None:
        categories = _decode_strings(root[categories_reference][:])
        codes = np.asarray(values, dtype=int)
        output = np.full(len(codes), None, dtype=object)
        valid = codes >= 0
        output[valid] = categories[codes[valid]]
        return output
    return _decode_strings(np.asarray(values))


def read_dataframe_group(
    root: h5py.File,
    group_name: str,
    columns: list[str] | None = None,
) -> pd.DataFrame:
    group = root[group_name]
    index_name = str(group.attrs.get("_index", "_index"))
    available = [name for name in group.keys() if name != "__categories"]
    requested = (
        available
        if columns is None
        else list(dict.fromkeys([index_name, *columns]))
    )
    data: dict[str, np.ndarray] = {}
    for name in requested:
        if name in group and isinstance(group[name], h5py.Dataset):
            data[name] = _read_series(root, group[name])
    if index_name not in data:
        raise ValueError(
            f"Missing dataframe index `{index_name}` in h5ad group `{group_name}`."
        )
    index = pd.Index(data.pop(index_name).astype(str), name=index_name)
    return pd.DataFrame(data, index=index)


@dataclass
class BackedCSR:
    group: h5py.Group

    @property
    def shape(self) -> tuple[int, int]:
        shape = self.group.attrs["shape"]
        return int(shape[0]), int(shape[1])

    def rows(self, start: int, stop: int) -> sp.csr_matrix:
        indptr = np.asarray(self.group["indptr"][start : stop + 1], dtype=np.int64)
        first, last = int(indptr[0]), int(indptr[-1])
        data = np.asarray(self.group["data"][first:last])
        indices = np.asarray(self.group["indices"][first:last], dtype=np.int32)
        return sp.csr_matrix(
            (data, indices, indptr - first),
            shape=(stop - start, self.shape[1]),
        )

    def iter_chunks(self, chunk_size: int) -> Iterator[tuple[int, int, sp.csr_matrix]]:
        for start in range(0, self.shape[0], chunk_size):
            stop = min(start + chunk_size, self.shape[0])
            yield start, stop, self.rows(start, stop)


@dataclass
class H5adInput:
    path: Path
    handle: h5py.File
    counts: BackedCSR
    obs: pd.DataFrame
    gene_names: list[str]
    matrix_source: str
    gene_name_source: str

    def close(self) -> None:
        self.handle.close()

    def __enter__(self) -> "H5adInput":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        self.close()


def open_h5ad(
    path: str | Path,
    counts_layer: str = "counts",
    obs_columns: list[str] | None = None,
    gene_name_column: str = "gene_name",
) -> H5adInput:
    source = Path(path).resolve()
    handle = h5py.File(source, "r")
    try:
        if counts_layer in handle.get("layers", {}):
            matrix_group = handle["layers"][counts_layer]
            matrix_source = f"layers/{counts_layer}"
        else:
            matrix_group = handle["X"]
            matrix_source = "X"
        if not isinstance(matrix_group, h5py.Group):
            raise ValueError("Only sparse CSR h5ad count matrices are supported.")
        encoding = matrix_group.attrs.get("encoding-type", "")
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        if encoding != "csr_matrix":
            raise ValueError(f"Expected CSR counts, found `{encoding}`.")
        obs = read_dataframe_group(handle, "obs", columns=obs_columns)
        var = read_dataframe_group(handle, "var", columns=[gene_name_column])
        if gene_name_column in var.columns:
            gene_names = var[gene_name_column].astype(str).tolist()
            gene_name_source = f"var/{gene_name_column}"
        else:
            gene_names = var.index.astype(str).tolist()
            gene_name_source = "var index"
        counts = BackedCSR(matrix_group)
        if counts.shape != (len(obs), len(gene_names)):
            raise ValueError("h5ad matrix, obs, and var dimensions do not agree.")
        return H5adInput(
            path=source,
            handle=handle,
            counts=counts,
            obs=obs,
            gene_names=gene_names,
            matrix_source=matrix_source,
            gene_name_source=gene_name_source,
        )
    except Exception:
        handle.close()
        raise

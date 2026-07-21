"""Interpretable Youth Score baselines and technical-confound diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import joblib
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import SGDClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .features import TokenizerState, log_normalize_counts
from .utils import read_json, write_json


@dataclass
class BaselinePredictions:
    logits: dict[str, np.ndarray]


def _fit_logistic_1d(values: np.ndarray, labels: np.ndarray, weights: np.ndarray) -> LogisticRegression:
    model = LogisticRegression(solver="lbfgs", max_iter=2_000, random_state=20260717)
    model.fit(values.reshape(-1, 1), labels, sample_weight=weights)
    return model


def _signature_components(
    normalized_train: sp.csr_matrix,
    train_cells: pd.DataFrame,
    gene_count: int = 50,
) -> tuple[np.ndarray, np.ndarray]:
    donor_rows = []
    donor_labels = []
    for mouse_id, group in train_cells.groupby("mouse_id"):
        indices = group["_matrix_row"].to_numpy(dtype=int)
        donor_rows.append(np.asarray(normalized_train[indices].mean(axis=0)).ravel())
        donor_labels.append(int(group["label"].iloc[0]))
    donor_matrix = np.vstack(donor_rows)
    donor_labels_array = np.asarray(donor_labels)
    young = donor_matrix[donor_labels_array == 1]
    old = donor_matrix[donor_labels_array == 0]
    difference = young.mean(axis=0) - old.mean(axis=0)
    young_genes = np.argsort(difference, kind="stable")[-gene_count:]
    old_genes = np.argsort(difference, kind="stable")[:gene_count]
    return young_genes.astype(np.int32), old_genes.astype(np.int32)


def _signature_values(matrix: sp.csr_matrix, young: np.ndarray, old: np.ndarray) -> np.ndarray:
    young_mean = np.asarray(matrix[:, young].mean(axis=1)).ravel()
    old_mean = np.asarray(matrix[:, old].mean(axis=1)).ravel()
    return young_mean - old_mean


def _technical_features(cells: pd.DataFrame) -> np.ndarray:
    sex = cells["sex"].str.lower().map({"female": 0.0, "male": 1.0}).fillna(0.5).to_numpy()
    return np.column_stack(
        [
            np.log1p(cells["total_counts_matrix"].to_numpy(dtype=float)),
            np.log1p(cells["detected_genes_matrix"].to_numpy(dtype=float)),
            sex,
        ]
    )


class BaselineFold:
    """Fit and score the three baseline models for one leakage-safe fold."""

    def __init__(self, state: TokenizerState) -> None:
        self.state = state
        self.young_local: np.ndarray | None = None
        self.old_local: np.ndarray | None = None
        self.signature_model: LogisticRegression | None = None
        self.elastic_model: Pipeline | None = None
        self.technical_model: Pipeline | None = None

    def fit(
        self,
        counts: sp.csr_matrix,
        cells: pd.DataFrame,
        train_indices: np.ndarray,
        weights: np.ndarray,
    ) -> None:
        feature_indices = np.asarray(self.state.source_gene_indices, dtype=int)
        normalized = log_normalize_counts(counts[:, feature_indices])
        train_cells = cells.iloc[train_indices].copy()
        train_cells["_matrix_row"] = train_indices
        normalized_train_reference = normalized
        self.young_local, self.old_local = _signature_components(normalized_train_reference, train_cells)
        signature = _signature_values(normalized, self.young_local, self.old_local)
        labels = cells["label"].to_numpy(dtype=int)
        self.signature_model = _fit_logistic_1d(signature[train_indices], labels[train_indices], weights[train_indices])

        self.elastic_model = Pipeline(
            [
                ("scale", StandardScaler(with_mean=False)),
                (
                    "model",
                    SGDClassifier(
                        loss="log_loss",
                        penalty="elasticnet",
                        l1_ratio=0.5,
                        alpha=1e-4,
                        max_iter=3_000,
                        tol=1e-4,
                        random_state=20260717,
                        average=True,
                    ),
                ),
            ]
        )
        self.elastic_model.fit(normalized[train_indices], labels[train_indices], model__sample_weight=weights[train_indices])

        self.technical_model = Pipeline(
            [("scale", StandardScaler()), ("model", LogisticRegression(max_iter=2_000, random_state=20260717))]
        )
        technical = _technical_features(cells)
        self.technical_model.fit(technical[train_indices], labels[train_indices], model__sample_weight=weights[train_indices])

    def predict_logits(self, counts: sp.csr_matrix, cells: pd.DataFrame) -> BaselinePredictions:
        if any(model is None for model in (self.signature_model, self.elastic_model, self.technical_model)):
            raise RuntimeError("BaselineFold must be fitted before prediction.")
        feature_indices = np.asarray(self.state.source_gene_indices, dtype=int)
        normalized = log_normalize_counts(counts[:, feature_indices])
        signature = _signature_values(normalized, self.young_local, self.old_local)  # type: ignore[arg-type]
        return BaselinePredictions(
            logits={
                "gene_signature": self.signature_model.decision_function(signature.reshape(-1, 1)),  # type: ignore[union-attr]
                "elastic_net": self.elastic_model.decision_function(normalized),  # type: ignore[union-attr]
                "technical_only": self.technical_model.decision_function(_technical_features(cells)),  # type: ignore[union-attr]
            }
        )

    def save(self, directory: str | Path, all_gene_names: Sequence[str]) -> None:
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        joblib.dump(self.signature_model, directory / "signature_model.joblib")
        joblib.dump(self.elastic_model, directory / "elastic_model.joblib")
        joblib.dump(self.technical_model, directory / "technical_model.joblib")
        local_names = list(self.state.gene_names)
        write_json(
            directory / "baseline_state.json",
            {
                "feature_gene_names": local_names,
                "young_signature_gene_names": [local_names[i] for i in self.young_local],
                "old_signature_gene_names": [local_names[i] for i in self.old_local],
            },
        )

    @classmethod
    def load(cls, directory: str | Path, state: TokenizerState) -> "BaselineFold":
        directory = Path(directory)
        value = cls(state)
        metadata = read_json(directory / "baseline_state.json")
        name_to_local = {gene: i for i, gene in enumerate(state.gene_names)}
        value.young_local = np.asarray([name_to_local[g] for g in metadata["young_signature_gene_names"]], dtype=int)
        value.old_local = np.asarray([name_to_local[g] for g in metadata["old_signature_gene_names"]], dtype=int)
        value.signature_model = joblib.load(directory / "signature_model.joblib")
        value.elastic_model = joblib.load(directory / "elastic_model.joblib")
        value.technical_model = joblib.load(directory / "technical_model.joblib")
        return value

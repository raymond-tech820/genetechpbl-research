"""Probability calibration, donor aggregation, metrics, and model selection."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import expit
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    roc_auc_score,
)

from .constants import MODEL_COMPLEXITY_ORDER


@dataclass(frozen=True)
class TemperatureScaler:
    temperature: float = 1.0

    @classmethod
    def fit(cls, logits: np.ndarray, labels: np.ndarray, weights: np.ndarray | None = None) -> "TemperatureScaler":
        logits = np.asarray(logits, dtype=np.float64)
        labels = np.asarray(labels, dtype=np.float64)
        if weights is None:
            weights = np.ones_like(labels)
        else:
            weights = np.asarray(weights, dtype=np.float64)
        weights = weights / weights.sum()

        def objective(log_temperature: float) -> float:
            temperature = math.exp(log_temperature)
            probabilities = np.clip(expit(logits / temperature), 1e-7, 1 - 1e-7)
            losses = -(labels * np.log(probabilities) + (1 - labels) * np.log(1 - probabilities))
            return float(np.sum(weights * losses))

        result = minimize_scalar(objective, bounds=(-4.0, 4.0), method="bounded")
        return cls(temperature=float(math.exp(result.x)))

    def transform(self, logits: np.ndarray) -> np.ndarray:
        return np.clip(expit(np.asarray(logits, dtype=np.float64) / self.temperature), 1e-7, 1 - 1e-7)


def donor_balanced_weights(cells: pd.DataFrame) -> np.ndarray:
    labels = cells["label"].astype(int)
    donors_per_class = cells.groupby("label")["mouse_id"].nunique().to_dict()
    cells_per_donor = cells.groupby("mouse_id").size().to_dict()
    weights = np.asarray(
        [0.5 / (donors_per_class[label] * cells_per_donor[mouse]) for label, mouse in zip(labels, cells["mouse_id"])],
        dtype=np.float64,
    )
    return (weights / weights.mean()).astype(np.float32)


def select_decision_threshold(labels: Iterable[int], scores: Iterable[float]) -> float:
    """Choose a balanced-accuracy threshold using donor-level validation scores only."""

    labels_array = np.asarray(list(labels), dtype=int)
    scores_array = np.asarray(list(scores), dtype=float)
    if len(np.unique(labels_array)) != 2:
        return 0.5
    unique_scores = np.unique(scores_array)
    midpoints = (unique_scores[:-1] + unique_scores[1:]) / 2.0
    candidates = np.unique(np.concatenate(([0.5], midpoints)))
    ranked = []
    for threshold in candidates:
        predictions = (scores_array >= threshold).astype(int)
        ranked.append(
            (
                balanced_accuracy_score(labels_array, predictions),
                -abs(float(threshold) - 0.5),
                f1_score(labels_array, predictions, zero_division=0),
                -float(threshold),
                float(threshold),
            )
        )
    return max(ranked)[-1]


def aggregate_donor_scores(cell_predictions: pd.DataFrame) -> pd.DataFrame:
    required = {"model", "mouse_id", "label", "age_months", "youth_score"}
    missing = required - set(cell_predictions.columns)
    if missing:
        raise ValueError(f"Prediction table is missing required columns: {sorted(missing)}")
    aggregations = {
        "label": "first",
        "age_months": "first",
        "youth_score": "mean",
        "cell_id": "count",
    }
    if "predicted_age_months" in cell_predictions:
        aggregations["predicted_age_months"] = "mean"
    if "decision_threshold" in cell_predictions:
        aggregations["decision_threshold"] = "mean"
    donor = (
        cell_predictions.groupby(["model", "mouse_id"], as_index=False)
        .agg(aggregations)
        .rename(columns={"cell_id": "cell_count"})
    )
    return donor


def _safe_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    return float(roc_auc_score(labels, scores)) if len(np.unique(labels)) == 2 else float("nan")


def compute_binary_metrics(
    labels: Iterable[int],
    scores: Iterable[float],
    thresholds: Iterable[float] | float = 0.5,
) -> dict[str, float]:
    labels_array = np.asarray(list(labels), dtype=int)
    scores_array = np.asarray(list(scores), dtype=float)
    thresholds_array = np.asarray(thresholds, dtype=float)
    predictions = (scores_array >= thresholds_array).astype(int)
    return {
        "roc_auc": _safe_auc(labels_array, scores_array),
        "pr_auc": float(average_precision_score(labels_array, scores_array)),
        "balanced_accuracy": float(balanced_accuracy_score(labels_array, predictions)),
        "f1": float(f1_score(labels_array, predictions, zero_division=0)),
        "brier": float(brier_score_loss(labels_array, scores_array)),
    }


def stratified_bootstrap_auc(
    donor_predictions: pd.DataFrame,
    iterations: int = 2_000,
    seed: int = 20260717,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    groups = [group for _, group in donor_predictions.groupby("label")]
    if len(groups) != 2:
        return float("nan"), float("nan")
    estimates = []
    for _ in range(iterations):
        sample = pd.concat(
            [group.iloc[rng.integers(0, len(group), size=len(group))] for group in groups],
            ignore_index=True,
        )
        estimates.append(_safe_auc(sample["label"].to_numpy(), sample["youth_score"].to_numpy()))
    return float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975))


def summarize_models(donor_predictions: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for model, group in donor_predictions.groupby("model"):
        thresholds = group["decision_threshold"] if "decision_threshold" in group else 0.5
        metrics = compute_binary_metrics(group["label"], group["youth_score"], thresholds)
        low, high = stratified_bootstrap_auc(group)
        age_corr = float(group[["age_months", "youth_score"]].corr(method="spearman").iloc[0, 1])
        rows.append(
            {
                "model": model,
                **metrics,
                "roc_auc_ci_low": low,
                "roc_auc_ci_high": high,
                "age_spearman": age_corr,
                "donors": len(group),
            }
        )
    return pd.DataFrame(rows).sort_values(["roc_auc", "brier"], ascending=[False, True]).reset_index(drop=True)


def select_official_model(metrics: pd.DataFrame) -> str:
    eligible = metrics[metrics["model"].isin(MODEL_COMPLEXITY_ORDER)].copy()
    if eligible.empty:
        raise ValueError("No eligible Youth Score models were evaluated.")
    best_auc = eligible["roc_auc"].max()
    eligible = eligible[eligible["roc_auc"] >= best_auc - 0.02]
    best_brier = eligible["brier"].min()
    eligible = eligible[eligible["brier"] <= best_brier + 0.02].copy()
    eligible["complexity"] = eligible["model"].map(MODEL_COMPLEXITY_ORDER)
    return str(eligible.sort_values(["complexity", "model"]).iloc[0]["model"])


def direction_is_reversed(donor_predictions: pd.DataFrame) -> bool:
    means = donor_predictions.groupby("label")["youth_score"].mean()
    return bool(0 in means.index and 1 in means.index and means.loc[1] <= means.loc[0])

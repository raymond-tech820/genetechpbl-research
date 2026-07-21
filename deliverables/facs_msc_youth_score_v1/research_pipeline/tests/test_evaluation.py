from __future__ import annotations

import numpy as np
import pandas as pd

from youth_score.evaluation import (
    TemperatureScaler,
    compute_binary_metrics,
    select_decision_threshold,
    select_official_model,
)


def test_temperature_scaler_returns_probabilities() -> None:
    logits = np.array([-3.0, -1.0, 1.0, 3.0])
    labels = np.array([0, 0, 1, 1])
    scaler = TemperatureScaler.fit(logits, labels)
    probabilities = scaler.transform(logits)
    assert np.all((probabilities > 0) & (probabilities < 1))
    assert probabilities[0] < probabilities[-1]


def test_model_selection_prefers_simpler_near_tie() -> None:
    metrics = pd.DataFrame(
        [
            {"model": "gene_signature", "roc_auc": 0.80, "brier": 0.15},
            {"model": "elastic_net", "roc_auc": 0.81, "brier": 0.15},
            {"model": "gene_transformer", "roc_auc": 0.81, "brier": 0.16},
            {"model": "technical_only", "roc_auc": 0.99, "brier": 0.01},
        ]
    )
    assert select_official_model(metrics) == "gene_signature"


def test_validation_threshold_uses_donor_scores_without_degenerate_tie() -> None:
    labels = np.array([0, 1])
    ordered_scores = np.array([0.42, 0.58])
    threshold = select_decision_threshold(labels, ordered_scores)
    metrics = compute_binary_metrics(labels, ordered_scores, threshold)
    assert 0.42 < threshold <= 0.58
    assert metrics["balanced_accuracy"] == 1.0

    reversed_scores = np.array([0.58, 0.42])
    assert select_decision_threshold(labels, reversed_scores) == 0.5

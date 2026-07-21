from __future__ import annotations

import numpy as np
import pandas as pd

from youth_score.data import assign_donor_folds, gene_exclusion_reason, parse_age_months


def test_parse_age_months() -> None:
    assert parse_age_months("3m") == 3.0
    assert parse_age_months("24m") == 24.0


def test_gene_exclusions() -> None:
    assert gene_exclusion_reason("mt-Nd1") == "mitochondrial"
    assert gene_exclusion_reason("Rpl13") == "ribosomal"
    assert gene_exclusion_reason("Xist") == "sex_associated"
    assert gene_exclusion_reason("Kdm5d") == "y_chromosome"
    assert gene_exclusion_reason("Col1a1") == ""


def test_donor_folds_are_stratified() -> None:
    rows = []
    for label, prefix in [(1.0, "young"), (0.0, "old")]:
        for donor in range(5):
            for cell in range(3):
                rows.append({"mouse_id": f"{prefix}_{donor}", "label": label, "age_group": prefix})
    cells = pd.DataFrame(rows)
    folds = assign_donor_folds(cells, n_splits=5, seed=7)
    assert folds["outer_fold"].nunique() == 5
    assert folds.groupby("outer_fold")["label"].nunique().eq(2).all()
    assert not folds["mouse_id"].duplicated().any()


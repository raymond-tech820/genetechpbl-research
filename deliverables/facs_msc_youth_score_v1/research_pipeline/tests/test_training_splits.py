from __future__ import annotations

import pandas as pd

from youth_score.training import build_fold_indices


def test_train_validation_test_have_disjoint_donors() -> None:
    rows = []
    for label, age in [(1.0, 3.0), (0.0, 18.0)]:
        for fold in range(5):
            mouse = f"m_{int(label)}_{fold}"
            for cell in range(2):
                rows.append({"mouse_id": mouse, "label": label, "outer_fold": fold, "age_months": age})
    cells = pd.DataFrame(rows)
    split = build_fold_indices(cells, test_fold=0, seed=3)
    donor_sets = {name: set(cells.iloc[idx]["mouse_id"]) for name, idx in split.items()}
    assert donor_sets["train"].isdisjoint(donor_sets["validation"])
    assert donor_sets["train"].isdisjoint(donor_sets["test"])
    assert donor_sets["validation"].isdisjoint(donor_sets["test"])

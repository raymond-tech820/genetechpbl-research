from __future__ import annotations

import pandas as pd

from youth_score.external import infer_external_schema


def test_external_schema_inference() -> None:
    obs = pd.DataFrame(
        {
            "Age": ["Young", "Aged", "Young", "Aged"],
            "Condition": ["Control", "Control", "SOKM", "SOKM"],
            "donor_id": ["A", "B", "C", "D"],
        }
    )
    schema = infer_external_schema(obs)
    assert schema["age_column"] == "Age"
    assert schema["treatment_column"] == "Condition"
    assert schema["replicate_column"] == "donor_id"


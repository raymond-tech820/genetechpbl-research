import pandas as pd

from identity_score_v1.aggregation import aggregate_identity_scores


def _cell(arm: str, animal: str, value: float) -> dict:
    return {
        "cell_id": f"{arm}_{animal}",
        "age_raw": "old",
        "age_group": "aged",
        "treatment_raw": arm,
        "analysis_role": "SOKM" if arm == "Tg+/Dox+" else "control",
        "animal": animal,
        "known_animal": True,
        "state": "baseline",
        "sample": "sample",
        "batch": "batch",
        "msc_identity_core_primary": value,
        "msc_identity_core_rank": value,
    }


def test_reused_animal_labels_remain_distinct_unpaired_units() -> None:
    cells = pd.DataFrame(
        [
            _cell("Tg+/Dox+", "1", 0.1),
            _cell("Tg+/Dox-", "1", 0.2),
            _cell("Tg-/Dox+", "1", 0.3),
        ]
    )
    tables = aggregate_identity_scores(
        cells,
        dataset_id="test",
        species="Mus musculus",
        cell_population="MSC",
        unit_prefix="TEST",
        gene_coverage=1.0,
        expected_biological_units=3,
    )
    donors = tables["donors"]
    assert len(donors) == 3
    assert donors["biological_unit_id"].nunique() == 3
    assert set(donors["exact_treatment_arm"]) == {
        "Tg+/Dox+",
        "Tg+/Dox-",
        "Tg-/Dox+",
    }
    assert set(donors["cross_arm_pairing_status"]) == {"unpaired"}

import numpy as np
import pandas as pd
import scipy.sparse as sp

from identity_score_v1.scoring import prepare_identity_plan, score_chunk


def _identity_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_symbol": ["Identity"],
            "module": ["msc_identity_core"],
            "direction": [1],
            "included": [True],
        }
    )


def test_primary_direction_and_rank_scale_invariance() -> None:
    genes = ["Identity", *[f"Background{i}" for i in range(130)]]
    reference_means = np.linspace(0.1, 2.0, len(genes))
    plan, qc = prepare_identity_plan(
        genes,
        _identity_table(),
        _identity_table(),
        reference_means,
        n_control_sets=20,
        control_bins=8,
        seed=7,
    )
    counts = np.ones((2, len(genes)), dtype=np.float32)
    counts[0, 0] = 100
    sparse = sp.csr_matrix(counts)
    scores = score_chunk(sparse, plan, ["high", "low"])
    scaled = score_chunk(sparse * 9, plan, ["high", "low"])

    assert qc.iloc[0]["status"] == "pass"
    assert scores.loc[0, "msc_identity_core_primary"] > scores.loc[
        1, "msc_identity_core_primary"
    ]
    np.testing.assert_allclose(
        scores["msc_identity_core_rank"],
        scaled["msc_identity_core_rank"],
        rtol=0,
        atol=1e-7,
    )


def test_low_gene_coverage_fails_closed() -> None:
    identity = pd.DataFrame(
        {
            "gene_symbol": ["Present", "Missing"],
            "module": ["msc_identity_core", "msc_identity_core"],
            "direction": [1, 1],
            "included": [True, True],
        }
    )
    genes = ["Present", *[f"Background{i}" for i in range(130)]]
    plan, qc = prepare_identity_plan(
        genes,
        identity,
        identity,
        np.linspace(0.1, 2.0, len(genes)),
        minimum_gene_coverage=0.70,
    )
    assert plan.status == "insufficient_gene_coverage"
    assert qc.iloc[0]["gene_coverage"] == 0.5

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from youth_score.config import DatasetConfig, ModelConfig, TrainingConfig
from youth_score.scoring import YouthScoreEnsemble
from youth_score.training import run_cross_validation
from youth_score.utils import write_json


def test_cpu_two_epoch_end_to_end(tmp_path: Path) -> None:
    rng = np.random.default_rng(20260717)
    genes = 96
    cells_per_donor = 6
    records = []
    matrices = []
    donor_rows = []
    for label, age, prefix in ((1, "3m", "young"), (0, "18m", "old")):
        for donor_index in range(4):
            donor = f"{prefix}_{donor_index}"
            outer_fold = donor_index % 2
            donor_rows.append({"mouse_id": donor, "label": label, "outer_fold": outer_fold})
            values = rng.poisson(0.35, size=(cells_per_donor, genes)).astype(np.float32)
            values[:, :12] += label * rng.poisson(2.0, size=(cells_per_donor, 12))
            values[:, 12:24] += (1 - label) * rng.poisson(2.0, size=(cells_per_donor, 12))
            matrices.append(values)
            for cell_index in range(cells_per_donor):
                records.append(
                    {
                        "cell_id": f"{donor}_cell_{cell_index}",
                        "mouse_id": donor,
                        "age": age,
                        "age_months": 3.0 if label else 18.0,
                        "age_group": "young" if label else "old",
                        "label": label,
                        "sex": "male" if donor_index % 2 else "female",
                        "outer_fold": outer_fold,
                    }
                )

    counts = sp.csr_matrix(np.vstack(matrices))
    cells = pd.DataFrame(records)
    cells["total_counts_matrix"] = np.asarray(counts.sum(axis=1)).ravel()
    cells["detected_genes_matrix"] = np.diff(counts.indptr)
    gene_table = pd.DataFrame(
        {
            "gene_name": [f"Gene{index}" for index in range(genes)],
            "eligible_primary": np.ones(genes, dtype=bool),
        }
    )

    config = DatasetConfig(
        path=tmp_path / "synthetic.yaml",
        project_root=tmp_path,
        dataset_id="synthetic_cpu",
        role="primary",
        modality="facs",
        tissue="Synthetic",
        cell_type="synthetic cell",
        young_ages=("3m",),
        old_ages=("18m",),
        sensitivity_ages=(),
        expected_cells=len(cells),
        expected_young_cells=int(cells["label"].sum()),
        expected_old_cells=int((1 - cells["label"]).sum()),
        expected_young_mice=4,
        expected_old_mice=4,
        minimum_detected_genes=1,
        folds=2,
        seeds=(20260717,),
        feature_count=64,
        sequence_length=32,
        expression_bins=8,
        model=ModelConfig(
            d_model=32,
            heads=4,
            layers=1,
            feedforward=64,
            dropout=0.1,
            token_dropout=0.1,
            age_loss_weight=0.25,
        ),
        training=TrainingConfig(
            batch_size=16,
            learning_rate=3e-4,
            weight_decay=1e-2,
            max_epochs=2,
            patience=2,
            gradient_clip=1.0,
            num_workers=0,
        ),
    )
    bundle = config.processed_dir
    bundle.mkdir(parents=True)
    sp.save_npz(bundle / "counts.npz", counts)
    cells.to_parquet(bundle / "cells.parquet", index=False)
    gene_table.to_parquet(bundle / "genes.parquet", index=False)
    pd.DataFrame(donor_rows).to_csv(bundle / "donor_folds.csv", index=False)
    write_json(bundle / "manifest.json", {"dataset_id": config.dataset_id, "synthetic": True})

    output = run_cross_validation(config, device_request="cpu", force=True)
    selection = output / "selection.json"
    assert selection.exists()
    assert len(list(output.glob("folds/fold_*/transformer/seed_*/model.safetensors"))) == 2

    scorer = YouthScoreEnsemble(output, device="cpu", batch_size=16)
    first = scorer.score_counts(counts, gene_table["gene_name"].tolist(), cells["cell_id"].tolist())
    second = scorer.score_counts(counts, gene_table["gene_name"].tolist(), cells["cell_id"].tolist())
    assert np.allclose(first["youth_score"], second["youth_score"])
    assert np.all((first["youth_score"] >= 0.0) & (first["youth_score"] <= 1.0))


from __future__ import annotations

import torch

from youth_score.model import GeneTransformer, TransformerShape


def test_transformer_forward() -> None:
    shape = TransformerShape(
        vocabulary_size=34,
        expression_bin_count=8,
        sequence_length=12,
        d_model=16,
        heads=4,
        layers=1,
        feedforward=32,
        dropout=0.0,
        token_dropout=0.0,
    )
    model = GeneTransformer(shape)
    gene_ids = torch.randint(1, 33, (5, 13))
    expression_ids = torch.randint(1, 9, (5, 13))
    rank_ids = torch.arange(13).repeat(5, 1)
    mask = torch.ones((5, 13), dtype=torch.bool)
    logits, ages = model(gene_ids, expression_ids, rank_ids, mask)
    assert logits.shape == (5,)
    assert ages.shape == (5,)


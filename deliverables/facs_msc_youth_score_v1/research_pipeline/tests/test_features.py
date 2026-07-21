from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from youth_score.features import build_tokenizer_state, encode_named_counts, encode_tms_counts


def test_tokenizer_shapes_and_gene_overlap() -> None:
    rng = np.random.default_rng(4)
    dense = rng.poisson(1.0, size=(20, 30)).astype(np.uint32)
    dense[dense < 1] = 0
    counts = sp.csr_matrix(dense)
    genes = [f"Gene{i}" for i in range(30)]
    state = build_tokenizer_state(
        counts[:12], genes, np.ones(30, dtype=bool), feature_count=20, sequence_length=8, expression_bins=4, seed=1
    )
    encoded = encode_tms_counts(counts, state)
    assert encoded["gene_ids"].shape == (20, 9)
    assert np.all(encoded["gene_ids"][:, 0] == state.cls_token_id)
    external, overlap = encode_named_counts(counts, genes, state)
    assert overlap == 1.0
    np.testing.assert_array_equal(encoded["attention_mask"], external["attention_mask"])


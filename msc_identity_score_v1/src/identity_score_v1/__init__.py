"""Frozen mouse MSC Identity Score v1."""

from .aggregation import aggregate_identity_scores
from .pipeline import run_pipeline
from .scoring import IdentityPlan, prepare_identity_plan, score_chunk

__all__ = [
    "IdentityPlan",
    "aggregate_identity_scores",
    "prepare_identity_plan",
    "run_pipeline",
    "score_chunk",
]

__version__ = "1.0.0"

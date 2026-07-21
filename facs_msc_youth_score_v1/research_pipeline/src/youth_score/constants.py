"""Shared project constants."""

from pathlib import Path

EXPECTED_MATRIX_SHAPES = {
    "facs": (22_966, 110_824),
    "droplet": (20_138, 245_389),
}

ARCHIVE_NAMES = {
    "facs": "tabula-muris-senis-facs_bpcells.tar.zst",
    "droplet": "tabula-muris-senis-droplet_bpcells.tar.zst",
}

MATRIX_DIR_NAMES = {
    "facs": "tabula-muris-senis-facs_bpcells",
    "droplet": "tabula-muris-senis-droplet_bpcells",
}

OBS_NAMES = {
    "facs": "facs_obs.csv",
    "droplet": "droplet_obs.csv",
}

MODEL_COMPLEXITY_ORDER = {
    "gene_signature": 0,
    "elastic_net": 1,
    "gene_transformer": 2,
}

DEFAULT_ROOT = Path(__file__).resolve().parents[2]


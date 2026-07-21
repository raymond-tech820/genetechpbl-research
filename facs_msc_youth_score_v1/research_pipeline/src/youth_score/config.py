"""Configuration loading and validation."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True)
class ModelConfig:
    d_model: int = 128
    heads: int = 8
    layers: int = 2
    feedforward: int = 512
    dropout: float = 0.2
    token_dropout: float = 0.15
    age_loss_weight: float = 0.25


@dataclass(frozen=True)
class TrainingConfig:
    batch_size: int = 64
    learning_rate: float = 3e-4
    weight_decay: float = 1e-2
    max_epochs: int = 80
    patience: int = 12
    gradient_clip: float = 1.0
    num_workers: int = 0


@dataclass(frozen=True)
class DatasetConfig:
    path: Path
    project_root: Path
    dataset_id: str
    role: str
    modality: str
    tissue: str
    cell_type: str
    young_ages: tuple[str, ...]
    old_ages: tuple[str, ...]
    sensitivity_ages: tuple[str, ...]
    expected_cells: int
    expected_young_cells: int
    expected_old_cells: int
    expected_young_mice: int
    expected_old_mice: int
    minimum_detected_genes: int = 500
    folds: int = 5
    seeds: tuple[int, ...] = (20260717, 20260718, 20260719)
    feature_count: int = 4096
    sequence_length: int = 256
    expression_bins: int = 64
    model: ModelConfig = field(default_factory=ModelConfig)
    training: TrainingConfig = field(default_factory=TrainingConfig)

    @property
    def processed_dir(self) -> Path:
        return self.project_root / "data" / "processed" / "youth_score" / self.dataset_id

    @property
    def output_dir(self) -> Path:
        return self.project_root / "outputs" / "youth_score" / self.dataset_id


def _resolve_root(config_path: Path, raw_root: str) -> Path:
    root = Path(raw_root)
    if root.is_absolute():
        return root.resolve()
    return (config_path.resolve().parents[1] / root).resolve()


def load_dataset_config(path: str | Path) -> DatasetConfig:
    config_path = Path(path).resolve()
    with config_path.open("r", encoding="utf-8") as handle:
        raw: dict[str, Any] = yaml.safe_load(handle)
    model = ModelConfig(**raw.pop("model", {}))
    training = TrainingConfig(**raw.pop("training", {}))
    project_root = _resolve_root(config_path, raw.pop("project_root", "."))
    for key in ("young_ages", "old_ages", "sensitivity_ages", "seeds"):
        raw[key] = tuple(raw.get(key, ()))
    config = DatasetConfig(
        path=config_path,
        project_root=project_root,
        model=model,
        training=training,
        **raw,
    )
    if config.modality not in {"facs", "droplet"}:
        raise ValueError(f"Unsupported modality: {config.modality}")
    if set(config.young_ages) & set(config.old_ages):
        raise ValueError("Young and old age sets must be disjoint.")
    if config.model.d_model % config.model.heads != 0:
        raise ValueError("d_model must be divisible by heads.")
    return config


def load_yaml(path: str | Path) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


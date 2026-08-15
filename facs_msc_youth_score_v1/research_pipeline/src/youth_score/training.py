"""Donor-aware cross-validation and model training orchestration."""

from __future__ import annotations

import copy
import logging
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader

from .baselines import BaselineFold
from .config import DatasetConfig
from .data import load_prepared_bundle
from .evaluation import (
    TemperatureScaler,
    aggregate_donor_scores,
    direction_is_reversed,
    donor_balanced_weights,
    select_decision_threshold,
    select_official_model,
    summarize_models,
)
from .features import (
    TokenizerState,
    build_tokenizer_state,
    encode_tms_counts,
    load_encoded,
    save_encoded,
)
from .model import EncodedCellDataset, GeneTransformer, TransformerShape
from .utils import ensure_directory, runtime_info, set_global_seed, stable_hash, utc_now, write_json


LOGGER = logging.getLogger(__name__)


def _choose_validation_donors(cells: pd.DataFrame, test_fold: int, seed: int) -> set[str]:
    candidates = cells.loc[cells["outer_fold"].ne(test_fold), ["mouse_id", "label"]].drop_duplicates()
    selected: set[str] = set()
    for label, group in candidates.groupby("label"):
        donor = min(group["mouse_id"].astype(str), key=lambda value: stable_hash(value, seed + test_fold * 101))
        selected.add(donor)
    if len(selected) != 2:
        raise RuntimeError(f"Expected one validation donor per class, observed {selected}.")
    return selected


def build_fold_indices(cells: pd.DataFrame, test_fold: int, seed: int) -> dict[str, np.ndarray]:
    validation_donors = _choose_validation_donors(cells, test_fold, seed)
    test = cells["outer_fold"].eq(test_fold).to_numpy()
    validation = cells["mouse_id"].isin(validation_donors).to_numpy()
    train = ~(test | validation)
    result = {
        "train": np.flatnonzero(train),
        "validation": np.flatnonzero(validation),
        "test": np.flatnonzero(test),
    }
    donor_sets = {key: set(cells.iloc[idx]["mouse_id"]) for key, idx in result.items()}
    if donor_sets["train"] & donor_sets["validation"] or donor_sets["train"] & donor_sets["test"]:
        raise RuntimeError("Donor leakage was detected between training and held-out splits.")
    if donor_sets["validation"] & donor_sets["test"]:
        raise RuntimeError("Donor leakage was detected between validation and test splits.")
    for split, indices in result.items():
        if cells.iloc[indices]["label"].nunique() != 2:
            raise RuntimeError(f"Split {split} does not contain both age classes in outer fold {test_fold}.")
    return result


def _weights_for_split(cells: pd.DataFrame, indices: np.ndarray) -> np.ndarray:
    result = np.zeros(len(cells), dtype=np.float32)
    result[indices] = donor_balanced_weights(cells.iloc[indices])
    return result


def _device_from_request(requested: str) -> torch.device:
    if requested == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if requested == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but torch.cuda.is_available() is false.")
    return torch.device(requested)


@torch.inference_mode()
def _predict_transformer(
    model: GeneTransformer,
    encoded: dict[str, np.ndarray],
    labels: np.ndarray,
    ages_scaled: np.ndarray,
    weights: np.ndarray,
    indices: np.ndarray,
    batch_size: int,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    dataset = EncodedCellDataset(encoded, labels, ages_scaled, weights, indices)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False, num_workers=0)
    model.eval()
    all_indices: list[np.ndarray] = []
    all_logits: list[np.ndarray] = []
    all_ages: list[np.ndarray] = []
    use_amp = device.type == "cuda"
    for batch in loader:
        with torch.amp.autocast(device_type=device.type, dtype=torch.bfloat16, enabled=use_amp):
            logits, age_output = model(
                batch["gene_ids"].to(device),
                batch["expression_ids"].to(device),
                batch["rank_ids"].to(device),
                batch["attention_mask"].to(device),
            )
        all_indices.append(batch["index"].numpy())
        all_logits.append(logits.float().cpu().numpy())
        all_ages.append(age_output.float().cpu().numpy())
    return np.concatenate(all_indices), np.concatenate(all_logits), np.concatenate(all_ages)


def _validation_loss(
    model: GeneTransformer,
    encoded: dict[str, np.ndarray],
    labels: np.ndarray,
    ages_scaled: np.ndarray,
    weights: np.ndarray,
    indices: np.ndarray,
    config: DatasetConfig,
    device: torch.device,
) -> float:
    ordered, logits, age_output = _predict_transformer(
        model,
        encoded,
        labels,
        ages_scaled,
        weights,
        indices,
        config.training.batch_size,
        device,
    )
    local_weights = weights[ordered].astype(np.float64)
    local_weights /= local_weights.sum()
    binary = F.binary_cross_entropy_with_logits(
        torch.from_numpy(logits), torch.from_numpy(labels[ordered]), reduction="none"
    ).numpy()
    age = F.smooth_l1_loss(
        torch.from_numpy(age_output), torch.from_numpy(ages_scaled[ordered]), reduction="none"
    ).numpy()
    return float(np.sum(local_weights * (binary + config.model.age_loss_weight * age)))


def _train_transformer_seed(
    config: DatasetConfig,
    encoded: dict[str, np.ndarray],
    cells: pd.DataFrame,
    indices: dict[str, np.ndarray],
    state: TokenizerState,
    seed: int,
    directory: Path,
    device: torch.device,
    max_epochs_override: int | None = None,
    force: bool = False,
) -> pd.DataFrame:
    prediction_path = directory / "test_predictions.parquet"
    if prediction_path.exists() and not force:
        return pd.read_parquet(prediction_path)

    ensure_directory(directory)
    set_global_seed(seed)
    labels = cells["label"].to_numpy(dtype=np.float32)
    ages_scaled = (cells["age_months"].to_numpy(dtype=np.float32) / 30.0).astype(np.float32)
    train_weights = _weights_for_split(cells, indices["train"])
    validation_weights = _weights_for_split(cells, indices["validation"])
    shape = TransformerShape(
        vocabulary_size=state.vocabulary_size,
        expression_bin_count=state.expression_bin_count,
        sequence_length=state.sequence_length,
        d_model=config.model.d_model,
        heads=config.model.heads,
        layers=config.model.layers,
        feedforward=config.model.feedforward,
        dropout=config.model.dropout,
        token_dropout=config.model.token_dropout,
    )
    model = GeneTransformer(shape).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=config.training.learning_rate, weight_decay=config.training.weight_decay
    )
    train_dataset = EncodedCellDataset(
        encoded, labels, ages_scaled, train_weights, indices["train"]
    )
    generator = torch.Generator().manual_seed(seed)
    loader = DataLoader(
        train_dataset,
        batch_size=config.training.batch_size,
        shuffle=True,
        num_workers=config.training.num_workers,
        generator=generator,
        pin_memory=device.type == "cuda",
    )
    max_epochs = max_epochs_override or config.training.max_epochs
    total_steps = max(1, max_epochs * len(loader))
    warmup_steps = max(1, int(total_steps * 0.1))

    def learning_rate(step: int) -> float:
        if step < warmup_steps:
            return max(step, 1) / warmup_steps
        progress = (step - warmup_steps) / max(total_steps - warmup_steps, 1)
        return 0.5 * (1.0 + np.cos(np.pi * min(progress, 1.0)))

    scheduler = torch.optim.lr_scheduler.LambdaLR(optimizer, learning_rate)
    best_loss = float("inf")
    best_state: dict[str, torch.Tensor] | None = None
    patience = 0
    global_step = 0
    history: list[dict[str, float | int]] = []
    use_amp = device.type == "cuda"

    for epoch in range(max_epochs):
        model.train()
        epoch_loss = 0.0
        epoch_weight = 0.0
        for batch in loader:
            optimizer.zero_grad(set_to_none=True)
            with torch.amp.autocast(device_type=device.type, dtype=torch.bfloat16, enabled=use_amp):
                logits, age_output = model(
                    batch["gene_ids"].to(device),
                    batch["expression_ids"].to(device),
                    batch["rank_ids"].to(device),
                    batch["attention_mask"].to(device),
                )
                batch_weights = batch["weight"].to(device)
                binary = F.binary_cross_entropy_with_logits(logits, batch["label"].to(device), reduction="none")
                age = F.smooth_l1_loss(age_output, batch["age"].to(device), reduction="none")
                losses = binary + config.model.age_loss_weight * age
                loss = torch.sum(losses * batch_weights) / torch.sum(batch_weights)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), config.training.gradient_clip)
            optimizer.step()
            scheduler.step()
            global_step += 1
            epoch_loss += float(torch.sum(losses.detach() * batch_weights).cpu())
            epoch_weight += float(torch.sum(batch_weights).cpu())

        validation_loss = _validation_loss(
            model,
            encoded,
            labels,
            ages_scaled,
            validation_weights,
            indices["validation"],
            config,
            device,
        )
        history.append(
            {
                "epoch": epoch + 1,
                "train_loss": epoch_loss / max(epoch_weight, 1e-8),
                "validation_loss": validation_loss,
                "learning_rate": optimizer.param_groups[0]["lr"],
            }
        )
        if validation_loss < best_loss - 1e-5:
            best_loss = validation_loss
            best_state = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
            patience = 0
        else:
            patience += 1
            if patience >= min(config.training.patience, max_epochs):
                break

    if best_state is None:
        raise RuntimeError("Transformer training did not produce a checkpoint.")
    model.load_state_dict(best_state)
    model.to(device)
    model.save(directory)
    pd.DataFrame(history).to_csv(directory / "history.csv", index=False)

    validation_order, validation_logits, _ = _predict_transformer(
        model,
        encoded,
        labels,
        ages_scaled,
        validation_weights,
        indices["validation"],
        config.training.batch_size,
        device,
    )
    temperature = TemperatureScaler.fit(
        validation_logits,
        labels[validation_order],
        validation_weights[validation_order],
    )
    validation_frame = cells.iloc[validation_order][["mouse_id", "label"]].copy()
    validation_frame["youth_score"] = temperature.transform(validation_logits)
    validation_donors = validation_frame.groupby("mouse_id", as_index=False).agg(
        label=("label", "first"), youth_score=("youth_score", "mean")
    )
    decision_threshold = select_decision_threshold(
        validation_donors["label"], validation_donors["youth_score"]
    )
    test_order, test_logits, test_age = _predict_transformer(
        model,
        encoded,
        labels,
        ages_scaled,
        np.ones(len(cells), dtype=np.float32),
        indices["test"],
        config.training.batch_size,
        device,
    )
    predictions = cells.iloc[test_order][
        ["cell_id", "mouse_id", "age", "age_months", "age_group", "label", "sex", "outer_fold"]
    ].copy()
    predictions["model"] = "gene_transformer"
    predictions["seed"] = seed
    predictions["logit"] = test_logits
    predictions["youth_score"] = temperature.transform(test_logits)
    predictions["predicted_age_months"] = test_age * 30.0
    predictions["decision_threshold"] = decision_threshold
    predictions.to_parquet(prediction_path, index=False)
    write_json(
        directory / "calibration.json",
        {
            "temperature": temperature.temperature,
            "decision_threshold": decision_threshold,
            "best_validation_loss": best_loss,
            "epochs": len(history),
        },
    )
    return predictions


def _train_baselines_fold(
    counts,
    cells: pd.DataFrame,
    genes: pd.DataFrame,
    indices: dict[str, np.ndarray],
    state: TokenizerState,
    directory: Path,
    force: bool = False,
) -> pd.DataFrame:
    ensure_directory(directory)
    prediction_path = directory / "baseline_test_predictions.parquet"
    if prediction_path.exists() and not force:
        return pd.read_parquet(prediction_path)
    weights = _weights_for_split(cells, indices["train"])
    validation_weights = _weights_for_split(cells, indices["validation"])
    baseline = BaselineFold(state)
    baseline.fit(counts, cells, indices["train"], weights)
    all_logits = baseline.predict_logits(counts, cells).logits
    calibration: dict[str, float] = {}
    thresholds: dict[str, float] = {}
    frames = []
    for model_name, logits in all_logits.items():
        scaler = TemperatureScaler.fit(
            logits[indices["validation"]],
            cells.iloc[indices["validation"]]["label"].to_numpy(dtype=int),
            validation_weights[indices["validation"]],
        )
        calibration[model_name] = scaler.temperature
        validation_frame = cells.iloc[indices["validation"]][["mouse_id", "label"]].copy()
        validation_frame["youth_score"] = scaler.transform(logits[indices["validation"]])
        validation_donors = validation_frame.groupby("mouse_id", as_index=False).agg(
            label=("label", "first"), youth_score=("youth_score", "mean")
        )
        thresholds[model_name] = select_decision_threshold(
            validation_donors["label"], validation_donors["youth_score"]
        )
        frame = cells.iloc[indices["test"]][
            ["cell_id", "mouse_id", "age", "age_months", "age_group", "label", "sex", "outer_fold"]
        ].copy()
        frame["model"] = model_name
        frame["seed"] = -1
        frame["logit"] = logits[indices["test"]]
        frame["youth_score"] = scaler.transform(logits[indices["test"]])
        frame["predicted_age_months"] = np.nan
        frame["decision_threshold"] = thresholds[model_name]
        frames.append(frame)
    result = pd.concat(frames, ignore_index=True)
    result.to_parquet(prediction_path, index=False)
    baseline.save(directory, genes["gene_name"].astype(str).tolist())
    write_json(
        directory / "baseline_calibration.json",
        {
            model_name: {
                "temperature": calibration[model_name],
                "decision_threshold": thresholds[model_name],
            }
            for model_name in calibration
        },
    )
    return result


def _average_seed_predictions(predictions: pd.DataFrame) -> pd.DataFrame:
    transformer = predictions[predictions["model"].eq("gene_transformer")]
    others = predictions[~predictions["model"].eq("gene_transformer")]
    keys = ["cell_id", "mouse_id", "age", "age_months", "age_group", "label", "sex", "outer_fold", "model"]
    transformer_average = (
        transformer.groupby(keys, as_index=False, dropna=False)
        .agg(
            youth_score=("youth_score", "mean"),
            predicted_age_months=("predicted_age_months", "mean"),
            decision_threshold=("decision_threshold", "mean"),
            logit=("logit", "mean"),
        )
    )
    others = others.drop(columns=["seed"], errors="ignore")
    return pd.concat([transformer_average, others], ignore_index=True, sort=False)


def run_cross_validation(
    config: DatasetConfig,
    device_request: str = "auto",
    force: bool = False,
    max_epochs_override: int | None = None,
    folds_override: list[int] | None = None,
    seeds_override: list[int] | None = None,
    feature_mode: str = "primary",
) -> Path:
    """Run all donor-aware folds, baselines, calibration, and model selection."""

    if config.role == "sensitivity":
        raise ValueError("Sensitivity bundles are scored by a trained model and cannot be trained directly.")
    if feature_mode not in {"primary", "all"}:
        raise ValueError(f"Unsupported feature mode: {feature_mode}")
    bundle = load_prepared_bundle(config.processed_dir)
    cells = bundle.cells.copy().reset_index(drop=True)
    if cells["label"].isna().any():
        raise ValueError("Training bundles must contain only young/old labels.")
    output_dir = ensure_directory(
        config.output_dir
        if feature_mode == "primary"
        else config.output_dir / "all_genes_sensitivity"
    )
    final_path = output_dir / "selection.json"
    if final_path.exists() and not force and folds_override is None and seeds_override is None:
        return output_dir

    device = _device_from_request(device_request)
    folds = folds_override if folds_override is not None else list(range(config.folds))
    seeds = seeds_override if seeds_override is not None else list(config.seeds)
    all_predictions: list[pd.DataFrame] = []
    labels = cells["label"].to_numpy(dtype=np.float32)

    for fold in folds:
        LOGGER.info("Training %s outer fold %d", config.dataset_id, fold)
        fold_dir = ensure_directory(output_dir / "folds" / f"fold_{fold}")
        split_indices = build_fold_indices(cells, fold, config.seeds[0])
        write_json(
            fold_dir / "split.json",
            {
                split: sorted(cells.iloc[idx]["mouse_id"].unique().tolist())
                for split, idx in split_indices.items()
            },
        )

        state_path = fold_dir / "tokenizer.json"
        encoded_path = fold_dir / "tokenized.npz"
        if state_path.exists():
            state = TokenizerState.load(state_path)
        else:
            state = build_tokenizer_state(
                train_counts=bundle.counts[split_indices["train"]],
                all_gene_names=bundle.genes["gene_name"].astype(str).tolist(),
                eligible=(
                    bundle.genes["eligible_primary"].to_numpy(dtype=bool)
                    if feature_mode == "primary"
                    else np.ones(len(bundle.genes), dtype=bool)
                ),
                feature_count=config.feature_count,
                sequence_length=config.sequence_length,
                expression_bins=config.expression_bins,
                seed=config.seeds[0] + fold,
            )
            state.save(state_path)
        if encoded_path.exists():
            encoded = load_encoded(encoded_path)
        else:
            encoded = encode_tms_counts(bundle.counts, state)
            save_encoded(encoded_path, encoded)

        baseline_predictions = _train_baselines_fold(
            bundle.counts,
            cells,
            bundle.genes,
            split_indices,
            state,
            fold_dir / "baselines",
            force=force,
        )
        all_predictions.append(baseline_predictions)

        for seed in seeds:
            seed_predictions = _train_transformer_seed(
                config,
                encoded,
                cells,
                split_indices,
                state,
                seed,
                fold_dir / "transformer" / f"seed_{seed}",
                device,
                max_epochs_override=max_epochs_override,
                force=force,
            )
            all_predictions.append(seed_predictions)

    raw_predictions = pd.concat(
        [frame.dropna(axis=1, how="all") for frame in all_predictions],
        ignore_index=True,
    )
    oof = _average_seed_predictions(raw_predictions)
    if folds_override is None and oof.groupby("model")["cell_id"].nunique().min() != len(cells):
        raise RuntimeError("Full cross-validation did not produce exactly one OOF prediction per cell and model.")
    oof.to_parquet(output_dir / "oof_cell_scores.parquet", index=False)
    donor = aggregate_donor_scores(oof)
    donor.to_csv(output_dir / "oof_donor_scores.csv", index=False)
    metrics = summarize_models(donor)
    metrics.to_csv(output_dir / "model_metrics.csv", index=False)
    official = select_official_model(metrics)

    official_donors = donor[donor["model"].eq(official)]
    technical_donors = donor[donor["model"].eq("technical_only")]
    official_auc = float(metrics.loc[metrics["model"].eq(official), "roc_auc"].iloc[0])
    technical_auc = float(metrics.loc[metrics["model"].eq("technical_only"), "roc_auc"].iloc[0])
    male_only = official_donors[official_donors["mouse_id"].isin(cells.loc[cells["sex"].eq("male"), "mouse_id"])]
    age_3_18 = official_donors[official_donors["age_months"].isin([3.0, 18.0])]
    sensitivity = {
        "male_only_reversed": direction_is_reversed(male_only) if male_only["label"].nunique() == 2 else None,
        "age_3m_vs_18m_reversed": direction_is_reversed(age_3_18) if age_3_18["label"].nunique() == 2 else None,
        "technical_auc": technical_auc,
        "official_auc": official_auc,
    }
    confound_limited = (
        technical_auc >= official_auc - 0.02
        or sensitivity["male_only_reversed"] is True
        or sensitivity["age_3m_vs_18m_reversed"] is True
    )
    status = "confound_limited" if confound_limited else "internally_supported"
    selection = {
        "created_at": utc_now(),
        "dataset_id": config.dataset_id,
        "feature_mode": feature_mode,
        "official_model": official,
        "status": status,
        "selection_rule": "Highest donor OOF ROC-AUC; within 0.02 use lower Brier; remaining ties prefer simpler models.",
        "sensitivity": sensitivity,
        "runtime": runtime_info(),
        "config": {
            "dataset": asdict(config),
        },
    }
    selection["config"]["dataset"]["path"] = str(config.path)
    selection["config"]["dataset"]["project_root"] = str(config.project_root)
    write_json(final_path, selection)
    return output_dir

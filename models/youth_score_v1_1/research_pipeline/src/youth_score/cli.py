"""Command-line interface for the complete Youth Score workflow."""

from __future__ import annotations

import argparse
import json
import logging
import shutil
import sys
from pathlib import Path

import torch  # Import torch before h5py/anndata to avoid a Windows CUDA DLL load-order conflict.
import anndata as ad
import pandas as pd
import scipy.sparse as sp

from .config import load_dataset_config
from .data import load_prepared_bundle, prepare_tms_bundle
from .external import download_gse176206, validate_gse176206
from .reporting import generate_project_summary, generate_training_report, score_limb_droplet_sensitivity
from .scoring import YouthScoreEnsemble
from .training import run_cross_validation
from .utils import runtime_info


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIGS = {
    "scat": ROOT / "configs" / "scat_primary.yaml",
    "limb": ROOT / "configs" / "limb_secondary.yaml",
    "sensitivity": ROOT / "configs" / "limb_droplet_sensitivity.yaml",
    "external": ROOT / "configs" / "external_gse176206.yaml",
}


def _parse_int_list(value: str | None) -> list[int] | None:
    return None if not value else [int(item.strip()) for item in value.split(",") if item.strip()]


def command_doctor(_: argparse.Namespace) -> int:
    result = runtime_info()
    result["project_root"] = str(ROOT)
    result["free_disk_gb"] = round(shutil.disk_usage(ROOT).free / (1024**3), 2)
    result["tms_source_files"] = {
        name: path.exists()
        for name, path in {
            "facs_obs": ROOT / "data" / "tms_compact" / "facs_obs.csv",
            "droplet_obs": ROOT / "data" / "tms_compact" / "droplet_obs.csv",
            "facs_matrix": ROOT / "data" / "tms_compact" / "tabula-muris-senis-facs_bpcells.tar.zst",
            "droplet_matrix": ROOT / "data" / "tms_compact" / "tabula-muris-senis-droplet_bpcells.tar.zst",
        }.items()
    }
    print(json.dumps(result, indent=2))
    if not all(result["tms_source_files"].values()):
        return 1
    return 0


def command_prepare(args: argparse.Namespace) -> int:
    for config_path in args.config:
        config = load_dataset_config(config_path)
        output = prepare_tms_bundle(config, force=args.force)
        print(output)
    return 0


def command_train(args: argparse.Namespace) -> int:
    config = load_dataset_config(args.config)
    output = run_cross_validation(
        config,
        device_request=args.device,
        force=args.force,
        max_epochs_override=args.max_epochs,
        folds_override=_parse_int_list(args.folds),
        seeds_override=_parse_int_list(args.seeds),
        feature_mode=args.feature_mode,
    )
    print(output)
    return 0


def command_evaluate(args: argparse.Namespace) -> int:
    config = load_dataset_config(args.config)
    print(generate_training_report(config))
    return 0


def command_download_external(args: argparse.Namespace) -> int:
    keys = args.keys.split(",") if args.keys else None
    print(download_gse176206(args.config, keys=keys))
    return 0


def command_validate_external(args: argparse.Namespace) -> int:
    print(validate_gse176206(args.config, device=args.device, force=args.force))
    return 0


def command_sensitivity(args: argparse.Namespace) -> int:
    config = load_dataset_config(args.config)
    print(score_limb_droplet_sensitivity(config, args.model_dir, device=args.device))
    return 0


def _matrix_from_h5ad(path: Path) -> tuple[sp.csr_matrix, list[str], list[str]]:
    adata = ad.read_h5ad(path)
    matrix = adata.raw.X if adata.raw is not None else adata.X
    genes = adata.raw.var_names.astype(str).tolist() if adata.raw is not None else adata.var_names.astype(str).tolist()
    counts = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)
    return counts, genes, adata.obs_names.astype(str).tolist()


def command_score(args: argparse.Namespace) -> int:
    scorer = YouthScoreEnsemble(args.model_dir, device=args.device)
    input_path = Path(args.input)
    if input_path.is_dir():
        bundle = load_prepared_bundle(input_path)
        counts = bundle.counts
        genes = bundle.genes["gene_name"].astype(str).tolist()
        cells = bundle.cells["cell_id"].astype(str).tolist()
    elif input_path.suffix.lower() == ".h5ad":
        counts, genes, cells = _matrix_from_h5ad(input_path)
    else:
        raise ValueError("score input must be a prepared bundle directory or an h5ad file.")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    scores = scorer.score_counts(counts, genes, cells)
    if output.suffix.lower() == ".csv":
        scores.to_csv(output, index=False)
    else:
        scores.to_parquet(output, index=False)
    print(output)
    return 0


def command_run_all(args: argparse.Namespace) -> int:
    scat = load_dataset_config(DEFAULT_CONFIGS["scat"])
    limb = load_dataset_config(DEFAULT_CONFIGS["limb"])
    sensitivity = load_dataset_config(DEFAULT_CONFIGS["sensitivity"])
    for config in (scat, limb, sensitivity):
        prepare_tms_bundle(config, force=args.force)
    for config in (scat, limb):
        run_cross_validation(config, device_request=args.device, force=args.force)
        run_cross_validation(
            config,
            device_request=args.device,
            force=args.force,
            feature_mode="all",
        )
        generate_training_report(config)
    score_limb_droplet_sensitivity(sensitivity, limb.output_dir, device=args.device)
    if not args.skip_external:
        download_gse176206(DEFAULT_CONFIGS["external"])
        validate_gse176206(DEFAULT_CONFIGS["external"], device=args.device, force=args.force)
    generate_project_summary(ROOT)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="youth-score", description=__doc__)
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    commands = parser.add_subparsers(dest="command", required=True)

    doctor = commands.add_parser("doctor", help="Check the runtime, GPU, disk, and local TMS files.")
    doctor.set_defaults(function=command_doctor)

    prepare = commands.add_parser("prepare-tms", help="Create one or more compact TMS training bundles.")
    prepare.add_argument("--config", action="append", required=True)
    prepare.add_argument("--force", action="store_true")
    prepare.set_defaults(function=command_prepare)

    train = commands.add_parser("train", help="Run donor-aware cross-validation and model selection.")
    train.add_argument("--config", required=True)
    train.add_argument("--device", default="auto", choices=["auto", "cuda", "cpu"])
    train.add_argument("--force", action="store_true")
    train.add_argument("--max-epochs", type=int)
    train.add_argument("--folds", help="Comma-separated fold indices for a partial run.")
    train.add_argument("--seeds", help="Comma-separated seeds for a partial run.")
    train.add_argument(
        "--feature-mode",
        default="primary",
        choices=["primary", "all"],
        help="Use the primary exclusion mask or retain all genes for sensitivity analysis.",
    )
    train.set_defaults(function=command_train)

    evaluate = commands.add_parser("evaluate", help="Generate the English training report and plots.")
    evaluate.add_argument("--config", required=True)
    evaluate.set_defaults(function=command_evaluate)

    sensitivity = commands.add_parser("score-sensitivity", help="Score the Droplet Limb MSC sensitivity bundle.")
    sensitivity.add_argument("--config", default=str(DEFAULT_CONFIGS["sensitivity"]))
    sensitivity.add_argument("--model-dir", default=str(ROOT / "outputs" / "youth_score" / "limb_facs"))
    sensitivity.add_argument("--device", default="auto", choices=["auto", "cuda", "cpu"])
    sensitivity.set_defaults(function=command_sensitivity)

    download = commands.add_parser("download-external", help="Resume and verify GSE176206 downloads.")
    download.add_argument("--config", default=str(DEFAULT_CONFIGS["external"]))
    download.add_argument("--keys", help="Comma-separated dataset keys; defaults to both SOKM files.")
    download.set_defaults(function=command_download_external)

    validate = commands.add_parser("validate-external", help="Score and summarize GSE176206 SOKM data.")
    validate.add_argument("--config", default=str(DEFAULT_CONFIGS["external"]))
    validate.add_argument("--device", default="auto", choices=["auto", "cuda", "cpu"])
    validate.add_argument("--force", action="store_true")
    validate.set_defaults(function=command_validate_external)

    score = commands.add_parser("score", help="Score a prepared bundle or h5ad file with a selected ensemble.")
    score.add_argument("--model-dir", required=True)
    score.add_argument("--input", required=True)
    score.add_argument("--output", required=True)
    score.add_argument("--device", default="auto", choices=["auto", "cuda", "cpu"])
    score.set_defaults(function=command_score)

    run_all = commands.add_parser("run-all", help="Run preparation, training, reporting, sensitivity, and external validation.")
    run_all.add_argument("--device", default="auto", choices=["auto", "cuda", "cpu"])
    run_all.add_argument("--force", action="store_true")
    run_all.add_argument("--skip-external", action="store_true")
    run_all.set_defaults(function=command_run_all)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(asctime)s %(levelname)s %(message)s")
    try:
        return int(args.function(args))
    except Exception:
        logging.exception("Youth Score command failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())

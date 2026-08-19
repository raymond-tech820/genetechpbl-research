"""Command-line interface for Identity Score v1."""

from __future__ import annotations

import argparse

from .pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="identity-score-v1",
        description="Score and aggregate the frozen mouse MSC Identity Score v1.",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="Run cell scoring and animal aggregation.")
    run.add_argument("--config", required=True)
    run.add_argument("--input-h5ad")
    run.add_argument("--output-directory")
    run.add_argument("--force", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    output = run_pipeline(
        args.config,
        force=args.force,
        input_h5ad=args.input_h5ad,
        output_directory=args.output_directory,
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

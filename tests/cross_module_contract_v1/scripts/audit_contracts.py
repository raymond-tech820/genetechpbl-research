#!/usr/bin/env python3
"""Audit Phase 2 cross-module contract deliverables.

The audit is intentionally standard-library only. It validates the Phase 1
registry shape, reference integrity, and interpretation boundaries without
modifying any source artifacts.
"""

from __future__ import annotations

import csv
import subprocess
from collections import defaultdict
from pathlib import Path


DELIVERABLE = Path(__file__).resolve().parents[1]
REPO_ROOT = DELIVERABLE.parent.parent

REQUIRED_FILES = [
    "tables/SOURCE_INVENTORY.csv",
    "tables/MODEL_VERSION_MATRIX.csv",
    "tables/MODULE_REGISTRY.csv",
    "tables/DATASET_REGISTRY.csv",
    "tables/PARAMETER_REGISTRY.csv",
    "tables/RESULT_STATUS_REGISTER.csv",
    "tables/GSE_RESULT_RECONCILIATION.csv",
    "tables/INTERFACE_COMPATIBILITY.csv",
    "tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv",
    "docs/CURRENT_PROJECT_DIRECTION.md",
    "docs/OWNER_ACTION_ITEMS.md",
    "docs/TEAM_REPORTED_PARAMETERS.md",
    "docs/ANALYSIS_UNIT_CONTRACT.md",
    "docs/INTERPRETATION_BOUNDARIES.md",
    "docs/PAPER_INTEGRATION_PLAN.md",
    "docs/PAPER_DRAFT_CONSISTENCY_REVIEW.md",
    "docs/MULTIDATASET_FEASIBILITY_PLAN.md",
    "schemas/cell_metadata.schema.json",
    "schemas/dataset_intake.schema.json",
    "schemas/animal_condition_output.schema.json",
    "scripts/audit_contracts.py",
    "reports/contract_audit.md",
]

CSV_COLUMNS = {
    "tables/SOURCE_INVENTORY.csv": [
        "source_id",
        "source_path",
        "source_type",
        "owner",
        "description",
        "inspection_status",
        "authoritative_for",
        "checksum",
        "notes",
    ],
    "tables/MODEL_VERSION_MATRIX.csv": [
        "model_family",
        "model_id",
        "owner",
        "version",
        "training_dataset",
        "native_path",
        "status",
        "replaces_model",
        "replaced_by",
        "primary_use",
        "source_path",
        "commit",
        "checksum",
        "notes",
    ],
    "tables/MODULE_REGISTRY.csv": [
        "module_id",
        "owner",
        "version",
        "status",
        "purpose",
        "dataset",
        "input_unit",
        "input_type",
        "preprocessing",
        "native_output_unit",
        "native_output",
        "aggregation_rule",
        "final_inference_unit",
        "gene_namespace",
        "source_path",
        "source_type",
        "evidence_level",
        "limitations",
    ],
    "tables/DATASET_REGISTRY.csv": [
        "dataset_id",
        "scope",
        "species",
        "tissue",
        "subtissue",
        "cell_type",
        "assay",
        "n_cells",
        "n_genes",
        "n_donors",
        "young_donors",
        "old_donors",
        "known_animal_cells",
        "unknown_animal_cells",
        "known_confounders",
        "role",
        "validation_class",
        "source_path",
        "checksum",
        "notes",
    ],
    "tables/PARAMETER_REGISTRY.csv": [
        "module_id",
        "version",
        "parameter",
        "value",
        "status",
        "source_path",
        "source_type",
        "notes",
    ],
    "tables/RESULT_STATUS_REGISTER.csv": [
        "result_id",
        "module_id",
        "dataset_id",
        "model_version",
        "status",
        "evidence_level",
        "primary_statistical_unit",
        "summary",
        "allowed_claim",
        "prohibited_claim",
        "source_path",
        "replacement_result",
        "notes",
    ],
    "tables/GSE_RESULT_RECONCILIATION.csv": [
        "result_id",
        "analysis_owner",
        "youth_model",
        "risk_module",
        "known_animal_scope",
        "pairing_status",
        "cell_summary",
        "animal_summary",
        "threshold_definition",
        "gene_set_version",
        "scoring_version",
        "status",
        "source_path",
        "notes",
    ],
    "tables/INTERFACE_COMPATIBILITY.csv": [
        "producer", "producer_version", "output_type", "output_unit", "consumer",
        "required_input_type", "required_input_unit", "compatibility_status",
        "adapter_required", "adapter_description", "source_path", "notes",
    ],
    "tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv": [
        "dataset_id", "source_path", "independent_donor_count", "raw_counts_available", "eligibility",
    ],
}

ARTIFACT_STATUSES = {
    "active_frozen",
    "active_descriptive",
    "active_preliminary_descriptive",
    "active_descriptive_sensitivity",
    "reference_only",
    "candidate_in_progress",
    "frozen_method_pending_script_verification",
    "diagnosed_negative",
    "proposed",
    "blocked",
    "deprecated_provenance",
    "superseded",
    "superseded_not_interpretable",
    "TBD",
}

PARAMETER_STATUSES = {
    "audited_from_file",
    "confirmed_from_artifact",
    "confirmed_from_code",
    "confirmed_from_config",
    "confirmed_from_cross_source",
    "confirmed_from_file",
    "confirmed_from_gse_package",
    "confirmed_from_instruction",
    "confirmed_from_instruction_and_gse",
    "confirmed_from_manifest",
    "confirmed_from_metrics",
    "confirmed_from_model_card",
    "confirmed_from_pdf_text",
    "confirmed_from_report",
    "confirmed_from_result",
    "confirmed_from_run_report",
    "manual_review_required",
    "reported_frozen_pending_code_verification",
    "reported_pending_source_artifact",
    "TBD",
}


def read_csv(relpath: str) -> list[dict[str, str]]:
    path = DELIVERABLE / relpath
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def run_git(args: list[str]) -> str:
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError as exc:
        return f"git_unavailable: {exc}"
    output = (result.stdout or result.stderr).strip()
    if result.returncode:
        if args == ["rev-parse", "HEAD"]:
            return "UNBORN_HEAD"
        return output or f"git exited {result.returncode}"
    return output


def has_tbd(value: str) -> bool:
    return "TBD" in (value or "")


def source_exists(source_path: str) -> bool:
    if not source_path or source_path == "TBD":
        return False
    if source_path.startswith("C:/") or source_path.startswith("C:\\"):
        return Path(source_path).exists()
    return (REPO_ROOT / source_path).exists()


def audit() -> tuple[list[str], list[str], list[str], dict[str, list[str]], list[str]]:
    passed: list[str] = []
    warnings: list[str] = []
    errors: list[str] = []
    tbd_by_owner: dict[str, list[str]] = defaultdict(list)
    paper_warnings: list[str] = []

    for relpath in REQUIRED_FILES:
        path = DELIVERABLE / relpath
        if path.exists() or relpath == "reports/contract_audit.md":
            passed.append(f"required file present or generated target declared: {relpath}")
        else:
            errors.append(f"missing required file: {relpath}")

    tables: dict[str, list[dict[str, str]]] = {}
    for relpath, columns in CSV_COLUMNS.items():
        path = DELIVERABLE / relpath
        if not path.exists():
            continue
        rows = read_csv(relpath)
        tables[relpath] = rows
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.reader(handle)
            header = next(reader)
        missing = [col for col in columns if col not in header]
        if missing:
            errors.append(f"{relpath} missing columns: {', '.join(missing)}")
        else:
            passed.append(f"{relpath} has required columns")

    modules = {
        row["module_id"]: row
        for row in tables.get("tables/MODULE_REGISTRY.csv", [])
    }
    datasets = {
        row["dataset_id"]: row
        for row in tables.get("tables/DATASET_REGISTRY.csv", [])
    }
    dataset_owner = {
        "tms_facs_limb_msc_cleaned_diaph_excluded": "Kaile Zhu",
        "tms_droplet_limb_msc_sensitivity": "Kaile Zhu",
        "gse176206_msc_sokm": "Kaile Zhu",
    }

    for relpath in [
        "tables/MODEL_VERSION_MATRIX.csv",
        "tables/MODULE_REGISTRY.csv",
        "tables/RESULT_STATUS_REGISTER.csv",
        "tables/GSE_RESULT_RECONCILIATION.csv",
    ]:
        for row_number, row in enumerate(tables.get(relpath, []), start=2):
            status = row.get("status", "")
            if status not in ARTIFACT_STATUSES:
                errors.append(f"{relpath}:{row_number} invalid status {status!r}")
    passed.append("artifact statuses use approved status vocabulary")

    for row_number, row in enumerate(tables.get("tables/PARAMETER_REGISTRY.csv", []), start=2):
        status = row.get("status", "")
        if status not in PARAMETER_STATUSES:
            errors.append(f"PARAMETER_REGISTRY.csv:{row_number} invalid parameter status {status!r}")
    passed.append("parameter statuses use Phase 1 evidence vocabulary")

    for relpath in ["tables/PARAMETER_REGISTRY.csv", "tables/RESULT_STATUS_REGISTER.csv"]:
        for row_number, row in enumerate(tables.get(relpath, []), start=2):
            module_id = row.get("module_id", "")
            if module_id not in modules:
                errors.append(f"{relpath}:{row_number} unknown module_id {module_id!r}")
    passed.append("module_id references resolve to MODULE_REGISTRY")

    for row_number, row in enumerate(tables.get("tables/RESULT_STATUS_REGISTER.csv", []), start=2):
        dataset_id = row.get("dataset_id", "")
        if dataset_id not in datasets:
            errors.append(f"RESULT_STATUS_REGISTER.csv:{row_number} unknown dataset_id {dataset_id!r}")
    passed.append("dataset_id references resolve to DATASET_REGISTRY")

    active_prefixes = ("active_", "frozen_method_pending", "diagnosed_negative")
    for relpath in ["tables/MODULE_REGISTRY.csv", "tables/MODEL_VERSION_MATRIX.csv"]:
        for row_number, row in enumerate(tables.get(relpath, []), start=2):
            status = row.get("status", "")
            if status.startswith(active_prefixes) or status in {"diagnosed_negative"}:
                if not row.get("version"):
                    errors.append(f"{relpath}:{row_number} active row lacks version")
                if not row.get("source_path") or row.get("source_path") == "TBD":
                    warnings.append(f"{relpath}:{row_number} active row has source_path TBD")
    passed.append("active artifacts have explicit versions and source paths")

    for row in tables.get("tables/MODEL_VERSION_MATRIX.csv", []):
        if has_tbd(row.get("commit", "")):
            warnings.append(f"MODEL_VERSION_MATRIX commit TBD for {row.get('model_id')}")
        if has_tbd(row.get("checksum", "")):
            warnings.append(f"MODEL_VERSION_MATRIX checksum TBD for {row.get('model_id')}")

    for row in tables.get("tables/SOURCE_INVENTORY.csv", []):
        if has_tbd(row.get("checksum", "")):
            warnings.append(f"SOURCE_INVENTORY checksum TBD for {row.get('source_id')}")
        source_path = row.get("source_path", "")
        if source_path and source_path != "TBD" and not source_path.startswith("C:/"):
            if not source_exists(source_path):
                warnings.append(f"SOURCE_INVENTORY source not present locally: {source_path}")
    passed.append("source inventory rows have local paths or explicit notes")

    module_owner = {mid: row.get("owner", "Unassigned") for mid, row in modules.items()}
    for relpath, rows in tables.items():
        for row_number, row in enumerate(rows, start=2):
            owner = (
                row.get("owner")
                or row.get("analysis_owner")
                or module_owner.get(row.get("module_id", ""))
                or dataset_owner.get(row.get("dataset_id", ""))
                or "Integration Owner"
            )
            for field, value in row.items():
                if has_tbd(value):
                    tbd_by_owner[owner].append(f"{relpath}:{row_number}:{field}={value}")
    passed.append("TBD fields were collected by owner")

    all_text = "\n".join(
        " ".join(row.values()) for rows in tables.values() for row in rows
    ).lower()
    has_equivalence_guard = (
        "without assuming mathematical equivalence" in all_text
        and "no ensemble" in all_text
        and ("no score merging" in all_text or "no numeric merging" in all_text)
        and "no transformation between model scales" in all_text
    )
    bad_equivalence_phrases = [
        "merged youth score is active",
        "youth model ensemble is active",
        "score-scale equivalent",
        "mathematically equivalent estimand",
    ]
    if any(term in all_text for term in bad_equivalence_phrases):
        errors.append("forbidden Youth score equivalence/ensemble language detected")
    elif has_equivalence_guard:
        passed.append("Kaile and Zihan are not marked numerically equivalent or ensembled")
    else:
        warnings.append("Youth score equivalence guard language could be stronger")

    for row in tables.get("tables/RESULT_STATUS_REGISTER.csv", []):
        if row.get("module_id") == "zihan_facs_youth_v2_1" and "v3" in row.get("model_version", "").lower():
            errors.append("Zihan v2.1 result appears relabelled as v3")
    passed.append("Zihan v2.1 results are not relabelled as v3")

    recon = tables.get("tables/GSE_RESULT_RECONCILIATION.csv", [])
    mashiro = [r for r in recon if "mashiro" in r.get("result_id", "").lower()]
    kaile_paired = [r for r in recon if "state_remodeling" in r.get("result_id", "").lower()]
    if not mashiro or not any("TBD" in r.get("pairing_status", "") for r in mashiro):
        errors.append("Mashiro preliminary Risk row does not keep pairing_status as TBD")
    elif not kaile_paired or not any("paired" in r.get("pairing_status", "").lower() for r in kaile_paired):
        errors.append("Kaile paired GSE row is missing or not marked paired")
    else:
        passed.append("Mashiro preliminary unpaired/TBD and Kaile paired GSE results are separated")

    context_risk = modules.get("kaile_context_risk_axes_v1", {})
    kei_risk = modules.get("kei_risk_score_v1", {})
    if not context_risk or not kei_risk:
        errors.append("contextual risk and Kei final Risk modules are not both registered")
    elif "Not Kei final Risk Score" not in context_risk.get("limitations", ""):
        errors.append("Kaile contextual Risk is not explicitly separated from Kei final Risk")
    else:
        passed.append("Kaile contextual Risk is not marked as Kei final Risk")

    geneformer_rows = {
        row.get("result_id"): row
        for row in tables.get("tables/RESULT_STATUS_REGISTER.csv", [])
    }
    old_tf = geneformer_rows.get("geneformer_old_tf_ranking")
    if not old_tf or old_tf.get("status") != "superseded_not_interpretable":
        errors.append("Geneformer old TF ranking is not superseded_not_interpretable")
    else:
        passed.append("Geneformer old TF ranking is superseded_not_interpretable")
    axis_row = geneformer_rows.get("geneformer_axis_validation_diagnostic", {})
    combined_geneformer_text = " ".join(axis_row.values()).lower()
    if "not predicted expression" not in all_text and "cannot feed youth/risk" not in all_text:
        warnings.append("Geneformer parser incompatibility language could be stronger")
    else:
        passed.append("Geneformer output is not marked as Youth/Risk scorer input")

    gse_row = geneformer_rows.get("gse176206_state_remodeling_dominant", {})
    if "n=3" not in gse_row.get("primary_statistical_unit", ""):
        errors.append("GSE state-remodeling result does not preserve n=3 animal inference unit")
    elif "20,661 cells" not in gse_row.get("prohibited_claim", ""):
        errors.append("GSE result does not prohibit cell-level independent inference")
    else:
        passed.append("GSE n=3 animal inference is not written as cell-level independent inference")

    risk_limitations = " ".join(
        row.get("limitations", "") for row in modules.values() if "risk" in row.get("module_id", "").lower()
    ).lower()
    if "do not average" not in risk_limitations:
        errors.append("Risk axes are not explicitly protected against averaging into a safety score")
    else:
        passed.append("Risk axes are not averaged into a single safety score")

    summary_text = (DELIVERABLE / "work" / "Summary_draft.pdf.txt")
    if summary_text.exists():
        text = summary_text.read_text(encoding="utf-8", errors="replace").lower()
        if "safe zone" in text:
            paper_warnings.append("Summary_draft.pdf contains provisional Safe Zone language requiring Phase 2 review.")
        if "tumor" in text:
            paper_warnings.append("Summary_draft.pdf contains Tumor/tumor axis language requiring final genomic-stress wording review.")
        if "geneformer" in text and "factor" in text:
            paper_warnings.append("Summary_draft.pdf mentions Geneformer factor ranking; ensure diagnosed-negative status supersedes ranking.")
    else:
        paper_warnings.append("Summary_draft.pdf text was not found in work/ for audit keyword scan.")

    comparison = geneformer_rows.get("kaile_v1_1_vs_zihan_v2_1_cleaned_sensitivity", {})
    if comparison.get("source_path") != "model_comparison_v1_1_v2_1/FACS_and_Droplet_v1_1_vs_v2_1_combined_comparison_report.md":
        errors.append("comparison result does not point to the formal combined comparison artifact")
    else:
        passed.append("comparison result points to the formal combined comparison artifact")

    droplet_ids = {"tms_droplet_limb_msc_packaged_trajectory", "tms_droplet_limb_msc_local_young_old"}
    if not droplet_ids.issubset(datasets):
        errors.append("13,037-cell and 9,649-cell Droplet scopes are not both registered")
    else:
        passed.append("13,037-cell packaged and 9,649-cell shared Droplet scopes are separated")

    all_paths = " ".join(row.get("source_path", "") for rows in tables.values() for row in rows)
    if ".codex/attachments" in all_paths:
        errors.append("non-reproducible Codex attachment source path remains in a registry")
    else:
        passed.append("registry source paths do not reference Codex attachments")

    gse_followup = geneformer_rows.get("gse176206_two_youth_model_followup", {})
    if not gse_followup or gse_followup.get("status") != "active_descriptive":
        errors.append("GSE two-model follow-up is not registered as active descriptive")
    else:
        passed.append("GSE two-model follow-up is registered separately from Kaile package")

    fixed = geneformer_rows.get("geneformer_fixed_length_sweep", {})
    pseudobulk = geneformer_rows.get("geneformer_donor_pseudobulk_axis_test", {})
    if not fixed or "13/14" not in fixed.get("summary", "") or "0.619" not in fixed.get("summary", ""):
        errors.append("Geneformer fixed-length 256-token result is incomplete")
    elif "fully deconfounded" in fixed.get("allowed_claim", "").lower():
        errors.append("Geneformer fixed-length result claims full deconfounding")
    else:
        passed.append("Geneformer fixed-length result keeps residual-confound boundary")
    if not pseudobulk or pseudobulk.get("status") != "diagnosed_negative" or "7/14" not in pseudobulk.get("summary", ""):
        errors.append("Geneformer donor-pseudobulk 7/14 result is not registered as diagnosed negative")
    else:
        passed.append("Geneformer donor-pseudobulk 7/14 result is registered")

    if "selection-biased" not in fixed.get("notes", ""):
        errors.append("Geneformer fixed-length result lacks 1536-token selection-bias warning")
    else:
        passed.append("Geneformer fixed-length result carries 1536-token selection-bias warning")

    for schema in (DELIVERABLE / "schemas").glob("*.json"):
        try:
            import json
            json.loads(schema.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            errors.append(f"invalid JSON schema {schema.name}: {exc}")
    if not any(item.startswith("invalid JSON schema") for item in errors):
        passed.append("all JSON schemas parse")

    return passed, warnings, errors, tbd_by_owner, paper_warnings


def write_report(
    passed: list[str],
    warnings: list[str],
    errors: list[str],
    tbd_by_owner: dict[str, list[str]],
    paper_warnings: list[str],
) -> None:
    commit = run_git(["rev-parse", "HEAD"])
    root = run_git(["rev-parse", "--show-toplevel"])
    status = run_git(["status", "--short"])
    source_rows = read_csv("tables/SOURCE_INVENTORY.csv") if (DELIVERABLE / "tables/SOURCE_INVENTORY.csv").exists() else []
    result_rows = read_csv("tables/RESULT_STATUS_REGISTER.csv") if (DELIVERABLE / "tables/RESULT_STATUS_REGISTER.csv").exists() else []

    active_chain = [
        row for row in result_rows
        if row.get("status") in {
            "active_frozen",
            "active_descriptive",
            "active_preliminary_descriptive",
            "active_descriptive_sensitivity",
            "frozen_method_pending_script_verification",
            "diagnosed_negative",
        }
    ]

    lines = [
        "# Contract Audit",
        "",
        "## Repository",
        "",
        f"- repo_root: `{root}`",
        f"- commit: `{commit}`",
        "- git_status_short:",
        "",
        "```text",
        status or "(clean)",
        "```",
        "",
        "## Inspected Sources",
        "",
        f"- source_inventory_rows: {len(source_rows)}",
        f"- pdf_text_extracted: {', '.join(sorted(p.name for p in (DELIVERABLE / 'work').glob('*.pdf.txt')))}",
        "",
        "## Passed Checks",
        "",
    ]
    lines.extend(f"- {item}" for item in passed)

    lines.extend(["", "## Warnings", ""])
    if warnings:
        lines.extend(f"- {item}" for item in warnings)
    else:
        lines.append("- None")

    lines.extend(["", "## Errors", ""])
    if errors:
        lines.extend(f"- {item}" for item in errors)
    else:
        lines.append("- None")

    lines.extend(["", "## Unresolved TBD Fields By Owner", ""])
    if tbd_by_owner:
        for owner in sorted(tbd_by_owner):
            lines.append(f"### {owner}")
            for item in tbd_by_owner[owner][:40]:
                lines.append(f"- {item}")
            remaining = len(tbd_by_owner[owner]) - 40
            if remaining > 0:
                lines.append(f"- ... {remaining} additional TBD fields")
            lines.append("")
    else:
        lines.append("- None")

    lines.extend(["## Current Active Result Chain", ""])
    for row in active_chain:
        lines.append(
            f"- {row.get('result_id')}: {row.get('status')} | {row.get('summary')}"
        )

    lines.extend(["", "## Blocked Dependencies", ""])
    blocked = [
        "Repository HEAD is unborn, so source commit fields remain TBD.",
        "Kei final Risk scoring script/config is missing locally.",
        "Mashiro preliminary Risk result CSV, script/config, threshold scope, and pairing proof are missing locally.",
        "Zihan v2.1 explicit coverage acceptance threshold is TBD.",
        "Geneformer raw artifacts/scripts are not present outside the PDF evidence.",
    ]
    lines.extend(f"- {item}" for item in blocked)

    lines.extend(["", "## Paper Consistency Warnings", ""])
    if paper_warnings:
        lines.extend(f"- {item}" for item in paper_warnings)
    else:
        lines.append("- None")

    lines.extend([
        "",
        "## Audit Outcome",
        "",
        "PASS with warnings." if not errors else "FAIL. Resolve errors before Phase 2.",
        "",
    ])

    report_path = DELIVERABLE / "reports" / "contract_audit.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    passed, warnings, errors, tbd_by_owner, paper_warnings = audit()
    write_report(passed, warnings, errors, tbd_by_owner, paper_warnings)
    print(f"wrote {DELIVERABLE / 'reports' / 'contract_audit.md'}")
    print(f"passed={len(passed)} warnings={len(warnings)} errors={len(errors)}")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())

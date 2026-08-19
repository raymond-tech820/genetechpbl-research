# Cross-Module Contract and Dataset-Aware Intake Framework

This deliverable is a reproducibility and compatibility layer, not a new biological model, model retraining, score merger, or multi-dataset harmonization result.

> Align module outputs without assuming mathematical equivalence.

## Purpose

The project combines several analyses that use different inputs, scoring units, and inference units. This contract records those boundaries so that results can be interpreted together at the animal/donor level without treating different outputs as one numerical scale.

In particular, it distinguishes cell-level measurement from donor/animal-level inference, preserves each frozen scorer's native preprocessing path, and makes missing evidence visible as `TBD` rather than inferring values.

## What this release records

- Kaile v1.1 and Zihan v2.1 as independent active frozen Youth Score models.
- The official frozen-model sensitivity comparison, with formal provenance and no overall winner claim.
- Separate 16-donor packaged and 12-donor shared-comparison TMS Droplet scopes.
- GSE176206 state-remodeling evidence and the separate two-Youth-model biological follow-up.
- The distinction between Kaile contextual Risk modules, Kei's pending final Risk method, and Mashiro's preliminary visual result.
- Module interfaces, analysis units, interpretation boundaries, paper-integration guidance, and a paper-draft consistency review.
- Future multi-dataset intake criteria and a standard animal-condition output schema.

## Geneformer Scope

The currently registered Geneformer arm remains a methodological diagnostic and is not an input to the Youth or Risk scorers. New Geneformer progress is intentionally outside the scope of this release. Its updated evaluation criteria and any new supporting artifacts will be added in a later contract update after they are documented and reviewed.

## Directory Guide

- `tables/`: source, model, dataset, parameter, result, and interface registries.
- `docs/`: contracts, interpretation limits, owner action items, and manuscript planning.
- `schemas/`: future dataset-intake, cell-metadata, and animal-condition output specifications.
- `scripts/audit_contracts.py`: reproducibility and consistency checks.
- `reports/contract_audit.md`: latest audit outcome and unresolved dependencies.

## Important Boundaries

- Kaile and Zihan raw score magnitudes are not comparable and must not be merged or ensembled.
- Droplet transfer is cross-assay sensitivity, not independent external validation.
- GSE176206 uses animals/donors as the inference layer; cells are measurement units.
- Kaile contextual Risk is not Kei final Risk.
- Current work does not establish a clinical safety result, accelerated aging result, or validated Safe Zone.

## Verification

Run:

```text
python deliverables/cross_module_contract_v1/scripts/audit_contracts.py
```

The audit parses all registered CSV and JSON-schema files and checks the principal provenance and interpretation boundaries.

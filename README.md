# Genentech PBL - AI for Regenerative Biology (Stem Cells and Repair)

This repository contains model packages, historical implementations,
sensitivity tests, cross-module audits, and manuscript materials developed for
the Genentech AI for Regenerative Biology PBL project.

## Repository Structure

```text
models/       Current or maintained model packages
archived/     Superseded models retained for provenance
tests/        Sensitivity analyses, comparisons, and compatibility audits
reports/      Manuscript contribution packages
```

Assay and dataset names are intentionally omitted from top-level Youth Score
directory names. Complete training-data and version provenance remains inside
each package.

## Models

| Package | Contributor | Description |
|---|---|---|
| [`models/youth_score_v1_1`](./models/youth_score_v1_1/) | Kaile Zhu | Current cleaned-limb Youth Score package |
| [`models/youth_score_v2_1`](./models/youth_score_v2_1/) | Zihan Zhou | Maintained historical Youth Score baseline |
| [`models/youth_score_v3_1`](./models/youth_score_v3_1/) | Zihan Zhou | Current frozen Youth Score framework with M1 primary and M4 sensitivity models |
| [`models/geneformer`](./models/geneformer/) | Jia Qi Choy | Zero-shot Geneformer perturbation and embedding analysis |
| [`models/identity`](./models/identity/) | Kaile Zhu | Frozen knowledge-driven MSC Identity Score |
| [`models/risk`](./models/risk/) | Kei Hasegawa | Risk Score module |

See [`models/README.md`](./models/README.md) for model roles.

## Archived Models

Superseded packages are retained under [`archived/`](./archived/):

- [`archived/risk_score_archived`](./archived/risk_score_archived/)
- [`archived/youth_score_v1`](./archived/youth_score_v1/)
- [`archived/youth_score_v2`](./archived/youth_score_v2/)

Archived packages support provenance and historical reproducibility. They are
not recommended as current operational entry points.

## Tests and Cross-Module Audits

Evaluation resources are under [`tests/`](./tests/):

- [`tests/youth_score_v1_sensitivity`](./tests/youth_score_v1_sensitivity/)
- [`tests/youth_score_v1_1_vs_v2_1_comparison`](./tests/youth_score_v1_1_vs_v2_1_comparison/)
- [`tests/cross_module_contract_v1`](./tests/cross_module_contract_v1/)

Run the compatibility audit from the repository root:

```bash
python tests/cross_module_contract_v1/scripts/audit_contracts.py
```

## Manuscript Reports

Manuscript contribution packages and role allocation are documented in
[`reports/`](./reports/) and [`reports/README.md`](./reports/README.md).

## Interpretation Boundaries

- Cells are measurement units; donor or animal identifiers define the formal
  inference unit unless a package explicitly documents another design.
- Youth Score implementations may use different scoring and aggregation
  contracts. Their raw magnitudes must not be merged or treated as equivalent.
- Geneformer embedding shifts are not expression matrices and cannot be passed
  directly into Youth, Identity, or Risk parsers.
- Youth, Identity, and Risk represent complementary state axes.
- Internal validation, cross-assay transportability, and external biological
  application answer different questions and must be reported separately.
- No package in this repository establishes a universal aging clock, clinical
  safety score, or validated rejuvenation endpoint.

## Getting Started

1. Select the required package under `models/`.
2. Read its README and model card.
3. Confirm input unit, gene namespace, normalization, and coverage requirements.
4. Use only the frozen parser and model artifacts supplied by that package.
5. Report the complete repository path and internal model version.

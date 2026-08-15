# genetechpbl-research

# Genentech PBL Youth Score Models

This repository temporarily contains MSC Youth Score model deliverables and model-comparison outputs developed for the Genentech AI for Regenerative Biology PBL project.

| Model                                                                                 | Contributor | Description                                                                                                                                                                                                                                                                 |
| ------------------------------------------------------------------------------------- | ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`droplet_msc_youth_score_v1`](./droplet_msc_youth_score_v1/)                         | Zihan Zhou  | Droplet-trained model retained as a sensitivity/reference model rather than the current primary FACS model.                                                                                                                                                                 |
| [`facs_msc_youth_score_v1`](./facs_msc_youth_score_v1/)                               | Kaile Zhu   | Deprecated FACS v1 model, retained for provenance only. Use `facs_msc_youth_score_v1_1_cleaned_limb` for the current Kaile model.                                                                                                                                           |
| [`facs_msc_youth_score_v1_1_cleaned_limb`](./facs_msc_youth_score_v1_1_cleaned_limb/) | Kaile Zhu   | Current Kaile FACS model for comparison and further refinement, retrained after explicit exclusion of age-confounded diaphragm cells.                                                                                                                                       |
| [`facs_msc_youth_score_v2`](./facs_msc_youth_score_v2/)                               | Zihan Zhou  | Deprecated FACS v2 model, retained for provenance only.                                                                                                                                                                                                                     |
| [`facs_msc_youth_score_v2_1`](./facs_msc_youth_score_v2_1/)                           | Zihan Zhou  | Superseded FACS model and historical diaphragm-excluded baseline. Retained for provenance and reproducibility of the v1.1-v2.1 comparison; use `facs_msc_youth_score_v3_1` for the current Zihan model.                                                                     |
| [`facs_msc_youth_score_v3_1`](./facs_msc_youth_score_v3_1/)                           | Zihan Zhou  | Current frozen FACS model package. M1 is the primary deployable reference; post-ablation M4 is retained as a high-stringency sensitivity comparator. Includes training, formal 999-permutation, GSE176206 biological-application, and TMS Droplet transportability reports. |
| [`msc_identity_score_v1`](./msc_identity_score_v1/)                                   | Kaile Zhu   | Frozen knowledge-driven mouse MSC Identity Score v1 used for the 2026-08-06 GSE176206 handoff. Provides deterministic cell-level scoring, animal-within-exact-arm aggregation, QC, reference results, and reproducibility tests.                                                |
| [`geneformer_perturbation`](./geneformer_perturbation/)                               | Jia Qi Choy | Zero-shot Geneformer in-silico perturbation arm. Diagnoses that the Old-to-Young embedding axis is confounded with sequencing depth and not interpretable on this cohort; see folder for full methods log.                                                                  |
| [cross_module_contract_v1](./cross_module_contract_v1/)                           | Bowen Liu   | Cross-module reproducibility and compatibility layer. Defines native scoring paths, donor/animal inference units, model and dataset provenance, interface compatibility, interpretation boundaries, future multi-dataset intake, and standardized animal-condition outputs. |

## Cross-Module Integration Contract

The `cross_module_contract_v1/` directory documents the integration work that connects the Youth Score, Identity, Risk, and Geneformer modules without forcing them into one mathematical score. It is a reproducibility and compatibility layer, not a new biological model.

The contract records:

- native input and scoring paths for each module;
- the distinction between cell-level measurement and donor/animal-level inference;
- model, dataset, result, parameter, and source provenance;
- compatibility status and adapter requirements between module outputs;
- interpretation boundaries, unresolved `TBD` items, and owner action items;
- future multi-dataset intake criteria; and
- standardized animal-condition output fields for downstream integration.

Important boundaries:

- Kaile and Zihan Youth Score raw magnitudes are not mathematically equivalent and must not be merged or ensembled;
- current integration is performed through native outputs and explicit animal-condition summaries;
- Geneformer embedding shifts are not expression matrices and cannot be passed directly into Youth or Risk parsers; and
- the contract does not establish a universal aging clock, clinical safety score, or validated Safe Zone.

### Running the Contract Audit

From the repository root, run:

```bash
python cross_module_contract_v1/scripts/audit_contracts.py
```

The audit checks registry references, source provenance, CSV parsing, JSON schemas, model-role boundaries, and forbidden score-equivalence claims. The latest audit report is in:

```text
cross_module_contract_v1/reports/contract_audit.md
```

## MSC Identity Score

The frozen Identity module is available in:

```text
msc_identity_score_v1/
```

Identity Score v1 is a knowledge-driven expression-module score rather than a
trained classifier. It measures how strongly a cell expresses the frozen mouse
MSC identity program relative to expression-matched background genes. It is a
contextual companion to the Youth Score and must not be interpreted as a
clinical safety, potency, tumour-risk, or causal rejuvenation score.

The scorer operates on non-negative raw counts, applies library-size
normalization and `log1p`, and calculates a deterministic matched-background
score for each cell. Formal GSE176206 integration uses the median score within:

```text
age_group x exact_treatment_arm x animal_label
```

Animal labels are nested within age and exact treatment arm. Reused labels
across treatment arms are not global subject identifiers, so cross-arm
comparisons are unpaired. See the package README, frozen methods, statistical
unit notice, reference results, and tests for the complete scoring contract.

## Current Model Status

Zihan Zhou's current frozen model is:

```text
facs_msc_youth_score_v3_1/
```

Start with:

```text
facs_msc_youth_score_v3_1/README.md
```

Within v3.1, M1 is the operational primary model and M4 is a post-ablation exploratory sensitivity comparator. The package contains minimal frozen scoring artifacts plus the completed training and cross-assay evidence. It supports a FACS-relative young-like MSC transcriptional state score, not a universal aging clock or proof of technical independence.

Kaile Zhu's current cleaned-limb model remains:

```text
facs_msc_youth_score_v1_1_cleaned_limb/
```

These current models may be compared in future work using a separately frozen protocol. The existing comparison folder:

```text
model_comparison_v1_1_v2_1/
```

is a historical comparison of Kaile v1.1 with the superseded Zihan v2.1 model. It includes:

- FACS same-input frozen scorer application;
- Droplet cross-assay sensitivity;
- internal validation summary separating Kaile packaged donor OOF CV from Zihan nested LOMO;
- paired donor bootstrap uncertainty analysis;
- model concordance decomposition within and after age adjustment;
- technical-variable audits;
- FACS gene-signature overlap analysis.

Important interpretation boundary: the historical comparison supports directionally consistent young-old state separation for both evaluated models, but it does not prove that either model is technically independent or externally validated. Paired donor bootstrap does not establish statistically supported overall superiority of Kaile v1.1 over Zihan v2.1, although Kaile v1.1 shows more consistent descriptive age ordering.

The older FACS v1, v2, and v2.1 folders are preserved to keep prior analyses reproducible. See each model directory and comparison directory for detailed documentation, scoring files, validation results, and limitations.

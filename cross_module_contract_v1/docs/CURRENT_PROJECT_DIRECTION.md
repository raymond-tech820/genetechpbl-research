# Current Project Direction

This Phase 1 deliverable is a reproducibility and compatibility layer, not a new biological model.

Core principle:

> Align module outputs without assuming mathematical equivalence.

## Active Direction

The current paper direction should treat the project as a dataset-aware, multi-axis evaluation of partial reprogramming in MSCs. The main GSE176206 interpretation is `state_remodeling_dominant`: SOKM is associated with lower Kaile Youth Score, lower MSC identity, lower ECM/collagen programs, and redistribution away from baseline/Col11a1-like states toward reprogramming and reprogramming-cycling states.

This should not be rewritten as accelerated aging, clinical harm, tumourigenicity, therapeutic efficacy, or validated safety.

## Youth Models

There are two active frozen Youth models in scope:

- `kaile_facs_youth_v1_1_cleaned_limb`: raw single-cell counts are scored at cell level and summarized to donor/animal level. The confirmed native donor aggregation is arithmetic mean of cell `youth_score`.
- `zihan_facs_youth_v2_1`: raw single-cell counts are first aggregated within donor/sample to pseudobulk counts, then scored at donor/sample level.

The two native estimands are different:

`aggregate(score(cell)) != score(aggregate(expression))`

Allowed comparisons are direction, ranking, age association, SOKM response direction, standardized within-model effects, and technical sensitivity. The contract does not permit numeric score-scale equivalence, a Youth model ensemble, or merged Youth Score outputs.

The official Kaile v1.1 versus Zihan v2.1 sensitivity result is sourced from `model_comparison_v1_1_v2_1/FACS_and_Droplet_v1_1_vs_v2_1_combined_comparison_report.md`. It supports concordant Young-Old direction in cleaned FACS and the shared local Droplet subset, without proving an overall winner, technical independence, external validation, or a continuous aging clock. The two Droplet analyses use different donor scopes and must not be reported as one dataset size: Kaile's packaged trajectory has 13,037 cells and 16 donors, whereas the shared comparison subset has 9,649 cells and 12 donors.

## GSE176206

The current active GSE package is:

`GSE176206_MSC_Reprogramming_Evaluation_v1_revised_20260730`

The revised package uses known animals as the primary inference unit and treats cells as measurement units. Unknown-animal cells are descriptive sensitivity only. The package records 20,661 cells and 28,694 genes, with 18,686 known-animal cells and 1,975 unknown-animal cells.

The active interpretation is `state_remodeling_dominant`, not stable rejuvenation and not safety/harm evidence.

The standalone two-model GSE follow-up is a separate artifact from the Kaile state-remodeling package. Both frozen models show negative SOKM-minus-control direction in valid young and aged animals; their raw score magnitudes are not comparable. The aggregate effect includes state redistribution and selected within-state remodeling, not uniform accelerated aging.

## Risk Score

`kei_risk_score_v1` is registered as `frozen_method_pending_script_verification`.

The local master gene file is `risk_genes_housekeeping_revised.csv`. It has source columns `Dedifferentiation`, `Inflam`, and `Tumor`. The contract labels the three axes as identity-loss-associated, SASP-associated, and genomic-stress-associated, but the source header-to-final-label mapping still needs owner confirmation.

`mashiro_gse176206_kei_risk_rerun` remains preliminary descriptive. No exact result CSV, script/config, pairing proof, or threshold scope was found locally.

## Geneformer

`jia_geneformer_v1_10m_axis_diagnostic` is `diagnosed_negative`.

The Geneformer output is a scalar cosine-similarity shift toward a young centroid, not a predicted expression matrix. It cannot be sent into the Youth or Risk scorers. The old seven-factor perturbation ranking is superseded and not interpretable under the failed axis-validation diagnostic.

Fixed token length is a confirmed contributor to Geneformer depth sensitivity but not a complete solution: the reported 256-token same-cell experiment improves LOO to 13/14 with residual confound r=0.619, while donor pseudobulk separation is 7/14. The axis remains unsuitable for perturbation ranking.

## Phase Boundary

Phase 1 completed only the current registries, owner actions, and audit contract. Phase 2 items were intentionally not created:

- `INTERFACE_COMPATIBILITY.csv`
- schemas
- paper integration plan
- paper draft consistency review
- multidataset feasibility plan

All `TBD` fields are explicit owner-confirmation points, not silent failures.

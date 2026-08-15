# Contract Audit

## Repository

- repo_root: `C:/Users/LiuBW/Documents/PBL 2`
- commit: `UNBORN_HEAD`
- git_status_short:

```text
?? GSE176206_MSC_Reprogramming_Evaluation_v1_revised_20260730/
?? GSE176206_youth_score_biology_followup_report.md
?? Summary_draft.pdf
?? "animal_risk_summaries mashiro.png"
?? axis_validation_update.md
?? deliverables/
?? facs_msc_youth_score_v1_1_cleaned_limb/
?? facs_msc_youth_score_v2_1/
?? feedback.md
?? genetech_pbl_paper.pdf
?? model_comparison_v1_1_v2_1/
?? risk_genes_housekeeping_revised.csv
```

## Inspected Sources

- source_inventory_rows: 43
- pdf_text_extracted: GSE176206_MSC_Reprogramming_Evaluation_Report_EN.pdf.txt, Summary_draft.pdf.txt, genetech_pbl_paper.pdf.txt

## Passed Checks

- required file present or generated target declared: tables/SOURCE_INVENTORY.csv
- required file present or generated target declared: tables/MODEL_VERSION_MATRIX.csv
- required file present or generated target declared: tables/MODULE_REGISTRY.csv
- required file present or generated target declared: tables/DATASET_REGISTRY.csv
- required file present or generated target declared: tables/PARAMETER_REGISTRY.csv
- required file present or generated target declared: tables/RESULT_STATUS_REGISTER.csv
- required file present or generated target declared: tables/GSE_RESULT_RECONCILIATION.csv
- required file present or generated target declared: tables/INTERFACE_COMPATIBILITY.csv
- required file present or generated target declared: tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv
- required file present or generated target declared: docs/CURRENT_PROJECT_DIRECTION.md
- required file present or generated target declared: docs/OWNER_ACTION_ITEMS.md
- required file present or generated target declared: docs/TEAM_REPORTED_PARAMETERS.md
- required file present or generated target declared: docs/ANALYSIS_UNIT_CONTRACT.md
- required file present or generated target declared: docs/INTERPRETATION_BOUNDARIES.md
- required file present or generated target declared: docs/PAPER_INTEGRATION_PLAN.md
- required file present or generated target declared: docs/PAPER_DRAFT_CONSISTENCY_REVIEW.md
- required file present or generated target declared: docs/MULTIDATASET_FEASIBILITY_PLAN.md
- required file present or generated target declared: schemas/cell_metadata.schema.json
- required file present or generated target declared: schemas/dataset_intake.schema.json
- required file present or generated target declared: schemas/animal_condition_output.schema.json
- required file present or generated target declared: scripts/audit_contracts.py
- required file present or generated target declared: reports/contract_audit.md
- tables/SOURCE_INVENTORY.csv has required columns
- tables/MODEL_VERSION_MATRIX.csv has required columns
- tables/MODULE_REGISTRY.csv has required columns
- tables/DATASET_REGISTRY.csv has required columns
- tables/PARAMETER_REGISTRY.csv has required columns
- tables/RESULT_STATUS_REGISTER.csv has required columns
- tables/GSE_RESULT_RECONCILIATION.csv has required columns
- tables/INTERFACE_COMPATIBILITY.csv has required columns
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv has required columns
- artifact statuses use approved status vocabulary
- parameter statuses use Phase 1 evidence vocabulary
- module_id references resolve to MODULE_REGISTRY
- dataset_id references resolve to DATASET_REGISTRY
- active artifacts have explicit versions and source paths
- source inventory rows have local paths or explicit notes
- TBD fields were collected by owner
- Kaile and Zihan are not marked numerically equivalent or ensembled
- Zihan v2.1 results are not relabelled as v3
- Mashiro preliminary unpaired/TBD and Kaile paired GSE results are separated
- Kaile contextual Risk is not marked as Kei final Risk
- Geneformer old TF ranking is superseded_not_interpretable
- Geneformer output is not marked as Youth/Risk scorer input
- GSE n=3 animal inference is not written as cell-level independent inference
- Risk axes are not averaged into a single safety score
- comparison result points to the formal combined comparison artifact
- 13,037-cell packaged and 9,649-cell shared Droplet scopes are separated
- registry source paths do not reference Codex attachments
- GSE two-model follow-up is registered separately from Kaile package
- Geneformer fixed-length result keeps residual-confound boundary
- Geneformer donor-pseudobulk 7/14 result is registered
- Geneformer fixed-length result carries 1536-token selection-bias warning
- all JSON schemas parse

## Warnings

- MODEL_VERSION_MATRIX commit TBD for kaile_facs_youth_v1_1_cleaned_limb
- MODEL_VERSION_MATRIX commit TBD for zihan_facs_youth_v2_1
- MODEL_VERSION_MATRIX commit TBD for kaile_gse176206_state_remodeling_v1
- MODEL_VERSION_MATRIX commit TBD for kei_risk_score_v1
- MODEL_VERSION_MATRIX commit TBD for mashiro_gse176206_kei_risk_rerun
- MODEL_VERSION_MATRIX commit TBD for jia_geneformer_v1_10m_axis_diagnostic
- SOURCE_INVENTORY checksum TBD for geneformer_axis_length_pdf
- SOURCE_INVENTORY source not present locally: axis_length_report.pdf

## Errors

- None

## Unresolved TBD Fields By Owner

### Hasegawa Kei
- tables/MODEL_VERSION_MATRIX.csv:5:replaces_model=TBD
- tables/MODEL_VERSION_MATRIX.csv:5:replaced_by=TBD
- tables/MODEL_VERSION_MATRIX.csv:5:commit=TBD
- tables/MODEL_VERSION_MATRIX.csv:5:notes=Team-reported settings are registered but scoring script path is TBD.
- tables/MODULE_REGISTRY.csv:10:aggregation_rule=Aggregation and threshold scope are TBD pending script/config.
- tables/MODULE_REGISTRY.csv:10:final_inference_unit=animal-condition TBD
- tables/PARAMETER_REGISTRY.csv:67:notes=Exact matched-control scope is TBD.
- tables/PARAMETER_REGISTRY.csv:68:value=TBD
- tables/PARAMETER_REGISTRY.csv:68:status=TBD
- tables/PARAMETER_REGISTRY.csv:68:source_path=TBD
- tables/PARAMETER_REGISTRY.csv:68:source_type=TBD
- tables/PARAMETER_REGISTRY.csv:69:value=TBD
- tables/PARAMETER_REGISTRY.csv:69:status=TBD
- tables/PARAMETER_REGISTRY.csv:69:source_path=TBD
- tables/PARAMETER_REGISTRY.csv:69:source_type=TBD
- tables/PARAMETER_REGISTRY.csv:70:value=TBD
- tables/PARAMETER_REGISTRY.csv:70:status=TBD
- tables/PARAMETER_REGISTRY.csv:70:source_path=TBD
- tables/PARAMETER_REGISTRY.csv:70:source_type=TBD
- tables/PARAMETER_REGISTRY.csv:71:value=TBD
- tables/PARAMETER_REGISTRY.csv:71:status=TBD
- tables/PARAMETER_REGISTRY.csv:71:source_path=TBD
- tables/PARAMETER_REGISTRY.csv:71:source_type=TBD
- tables/PARAMETER_REGISTRY.csv:72:value=TBD
- tables/PARAMETER_REGISTRY.csv:72:status=TBD
- tables/PARAMETER_REGISTRY.csv:72:source_path=TBD
- tables/PARAMETER_REGISTRY.csv:72:source_type=TBD
- tables/RESULT_STATUS_REGISTER.csv:9:primary_statistical_unit=cell; animal-condition TBD
- tables/RESULT_STATUS_REGISTER.csv:9:summary=Kei parameters are registered from team instruction and master gene list was audited; scoring script/config and threshold scope are TBD.
- tables/RESULT_STATUS_REGISTER.csv:9:replacement_result=TBD

### Integration Owner
- tables/DATASET_REGISTRY.csv:3:subtissue=TBD
- tables/DATASET_REGISTRY.csv:3:n_genes=TBD
- tables/DATASET_REGISTRY.csv:4:subtissue=TBD
- tables/INTERFACE_COMPATIBILITY.csv:4:source_path=TBD
- tables/INTERFACE_COMPATIBILITY.csv:4:notes=Threshold scope and final script/config remain TBD.
- tables/INTERFACE_COMPATIBILITY.csv:11:producer_version=TBD

### Jia Qi Choy
- tables/SOURCE_INVENTORY.csv:43:checksum=TBD
- tables/MODEL_VERSION_MATRIX.csv:7:replaces_model=TBD
- tables/MODEL_VERSION_MATRIX.csv:7:replaced_by=TBD
- tables/MODEL_VERSION_MATRIX.csv:7:commit=TBD
- tables/RESULT_STATUS_REGISTER.csv:11:replacement_result=TBD
- tables/RESULT_STATUS_REGISTER.csv:14:replacement_result=TBD
- tables/RESULT_STATUS_REGISTER.csv:15:replacement_result=TBD

### Kaile Zhu
- tables/MODEL_VERSION_MATRIX.csv:2:replaced_by=TBD
- tables/MODEL_VERSION_MATRIX.csv:2:commit=TBD
- tables/MODEL_VERSION_MATRIX.csv:4:replaces_model=TBD
- tables/MODEL_VERSION_MATRIX.csv:4:replaced_by=TBD
- tables/MODEL_VERSION_MATRIX.csv:4:commit=TBD
- tables/DATASET_REGISTRY.csv:5:tissue=TBD
- tables/DATASET_REGISTRY.csv:5:subtissue=TBD
- tables/DATASET_REGISTRY.csv:5:checksum=TBD
- tables/RESULT_STATUS_REGISTER.csv:2:replacement_result=TBD
- tables/RESULT_STATUS_REGISTER.csv:4:replacement_result=TBD
- tables/RESULT_STATUS_REGISTER.csv:7:replacement_result=TBD
- tables/GSE_RESULT_RECONCILIATION.csv:4:pairing_status=paired in GSE package for contextual modules; Kei/Mashiro final Risk pairing is TBD

### Liu Bowen
- tables/RESULT_STATUS_REGISTER.csv:5:replacement_result=TBD

### Mashiro Yasuda
- tables/MODEL_VERSION_MATRIX.csv:6:version=TBD
- tables/MODEL_VERSION_MATRIX.csv:6:replaces_model=TBD
- tables/MODEL_VERSION_MATRIX.csv:6:replaced_by=TBD
- tables/MODEL_VERSION_MATRIX.csv:6:commit=TBD
- tables/MODEL_VERSION_MATRIX.csv:6:notes=No numeric result CSV or script was found; pairing is TBD.
- tables/MODULE_REGISTRY.csv:11:version=TBD
- tables/MODULE_REGISTRY.csv:11:preprocessing=Reported Kei finalized settings; exact script/config TBD
- tables/MODULE_REGISTRY.csv:11:aggregation_rule=Pairing_status is TBD; current analysis described as unpaired.
- tables/MODULE_REGISTRY.csv:11:final_inference_unit=animal TBD
- tables/PARAMETER_REGISTRY.csv:78:version=TBD
- tables/PARAMETER_REGISTRY.csv:79:version=TBD
- tables/PARAMETER_REGISTRY.csv:80:version=TBD
- tables/PARAMETER_REGISTRY.csv:80:value=TBD
- tables/PARAMETER_REGISTRY.csv:80:status=TBD
- tables/PARAMETER_REGISTRY.csv:80:source_path=TBD
- tables/PARAMETER_REGISTRY.csv:80:source_type=TBD
- tables/PARAMETER_REGISTRY.csv:81:version=TBD
- tables/PARAMETER_REGISTRY.csv:81:value=TBD
- tables/PARAMETER_REGISTRY.csv:81:status=TBD
- tables/PARAMETER_REGISTRY.csv:81:source_path=TBD
- tables/PARAMETER_REGISTRY.csv:81:source_type=TBD
- tables/PARAMETER_REGISTRY.csv:82:version=TBD
- tables/PARAMETER_REGISTRY.csv:82:value=TBD
- tables/PARAMETER_REGISTRY.csv:82:status=TBD
- tables/PARAMETER_REGISTRY.csv:82:source_path=TBD
- tables/PARAMETER_REGISTRY.csv:82:source_type=TBD
- tables/PARAMETER_REGISTRY.csv:83:version=TBD
- tables/PARAMETER_REGISTRY.csv:83:value=TBD
- tables/PARAMETER_REGISTRY.csv:83:status=TBD
- tables/PARAMETER_REGISTRY.csv:83:source_path=TBD
- tables/PARAMETER_REGISTRY.csv:83:source_type=TBD
- tables/RESULT_STATUS_REGISTER.csv:10:model_version=TBD
- tables/RESULT_STATUS_REGISTER.csv:10:primary_statistical_unit=animal TBD; current analysis described as unpaired
- tables/RESULT_STATUS_REGISTER.csv:10:replacement_result=TBD
- tables/RESULT_STATUS_REGISTER.csv:10:notes=Exact CSV, script, threshold, and pairing are TBD.
- tables/GSE_RESULT_RECONCILIATION.csv:5:known_animal_scope=TBD
- tables/GSE_RESULT_RECONCILIATION.csv:5:pairing_status=TBD; current analysis described as unpaired
- tables/GSE_RESULT_RECONCILIATION.csv:5:cell_summary=exact cell-level values TBD; gene coverage reported >95% for all axes
- tables/GSE_RESULT_RECONCILIATION.csv:5:threshold_definition=TBD: per-axis vs joint; age-specific vs pooled; donor-specific vs group-level
- tables/GSE_RESULT_RECONCILIATION.csv:5:scoring_version=Kei reported Scanpy score_genes settings; script/config TBD

### TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:dataset_id=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:source_path=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:owner=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:species=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:tissue=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:subtissue=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:cell_type=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:assay=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:n_cells=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:n_genes=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:independent_donor_count=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:donors_per_age_group=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:sex_balance=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:age_sex_confounding=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:dataset_age_confounding=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:assay_age_confounding=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:tissue_compatibility=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:msc_annotation_compatibility=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:raw_counts_available=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:donor_resolved_metadata=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:paired_sample_metadata=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:genome_build=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:gene_identifier_namespace=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:transcript_complexity_by_age=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:kaile_feature_coverage=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:zihan_signature_coverage=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:kei_risk_coverage=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:eligibility=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:exclusion_reason=TBD
- tables/DATASET_CANDIDATE_INTAKE_TEMPLATE.csv:2:notes=TBD

### Zihan Zhou
- tables/MODEL_VERSION_MATRIX.csv:3:replaces_model=TBD
- tables/MODEL_VERSION_MATRIX.csv:3:replaced_by=TBD
- tables/MODEL_VERSION_MATRIX.csv:3:commit=TBD
- tables/PARAMETER_REGISTRY.csv:33:value=TBD
- tables/PARAMETER_REGISTRY.csv:33:status=TBD
- tables/PARAMETER_REGISTRY.csv:33:source_path=TBD
- tables/PARAMETER_REGISTRY.csv:33:source_type=TBD
- tables/RESULT_STATUS_REGISTER.csv:3:replacement_result=TBD

### Zihan Zhou and Kaile Zhu
- tables/RESULT_STATUS_REGISTER.csv:13:replacement_result=TBD

## Current Active Result Chain

- kaile_v1_1_internal_limb_facs_cv: active_frozen | Gene-signature donor ROC-AUC 0.9792, Brier 0.0547, age Spearman -0.8767 across 14 donors.
- zihan_v2_1_nested_lomo: active_frozen | Primary nested LOMO AUC 0.979, all-age rho -0.633, old-only rho 0.655, library-size rho -0.521.
- kaile_v1_1_limb_droplet_sensitivity: active_descriptive_sensitivity | Donor age Spearman -0.838 in Kaile packaged 16-donor TMS Droplet trajectory sensitivity analysis.
- kaile_v1_1_vs_zihan_v2_1_cleaned_sensitivity: active_descriptive_sensitivity | Both models preserve Young-Old direction in cleaned FACS and the 12-donor local Droplet subset; score units and estimands differ. Bootstrap intervals do not establish an overall winner.
- kaile_gse176206_msc_external_youth_response: active_descriptive | Young control minus aged control +0.0359; aged SOKM minus aged control -0.1486; minimum gene overlap 0.942.
- gse176206_state_remodeling_dominant: active_descriptive | Aged SOKM reduced Youth Score, MSC identity, and ECM/collagen; baseline/Col11a1 states decreased while reprogramming states increased; global cell-cycle activation is method-sensitive.
- gse176206_contextual_risk_axes: active_descriptive | Contextual identity-loss risk increased, inflammation/SASP risk decreased, and abnormal-proliferation/DDR was method-sensitive in the aged integrated assessment.
- kei_risk_score_method_freeze: frozen_method_pending_script_verification | Kei parameters are registered from team instruction and master gene list was audited; scoring script/config and threshold scope are TBD.
- mashiro_preliminary_gse176206_risk_rerun: active_preliminary_descriptive | Preliminary visual: gene coverage reported >95%; each point is one animal; SOKM lower Dedifferentiation and Inflammation; Genomic Stress slight or mixed.
- geneformer_axis_validation_diagnostic: diagnosed_negative | Raw donor LOO 11/14 p=0.050; after detected-gene regression 8/14 p=0.531; matched-depth band 10/14 p=0.094.
- gse176206_two_youth_model_followup: active_descriptive | Both separately run frozen models show negative SOKM-minus-control direction in all valid young and aged animals; aggregate decrease reflects both state redistribution and selected within-state remodeling.
- geneformer_fixed_length_sweep: active_descriptive | 256-token same-cell experiment retained all 815 cells and dropped 0 old cells: fixed LOO 13/14, p=0.001, confound r=0.619; variable-length control 11/14, p=0.050, r=0.691.
- geneformer_donor_pseudobulk_axis_test: diagnosed_negative | Donor pseudobulk axis separation = 7/14.

## Blocked Dependencies

- Repository HEAD is unborn, so source commit fields remain TBD.
- Kei final Risk scoring script/config is missing locally.
- Mashiro preliminary Risk result CSV, script/config, threshold scope, and pairing proof are missing locally.
- Zihan v2.1 explicit coverage acceptance threshold is TBD.
- Geneformer raw artifacts/scripts are not present outside the PDF evidence.

## Paper Consistency Warnings

- Summary_draft.pdf contains Tumor/tumor axis language requiring final genomic-stress wording review.
- Summary_draft.pdf mentions Geneformer factor ranking; ensure diagnosed-negative status supersedes ranking.

## Audit Outcome

PASS with warnings.

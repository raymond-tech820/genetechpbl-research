# Owner Action Items

These items collect all material `TBD` fields and owner confirmations required before Phase 2 or manuscript integration.

## Kaile Zhu

- Confirm the source commit or immutable release identifier for `facs_msc_youth_score_v1_1_cleaned_limb`.
- Confirm whether the GSE176206 condition-level validation output should be treated as superseded by the revised known-animal paired package, or retained as a separate descriptive validation layer.
- Confirm the formal comparison artifact's author/release identifier if one exists; the local `model_comparison_v1_1_v2_1/` package is now registered as the primary comparison evidence.
- Confirm whether any additional Kaile contextual modules in the GSE package should be registered in Phase 2 interface compatibility.
- Confirm that contextual `risk_identity_loss`, `risk_inflammation_sasp`, and `risk_abnormal_proliferation_ddr` remain separate from Kei final Risk Score.
- Confirm the packaged 16-donor Droplet trajectory gene count and exact subtissue label, and confirm it remains distinct from the local 12-donor shared comparison subset.
- Confirm the GSE176206 MSC dataset tissue/subtissue labels to use in dataset registry.
- Confirm whether the source h5ad checksum recorded inside `scoring_manifest.json` should be copied into `DATASET_REGISTRY.csv` or left as a non-local source checksum note.

## Zihan Zhou

- Confirm whether Zihan v3 will be frozen and added later. No independent v3 artifact was found locally.
- Provide the source commit or immutable release identifier for `facs_msc_youth_score_v2_1`.
- Confirm whether v2.1 has an explicit minimum acceptable `gene_coverage` or `weighted_coverage` threshold for deployment.
- Confirm whether v2.1 has any completed external validation artifact, or whether external use should remain sensitivity/transportability only.
- Confirm how future v3, if any, should replace or coexist with v2.1 in the model matrix.

## Hasegawa Kei

- Provide the exact final scoring script or config for `kei_risk_score_v1`.
- Confirm that the final axis labels are identity-loss-associated, SASP-associated, and genomic-stress-associated.
- Confirm whether the source CSV column `Tumor` should be renamed or mapped to genomic-stress-associated in public-facing materials.
- Confirm whether the 95th percentile cutoff is computed separately per axis or jointly.
- Confirm whether the cutoff is computed separately by age group or across all controls.
- Confirm whether matched Control means donor-specific, group-level, age-matched, condition-matched, or another rule.
- Confirm whether all controls are pooled.
- Confirm whether the primary animal-level summaries are median, p95, high-risk fraction, or all three.
- Provide the Scanpy environment proof for version 1.12.2 and the final random seed handling.

## Mashiro Yasuda

- Can Control and SOKM be paired by the same animal ID?
- Where is the exact result CSV behind `animal_risk_summaries mashiro.png`?
- Where is the scoring script/config used for the preliminary rerun?
- Was the 95th percentile cutoff calculated per axis?
- Was the cutoff calculated by age group or across all controls?
- What animal-level summaries are primary: median, p95, high-risk fraction, or all three?
- Provide exact gene-coverage values for all three axes rather than the current reported `>95%`.

## Jia Qi Choy

- Provide the Geneformer scripts, ortholog mapping table, token dictionary/version manifest, and raw output tables if they exist locally outside this repository.
- Confirm that the old seven-factor ranking is superseded and should not be used as a formal result.
- Provide the raw fixed-length sweep CSV and fixed-length script/config for the reported 256-token result.
- Provide donor-pseudobulk diagnostic script/config, raw output, repository commit/checksum, and the missing `axis_length_report.pdf`.
- Obtain advisor confirmation for the proposed criterion: any future depth-robust model must demonstrate donor-level separation that survives explicit evaluation of detected-gene count (`proposed_pending_advisor_confirmation`).
- Confirm that Geneformer output should remain incompatible with Youth/Risk parser input unless a future predicted-expression interface is created.

## Liu Bowen / Integration Owner

- Keep Phase 2 blocked until Phase 1 is reviewed and approved.
- Do not create a Youth model ensemble or merged score.
- Keep `TBD` visible in all registries until the responsible owner resolves it.
- Populate the `DATASET_CANDIDATE_INTAKE_TEMPLATE.csv` placeholders only when a candidate dataset is proposed; its `TBD` cells are template fields, not registered data.
- When Phase 2 starts, produce interface compatibility and schema files without changing source model artifacts.
- Review `Summary_draft.pdf` for provisional `Safe Zone`, `Tumor` axis, and Geneformer ranking language before manuscript integration.

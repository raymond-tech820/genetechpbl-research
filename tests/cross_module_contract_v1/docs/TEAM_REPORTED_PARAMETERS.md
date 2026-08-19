# Team-Reported Parameters

This file replaces non-reproducible Codex attachment paths. These records are reported settings or governing decisions that still require local code/config verification; none are upgraded to `confirmed_from_code`.

| owner | module | parameter | reported value | current verification status | unresolved question |
|---|---|---|---|---|---|
| Integration owner | cross_module_youth_contract_v1 | governing principle | Align module outputs without assuming mathematical equivalence | reported_frozen_pending_code_verification | Confirm in the next reviewed contract release. |
| Integration owner | cross_module_youth_contract_v1 | cross-model operations | Compare direction, donor ranking, age association, SOKM response direction, standardized within-model effects, technical association, and pathway concordance; do not merge or rescale raw scores | reported_frozen_pending_code_verification | Owner approval pending. |
| Kei | kei_risk_score_v1 | normalization | normalize_total target_sum=10000, then log1p | reported_frozen_pending_code_verification | Provide final scoring script/config. |
| Kei | kei_risk_score_v1 | score_genes settings | ctrl_size=50, n_bins=25, use_raw=False, seed=42, Scanpy 1.12.2 | reported_frozen_pending_code_verification | Provide environment proof and exact script. |
| Kei | kei_risk_score_v1 | threshold | 95th percentile of matched Control cells | reported_frozen_pending_code_verification | Confirm per-axis vs joint, age scope, pooling, and paired-control rule. |
| Mashiro | mashiro_gse176206_kei_risk_rerun | preliminary direction | SOKM lower Dedifferentiation and Inflammation; Genomic Stress slight or mixed; gene coverage reported >95% | reported_pending_source_artifact | Provide numeric CSV, script/config, pairing proof, and exact coverage. |
| Jia Qi | geneformer_fixed_length_sweep | 256-token follow-up | 815 cells retained, 0 old cells dropped, LOO 13/14, p=0.001, confound r=0.619; variable-length control 11/14, p=0.050, r=0.691 | reported_pending_source_artifact | Provide raw sweep CSV, script/config, and repository checksum. |
| Jia Qi | geneformer_donor_pseudobulk_axis_test | donor validation | donor pseudobulk axis separation 7/14 | reported_pending_source_artifact | Provide raw diagnostic output, script/config, and advisor acceptance decision. |
| Jia Qi | future depth-robust model | proposed acceptance criterion | donor-level separation must survive explicit evaluation of detected-gene count | proposed_pending_advisor_confirmation | Advisor confirmation is required before future model selection. |


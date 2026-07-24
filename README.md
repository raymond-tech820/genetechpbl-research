# genetechpbl-research

# Genentech PBL Youth Score Models

This repository temporarily contains MSC Youth Score model deliverables developed for the Genentech AI for Regenerative Biology PBL project.

| Model | Contributor | Description |
|---|---|---|
| [`droplet_msc_youth_score_v1`](./droplet_msc_youth_score_v1/) | Zihan Zhou | Droplet-trained model retained as a sensitivity/reference model rather than the current primary FACS model. |
| [`facs_msc_youth_score_v1`](./facs_msc_youth_score_v1/) | Kaile Zhu | Deprecated FACS v1 model, retained for provenance only. Use `facs_msc_youth_score_v1_1_cleaned_limb` for the current Kaile model. |
| [`facs_msc_youth_score_v1_1_cleaned_limb`](./facs_msc_youth_score_v1_1_cleaned_limb/) | Kaile Zhu | Current Kaile FACS model for comparison and further refinement, retrained after explicit exclusion of age-confounded diaphragm cells. |
| [`facs_msc_youth_score_v2`](./facs_msc_youth_score_v2/) | Zihan Zhou | Deprecated FACS v2 model, retained for provenance only. Use `facs_msc_youth_score_v2_1` for the current Zihan model. |
| [`facs_msc_youth_score_v2_1`](./facs_msc_youth_score_v2_1/) | Zihan Zhou | Current Zihan FACS model for comparison and further refinement, retrained on diaphragm-excluded TMS FACS limb-muscle MSC pseudobulk data. |

The current active comparison should focus on Kaile Zhu's `facs_msc_youth_score_v1_1_cleaned_limb` and Zihan Zhou's `facs_msc_youth_score_v2_1`. The older FACS v1/v2 folders are preserved only to keep prior analyses reproducible.

See each model directory for its detailed documentation, scoring files, validation results, and limitations.

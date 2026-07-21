# genetechpbl-research

# Genentech PBL Youth Score Models

This repository temporarily contains three MSC Youth Score models developed for the Genentech AI for Regenerative Biology PBL project.

| Model | Contributor | Description |
|---|---|---|
| [`droplet_msc_youth_score_v1`](./droplet_msc_youth_score_v1/) | Zihan Zhou | An early model trained on TMS Droplet limb-muscle MSC data. Because the training samples are substantially unbalanced, this model is retained mainly as a reference. |
| [`facs_msc_youth_score_v1`](./facs_msc_youth_score_v1/) | Kaile Zhu | A finalized Youth Score model developed using FACS MSC data. |
| [`facs_msc_youth_score_v2`](./facs_msc_youth_score_v2/) | Zihan Zhou | A finalized Youth Score model trained on TMS FACS limb-muscle MSC pseudobulk data, with nested donor-level validation and cross-assay sensitivity analysis. |

The two finalized models—Kaile Zhu's FACS MSC Youth Score v1 and Zihan Zhou's FACS Limb MSC Youth Score v2—were developed independently. A direct comparison between them has not yet been completed.

See each model directory for its detailed documentation, scoring files, validation results, and limitations.

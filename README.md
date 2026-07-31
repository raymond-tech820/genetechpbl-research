# genetechpbl-research

# Genentech PBL Youth Score Models

This repository temporarily contains MSC Youth Score model deliverables and model-comparison outputs developed for the Genentech AI for Regenerative Biology PBL project.

| Model | Contributor | Description |
|---|---|---|
| [`droplet_msc_youth_score_v1`](./droplet_msc_youth_score_v1/) | Zihan Zhou | Droplet-trained model retained as a sensitivity/reference model rather than the current primary FACS model. |
| [`facs_msc_youth_score_v1`](./facs_msc_youth_score_v1/) | Kaile Zhu | Deprecated FACS v1 model, retained for provenance only. Use `facs_msc_youth_score_v1_1_cleaned_limb` for the current Kaile model. |
| [`facs_msc_youth_score_v1_1_cleaned_limb`](./facs_msc_youth_score_v1_1_cleaned_limb/) | Kaile Zhu | Current Kaile FACS model for comparison and further refinement, retrained after explicit exclusion of age-confounded diaphragm cells. |
| [`facs_msc_youth_score_v2`](./facs_msc_youth_score_v2/) | Zihan Zhou | Deprecated FACS v2 model, retained for provenance only. Use `facs_msc_youth_score_v2_1` for the current Zihan model. |
| [`facs_msc_youth_score_v2_1`](./facs_msc_youth_score_v2_1/) | Zihan Zhou | Current Zihan FACS model for comparison and further refinement, retrained on diaphragm-excluded TMS FACS limb-muscle MSC pseudobulk data. |

## Active Model Comparison

The current active comparison is between Kaile Zhu's `facs_msc_youth_score_v1_1_cleaned_limb` and Zihan Zhou's `facs_msc_youth_score_v2_1`:

```text
model_comparison_v1_1_v2_1/
```

Start with:

```text
model_comparison_v1_1_v2_1/FACS_and_Droplet_v1_1_vs_v2_1_combined_comparison_report.md
```

The comparison folder includes:

- FACS same-input frozen scorer application;
- Droplet cross-assay sensitivity;
- internal validation summary separating Kaile packaged donor OOF CV from Zihan nested LOMO;
- paired donor bootstrap uncertainty analysis;
- model concordance decomposition within and after age adjustment;
- technical-variable audits;
- FACS gene-signature overlap analysis.

Important interpretation boundary: the comparison supports directionally consistent young-old state separation for both active models, but it does not prove that either model is technically independent or externally validated. Paired donor bootstrap does not establish statistically supported overall superiority of Kaile v1.1 over Zihan v2.1, although Kaile v1.1 shows more consistent descriptive age ordering.

The older FACS v1/v2 folders are preserved only to keep prior analyses reproducible. See each model directory and the comparison directory for detailed documentation, scoring files, validation results, and limitations.

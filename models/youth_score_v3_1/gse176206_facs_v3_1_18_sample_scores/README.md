# GSE176206 FACS Youth Score v3.1: 18-Sample M1/M4 Scores

## Purpose

This folder supplies the machine-readable donor-level output requested for the frozen GSE176206 application. The table has one row per:

`age group x exact treatment arm x animal label`

The 18 rows comprise two age groups, three exact treatment arms, and three nested animal labels per age-treatment arm.

## Files

- `gse176206_m1_m4_18_sample_scores.csv`: wide-format M1/M4 score and QC table.
- `generation_audit.csv`: row-count, uniqueness, score-equivalence, and coverage checks.

## Biological Unit

`biological_unit_id` is the authoritative key:

`age_group|exact_treatment_arm|animal_label`

Numeric animal labels are nested within age and exact treatment arm. For example, animal `1` under `Aged|Tg+/Dox+` and animal `1` under `Aged|Tg+/Dox-` are distinct biological units. The treatment arms are therefore unpaired and must not be joined by numeric animal label alone.

## Score Columns

- `m1_calibrated_score`: frozen M1 primary operational score.
- `m4_calibrated_score`: frozen M4 high-stringency sensitivity score.
- `m1_*` and `m4_*`: model identifiers, versions, roles, module scores, raw scores, coverage, and coverage-pass fields.

Higher calibrated values indicate a more young-like state relative to the FACS v3.1 training calibration. Scores were generated with frozen genes, weights, training means and standard deviations, calibration centers, orientation, and coverage thresholds. GSE176206 was not used for retraining or recalibration.

## QC Fields

- `n_cells` and `all_gene_library_size` describe each arm-level pseudobulk unit.
- `gene_coverage` is the fraction of frozen signature genes available.
- `weighted_coverage` is the fraction of total frozen signature weight available.
- young-high and old-high weighted coverage are reported separately.
- `sample_qc_pass` requires positive cell/library counts, finite M1/M4 scores, and both frozen coverage gates to pass.

All 18 rows pass the frozen M1 and M4 coverage gates and numeric sample QC. `generation_audit.csv` confirms exact calibrated-score equivalence to the original frozen scoring output.

## Provenance

Source expression object:

`data_geo/processed/GSE176206_msc_sokm.h5ad`, `layers['counts']`

Authoritative frozen score source:

`outputs/facs_v3_1/external_validation/gse176206/gse176206_arm_scores_m1_m4.csv`

This deliverable reshapes the existing frozen output from 36 model-specific rows into 18 biological-unit rows. It does not recompute or alter any score.

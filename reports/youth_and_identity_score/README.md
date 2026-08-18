# My Youth and Identity LaTeX Manuscript Package

This folder is an English-only manuscript contribution for My Youth Score and MSC Identity Score. The manuscript author is Kaile Zhu.

## Contents

- `main.tex`: two-column manuscript source.
- `references.bib`: bibliography used by the manuscript.
- `main.pdf`: compiled review PDF.
- `figures/`: all figures referenced by `main.tex`.
- `source_data/`: frozen result tables used to generate the new validation and Identity figures.
- `make_figures.py`: reproducible figure-generation script.

## Layout contract

The body remains in two-column format. Figures and tables span the full text width through controlled `strip` blocks rather than unconstrained floating environments. This keeps each visual next to its substantive discussion and prevents figures from accumulating on float-only pages. Figure 1 uses a vertically stacked A/B arrangement to improve legibility and balance the validation page; the remaining figures retain side-by-side panels where that arrangement best supports direct comparison.

The compiled review PDF is five pages with the following narrative flow:

1. Youth Score methods;
2. Youth validation, Table 1, and the exact-arm GSE176206 result setup;
3. GSE176206 figures and Table 2, Youth interpretation, and the start of Identity methods;
4. Identity methods and results, Table 3, Identity interpretation, and Figure 3;
5. the integrated framework, limitations, and references.

This layout requires the LaTeX `cuted` package, which is included in standard TeX Live installations.

## Scientific scope

This contribution covers only:

1. My Youth Score;
2. the frozen MSC Identity Score; and
3. their integrated interpretation and limitations.

The independently developed Youth model in the shared repository should be described by its authors. Risk scoring is outside this package.

## Compile

Run from this folder:

```powershell
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

To regenerate the new figures before compiling:

```powershell
python make_figures.py
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The two original GSE176206 Youth figures are copied from the finalized 18-unit exact-arm package. The other two data figures are regenerated from the frozen CSV files included under `source_data/`.

## Interpretation boundary

A lower Youth Score means reduced similarity to the young TMS FACS Limb_Muscle MSC reference. It is not evidence of accelerated aging. Identity is an independent explanatory axis and is not a mathematical correction to Youth. The frozen Identity panel is a composite of murine skeletal/stromal-lineage-associated genes whose expression is subset- and state-dependent; it is not a set of universally expressed or MSC-exclusive markers.

## Cross-module handoff

Cross-module integration should occur only after module-native animal-level aggregation. Join the frozen Identity output on `biological_unit_id`. If another module does not provide that field, construct and validate the equivalent composite key (`dataset_id`, `age_group`, `exact_treatment_arm`, `animal_label`); never join on `animal_label` alone. Preserve Youth Score, primary Identity Score, rank-based Identity Score, and all module-specific quality-control fields as separate named outputs on their native scales. Do not rescale one module with another, impute a missing module output from a different score, or collapse the axes into a scalar composite without a separately validated integration model.

For export to `cross_module_contract_v1`, use the following field mapping:

- `dataset_id` -> `dataset_id`;
- `animal_label` -> `animal_id`;
- `age_group` -> `age_group`;
- `exact_treatment_arm` -> `condition` (and, if required by the integrator, `treatment`);
- each Identity metric column -> one long-form `score_name` / `score_value` row;
- `n_cells`, `identity_gene_coverage`, and `qc_status` -> the corresponding contract fields;
- `cross_arm_pairing_status` -> `pairing_status`.

Do not derive `paired_group_id` from repeated numeric animal labels across arms. Populate `analysis_id`, `module_id`, `module_version`, `native_scoring_path`, `evidence_level`, and `source_file` from the shared frozen-module registry at integration time.

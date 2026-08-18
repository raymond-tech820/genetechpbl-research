# Youth Score Manuscript TeX Report

This folder is a self-contained LaTeX delivery of the Youth Score manuscript module.

## Contents

- `Youth_Score_Manuscript_Module.tex`: portable two-column LaTeX manuscript.
- `Youth_Score_Manuscript_Module.pdf`: verified compiled manuscript.
- `Youth_Score_Manuscript_Module.md`: canonical Markdown content used for the consistency audit.
- `figures/`: the six raster figures referenced by the TeX source.
- `CONTENT_AUDIT.md`: text, table, formula, figure, compilation, and visual-QA record.
- `MANIFEST.sha256`: SHA-256 checksums for every delivered file.

## Compile

Run from this folder:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error Youth_Score_Manuscript_Module.tex
```

The document uses standard TeX Live packages: `geometry`, `newtxtext`, `newtxmath`, `microtype`, `graphicx`, `booktabs`, `array`, `siunitx`, `caption`, `subcaption`, `cuted`, `tabularx`, `placeins`, `natbib`, and `hyperref`.

## Portable Table Contract

All 12 tables are non-floating single-column blocks with fixed `\columnwidth`, explicit font sizes, and `tabularx` or `tabular*` column allocation. No table uses `table`, `table*`, or `resizebox`. Wide source tables are divided into labeled panels under one table number, preserving every source field and value. This keeps each table independent when copied into a larger two-column manuscript.


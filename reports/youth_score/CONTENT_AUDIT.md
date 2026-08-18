# Youth Score Manuscript Delivery Audit

## Audit Scope

Canonical content source: `Youth_Score_Manuscript_Module.md`

Delivered typeset source: `Youth_Score_Manuscript_Module.tex`

The audit covered manuscript prose, section order, formulas, table fields and values, figures, citations, compilation, and rendered PDF layout.

## Content Consistency

- The abstract, Introduction, Methods, Results, Discussion, Limitations, Conclusion, Appendices, and reference content were checked against the canonical Markdown source.
- Key scientific claims and interpretive boundaries are unchanged.
- All Markdown formulas are retained in the TeX source. TeX contains additional repeated inline labels only where a wide table was divided into panels.
- All 12 Markdown tables are represented by Tables 1-12 in the TeX manuscript.
- The Markdown tables contain 209 numeric cells. Every source numeric cell is present in TeX with the same value and at least the same multiplicity. Additional repetitions are restricted to panel keys such as age, model, or stability threshold.
- Tables 2 and 3 stack their original text fields vertically for single-column readability. No source field wording or data value was removed.
- Bibliographic citations are rendered with `natbib`; the two full references are embedded in the TeX source.

## Figure Verification

All six Markdown figures are included in TeX and copied byte-for-byte into `figures/`:

| Figure file | Pixel dimensions | SHA-256 |
|---|---:|---|
| `manuscript_stability_tradeoff.png` | 2850 x 1740 | `305c5133d21ac71f8f4871b8f425d798ad077007cfc3cd589a4d4b1d7458dc9f` |
| `frozen_signature_module_composition.png` | 1000 x 800 | `97054ac8c5bd0d0ebf4257f03dff9437f1268f495901816957b1aeadc8e39091` |
| `training_nested_lomo_scores_by_age.png` | 1500 x 800 | `2d1d3a8cf573187ae79b2085f959016880b2464391442fb98bb92ee44966f16e` |
| `formal_permutation_abs_all_age_rho_null.png` | 1500 x 800 | `6c9dfd84af021f67704c1e044984242a0e1bc4366721e6162dcbb8245dbdef55` |
| `gse176206_scores_by_exact_treatment.png` | 1440 x 800 | `47e09eb078f0b59327e2baec2b37a9bc12604e6edfaeebc43c3ecce80ca7c8eb` |
| `droplet_scores_by_age.png` | 1200 x 880 | `8ce321c9a3ed65e96d7c6ff483951cbaec2b2cb15992aceae89205c9fea93ff4` |

## Table Stability

- 12 of 12 tables use non-floating `minipage` containers.
- 0 tables use `table` or `table*` floats.
- 0 tables use `resizebox` or implicit font scaling.
- Every table has an explicit font size and a fixed `\columnwidth` layout.
- Wide tables retain one caption and one table number while using multiple panels.

## Compilation And Visual QA

- Compiled with `latexmk -pdf -interaction=nonstopmode -halt-on-error`.
- No LaTeX errors, undefined references, undefined citations, or overfull boxes were present in the final log.
- All eight rendered pages were inspected at 120 dpi.
- No clipped text, overlapping elements, broken tables, missing figures, abnormal blank pages, or unreadable glyphs were observed.
- Page numbering, section order, table numbering, figure numbering, and references were intact.


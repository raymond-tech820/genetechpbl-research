# Tests and Comparisons

This directory contains sensitivity analyses, historical model comparisons,
and cross-module compatibility audits. These packages are evaluation resources,
not current deployable models.

| Package | Purpose |
|---|---|
| [`youth_score_v1_sensitivity`](./youth_score_v1_sensitivity/) | Historical Youth Score v1 sensitivity and reference analysis |
| [`youth_score_v1_1_vs_v2_1_comparison`](./youth_score_v1_1_vs_v2_1_comparison/) | Historical comparison of the v1.1 and v2.1 Youth Score implementations |
| [`cross_module_contract_v1`](./cross_module_contract_v1/) | Reproducibility, interface-compatibility, provenance, and interpretation-boundary audit |

The historical model comparison does not establish universal superiority,
technical independence, or external validation of either model. Consult its
reports for the applicable evidence boundaries.

Run the cross-module contract audit from the repository root:

```bash
python tests/cross_module_contract_v1/scripts/audit_contracts.py

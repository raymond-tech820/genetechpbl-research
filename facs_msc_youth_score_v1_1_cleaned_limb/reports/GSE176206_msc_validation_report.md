# GSE176206 External Validation Report

Generated: 2026-07-24T01:09:12.580491+00:00

- Source: `data\external\GSE176206\GSE176206_msc_sokm.h5ad`
- Model: `gene_signature` from `outputs\youth_score\limb_facs`
- Minimum model-gene overlap: `0.942`
- TMS/external reference-expression Spearman: `0.520`
- Cell-identity evidence: `source_file_declared`
- External replication status: `replicated`
- Young-control sanity direction: `True`
- Aged-SOKM rejuvenation direction: `False`

## Donor/replicate-aggregated conditions

| external_age_group   | external_treatment   |   youth_score |   replicates |   cells |   ci_low |   ci_high |
|:---------------------|:---------------------|--------------:|-------------:|--------:|---------:|----------:|
| aged                 | SOKM                 |        0.344  |            3 |    8181 |   0.2958 |    0.4171 |
| aged                 | control              |        0.4925 |            3 |    2012 |   0.484  |    0.5045 |
| young                | SOKM                 |        0.3658 |            3 |    6811 |   0.3273 |    0.4013 |
| young                | control              |        0.5284 |            3 |    1682 |   0.5116 |    0.5393 |

The external files were not used for feature selection, fitting, early stopping, or calibration. Confidence intervals use 2,000 biological-replicate bootstrap draws only when at least two identifiable replicates are available in each contrasted condition. A `replication_limited` result is descriptive because cell-level variation is not a substitute for biological replication. These results measure expression similarity to the TMS young reference and do not establish safety, preserved cell identity, causality, or treatment efficacy.

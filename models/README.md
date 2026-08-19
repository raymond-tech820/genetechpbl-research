# Models

This directory contains current or maintained model packages. Each package
retains its own scoring contract, model card, version provenance, validation
evidence, and interpretation boundaries.

| Package | Contributor | Role |
|---|---|---|
| [`youth_score_v1_1`](./youth_score_v1_1/) | Kaile Zhu | Current cleaned-limb Youth Score package; internally documented as the v1.1 release |
| [`youth_score_v2_1`](./youth_score_v2_1/) | Zihan Zhou | Maintained historical baseline for reproducibility and model comparison |
| [`youth_score_v3_1`](./youth_score_v3_1/) | Zihan Zhou | Current frozen Youth Score framework; M1 is primary and M4 is a high-stringency sensitivity model |
| [`geneformer`](./geneformer/) | Jia Qi Choy | Zero-shot Geneformer perturbation and embedding analysis |
| [`identity_score`](./identity_score/) | Kaile Zhu | Frozen knowledge-driven MSC Identity Score |
| [`risk_score`](./risk_score/) | Kei Hasegawa | Risk Score module |


Model outputs are not automatically interchangeable. Youth, Identity, Risk,
and Geneformer outputs have different physical meanings, normalization
contracts, and inference units. Read the package README or model card before
application or comparison.

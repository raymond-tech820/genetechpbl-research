# FACS Youth Score v3.1 M1-M4 External Validation Synthesis

## Final role assignment

The evidence most closely matches frozen Scenario A. M1 remains the primary deployable reference. M4 remains a post-ablation sensitivity comparator and is not upgraded.

This is a role assignment, not a claim that M1 and M4 are mathematically equivalent or that M1 is universally superior.

## Evidence hierarchy

### FACS statistical evidence

The fully nested sex-stratified complete-age permutation test was significant for the predefined M1 pipeline (`p=0.023`) and for exploratory M4 (`p=0.010`). This establishes stronger age-label association than the practical randomized-label null, not technical independence or aging-clock validity.

### GSE176206 biological external application

Both models passed total and module coverage gates. In both independent control arms, baseline Young scores exceeded Aged scores. For both ages and both models, the prespecified median SOKM-minus-control contrast was negative against each control arm.

The GSE units were unpaired across treatment arms. Equal numeric animal labels were not treated as the same animal. M1's Aged SOKM-versus-Tg+/Dox- mean difference was 0.0042 while its primary median difference was -0.0148; therefore the treatment result is directional median evidence, not uniform agreement across summaries.

M1-M4 Spearman agreement was 0.930 overall and 0.907 after removing age-by-treatment centers. State-restricted contrasts were heterogeneous, so aggregate SOKM shifts cannot be described as every state moving uniformly.

Formal Youth-Identity and Youth-Risk coupling remained not assessable because no standalone frozen scorers existed before outcome analysis.

### TMS Droplet cross-assay transportability

Both models had complete signature coverage, AUC 1.0, negative all-age and old-only correlations, and Young medians above Old medians. M1 and M4 overall donor Spearman agreement was 0.965; exact-age-residualized agreement was 0.846.

M4 had slightly more negative age-correlation point estimates, but every prespecified paired-bootstrap interval for M4-minus-M1 included zero. This does not support upgrading M4 from its exploratory role.

Droplet score correlations with all-gene library size, filtered-DGE library size, TMM effective library size, and cell count were small. Detected-gene correlations were moderate (M1 0.336; M4 0.308). FACS nested scores remained substantially associated with library size and cell count, so technical independence has not been demonstrated.

## Combined evidence table

```
 model                         frozen_role signature_genes facs_nested_auc
    M1        primary_deployable_reference              58               1
    M4 post_ablation_exploratory_secondary              29               1
 facs_all_age_rho facs_old_only_rho facs_raw_library_rho facs_cell_count_rho
       -0.8063484        -0.1091089           -0.6923077          -0.5566560
       -0.8438530        -0.3273268           -0.7274725          -0.5588562
 permutation_abs_rho permutation_z_null permutation_empirical_p
           0.8063484           1.952626                   0.023
           0.8438530           2.176625                   0.010
 gse_min_module_weighted_coverage gse_young_higher_both_controls
                        0.8702639                           TRUE
                        0.9224784                           TRUE
 gse_sokm_lower_both_controls_young gse_sokm_lower_both_controls_aged
                               TRUE                              TRUE
                               TRUE                              TRUE
 droplet_min_module_weighted_coverage droplet_auc droplet_all_age_rho
                                    1           1          -0.8079193
                                    1           1          -0.8297550
 droplet_old_only_rho droplet_all_gene_library_rho
           -0.6617241                  -0.02797203
           -0.7006490                  -0.03496503
 droplet_effective_library_rho droplet_cell_count_rho
                   -0.04895105             -0.0979021
                   -0.05594406             -0.1048951
 droplet_detected_genes_rho
                  0.3356643
                  0.3076923
```

## M1-M4 agreement table

```
   dataset                      scope  n  spearman
 GSE176206                    overall 18 0.9298246
 GSE176206 age_treatment_residualized 18 0.9071207
   Droplet                    overall 12 0.9650350
   Droplet     age_group_residualized 12 0.9300699
   Droplet     exact_age_residualized 12 0.8461538
```

## Scientific conclusion

M1 and M4 detect a largely shared MSC youth-associated transcriptional axis that survives the GSE biological context and the TMS assay shift. The stricter M4 threshold changes effect magnitudes modestly but does not add a clearly separable external signal under the frozen comparisons. M1 therefore remains the operational primary model, while M4 is retained as a high-stringency sensitivity analysis.

The completed evidence supports non-random age-label association, external directional consistency, perturbation sensitivity, and within-TMS cross-assay transportability. It does not establish a universal MSC aging clock, continuous trajectory validity, technical independence, rejuvenation, treatment safety, or formal Identity/Risk coupling.

## Provenance

GSE and Droplet analyses used only processed inputs frozen in the protocol. No external result changed model genes, weights, standardization, calibration, coverage thresholds, or orientation.

## Linked Figures

The GSE176206 and Droplet figures supporting this synthesis are stored beside this report and embedded in reports 01 and 02.

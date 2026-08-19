# FACS v2 Youth Score Cross-Assay Sensitivity on TMS Droplet Limb Muscle MSC

## Scope

This is a within-TMS cross-assay sensitivity / transportability analysis, not independent external validation.

The FACS v2 frozen models were applied to Droplet mouse-level pseudobulk without changing signature genes, weights, training means/SDs, calibration centers, thresholds, or gene filtering. Droplet data were used only for scoring and evaluation.

## Droplet Input Audit

- Counts: `data_droplet/processed/tms_limb_msc_pseudobulk_counts.rds`
- Metadata: `data_droplet/processed/tms_limb_msc_pseudobulk_metadata_labeled.csv`
- Matrix dimensions: 20138 genes x 12 mice
- Mice: 12
- Age composition: 18m=4; 21m=2; 24m=4; 3m=2
- Sex x age_group mice:
       
        female male
  Old        4    6
  Young      2    0

## Signature Coverage

                                  model training_assay       module
     FACS_v2__factorial_medium_original        FACS_v2     old_high
     FACS_v2__factorial_medium_original        FACS_v2   young_high
      FACS_v2__factorial_large_original        FACS_v2     old_high
      FACS_v2__factorial_large_original        FACS_v2   young_high
  FACS_v2__factorial_stability_selected        FACS_v2     old_high
  FACS_v2__factorial_stability_selected        FACS_v2   young_high
 FACS_v2__factorial_medium_equal_weight        FACS_v2     old_high
 FACS_v2__factorial_medium_equal_weight        FACS_v2   young_high
                   FACS_v2__age_only_de        FACS_v2     old_high
                   FACS_v2__age_only_de        FACS_v2   young_high
                     Droplet_v1__Medium     Droplet_v1   old_module
                     Droplet_v1__Medium     Droplet_v1 young_module
        Droplet_v1__Medium_equal_weight     Droplet_v1   old_module
        Droplet_v1__Medium_equal_weight     Droplet_v1 young_module
 signature_genes covered_genes gene_coverage weighted_coverage
              50            50             1                 1
              50            50             1                 1
             100           100             1                 1
             100           100             1                 1
              28            28             1                 1
              29            29             1                 1
              50            50             1                 1
              50            50             1                 1
              50            50             1                 1
              50            50             1                 1
              50            50             1                 1
              50            50             1                 1
              50            50             1                 1
              50            50             1                 1

## Model Summary on Droplet

                                  model training_assay test_assay n_mice
                     Droplet_v1__Medium     Droplet_v1    Droplet     12
        Droplet_v1__Medium_equal_weight     Droplet_v1    Droplet     12
                   FACS_v2__age_only_de        FACS_v2    Droplet     12
      FACS_v2__factorial_large_original        FACS_v2    Droplet     12
 FACS_v2__factorial_medium_equal_weight        FACS_v2    Droplet     12
     FACS_v2__factorial_medium_original        FACS_v2    Droplet     12
  FACS_v2__factorial_stability_selected        FACS_v2    Droplet     12
 auc_young_vs_old all_age_rho old_only_rho young_median    old_median
                1  -0.9389333   -0.8952738    1.0000000 -2.395545e-16
                1  -0.9389333   -0.8952738    1.0008482 -1.606882e-02
                1  -0.8734263   -0.7784989    0.4752532  3.348226e-01
                1  -0.8515906   -0.7395740    0.4888861  3.612856e-01
                1  -0.8515906   -0.7395740    0.4900510  3.465016e-01
                1  -0.8734263   -0.7784989    0.4797172  3.293127e-01
                1  -0.8297550   -0.7006490    0.4820777  3.265880e-01
 young_minus_old_median raw_library_rho effective_library_rho
              1.0000000       0.1328671             0.0979021
              1.0169171       0.1328671             0.0979021
              0.1404306      -0.1608392            -0.1818182
              0.1276005      -0.2097902            -0.2377622
              0.1435494      -0.2097902            -0.2377622
              0.1504046      -0.2377622            -0.2587413
              0.1554897      -0.0979021            -0.1118881
 detected_genes_rho cell_count_rho min_gene_coverage min_weighted_coverage
          0.1818182     0.06293706                 1                     1
          0.1818182     0.06293706                 1                     1
          0.1678322    -0.23076923                 1                     1
          0.1398601    -0.27972028                 1                     1
          0.1398601    -0.27972028                 1                     1
          0.1188811    -0.30769231                 1                     1
          0.2797203    -0.16783217                 1                     1

## Primary FACS v2 vs Droplet v1

                                             comparison facs_genes
 FACS_v2_factorial_medium_original_vs_Droplet_v1_Medium        100
 droplet_v1_genes overlap_genes union_genes    jaccard
              100             7         193 0.03626943
 score_spearman_on_droplet
                 0.7762238

## Interpretation

FACS v2 primary gene coverage on Droplet was 1.000 and weighted coverage was 1.000.
FACS v2 primary on Droplet: AUC 1.000, all-age rho -0.873, old-only rho -0.778, raw-library rho -0.238, young-minus-old median 0.150.
Droplet v1 Medium on Droplet: AUC 1.000, all-age rho -0.939, old-only rho -0.895, raw-library rho 0.133, young-minus-old median 1.000.

Interpretation should separate cross-assay direction from technical independence. Directional transfer can be considered supportive only if coverage is adequate and young/old ordering is preserved; persistent library-size association remains a limitation.

## Outputs

- `cross_assay_droplet_scores_all_models.csv`
- `cross_assay_droplet_model_summary.csv`
- `cross_assay_droplet_signature_module_coverage.csv`
- `cross_assay_droplet_missing_signature_genes.csv`
- `facs_v2_vs_droplet_v1_primary_mouse_ranking.csv`
- `facs_v2_vs_droplet_v1_primary_overlap_and_score_correlation.csv`
- `figures/cross_assay_droplet_scores_by_age.png`
- `figures/cross_assay_droplet_score_vs_library_size.png`
- `figures/facs_v2_vs_droplet_v1_medium_score_correlation.png`

# FACS v2 vs Droplet v1 Pathway-Level Concordance

## Scope

This analysis asks whether the FACS v2 Medium and Droplet v1 Medium signatures, despite low gene-level overlap, enrich similar directional biological programs on a shared feature background.

No model parameters were changed. This is interpretive post-hoc annotation of frozen signatures.

## Gene Set Source

MSigDB Hallmark via `msigdbr` was attempted but Zenodo downloads repeatedly failed with partial-file errors. g:Profiler API was also attempted but did not connect reliably. To avoid fabricating pathway results, this completed run uses local Bioconductor annotation only:

       package version
 AnnotationDbi  1.74.0
  org.Mm.eg.db  3.23.0
         GO.db  3.23.1
   reactome.db  1.96.0

Collections used: GO Biological Process and Reactome.

## Unified Gene Universe

FACS filtered genes: 12859
Droplet filtered genes: 12735
Unified universe intersection: 11945

## Directional Module Sizes

                module n_genes
    FACS_v2_young_high      38
      FACS_v2_old_high      46
 Droplet_v1_young_high      31
   Droplet_v1_old_high      32

## Concordance Summary

               category n_terms
 not_shared_significant    2433

## Top Shared Significant Terms

No pathways were significant at FDR < 0.1 in both models under the local GO:BP/Reactome ORA definition.

## Top Terms By Module

              module   source       term_id
 Droplet_v1_old_high    GO:BP    GO:0006335
 Droplet_v1_old_high    GO:BP    GO:0034331
 Droplet_v1_old_high    GO:BP    GO:0010958
 Droplet_v1_old_high    GO:BP    GO:1903789
 Droplet_v1_old_high    GO:BP    GO:0045217
 Droplet_v1_old_high    GO:BP    GO:0043954
 Droplet_v1_old_high    GO:BP    GO:0060078
 Droplet_v1_old_high    GO:BP    GO:0014850
 Droplet_v1_old_high REACTOME R-MMU-8955332
 Droplet_v1_old_high    GO:BP    GO:2000463
 Droplet_v1_old_high    GO:BP    GO:0045840
 Droplet_v1_old_high    GO:BP    GO:0051955
 Droplet_v1_old_high    GO:BP    GO:0003382
 Droplet_v1_old_high    GO:BP    GO:0015804
 Droplet_v1_old_high    GO:BP    GO:0089718
 Droplet_v1_old_high    GO:BP    GO:0071825
 Droplet_v1_old_high    GO:BP    GO:0098815
 Droplet_v1_old_high    GO:BP    GO:0006334
 Droplet_v1_old_high    GO:BP    GO:0051932
 Droplet_v1_old_high    GO:BP    GO:0042391
                                                                 term_name
                              DNA replication-dependent chromatin assembly
                                                 cell junction maintenance
                    regulation of amino acid import across plasma membrane
                          regulation of amino acid transmembrane transport
                                            cell-cell junction maintenance
                                            cellular component maintenance
                             regulation of postsynaptic membrane potential
                                               response to muscle activity
 Mus musculus: Carboxyterminal post-translational modifications of tubulin
                  positive regulation of excitatory postsynaptic potential
                           positive regulation of mitotic nuclear division
                                        regulation of amino acid transport
                                             epithelial cell morphogenesis
                                              neutral amino acid transport
                                  amino acid import across plasma membrane
                                        protein-lipid complex organization
                           modulation of excitatory postsynaptic potential
                                                       nucleosome assembly
                                          synaptic transmission, GABAergic
                                          regulation of membrane potential
 overlap_n adjusted_p_value enrichment_ratio
         2        0.1082306       106.651786
         3        0.1082306        22.396875
         2        0.1082306        53.325893
         2        0.1082306        53.325893
         2        0.1082306        46.660156
         3        0.1082306        15.997768
         3        0.1082306        15.133024
         2        0.1082306        41.475694
         2        0.1164390        49.770833
         2        0.1577600        32.459239
         2        0.1814439        25.743534
         2        0.1814439        24.885417
         2        0.1814439        23.330078
         2        0.1814439        22.623106
         2        0.1814439        20.737847
         2        0.1814439        20.177365
         2        0.1814439        20.177365
         2        0.1814439        19.646382
         2        0.1814439        19.646382
         4        0.1814439         5.677281

## Interpretation Guardrails

- ORA is based on small 50-gene directional modules, so FDR-significant overlap can be sparse.
- Directional concordance is assessed separately for young-high and old-high modules; mixed-direction signatures were not pooled.
- A lack of shared significant pathways does not negate score-level transportability; it means pathway-level mechanism agreement is not strongly supported by this ORA setup.
- Hallmark/SASP/cell-cycle curated collections remain useful future additions if a stable local GMT source is provided.

## Outputs

- `pathway_ora_results_all_modules.csv`
- `pathway_direction_concordance_matrix.csv`
- `pathway_shared_significant_concordance.csv`
- `pathway_concordance_summary.csv`
- `pathway_top20_by_module.csv`
- `pathway_term_catalog_used.csv`
- `pathway_signed_ora_concordance.png`

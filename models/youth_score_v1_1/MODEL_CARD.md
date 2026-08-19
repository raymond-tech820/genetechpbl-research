# Model Card: corrected FACS Limb Muscle MSC Youth Score

## Intended use

The model scores mouse Limb Muscle MSC expression profiles by similarity to the young (3-month) TMS FACS reference. Higher scores indicate more young-like expression within this reference system.

## Training and validation design

- Input: raw cell-by-gene counts.
- Reference cohort: 815 TMS FACS Limb Muscle MSCs from 14 mice after explicit diaphragm exclusion.
- Internal validation: five donor-grouped, age-stratified outer folds. No mouse contributes cells to more than one of train, validation, or test within a fold.
- Official scorer: fold-ensemble gene-signature model, selected by donor-level out-of-fold ROC-AUC and Brier score.
- Technical-only model: diagnostic control using library size, detected genes, and sex; not a formal candidate scorer.

## Preprocessing and output

For direct cell scoring, the code performs library-size normalization to 10,000 counts followed by `log1p`; do not supply pre-normalized or log-transformed inputs. The output `youth_score` is the mean calibrated `P(young)` over fold-specific models. For biological inference, aggregate scores within true mouse/donor or use the prescribed pseudobulk workflow; cells are not independent biological replicates.

## Limitations

- The cohort has only 14 biological donors, and age, sex, assay, and technical covariates remain partly confounded.
- Cross-platform scores should be interpreted primarily as within-dataset contrasts, not as directly comparable absolute probabilities.
- The model does not establish cellular identity preservation, causal rejuvenation, clinical age, treatment efficacy, or safety.
- The external SOKM result is specific to this transcriptional score and should not be interpreted as evidence of broad treatment harm.

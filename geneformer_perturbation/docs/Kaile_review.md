# Kaile: three-model Youth Score on TMS FACS SCAT + Limb Muscle MSC — Review

## 1. What was built (from the three markdown files in the zip)

- **Data.** Two independent FACS MSC populations: SCAT adipose MSC (1,570 cells, 7 young + 8 old mice) and Limb Muscle MSC (935 cells, 6 young + 8 old). Young = 3 mo; Old = 18/24 mo.
- **Models compared.** (i) Gene signature (bioinformatic score → logistic), (ii) Elastic Net on gene expression, (iii) Gene-token Transformer with expression-bin and rank tokens + auxiliary age head, (iv) Technical-only logistic (total counts, detected genes, sex) as a diagnostic control — not eligible to be the official model.
- **Validation.** Mouse-grouped, age-stratified 5-fold outer CV with inner validation for calibration/early stopping/threshold — cells from the same mouse never split. Donor-level pooled out-of-fold predictions → donor-level ROC-AUC.
- **Sensitivity.** Applied frozen Limb-FACS model to independently processed TMS Droplet Limb MSC.
- **External.** GSE176206 partial-reprogramming Adipo + MSC — Young/Aged controls and Aged SOKM.
- **Scoring interface.** `.h5ad` raw counts in → per-cell `youth_score`, `predicted_age_months`, `gene_overlap`, `qc_status` out. Coverage QC gate at 0.70.

## 2. What is genuinely strong

1. **Model class comparison is properly structured.** Three model families run through the *same* mouse-grouped stratified 5-fold protocol, judged on donor-level AUC, with a Technical-only baseline as a check that the model isn't just riding sequencing depth. This is the right experimental design.
2. **Technical-only baseline is a big deal.** AUCs of 0.696 (SCAT) and 0.688 (Limb) from depth + detected genes + sex alone quantify how much of the "youth signal" is confoundable. That the formal models beat this by 0.20–0.30 AUC is what makes the AUCs interpretable.
3. **Honest Transformer read.** "Learned an age signal but did not outperform simpler models" is exactly the correct conclusion at ~15 donors — no p-hacking, no overselling. That is a hard call to make when you built the fancier thing.
4. **External CI interpretation is precise.** The Session 4 note explicitly says +0.030 (95% [+0.0135, +0.0453]) is the uncertainty of the *estimate itself*, not a between-Adipo-and-MSC contrast. This distinction is routinely muddled in this literature; getting it right at this stage is very good.
5. **Cross-modality direction preserved.** ρ = −0.856 at the donor level from FACS-trained Limb model on Droplet Limb is meaningful *directional* evidence that the model isn't tied to the training platform.
6. **Correct framing of the SOKM negative result.** "Our model didn't detect an SOKM-driven trend" ≠ "SOKM is harmful" — the report keeps these separate.

## 3. What is fragile — in decreasing order of impact

### (a) AUC = 1.000 for the Limb Gene signature is a warning sign, not a strength

Perfect ranking of 14 held-out mice is possible with a genuinely strong signal at this cohort size, but a report should treat it as a **saturation ceiling**, not evidence that the model is uniformly excellent. Two things to add before shipping:

1. **Gene-set stability across the 5 folds** (Jaccard between fold-selected signatures) — if the "best gene signature" is very different in each fold but each fold happens to work, that's overfitting-to-donor-batch.
2. **Which donors are ranked adjacent to each other** — check whether the AUC = 1 folds are separating age but not, say, plate/library/sex.

### (b) Bootstrap CIs on 3 donors per external condition are architecturally weak

In GSE176206, MSC has 3 donors per condition. A donor-level bootstrap resamples 3 out of 3 with replacement — the number of distinct bootstrap draws is small (C(3+3−1,3) = 10 combinations before considering permutations), so the 95% donor-bootstrap interval [+0.0135, +0.0453] is optimistic. **Report the interval alongside a Wilcoxon or permutation p-value at n=3 vs n=3**, or a min–max range across the 3! = 6 donor-swap resamples, so the reader can see the small-N ceiling. This is not a fatal issue — the direction is right — but the current CI implies more precision than 3 donors support.

### (c) `predicted_age_months` and `youth_score` come from different models

The spec says `predicted_age_months` **always** comes from the Transformer, while `youth_score` may come from Elastic Net (SCAT) or Gene signature (Limb). Two consequences:

1. A user reading a CSV where `youth_score = 0.9` and `predicted_age_months = 20` gets a contradictory signal from a model committee that never actually voted together.
2. If the Transformer was rejected on AUC grounds, why is its auxiliary head still shipped as a first-class output? **Either use the same model class for both outputs**, drop `predicted_age_months` from the primary release, or explicitly label it "auxiliary — do not use as primary readout" in the column name.

### (d) Binary Young/Old collapses the age axis unnecessarily

TMS has 1, 3, 6, 12, 18, 21, 24, 30 mo mice for Droplet Limb (per Zihan's droplet trajectory table); Kaile's FACS subset was collapsed to Young=3m vs Old=18/24m. Two costs: middle ages carry no label information, and the "youth score" trained on the extremes will be poorly calibrated in the middle where it is most biologically interesting. **Consider training a continuous age regression** (elastic-net Cox or Gaussian regression, or ordinal logistic) — the same signature/EN framework accepts a continuous target with no architectural change.

### (e) 70% gene-overlap QC threshold is loose for a signature-based model

Missing 30% of the model's genes is a lot when the model is a weighted 100-gene score. The signature is not robust to arbitrary 30% loss — some genes carry more weight than others. **Compute the weighted-coverage metric that Zihan's parser already uses** (sum of `|w_g|` present / sum of all `|w_g|`) and gate on that too, at ≥ 0.85.

### (f) Cross-modality Droplet trajectory has a 24m vs 30m inversion

The reported means are 1m = 0.448, 3m = 0.273, 18m = 0.215, 21m = 0.162, 24m = 0.112, 30m = 0.142. The 1m→24m stretch is monotonic-decreasing as expected, but 30m scores *higher* than 24m. The report's own text acknowledges this ("30m is slightly higher than 24m and the trajectory is not strictly monotonic"), and with only 2 mice per age this may just be donor noise. Still: with 2 mice at 30m the "monotonic aging axis" claim can only be defended if the 30m point comes with a mouse-level CI that overlaps 24m. **Either attach a donor-level uncertainty bar at each age or downweight 30m in the summary trajectory rho.** More broadly, the 1m point being *higher* than the 3m training-anchor young is worth flagging as evidence that the model extrapolates outside its training range and should not be interpreted as a probability at ages below 3m.

### (g) Reproducibility surface is very thin

The zip contains three markdown files; no code, no model files, no example input, no versions, no seed. Compared to Zihan's repo (17 scripts + tables + figures + models + parser + calibration JSON), Kaile's deliverable is not independently runnable. This is the most consequential gap for review. **Ship at minimum: (i) the signature genes and weights CSV for both SCAT and Limb, (ii) the Elastic Net coefficients, (iii) fold-level CV predictions for reproducing the AUC, (iv) the scoring script.** The GitHub repo pattern Zihan used is a good template.

### (h) Transformer specification is missing

No parameter count, no learning rate, no epochs, no early-stopping criterion, no seed, no loss function is documented. Given the Transformer is a headline comparison and the *negative* conclusion ("did not outperform simpler models") is being reported, the Transformer configuration is exactly what a reader needs to trust that negative. **Document the Transformer as thoroughly as the Gene signature is documented.**

### (i) Ensemble-of-5-folds "official" model is ambiguous

Section 3 of the input/output doc says the scorer "averages predictions across the five cross-validation folds". That means the deployed model is a bag of 5 models, each with different HVGs and calibration. This is reasonable, but then the folds themselves are the model — losing any one fold changes the score. **Ship all five fold artifacts**, and consider retraining a single final model on the full cohort (with calibration held out) for a single deployable object.

### (j) SCAT model applied to Adipo GSE176206 is a documented cross-cell-state check, but the Adipo n=1 pooled library makes it non-inferential

The report says this. The recommended presentation is to keep Adipo out of any bar/CI plot entirely and only mention it as a text-level directional check — otherwise readers will read the +0.050 as a real point estimate with power behind it.

## 4. Smaller notes

- The Session 4 update references three images pinned to Windows local paths (`C:\Users\HUAWEI\...`) — will render as broken in any shared markdown viewer. Move images into a `figures/` subfolder next to the markdown and use relative paths, matching what Zihan did.
- "Youth Score is defined as fold-calibrated P(young)" — say which calibrator (Platt vs isotonic vs raw sigmoid) and show a reliability diagram at least for the Elastic Net on SCAT.
- Coverage metric in the output CSV is "minimum fraction across the five folds" — this is a conservative reading; also report the *mean* fold coverage so users see whether one fold is dragging the QC down.
- The auxiliary Transformer's 64 expression bins per gene × up to 4,096 HVGs × identity embeddings is a large model for 15 donors. It's fine as an experiment but should not be quietly retained as an official output surface after being rejected on primary AUC.

## 5. Recommended next steps

1. **Ship code, model weights, and fold-level predictions** as a repo (Zihan's `deliverables/…` layout is a good template).
2. **Report per-fold gene-set overlap** for the Limb Gene signature to sanity-check the AUC = 1.000.
3. **Reframe the external CI** — pair the donor-bootstrap CI with a permutation p at n=3 vs n=3 so the small-N ceiling is visible.
4. **Either drop `predicted_age_months` from the primary output or align it to the official model class** for the tissue.
5. **Attach mouse-level CIs to the Droplet age trajectory** so the 24m/30m non-monotonicity is either quantified or explained, and flag 1m as extrapolation below the 3m training anchor.
6. **Add continuous age regression** as a comparator against the binary Young/Old model.

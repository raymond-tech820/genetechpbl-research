#!/usr/bin/env python3
"""Step 10 candidate signature compression for the Limb_Muscle MSC Youth Score."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


STABILITY_PATH = Path("outputs/stability/gene_reliability_scores.csv")
OUT_DIR = Path("outputs/signature")
MODELS_DIR = Path("models")

TOP_N_PER_DIRECTION = 100
CORE_N_PER_DIRECTION = 50


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    MODELS_DIR.mkdir(parents=True, exist_ok=True)

    rel = pd.read_csv(STABILITY_PATH)
    required = {
        "gene",
        "adjusted_logFC",
        "age_rho",
        "pi_LOMO",
        "pi_depth",
        "pi_sex",
        "passes_step9_initial_reliability",
        "continuous_age_trend_compatible",
    }
    missing = required.difference(rel.columns)
    if missing:
        raise ValueError(f"Missing required columns in {STABILITY_PATH}: {sorted(missing)}")

    rel = rel.copy()
    rel["direction"] = rel["adjusted_logFC"].map(lambda x: "young_high" if x > 0 else "old_high")
    rel["r_g"] = rel["pi_LOMO"] * rel["pi_depth"] * rel["pi_sex"]
    rel["q_g"] = rel["adjusted_logFC"].abs() * rel["age_rho"].abs() * rel["r_g"]
    rel["abs_logFC"] = rel["adjusted_logFC"].abs()
    rel["abs_age_rho"] = rel["age_rho"].abs()
    rel["candidate_pool"] = rel["passes_step9_initial_reliability"].astype(bool)

    candidate = rel[rel["candidate_pool"]].copy()
    candidate = candidate.sort_values(
        ["direction", "q_g", "abs_logFC", "abs_age_rho", "gene"],
        ascending=[True, False, False, False, True],
    )

    young = candidate[candidate["direction"] == "young_high"].sort_values(
        ["q_g", "abs_logFC", "abs_age_rho", "gene"],
        ascending=[False, False, False, True],
    )
    old = candidate[candidate["direction"] == "old_high"].sort_values(
        ["q_g", "abs_logFC", "abs_age_rho", "gene"],
        ascending=[False, False, False, True],
    )

    ranking_cols = [
        "gene",
        "direction",
        "adjusted_logFC",
        "FDR",
        "age_rho",
        "abs_logFC",
        "abs_age_rho",
        "pi_LOMO",
        "pi_depth",
        "pi_sex",
        "r_g",
        "q_g",
        "LOMO_sign_rate",
        "max_LOMO_delta",
        "low_depth_abs_effect_size_change",
        "sex_effect_old_male_minus_old_female_logcpm",
        "sex_effect_ratio_vs_age",
        "stability_weighted_abs_logFC",
    ]
    existing_ranking_cols = [col for col in ranking_cols if col in rel.columns]

    candidate[existing_ranking_cols].to_csv(OUT_DIR / "step10_ranked_reliable_candidates.csv", index=False)
    young[existing_ranking_cols].to_csv(OUT_DIR / "step10_young_high_ranked_candidates.csv", index=False)
    old[existing_ranking_cols].to_csv(OUT_DIR / "step10_old_high_ranked_candidates.csv", index=False)

    top_young = young.head(TOP_N_PER_DIRECTION).copy()
    top_old = old.head(TOP_N_PER_DIRECTION).copy()
    compressed = pd.concat([top_young, top_old], ignore_index=True)
    compressed["selected_top_n_per_direction"] = TOP_N_PER_DIRECTION
    compressed[existing_ranking_cols + ["selected_top_n_per_direction"]].to_csv(
        OUT_DIR / "step10_compressed_signature_candidates_top100_each.csv",
        index=False,
    )

    core = pd.concat([young.head(CORE_N_PER_DIRECTION), old.head(CORE_N_PER_DIRECTION)], ignore_index=True)
    core["selected_core_n_per_direction"] = CORE_N_PER_DIRECTION
    core[existing_ranking_cols + ["selected_core_n_per_direction"]].to_csv(
        OUT_DIR / "step10_core_signature_candidates_top50_each.csv",
        index=False,
    )

    model_signature = core.copy()
    model_signature["base_weight_abs_logFC"] = model_signature["abs_logFC"]
    model_signature["signed_weight"] = model_signature["adjusted_logFC"].clip(-3, 3) * model_signature["r_g"]
    model_signature["module"] = model_signature["direction"].map(
        {"young_high": "young_module", "old_high": "old_module"}
    )
    model_cols = [
        "gene",
        "module",
        "direction",
        "adjusted_logFC",
        "age_rho",
        "r_g",
        "q_g",
        "signed_weight",
        "FDR",
        "pi_LOMO",
        "pi_depth",
        "pi_sex",
    ]
    model_signature[model_cols].to_csv(
        MODELS_DIR / "limb_msc_general_youth_score_v1_signature_candidates_step10.csv",
        index=False,
    )

    summary = pd.DataFrame(
        [
            {
                "total_genes_in_reliability_table": len(rel),
                "genes_passing_step9_reliability": len(candidate),
                "young_high_passing": len(young),
                "old_high_passing": len(old),
                "top_n_per_direction_exported": TOP_N_PER_DIRECTION,
                "core_n_per_direction_exported": CORE_N_PER_DIRECTION,
                "compressed_top100_total": len(compressed),
                "core_top50_total": len(core),
                "q_g_formula": "|adjusted_logFC| * |age_rho| * r_g",
                "r_g_formula": "pi_LOMO * pi_depth * pi_sex",
            }
        ]
    )
    summary.to_csv(OUT_DIR / "step10_signature_compression_summary.csv", index=False)

    def fmt_top(frame: pd.DataFrame, n: int = 15) -> str:
        top = frame.head(n)[["gene", "adjusted_logFC", "age_rho", "r_g", "q_g", "FDR"]].copy()
        top[["adjusted_logFC", "age_rho", "r_g", "q_g", "FDR"]] = top[
            ["adjusted_logFC", "age_rho", "r_g", "q_g", "FDR"]
        ].round(4)
        return top.to_markdown(index=False)

    report = f"""# Step 10: Candidate Signature Compression

## Inputs

- Reliability table: `{STABILITY_PATH}`
- Guidance basis: Step 9 instruction plus `outputs/eda/guidance_review_for_next_steps.md`

## What This Step Did

1. Started from genes passing Step 9 initial reliability.
2. Computed the combined reliability score:

```text
r_g = pi_LOMO * pi_depth * pi_sex
```

3. Computed the candidate ranking score:

```text
q_g = |adjusted_logFC| * |age_rho| * r_g
```

4. Split candidates into `young_high` and `old_high` using the sign of the sex-adjusted age logFC.
5. Ranked each direction separately by `q_g`.
6. Exported top 100 per direction as a compressed candidate pool.
7. Exported top 50 per direction as a smaller core candidate pool for the first Youth Score construction pass.

## Summary

- Genes in reliability table: {len(rel)}
- Genes passing Step 9 reliability: {len(candidate)}
- Young-high passing genes: {len(young)}
- Old-high passing genes: {len(old)}
- Top candidates exported per direction: {TOP_N_PER_DIRECTION}
- Core candidates exported per direction: {CORE_N_PER_DIRECTION}

## Top Young-High Candidates

{fmt_top(young)}

## Top Old-High Candidates

{fmt_top(old)}

## Important Interpretation

This step does not yet define the final Youth Score. It compresses the candidate space into ranked young-high and old-high modules using sex-adjusted effect size, continuous-age trend strength, and robustness from Step 9. Final acceptance still requires score-level validation against age, old-sex separation, library size, cell count, leave-one-mouse-out score stability, and permutation/random controls.

## Outputs

- `outputs/signature/step10_ranked_reliable_candidates.csv`
- `outputs/signature/step10_young_high_ranked_candidates.csv`
- `outputs/signature/step10_old_high_ranked_candidates.csv`
- `outputs/signature/step10_compressed_signature_candidates_top100_each.csv`
- `outputs/signature/step10_core_signature_candidates_top50_each.csv`
- `outputs/signature/step10_signature_compression_summary.csv`
- `models/limb_msc_general_youth_score_v1_signature_candidates_step10.csv`
"""
    (OUT_DIR / "step10_signature_compression_report.md").write_text(report, encoding="utf-8")

    print("Step 10 complete")
    print(summary.to_string(index=False))
    print("\nTop young-high:")
    print(young[["gene", "adjusted_logFC", "age_rho", "r_g", "q_g", "FDR"]].head(10).to_string(index=False))
    print("\nTop old-high:")
    print(old[["gene", "adjusted_logFC", "age_rho", "r_g", "q_g", "FDR"]].head(10).to_string(index=False))


if __name__ == "__main__":
    main()

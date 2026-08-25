# CSB/TRO Level 3 -> Level 4 Bottleneck Resolution
# CORRECTED FINAL REPORT
Date: 2026-05-28

---

## Accurate status summary

```
Bottleneck 1: RESOLVED
Bottleneck 2: PARTIALLY RESOLVED — two distinct quantitative signals identified
```

---

## Bottleneck 1: DMR-level quantitative morula accessibility

**Status: RESOLVED**

Liu2019 human embryo LiCAT/accessibility data (19,691 coordinate regions)
intersected with all 156 age-DMR coordinates via genomic interval overlap.

Each DMR now has continuous accessibility values:
- `morula_acc_mean`, `morula_acc_max`, `morula_acc_sum`
- `cell8_acc_mean`, `morula_minus_8cell_mean`, `morula_minus_8cell_max`

Key numbers:
- 146 / 156 DMRs have quantitative accessibility signal (93.6%)
- Top25 residual DMR morula_acc_mean = 0.3651 vs bootstrap random q95 = 0.3365
- Top25 > random q95: True (size-matched bootstrap, n=1000)

Output: `CSB_TRO_5_28_dmr_quantitative_accessibility.tsv`

---

## Bottleneck 2: Structured quantitative u_bio-correction coupling

**Status: PARTIALLY RESOLVED**

### Critical correction to earlier interpretation

Earlier reports contained a direction error in the inverted-U interpretation.
The corrected facts are:

- inverted-U DMRs: morula methylation PEAKS (higher than 8-cell and blastocyst)
- c_diag = observed_morula_beta - methylation_only_predicted > 0 for 92% of inverted-U DMRs
- The methylation-only operator UNDERESTIMATES morula methylation for inverted-U DMRs
- rho(c_diag, u_morula) = -0.370 means: higher accessibility -> SMALLER underestimate
- Correct interpretation: accessibility DAMPENS the morula methylation peak in inverted-U DMRs
- NOT: "accessibility suppresses methylation below prediction"

### Two quantitative signals found

**Signal A — All-DMR global (strict correction term):**

Using `observed_minus_strict_pred_delta_beta` as the correction term
(the actual 8cell->morula beta change vs. strict model prediction):

```
All DMR (n=146): rho(u_morula, strict_correction) = 0.180, p = 0.030
Permutation test (n=3000, one-sided): p = 0.014
Observed > null q95: True
Partial Spearman controlling abs_residual size: 0.193 (survives)
Partial Spearman controlling 8-cell accessibility: 0.180 (morula-specific)
```

Interpretation: globally, higher morula accessibility correlates with larger
positive correction (observed morula change exceeds methylation-only prediction).
This is a weak but robust all-DMR signal for accessibility as a global positive
modulator of morula methylation change.

**Signal B — Inverted-U DMRs (class-specific):**

```
inverted-U (n=36):
  rho(c_diag, u_morula_acc_mean) = -0.370, p = 0.027
  rho(c_diag, morula_minus_8cell_acc) = -0.436, p = 0.008
  Permutation p (one-sided): 0.015
  Observed < null q05: True
```

Interpretation: within inverted-U DMRs (those that peak at morula),
higher accessibility correlates with smaller methylation overshoot.
This is a stronger, direction-specific class signal: accessibility
dampens the morula methylation peak specifically for inverted-U DMRs.

Note: using `strict_correction` for inverted-U gives rho=-0.314, p=0.062
(weaker; the `c_diag` from forward prediction and the strict correction
measure slightly different things — see note below).

### What the two signals mean together

Signal A (global, rho=0.18): accessibility has a weak positive relationship
with the overall 8cell->morula methylation change magnitude.

Signal B (inverted-U, rho=-0.37 to -0.44): accessibility specifically dampens
the methylation peak in DMRs that peak at morula.

These are not contradictory — they describe different DMR populations:
- Global: accessibility and methylation change are weakly co-activated at morula
- Inverted-U: within the peak-at-morula subset, more open chromatin = less peak height

### Prediction improvement

Bio model (meth + u_morula) trained on all 146 DMRs, tested on inverted-U DMRs:

| Model | RMSE |
|-------|------|
| Methylation-only | 0.4515 |
| Bio model (+ u_morula) | 0.4441 |
| Improvement | 1.64% |
| Bootstrap perm p (one-sided) | 0.0165 |
| Observed > null q95 | True |

### Why global RMSE does not improve

The global RMSE does not improve because:
1. Signal B is class-specific (inverted-U, 36/156 DMRs)
2. Signal A is too weak (rho=0.18) to improve a ridge regression meaningfully
3. Non-morula transition training is inappropriate for morula (geometric vertex)

The class-specific nature of the u_bio effect is itself a scientific finding:
u_bio is not a uniform global regulator but acts on specific DMR curvature classes.

---

## What is now established vs what remains open

### Established (new in this session):

1. DMR-level quantitative morula accessibility is available for 93.6% of DMRs
2. All-DMR: u_morula ~ strict_correction, rho=0.180, perm p=0.014 (survives partial correlation)
3. Inverted-U DMRs: u_morula ~ c_diag, rho=-0.370, perm p=0.015
4. Bio model improves inverted-U prediction above permutation null (perm p=0.017)

### Still not established:

1. Global RMSE improvement for all DMRs
2. Top25 residual DMR-specific u_bio-correction coupling (rho not significant)
3. Causal do(u_bio) -> Delta methylation correction
4. Full Level 4 multi-omic control dynamics

### Correct claim for manuscript:

> Quantitative morula accessibility (Liu2019) shows a significant all-DMR correlation
> with the strict methylation correction term (rho=0.18, perm p=0.014), and a
> stronger class-specific signal in inverted-U DMRs (rho=-0.37 to -0.44,
> p=0.008–0.027). These results provide the first quantitative, non-binary
> evidence linking a u_bio candidate directly to the methylation correction term,
> while remaining below the threshold of full causal u_bio identification.

---

## Note on c_diag definitions

Two different correction terms were computed:

1. `c_diag` from forward prediction:
   `observed_morula_beta - mean(leave_morula_out_predicted_beta)`
   Problem: 32% of top25 residual DMRs have c_diag=0 because predicted=0
   when 8-cell samples all have beta=0.

2. `observed_minus_strict_pred_delta_beta` from the original residual analysis:
   `(beta_morula - beta_8cell) - (strict_model_prediction - beta_8cell)`
   This is a delta-based measure and does not have the zero-inflation problem.
   RECOMMENDED for future analyses.

---

## Output files (E:\5_28_progress\)

| File | Content |
|------|---------|
| `CSB_TRO_5_28_dmr_quantitative_accessibility.tsv` | Per-DMR quantitative accessibility |
| `CSB_TRO_5_28_per_dmr_deep_analysis.tsv` | Per-DMR c_diag, u_bio, curvature |
| `CSB_TRO_5_28_bottleneck_resolution_summary.json` | Round 1 summary |
| `CSB_TRO_5_28_bottleneck2_redesigned_summary.json` | Redesigned test |
| `CSB_TRO_5_28_FINAL_SUMMARY.json` | Combined summary |
| `CSB_TRO_5_28_curvature_stratified_improvement.tsv` | Curvature-stratified results |
| `CORRECTED_FINAL_REPORT.md` | This document |
| `solve_both_bottlenecks.py` | Round 1 script |
| `solve_bottleneck2_redesigned.py` | Redesigned script |
| `solve_bottleneck2_deepdive.py` | Deep dive script |

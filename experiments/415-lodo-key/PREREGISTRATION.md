# 415 — preregistration (written 2026-09-22, before any LODO-key number was computed)

## Question
410's key `k410_p = mean_{i<j} ½(<W d_i,d_j>+<W d_j,d_i>) / |W m|²` estimates the mse-optimal per-target scale of the
whitened donor BASE. The shipped bulk row also carries a CD4 residual and an H1-union residual. Does a leave-one-donor-out
key that estimates `<P_p,T_p>` for the FULL production-composition row P rank/scale targets better than 410's key, on the
same production-composition fold predictions?

## Objects (PLOLO, 162's cached folds; three cancer folds + kolf2 reported separately)
- `P` = 228's `build_final(m=None, a=1)`: `whitened = Operator.apply(ALPHA·pre)` (γ 0.65, λ 0, fold operator with its
  on-axis preserve), `+ merged(frac_cd4·CD4, frac_u·H1union)` with CD4 f = 0.34 (all-rows sizing) and H1 union
  (`gse295214_h1`, `h1_2025`) g = 0.28 (carrying sizing), production's post-add preserve, panel coordinates zeroed.
- `T` = the held library's truth row. Everything is evaluated on the fold's `scored` columns, `live` targets
  (non-zero P and T rows), exactly as check410.
- Oracle `s*_p = <P_p,T_p>/|P_p|²`.

## Keys compared (all read only source-donor rows, CD4 and H1 residual rows; never T)
- **K410**: `key410.shared_energy_key(targets, fold sources, fold op, wcols, cols=scored)` — unchanged from 410.
- **LODO-V1 (PRIMARY, the mandate's literal construction)**: `key_p = mean_j <P^{−j}_p, d_{j,p}> / |P_p|²`, mean over
  source donors j carrying p; `d_j` = donor j's raw bank row (410's `donor_rows`, uncentred); `P^{−j}` = P rebuilt
  linearly without donor j: the SAME global constants as P (on-axis preserve c, frac_cd4, frac_u, post-add preserve s),
  base = `c·W(ALPHA·pre^{−j})` where `pre^{−j}` is the mean of the other carrying sources' prepared rows (the rows
  `build_baseline` averages; k562 PM-blended as in build_fold). Linear LOO with fixed constants is the unbiased
  construction under donor exchangeability (E<P^{−j},d_j> = E<P,T>).
- LODO-V2 (secondary): as V1 but `d_j` = donor j's prepared row (centred, ×EFFECT_SCALE/raw_median; the row that enters
  `pre`), so donors mix on a common scale and d_j matches the truth's centring.
- LODO-V3 (secondary): V2 with P^{−j} fully rebuilt through the production pipeline on the reduced donor set
  (preserve constants and residual fractions recomputed per j).
- Decomposition row: LODO-V2 with the residual term dropped (base only) — isolates the residual-coverage effect from the
  row-preprocessing effect.
- Carried set for every LODO key = targets carried by ≥ 2 source donors (identical to K410's), non-carried get w = 1.

## Endpoints (per fold)
1. Spearman(key, s*) over carried targets.
2. Relative aggregate transfer cosine gain `agg(w)/agg(1) − 1`, `agg(w) = Σ w_p<P_p,T_p> / sqrt(Σ w_p²|P_p|² · |T|²)`,
   `w = clip(key/median(key over carried), lo, hi)`, at clip(0.3, 3) and clip(0.3, ∞). (Energy-preserving normalisation
   is a global constant and does not change the cosine.)
3. 1,000× target bootstrap (seed 415): CI of each gain and of the PAIRED difference LODO − K410.
4. Negative control: the LODO-V1 weights permuted among carried targets (mean of 200 permutations).

## Known answers (run aborts if any fails)
- KA-1: K410 on the ORIGINAL fold prediction `Operator.apply(pre, wcols, 0.65, 0)` reproduces check410.json
  (Spearman and gains at clip(0.3,3), clip(0.3,∞) on all four folds to ≤ 1e-9).
- KA-2: `build_baseline` of the re-prepared source rows equals the cached fold `pre` (float32 cache, ≤ 1e-5 relative).
- KA-3: the production composition reproduces 228's KA rows `base222` (a23 cos, pds; 223 focus223) to ≤ 1e-12.
- KA-4 (algebraic): on the original prediction with no residual and 2-source folds, LODO with d = the rows W acts on is
  identical to the 410 formula applied to those rows (checked numerically on the prepared rows, ≤ 1e-9 relative).

## Decision rule (fixed now)
- "Beats 410": LODO-V1 gain at clip(0.3,3) > K410 gain at clip(0.3,3) on ALL THREE cancer folds (point estimates),
  both evaluated on the production-composition P. Then — and only then — the A/B/C panel artifact is written for
  clip(0.3,3), and for clip(0.3,∞) only if V1 also wins at ∞ on all three cancer folds.
- Secondary variants are reported, not selected: if only a secondary wins, the recommendation is a preregistered
  follow-up, not an artifact.
- kolf2 is reported separately (different lineage; production operator) and does not enter the rule.

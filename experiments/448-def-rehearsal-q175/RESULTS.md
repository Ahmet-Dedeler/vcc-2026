# 448 — October rehearsal of the current recipe at q 1.75 on the essential-free, realistic-coverage panel (efr)

**Date:** 2026-09-25 00:5x–02:31Z (coordinator). Builder `b6c3f131` (447's widened bound), driver `scripts/def_pipeline.py`,
bundle `experiments/238-essential-free-panel/work/ef-bundle-r`. No scoring, no submission. Log: `work/run448.log`.

## Why
447 widened `VCC_BULK_CTX_SCALE`'s bound to [0.7, 1.4] because 438's level file asked for factors above 1.2 on this panel
at q ≥ 1.5. Nobody had built and re-levelled a D/E/F-like panel end to end at a factor above 1.2. This does that, with the
recipe the board now carries (433k + q 1.75 = 446b) and 438's October stack default (H1 ×2, CD4 ×1).

## What ran
1. Driver setup stages (preflight → identity → CD4 extract → census → the day's whitener → the day's 410 key → bare-env
   pins). All rc 0.
2. Stack dump H1 ×2 / CD4 ×1 on the day's operator: `work/efr/dump_h1x2_cd4x1.npz` sha `16cdf3d3e3d2e8f8`.
3. lvl0: the un-levelled LOAD build at q 1.75 (with SHIPPED_ENV's coarse MDE and zero repair). `ntclevel417.py --sd 0.01`.
4. lvl1: rebuilt at ntclevel's factors. ntclevel again.
5. lvl2: one correction pass (fixed exponent 1.86 in `run448.sh`). ntclevel again.
6. Builder restored by the trap: `b6c3f131` (identical).

## Result

| build | factors D / E / F | ρ D | ρ E | ρ F |
|---|---|---:|---:|---:|
| lvl0 | 1 / 1 / 1 | 0.6706 | 0.6206 | 0.5911 |
| lvl1 | 1.2084 / 1.2562 / 1.2871 | 0.9595 | 0.9566 | **1.0568** (past the cliff) |
| lvl2 | 1.2218 / 1.2721 / 1.2356 | **0.9798** | **0.9798** | 0.9463 |
| **lvl3** | 1.2218 / 1.2721 / **1.2514** (secant) | 0.9798 | 0.9798 | **0.9784** |

EV target ρ* = 0.9793 (sd 0.01, q = 0.31 / 2.09). Audits: zero-gene 0 on lvl0 and lvl1; `mde_lambda` median 0.951 on lvl1
(A/B/C's 433a build: 0.956 / 1 / 1), so the coarse MDE delivers on this panel.

1. **The bound change was needed and works.** Every factor the panel needs is in [1.2, 1.29]; the old [0.8, 1.2] would
   have refused all three contexts at q 1.75. The build path at factors > 1.2 is clean (no refusals, zero audit 0).
2. **One ntclevel pass is not enough on D/E/F.** ntclevel converts its claim into a factor with a square root (ρ ∝ f²).
   The render is not quadratic: ρ ∝ f^e with e = 1.89 / 1.90 / 2.30 between f = 1 and lvl1. D and E came in 2 % short.
   F went 5.7 % past its cliff: shipping lvl1 would have put F above its cliff, which is the expensive side of the
   asymmetric loss (1.78 vs 0.31 per unit).
3. **F is expansive near the top.** Between lvl1 and lvl2 F's local exponent is 2.70. The fixed-exponent correction
   (1.86) put F at 0.946, 3.4 % below target. That's the safe side (≈ 0.0003 overall left on the table by the cliff model),
   but it is not converged.
4. **Per-context secant steps fix it.** `relevel_iter.py` (rewritten): for each context, take the two builds nearest ρ*,
   bracketing when possible, read the local log-log exponent and solve for ρ*. From lvl0 + lvl1 it retrodicts lvl2's D
   and E to ±0.0001 in ρ (0.9797 predicted, 0.9798 read). From all three it declares D and E converged and puts F at
   1.2514 (bracketed between 1.2356 and 1.2871).

5. **lvl3 (03:28Z, `run448b.sh`): one bracketing secant step lands F at ρ 0.9784** (target 0.9793, −0.09 %), inside the
   1 % tolerance. D and E, rebuilt at the same factors, read the same ρ to 4 decimals (the render is deterministic).
   `relevel_iter.py` on all four builds prints CONVERGED. The loop took two levelled builds for D/E and three for F.

## For the day (written into DEF-RUNBOOK §5.-1)
- Budget 3–4 builds for the re-level, about 25 min each on this machine.
- The loop: lvl0 (unit factors) → ntclevel → lvl1 → ntclevel → `relevel_iter.py` → lvl2 → … Stop when every context is
  within 1 % of ρ* in ρ (0.5 % in scale) and none is above 1.
- Never ship a build with any ρ > 1. If the loop has to stop early, ship the last build whose contexts are all below
  their cliffs.

## Files
- `run448.sh` (lvl0–lvl2), `run448b.sh` (lvl3), `relevel_iter.py` (the secant step).
- `work/efr/448-lvl{0,1,2}.ntclevel.json`, `work/efr/*.audit428.json`, `work/efr/rc_audit_lvl1.jsonl`.
- The h5ad builds were deleted after their reads (rehearsal only; nothing here is a submission).

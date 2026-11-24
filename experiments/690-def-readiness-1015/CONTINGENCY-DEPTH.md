# Depth contingency: D/E/F released at another depth (690, 2026-10-15)

Class: day-of procedure. No recipe value changes. Coordinator's question: if D/E/F arrive at full depth (×1.47, Arc's
29,482 UMI/cell) instead of A/B/C's ~20k, what does the day do, and is `DEF_DAY_ALLOW_DEPTH=1` safe?

## Likelihood
Low. Arc downsampled every released A/B/C context to one ~20k median (20,109 / 19,946 / 20,034; 641), and the final
bundle is served by the same download handler (`getFinal=true`, 680's NEWS). It still has to be planned, because 641's
row turns the case into a preflight FAIL, and a FAIL stops AUTO_R0.

## What changes with depth, and what does not
- **Refitted from the day's own controls, so depth-adaptive:** every per-gene quantity (control means, CPM gates,
  b_ctrl), the NTC-proxy re-level target ρ* (ntclevel417 runs on the day's controls), library sizes. The emitter draws
  each cell's library from the controls, so predictions come out at the controls' depth, which is the truth's depth.
- **Fixed in CPM, so scale-free:** `--min-cpm 10`, the scorer's own > 5 CPM gene filter, the call-set CPM block (K 2,200).
- **Depth-sensitive, and the real risk:** statistical power. At ×1.47 counts per cell, Wilcoxon calls more genes in
  BOTH the truth (n_real grows) and our cells. fid (sign × min(1, n_pred / n_real)), jac and reach depend on how those two
  grow together, and the call-set size K was tuned at 20k (521: K is regime-sensitive, +0.0048 at H1 f 0.35, −0.015 at
  f 1). The MDE scales ∝ depth^−1/2 (−18 % at ×1.47), beyond 554's 10 % dose tolerance in the conservative direction.
- **Downsampling our controls to 20k is worse:** the predictions would be emitted at 20k while the truth sits at ×1.47.
  Per-cell noise would mismatch and DE power fall below the truth's, on every member.

## Decision for the day (proposed; the coordinator decides)
1. **Uniform shift** (all three medians agree within 2 %, ratio in [0.6, 1.8]): build at native depth with
   `DEF_DAY_ALLOW_DEPTH=1`. R0 is not the last entry, so an unpriced K costs at most R0's DE members; a priced K shift
   (step 3) can follow in a later entry, one variable.
   - 691 does this unattended only if armed with `AUTO_ALLOW_DEPTH691=1`. `release691.py overrides` grants exactly this
     case and nothing else. Preflight re-runs with the override and must pass. NEXT.txt and the log record it.
   - Tested: a uniformly thinned def-bundle (×0.80 on D/E/F) went preflight FAIL → override → preflight rc 0 → AUTO_R0
     gate (`691/work/testrecords/v3_depth`). The non-uniform 690 mock (0.84 / 1.01 / 0.80) is refused.
2. **One context off alone**: STOP, unchanged (a damaged or mislabelled file). The helper refuses it.
3. **Price K at the new depth before the last entry.** An instrument exists: the 2025 H1 training file
   (`data/external/vcc2025/train/adata_Training.h5ad`, median **54,208 UMI**, q25/q75 42.9k/67.6k, measured 10-15) is at
   full depth, while our h1-test harness is downsampled to ~20k (19,968). Downsample H1 (truth and controls, binomially) to
   20k and to 20k × f. Then price the recipe's DE on the real emitter (run496, with production's zero repair) at K 2,200
   and one or two alternatives at both depths, with cell-eval2 (cloud: `cloud_score.sh`, ~$1.50 a run). It only runs if
   the case happens, after R0, in a machine window.
4. **Spread** (554, q75/q25 > 3.0): the card's own rule. ≤ 3.5 in every context → `DEF_DAY_ALLOW_WIDE_LIBS=1` (armed by
   `AUTO_ALLOW_WIDE691=1`); above 3.5 anywhere → STOP. 554's postbuild dose guard still checks the build.

## Arming
Both flags are OFF in the armed watcher (pid in PLAN Notes). Arming them is the coordinator's decision: re-launch with
`AUTO_R0=1 AUTO_ALLOW_DEPTH691=1 AUTO_ALLOW_WIDE691=1` (the watcher runs from a frozen copy, so a kill + re-arm is safe).

# 477 — 480's donor-consensus rider through 490's exact emitter reconstruction: nmae −0.003 reproduces, but the mixed form fails the preregistered rule. Not GO.

2026-09-26 23:40Z, third session (at the coordinator's request). Analytic only (no emitter hook, no build, no
`cell-eval2` CLI; cell-eval2's own `prepare_de` + `dispatch_de_metrics` through 490's `score_de.py`). Preregistered in
`PREREGISTRATION.md` (commit `e2413582`). Script `price477.py`, output `work/price477.{json,log}`. 12 s, peak RSS 1.2 GB.

## Known answers (all pass)
- 490's per-target reconstruction known answers (kf, mde_fu/fd, mde_n_up/dn, flip, n_block == each run's diag) on
  85/85 targets × 3 seeds.
- The delivered tables reproduce the harness's metrics JSON, |d| ≤ 1e-9.
- k 0 returns Δ = 0 exactly on every member.

## Result (H1, 85 targets; Δ vs the shipped coarse-MDE construction; overall-equiv with 475's slopes; 3,000-draw target bootstrap)
| arm | nmae (raw) | reach | fid | overall-equiv seed 0 / 1 / 2 [seed-0 CI] | P(≤ 0) per seed | sign flips |
|---|---|---|---|---|---|---|
| mixed k 0.1 | −0.0024 | +0.0073 / −0.0061 / −0.0061 | +0.0001 | +0.0021 / −0.0004 / −0.0004 | 0.00 / 0.60 / 0.60 | 2.2 % |
| **mixed k 0.2** | **−0.0030** (3/3) | **+0.0033 / −0.0029 / −0.0031** | +0.0001 | **+0.0015 / +0.0004 / +0.0003** [−0.0014, +0.0049] | 0.17 / 0.43 / 0.38 | 3.5 % |
| mixed k 0.3 | −0.0018 | +0.0033 / −0.0029 / −0.0071 | −0.0001 | +0.0010 / −0.0001 / −0.0009 | 0.29 / 0.53 / 0.70 | 4.4 % |
| magnitude-only k 0.1 | −0.0013 (3/3) | 0 | 0 | +0.0004 / +0.0004 / +0.0004 [−0.0001, +0.0009] | 0.08 / 0.06 / 0.06 | 0 |
| magnitude-only k 0.2 | −0.0008 | 0 | 0 | +0.0002 | 0.3–0.4 | 0 |

- **The nmae part reproduces 480** (mixed k 0.2: −0.0030 raw on every seed, against 480's −0.0030 from adding to the
  delivered lfc directly). Carrying it through the carriers' means, with the ballast as delivered, loses almost nothing.
- **The mixed form's sign flips (3.5 % of cells) move reach by ±0.003 raw, with a sign that depends on the emission
  seed** (seed 0 +0.0033, seeds 1–2 −0.003). The flipped genes sit near the top of the reach walk in some realisations
  and not others. That seed noise is as large as the nmae gain in overall units, so the per-seed intervals all span 0.
- **The magnitude-only form is clean but small:** nmae −0.0013, no flips, fid/reach/jac exactly 0, overall-equiv
  +0.0004 on every seed. On the board, if nmae prices as a pure magnitude lever (×1 … ×2.3), that is about
  +0.0004 … +0.0009 overall.

## Against the preregistered GO rule (mixed k 0.2) → NOT GO
| item | rule | measured | pass |
|---|---|---|---|
| 1 | overall-equiv > 0 on 3/3 seeds, bootstrap lower bound > 0 on ≥ 2 | +0.0015 / +0.0004 / +0.0003; lower bounds −0.0014 / −0.0034 / −0.0025 | **no** |
| 2 | nmae ≤ −0.0015 raw | −0.0030 | yes |
| 3 | non-essential overall-equiv > 0 | +0.0018 / +0.0004 / +0.0003 | yes |

**Reading.**
- Close the mixed rider at this construction: its sign component is unpriceable (calibration rules 3–4) and
  seed-fragile here.
- The magnitude-only k 0.1 form is the only clean piece, worth ≤ +0.001 on the board. It is a rider to include only if
  an emitter hook for arithmetic-mean multipliers gets built for another reason.
- Not priced here: class capacity. 426/435 found adds above the MDE's dose leak calls, and 0.1·|cons| exceeds the MDE
  cap on strong donor genes. Only an emitter run measures that.

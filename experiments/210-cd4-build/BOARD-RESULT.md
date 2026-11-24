# `cd4f034` board result — **+0.0130557**, rank 7, the largest single-arm gain in the project

Entry `QSe3JlNUYZ2BzvlcsL14`, model name `The Ribosome Has Filed A Grievance`, **empty
description**, submitted at the 2026-09-04 00:00 UTC slot reset (armed to fire itself).
**Overall `0.2049195830` against `rc_callset`'s `0.1918638367` — `+0.0130557`, rank 7.**

Preregistration: `PREREGISTRATION.md` + `PREREGISTRATION-ADDENDUM.md`, both committed before the
slot. One variable against the incumbent: `VCC_CD4_TARGET_F` **0.10 → 0.34**.

| member | `rc_callset` | `cd4f034` | Δ scaled | Δ overall | rc raw | cd raw | Δ raw |
|---|---:|---:|---:|---:|---:|---:|---:|
| `pds` | 0.745871 | 0.737800 | −0.008071 | −0.001345 | 0.83712 | 0.83347 | −0.00365 |
| **`mse`** | 0.120057 | **0.159276** | **+0.039220** | **+0.006537** | 0.87367 | **0.83635** | −0.03732 |
| `nmae` | 0.094290 | 0.104861 | +0.010571 | +0.001762 | 0.94326 | 0.93687 | −0.00639 |
| `fid` | −0.029268 | −0.020127 | +0.009141 | +0.001524 | 0.50370 | 0.50642 | +0.00272 |
| **`reach`** | 0.192898 | **0.219300** | **+0.026402** | **+0.004400** | 0.24983 | **0.27321** | +0.02338 |
| `jac` | 0.027335 | 0.028408 | +0.001072 | +0.000179 | 0.04043 | 0.04082 | +0.00039 |
| | | | **+0.078334** | **+0.013056** | | | |

**Predicted +0.00437 conservative / +0.00572 calibrated. Realised +0.0130557 — 3.0× / 2.3×.**

## The falsifiers

| # | test | result |
|---|---|---|
| 1 | raw `mse` must improve on 0.8738045, else the 096/197 identity does not transfer at dose and the CD4 dose family is retired | **PASS**, 0.83635 |
| 2 | raw `mse` in [0.836, 0.852]; the addendum revised the point prediction to **≈0.8333** on the measured `A` shortfall | **PASS at the very edge**, 0.83635 — 0.003 above the revised point, inside the original band |
| 3 | `pds` must not fall > 0.004 raw | **PASS**, −0.00365 (predicted −0.0011, so real but inside) |
| 4 | raw `nmae` must not rise > 0.002 | **PASS**, it **fell** 0.00639 |
| 5 | `fid`/`reach`/`jac` must not move > 0.004 scaled | **FIRES — all three, all UPWARD**: +0.0091 / +0.0264 / +0.0011 |

**Falsifier 5 fired in the favourable direction and that is the most informative thing here.**
`PREREGISTRATION-ADDENDUM.md` §3 classified the arm shape-only on the a14 → a15 board pair, whose
unpriced channels moved **+0.0005 overall favourably**. This arm's moved **+0.0061 overall
favourably — 12× that.** The classification was right in direction and badly understated in size;
the +0.0089 amplitude floor would have vetoed a +0.0131 arm.

## What it proves

**The transfer cosine rose ×1.249.** By 197's one-sided inversion, `cos ≥ (1.023611 − raw)/2.074495`
takes our floor from **0.07228 → 0.09027**. This is the **first arm in the project's history to move
the direction channel on the board**, and it moved five of six members at once — everything except
`pds`, which paid the −0.0037 raw the sweep predicted it would.

**And the board converts this lever ~1.8× better than PLOLO said.** 212 measured **+13.83%** PLOLO
cosine at `f = 0.34`; the board delivered **+24.9%**. That ratio is now the most important open
number for experiment 216 — see `experiments/216-cd4-hard/INBOX-2`.

## Position

**Rank 7. Overall 0.2049196.** Gap to `cqawesome` (0.22082) **+0.01590**, from +0.04383 at
2026-09-03 00:30 UTC. **Three slots, +0.0308331: a15 0.1740865 → 0.2049196 in two UTC days.**

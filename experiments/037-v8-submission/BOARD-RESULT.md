# v8 board result: the mechanism was exact, the magnitude was 8x over-predicted

**Submitted** 2026-08-22 13:48 UTC as `We regret to inform the genes`.
**Overall 0.13689 → 0.13771 (+0.00082).** Rank 2 of 122.

| metric (scaled) | v6 | v8 | delta |
|---|---|---|---|
| `pds` | 0.7013719706431124 | **0.7013719706431124** | **0** |
| `mse` | 0 | 0 | 0 |
| `nmae` | 0.06517348739698008 | 0.06517740507579672 | +0.0000039 |
| `reach` | 0.08987 | 0.09460 | **+0.00473** |
| `jac` | 0.001539 | 0.002692 | +0.00115 |
| `fid` | −0.036603 | −0.037569 | −0.00097 |

## What worked exactly

**`pds` is identical to all sixteen significant figures, and `nmae` to six decimals.** The
exact-total quasi emitter does what it was built to do: the pseudobulk is unchanged, so the
three metrics that read it cannot move. Every point of the board delta is attributable to the
emission change and nothing else. That is a cleaner controlled experiment than a leaderboard
usually permits.

## What did not

Predicted **+0.007**, delivered **+0.00082** — a factor of 8.5.

| | h1-fast predicted | board delivered | ratio |
|---|---|---|---|
| `reach` raw | 0.1669 → 0.1911 (+14.5%) | 0.1526 → 0.1568 (+2.8%) | **5.2x** |
| `fid` raw | +0.0045 | −0.0003 | wrong sign |
| `jac` raw | −0.004 | +0.0012 | wrong sign |

The direction on `reach` was right and it is the largest single component. `fid` and `jac` both
came back with the opposite sign to the prediction, though both are small enough to be noise on
a 40-perturbation harness.

## The diagnosis, and it invalidates a rule I wrote yesterday

I claimed `h1-fast` governs anything routed through `n_real` "because it has 400 reference cells,
as the competition does". That is half the story and the wrong half was load-bearing:

| harness | perturbations | cells/pert | **control cells** |
|---|---|---|---|
| h1-fast | 40 | 400 | **4,000** |
| somatic HCT116 | 132 | **150** | 18,400 |
| **h1-test** | **85** | **400** | **18,400** |
| **competition** | 300 | 400 | 18,400 |

`n_real` depends on the Wilcoxon's power, which depends on **both** groups. `h1-fast` matches on
the perturbed side and is **4.6x short on controls**; the somatic harness matches on controls and
is 2.7x short on the perturbed side. **Neither matches the competition's DE regime**, and the
emission lever acts precisely through that regime — which is why both over-predicted it, in the
same direction, by roughly the same factor.

**`h1-test` matches on both axes and has been sitting unused all along.** It is the harness that
should have priced this change. It costs ~100 minutes per configuration rather than 7, which is
exactly why it was skipped.

## What this is worth beyond the 0.0008

The clean freeze on `pds`/`mse`/`nmae` turns v8 into a **calibration measurement**: the
DE-metric transfer ratio from `h1-fast` to the board is now measured at roughly **0.2**, not
assumed. Any future emission or DE-calling change screened on `h1-fast` should be discounted by
about 5x before deciding whether it is worth a slot — or screened on `h1-test` instead.

The earlier log entry already suspected this shape ("the harness overstates `fid` gains by
~2.75x, because its baseline sits much lower"). The correction is the same one, larger, and now
attributable to a specific missing 14,400 control cells.

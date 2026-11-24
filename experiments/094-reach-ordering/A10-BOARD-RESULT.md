# a10 board result — the emission frontier is flat, and the falsifier fired

Entry `lwHUEPPfuUAwvrfwrkne`, submitted 2026-08-28 13:35 UTC. **Overall `0.1541841`** against
a8's `0.1543607` — **`-0.0001766`**. Rank unchanged at 16.

Config: a8 base with `--tail-spread 0.25 -> 0.15` and `--overdispersion 4.0 -> 16.0`.

| member | a8 | **a10** | Δ |
| --- | ---: | ---: | ---: |
| `pds` | 0.736560993522496 | 0.736560993522496 | **+0.000000** |
| `mse` | 0 | 0 | 0 |
| `nmae` | 0.094301 | 0.094338 | +0.000037 |
| `fid` | −0.028694 | −0.027897 | +0.000797 |
| `reach` | 0.124349 | 0.123595 | **−0.000754** |
| `jac` | −0.000353 | −0.001493 | **−0.001140** |
| **overall** | **0.1543607** | **0.1541841** | **−0.000177** |

## The falsifier

`docs/CURRENT_STATE.md` pre-registered: *"if a10 scores below a8's 0.1543607, then `reach`
magnitude beats the `n_pred` charge on the challenge, orion was right, and the next move is
higher spread, not lower."*

**a10 scored below a8. The falsifier fired.** But it fired by `-0.00018`, and the honest
reading is not "orion was right" — it is that **the difference is too small to have a sign we
should trust.** Moving `n_pred` from 3,836 to 4,757 — a 24% change in the call set, driven by
a 4x change in `phi` and a 40% change in `spread` — moved the overall score by less than two
parts in ten thousand.

## What this actually establishes

**The (spread, phi) emission surface is flat near our operating point.** Two harnesses
disagreed about which side was better (orion ranked a8 first, h1-test ranked a10 first); the
board says the question was not worth asking, because the surface has no gradient there.

The prediction that h1-test's larger `n_real` made it the more trustworthy ranking was
**not** vindicated: it predicted a10 > a8 and a10 came in marginally below. Nor was orion's
ranking vindicated in magnitude — it predicted a8 ahead by `+0.0097` scaled and the gap was
`+0.0002`. **Both harnesses over-stated a difference the board could barely measure.**

The one thing that reproduced exactly, for the third consecutive submission, is
`pds` Δ = `0.000000` under a large emission change. The exact-total repair holds.

## Consequence for how the remaining slots are spent

Nine experiments (088-095) have now been spent on the emission channel. It produced one real
gain — a8's `+0.0041` from over-dispersion — and has now returned two consecutive results
inside noise (`a2-prune` `+0.00013`, `a10` `-0.00018`).

**The emission channel is exhausted.** `spread`, `phi`, `n_pred`, the whitening basis and the
call-set geometry have all been tuned to the point where the board cannot resolve further
moves. Remaining headroom is in the two members we have never moved — `mse`, which is at a
hard 0, and the model quality underneath `reach` — not in how the counts are emitted.

**Do not spend another submission on an emission parameter without a mechanism that predicts
a move larger than `0.001`.**

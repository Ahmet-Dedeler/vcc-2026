# a13 board result — the mse lever works and the family is dead

Entry `RnN67G6OWerJ1dTMANzQ`, 2026-08-29. **Overall `0.1476855`** against a12's `0.16798741`
— **`−0.020302`**. The largest single-arm loss this project has recorded.

| member | a12 | **a13** | Δ | gate predicted |
| --- | ---: | ---: | ---: | ---: |
| **`mse`** | 0.091163 | **0.145924** | **+0.054762** | +0.0150 / +0.0214 |
| `pds` | 0.737885 | 0.709106 | **−0.028779** | +0.0044 / −0.0086 |
| **`nmae`** | 0.093974 | **−0.029332** | **−0.123305** | +0.0026 / −0.0145 |
| `fid` | −0.031473 | −0.036758 | −0.005285 | −0.0011 / +0.0057 |
| `reach` | 0.116790 | 0.097121 | −0.019668 | unreadable |
| `jac` | −0.000413 | 0.000051 | +0.000464 | +0.0006 / +0.0000 |
| **overall** | **0.167987** | **0.147686** | **−0.020302** | **+0.0215 / +0.0040** |

Raw: `mse` **0.901349 → 0.849307**. `nmae` **~0.9441 → 1.018502**.

## What worked, and it worked better than predicted

**The shrinkage is the strongest `mse` lever this project has found.** Raw `mse` fell
`0.9013 → 0.8493`, a `−0.052` move against a `−0.013` model — **4x the prediction**, and it puts
us within `0.022` of the leader's `0.8270` on that member. `score_mse` `+0.0548` is larger than
any single-member gain we have ever recorded, a11's `mse` included.

## What killed it

**`nmae` at raw `1.0185` is past that member's own baseline anchor (~1.0013) — worse than
predicting nothing at all**, exactly the failure the `gamma = 1.0` default was demoted for.
It happened at `gamma = 0.25`, the demoted-to value, four times gentler.

The mechanism was understood in advance and mis-sized by an order of magnitude: `nmae` scores
absolute fold-change error over **the reference's** DE genes, and we shrink by **donor
agreement**, which is uncorrelated with which genes the reference calls DE. So the shrinkage
lands squarely on genes that are scored for fold-change accuracy.

**Both harnesses under-called it by ~8x.** h1-test said `+0.0026`, orion said `−0.0145`; the
board said `−0.1233`. Rule 1's "take the worse surface" was applied and was still 8x optimistic.
`pds` was also mis-called — the gate had it at `+0.0044 / −0.0086` and the board delivered
`−0.0288`.

## The verdict, and it is not "tune gamma"

Even granting the full `mse` gain and zeroing the `nmae` cost entirely, the remaining members
still net out at `(+0.0548 − 0.0288 − 0.0053 − 0.0197 + 0.0005)/6 = +0.0003` — **inside noise.**
A gentler `gamma` scales the gain down with the costs. **The donor-agreement shrinkage family is
closed**: it buys `mse` and pays it back across `pds`, `nmae` and `reach` at roughly one-to-one.

That is the board-scale confirmation of what three DE gates had already suggested at harness
scale: within this family the `mse` objective is not merely an incomplete view of the total, it
is very nearly its negative.

## What it cost, and what it bought

It cost a slot and dropped the displayed entry from `0.16799` to `0.14769` — **the leaderboard
shows a team's latest entry, not its best**, which had been assumed the other way and is worth
knowing. Under latest-wins the risky arm of a pair should go first, so a slot remains to restore.

What it bought is worth more than the slot: **`mse` is now a solved, measured, movable member**
— we know exactly how to reach `0.849` and roughly what `0.827` would take — and the shrinkage
route to it is definitively closed. The remaining question is no longer *whether* the metric can
be moved but *which* member to move, and the gap decomposition now answers that.

# 217 CANCELLED before submission — experiment 216 measured the answer first, and it is a VETO

**2026-09-04 01:30 UTC. No slot was spent. The build was killed mid-emission.**

## What happened

216's dose sweep finished while 217's build was still running, and it priced `f = 0.80` **directly**
rather than by extrapolation, using **188's exact `nmae` identity** instead of the proxy 207 used:

| operating point | ovr `pds` | ovr `mse` | ovr `nmae` | **TOTAL** | vs `cd4f034` |
|---|---:|---:|---:|---:|---:|
| `f = 0.28` | −0.00009 | +0.00414 | +0.000078 | +0.00413 | −0.00024 |
| **`f = 0.34` ← `cd4f034`** | −0.00068 | +0.00521 | −0.000152 | **+0.00437** | — |
| `f = 0.40` | −0.00135 | +0.00612 | −0.000683 | +0.00408 | −0.00029 |
| `f = 0.55` | −0.00335 | +0.00779 | −0.002279 | +0.00216 | −0.00221 |
| **`f = 0.80`** | **−0.00688** | +0.00912 | **−0.004298** | **−0.00205** | **−0.00642** |
| `f = 1.00` | −0.00960 | +0.00938 | −0.005256 | −0.00547 | −0.00984 |

**`cd4f034` is the argmax of the exact six-member price over the entire 31-point grid and all three
forms.** Its own neighbours are −0.00024 and −0.00029. 216 reproduces 207's three published
selection values exactly, so this is the same instrument with the proxy replaced.

## Why my prediction was wrong, and it is a real lesson

The addendum predicted **+0.004 to +0.017** by scaling `cd4f034`'s realised member deltas: positives
by the cosine ratio, `pds` cost by ×2–×4. **The `pds` half was roughly right** (216: −0.00688 at
`f = 0.80` against my −0.003 to −0.007). **The error is that I scaled all five positive members
together, and they do not move together.**

**`nmae` has a turning point and goes NEGATIVE with dose** — +0.000191 at `f = 0.20`, −0.000152 at
0.34, **−0.004298 at 0.80.** And §8.4's own `nmae` identity says exactly why: `nmae` is optimal at
`k = 1` with slope **−(2q−1) ≈ −0.074 below it and +1 above it**. Past the optimum, raising the dose
pushes DE-pool magnitudes over their true values and the +1 slope takes over. **I had the mechanism
written down in the same document and still scaled `nmae` as if it were a positive member.**

**Corrected rule, for §8.4 and for every future direction arm:** a direction gain propagates to five
members **near the optimum**, but `nmae` turns over first and `pds` turns over at `f ≈ 0.20`
(219). **Do not scale a member past a turning point you have already measured.**

## The decision-tree branch this is

`DECISION-TREE.md` named branch **D** (overall falls → revert to `cd4f034`, do not spend a slot to
restore rank, re-point 216 at the *form*). That branch fired **before the experiment ran**, from a
better local measurement rather than from a board result — which is the cheaper way for it to fire
and exactly what running the sweep in parallel with the build was for.

216 also closed the form: the library and 1/k routes price at **−0.051 to −0.137**, and the
published atlas is **byte-for-byte the object we already ship**, so the *Cell* date coincidence is
confirmed as coincidence and the dose cannot be what the leaders did (the lever's absolute ceiling
is **+32.4%** cosine against their measured **+106%**).

## Slot status
**The second 2026-09-04 slot is UNSPENT.** Per `AGENTS.md` §6, a slot must buy information toward
first place; a −0.0064 arm buys neither score nor information we do not already have. It is held for
218's candidate or for whatever 220/221 surface.

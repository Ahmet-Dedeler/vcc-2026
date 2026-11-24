# Pre-decided: what happens after `cd4f080`'s board result

Written 2026-09-04 00:50 UTC, **before the build finished and before the result exists**, so the
next move is not chosen after seeing the number. `PREREGISTRATION.md` §6 already names the branch;
this fills in the rest of the tree.

Incumbent `cd4f034` = **0.2049196**, rank 7. Leader `cqawesome` = **0.22082**.
Predicted `cd4f080` = **0.219 – 0.222**.

## The branches

**A — lands at or above ~0.215, `pds` cost within falsifier 3 (< 0.012 raw).**
The dose-response is still climbing at 207's PLOLO maximiser and the board converts it ~1.8×.
→ **Next slot: `f = 1.4`.** Not 1.2 — three board points (0.10, 0.34, 0.80) will pin the curvature,
and if it is still near-linear in `ln f` the information is in going further, not in halving the
step. 216's extended sweep should be in by then; **if 216 has located a PLOLO turnover below 1.4,
take 216's number instead.** Re-check `pds` against the realised cost at 0.80, not against 207's.

**B — lands positive but under ~0.212, i.e. the gain per e-fold is decaying.**
The curve is bending. → **Next slot: `f = 0.55`**, to locate the interior optimum between 0.34 and
0.80 rather than pushing past it. Two points bracket it; the third finds it.

**C — falsifier 3 fires (`pds` falls > 0.012 raw) but overall is still up.**
`pds` is the binding constraint, exactly as 207's guard suspected and the realised −0.00365 at
`f = 0.34` (3.3× its prediction) hinted. → **Next slot: `f = 0.55`**, and hand 216 the realised
`pds`-vs-dose curve so its sweep is re-guarded on board numbers instead of PLOLO ones.

**D — overall falls.**
The dose-response turned over between 0.34 and 0.80. → **Revert to `cd4f034` as the shipped
configuration**, do **not** spend a slot restoring it (rank is informational; `AGENTS.md` §6), and
re-point 216 at the *form* — the library form, the operator re-fit with CD4 in its basis, and the
published-atlas question — since the dose axis would then be exhausted.

**E — falsifier 1 fires (raw `mse` does not improve at all).**
The 096/197 identity does not transfer at dose. → **Stop the ladder entirely.** That would mean the
`f = 0.34` result was not a cosine gain, which would also invalidate §8.4's multiplier and force a
re-read of `cd4f034` itself. Highest-information failure on the tree; treat it as a finding, not a
setback.

## What does NOT change on any branch

- `VCC_ALPHA` stays **0.982** (207 §6.2 measured α\* = 0.982; 211 re-solved 0.9808/0.9879, all inside
  noise). Never re-solve it in the same arm as a dose change — that is two variables.
- The emitter stays `rc_callset` and the call set stays 192's union. Both are board-proven and
  orthogonal to this axis.
- **Never re-submit an already-scored package to restore rank** (`AGENTS.md` §6). If a dose point
  loses, the displayed row drops and we leave it.
- Every subsequent dose arm gets its own preregistration with all six members priced under §8.4's
  rule — **never `mse` alone**, which is the convention that cost 4.5×.

## The standing risk, stated once
Three board points on one parameter is a thin basis for a ladder, and the whole thing rests on a
single measured conversion factor (board ≈ 1.80 × PLOLO on cosine). **If branch C or D fires, the
factor is not general and every extrapolation in `AGENTS.md` §8.4 and in 216's mandate needs
re-deriving from the new points.** That is the honest failure mode and it is cheap: it costs one
slot and returns a corrected constant that every future direction arm uses.

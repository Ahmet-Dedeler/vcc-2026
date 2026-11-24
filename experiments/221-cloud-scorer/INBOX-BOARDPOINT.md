# INBOX — a SECOND board point landed mid-run, and it contradicts 216's transfer factor

`h1union` published at **0.2118310, rank 3** (entry `iUziqE1gewfcPKzYUDaB`), against `cd4f034`'s
**0.2049196**. **Δ = +0.0069114, all six members up.** Full per-member table in
`experiments/222-h1-union-build/BOARD-RESULT.md`.

**This doubles your calibration data and it hands you a contradiction that is worth more than the
score.** On the matched energy-weighted-aggregate functional (`cos_agg`, the one 197's inversion
actually returns):

| arm | PLOLO rel `cos_agg` | board rel cosine | factor |
|---|---:|---:|---:|
| `cd4f034` — 216 §4.4 | +17.98% | +24.89% | **1.384** |
| `h1union` — 222, 3 cancer folds | **+10.3065%** | **+10.3045%** | **1.000** |

**Same family, same instrument, two factors.** One of these readings is wrong, or the factor is
conditional on something neither experiment controlled. Hypotheses to separate, none tested:

1. **Fold set.** 222 excluded `kolf2` from every headline (it is pluripotent, like the added H1
   source, and gains 3–5× more). Did 216's ladder aggregate over four folds or three? If four,
   recompute 216's row on three and see whether 1.384 moves toward 1.0.
2. **Functional/tag convention.** 216 explicitly diagnosed this trap — *"212's macro per-target
   cosine vs 197's energy-weighted aggregate"* — and it cost it four Phase-0 gates. **Check whether
   it re-entered its own trap when computing the 1.384 row.** Gate (`min_cpm`), tag, and target-
   coordinate exclusion must match between the two sides.
3. **Genuinely conditional.** The factor may depend on where the added source sits relative to the
   stack already held — a first source into empty rows versus a fourth estimate onto occupied ones.
   If so it is not a constant and the calibration must be expressed as a function.

**Your deliverable was already the six-member local→board transfer table. This makes it the
highest-value artefact in the project**, because every future arm is sized through it and we now
know the current number is unreliable by ~38%. **Add both board points as known-answer rows and do
not report a single scalar factor unless the two reconcile** — a factor that is conditional, stated
as a function with its condition named, is worth far more than a wrong constant.

Nothing else in your mandate changes. Keep the reaper results and the housekeeping proof.

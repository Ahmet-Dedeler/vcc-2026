# ADDENDUM to the `cd4f080` preregistration — written 2026-09-04 01:15 UTC, BEFORE the build finished and BEFORE submission

The preregistration is **not edited**. This records information that arrived from experiment 219
*after* it was committed and *before* the measurement exists, together with the revised prediction.
Acting on a better prior before an experiment is legitimate; revising an analysis after seeing its
result is not, and this is the former.

## 1. The multiplier I built the prediction on was double-counted

219 factored §8.4's 2.51× exactly as **1.997 (propagation) × 1.256 (transfer residual)** and showed
that 207's predicted `mse` had **already** passed the PLOLO cosine through `k_cos = 1.2598` — so the
1.256 *is* Error 2, and multiplying by 1.80 charged it twice. **The honest end-to-end multiplier is
2.51 (range 2.00–2.56), not 4.5.** `AGENTS.md` §8.4 is corrected.

**This does not weaken the case for the arm** — the preregistration's prediction was built from
`cd4f034`'s *realised* member deltas scaled by the PLOLO cosine ratio, not from the multiplier — but
it does mean the framing "the instrument understates by 4.5×" was wrong and is withdrawn.

## 2. The propagation factor FALLS with dose, and that widens the range a lot

219's decisive addition: `pds`-per-unit-cosine is **+0.2797 at `f = 0.10`** and **−0.0748 at
`f = 0.34`** — **the `pds` leg turns over at `f ≈ 0.20`.** So `pds` cost cannot be scaled linearly to
`f = 0.80`, which is what the preregistration did (it assumed ×2 to ×4).

Re-deriving from `cd4f034`'s realised deltas, with the cosine gain ×2.28 and the `pds` cost scaled
several ways:

| `pds` cost scaling | `pds` at `f = 0.80` | net vs `cd4f034` | predicted board |
|---|---:|---:|---:|
| linear ×2.28 | −0.00306 | **+0.0167** | 0.2216 |
| quadratic ×5.2 | −0.00698 | **+0.0128** | 0.2178 |
| cubic ×11.9 | −0.01589 | **+0.0039** | 0.2088 |
| quartic ×27 | −0.03620 | **−0.0165** | 0.1885 |

**Revised prediction: +0.004 to +0.017 over `cd4f034`, board 0.209 – 0.222**, against the
preregistration's narrower +0.0143 to +0.0170. **The arm still clears its gate under three of four
scalings and is roughly break-even under the fourth.** It stays.

## 3. And the bar it is judged against moved in our favour
219: **§8.2's "+31% of cosine for a single direction arm" was itself an `mse`-only bar and is
2.0–2.3× too high.** Corrected: **+11.6–13.6%** against the amplitude floor. `f = 0.80` projects to
**+57%** board cosine — comfortably over the corrected bar rather than marginally over the old one.

## 4. What I am NOT changing, and why
The obvious defensive move is to drop to `f = 0.55`. I am not taking it:
- three of four `pds` scalings clear the gate, and the worst is ≈ break-even, so the downside is
  bounded and the upside reaches the leader;
- the branch tree in `DECISION-TREE.md` **already names `f = 0.55` as the response** if falsifier 3
  fires, so that option is preserved rather than spent;
- a widely-spaced third point (0.10, 0.34, 0.80) pins the **curvature**, which is exactly what
  §2's turnover question needs and what a timid step would not deliver.

**Falsifier 3 (`pds` falls > 0.012 raw) is now the most likely to fire**, and 219 has given it a
mechanism rather than a hunch. If it does, that is the turnover located — a result, not a loss.

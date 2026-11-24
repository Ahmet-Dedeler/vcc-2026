# INBOX 4 — 228 just landed a constraint your six-member price must be read against. `nmae` is the gate.

228 finished: **224's per-gene reliability gain is a NO-GO**, and it died in a way that bears
directly on `ARM6`.

**The gain survived the scorer read-back on the channels you are reporting** — **+4.55%** relative
aggregate cosine, **Δ`pds` +0.00484**, positive **3/3** cancer folds. **And its six-member composite
is −0.00702.** `nmae` alone contributed **−0.01043**.

**The mechanism.** 188's `nmae` splits as matched `M` (= sign-matched `min(|pred|,|truth|)`) over
unmatched `U`. **`M` cannot collect a magnitude rebalance; `U` collects all of it.** 224's gain moved
|t|-weighted **sign accuracy by +0.00027 — nothing** — so 93% of the magnitude it pushed onto the
reference's DE pool was unmatched: efficiency `e = ΔM/ΔU = 0.072` against break-even **0.5**.

**So "3/3 cosine, 3/3 `pds`" is NOT sufficient to clear.** It was true of 228's arm too.

## What I need from your `RESULTS.md`

1. **The full six-member composite for `ARM6`, with `nmae` explicit**, on the three cancer folds
   against the shipped `h1union`. If you have it, lead with it. **A cosine/`pds` headline without
   `nmae` is not a price** — that is exactly how 228's arm looked one step before it died.
2. **State whether `ARM6` moves SIGNS or MAGNITUDES.** I expect **signs**, and the precedent is
   strong: `h1union` was also a donor-estimate addition and **its `nmae` went UP +0.00588** with all
   six members rising. A source adds *information about direction*; a per-gene rescale does not.
   **But measure it — report |t|-weighted sign accuracy with and without the arm**, the same statistic
   228 used. If `ARM6`'s sign accuracy moves and its `nmae` holds, that is the positive control that
   distinguishes it from 228's arm, and it is worth one run.
3. **If `nmae` falls materially, say so and mark the arm NOT recommended** regardless of cosine and
   `pds`. We have two slots tomorrow and no obligation to fill both. An unspent slot beats a
   composite we misread.

Nothing else changes. The registry work and its known-answer hashes stand on their own merits
whatever the arm does, and **INBOX-1 still asks you to commit the builder** once your gates pass.

# INBOX 2 — 230 just landed a number that bears directly on your verdict. Read it before you write one.

230 finished (`experiments/230-toll-magnitude/RESULTS.md`, 27 known-answer rows, all pass). Two of
its results change how your arm must be read.

## 1. The toll is MIS-KEYED — charge 0. This helps you.

Rule 13's `−0.0146` `reach` toll is keyed on the wrong test. Its own calibration points (a21, a22)
are **digest-identical to the incumbent on 900/900 blocks**, so by the table's own test they belong
in the **zero-toll** row; the row the toll is charged to **contains zero board points**. `AGENTS.md`
rule 13 is corrected: **charge no `reach` toll.** If you were going to report your arm both ways,
you no longer need to — **report it at toll 0**, and say you are following 230.

## 2. 🔴 The noise floor on YOUR channel — this is the hard part

> **The unexplained spread of board Δ`reach` among DIGEST-IDENTICAL arms is
> sd 0.027011 scaled = 0.004502 on overall — 2.4× the entire slot floor.**

**219 says every point of Arm A is in `reach`.** So your headline range of **+0.0046 … +0.0088** is
**1.0 to 2.0 standard deviations of the noise on the very channel it lives in.** At the low end the
arm is indistinguishable from a re-roll of the same submission.

**And note the asymmetry, because it is the whole point:** `pds` and `mse` are *deterministic* on
this board — `rc_hybrid` reproduced both **to sixteen significant figures**. **The noise is specific
to the ordering channel.** An arm that delivers cosine/`pds` is measurable at a precision an
ordering arm simply is not.

## What I need in your verdict

1. **State your price in units of that noise**, not only against the floor: "+X, i.e. Y sd of the
   0.004502 board `reach` noise." That single number is what decides the slot.
2. **If the arm's gain is entirely in `reach` and lands under ~2 sd, say it does not clear** — even
   if it beats the floor arithmetically. A floor is not a noise model, and 230 has now measured the
   noise. **This is the same call I just made against 208's Hamilton (1.16× floor, margin ≈1/15th of
   the noise): live candidate, not cleared, no slot.** Apply it consistently to Arm A.
3. **If any part of the gain is in cosine or `pds`, price and report that part separately** — it is
   measurable where the `reach` part is not, and it may be the shippable half of the arm even if the
   whole does not clear.
4. **Phase 0's reconciliation against 220's bound still governs.** If the implied raw `reach`
   exceeds the field best, the arm is refuted on arithmetic regardless of everything above.

**"Real but unmeasurable at this precision" is a complete and honest verdict**, and it is more
useful than a number that beats a floor we now know is narrower than the noise.

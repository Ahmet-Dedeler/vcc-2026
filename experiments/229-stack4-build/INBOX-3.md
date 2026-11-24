# INBOX 3 — on your "the pre-registered arm is a measured regression". Do NOT build a package you already know is bad.

You report that **a single-union build cannot express 225's construction** and that the literally
pre-registered arm regresses against the incumbent. **That is a real and valuable finding** — it is a
gap between the instrument that priced the arm and the path that ships it, and finding it is exactly
what a build step is for. Take the credit for it in `RESULTS.md`.

My mandate told you to *"build the pre-registered one anyway"* if you disagreed with `g`. **That
instruction was about the DOSE, and it does not apply here.** It exists to stop an agent tuning a
dose against the instrument that scores it. It was never meant to make you ship a construction you
have measured as a regression. **Do not build a package you already know is worse than what is on
the board.**

## What I want instead

1. **Characterise the gap precisely.** Say what 225's construction does that a single union with one
   `g` cannot — per-source doses? sequential residuals? separate carrying-row normalisation per
   source? — and **which of those differences carries the regression.** That is the durable result.
2. **Then find the best arm the builder CAN express**, and pre-register it **before** you build it,
   in a file, with a one-line reason.
3. **Choose that arm on 225's measurements, not on your own build's score.** 225's harness is
   independent of the builder and reproduces 223 and 218 exactly; your build is the thing being
   tested. **Selecting a construction by what your own build scores is fitting the arm to the
   instrument** — that is the failure the original instruction was guarding against, and it applies
   with more force here, not less.
4. **`dld1_lowmoi` alone is the obvious fallback** — 223 priced it at **+0.00228**, 225 reproduced
   that to the last digit through an independent harness, it is a **single** source so a union of
   one expresses it exactly, and it is one of the few sources whose October coverage **improves**
   (67 → 73 on the pessimistic draw, §8.3). **A clean +0.00228 that we can actually build beats a
   +0.00421 we cannot.**
5. **If NOTHING expressible clears the gate, say so and build nothing.** An unspent slot is worth
   more than a regression, and we have two slots tomorrow. **`NOT BUILDABLE` with a clear reason is
   a complete and successful outcome for this experiment.**

## Still required either way

The **registry change itself** is valuable independent of which arm ships — it is the thing standing
between us and shipping any new source, and `docs/NEXT-BUILD-SPEC.md` §3's two known-answer hashes
still govern it. **Land the registry with its gates passing and commit it** (INBOX-1), plus the four
inert October sources (INBOX-2), **even if the answer on the arm is `NOT BUILDABLE`.**

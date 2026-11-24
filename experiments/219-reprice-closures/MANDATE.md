# 219 — Re-price EVERY closure under §8.4's corrected rule. We under-priced direction arms 4.5×.

You are an audit subagent, in the tradition of 181 and 198. Work ONLY in
`experiments/219-reprice-closures/`. Read-only elsewhere. No `vcc`, no builds, no `cell-eval2`.
Peak RSS ≤ 10 GB. `uv run python`.

## READ FIRST
`AGENTS.md` **§8.4** (the correction), §8.1, §8.2, §11.1 rules 1–13 **with the rule-13 table**.
Then `experiments/181-audit-0902/RESULTS.md` and `experiments/198-closure-reaudit/RESULTS.md` —
**the two prior audits. Their method is your method; do not repeat their findings.**

## The correction, and why it invalidates a class of verdicts

`cd4f034` predicted **+0.00437** and delivered **+0.0130557**. The post-mortem (§8.4) locates two
systematic conventions, both present in essentially every direction price this project has computed:

1. **A direction gain propagates to FIVE members, not one.** The four we carried at ~zero —
   `nmae`, `fid`, `reach`, `jac` — delivered **60%** of the realised gain, because 187 proved `nmae`
   moves only through `q_w`, `reach` is `f(sign accuracy, ordering)`, and `fid` is our call set's
   precision. **All are sign-accuracy members and cosine sets sign accuracy.** `mse`-only pricing
   understates by **2.51×**.
2. **PLOLO understates the board's cosine response by 1.80×** (212 measured +13.83% at `f = 0.34`;
   the board delivered +24.9%).

**Combined ≈ 4.5×.** Any NO-GO whose margin was inside that factor was decided by a convention, not
by the data.

## Work

### Phase 0 — calibrate the correction, do not inherit it (mandatory, abort on failure)
- Re-derive the 2.51× and 1.80× from the artifacts, not from §8.4's prose:
  `cd4f034`'s board row against `rc_callset`'s, 212's PLOLO forward measurement, and 197's
  inversion `cos ≥ (1.023611 − raw)/2.074495` giving 0.07228 → 0.09027.
- Reproduce **one** number from each audit you are extending, so your instrument is comparable to
  theirs: 181's and 198's own known-answer rows are listed in their RESULTS.
- Reproduce PLOLO's a15 `pds` **0.795784** / a23 **0.789121** / mass-fix **+0.004014**.
- **State the correction's uncertainty.** It rests on **one** board point. Give the multiplier a
  range, and carry that range through every re-price rather than a point estimate. An audit that
  turns one measurement into a fleet of confident reversals is worse than no audit.

### Phase 1 — enumerate and classify every closed arm
Sweep every `experiments/*/RESULTS.md` with a terminal verdict. For each, record: the verdict, the
**margin** against its gate, **which members it priced**, **which it carried at zero**, and **which
instrument** it used. Then classify:

- **Class A — re-priceable**: a direction/profile arm priced on `mse` and/or PLOLO with DE members
  carried at zero, whose margin is inside the corrected multiplier.
- **Class B — unaffected**: negative by more than the multiplier, or priced on all six already, or
  closed by a **ratio** argument (212's gamma closure is the canonical case — the operator trades
  cosine for retrieval at 96.2% of the board's own rate, and a ratio is instrument-independent).
- **Class C — closed structurally**: closed by a fact, not a number (205: the artifacts do not
  exist; 194: 92% of targets carry one guide so guide-averaging has no axis; 184: the winning law is
  sign-preserving and absent in CRISPRi).

**Report the full table.** The classification is the deliverable, more than any individual reversal.

### Phase 2 — re-price Class A, with the DE members ON
For each Class A arm, re-price using §8.4's rule and **say explicitly what it becomes**. Candidates I
expect to matter, and you should not limit yourself to them:
- **193** (operator eigenvalue map): best composed arm **+0.000512** against +0.0019.
- **189 / 190** (ordering keys): +0.00408 ± 0.00424 and +0.00137 — **and note that 214 measured a
  77× tightening of the `reach` instrument under `rc_callset` (three-seed sem 0.00009 against 200's
  0.00695), which is a separate reason these two are re-openable.** Coordinate: experiment 220 owns
  the ordering keys specifically — hand them over rather than duplicating.
- **208 / 214**'s riders, whose floors turned on the rule-13 class table.
- Anything in 200's ten-component ablation, **noting 214's finding that 200 and 208 both priced on
  a renderer that has not shipped since a23** — Hamilton is 1.84× larger on the shipped renderer,
  PerturbMatch 1.77× smaller. **Every row of that ablation inherits the mis-scaling.**

### Phase 3 — rank what to actually re-run
Output a ranked shortlist: arm, corrected six-member estimate with its range, what would have to be
re-measured to confirm it, and the cost of doing so. **Rank by (corrected estimate × confidence) ÷
cost**, and be explicit that a re-price is a *hypothesis*, not a result — nothing goes near a slot
without a fresh measurement and a preregistration.

## Falsifiers
1. **The correction must not reverse a verdict whose own evidence contradicts it.** For each proposed
   reversal, check whether the original experiment measured the DE members directly (some did) — if
   it did, the 2.51× does not apply and you must say so.
2. **A reversal that depends on the multiplier's upper end is not a reversal.** Flag those separately.
3. **Do not re-open anything closed by a permutation control that tied or beat the real arm**
   (204's donor weighting, 209's subspace deflation, 190's shuffled regulon). Those are dead
   regardless of scale.

## Deliverable
`RESULTS.md`: the Phase-0 re-derivation with the multiplier's range, the full classification table,
the Class-A re-prices, and the ranked shortlist. **"Nothing reverses" is a full result** — it would
mean the eleven direction closures are real and the remaining time belongs to the dose ladder and
to acquisition.

## Territory
216 owns the CD4 dose/form. 217 is building `f = 0.80`. 218 owns donor independence and acquisition.
220 owns the ordering keys. **You own the audit.** Do not re-measure anything yourself beyond what
Phase 0 needs — your job is to say what is worth re-measuring.

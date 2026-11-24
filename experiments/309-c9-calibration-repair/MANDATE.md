# 309 — Claim C9: repair the two calibration defects that aborted C8, and re-test the SAME hypothesis at the SAME α. The bias repair works AGAINST the effect, so a surviving result is stronger, not weaker.

You are a research subagent. Work in `experiments/309-c9-calibration-repair/`. **Edit nothing outside it** except
appending C9's preregistration and result to `docs/CONFIRMATION-LEDGER.md` under a heading you create,
`## Claim C9 — 309 (calibration repair of C8; no reserve)`. No `vcc`, no `cell-eval2`, no build, nothing under
`submissions/` or `data/`, no download. `uv run python`. Peak RSS under ~4 GB. **Time box: 3 hours.**

## Why, and the honest framing you must preserve
307 measured the aggregate candidate-repression statistic with the construct as the unit and got
**S = +0.786 [+0.387, +1.166]** (n = 34), joint-null p at the floor 1.0e-5, surviving a real-spacer null
(p = 0.0015), a wrong-recipient permutation (p = 2.0e-4) and a partner-construct specificity test (the same
candidate genes read at wrong-recipient level, +0.041 [−0.295, +0.369], in the cells of the partner construct that
knocks down the SAME intended gene with DIFFERENT protospacers). **But 2 of 7 preregistered known-answer rows
failed and 307 correctly aborted to `UNRESOLVABLE-AT-THIS-RESOLUTION` rather than loosen its own rule:**
- **KA-4** (no-signal row): a construct's candidate genes read slightly repressed in *other* constructs' cells —
  mean **+0.0683** against a `≤ 0.050` rule (post-hoc at 5,000 derangements: +0.0717 ± 0.0022). Real, and **9 % of
  the observed +0.786**.
- **KA-6** (matching balance): one axis of nine missed — promoter-window count SMD **+0.2905** against `≤ 0.250`.

**These are defects of the null/matching machinery, not of the hypothesis.** Repairing a demonstrably
miscalibrated null and re-testing is legitimate; silently relaxing a threshold is not, and 305's coordinator
decision (`experiments/305-t7-expression-test/COORDINATOR-DECISION-2026-09-07.md`) is binding precedent that you
must not do the latter. **Three rules bind you:**
1. **Same hypothesis, same α (0.05, one-sided), same ONE verdict cell, same statistic family as C8.** You may not
   redefine the effect to make it easier to detect.
2. **C8's bounded estimate stands and is reported alongside C9 whatever happens.** C9 does not replace or retract
   it. Report explicitly that this is a **second test of the same hypothesis on the same data** and say what that
   costs in interpretation (it is not an independent replication; state the multiplicity plainly).
3. **If a known-answer row fails again, abort again.** The abort clause is not weakened because it fired once.

**State this in the preregistration before running:** the KA-4 repair removes a positive bias and therefore
**shifts the estimate DOWN by roughly 9 %**; it makes passing harder. The KA-6 repair may move it either way. You
are not permitted to choose between repair variants after seeing their effect on the primary.

## The two repairs (specify exactly, in `PREREGISTRATION.md`, before computing the primary)
- **KA-4 / residual bias.** Either (a) make the statistic bias-immune by centring each construct's score on its own
  derangement distribution (score minus the mean of that construct's wrong-recipient scores), or (b) subtract the
  estimated construct-independent bias with a propagated interval. **Choose ONE and justify it from the
  diagnostic, not from its effect on the answer.** Whichever you choose, the new KA-4 must be recomputed on the
  repaired statistic and must meet the original `≤ 0.050` rule.
- **KA-6 / matching.** Add **promoter-window count** to the matched axes (it is the axis that missed) and re-run
  the joint matching. Report the full nine-axis SMD table for the repaired pools; every axis must meet `≤ 0.250`.
  If adding the axis makes the matched pools too small for some construct, say so and report the resulting
  attainable resolution rather than dropping constructs silently.

## Carry forward from 307 (reuse its code and artifacts; do not re-derive)
`experiments/307-offtarget-enrichment/` — `PREREGISTRATION.md`, `RESULTS.md`, and `work/`: `lib307.py`,
`seq307.py`, `pool307.py` (**the 209,070-spacer rescan, already corrected for the Horlbeck `="SYMBOL"`
Excel-escaping bug that had inflated the per-spacer L=12 hit rate 0.107 → 0.597 — inherit the CORRECTED pool and
carry a known-answer row proving you did**), `enrich307.py`, `enrichment307.csv`, `ka307.py`, `results307.json`.
Also 305's saved integer sums (`sums305_<split>.npz`) and 301's rule. **Re-run 307's seven known-answer rows
unchanged** plus the two repaired ones; report all nine.

## Keep, unchanged, as declared secondaries (no verdict)
The real-spacer resampling null; the wrong-recipient permutation; the partner-construct specificity test; the
failing-vs-working contrast (301's original T7 question); the 31 non-targeting constructs run identically and
labelled as candidate-bearing observations, not negative controls.

## Deliverable
`PREREGISTRATION.md` (hashed into `work/PREREGISTRATION.sha256`, appended to the ledger as C9 **before** any
repaired primary value is computed), `RESULTS.md`, `work/enrichment309.csv`, and one line:
**C9 (aggregate candidate-repression enrichment, construct as unit, calibration-repaired): PASS / FAIL /
UNRESOLVABLE-AT-THIS-RESOLUTION (S … [lo, hi]; p …; nine KA rows …/9; C8's estimate was S = +0.786
[+0.387, +1.166]).** State in the abstract what a positive does and does not establish: it would NOT establish
dCas9-KRAB off-target binding, because expression here is observational and the partner-construct test narrows
the explanation to something carried by the protospacer sequence without separating seed-mediated binding from a
construct artefact correlated with spacer identity. No model, no submission, no production change. Then stop.

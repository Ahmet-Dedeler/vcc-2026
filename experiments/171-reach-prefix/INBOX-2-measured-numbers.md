# INBOX 2 for experiment 171 — numbers already measured tonight, 2026-09-02 06:55 UTC+3

Two sibling experiments have produced numbers that bear directly on your mandate. Use them; do not
re-derive them from scratch, but **do** reproduce the known-answer row before trusting your own
pipeline.

---

## 1. h1-test is CONFIRMED to reproduce our board DE members — your instrument is sound

Experiment 169 is running the six-member table on `h1-test` right now and its known-answer row for
the shipped a15 arm reads:

| arm | n_pred | reach f1 | precision | fid raw | jac raw | recall | nmae | pds | **sign acc** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `KA_a15_raw` | 6034 | **0.18283** | 0.50701 | 0.50701 | 0.06824 | 0.66852 | 0.79944 | 0.53435 | **0.53435** |

Compare to our **live board** raw values: `reach` **0.17422**, `fid` **~0.505**. h1-test reads
0.18283 and 0.50701. **The surface reproduces both to within about 0.009 raw**, which is inside the
0.014 AGENTS.md already records for reach. So for this experiment specifically, h1-test is a
trustworthy instrument and rule 9's "unreadable" caveat is weaker than usual — **but still report
deltas as your primary quantity.**

## 2. THE DEFECT, MEASURED: our call set has sign accuracy 0.534

That last column is the finding. **Of the genes we call differentially expressed, 53.4% have the
correct sign.** Chance is 0.50. Precision against the reference's DE set is 0.507.

This is the numerical form of the pattern in INBOX 1: our profile ranks 5th in the world on `pds`
while the calls we derive from it are **barely distinguishable from coin flips**. A call set at
53.4% sign accuracy cannot sustain a 90%-pure prefix for any depth, which is exactly why `k*`
collapses and `reach` reads 192nd.

**This reframes your task sharply.** The question is not only "which genes go in the significant
block and in what order" — it is **"why is the sign of our own called genes almost uninformative
when the profile those calls come from is world-class?"** Candidates worth separating:

- the calls are dominated by genes where the predicted effect is ~0 (157's 74–79% zero-effect
  finding) so their sign is numerical noise;
- the emission's cell-to-cell variance decides significance nearly independently of the predicted
  effect, so the call set is close to a random draw over expressed genes;
- the sign is being read off a quantity that the renderer or the mass fix has altered.

Distinguish these — a `sign accuracy vs |predicted effect| decile` table over the call set would
separate the first two immediately, and it is cheap.

## 3. A live contradiction you should be aware of but NOT resolve

Experiment 161 (just closed NO-GO on its own question) reports in passing that a zero-learning
per-entry re-shaping of the donor consensus, `C ← sign(C)·|C|²`, is worth **+0.0118 worst-line /
+0.0206 mean raw pds** across five lines in the delivered space.

That is **magnitude sharpening at q = 2**, which experiment 168d closed as a **TERMINAL NO-GO**
because its gain did not survive on h1-test (GLOBAL q=1.15 read +0.00112 on h1-test against
+0.01145 on Orion, a 10.2× overstatement). 161's five lines are K562/RPE1/Jurkat/HCT116/HEK293T —
**none of them h1-test, and two of them are the discredited Orion pair.**

The most likely reading is that 161 re-measured 168d's Orion overstatement on three more lines
without knowing it, which would make it *corroboration of the instrument rule* rather than a
contradiction. **Experiment 169 owns this question — do not spend time on it.** It is here only so
that if your own work touches magnitude, you know the family is closed on the surface you are
using and you should not re-derive a gain from it.

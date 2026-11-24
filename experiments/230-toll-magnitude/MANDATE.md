# 230 — Is rule 13's `reach` toll a STEP or does it SCALE with the profile change? It gates a live +0.0022 arm and a whole class.

You are a research subagent. Work in `experiments/230-toll-magnitude/`. **No `vcc`, no submission
build, no production change. No `cell-eval2` while another agent runs one — check first.** Peak RSS
≤ 10 GB. `uv run python`. Four other agents are live (221, 226, 228, 229); write only inside your
own directory.

## The question, and why it is unresolved rather than unasked

`AGENTS.md` §11.1's rule 13 charges a **single-valued** `reach` toll, keyed on the per-gene
panel-total digest:

| class | toll |
|---|---:|
| emission-law **profile-CHANGING** | **−0.0146 scaled** |
| emission-law profile-PRESERVING | 0 |
| profile-class | 0 |

**That table treats the toll as a STEP: any profile change at all costs the full −0.0146.** Nobody
has tested whether it instead **scales with the size of the change** — and the two arms it was
calibrated on both moved the profile *enormously* (a21 took `n_pred` **6,034 → 1,913** and `nmae`
**−0.064 scaled**).

## What it currently blocks — this is not hypothetical

219 §2.2 re-priced **208's Hamilton rider** and found that 208 and 200 had both measured on
`VCC_MASS_FIX=minnorm`, a renderer **only a23 ever shipped** (a23 scored −0.00073 and was
abandoned). On the renderer that actually ships, 214 measured the factor at **1.84×**:

| | value | vs the +0.0019 floor |
|---|---:|---|
| 208's original rule-13 composite | +0.001197 | VETO |
| **corrected to the shipped renderer** | **+0.002203** | **CLEARS, 1.16×** |
| six-member, corrected | **+0.003909** | |
| …with the **−0.00243** emission-law toll charged | **−0.000227** | fails |

**So the arm lives or dies on the toll alone.** 208 proved it is not pseudobulk-preserving (900/900
blocks differ), so **the letter of the table charges it** — but Hamilton moves the challenge tilt
from **−0.01253 to +0.00030 log2**, i.e. it *removes* a systematic rather than adding one, and that
is a **0.03%** profile change against calibration points that were catastrophic.

219 flagged this and was **explicitly forbidden** from testing it (*"any re-measurement of a Class-A
arm — out of mandate"*). **You are the follow-up.**

## Work

### Phase 0 — known-answer rows (mandatory)
Reproduce: the toll's two calibration points (**a21, a22**) and their board deltas; 208's original
**+0.001197**; 214's **1.84×** renderer factor; and `rc_hybrid`'s digest result (it is the arm that
refuted the naive profile-PRESERVING row). If these do not reproduce, stop.

### Phase 1 — fit the toll against profile-change MAGNITUDE
**We have many board points and the toll has only ever been read as a class constant.** For every
scored arm with a known emission law — at least **a21, a22, `rc_hybrid`, `rc_callset`, `cd4f034`,
`h1union`**, plus any others in `docs/BOARD-LEDGER.md` — compute:

1. a **magnitude** of profile change (the per-gene panel-total digest gives the test; you need a
   *size*: total-variation or relative-L1 on the per-gene panel totals, and the challenge tilt in
   log2, are the two obvious candidates — report both);
2. the arm's **realised board `reach` delta**, isolated from its other channels as far as the ledger
   allows.

**Then fit `reach` delta against magnitude, and test STEP vs LINEAR vs SATURATING.** State the
evidence for each and which the data actually supports. **A step function and a steep linear are
distinguishable precisely because Hamilton sits at 0.03% while the calibration points sit at
catastrophic** — the leverage is real.

### Phase 2 — price Hamilton under the supported law
If the toll scales, price 208's rider under the fitted law and report the six-member composite with
a confidence interval. **If it is a step, say so and Hamilton is dead** — that is a complete result
and it retires a candidate cleanly.

### Falsifiers
1. **Out-of-sample.** Fit the law leaving out one calibration arm; predict it. A law that cannot
   retrodict a held-out board point is not a law.
2. **Confounding.** a21/a22 changed the **call set** as well as the profile (`n_pred` 6,034 →
   1,913). **Separate the call-set effect from the profile-change effect, or state plainly that you
   cannot** — if the −0.0146 is really a call-set toll misattributed to the profile, that is a
   bigger finding than Hamilton and it would re-open every arm rule 13 has ever charged.
3. **Do not fit and price on the same points** without saying so.

### Reporting
**Do not convert PLOLO to board with a scalar** — the factor is **not identified** (1.384 on
`cd4f034`, 1.000 on `h1union`, matched `cos_agg`). **And only invert a board row whose `score_mse`
is strictly positive** — 15 of the 55 teams above 0.15 are clamped at exactly 0 and 197's inversion
returns cosines for them that were never measured. Both traps have already cost this project real
errors today.

## Deliverable
`RESULTS.md`: the reproduced known-answer rows, the fitted law with its evidence and out-of-sample
test, the call-set/profile separation, Hamilton's corrected price, and a one-line verdict. **Either
answer is valuable**: a step retires Hamilton and confirms rule 13 as written; a scaling law
re-opens Hamilton **and every other arm rule 13 has charged**, which is much the larger result.

**Propose no slot. Build nothing.** The coordinator decides.

## Territory
221 calibration · 226 the operator's fitting rows · 228 the per-gene gain's PLOLO price ·
229 the source registry and build. **You own rule 13's toll and nothing else.**

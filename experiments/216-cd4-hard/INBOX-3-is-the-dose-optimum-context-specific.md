# INBOX 3 for 216 — report the dose curve PER CONTEXT. It is the one D/E/F risk nobody owns.

Coordinator (lead), 2026-09-04 01:10 UTC. Small addition to your sweep, large consequence for
October. Costs you nothing extra if you stratify what you are already computing.

## The risk

`cd4f034` earned **+0.0130557** on the board and the dose ladder is now the main line. **But every
number we have is an A/B/C aggregate, and the only entry that counts is on D/E/F**, whose lineages
we will not know until Oct 22 and whose panel is a different 300 genes.

207 established the gain is **lineage-generic** — +7.6% to +15.2% on four unrelated non-T lines —
which is why `AGENTS.md` §8.3's decision table ships CD4 *whatever* D/E/F's lineages turn out to be.
**That result is about whether the gain exists. It is not about where the optimum sits.**

A lineage-generic gain is entirely compatible with a **context-specific optimal dose**: context A is
**JURKAT** (041, unanimous across all eight feature/rank variants), a T-lymphoblastic line, and the
donor is a primary **CD4⁺ T-cell** atlas. If A's optimum sits above B's (HELA, cervix) and C's
(CAL33, head-and-neck squamous), then the aggregate optimum we are about to ship is a **weighted
average of three different optima**, and on three *new* lineages it could sit anywhere in that
spread — or outside it.

This is exactly the shape of failure that **`abc-knee` already produced once**: a per-context
quantity fitted on A/B/C, shipped as a global, and falsified on the board (Δ`mse` −0.00226 against a
predicted +0.005). We should not walk into the mirror image of it in the dark window.

## What to add — stratification, not new work

You are already sweeping `f` on PLOLO's four folds in the challenge frames. **Report every curve
per context as well as pooled:**

1. **The per-context cosine-vs-`f` curve**, with each context's own maximiser. Report the **spread**
   between the three maximisers explicitly — that number is the deliverable here.
2. **The per-context `pds` cost curve**, since 217's realised `pds` cost at `f = 0.34` came in
   **3.3× larger** than 207 predicted and `pds` is now the binding constraint.
3. **Whether A's maximiser is the highest.** If it is, the lineage story survives as a *dose* effect
   even though 207 refuted it as an *existence* effect, and D/E/F should ship a **lower** dose than
   A/B/C's aggregate. If A's is *not* the highest, lineage is irrelevant to the dose too and the
   aggregate transfers cleanly — which is the better outcome and worth stating plainly.
4. **The loss from shipping the pooled optimum to a context whose own optimum differs.** Evaluate
   each context's curve at the other two contexts' maximisers. If the curve is flat near its top,
   the whole risk evaporates and we ship the aggregate with confidence; if it is sharp, October
   needs a rule.

## Why this decides an October rule, not just a number

`docs/DEF-RUNBOOK.md` §3's decision table currently pre-decides `VCC_CD4_TARGET_F` as a **single
global value**, with the note that a board-validated dose supersedes `0.1`. **If your per-context
spread is wide and the curve is sharp, that row needs a conditional** — computable from the D/E/F
controls alone, the way 187 wrote the α-knee rule with a numeric trigger rather than a judgement.
If the spread is narrow, the row stays a constant and we have retired a risk for the cost of a
`groupby`.

**Either answer is worth having, and only you are positioned to produce it** — 217 is building a
single global dose and 218 is searching for new sources. Do not re-derive anything; just stratify.

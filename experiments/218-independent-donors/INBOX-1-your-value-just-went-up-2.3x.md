# INBOX 1 for 218 — the bar you are judged against was wrong. Your target just re-priced 2.3× upward.

Coordinator (lead), 2026-09-04 01:15 UTC. From `experiments/219-reprice-closures/RESULTS.md`.
Verify before relying on it.

## What changed

`AGENTS.md` §8.2 carried a bar of **"+31% of cosine for a single direction arm to clear the gate"**.
219 showed **that bar was itself computed `mse`-only** and is therefore **2.0–2.3× too high**.
Corrected: **+11.6–13.6%** against the amplitude floor, **+2.5–2.9%** shape-only. (§8.4 is also
corrected — my 4.5× was a double-count; the honest multiplier is **2.51**, range 2.00–2.56.)

**The direct consequence is about your mandate.** 205 priced a fifth panel-covering donor line at
**+0.0062 overall** and concluded it sat **below** the single-arm bar. Under the corrected bar the
same line prices at **+0.0140 to +0.0163** — from **1.4× below** to **1.6–1.8× above**.

**166's feasibility verdict is unchanged** — no fifth genome-wide CRISPRi line covering the panel
exists. **What changed is the value of finding one**, and it is now 2.3× larger. That is exactly the
asymmetry your mandate was built on: 166 searched under the criterion *"genome-wide human Perturb-seq
covering the panel"*, which **would have missed the CD4 atlas**, and CD4 is our single best
component — it delivered **+0.0131 on the board** yesterday.

## What to do with it

- **Do not lower your standards to find something.** A source that is not genuinely independent is
  worth nothing; 204 and 209 both died on permutation controls that tied or beat the real arm.
- **Do raise how far you are willing to look.** A candidate that gets you to **+0.014** justifies
  substantially more acquisition and processing effort than one worth +0.006. Primary cells,
  iPSC-derived and differentiated systems, organoids, CRISPRa, base editing, non-10x chemistries —
  all in scope, and a partial-coverage source is now worth building if its independence is high.
- **Your Phase 0 matters more, not less.** The independence-vs-contribution regression on the ten
  sources we already hold is still the gate: if independence does not predict contribution, none of
  this pays and the honest answer is to say so.

## One number to carry
The corrected shape-only bar is **+2.5–2.9% of cosine**. A source that adds even a few percent of
*independent* cosine now clears a gate it would have failed a day ago. **Report every candidate's
projected cosine contribution against the corrected bar, not the old one**, and state which bar you
are using in every table — the old number is still written in several experiment RESULTS files.

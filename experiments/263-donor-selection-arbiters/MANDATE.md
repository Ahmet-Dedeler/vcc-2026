# 263 — H2 of 262: choose the donor per target where the donors disagree, with a BUILD-TIME-OBSERVABLE arbiter. Screen every arbiter we can observe against the oracle ceiling. Measurement only, ≤ 3 hours. DO NOT BUILD, DO NOT SUBMIT.

You are a research subagent. Work in `experiments/263-donor-selection-arbiters/`. **Edit nothing outside it.** No build,
no `vcc`, no `cell-eval2`, nothing under `submissions/`. `uv run python`. One PLOLO fold per process, peak RSS under
~8 GB (262's pilot staged three folds in one process at 8.85 GB — do not). **Time box: write up after ~3 hours.**

## Why (R262 §5–§8)
`cd4svd5` loses `pds` on a fold-consistent set of fully-supported targets whose whitened donor CONSENSUS points away
from the held line's truth: 17.5 % of the panel scores below chance (53 % of the deficit) and 66 % sits in the marginal
band (base·truth cosine 0–0.05) where **+0.01 cosine buys +0.095 `pds`** (R262 §4.2). The equal-weight mean of two
disagreeing donor rows is a compromise pointing at neither: on the marginal band the BEST single donor's agreement with
the held line is **1.7–2× the mean's**, and the oracle choice would move 20–50 % of the band past the 0.05 saturation
point (R262 §7.2). But the one arbiter 262 tried — the CD4 residual's agreement with each donor — **picks the wrong
donor**, and no observable bin finds the below-chance group (§7.3). **Question: is there ANY build-time-observable
arbiter that captures a material fraction of the oracle headroom, with the controls below half the arm?** Production
weights the three donors 1 : 1 : 1 (RB decision table: equal weighting was best under union whitening — that was a
CONTEXT-level choice; this is PER-TARGET).

## READ FIRST
`docs/READ-FIRST.md` (v3.1; §3 pricing rule and the marginal-band reading rule), `experiments/262-where-pds-is-lost/RESULTS.md`
§4–§8 and its `diag262.py` / `pilot262.py` (REUSE: the per-target arrays, the covariates, the `cd4svd5` construction from
261's `cd4svd5_parts`, the fold caches — byte-copy `work/plolo/`), `experiments/162-retrodicting-validator/plolo.py`
(which donors each fold whitens; `kolf2` is the ONLY fold with all three production donors — use it for any
three-donor arbiter and report it apart), `experiments/246-direction-hypotheses/lib246.py`, `AGENTS.md` §2.1, §4.2
(donor weights), §8.1b, §11.1 rule 10.

## Construction
On each fold, the whitened base row for target t is the equal-weight mean of the fold's available donor rows
(two on the cancer folds, three on `kolf2`). Replace it, per target, by a selection or a weighting of those rows, then
add the residuals exactly as `cd4svd5` does (CD4 deflated k = 5, re-solved at f 0.34; the H1-pair union at g 0.28).
Everything downstream unchanged. Score `pds`, `cos_agg`, macro, per-target `pds` (the scorer's definition, as 262).

**Known-answer / headroom rows first:** (R1) equal weights reproduce `cd4svd5` at 0.0e+00 (262 KA-1); (R2) the ORACLE
selection (per target, the donor row with the higher cosine to the held truth — never an arm, the CEILING) and its
`pds`/cosine; (R3) 262's CD4-as-arbiter row reproduces §7.2's numbers.

## Arbiters to screen (each with its gating: apply only where donor-donor cosine < τ, τ ∈ {0, 0.1, always})
1. CD4 residual row → donor with the higher cosine (262's; the control that must reproduce and must NOT work).
2. H1-union row → donor with the higher cosine (only where the union carries the target; else mean).
3. **Third donor (on `kolf2` only, three donors available): pick the two that agree, drop the odd one out** — the
   production-realisable "majority" rule; and the pairwise version on the cancer folds via the held-out donor's
   *other* fold appearances if the caches allow (say so if not).
4. Donor response strength: the donor with the larger whitened row norm / the larger DE count in its own line
   (from the donor tables' per-target statistics if on disk; else norm).
5. Donor self-consistency: for each donor, the cosine of its row for t to its rows for t's nearest-neighbour targets
   (kNN in the OTHER donor's space) — a donor whose response is coherent with its neighbours is trusted.
6. Continuous weighting instead of selection: w ∝ softmax(β · arbiter score), β ∈ {1, 3}, for the two best arbiters.
Controls for every arbiter: random per-target choice (seeded, 3 seeds), the REVERSE of the arbiter, and the mean.
Report per fold: Δ`pds` vs `cd4svd5`, 3/3 sign count, LOO σ, Δ`cos_agg`, the fraction of the oracle headroom captured,
how many marginal-band targets crossed 0.02 and how many fell below 0, leave-top-3-out retention.

## Price (READ-FIRST §3: coherent pessimistic / point / optimistic; the class factor for a base re-direction is
UNMEASURED — use the re-shape band as pessimistic, k_pds 0.85 / k_cos 3.76, and the re-dose band as point, 1.22 / 1.36;
`nmae` by the identity ×1.5 on a loss; amplitude leg if the delivered ratio moves > 0.2 %). CLEARED = pessimistic
≥ +0.0019 over `cd4svd5` AND retention ≥ 50 % AND `pds` > 0 on 3/3 folds at ≥ 2 LOO-σ AND every control below half.

## Deliverable
`RESULTS.md`: R1–R3; the arbiter × gating table; the controls; the oracle-headroom fractions; the best arbiter in
full; the builder variable it would need (none exists — name it, e.g. `VCC_DONOR_SELECT=<arbiter>:<tau>`, and the
one-place change at the base-merge step, do not make it); one line **CLEARS / DOES NOT CLEAR**; if it clears,
`PREREGISTRATION-DRAFT.md` in 249's form. "No observable arbiter captures more than a tenth of the oracle headroom"
is a valid result that closes H2. No model name, no description, no submission. Then stop.

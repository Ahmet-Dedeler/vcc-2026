# 267 — Reference reliability: how much of the donor-selection oracle's +0.046 `pds` is REPRODUCIBLE biology? The replicate-validated oracle advantage. Feasibility stage, ≤ 4 hours. Measurement only. No build, no submission.

You are a research subagent. Work in `experiments/267-replicate-validated-oracle/`. **Edit nothing outside it.** No
build, no `vcc`, no `cell-eval2`, nothing under `submissions/` or `data/harness/`, nothing written under `data/`. `uv run
python`. Raw-cell extractions are the heavy step: one line at a time, streaming/backed reads, peak RSS under ~10 GB;
three sibling subagents (268–270) run alongside at ≤ 6 GB each. **Time box: write up after ~4 hours.** Programme:
`docs/RESEARCH-PROGRAMME-2026-09.md` (read it first; its rules bind you: name the experimental unit, refit inside the
split, strong baselines, paired bootstraps).

## Why
263 measured that choosing, per target, the whitened donor row with the higher cosine to the held truth lifts `pds` by
**+0.046** on the three cancer folds (`experiments/263-donor-selection-arbiters/RESULTS.md` §0.1), 3.7× the deflation,
65 % of it on 262's below-chance targets. But that oracle chooses WITH the held reference and is scored AGAINST the same
reference: part of the advantage may be selection against the reference's own noise, not transferable biology. 159
already reports large uncertainty in target-specific responses at 150 cells (`experiments/159-transfer-ceiling/RESULTS.md`
§1.2–1.3). **Question: choose the donor on one independent half of the held line's reference and score on the other —
how much of +0.046 survives, overall, per fold, per 262 group, per reliability tercile?** Outcomes: survives → real
response structure worth learning to select; disappears → response ESTIMATION must improve before any selector;
survives for particular target classes → a concrete biological hypothesis.

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; 263's `RESULTS.md` §1–§3 and `lib263.py` (the selection path, norm matching, the
oracle); 262's `RESULTS.md` §4–§5 (groups, covariates; its per-target arrays); `experiments/162-retrodicting-validator/plolo.py`
(how each fold's TRUTH is built: which raw data, which controls, pseudobulk → log effect; `fold_cache_status`);
`scripts/build_orion_harness.py` and `scripts/build_donor_table.py` (raw cell layouts for Orion HCT116 / HEK293T under
`data/external/orion/` and K562-GWPS under `data/external/scperturb/`, guide and batch columns; 260/266's
`work/extract_and_de.py` / `extract266.py` for the K562 streaming pass, 1.7 min); 159 §1 (how sampling noise enters the
residual space); `AGENTS.md` §2.1, §8.1b.

## Stage 1 (this mandate)
1. **Inventory (write it down first):** for each held line (K562-GWPS, Orion HCT116, Orion HEK293T) the per-target cell
   counts, the number of GUIDES per target, batch / gemgroup structure, control cell counts and their batches. Decide the
   split unit per line in this order and say which applies to how many targets: (a) by guide (≥ 2 guides: guide sets
   A/B), (b) by batch, (c) random cell halves. Controls are split the same way, so each half has its OWN controls.
2. **Half-truths in the instrument's format:** build truth_A and truth_B per target exactly as `plolo.py` builds the
   fold truth (same pseudobulk, same log effect vs the half's own controls, same gene axis and residual centring). **KA:**
   the truth from the UNION of halves (with the union of controls) reproduces the fold cache's truth to the instrument's
   precision (state the max |Δ| and why any residual difference exists, e.g. cell-count rounding).
3. **Reliability of the reference itself:** per target, cosine(truth_A, truth_B) in the residual space; its distribution
   per line; its relation to 262's groups (anti-aligned / marginal / saturated) and to cell count. This number is a
   scientific result on its own: how reproducible is a 150-cell Perturb-seq response on this panel?
4. **The replicate-validated oracle:** per target, choose the donor on half A (cosine of each whitened donor row to
   truth_A, exactly 263's `oracle` rule, norm-matched), then score `pds` on half B — with half-B truths of ALL targets
   as the distractor set, so the rank statistic is the scorer's. Compare: (i) the equal-weight mean (`cd4svd5`'s base)
   scored on half B; (ii) the in-sample oracle (chosen on B, scored on B); (iii) the cross-validated oracle (chosen on A,
   scored on B) and its mirror (B→A); (iv) random choice (3 seeds); (v) 263's best observable arbiter for reference.
   Report per fold and pooled: Δ`pds` vs the mean, paired target bootstrap (≥ 10,000), the fraction of (ii) that (iii)
   retains, and the same by 262 group and by reliability tercile of item 3.
5. **Sanity on the split**: scoring the equal-weight mean on half A vs half B must agree to within the half-truth noise;
   report it. If guide-based halves exist, also report whether the cross-validated advantage differs between guide
   splits and random-cell splits (a guide-specific effect would show there).

## Deliverable
`RESULTS.md`: the inventory; the KA; the reliability distribution (with a figure-free table); the oracle table with
bootstraps; the three-outcome reading (survives / disappears / class-specific) with the numbers that decide it; the
named reserve — **the `orion_hek293t` half-B truths are the programme's reserved confirmation data: save them under
`work/reserve/` with a sha256 and DO NOT use them for any selection in this experiment beyond item 4's single scoring;
log that scoring in `docs/CONFIRMATION-LEDGER.md` is the coordinator's job — just report the file and sha**. One line:
**SURVIVES / DISAPPEARS / CLASS-SPECIFIC (state the retained fraction)**, and what stage 2 would be. No model name,
no description, no submission. Then stop.

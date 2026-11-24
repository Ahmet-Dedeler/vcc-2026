# MANDATE — experiment 164: the call-set artefact and `fid` precision

## CONTINUATION NOTICE — read this block first

An earlier agent started this experiment and was stopped mid-run for runtime/cost reasons. **Its
state is in `experiments/164-callset-artefact/HANDOFF.md` — read that file in full before
anything else.** It records the scripts written (`blocks.py` ran cleanly; `realize164.py` was
launched in the background), every number measured, two findings that *reframe the thesis*, the
planned next steps, and five traps. Do not redo its work; continue from it.

Two facts the HANDOFF could not know:
1. **`realize164.py` DIED.** Its background process (PID 94461) is gone, `realize164.log` is 0
   bytes, and no `de164.*.npz` outputs exist. It was killed when the launching shell closed, or
   crashed silently. Rerun it in the foreground (or with `nohup` and stderr captured) and confirm
   the known-answer gates the HANDOFF lists: `a15` must reproduce 157's median `n_pred` 6035 /
   `reach` 0.1828; `armC` (a8 flat) must reproduce `reach` 0.2029. If they do not reproduce,
   stop and debug before measuring anything new.
2. The HANDOFF's reframing — **`fid`'s direction comes from the lfc (mean per-cell CPM), not
   the Wilcoxon z, so spread/phi arms change WHICH genes are called and never a gene's scored
   DIRECTION** — supersedes the original thesis statement below wherever they conflict. Its
   corollary (the sub-gate 5-10 CPM band's scored direction is a pure normalisation artefact,
   `lfc_0(g) ≈ log2(pb_cpm_g / cell_cpm_g)`, computable from the controls alone) is the single
   most important unmeasured quantity. Measure it first (`zero_lfc.py` as the HANDOFF planned).

Everything below is the original brief, kept verbatim so the goal and constraints are unchanged.

---

You are working in ~/Code/virtual-cell-challenge-2026 on the Arc Virtual Cell Challenge 2026. Read, in this order: `AGENTS.md` (§4 metric suite, §6 submission mechanics, §11.1 rules 1-12, the §6.4 harness warnings), `experiments/157-reach-order/RESULTS.md` §5 in full, and `experiments/156-delivered-geometry/{RESULTS.md,A23-BOARD-RESULT.md}`.

## HARD CONSTRAINTS
- **Do NOT submit.** No `vcc submit`. You MAY NOT create or upload a challenge package.
- **Do NOT run `cell-eval2`** (25 GB RSS, one at a time, machine shared).
- Write only inside `experiments/164-callset-artefact/`. Do not modify existing files elsewhere. If you need a builder change, write the patch as a `.diff` inside your directory and describe it — do not apply it.
- `uv run python`, never bare `python`. Peak RSS under 10 GB.
- Other agents own: validator calibration (162), the cross-line residual predictor (161). You own the CALL SET and the emitted p-value law. Do not touch donor weighting, whitening, the effect matrix, or amplitude.

## The thesis you are testing
Because `n_pred` (~4,400) far exceeds the challenge's `n_real` (~240, settled by the a21 board result), `de_wilcoxon_direction_fidelity_yield_raw` reduces to **precision on our own call set**: `fid = k / max(n_pred, n_real) = k / n_pred`. Our board raw `fid` is **0.5022** — chance. The board maximum is 0.5175 (`The_Way`) and the leader cqawesome is *below* us at 0.4877. **Nobody in this competition calls a DE gene's direction better than a coin flip.**

Experiment 157 §5b measured, in closed form on the real challenge controls (`data/challenge/context_{A,B,C}.h5ad`, so NOT harness-limited and NOT subject to AGENTS.md rule 9), that the shipped a15 emission law calls **28.8-31.4%** of gated genes significant at **zero predicted effect**. The cause is the tier spread policy: a15 ships `--tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 --overdispersion 4.0`, and on the ~1,900 tail-block genes the combination `spread 1.0 / phi 4` (gamma shape 0.25) sends most of a low-expression gene's 400 cells to zero, so the rank-sum reports a large, highly significant **DOWN** shift at zero true effect: **73.9-78.9%** of that block called significant, **100% DOWN**, median `|z| = 7.4` against the keep block's 1.7. *(HANDOFF correction: that "100% DOWN" is the z sign; the SCORED direction is the lfc sign, and 41-59% of that block is the sub-gate 5-10 CPM band where the predicted effect is exactly zero.)*

**If ~30% of our call set carries no model signal, it dilutes precision toward 0.5.** Removing it should raise `fid`; those same genes sit at the very front of `reach`'s confidence order (79.7% of the top ten on HCT116), so removing them should raise `reach`; and they are false positives, so `jac` should rise too. Those are three members where we rank 179th-207th of 578.

## Your task
**Decompose our precision, then price the arms that remove the artefact — on all six members.**

1. **Decompose `fid`.** Split our call set into (a) genes called because of the zero-effect emission artefact and (b) genes called because of real predicted effect, and measure the directional precision of each separately. The artefact set is identifiable at build time — 157's `zero_effect.py` computes it in closed form from the control cells alone. Reconcile your decomposition with the board's raw `fid` 0.5022: if artefact genes are ~30% of calls at ~0% skill, what precision must the signal genes have, and does that reconcile? Report the arithmetic explicitly. This is the load-bearing calculation; if it does not reconcile, say so and stop.
2. **Price the arms.** `src/vcc2026/donor_transfer.py:emit_counts(..., overdispersion=phi)` accepts **a scalar or one value per gene**. So a per-gene `phi` is already supported and needs no new machinery. Candidate arms, all one-variable:
   - **A**: tail block `phi = 0`, keep block `phi = 4` (a per-gene phi vector). 157 says this takes the mixture's zero-effect call rate from 28.8-31.4% to 14.6-20.7%.
   - **B**: `--tail-spread 0.5` instead of 1.0 (a15 currently 1.0).
   - **C**: revert to a8's flat `spread 0.25, phi 4` — **this one has a known board answer**, so it is your known-answer row: a8 scored raw `reach` 0.1890 and raw `fid` 0.5021 against a15's 0.1811 / 0.5022, and a15 beat a8 by +0.0197 overall via `mse` and `nmae`. Any instrument that cannot reproduce the sign of the a8→a15 `reach` and `mse` moves is not trustworthy for this family.
   - **D**: leave the emission alone and simply drop the artefact genes from the call set post hoc — NOT expressible in the submission format (we emit counts, not calls); treat it as a bound, not an arm, and say so.
3. **The binding constraint, which you must respect.** The tail block was introduced (a12) to raise emitted dispersion toward the metric's `#247` correction cap and thereby feed `mse`'s jackknife credit — `experiments/096-mse-floor/MSE-MECHANISM.md` is the reference. Removing it may cost `mse`. **Price `mse` explicitly for every arm.** An arm priced on `fid`/`jac`/`reach` alone is exactly the a19 failure mode (-0.01352) and will be rejected. The HANDOFF's suspicion — a15 sits ON the ρ cliff, so any `claim` reduction raises raw `mse` by ≈ Δclaim/0.229 with no amplitude compensation available — must be confirmed or refuted numerically, reusing `experiments/096-mse-floor/mech_claim_grid2.py` as the HANDOFF describes.
4. Use `h1-test` as the primary surface — 157 §1 showed it reproduces two board `reach` points to within 0.014 raw while repaired HCT116 misses a15's by 2.6x. **Beware: `data/donors/h1_2025` covers 85/85 of h1-test's targets — pure oracle leakage. Use `k562_gwps` + `orion_hek293t` only.** `debench.Surface(n_real_target=N)` uniform power-matching is the WRONG instrument for `reach`; scale each perturbation's pool proportionally instead.

## Deliverable
`experiments/164-callset-artefact/RESULTS.md`. Then a final message containing:
1. the precision decomposition and whether it reconciles with board `fid` 0.5022;
2. a six-member table for arms A/B/C with the known-answer row C reproducing a8→a15's signs;
3. a single recommended arm with the **exact** build invocation (env vars + CLI flags, diffed against a15's: `VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.10`, `--min-cpm 10.0 --overdispersion 4.0 --tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 --whitener-id union_pruned --whitener data/donors/whitener_union_pruned.npz --whitener-metadata data/donors/whitener_union_pruned
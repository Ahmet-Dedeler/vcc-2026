# 276 — The conservation map: for every panel target, how conserved is its response DIRECTION across cell lines and datasets on disk, replicate-validated — and does conservation predict transfer? Stage 2 of the organising question. ≤ 6 hours, feasibility first. Measurement only. No build, no submission.

You are a research subagent. Work in `experiments/276-conservation-map/`. **Edit nothing outside it.** No build, no
`vcc`, no `cell-eval2`, nothing under `submissions/` or `data/` (extractions and caches under `work/`; no download —
everything is on disk). `uv run python`. One raw-cell extraction at a time, peak RSS under ~10 GB. **Time box: write up
after ~6 hours; a partial map with its retrodictions is a valid result.** `docs/RESEARCH-PROGRAMME-2026-09.md` binds you.

## Why (FINDINGS-2026-09 F1–F4; R267, R273, R270 §4)
The programme's organising question — which perturbation responses are conserved across cell states, which depend on
the recipient, and when is a transfer unreliable — has so far been asked one dataset at a time. The reliability law
(F1) says a single 150-cell response barely resolves a direction; the Flex-400 regime (273) resolves them; the
selectors failed because no on-disk signal knows the held line (F2). **Nobody has yet built the per-target
CONSERVATION statistic across ALL the lines and datasets on disk and asked whether a target's conservation — a
property of the target, measurable before any held line is seen — predicts how well its response transfers.** If it
does, that is both a biological result (conserved vs state-dependent perturbations on this panel) and, unlike a per-
target selector, a prior that could inform how a residual is dosed per target class; if it does not, the map is still
the scientific deliverable.

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; `docs/FINDINGS-2026-09.md` (F1, F2, F4 and §3's table); 270's `RESULTS.md` §4 (the
on-disk multi-context table: CMap 38 lines × 148 panel targets across two reagents — `cmap270.py`, `feature270_cmap.csv`;
Jiang 6 lines × 218 targets; Nadig HepG2/Jurkat; Replogle K562/RPE1; Feng iPSC; HSPC; Tahoe is chemical, note it only)
and 176's platform-agreement result (the CMap known answer); 267's `inv267.py` / `truth267.py` and 271's `halfrows271.py`
(split-half reliability per target on our three lines — REUSE their half-banks for the replicate validation);
262's covariates (`cov262-*.npz`: base·truth cosine, cross-line agreement, donor-donor agreement, CD4 DE count) and
273's per-target `pds` on h1-test (`work/…pertarget…`); `scripts/build_donor_table.py` (how a line's responses are
turned into effect rows — use the same object everywhere); AGENTS.md §2.1 (the panel's design: non-essential).

## Stage 1 (this mandate)
1. **Inventory (write it first):** for each of the 300 panel targets, every (dataset, line) that measured its genetic
   perturbation on disk, with the unit (cells / guides / replicates / bulk), depth, and perturbation type; the leak
   surface of each source against the challenge lines and against our donors (which are themselves sources here —
   say so). Table: targets × number of independent lines / datasets.
2. **The conservation statistic (freeze it before computing):** for target t, the mean pairwise cosine of its response
   direction across all (dataset, line) pairs, in a common residual space (each source centred on its own across-target
   mean; genes on the common axis), with (a) a raw version and (b) a **reliability-corrected** version that divides by
   √(r_i r_j) where the split-half reliabilities exist (267/271/273's halves; for sources without cells, state the
   reliability as unknown and report (a) only); the number of pairs; a paired bootstrap over pairs; the null from
   target-shuffled pairings. Report the distribution over targets and the top / bottom deciles by name.
3. **Replicate validation:** conservation computed from half A of every source that has halves, and again from half B —
   the split-half reliability OF THE CONSERVATION SCORE itself (a target-level Spearman); a score that does not
   reproduce across halves is noise and must be labelled so.
4. **Does conservation predict transfer?** Against (i) `cd4svd5`'s per-target `pds` and base·truth cosine on the three
   PLOLO folds (262's arrays) and (ii) the incumbent's per-target `pds` on h1-test (273): Spearman with bootstrap; AUC for
   262's anti-aligned group; and the same after removing the held line from the conservation score (LEAVE-THE-HELD-LINE-
   OUT — the honest version: conservation among the OTHER sources only). If the LOLO version predicts transfer, that
   is the biological result. Baselines: 262's donor-donor agreement (the two-donor version of the same idea), CD4 DE
   count, expression level of the target.
5. **What conservation is made of (small):** for the top and bottom deciles, the DepMap co-essentiality cluster,
   expression breadth across the lines, and whether the conserved targets are the ones the whitener's top decile
   captures (the shared programme) — descriptive, no model.

## Deliverable
`RESULTS.md`: the inventory; the statistic's definition (frozen) and its distribution with intervals; the split-half
reliability of the score; the prediction table (raw, reliability-corrected, LOLO) with baselines; the decile
descriptions; one line: **CONSERVATION REPRODUCES ACROSS HALVES: YES / NO (ρ …); PREDICTS TRANSFER UNDER LOLO: YES /
NO (ρ …, AUC …)**; and the stage-2 recommendation (a target-class prior for the residual dose — priced later, not
here — or the closure at this construction). Save the per-target table as `work/conservation276.csv` (it is a
deliverable of the scientific track). No price, no model name, no submission. Then stop.

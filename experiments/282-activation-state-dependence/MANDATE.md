# 282 — State dependence within one individual: how conserved is a perturbation response across activation states (Rest → Stim8hr → Stim48hr) of the SAME donors' CD4 T cells, replicate-validated by donor? The organising question's cleanest test on disk. Preregistered analysis of 280's extracted bank, ≤ 4 hours. No build, no submission.

You are a research subagent. Work in `experiments/282-activation-state-dependence/`. **Edit nothing outside it.** No
build, no `vcc`, no `cell-eval2`, nothing under `submissions/` or `data/`; no download — 280 already extracted the rows
(`experiments/280-gwcd4i-donor-replicates/work/bank280_X.npy` 1.2 GB, `bank280_rows.npy`, `obs280.npz`, `axis280.npz`,
`extra280_*`; read them backed / memmapped). `uv run python`. Peak RSS under ~8 GB. **Time box: 4 hours.**
`docs/RESEARCH-PROGRAMME-2026-09.md` binds you: preregister before looking; name the unit; strong baselines.

## Why (F8–F10; R280 §2–§6)
Every cross-line comparison in the programme confounds cell state with genotype, culture history and technical
platform. GWCD4i measures the same 4 donors' CD4 T cells in three activation states, in one lab, one chemistry (10x
Flex), ~500 cells per (target, state). 280 established the replicate ceiling across donors within a state (+0.054) and
that the donor is a marginal unit. **So the across-STATE conservation of a target's response, measured against the
across-DONOR reproducibility within a state as the ceiling, is the purest available estimate of how much of a
perturbation response depends on cell state — and, per target and per family, which responses are state-dependent.**
This is the organising question asked within one person.

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; `docs/FINDINGS-2026-09.md` F8–F10; 280's `RESULTS.md` §1 (what was extracted:
donors, conditions, guides, panel targets; the residual space; KA rows), §2 (the across-donor statistic and the null
construction — REUSE `extract280.py` / `desc280.py` / `bench280.py`'s functions and the exact residual geometry), §7
(files); 279's `work/znf279.csv` (family labels and the cross-line conservation of the same targets — the comparison
axis); 276's `DEFINITION-276.md`.

## Protocol (write `PREREGISTRATION.md`, hashed, before any across-state number)
1. **Objects:** per (donor pair, state, target) residual effect rows exactly as 280 built them (each state centred on
   its own across-target mean, its own controls). States: Rest, Stim8hr, Stim48hr. Targets: the 297 panel targets
   present (and, if `extra280` carries non-panel targets, those too, labelled).
2. **Statistics, frozen:** (S1) across-state cosine per target: cos(row_state_i, row_state_j) with BOTH sides from the
   same donor pair (so donor variation is held fixed), for the three state pairs; (S2) the ceiling: across-donor
   cosine within a state (280's +0.054 family; recompute it here as the KA); (S3) the state-dependence index
   `SD = 1 − S1 / S2` per target (with S2 > 0), and its paired bootstrap; (S4) the family / class contrast: SD by HGNC
   family and by 279's cross-line conservation tercile — prediction (frozen): targets with low cross-line conservation
   have high within-individual state dependence (Spearman(SD, C_cross) < 0, interval excluding 0). **Null:** target-
   shuffled state pairings; a "same state, different donor pair" comparator to show S1 < S2 is not a donor artefact.
3. **The shared vs state-specific decomposition** (as 268 did across lines, now within individuals): per target, the
   fraction of response energy shared across the three states vs state-specific-and-reproducible (across donor
   pairs) vs noise; pooled and by family; the effective rank of the reproducible state-specific part.
4. **Descriptive (post-hoc, labelled):** the 20 most state-dependent targets and the 20 least; whether the
   state-dependent responses are the immune-programme genes (activation, cytokine, proliferation — say how you label
   them) or scattered.

## Deliverable
`PREREGISTRATION.md` (hashed), `RESULTS.md`: the KA (S2 reproduces 280), the S1/S2/SD tables with intervals and nulls,
the family/tercile contrast with its pass rule, the decomposition, the descriptive tail, and one line:
**STATE DEPENDENCE WITHIN INDIVIDUALS: median SD … [lo, hi]; PREDICTED BY CROSS-LINE CONSERVATION: YES / NO (ρ …)**.
Write `work/state282.csv` (per target). No model, no submission. Then stop.

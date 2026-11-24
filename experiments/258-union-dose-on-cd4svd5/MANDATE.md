# 258 — The six-source union on top of `cd4svd5`: is there a dose or a source subset where its `pds`+`mse` gain beats its `nmae` cost? Six-member price on the instruments the board just calibrated. Measurement only. DO NOT BUILD, DO NOT SUBMIT.

You are a research subagent. Work in `experiments/258-union-dose-on-cd4svd5/`. **Edit nothing outside it.** No build,
no `vcc`, no `cell-eval2`, nothing under `submissions/`. `uv run python`. One PLOLO fold at a time, peak RSS under
~10 GB, never alongside another ≥ 15 GiB job (check `ps`; a sibling subagent, 257, may be running PLOLO folds too —
that is fine, but stagger anything above ~6 GB). Record peak RSS per script.

## Why (2026-09-05 00:28Z board point, `experiments/237-gate-validation/BOARD-RESULT.md`)
`ARM6` — 229's `stack6`, the six-source union at g = 0.35 in place of the H1 pair, on the `h1union` base — scored
**0.2081907, −0.0036403** vs `h1union`. But by member it **bought `pds` +0.00248 overall (raw +0.00666, 1.9× the
`h1union` step itself) and `mse` +0.00088 (raw 0.81227 vs 243's point prediction 0.81244: the new-source class rule
confirmed, realised k_cos 0.962)** and **paid `nmae` −0.00716** (the 216/188 exact identity predicted −0.004783 —
right sign, 1.5× magnitude; the identity is now VALIDATED at scale and is a conservative bound). The four DE-three
members were flat. So the union's *directions* are worth +0.0034 and its *magnitude on the Wilcoxon-significant
genes* costs twice that. **Question: is there a dose g, or a subset of the six sources, where the gain exceeds the
identity's `nmae` leg (×1.5) by at least the +0.0019 gate — measured on top of `cd4svd5`, the incumbent (its CD4
deflation changes the base the union rides on)?** 229's ARM3 (H1 pair + `dld1_lowmoi`) read +0.00067 and 248's
DLD-1 rider +0.00216 on the old base and the old, un-calibrated factors — re-price them here.

## READ FIRST
`docs/READ-FIRST.md` (v2.3), the two `BOARD-RESULT.md` of 237 and 249, `experiments/229-stack4-build/six229.py` and
`RESULTS.md` (the six-member composite and the registry of the six sources), `experiments/243-transfer-factor/RESULTS.md`
§4.4 (class rule) + `work/arm6_243.json`, `experiments/248-second-residual/RESULTS.md` (the DLD-1 rider and why it was
NOT BUILDABLE as a second residual — a single union only), `experiments/216-cd4-hard/` and `188-magnitude-calibration/`
(the `nmae` identity), `experiments/249-cd4-deflate/stage249.py` (how the `cd4svd5` base is staged), `AGENTS.md`
§8.1b, §11.1 rule 10, §14.

## The instrument and the price recipe (put it in `price258.py`; retrodict before you predict)
Six-member composite as `six229.py` does it, with the base changed to `cd4svd5` (stage the CD4 residual with
`VCC_CD4_DEFLATE_K=5` exactly as 249 did; KA: the staged effects reproduce 249's `a28d5be2…` — cite 255 KA 2 if you
re-derive it) and the factors calibrated by tonight's two points:
- `mse` leg: 197's identity at the arm's own delivered terms (243 §5.2, measured cosine, never the floor), with
  **k_cos = 0.96** for a new source into empty rows (ARM6 realised 0.962; band 0.94–1.00).
- `pds` leg: 2.2138 × k_pds × Δ`pds`_3cancer / 6 with **k_pds at 0.7 (conservative, 243's class value) and at the
  ARM6-realised value** — compute the latter yourself from BR237 (`pds` +0.00248 overall vs the arm's PLOLO Δ`pds`
  in `six229.json`) and report it; it is well above 0.7.
- DE-three: 0.24 × board Δcos (243 §4.3; ARM6 realised +0.00017 total — check the recipe against that).
- **`nmae`: the exact identity, ×1.0 (point) and ×1.5 (ARM6's realised ratio, conservative).**
- **Retrodiction rows, must pass first:** (R1) `stack6` g = 0.35 on the `h1union` base reproduces 229's `six229.json`
  row to the digit (`nmae_wald89` 1.0029112, `overall_nmae_vs_incumbent` −0.0047829, the `pds`/cosine numbers);
  (R2) its six-member price with the factors above lands on **−0.0036 ± 0.0010** (the board), and by member on
  `pds` +0.0025 ± 0.0008 / `mse` +0.0009 ± 0.0003 / `nmae` −0.0072 ± 0.0025 — state which factor set does it;
  (R3) `cd4svd5` vs `h1union` on the same recipe lands on +0.0043 ± 0.0008 (BR249; use k_cos 3.76 / k_pds 1.22 for
  that re-shape step as 257 does).

## The screen (all on top of `cd4svd5`, three cancer folds, per-fold rows, LOO σ)
1. **Dose:** the six-source union at g ∈ {0.10, 0.15, 0.20, 0.28, 0.35} replacing the H1 pair (229's construction).
2. **Subsets at the best dose and at 0.28:** the H1 pair (the incumbent's own union — the zero row), pair + `dld1_lowmoi`
   (ARM3), pair + each single other source, the four-source `stack4`, all six.
3. For every cell: the six-member increment over `cd4svd5` at conservative (k_pds 0.7, `nmae` ×1.5) and point factors;
   delivered amplitude via 243's stream if a full staging is cheap, else 246's delivered ratio, and 243's amplitude
   leg when it moves > 0.2 %; leave-top-3-out retention for the best cell (229 flagged ARM6's `pds` as concentrated).

## Deliverable
`RESULTS.md`: the three retrodiction rows first (abort and report if R1 fails; if R2 misses, say by how much and
which member); the dose and subset tables; the best cell; **whether ANY cell clears +0.0019 overall over `cd4svd5` at
the conservative factors**; the exact builder configuration it would need (`VCC_UNION_SOURCES`, `VCC_H1_UNION_G`, on
top of `VCC_CD4_DEFLATE_K=5`, α 0.982, f 0.34, `rc_callset` — all existing flags; if it needs a flag that does not
exist, name it and stop); and one line **CLEARS / DOES NOT CLEAR**. If it clears, `PREREGISTRATION-DRAFT.md` in
249's form. A well-supported negative — "no dose of the six-source union clears on the validated identity" — is a
valid and useful result: it closes the family. No model name, no description, no submission. Then stop.

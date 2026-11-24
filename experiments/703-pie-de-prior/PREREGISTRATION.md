# 703 — PIE as a per-target DE sign prior on the CPM block, prior-only (preregistration)

Written 2026-10-14 ~23:15Z by the science session [bf003b]. At this time no PIE prediction has been read, the artifact and
package do not exist, and no 10-15 entry exists. Slot: assigned by the coordinator [400b6e] after 702.

## The arm (one variable vs ALLSW, `658allsw`)
- **703piept** = ALLSW's exact env + `VCC_SAMELINE=<pie_sameline.npz>:<sha256>:0:pt:noprog`.
  - **w = 0:** `sameline_blend` returns the matrices unchanged, so the bulk rows and design rows are ALLSW's, and pds/mse and
    the design are untouched.
  - **pt:** for each covered target, the DE block prior is PIE's sign wherever the artifact's `signs` are non-zero, else the
    context default (`sl_prior[t] = where(signs != 0, signs, de_sign_full)`, builder 451).
  - **noprog:** the default stays CLEAN4 q80.
- Artifact: `make_pie_sameline700.py --seen work/pie_train_perts.json --sign-pde 0.5`. Covered targets are the SEEN ones,
  i.e. target ∈ PIE's pinned train-split perturbations (154 / 300 on A/B/C), a rule computable for any panel on the day.
  Signs = sign(PIE lfc_pred) where PIE p_de ≥ 0.5 (PIE's own DE threshold; not tuned), else 0, which falls back to CLEAN4.
  Per context: A = jurkat, B = CVCL_0030 (HeLa), C = CVCL_1108 (CAL-33), each with its own control mean and Cellosaurus text.
- Expected config diff vs ALLSW: the same-line keys only (`sameline` record: path/sha/w/pt/noprog), plus the audit block.
  pds/mse must equal ALLSW's to ≤ 1e-6 raw.

## Class and noise
**Transfer** (a better sign prediction). It is a prior-only DE arm that keeps pass 1, which consumes no random draw
differently (468), so the board read is near-exact (≈ 1e-5 overall; AGENTS rule 13).
Held-out reads (PIE-clean, from `eval700_h1.py` and `eval700_k562.py`): sign agreement with the truth's lfc on OUR call set
(head | CPM block, adjudicable genes) for "PIE where p_de ≥ 0.5, else CLEAN4" vs CLEAN4, on h1-test (seen subset and all).
On K562 the bulk delta's sign is the proxy, since the GWPS table carries deltas, not DE calls. Sizing rule (625): ≈ +0.002
overall per +0.01 of gate-sign gain on route genes.

## Primary statistic and bins (board, vs ALLSW's own entry)
**Primary: Δ(fid + reach + jac + nmae) in overall units** (the four DE members, scaled, / 6).
- **TRANSFERS, ≥ +0.0008:** adopt `pt` with PIE signs for October. The day-of path needs PIE inference on D/E/F + the artifact
  + VCC_SAMELINE in def_day, which is a def_day pass-through and a rehearsal by 10-21 ([17adab]).
- **WEAK, +0.0002 … +0.0008:** adopt only if h1-test's held-out sign gain is positive on the seen subset.
- **NULL, −0.0002 … +0.0002:** not adopted.
- **HARM, ≤ −0.0002:** not adopted. Record the held-out over-read; reach is the likely loser, as in 458d.
- Falsifiers (stop and diagnose before reading): pds moving by more than one rank flip (3.716e-6 raw = 1/(3·300·299);
  ALLSW itself differs from 498p2g by one such flip from documented rounding swaps), or mse by more than 1e-6 raw; a config
  diff beyond the same-line keys. (Amended 10-14 23:10Z, before any PIE number, on the coordinator's ALLSW note.)

## Predictions
Filled in by an amendment after the held-out reads and before the build: the held-out sign-agreement gains and the
625-rule board prediction. Tries: 1 board arm on this switch, 1 threshold (0.5, not tuned).

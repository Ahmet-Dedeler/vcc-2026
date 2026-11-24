# 162 HANDOFF — panel-restricted LOLO validator (PLOLO). Stopped mid-sweep.

## 1. Scripts (all in `experiments/162-retrodicting-validator/`)

- **`board.py`** — RAN CLEAN. Ground-truth board record, transcribed from primary sources
  (each entry carries a `src` field). `python board.py` prints per-arm, per-member realised
  deltas and asserts `mean(6 member deltas) == overall delta`: **PASS for all 12 arms**, so
  the transcription is internally consistent. Also classifies each arm's channel
  (geometry / emission / mixed) — the emission arms are the validator's known-answer zeros.
- **`plolo.py`** — RAN CLEAN. Fold builder + renderer. 4 folds, truth = held donor library's
  own rows over the real 300 panel targets; sources = production donors minus the held line;
  leak-free operator per fold (062 held-safe trio, PRODUCTION `union_pruned` for kolf2 since
  kolf2 is not a production donor). Ports production `mass_consistent_effect` and
  `profile_from_effect` verbatim, renders and reads back like the scorer, and can render in
  either the **library frame** or the real **challenge A/B/C frames**.
- **`run_arms.py`** — RAN, ~35% complete when stopped. Sweeps 10 reconstructed arm configs
  (exact `alpha` / `gamma` / `lam` / `cd4` / `mass_fix` / `whitener_id` read out of the
  shipped `submissions/*.meta.json`) × 4 folds × 5 spaces. Resumable: it appends to
  `work/arms.json` and skips `(fold, arm, space)` already present.
  **`retrodict.py` was NOT written yet** — that is the deliverable table's generator.

## 2. Measured so far

Folds (from `work/fold-*.npz` audits): panel coverage / scored genes / op sha —
`orion_hct116` 290/300, 17806, `7cecbf7f`; `orion_hek293t` 294/300, 17806, `97e486e0`;
`k562_gwps` 269/300, **7448**, `13d988fc`, **pm_supported_targets 0**;
`kolf2` 273/300, 17806, `e2f65ecf` (production). Donor↔truth target overlap is recorded per
fold in the audit (that is same-target-different-line, not leakage).

Known-answer rows, both PASS: `a15-dup` reproduced `a15` bit-for-bit; `intended` beats
`library/delivered` (0.803937 vs 0.789536 on `orion_hct116`), same sign as 156's
+0.00324/+0.00382. **Challenge-frame reproduction is the strong one**: mass ratio
1.00178/1.00174/1.00171 and render residual **2.88 / 2.90 / 2.87 %** on A/B/C, against
`A23-PREDICTION.md`'s independently measured 2.83–3.35 %.

**a23 (mass_fix) retrodiction, the only complete one:** on `orion_hct116`, challenge-frame
mean raw pds 0.787706 → 0.791823 = **+0.004117**; library frame **+0.010607**. Board raw pds
moved **+0.00075** (scaled +0.001737). So sign correct, magnitude over-predicted **5.5×**
(challenge frame) / **14×** (library frame). Challenge frame is the better instrument.
`sgn_w` fell 0.492 → 0.417 under a23 while board `nmae` was flat — treat `sgn_w` as unvalidated.

## 3. Next step

Finish the sweep (`uv run python run_arms.py`, resumes), then write `retrodict.py`: pair each
arm with its base per `board.ARMS`, convert Δraw pds via the board map (span 0.45208), fit
Δ`mse` against Δ`sse` with leave-one-arm-out, and emit the predicted-vs-realised table + RMS
residual. `a13` (shrinkage) and the four emission arms are not reconstructible in this
pipeline; the emission arms enter as exact zeros (board confirms a21/a22 pds/mse = 0.0).

## 4. Traps hit

1. **`reach`/`fid`/`jac` are NOT priceable by this design.** The refde tables
   (`experiments/098-nmae/work/refde_*.npz`, `099-donors/work/refde_hek293t.npz`) exist but
   cover only the harness's 85/132 targets — a donor library gives pseudobulk effects, never
   cells, so a panel-restricted validator has no reference DE pool. This is the answer to
   "does it get `reach`'s sign right": it cannot produce a `reach` number at all.
2. E062's `enforce_memory_limit` reads **peak** RSS (`ru_maxrss`), so loading a 2.56 GB
   float64 operator permanently poisons later folds' guards. Fix used: pre-build every fold
   in its own process (cached), then run arms.
3. `var/_index` in `data/challenge/context_*.h5ad` is a nullable-string **group**
   (`_index/values`), not a dataset.
4. `kolf2`'s library-frame CPM gate passes 18,106 genes (its `b_ctrl` covers a gene subset
   and renormalising inflates CPM). Another reason to read only the challenge-frame rows.

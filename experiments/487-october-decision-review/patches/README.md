# 487 proposed patch (NOT applied; the coordinator decides)

## `in_mix_sameline.diff` — `scripts/def_sameline.py` (base `6024466d`, sha256 prefix `41d76266b3259d73`)

**What it does.** A new opt-in flag, `--in-mix-weak`. A matched context whose covered rows all come from a production
donor or a same-lab sibling (`IN_MIX_BANKS = k562_gwps, k562_essential, orion_hct116, orion_hek293t`) gets:
- `w_<c>` = `--w-weak` (0.5) instead of `--w-strong` (2);
- all-zero per-target signs. The builder's `pt` then falls back to the context's default prior on every gene
  (`np.where(_s != 0, _s, de_sign_full)`, `make_submission_pm_cd4.py` ~l. 5527).

**Why.** RESULTS.md §1. For an in-mix bank:
- the direction lift is ×1.16–1.29 with a production-like whitened composition, against Pan's ×4.37 on H1;
- it is flat from w 0.5 to w 2;
- it sits below the weakest bank 450f priced, where w 2 ≈ 0 with mse worse than base;
- the bank's own signs beat CLEAN4 by −0.017 (K562, VIPerturb truth) and +0.038 (HEK293, XuCao truth).

**Known answers (run 2026-09-26 ~22:45Z on scratch copies with ROOT pointed at the repo; outputs in the session
scratchpad):**

| run | result |
|---|---|
| flag off, 465's inputs (ef-bundle-r, 465's identity, `--programme-source h1_2025_full`) | sha256 `b4eb2a6d…` = 465's recorded `sameline.npz`, **byte for byte** |
| flag on, same inputs (D = Jurkat → `nadig_jurkat`, not in-mix) | every array identical (NaN-aware); only `provenance` gains `in_mix_rows` / `in_mix_rule` |
| `--force D=HCT116_LARGE_INTESTINE`, flag off vs on | exactly two arrays differ: `w_D` 2.0 → 0.5 and `signs_D` (5,171,487 non-zero → 0); rows, targets, programme identical |

Peak RSS 3.3 GB (streams Orion's 2.7 GB bank for its programme).

**To adopt:**
1. Apply the patch.
2. Pass `--in-mix-weak` in `def_day.sh`'s stage 2. Its default-off form leaves every rehearsal's artifact unchanged.
3. Note in DEF-RUNBOOK's decision line.

**A builder follow-up the patch cannot do alone.** With the programme default on (`DEF_DAY_PROGRAMME_SOURCE`, flags
`pt`), the prior gate hands CLEAN4 only to weak targets *without* a same-line entry (`weak_kept_pt`). The builder still
adds an `sl_prior` entry for a covered target with all-zero signs, so an in-mix context's covered weak targets keep the
programme, not CLEAN4. For a covered target WITH signs, keeping its same-line prior over CLEAN4 is 468's deliberate
choice. It is defensible for a new-information bank; for an in-mix one it is not.

The one-line fix, in `make_submission_pm_cd4.py`'s pt loop (~l. 5525):

    if SAMELINE["pt"]:
        for _k, _ti in enumerate(...):
            _s = ...
            if np.any(_s != 0.0):          # 487: a target with no same-line signs keeps the gate's decision
                sl_prior[...] = np.where(_s != 0.0, _s, de_sign_full)

Its known answers:
- KA-1 (unset rebuilds 446b) and KA-2 (458d) are untouched: no pt entries there, or all non-zero.
- 465dr's artifact changes only if some covered target had all-zero signs. Check that before rebuilding.

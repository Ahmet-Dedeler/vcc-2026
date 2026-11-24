# 460dr — the October path with the LODO stack key, vs 453dr (preregistration, 2026-09-26 10:22Z)

Written before 446bj2 (446b + the LODO key, 09-27 s1) is scored, before 460dr's package exists, and before any 09-27…
09-30 entry exists.

## The arm
- 460dr = 453dr (`scripts/def_day.sh` on 180's relabelled A/B/C bundle, q 1.75, `auto`, `pt:noprog`) with
  `DEF_DAY_LODO=1`. The stack dump is scaled by 415's LODO-V1 key, re-fitted on the day's operator and CD4 banks under
  October rules (`d9a6b4c4`; w Pearson 0.9999988 with 446bj2's `8c82bef3`, max |Δw| 0.0067).
- Expected config diff vs 453dr: exactly `bulk_load`, `bulk_load_sha256`, `bulk_ctx_scale`. The re-level re-converges on
  the new dump. `gate460.py` checks this.
- It is the October code path's known answer for the key. It shows that the day-of fit, the second pin and the
  dump-only export reproduce on the board what 446bj2 measured on the A/B/C path.

## When it goes
Only if **key_ok** on 09-27: 446bj2 − 446b ≥ +0.001 overall with pds Δ ≥ 0 (433's PREREGISTRATION-0926 bins). Otherwise
the package is never submitted: the key stays out of October and there is nothing to confirm.

## Prediction and bins (460dr − 453dr, raw)
- The bulk key differs from 446bj2's by < 1 % in w, and the stacks differ only by CD4 ×1 (453) vs ×1.5 (446b).
- **Prediction:** pds and mse move in the same direction as 446bj2 − 446b, at 0.7–1.3× its size. DE members stay within
  ballast noise (nmae |Δ| ≤ 0.0015, reach |Δ| ≤ 0.010).
- **Bins:**
  - pds and mse both in [0.7, 1.3] × 446bj2's Δ, same sign: **the October LODO path reproduces the board key.** Set
    `DEF_DAY_LODO=1` as the day's default (DEF-RUNBOOK).
  - Same sign but outside the band: it transfers at a different size. Keep it on if overall ≥ +0.001, and record the ratio.
  - Opposite sign on pds or mse: a path defect. Keep `DEF_DAY_LODO` off until diagnosed.
- Name: a fresh joke at arm time. Description: empty.

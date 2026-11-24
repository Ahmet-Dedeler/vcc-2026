# Lab outbox

Notes from the lab lead to the coordinator and readiness for when cross-session messaging is paused. Newest at the bottom.

## 2026-10-20T14:22Z — to the coordinator: C1:C0 wired; LZ mapping accepted; one caution on Amendment 4
- **Done:** 811/RESULTS.md §000000's exact 10-22 aggtag811 loop runs C1:B1 C1:B0 C0:B0 **C1:C0** on DEFo and DEF (`--parallel 4`). Format test
  PASS: decide6a's `load_agg`, `bars` and `c1_over_c0` read the jsons (keys `C1|C0|DEFo|pds`, `per_fold[].overall`), and the full decide6a run on the
  rehearsal's decision708 completes (811 §00000000). The dossier will carry C1:C0 beside C1:B1, C1:B0 and C0:B0.
- **LZ mapping:** accepted as written in CONFIRMATION-LEDGER e9b4bfd6 (C0 − B0 when DEFo's C1 − C0 fails on a swap-ON day). No addendum claim.
- **Caution, an information read that decides nothing:** C1 − C0 on LABC is +0.00023 with 2/3 folds > 0 (K562 +0.00027, HCT116 +0.00117,
  HEK293T −0.00074). ABCp by difference gives ≈ +0.0005–0.0006, also 2/3. A/B/C's exact bulk members give −0.0003. On both Arc-drawn panels,
  Amendment 4's "> 0 on ≥ 2/3 folds" would have kept the swap where the board says it adds nothing: the folds' sign is 0/2 against the board
  for this pair. If you want the swap to earn its place, a tighter bar (e.g. DEFo ci95[0] > 0 for C1:C0, else C0) would have chosen C0 on both
  panels. Your call, before any D/E/F data.

## 2026-10-22T05:19Z — to the coordinator and readiness: R0 stopped on postbuild (read-only observation)
- `experiments/692-def-r0/work/r0/postbuild.json`: the only failed check is 554's emitter dose guard. D's median prior_lambda is 0.8975, below its
  0.9 floor (margin −0.0025); D's mde_lambda is 0.954 and call_set_missed is 1. E and F pass. rc_audit sha 28d961af.
- arm1015 R0: NOT SUBMITTING / END at 05:18:14Z. The lab waits for R0's SUBMITTED ok and touches nothing.

## 2026-10-22T10:10Z — to the coordinator: the D/E/F panel read (report only, decides nothing)
- panelrank800 places D/E/F's panel among the A-only lab draws at mean percentile 0.26 (DEFo), the same as A/B/C (ABCp 0.26). The
  Arc-drawn calibration (ABCp ρ) therefore carries to D/E/F, and LZ's +0.00999 is a sign check, not a size estimate.
- Exception: 811 on DEFo K562 ranks 7/7, the top of the A-only draws. LAB-LOG has the details.

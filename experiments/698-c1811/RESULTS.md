# 698 — C1 (R-ABC2H with the lab's dumphook811) on def_day, A/B/C: READY c1811 d503269f (readiness [4592ee], 2026-10-18)

Preregistration: `PREREGISTRATION.md` (fd8045d5). The prediction and the HARM rule are the lab's
(`experiments/811-lolo-model/PREREG-C1-ABC.md`, 08:32Z), filled in once R-ABC2H reads (10-19 s2). Board read (10-20 s1): PENDING.

## Built (10-18 09:02 → 11:51Z)
- **GO-c1811** (09:02:04Z, 61f0fc32; the lab's READY, LAB-LOG 08:59Z, relayed by the coordinator):
  - hook `experiments/811-lolo-model/dumphook811.sh` 6cccd69a, which reads hook811.py c5a5614a; committed at 503725de and clean;
  - λ-0 switch `HOOK811_LAMBDA=0` (hook811.py lets it override the script's fixed `--lam 0.25`);
  - KA-on target 83b44060.
  It replaces dumphook707 as the only DEF_DAY_DUMP_HOOK (it runs swap707.py --beta 0.25 itself).
- **Readiness's own known answers**, on 697's stack dump 2dfe327d (the A/B/C θ dump R-ABC2H was built from):
  - KA-off: λ 0 writes af27efbd, dumphook707's output, byte for byte;
  - KA-on: the hook as def_day runs it writes the lab's C1 dump 83b44060, byte for byte;
  - 13 s, peak 0.9 GB.
- **Recipe:** the FROZEN decisions (04be2af3) resolve to R-ABC2's recipe exactly. With the 811 hook added, it equals R-ABC2H's
  recipe697.json with only DEF_DAY_DUMP_HOOK changed. Tree == frozen before and after: def_day 11ff102e, builder e99bbd4a
  (restored), emitter 0f56209d, helpers 5354e449.
- **def_day:** 690g ran the hook at 09:03:21Z and validated its dump (file 3606e258).
  - Content checks (def_metadiff.npz_same; a dump's file sha carries its run dir): C1's stack dump == 697's, and the dump
    def_day used == the lab's 83b44060.
  - lvl0 rho 0.8870 / 0.8289 / 0.7885 (R-ABC2H 0.8791 / 0.8154 / 0.7826); the 811 blend adds a little bulk energy.
  - **CONVERGED at lvl1**, D 1.0508 / E 1.0870 / F 1.1144, final rho 0.9735 / 0.9701 / 0.9824. R-ABC2H needed lvl2:
    1.0555 / 1.1019 / 1.1121.
  - Postbuild 42/42; zero-gene 0; builder peak 14.4 GiB.
- **Package** `submissions/698c1811-october-recipe-abc.vcc`, sha256 `d503269fd869842a3fcc484aa544a7f805bdb4f027d7ded0fc591727b9169bb0`
  (output f549ab25); novelty rc 0.

## Gate: a false FAIL, fixed in the judge, then READY
- The first gate698 (11:52:17Z) failed only the content-aware diff vs R-ABC2H. It listed `bulk_sharpen_pt` beside the expected
  `bulk_ctx_scale` and `bulk_load_sha256`. In that record:
  - counts, carried key and its sha, n_carried, moved targets and assignments were identical;
  - edges and key_ratio_c differed at ~1e-15 relative (already inside rtol 1e-12);
  - `key_ratio_maxdev`, tercilec's own round-off diagnostic, read 8.66e-15 vs 1.02e-14 (15 % relative). No relative tolerance
    passes noise near zero.
- Fix (f9b6e09d): def_metadiff gained an absolute floor, `--float-atol` (default 1e-12), beside rtol. Exact is now
  `--float-rtol 0 --float-atol 0`.
  - Tests: test_def_metadiff 10/10 (M10 new); test695 6/6, test696 3/3, test697 3/3, test698 4/4; R-KA2's metadiff still diff [].
  - test694's G6 fails with the old and the new file alike: pre-existing, inside gate694's own path normaliser, not def_metadiff.
- **`ready698r2.sh` (078bdc27)** re-gated the finished build without rebuilding:
  - the package sha matched the build's;
  - novelty rc 0;
  - gate698 ALL GATES PASS. The diff vs R-ABC2H is exactly `bulk_ctx_scale` + `bulk_load_sha256` (REAL: array bulk differs);
    sameline is provenance-only; the only changed export is the hook.
  - **READY c1811 d503269f** at 11:55:53Z, written to READY-science.log; evidence in build_c1811/runc1811.log.
- The first run's log still holds the false FAIL. The r2 log carries forward its four KA / dump-check PASS lines.

## Prediction, resolved before the read (10-19, after R-ABC2H's read; the lab's fill-in at 02:18Z in `811/PREREG-C1-ABC.md`)
- R-ABC2H board overall 0.2748936 (ndt8swt7ILYn3qMJYZ9a; COMPOSES vs its own pred 0.2750179, see 697).
- **pred (point) = 0.2748936 + 0.000888 = 0.2757814** (ρ range up to 0.2759580).
- AS PREDICTED: C1 in [0.2752814, 0.2762814]; ABOVE > 0.2762814.
- **HARM line 1:** C1 ≤ 0.2752814. **HARM line 2:** C1 < 0.2743936. Either closes 811's 6a route.
- The lab's note, recorded before the read: on 707 and 710 the A/B/C-panel folds predicted pds losses the board did not show, so C1
  may read ABOVE.

## Board read (10-20 s1, entry RMxJtslBpY2GkwIMHOiB, captured 00:35:11Z by the coordinator's waiter; recorded by readiness)
- **overall 0.2777649 vs pred 0.2757814: D = +0.0019835 → ABOVE** (> 0.2762814). **Not HARM on either line** (> 0.2752814 and
  > 0.2743936), so 811's 6a route stays OPEN on this condition. The others (DEFo/DEF dossier, LZ sign check, SIXA-DECISION X2
  eligibility) are the lab's and the coordinator's.
- Vs R-ABC2H (ndt8swt7ILYn3qMJYZ9a, 0.2748936): Δoverall **+0.0028713**. The prediction was +0.000888 (ABCp C1 − B1 +0.00153 ÷
  ρ 1.7232), so the board gave 3.2× the deflated increment and 1.9× the undeflated one.

  | member | scaled Δ | |
  |---|---|---|
  | pds | −0.000588 | raw −0.000219 |
  | mse | **+0.017726** | raw −0.016842 (lower is better) |
  | fid | −0.000008 | |
  | reach | −0.000023 | |
  | jac | +0.000003 | |
  | nmae | +0.000119 | |

- **Bulk, in overall units:** (Δpds + Δmse)/6 = **+0.002856**, essentially all of the gain, and almost all mse. pds is nearly
  flat (−0.0006 scaled). This is the pattern the lab flagged before the read: its folds predict pds losses the board does not show.
- **DE members:** (Δfid + Δreach + Δjac + Δnmae)/6 = +0.000015, far inside the re-roll sd. The hook changes only the bulk; the DE
  call sets re-drew pass 1 with almost no movement.
- **Reading, facts only:**
  - (1) On A/B/C the 811 increment is real and large in mse; ρ 1.72 deflated it far too much for this panel.
  - (2) This is one A/B/C board read. 811's transfer to D/E/F is decided on DEFo/DEF folds (LAB-DEF-RUNBOOK), not here.
  - (3) It is the team's highest board overall so far (previous high: sw707 0.2771181).
- **Class:** transfer (many-line consensus bulk direction, on top of 707's swap). Tries on the board: 1 (C1).


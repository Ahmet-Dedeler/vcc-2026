# 673 — W1 (R-ABC2 + DEF_DAY_H1_STACK=x2.5) on def_day, A/B/C: built by readiness, READY w1 1337cc9b after the coordinator's re-gate

Preregistration: `PREREGISTRATION.md` (coordinator, 5c4f4bc5 + 1182d686): pred 0.2744673; BELOW (D ≤ −0.0005) blocks any 708 H1 move
above ×2 (route: wait708 re-armed with WAIT708_R710=HARM before the release, OCT22-CHECKLIST 4b). Board read (10-21 s1): **0.2746084, COMPOSES** (D +0.0001410 vs pred); wait708 stays at R710=POSITIVE.

## Built (readiness [1efb7d], w1_673.sh 8fe13218; 10-20 07:04 → 11:15Z, right after C0)
- **Recipe:** the frozen decisions (04be2af3) resolve to R-ABC2's recipe exactly, with no DEF_DAY_H1_STACK (the inert-x2 KA: h1_stack
  DECIDED x2 = no export; R-KA2 877e7c28 rebuilt R-ABC2 byte for byte from it). + DEF_DAY_H1_STACK=x2.5 is the only addition. Tree ==
  frozen before and after; builder restored e99bbd4a.
- **def_day:** H1 stack x2.5 → dump `dump_h1x25_cd4x1_s0.25.npz` with VCC_H1_UNION_G 0.7 (checked in the dump's own env).
  - lvl0 rho 0.8881 / 0.8202 / 0.7902 (R-ABC2 0.8673 / 0.7994 / 0.7766). lvl1 0.9745 / 0.9700 / 0.9912.
  - **CONVERGED at lvl2:** D 1.0501 / E 1.0927 / F 1.1070, final rho 0.9745 / 0.9700 / 0.9787.
  - Postbuild 42/42, zero-gene 0, peak 14.0 GiB.
- **Package** `submissions/673w1-october-recipe-abc.vcc`, sha256 `1337cc9b489121e0782616f2ebccd26b6b7acc18e0bd8cbf4360480d90ac40e1`
  (output 621371bd); novelty rc 0.

## Gate: gate673 FAILED on its own spec; the coordinator re-gated the finished build (81b3f314)
- gate673 allowed {bulk_ctx_scale, bulk_load, bulk_load_sha256}. The diff vs R-ABC2 also listed **bulk_sharpen_pt** → NOT READY at
  11:15:45Z. That NOT READY also stopped W2's driver (674).
- **What moved:** only tercilec's data-derived scale fields (edges, key_ratio_c, key_ratio_maxdev). Every behavioural field is
  identical: mode, q, counts [99, 103, 98], carried key + sha, n_carried 294, the per-target moves (AKAP1, PIM2, SLC25A13), and the
  floor/tie/non-carried counts. Each target is sharpened exactly as in R-ABC2.
- **Readiness's check of the cause:** edges and key_ratio_c all scale by ONE common factor, **0.97733018428809 (to 1e-15)**:
  0.19650/0.20106, 0.63821/0.65302, 0.39300/0.40211. That is a pure re-normalisation. tercilec's key ratio is normalised on the
  stack dump, which the H1 dose rewrites by design. Neither the coordinator's spec nor gate673 anticipated it: a spec gap, not a
  bad build.
- **The coordinator's `regate673.py` (81b3f314):** it accepts bulk_sharpen_pt only if the behavioural fields are identical and only
  the data-derived ones moved, with exactly one FAIL row (the config diff) and every other row PASS. **READY w1 1337cc9b** at
  11:32:45Z into `experiments/670-coord-1015/work/regate673.ready.log`, the 10-21 s1 waiter's ready-log (arm1015 1021s1 pid 92496,
  name 'social media isnt social at this point, its just media.').
- **For any later weights build** (a 6a X1 with an H1 or CD4 move): expect bulk_sharpen_pt's scale fields to move. Judge it as
  regate673 does: behavioural fields identical and the scale fields moving by one factor. OCT22-CHECKLIST 6a (a) was fixed by the
  coordinator for this.

## Board read (10-21 s1): 0.2746084, COMPOSES (entry WULHAm2fytPtmAb3ryVN, rank 88, submitted 00:01:52Z)

Read by the coordinator's `read695 --w1` at 00:53:43Z (w1route log `experiments/670-coord-1015/work/w1route.log`; record
`data/board-entries/WULHAm2fytPtmAb3ryVN.json`; ledger b7c0d14c). Pred 0.2744674 (R-ABC2 0.2739849 + 710's H1 ×2 → ×2.5 step
+0.0004824), read 0.2746084, **D +0.0001410**: inside the COMPOSES bin (BELOW needed D ≤ −0.0005). W1 − R-ABC2 = +0.0006235 overall.

| scaled member | W1 | R-ABC2 | W1 − R-ABC2 |
|---|---|---|---|
| pds | 0.8095617 | 0.8091907 | +0.0003710 |
| mse | 0.3841555 | 0.3825495 | +0.0016060 |
| fid | 0.0055074 | 0.0054189 | +0.0000886 |
| reach | 0.2354233 | 0.2340327 | +0.0013906 |
| jac | 0.0284677 | 0.0284616 | +0.0000062 |
| nmae | 0.1845346 | 0.1842562 | +0.0002784 |

- **Class:** transfer (the H1 donor weight in the stack). This is its A/B/C board read. The held-out read is 708's day-of re-pricing on
  D/E/F's own panel, which decides whether H1 moves above ×2; W1 only guards that move.
- **Route:** no block. w1route left wait708 (pid 5079) at R710=POSITIVE, armed line 2026-10-20T15:37:31Z; the coordinator asked to
  keep it so.
- **Noise:** this changes the bulk, so pass 1 is re-drawn: re-roll class, sd ≈ 0.00025 overall (AGENTS 11.1 rule 13). The +0.00062
  over R-ABC2 is about 2.5 sd, and almost all of it is bulk (mse +0.0016, pds +0.0004) plus reach +0.0014. The read
  sits +0.00014 above its own prediction, under 1 sd: W1 composes with R-ABC2 as 710 predicted.
- **Tries:** one W1 arm on the board; W2 (H1 ×2.5 + CD4 ×1.5) was preregistered but not built (674).

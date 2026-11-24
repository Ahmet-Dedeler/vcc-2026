# 549 — arm 446bj2zcph built and packaged (results, 2026-09-30T22:03Z)

**Arm:** 446bj2z + 528's call-propensity block (446bj2zcp's A/B/C artifact `04fb24ec…`) + 496's CLEAN4 uncalled-gene hook
(`VCC_DE_MDE_SIGN=clean4:0.05`). Why and launch rule: `PLAN.md`. Built early (09-30 21:08Z, before 446bj2zcp's read) so
that 10-05/10-06, which are crowded, need no build for it; it is slotted only if 446bj2zcp reads ≥ PARTIAL.

## Build (`build549.sh`, 21:08–21:49Z, builder 4daae237 / emitter 0f56209d == HEAD)
- Artifact `cph549.h5ad` 11,960,933,560 bytes, sha256 `d8ca81bc031b2be8…`.
- `[446bj2zcph] PASS: exactly ['de_mde_sign', 'de_sign_emitter_sha256', 'rc_callset_block']` vs 446bj2z; zero-gene audit 0
  genes in 0 blocks; `CLIFF OK {"A": 0.9602, "B": 1.0281, "C": 1.0784}`, identical to 446bj2zcp's (a DE-only arm).
- 537's driver was paused (SIGSTOP, idle at the time) for the build and resumed at 22:01:50Z by the EXIT trap.

## Packaging (by hand, `prep549.sh`, 22:01–22:03Z)
- `run_arm_generic.sh`'s prep gate needs vm_stat free+inactive ≥ 23 GiB; from 21:46Z the Mac sat at 18.8–20.5 GiB because the
  compressor held ~18.5 GiB of other apps' pages (browser renderers), which free+inactive does not count, while
  memory_pressure reported 58–62 % free and swap was unused (diagnosis: peer session "Explore contribution areas", 22:00Z).
  The gate would have timed out at ~23:58Z, into the 00:01Z upload.
- The waiting `run_arm_generic.sh` (pid 93078, in `wait_gate`; builder already restored at 21:48:56Z; the script has no
  trap) was stopped; `build549.sh` exited "no package" and resumed 537's driver. `prep549.sh` then ran the arm script's
  exact prep lines into the same run log, with a `vcc prep --dry-run` first and the same gate minus the free+inactive term
  (pressure 58 % ≥ 40, no builder, no prep, disk ≥ 20 GiB). Peak memory footprint of `vcc prep`: 24.0 GB; no swap.
- **Package** `submissions/446bj2zcph-…-zs250-cph.vcc`, 2,771,148,800 bytes, **sha256 `de1e4ab950f08990…`**; content identity
  `859d8b22b27a14ff…`; novelty rc 0; meta copied beside it (`5044276a…`).
- **550's checks:** `config_keys550.py` vs 446bj2zcp's meta: exactly `de_mde_sign` (PASS). Bulk identity vs
  `callprop535.h5ad` (per-(context, target) summed counts): see `work/bulk_identity_vs_cp.log` (appended below).

### Bulk identity vs 446bj2zcp (`550/repro/bulk_identity550.py`, 22:03Z; `work/bulk_identity_vs_cp.log`, `diffgroups549.py`)
- **Not byte-identical, and not a bulk leak:** 895 of 900 (context, target) groups have identical summed counts; in 5
  (A MAPK12, B AKT2, B GRAMD1A, C STK3, C ZBTB48) exactly ONE count moved between two low-count genes of the group
  (+1 / −1; group totals equal): VPS9D1 30/29 & PLD2 50/51; ZNF365 114/113 & FRRS1 98/99; CRYBG3 65/64 & PDZD7 55/56;
  PRKAR2B 57/58 & FAM184A 49/48; ANKRD44 63/64 & SLC31A2 104/103 (cph / cp; group totals 12.2–23.6 M counts).
- Mechanism: the hook adds only CLEAN4-steered whole-count moves on uncalled genes (the only difference between the two
  arms); a move that lands on a neighbouring gene of the same cell keeps the cell's total and swaps one count between genes.
- Size: 10 gene values move by 1–3.5 % in 5 of 900 pseudobulks; expected pds and mse change far below the 1e-6 raw
  falsifier of the DE-only bins. 496a (the same hook on 446bj2, board read 10-05) measures it directly on the board.
- Decision: 446bj2zcph stays a candidate. Its PREREGISTRATION must state this expected tiny bulk move and keep the 1e-6
  falsifier; the read uses 537's bins.

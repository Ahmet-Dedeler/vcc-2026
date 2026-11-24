# 693 — why DNMAG doubles a build, and a byte-identical faster design (2026-10-15, session [17adab])

Class: engineering of the shipped emitter pass (no output may change). Asked by the coordinator after 690 §11 (the DN build
ran ≥ 139 min in its build stage vs 630arm's 74 min). Stage 1 (now): profile + design. Stage 2 only if DN reads TRANSFERS or
PARTIAL: a new module (dn_mag v2, new sha) whose KA is byte-identical outputs, plus an UNAPPLIED builder pin-bump patch.

## Profile (`prof693.py`; real H1 blocks from 630's fixtures, the shipped order zs → UP b_k100_q20r3p10 → DN d_m200)
Six blocks of seed 0, unprofiled: **6.04 s per block** (5.37–7.47), peak RSS 3.9 GB. A production build has 900 blocks
(300 targets × 3 contexts): ~90 min of DNMAG per build, which matches the DN build's overrun. Per block: ~1,650 called
DOWN genes selected, ~1,640 moved, 40–52k counts moved (quantum 5), 0 guard rounds, de_stats twice (0.048 s each, negligible).

cProfile of the slowest block (ASF1A, 8.9 s under the profiler):

| where | calls | cumulative s | share |
|---|---:|---:|---:|
| `rank_mag._take` (the sub-gate payback: one call per move, each an argsort over the ballast row's ~8k sub-gate genes) | 51,652 | 3.0 | 34 % |
| `rank_mag.ranks_of` on 1-element arrays (2 per carrier move + 2 per ballast row per move, in a Python loop) | 853,166 | 2.5 | 28 % |
| `dn_pass`'s own Python loop | 1 | 2.7 | 31 % |
| `de_stats` | 2 | 0.1 | 1 % |

## Design of v2 (byte-identical by construction; the KA decides)
1. **Payback by per-row heaps instead of argsorts.** `_take(row, need)` removes one count from each of the `need`
   largest current holdings, ties broken by the lowest index (stable argsort of −value over live genes), repeating sweeps
   while `need` remains. Sequential takes are NOT one aggregated take: [3, 1] taken 1 then 1 gives [2, 0], while one take of
   2 gives [1, 1]. So v2 replays each ballast row's takes in their original order, on a heap keyed (−holding, index) per
   ballast row. One take = pop `need` entries, decrement, push back the ones still > 0: the same genes in the same order.
   Rows are independent (a ballast row's sub-gate counts change only through its own takes, and the carriers only receive
   them), so per-row replay equals the global order. If a take ever needs more genes than are live (never at quantum ≤ 5
   with ~8k sub-gate genes), v2 falls back to `_take` itself.
2. **Batched rank lookups.** The carrier's two lookups (x−k, x) become one 2-element call. The ballast loop becomes one
   call over all eligible rows (current and +k). The choice is argmin over eligible rows in `bal` order, the first minimum,
   which is what the strict `<` loop picks. Eligibility is unchanged (`sub_bal ≥ k`, `cw[b] < cw[a]`).
3. **Unchanged:** selection, ordering, the capU rule, the quantum fallback to 1, the carrier top-holder rule, the
   invariant asserts, the guard rounds and the diagnostics. The asserts stay as the fail-safe.
4. **Risk to exactness:** the same float64 `log1p(x / s)` evaluated in a longer array. numpy's float64 log1p on arm64 is
   elementwise libm, and searchsorted is exact. The KA catches any ULP difference.

**Expected:** about −3 s (payback), −1.8 s (lookups) and −1 to −1.5 s of loop overhead per profiled block, i.e. ~1.5–2 s
per real block, 3–4× faster: DNMAG ≈ +20–30 min per build instead of ~+90.

**KA for stage 2 (the coordinator's conditions):**
- byte-identical blocks and diagnostics vs v1 on test_dn_mag631's 3 seeds × {DN-only, UPDN} × 85 targets;
- one real A/B/C block vs v1;
- then an UNAPPLIED builder patch that bumps DE_DNMAG_MODULE_PATH/SHA to v2, so the swap is one reviewed step.

## Timing for the day (690 §6 / §11)
With DNMAG ON at v1, a def_day level build gains ~90 min, so R0 (dump + 3 level builds) takes ~4–5 h instead of 1 h 41 m.
At v2 (if built), ~2–2.5 h.

## Stage 2 (GO after DN read PARTIAL, 10-15 ~01:50Z)
- **v2 written** (`dn_mag2.py`, sha256 `b9d6cfcf17ad783c…`): v1 plus the three exact rewrites above (a 2-element carrier
  lookup, one batched ballast lookup with the first minimum, per-row heap replay of the payback).
- **Smoke KA** (seed 0, UPDN, 3 blocks): v2 == v1 byte for byte (blocks and diagnostics), v1 17.7 s, v2 2.9 s, **×6.19**.
  So ~1 s per block, i.e. DNMAG adds ~15 min to a 900-block build instead of ~88.
- **Full KA running** (`ka693.py`, one serial ~4 GB stream, the light budget beside a running build): 3 seeds ×
  {UPDN, DN} × 85 blocks, records `work/ka693_s<seed>_<mode>.json`.
- **Pin bump 693p** (UNAPPLIED; `patch693p.py`, `test693p.py` 8/8 on synthetic records): installs src/vcc2026/dn_mag2.py,
  points the builder's module path, sha and both imports at it, and moves def_day's 690a meta check to v2. It refuses
  unless the six real KA records exist and name this v2. v1 stays on disk.
- **690h** (UNAPPLIED; `patch690h.py`, `test690h.py` 9/9): the day-of DNMAG delivery report (per context mean dn_frac
  ≥ 0.10 and mean dn_moved ≥ 200, from the final rc_audit dump; REPORT, never fatal; `690h: dnmag LOW_CONTEXTS …`;
  a warning before DONE). On the DN build's real dump: A 0.193 / 1,512, B 0.219 / 1,488, C 0.185 / 1,541, LOW none.
  691's watcher copies the lines into NEXT.txt after R0, and 692 logs them; both re-armed (pids in PLAN).
- The stray `scripts/def_pipeline.py.orig` was == cb9e432a's committed def_pipeline.py (a 09-27 patch backup): deleted.

## Stage 2 verdicts (10-15 03:24Z): every known answer PASS; 693p cleared for apply
- **H1 KA, 6/6 PASS** (`work/ka693_all.log`, `work/ka693_s*_*.json`): v2 == v1 byte for byte (blocks and diagnostics) on
  510 real blocks, 0 fail-safes. UPDN (the shipped order) ×6.72–6.78: v1 ~560 s vs v2 ~83 s per 85 blocks. DN-only ×5.0–5.1.
- **test693p --real-ka 6/6** on the six real records.
- **Real A/B/C blocks, PASS** (`abc693.py`, `work/abc693.json`; context A: STK3 (largest K, most moved), PHF23, CHKB (fewest
  moved, 765)). For each, v1 on 630arm's written block (DNMAG's exact input; n_car from the DN build's rc_audit; scorer
  controls rebuilt from context_A as the builder makes them) **equals the DN build's written block**, so the input
  reconstruction is faithful and this is the real day's object. Then v2 == v1. Peak 4.07 GB.
  - The first run died on a reader bug, not a mismatch: the restored `var/_index` is an encoded string group. Fixed by using
    restore_prediction_h5ad's `strings()`.
  - The second run picked STK3 twice (2 distinct blocks), so the selection was made distinct and the check re-run.
  - Both restores were deleted after each run with their sha256 logged (`work/chain693.log`; regenerable from the .vcc).
- Day-of consequence: with 693p applied, DNMAG adds ~15 min per 900-block build (≈ 900 × 1 s) instead of ~88 min.
  R0 ≈ 2 h with DNMAG ON.

## Build-level known answers on the v2 tree (ka693b, 10-15; re-run after the 10-15 sleep/shutdown killed the first attempt)
- **KAU2 PASS** (05:02Z, again 22:02Z): all knobs unset on builder e99bbd4a (dn_mag2 pinned) == 446bj2z's ec3a5ae6.
- **DNV2 PASS** (23:39Z): dn671's DN build (630arm + VCC_DE_DNMAG=d_m200) rebuilt on the v2 tree is recognised by the novelty
  checker as the already-scored DN entry wzaVLFTmlAMW3Mfw0dD8 (the ledger judges the prediction, not our own checks).
  Whole run 22:02→23:39Z = 97 min, against the v1 DN build stage alone of 162 min.
- So dn_mag v2 is proven byte-identical by three judges: 510 H1 blocks, 3 real A/B/C blocks, and the board-scored DN package.

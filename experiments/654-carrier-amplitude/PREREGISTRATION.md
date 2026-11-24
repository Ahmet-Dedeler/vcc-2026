# 654 — DE carrier design amplitude ×0.8 / ×1 / ×1.25 on H1's real emitter, by 631's regions (PREREGISTRATION)

2026-10-02, council session [7772dc]. Approved by the coordinator (low priority; range 654–655). [d4cd5f]: no overlap with v7.
Written and committed BEFORE any 654 run exists. Class: **scorer mechanics** (how the DE carriers are written; the bulk is
pinned by ballast). Nothing is fitted to A/B/C.

## Why
- 652 F6: `VCC_ALPHA` reaches the DE channel only as the carriers' design amplitude. The re-level scales only the ballast's
  bulk rows.
- 410 §7 (09-22 construction: no coarse MDE, no zs) measured carriers ×1.25 / ×1.5 / ×2 with the bulk pinned:
  - nmae +0.0025 / +0.0123 / +0.0287 raw;
  - reach −0.0105 / −0.0137 / −0.0107 raw;
  - fid ~0.
- Nothing below ×1 has been measured. Today's carrier amplitude was never tuned in the ballast era.

## Arms (one variable: the DE design `eff_de_native` × s; bulk = `ps_bulk(eff)` unchanged)
- `x080` (s 0.8; new design in `run654_496.py`, a copy of run496 that adds ONE design branch);
- `base` (s 1);
- `x125` (s 1.25; run496's existing known-answer design).

Construction = 446bj2z's on H1, as 537 runs it: `run654.py cpm <design> <seed> --bulk ps --mde-c 0.1 --mde-cap 0.1
--mde-min-count 0.3 --mde-coarse --zs dn_k250_zs`, with the production emitter 0f56209d (asserted). Seeds 0, 1, 2. Read the
zs stats (production).

## Known answer (aborts if it fails)
- `run654.py cpm base 0 …` == 537's `work/h1_cpm/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse{,_dnk250zs}.npz` bit for
  bit (q, p, lfc, sums), pre-zs and zs.
- Seeds 1–2 base: 630's stored `work/h1/stats-base_s{1,2}_…_dnk250zs.npz`, the same construction (the base used by 631's
  decomposition), are re-checked by running `run654.py cpm base 1` and comparing.
- The copies' diffs against run496 / run537 are recorded and must contain only the design branch, the import and the work
  path.

## Measures (paired per seed vs base; 521's four regimes; lib521 members with the realised order, as score625d / score612)
1. The six members.
   - pds and mse are pinned by construction. The check: per-gene block totals (`sums`) are identical across arms;
     otherwise the arm is a bulk leak and is void.
   - fid / jac / reach / nmae raw Δ, and overall-equivalent Δ (lib521.overall_delta: 192's INV) at f 1 / 0.5 / 0.35 / cmp.
2. 631's regions: head (400 by |key|) / CPM block (2,200) / uncalled, each × DELIVERED direction (UP / DOWN).
   - Per region at f 1 and f 0.35: the Δ nmae share, and the mean delivered |lfc| over the region's pool genes ([d4cd5f]'s
     ask). Region masks come from lib521.load_state, as in decomp631.
3. Calls per target (q < 0.05): the call set should be unchanged, since CallSet(400, 2200, 1.9) is realised.

## Predictions (registered now)
- **P1, retrodiction of 410 §7's sign on today's construction:** ×1.25 vs ×1 overall-eq at f 0.35 < 0, with nmae worse (> 0
  raw) and reach worse. If this fails, 410 §7 does not carry to today's construction, and P2's basis is gone.
- **P2, the sign bet:** ×0.8 vs ×1 overall-eq at f 0.35 > 0 (×0.8 BETTER), small: +0.0003 … +0.0015, with nmae −0.0005 …
  −0.0025 raw and reach ≥ 0.
  - Reasoning: 410 §7's nmae curve is convex and rising at ×1 (slope ≈ +0.01 per unit), and the point was never tuned in
    the ballast era.
  - **Against my bet:** 631 says called DOWN genes are under-delivered ~2×, and the coordinator expects a global ×0.8 to hurt
    there.
- **P3 ([d4cd5f], from 631 + 484/490):** the mass-balance kf (~0.60) pins head DOWN genes' delivered mean, because the
  ballast tops each gene up to the bulk. So amplitude barely reaches them. Test:
  - mean delivered |lfc| on head DOWN moves by < 5 % relative between ×0.8 and ×1.25;
  - it moves by > 10 % relative on block UP or on uncalled UP.

## Bins (overall-eq Δ ×0.8 − ×1, needed at BOTH f 0.35 and f 0.5, all three seeds the same sign)
- **≥ +0.0010:** a candidate knob. Propose a builder knob for the carrier amplitude, decoupled from α (default unset =
  today), a dose curve (×0.7 / ×0.9), and 614's K562 real-emitter read before any October use (11.1c.2). Not an October
  switch by itself.
- **(0, +0.0010):** inconclusive at H1's resolution. Recorded; no proposal.
- **≤ 0:** amplitude is at or past its optimum below ×1 too. 410 §7's "DE amplitude is not the lever" extends to ×0.8.
  Closed.
- Whatever the bin, the region table is reported. Where amplitude reaches is the input [d4cd5f] asked for v7.

## Resources
- Each run goes through `experiments/602-council-gate/gate602.sh` (one ticket per run, max ~8 GB), serially.
- 7 runs: base s0 (KA); x080 s0–2; x125 s0–2; base s1 (re-check). Analysis < 3 GB.
- Outside the coordinator's windows (10-04 05:00Z – 10-05 18:00Z, 10-06 06:00–18:00Z, 10-07 06:00–12:00Z).

## Amendment A1 (10-02 18:05Z, commit b1590abd; written AFTER seed 0's as-registered read and BEFORE any zero-repair run)
**What seed 0 showed.** The as-registered harness (537's command) runs WITHOUT production's zero repair:
- production sets `VCC_DE_ZERO_REPAIR=-5`, and run496's `--zero-repair` defaults to 0 = off;
- at ×1.25, 10 gated cells per seed-0 panel reach a zero block total (delivered lfc ≈ −35; base has 1, ×0.8 has 0);
- 3 of them sit in head-DOWN f 0.35 pools and carry most of ×1.25's nmae loss there (+31 error each). Production would floor
  them near log2 −5.
So the as-registered ×1.25 read is not production-faithful. (A collapse that production also repairs is part of that loss,
not the lever's.)

**Added arms (the production-faithful construction):** `--zero-repair -5` on base and x125, seeds 0–2 (6 runs).
- `zero_repair` is an exact no-op on a block with no zero-total gated gene (it returns the input object), so x080_zr ==
  x080 on every seed whose stats have no lfc < −30. That is checked per seed from the stats, and an x080_zr run is added for
  any seed that has one.
- The function body is identical to production's `src/vcc2026/zero_repair.py`; only the repair's own RNG seed formula
  differs.

**Structural KA:** base_zr s<k> equals base s<k> on every target except the one(s) with a zero-total gated gene (q, p, lfc,
sums).

**Reading:**
- The preregistered bins (P2) and P1 are read on BOTH constructions, reported side by side.
- The production-faithful (zr) read governs the conclusion about the production lever. The as-registered read stays
  reported, as registered.
- The bins and predictions are unchanged.

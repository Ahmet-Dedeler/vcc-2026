# 630 — Rank-neutral UP magnitude on the REAL production emitter: H1 nmae −0.0038 raw with fid/jac/reach/pds/mse pinned (v1); v2 running

> **VOID ROWS (2026-10-02T19:10Z; [21dd4c]'s 616 audit, 7dac7a98, `experiments/616-clean-orion-folds/AUDIT-orion-dependence.md`):** the thinned Orion gate tables' DE truth was ~98 % selection artefact (HEK293T: 37,748 calls on the thinned table vs 774 on clean raw cells; 13 / 132 targets scoreable). In this experiment, 613e's HCT116 / HEK293T premise rows (UP-sign share, footprint) are VOID. 630's prices are all H1; the premise reading rests on H1 + K562.

2026-10-02, council session [d4cd5f], council portfolio #1 (the design-sweep pivot). Real production emitter
(rc_emit_k 0f56209d, 446bj2z's DE construction: CLEAN4 q80 block prior, coarse MDE, ballast, zero repair, zs).
Scripts: `run630.py` (dumps the final blocks via run496 with outputs redirected here), `rnup630.py` (the pass + six-member
scoring with evaluate323), `price630.py` (521's regimes, 537's machinery). Every heavy run goes through
`experiments/602-council-gate/gate602.sh`: ~7 GB per dump, ~4 GB per scoring run.

## Why
nmae reads per-gene arithmetic-mean CPM. Calls, fid, jac and reach read Wilcoxon ranks. pds and mse read per-gene block
totals. For uncalled genes we sign UP, three independent analytic pricings on H1 (603, 610, 632) found the L1-optimal
delivered |lfc| above today's MDE add. But raising it through the emitter leaks calls (426/435). DOWN genes HURT under the
same rule (603: +0.024; 632: +0.0035), so UP only.

## Known answers
- run630 seeds 0/1 reproduce 435's six members: nmae equal to 1e-16, the other five exact.
- rnup630 re-applies the production zs pass (src/vcc2026/rank_mag.py) to each dumped pre-zs block. It reproduces run630's
  zs stats (sums, jk, p, q, lfc) bit for bit on all 85 targets, and the six members match its metrics json exactly.
- The strict variant (`up_k100`: recipients must already sit above the control max) is an exact no-op: zero capacity.

## The capacity wall (measured)
- STRICT U-invariance has ZERO capacity. On real blocks, none of the ~3,700–4,000 uncalled UP genes per target has
  any carrier cell above the control max. The top carrier sits a median 116 controls below it (q10 45, q90 245).
- So one carrier per gene is promoted to just above the control max. ΔAUC for that gene ≈ 116/(400·18,400) ≈ 1.6e-5.
- The ballast is the count donor, with matched sub-gate library compensation (carrier gives a sub-gate count → the same
  sub-gate gene on the donor ballast row). Per-gene totals of every gene and every library are exact, so untouched genes
  keep identical values and U.
- The binding limit is the DONOR: the deep ballast cells cannot give up a count of most genes without crossing a
  neighbouring control value (about 2,600 of ~3,900 genes have no in-interval donor room). More recipients do not help
  (R 1 → 20 moves 1,018 → 1,245 genes).
- Ballast holds only ~20 % of each gene's block total (carriers 80 %), so adding counts without a donor would move
  pds/mse visibly. Not an option.

## v1 result (promote-one, ballast donors inside their intervals): six members, H1 (f 1), seeds 0/1

| variant | Δnmae raw | Δfid | Δjac | Δreach | Δpds, Δmse | delivered / asked | calls changed | max ΔAUC | max \|Δlog10 p\| |
|---|---|---|---|---|---|---|---|---|---|
| up_k100_p (K 1.0) | **−0.00383 / −0.00373** | 0 / +2e-6 | 0 / −0 | 0 / 0 | 0 / 0 exactly | 0.112 / 0.113 | 0 / 1 (MED19) | 5.8e-5 | 0.004 |
| up_k125_p (K 1.25) | −0.00237 / −0.00241 | 0 / +2e-6 | 0 | 0 | 0 | 0.073 / 0.074 | 0 / 1 | 5.9e-5 | 0.004 |

At 521's regimes (`price630.py`, 2 seeds, overall-equivalent [target-bootstrap 95 % CI]):

| variant | f 1 | f 0.5 | **f 0.35** | cmp |
|---|---|---|---|---|
| up_k100_p | +0.0010 [−0.0001, +0.0022] | +0.0020 [+0.0008, +0.0031] | **+0.0019 [+0.0007, +0.0031]** (nmae −0.0070) | +0.0021 |
| up_k125_p | +0.0007 | +0.0018 | +0.0018 | +0.0020 |

Only nmae moves. K 1.25 is worse than K 1.0 at this capacity, because large asks drain the recipients' budgets first.

## Variants tried and rejected
- "relaxed" (ballast donors give everything, crossing controls freely): delivers ~17 %, but 0–2 call flips per target
  and |Δlog10 p| up to 0.2. Not rank-neutral. Rejected.
- R > 1 recipients per gene: capacity barely moves (donor-limited); ΔAUC rises to 2–6e-4. Not used.

## v2 RESULT (the coordinator's bar): six members, H1, real emitter

| variant | Δnmae raw f 1 (s0) | fid / jac / reach / pds / mse | genes moved / selected | delivered | call flips | max ΔAUC | max \|Δlog10 p\| |
|---|---|---|---|---|---|---|---|
| s_k100 (in-interval donors) | −0.0059 | all exactly 0 | 963 / 3,693 | 0.143 | 0 | 9.9e-5 | 0.010 |
| **b_k100 (budgeted crossings)** | **−0.0066** | **all exactly 0** | 1,584 / 3,693 | 0.254 | 0 (guard excluded 0.19/target) | 1.0e-4 | 0.013 |
| b_k125 | −0.0029 | all 0 | 1,619 / 3,811 | 0.173 | 0 | 1.0e-4 | 0.014 |
| b_k150 | −0.0018 | all 0 | 1,606 / 3,831 | 0.126 | 0 | 1.0e-4 | 0.014 |

At 521's regimes (`price630.py`, seeds 0/1; overall-equivalent [target-bootstrap 95 % CI]):

| variant | f 1 | f 0.5 | **f 0.35** | cmp |
|---|---|---|---|---|
| **b_k100** | +0.0015 [+0.0001, +0.0030] | +0.0027 | **+0.0028 [+0.0013, +0.0043]** (nmae −0.0103) | +0.0030 |
| s_k100 | +0.0016 | +0.0025 | +0.0025 [+0.0010, +0.0038] | +0.0029 |
| b_k125 | +0.0008 | +0.0023 | +0.0024 | +0.0026 |

- The gain grows in the board's thinner-tail regimes (f 0.35–0.5), as the analytic 603 predicted.
- Larger K is worse at this capacity: big asks drain the recipients' sub-gate budgets. So the board dose pair is
  **K 0.75 / K 1.0** (b_k075 queued on 3 seeds; seed 2 of K 1.0 running).
- **Builder path (UNAPPLIED, for the coordinator):** `src_rank_up.py` (→ `src/vcc2026/rank_up.py`), `patch_builder630.py`
  (VCC_DE_UPRN, after zs; 4 anchors, applies cleanly to a copy of the live builder) and `test_rank_up630.py`: the module
  reproduces rnup630's b_k100 stats bit for bit (**PASS on all 85 seed-0 targets**, 2026-10-02 12:20Z). Prereg: `PREREGISTRATION-DRAFT.md`.

## v2 design (the coordinator's bar, 2026-10-02)
The bar is: 0 call flips (recomputed BH calls; a guard re-runs the block without any gene whose call flipped),
|ΔAUC| ≤ 1e-4 per gene (net), pds/mse ≤ 1e-6, all six members measured, reach reported separately.

`up_pass2`:
- per gene, one promoted carrier;
- a net U budget: the promotion's upward U change is offset by ballast donors crossing controls downward, within
  |ΔAUC| ≤ 1e-4, otherwise the gene is reverted;
- `s_k*`: donors stay inside their intervals. `b_k*`: budgeted crossings.
Smoke test on 3 blocks: b_k100 delivers 21–26 % (s_k100: 10–18 %), with 0 call flips, max ΔAUC exactly 1e-4, max
|Δlog10 p| 0.008, and per-gene totals exact. Seeds 0/1/2 × {s_k100, b_k100, b_k125, b_k150} are queued on the gate.

## Closures it must escape (the coordinator's condition 2)
- **426/435 (MDE doses above c 0.1 leak calls):** the leak is a rank effect. This pass moves no call (guarded) and keeps
  |ΔAUC| ≤ 1e-4, so the leak cap does not bind.
- **426b (prior-sign MDE over-fills the classes):** no class fill. The pass acts on the final block after emission, zero
  repair and zs, and never changes a delivered sign (UP genes only get larger).
- **434/480 (no learned / class-based magnitude rule beats the MDE; the next nmae points are signs):** 480 priced
  both-direction class rules under leak-capped delivery. Here the rule is UP-only (DOWN measured harmful) at K·MDE,
  above the leak cap. It is consistent with 480 that the larger prize needs better signs: the hybrid key/identity sign
  source (612/622) is a separate switch that plugs into the same pass as the UP selector.
- **Up/down asymmetry:** measured. DOWN is not touched.

## Class (11.1c)
Scorer mechanics: delivery only, computed from each context's own controls (control max, intervals, ref_lin) and the
block itself; nothing is fitted to A/B/C. K is one global dose. Its signs are today's delivered signs, so it does not make
the model better at biology.

## Builder readiness (2026-10-02 ~12:25Z; UNAPPLIED, for the coordinator)
- `patch_builder630.py --root <repo> [--check]`: order 538 → 554 → 590 → 598 → 612 → **630**. It refuses without 612's
  mark and refuses a second apply. 4 anchored insertions in `scripts/make_submission_pm_cd4.py` (knob, per-context setup,
  per-block call after zs, meta `de_uprn`); it installs `src/vcc2026/rank_up.py` = `src_rank_up.py`, sha pinned
  (5ae935b9…).
- **`test630.py`: 27/27 PASS** on temp trees built with 598's own canonical post-bundle construction (test598.canon_tree).
  Refused out of order; 598 → 612 → 630 applies, compiles and equals `patched_text`; the diff is only insertions (44
  lines = the patch's); the knob grammar is validated (K in [0.5, 3]; needs VCC_DE_MAG); the sha chain is consistent.
- **KA plan:**
  - KA0: knob unset on the post-630 tree rebuilds the base byte for byte (the coordinator's ec3a5ae6 row; no code runs
    when unset).
  - KA1: `test_rank_up630.py <seed> b_k100`: the installed module reproduces the H1 pricing run bit for bit (**PASS on all
    85 targets**, 12:20Z; `work/ka1-s0.log`).
  - KA2: one K 1.0 artifact whose config diff vs the base is exactly `de_uprn` (rule 10, def_postbuild
    `--allow-config-diff de_uprn`), with pds/mse equal to the base to ≤ 1e-6 raw in its local check.

## Informational: relaxing the U bar does not help (m_k100, seed 0, 12:22Z)
Per-gene sign-aware z-margin bounds ([21dd4c]'s idea): each gene may move its Wilcoxon z anywhere keeping |z| ≥ 0.5
below its target's BH critical z, capped at |ΔAUC| ≤ 1e-3, i.e. 10× the shipped bar. Delivered fraction 0.313 (vs
b_k100's 0.254), but FEWER genes moved (1,323 vs 1,584: the first genes take more budget). Δnmae −0.0064 (b_k100
−0.0066), f 0.35 +0.0024 (b_k100 s0 +0.0027), max |Δlog10 p| 0.10, 0 call flips. **The coordinator's 1e-4 bar costs
nothing; b_k100 under 1e-4 stands.** Side fact: on uncalled genes the emitter's Wilcoxon z and the delivered lfc agree in
sign only 50 % of the time (z ≈ −1 on average even where lfc > 0), so the promotion moves z toward 0 for most UP genes.

## October path (TODO after the board read; deliberately deferred)
A `DEF_DAY_DE_UPRN` pass-through is needed only if the dose pair reads TRANSFERS or PARTIAL (10-08 at the earliest). Pattern:
502's `DEF_DAY_DE_MAG`, which touches `scripts/def_day.sh` (validate the variant; check the builder has the knob; export
`VCC_DE_UPRN` for the level builds; write it to def_day.final.env; add `de_uprn` to the config-diff ALLOW list),
`scripts/def_day_helpers.py riders_match` (the final meta must carry the rider), `scripts/def_pipeline.py` (a config row
`de_uprn.variant`) and `scripts/october_recipe.py` (switch → env, BUILDER_NEEDS, export line). It must compose after
598/612/630 on the canonical post-bundle tree (test598.canon_tree), with an inert KA (unset == today's def_day output).
The pass needs nothing but each context's controls and the block, so D/E/F's day-of builds compute it unchanged.

**630dd (UNAPPLIED, built by [4d0207], commit 626949ef):** `patch630dd.py` (sha 4db2a0c5…) covers the def_day.sh part:
six anchors on DE_MAG's original lines; V refuses without DEF_DAY_DE_MAG, outside the builder grammar, or without the
630 mark; C requires meta `de_uprn.variant` == request and `module_sha256` == d8c39c3f (regenerate if 607 changes the
module); unset knob is inert. test630dd 39/39. 606's chain 641 → 598 → 612 → 612dd → 630 → 630dd → 642 → 604 → 644
ALL PASS with 612dd / 630dd in either order (def_day 7f78bd04). Still open, only if the board read adopts it:
def_pipeline.py's config row and october_recipe.py's switch.
**630dd v2** (28bc1e46, patch630dd ce9f0fe0): module constant 7e6b3e52; the postbuild step REPORTS per context (targets,
mean ur_frac_got / ur_n_moved, guard_failed / invariant_failed counts, no_pos_donor genes, the prereg capacity floor)
to the log and `<work>/uprn_audit.json`, never fatal. test630dd 43/43. 606 re-pinned: def_day f9552b6b, builder
2c6cfe7d, rank_up 7e6b3e52, ALL PASS.

## A/B/C capacity pre-check, v2 b_k100 (12:32Z): B FAILS the n_moved gate
`capacity630.py` + `work/capacity_chain.sh` restore 446bj2z's package (it carries zs) and run the production module on
60 real targets per context, with each context's own scorer-order controls. The restored h5ad is deleted after
(`work/capacity630_b_k100.json`).

| context | delivered (gate ≥ 0.127) | genes moved / selected (gate ≥ 500) | ballast K | binding constraint |
|---|---|---|---|---|
| A (Jurkat) | 0.170 ✓ | 772 / 2,996 ✓ | ~9 | "no candidate" 2,223 |
| B (HeLa) | 0.134 ✓ | **400 / 2,897 ✗** | ~9 | "no candidate" 2,497 |
| C (CAL33) | 0.213 ✓ | 1,017 / 3,317 ✓ | ~9 | "no candidate" 2,299 |
| (H1, reference) | 0.254 | 1,584 / 3,693 | 3 | |

A/B/C deliver ~69 % of H1's fraction, so the board prediction scales down accordingly. **Diagnosis** (replayed on H1
blocks): on the UNMODIFIED block almost every selected gene is servable (promotion needs a median 4 counts; donors have
11–64; the top carrier holds ~100 sub-gate counts). "No candidate" arises because the first, largest asks drain the few
shared top carriers' budgets.

## v3: load-spreading (12:35Z)
`up_pass2(..., quantum, rounds, rec_pool)`. Each gene receives at most `quantum` counts per round, over `rounds` passes
(water-filling: the lfc gained, log2(1 + k·cw/A), is concave in counts, so the first counts per gene are worth most).
A promoted gene's recipient is chosen among the `rec_pool` cheapest-to-promote carriers by most remaining budget. Defaults
reproduce v2 bit for bit (re-verified after removing an early exit added at 12:20Z that had changed v2's semantics; every
priced v2 / K 0.75 result and KA1 predate it). Smoke on 3 H1 blocks, b_k100, vs v2:

| setting | genes moved | delivered | Σ Δlfc on selected genes |
|---|---|---|---|
| v2 | 1,390–1,613 | 0.21–0.26 | 262–353 |
| q20 r3 p10 | 2,648–2,919 | 0.33–0.38 | 331–376 |
| q40 r3 p10 | 2,256–2,632 | 0.34–0.39 | 335–402 |

0 call flips, ΔAUC ≤ 1e-4 throughout. Full 3-seed six-member scoring of q20r3p10 and q40r3p10 is queued, and the A/B/C
capacity check will re-run with the winner. The module (`src_rank_up.py`, sha 8e1f05a6) carries v3 (AST-equal to the
driver); the knob grammar accepts `b_k100_q20r3p10`; test630 31/31 PASS.

## v3 on 3 seeds, and the dose pair (13:18Z; `price630.py`, `work/price630-dosepair.log`, `work/h1/price630.json`)
Overall-equivalent at 521's regimes, seeds 0/1/2 [target-bootstrap 95 % CI]. fid / jac / reach are exactly 0 for every
variant (0 call flips); pds/mse pinned check passes (per-gene sums equal).

| variant | f 1 | f 0.5 | **f 0.35** | cmp | nmae raw f 0.35 |
|---|---|---|---|---|---|
| b_k100 (v2) | +0.0016 [+0.0001, +0.0031] | +0.0028 | +0.0030 [+0.0014, +0.0044] | +0.0030 | −0.0107 |
| b_k075 (v2) | +0.0019 | +0.0026 | +0.0025 | +0.0028 | −0.0090 |
| **b_k100_q20r3p10 (v3)** | +0.0023 [+0.0005, +0.0042] | +0.0030 | **+0.0029 [+0.0015, +0.0044]** | +0.0034 | −0.0106 |
| **b_k075_q20r3p10 (v3)** | +0.0023 [+0.0006, +0.0041] | +0.0031 | **+0.0030 [+0.0015, +0.0044]** | +0.0034 | −0.0108 |
| b_k100_q40r3p10 (v3) | +0.0015 | +0.0029 | +0.0027 | +0.0029 | −0.0099 |

- v3 (q20 r3 p10) gains over v2 at f 1 (+0.0007) and cmp (+0.0004), ties at f 0.35. q40 is worse than q20.
- **K 0.75 and K 1.0 are indistinguishable under v3** (≤ 0.0001 apart at every regime). A board pair of them would
  read the same number twice. The dose pair needs a dose H1 separates: K 0.5 / K 1.5 under v3 are queued.
- A/B/C capacity under v3 b_k100_q20r3p10 (60 targets/context, 446bj2z's real blocks): A 0.284 / 1,550 moved,
  B 0.295 / 970, C 0.325 / 1,985 (frac_got mean / genes moved per target). **All three pass the gate** (≥ 0.127,
  ≥ 500). v2 had failed B (400). Build cost: ~2.8 s per target (157–185 s per 60), i.e. ~40 min added to a
  900-target build.
- KA1 for v3 (module == pricing run, seed 0, 85 targets): PASS (`work/ka1v3-s0.log`).

## 605 cloud-review fix: negative-gain donors (13:15Z; `test_neggain630.py` 11/11)
**Finding** (Grok via [4d0207], 605/REVIEW.md). One count of gene j moved from ballast row b to carrier a changes the
mean CPM by cw[a] − cw[b], with cw = 1/(n·library/1e6). That is positive only if the donor's library is LARGER than the
carrier's. The pre-fix ask clamped its denominator at 1e-300 and `ZS.fill` took from every donor. With shallow donors,
the mean could fall and an UP lfc flip DOWN while libraries, per-gene totals and |ΔAUC| stayed green.
Synthetic case (8 carriers at 22k, 32 ballast at 18k): +0.00144 → −0.286 with 12,736 counts moved.

**Exposure on production blocks: none MEASURED (an empirical claim, not structural; [4d0207]/607).** The filter is
correct for any block. It is idle only while every ballast library exceeds every carrier library, and that is measured, not
guaranteed: on H1 85/85 here, and on A/B/C by capacity630's `neg_donor_rows` counter in the chain605fix rerun.
Ballast rows are ~995k counts and carriers ~20k (H1 seed-0 blocks: median
ballast library 995,000, minimum over targets 994,978; carrier maximum 24,917). On 85/85 targets no ballast row's
library is at or below any carrier's, so every move had positive gain. The fixed module's KA (bit-identical to every
priced H1 run, seeds 0–2 × b_k100 / b_k075_q20r3p10 / b_k100_q20r3p10) and the A/B/C capacity rerun with the new
counters are running (`work/chain605fix.sh`).

**Fix** (module and driver, AST-equal):
- donors restricted to `cw[bal] < cw[a]`; a gene with none is skipped (`n_no_pos_donor`);
- the ask divides by the eligible donors' largest gain. When any donor is eligible this equals the old denominator,
  so the fix is a no-op on all-deep blocks (test N3: bit for bit);
- the gene is reverted if its mean falls or its lfc sign changes (`n_mean_revert`);
- new counters `ur_n_neg_donor_rows`, `ur_n_no_pos_donor`, `ur_n_mean_revert`.

**Guard fails safe.** If calls still flip after 5 guard rounds, that target keeps its zs block and the meta records
`ur_guard_failed = 1`; nothing raises. A raise on D/E/F day would abort the build.

Module sha after this fix: d8c39c3f (was 8e1f05a6; pre-fix copy kept at `work/src_rank_up_prefix605.py`).

**Invariant asserts fail safe too** ([7772dc]'s review of 630, `600/review630_7772dc.py`, feb07f20: 1,054 synthetic
runs, 0 fails, all five claims hold; the fix is a no-op on deep-ballast blocks and needed on hostile ones, where v2/v3 also
raised at ballast_k = 0). `apply()` wraps the pass. Any AssertionError (libraries, per-gene totals, negative counts,
the U record, rank_mag's helpers) returns the input block with `ur_invariant_failed = 1`. It is bit-identical when
nothing fires, and all three paths return the same 19 keys (`DIAG_KEYS`). test_neggain630 14/14 (N5-N7 added),
test630 31/31.
**FINAL module sha 7e6b3e52e2ad374b7b6e21e92327fcf65771433693eb3846c09c0189f3f2959f.**
The bit-for-bit KA of the fixed module against all nine priced H1 runs was re-started on this module (gate ticket 52,
`work/chain605fix.log`). The first two KAs on d8c39c3f had passed (`work/chain605fix.d8c39c3f-partial.log`).

**mse is pinned by a cap, not by construction** ([7772dc] (b), [4d0207]). The pass raises the prediction's
jackknife term (H1 seed 0: jk 37 → 51 per target, median +36 %). mse stays bit-identical only because cell-eval2's
per-row cap min(C_pred, k·C_real) (`metrics/delta.py` mse_unbiased_capped, #247) already binds on ballast blocks, and
the pass only raises C_pred. KA2's |Δmse| ≤ 1e-6 on A/B/C is the gate (the prereg's falsifier is at 1e-4 raw).

## Capacity ceiling: the oracle (13:30Z; `oracle630.py`, `work/oracle630.log`)
Each selected gene's delivered lfc set to its full target t = max(l, K·MDE), with the same selection as the module and
everything else from the zs base. This is the upper bound for any rank-neutral UP pass: only nmae can move. Overall-equivalent,
3 seeds [95 % CI]:

| | f 1 | f 0.5 | **f 0.35** | cmp |
|---|---|---|---|---|
| oracle K 0.5 | +0.0016 | +0.0021 | +0.0019 | +0.0022 |
| oracle K 0.75 | +0.0030 | +0.0040 | +0.0039 | +0.0046 |
| **oracle K 1.0** | +0.0029 [−0.0003, +0.0062] | +0.0060 | **+0.0060 [+0.0034, +0.0087]** | +0.0070 |
| oracle K 1.5 | **−0.0101** (overshoot) | +0.0044 | +0.0077 | +0.0048 |
| actual v3 K 1.0 (delivers ~⅓) | +0.0023 | +0.0030 | +0.0029 | +0.0034 |

- **Full delivery at K 1.0 doubles the gain at board-like regimes** (+0.0060 vs +0.0029 at f 0.35; +0.0070 vs +0.0034
  cmp). It matches the analytic 603 (+0.0060). The capacity headroom is ≈ +0.003 overall-equivalent on H1, ×1.2–1.6 on
  the board (611's delivery multipliers). Engineering more capacity is worth it.
- K 1.0 is the robust target. K 1.5 overshoots at f 1 (−0.010), and K 0.75 leaves a third of the gain.
- At K 1.0, today's pass sits between oracle K 0.5 and oracle K 0.75. `diagcap630.py` tags which branch binds.

## October premise reading (11.1c.2): signs transfer; the footprint question goes to the board (13:32Z)
[21dd4c]'s 613e (`613/block613e.py`, committed before running; `613/work/block613e.out`) measured uncalled gate genes with
UP delivered signs on 636's surfaces at 521's pools:
- **Signs transfer.** The real-UP share is HCT116 0.70 / 0.70 (f 1 / f 0.35), HEK293T 0.73 / 0.79, K562 0.73 / 0.89,
  against H1's 0.54 / 0.58. Headroom: median |T|/add 4.7–8.5.
- **Footprint is small on the folds.** Uncalled-UP gate cells per target at f 1: HCT116 8, HEK293T 38, K562 1, H1 189.
  The folds keep 20–42 % of gate mass uncalled; H1 keeps 55–69 %.

**CORRECTED (13:34Z, after [21dd4c]'s reply).** 613's cell counts required add > 0, and the Orion folds' low depth
(libraries 12–14k vs ~20k on H1 and A/B/C) zeroes the add on many genes. [21dd4c] retracted "5–25× smaller". In nmae's
own units (613e2), the share of each target's Σ|T| on uncalled UP-signed gate genes:

| surface | f 1 | f 0.35 |
|---|---|---|
| H1 | 0.25 | 0.23 |
| HCT116 | 0.185 | 0.10 |
| HEK293T | 0.184 | 0.11 |
| K562 | 0.134 | 0.08 |

So the footprint is 1.3–1.9× smaller than H1's at f 1 and ~2–3× smaller at f 0.35. Part of that is the Orion folds'
consensus key: each Orion line's sister (HCT116 ↔ HEK293T) is a same-platform donor, so the head captures more of the
gate.

Two of my arguments were wrong and are withdrawn:
- (1) 433a's coarse MDE add has `mde_scope` "all" (block genes too), so its 1.59× board ratio does not isolate the
  uncalled footprint.
- (2) The K562 fold's key donors are the two Orion lines; there is no K562 donor. Same-line donors exist in CLEAN4's
  banks, not in the key.

**Orion caveat (14:19Z; [4d0207]'s 609b, 7b419f72):** on the Orion folds the LOW-expression uncalled genes skew UP (truth
0.83 HCT116 / 0.68 HEK293T). An expression-only "bottom decile → UP" rule earns +0.26 / +0.15 there; this is the mirror of
the high-expression DOWN skew (CPM normalisation of shallower perturbed cells inflates low counts). So the HCT116 /
HEK293T rows above (UP-sign share 0.70–0.79, footprint 0.10–0.19) are artifact-inflated and are NOT readings. **Valid
surfaces for UP-pass and UP-sign levers: H1, Broad-16 and K562 only** (K562: real-UP share 0.73 / 0.89, footprint
0.13 / 0.08 of |T| mass). On clean H1 the expression-only control is flat (0.580 = CLEAN4). Any Orion row shown later
sits next to an expression-only row.

**Expectation for A/B/C: between H1-like and ~2× smaller** than the H1 f 0.35 price × 611's delivery multiplier. The 630
board read settles it. If it reads below 0.3× the H1 price, the footprint is smaller still, and October's expectation
for 630 shrinks accordingly.

## v3 dose bracket and the capacity diagnosis (13:43Z; `work/price630-dose.log`, `diagcap630.py` → `work/diagcap630.log`)
**Dose (v3, 3 seeds, overall-equivalent):**

| K | f 1 | f 0.5 | f 0.35 | cmp |
|---|---|---|---|---|
| 0.5 | +0.0017 | +0.0022 | +0.0020 | +0.0023 |
| 0.75 | +0.0023 | +0.0031 | +0.0030 | +0.0034 |
| **1.0** | +0.0023 | +0.0030 | **+0.0029** | +0.0034 |
| 1.5 | +0.0020 | +0.0027 | +0.0025 | +0.0031 |

The curve is flat from K 0.75 to 1.0 and falls off on both sides by ≤ 0.001. Above K 0.75 the pass is capacity-
saturated: a larger K asks for more than it can deliver. **Decision (the coordinator's rule): ONE board arm,
`b_k100_q20r3p10`.** The freed slot goes to 524dr2. A dose pair becomes worth two slots only once capacity grows.

**Why delivery stops at ~⅓ (12 seed-0 H1 targets):** the binding constraint is the library-compensation currency.
- Carriers hold only ~44k sub-gate counts per block, 0.55 % of their library.
- The pass spends 93 % of them (41k), and 76 % of carriers end fully drained.
- Per target: k = 0 events (drained recipient) ≈ 5,050; recipient budget below the promotion cost ≈ 950; U-budget
  reverts ≈ 2,100; donor shortfall ≈ 26; no positive-gain donor 0; cmax reverts 0.
- So flags a / z / r (the v4 screen) can only re-allocate an almost exhausted pool. Real capacity needs a new
  compensation currency.

Candidate: counts of DOWN genes parked on zs's mass-holder carriers above the control max, paid carrier → ballast. That
keeps the library and per-gene totals exact and lowers the DOWN gene's mean, the direction zs is capacity-bound on.
Its size is measured next (`currency630.py`).

## Knob grammar tightened (607's MINOR; builder patch only, module 7e6b3e52 unchanged)
`VCC_DE_UPRN` must match `(s|b)_k(\d{2,3})(?:_q([1-9]\d*)r([1-9]\d*)p([1-9]\d*))?` with K in [0.5, 3]. A zero or
zero-led quantum / rounds / pool is refused: `b_k100_q20r0p10` used to be accepted and moved nothing. test630 35/35.
The module's own parse_spec stays permissive (production reaches it only through the knob). It gets the same check in the
next module revision.

## Design decision: the fail-safe catches every AssertionError (607b MINOR, kept on purpose)
607b (Grok via [4d0207]) confirmed 7e6b3e52 is bit-identical to d8c39c3f when nothing fires. It noted that `apply()`
also catches the precondition ("base lfc does not match the block") and any future programming assert, so a systematic
bug would silently skip the pass on every target. The catch-all is kept, for three reasons:
- every caught path returns the input block unchanged, so the pass can be lost but a wrong block is never written;
- a systematic failure fires on every target, and the coordinator's rule (any FAILSAFE trigger in the KA630 build = KA
  FAIL) catches it on A/B/C first;
- on D/E/F day a re-raise would abort the whole build and force the same outcome (no pass) later.
630dd is asked to make it loud: "630 EFFECTIVELY OFF in context <c>" when fail-safe triggers exceed 10 % of a
context's targets. Candidates for the next module revision, batched into one sha change: an allow-list re-raising
programming asserts, parse_spec's q/r/p ≥ 1 check, and trying the next recipient candidate.
606 pin with the tightened knob: builder fe4678c3, def_day 7075d83d (630dd mirrors the grammar, test630dd 50/50).

## KA verdict on the FINAL module 7e6b3e52 (14:26Z; `work/chain605fix.log`, `work/capacity_chain_v3fix.log`)
- **9/9 bit-for-bit PASS.** `src_rank_up.apply` (module 7e6b3e52) reproduces every priced H1 run (seeds 0/1/2 ×
  b_k100_q20r3p10, b_k075_q20r3p10, b_k100; 85 targets each). So the pre-fix code never used a negative-gain donor, never
  lowered a mean and never hit a fail-safe on H1. Every priced number above stands unchanged.
- **A/B/C (446bj2z's real blocks, 60 targets/context):** neg-donor rows 0, no-pos-donor 0, mean reverts 0, guard
  failures 0 in all three contexts. Per-target frac_got and n_moved are identical to the pre-fix run on 60/60 targets
  per context. This answers the coordinator's exposure KA on A/B/C: none, measured.
- Capacity unchanged: A 0.284 / 1,550, B 0.295 / 970, C 0.325 / 1,985 (gate ≥ 0.127 / ≥ 500 passes everywhere).
- Peak RSS of the KA chain 6.05 GB (one gate602 ticket, v2 rules).

## v4 (allocation flags) NULL, and the DOWN-holder currency is ballast-poor (14:44Z)
**v4** (`up4_630.py`, `work/price630-v4.log`, 3 seeds; the flags-off KA == v3 bit for bit on every seed). Overall-
equivalent at f 0.35 / cmp:
- v3 +0.0029 / +0.0034;
- a (above carrier must have budget) identical;
- z (zero-cost donor room after the budget) +0.0029 / +0.0034;
- r (re-select a drained recipient) +0.0029 / +0.0034;
- azr +0.0030 / +0.0034.

Every change is ≤ 0.0001 at every regime. Re-allocating an almost exhausted currency cannot raise capacity, as the
diagnosis predicted. **Closed.**

**DOWN-holder currency** (`currency630.py`, `work/currency630.log`, 12 seed-0 H1 targets): carriers hold ~39k counts of
DOWN genes above those genes' control max, about the size of the sub-gate pool. But the ballast rows can receive almost
none at zero rank cost: paired zero-cost currency is ~670 counts per block, 0.02× the sub-gate pool. **The binding side
is the ballast's zero-cost room, not the carriers.** v5 (`up5_630.py`) is correct but can add only ~2 % capacity. Its
H1 screen runs as the negative control. `currency630b.py` repeats the measurement over ALL kept genes, the last check
before concluding that untested (sub-gate) genes are the only free currency.
**currency630b (all kept genes, not only DOWN holders; `work/currency630b.log`):** the carrier side is large (~92k counts
of zero-cost downward room per block), but **no ballast row sits above any gene's control max** (0 genes), and the
ballast's in-interval upward room pairs only ~1,056 counts per block (0.02× the sub-gate pool, all on called genes).
**Untested (sub-gate) genes are the only zero-cost currency.** The remaining capacity route spends rank budget instead
of counts: v6 (library drift, `up6_630.py`), probed on H1 now.

## 630 stacked on the hooks that ship (15:44Z; seed 0; [21dd4c]'s KA'd dumps h1_c4kgd / h1_idkg3d)
v3 `b_k100_q20r3p10` on each construction's own zs base. The zs KA reproduces each run's metrics, 0 calls change, and
fid / jac / reach are exactly 0. Overall-equivalent:

| base hook | f 1 | f 0.5 | f 0.35 | cmp | frac_got |
|---|---|---|---|---|---|
| plain CLEAN4 hook (the price above, seed 0) | +0.0022 | +0.0029 | +0.0027 | +0.0034 | 0.33 |
| **c4kg** (KEYGATE 0.5 + CLEAN4; applies 10-06) | +0.0024 | +0.0032 | +0.0030 | +0.0037 | 0.37 |
| **idkg3** (v3 identity, key-gated) | +0.0028 | +0.0030 | +0.0029 | +0.0038 | 0.33 |

- 630 is ADDITIVE to the shipping hook and compounds mildly with better uncalled-gene signs (f 1 +0.0024 to +0.0028;
  cmp up to +0.0038). This agrees with the earlier idkgd stack and with 632's prediction that better UP signs make
  the extra magnitude worth more.
- **idall3** (the ungated v3 identity hook that replaced KGID, [21dd4c]; 3 seeds, all three dumps KA'd byte-identical):
  f 1 +0.0027 / f 0.5 +0.0029 / f 0.35 +0.0028 / cmp +0.0035, against the plain hook's 3-seed +0.0023 / +0.0030 / +0.0029 /
  +0.0034. **630 is additive on the likely October hook too**, slightly stronger at f 1. 0 calls changed.

## v5 (DOWN-holder currency) NEGATIVE, as the currency measurement predicted (15:44Z; `work/price630-v5.log`, 3 seeds)
Overall-equivalent f 1 / f 0.35 / cmp:
- v3 +0.0023 / +0.0029 / +0.0034;
- d15 (lfc floor ×1.5) +0.0019 / +0.0028 / +0.0034;
- d25 (×2.5) +0.0018 / +0.0027 / +0.0034;
- d15c / d25c identical to d15 / d25: every paying gene was called.

The payer found only ~210–310 DOWN counts per target, on ~19 genes, against ~32k sub-gate counts paid. Its availability
estimate (an upper bound) over-promised, so 125–160 genes per target hit a short payment and were reverted: lost
capacity. The extra DOWN magnitude on those called genes costs nmae slightly at f 1. 0 calls changed. **Closed.**
Capacity routes tried: v4 (allocation) NULL; v5 (DOWN currency) negative; currency630b says no other zero-cost
currency exists. The remaining route is v6 (library drift), running.

## v6 probe (library drift as compensation): not rank-neutral, small gain. CAPACITY LINE CLOSED (16:00Z)
`up6_630.py`, `work/price630-v6.log`, 3 seeds. With unlimited currency, recipients' libraries drift by up to 2.8 %.
- Delivery rises from 0.33 to 0.55 of the ask.
- **9,090–9,110 genes per target exceed the 1e-4 |ΔAUC| bar** (worst 0.0115), and ~20 calls per target flip despite
  the guard. The guard can only exclude UP-selected genes; the flips are compressed bystanders.
- Overall-equivalent f 1 / f 0.5 / f 0.35 / cmp: v6 +0.0022 / +0.0034 / +0.0036 / +0.0039, against v3 +0.0023 /
  +0.0030 / +0.0029 / +0.0034.
- fid / jac / reach move (reach −0.0003 … −0.0015 raw), and 611 says H1 cannot price reach for call-set changes.
So library drift buys ≤ +0.0007 on H1 by breaking the contract that makes 630's board read exact and its transfer
argument mechanical. Not a candidate. A capped drift would buy proportionally less.

**Capacity conclusion.** Every route has been tried:
- allocation (v4): NULL;
- DOWN-holder currency (v5): negative;
- any other zero-cost currency (currency630b): none exists;
- library drift (v6): breaks rank-neutrality.

The UP pass is bounded by the carriers' untested-gene pool (0.55 % of their library), at about a third of its ask.
**630 ships as v3 b_k100_q20r3p10 (module 7e6b3e52).** The only remaining lever is upstream: an emitter that writes
carriers with a larger untested-gene share. That is an emission redesign and is not pursued now.
**Closed on paper too: a reservoir in the target gene's own column.** The own target is excluded from every metric's gene
set, but every bulk profile normalises by the sum over ALL genes, target included (BRIEF: b_{p,g} = log(1 + TS·P_g/Σ_g'
P_g')). A target-gene reservoir of share r therefore shifts every other gene's predicted bulk value by ≈ −r, and the bulk
designer cannot pre-compensate, because the non-target shares must sum to 1 − r. At r = 0.5 % (≈ 70k counts, 1.6× today's
pool) mse pays ≈ r² = 2.5e-5 per gene, ≈ 2.5 % of the panel's mean effect size: about −0.004 overall, more than the UP gain.
At r = 0.1 % the cost and the gain are both negligible. Not pursued.

**KA630 audit tool** ([4d0207], 26c8eddb): `ka630_audit.py <build>.rc_audit.jsonl` runs 630dd's own report snippet
(extracted verbatim from patch630dd at run time) and exits 0 only on `FAILSAFE_TOTAL 0`; test_ka630_audit 6/6. The KA630
build needs `VCC_RC_AUDIT_DUMP=<path>` (remove a stale file first) plus `VCC_DE_UPRN`. **Capacity floor, now gated**
(d4cd5f, 2026-10-02 17:03Z): exit 0 also requires every context in `--contexts` (default A,B,C) to be present with mean
frac_got ≥ 0.127 and mean moved ≥ 500 per target (`FLOOR FAIL <ctx>` otherwise). A context below the floor means the build
is re-diagnosed, not submitted. Reference, 446bj2z v3 blocks (60 targets/context): A 0.284 / 1,550, B 0.295 / 970,
C 0.325 / 1,985, about 2× margin. test_ka630_audit 11/11.

## Depth dependence of delivery levers ([21dd4c]'s H1@11k control, 614 `work/h1t614.json`; information for D/E/F)
h1-test binomially thinned to 11.3k median library (K562's), 3 seeds, f 1, raw Δnmae:
- zs: −0.00337 at 11k vs −0.00183 at 20k. Magnitude reweighting GROWS at lower depth.
- coarse MDE: −0.00648 at 11k vs −0.00966 at 20k. Add-based delivery shrinks to ~⅔.

The depth factor is lever-specific. 630 / 631 have not been priced at 11k. Arc's 2026 post says the released files are
downsampled to a ~20k median (FINDINGS §0), so D/E/F should sit at A/B/C's depth and no rescale is planned. If D/E/F's
control median library comes out far from 20k on Oct 22, re-price 630 / DN on a thinned H1 before the October freeze.

## Next-revision candidate found by 631's KA: float-noise asks
For a selected gene already at or above its target (t = l), A_t − A should be 0 but can be a tiny positive number from
rounding, and the pass then promotes a carrier for it (in 631's K = 0 arm, ~425 genes and ~1,850 counts per target). In
v3 the same can happen for UP genes with l ≥ K·MDE. That behaviour is inside every priced and KA'd 630 number, so nothing
changes for 10-06. The next module revision should require dA > 1e-9·A (batched with the other candidates: the
programming-assert allow-list, parse_spec's q/r/p ≥ 1, next-recipient retry).
**Rule (coordinator):** no 630 module revision before 630's board read. The read must be of exactly what was priced and
KA'd. Reconsider the dA guard (and the other candidates) after the read, and only if its re-price moves a number by more
than noise.

## 10-07 build driver: `build_630arm.sh` (UNARMED; the coordinator launches it after 10-06 APPLY DAY step 8)
446bj2z's recipe (as build549.sh: LODO bulk load 00ccfcba, MDE 0.1/0.1/0.3 coarse, zero repair −5, re-level scales
A 1.0962 / B 1.1461 / C 1.1681, zs dn_k250_zs, q 1.75) + `VCC_DE_UPRN=b_k100_q20r3p10` + `VCC_RC_AUDIT_DUMP`, through
run_arm_generic.sh.
- **Pre-build gates:**
  - builder and emitter == HEAD;
  - the builder carries 630's mark; rank_up.py == 7e6b3e52;
  - no other build or vcc prep; disk ≥ 35 GiB; free+inactive ≥ 23 GiB;
  - CAP 10-07 18:00Z (env CAP).
- **Post-build gates:**
  - config diff vs 446bj2z exactly `de_uprn` (+ `de_sign_emitter_sha256`, since the emitter moved after 446bj2z);
  - CLIFF OK;
  - the zero-gene audit's "0 in 0 blocks";
  - novelty rc 0;
  - a package;
  - `ka630_audit.py` rc 0 (FAILSAFE_TOTAL 0 plus the capacity floor in A / B / C).
- Terminal line: `READY 630arm <pkg sha256>` or `NOT READY <reason>`. Writes `work/build630arm/env.json`.
- DRY=1 prints the plan and gates (today: NOT READY, no 630 mark yet, as expected). `bash -n` clean; BSD grep
  (`/usr/bin/grep -E`) throughout.

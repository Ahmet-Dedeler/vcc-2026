# 631 — PREREGISTRATION (FROZEN 2026-10-02T18:50Z): the DOWN magnitude rider `VCC_DE_DNMAG=d_m200`, read against 630's entry

Drafted by [d4cd5f] before any entry with this knob exists (draft 70eb605a), and frozen 2026-10-02T18:50Z with the coordinator's base and
day (NEXT-SLOT-DECISION). Bins, falsifiers and gates were fixed at the draft and are unchanged.
- **Base:** 630's scored entry: b_k100_q20r3p10 on the best base at 630's arming.
- **Day:** the first free slot after 630's and idall3's reads, before 625.
- **Apply:** patch_builder631 is on the 10-08+ apply list, with a 606-style composition re-run on the applied 10-06 tree
  plus KA0.

**Expected bin: PARTIAL, not TRANSFERS (the coordinator's note, written down before the read).** H1's f 0.35 Δnmae is
−0.0030 raw. A/B/C's capacity is expected at about ½ H1's (630's A/B/C capacity was 0.70× H1's for zs; DN's A/B/C numbers
are below once measured). With 611's ×1.2–1.6 multiplier the centre is ≈ −0.002 raw, in PARTIAL. A PARTIAL read is the
expected outcome, not a disappointment. TRANSFERS needs the footprint and multiplier at the top of their ranges.

**654 interaction:** if 654's ×1.25 carrier amplitude is adopted before this read, it changes the base. DN's ×2 target is
relative to the delivered lfc, so on a ×1.25 base DN aims higher and may overshoot. The stacking test (DN on 654's ×1.25
zero-repair dumps) is recorded below when it lands. If DN adds less than half its standalone H1 gain on that base, the
predicted range is halved before the read. The bins stay; the expectation moves to NULL/PARTIAL.

## Arm (ONE variable vs the base)
- **Base: 630's own scored entry** (VCC_DE_UPRN=b_k100_q20r3p10 on the zs-carrying production base).
- **DN = base + `VCC_DE_DNMAG=d_m200`.** Module `src/vcc2026/dn_mag.py` = experiments/631-nmae-headroom/src_dn_mag631.py
  (sha recorded at apply). Builder patch `patch_builder631.py` (UNAPPLIED; order … → 630 → 631; test631 31/31).
- Expected config diff vs the base: exactly `de_dnmag` (variant + module sha).
- **(c) Unset is byte-identical to the 10-06 builder:** KA0 (unset rebuilds the base byte for byte) before any DN build.

## Mechanism and class
Called DOWN genes are delivered at about half the right magnitude: the mass-balance shrink and mass holders keep typical
DOWN genes under-delivered while the strongest sit near the floor (631, corrected by 654). The pass lowers each called
DOWN gene's arithmetic mean at FIXED per-gene totals: counts move top carrier → the ballast row with the fewest
crossings, net |ΔAUC| ≤ 1e-4 per gene, and the carrier's library is paid back with untested-gene counts from the ballast,
so libraries and totals are exact. **Class: scorer mechanics (delivery).** The dose m = 2 is fitted on H1 (proxy631's
cross-fit). There is no validated held-out reading: 631f's fold model failed its retrodiction, and 614's K562 harness
failed its trust gate. Under 11.1c.2 an October switch therefore needs a written reason. Candidate reason: the mechanism
carries, because the scorer and the emitter are the same on D/E/F, and cancer folds descriptively want MORE DOWN
magnitude than H1, not less.

## Known answers and gates (all before the slot)
- **KA1:** the module reproduces the H1 pricing runs bit for bit (`test_dn_mag631.py`, 3 seeds × DN-only / UPDN;
  `work/chain_ka631.log`).
- **(d) Fail-safe count:** the build's audit must show dn_invariant_failed = 0 and dn_guard_failed = 0 in every context
  (the ka630_audit style). Any trigger means a KA FAIL.
- **Capacity floor** (A/B/C, `capacity631.py` on 446bj2z's real blocks, 60 targets/context): each context's mean
  dn_frac ≥ 0.10 and dn_moved ≥ 200 genes/target, else re-diagnose, don't submit. (H1 after 630: dn_frac 0.28; A/B/C
  numbers to be filled in from `work/capacity_chain631.log` before freezing.)

## H1 price (definitive, shipping semantics; `work/price-dn2.log`; 3 seeds; DN on top of 630 v3)
| | f 1 | f 0.5 | f 0.35 | cmp |
|---|---|---|---|---|
| Δnmae raw, UPDN − v3 | −0.0005 | −0.0012 | −0.0030 | −0.0058 |
| overall-eq, UPDN − v3 | +0.0000 | +0.0001 | +0.0008 | +0.0015 |
| Δreach raw, UPDN − v3 | −0.0007 | −0.0014 | −0.0001 | −0.0005 |
fid and jac exactly 0. DN delivery after 630: 0.28 of the ×2 need.

## Board prediction (DN vs 630's entry)
H1 f 0.35 Δnmae −0.0030 raw, times 611's delivery multiplier (1.2–1.6), times the footprint factor for called DOWN genes
on A/B/C (unmeasured; take 0.5–1):
**Δnmae raw −0.0018 … −0.0048**, centre ≈ −0.003 (≈ +0.0008 overall). Every other member pinned except reach (below).

## (a) Bins (primary statistic: DN's board Δnmae raw vs 630's entry)
Overall equivalents use 192's slope (scaled nmae per raw unit 1.64787, / 6).
- **≤ −0.0025 (≥ +0.0007 overall): TRANSFERS.** Adopt for A/B/C arms. October adopts via a def_day pass-through (after
  10-08), with the written reason above.
- **−0.0025 … −0.0008: PARTIAL.** Adopt for A/B/C. October adopts only if the build's capacity floor holds on D/E/F's
  own blocks, flagged.
- **> −0.0008 … +0.0008: NULL.** Do not adopt; record the H1 → board ratio for DOWN-magnitude levers.
- **> +0.0008: HARM.** Wrong sign. Close the route and diagnose (footprint, regime).

## (b) Falsifiers (any one: no adoption, diagnose)
- |Δfid| or |Δjac| > 1e-6 raw (a call changed; the guard should prevent it);
- **|Δreach| > 0.0014 raw** (the H1 maximum; reach moves because called genes' p shift within the U bar);
- |Δpds| > 1e-6 raw; |Δmse| > 1e-6 raw. mse is pinned by cell-eval2's #247 cap as for 630: measured, not structural.
  Report the exact value.

## Try count (11.1c.4)
- DOWN magnitude forms tried on H1: m = 2 only (proxy631 cross-fit chose it on every half and regime).
- Delivery forms: DN, UPDN, v7r3, v7r5, C5, v7+DN (running).
- Instruments: 631f (FAILED retrodiction), 614 (FAILED trust gate).
- Board: 0.

## Stacking note
654 (7772dc) finds a global carrier amplitude ×1.25 (with zero repair) also helps on H1 (one seed, +0.0011 at f 0.35).
It moves calls (+26 per target), so it is a call-set change. DN's ×2 target is relative to the delivered lfc, so DN on a
×1.25 base may overshoot. The stacking test (DN on 654's ×1.25 zr dumps) must read before both are adopted together.

## Addendum (2026-10-02T20:25Z, before any stacking output exists): stacking-test prediction on 654's ×1.25 + zero-repair base
Arms (`work/chain_a125.sh`, H1 seeds 0–2, priced against the ×1.25 zs base):
- 630 v3 on ×1.25: overall-eq f 0.35 **+0.0024 … +0.0032**. Additive: it acts on uncalled UP genes, and ×1.25 moved
  them little (654: uncalled UP share −0.0014).
- (630 + DN) − 630 on ×1.25: **+0.0002 … +0.0008** at f 0.35 (standalone on ×1: +0.0008). Expected smaller: ×1.25 already
  lifted head-DOWN delivered |lfc| by ~22 %, and DN's ×2 target is relative to it.
- DN alone on ×1.25: +0.0002 … +0.0008 at f 0.35.
Rule (already frozen): if DN adds < ½ its standalone gain (< +0.0004 at f 0.35) on ×1.25, the predicted board range for
DNMAG halves before its read whenever ×1.25 is in the base.

## Gate readings (2026-10-02T20:27Z; appended, bins unchanged)
- **KA1 (module vs pricing run), UPDN:** `src_dn_mag631.apply` after `src_rank_up.apply` reproduces rnup630's
  `b_k100_q20r3p10_D200` stats bit for bit on 85/85 targets, seeds 0 / 1 / 2 (module sha 0a8cc560…;
  `work/ka631-updn-s*.log`). PASS.
- **KA1, DN alone (clean re-run, `work/chain_dn3.log`): PASS bit for bit, 85/85 × 3 seeds.** Clean DN-alone price
  (`work/price-dn3.log`), overall-eq f 1 / f 0.5 / f 0.35 / cmp: +0.0001 / +0.0002 / +0.0007 / +0.0015 (the contaminated
  first run read +0.0002 / +0.0004 / +0.0008 / +0.0018); delivery 0.128; UP moves 0.
- (History) KA1, DN alone: the first run FAILED because the driver's K = 0 arm was not a true no-UP arm: float-noise asks promoted
  ~425 carriers per target. The driver is fixed (K = 0 now selects nothing), and the clean re-run plus its KA are running
  (`work/chain_dn3.sh`). This affects only the informational DN-alone price; the frozen prereg uses UPDN − v3.
- **Capacity floor (A/B/C, 446bj2z's real blocks, 630 v3 then DN d_m200, 60 targets/context; `work/capacity_chain631.log`):**
  - A: dn_frac 0.196, 1,529 genes moved/target;
  - B: 0.220, 1,483;
  - C: 0.185, 1,536.
  All pass the floor (≥ 0.10 and ≥ 200), at 0.66–0.79× H1's 0.28. **Fail-safe count 0 in every context**; 630's own
  counters stay at 0.

## Stacking-test reading (2026-10-02T21:34Z; `work/price-a125.log`, 3 seeds; zs KA true on every seed): the interaction rule TRIGGERS
On 654's ×1.25 + zero-repair base, overall-eq f 1 / f 0.5 / f 0.35 / cmp:
- 630 v3: +0.0022 / +0.0033 / +0.0033 / +0.0037. Registered +0.0024 … +0.0032 at f 0.35: additive, at the top.
- 630 + DN: +0.0017 / +0.0028 / +0.0030 / +0.0046, so **DN's increment is −0.0005 / −0.0005 / −0.0003 / +0.0009**.
  Registered +0.0002 … +0.0008 at f 0.35 ✗.
- DN alone: −0.0005 / −0.0006 / −0.0004 / +0.0008. Seed 0 started before the driver's K = 0 fix (414 spurious
  promotions/target); seeds 1–2 are clean. All three read alike (−0.00038 / −0.00042 / −0.00042 at f 0.35).

**DNMAG and ×1.25 are SUBSTITUTES.** Both raise DOWN magnitude, and DN's ×2 target is relative to the deeper ×1.25
delivery, so it overshoots. By the frozen rule, DN adds < ½ its standalone gain on ×1.25 (it is negative). **DNMAG must
not be read, or adopted, on a base that carries ×1.25.**

Choice for the coordinator (one of the two):
- ×1.25 + 630: ≈ +0.0052 (f 0.35) / +0.0049 (cmp) over the plain base. A call-set change (+26 calls/target; 611's
  call-set rule applies).
- ×1 + 630 + DN: +0.0037 / +0.0049. Pure delivery (fid / jac exactly 0, reach ≤ 0.0014 raw), an exact read.

## STATUS (2026-10-02T21:35Z, the coordinator's decision): FALLBACK
Path (a) was chosen: α ×1.25 + 630 (larger at f 0.35 / 0.5, one setting with no code change; 654's fid / jac ~0 under the
call-set rule). Its conditions: 618's held-out read confirms ×1.25, and the α-path KA passes, including 554's dose guard on
A/B/C. **DNMAG on a ×1 base is the FALLBACK.** It is read only if α fails those conditions or reads NULL / HARM. Only one
of the two is ever read. This prereg, patch_builder631 and the module stay frozen and ready.
Queue (NEXT-SLOT-DECISION): 630 → 524dr2 → ALLSW → α → idall3 → 625.

## STATUS AMENDED (2026-10-02T21:36Z, the coordinator, before any 618 read): a LIVE candidate again
The pick between (a) α ×1.25 + 630 and (b) ×1 + 630 + DN is made on 618's fresh-H1 surface. Both stacks are priced
there, and the higher minimum over {f 0.35, f 0.5} wins; a tie within 0.0005 goes to (b), which has no capacity risk.
(a) still needs 554's dose guard. Only the winner is ever read on the board.

## STATUS (2026-10-03T17:04Z): (b) CHOSEN on 618's fresh-H1 surface ([21dd4c], 08ebeed1). DNMAG is the live arm.
618, 3 seeds, overall-eq vs plain (f 1 / 0.5 / 0.35 / cmp):
- (b) 630 + DN: +0.0031 / +0.0042 / +0.0040 / +0.0042;
- 630 alone: +0.0028 / +0.0033 / +0.0027 / +0.0030;
- (a) ×1.25 + 630: +0.0036 / +0.0036 / +0.0037 / +0.0052.
The coordinator's rule (higher min over {f 0.35, f 0.5}; ties within 0.0005 to (b)) picks (b): min difference −0.00044.
**DN's increment over 630 on 618: +0.0003 / +0.0009 / +0.0013 / +0.0012**, at or above H1's (+0.0000 / +0.0001 / +0.0008 /
+0.0015). That is a second, fresh surface where it adds. rnup630 ran unchanged: zs KA on all 126-block runs, numeric
block order and target assert held, 0 calls moved. 618 also found ×1.25 is not bulk-neutral in that harness (price630's
sums assert).
Arm, bins, falsifiers and gates as frozen above. The board read is against 630's entry, in the coordinator's queue.

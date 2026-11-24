# 631 — where our nmae loss sits: called DOWN genes are under-magnitude (×2 is worth −0.023 / −0.025 raw on H1)

> **VOID ROWS (2026-10-02T19:10Z; [21dd4c]'s 616 audit, 7dac7a98, `experiments/616-clean-orion-folds/AUDIT-orion-dependence.md`):** the thinned Orion gate tables' DE truth was ~98 % selection artefact (HEK293T: 37,748 calls on the thinned table vs 774 on clean raw cells; 13 / 132 targets scoreable). In this experiment, the Orion rows of 631f are VOID (they were descriptive only). Every H1 number (decomp, proxy, cap, DN, UPDN, v7) is unaffected.

2026-10-02, council session [d4cd5f] (631 is in its 630–639 range). Started 16:13Z. Class: **scorer mechanics** (how DOWN
magnitude is written for called genes; nothing fitted to A/B/C), with a held-out dose check owed before any October use.

## Why
On the board's best entries our nmae is scaled 0.168 against 0.247–0.288 for teams #2–#5 (~0.06 raw), the largest
member gap; fid is second (0.006 vs 0.033–0.083). 630 covers uncalled UP genes (−0.009 … −0.015 raw expected).
`decomp631.py` decomposes the rest.

## Decomposition (`decomp631.py`, `work/decomp631.out`; H1, 446bj2z's zs base, 3 seeds, 521's pools)
Raw nmae, negative = better. "m* cross-fit" = one global factor on the region, chosen on odd targets and scored on even
targets and vice versa.

| region × delivered direction | share of nmae (f 1 / 0.35) | \|T\| oracle at our signs | m* cross-fit (m chosen) | perfect signs |
|---|---|---|---|---|
| **head DOWN** | 0.120 / 0.161 | −0.111 / −0.149 | **−0.017 / −0.022 (m 2, both halves)** | −0.007 / −0.007 |
| head UP | 0.017 / 0.018 | −0.008 / −0.006 | +0.000 | −0.007 |
| block DOWN | 0.159 / 0.198 | −0.040 / −0.048 | −0.006 / +0.009 (unstable: 0.25 vs 3) | −0.033 / −0.032 |
| block UP | 0.042 / 0.048 | −0.004 / −0.006 | −0.000 / +0.001 | −0.019 / −0.015 |
| uncalled UP (630's) | 0.227 / 0.184 | −0.020 / −0.022 | −0.001 / −0.006 | −0.057 / −0.031 |
| uncalled DOWN | 0.336 / 0.280 | −0.077 / −0.054 | +0.001 / −0.003 | −0.066 / −0.043 |

- **Called head DOWN genes have the right signs but far too little magnitude.** Perfect signs buy only −0.007, while
  |T| at our signs buys −0.11 / −0.15. A single factor m = 2 is chosen on both halves at both regimes.
- Block and uncalled regions are sign-limited (perfect signs −0.03 … −0.07), not magnitude-limited. Their magnitude
  levers are unstable or ~0, which matches 603 / 632's "DOWN hurts" for uncalled genes.

## A head-free selection carries it (`proxy631.py`, `work/proxy631.out`)
The builder module sees only the block and the controls, not the design head.

| selection | genes/target | head overlap | m = 2 (f 1 / 0.35) | cross-fit |
|---|---|---|---|---|
| head DOWN (design) | 201 | 1.00 | −0.0168 / −0.0223 | same (m 2) |
| called DOWN, top 200 by \|l\| | 200 | 0.96 | −0.0176 / −0.0231 | same (m 2) |
| **called DOWN, all** (q < 0.05, l < 0) | 1,666 | 0.12 | **−0.0227 / −0.0251** | same (m 2) |
| called DOWN, \|l\| ≥ 0.3 | 79 | 0.99 | −0.0198 / −0.0243 | m 2–3 |

**Doubling the delivered magnitude of every called DOWN gene is worth −0.023 / −0.025 raw nmae on H1 at full delivery**,
≈ +0.006 / +0.007 overall-equivalent. On the board, ×1.2–1.6 (611) would put it at ≈ +0.008 … +0.011, about half our
nmae gap to #2–#5. fid / jac / reach are unchanged as long as calls are.

## Delivery is the constraint (`cap631.py`, `work/cap631.log`; 12 seed-0 targets, head DOWN, |ΔAUC| ≤ 1e-4)
The emitter parks a DOWN gene's bulk total on mass-holder carriers above the control max, and their CPM inflates the
arithmetic mean. Lowering the mean at a fixed total means moving counts carrier → ballast, the mirror of 630.
- Need for m = 2 on head DOWN: ~23k counts per block.
- Carriers can give ~12k with no rank change (holders down to the minimum above the max, plus in-interval room).
- Ballast zero-cost room: ~860.
- With the per-gene U budget (1e-4, upward crossings) ~19 % of each gene's need is reachable.
- The library payback comes from the ballast's own sub-gate counts (~31k per block), which is ample.

The relaxed-budget sweep (1e-3, 3e-3; head and all called DOWN) is running. Reach orders called genes by our adjusted
p, so a relaxed U budget is not free there. It needs pricing on all six members.

## 631f — held-out dose check: the fold model FAILS its retrodiction, so no verdict (16:20Z; `dose631f.py`, gate and rule committed first)
The fold model's h1 row: Δ(m = 2) +0.0153 (f 1) / +0.0046 (f 0.35) raw, argmin m 1 / 1.5. The real emitter says
m = 2 is best (−0.0227 / −0.0251). **Gate FAILED**, so every fold row is descriptive only and the K562 decision is not
binding.

Mechanism of the failure: the fold model's base magnitudes come from BETA, the H1 delivered-to-key ratio over ALL gate
genes. It lacks the real emitter's specific DOWN under-delivery (mass holders inflating called DOWN genes' means), and
that under-delivery is exactly what the lever corrects. A fold model cannot price a delivery defect it does not contain.

**Descriptive** (Δnmae raw at the argmin, relative to the fold base, f 0.35): K562 argmin m 4 (m 3 −0.108); Jurkat
(not clean) m 3 (−0.118); Orion (not clean for this) m 3. Every cancer surface wants ~2× MORE called-DOWN magnitude,
relative to the same model, than H1 does (m 2.5–3 vs H1's 1–1.5). That is consistent with 633: cancer lines' signs are
more often right on these genes.

**Status:** the H1 real-emitter lever (×2) is real on H1. October transfer needs an instrument that contains the real
emitter, ideally the emitter run on a clean cancer surface (K562), which no harness here does yet. Until then the dose
for October is unvalidated (11.1c.2).

## Delivered-level H1 price of the DOWN pass (17:19Z; predictions registered first: PREREGISTRATION-H1.md, 9251902c)
`dn631.dn_pass` (m = 2 target, net |ΔAUC| ≤ 1e-4, ballast sub-gate payback), 3 seeds, `work/price-dn.log`.

| | DN (DOWN only) | registered | UPDN (630 v3 then DN) | registered |
|---|---|---|---|---|
| delivery (dn_frac) | 0.138 | 0.05–0.15 ✓ | **0.280** | DN ± 0.02 ✗ (synergy) |
| Δnmae raw f 1 | −0.0008 | −0.002 … −0.005 ✗ (smaller) | −0.0089 | v3 + DN ± 0.0005 ✓ |
| Δnmae raw f 0.35 | −0.0031 | −0.002 … −0.006 ✓ | −0.0136 | v3 + DN ✓ (−0.0137) |
| overall-eq f 0.35 | +0.0008 | +0.0006 … +0.0017 ✓ | **+0.0037** (v3 +0.0029) | ✓ |
| overall-eq cmp | +0.0018 | — | **+0.0049** (v3 +0.0034) | — |
| fid / jac | 0 / 0 | exactly 0 ✓ | 0 / 0 | ✓ |
| reach raw | −0.00005 … 0 | exactly 0 ✗ (tiny) | −0.0001 … −0.0014 | exactly 0 ✗ |

- **The DOWN lever as delivered is real but modest**: +0.0008 overall-eq at f 0.35 alone, +0.0008 … +0.0015 on top of
  630, and ~0 at f 1, where the gate's weaker genes dominate.
- **Synergy, unpredicted:** after 630's UP pass, the DOWN pass delivers twice as much (0.28 vs 0.14). 630 moves sub-gate
  counts carrier → ballast, and that refills the ballast's payback pool. The prices still add.
- **reach is not exactly 0** (called genes' p shift within the U bar, which reorders reach's walk). The move is small
  (≤ 0.0014 raw, f 0.5). 611 gives no validated reach band for delivery levers, so the board read would carry it.
- **The coordinator's ≤ 25 % rule:** DN alone delivers 14 % (UPDN 28 %) of the ×2 need. Whether an emitter redesign is
  worth it this late is decided jointly. v7 (paired payments) is pricing.

## FINDINGS-ready summary (for the coordinator's FINDINGS §0 / §2)
**631 ([d4cd5f], 2026-10-02): called DOWN genes are delivered at about half the right magnitude.**
- **Size of the defect:** on H1's real emitter, one cross-fitted factor ×2 on every called DOWN gene (q < 0.05,
  delivered l < 0; no design input needed) is worth −0.023 / −0.025 raw nmae at full delivery (f 1 / f 0.35). The
  signs are right there: perfect signs buy only −0.007.
- **Mechanism (CORRECTED 18:04Z by [7772dc]'s 654, seed 0; my "kf pins head DOWN" prediction FAILED):** design amplitude
  DOES reach head-DOWN delivered magnitude, super-linearly (−15 % at ×0.8, +39 % at ×1.25, f 1), and ×1 is the optimum
  of a GLOBAL amplitude. ×1.25's nmae loss is almost all head DOWN: the strongest genes overshoot into the zero-count
  floor before the typical ones lift. So the defect is the SHAPE of the delivered DOWN magnitudes, with the typical
  called DOWN gene under-delivered. A global design knob cannot fix it, and that is what 484/490's null meant. The
  original text below is kept for the record; its "cannot move" claim is wrong.
  - the ballast tops each gene up to the bulk total;
  - the carriers' DE mean follows the DE design, where the mass-balance kf (~0.60) deliberately keeps DOWN genes less
    down, because freed mass would inflate every other gene's per-cell CPM and create false calls (610's coupling);
  - so design-side fixes cannot move the delivered DOWN mean (484/490's "protect" was null);
  - only moving counts carrier → ballast at fixed totals can, and that needs an invisible mass sink for the carriers'
    library payback.
- **Delivered price** (`dn631.dn_pass`; net |ΔAUC| ≤ 1e-4, 0 call flips, libraries/totals exact; payback from the
  ballast's untested genes), H1, 3 seeds:
  - alone: delivery 0.14 of the ×2 need; +0.0002 / +0.0008 / +0.0018 overall-eq (f 1 / f 0.35 / cmp);
  - after 630: delivery 0.28; +0.0023 / +0.0037 / +0.0049, against 630's +0.0023 / +0.0029 / +0.0034.
- **Synergy, a mechanism:** 630 pays its UP receipts by moving untested-gene counts carrier → ballast, which refills
  exactly the pool the DOWN pass pays back from. The two passes are mirror images that share one currency, so running
  both doubles the DOWN pass's delivery.
- **Both passes are bounded by untested-gene currency.** Carriers hold ~0.55 % of their library in it; the ballast holds
  its own pool. Capacity beyond that needs either pairing (v7) or an emitter that writes a larger invisible sink.
- **Transfer:** scorer mechanics. Its dose is fitted on H1. A fold-model held-out check failed its retrodiction (the
  fold model cannot contain the emitter's delivery defect); descriptively, cancer folds want ~2× H1's DOWN magnitude.
  A real-emitter cancer-line read (614) is needed before October.
- **Why a targeted pass, not an amplitude knob:** 654 shows a global ×1.25 loses −0.0041 overall at f 0.35, nmae +0.0147
  from head-DOWN overshoot. proxy631 locates the gain in the moderate called DOWN genes: the top ~80–100 by |l| per target
  carry most of it, and the ~1 gene/target with |l| ≥ 1 is better left alone (cross-fit +0.004 / +0.013). The DOWN pass
  targets ×2 in delivered-lfc space per gene and stops at each gene's own target.
- **Plan (coordinator, NEXT-SLOT-DECISION):** emitter-level work only if v7 lifts delivery well beyond 0.28. Otherwise
  DN becomes a 10-08+ rider on 630, read against 630's own entry with raw-unit bins (only nmae moves; reach beyond
  0.0014 raw is a falsifier), after idall3 and before 625, with 614's read if 614 passes.

## 614 (the K562 real-emitter DE harness) FAILED its trust gate (17:47Z; [21dd4c], 4db74d02)
- zs retrodicts (−0.0041 nmae).
- The coarse MDE comes out at ~1/22 of its board effect and unresolved across seeds; the call-set pair's fid is noise.
- Cause: a low-power truth (median 5 DE genes per target at ≤ 250 cells and 11k UMI). Levers that act on weakly detected
  genes have almost nothing gated to act on; 630's UP pass is in that class.
**K562 rows are descriptive only, and 614 is not used to price 630 or 631.** The DE channel still has no cancer-line
held-out instrument (534's gap stands). 630's and 631's October cases rest on H1 plus their A/B/C board reads. [21dd4c]'s
H1@11k control will show whether depth alone shrinks these levers, which matters if D/E/F's depth differs from A/B/C's.
**Floor check ([7772dc]'s question, seeds 0–1):** the DOWN pass creates no new extreme lfc. Cells below −5 stay at 2 per
seed in DN and UPDN, the base's two zero-total cells, which lie outside the pools and are never touched. Only 4–6 changed
cells per seed started below −2, and the largest |Δlfc| is 1.6–2.2. Totals are fixed, so no zero total can arise:
production's zero repair (VCC_DE_ZERO_REPAIR=−5, absent from the 537 harness) is an exact no-op for these arms, and the
harness-vs-production difference 654 found does not touch 631's prices.

## v7 (630's UP receipts paid with budgeted called-DOWN counts): delivered-level H1 price (18:12Z; predictions registered first)
`up7_630.py`, `work/price-v7.log`, 3 seeds. ~15.6k of ~48k payment counts per target came from DOWN genes (32 %).
Overall-equivalent:

| arm | f 1 | f 0.5 | f 0.35 | cmp | registered |
|---|---|---|---|---|---|
| v3 (630) | +0.0023 | +0.0030 | +0.0029 | +0.0034 | |
| C5 (v3, 5 rounds) | +0.0023 | +0.0030 | +0.0029 | +0.0034 | ±0.0003 of v3 ✓ |
| **V7r3** | **+0.0010** | +0.0038 | **+0.0043** | +0.0043 | f 0.35 +0.0000 … +0.0008 over v3 ✗ (above: +0.0014) |
| V7r5 | +0.0010 | +0.0038 | +0.0043 | +0.0044 | +0.0002 … +0.0015 over C5 ✓ (+0.0014) |
| UPDN (630 + DN, for reference) | +0.0023 | +0.0030 | +0.0037 | +0.0049 | |

- fid / jac exactly 0; reach ~0 (+0.000005), cleaner than UPDN's.
- Floor check: median DOWN Δlfc −0.11; at most one new cell below −5 per seed; totals fixed.
- **Reading:** pairing pays at f 0.35 / 0.5, the regime 611 validates for delivery levers (+0.0014 over 630, +0.0006 over
  UPDN). It loses at f 1 (−0.0013): payments go to whichever called DOWN genes the recipient carrier holds, including
  weak ones where ×2 overshoots the full-power gate. It does not lift DOWN delivery "well beyond 0.28", so by the
  coordinator's rule there is no emitter-level work.
- Last combination, registered now: **v7 then DN** (`b_k100_q20r3p10_D200_P200`). Prediction: f 0.35 +0.0048 …
  +0.0052 (v7 + DN-alone's +0.0005 … +0.0008), f 1 +0.0010 … +0.0013, DN delivery after v7 ≥ 0.20, fid / jac 0.

## v7 + DN: overshoot, rejected (20:25Z; registered first; `work/price-v7dn.log`, 3 seeds)
| arm | f 1 | f 0.5 | f 0.35 | cmp |
|---|---|---|---|---|
| v3 | +0.0023 | +0.0030 | +0.0029 | +0.0034 |
| UPDN (the frozen rider) | +0.0023 | +0.0030 | +0.0037 | +0.0049 |
| v7 | +0.0010 | +0.0038 | +0.0043 | +0.0043 |
| **v7 + DN** | **−0.0008** | +0.0020 | +0.0032 | +0.0039 |

Registered prediction: f 0.35 +0.0048 … +0.0052 ✗; f 1 +0.0010 … +0.0013 ✗. v7's DOWN payments already deepen called
DOWN genes (median Δlfc −0.11), and DN then targets ×2 of that DEEPER lfc. Its relative target compounds into overshoot
(nmae +0.0025 raw at f 1). **Lesson: DN's target is relative to the delivered lfc it finds, so it must not stack on
anything that already raised DOWN magnitude without re-pricing.** That is exactly the ×1.25 question (stacking test
running). Switch rule: v7 + DN does not beat UPDN without an f 1 loss, so UPDN stays the rider.

## Stacking test on 654's ×1.25: DN and ×1.25 are substitutes (2026-10-02T21:34Z)
630 stays additive on ×1.25 (+0.0033 at f 0.35). DN on top overshoots (increment −0.0005 f 1, −0.0003 f 0.35, +0.0009
cmp), because its ×2 target is relative to the deeper ×1.25 delivery. This is the same compounding as v7 + DN. So
DNMAG and ×1.25 are alternatives: ×1.25 + 630 ≈ +0.0052 at f 0.35 (a call-set change) vs ×1 + 630 + DN +0.0037 (pure
delivery); they tie at cmp. The choice is the coordinator's (DNMAG prereg, stacking reading).

**Clean DN-alone (2026-10-02 ~21:35Z, driver K = 0 fixed; module KA 3/3 bit for bit):** +0.0001 / +0.0002 / +0.0007 /
+0.0015 overall-eq (f 1 / 0.5 / 0.35 / cmp), delivery 0.128. The first DN-alone run carried ~425 spurious UP promotions per
target and read slightly higher (+0.0008 at f 0.35). UPDN, the rider's priced form, was never affected.

## 618 picks 630 + DN (2026-10-03T17:04Z; [21dd4c], 08ebeed1)
On 618's fresh-H1 surface (3 seeds), DN adds +0.0003 / +0.0009 / +0.0013 / +0.0012 over 630 (f 1 / 0.5 / 0.35 / cmp). The
coordinator's pick rule chooses (b) 630 + DN over (a) ×1.25 + 630, so DNMAG is the live A/B/C arm, read against 630's
entry under its frozen prereg.

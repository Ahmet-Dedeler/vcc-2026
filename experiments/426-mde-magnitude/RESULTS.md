# 426 — MDE-structured magnitude on the DE channel: a pure-nmae lever, nmae −0.0096 on H1 at 5 seeds, in production as `VCC_DE_MDE`

Status: H1 complete 2026-09-23 04:25Z (PI). Production patch committed `8dc263b3` with known answers. Board arm 426a
(417e + `VCC_DE_MDE=0.1,0.1,0.3`) building for the 09-24 slot 2; preregistration in `PREREGISTRATION.md`.

## Mechanism (from 425 §2)
- nmae is scored only on the reference's significant genes. Given significance, |lfc| is bounded below by the minimum
  detectable effect, MDE ≈ 0.24/sqrt(mean control count per cell).
- A predicted magnitude x < |l| changes the expected L1 error by −(2a − 1)·x, where a is the gene's sign accuracy. So
  adding magnitude up to about the MDE is a sure gain whenever a > 1/2.
- Our delivered magnitudes are ~0.17× the truth, far below the MDE.
- **Implementation (A-channel only).** The K-class emitter's prior mode sets each gene's arithmetic mean A through its
  class split at a fixed per-gene total, then designs the Wilcoxon |z| per class.
  - The add, sign(target) × min(c/sqrt(count), cap) (0 below min_count), enters as an extra steering move on every
    gated gene (`mde_add` in the prior spec).
  - It is mass-balanced by one factor on the heavier side, because the gated genes' arithmetic mass is fixed by the
    class levels.
  - It gets its own feasible dose after the prior's dose is fixed.
  - Per-gene totals (the bulk) are untouched. The |z| targets are unchanged; the realised |z| drifts only where the new
    split clips the design.

## Harness and known answers
- `run426.py` is 416's `run416.py` (ballast, ratio mode, jk boost) with `--mde-*` flags; `k_emit426.py` is 416's emitter
  plus the `mde_add` block.
- **KA0:** knob off reproduces 416's `base_s0_bps` (the shipped 410b construction). Every stats array is bit-equal;
  members agree to 1 ulp (the scorer's polars aggregation, as in 419).
- **Production (`ka426prod.py`, `work/known-answer-prod.json`):**
  - KA1: the patched `src/vcc2026/rc_emit_k.py` with the knob off equals 416's base, bit-equal.
  - KA2: with the knob on, it equals `k_emit426`'s run, bit-equal.
- `test_mde426.py`: the builder knob parses (and is None unset), with 4 refusals.

## Results (H1, Δ vs the same-seed base; board slopes for overall)

**The chosen arm, c 0.1, cap 0.1, min_count 0.3, 5 seeds:**

| member | Δ raw (mean ± sd) | → overall |
|---|---:|---:|
| nmae | **−0.00962 ± 0.00016** | **+0.00265** |
| fid | −0.00030 ± 0.00004 | −0.00017 |
| jac | −0.00015 ± 0.00002 | −0.00007 |
| reach | −0.0019 ± 0.0023 | −0.00035 |
| pds, mse | 0 (to 5 decimals) | 0 |
| **total** | | **+0.00206 ± 0.00041** (non-reach +0.0024) |

**The menu (seed 0 unless marked):**

| variant | Δnmae | Δfid | Δjac | Δreach | unsuppressible calls/target | clip frac | ovr |
|---|---:|---:|---:|---:|---:|---:|---:|
| c 0.1 all, uncapped | −0.0119 | +0.0026 | −0.0073 | −0.030 | 2,315 | 0.53 | −0.0042 |
| c 0.05 all, uncapped | −0.0072 | −0.0001 | −0.0001 | −0.011 | 0.8 | 0.18 | −0.0003 |
| c 0.1 block only | −0.0012 | −0.0002 | −0.0001 | −0.004 | 0.6 | 0.05 | −0.0006 |
| c 0.05 non-block only | −0.0050 | 0.0000 | 0.0000 | −0.006 | 1.8 | 0.21 | +0.0002 |
| **c 0.1 cap 0.1 min 0.3** | **−0.0098** | −0.0003 | −0.0002 | −0.005 | 3.4 | 0.28 | +0.0014 (5-seed +0.0021) |
| c 0.05 cap 0.1 min 0.3 | −0.0065 | −0.0001 | −0.0001 | −0.006 | 0.6 | 0.06 | +0.0005 |
| c 0.1 cap 0.1 min 0.3, concordant only | −0.0106 | +0.0063 | +0.0011 | −0.005 | 522 | 0.31 | (extra calls: not transferable) |
| base (knob off) | 0 | 0 | 0 | 0 | 0.6 | 0.035 | 0 |

**Reading.**
- **The analytic prize is delivered.**
  - Analytic, at the delivered LFC with caps: −0.0107.
  - Through the real emitter: −0.0098 on one seed and −0.0096 over five, with every other member pinned.
- **Uncapped low-count adds break the z design.** Genes below ~0.3 counts/cell need lopsided class splits, which clip
  the |z| targets. Complement genes are then called ("unsuppressible", thousands per target) and jac and reach fall.
  - The cap (0.1 log2) and the count floor (0.3) remove this, leaving 3.4 extra calls per target (base 0.6).
- **Block-only is small** (−0.0012): the CPM block holds the most expressed genes, where the add is small.
- **The reach cost is inside noise:** −0.0019 ± 0.0023 over 5 seeds. It moves on ~18 of 82 targets, in discrete
  small-n_real jumps.
  - Its plausible source is the reach order among our calls. Spearman of log p over shared calls vs base is 0.95,
    because the |z| design clips slightly more under the new splits.
  - The concordant restriction (add only where the LFC sign agrees with the |z| sign, i.e. not on prior-flipped
    genes) raises that Spearman to 0.99, but it creates extra calls, so it is not the shipped variant.
- **Signs:** 0.17 % of delivered LFC signs change (78 of ~221k calls).

## Dose and sign follow-ups (04:25–04:55Z)
- **The dose is at the emitter's edge.**
  - c 0.15 / cap 0.15 / min 0.3 gives nmae −0.0101 on seed 0, barely beyond c 0.1's −0.0098. But it calls 4,060
    unsuppressible genes per target (clip 0.52), so jac −0.0048, fid −0.0028: **−0.0020 overall**.
  - The class split cannot move A further without breaking the |z| design. **c 0.1 / cap 0.1 / min 0.3 stays the
    shipped dose.**
- **The CLEAN4 prior's signs on non-block genes (`--mde-prior-sign`) do not fit the class structure.**
  - About 2,200 extra flips per target shift enough mass that one H1 block over-filled a class by 8,324 counts against
    14,430 sub-gate counts. The per-cell repair stalled: "did not meet the per-cell totals".
  - With an over-fill guard at 0.5× the sub-gate mass, the MDE dose collapses: λ median 0.18, nmae −0.0004. Shelved.
  - The analytic prize (−0.0128 vs −0.0098) needs another delivery: fewer, confident flips, or flips paid inside the
    class.
- **Prior signs through the A-channel are closed, even with an over-fill guard at 0.55.**
  - All covered genes: 2,220 flips per target; the MDE dose collapses (λ median 0.21, nmae −0.0005).
  - Confident genes only (|P−0.5| > 0.1): 184 flips per target, λ 0.76, nmae −0.0072, which is worse than the plain
    MDE. It also adds 103 unsuppressible calls per target, which carry an untransferable fid +0.0022.
  - Sign changes on uncalled genes do not fit the class structure. Delivering better signs needs a carrier-design
    route (427e), with its board risk.
- **Agreement gate (`--mde-agree`, `priorsign426.py` 426e): no clear gain.** The add is dropped where our sign disagrees
  with the covered CLEAN4 prior. On those genes our sign is right 0.45 of the time on H1, so magnitude there costs.
  - Analytic prize: −0.0121 vs −0.0098.
  - Through the emitter, 3 seeds: nmae −0.0107 (−0.0011 better than plain), but reach −0.0062 ± 0.0019 vs plain
    −0.0031 ± 0.0022, and fid −0.0004. Unsuppressible calls 0.7 and clip 0.16 are cleaner than plain's.
  - Total +0.0015 vs plain +0.0019. **Not shipped.** Revisit if the board shows reach is insensitive to the MDE.
- **Relative targets (`--mde-rel`, 426f) change nothing on H1.**
  - Non-block MDE targets taken from the unsteered split's delivered A instead of the ideal A0 give nmae −0.00972 vs
    −0.00976. The over-fill distribution is the same (median 0.48, max 0.60).
  - Even H1 has a few blocks at λ = 0: the per-target median is 1, the min 0.
  - The single-count-polish hypothesis for A/B/C's 0.84 dose is not supported on H1. The per-block A/B/C diagnostics
    (431's `VCC_RC_AUDIT_DUMP` in the next build) will say what binds.
- **Context A gets no MDE: diagnosed with the per-block dump (431, 429a's build) and fixed on H1 (426h).**
  - 429a's dump: context A `mde_lambda` median **0.000** (92 % of blocks at 0, none at 1); B median 0.935. So 426a's
    board read is effectively B and C only.
  - Cause: the MDE on non-block genes shifts ~6,400 counts between library classes (H1, `mde_dclass_*`).
    - `split_for_A`'s single-count polish takes its fine steps from the MIDDLE class for both up and down genes.
      Class 2 is drained (−6,415 net) while classes 1 and 3 gain.
    - Context A has only ~77 sub-gate counts per cell, so ~4k of slack under the 6,000-count capacity margin
      (H1 ~8.6k). The check fails at any positive dose.
  - **Fix (`mde_coarse`):** non-block MDE moves use only the extreme-class channel, whole counts, no polish, targets
    relative to the unsteered split (`split_coarse_extremes`). An A-balanced set of moves is then count-balanced.
  - H1, 3 seeds: per-class shift median **~390** (was ~6,800; max ~14k on a few blocks); nmae **−0.00953 ± 0.00021**
    (plain −0.0096); fid −0.0002, jac −0.0001, reach −0.0036 ± 0.0018; total +0.0018. Same prize, 17× less class
    pressure.
- **A/B/C delivery (the 426a build):** `mde_lambda` median 0.84 (min 0). The prior's class-capacity check binds on many
  A/B/C blocks, against none on H1.
- **Over-fill census for the shipped variant** (c 0.1, seed 0; `mde_overfill_frac` in the diag): over-fill per class /
  sub-gate mass is median 0.48, q90 0.51, max 0.60, and every block placed.
  - The one stall seen was 0.58, with the flips.
  - Production builds run without a guard (the committed emitter). The `mde_overfill` option exists in `k_emit426.py`
    only, and the same c 0.1 run reproduces bit-for-bit without it.
- **423 (worker) is complementary.** Its exact-U, zero-sum carrier↔carrier pass on down-regulated "mass-holder" genes
  gives nmae −0.0026 on H1 with the other five members bit-identical.
  - It works inside the control intervals after emission; the MDE works on the class split before the |z| design.
  - A stacked test (MDE + zs pass) comes after 426a's board read.

## Next (after the board reads 426a)
- If nmae transfers (≤ −0.008 on the board):
  - the MDE rides every later arm;
  - then the stacked MDE + 423 zs pass on H1;
  - then a delivery of the prior signs that fits the class structure (confident flips only, or flips paid inside the
    class).
- **October:** the add needs only the context's own controls (ref_lin and the library median), so it transfers to
  D/E/F as it is. Add a DEF-RUNBOOK row after the board read.

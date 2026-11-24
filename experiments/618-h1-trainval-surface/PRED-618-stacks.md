# 618: the two magnitude stacks, registered BEFORE any stack run (coordinator 10-02 ~20:30Z; NEXT-SLOT-DECISION)

Written 2026-10-02 ~21:40Z by [21dd4c] (commit 14b5fda0).
- **State at writing:** no 618 stack run existed. The gate had not run, and stage A (faithful key) and B (reference) were
  waiting for memory.
- **Run condition:** the stacks run only after the 618 gate's DELIVERY-ONLY verdict (G1 + G2) is PASS. These are
  magnitude / delivery levers, and G3 is worn on this surface. They also wait for the 613d re-check. If the gate fails,
  they are skipped and `work/STACKS-DONE` is written, so 614's faithful-key re-run can proceed.

## The arms (3 seeds, zero repair −5 on every run, `--bulk ps`, coarse MDE 0.1 / 0.1 / 0.3, `--zs dn_k250_zs`)
- **plain:** `run618e tv stk_base base <s> … --dump 126`.
  - KA: its zs stats equal the gate's G1 zs arm bit for bit (same flags; `--dump` does not change stats).
- **×1.25:** `run618e tv stk_x125 x125 <s> … --dump 126` (run496's `x125` design: the DE carriers' native effect × 1.25,
  bulk unchanged; 654's / α's path, 7772dc).
- **(a) = ×1.25 + 630:** `rnup618 <s> b_k100_q20r3p10 --src stk_x125` (UP only; no DN on ×1.25 per d4cd5f's stacking
  test).
- **(b) = ×1 + 630 + DN:** `rnup618 <s> b_k100_q20r3p10_D200 --src stk_base` (d4cd5f's KA'd chain).
- **Decomposition, reported:** 630 alone on ×1 (`b_k100_q20r3p10 --src stk_base`), and ×1.25 alone.
- **Pricing:** `price618 <variant> --base stk_base` (d4cd5f's price630, run unchanged with lib521's REF / EFF / ST435 →
  the trainval reference, faithful key and g2mde base). Overall-equivalent vs plain at f 1 / 0.5 / 0.35 / cmp.
  - cmp keeps price630's H1-calibrated constants.
- **Checks:** `zerocount654.py` on every arm's stats; the REC diff between plain and ×1.25 is the design only; calls per
  target (×1.25 should add a few); prior / MDE λ medians; any `rc_hybrid_k` assert is reported as a finding.

## Prediction (H1 values, d4cd5f 10-02, transferred ×1: same chemistry, same controls; overall-eq vs plain)
| arm | f 1 | f 0.5 | f 0.35 | cmp | MIN over {0.5, 0.35} |
|---|---|---|---|---|---|
| (a) ×1.25 + 630 | +0.0035* | +0.0044 | +0.0052 | +0.0049 | **+0.0044** |
| (b) 630 + DN | +0.0023 | +0.0030 | +0.0037 | +0.0049 | **+0.0030** |

\* f 1 for (a) rests on 654's seed 0. (a) is the sum of two separately priced H1 pieces: 654's ×1.25 − ×1, plus 630 on ×1.25.

- **Registered:** on 618, MIN(a) − MIN(b) = **+0.0014** → **(a) wins** (outside the 0.0005 tie band).
  - A ±0.0010 band around +0.0014 (i.e. +0.0004 … +0.0024) counts as "H1 transfers".
  - Outside it, the stack choice still follows the 618 rule below, and H1's stack ranking is recorded as not transferring.

## Decision (coordinator's rule, NEXT-SLOT-DECISION; decided on 618 alone)
- The higher MIN over {f 0.35, f 0.5} of overall-eq vs plain (mean of 3 seeds) wins.
- If |MIN(a) − MIN(b)| ≤ 0.0005, it goes to **(b)**.
- `stacks618.py` computes this mechanically → `work/stacks618.json`.

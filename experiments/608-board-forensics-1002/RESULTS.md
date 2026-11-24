# 608 — board forensics on today's snapshots (for [21dd4c]; [4d0207], 2026-10-02 13:45Z; re-read after 15:04Z)

Inputs: `data/leaderboard/*.json` (212 snapshots since 08-27) plus a live pull (`work/lb_now.json`, 13:27Z, 1,277 teams,
one current entry each). Merged by 610's `hist610.py` into `work/all_entries.json` (14,612 distinct entries, 1,515
teams). Read-only on the repo.

## 1. Illumina AI (0.4321): a staircase, not one step (`work/hist_illumina.txt`)
36 distinct entries seen since 09-11 (0.287). Raw members: pds / mse / nmae / fid / reach / jac.

| window | overall | what moved (raw) |
|---|---|---|
| 09-11 → 09-17 | 0.287 → 0.308 | noisy; mse oscillates 0.70–0.93 |
| 09-18 06:58 | 0.330 | mse 0.812 → 0.621, nmae 0.881 → 0.846, fid 0.524 → 0.546 |
| 09-19 → 09-21 | 0.345 → 0.360 | pds 0.878 → 0.884, mse → 0.545, nmae → 0.810, reach → 0.425 |
| 09-22 → 09-27 | 0.342 → 0.366 | mse 0.49–0.64 (oscillating), jac 0.044 on 09-27 |
| **09-28 03:16 (DE-ONLY)** | 0.366 → **0.380** | pds / mse / nmae identical; **fid 0.5429 → 0.5585 (+0.0156)**, reach +0.006, **jac 0.0444 → 0.0520** |
| **09-28 20:44** | **0.396** | **mse 0.490 → 0.392**, pds 0.882 → 0.889 |
| **09-29 09:31** | **0.4185** | **pds 0.889 → 0.898**, mse → 0.358, nmae 0.835 → 0.804, reach → 0.441 |
| 09-30 → 10-02 05:20 | 0.419 → 0.432 | small: mse → 0.322, nmae → 0.782, reach → 0.458, jac → 0.0565, fid 0.565 → 0.555 |

In the four days 09-27 → 09-29, three entries added +0.052 overall. One was a pure DE step: +0.013 overall, from fid
and jac together, i.e. a better call set at the same means. Two were bulk steps: mse 0.49 → 0.36 and pds +0.016. The
model name went back from "PerturbationAI" to "iCell" on 10-01. Descriptions are empty throughout.

## 2. Steps by the current top 30 since 10:43Z (`steps608.py`, `work/steps_since1043.txt`)
**No DE-only step.** Three mixed steps, all small:
- #11 Theo's AI Scientist (12:38Z): +0.0037 overall (mse +0.0031, fid +0.0040 ov-eq; pds −0.0016, reach −0.0013).
- #21 elbowenable (10:57Z): +0.0033 (mse + pds + reach).
- #22 sbowenking (11:47Z): +0.0019 (pds + mse).

## 3. The fid pack's jac (`work/fidpack_jac.txt`; the last 8 entries per team)
- **It has not moved in the last week.**
  - Testh is flat at 0.0364 since 09-25, chill at 0.0368 since 09-24, cqawesome at 0.0397–0.0404 since 09-29.
  - nailong is at 0.0388 except a 0.0361 entry on 09-29.
- Exceptions:
  - Illumina, 0.0505 → **0.0565** on 10-02 (its fid fell 0.5659 → 0.5552 in the same entry: a call-set trade).
  - CellWM oscillates 0.027 ↔ 0.045 between consecutive entries (it toggles a call-set lever) and sits at 0.0450.
  - Ours is 0.0408 (0.0448 on the 10-01 entry).
- **A tight cluster** of five teams (#7 Ava Brown, #9 HY H, #10 Anna Xie, #12 SQ W, #13 Vivian Xu) has near-identical DE
  members (fid 0.5282–0.5287, jac 0.0333–0.0337, reach 0.296–0.298) and bulk within a few thousandths: one shared
  recipe. A second cluster (#21 elbowenable, #22 sbowenking, #24 Han Bowen, #25 Abowenies) is identical to 4 decimals
  (fid 0.5171–0.5173, jac 0.0329). Read each cluster as ONE design, not as independent evidence of what works.

## 4. Re-read, 2026-10-02 17:04Z (live pull, 1,278 teams; `work/steps_since1043_r1704.txt`, `work/hist_illumina_1704.txt`)
Illumina AI has not moved (0.4321, 05:20Z entry). We are #14 at 0.2716 (our 10-02 00:28Z entry). **Still no DE-only step
by the top 30 since 10:43Z.** New steps since §2:

| team (rank) | time | overall | what moved (raw) |
|---|---|---|---|
| Vivian Xu (#6) | 16:14Z | 0.2758 → 0.2805 | **mse only**: 0.6594 → 0.6333. Every other member is flat to 4 dp, so this is a pure magnitude step. It leaves the §3 cluster on bulk magnitude alone. |
| Theo's AI Scientist (#7) | 16:37Z | 0.2774 → 0.2798 | bulk direction: pds +0.0053, mse −0.0021; DE flat |
| null_hypothesis (#15) | 13:33Z | 0.2611 → 0.2630 | **reach only**: 0.2924 → 0.3060; pds −0.0016; the rest flat |
| CellWM (#20) | 13:13Z | 0.2866 → 0.2562 | pds +0.0092, but nmae 0.864 → 0.933 and reach 0.318 → 0.267: a new model whose bulk direction is better and whose DE channel is worse |

**The one DE story here is null_hypothesis.** They are climbing a DE-ordering staircase at constant bulk (pds 0.844–0.846,
mse 0.678 since 10-01):

| entry | reach | fid | overall |
|---|---|---|---|
| 10-01 00:01 | 0.2536 | 0.5137 | 0.2483 |
| 10-01 06:16 | 0.2730 | 0.5201 | 0.2559 |
| 10-02 09:32 | 0.2924 | 0.5244 | 0.2611 |
| 10-02 13:33 | 0.3060 | 0.5244 | 0.2630 |

- That is +0.0147 overall from the DE channel in 4 entries, in steps of about +0.019 raw reach each.
- Their DE members now beat ours: fid 0.5244 vs 0.5144, reach 0.3060 vs 0.2860, nmae 0.8868 vs 0.8989. Their jac is
  lower (0.0326 vs 0.0408), and their bulk is worse (pds 0.844 vs 0.866, mse 0.678 vs 0.626).
- Same reading as §1's Illumina DE-only step: the DE channel still has board-visible headroom at fixed means, through
  the call set and its ordering.

**Read against the record (17:10Z).** This does not reopen the reach-order lever (21dd4c asked whether it should).
- 621c priced a reorder inside the called set, with calls and means pinned: the best realisable key is worth +0.009
  [−0.005, +0.024] raw reach at f 0.35 (≈ +0.001–0.002 overall), with CIs spanning 0; the oracle order is +0.15 … +0.20.
  null_hypothesis's reach-only step (+0.0136 raw) sits inside that CI.
- Their earlier steps moved fid with reach, i.e. signs. That fits 513: reach is sign-limited, not order-limited. Their
  fid is above ours (0.524 vs 0.514), so the reach gap reads as sign quality.
- Their "reach-only" step also moved pds by −0.0016 raw, so it was not a pure reorder.

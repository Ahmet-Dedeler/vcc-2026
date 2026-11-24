# PREREGISTRATION: C1 on A/B/C, the board read of 811's increment on top of R0 + 707's swap

Written 2026-10-18 08:32Z by the clock. This is BEFORE R-ABC2H's board read (10-19 s2) and BEFORE C1 is built (readiness, tag `c1811`,
10-19 ~12:00–17:00Z) or read (10-20 s1). Only R-ABC2H's board number gets filled in later. Nothing else here may change.

## What is read
- **C1** = R-ABC2H's recipe (R0 + 707's strong-tier swap via dumphook707) plus 811's c0_l25 many-line blend, applied after the
  swap.
- Built from the lab's UNAPPLIED `experiments/811-lolo-model/dumphook811.sh` (sha recorded in the READY note), via
  DEF_DAY_DUMP_HOOK on def_day.
- Compared against R-ABC2H, the same recipe without 811.

## The prediction (fixed now)
- **Increment source:** 811's information read of C1 − B1 on ABCp, **+0.00153** [−0.00096, +0.00428] 2/3 (composite811
  `6ec25788…`, `811/RESULTS.md` §000000). LABC +0.00142 is reported beside it.
- **Deflation (fixed now):** ρ = **1.7232**, the lab's bulk-direction calibration on ABCp (`800-checker/ka800.json`,
  `bulk-direction|ABCp`, from 707's board point). The range uses ρ ∈ [1.44, 1.72] (ABCp key / direction).
- **pred (point) = R-ABC2H's board overall + 0.00089**
- pred (ρ range) = R-ABC2H + [0.00089, 0.00106]. The increment's own CI ÷ ρ = [-0.00056, +0.00248]
  (information only).
- **Bins, on D = C1's board overall − pred (point):**
  - |D| ≤ 0.0005 → AS PREDICTED;
  - D > +0.0005 → ABOVE;
  - D ≤ −0.0005 → **HARM**.
- **Second HARM condition:** C1's board overall < R-ABC2H's board overall − 0.0005.
- **Either HARM condition closes 811's 6a route for the final** (runbook Amendment 3), whatever DEFo reads. A non-HARM read
  does not qualify 811 by itself.
- **Noise note:** the increment is bulk-only. Through ballast K, the DE members carry the same-law re-roll sd ≈ 0.00025 overall
  (AGENTS rule 13). That is inside the ±0.0005 bins and is not charged separately.

## Fill-in (after R-ABC2H reads, before C1 is read; filled 2026-10-19 02:18Z by the clock, the only edit this file allows)
- R-ABC2H board overall: **0.2748936** (entry `ndt8swt7ILYn3qMJYZ9a`, captured 2026-10-19T01:49:28Z,
  `data/board-entries/ndt8swt7ILYn3qMJYZ9a.json`; production's own pred was 0.2750179)
- **pred (point) = 0.2748936 + 0.000888 = 0.2757814** (ρ range: 0.2757814 … 0.2759580)
- AS PREDICTED: C1 in [0.2752814, 0.2762814]; ABOVE: C1 > 0.2762814
- **HARM line 1:** C1 ≤ pred − 0.0005 = **0.2752814**
- **HARM line 2:** C1 < R-ABC2H − 0.0005 = **0.2743936**
- (Note, not a change: across 707 and 710 the A/B/C-panel folds predicted pds losses the board did not show. If that holds
  here, C1 may read ABOVE. Recorded before the read.)

## Read (2026-10-20 02:37Z by the clock; from the coordinator [3c7a7e] and BOARD-LEDGER)
- **C1 board overall: 0.2777649** (entry `RMxJtslBpY2GkwIMHOiB`, 10-20 s1). **D = C1 − pred = +0.0019835 → ABOVE.** Not HARM: C1 − R-ABC2H =
  +0.0028713 ≥ −0.0005.
- Increment over R-ABC2H +0.0028713, all bulk: Δmse +0.01773 scaled, Δpds −0.00059 scaled, DE +0.000015 (re-roll).
- **The board gave 811's increment 1.9× ABCp's checker read** (+0.00153), i.e. an implied ρ = 0.00153 / 0.0028713 = **0.533**.
  The lab's deflation (ρ 1.72, from 707) does not hold for this construction. One read: it counts against that deflation, and
  it does not by itself recalibrate ABCp. It matches the pds-pessimism pattern logged after 707 and 710: the board's pds cost
  (−0.00059 scaled) is far smaller than the folds' literal-pds leg predicted.
- 811's 6a route stays open on the swap-ON branch (Amendment 3's HARM did not fire).

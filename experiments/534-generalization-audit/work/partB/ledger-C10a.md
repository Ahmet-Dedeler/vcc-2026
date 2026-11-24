### Outcome C10a (the single R6 access, 2026-09-29, by the 534 Part B worker): the October bulk recipe keeps 70 % of its A/B/C gain on never-used panels; every milestone TRANSFERS; the LODO key is the weak switch

*Draft for the owner to append under "Claim C10a" in `docs/CONFIRMATION-LEDGER.md`. Numbers from
`experiments/534-generalization-audit/work/partB/results534B.json`; write-up in `experiments/534-generalization-audit/RESULTS-B.md`.*

- **Access.** R6 was drawn once by `draw534.py` at 2026-09-29T15:50:51Z, after the efp known answers had passed. The draws
  and their sha256 are in `work/partB/draws.json`:

  | panel | rule | seed | pool (essential-free) | sha256 prefix | K562 / HCT116 / HEK293T |
  |---|---|---|---|---|---|
  | S1 | efr | 534101 | 5,599 | 93237a3a | 300 / 300 / 300 |
  | S2 | efp | 534201 | 14,016 | 835d4b32 | 137 / 291 / 293 |
  | S3 | efr | 534102 | 5,184 | d9112014 | 300 / 300 / 300 |
  | S4 | efp | 534202 | 13,416 | 88a7802a | 126 / 293 / 293 |

  - Every panel has 300 targets, 0 essential, and 0 overlap with the 9 panel files on disk (7 distinct panels, 1,829
    excluded targets) or with each other.
  - All four panels were priced: 12 held-out cancer folds, 16 fresh operators. Each operator passed `_check_axis` on its
    own panel, and its axis differs from A/B/C's and from all 16 tuned-on (panel × donor set) axes.
  - The scores were read once, by `aggregate534.py`: 4,000 draws, seed 534. No fit, selection, dose or debugging step
    read an R6 score.
- **Known answers (efp, before the draw): all passed with |d| = 0.**
  - 438's stored rows reproduce exactly: M0–M3 and C-a, design and realised, both levelled twins, 3/3 folds. So do
    493's M4 row, 546o's M5 row and 522's `pctPROD` row for M6.
  - The fold caches are bit-identical, and so are the full-panel energies and level factors.
  - 72 prices from 438's own pricer match exactly; the H1 ×2 half on efp is +0.002545904997, as recorded.
- **Result.** 438's levelled route, k = 1, pooled over R6's 12 folds, 95 % paired target-bootstrap intervals. The
  A/B/C denominators are the record's values.

  | step | R6 pooled [95 % CI] (folds positive) | verdict | retention vs A/B/C |
  |---|---|---|---|
  | M0→M1, q 1.25 | +0.00652 [+0.00593, +0.00727] (12/12) | TRANSFERS | 0.55 |
  | M1→M2, H1 ×2 | +0.00099 [+0.00037, +0.00150] (8/12) | TRANSFERS | 0.42 |
  | M2→M3, q 1.5 | +0.00458 [+0.00403, +0.00526] (12/12) | TRANSFERS | 0.53 |
  | M3→M4, q 1.75 | +0.00282 [+0.00232, +0.00342] (12/12) | TRANSFERS | 0.61 |
  | M4→M5, LODO key | +0.00099 [+0.00033, +0.00181] (7/12) | TRANSFERS | 0.30 [0.10, 0.55] |
  | M5→M6, percentile cap | +0.00694 [+0.00517, +0.00828] (10/12) | TRANSFERS | 3.31 |

  - **Cumulative M0 → M6: +0.0229 [+0.0204, +0.0255] against +0.0327 on A/B/C, retention 0.70 [0.62, 0.78].** Against
    the exact designs priced on A/B/C's own panel (supplementary) it is 0.68 [0.61, 0.76].
  - By panel: S1 +0.0305, S2 +0.0250, S3 +0.0221, S4 +0.0138.
  - Calibrated board-unit projection: +0.0182 [+0.0157, +0.0206], against +0.0234 for the same designs on A/B/C.
- **Rules (PREREGISTRATION-B §7).**
  - No switch FAILS or is FRAGILE, so nothing is a switch-off candidate.
  - Tuned-on check: 1 of 3 October switches reads below half of its real/efp/efr mean. The LODO key reads +0.00099
    against half-mean 0.00147; q 1.75 and the percentile cap do not. So real/efp/efr are **not over-read** by the rule.
  - Controls:
    - C-a (CD4 ×1.5 on M4): −0.00087 [−0.00165, −0.00006]. 438's CD4 ×1 decision stands.
    - C-b (the LODO key permuted over targets): −0.0102 [−0.0116, −0.0085], 0/12 folds.
  - Surprise rule: no sign flip on M3–M6, and cumulative retention is above 0.3, so this is not a surprise.
- **For the coordinator (information only, nothing changes).**
  - The LODO key is the weakest October switch off A/B/C.
    - It is the one switch whose R6 value sits below half its tuned-on value.
    - It is negative on S4: −0.00060, 0/3, relative cosine −1.9 %.
  - The percentile cap carries most of R6's October value: +2.62 % relative cosine on A/B/C, +18.9 % on R6. It is
    ≈ 0 on S4, the panel with the weakest donor transfer (ret 0.21–0.30).
- **Deviation.** The owner set the gate floor to 27 GiB (`work/partB/GATE-OVERRIDE.json`, `swapguard534.sh`), because
  §9's 33 GiB free+inactive was unreachable (25.3–29.3 GiB) with the desktop apps open. No design, seed or reading rule
  moved.
- **C10b remains:** one more access, after the 10-05 read, for the final October recipe's exact bulk configuration.
  The four R6 operators per panel were deleted after pricing to free ~19 GB. They are deterministic re-fits (seed 0), and
  their sha256 values and exact commands are in `work/partB/state534.json` for C10b.

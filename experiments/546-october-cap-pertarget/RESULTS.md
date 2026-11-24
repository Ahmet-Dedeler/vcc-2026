# 546 — the October bulk candidate as a whole (LODO key + percentile cap + per-target sharpening): positive on 9/9 October-like folds

2026-09-27 05:12Z, third session. Exploratory, `comp546o.py`.
- lib438, 438's October-like panels, the October stack h2c1, the fold's own LODO fit, q 1.75 base, realised 16.
- Known answers: the q-1.5 cap-3 row == 438's `h2c1|V1_5_3+sh1.5`, and fixed cap 6 == 540o's row, |d| 0 on 9/9.

## Result (relative cos, mean of three cancer folds, folds positive; Δpds literal vs the percentile-cap row)
| panel | percentile cap (p96.647) vs cap 3 (the current October default) | + P2g on top (carried-only terciles, 493h) | + S_P2 on top | total, cap + S_P2 vs cap 3 |
|---|---|---|---|---|
| efp | +19.43 % (3/3) | +0.60 % (3/3), Δpds +0.0068 | +0.71 % (3/3), Δpds +0.0056 | +20.27 % |
| efr | +11.58 % (3/3) | +0.83 % (3/3), Δpds +0.0041 | **+2.56 % (3/3)**, Δpds +0.0062 | +14.44 % |
| real | +14.06 % (3/3) | +0.38 % (3/3), Δpds +0.0103 | +0.24 % (3/3), Δpds +0.0110 | +14.33 % |

- For comparison, fixed cap 6 gives +12.7 / +5.1 / +8.2 % (the percentile cap adds ~6 points; 543).
- The per-target arms still add on top of the heavier cap on every fold. Their direction increment is smaller, since
  the cap already moved amplitude onto strong, consistent rows. Their pds gain is intact: literal +0.004 … +0.011 on
  PLOLO, ≈ +0.0015 … +0.004 raw on the board at k 0.34.
- S_P2 and P2g are close on efp and real, and S_P2 is clearly ahead on efr. That matches 493g: S_P2 is the robust
  October partner, with no ties to handle.

**October bulk candidate for the coordinator (after 479 and 495 read on A/B/C):** `DEF_DAY_LODO=1` with
`def_lodo_key.py`'s cap at the 96.647th percentile of carried key / median (= cap 6 on A/B/C, so 495 calibrates it),
plus the per-target arm that wins 479 (S_P2 by default).

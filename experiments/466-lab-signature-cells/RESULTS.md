# 466 — where a perturbation programme's lab signature lives: with single cells, by target class. Analysis only.

2026-09-26 15:52Z → 16:05Z (`chain466.sh`, after 465). Parts: HCT116 4 min (30.9 GB peak), HEK293T 5 min (32.5 GB),
K562 3 min (4.6 GB), summary.

## Design
- Lines: Xaira's Orion HCT116 and HEK293T (one lab, genome-wide CRISPRi, single cells from parquet) and Weissman's K562
  genome-wide screen (scPerturb single cells). 7,976 genes all three measure.
- Targets: 1,500 non-essential (outside the Replogle K562/RPE1 essential sets) and 500 essential. Each has ≥ 50
  guide-passing cells in every line and was sampled with seed 466 from 7,797 candidates.
- Every cell of a target and of the controls goes to half A or B at random. Programme = per-gene P_down over the class's
  targets, each half against its own control half (disjoint: no shared noise).
- Readings: within-line (A vs B); same lab, other line (HCT116 vs HEK293T, A vs B both ways); cross-lab (each Xaira line
  vs K562).

## Result (Pearson r of P_down − 0.5)
| targets | genes | within-line (HCT116 / HEK293T / K562) | same lab, other line | other lab (HCT116×K562 / HEK293T×K562) |
|---|---:|---|---:|---|
| non-essential, ≥ 5 CPM | 7,484 | 0.127 / 0.512 / 0.140 | **0.182** | 0.059 / 0.094 |
| non-essential, ≥ 50 CPM | 2,378 | 0.105 / 0.528 / 0.148 | **0.165** | 0.050 / 0.076 |
| essential, ≥ 5 CPM | 7,484 | 0.630 / 0.660 / 0.826 | 0.346 | **0.420 / 0.406** |
| essential, ≥ 50 CPM | 2,378 | 0.721 / 0.775 / 0.900 | 0.428 | **0.515 / 0.484** |

## Reading
- **Non-essential (weak) responses: the programme is the lab's.** Two lines from one lab agree 2–3× more than lines
  across labs (0.17–0.18 vs 0.05–0.09).
- **Essential (strong) responses: the programme is biology.** Other labs' lines agree at least as well as the same lab's
  (0.41–0.52 vs 0.35–0.43). This confirms 462i on a different design, with cells.
- The smoke test's lesson holds: at 1–2 cells per target, "down" tracks expression level alike in every line and
  inflates every r. The ≥ 50-cell and ≥ 50-CPM floors are what make this readable. The class contrast is the same at
  both CPM floors.
- For October: Arc's panel is non-essential, the regime where a same-lab programme can carry across lines. This
  supports 458d's premise; it does not replace the board read.
- HEK293T's within-line reliability on non-essential targets (0.51) is far above HCT116's and K562's (0.13–0.15).
  Its screen has more cells per target (223 batches).

| step | Δ pooled real/efp/efr (k = 1, board units) | 95 % CI | folds positive (of 9) | verdict | Δ on `abc` [95 % CI] | the board's A/B/C read |
|---|---:|---|---|---|---|---|
| normrel: normrel -1.5 with its cap (per-target donor weights vs equal weights) (B0 − B2) | **+0.00112** | [+0.00041, +0.00182] | 8/9 | **HOLDS** | +0.00116 [+0.00046, +0.00193] (3/3) | +0.00257 (312 +0.0016245 (off -> plain -1.5), then 319 +0.0009408 (plain -> cap)) |
| cap: the cap form of normrel (cap vs plain at -1.5) (B0 − B1) | **-0.00011** | [-0.00049, +0.00028] | 3/9 | **FAILS** | +0.00028 [-0.00017, +0.00074] (3/3) | +0.00094 (319 +0.0009408) |
| partial: the CD4 partial rows at g 0.5 (vs none) (B0 − B3) | **-0.00026** | [-0.00083, +0.00026] | 4/9 | **FAILS** | +0.00094 [+0.00017, +0.00173] (3/3) | +0.00146 (319p +0.0011593 (g 0.34), then 319p2 +0.0002969 (g 0.34 -> 0.5)) |

| step | real | efp | efr | abc | pooled, unlevelled |
|---|---:|---:|---:|---:|---:|
| normrel | +0.00038 (3/3) | +0.00154 (3/3) | +0.00145 (2/3) | +0.00116 (3/3) | +0.00112 |
| cap | -0.00007 (1/3) | +0.00023 (2/3) | -0.00049 (0/3) | +0.00028 (3/3) | -0.00005 |
| partial | -0.00033 (1/3) | +0.00013 (2/3) | -0.00059 (1/3) | +0.00094 (3/3) | -0.00021 |

## POST-HOC — the same map on each pair's own shared genes

- 315 cross-line pairs; shared genes per pair: median 9517, min 4185, max 18533 (the frozen map uses 2,994 for every pair)
- `C_raw` wide: mean +0.0201, median +0.0143, share <= 0 8.3 %; rho with the frozen `C_raw` **+0.811**

| predictor (wide axis) | rho with `pds`, pooled | frozen-axis value |
|---|---|---|
| `C_lolo` (held line out) | **+0.457** [+0.401, +0.510] (n 850) | +0.387 [+0.326, +0.444] |
| baseline: donor-donor agreement | **+0.396** [+0.331, +0.456] (n 777) | +0.307 [+0.239, +0.372] |
| `C_out` (held line out AND every incumbent source) | **+0.253** [+0.184, +0.316] (n 713) | +0.198 [+0.127, +0.263] |

paired dRho (`C_lolo` - donor-donor), wide axis: **+0.052** [-0.011, +0.117] on 777 rows (frozen axis +0.066 [+0.004, +0.133])

## POST-HOC — the prediction on two disjoint halves of the panel

| half | n rows | `C_lolo` (frozen axis) | `C_lolo` (wide axis) | donor-donor (wide) |
|---|---:|---|---|---|
| 0 | 407 | +0.392 [+0.306, +0.474] | +0.454 [+0.374, +0.525] | +0.367 [+0.273, +0.455] |
| 1 | 443 | +0.380 [+0.295, +0.453] | +0.458 [+0.374, +0.530] | +0.417 [+0.325, +0.494] |

peak RSS 1.20 GiB

# 429 table (H1, production bulk `--bulk ps`, seeds [0, 1, 2])

## Absolute (seed mean)

| arm | pds | mse | nmae | fid | reach | jac |
|---|---:|---:|---:|---:|---:|---:|
| base | 0.89486 | 0.92529 | 0.95232 | 0.50920 | 0.27382 | 0.07480 |
| base + zs | 0.89486 | 0.92529 | 0.95000 | 0.50920 | 0.27382 | 0.07480 |
| MDE | 0.89486 | 0.92529 | 0.94272 | 0.50892 | 0.27077 | 0.07465 |
| MDE + zs | 0.89486 | 0.92529 | 0.94093 | 0.50892 | 0.27077 | 0.07465 |
| MDE + repair | 0.89486 | 0.92529 | 0.94272 | 0.50892 | 0.27077 | 0.07465 |
| MDE + repair + zs | 0.89486 | 0.92529 | 0.94093 | 0.50892 | 0.27077 | 0.07465 |

## Δ raw (mean ± sd over seeds), overall = Σ slope·Δ / 6

| contrast | Δpds | Δmse | Δnmae | Δfid | Δreach | Δjac | overall |
|---|---:|---:|---:|---:|---:|---:|---:|
| base + zs − base (zs alone) | 0 | 0 | -0.00232 ± 0.00005 | 0 | 0 | 0 | **+0.00064 ± 0.00001** |
| MDE − base | 0 | 0 | -0.00960 ± 0.00022 | -0.00028 ± 0.00004 | -0.00305 ± 0.00220 | -0.00014 ± 0.00001 | **+0.00185 ± 0.00039** |
| MDE + zs − base (the stack) | 0 | 0 | -0.01139 ± 0.00020 | -0.00028 ± 0.00004 | -0.00305 ± 0.00220 | -0.00014 ± 0.00001 | **+0.00234 ± 0.00039** |
| MDE + repair + zs − base (the stacked arm) | 0 | 0 | -0.01139 ± 0.00020 | -0.00028 ± 0.00004 | -0.00305 ± 0.00220 | -0.00014 ± 0.00001 | **+0.00234 ± 0.00039** |
| MDE + zs − MDE (zs on top of the MDE) | 0 | 0 | -0.00179 ± 0.00007 | 0 | 0 | 0 | **+0.00049 ± 0.00002** |
| MDE + repair − MDE (repair on H1) | 0 | 0 | +0.00000 ± 0.00000 | 0 | 0 | 0 | **-0.00000 ± 0.00000** |
| MDE + repair + zs − MDE + zs (repair on H1, with zs) | 0 | 0 | +0.00000 ± 0.00000 | 0 | 0 | 0 | **-0.00000 ± 0.00000** |

nmae interaction (MDE+zs) − MDE − zs: +0.00053 (per seed +0.00055, +0.00055, +0.00050); zs keeps 0.77 of its stand-alone nmae gain on top of the MDE

## nmae paired target bootstrap (seed-averaged per-target Δ, 20,000 resamples)

| contrast | targets | mean | 95% CI | p(≥0) | better / worse | worse by > 0.005 | top-3 share |
|---|---:|---:|---|---:|---|---:|---:|
| base + zs − base (zs alone) | 58 | -0.00232 | [-0.00433, -0.00073] | 0.0006 | 32 / 21 | 0 | 0.60 |
| MDE − base | 58 | -0.00960 | [-0.01350, -0.00577] | 0.0000 | 45 / 12 | 9 | 0.23 |
| MDE + zs − base (the stack) | 58 | -0.01139 | [-0.01587, -0.00704] | 0.0000 | 46 / 12 | 8 | 0.25 |
| MDE + repair + zs − base (the stacked arm) | 58 | -0.01139 | [-0.01590, -0.00711] | 0.0000 | 46 / 12 | 8 | 0.25 |
| MDE + zs − MDE (zs on top of the MDE) | 58 | -0.00179 | [-0.00338, -0.00052] | 0.0006 | 26 / 25 | 0 | 0.61 |
| MDE + repair + zs − MDE + zs (repair on H1, with zs) | 58 | +0.00000 | [+0.00000, +0.00000] | 1.0000 | 0 / 0 | 0 |  |

## zs verify census (423's verify, per target block, all seeds)

- **base + zs**: n_genes_U_changed_sum 0, n_touched_U_changed_sum 0, n_untouched_U_changed_sum 0, call_flips_sum 0, sign_flips_sum 0, sums_equal_all True, libs_equal_all True, ballast_equal_all True, rp423_equal_all True, blocks 255, n_moved_mean 123.4, counts_moved_mean 2.758e+04, frac_massw_mean 0.02056, n_genes_tie_changed_mean 82.73, p_logmax_max 0.001108, sign_guard_sum 0, sub_budget_exhausted_cells_mean 19.96
- **MDE + zs**: n_genes_U_changed_sum 0, n_touched_U_changed_sum 0, n_untouched_U_changed_sum 0, call_flips_sum 0, sign_flips_sum 0, sums_equal_all True, libs_equal_all True, ballast_equal_all True, rp423_equal_all True, blocks 255, n_moved_mean 98.66, counts_moved_mean 2.048e+04, frac_massw_mean 0.007394, n_genes_tie_changed_mean 55.66, p_logmax_max 0.001482, sign_guard_sum 0, sub_budget_exhausted_cells_mean 28.43
- **MDE + repair + zs**: n_genes_U_changed_sum 0, n_touched_U_changed_sum 0, n_untouched_U_changed_sum 0, call_flips_sum 0, sign_flips_sum 0, sums_equal_all True, libs_equal_all True, ballast_equal_all True, rp423_equal_all True, blocks 255, n_moved_mean 98.66, counts_moved_mean 2.048e+04, frac_massw_mean 0.007394, n_genes_tie_changed_mean 55.66, p_logmax_max 0.001482, sign_guard_sum 0, sub_budget_exhausted_cells_mean 28.43

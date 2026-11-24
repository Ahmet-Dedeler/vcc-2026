## Known-answer row — the cloud scorer against the frozen 095 h1-test artifact

`experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad`,
scored on 2026-08-27 by `cell-eval2` 0.16.0, re-scored today on the cloud worker.
All ten statistics of `agg_results.csv`'s `mean` row, not just `mse`.

| metric | 2026-08-27, on disk | cloud worker | rel |
|---|---:|---:|---:|
| `de_wilcoxon_direction_fidelity_yield_raw` | 0.5046621497798164 | 0.5046621497798164 | 0.00e+00 |
| `de_wilcoxon_direction_reach_raw` | 0.20045827672164832 | 0.20045827672164832 | 0.00e+00 |
| `de_wilcoxon_lfc_nmae` | 0.9667252780145132 | 0.9667252780145129 | 3.33e-16 |
| `de_wilcoxon_sig_jaccard` | 0.06359002047099205 | 0.06359002047099205 | 0.00e+00 |
| `expr_distance_unbiased` | 0.0018008600704525074 | 0.0018008600704525074 | 0.00e+00 |
| `expr_mse_unbiased` | 0.003700950826313404 | 0.003700950826313404 | 0.00e+00 |
| `expr_mse_unbiased_capped` | 0.003700950826313404 | 0.003700950826313404 | 0.00e+00 |
| `expr_mse_unbiased_capped_norm` | 2.0551018299736388 | 2.0551018299736388 | 0.00e+00 |
| `expr_real_mass_ratio` | 1.0 | 1.0 | 0.00e+00 |
| `pds_cosine` | 0.8655462184873949 | 0.8655462184873949 | 0.00e+00 |

**Worst relative difference over all ten: `3.331e-16` — PASS.**

## Gate — the board's raw deltas, derived from its scaled ones

| rung | member | Δscaled / slope | Δraw recorded | rel |
|---|---|---:|---:|---:|
| `a15`→`rc_callset` | `reach` | +0.068659 | +0.068657 | 3.2e-05 |
| `rc_callset`→`cd4f034` | `pds` | -0.003646 | -0.003650 | 1.2e-03 |
| `rc_callset`→`cd4f034` | `mse` | -0.037316 | -0.037320 | 1.0e-04 |
| `rc_callset`→`cd4f034` | `nmae` | -0.006394 | -0.006390 | 6.3e-04 |
| `rc_callset`→`cd4f034` | `fid` | +0.002717 | +0.002720 | 1.1e-03 |
| `rc_callset`→`cd4f034` | `reach` | +0.023381 | +0.023380 | 4.9e-05 |
| `rc_callset`→`cd4f034` | `jac` | +0.000394 | +0.000390 | 9.2e-03 |

Worst relative disagreement **9.2e-03** over 7 recorded raw deltas.

## h1-test, real `cell-eval2` 0.16.0, RAW `agg_results.csv` means

| candidate | pds | mse | nmae | fid | reach | jac |
|---|---:|---:|---:|---:|---:|---:|
| `ka095` | 0.865546 | 2.055102 | 0.966725 | 0.504662 | 0.200458 | 0.063590 |
| `a15` | 0.873810 | 1.015978 | 0.960014 | 0.507478 | 0.166269 | 0.064039 |
| `a15_seed1` | 0.874510 | 1.016016 | 0.960063 | 0.508254 | 0.162949 | 0.064143 |
| `rc_hybrid` | 0.873810 | 1.015978 | 0.960046 | 0.508602 | 0.229541 | 0.064304 |
| `rc_callset` | 0.873810 | 1.015978 | 0.960046 | 0.495058 | 0.226704 | 0.074204 |
| `cd4f034` | 0.873249 | 1.018538 | 0.958995 | 0.495303 | 0.225905 | 0.074199 |

### The h1-test noise floor — `a15` against its re-seeded twin `a15_seed1`

Identical config and identical effect matrix; `library_pick_seed` 0 → 1, which also re-draws each target's 400 control library sizes, so this is the seed noise on all six members and not only on the DE ones.

| | pds | mse | nmae | fid | reach | jac |
|---|---:|---:|---:|---:|---:|---:|
| |Δ| raw seed noise | 0.000700 | 0.000039 | 0.000048 | 0.000777 | 0.003319 | 0.000104 |

## THE TRANSFER TABLE — board Δraw / h1-test Δraw, per member, per rung

| rung | member | h1-test Δraw | board Δraw | **transfer** | readable? |
|---|---|---:|---:|---:|---|
| `a15`→`rc_hybrid` | `pds` | +0.000000 | -0.000000 | — | NO — |Δ| ≤ noise 0.000700 |
| `a15`→`rc_hybrid` | `mse` | +0.000000 | +0.000000 | — | NO — |Δ| ≤ noise 0.000039 |
| `a15`→`rc_hybrid` | `nmae` | +0.000032 | -0.000031 | — | NO — |Δ| ≤ noise 0.000048 |
| `a15`→`rc_hybrid` | `fid` | +0.001124 | +0.000190 | +0.169× | yes |
| `a15`→`rc_hybrid` | `reach` | +0.063272 | +0.048654 | +0.769× | yes |
| `a15`→`rc_hybrid` | `jac` | +0.000265 | +0.000134 | +0.507× | yes |
| `rc_hybrid`→`rc_callset` | `pds` | +0.000000 | +0.000000 | — | NO — |Δ| ≤ noise 0.000700 |
| `rc_hybrid`→`rc_callset` | `mse` | +0.000000 | -0.000000 | — | NO — |Δ| ≤ noise 0.000039 |
| `rc_hybrid`→`rc_callset` | `nmae` | -0.000000 | -0.000000 | — | NO — |Δ| ≤ noise 0.000048 |
| `rc_hybrid`→`rc_callset` | `fid` | -0.013544 | +0.000146 | -0.011× | yes |
| `rc_hybrid`→`rc_callset` | `reach` | -0.002837 | +0.020005 | — | NO — |Δ| ≤ noise 0.003319 |
| `rc_hybrid`→`rc_callset` | `jac` | +0.009900 | +0.010119 | +1.022× | yes |
| `rc_callset`→`cd4f034` | `pds` | -0.000560 | -0.003646 | — | NO — |Δ| ≤ noise 0.000700 |
| `rc_callset`→`cd4f034` | `mse` | +0.002561 | -0.037316 | -14.573× | yes |
| `rc_callset`→`cd4f034` | `nmae` | -0.001052 | -0.006394 | +6.080× | yes |
| `rc_callset`→`cd4f034` | `fid` | +0.000244 | +0.002717 | — | NO — |Δ| ≤ noise 0.000777 |
| `rc_callset`→`cd4f034` | `reach` | -0.000799 | +0.023381 | — | NO — |Δ| ≤ noise 0.003319 |
| `rc_callset`→`cd4f034` | `jac` | -0.000005 | +0.000394 | — | NO — |Δ| ≤ noise 0.000104 |

## The CD4 rung on the targets the frozen bank actually covers

The `complete mean3` CD4 policy covers **4 of h1-test's 85 targets** against **276 of the challenge panel's 300**. The whole-surface delta is diluted ~19.6× by construction; these four are where the arm is even defined.

| target | Δpds | Δmse | Δnmae | Δfid | Δreach | Δjac |
|---|---:|---:|---:|---:|---:|---:|
| `ANKZF1` | +0.00000 | +nan | -0.05332 | +0.00693 | +0.00000 | -0.00039 |
| `MED15` | -0.03571 | +nan | -0.00258 | +0.01077 | -0.00154 | -0.00038 |
| `MED25` | +0.00000 | +nan | -0.00509 | +0.00231 | -0.00267 | +0.00033 |
| `ZNF32` | -0.01190 | +nan | +nan | +0.00115 | +0.00000 | +0.00000 |
| *the other 81, mean* | +0.00000 | +nan | -0.00000 | -0.00000 | -0.00079 | +0.00000 |

### The same rung, coverage-corrected

Only covered targets move (proven: 81/85 blocks digest-identical, and the *other 81* row above is zero to five decimals). A per-target mean over a whole surface is therefore the per-COVERED-target effect times that surface's coverage — **0.047 on h1-test, 0.920 on the challenge panel**. Dividing both sides out makes them comparable.

| member | h1-test, per covered target | board, per covered target | transfer |
|---|---:|---:|---:|
| `pds` | -0.01190 | -0.00396 | +0.333× |
| `nmae` | -0.02033 | -0.00695 | +0.342× |
| `fid` | +0.00529 | +0.00295 | +0.558× |
| `reach` | -0.00105 | +0.02541 | -24.117× |
| `jac` | -0.00011 | +0.00043 | -3.904× |

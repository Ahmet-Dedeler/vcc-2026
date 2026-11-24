# a15 board result — CD4 residual transfer confirmed

Entry `KpoZt5Rie0EjgMaZQECE`, submitted 2026-08-30 20:14 UTC.
Model name: `The Fern Filed An Appeal`.

**Overall `0.1740865087`**, rank 9 of 482 when read at 2026-08-30
20:50 UTC. Against a14's `0.1704242443`, this is **`+0.0036622644`**.

| member | a14 | a15 | delta |
| --- | ---: | ---: | ---: |
| `score_pds` | 0.738614 | 0.745871 | **+0.007256** |
| `score_mse` | 0.110615 | 0.120057 | **+0.009441** |
| `score_nmae` | 0.091972 | 0.094239 | **+0.002266** |
| `score_fid` | -0.032586 | -0.030397 | **+0.002190** |
| `score_reach` | 0.114542 | 0.115368 | **+0.000826** |
| `score_jac` | -0.000612 | -0.000619 | -0.000007 |

Raw metrics:

| member | a14 | a15 | direction |
| --- | ---: | ---: | --- |
| `pds_cosine` | 0.832702 | 0.835987 | improved |
| `expr_mse_unbiased_capped_norm` | 0.882775 | 0.873805 | improved |
| `de_wilcoxon_lfc_nmae` | 0.945764 | 0.944251 | improved |
| `de_wilcoxon_direction_fidelity_yield_raw` | 0.502881 | 0.503536 | improved |
| `de_wilcoxon_direction_reach_raw` | 0.180615 | 0.181173 | improved |
| `de_wilcoxon_sig_jaccard` | 0.030081 | 0.030080 | effectively flat |

## Frozen falsifiers

1. Raw nmae rises by more than `0.004`: **PASS** — it fell by
   `0.001513`.
2. Scaled PDS falls at all: **PASS** — it rose by `0.007256`.
3. Raw mse improves by less than `0.008`: **PASS** — it improved by
   `0.008971`.
4. Challenge-control `n_pred` changes by more than 5%: **PASS** before
   submission — measured `-0.45%`.

All six scaled members moved favorably except Jaccard, whose change was only
`-0.0000066`. This is a successful independent board validation of the
CD4-derived target-specific residual at frozen effective fraction
`f=0.09997`.

The submitted package SHA-256 is
`55f7bd41bfc4d751bcc125d42e25b40c7aed474589257c6db4b575afe26aac59`.

# a14 board result — the per-context cliff was real, with one directional falsifier

Entry `0NvkHnMPzLFaJl8phAxS`, submitted 2026-08-30 00:00 UTC.
Model name: `The Houseplant Has Notes`.

**Overall `0.1704242443`** against a12's `0.16798741`: **`+0.00243683`**.
It was rank 13 of 449 when read.

| member | a12 | a14 | delta |
| --- | ---: | ---: | ---: |
| `score_pds` | 0.737885 | 0.738614 | +0.000729 |
| `score_mse` | 0.091163 | 0.110615 | **+0.019452** |
| `score_nmae` | 0.093974 | 0.091972 | **-0.002002** |
| `score_fid` | -0.031473 | -0.032586 | -0.001113 |
| `score_reach` | 0.116790 | 0.114542 | -0.002248 |
| `score_jac` | -0.000413 | -0.000612 | -0.000199 |

Raw `mse` moved `0.901349 -> 0.882775`, very close to the registered
`0.8813` prediction and inside the `[0.86, 0.90]` gate.

## Falsifiers

1. Raw `mse` outside `[0.86, 0.90]`: **PASS**.
2. `score_mse` below a12's `0.0912`: **PASS**.
3. `nmae` or `pds` falls: **PARTIAL FAIL**. `pds` rose, but `score_nmae`
   fell by `0.002002`.

## Interpretation

The per-context cap correction is real and delivered 92.6% of its predicted
`score_mse` gain (`+0.01945` observed versus `+0.0210` predicted). The total
gain under-delivered because all four DE members moved slightly negative.

This closes the a12 amplitude correction successfully, but invalidates the
simple monotone-amplitude prediction for `nmae`. It does not by itself reject
a15's CD4 residual, whose preregistered mechanism is directional rather than
another global-amplitude change. a15 must still be judged against its own
falsifiers before using the second daily slot.

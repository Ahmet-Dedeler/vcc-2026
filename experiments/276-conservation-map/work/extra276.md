## POST-HOC 1 — is it the prediction's own self-consistency?

| conservation score | sources in the map | median pairs / target | rho with `pds` (pooled, 3 cancer folds) | rho with base·truth cosine |
|---|---:|---:|---:|---:|
| C_lolo (held line out) | 24 | 25 | +0.387 [+0.326, +0.444] | +0.388 [+0.328, +0.447] |
| C_out (held line out AND every source the incumbent uses) | 16 | 3 | +0.198 [+0.127, +0.263] | +0.225 [+0.155, +0.297] |
| C_out, no CD4 pairs only | 21 | 10 | +0.370 [+0.309, +0.430] | +0.384 [+0.323, +0.443] |

## POST-HOC 2 — the mechanism, by band of `C_out` (no incumbent source, held line out)

| `C_out` band | n | mean `pds` | share anti-aligned | median `C_out` |
|---|---:|---:|---:|---:|
| bottom decile | 72 | 0.710 | 29.2 % | -0.0186 |
| 10–30 % | 144 | 0.760 | 18.8 % | -0.0038 |
| 30–50 % | 141 | 0.837 | 13.5 % | +0.0058 |
| 50–70 % | 144 | 0.784 | 19.4 % | +0.0165 |
| 70–90 % | 141 | 0.819 | 16.3 % | +0.0337 |
| top decile | 71 | 0.914 | 4.2 % | +0.0697 |
| all | 713 | 0.802 | 17.0 % | +0.0095 |

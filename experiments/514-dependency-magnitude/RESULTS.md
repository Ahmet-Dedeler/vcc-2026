# 514 — Line-specific dependency does not predict cross-line response size on non-essential targets (NEGATIVE)

2026-09-28 12:45Z, coordinator. Analysis only, public data only (no challenge truth, no build, no `vcc`). Script
`dep514.py` (peak < 1 GB, 2 s). DepMap 24Q4 Chronos subset from 360 (`CRISPRGeneEffect_subset.csv`, CC BY 4.0).

## Question
360 recorded that no panel target is strongly essential in A or C and concluded dependency is "not a response-strength
proxy on this panel", without testing it. The idea tested here: a target that is more essential in the recipient line
than in the donors (A has 15 panel targets with Chronos below the donor mean by > 0.3, C has 3) might respond more
strongly there, which would be a per-(target, context) magnitude correction known on the day for any D/E/F line in DepMap.

## Result
Response energy E = Σ δ² over non-target genes of each donor row (the scorer's log-pseudobulk units).
| test | n | r |
|---|---:|---:|
| K562 within-line, log E vs Chronos, all targets | 9,274 | −0.354 |
| K562, non-essential (Chronos > −0.5) | 7,650 | −0.062 |
| HCT116 within-line, all | 16,537 | −0.266 |
| HCT116, non-essential | 15,050 | −0.046 |
| **HCT116 vs K562: log(E_hct/E_k562) vs ΔChronos, all shared targets** | 8,781 | −0.041 (slope −0.13) |
| same, both non-essential | 7,051 | **+0.001** |
| same, \|ΔChronos\| > 0.3 | 1,631 | −0.107 (slope −0.17) |

Dependency predicts response size only through essential genes. Across lines it explains ≤ 1 % of the variance of the
energy ratio, and nothing among non-essential targets (the challenge panels are essential-free by design, §8.1g).

## Reading
Closed at its construction: a Chronos-based per-(target, context) magnitude correction. The cross-line energy ratio is
dominated by line-independent and technical factors (the LODO key already models the donor-agreement part). Do not
re-propose without a different response-size measure (e.g. DE count rather than energy) and an essential-containing panel.

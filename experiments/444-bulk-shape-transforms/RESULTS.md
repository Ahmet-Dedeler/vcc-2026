# 444 — Is the power law the right bulk shrinker? Yes: soft/hard thresholds and hybrids lose. (2026-09-24, coordinator)

Instrument: `comp444.py`, 420's PLOLO harness on the current bulk (stack|K3), each transform swapped in for 169's
`sharpen` inside `lib420.design_and_realised` (gated row norm restored in every variant). **Known answers: pow 1.5 and
pow 1.75 through this path reproduce comp420's `stack|K3+sh1.5` / `+sh1.75` rows with |d| = 0 on all four folds.**
~40 s per fold, peak RSS 5.7 GiB.

Rationale tested: under a sparse effect prior the posterior mean is closer to soft thresholding than to a power law.

| transform (vs pow 1.5; realised KMAX 16: rel cos / Δ β-emulated pds) | k562_gwps | orion_hct116 | orion_hek293t | kolf2 |
|---|---|---|---|---|
| pow 1.75 | +3.31 % / +0.0053 | +3.14 % / +0.0096 | +3.62 % / +0.0116 | +0.80 % / +0.0018 |
| soft τ = q25 of the row's gated \|h\| | −10.4 % / −0.017 | −9.0 % / −0.022 | −10.2 % / −0.026 | −4.8 % / −0.006 |
| soft q50 | −6.4 % / −0.006 | −6.3 % / −0.011 | −6.3 % / −0.015 | −5.3 % / −0.005 |
| soft q70 | −3.8 % / +0.004 | −5.2 % / −0.001 | −3.6 % / −0.001 | −8.6 % / −0.006 |
| soft q85 | −3.2 % / +0.008 | −6.1 % / +0.010 | −2.8 % / +0.017 | −13.6 % / −0.003 |
| hard q50 | −15.0 % / −0.026 | −13.2 % / −0.028 | −14.7 % / −0.035 | −7.4 % / −0.005 |
| hard q85 | −16.8 % / −0.013 | −17.0 % / −0.022 | −15.4 % / −0.019 | −16.4 % / −0.015 |
| pow 1.25 + soft q25 | −4.1 % / −0.006 | −3.8 % / −0.009 | −4.1 % / −0.012 | −2.2 % / −0.002 |
| pow 1.25 + soft q50 | −1.8 % / +0.001 | −2.4 % / −0.002 | −1.8 % / −0.003 | −3.6 % / −0.002 |
| pow 1.5 + soft q25 | +0.7 % / +0.002 | +0.5 % / +0.002 | +0.8 % / +0.002 | −0.4 % / −0.000 |
| pow 1.5 + soft q50 | +1.4 % / +0.005 | +0.7 % / +0.006 | +1.6 % / +0.008 | −2.5 % / −0.001 |

## Verdict
The power law is the right family: every threshold loses cosine on every fold, and the best hybrid is below pow 1.75
and negative on kolf2 (the fold that predicted the board for shape arms). The next bulk shape step is the q dose itself
(q 1.75 levelled), decided by 433k's board read of the q 1.25 → 1.5 step. No builder change.

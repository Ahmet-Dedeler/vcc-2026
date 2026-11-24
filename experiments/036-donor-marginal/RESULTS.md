# 036 — What one more donor line is worth, and what whitening did to that number

**Date:** 2026-08-22 · 4-line LODO (K562-gwps, Orion-HCT116, Orion-HEK293T, H1-2025), 226
perturbations shared by all four, whitening at gamma = 0.65. One whitener fitted **per
predictor subset**, so a k-donor configuration is never handed a (k+1)-donor covariance.

## Why this needed measuring

Every calibration avenue is now closed — `mse` needs a cosine we trade away by design (033),
`jac` is 1.31x chance at gene identity (014), `nmae` has no subgroup below 50% sign accuracy so
nothing is worth silencing (035), per-perturbation alpha is worth 0.25% (012). The binding
constraint is directional accuracy, and only two things have ever moved it: **more donor lines**
and **better geometry**. Geometry is exhausted (030f). So: is a fourth line worth a day?

## Result

| donors | folds x combos | macro PDS | marginal |
|---|---|---|---|
| 1 | 12 | 0.7598 | — |
| 2 | 12 | 0.8119 | **+0.0520** |
| 3 | 4 | 0.8365 | **+0.0246** |
| 4 | — | — | **+0.0117 (extrapolated)** |

Decay ratio 0.47 — each donor is worth about half the previous one. Extrapolating one step and
applying the LODO→board ratio calibrated on v3→v6 (+0.070 LODO produced +0.088 board raw `pds`,
ratio 1.26): a fourth line is worth roughly **+0.015 board raw `pds` ≈ +0.033 scaled ≈ +0.005
overall**.

## The number that reframes the whole project

**The second donor is worth +0.0520 here. Before whitening, experiment 010 measured it at
+0.0166** — and with equal weights it was *negative*.

Whitening did not merely add its own +0.070. It made every additional donor roughly **three
times more valuable**, because averaging's only cost was the smoothing whitening removes. The
two levers are multiplicative, not additive, and that changes the value of data acquisition
retrospectively: donor lines that were not worth chasing in the pre-whitening regime are worth
chasing now.

## What it means for Feng

The one genome-scale candidate in reach is Feng et al. 2025 — 7,226 genes targeted across 34
iPSC lines, ~2M cells, on Figshare. At +0.005 overall it is worth a day, with two discounts:

- **Tissue mismatch.** It is pluripotent where the contexts are somatic, the same mismatch that
  made H1 a weak donor. H1 is one of the four lines averaged into the +0.0117 above, so the
  estimate already contains one pluripotent line — but a somatic fourth line would be worth
  more and a second pluripotent one less.
- **Weak external evidence.** Bio_NUDT runs `gwps+feng` and sits at raw `pds` 0.665 against our
  0.816, though their model string also says `alpha=2 probe`, so the comparison is not clean.

Worth doing, after the submission rather than before it.

# 004 — Real K562-gwps → H1 donor transfer

**Date:** 2026-08-21 · Donor: Replogle K562 genome-wide (9,787 perturbations × 8,248 genes).
Target: `h1-test` harness (85 perturbations, 18,080 genes), 82 covered by the donor.

## The operating point

| measurement | centered cos |
|---|---|
| K562 → H1, full 18,080-gene axis | **+0.0727** |
| K562 → H1, donor-covered genes only (7,581) | **+0.0949** |
| K562 → RPE1 (experiment 003, for comparison) | +0.128 |
| cos(donor mean response, H1 mean response) | +0.161 |

**Positive for 100% of the 82 covered perturbations.** Small, but utterly consistent — there is
no subset where transfer fails.

H1 is a harder target than RPE1: embryonic stem cells versus a somatic line. The 2026 contexts
are described as six lines "from different tissues", which are more likely somatic/immortalised
than embryonic, so **0.073 is probably a pessimistic bound and 0.128 a more representative one.**

Gene-axis dilution costs ~0.022: the donor covers 41.9% of the full axis (71% of the >5 CPM
DE-scored set), and the uncovered genes are predicted unchanged.

Adding the mean response `m` changes centered cos by +0.0002 — confirming experiment 003's
point that `m` is orthogonal to what earns score.

## Negative result 1 — low-rank denoising hurts, monotonically

Hypothesis: gwps has a median of 178 cells/perturbation, so effect vectors are noisy;
perturbation responses are believed low-dimensional; truncating the SVD should strip noise.

| rank | var explained | transfer cos | vs full |
|---|---|---|---|
| 5 | 0.124 | 0.0392 | −0.033 |
| 50 | 0.213 | 0.0443 | −0.028 |
| 200 | 0.304 | 0.0461 | −0.027 |
| 800 | 0.514 | 0.0519 | −0.021 |
| **full** | 1.000 | **0.0727** | — |

Every truncation loses. And the spectrum is genuinely flat — **800 components explain only 51%
of variance.** The transferable signal is high-dimensional and spread across the spectrum, not
concentrated in a few shared programs.

**This matters beyond the experiment.** Any architecture with a bottleneck — autoencoders,
latent-space perturbation models, low-rank factorisations — destroys exactly the signal that
transfers. It is a mechanistic explanation for the published finding that simple baselines beat
deep models on perturbation prediction.

## Negative result 2 — per-gene transfer weights hurt

Hypothesis: some genes' responses are portable (core machinery), others context-specific. If that
is a property of the *gene*, weights fitted on one donor→target pair should help another.
Fitted per gene as the least-squares slope `w_j = Σ(δ_target·δ_donor)/Σδ_donor²`.

| | cos |
|---|---|
| unweighted baseline | **0.0949** |
| w fitted on H1 train+val, eval H1 test (same line, held-out perts) | 0.0710 |
| w fitted on K562→RPE1, eval K562→H1 (unseen line) | 0.0672 |
| same weights, gene labels shuffled (null) | 0.0578 |
| binary gate, keep w > q25 / q50 / q75 | 0.0884 / 0.0858 / 0.0696 |

The weights beat their own shuffled null (0.067 vs 0.058), so they carry real information — but
they are **worse than uniform**. Cosine is a statement about direction, and per-gene rescaling
learned from noisy data distorts direction more than it corrects magnitude. Binary gating fails
the same way.

## Conclusion

Two independent attempts to improve the donor estimate — spectral truncation and per-gene
reweighting — both lose to the raw delta. **Uniform, full-rank donor transfer is the estimator to
beat**, and the remaining headroom is in calibration (experiment 002), not in transforming `s`.

Do not spend more time on transforms of the donor matrix without a mechanism that explains why
these two failed.

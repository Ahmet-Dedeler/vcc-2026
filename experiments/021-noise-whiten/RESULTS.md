# 021 — Whitening by the noise covariance loses to whitening by the signal covariance

**Date:** 2026-08-21 · LODO, 259 VCC-panel perturbations, truth = raw held-out effects.

## The theory that motivated it

Experiment 017 whitens with `C_signal = B^T B / n`, the covariance of donor **effect**
vectors. That is not what estimation theory asks for. Writing a prediction as `P = T + e`, the
optimal retrieval similarity is Mahalanobis in the **error** covariance:

```
similarity(P, T) ~ P^T Sigma_e^{-1} T
```

which coincides with signal whitening only when the error is proportional to the signal. The
measured optimum sitting at `gamma = 0.65` rather than 1.0 looked like a hint that the right
object had not been used yet.

`Sigma_e` is directly estimable and target-blind: two donors measuring the same perturbation
differ by exactly what we are calling error, so `Sigma_e ≈ cov((d_A − d_B)/2)` over the
thousands of perturbations both predictor lines carry.

## Result

| basis | best gamma | macro PDS |
|---|---|---|
| **signal** (donor effects) | 0.65 | **0.7990** |
| noise (donor half-differences) | 0.3 | 0.7548 |
| ratio (noise-whiten then signal-colour) | 0.3 | 0.7576 |
| none | — | 0.7210 |

Noise whitening helps a little at small gamma and then **plateaus flat** (0.7310, 0.7304,
0.7303 at gamma 0.65, 0.85, 1.1) — it saturates rather than peaking, which says it is
removing one small thing and then doing nothing. The signal basis is better by 0.044.

## Why, and it is worth stating

The theory assumes the error is *noise*. It is not. The difference between two donor lines is
dominated by **genuine line-specific biology**, and line-specific biology loads on the same
dominant programmes as the shared response — growth rate, ribosome content, stress tone.
Whitening by it therefore removes signal along with error.

The "ratio" arm is the sharper version of the same finding. Noise-whitening then
signal-colouring is the LDA / generalised-eigenproblem direction: up-weight directions with
high signal-to-noise. It loses monotonically. So the data says the high-variance directions
have high SNR *and* low discriminative value — they are estimated accurately and shared by
every perturbation.

**Discriminability is not signal-to-noise.** That is the reason a metric like `pds` rewards an
operator that deliberately down-weights the best-measured part of the prediction, and it is
why the naive "improve the estimate" instinct keeps failing on this task.

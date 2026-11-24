# 160 — the control cells DO carry line-specific response structure, and it is unusable

**Date:** 2026-09-02 · two held Orion lines, exact-space pseudobulk truth, deterministic
(seed `20260902`, no emitter, no `cell-eval2`).

## Verdict

**NO-GO, and the closure is mechanistic.** A line's own control-cell variation captures **8-11x**
more of its target-specific response residual than another line's does and **25-340x** more than
an envelope-matched random subspace — a real, large, context-specific signal. It cannot be
monetised: **every way of using it costs more `pds` than the cosine it buys**, because the
capture is concentrated in a single direction that is the line's *mean* response, which is
exactly the hub `pds` punishes.

This is the fifth independent route into the same wall (155 per-gene regression, 067 variance
transport, 017 basal modulation, Heisenberg's board arm adding context-shared co-expression).
Do not reopen control-cell geometry without a mechanism that separates the line's mean response
from its target-specific structure — this experiment shows the control cells do not.

## 1. The hypothesis, and why it was worth testing

Experiment 159 measured that ~70% of a target's response energy is line-specific and that
knowing the **context** is worth ~3x knowing the **target** (same-line cos 0.2627 vs best
cross-line single donor 0.0897). Production spends all its modelling on the target axis and
none on the context axis. The only context information the challenge gives is 18,400
unperturbed cells, so the geometric hypothesis: *a perturbation pushes the cell along directions
the population already varies in* — unperturbed populations contain cells at different
cell-cycle, stress and metabolic states, which are the programmes a knockdown engages. If true,
the control cells hand us a context-specific basis for free and the transfer problem collapses
from 18,533 gene directions to ~50 loadings.

## 2. Containment — strongly positive

Control cells in the scorer's per-cell space (`log1p(TS·p)`), cell-centred, top-k right singular
vectors. True effects split into shared mean `m` and target-specific residual `R = E − m`; only
`R` can earn score, so `R` is what must be contained. Fraction of `‖R‖²_F` inside the span:

| | k=5 | k=20 | k=50 | k=200 |
| --- | ---: | ---: | ---: | ---: |
| **HCT116** own control basis | **0.0653** | 0.0938 | **0.1103** | 0.1284 |
| HCT116 *other line's* control basis | 0.0081 | 0.0107 | 0.0145 | 0.0249 |
| HCT116 random, per-gene sd envelope matched | 0.0004 | 0.0018 | 0.0044 | 0.0174 |
| HCT116 oracle (top-k PCs of `R` itself) | 0.1238 | 0.2665 | 0.4978 | 1.0000 |
| **HEK293T** own control basis | **0.1392** | 0.1597 | **0.1719** | 0.1837 |
| HEK293T *other line's* basis | 0.0128 | 0.0164 | 0.0202 | 0.0302 |
| HEK293T random, envelope matched | 0.0004 | 0.0017 | 0.0042 | 0.0168 |
| HEK293T oracle | 0.1863 | 0.3116 | 0.5257 | 1.0000 |

At k=5 the own basis is **8.0x / 10.9x** the other line's and **153x / 339x** the random null.
The `other` column is the load-bearing control: it rules out generic co-expression structure.
The envelope-matched random null rules out an expression-concentration artefact.

**Confounds tested and cleared.** Regressing log library size out of the control cells before
building the basis changes capture by `0.0000 / 0.0000`; adding mitochondrial fraction,
`−0.0001 / −0.0007`. A basis built on one disjoint half of the control cells gives
`0.1045 / 0.1645` against the full `0.1091 / 0.1698`, and the two half-bases overlap at
`0.492 / 0.383` of their dimension — so the basis is reproducible and the capture is not a
technical-covariate artefact.

## 3. The payoff — negative on every variant

Projection applied to the a23 delivered effect, then mass-fixed and rendered (experiment 156:
scoring the pre-render effect measures a different object), scored against exact-space truth.

| arm | HCT116 Δpds | HCT116 Δcos | HEK293T Δpds | HEK293T Δcos |
| --- | ---: | ---: | ---: | ---: |
| `project` onto own k=50 | **−0.21096** | −2.2% | **−0.17580** | −77.8% |
| `blend w0.5` own k=50 | −0.01637 | +9.9% | −0.00931 | +3.9% |
| `amplify g1.0` own k=50 | −0.01671 | +10.2% | −0.00948 | +5.2% |
| `blend w0.5` **other** k=50 | +0.00150 | +1.7% | −0.00173 | +1.4% |
| drop PC1, `blend w0.5` k=50 | −0.01561 | +9.1% | −0.00798 | +5.8% |
| basis ⊥ panel-mean, `blend w0.5` | −0.01677 | +3.4% | −0.00948 | +7.6% |

**Hard projection is catastrophic** (−0.18 to −0.27 pds): forcing 132 predictions into a ≤200-dim
space makes them mutually similar, which is hubness by construction. **Blend and amplify buy
cosine and sell `pds`**, and at a23's measured board exchange rate (+14.6% cosine bought
+0.00174 pds and +0.00153 mse) a −0.016 pds loss is roughly 5x larger than the cosine is worth.
The `other`-basis rows are +0.0015 / −0.0017 — inconsistent in sign, i.e. noise.

## 4. Why, and it is one number

**On HEK293T, dropping PC1 alone takes residual capture from 0.1698 to 0.0525.** The
context-specific capture is essentially *one direction*, and the own basis captures the **shared
mean response** at `0.710 / 0.888` (k=50) against the residual's `0.110 / 0.172` — 5-6x better.

So the control cells' apparent context-specificity is concentrated in their top principal
direction, and that direction aligns with the line's *mean* perturbation response. That is the
shared programme, which is the metric's zero for `pds` and the hub the whitening operator exists
to suppress. Amplifying it raises cosine — because the shared programme is partly correct — and
destroys discrimination. Heisenberg's board arm is the same experiment run live: adding
context-shared co-expression moved pds `0.4509 → 0.3974`.

**The target-specific component of the response is not recoverable from unperturbed variation.**
That is the honest form of the negative, and it is a statement about biology rather than about
our pipeline.

## 5. A vacuous control, recorded per rule 10's corollary

The design included a `shuffle_targets` control — deal the residual rows to the wrong targets and
re-measure capture — intended as the sharpest test of target-specificity. It returned **exactly**
`own` to four decimals on both lines (`0.1091 / 0.1698`), because `‖R Vᵀ‖²/‖R‖²` is a sum of
independent per-row terms and is therefore **invariant to any permutation of R's rows by
construction**. It is a theorem, not evidence, and it would have looked like a passed control.
The control that can actually fail is the row-row cosine of the *projections* (§4), and the
`other`-basis column.

This is the second vacuous check in two days (the first: verifying the mass fix by testing
`Σ expm1(b + log1p(TS·p) − b) = TS`, which is an identity). Both were caught only by asking what
the check would read if the thing being tested were absent.

## Files

`subspace.py` (containment + nulls) · `transfer.py` (confounds + payoff) · `work/*.json`.

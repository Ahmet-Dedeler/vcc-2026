# Why the dual-moment emitter fails: the reach overwrite law

**Status:** exploratory post-hoc analysis, 2026-08-31. Uses only source-derived
`effects.npz` and the **already-unblinded** scored frames of terminal
experiments 149 and 150. No sealed data was opened and no gate was changed.
Nothing here licenses a build or a submission. **Any hypothesis derived from
this document must be preregistered before it is tested on a new surface.**

Reproduce with `analysis/reach_signflip_diagnostic.py`.

## The question

Experiments 148-152 carried one fixed lambda-1 dual-moment emitter across four
held surfaces. It gained `+0.009707` on H1 and `+0.012919` on HEK293T, then lost
`-0.005045` / `-0.004344` / `-0.003235` on HCT116 under uniform, target-gated and
gene-gated placement. Every arm left per-target PDS exactly unchanged and every
arm improved NMAE and Jaccard. The whole swing lives in `reach`.

The obvious reading was "the transform is context-dependent." That reading is
wrong, and the truth is more useful.

## Finding 1 — the emitter overwrites the incumbent's ranking

Regress per-target `reach_delta` on the incumbent's own per-target reach:

| context | pearson | fit | break-even |
| --- | ---: | --- | ---: |
| HEK293T | `-0.731` | `delta = -0.9632 * inc + 0.01663` | `+0.0173` |
| HCT116 | `-0.575` | `delta = -0.4784 * inc + 0.00781` | `+0.0163` |

Equivalently, regressing the **candidate's** reach on the incumbent's gives a
slope of `+0.037` on HEK — the candidate's reach is very nearly independent of
what the incumbent achieved. The emitter does not refine the incumbent's
directional ranking; it **replaces** it with a ranking of its own roughly fixed
quality.

The two break-even points agree to within 6% across two unrelated cell lines.
This is one law, not two behaviours.

## Finding 2 — the apparent context-dependence is a distribution difference

Both contexts show the same monotone quartile pattern:

| incumbent-reach quartile | HEK delta | HCT delta |
| --- | ---: | ---: |
| Q1 (lowest) | `+0.0130` | `+0.0192` |
| Q2 | `+0.0258` | `+0.0077` |
| Q3 | `+0.0076` | `-0.0052` |
| Q4 (highest) | `-0.0317` | `-0.1835` |

What differs is the input, not the response:

| | HEK293T | HCT116 |
| --- | ---: | ---: |
| incumbent reach, mean | `0.0134` | `0.1008` |
| incumbent reach, median | `0.0000` | `0.0185` |
| targets at exactly zero | `68/132` | `40/132` |

HEK's incumbent barely reaches at all, so replacement is an upgrade. HCT's
incumbent reaches well, so replacement is destruction. Same law, different
starting point.

## Finding 3 — the gate that would fix it is unobtainable

Applying the transform only where the incumbent reaches shallowly would work.
That requires knowing per-target incumbent reach on the **target** context.

It does not transfer:

- `spearman(inc_reach_HEK, inc_reach_HCT) = -0.033` over 120 shared targets
- `P(HCT reach = 0 | HEK reach = 0) = 0.290`, base rate `0.292` — **zero
  information**

And no source-side feature predicts it either. Across seven features
(effect norm, relaxation norm and relative magnitude, incumbent-direction
cosine, relaxation-incumbent cosine, sign-agreement fraction, top-1% mass
concentration) the largest `|rho|` against incumbent reach was `0.16` on HEK and
`0.07` on HCT; against `reach_delta`, `0.11`.

Gating HCT116 on HEK293T's incumbent reach — the honest cross-context version —
moves the reach delta from `-0.0377` to `-0.0277` at its best threshold. Still
deeply negative.

## Finding 4 — even a perfect oracle does not clear the bar

A greedy per-target oracle that knows the held HCT116 outcomes and applies the
transform only where it improves the projected average:

```
incumbent projected avg              +0.154002
uniform candidate                    +0.150051   (-0.003951)
PERFECT ORACLE (59/120 targets)      +0.161209   (+0.007207)

experiment 152's gate required                    >= +0.006000
```

The entire family had at most `+0.0072` available, and only to a selector with
perfect foresight. Experiment 152's preregistered gate sat at 83% of a
theoretical maximum that requires the held answer.

## Conclusion

The dual-moment relaxation direction is closed, and now for a demonstrated
reason rather than an accumulation of failures:

1. it replaces rather than refines the incumbent's directional ranking;
2. whether that helps is decided by the incumbent's own per-target reach;
3. that quantity does not transfer across contexts and is not predicted by any
   source feature tested;
4. the perfect-information ceiling is `+0.0072`, barely above the gate that was
   already missed.

**Do not** build further variants that re-weight, re-place, or re-dose this
transform. A useful successor would have to *combine* the incumbent's ranking
with the emitter's rather than substitute one for the other — and it would need
a build-time confidence signal that this analysis has not found.

## Caveats

The regressions are post-hoc on 120-132 targets per context, and the two
contexts' fits differ in slope (`-0.963` vs `-0.478`) even though their
break-even points agree. The oracle is greedy and coordinate-wise, so it is a
lower bound on the true optimum, not an exact ceiling. `reach` is heavily
zero-inflated (68 and 40 targets at exactly zero), which inflates Pearson
relative to Spearman — both are reported for that reason.

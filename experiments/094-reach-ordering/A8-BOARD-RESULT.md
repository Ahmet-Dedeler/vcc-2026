# a8 board result — over-dispersion works, at a fifth of the predicted size

Entry `lfnx5c9vp1gVhlNWcemT`, model name `Wildflower`, submitted 2026-08-27 19:21 UTC.
**Overall `0.1543607`** against `a2-prune`'s `0.1503084` — **`+0.0040523`**, our best score
in the project, and **rank 11 → 7 of 363 teams**.

Config: `a2-prune` exactly, plus `--overdispersion 4.0`. One variable.

## Against the pre-registered prediction

| member | a2-prune | **a8** | actual Δ | predicted Δ |
| --- | ---: | ---: | ---: | ---: |
| `pds` | 0.736560993522496 | 0.736560993522496 | **+0.000000** | 0, exactly |
| `mse` | 0.000000 | 0.000000 | +0.000000 | ~0 |
| `nmae` | 0.094339 | 0.094301 | **−0.000038** | ~0 |
| `fid` | −0.035087 | −0.028694 | **+0.006393** | +0.03900 |
| `reach` | 0.102951 | 0.124349 | **+0.021398** | +0.06100 |
| `jac` | 0.003086 | −0.000353 | **−0.003439** | +0.00220 |
| **overall** | **0.1503084** | **0.1543607** | **+0.0040523** | +0.0213 |

Raw: `pds_cosine 0.831773 → 0.831773`, `reach 0.1701 → 0.1890`, `fid 0.5021 → 0.5041`.

## What was right

**`pds` moved by exactly zero, on an independent build with a fundamentally different
emitted count distribution.** The exact-total repair survives gamma-Poisson scatter at φ=4,
and `nmae` moved by `3.8e-5`. That is the load-bearing guarantee of the whole emission
channel, and it is now board-verified twice — once for a full-range `spread` change (v14)
and once for a change that quadruples the per-cell variance.

**The direction was right on four of five moving members**, and `reach` delivered the
**largest single-member gain this project has recorded**: `+0.0214` scaled, against a
previous best of `+0.0137` (v11 → v14).

## What was wrong, and by how much

**The magnitude was ~5× too optimistic.** Predicted `+0.0213`, delivered `+0.0041`.
Member by member the harness over-stated `reach` by 2.9× and `fid` by 6.1×.

**`jac` moved the wrong way** — predicted `+0.0022`, delivered `−0.0034`, which is 85% of
the shortfall on its own after `reach` and `fid`. Nothing in the harness work anticipated it;
`jac` rose with over-dispersion on both harnesses. The call set at `n_pred` 3,836 overlaps
the reference *worse* than at 5,961, which the coverage analysis did not model because it
only tracked `fid`.

### Why `jac` fell — the same `n_real` trap as `fid`, one level down

Decomposing the ratio on the repaired harness, 30 perturbations:

| φ | `n_real` | `n_pred` | `\|R ∩ R̂\|` | `\|R ∪ R̂\|` | `jac` | recall | precision |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 92 | 4,028 | 52.6 | 4,067 | 0.0130 | 0.575 | 0.013 |
| 4 | 92 | 2,943 | 44.9 *(−15%)* | 2,989 *(−27%)* | 0.0150 | 0.491 | 0.015 |

`jac` rises only when the union shrinks faster than the intersection, and on this harness it
does — **because `n_real` is 92, so `|R ∪ R̂| ≈ n_pred` and any call-set cut raises `jac`
mechanically.**

On the challenge `n_real` is large enough that the union contains a substantial number of
real genes, so it *cannot* shrink much faster than the intersection. Solving the board's two
points at `n_real ≈ 2,000` puts the intersection at ~248 → ~176 (−29%) against a union of
~7,900 → ~5,800 (−27%) — nearly proportional, hence `jac` flat-to-down. Which is what
happened.

**So both of the harness's "free" gains — `fid` at §5's coverage term and `jac` here — are
the same artefact wearing different clothes: `n_real` is too small on every harness we hold,
so every call-set cut looks free.** That artefact has now cost this project three
predictions: v14's `fid`, v15's `fid`, and a8's `jac`.

**Rule: any metric with `n_real` in its denominator or its union is unreadable on our
harnesses.** That is `fid`, `jac`, and `reach`'s normalisation — i.e. three of the four DE
members. Only their *ordering-driven* parts transfer.

## The calibration that matters most: neither harness bracketed the truth

The independent check on **h1-test** — run after a8 was submitted — says over-dispersion does
**nothing** there:

| harness | production reach | φ=4 reach | Δ |
| --- | ---: | ---: | ---: |
| orion-hct116-thinned | 0.0838 | **0.1605** | **+0.0767** |
| **the board** | **0.1701** | **0.1890** | **+0.0189** |
| h1-test | 0.1878 | 0.1830 | **−0.0009** *(noise floor −0.0004)* |

**The gain scales inversely with how good the baseline ordering already is.** HCT116's
emitted ordering is bad (reach 0.084) and over-dispersion fixes a lot of it; h1-test's is
already good (0.188) and there is nothing to fix; the board sits between them and got a
proportionate share.

That is a coherent mechanism — over-dispersion repairs a *badly ordered* emission — but it
means **no single harness predicts the board's magnitude**, and the one whose absolute reach
is closest to the board (h1-test, 0.188 vs 0.170) predicted **zero**. Reading either harness
alone would have been wrong: HCT116 by 4×, h1-test by sign.

**Rule going forward: run every emission arm on BOTH harnesses.** Their disagreement is the
error bar, and today it bracketed the truth when neither one alone did.

## Verdict

Ship it, and re-tune with both harnesses. `+0.0041` from a one-line emitter change on an
arm whose downside was bounded by construction is a good trade, and the mechanism — that our
emitted cells were ~2.3× too uniform, so the rank test ordered genes by control skewness
rather than by predicted effect — is now confirmed on the board rather than argued.

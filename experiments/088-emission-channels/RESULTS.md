# 088 results — the variance channel, and mse settled

Date: **2026-08-23**. 69 arms on h1-test and Orion HCT116. Production operator
unless stated: `g=0.65, q=1.0, a=1.08, cpm=5.0`, union axis, challenge panel
excluded.

## 1. mse is reachable, and catastrophically unprofitable. Closed.

This is the headline. Turning whitening **off** produces the first nonzero
scaled `mse` in this project's history:

| config | raw mse | scaled mse | pds | avg6 |
| --- | ---: | ---: | ---: | ---: |
| `g=0.00, a=0.392` | `0.9570` | **`+0.03363`** | `0.4189` | `0.0890` |
| `g=0.00, a=0.250` | `0.9815` | `+0.00789` | `0.4187` | `0.0871` |
| `g=0.65, a=1.08` (production) | `1.9678` | `0` | `0.8063` | **`0.1892`** |
| `g=0.95`, any alpha tested | `≥0.9916` | `0` | `0.63-0.67` | `0.119-0.121` |

`mse` costs **`-0.387` pds to buy `+0.034` mse** — an 11:1 trade against us, and
`avg6` more than halves. No intermediate gamma scores it at any alpha tested.

This retires the question three ways at once. `STRATEGY.md`'s "structurally
unreachable" was wrong — it is reachable. My own retracted cosine derivation
(`92a9c3e`) is moot — the normalization never needed to be known. And "we have
no method" is wrong — we have one, it is simply a terrible trade. **Do not spend
further compute on mse.**

## 2. The alpha response is a cliff, not a curve

At `g=0.65`, between `a=0.9` and `a=1.08`:

| alpha | reach | jac | fid | nmae |
| ---: | ---: | ---: | ---: | ---: |
| 0.90 | **`-0.0643`** | `0.0124` | `+0.0489` | `0.0763` |
| 1.08 | **`+0.1845`** | `0.0938` | `-0.0243` | `0.0752` |

reach moves `0.25`, jac jumps 8x, fid flips sign, across an 18% alpha change.
That is a threshold effect, almost certainly the `min_cpm=5.0` gate with
`gate_unexpressed=True`: below some magnitude, predicted effects stop clearing
the gate and genes leave the called set entirely.

**Alpha has never been swept above `1.08` on this operator.** The only grid in
the repo with `a>1.08` is exp031 at `g=0.85` on the *inter* axis, which predates
the union operator. And alpha was tuned for `pds`, which is scale-invariant — so
it was effectively chosen by the one metric that cannot see it, while the four
DE metrics can. Experiment 088b tests `a ∈ {1.2, 1.4, 1.7, 2.0}` plus
`cpm ∈ {2, 10}` at fixed alpha as a mechanism check.

## 3. Variance channel: positive, but not yet attributable

Best arm `k=1000, spread=0.15`: **`+0.00075` avg6** over the shipped
`k=0, spread=0.25`, `t=2.09`, all three seeds positive
(`+0.00062 / +0.00105 / +0.00057`).

Gates 1, 2 and 4 pass. `pds` holds at `0.806-0.809` and `nmae` at `0.075` across
all 48 arms, so mean-preservation is intact and there is no implementation bug.
Gate 3 passes on a reading: five of six neighbours are positive, and the one
that is not (`spread=0.10`) is the systematically worst column at every k.

**The design has a hole.** The baseline is `k=0` at `spread=0.25`; the winner is
`k=1000` at `spread=0.15`. Those differ in two ways — confidence weighting *and*
less noise — and there is no `k=0, spread=0.15` control, so the gain cannot be
attributed. The decomposition makes this worse rather than better:

| metric | baseline | best arm | delta |
| --- | ---: | ---: | ---: |
| jac | `0.0939` | `0.0987` | **`+0.0048`** |
| reach | `0.1828` | `0.1848` | `+0.0020` |
| fid | `-0.0245` | `-0.0278` | `-0.0033` |

The gain is mostly **jac**, and less noise trivially raises jac by leaving more
genes above the significance cut. That needs no confidence weighting at all. The
missing control is running.

Note also that this is **5x smaller than exp027 advertised** (`+0.004` overall),
and h1-test overstated reach by 4.3x on v9, so the board delta may be smaller
again.

## 4. basal_beta refutes its own hypothesis. Closed.

`run.py:113`'s comment proposed it as the untested `nmae` lever. It moves nmae
the wrong way:

| b | avg6 | nmae | pds |
| ---: | ---: | ---: | ---: |
| 0.05 | `0.18932` | `0.0708` | `0.8119` |
| 0.10 | `0.18939` | `0.0652` | `0.8153` |
| 0.25 | `0.18688` | **`0.0402`** | **`0.8267`** |

It is a pds-for-nmae trade. We are already first on pds, so it is the wrong
direction for us at every setting tested.

## 5. h1-test is a valid instrument for reach

Raw reach `0.23-0.24` against the board's `0.169` and HCT116's `0.015`. The risk
that h1-test was as reach-blind as HCT116 is retired; the arms above are
measured on a surface that can actually see the metric.

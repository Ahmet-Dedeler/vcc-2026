# 007 — The emission model is worth more than the transfer model

**Date:** 2026-08-21 · Harness `h1-fast`, real K562-gwps donor transfer, official anchors.

## Result

| α | emitter | pds | mse | fid | reach | jac | nmae | **est** |
|---|---|---|---|---|---|---|---|---|
| 0.7 | **even** | **0.756** | 0 | −0.003 | 0.076 | **0.127** | 0.054 | **0.168** |
| 0.85 | even | 0.754 | 0 | −0.000 | 0.058 | 0.128 | 0.042 | 0.164 |
| 1.0 | even | 0.753 | 0 | 0.002 | 0.062 | 0.129 | 0.027 | 0.162 |
| 0.4 | even | 0.756 | 0 | −0.008 | 0.057 | 0.123 | 0.044 | 0.162 |
| 1.0 | nb | 0.701 | 0 | −0.027 | 0.126 | 0.074 | −0.107 | 0.128 |
| 0.7 | nb | 0.694 | 0 | −0.160 | 0.100 | 0.051 | 0.052 | 0.123 |
| 0.4 | nb | 0.708 | 0 | −0.453 | 0.081 | 0.011 | 0.062 | 0.068 |

Effect of switching emitter at matched α:

| α | NB → even | Δest | fid |
|---|---|---|---|
| 0.4 | 0.068 → 0.162 | **+0.094** | −0.453 → −0.008 |
| 0.7 | 0.123 → 0.168 | **+0.045** | −0.160 → −0.003 |
| 1.0 | 0.128 → 0.162 | **+0.035** | −0.027 → +0.002 |

**Live board for comparison: Vivai 0.1341, IsoCell 0.1207, Sidechain 0.0788.**

## Why it works — the decoupling

Experiment 005 exposed a conflict that one parameter cannot resolve:

- `fid = k / max(n_pred, n_real)` needs **many** genes to clear the Wilcoxon/BH bar, which with
  realistic NB noise needs α ≥ 1.0
- `nmae = Σ|lfc_hat − lfc| / Σ|lfc|` punishes **overshooting** magnitude, and wants α ≈ 0.7

At α = 1.5 `nmae` collapses to −0.448; at 2.5 to −1.211. At α = 0.4 `fid` collapses to −0.453.

The resolution is that **significance and magnitude are separate functionals of the same counts**:
`nmae` reads the per-cell CPM *mean*, which α sets; significance reads the Wilcoxon on *ranks*,
which depends on mean **and** cell-to-cell variance. The even emitter drives that variance to the
integer floor (within-perturbation library cv **0.0000**, per-gene sd 0.285 vs 1.34 for NB), so
the DE test reaches full power at whatever α `nmae` wants.

Both metrics are now near their best simultaneously: `fid` −0.003 (raw 0.513 against a 0.5135
baseline — exactly at parity) and `nmae` +0.054.

## The unexpected win: `jac`

`de_wilcoxon_sig_jaccard` went **0.051 → 0.127** at α=0.7, and is **0.127 vs the leader's 0.008 —
about 16×**. I expected the opposite: `jac` divides by the union, so over-calling should hurt.
It does not, because the even emitter's extra calls are not random — the DE test now detects the
genuinely-shifted genes it previously missed for lack of power, so intersection grows faster than
union. This is the single largest per-metric edge we have over the field.

## CPM gate: a wash

| α | cpm20 | cpm5 |
|---|---|---|
| 0.7 | 0.1683 | 0.1685 |
| 1.0 | 0.1621 | 0.1575 |

The 5→20 CPM difference at α=0.7 is **+0.0001**, far below the harness's resolution floor
(`jac` seed noise alone is 6.7% of span). Keeping cpm20; the choice does not matter.

## Credit and caveat

The even-spread idea comes from **Sidechain** (live rank 5, public code at
`github.com/saberhq/sidechain`), who measured an even-spread null at −0.01 against −0.30 for a
resampled-control null. Our contribution is quantifying *why* it works and combining it with
α set for `nmae`.

**Caveat:** this is 40 perturbations on one harness. `pds` saturates at that panel size, so the
0.756 is optimistic and the true competition value is unknown. The `est` column is a projection,
not a prediction — treat the *ranking* of configs as the finding, not the absolute number.

# 159 — the cross-line transfer ceiling: the paradigm is at ~30%, and the +0.0091 was an artefact

**Date:** 2026-09-02 · deterministic (seed `20260902` everywhere) · **no `cell-eval2` was run, no
submission was built, no existing file was modified.** Every number below is a profile-space
pseudobulk statistic computed from cell-level h5ad or from `data/donors/*/delta.npy`.

---

## Verdict

Two findings, and the second overturns the brief's premise.

1. **The disattenuated cross-line transfer is real and large, not marginal.** Only ~30% of a
   target's target-specific response energy is shared across cell lines (κ² ≈ 0.34), so ~70% is
   line-specific — but 30% of a real signal is not the same thing as "barely above the decoy
   level".
2. **The brief's `+0.1881 same / +0.1790 different / +0.0091 excess` is a property of the
   thinned Orion harnesses, not of cross-line biology.** Those harnesses put **86.5% of their
   true effect energy into a single target-independent vector** that the same two lines' donor
   tables do not have. Measured on the same two cell lines, the same 120 targets and the same
   17,986 genes, with all their cells: **same `+0.0317`, different `+0.0006`, ratio `54×`** —
   not `1.05×`. The excess is not 5% above the decoy level; it is 5,000%.

**Production sits at 31.1% (HCT116) / 28.5% (HEK293T) of the infinite-donor
context-independent ceiling**, measured on the actual challenge panel. Not 90%. But the
headroom is expensive: the ceiling is approached as `1/k` in the number of *independent cell
lines*, and reaching even half of it needs about **ten** lines against the four we have.

---

## 0. Known-answer row — reproduced exactly (`known.py`)

Required before anything else. Reads 156's frozen exact-space truth read-only.

| quantity | brief | this experiment |
|---|---:|---:|
| same target, HCT116↔HEK293T, per-target cosine | +0.1881 | **+0.1881** |
| median / sd | 0.1965 / 0.1031 | **0.1965 / 0.1035** (sd is `ddof=1`) |
| different target, cross-line | +0.1790 | **+0.1790** |
| target-specific excess | +0.0091 | **+0.0091** |
| within-HCT116, different target | +0.2476 | **+0.2476** |
| within-HEK293T, different target | +0.4175 | **+0.4175** |
| energy split, shared / line-specific | 60.1% / 39.9% | **60.1% / 39.9%** |
| context-independent ceiling, cos to each line | 0.7700 | **0.7676** |

120 shared targets, 18,413 scored genes (18,533 minus the 120 panel target coordinates), 150
cells per target, 18,400 control cells per line. The drop set does not matter: dropping the
union of both panels (144 coords) gives the identical numbers to four decimals; dropping
nothing gives same +0.1910 / excess +0.0121.

**No bug on either side.** Everything after this is about what those numbers mean.

---

## 1. Split-half reliability, as a function of cells per target (`reliability.py`, `diag.py`)

### 1.1 Design

Each target's cells are randomly partitioned into K shards (`shards.py`, seed `20260902`;
10 shards × 15 cells on Orion, 20 × 20 on H1), so any two **disjoint** cell subsets of size
`m = j·shard` can be pseudobulked without re-reading the h5ad. Control cells are sharded the
same way, so each half is scored against a **disjoint half of the control** — sharing one
control between the halves correlates their noise. `r(m) = cos(e¹, e²)`, averaged over 24
disjoint pairings and over targets.

Two spaces, and the distinction is load-bearing:

- `full` — the raw effect vector, what `pds_cosine` and `expr_mse` read.
- `resid` — the effect minus the mean over targets **of that side at that depth**. This is the
  target-specific component, and it is the only part that can earn score (003 §"the
  decomposition that matters": the mean response *is* the metric's zero).

A first pass centred the residual on the full-depth mean and produced a **non-monotone** `r(m)`
(0.0178 → 0.0143 → 0.0153 → 0.0195 → 0.0218 on HCT116). That is the signature of the
depth-dependent `log1p` bias: a shallow pseudobulk carries a systematic offset that is common
to every target and to both halves, so it inflates the cosine. Per-depth, per-side centring
removes it exactly and the curve becomes monotone. **The published residual numbers are the
corrected ones.**

### 1.2 Measured curves

`r(m)` per side, disjoint controls, per-target cosine averaged over targets:

| m cells/side | HCT116 full | HCT116 resid | HEK293T full | HEK293T resid | H1-test full | H1-test resid |
|---:|---:|---:|---:|---:|---:|---:|
| 15 / 20 | 0.0402 | 0.0067 | 0.0761 | 0.0131 | 0.0858 | 0.0804 |
| 30 / 40 | 0.0697 | 0.0119 | 0.1460 | 0.0279 | 0.1347 | 0.1329 |
| 45 / 60 | 0.1038 | 0.0174 | 0.2077 | 0.0396 | 0.1724 | 0.1727 |
| 60 / 80 | 0.1357 | 0.0243 | 0.2608 | 0.0505 | 0.2040 | 0.2061 |
| 75 / 100 | 0.1637 | 0.0286 | 0.3070 | 0.0637 | 0.2292 | 0.2331 |
| 200 | — | — | — | — | 0.3178 | **0.3300** |

Distribution, not just the mean — HCT116 `resid` at m = 75: p10 −0.011, median 0.021, p90 0.086.
HEK293T: p10 −0.007, median 0.048, p90 0.150. A third of targets are at or below zero.

**The control's contribution is measured, not assumed.** In the `full` space, using the full
control for both halves instead of disjoint halves inflates `r` by **2.4%** (HCT116 0.1637 →
0.1676) — small, because 18,400 control cells against 150 target cells. In the `resid` space
the two variants agree to **every printed digit**, because centring on the across-target mean
removes any vector common to all targets of a side, and the control profile is exactly such a
vector. **The residual space is immune to control noise by construction.**

### 1.3 Extrapolating 150 → 400, and the check that licenses it

`h1-test` and `h1-full` carry **exactly 400 cells per target — the challenge's own count** — so
the extrapolation the 150-cell Orion surfaces need can be tested against a measured value.
Fit on `m ≤ 100`, predict the measured `m = 200`:

| surface | law | predicted | measured | ratio |
|---|---|---:|---:|---:|
| h1-test `full` | `1/r−1 = a/m + c` | 0.3632 | 0.3178 | ×1.143 |
| h1-test `full` | `1/r−1 = a/m^p` | **0.3304** | 0.3178 | **×1.039** |
| h1-test `resid` | linear | 0.3719 | 0.3300 | ×1.127 |
| h1-test `resid` | **power** | **0.3451** | 0.3300 | **×1.046** |

The pure-`1/n` Poisson law overshoots by 13–14%; the power law (fitted exponent `p = 0.70` on
H1, `0.93–1.05` on Orion — `log1p` is not variance-stabilising at these depths) is within 4–5%.
Both are reported; the spread between them is the honest error bar.

| | HCT116 | HEK293T | h1-test |
|---|---:|---:|---:|
| ρ(150), `full` | 0.262–0.276 | 0.462–0.480 | 0.281–0.343 |
| ρ(150), `resid` | **0.050–0.055** | **0.110–0.124** | 0.289–0.338 |
| ρ(400), `full` | 0.43–0.49 | 0.72–0.77 | 0.435–0.469 |
| ρ(400), `resid` | **0.099–0.123** | **0.255–0.273** | **0.481 (measured, not fitted)** |

**A 400-cell reference is still only 10–48% reliable on the target-specific component.** That
is the single most important number for interpreting every cosine anyone quotes, ours or the
board's.

### 1.4 A second estimator, because the first one is a ratio-of-expectations

`diag.py` re-derives ρ from moments only: `‖d‖² = E‖mean‖² − E‖diff‖²/4`, so
`ρ(2m) = ‖d‖²/E‖mean‖²`. Agreement with Spearman–Brown on the cosine estimator: **HCT116 ×1.037,
HEK293T ×1.069**. (On h1-test the two differ ×1.43, because the moment estimator is
energy-weighted and H1's effect sizes are heavy-tailed; the cosine estimator is the one that
matches how `pds` averages, and is used throughout.)

The same routine returns the energy budget, which is where §2 comes from:

| surface (per target, at its own depth) | observed | = shared | + target-specific | + noise |
|---|---:|---:|---:|---:|
| HCT116, 150 cells, 18,401 genes | 94.42 | **24.99** | 3.90 | 65.54 |
| HEK293T, 150 cells | 119.34 | **52.75** | 7.92 | 58.66 |
| h1-test, 400 cells, 17,995 genes | 46.69 | **1.76** | 30.47 | 14.47 |

---

## 2. The Orion harnesses carry a cell-pool blob, and it is what the brief measured

**86.5% (HCT116) and 86.9% (HEK293T) of the harnesses' *true* effect energy sits in a single
target-independent vector.** On h1-test the same statistic is **5.4%**. That is not a subtle
disagreement between surfaces; it is a different object.

`reconcile.py` walks from the donor table to the harness one variable at a time. Same two cell
lines throughout.

| configuration | n | genes | same | diff | ratio | within-A | within-B |
|---|---:|---:|---:|---:|---:|---:|---:|
| donor, 3,000 random targets, Orion 38.6k axis | 3000 | 35,605 | +0.0296 | +0.0018 | 16.1 | +0.0023 | +0.0089 |
| donor, 3,000 random targets, challenge axis | 3000 | 15,267 | +0.0308 | +0.0019 | 16.0 | +0.0024 | +0.0093 |
| donor, **the harness's 120 targets**, Orion axis | 120 | 38,486 | +0.0299 | +0.0006 | 46.7 | +0.0018 | +0.0021 |
| donor, **the harness's 120 targets, challenge axis** | 120 | 17,986 | **+0.0317** | **+0.0006** | **54.2** | +0.0018 | +0.0021 |
| **harness, 120 targets, challenge axis** | 120 | 17,986 | **+0.1881** | **+0.1790** | **1.05** | **+0.2476** | **+0.4175** |
| harness, 120 targets, donor's control substituted | 120 | 17,986 | +0.8202 | +0.8178 | 1.00 | +0.8736 | +0.8677 |

`‖mean effect‖` over the same 120 targets and 17,986 genes: **harness 4.93 / 7.16, donor
0.60 / 0.52** — an 8–14× inflation, and the two directions barely agree (cos +0.05 / +0.21).

Four candidate mechanisms were tested and three are ruled out:

- **not the gene axis** — the donor table restricted to the challenge axis is unchanged
  (within-line +0.0018 vs +0.0018);
- **not the target set** — 3,000 random targets and the harness's 120 give the same picture;
- **not the pseudobulk depth ratio** (`blob.py`) — sweeping the control from 920 to 18,400 cells
  moves `‖mean effect‖` only 6.30 → 5.00 on HCT116, and sweeping the target from 15 to 150 cells
  moves it 5.38 → 5.00. The blob is **invariant to both depths**, so it is not a `log1p`
  sampling bias;
- **not gene-axis renormalisation** — inverting the donor table to profiles, restricting to the
  18,533-gene axis and *renormalising* there (exactly what the harness pseudobulk does) moves
  `‖mean‖` 0.586 → 0.686 and within-line cosine 0.0020 → 0.0023. Nothing.

What is left is the harnesses' **cell selection**: `build_orion_harness.py --sampling greedy`
fills each perturbation with 150 cells and the control pool separately. This is the same family
as `docs/HARNESS-DEPTH-CONFOUND.md`, and **the depth gate does not catch it** — the thinned
build passes at ratio 0.9999 and still carries the blob. `data/harness/orion-hct116-thinned.h5ad`
and `orion-hek293t-thinned.h5ad` have no batch column, so the offset cannot be regressed out
from the harness alone.

> **Consequence for this repo, stated plainly.** Any statistic on the thinned Orion harnesses
> that reads the *un-centred* effect — absolute cosine, `pds` decoy level, the row-row cosine of
> the truth, "shared vs target-specific energy" — is distorted. Statistics computed on
> **mean-removed** effects, and all *relative* comparisons between arms, are not: `reconcile.py`
> §D shows harness and donor agree in kind once each line's mean is removed (+0.0166 vs +0.0311,
> a factor 1.9 attributable to the harness's UMI thinning and the donors' extra cells). 156's
> and 099's arm-to-arm deltas are in the second category. This is a *new* pathology on top of
> rule 9, and it should go in `docs/`.

---

## 3. Disattenuation (`crossline.py`)

### 3.1 The estimator and its assumptions

With `d_obs = d + ε`, ε zero-mean, independent of `d`, independent between lines and between
disjoint cell sets:

```
E cos(d_A + eps_A, d_B + eps_B)  ~=  cos(d_A, d_B) * sqrt(rho_A * rho_B)
```

Two consequences, and the first corrects the brief:

1. **The attenuation multiplies the same-target and the different-target cosine identically**,
   because it lives entirely in the denominators and the noise does not know which pair it is
   in. So the *ratio* same/diff is attenuation-invariant, and the excess and the decoy level
   disattenuate by the same factor. The brief's expectation that "noise suppresses the EXCESS
   far more than the baseline" is **true in the residual space** (where the target-specific
   vector's own low SNR is what attenuates) and **false in the full-vector space**.
2. **That invariance is testable, and it was tested.** Sub-sampling the same cells sweeps
   `sqrt(ρ_A ρ_B)` from 0.053 to 0.347 — a factor 6.5:

| m cells/target | 15 | 30 | 45 | 60 | 75 | 90 | 105 | 120 | 135 | 150 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| same (full) | 0.0242 | 0.0466 | 0.0692 | 0.0928 | 0.1102 | 0.1292 | 0.1462 | 0.1599 | 0.1758 | 0.1881 |
| diff (full) | 0.0222 | 0.0427 | 0.0647 | 0.0875 | 0.1049 | 0.1224 | 0.1392 | 0.1521 | 0.1669 | 0.1790 |
| **ratio** | 1.085 | 1.093 | 1.069 | 1.062 | 1.050 | 1.056 | 1.050 | 1.051 | 1.053 | **1.051** |
| `sqrt(rA·rB)` | 0.053 | 0.099 | 0.141 | 0.184 | 0.220 | — | — | — | — | 0.347 |
| resid same | 0.0025 | 0.0048 | 0.0063 | 0.0077 | 0.0081 | 0.0106 | 0.0118 | 0.0135 | 0.0157 | 0.0166 |
| resid same / attenuation | 0.279 | 0.271 | 0.258 | 0.234 | 0.202 | — | — | — | — | **0.207** |

The ratio is flat to ±0.02 over a 6.5× change in attenuation. **The disattenuation model holds
on this data**, which is what licenses reading the donor-pair matrix in §4 through the ratio.

Where it is weakest: ρ varies across targets and is *correlated* across lines (a strong
perturbation is strong in both), so `E[sqrt(ρ_At ρ_Bt)]` on the diagonal exceeds
`E[sqrt(ρ_At)]·E[sqrt(ρ_Bu)]` off it. That biases the observed excess **up**, not down — which
is the residual decline from 0.279 to 0.207 in the last row.

### 3.2 Full depth, target-cluster bootstrap (2,000 draws, 120 target clusters)

| space | attenuation | same obs → disatten | diff obs → disatten | **excess obs → disatten** |
|---|---:|---|---|---|
| `full` | 0.3469 | +0.1881 → **+0.5422** [+0.5064, +0.5736] | +0.1790 → +0.5160 [+0.4794, +0.5483] | +0.00908 [+0.00395, +0.01492] → **+0.02617 [+0.01104, +0.04303]** |
| `resid` | 0.0808 | +0.0166 → **+0.2050** [+0.0665, +0.3348] | −0.0001 → −0.0016 | +0.01670 [+0.00513, +0.02801] → **+0.20661 [+0.06528, +0.33557]** |

**The `full` row is contaminated by §2's blob and should not be quoted.** The `resid` row is
the honest one: the disattenuated cross-line transfer of the target-specific programme on this
surface is **τ ≈ 0.21, 95% CI [0.065, 0.336]**, and the depth trend says the true value is at
the low end of the point estimate rather than above it.

### 3.3 Three independent routes to τ, and they agree

| route | τ (cross-line cosine of the true target-specific effect) |
|---|---:|
| harness residual split-half disattenuation (§3.2) | **0.207** [0.065, 0.336] |
| donor single-pair / same-line control ratio: `cos(K562←HCT116) / cos(K562←K562-ess)` = 0.0897 / 0.2627 | **0.341** |
| `κ²` from the k→∞ consensus extrapolation (§5), `κ = 0.591` | **0.349** |

Call it **τ = 0.30 ± 0.10**, i.e. `κ = sqrt(τ) ≈ 0.55`: **about 30% of a target's
target-specific response energy is line-independent and about 70% is line-specific.** The
harness route reads low, consistently with its UMI thinning.

---

## 4. The line-pair matrix (`donorpairs.py`)

All ten donor tables, every pair with ≥30 shared perturbations, gene axes intersected, panel
target coordinates dropped, capped at 3,000 targets. `RATIO` is the attenuation-invariant
column; `resid null` is the same statistic with the target correspondence permuted.

| line A | line B | n pert | n gene | same | diff | excess | **RATIO** | within-A | within-B | resid same | resid null |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **k562_essential** | **k562_gwps** | 2003 | 6402 | +0.2712 | +0.0329 | +0.2384 | 8.26 | +0.0383 | +0.0429 | **+0.2857** | +0.0055 |
| h1_2025 | orion_hek293t | 247 | 17454 | +0.0265 | +0.0002 | +0.0264 | 165.2 | +0.0279 | +0.0106 | +0.0275 | +0.0022 |
| h1_2025 | orion_hct116 | 244 | 17457 | +0.0335 | +0.0014 | +0.0321 | 23.2 | +0.0280 | +0.0026 | +0.0332 | −0.0008 |
| h1_2025 | k562_gwps | 237 | 7359 | +0.0639 | +0.0029 | +0.0611 | 22.4 | +0.0293 | +0.0078 | +0.0653 | −0.0019 |
| **orion_hct116** | **orion_hek293t** | 3000 | 35584 | +0.0306 | +0.0016 | +0.0290 | **19.5** | +0.0025 | +0.0067 | +0.0300 | −0.0032 |
| k562_gwps | orion_hek293t | 3000 | 5642 | +0.0289 | +0.0016 | +0.0273 | 17.6 | +0.0061 | +0.0055 | +0.0286 | +0.0002 |
| k562_gwps | orion_hct116 | 3000 | 5662 | +0.0321 | +0.0020 | +0.0300 | 15.7 | +0.0063 | +0.0035 | +0.0321 | +0.0003 |
| cd4_rest | cd4_stim48hr | 276 | 9839 | +0.2249 | +0.0221 | +0.2028 | 10.2 | +0.0317 | +0.0276 | +0.2078 | −0.0001 |
| k562_essential | orion_hek293t | 1807 | 6609 | +0.0941 | +0.0114 | +0.0827 | 8.22 | +0.0356 | +0.0175 | +0.0983 | +0.0006 |
| k562_essential | orion_hct116 | 1800 | 6617 | +0.0999 | +0.0123 | +0.0876 | 8.11 | +0.0357 | +0.0248 | +0.1032 | +0.0008 |
| k562_essential | nadig_jurkat | 1916 | 5856 | +0.1552 | +0.0229 | +0.1323 | 6.77 | +0.0364 | +0.0444 | +0.1573 | +0.0040 |
| k562_gwps | nadig_jurkat | 2252 | 5373 | +0.1495 | +0.0241 | +0.1254 | 6.20 | +0.0408 | +0.0440 | +0.1510 | +0.0051 |
| nadig_jurkat | orion_hct116 | 2047 | 6731 | +0.0786 | +0.0129 | +0.0658 | 6.11 | +0.0420 | +0.0196 | +0.0778 | +0.0005 |
| nadig_jurkat | orion_hek293t | 2055 | 6723 | +0.0804 | +0.0133 | +0.0671 | 6.06 | +0.0420 | +0.0142 | +0.0809 | +0.0070 |
| orion_hek293t | rpe1 | 1993 | 6780 | +0.0692 | +0.0149 | +0.0544 | 4.65 | +0.0140 | +0.1726 | +0.0705 | +0.0036 |
| nadig_jurkat | rpe1 | 2123 | 5544 | +0.1259 | +0.0330 | +0.0929 | 3.82 | +0.0420 | +0.1890 | +0.1215 | +0.0027 |
| k562_essential | rpe1 | 1875 | 5598 | +0.1333 | +0.0354 | +0.0979 | 3.76 | +0.0376 | +0.2067 | +0.1342 | +0.0033 |
| k562_gwps | rpe1 | 2187 | 5220 | +0.1244 | +0.0358 | +0.0886 | 3.47 | +0.0398 | +0.1903 | +0.1164 | +0.0043 |
| orion_hct116 | rpe1 | 1986 | 6788 | +0.1200 | +0.0365 | +0.0836 | 3.29 | +0.0189 | +0.1737 | +0.1242 | +0.0072 |
| kolf2 | orion_hek293t | 268 | 17838 | −0.0218 | −0.0247 | +0.0028 | 0.89 | +0.2423 | +0.0068 | +0.0033 | +0.0008 |
| k562_gwps | kolf2 | 249 | 7333 | −0.0087 | −0.0151 | +0.0064 | 0.57 | +0.0049 | +0.2836 | +0.0089 | +0.0011 |

(13 further pairs in `work/donorpairs.json`.)

Readings:

- **The same-line positive control is the anchor.** K562-gwps ↔ K562-essential — one cell line,
  two libraries, two experiments — reads `+0.2712` where every genuine cross-line pair reads
  `+0.03` to `+0.16`. A single donor's measurement of its **own** line is 2–9× more informative
  than another line's, at comparable cell counts. That difference is the size of the
  context-conditional gap.
- **Every cross-line pair beats its permutation null by 10–70×** (`resid null` is 0.000–0.007
  against `resid same` 0.003–0.29). The transfer is unambiguously real everywhere.
- **`kolf2` is the exception and it is anti-aligned**, reproducing 099 §3.3 (κ −0.0027/+0.0035)
  and 057 from a third statistic. Its within-line cosine is +0.24–0.28, i.e. **`kolf2` has its
  own blob** — it is the one donor table with the harness's pathology.
- **`rpe1` also has a blob** (within-line +0.17–0.21), which deflates every `rpe1` ratio in the
  table. Its 003-era number (K562↔RPE1 transfer +0.1416, centred +0.1275) is reproduced here as
  `resid same +0.1164` on 2,187 shared perturbations — 003's measurement stands.
- The `orion_hct116 ↔ orion_hek293t` row is the same pair the brief measured. **Ratio 19.5 on
  3,000 targets, 54.2 on the harness's 120 — against the harness's 1.05.**

Not disattenuable from their own data (no cell-level split is possible on a donor table), which
is exactly why the ratio column is the primary one.

---

## 5. The ceiling, by construction (`ceiling.py`, `place.py`)

### 5.1 What the ceiling is

A context-independent model emits one effect matrix and applies it to every context. Given
infinitely many donor lines, the best such matrix is the **across-line mean effect** for each
target — the conditional expectation of a line's effect given only the target's identity. So
the ceiling needs no model: build the equal-weight consensus of `k` other lines, score it, and
extrapolate `k → ∞`.

Under "shared component + independent per-line component + independent noise" the consensus
obeys `1/cos²(k) = A + B/k` exactly. Measured, that fit has **R² = 0.99998**.

Everything is mean-removed per line first: the across-target mean *is* the metric's zero, so no
predictor may be credited for it.

### 5.2 On the six-line intersection (1,592 targets **including common-essential genes**, 5,206 genes)

| reference | k=1 | k=2 | k=3 | k=4 | k=5 | **k→∞** |
|---|---:|---:|---:|---:|---:|---:|
| orion_hct116, cos | 0.1065 | 0.1382 | 0.1579 | 0.1718 | 0.1822 | **0.253** |
| orion_hct116, pds | 0.7861 | 0.8207 | 0.8375 | 0.8492 | 0.8559 | — |
| k562_gwps, cos | 0.1132 | 0.1492 | 0.1728 | 0.1900 | — | **0.303** |
| nadig_jurkat, cos | 0.1232 | 0.1602 | 0.1841 | 0.2013 | 0.2146 | — |

Same-line positive control on this surface: predict K562-gwps from K562-essential → **pds
0.9468, cos 0.2627**; predict it from HCT116 → pds 0.7839, cos 0.0897.

That control also gives a second, independent read on donor reliability: if two libraries of the
same line have true cosine 1, then `sqrt(ρ_gwps·ρ_ess) = 0.2627`, so `ρ ≈ 0.26`, and
`κ = cos_∞/sqrt(ρ) = 0.303/0.513 = 0.591`, hence `τ = κ² = 0.349` — the third route in §3.3.

### 5.3 On the **actual challenge panel** — the number that matters

`place.py`. Reference = the Orion donor table for the held line; prediction = 156's frozen
production artifact (`104-.../work/<line>/frozen.npz`, built from the *legal* donor pair, so
scoring against the held line is leakage-free); four independent donor lines
(other-Orion, K562-gwps, KOLF2, CD4-T-as-one-line). 116 targets × 6,535 genes.

`rpe1`, `nadig_jurkat` and `k562_essential` cannot appear here: they carry **0 of the 300**
challenge targets (099 §1.1 — the 2026 panel excludes common-essential genes). The
challenge-panel cosines are 3–5× smaller than §5.2's, and **two things are confounded in that
gap**: the panel excludes the common-essential genes whose strong, housekeeping-driven responses
transfer best (003 §"best-transferring perturbations"), and the donors that survive the panel
filter include two weak ones (`kolf2`, which §4 shows is anti-aligned and carries its own blob,
and CD4-T). This experiment does not separate them, and the k→∞ fit is run on the surface as it
is, which is the conservative choice.

| arm | pds | cos | sse* | row-row |
|---|---:|---:|---:|---:|
| **production (shipped effect matrix)** | **0.8262** | **0.03814** | 0.998493 | −0.00842 |
| consensus of 1 donor | 0.6501 | 0.02239 | 0.999393 | −0.00868 |
| consensus of 2 donors | 0.6872 | 0.03108 | 0.998947 | −0.00867 |
| consensus of 3 donors | 0.7044 | 0.03753 | 0.998510 | −0.00867 |
| consensus of 4 donors (all we have) | 0.7214 | 0.04278 | 0.998082 | −0.00866 |
| **k → ∞ ceiling** | — | **0.12262** | 0.98497 | — |
| oracle (the reference's own truth) | 1.0000 | 1.00000 | 0.000000 | — |

HEK293T reads the same: production cos 0.03473, 4-donor 0.04357, **k→∞ 0.12182**.

Robustness of the extrapolation: sub-fits on `k∈{1,2,3,4} / {2,3,4} / {3,4}` give **0.1226 /
0.1317 / 0.1346**; a 300-draw target-cluster bootstrap gives **median 0.1225, 95% CI [0.0845,
0.1610]** (HEK293T: 0.1248, [0.0839, 0.1608]).

> **`production is at 31.1% (HCT116) and 28.5% (HEK293T) of the context-independent ceiling`**,
> and at **89.1% / 79.7% of the best consensus that can actually be built from the four donor
> lines that carry the panel.**

### 5.4 The pds-vs-transfer-level curve

Built on the same surface: `P_t = τ·s_t + sqrt(1−τ²)·z_t`, `z_t` a random combination of the
**other** targets' truths orthogonalised against `s_t`, so the decoy geometry is the real one
and only the target-specific information is dialled.

| cos to truth | 0.000 | 0.020 | 0.030 | **0.038** | **0.043** | 0.050 | **0.072** | 0.100 | **0.123** | 0.150 | 0.200 | 0.250 | 0.300 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **pds (HCT116 surface)** | 0.5043 | 0.5730 | 0.6096 | **0.6374** | **0.6538** | 0.6765 | **0.7430** | 0.8199 | **0.8705** | 0.9175 | 0.9723 | 0.9920 | 0.9981 |
| pds (HEK293T surface) | 0.4998 | 0.5594 | 0.5879 | 0.6091 | 0.6231 | 0.6414 | 0.6960 | 0.7650 | 0.8112 | 0.8571 | 0.9200 | 0.9567 | 0.9766 |
| `sse*` = 1−cos² | 1.0000 | 0.99960 | 0.99910 | 0.99856 | 0.99815 | 0.99750 | 0.99482 | 0.99000 | 0.98487 | 0.97750 | 0.96000 | 0.93750 | 0.91000 |
| marker | baseline | | | **PRODUCTION** (on-surface) | 4-donor consensus | | **PRODUCTION** (board-equivalent) | | **k→∞ CEILING** (on-surface) | | | | |

**The curve validates itself against the board.** Our board-implied aggregate cosine is 0.0722
(brief: raw `expr_mse` 0.8738 through `raw_mse ≈ 1.0237 − 2.077·cos`). The curve at cos = 0.072
predicts raw `pds` **0.743**; the board reads our raw `pds` at **≈0.73–0.74**. A synthetic
construction on a held cell line reproduces a live board number to 0.01 — that is the
known-answer row for §5.4.

Converting the on-surface ceiling to board units by the ratio our own single prediction spans
between the two surfaces (0.0722 / 0.03814 = ×1.893, a ratio in which the surface factor
cancels because it is the *same* prediction):

| quantity | production | leader | **context-independent CEILING** |
|---|---:|---:|---:|
| aggregate cos (board scale) | 0.0722 | 0.119 | **0.232** [0.160, 0.305] |
| min raw `expr_mse`, board map | 0.8738 | 0.7767 | **0.542** |
| raw `pds`, curve read at that cos | 0.743 | 0.850 | **~0.98** |
| raw `pds`, actually observed on the board | ≈0.728 | 0.7835 | — |
| **% of ceiling** | **31%** | **51%** | 100% |

Two honest caveats on that table. **The curve tracks the board at our own point (0.743 predicted
vs ≈0.728 observed) and overshoots at the leader's by 0.07** (0.850 vs 0.7835), so the ~0.98 at
the ceiling is an upper reading and the real curve flattens sooner. And `sse* = 1−cos²` and the
board map disagree (0.946 vs 0.542) because the scored `mse` carries a jackknife correction and
per-perturbation normalisation the proxy does not; both are quoted, and the board map — fitted
to real board entries — is the one to use.

### 5.5 How expensive the remaining 69% is

From `1/cos²(k) = A + B/k` with `A = 66.5`, `B = 1,920` (`B/A = 28.9`):

| target | independent cell lines needed |
|---|---:|
| 50% of the ceiling | **9.6** |
| 70% of the ceiling | **27.7** |
| 90% of the ceiling | **123** |

And inverting the same fit for what our arms are *worth* in donor-equivalents:

| arm | cos (on-surface) | effective independent donor lines |
|---|---:|---:|
| production | 0.03814 | **3.09** |
| 4-donor equal consensus | 0.04278 | 4.00 |
| board leader (board 0.119, rescaled) | 0.0629 | **10.3** |

**Production's whitening + PerturbMatch + CD4 stack is worth about three ideal equal donors,
which is what it uses.** The leader is operating at an effective ten. This is the sharpest
quantitative statement in the experiment about the gap to first place, and it says the gap is
*not* a modelling subtlety — it is roughly a factor of three in effective independent
biological views.

**One thing production is already doing well, and it should not be lost.** At the same cosine
(0.038) the isotropic-error curve gives pds **0.637**; production gets **0.826**. That +0.19 is
the whitening geometry of AGENTS.md §4.2 and 156's mass fix, and it means production converts
cosine into retrieval far more efficiently than a raw consensus. The 4-donor consensus has a
*higher* cosine (0.0428) and a *lower* pds (0.7214) than production. On `pds` specifically,
production is much closer to the ceiling than the cosine ratio suggests.

---

## 6. The answer to the strategic question

**Production is at about 30% of the donor-transfer paradigm's ceiling, not 90% — but the
paradigm's ceiling is itself low, and the brief's evidence that it had been reached was an
artefact of the measurement surface.** The `+0.0091` excess is not a fact about cross-line
biology: the two thinned Orion harnesses put 86.5% of their true effect energy into a
target-independent cell-pool vector that the same two lines' donor tables do not have, and once
that vector is removed the same comparison reads **same +0.0317 against a decoy level of
+0.0006, a ratio of 54 rather than 1.05**, and every one of the 33 donor line pairs beats its
own permutation null (by 10–70× for the genuine cross-line pairs, by only 2–8× for the five
`kolf2` pairs, which are the exception in every other respect too). What is genuinely true, and
what all three independent estimators agree on, is that **only about 30% of a target's
target-specific response energy is line-independent** (τ = 0.30 ± 0.10, κ ≈ 0.55) — so a
context-independent matrix is structurally capped at an aggregate cosine of ≈0.23 on the board's
scale against ≈0.072 achieved, i.e. raw `expr_mse` ≥ 0.54 and raw `pds` ≲ 0.98, and it is
provably incapable of expressing the remaining ~70%. The catch is that the ceiling is approached
as `1/k` in the number of *independent cell lines*: reaching half of it needs about ten lines and
we have four that carry the panel (the common-essential libraries carry zero of the 300 by
construction of the panel), so the honest read is that **the cheap part of donor transfer is
close to spent — production is already worth about three ideal donors and 89% of the best
four-donor consensus buildable today — while the expensive part is a data-acquisition problem,
not a modelling one.** The class of information that would exceed the ceiling is
**context-conditional response**: a model whose predicted effect for target *t* is a function of
the context's basal profile as well as *t*, which is the only object that can reach the ~70% of
energy the line-specific component holds. The size of that prize is measurable and it is large —
the same-line control (K562 predicted from a second K562 library) reads cos **0.2627** where the
best cross-line single donor reads **0.0897**, so knowing the *context* is worth roughly 3× what
knowing the *target* is worth, on identical cell counts. That is where the next two months
should go, and the first thing to build is not another donor but the **cross-line residual
predictor**: given (basal profile of line X, target *t*, consensus effect), predict the
line-specific residual. Every one of the ~40 NO-GO experiments improved the consensus; none of
them attempted the residual.

---

## 7. What is surface-limited, and what is not (AGENTS.md §11.1 rule 9)

**Not computed at all, deliberately:** `fid`, `jac`, and `reach`. All three carry `n_real` in a
denominator or a union and are unreadable on every harness we hold. Nothing in this document
should be read as a claim about them.

**Surface-limited, do not quote as absolutes:**

- Anything from the thinned Orion harnesses that reads the *un-centred* effect (§2). That
  includes the brief's `+0.1881 / +0.1790 / +0.0091`, the 0.2476 / 0.4175 within-line cosines,
  the 60.1/39.9 energy split and the 0.7676 "ceiling".
- The `full`-space row of §3.2.
- Absolute `pds` on any surface here: the challenge has 300 decoys, these surfaces have 113–300,
  and `pds` depends on the decoy count. Only the §5.4 curve, which holds the panel at 116–300
  and was validated against a live board `pds`, should be read for level.
- ρ(400) for the Orion lines is an **extrapolation** from a 75-cell measurement (5.3× beyond the
  data). It is bracketed by two laws whose predictions differ by 13% on the one surface where a
  400-cell value is directly measurable (h1-test), and h1-test's own 400-cell residual
  reliability (0.481) is quoted in preference wherever a single number is needed.
- The board-scale conversion (×1.893) assumes our own prediction's surface ratio applies to the
  ceiling and to the leader. It cancels the surface factor for *our* point by construction; for
  the leader's it is an assumption.

**Not surface-limited:**

- The known-answer reproduction (§0) — exact.
- The attenuation-invariance of same/diff, demonstrated over a 6.5× attenuation sweep (§3.1).
- The donor-pair matrix ratios (§4) and their permutation nulls.
- The `1/cos²(k) = A + B/k` consensus law (R² = 0.99998, three sub-fits, target bootstrap).
- The three-route agreement on τ (§3.3), which spans two surfaces and two unrelated estimators.

**Not run:** `shards.py h1full` (255 targets × 400 cells) was started and did not complete; it
would only have added statistical weight to the 400-cell reliability anchor that h1-test already
provides, and nothing here depends on it.

---

## 8. Files

| file | what it does |
|---|---|
| `known.py` | §0 — reproduces the brief's row from 156's frozen truth, three drop-set variants |
| `shards.py` | per-target/per-control count shards, streamed with h5py, seed 20260902 |
| `reliability.py` | §1 — `r(m)` curves, both spaces, both control variants, the depth-law fits and the h1 extrapolation check |
| `diag.py` | §1.4 — moment-based ρ, and the shared / target-specific / noise energy budget |
| `crossline.py` | §3 — the depth sweep, the disattenuation, the target-cluster bootstrap |
| `reconcile.py` | §2 — harness vs donor table, one variable at a time |
| `blob.py` | §2 — control-depth and target-depth sweeps of the harness blob |
| `donorpairs.py` | §4 — the line-pair matrix and its permutation nulls |
| `ceiling.py` | §5.2 — consensus-of-k on the six-line intersection, same-line control, synthetic curve |
| `place.py` | §5.3–5.5 — production on the challenge panel, k→∞ fit, the marked pds curve |
| `work/*.json` | every number above |

`work/shards-*.npz` (311 MB, gitignored) are a regenerable cache; `shards.py` rebuilds them in
~10 s per surface and skips any that exist.

```bash
cd experiments/159-transfer-ceiling
uv run python known.py                                   # the known-answer row, first
uv run python shards.py hct116 hek293t h1test            # ~30 s, peak RSS 2.7 GB
uv run python reliability.py && uv run python diag.py
uv run python crossline.py && uv run python reconcile.py && uv run python blob.py
uv run python donorpairs.py && uv run python ceiling.py && uv run python place.py
```

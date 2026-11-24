# 173 — sign repair: set whitening-flipped signs to the un-whitened consensus. **NO-GO** (net ≈ +0.002, profile-class gate +0.0089)

Date: **2026-09-02**, coordinator-run. Measurement only. No `cell-eval2`, no build, no
submission. Surface: `h1-test`, a15 emission law, donors `k562_gwps` + `orion_hek293t` (the 095
artifact's; `h1_2025` never enters). Scripts: `signfix173.py` (delivered-space `pds`/cos, 168d's
exact pipeline), `realize173.py` (169's realisation pipeline with the transform swapped),
`six173.py` (169's six-member reader). Peak RSS ~4 GB.

## The finding that motivated it

172 reported that whitening costs 0.011 pooled sign accuracy on the reference DE pool while buying
`pds`. Decomposed here against the un-whitened consensus `cons` (172's `h1test-feats.npz`):

| on the reference DE pool (65,044 pairs) | share | H sign right | cons sign right |
|---|---:|---:|---:|
| `sign(H) == sign(cons)` | 83.0% | 0.5532 | 0.5532 |
| **whitening FLIPPED the sign** | **15.7%** | **0.4556** | **0.5444** |
| `cons == 0` (no donor coverage) | 1.3% | 0.5111 | — |

**On the 15.7% of pool genes where whitening flipped the sign, `H` is right 45.6% of the time —
worse than chance — and the raw consensus is right 54.4%.** `cons` beats `H` on sign at every
`|cons|` quantile; the flips concentrate where `|cons|` is small (35.1% of the bottom quartile,
1.4% of the top decile) and carry only **1.6% of `H`'s energy**. In the top `|H|` decile the flipped
genes are confidently wrong (28% right) but are 0.5% of that decile.

## Arms

`H' = H` with the sign set to `cons`'s wherever they disagree, magnitude kept (`flip`); zeroing
those entries instead (`zero`); and a **null** that flips the same per-row count of random gated
entries (`permflip`). All through the production `mass_consistent_effect` (`minnorm`), rendered,
emitted with a15's exact law (`tail 1.0 / keep_top 8000 / keep_spread 0.25 / φ 4`), two matched
seeds on base and flip.

## Known-answer rows

`KA_a15_raw`: `reach` **0.18283** (want 0.18283) PASS; median `n_pred` 6034 against a stored 6035
— 169's own KA row also reads 6034, so this is the stored value's rounding, not a drift. Recovery
round trip `|base_mf − KA|` = 0.0015 against a re-seed spread of 0.016, PASS. `signfix173.py`
reproduces 168d's h1-test base (`pds` 0.86134).

## Results

**Delivered space (`signfix173.py`, 168d's pipeline, 17,995 scored genes):**

| arm | pds | Δ pds | cos | Δ cos |
|---|---:|---:|---:|---:|
| base H | 0.86134 | — | 0.03531 | — |
| **flip** | 0.85714 | **−0.00420** | 0.03806 | **+0.00275 (+7.8%)** |
| zero | 0.86008 | −0.00126 | 0.03706 | +0.00174 |
| **null: random flips, same count** | 0.82745 | **−0.03389** | 0.02560 | −0.00971 |

**DE members (`six173.py`, h1-test, a15 emission), deltas vs `base_mf`:**

| arm | Δ prec / fid raw | Δ nmae raw | Δ jac ovr | Δ reach raw f=1 | call-set sign acc |
|---|---:|---:|---:|---:|---:|
| flip_mf | **+0.00313** | **−0.00417** | −0.00016 | −0.02019 | 0.5454 |
| flip_mf_b (twin) | **+0.00315** | **−0.00447** | +0.00009 | −0.00897 (vs its twin) | 0.5457 |
| zero_mf | +0.00045 | −0.00221 | −0.00028 | −0.00303 | 0.5390 |
| permflip_mf (null) | −0.00195 | **+0.01129** | +0.00007 | −0.01234 | 0.5220 |
| seed spread (base vs base_b) | 0.00005 | 0.00039 | — | **0.01628** | 0.00007 |

`fid` and `nmae` are resolved to better than a tenth of their deltas across matched seeds. `reach`
is not: matched-seed deltas are −0.020 and −0.009 against a seed spread of 0.016.

## Pricing, all six, in overall board units

| member | estimate | basis |
|---|---:|---|
| `fid` | **+0.00175** | precision +0.0031, yield cap saturated on the challenge (169's convention) |
| `nmae` | **+0.00116** | −0.0043 raw / span 0.60 / 6 |
| `jac` | −0.0001 | — |
| `mse` | ≈ **+0.0022** | cos +7.8% → ΔX through 096's mechanism (`X` 0.01885, `sumD` 0.2286); not directly measured |
| `pds` | ≈ **−0.0016** | −0.0042 raw on h1-test; 165 puts h1-test's profile-channel fidelity at ~1.5× |
| `reach` | ≈ −0.0015, **unresolved** | −0.0145 mean matched-seed × 0.393 transfer; inside seed noise; 165 says h1-test is blind to reach on *profile* changes anyway |
| **net** | **≈ +0.002** | |

**This is a profile change, so rule 13's +0.0089 gate applies. NO-GO.**

## Mechanism, and why `reach` does not respond to +0.014 pooled sign accuracy

The sign gain lives in the low-`|cons|`, low-`|H|` entries — the ones with the weakest emitted
significance, i.e. the **back** of the p-value order. `reach` reads only the front, where flips are
0.5–1.7% of genes; so the pool-wide sign gain is invisible to it, while the mass fix's small
redistribution jitters the front and shows up as noise. This is 171's finding from the other side:
**`reach` is bounded by front purity, and nothing that acts on the tail can move it.**

The `pds` cost is the same trade as gamma. Whitening decorrelates the 85 effect vectors partly by
flipping small entries into low-variance directions; putting them back re-correlates the vectors
and returns some hubness. The exchange rate here (−0.054 raw `pds` per +100% cos) is essentially
097's for gamma on h1-test (−0.068). **Sign repair is surgical partial un-whitening, and it lands on
the same curve.** Its one advantage is that it is a mechanism (transfers to D/E/F) rather than a
fitted exponent; its disadvantage is that, like gamma, it nets ~+0.002.

## What survives

- **The measurement**: whitening makes 15.7% of the reference's DE genes *anti*-informative on
  sign (45.6%). This is the mechanism behind 172's "whitening costs 0.011" and 097's `kappa` loss.
- **A composable +0.002.** If a future arm buys `pds` back (the only thing sign repair and gamma
  cost), both become shippable riders. Neither is shippable alone.
- The `zero` variant costs a third of the `pds` for two-thirds of the cos gain; if a rider is ever
  wanted, price it too.

No slot spent. a15 unchanged.

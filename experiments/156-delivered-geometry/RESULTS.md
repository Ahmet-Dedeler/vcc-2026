# 156 — the effect the scorer reads is not the effect we designed

**Date:** 2026-09-02 · two held Orion lines, exact-space pseudobulk truth, deterministic
(no emitter, no seeds).

## Verdict

**The renderer silently re-injects the shared component that whitening exists to remove**, and
fixing it is worth **+13.6% / +56.0%** aggregate transfer cosine and **+0.00035 / +0.00278** raw
pds, with **every one of the six members moving in our favour**. Shipped as `a23`
(`VCC_MASS_FIX=minnorm`); prediction and falsifiers pre-registered in `A23-PREDICTION.md`.

A second, independent consequence: **every effect-level screen in this project since experiment
017 has scored an object the board never sees.** `identity/intended` beats
`identity/delivered` by `+0.00324 / +0.00382` raw pds, so the harness and the board have been
looking at different geometry all along. That is a fact about the instrument, not about a arm.

## 1. The mechanism

`profile_from_effect` renders a log-space effect as `p ∝ expm1(b_ctrl + gated)` and normalises
`p` to sum 1. The scorer reads back `δ̂ = log1p(TS·p) − b_ctrl`. The normalisation is the
identity **iff**

```
S(h) := sum_g clip(expm1(b_ctrl + h_g), 0)  ==  TS      (TS = 5e4)
```

`b_ctrl` is itself built as `log1p(TS · totals/totals.sum())`, so `S(0) = TS` **exactly** — I
verified `sum(expm1(b_ref)) = 50000.0` and the control round-trip `‖b_recov − b‖ = 0.0` on both
surfaces. Every *reference* perturbation satisfies it too, so **the truth lives on that
manifold.** Our effect does not: the gate, the alpha scaling, the whitening and the CD4 addition
all break it, and the measured mass ratio is **1.00181 / 1.00168**.

So the renderer applies its own correction — a division in the linear space, which in log space
is a shift concentrated on the highest-expressed genes and **nearly the same direction for every
target**:

| | HCT116 | HEK293T |
| --- | ---: | ---: |
| render residual, as % of row norm | **3.96%** | **4.53%** |
| row-row cosine of the *residual* | **0.303** | — |
| row-row cosine of the *intended* effect | **0.0040** | 0.0035 |
| row-row cosine, delivered | 0.00564 | 0.00480 |

The contamination is **75× more mutually-aligned than the signal it is added to.** `pds` is a
retrieval metric and hubness is exactly what it punishes; whitening's entire job (AGENTS.md
§4.2) is to make those rows mutually dissimilar. The renderer was putting the blob back.

## 2. Why scaling cannot fix it, and what can

`expm1` is convex and `S(0) = TS`, so `S` has **no second root along the ray** `{γ·h}` — there is
no `γ > 0` that restores the mass, and the correction *must* move off the ray. But it is one
scalar constraint, so its direction is free. With `dS/dh_g = exp(b_g + h_g) =: w_g`, minimising
`Σ ε_g²/u_g` subject to `⟨ε, w⟩ = −(S − TS)` gives `ε = λ·u·w`, a one-parameter family in the
weight `u`. `u = 1` is minimum-L2 (`minnorm`); `u = 1/w` is a uniform log shift (`const`).

The correction is confined to the CPM gate, so a gene the gate zeroed stays at exactly zero.

## 3. Result — production-exact, both surfaces

Effects gated at a15's shipped `min_cpm = 10.0`; the screen calls the **production**
`mass_consistent_effect` directly, so screen and build cannot drift. `sse*` is the uncorrected
`mse` ratio at its own best global amplitude, i.e. `1 − cos_agg²`.

| arm | HCT116 Δpds | HEK Δpds | HCT Δcos | HEK Δcos | HCT sse* | HEK sse* |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| delivered (production) | — | — | — | — | 0.999478 | 0.999919 |
| `intended` (INFEASIBLE) | +0.00098 | +0.00237 | +13.4% | +61.4% | 0.999339 | 0.999815 |
| **`massfix-minnorm`** | **+0.00035** | **+0.00278** | **+14.6%** | **+67.2%** | **0.999326** | **0.999803** |
| `massfix-const` | −0.00555 | −0.00312 | −4.7% | −18.8% | 0.999523 | 0.999943 |

After the fix the read-back is the identity to **1e-12**, for a correction of **1.76% / 2.29%**
of the row norm — against the renderer's implicit **3.96% / 4.53%**.

### Three internal controls, and they are the reason to believe this

1. **`const` satisfies the identical scalar constraint and LOSES** 0.0056/0.0031 pds and
   4.7%/18.8% cosine. Same constraint, different direction, opposite sign. So this is geometry,
   not normalisation bookkeeping. It is kept in the shipped code as that control.
2. **`minnorm` beats the infeasible `intended` upper bound on HEK293T** (+0.00278 vs +0.00237,
   cosine 0.01311 vs 0.01265). The feasible object is *closer to the truth* than the one we
   designed — which is what you expect if the truth's own manifold is the constraint set.
3. **The placement does not matter beyond "not uniform" (§5).** Seven correction directions all
   land at +14% / +65% cosine. A tuned parameter would not behave like that.

## 4. All six members (AGENTS.md 11.1 rule 12)

DE side is deterministic on the HCT116 cached table, `pred_lfc` recomputed from the delivered
profile. Our emitter gives all 400 cells the same library, so a gene's per-cell CPM mean is
exactly `1e6·p_g` and the arm-to-arm lfc difference is exactly `log2(p_new/p_old)` — no proxy.

| | production | `massfix-minnorm` | `massfix-const` |
| --- | ---: | ---: | ---: |
| sign accuracy (reference DE pool) | 0.5311 | **0.5417** (+0.0106) | 0.5372 |
| head purity (top decile) | 0.6751 | **0.6844** (+0.0093) | 0.6717 |
| reach (keyed on \|effect\|) | 0.09069 | **0.09668** (+0.00598) | 0.09047 |
| nmae (lower better) | 0.9501 | **0.9479** (−0.0023) | 0.9505 |
| sign flips vs production | — | **318** | 165 |
| …of which the NEW sign is correct | — | **68.9%** | 70.9% |

**The 68.9% is the sharpest number in the experiment.** The fix moves the predicted sign on 318
(perturbation, gene) pairs in the reference's DE pool, and more than two thirds land on the
*correct* side against 0.500 chance. It is removing a systematic bias, not reshuffling noise.

`fid` and `jac` are not counted: `fid`'s cohort is our call set rather than the reference pool,
and both carry `n_real` (rule 9). `reach` here is keyed on `|effect|`, not the shipped p-value
order, so it prices the ordering's sign content only. Both flagged rather than laundered.

## 5. Where to put the correction — flat, and that is the finding (`weighted.py`)

| weight `u` | worst-surface Δpds | Δcos HCT | Δcos HEK |
| --- | ---: | ---: | ---: |
| `1/(|e|+0.001)` | +0.00075 | +13.7% | +64.0% |
| `|e|+0.01` | +0.00069 | +14.9% | +69.9% |
| **`1` (minnorm, shipped)** | **+0.00035** | **+14.6%** | **+67.2%** |
| `w` | −0.00052 | +14.9% | +64.7% |
| `w²` | −0.00029 | +14.7% | +62.1% |
| `1/w` (`const`) | −0.00555 | −4.7% | −18.8% |

Every non-`const` member gives the same cosine gain to within a point, and the pds spread
(0.0013) is smaller than the between-surface disagreement. **There is no placement to tune.**
`minnorm` is kept because it is the parameter-free, principled member. Do not spend a slot on
another member of this family.

## 6. What this closes on the competitors' "anti-share" family

Five teams report gains from removing the shared component of the 300 predicted rows
(AgingLab_BI `L ← L − mean_t(L)`, pds 0.7573→0.7961 through the real scorer; Heisenberg rank-1
deflation, 0.6518→0.6829; Zhou Li's row-Gram Tikhonov inverse; Cameron S; Lingkuan Meng).

**Measured here, that lever is spent for us:** `centroid` alone, pre-render, is worth exactly
**+0.00000** on HCT116 and +0.00098 on HEK293T, and an independent audit of our panel matrix
puts its shared energy at 0.35% and mean pairwise cosine at +0.0009 — against the 0.0747
Heisenberg's operator was fixing. Our donors are already mean-removed by
`Donor.decompose()` and production runs `beta = 0`.

**And it explains a competitor's published board failure.** Duy reports cross-target centring
measured at "+0.041 scaled perturbation discrimination post-hoc on the pseudobulk delta" but "a
seed-matched local render against an identical uncentred control returns only −0.002". That is
this experiment's `centroid only (pre-render) = +0.00000` exactly: **the renderer undoes a
pre-render geometry transform.** Geometry applied *after* the mass constraint is imposed does
stack (`massfix + centroid d1.0`, worst-surface +0.00237), which is the honest single-variable
follow-up if a23 confirms — held deliberately, since it trades cosine on HCT116 and would be a
two-variable arm.

## 7. Caveats, stated rather than buried

- Both frozen baselines were built with **leakage-safe** operators (`held_orion_hct116` sha
  `7cecbf7f`, `held_orion_hek293t` sha `97e486e0`), not production's `union_pruned`. Per the
  standing caveat in `docs/LOG.md`, a held-safe operator and the shipped one are different
  objects; the *mechanism* here is a property of `expm1`/`log1p` and the gate and is
  operator-independent, but the magnitudes are extrapolations.
- Absolute cosine on these surfaces (0.008–0.025) is far below the board-implied 0.072 because
  the harness truth carries 150 cells per target against the challenge's 400. Only **relative**
  changes are claimed to transfer.
- The two surfaces disagree by 4× on the size of the cosine gain (+13.6% vs +56.0%). The
  prediction uses the worst surface.

## Files

`truth.py` (exact-space pseudobulk truth, both surfaces) · `screen.py` (render/invert + the
geometry family) · `massfix.py` (the mass constraint; calls production) · `de_price.py` (the DE
members) · `weighted.py` (§5) · `A23-PREDICTION.md`.

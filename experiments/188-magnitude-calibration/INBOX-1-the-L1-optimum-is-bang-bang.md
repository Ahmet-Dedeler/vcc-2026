# INBOX 1 for 188 — what the L1 optimum actually looks like, and where your value really is

Coordinator (lead), 2026-09-03 01:35 UTC. Reproduce, do not trust. Two corrections to the mandate
and one reframing that changes what you should be optimising.

## 1. The exact `nmae` functional, from the spec (§7 of `docs/reference/vcc2026-metrics-BRIEF.md`)

```
NMAE_p = Σ_{g∈S_p} |lfc_hat − lfc| / Σ_{g∈S_p} |lfc|
```
`S_p` is the **reference's** significant set, target gene removed, gate ≥ 10 genes. So the loss is a
**|t|-weighted L1 on log2 fold change**, and predicting zero everywhere gives exactly 1.000.

Substituting `p̂ = k·t` for correct signs and `−k·t` for wrong ones recovers the 02:15 identity
`nmae = q|k−1| + (1−q)(k+1)` **with `q` the |t|-WEIGHTED sign accuracy**, not the plain one. Use the
weighted `q` everywhere; the two differ and the identity is only exact for the weighted one.

## 2. The mandate's "fit the conditional median" is right, but the shape it implies is NOT a smooth
   shrinkage — it is closer to bang-bang, and that changes the estimator you should fit

Under L1 the optimal prediction for gene `g` is the **conditional median of `t_g` given everything we
observe**. `t_g` is a signed scalar and our knowledge of its sign is nearly a coin flip: with a
per-gene correct-sign probability `q_g`, the conditional distribution of `t_g` is a two-sided mixture
and its median sits at the `(1 − 1/(2q_g))` quantile of `|t|`. At `q_g = 0.537` that is the
**0.069 quantile** — i.e. the L1-optimal magnitude is *tiny*, and it only becomes the full magnitude
as `q_g` rises well above 0.5.

**Consequence you must handle:** a naive isotonic fit of `median(|t|) ~ |h|` on *unsigned* magnitudes
will not find this, because it throws away the sign uncertainty that drives the answer. **Fit the
conditional median of the SIGNED `t` given the signed prediction**, or equivalently fit `q_g` as a
function of observables and then map `(q_g, |h|) → the quantile`. Report both and say which you used.

## 3. The reframing: `nmae` is nearly closed by arithmetic, so your value is in `pds`

From the identity, the best any magnitude map can do on `nmae` at our `q` is `2(1 − q) = 0.926`
raw, against our board 0.944. That is **+0.030 scaled ≈ +0.005 overall, and that is the ceiling for
the entire magnitude family on this member** — it is not a target, it is a budget. Three of the six
members are `q`-bound by formula and `nmae` is the most tightly bound of them.

**So do not spend the run maximising `nmae`. Your job is to harvest 169's `pds` while holding `nmae`
at or better than neutral.** 169's arm is worth **+0.019452 scaled `pds` = +0.00324 overall** on the
board-validated instrument, positive on all four PLOLO folds including the leak-free one. Nothing
else in this project's history is close. It failed by **0.001** on the composed gate, and only
because `nmae` went the wrong way.

## 4. Which implies a specific shape, and it is not the conditional median either

169 §6.2 measures where the `nmae` damage lands: on h1-test **9% of gated genes already have
`|h| > |r|`** and those carry **22.6% of the L1 error**. The identity charges them slope **+1** while
paying only `−(2q−1) = −0.074` on the other 91%. A power law is expansive **everywhere**, so it buys
`pds` contrast in the middle of the magnitude distribution and pays a 13:1 penalty at the top.

**The contrast that `pds` rewards comes from the middle of the distribution; the `nmae` penalty comes
from the top.** They are separable. So the shape to test is an **S-curve**: expansive over the bulk,
flat or compressive above the point where `P(|h| > |r|)` becomes material — a winsorised sharpening
rather than a monomial. Add it as an explicit arm alongside the fitted conditional median, with the
knee located from the donor folds (the quantile of `|h|` at which `P(|h| > |r|)` crosses, say, 0.25),
not tuned on the endpoint.

Concretely, three arms worth pricing, in increasing order of assumption:
- **A. winsorised power** — `|h|^q` for `|h|` below the knee, identity above it. Nearest neighbour of
  169's arm; isolates whether the top tail is the entire `nmae` cost.
- **B. fitted conditional median** of signed `t` given signed `h`, monotone, LOLO. The principled one.
- **C. A ∘ B** only if both are individually clean.

Arm A is the cheap decisive test: if it keeps most of 169's `+0.0195` `pds` and cuts the `nmae` cost
by the 22.6% share the saturated genes carry, the arm clears the gate and the family reopens. If it
keeps the `pds` but not the `nmae` relief, then the `pds` gain lives in the same tail as the cost and
the family is closed for a reason, which is also worth knowing.

## 5. One number to carry, and one warning

`nmae`'s span: 169 §7 used the conservative 0.60 (AGENTS.md range 0.57–0.63) while "our own two board
points solve to 0.6674, which would shrink every `nmae` entry by 10%". The board-fitted slope over
612 entries is `dScaled/dRaw = −1.6532`, i.e. **span 0.6049** — use that, and say what it changes.

**Warning (audit 181 rule 1):** a13 is this family's only board point and it lost **−0.0203 overall**
with raw `nmae` 1.0185, and **both harnesses were off by 0.075 raw and in the wrong direction on it**.
If a13 is not reconstructible, say so and label your `nmae` number unvalidated rather than working
around it.

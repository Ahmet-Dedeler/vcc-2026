# a12 audited on the shipped artifact — `claim` DID reach the cap; the miss was the model's own arithmetic

Measured on `submissions/a12-tier8000-cap.vcc` itself — the exact file that produced the board
number — against `data/challenge/context_{A,B,C}.h5ad`. Nothing here is a re-emission of a
policy. `a11`, `a13` and the previously-cached `a8`/`a10` were extracted the same way, so four
shipped artifacts now have their computable terms measured on the real panel.

---

## Headline, in one paragraph

**`claim` reached 99.75 % of the `#247` cap on the built file** — the cap binds on 720 of 900
rows and the whole shortfall to the cap is 0.25 %. So the emission did what it was designed to
do, and **the shortfall is not in `claim`.** It is in how the model *evaluated* the identity:
`MSE − 1 = (A − 2X + Kc − ρ·claim)/ΣD` was scored by applying `min(claim, B)` to the
**context means**, while the scorer forms `#348`'s budget **inside each context's own scoring
run**. At a12's raised dispersion, context C's `claim` (0.265) overran its own budget (0.2496)
and `ρ_C` fell to **0.9418** — invisible to the mean form, because the mean `claim` (0.245917)
sat just under the mean budget (0.245927). That single term is worth **+0.0218 raw**, against a
total prediction error of **+0.0131**.

With the branch taken per context and `(X, ΣD)` re-solved the same way, **the three arms whose
`claim` is their own measured `Ĉ_p` — a8, a10, a11 — fit a single `(X, ΣD)` to `≤ 0.0001` raw**,
a three-point fit on two unknowns. a12, the one arm whose `claim` is the *reference's* withheld
`C_p`, lands at **−0.0028**, which is inside the `C_p` proxy's own measured uncertainty (below).

```
X = 0.019240        ΣD = 0.229141        deconfounded cosine = 0.07556
```

against the published `0.018846 / 0.228572 / 0.0740` (v11+v12) and the re-anchored
`0.018089 / 0.235508 / 0.0701` that a12 was predicted on. **The re-anchored solve was the
contaminated one**; the correction moves the cosine *up*, not down.

---

## 1. Did `claim` reach the cap? Yes — measured, per row

`mech_a12_audit.py`, on the shipped file. `C_p` is the ntc-construct proxy that produced the
published cap `K = 0.246529`.

| ctx | `Ĉ_p` mean | `Ĉ_p/C_p` mean | p5 | min | cap binds | `claim_c` | `K_c` | % of cap |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 17.110 | 1.176 | 1.141 | 1.091 | **300/300** | 0.235434 | 0.235434 | **100.00 %** |
| B | 14.717 | 0.996 | 0.974 | 0.958 | 120/300 | 0.237291 | 0.239126 | 99.23 % |
| C | 17.810 | 1.088 | 1.055 | 1.039 | **300/300** | 0.265026 | 0.265026 | **100.00 %** |
| mean | 16.546 | 1.087 | | | **720/900** | **0.245917** | 0.246529 | **99.75 %** |

**The earlier assertion — "the tail policy over-disperses past the cap so the cap binds
regardless" — is true in A and C and false in B**, where the emitted dispersion lands at 0.996×
the proxy and 180 rows sit below it. The cost of that is 0.001835 of `claim` in context B,
i.e. **0.77 % of that context's cap and 0.25 % of the panel's**. It is real, it is measured, and
it is not the miss.

Using the artifact's own `A`, `B`, `claim` instead of the grid's assumed values moves the
prediction by **+0.000372** — and it moves it *away* from the board, so `claim` accounts for
**−3 %** of the error.

## 2. Where the miss is: the attribution ladder

`mech_a12_audit.py`, each step swapping exactly one thing:

| raw `mse` | Δ | step |
|---:|---:|---|
| 0.888247 | | **as registered** (assumed `A = s²A_a8`, `claim` at the cap, `min` of the context MEANS, a8+a11 solve) |
| 0.888619 | +0.000372 | + `A`, `B`, `claim` measured on the shipped artifact |
| **0.910448** | **+0.021830** | + **`#348`'s budget evaluated PER CONTEXT** — the model's error |
| 0.900182 | −0.010266 | + `(X, ΣD)` re-solved with the same correction applied to the anchor points |
| **0.901349** | **+0.001167** | **BOARD** |

So of the `+0.0131` miss: **`claim` −3 %, the per-context branch +167 %, the re-solve −78 %,
residual +9 %.** The two large terms are the same defect seen twice — once in a12's own
numerator and once in the a11 anchor that `(X, ΣD)` had been fitted to.

**The counterfactual.** Had the branch been taken per context at prediction time — same assumed
terms, same two anchor points, no new information — a12 would have been pre-registered at raw
**0.8960 / `score_mse` 0.0967**, against the board's 0.9013 / 0.0911. Falsifier 2 would still
have fired, but on an error of 0.0053 rather than 0.0131.

### Why the mean form hides it

`min(claim, B)` is concave, so `mean_c min(claim_c, B_c) ≤ min(mean claim, mean B)` always, with
equality only while every context is on the same side of its own cliff. Measured gap:

| arm | `mean_c min(claim_c, B_c)` | `min(mean claim, mean B)` | gap | `ρ` in the worst context |
|---|---:|---:|---:|---:|
| a8 | 0.166762 | 0.166762 | 0.000000 | 1.0000 |
| a10 | 0.169139 | 0.169139 | 0.000000 | 1.0000 |
| a11 | 0.166964 | 0.169060 | **−0.002097** | 0.9775 (ctx C) |
| **a12** | **0.240776** | 0.245917 | **−0.005141** | **0.9418 (ctx C)** |
| a13 | 0.240080 | 0.245925 | −0.005845 | 0.9338 (ctx C) |

The gap is exactly zero on the two arms far from the cliff and grows with every step toward it.
**a12's design — put the panel on the ρ cliff — is precisely the design that makes the two forms
disagree**, so the model was least trustworthy on the arm built to its own recommendation.

### The board now discriminates the aggregation, independently

Three hypotheses, same two free parameters, least squares over all four board points:

| hypothesis | `X` | `ΣD` | cos | max residual |
|---|---:|---:|---:|---:|
| **per-context `min(claim_c, B_c)`** | 0.019086 | 0.229996 | 0.07482 | **0.0014** |
| `min` of the context means (what 096 did) | 0.017425 | 0.239255 | 0.06697 | 0.0068 |
| `ρ = 1` always (one pooled 900-row budget) | 0.017405 | 0.239366 | 0.06688 | 0.0066 |

Per-context fits **4.7× tighter with no extra parameter**. That is a second, independent
confirmation of the per-context reading, which until now rested on the `anchor_version` string
and the We527 two-raws-one-score proof (`MSE-MECHANISM.md` §(d)).

## 3. The three-point re-solve, and a better partition than three-vs-two

Asked for: re-solve on a8, a11, a12. Done — but the honest partition is not "how many points"
but **which points carry `C_p` risk.** a8, a10 and a11 all claim their own measured `Ĉ_p`; the
`#247` cap never binds on them, so their terms are exact. a12 is the first arm in the project
whose `claim` is the *reference's* withheld `C_p`.

| points | `X` | `ΣD` | cos | residuals |
|---|---:|---:|---:|---|
| a8+a11 (the pair a12 was predicted on, corrected) | 0.019242 | 0.229169 | 0.07557 | a10 −0.0001, a12 −0.0028 |
| a10+a11 | 0.019239 | 0.229110 | 0.07556 | a8 +0.0001, a12 −0.0028 |
| a8+a12 | 0.018976 | 0.230635 | 0.07428 | a10 −0.0001, a11 +0.0024 |
| a11+a12 | −0.001129 | −0.098884 | — | **degenerate** (a8 +0.19) |
| LS a8+a11+a12 (**the requested three-point**) | 0.019087 | 0.230010 | 0.07482 | a8 +0.0000, a11 +0.0014, a12 −0.0012 |
| LS all four | 0.019086 | 0.229996 | 0.07482 | ≤ 0.0014 |
| **LS a8+a10+a11 (`C_p`-free)** | **0.019240** | **0.229141** | **0.07556** | **≤ 0.0001** |

**How much do the constants move?** From the re-anchored solve a12 was predicted on
(`0.018089 / 0.235508`) to the `C_p`-free solve: `X` **+6.4 %**, `ΣD` **−2.7 %**, cosine
**+7.8 %** (0.0701 → 0.0756). From the original v11+v12 solve (`0.018846 / 0.228572`): `X`
**+2.1 %**, `ΣD` **+0.2 %**, cosine **+2.1 %**. The original two-point solve was closer to right
than the re-anchoring was, because v11 and v12 both sit far from the cliff and the re-anchoring
pulled in a11, which does not.

**Is three points on a moving family enough?** For `(X, ΣD)`, yes, and the evidence is that
three *cliff-free* points over-determine them and agree to `0.0001`. What three points cannot do
is constrain a **fourth** unknown, and a12 introduced one — the reference's `C_p`. a13
introduces a fifth, `κ` (§6). **The rule that falls out is not "get more points"; it is: a board
point constrains the unknowns that were already in the model when it was scored, and an arm that
introduces a new unknown cannot be predicted from the arms that preceded it.**

`a11+a12` alone returns a negative `ΣD`: the two are 0.16 apart in `s` and both sit on the
cliff, where the curve is nearly linear in `s`, so the pair is ill-conditioned. **Two points are
only informative when they differ in something the identity is stiff in** — a8/a10 differ in
`claim` alone and pin `ΣD` on their own.

## 4. The `C_p` proxy, measured for the first time

`mech_cp_bias.py`. The cap is `min(Ĉ_p, C_p)` with `C_p` the reference's own jackknife, which
096 substitutes with 400 CONTROL cells and calls "a mild UNDER-estimate, conservative in the
direction that matters." That was an argument. `data/harness/h1-test.h5ad` has the challenge's
**exact geometry** — 85 perturbations × 400 cells plus 18,400 controls — so the same substitution
can be made where the real groups are visible.

| surface | depth ctrl/pert | real `C_p` (perturbation groups) | 400-cell control draws | bias | across-pert CV |
|---|---:|---:|---:|---:|---:|
| **h1-test** (400 cells, 20k UMI) | 1.016 | 14.2523 ± 0.7426 | 13.9796 ± 0.2511 | **1.0195** | **0.0521** |
| orion-thinned (150 cells, 12k) | 1.000 | 67.0177 ± 12.19 | 63.2502 ± 4.07 | 1.0596 | 0.1818 |

**The proxy is an under-estimate, as argued — by 1.95 % on the matched-geometry surface.** That
is the wrong direction to explain a12's residual: a higher true `C_p` means a higher cap and
*more* claim. The *second* effect runs the other way and was never modelled: `min(·, C_p)` is
concave, so replacing a spread `C_p` by its mean **always over-states `claim`**.

a12's raw re-priced under each, per-context branch, `(X, ΣD)` from the `C_p`-free solve so the
sensitivity is not circular:

| `C_p` assumption | `claim` | worst `ρ_c` | raw pred | vs board |
|---|---:|---:|---:|---:|
| proxy as used | 0.245917 | 0.9418 | 0.898551 | −0.0028 |
| × 1.0195 (h1's measured level bias) | 0.249462 | 0.9238 | 0.890602 | −0.0107 |
| proxy, CV 0.052 (h1's measured spread) | 0.244545 | 0.9429 | 0.904073 | **+0.0027** |
| both | 0.248254 | 0.9265 | 0.894727 | −0.0066 |

**a12's remaining −0.0028 sits inside that band**, which spans −0.011 to +0.003. The two
corrections are of similar size and opposite sign, and neither is currently pinned on the
challenge panel. **This, not `X`, is the honest error bar on any arm whose `claim` is capped.**

## 5. What a12's tiered emission did NOT do

`mech_a12_profile.py`. a10's profile re-rendered at `s = 1.0108/1.0796` against the shipped a12:
`A` ratio 1.0004–1.0006, `B` ratio 1.0000–1.0002, per-row cosine median 0.99973, min 0.99906.
a11 (amplitude only, emission unchanged) reads a *larger* residual — cosine median 0.99936 — so
the tiered spread did not move the profile more than a plain α rescale does.

⚠️ **That check has no power and must not be quoted as if it did.** The residual is
`‖r‖/‖dĥ‖ = 2.3 %`, and `X` lives in a cosine of 0.076, so a residual of only **1.0 %**
perfectly aligned with the truth would move `X` by **12.8 %**. A whole-vector cosine of 0.9997
is consistent with `X` being anywhere in a ±30 % band. The measurement rules out a gross
profile corruption and nothing finer — exactly the failure `RESEARCH-METHOD` §2 records as
"a validation statistic chosen at the wrong scale will pass."

## 6. The corrected operating point — a12 shipped 2.9 % too hot

With the branch per context, the `ρ` cliff is a **per-context** object and there are three of
them. Under a12's emission (`claim_c` flat in `s`, measured at 0.01 %):

| ctx | `claim_c` | `B_c(s=1)` | `s_c` | `VCC_ALPHA` at that cliff |
|---|---:|---:|---:|---:|
| A | 0.235434 | 0.284448 | **0.9098** | **0.9822** |
| B | 0.237291 | 0.272386 | 0.9334 | 1.0077 |
| C | 0.265026 | 0.284676 | 0.9649 | 1.0417 |
| mean form (what 096 solved) | 0.245917 | 0.280503 | 0.9363 | **1.0109** ← a12 shipped here |

Below every cliff the slope is `2s·mean(A_c − B_c) − 2X`, dominated by `−2X`, so the curve falls.
The instant the **lowest** cliff is crossed, that context's slope becomes `2s·A_c` and the mean
turns up immediately. **The minimum is at `min_c s_c`, not at the mean cliff:**

| `s` | `VCC_ALPHA` | raw `mse` | `score_mse` | overall vs a12 |
|---:|---:|---:|---:|---:|
| **0.9097** | **0.9821** | **0.8813** | **0.1121** | **+0.0035** |
| 0.9363 | 1.0108 (a12 as shipped) | 0.8998 | 0.0927 | +0.0003 |
| 0.9649 | 1.0417 | 0.9394 | 0.0510 | −0.0067 |
| 1.0000 | 1.0796 | 1.0186 | 0.0000 | −0.0152 |

The model's own recommendation of `VCC_ALPHA = 1.0108` was the **mean-form** argmin. The
per-context argmin is **0.982**, worth about **+0.0035 overall on `mse` alone** against a12 as
shipped, less roughly `−0.0001` for a 2.8 % amplitude cut at the α slopes in §8. ⚠️ The cliff
positions move with `C_p`: at h1's measured +1.95 % bias every `s_c` rises by ~1 %, and the
optimum with it. Treat `0.982 ± 0.010` as the band, not `0.982`.

## 7. a13, forecast from the shipped artifact

`mech_a13_forecast.py`, on `submissions/a13-shrink.vcc`. **This settles the check
`A13-PREDICTION.md` explicitly left open** ("the `claim`-stays-at-the-cap check is per-row and
has not landed"):

| ctx | `A_c` | `B_c` | `claim_c` | `Ĉ_p/C_p` | cap binds | `ρ_c` | `num_c` |
|---|---:|---:|---:|---:|---:|---:|---:|
| A | 0.252341 | 0.250336 | 0.235434 | 1.177 | 300/300 | 1.0000 | 0.021991 |
| B | 0.240774 | 0.238729 | 0.237314 | 0.997 | 130/300 | 1.0000 | 0.008688 |
| C | 0.249549 | 0.247493 | 0.265026 | 1.088 | 300/300 | **0.9338** | 0.007976 |
| mean | 0.247555 | 0.245519 | 0.245925 | | 730/900 | | **0.012885** |

**a13's `claim` is at 99.75 % of the cap — identical to a12's — and its margin was not consumed
by the shrinkage** (730 rows at the cap against a12's 720). `--effect-scale 2.4004` also
restored the amplitude almost exactly: `√(A_a13/A_a8) = 0.9354` against a12's 0.9365.

**The consequence is the important part: a13's computable numerator (0.012885) is 0.8 % LARGER
than a12's (0.012782).** Every calibration lever 096 identified is already spent. **a13 is a
pure direction bet**, and its whole `mse` case rides on one quantity no surface we hold measures
on the challenge panel:

```
   kappa = X_a13 / X_prod  ,   and  raw = 1 + (0.012885 - 2*kappa*0.019240)/0.229141
```

| cosine multiplier `c` | `κ` | raw `mse` | `score_mse` | overall (mse leg + rule-1 DE) | vs a12 |
|---|---:|---:|---:|---:|---:|
| 1.000 — no directional gain at all | 0.9354 | 0.8992 | 0.0934 | 0.16547 | **−0.0025** |
| 1.186 — per-gene half, h1-test | 1.1092 | 0.8700 | 0.1241 | 0.17059 | +0.0026 |
| 1.227 — per-gene h1 × per-pert orion (**both rule-1 minima**) | 1.1480 | 0.8634 | 0.1309 | **0.17174** | **+0.0038** |
| 1.357 — both halves at h1 | 1.2689 | 0.8431 | 0.1523 | 0.17530 | +0.0073 |
| 1.490 — what the published 0.8212 requires | 1.3939 | 0.8222 | 0.1744 | 0.17898 | +0.0110 |
| 1.836 — per-gene half, orion | 1.7169 | 0.7679 | 0.2315 | 0.18849 | +0.0205 |

`c` is the cosine multiplier; `κ = 0.9354 · c`. The DE column applies the 097 gate's non-`mse`
net for a13 vs a12 at **rule 1** (`−0.0174` scaled on orion, `reach` dropped as unreadable),
i.e. `−0.0029` overall.

**a13 should now be expected between 0.170 and 0.175 overall, with a real downside below a12.**
The published `0.1866` needs `c = 1.49`, above every multiplier either harness measured for the
per-gene half alone and reachable only if the per-perturbation half contributes another 1.26×.
The published raw `0.8212` is not retired by this — it is exactly the `c = 1.49` row — but it
was quoted as a point estimate when it is the top of a range whose bottom is negative.
Add the `C_p` band of §4 (`−0.011 / +0.003` raw) on top; a13's `claim` is capped too.

## 8. The DE ledger: not biased, and two members are not unmodelled at all

The brief asked whether the DE side is systematically biased in our favour. **On this evidence,
no.**

**(a) Most of a12's favourable DE surprise is inside a measured noise floor.** The 097 gate ran
same-policy re-seeded twins and they differ by `0.0207` / `0.0245` raw on `reach` — **0.023–0.028
scaled**. a12's entire `reach` gain is `+0.0046` scaled: **five times below the floor.** `jac`
(+0.0025) and `fid` (+0.0001) are smaller still. Of the `+0.0107` scaled of favourable DE
surprise on a12, `+0.0066` is `reach` and is unreadable.

**(b) The count is 4-for-6 across two arms, which is nothing.** a11: `reach` +0.0100 favour,
`nmae` +0.0177 favour, `fid` −0.0036 against, `jac` −0.0014 against. a12: `fid`, `reach`, `jac`
all favour, `nmae` −0.0030 against. Binomial 4/6 is `p = 0.34`.

**(c) `nmae` is not "the unmodelled cost" — it is a monotone amplitude response, and the board
has now measured it twice in opposite directions.** The two arms are pure `VCC_ALPHA` moves of
opposite sign, so they are a two-point derivative:

| member | a11 (`Δα = −0.2416`) | a12 (`Δα = +0.1728`) | slope per unit α | |
|---|---:|---:|---:|---|
| `pds` | +0.003328 | −0.002004 | −0.0138 / −0.0116 | **consistent, 19 % apart** |
| `nmae` | +0.002644 | −0.003008 | −0.0109 / −0.0174 | **consistent, 60 % apart** |
| `reach` | −0.011375 | +0.004569 | +0.0471 / +0.0264 | consistent, but both inside the twin floor |
| `jac` | −0.001414 | +0.002494 | +0.0059 / +0.0144 | consistent, both tiny |
| `fid` | −0.003638 | +0.000062 | +0.0151 / +0.0004 | 38× apart — not a slope |

**`pds` and `nmae` both fall as amplitude rises, on both board points, at slopes that agree to
within a factor of 1.6.** 096 modelled both as invariant ("the exact-total repair holds them
under emission changes") — true of *emission*, false of *amplitude*, and a11's prediction of
`nmae −0.0150` had the sign wrong. So the claim that should replace "`nmae` is the unmodelled
cost" is: **`nmae` and `pds` are amplitude-priced members with board-measured slopes of about
`−0.012` and `−0.014` scaled per unit `VCC_ALPHA`, and 096's invariance assumption for them is
retracted.** ⚠️ Two points, one family, and a per-gene shrinkage is not a global α move — use
these to price a12→a12′ amplitude trims, not a13's shrinkage.

## 9. Status of the decomposition — re-registered, not rescued

Falsifier 2's pre-registered consequence was that the point estimates are retired. **They are
retired as they stood, and the corrected model is registered as a new claim rather than as a
patch.** The distinction matters and it is checkable:

* **No free parameter was acquired.** The unknowns are the same two, `X` and `ΣD`. What changed
  is that a concave function is now evaluated where the scorer evaluates it. The corrected form
  has *fewer* approximations, not more.
* **It fits out of sample.** Three `C_p`-free board points over-determine two unknowns and agree
  to `≤ 0.0001` raw. The mean form fits the same four points 4.7× worse with the same two
  parameters.
* **It is not the thing that failed.** The identity (★) was and remains arithmetic, verified to
  `8.9e-16` against `cell_eval2`'s own functions.

**What it may be trusted for, plainly:**

| use | verdict |
|---|---|
| **ranking** arms within the calibration family | **yes** — unchanged, and never was in doubt |
| **predicting the level** of an arm whose terms are MEASURED on the built artifact, per context, with `claim` uncapped | **yes, to ±0.001 raw** — three points, three amplitudes, three claims |
| the same, with `claim` **at the cap** | **±0.010 raw**, set by the `C_p` proxy, not by the model |
| predicting an arm from **assumed** terms (`A = s²A`, `claim = K`) before it is built | **only to ±0.005** — that is what the a12 counterfactual measured |
| predicting an arm that **changes the direction** (a13, any shrinkage) | **no.** `κ` is a new unknown, unmeasured on the challenge panel and bracketed only by two harnesses that disagree 1.19× vs 1.84× |
| `pds`, `nmae`, `fid`, `reach`, `jac` | **no** — 096 never modelled them; §8 replaces the invariance assumption for two of them with a two-point board slope |

**The rule this cost, stated so it transfers:** an aggregate of a **non-linear** per-context
quantity is not the non-linear function of the aggregate. Every term in (★) is linear and
averages cleanly except `min(claim, B)`, and that one term carried the entire error.
Symptomatically, the model was *most* wrong on the arm built to sit exactly where the
non-linearity lives — a partial objective's optimum is the point where its own approximation is
weakest, which is the same shape as `RESEARCH-METHOD` §2's "a single-member objective can point
hardest at the worst arm."

## 10. Measured / derived / assumed

**MEASURED** (on the shipped artifacts + `data/challenge/context_*.h5ad`, through `cell_eval2`'s
own `_across_pert_budget`, `_row_weights`, `jackknife_correction`, `bulk_lognorm_means`):

* `A`, `B`, `Kc`, `Ĉ_p` per row and per context for **a8, a10, a11, a12, a13** — 4,500 rows.
* `claim` and the per-row cap margin under the ntc `C_p` proxy; the 720/900 binding count.
* The `C_p` proxy's own bias and across-perturbation CV on both harnesses (h1-test at the
  challenge's exact geometry).
* Board raws: a8 `1.362666960949` (full precision, `data/leaderboard/lb_20260828T1325.json`),
  a10 `1.3524148`, a11 `0.9099893`, a12 `0.9013492`.
* The profile residual between a12 and a10-re-rendered (§5), and the `X` sensitivity it permits.
* The `α` slopes of all five other members from the two board amplitude moves (§8).

**DERIVED** (arithmetic on the above): `X = 0.019240`, `ΣD = 0.229141`, cosine `0.07556`, all
residual tables, the per-context cliff positions, the corrected operating point, the attribution
ladder, a13's raw as a function of `κ`.

**ASSUMED** (each load-bearing, each falsifiable):

1. **`ΣD_c` and `X_c` are equal across the three contexts.** Unchanged from 096 §(e) and still
   unmeasurable. Under per-context scoring the board's pooled raw is a mean of three ratios; if
   the `ΣD_c` differ, the recovered `ΣD` is an effective weighted mean whose weights differ
   between a `claim` move and an amplitude move. This is the leading suspect for the residual
   ±0.0014 in the four-point fit.
2. **The board's reported `expr_mse_unbiased_capped_norm` is the mean of three per-context
   values.** §2 shows the board prefers this 4.7× over the alternatives, which is evidence, not
   proof.
3. **`X(s) = s·X`** — the direction does not rotate with `VCC_ALPHA`. Supported to the extent
   `A(s)/A(1) = s²` holds (measured, 0.02–0.2 %); §5 shows a cosine check cannot confirm it at
   the scale that matters.
4. **`claim` is flat in `s`** (measured at 0.01 % in 096; not re-measured here).
5. **`C_p` = the ntc-construct control proxy.** Now bracketed rather than trusted: §4.
6. **a13's `κ`.** Assumed only as a range; nothing here estimates it.
7. **The 097 DE gate's non-`mse` net transfers to the board.** Rule 1 is applied (worst surface)
   and `reach` is dropped; the remaining members are still harness numbers.

## 11. Pre-registered, for the next board points

1. **a13's raw `mse` lands in `[0.86, 0.90]`** if the per-gene multiplier transfers at the
   rule-1 rate, and in `[0.895, 0.905]` if it does not transfer at all. **Falsified outside
   `[0.76, 0.91]`**, which is the full `κ` range above plus the `C_p` band. A landing at
   `≈ 0.899` is the informative one: it would say the donor-agreement direction gain is **zero**
   on the challenge panel, and that both harnesses' positive readings are model-surface
   artefacts.
2. **a13's overall lands in `[0.165, 0.176]`.** Falsified above 0.180.
3. **`nmae` and `pds` move against `VCC_ALPHA` on the next amplitude arm**, at
   `−0.012 ± 0.004` and `−0.014 ± 0.005` scaled per unit α. Falsified on a sign flip.
4. **a12 rebuilt at `VCC_ALPHA = 0.982`, emission unchanged, scores raw `0.881 ± 0.010`** and
   `score_mse ≈ +0.112`. Falsified outside `[0.87, 0.90]`. This is the cheapest remaining
   `mse` gain identified and it is one environment variable.
5. **The three `C_p`-free arms continue to fit one `(X, ΣD)` to `≤ 0.002`** as more points
   arrive. Falsified by any uncapped arm missing by more than that — which would put the
   per-context-`ΣD` assumption (1) rather than the model in the dock.

## Files

| file | what it does |
|---|---|
| `mech_a12_audit.py` | per-context terms for every shipped artifact; the branch table; the attribution ladder; every `(X, ΣD)` solve |
| `mech_cp_bias.py` | the `C_p` proxy measured against real perturbation groups on both harnesses; a12 re-priced |
| `mech_a12_profile.py` | did the tiered emission move the profile; and why the cosine check has no power |
| `mech_a13_forecast.py` | a13's measured terms and its raw as a function of `κ`; the two-point α slopes |
| `work/{a11,a12,a13}-stats.npz` | group sums + jackknife streamed from the shipped `.vcc` files |
| `work/{a12_audit,cp_bias,a12_profile,a13_forecast}.log` | the runs behind every number above |

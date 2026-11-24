# 188 — L1-optimal magnitude calibration. **The hypothesis is FALSIFIED, and 169's `nmae` surface question is RESOLVED against the arm.**

**Date:** 2026-09-03 · measurement only. No `cell-eval2`, no build, no `vcc prep`, no submission,
nothing written to `submissions/`, no file outside `experiments/188-magnitude-calibration/`
modified. Peak RSS **7.45 GB** (single process, one-off whitening-operator load), one job at a
time.

---

## Decision, in one line

**NO-GO on every arm in this experiment, and the fitted map is worse than the power law it was
supposed to replace.** The premise — "a power law is the wrong magnitude map; fit the conditional
median and the arm becomes correctly specified" — fails on a measurement that was never taken
before: **our predicted per-gene magnitude carries essentially no information about the true
magnitude.** Over a 36× range of `|h|` the conditional median of `|r|` moves **+9.5 %**
(0.04454 → 0.04879) and per-gene sign accuracy moves **1.8 pp** (0.4996 → 0.5176). The map the
data supports is therefore nearly flat, and a nearly-flat map is a catastrophic `pds` arm
(**−0.1497 scaled**, measured). Every fitted variant — conditional median, budget-constrained
quantile, `|r|`-median, DE-restricted, in both spaces and both parametrisations — is **strictly
dominated** by `sign(h)|h|^q` on the two endpoints that decide the arm.

**Second result, larger than the first.** 169 §11 named the surface question as "worth more than
another exponent sweep": *both* `nmae` harnesses missed a13 by 0.075 raw in the wrong direction,
and that gap was blocking the biggest `pds` lever on the board. **It is now closed, and it closes
against the arm.** 169's 56× disagreement between HCT116 (`+0.00082`) and h1-test (`+0.04640`) is
a **DE-pool-composition artefact**, not instrument bias: `nmae` is paid on the reference's
*significance*-selected DE set, which on h1-test is reproduced **0.7219** by `|lfc|·√cpm`, 0.6406
by `|h|` and only **0.3796** by `|lfc|`, with a **median size of 89** — matching audit 181's
board-implied `n_real` median of 90–110. Priced on a pool built that way, PLOLO — a third,
independent surface, four donor lines, deep-pseudobulk truth — reads the arm at **+0.02392 raw
`nmae`** against h1-test's **+0.04640**: **same sign, within 2×.** The shallow-pool readings agree;
the deep-pool readings (HCT116 at BH 0.35, 767 genes/row; my own first-pass top-2000 pool) all read
5–20× smaller. **The `nmae` branch of 169 §8's gate table that said GO is the branch computed on
the wrong pool.**

**a15 stays exactly as it is. Zero flags change.**

---

## 0. Phase 0 — known-answer rows. **Eleven gates, all passing, several bit-exact.**

| gate | expected | measured | script |
| --- | ---: | ---: | --- |
| identity map reproduces 169's PLOLO a23 rows, chalA, 4 folds | — | **worst \|Δpds\| = 0.00e+00** | `budget188.py`, `frontier188.py` |
| PLOLO Δ`pds`, `h` space, q = 1.15, A/B/C mean | +0.019452 | **+0.019452** | `frontier188.py` |
| PLOLO Δ`pds`, `lfc` space, q = 1.15 | +0.017631 | **+0.017631** | `frontier188.py` |
| a15 → a23 mass-fix PLOLO price | +0.004014 | **not re-measured** — my identity row reproduces the *a23* endpoint of that delta bit-exactly on all four folds; the a15 endpoint was not rebuilt, so this row is inherited from 162, not independently confirmed here | — |
| `claim` flat 0.25/φ0 | 0.11860 | **0.11860** (−0.00 %) | `claim188.py` |
| `claim` flat 1.00/φ0 | 0.19216 | 0.19438 (+1.16 %) | `claim188.py` |
| `claim` flat 0.25/φ4 (a8) | 0.16478 | 0.16652 (+1.06 %) | `claim188.py` |
| **`claim` a15 as shipped** | **0.24581** | **0.24581** (+0.00 %) | `claim188.py` |
| `A`/`B` round trip at s = 1 | 0.283614 / 0.281236 | **exact, +0.000 %** | `claim188.py` |
| **a15 raw `mse` / `score_mse`** | **0.88228 / 0.11112** | **0.88228 / 0.11112** | `claim188.py` |
| 169's Δ`claim` / Δ`mse` overall at q = 1.15 | −0.00001 / +0.000037 | **−0.00001 / +0.000037** | `claim188.py` |
| 169's Δ`claim` / Δ`mse` overall, lfc arm | +0.00010 / +0.000006 | **+0.00010 / +0.000006** | `claim188.py` |
| **h1-test emitted Δ`nmae`, `h` space, 5 matched seeds** | **+0.04640** | **+0.04640** | `six188.py` |
| **h1-test emitted Δ`nmae`, `lfc` space, 5 seeds** | **+0.00860** | **+0.00860** | `six188.py` |
| **re-emission of 169's q = 1.15 arm** (`KA_q1.15_mf`, seed 90 000) | bit-identical to `de169.q1.15_mf.npz` | **0.000e+00 on `real_lfc`, `real_padj`, `pred_lfc`, `pred_padj`, `pred_p`** | `realize188.py` |
| identity (I) vs `debench.Surface.nmae()` on 11 tables | < 1e−12 | **0.0e+00 on all 11** | `identity188.py` |

One row that did **not** reproduce on the first attempt and is reported as a correction rather
than dropped: pricing `mse` against the **ungated** v11 effect instead of the gated one charges
the arm `−0.00254` raw `mse` for zeroing the renderer's off-gate residue — **10× the arm itself**.
169 §3 used the gated `q = 1.00` row for exactly this reason. Both rows are in `claim188.log`; the
ungated one is the anchor for 164's `claim = 0.24581`, the gated one is the baseline for every
delta.

---

## 1. Phase 1 — the fitted map, and the estimand the mandate asked for is not the L1 optimum

### 1.1 A correction to the estimand, derived not chosen

The mandate says "the conditional median of the true magnitude given the predicted one". That is
not the L1 optimum. Writing our emitted value as `x·s` with `s = sign(p)` fixed and `x ≥ 0` free,

```
|x·s − r| = |x − r·s|          (|s| = 1)
```

so the loss is `E|x − y|` with **`y = r·sign(p)`** — the truth *projected onto our own sign* — and
the minimiser is `m(p) = max(0, median(y | p))`. The distinction is the whole point: `median(|r|)`
is positive on every gene, so it hands a large magnitude to genes whose sign we get right 50 % of
the time, which is exactly the population 02:15's landscape charges at slope +1. `median(r·sign p)`
encodes sign accuracy automatically: where `P(sign correct) ≤ ½` it returns ≤ 0 and says *emit
nothing*. Both are fitted; the `|r|`-median is the `map_abs` row throughout and it is the worst
arm in the entire screen (**−0.1497 scaled `pds`**).

Fitted leave-one-line-out on the four PLOLO folds (162's validator, whose truth is the held donor
library's own effect matrix over the real 300-target panel), binned-median + PAVA, 60 quantile
bins, ~2.5 M (target, gene) pairs per fold, on `scored ∩ challenge gate` = 5,767–6,655 columns.
Also fitted: the **L2-budget-constrained** L1 optimum, i.e. the first-order condition
`F(x) = ½ + μ·x` — the conditional `½ + μ|p|` **quantile** of `y`, with `μ` solved by bisection so
the map hits the row norm itself.

### 1.2 The shape, which is the finding

Pooled over all four folds, `h` space, per quantile bin of `|p|`:

| pctile of \|p\| | \|p\| | `m_signed` | `m/x` | `m_abs` = median \|r\| | `m_de` |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.8 | 0.00043 | 0.00000 | 0.000 | 0.04439 | 0.00499 |
| 15.8 | 0.00815 | 0.00020 | 0.025 | 0.04453 | 0.03021 |
| 50.8 | 0.02858 | 0.00045 | 0.016 | 0.04481 | 0.04251 |
| 80.8 | 0.05691 | 0.00114 | 0.020 | 0.04568 | 0.04877 |
| 95.8 | 0.09723 | 0.00261 | 0.027 | 0.04791 | 0.05831 |
| **99.2** | **0.14187** | **0.00722** | **0.051** | **0.05480** | **0.07568** |

**Convex over the bulk in the sense that matters** (`m/x` roughly doubles from the middle to the
head), **hard-zero below the 13th percentile**, and — the number that kills the family —
**`m_abs` is essentially FLAT**: the conditional median of the true magnitude rises **+23 %** over
a **330× range** of our own predicted magnitude (and, on the DE band, **+10 %** across the ten
deciles, 0.10006 → 0.11004). The map crosses the identity (after row-L2 restore) at the **0–20th
percentile** of `|p|` (median over rows: 0.0 / 12.0 / 13.3 / 19.9 by fold), and moves **20–25 %**
of the gated L1 mass up in `h` space (27–34 % in `lfc` space) and the rest down.

Directly measured, per decile of `|p|`, pooled over four folds:

| decile | \|p\| mid | `P(sign ok)` | median \|r\| | median `r·s` | \|p\|/\|r\| med | `P(sign ok \| DE)` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.00255 | 0.4996 | 0.04454 | −0.00006 | 0.051 | 0.5010 |
| 5 | 0.02473 | 0.5029 | 0.04473 | +0.00044 | 0.553 | 0.5049 |
| 9 | 0.06359 | 0.5084 | 0.04603 | +0.00137 | 1.388 | 0.5155 |
| **10** | **0.09249** | **0.5176** | **0.04879** | **+0.00304** | **2.023** | **0.5358** |

> **`|p|` is a sign-confidence variable, not a magnitude variable.** It moves sign accuracy by
> 1.8 pp and the true magnitude by 9.5 %. That is why `median(r·sign p)` is such a violent
> function of `|p|` (0 → 0.0072) while `median(|r|)` is nearly constant: under independence
> `median(r·s) = Q_{|r|}(1 − 1/(2q))`, so `dτ/dq = 1/(2q²) ≈ 2` — **the estimator amplifies a
> 1 pp error in sign accuracy into a 2 pp quantile-level error near the bottom of a steep
> quantile function.** The conditional-median map is therefore an ill-conditioned readout of a
> nearly-absent signal, and its non-monotonicity in `m/x` (0.025 → 0.014 → 0.051) is that noise.

One more number worth recording because 169 §6.2 flagged the inherited one as unreproducible:
**median `|h|/|r|` over `scored ∩ challenge gate` is 0.4727–0.7953** across the four folds. That is
a *third* population, different again from 098's 0.138 and 169's 0.2825/0.1266; the low values in
the literature are a **selection effect of the DE pool** (large `|r|`), not evidence of heavy
per-gene shrinkage. On the DE band our median ratio is 0.05 in the bottom decile and **2.02 in the
top** — we already overshoot at the head.

### 1.3 The held-out screen — the fitted map is dominated

4 folds × 3 challenge frames, delivered through the production renderer with `VCC_MASS_FIX=minnorm`
(156: an effect-level score measures an object the board never sees). `nmae*` on the **corrected**
pool (§3): the reference's top-89 by `|lfc|·√cpm` per row. All arms preserve the row L2 norm
exactly (`rA = 1.00000`), so `mse` is neutral by construction.

| arm | Δ`pds` raw | Δ`pds` scaled | Δ`nmae*` | `dU/D` | `e` | `pds` ovr | `nmae` ovr | SUM(2) | keep-top ov |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| identity | 0 | 0 | 0 | 0 | — | 0 | 0 | 0 | 1.00000 |
| **q1.15** | +0.008794 | **+0.019452** | +0.02392 | +0.01709 | −0.20 | +0.003278 | −0.006643 | −0.003366 | 1.00000 |
| **q1.15L** | +0.007971 | **+0.017631** | +0.01766 | +0.00883 | −0.50 | +0.002971 | −0.004905 | −0.001934 | 0.99889 |
| q1.25 | +0.014143 | +0.031284 | +0.05081 | +0.03819 | −0.17 | +0.005271 | −0.014113 | −0.008842 | 1.00000 |
| `map_median` | −0.001537 | −0.003400 | −0.00181 | −0.03230 | 0.47 | −0.000573 | +0.000502 | −0.000071 | 0.96611 |
| `map_abs` | −0.067699 | **−0.149750** | +0.00247 | +0.00594 | 0.29 | −0.025233 | −0.000687 | −0.025919 | 0.98544 |
| `map_de` | −0.050539 | −0.111792 | −0.00197 | −0.00253 | 0.11 | −0.018837 | +0.000549 | −0.018288 | 0.97093 |
| `map_budget` | +0.006424 | +0.014211 | +0.00252 | −0.03050 | 0.54 | +0.002395 | −0.000701 | +0.001694 | 1.00000 |
| `map_debudget` | −0.053371 | −0.118056 | −0.00132 | +0.00089 | 1.24 | −0.019892 | +0.000366 | −0.019527 | 0.98505 |
| `nmap_median` (scale-free) | −0.007354 | −0.016267 | +0.00247 | −0.04409 | 0.53 | −0.002741 | −0.000686 | −0.003427 | 0.95765 |
| `nmap_budget` (scale-free) | +0.002425 | +0.005365 | +0.01249 | −0.06989 | 0.59 | +0.000904 | −0.003469 | −0.002565 | 1.00000 |

Full table (32 arms, both spaces, both parametrisations) in `work/frontier188.log`.

**`map_budget` is the one row that looks like a win, and it is not one.** Its `pds` is
**negative on `orion_hek293t`** (+0.01267 / +0.00917 / **−0.00399** / +0.00784 across folds) where
every power arm is positive on all four; and its shape is smooth `m/x` 0.43 → 0.65 across the whole
bulk with a **single top-bin jump to 2.11**, i.e. the τ → 0.9995 saturation of the first-order
condition rather than a shape. Its scale-free twin `nmap_budget`, the only version transferable to
another surface, reads **+0.07453 raw `nmae`** on h1-test — the worst arm measured anywhere in this
experiment. `map_budget` is also **not priceable on `mse`** (096's cached object is v11's
*delivered* effect at a different amplitude; a power arm with the row norm restored is exactly
invariant to that rescale, a raw-axis map is not), so under rule 12 it cannot be shipped even if it
were consistent. In the `lfc` space the budget solve is degenerate on 3 of 4 folds (bisection hits
the τ cap, `μ → 1.7e7`) and its output is the cap map, not the first-order condition.

### 1.4 The extra covariate — control expression earns a little, donor agreement is board-refuted

`tilt_θ: h_g → h_g·(cpm_g/median)^(−θ)`, row L2 restored, held out fold by fold:

| arm | Δ`pds` scaled | per fold | Δ`nmae*` | SUM(2) | keep-top ov |
| --- | ---: | --- | ---: | ---: | ---: |
| `tilt−0.10` | +0.004333 | +0.00102 / +0.00272 / +0.00262 / +0.00149 | −0.00049 | **+0.000866** | 0.99671 |
| `tilt+0.20` | −0.013386 | — | +0.00739 | −0.004308 | 0.99360 |
| **`q1.15_tilt−0.10`** | **+0.023373** | +0.01128 / +0.01215 / +0.01387 / +0.00496 | +0.02421 | −0.002787 | 0.99713 |

**It earns its complexity on `pds` and nowhere near enough to matter.** A −0.10 expression tilt on
top of q = 1.15 adds **+0.0039 scaled `pds`** (+20 % of the arm) at *equal* `nmae*` cost, positive on
all four held-out lines. `tilt−0.10` alone is the only arm in the screen with a positive two-member
sum and a consistent sign, at **+0.00087** — less than half the shape-only gate, before
`reach`/`fid`/`jac` are even measured, and it breaks the `|effect|` ordering (keep-top overlap
0.9967) so it is a two-variable arm. Recorded as a live-but-sub-threshold lead, not an arm.

**Donor agreement is the other covariate the mandate names, and the board has already tested it:
that is a13.** It lost **−0.0203 overall** with raw `nmae` **1.0185**. It is *not* a clean known
answer for this family — `docs/LOG.md` records a13 as `--effect-scale 1.4025` **and** `gamma 0.5`
**and** the per-gene agreement shrinkage, i.e. an amplitude change *and* an operator change *and*
the magnitude map, and its raw `mse` moved 0.9013 → 0.8493. 162 records it as not reconstructible;
**it still is not**, and reconstructing it would need per-donor effect tables that `build_fold` does
not cache. Stated as a limitation (§7), not worked around.

---

## 2. Phase 2 — the six-member price, on the instrument each member is validated on

Rule 12. `pds` from PLOLO × 1.011/6; `mse` from 096's `claim` route at the `rX = 1` null;
`reach`/`fid`/`jac`/`nmae` from h1-test emitted tables with **five matched seeds** and 169 §4's
conventions verbatim (rule 9 handled: `fid` read as precision, `jac` rebuilt on the board's
`n_real`/`n_pred`, `reach` transferred by 157's a8↔a15 factor 0.3930).

### 2.1 The arm 169 wants resurrected — `h` space, q = 1.15

| member | instrument | Δ overall |
| --- | --- | ---: |
| `pds` | PLOLO, 4 folds × 3 frames | **+0.003278** |
| `mse` | 096 `claim`, `rX = 1` | +0.000037 |
| `reach` | h1-test, 5 matched seeds (+0.00848 raw, sem 0.00387) | +0.000627 |
| `fid` | h1-test precision | −0.000032 |
| `jac` | h1-test recall → board counts | −0.000419 |
| **`nmae` — h1-test, emitted, 5 seeds (+0.04640 raw)** | board-validated pool shape | **−0.012862** |
| `nmae` — PLOLO, wald-89 pool (+0.02392 raw) | pool-matched, independent | −0.006643 |
| `nmae` — HCT116 @ BH 0.35 (+0.00082 raw, 169) | deep pool, disqualified (§3) | −0.000229 |
| **SUM, h1-test `nmae`** | | **−0.009371** |
| **SUM, PLOLO `nmae`** | | **−0.003152** |
| SUM, HCT116 `nmae` | | +0.003262 |

### 2.2 `lfc` space, q = 1.15 — 169's best candidate, which missed by 0.001

| member | Δ overall |
| --- | ---: |
| `pds` | +0.002971 |
| `mse` | +0.000006 |
| `reach` | +0.000633 |
| `fid` | +0.000002 |
| `jac` | −0.000270 |
| **`nmae` h1-test (+0.00860 raw)** | **−0.002383** |
| `nmae` PLOLO wald-89 (+0.01766 raw) | −0.004905 |
| `nmae` HCT116 (+0.00054 raw) | −0.000149 |
| **SUM, h1-test** | **+0.000959** |
| **SUM, PLOLO** | **−0.001563** |
| SUM, HCT116 | +0.003193 |

169 reported `+0.000930` on the h1-test branch; I reproduce `+0.000959`. **The third surface moves
it below zero.** The arm no longer "misses by 0.001" — under the pool-matched reading it is
negative.

### 2.3 The fitted maps

| arm | `pds` ovr | `mse` ovr | `reach` ovr | `fid` ovr | `jac` ovr | `nmae` ovr (h1) | SUM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `nmap_median` | −0.002741 | +0.000023 | +0.001009 | −0.000224 | −0.001656 | −0.009691 | **−0.013280** |
| `nmap_budget` | +0.000904 | +0.000116 | +0.001220 | −0.000555 | −0.002543 | −0.020657 | **−0.021516** |

One seed each on the DE members; the arms are 3–10× outside the shape-only gate so a seed ladder
would not change the verdict. `nmap_median` also carries `rA = 1.0174` on 096's object — a 1.7 %
amplitude ride-along, i.e. a second variable, reported rather than hidden.

### 2.4 `nmae` from the identity, and how it compares to the surfaces

169 §5's identity, re-gated here against `debench` to **0.0e+00 on all eleven tables**:

```
nmae_p = 1 + (U_p − 2 M_p)/D_p ,  U = Σ|h| , M = Σ min(|h|,|r|)·1[sign ok] , D = Σ|r|
```

so a shape arm is priced by exactly two numbers and the marginal refund rate `e = dM/dU`;
**`nmae` improves iff `e > ½`, and that threshold is not fitted.**

| surface / pool | arm | `dU/D` | `dM/D` | `e` | Δ`nmae` |
| --- | --- | ---: | ---: | ---: | ---: |
| h1-test, reference BH<0.05 (median 89) | q1.15 | **+0.03828** | −0.00402 | −0.11 | +0.04633 |
| h1-test, same | q1.15L | **+0.00221** | −0.00317 | −1.44 | +0.00855 |
| PLOLO, wald-89 | q1.15 | **+0.01709** | −0.00341 | −0.20 | +0.02392 |
| PLOLO, wald-89 | q1.15L | **+0.00883** | −0.00441 | −0.50 | +0.01766 |
| PLOLO, top-2000 \|lfc\| *(wrong pool)* | q1.15 | **−0.01465** | −0.00879 | +0.60 | +0.00293 |
| PLOLO, top-2000 \|h\| | q1.15 | −0.01521 | −0.00930 | +0.61 | +0.00339 |

**The identity and the surfaces do not disagree — they never did.** All six rows satisfy
`Δnmae = dU − 2dM` exactly by construction, and every apparent surface disagreement is a
disagreement about `dU`, which is a property of *which genes are in the pool*, not of the arm.
`e < ½` on **every** row: even where sharpening *removes* L1 from the pool it removes refund faster.
**Within the shape family `nmae` opposes concentration structurally.** The `dU` term is also the
whole explanation of the 5.4× gap between the `h`-space and `lfc`-space arms on h1-test: +0.0383 vs
+0.0022, a 17× difference, with near-identical `dM`.

**Mechanism, stated once because it governs the whole family:** the reference's DE pool is
*expression-biased* (a Wilcoxon has more power on well-expressed genes — measured: the pool's median
control CPM is **76.99** against the gate's **47.68**, 1.6×), and the builder's effect space
`log1p(TS·p)` is *expression-weighted*. So any transform that concentrates magnitude in `h` dumps
its L1 **exactly onto the scored pool**. That is why `pds` and `mse` gain in that space while
`nmae`, `reach`, `fid` and `jac` are all paid there.

### 2.5 The gate-width hypothesis — tested, **falsified**

169 priced `nmae` on an artifact gated at `--min-cpm 5.0` (10,783/18,080 genes) while production
gates at `--min-cpm 10.0` (9,135/8,668/9,127 of 18,533). Restricting the transform's support to the
sub-gate and passing the 5–10 CPM band through **unchanged** (so only one variable moves):

| arm | transform support | median n_pred | Δ`reach` raw (5 seeds) | Δ`nmae` raw | SUM DE ovr |
| --- | --- | ---: | ---: | ---: | ---: |
| `q1.15_mf` (169) | 10,783 (CPM>5) | 5,820 | +0.00848 ± 0.00387 | **+0.04640** | −0.012687 |
| `g10_q1.15` | 9,598 (CPM>10) | 5,871 | +0.00989 ± 0.00299 | **+0.04541** | −0.012112 |
| `g20_q1.15` | 8,047 (CPM>20) | 5,920 | +0.01029 ± 0.00276 | **+0.04278** | −0.011539 |

**Matching production's own gate width removes 2.1 % of the cost; a 25 %-narrower gate removes
7.8 %.** The h1-test `nmae` reading is a reading of the production arm. This closes the last
available "the harness isn't measuring what we'd ship" objection.

---

## 3. Phase 3 — the `pds` sign contradiction. **RESOLVED: it is the SPACE, localised to an 11 % gene band.**

169 §6.4: the same lfc-space arm reads **+0.017631 scaled on PLOLO** and **−0.02339 raw on the
h1-test DE table**. Three things differ at once — space, decoy pool, gene set. Varied one at a time:

**(a) Decoy pool: not the cause.** Subsampling PLOLO's 269–294 decoys to **85** (40 draws) leaves
every delta positive and in fact *larger*: q1.15 +0.01943 → **+0.02838** scaled, q1.15L +0.01786 →
**+0.02505**.

**(b) Space: the cause, demonstrated on one table.** Converting **both** sides of the same h1-test
DE table into the scorer's `h` space with the harness's own control CPM
(`h = log1p(TS·10⁻⁶·cpm·2^lfc) − log1p(TS·10⁻⁶·cpm)`, exact and invertible) — same 85 decoys, same
10,694 genes, same prediction:

| arm | `pds` lfc | Δ lfc | `pds` **h** | Δ **h** | Δ h scaled |
| --- | ---: | ---: | ---: | ---: | ---: |
| q1.00 | 0.80182 | 0 | 0.85560 | 0 | 0 |
| q1.15 | 0.80658 | +0.00476 | 0.86218 | +0.00658 | **+0.01456** |
| **q1.15L** | 0.77843 | **−0.02339** | **0.86064** | **+0.00504** | **+0.01115** |
| q1.25L | 0.75938 | −0.04244 | 0.86331 | +0.00770 | +0.01704 |

`debench.Surface.pds()` is reproduced exactly (0.80182 / 0.80658 / 0.77843 / …), so the negative
number is real — **and it is a log2-FC-space statistic that has never been validated against a
board `pds` point.** `pds_cosine` is defined on `bulk_lognorm_means` differences, i.e. `h`. In that
space **h1-test and PLOLO agree on sign, on ordering, and within 25 % on size** for both arms.

**(c) Where the negative lives.** Sweeping an expression floor on the same table:

| CPM floor | n genes | Δ q1.15L, lfc | Δ q1.15L, **h** |
| ---: | ---: | ---: | ---: |
| 5 | 10,686 | **−0.04926** | +0.01053 |
| **10** | **9,516** | **+0.01580** | +0.01425 |
| 20 | 7,966 | +0.00496 | +0.01239 |
| 100 | 2,816 | −0.00062 | +0.00031 |

**Dropping the 5–10 CPM band flips the sign.** That band is 11.0 % of genes and carries **21.4 % of
the reference's total `|lfc|` norm** but only **9.9 % of its `|h|` norm** — a >2× loading
difference, and it is where a ~500-cell reference mean is noisiest. `sharpen_lfc` concentrates our
magnitude there by construction (169 §6.4), so the lfc statistic scores us against low-count noise
while the scorer's own statistic weights that band at less than half. The `h`-space reading is
stable at every floor.

**Verdict: PLOLO was right, the contradiction is closed, and it was never about the decoy pool.**
Also closed as a by-product: on PLOLO the lfc-space `pds` does **not** invert (+0.02250 scaled),
because PLOLO's truth in that band is a deep donor pseudobulk rather than a 500-cell mean.

---

## 4. Phase 4 — the composed slot gate

**Gate 0 — eligibility.** Effect matrix only; no emitter change, no count law, no per-cell library
change. PLOLO eligible. **PASS.**

**Gate 2 — amplitude.** Row-wise L2 preserved by construction; `rA` = 0.97869 (q1.15) / 0.97810
(q1.15L) on 096's object, `VCC_ALPHA` untouched. **SHAPE-ONLY → threshold +0.0019.**

**Gate 3 — rule 13's `reach` toll.** Profile class, so the −0.0146 emission toll does **not**
apply (audit 181 amendment; 169 falsified it here and this experiment reproduces the measurement at
+0.00848 ± 0.00387 raw over five matched seeds, and finds the same for the gate-width variants).

**Gate 4 — 162's budget bar, +0.0126 scaled `pds` for a shape arm.** `q1.15` **+0.019452 PASS**;
`q1.15L` **+0.017631 PASS**; `q1.15_tilt−0.10` **+0.023373 PASS**; `map_budget` +0.014211 PASS but
fold-inconsistent; every other fitted map **FAIL**.

**The composite:**

| arm | `nmae` from | predicted Δoverall | vs +0.0019 |
| --- | --- | ---: | --- |
| q1.15 (`h`) | **h1-test (board-validated pool shape)** | **−0.009371** | **FAIL** |
| q1.15 (`h`) | PLOLO wald-89 (pool-matched, independent) | **−0.003152** | **FAIL** |
| q1.15 (`h`) | HCT116 @ BH 0.35 (deep pool — §3 disqualifies it) | +0.003262 | pass |
| **q1.15L** | **h1-test** | **+0.000959** | **FAIL** |
| **q1.15L** | **PLOLO wald-89** | **−0.001563** | **FAIL** |
| q1.15L | HCT116 | +0.003193 | pass |
| `nmap_median` | h1-test | −0.013280 | FAIL |
| `nmap_budget` | h1-test | −0.021516 | FAIL |
| `map_budget` | PLOLO only; `mse`/`reach`/`fid`/`jac` **unpriced** | — | **INELIGIBLE (rule 12)** |
| `tilt−0.10` | PLOLO only, 2 members | +0.000866 | FAIL |

### NO-GO. `a15` stays exactly as it is:

```
VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.10 VCC_MASS_FIX=minnorm \
  --min-cpm 10.0 --overdispersion 4.0 --tail-spread 1.0 --keep-top 8000 \
  --keep-spread 0.25 --whitener-id union_pruned \
  --whitener data/donors/whitener_union_pruned.npz \
  --whitener-metadata data/donors/whitener_union_pruned.json
```

**Zero flags change. There is no recommended arm, no builder invocation, and no preregistration**
— the deliverable's conditional ("if it clears") is not met. Nothing was built.

**What changed against 169.** 169 left the family alive on one branch: "five of the six members
clear comfortably; the sixth flips the answer depending on which of two harnesses prices it." A
third surface now exists, its pool is matched to the scorer's own DE call, and **it lands with
h1-test.** The HCT116 branch of 169 §8 is closed: HCT116's reading is a deep-pool reading
(767 genes/row at BH 0.35) and the deep-pool readings are the ones that invert the sign of `dU`.
The correct summary is no longer "one of two harnesses says no" but **"the two shallow-pool
instruments agree it costs 2–7× the entire `pds` gain, and the one that disagrees is measuring a
different pool."**

---

## 5. What is now known that was not, in five lines

1. **Our per-gene magnitude is a sign-confidence variable, not a magnitude variable.** `median(|r|)`
   moves +9.5 % over a 330× range of `|p|`; sign accuracy moves 0.4996 → 0.5176. Any "calibrate the
   magnitude" arm is therefore a magnitude-*discarding* arm in disguise, which AGENTS.md §12 already
   lists as closed — this is the measurement of *why*.
2. **`nmae`'s DE pool is expression-biased and the builder's space is expression-weighted**, so
   concentrating in `h` dumps L1 exactly on the scored pool (`dU/D` +0.0383 vs the lfc arm's
   +0.0022). Any future `pds` lever should be checked against `dU` before it is built.
3. **The refund rate `e = dM/dU` is below ½ on every surface, every pool and every arm measured**
   (−1.44 … +0.61). `nmae` opposes concentration structurally, not incidentally.
4. **The `nmae` surface disagreement was pool composition.** Shallow (89) → large cost; deep
   (767–2,000) → small cost. `|lfc|·√cpm` reproduces the scorer's own call at 0.7219; `|lfc|` at
   0.3796. **Any future `nmae` screen must state its pool's size and selection key.**
5. **`debench.Surface.pds()` is a log2-FC statistic and must not be used to size a `pds` arm.** On
   h1-test it inverts the lfc-space arm's sign, and the inversion is carried entirely by the 5–10
   CPM band, which holds 21.4 % of the lfc norm and 9.9 % of the `h` norm.

---

## 6. Files

| file | what it produces |
| --- | --- |
| `mapfit.py` | the shared definitions — the derivation of `y = r·sign(p)`, PAVA, the binned-median fit, the L2-budget-constrained quantile fit, and the `h ↔ lfc` conversions |
| `prep188.py` | freezes the four PLOLO (prediction, truth) pairs and the three challenge frames |
| `shape188.py` | §1.2 — the fitted maps, LOO and pooled, both spaces, three estimands, and the shape report |
| `normfit188.py` | the scale-free (`u = |p|/rms`) refit, needed to carry a map to another surface; `maps_norm.json` |
| `budget188.py` | the first held-out screen, frame A, and the decile table of what `|p|` carries |
| `frontier188.py` | §1.3 — the consolidated 32-arm frontier, 4 folds × 3 frames, corrected pool; the exact 169 known-answer rows |
| `tilt188.py` | §1.4 — the expression-tilt family and the L1-norm control |
| `realize188.py` | the h1-test emissions: the bit-identical re-emission of 169's arm, the gate-width ladder, the map arms |
| `six188.py` | §2 — `reach`/`fid`/`jac`/`nmae` in board units, matched-seed ladders |
| `claim188.py` | §0, §2 — `mse` via 096's `claim` route, four known-answer gates plus the gated/ungated baseline correction |
| `identity188.py` | §2.4 — identity (I) on every surface, `U`/`M`/`D`/`e` decomposition |
| `pool188.py` | §2.4, §3 — what kind of set the reference's BH pool is, and the arms re-priced under each construction |
| `phase3.py` / `phase3b.py` | §3 — space vs pool vs gene set, and the expression-floor localisation |
| `de188.*.npz` | the h1-test DE tables (12 arms) |
| `work/*.log`, `work/*.json` | every table above, reproducible by re-running the script named beside it |

---

## 7. What I did NOT measure, and the limitations that matter

| not measured | why | can it move the verdict? |
| --- | --- | --- |
| **a13, retrodicted** | still not reconstructible (162's finding holds): it needs per-donor agreement tables `build_fold` does not cache. **And it would not be a clean row anyway** — `docs/LOG.md` records a13 as `--effect-scale 1.4025` **+** `gamma 0.5` **+** the per-gene agreement map, i.e. an amplitude change and an operator change riding with the magnitude map, with raw `mse` moving 0.9013 → 0.8493. | **This arm class has NO board point at all.** Under audit 181 rule 1 my `nmae` numbers are therefore **unvalidated against the board**, and the verdict rests on two independent local surfaces agreeing after their pools are matched — not on either one's absolute level. Stated as the central limitation. |
| the challenge's own `n_real` distribution | needs the withheld reference. Audit 181's board-implied **median 90–110** is the best available and it matches h1-test's 89. | It is the input that decides which pool depth is right; if the challenge's pool were really the deep one, the HCT116 branch would revive. I could not close this. |
| `map_budget` (raw-axis) on `mse`, `reach`, `fid`, `jac` | 096's cached object is v11's *delivered* effect at a different amplitude, and a raw-axis map is not invariant to that rescale; its scale-free twin is measurable and reads −0.0215 overall. | No — it is fold-inconsistent on `pds` (negative on `orion_hek293t`) and its gain is one saturated quantile bin. |
| `tilt−0.10`'s `reach`/`fid`/`jac` | its two measured members sum to +0.00087, under half the gate, so a 20-minute emission could not change the verdict. | No, but it is the one live lead (§8). |
| donor agreement as the map's covariate | a13 is the board's own test of it and it lost −0.0203; reconstructing it is out of scope here. | Recorded as board-refuted-with-confounds rather than untested. |
| `cell-eval2`, any build, any submission | hard constraints. Never run. | — |

**Three of my own controls failed and are reported as failures.** (i) My first `nmae` pool was the
top-2000 by `|lfc|` — the wrong construction *and* the wrong size; it inverts the sign of `dU` and
made the fitted `budget` map look like a win. It was caught only because it disagreed with h1-test,
which forced the pool measurement in §3. (ii) My first `mse` baseline was the ungated v11 effect,
which charged the arm 10× its own size; caught by failing to reproduce 169's `Δclaim = −0.00001`.
(iii) My first `pds` run applied the frame-A gate in frames B and C, reading 0.8 % low on q1.15;
caught by the A/B/C known-answer row, which then reproduced `+0.019452` exactly.

---

## 8. The one thing worth doing next

Not another map. The measurement that would move the campaign is in §2.4's `dU` term:
**`nmae`, `reach`, `fid` and `jac` are all paid on an expression-biased pool that `|lfc|·√cpm`
reproduces at 0.7219, and `pds` is not.** Every `pds` lever this project has found concentrates
magnitude in a space that loads that pool. The unexplored move is a norm-preserving transform whose
`dU` on a *predicted* pool proxy is **≤ 0 by construction** — i.e. buy `pds` from coordinates the
reference is unlikely to call DE. `tilt−0.10` is the crudest possible version of that idea and it
already reads `dU +0.0052` with `e = 0.55` and a consistent +0.0043 scaled `pds`; a version that
targets the pool proxy directly is the natural next arm, and it belongs with whoever owns the DE
call set (185), not with the magnitude shape.

Two sibling experiments landed while this one ran and both point the same way. **189** finds a
per-gene cross-line reliability prior that is genuinely orthogonal to `|H|` (spearman +0.031) and
genuinely orders sign accuracy (0.5873 → 0.6694 inside the first `|H|` decile of h1-test's DE pool)
and still buys nothing, because `reach`'s 90 % purity bar does not move. **190** closes the
regulatory-prior ordering key with an oracle ceiling 13× below the ship floor. Together with §1.2's
measurement — that our own magnitude carries 1.8 pp of sign information and 9.5 % of magnitude
information — the three results say the same thing from three directions: **the binding constraint
is per-gene directional accuracy on the reference's own DE pool, and no reweighting of what we
already have reaches it.**

Second: **every `nmae` number in this repo's history should be re-read with its pool size and
selection key attached.** 098 §4's gate-depth ladder, 169 §6.1's, and 164's are all the same
statistic evaluated on pools that differ by 20× in depth, and §2.4 shows the depth flips the sign
of the term that dominates the answer.

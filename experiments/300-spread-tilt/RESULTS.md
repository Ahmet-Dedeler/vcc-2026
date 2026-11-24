# 300 — RESULTS: is the whitener's rank-retrieval gain paid by the SPREAD of its training rows?

**Preregistration** `PREREGISTRATION.md`, sha256
`4a8ef74016b0c3202c89dfab840cbfa94c3cb9e4ef40545c583886f1cbf2a381`, hashed 2026-09-06T20:11:40Z into
`work/PREREGISTRATION.sha256` — **before any operator was fitted and before any new fold number existed**. Ended
inside the 3.5-hour box. No build, no `vcc`, no `cell-eval2`, nothing under `submissions/` or `data/`; nothing
written outside `experiments/300-spread-tilt/`.

## 0. Headline

**SPREAD MECHANISM: P1 monotone YES / saturating NO (and the shape-only band ends before the curve turns over) →
FAILS its frozen rule; P2 interaction −0.00023 [−0.00089, +0.00043] → NOT SEPARATED FROM ZERO (mechanism not
confirmed, not refuted); P3 collinear YES (Spearman +0.721 [+0.649, +0.760]) but SMALLER, −0.00072 [−0.00125,
−0.00019] below the matched tilt → FAILS; best shape-only arm `tilt:d1.0` Δ`pds` +0.00094 [+0.00016, +0.00173],
3/3, retention +47.8 %, controls ALL THREE PASS (mirror −0.00105, random seeds +0.00025 / −0.00023), priced +0.00140 pessimistic vs +0.0019
→ DOES NOT CLEAR; FLEX +0.001401 [−0.000420, +0.003361].**

Nothing here is a candidate. The one live operator thread from 299 is now measured as a mechanism and the answer
is a **qualified negative with a sharp boundary**: the reversed-reliability tilt is real, monotone and 3/3
positive at every dose, but

1. it **does not saturate** — Δ`pds` is still rising at the largest *realisable* dose (d = 2.0), so there is no
   plateau inside the band to sit on; and
2. the **shape-only band edge is measured at d ≈ 1.1–1.2** (d = 1.0 delivers +0.233 %, d = 1.25 delivers
   +0.528 % — just outside), so the largest shape-only dose is exactly 299's `wgt:mildrev`. And the arm's price
   crosses the +0.0019 shape-only gate only at a dose that has already left the band: +0.00140 at d = 1.0
   (0.74 × the +0.0019 gate), +0.00180 at d = 1.25 — which would be 0.95 × that gate but is priced against the
   ×4.7 harder +0.0089 one, where it is 0.20 ×. **There is no dose of this arm that clears, and the reason is
   structural: the price and the amplitude rise together, and the amplitude wins.**

And the mechanism named in the mandate — "the metric is paid by the SPREAD of the donor rows over many
directions" — is **only half right as stated**. The direct spread manipulation (arm C) is strongly collinear with
the reliability tilt per target (ρ +0.721) and at the row level the two ranks agree at ρ +0.70…+0.81, yet C
delivers only 23 % of A's Δ`pds`, and the gap is significant. Reliability is **not** merely a proxy for leverage
on the bank's leading 50 directions: it carries something that top-50 leverage does not.

**Ordering noise.** READ-FIRST §3: the digest-identical ordering-noise SD is **0.004502 overall** — three times
299's whole +0.0014 effect and three times the largest number in this file. Every interval below is a paired
target bootstrap on PLOLO folds, not a board prediction.

## 1. Known-answer rows (every one passed; nothing was priced before they did)

| id | claim | tolerance | measured |
|---|---|---|---:|
| KA-0 | `E062.build_baseline` on the rebuilt per-donor rows, float32, == the fold cache's `pre` | 0.0e+00 | **0.0e+00** on 3/3 folds |
| KA-3 | `cd4svd5` through this file's plumbing == 263's row on (`cos_agg`, `a23_cos`, `pds`) | ≤ 1e-9 | **0.0e+00** on 3/3 folds |
| KA-4 | `cd4svd5` per-target `pds` == 263's `pertarget263-<fold>.npz` | 0.0e+00 | **0.0e+00** on 3/3 folds |
| **KA-W1** | A at d = 1.0 per-fold Δ`pds` == 299's `wgt:mildrev` (+0.00032366 / +0.00139999 / +0.00110670) | ≤ 1e-9 | **0.0e+00** on 3/3 folds |
| **KA-W2** | A at d = −1.0 per-fold Δ`pds` == 299's `wgt:mild` (−0.00038377 / −0.00136022 / −0.00140853) | ≤ 1e-9 | **0.0e+00** on 3/3 folds |
| **KA-B** | the reproduced 062 prune's kept rows == the production operator's 7,784 / 10,807 / 10,803 | exact | **7,784 / 10,807 / 10,803** — exact |
| R2 | the BR249 retrodiction, `cd4svd5` over `h1union` on 261's POINT set | 0.0e+00 vs 261 | **0.0e+00** — `compare` max\|d\| 0.0e+00, POINT total \|d\| 0.0e+00 (mine +0.004722217, 261 +0.004722217; board realised +0.0042634) |
| KA-F0 | `DampedOperator` + preserve == `plolo.Operator.apply` on the Flex frame | 0.0e+00 | **0.0e+00** |
| KA-F1 | `flex:prod` == 273's saved `INCUMBENT` (`pds`, own cosine, `cos_agg`) | 0.0e+00 | **0.0e+00** (0.883333 / 0.05206 / 0.07107) |

**KA-W1 and KA-W2 reproduce 299's post-hoc rows bit-for-bit through independently written plumbing** — not to
1e-9, to zero — so every number below is on the same footing as 299's.

**KA-W3 was not run** (see §9): d = 1.8 was queued behind the P1/P2/P3 arms and the box closed first. Its two
neighbours (d = 1.5, d = 2.0) are measured on 3/3 folds and bracket it.

**KA-B is also the residual difference from production.** The row SET is reproduced exactly per bank — the same
`default_rng(0)` 12,000-row cap visiting the donors in production order (donors under the cap consuming no random
numbers) followed by the same common-essential drop. What is *not* the production operator is the bank count: a
PLOLO fold operator is fitted on the fold's TWO banks, the production `whitener_union_pruned` on THREE. On the
Flex frame (§5) all three banks are used, and there the production operator is carried directly as `flex:prod`.

## 2. The design

Every arm is fitted at **the incumbent's exact row count, donor proportions, axis, γ = 0.65 and λ = 0** — all rows
of both fold banks teach every arm — except arm B, whose row set is 062's production prune. The only thing that
changes between an arm and the incumbent is the WEIGHT each training row carries in the covariance, and the weight
is a bounded rank ramp:

    w(d) = 1 + d · (0.5 − q),        q = the row's rank quantile in [0, 1] of the arm's axis

Σ(0.5 − q) = 0 exactly over a rank quantile, so **Σw = n exactly**: same number of rows, same donor mixture, same
coverage, at every dose. The mirror of a dose is d → −d (the same weight distribution reflected about 1). Three
axes: **A `tilt`** — q = u, the rank of the row's Spearman–Brown-corrected split-half reliability (d > 0 tilts
AWAY from the reproducible rows); **B `prune`** — the same ramp on 062's production row set, d = 0 being the prune
alone; **C `spread`** — q = v, the rank of the row's energy share on the top 50 eigenvectors of the fold's own
cached (un-weighted) operator, d > 0 up-weighting LOW-leverage rows, reliability not used.

**One dose of the frozen ladder is not realisable.** w_min = 1 − d/2, so **d > 2 gives negative row weights** and
no valid covariance weighting exists. `tilt:d2.5` (and its mirror, and `prune:d2.5`) failed on all three folds
with `ValueError: array must not contain infs or NaNs` inside the eigendecomposition, in 1.0–1.3 min each. The
ladder therefore tops out at **d = 2.0**, where the least-reliable row's weight is exactly 0 and the most
reliable's is 2. This is a property of the preregistered ramp, not of the machine; it is reported, not repaired.

## 3. The arm × fold table (everything over `cd4svd5`, three cancer folds, PLOLO 162)

| arm | k562_gwps | orion_hct116 | orion_hek293t | mean Δ`pds` | pooled 95 % | 3/3 | `cos_agg` % | Δ`nmae89` | amplitude | LT3O | PESSIMISTIC | gate | verdict |
|---|---:|---:|---:|---:|---|:--:|---:|---:|---:|---:|---:|---:|:--:|
| `tilt:d-1.0` (mirror) | −0.00038 | −0.00136 | −0.00141 | **−0.00105** | [−0.00184, −0.00026] | 0/3 | −2.118 | −0.00021 | +0.304 % | +140.8 % | −0.00533 | +0.0019 | DOES NOT CLEAR |
| `tilt:d0.5` | +0.00012 | +0.00068 | +0.00068 | **+0.00049** | [+0.00009, +0.00089] | 3/3 | +0.940 | +0.00021 | +0.469 % | +46.9 % | +0.00076 | +0.0019 | DOES NOT CLEAR |
| **`tilt:d1.0`** | +0.00032 | +0.00140 | +0.00111 | **+0.00094** | [+0.00016, +0.00173] | **3/3** | +1.901 | +0.00011 | +0.233 % | +47.8 % | **+0.00140** | +0.0019 | DOES NOT CLEAR |
| `tilt:d1.25` *(POST-HOC, §4)* | +0.00023 | +0.00160 | +0.00126 | **+0.00103** | [+0.00005, +0.00204] | 3/3 | +2.426 | +0.00018 | +0.528 % | +39.7 % | +0.00180 | +0.0089 | DOES NOT CLEAR |
| `tilt:d1.5` | +0.00032 | +0.00175 | +0.00138 | **+0.00115** | [−0.00004, +0.00235] | 3/3 | +3.038 | +0.00035 | +0.867 % | +34.2 % | +0.00223 | +0.0089 | DOES NOT CLEAR |
| `tilt:d2.0` | +0.00078 | +0.00207 | +0.00165 | **+0.00150** | [−0.00021, +0.00320] | 3/3 | +4.584 | +0.00088 | +0.751 % | +30.8 % | +0.00296 | +0.0089 | DOES NOT CLEAR |
| `prune:d0.0` (the prune alone) | +0.00093 | −0.00155 | −0.00446 | **−0.00169** | [−0.00457, +0.00118] | 1/3 | −2.915 | −0.00115 | −3.350 % | +198.5 % | −0.00858 | +0.0089 | DOES NOT CLEAR |
| `prune:d1.0` | +0.00097 | −0.00065 | −0.00326 | **−0.00098** | [−0.00398, +0.00205] | 1/3 | −1.478 | −0.00088 | −2.602 % | +283.9 % | −0.00473 | +0.0089 | DOES NOT CLEAR |
| `spread:d1.0` | +0.00038 | +0.00080 | −0.00052 | **+0.00022** | [−0.00068, +0.00113] | 2/3 | +1.410 | +0.00125 | +0.676 % | −166.5 % | +0.00057 | +0.0089 | DOES NOT CLEAR |
| `spread:d2.0` | +0.00125 | +0.00076 | −0.00301 | **−0.00033** | [−0.00249, +0.00181] | 2/3 | +3.286 | +0.00075 | +0.398 % | +501.3 % | +0.00153 | +0.0019 | DOES NOT CLEAR |
| `rand:d1.0:s300001` (control) | +0.00059 | +0.00034 | −0.00017 | **+0.00025** | [−0.00038, +0.00087] | 2/3 | **−0.377** | +0.00132 | +0.221 % | −52.2 % | −0.00128 | +0.0019 | DOES NOT CLEAR |
| `rand:d1.0:s300002` (control) | −0.00008 | +0.00038 | −0.00098 | **−0.00023** | [−0.00085, +0.00040] | 1/3 | **−0.120** | +0.00005 | +0.436 % | +260.2 % | −0.00026 | +0.0019 | DOES NOT CLEAR |
| `tilt:d2.5`, `prune:d2.5`, `tilt:d-2.5` | — | — | — | **not realisable** (negative weights) | — | — | — | — | — | — | — | — | — |

Paired target bootstrap: 265's `boot`/`boot_folds`, B = 20,000, a generator seeded per arm name so an arm's
interval does not move when other arms are added. Price: 261's recipe on the three COHERENT scenarios (the
pessimistic column takes, per leg, the end of its measured band that LOWERS the total for that cell's sign).
Gate: +0.0019 if the delivered amplitude moves ≤ 0.5 % (shape-only), else +0.0089.

**Where the `pds` moves (262's groups, contribution to the fold-mean Δ`pds`, and the 0.02-knee):**

| arm | below-chance (n 139) | marginal (n 527) | saturated (n 187) | knee up / down / net |
|---|---:|---:|---:|---:|
| `tilt:d1.0` | +0.000205 | **+0.000724** | +0.000014 | 3 / 4 / **−1** |
| `tilt:d2.0` | +0.000578 | **+0.000905** | +0.000018 | 12 / 9 / **+3** |
| `spread:d1.0` | +0.000101 | +0.000129 | −0.000010 | 4 / 3 / **+1** |
| `prune:d0.0` | +0.000326 | **−0.001983** | −0.000035 | 11 / 20 / **−9** |

Every arm moves `pds` in the marginal band, which is where 262 says the metric is decided; the tilt's gain and
the prune's loss are both marginal-band effects, and neither is carried by knee crossings.

**All three controls on the primary arm PASS**, and they fail in different ways, which is the point. The
**mirror** (`tilt:d-1.0`, the same weight distribution reflected) is −0.00105 — the arm beats it by +0.00199
against a required +0.00047, and its sign is reversed, so the effect is in the *direction* of the tilt and not in
the fact that the rows were reweighted. The **two random-weight controls** (the same ramp assigned to rows at
random) read +0.00025 (seed 300001) and −0.00023 (seed 300002) — the arm beats them by +0.00069 and +0.00117
against the same +0.00047 — and, more tellingly, both random controls **lose** `cos_agg` (−0.377 % and −0.120 %)
while the tilt gains it (+1.901 %). The two seeds straddle zero, which is what an inert control should do.
Perturbing the covariance with an arbitrary bounded ramp is not what produces the tilt's effect. The primary arm
therefore **passes the 3/3 and all three control conditions and fails only price and retention** — the same two
narrow misses 299 reported, now with the controls 299 never ran.

## 4. The three preregistered predictions

### P1 — dose. Monotone YES, saturating NO, and the band closes first. **FAILS its frozen rule.**

| dose d | Δ`pds` | pooled 95 % | 3/3 | `cos_agg` % | delivered | shape-only? |
|---:|---:|---|:--:|---:|---:|:--:|
| 0.5 | +0.00049 | [+0.00009, +0.00089] | 3/3 | +0.940 | +0.469 % | **yes** |
| 1.0 | +0.00094 | [+0.00016, +0.00173] | 3/3 | +1.901 | +0.233 % | **yes** |
| 1.25 *(POST-HOC)* | +0.00103 | [+0.00005, +0.00204] | 3/3 | +2.426 | **+0.528 %** | **no — the edge** |
| 1.5 | +0.00115 | [−0.00004, +0.00235] | 3/3 | +3.038 | +0.867 % | no |
| 2.0 | +0.00150 | [−0.00021, +0.00320] | 3/3 | +4.584 | +0.751 % | no |
| 2.5 | — | — | — | — | — | **not realisable** |

The curve is **monotone non-decreasing across all five measured doses** and 3/3 folds positive at every one of
them, and it is still rising at d = 2.0 — `cos_agg` rises almost linearly in the dose (+0.94, +1.90, +2.43,
+3.04, +4.58 %) while Δ`pds` rises sub-linearly (+0.00049, +0.00094, +0.00103, +0.00115, +0.00150), i.e. the
*cosine* return is linear in the tilt and the *`pds`* return is compressing, but neither has turned over.
**The mandate's P1 predicted saturation before the amplitude leaves the band; the measurement shows the opposite
order — the amplitude leaves the band first, and nothing saturates.**

The frozen rule needed "the largest shape-only dose's Δ`pds` above d = 1.0's with the paired interval above 0".
The largest shape-only dose **is** d = 1.0, so the second clause cannot be satisfied and **P1 does not pass**.

**The band edge was then measured.** With ~70 minutes of the box left and every preregistered arm complete, one
POST-HOC dose, d = 1.25, was added on all three folds for the single purpose of locating where the amplitude
leaves the 0.5 % band — a quantity the frozen P1 rule itself refers to ("the largest shape-only dose") and which
no preregistered dose pins down. It is labelled POST-HOC everywhere it appears. It delivers **+0.528 %**, i.e.
*just* outside, so **the edge sits at d ≈ 1.1–1.2 and the largest shape-only dose is d = 1.0**, confirmed rather
than inferred. (The delivered amplitude is not otherwise a clean function of the dose — +0.469 % at d = 0.5,
+0.233 % at d = 1.0, +0.867 % at d = 1.5, +0.751 % at d = 2.0; it is a mean over folds of per-fold delivered
norms that differ in sign — so the edge is bracketed, not solved. It is bracketed tightly enough to settle the
question.)

**And d = 1.25 is the number that closes this thread.** Its pessimistic price is **+0.00180** — **0.95 ×** the
+0.0019 shape-only gate, the closest any arm in 299 or 300 has come to clearing anything. It does not get that
gate: at +0.528 % delivered it is amplitude-changing and is judged against +0.0089, where +0.00180 is 0.20 ×.
Inside the band the arm prices +0.00140 (0.74 ×); outside it the larger doses price +0.00223 and +0.00296 against
+0.0089 (0.25 × and 0.33 ×). Price and amplitude are driven by the same quantity — the cosine gain — so **buying
enough price to clear the shape-only gate necessarily buys enough amplitude to lose it.** Retention falls with
the dose as well (+47.8 % → +39.7 % → +34.2 % → +30.8 %): the gain concentrates as it grows. **The dose axis of
this arm is closed, and it is closed structurally, not by a near miss.**

### P2 — the essential-gene prune. Premise CONFIRMED, interaction NOT SEPARATED FROM ZERO.

The premise is measured directly and holds strongly: the rows 062's prune drops are the **most reproducible rows
in the bank** — median reliability +0.2729 / +0.1422 / +0.1474 for the essential rows against +0.0706 / +0.0329 /
+0.0546 for the rest (**2.7–4.3 ×**), mean reliability rank-quantile 0.741 / 0.690 / 0.601 against 0.438 / 0.479 /
0.489 (`work/prune_reliability300.json`). The prune is a coarse version of the same tilt.

| reading | value | 95 % |
|---|---:|---|
| the prune alone, `prune:d0.0` | −0.00169 | [−0.00457, +0.00118] (1/3 folds) |
| the tilt on the un-pruned banks, `tilt:d1.0` | +0.00094 | [+0.00016, +0.00173] (3/3) |
| the tilt on top of the prune, `prune:d1.0` | −0.00098 | [−0.00398, +0.00205] (1/3) |
| **the letter of the rule**, (B at d*) − (A at d*) | **−0.00193** | [−0.00482, +0.00093] |
| **the interaction**, [B(d\*) − B(0)] − [A(d\*) − A(0)] | **−0.00023** | [−0.00089, +0.00043] |

The tilt buys **+0.00071** on top of the prune against **+0.00094** on the un-pruned banks. The point estimate has
the sign the mechanism predicts — tilting gains less once the most reproducible rows are already gone — but the
interval covers zero comfortably, so **P2's mechanism prediction is neither confirmed nor refuted**: at this
effect size the interaction would need roughly an order of magnitude more targets to resolve.

**A second, unexpected result sits in this arm.** 062's essential-gene prune, reproduced exactly and evaluated on
PLOLO, **loses** `pds` (−0.00169, 1/3 folds) and 2.9 % of `cos_agg`, and moves the delivered amplitude −3.35 %.
On 273's Flex instrument the production pruned operator **gains** (+0.0021, §5). Both intervals cross zero
([−0.00457, +0.00118] and [−0.001821, +0.006162]), so this is an instrument disagreement between two
non-significant readings rather than a contradiction — but it does mean R299 §6's "062's prune is itself worth
+0.0021, another row-composition effect in the opposite direction from reliability" **is not supported on the
PLOLO folds and is not significant on the instrument it was read on.** The prune's value to production is not
established by either measurement.

### P3 — the direct spread manipulation. Collinear YES, equal NO. **FAILS its frozen rule.**

| reading | value | 95 % |
|---|---:|---|
| **Spearman(per-target Δ`pds` of C at d = 1, of A at d = 1)**, pooled over 853 targets | **+0.721** | [+0.649, +0.760] |
| C's Δ`pds` at d = 1.0 | +0.00022 | [−0.00068, +0.00113] (2/3 folds) |
| A's Δ`pds` at the matched dose | +0.00094 | [+0.00016, +0.00173] (3/3) |
| **C − A, paired** | **−0.00072** | **[−0.00125, −0.00019]** |

The two arms move the **same targets in the same direction** — ρ +0.721, the tightest correlation in this file —
and at the row level the two ranks they weight by agree at Spearman **+0.8125** (k562_gwps) and **+0.6991**
(orion_hek293t) (`work/uv300-orion_hct116.json`; the top 50 eigenvectors carry 24.4 % of the operator's trace, and
the essential rows sit high on both ranks: mean v 0.750 / 0.589 vs 0.436 / 0.490, mean u 0.741 / 0.601 vs 0.438 /
0.489). So the mandate's collinearity clause is satisfied and the two axes are largely one axis.

But the frozen rule also required C to gain **at least as much** as A, and it does not: C delivers 23 % of A's
Δ`pds` and the paired difference is **entirely below zero**. C is 2/3 folds, not 3/3, and its `cos_agg` gain
(+1.41 %) is 74 % of A's for a quarter of the `pds`. **P3 does not pass.**

C's dose curve is also **not monotone**, unlike A's: `spread:d2.0` reads **−0.00033** [−0.00249, +0.00181]
(2/3 folds) against `spread:d1.0`'s +0.00022, driven by a −0.00301 loss on `orion_hek293t`, even though its
`cos_agg` keeps rising (+1.41 % → +3.29 %). Doubling the leverage tilt buys more cosine and loses `pds`; doubling
the reliability tilt buys both. That is a second, independent way in which the two axes are not the same object.

The reading that survives is sharper than either alternative the mandate offered. It is **not** the case that
"the reliability axis was a proxy" (C would then match A), and it is **not** the case that reliability is acting
on its own (the two are 72 % collinear per target). Up-weighting rows with low leverage on the bank's leading 50
directions reproduces most of the *direction* of the reliability tilt and about a quarter of its *size*:
**spread over the leading directions is part of what reliability measures, but the useful part of the reliability
axis is not captured by leverage on the top 50 eigenvectors.** The tail of the spectrum — where 272 located the
`pds` — is exactly where a top-50 leverage rank is blind.

## 5. The Flex reference (273's union-OFF Flex-400 instrument, 85 targets) WITH an interval

Matched refits on the same three donor banks and the same axis, everything else 273's incumbent verbatim (base +
H1-CD4 deflated k = 5 re-solved at f = 0.34, union OFF, scored with the scorer's own per-target statistic).

| arm | `pds` | own cosine | `cos_agg` | delivered | base ~ prod | Δ`pds` vs `flex:full` | 95 % (paired, 85 targets) | Δ own cos | amplitude |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|
| `flex:prod` (= 273's INCUMBENT, KA-F1 0.0e+00) | 0.883333 | 0.05206 | 0.07107 | 3.8711 | 1.0000 | +0.002101 | **[−0.001821, +0.006162]** | −0.000600 | ×0.9829 |
| `flex:full` (un-pruned refit, the reference) | 0.881232 | 0.05266 | 0.07003 | 3.9383 | 0.9651 | — | — | — | — |
| **`flex:tilt:d1.0`** | **0.882633** | 0.05323 | 0.07169 | 3.9269 | 0.9660 | **+0.001401** | **[−0.000420, +0.003361]** | +0.000577 | ×0.9971 |
| **`flex:spread:d1.0`** | **0.882213** | 0.05314 | 0.07094 | 3.9270 | 0.9672 | **+0.000980** | **[−0.000840, +0.003081]** | +0.000478 | ×0.9971 |

**The Flex instrument agrees with PLOLO in sign and in ordering and resolves neither.** The tilt reads +0.001401
(PLOLO +0.00094), the spread arm +0.000980 (PLOLO +0.00022), the production prune +0.002101 (PLOLO −0.00169) —
and **all three intervals contain zero**, because 85 targets is not enough to resolve a 0.001 `pds` effect. 299
reported the tilt's Flex number without an interval and called it "a sign-and-size agreement, not an independent
confirmation"; with the interval now computed, that caution is quantified: the Flex reading is compatible with
the PLOLO reading and with zero alike. It is not evidence that the tilt transfers.

The one thing the interval does settle is the ordering: on Flex the tilt is above the spread arm, as on PLOLO,
and the gap between them (+0.00042) is far smaller than either interval.

## 6. The mechanism row (272's statistic, fold `orion_hct116`)

272's statistic with an INDEPENDENT noise estimate built from the donors' batch-parity half differences: the
share of eigen-directions on which the noise estimate n_i exceeds the eigenvalue λ_i (s_i < 0.5, and s_i = 0 for
a direction the noise estimate swamps entirely), and the share of the fold's WHITENED PREDICTION ENERGY those
directions carry.

| operator | Δ`pds` | directions s < 0.5 | s = 0 | share of the WHITENED prediction energy on them | raw-energy share | Σn / Σλ | energy-weighted mean s |
|---|---:|---:|---:|---:|---:|---:|---:|
| incumbent (cached fold operator) | 0 | 93.0 % | **76.8 %** | **86.8 %** | 71.2 % | 0.855 | **0.1530** |
| `spread:d1.0` (arm C) | +0.00022 | 93.1 % | 77.2 % | **87.0 %** | 72.4 % | 0.926 | 0.1503 |
| `tilt:d2.0` (arm A, the largest realisable dose) | +0.00150 | 92.6 % | **77.8 %** | **88.1 %** | 74.0 % | 0.972 | **0.1332** |
| *299's `sel:q0.5` (selection TOWARD reliability, for contrast)* | *−0.06748* | *90.6 %* | *78.5 %* | *99.9 %* | *43.3 %* | *0.747* | *0.0007* |

The incumbent row reproduces R272 §0.10 exactly (76.8 % / 86.8 % / Σn / Σλ 0.855 / mean s 0.1530), so the
statistic is the one 272 and 299 computed.

The two arms are read at their own best doses (A at d = 2.0, C at d = 1.0), not at a matched dose, so the
comparison below is between each arm's largest measured effect and the incumbent.

**Both arms move the operator the same way, and by amounts that track their Δ`pds`.** Tilting away from the
reproducible rows puts MORE of the whitened prediction energy on directions an independent noise estimate scores
as noise (86.8 % → 88.1 %, **+1.3 pp**) and lowers the energy-weighted signal share (0.1530 → 0.1332, −0.0198);
the direct spread manipulation moves the same two numbers the same way and about a sixth as far (**+0.2 pp**,
0.1503, −0.0027). The ratio of the two arms' movements on the energy statistic (≈ 6.5 ×) and on the signal share
(≈ 7.3 ×) brackets the ratio of their Δ`pds` (≈ 6.8 ×) — two points, so a consistency check and not a fitted
relation, but it is the check the mechanism predicts, and it is the check that would have failed had the `pds`
come from somewhere else.

This is 272's finding read forward rather than backward. 272 showed that the incumbent's `pds` is *paid in*
directions that are mostly sampling noise; 299 showed that fitting the operator preferentially on the rows that
reproduce moves the statistic the wrong way (99.9 % of the whitened energy on noise-dominated directions, signal
share 0.0007) and destroys retrieval; 300 now shows that **tilting the training rows AWAY from the reproducible
ones moves the statistic further in the paying direction and buys `pds` in proportion.** The mechanism is
coherent across all three experiments and its sign is now established in both directions from the incumbent.
What it does not do is clear a gate: the entire span of this statistic between the incumbent and the largest
realisable dose is worth +0.00150 `pds` — a third of the ordering-noise SD (0.004502).

## 7. Two supporting measurements

**(a) The two axes are one axis at the row level, but only mostly.** Spearman(u, v) = **+0.8125** on `k562_gwps`
and **+0.6991** on `orion_hek293t` — reproducible rows are the rows that concentrate on the operator's leading
directions. This is the cleanest statement of the spread reading available in this experiment, and it is measured
on the operator's own eigenbasis, not inferred.

**(b) 062's prune is a reliability tilt in disguise.** The dropped common-essential rows are 2.7–4.3 × more
reproducible than the rows kept, and sit at rank quantile 0.60–0.74 on both the reliability and the leverage axis.
Whatever the prune does to the operator, it does it along the same axis this experiment doses continuously.

## 8. What this settles, and what it does not

* **Settled — the dose axis of the reversed tilt is closed, and closed for a structural reason.** Monotone, 3/3
  at every one of five doses, no saturation, and the shape-only band edge measured at d ≈ 1.1–1.2. The best
  shape-only member is exactly 299's `wgt:mildrev` at +0.00140 pessimistic = 0.74 × its gate. The next dose up
  (d = 1.25, post-hoc) prices +0.00180 = **0.95 ×** the shape-only gate but delivers +0.528 % and is therefore
  judged against +0.0089, where it is 0.20 ×; the larger doses reach at most 0.33 ×, with retention falling from
  47.8 % to 30.8 %. Price and amplitude are both driven by the cosine gain, so **the dose that would pay for the
  shape-only gate is always a dose that no longer qualifies for it.** No dose of this arm clears, and the miss
  gets worse with dose, not better.
* **Settled — the reliability axis is not merely leverage on the leading directions.** A direct spread
  manipulation on the top-50 eigenbasis is 72 % collinear with the tilt per target and 70–81 % collinear at the
  row level, yet delivers a quarter of the `pds`, with the paired difference entirely below zero. The mandate's
  P3 disjunction ("if P3 holds the reliability axis was a proxy; if the tilt gains and P3 does not, reliability
  itself matters") resolves to the second branch **with the first branch's collinearity intact** — an outcome the
  disjunction did not anticipate, and the more informative one: the axis is spread-like, but the operative part
  of it lives outside the leading 50 directions, which is exactly where 272 located the `pds`.
* **Not settled — the prune interaction.** The interaction has the predicted sign (−0.00023) and an interval that
  covers zero by a factor of four. P2 is unresolved at this sample size, not answered.
* **Corrected — R299 §6's prune side-observation.** "062's essential-gene prune is worth +0.0021 on the Flex
  instrument, ≈ 1.5 × the reliability tilt" is now read with an interval ([−0.001821, +0.006162]) and on a second
  instrument (PLOLO, −0.00169 [−0.00457, +0.00118], 1/3 folds, −3.35 % amplitude). The prune's `pds` value is
  **not established on either instrument** and the two point estimates have opposite signs.
* **Not tested.** Doses above 2.0 (the negative-weight wall IS 2.0); the mirror of every A dose beyond d = −1.0; the mirrors of the C doses; arm C's own
  mechanism at doses other than 1.0 and 2.0; a leverage rank built on a deeper slice of the spectrum than the top
  50 (the obvious next construction, and the one this result points at); any per-target use of either axis, which
  271 and 263 already closed on their own constructions.

## 9. Provenance, cost and declared deviations

Inputs, all read-only: PLOLO's fold caches (byte copies of 263's in `work/plolo/`), 299's per-row reliabilities
(`work/rel/<line>.npz`, **copied and reused, not recomputed**), the three donor banks, 271's saved half rows
(272's noise covariance), 261's and 257's pricing constants, 265's bootstrap, 273's Flex instrument, 062's
essential list. Reserves and R3 were not touched.

Files written (all inside `experiments/300-spread-tilt/`): `PREREGISTRATION.md`, `work/PREREGISTRATION.sha256`,
`run300.py`, `sched300.py`, `price300.py`, `flex300.py`, `mech300.py`, `uv300.py`, `endgame.sh`,
`work/sweep300-<fold>.json`, `work/pertarget300-<fold>.npz`, `work/identity300-<fold>.npz`,
`work/pertarget300-flex.npz`, `work/price300.{json,txt}`, `work/flex300.{json,log}`, `work/mech300-<fold>.{json,log}`,
`work/uv300-<fold>.json`, `work/prune_reliability300.json`, `work/sched300.log`, the run logs.

Cost: one operator eigendecomposition (17,942² or 17,806², `scipy` `eigh(driver='evd', overwrite_a=True)`) per
arm per fold, 7–10 min contended; staging a fold 9–11 s; scoring an arm 1–2 s. Peak `ru_maxrss` per process
10.6–11.1 GiB (that figure includes the memmapped donor-bank pages the Gram pass touches). `memory_pressure`
stayed at **77–85 % free** throughout.

**Declared deviations from the mandate, every one of them:**

1. **`tilt:d2.5` is not realisable** and was not measured (§2) — the preregistered ramp gives negative row weights
   for d > 2. Its mirror and `prune:d2.5` are unrealisable for the same reason.
2. **The full mirror ladder was not run.** Only d = −1.0 (which is both KA-W2 and the primary arm's mirror
   control) was measured, on 3/3 folds. The mirrors of d = 0.5, 1.5, 2.0 and of the C doses were queued behind
   the P1/P2/P3 arms in the frozen priority order of `PREREGISTRATION.md` §7 and the box closed first.
3. **KA-W3 (d = 1.8 == 299's `wgt:steeprev`) was not run**, for the same reason. KA-W1 and KA-W2 both reproduced
   at 0.0e+00, in both signs.
4. **Arm B was measured at two doses (0 and 1.0), not five**; arm C at two (1.0 and 2.0), plus its two mirrors
   dropped. P2 and P3 are read at d* = 1.0, which is the best shape-only A dose and therefore the dose the frozen
   rules name.
5. **The random-weight control was run at both preregistered seeds (300001, 300002) on 3/3 folds**, and the
   mirror control at d = −1.0 on 3/3 folds. All three control conditions were evaluated and all three PASS.
6. **One POST-HOC dose, d = 1.25, was added after the preregistered arms were read** (§4), on 3/3 folds, for the
   single purpose of locating the shape-only band edge that the frozen P1 rule refers to. It is labelled POST-HOC
   in every table it appears in, it was not given its own mirror or random control, and no other dose, arm or
   rule was changed after the hash. `price300.py`'s `A_DOSES` tuple was extended to include it so that the P1
   table prints the full ladder.
7. **Three heavy processes ran concurrently for ~40 min** — the two PLOLO fold processes the mandate allows, plus
   the Flex-instrument process, which stages 273's Flex frame rather than a PLOLO fold. The memory gate held with
   large margin throughout (77–85 % free against the 40 % floor; three processes at ≤ 11.1 GiB against 64 GiB).
   At no point were more than two PLOLO fold processes in flight.
8. **`flex300.py` had a bug on its first run** (an arm-kind `split` index) that silently skipped arm C on the Flex
   frame; it was fixed, the reference refit was recomputed and arm C measured. `flex:prod` reproduces 273's
   INCUMBENT at 0.0e+00 on both runs and `flex:full`/`flex:tilt:d1.0` are identical across them.
9. **The last three cells landed inside the final 20 minutes** (`spread:d2.0` on the third fold, the
   random-weight control on all three, and the mechanism row); they are in the tables above and in
   `work/price300.txt` / `work/mech300-orion_hct116.json`. For ~10 minutes of that window three fold processes
   were briefly in flight (two random-control folds started by the scheduler plus the mechanism process) while
   `memory_pressure` read 72–82 % free; peak `ru_maxrss` was 11.04 GiB.

**No arm cleared the frozen gate. Nothing here justifies a slot, and no build was started, designed or prepared.**

---

**SPREAD MECHANISM: P1 monotone YES, saturating NO, band edge measured at d ≈ 1.1–1.2 → FAILS its frozen rule /
P2 interaction −0.00023 [−0.00089, +0.00043], mechanism NOT SEPARATED FROM ZERO / P3 collinear (ρ +0.721 [+0.649,
+0.760]) but smaller by −0.00072 [−0.00125, −0.00019] → FAILS; best shape-only arm `tilt:d1.0` Δ`pds` +0.00094
[+0.00016, +0.00173], 3/3 YES, retention +47.8 %, controls ALL THREE PASS (mirror −0.00105, random seeds +0.00025 / −0.00023), priced
+0.00140 vs +0.0019 → DOES NOT CLEAR; FLEX +0.001401 [−0.000420, +0.003361].**

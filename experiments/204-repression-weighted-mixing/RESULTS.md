# 204 — the 3.3× repression spread is real, `a_i` is not identified, and equal weighting is the **maximum** of the transfer cosine

**Date:** 2026-09-03 · deterministic (seed `20260903` everywhere) · **no `vcc` command was run, no
`cell-eval2` was run, no submission was built, nothing in `submissions/` was touched, and nothing
outside `experiments/204-repression-weighted-mixing/` was written or modified.** Peak RSS
**4.73 GiB** (measured on the heaviest fold, `work/rss_probe.py`) against a 10 GB cap.

---

## Verdict

**NO-GO on every donor weighting derived from achieved repression, and the negative is stronger
than a null: on the channel that AGENTS.md §8.2 says is the entire remaining gap — direction, i.e.
the transfer cosine — the incumbent equal weighting is the MAXIMUM over all 29 alternatives
tested, and 24 of the 29 are below it at 95%.**

Four things were measured, in this order:

1. **194's repression table reproduces exactly by an independent route, and it is bigger on the
   panel than 194 reported.** A closed form on the shipped pseudobulk tables reproduces 194's
   per-cell measurement to **2e-7** on all three donors. Restricted to the **300-target panel** —
   the population production actually predicts, which nobody had measured — the median residual is
   **0.135 / 0.298 / 0.504**, a **3.74×** spread (194's 3.3× was on its extras-padded row set), and
   the share of panel targets clearing the challenge's own >80% criterion is **62.2% / 17.4% /
   1.5%**.

2. **De-attenuating a donor's TABLE is an exact no-op in production — measured, 4.1e-8 relative.**
   `prepare_main_source` divides every donor by its own median centred row norm, and every term in
   it is degree-1 homogeneous in `delta`. So production *already* performs "de-attenuate then
   average", using the median row norm as its attenuation estimate. The only place `a_i` can act
   is as a **relative weight** between donors, which is what every arm here moves.

3. **`a_i` is not identified on the panel population.** Two routes were run as the mandate asked
   and they disagree: the self-gene residual gives `(1, 0.81, 0.57)`, the 3-donor shared-target
   Gram gives `(1, 1.44, 1.25)`. The disagreement was then **located**, not averaged away: a
   one-factor (rank-1) test on 4–6 lines rejects the single-shared-component model wherever it has
   power (RMS residual **19.7–22.4%** of mean |ρ|, |z| up to 11), its built-in same-line positive
   control fires at **z = +11.2**, and the largest structural excess is the **Orion↔Orion pair**
   (z = +5.1 to +7.6) — exactly the artefact 194 predicted. Every estimator that must use the
   Orion↔Orion product (the 3-donor triple has zero degrees of freedom and cannot avoid it) says
   Orion is the higher-amplitude platform; every estimator that avoids it (a fourth line as
   instrument, or a ≥5-line fit) says the opposite. Across four independent fits `a_rel(HEK293T)`
   spans **0.504 to 1.254**, a 2.5× range with barely overlapping intervals.

4. **It does not matter, because the pipeline's response to the weight vector has almost no
   component along the attenuation direction.** A quadratic fit of Δ`pds` on the two weight
   coordinates reaches **R² = 0.997**, and the K562 share carries **partial R² 0.978** against the
   HCT116:HEK293T split's **0.046**. The two Orion lines are the same platform, the same chemistry
   and the same depth, differing essentially only in achieved repression (0.298 vs 0.504) — and the
   fitted optimum of their split is **0.92**, i.e. equal, against the **1.43** attenuation predicts.
   **On the one axis where the hypothesis is cleanly testable, it is falsified.**

**And the model was not merely un-supported, it was rejected with power.** The minimum-variance
model's own numbers predict optimal weighting should remove **15.4–27.1%** of the consensus noise
variance and raise the transfer cosine by **+8.5% to +16.8%** relative — worth **+0.0025 to
+0.0049 overall through `mse` alone**, which would clear rule 13's gate. Scored on the model's own
prescribed weights the measurement is **−1.50%**, 95% CI **[−2.38%, −0.67%]**. Wrong sign, and the
interval excludes the prediction by 7×.

**The permutation falsifier ties the best real arm, exactly as in 194.** `perm_a_5` — the
attenuation vector with the two Orion lines swapped — reads Δ`pds` **+0.00074** against the
mandate's headline arm `a_self`'s **+0.00019**, and `p(Δ>0)` **0.734** against **0.561**.

**Best priced arm: +0.00017 overall** (`deatten_slope`, a 1.15:0.88:0.98 mix) against rule 13's
shape-only ship floor of **+0.0019** — **11× below**. `a_self` (`w ∝ a_i`) prices at **−0.00036**,
`deatten` (`w ∝ 1/a_i`) at **−0.00137**, `minvar` (`w ∝ a_i/σ_i²`) at **−0.00242**. Nothing is
shippable and no slot is requested.

**AGENTS.md §4.2's equal-donor-weight finding survives the measurement it predates**, and it now
rests on a mechanism (the score responds to one degree of freedom, and it is not attenuation)
rather than only on a tuned-share sweep.

---

## 0. Known-answer rows (rule 10's corollary) — five, all pass

| # | what | expected | measured | source |
|---|---|---|---|---|
| 1 | 194's achieved-repression table, by an **independent** closed form on the shipped tables | `0.15515195` / `0.29975643` / `0.51016574` | `0.15515196` / `0.29975641` / `0.51016574` | max abs Δ **1.6e-07 / 1.2e-07 / 1.9e-07**, corr 1.000000 |
| 2 | PLOLO a15 `pds`, `k562_gwps` fold, equal weights | 0.79578409 (194 §0.3) | **0.795784** | `work/score204.log` |
| 2b | PLOLO a23 `pds`, same fold — and `q`, `q_w` | 0.78912131 / 0.50164311 / 0.50487994 | **0.789121 / 0.501643 / 0.504880** | same |
| 3 | the a15 → a23 mass-fix price, 4 folds, pre-committed unweighted mean | **+0.004014** scaled (162 §2) | **+0.004014** (+0.01%) | `analyse204.py` |
| 4 | the weighted mixer at `w = (1,1,1)` reproduces 162's own cached `pre` | 0 | max abs dev **4.1e-08** (float32 cache round-trip) | `audit.max_abs_dev_from_162_pre` |
| 5 | the Gram/triple attenuation estimator recovers a **planted** attenuation | 0.500 / 0.800 | **0.5063 (+1.26%) / 0.7993 (−0.09%)** | `atten2_204.py` |

Row 5 is the one that makes Phase 1 readable at all: the estimator that returns "Orion is 1.4×
K562" has been shown a right answer and returned it to ~1%. Its disagreement with the residual
route is therefore not an implementation bug.

### Rule 11 — a path and a sha for every operator loaded

| path | sha256 (8) | shape | used by |
|---|---|---|---|
| `data/donors/whitener_union_pruned.npz` | **`e2f65ecf`** | `evals[17942]`, `evecs[17942,17942]`, basis `[29394, 17942]` (from its own `.json`) | `kolf2` fold (production) |
| `experiments/062-cd4-post-whiten/cloud-output/whiteners/held_k562_gwps.npz` | `13d988fc` | `[17806,17806]` | `k562_gwps` fold |
| `.../held_orion_hct116.npz` | `7cecbf7f` | `[17942,17942]` | `orion_hct116` fold |
| `.../held_orion_hek293t.npz` | `97e486e0` | `[17942,17942]` | `orion_hek293t` fold |
| *(present on disk, NOT loaded)* `whitener_union.npz` `5608e807`, `whitener_union_full.npz` `965d3cef`, `whitener_full.npz` `dc6bd2ce` | | | listed so the near-name confusion of rule 11 is closed by hash, not by filename |

**The mandate's "29,394 × 17,942" is the BASIS the operator was fitted on** (`basis` in
`whitener_union_pruned.json`), not the stored array; the stored operator is `17,942²`. Both are
recorded above so the next agent does not have to re-resolve it.

Frozen inputs (sha256, first 16): `k562_gwps/delta.npy` `0c9422280fed7afb`, `b_ctrl.npy`
`d06cf51d51cde1ae`; `orion_hct116/delta.npy` `72c7a01dc2869bd3`, `b_ctrl.npy` `b3044d439cfac932`;
`orion_hek293t/delta.npy` `d95202d984d9322d`, `b_ctrl.npy` `e04234cffbe10273`. The first two match
194's recorded source hashes.

### The fold record (162's, unchanged)

| fold (held) | sources | operator | panel targets | scored genes | vote histogram |
|---|---|---|---:|---:|---|
| `k562_gwps` | hct116, hek293t | `held_k562_gwps` `13d988fc` | 269/300 | 7,448 | 2:259, 1:5, 0:5 |
| `orion_hct116` | k562, hek293t | `held_orion_hct116` `7cecbf7f` | 290/300 | 17,806 | 2:259, 1:31 |
| `orion_hek293t` | k562, hct116 | `held_orion_hek293t` `97e486e0` | 294/300 | 17,806 | 2:259, 1:34, 0:1 |
| **`kolf2`** (leak-free) | **all three** | **`whitener_union_pruned` `e2f65ecf`** | 273/300 | 17,806 | 3:240, 2:28, 1:5 |

`debench.Surface.pds()` was never called (188 Phase 3). `h1_2025` was never used as a truth or a
fold; it appears only as an *instrument* in Phase 1, which touches no score.

---

## 1. The repression table, reproduced independently and extended to the panel

**How residual expression is computed** (194's definition, stated exactly so it can be checked):
for target `T`, over `T`'s own perturbed cells,

```
r_T = sum_c obs_c(T) / sum_c exp_c(T),      exp_c(T) = p_ctrl(T) * depth_c
```

`obs_c(T)` is that cell's count of **its own target gene only**; `depth_c` its library size;
`p_ctrl(T)` the target gene's fraction of the pooled control profile. Cells whose target gene is
off the donor's gene axis are excluded (`r_T = NaN`). No knockdown gate, no `KD_MIN_EXPECTED`
threshold — every perturbed cell of the target contributes.

Because numerator and denominator are both **pooled sums**, that ratio has a closed form on the
shipped pseudobulk table and needs no cell streaming at all:

```
delta[T,g] = log1p(TS*f_T(g)) - log1p(TS*f_ctrl(g)),   b_ctrl[g] = log1p(TS*f_ctrl(g))
=>  r_T = expm1(delta[T,T] + b_ctrl[T]) / expm1(b_ctrl[T])
```

That is the independent route. It reproduces 194's 10-million-cell streaming pass to **2e-7**
(known-answer row 1) — two entirely different code paths on the same quantity.

| | K562-gwps | Orion HCT116 | Orion HEK293T |
|---|---:|---:|---:|
| **194's restricted row set** (its own number, reproduced) | | | |
| rows / with a readout | 1,569 / 1,230 | 790 / 785 | 794 / 793 |
| median residual | **0.15515** | **0.29976** | **0.51017** |
| clearing >80%, over all rows | **47.6%** | 21.4% | 4.4% |
| clearing >80%, over rows with a readout | 60.7% | 21.5% | 4.4% |
| **the 300-target PANEL** (new here) | | | |
| panel targets with a readout | **216 / 259** | 259 / 259 | 259 / 259 |
| **median residual** | **0.1348** | **0.2983** | **0.5039** |
| implied `a = 1 − r` | **0.865** | **0.702** | **0.496** |
| clearing >80% (of those with a readout / of all 259) | **74.5% / 62.2%** | 17.4% | 1.5% |
| **the full shipped table** | | | |
| rows / with a readout | 9,787 / 7,617 | 17,767 / 17,265 | 17,977 / 17,710 |
| median residual (p10 / p90) | 0.1651 (0.048 / 0.577) | 0.2787 (0.000 / 0.813) | 0.4919 (0.000 / 0.907) |

**Two corrections for the record, both small and both definitional.**

* 194 §1's "**47.6%** of K562 targets achieve >80% knockdown" counts the 339 rows with **no
  readout as failures**. Among rows that can be read it is **60.7%**. Both are reproduced exactly;
  the other two donors' figures are insensitive (785/790 and 793/794 rows are readable).
* On the **panel** the spread is *larger* than the mandate's premise, not smaller: median residual
  0.1348 → 0.5039 is **3.74×**, and only **1.5%** of HEK293T's panel targets clear >80% against
  **74.5%** for K562. §2.1's selection criterion ("sufficiently high basal expression of the target
  gene") is why panel targets read better than the extras.

---

## 2. The algebraic fact that decides what an arm can even be

`prepare_main_source` (`experiments/062-cd4-post-whiten/run.py:239-297`) computes

```
effects_i = (delta_i - mean_over_all_rows(delta_i)) * (EFFECT_SCALE / raw_median_norm_i)
```

and `raw_median_norm_i` is the median centred row norm of that same matrix. Every term is degree-1
homogeneous in `delta_i`, so **scaling a donor's entire table by any constant leaves `effects`
unchanged.** Measured, not argued (`homog204.py`): `h1_2025`'s table scaled by 2.5 gives
`raw_median_norm` ratio **2.5000000032** and `effects` identical to **4.1e-08 relative**.

Three consequences, and they reshape the mandate's model:

1. **Production already implements "de-attenuate then average"** — with the donor's own median row
   norm as the attenuation estimate. The mandate's de-attenuation applied at the table level is an
   exact no-op.
2. So the live question is not *whether* to de-attenuate but *whether the median row norm is the
   right attenuation estimate*. It is not obviously wrong: in the production frame on the official
   axis the three donors' median row norms are **6.721 / 7.027 / 6.921** — within **5%** of each
   other — while their measured repression differs **1.74×**. The normalisation flattens amplitude
   completely and is blind to repression.
3. The only expressible arm is therefore a **relative weight** inside `build_baseline`. Every arm
   below is exactly that, and `w = (1,1,1)` reproduces the incumbent expression identically.

**One property deliberately preserved (039).** `build_baseline` divides by the summed weight of the
donors that carry the *target*, while a donor contributes zero on a gene it does not *measure* — so
a gene seen by one donor is emitted at a fraction of full magnitude, an implicit coverage
shrinkage 039 measured at **+0.019 macro PDS**. Re-weighting necessarily moves that fraction (K562
carries 7,448 of the official axis against Orion's 17,678/17,756). That is a real consequence of
re-weighting, not a bug, and it is one reason the `pds` surface responds to the K562 share at all
(§4.3).

---

## 3. Phase 1 — `a_i` from independent routes, and why they disagree

All routes are computed in the **production frame** (the object `build_baseline` averages) on a
fixed shared set: **2,259 targets = 259 panel + a seeded draw of 2,000 non-panel targets carried by
all three donors**. Gene coverage is controlled by restricting every inner product to the **7,310
columns all three donors measure** — otherwise a ratio compares a 7,448-gene donor with a
17,700-gene one, which is 039's lever smuggled into this one.

### 3.1 The routes

| route | what it uses | `a_rel` (K562 = 1) | HCT116 | HEK293T |
|---|---|---|---:|---:|
| **A — self-gene residual**, panel scope | the target's own transcript | | **0.811** | **0.573** |
| A — self-gene residual, shared scope | | | 0.819 | 0.571 |
| **B — 3-donor Gram (triple)**, panel, gene-controlled | `a_i = √(G_ij G_ik / G_jk)` | | **1.414** | **1.288** |
| B — triple, shared set, gene-controlled | | | 1.305 | 1.173 |
| B — triple, **not** gene-controlled | | | 1.694 | 1.524 |
| **B′ — naive slope on the others' consensus** | `⟨d_i, m_−i⟩/‖m_−i‖²` | | 1.401 | 1.335 |
| **B — fourth-line instrument**, `kolf2` (n=240, all panel) | `G_iX / G_jX` | | 0.639 `[0.16, 1.30]` | 0.436 `[−0.03, 1.03]` |
| B — instrument, `h1_2025` (n=69, 20 panel) | | | 0.644 `[0.49, 0.85]` | 0.513 `[0.44, 0.60]` |
| B — instrument, `nadig_jurkat` (n=469, 0 panel) | | | 0.485 `[0.42, 0.55]` | 0.466 `[0.41, 0.52]` |
| B — instrument, `rpe1` (n=458, 0 panel) | | | 1.221 `[1.04, 1.45]` | 0.613 `[0.54, 0.69]` |
| B — instrument, `k562_essential` **(same cell line — biased, shown as the control)** | | | 0.363 `[0.32, 0.41]` | 0.313 `[0.28, 0.35]` |

Intervals are 2,000-draw target bootstraps. **Route A and route B disagree, and they disagree in
opposite directions.** Route B′ agrees with B in direction and is attenuated relative to it, which
is the errors-in-variables bias its regressor predicts — so B′ is consistent with B and adds no
independent information.

### 3.2 The blob control — it is already deflated, and deflating it further changes nothing

The target-independent mean vector, as a share of table energy, in the production frame on the
shared set: **K562 0.044%, HCT116 0.071%, HEK293T 0.185%**. These are far below 194's
0.81%/0.52%/1.00% because `prepare_main_source` already subtracts each bank's own gene mean over
**all** its rows. Cross-donor blob alignment: `cos(HCT116, HEK293T) = +0.218` against
`cos(K562, HCT116) = −0.107` and `cos(K562, HEK293T) = −0.101` — the two Orion lines' blobs point
the same way and K562's points away, which is the technical axis 194 named.

Subtracting each donor's blob before estimating moves `a_rel` by **under 1%** on every route
(triple, shared set, all official columns: **1.6941/1.5236 → 1.6818/1.5143**, i.e. −0.7%/−0.6%;
naive slope 1.4008/1.3353 → 1.3942/1.3288; instruments unchanged to 3 decimals). **Reported with
and without, as the mandate requires; the blob is not what drives the disagreement.**

### 3.3 What DOES drive it — a one-factor test with a positive control

The mandate's model `d_i = a_i s + ε_i` with `ε` independent across lines predicts the off-diagonal
cross-line correlation matrix is **rank one** (`ρ_ij = c_i c_j`). With four or more lines that is
over-identified and can be rejected. `k562_essential` is the **same cell line** as `k562_gwps`, so
it must show a large positive residual — the test's built-in positive control.

| block | lines | n targets (panel) | rank-1 RMS resid / mean\|ρ\| | verdict | `a_rel` HCT116 | `a_rel` HEK293T |
|---|---|---:|---:|---|---:|---:|
| `panel_block` | k562, hct116, hek293t, **kolf2** | 240 (240) | **11.1%** `[0.0006, 0.0037]` | not rejected (2 df, weak instrument) | 1.438 `[0.82, 2.19]` | 1.254 `[0.74, 1.88]` |
| `h1_block` | k562, hct116, hek293t, **h1_2025** | 69 (20) | 7.8% | not rejected (n = 69) | 0.656 `[0.46, 0.98]` | 0.581 `[0.48, 0.76]` |
| `five_no_sameline` | + rpe1, nadig_jurkat | 445 (0) | **19.7%** | **REJECTED**, \|z\| to 11.0 | 0.880 `[0.76, 1.03]` | 0.679 `[0.60, 0.78]` |
| `replogle_block` | + k562_essential | 380 (0) | **22.4%** | **REJECTED**, \|z\| to 11.2 | 0.623 `[0.54, 0.72]` | 0.504 `[0.44, 0.58]` |

**The positive control fires**: `k562_gwps ↔ k562_essential` residual **+0.0495, z = +11.2**, the
largest in its block. The test is sensitive.

**And the second-largest structural excess is the Orion pair.** `orion_hct116 ↔ orion_hek293t`
residual **+0.0306, z = +5.14** in the five-line block and **+0.0578** in the six-line block —
i.e. the two X-Atlas lines agree with each other **beyond what any single shared component
allows**, by more than a same-cell-line pair does in relative terms. This is 194 §3's
target-independent-blob finding generalised: the excess is a **subspace**, not one vector, which is
why blob deflation (§3.2) does not remove it.

**That is the whole mechanism of the route disagreement.** The 3-donor triple has **zero degrees of
freedom** — with three lines the Gram is exactly identified and the estimator is *forced* to use
`G(hct116, hek293t)`. Inflating that one entry pushes both Orion amplitudes up and K562's down,
which is precisely the reported `(1, 1.41, 1.29)`. Every route that avoids that product — a fourth
line as instrument, or a ≥5-line fit that can average the excess away — returns `(1, ~0.6, ~0.5)`,
in agreement with the self-gene residual's `(1, 0.81, 0.57)`.

**But the panel is where the answer is needed, and there it cannot be settled.** The only
independent line with real panel coverage is `kolf2`, and it is too noisy: its interval on
`a_rel(HCT116)` is **[0.16, 1.30]**. The four independent estimates of `a_rel(HEK293T)` span
**0.504 → 1.254**, a 2.5× range with barely-overlapping intervals, and the estimate depends on which
fourth line and which target population is used.

> **The Phase-1 result, stated as the mandate asked it to be: the two routes disagree, the
> disagreement is traced to a measured excess in the Orion↔Orion inner product (z = +5.1 to +7.6),
> and `a_i` is NOT identifiable as a scalar donor property on the 300-target panel with the lines
> that exist. It is not averaged away here, and it should not be averaged away later.**

A further caution for whoever revisits this: the two populations do not merely differ in precision.
The essential-gene population (`replogle_block`, 0 panel targets) and the panel population
(`panel_block`) give `a_rel(HEK293T)` of 0.504 `[0.44, 0.58]` and 1.254 `[0.74, 1.88]`. §2.1 says
the panel was deliberately built to be diverse and *non*-common-essential, so an attenuation
calibrated on essential-gene screens is being read on the wrong population.

---

## 4. Phase 2 — the weightings, scored

162's instrument, byte-unchanged: four PLOLO folds, challenge frames A/B/C, shipped geometry
(`alpha 0.982, gamma 0.65, lam 0, cd4_f 0.10`), a23 render for `pds`/`cos`, a15 render for `q`/`q_w`
(194 §0.3: the mass repair puts a 1e-3 numerical floor under the sign channel). **The only thing
that varies across arms is the donor weight vector.** The truth is the held-out donor's own shipped
table on every arm; the whitening operator is the fold's frozen one on every arm and is **never
refitted** per weighting.

Absolute equal-weight baseline `pds`: **0.78912 / 0.79182 / 0.77738 / 0.61013** (269 / 290 / 294 /
273 panel targets). Intervals are a 4,000-draw paired target-cluster bootstrap resampling targets
within each fold (194 §5.2's protocol).

### 4.1 Every arm, sorted by the priced total

`ovr` columns are board-overall. `mse` is priced through 197's identity from the cosine (§6);
`TOT|abs` is the same total under the conservative absolute-cosine transfer instead of the relative
one. Weights shown normalised to sum 3.

| arm | w (k562, hct116, hek293t) | Δ`pds` | 95% CI | p>0 | Δcos | Δcos 95% CI | ovr `pds` | ovr `mse` | **TOTAL** | TOT\|abs |
|---|---|---:|---|---:|---:|---|---:|---:|---:|---:|
| `deatten_slope` | (1.146, 0.878, 0.976) | +0.00076 | [−0.00046, +0.00202] | 0.889 | −0.000086 | [−0.000155, −0.000015] | +0.00028 | −0.00011 | **+0.00017** | +0.00025 |
| `k13` (sweep) | (1.182, 0.909, 0.909) | +0.00082 | [−0.00058, +0.00220] | 0.876 | −0.000116 | [−0.000194, −0.000034] | +0.00031 | −0.00015 | +0.00015 | +0.00026 |
| `k115` (sweep) | (1.095, 0.952, 0.952) | +0.00049 | [−0.00025, +0.00124] | 0.904 | −0.000045 | [−0.000087, −0.000001] | +0.00018 | −0.00006 | +0.00012 | +0.00017 |
| `k15` (sweep) | (1.286, 0.857, 0.857) | +0.00104 | [−0.00117, +0.00321] | 0.826 | −0.000228 | [−0.000348, −0.000102] | +0.00039 | −0.00030 | +0.00008 | +0.00030 |
| **`equal` (incumbent)** | (1.000, 1.000, 1.000) | 0 | — | — | **0** | — | 0 | 0 | **0** | 0 |
| `h085` (sweep) | (1.053, 1.053, 0.895) | +0.00003 | [−0.00080, +0.00087] | 0.523 | −0.000058 | [−0.000107, −0.000007] | +0.00001 | −0.00008 | −0.00007 | −0.00001 |
| `h115` (sweep) | (0.952, 0.952, 1.095) | −0.00025 | [−0.00102, +0.00052] | 0.261 | −0.000008 | [−0.000053, +0.000036] | −0.00009 | −0.00001 | −0.00010 | −0.00010 |
| **`perm_a_5` (falsifier)** | (1.246, 0.722, 1.032) | **+0.00074** | [−0.00175, +0.00329] | **0.734** | −0.000294 | [−0.000430, −0.000158] | +0.00028 | −0.00039 | −0.00011 | +0.00017 |
| `h07` (sweep) | (1.111, 1.111, 0.778) | −0.00006 | [−0.00177, +0.00168] | 0.462 | −0.000210 | [−0.000317, −0.000100] | −0.00002 | −0.00028 | −0.00030 | −0.00010 |
| **`a_self`** (`w ∝ a_i`, route A) | (1.246, 1.032, 0.722) | **+0.00019** | [−0.00214, +0.00261] | **0.561** | −0.000327 | [−0.000470, −0.000182] | +0.00007 | −0.00043 | **−0.00036** | −0.00005 |
| `k085` (sweep) | (0.895, 1.053, 1.053) | −0.00097 | [−0.00180, −0.00012] | 0.013 | −0.000002 | [−0.000053, +0.000046] | −0.00036 | −0.00000 | −0.00036 | −0.00036 |
| `h13` (sweep) | (0.909, 0.909, 1.182) | −0.00079 | [−0.00222, +0.00061] | 0.140 | −0.000063 | [−0.000147, +0.000020] | −0.00030 | −0.00008 | −0.00038 | −0.00032 |
| `a_inst` (`w ∝ a_i`, h1_2025 instrument) | (1.391, 0.895, 0.714) | +0.00039 | [−0.00267, +0.00346] | 0.601 | −0.000431 | [−0.000604, −0.000252] | +0.00014 | −0.00057 | −0.00043 | −0.00001 |
| `minvar_self` (`w ∝ a_i/σ_i²`, a from route A) | (0.944, 1.171, 0.885) | −0.00089 | [−0.00215, +0.00033] | 0.081 | −0.000089 | [−0.000161, −0.000017] | −0.00033 | −0.00012 | −0.00045 | −0.00037 |
| `a_slope` (`w ∝ a_i`, route B triple) | (0.863, 1.125, 1.012) | −0.00117 | [−0.00232, +0.00001] | 0.026 | −0.000031 | [−0.000100, +0.000037] | −0.00043 | −0.00004 | −0.00048 | −0.00045 |
| `perm_a_1` (falsifier) | (1.032, 0.722, 1.246) | −0.00037 | [−0.00276, +0.00216] | 0.390 | −0.000303 | [−0.000442, −0.000162] | −0.00014 | −0.00040 | −0.00054 | −0.00025 |
| `share_2_1_1` (010/019) | (1.500, 0.750, 0.750) | +0.00034 | [−0.00341, +0.00400] | 0.582 | −0.000534 | [−0.000732, −0.000325] | +0.00013 | −0.00071 | −0.00058 | −0.00007 |
| `perm_a_3` (falsifier) | (1.032, 1.246, 0.722) | −0.00100 | [−0.00331, +0.00123] | 0.186 | −0.000339 | [−0.000479, −0.000202] | −0.00037 | −0.00045 | −0.00082 | −0.00050 |
| `h15` (sweep) | (0.857, 0.857, 1.286) | −0.00164 | [−0.00380, +0.00052] | 0.072 | −0.000178 | [−0.000305, −0.000049] | −0.00061 | −0.00024 | −0.00085 | −0.00068 |
| `k07` (sweep) | (0.778, 1.111, 1.111) | −0.00245 | [−0.00418, −0.00069] | 0.003 | −0.000080 | [−0.000190, +0.000024] | −0.00091 | −0.00011 | −0.00102 | −0.00094 |
| `h05` (sweep) | (1.200, 1.200, 0.600) | −0.00110 | [−0.00423, +0.00205] | 0.238 | −0.000619 | [−0.000820, −0.000420] | −0.00041 | −0.00082 | −0.00123 | −0.00063 |
| **`deatten`** (de-attenuate then average, route A) | (0.763, 0.921, 1.316) | **−0.00291** | **[−0.00541, −0.00048]** | 0.009 | −0.000216 | [−0.000367, −0.000069] | −0.00109 | −0.00029 | **−0.00137** | −0.00117 |
| `perm_a_4` (falsifier) | (0.722, 1.032, 1.246) | −0.00329 | [−0.00569, −0.00094] | 0.002 | −0.000176 | [−0.000327, −0.000037] | −0.00123 | −0.00023 | −0.00146 | −0.00129 |
| `perm_a_2` (falsifier) | (0.722, 1.246, 1.032) | −0.00333 | [−0.00558, −0.00104] | 0.003 | −0.000185 | [−0.000324, −0.000049] | −0.00124 | −0.00025 | −0.00149 | −0.00131 |
| **`minvar`** (`w ∝ a_i/σ_i²`, route B) | (0.618, 1.208, 1.174) | **−0.00532** | **[−0.00827, −0.00228]** | 0.001 | −0.000328 | [−0.000520, −0.000146] | −0.00198 | −0.00043 | **−0.00242** | −0.00211 |
| `share_1_2_2` = `k05` | (0.600, 1.200, 1.200) | −0.00572 | [−0.00883, −0.00254] | 0.000 | −0.000367 | [−0.000569, −0.000175] | −0.00213 | −0.00049 | −0.00262 | −0.00227 |
| `share_4_1_1` (010's optimum) | (2.000, 0.500, 0.500) | −0.00352 | [−0.00991, +0.00298] | 0.144 | −0.001429 | [−0.001778, −0.001060] | −0.00131 | −0.00189 | −0.00320 | −0.00184 |
| **`a_self_pt`** (per-target `w ∝ 1−r_it`) | per target | **−0.02162** | [−0.03057, −0.01274] | 0.000 | −0.002066 | [−0.002798, −0.001387] | −0.00806 | −0.00273 | **−0.01079** | −0.00881 |
| **`deatten_pt`** (per-target `w ∝ 1/(1−r_it)`) | per target | **−0.02290** | [−0.03181, −0.01420] | 0.000 | −0.001999 | [−0.002732, −0.001301] | −0.00854 | −0.00264 | **−0.01118** | −0.00926 |
| *194's `kdrescale`* (per-target de-attenuation at the TABLE level, cited) | — | *+0.00068 ± 0.00288* | — | *0.808 @λ0.25* | *−0.00009* | — | — | — | *null* | — |

### 4.2 Equal weighting is the maximum of the transfer cosine

**Every one of the 29 alternatives has a lower macro cosine than `equal`, and 24 of them are below
it at 95%.** The best arm on `pds` (`deatten_slope`) still loses **−0.39%** relative cosine; the
attenuation arm `a_self` loses **−1.50%**; the minimum-variance arm loses **−1.50%**;
`share_4_1_1` loses **−6.53%**.

Fitting the surface (quadratic in the two weight coordinates, 28 scalar arms) gives the argmax
directly:

| endpoint | R² | argmax K562 / geomean(Orion) | argmax HCT116 : HEK293T | peak height over `equal` |
|---|---:|---:|---:|---:|
| Δ`pds` | **0.9969** | **1.559** | **0.924** | +0.00104 raw |
| Δcos | **0.9813** | **0.947** | **0.974** | **+0.000003** |

**On the cosine — the channel §8.2 says is the entire remaining gap — the fitted optimum is
`(0.95, 1.02, 1.02)`, i.e. the incumbent, and the peak height above the incumbent is 3e-6.**

### 4.3 The one degree of freedom the score actually responds to

| endpoint | partial R², **K562 share** | partial R², **HCT116:HEK293T split** |
|---|---:|---:|
| Δ`pds` | **0.978** | **0.046** |
| Δcos | **0.930** | **0.133** |

The K562 share carries essentially everything. And the split is the axis on which the attenuation
hypothesis makes its cleanest prediction: **the two Orion lines are the same platform, the same
chemistry, the same depth and the same builder, differing essentially only in achieved repression
(0.298 vs 0.504, a 1.69× ratio in `1−r`: 0.702 vs 0.496 → predicted weight ratio 1.43)**. The
measured argmax of that ratio is **0.924**, and the dedicated sweep along it reads

| HEK293T share | ×0.5 | ×0.7 | ×0.85 | **×1.0** | ×1.15 | ×1.3 | ×1.5 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Δ`pds` | −0.00110 | −0.00006 | +0.00003 | **0** | −0.00025 | −0.00079 | −0.00164 |
| Δcos | −0.000619 | −0.000210 | −0.000058 | **0** | −0.000008 | −0.000063 | −0.000178 |

**Down-weighting the 49%-repression donor to exactly the level its repression prescribes (×0.70)
is worth −0.00006 raw `pds` and −0.000210 cosine.** With platform, depth and chemistry held, the
attenuation hypothesis has nothing left to explain and it explains nothing.

### 4.4 The K562-share ridge, and where it comes from

| K562 share | ×0.5 | ×0.7 | ×0.85 | **×1.0** | ×1.15 | ×1.3 | ×1.5 | ×2 | ×4 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Δ`pds` (4-fold mean) | −0.00572 | −0.00245 | −0.00097 | **0** | +0.00049 | +0.00082 | +0.00104 | +0.00034 | −0.00352 |
| Δ`pds`, **`kolf2` fold only** (leak-free, 3 sources, shipped operator) | −0.00265 | −0.00110 | −0.00064 | **0** | −0.00030 | −0.00057 | −0.00034 | −0.00089 | −0.00304 |

**The ridge exists only on the two-source folds.** On `kolf2` — the only fold that uses the exact
shipped operator *and* the exact shipped three-source mixture, and the only leak-free one — **every
K562 up-weighting is negative**. The mean is the pre-committed predictor and is what §4.1 reports,
but the decomposition says the small `pds` ridge is a property of a 2-donor consensus, not of
production's 3-donor one. That is consistent with 019's own curve, which was measured on 2-donor
LODO and found the optimum moving *toward* equal as the geometry improved.

---

## 5. Phase 3 — the falsifiers

### 5.1 Permuted attenuation — **the control is not beaten** (mandate falsifier 1)

All five non-identity permutations of `a = (0.845, 0.700, 0.490)` across the three donors:

| arm | w | Δ`pds` | p(Δ>0) |
|---|---|---:|---:|
| **`a_self`** (the real vector) | (1.246, 1.032, 0.722) | **+0.00019** | **0.561** |
| **`perm_a_5`** (the two Orion lines swapped) | (1.246, 0.722, 1.032) | **+0.00074** | **0.734** |
| `perm_a_1` | (1.032, 0.722, 1.246) | −0.00037 | 0.390 |
| `perm_a_3` | (1.032, 1.246, 0.722) | −0.00100 | 0.186 |
| `perm_a_4` | (0.722, 1.032, 1.246) | −0.00329 | 0.002 |
| `perm_a_2` | (0.722, 1.246, 1.032) | −0.00333 | 0.003 |

`perm_a_5` holds K562's weight **identical to the real arm** and swaps only the two Orion lines —
and it scores **4× better** with a higher `p(Δ>0)`. This is 194's trap exactly (`kd_perm@0.25`
p = 0.947 against the real arm's 0.907), and it has the same reading: the arm's apparent gain lives
entirely in the K562 share, which permutations preserve, and the attenuation ordering *within* the
Orion pair is not merely uninformative — reversing it does better.

### 5.2 The truth was held fixed (mandate falsifier 2)

194 measured that rebuilding the truth with the same estimator swings `kd_post` from −0.00156 to
+0.00433, a 0.006 swing from the design choice alone. Here the truth is the held-out donor's own
shipped table on every arm — a weighting cannot apply to a single donor, so the trap is
structurally absent. Two further quantities were held fixed for the same reason, and both would
have been second variables: the **whitening operator** (fitted on the equal-weight donor basis, and
never refitted per arm), and the **PerturbMatch blend** (`PM_BLEND = 0.5` on K562's shared columns).

### 5.3 Both surfaces and the leak-free fold (mandate falsifier 3)

All four folds are reported per arm in `work/six204.json`; §4.4 breaks out `kolf2`. **No absolute is
quoted from the Orion pair** — the Orion folds appear only as differences against a matched `equal`
arm on the same row set, and the one place their behaviour is decisive (the K562 ridge) is
explicitly attributed to them and contrasted with `kolf2`.

### 5.4 The model's own predicted effect size — it is rejected, not merely unsupported

A negative is only worth something if the instrument could have seen the positive. The
minimum-variance model's own measured inputs give the answer (`power204.py`):

| gene set / target scope | optimal weights remove … of consensus noise variance | predicted relative cosine gain | predicted overall via `mse` alone |
|---|---:|---:|---:|
| common genes / shared | 15.4% | **+8.5%** | +0.00247 |
| common genes / panel | 19.4% | **+11.2%** | +0.00323 |
| all official / shared | 17.4% | +9.8% | +0.00285 |
| all official / panel | 27.1% | **+16.8%** | +0.00486 |

(`Var_equal = Σv_i/k²`, `Var_opt = 1/Σ(1/v_i)`, `v_i = N_i/A_i²` from the Gram decomposition;
`cos` gain `= √(1/(1−f(1−c_pred²))) − 1` with `c_pred = 0.130` implied by PLOLO's own macro cosine
0.021885 and the panel block's mean single-donor reliability 0.169. `c_pred² ≪ 1`, so the gain is
insensitive to it: `c_pred ∈ [0.12, 0.30]` moves +8.5% to +8.9%.)

**Measured, on the model's own prescribed weights (`minvar`): −1.50% relative, 95% CI
[−2.38%, −0.67%].** Wrong sign, and the interval excludes the smallest prediction by more than 7×.
The same arm is **−0.00532 raw `pds`, p = 0.001** — both channels reject it independently.

The reason is 022's, unchanged and now with a `z`: `N_i` is not independent noise. §3.3 measures the
Orion pair sharing structure beyond the single factor at **z = +5.1 to +7.6**, so down-weighting a
donor because it carries more `N` does not remove a proportional amount of consensus noise. This is
the fourth time this repo has measured that "improve the estimate" instincts fail here on the same
mechanism (AGENTS.md §4.2 point 4; exp 022; exp 194 §6), and the first time the shared component of
the residual has been given a significance.

### 5.5 The per-target route fails 30× worse, and fails where 194 predicted

`a_self_pt` and `deatten_pt` apply the SAME correction at per-target resolution, using each donor's
own measured `r_it` (floor 0.25, 194's floor). They lose **−0.0216 / −0.0229** raw `pds`, and the
loss is concentrated exactly where 194's mechanism says it must be:

| fold (held) | sources | `a_self_pt` Δ`pds` | `deatten_pt` Δ`pds` |
|---|---|---:|---:|
| `k562_gwps` | both Orion lines | −0.00209 | +0.00072 |
| `orion_hct116` | K562 + HEK293T | **−0.04622** | **−0.05181** |
| `orion_hek293t` | K562 + HCT116 | **−0.03569** | **−0.03911** |
| `kolf2` | all three | −0.00250 | −0.00141 |

The two folds that mix **across platforms** are destroyed; the fold whose two sources are both
Orion is barely touched. Per-target de-attenuation is per-target donor *selection*, and the
selector is a readout 194 measured to be inverted on Orion (`corr(proj, residual)` −0.232 /
−0.159). It is also the same failure 010 recorded for basal-similarity donor weighting: a
per-target signal "concentrates weight confidently on the wrong donor". Note that 194's own
table-level `kdrescale` was a null (+0.00068 ± 0.00288), not a loss — the difference is that at the
table level the row rescale is largely absorbed by the median-norm normalisation of §2, whereas in
the mixing frame it acts at full strength.

---

## 6. Rule 12 — all six members

The arm class is **profile-changing, shape-only**: `alpha`, the emission law, the spread policy and
`keep_top` are untouched; only the pre-whitening effect matrix moves. Rule 13's **−0.0146 `reach`
toll does not apply** — it is an EMISSION-LAW toll (audit 181; 169 measured +0.0006 for profile
arms; the `rc_hybrid` board point confirmed the class boundary). Gate: **+0.0019**.

| member | instrument | value (best arm `deatten_slope` / mandate arm `a_self`) | verdict |
|---|---|---|---|
| `pds` | PLOLO, 4 folds, challenge frames, a23, paired target bootstrap | **+0.00076 [−0.00046, +0.00202]** / **+0.00019 [−0.00214, +0.00261]** raw → **+0.00028 / +0.00007** overall | inside the noise floor (162's LOO RMS 0.001494 scaled ≈ 0.00067 raw) |
| `mse` | 096/197's identity, `claim` route. Only `X` moves: `d raw_mse/d cos = −2√(q(1+ε)) = −2.0826` → **+0.3653 overall per unit cosine**, i.e. §8.2's **+0.00365 per +0.01** re-derived | Δcos **−0.000086 / −0.000327** raw (−0.39% / −1.50% relative) → **−0.00011 / −0.00043** overall (relative transfer); **−0.00003 / −0.00005** (absolute transfer) | **negative on every arm**; equal is the maximum |
| `nmae` | 187 §5's exact identity. With `alpha` fixed, `K_lo`/`K_hi`/`W_hi` are fixed and only `q_w` lives: `d raw_nmae = −2·W_hi·d q_w`, `W_hi = 0.0153` | Δ`q_w` **−0.000038 / −0.000076** → Δraw_nmae **+1.2e−6 / +2.3e−6** → **−3.3e−7 / −6.5e−7** overall | **pinned.** Largest \|Δ`q_w`\| over all 30 arms is **1.04e−3** (`a_self_pt`), i.e. \|Δ overall\| ≤ **8.9e−6**; over the 28 scalar-weight arms it is 2.6e−4, i.e. ≤ 2.2e−6 |
| `reach` | ordering half only: Δ`q` (a15 render). Not locally priceable (rule 13: 2/6 signs, worse than chance; rule 9: `n_real`) | Δ`q` **+0.00067 / +0.00102** | unpriced; see the note below |
| `fid` | ordering half: same Δ`q`. Yield half unreadable (rule 9) | Δ`q` as above | unpriced |
| `jac` | call set = `keep_top` on the delivered profile; `keep_top` is not touched and 192 owns the call set | — | unmeasured, mechanism stated |

**Two disclosures.**

1. **The h1-test `nmae` route (188) was not run**, and this is a deliberate limitation, not an
   oversight. It would require rebuilding the h1-test artifact from re-weighted donor tables for
   each arm × 5 seeds — a build of a scale this mandate's 10 GB cap and the concurrent production
   build do not allow. Two things make it a poor instrument for *this* arm class anyway: **h1-test
   carries only `k562_gwps` + `orion_hek293t`**, so a three-donor weighting is not even expressible
   on it; and 188's route prices arms that change the *emission law*, which this class does not.
   The 187 identity above bounds the member at **|Δ| ≤ 2.2e−6 overall across all 30 arms**, which
   is four orders of magnitude below the gate.
2. **The `mse` identity holds `A` (delivered energy) fixed.** Checked, not assumed: the delivered
   median row norm relative to `equal` stays inside **±1.5%** for every arm near the incumbent
   (`deatten_slope` 0.987–1.008, `k13` 0.991–1.010, `a_self` 1.002–1.024). It reaches **1.131** for
   `share_4_1_1` and **1.140** for `deatten_pt`; for those arms the amplitude term is real and my
   pricing **understates** their loss, which is conservative in the direction that matters since
   they already lose.

**An observation for the ordering channel, explicitly NOT priced.** Δ`q` — unweighted per-gene sign
accuracy, the ordering half of `reach`/`fid` — is *positive* for every K562-heavy arm and grows
monotonically with the K562 share: `k15` +0.00114, `a_inst` +0.00161, `share_4_1_1` +0.00352, while
Orion-heavy arms go the other way (`minvar` −0.00138). §8.1 rule 1 forbids pricing an ordering
signal here, and rule 13 says `reach` is not locally priceable at all (2/6 signs). It is recorded
because 172 measured **+0.052 raw `reach` from ordering alone** on our own predictions, and the sign
ordering is 185/189/190's territory. **A K562-heavy mix trades cosine for sign ordering** —
`share_4_1_1` is −6.53% cosine for +0.0035 `q` — which is the same trade AGENTS.md §4.2 records for
whitening, in miniature. It is not a recommendation: on every priced member that arm loses.

---

## 7. Verdict in overall units

| arm | what it is | overall (relative cosine transfer) | overall (absolute transfer) | vs the +0.0019 gate |
|---|---|---:|---:|---|
| `deatten_slope` | best measured | **+0.00017** | +0.00025 | **11× below** |
| `k13` | best tuned share | +0.00015 | +0.00026 | 13× below |
| `a_self` | **the mandate's `w ∝ a_i`** | **−0.00036** | −0.00005 | **negative** |
| `deatten` | **de-attenuate then average** | **−0.00137** | −0.00117 | **negative, CI excludes 0** |
| `minvar` | **`w ∝ a_i/σ_i²`** | **−0.00242** | −0.00211 | **negative, CI excludes 0** |
| `share_4_1_1` | 010's tuned optimum | −0.00320 | −0.00184 | negative |
| `a_self_pt` / `deatten_pt` | per-target | −0.01079 / −0.01118 | −0.00881 / −0.00926 | negative |
| best permutation falsifier | `perm_a_5` | −0.00011 | +0.00017 | ties the best real arm |

**"Equal weighting is optimal even given the 3.74× attenuation spread" is the result**, and it is
now supported on a mechanism rather than a fit:

* the spread is real and **larger on the panel than reported** (3.74×, 62.2% vs 1.5% clearing >80%);
* production's per-donor median-norm rescale **already absorbs any global attenuation exactly**
  (4.1e-8), so only relative weights are expressible;
* `a_i` **is not identified** as a scalar on the panel population — four independent estimates span
  2.5× and the one-factor model they assume is rejected at |z| up to 11 wherever it has power;
* the score's response to the weight vector is **97.8% a K562-share effect** and the attenuation
  hypothesis's clean prediction (the Orion split) is falsified on its own axis;
* the model's own predicted gain (**+8.5% to +16.8% cosine**) is excluded by the measurement
  (**−1.50%, CI [−2.38%, −0.67%]**);
* the permutation control ties the best real arm.

**No slot is requested and none should be spent on this.**

---

## 8. What this closes, and what it hands on

**Closed.** Donor mixing is not a lever on the direction channel, and the closure now covers the
one case AGENTS.md §4.2's equal-weight finding did not: a weight **derived** from a measured donor
property rather than tuned. §4.2's four numbered findings stand. The sentence in §4.2 that equal
weighting is best can now carry a second clause: *and it remains best when the donors' achieved
repression is measured and used to derive the weights, because the score does not respond along
that axis.* 194's hand-off is answered: the 3.3× (3.74× on the panel) repression difference is a
real and previously unrecorded quality difference between the three tables, and it is **not** a
reason to stop weighting them equally.

**Three things handed on.**

1. **The Orion↔Orion excess now has a significance and a shape.** 194 measured it as a
   target-independent blob (1.00% → 10.45% under knockdown weighting). This experiment shows the
   excess survives blob deflation and is therefore a **subspace**, and sizes it against a
   single-factor null at **z = +5.14 (5-line) / +7.56 (6-line)** — larger in relative terms than a
   same-cell-line pair. Anything that reads cross-line agreement between HCT116 and HEK293T as
   evidence about shared biology is reading that subspace. This includes 159's κ² route-2 estimate,
   whose HCT116↔HEK293T entry feeds the effective-donor-count arithmetic §8.2 uses.
2. **`a_i` is population-dependent, and the panel is the wrong population to have calibrated on.**
   The essential-gene population gives `a_rel(HEK293T) = 0.504 [0.44, 0.58]`; the panel population
   gives `1.254 [0.74, 1.88]`. §2.1 says the panel was built to be diverse and non-common-essential.
   Any future donor-level calibration must be estimated on panel-like targets or it is measuring a
   different quantity.
3. **The cosine is measurable to ±0.000044 absolute on PLOLO with a paired target bootstrap**
   (`_cos_per_target` in `scores204.json`, 4,000 draws; the tightest arm intervals, `k115` and
   `h115`, have half-width 4.3e−5 and 4.5e−5), which is **±0.20% relative**, i.e.
   **±0.00016 board cosine ≈ ±0.00006 overall** through §8.2's exchange rate. Given §8.2's finding
   that everything left is direction, that is a usable local instrument for direction arms — sharper
   for this purpose than PLOLO's `pds` channel (LOO RMS 0.00067 raw ≈ 0.00025 overall), and it is
   what turned this experiment's ambiguous `pds` ridge into a clean answer. **A direction arm should
   be priced on Δcos first and Δ`pds` second.**

**Not pursued, per the mandate's territory clause:** 193's basis composition, which it hands on as
the live lever on the operator, is noted and left alone.

---

## 9. Files

```
repress204.py   Phase 0: the closed-form repression table from the shipped pseudobulk, checked
                against 194's per-cell _kd_target.npy; operator + input sha table (rule 11)
prep204.py      production-frame donor effects (062's prepare_main_source), cached per target set
atten204.py     Phase 1: route A (self-gene residual), route B (3-donor Gram triple),
                route B' (naive slope), fourth-line instruments, blob geometry, with/without
                blob deflation
atten2_204.py   Phase 1b: the same with GENE COVERAGE controlled + the planted-attenuation
                known-answer row + target bootstraps
factor204.py    Phase 1c: one-factor (rank-1) test of the shared-component model on 4-6 line
                blocks, with the same-cell-line positive control
factor2_204.py  Phase 1d: target bootstrap of the rank-1 residuals, tetrads, and both the
                reliability (c) and the amplitude (a = c*||d||) scales
homog204.py     the scale identity: de-attenuating a donor table is an exact no-op (4.1e-8)
score204.py     Phase 2: PLOLO with a weighted `build_baseline`; per-target pds AND cosine
analyse204.py   the pds screen + the a15->a23 mass-fix known-answer row
six204.py       Phase 2/3: cosine bootstrap, six-member pricing, the mechanism regression
power204.py     Phase 3: the gain the minimum-variance model predicts for itself
work/           logs, scores204.json (120 fold x arm rows), *.json for every stage,
                plolo/ (162's fold caches, copied read-only)
```

**Memory.** Peak RSS **4.73 GiB** on the heaviest fold, measured with
`resource.getrusage(RUSAGE_CHILDREN)`. One earlier design loaded the 2.57 GiB operator and the CD4
banks *before* the donor streaming passes and tripped 062's own guard at **9.005 GiB**; the guard
reads *peak* RSS, so the fix is ordering (stream first, allocate the operator last) plus one process
per fold. Both are in `score204.py` with the comment.

**Hygiene.** Nothing outside this directory was written. `data/donors/` was read only. No `vcc`, no
`cell-eval2`, no build, no `submissions/`. `debench.Surface.pds()` was never called. The two Orion
folds are never quoted as absolutes. Every arm is deterministic: re-scoring `equal` on the
`orion_hek293t` fold in a fresh process reproduced `pds 0.777375` and `max_abs_dev_from_162_pre`
5.8e−08.

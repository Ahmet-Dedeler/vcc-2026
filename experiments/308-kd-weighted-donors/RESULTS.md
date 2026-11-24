# 308 — RESULTS: weighting donor rows by each donor's OWN measured on-target knockdown

**308 (|kd|-weighted donors): DOES NOT CLEAR** (best arm `soft:abs:p1.0:dis` `pds` **+0.00077
[−0.00114, +0.00283]** on **3/3** folds; vs the reliability control `relw` **+0.00447 [−0.00213, +0.01102]**
— not separated, and vs 097's norm control `normrel1.5` **−0.00182 [−0.00484, +0.00106]** — behind it;
`nmae` identity **+0.00074** (a cost, on the arm's own `pds` gain); priced **−0.00145 pessimistic /
−0.00109 point / +0.00003 optimistic vs +0.0089**). **0 of 16 arms clear; 51 cells screened.**
And the target-set-matched controls close it: on exactly the arm's own rows, entropy-matched **random**
weights read **+0.00101 / +0.00054 / +0.00100** — two of three *larger* than the arm.

**Date 2026-09-07**, HEAD `89d38f7`. Analysis only: **no build, no submission, no `vcc`, no `cell-eval2`,
no download**, nothing under `submissions/`, nothing written outside
`experiments/308-kd-weighted-donors/` (write audit at §10). `uv run python` throughout,
`PYTHONDONTWRITEBYTECODE=1`; heavy jobs serial. **Peak RSS 5.12 GiB** (fold sweeps), 8.68 GiB (the
feature pass's memmap scan of the two 2.7 GB Orion `delta.npy` files), well under the 10 GB cap.
Fold caches are byte copies of 262's (`work/plolo/SHA256SUMS.txt`, `diff` against 262's list empty).
**Every known-answer row passes — 64 of 64: 23 on the feature and the scorer, 3 × 12 on the PLOLO
scoring path, 4 on the h1-test instrument, 1 retrodiction — all at 0.0e+00 except three at ≤ 5e-12.**

---

## 0. The answer in eight lines

1. **The feature exists, is build-time observable, and is clean.** |kd|(t, d) — target *t*'s own gene
   coordinate in donor screen *d*'s own measured response to perturbing *t*, in *d*'s own
   `log1p(50,000·p)` units — is defined for **226 / 290 / 294 of the 300 panel targets** in K562-GWPS /
   Orion-HCT116 / Orion-HEK293T. It is read from one donor bank's `delta.npy` and nothing else: no truth,
   no fold, no operator, no held bank (§2, KA-D at 0.0e+00). It is *absent from the incumbent's base* by
   construction (`E062.prepare_main_source` zeroes every panel column), so it is genuinely new information.
2. **The set the arm can act on is large**, so nothing below is a power problem: **106 of 300 panel
   targets disagree** between donors at 296's library-relative cut (67 at 293's absolute cut); per PLOLO
   fold **30 / 87 / 93** (rel) and **20 / 63 / 49** (abs); on the h1-test Flex reference **34 / 85** and
   **21 / 85**.
3. **The hard gate (arm 1) loses.** `gate:rel` **−0.00863 [−0.01413, −0.00335]**, 0/3 folds;
   `gate:abs` **−0.00320 [−0.00735, +0.00077]**, 1/3. Arm 3 is *identical* to arm 1 for a hard gate — a
   structural fact, not a result: a gate can only act where donors disagree (§4).
4. **The soft weight (arm 2) loses monotonically in its parameter**, on both cuts: at the relative cut
   −0.0023 / −0.0060 / −0.0124 / −0.0189 for p = 0.5 / 1 / 2 / 4; at the absolute cut −0.0004 / −0.0028 /
   −0.0069 / −0.0124. The dose curve points the wrong way at every rung.
5. **The one cell with 3/3 positive folds is `soft:abs:p1.0:dis` at +0.00077 [−0.00114, +0.00283]** — an
   interval through zero, retention −151 % (three targets carry it), and **indistinguishable from randomly
   re-weighting the same rows** (§6). It is the act of perturbing those rows, not |kd|.
6. **The decisive control says the mechanism is refuted at this construction.** |kd| does not beat
   reliability weighting by a separated margin anywhere, and it is *behind* 097 route (a)'s
   `normrel^−1.5` on every arm (`gate:abs − normrel1.5` **−0.00579 [−0.01042, −0.00132]**). `normrel1.5`
   is the only positive row in the whole screen (+0.00259 [+0.0000, +0.0053], 3/3, priced +0.00104
   pessimistic = 0.55× the shape-only gate) — and it too does not clear.
7. **Why, measured, not argued.** On the disagreeing set the higher-|kd| donor's row is *not* the one
   closer to the held truth: dropping the failing donor moves **more targets DOWN than UP in all six
   (cut × fold) cells, by 1.5–2.4×** (rel/HEK293T 25 up vs 53 down, mean Δ`pds` −0.031), and the mean
   Δ`pds` on the set is negative in five of the six (the exception is abs/K562, +0.010 on 20 targets).
   And |kd| is **not** a norm
   proxy (Spearman of the log |kd| ratio with the log row-norm ratio +0.03 / +0.03 / −0.005) and only
   weakly a reliability proxy (+0.04 … +0.19): it carries information the existing keys do not, and that
   information does not help on this metric.
8. **The scorer-exclusion distinction the mandate asked me to verify holds, and it is not what killed
   this.** `_drop_on_target` / `_exclusion_cols` remove each row's own-target coordinate from every member
   and `pds` drops all 300 panel genes — verified numerically in the installed `cell_eval2` 0.16.0
   (KA-C). 097 route (d) is closed *as a prediction target*; using |kd| as a *weight* was open, was
   measured here, and is now closed at this construction too.

**One line: 308 (|kd|-weighted donors): DOES NOT CLEAR (best arm `soft:abs:p1.0:dis` `pds` +0.00077
[−0.00114, +0.00283] on 3/3 folds; vs reliability control `relw` +0.00447 [−0.00213, +0.01102], not
separated, and behind `normrel1.5` by −0.00182 [−0.00484, +0.00106]; `nmae` identity +0.00074 (a cost);
priced −0.00145 pessimistic vs +0.0089).**

---

## 1. Known-answer rows — 64/64 PASS (`work/ka308.{py,json,log}`, the sweeps' own `ka` blocks)

### 1a. The scoring path (`work/sweep308-<fold>.json`, 12 rows per fold, all three folds)

| row | required | measured |
|---|---|---|
| **KA-0** `E062.build_baseline` over the rebuilt per-donor rows, float32, == the fold cache's `pre` | 0 | **0.0e+00** on all three folds |
| **KA-0b** votes == per-donor availability | 0 | **0.0e+00** ×3 |
| **KA-1** `cd4svd5` == 262 (= 257 k5_s1 = 261 R1 = 271) on `cos_agg` / macro / `pds` | 1e-9 | **0.0e+00** ×3 — `pds` **0.815384 / 0.801659 / 0.780978** exactly as the mandate specifies |
| **KA-1b** `cd4svd5` per-target `pds` array bit-identical to 262's | 0 | **0.0e+00** ×3 |
| **KA-1c** `cd4svd5` `nmae_wald89` == 271's | 0 | **0.0e+00** ×3 |
| **KA-2** `mean_explicit` (the equal-weight mean rebuilt through the selection path) == `cd4svd5` | 1e-9 | **0.0e+00** ×3 |
| **KA-2b** `soft:*:p0` (equal weights through the *weight* path) == `cd4svd5` | 1e-9 | **5.0e-12 / 2.9e-12 / 1.9e-12** |
| **KA-C** each donor's production row == 271's independently built half-builder row (float32) | 0 | **0.0e+00** ×2 per fold |
| **KA-3** `relw` == 271's published `relw` on `cos_agg` / macro / `pds` | 1e-9 | **0.0e+00** ×3 (0.808804 / 0.799049 / 0.779066) |
| **KA-3b** `relw_glob` == 271's | 1e-9 | **0.0e+00** ×3 |

### 1b. The feature and the scorer (`work/ka308.json`)

| row | required | measured |
|---|---|---|
| **KA-A1** 293 §4's per-target `kd big` / `kd small`, all 16 two-construct targets, re-derived by this code from `construct293.csv` | printed 3 d.p. | **16/16, max abs 0.0e+00** |
| **KA-A2/3/4** 293 §5b's Arc-2025 panel rows from 293's `kd_all` | 1.485 / 18.7 % / 4.3 % | **1.48485 / 18.667 % / 4.333 %** |
| **KA-B** the chunked gather == a direct element-wise memmap read, all three banks | 0 | **0.0e+00** (226 / 290 / 294 pairs) |
| **KA-C1** `cell_eval2.metrics.delta` defines `_exclusion_cols` and `_drop_on_target` | present | present (`cell_eval2` 0.16.0) |
| **KA-C2** `run._ONTARGET_EXCLUDING_FUNCS` covers `mse_unbiased_capped`, `de_sig_jaccard`, `de_lfc_nmae` | yes | yes (also `mse_unbiased`, `distance_unbiased`) |
| **KA-C3** `metrics.direction._exclude_on_target` (fidelity / reach) | present | present |
| **KA-C4** `metrics.discrimination` takes `exclude_target_gene` with `exclusion_scope="panel"` | present | present |
| **KA-C5** **numeric**: `_exclusion_cols(["G1","G2"], ["G0","G1","G2","G3"], target_gene_map={"G1":"G1","G2":"G2"})` | row 0 → col 1, row 1 → col 2 | `{'0': [1], '1': [2]}` — **the own-gene column, exactly** |
| **KA-C5b** with an EMPTY map the resolver still matches each label to its own symbol and excludes it | same 2 columns | same 2 columns — the exclusion is the default |
| **KA-C5c** an unresolvable map is REFUSED, not ignored | raises | `ValueError: "…exclusion is part of the metric and cannot be switched off"` |
| **KA-D1/2** LEAKAGE: the two surviving donors' |kd| computed in a process that never opens the held bank == the panel-wide array, identical NaN masks | 0 | **0.0e+00** for all three held lines |
| **KA-E** the K562-GWPS bank-derived |kd| vs 296's independently cell-derived `kd_single` on 218 shared panel targets | Spearman > +0.8 | **Spearman +1.000, Pearson +1.000**, median \|kd\| 0.762 vs 0.762, < 1.0 70.2 % vs 296's published 70.1 % |

### 1c. The h1-test instrument and the price

| row | required | measured |
|---|---|---|
| **KA-F** my h1-test `mean:rebuilt` / `relw` / `relw:rev` / `ebrow` == 273's own published rows on `pds`, `cos`, `cos_agg` | 0 | **0.0e+00** on all four (`mean:rebuilt` 0.882213, `relw` 0.876611, `relw:rev` 0.842577, `ebrow` 0.882213) |
| **KA-G** the retrodiction, built inside this experiment (`stack308.py`'s `STACK` + the sweep's `cd4svd5`) through 261's recipe: BR249 | 257's R1 to 1e-6 | Δ`pds` **+0.012263** (257 +0.012263), agg **−0.4617 %** (257 −0.4617 %), Δ`nmae89` **−0.00138** (257 −0.00138), max abs **0.0e+00**; priced POINT `pds` +0.00552, `mse` −0.00083, DE3 −0.00035, `nmae` +0.00038 = **+0.00472** vs 257's R1 +0.00472 and the **board's +0.0042634** |

---

## 2. The feature: provenance, distribution, missingness (`kd308.py`, `work/kd308_feature.{npz,json}`)

**Definition.** `kd(t, d) = delta_d[row of t, column of gene t]`, where `delta_d[p, g] =
log1p(50,000·pert_frac) − log1p(50,000·ctrl_frac)` is donor screen *d*'s own pseudobulk perturbation
table (`scripts/build_donor_table.py:10`, `scripts/build_orion_donor.py:112-116`). That is **exactly
293's `|kd|` scale**, which is why KA-E reads Spearman +1.000 against 296's independent extraction.

**Provenance of every array, stated as the mandate requires.**

| array | source | touches held truth? |
|---|---|---|
| `kd`, `abs_kd`, `measured` | `data/donors/<d>/{delta.npy, genes.txt, perts.txt}` — one donor bank at a time | **no** |
| `cut_rel` (296's C2) | the same banks: the median \|kd\| over every perturbation whose own gene is on the bank's axis | **no** |
| `cut_abs` (293's C1) | the constant 1.0 | **no** |
| `cells` (depth covariate) | `experiments/278-depth-inventory/work/depth278.csv`, read-only | **no** |
| donor rows, operator, truth | 262's fold caches (byte copies) — the *incumbent's* objects, unchanged | truth only as the scoring target, never as a feature input |

**The leakage row.** KA-D recomputes the feature for the two surviving donors of each fold in a process
that never opens the held line's bank: bit-identical to the panel-wide array, identical missingness,
on all three folds. The production prepared rows set `effects[:, panel_positions] = 0.0`, so this
coordinate is not in `cd4svd5`'s base at all — it had to be read from the bank, and it was.

**Distribution and missingness (300 panel targets).**

| donor | perturbed | own gene on the bank axis | **feature defined** | median \|kd\| (panel) | \|kd\| < 1.0 | < 0.5 | lfc > 0 | library-wide median \|kd\| (**296 C2 cut**) | library \|kd\| < 1.0 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `k562_gwps` (8,248 genes) | 269 | 231 | **226** | 0.746 | 71.2 % | 25.2 % | 0.0 % | **0.6571** (7,617 self-measured perts) | 76.5 % |
| `orion_hct116` (38,606) | 290 | 300 | **290** | 0.476 | 92.1 % | 52.4 % | 0.3 % | **0.1563** (17,767) | 95.4 % |
| `orion_hek293t` (38,606) | 294 | 300 | **294** | 0.363 | 99.3 % | 75.5 % | 0.3 % | **0.1476** (17,977) | 99.5 % |

**296's warning reproduces exactly**: the absolute cut is *not* comparable across libraries. On Orion's
38,606-gene axis the `log1p` CPM scale is compressed and 293's |lfc| ≥ 1.0 passes only **8 %** (HCT116)
and **1 %** (HEK293T) of panel targets against **29 %** in K562-GWPS. Both cuts are carried below; the
library-relative one is the transferable key and gets the wider disagreeing set.

**Depth (programme rule 9).** Per-target cell counts from 278: median 174 / 153 / 202 cells,
Spearman(cells, |kd|) **+0.04 / −0.11 / −0.16** — |kd| is essentially not a depth statistic. The arms are
also re-read inside depth terciles in §7.

### 2a. The disagreeing set — the bound on everything below

| cut | panel-wide (of 300) | k562 fold | hct116 fold | hek293t fold | h1-test (of 85) |
|---|---:|---:|---:|---:|---:|
| **abs** (293, \|kd\| ≥ 1.0) | **67** disagree, 231 all-fail, 2 all-work | 20 | 63 | 49 | 21 |
| **rel** (296 C2, per library) | **106** disagree, 24 all-fail, 170 all-work | 30 | 87 | 93 | 34 |

The set is a third of the panel at the relative cut; the negative below is not a power failure. The
**gate's fallback rate** (all carrying donors below the cut, so the incumbent's row stands) is 237 / 15
per fold at the abs / rel cut on the K562 fold, i.e. the absolute cut is nearly inert on Orion donors.

---

## 3. What was run — 51 cells × 3 folds + 47 on h1-test + 8 matched cells × 3 folds

Every row rides `cd4svd5`'s downstream verbatim (lib263: whitening through `plolo.Operator.apply`,
CD4 deflated k = 5 re-solved at f = 0.34, the H1-pair union re-solved by the carrying-row rule at
g = 0.28, one add, one preserve), scored by `lib246.score` with 188's exact `nmae` identity. Weighted
rows are **norm-matched** to the cached mean row (263's `:nm`, the shape-only reading); `:lit` is the
literal amplitude-changing form.

* **Arm 1 `gate:<cut>`** — drop donor *d* for target *t* when |kd|(t,d) < cut_d, renormalise over the
  survivors; if all fail, fall back to the equal-weight mean.
* **Arm 2 `soft:<cut>:p<P>`** — w_j ∝ (max(|kd|_j, 1e−3) / cut_j)^P, P ∈ {0.5, 1, 2, 4} (P = 0 is the
  incumbent, a KA). Standardising by the donor's own cut makes the weight scale-free across libraries.
* **Arm 3 `<arm>:dis`** — the same restricted to the disagreeing set, incumbent elsewhere.
* **Controls** — `cd4svd5`, `mean_explicit`, `relw` (271's reliability weights, w ∝ R/(1−R)),
  `relw_glob` (271's donor-level scalar reliabilities, the 019/142 class), `normrel1.5`
  (097 route (a) a6, w ∝ ‖row_j(t)‖^−1.5 relative to the donor's own median), `:shuf` (|kd| permuted
  across targets within donor, 3 seeds), `:rand` (weights with the arm's realised entropy, permuted over
  donors, 3 seeds), `:rev` (weight the FAILING construct more).

---

## 4. The PLOLO table (three cancer folds, paired target bootstrap B = 20,000, seed 20260907)

Incumbent `pds` 0.815384 / 0.801659 / 0.780978, `nmae` identity 0.981452.
`nmae` is reported as the **change in the raw identity — positive is a COST**.

| cell | pooled Δ`pds` [95 % CI] | per-fold | folds > 0 | Δcos_agg % | **Δ`nmae`** | pessimistic price |
|---|---|---|---:|---:|---:|---:|
| `mean_explicit` / `soft:*:p0` (KA) | +0.00000 [+0.00000, +0.00000] | 0 0 0 | — | +0.00 | +0.00000 | — |
| **`gate:abs`** (= `gate:abs:dis`) | −0.00320 [−0.00735, +0.00077] | +0.0007 −0.0043 −0.0061 | 1/3 | −4.49 | **+0.00125** | −0.01243 |
| **`gate:rel`** (= `gate:rel:dis`) | **−0.00863 [−0.01413, −0.00335]** | −0.0086 −0.0075 −0.0098 | 0/3 | −3.20 | −0.00070 | −0.01157 |
| `soft:abs:p0.5` | −0.00035 [−0.00342, +0.00278] | −0.0015 −0.0009 +0.0014 | 1/3 | −0.73 | +0.00063 | −0.00286 |
| `soft:abs:p1.0` | −0.00277 [−0.00762, +0.00204] | −0.0030 −0.0048 −0.0005 | 0/3 | −2.10 | +0.00174 | −0.00836 |
| `soft:abs:p2.0` | −0.00685 [−0.01353, −0.00022] | | 0/3 | −4.70 | +0.00319 | −0.01797 |
| `soft:abs:p4.0` | −0.01243 [−0.02059, −0.00429] | | 0/3 | −8.26 | +0.00357 | −0.02975 |
| **`soft:abs:p1.0:dis`** (best) | **+0.00077 [−0.00114, +0.00283]** | +0.0008 +0.0013 +0.0002 | **3/3** | −0.71 | **+0.00074** | **−0.00145** |
| `soft:abs:p2.0:dis` | +0.00024 [−0.00253, +0.00306] | | 1/3 | −1.75 | +0.00109 | −0.00447 |
| `soft:rel:p0.5` | −0.00233 [−0.00484, +0.00022] | | 0/3 | −0.78 | −0.00142 | −0.00327 |
| `soft:rel:p1.0` | −0.00596 [−0.01026, −0.00166] | | 0/3 | −2.75 | −0.00074 | −0.00961 |
| `soft:rel:p2.0` | −0.01236 [−0.01846, −0.00625] | | 0/3 | −6.54 | −0.00016 | −0.02143 |
| `soft:rel:p4.0` | −0.01892 [−0.02639, −0.01142] | | 0/3 | −10.39 | +0.00051 | −0.03359 |
| `soft:rel:p1.0:dis` | −0.00270 [−0.00578, +0.00035] | | 0/3 | −1.01 | −0.00095 | −0.00378 |
| `soft:rel:p2.0:dis` | −0.00519 [−0.00948, −0.00095] | | 0/3 | −2.11 | −0.00082 | −0.00745 |
| **CONTROLS** | | | | | | |
| `relw` (271, the decisive control) | −0.00370 [−0.01030, +0.00293] | −0.0066 −0.0026 −0.0019 | 0/3 | −3.23 | +0.00490 | −0.01117 |
| `relw_glob` (019/142 class) | −0.00189 [−0.00625, +0.00243] | | 1/3 | −3.30 | +0.00469 | −0.01298 |
| **`normrel1.5`** (097 a6) | **+0.00259 [+0.00000, +0.00527]** | +0.0006 +0.0032 +0.0040 | **3/3** | +0.12 | **−0.00095** | **+0.00104** |
| `gate:abs:rev` | −0.00374 [−0.00746, −0.00003] | | 0/3 | −3.25 | −0.00043 | −0.00942 |
| `gate:rel:rev` | −0.00019 [−0.00505, +0.00486] | | 2/3 | −3.67 | +0.00213 | −0.00936 |
| `gate:abs:shuf:s0/1/2` | −0.0034 / −0.0089 / −0.0058 | | 0/3 each | −4.2 … −4.5 | +0.0013 / +0.0035 / +0.0009 | −0.012 … −0.016 |
| `gate:rel:shuf:s0/1/2` | −0.0102 / −0.0090 / −0.0071 | | 0/3 each | −4.2 … −5.0 | −0.0001 / +0.0025 / −0.0005 | −0.013 … −0.017 |
| `gate:abs:rand:s0/1/2` | −0.0019 / −0.0020 / −0.0030 | | 0–1/3 | −3.1 … −4.3 | +0.0010 / −0.0001 / +0.0010 | −0.008 … −0.012 |
| `gate:rel:rand:s0/1/2` | −0.0029 / −0.0011 / −0.0047 | | 1–2/3 | −3.3 … −4.1 | ≈ 0 | −0.009 … −0.012 |
| `gate:rel:lit` / `soft:rel:p1.0:lit` / `gate:abs:lit` | −0.0098 / −0.0062 / −0.0040 | | 0/3 each | −6.4 / −3.0 / −5.0 | +0.0027 / +0.0012 / **+0.0048** | −0.029 / −0.012 / −0.020 |

Two structural facts in that table, both worth stating plainly:

* **Arm 3 ≡ arm 1 for a hard gate** (`gate:*` and `gate:*:dis` are bit-identical on every fold): a gate
  can only change a row where at least one donor passes and at least one fails, which *is* the
  disagreeing set. Arm 3 is a real, separate arm only for the soft weight.
* **`gate:rel`'s reverse control beats it.** `gate:rel` −0.0086 against `gate:rel:rev` −0.0002. Keeping
  the *failing* construct's row costs less than keeping the working one. That is the sign the mechanism
  forbids.

### 4a. The decisive contrasts (`arm − control`, paired over targets, pooled)

| contrast | pooled | folds > 0 |
|---|---|---:|
| `gate:abs` − `relw` | +0.00050 [−0.00660, +0.00764] | 1/3 |
| `gate:abs` − `relw_glob` | −0.00131 [−0.00590, +0.00336] | 2/3 |
| **`gate:abs` − `normrel1.5`** | **−0.00579 [−0.01042, −0.00132]** | 1/3 |
| `gate:rel` − `relw` | −0.00493 [−0.01323, +0.00315] | 0/3 |
| `gate:rel` − `normrel1.5` | −0.01122 [−0.01718, −0.00537] | 0/3 |
| `soft:abs:p1.0:dis` − `relw` | +0.00447 [−0.00213, +0.01102] | 3/3 |
| `soft:abs:p1.0:dis` − `relw_glob` | +0.00266 [−0.00115, +0.00661] | 2/3 |
| **`soft:abs:p1.0:dis` − `normrel1.5`** | **−0.00182 [−0.00484, +0.00106]** | 1/3 |

**Reading.** No arm beats the reliability control by a bootstrap-separated margin (every interval spans
zero), and every arm is at or behind the *norm* control, decisively so for the gates. **|kd| adds nothing
over reliability, and it is worse than the key 097 already had.** By the mandate's own test — "if |kd|
adds nothing over reliability, the mechanism above is refuted and you say so plainly" — the mechanism is
refuted at this construction.

---

## 5. The `nmae` leg (the killer, BR237/ARM6)

Reported for every cell above as Δ of 188's exact identity concatenated over the three folds (positive =
cost). The best arm `soft:abs:p1.0:dis` costs **+0.00074** — small, but it is a cost on the arm whose
whole `pds` gain is +0.00077, which is why its six-member price is negative even where `pds` is nominally
positive. The literal (amplitude-changing) readings are worse: `gate:abs:lit` **+0.00484**, i.e. the same
order as ARM6's board-realised −0.00716 leg. `normrel1.5` is the only cell with an `nmae` *gain*
(−0.00095) — consistent with 097's finding that the direction-only leg is where its 1.035× lives.
**No arm in this experiment buys `pds` and keeps `nmae`.**

---

## 6. The only 3/3 cell, killed by its own matched control (`matched308.py`, `work/matched_analysis.log`)

`soft:abs:p1.0:dis` is the single cell with Δ`pds` > 0 on all three folds. The §4 controls are
**unrestricted** (they re-weight all ~259 two-donor rows) while this arm touches only 20 / 63 / 49 rows,
so comparing them is unfair to the control. Re-run with every control acting on **exactly the arm's own
row set** (the real disagreeing gate held fixed; only the weight inside it changes):

| row (matched to the arm's own 20/63/49 rows) | pooled Δ`pds` [95 % CI] | folds > 0 | Δ`nmae` |
|---|---|---:|---:|
| **`soft:abs:p1.0:dis`** (the arm) | **+0.00077 [−0.00115, +0.00285]** | 3/3 | +0.00074 |
| `:rev` (weight the failing construct) | −0.00112 [−0.00318, +0.00085] | 1/3 | −0.00031 |
| `:shuf:s0` / `s1` / `s2` (|kd| permuted within donor) | −0.00080 / +0.00017 / +0.00052 | 1/3, 1/3, 2/3 | +0.0005 … +0.0008 |
| **`:rand:s0` / `s1` / `s2`** (entropy-matched random) | **+0.00101 / +0.00054 / +0.00100** | 2/3, 1/3, 2/3 | +0.0000 … +0.0004 |
| arm − `:rand:s0` / `s1` / `s2` | −0.00025 / +0.00023 / −0.00023, all intervals through zero | | |

**Two of three entropy-matched random controls are larger than the arm.** The +0.00077 is what you get
from softening the equal-weight mean on those rows *at all*; |kd| contributes nothing identifiable. The
reverse control being only −0.0011 below (and positive on one fold) says the same. Retention (leave-top-3-out)
is **−151 %**: three targets carry the whole thing.

---

## 7. Depth conditioning (programme rule 9) — no depth stratum rescues it

`gate:<cut>` re-read inside terciles of the minimum cells/target over the fold's donors (278's counts),
on the disagreeing set only:

| cut | fold | T1 shallow | T2 | T3 deep |
|---|---|---:|---:|---:|
| abs | k562 | −0.0505 (7) | +0.0199 (7) | +0.0692 (6) |
| abs | hct116 | −0.0692 (21) | +0.0233 (21) | −0.0138 (21) |
| abs | hek293t | −0.0719 (17) | +0.0261 (16) | −0.0607 (16) |
| rel | k562 | −0.0928 (10) | −0.0978 (10) | −0.0439 (10) |
| rel | hct116 | −0.0049 (29) | −0.0430 (29) | −0.0268 (29) |
| rel | hek293t | −0.0040 (31) | −0.0356 (32) | −0.0534 (30) |

At the **relative** cut the arm is negative in **9 of 9** strata. At the absolute cut the shallow tercile
is the worst on all three folds (−0.050 … −0.072) and the mid tercile positive on all three — a pattern
worth one line but not a result: these are 18 post-hoc cells of ≤ 21 targets each on truths whose own
per-target readings 162/267 label noise-dominated. No conditioning turns the arm positive.

---

## 8. The h1-test Flex-400 reference (273's instrument, union OFF, 85 targets)

Built through 273's own item-(b) path; `mean:rebuilt`, `relw`, `relw:rev`, `ebrow` reproduce 273's
published rows at **0.0e+00** (KA-F). Incumbent `pds` **0.882213**, `cos_agg` 0.070681. Disagreeing set
21 (abs) / 34 (rel) of 85; 80 of 85 targets carry all three donors here.

| arm | Δ`pds` | paired 95 % CI | Δcos_agg % |
|---|---:|---|---:|
| `gate:abs` (= `:dis`) | **−0.00140** | [−0.00910, +0.00560] | −3.13 |
| `gate:rel` (= `:dis`) | **−0.01919** | [−0.03810, −0.00434] | −5.89 |
| `soft:abs:p0.5` / `p1.0` / `p2.0` / `p4.0` | −0.0032 / −0.0085 / −0.0172 / −0.0268 | `p1.0` [−0.01709, −0.00042] | −0.1 … −10.3 |
| `soft:abs:p1.0:dis` / `p2.0:dis` (the PLOLO best arm) | **+0.00042** / −0.00182 | [−0.00364, +0.00518] | +0.41 / −0.09 |
| `soft:rel:p0.5` / `p1.0` / `p2.0` / `p4.0` | −0.0048 / −0.0104 / −0.0198 / −0.0333 | | −2.2 … −20.2 |
| `soft:rel:p1.0:dis` / `p2.0:dis` | −0.0063 / −0.0116 | | −1.7 / −4.5 |
| `relw` (control) | −0.00560 | [−0.02493, +0.01513] | +0.89 |
| `normrel1.5` (control) | −0.00420 | [−0.01457, +0.00476] | +3.81 |
| `relw:rev` (271's reversed reliabilities) | −0.03964 | [−0.07479, −0.00602] | −19.78 |

**The second instrument agrees**: every |kd| arm is at or below zero, the relative-cut gate is decisively
negative (its interval excludes zero), and the best PLOLO arm reads +0.0004 — inside its own noise on a
reference that resolves every target. The transfer cosine reading the mandate asked for: every |kd| arm
moves `cos_agg` **down**, by −0.1 % (`soft:abs:p0.5`) to −20.2 % (`soft:rel:p4.0`), on the instrument
where per-target cosine is measurable; only `normrel1.5` moves it up (+3.8 %) and it still loses `pds`.

---

## 9. The price (AGENTS §8.4, `price308.py`, `work/price308.{json,log}`)

261's recipe imported verbatim (five-member propagation; coherent-pessimistic convention: per leg the end
of its measured band that *lowers* the total for that cell's sign; class factors k_pds {0.85, 1.22},
k_cos {1.36, 3.76}, DE-three {0.0425, 0.24}, identity ×1.5 on a loss / ×0 on a gain). The retrodiction
row (KA-G) reproduces BR249 at +0.00472 (board +0.0042634) before anything is priced.

**Price stated BEFORE the verdict sentence:**

| cell | PESSIMISTIC | POINT | OPTIMISTIC | class | its gate | ×gate | retention |
|---|---:|---:|---:|---|---:|---:|---:|
| **`soft:abs:p1.0:dis`** (best arm) | **−0.00145** | −0.00109 | +0.00003 | amplitude (dlv +0.87 %) | +0.0089 | **−0.16×** | −151 % |
| `soft:abs:p0.5` | −0.00286 | −0.00262 | −0.00142 | amplitude | +0.0089 | −0.32× | 841 % |
| `gate:abs` | −0.01243 | −0.01130 | −0.00380 | amplitude | +0.0089 | −1.40× | 198 % |
| `gate:rel` | −0.01157 | −0.01070 | −0.00455 | amplitude | +0.0089 | −1.30× | 143 % |
| `soft:rel:p4.0` | −0.03359 | −0.03130 | −0.01275 | amplitude | +0.0089 | −3.77× | 124 % |
| *(control)* `relw` | −0.01117 | −0.00980 | −0.00422 | amplitude | +0.0089 | −1.26× | 223 % |
| *(control)* **`normrel1.5`** | **+0.00104** | +0.00184 | +0.00186 | **shape-only** (dlv +0.46 %) | +0.0019 | **0.55×** | 30 % |

**0 of 16 arms clear; 51 cells screened on PLOLO (16 arms + 35 controls/variants), 47 on h1-test, 8
matched cells** (each beside its own incumbent row). The best arm prices
**−0.00145 pessimistic against the +0.0089 amplitude gate** (its optimistic scenario, +0.00003, is still
below the shape-only gate it does not qualify for). 097 route (a)'s `normrel^−1.5` is the only positive
price in the file at 0.55× the shape-only gate — a reproduction of 097's +0.0008…+0.0009 reading on
today's base, not a candidate, and not this experiment's arm.

**Multiplicity.** 51 PLOLO cells (16 arms + 35 controls/variants) × 3 folds, 47 h1-test cells, 8 matched
cells, 18 post-hoc depth cells. One arm reached 3/3 folds with an interval through zero; under 51 cells
that is the expected count, and its matched control (§6) explains it without |kd|.

---

## 9a. The distinction the mandate asked me to state, and where it now stands

097 route (d) closed the on-target knockdown **as a prediction target**: the scorer removes each row's own
gene from every member (`_drop_on_target` / `_exclusion_cols`; `pds` drops all 300 panel genes at
`exclusion_scope="panel"`). I verified that claim myself, numerically, in the installed `cell_eval2`
0.16.0 (KA-C5: the resolver maps row *i* to exactly its own gene's column; KA-C5c: an unresolvable map is
refused with "exclusion is part of the metric and cannot be switched off"). **That closure is real and it
is not what this experiment tested** — a weight on which donor to believe changes the *direction* of the
17.9 k scored coordinates, which the exclusion does not touch.

**But the family is not as open as the mandate assumed, and I should flag it:** 097 route (a) already
screened `a6 selfkd` — "on-target knockdown depth: `max(0, −s_n(p)[gene p])`, a pure QC signal"
(`experiments/097-direction/w_arms.py:35`) — as a per-perturbation donor weight at exponents −4 … +4 plus
a hard argmax, and recorded **0.981 (h1-test) / 1.027 (orion), "dead cross-surface"**. That is the same
quantity in the same role, on the a15-era pipeline and the older harnesses. **308 is a re-measurement of
that family on today's incumbent with better instruments, a full six-member price, and matched controls —
and it reproduces 097's verdict with the sign resolved: not merely inconsistent across surfaces, but
negative on both of today's instruments.**

---

## 10. Limits — read before quoting

1. **Scope of the closure (programme rule 6).** What is closed is *this construction*: |kd| read as a
   bulk pseudobulk lfc on the target's own gene column, used as a hard gate or a one-parameter monotone
   soft weight over the **three production donor screens**, on `cd4svd5`'s base, scored on three PLOLO
   cancer folds and the h1-test Flex-400 reference. Not closed: |kd| as an input to something that is not
   a donor weight; |kd| for donors that do not exist yet; a per-cell (rather than pseudobulk) efficiency.
2. **The absolute cut is nearly inert on the Orion donors** (1–8 % pass) because their 38,606-gene axis
   compresses the `log1p` scale — 296's point. The relative cut is the informative one, and it is the one
   that loses hardest.
3. **The measurement is `pds`-first.** PLOLO's truths are 150–180-cell references whose per-target
   readings are noise-dominated (R267 §3); only totals and group aggregates carry signal, which is why
   every number here is a fold total or a bootstrap over targets. The h1-test reference resolves targets
   and agrees.
4. **A donor's |kd| is a property of that donor's own construct in that own cell line.** Nothing here
   measures the *challenge* library's constructs — no source on disk carries the 2026 constructs
   (READ-FIRST §6). F14/293's term binds the truths; this experiment could only ever act on the donors.
5. **`soft:abs:p1.0:dis` is classified amplitude-changing** (delivered-norm ratio +0.87 %, above the
   0.5 % shape bound) even though its rows are norm-matched: the global preserve moves when the row set
   changes. Its shape-only reading does not exist as built.
6. **The depth terciles in §7 are post-hoc** and small (6–32 targets).
7. **`relw` here is 271's estimator on 271's half rows** (read-only), reproduced at 0.0e+00; it is the
   published reliability control, not a re-derivation from cells.

---

## 11. What this means for the programme

* **Nothing to build. No slot. No candidate.** 0/16 arms clear; the best prices −0.00145 pessimistic
  against +0.0089.
* **A finding for FINDINGS, and it is about F14's term.** F14/293/296 established that the guide
  construct — largely *guide efficacy* — is the dominant irreducible term of a perturbation response, and
  that a failed construct still yields a reproducible transcriptome. 308 measures the competition-facing
  corollary and it is negative: **a donor's own measured on-target knockdown does not identify which
  donor's response transfers to a held cell line.** On the disagreeing set the higher-|kd| donor's row is
  the one *farther* from the held truth 1.5–2.4× as often as it is nearer, in all six cut × fold cells.
  The construct term is a ceiling on the *truths*; it is not a usable selector among *donors*.
* **It also isolates what reliability weighting was and was not blind to.** The mandate's premise was that
  reliability keys are structurally blind to a failed construct. That premise is *correct* — |kd| is
  nearly orthogonal to both the split-half reliability (ρ +0.04 … +0.19) and the row norm (ρ ≈ 0.00) — and
  the conclusion still does not follow: the information |kd| carries beyond those keys does not move the
  transfer. Adding a non-redundant, mechanism-motivated key to the donor-weighting family did not rescue
  it. Together with 016, 097(a), 262/263/264, 271, 299 and 300, **per-target donor re-weighting on the
  three production donors is now closed on reliability, agreement, norm, arbiter, spread AND measured
  on-target knockdown.** The equal-weight mean survives every key tried.
* **The one open thread it leaves** is not a knob: 263's oracle ceiling (+0.046 `pds` for the same
  drop-one operation) is still there and still needs an arbiter that knows the held line. |kd| is one more
  build-time-observable signal that is at chance for that job.

---

## Files

`kd308.py` (the feature, sha256 `d18b2ca066309e98fe994274f5191bc006d59cac948a0fe911f43f4b6af2b268`) ·
`lib308.py` (the weight constructions, frozen before any fold was scored — sha256
`485c0018…` at the three PLOLO sweeps, recorded in each `sweep308-<fold>.json`;
`8bc01663a47271af0ce1d7708ef887c59d8e0aaf68e2431072cbd9e4e90fb061` after the h1-test target-list helper
`kd_for_targets` was appended, which changed no weight construction) · `ka308.py` (23 KA rows) ·
`sweep308.py` (the three folds) · `stack308.py` (the retrodiction row) · `matched308.py` (the matched
controls) · `flex308.py` (the h1-test instrument) · `analyse308.py` · `price308.py` ·
`work/kd308_feature.{npz,json}` · `work/ka308.{json,log}` · `work/sweep308-<fold>.{json,log,nohup}` ·
`work/pertarget308-<fold>.npz`, `work/identity308-<fold>.npz`, `work/weights308-<fold>.npz` ·
`work/stack308-<fold>.json`, `work/stackid308-*.npz` · `work/matched308-<fold>.json`,
`work/matched_analysis.{py,json,log}` · `work/flex308.{json,log}`, `work/flex308_pertarget.npz` ·
`work/analyse308.{json,log}` · `work/price308.{json,log}` · **`work/kd308.csv`** (1,706 rows =
target × donor × fold: `kd_lfc`, `abs_kd`, `measured`, `cells`, both cuts and their pass calls, the
`gate:abs` / `gate:rel` / `soft:rel:p1` weights, the disagreement flags, and the per-target `pds` under
`cd4svd5` / `gate:abs` / `gate:rel` / `soft:rel:p1.0` / `relw`) · `work/plolo/` (byte copies of 262's
fold caches + `SHA256SUMS.txt`).

**Write audit**: `find experiments data scripts src docs submissions -newer MANDATE.md` returns only
`experiments/308-kd-weighted-donors/` plus the two sibling agents' trees (306, 307) and the
coordinator's `docs/` — nothing written by this agent outside its own directory.

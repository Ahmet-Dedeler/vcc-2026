# 209 — the Orion↔Orion excess IS the long-gene axis, and deflating it LOSES to its own permutation control on the direction channel

**Date:** 2026-09-03 · deterministic (seed `20260903` everywhere; `deflate175.py`'s own seed `0`
where 175's row is reproduced) · **no `vcc` command was run, no `cell-eval2` was run, no submission
was built, nothing in `submissions/` was touched, and nothing outside
`experiments/209-orion-subspace/` was written or modified.** Peak RSS **9.53 GiB**
(`subspace209.py`, the 7,310² eigendecomposition stage) against the 12 GB cap; every scoring
process peaked at **4.67–4.94 GiB**. `debench.Surface.pds()` was never called (188 Phase 3).

---

## Verdict

**NO-GO, and the negative is the strong kind: on the transfer cosine — the channel AGENTS.md §8.2
says is 110% of the remaining gap — the maximiser over all 36 arms is a PERMUTATION CONTROL, and
every real deflation arm is negative with a 95% interval excluding zero.**

Five findings, in the order they were measured.

1. **The excess is real, and it is rank 2.** By a criterion declared before any endpoint was seen
   (family-wise: the 95th percentile of the largest eigenvalue over 200 target-permutations), the
   subspace the two Orion donors share beyond what each shares with K562 has dimension **r = 2** on
   the 3-donor-shared columns and **r = 2** on the wider Orion column set. It carries **9.9%
   (HCT116) and 18.6% (HEK293T)** of each donor's production-frame energy on the estimation rows,
   and **7.4% / 14.1%** on the 259 panel targets the estimate never saw, against **1.0% / 0.7%** for
   K562. So 204's hand-off — *"the Orion↔Orion excess is a subspace, not the blob"* — is confirmed
   in its literal form: it is not rank 1.

2. **But it is not distinct from the blob either. It IS the long-gene axis, at rank 2.**
   **82.6% of HCT116's own top PC and 83.3% of HEK293T's own top PC lie inside the two-dimensional
   excess subspace** (leading principal cosine **0.947** against the rank-2 long-gene basis).
   Deflating the long-gene axis first collapses the excess's leading eigenvalue to **28.1%** of its
   original value, and the surviving directions score **worse than the identity on every channel**:
   Δcos **−0.000451** (−2.06%, CI excludes 0), and on the h1-test sign pool **−0.0022**, against the
   long-gene axis's own **+0.0060**. The mandate's Phase-1 falsifier fires.

3. **Production's whitener has already taken it.** The subspace's share of the equal-weight
   consensus's energy goes **6.8% → 0.09%** (C) and **10.7% → 0.12%** (WIDE) through the shipped
   operator; the rank-2 long-gene basis goes **12.5% → 0.14%**, reproducing AGENTS.md's
   9.8% → 0.11% for the rank-1 axis on an independently constructed basis. A random 2-D subspace
   sits at 0.017% → 0.014%. **There is ~0.1% of delivered energy left to take, and taking it costs
   more than it is worth.**

4. **Every real deflation arm loses cosine, monotonically in strength and in dimension, and the
   permutation control wins.** The three random 2-D controls in the same column span read Δcos
   **+0.000065 / +0.000065 / +0.000051** (+0.23% to +0.30% relative); the real rank-2 arm reads
   **−0.001195 (−5.46%)**, and the sweep to r = 20 reaches **−0.002008 (−9.17%)**. **The best real
   arm on the primary endpoint is `subC_l0.25` at −1.20%, still below every permutation control.**

5. **The mechanism claim is falsified with its own instrument.** The mandate's argument was that
   de-correlating two donors that are secretly one platform should move our effective donor count
   toward 3. By 159's `1/cos²(k) = 66.5 + 1920/k` on 197's a15 anchor, `k_eff` goes **3.09 → 2.73**
   (rank-2, full strength), **→ 2.50** at r = 20, and **→ 2.94** for the long-gene axis alone. The
   only arms that move `k_eff` up are the permutation controls (**3.09 → 3.11**), by nothing.

**Best priced arm: `blob2_l0.5` at +0.00086 overall** — and it is not this experiment's arm, it is
175's rank-1 deflation generalised to rank 2, whose family 178 closed. It sits **2.2× below** rule
13's +0.0019 shape-only gate, and its `pds` gain is not resolvable (95% CI `[−0.00022, +0.00657]`).
The mandate's own arm, `subC` at full strength, prices at **−0.00603**. Nothing is shippable and no
slot is requested.

**And a second closure falls out.** 175's per-donor rank-1 deflation — handed to 178 as a fourth
switch, and closed there on composition grounds rather than its own merits — is priced here on the
board-validated instrument for the first time: **Δcos −2.34% to −2.66%, overall −0.00030 (Orion
only) / −0.00052 (all three donors)**. It buys **+0.0073** sign accuracy on the h1-test pool and
pays for it in direction. That is the same `pds`/cosine-versus-sign trade AGENTS.md §4.2 records for
whitening, and 175 §"Reading" predicted it in words; it is now measured in overall units.

---

## 0. Known-answer rows — nine, all pass (`phase0_209.py`, `analyse209.py`)

| # | what | expected | measured | source |
|---|---|---|---|---|
| 1 | 204's one-factor test, Orion↔Orion excess, **five**-line block | z = +5.14 | **z = +5.1353**, resid **+0.030643** | bit-identical to 204's `factor2_204.json` |
| 2 | same, **six**-line block | z (see correction below) | **z = +6.8133**, resid **+0.057776** | bit-identical to 204's stored artefact |
| 3 | **positive control**: `k562_gwps ↔ k562_essential`, same cell line | z = +11.2 | **z = +11.1687**, resid **+0.049514** | 204 §3.3 |
| 4 | 175: cos(top PC HCT116, top PC HEK293T) | 0.872 | **0.87165** | 175 addendum / harness warning |
| 5 | 175: HEK293T top-1 variance share | 19.7% | **19.653%** | same |
| 6 | 175: consensus sign accuracy, no deflation | 0.5697 | **0.5697** | same |
| 7 | 175: consensus sign accuracy, per-donor rank-1 | 0.5770 | **0.5770** (k ≥ 2 removes signal: 0.5683, 0.5650, 0.5567) | same |
| 8 | PLOLO a15 `pds` / a23 `pds`, `k562_gwps` fold, equal weights | 0.79578409 / 0.78912131 | **0.79578409 / 0.78912131** | 162, 204 §0 |
| 8b | — its `q` / `q_w` | 0.50164311 / 0.50487994 | **0.50164311 / 0.50487994** | 204 §0 |
| 9 | the a15 → a23 mass-fix price, 4 folds, pre-committed unweighted mean | +0.004014 scaled | **+0.0040145** (+0.01%) | 162 §2 |

Plus two internal identities: the `equal` arm reproduces 162's cached `pre` to **4.1e-08 / 5.5e-08 /
5.8e-08 / 5.7e-08** across the four folds (float32 cache round-trip, matching 204's 4.1e-08), and
`subCx_r2` — the rank-2 arm re-scored off an independently rebuilt 20-direction basis — reproduces
`subC_l1` **exactly on all four folds** (principal cosines of the two bases: **1.000000, 1.000000**).

> **A correction to 204 §3.3 / §8, for the record.** 204's RESULTS.md reports the Orion↔Orion excess
> as "z = +5.14 (5-line) / **+7.56** (6-line)" and the mandate inherits it as "z = +5.1 to +7.6".
> **204's own stored `work/factor2_204.json` contains no z of +7.56 for that pair in any block.** The
> six-line (`replogle_block`) value is **+6.813**, which this experiment reproduces to sixteen
> digits along with the residual +0.05778 that 204 quotes correctly. The range should read
> **+5.1 to +6.8**. Nothing downstream depends on it — the excess is significant either way, and the
> positive control at +11.17 is unaffected — but the number is wrong where it is written down.

### Rule 11 — a path and a sha256 for every operator

| path | sha256 (8) | stored shape | used by |
|---|---|---|---|
| `data/donors/whitener_union_pruned.npz` | **`e2f65ecf`** | `evals[17942]`, `evecs[17942,17942]`; its own `.json` records `basis: [29394, 17942]`, `basis_prune: essential` | `kolf2` fold (production) + the whitener-survival test |
| `experiments/062-cd4-post-whiten/cloud-output/whiteners/held_k562_gwps.npz` | `13d988fc` | `[17806,17806]` | `k562_gwps` fold |
| `.../held_orion_hct116.npz` | `7cecbf7f` | `[17942,17942]` | `orion_hct116` fold |
| `.../held_orion_hek293t.npz` | `97e486e0` | `[17942,17942]` | `orion_hek293t` fold |
| *present on disk, **NOT** loaded:* `whitener_union.npz` `5608e807`, `whitener_union_full.npz` `965d3cef`, `whitener_full.npz` `dc6bd2ce` | | | listed so rule 11's near-name confusion is closed by hash |

**The mandate's `29,394 × 17,942` is the BASIS the operator was fitted on**, recorded in
`whitener_union_pruned.json`; the stored array is `17,942²`. 193's correction is confirmed against
the file, and both numbers are above so the next agent does not re-resolve it.

Donor inputs (sha256, first 16): `k562_gwps/delta.npy` `0c9422280fed7afb`, `b_ctrl.npy`
`d06cf51d51cde1ae`; `orion_hct116/delta.npy` `72c7a01dc2869bd3`, `b_ctrl.npy` `b3044d439cfac932`;
`orion_hek293t/delta.npy` `d95202d984d9322d`, `b_ctrl.npy` `e04234cffbe10273`. All three match
194's and 204's recorded hashes.

### The fold record (162's, unchanged)

| fold (held) | sources | operator | panel targets | equal-weight `pds` (a23) |
|---|---|---|---:|---:|
| `k562_gwps` | hct116, hek293t | `held_k562_gwps` `13d988fc` | 269/300 | 0.789121 |
| `orion_hct116` | k562, hek293t | `held_orion_hct116` `7cecbf7f` | 290/300 | 0.791823 |
| `orion_hek293t` | k562, hct116 | `held_orion_hek293t` `97e486e0` | 294/300 | 0.777375 |
| **`kolf2`** (leak-free) | all three | `whitener_union_pruned` `e2f65ecf` | 273/300 | 0.610128 |

All four reproduce 204 §4's absolutes exactly.

---

## 1. Phase 1 — what the excess subspace is

### 1.1 The estimator, declared before the endpoint (`subspace209.py` docstring)

* **Frame:** `E062.prepare_main_source`'s output — the object `build_baseline` averages.
* **Rows:** the **2,000 NON-PANEL** targets shared by all three donors (204's seeded draw). The 259
  panel targets are excluded, so no fold's scored rows enter the estimate.
* **Columns:** `C` = the **7,310** genes all three donors measure (the K562 regression is only
  defined where K562 can see); a `WIDE` variant on the **17,770** columns either Orion donor
  measures, where the K562 residual off `C` is the raw row.
* **Residual:** `ε_i(t) = d_i(t) − β_i(t)·d_k(t)`, `β = ⟨d_i,d_k⟩/‖d_k‖²`, **per target**. K562 is
  the instrument for the shared factor; no global attenuation constant is needed, which matters
  because 204 proved none is identified.
* **Statistic:** `S = Sym(ε_h^T ε_e / n)` — a **cross**-covariance, so independent per-donor noise
  enters only at `O(1/√n)` and only shared residual structure accumulates.
* **Dimension:** `r = #{i : λ_i > λ*}` with `λ*` the **95th percentile of the largest eigenvalue of
  S over 200 permutations of the target correspondence** — a Westfall-Young family-wise criterion,
  a property of the estimate alone, and it never touches a score.

Computed exactly through the row space (`rank(ε) ≤ n ≪ p`), so the 7,310² and 17,770² eigenproblems
are 4,000² ones and the permutation null needs no re-factorisation.

### 1.2 The answer: r = 2, and it is big

| | `C` (7,310 cols) | `WIDE` (17,770 cols) |
|---|---:|---:|
| **admitted dimension `r`** | **2** | **2** |
| threshold `λ*` (95th pct of the permuted max) | 0.37156 | 0.57217 |
| permuted max eigenvalue: median / max over 200 draws | 0.2427 / 0.4074 | 0.3870 / 0.6416 |
| eigenvalues 1–6 | **1.4964, 0.4243**, 0.2572, 0.1332, 0.0880, 0.0822 | **2.1401, 0.5920**, 0.3862, 0.2057, 0.1383, 0.1298 |
| energy share, `k562_gwps` (estimation / panel) | 1.66% / **0.98%** | 1.18% / **0.70%** |
| energy share, `orion_hct116` (estimation / panel) | 9.94% / **7.43%** | 8.99% / **7.07%** |
| energy share, `orion_hek293t` (estimation / panel) | 18.62% / **14.11%** | 17.72% / **14.01%** |

The gap between the two Orion donors and K562 — **7–14× on the panel rows the estimate never saw** —
is the excess, and it is not an artefact of the fit: the panel is a clean holdout for it.

**Dimension is not knife-edge but it is not generous either.** The declared family-wise criterion
admits 2; a laxer criterion (above the *median* permuted top eigenvalue) admits 3. The eigenvalue
spectrum is close to symmetric about zero beyond the top few (`Σλ⁺/Σ|λ| = 0.544`), which is what a
noise spectrum looks like. **Dimension is swept as an arm anyway (r = 1, 2, 5, 10, 20) and the
endpoint is monotone in it — in the wrong direction (§2.2).**

### 1.3 The blob overlap — the mandate's mandatory check, and it fires

The long-gene / nuclear-fraction axis was rebuilt independently in the production frame on the same
estimation rows: HCT116 top-1 **8.00%** of variance, HEK293T top-1 **19.99%**, and
**cos(top PC HCT116, top PC HEK293T) = 0.8387** against 175's donor-axis 0.872 — the same object,
measured in a different frame. Its top loadings are the documented gene class:
HEK293T `AUTS2, RORA, GRIP1, EFNA5, EXT1, MECOM, PRKCE, …, NRXN3, WWOX, GPC6, RAD51B, …, PRKG1`;
HCT116 `RGS6, ASTN2, MCC, MAML3, AHNAK2, IMMP2L, FHIT, MAGI1, …` — 0.3–1.7 Mb intron-rich genes on
both.

| overlap of the excess subspace with … | `C` | `WIDE` |
|---|---:|---:|
| `‖P_S · pc1(HCT116)‖²` | 0.557 | **0.826** |
| `‖P_S · pc1(HEK293T)‖²` | 0.511 | **0.833** |
| principal cosines vs the rank-2 long-gene basis | 0.762, 0.183 | **0.947**, 0.258 |
| share of the subspace's own energy the blob explains | 0.307 | **0.482** |
| `‖P_S · (target-independent mean)‖²`, HCT116 / HEK293T / K562 | 0.314 / 0.494 / 0.005 | 0.437 / 0.729 / 0.004 |

**On the columns the deflation actually has to act on (`WIDE`), 83% of each Orion donor's dominant
direction lies inside the excess subspace.** 204's §3.2 finding that blob deflation moves *its own*
numbers by <1% is correct and is not in conflict: 204 was deflating the **target-independent mean
vector**, which carries only 0.04–0.22% of table energy here. The **long-gene axis** is a different
object — a *between-perturbation* direction, not a mean — and it is 83% of the excess.

### 1.4 Production's whitener has already removed it (`subspace209.py`, kolf2 fold, `e2f65ecf`)

Energy share of the equal-weight three-donor consensus, before and after the shipped operator at
`gamma 0.65, lam 0`:

| subspace | raw consensus | **whitened consensus** | ratio |
|---|---:|---:|---:|
| excess, `C`, r = 2 | 6.76% | **0.092%** | 74× |
| excess, `WIDE`, r = 2 | 10.71% | **0.124%** | 86× |
| the rank-2 long-gene basis | 12.45% | **0.137%** | 91× |
| a **random** 2-D subspace, same column span | 0.017% | 0.014% | 1.2× |

The rank-2 long-gene row reproduces AGENTS.md's **9.8% raw → 0.11% whitened** for the rank-1 axis on
a basis built here from scratch. **What is left after whitening is 0.09–0.14% of delivered energy,
about 9× a random 2-D subspace's share on a 0.1% base.** That is the whole prize a pre-whitening
deflation can claim, and §2 measures what claiming it costs.

### 1.5 Blob-first re-estimation — is anything left that is not the blob?

Deflating the rank-2 long-gene basis from both Orion donors and re-running the entire estimator:

| | before | after the blob |
|---|---:|---:|
| admitted `r` | 2 | 3 |
| top eigenvalue | 1.4964 | **0.4209 (28.1% of before)** |
| eigenvalues 2–4 | 0.4243, 0.2572, 0.1332 | 0.3070, 0.1398, 0.0935 |
| principal cosines against the original subspace | — | 0.976, 0.352 |

**Something survives** — a 3-D residual excess, still above its own permutation threshold, and it
overlaps the original subspace only partly. **It is that surviving object that §2 and §4 price, and
it is negative on every channel.**

---

## 2. Phase 2 — the deflation sweep, cosine primary

162's instrument, byte-unchanged: four folds, challenge frames A/B/C, shipped geometry
(`alpha 0.982, gamma 0.65, lam 0, cd4_f 0.10`), a23 render for `pds`/`cos`, a15 render for `q`/`q_w`.
The only thing that varies across arms is the projection applied to the two Orion donors' prepared
effects before `build_baseline`. Intervals are a **4,000-draw paired target-cluster bootstrap**
resampling targets within each fold (194 §5.2's protocol, seed `20260903`), on the per-target rank
statistic and on the per-target cosine separately.

Baseline 4-fold macro cosine: **0.021885**. 204 §8 measured this instrument at **±0.20% relative
≈ ±0.00006 overall**.

### 2.1 The strength sweeps

**`subC` — the excess subspace, `C`, r = 2 (the mandate's arm):**

| | equal | ×0.25 | ×0.5 | ×0.75 | ×1.0 |
|---|---:|---:|---:|---:|---:|
| **Δcos** | 0 | **−0.000263** | **−0.000552** | **−0.000862** | **−0.001195** |
| Δcos, relative | 0 | −1.20% | −2.52% | −3.94% | **−5.46%** |
| Δ`pds`, 4-fold | 0 | +0.00092 | −0.00025 | −0.00452 | −0.01204 |
| Δ`pds`, **`kolf2`** (leak-free) | 0 | −0.00040 | −0.00076 | −0.00073 | −0.00107 |
| Δ`q_w` | 0 | −0.00035 | −0.00062 | −0.00092 | −0.00118 |
| **`k_eff`** | **3.09** | 3.01 | 2.92 | 2.83 | **2.73** |

**`subW` — the same on the wider Orion column set:** Δcos 0 → −0.000305 → −0.000630 → −0.000971 →
**−0.001318 (−6.02%)**; `k_eff` 3.09 → **2.70**. Its 4-fold Δ`pds` collapses to **−0.01900** while
its `kolf2` Δ`pds` rises to **+0.00352** — the sharpest fold disagreement in the screen, and §4.3
explains it.

**`blob2` — the long-gene axis alone (the control):** Δcos 0 → −0.000114 → −0.000238 → −0.000369 →
−0.000501 (−2.29%); Δ`pds` 4-fold 0 → +0.00218 → **+0.00314** → +0.00207 → −0.00108; `k_eff` → 2.94.

### 2.2 Dimension is monotone, in the wrong direction

| r | 1 | **2** | 5 | 10 | 20 |
|---|---:|---:|---:|---:|---:|
| Δcos | −0.001166 | **−0.001195** | −0.001615 | −0.001823 | **−0.002008** |
| relative | −5.33% | −5.46% | −7.38% | −8.33% | **−9.17%** |
| Δ`pds`, 4-fold | −0.01281 | −0.01204 | −0.01588 | −0.01691 | −0.01705 |
| `k_eff` | 2.74 | 2.73 | 2.61 | 2.55 | **2.50** |

**Adding directions makes it worse.** The r = 2 admitted by the pre-declared criterion is close to
the least-bad rank, and even r = 1 is −5.33%.

### 2.3 The two maximisers, kept separate (062's mistake, not repeated)

| endpoint | maximiser over all 36 arms | value | reading |
|---|---|---:|---|
| **Δcos (PRIMARY)** | **`permW_s1` — a permutation control** | **+0.000065 (+0.30%)** | doing nothing is better than every real arm |
| Δ`pds` | `blob2_l0.5` — the long-gene axis at half strength | +0.00314, CI `[−0.00022, +0.00657]`, p = 0.966 | not resolvable, and it is 175's arm, not this one |
| best **real** arm on Δcos | `subC_l0.25` | −0.000263 (−1.20%), CI `[−0.000386, −0.000137]` | negative, interval excludes 0 |

The `pds` and cosine optima disagree by more than the whole span of the screen — **`pds` peaks at a
blob deflation of half strength while cosine peaks at zero deflation.** 197 §3 measured why: `pds`
is a rank statistic and our whitening already converts direction into retrieval at 1.5× the cohort
median, so `pds` is the member least sensitive to what an arm does to direction. **Had this arm been
selected on `pds`, as 062's CD4 dose was, it would have shipped a −5% cosine change as a gain.**

### 2.4 `k_eff` — the mechanism claim, falsified

By 159 §5.5's `1/cos²(k) = 66.5 + 1920/k` on 197 §4.1's a15 board↔surface anchor (×2.078):

| arm | implied board cosine | **`k_eff`** |
|---|---:|---:|
| incumbent (`rc_hybrid`, measured) | 0.079249 | **3.09** |
| permutation controls (3 seeds, WIDE) | 0.079484 / 0.079483 / 0.079433 | **3.11 / 3.11 / 3.11** |
| `blob2_l0.5` | 0.078386 | 3.02 |
| **`subC_l1.0`** (rank 2, full) | **0.074923** | **2.73** |
| `subW_l1.0` | 0.074477 | 2.70 |
| `subCx_r20` | 0.071979 | **2.50** |

The mandate's hypothesis was that if the Orion pair is partly one line our true `k` is nearer 2, so
removing the shared excess should move `k_eff` toward 3. **It moves it toward 2.5.** Removing the
shared excess removes shared *signal* faster than it removes shared *artefact* — which is exactly
what §1.3 predicts, since 83% of the direction being removed is the axis the whitener has already
handled and the remainder is real biology the two Orion lines genuinely share.

---

## 3. Phase 3 — six members and the gate (rule 12)

**Arm class: profile-changing, shape-only.** `alpha`, the emission law, the spread policy,
`keep_top` and the per-cell library sizes are all untouched; only the pre-whitening effect matrix
moves. Rule 13's **−0.0146 `reach` toll does not apply** — it is an EMISSION-LAW toll (audit 181;
169 measured +0.0006 for profile arms; the `rc_hybrid` board point fixed the class boundary). **Gate
+0.0019**, and **+0.0089** for the one arm that changes amplitude (§3.2).

| member | instrument | `subC_l1.0` (the arm) | `blob2_l0.5` (best priced) | verdict |
|---|---|---|---|---|
| `cos`→`mse` | 197's identity through 096's `claim` route; `d raw_mse/d cos = −2.0826` ⇒ **+0.3653 overall per unit cosine** (§8.2's +0.00365 per +0.01, re-derived) | Δcos **−0.001195** (−5.46%) → **−0.00153** overall | Δcos −0.000238 (−1.09%) → **−0.00031** | **negative on every real arm** |
| `pds` | PLOLO, 4 folds, a23, paired target bootstrap; `Δ/0.452/6×1.011` | **−0.01204** `[−0.02020, −0.00343]`, p = 0.003 → **−0.00449** | +0.00314 `[−0.00022, +0.00657]`, p = 0.966 → **+0.00117** | arm negative and resolvable; control positive and not |
| `nmae` | 187 §5's exact identity, `alpha` fixed: `d raw_nmae = −2·W_hi·d q_w`, `W_hi = 0.0153` | Δ`q_w` −0.00118 → **−0.0000100** overall | Δ`q_w` −0.00012 → −0.0000010 | **pinned.** Largest \|Δ`q_w`\| over all 34 Orion-only arms is 1.52e−3 ⇒ \|Δoverall\| ≤ **1.3e−5** |
| `reach` | ordering half only: Δ`q` (a15 render) and the h1-test sign pool (§4.1). Not locally priceable (rule 13: 2/6 signs; rule 9: `n_real`) | Δ`q` **−0.00290**; h1 pool **+0.0058** | Δ`q` +0.00005; h1 pool **+0.0060** | unpriced, and the two instruments disagree in sign — see below |
| `fid` | ordering half: same Δ`q`. Yield half unreadable (rule 9) | as above | as above | unpriced |
| `jac` | the call set is `keep_top` on the delivered profile. `keep_top` is untouched, but a profile-class arm *does* move which genes clear the bar, so this is **not** exempt by 192's mechanism | — | — | **unmeasured**, see §3.3 |

**Total, `alpha` re-solved (§3.2):** `subC_l1.0` **−0.00603**, `blob2_l0.5` **+0.00086**,
best real arm `subC_l0.25` **+0.00000**.

### 3.1 The `reach`/`fid` ordering half: two instruments, two answers

The h1-test sign pool (175's instrument, `real_sig & adj`, n = 39,945) says the deflation **helps**
ordering; PLOLO's `q` on the delivered challenge-frame profile says it **hurts** (−0.0029). They are
not measuring the same object: the pool statistic is read off the *raw donor consensus* before the
operator, the CD4 residual and the render, and 172 measured +0.052 raw `reach` available from
ordering alone on predictions whose profile is pinned. Rule 13 says `reach` is not locally priceable
at all (2/6 signs, worse than chance), so **neither number is priced**, and the disagreement is
recorded rather than resolved. It does not change the verdict: on the two members that *are*
priceable the arm is negative by 3.2× the gate.

### 3.2 Amplitude, and what re-solving `VCC_ALPHA` costs (`price209.py`)

Deflation with the donor-level norm restoration holds the **delivered** median row norm to within
**1.2%** (`subC_l1.0`: s = 0.9883) and **1.5%** (`blob2_l0.5`: s = 0.9854) — inside 204 §6's ±1.5%
band, so both stay in the shape-only class. The one arm without norm restoration (`subC_l1_raw`)
moves it **−2.7%** and is charged the +0.0089 amplitude gate; it prices at −0.00592 and the
distinction is moot.

`VCC_ALPHA` is degree-(−1) homogeneous in the delivered amplitude, so the re-solve is a division:
**0.982 → 0.9937** for `subC_l1.0`, **0.982 → 0.9966** for `blob2_l0.5`. 187's global-α optimum is
a property of the shipped *shape*, so re-solving costs nothing to find and restores `A` exactly.

**If `VCC_ALPHA` is NOT re-solved**, a15 sits on the `B` leg of `min(claim, B)` in two of three
contexts (ρ = 0.992 / 0.940 / 0.881), so shrinking the prediction's across-perturbation spread
shrinks the credited jackknife term. Charged explicitly:

| arm | s | α re-solve | Δ`mse` overall, α re-solved | Δ`mse` overall, α fixed | extra cost of not re-solving |
|---|---:|---:|---:|---:|---:|
| `subC_l1.0` | 0.9883 | 0.9937 | −0.00153 | −0.00180 | **−0.00027** |
| `blob2_l0.5` | 0.9854 | 0.9966 | −0.00031 | −0.00066 | **−0.00036** |
| `subC_l1_raw` | 0.9731 | 1.0092 | −0.00152 | −0.00214 | −0.00063 |
| permutation controls | 1.0000 | 0.9820 | +0.00008 | +0.00008 | −0.0000001 |

**Saturation.** `claim` is a functional of the emitted counts' jackknife noise; this arm class
touches neither the emission law, the spread policy, `keep_top`, nor the per-cell library sizes, so
**`claim` does not move.** We ship it at 100.0% of `PRED_TRACE_CAP_K`, so there is no headroom
above; what moves is `B`, and it moves **down** by `s²` — i.e. the arm can only reduce the credited
term. That is the −0.00027 to −0.00036 charged above, and re-solving α removes it.

### 3.3 What is NOT measured, and why (rule 12's narrow exemption, honestly)

**`nmae`/`reach`/`fid`/`jac` were not measured on `h1-test` with a significance-selected pool, and
this is a structural limitation, not an omission.** Two reasons, in order of force:

1. **The surface cannot see the arm.** `h1-test`'s legal donors are `k562_gwps` + `orion_hek293t`
   only (AGENTS.md §6.4 trap 1 — `h1_2025` is pure oracle leakage there). **It holds out HCT116,
   one of the two lines whose shared excess is being removed.** On that surface the arm degenerates
   from "de-correlate two donors" to "remove a fixed 2-D projection from one donor", which is a
   different arm with a different mechanism. This is measured, not asserted: on the h1-surface
   donor pair the rank-2 deflation reads **+0.0057** sign accuracy against **+0.0058** on the full
   three-donor consensus, i.e. essentially all of the pool-level effect is the single-donor
   projection and none of it is the de-correlation. **Reporting a `reach` number from that surface
   would be reporting an unrelated arm's number.** 200 caught exactly this failure mode on the CD4
   residual today.
2. Producing the four DE members requires rebuilding the h1-test artifact from deflated donor
   tables for each arm × seeds — a count build this mandate's memory cap and the concurrent
   production build do not allow, and the numbers would then be unreadable anyway under rule 9
   (`n_real` is 222 on h1-test; the challenge is in no such regime).

**`jac` is unmeasured and is NOT exempt.** 204 pinned it by "`keep_top` is not touched"; that
mechanism holds for call-set and spread arms, not for a profile arm, because moving the delivered
profile moves which genes clear the Wilcoxon/BH bar. It is stated as unmeasured. It does not change
the verdict: `jac`'s scaled span is 0.34–0.39 and we are already *ahead* of the leader on it
(197 §4.4), so a member we lead cannot rescue an arm that is −3.2 gates down on the two members
that are priceable.

---

## 4. Phase 4 — the four falsifiers

### 4.1 Permutation control — **it beats every real arm** (mandate falsifier 1)

Random orthonormal 2-D subspaces drawn inside the *same column span*, three seeds each, at full and
half strength, deflated from the same two donors with the same norm restoration:

| arm | Δcos | 95% CI | relative | Δ`pds` (4-fold) | overall |
|---|---:|---|---:|---:|---:|
| **`permW_s1` (random, WIDE)** | **+0.000065** | [+0.000045, +0.000086] | **+0.30%** | −0.00008 | **+0.00005** |
| `permW_s2` | +0.000065 | [+0.000047, +0.000084] | +0.30% | +0.00002 | +0.00009 |
| `permW_s3` | +0.000051 | [+0.000033, +0.000070] | +0.23% | −0.00011 | +0.00002 |
| `permC_s1/2/3` (random, C) | +0.000000 / +0.000004 / −0.000005 | — | ±0.02% | +0.00006 / −0.00000 / +0.00005 | +0.00002 / +0.00000 / +0.00001 |
| **`subC_l0.25`** (best real arm) | **−0.000263** | [−0.000386, −0.000137] | **−1.20%** | +0.00092 | +0.00000 |
| **`subC_l1.0`** (the arm) | **−0.001195** | [−0.001682, −0.000691] | **−5.46%** | −0.01204 | **−0.00603** |

**204's repression falsifier *tied* its best real arm; this one *beats* it.** The three WIDE controls
are all positive with intervals excluding zero, and the real arm is 18× further from them in the
other direction. Two readings, and both matter:

* The +0.30% from a random 2-D deflation is **statistically resolvable and economically nil** —
  +0.00005 overall, 38× below the gate. It is the norm restoration rescaling the donors by a hair,
  not a lever, and it should not be mistaken for one.
* **A deflation that removes nothing in particular does not hurt. A deflation that removes the
  Orion-shared excess hurts a great deal.** That is the cleanest statement of the negative: the
  excess is not noise the consensus would be better off without.

### 4.2 Blob control — the gain does not survive it, on either channel (mandate falsifier 2)

Two arms answer this, and they agree:

| arm | Δcos | Δ`pds` (4-fold) | h1 sign pool | overall |
|---|---:|---:|---:|---:|
| `blob2` — the long-gene axis alone, rank 2 | −0.000501 | −0.00108 | **+0.0060** | −0.00104 |
| `blob2_l0.5` — the same, half strength | −0.000238 | **+0.00314** | — | **+0.00086** |
| **`subCblobDefl`** — ONLY the blob-orthogonal excess (r = 3, re-estimated after the blob) | **−0.000451** | +0.00011 | **−0.0022** | **−0.00054** |
| `blob2_then_subC` — the blob and then the residual excess | −0.000940 | −0.00118 | — | −0.00165 |

**Every positive number in this screen belongs to the long-gene axis. The part of the excess that is
NOT the long-gene axis is negative on the cosine (−2.06%, interval excluding zero), negative on the
h1 sign pool (−0.0022 against the axis's +0.0060), and null on `pds`.** Composing the two makes both
worse. The mandate's falsifier 2 is unambiguous: **the arm does not survive deflating the long-gene
axis first, and the arm is therefore the blob.**

And the blob is production's whitener's business, not a donor-side arm's: §1.4 measures 91× of it
removed already, and 178 recorded the same conclusion from the other direction ("deflation acts only
via the whitener refit and the sign source").

### 4.3 All four folds, including `kolf2`, and no absolute from the Orion pair (mandate falsifier 3)

Δ`pds` per fold. **No absolute is quoted from `orion_hct116` or `orion_hek293t` anywhere in this
document** — they appear only as differences against a matched `equal` arm on the same row set.

| arm | `k562_gwps` | `orion_hct116` | `orion_hek293t` | **`kolf2` (leak-free)** | 4-fold mean |
|---|---:|---:|---:|---:|---:|
| `subC_l0.25` | −0.00072 | +0.00176 | +0.00304 | **−0.00040** | +0.00092 |
| `subC_l1.0` | −0.01370 | −0.01813 | −0.01525 | **−0.00107** | −0.01204 |
| `subW_l1.0` | +0.00029 | −0.04184 | −0.03797 | **+0.00352** | −0.01900 |
| `blob2_l0.5` | +0.00068 | +0.00298 | +0.00797 | **+0.00092** | +0.00314 |
| `permW_s1` | +0.00010 | +0.00000 | −0.00056 | **+0.00012** | −0.00008 |

**Two structural biases in the fold set, stated so they are not read as findings.**

* On the two Orion folds the **truth carries the artefact being removed** (the axis is 19.7% of
  HEK293T's own between-perturbation variance), so deflating it from the *source* removes a
  component the truth genuinely has. That biases those folds against the arm on cosine — and
  indeed the two Orion folds hold the largest cosine losses.
* Those same two folds hold the largest `pds` *gains* for the blob arms, which is the hubness
  mechanism of AGENTS.md §4.2 in miniature: a direction shared by the truth and by every decoy
  raises similarity to everything, so removing it from the prediction buys discrimination while
  costing bulk cosine.
* **`kolf2` is the only fold with the shipped operator, the shipped three-source mixture and no
  leakage, and on it the mandate's arm is negative at every strength** (−0.00040 → −0.00107). The
  4-fold mean is the pre-committed predictor and is what §2 reports; `kolf2` is reported alongside
  every arm so the decomposition is visible.

### 4.4 The K562 sanity check — asserted, and it holds (mandate falsifier 4)

Every Orion-only arm carries a hard assertion in `score209.py`: the K562 source's prepared effects,
after the arm is applied, must be **bit-identical** to the un-deflated build. Measured
`k562_max_abs_dev` = **0.0 exactly** on **all four folds × all 34 Orion-only arms** (136 rows), including the
`k562_gwps_pm50` PerturbMatch-blended source. The assertion is live and it fired once, correctly, on
the one arm that deflates all three donors by design (175's `rank1_175`, dev 1.10 / 1.02 / 1.03) —
which is why that arm is reported separately and is exempt by construction.

Orthonormality of every deflation basis: max `|VVᵀ − I|` ≤ **6.7e−15**.

---

## 5. Verdict in overall units

| arm | what it is | overall (α re-solved) | overall (α fixed) | vs the gate |
|---|---|---:|---:|---|
| **`permW_s1/2/3`** | **the permutation controls** | **+0.00005 / +0.00009 / +0.00002** | same | 21–95× below +0.0019 |
| `blob2_l0.5` | 175's arm at rank 2, half strength — **the best priced arm, and not this experiment's** | **+0.00086** | +0.00050 | **2.2× below +0.0019** |
| `blob2_l0.25` | same, quarter strength | +0.00066 | +0.00051 | 2.9× below |
| `subC_l0.25` | **the mandate's arm, weakest strength** | **+0.00000** | −0.00008 | **nil** |
| `subW_l0.25` | wide variant, weakest strength | +0.00023 | +0.00017 | 8× below |
| `rank1_orion` / `rank1_175` | **175's rank-1 deflation, priced for the first time** | **−0.00030 / −0.00052** | −0.00086 / −0.00112 | **negative** |
| `subCblobDefl` | **the blob-orthogonal excess alone** | **−0.00054** | −0.00057 | **negative** |
| `blob2` | the long-gene axis, full strength | −0.00104 | −0.00160 | negative |
| **`subC_l1.0`** | **the mandate's arm, full strength** | **−0.00603** | −0.00630 | **3.2× the gate, negative** |
| `subW_l1.0` | wide variant, full strength | −0.00878 | −0.00917 | negative |
| `subCx_r20` | dimension pushed to 20 | **−0.00894** | −0.00923 | negative |

**Both of the mandate's stated full results are delivered, and they point the same way:**

> **"The excess is the blob at higher rank"** — measured: 83% of each Orion donor's dominant
> direction lies inside the rank-2 excess; deflating the axis first collapses the excess's leading
> eigenvalue to 28%; and production's whitener already removes 91× of it (12.5% → 0.14% of consensus
> energy). The residual that is genuinely blob-orthogonal is **negative on the cosine (−2.06%),
> negative on the h1 sign pool (−0.0022), and worth −0.00054 overall.**
>
> **"Deflation ties its permutation control"** — measured, and it is worse than a tie: **the
> permutation control is the maximiser of the primary endpoint over all 36 arms**, and the real arm
> is −5.46% relative cosine against the control's +0.30%.

**No slot is requested and none should be spent on this.**

---

## 6. What this closes, and what it hands on

**Closed.**

1. **The donor–donor correlation structure is not a lever.** 204 closed scalar re-weighting (equal
   weights are the maximum of the transfer cosine over 29 alternatives). This closes the subspace
   generalisation it handed on: a projection does no better than a coefficient, and both lose to
   doing nothing. The Orion↔Orion excess is real, significant (z = +5.1 to +6.8), large (7–14% of
   panel energy), rank-2 — **and 83% of it is a sequencing artefact production already removes.**
   AGENTS.md §4.2's equal-donor-weight finding now survives a *shape* attack as well as a *weight*
   attack.
2. **175's rank-1 deflation is closed on its own merits.** 178 retired it on composition grounds
   (*"gains compose sub-additively, costs do not"*). Priced here on the board-validated instrument it
   is **−0.00030 to −0.00052 overall**, with Δcos −2.34%/−2.66% and intervals excluding zero. Its
   +0.0073 sign gain is real and is the ordering channel; it costs more direction than the ordering
   is worth, which is §8.1 rule 1 exactly.
3. **"Our true `k` is nearer 2" is falsified as an operational claim.** Removing the shared platform
   component moves `k_eff` **3.09 → 2.73** (→2.50 at r = 20), not toward 3. The Orion pair is not
   "secretly one line" in any sense the consensus can exploit: what they share beyond K562 is 83%
   artefact the whitener has taken and 17% biology worth keeping.

**Handed on — three things.**

1. **A correction to 204 §3.3 and §8, and to anything that inherited it.** The six-line
   Orion↔Orion z is **+6.81**, not +7.56; 204's own `factor2_204.json` contains no +7.56 for that
   pair in any block. Reproduced bit-identically here. The residual (+0.0578) and the five-line z
   (+5.135) that 204 also quotes are correct.
2. **The whitener-survival measurement is a cheap general test, and it should be run before any
   future donor-side shape arm.** `energy share in the raw consensus → in the whitened consensus`,
   against a random subspace of the same dimension, costs one operator application and tells you
   whether production has already taken what you are about to remove. Here it read
   **12.5% → 0.14% against a random 0.017% → 0.014%**, i.e. 91× removed, and it predicted the
   negative before any score was computed. **`subspace209.py`'s last stage does it in ~40 s.**
3. **The `pds`-versus-cosine optimum split is now measured on a fifth arm family, and it is the
   widest yet.** `pds` peaks at a *half-strength blob deflation* (+0.00314) while cosine peaks at
   *no deflation at all*; the two prescriptions differ by the entire span of the screen. 204 §8
   recommended pricing a direction arm on Δcos first; this experiment is the case where taking
   that advice and ignoring it give opposite verdicts. **Any future donor-side arm must report
   both maximisers separately or it will ship a cosine loss as a `pds` gain.**

**Not pursued, per the mandate's territory clause.** The whitening operator was held frozen on every
arm and never refitted — refitting it on deflated donor tables is a second variable and a different
arm (178 identified it as the channel through which deflation could still act). 193 owns the
operator; 194 owns the donor tables; 207 owns the CD4 residual; 208 owns emission riders. None were
touched.

---

## 7. Files

```
phase0_209.py    Phase 0: 204's one-factor test + positive control, 175's rank-1 rows, and
                 rule 11's operator/donor sha table.  Aborts on any failure.
subspace209.py   Phase 1: the excess-subspace estimator with its pre-declared dimension
                 criterion, the permutation null, the energy shares, the long-gene axis and
                 its loadings, the blob overlaps, the whitener-survival test, and the
                 blob-first re-estimation.
ext209.py        the 20-direction basis for the dimension sweep; its top-2 reproduces
                 subspace-C to principal cosines 1.000000
score209.py      Phase 2: PLOLO with a per-donor projection; the K562 bit-identity assertion
sign209.py       the ordering channel: h1-test sign-pool accuracy, both the three-donor
                 consensus and the two donors h1-test actually carries
analyse209.py    the bootstrap, the six-member pricing, the sweeps, the two maximisers
price209.py      Phase 3: the amplitude term, the `VCC_ALPHA` re-solve, and saturation
work/            scores209.json (144 fold x arm rows: 36 arms x 4 folds), subspace209.json, sign209.json,
                 analyse209.json, price209.json, phase0_209.json, subspace-*.npz,
                 blob-axis.npz, logs, plolo/ (162's fold caches, copied read-only)
```

```bash
cd experiments/209-orion-subspace
uv run python phase0_209.py                       # known-answer rows first; aborts on failure
uv run python subspace209.py && uv run python ext209.py
for f in k562_gwps orion_hct116 orion_hek293t kolf2; do uv run python score209.py --folds $f; done
uv run python sign209.py && uv run python analyse209.py && uv run python price209.py
```

**Memory.** Peak RSS **9.53 GiB**, in `subspace209.py`'s eigen/whitener stage (the 2.57 GiB
production operator plus three 2,259 × 18,533 float64 donor frames). The scoring processes peak at
**4.67–4.94 GiB** with 204's ordering (stream the donors first, allocate the operator last) and one
process per fold.

**Hygiene.** Nothing outside this directory was written or modified. `data/donors/` and every other
experiment's `work/` were read only — 204's cached `prepare_main_source` outputs and 162's fold
caches were **copied into this experiment's own `work/`** and `prep204.WORK` / `plolo.WORK` were
repointed before any call, so no cache miss could write into another agent's directory. No `vcc`, no
`cell-eval2`, no build, no `submissions/`. `debench.Surface.pds()` was never called. Every arm is
deterministic: `subCx_r2` re-derived from an independently rebuilt basis reproduces `subC_l1`'s
`pds` to all printed digits on all four folds.

# 194 — re-estimating the donor tables: the plain pseudobulk difference wins, and the two mechanisms that make it win are measurable

**Date:** 2026-09-03 · deterministic (seed `20260903` everywhere) · **no `cell-eval2` was run, no
`vcc` command was run, no submission was built, nothing in `submissions/` was touched, and no
existing file in `data/donors/` was modified.** Everything written is inside
`experiments/194-donor-reestimation/` except three new tables in `data/donors/194-debias-*/`,
which the mandate authorises. Peak RSS **13.0 GiB** against a 14 GiB cap (three runs were killed
by the guard at 13.3/14.1/17.8 GiB and re-engineered; §7).

---

## Verdict

**NO-GO on every alternative estimator tested. The 2026-08-22 plain pseudobulk difference is the
best of eleven candidates on cross-line transfer, and this closes the estimator as a lever.**

Not because the alternatives do nothing — two of them demonstrably do exactly what the
literature says they should — but because of a two-part mechanism that this experiment measured
directly:

1. **Cell-level knockdown weighting is real on K562 and it is a MAGNITUDE correction.** On the 420
   K562-gwps targets whose own gene is measured well enough to read, down-weighting escaper cells
   moves the effect vector **along its own direction** by `+0.0275 ± 0.0039` (7.1σ), amplifies it
   by 2.5–8.3%, and the amplification **correlates +0.258 with the target's own measured residual
   expression** — worse knockdown, bigger correction, exactly as the escaper-dilution model
   predicts. Both matched falsifiers are flat (`kd_perm` `−0.0006 ± 0.0007`, correlation `−0.027`).
   **But `pds_cosine` is exactly scale-invariant per row (AGENTS.md §4.1) and
   `prepare_main_source` rescales each bank by its own median centred norm, so the entire
   recoverable part of the correction lands in the one channel production is already invariant to.**
   What is left over is `1 − cos = 1.2%` of the row, and it is noise: it does not move `pds`.
2. **On Orion the same weighting does not identify escapers at all, and it inflates the
   target-independent blob.** The correlation between amplification and measured residual flips
   sign (HCT116 `−0.232`, HEK293T `−0.159`), and the share of table energy in the single
   target-independent mean vector goes **1.00% → 6.67% (`kd_post`) → 10.45% (`kd_gate80`)** on
   HEK293T while the permuted control stays at 0.94%. That artefact is why the HCT116↔HEK293T
   cross-line cosine appears to rise 63% (`+0.0366 → +0.0596`) while the cross-platform
   K562↔Orion cosines do not move at all (`+0.0276 → +0.0282`, `+0.0271 → +0.0256`) — the two
   Orion lines are agreeing better about a shared library artefact, not about biology. Growing
   that vector is the harness pathology of AGENTS.md §6.4, applied to ourselves.

**κ² does not move.** Recomputed through 159's own `donorpairs.pair`, τ = κ² is `0.0951` on the
incumbent estimator and `0.0815–0.0978` across all eight variants — the knockdown arms *lower* it,
and the best (`qc`, `0.0978`) is inside the noise of a 250-target anchor. **No estimator recovers
more of the line-independent component than the plain difference does.**

**Two incumbent defects were found and priced. Both are real and both are negligible.** (a) The
Orion builder's control regex `^neg` matches the gene **NEGR1**, so 154 HCT116 and 233 HEK293T
*perturbed* cells sit in the control pool and NEGR1 has no row in either shipped table; correcting
it moves the tables by 0.09% and `pds` by `−0.00000`. (b) `log1p` is concave and the control pool
carries **424–1,086×** more counts than the median target, so the incumbent's rows carry a
systematic, target-independent negative offset on low-expressed genes; correcting it exactly
(full tables, shipped row sets, `data/donors/194-debias-*/`) gives `Δpds = +0.00033`, 95% CI
`[−0.00088, +0.00155]`, 2 of 4 folds positive. NO-GO.

**Raising `min_cells` is strictly harmful**, monotonically: `20 → 50 → 100 → 200` gives
`Δpds 0.00000 → −0.00113 → −0.01678 → −0.08317`. A noisy row from a donor beats falling back to
the other donors.

**Guide-level aggregation and every robust-over-guides estimator (candidates 2 and 3) are
structurally unavailable here, which had never been checked.** In all three donors **92% of
targets carry exactly one guide construct** (K562-gwps 92.01%, both Orion lines 91.84%; mean
1.08 constructs per target, max 5). There is no guide axis to average, trim, or MODZ over. Relatedly,
Orion's `pass_guide_filter` — the only guide QC in either builder — **passes 100.0% of
7,943,468 cells** in both lines, i.e. the incumbent's one nod to guide QC filters nothing.

---

## 0. Known-answer rows (rule 10) — four, all reproduced

### 0.1 The incumbent rebuilt, bit-identically (mandate Phase 0)

`scripts/build_donor_table.py --donor k562_gwps` re-run into
`work/repro/k562_gwps`, against the shipped `data/donors/k562_gwps`:

| artefact | shipped sha256 | rebuilt sha256 | |
|---|---|---|---|
| `delta.npy` | `0c942228…e296aa9` | `0c942228…e296aa9` | **bit-identical** |
| `b_ctrl.npy` | `d06cf51d…c5d81b957` | `d06cf51d…c5d81b957` | **bit-identical** |
| `perts.txt`, `meta.json` | — | — | **byte-identical (`diff` empty)** |

Source hashes recorded in `work/sha256-sources.txt`:
`ReplogleWeissman2022_K562_gwps.h5ad` `d3269d9c86…8278f104`,
`orion/metadata/gene_metadata.parquet` `33900a0dba…1f2bf0fdbf`,
`orion_hct116/delta.npy` `72c7a01dc2…2afde4ec2b33`,
`orion_hek293t/delta.npy` `d95202d984…5d064d7b1`.

**The build is deterministic.** Every number below is therefore a measurement of an estimator, not
of run-to-run variance.

### 0.2 The re-implementation reproduces the incumbent (`accumulate.py`, scheme `v0`)

The streaming accumulator used for all eleven variants must reproduce the shipped tables on the
rows it shares, or nothing it produces can be read:

| donor | rows compared | `b_ctrl` | max abs Δ on `delta` | as fraction of median row norm |
|---|---:|---|---:|---:|
| `k562_gwps` | 1,569 | **bit-identical** | 4.77e-07 | **6.8e-08** |
| `orion_hct116` | 790 | **bit-identical** | 4.77e-07 | **6.0e-08** |

The residual is float32-vs-float64 accumulation order and is 4 orders of magnitude inside the
1e-6 tolerance the mandate allows.

### 0.3 The scoring path reproduces 162's arms (`score194.py --variants shipped`)

Scoring the incumbent tables through the PLOLO fold `k562_gwps`, challenge frames A/B/C, against
162's own recorded endpoints:

| arm | member | 162 | 194 | Δ |
|---|---|---:|---:|---:|
| a15 | `pds` | 0.79578409 | 0.79578409 | **0.00e+00** |
| a15 | `sgn` (`q`) | 0.50164311 | 0.50164311 | **0.00e+00** |
| a15 | `sgn_w` (`q_w`) | 0.50487994 | 0.50487994 | −5.2e-12 |
| a23 | `pds` | 0.78912131 | 0.78912131 | **0.00e+00** |
| a23 | `cos` / `rr` / `sse` | — | — | ≤ 5.4e-11 |
| a23 | `sgn` | 0.46878748 | 0.46777958 | **−1.0e-03** |

**A finding that changes how the sign channel must be read.** `pds`, `cos`, `rr` and `sse`
reproduce to float noise under both render laws, but `sgn`/`sgn_w` reproduce **exactly** under
a15 and drift by `1e-3` under a23. a23's exact-mass repair (`mass_consistent_effect`) leaves a
large fraction of delivered genes within ~1e-11 of zero, so a BLAS-level perturbation flips
thousands of signs. **`q` and `q_w` are therefore reported through the a15 (no-repair) render in
every table below; `pds` is reported through the shipped a23 render.** Reading sign accuracy
through a23 would put a 1e-3 numerical floor under a channel whose real arm-to-arm differences
are 1e-4.

### 0.4 159's κ² route reproduces (`kappa194.py --shipped`)

Using 159's own `donorpairs.pair`, imported not reimplemented, on the shipped tables:

| pair | 159 §4 `resid same` | 194 |
|---|---:|---:|
| `k562_essential ↔ k562_gwps` (same line, two libraries) | +0.2857 | **+0.2857** |
| `k562_gwps ↔ orion_hct116` | +0.0321 | +0.0311 |
| `k562_gwps ↔ orion_hek293t` | +0.0286 | +0.0269 |
| `orion_hct116 ↔ orion_hek293t` | +0.0300 | +0.0286 |

The anchor is exact; the three 3,000-target rows differ only through which 3,000 of >9,000 shared
targets the shared RNG drew in a different pair order. **All variant-vs-variant comparisons in §5
are exact by construction** — those tables have fewer than 3,000 rows, so no subsampling occurs.

> **A correction for the record.** 159 §3.3 quotes route 2 as `0.0897 / 0.2627 = 0.341`. Neither
> `0.0897` nor `0.2627` appears in 159's own §4 pair matrix, whose corresponding entries are
> `+0.0321` and `+0.2857`, giving **0.112**. This experiment reproduces the §4 table, not the
> §3.3 headline. τ from this route is therefore ~0.11, not 0.34 — which moves 159's
> `τ = 0.30 ± 0.10` consensus down if route 2 is reweighted, and is worth a look by whoever next
> touches the ceiling estimate. It does not change anything here: **the comparison across
> estimators uses one protocol throughout.**

---

## 1. What the incumbent estimator actually is, and what the data allows

`scripts/build_donor_table.py` / `build_orion_donor.py`: sum a target's counts, sum the control's,
normalise each to `TS = 5e4`, `log1p`, subtract. `min_cells = 20`. No knockdown weighting, no
guide QC on K562-gwps, `pass_guide_filter` on Orion.

Measured facts about the substrate, all new here (`prepass.py`, `accumulate.py --stage obs`):

| | K562-gwps | Orion HCT116 | Orion HEK293T |
|---|---:|---:|---:|
| cells | 1,989,578 | 3,409,169 | 4,534,299 |
| control cells | 75,328 | 165,931 | 219,071 |
| cells/target (median) | 184 | 151 | 202 |
| batches / gem groups | 274 | 109 | 223 |
| **targets with exactly 1 guide construct** | **92.01%** | **91.84%** | **91.84%** |
| guides per target (mean / max) | 1.08 / 5 | 1.08 / 5 | 1.08 / 5 |
| `pass_guide_filter` pass rate | n/a (no filter) | **1.0000** | **1.0000** |
| cells whose own target gene is measured | 73.9% | 95.13% | 95.17% |
| **expected counts of the target's own gene, median** | **0.376** | **0.106** | **0.200** |
| cells with expected counts ≥ 0.5 | 40.2% | 27.3% | 35.5% |
| targets with a usable per-cell readout (median expected ≥ 0.5) | 511 / 1,569 | 399 / 790 | 505 / 794 |
| control median depth vs perturbed | 1.000 | 17,312 / 17,291 | matched |
| **target-level residual expression (median)** | **0.155** | **0.300** | **0.510** |
| targets achieving >80% knockdown | **47.6%** | 21.4% | **4.4%** |

Three consequences that constrain the whole experiment and were not known before it:

- **Candidates 2 and 3 of the mandate (guide-level aggregation; robust location over guides —
  trimmed mean, median-of-guides, within-donor MODZ) cannot be built.** There is no guide axis:
  92% of targets have one construct. This is a property of both screens, not a choice.
- **A per-cell knockdown readout carries at most ~1 expected count.** At `exp ≈ 1`, the posterior
  `P(knocked down | obs)` runs 0.90 / 0.64 / 0.26 for `obs = 0 / 1 / 2` — real discrimination, but
  one weak bit per cell.
- **The three donors' repression differs 3.3×.** HEK293T leaves 51% of the target transcript
  standing and only 4.4% of its targets clear the challenge's own >80% criterion. That motivated a
  candidate the mandate did not list — per-target de-attenuation (`kdrescale`, §5).

Also found: the control regex `^neg` in `build_orion_donor.py` matches **NEGR1**. Its 154 (HCT116)
/ 233 (HEK293T) perturbed cells are pooled into the control and **NEGR1 has no row in either
shipped Orion table**. Priced in §5 (`ctrlfix`): 0.09% of the table, `Δpds −0.000004`.

---

## 2. The eleven estimators, and why they are all one accumulation pass

Every standard alternative in this literature is the incumbent's own functional with a per-cell
weight in front of the sum, so K variants are K weighted sums over one stream (`schemes.py`):

| id | mandate candidate | weight |
|---|---|---|
| `v0` | incumbent | `1` |
| `cellmean` | 5 | `1/depth` (mean of relative profiles, not of counts) |
| `qc` | 5 | `1[depth ∈ p2..p99.5, mito ≤ p97.5, genes ≥ p2]` |
| `batchmatch` | 4 | `f_ctrl(batch)/f_T(batch)`, capped 20× |
| `kd_soft` | 1 | `clip(1 − (obs+1)/(exp+1), 0, 1)` |
| `kd_post` | 1 | `P(knocked down | obs)`, two-component Poisson, prior 0.8, residual 0.2 |
| `kd_gate50` / `kd_gate80` | 1 | `1[posterior ≥ 0.5 / 0.8]` — the challenge's own >80% criterion |
| `kd_perm` | **falsifier** | `kd_post` permuted within target |
| `kd_perm_depth` | **falsifier** | `kd_post` permuted within (target × depth decile) |
| `combo` | 1+4+5 | `kd_post · qc · batchmatch` |
| `debias` / `debias2` | new | post-hoc second-order `log1p` bias correction, 1× / 2× |
| `ctrlfix` | new | control pool with NEGR1 removed |
| `kdrescale` (+`_perm`) | new | per-target de-attenuation by `1/(1−r_target)`, floor 0.25 |

Doses `scheme@λ` are exact, not approximate: because every scheme is normalised to mean weight 1
per target, the sums under `w = (1−λ) + λ·w_scheme` are exactly `(1−λ)·S_v0 + λ·S_scheme`. `λ > 1`
is **not expressible** as a cell weighting (weights would go negative) and is refused.

**A design bug the depth audit caught, recorded because the same trap is easy to fall into again.**
The first version gated knockdown applicability *per cell* on `exp = p_ctrl · depth`. That leaves
the shallow cells of an applicable target un-weighted and the deep ones down-weighted, i.e. it
builds a library-depth preference into the estimator — the confound that broke
`data/harness/orion-hct116.h5ad`. It showed up as `kd_perm` shifting mean weighted depth to
**0.947×**, as much as `kd_soft` itself. Applicability is now decided **per target** on the
target's median expected count, and every scheme's depth shift is clean:

| | `v0` | `kd_post` | `kd_gate80` | `kd_perm` | `kd_perm_depth` | `cellmean` |
|---|---:|---:|---:|---:|---:|---:|
| K562 mean weighted depth / baseline | 1.0000 | 0.9992 | 0.9912 | **1.0002** | 0.9994 | 0.7761 |
| HCT116 | 1.0000 | 0.9920 | 0.9797 | **1.0002** | 0.9928 | 0.7960 |
| HEK293T | 1.0000 | 0.9874 | 0.9778 | **1.0000** | 0.9884 | 0.7809 |

`cellmean`'s 0.78 is intended (that is what `1/depth` weighting means) and is reported, not fixed.

---

## 3. Does each estimator do what it claims? (mechanism, before transfer — `diag194.py`)

The decisive test for knockdown weighting needs no transfer score. Escaper cells carry a control
profile, so removing them must make the effect **bigger along its own direction**:
`proj = ⟨d_w − d_0, d_0⟩ / ⟨d_0, d_0⟩` must be positive, and the permuted control must be flat.
A second, independent prediction: worse-repressed targets have more escaper mass, so `proj` must
**correlate positively with the target's own measured residual expression**.

Statistics over the rows each variant actually changes (33–60% of targets — a median over all rows
is a median over rows the arm never touched):

### K562-gwps (420 changed rows) — **the mechanism is confirmed**

| variant | ampl. | rel. change | **proj** | ±sem | cos to `v0` | **corr(proj, residual)** |
|---|---:|---:|---:|---:|---:|---:|
| `kd_post` | 1.0253 | 0.162 | **+0.0275** | 0.0039 | 0.9877 | **+0.258** |
| `kd_gate80` | 1.0830 | 0.347 | **+0.0467** | 0.0048 | 0.9479 | **+0.321** |
| `kd_soft` | 1.0869 | 0.361 | +0.0380 | 0.0037 | 0.9449 | +0.300 |
| `kd_post@0.5` | 1.0091 | 0.081 | +0.0135 | 0.0019 | 0.9968 | +0.255 |
| **`kd_perm`** | 1.0100 | 0.138 | **−0.0006** | 0.0007 | 0.9906 | **−0.027** |
| **`kd_perm_depth`** | 1.0134 | 0.158 | **+0.0031** | 0.0009 | 0.9878 | +0.151 |
| `cellmean` | 1.1178 | 0.496 | +0.0056 | 0.0008 | 0.8967 | +0.091 |
| `qc` | 1.0136 | 0.174 | −0.0025 | 0.0003 | 0.9852 | −0.014 |
| `batchmatch` | 1.0713 | 0.391 | −0.0020 | 0.0003 | 0.9309 | +0.008 |
| `debias` | 0.9993 | 0.031 | −0.0011 | 0.0001 | 0.9995 | +0.094 |

`kd_post` is 7.1σ from zero and 40× its own permuted control; the dose is exactly linear
(`@0.5` gives half of `@1.0`); and the correlation with measured repression has the predicted
sign in the real arm and not in the permuted one. **Escaper dilution is real and this estimator
finds it.** `kd_perm_depth` at `+0.0031` says ~11% of the amplification is depth-mediated and 89%
is knockdown-specific.

### Orion (254 / 311 changed rows) — **the mechanism inverts**

| | HCT116 `proj` | ±sem | corr | HEK293T `proj` | ±sem | corr |
|---|---:|---:|---:|---:|---:|---:|
| `kd_post` | **−0.0046** | 0.0057 | **−0.232** | **−0.0103** | 0.0122 | **−0.159** |
| `kd_gate80` | −0.0019 | 0.0084 | −0.201 | −0.0110 | 0.0187 | −0.176 |
| `kd_perm` | +0.0002 | 0.0019 | +0.130 | +0.0011 | 0.0022 | +0.008 |
| `cellmean` | −0.0510 | 0.0027 | −0.047 | −0.0727 | 0.0031 | +0.007 |
| `batchmatch` | −0.0120 | 0.0018 | +0.003 | −0.0427 | 0.0087 | −0.067 |

No amplification along the effect's own direction, and the correlation with measured repression is
**negative** — the opposite of the prediction. Meanwhile the rows move a great deal
(`rel change` 0.32–0.95, `cos to v0` down to 0.75). **On Orion the readout is selecting Poisson
noise, not escapers**, and the most likely reason is in the table in §1: HEK293T's median target
retains 51% of its transcript, so there is no bimodal escaper population to separate — repression
is partial in all cells, and selecting cells with few target counts selects shallow sampling.

### And what the Orion arms are actually doing (`work/geometry.json`)

Share of table energy in the single target-independent mean vector — the metric's own zero, and
the vector AGENTS.md §6.4 warns about:

| donor | `v0` | `kd_post` | `kd_gate80` | `kd_perm` | `batchmatch` | `cellmean` | `debias` |
|---|---:|---:|---:|---:|---:|---:|---:|
| `orion_hek293t` | **1.00%** | **6.67%** | **10.45%** | 0.94% | **16.88%** | 0.66% | 0.84% |
| `orion_hct116` | 0.52% | 0.99% | 1.53% | 0.52% | 0.62% | 0.57% | **0.41%** |
| `k562_gwps` | 0.81% | 0.86% | 0.90% | 0.80% | 0.76% | 0.86% | 0.88% |

`kd_post` multiplies HEK293T's blob by **6.7×** and `batchmatch` by **16.9×**, while the permuted
control leaves it exactly where it was and K562 is untouched by any of them. That is the whole
explanation of §5's `kd`/`batchmatch` losses and of §6's illusory Orion-pair gain.

---

## 4. Effect-norm distributions and retention (mandate reporting requirement)

| donor | `v0` effect norm p10 / median / p90 | cells retained | cells reweighted |
|---|---|---|---|
| `k562_gwps` | 5.01 / **6.84** / 11.75 | `qc` 95.00%, `kd_gate80` 83.8% of applicable | `kd_post` 33.2% of perturbed cells (511 targets) |
| `orion_hct116` | 5.29 / **7.79** / 13.42 | `qc` 94.91%, `kd_gate80` 77.1% | `kd_post` 48.0% (399 targets) |
| `orion_hek293t` | 4.66 / **7.01** / 11.87 | `qc` 94.68%, `kd_gate80` 64.1% | `kd_post` 59.7% (505 targets) |

Median effect norm of the changed rows, relative to `v0`: `kd_post` ×1.025 / ×1.052 / ×1.108,
`kd_gate80` ×1.083 / ×1.180 / ×1.384, `cellmean` ×1.118 / ×1.071 / ×1.056, `debias` ×0.9993 /
×1.0004 / ×0.9999, `kdrescale` ×1.25 / ×1.59 / ×2.17 (median norm 6.84→8.57, 7.79→12.35,
7.01→15.18 — by construction, since it divides by the achieved repression).

---

## 5. Transfer scores — PLOLO, four folds, challenge frame (`score194.py`)

162's instrument, unchanged: hold each donor library out as truth over the real 300-target panel,
rebuild the effect matrix from the remaining libraries with a whitening operator fitted **without**
the held line, apply the shipped geometry (`alpha 0.982, gamma 0.65, lam 0, cd4_f 0.10`), render
through production's `profile_from_effect` into challenge frames A/B/C, read back as the scorer
does. **`debench.Surface.pds()` was never used** (188 Phase 3). The scored line is held out of
everything including the whitening basis; `h1_2025` and `h1-test` were never touched.

**The truth is held at the incumbent table and only the SOURCES vary.** This is not a detail — see
§5.3, where varying both moves the same arm by 0.006 `pds`.

### 5.1 The screen, all three donors swapped, Δ vs a matched `v0` on the same row set

| variant | k562 fold | hct116 fold | hek293t fold | kolf2 fold | **mean Δpds** | fold sd | Δcos | Δ`q` | Δ`q_w` | Δ`sse` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `v0` | — | — | — | — | **0.00000** | — | — | — | — | — |
| `kd_post@0.25` | +0.00107 | +0.00064 | +0.00082 | +0.00091 | **+0.00086** | 0.00018 | +0.00009 | +0.00001 | +0.00002 | −5e-06 |
| `kd_post@0.4` | +0.00155 | −0.00015 | +0.00030 | +0.00198 | +0.00092 | 0.00100 | +0.00009 | +0.00021 | +0.00003 | −6e-06 |
| `kd_gate80@0.25` | +0.00165 | +0.00135 | −0.00042 | +0.00136 | +0.00098 | 0.00094 | +0.00008 | +0.00021 | +0.00016 | −4e-06 |
| **`kd_perm@0.25`** (falsifier) | +0.00026 | +0.00019 | +0.00075 | +0.00097 | **+0.00054** | 0.00038 | −0.00002 | −0.00030 | −0.00005 | +1e-06 |
| `kd_post@0.1` | +0.00041 | +0.00039 | +0.00035 | +0.00035 | +0.00038 | 0.00003 | +0.00005 | −0.00004 | +0.00001 | −3e-06 |
| `kdrescale` | +0.00339 | −0.00115 | −0.00239 | +0.00287 | +0.00068 | 0.00288 | −0.00009 | +0.00048 | +0.00028 | +6e-05 |
| `kdrescale@0.25` | +0.00119 | −0.00030 | −0.00053 | +0.00063 | +0.00025 | 0.00080 | −0.00001 | +0.00003 | −0.00010 | +7e-06 |
| `debias` | +0.00280 | −0.00110 | −0.00143 | +0.00057 | +0.00021 | 0.00194 | −0.00002 | −0.00016 | +0.00036 | +4e-06 |
| `kd_perm` (falsifier) | −0.00140 | −0.00147 | +0.00021 | +0.00340 | +0.00018 | 0.00228 | −0.00045 | −0.00041 | +0.00024 | +2e-05 |
| `debias2` | +0.00629 | −0.00214 | −0.00468 | +0.00103 | +0.00013 | 0.00473 | −0.00005 | −0.00037 | +0.00065 | +7e-06 |
| `ctrlfix` | −0.00001 | +0.00000 | +0.00000 | +0.00000 | **−0.00000** | 0.00001 | −0.00000 | −0.00000 | +0.00000 | −0e+00 |
| `kd_perm_depth@0.25` | −0.00096 | +0.00059 | −0.00079 | +0.00054 | −0.00015 | 0.00083 | −0.00004 | +0.00010 | −0.00005 | +2e-06 |
| `qc@0.5` | +0.00026 | −0.00286 | −0.00332 | +0.00145 | −0.00111 | 0.00234 | −0.00005 | −0.00022 | +0.00001 | +4e-06 |
| `kd_gate50` | +0.00141 | −0.00778 | −0.00511 | +0.00646 | −0.00126 | 0.00643 | −0.00040 | +0.00025 | +0.00005 | +3e-05 |
| `kd_post` | +0.00100 | −0.00615 | −0.00604 | +0.00496 | −0.00156 | 0.00549 | −0.00020 | +0.00044 | +0.00017 | +1e-05 |
| `kd_perm_depth` | −0.00700 | −0.00060 | −0.00652 | +0.00261 | −0.00288 | 0.00468 | −0.00063 | −0.00009 | −0.00031 | +4e-05 |
| `qc` | −0.00128 | −0.00786 | −0.01008 | +0.00351 | −0.00393 | 0.00621 | −0.00035 | −0.00018 | +0.00005 | +2e-05 |
| `cellmean` | −0.01911 | −0.00217 | +0.00331 | −0.00162 | −0.00490 | 0.00979 | −0.00174 | +0.00068 | −0.00016 | +9e-05 |
| `kd_soft` | +0.00147 | −0.01270 | −0.01677 | +0.00388 | −0.00603 | 0.01024 | −0.00083 | −0.00004 | +0.00005 | +5e-05 |
| `kd_gate80` | −0.00448 | −0.01109 | −0.01407 | +0.00423 | −0.00635 | 0.00812 | −0.00111 | +0.00073 | +0.00013 | +7e-05 |
| `batchmatch` | −0.00992 | −0.00617 | −0.00866 | −0.00598 | −0.00768 | 0.00193 | −0.00176 | −0.00035 | −0.00036 | +7e-05 |
| `combo` | −0.00945 | −0.01706 | −0.02501 | −0.00213 | **−0.01341** | 0.00984 | −0.00216 | −0.00047 | −0.00048 | +9e-05 |

Absolute `v0` baseline: `pds` 0.78977 / 0.79214 / 0.77770 / 0.60982 (269 / 290 / 294 / 273 panel
targets), `q` 0.50170 / 0.48755 / 0.46847 / 0.47222, `sse` 0.99918 / 0.99944 / 0.99960 / 0.99994.

### 5.2 Falsifier 1 — the permuted-weight control is NOT beaten

PLOLO's own LOO RMS is **0.001494 scaled `pds`** (162), i.e. ≈0.00067 raw. Paired
target-cluster bootstrap, 4,000 draws, resampling targets within each fold:

| arm | Δpds | 95% CI | p(Δ>0) |
|---|---:|---|---:|
| `kd_post@0.1` | +0.00038 | [−0.00012, +0.00089] | 0.930 |
| `kd_post@0.25` | +0.00086 | [−0.00040, +0.00210] | 0.907 |
| `kd_post@0.4` | +0.00092 | [−0.00108, +0.00290] | 0.817 |
| **`kd_perm@0.25`** | **+0.00054** | [−0.00011, +0.00127] | **0.947** |
| `kd_perm_depth@0.25` | −0.00015 | [−0.00092, +0.00062] | 0.350 |
| `kd_gate80@0.25` | +0.00098 | [−0.00090, +0.00283] | 0.847 |
| `kdrescale@0.25` | +0.00025 | [−0.00032, +0.00082] | 0.808 |

**Not one interval excludes zero, and the permuted falsifier has a HIGHER `p(Δ>0)` (0.947) than
the real arm it controls (0.907).** The knockdown-specific increment at the best dose is
`+0.00086 − 0.00054 = +0.00032` raw `pds`, half the instrument's own LOO RMS. The small positive
trend at small `λ` is a generic mild-reweighting effect, not knockdown: it survives permutation.

At the shipped mapping (`predicted overall = PLOLO Δpds / 6 × 1.011`, and raw→scaled `pds`
≈ ÷0.45), the best arm is worth **+0.0003 overall** against rule 13's shape-only ship floor of
**+0.0019** — 6× below it. **Nothing here is shippable, and nothing here is close.**

### 5.3 Falsifier 2 — why the truth had to be held fixed

Rebuilding the TRUTH with the same estimator as the sources (secondary arm, 3 folds):

| variant | Δpds, truth = incumbent | Δpds, truth = variant | swing |
|---|---:|---:|---:|
| `kd_post` | **−0.00156** | **+0.00433** | **+0.0059** |
| `kdrescale` | +0.00068 | −0.00007 | −0.0008 |
| `debias` | +0.00021 | −0.00038 | −0.0006 |
| `cellmean` | −0.00490 | −0.01102 | −0.0061 |

`kd_post` flips from a loss to a gain of 0.0043 purely from letting both sides of the comparison
be built by the same estimator. That is 4× the arm's own effect and 6× PLOLO's noise floor. **An
estimator screen that rebuilds its own truth measures agreement about a shared artefact.** Any
future donor-estimator work in this repo must hold the truth fixed.

### 5.4 Falsifier 3 — does it survive whitening? (022's question)

Every number in §5.1–5.3 is measured **through the production operator**, never before it: 062's
held-out whitener for the three donor folds, the shipped `whitener_union_pruned` for the leak-free
`kolf2` fold, followed by the CD4 residual, the renderer and the scorer's read-back. There is no
pre-whitening arm in this experiment, so 022's failure mode ("it looked good before the operator")
cannot apply. The `kolf2` fold — the one that uses the exact shipped operator and the exact shipped
three-source mixture, and is leak-free because `kolf2` contributed no row to any production
operator — behaves like the others: no arm clears the noise floor.

### 5.5 `min_cells` (candidate 5) — the incumbent's 20 is not just acceptable, it is right

| `min_cells` | k562 | hct116 | hek293t | kolf2 | **mean Δpds** |
|---:|---:|---:|---:|---:|---:|
| **20 (shipped)** | 0.78977 | 0.79214 | 0.77770 | 0.60982 | **0.00000** |
| 50 | 0.77863 | 0.79413 | 0.77379 | 0.61833 | −0.00113 |
| 100 | 0.74746 | 0.78820 | 0.75434 | 0.61233 | **−0.01678** |
| 200 | 0.67673 | 0.72026 | 0.67148 | 0.56829 | **−0.08317** |

Monotone and large. A 20-cell pseudobulk row is worth keeping: the alternative is falling back on
the other donors, or on the mean response, and that is much worse than a noisy measurement.

---

## 6. κ² and the effective donor count (`kappa194.py`)

159's route 2, one protocol, all estimators. `k562_essential` (same cell line, different library,
covers **0/300** panel targets) is the same-line anchor; the extras pool was made donor-independent
so gwps and essential share 250 rows rather than 38.

| estimator | same-line anchor | k562↔HCT116 | k562↔HEK293T | **τ = κ²** | κ | HCT116↔HEK293T |
|---|---:|---:|---:|---:|---:|---:|
| shipped full tables | +0.2857 | +0.0311 | +0.0269 | 0.1016 | 0.319 | +0.0286 |
| **`v0` (matched rows)** | **+0.2878** | +0.0276 | +0.0271 | **0.0951** | **0.308** | +0.0366 |
| `qc` | +0.2776 | +0.0277 | +0.0265 | 0.0978 | 0.313 | +0.0330 |
| `debias` | +0.2886 | +0.0276 | +0.0273 | 0.0951 | 0.308 | +0.0359 |
| `kd_perm` (falsifier) | +0.2788 | +0.0272 | +0.0254 | 0.0942 | 0.307 | +0.0341 |
| `kdrescale` | +0.2968 | +0.0282 | +0.0277 | 0.0942 | 0.307 | +0.0367 |
| `batchmatch` | +0.2578 | +0.0244 | +0.0239 | 0.0939 | 0.306 | +0.0300 |
| `cellmean` | +0.2775 | +0.0256 | +0.0235 | 0.0885 | 0.298 | +0.0369 |
| `kd_post` | +0.3049 | +0.0282 | +0.0256 | 0.0882 | 0.297 | **+0.0596** |
| `kd_gate80` | +0.3047 | +0.0275 | +0.0222 | 0.0815 | 0.286 | **+0.0717** |

**κ² does not rise for any estimator.** The best (`qc`, 0.0978) is +2.8% relative on a 250-target
anchor; the knockdown arms lower it by 7–14%. The effective donor count is unchanged: production
still stands on three panel-covering lines, and 159's `1/k` law still says ~10 are needed.

**The one row that looks like a win is the trap.** `kd_post`/`kd_gate80` raise the same-line
anchor by 5.9–6.6% — that is escaper removal working, exactly as §3 predicted, on two libraries of
one cell line — and raise the HCT116↔HEK293T cosine by **63%/96%**. But the two cross-platform
pairs move by `+0.0006` and `−0.0015`. Cross-line agreement that appears only between the two
X-Atlas lines and only for the arms that inflate the blob 6.7–10.5× (§3) is the long-gene
technical axis, not shared biology. **An experiment that had measured only the Orion pair would
have reported a large, reproducible, entirely artefactual gain.**

---

## 7. Rule-12 accounting: all six members

A donor-table change moves the profile, so no member is exempt by construction. For the best arm
(`kd_post@0.25`) and the delivered arm (`full:debias`):

| member | instrument | value | verdict |
|---|---|---|---|
| `pds` | PLOLO, 4 folds, challenge frame, bootstrap CI | `+0.00086 [−0.00040, +0.00210]` raw; `+0.00033 [−0.00088, +0.00155]` for `full:debias` | inside the noise floor; ~+0.0003 overall vs a +0.0019 ship floor |
| `mse` | `sse = 1 − cos²θ = min_α mse`, measured per arm | `Δsse ≤ 9e-05` on every arm; transfer cosine 0.0219 → 0.0219 | **pinned.** §4.1b: `mse` is clamped at 0 until cos > 0.105; we are at 0.022 and no arm moves it by more than 0.4% relative |
| `nmae` | `q_w` (187 §5: `nmae` moves only through `q_w`), a15 render | `Δq_w = +0.00002` (`kd_post@0.25`), `+0.00024` (`full:debias`) | unmoved; ≈93% of the `nmae` deficit is direction (187) and the direction did not move |
| `fid` | not locally priceable (rule 13: 2/6 signs; rule 9: `n_real`) | `Δq = +0.00001 / +0.00012` | ordering half unmoved; yield half unmeasurable |
| `reach` | not locally priceable (rule 13: 2/6 signs, worse than chance) | `Δq = +0.00001 / +0.00012` | ordering half unmoved |
| `jac` | not locally priceable (rule 13) | call set is not touched by a donor-table change at fixed `keep_top` | 192 owns the call set |

Since every priceable member is inside its own noise floor and the two unpriceable DE members read
their ordering off a per-gene sign accuracy that moved by `1e-5`, **no arm here can be argued onto
the board.** No slot is requested and none should be spent on this.

---

## 8. Deliverable

`data/donors/194-debias-{k562_gwps,orion_hct116,orion_hek293t}/` — full-row, byte-compatible
tables (same gene axis, same target list, same dtype and layout as the incumbent), built by
`fulltables.py` as an exact closed-form correction of the shipped tables plus the per-target total
counts. `delta.npy` sha256 `ce64b851…`, `b5f57465…`, `210db298…`; each `meta.json` records the
sha256 of the shipped table it derives from.

Table properties (per donor):

| donor | N_control / N_median_target | median relative change | median amplification |
|---|---:|---:|---:|
| `k562_gwps` | **424×** | 0.0316 | 0.99933 |
| `orion_hct116` | **1080×** | 0.0237 | 1.00037 |
| `orion_hek293t` | **1086×** | 0.0175 | 0.99951 |

Transfer (per PLOLO fold, all three source tables swapped at once, shipped truth):

| fold (held line) | shipped `pds` | `194-debias` `pds` | Δ |
|---|---:|---:|---:|
| `k562_gwps` | 0.78912 | 0.79307 | +0.00395 |
| `orion_hct116` | 0.79182 | 0.79086 | −0.00097 |
| `orion_hek293t` | 0.77738 | 0.77520 | −0.00217 |
| `kolf2` | 0.61013 | 0.61065 | +0.00052 |
| **mean** | **0.74211** | **0.74245** | **+0.00033** |

Scored against the **shipped** tables with the shipped row sets (no restricted-row caveat), 4
PLOLO folds: `Δpds +0.00033`, 95% CI `[−0.00088, +0.00155]`, `p(>0) = 0.702`, 2 of 4 folds
positive, `Δsse +0.000001`, `Δq +0.00012`, `Δq_w +0.00024`.

**These tables are correct and they are not an improvement. Do not swap them into production.**
They are shipped so the correction is on disk and reproducible, and because the ratio column is
worth knowing: the control pool outweighs a median target's counts by three orders of magnitude,
which is where the incumbent's small shared offset comes from.

---

## 9. What this closes, and the one thing it opens

**Closed.** The donor estimator is not a lever on the direction channel. Eleven candidates —
including every one the mandate ranked, and three it did not — either fail, or succeed at
something the pipeline is invariant to, or succeed at inflating an artefact. Two of the mandate's
five candidates (guide-level aggregation, robust location over guides) are **structurally
unavailable**: 92% of targets in all three donors carry a single guide
construct (92.01% / 91.84% / 91.84%), so there is no guide axis to aggregate over. `min_cells` is at its optimum. Combined
with 166 (no further panel-covering line exists) and 184 (the line-specific residual is not
predictable), **the direction channel has no remaining cheap lever from the donor side.**

**Open, and it is not an estimator question.** The three donors' *achieved repression* differs
3.3× (median residual expression K562 0.155, HCT116 0.300, HEK293T 0.510), and only **4.4%** of
HEK293T targets clear the challenge's own >80% knockdown criterion against **47.6%** for K562.
That is a large, previously unrecorded quality difference between the three tables production
weights equally. `kdrescale` — the obvious per-target de-attenuation — was tried here and is
inside the noise (`+0.00068 ± 0.00288`, falsifier `−0.00187`), so the naive correction does not
collect it. But the fact itself belongs in the donor-mixing conversation, which is 193's and the
coordinator's territory, not this experiment's: a donor whose guides only half-repress is
contributing a systematically attenuated *and* differently-shaped row, and the equal-weight
consensus (010/019) was established before anyone measured that.

---

## 10. Files

```
prepass.py        metadata-only pass over Orion (6 s/line): guides per target, batch
                  composition, depth/mito distributions, pass_guide_filter rate
schemes.py        the eleven per-cell weight schemes + the two permutation falsifiers,
                  with the per-target applicability gate and the depth audit
accumulate.py     two streaming stages (obs, sums); Arrow record-batch reader; buffered
                  flush; `v0` reproduces the incumbent to 6e-8 of a row norm
tables.py         sums -> donor tables; exact dose blending; log1p debias; ctrlfix;
                  per-target knockdown de-attenuation
fulltables.py     the full-row `debias` tables written into data/donors/194-debias-*
diag194.py        mechanism test: projection on v0, amplification, corr with residual
score194.py       PLOLO scoring; a23 and a15 render laws; per-target pds; min_cells gate
kappa194.py       159's donorpairs.pair on variant tables -> tau = kappa^2
crossline194.py   blob share, PC1 share, PC1 alignment to the incumbent
analyse.py        assembles work/scores*.json into the tables above
work/             logs, sums, variant tables, scores*.json, diag.json, kappa.json,
                  geometry.json, sha256-sources.txt
```

**Memory.** Peak RSS 13.0 GiB against the 14 GiB cap. Three earlier designs were killed by the
guard and re-engineered, and the reasons are worth recording because they are Orion-specific:
each Orion parquet file is a **single row group**, so pyarrow must decompress ~1.9 GiB before
`iter_batches` can slice it and the reader alone peaks at 5.6 GiB; converting the two list columns
to pandas materialises one numpy array per cell and reached **27 GiB**; and holding one
`PreparedSource` per (variant, donor, fold) reaches 8 GiB across a screen. The fixes — Arrow flat
buffers, a 6,000-cell buffered flush, row-by-row accumulation instead of fancy-index temporaries,
a 790-target Orion row set, and evicting the variant side of the preparation cache — are in the
code with comments.

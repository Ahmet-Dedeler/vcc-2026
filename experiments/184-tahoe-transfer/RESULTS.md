# 184 — Tahoe-100M as a cross-line transfer substrate: the operator is real, and it does not exist in our CRISPRi data

**Date:** 2026-09-03 · deterministic (seed `20260903` everywhere) · **no `cell-eval2` was run, no
submission was built, no production code was touched, and nothing outside
`experiments/184-tahoe-transfer/` was written except the downloads in `data/external/tahoe/`
(1.2 GB, budget was 120 GB).** Peak RSS 5.3 GB against a 12 GB cap.

---

## Verdict

**The Phase-2 gate FAILS as written, and the reason it fails is the most informative thing in the
experiment. Do not proceed to Phase 3.**

1. **A basal-conditioned cross-line transfer operator is real, large and easy to fit — on
   Tahoe.** With 49 training lines the recipient's DMSO profile predicts its line-specific
   response residual at cosine **+0.4635** (sem 0.017 across the 50 held lines), against a
   context-free control of **+0.0316**, a line-permutation null of **+0.0389** (sem 0.014, not
   distinguishable from zero) and a leaky oracle ceiling of **+0.5224**. The arm reaches **89% of
   the ceiling** and sits ~19 combined sem above its null.
2. **It grows steeply in the number of training lines** — `+0.2293` at k=3, `+0.4635` at k=49 —
   exactly the signature the brief said would mean "we have found the missing information."
3. **It is entirely an AMPLITUDE law, not a DIRECTION law.** A per-gene gain is sign-preserving
   by construction, so the arm that wins cannot move sign accuracy at all; the only family that
   can (a low-rank rotation) buys **+0.0025 sign / +0.0063 sign on the weak-consensus half** at
   k=49, against a gate of +0.010. Best sign delta anywhere in the run: **+0.0037** (k=20).
4. **And the law does not exist in our CRISPRi donors.** Running the identical measurement on the
   same five lines 161 used: arm **+0.0665**, context-free control **+0.0667**, line-permutation
   null **+0.0663**. Arm ≡ control ≡ null to three decimals. The recipient-specific content is
   **zero**, while the oracle says +0.1126 of diagonal structure is there to be had.
5. **That CRISPRi null is not a noise floor.** Degrading Tahoe to CRISPRi-like reliability
   (ρ 0.281 → 0.047, a 6× degradation) shrinks the Tahoe effect but never removes its
   specificity: at ρ = 0.047, k=5, the arm reads **+0.0441 against a null of +0.0044** — a 10×
   ratio at a reliability *below* our worst donor table, where CRISPRi reads 1.0×.
6. Model-free confirmation, no ridge and no basis: the correlation between the line-specific part
   of the empirical per-gene gain and the line-specific part of log basal expression is
   **+0.666 on Tahoe** (and **+0.640** after the 6× noise degradation — the statistic is
   noise-insensitive because the gain averages over 1,137 drugs) against **+0.166 on CRISPRi**.
   45% of the line-specific gain variance explained versus 3%.

**So: the five CRISPRi failures were not sample-size limited.** 161 measured
`cos(orth, TRUE orth) = 0.0592/0.0606` with five lines and called the sample size the honest
limit of its NO-GO. This experiment reproduces that number from a different model class
(**+0.0665**), shows what the same estimator does when the information *is* present (**+0.4635**),
and shows it still does it at CRISPRi noise levels. The limitation was the data's content, not
its k. That closes the line-specific-residual route on CRISPRi donor tables for good, in the
strong form 161 could not reach.

---

## 0. Known-answer rows (rule 10) — all reproduced

### 0.1 Phase 0, before touching Tahoe (`known.py`)

Reads 156's frozen exact-space truth and `data/donors/orion_*` read-only.

| quantity | 159 | this experiment |
|---|---:|---:|
| harness route, same target, HCT116↔HEK293T | +0.1881 | **+0.1881** (med +0.1967, sd 0.1035) |
| harness route, different target | +0.1790 | **+0.1790** |
| harness route, excess | +0.0091 | **+0.0091** |
| within-HCT116 / within-HEK293T | +0.2476 / +0.4175 | **+0.2477 / +0.4176** |
| donor route, 120 targets, challenge axis, same | +0.0317 | **+0.0317** |
| donor route, different | +0.0006 | **+0.0006** |
| donor route ratio | 54.2 | **54.2** |

Both of the brief's alternatives reproduce exactly. Instrument calibrated.

### 0.2 Rows carried *inside* every later screen

| row | expected | measured |
|---|---|---|
| oracle = the truth itself | cos 1, sign 1 | **1.0000 / 1.0000** |
| zero prediction | sign 0.5 | **0.5000** |
| sign-flipped truth | sign 0.0 | **0.0000** |
| drug-label permuted consensus (`M/perm-drug`) | sign 0.5, cos 0 | **0.5013 / +0.0029** (k=49) |
| single global rescale of the consensus (`scale-ORACLE`) | cos_orth exactly 0 | **+0.0001** |

**One of these caught a real bug.** The first centring implementation subtracted the mean within
each measurement's own plate-set group; 156 of the 1,137 treatments sit on more than one plate and
some plate sets have a single member in a line, so those rows were centred against themselves and
annihilated to exactly zero. The oracle row read **0.98 instead of 1.00** and aborted the run.
Centring was replaced by subtracting the mean of the per-`(line, plate)` blobs, each estimated
from that plate's ~94 drugs, which cannot annihilate a row. Without the oracle row this would
have silently deleted 2% of every fold.

---

## 1. Acquisition — the smallest sufficient representation (Phase 1)

The 100M-cell atlas is not the object needed and was never downloaded. Four routes were checked
in the brief's order:

| route | what it is | size | used |
|---|---|---:|---|
| `gs://arc-institute-virtual-cell-atlas/tahoe100M/2025-02-25/h5ad/` | 14 per-plate cell-level h5ad | **336 GB** | no — over budget, and the wrong object |
| HF `tahoebio/Tahoe-100M` `data/train-*` | 3,388 cell-level parquet shards | **338 GB** | no |
| HF `tahoebio/Tahoe-100M` `metadata/pseudobulk_differential_expression/` | 1,026 shards, DESeq2 per (drug, dose, line, plate) × 62,710 genes | **92 GB** | no — one 91 MB shard downloaded as a schema probe (`pdde-00000.parquet`, 3,986,181 rows, 16 columns) |
| **HF `tahoebio/tahoe-de-rhaister`** | the Rhaister training tables | | **yes** |

`tahoe-de-rhaister` is the right object and nothing in the brief's list pointed at it — it is a
*separate* repo from `tahoebio/Tahoe-100M`, found by searching HF for Tahoe-derived datasets.
Downloaded, 1.2 GB total:

| file | rows / shape | bytes |
|---|---|---:|
| `cell_eval/plate_plate{1..14}.parquet` | one row per (cell_line, treatment, plate), 2,000 gene columns of pseudobulk delta `mean(treated) − mean(plate-matched DMSO)` | 1,065 MB |
| `zeroshot/control_expression.parquet` | 43,897,000 rows = 700 (line, plate) × 62,710 genes of DMSO-control mean expression + library size + control cell count | 186 MB |
| `zeroshot/cell_line_centroids_dmso{,_hvg}.parquet` | per-line DMSO centroids | 11 MB |
| `metadata/pseudobulk_differential_expression/train-00000-of-01026.parquet` | schema probe only, unused | 91 MB |

### 1.1 The tensor (`build.py` → `work/tensor.npz`)

```
delta        (1137, 50, 2000) float32   sha256 68bf46f4e9e47aa00d2627e1d177bb257ed5ea65af82dd78e61d5f75c0437992
basal_line   (50, 2000)   float32       basal_plate (50, 15, 2000)
nrep, group  (1137, 50)                 ctrl_cells (50, 15)
```

- **1,137 treatments × 50 cell lines × 2,000 genes, 56,827 of 56,850 grid cells populated —
  99.96% complete.** Every drug is measured in essentially every line, because Tahoe's MOSAIC
  design pools all 50 lines in the same well, so line and treatment are not confounded by batch.
  Nothing in `data/donors/` is remotely this dense: our best pair is 3,000 shared perturbations
  across 2 lines.
- Treatments are (compound, concentration) pairs; 384 compounds × 3 doses.
- **7,800 (line, treatment) pairs are measured on 2–14 plates.** These replicates are the only
  same-line same-perturbation repeat structure in the dataset and they carry the reliability
  anchor (`work/reps.npz`, sha `408aad65…`).
- Scale check, run rather than assumed: `‖delta‖` per line is 0.86–1.15 and correlates
  −0.36 with that line's control library size, while `‖basal‖` is 6.7–17.7 and correlates
  strongly — so **the deltas are library-size normalised and the control means are not**. Basal
  is therefore converted to a library fraction before use.

---

## 2. Anchors — what Tahoe's cross-line structure looks like next to ours (`anchor.py`)

Everything is centred per `(line, plate)`: the line's mean response across drugs is the metric's
zero (159 §5.1) *and* the plate-matched DMSO's own noise enters every delta on that plate with the
same sign — the exact shape of the 86.5% blob that made 159's brief misread the Orion harnesses.

| quantity | Tahoe (this experiment) | our CRISPRi donors (159) |
|---|---:|---:|
| same-line, same-perturbation reliability | **+0.2808** (7,641 replicate pairs; +0.3984 uncentred) | ρ ≈ 0.26 (K562 two libraries); residual ρ 0.05–0.12 at donor depth |
| cross-line, same perturbation | **+0.2418** | +0.0306 (HCT116↔HEK293T) |
| cross-line, different perturbation (decoy) | **+0.0027** | +0.0016 |
| ratio | **88.8** | 19.5 |
| within-line, different perturbation (blob check) | **+0.0116** | +0.0025 / +0.0067 |
| disattenuated cross-line transfer τ | **0.72 – 0.86** | **0.30 ± 0.10** |
| consensus law `1/cos²(k) = A + B/k` | A 4.92, B 12.48, R² 0.990 | A 66.5, B 1,920 |
| `cos(k→∞)` | 0.4507 (0.8504 disattenuated) | 0.1226 |
| lines needed for half the ceiling | **0.84** | **9.6** |

Two readings, both load-bearing for how the rest is interpreted:

- **The surface is clean.** Within-line different-perturbation cosine is +0.0116 and the decoy
  level is +0.0027, so this dataset does not carry the harness pathology of 159 §2. Uncentred and
  centred replicate reliabilities differ by only 0.40 → 0.28, which is the size of the plate blob:
  real but small.
- **Tahoe is a much easier transfer problem than ours.** 72–86% of a drug's drug-specific
  response energy is line-independent, against ~30% for CRISPRi knockdowns; the consensus
  saturates at `k < 1` rather than needing ten lines. Drug responses at these doses are largely a
  shared MOA/stress signature. **So Tahoe is a substrate for measuring whether an operator is
  learnable, not a source of donor rows, and no absolute level here transfers to the board.**

---

## 3. Phase 2 — leave-one-LINE-out over all 50 lines (`phase2.py`, `work/phase2-main.json`)

50 held lines × k ∈ {3, 5, 10, 20, 40, 49} × 4 draws (1 at k=49) = 1,050 folds, 53 min, peak RSS
5.3 GB. Everything is fitted on training lines only: the rank-64 response basis, the basal
standardisation, the basal embedding (4 PCs), the gains, the operator, and the ridge parameter
(chosen by leave-one-*training*-line-out inside the training set).

**Arms.** `M` = mean of the training lines' centred responses (production's estimator,
transplanted). `gain-mean` = `M` reweighted per gene by a gain that is the same for every
recipient — 161's lesson, since a flexible learner pointed at this objective converges on a
context-free per-entry reshaping, so that reshaping is the control the context arm must beat.
`gain-basal` = the gain is a linear function of the recipient's basal log-expression of that gene
(`[1, ḡ, z, ḡ·z]`, two-stage: per-(line, gene) least-squares gain, then regressed on basal).
`gain-basal-clip` = the same, clipped at zero, because a recipient can attenuate or amplify a
response but not invert it. `bilinear` = the low-rank operator `Δ̂_L = M + (W₀ + Σ u_L[k] W_k)`
applied in the rank-64 response basis, parameterised as a *correction* to `M` so that `W = 0`
reproduces the baseline exactly and the ridge shrinks towards it. `-nc` suffix = the same operator
with the context block removed. `ORACLE` rows are fitted on the held line's own truth and are
ceilings, not arms.

**Endpoints.** `sign` = per-row sign accuracy on that row's own top-200 |truth| genes, macro over
rows then over lines (rule 8, never pooled). `sign_weak` = the same on the half of those genes
where the consensus is *weakest* — 174's framing, and the half where a context model would have to
earn its keep. `cos` = per-row cosine. `cos_orth` = 161's endpoint, the cosine between the arm's
and the truth's components **orthogonal to `M`** — a rescaling of the consensus cannot fake it, so
it measures line-specific information directly.

### 3.1 `cos_orth` — the line-specific residual, and the learning curve in k

| arm | k=3 | k=5 | k=10 | k=20 | k=40 | k=49 |
|---|---:|---:|---:|---:|---:|---:|
| **`gain-basal-clip`** | **+0.2293** | **+0.3159** | **+0.3891** | **+0.4244** | **+0.4557** | **+0.4635** |
| `gain-basal` | +0.1164 | +0.2167 | +0.3389 | +0.4021 | +0.4363 | +0.4439 |
| `bilinear` | −0.0003 | +0.0999 | +0.1475 | +0.1647 | +0.1848 | +0.2062 |
| `gain-mean` (context-free control) | +0.0026 | +0.0113 | +0.0221 | +0.0435 | +0.0351 | +0.0316 |
| `bilinear-nc` (context-free control) | +0.0376 | +0.0453 | +0.0374 | +0.0349 | +0.0263 | +0.0361 |
| `gain-basal/perm-line` (**null**) | +0.0195 | +0.0269 | +0.0316 | +0.0188 | +0.0089 | +0.0389 |
| `bilinear/perm-line` (**null**) | +0.0022 | +0.0239 | +0.0166 | −0.0067 | +0.0093 | −0.0028 |
| `bilinear/perm-drug` (**null**) | −0.0000 | +0.0011 | +0.0021 | +0.0022 | +0.0030 | +0.0031 |
| `M/perm-drug` (**null**) | +0.0000 | +0.0004 | +0.0006 | +0.0010 | +0.0003 | +0.0005 |
| `gain-ORACLE` (ceiling, leaky) | +0.3896 | +0.4217 | +0.4595 | +0.4809 | +0.5148 | +0.5224 |

**This is the headline.** The arm doubles from k=3 to k=49, both nulls stay flat at zero, and the
context-free control stays flat at 0.03. At k=49 the arm holds 89% of the diagonal oracle.

Spread across the 50 held lines (sd, so sem = sd/√50): `gain-basal-clip` 0.0962 / 0.1009 / 0.1184
at k = 3 / 5 / 49; `gain-basal/perm-line` 0.0348 / 0.0583 / 0.1016. The k=49 null's own sem is
0.014, so it is not distinguishable from zero, and the arm sits ~19 combined sem above it.

### 3.2 `cos` — delta against `M`, with the count of held lines that improve

| arm | k=3 | k=5 | k=10 | k=20 | k=40 | k=49 |
|---|---|---|---|---|---|---|
| **`gain-basal-clip`** | +0.0194 (33/50) | +0.0871 (48/50) | +0.1241 (48/50) | +0.1373 (49/50) | +0.1437 (49/50) | **+0.1467 (49/50)** |
| `gain-basal` | −0.1875 (0/50) | −0.0703 (8/50) | +0.0444 (39/50) | +0.0974 (49/50) | +0.1062 (49/50) | +0.1087 (49/50) |
| `bilinear` | −0.3262 (0/50) | +0.0044 (26/50) | +0.0239 (46/50) | +0.0283 (47/50) | +0.0310 (48/50) | +0.0379 (48/50) |
| `gain-mean` | −0.0760 (0/50) | −0.0355 (4/50) | −0.0080 (8/50) | +0.0017 (28/50) | +0.0010 (30/50) | +0.0008 (31/50) |
| `gain-basal/perm-line` (null) | −0.2691 (0/50) | −0.2333 (0/50) | −0.1849 (0/50) | −0.1687 (0/50) | −0.1749 (0/50) | −0.1643 (1/50) |
| `bilinear/perm-drug` (null) | −0.3197 (0/50) | −0.0516 (0/50) | −0.0306 (0/50) | −0.0324 (0/50) | −0.0393 (0/50) | −0.0524 (0/50) |
| `gain-ORACLE` (ceiling) | +0.1627 (50/50) | +0.1661 (50/50) | +0.1721 (50/50) | +0.1738 (50/50) | +0.1812 (50/50) | +0.1838 (50/50) |

Absolute levels: `M` runs 0.3195 → 0.4599, `gain-basal-clip` 0.3389 → 0.6066, oracle 0.4822 →
0.6437. **The arm captures 80% of the oracle's cosine gain and wins on 49 of 50 lines**, while its
line-permutation null loses 0.16.

### 3.3 `sign` and `sign_weak` — the gate endpoint, delta against `M`

| arm | metric | k=3 | k=5 | k=10 | k=20 | k=40 | k=49 |
|---|---|---|---|---|---|---|---|
| `bilinear` | sign | −0.2052 (0/50) | −0.0048 (9/50) | **+0.0036 (37/50)** | **+0.0037 (35/50)** | +0.0019 (31/50) | +0.0025 (33/50) |
| `bilinear` | sign_weak | −0.1093 (1/50) | −0.0027 (20/50) | **+0.0080 (41/50)** | **+0.0085 (41/50)** | +0.0050 (37/50) | +0.0063 (37/50) |
| `bilinear-nc` | sign | −0.0134 | −0.0263 | +0.0005 | +0.0005 | +0.0002 | +0.0009 |
| `gain-basal-clip` | sign | −0.0293 | −0.0207 | −0.0154 | −0.0126 | −0.0130 | −0.0133 |
| `gain-basal` | sign | −0.0587 | −0.0414 | −0.0307 | −0.0251 | −0.0260 | −0.0266 |
| `gain-ORACLE` | sign | +0.0004 | −0.0000 | −0.0001 | −0.0001 | −0.0003 | −0.0003 |
| `bilinear/perm-line` (null) | sign | −0.1901 | −0.0156 | −0.0153 | −0.0207 | −0.0256 | −0.0333 |
| `bilinear/perm-drug` (null) | sign | −0.2010 | −0.0317 | −0.0204 | −0.0228 | −0.0288 | −0.0376 |

Absolute `sign` for `M`: 0.7012 → 0.7984. `sign_weak` for `M`: 0.6071 → 0.7175.

**The gate is +0.010 absolute sign accuracy. The best number in the whole run is +0.0037**
(`bilinear`, k=20), reaching +0.0085 on the weak-consensus half. It beats both nulls comfortably
(the nulls are −0.02 to −0.04) and it is positive on 35–41 of 50 lines, so it is not noise — it is
just small, and it does **not** grow with k past k=10.

**Why, structurally.** `gain-ORACLE` — a per-gene gain fitted on the held line's *own truth* —
moves sign accuracy by **−0.0003**. A diagonal gain is sign-preserving wherever it is positive, so
the entire diagonal family, oracle included, is worth **zero** on this endpoint by construction.
Only a rotation can move sign, and the rotation is where the information is thin. This is not a
tuning failure; it is a property of the decomposition. The line-specific residual on Tahoe is
**an amplitude object**, and the gate asks about direction.

Rank and context-dimension sweep at k=49 (4 held lines): `cos_orth` for `bilinear` is 0.129
(r=32) / 0.139 (r=64) / 0.143 (r=128), and 0.139 (4 basal PCs) / 0.145 (8). Rank 64 / 4 PCs was
used throughout. Composing the operator on top of the gain-corrected consensus was tried and
dropped: it inherits the gain's `cos_orth` while making the line-permutation null uninterpretable
(arm 0.4055 against a null of 0.3902, because the null permutes only the operator's context and
not the gain's).

---

## 4. Mechanism — it is not a detection floor (`mechanism.py`)

The deflationary reading of a basal-conditioned gain is that a gene the recipient barely expresses
cannot show a large delta, so "predicting" it is predicting a zero. That would still be usable on
the challenge — we are given the recipient's control cells — but it is not recipient-conditional
biology. Stratifying the gene axis by the **minimum** basal expression across all 50 lines, at
k=49:

| stratum | genes | `gain-basal-clip` | `gain-mean` (ctx-free) | `perm-line` (null) | `gain-ORACLE` |
|---|---:|---:|---:|---:|---:|
| all genes | 2,000 | +0.4635 | +0.0316 | +0.0040 | +0.5224 |
| min basal ≥ median | 1,000 | +0.4707 | +0.0236 | −0.0255 | +0.5259 |
| min basal ≥ p75 | 500 | +0.4815 | +0.0096 | +0.0101 | +0.5307 |
| **min basal ≥ p90** | **200** | **+0.4817** | +0.0008 | +0.0305 | +0.5185 |

**The effect is strongest, not weakest, on the 200 genes that are well expressed in every one of
the 50 lines**, where no detection floor can operate, and there the context-free control collapses
to +0.0008. The law is: *the recipient's expression level of a gene sets how far that gene moves*,
among genes every line expresses. A single global per-line rescale (`scale-ORACLE`) scores
`cos_orth` +0.0001 — the effect is per-gene, not per-line amplitude.

---

## 5. The same measurement on our CRISPRi donors (`crispri.py`)

Not Phase 3: no operator is transferred, no `h1-test` surface is touched, no `pds` is computed and
nothing is built. This is the k = 2–5 end of the learning curve **measured on our own data**
instead of extrapolated, using the identical arms, nulls and endpoints. Donor tables read-only
from `data/donors/` and `experiments/166-donor-acquisition/donors/`; the across-target mean is
removed per line; each target's own gene coordinate is dropped.

`core5` = K562-gwps / RPE1 / Jurkat / HCT116 / HEK293T — **161's exact surface**, 1,883 shared
targets × 4,980 genes. `core6` adds 166's `nadig_hepg2`: 1,758 targets × 4,938 genes.

### 5.1 `cos_orth`, `core5`

| arm | k=2 | k=3 | **k=4** |
|---|---:|---:|---:|
| `gain-basal-clip` | +0.0232 | +0.0625 | **+0.0665** |
| `gain-basal` | +0.0141 | +0.0613 | +0.0662 |
| `gain-mean` (context-free control) | +0.0617 | +0.0652 | **+0.0667** |
| `gain-basal/perm-line` (**null**) | +0.0467 | +0.0646 | **+0.0663** |
| `gain-ORACLE` (ceiling) | +0.1064 | +0.1098 | +0.1126 |
| `M/perm-drug` (null) | +0.0037 | +0.0040 | +0.0059 |

`core6` at k=5 reads the same: arm **+0.0694**, context-free **+0.0696**, null **+0.0695**,
oracle +0.1224.

**Arm = context-free control = wrong-recipient null, to three decimals, on both surfaces.** The
recipient-specific content of the basal-conditioned gain on CRISPRi donor tables is zero. The
+0.066 that all three share is the context-free per-entry reshaping 161 already identified and
experiments 030/168/169 already own.

**This reproduces 161 from a different model class.** 161 reported `cos(orth, TRUE orth)` of
0.0592 (consensus × context) and 0.0606 (target-side only) from a 27-feature gradient-boosted
model on the same five lines; a two-stage closed-form diagonal gain lands at 0.0665/0.0667. Two
unrelated estimators, same number.

Sign and cosine tell the same story: on `core5` at k=4, `gain-basal-clip` moves sign by −0.0011
and cosine by +0.0092 — and `gain-mean`, which cannot see the recipient at all, moves cosine by
+0.0093. Every gain on CRISPRi is context-free.

### 5.2 Side by side at matched k

| | Tahoe, k=3 | Tahoe, k=5 | CRISPRi, k=3 | CRISPRi, k=5 (core6) |
|---|---:|---:|---:|---:|
| arm `cos_orth` | **+0.2293** | **+0.3159** | +0.0625 | +0.0694 |
| context-free control | +0.0026 | +0.0113 | +0.0652 | +0.0696 |
| line-permutation null | +0.0195 | +0.0269 | +0.0646 | +0.0695 |
| **arm ÷ null** | **11.8×** | **11.7×** | **0.97×** | **1.00×** |
| oracle ceiling | +0.3896 | +0.4217 | +0.1098 | +0.1224 |

Production has three donor lines. At k=3 Tahoe already shows a 12× separation from its null and
CRISPRi shows none. **Whatever is missing from our CRISPRi transfer problem, it is not lines.**

---

## 6. Is the CRISPRi null a noise floor? No (`noisefloor.py`)

Tahoe pseudobulks are 3–5× more reliable than our donor tables. If the same law were present in
CRISPRi and merely buried, deeper references would reveal it, and that is a completely different
strategic conclusion. Discriminated by degrading Tahoe: with observed energy split as signal
`S = ρV` and noise `(1−ρ)V`, adding independent noise of energy `S(1/ρ' − 1/ρ)` moves the
reliability from ρ to ρ'. The added noise carries each gene's own residual sd, so only the SNR
changes. Achieved reliability verified by re-measuring the replicate-plate cosine after noising.

| ρ target | ρ achieved | k | `M` cos | arm `cos_orth` | ctx-free | **null** | oracle | **arm ÷ null** |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.281 (native) | 0.281 | 5 | 0.365 | +0.3168 | −0.0055 | +0.0348 | +0.4325 | 9.1× |
| 0.281 | 0.281 | 49 | 0.460 | +0.4635 | +0.0316 | +0.0114 | +0.5224 | 40.7× |
| 0.120 | 0.115 | 5 | 0.126 | +0.1043 | +0.0135 | +0.0107 | +0.1478 | 9.7× |
| 0.120 | 0.115 | 49 | 0.210 | +0.1828 | +0.0362 | +0.0013 | +0.2148 | 141× |
| 0.080 | 0.076 | 5 | 0.083 | +0.0676 | +0.0117 | +0.0057 | +0.1023 | 11.9× |
| 0.080 | 0.076 | 49 | 0.158 | +0.1362 | +0.0339 | +0.0075 | +0.1619 | 18.2× |
| 0.050 | 0.047 | 5 | 0.059 | +0.0441 | +0.0104 | +0.0044 | +0.0712 | 10.0× |
| 0.050 | 0.047 | 49 | 0.114 | +0.0972 | +0.0298 | +0.0081 | +0.1174 | 12.0× |

Our Orion donors sit at residual reliability 0.05–0.12 (159 §1.3). **At ρ = 0.076, k=5 — inside
that band — noise-matched Tahoe reads arm +0.0676 against a null of +0.0057, a 12× separation.
CRISPRi at the same k and the same band reads 1.0×.** Note also that the *oracle* ceilings match
once noise is matched (Tahoe +0.102 at ρ 0.076 / k=5 against CRISPRi +0.113 at k=4), so the two
surfaces have comparable amounts of diagonal structure to find. The difference is entirely in
whether the recipient's basal profile predicts it.

---

## 7. The model-free version (`gaincorr.py`)

No ridge, no basis, no arm. Take the empirical per-gene gain `g_L[j]` (least-squares coefficient
of line `L`'s response on the other lines' consensus), remove row and column means from the
`(line, gene)` matrix, do the same to log basal, and correlate what is left.

| surface | lines | genes | pooled r | per-line r, mean [min, max] |
|---|---:|---:|---:|---|
| **Tahoe** | 50 | 2,000 | **+0.666** | +0.685 [+0.055, +0.780] |
| Tahoe, noised to ρ 0.115 | 50 | 2,000 | +0.673 | +0.689 |
| Tahoe, noised to ρ 0.076 | 50 | 2,000 | +0.668 | +0.684 |
| Tahoe, noised to ρ 0.047 | 50 | 2,000 | +0.640 | +0.662 |
| **CRISPRi `core5`** | 5 | 4,980 | **+0.166** | +0.186 [+0.091, +0.325] |
| CRISPRi `core6` | 6 | 4,938 | +0.167 | +0.183 [+0.077, +0.309] |

The statistic is essentially **noise-invariant** — the gain averages over 1,137 drugs, so a 6×
degradation of each individual measurement barely moves it. **45% of the line-specific per-gene
gain variance is explained by the recipient's basal profile for drug responses; 3% for CRISPRi
knockdowns.** That is the whole experiment in one row, and it does not depend on any model choice
in §3.

---

## 8. The gate, stated plainly, and the one concern with it

**Gate:** `T` beats `M` by ≥ +0.010 absolute sign accuracy, macro over held lines, beating both
nulls. **Measured best: +0.0037** (`bilinear`, k=20; +0.0085 on the weak-consensus half). Both
nulls are beaten, the win is on 35–41 of 50 lines, and it is flat in k past k=10. **The gate
fails. Phase 3 was not run and no operator was transferred to CRISPRi.**

The concern, stated once and not used to argue the gate away: the gate is on the one endpoint the
winning mechanism cannot move. `gain-ORACLE` — fitted on the held line's own truth — moves sign
accuracy by −0.0003, so **the entire diagonal family has a sign-accuracy ceiling of zero by
construction**, and the arm that does pass everything else (`gain-basal-clip`: `cos_orth` +0.4635
against a +0.0389 null that is itself indistinguishable from zero, `cos` +0.1467 on 49/50 lines)
is diagonal. On cosine and on the line-specific residual the gate would have passed by a wide
margin.

That does not change the recommendation, because §5–§7 make the Phase-3 question moot: the object
Phase 3 would transfer does not exist in the receiving data. An amplitude operator learned on
Tahoe would be applied to CRISPRi donor effects whose line-specific gains correlate with basal at
r = 0.17, where arm ≡ null. The right reading is not "the operator failed a threshold" but
**"the operator is real, it is an amplitude law, and CRISPRi knockdown responses do not obey it."**

---

## 9. What this closes, and what it leaves open

**Closed.** The line-specific residual of a CRISPRi knockdown response is not predictable from the
recipient's basal profile, and the failure is not a sample-size failure. It joins 155, 067, 017,
104, 160, 161 and 174 as another route into 159/160's wall, and is the first with the statistical
power to distinguish "no signal" from "not enough lines". 161's own stated limit — *"this experiment cannot
rule out a residual signal that 400-cell references would reveal"* — is now partly answered from
the other side: at CRISPRi reliability the estimator finds a 12× null separation when the signal
is present, so its silence on CRISPRi is informative. **Do not fund another cross-line residual
predictor on donor-table basal features.**

**Open, and worth someone's time.**

1. **The amplitude law is real for chemistry and the challenge has three chemical-free contexts —
   but it is also the shape of what 187 owns.** `gain-basal-clip` is a per-gene multiplicative
   rescaling of the emitted effect driven by the recipient's control-cell expression. On the
   challenge that would move `mse` and `nmae` (which read the profile) and leave `reach`/`fid`/`jac`
   untouched, since it is sign-preserving. Its CRISPRi recipient-specificity is zero, so it should
   **not** be shipped as a context adapter — but the *context-free* half of it (`gain-mean`,
   `cos_orth` +0.0667 on core5) is the same object as 030/168/169's sharpening, which is already
   priced elsewhere. Nothing here is proposed for the builder.
2. **`τ` differs by 2.5× between chemical and genetic perturbation** (0.72–0.86 vs 0.30). If a
   competitor's Tahoe-derived pretraining (PertMind, `luka`, board `reach` 0.2908) is transferring
   anything, it is transferring a *shared* response manifold, not a context adapter — Tahoe's
   context adapter is an amplitude law that does not survive the jump to knockdowns. That is a
   concrete, testable claim about the leader's mechanism, and 186 owns competitive intel.
3. **Untouched by budget:** the 92 GB `pseudobulk_differential_expression` table carries DESeq2
   `padj` on all 62,710 genes, so "the held line's own *significant* genes" could be defined by an
   actual test rather than by top-|truth|; and the `pdex` long table (40.6 GB) carries
   Mann–Whitney p-values. Neither would change §5–§7, which are model-free.

---

## 10. Limits — what is surface-limited here (rule 9)

- **`fid`, `jac`, `reach` were not computed and nothing here should be read as a claim about
  them.** This surface is pseudobulk: there are no cells, no Wilcoxon and no call set, and all
  three carry `n_real` in a denominator or a union.
- **`pds` was not computed either**, on any surface. The endpoints are cosine, orthogonalised
  cosine and per-perturbation sign accuracy.
- **No absolute level from Tahoe transfers to the board.** Its τ is 0.72–0.86 against our 0.30,
  its baseline sign accuracy is 0.80 against our incumbent 0.537, and its gene axis is a 2,000-gene
  analysis set, not the challenge's 18,533. Only the *contrasts* — arm vs null, arm vs
  context-free, and the shape of the k-curve — are claimed to carry.
- **`sign` is measured against a noisy truth** (ρ = 0.28 native, 0.05–0.12 on CRISPRi), which
  attenuates every arm's sign accuracy. It attenuates arms and nulls alike, so deltas are readable,
  but the absolute levels are not.
- **The k-curve's basis and operator are refitted at every k**, so the curve is honest about the
  operator's k-dependence; the basal *standardisation* at k=3 is estimated from 3 lines and is
  correspondingly noisy, which works against the low-k end of the curve.
- **The CRISPRi comparison uses global centring** (donor tables have no plate structure) while
  Tahoe uses per-`(line, plate)` centring. Global centring on Tahoe was run as a variant at k=49
  (`work/phase2-global.json`) and changes nothing material: `gain-basal-clip` `cos_orth` **+0.4554**
  against a null of **−0.0182** and a context-free control of **+0.0210**, oracle +0.5372;
  `bilinear` `cos_orth` +0.2009. The one number that moves is the gate endpoint, and it moves
  *down*: `bilinear` sign delta **+0.0010** under global centring against +0.0025 under plate
  centring. The gate fails either way.
- **`nadig_hepg2` and `frangieh_melanoma`** were built by 166 and carry 0/300 and 5/300 panel
  targets; `nadig_hepg2` is used here only as a sixth *training* line, which is exactly what 166
  built it for.

---

## 11. Files

| file | what it does |
|---|---|
| `known.py` | §0.1 — Phase 0, reproduces 159's harness and donor rows from frozen artifacts |
| `build.py` | §1.1 — the tensor, the replicate store, the per-(line, plate) blobs, shas |
| `common.py` | loading and the centring convention |
| `anchor.py` | §2 — reliability, transfer level, decoy level, consensus law, τ |
| `phase2.py` | §3 — the arms, the nulls, the k-curve, all four endpoints |
| `report.py` | §3 tables |
| `mechanism.py` | §4 — the expression-floor stratification |
| `crispri.py` | §5 — the same measurement on our donor tables |
| `noisefloor.py` | §6 — the reliability-matched discrimination |
| `gaincorr.py` | §7 — the model-free correlation |
| `work/*.json` | every number above |

```bash
cd experiments/184-tahoe-transfer
uv run python known.py                       # Phase 0 first; abort if it does not pass
uv run python build.py                       # ~20 s, peak RSS 7.4 GB
uv run python anchor.py                      # ~25 s
uv run python phase2.py --ks 3,5,10,20,40,49 --draws 4 --rank 64 --nctx 4 \
    --out work/phase2-main.json              # 53 min, peak RSS 5.3 GB
uv run python report.py work/phase2-main.json
uv run python mechanism.py                   # ~4 min
uv run python crispri.py                     # ~50 s
uv run python noisefloor.py                  # ~9 min
uv run python gaincorr.py                    # ~1 min
uv run python phase2.py --ks 49 --draws 1 --rank 64 --nctx 4 --centring global \
    --out work/phase2-global.json            # §10 robustness, ~7 min
```

`work/tensor.npz` (394 MB) and `work/reps.npz` (114 MB) are regenerable from
`data/external/tahoe/` in 20 s.

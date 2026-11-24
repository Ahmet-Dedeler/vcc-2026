# 280 — PREREGISTRATION: the biological donor as a replicate unit in GWCD4i

**Written 2026-09-06, 07:05Z (10:05 Ankara), before any expression value from any GWCD4i cell-level or
pseudobulk object had been read.** Everything below is fixed by this document; the hash of this file is
recorded in `work/prereg280.sha256` and quoted in `RESULTS.md`. Facts used to write it come from
**metadata only** — the S3 key listing, `data_sharing_readme.md`, `sample_metadata.suppl_table.csv`,
`QC_summaries_per_sample_lane.csv`, the GEO `filelist.txt`, and the `.obs`/`.var`/`X.indptr` groups of
`GWCD4i.pseudobulk_merged.h5ad` read over HTTP range requests (`work/obs280.npz`). **No `X` value was read.**

---

## 0. The experimental unit (programme rule 1)

| unit | levels | note |
|---|---|---|
| **biological donor** (the unit this experiment adds) | 4: `CE0008162`=**D1** (29 F, Hispanic), `CE0010866`=**D2** (22 F, African American), `CE0008678`=**D3** (23 F, African American), `CE0006864`=**D4** (34 M, Asian) | person-to-person; no fold on disk has one |
| culture condition / cell state | 3: `Rest`, `Stim8hr`, `Stim48hr` | a within-donor state contrast |
| **10x run (processing batch)** | 2: `CD4i_R1`, `CD4i_R2` | **CONFOUNDED with donor in Rest and Stim8hr** (see §0.1) |
| 10x lane | 284 lanes | collapsed in the pseudobulk object |
| **guide construct** | 1–2 per (target, donor, condition) | F1: a real unit of variance (+0.014 Flex) |
| cell | median **101** cells per (panel target, donor, condition); **482** pooled over 4 donors | q10 16 / q90 294 |

### 0.1 The confound that decides which condition is primary — declared before any result

From `sample_metadata.suppl_table.csv` and the pseudobulk `obs`:

| condition | D1 | D2 | D3 | D4 | donor ⟂ run? |
|---|---|---|---|---|---|
| Rest | R1 | R1 | R2 | R2 | **NO** — the {D1,D2} vs {D3,D4} half is exactly the run split |
| Stim8hr | R1 | R1 | R2 | R2 | **NO** — same |
| **Stim48hr** | **R2** | **R2** | **R2** | **R2** | **YES** — all four donors in one 10x run |

**Therefore `Stim48hr` is the PRIMARY condition for every donor claim.** Rest and Stim8hr are reported as
secondary, and within them the pairing {D1,D2} vs {D3,D4} is labelled **run-confounded** while {D1,D3} vs {D2,D4}
and {D1,D4} vs {D2,D3} are run-balanced (one R1 donor and one R2 donor on each side). Prediction **P4** below is
the falsifiable form of this.

---

## 1. Acquisition plan (fixed before downloading anything)

**Rejected — the cell-level files.** `D*_*.assigned_guide.h5ad` are **118.6–172.8 GB each**
(D1_Rest 142,828,875,662 B; the four Rest files sum to **540.2 GB**). One file alone is 2–3× the mandate's
60 GB cap. The GEO route (`GSE314342_RAW.tar`, 171,154,268,160 B = 284 CellRanger
`sample_filtered_feature_bc_matrix.h5`; Rest alone 54.58 GB over 94 files) is under the cap only for one
condition, carries **no published guide assignment**, and would require re-deriving it. **Not acquired.**

**Chosen — a range-request extraction from `GWCD4i.pseudobulk_merged.h5ad`.** The object is 44,566,657,140 B,
`X` is **CSR, uncompressed** (`data` float64, `indices` int64, `indptr` int64; 278,684 rows × 18,129 genes;
2,782,896,823 nnz), so individual rows are addressable by byte range. The rows needed are

* every **targeting** pseudobulk whose `perturbed_gene_name` is one of the 300 panel targets (6,477 rows), and
* every **non-targeting** pseudobulk (11,018 rows) — the per-(donor, condition) control pool,

= **17,495 rows, 177,782,218 nnz ⇒ 2.84 GB transferred**, stored locally as a dense `float32` bank of
17,495 × 18,129 = **1.27 GB**. Peak disk in `work/` stays **under 2 GB**; nothing else is downloaded.
Terms: **MIT licence** (`https://opensource.org/licenses/MIT`), stated on the CZI Virtual Cells Platform page
and in the bucket's `data_sharing_readme.md`. Bucket `s3://genome-scale-tcell-perturb-seq/marson2025_data/`,
anonymous HTTPS. Every retained file's sha256 is recorded.

**Optional secondary (only if the primary finishes inside the time box):** the same range-request extraction of
the panel rows from `GWCD4i.DE_stats.by_donors.h5mu` (16,866,278,447 B) — the authors' own DESeq2 fits computed
**within each donor pair** — as an independent estimator for Q1. Declared here so that using it is not a
post-hoc choice; if it is not reached, that is stated.

---

## 2. The response object (267's residual space, arithmetic fixed here)

For a **side** = a set of pseudobulk rows (a choice of donors × guides for one target, in one condition):

1. `S` = elementwise sum of the side's raw UMI count rows (18,129 genes); `tot = S.sum()`.
2. `p = log1p(TS * S / tot)` with **`TS = 5e4`** — `build_donor_table.py` / `truth267.py` arithmetic.
3. `b` = the same transform applied to the summed **non-targeting** rows of **the same donor(s) and the same
   condition** (so each side's own donor basal state is removed, and the statistic is a *response*, not a
   basal difference).
4. `delta = p − b`.
5. **Residual centring:** `r = delta − mean(delta over the 297 panel targets built the same way on that same
   side)`; then the **300 panel gene positions are set to 0** (`pds_cosine` removes all 300).
6. Scored columns = the official 18,533-gene axis (`data/challenge/gene_names.csv`) ∩ the 18,129 CD4 genes,
   minus the panel positions. All cosines are on those columns.

The **un-centred** ("full") cosine of step 4 is reported alongside every centred number, exactly as
`truth267.py` does, so the two conventions are never mixed.

**Cell-count matching** is done by **multinomial thinning of the pseudobulk UMI counts** to a common per-side
total, applied *before* steps 2–5 (so the centring is computed on the thinned rows too). Thinning equalises
**counting** noise at equal depth; it does **not** reproduce cell-to-cell biological variance, and every
thinned number is labelled `counting-matched`, never "cell-matched", except where the sides genuinely carry
equal numbers of cells by construction (the 2×2 of §3.1, where matching is by cell count band **and** thinning).

Seed for every random operation: **`20280906`**, one `numpy.random.default_rng` per statistic, recorded.

---

## 3. Q1 — how reproducible is a deep Flex response across human donors?

### 3.1 The 2×2 that isolates the donor (primary)

For each panel target `t` present with 2 guides in ≥ 2 donors, and each condition, form the four classes of
pairs over (donor, guide) pseudobulks, both sides counting-matched to a common UMI total:

| class | pair | variance it contains |
|---|---|---|
| **D** | same donor, same guide, binomial(0.5) count split | counting noise only (**floor**) |
| **A** | **same donor, different guide** | construct + cells + counting |
| **B** | **different donor, same guide** | donor + cells + counting |
| **C** | **different donor, different guide** | donor + construct + cells + counting |

* **Primary pass rule — "DONOR IS A REAL UNIT":** `Δ_donor = median_t cos_A(t) − median_t cos_C(t) > 0` with a
  **paired target bootstrap** (10,000 resamples over targets, the unit resampled is the target) whose 95 %
  interval **excludes 0**, computed in **Stim48hr**. Construct is held constant as a nuisance present on both
  sides; only the donor differs between A and C.
* Secondary readings, same machinery: `D − B` (donor beyond counting noise), `D − A` (the construct's own
  contribution — an independent replication of F1's construct finding on a 4th surface, in primary human cells).

### 3.2 The mandate's donor halves (headline reliability)

For each condition and each of the three **disjoint** pairings — {D1,D2}|{D3,D4}, {D1,D3}|{D2,D4},
{D1,D4}|{D2,D3} — pool all cells and all guides on each side (~200 cells/side) and report the per-target
cosine. Headline = **median across-donor residual cosine over the 297 panel targets, pooled over the three
pairings, in Stim48hr**. Comparator at the same per-side cell count: class **A** of §3.1 (~50 cells/side) and
class **D** (the counting floor), both reported with their measured cells/side so the mismatch is auditable,
plus the thinned version of §3.2 brought down to A's per-side total.

### 3.3 Predictions (stated before looking; each is falsifiable)

* **P1.** The Stim48hr across-donor median (§3.2) lands in **[+0.08, +0.22]** — above F1's Orion 150–200-cell
  regime (+0.028) and **below** the H1 Flex lane split at a comparable per-side cell count
  (+0.243 at 200/200; +0.114 at 80/80). Direction claimed: **below +0.243**.
* **P2.** `Δ_donor = median cos_A − median cos_C > 0` (donor is a real unit), and the margin is **small** —
  point estimate in **[+0.005, +0.06]** — because F1 already shows the construct is itself a real unit
  (+0.014) and because at ~50 cells/side counting noise dominates both classes.
* **P3.** Across-donor cosine rises with 276's conservation: **Spearman(across-donor cosine, `C_raw`) > 0**,
  predicted **+0.15 to +0.40**, on the 297 panel targets joined to `work/conservation276.csv`.
* **P4.** In **Rest** and **Stim8hr** the run-confounded pairing {D1,D2}|{D3,D4} reads **lower** than the mean
  of the two run-balanced pairings; in **Stim48hr** the three pairings agree (all differences inside their
  paired bootstrap intervals).
* **P5 (null).** The label-shuffled cosine — side A's row for target `t` against side B's row for a random
  `t′ ≠ t` — is **0.00 ± 0.02**. If it is not, the residual centring is wrong and Q1 is void.

---

## 4. Q2 — replicates vs cells (programme rule 7's depth/replicate axis)

**Leave-one-donor-out.** With 4 donors the largest honest training set is 3, so the curve is `k ∈ {1,2,3}`
training donors, not "four-donor mean vs one" — that framing is recorded here as unavailable and is replaced.

For each held-out donor `h` and each target `t`: truth = `h`'s own residual row (all of `h`'s guides and cells).
Prediction = the **mean of the residual rows** of a size-`k` subset of the other three donors, averaged over all
`C(3,k)` subsets. Score = `cos(pred, truth_h)`, averaged over the 4 choices of `h`.

* **Full-cells curve:** every training donor contributes all its cells.
* **Counting-matched curve (the one that answers the question):** each of the `k` training donors is
  multinomially thinned to `1/k` of the per-donor UMI budget, so the **total** counts behind the prediction are
  equal at `k = 1, 2, 3`. `Δcos = cos(k=3) − cos(k=1)` at matched total counts is the reported quantity, with a
  paired target bootstrap.
* Reported per condition, and split by 276's `C_raw` terciles.
* **P6.** At matched total counts `Δcos = cos(k=3) − cos(k=1) > 0` — spreading a fixed budget over more donors
  beats concentrating it in one — with a point estimate in **[+0.01, +0.08]**. If the interval contains 0, the
  answer is "replicates buy nothing beyond cells at this budget", which is a valid result.
* **P7.** The full-cells curve rises monotonically in `k` (trivially, 3× the cells) — this is a sanity check,
  not a finding; if it does not rise, the pipeline is broken.

---

## 5. Known-answer rows (every screen carries one; abort on failure)

| row | required | why it can fail |
|---|---|---|
| **KA-1** | `n_cells` summed over a panel target's **`keep_for_DE`** guides across all 4 donors in condition `c`, from the pseudobulk `obs`, **equals** `DE_stats.n_cells_target` on the local `experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad` (`c355f535…`) for that (target, condition). **Exact, max\|Δ\| = 0** | wrong keep-filter, wrong donor set, wrong target curation |
| **KA-2** | the **pooled four-donor** Rest residual row agrees with the on-disk `data/external/cd4t/cloud-compact/cd4_compact/cd4_rest` bank row for the same target. The bank's `delta` is `ln(2) × DESeq2 log2FC`, a *different estimator* from a log1p-CPM pseudobulk difference, so exactness is **not** expected. **Pass rule fixed here: median per-target cosine ≥ 0.50** over the panel targets with ≥ 400 pooled cells, and > 0 for ≥ 95 % of them. Max \|Δ\| reported with its explanation | a wrong gene axis, a wrong control pool, a sign error |
| **KA-3** | the join to `experiments/276-conservation-map/work/conservation276.csv` reproduces `C_raw` for every matched target, max\|Δ\| = 0 (a regression check) | a bad merge |
| **KA-4** | **P5**'s label-shuffled null ≈ 0 | broken residual centring |
| **KA-5** | 278 §1.1's CD4 depth row: pooled-4-donor cells per (panel target, condition) has median 482 (all guides) — reconcilable with 278's DE-table medians 512 / 490 / 495.5 | a wrong row selection |

---

## 6. What this experiment will NOT do

No model, no build, no submission, nothing under `submissions/` or `data/`, no production file touched, no
`cell-eval2`, no `vcc`. Nothing here is priced for a slot; a positive Q1/Q2 is a **scientific** result about the
replicate ceiling of a Flex-chemistry primary-cell reference, and any competition use would need its own
measured transfer (programme, "two evaluation tracks").

**Descriptive, post-hoc, labelled as such:** which targets are donor-variable, and whether donor variability
tracks 276's low-conservation tail. These are hypothesis-generating and carry no pass rule.

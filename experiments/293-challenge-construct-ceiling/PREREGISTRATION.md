# 293 — PREREGISTRATION: the guide-construct ceiling in the challenge's own regime

Written **before any expression value was read**. Only `obs` metadata (`guide_id`, `target_gene`, `batch`)
and `pert_counts_*.csv` have been read at the time of hashing; `work/inventory293.json` (metadata only) is
the sole numeric input so far. Hashed into `work/PREREGISTRATION.sha256` before any cosine is computed.

Timestamp of writing: 2026-09-06, experiment `experiments/293-challenge-construct-ceiling/`.

---

## 0. The question

The 2026 challenge scores, per target, **400 cells of ONE dual-guide CRISPRi construct** in six somatic
lines, zero-shot. 285 (F14) measured on five *public* CRISPRi libraries that the guide construct is the
dominant irreducible term of a perturbation response (ρ_guide 0.09 K562-GWPS / 0.17 Jurkat / 0.41 HepG2 /
0.02–0.06 CD4 / 0.76–0.81 Broad-16 Cas9). None of those libraries is Arc's. The 2026 panel shares only
25/300 targets with the 2025 H1 library, so no source on disk carries the 2026 constructs and the construct
ceiling binds every zero-shot entry.

**What is that ceiling in the challenge's own regime — Arc's own dual-guide library, H1 hESC, 10x Flex,
~400 cells, the same knockdown-validated construct design the 2026 truths use?**

273 read it on the 5 two-construct targets of the *test* split (median +0.0135; n = 5). This experiment
reads all 17 two-construct targets across the three 2025 splits on one matched design.

---

## 1. Inventory (already measured, metadata only — `work/inventory293.json`, `work/inv293.py`)

A **construct** is the dual-guide vector id: `GENE_Px_A|GENE_Px_B` → `GENE_Px`. `GENE_P1P2_A|GENE_P1P2_B`
is ONE construct. Parser check: constructs = 31 NTC + (targets + multi-construct extras) in every split.

| split | cells | targets | constructs | NTC constructs | NTC cells | Flex lanes | ≥2-construct targets |
|---|---:|---:|---:|---:|---:|---:|---:|
| train | 221,273 | 150 | 189 | 31 | 38,176 | 48 | **8** |
| validation | 98,927 | 50 | 85 | 31 | 38,176 | 48 | **4** |
| test | 170,846 | 100 | 136 | 31 | 38,176 | 48 | **5** |

n = **17** (13 as the mandate expected + validation's 4). The 38,176 non-targeting cells and the 31 NTC
constructs are the **same pool in all three splits**. Every construct spans all 48 lanes except `DPF2_P1`
(42 lanes, 95 cells). Every target gene is present on the 18,080-gene axis in every split.

The 17 targets and their two constructs (cells):

| split | target | construct 1 | cells | construct 2 | cells |
|---|---|---|---:|---|---:|
| train | CAST | CAST_P1 | 894 | CAST_P2 | 1,358 |
| train | DAXX | DAXX_P1 | 350 | DAXX_P2 | 982 |
| train | HIRA | HIRA_P1 | 2,888 | HIRA_P2 | 519 |
| train | KDM2B | KDM2B_ENST00000377069.4 | 1,148 | KDM2B_ENST00000377071.4 | 994 |
| train | KLF10 | KLF10_P1 | 449 | KLF10_P2 | 1,483 |
| train | PRCP | PRCP_P1 | 3,415 | PRCP_P2 | 916 |
| train | STAG2 | STAG2_P1 | 662 | STAG2_P2 | 911 |
| train | TMSB4X | TMSB4X_P1 | 2,499 | TMSB4X_P2 | 2,261 |
| validation | ARPC2 | ARPC2_P1 | 491 | ARPC2_P2 | 489 |
| validation | IGF1R | IGF1R_P1 | 1,614 | IGF1R_P2 | 442 |
| validation | SMARCB1 | SMARCB1_P1 | 1,284 | SMARCB1_P2 | 422 |
| validation | ZNF581 | ZNF581_P1 | 1,393 | ZNF581_P2 | 1,109 |
| test | DPF2 | DPF2_P1 | 95 | DPF2_P2 | 1,049 |
| test | GRB10 | GRB10_P1 | 1,938 | GRB10_P2 | 2,136 |
| test | MED15 | MED15_P1 | 1,459 | MED15_P2 | 2,178 |
| test | PCCB | PCCB_P1 | 2,358 | PCCB_P2 | 393 |
| test | SOX4 | SOX4_P1 | 2,785 | SOX4_P2 | 558 |

16 of 17 are a **two-promoter pair (P1 / P2)**; KDM2B is **two transcript ids**. Expected per-`n` coverage
(needs ≥ 2n in the larger construct and ≥ n in the smaller): n = 25 → 17, n = 50 → 17, n = 100 → 16
(DPF2 out), n = 200 → 12 (DPF2, DAXX, KLF10, IGF1R, SMARCB1 out — smaller construct < 200 for the last
three… stated as an expectation, the run reports the realised set).

---

## 2. Known-answer rows (abort if any fails; tolerance stated per row)

* **KA-1** — 273's lane-parity 200/200 median on `data/harness/h1-test.h5ad`: **+0.242920** (exactly
  `0.24291977562969289`), 85 targets, batch parity of the sorted lane index, disjoint control halves
  (9,208 / 9,192 by NTC-construct parity), residual = minus the mean over the 85 targets of that side,
  17,995 scored columns. Required at **0.0e+00** with this experiment's code.
* **KA-2** — 273's 5-target construct row from the same file, 273's construction (all cells of each
  construct, centred on the 85-target FULL-side mean): median **+0.013512255772155394**; per target
  DPF2 −0.03155149708284103, GRB10 +0.013512255772155394, MED15 −0.09916872748565872,
  PCCB +0.027145874176245796, SOX4 +0.35394153008771090; construct cell counts 37/363, 191/209, 150/250,
  351/49, 338/62. Required at **0.0e+00**.
* **KA-3** — integer identity: for every target, (lane-half A counts) + (lane-half B counts) equals the
  full-target count vector exactly, and the cell counts add (`0` absolute difference on int64 sums).
* **KA-4** — the construct parser applied to `adata_Test.h5ad`'s obs reproduces 273's five multi-construct
  targets and their two guide ids exactly, and every `h1-test.h5ad` barcode is found in the test source
  with a matching lane suffix.

If KA-1 or KA-2 misses, the run aborts and RESULTS.md reports the failure instead of the ceiling.

---

## 3. The matched design (FROZEN)

### 3.1 Response object (285 §3 / 273 / 159, verbatim `lib285.py`)

`prof(S) = log1p(50,000 · S / ΣS)` on summed integer counts. `delta = prof(S) − prof(control side)`.
Residual = `delta − m_side`, where **`m_side` is rebuilt for every (split, n, draw, side)** from a
**matched-n sample of every target of that split** (285 §3's first artefact: centring a sub-sampled side on
a full-depth mean inflates the cosine). Per-side centring (side A's pool from side-A cells, side B's from
side-B cells) — a shared centring vector would inject a common component into both residuals.

Scored columns = 18,080 minus that split's own target-gene columns (150 / 50 / 100 removed → 17,930 /
18,030 / 17,980). Controls = the split's 38,176 non-targeting cells, split into **disjoint halves by NTC
construct parity** (273's convention), side A ↔ parity 0.

Two readouts, both from `lib285.views` verbatim:
* **residual cosine** per target (`cos_resid`);
* the scorer's own per-target **`pds`** (218's `_per_target`, == `cell-eval2`'s `pds_cosine` to 5.6e-17 by
  273 KA-2), computed with **the split's other targets as distractors** — side A/B matrices are the
  contrast's rows for the two-construct targets stacked on the matched-n centring-pool rows of every other
  target; the two-construct rows' `pds` is read out. `pds` is the mean of both directions.

### 3.2 The three contrasts, at n ∈ {25, 50, 100, 200} cells per side, 8 draws averaged

For each two-construct target, `C_big` = the construct with more cells, `C_small` = the other. One random
permutation of `C_big`'s cells per draw supplies side A to all three contrasts (285's `k562_guide.py`
design; the tightest paired comparison):

* **(a) cell-sampling reference** — `C_big` ranks [0, n) vs ranks [n, 2n). Requires |C_big| ≥ 2n.
* **(b) construct contrast** — `C_big` ranks [0, n) vs n cells drawn without replacement from `C_small`.
  Requires |C_big| ≥ 2n and |C_small| ≥ n (same cohort as (a), so (a) and (b) are paired per target).
* **(d) lane halves** — n cells from `C_big`'s even-lane cells vs n from its odd-lane cells (parity of the
  sorted 48-lane index, 273's convention). Requires each parity group of `C_big` ≥ n. Reported on its own
  cohort.

Side-identity note: `C_big` is P1 for some targets and P2 for others, but the centring vectors are built
from each target's **pooled** cells and so carry no construct identity — 285 §3's second artefact (a
side whose identity varies makes its own centring a mixture) cannot arise. A **sensitivity row** re-runs
(b) with the fixed assignment construct-1 → side A, construct-2 → side B and must not move the pooled
median by more than the draw-to-draw spread.

### 3.3 The statistics

* **ρ_guide = median(b) / median(a)** at the same n, on the paired cohort. Per target: `b_t / a_t`
  (reported, with the median of ratios alongside the ratio of medians). **Paired target bootstrap, 20,000
  draws** (`lib285.boot_median`, `boot_median_diff`, and the ratio-of-bootstrapped-medians used by
  `k562_guide.py`) for: (a), (b), (d), (a)−(b), and ρ_guide.
* **Spearman–Brown**: `SB(r, k) = k·r / (1 + (k−1)·r)`. Extrapolate (a) from each n to 400/400
  (k = 400/n) and report the implied between-construct cosine at the challenge's budget as
  **ρ_guide × SB(a_n, 400/n)**, together with the directly SB-extrapolated (b).
* **Scorer statistic**: mean and median per-target `pds` for (a) vs (b) at each n, against the 0.5 chance
  level of the rank statistic.
* **On-target knockdown** per construct: `prof(S_c)[own gene] − prof(control)[own gene]` at FULL depth
  (all of that construct's cells) and at n = 100. This column is excluded from every scored statistic.
* **Weaker-construct covariate** (285 §6b's stratification as a per-target covariate): Spearman rank
  correlation across the 17 targets between ρ_guide (and (b)) and (i) the weaker construct's knockdown,
  (ii) the mean of the two knockdowns, (iii) |knockdown difference|, (iv) the target's response magnitude
  (full-depth residual row norm), (v) min(cells) per construct.

### 3.4 What would be reported as the headline, whichever way it comes out

One line: `CHALLENGE-LIBRARY CONSTRUCT CEILING (H1, Flex): ρ_guide <point> [lo, hi] on <n> targets;
between-construct cosine at 400/400 (SB) <x> vs same-construct <y>; scorer pds (a) <p> → (b) <q>;
knockdown explains <r> of it.` **A high ρ_guide (constructs agree) is as publishable as a low one** — it
would say Arc's validated dual-guide library escapes the public-library construct term, and that 273's
n = 5 read was small-sample. No model, no build, no submission either way.

---

## 4. Frozen predictions (stated before any cosine)

Recorded so the result can falsify something. These are **predictions, not gates** — the deliverable is
the measured ceiling with its interval either way.

* **P-1.** ρ_guide at n = 100 will be **below 0.5** — i.e. the construct term is real in Arc's library too.
* **P-2.** ρ_guide will be **above K562-GWPS's 0.09**, because Arc's constructs are knockdown-validated
  dual-guide vectors while GWPS's are single guides of unscreened efficiency. Point guess 0.2–0.5.
* **P-3.** (d) lane halves will read within ±0.03 of (a) at every n (285: batch is free; 273: Spearman
  +0.994 between batch and random halves on this line).
* **P-4.** SOX4 will remain the top target by (b) and MED15 near the bottom (273's ordering), but with
  17 targets the pooled median will differ from 273's +0.0135 by more than its own bootstrap width.
* **P-5.** The scorer's `pds` for (b) will be **closer to chance (0.5)** than for (a) at every n, and the
  gap will be larger than the (a)−(b) cosine gap suggests, because `pds` is a threshold statistic
  (R262 §4.2's knee).
* **P-6.** On-target knockdown of the weaker construct will correlate positively with ρ_guide
  (Spearman > 0), but will not explain most of the between-target spread (|ρ| < 0.6).

---

## 5. Constraints obeyed

Analysis only. Reads: the three 2025 h5ad files (streamed in CSR row blocks, one file at a time),
`data/harness/h1-test.h5ad`, `pert_counts_*.csv`, `experiments/285-reliability-benchmark/work/lib285.py`,
`experiments/273-flex-reference-instrument/{lib273.py,truth273.py}` (read-only, imported/adapted).
Writes: only inside `experiments/293-challenge-construct-ceiling/`. No `vcc`, no `cell-eval2`, no
build, no download, nothing under `submissions/` or `data/`. Peak RSS target < 6 GiB (cells stored as
uint16, capped at 600 per target for the centring/distractor pool and 1,600 per two-construct construct;
the cap is recorded and asserted). Time box 3 h.

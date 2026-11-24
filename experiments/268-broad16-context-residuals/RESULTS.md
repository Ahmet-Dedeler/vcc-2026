# 268 — Multi-context biology on Broad-16: the reproducible context-residual map (stage 1, feasibility)

Date: 2026-09-05. Measurement only; no build, no `vcc`, no `cell-eval2`, nothing written outside this directory
(the Figshare single-cell archive members were fetched into `work/sc/`, never under `data/`).

## 0. Reserved confirmation lines (declared before any residual was computed)

`work/reserve.json` (sha256 `1ce639b509d353708feded808914c63233bcee03513ace1c3f47cb93d071b9ce`, written 16:20Z, before the
first response vector existed). Rule: sort the 16 lines by name; reserve positions 4, 8, 12, 16.

| reserved line | batch | lineage |
|---|---|---|
| `JVE127` | apollo34 | Bowel |
| `MG63` | apollo34 | Bone |
| `SLR23` | apollo32 | Kidney |
| `UMRC3` | apollo35 | Kidney |

Every batch of a reserved line also contains at least one working line (apollo34: ONE58, T3M5; apollo32: SKGII, U343;
apollo35: KMRC20, SNU8), so a stage-2 confirmation on a reserved line is not confounded with an unseen batch. The reserved
lines' single-cell files were NOT downloaded; their pseudobulk columns in the on-disk CSV were not read by any analysis
script (`load268.py` asserts on them). Working lines (12): C4I, HCT15, HKGZCC, KMRC20, KYSE140, KYSE450, ONE58, RPMI7951,
SKGII, SNU8, T3M5, U343.

## 1. Data, experimental units, and the parse (known-answer row)

**What is on disk.** The Broad/DepMap 16-line Perturb-seq release (preprint DOI `10.64898/2026.08.24.746802`, Figshare
article 33273600, CC BY 4.0). `data/external/broad16line/{data,metadata}.zip` hold only gene sets, gene locations and the
line metadata — no responses. The responses used by 111 are its selective extraction
`data/external/depmap26q1/extracted/processed/pseudobulk_sum_by_ko.csv`: **one pseudobulk SUM per (line, KO)**, 34,648
genes × 1,593 columns, 7 missing pairs (C4I/PSMA1, C4I/RRM1, HKGZCC/KIF11, HKGZCC/RRM1, KMRC20/RRM1, U343/KIF11,
U343/RRM1), NaN where a gene is absent from a line's feature set (2,586–5,659 genes per line; 24,979 genes present in all
16 lines). No guide, replicate or cell unit exists in that file. `work/cache_pb.py` reproduces every count in 111's
inventory (16 lines, 6 batches, 11 lineages, 100 targets, 10 olfactory-receptor (OR) controls, 87 complete, 1,593 columns).

**What was fetched to get units (HTTP range requests into `work/`, CRC32-verified, nothing under `data/`):** from
`single_cell_data.zip` (3.17 GB) the 10x `filtered_feature_bc_matrix` h5 and the guide-UMI table (`*_CRISPR.filt`) of
each of the **12 working lines only** (2.3 GB compressed, 8.3 GB on disk; the four reserved lines were not fetched), and from
`processed.zip` the authors' per-cell table `processed/all_cells_table.csv` (375,074 cells: `pass_qc`, `assigned_grna`,
`assigned_ko`, `arm_trunc`, `arm_gain`). Design: 100 Cas9 knockout targets (90 analysis + 10 OR controls), **two guides per
target** (`_v1_`, `_v2_`), one guide per cell, 16 lines in 6 assay batches (apollo31–38), 11 lineages.

**Parse known-answer: the on-disk pseudobulk is reproduced EXACTLY.** `pseudobulk_sum_by_ko[line, KO]` = sum of the h5
counts over cells with `pass_qc` and neither `arm_trunc` nor `arm_gain`, grouped by `assigned_ko`, duplicate gene symbols
summed (the authors' `anndata_groupby(..., sum)`): **1,193 / 1,193 (line, KO) columns identical on the 12 working lines**
(`work/halves/<line>.json`). Rules that do not reproduce it, for the record: `pass_qc` alone (cell counts match the
authors' `n_cells_per_ko` exactly but 2/100 columns identical, median L1 deviation 0.22 — the arm-event cells are in the
QC count but not in the pseudobulk); ratio callers on the raw guide UMIs (median L1 deviation 0.29–0.39); a per-guide
mixture caller (0.29). Every downstream number therefore stands on the authors' own cell assignment.

**Units used below.** *Cell halves*: each (line, KO)'s cells split at random into two halves (seed 268); each half's
response is measured against the SAME half of that line's 10 OR-control KOs, so the two halves share no cell on either
side. *Guide halves*: the `v1` cells vs the `v2` cells, each against the OR controls' `v1`/`v2` cells. *Targets*: 82 of the
87 complete targets have ≥ 6 cells on every working line (both halves ≥ 3); the five excluded are AURKB, MAD2L1, MDM2,
TOP2A, UBE2D3 — the most lethal knockouts, so the strongest-response tail is under-represented. *Gene axis*: genes present
in all 12 working lines, the 100 target genes removed, mean control CPM ≥ 5 over the lines used for the fit → **12,199
genes** (7,628 at CPM ≥ 20, the robustness run). *Response*: `log1p(5e4·p/Σp)` (the scorer's pseudobulk frame) minus the
same-half control log-profile. *Noise floor*: the cross-half method below needs no model — the noise energy is measured
as the part of a half's energy that the other half does not reproduce; the OR-control KOs, passed through the identical
pipeline against the other nine controls, are the empirical null in the same units.

| line | batch | lineage | h5 cells | cells kept | KOs | cells/KO median | min | max | control cells (10 OR KOs) | UMI/cell | targets < 6 cells | KA exact |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| C4I | apollo33 | Cervix | 20682 | 5442 | 98 | 55 | 1 | 104 | 652 | 15507 | 3 | 98/98 |
| HCT15 | apollo31 | Bowel | 14305 | 4270 | 100 | 43 | 3 | 76 | 495 | 50005 | 2 | 100/100 |
| HKGZCC | apollo38 | Biliary Tract | 40262 | 9911 | 98 | 96 | 2 | 224 | 1052 | 31881 | 3 | 98/98 |
| KMRC20 | apollo35 | Kidney | 28936 | 8771 | 99 | 90 | 2 | 185 | 1004 | 12611 | 3 | 99/99 |
| KYSE140 | apollo31 | Esophagus/Stomach | 13463 | 6278 | 100 | 62 | 6 | 118 | 727 | 46714 | 0 | 100/100 |
| KYSE450 | apollo31 | Esophagus/Stomach | 24964 | 12833 | 100 | 129 | 13 | 232 | 1462 | 37502 | 0 | 100/100 |
| ONE58 | apollo34 | Pleura | 25284 | 11535 | 100 | 119 | 2 | 244 | 1244 | 20013 | 1 | 100/100 |
| RPMI7951 | apollo31 | Skin | 11727 | 5658 | 100 | 58 | 1 | 97 | 650 | 51585 | 2 | 100/100 |
| SKGII | apollo32 | Cervix | 11074 | 5208 | 100 | 52 | 2 | 112 | 576 | 76867 | 5 | 100/100 |
| SNU8 | apollo35 | Ovary/Fallopian Tube | 37795 | 8753 | 100 | 87 | 4 | 217 | 1014 | 11723 | 3 | 100/100 |
| T3M5 | apollo34 | Head and Neck | 18586 | 4827 | 100 | 47 | 5 | 84 | 552 | 30006 | 1 | 100/100 |
| U343 | apollo32 | CNS/Brain | 8613 | 4122 | 98 | 42 | 5 | 83 | 468 | 28572 | 2 | 98/98 |

87,608 cells kept in total. Depth differs 6.6× across lines (11.7k–76.9k UMI per cell), and cells per KO 3×: the noise
floor is line-specific, which is exactly why it is measured per (line, target) rather than assumed.

## 2. Variance decomposition — the core table

**Estimator (`work/decomp268.py`).** For target *t* with half-responses `R_A[t,l]`, `R_B[t,l]` over the 12 lines: shared
`m = mean_l R`, residual `r = R − m`. The total energy of a half is `E = Σ_l ‖R‖²` (averaged over the two halves); the
cross-half products are unbiased for the true energies because the halves' noises are independent:
`X_shared = 12·⟨m_A, m_B⟩`, `X_line = Σ_l ⟨r_A[l], r_B[l]⟩`, `noise = E − X_shared − X_line` — and the identity
`Σ_l⟨R_A,R_B⟩ = X_shared + X_line` holds exactly. The half-depth noise is converted to the full pseudobulk's depth by
halving it (checked: the measured `(E_full − X)/(E_half − X)` is **0.4996** per target, median; expected 0.5). The
line-specific part is further split into the component parallel to the target's shared vector (amplitude modulation) and
the orthogonal part (a different direction).

**Known-answer rows.** (i) Synthetic tensors with known energies (87 × 12 × 2,000): truth 0.60/0.25/0.15 → estimate
0.599/0.251/0.150; truth 0.20/0.10/0.70 → 0.200/0.100/0.700; pure noise → 0.000/−0.000/1.000 with 6.5 % of (t,l)
residual cosines above the permutation null's 95th percentile (5 % expected). (ii) The OR-control null (below).

| construction (pooled over targets, energy-weighted, at FULL pseudobulk depth) | shared | line-specific & reproducible | of which direction-changing | of which amplitude only | noise |
|---|---:|---:|---:|---:|---:|
| **cell halves, 82 targets, 12,199 genes (CPM ≥ 5)** | **0.091** | **0.115** | 0.099 | 0.019 | **0.794** |
| same, per-target median (IQR of the line-specific fraction 0.040–0.138) | 0.031 | 0.080 | 0.074 | — | 0.880 |
| same at the measured HALF depth (what one split half sees) | 0.051 | 0.064 | — | — | 0.885 |
| guide halves (v1 vs v2), 78 targets with both guides on all 12 lines | 0.080 | 0.091 | 0.076 | 0.013 | 0.829 |
| **OR-control null** (10 control KOs, each vs the other nine, same pipeline) | 0.019 | **0.028** | 0.027 | 0.003 | 0.953 |
| robustness: CPM ≥ 20 (7,628 genes) | 0.109 | 0.129 | 0.107 | 0.021 | 0.762 |
| robustness: CPM ≥ 20, OR-control null | 0.024 | 0.033 | — | — | 0.944 |

Readings, all measured:

- **A reproducible line-specific component exists and is as large as the shared component**: 11.5 % of the full-depth
  response energy pooled (median 8.0 % per target) against 9.1 % shared (median 3.1 %); it exceeds the shared part on
  63 of 82 targets. It is **four times the OR-control null** (2.8 %), which bounds the target-independent artefacts a
  KO-vs-controls comparison carries in this design (guide-population differences, the cut itself, line-level QC).
- **It is a different direction, not a different dose**: 86 % of the reproducible line-specific energy is orthogonal to
  the target's shared response vector (0.099 of 0.115); amplitude modulation is 2 % of response energy.
- **It survives the guide unit**: across guides (each guide its own cells and controls) the line-specific fraction is
  0.091 vs 0.115 across cell halves, i.e. 79 % of it is not guide-specific (Spearman across targets 0.87); the shared part
  keeps 88 %.
- **Per (target, line) it is weak and diffuse**: the split-half cosine of the residual has median 0.046 (permutation null
  95th percentile 0.096); 26.2 % of the 984 (t,l) pairs exceed the null's 95th percentile (OR controls: 7.5 %; direction-
  only part 23.2 % vs 6.7 %); 72/82 targets have at least one such line and 34/82 at least three (chance ≈ 0.6 per target).
  The target-independent line offset (`mean_t r[t,l]`) is 10 % of the line-specific energy.
- **Strong responders carry both parts**: Spearman(shared, line-specific) = +0.71 across targets. The top of the table is
  MTOR, MYC, BRD4, SMG6, XRN1, HDAC3, SLC39A10 (shared 0.16–0.37, line-specific 0.14–0.27, noise < 0.65); the bottom
  (ADAT2, SLC20A1, TWNK, STEAP3, UBE2C, N6AMT1, TXN) is noise to within ±0.02 on both. Full table: appendix.
- **Batch structure is present**: the cross-half residual similarity between lines of the same batch is +0.020 vs −0.106
  between batches (9 within-batch pairs; permutation of batch labels over lines, p = 0.001; CPM ≥ 20: p = 0.0015). Split
  halves within one experiment reproduce batch effects as faithfully as biology, so part of the 11.5 % is line-level
  technical structure that the OR null (target-independent by construction) does not bound. Biology and batch are only
  separable with a line measured twice, which this set does not contain.

## 3. Is the line-specific part programme-level? Mostly no

SVD of the residual unfolded to (82·12) × 12,199, fitted on half A; each component's reproducible energy is its cross-half
energy `⟨r_A v_k, r_B v_k⟩`, which sums exactly to `X_line`; the null is the same product with the rows of half B
permuted (200 draws).

| matrix | components for 50 / 80 / 90 % of the reproducible energy | top 1 / 2 / 5 / 10 / 20 cumulative share | components above the permutation p99 (first 60) |
|---|---|---|---|
| **residual (line-specific)** | **27 / 167 / 370** | 0.094 / 0.229 / 0.311 / 0.364 / 0.449 | 58 of 60 |
| residual with the line offset removed | 27 / 189 / 391 | 0.103 / 0.231 / 0.293 / 0.363 / 0.454 | 60 of 60 |
| **shared response (82 × G)** | **3 / 10 / 19** | 0.338 / 0.430 / 0.686 / 0.814 / 0.913 | 58 (reference: low-rank) |
| OR-control residual (noise null) | 35 / 74 / 90 (of a total that is itself ≈ 0) | 0.204 / 0.248 / — | 10 |
| robustness CPM ≥ 20, residual | 22 / 138 / 327 | 0.112 / 0.264 / — | 59 |

**Effective rank**: the reproducible line-specific residual needs **27 components for half of its energy** where the shared
response needs 3, and 58 of the first 60 components carry reproducible energy — it is a genuine, high-rank object, not a
few programmes. Its two largest components (23 % together) are, in words (gene loadings and Hallmark enrichment of the
component; `work/decomp268_cpm5.json`):

- **C2 (13.5 %, cross-half score correlation 0.40)** — E2F targets, G2M checkpoint, MYC targets (MKI67, CENPF, TPX2, RRM2,
  CCNB1) against NEAT1, GDF15, CDKN1A, SQSTM1, FBXO32: the **proliferation-versus-p53/stress-arrest axis**, i.e. how far a
  line arrests. Strongest in HKGZCC, KMRC20, T3M5 and on the strongest targets (POLR2D, MYC, SLC39A10, BRD4, USP37).
- **C1 (9.4 %, correlation 0.25)** — G2M/mitotic-spindle/E2F genes up against oxidative phosphorylation, ribosomal and
  S100/LGALS genes down (KMRC20, HCT15, HKGZCC; POLR2D, USP37, MYC, TRAIP, PKMYT1): a mitotic-accumulation versus
  metabolic/differentiation axis.
- C3 (2.2 %) is a **line offset** (64 % of its score variance is between lines: C4I and U343): keratinocyte/S100A7/PI3
  genes against IFIT2/ATF3/GADD45A (interferon, TNFα-NFκB). C4 (2.7 %): mitotic (AURKA, CENPF, BUB1) vs OXPHOS. C5 (3.3 %):
  TNFα-NFκB, p53 pathway, apoptosis, hypoxia (GADD45A, DDIT3, PPP1R15A, CXCL3) against DNA replication (MCM6, CCNE2, POLA1)
  in U343, HKGZCC, SNU8. C6 (3.3 %, 53 % between-line): GADD45A/B, SLC7A11, DDIT3, IRF1 stress genes against mitotic genes
  in T3M5, U343, C4I.

So the only programme-level structure in the line-specific part is **the line-dependent balance of proliferation arrest
versus p53/stress response inside the programme that also dominates the shared response** (the shared response's first
component is 34 % of its energy and is the same E2F/G2M axis). Everything else is spread over dozens of weak, reproducible
components — target-specific line-specific biology, or batch, at a rank this data cannot resolve into programmes.

## 4. Does it relate to basal state? Leave-one-line-out, with the three controls

**Protocol (`work/lolo268.py`, `work/lolosum268.py`).** 12 folds. Inside each fold everything is fitted on the 11 training
lines: the gene filter, the centring of the control profiles, the basal programmes (PCA of the 11 control log-profiles;
`K = 2` components hold 38 % of their variance, `K = 3` 50 %), the residual's programme basis (SVD of the training
interaction residual, `J = 5`). The held line's truth is its FULL-data residual from the training mean; its reproducible
part is read from the held line's own two cell halves (each against the training halves' mean), whose mutual cosine is
the **ceiling** (median 0.045, mean 0.061 — the per-(t,l) residual is barely reproducible at this depth, §2).

Arms, each predicting the held residual `r̂[t]` (few parameters: ridge on standardised basal scores, λ = 1):
(a) **target-only** = the shared response, `r̂ = 0`, the strong baseline; (b) **recipient-only** = a target-independent line
offset regressed on basal scores; (c1) **target × basal, per-gene ridge** (`K` coefficient vectors per target); (c2) **target ×
basal, programme-level** (`K` coefficients per target and residual programme, `J = 5`); (nn) the residual of the **nearest
basal line**; and the **oracle line offset** (the held line's true mean residual — uses held data; reference only).
Controls: **WRONG recipient** — the same fitted model driven by each other training line's basal scores (or residual), mean
over the 11; **SHUFFLED labels** — basal scores permuted across the training lines before fitting, predicted with the
correct held scores (5 permutations). Three readouts per target: cosine of `r̂` to the held residual; the change in cosine
of the FULL prediction `m + r̂` to the held response over target-only (mean target-only cosine 0.113); the change in a
`pds`-like within-line retrieval (rank of the held truth among the 82 by cosine; target-only 0.766). Paired bootstrap
over targets (each target carries its 12 folds; 4,000 draws).

| arm (K = 2) | cos(pred. residual, held residual) [95 % CI] | folds > 0 | Δ cos(full response) [95 % CI] | folds > 0 | Δ retrieval [95 % CI] | folds > 0 |
|---|---|---|---|---|---|---|
| (a) target-only | 0 | — | 0 | — | 0 | — |
| (b) recipient-only offset | +0.0052 [+0.0026, +0.0081] | 7/12 | −0.0011 [−0.0015, −0.0007] | 3/12 | −0.0046 [−0.0078, −0.0015] | 3/12 |
| **(c1) target × basal, ridge** | **+0.0181 [+0.0134, +0.0229]** | 9/12 | **−0.0071 [−0.0093, −0.0050]** | 1/12 | **+0.0003 [−0.0098, +0.0105]** | 6/12 |
| **(c2) target × basal, programmes** | **+0.0174 [+0.0110, +0.0240]** | 8/12 | **+0.0000 [−0.0017, +0.0017]** | 5/12 | **−0.0010 [−0.0109, +0.0097]** | 6/12 |
| (nn) nearest basal line's residual | +0.0353 [+0.0306, +0.0401] | 8/12 | −0.0494 [−0.0603, −0.0387] | 0/12 | −0.0596 [−0.0810, −0.0387] | 0/12 |
| oracle line offset (held data) | +0.2968 [+0.2909, +0.3026] | 12/12 | +0.1690 [+0.1572, +0.1808] | 12/12 | −0.0192 [−0.0392, +0.0025] | 4/12 |
| (b) WRONG recipient | +0.0003 [+0.0002, +0.0005] | 7/12 | −0.0042 [−0.0047, −0.0036] | 0/12 | −0.0061 [−0.0087, −0.0034] | 1/12 |
| (b) SHUFFLED labels | −0.0006 [−0.0015, +0.0003] | 4/12 | −0.0010 [−0.0013, −0.0007] | 1/12 | −0.0007 [−0.0020, +0.0004] | 5/12 |
| (c1) WRONG recipient | +0.0012 [+0.0006, +0.0018] | 9/12 | −0.0295 [−0.0343, −0.0248] | 0/12 | −0.0355 [−0.0459, −0.0249] | 0/12 |
| (c1) SHUFFLED labels | −0.0073 [−0.0091, −0.0055] | 2/12 | −0.0124 [−0.0144, −0.0105] | 0/12 | −0.0104 [−0.0159, −0.0048] | 1/12 |
| (c2) WRONG recipient | +0.0007 [+0.0003, +0.0012] | 9/12 | −0.0083 [−0.0101, −0.0068] | 0/12 | −0.0216 [−0.0297, −0.0131] | 0/12 |
| (c2) SHUFFLED labels | −0.0007 [−0.0036, +0.0022] | 5/12 | −0.0024 [−0.0031, −0.0017] | 1/12 | −0.0053 [−0.0095, −0.0012] | 2/12 |
| (nn) WRONG recipient (a random other line's residual) | +0.0112 [+0.0098, +0.0130] | 12/12 | −0.0624 [−0.0748, −0.0507] | 0/12 | −0.0872 [−0.1043, −0.0697] | 0/12 |

Paired differences that decide the question (per-target means, bootstrap over targets; `P` = P(Δ ≤ 0)):

| comparison | residual cosine Δ | Δ cos(full) | Δ retrieval |
|---|---|---|---|
| (c1) − (c1) WRONG recipient | +0.0169 [+0.0121, +0.0218], P 0.000 | +0.0224 [+0.0186, +0.0263], P 0.000 | +0.0358 [+0.0249, +0.0467], P 0.000 |
| (c1) − (c1) SHUFFLED labels | +0.0254 [+0.0201, +0.0309], P 0.000 | +0.0053 [+0.0033, +0.0072], P 0.000 | +0.0107 [+0.0007, +0.0211], P 0.018 |
| **(c1) − (a) target-only** | +0.0181 [+0.0135, +0.0231], P 0.000 | **−0.0071 [−0.0093, −0.0050], P 1.000** | **+0.0003 [−0.0097, +0.0110], P 0.48** |
| (c2) − (c2) WRONG recipient | +0.0167 [+0.0103, +0.0234], P 0.000 | +0.0083 [+0.0061, +0.0107], P 0.000 | +0.0207 [+0.0110, +0.0309], P 0.000 |
| (c2) − (c2) SHUFFLED labels | +0.0182 [+0.0104, +0.0261], P 0.000 | +0.0024 [+0.0006, +0.0043], P 0.005 | +0.0044 [−0.0053, +0.0151], P 0.20 |
| **(c2) − (a) target-only** | +0.0174 [+0.0108, +0.0238], P 0.000 | **+0.0000 [−0.0016, +0.0017], P 0.51** | **−0.0010 [−0.0109, +0.0099], P 0.58** |
| (b) − (b) WRONG recipient | +0.0049 [+0.0022, +0.0077], P 0.000 | +0.0030 [+0.0024, +0.0036], P 0.000 | +0.0014 [−0.0020, +0.0047], P 0.21 |
| (nn) − (a) target-only | +0.0353, P 0.000 | −0.0494 [−0.0606, −0.0389], P 1.000 | −0.0596 [−0.0809, −0.0376], P 1.000 |

`K = 3` changes nothing that matters (c1 +0.0204 / −0.0082 / −0.0004; c2 +0.0183 / −0.0001 / −0.0019; `work/lolo268_K3.json`).

**Reading, plainly.** The basal-conditioned models learn something real: their predicted residual has cosine +0.017 to
+0.018 with the held line's residual, against +0.001 for the same model handed the wrong recipient's basal state and
−0.007/−0.001 for shuffled line labels (P ≤ 0.000 on both controls, 8–9 of 12 folds positive), i.e. about 30 % of the
split-half ceiling (0.045–0.061). **But it is not a usable prediction**: adding it to the shared response does not improve
the prediction of the held line's response (c1 −0.007, c2 0.000 in cosine) nor within-line retrieval (+0.0003 and −0.001,
both ±0.010), so at this construction **(c) does not beat (a)**; the recipient-only offset and the nearest-basal-line rule
are net harmful; and even the oracle line offset — the held line's true mean residual, unavailable to any model — buys
+0.169 in cosine but nothing in retrieval (−0.019), because a target-independent offset is the one component a
`pds`-type metric cannot reward. Context learning of the residual DIRECTION is demonstrated at a tiny effect size;
context learning that improves a transfer prediction is **not demonstrated on this data**.

## 5. The bridge to the panel

**Overlap of the 90 Broad targets with the 300-target VCC panel: 0** (`data/challenge/pert_counts.csv`, verified;
111 recorded the same). No per-target comparison with 262's anti-aligned group is possible; the comparison the mandate
asked for is not defined on this data. The only bridge is programme-level: the dominant reproducible line-specific
axis here (proliferation/E2F–G2M versus p53/stress arrest, §3) is a line-dependent balance inside the most shared
programme, which is the kind of component our whitening operator suppresses as a hub rather than predicts.

## 6. Verdict and what stage 2 would be

**REPRODUCIBLE LINE-SPECIFIC COMPONENT: YES (fraction 0.115 of full-depth response energy pooled, median 0.080 per target;
4× the OR-control null's 0.028; 86 % direction-changing; 79 % of it survives across guides; part of it is batch, p = 0.001);
PROGRAMME-LEVEL: NO (rank 27 components for half of its reproducible energy vs 3 for the shared response; the two largest
components, 23 %, are the proliferation-vs-p53-arrest balance of the shared programme); BASAL-STATE-PREDICTABLE under
LOLO: NO (detectable — residual cosine +0.018 vs +0.001 wrong recipient — but zero gain on the full response and on
retrieval; (c) does not beat (a)).**

Stage 2 is **not earned at this construction**, and the closure is at its scope: linear conditioning of a target's residual
on 2–3 basal principal programmes (per-gene ridge or a rank-5 residual basis), fitted on 11 lines, does not produce a
transfer prediction better than the shared response on 12 cancer lines × 82 Cas9 knockouts. Three things this does NOT
close: (i) a targeted hypothesis on the one programme-level component that exists — that a line's position on the
arrest-versus-p53 axis (residual C2) follows its TP53 status / basal E2F–p53 activity (the archive carries
`processed/p53_cell_line_properties.csv`); a 1–2-parameter test, cheap, with the four reserved lines as the preregistered
confirmation; (ii) more contexts — 11 training lines fit 2–3 basal parameters, and rule 7 says buy independent contexts
before capacity; (iii) 224's per-gene reliability gain, which remains the baseline any stage-2 model must beat. The
reserved lines were not accessed; no confirmation entry is needed for this stage (`docs/CONFIRMATION-LEDGER.md` should
record the four names and the sha of `work/reserve.json`).

## 7. Honest limits

1. Split halves inside one experiment reproduce batch and line-level technical effects as well as biology; the batch
   test (§2) shows such structure exists. The OR-control null bounds only the target-independent part (2.8 %). A biology-
   only number needs a line measured in two experiments, which this set lacks.
2. Cells per KO are small (median 42–129 per line, halves of 21–65), so the per-(t,l) residual reproducibility is low
   (median cosine 0.046) and per-target fractions carry ±0.03–0.05 of sampling error; the pooled, energy-weighted numbers
   are the reliable ones, and they are dominated by the strong responders (Spearman shared/line +0.71).
3. Cas9 knockouts with the authors' arm-truncation/gain cells excluded, 12 cancer lines of 11 lineages, 82 targets with
   zero panel overlap, a 12,199-gene axis; the five most lethal targets were excluded by the 6-cell rule. The CPM ≥ 20 run
   moves every fraction by ≤ 0.02 in the same direction.
4. "Direction-changing" is orthogonality to the target's own shared vector, not to a programme; a line that runs the
   shared programme with a different gene composition counts as direction-changing.
5. Peak RSS was 6.7 GB for ~30 s on the largest line's build (KYSE450); every analysis step stayed under 2 GB.
6. The p53 hypothesis of §6 was not tested; the LOLO used only generic basal principal programmes.

## 8. Files, provenance, resources

- Reserve: `work/reserve.json` (sha256 `1ce639b5…`). Inventory and parse: `work/cache_pb.py` → `work/pb_cache.npz`,
  `work/inventory.json`; `work/nanaudit.py`.
- Remote-archive tooling (range requests, CRC-verified): `work/zipdir.py`, `work/fetchzip.py`; fetched members in
  `work/sc/` (12 working lines only) and `work/proc/all_cells_table.csv`; the authors' code consulted:
  `work/authors_data_utils.py`, `work/authors_constants.py` (GitHub `broadinstitute/perturb-seq-depmap-public`).
- Build: `work/build_halves2.py` (rule `EXCLUDE_ARM=1`), driver `work/run_halves2.sh`, log `work/halves2.log`, outputs
  `work/halves/<line>.{npz,json}` (per-KO sums: full, halves A/B, guides v1/v2; cell counts; per-column KA). Earlier
  non-reproducing rule tests are recorded in `work/build_halves.py` and `work/assign_rules.py` (kept for the record).
- Analysis: `work/load268.py`, `work/decomp268.py` (→ `work/decomp268_cpm{5,20}.json`, `work/per_target_decomp_cpm{5,20}.csv`,
  `work/decomp268_arrays_cpm{5,20}.npz`, logs `work/decomp{5,20}.log`), `work/lolo268.py` + `work/lolosum268.py`
  (→ `work/lolo268_K{2,3}.{csv,json}`, logs `work/lolo_K{2,3}.log`), `work/units268.py` (→ `work/units268.md`),
  `work/appendix_table.md`.
- Wall clock ≈ 3 h 10 min; download 2.3 GB at 20–36 MB/s; disk 8.3 GB in `work/sc/`; no `cell-eval2`, no build, no `vcc`.

## Appendix — per-target decomposition (82 usable targets, cell halves, CPM ≥ 5, fractions at full depth, sorted by the reproducible line-specific fraction)

| target | shared | line-specific (reproducible) | of which direction-changing | noise | line-specific, guide halves | lines with reproducible residual (of 12) | lines with reproducible direction residual | split-half cosine of the response (median over lines) | min cells |
|---|---|---|---|---|---|---|---|---|---|
| SLC39A10 | 0.164 | 0.272 | 0.165 | 0.563 | 0.199 | 10 | 6 | 0.210 | 19 |
| MTOR | 0.278 | 0.240 | 0.217 | 0.482 | 0.222 | 11 | 10 | 0.368 | 31 |
| UBE3A | 0.052 | 0.238 | 0.122 | 0.710 | 0.247 | 3 | 3 | 0.043 | 21 |
| MYC | 0.277 | 0.234 | 0.165 | 0.489 | 0.225 | 11 | 9 | 0.271 | 13 |
| USP37 | 0.059 | 0.226 | 0.108 | 0.715 | n/a | 7 | 6 | 0.118 | 6 |
| BRD4 | 0.235 | 0.226 | 0.217 | 0.539 | n/a | 11 | 11 | 0.326 | 18 |
| HDAC3 | 0.165 | 0.215 | 0.201 | 0.621 | 0.185 | 10 | 9 | 0.230 | 26 |
| SMG6 | 0.343 | 0.192 | 0.171 | 0.466 | 0.184 | 9 | 7 | 0.343 | 21 |
| AP2M1 | 0.060 | 0.184 | 0.182 | 0.756 | 0.155 | 9 | 9 | 0.146 | 30 |
| TFRC | 0.056 | 0.176 | 0.155 | 0.768 | 0.181 | 4 | 4 | 0.073 | 16 |
| RAB7A | 0.099 | 0.166 | 0.144 | 0.735 | 0.130 | 5 | 5 | 0.101 | 45 |
| RHOA | 0.087 | 0.155 | 0.149 | 0.758 | 0.188 | 6 | 6 | 0.103 | 29 |
| ATP5PO | 0.095 | 0.154 | 0.125 | 0.751 | 0.130 | 4 | 2 | 0.074 | 34 |
| WDR89 | 0.016 | 0.153 | 0.139 | 0.831 | n/a | 2 | 2 | 0.042 | 24 |
| SEC23IP | 0.064 | 0.147 | 0.121 | 0.789 | 0.126 | 6 | 5 | 0.068 | 48 |
| STXBP3 | 0.035 | 0.142 | 0.109 | 0.822 | 0.084 | 2 | 2 | 0.040 | 38 |
| XRN1 | 0.372 | 0.142 | 0.142 | 0.486 | 0.115 | 9 | 9 | 0.358 | 43 |
| CAB39 | 0.100 | 0.142 | 0.146 | 0.759 | 0.172 | 5 | 5 | 0.121 | 42 |
| MAPK1 | 0.049 | 0.140 | 0.131 | 0.810 | 0.084 | 3 | 3 | 0.073 | 38 |
| ADAR | 0.131 | 0.139 | 0.109 | 0.730 | 0.105 | 5 | 2 | 0.107 | 12 |
| IER3IP1 | 0.021 | 0.139 | 0.094 | 0.840 | 0.012 | 4 | 3 | 0.012 | 43 |
| POLR2D | 0.169 | 0.137 | 0.131 | 0.694 | 0.110 | 8 | 6 | 0.188 | 8 |
| KCTD5 | 0.054 | 0.131 | 0.126 | 0.815 | 0.109 | 5 | 4 | 0.084 | 53 |
| PPP2R2A | 0.072 | 0.130 | 0.130 | 0.798 | 0.146 | 4 | 4 | 0.083 | 68 |
| RNF31 | 0.033 | 0.126 | 0.108 | 0.841 | 0.131 | 2 | 2 | 0.064 | 35 |
| CCNE1 | 0.039 | 0.123 | 0.062 | 0.838 | 0.102 | 1 | 1 | 0.038 | 26 |
| MIB1 | 0.021 | 0.122 | 0.116 | 0.857 | 0.093 | 5 | 4 | 0.081 | 41 |
| NRBP1 | 0.145 | 0.119 | 0.131 | 0.735 | 0.106 | 7 | 7 | 0.213 | 19 |
| PRKRA | 0.024 | 0.116 | 0.103 | 0.859 | 0.096 | 4 | 3 | 0.046 | 20 |
| GET4 | 0.027 | 0.115 | 0.091 | 0.857 | 0.120 | 2 | 2 | 0.057 | 33 |
| C3orf38 | 0.037 | 0.114 | 0.085 | 0.849 | 0.074 | 2 | 2 | 0.052 | 21 |
| NFKB1 | 0.031 | 0.109 | 0.106 | 0.859 | 0.106 | 3 | 3 | 0.066 | 49 |
| MIPEP | 0.023 | 0.103 | 0.100 | 0.874 | 0.116 | 2 | 2 | 0.030 | 60 |
| SUZ12 | 0.066 | 0.100 | 0.112 | 0.834 | 0.053 | 3 | 3 | 0.112 | 25 |
| CLPB | 0.044 | 0.099 | 0.091 | 0.857 | 0.039 | 4 | 2 | 0.079 | 31 |
| TRAF2 | 0.030 | 0.093 | 0.090 | 0.876 | 0.118 | 5 | 5 | 0.061 | 32 |
| METTL1 | 0.018 | 0.091 | 0.081 | 0.891 | 0.090 | 2 | 2 | 0.046 | 44 |
| PPIL1 | 0.029 | 0.088 | 0.082 | 0.884 | 0.064 | 2 | 2 | 0.055 | 29 |
| YTHDF2 | 0.080 | 0.087 | 0.093 | 0.833 | 0.108 | 3 | 3 | 0.065 | 38 |
| DCLRE1B | 0.018 | 0.084 | 0.075 | 0.898 | 0.055 | 1 | 1 | 0.035 | 34 |
| CDK4 | 0.004 | 0.083 | 0.088 | 0.913 | 0.093 | 2 | 2 | 0.041 | 43 |
| UBE2S | 0.036 | 0.078 | 0.067 | 0.886 | n/a | 2 | 1 | 0.039 | 35 |
| UBA5 | 0.065 | 0.076 | 0.080 | 0.859 | 0.065 | 2 | 2 | 0.088 | 59 |
| PIK3CA | 0.044 | 0.074 | 0.063 | 0.882 | 0.064 | 3 | 1 | 0.042 | 35 |
| SMCHD1 | 0.009 | 0.073 | 0.072 | 0.918 | 0.066 | 1 | 1 | 0.053 | 50 |
| SDHB | 0.027 | 0.072 | 0.074 | 0.901 | 0.020 | 2 | 2 | 0.045 | 49 |
| TRAIP | 0.050 | 0.067 | 0.062 | 0.883 | 0.081 | 3 | 3 | 0.046 | 12 |
| PKMYT1 | 0.065 | 0.066 | 0.067 | 0.870 | 0.051 | 1 | 0 | 0.088 | 15 |
| METAP1 | 0.038 | 0.062 | 0.069 | 0.899 | 0.053 | 1 | 1 | 0.054 | 39 |
| RAB18 | 0.025 | 0.062 | 0.061 | 0.912 | 0.031 | 1 | 1 | 0.039 | 41 |
| CCDC51 | 0.063 | 0.059 | 0.056 | 0.878 | 0.057 | 3 | 3 | 0.035 | 20 |
| UFL1 | 0.024 | 0.059 | 0.056 | 0.917 | 0.017 | 3 | 3 | 0.026 | 42 |
| WDR48 | 0.031 | 0.058 | 0.057 | 0.910 | 0.048 | 0 | 0 | 0.065 | 35 |
| SLC39A9 | 0.003 | 0.056 | 0.060 | 0.940 | 0.105 | 3 | 3 | 0.053 | 62 |
| MTPAP | 0.169 | 0.055 | 0.053 | 0.776 | -0.002 | 2 | 1 | 0.110 | 47 |
| TYMS | 0.032 | 0.051 | 0.044 | 0.917 | 0.055 | 2 | 3 | 0.057 | 48 |
| PITRM1 | 0.020 | 0.049 | 0.048 | 0.931 | 0.021 | 2 | 2 | 0.027 | 42 |
| LRP8 | 0.014 | 0.045 | 0.043 | 0.941 | -0.000 | 1 | 1 | 0.032 | 38 |
| GPX4 | 0.028 | 0.045 | 0.035 | 0.928 | 0.032 | 1 | 1 | 0.025 | 19 |
| SCYL1 | 0.024 | 0.043 | 0.049 | 0.933 | 0.083 | 4 | 3 | 0.044 | 27 |
| NHLRC2 | 0.002 | 0.042 | 0.040 | 0.957 | 0.037 | 0 | 0 | 0.021 | 51 |
| ATF4 | 0.028 | 0.040 | 0.038 | 0.932 | 0.025 | 1 | 1 | 0.030 | 30 |
| RAD51 | 0.053 | 0.039 | 0.042 | 0.908 | 0.029 | 4 | 3 | 0.058 | 6 |
| DNAJC13 | 0.029 | 0.036 | 0.048 | 0.934 | 0.012 | 1 | 1 | 0.033 | 41 |
| IBA57 | 0.009 | 0.033 | 0.032 | 0.958 | 0.033 | 1 | 1 | 0.018 | 38 |
| CD46 | 0.008 | 0.033 | 0.026 | 0.958 | 0.018 | 0 | 0 | 0.029 | 42 |
| EEFSEC | 0.005 | 0.033 | 0.033 | 0.963 | -0.026 | 0 | 0 | 0.043 | 39 |
| BCL2L2 | 0.012 | 0.031 | 0.027 | 0.956 | 0.006 | 0 | 0 | 0.024 | 43 |
| UBE2H | 0.047 | 0.030 | 0.048 | 0.924 | 0.049 | 2 | 2 | 0.070 | 27 |
| SLC25A3 | 0.019 | 0.029 | 0.031 | 0.951 | 0.001 | 1 | 1 | 0.012 | 44 |
| KIF2C | 0.015 | 0.025 | 0.020 | 0.960 | 0.023 | 1 | 1 | 0.008 | 37 |
| COX5B | 0.022 | 0.024 | 0.031 | 0.954 | 0.025 | 0 | 0 | 0.033 | 40 |
| FBXO7 | 0.014 | 0.024 | 0.015 | 0.962 | 0.053 | 1 | 1 | 0.035 | 45 |
| UBE2G2 | 0.023 | 0.023 | 0.024 | 0.954 | 0.032 | 0 | 0 | 0.025 | 42 |
| SPRING1 | 0.017 | 0.022 | 0.033 | 0.961 | 0.022 | 1 | 1 | 0.032 | 50 |
| TXN | 0.020 | 0.018 | 0.022 | 0.961 | 0.033 | 1 | 1 | 0.025 | 43 |
| N6AMT1 | 0.017 | 0.013 | 0.015 | 0.970 | 0.013 | 0 | 1 | 0.025 | 29 |
| UBE2C | 0.024 | 0.006 | 0.004 | 0.969 | 0.024 | 2 | 1 | 0.004 | 37 |
| STEAP3 | 0.008 | 0.004 | 0.001 | 0.987 | -0.002 | 0 | 0 | 0.021 | 33 |
| TWNK | 0.015 | -0.012 | -0.002 | 0.997 | -0.060 | 1 | 1 | 0.020 | 40 |
| SLC20A1 | 0.002 | -0.014 | -0.016 | 1.012 | -0.040 | 0 | 0 | -0.008 | 46 |
| ADAT2 | 0.001 | -0.018 | -0.011 | 1.017 | -0.025 | 1 | 2 | -0.007 | 71 |

# 289 — PREREGISTRATION

Written 2026-09-06, **before any state-feature LOLO number, before any P1/P2 statistic, before any gating number,
and before any per-line state-feature value was computed or inspected.** The only number that exists at the time of
writing is the known-answer row KA-1 (§7), which is a bit-identical re-run of 268's published arm and contains no new
information.

Frozen by sha256 into `work/PREREGISTRATION.sha256`. Anything added after the hash is labelled **POST-HOC** in
`RESULTS.md`.

---

## 1. Question and unit

**Question.** 268 reported that the reproducible line-specific component of a Broad-16 perturbation response is *not
basal-state-predictable*: a ridge on **two basal expression PCs** predicted the held line's residual with cosine
+0.0181 but gained nothing on the full response (Δcos −0.0071 per gene, +0.0000 programme-level; R268 §4). The source
paper (Maffa et al., *Dissecting context-dependent cancer vulnerabilities using Perturb-seq*, bioRxiv
10.64898/2026.08.24.746802, 2026-08-24) predicts line-specific responses from cell **state** in the mechanistic sense
(p53 activation after UBE3A/MDM2 knockout only in the HPV-E6 cervical lines SKGII and C4I; IER3IP1 dependence from a
hemizygous 18q/11q co-deletion). **Does the paper's kind of state feature predict the residual where expression PCs did
not, on our data, our statistic, our controls?**

**Experimental units, named.**
- **Cell line** (12 working Broad-16 lines) — the unit of the leave-one-line-out (LOLO) split and of every permutation
  in §5 and §6. n = 12; every interval that resamples lines is a 12-point interval and is reported as such.
- **Target** (Cas9 knockout, 82 usable targets by 268's ≥ 6-cell rule) — the unit of the paired bootstrap in §4
  (4,000 draws; each target carries its 12 folds), exactly as 268.
- **Cell** — the unit of the split-half ceiling only (268's halves A/B, seed 268).
- Not available on this data: guide-level replicates for the LOLO (used only in 268's robustness), a line measured
  twice (so line-level batch and biology are not separable; 268 §7.1).

**Reserved lines.** JVE127, MG63, SLR23, UMRC3 are SPENT (claim C2, 283) and are **not read**. `load268.py`'s assertion
is kept unmodified. Every number below is on the 12 working lines: C4I, HCT15, HKGZCC, KMRC20, KYSE140, KYSE450, ONE58,
RPMI7951, SKGII, SNU8, T3M5, U343.

---

## 2. The frozen state features

All are build-time observable for a recipient line: none is derived from any perturbation response. Feature values are
computed **after** this file is hashed.

### F-p53 — p53 functional status (DepMap Chronos gene effects)

Source: `experiments/269-utility-router-feasibility/work/depmap/CRISPRGeneEffect.csv` (DepMap 24Q2 Public, figshare
25880521, sha256 `d1551491…` recorded by 269), rows = ArxSpan ids from
`data/external/broad16line/metadata/cell_line_metadata.csv`. All 12 working lines are present (verified before this
file was written; a coverage check, not a value).

**Sign convention (stated explicitly because the mandate's two clauses conflict).** In DepMap Chronos, a line with
*functional* p53 (a) depends on MDM2, because MDM2 restrains p53 (strongly negative MDM2 effect), and (b) *gains*
from losing TP53 (positive TP53 effect). The mandate's clause 1 says "strong MDM2 dependency and a non-negative TP53
effect is p53-wild-type" (consistent with this) and then "a TP53 knockout benefit is p53-deficient" (the opposite).
I freeze the biologically and conventionally correct reading: **MDM2 dependency + TP53 knockout benefit ⇒ p53
functional (wild-type)**, and I record the conflict here rather than silently choosing.

- **Binary**: `p53_wt = 1` iff `MDM2_effect ≤ −0.5` **and** `TP53_effect ≥ 0`; else `p53_wt = 0` (p53-deficient).
- **Continuous**: `p53_score = TP53_effect − MDM2_effect` (higher = more p53-functional).
- **Degeneracy contingency, frozen now**: if `p53_wt` puts fewer than 2 of the 12 lines in either class, the binary is
  reported as degenerate and `p53_score` takes its place in every arm that uses "p53 status"; both are always reported.
- **Agreement check (reported, not used to define the feature)**: DepMap 24Q2
  `OmicsSomaticMutationsMatrixDamaging.csv` TP53 column > 0 (damaging TP53 mutation ⇒ p53-deficient), and
  `OmicsSomaticMutationsMatrixHotspot.csv` TP53 > 0. Agreement of the dependency call with the mutation call is
  reported per line.

### F-HPV — HPV-E6 status

`hpv = 1` for **SKGII** and **C4I**, `0` for the other ten. Source: the source paper's own statement as quoted in
`experiments/287-novelty-audit/NOVELTY-AUDIT.md` row 10 and §0.3 — "*p53 pathway activation on UBE3A and MDM2 only in
the HPV-E6 cervical lines SKGII and C4I*" (287's paraphrase of [BROAD16]; the preprint was not re-read here, and the
provenance is recorded as second-hand).

**Confound declared in advance**: SKGII and C4I are the *only* two Cervix lines among the 12 working lines, so on this
data `hpv` is numerically identical to a "Cervix lineage" indicator. Any effect attributed to HPV-E6 here is not
separable from cervical lineage. This is stated in every reading of an HPV result.

### F-E2F, F-P53H — basal hallmark programme scores

MSigDB hallmark sets from `data/external/broad16line/metadata/hallmark_symbols.gmt`.
- `e2f_g2m` = HALLMARK_E2F_TARGETS ∪ HALLMARK_G2M_CHECKPOINT.
- `p53_path` = HALLMARK_P53_PATHWAY.

Score of line *l* = (mean over the set's genes of *l*'s **basal** control log-profile) − (mean over **all** kept genes
of that profile). The basal profile is the 10 OR-control knockouts summed, `logprof` as in `load268.py`; the gene axis
and the CPM ≥ 5 filter are refit **inside the 11 training lines of each fold**, exactly as 268 fits its PCs
(`lolo268.py` line 23, `filter_lines=train`). Standardisation (mean, sd) is also fitted on the 11 training rows.

### F-CN — arm/region-level copy number (the paper's 18q/11q feature)

Acquisition (≤ 30 min, completed before this file was hashed; coverage verified, values not inspected): DepMap 24Q2
Public segment tables, figshare article 25880521, downloaded to `work/depmap/`:

| file | bytes | sha256 |
|---|---:|---|
| `OmicsCNSegmentsProfile.csv` | 31,571,033 | `4ad1d5dbabf804aa81007bc30bfac2a655950c60be5614e8565985d2697b797c` |
| `OmicsAbsoluteCNSegmentsProfile.csv` | 29,622,816 | `097bb2b557aacb02d675a2bc47173019ab7fcd3e2c32348da7a5f7decc66bfb3` |
| `OmicsProfiles.csv` | 214,983 | `e0125eebaca48c4f084218e372ca2da70a2de3ab3f6756462f3daecc258ea4aa` |
| `OmicsSomaticMutationsMatrixDamaging.csv` | 135,585,447 | `4f2e80d9fbf282600920dd16a9ad4782fb565ad9406cf1fd1dc5a8137d7fee46` |
| `OmicsSomaticMutationsMatrixHotspot.csv` | 3,874,822 | `a767624f6ee1997d2d6644f37796bbe486f58687c4c297ccd1972050b9c73468` |

Total 200.9 MB (≤ 1 GB). All 12 working lines have a CN-segment profile and a mutation row (coverage checked; no value
read). Profile → ModelID via `OmicsProfiles.csv`; if a model has more than one profile with segments, the **wgs**
profile is preferred, else the first by `ProfileID` sort.

Region coordinates come from `data/external/broad16line/data/gene_locations.csv` (Ensembl karyotype bands, hg38): a
region is the [min gene start, max gene end] over genes whose `Chromosome/scaffold name` and `Karyotype band` match.
Regions, chosen to match the paper's "hemizygous 18q21.1-32 / 11q23-4 co-deletion" as quoted by 287:
- `cn_18q21` — chr18, bands beginning `q21`.
- `cn_11q23_24` — chr11, bands beginning `q23` or `q24`.
- Secondary, reported: `cn_18q` and `cn_11q` (whole arms; all bands beginning `q`).

Value = **log2 of the length-weighted mean `SegmentMean`** over segments overlapping the region, from
`OmicsCNSegmentsProfile.csv` (relative CN; 1.0 = ploidy-neutral). Binary
`codel = 1` iff the length-weighted mean relative CN of **both** `cn_18q21` and `cn_11q23_24` is `< 0.85`.
`OmicsAbsoluteCNSegmentsProfile.csv` is used only as a reported cross-check of the calls.

### F-lin — lineage

`OncotreeLineage` from `data/external/broad16line/metadata/cell_line_metadata.csv`. **Used only** as (nnL) "nearest
same-lineage line's residual" and as a permutation label. Never a one-hot in a ridge with 11 training rows.

### F-depth — sequencing depth (NUISANCE, never a predictor)

From `data/external/depmap26q1/extracted/processed/cell_quality_covariates.csv`: `avg_umis_per_ko` (UMIs per cell) and
`total_umis_per_ko` (UMIs per perturbation — the variable of the source paper's power law, r = 0.85, elbow ≈ 5 M, and
of programme rule 9), plus `n_cells_per_ko`. Reported next to every fold's readout; never entered into any model.

---

## 3. Arms (frozen)

The fold loop, the ridge (λ = 1 on standardised features), the residual programme basis (J = 5, SVD of the training
interaction residual), the gene filter, the centring, the split-half ceiling and the three readouts are **exactly**
`lolo268.py`'s; only the feature matrix `Z` changes. Every fitted object is refit inside the 11 training lines
(programme rule 3).

Reference arms, re-run unchanged: **(a)** target-only (`r̂ = 0`), **(c1)** 2 basal PCs / per-gene ridge, **(c2)** 2 basal
PCs / programme-level, **(nn)** nearest basal line's residual, **oracle line offset** (uses held data; reference only).

New arms, each in a **per-gene** (`.g`, 268's c1 form) and a **programme-level** (`.p`, 268's c2 form) variant:

| arm | features |
|---|---|
| **S1** (headline) | `p53_wt`, `hpv`, `e2f_g2m`, `p53_path`, `cn_18q21`, `cn_11q23_24` |
| **S1n** | S1 without the two CN features (the mandate's core four) |
| **S2-\*** | each single feature alone: `p53_wt`, `hpv`, `e2f_g2m`, `p53_path`, `cn_18q21`, `cn_11q23_24`, `p53_score`, `codel` |
| **S3** | S1 + the 2 basal PCs (8 features) |
| **S4** | S1 with `p53_wt` replaced by the DepMap **damaging-mutation** p53 call |
| **(nnL)** | nearest **same-lineage** training line's residual; defined only for held lines that have a same-lineage partner among the 12 working lines — evaluated on those folds only, with the fold count reported |

**Controls for every new arm** (268's): **WRONG recipient** — the fitted model driven by each other training line's
feature vector, mean over the 11; **SHUFFLED labels** — feature rows permuted across the 11 training lines before
fitting, then predicted with the correct held feature vector, 5 permutations (rng seed 268).

**Readouts** (268's three): cosine of the predicted residual to the held residual; Δ cosine of the FULL prediction
`m + r̂` to the held response over target-only; Δ within-line retrieval. Paired bootstrap over **targets**, 4,000
draws, each target carrying its 12 folds.

### Pass rule (frozen)

An arm **predicts the line-specific residual from state** iff

1. its predicted-residual cosine beats **both** controls — paired difference (arm − WRONG recipient) and
   (arm − SHUFFLED labels) with the 95 % bootstrap interval entirely above 0; **and**
2. its **full-response Δcos is > 0 with the 95 % interval above 0** on the **programme-level** (`.p`) readout.

Headline verdict is read off **S1**. S1n / S2 / S3 / S4 / (nnL) are reported under the same rule.

---

## 4. Multiplicity

Cells screened, counted in advance: 6 feature sets (S1, S1n, S3, S4 + 8 single-feature S2 arms = 12 feature sets)
× 2 model forms (`.g`, `.p`) = **24 arm-cells**, plus (nnL) = 25, each with 3 readouts. No per-cell p-value is read as
a discovery on its own; the headline is S1 alone, fixed before any number. The remaining cells are descriptive and
labelled as a screen of 25.

---

## 5. P1 — the paper's p53/HPV prediction, on our residual object

**Claim (frozen).** After **UBE3A** and **MDM2** knockout, the HALLMARK_P53_PATHWAY score of the *line-specific
residual* is higher in the HPV-E6 lines {SKGII, C4I} than in (i) the p53-wild-type non-HPV lines and (ii) the
p53-deficient lines.

**Construction.** Full-data (not half) responses on the 12 working lines, 268's `build_responses` with
`min_cpm = 5` fitted on all 12 lines and the 100 Broad target genes removed. Residual `r[t,l] = R[t,l] − mean_l R[t,l]`.
Statistic `S(t,l)` = mean of `r[t,l]` over HALLMARK_P53_PATHWAY genes − mean of `r[t,l]` over all kept genes.
Inclusion: a (target, line) cell needs ≥ 6 assigned cells (268's rule); cells per (target, line) are reported. MDM2 is
**not** among 268's 82 usable targets (it fails the ≥ 6-cell rule on at least one line), so its included-line set may be
smaller than 12; the permutation is over the included lines only.

**Test.** Primary statistic = the mean of the two targets' `S` per line. Exact permutation over **all line pairs**
(C(n,2) with n = the number of included lines): one-sided p = fraction of pairs whose 2-line mean ≥ the observed
{SKGII, C4I} mean.
**P1 PASSES** iff p ≤ 0.05 **and** both orderings hold (HPV mean > p53-WT-non-HPV mean, and HPV mean > p53-deficient
mean). Per-target results are reported separately as secondary.

**Specificity check (frozen)**: the identical statistic on the E2F/G2M direction, where no HPV effect is predicted.

---

## 6. P2 — the paper's IER3IP1 co-deletion prediction

**Claim (frozen).** Lines carrying the 18q21/11q23-24 co-deletion show a **stronger** IER3IP1 knockout response than
lines without it.

**Readouts.** (a) *Fitness*: DepMap Chronos IER3IP1 gene effect (more negative = stronger dependency). (b)
*Transcriptome, ours*: the norm of the line's IER3IP1 full response divided by that line's median response norm over
the 82 usable targets (a per-line scale/depth normalisation), and the same for the line-specific residual.

**Test.** One-sided exact Mann–Whitney between `codel = 1` and `codel = 0` lines in the predicted direction, plus a
one-sided Spearman of each readout against the continuous co-deletion score `(cn_18q21 + cn_11q23_24)/2` (more
negative = more deleted). **P2 PASSES** iff the fitness readout is in the predicted direction at p ≤ 0.05 **and** the
transcriptome readout is in the predicted direction (any p). If fewer than 2 lines are called `codel = 1`, the group
test is declared **UNTESTABLE** and only the continuous correlations are reported.

---

## 7. Per-target gating (frozen)

Per fold, two fixed directions on the fold's gene axis: `u_P53` (HALLMARK_P53_PATHWAY membership indicator, mean-
centred over kept genes, unit-normed) and `u_E2F` (E2F_TARGETS ∪ G2M_CHECKPOINT, same construction). Nothing is fitted
in them.

For each target *t* and direction *u*: training projections `y[l] = ⟨r_tr[t,l], u⟩` over the 11 training lines;
predictors — for `u_P53` the 2-vector (`hpv`, `p53_wt`), for `u_E2F` the scalar `e2f_g2m` — standardised on the
training rows; ridge λ = 1 (no intercept, `y` centred on the training rows); predict the held line's projection and
take its **sign**. A fold is correct iff the predicted sign equals the sign of the held line's true projection.

- Per-target accuracy over 12 folds. **"Beats chance" = accuracy ≥ 9/12** (one-sided binomial p ≤ 0.073 at 0.5).
- Reported: the fraction of the 82 targets meeting it, next to (i) the same fraction under the SHUFFLED-label control
  (5 permutations of the training rows' features) and (ii) the binomial expectation 0.0730.
- Pooled: mean accuracy over targets × folds with a paired bootstrap over targets against the shuffled control.

---

## 8. Depth conditioning (programme rule 9)

Every fold's readout is reported next to the held line's `avg_umis_per_ko`, `total_umis_per_ko` and `n_cells_per_ko`.
Spearman of the per-fold S1 gain (both `cos_res` and `Δcos_full`) with each depth variable over the 12 folds is
reported with the n = 12 caveat, together with the deep-half / shallow-half medians (split at the median
`total_umis_per_ko`). No pass rule is conditioned on depth; depth never enters a model.

---

## 9. Known-answer rows

- **KA-1 (already run, contains no new information).** `work/ka289.py` is `lolo268.py` with the output paths
  redirected into `work/`; K = 2. It must reproduce R268 §4 exactly. **Result: PASS at 0.0e+00** — all 788 numeric
  entries of the summary JSON and all 984 × 45 CSV cells identical to
  `experiments/268-broad16-context-residuals/work/lolo268_K{2}.{json,csv}`; (c1) cos_res +0.0181 [+0.0136, +0.0229],
  Δcos_full −0.0071, Δretrieval +0.0003; ceiling median 0.0455. Abort condition (did not fire): any deviation.
- **KA-2 (instrument check for F-p53, to run).** Across all 1,150 DepMap 24Q2 models, the TP53 and MDM2 Chronos gene
  effects must be **negatively** correlated (Pearson ≤ −0.3), and models with `MDM2 ≤ −0.5` must have a higher mean
  TP53 effect than models with `MDM2 > −0.5`. If either fails, the p53 dependency proxy is reported as unvalidated and
  the mutation-based call (S4) becomes the primary p53 feature.
- **KA-3 (feature-pipeline check, to run).** Re-running the S1 fold loop with the feature matrix replaced by the two
  basal PCs must reproduce KA-1's (c1)/(c2) rows to ≤ 1e-12, proving the generalised loop is 268's loop.
- **KA-4 (hallmark parse).** HALLMARK_P53_PATHWAY, HALLMARK_E2F_TARGETS and HALLMARK_G2M_CHECKPOINT must each have
  200 symbols in the `.gmt` (MSigDB's published size); the number intersecting the fold gene axis is reported.

---

## 10. What this experiment will NOT do

No build, no `vcc`, no `cell-eval2`, no write under `data/` or `submissions/`, no read of the four reserved lines, no
proposal of a model for the competition. A negative is a result: if S1 fails the §3 pass rule, F5's amended headline
("not predictable from basal expression by a linear model at 11 training lines") is extended to the paper's own state
features **at this construction and this n** (12 lines, 82 targets, 11 training rows), and nothing more.

# 304 — the interrupted F15 audit, completed

> Coordinator review: `COORDINATOR-REVIEW-2026-09-07.md` governs integration.
> The numerical audit reproduces, but broad depth-law refutation, causal
> state explanations and exhaustive novelty claims are not adopted.
> The original worker report below is retained as the audit record.

**Result in one sentence.** F15's cell-line reading survives and is in fact stronger
than 302 wrote it (both arms are H1/WA01, not "H1/HUES"), but **nothing in F15 is
NEW**, its "same state" premise is wrong (the overexpression arm is 7 days in a
**differentiation** medium), and three of its four headline by-products do not
survive: the **depth law is a cross-target association that the same data refute
within targets**, the **geometry result is mostly a target-shuffled baseline**, and
the **isoform term is the source paper's own biology, not a construct ceiling**.

**Scope kept.** No reserve was opened: the spent Joung atlas `X` was not re-read and
`pseudo_bulk_perturbation_database.h5ad` (R5, 303's) was never touched. Every number
below comes from 302's saved `work/sums302.npz` + the already-used `gse295214_h1`
bank, or from primary text/source code. Nothing outside
`experiments/304-f15-audit/` was written; no commit, no publication, no submission,
no builder, no subagent. Peak RSS **1.34 GiB** (cap 6 GiB). One 200-draw permutation
and one 4-level shard ladder, 2.1 s wall.

Deliverables: **`NOVELTY-AUDIT.md`** (21 rows: exact local claim + location, primary
source with DOI and date, what the source shows, verdict, proposed wording),
**`work/refs304.md`** (sources, including the two that could **not** be opened),
**this file**. Audit code and outputs: `work/audit304.py`, `work/audit304.json`,
`work/audit304.log`.

---

## 1. Known-answer row

Everything computed here rests on one row, run first:

| row | reading | verdict |
|---|---|---|
| **KA-1** | rebuild the ORF residual from `sums302.npz` + the `gse295214_h1` bank and recompute 302's per-TF statistic: max \|ΔS_g\| over 175 scored TFs vs `symmetry302.csv` = **4.85e-07**, against a CSV stored at `%.6f` (rounding bound 5.0e-07) | **PASS** |
| **KA-1b** | the rebuilt median S over the 148 primary TFs vs `stats302.json`'s full-precision **−0.0006752075** | **PASS at 0.0e+00** |
| **KA-1c** | my within-target ladder at k = 8 shards/side reproduces 302's five published per-bin constants `sblaw302.json` to the digit: **6,942 · 11,127 · 9,562 · 7,741 · 9,592** | **PASS** |

So the re-analysis below is 302's own arithmetic on 302's own objects, not a
different pipeline.

## 2. Provenance (mandate item 1) — sent as the early note

- **Both arms are H1 hESC.** Joung et al. 2023 STAR Methods: *"We chose H1 hESCs for
  the TF Atlas because it is a commonly used stem cell line"*; key-resources table
  `H1 hESCs · WiCell · WA01`. **HUES66 was used only for the iNP/iCM sub-studies**,
  which are not in the atlas. GSE295214's SOFT file: `cell line: H1 hESC`,
  `genotype: CRISPRi`, `10x 5' scRNA-seq`. **The planned same-line interpretation is
  NOT invalidated — it is under-stated.** Correct "H1/HUES hESCs" to "H1 (WA01)".
- **The same state premise is wrong.** The atlas is *"cultured cells in STEMdiff APEL
  media for **7 days**"* after transduction — a directed-**differentiation** arm, which
  is the paper's entire subject. The CRISPRi arm's medium and days-post-transduction
  are **not established**: biorxiv.org returned HTTP 429/1015 for the JATS source and
  the HTML full text inside the time box. F15's "one line, one state" must become
  "one line, **two assay states**, one of them unread".
- **Platform and DOI.** SHARE-seq combinatorial indexing, >1.1 M profiles at 3,761
  UMI/cell; the reserve file's 1,145,823 cells is the **full** set, not the authors'
  671,453-cell representation-balanced subsample. The DOI carried in the 301/304
  chain (`10.1016/j.cell.2023.01.026`) is **wrong**; the paper is
  **`10.1016/j.cell.2022.11.026`**, Cell 186(1):209–229.e26, 5 Jan 2023.
- **The knockdown side is a third-party re-processing**, not GEO counts: the
  `dataset == "Pan2025_hESC"` slice of `xingjiepan/SCMG_data` (rev `1d3f28d9…`),
  effect `log1p(CPM1e4 perturbed) − log1p(CPM1e4 control)` with a `measure_mask`.
- **For 303 (Nakatake/R5).** The SCMG database's hESC overexpression rows are built
  by `SCMG_scripts@aa848b9d …/pseudobulk_perturbation/hESC_TF_screen.ipynb` from a
  single flat file `TF_hESC_screen.csv`: columns are TF symbols, **the control column
  is literally `Emerald`**, values are **log10** (`10 ** df`), rescaled to 1e4, genes
  with any value `< 1e-6` dropped, four symbols hand-renamed
  (`H2AFZ→H2AZ1`, `HIST2H3C→H3C14`, `SSX6→SSX6P`, `T→TBXT`), `perturbation_sign = 1`.
  It is a **bulk, per-transgene summary matrix**, not single cells and not per-guide.
  Which instrument produced it is **UNRESOLVED**: Nakatake 2020 (Cell Reports
  31:107655) is not in PMC/Europe PMC and cell.com returned 403, and the *Scientific
  Data* DOI in the mandate (`10.1038/s41597-021-00806-1`) **does not resolve** in
  Crossref or Europe PMC. Secondary sources give the line as the male **SEES3** hESC
  line with readout at **48 h** — **unverified**. If that holds, a Nakatake-vs-Pan
  test is **cross-line and cross-timescale** (SEES3/48 h vs H1/CRISPRi), not a
  second instance of F15's design, and C7 should be worded accordingly.
  Separately: the same SCMG database already contains a **`Joung_TFScreen_HS_2023`**
  pseudobulk notebook — an OE-vs-KD test inside one processing pipeline is available
  and is a cleaner comparison than 302's cross-pipeline one.

## 3. The three measured corrections

### 3a. The depth law is refuted by 302's own shards (item 4)

302 did **not** downsample the same TFs. `work/sblaw302.json` bins **different TFs**
by cell count and derives `c = n(1−ρ)/ρ` per bin; the quoted "[6,942, 11,126]" is the
**min–max across bins**, not an interval, and there is no fit and no residual.

Running the missing within-target ladder on 302's 16 saved shards (nested subsets,
**same TFs**, k shards per side):

| k (shards/side) | median cells in 2k shards | split-half r | ρ_full(2k) | implied c |
|---:|---:|---:|---:|---:|
| 1 | 59.9 | +0.0422 | +0.0810 | 680 |
| 2 | 121.4 | +0.0337 | +0.0653 | 1,738 |
| 4 | 241.9 | +0.0264 | +0.0515 | 4,458 |
| 8 | 474.5 | +0.0246 | **+0.0480** | **9,404** |

Reliability **does not improve** as the same TF's cells go 60 → 475; c rises almost
proportionally to n, which is what a *non*-Spearman–Brown curve looks like. The same
holds inside **every one of 302's five bins** (3,000–20,000 cells: c = 407 → 1,045 →
3,574 → 9,592; 100–200: 834 → 1,457 → 4,459 → 6,942). The k = 8 column *is* 302's
published table (KA-1c), so the cross-target ladder and the within-target ladder are
computed from the same objects and disagree.

Two further defects: cells-per-TF in a pooled 7-day differentiation screen is itself
a phenotype (Joung et al. down-sampled *"to ensure even representation"*), so the
cross-target association is confounded with the TF's effect on proliferation; and the
"≈ 3,100 cells" figure solves ρ_**full** = 0.243 when 0.243 is the Flex reference's
**half-sample** number — like-for-like it is **6,139**, and since c is not identified,
neither number is supported.

### 3b. The geometry result is mostly a shuffled-pairing baseline (item 6)

302's evidence that ρ(\|S\|, `gen_pc20`) = **+0.3823** "is not a measurement artefact"
is a set of partials that never break the TF pairing. The correct control is the one
302 already used for C6a — a derangement of the pairing:

| statistic | reading |
|---|---|
| observed ρ(\|S\|, gen_pc20), n = 142 | **+0.3823** (reproduces 302 exactly) |
| target-shuffled ρ(\|S_perm\|, gen_pc20), 200 derangements | **+0.293 ± 0.067**, 95 % [+0.167, +0.427] |
| draws ≥ observed | **10 % of 200**; z = **+1.32** |
| ρ(gen_pc20, median over other TFs' ORF rows of \|cos\|) | **+0.517** (mean +0.556, sd +0.582, q90 +0.561) |
| ρ(\|S\|, gen_pc20) partialling out that alignability | **+0.126 … +0.165** |
| ρ(gen_pc20, ·) for 302's chosen covariates | reliability +0.079 · kd_mag +0.315 · orf_mag −0.069 · log cells +0.011 |

≈ 77 % of the association survives pairing an ORF response with a **different TF's**
knockdown row. `gen_pc20` predicts a knockdown row's alignment to the **whole ORF
ensemble**; the pairing-specific remainder is +0.13…+0.17 and is not tested against
its own null. 302's partials could not have moved — its covariates are nearly
orthogonal to `gen_pc20`. Multiplicity is not the issue (≈ 18 descriptive
correlations; +0.38 at n = 142 is p ≈ 3e-6 against independence) — **independence is
the wrong null**.

### 3c. The excluded "strong mirror" is an assumption, not a result (item 3)

`strat302.py:18` substitutes `RHO_KD = 0.243` — 273's **half-sample** split-half
figure from a *different* dataset (H1 Flex, 200/200 cells) — into
`sqrt(ρ_OE_full · ρ_KD)`, while the ORF side is Spearman–Brown-doubled to full
sample first. Two sides, two conventions. Recomputed on the KA-1 rebuild:

| ρ_CRISPRi (assumed) | ceiling | implied TRUE median cosine [95 %] |
|---|---:|---|
| 0.243 — 302's substitution (a half-sample figure used as full) | 0.1080 | −0.0062 [−0.0435, +0.0150] |
| **0.391 = SB(0.243)** — the like-for-like carry | **0.1370** | **−0.0049 [−0.0343, +0.0118]** |
| 0.10 — a poor bank | 0.0693 | −0.0097 [−0.0678, +0.0234] |
| 0.05 — a poor bank | 0.0490 | **−0.0138 [−0.0958, +0.0330]** |
| 0.60 — a good bank | 0.1698 | −0.0040 [−0.0277, +0.0095] |

At ρ_CRISPRi = 0.05 the same data are compatible with a true median mirror of
**−0.10**. The exclusion is carried entirely by an unmeasured number.
Two further construction faults: the "[95 %]" columns are the bootstrap interval of
**median S** divided by a **constant** built from the **median** reliability — the
reliability's own error (each per-TF value is a mean of 8 noisy split-halves, median
+0.025, 5.4 % ≤ 0) is not propagated, and a median of ratios is not a ratio of
medians. And the residualisation removes any mirror **shared across TFs** from both
sides before the cosine is taken, so the statistic can only see TF-*specific*
mirroring — the shared "dosage-monotone" component the programme predicted is
centred out by construction.

## 4. What survives

- **The line identity** (both arms H1/WA01) — stronger than written.
- **The measurement itself**: median S = −0.0007 [−0.0047, +0.0016] against a
  shuffled null of +0.0008 ± 0.0017, z −0.83, 16 TFs below −0.05 vs 6 above +0.05.
  Reproduced here at 0.0e+00.
- **The honest narrowing of F13** (`R302:174-176`): the sign travels, the ordering
  claim does not, and on the signed statistic there is nothing (+0.058, interval
  across zero). §3b makes this stronger.
- **The negative as a bounded negative**: "no mirror is detected at this depth" is
  supportable. "Orthogonal" is not — a non-significant median is failure to reject,
  and 302 never states an equivalence bound. The defensible form is in
  `NOVELTY-AUDIT.md` row B2.
- **Prior art**: that gain- and loss-of-function readouts are not mirror images is
  **KNOWN** (Tian et al. 2021, *Nat Neurosci*, `10.1038/s41593-021-00862-0`:
  *"Overall, there was little overlap between CRISPRi and CRISPRa hits, consistent
  with our previous study comparing parallel CRISPRi and CRISPRa screens in K562
  cells"*; Nishiyama et al. 2013, *Sci Rep* 3:1390, `10.1038/srep01390`, the same
  laboratory as the Nakatake resource, on the asymmetry of TF repression vs
  induction). Tian's statement is about **survival-phenotype hit overlap**, not a
  transcriptome cosine — so 302's *quantity* (a transcriptome-wide cosine for 148
  matched TFs with a measured reliability bound) is **EXTENDED**, not KNOWN, and
  not NEW either: it is a new measurement of an established qualitative fact.

## 5. Verdict tally (details in `NOVELTY-AUDIT.md`)

21 rows: **0 NEW · 2 EXTENDED · 5 KNOWN (+1 retained) · 12 CONTRADICTED · 1
UNRESOLVED**. F15 contains no row that is NEW.

## 6. The three most consequential corrections

1. **Delete the Spearman–Brown depth law for ORF overexpression (F15 by-product (i),
   `FINDINGS §F15`, `READ-FIRST:60`).** 302 compared **different TFs in cell-count
   strata**; on the same TFs' own cells reliability is flat-to-declining over an 8×
   range (ρ_full +0.081 → +0.048 from 60 to 475 cells) and the implied constant rises
   from 680 to 9,404. "A single constant c ≈ 9,600 [6,942, 11,126]" and "≈ 3,100 cells
   per TF" should both go; the surviving statement is the cross-target association
   (ρ(reliability, cells) +0.678), which is confounded with the screen's own
   proliferation selection. **This also removes F15's claim to extend F14's law to a
   second modality** — F14's law rests on real cell ladders (285), 302's does not.
2. **Withdraw "it is not a measurement artefact" from the geometry by-product
   (F15 by-product (iii), `FINDINGS:382`, `READ-FIRST:60`).** Against a
   target-shuffled baseline of +0.293 ± 0.067 the matched-pair +0.382 is z +1.32 with
   10 % of 200 draws above it; `gen_pc20` predicts a knockdown row's alignment to the
   **whole** ORF ensemble (ρ +0.52) and the pairing-specific part is +0.13…+0.17.
   The defensible claim is "spectrally central knockdown rows align more with
   everything", not "geometry predicts the size of the cross-modality overlap".
3. **Make the excluded mirror conditional, and fix "one line, one state".** The
   "strong mirror excluded / −0.05 compatible" bound holds only if the
   `gse295214_h1` bank's reliability is ≥ ≈ 0.24; at 0.05 the data admit −0.10. And
   the two arms are the same **line** but not the same **state** — the ORF arm is
   7 days in STEMdiff APEL differentiation medium (Joung et al. 2023), which is a
   more likely explanation of a null cosine than "biological orthogonality".
   *(Runner-up, narrow: `READ-FIRST:60`'s "the term is knockdown failure, not
   off-targeting" over-reads 301's P = 0.68 at n = 34 — see `NOVELTY-AUDIT.md` G1;
   305 owns the expression test.)*

## 7. Paths written

```
experiments/304-f15-audit/RESULTS.md            (this file)
experiments/304-f15-audit/NOVELTY-AUDIT.md      (21 audit rows)
experiments/304-f15-audit/work/refs304.md       (sources, incl. 2 unreachable)
experiments/304-f15-audit/work/audit304.py      (KA-1 + items 3, 4, 6)
experiments/304-f15-audit/work/audit304.json    (all audit numbers)
experiments/304-f15-audit/work/audit304.log
experiments/304-f15-audit/work/hESC_TF_screen.ipynb   (SCMG source notebook, fetched)
experiments/304-f15-audit/work/scmg_scripts_tree.json (SCMG_scripts file tree @aa848b9d)
experiments/304-f15-audit/work/joung_pmc.html|.txt    (Joung 2023 full text)
experiments/304-f15-audit/work/PMC8254803.html|.txt   (Tian 2021 full text)
experiments/304-f15-audit/work/PMC3589720.html|.txt   (Nishiyama 2013 full text)
```

Nothing outside this directory was modified. `docs/FINDINGS-2026-09.md`,
`docs/CONFIRMATION-LEDGER.md`, `docs/READ-FIRST.md` and
`experiments/294-writeup/` carry the wording this audit corrects and are the
**coordinator's** to change.

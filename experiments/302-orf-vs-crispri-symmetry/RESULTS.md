# 302 — claim C6 on reserve R4: does ORF overexpression reverse CRISPRi knockdown in the same cell line?

**C6a (overexpression reverses knockdown, hESC): FAIL — VOID by the frozen KA-E clause**
(median S **−0.0007** [−0.0047, +0.0016]; share < −0.05 **10.8 %** (16) vs > +0.05 **4.1 %** (6);
null **+0.0008 ± 0.0017**, z −0.83);
**C6b (conservation predicts mirroring): FAIL — VOID** (ρ **−0.1143** [−0.2846, +0.0673]; partial −0.1057);
**C6c (geometry sign): FAIL — VOID** (ρ **+0.0583** [−0.1106, +0.2278]).

**What the measurement says, in one sentence:** in H1/HUES hESCs, over 148 transcription factors measured
both ways, an ORF overexpression response and a CRISPRi knockdown response are **orthogonal, not opposed**
(median cosine −0.0007, indistinguishable from the TF-shuffled null at +0.0008), and the reason is measured
rather than assumed — the ORF side's split-half reliability is **+0.025** (Spearman–Brown full-sample
**+0.048**), which caps any cosine at **0.108**, so the data exclude a true mirror stronger than about
**−0.04** at the median but cannot resolve a weaker one.

Reserve **R4 is SPENT** (`JoungZhang2023_atlas.h5ad`, sha256 `8bd00eac…`, one streaming pass, 77 s, peak RSS
7.04 GiB). Preregistration `PREREGISTRATION.md` sha256 `69c2a92415b9e61149bdbee6eb3e61c61a8b23711954d237fcf8bd19d341b9e0`,
hashed 2026-09-07T03:33:19Z and appended to `docs/CONFIRMATION-LEDGER.md` **before** the first byte of `X`.
No build, no `vcc`, no `cell-eval2`, nothing under `submissions/` or `data/`, no download. Panel coverage
0/300 — nothing competition-facing was selected, tuned or checked, and nothing is proposed for the model.

---

## 1. Known-answer rows, in the order they ran

| row | reading | verdict |
|---|---|---|
| KA-0 | atlas sha256 `8bd00eac883b469a1e82a2beca18f647e8023904920f233dbb8568e13503fc22` | **MATCH** |
| KA-A | GFP = code 3367, **49,539** cells; mCherry = code 3368, **78,274** cells (301's counts) | **PASS** |
| KA-B | **249** TFs in both the resolved ORF map and the H1 bank; **233** in 279's 17,106-gene held map; **148** with ≥ 100 ORF cells — the coordinator's metadata check, to the digit | **PASS** |
| KA-C | the H1 bank reproduces 279's `Ccross` on all **26** cross-line pairs involving `gse295214_h1`, **15,047** (pair, gene) cosines, max abs diff **0.0e+00**, identical gene sets, through 279's own `load_unit` | **PASS** |
| KA-D | 301's `joung_orf_map.csv` sha256 `a6b60664…`; 3,136 / 3,369 codes, 1,647 symbols, 99.4836 % of cells | **PASS** |
| **KA-E** | computed library size vs `obs['ncounts']`: max abs diff **9.96e+03** | **FAIL** (see §2) |
| KA-F | GFP-vs-mCherry no-signal row against the 148 knockdown rows: median \|cos\| **0.0097**, max **0.0461**, **0 %** above 0.05 | **PASS** |
| **KA-G** | the positive control: ORF residual's own-gene coordinate median **+0.3368**, **95.1 %** positive (144 TFs); CRISPRi residual's own-gene coordinate median **−0.6998**, **100 %** negative | **PASS** |
| **KA-H** | `X` integrality: max fractional part **0.5** | **FAIL** (see §2) |

KA-G is the row that matters for trust in the pipeline: **overexpressing TF *X* raises *X*'s own transcript
in 95.1 % of the cohort, and knocking it down lowers it in 100 %**. That is an independent expression-level
confirmation of 301's ORF-id resolution (which was built from barcodes and a supplementary table, never from
expression), of the cell→symbol assignment, of the gene axis and of the profile arithmetic.

## 2. KA-E and KA-H: the file's own count columns are inconsistent, and X is not integer-valued

`X` holds **half-integers** (max fractional part exactly 0.5), and the per-cell row sums do not equal
`obs['ncounts']`. Diagnosed from **free metadata plus the sums already collected** — `X` was not reopened
(`work/ka_e_diagnosis.json`):

- `sum(obs['ngenes'])` = `sum(var['ncells'])` = **2,673,739,132** = exactly `X`'s nnz. The sparsity structure
  I read is the file's, cell for cell and gene for gene.
- `sum(obs['ncounts'])` = **4.18355e9** but `sum(var['ncounts'])` = **4.04904e9** — the file's own two count
  columns disagree by **3.2 %**, and `var['ncounts']` is itself fractional. So the mismatch predates this
  experiment.
- Per scored group, computed / `obs['ncounts']` runs **0.882 … 1.182** (median 0.942), and it is
  **batch-keyed**: 0.940 for the 134 groups that are 100 % batch 0, 1.142 / 1.149 for GFP / mCherry (100 %
  batch 1), Pearson(batch-1 fraction, ratio) **+0.910**.

The reading: `obs['ncounts']` was recorded at a different processing stage from the `X` on disk, batch by
batch, and `X` carries fractional (EM- or multimapping-assigned) counts. **KA-E as I froze it — "max abs
diff 0" — tests a property this file does not have**, and KA-H's integer assumption is likewise wrong for
this matrix. Nothing in the claim's arithmetic uses `obs['ncounts']`: the profile
`log1p(5e4·S_s/L_s)` takes both `S_s` and `L_s` from the same matrix, so it is internally consistent.

Nevertheless the preregistration says, in writing: *"If KA-E or KA-G fails … the three verdicts are recorded
as **VOID**, not PASS/FAIL, and the numbers are reported as descriptive only."* That clause fires. **The three
verdict cells are VOID.** Every number below is reported, with the rules' verdicts stated as they would have
read, and the reserve is spent either way. R4 cannot be re-opened to fix the KA row.

## 3. The design as executed

- **Axis: 15,710 genes** = the H1 bank's 18,533 genes ∩ its row-constant `measured_mask` (16,268; `delta` is
  exactly 0 elsewhere) ∩ the atlas's 37,528 distinct `var` symbols.
- **Scored set T: 175 TF symbols** (≥ 50 pooled ORF cells and an H1 bank row), 419 ORF codes,
  **212,329 ORF cells** + 127,813 control cells = 340,142 cells accumulated. **Primary cohort T100 = 148**
  (≥ 100 cells; median 474, q10 156, max 14,339); secondary = 27 (50–99 cells).
- **ORF residual** `a_s = log1p(5e4·S_s/L_s) − mean_{t∈T}(…)`; the GFP+mCherry control profile cancels
  exactly, which is why the controls sitting entirely in batch 1 (while 94.1 % of scored cells are batch 0)
  cannot enter `S_s`. **CRISPRi residual** `b_s` = the bank's `delta` row minus the mean over the same 175
  rows (279's `load_unit` arithmetic). Own column deleted from both for the 144 / 148 TFs that have one.
- **`S_s = cos(a_s, b_s)`**; ORF split-half reliability from 16 seeded cell shards × 8 balanced partitions.
- One streaming pass over `X` (58 row blocks of 20,000 cells, 2.674e9 nonzeros), 77 s, peak RSS 7.04 GiB.

## 4. The three cells

| cell | statistic | reading | rule | verdict |
|---|---|---|---|---|
| **C6a** | median `S_s`, n = 148 | **−0.0007** [−0.0047, +0.0016]; 52.7 % negative; share < −0.05 **10.8 %** (16) vs > +0.05 **4.1 %** (6) | median < 0 with the interval below 0 **AND** the tail asymmetry | **FAIL** (condition 2 holds, condition 1 does not) → VOID |
| **C6b** | Spearman(`C_held`, −`S_s`), n = 142 | **−0.1143** [−0.2846, +0.0673]; partial given both magnitudes + reliability **−0.1057** | ρ > 0 with the interval above 0 | **FAIL** (sign inverted, interval spans 0) → VOID |
| **C6c** | Spearman(`gen_pc20`(H1), −`S_s`), n = 142 | **+0.0583** [−0.1106, +0.2278]; partial +0.0653 | ρ > 0 with the interval above 0 | **FAIL** (right sign, interval spans 0) → VOID |

**The null matters here.** 200 TF-shuffled pairings give median **+0.0008 ± 0.0017**; the observed −0.0007 sits
**z = −0.83** from it, with 19 % of null draws at or below the observation. The frozen C6a rule did not
require beating the null, and the mock run (a synthetic atlas with no signal, run before the access) showed
the null median is *not* exactly 0 for this construction — so the null, not the sign of the median, is the
reference a reader should use. On that reference C6a is a clean zero.

**Baselines for C6b/C6c** (reported, never in the rules), Spearman against −`S_s` on the primary cohort:
knockdown magnitude −0.0128 [−0.1821, +0.1577]; 290's `mag` −0.0381; ORF magnitude −0.0467; ORF reliability
+0.0799; `r_bar` −0.0455; ORF cells −0.0003; ORF UMI/cell −0.1032; isoform count −0.0004; map-wide `gen_pc20`
−0.0358; `breadth`(H1) +0.0143; map-wide `breadth` −0.0249; batch-1 fraction +0.1002. **Nothing predicts the
sign of the alignment.**

## 5. Why: the ORF side does not resolve a direction

| cells per TF | n | median cells | ORF split-half `r` | Spearman–Brown ρ_full | implied `c` (cells at ρ = 0.5) |
|---|---:|---:|---:|---:|---:|
| 100–200 | 28 | 154 | +0.0109 | +0.0216 | 6,942 |
| 200–500 | 49 | 321 | +0.0142 | +0.0280 | 11,126 |
| 500–1,000 | 30 | 633 | +0.0320 | +0.0621 | 9,562 |
| 1,000–3,000 | 22 | 1,636 | +0.0956 | +0.1745 | 7,741 |
| 3,000–20,000 | 19 | 6,222 | +0.2449 | +0.3934 | 9,592 |

`ρ(reliability, cells)` **+0.678**, `ρ(reliability, total UMIs)` +0.666, `ρ(reliability, UMI/cell)` −0.277 —
a clean Spearman–Brown ladder over a 40× cell range with a **single constant c = 9,562 cells
[6,942, 11,126]**. To reach the H1 Flex reference's +0.243 an ORF response needs **≈ 3,100 cells**; the atlas's
median resolved symbol has **110** (the cohort's median 474). 5.4 % of the cohort reads reliability ≤ 0.

**Attenuation.** `ρ(|S|, ORF reliability)` = **+0.469** [+0.331, +0.583] — the size of the alignment is set by
how well the ORF side is measured, as attenuation predicts. With the ORF side at ρ_full 0.048 and the CRISPRi
side assumed at **0.243** (273's H1 Flex figure — the `gse295214_h1` bank's own per-target reliability is not
measured here), the ceiling on any cosine is **0.108**:

| stratum | n | ORF `r` | ceiling | median `S` | implied TRUE mirror cosine [95 %] |
|---|---:|---:|---:|---:|---|
| all | 148 | +0.025 | 0.108 | −0.0007 | **−0.006 [−0.044, +0.015]** |
| reliability T1 | 49 | +0.006 | 0.052 | +0.0009 | +0.018 [−0.094, +0.096] |
| reliability T2 | 49 | +0.025 | 0.108 | −0.0009 | −0.008 [−0.037, +0.068] |
| reliability T3 | 50 | +0.145 | 0.248 | −0.0106 | −0.043 [−0.094, +0.040] |
| `r` ≥ 0.20 | 21 | +0.280 | 0.326 | −0.0180 | −0.055 [−0.196, +0.075] |
| `r` ≥ 0.30 | 10 | +0.408 | 0.375 | +0.0032 | +0.008 [−0.179, +0.179] |

The better-measured terciles drift negative (T1 +0.0009 → T3 −0.0106; `r` ≥ 0.20 −0.0180) but **no stratum's
interval excludes 0**, and the trend breaks in the top stratum and in the depth strata (`cells ≥ 2000`,
n = 27, reads **+0.0101**). Read it as: *a weak reversal of order −0.05 in true cosine is compatible with
these data; a strong one is not; and this measurement cannot separate the two.*

**The two-construct ceiling, in a new modality.** Two ORF splice-isoform constructs of the same TF, in the
same cells and the same batch, agree at a median cosine of **+0.021** over **365 pairs** (70 TFs, 219
isoform-level rows). That is the same near-zero construct term 296 measured for CRISPRi libraries once a
cross-gene null was subtracted (≈ 0.03–0.08) — with the honest caveat that TF splice isoforms are genuinely
different proteins, so this mixes a real biological difference with measurement noise and is a **lower**
bound on the construct term.

## 6. Secondaries (declared) and descriptive tails (post-hoc, labelled)

- **|S| (symmetry of either sign).** Median |S| **0.0163**. `ρ(|S|, gen_pc20(H1))` = **+0.3823
  [+0.2261, +0.5192]**, and the partial given ORF reliability, both magnitudes and log cells is **+0.3831** —
  unchanged. `ρ(|S|, C_held)` −0.1097 [−0.2825, +0.0753]; `ρ(|S|, breadth)` +0.0992 (partial +0.1713);
  `ρ(|S|, kd_mag)` +0.0429 (partial −0.0031); `ρ(|S|, reliability)` +0.4690 (partial −0.0155).
  **So F13's geometry feature predicts how much the two modalities' responses share a subspace, and it is not
  a measurement artefact — but it does not predict whether they oppose** (C6c, on the signed statistic,
  +0.058 with the interval across 0). Size, not sign.
- **Secondary cohort (50–99 cells), n = 27:** median `S` **+0.0024** [−0.0121, +0.0113], 40.7 % negative.
- **Per-isoform:** 219 rows over 70 TFs, median `S` **+0.0029**; a TF's isoform-level `S` values scatter by a
  median sd of **0.0083** against a cohort sd of 0.0418 (they share the same CRISPRi row, so this bounds the
  ORF-side construct contribution to `S`, not the reproducibility of a biological mirror).
- **Disattenuated median `S`** (ORF side only, `S/sqrt(r)`): **−0.0070**. Positively aligned fraction 47.3 %.
- **Most mirrored** (post-hoc, labelled): `ELOF1` −0.232 (7,994 cells, r +0.28), `DAXX` −0.174, `MYC` −0.117,
  `ZFP90` −0.117, `NANOG` −0.108, `ZNF281` −0.096, `TGIF1` −0.085, `SUPT4H1` −0.076, `KLF16` −0.074,
  `POU5F1` −0.067, `NHP2` −0.065, `YAP1` −0.064, `KLF4` −0.064, `TEAD2` −0.059, `TCF4` −0.057, `ZNF544`
  −0.051, `GTF3A` −0.043, `TCF3` −0.039, `PA2G4` −0.039, `YY1` −0.034. The hESC pluripotency network
  (`NANOG`, `POU5F1`/OCT4, `MYC`, `KLF4`, `TCF3`, `TEAD2`/`YAP1`) and core transcription machinery
  (`ELOF1`, `NHP2`, `SUPT4H1`, `GTF3A`) sit on the mirrored side — the dosage-monotone shape the programme
  predicted, present as a tail rather than as a cohort median.
- **Most aligned** (post-hoc, labelled): `MAX` +0.134 (8,305 cells, r +0.30), `SATB1` +0.097, `SUB1` +0.092,
  `TP53` +0.057, `GABPB1` +0.052, `NFYC` +0.051, `MAZ` +0.042, `ZNF574` +0.041, `RUNX1T1` +0.039,
  `HMGN3` +0.038. `MYC` mirrors while its obligate partner `MAX` aligns — a description, not a finding.

## 7. What this does and does not establish

**Establishes.** (i) In one cell line, across 148 TFs, overexpression and knockdown responses are orthogonal
at the cohort median, and the tail asymmetry (16 TFs below −0.05 vs 6 above +0.05) is the only sign of the
predicted reversal. (ii) Cross-line conservation does **not** predict mirroring (C6b's point estimate has the
wrong sign). (iii) F13's `gen_pc20` predicts the **magnitude** of cross-modality overlap (+0.38, partial
+0.38) but not its **sign** (+0.06, interval across 0) — F13 is narrowed again, in the same direction 291 and
295 narrowed it: the sign of the geometry rule travels, the ordering claim does not, and here even the sign
applies only to |S|. (iv) A Spearman–Brown depth law for ORF overexpression in hESCs with c ≈ 9,600 cells,
extending 285's F14 law to a second perturbation modality. (v) A ≈ +0.02 two-construct ceiling in that
modality.

**Does not establish.** The three verdicts are VOID by the frozen clause, so none of C6a–C6c is a decided
preregistered claim. The negative is bounded by measurement, not by biology: at these cell counts a true
mirror cosine of −0.05 would have been invisible. ρ_CRISPRi is assumed, not measured, so every disattenuated
number moves with it. One line, one state, two batches, no biological donor. Panel coverage is 0/300 and
nothing here touches the competition track.

**What would decide it.** ≈ 3,000 ORF cells per TF (a 6× deepening of the atlas's cohort, or restriction to
the 19 TFs already above 3,000 cells plus new data), a measured per-target reliability for the
`gse295214_h1` CRISPRi bank so the ceiling is not assumed, and a second overexpression library so the
construct term can be separated from the isoform difference.

## 8. Files

`PREREGISTRATION.md` (`69c2a924…`) · `work/PREREGISTRATION.sha256` · `meta302.py` → `work/meta302.json`,
`work/cohort302.csv`, `work/axis302.txt` · `ka302.py` → `work/ka302_C.json` (KA-C) · `lib302.py` (frozen
design) · `mock302.py` (the synthetic-atlas dry run executed **before** the access; `work/*_mock.*`) ·
`access302.py` → `work/access302.json`, `work/access302.log`, `work/sums302.npz`
(sha256 `01c2a69669ff17b694b349f1b50783b0caeec686a8781fea9e9b700f76d51671`, the only object that lets a later
reader re-derive anything without re-opening the spent reserve) · `analyse302.py` → `work/stats302.json`,
**`work/symmetry302.csv`** (sha256 `9c15824aed80d984b42a548ab5f8b46c4929adec1752038c6bf1b0bf043ace9f`; per TF:
cells, UMI/cell, isoforms, ORF reliability, `S_g`, `C_held`, `kd_mag`, `mag290`, `gen_pc20` H1 and map-wide,
`breadth`, ORF magnitude, `r_bar`, batch-1 fraction, own-gene coordinates), `work/isoform302.csv` ·
`desc302.py` → `work/desc302.json` · `strat302.py` → `work/strat302.json` · `work/sblaw302.json` ·
`work/absS302.json` · `work/ka_e_diagnosis.json`.

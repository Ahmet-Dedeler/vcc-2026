# 303 — PREREGISTRATION **DRAFT** for claim C7 (reserve R5)

**Status: DRAFT. NOT FROZEN. R5 IS UNREAD.**
This file is written before any Nakatake `X` or `layers/control` value has been
read. It is not a release. The coordinator alone freezes it (by renaming/adopting
it as `PREREGISTRATION.md`, hashing it, appending C7 to
`docs/CONFIRMATION-LEDGER.md`, and writing `work/RELEASE.json`).
`experiments/303-*/access303.py` refuses to open the reserve without that release.

Every number below is from metadata, measurement masks, the already-used
`Pan2025_hESC` rows, saved artefacts of 276 / 279 / 290 / 302, and the data
producer's own source code. None is an R5 expression value.

---

## 0. What changed from the interrupted draft, and why

The original `MANDATE.md` assumed the two libraries sit **in the same cell line**
("hESC, the map's `gse295214_h1` bank"), that the Nakatake side is **RNA-seq**,
and that **bulk depth removes the reliability ceiling**. The pre-access checks
falsify or fail to support all three (§2, §7). The claim is therefore restated as
a **cross-study, cross-line** comparison, its estimand is the **observed**
alignment of two differently-measured responses, and no latent biological
equivalence is promised from a null. Six further amendments (co-centring, the
shuffled-pairing null, Holm, the source-excluded conservation feature, the
OE-randomised geometry null, and the demotion of the 302 comparison to
descriptive) are listed in `PREACCESS.md` §4.

## 1. The question

Over transcription factors measured **both** as a bulk transgene overexpression
(Nakatake et al. 2020, SEES3 hESC, 48 h Dox) and as a single-cell CRISPRi
knockdown (Pan et al. 2025, H1 hESC, 10x Perturb-seq pseudobulk), is the
observed overexpression response the **reverse** of the observed knockdown
response, does **cross-line conservation** predict how mirrored a TF is, and does
**response geometry** (F13's `gen_pc20`) predict the alignment's **sign** and its
**size**?

## 2. The estimand, stated honestly

`S_i` is the cosine between two **observed** response estimates from **two
different studies, two different hESC lines, two different assays and two
different perturbation modalities**, after both are mapped into one processing
convention. It is **not** a repeat measurement of one intervention. A value near
zero is consistent with (a) genuine orthogonality, (b) a real mirror destroyed by
line, protocol, assay, duration or construct differences, or (c) attenuation by
measurement error on either side — and **this design cannot separate them**,
because the Nakatake side in this file carries **no replicate structure** from
which a reliability could be estimated (§7, risk U3). No disattenuated cosine, no
"true mirror" bound and no reliability ceiling will be reported.

## 3. Frozen inputs (all hashed in `work/meta303.json`)

| object | value |
|---|---|
| source | `data/external/gse295214/pseudo_bulk_perturbation_database.h5ad`, 4,425,260,319 B, sha256 `0754d620…` (= 115's pin) |
| R5 rows | `obs.dataset == 'Nakatake2020_hESC_TF_screen'`, **510 rows, contiguous 184…693**, ordered index sha256 `85abd9a4…`, symbol list sha256 `a76bb7ec…`, `obs.id` sha256 `e1ba1c32…` |
| KD rows | `experiments/115-gse295214-h1/work/compact-donor` (`delta.npy` 2,354 × 18,533) |
| **cohort** | **100 TFs** = R5 symbols ∩ bank `perts.txt`, sorted; `work/COHORT.tsv`, sha256 `fca12254…` |
| **gene axis** | **15,309** canonical genes; `work/AXIS.txt`, sha256 `ec1a5719…` |
| conservation | `C_excl` = mean of 279's `Ccross` over the **338** cross-line pairs **not** involving `gse295214_h1` (`held279.npz`); finite for 97/100 |
| geometry | `gen_pc20` = `persrc290.npz['pc20']` row for source `gse295214_h1` |

**Rows are TF targets; columns are measured genes.** The cohort is a set of
*rows* on both sides; the axis is a set of *columns*.

**Axis definition.** A canonical gene enters the axis iff it is measured in
**every** cohort Nakatake row (after 115's duplicate-collapse rule: a canonical
gene is measured only if *every* source column mapping onto it is measured) **and**
in **every** cohort bank row. Both masks are row-constant (measured, §7 C4), so
the axis is the same for every TF and **no padded zero is ever treated as an
expression value**. Coordinates outside a row's support are exactly 0 in 115's
convention and are excluded by construction, not by imputation.

## 4. The objects (ONE computation from the single access)

Both sides are put in the bank's convention by 115's transform, verbatim
(reproduced bit-identically on `Pan2025_hESC`, §6 KA-L1):

```
perturbed_log = layers/control + X
*_cpm         = expm1(max(*_log, 0))                     # sums to 10,000
collapse      = sum duplicate source symbols in LINEAR cpm space
delta         = log1p(5*collapsed_perturbed) - log1p(5*collapsed_control)
```

- `OE_i` = the Nakatake row of TF *i*, transformed, restricted to the axis.
- `KD_i` = the bank `delta` row of TF *i*, restricted to the axis.
- **Primary (CO-CENTRED):** `OEr = OE - mean_{cohort} OE`, `KDr = KD - mean_{cohort} KD`
  — **both sides centred on the SAME 100 TFs**.
- `S_i = cos(OEr_i, KDr_i)` with the own-transcript coordinates of **both**
  members of the pair removed (98/100 cohort TFs have an own column on the axis;
  a shuffled pair drops two coordinates, a matched pair one — declared, so the
  matched and shuffled statistics are computed by the same code path).

**Labelled sensitivities, never a verdict:**
(S1) bank-wide centring — OE on all 510 Nakatake rows, KD on all 2,354 bank rows
(the interrupted draft's convention);
(S2) **uncentred**, control-referenced `S⁰_i = cos(OE_i, KD_i)` — so that any
conclusion created by centring can be recognised;
(S3) cohort-centred but own column retained.

## 5. The four verdict cells

Let `N` be the **target-shuffled pairing null**: 2,000 derangements of the
OE→KD assignment, each recomputed through the identical code path. No statistic
is tested against 0; each is tested against `N`.

| cell | question | statistic | rule (one-sided) |
|---|---|---|---|
| **C7a** | antisymmetry exists | median `S`, and tail asymmetry `share(S<−0.05) − share(S>+0.05)` | median `S` **below** `N`'s medians (p_less) **AND** tail asymmetry **above** `N`'s (p_greater); PASS iff **both** survive Holm at FWER 0.05 |
| **C7b** | conservation predicts mirroring | Spearman(`C_excl`, −`S`), n = 97 | **above** `N`'s (p_greater), Holm |
| **C7c** | geometry predicts the SIGN | Spearman(`gen_pc20`, −`S`) | **above** `N`'s (p_greater), Holm |
| **C7d** | geometry predicts the SIZE | Spearman(`gen_pc20`, \|`S`\|) | **above** `N`'s (p_greater), Holm |

**Multiplicity.** Holm step-down over the four cells at familywise 0.05; C7a's
two conditions are a conjunction inside one cell, entered at the larger of its
two p-values. Unadjusted paired **TF-bootstrap** 95 % intervals (20,000 draws)
are reported for every cell as **descriptive** readings, never as the rule.

**Why the null and not zero.** A centred cosine over 100 correlated rows has no
guaranteed zero median (302's mock showed the same), and `gen_pc20` is computed
**on the KD row itself**, so it is coupled to the KD side of every cosine. In `N`
the KD rows, the KD features, `C_excl` and the centring are all **held fixed**
and only the OE→target matching is randomised: `N` therefore measures exactly
"how much of this could a generic chance alignment produce?". Partialling out
magnitudes does not remove that coupling and is reported separately.

## 6. Known-answer rows

**LOADER rows — already run, all PASS (`work/ka303.json`); a failure is an
ABORT, not a VOID.**

| row | reading | verdict |
|---|---|---|
| KA-L1 | `delta.npy` for 20 seeded `Pan2025_hESC` targets through this loader: max abs diff **0.0e+00**, bit-identical (tol 1e-6, justified by float32 storage) | **PASS** |
| KA-L2 | `b_ctrl.npy` from the first Pan2025 control row: max abs diff **0.0e+00** | **PASS** |
| KA-L3 | `measured_mask.npy` rows identical for the same 20 | **PASS** |
| KA-L4 | the bank's `rows.tsv` source-row map and row order reproduce from `obs` alone (2,354 rows) | **PASS** |
| KA-L5 | Pan2025 control rows identical to the reference row on the sample | **PASS** |
| KA-C | 279's `Ccross` for the **26** pairs involving `gse295214_h1`, **15,047** (pair, gene) cosines: max abs diff **0.0e+00** (reproduces 302's KA-C to the digit) | **PASS** |
| KA-W | importing 276 / 279 wrote nothing outside `experiments/303-*/work` (0 files changed) | **PASS** |
| KA-M | cohorts 100 / 45 / 13 vs the coordinator's 100 / 58 / 13 — **one discrepancy, reported, not forced** (§7 C2) | **REPORTED** |

**DATA rows — computed at access time; they CONDITION THE INTERPRETATION and
never the verdict** (302's KA-E is the reason this separation is explicit):

- **KA-D1 (CPM gate).** Every Nakatake row's reconstructed control and perturbed
  CPM10k sums lie within ±1 of 10,000. *Expected to hold from the producer's code
  (§7 C5); if it fails, the transform's premise is wrong and the verdicts are
  reported as **VOID**, exactly as 302's KA-E clause fired.*
- **KA-D2 (shared control).** All 510 Nakatake `layers/control` rows are
  identical. *Expected: the producer used a single `Emerald` GFP sample as the
  control for every row (§7 C5). If false, the shared-control reading in §8 is
  withdrawn; the verdicts stand.*
- **KA-D3 (transgene detected).** The overexpressed TF's own transcript is up in
  its own Nakatake row: share positive, median own-coordinate. Reported next to
  the CRISPRi side's own coordinate (302 read +0.337 / 95.1 % and −0.700 / 100 %
  on its libraries). *Descriptive; a low share bounds interpretation, it does not
  void a verdict, because a PiggyBac transgene need not map to the endogenous
  probe (§7 U4).*
- **KA-D4 (no-signal row).** The distribution of `N` printed next to the observed
  statistic for every cell.
- **KA-D5 (finiteness/type).** `delta` all finite, float32; mask boolean and
  row-constant; axis order equal to `work/AXIS.txt`.

## 7. Declared secondaries (reported whatever they show, no verdict)

1. **Cross-library OE agreement on the 45 TFs** that also carry 302's `S_g`
   (the coordinator's "58" is against 302's 175-row *scored* set; only 45 are in
   302's 148-row *primary* cohort — §C2): `cos(Joung ORF residual rebuilt from
   `sums302.npz`, Nakatake residual)` and Spearman(`S`, 302's `S_g`).
   **DESCRIPTIVE and conditional on 302's unresolved processing assumptions.**
   Two overexpression libraries in **different lines** (H1 vs SEES3), different
   modalities (pooled scRNA ORF vs arrayed bulk transgene) and different
   durations are **not** repeat measurements of one intervention: this agreement
   is **never** used as a plug-in reliability, a disattenuation factor or a
   construct-noise ceiling.
2. `S` restricted to TFs whose own transcript is up (KA-D3), and to the 13 TFs
   the database also carries as Joung rows.
3. Baselines for C7b–C7d, Spearman against the same target: KD magnitude, OE
   magnitude, `C_excl` vs coupled `C_held`, `r_bar`, `breadth`, `mag` (290).
   Partial Spearman given both magnitudes reported alongside.
4. **Family-aware robustness**: an HGNC gene-group block bootstrap over the
   cohort. With n = 100 and mostly singleton groups this **does not support a
   formal family-held adjustment**; it is reported as a robustness reading and
   labelled as such.
5. The three centring sensitivities S1–S3 of §4.
6. Post-hoc tails: top 15 each way, tail shares at ±0.10 — labelled post-hoc.

## 8. Deliverables

`work/nak303.npz` (the sufficient arrays, hashed, written once),
`work/symmetry303.csv` (per TF: coverage, own-transcript delta, `S`, `|S|`,
`C_excl`, `C_held`, KD mag, OE mag, `gen_pc20`, `breadth`, 302's `S_g` where
present), `RESULTS.md`, and one line reporting the four cells with their
permutation p-values, Holm decisions, bootstrap intervals and the null.

Nothing competition-facing is selected, tuned or proposed. Panel overlap of the
cohort is reported for the record only.

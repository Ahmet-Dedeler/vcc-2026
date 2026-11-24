# 303 — PREACCESS report (reserve R5 UNREAD)

Written 2026-09-07 by the 303 worker under `HANDOFF-MANDATE.md`.
**No byte of `X` or `layers/control` for the 510 `Nakatake2020_hESC_TF_screen`
rows has been read.** Reads performed: `obs` / `var` / `uns`,
`layers/measure_mask` (free per the R5 declaration), and `X` / `layers/control`
for **20 seeded `Pan2025_hESC` rows + 1 control row** (explicitly free — the
already-used bank source). A monkeypatched `h5py.Dataset.__getitem__` spy
confirms the metadata path touches **0** reserve value slices (`test303.py` H3).

Artefacts: `work/meta303.json`, `work/R5-ROWS.tsv`, `work/COHORT.tsv`,
`work/AXIS.txt`, `work/ka303.json`, `work/test303.json`,
`work/producer_hESC_TF_screen_source.py`, `work/source.sha256`.

---

## 1. What exists, and what was never launched

At takeover the directory held only `MANDATE.md` (an unexecuted draft) and
`HANDOFF-MANDATE.md`. There was **no** preregistration, no code, no `work/`, no
access marker, no C7 ledger entry and no release. The coordinator's recovery note
in `docs/CONFIRMATION-LEDGER.md` (R5, 2026-09-07) says the same. **Nothing had
been run; nothing had to be undone.** This phase created the first artefacts.

The reserve is intact: `access303.py --dry-run` reaches the release step and
stops (`work/ACCESS-MARKER.json` absent, `work/RELEASE.json` absent, refusal
code 2 without `--dry-run`).

---

## 2. CONFIRMED metadata (measured, hashed)

| fact | reading |
|---|---|
| source | 4,425,260,319 B, sha256 `0754d6206fe045744dfd9b21eb37c8e3435152c76b158b9810a84c64a029c0ac` — **equal to 115's pin** in `compact-donor/meta.json` |
| shape | `X`, `layers/control`, `layers/measure_mask` all dense float32 **20,345 × 18,108**; `uns` **empty**; `obsm`/`varm`/`obsp`/`varp` empty |
| obs columns | `condition`, `perturbed_gene` (ENSG), `perturbation_sign`, `perturbed_gene_name`, `dataset`; index `id` |
| **R5 rows** | **510**, `dataset == 'Nakatake2020_hESC_TF_screen'`, **contiguous 184…693**, all `condition == 'hESC_TF_screen'`, all `perturbation_sign == +1`, `obs.id == 'hESC_TF_screen_<SYMBOL>'` for all 510 |
| **R5 row identity** | ordered row-index sha256 `85abd9a4c04591313ed66895157cfb9eab1eaf78b1c0fc4f974b5981f871f33b`; symbols sha256 `a76bb7ec0ff7aea94e7068b9176af215b2554ceca5ecba5240dc1e949a2fddb4`; `obs.id` sha256 `e1ba1c32…`; ENSG sha256 `c8481175…`; **510 distinct symbols, 0 duplicates** |
| **R5 measure_mask** | **row-constant**; **15,794 / 18,108 source columns = 0.87221** (the declaration's "0.872" reproduces); after 115's duplicate-collapse rule, **15,324** canonical coordinates, also row-constant |
| bank rows in the file | 2,354 `Pan2025_hESC` rows, contiguous 17,991…20,344; ordered index sha256 `77ff6348…`; mask row-constant on the canonical axis (16,268) |
| **cohort (rows)** | **100** = R5 symbols ∩ bank `perts.txt`; sha256 `fca12254412534255a618dc07aa9f2da3fbdcef05b4c00de6fe900efed0a07ff` — **the coordinator's 100 reproduces exactly** |
| **db-Joung overlap** | **13** of the database's own 28 `Joung2023_hESC_TF_screen` rows — **reproduces exactly**: `CDX2 EOMES HNF4A KLF1 KLF2 KLF4 KLF5 MSC MSGN1 PAX2 PAX5 PTF1A TSHZ1` |
| **302 overlap** | **45**, not 58 — see §3 D1 |
| **gene axis (columns)** | **15,309** canonical genes = measured in every cohort R5 row **and** every cohort bank row; sha256 `ec1a5719fc5a0bb624e15c13d1a3787a4d89d9c59af31c3dc5e447ef7ecdda20` |
| own-target column | present on the axis for **98 / 100** cohort TFs |
| competition exposure | cohort ∩ the 2026 panel = **1** (`TEAD4`); all 510 R5 rows ∩ panel = 13. Nothing is selected, tuned or proposed for the model |

**Loader known-answer rows, all PASS** (`work/ka303.json`, peak RSS 4.06 GiB, 22 s):
`delta.npy` for 20 seeded Pan2025 targets reproduces **bit-identically**
(max abs diff **0.0e+00**, frozen tolerance 1e-6 justified by float32 storage);
`b_ctrl.npy` **0.0e+00**; `measured_mask.npy` rows identical; the bank's
`rows.tsv` source-row map and row order reproduce from `obs` alone for all 2,354
rows; Pan2025 control rows identical to the reference row on the sample; **279's
`Ccross` for the 26 pairs involving `gse295214_h1`, 15,047 (pair, gene) cosines,
max abs diff 0.0e+00** — reproducing 302's KA-C to the digit. Importing 276/279
changed **0** files outside `experiments/303-*/work` (measured by a size+mtime
snapshot, not asserted).

**These are LOADER REGRESSION rows, not biological validation.** They establish
that this code reproduces 115's and 279's arithmetic; they say nothing about
whether the Nakatake rows mean what the draft assumed.

---

## 3. UNCONFIRMED / CORRECTED provenance

The draft's premises were checked against the **producer's own source code**
(`xingjiepan/SCMG_scripts` @ `aa848b9d…`, notebook
`manuscript2025/datasets/pseudobulk_perturbation/hESC_TF_screen.ipynb`, sha256
`2461207bda93b43b15e2b35ccf8e77918d70c9d0cf98c78c72afeeddb2cda455`, **outputs
stripped**, saved as `work/producer_hESC_TF_screen_source.py`) and against the
primary paper (Nakatake et al., *Cell Reports* 31:107655, 2020).

**C1 — the matrix meaning is now established from the producer's code.**
The producer reads a gene × sample CSV of **log10** values, un-logs it
(`10 ** df`), renormalises each **sample** to CPM10k, drops genes not detected in
every sample (`min > 1e-6`), maps symbols → ENSG keeping the **first** of each
duplicate, renames four columns (`H2AFZ→H2AZ1`, `HIST2H3C→H3C14`, `SSX6→SSX6P`,
`T→TBXT`), then places each column on the standard 18,108-gene axis with
`NaN → 0`, renormalises to CPM10k **again** over that axis, takes `log1p`, and
stores `X = log1p(CPM10k perturbed) − log1p(CPM10k control)` with
`layers/control = log1p(CPM10k control)` and
`layers/measure_mask = ~isna(perturb)`. **This is the identical representation
115 assumes for `Pan2025_hESC`** — so 115's transform is applicable, and the CPM
gate (KA-D1) should hold by construction. Log base: the *stored* values are
natural-log `log1p`; the *input file* was base-10 and is un-logged first.

**C2 — the control is ONE sample, shared by all 510 rows.**
`control_exp = tf_screen_df['Emerald']` — a single Emerald-GFP transgene sample.
Consequences: (i) `layers/control` should be identical across all 510 rows
(KA-D2); (ii) **the control carries no replication**, so its measurement noise is
a component common to every row — largely removed by the cohort centring, fully
present in the uncentred secondary S2; (iii) **510 = the paper's 511 transgene
inductions minus Emerald** — consistent, and the count is not forced.

**C3 — the cell line is NOT the bank's line. This falsifies the draft.**
Nakatake et al. used **SEES3** hESCs (male; Akutsu et al. 2015) with a
Dox-inducible PiggyBac vector in a common rtTA3G host line. The `gse295214_h1`
bank is **H1** hESC CRISPRi Perturb-seq (115's `SOURCE_MANIFEST.json`: "H1 hESC",
"CRISPRi Perturb-seq, 10x 5-prime"). The draft's "the same line (hESC, the map's
`gse295214_h1` bank)" is **wrong**: "hESC" is not an identity. C7 is a
**cross-study, cross-line** comparison and is stated as such.

**C4 — the assay is not established as RNA-seq.**
The paper reports "2,174 transcriptome datasets — RNA-seq **or** microarrays" at
48 h of Dox. The database's input CSV (`TF_orf_screen/TF_hESC_screen.csv`) is not
in the repository and carries no platform field; `uns` is empty. **Whether the
510 profiles are RNA-seq, microarray, or a mixture is UNKNOWN**, and the log10
input is as consistent with a microarray table as with a log-TPM table. The
draft's "profiled by RNA-seq" is a hypothesis.

**C5 — duration and modality differ.** Nakatake: 48 h transgene induction,
arrayed, bulk, one profile per TF. Pan2025: pooled CRISPRi Perturb-seq
pseudobulk over ~1.3 M cells. Different perturbation direction, different
duration, different measurement unit.

**D1 — the "58" cohort does not reproduce; **45** does.**
`Nakatake ∩ 302's PRIMARY cohort (T100, n = 148)` = **45**.
`Nakatake ∩ 302's full scored set (T, n = 175)` = **58**.
The coordinator's 58 is against the *scored* set; only **45** TFs carry 302's
primary-cohort `S_g`. Reported, not forced: the cross-library secondary is
declared on **45** with the 58 stated.

---

## 4. Design amendments made BEFORE any outcome was viewed

These are corrections to an interrupted draft, not post-hoc repairs; no C7
statistic has been computed on any R5 value.

| # | amendment | why |
|---|---|---|
| A1 | estimand restated as **observed cross-study, cross-line** alignment; no latent-equivalence claim from a null | C3, C4, C5 |
| A2 | **no disattenuation, no reliability ceiling, no "true mirror" bound** | the file has one profile per TF and one control sample: **no replicate structure exists** to estimate an OE-side reliability (unlike 302, which had cells to shard). 302's ceiling arithmetic is not transferable |
| A3 | both sides **co-centred on the same 100-TF cohort**; the bank-wide (510-row / 2,354-row) convention demoted to a labelled sensitivity; the **uncentred** control-referenced cosine added as a secondary | the draft would have compared a 510-row-centred OE vector to a 2,354-row-centred KD vector; centring can create the effect being tested |
| A4 | axis = coordinates **genuinely measured on both sides for every cohort row** (15,309); no padded zero used as expression | 115 writes exact 0 outside a row's support |
| A5 | own-target coordinate removed for **both** members of every pair, matched and shuffled | otherwise the null keeps a coordinate the observed statistic drops |
| A6 | every cell tested against a **target-shuffled pairing null** (2,000 derangements), not against 0; the null's median is measured, not assumed | a centred cosine over correlated rows has no guaranteed zero median |
| A7 | **Holm** familywise 0.05 across C7a–C7d; unadjusted bootstrap intervals reported as descriptive | four verdict cells |
| A8 | conservation feature = **`C_excl`**, 279's `Ccross` averaged over the **338** cross-line pairs that exclude `gse295214_h1` (finite for 97/100 cohort TFs; Spearman with the coupled `C_held` +0.976 on the cohort). The coupled `C_held` is reported beside it, labelled | C_held shares data with the KD side of every cosine |
| A9 | C7d's null **preserves the KD rows and their `gen_pc20`** and randomises only OE→target matching | `gen_pc20` is measured on the KD row itself; partialling magnitudes does not remove that coupling |
| A10 | the 302 cross-library comparison is **descriptive**, on 45 TFs, and is never a plug-in reliability or a construct-noise ceiling | different lines, modalities and durations; 302's own processing assumptions are unresolved (its KA-E VOID) |
| A11 | LOADER KA rows **abort** before access; DATA KA rows **condition interpretation** and can VOID but never silently pass | 302's KA-E froze a property its file did not have |
| A12 | family-aware robustness declared as an HGNC block bootstrap and explicitly stated **not** to support a formal family-held adjustment at n = 100 | rule 5 / mandate item 6 |

---

## 5. Unresolved risks (recorded, not resolved)

- **U1 — platform unknown** (C4). If the 510 profiles are microarray, dynamic
  range compression and probe-level effects apply and the "bulk depth removes the
  ceiling" premise fails for a second reason. Not resolvable from this file.
- **U2 — line difference is a confound, not a nuisance** (C3). A null S is
  compatible with a real mirror abolished by SEES3-vs-H1 biology. C7a's negative
  would be a statement about *observed* responses only.
- **U3 — no OE-side reliability is estimable.** One profile per TF, one control
  sample. Nothing in this design bounds attenuation. Any claim of the form "a
  true mirror of order X is excluded" is **not available** here.
- **U4 — the transgene need not map to the endogenous probe.** KA-D3's own-gene
  check may read low for reasons of construct design (codon-optimised ORF, 3′-UTR
  loss), so a low share bounds interpretation without voiding a verdict.
- **U5 — independence traced across every reader on disk, still not a proof.**
  Every file in the repository that names this h5ad was inspected: 115's
  `acquire.py` (download) and `extract_compact_donor.py` (selects
  `dataset == 'Pan2025_hESC'`, refuses if the count is not 2,354); 115's three
  manifests/receipts; 117 (reads only the `compact-donor` artefacts, never the
  h5ad); and `E166`/`E218`'s `census.json`, which recorded shape, `obs` columns
  and gene-axis overlap only — **no value read**. 276/279's `CANDIDATES` registry
  reaches this source exclusively through the `gse295214_h1` compact-donor bank.
  So no loader on disk reads the Nakatake rows. This is a trace of the code that
  exists, not a proof that no earlier run ever touched them.
- **U6 — one control sample is a single point of failure.** If the Emerald
  profile is atypical, every row inherits the same distortion; the cohort
  centring removes its *shared* part, not its interaction with per-TF support.
- **U7 — 45, not 58** (D1) reduces the cross-library secondary's power.
- **U8 — the 100-TF cohort was fixed by metadata alone**, but the axis (15,309)
  depends on both masks; a coordinator who changes the cohort changes the axis
  hash and invalidates any release bound to it.

---

## 6. The fail-closed execution path (built, tested, unarmed)

`access303.py` is the only code permitted to open the reserve. It refuses unless
`work/RELEASE.json` exists **and** binds, reproducing live: the source sha256,
the R5 ordered row-index sha256, the R5 symbol-list sha256, the axis sha256, the
**final** `PREREGISTRATION.md` sha256, and the sha256 of `lib303.py`,
`stats303.py`, `access303.py`. **303 never creates that file and none exists.**
Before the first byte of `X` it writes and `fsync`s `work/ACCESS-MARKER.json`
(O_EXCL); the marker's presence refuses every later run, *including* one where
`work/nak303.npz` is missing — an interrupted read spends the reserve and is not
resumable. The single pass saves the sufficient arrays (`delta` on the frozen
axis, mask, `b_ctrl`, per-row CPM sums, row/symbol/axis lists) with hashes to
`work/nak303.npz`; every statistic afterwards reads only that file.

**18 / 18 tests pass** (`work/test303.json`, 0.4 s):
synthetic **zero** (median −0.0008), **exact mirror** (−1.000000), **exact
alignment** (+1.000000), a **planted −0.30 mirror recovered at −0.2999**, the
shuffled null measured (median +0.00113, sd 0.00679) rather than assumed, and
symmetric own-column removal; refusals for **wrong axis**, **non-binary mask**,
**nonfinite values**, **CPM10k sum drift**, **a value read of a reserve
dataset** (with mask-only reads still allowed), **no release**, **a tampered
release hash**, and **a second access after a marker**; determinism under a fixed
seed; Holm behaving over four cells; and the preflight-touches-no-reserve spy.

---

## 7. Verdict for the coordinator

**READY_FOR_REVIEW.** Every precondition this phase can test passes; the two
substantive findings are that the reserve is in a **different cell line from the
knockdown bank** and that the **58-TF cross-library cohort is 45**, both of which
change the claim rather than block it. The reserve stays closed: no release JSON
exists, and none may be created by 303.

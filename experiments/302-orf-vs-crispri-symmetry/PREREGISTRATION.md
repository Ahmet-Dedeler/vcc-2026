# 302 — preregistration of claim C6: does overexpression reverse knockdown in the same cell line?

Frozen 2026-09-07, **before the first byte of `data/external/scperturb/JoungZhang2023_atlas.h5ad`'s `X` was
read**. Everything below was written from `obs` / `var` / `uns` metadata, from 301's id map, and from the on-disk
`gse295214_h1` CRISPRi bank — all free reads under the R4 access rule. `X` is opened exactly once, by
`access302.py`, after this file's sha256 is recorded in `docs/CONFIRMATION-LEDGER.md`.

Reserve R4 = the single file `data/external/scperturb/JoungZhang2023_atlas.h5ad`, sha256
`8bd00eac883b469a1e82a2beca18f647e8023904920f233dbb8568e13503fc22`, 5,804,454,475 B, 1,145,823 hESC cells ×
37,528 genes, ORF overexpression. One computation; the result is appended to the ledger whatever it is; the file
is then SPENT. R1–R3 are spent and are not touched. Nothing competition-facing: panel coverage is 0/300.

---

## 1. The question

Every measurement in this programme so far compares the SAME perturbation modality (CRISPRi knockdown) across
contexts. R4 allows the reverse question inside ONE context: **in hESCs, is a transcription factor's
overexpression response the reverse of its knockdown response, and what predicts how strongly it mirrors?**

Joung et al. 2023 overexpressed TF ORFs in hESCs (SHARE-seq, 1,145,823 cells); the map's `gse295214_h1` bank
(Pan et al. 2025, CRISPRi in H1 hESCs, 2,354 perturbations) knocked down 249 of the same TFs. The programme's
prediction (F9, F13) is that the mirrored TFs are the conserved, core-machinery-like, spectrally central
responses, and that recipient-specific responses are not mirrored.

## 2. Objects (all frozen here; one computation)

**Gene axis `A` — 15,710 genes.** The genes that are (i) in the H1 bank's `genes.txt` (18,533), (ii) flagged in
its `measured_mask.npy` (16,268; the mask is row-constant — verified min = max = 16,268 over all 2,354 rows, and
`delta` is exactly 0.0 on unmeasured entries, so they are structural zeros, not measurements), and (iii) present
in the atlas's `var/gene_symbol` (37,528 symbols, all distinct — 0 duplicates). Order: the bank's `genes.txt`
order. Restricting to the bank's measured genes matters: with the structural zeros left in, each TF's cosine
would be attenuated by a TF-specific factor `|a_measured| / |a_full|`.

**Scored set `T` — 175 TF symbols.** Symbols that (a) 301's `joung_orf_map.csv` (sha256
`a6b60664943abc28f3a6a7e9ef257a881adb097ca0bacc363a25f6ce29323211`, 3,136 / 3,369 codes resolved, 99.4836 % of
cells) resolves, (b) have ≥ 50 ORF cells pooled over all their isoform codes, and (c) appear in the H1 bank's
`perts.txt`. 419 ORF codes, 212,329 cells. The 233 unresolved codes (5,917 cells) are dropped, never guessed.

**Primary cohort `T100` — 148 symbols** (≥ 100 pooled ORF cells). Secondary cohort: the other 27 (50–99 cells).

**ORF profile.** For a symbol `s`, `S_s` = the summed raw counts over its cells restricted to `A`, and `L_s` =
the summed total counts of the same cells over ALL 37,528 columns. `p_s = log1p(5e4 · S_s / L_s)`.
Control profile `p_ctrl` from the pooled GFP (code 3367, 49,539 cells) + mCherry (code 3368, 78,274 cells) cells,
same formula.

**ORF response** `r_s = p_s − p_ctrl`. **ORF residual** `a_s = r_s − mean_{t∈T} r_t` (the map's convention:
minus the bank's mean over its scored perturbations). Note the control profile **cancels exactly**:
`a_s = p_s − mean_{t∈T} p_t`. This is stated in advance because the two control ORFs sit entirely in batch 1
while 94.1 % of the scored cells sit in batch 0 — that confound cannot enter `S_s`, only the descriptive
absolute responses.

**CRISPRi residual.** `d_s` = the H1 bank's `delta.npy` row for `s` restricted to `A`;
`b_s = d_s − mean_{t∈T} d_t` (279's `load_unit` path, same arithmetic).

**Own column.** For each `s` present on `A` (144 of `T100`, 26 of the secondary cohort), the coordinate of `s`
itself is deleted from both `a_s` and `b_s` before normalising.

**The statistic.** `S_s = cos(a_s, b_s)`. Negative = overexpression reverses knockdown.

**Reliability.** Cells are assigned to 16 shards by `numpy.random.default_rng(20260907)`; 8 distinct balanced
8/8 partitions of the shards are drawn with the same generator. For draw `d`, each half's profile is built with
the identical formula from that half's shard sums and residualised by the mean over `T` of the SAME half's
profiles; the own column is removed. `rel_s` = mean over the 8 draws of `cos(a_s^{A,d}, a_s^{B,d})`.

**Covariates, reported next to every cell, never conditioned on (programme rule 9):** ORF cells per TF (T100
median 474, q10 156, q90 4,186, min 100, max 14,339), ORF UMI/cell (cohort median 3,502; controls 2,503), the
ORF isoform count (median 2, max 17), `rel_s`, and the batch-1 fraction (T100 median 0.0000, mean 0.0586).

## 3. The claims

Verdict cells: **three** (C6a, C6b, C6c). Everything else is secondary and labelled. All three are read on
`T100`; the rules are unconditioned, as declared.

**C6a — antisymmetry exists.** median `S_s` over `T100` < 0 with the paired TF-bootstrap 95 % interval (20,000
draws, seed 20260907, TFs resampled with replacement) entirely below 0, **AND** #{`S_s` < −0.05} > #{`S_s` >
+0.05}. **PASS iff both hold.** Null: 200 TF-shuffled pairings (the residual sets are fixed; only the pairing
between `{a_s}` and `{b_s}` is permuted, seed 20260907) — expected median ≈ 0, and it is exactly 0 in
expectation because both sides are centred over the same `T`.

**C6b — conservation predicts mirroring.** Spearman(`C_held`, −`S_s`) over the `T100` TFs with a finite
`C_held` in `experiments/279-znf-state-dependence/work/znf279.csv` (**n = 142**), paired TF bootstrap 10,000
draws. **PASS iff ρ > 0 and the 95 % interval is entirely above 0.** Baselines reported, not in the rule:
knockdown magnitude (290's `mag` for `gse295214_h1`, and `‖d_s‖`), ORF response magnitude `‖r_s‖`, the ORF
reliability `rel_s`, `r_bar` from `znf279.csv`; and the partial Spearman given both magnitudes and `rel_s`.

**C6c — F13's sign in a new modality.** Spearman(`gen_pc20` of the H1 knockdown row, −`S_s`) over the same 142,
where `gen_pc20` is `persrc290.npz`'s `pc20` row for source `gse295214_h1`. **PASS iff ρ > 0 and the 95 %
interval is entirely above 0.** Second reading (secondary): the map-wide aggregate `gen_pc20` over all 28
sources, and `breadth` likewise.

**Declared secondaries** (labelled, no verdict): the same three with |`S_s`| (symmetry of either sign); the
50–99-cell cohort (n = 27); per-isoform `S` for the 70 `T100` TFs with ≥ 2 isoform codes at ≥ 50 cells each, and
the isoform agreement as a reliability reading; the fraction of TFs whose ORF response is POSITIVELY aligned with
their knockdown response, and the top 20 each way described post-hoc and labelled; `S_s` disattenuated by
`sqrt(rel_s)`; the association of `S_s` with each covariate.

## 4. Known-answer rows

Computed and printed **before any claim number**. KA-A … KA-D need no `X` and have already run.

| row | statement | status |
|---|---|---|
| **KA-A** | the atlas's control-code cell counts equal 301's: GFP = code 3367, 49,539 cells; mCherry = code 3368, 78,274 cells | **PASS** (`work/meta302.json`) |
| **KA-B** | the coordinator's metadata check reproduces: 249 TFs in both the resolved ORF map and the H1 bank; 233 of them in 279's 17,106-gene held map; 148 with ≥ 100 ORF cells | **PASS** (249 / 233 / 148) |
| **KA-C** | the H1 bank's rows reproduce 279's `Ccross` for every cross-line pair involving `gse295214_h1` at 0.0e+00 through 279's own loader | **PASS** (26 pairs, 15,047 (pair, gene) cosines, max abs diff 0.0e+00; `work/ka302_C.json`) |
| **KA-D** | the ORF map's TF list is 301's, by sha256 `a6b60664943abc28f3a6a7e9ef257a881adb097ca0bacc363a25f6ce29323211`; 1,647 symbols, 99.4836 % of cells covered | **PASS** |
| **KA-E** (ABORT) | inside the single pass, the computed library size (sum over all 37,528 columns) equals `obs['ncounts']` for every cell read — max abs diff 0 | run in `access302.py` |
| **KA-F** | a controls-vs-controls split is a no-signal row: `v = log1p`-profile(GFP) − `log1p`-profile(mCherry) on `A`, own column removed per TF, reads **median over `T100` of \|cos(v, b_s)\| < 0.05**; the whole distribution and the max are reported | run in `access302.py` |
| **KA-G** (ABORT) | the pipeline's positive control: overexpression raises the TF's own transcript and CRISPRi lowers it. Over the 144 `T100` TFs whose own gene is on `A`: median of the ORF residual's own-gene coordinate **> 0** with **≥ 80 %** positive, and median of the CRISPRi residual's own-gene coordinate **< 0** with **≥ 80 %** negative | run in `access302.py` |
| **KA-H** | `X` holds non-negative integers (max fractional part 0 over every value read) | run in `access302.py` |

If **KA-E** or **KA-G** fails, the claim statistics are still computed (the file is spent by the access either
way) but the three verdicts are recorded as **VOID**, not PASS/FAIL, and the numbers are reported as descriptive
only.

## 5. What this cannot show

(i) Modality and sign: the comparison is exactly the point here, so no re-orientation is applied — but a
negative median is a statement about hESCs and this pair of libraries, not a general law. (ii) One line, one
state, two batches, no donor unit. (iii) Panel coverage 0/300: nothing competition-facing is selected, tuned or
checked, and nothing is proposed for the model. (iv) The CRISPRi side has no cell-level data here, so only the
ORF side carries a measured reliability; `S_s` is attenuated by both sides' noise and the disattenuation is
one-sided and labelled. (v) 233 codes (0.52 % of cells) are unresolved and dropped. (vi) Depth: cells per
symbol span 100 to 14,339 — programme rule 9 is honoured by reporting, not by conditioning the rules, and that
choice is frozen here.

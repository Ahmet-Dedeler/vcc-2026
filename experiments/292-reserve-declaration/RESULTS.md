# 292 — RESERVE R3 declared: 8 untouched human genetic-perturbation datasets on disk, 174 held-map genes pooled, 1 new line, from metadata only

**RESERVE R3 DECLARED: 8 datasets, 8 human genetic-perturbation sets, 174 perturbations in the held map (pooled;
152 in the 7-file core), 1 new line; NOTHING COMPUTED BEYOND METADATA.**

**Date 2026-09-06, 15:38–16:05Z (~27 min of a 2 h box).** Deterministic; nothing random was drawn. **No build,
no `vcc`, no `cell-eval2`, no submission, no download, nothing written under `data/` or `submissions/`.** All
writes are inside `experiments/292-reserve-declaration/` except the R3 entry appended to
`docs/CONFIRMATION-LEDGER.md` (a table row + a section, no existing line changed). Sibling 285 was neither read
nor written. `uv run python` throughout. **No `X`, `layers`, `obsm`, `varm` or `obsp` was opened on any file** —
h5ad reads are `h5py` on `obs`, `var`, `uns` only; no response, cosine, conservation value or feature exists in
this experiment's output. Declaration: `RESERVE.md`; ledger: `docs/CONFIRMATION-LEDGER.md` "Reserve R3".

| known-answer row | required | measured | verdict |
|---|---|---|---|
| **KA-1** every scPerturb file's `n_cells` / `n_labels` reproduce 278's independent obs-only scan | equality | **50/50 files, exact** | PASS |
| **KA-2** 279's 17,106 held genes and 290's feature-store genes are one object | same length, same order | **element-wise identical** — the mandate's overlaps (a) and (c) are one number | PASS |
| **KA-3** the 300 panel targets in `conservation276.csv` equal `data/challenge/pert_counts.csv` | equality | **300, equal** | PASS |
| **KA-4** the sealed HCASMC asset is byte-identical to 126's acquisition | `b87452c2…` (RB:296) | **`b87452c2db1f…d9b9`**, 1,531,381,676 B (sha256 only; not decompressed, header not read) | PASS |

---

## 0. Verdict in eight lines

1. **42 of the 104 inventoried datasets on disk are UNUSED** — no code reads them for anything but an obs/var
   census (`work/inventory292.csv`, `evidence` column names the reader for every USED row). Twelve of the 50
   scPerturb files feed banks in 276's source inventory (ten of them among the 28 map banks); the other 38 have
   never had an outcome computed.
2. **Eight of those 38 are human, genetic, transcriptome-wide, deep enough and controlled** → **Reserve R3**:
   Adamson 2016 10X010 (K562, CRISPRi, 82 genes ≥ 50 cells), Datlinger 2017 (Jurkat, Cas9 KO, 28, 3 guides/gene),
   Datlinger 2021 (Jurkat, Cas9 KO, 20 at 1,702 cells), Dixit 2016 (K562, CRISPRi, 10 at 2,435 cells),
   Papalexi 2021 ECCITE (THP-1, Cas9 KO, 23, 4 guides/gene), Papalexi 2021 arrayed (THP-1, 7, one gene per well),
   **Shifrut 2018 (primary human CD8⁺ T cells, Cas9 RNP, 18, 2 donors, 39,813 control cells)**, and — RESTRICTED
   — Tian 2019 day-7 neuron (24; **same study** as the map bank `tian2019_ipsc`). 1,575,574,586 B, 427,073 cells,
   85,859 control cells, sha256 for all eight in `work/sha256.txt`.
3. **Pooled coverage: 181 distinct genes at ≥ 50 cells, 174 of them in 279's 17,106-gene held map, 7 on the
   challenge panel** (152 / 6 for the 7-file core without Tian). 212 (gene × dataset) rows.
4. **One new line.** K562, Jurkat, THP-1 and iPSC-neuron repeat the map's 21 lines (each still adds an
   independent STUDY, and Papalexi adds a new MODALITY in THP-1, whose map bank `wessels_thp1` is Cas13).
   Shifrut's primary human CD8⁺ T cells are the only context outside the 21.
5. **Comparability is measured, not assumed**: each reserve file's `var` covers **2,636–2,739 of the map's
   2,739 response-axis genes** (Shifrut and Tian cover all 2,739).
6. **The access rule is frozen** (§3 of `RESERVE.md`, copied into the ledger): outcomes once, for a
   preregistered claim written into the ledger first, by a named experiment; metadata free; no exploratory look;
   files spent independently; two riders — independence is a property of the CELLS not of the freeze moment
   (R1's relabelling), and depth-condition the pass rule (programme rule 9; cells/gene span 150–2,435).
7. **What R3 cannot do** is stated in the ledger: it cannot measure conservation *within* itself (only 3 genes —
   `EGR1`, `IRF1`, `SOCS1` — appear at ≥ 50 cells in two distinct reserve studies, 0 in three), half of it is
   Cas9 knockout rather than the map's CRISPRi, per-file n is 7–82 so it supports ONE pooled claim with the file
   as the bootstrap block, and its 6-target panel overlap makes it useless and off-limits for anything
   competition-facing.
8. **The sealed IGVF HCASMC asset is NOT part of R3** and was not touched beyond re-digesting it (KA-4): it
   stays under `docs/DEF-RUNBOOK.md` §5.-3's older, competition-facing rule.

---

## 1. How USED / UNUSED was decided

A dataset is **USED** iff some experiment computed an outcome from it. Raw name-grep does not separate that from
inventory noise: every scPerturb file appears in ≥ 14 files on disk because 045/076/166/218/225/270/278 all
enumerate the directory. The discriminator is *code that reads the file*: `rg` over `*.py`/`*.sh` returns readers
for only 15 of the 50 files, and each was then checked for whether it opens `X`:

- **builds a response table (USED):** `scripts/build_donor_table.py` (Replogle ×3), `E166/build_extra.py`
  (Frangieh), `E225/build225.py` lines 60–67 (Sunshine, Tian 2021 ×2, Tian 2019 iPSC, Xu, Norman, Gasperini
  lowMOI, Wessels). Twelve files → twelve banks, all inside 276's inventory; ten of them are among the 28 map
  banks of `E279/build279.json`.
- **metadata only (does NOT spend the data):** `E166/shared.py` `labels()` reads one obs column and nothing else
  — this is the whole of 166's contact with Datlinger 2017, Dixit 2016 and Papalexi ECCITE (166 built donor
  tables for `nadig_hepg2` and `frangieh_melanoma` only, R166's own table); `E218/census166_copy.py` states in
  its docstring that it reads obs codes and the var index through h5py and never opens `X`; `E225/census225.py`
  counts panel overlap from obs; E270's scan and `E278/scan278_h5ad.py` are obs-only by construction.
  `E278/build278.py` contains no `.X`, `read_h5ad`, `h5py` or `layers` reference at all (grep returns nothing) —
  it consumes 278's obs counts.

So the three files that looked most at risk of being spent (Datlinger 2017, Dixit 2016, Papalexi ECCITE) have had
their label vocabulary counted and nothing else. That is the property R3 protects.

## 2. Resolving labels to genes (metadata, frozen grammars)

`work/parse292.py` carries one frozen grammar per dataset, written after reading label vocabularies only:
`<GENE>_p<plasmid>` (Adamson), `Tcrlibrary_<GENE>_<n>` (Datlinger 2017), `<GENE>_<n>` (Datlinger 2021),
`p-sg<GENE>-<n>` (Dixit), `<GENE>g<n>` (Papalexi pooled), plain `<GENE>` (Papalexi arrayed, Tian),
`ES.sg<n>.<GENE>` (Shifrut), `sg<GENE>[-<n>]` (McFarland), `<GENE>_<strand>_<TSS>-P1P2` (TAP-seq), and explicit
NOT-RESOLVABLE for Joung (numeric ORF ids), Gasperini at-scale/high-MOI (guide-combination strings) and Xie
(chromosome coordinates). A resolved symbol counts iff it is in 279's 17,106 held genes ∪ the 300 panel; cells
for a gene are summed over its guides. Everything downstream of the h5ad scan runs off 278's and 292's JSON
artefacts, so it re-runs at < 0.1 GiB.

Two label traps this caught, both worth recording: **`CuiHacohen2023` resolves 31 "genes"** (`IL2`, `IFNG`, …)
but is a MOUSE cytokine-treatment experiment — symbol matching alone would have put it in the reserve; and
**`GasperiniShendure2019_atscale`'s 93 single-token labels** are 4–6-cell residue of a high-MOI enhancer screen,
not gene perturbations.

## 3. Why the other 30 unused files are not in the reserve

10 mouse · 7 chemical · 3 not a transcriptome (Frangieh protein — which is also the SAME cells as a used bank —
and the two Papalexi protein files, 4–24 features) · 3 enhancer/non-gene targets · 3 too small or too shallow
(Adamson 10X001 7 genes, Adamson 10X005 median 3 cells/gene, McFarland 1 gene) · 2 identity unresolvable on disk
· **2 targeted readout, excluded on a measured number**: the two Schraivogel TAP-seq screens are real K562
CRISPRi with 9–10 gene targets at ≥ 50 cells and 4 guides each, but their 3,185 / 4,191 measured genes cover only
**30 / 28** of the map's 2,739 response-axis genes, so no comparable response vector exists. Full table with
reasons: `RESERVE.md` §5.2.

**One recoverable exclusion.** `JoungZhang2023_atlas.h5ad` (sha256 `8bd00eac883b…fc22`, 5.8 GB) is **1,145,823
hESC cells × 3,369 TF ORFs** — the largest untouched human genetic set on disk — but `obs.perturbation` is a bare
numeric ORF id that no `obs`/`var`/`uns` column resolves, and its modality is ORF overexpression, not the map's
knockdown. A ≤ 100 KB supplementary-table download would resolve the ids and would not spend a reserve; this
mandate forbids downloads, so it is recorded, not opened.

## 4. Suggested first claim (not preregistered here — the coordinator writes it)

One pooled, out-of-map generalisation test of both live claims at once over the 174 genes: per reserve gene,
compute its response and its agreement with the map's consensus response for that gene in the 28 banks' other
lines; test (i) that 279's `C_held` predicts that agreement (F8 off-map, 291's C4 shape) and (ii) that 290's
`gen_pc20` predicts it with F13's sign; baselines `M_ann`, magnitude, the gene's own DE count, a shuffled-gene
null; file as the bootstrap block; cells-per-gene a declared covariate; modality a declared stratum. One
preregistration, one access, all eight files. Anything smaller wastes it.

## 5. Process breach, recorded

The mandate capped peak RSS at ~6 GB. `work/scan292.py` peaked at **11.5 GiB** (`work/scan292_ka.json`) because
it materialised whole obs string columns as NumPy object arrays for the two largest files (Joung atlas,
1,145,823 cells; Gasperini at-scale, 203,496 distinct labels). It completed, `memory_pressure` read 752,283 free
pages afterwards and no sibling job failed, but the cap was exceeded; if that scan is ever re-run it must go
per-file in a subprocess. Every later step peaks under 0.1 GiB. No other rule was bent: no build, no `vcc`, no
`cell-eval2`, no download, no write outside the experiment directory except the ledger append.

## 6. Files

`RESERVE.md` (the declaration, the inventory summary, the access rule, what it cannot confirm) ·
`work/inventory292.csv` (104 rows over every dataset on disk, USED/UNUSED + evidence) · `work/sha256.txt`
(55 digests) · `work/reserve292.json` (the eight files' paths, digests, sizes) · `work/scan292.py`,
`work/scan292.json`, `work/scan292_ka.json` (the obs/var/uns scan and its KA rows) · `work/parse292.py`,
`work/perts292.json` (the frozen label grammars and the per-dataset gene resolution) · `work/inventory292.py` ·
`work/varaxis292.json`, `work/varoverlap292.json` (measured-gene overlap with the map's response axis) ·
`work/reserve_totals292.txt` · `work/sha256_scperturb.txt`, `work/sha256_other.txt`.

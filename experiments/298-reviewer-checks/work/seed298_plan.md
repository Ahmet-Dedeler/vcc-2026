# 298 item 4 (T7) — the seed scan, specified but NOT RUN: what is missing and what it would take

**Verdict: NOT TESTABLE ON DISK.** The 2025 library's spacer sequences are not on this machine and
the mandate does not permit downloading them (its one download allowance is a promoter-sequence
table). The match logic is written out below so a future run needs no re-derivation.

## What was checked, and what each source actually holds

| source | what it holds | carries a spacer? |
|---|---|---|
| `data/external/vcc2025/{train,validation,test}/adata_*.h5ad` `obs/guide_id` | 410 categories over the three splits, all of the form `GENE_P1P2_A\|GENE_P1P2_B` — protospacer **names** (`P1`/`P2` = the two protospacers of a target, `_A`/`_B` = the dual-guide construct's two positions). **0/410 are ACGT strings.** `obs` = `[_index, batch, guide_id, target_gene]`; `uns` is empty; `var` carries only gene ids and QC columns. | no |
| `data/external/vcc2025/` (rest) | `gene_names.csv`, three `pert_counts_*.csv` (target, n_cells, median UMI). | no |
| `data/challenge/` | three 2026 context h5ads, `gene_names.csv`, `manifest.json`, `pert_counts.csv`, `vcc_2026_controls.zip`. | no |
| `.../cd4_cell2026/suppl_tables_2026-08-24/sgrna_library_metadata.suppl_table.csv` | the **GWCD4i primary-CD4** CRISPRi library. Different library. | not this one |
| `data/external/gse249595/GSE249595_guides_reference.fa.gz` | GSE249595's guide reference. Different study. | not this one |
| `data/external/myllia-a549/GSE337804_sgRNA_lib.csv.gz` | the Myllia A549 library. Different study. | not this one |
| `data/external/broad16line/metadata/gene_locations_biomart.csv` | 86,402 rows: Ensembl gene id, karyotype band, chromosome, gene start/end (bp), gene name. **Coordinates only — no bases.** | n/a (this is the promoter side) |
| anywhere under `data/` | no `.fa`, `.fasta`, `.2bit` or any genome/promoter **sequence** file. | n/a |

So **both** halves of the scan are missing, and only the promoter half is downloadable under this
mandate. Stopping here is the mandate's own instruction ("if no spacer sequence is on disk, stop
item 4 at the plan and say so").

## The recipe, fully specified (Hartman 2026's, as far as this data can go)

1. **Arms.** The 34 arms of the 17 two-construct targets, in `work/seed298.csv`
   (`target, split, arm, construct, cells, kd_lfc, abs_lfc, failing`), taken from
   `experiments/293-challenge-construct-ceiling/work/construct293.csv`
   (`challenge_20k` regime, one row per target, columns `c_big`/`c_small`, `kd_big`/`kd_small`).
   **Failing = `|kd_lfc| < 1.0`: 13 of 34 arms** (ARPC2_P2, CAST_P1, CAST_P2, DPF2_P2, GRB10_P2,
   HIRA_P1, IGF1R_P2, KDM2B_ENST00000377071.4, KLF10_P1, MED15_P2, PCCB_P2, SMARCB1_P2, ZNF581_P2).
   These 13 are the scan's numerator population; **13**, not 34, is the denominator that a
   "how many failing spacers have a candidate seed off-target" count should be read against.
2. **Spacers (MISSING).** Each `GENE_Pn_A|GENE_Pn_B` id must be resolved to its two 19–20 nt
   protospacers. Source would be the VCC 2025 library table (Arc's release) or the paper's
   supplement. Two spacers per construct → 26 spacers over the 13 failing arms.
3. **Promoter windows.** TSS ± 500 bp per gene. `gene_locations_biomart.csv` gives
   chromosome + gene start/end; the TSS is `start` on the `+` strand and `end` on the `−` strand —
   **the strand column is absent from this file** and must come from the same BioMart query
   (`strand`) or an Ensembl GTF. CRISPRi's active window is roughly TSS −50 … +300, so a
   TSS ± 500 window is the conservative superset.
4. **Promoter sequence (MISSING, downloadable).** The ~20,000 × 1,001 bp windows are ≈ 20 MB of
   bases (≈ 6 MB gzipped) — inside the mandate's 50 MB cap. Ensembl REST/BioMart `sequence` with
   `upstream_flank=500&downstream_flank=500` on the gene set, or a UCSC `hgTables` promoter
   extraction, both public; the run must record the URL, the assembly (GRCh38) and the sha256.
5. **Seed match.** For every failing spacer, take the PAM-proximal seed of length L ∈ {9, 10, 11, 12}
   (the 3′ L nt of the protospacer) and require an exact match of the seed in either orientation
   (the spacer's seed, or its reverse complement) inside a promoter window, with a canonical
   `NGG` PAM immediately 3′ of the matched protospacer-orientation site. Scan with a suffix-array
   or a simple 4^L-bucket hash over the concatenated windows; 26 spacers × 4 lengths × 2 strands is
   trivial compute. Exclude hits in the arm's **own** target's promoter. Record every (spacer,
   gene, L, strand, offset-from-TSS) hit; a candidate off-target = a gene with ≥ 1 hit at the
   longest L that yields any hit.
6. **Expression test.** For each failing construct with ≥ 1 candidate, read that candidate gene's
   own column in the construct's reproducible residual response and test whether it is knocked
   down: the construct's own pseudobulk lfc for the candidate gene vs the distribution of that
   gene's lfc over the other 2025 constructs (a per-gene z or an empirical rank). The response
   object and the residual convention are already built —
   `experiments/293-challenge-construct-ceiling/work/agg293.py` + `lib293.py` (which imports 285's
   `lib285.prof/delta/resid`); the 2025 h5ads carry the cells. **This step needs no download.**
7. **Report.** Count of failing constructs with ≥ 1 candidate seed off-target (of 13), and of those,
   how many show the candidate gene knocked down in their own column. A null: the same scan run on
   the 21 **non-failing** arms' spacers, which should show candidates at the same rate if seed
   off-targeting is not what distinguishes a failing construct.

## What is already known without the scan

- 296 measured the relevant negative directly: failing constructs' responses are reproducible in
  all five public libraries but carry **no shared vector signature**; and the published off-target
  rate for a validated CRISPRi library is 0.8–4.7 % (297, T7). A seed hit in a promoter is
  therefore expected to be rare, and a 13-arm scan is underpowered to estimate that rate — it can
  find an existence proof, not a rate.
- 293's own reading of these 34 arms is that the two-construct term is **knockdown failure**, not
  off-target action: both constructs knocking down gives +0.687 [+0.197, +1.067], the weaker
  failing −0.003.

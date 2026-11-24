# Perturb-ME: scalable mechanism discovery from phenotype-enriched genome-wide screens

Processed data accompanying the manuscript *"Perturb-ME: Scalable mechanism discovery from
phenotype-enriched genome-wide screens"* (Wang, Gu, Frangieh et al.).

A genome-scale CRISPR knockout screen in the A375 melanoma line, FACS-enriched on surface
HLA-A/B/C (bottom and top 5%) and profiled by Perturb-CITE-seq, giving matched single-cell
RNA, surface-protein and sgRNA measurements for phenotype-associated perturbations.

---

## What is and is not in this record

**This record holds processed and derived data**: analysis-ready single-cell objects, the
fitted regulatory model, bulk-screen result tables, and the agentic-interpretation inputs and
outputs. Everything a reader needs to reproduce the figures and re-examine the conclusions.

**Raw sequencing reads are not here.** FASTQ files for all single-cell libraries (gene
expression, CITE-seq, cell hashing, CROP-seq guide capture) and for the bulk sgRNA
enrichment sequencing are deposited separately at the NCBI Gene Expression Omnibus; the
accession will be reported in the paper's Data Availability statement.

---

## Contents

### `01_single_cell/`

Two AnnData (`.h5ad`) objects, `X` stored as a sparse float32 matrix.

| File | Cells | Description |
|---|---|---|
| `PerturbME_A375_HLA_all_cells.h5ad` | 353,091 | all QC-passing cells, all sgRNA multiplicities |
| `PerturbME_A375_HLA_MOI1.h5ad` | 277,152 | sgRNA-singlet cells; the set used for the regulatory model |

**Features** (`var`, 18,683 total): 18,670 gene-level RNA features plus 13 CITE-seq
surface proteins, prefixed `CITE-`. The `modality` column separates them.

> ⚠️ **`X` is not on a single scale.** The RNA block (`var.modality == 'RNA'`) holds **raw
> integer UMI counts**. The 13 CITE-seq columns (`var.modality == 'protein'`) hold
> **IgG-control-normalised values**, `log((antibody UMI + 1) / (isotype control UMI + 1))`,
> as described in Methods — they are already normalised and must not be re-normalised
> together with the RNA counts. The same statement is carried in `adata.uns['X_contents']`.

**Cell annotations** (`obs`):

| Column | Meaning |
|---|---|
| `Condition` | sorted population — `Low`, `High`, `Control` |
| `Sample` | 10x channel (20 channels: 9 HLA-low, 9 HLA-high, 2 control) |
| `MOI` | number of distinct perturbation targets detected in the cell |
| `sgRNA` | assigned guide identity |
| `cell_cycle_phase` | Tirosh et al. G1/S/G2M call |
| `sequencing_run` | 10x run covariate used in the linear model |

**Embeddings** (`obsm`): `X_pca` (50 PCs) and `X_umap`, the coordinates used in Fig. 1c–d.

Loading:

```python
import anndata as ad
adata = ad.read_h5ad('PerturbME_A375_HLA_MOI1.h5ad')
print(adata.uns['X_contents'])                    # what X holds, per modality
rna  = adata[:, adata.var.modality == 'RNA']      # raw integer UMI counts
cite = adata[:, adata.var.modality == 'protein']  # already IgG-normalised
```

The kNN graph and the log-normalised matrix are not stored — both are re-derivable from the
counts in a few lines, and omitting them keeps these files small enough to download.

### `02_regulatory_model/`

| File | Contents |
|---|---|
| `elasticnet_beta_matrix.csv.gz` | EM-refined ElasticNet coefficients, 18,683 features × 231 covariates. This is the matrix behind Fig. 2b and everything downstream. |
| `target_modules.csv` | K-means partition of the β columns into the seven target modules (#0–#6) |
| `gene_programs.csv` | K-means partition of 1,998 responsive RNA features into the nine gene programs (#0–#8) |
| `module_program_edges.csv` | mean β per (module × program) pair — the edges in Fig. 2c and the Fig. 2e Sankey |
| `target_module_representatives.csv`, `gene_program_representatives.csv` | the eight representative genes shown per module/program in Fig. 2c |
| `Fig2d_curated_gene_reconciliation.csv` | see the note below |

Covariates in the β matrix: 221 gene-target perturbations, 2 pooled control categories
(`HIGH_CTL`, `LOW_CTL`) and 8 technical or biological covariates (cell-cycle `G1`/`G2M`/`S`,
sequencing run `TENX_1`–`3`, `UMI`, `cells_per_target`).

> **Two things worth knowing before you read these tables.**
>
> **1. `target_modules.csv` has 224 rows, not 221.** K-means was run across the β columns,
> which include the two pooled control categories and the `G1` cell-cycle covariate alongside
> the 221 gene targets. The partition is reproduced here exactly as it was fitted, with a
> `covariate_type` column labelling each row, so the 221 gene targets the manuscript refers
> to can be selected directly.
>
> **2. Fig. 2d gene lists are curated, and intentionally differ from Fig. 2b/2c.** The module
> membership in Fig. 2b/2c is purely algorithmic (the K-means partition in
> `target_modules.csv`). The gene lists drawn in Fig. 2d, and the module descriptions in the
> main text, were curated from the Co-Scientist reports and add literature-derived
> regulators that are not members of the fitted model — for example EZH2, SUZ12, PTPN2,
> PTPN11, SOCS1 and IRF9. `Fig2d_curated_gene_reconciliation.csv` classifies every Fig. 2d
> gene as an algorithmic module member, a bulk-screen hit outside the model, or a
> literature-only addition, so the two conventions can be told apart gene by gene.

### `03_bulk_screen/`

MAGeCK v0.5.9.2 output for the bulk sgRNA enrichment screen (9 libraries: 4 input
replicates, 1 unsorted control, 2 HLA-high, 2 HLA-low).

- `mageck_sgRNA_counts.txt`, `mageck_sgRNA_counts_normalized.txt` — guide-level count tables
- `mageck_<comparison>.gene_summary.txt` / `.sgrna_summary.txt` for `Low_vs_Input`,
  `High_vs_Input`, `CTL_vs_Input`, `Low_vs_CTL`, `High_vs_CTL`, `High_vs_Low`

Positive enrichment in HLA-low indicates a perturbation whose knockout *decreases* surface
HLA-I; positive enrichment in HLA-high indicates one whose knockout *increases* it.

### `04_cosci/`

- `CoScientist_prompt.docx` — the exact input given to Google DeepMind's Co-Scientist
- `CoScientist_response.pdf` — its full response, including the hypothesis tournament

### `05_reference/`

- `CITE_antibody_panel.csv` — the 13-antibody CITE-seq panel

### `MD5SUMS.csv`

Checksums for every file in this record.

---

## Reproducing the figures

The full analysis codebase is released at <https://github.com/Genentech/PerturbME> and
archived at Zenodo with its own DOI; both are cited in the paper's Code Availability
statement.

The Fig. 2b/2c module and gene-program partition can be reproduced from
`02_regulatory_model/elasticnet_beta_matrix.csv.gz` alone, without the single-cell objects,
by `notebook/20260519_Figure2B_Rerun_BetaCorr.ipynb` in that repository: retain RNA rows in
which fewer than 25% of coefficients have |β| ≤ 0.05 and covariate columns in which fewer
than 80% do (giving 1,998 features × 224 covariates), then run K-means with k = 7 on the
target-target β correlation and k = 9 on the feature-feature β correlation
(`random_state=3`, `n_init=10`).

Note that the 224 clustered covariates are the 221 gene-target perturbations plus the two
pooled controls and the `G1` cell-cycle covariate; `target_modules.csv` labels each row
accordingly.

## Citation

Please cite the manuscript and this record's DOI.

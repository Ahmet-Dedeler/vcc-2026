# 047 — Direct primary-CD4 donor extraction and evaluation

## Goal

Evaluate the public genome-scale primary CD4+ T-cell CRISPRi atlas as a
direct perturbation-response donor without bringing its 16.8 GB DE object onto
the laptop. The source covers 291/300 VCC targets across at least one of
`Rest`, `Stim8hr`, and `Stim48hr`, making it the highest-value new direct donor
found so far.

This experiment has two stages:

1. On an ephemeral cloud disk, stream only VCC-panel rows and official VCC
   genes from `GWCD4i.DE_stats.h5ad`.
2. Run a small effect-level leave-one-donor-out (LODO) screen from the compact
   outputs. The most informative fold holds out K562, the available
   hematopoietic-lineage analogue.

It does not run `cell-eval2`, emit single-cell counts, fit a full covariance,
or submit anything.

## Source and schema assumptions

Authoritative source:

```text
s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad
```

Expected shape is approximately 33,983 perturbation-condition rows by 10,282
measured genes. Required fields are:

- `.obs['target_contrast_gene_name']`
- `.obs['target_contrast']`
- `.obs['culture_condition']`
- `.var['gene_name']`
- `.var['gene_ids']`
- `.layers['log_fc']`, `zscore`, `lfcSE`, `adj_p_value`, `p_value`, and
  `baseMean`

The source documents `log_fc` as DESeq2 **log2 fold change**. This repository's
donor effects use natural-log units, so extraction defines:

```text
delta = ln(2) * log_fc
```

Non-finite source log fold changes are retained in `log_fc.npy`, marked in
`valid_log_fc.npy`, and replaced by zero only in donor-compatible `delta.npy`.

`b_ctrl.npy` is required by the existing `Donor` file contract, but this DE
object has no control-only pseudobulk. The extractor therefore derives a
normalized proxy from median `baseMean`. It is acceptable for loading and
rough expression gating; it must not be used for basal cell-line matching or
count-noise/empirical-Bayes calibration. `meta.json` records this restriction.

### Symbol mapping

- VCC target rows use exact source symbols first.
- A unique case-insensitive match is allowed second.
- A unique Ensembl-ID match through `.var` is allowed last.
- No fuzzy or hand-curated alias mapping is performed.
- Duplicate measured-gene symbols are collapsed by selecting the source
  feature with the largest median `baseMean` over selected panel rows.
- Duplicate perturbation rows are resolved deterministically by source QC:
  on-target significance, absence of low-expression/off-target flags, guide
  count, then cell count. Every candidate and decision is written to
  `target_row_map.csv`.

## Cloud-only extraction

Use an instance with at least 40 GB free ephemeral disk and 8 GB RAM. The
extractor's default RSS guard is 4 GiB; expected live RSS is far lower because
layers are copied in 32-row blocks.

Example on AWS:

```bash
aws s3 cp --no-sign-request \
  s3://genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad \
  /mnt/vcc/GWCD4i.DE_stats.h5ad

uv run python experiments/047-cd4-direct/extract_cd4.py \
  --input /mnt/vcc/GWCD4i.DE_stats.h5ad \
  --output /mnt/vcc/cd4_compact
```

Copy only `/mnt/vcc/cd4_compact/` back to durable storage. Do not copy the
16.8 GB source object to this Mac.

The three donor directories are:

```text
cd4_rest/
cd4_stim8hr/
cd4_stim48hr/
```

Each has `delta.npy`, `b_ctrl.npy`, `perts.txt`, `genes.txt`, source statistic
matrices, `valid_log_fc.npy`, `obs.csv`, and `meta.json`.

Before transferring outputs, inspect `/mnt/vcc/cd4_compact/manifest.json`.

## Extraction validation gates

The script fails unless all of these hold by default:

- all six documented layers exist and match the `.obs × .var` shape;
- at least 275 VCC targets map in each condition;
- at least 290 VCC targets map in the union of conditions;
- at least 9,500 official VCC genes map;
- at least 98% of selected `log_fc` values are finite;
- perturbation and gene labels are unique;
- every donor-compatible `delta.npy` value is finite;
- `baseMean` can produce a nondegenerate compatibility profile;
- peak process RSS remains below 4 GiB.

Expected values from the public metadata are 283 Rest, 287 Stim8hr, 284
Stim48hr, and 291 targets in their union. Large deviations mean the schema,
source revision, or symbol mapping changed; do not weaken thresholds until the
cause is understood.

## Leak-free LODO screen

After the compact outputs are available locally or on the cloud worker:

```bash
uv run python experiments/047-cd4-direct/evaluate_lodo.py \
  --cd4-root /path/to/cd4_compact
```

The evaluator:

- memory-maps every `delta.npy`;
- uses only VCC perturbations shared by all three existing broad donors and
  all three CD4 conditions;
- defines one fixed gene axis from the official VCC genes shared by every
  bank;
- removes **all 300 VCC target genes before any centering, normalization, or
  combination**, preventing target-coordinate leakage;
- holds out K562-gwps, Orion-HCT116, and Orion-HEK293T in turn;
- compares the existing two-donor consensus with Rest, Stim8hr, Stim48hr, and
  a three-condition CD4 mean at fixed blend fractions;
- reports raw-centered and mean-absolute gene-normalized variants;
- selects one candidate using HCT116 and HEK293T only, then reports its K562
  result as a leakage-free lineage-analogue check.

The realized fixed axis is 6,406 non-panel genes. The original 6,500-gene
sanity threshold was an estimate made before extraction; the default gate is
6,300 so it accepts the measured axis while still failing on a material mapping
loss.

Outputs are `results.json` and `results.csv`. This is PDS-only and
effect-level; a winner still needs the established union-whitening and
count-level validation pipeline before it can affect a submission.

## Decision thresholds

Treat CD4 as a serious production-donor candidate only if:

1. the HCT116/HEK293T-selected configuration improves held-out K562 by at
   least **+0.003 raw PDS**;
2. its three-fold macro delta is positive and no non-K562 fold loses more than
   **0.002**;
3. the direction of the K562 gain agrees between raw and gene-normalized
   screens, or there is a predeclared reason to prefer one;
4. the result survives a second evaluation using the production union
   whitening operator fitted without the held-out donor; and
5. any final full-metric test is run one process at a time under the repository
   memory guard.

Below +0.002 is treated as unresolved because it is in the range of ordinary
fold/seed variation. A K562-only winner selected directly on K562 is
exploratory evidence, not a deployment decision.

## Lightweight checks

The included tests use a tiny synthetic H5AD and never touch the public atlas:

```bash
uv run python experiments/047-cd4-direct/test_cd4_direct.py
```

No VCC submission is made by either script.

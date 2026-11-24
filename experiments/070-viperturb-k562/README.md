# 070 — source-wide VIPerturb K562 assay-balanced slot

Status: **corrected evaluation completed; rejected for production**.

Evaluation-only retry `20260823T070153Z` completed from the preserved compact.
Macro PDS changed by `-0.0069624664`, with both held folds negative and a
bootstrap 95% interval of `[-0.0108373282, -0.0034542003]`. See
`RESULTS.md`; no count-level confirmation or challenge submission is allowed.

## Question

Does an independent K562 CRISPRi assay improve the frozen experiment-064 arm-D
direction when it is added without changing the selected whitening or CD4
configuration?

The official source is Zenodo record `18460279`,
`genome_wide_filtered.rds`:

- expected bytes: `3,611,371,127`;
- expected MD5: `ac52c41233fc7a5a5c3810e0a075007f`;
- manifest: `genome_wide_manifest.txt`;
- expected manifest MD5: `c095c8f2d16136753c22ef5d9b30a22b`.

The manifest contains exactly 124 VCC-panel genes marked
`is_in_filtered_obj=TRUE`. These are retained even when cell count is low.
Rows below 30 cells are explicitly reported from both the manifest and the
loaded Seurat metadata. The current official manifest has 22 such VCC targets,
including one with only 8 cells.

The manifest also repeats its `NO-TARGET` control row once in each of the
three source bins. Those rows describe one control class, not three
perturbations. The extractor removes recognized control rows before enforcing
the unique, frozen 6,724-perturbation source catalogue; duplicate non-control
genes still fail closed.

A no-RDS preflight against the existing donor label files predicts 98 targets
will satisfy the incumbent-K562/PerturbMatch/VIPerturb/full-gallery
intersection, above the 70-target gate. The extractor and evaluator still
recompute this from actual output rather than trusting the preflight.

## Frozen comparison

Valid held folds are Orion HCT116 and Orion HEK293T. The full gallery is the
289 challenge targets present in both held donors.

The baseline is the full-gallery extension of experiment-064 arm D:

- 50:50 incumbent/PerturbMatch K562, with incumbent fallback;
- the non-held Orion donor;
- the existing candidate-independent gamma-`0.65` experiment-062 whitener;
- post-whitening CD4 mean3 `gene_norm`, fraction `0.05`, applied only where all
  three CD4 conditions contain the target.

The candidate changes only targets available in incumbent K562,
PerturbMatch, and filtered VIPerturb. On official non-panel genes shared by
all three K562 views, its K562 slot is:

```text
0.5 * (0.5 * incumbent + 0.5 * PerturbMatch) + 0.5 * VIPerturb
```

Every other target and coordinate keeps the baseline value. Candidate
whitening and CD4 use the baseline arm's fixed preservation scalars, making
unmodified target rows exactly identical after the full pipeline.

## Source-wide extraction and normalization

`extract_viperturb.R` reads the Seurat RDS directly. It:

- identifies the RNA assay and prefers a single counts layer;
- falls back to a single data layer only when counts are absent, recording the
  resulting effect definition;
- detects and audits the perturbation metadata column;
- aggregates controls and **every perturbation marked as present in the
  filtered source object**;
- selects only genes on the official VCC axis;
- computes grouped pseudobulks in 256-perturbation chunks;
- computes the per-gene center, per-gene variance, and median centered norm over
  the complete filtered source catalogue;
- only after freezing those source-wide statistics, retains the 124 VCC-panel
  rows for evaluation;
- emits the raw panel effects, source-wide center/variance/scale, the old
  panel-centered center/scale for one diagnostic ablation, and audit files.

It never writes a MatrixMarket, H5AD, h5Seurat, or whole-cell matrix
intermediate. Python rejects a compact artifact unless its metadata and hashes
prove that normalization used more rows than the retained panel and that every
panel row belongs to the source catalogue.

The primary candidate uses the source-wide normalization. The panel-centered
version reproduces the old normalization mistake as a labeled diagnostic
ablation and can never be selected for shipping.

## Null controls

Both nulls reuse each fold's frozen arm-D output, whitener, CD4 preservation
scale, truth axis, and exact PDS rank projection. They never fit or select on
held truth.

- **Target-label null:** permute VIPerturb rows only among modified targets,
  stratified by deterministic observed-cell-count quartile.
- **Gene-coordinate null:** apply one shared coordinate permutation per
  replicate within deterministic source-variance deciles on the shared
  non-panel triplet genes, before the frozen whitening operator.

The default is 1,000 deterministic replicates per null and is configurable
with `--null-reps`. Nulls run serially/in small bounded chunks and write their
full macro-gain arrays as `.npy` files plus p95 and empirical upper-tail
p-values in `results.json`.

The PDS null acceleration is exact: prediction row norms are omitted only
inside the rank calculation because a positive row-wise denominator cannot
change the ordering against the 289 truth rows. The observed candidate is
checked against the direct full pipeline before nulls run.

## Stopping gate

The candidate advances only if all conditions hold:

1. at least 70 full-gallery targets are actually modified;
2. both held-fold PDS deltas are nonnegative;
3. macro PDS gain is at least `+0.003`;
4. the paired target-cluster bootstrap 95% lower bound is positive.
5. observed macro gain is strictly above the target-label-null p95;
6. observed macro gain is strictly above the gene-coordinate-null p95;
7. after excluding modified targets with fewer than 30 cells in either the
   manifest or extracted metadata, the mean gain remains positive.

The low-cell condition is statistically defined only when at least 20
non-low-cell modified targets remain. Fewer than 20 fails closed. The report
also includes low/non-low mean deltas and the full-gallery delta obtained by
resetting low-cell modified rows to arm D. Failure of any gate means stop.

## Cloud runner

The launcher requires explicit clearance and an already-running EC2 instance.
It never creates or starts an instance. The runner holds the shared
`/work/vcc-heavy.lock`, rejects overlapping heavy processes, requires at least
24 GiB free disk, downloads the RDS only on the worker, checks byte count and
MD5, and runs extraction/evaluation serially under memory guards.

R is installed from conda-forge binaries in a pinned micromamba environment:
R `4.4.3`, Seurat `5.3.1`, and SeuratObject `5.1.0`. This avoids source
compilation and system-library drift. On Amazon Linux the bootstrap installs
only missing utilities one at a time, excludes all curl packages from those
transactions, verifies that an existing `curl-minimal` package remains
installed, and includes the repository root `README.md` in the source archive.

No cloud command should be run without separate approval:

```bash
I_HAVE_EXPLICIT_CLEARANCE=YES \
INSTANCE_ID=<already-running-worker> \
bash experiments/070-viperturb-k562/aws/launch.sh
```

After a successful extraction, an evaluation-only retry can reuse the
preserved compact and all existing inputs without downloading or reading the
RDS:

```bash
I_HAVE_EXPLICIT_CLEARANCE=YES \
INSTANCE_ID=i-08d90ae88c7b210ed \
BASE_RUN_ID=20260823T063300Z \
bash experiments/070-viperturb-k562/aws/evaluate_only.sh
```

The retry stages only the current code, verifies the preserved input paths,
runs the synthetic suite, and writes results beneath a fresh
`evaluation-retry/<timestamp>/` S3 prefix. It does not invoke R, reference the
RDS, or submit anything.

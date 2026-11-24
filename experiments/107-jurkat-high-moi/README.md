# 107 — high-MOI Jurkat target deconvolution

Status: **source identified; acquisition resumed; identifiable modeling
protocol frozen before transcriptome coordinate values are read**

Date: **2026-08-30**

## Question

Can the genome-scale Jurkat CRISPRi CROP-seq screen GSE249595 supply a
high-confidence, narrow response channel for the VCC predictor after jointly
deconvolving the multiple perturbations carried by each cell?

This is motivated by two source-only facts:

1. prior basal-profile work independently identified challenge context A as
   Jurkat-like;
2. the released guide library covers 299 of the 300 VCC perturbation targets,
   whereas the existing Nadig Jurkat donor covers none of them.

Neither fact authorizes direct use. GSE249595 is a high-MOI combinatorial
screen, so averaging all cells carrying one target would confound that target
with several co-perturbations.

## Frozen source

- GEO series: `GSE249595`
- parent series: `GSE247601`
- BioProject: `PRJNA1049794`
- assay: Jurkat E6 `dCas9-KRAB` genome-scale CRISPRi CROP-seq
- raw archive:
  `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249595/suppl/GSE249595_RAW.tar`
- guide reference:
  `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249595/suppl/GSE249595_guides_reference.fa.gz`
- GEO file inventory:
  `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE249nnn/GSE249595/suppl/filelist.txt`
- official supplementary table:
  `https://static-content.springer.com/esm/art%3A10.1038%2Fs41556-025-01626-9/MediaObjects/41556_2025_1626_MOESM3_ESM.xlsx`
- raw archive size: `3,753,820,160` bytes

The archive contains 16 channels, each with transcriptome, guide, and hash
Matrix Market triples. The guide FASTA has 83,401 records, including 250
non-targeting controls. Exact guide-library overlap is 299/300 VCC targets;
`ZNF714` is absent. The source paper reports more than 586,000 retained cells,
a median of 13 called guides per cell, about 400 cells per target, a 10%
untreated / 90% activated pool, deMULTIplex hash assignment, and binomial
guide calling at BH-adjusted `p <= 0.05`.

The small official metadata files were downloaded and hashed before the raw
archive. Supplementary Table 1 contains 18,466 populated rows and 18,465
unique symbols. Table 2 contains 374 gene rows plus one prose footnote in the
symbol column; only rows with a non-null integer Entrez ID are gene records.
Table 6 contains seven hash oligos and explicitly maps them to untreated or
activated conditions. The FASTA parses to 18,601 labels including
controls/composite labels and independently confirms the 299/300 VCC overlap.

## Scientific limitation

The experiment used a very high multiplicity of infection, with a reported
median of 13 guides per retained cell. Effects must therefore be estimated jointly.
Target-carrier pseudobulks, nearest-guide assignment, and arbitrary
single-guide filtering are forbidden.

The transcriptome matrices expose 20,606 features and overlap 18,378 of the
18,533 VCC genes, but the assay deliberately enriched only 374 response genes.
After symbol alignment and removal of all 300 VCC target coordinates, the
source-model input physically contains only the 349 high-sensitivity,
PDS-visible response genes. The original 374-coordinate symbols, VCC mappings,
exportability mask, and selected original indices are retained as provenance,
but excluded coordinates cannot enter response PCA, ridge fitting, held-fold
weights, lambda selection, coefficient fitting, or stability calculations.
The compatibility-format frozen model restores the original 374-coordinate
axis only after fitting and fills all 25 excluded coordinates with exact zeros.
Shallow non-panel transcriptome coordinates are not silently promoted to a
full-transcriptome donor.

One exportable response gene, `BTBD11`, is present on the official and held
score axes but absent from every incumbent donor gene axis, so it is not a
column of any incumbent whitening operator. It remains in held evaluation
instead of being dropped optimistically. The source integration requires the
incumbent and source values to be exact zero on that coordinate and applies
the supported-principal whitening transform only to the remaining response
coordinates.

## Acquisition gate

Acquisition and schema inspection may:

- download only the official files above;
- verify exact sizes and cryptographic hashes;
- list archive members and Matrix Market dimensions;
- inspect guide, hash-label, and feature metadata;
- count guide/cell/condition support without reading expression values;
- write only under `experiments/107-jurkat-high-moi/` and
  `data/external/gse249595/`.

It may not:

- build target effects;
- open VCC panel outcomes or KOLF outcomes;
- run `cell-eval2`;
- modify production code;
- build or submit a challenge artifact.

## Frozen estimator protocol

The original guide-held-out ridge proposal was rejected before expression was
opened: withholding one guide makes its coefficient unestimable, while
retaining its carrier cells leaks its effect into co-guides. The replacement
protocol is frozen in `preregistered_protocol.json`.

Assignments reproduce the source methods as closely as the released artifacts
permit:

- cell keys are `(channel, barcode)`;
- cell calling is a streaming Python port of `DropletUtils::emptyDrops`
  1.31.2 at commit `b423c0fafe7271b173317b407a8538bbbcbede5c`,
  with `lower=100`, `niters=10000`, `test.ambient=FALSE`,
  `ignore=NULL`, `alpha=Inf`, `round=FALSE`, `SerialParam`, and
  `retain=barcodeRanks(..., lower=100)$knee`; the accepted set is
  `FDR <= 0.01`;
- the port preserves the pinned multinomial/Good-Turing/knee/BH semantics but
  has an explicit Python RNG contract rather than claiming bit identity with
  R's `dqrng`: NumPy 2.5.2 `PCG64`, seed
  `SeedSequence([1072026, channel, iteration])`, one serial categorical
  prefix stream per Monte Carlo iteration;
- hash singlets use `deMULTIplex` 1.0.2 at commit
  `9883367c16d374962955e8eeb9fab720a64d7a57`, exactly three
  quantile-sweep rounds over `q=0.01,0.03,...,0.99`, no negative-cell
  rescue, and exactly one known label. Each round removes only `Negative`
  cells; surviving singlets and doublets are reclassified in the next round.
  A round with no new negatives ends the state machine early;
- guide calls use
  `scipy.stats.binom.sf(observed_umi - 1, total_guide_umi, channel_probability)`
  and BH correction within each cell over exactly all 83,401 guides,
  including zero observations as `p=1`, at `q <= 0.05`;
- all 83,401 guide calls contribute to total-guide burden, the pinned 250
  non-targeting records contribute to NTC burden, and only guides mapping
  unambiguously to one target are collapsed to binary target presence.

Before any real transcriptome coordinate is parsed, synthetic tests must show
that the streaming port matches the pinned R implementation exactly for
totals, the knee, ambient probabilities, observed log-probabilities and BH,
and matches independent exact enumeration for the Monte Carlo statistic.
Changing chunk size must not change a call.

The deterministic R-oracle gate passed on 2026-08-30 against the exact pinned
`DropletUtils 1.31.2` and `edgeR 4.99.0` source builds. The receipt is
`R_ORACLE_RECEIPT.json`; `r_oracle.R` is synthetic-only and is not part of the
16-channel production pass.

The primary estimator uses only untreated cells. Channel effects are removed
by deterministic within-channel centering of responses and every design
column within each training or held partition; no held-channel one-hot
coefficient is ever estimated or transferred:

```text
within-channel centered log1p(10,000 * count / full-transcriptome-library-size)
  ~ centered log1p(library size)
  + non-targeting-guide burden
  + total called-guide burden
  + target-presence indicators
```

It is a target-level, rank-64 multiresponse ridge. One penalty from
`{1, 3, 10, 30, 100}` is selected by four channel-blocked folds, each holding
out four complete channels. Nuisance slopes are fitted on training channels
and residualize both the response and the original binary target design;
residualizing only the response is forbidden. Target scales and the response
basis are training-only.

Before any response basis or coefficient is fitted, `model_input.py` computes
and serializes a deterministic design-only audit from the target design,
nuisance matrix, guide evidence, and channel labels. It contains every fold's
support, nuisance aliases, nullspace, projection diagonal, effective rank, and
condition number. `fit_source.py` binds that artifact to the exact input arrays,
validates its nullspace and diagnostics, and reuses it rather than repeating
the expensive sparse nullspace decomposition.

After fitting, coefficients are mapped back to original binary-target units
and projected onto the identifiable row space of the nuisance-residualized
training design. Because total called-guide burden is retained as a nuisance
and cells carry many simultaneous guides, these rows are **identifiable
high-MOI conditional target contrasts**: minimum-norm row-space components
conditional on the observed guide-burden and co-target distribution. They are
not claimed to be uniquely identified singleton marginal effects. External
Nadig, H1, and HEK293T gates decide whether these conditional contrasts are
useful for transfer.

Singular directions with
`sigma / sigma_max <= 1e-6` are discarded. A target is transferable only when
its support gates pass and the corresponding diagonal of the identifiable
projection is at least `0.99`; otherwise its source effect is unsupported and
integration falls back bit-exactly to the incumbent. This projection is
invariant to arbitrary coefficient shifts in the design nullspace, including
the high-MOI guide-burden gauge. Post-fit population-mean subtraction is
forbidden.

Rank-deficiency checks are training-only. Every unsupported, zero-variance or
aliased target, the projection diagonal, effective rank, condition number and
held-cell unmodelled-target burden are recorded in the design audit before
response fitting and propagated into the frozen-model report.
Activated cells are a diagnostic replicate only; they cannot be selected,
averaged, or used to rescue the primary result in experiment 107.

Eligibility is frozen before support is read: at least three distinct called
guides, at least two guides represented in three untreated cells each, at
least 20 untreated carrier cells total, and support across at least four
channels. The thresholds are recomputed using training channels only in each
fold and using all 16 channels only after lambda selection. Expression-QC
retention is reported by guide, target, channel and condition. No
target-specific hand editing is allowed.

## Validation boundary

The first integration is a partial-coordinate residual channel over the
high-sensitivity response genes only. Unsupported coordinates and targets are
exact incumbent fallbacks.

Nadig Jurkat is a **reused development benchmark, not an independent
confirmation set**: committed repository code already used Nadig as a teacher
donor before Experiment 107. Its results can reject this source channel and can
measure incremental value, but Nadig alone cannot authorize promotion.
Promotion additionally requires the frozen H1 and HEK293T PDS and directional-
credit gates. Metadata-only intersection yields 2,194 shared non-VCC targets
and 236 response coordinates after target masking. Targets are SHA-256 sorted,
the first 1,800 are split into six disjoint 300-target galleries, and one
gallery is held out at a time.
For every gallery, its 300 target coordinates and the VCC 300 target
coordinates are removed from source effects, incumbent predictions, truth,
and the operator before any fit or score. The GSE249595 source model is fit
once without Nadig expression or gallery scores. Each strict gallery
reconstructs the incumbent whitener from allowed incumbent-donor rows only;
the a14 gamma, regularization, donor weights, gene-normalization choice and
norm-preservation rule are frozen and are never tuned on Nadig.

The source is not inserted as an ordinary fourth donor: the production donor
model would normalize a 349-coordinate source as though it were full-axis
evidence, and dense whitening would violate exact fallback. The only allowed
integration is a post-whitening masked residual on supported coordinates,
with `f in {0, 0.05, 0.10, 0.20}` and bit-exact incumbent values elsewhere.
The source vector is zero-padded and transformed only through the
supported-by-supported principal submatrix of that gallery's frozen incumbent
whitener, then masked again; GSE249595 rows never refit the operator.
The first screen applies the same candidate to all held contexts. Routing only
to context A is forbidden until a separate cross-context routing protocol is
preregistered and validated; basal identity alone has failed as a selector in
this project before.

Promotion requires:

- pure source PDS at least `0.50` on every Nadig gallery and `0.55` macro;
- nonnegative nested blend PDS delta on every Nadig gallery;
- macro PDS delta at least `+0.003`;
- nonnegative PDS and direction-credit deltas on both promotion-authorizing
  surfaces, H1 and HEK293T; HCT116 remains diagnostic-only;
- target-block bootstrap 95% lower bound above zero for Nadig macro PDS;
- observed gains above 99 whole-target permutation nulls;
- no single target cluster contributing more than half the gain;
- exact fallback outside supported target/response coordinates, verified from
  incumbent/candidate arrays or from a hash-pinned NPZ evidence artifact rather
  than accepted as a caller Boolean.

Before any Experiment-107 Nadig outcome array is opened, the strict validator
now requires `work/strict-validation-inputs.json`. It must have schema version
1, the frozen pre-outcome status, and exactly the artifact roles listed in
`strict_validation_contract.json`. Each role maps to exactly two unique,
repo-relative files: a role-owned NPZ followed by its JSON manifest. The role
manifest binds the NPZ, every array's dtype/shape/C-order/logical SHA-256, every
transitive input file, and its parent roles. Shared donor, PerturbMatch, CD4,
challenge-axis, source-model, and implementation files therefore remain fully
bound without violating the top-level no-duplicate-path rule. The pretruth plan
binds:

- the strict-validation implementation itself;
- the frozen source manifest and the exact current production NPZ schema;
- every file in the exact input manifest;
- the ordered Nadig target and gene axes; and
- the opaque SHA-256 and byte count of `delta.npy`, without parsing it.

Outcome authorization reconstructs every binding and rejects stale or
reordered artifacts. It then creates
`work/nadig-outcome-access-receipt.json` with exclusive-create semantics
**before** NumPy opens `delta.npy`. An existing receipt is a hard replay
failure; the command never overwrites it.

Build preparation is split deliberately:

```bash
# Metadata only. Reconstructs all six Nadig galleries plus H1 and HEK293T,
# checks frozen hashes/counts, and never opens or hashes Nadig delta.npy.
uv run python experiments/107-jurkat-high-moi/strict_artifact_builder.py check

# Expensive: eight sequential eigendecompositions, then exact a14 incumbents
# and all four frozen integration fractions. The explicit confirmation is
# required and outputs stay under this experiment.
uv run python experiments/107-jurkat-high-moi/strict_artifact_builder.py build \
  --confirm-heavy-eigendecomposition
```

H1 uses all three production donors. HEK293T is strict LODO and uses only
K562-GWPS plus Orion HCT116. HCT116 is not a strict promotion input, so no
unbound diagnostic artifact or metric is accepted.

Acquisition success is not a modeling success, and a successful deconvolution
is not submission authorization.

## Acquisition tooling

All acquisition artifacts stay under this experiment:

```text
experiments/107-jurkat-high-moi/
  downloads/   # ignored; partial and completed source files plus receipt
  audit/       # ignored; JSON reports and optional support tables
```

The manifest in `ACQUISITION_MANIFEST.json` pins the three small metadata
files by exact byte count and source SHA-256. NCBI does not publish a
cryptographic checksum for the raw tar in the frozen source inventory, so the
raw archive is pinned by official URL and exact byte count. The downloader
computes a local SHA-256 after transfer and records it in
`downloads/download_receipt.json`; every receipt entry is bound to the exact
manifest hash, asset URL, filename, size, and expected source checksum. Reports
must describe the raw hash as transport-observed, not as an independently
source-pinned checksum.

First validate the frozen manifest/provenance contract without downloading:

```bash
uv run python experiments/107-jurkat-high-moi/audit_archive.py provenance
```

Preview the small metadata acquisition:

```bash
experiments/107-jurkat-high-moi/download.sh metadata --dry-run
```

Download and verify only the small metadata:

```bash
experiments/107-jurkat-high-moi/download.sh metadata
uv run python experiments/107-jurkat-high-moi/audit_archive.py \
  provenance --require-metadata
```

The 3.75 GB archive cannot be selected implicitly. Its command requires an
explicit large-download acknowledgement:

```bash
experiments/107-jurkat-high-moi/download.sh raw \
  --confirm-large-download
```

Interrupted transfers remain as `downloads/GSE249595_RAW.tar.part`; rerunning
the same command resumes with HTTP Range only when its sidecar contains the
same manifest binding and the remote ETag or Last-Modified validator is
unchanged. If the server exposes no stable validator, an interrupted transfer
must be removed and restarted. Curl ignores user configuration, permits only
HTTPS redirects, enforces the expected maximum size, and uses an acquisition
lock. A wrong-size file or source-pinned metadata hash mismatch is never
renamed into place. The unpinned raw archive is never reused without its
manifest-bound receipt, and ambiguous final-plus-partial state fails closed.

After acquisition, run the gates in this order:

```bash
# Receipt bindings, metadata hashes, archive byte count, and a freshly
# recomputed local raw SHA-256.
uv run python experiments/107-jurkat-high-moi/audit_archive.py \
  provenance --require-metadata --require-raw --verify-raw-hash

# Exact 144-member inventory, safe flat paths, official member sizes, and the
# complete 16 x 3 modalities x 3 roles grid. Does not extract the tar.
uv run python experiments/107-jurkat-high-moi/audit_archive.py inventory

# Matrix Market headers and feature axes only. Transcriptome coordinates are
# never scanned.
uv run python experiments/107-jurkat-high-moi/audit_archive.py schema

# Confirm that transcriptome/guide/hash barcode files are byte-identical
# within each channel and that the canonical axis has exactly 6,794,880
# nonempty one-field rows. This reads barcode members but no matrix values.
uv run python experiments/107-jurkat-high-moi/audit_archive.py schema \
  --verify-barcode-hashes

# Stream only guide and hash-label coordinate matrices. This requires positive
# integer values, exact declared entry counts, duplicate-free canonical
# row-major or column-major order, bounded sparse coordinates, and writes
# compressed support tables. It never reads transcriptome coordinate values.
uv run python experiments/107-jurkat-high-moi/audit_archive.py schema \
  --verify-barcode-hashes --scan-assignments

# Expensive integrity gate: read every gzip member to EOF so CRC checks run.
# This includes the transcriptome payload and is intentionally explicit.
uv run python experiments/107-jurkat-high-moi/audit_archive.py verify-gzip \
  --confirm-full-stream
```

The schema audit is sparse by construction:

- it reads Matrix Market dimensions before allocating support counters;
- it refuses a declared or observed guide-plus-label total above 200 million
  entries and aborts on the first entry beyond a matrix declaration;
- it uses fixed-width arrays proportional to feature and barcode axes;
- it processes one assignment matrix at a time;
- it never calls `scipy.io.mmread`, creates a dense matrix, or extracts the
  archive;
- transcriptome coordinate scanning is hard-disabled in the manifest;
- support tables are built in a staging directory and published only after all
  matrices pass, together with `support/COMPLETE.json`;
- every inventory/schema/gzip report is bound to a freshly recomputed raw
  archive SHA-256, the manifest SHA-256, and the provenance SHA-256.

## Acquisition/schema promotion gates

Model implementation remains blocked until all of these pass:

1. the frozen provenance fields and official URLs match;
2. all small metadata files match their source-pinned byte counts and SHA-256;
3. the raw tar matches the official 3,753,820,160-byte size and its acquisition
   receipt;
4. the archive has exactly 144 safe regular-file members with no missing,
   duplicate, nested, or extra paths;
5. every gzip member reaches EOF without CRC failure;
6. each channel has transcriptome, guide, and label matrices with exactly
   6,794,880 barcode columns and respectively 20,606, 83,401, and 7 rows;
7. feature axes are identical across all 16 channels;
8. all three barcode files are byte-identical within every channel;
9. every barcode axis has exactly 6,794,880 nonempty identifiers;
10. every guide and label coordinate is in bounds, positive-integer, unique,
    canonically ordered, and the observed entry count equals the Matrix Market
    declaration;
11. guide/cell/hash support reports and their input-bound completion manifest
    exist before any QC threshold is chosen.

Passing these gates authorizes only the frozen source model and held-context
screen described in `preregistered_protocol.json`. It does not authorize a
challenge build, submission, context-A-only route, count-level score, or any
post-result threshold/grid change.

## Lightweight tests

Synthetic tarballs exercise authenticated resume state, final/partial
ambiguity, stale receipt rejection, raw tamper detection, path confinement,
the exact member grid, barcode row counts, Matrix Market dimensions and entry
bounds, negative/duplicate coordinate rejection, staged-output cleanup,
support reports, path-traversal rejection, provenance mismatch rejection, and
gzip EOF/CRC validation:

```bash
uv run python -m unittest \
  experiments/107-jurkat-high-moi/test_acquisition.py

bash -n experiments/107-jurkat-high-moi/download.sh
uv run ruff check \
  experiments/107-jurkat-high-moi/acquire.py \
  experiments/107-jurkat-high-moi/audit_archive.py \
  experiments/107-jurkat-high-moi/test_acquisition.py \
  experiments/107-jurkat-high-moi/test_protocol.py

uv run python experiments/107-jurkat-high-moi/test_protocol.py

uv run python -m unittest \
  experiments/107-jurkat-high-moi/test_strict_nadig_validation.py

uv run python -m unittest discover \
  -s experiments/107-jurkat-high-moi -p 'test_*.py'

uv run ruff check \
  experiments/107-jurkat-high-moi

uv run python -m py_compile \
  experiments/107-jurkat-high-moi/strict_nadig_validation.py \
  experiments/107-jurkat-high-moi/test_strict_nadig_validation.py
```

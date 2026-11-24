# 111 — DepMap context-program gain

Status: **source-only complete; NO-GO; held truth remains unopened**

Date: **2026-08-30**

## Question

Can a target-independent response-program operator learned from the same
perturbations across 16 cancer contexts use a supplied control profile to
improve an a14-shaped prediction without giving back its retrieval geometry?

The released Broad/DepMap screen contains 100 Cas9 knockout targets across 16
cell lines and 11 lineages. The target set has zero overlap with the 300-gene
VCC 2026 panel. It cannot be a direct prediction donor. Its distinctive asset
is the repeated target-by-context response tensor.

The proposed model family and all promotion thresholds are frozen in
`preregistered_protocol.json`. Acquisition and source-only nested validation
are implemented here. There is no held-truth or count-level evaluator.

## Prior-family boundary

- Experiment 042 rejects adding mismatched auxiliary perturbations as ordinary
  whitening rows.
- Experiment 046's apparent context adapter lost to target-label permutations
  and wrong-context maps.
- Experiment 048 identifies the surviving effect as generic source geometry,
  not context adaptation.
- Experiment 072 shows why low target overlap cannot authorize a direct donor.
- Experiment 082 rejects a linear ESM2 target-kernel residual.
- Experiment 080 already owns the context-independent listwise metric-adapter
  family and remains unrun.

Experiment 111 is narrower and different: it permits only multiplicative gains
on target-specific response-program coefficients, predicted from a context's
control expression. Additive context vectors, target embeddings, target
identity, gene identity, and ordinary covariance pooling are forbidden.

## Frozen sources

`SOURCE_MANIFEST.json` pins:

- preprint DOI `10.64898/2026.08.24.746802`;
- Figshare data article `33273600`, including every top-level file id, name,
  size, and published MD5;
- Figshare code article `33282171` and its 234,195-byte archived source;
- GitHub commit `ccb9ee5c37fb5a117cd76d8e9a5a5217ce12338d`;
- the GitHub commit archive's byte count, MD5, SHA-256, and normalized tree
  SHA-256;
- the byte-identical normalized tree SHA-256 shared by the GitHub and Figshare
  code archives, using an explicit path/NUL/content/NUL canonicalization;
- the complete processed ZIP central-directory offset, size, entry count, and
  SHA-256;
- the exact CRC32, compressed size, uncompressed size, flags, compression
  method, and local-header offset of the five authorized processed members;
- the exact metadata ZIP and `cell_line_metadata.csv` hashes;
- the exact 16 cell lines, 11 lineages, six assay batches, 100 targets, ten
  olfactory controls, seven absent pseudobulk context-target pairs, and zero
  VCC-panel overlap.

The full `processed.zip` SHA-256 is intentionally not claimed: computing it
would require the forbidden 3.23 GB download. Its Figshare-published MD5 is
bound into every selective-extraction receipt, while each extracted member
receives a locally computed SHA-256 and CRC32 bound to the exact manifest and
central directory.

## Authorized members

Only these data members may be published under
`data/external/depmap26q1/extracted/`:

```text
metadata/cell_line_metadata.csv
processed/pseudobulk_sum_by_ko.csv
processed/lfc_matrix.csv
processed/fdr_matrix.csv
processed/crispr_table.csv
processed/cell_quality_covariates.csv
```

The acquisition code cannot select `single_cell_data.zip`,
`h2h_data_formatted.zip`, `downloaded.zip`, `data.zip`, the full
`processed.zip`, or any unlisted processed member.

## Acquisition design

The large processed archive is never downloaded as a whole.

1. Fetch the last 131,072 bytes with an HTTP Range request.
2. Require status `206` and an exact `Content-Range` total of
   `3,233,927,690`.
3. Parse the EOCD and require the pinned central-directory offset, size, and
   entry count.
4. Range-fetch the central directory and require its pinned SHA-256.
5. Reject unsafe, duplicate, encrypted, symlink, malformed, or unexpected
   selected entries.
6. For each authorized member, validate its local header and stream only its
   compressed byte interval in chunks no larger than 8 MiB.
7. Inflate into a temporary file while computing CRC32, SHA-256, and exact
   uncompressed byte count.
8. Publish atomically only after every check succeeds.
9. Bind the receipt to the source-manifest SHA-256, Figshare article/file
   metadata, processed central-directory SHA-256, selected-member metadata,
   and extracted-file SHA-256.

`metadata.zip` is only 1.84 MB, so it is downloaded completely and checked
against both its pinned MD5 and SHA-256. Every member path is audited, but only
`metadata/cell_line_metadata.csv` is extracted.

HTTPS-only redirects, exact response sizes, output confinement, a nonblocking
acquisition lock, atomic publication, and stale-receipt rejection are
mandatory. A server that ignores Range and returns `200` fails immediately
without reading the response body.

## Commands

Local manifest, panel-overlap, and protocol checks:

```bash
uv run python experiments/111-depmap-context-program-gain/acquire.py manifest
```

Read-only live provenance and ZIP-directory audit. This transfers only API
JSON, two sub-250 KB code archives, the 1.84 MB metadata archive, and ZIP
directory/range metadata—not selected expression members:

```bash
uv run python experiments/111-depmap-context-program-gain/acquire.py remote
```

Preview the selective acquisition:

```bash
experiments/111-depmap-context-program-gain/download.sh acquire --dry-run
```

The approximately 300 MB compressed selective transfer requires the exact
separate execution amendment and an explicit acknowledgement. The original
preregistration remains unchanged and closed:

```bash
experiments/111-depmap-context-program-gain/download.sh acquire \
  --confirm-selective-download
```

Verify the authenticated local extraction and receipt:

```bash
uv run python experiments/111-depmap-context-program-gain/acquire.py local
```

`local` verifies the receipt and streams every CSV through its frozen schema:

- exact 16-line/11-lineage/six-batch metadata;
- exact 1,600-row cell-line-by-target CRISPR table;
- exact 100 targets and ten olfactory controls;
- exact 1,593 pseudobulk columns and seven missing pairs;
- exact 1,440 SCEPTRE LFC/FDR columns;
- unique, nonempty response-gene axes;
- finite nonmissing values, nonnegative integral counts, and FDR values in
  `[0,1]`;
- source-native empty-field missing values, with no entirely missing response
  rows;
- bit-identical LFC/FDR response axes and missing-value masks;
- exact zero overlap with the pinned 300-target challenge panel.

The exact 70-row olfactory-control metadata exception in `crispr_table.csv` is
bound by `ACQUISITION_SCHEMA_RECOVERY.v1.json`; no analysis-target exception is
permitted.

Source-only preflight and execution:

```bash
uv run python experiments/111-depmap-context-program-gain/source_model.py \
  preflight

uv run python experiments/111-depmap-context-program-gain/source_model.py \
  run --confirm-source-only-fit I_CONFIRM_EXPERIMENT_111_SOURCE_ONLY_V1
```

## Lightweight validation

The synthetic suite uses small in-memory ZIP files and fake range sources. It
does not contact Figshare, open donor matrices, read held truth, or fit a
model.

```bash
uv run python -m unittest \
  experiments/111-depmap-context-program-gain/test_acquisition.py

uv run ruff check \
  experiments/111-depmap-context-program-gain/acquire.py \
  experiments/111-depmap-context-program-gain/test_acquisition.py

uv run ruff format --check \
  experiments/111-depmap-context-program-gain/acquire.py \
  experiments/111-depmap-context-program-gain/test_acquisition.py

bash -n experiments/111-depmap-context-program-gain/download.sh
python3 -m json.tool \
  experiments/111-depmap-context-program-gain/SOURCE_MANIFEST.json >/dev/null
python3 -m json.tool \
  experiments/111-depmap-context-program-gain/preregistered_protocol.json \
  >/dev/null
```

The separately frozen execution amendment authorized the exact selective
acquisition and source-only fit. The authenticated files were published under
`data/external/depmap26q1/extracted`; the immutable source-only outcome is under
`work/source-only-v1`.

The nested source result is a **NO-GO**:

- macro raw PDS delta `+0.00490196`, below the frozen `+0.010` gate;
- context-label complete-selection null p95 `+0.00680828`;
- target-bootstrap lower 95% `-0.00095316`;
- the gain is concentrated in one batch and the shallow-effect half.

See `RESULTS.md`. No candidate or held-evaluation proposal was created, and no
held truth, `cell-eval2`, counts, production build, submission, commit, or push
was used.

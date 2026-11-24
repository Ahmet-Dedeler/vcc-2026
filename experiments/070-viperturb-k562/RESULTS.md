# 070 — VIPerturb K562 results

Date: 2026-08-23

## Status

**Corrected evaluation completed; reject VIPerturb for production.**

The evaluation-only retry `20260823T070153Z` reused the authenticated compact
from run `20260823T063300Z` and exited successfully without rereading the RDS.
The frozen candidate reduced full-gallery PDS in both valid held-donor folds:

| held donor | arm-D baseline | VIPerturb candidate | delta |
| --- | ---: | ---: | ---: |
| Orion HCT116 | 0.8023716840 | 0.7968449635 | -0.0055267205 |
| Orion HEK293T | 0.7802287582 | 0.7718305459 | -0.0083982122 |
| **macro** | **0.7913002211** | **0.7843377547** | **-0.0069624664** |

The paired target bootstrap 95% interval was
`[-0.0108373282, -0.0034542003]` with zero positive bootstrap probability.
Resetting every modified low-cell target still left a negative
`-0.0046015955` full-gallery delta. The candidate failed both-fold,
macro-gain, confidence-bound, null, and low-cell-sensitivity gates. Stop this
branch; do not run a count-level confirmation or use VIPerturb in a challenge
submission.

The retry peaked at `4.5166 GiB` process-tree RSS, all 16 synthetic tests
passed, and the downloaded artifacts verified against
`cloud-output/20260823T063300Z/evaluation-retry/20260823T070153Z/artifact-sha256.txt`.
No challenge submission was made.

Run `20260823T063300Z` completed the full extraction and wrote an authenticated
16 MB compact. Extraction peaked at `48.484 GiB` process-tree RSS
(`47.246 GiB` at the first observed peak) and the failed evaluator peaked at
`2.619 GiB`. Evaluation stopped before the biological comparison because the
stored panel scale was inconsistent with the raw panel rows.

The cause had two parts. R's `panel_effects - panel_center` recycled the
18,111-element center through the column-major 124 x 18,111 matrix rather than
subtracting it by column. That produced `12.50221075846`; jsonlite's default
four-decimal serialization stored `12.5022`. The canonical column-wise value
from the authenticated effects and center is `12.26301798025`.

The extractor now uses `sweep(..., MARGIN = 2)` and writes JSON numbers at full
precision. The Python loader accepts the preserved compact only because its
stored value matches the exact legacy recycling formula within the known
four-decimal serialization bound, then uses the recomputed canonical value.
Values matching neither the canonical nor legacy formula still fail closed.
An evaluation-only retry launcher reuses the existing compact and donor inputs
without invoking R or accessing the RDS.

The corrected cloud run `20260823T054706Z` failed before expression extraction
because the official manifest repeats its `NO-TARGET` control row once for
each of `partA`, `partB`, and `partC`. The RDS itself is a Seurat v5 object
whose RNA assay is an `Assay5` with one `counts` layer; `object[[]]` has the
columns `orig.ident`, `nCount_RNA`, `nFeature_RNA`, `nCount_GDO`, `guide`,
`gene`, `sample`, and `nFeature_GDO`. Its `gene` column contains 6,724
perturbations plus one `NO-TARGET` control label. The extractor had treated
all 6,727 filtered manifest rows as unique perturbations and therefore stopped
before calling `readRDS()`.

The minimal correction now removes recognized control rows before validating
the 6,724-entry perturbation catalogue. It still fails closed on missing or
ambiguous controls, empty labels, any duplicated non-control perturbation, an
unexpected perturbation count, or any non-control metadata label absent from
the frozen manifest catalogue. A base-R regression fixture reproduces the
three repeated manifest controls and the observed eight-column Seurat
metadata schema.

The previously started run `20260823T052157Z` on worker
`i-08d90ae88c7b210ed`, with durable output prefix
`s3://<your-s3-bucket>/output/exp070-viperturb-k562/20260823T052157Z/`,
was packaged before the source-wide correction. Any output from it is
superseded and must not be used for promotion or submission. No cloud worker
was launched, restarted, or stopped during this diagnosis; the already-running
failed worker was queried only to inspect the retained RDS schema.

The corrected protocol now normalizes across all perturbations in the filtered
VIPerturb source object before retaining the 124 panel rows. It adds a
panel-centered diagnostic ablation, target-label and gene-coordinate nulls,
and a fail-closed low-cell sensitivity gate. The frozen arm-D fallback remains
exact for unsupported targets. Real biological results must not be inferred
from synthetic tests.

An earlier launch, `20260823T045858Z`, failed before tests, data download, or
scientific computation because Amazon Linux's installed `curl-minimal`
conflicted with an unnecessary request to install the full `curl` package.
The host bootstrap now retains the working system curl and installs only the
missing utilities. The rerun uses a fresh per-run work directory.

## Local verification

- Synthetic harness: `16 synthetic tests passed`.
- Python self-reported peak RSS: `0.046 GiB`.
- Guarded process-tree peak RSS: `0.0056 GiB` in the latest local run.
- Memory limit: `1.0 GiB`; return code `0`, no kill.
- The exact mirrored 16 MB compact passes the patched loader at `0.102 GiB`
  peak process-tree RSS. It is classified as
  `legacy_r_recycling_jsonlite_four_decimal`, preserves stored `12.5022` for
  audit, and uses canonical `12.263017980250217` for the diagnostic ablation.
- Live-manifest regression: 18,888 rows, 6,724 non-control source
  perturbations, one control label (`NO-TARGET`) represented by three manifest
  rows, 124 VCC-panel targets, and 22 low-cell panel targets.
- Ruff lint: passed.
- Ruff formatting check: passed.
- Python compilation: passed for all four Python files.
- AWS shell syntax: `bash -n` passed for both scripts.
- R source-wide/static contract and delimiter checks: passed. A pure-base-R
  schema regression is wired into the cloud runner before any RDS download;
  local execution remains unavailable because this Mac has no R runtime.

The live official manifest was streamed without saving the RDS and passed:

- 18,888 manifest rows;
- 124 filtered VCC targets;
- 22 targets below 30 cells;
- minimum manifest cell count 8;
- 289-target held gallery;
- 98 expected incumbent-K562/PerturbMatch/VIPerturb shared targets, above the
  70-target gate.

## Unresolved environment risk

`Rscript` and micromamba are not installed on the local Mac, so the corrected R
extractor and complete conda environment could not be executed locally. The
runner uses pinned conda-forge Linux binaries for R 4.4.3, Seurat 5.3.1, and
SeuratObject 5.1.0, whose individual Linux builds were verified available.
The failed worker's retained RDS verified the exact object, assay, layer, and
metadata schema, but the post-fix extractor has not been run against that RDS.
The extractor handles single or multiple counts/data layers, audits its
selected metadata column, processes source groups in bounded chunks, and fails
closed on missing source perturbations/cells, ambiguous mappings, insufficient
gene overlap, wrong panel coverage, or invalid normalization scales.

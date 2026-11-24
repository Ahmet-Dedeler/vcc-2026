# 324 — Are independent downstream-residual labels identifiable from existing metadata?

September 9, 2026. Solo Codex, following 323. **Metadata-only feasibility, not
model training, response measurement, an independent confirmation, or a
literature novelty claim.** Commit this protocol and the implementation before
running the real inventory. Source headers and earlier results were inspected
to design this test; the crossed counts below have not been computed.

## Decision and hypothesis

Can an already inventoried source support the proposed new ingredient:
repeatable, context-specific downstream-response labels across independent
experiments, without renaming a scalar reliability feature or basal ridge?

The hypothesis is that the nominal many-study/many-line coverage conceals too
little **crossed** replication to identify that ingredient. A negative closes
this proposed label construction on the inspected metadata, not all response
models or all publicly available data.

## Authorized inputs

Only these existing files may be opened by the inventory:

- `experiments/270-data-reconnaissance/work/gpsadb_metadata_browse.csv`
  (target, line, accession, method, sample-count and annotation metadata).
- `data/challenge/pert_counts.csv` (the target list, not expression).
- `experiments/285-reliability-benchmark/work/inv285.json` (existing metadata
  inventory; use only the explicitly named metadata fields).
- `experiments/268-broad16-context-residuals/work/units268.md` (existing
  12-working-line batch/count table).

Hash-bind them, this protocol, and the implementation in `inputs.json`.
No raw HDF5, response/DE matrices, pseudobulk arrays, reserve contents, source
test targets, new data acquisition, or model calls. The old R1–R5 files and
132's sealed 410-target outcome stay unopened. An accession is not proof of
an independent publication, sample pool, control pool, or guide construct.

## Frozen inventory definitions

1. Reproduce 270's raw metadata totals: 7,665 dataset groups, 4,623 accessions,
   96 panel targets, 329 panel groups and 50 panel targets with at least two
   accessions. Abort if these known-answer counts fail; do not silently edit
   the expected numbers.
2. Report the historical raw line-name convention and a second convention:
   use Cellosaurus ID when present, otherwise an exact stripped line name.
   Do not silently merge unregistered aliases or infer tissue identity.
3. Count targets with two accessions and two contexts; separately count
   **same-context** multi-accession target/context pairs.
4. A crossed rectangle is one target, two context keys, and the same two
   accessions observed in both contexts. Deduplicate dataset groups into
   edges first; multiple comparisons within one accession are not studies.
   Report all-method and **same-method** rectangles separately. Preserve
   shRNA, siRNA, knockout and CRISPRi as distinct methods.
5. Report the same-method, stable-Cellosaurus-ID-only subset as a conservative
   identity sensitivity. Report exact groups supporting each panel rectangle
   and each panel same-context pair for later source-provenance checks.
6. Compare with the existing single-cell inventory at the actual independent
   unit: Replogle's two K562 experiments are one publication and no shared
   panel targets; Nadig's two lines, Broad's 12 working lines, and CD4's four
   donors/three states do not themselves create independent publications.
   Do not open response arrays to strengthen this comparison.

Counts above are descriptive, not effect estimates or statistical power.
Two-by-two crossing is necessary for a direct matched context contrast across
the same two studies, not necessary for every possible prediction design.
No arbitrary count threshold will be sold as sufficient power.

## What could earn a response-label pilot

Metadata must first identify genuinely independent experimental/sample/control
and construct units for the same target/context/method. Publication provenance
must be resolved. A crossed contrast needs those units in both contexts, not
just a target observed once in each of two different lines.

Even qualifying metadata only earns a separately registered **development**
repeatability measurement. It cannot establish repeatability, six-metric
benefit, construct equivalence, generalization, or an untouched confirmation.
If provenance is absent, classify it as unresolved; do not impute independence.
No production candidate or response fit is authorized by this inventory.

## Design constraint for a later residual label

With independent noisy responses `y_A = s + r + e_A` and `y_B = s + r + e_B`,
subtracting the **same noisy shared estimate** `s_hat = s + u` makes
`E[(y_A-s_hat) dot (y_B-s_hat)] = ||r||² + E[||u||²]` under independent
zero-mean errors. Shared control noise produces an analogous artifact.
This is a known measurement-error issue, not a claimed new biological finding.

A later cross-study residual must use independent shared-estimate/control
views (or an explicitly validated correction), exclude the recipient from
every training example/operator, and include wrong-target and wrong-recipient
controls. Add a small deterministic synthetic test demonstrating both the
false-positive shared-noise case and its independent-view repair. Do not use
this synthetic demonstration as empirical evidence about 268.

## Novelty boundaries

- 132 replaced direct same-target means with DepMap PCA/ridge: failed.
- 161 already tested a context residual; its target-only consensus reshaping
  is a mandatory baseline, not a new proposal.
- 224's small per-gene diagonal **succeeded at effect level**; higher-capacity
  and recipient-conditioned models did not improve its result. Do not call
  the whole experiment a failure.
- 268 already measured half/guide reproducibility and basal-conditioned
  residuals on 12 working lines. Its line/batch nesting remained a limitation.
- 269 exposed held-line leakage through other folds' donors and failed its
  utility routers. A new output model cannot reuse those leaky label objects.
- 283 spent its reserved lines; its preregistered failed/positive claims stand.
- 285 measured guide/donor/batch effects and two K562 experiments; it did not
  turn one publication into two independent studies.
- 270 already measured GPSAdb top-gene sign stability as a scalar feature.
  Here only crossed metadata and identifiability are new, not that statistic.

## Deliverables and operating bound

`work/results.json`, supporting CSVs, `RESULTS.md`, and unit tests. Stdlib-only
inventory; one local process, under 1 GiB RAM and 25 MiB new disk. No scorer
metrics are measured: all six are explicitly **unmeasured**, because no
prediction changes or response values are read. Local commit only, no push.

# 303 — Pre-access recovery and validation of the second-library test

The original `MANDATE.md` is an unexecuted draft. This mandate supersedes its
execution authority. The user asked Codex to take over on 2026-09-07.

## This phase is PRE-ACCESS ONLY

Do NOT read the Nakatake rows of `X` or `layers/control`, directly or through
an eager AnnData load. Do NOT read the spent Joung atlas X or any other reserve.
Do not compute a C7 outcome. You may inspect obs/var/uns and measurement masks,
source code / methods, the already-used Pan2025 rows, existing banks, and saved
302 artifacts. Keep the reserve closed even if all your checks pass.

Own only `experiments/303-second-overexpression-library/`. The coordinator
alone edits the shared ledger and freezes/releases C7. No submissions, builders,
production changes, commits, publication, or new agents. 304 owns the independent
novelty/F15 audit; 305 owns the T7 expression test. Do not duplicate their work.
Use `apply_patch` for authored files, `uv run python`, <6 GiB RSS, no large
download. Time box 60 minutes for this phase.

## Immediate deliverable

Write `PREACCESS.md`, `PREREGISTRATION-DRAFT.md`, metadata/input manifests,
a guarded single-access runner and synthetic/preflight tests. Your final
`RESULTS.md` must say **PREACCESS ONLY — R5 UNREAD — NO C7 RESULT**. Report
`READY_FOR_REVIEW` or the precise precondition failure, not scientific PASS/FAIL.
The coordinator will inspect and authorize any later phase explicitly.

Read READ-FIRST, AGENTS, the research programme, the R5 ledger and original
303 draft, 302's RESULTS / PREREGISTRATION / lib302 / ka302, and experiment
115's **extract_compact_donor.py** (not just acquire.py).

## Checks and corrections required BEFORE freezing C7

1. Reconstruct exactly what exists and what was never launched. Confirm metadata
   cohorts 100 / 58 / 13, reporting any discrepancy rather than forcing counts.
   Record an explicit ordered row-index list for R5 and gene axis with hashes.
   Do not confuse target genes (rows) with measured genes (columns).
2. Establish the matrix meaning from the ORIGINAL producer's source code and
   methods: X vs control, log base, normalization, sign handling, missing probes,
   mask, duplicated symbols, platform, induction time and original cell line.
   hESC is not an identity. The draft's claims of RNA-seq, same line, and no
   reliability ceiling are hypotheses, not authority. You may check primary
   methods, but never source-notebook outputs with reserve expression results.
3. On already-used Pan2025 rows, reproduce the actual bank loader and its
   delta.npy, b_ctrl.npy and measured_mask for a fixed seeded 20-target sample
   (tolerance justified by stored dtype). Reproduce 302's KA-C / 279's Ccross
   through the existing path. Distinguish loader regression from biological
   validation. Abort before access if either fails. Ensure imports have no
   writes to another experiment.
4. Draft a scientifically valid C7 retaining the four planned questions:
   antisymmetry; conservation vs mirroring; geometry vs signed alignment;
   geometry vs absolute alignment. Use observed cross-study responses as the
   estimand unless actual identical cell lines/protocols are established.
   Do not promise latent biological equivalence from a null result.
5. Center both study objects on the SAME predeclared overlapping TF cohort
   for the primary residual statistic; do not silently compare a 510-row
   centered OE vector to a 2354- or 175-row-centered KD vector. Carry the old
   draft's bank-wide convention only as a labelled sensitivity. Use genuinely
   measured coordinates; no padded zeros as measured expression. Remove the
   own-target coordinate, including for shuffled pairings of two different TFs.
   Include the direct control-referenced (unresidualized) effect as a secondary
   so centering-induced conclusions can be recognized.
6. Require the matching effect to beat a target-shuffled pairing null, not only
   zero. Preserve the draft's median/tail statistics but do not assume a centered
   null has median zero. Four tests require explicit multiplicity control
   (Holm familywise 0.05); report unadjusted intervals as descriptive too.
   Add family-aware robustness and matched-null checks where metadata supports
   them; say when sample size does not support an adjustment.
7. C_held and the KD geometry share data with the alignment. For conservation,
   prefer a frozen source-excluded conservation feature recomputed from saved
   map rows with gse295214_h1 excluded, recording all sources. If unavailable,
   label the original C_held statistic as coupled/descriptive, not independent.
   For geometry-versus-|S|, preserve KD features in permutations and randomize
   OE-to-target matching to test whether geometry explains generic chance
   alignment. Partial magnitudes alone do not eliminate that coupling.
8. Cross-library comparison to 302 is DESCRIPTIVE and conditional on 302's
   unresolved processing assumptions. A different study is not a repeat
   measurement of an identical intervention; never treat that agreement as a
   plug-in reliability for disattenuation or a construct-noise ceiling.
9. Create synthetic zero, exact-mirror and exact-alignment tests; wrong-axis,
   missing-mask and failed-known-answer refusal tests; deterministic seeds,
   finite-value/type validation and a test ensuring metadata/preflight mode
   never reads an R5 X/control slice.
10. Build fail-closed execution: the single reserve read must require a
    coordinator-created release JSON binding the source hash, row-index hash,
    final preregistration hash and code hashes. Record a durable access-start
    marker BEFORE reading, and refuse a second access even after interruption.
    Save the necessary sufficient arrays once with hashes; later computation
    uses only the saved arrays. A marker is not permission to resume a partial
    source read. No release JSON exists in this phase; do not create one.

In PREACCESS.md, clearly separate confirmed metadata, unconfirmed provenance,
design amendments to the interrupted draft and unresolved risks. These are
changes BEFORE any outcome has been viewed, not post-hoc repairs. End this
phase for coordinator review rather than waiting in a polling loop.

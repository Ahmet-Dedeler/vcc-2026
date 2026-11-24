# 087 — STRING target-kernel residual

**Verdict: NO-GO before held-panel evaluation.**

The implementation follows the frozen README and experiment 082's corrected
infrastructure:

- reconstructs experiment-064 arm D independently for held HCT116, HEK293T,
  and KOLF roles;
- freezes the exact experiment-082 60/20/20 target lists and hashes before the
  first fit;
- keeps development-target coordinates out of development preprocessing and
  uses the pinned 17,806-gene final score axis;
- loads only the pinned 21,688 × 64 STRING HDF5 artifact;
- represents missing STRING targets as literal zero rows with an explicit
  availability mask and exact arm-D fallback;
- searches ranks 8/16/32/64, the frozen ridge grid, the frozen residual-fraction
  grid, and mandatory identity;
- repeats complete selection for stratified STRING-row permutation and
  equal-dimensional Gaussian nulls;
- serializes and hashes real/null final weights and source-query clusters
  before authenticating held-panel truth;
- evaluates the complete 290-target HCT116 and 294-target HEK293T galleries,
  rather than the legacy 259-target shared intersection, and keeps KOLF's
  complete 273-target gallery outside the primary macro, bootstrap, and null
  thresholds while requiring its independent nonnegative confirmation;
- retains the 1 GiB local synthetic guard and fail-closed AWS launcher/worker.

Local validation on 2026-08-29:

- 23 synthetic tests passed;
- peak process-tree RSS: see `local-test-memory.json`;
- Python compilation passed;
- shell syntax and JSON parsing passed;
- dependency SHA-256 checks passed;
- the pinned STRING artifact was verified as 21,688 root datasets × 64
  float64 features with 298/300 panel coverage; missing targets are exactly
  `EEF1A2` and `KIF21B`;
- scoped Ruff checks passed with the inherited executable/import-layout and
  redundant-cast rules excluded.

No donor-truth score, count generation, `cell-eval2` run, challenge build, or
submission has completed.

Two 2026-08-30 cloud attempts failed before held-truth staging:

- `20260830T004152Z` stopped in the cloud synthetic test because the source
  archive omitted the RPE1 and Nadig-Jurkat `genes.txt`/`perts.txt` files used
  by the real-axis contract test.
- `20260830T005453Z` passed that test after the archive fix, then stopped at
  source staging because the expected SHA-256 for
  `cd4_stim8hr/perts.txt` contained a one-character transcription error
  (`...6cf1...` instead of the independently reproduced
  `...6cf3...`).

Both failures occurred before cache construction, model fitting, or any held
HCT116, HEK293T, or KOLF perturbation value was used. The archive omission now
has a regression assertion, and the corrected CD4 hash was verified against
the local file, the exact S3 tarball, and experiment 085's independent frozen
manifest.

## Authoritative cloud run

Run `20260830T011413Z` completed successfully on the isolated
`r7i.8xlarge` worker after the two fail-closed staging repairs above.

- all 26 synthetic and contract tests passed;
- source/query staging passed every leakage invariant;
- the complete frozen grid contained 80 non-identity candidates plus identity;
- **zero** non-identity candidates satisfied the inner all-cell nonnegative
  selection rule;
- the best non-identity candidate was
  `rank=64, ridge=0.1, beta=0.02`, with inner macro PDS delta
  **`-0.004302`** and minimum cell delta **`-0.007443`**;
- identity was therefore selected;
- the untouched non-panel pretest was exactly identity (`delta=0`);
- both complete-selection null p95 values were `0`;
- no outer HCT116, HEK293T, or KOLF numerical value was indexed into fitting,
  selection, or prediction, and no held-panel score was computed.

Procedural caveat: staging did hash the held donor files and the shared loader
memory-mapped their arrays before the pre-gate stopped. Thus the literal claim
that held-truth files were “not opened” would be false. This is a provenance
boundary defect, not demonstrated numerical leakage: the candidate selected
identity without reading held values into any computation.

The preregistered pre-gate failed on selected-nonidentity, inner
nonnegativity, positive non-panel test gain, and both null comparisons. No
held-panel evaluation, count generation, `cell-eval2` run, challenge build, or
submission followed.

All returned artifacts were downloaded from S3 and verified against
`cloud-output/20260830T011413Z/artifact-sha256.txt`. Evaluation peak
process-tree RSS was `0.390 GiB`; exit code was `0`. The worker was stopped
after confirming no other VCC process or service was active.

Machine-readable result:
`cloud-output/20260830T011413Z/evaluation/results.json`.

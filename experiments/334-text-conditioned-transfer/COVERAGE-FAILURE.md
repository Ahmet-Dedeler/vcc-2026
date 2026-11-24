# 334 original run — text-coverage stop, before fitting

Original freeze **ba649342**, September 9, 2026. The guarded command
exited 1 before the first model fit. Ten pure tests had passed, and all
five historical raw-loader/centered-array checks passed in the run.
No prediction, model, terminal accuracy result or production output was
written. The process group is gone.

Failure: `representations()` found **971/1,024** usable target texts,
below the frozen 95% requirement, and **482/512** usable gene texts.
The target failure is preserved; do not call the original run a pass.
The memory guard's 0.01-GiB sampled peak missed the short execution and
is not an adequate peak-memory measurement.

Metadata-only diagnosis: 50 missing target descriptions are exact
`Gene Symbol KEY` headers with no description body; three primary keys
are absent. The combined NCBI/UniProt dictionary in the same verified
archive fills those 50 bodies under the same exact-key/header check.
It also fills 25 response-gene bodies. Fallback coverage is therefore
1,021/1,024 targets and 507/512 genes. No response accuracy was read.
This source-augmentation possibility motivates a disclosed recovery,
not a retrospective change to the original feature artifact.

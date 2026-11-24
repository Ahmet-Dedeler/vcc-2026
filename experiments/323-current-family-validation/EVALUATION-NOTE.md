# Evaluation implementation freeze — September 9, 2026

Written after the construction gates passed, before new reference statistics
or model scores in 323.

The replacement operator and model composition are built. Source replay matches
170 exactly and matches 319's two A/B/C pre-whitening hashes. This is an
engineering result, not a held-reference outcome.

`evaluate323.py` evaluates actual emitted count blocks without serializing and
re-reading whole prediction files. It uses production's pre-sorted-control
Wilcoxon kernel and the installed cell-eval2 0.16.0 PDS, jackknife/MSE,
DE-schema preparation and four DE-metric functions. These sufficient statistics
are not predicted profiles used in place of emitted counts. Exact integer block
totals and an ordered count digest are retained.

Before any new arm, reconstruct 221's historical rc_callset effect/emission:
the combined panel-total digest must equal its existing build record. The six
raw aggregates must match its recorded full-scorer row within absolute
tolerances PDS 1e-12, MSE 1e-10, and each DE metric 1e-6. Tolerances allow
floating-point reductions but not a rank step or changed call set. Failure
stops the fast evaluator; record the failure and diagnose without changing
candidate settings or selecting targets.

The historical known-answer adapter explicitly uses minnorm and is only a
known-answer row. Current arms explicitly use the production unset mass fix.
All 85 targets and all emitted genes are retained, including empty reference
DE pools according to cell-eval2's own conventions.

No usable H1 real-bundle anchors were found among retained local artifacts;
a read-only listing of the historically named S3 output prefix was also empty.
Do not manufacture an overall score by applying A/B/C anchors to H1. Until
H1-specific anchors are independently recovered or rebuilt, overall-sign
qualification is **unassessed**, not passed. Raw six-member comparisons,
negative-control checks, and paired per-target uncertainty can still be
completed and used to describe this reference. This limitation does not change
the original qualification rule.

The evaluation code, supporting package files, historical known-answer inputs,
and composed effects are hash-pinned in `evaluation-inputs.json` and committed
before execution. No new primary-data download or reserve read.

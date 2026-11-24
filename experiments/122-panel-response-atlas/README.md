# Experiment 122 — Panel Response Atlas

Phase 1 builds a source-only, PDS-null residual for the 300 VCC targets on the
other 299 panel genes.

## Scope

The only effect sources are:

- Replogle K562-GWPS;
- Orion HCT116 and HEK293T, collapsed into one Orion family; and
- the three CD4 compact conditions, collapsed into one CD4 family.

A549, Pan H1, Experiment 116 DLD1 values, VCC held harnesses/outcomes, challenge
truth, `cell-eval2`, count generation, challenge builds, and submissions are
outside this phase.

The preregistration and source authorization were hash-frozen before any
authorized effect matrix was loaded. The first real run later exposed a schema
bug: Orion has duplicate non-panel gene labels. The run stopped before writing
results. `IMPLEMENTATION_AMENDMENT_001.json` records the narrow repair allowing
duplicate non-panel labels while continuing to reject duplicate panel labels;
no candidate, metric, gate, or null changed.

## Method

Each source is reduced to a 300 × 300 target-by-panel-response matrix using
exact gene symbols. Target-self entries are set to exact zero.

Correlated sources are collapsed before validation:

- Orion HCT116 + HEK293T: coordinate-wise median;
- CD4 Rest + Stim8hr + Stim48hr: coordinate-wise median, requiring support in
  at least two conditions and the frozen Experiment 118 target-QC rules;
- K562-GWPS: one family.

Rows are normalized by their family-specific median absolute nonzero effect.
Candidate consensus is computed with equal family weight, then rescaled by the
median contributing family row scale. The mean residual and nonnegative
confidence/ranking signal are separate arrays.

Selection uses leave-one-biological-family-out validation. Its primary
statistics are target-equal signed credit and top-20%-by-confidence,
truth-magnitude-weighted sign precision. A nonidentity candidate must have
positive signed credit and at least 0.5 sign precision in every held family.
Two 200-draw complete-selection nulls rerun the full candidate selection after
coherent target-label or response-label derangement.

## Reproduction

The committed-style workflow is deliberately overwrite-refusing:

```bash
uv run python experiments/122-panel-response-atlas/test_run.py
uv run python experiments/122-panel-response-atlas/freeze_plans.py --verify
uv run python experiments/122-panel-response-atlas/run.py --audit
```

The completed immutable result is under
`work/source-only-phase1-v1/`. Re-running without choosing a new output path
will fail rather than overwrite it.

## Main artifacts

- `work/source-only-phase1-v1/selected-panel-residual.npy`: 300 × 300 mean
  residual.
- `work/source-only-phase1-v1/selected-panel-confidence.npy`: separate 300 ×
  300 ranking signal.
- `work/source-only-phase1-v1/selected-full-residual.npy`: 300 × 18,533
  expansion with every non-panel coordinate exactly zero.
- `work/source-only-phase1-v1/selected-full-confidence.npy`: matching expanded
  confidence signal.
- `work/source-only-phase1-v1/results.json`: complete candidate, fold, null,
  gate, and invariant report.
- `work/source-only-phase1-v1/artifact-manifest.json`: SHA-256 and size for all
  primary result files.

See `RESULTS.md` for the concise outcome.

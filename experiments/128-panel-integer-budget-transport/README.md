# Experiment 128 — Panel integer-budget transport

Status: preregistered source/mechanics audit. No held outcome, `cell-eval2`,
count artifact, challenge build, or submission is authorized.

Experiment 124 correctly showed that naive panel-only effect addition changes
non-panel abundance during whole-profile normalization, but its proposed repair
was not exactly PDS-invisible after production integerization. Production quasi
emission first computes each gene's pseudobulk total as
`floor(400 * L * profile_g)`. Redistributing continuous panel mass can therefore
change the sum of the panel's independently floored totals, changing the
all-gene denominator used by PDS.

This experiment fixes both identified problems:

1. Experiment 122 residuals remain in their declared
   `log1p(50000 * relative_abundance)` units. The continuous transformed value is
   exactly
   `clip(expm1(log1p(50000*q0_j) + eta*r_j), 0, inf) / 50000`.
2. After the production per-target median library `L` is known, the baseline
   integer pseudobulk totals are computed exactly as production does. Every
   non-panel coordinate, target-self, and ineligible panel coordinate is copied.
   The existing mutable-panel integer budget is redistributed with deterministic
   largest remainder allocation, with canonical panel order breaking ties.

Because the grand emitted total and every non-panel integer total are unchanged,
the scorer-normalized non-panel linear profile, log1p profile, and effect are
bit-identical before PDS removes the panel.

The fixed doses `0.10`, `0.25`, `0.50`, and `1.00` are invariant diagnostics
only. They cannot select a dose or authorize production.

Run only after `IMPLEMENTATION_FREEZE.json` exists:

```bash
uv run python experiments/128-panel-integer-budget-transport/test_run.py
uv run python experiments/128-panel-integer-budget-transport/run.py --audit
```

Immutable outputs are written to `work/integration-audit-v1/`.

# Experiment 124 — Panel-simplex transport

Status: **NO-GO after independent audit**. No held outcome, scorer, count
build, or challenge submission was used.

Experiment 122 produced a source-validated residual supported only on the 300
VCC panel genes. That support property is sufficient for PDS invisibility when
effects are compared directly, but production converts log effects to relative
abundance by exponentiating and normalizing the whole 18,533-gene profile.
A naive panel-only log-effect addition can therefore change every non-panel
abundance through the normalization denominator.

This experiment tests two fixed integrations on the exact staged a15 effects
and the released A/B/C control profiles:

1. `naive_log_effect`: add the Experiment 122 residual before the production
   `profile_from_effect` normalization. This is a diagnostic and is expected
   to fail exact non-panel identity.
2. `panel_simplex`: start from the production baseline relative-abundance
   profile, multiply only eligible off-diagonal panel genes by
   `exp(eta * residual)`, and renormalize only those mutable panel genes back
   to their original total mass. Non-panel genes, the target-self coordinate,
   total mass, and nonnegativity are then preserved by construction.

The fixed audit doses are `eta = 0.10, 0.25, 0.50, 1.00`. They are not selected
using the released controls, which contain no perturbation outcome.

Run:

```bash
uv run python experiments/124-panel-simplex-transport/test_run.py
uv run python experiments/124-panel-simplex-transport/run.py --audit
```

The immutable output is written to `work/integration-audit-v1/`, but its
claimed exact post-emission PDS invariant is invalid. See `RESULTS.md`.

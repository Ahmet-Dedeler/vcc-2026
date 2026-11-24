# Experiment 125 — Hierarchical panel posterior

Status: source-only nested evaluation. This directory does not access held
challenge outcomes, Pan, DLD1, H1, the Orion held harness, `cell-eval2`,
count generation, submission construction, submission, git commit, or push.

This experiment compares a low-capacity hierarchical Gaussian
empirical-Bayes/random-effects estimator with Experiment 122's frozen
`equal_mean_agreement` estimator. It uses the exact three already-collapsed
biological families from Experiment 122:

- Replogle K562;
- Orion, with HCT116 and HEK293T already collapsed into one family;
- CD4, with Rest, Stim8hr, and Stim48hr already collapsed into one family.

The latent effect is the posterior mean. Ranking confidence is the posterior
sign confidence `abs(2 * Phi(mean / sd) - 1)` and is never substituted for the
effect. Four posterior settings are crossed with two frozen dose-1 panel
transports: Experiment 124's `q * exp(r)` simplex transport and the exact
production-log1p tilt `expm1(log1p(50000*q)+r)/50000`. Both clip negative
pre-normalization abundance only where mathematically necessary, renormalize
only the mutable panel slice, and must pass exact positivity, target-self,
panel-mass, and PDS-visible-coordinate invariants.

Every outer leave-one-family-out evaluation is nested across five whole-target
folds. For each outer test-target fold, candidate selection uses only the other
four target folds. Each inner validation fold is predicted with hyperparameters
estimated after excluding both the outer test targets and the inner validation
targets. The held biological family is never used for fitting or selection.

The complete procedure is repeated for 200 target-label derangements and 200
within-expression-bin response-label derangements. Uncertainty uses 2,000
target-cluster bootstrap draws.

Commands:

```bash
uv run python experiments/125-hierarchical-panel-posterior/freeze_plans.py
uv run python experiments/125-hierarchical-panel-posterior/test_run.py
uv run ruff check experiments/125-hierarchical-panel-posterior
uv run python experiments/125-hierarchical-panel-posterior/run.py
uv run python experiments/125-hierarchical-panel-posterior/verify.py
```

All outputs are written below
`experiments/125-hierarchical-panel-posterior/work/` and refuse overwrite.

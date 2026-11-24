# 328 — Numerical continuation, not a model or selection change

September 9, 2026, after source v1. Preserve v1's **SOURCE_PILOT_NO_GO**.
The coordinator observed the joint training fit miss tolerance at iteration
500 before reviewing terminal source results. The terminal gate status and
numerical source comparisons have now been reviewed; this is explicitly a
post-result amendment, not a retroactive preregistration.

The failure is numerical: joint mean training log-likelihood is still increasing
by 1.56e-6 per iteration, versus the fixed 1e-6 stopping threshold. Its 45-minute
resource allowance was not exhausted (three fits took 355 s).

## Allowed operation

- Same source archive, labels, rows, gene identities, conditions and complete cases.
- Same 99/8 components, variances, likelihood, optimizer update and uniform-start
  trajectory. **Only continue a fit whose v1 convergence flag is false**,
  from its recorded last mixture weights.
- Same mean-improvement threshold **1e-6**, monotonicity tolerance **1e-9**.
  At most **2,000 total EM iterations**, with at least ten additional iterations
  before stopping. Already-converged control weights remain bit-identical.
- No component pruning, regularizer, SE adjustment, dose, target or gene change.
  The held target values cannot influence iteration count or fitting.
- Before continuation, reproduce every v1 training likelihood endpoint and every
  v1 per-target held density to 1e-10. Freeze the recovery code, this note and
  exact v1 files in `recovery-inputs.json`, and commit before running.
- Evaluate revised weights on the same development targets, report the same
  bootstrap and identical scientific gates. Keep v1 and recovery separate.

The revised numerical gate passes only if every fit satisfies the unchanged
tolerance. Otherwise stop; no further solver extension is authorized here.
Even a numerical and source-gate pass earns only a separately frozen H1
count-stage test, not a package. This is repeated source-development use,
not independent confirmation.

Bound: one process, ≤4 GiB, ≤15 additional minutes, no new response acquisition
or H1 outcome access. Outputs only under `work/recovery/`.

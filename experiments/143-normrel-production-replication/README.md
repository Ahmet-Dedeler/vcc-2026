# Experiment 143 — normrel production replication

This experiment is a one-candidate, fail-closed test of a leakage-safe
development-panel analogue of Experiment 097's target-specific
`normrel^-1.5` donor weighting. It is intentionally limited to the
already-defined 613-target development panel and frozen 4,992-gene response
axis.

The candidate is fixed before any development truth access:

```text
k[d,t] = ||centered donor response[d,t] on 4,992 genes||
         / median over the same 613 rows of that 4,992-gene norm[d]

floor[group] = percentile_10(k[d,t] > 0 over available group slots)
raw_weight[d,t] = max(k[d,t], floor[group]) ** -1.5
weight[d,t] = raw_weight[d,t] / sum_available(raw_weight[:,t])
```

The original K562-GWPS value of `k` controls the complete frozen 50:50
K562/PerturbMatch composite slot. PerturbMatch is not a fourth donor. There is
no alternate exponent, floor, proxy, or blend.

`prepare` is source-only. It validates a separately supplied
`SOURCE_ONLY_PREPARE_AUTHORIZATION.v1.json`, validates every frozen upstream
hash, imports Experiment 142 only after checking its exact source hash, and
materializes:

- Experiment 142's exact folds, target blocks, and 20,000 bootstrap draws;
- the already-frozen 613-row pooled-centered source matrices from Experiment
  142, without reopening any raw donor or held-outcome file;
- equal and candidate donor weights;
- exact equal-weight and weighted effects in Experiment 142's frozen
  PerturbMatch-aware, pruned-whitener geometry;
- a copied, hash-bound development truth ledger without reading truth arrays.

PDS and signed credit are invariant to the positive global alpha, global donor
scale, and post-whitening preserve scalars omitted from this geometry.
None of the 613 development targets has complete support in all three CD4
conditions, so this panel cannot measure the challenge-specific interaction
with a15's CD4 residual. That interaction requires a separate source-only
challenge-panel audit only if this development gate passes.

`develop` requires a different, separately supplied
`DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json`. It writes the persistent outcome
access lock before validating or loading any truth bytes, then loads only the
four 613-target development matrices through the bound Experiment 142/132
loader. It computes full-gallery PDS and signed-credit deltas, three-family and
five-fold cells, and the frozen paired block bootstrap.

The candidate passes only if both metrics satisfy all frozen raw gates:

- macro delta at least `0.003`;
- paired-bootstrap 95% lower bound strictly above zero;
- every family and every family-by-fold cell nonnegative;
- positive-gain family concentration at most `0.5`;
- positive-gain target-block concentration at most `0.5`.

Any failure is `TERMINAL_NO_GO`. This experiment has no null, sealed
confirmation, count-generation, `cell-eval2`, challenge-build, or submission
stage.

Commands, after independent authorization files have been supplied:

```bash
.venv/bin/python experiments/143-normrel-production-replication/run.py audit
.venv/bin/python experiments/143-normrel-production-replication/run.py prepare
.venv/bin/python experiments/143-normrel-production-replication/run.py develop
```

Implementation checks only:

```bash
.venv/bin/python -m unittest \
  experiments/143-normrel-production-replication/test_run.py
```

Neither `prepare` nor `develop` is run as part of implementation.

This is not an exact replay of Experiment 097's full-donor-panel key: strict
410-target isolation prevents recomputing that denominator. The frozen
`normrel613` analogue keeps numerator and denominator on the same allowed axis
and rows; any positive result would still require independent confirmation.

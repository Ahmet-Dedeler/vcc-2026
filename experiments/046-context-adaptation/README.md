# 046 — Supervised target-context adaptation: no-go

This experiment asks whether a small amount of public, target-context
Perturb-seq data can teach a generic donor consensus how that context responds.

For each of K562-gwps, Orion-HCT116, and Orion-HEK293T:

1. Average the other two donors.
2. Fit a low-capacity map to the held target on 600 common, non-VCC
   perturbations.
3. Select hyperparameters with five perturbation-held-out folds.
4. Evaluate once on 259 shared VCC-panel perturbations.

Candidate maps are scalar and diagonal shrinkage, randomized-PCA ridge maps,
and residual-program regressions. No full gene covariance or single-cell scorer
is used. All target coordinates from both the training and VCC panels are
excluded before fitting a map.

## Verdict

The apparent gain was not target-context adaptation. A paired rank-32 ridge map
improved macro PDS by +0.0421, but permuting target perturbation labels improved
it by +0.0566 and wrong-context maps by as much as +0.0579 macro. The gain is a
generic geometry effect from the target effect subspace, not learned
perturbation-specific context biology.

Do not apply the generated Jurkat artifact. It is retained only for audit
reproducibility and marked as quarantined.

The strict gate for fitting an analogous generic-to-Nadig-Jurkat artifact is:

- target-specific CV selection improves every outer donor;
- donor-nested selection improves every outer donor; and
- both macro gains are at least 0.002 raw PDS.
- the paired map beats a target-label permutation and every wrong-context map
  in each outer donor.

The final condition failed in all three donors.

Run under the repository memory guard:

```bash
python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 4 --min-free-percent 35 --poll-seconds 2 -- \
  uv run python experiments/046-context-adaptation/run.py
```

Negative controls:

```bash
python3 scripts/run_with_memory_guard.py \
  --max-rss-gib 2 --min-free-percent 35 --poll-seconds 2 -- \
  uv run python experiments/046-context-adaptation/audit_controls.py
```

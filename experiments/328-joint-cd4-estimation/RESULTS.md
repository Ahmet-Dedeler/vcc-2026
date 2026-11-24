# 328 — Source pilot v1: positive held-target density, numerical NO-GO

September 9, 2026. Original registered status: **SOURCE_PILOT_NO_GO**.
No H1 outcomes, count scoring, package, submission, new reserve or paid resource.

The joint covariance mixture predicts the held CD4 source observations better
than its two controls. However, the joint EM fit reached its 500-iteration
limit at mean log-likelihood improvement **1.56e-6**, above the registered
**1e-6** threshold. It is not recorded as a scientific pass.

## Registered source result

| comparison | mean log-density delta (nats/triple, target-macro) | paired target 95% CI | fraction of targets favoring joint |
|---|---:|---|---:|
| joint minus diagonal | +.0132563 | [+.0102059, +.0167506] | .704 |
| joint minus shuffled_joint | +.0222034 | [+.0140805, +.0354749] | .732 |

All three held-target gates pass; the all-fits-converged gate fails.
Joint: 500 iterations, not converged. Diagonal: 143, converged.
Shuffled joint: 500, converged (last improvement just below the threshold).
The results are in `work/results.json`; do not overwrite or relabel them.

This is predictive density of observed DESeq2 coefficients conditional on
reported SE, **not latent biological accuracy or competition score**.
Related perturbation targets, one source study, fixed 256 genes and three
conditions limit generalization. Target shuffling preserves gene identity,
the (effect, SE) pairing and every complete-case marginal.

## Source and engineering checks

- Metadata selection frozen at **8c7c4f0a** before reading its effect layers.
- 10,737 unambiguous complete-condition, non-panel target IDs available;
  fixed first 1,000 hash-ordered IDs for fitting and next 250 for development.
- Both A/B/C and H1 targets excluded by symbol/case-fold and mapped Ensembl ID;
  their gene coordinates also excluded.
- **320,000/320,000 complete triples**: 256,000 fit, 64,000 development.
  Every selected target has all 256 triples; no support exclusion occurred.
- Each condition's fixed first compact-bank row exactly reproduces parent
  `log_fc`, `lfcSE` and float32 natural-log delta. Seven synthetic tests pass.
- Extraction: 6.86 s, sampled process-tree peak .19 GiB.
  Three fits plus report: 355.15 s, peak 2.40 GiB. All exited zero.
- Source archive SHA-256:
  `ad14bdd98ae8c89303a96bec85271657bf8feebe56348dab582e919a7749fff3`.
- Input manifest SHA-256:
  `d84c3a137b5964d37316593f86b877bb2f6ec3d4bb131ed89a6b21e95908e54b`.

## Numerical follow-up is separate

`NUMERICAL-RECOVERY.md` defines a solver-only continuation, motivated by the
training convergence failure. The original data, model/covariance grid,
likelihood, target split and held-target gates do not change. This is a
**post-result numerical amendment**: v1's gate status and source outcomes
were inspected, and any revised development report is not an untouched
confirmation. No H1 count test is authorized by v1 itself.

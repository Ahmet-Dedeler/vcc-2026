# 085 — Frozen PM inclusion gate on untouched KOLF

Date: **2026-08-23**

Status: **completed; NO-GO at the source-only training gate.**

No KOLF score was computed. The fixed training gate selected the exact
`always_D` identity fallback, so the runner stopped before the durable freeze,
global access marker, or first KOLF file authentication. The experiment's
one-time KOLF evaluation therefore remains unused. No count artifact or
challenge submission was created.

## Result

The run used the exact source tree pinned by the preregistration at commit
`2d5cd1dd9382af89ef43b61ec98140b999887b26`.

| rule | HCT116 delta | HEK293T delta | macro delta | eligible |
| --- | ---: | ---: | ---: | --- |
| `always_D` | 0 | 0 | **0** | yes |
| best nonidentity: `high:q=0.20` | +0.000432 | -0.001255 | -0.000411 | no |
| `low:q=0.20` | -0.002674 | -0.000123 | -0.001399 | no |

The full-selection permutation test for the selected identity rule gave
`p = 1.0` with null p95 `0.001810`; the identity fallback was selected in
71.4% of 5,000 null selections. No nonidentity rule had nonnegative gain in
both source donors, and none met the preregistered `+0.0015` macro threshold.

All experiment-064 arm reconstructions were exact (`maximum_arm_c_pds_error =
0`, `maximum_arm_d_pds_error = 0`) on both held donors. Peak RSS was
`4.043 GiB`.

Machine-readable artifacts:
`cloud-output/20260830T003317Z/results.json`.

## Operational provenance

The first 2026-08-30 launch attempt stopped locally before dispatch because
the historical S3 location for the two held-whitener artifacts was absent.
The local copies matched the frozen artifact manifest exactly and were
restored to that original location. A fresh run then completed.

Current HEAD could not be used directly: three shared source files had moved
since the experiment was frozen, and the current builder also permits ambient
`VCC_*` variables to change the arms. The launch therefore came from a detached
checkout of the exact pinned commit rather than changing expected hashes.

## Implemented frozen path

- One feature only: midrank percentile of the exact arm-C/arm-D angular
  displacement on each frozen whitening-operator axis.
- One model family only: low/high monotone stumps at
  `q ∈ {0.20, 0.30, 0.40, 0.50, 0.60}`, plus mandatory always-D fallback.
- Training and all candidate selection use only the 221-target HCT116 and
  HEK293T experiment-064 outcomes.
- A 5,000-replicate joint target-label permutation reruns the complete
  candidate selection before KOLF can be opened.
- Production builds exact three-donor gamma-0.65 arms C and D for all 300
  targets. Percentile ranks are computed only among the 248 joint
  K562/PerturbMatch-supported targets; unsupported targets are forced to exact
  D.
- Arm C, arm D, assignments, gated queries, source counts, and source-only
  clusters are atomically serialized and hashed, then every sealed file is
  durably uploaded with no-overwrite semantics before KOLF access is claimed.
- The KOLF access marker is one atomic global S3 claim, not merely a per-run
  object. Its content is read back for SHA-256 verification before any KOLF
  file is authenticated, preventing a cross-worker double evaluation.
- KOLF is then restricted to the fixed complete 273-target gallery and fixed
  17,806-gene non-panel score axis. Its all-ones `b_ctrl.npy` is hashed for
  provenance but never loaded as a feature.
- The only post-KOLF operations are fixed scoring, the preregistered
  source-cluster bootstrap, the fixed-count matched-routing null, source-strata
  summaries, and promotion checks. There is no fit, feature choice, threshold
  choice, or fallback change after KOLF opens.

## Verification completed

- All **57 non-KOLF dependency artifacts** passed their pinned SHA-256 checks,
  including the two held experiment-064 operators, the production operator,
  donor/PM/CD4 banks, challenge axes, frozen helper code, and exp064 outputs.
- Experiment-064's result contract was checked: 221 targets, both directional
  folds, 17,942-gene operators, 17,806-gene retrieval axes, panel exclusions,
  candidate-independent operators, untransformed truth, and no retuning.
- **24 synthetic checks passed** under a 1 GiB process-tree guard.
- Peak synthetic-test RSS: **0.260 GiB**.
- Peak non-KOLF hash-preflight RSS: **0.109 GiB**.
- Python compilation, Ruff checks/formatting, JSON parsing, shell syntax, and
  executable-bit checks pass.
- Final adversarial review caught and fixed two fail-open risks before sealing:
  a cross-worker race in the one-time S3 marker and non-finite vectors being
  able to receive a spurious perfect PDS rank.
- The CLI exposes no feature, direction, quantile, threshold, seed, null-size,
  cluster, or operator-hash override.
- The AWS path requires explicit clearance, an already-running SSM-online
  worker with at least 64 GiB nominal and 48 GiB available memory, a global
  heavy-process lock, a 44 GiB systemd memory ceiling with swap disabled, no
  conflicting process/service, fresh work and S3 prefixes, pinned source
  archive and package versions, and no prior KOLF access marker. It contains
  no instance-create/start or VCC submission command.

## Original launch recommendation

The pre-run recommendation was **GO for exactly one cloud run after explicit
approval**. The completed run demonstrated that the fail-closed design worked:

1. if the HCT116/HEK293T training gate fails, it writes a rejection result and
   never authenticates or opens KOLF;
2. if the training gate passes, it freezes and hashes the complete production
   routing first, then consumes the one-time KOLF evaluation;
3. promotion still requires `ΔPDS ≥ 0.0015`, bootstrap lower 95% above zero,
   matched-routing `p ≤ 0.05`, both routes represented, exact D for unsupported
   targets, and all provenance/axis checks.

The empirical question was resolved before held-truth access: the source-only
feature does not support routing away from PerturbMatch. Close this rule
family unless a new external feature or mechanism is introduced.

KOLF is not globally unseen within the repository: prior disjoint experiments
used it for other hypotheses. The defensible claim is narrower and explicit:
experiment 085's gate, candidate family, and promotion criteria were frozen
without reading KOLF outcomes for this PM-routing question.

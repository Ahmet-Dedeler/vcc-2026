# 331 — TabICL synthetic-only local execution feasibility

September 9, 2026. This is an engineering feasibility test, not a biological
experiment, a response-model selection or a VCC score prediction.
It follows the read-only official-source check in
`docs/sprint-2026-09-09/MODEL-FEASIBILITY.md`.

## Fixed materials

- TabICL 2.2.0, official PyPI wheel:
  SHA-256 `c80b8ec2719be26bc0bf1dd758a2361f58019d8f7877038b785529a425789ce4`,
  310,448 bytes, published September 2, 2026.
- Official `jingang/TabICL` model revision
  `4dcd344ece2c00be9e831fdd35bed57b5ad83e19`,
  `tabicl-regressor-v2-20260212.ckpt`, **114,324,594 bytes**,
  SHA-256 `0db9cb538f114e79026bf08f45f41ad8dd7ad2de2aaca9a5ca8cd3bd9748ae7a`.
  Exact byte size/hash verified before any model load.
- Code and model card declare BSD-3-Clause; retain source attribution.
  Main code snapshot inspected:
  `soda-inria/tabicl@0dbff3ec8fc68c123c87af77b0ea8b25cd2d23f3`.
- Resolve dependencies into a hash-pinned requirements lock, install only in
  this experiment's `.venv`. Never edit the project's environment or `uv.lock`.
  Pin code, protocol, dependency lock and root `uv.lock` hash before inference.

The checkpoint is synthetically pretrained, not a biological donor. No
challenge, H1 or public perturbation response is loaded by this test.
Checkpoint download is public and free. No API model call, account grant,
paid service, cloud, telemetry or automatic upload.

## Fixed workload

CPU only for this first check, four Torch threads, float32/no AMP, one
estimator, batch size one, normalization `none`, no KV cache, seed 20260909.
No default automatic checkpoint download during inference.
The MPS backend is not tested or claimed operational by this protocol.

Two independent deterministic synthetic problems, 512 and 2,048 training
rows, 256 test rows, 12 uniform [0,1] features:

`y = 10*sin(pi*x0*x1) + 20*(x2-.5)^2 + 10*x3 + 5*x4 + N(0,.1^2)`.

For each size:

1. Fit once and predict mean/median; repeat mean prediction to measure
   reproducibility and warm runtime.
2. Fit the same model on a fixed permutation of training labels and predict.
3. Report R2 relative to the training-mean null, MSE, timing and peak process
   memory. The shuffled fit checks that supplied labels affect predictions.

Engineering pass requires finite, correctly shaped outputs; maximum repeat
mean difference <=1e-6; positive R2 versus the training-mean null on both
sizes; real-label MSE below shuffled-label MSE; and the registered resource
bound. This is a synthetic sanity gate, not an accuracy guarantee or a
statistical biological result.

## Bound and decision

One heavy process, repository launch gates, <=12 GiB RSS and <=4 GiB scratch.
Setup plus measurements <=45 minutes; each size <=15 minutes. Record failures
and stop, not a device/ensemble/normalization sweep. Use no thread/subagent.

A pass only permits consideration of a separately registered biological
pilot with a genuinely distinct estimand and strict held-context controls.
It does not reopen 329 or 330, and does not authorize packaging or submission.

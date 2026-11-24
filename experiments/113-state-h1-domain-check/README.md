# 113 — STATE H1 domain check

Status: **complete — STATE is near random on H1 and does not justify further retraining**.

Date: **2026-08-30**.

## Question

Did experiment 106 reproduce a checkpoint that learned perturbation-specific
signal at all, or is its severe Orion failure only an out-of-context transfer
failure?

The fixed experiment-106 `state_sm` step-20,000 checkpoint is evaluated on the
85-target H1 test panel. H1 is the same biological context as the 2025
competition training data but uses held perturbations, so this is a
target-generalization sanity check rather than evidence for 2026
context-generalization.

## Frozen inputs

- experiment-106 immutable `state_sm` step-20,000 checkpoint;
- `data/harness/h1-test.h5ad`, using only its 18,400 non-targeting rows during
  inference;
- experiment-100's already-frozen H1 identity prediction as the incumbent
  comparator;
- experiment-098's H1 reference-DE bundle and `data/donors/h1_2025` truth,
  opened only by the separate score command.

## Inference

- exact STATE 18,080-gene order;
- four deterministic 128-cell control banks;
- `log1p(raw counts)` input, matching the untransformed log-count support
  matrices used for training;
- one homogeneous target per forward pass;
- clip model outputs to `[0, 14]`, exponentiate, average in count space, then
  construct a 50,000-count log-pseudobulk effect;
- preserve each bank separately and also aggregate all 512 predicted cells in
  count space;
- set the union of the 300 VCC-2026 panel genes and all 85 H1 test target
  coordinates to exact zero before scoring.

No perturbation-expression row is read during inference.

## Fixed diagnostics and decision

Report:

1. pure STATE PDS, cosine, and uncorrected effect NMSE;
2. the same metrics for the frozen incumbent;
3. four individual-bank PDS values and their aggregate;
4. directional credit on the fixed experiment-098 H1 DE gate;
5. the oracle directional-credit headroom obtainable by replacing only
   incumbent/STATE sign disagreements after matching each STATE LFC row to the
   incumbent row norm.

Interpretation is binary:

- If pure STATE PDS is materially above random (`>0.55`) and perturbation
  labels beat a frozen row-permutation null, experiment 106 learned an
  in-domain perturbation signal and the Orion result is evidence of failed
  context transfer.
- Otherwise the current checkpoint/reproduction is invalid as a scientific
  basis for experiment 110, and no RPE1-out or Jurkat-out retraining is
  justified.

The check is diagnostic only. It authorizes no challenge build, submission,
count-level scorer, checkpoint selection, or post-result tuning.

## Commands

```bash
STATE_USE_MPS=1 PYTORCH_ENABLE_MPS_FALLBACK=0 \
  data/external/state-vcc/state/.venv/bin/python \
  experiments/113-state-h1-domain-check/run.py freeze

uv run python experiments/113-state-h1-domain-check/run.py score
```

# Experiment 118 Pan confirmation — one-shot evaluator

Status: **GO for independent review; Pan values unopened by this evaluator**.

This is a separate confirmation stage for the already frozen
`unweighted_eta1p00` residual from Experiment 118. It has no candidate grid,
no tuning, and no reselection. Its only outcome source is the hash-pinned
Experiment 115 Pan2025 H1 compact donor.

## Safety boundary

`authenticate`:

- verifies every Experiment 118 source-input and frozen-output hash;
- requires the frozen `selection.json` hash
  `3d849a0fb1227b3e968861fb891d954513a73de12c80130a2b20059f880028d4`;
- verifies the Experiment 115 compact receipt, artifact manifest, axes, control,
  delta, and measured-mask hashes;
- parses only JSON and text axes;
- hashes binary files as opaque bytes and never calls `numpy.load`;
- requires both the permanent access lock and result directory to be absent.

`run` repeats authentication, requires the exact confirmation token, creates
and fsyncs an exclusive permanent access lock, and only then enters the single
Pan-value-loading phase. The lock is never removed, including after failure.
Results are written once with exclusive creation, fsynced, hash-manifested,
and made read-only.

No command in this evaluator generates counts, invokes `cell-eval2`, builds or
submits a challenge artifact, commits, or pushes.

## Frozen candidate and exact intersection

The only permitted candidate is `unweighted_eta1p00`.

Pan target rows are matched to the canonical 300-target axis by exact symbol.
Response coordinates use the exact canonical panel-gene intersection and the
Pan row's measured mask. The destination target's own coordinate is always
excluded. No aliases, case folding, synonym mapping, imputation, or fallback
candidate are allowed.

## Linear-mean application

The Experiment 118 residual is **not exponentiated**. Its numeric values are
interpreted directly as additive units on the Pan TS=50,000 linear mean:

```text
control_linear = expm1(Pan b_ctrl)
candidate_linear = control_linear + intersected_centered_residual
```

Exponentiation is used only to reconstruct Pan's released truth:

```text
truth_linear = expm1(Pan b_ctrl + Pan delta)
truth_delta_linear = truth_linear - control_linear
```

The intersected residual is recentered over mutable off-diagonal panel
coordinates. No dose scaling, clipping, or positivity repair is allowed.
Negative candidate means are an automatic no-go.

The run must prove:

- exact non-panel identity;
- exact target-diagonal identity;
- exact panel mass conservation under compensated summation;
- exact identity of all PDS-visible normalized linear and log1p coordinates.

## Metrics and gates

Each of the 64 expected exact-symbol Pan/VCC target intersections receives
equal weight.

- Signed credit per target:
  `sign(truth_delta) @ residual / sum(abs(truth_delta))`.
- Top-20% precision: rank eligible responses by absolute residual, select
  `ceil(20%)`, score sign agreement weighted by absolute truth magnitude.

Eligibility is Pan-measured, finite, nonzero-truth, off-diagonal panel genes.

The evaluator freezes 20,000 whole-target bootstrap draws and evaluates the
already frozen 100 target-label and 100 within-expression-bin response-label
permutations without reselection. GO requires every threshold in
`PAN_CONFIRMATION_PROTOCOL.v1.json`, including positive credit, credit
bootstrap lower bound above zero, precision above 0.55, precision bootstrap
lower bound above 0.50, and both observed metrics strictly above both null
95th percentiles.

## Commands

Truth-safe authentication:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -B \
  experiments/118-pds-null-panel-residual/evaluate_pan_confirmation_v1.py \
  authenticate
```

Synthetic and authentication-safety tests:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -B -m unittest \
  experiments/118-pds-null-panel-residual/test_pan_confirmation_v1.py
```

One-shot value access, for an independent reviewer only:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -B \
  experiments/118-pds-null-panel-residual/evaluate_pan_confirmation_v1.py \
  run --confirm-pan-access I_CONFIRM_EXPERIMENT_118_PAN_CONFIRMATION_V1
```

Do not run the final command twice. The permanent lock intentionally makes a
second execution impossible.


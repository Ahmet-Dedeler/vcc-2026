# Experiment 130 — TeloHAEC panel confirmation

Status: **pre-value freeze complete; biological value access locked**

This experiment preregisters one independent, 65-target TeloHAEC confirmation
of a future exact Experiment 125 panel candidate. It is a distinct mechanism
from Experiment 105's failed broad direct-donor transfer: only off-diagonal VCC
panel responses are evaluated.

The current phase is metadata-only. It may inspect file structure, labels,
guide assignments, and already-authorized cell counts, but it must not read
TeloHAEC `X`/layer expression values or any candidate biological arrays.

Two-stage boundary:

1. This directory freezes the source receipt, axes, statistics, nulls,
   bootstrap, Experiment 128 transport dependency, implementation, and gates.
2. A later `CANDIDATE_BINDING_AMENDMENT.json` must validate before any
   TeloHAEC panel value is read. It must bind an immutable Experiment 125
   `GO` result and the exact modal candidate selected by the frozen rule,
   including residual, confidence, support, axes, and hashes. The candidate
   must use `production_log1p` transport.

No candidate binding exists now. `run.py --check-lock` and the synthetic tests
prove that value access fails before any protected file is opened.

Metadata-only commands:

```bash
uv run python experiments/130-telohaec-panel-confirmation/run.py --metadata-preflight
uv run python experiments/130-telohaec-panel-confirmation/test_run.py
```

The future value command is intentionally unusable until the amendment exists:

```bash
uv run python experiments/130-telohaec-panel-confirmation/run.py --evaluate
```

Any scientific or implementation repair after the first TeloHAEC value read
invalidates Experiment 130. It must receive a new experiment number.

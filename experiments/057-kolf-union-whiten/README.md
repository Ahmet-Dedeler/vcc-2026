# 057 — KOLF under source-only union whitening

This experiment asks one narrow question: does experiment 055's already-fixed
raw KOLF rule (`fraction=0.1`) still help after the production-style union
whitener?

Protocol:

- target-role holdout: K562-GWPS;
- predictive sources and whitening basis: Orion HCT116 + Orion HEK293T only;
- auxiliary prediction: KOLF2.1J raw effect, fixed fraction `0.1`;
- target gallery: every VCC-panel perturbation shared by all four banks;
- whitening axis: source-gene union on the official VCC axis, after removing
  all 300 VCC target-gene coordinates;
- primary gamma: `0.65`, fixed before opening K562 effect values;
- diagnostics: `0`, `0.55`, `0.75`, and `0.85`, never used as selectors;
- score axis: every non-panel K562-observed VCC gene. The exact experiment-055
  four-bank shared axis is also reported as a checksum;
- query-fold summaries partition query rows only. Every query is always ranked
  against the complete shared-target gallery.

The operator uses a deterministic production-sized sample of 6,000 background
perturbations per source donor. K562 and KOLF contribute zero rows to its
covariance. The unwhitened shared-axis result must reproduce experiment 055
within `1e-9` or the run fails.

Real-data computation belongs on AWS. Local work is limited to the synthetic
tests:

```bash
uv run --frozen pytest -q experiments/057-kolf-union-whiten/test_run.py
```

No count emission, `cell-eval2`, submission build, or challenge submission is
performed.

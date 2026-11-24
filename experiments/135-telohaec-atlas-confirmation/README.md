# Experiment 135 — TeloHAEC confirmation of the frozen panel atlas

Status: **independently reviewed; ready for immutable pre-value freeze**

Experiment 122 selected one source-only 300 × 300 panel-response atlas,
`equal_mean_agreement`, from K562, Orion, and CD4. Experiment 128 later
established an exact integer-budget transport that changes only PDS-excluded
panel coordinates while preserving every non-panel emitted total, target-self
total, panel budget, and grand total.

This experiment performs one independent biological confirmation of that
exact frozen candidate on the 65 TeloHAEC panel targets whose panel responses
have never been aggregated, inspected biologically, or used for selection.
Experiment 105 previously streamed every raw count chunk while excluding all
300 panel coordinates before aggregation; this is a fresh semantic outcome,
not a claim that the source bytes were never read. There is no candidate grid
and no dose selection:

- candidate: Experiment 122 `equal_mean_agreement`;
- comparator: the exact staged-a15 effect array;
- transport: Experiment 128 integer-budget transport;
- dose: `1.0`;
- target aggregation: median of all eligible guide effects;
- controls: pooled pseudobulk of all eligible control cells, matching Experiment 128;
- metrics, nulls, bootstrap seeds, and promotion gates: fixed in
  `preregistered_protocol.json`.

The first semantic pass over the H5AD reads only eligible control-cell rows.
It maps the exact staged-a15 effects onto the collapsed TeloHAEC axis, replays
all 300 target-specific production libraries from one `default_rng(0)`,
integerizes the 65 incumbent rows, then freezes control-abundance response
bins and both unrestricted exchangeable null banks. Only after those artifacts
are hashed does a second pass aggregate targeting-guide panel responses.

For every target, the baseline is the exact staged-a15 incumbent in the same
TeloHAEC control and integer-library frame. The candidate is that same row
after the exact Experiment 128 transport. The two primary endpoints are
candidate-minus-incumbent signed credit and candidate-minus-incumbent
top-20%-magnitude sign precision. Experiment 122 support is provenance-only:
every finite residual coordinate follows Experiment 128 mutability semantics.

The target-label null permutes the 65 rows and conjugates the same unrestricted
permutation onto their 65 response labels, preserving the source diagonal at
the destination self coordinate. The response-label null is row-specific: it
fixes self and randomly permutes all other coordinates inside frozen pooled
control-abundance bins. Target cell totals and held-out-guide-versus-remaining
guide diagnostics are always summarized. After target access, the evaluator
releases aggregate results and hashes only—not truth, per-target metrics, null
samples, or bootstrap samples. The output directory is overwrite refusing.

A pass authorizes only a separate production implementation and full artifact
audit. It does not authorize a challenge build or submission.

Pre-value checks:

```bash
uv run python experiments/135-telohaec-atlas-confirmation/test_run.py
uv run python experiments/135-telohaec-atlas-confirmation/freeze.py
```

One-shot value access:

```bash
python scripts/run_with_memory_guard.py \
  --max-rss-gib 4 --min-free-percent 30 --poll-seconds 1 -- \
  uv run python experiments/135-telohaec-atlas-confirmation/run.py \
    --evaluate \
    --confirm-value-access I_CONFIRM_EXPERIMENT_135_TELOHAEC_PANEL_V1
```

Forbidden: changing the candidate, dose, comparator, support rule, metric,
null, seed, threshold, or implementation after the first targeting-guide
panel response is aggregated;
using HCASMC, DLD1, Pan, KOLF, H1/Orion held outcomes, challenge truth,
`cell-eval2`, count generation, a challenge build, or submission.

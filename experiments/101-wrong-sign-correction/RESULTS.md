# 101 — cross-fitted selective wrong-sign correction

**Verdict: NO-GO.**

The experiment used the preregistered rank-16 signed ridge grid, strict outer-donor holdout, four whole-target folds, 2,000 target-block bootstrap draws, and 20 complete target-label permutation nulls. No count-level scorer, AWS worker, challenge artifact, submission, or commit was used.

## Selected configurations and outer results

| surface | selected | linear NMAE identity | flip | zero | flip credit gain | PDS delta | tail error | selected pairs |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| H1 | r=0.1, tau=0.50, lambda=1 | 0.93062 | 0.93059 | 0.93062 | +0.00003 | +0.00014 | 0.44824 | 26 |
| HCT116 | identity | 0.87011 | 0.87011 | 0.87011 | +0.00000 | +0.00000 | n/a | 0 |
| HEK293T | identity | 0.99512 | 0.99512 | 0.99512 | +0.00000 | +0.00000 | n/a | 0 |

## Macro uncertainty and null

- Flip linear-credit gain: +0.00001 [-0.00005, +0.00006].
- Zero-tail linear-credit gain: -0.00000.
- Flip minus zero: +0.00001 [-0.00001, +0.00003].
- Flip PDS delta: +0.00005 [-0.00065, +0.00075].
- Selected-surface weighted tail error (1/3 surfaces: h1): 0.44824 [0.34604, 0.84519].
- Permutation null p50 / p95: +0.00000 / +0.00055; empirical one-sided p=0.0952.

## Preregistered decision

| condition | result |
| --- | --- |
| linear gain >= 0.003 on every surface | **FAIL** |
| macro bootstrap lower bound above zero | **FAIL** |
| tail error > 60% on every surface | **FAIL** |
| PDS nonnegative on every surface | PASS |
| observed gain above permutation p95 | **FAIL** |
| flip beats matched zero-tail control | PASS |
| overall | **NO-GO** |

## Interpretation

The proposed mechanism is not supported. H1 selected only 26 pairs across 10 targets, and just 44.8% of their metric weight was actually wrong-sign. In other words, the selected tail was more often correct than wrong, so reversal was not a reliable repair. HCT116 and HEK293T found no nested PDS-safe candidate with positive signed-credit gain and therefore selected identity.

The tiny positive H1 result is unresolved: the macro gain is roughly 270 times smaller than the preregistered 0.003 per-surface threshold, its bootstrap interval crosses zero, and it falls below the 20-draw permutation-null p95. Although flip is numerically above the matched zero control, that difference is also tiny and its interval crosses zero.

See `work/summary.json` for the complete machine-readable result and each surface's `selection.json`, `freeze-audit.json`, `metrics.json`, and `per_target_metrics.npz` for the full selection, isolation, bootstrap, and null evidence.

Inner selection intentionally used an unwhitened, median-norm-scaled base because some available precomputed whiteners contain the outer held donor. Outer evaluation retained experiment 100's exact surface-specific whitened base.

## Provenance limitation

The freeze audit recorded hashes for predictions, selections, and null-model arrays, but the scorer did not enforce those hashes before opening held truth and did not freeze source-code or whitener hashes. `run.py` was changed after H1/HCT116 scoring to repair the HEK293T identity-null crash. The change did not alter the selected H1 model or either identity result, and every stored number reproduces from the frozen arrays, but code immutability across all three score stages cannot be proven retroactively. This weakens provenance, not the NO-GO: the corrected conservative null threshold is larger and the result still fails every substantive promotion gate.

Each permutation draw also creates one deterministic mapping per whole-target fold after that fold's targets are removed, rather than one global mapping reused across folds. This is leakage-safe but should be made explicit in any reuse.

## Reproduction

```bash
uv run python experiments/101-wrong-sign-correction/test_run.py
uv run python scripts/run_with_memory_guard.py --max-rss-gib 28 -- uv run python experiments/101-wrong-sign-correction/run.py prepare
uv run python scripts/run_with_memory_guard.py --max-rss-gib 28 -- uv run python experiments/101-wrong-sign-correction/run.py freeze --surface h1
uv run python scripts/run_with_memory_guard.py --max-rss-gib 28 -- uv run python experiments/101-wrong-sign-correction/run.py score --surface h1
# Repeat freeze/score sequentially for hct116 and hek293t.
uv run python experiments/101-wrong-sign-correction/summarize.py
```

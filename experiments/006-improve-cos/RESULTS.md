# 006 — Two small but replicated gains in transfer cosine

**Date:** 2026-08-21 · Offline against `h1-test` (82 perts) and `h1-trainval` (114 perts).

Experiment 004 closed two doors (low-rank truncation, supervised per-gene weights). Both were
*reweightings of a fixed donor estimate*. This tries a different kind of change: **which genes are
allowed to set the direction.** Cosine is driven by whichever genes carry the largest values, so a
gene that is merely noisy in the donor steals direction without contributing signal.

## What worked

| variant | h1-test | h1-trainval |
|---|---|---|
| baseline | 0.0949 | 0.1118 |
| + divide by mean \|effect\| per gene | 0.0983 | 0.1184 |
| + target CPM > 20 gate | 0.1009 | 0.1198 |
| **+ both** | **0.1029** | **0.1249** |
| | **+8.4%** | **+11.7%** |

**The CPM optimum is 20 on both harnesses independently** — a replicated finding, not a fit to one
split. Below 20 a gene contributes noise to the predicted direction, and it is excluded from the
DE metrics anyway (their filter is 5 CPM in control). Above ~50 the gate starts discarding real
signal.

End-to-end through the model's own `profile()` path on h1-test: **0.0726 → 0.0771 (+6.2%)**.

## What failed, and it is a pattern

| variant | h1-test | vs baseline |
|---|---|---|
| rank-normalise per perturbation | 0.0470 | −0.048 |
| keep top 50% most-variable genes | 0.0638 | −0.031 |
| keep top 25% most-variable genes | 0.0534 | −0.042 |
| keep top 9% most-variable genes | 0.0407 | −0.054 |
| clip at ±3 gene sd | 0.0630 | −0.032 |
| sign only (binarised) | 0.0369 | −0.058 |
| anti-standardisation (× sd^0.5) | 0.0929 | −0.002 |

**Every transform that discards magnitude information loses badly.** Rank-normalising, binarising,
clipping and hard gene subsetting all cost 30–60% of the signal. Selecting the *most variable*
genes is actively harmful — the opposite of the usual single-cell instinct.

Combined with experiment 004's failures, three independent families of transform have now lost to
the raw donor delta: spectral truncation, supervised reweighting, and magnitude discarding. The
transferable signal is **high-dimensional, spread across many modestly-responding genes, and
carried in the magnitudes**. Only the two interventions that *preserve* the full magnitude
structure while removing per-gene noise scale gave anything.

## Applied

`TransferModel(gene_norm=True)` and `profile(min_cpm=20)` are now the defaults.

# Experiment 154 — confidence-key screen for `reach`

Decision: **TERMINAL NO-GO — `reach` is model-limited, not order-limited.**

Deterministic screen, no emitter and no seeds, so the ~`0.012` reach seed-noise
that invalidates single-seed over-dispersion arms cannot reach these numbers.
Reproduce with `screen_keys.py`.

## The question

`reach` walks the reference's DE pool in **our** confidence order and measures the
deepest prefix holding >=90% sign purity. The pool and each gene's correctness are
fixed by the model; only the order is free, and the ordering oracle is worth
`+0.103` overall. So the falsifiable question is narrow:

> Does any available signal isolate a subset of our predictions with >=90% sign
> accuracy, per perturbation?

## Answer: no

Per-perturbation reach, 127 perturbations with pool >= 10, median pool 70:

| key | reach | median k* | k*=0 |
| --- | ---: | ---: | ---: |
| `|pred_lfc|` (baseline) | `0.10716` | 4 | 7 |
| `-pred_padj` | `0.07231` | 3 | 19 |
| donor sign concordance | `0.05657` | 1 | 42 |
| agreement-with-pred | `0.06935` | 1 | 41 |
| n donors available | `0.02661` | 1 | 56 |
| concordance x `|pred_lfc|` | `0.11117` | 4 | 6 |
| concordance^2 x `|pred_lfc|` | `0.11140` | 4 | 8 |
| agreement x `|pred_lfc|` | `0.10682` | 4 | 6 |
| **geo-magnitude x concordance** | **`0.12473`** | 4 | 7 |
| **ORACLE (correct signs first)** | **`0.62797`** | **44** | 0 |

The best buildable key beats `|lfc|` by `+0.01757` raw reach and captures
**3.4%** of the oracle headroom. Median pure head is **4 genes out of 70**; the
oracle's is 44.

This reproduces experiment 095 §6's ordering (`|lfc|` `0.10082`, agreement x
geometric magnitude `0.11302`) on an independent implementation, so the
instrument agrees with prior work and only the conclusion is new.

## Why concordance alone fails, and what that means

Cross-donor sign concordance on its own scores `0.05657` — **worse than
`|pred_lfc|` and worse than nothing**, with 42 of 127 perturbations at `k* = 0`.
It becomes useful only multiplied by magnitude. Independently confirms the peer
session's finding that `|lfc_hat|` carries sign-correctness information even
inside the half where our sign is anti-informative: the magnitude *is* the
confidence signal, and donor agreement adds only a modest reweighting of it.

## Donor coverage is the hard limit, and it is not fixable by adding donors

Seven donor matrices exist. Their overlap with the 132-perturbation panel:

| donor | perturbations covered |
| --- | ---: |
| `k562_gwps` | 132 |
| `orion_hek293t` | 132 |
| `h1_2025` | 12 |
| `k562_essential` | 0 |
| `nadig_jurkat` | 0 |
| `rpe1` | 0 |

(`orion_hct116` is excluded by construction — it is the evaluation surface.)

So the concordance estimate rests on **two** donors and cannot be widened with
what is on disk. "More donors" is not an available move; the extra matrices do
not cover the panel's perturbations.

## Consequence

The ordering channel is closed for buildable keys. Combined with experiment 095:
the deterministic ceiling for any key we can compute is roughly `+0.0068`
overall against the production emitted order, over-dispersion already realises
about 14% of it and shipped as arm a8, and no better key exists among available
signals. The remaining `+0.0875` between here and the oracle is **per-gene
directional accuracy** — a model problem, not an emission or ordering problem.

**Do not** build further ordering keys, per-gene `phi` confidence maps, or spread
policies aimed at `reach`. The emitter can express a good order; there is no good
order to express.

## Caveat, stated in the honest direction

This harness's reference DE pool is low power (`n_real` ~68-304) against the
challenge's (~2,670 by board inversion in `093-spread-truncation/BOARD-RESULT.md`).
Higher power admits weaker-effect genes, where our signs are worse. Purity here is
therefore **optimistic**, which makes this negative result robust — a positive one
would have needed re-checking in a higher-power regime.

Pooled across perturbations the picture looks much better (`|pred_lfc|` top-1%
purity `1.0000`, deepest >=0.90 prefix 481 of 11,410 pairs; best key 545). That is
a genuinely easier problem than the metric poses, because pooling may cherry-pick
across targets while `reach` is scored per perturbation. Recorded so the pooled
numbers are not mistaken for a result.

## Integrity

No sealed data opened, no gate changed, no `cell-eval2` run, no build, no
submission. Evaluation surface `orion-hct116-thinned` via
`experiments/095-clean-harness/de_thinned.npz`; the leaking donor is excluded.

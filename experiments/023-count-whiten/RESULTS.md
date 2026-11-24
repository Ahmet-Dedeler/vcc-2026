# 023 — Count-level: what whitening costs on the four DE metrics, and the CPM gate

**Date:** 2026-08-21 · real `cell-eval2` on `h1-fast` (40 perturbations, 2025 H1 hESC, its own
measured anchors). `pds` saturates on a 40-perturbation panel and is not read here; `jac` on
this harness has 6.7% seed noise against the competition's 0.8–2.2% and lied once, so it is
read as directional only. This run exists for **fid, reach, nmae**.

## Whitening: the DE cost is small, and the proxy was wrong about `reach`

At alpha = 1.08, CPM gate 5:

| | fid | reach | nmae | jac | mse |
|---|---|---|---|---|---|
| gamma = 0 | 0.5178 | 0.1538 | 0.9397 | 0.0728 | 1.385 |
| gamma = 0.65 | 0.5121 | **0.1532** | 0.9531 | 0.0747 | 1.467 |

Converted to the official scale: `fid` −0.019, `nmae` −0.022, `jac` +0.005 — about **−0.006
overall**, against a `pds` gain worth +0.017 to +0.023. Net clearly positive.

**`reach` is flat (0.1538 → 0.1532).** Experiment 018's proxy predicted a 39% collapse. It was
wrong, and the reason is instructive: the proxy ordered the reference's DE genes by
`|predicted effect|` and used "top-K by |true effect|" as the reference's set, but the real
metric orders by *adjusted p-value from a Wilcoxon on emitted counts*, which under a
minimum-variance emitter is `|effect|` measured **in units of the control's own cell-to-cell
spread**. That per-gene rescaling is exactly the thing whitening perturbs, and it happens to
absorb most of the damage. Proxy `reach` should not be trusted even directionally.

## The CPM gate: 5 beats 20, consistently, and it is not small

| alpha | gamma | gate 20 | gate 5 | delta avg |
|---|---|---|---|---|
| 0.7 | 0 | 0.2131 | 0.2184 | +0.0053 |
| 1.08 | 0 | 0.2177 | 0.2245 | +0.0068 |
| 1.6 | 0 | 0.2172 | 0.2262 | +0.0090 |
| 0.7 | 0.65 | 0.2176 | 0.2306 | +0.0130 |
| **1.08** | **0.65** | 0.2247 | **0.2362** | **+0.0115** |
| 1.6 | 0.65 | 0.2005 | 0.2095 | +0.0090 |

Six for six. At the shipped operating point the gain is carried by `reach`
(0.1416 → 0.1532) and `nmae` (0.9635 → 0.9531), worth roughly **+0.006 overall** on the
official scale.

The mechanism is the one the flag's own docstring gave and which v1/v4b ignored: the scorer
tests every gene above **5** CPM in the reference control. Gating at 20 makes us predict *no
change* on ~1,900 genes the scorer does test — genes that then enter `reach`'s pool with an
arbitrary sign and enter `jac`'s union for free.

## alpha: 1.08 is the optimum and it is not flat

At gamma = 0.65, gate 5: avg 0.2306 (alpha 0.7), **0.2362 (1.08)**, 0.2095 (1.6). `nmae` is
the metric that decides it — 0.9623 / 0.9531 / **1.0580**. Pushing alpha to 1.6 buys `reach`
(0.1532 → 0.1756) and destroys `nmae`, exactly the landmine logged on 2026-08-21. The alpha
already chosen for v5 on magnitude-parity grounds lands on the measured optimum.

# 167 — a target's own expression does NOT predict whether its effect is real

**Date:** 2026-09-02 · three surfaces with truth (both thinned Orion lines + h1-test),
deterministic (seed `20260902`), no emitter, no `cell-eval2`.

## Verdict

**NO-GO, and the null is more useful than the arm would have been.** A panel target's own
expression in the recipient context carries **no** information about whether its perturbation
response is real: `Spearman(own CPM, split-half reliability)` is **+0.122 / +0.046 / +0.057**
across the three surfaces, and targets below 5 CPM are as reliable as targets above
(`+0.142` vs `+0.159`, `+0.266` vs `+0.310`, and on h1-test the low ones are *higher*,
`+0.499` vs `+0.308`).

## 1. The observation that motivated it

Production builds ONE context-independent effect matrix and applies it to A/B/C. But a panel
target's own expression swings hard between contexts (measured on the real challenge controls):

| | median | q90 | frac > 2x | frac > 4x |
| --- | ---: | ---: | ---: | ---: |
| max/min across A/B/C, log2 | **1.411** | 2.871 | **0.690** | **0.300** |

and the tail is extreme — `VIM` is 1626 / 6132 / **2.02** CPM in A/B/C, `S100A11` is 0.12 / 525 /
2353, `KIF21B` 107 / 0.22 / 2.00, `PSMB9` 201 / 0.80 / 224. Below 10 CPM: 8 targets in A, **31**
in B, **27** in C.

CRISPRi represses transcription, so a gene that is already silent should not be repressible and
its downstream effect should be ~0 — yet we predict a full-magnitude effect for it. That is a
context-specific error made by construction, and being a mechanism rather than a fit it would
have transferred to the unseen D/E/F panel.

## 2. Why the obvious test would have lied, and what was measured instead

Correlating expression with `||true effect||` is confounded: at ~150 cells per target a NULL
perturbation still has a large measured effect norm made entirely of sampling noise, so a silent
target and a noisy one are indistinguishable. Indeed `Spearman(own CPM, ||effect||)` is
`+0.053 / +0.041 / +0.070` — uninformative, and it would have been read as "no relationship"
for the wrong reason.

The endpoint here is **split-half reliability**: each target's cells split into two disjoint
halves, each pseudobulked against a **disjoint half of the control** (sharing one control between
halves correlates their noise), effects taken in the scorer's `log1p(TS·p)` space with all
panel-target coordinates dropped, and the cosine between the halves. That is the reproducible
fraction of the effect and it goes to zero for a target with no real response regardless of how
large its noisy norm is.

## 3. Result

| own-CPM quintile (HEK293T) | n | median CPM | mean reliability |
| --- | ---: | ---: | ---: |
| 0.3 – 19.9 | 27 | 10.5 | 0.3017 |
| 19.9 – 34.9 | 26 | 26.5 | 0.2997 |
| 34.9 – 63.0 | 26 | 46.8 | 0.3149 |
| 63.0 – 105.2 | 26 | 79.5 | 0.3156 |
| 105.2 – 734.1 | 27 | 177.4 | 0.3086 |

Dead flat across a 17x expression range. HCT116 is 0.133 / 0.169 / 0.153 / 0.148 / 0.184 over
8.3 → 141 CPM; h1-test is non-monotone. **Nine HCT116 targets sit below 5 CPM (minimum 0.02) and
carry normal reliability**, so even near-silent targets have real, reproducible responses.

## 4. Why — and this is the part worth keeping

The on-target knockdown is **roughly constant in the target's expression**: median on-target
log2FC per quintile is −1.71 / −1.64 / −1.55 / −2.01 / −1.93 on HCT116 and −0.94 / −0.84 / −0.78
/ −0.84 / −1.10 on HEK293T, with `Spearman(own CPM, on-target log2FC)` = −0.164 / −0.104 / +0.119
— no consistent sign across surfaces.

That is the mechanism. CRISPRi acts at the promoter, so it delivers a roughly **fixed relative**
repression whatever the baseline abundance; and the size of the downstream response is set by the
target gene's **function**, not by how much of its transcript is present. Absolute expression
therefore has no bearing on whether the perturbation does something.

**This is a positive result for the current architecture.** It says the context-independent
effect matrix is *justified in its magnitude dimension*, and it independently explains
experiment 158's otherwise surprising finding that an oracle per-target **magnitude** — the true
`||delta_p||` handed over for free — is worth `+0.000011` on `mse`. There is no per-target
magnitude signal to exploit because relative knockdown is expression-invariant.

## 5. Scope and what is not claimed

- The surfaces' targets are the Orion 132 and h1 85, not the challenge 300. Their expression
  distributions are comparable (harness median 40–111 CPM, challenge 43–58) and the harness does
  reach down to 0.02 CPM, but the very extreme challenge cases (a target at 0.12 CPM in one
  context and 2353 in another) are represented by only 1–2 rows per context. **For those one or
  two rows the question is not settled by this experiment**, and it is worth at most a rounding
  error on a 300-row panel.
- Reliability is a property of the *reference*, not of our prediction. A high-reliability target
  can still be one we predict badly; this experiment says nothing about our accuracy.
- Reliability levels differ 2x between surfaces (0.157 HCT116, 0.308 HEK293T) at the same cell
  count, which is consistent with the 26.4% / 44.1% target-independent blob those two harnesses
  carry (`AGENTS.md` §6.4). **The levels are not comparable across surfaces; only the slopes
  within a surface are read here**, and all three slopes are flat.

## Files

`reliability.py` · `work/reliability.json`.

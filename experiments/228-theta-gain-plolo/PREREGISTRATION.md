# 228 — PREREGISTRATION (written before the first fit)

**Date 2026-09-04.** Committed before `gain228.py` produced a number.

## The object

224 §10's specified successor: take 224's **per-gene diagonal gain** — 189's cross-line
reliability prior `theta`, recovered by a regression estimator — and wire it in as a **per-gene
multiplier on the production effect matrix**, then price it on **all six members through PLOLO**
in the challenge frame.

## Wiring (fixed before any measurement)

```
pre  ──(× ALPHA)──> Operator(gamma=0.65, lam=0) ──> whitened
whitened ──(× m, then median-row-norm preserve over wcols)──> whitened'
whitened' + merged(CD4 @ f=0.34, H1union @ g=0.28)  ──> final   (panel coords zeroed)
final ──> deliver (a23 / minnorm render, challenge frames A/B/C) ──> scored
```

`m_j = (g_j / median g)^t`. `t = 0` is the shipped stack **exactly** and is the known-answer row.
The preserve rescale makes the arm amplitude-neutral at fixed `ALPHA` by construction; the
`VCC_ALPHA` interaction is then tested as a **separate** axis on top.

## Gain fit (fixed before any measurement)

- **Leave-one-donor-LINE-out**, inner, inside the fit pool. The PLOLO held line is excluded from
  the fit pool by cell-line GROUP, not by table name.
- **NON-PANEL targets only** (224 §7 — the panel surfaces read ~0 from power, not substance).
- Fitted on **whitened** inner pairs with the fold's own operator (224 §6's `g_w`).
- Shrinkage selected on a held inner LINE **and** a disjoint target block (224 §3.2's bug).

## Primary endpoint and gate

- **PRIMARY**: six-member PLOLO price, **marginal over the shipped stack** (`base222` =
  three donors + PM50 + CD4 `f = 0.34` + H1 union `g = 0.28`, board `0.2118310`), on the
  **three cancer folds** (`k562_gwps`, `orion_hct116`, `orion_hek293t`). `kolf2` reported
  **separately** and excluded from the headline.
- **GATE** (rule 13, shape-only class — profile-class arm, `reach` toll 0 per 169, to be
  *verified* with the per-gene panel-total digest, not assumed):

  ```
  price > +0.0019   AND   d_pds >= 0 on 3/3 cancer folds
  ```

- **Delta-pds-first filter** (216 INBOX-1): an arm that buys `cos_agg` while losing `pds` is
  ranked below one that does not, whatever its composite.
- **NO scalar PLOLO->board conversion.** The factor is not identified (1.384 / 1.000). Board
  points are quoted as a RANGE over `k_cos in [1.000, 1.384]` and always labelled as such.

## Falsifiers (mandatory, all four, through THIS instrument)

1. **gene-shuffle** — `g` permuted over genes. Expect NEGATIVE (224: −7.2%).
2. **moment-matching** — `sd(Y)/sd(P)` per gene, same estimator without target correspondence.
   Expect NEGATIVE (224: −3.5%).
3. **label-derangement** — `g` fitted with the target correspondence deranged. Expect strongly
   NEGATIVE (224: −25.7%).
4. **wrong-recipient** — 224's R4a shape: a 12-bin gain conditioned on the RECIPIENT's basal
   profile, right recipient vs wrong recipient. Expect the wrong recipient to be **no worse**
   than the right one and the whole rung to sit below the gate (224: +0.41%, wrong > right).
   The per-gene diagonal is recipient-independent by construction, so its own analogue is the
   **decoy/`pds` pair**: `pds` must RISE, and the off-diagonal (decoy) cosine must not rise as
   fast as the same-target cosine.

If a null does not reproduce its sign, the wiring is wrong and the arm is not read.

## Known-answer rows (the run aborts if any fails)

| # | gate | expected |
|---|---|---|
| KA-1 | `t = 0`, `cd4@0.34` | 223 `focus223` `BASE:cd4` per fold, to 1e-12 |
| KA-2 | `t = 0`, `cd4@0.34 + h1pair@0.28` | 223 `focus223` `BASE:base222` per fold, to 1e-12 |
| KA-3 | gain vs 189's `theta` | `r > 0` and of the order 224 §7 reports (+0.81/+0.90) |
| KA-4 | fit ∩ eval targets | exactly 0 |
| KA-5 | `t = 0` digest | per-gene panel-total digest identical to the shipped stack's |

## Leak audit

- Fit pool excludes the held line **by group**.
- `h1_2025` covers 85/85 of h1-test and is inside the shipped union — it is a production INPUT
  on every fold, never a truth, and is never in the fit pool.
- Fit targets are non-panel; eval targets are panel. Overlap asserted 0.
- `debench.Surface.pds()` is never called and sizes nothing (188 Phase 3).

## Stop rule

If the arm reads `<= 0` on `cos_agg` **or** `pds` on the three cancer folds at `t = 1`, report the
negative and do not sweep further doses looking for a positive. A well-supported negative is the
deliverable; 156 exists precisely because effect-level gains vanish at the read-back.

## Out of scope, by mandate

No `vcc`. No submission package. No production file touched. No `cell-eval2`. No model name.
Nothing written outside `experiments/228-theta-gain-plolo/`.

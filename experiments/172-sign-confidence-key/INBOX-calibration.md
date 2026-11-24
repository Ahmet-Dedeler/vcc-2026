# INBOX for 172 — calibration measured after your launch, 2026-09-02 08:05 UTC+3

Same surface and instrument as your brief (a15 h1-test DE table + `H` from the 095 artifact,
`Surface.reach(key=…)`, `|H|` = 0.23533 reproduced). Two things you need before designing the model.

## 1. The learnable oracle is 0.621, not 0.646 — and that is good news

| oracle | reach |
|---|---:|
| full: emitted sign correct, perfect order | 0.64553 |
| **H-oracle: DESIGN sign `sign(H)` correct, perfect order** | **0.62119** |
| H-oracle restricted to emission-consistent genes | 0.62166 |
| `|H|` (incumbent key) | 0.23533 |

Over the reference pool, **94.5% of emitted signs equal the design sign** — emission barely corrupts
the sign, so the "lucky coin flip" component is small and **96% of the ordering headroom is a
property of `H` itself, hence learnable in principle from features of the transfer.** Your target
label should be **`sign(H) == sign(real_lfc)`** (design-sign correctness), not the emitted sign — it is
what the features can explain, and it is available in every LODO fold without an emitter.

## 2. The signal is a thin excess over chance, and `|H|` is a weak but monotone predictor of it

Design-sign accuracy over the reference pool: **0.537** overall. By `|H|` decile (d1 = smallest):

| decile | `\|H\|` range | P(emission flips sign) | **P(design sign correct)** |
|---|---|---:|---:|
| d1 | 0 – 0.0062 | 0.363 | 0.510 |
| d2 | 0.0062 – 0.0132 | 0.128 | 0.514 |
| d3 | 0.0132 – 0.0208 | 0.037 | 0.516 |
| d4 | 0.0208 – 0.0290 | 0.011 | 0.524 |
| d5 | 0.0290 – 0.0384 | 0.007 | 0.526 |
| d6 | 0.0384 – 0.0490 | 0.001 | 0.528 |
| d7 | 0.0490 – 0.0621 | 0.001 | 0.533 |
| d8 | 0.0621 – 0.0794 | 0.001 | 0.545 |
| d9 | 0.0794 – 0.1073 | 0.000 | 0.554 |
| **d10** | **0.1073 – 1.33** | 0.000 | **0.624** |

Three consequences:

- **The bottom half of the pool is at chance (0.51–0.53).** No ordering of those genes can help;
  the whole game is identifying the correct-sign genes in the upper half and putting them first.
  A key that can *also* push confidently-wrong genes to the back is worth as much as one that pulls
  correct ones forward — consider modelling P(wrong) explicitly.
- **`|H|` alone reaches only 0.624 at its top decile.** The metric's purity bar is 0.9. Whatever
  features you add must lift the top 5% well above 0.72 (the current top-5% figure). If after honest
  LODO the top-5% accuracy does not exceed ~0.75, the information is not in these features and the
  experiment should say so — that is a real result.
- **Emission sign-flips are confined to d1–d2** (`|H|` < 0.013). Below that the emitted sign is noise
  regardless of the key. A schedule that delivers your key need not care about ordering inside d1–d2;
  it only needs to keep them *behind* everything else, which the metric's "significant calls first"
  rule does automatically if those genes are not called.

## 3. What the leader's number implies

`cqawesome` holds raw reach 0.2606 with `pds` 0.868 (ours 0.837). A modestly better direction
(perhaps 0.56–0.60 design-sign accuracy instead of our 0.537) delivered through a p-value order that
actually tracks it could plausibly produce 0.26 with **no** confidence model at all. So do not assume
the leader has solved your problem; the safer reading is that they have solved 171's. Your job is the
part beyond that.

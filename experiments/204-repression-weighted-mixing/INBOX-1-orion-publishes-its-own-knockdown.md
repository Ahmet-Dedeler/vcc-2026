# INBOX 1 for 204 — Orion publishes its own knockdown efficiency, and it confirms 194's measurement

Coordinator (lead), 2026-09-03 19:20 UTC+3. External datum, not a measurement of mine. Verify it
against the source before using it.

## The dataset card states the numbers directly

`hf datasets info Xaira-Therapeutics/X-Atlas-Orion` returns, in the description:

> *"eight million HCT116 and HEK293T cells, each deeply sequenced to a median of 16,000 unique
> molecular identifiers (UMIs) per cell. **The median on-target knockdown efficiency is 75.4% in
> HCT116 cells and 51.5% in HEK293T cells**"*

Against 194's independent measurement from the cells themselves (median self-gene residual
**0.300** and **0.510**, i.e. **70.0%** and **49.0%** knockdown):

| donor | Orion's published median knockdown | 194's measured | agreement |
|---|---:|---:|---|
| `orion_hct116` | **75.4%** | 70.0% | within 5.4 pp |
| `orion_hek293t` | **51.5%** | 49.0% | within 2.5 pp |

**Two things follow.**

1. **194's repression measurement is externally validated**, which is the premise your Phase 0 was
   asked to reproduce. Use the published figures as a **third route** alongside your self-gene
   residual and your shared-target regression slope. If your slope route disagrees with *both* the
   published number and 194's residual, the slope route is what is wrong.
2. **The HEK293T attenuation is a property of the screen, published by the people who ran it, not an
   artefact of our pipeline or of the long-gene blob.** That materially strengthens the mandate's
   premise: we are averaging a donor whose perturbations achieve ~51% knockdown with one achieving
   ~75% and one (K562) achieving ~85%, at equal weight. The blob control in your Phase 1 is still
   required — it protects the *estimate* — but the underlying 3.3× spread is no longer in question.

**One caution.** "Median on-target knockdown efficiency" as Orion defines it need not equal
`1 − median(self-gene residual)` as 194 defines it — the two agree to 2.5–5.4 pp, which is close
enough to corroborate and not close enough to substitute. **State which definition each of your
`a_i` uses**, and do not silently mix them.

Note also that Replogle's K562-gwps has no comparable published figure, so K562's `a_i` must come
from our own measurement on all three routes consistently — otherwise the comparison across donors
inherits a definitional difference, which is exactly the failure mode that would make a derived
weight look like a tuned one.

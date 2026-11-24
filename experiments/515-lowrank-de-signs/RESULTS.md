# 515 — Low-rank / Wiener denoising of a donor row lifts DE-pool sign accuracy by ~+0.01 on H1 (small; not queued)

2026-09-28 12:50Z, coordinator. Analysis only (`lowrank515.py`, 47 s, one SVD of 9,705 × 7,101 K562-GWPS rows).

## Question
The DE channel is sign-limited (513, 480). Our carriers use the WHITENED donor composition (γ 0.65), which is right for
pds (it up-weights target-specific directions) but is the opposite of denoising a weak gene's sign. Does projecting a
donor row onto K562-GWPS's response modules (basis learned on non-h1-test perturbations, no leakage) give better signs on
the recipient's DE pool?

## Result (H1 pool genes: reference q < 0.05, finite, non-zero, target excluded; K562 rows for 82/85 targets)
| source | all | pool < 20 | 20–100 | ≥ 100 |
|---|---:|---:|---:|---:|
| raw K562 row | 0.572 | 0.755 | 0.597 | 0.571 |
| top-10 projection | 0.574 | 0.446 | 0.516 | 0.575 |
| top-100 | 0.588 | 0.496 | 0.507 | 0.589 |
| top-300 | 0.589 | 0.568 | 0.549 | 0.589 |
| top-1000 | 0.586 | 0.597 | 0.572 | 0.586 |
| Wiener shrink s²/(s²+n²) | 0.581 | 0.734 | 0.613 | 0.580 |
| K562 programme (mean row) alone | 0.552 | 0.489 | 0.495 | 0.553 |

## Reading
- Low-rank projection helps the big pools (+0.018) and destroys the small ones: small pools are carried by specific genes
  (cis neighbours, direct targets) that module projections smear. The Wiener shrink keeps them and gains +0.009 overall,
  +0.016 on the 20–100 band.
- The size of the effect (+0.01 of sign accuracy ≈ −0.006 … −0.008 raw nmae on H1 by 480's rate) is the class of the
  consensus rider (477/480), not a new lever. Not queued. If a DE-sign build is ever made for another reason, the Wiener
  form of the donor rows is the one to use for the sign source (never the hard projection).

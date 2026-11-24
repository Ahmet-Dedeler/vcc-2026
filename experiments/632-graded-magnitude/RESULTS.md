# 632 — Analytic H1 price: identity-boosted signs + L1-optimal graded magnitudes on UNCALLED UP genes, Δnmae −0.0092 raw (if delivered rank-neutrally)

2026-10-02 ~13:05Z, council session [d4cd5f]; council #1 (the rank-neutral nmae lever), analytic stage. `grade632.py` (renumbered from 640: 640–649 belong to [537757]),
`work/grade632.json`; 10 s, 3 GB.

## Method
- Base: 446bj2z's DE construction on H1 (435's coarse + zs stats, seeds 0/1/2), base nmae 0.9409.
- Uncalled genes only (called genes keep their delivered lfc; their sign and size are tied to the call design).
- Each gene's direction comes from a score S, with magnitude m* = F⁻¹((2a − 1)/(2a)), the expected-L1 optimum. a and F
  (the gate-conditional |l| distribution) are taken per decile of |S| and CROSS-FITTED on the other half of the targets
  (odd/even, both ways).
- Called genes, p-values and order are held fixed. This is what a rank-neutral delivery would allow, since the leak cap
  of 426/435 is a rank effect.

## Result (Δ raw nmae vs the delivered lfc; seeds 0/1/2)

| signs on uncalled genes | all directions | UP-direction genes only | DOWN-direction only |
|---|---|---|---|
| z(CLEAN4) + 1.0·z(identity) (639) | −0.0057 / −0.0056 / −0.0059 | **−0.0092 / −0.0091 / −0.0092** | +0.0035 / +0.0035 / +0.0033 |
| z(CLEAN4) + 0.5·z(identity) | −0.0049 / −0.0048 / −0.0051 | −0.0083 / −0.0082 / −0.0083 | +0.0034 |
| CLEAN4 alone | −0.0041 / −0.0039 / −0.0040 | −0.0061 / −0.0059 / −0.0060 | +0.0020 |
| delivered signs, graded magnitude only | +0.0032 / +0.0023 / +0.0034 | −0.0025 | +0.0057 |

## Reading
- **The prize is the sign source times UP magnitude.** Magnitude alone on today's signs HURTS (+0.003), as 423/480
  found. Better signs (CLEAN4 + identity) raise the L1-optimal magnitude above the MDE on confident genes, which
  today's delivery cannot reach without leaking calls. UP genes can be raised rank-neutrally (inflate cells above the
  control max). DOWN genes already carry the zs / MDE treatment, and re-grading them hurts, so leave them.
- The H1 → board magnitude multiplier is 1.2–1.76 (611 / 433a), so −0.0092 raw ≈ −0.011 … −0.016 raw on the board ≈
  **+0.003 … +0.0045 overall**, with every other member pinned if the delivery is truly rank-neutral.
- Independent agreement: [21dd4c]'s 610 dosecurve (HYBRID key/CLEAN4 signs on uncalled genes) prices −0.005 … −0.011.
- Next: (1) fold the per-target key (21dd4c's HYBRID) into the CLEAN4 + identity score; (2) the real-emitter
  rank-neutral UP pass, under the coordinator's 4 conditions (KA, six members, 3 seeds, H1-full + h1r7; closures
  426/435/426b/434/480 addressed; U/calls/pds/mse pinned; a 2-point dose prereg).

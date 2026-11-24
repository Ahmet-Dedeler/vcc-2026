# 712: is A/B/C's light key tail a property of Arc's panel selection? (preregistration)

Written 2026-10-16 ~07:58Z by the science session [d80a10], before any number below exists.

## Why
The fold cap `hi` (P 96.647 of the LODO key / median over carried targets) is 5.2–7.0 on A/B/C's folds. It is 7.8–19.9 on the
lab's L panels, which were matched to A/B/C on coverage, strength mix and A/B/C basal expression, and 8.8–17.2 on our
real / efp / efr draws (708 RESULTS). The lab calls this an open question with a direct October consequence. If Arc's
selection makes light tails, D/E/F (Arc-selected) is likely light too. Then 708's D/E/F prices read like A/B/C's, and the
same-panel ρ (key 1.25) applies. If A/B/C is just unusual, D/E/F's tail is open, and only the day-of `hi` (708 prints it)
tells. A second Arc-selected panel exists on disk: the 2025 VCC H1 panel (train + validation + test, ~300 targets).

## Instrument (light; no operator fits)
- Proxy for the key tail: each target's donor strength in a fold line h (K562-GWPS, Orion HCT116, Orion HEK293T). Strength
  = the L2 norm of the h bank's delta row over the bank's genes, excluding the panel's own genes.
- Tail ratio per panel × fold = P 96.647(strength) / median(strength) over the covered targets, the same order statistic as
  `hi`. Also reported: P 90 / median, and the covered count.
- Panels:
  - Arc: A/B/C 2026 (`data/challenge/pert_counts.csv`) and H1 2025 (`data/external/vcc2025/*/pert_counts_*.csv`, the union).
  - Ours: 199's real (`def-bundle-r`), 238's efp / efr.
  - The lab's: L1–L6 (`experiments/801-lab-folds/work/L<k>/bundle/pert_counts.csv`). L7 and L8 are SEALED and are never read.

## Validity gate (before any reading)
The proxy must track the real `hi` where both are known: the 9 real / efp / efr folds (522's `hi546`), with Spearman
≥ 0.6. The proxy must also put A/B/C's three folds below every one of those 9 folds' proxy values on the same line. If
either fails, the result is "no read". The proxy then does not see the key's tail, and the question goes back to the lab
(fold operators on the H1 2025 panel).

## Readings (only if valid)
- **ARC-LIGHT:** H1 2025's tail ratio sits at or below A/B/C's on ≥ 2 of 3 lines and below every draw of ours or the
  lab's on those lines. Arc's selection makes light tails, so expect D/E/F light.
- **A/B/C-UNUSUAL:** H1 2025's tail ratio sits inside the range of our and the lab's draws on ≥ 2 of 3 lines. D/E/F's
  tail is open, and 708's day-of `hi` decides which ρ applies.
- Anything else: MIXED, reported as is.

## Class and use
Diagnostic only. No recipe switch reads it. It informs how 708's D/E/F read is interpreted (the information column) and the
writeup. Tries: one proxy and one statistic, fixed here.

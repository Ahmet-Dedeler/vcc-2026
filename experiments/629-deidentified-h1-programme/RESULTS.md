# 629 — Arc's H1 programme does not transfer, and identity erosion is NOT why (negative; closed)

> **AUDIT 616 (7dac7a98, [21dd4c], 10-02): the thinned Orion gate tables' DE truth was ~98 % selection artefact.** Here the HCT116 / HEK293T columns are VOID. K562 (H1 programme −0.064 … −0.023 vs CLEAN4) and the Broad-16 (0/12 lines) stand, so the closure stands.

2026-10-02 ~16:20Z, council session [7772dc] (territory 620–629). Analysis only (`deid629.py`, `work/deid629.json`), < 2 GB,
20 s.

## Question
Arc's own 2025 H1 programme (458d's artifact: P_down per gene from Arc's H1 screen, the same lab and chemistry as A–F)
is +0.087 gate sign on H1 itself but far below CLEAN4 on other lines (613: −0.17 / −0.19 / −0.09). 458d hurt on the
board. If a line's programme were a shared core plus its own identity erosion (622 / 627 / 628), removing H1's ESC
identity term would recover the core. Arc's lab programme would then become a sign source for D/E/F.

## Method
- H1 identity: I_H1 = log1p CPM_H1 − mean log1p CPM over A/B/C (flex_refs604, Flex). zP = (P_H1 − ½)/sd.
- De-identified: zP − β·z(I_H1), with β the OLS slope over H1's genes.
- Controls:
  - an expression-detrended version, zP − β_E·z(log1p CPM_H1) (609's caution);
  - H1 restricted to low-|I| genes (|z(I_H1)| < 0.5).
- Targets:
  - 636's HCT116 / HEK293T / K562 gate surfaces (I_s as in 639);
  - the Broad-16's 12 working lines (same-lab I_s; the reserved 4 never read).
- Statistic: paired gate-sign Δ vs CLEAN4 with a target bootstrap.

## Result
corr(zP, z(I_H1)) = +0.16 (β +0.10): H1's programme carries only a modest identity term.

| source, Δ vs CLEAN4 | HCT116 | HEK293T | K562 | Broad-16 mean (lines > 0) |
|---|---:|---:|---:|---:|
| H1 raw | −0.149 | −0.154 | −0.064 | −0.109 (0/12) |
| H1 de-identified | −0.144 | −0.152 | −0.055 | −0.109 (0/12) |
| H1 expression-detrended | −0.159 | −0.174 | −0.064 | −0.108 (0/12) |
| H1, low-\|I\| genes only | −0.109 | −0.128 | −0.023 | −0.111 (0/12) |
| CLEAN4 + 0.5·H1 de-identified | +0.002 | −0.010 | −0.006 | −0.007 (2/12) |
| 639 rider (CLEAN4 + 0.5·I_s) | +0.037 | +0.017 | −0.005 | +0.009 (10/12) |
| rider + 0.5·H1 de-identified, vs the rider | −0.012 | −0.016 | +0.003 | −0.008 (2/12) |

**Where the H1 programme fails.** By I_H1 tercile, the raw programme is WORST on the genes H1 expresses BELOW the cancer
lines: low-I tercile −0.166 / −0.200 / −0.101 on HCT116 / HEK293T / K562, against mid −0.107 / −0.140 / −0.023 and high
−0.150 / −0.116 / −0.042. On other lines its accuracy is 0.51–0.54, close to a coin flip, against CLEAN4's 0.60–0.69.
It is uninformative, not reversed.

## Reading
1. **The H1 programme's failure to transfer is not identity erosion.** Removing H1's identity term, its expression trend,
   or its identity-extreme genes leaves it 0.06–0.17 below CLEAN4 everywhere: 0/12 Broad lines, all three folds.
2. It is uninformative on other lines (near 0.5), and worst on genes that are low in H1. Those genes' signs in H1's own
   screen are plausibly set by H1-specific mass changes (610's CPM coupling), which do not transfer.
3. **Closed:** Arc's 2025 H1 programme in any form tested (raw, de-identified, detrended, restricted, or blended at 0.5
   into CLEAN4 or the identity rider) is not a sign source for D/E/F. Adding it to the identity rider costs 0.008–0.016.
   This confirms 458d / 613. The identity principle explains the TARGET line's own deviations, which the rider uses; it
   does not explain why another line's whole programme fails to carry over.

# 628 — CLEAN4's banks carry their own lines' identity erosion; removing it does NOT transfer better (negative)

> **AUDIT 616 (7dac7a98, [21dd4c], 10-02): the thinned Orion gate tables' DE truth was ~98 % selection artefact.** Here the hct116 / hek293t rows are VOID (de-identification Δ there was −0.009 / +0.006). The h1 (−0.002) and K562 rows stand, so the negative lever stands. The bank-erosion science reads the banks' own data and is unaffected.

2026-10-02 ~13:55Z, council session [7772dc] (territory 620–629). Analysis only (`deid628.py`, `work/deid628.json`), < 1 GB,
under a second.

## Question
Under perturbation, a line's identity genes (those expressed above other lines) tend to go down (639, 622, 627). CLEAN4
is the mean P_down of four CRISPRi banks: K562 GWPS, K562 essential, Nadig Jurkat and Nadig HepG2. If the principle holds
inside each bank, CLEAN4 carries the BANK lines' identity erosion, with K562 counted twice. For a new line that component
should be noise. Does removing it make CLEAN4 transfer better leave-one-line-out?

## Method
- Per bank b: I_b = b_ctrl_b − the mean of 3′ reference control profiles (log1p CP50k).
  - References: K562, Jurkat, HepG2, RPE1, HCT116, HEK293T, TeloHAEC and melanoma, one profile per line.
  - Excluded: b's own line, and on each evaluation surface s, s's own line and its same-study sibling. The surface's own
    expression therefore never enters its prior.
- β_b = OLS slope of (P_b − ½) on z(I_b). P_b' = P_b − β_b·z(I_b). CLEAN4' = the LOLO mean of P_b'.
- Scored on 636's gate surfaces with 636's paired target bootstrap. The 639 rider (z(C − ½) + 0.5·z(I_s)) is scored on
  both arms with the same I_s.

## Result 1 (science, positive): the identity principle holds INSIDE the banks
Slope of a bank's P_down on its own line's identity, against the H1 surface's reference set:

| bank | β | r | n genes |
|---|---:|---:|---:|
| K562 GWPS | +0.0049 | +0.159 | 6,346 |
| K562 essential | +0.0227 | +0.278 | 6,346 |
| Nadig Jurkat | +0.0121 | +0.135 | 6,346 |
| Nadig HepG2 | +0.0007 | +0.006 | 6,346 |

Other reference sets give the same picture (K562 GWPS r +0.20 on the Orion surfaces' set). In three of the four banks,
the genes a line expresses above other lines are called DOWN more often in that line's own screen. That includes the
genome-wide K562 GWPS screen, a different lab, year and design from every surface in 639/622/627. HepG2 shows no
effect.

## Result 2 (lever, negative): de-identified CLEAN4 does not transfer better

| surface | CLEAN4' − CLEAN4 | rider on CLEAN4' − rider on CLEAN4 | (639 check: rider − CLEAN4) |
|---|---|---|---|
| h1 (Flex) | −0.0020 [−0.0043, +0.0004] | −0.0008 [−0.0026, +0.0009] | +0.0140 |
| hct116 | **−0.0091 [−0.0126, −0.0058]** | **−0.0047 [−0.0088, −0.0011]** | +0.0387 |
| hek293t | **+0.0059 [+0.0043, +0.0076]** | +0.0020 [+0.0005, +0.0035] | +0.0174 |
| k562 (low power) | −0.0092 [−0.0193, +0.0006] | −0.0037 [−0.0106, +0.0031] | −0.0049 |

The effect is mixed, and the mean is slightly negative. **Closed: do not de-identify CLEAN4.**
Per [4d0207]'s 609 caution, read h1 first: on the clean surface the effect is null-to-negative (−0.0020 [−0.0043, +0.0004];
on the rider −0.0008). The Orion rows (hct116 −0.009, hek293t +0.006) are tilted by those folds' 72–74 % DOWN block gate
genes. They point in opposite directions, so they do not change the reading.

A likely mechanism: "identity" measured against a reference mean is partly SHARED across cancer lines. A gene K562
expresses above RPE1, TeloHAEC and the rest can also be high in HCT116, so its down-tendency transfers. The bank
component is not pure noise for a new line. What carries is the target line's OWN identity, which is the 639 rider.

## Class / use
- Science: one more independent dataset family for SCIENCE.md (sent to [4d0207]). The sign effect appears in three public
  CRISPRi banks' own screens, in Broad's Cas9 KO screen (622g), in Arc's H1 and in Orion's lines (639).
- No production change.

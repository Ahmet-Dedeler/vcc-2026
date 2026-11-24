# 627 — Identity across perturbation modality: Tahoe drugs show a continuous identity drift (partly ambient); Broad KO does not, but its gate-gene SIGN effect replicates

2026-10-02, council session [7772dc] (science item from [4d0207]). Analysis only, on 184's cached Tahoe tensor
(tahoebio/tahoe-de-rhaister: 1,137 drug-dose treatments × 50 cancer lines × 2,000 genes, DMSO-subtracted, library-
normalised, with per-plate DMSO basal profiles). `tahoe627.py`, `ambient627.py`; `work/*.json`. Peak RSS 1.7 GB.

## Design
- **Shared-noise trap removed.** delta subtracts the plate's own DMSO, so an identity axis from the same controls would
  carry the opposite-signed sampling noise. Every (drug, line)'s identity is therefore built only from plates that
  pair was NOT measured on (leave-plate-out): I_L = log1p CPM_L − mean of the other 49 lines, on genes expressed in L.
- **The drug's own signature removed.** r = delta − the drug's mean across all lines, i.e. the line-specific deviation.

## Result A: the line-specific deviation follows −identity
- 56,777 (drug, line) pairs. Mean Spearman(r, −I) is **+0.147**, positive in 70.1 % of pairs. The wrong-line-identity
  null gives +0.013 (54.2 %). The top-100 |r| genes' sign agrees with −I at 0.584.
- **It grows with response strength.** By |delta| quintile: +0.085 / +0.092 / +0.124 / +0.204 / +0.231; top-100 sign
  0.541 → 0.627. Nulls stay at +0.009 … +0.019.
- **Dose:** 0.05 µM +0.112, 0.5 µM +0.232, 5 µM +0.097. The non-monotone shape is plausibly death or selection at
  the top dose.

## Result B: the generic response, a line's mean over a plate's drugs
Against −I from the other plates: +0.206 (77 % of 600 line × plate cells positive), against a null of +0.064.

## Confound: ambient RNA (Tahoe pools all 50 lines per well)
- Extra soup in treated cells pulls each line toward the pooled mean, which mimics erosion. Its fingerprint is that
  genes ABSENT in a line (CPM < 5) but expressed elsewhere (> 50 CPM, ≥ 10×) gain reads. Biology predicts ~0 there.
- `ambient627.py`, 49,955 pairs, median 30 such genes per line, by strength quintile:

  | quintile | absent-gene rise | identity-gene (top-decile I) drop | low-I decile |
  |---|---:|---:|---:|
  | Q1 | +0.0015 | −0.0029 | +0.0020 |
  | Q2 | +0.0019 | −0.0040 | +0.0028 |
  | Q3 | +0.0022 | −0.0081 | +0.0033 |
  | Q4 | +0.0036 | −0.0171 | +0.0057 |
  | Q5 | +0.0064 | −0.0386 | +0.0133 |

- **An ambient component exists** and grows with cytotoxicity, but **it is too small to explain the drop.** In log1p
  units a soup increase Δa raises an absent gene by ≈ 0.5·Δa (soup ≈ 0.5–1 CP10k against ~0.05 in the line) and lowers
  a high-identity gene by ≲ Δa. So soup predicts |identity drop| ≲ 2 × absent rise. Observed: 4.5× overall, 6× in Q5.
  Roughly half to two-thirds of the erosion is beyond ambient. This is semi-quantitative: it rests on the soup model and
  on delta's units.

## Result C: Tahoe's statistic on CRISPR KO does NOT replicate (`broad627.py`, `work/broad627.json`, 2026-10-02 ~13:05Z)
This is the same test on the Broad-16 Cas9-KO screen (12 working lines; the 4 reserved lines are never parsed). Each line
is screened alone with its own OR-control guides, so there is no cross-line soup.
- **Shared-control trap.** The 10 OR control KOs are split by name, 5 + 5. The response is computed against half A
  (lfc = log2((pb + 1)/(ctrl_A + 1)) on CPM pseudobulks). I_L is built from half B against the other 11 lines' half B.
- r = lfc − the target's mean over the 12 lines, over 87 targets common to all lines.
- Mean Spearman(r, −I) over 1,044 (line, target) pairs is **−0.019**, positive in 33.7 %. The wrong-line null gives
  −0.003 (46.8 %).
- By strength quintile: −0.024, −0.006, −0.014, −0.023, −0.028 (nulls −0.007 … +0.008). It does not grow with strength.

**The continuous line-specific deviation does not follow −identity in KO.** If anything it leans slightly the other way.

## Result D: the gate-SIGN effect in KO is real and survives control halves (`halves622g.py`, `work/halves622g.json`)
622g's sign result was re-read with I from control half B and the sign truth recomputed against half A, on the authors'
fdr < 0.05 gates:

| | identity alone (mean over 12 lines) | z(C4) + 0.5 z(I) − CLEAN4 |
|---|---:|---:|
| halves | 0.5436 | +0.0093 |
| full controls (622g's construction) | 0.5440 | +0.0091 |

- 622g is not shared-control noise.
- 622m (`../622-identity-track/block622.py`) adds two facts:
  - the gain sits on gate genes OUTSIDE the top-2,200-expressed block;
  - most of it survives regressing expression level out of I (+0.0052 of +0.0081).

## Reading (revised 2026-10-02 ~13:35Z after C and D)
- **What generalises across modality is the gate-gene SIGN effect, not Tahoe's continuous statistic.**
  - Identity predicts the direction of called DE genes beyond CLEAN4:
    - CRISPRi in Arc H1, Orion HCT116 / HEK293T and K562 (636/639);
    - Cas9 KO on 12 Broad lines (+0.009, survives halves);
    - Tahoe drugs in 50 lines (top-100 |r| sign 0.584).
  - The **continuous line-specific deviation** follows −identity in Tahoe (+0.147, partly ambient), but **not in the
    Broad KO screen** (−0.019).
  - So "perturbations push lines toward a common state" is a Tahoe (drug, pooled-well) finding, and a soup-sensitive
    one. In KO the identity signal lives on the called genes' signs and is not a whole-transcriptome drift.
- **Manuscript claim** (for [4d0207]'s SCIENCE.md, replacing the earlier wording):
  - A line's control-derived identity predicts the sign of its called DE genes beyond public-bank priors.
    - It replicates across CRISPRi (4 surfaces) and Cas9 KO (12 lines, a different lab and chemistry).
    - It is not an expression-level or shared-control artifact.
  - In drug screens the whole-transcriptome deviation also tracks identity, but pooled designs inflate it with
    ambient RNA, and it does not replicate in a single-line KO screen.
  - Do not claim "identity erosion across modalities" as a continuous drift.
- **The effect strengthens with response strength in Tahoe only.** The Broad KO continuous statistic shows no trend.
- **For VCC:**
  - The identity source is transfer-class biology for SIGNS.
  - Its delivery route is unresolved: 622m puts the Broad-16 gain on uncalled genes, against 639d's Orion block
    localisation.
  - The delivered value is open (622 CORRECTION; 625d / the board).

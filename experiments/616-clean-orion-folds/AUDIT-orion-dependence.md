# 616 AUDIT — every reading and decision that used the thinned Orion DE truths

[21dd4c], 2026-10-02 ~19:35Z. Requested by the coordinator (19:06Z): list every instrument reading and decision that
used Orion fold rows and, for each, what remains (H1 / Broad-16 / K562 / nothing). Nothing has been recomputed: S2 is
HELD behind `work/RELEASE-S2`.

## The finding this audit rests on (616 S1; S0 known answers all PASS, `work/ka616.json`)
480's gate rules, with the same 132 targets and 150 cells per target. Only the truth's cell selection differs: the old
greedy selection vs random cells over all batches with 615's P′ control pool.

| HEK293T gate table | thinned (`480/work/gate-orion_hek293t.npz`) | **clean** (`work/gate-orion_hek293t_clean.npz`) |
|---|---|---|
| downstream DE calls | 37,748 | **774 (2 %)** |
| median calls per target (q25 / q75) | 141 (42 / 390) | **1 (0 / 3)** |
| targets with ≥ 10 calls (the ones the instruments score) | 125 / 132 | **13 / 132** |
| genes called in ≥ 50 % of targets (artefact signature) | 44 | **0** |
| DOWN share of calls | 0.562 | 0.565 |

- **HCT116 (landed 19:25Z):** same 132 targets, 150 cells.
  - Downstream calls 11,309 → **712 (6 %)**.
  - Median per target 68 → **2** (q75 113 → 4).
  - Targets ≥ 10 calls: 127 → **18 / 132**.
  - Genes in ≥ 50 % of targets: 12 → **0**.
  - DOWN share 0.598 → 0.643.
  - Library median 12.1k in both.
  - The same conclusion holds: **no usable DE truth on either Orion fold.**
- 615's deeper HEK truth (220 cells) agrees: median 1, 18 / 100 targets ≥ 10.
- **Conclusion for HEK293T (and HCT116, below): the Orion folds' DE truth carries no usable signal.** About 98 % of what every instrument
  graded against on that fold was the builder's selection artefact. A clean truth leaves ~13 scoreable targets, which is
  less power than K562's fold.

## NOT affected (checked)
- **The Orion DONOR tables** (`data/donors/orion_*`, `scripts/build_orion_donor.py`) pool every control and every cell
  over all batches. That is 615's configuration P, whose split-half blob share is 0.064 against h1-test's 0.049. So:
  - production's Orion donors are clean;
  - PLOLO's Orion folds (bulk truths from donor tables) are clean;
  - the bulk-channel October decisions resting on PLOLO are not touched: **q, lodo, lodo_hi_pct, sharpen_pt,
    cd4_stack**.
  - Small nit: the donor builder's `^neg` control regex would also count NEGR1-targeted cells as controls (~200 of ~219k
    controls), which is negligible.
- **h1-test, 480's K562 gates, the Broad-16 / Broad-12 lines, 614 / 615.**
- **de_mag** (H1 and board), and **K16 / K20** (611, H1-based).

## A. Instruments that read the thinned Orion DE truths (480's `gate-orion_*.npz` or the harness h5ads)
| # | instrument / reading | the Orion rows | status | what remains |
|---|---|---|---|---|
| 1 | **480 fold480**, the "cancer folds as a second instrument" for DE-magnitude constructions | HCT116 + HEK293T rows, including the "inversion via the same-study sibling" | **VOID** | H1 (calibrated); K562 (71 / 267 targets ≥ 10, low power) |
| 2 | **636 LOLO sign instrument** (gate-gene sign accuracy per source vs CLEAN4) | hct116 / hek293t surfaces (CLEAN4 0.69 / 0.66; "Nadig Jurkat best on Orion") | **VOID** | h1, k562; Broad-16 (independent-lab sign accuracy) |
| 3 | **639 identity axis LOLO** ("transfers": h1 +0.016, hct116 +0.015, hek293t +0.024) | 2 of its 3 positive lines | **VOID** (639d's block localisation was already retracted by 609) | h1 +0.016; Broad-16 +0.009; K562 null |
| 4 | **613 block-prior instrument**: R1 458d retrodiction (rule: 2 of 3 cancer folds), R2 392 | hct116 (−0.170 / −0.026), hek293t (−0.194 / −0.038) | **VOID for 2 of 3 folds.** The 2-of-3 rule can no longer be met, so **613 is no longer a validated block instrument by its own rule** | K562 alone still agrees with both board points: 458d gate −0.085, null −0.008; 392 null +0.004. One low-power fold. H1 is not part of the rule |
| 5 | 613 (B) identity null-block ("below CLEAN4 on 2 of 3 cancer folds", which restricted ID05 to TRANSFERS) | hct116 −0.014, hek293t −0.008 | **VOID** | K562 +0.0001 (not below). The TRANSFERS-only restriction now rests on 609 (clean surfaces: the block share is small) |
| 6 | 613c bank composition ("Nadig Jurkat alone beats CLEAN4 on both Orion lines") | both rows | **VOID**; no decision used it ("not actionable") | Broad-12 (613f / 613g decided "CLEAN4 stays"), H1, K562 |
| 7 | **613d hook add-rule instrument**: per-fold c4h − base, and **the registered 496a prediction** (cancer-fold mean at f 0.35 +0.00017 = HCT116 +0.00001, HEK −0.00076, K562 +0.00127) | 2 of the 3 folds in the registered mean | **VOID for HCT116 / HEK.** The registered prediction is now effectively K562 (+0.0013) averaged with two void rows | H1 real emitter: the instrument's 2 out-of-sample PASSES are on clean H1 (idkg3, idall3); K562 row |
| 8 | 613d3 depth rows, 613d4 KGID decomposition, **613d5 / 613d6 ungated forms** (idall3 − c4h: HCT116 −0.0016 / −0.0011 / −0.0020, HEK −0.0039 / −0.0026 / −0.0011) | Orion rows | **VOID** | H1 (ungated v3 −0.0016 / −0.0023 / −0.0023, out-of-sample PASS on the real emitter); K562 opposes (+0.0005 / +0.0016 / +0.0024, weakened CLEAN4); Broad-12 (612 addendum: neutral on KO lines, 10–15 % K562-like tail) |
| 9 | 613e / 613e2: **630's premise** (uncalled UP-signed gate mass share) | Orion rows | **VOID** | H1 (0.25), K562 (0.134) |
| 10 | 612 ds612b / ds612c fold model; **KGID's regret table** (max regret 0.0032, HEK vs CLEAN4); "the plain hook pays where CLEAN4 fits (Orion +0.015)" | Orion rows | **VOID.** The fold model was already void for hooks (612 A1) | H1 real emitter; K562 (−0.010, the stress case) |
| 11 | **516 DE sign sources**: "≈ 0 on every cancer fold" (HCT116 +0.0033, HEK −0.0002), which is part of **programme_source = off** | 2 of 3 cancer folds | **VOID** | K562 only (+ H1, the programme's own line) |
| 12 | 521 strength → n_real predictors (r on Orion folds) | HCT116 / HEK rows | **VOID** (the Orion n_real WAS the artefact) | H1; K562 (r +0.43) |
| 13 | 621 supervised DE (P(down \| gate) shift; membership AUROC 0.79 / 0.74 on Orion) | Orion rows | **VOID** | h1 0.68, k562 0.74 |
| 14 | 622 identity track (Orion rows +0.003 … +0.042; "public 3′ refs harmful on hek293t") | Orion rows | **VOID** | h1; Broad lines where used |
| 15 | 625 module readout (Orion \|x_M\| already flagged inflated) | Orion rows | **VOID** | h1 |
| 16 | 628 de-identified CLEAN4 (Orion rows already flagged as tilted by the DOWN block) | Orion rows | **VOID** | K562, h1, Broad |
| 17 | 629 de-identified H1 programme | HCT116 / HEK rows | **VOID** | K562; Broad-16 (0 / 12 lines) — the NEGATIVE stands |
| 18 | 631 nmae headroom ("every cancer surface wants ~2× more called-DOWN magnitude"; Orion already marked "not clean for this") | Orion rows | **VOID** | K562 (+ H1) |
| 19 | 633 dose folds (Orion d20 rows K 1.5–2, descriptive) | Orion rows | **VOID**; not binding | K562 d20 (binding) |
| 20 | 638 responder programme (Orion Δ gate sign −0.099 / −0.205) | Orion rows | **VOID** | h1 −0.013, k562 +0.001: the NEGATIVE stands |
| 21 | 609 identity split (diagnosed the Orion block artefact) | its Orion rows ARE the artefact it named | consistent; its conclusions sit on H1 / Broad-16 / K562 | — |
| 22 | 610 decile610 (CLEAN4 null-gene sign precision by expression rank: HCT116 0.543 → 0.501, HEK 0.550 → 0.494) | Orion rows | **VOID** | h1 0.518 → 0.511, K562 0.514 → 0.500: same conclusion |
| 23 | 513 pure-head subset (sign accuracy on H1's pools; the Orion lines appear only as DONORS) | donors only | not affected (donor tables are clean) | — |

Pre-council (031–194) experiments that read `orion-*-thinned.h5ad`, `orion-hct116.h5ad` or `-depthmatched` directly:
031, 033, 035, 043, 049, 054, 067, 092, 093, 095–100, 102, 104, 106, 117, 149–153, 155–157, 159, 160, 167, 192, 194,
273.
- FINDINGS H.4 (09-02) already restricted them to "relative direction only".
- **New:** their DE-channel truths (Wilcoxon calls: fid / jac / reach / nmae) were about 98 % artefact, so DE-channel
  conclusions drawn on those surfaces should be treated as void, including relative ones.
- Bulk mean-removed comparisons (H.4: 156, 099) are less exposed, because the artefact is mostly a shared offset.
- Not re-audited line by line here.

## B. Decisions and open calls touched
| decision / call | status | what used Orion | what remains | note for the re-ruling |
|---|---|---|---|---|
| **mde_sign** (PENDING; 496a board read 10-05) | open | 613d's registered prediction (+0.00017) averages two void rows with K562 | H1 real emitter (c4h − base −0.0015 / −0.0022 / −0.0021, clean); K562 +0.0013 | The 10-05 read no longer tests "613d on held-out lines". It can only test the K562-dominated composite. The registered prediction is untouched; a clean-fold version can't be computed while S2 is held, and with ~13 scoreable targets it would be power-limited |
| **idall3** (prereg FINAL, main queue slot 4) | open | 613d5 / 613d6 Orion rows (favourable) | H1 out-of-sample PASS (real emitter, −0.0019 … −0.0022 vs c4h); Broad-12 neutral with a 10–15 % K562-like tail; K562 against | Its held-out case is H1 + Broad-12. Both favourable Orion rows are gone, and the K562 stress row opposes |
| KGID | retired 15:25Z | 612 regret table, 613d4 | — | nothing to re-rule |
| **programme_source = off** (DECIDED) | decided | 516's cancer folds | K562 only | "off" stays the conservative ruling; the held-out line now reads "K562 only" |
| **604 v3 construction** ("identity-only genes get NO prior: their signs were right 0.27–0.34 on HEK293T") | shipped in v3 | HEK fold | none (that rule was HEK-only) | The exclusion is conservative (fewer genes hooked), but its stated evidence is void. Re-test the identity-only extension on H1 / Broad-16 before calling it closed |
| **639 "identity axis transfers"** (POSITIVE) | standing | hct116 / hek293t rows | h1 +0.016, Broad-16 +0.009 | still positive, on 2 independent surfaces instead of 4 |
| **613 block instrument = "validated (2 board points)"** | standing | 2 of 3 folds | K562 alone agrees with both board points | Downgrade to "consistent on K562 only; not validated by its own 2-of-3 rule". Anything priced on it (ID05's block arm, CLEAN4-vs-source block readings) needs K562 + Broad-16 rows |
| 630 UP pass (final 7e6b3e52, freeze 10-05 18:00Z) | standing | 613e / 613e2 Orion mass shares; 633's Orion d20 rows | H1 (d4cd5f's 3 seeds), K562 (633 binding row) | its binding evidence was already H1 + K562 |
| ID05 block prior → TRANSFERS-only | decided | 613(B), 639d | 609 (clean: the block share is small) | unchanged |
| CLEAN4 composition (stays) | decided | 613c, 628 | Broad-12 (613f / 613g) | unchanged |
| 516 / 621 / 622 / 625 / 629 / 638 NEGATIVES | closed | Orion rows | H1 / K562 / Broad rows | all still NEGATIVE without Orion |

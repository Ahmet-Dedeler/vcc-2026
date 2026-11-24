# 651 — a target's RNA-content shift transfers across lines, but it does not set the filler's signs (negative; closed)

2026-10-02 16:43–16:54Z (commit times), council session [7772dc]. Analysis only (`lib651.py`, `run651.py`; `work/lib651.out`,
`work/run651.out`, `work/targets651.csv`, `work/beta651.npz`). Peak RSS 4.7 GB, 11 s. The PLAN and its bar were committed
before steps 2–3 ran (47b018f1).

## Question
The pack's fid edge is filler sign precision: p_null ≈ 0.52–0.525 against our ≈ 0.508 on A/B/C, and still unexplained.
- 610 cpmcouple closed the uniform per-target route: DE-mass loss lifts every other gene's CPM together.
- 651 tested the gene-specific route. Gene g's per-cell CPM share moves with library size along a slope β_g, measured on the
  line's own controls. If target t's cells shift in RNA content by Δ_t, the truth's null-gene lfc carries β_g Δ_t, and
  sign(β_g Δ_t) varies across genes within a target.
- The decisive part was predicting Δ_t without the truth, from K562 GWPS's library shift for the same target.

## Results (h1-test truth, 85 targets; block null genes = the top 2,200 control-CPM gate genes not called by the truth)
**Step 1 (descriptive).**
- Perturbed cells do shift in library size. Δ_t has sd 0.050 across targets against a 400-cell se of 0.017, and 14 % of
  targets lie beyond 3 se.
- The extremes are the DE-mass targets: MED19 −0.29, HSP90B1, CSNK1E, EP300 ≈ −0.10, and ZBTB39 +0.17.
- The batch mix is balanced (TV 0.021).

**Step 2 (oracle Δ_t, descriptive).**
- β_g is gene-specific (relative slope sd 0.12, uncorrelated with expression at r 0.012).
- The truth's null genes do not move along it. The per-target slope of null lfc on β_g tracks Δ_t at only Spearman 0.13,
  with an OLS ratio of 0.016 where the model predicts 1.
- sign(β_g Δ_t) alone gives p_null 0.511, against CLEAN4's 0.516.
- Combined with CLEAN4, it gains +0.0027 [−0.0005, +0.0059] at λ 1 (+0.0032 at λ 0.5), under the plan's +0.003 oracle bar.

**Step 3 (decisive; K562 GWPS Δ̂_t = mean log UMI vs the same gem group's NTCs; 82 of 85 targets covered).**

| λ | p_null (CLEAN4 0.5169) | Δ vs CLEAN4 | 95 % CI |
|---|---:|---:|---|
| 0.5 | 0.5176 | +0.0007 | [−0.0015, +0.0029] |
| **1 (decides)** | 0.5167 | **−0.0002** | [−0.0031, +0.0025] |
| 2 | 0.5144 | −0.0025 | [−0.0063, +0.0013] |

**Verdict against the bar:** Spearman 0.39 ≥ 0.30 passes. Δ ≥ +0.005 fails. CI lower > 0 fails. **CLOSED.**

## Reading
1. **New transfer fact (science).** A target's effect on its cells' total RNA (library size) is shared across lines: K562
   GWPS predicts H1's per-target shift at Spearman 0.39 over 82 targets, across lab, chemistry and cell type.
2. **The shift does not reach the filler's signs.** Perturbed cells that lose library do not re-compose along the controls'
   cross-sectional size axis (ratio 0.016). The shift is mostly the DE genes' own lost mass, which 610 showed acts
   uniformly. Both library-coupling routes (uniform, 610; gene-specific, 651) are closed as filler sign sources.
3. The pack's p_null edge stays unexplained. Routes closed so far:
   - null skew from controls (610);
   - DE-mass coupling (610);
   - the identity filler (+0.001 … +0.003);
   - CLEAN4 bank subsets (613f/g);
   - expression deciles (610);
   - library-size composition (651);
   - per-target sizing (521).

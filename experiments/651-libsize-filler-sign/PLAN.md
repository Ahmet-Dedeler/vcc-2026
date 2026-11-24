# 651 — RNA-content (library-size) composition as a filler sign source for fid

2026-10-02 16:51Z (commit 47b018f1; the time was corrected to the commit time, the bar is unchanged), council session [7772dc]. 651 is one of the coordinator's reserved 651–653 (640–649 belongs to
[537757]). Analysis only, < 4 GB, h1-test plus K562 GWPS obs. No builder, slot, `vcc` or git outside this directory.

## How this differs from the closed 610 cpmcouple
- 610 cpmcouple (FINDINGS §2, closed) tested a **uniform per-target** offset: when a target's DE genes lose library mass,
  every other gene's CPM rises together. Null-gene signs follow it (Spearman −0.86), but our key cannot predict the mass
  change, and even the oracle sign rule adds only +0.0014 p_null on our calls.
- 651 tests a **gene-specific** term: β_g × Δ_t.
  - β_g is the slope of gene g's per-cell CPM on log library size across the line's own CONTROL cells (cell size,
    cell cycle, Flex probe saturation). It is measurable on the day from D/E/F's controls.
  - Δ_t is target t's mean log-library shift against the controls.
  - If perturbed cells shift in RNA content, the truth's null-gene lfc carries β_g Δ_t / ln 2, so its signs differ from
    gene to gene within a target. A uniform offset cannot capture that.
- Step 3 is where cpmcouple died: predicting the per-target quantity without the truth. Here the decisive test is
  whether **K562 GWPS's per-target library shift for the same target** (obs `UMI_count`, within gem group) predicts it.

## Measured before this plan (step 1, `lib651.py`, `work/lib651.out`)
- h1-test controls: median L 20,190, sd log L 0.346.
- Per-target Δ_t = mean log L(t) − mean log L(ctrl): sd 0.050 across 85 targets, against a 400-cell sampling se of 0.017.
  14 % of targets exceed 3 se. Extremes are MED19 −0.29, HSP90B1 −0.11, CSNK1E −0.11, EP300 −0.10 and ZBTB39 +0.17,
  PHF12 +0.09: the same targets cpmcouple flagged for DE-mass change, so Δ_t is partly the DE mass itself.
- Batch mix of controls vs perturbed cells: TV distance 0.021 (balanced).

## Definitions (fixed now)
- **Gate genes:** the real DE table's genes (`experiments/213-transport-law/work/de_real_h1.npz`, 10,779 genes, control CPM
  > 5).
- **β_g:** OLS slope of per-cell CPM (counts / L × 1e6) on centred log L over the 18,400 h1-test controls, divided by the
  gene's mean control CPM. This is the relative slope d log(mean CPM) / d log L. **b** = z(β) over the gate genes.
- **Block null genes of target t:** the top 2,200 gate genes by control mean CPM (our CPM block), minus t's own gene,
  restricted to adjudicable genes (finite, non-zero lfc) that the truth does NOT call (padj ≥ 0.05).
- **p_null(t):** the fraction of t's block null genes whose predicted sign equals sign(truth lfc). The panel value is the
  mean over targets. Δ is against the baseline, with a paired target bootstrap (4,000 draws, seed 651).
- **Baseline:** 391's CLEAN4 sign (`clean4_prior_challenge.npz`: P_down > 0.5 → down, else up; uncovered → down).
- **Combined rule:** up-score u_tg = −z(P_C4,g − ½) + λ · b_g · Δ̂_t / s_Δ. Uncovered genes get P_C4 = ½ + 1e-6 (a tie
  broken DOWN, as in 391). s_Δ is the sd of Δ̂ across the predictor's own targets. λ = 1 decides; λ 0.5 / 2 are
  descriptive. Predicted sign = sign(u), with ties going down.

## Step 2 (oracle feasibility, descriptive; does not decide)
- Δ_t is h1-test's observed shift (truth-side, so an oracle).
- Report the per-target regression of null-gene lfc on b_g (slope vs Δ_t), p_null of sign(b_g Δ_t) alone, and the combined
  rule.
- **If the ORACLE combined rule gains < +0.003 p_null over CLEAN4, the family closes here.** A predictor cannot be
  expected to beat its own oracle.

## Step 3 (DECISIVE): Δ̂_t from K562 GWPS, no H1 fitting
- **Δ̂_t:** for each target, the mean over its K562 GWPS cells of log UMI_count minus the mean log UMI_count of the
  non-targeting cells in the same gem group (`batch`). Only targets that GWPS covers are scored.
- **GO only if ALL three hold at λ = 1:**
  1. Spearman(Δ̂_t, h1-test's observed Δ_t) ≥ 0.30 over the covered targets;
  2. Δp_null (combined rule − CLEAN4) ≥ **+0.005**;
  3. its bootstrap 95 % lower bound > 0.
- **Otherwise CLOSED**, with the mechanism written down.
- **What GO buys.** A proposal to the coordinator: a day-of filler sign rider (β_g from D/E/F controls × Δ̂_t from K562 GWPS),
  class TRANSFER. It needs a held-out read before any October use: Orion HCT116 / HEK293T relative only, or the Broad-16
  with the KO caveat. Size if it holds: +0.005 p_null ≈ +0.0048 raw fid at f 0.35 ≈ +0.0027 overall (÷6 after scaling by
  the span 0.30).

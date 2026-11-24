# 285 — The apples-to-apples reliability benchmark: **the guide construct, not the cell count, is what a second measurement of the same perturbation disagrees about.** Cells and UMI depth are ONE removable axis (total UMIs; Spearman–Brown holds); the construct, the donor and the line are not removable by measuring more.

**The one line the mandate asks for.**

> **AT MATCHED CELLS (n = 25 per side, residual cosine): same-guide cell split +0.0063 (K562-GWPS, 748
> two-construct perturbations) / +0.0416 (Nadig Jurkat, 124) / +0.1622 (Nadig HepG2, 70) / +0.0391 (h1-test,
> 85 targets) — thinned counts reproduce the SAME curve at equal total UMIs (h1: 200 cells at ½ depth
> +0.1433 vs 100 cells at full depth +0.1354, Δ +0.0078 [−0.0068, +0.0119]) — different guide +0.0006
> [−0.0002, +0.0014] (K562, 48 % of targets ≤ 0) / +0.0070 (Jurkat) / +0.0659 (HepG2) / +0.0051 (h1) —
> different donor +0.0044 (CD4, matched 2×10⁵ UMI budget ≈ 16 cells/side) — different batch +0.0134 vs
> +0.0120 for cells (K562 gemgroups, i.e. NO batch term; the one exception is Orion HEK293T, 0.0282 vs
> 0.0706; rho 0.42–0.47, flat across a 5x budget ladder) — different study, same line +0.0479 vs +0.0605 within study (K562 GWPS vs the essential screen,
> 1,765 shared perturbations).
> **THE DOMINANT REMOVABLE TERM AT 400 CELLS: NEITHER — cells and depth are the same axis.** At equal total
> UMIs the two are interchangeable to within ±0.008, so there is one removable term (the UMI budget), and at
> 400 cells of one construct it is no longer the binding one: the construct's own ceiling ρ = 0.06–0.41 is
> reached at **n = 76–330 cells per side** on every source where the unit can be measured.
> **DEPTH LAW: measured = Spearman–Brown.** h1-test, per target on a fixed cohort of 85: SB from n = 25
> predicts +0.2455 at n = 200, measured **+0.2392**, measured − SB **−0.0063 [−0.0296, +0.0326]**; K562
> +0.0369 predicted vs +0.0308 measured; Nadig Jurkat +0.1016 vs +0.0973; Nadig HepG2 +0.2954 vs +0.2775.
> The one systematic deviation is a slight shortfall of the measured value above the base rung (≤ 0.018).

**Date 2026-09-06, 13:43–15:52Z (2 h 09 min of the 7-hour box).** Deterministic (bootstrap seed 20260906, 20,000 draws; 8 draws per ladder
rung, 4 per CD4 budget). `uv run python`, `PYTHONDONTWRITEBYTECODE=1` throughout. **No model, no build, no
submission, no `vcc`, no `cell-eval2`, no price, no model name; nothing written under `submissions/`,
`data/`, `scripts/`, `src/` or any sibling experiment; everything this run created is inside
`experiments/285-reliability-benchmark/`.** Preregistration `PREREGISTRATION.md`, sha256
`11f4414edac4a84e15749d84de05a05afe17f071b0da28a4c4b838d2571ca404`, hashed at **13:54:51Z**
(`work/PREREGISTRATION.sha256`) — at that moment only file *metadata* had been read (`work/inv285.json`,
the remote CD4 object's `obs` and `X` layout); no count, profile, cosine or `pds` had been computed on any
source. **The reserves were not touched**: `experiments/267-*/work/reserve/` was never opened and none of
JVE127 / MG63 / SLR23 / UMRC3 appears anywhere in this directory.

**Peak RSS** per job: 12.02 GiB (`pds285.py`), 11.84 (`k562_guide.py`), 10.42 (`k562_kd.py`), 9.27 (`h1_285.py`), 8.21 (`broad`),
8.03 (`csplit`), 6.2 (`cd4`), 6.09 (`k562_pass gwps`), 5.79 (`k562_analyse`), ≤ 7.7 (Nadig). **Three jobs
exceeded the mandate's ~10 GB, by 0.4–2.0 GB** (each holds a 0.8–2.8 GB `uint16` cell bank *and* its
`float32` copy); they ran one at a time with `memory_pressure` at 73–85 % free throughout. Peak disk in
`work/`: **2.5 GB** (cap 20 GB). Bytes transferred from the network: **≈ 0.78 GB** (cap 40 GB), all of it
metadata and layout probes of one CD4 object — see §7.

---

## 0. Known-answer rows — all pass

| row | required | measured | verdict |
|---|---|---|---|
| **KA-1** 267's batch-split residual-cosine medians from its cached `work/half/truth-<fold>.npz` `rel\|batch\|cos_resid` | 0.0e+00 | K562-GWPS **+0.0471**, HCT116 **+0.0276**, HEK293T **+0.0282** | **PASS** |
| **KA-2** 273's h1-test lane-parity 200/200 residual median, RECOMPUTED here from `data/harness/h1-test.h5ad` with this experiment's own code | \|Δ\| ≤ 1e-4 | **+0.242920** vs 273's **+0.242920**, \|Δ\| **0.00e+00** | **PASS** |
| **KA-2b** 159's `resid_disjoint` mean at 200 cells/side (cache read) | 0.0e+00 | **0.330039** | **PASS** |
| **KA-3** 280's Stim48hr donor-half median recomputed from `bank280_X.npy` through `lib280` (3 disjoint 2v2 pairings, 280 targets) | \|Δ\| ≤ 1.5e-3 vs +0.0539 | **+0.0534** | **PASS** |
| **KA-4** the union of two halves reproduces the source's full row (integer identity) | max \|Δ\| = 0 | h1-test lane halves **0.0**; K562 nested rank blocks vs the full sums on the 518 / 766 short-cell targets **0.0** | **PASS** |
| **KA-5** `pds` of a view against itself | 1.0 on every target | **1.000000** (min) | **PASS** |
| **KA-6** 268's Broad-16 halves satisfy `A + B == full` and `v1 + v2 ≤ full` on all 12 lines | 0 / ≤ 0 | **0.0 / 0.0** | **PASS** |
| **KA-probe** the CD4 object's `indptr` reconstructed as `[0, cumsum(obs/n_genes_by_counts)]` | total == `X/data` length AND 0 diff at 40 scattered rows | **11,417,773,379 == 11,417,773,379**; max \|Δ\| **0** | **PASS** |
| **KA-7** (extraction) | — | **NOT REACHED — the extraction was not admissible (§7)** | n/a |

The scorer's statistic is 218's `_per_target`, extracted from `resid218.py` **as source text and exec'd on
its own** (sha256 of the extracted function `bf12ede2…`), so none of 218's module-level code runs here; 273's
KA-2 proved that function equals `cell-eval2` 0.16.0's `pds_cosine` on a cloud-scored artifact to 5.6e-17.

---

## 1. The answer in twelve lines

1. **Cells and UMI depth are one axis, not two.** At equal *total UMIs per side*, halving the depth and
   doubling the cells changes the agreement by nothing measurable: h1-test **+0.0078 [−0.0068, +0.0119]**
   (200 cells at ½ depth vs 100 at full), K562 **+0.0012 [−0.0004, +0.0033]** and **−0.0017 [−0.0046,
   +0.0015]**. So "more cells" and "deeper sequencing" are the same purchase, and the correct budget
   variable is total UMIs per condition.
2. **On that axis the law is Spearman–Brown**, on all four sources with a ladder. Measured − SB is
   −0.0063 [−0.0296, +0.0326] (h1, 25 → 200), −0.0061 [−0.0121, +0.0036] (K562, 25 → 100),
   −0.0043 [−0.0101, +0.0058] (Jurkat), −0.0179 [−0.0371, −0.0004] (HepG2). **This is the corrected
   replacement for 280 §4**: 280's curve thinned counts out of a pseudobulk and never varied the cells; here
   both are varied separately, and they agree.
3. **The guide construct is the dominant irreducible term.** At 25 cells per side on **748** two-construct
   K562-GWPS perturbations, two constructs of the same gene agree at **+0.0006 [−0.0002, +0.0014]** (48 % of
   perturbations ≤ 0) while two halves of ONE construct's cells agree at **+0.0063 [+0.0055, +0.0071]** —
   `CELLS − GUIDE` **+0.0058 [+0.0045, +0.0069]**. In the scorer's own statistic the same swap takes mean
   `pds` from **0.633 to 0.509** (chance = 0.500).
4. **Replicated on four more datasets** (`CELLS − GUIDE`, matched n): Nadig Jurkat **+0.0345 [+0.0269,
   +0.0452]** on 124 targets (`pds` 0.800 → 0.560), Nadig HepG2 **+0.0963 [+0.0674, +0.1447]** on 70,
   h1-test **+0.0096 [−0.0028, +0.0536]** on 5, Broad-16 **+0.0022 [+0.0009, +0.0032]** pooled over 12 lines
   (positive on 9/12), Orion/K562 in 267's cache (guide halves at ~175 cells/side read +0.010 / −0.001 /
   +0.004 against random halves at 77–100 cells of +0.045 / +0.029 / +0.071).
5. **Expressed as a budget**: swapping the construct costs the same as throwing away **90 % of the cells**
   in K562 (25 → 2.4 cells' worth), **84 %** in Jurkat, **64 %** in HepG2, **66 %** in h1-test, **97 %** in
   CD4 — and only **20 %** in Broad-16, the one **Cas9-knockout** dataset here (all the others are CRISPRi).
6. **The batch is free almost everywhere.** Lane / gemgroup halves at matched n equal random cell halves:
   h1 lanes ρ = 0.94–1.00, K562 gemgroups ρ = 0.87–1.12, Nadig ρ = 0.99–1.03, and 267's cached K562 /
   HCT116 rows ρ = 1.05 / 0.94. **The single exception is Orion HEK293T** — and, laddered over a 5× UMI
   budget from 267's cached sums (`orion2_285.py`), it is **flat and irreducible**: ρ = **0.467 [0.282,
   0.697] / 0.473 [0.300, 0.705] / 0.421 [0.288, 0.647]** at 2×10⁵ / 5×10⁵ / 10⁶ UMIs per side, every
   interval below 1, against K562 (1.06 / 1.01 / 0.91) and HCT116 (0.91 / 0.91 / 0.74) on the same ladder
   with the same code. **A batch term, where one exists, behaves like the construct: more measurement does
   not shrink it.**
7. **A different STUDY in the same line is a mild term.** K562-GWPS vs the K562 essential screen (two
   independent Replogle 2022 experiments, 1,765 shared perturbations, matched cells): across-study
   **+0.0479** vs within-study **+0.0605** at n = 25 (ρ **0.791 [0.751, 0.820]**), **+0.0714** vs **+0.0942**
   at n = 50 (ρ **0.758 [0.718, 0.788]**). The scorer's statistic loses only 0.018 `pds` (0.843 → 0.825).
   **A whole independent experiment costs less than a second guide against the same gene.**
8. **A different LINE costs about half.** Broad-16, 66 line pairs at a matched UMI budget: **+0.0051** vs
   within-line cells **+0.0114** (ρ 0.45); at the larger budget **+0.0078** vs **+0.0223** (ρ 0.35).
9. **CD4's donor is as small as its guide, and neither grows with the budget.** Matched-budget rows
   (Rest / Stim8hr / Stim48hr): different donor **+0.005 / +0.004 / +0.004** at 1×10⁵ UMIs and
   **+0.007 / +0.008 / +0.005** at 4×10⁵ — flat — while the same-construct reference rises from +0.11 to
   +0.28 over the same range. Disattenuated, ρ_donor = **0.02–0.07** and ρ_guide = **0.02–0.06**.
   **Across-donor agreement is not the challenge's within-reference repeatability, and it is not a depth
   problem: more UMIs do not move it.**
10. **What DOES reproduce in CD4 is the construct.** Same donor **and** same guide, activation state changed:
    **+0.037 / +0.065 / +0.125** at the three budgets, ρ = **0.49–0.64** — ten times the donor and guide
    terms. This independently reproduces 282's ordering (a donor × guide interaction, not a gene effect)
    with fixed side conventions and a matched budget.
11. **The two readouts agree on the ordering and differ on the scale.** The scorer's `pds` compresses:
    a guide swap that takes the residual cosine from +0.0416 to +0.0070 takes `pds` from 0.800 to 0.560
    (Jurkat), and the across-study swap that costs 21 % of the cosine costs 2 % of `pds`. Where they
    disagree: **`pds` says the batch and the study are almost free and the construct is nearly fatal**,
    which is the same statement the cosine makes, amplified by the rank statistic's knee.
12. **Reproducibility is a property of the response, not of the budget.** Terciles of full-depth response
    magnitude, at fixed n: h1-test **+0.044 / +0.140 / +0.506**, K562 **+0.012 / +0.030 / +0.054**. Between
    sources the same thing: the K562 *panel* targets (non-essential by design) read +0.0095 at n = 25 where
    the K562 *essential* shared targets read +0.0605 at the same n in the same file.

---

## 2. Inventory (`work/inv285.json`, `work/inv285.log`; metadata only, read before the preregistration)

| source | shape | targets | cells/target | units present | rungs used |
|---|---|---|---|---|---|
| h1-test `data/harness/h1-test.h5ad` | 52,400 × 18,080 csr | 85 (+18,400 NTC) | **400 exact** | 48 lanes; 5 two-construct targets | 25 / 50 / 100 / 200 |
| K562-GWPS `…/ReplogleWeissman2022_K562_gwps.h5ad` | 1,989,578 × 8,248 **dense f32, gzip** | 9,866 (272 on the panel) | median 178 (panel 171.5) | 267 gemgroups; **775 two-construct** perturbations (8 on the panel); 75,328 NTC in 514 guides | 25 / 50 / 100 |
| K562 essential `…_K562_essential.h5ad` | 310,385 × 8,563 dense f32 | 2,057 — **2,053 shared with GWPS, 0 on the panel** | median 121 | 48 gemgroups | 25 / 50 |
| Nadig Jurkat / HepG2 `…_raw_singlecell.h5ad` | 262,956 × 8,882 / 145,473 × 9,624 **dense, uncompressed** | 2,393 each, **0 on the panel** | 83 / 45 | 55 / 56 gem groups; 135 multi-guide targets | 25 / 50 |
| Broad-16 (268's `work/halves/<line>.npz`) | 12 lines × 97 shared targets × 34,648 genes | 97 | mean 75 | Cas9 KO, **two guides `_v1`/`_v2`** per target; 12 lines | UMI budgets 2×10⁵ / 5×10⁵ |
| CD4 GWCD4i (280's bank) | 17,495 × 18,129 pseudobulks | 300 panel present | 54 per (guide, donor, condition) row | 4 donors × 3 states × ≥2 guides | UMI budgets 1 / 2 / 4 ×10⁵ |
| Orion HCT116 / HEK293T, K562 (267's cache) | — | 290 / 294 / 269 | 70–100 per half | batch, random, guide halves | as cached |

Two facts settled here that the mandate's table left open: **K562-essential shares 2,053 perturbations with
GWPS and NONE of them is on the 300-target panel** (so the study row is measured on the shared non-panel
set, and is a property of the library and the line, not of the panel); and **Nadig covers 0 panel targets**
in both lines (DepMap common-essential genes; the panel is non-essential by construction, AGENTS.md §2.1).

---

## 3. The response object, the two readouts, and TWO centring corrections found on the way

`r_V(t) = log1p(50,000 · S/ΣS) − (the mean of that statistic over the source's targets, on the same side)`,
target-gene columns removed. Because the control profile is identical for every target of a view, it cancels
in the centring: **the residual rows do not depend on the control at all**, which is why guide / donor /
batch rows that have no matching control pool are still comparable. The scorer's `pds` *does* use the
control, and uses each view's own.

Two artefacts were found by known-answer-style sanity checks and are recorded because both are easy to
repeat:

* **Centring a thinned or sub-sampled side on a FULL-DEPTH mean inflates the cosine.** The first K562 run
  centred every rung on the 172-cell panel mean and returned 1/8-depth agreement (**+0.0667**) ABOVE
  full depth (**+0.0190**) — physically impossible. The depth artefact is a vector shared by both sides.
  Fixed by rebuilding the centring vector for every (n, draw) and every thinning fraction from a matched-n,
  matched-depth sample of every target. Every number in this file uses the fixed convention.
* **A side whose IDENTITY varies across targets makes its own centring vector a mixture, and the mixture
  anti-correlates the two sides.** The first CD4 run assigned "side 0" to whichever donor had more cells;
  centred per side, the four donors' deviations (which sum to zero) put ±½(δ_A − δ_B) into the two sides
  with opposite signs and the donor row read **−0.05 with 70 % of targets ≤ 0**. The second run centred both
  sides on one pooled reference instead, which left each side's own donor deviation in the residual and made
  the *guide* row (both sides in one donor) read ABOVE the count-split ceiling. The correct construction —
  used here — is **a fixed (donor, guide-slot, condition) triple on both sides plus per-side centring**.

---

## 4. The depth law, corrected (`analyse285.py` §1; this replaces 280 §4)

Fixed cohort per source (targets present at every rung), 8 draws averaged per target, medians:

| source | n = 25 | 50 | 100 | 200 | SB from n = 25 → largest | measured − SB |
|---|---:|---:|---:|---:|---:|---|
| **h1-test** (85 targets) | **+0.0391** | +0.0641 | +0.1354 | **+0.2392** | +0.2455 | **−0.0063 [−0.0296, +0.0326]** |
| **K562-GWPS panel** (95) | **+0.0095** | +0.0187 | **+0.0308** | — | +0.0369 | **−0.0061 [−0.0121, +0.0036]** |
| **Nadig Jurkat** (336) | +0.0535 | **+0.0973** | — | — | +0.1016 | −0.0043 [−0.0101, +0.0058] |
| **Nadig HepG2** (176) | +0.1733 | **+0.2775** | — | — | +0.2954 | −0.0179 [−0.0371, −0.0004] |

Read downwards instead (SB from the largest rung) the agreement is even closer: h1 predicts +0.0378 at
n = 25 against +0.0391 measured (**+0.0013 [−0.0065, +0.0061]**) and +0.1358 at n = 100 against +0.1354
(**−0.0004**). **Spearman–Brown in the number of cells is the depth law**, with a slight (≤ 0.018)
optimism when extrapolating upward from the smallest rung.

**Cells vs UMIs at a fixed budget** (§1 line 1) — the pairs are frozen in the preregistration:

| source | more cells at half depth | fewer cells at full depth | difference (same total UMIs) |
|---|---:|---:|---|
| h1-test | 200 cells × ½ = **+0.1433** | 100 cells × 1 = +0.1354 | **+0.0078 [−0.0068, +0.0119]** |
| K562-GWPS | 50 × ½ = +0.0124 | 25 × 1 = +0.0120 | +0.0012 [−0.0004, +0.0033] |
| K562-GWPS | 100 × ½ = +0.0170 | 50 × 1 = +0.0224 | −0.0017 [−0.0046, +0.0015] |

The full thinning ladders are monotone in total UMIs and lie on the same curve as the cell ladder
(h1 at n = 200: f = 1 / ½ / ¼ / ⅛ → +0.2392 / +0.1433 / +0.0802 / +0.0406; K562 at n = 50 → +0.0187 /
+0.0107 / +0.0056 / +0.0038).

---

## 5. The benchmark table (`work/bench285.csv`, 136 rows; medians of the per-target residual cosine, 8 draws)

Headline rows only; the csv carries every cell with its interval, share ≤ 0 and target count.

| source | unit | cells/side | n targets | median cos | 95 % CI | share ≤ 0 | mean `pds` | `cos_agg` |
|---|---|---:|---:|---:|---|---:|---:|---:|
| **h1-test** (H1 hESC, Flex) | a cells, same guide | 25 | 85 | +0.0391 | [+0.0216, +0.0592] | 0.024 | 0.839 | — |
| h1-test | a cells | 100 | 85 | +0.1354 | [+0.0781, +0.1957] | 0.000 | 0.936 | — |
| h1-test | a cells | 200 | 85 | +0.2392 | [+0.1426, +0.3312] | 0.000 | 0.967 | — |
| h1-test | **d lane** | 25 / 50 / 100 | 85 | +0.0343 / +0.0728 / +0.1293 | — | ≤ 0.012 | — | — |
| h1-test | **b guide construct** | 37 | 5 | **+0.0051** | — | 0.40 | — | — |
| h1-test | a cells, one construct | 37 | 5 | +0.0147 | — | 0.00 | — | — |
| **K562-GWPS**, panel | a cells (fixed cohort) | 25 / 50 / 100 | 95 | +0.0095 / +0.0187 / +0.0308 | [+0.0077,+0.0114] … | 0.08 / 0.07 / 0.02 | 0.696 / 0.775 / 0.857 | 0.031 / 0.045 |
| K562-GWPS, panel | **d gemgroup** | 25 / 50 / 100 | 255 / 217 / 81 | +0.0134 / +0.0229 / +0.0286 | — | 0.03–0.06 | 0.736 / 0.817 | 0.031 / 0.046 |
| **K562-GWPS, 748 two-construct perturbations** | **b guide construct** | 25 | 748 | **+0.0006** | [−0.0002, +0.0014] | **0.476** | **0.509** | 0.010 |
| K562-GWPS, the same 748 | a cells, one construct | 25 | 748 | **+0.0063** | [+0.0055, +0.0071] | 0.217 | **0.633** | 0.027 |
| K562-GWPS, the same 653 | b guide / a cells | 50 | 653 | +0.0006 / +0.0102 | — | 0.48 / 0.18 | 0.510 / 0.689 | 0.006 / 0.041 |
| **K562 GWPS vs essential** | a cells within study | 25 | 1,765 | +0.0605 | [+0.0561, +0.0655] | — | 0.843 | 0.119 |
| K562 GWPS vs essential | **e different STUDY, same line** | 25 | 1,765 | **+0.0479** | [+0.0440, +0.0513] | — | **0.825** | 0.117 |
| K562 GWPS vs essential | a cells / e study | 50 | 1,163 | +0.0942 / +0.0714 | — | — | 0.866 / 0.849 | 0.188 / 0.179 |
| **Nadig Jurkat** | a cells / d gemgroup | 25 | 534 / 528 | +0.0623 / +0.0626 | — | — | 0.857 / 0.862 | 0.132 / 0.134 |
| Nadig Jurkat | **b guide** vs a cells, one construct | 25 | 124 | **+0.0070** vs +0.0416 | — | — | **0.560** vs 0.800 | 0.020 vs 0.118 |
| **Nadig HepG2** | a cells / d gemgroup | 25 | 507 / 494 | +0.1841 / +0.1801 | — | — | — | — |
| Nadig HepG2 | **b guide** vs a cells, one construct | 25 | 70 | **+0.0659** vs +0.1622 | — | — | — | — |
| **Broad-16**, 12 lines pooled (Cas9 KO) | a cells / **b guide** / **f LINE** | T = 2×10⁵ UMIs | ~1,000 rows | +0.0114 / **+0.0092** / **+0.0051** | [+0.0102,+0.0128] … | — | — | — |
| **CD4 GWCD4i** Stim48hr | a* count split (Poisson only) | T = 2×10⁵ (~16 cells) | 55 | +0.1749 | — | 0.000 | — | — |
| CD4 Stim48hr | **b guide** / **c donor** / **s state** | T = 2×10⁵ | 35 / 91 / 92 | **+0.0039 / +0.0044 / +0.0731** | — | 0.37 / 0.31 / 0.00 | — | — |
| CD4 Stim48hr | b / c / s | T = 4×10⁵ (~33 cells) | — / 51 / 44 | — / **+0.0054** / +0.1248 | — | 0.35 / 0.00 | — | — |
| **267 cache** K562 / HCT116 / HEK293T | a random cells | 86 / 77 / 100 | 269 / 290 / 294 | +0.0450 / +0.0293 / **+0.0706** | — | 0.06 / 0.30 / 0.15 | — | — |
| 267 cache, the same | **d batch** | 89 / 74 / 100 | " | +0.0471 / +0.0276 / **+0.0282** | — | 0.06 / 0.32 / 0.34 | — | — |
| 267 cache, the same | **b guide** | 175 / 172 / 196 | 8 / 10 / 10 | +0.0103 / −0.0010 / +0.0042 | — | 0.50 | — | — |

---

## 6. The decomposition: each unit's OWN ceiling (`analyse285.json` `rho`, `n_stop`)

The observed agreement of any unit at n cells per side is `ρ_unit · a(n)`, where `a(n)` is the cell-sampling
attenuation; the same-construct cell split has `ρ = 1`, so **the cell-split row IS `a(n)`** and
`ρ_unit = cos_unit(n) / cos_cells(n)` at the same n. `ρ` is what more measurement can never change.

| unit | source (n or budget) | ρ (95 % CI) | equivalent to keeping … of the budget |
|---|---|---|---|
| **d batch / lane / gemgroup** | h1-test (25 / 50 / 100) | **0.94** [0.72, 1.10] / **0.98** [0.78, 1.11] / **1.00** [0.86, 1.07] | ~100 % — free |
| d gemgroup | K562-GWPS (25 / 50 / 100) | 1.12 / 1.00 / 0.87 | ~100 % — free |
| d gemgroup | Nadig Jurkat / HepG2 (25) | 1.02 [0.95, 1.07] / 0.99 [0.96, 1.04] | ~100 % — free |
| d batch | 267 K562 / HCT116 (89 / 74 cells) | 1.05 / 0.94 | ~100 % — free |
| d batch, laddered | K562 (2×10⁵ / 5×10⁵ / 10⁶ UMIs) | 1.06 [0.83, 1.18] / 1.01 / 0.91 | ~100 % — free |
| d batch, laddered | Orion HCT116 (same ladder) | 0.91 [0.51, 1.30] / 0.91 / 0.74 | ~100 % — free |
| **d batch, laddered** | **Orion HEK293T (same ladder)** | **0.467 [0.282, 0.697] / 0.473 [0.300, 0.705] / 0.421 [0.288, 0.647]** | **~45 %** — the one real batch term, and FLAT in the budget |
| **e different STUDY, same line** | K562 GWPS vs essential (25 / 50) | **0.791 [0.751, 0.820]** / 0.758 [0.718, 0.788] | **78 %** |
| **f different LINE** | Broad-16 (2×10⁵ / 5×10⁵ UMIs) | **0.45 / 0.35** | **45 %** |
| **s different STATE** (same donor + guide) | CD4 (1 / 2 / 4 ×10⁵) | **0.53 / 0.49 / 0.53** | **52 %** |
| **b guide construct** | Broad-16, **Cas9 KO** (2×10⁵ / 5×10⁵) | **0.81 / 0.76** | 80 % |
| b guide construct | Nadig HepG2 (25) | **0.41 [0.20, 0.55]** | 36 % |
| b guide construct | h1-test (37) | 0.35 [−1.72, +1.43] — 5 targets | 34 % |
| b guide construct | Nadig Jurkat (25 / 50) | **0.17 [0.09, 0.29]** / 0.19 [0.03, 0.44] | 16 % |
| **b guide construct** | **K562-GWPS, 748 perturbations (25 / 50)** | **0.090 [−0.038, +0.232]** / 0.062 | **9.5 %** |
| b guide construct | CD4 (1 / 2 ×10⁵) | **0.02–0.06** | 3 % |
| b guide construct | K562-GWPS **panel targets only**, 8 targets | −0.78 [−13.2, +8.2] — uninformative | — |
| **c different DONOR** | CD4 (1 / 2 / 4 ×10⁵) | **0.02–0.07** | 3 % |

**The n at which more cells of ONE construct stop helping** (solve `a(n) = ρ` through Spearman–Brown from
the measured base rung): **89 cells/side** (Nadig Jurkat, ρ 0.19), **76** (HepG2, ρ 0.41), **328** (h1-test,
ρ 0.35). Below the construct's ceiling the budget still buys agreement; above it, a second 400-cell
measurement of the *same* construct keeps improving while a measurement of the *same gene* through another
construct does not.

**Calibrating the CD4 bank** (`csplit_285.py`). CD4 has no cells, so its cell-sampling reference is a
molecular count split, which is Poisson-only and therefore optimistic. On the two sources that have both,
at matched total UMIs, the **count-minus-cell offset is almost constant**: +0.0490 / +0.0445 / +0.0459 /
+0.0443 (h1 at 25/50/100/200) and +0.0398 / +0.0415 / +0.0416 (K562 at 25/50/100), mean **+0.0438**
(the *ratio* runs 0.25 → 0.85 with n, so the additive form is the stable one). The CD4 ρ's above use the
count split minus 0.0438.

### 6b. Is the guide term knockdown EFFICIENCY? Partly — and better guides do not remove it (`k562_kd.py`)

If two constructs disagree only because one of them knocks the gene down weakly, the guide term is a dose
artefact that guide selection removes. Measured on the 539 K562-GWPS two-construct perturbations whose own
gene is on the axis: each construct's **on-target log delta at its own gene column** (a column excluded from
every scored statistic) has median **−0.320**, the *weaker* construct of each pair **−0.119**, and both
constructs are negative on **79 %** of pairs. Stratified by the weaker construct's knockdown, at n = 25:

| tercile of the weaker construct's knockdown | n | weaker-construct on-target Δ | guide | cells, one construct | ρ |
|---|---:|---:|---:|---:|---|
| T1 weakest (no detectable knockdown) | 180 | **+0.015** | +0.0009 [−0.0010, +0.0022] | +0.0072 | 0.124 [−0.136, +0.320] |
| T2 | 179 | −0.119 | +0.0003 [−0.0008, +0.0018] | +0.0064 | 0.041 [−0.125, +0.293] |
| **T3 strongest** | 180 | **−0.300** | **+0.0021 [+0.0003, +0.0031]** | +0.0089 | **0.236 [+0.029, +0.416]** |

**Reading.** Efficiency matters — ρ roughly doubles from the weakest to the strongest tercile — but at the
strongest knockdown the two constructs still agree at only **ρ ≈ 0.24**, and the improvement is not specific
to the guide unit: knockdown strength predicts the *cell-split* agreement exactly as strongly
(Spearman **+0.108**) as it predicts the *guide* agreement (**+0.110**), i.e. what the tercile is really
selecting is response magnitude (§9), not off-target cleanliness. **Selecting effective guides raises the
construct ceiling but does not abolish it.**

---

## 7. The bounded CD4 cell-level extraction: **NOT ADMISSIBLE — the plan, the layout and the exact reason**

Object probed: **`D2_Stim48hr.assigned_guide.h5ad`**, `137,554,355,000 B`, over HTTP range requests with
280's `remote_h5.RemoteFile`.

**What the probe established** (`work/probe285.py`, `probe2_285.py`, `probe3_285.py`, `filt285.py`,
`layout285.py`, `findbase285.py`, `findbase2_285.py`; ≈ 0.75 GB transferred in total, nothing but metadata):

* `X` **is** `csr_matrix`, **`nfilters 0`, no compression**, `data` float32 and `indices` **int64**, both
  chunked `(7245,)`; 2,863,571 cells × 18,130 genes, 11,417,773,379 nnz, ~3,986 nnz/cell (12 B/nnz).
* `obs`: `perturbed_gene_name` (12,084 labels, **290 of the 300 panel targets present, 48,055 panel cells,
  median 99.5/target**), `guide_id` (22,864), `guide_group` (`targeting single sgRNA` 2,062,521 /
  `multi sgRNA` 508,111 / `no sgRNA` 292,939), `guide_type` (targeting 2,775,412 / non-targeting 88,159),
  `lane_id` (24), `low_quality` (18,252 true), `total_counts`, `n_genes_by_counts`.
* The QC-clean single-guide panel cohort is **47,842 cells over 290 targets, 239 of them with two
  constructs**, median 99 cells/target and a median smaller-construct size of 19 cells.
* **The plan**: 47,842 panel + 15,000 non-targeting = **62,842 rows, 244,153,523 nnz, 2.93 GB payload,
  a 1.95 GB kept sparse bank** — comfortably inside the 25 GB admissibility bar and the 20 GB disk cap.
* `indptr` = `[0, cumsum(obs/n_genes_by_counts)]` **exactly** (KA-probe: the total equals `X/data`'s length
  to the byte and 40 scattered rows agree at 0), so every row's *element* range is known.

**Why it was still not admissible.** Turning an element range into a byte range needs the chunk layout, and
**this object's chunks are NOT laid out linearly** — unlike 280's pseudobulk object, where `layout280.py`
verified `offset = base + i · chunk_bytes` and the whole extraction followed. An 8-byte-aligned byte-pattern
scan over the file (int64 values in [0, 18130) vs float32 positive integers, 80 probes) shows `X/data` and
`X/indices` **interleaved in multi-gigabyte runs**, and `get_chunk_info(0)` on a 1.57 M-chunk index did not
return in 10 minutes over HTTP. The only remaining route is h5py's own chunk walk, and it was benchmarked:
**1.53 s per scattered CSR row** (30 rows in file order, 40 range requests, 42 MB) — **≈ 27 hours
single-process for 62,842 rows**, against a 7-hour box. Reading `X/indptr[:]` through the same path had
already run > 20 minutes before being abandoned.

**Recorded, per the mandate**: bytes transferred **≈ 0.78 GB** (probe2 318 MB / 71 requests, probe3 287 MB /
69, the filter check 36 MB / 9, the row benchmark 50 MB, the two byte-pattern scans ≈ 7 MB in ~300 8 KB
probes, the abandoned first probe ≈ 78 MB), plus one 23 MB bulk read that returned garbage and is itself the
proof of the non-linear layout, **kept on disk**: `work/indptr_D2_Stim48hr.npy`
(22.9 MB), `work/probe2_D2_Stim48hr.npz` (51.4 MB, the object's `obs`), `work/plan_D2_Stim48hr.npy` (0.5 MB).
**No expression value was ever downloaded**, no second donor was attempted, and the CD4 rows in this file
are the bank-only units, as the preregistration's fallback directs.

*If a future run wants those cells*: the plan is computed and stored; what it needs is either the chunk-offset
map (one pass of `get_chunk_info` over 1.57 M chunks, feasible on a local copy but not over HTTP) or a
different object — the per-lane GEO CellRanger matrices (`GSE314342_RAW.tar`, 284 files, 0.04–0.76 GB each)
are small enough to fetch individually, at the cost of re-deriving the guide assignment (280 §1.1 rejected
that route for the same reason).

---

## 8. Where the two readouts disagree

The scorer's `pds` is a rank statistic with a knee (273 §3), so it compresses everything above the knee and
expands everything at it. Measured here on the same views:

| contrast (matched n) | residual cosine | mean `pds` | `cos_agg` |
|---|---|---|---|
| K562 748 perts, cells one construct → **other construct** (n = 25) | +0.0063 → **+0.0006** (−90 %) | 0.633 → **0.509** (−0.124) | 0.027 → 0.010 |
| Nadig Jurkat 124, cells one construct → **other construct** (n = 25) | +0.0416 → **+0.0070** (−83 %) | 0.800 → **0.560** (−0.240) | 0.118 → 0.020 |
| K562 1,765 shared, within study → **other study** (n = 25) | +0.0605 → +0.0479 (−21 %) | 0.843 → 0.825 (**−0.018**) | 0.119 → 0.117 |
| K562 panel, random cells → **gemgroup halves** (n = 25) | +0.0120 → +0.0134 | 0.731 → 0.736 | 0.031 → 0.031 |
| Nadig Jurkat, random cells → **gemgroup halves** (n = 25) | +0.0623 → +0.0626 | 0.857 → 0.862 | 0.132 → 0.134 |

They agree on every sign and on the ordering. They disagree on scale exactly where 273 said they would: a
21 % loss of cosine that stays above the knee costs 2 % of `pds`, while an 83 % loss that crosses it costs
30 % of `pds` and lands at chance.

---

## 9. Strata (`analyse285.json` `strata`)

Terciles of the full-depth residual row norm, at a fixed n:

| source | T1 weak | T2 mid | T3 strong |
|---|---:|---:|---:|
| h1-test, n = 100 | +0.0442 (29) | +0.1402 (28) | **+0.5063 (28)** |
| K562-GWPS panel, n = 50 | +0.0116 (75) | +0.0297 (74) | +0.0540 (74) |

And between cohorts inside one file: the K562 **panel** targets read +0.0095 at n = 25 while the K562
**essential** shared targets read +0.0605 at the same n with the same code — a 6× difference that has
nothing to do with the budget. This is the same fact 280 reported (magnitude, not cells, predicts
reproducibility) measured with cells held fixed.

---

## 10. Sensitivity and multiplicity

* **Leave-one-lane-out** (h1-test, n = 100, 48 lanes dropped one at a time): median-of-medians **+0.1385**,
  range **[+0.1221, +0.1525]** against the all-lane value +0.1354 — no lane drives the ladder.
* **Leave-one-donor-pair-out / leave-one-condition-out** (CD4): every donor row above is the mean of the
  **6 fixed donor pairs**, every state row the mean of the **fixed condition pairs across 4 donors**; the
  per-pair rows are in `work/cd4_285.json`, and the spread across pairs is smaller than the interval.
* **Per-line** (Broad-16): `CELL − GUIDE` is positive on **9/12** lines at 2×10⁵ and 8/12 at 5×10⁵; the
  pooled interval excludes zero.
* **Multiplicity**: **136 rows** in `bench285.csv` (source × unit × rung) plus **44** ρ cells, 12
  budget-equivalents, 9 Orion/K562 budget-ladder cells and 3 knockdown terciles. The unit list, the rungs, the thinning fractions and the frozen contrast pairs were
  all fixed in `PREREGISTRATION.md` before any expression value was read; nothing was added after a look,
  and the two design changes made during the run (§3) were *corrections of artefacts detected by
  impossibility*, applied uniformly to every row, not selections among outcomes.

---

## 11. Honest limits

* **The guide unit is measured on the targets that happen to carry two constructs**, which is 1–10 % of any
  library and is not a random sample of targets (in K562 the 775 two-construct perturbations are mostly not
  panel targets). The direction of the result is consistent across five datasets and two chemistries, but a
  library that deliberately paired constructs for *hard* targets would bias ρ_guide down. §6b bounds one
  version of that worry directly — at the strongest on-target knockdown ρ is still only 0.24 — but it cannot
  exclude that a library built from only well-validated constructs would read higher still.
* **ρ is a ratio of two medians, and its denominator is small on the weak-response sources.** The
  `CELLS − GUIDE` differences (tight, all reported) are the robust statistic; ρ is the interpretable one,
  and on K562 its interval [−0.04, +0.23] is wide even with 748 perturbations. The h1-test ρ (5 targets) and
  the K562 *panel* ρ (8 targets) are uninformative and are reported only for completeness.
* **Broad-16 is Cas9 knockout, not CRISPRi, and its cell-split reference mixes the two constructs**
  (268 built `A`/`B` from all assigned cells). Both differences push in the direction of the high ρ_guide
  it reports; the contrast with the four CRISPRi datasets is a hypothesis about the perturbation modality,
  not a controlled comparison.
* **CD4 has no cells here.** Its cell-sampling reference is a count split corrected by an offset measured on
  two *other* sources; the correction is stable (+0.040…+0.049 across 7 cells) but it is a transfer.
  The CD4 donor and guide rows themselves need no correction — they are flat in the budget, which is the
  claim that matters.
* **Matching is on cells where cells exist and on total UMIs where they do not.** §4 shows the two are
  interchangeable to ±0.008 on the two sources where both were measured; that is the licence for the CD4
  and Broad-16 rows, and it is an assumption there.
* **Nothing here is measured on A/B/C or on the challenge's own reference.** Every source is a public
  screen; h1-test is the challenge's chemistry and cell count but an embryonic line (273's warning).
* **Two jobs peaked at 11.8 and 12.0 GiB**, 1.8–2.0 GB above the mandate's ~10 GB, in the step that holds a
  large `uint16` cell bank and its `float32` copy at the same time. They ran serially with the machine at
  73–85 % free; no sibling job was affected as far as `memory_pressure` shows.
* **`pds` is reported only where the cohort is large enough for a rank statistic** (≥ ~100 targets); the
  5- and 8-target guide rows carry no `pds`.

---

## 12. Files, memory, provenance

| file (in `work/`) | what | peak RSS / time |
|---|---|---|
| `inv285.py` → `inv285.{json,log}` | the metadata inventory (§2), read before the preregistration hash | 2.4 GB / 4 s |
| `lib285.py` | the response object, the two readouts, the nested ladder, the bootstrap; 218's `_per_target` extracted verbatim (`bf12ede2…`) | — |
| `h1_285.py` → `h1_285.{npz,json,log}` | h1-test: KA-2/2b/4/5, the ladder, thinning, lanes, guides, leave-one-lane-out | 9.27 GiB / 41 s |
| `k562_pass.py gwps\|essential` → `k562_*_{panelX.npy,meta.npz,study_blocks.npy,study_full.npy}` | the two streaming passes (KA-4a) | 6.09 / 4.04 GiB; 226 / 117 s |
| `k562_analyse.py` → `k562_285.{npz,json,log}` | KA-1, the K562 ladder, thinning, gemgroups, guides, the across-study rows | 5.79 GiB / 43 s |
| `k562_guide.py` → `k562_guide.{npz,json,log}` | the guide row on **748** two-construct perturbations, both readouts | 11.84 GiB / 168 s |
| `csplit_285.py` → `csplit_285.{json,log}` | the count-split vs cell-split calibration (§6) | 8.03 GiB / 19 s |
| `cd4_285.py` → `cd4_285.{npz,json,log}` | KA-3 and the CD4 bank units at three UMI budgets | 6.21 GiB / 184 s |
| `broad_285.py` → `broad_285.{json,log}` | KA-6, Broad-16 cells vs guides vs lines at matched budgets | 8.21 GiB / 52 s |
| `nadig_285.py jurkat\|hepg2` → `nadig_*.{npz,json,log}` | Nadig cells / gemgroups / guides | 7.73 / 6.10 GiB; 11 / 6 s |
| `orion_285.py` → `orion_285.{json,log}` | 267's cached batch / random / guide rows, paired | < 1 GiB / 2 s |
| `pds285.py` → `pds285.{json,log}` | the scorer's statistic on the headline contrasts (§8) | 12.02 GiB / 21 s |
| `k562_kd.py` → `k562_kd.{npz,json,log}` | the on-target-knockdown stratification of the guide term (§6b) | 10.42 GiB / 148 s |
| `orion2_285.py` → `orion2_285.{json,log}` | the batch-vs-cells ladder on 267's cached Orion / K562 sums (mmapped) | 1.94 GiB / 27 s |
| `analyse285.py` → `analyse285.{json,log}`, **`bench285.csv`** | the depth law, cells-vs-UMIs, the strata, the ρ table, `n_stop` | < 1 GiB / 3 s |
| `probe285.py`, `probe2_285.py`, `probe3_285.py`, `filt285.py`, `layout285.py`, `findbase285.py`, `findbase2_285.py`, `bench285.py` → `probe*_D2_Stim48hr.{json,log}`, `indptr_D2_Stim48hr.npy`, `plan_D2_Stim48hr.npy`, `layout2_D2_Stim48hr.json` | the CD4 extraction probe and its verdict (§7) | < 1 GiB; ≈ 0.75 GB transferred |
| `equiv285.{txt,json}` | the budget-equivalent of each unit swap (§1 line 5) | — |

**Write audit.** `find experiments data scripts src docs submissions -newer work/.marker-start -type f`
returns, outside this directory, only the sibling agents' own trees (284, 286–291), the coordinator's
leaderboard snapshot `data/leaderboard/2026-09-06T1528Z.json`, `docs/`, and one
`experiments/265-*/__pycache__/boot265.cpython-312.pyc` — which no script here can have produced (nothing in
this directory imports `boot265`, and every run set `PYTHONDONTWRITEBYTECODE=1`).

Read-only inputs: `data/harness/h1-test.h5ad`, `data/external/vcc2025/test/adata_Test.h5ad` (obs only),
`data/external/scperturb/ReplogleWeissman2022_K562_{gwps,essential}.h5ad`,
`data/external/nadig2025/{jurkat,hepg2}_raw_singlecell.h5ad`, `data/challenge/pert_counts.csv`,
`experiments/267-replicate-validated-oracle/work/half/truth-<fold>.npz` (**not** `work/reserve/`),
`experiments/268-broad16-context-residuals/work/halves/*.npz`,
`experiments/273-flex-reference-instrument/work/reliability273.json`,
`experiments/280-gwcd4i-donor-replicates/work/{lib280.py, remote_h5.py, obs280.npz, bank280_*.npy, axis280.npz}`,
`experiments/218-independent-donors/resid218.py` (source text only), and the remote object
`s3://genome-scale-tcell-perturb-seq/marson2025_data/D2_Stim48hr.assigned_guide.h5ad` (metadata only).
No builder import, no `plolo` fold, no operator, no `vcc`, no `cell-eval2`, no submission, no price.

---

## 13. Verdict

**The decomposition, at the challenge's own budget: the removable term is one axis (total UMIs, obeying
Spearman–Brown), and it stops being the binding constraint between 76 and 330 cells per side. What remains
is the guide construct.** Two constructs against the same gene, measured at the same number of cells in the
same line, agree at ρ = 0.09 (K562, 748 perturbations), 0.17–0.19 (Jurkat), 0.41 (HepG2), 0.02–0.06 (CD4),
0.76–0.81 (Broad-16, Cas9 KO); in the scorer's own statistic the swap moves mean `pds` from 0.633 to 0.509
and from 0.800 to 0.560 — to chance. Against that, a different processing batch is free (ρ ≈ 1 on five
datasets; the one exception is Orion HEK293T at ρ 0.42–0.47, flat across a 5× budget ladder), a whole independent experiment in the same line costs
21 % of the cosine and 2 % of `pds` (ρ 0.79), a different line costs 55–65 % (ρ 0.35–0.45), and in CD4 the
donor is indistinguishable from the guide (ρ 0.02–0.07) while the same donor + same guide across an
activation state change keeps half (ρ 0.53). **So the challenge's truth — 400 cells of one construct per
target — is a well-measured quantity, and the thing it measures is that construct's response, not the
gene's.** Measurement only. No arm, no dose, no price, no build, no submission.

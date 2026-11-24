# 266 — 260's NO-GO revisited narrowly: batch-matched controls do NOT change K562-GWPS's real Wilcoxon DE-set sizes on the panel

**Date** 2026-09-05, 12:37–13:00 UTC (≈ 23 min wall clock of the 60; machine clock is UTC+3). **Measurement only.** No harness, no `cell-eval2`, no
build, no `vcc`, nothing under `submissions/` or `data/harness/`, no file outside this directory modified.
Scripts and outputs in `work/`: `extract266.py` (one streaming pass, 1.7 min, peak RSS **3.22 GiB**, `extract266.log`
/ `.json`; caches `counts266.npy` 2.94 GB float32 raw counts of all 75,328 control cells + the 13,654 cached seed-0
cells of 260's 50 candidates, `rows266.npz`), `pilot266.py` (450 Wilcoxon runs, 8.4 min on 6 forked workers,
`/usr/bin/time -l` max RSS **3.91 GiB** for the largest worker — that figure includes the memmap's file-backed pages,
which the six workers share; `memory_pressure` stayed at 70 % free throughout; `pilot266.log` / `.json`,
`runs266.csv`, `de_tables266.npz` = per-run p-values and LFCs), `analyse266.py` → `tables266.md`, `analysis266.json`,
`analyse266.log`. Both stages under the ~6 GB ceiling; the float64 CPM temporaries that took 260 to 8.25 GiB are gone
(CPM is float32 in place, per 512-gene chunk). The `.npy`/`.npz` caches in `work/` are git-ignored by the repo's `.gitignore:11-12`; `counts266.npy` regenerates from the raw file in 1.7 min via `extract266.py`, and the tracked `runs266.csv` / `analysis266.json` / `tables266.md` carry every number below.

---

## VERDICT LINE

**NO — batch-matched control selection does not change 260's coverage verdict.** With controls drawn from each
target's own gemgroups in proportion (composition TV distance 0.39 → 0.005) and additionally UMI-matched within
gemgroup (control-pool UMI CV 0.440 → 0.378, depth-ratio outliers 23/150 → 3/150), the 50 best-covered candidates'
median significant-gene count is **6.0 / 5.7 / 6.3** (pooled / gemgroup-matched / gemgroup + depth-matched,
seed-averaged) against 260's 6, and **17 / 17 / 18 of 50 (34–36 %)** clear `nmae`'s 10-gene gate against 260's 17;
the reopening thresholds were **≥ 20** and **≥ 60 %**. The per-target paired change has median **0** under both
matched schemes (up/down 16/18 and 21/12 of 50), the seed spread of the matched schemes is the same size as the
control-draw noise of the pooled scheme itself, and where sets do change the "new" genes are mostly the pooled
ordering's next genes at q ≈ 0.05–0.07. The known-answer row (260's exact control draw) reproduces all 50 of 260's
counts to the digit. **260's NO-GO stands on a cleaner measurement: the sparse DE sets are the panel's biology in
K562, not a control-pool artefact.**

---

## 1. Design

- **Targets and cells**: 260's 50 candidates (`experiments/260-raw-de-diagnostic-feasibility/work/candidates_de_power.csv`),
  the SAME perturbed cells 260 used (its seed-0 rows cached in `work/extract_rows.npz`, 200–400 cells each, cap 400).
- **Recipe**: 260's `extract_and_de.py` arithmetic verbatim — per-cell CPM as float32 counts / float32 row total × 1e6
  (bitwise 260's values), genes with control mean CPM > 5 on THAT run's control pool (8,248 / 8,248 in all 450 runs;
  minimum control-mean CPM over all runs 5.28), two-sided Mann–Whitney U with the normal approximation, tie-corrected
  (`tc`, the scorer's convention, `cell_eval2/de_compute.py:871-872`) and uncorrected (`uc`), BH per target, α 0.05,
  target gene excluded, `gate` = significant with finite non-zero LFC. The ranks are computed exactly (vectorised
  midranks over tie runs, `pilot266.py`), so R1 and the tie term are 260's to the last bit; z, p and BH use 260's
  formulas unchanged. Means for the LFC are float64 here (260: float32) — sign and finiteness only.
- **Three control schemes, 18,400 controls each, three seeds each** (`runs266.csv`, 50 × 3 × 3 = 450 runs):
  - `pooled`: uniform over all 75,328 controls; **seed 0 = 260's exact draw (the KA row)**, seeds 1–2 fresh draws.
  - `gem`: for each target, its cells' gemgroup counts (`obs.batch`, 267 gemgroups; a target's 200–400 cells lie in
    131–207 gemgroups, median 167, 1–3 cells each) → controls apportioned to those gemgroups in proportion
    (Hamilton rounding, capped by availability, shortfall redistributed; every draw reached 18,400), drawn without
    replacement within gemgroup.
  - `gemdepth`: same apportionment; within each gemgroup the controls are the nearest-UMI unused matches to a
    jittered quantile grid of the target's own cells' UMI in that gemgroup (greedy nearest neighbour, random order).
- **Diagnostics per run**: control-pool UMI CV, per-target depth ratio (target median UMI / control median UMI,
  `check_harness_depth.py`'s band [0.95, 1.05]), per-cell zero-fraction gap (pert − ctrl, the quantity the rank
  test reacts to, `docs/HARNESS-DEPTH-CONFOUND.md`), total-variation distance between the target's and the pool's
  gemgroup composition, number of gemgroups in the pool.

## 2. Known-answer row — pooled seed 0 reproduces 260 on all 50 targets

| field | identical / 50 |
|---|---:|
| `n_real_tc` | **50** |
| `n_real_uc` | 50 |
| `gate_tc`, `gate_uc`, `cells`, `median_umi` | 50 each |
| `frac_neg_tc` (|Δ| < 1e-9) | 50 |
| `target_knockdown_lfc` (|Δ| < 1e-4, float64 means) | 50 |

Control median UMI 11,406 (260: 11,406). Spot-checked before the run on DCP2 346, NABP2 9, MED15 0, PPIF 24
(`n_real_tc` and `n_real_uc` both to the digit). The exact-rank implementation is therefore 260's test.

## 3. The three-scheme table (`work/tables266.md` has the full per-target table; all 50 rows below)

Per target `n_real` (tie-corrected, target gene excluded); the last two columns are the seed-mean differences:

| target | cells | pooled s0 (=260) | pooled s1 | pooled s2 | gem s0 | s1 | s2 | gem+depth s0 | s1 | s2 | gem − pooled | gem+depth − pooled |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| DCP2 | 241 | 346 | 389 | 387 | 377 | 379 | 370 | 315 | 314 | 313 | +1.3 | −60.0 |
| ANP32B | 260 | 134 | 154 | 136 | 152 | 147 | 151 | 179 | 174 | 176 | +8.7 | +35.0 |
| DNTTIP1 | 267 | 113 | 109 | 121 | 120 | 132 | 125 | 120 | 120 | 120 | +11.3 | +5.7 |
| FCGRT | 224 | 103 | 93 | 97 | 85 | 90 | 92 | 107 | 105 | 107 | −8.7 | +8.7 |
| HMGXB3 | 342 | 85 | 92 | 85 | 95 | 91 | 96 | 101 | 102 | 102 | +6.7 | +14.3 |
| PHF23 | 261 | 80 | 73 | 72 | 71 | 70 | 79 | 88 | 88 | 90 | −1.7 | +13.7 |
| EPC2 | 259 | 77 | 76 | 80 | 83 | 82 | 77 | 79 | 81 | 79 | +3.0 | +2.0 |
| MED13 | 268 | 65 | 60 | 57 | 68 | 63 | 58 | 72 | 69 | 67 | +2.3 | +8.7 |
| PCGF1 | 260 | 64 | 64 | 63 | 68 | 75 | 69 | 68 | 69 | 66 | +7.0 | +4.0 |
| FOXN2 | 232 | 63 | 64 | 65 | 66 | 67 | 69 | 71 | 70 | 70 | +3.3 | +6.3 |
| PAN2 | 215 | 59 | 55 | 57 | 56 | 57 | 53 | 60 | 60 | 60 | −1.7 | +3.0 |
| GIGYF2 | 251 | 54 | 49 | 51 | 60 | 56 | 61 | 66 | 66 | 64 | +7.7 | +14.0 |
| ZNF462 | 300 | 50 | 63 | 55 | 61 | 61 | 57 | 63 | 66 | 64 | +3.7 | +8.3 |
| E2F1 | 226 | 31 | 26 | 26 | 23 | 24 | 25 | 30 | 29 | 30 | −3.7 | +2.0 |
| FOXO4 | 293 | 30 | 25 | 29 | 24 | 28 | 21 | 27 | 26 | 27 | −3.7 | −1.3 |
| ZNF507 | 279 | 25 | 25 | 26 | 24 | 23 | 23 | 21 | 23 | 21 | −2.0 | −3.7 |
| PPIF | 256 | 24 | 20 | 19 | 21 | 17 | 20 | 20 | 20 | 21 | −1.7 | −0.7 |
| GPATCH8 | 284 | 9 | 10 | 8 | 9 | 6 | 8 | 9 | 9 | 9 | −1.3 | 0.0 |
| NABP2 | 400 | 9 | 9 | 10 | 7 | 6 | 6 | 7 | 7 | 7 | −3.0 | −2.3 |
| NAIF1 | 264 | 8 | 8 | 8 | 8 | 8 | 8 | 8 | 9 | 8 | 0.0 | +0.3 |
| ZHX3 | 261 | 7 | 7 | 6 | 7 | 7 | 7 | 8 | 8 | 8 | +0.3 | +1.3 |
| API5 | 206 | 7 | 7 | 7 | 11 | 10 | 7 | 7 | 7 | 7 | +2.3 | 0.0 |
| MYRF | 400 | 7 | 4 | 4 | 2 | 3 | 2 | 3 | 3 | 3 | −2.7 | −2.0 |
| TRIM35 | 382 | 6 | 6 | 6 | 6 | 6 | 5 | 5 | 5 | 5 | −0.3 | −1.0 |
| THAP7 | 299 | 6 | 7 | 7 | 8 | 7 | 8 | 8 | 8 | 8 | +1.0 | +1.3 |
| ZNF821 | 336 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 0.0 | 0.0 |
| RABEP2 | 254 | 5 | 6 | 5 | 6 | 5 | 6 | 5 | 5 | 5 | +0.3 | −0.3 |
| SNRK | 224 | 5 | 4 | 4 | 4 | 7 | 4 | 10 | 10 | 10 | +0.7 | +5.7 |
| POR | 260 | 5 | 6 | 7 | 5 | 4 | 5 | 4 | 4 | 4 | −1.3 | −2.0 |
| BRPF1 | 230 | 4 | 4 | 4 | 4 | 4 | 4 | 3 | 3 | 3 | 0.0 | −1.0 |
| MBD4 | 200 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 3 | 0.0 | 0.0 |
| ZNF526 | 212 | 3 | 3 | 3 | 4 | 3 | 3 | 3 | 3 | 3 | +0.3 | 0.0 |
| ZNF282 | 362 | 3 | 3 | 5 | 3 | 3 | 3 | 4 | 4 | 4 | −0.7 | +0.3 |
| ZNF251 | 219 | 3 | 3 | 3 | 2 | 3 | 2 | 3 | 3 | 3 | −0.7 | 0.0 |
| UROS | 217 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 0.0 | 0.0 |
| ZKSCAN1 | 323 | 2 | 2 | 1 | 2 | 1 | 1 | 7 | 7 | 6 | −0.3 | +5.0 |
| TFCP2 | 330 | 2 | 2 | 3 | 2 | 2 | 2 | 2 | 2 | 2 | −0.3 | −0.3 |
| GLI3 | 324 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 0.0 | 0.0 |
| ZNF530 | 264 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0.0 | 0.0 |
| ANKZF1 | 252 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0.0 | 0.0 |
| DTNBP1 | 327 | 1 | 1 | 1 | 0 | 2 | 1 | 2 | 1 | 1 | 0.0 | +0.3 |
| IFNGR2 | 295 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0.0 | 0.0 |
| SMYD2 | 205 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0.0 | 0.0 |
| RNF141 | 314 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 0.0 | 0.0 |
| TWF2 | 284 | 0 | 1 | 1 | 0 | 0 | 0 | 0 | 0 | 0 | −0.7 | −0.7 |
| CYTH2 | 203 | 0 | 2 | 0 | 0 | 0 | 0 | 2 | 2 | 2 | −0.7 | +1.3 |
| ZNF514 | 300 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.0 | 0.0 |
| SMARCC2 | 272 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.0 | 0.0 |
| ZNF543 | 211 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.0 | 0.0 |
| MED15 | 306 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0.0 | 0.0 |

**Medians and gate fractions, every scheme × seed** (50 candidates; reopening thresholds: median ≥ 20, gate ≥ 60 % = 30/50):

| scheme | seed | median `n_real` | IQR | ≥ 10 (`nmae` gate) | ≥ 20 | ≥ 30 | ≥ 89 | zero | median paired Δ vs pooled s0 | up / down / tie vs pooled s0 |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---|
| pooled | 0 (=260) | **6.0** | [2, 45] | **17 (34 %)** | 17 | 15 | 4 | 6 | 0 | 0 / 0 / 50 |
| pooled | 1 | 6.0 | [2, 43] | 18 (36 %) | 17 | 13 | 5 | 4 | 0 | 11 / 12 / 27 |
| pooled | 2 | 6.0 | [1, 46] | 18 (36 %) | 16 | 13 | 4 | 5 | 0 | 13 / 14 / 23 |
| gem | 0 | 6.0 | [2, 48] | 18 (36 %) | 17 | 13 | 4 | 7 | 0 | 14 / 12 / 24 |
| gem | 1 | 6.0 | [2, 49] | 18 (36 %) | 16 | 13 | 5 | 6 | 0 | 13 / 13 / 24 |
| gem | 2 | 5.5 | [1, 46] | 17 (34 %) | 17 | 13 | 5 | 6 | 0 | 10 / 15 / 25 |
| gem+depth | 0 | **6.5** | [2, 52] | **18 (36 %)** | 17 | 14 | 5 | 5 | 0 | 19 / 10 / 21 |
| gem+depth | 1 | 6.5 | [2, 52] | 18 (36 %) | 17 | 13 | 5 | 5 | 0 | 19 / 10 / 21 |
| gem+depth | 2 | 6.0 | [2, 52] | 18 (36 %) | 17 | 14 | 6 | 5 | 0 | 18 / 10 / 22 |

Seed-averaged per target: median of the seed-mean `n_real` **6.0 / 5.7 / 6.3**, IQR [2, 46] / [2, 48] / [2, 52];
seed-mean ≥ 10: **17 / 17 / 18 (34 / 34 / 36 %)**; ≥ 20: 17 / 16 / 17; median paired Δ vs the pooled seed-mean 0 for
both matched schemes, up/down 16/18 (gem) and 21/12 (gem+depth); paired log₂((n+1)/(n_pooled+1)) median 0.000 for
both, mean −0.063 (gem) / +0.057 (gem+depth). Uncorrected-variance sensitivity (`n_real_uc`): median 4 in every
scheme and seed.

**Reading.** The seed spread inside a scheme (pooled s1/s2 vs s0: DCP2 346 → 389/387, ANP32B 134 → 154/136,
ZNF462 50 → 63/55) is control-draw noise of about ±10 % on the large sets and ±1–2 genes on the small ones, and the
between-scheme differences are the same size and unsigned: gem − pooled runs −8.7 … +11.3, gem+depth − pooled
−60 … +35 with the two largest sets moving in opposite directions (DCP2 −60, ANP32B +35). No target crosses from the
single digits into the board's regime: the 33 candidates under the 10-gene gate in 260 read 0–11 under every
scheme (API5 hits 11/10/7 under gem, SNRK 10/10/10 under gem+depth — the only gate crossings, both by one gene), and
the 17 above it stay above it. The three deepest targets (NABP2, MYRF, TRIM35 at 382–400 cells) read 7/3/5 under
gem+depth against 9/7/6 pooled.

## 4. Diagnostics — the matching did what it was meant to, and the counts did not follow

Median over 50 targets × 3 seeds:

| scheme | ctrl-pool UMI CV | ctrl UMI mean | depth ratio median [min, max] | outside [0.95, 1.05] / 150 | zero-fraction gap (pert − ctrl) | gemgroups in pool | gemgroup TV distance target vs pool |
|---|---:|---:|---|---:|---:|---:|---:|
| pooled | 0.440 | 11,719 | 1.0135 [0.886, 1.128] | 23 | −0.0007 | 267 | 0.390 |
| gem | 0.441 | 11,742 | 1.0102 [0.893, 1.117] | 21 | −0.0012 | 167 | **0.005** |
| gem+depth | **0.378** | 11,656 | 0.9903 [0.927, 1.037] | **3** | +0.0013 | 167 | **0.005** |

- The gemgroup-matched pool has the target's gemgroup composition to within TV 0.005 (pooled: 0.39), i.e. batch
  structure is removed as far as a control pool can remove it — and its UMI CV is unchanged (0.441 vs 0.440), because
  a target's 200–400 cells span ~167 gemgroups at 1–3 cells each, so "its gemgroups" are two thirds of the screen and
  carry the screen's whole depth range. The depth heterogeneity the mandate names is therefore NOT a property of
  pooling across gemgroups; it is inside every gemgroup-proportional pool as well.
- Depth-matching within gemgroup does cut the CV (0.378) and brings every target inside the acceptance band
  (3/150 outside vs 23/150 pooled; the pooled per-target IQR is [0.996, 1.028], so 260's ratio-1.0035 reading was a
  median that hid a ±10 % per-target scatter). The zero-fraction gap, the quantity the rank test reacts to, is
  already ≈ 0 under pooling (−0.0007) and stays ≈ 0 (+0.0013). Nothing was there for the matching to repair.
- The per-target change under gem+depth is only weakly related to the pooled depth ratio (Spearman −0.28 with the
  ratio, +0.25 with the pooled zero gap, +0.25 with the pooled `n_real`): targets whose controls were slightly too
  deep gain a little, slightly too shallow lose a little — a marginal-gene reshuffle, not suppressed signal.

## 5. DE-set overlap on the five best-covered targets (+ DCP2, the largest set), seed 0

| target | cells | pooled | gem | gem+depth | ∩ all three | J(pooled, gem) | J(pooled, gem+depth) | new in gem+depth | their pooled p-rank median [max] | their pooled q median | lost from pooled | ρ(−log p) pooled vs gem+depth |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| NABP2 | 400 | 9 | 7 | 7 | 6 | 0.60 | 0.60 | 1 | 11 [11] | 0.059 | 3 | +0.83 |
| MYRF | 400 | 7 | 2 | 3 | 2 | 0.29 | 0.43 | 0 | — | — | 4 | +0.83 |
| TRIM35 | 382 | 6 | 6 | 5 | 5 | 1.00 | 0.83 | 0 | — | — | 1 | +0.90 |
| ZNF282 | 362 | 3 | 3 | 4 | 3 | 1.00 | 0.75 | 1 | 6 [6] | 0.054 | 0 | +0.91 |
| HMGXB3 | 342 | 85 | 95 | 101 | 77 | 0.82 | 0.74 | 22 | 107 [355] | 0.065 | 6 | +0.94 |
| DCP2 | 241 | 346 | 377 | 315 | 268 | 0.82 | 0.69 | 45 | 451 [1,831] | 0.069 | 76 | +0.89 |

(p-rank k = the gene was the k-th smallest pooled p-value of the 8,248 tested.) The "new" genes are the pooled
ordering's next genes at a looser bar: NABP2's one new gene was pooled rank 11 at q 0.059, ZNF282's rank 6 at
q 0.054, HMGXB3's 22 new genes sit at median rank 107 against a pooled set of 85 with median pooled q 0.065, DCP2's
45 at median rank 451 against 346 (q 0.069) with a tail out to rank 1,831 — genuinely different genes exist only in
that tail, and DCP2 simultaneously LOSES 76 of its pooled genes under the depth-matched pool. Gene rankings agree
between schemes at ρ +0.83 … +0.94 over all 8,248 genes. The deepest three targets have 2–7 genes under any scheme
(new genes for HMGXB3: SEC62, RPL34, HEXB, DCP2, HIST1H2BC, EHMT2, RPS18, GGCT, KLHDC10, MTDH, MYC, PHF19, …).

## 6. What this settles and what it does not

- **Settled**: the low DE counts in 260 §2c are not a control-pool artefact of pooling across gemgroups, nor of depth
  heterogeneity inside the pool. Removing the batch mismatch entirely (TV 0.005) and most of the depth mismatch
  (every target inside the acceptance band) leaves the median at 6 and the gate fraction at 34–36 %. Together with
  260's strong-target known-answer row (RPS3 1,205, HSPA5 425, LMO2 342 genes on the SAME pooled controls), the
  reading is unchanged: K562 answers these knockdowns with a handful of significant genes because the panel is
  non-essential in K562 by design (`plolo.py:16-18`; panel ∩ K562-essential = 0), and the scorer's recipe finds
  that faithfully.
- **Not tested (time box)**: a tighter depth match at a SMALLER control count (the greedy matching at 18,400 is
  capped by ~110 controls per gemgroup per 1–3 target cells, which is why the CV falls only to 0.378). In the
  n₂ ≫ n₁ regime the test's power is saturated in n₂ (z ∝ √(n₁n₂/N) moves 15.7 → 15.4 from 18,400 to 5,000
  controls), so a smaller pool can only add calls through a higher AUC, and the depth-matched scheme already moved
  the AUC-relevant quantity (zero-fraction gap) from −0.0007 to +0.0013 for no median gain. A pseudo-replicate /
  mixed-model DE (batch as a covariate) would be a different test from the scorer's and is out of scope: the board
  scores Wilcoxon on a pooled 18,400-control reference, exactly the pooled scheme.
- **Interpretation caveat**: the gem+depth scheme's control pool is slightly DEEPER than the targets (median ratio
  0.990), a bias of the greedy nearest-neighbour under availability caps; it is inside the band and the seed
  spread (0/1 targets change gate status across seeds) shows it does not drive the reading.

**Does batch-matched control selection change 260's coverage verdict (median ≥ 20 and ≥ 60 % clearing the
10-gene gate would reopen it)? NO** — median 5.7–6.3 and 34–36 % under every matched scheme and seed, with the
known-answer row at 50/50.

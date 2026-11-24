# 260 — Feasibility of a K562-GWPS raw-count diagnostic scored on the REAL Wilcoxon DE sets

**Date** 2026-09-05, 02:46–03:45 UTC (59 min of the 90). **Measurement only.** No harness written, no
`cell-eval2`, no build, no `vcc`, nothing under `submissions/` or `data/harness/`, no file outside this
directory modified. Scripts and outputs are in `work/` (`inspect_h5.py`, `measure_gwps.py`,
`extract_and_de.py`, `challenge_gate.py`, `summarise_de.py`, `ka_strong_targets.py`; `gwps_numbers.json`,
`panel_targets_in_gwps.csv`, `panel_targets_coverage.csv`, `candidate_targets_ge200.csv`,
`de_power_panel_targets.csv`, `de_power_curve.csv`, `de_power_summary.json`, `candidates_de_power.csv`,
`summarise_de.log`, `ka_strong_targets.csv`, `extract_rows.npz`).

**RSS.** `measure_gwps.py` 0.53 GB; `challenge_gate.py` < 1 GB; `extract_and_de.py` **peaked at 8.25 GiB**
(`ru_maxrss`, transient at the control CPM-normalisation step: two float64 temporaries of 1.2 GB on top of
the 2.2 GB extracted block), steady state 5.3 GiB during the Wilcoxon loop. **That peak breached the
mandate's ~6 GB ceiling for a few seconds**; `memory_pressure` stayed at 71 % free (the sibling 259's folds
run at 4.96–5.26 GiB, R258 §9). The fix for any re-run is float32 CPM in place — `ka_strong_targets.py`,
written that way, peaked at 3.75 GiB on the same recipe; nothing else in this experiment came near the
ceiling. Reported here rather than hidden.

---

## VERDICT LINE

**NO-GO — on DE coverage, measured, not on cost.** The raw file is on disk and streams in 1.7 min, the
depth ratio passes without repair, the K562-free model already exists (PLOLO's fold), and the whole thing
would cost ≈ 1 agent-day and ≈ $2. But the scorer's own Wilcoxon recipe run on the actual cells gives the
panel targets a **median of 5 significant genes** in K562; only **85/272 (31 %)** clear `nmae`'s 10-gene
gate, **17/50** of the best-covered candidates, and 400 cells does not change that (§2c). The mandate's
400-cell / 20 k-UMI recipe is unmeetable anyway (3 targets at 400 cells; depth 0.56× the challenge's, and
thinning only goes down). A surface whose DE sets are an order of magnitude sparser than the board's, on a
cohort a tenth the size, cannot bound what the top-89 proxy misses on the board. §5 has the full reasoning
and what survives.

---

## 1. Availability

| item | value | source |
|---|---|---|
| raw cells on disk | **yes**: `data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad`, 8,805,466,154 B — byte-equal to `scripts/scperturb_files.tsv` row 50; **no download needed** | `ls -la data/external/scperturb/`; TSV row 50 |
| storage | `X` is a **dense** gzip-4 float32 dataset (1,989,578 × 8,248), chunks (3,886 × 33), 65.6 GB uncompressed; `obs` 24 columns incl. `perturbation`, `gene`, `batch`, `ncounts`, `UMI_count`, `ngenes` | `work/inspect_h5.py`; `scripts/build_donor_table.py:76` |
| counts are raw | first 3,886-row block: all entries integral, min 0, max 1,647; row sums equal `obs.ncounts` exactly (not `UMI_count`, a different column) | `work/gwps_numbers.json` `x_integral_block`, `rowsum_eq_ncounts` |
| read rate | one 3,886-row block in 0.20 s; the **full pass measured 1.7 min** (66,869 selected rows kept, 2.55 GiB peak) | `gwps_numbers.json` `block_read_s`; `work/extract_and_de.log` |
| the compact donor | `data/donors/k562_gwps/delta.npy` (9,787 × 8,248) is **log-pseudobulk effect vs control** (`b_p − b_ctrl`, TS 5e4), not cells; 75,328 control cells; ≥ 20 cells per group; median 179 cells/group | `data/donors/k562_gwps/meta.json`; `scripts/build_donor_table.py:1-22` |
| depth | median UMI/cell **11,304** (control 11,372, perturbed 11,301); percentiles 5/25/50/75/95 = 4,032 / 7,994 / 11,304 / 14,684 / 20,823; median 3,215 genes detected/cell | `gwps_numbers.json` |
| challenge depth | context A controls: 18,400 cells, median UMI **20,109** | `work/challenge_gate.py` on `data/challenge/context_A.h5ad` |
| gene coverage | 8,248 genes, **7,679 on the official 18,533 axis (41.4 %)**; 231/300 panel targets are among the 8,248 (the other 69 targets' own-gene knockdown is invisible; irrelevant to scoring, which drops the target gene) | `gwps_numbers.json`; `build_donor_table.py:18` |
| DE gate regime | K562: **8,248 / 8,248** genes pass the scorer's control-mean-CPM > 5 gate (the axis is pre-filtered to expressed genes, so CPM is inflated ~2.2×); challenge context A: **9,929 / 18,533** pass | `work/extract_and_de.log`; `work/challenge_gate.py` |
| batch structure | 267 gemgroups, 860–10,042 cells each; per-gemgroup median UMI 2,854–22,130; control fraction flat **3.1–4.8 %** in every gemgroup (pooled screen: every gemgroup carries controls and all guides) | `gwps_numbers.json` |
| depth ratio | panel targets' median UMI / control median: **1.0035**, IQR [0.969, 1.028] — inside `check_harness_depth.py`'s PASS band [0.95, 1.05] with NO repair (Orion was 0.725, `docs/HARNESS-DEPTH-CONFOUND.md`) | `gwps_numbers.json` `panel_depth_ratio_median_iqr`; `scripts/check_harness_depth.py:26` |
| chemistry | GWPS is 10x 3′ (Replogle 2022, outside the repo); the 2025 H1 / 2026 references are 10x Flex (`scripts/make_harness.py:5`) — a regime difference no downsampling removes | as cited |

**What the 20 k / 400 recipe does here.** `make_harness.py` thins at `rate = target_umi / median` and returns
the matrix unchanged when `rate ≥ 1` (`thin_counts`, `:52-60`); at 11,304 the rate is 1.77, so the harness
would sit at the file's native depth, 0.56× the challenge's. It drops any perturbation with < 400 cells
(`:91-92`); §2 shows that leaves 3 targets.

## 2. Panel overlap, cell counts, and DE power

### 2a. Overlap and cell counts (`work/measure_gwps.py`, `work/panel_targets_in_gwps.csv`)

- **272 / 300** panel targets are in GWPS (missing 28: ABCD1, CAPRIN2, EPHB2, HEATR6, KIF21B, MAPK7, MLKL,
  NICN1, PARP3, PBLD, PDK3, PHF11, PSMB9, RAB11FIP5, SEMA4F, SLC44A1, SNN, SOCS5, SPATA20, TAF4, TBC1D19,
  TDRD7, TESK2, TMEM104, TMEM79, TTBK2, TXNDC16, ZC2HC1A); the donor carries 269 (≥ 20 cells,
  `make_submission_pm_cd4.py:5835`).
- Per-target cells on the 272: **min 13, Q1 126, median 171.5, Q3 225, max 1,745.**

| cells ≥ | 100 | 150 | 200 | 250 | 300 | **400** |
|---|---:|---:|---:|---:|---:|---:|
| panel targets | 223 | 168 | 95 | 55 | 23 | **3** |

- Whole file: median 178 cells/group; 6,086 / 9,866 groups ≥ 150; 459 ≥ 400. GPT's prior ("many ~150, few
  400") is confirmed on the panel: 400 is met by 3 targets, 150 by 168.
- `k562_essential` (`data/donors/k562_essential/perts.txt`, 2,021 names) ∩ panel = **0** — the essential
  screen cannot serve; only GWPS can.

### 2b. Coverage of the same targets by every other production input (`work/panel_targets_coverage.csv`)

Perts files of each registry bank (`make_submission_pm_cd4.py:840-1000` paths; `data/donors/orion_*`;
`data/external/perturbmatch`), CD4-complete from `experiments/249-cd4-deflate/work/stage/cd4-off.npz`
(`target_condition_count == 3`, 276 targets):

| K562 cells ≥ | n | both Orion | six-source union | H1 pair | PerturbMatch | CD4-complete | union ∧ Orion ∧ CD4 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 150 | 168 | 161 | 100 | 35 | 154 | 159 | **91** |
| 200 | 95 | 90 | 56 | 24 | 88 | 90 | **50** |
| 250 | 55 | 55 | 36 | 16 | 50 | 53 | **35** |
| 300 | 23 | 23 | 14 | 3 | 20 | 23 | 14 |

(Panel-wide: six-source union 176/300, H1 pair 77, `dld1_lowmoi` 67, `igvf_cardiomyocyte` 87,
`gse295214_h1` 64, `h1_2025` 20, `myllia_a549` 14, `jiang_5line_noK562` 9, PerturbMatch 251, Orion 290/294.)

**A 32–40-target diagnostic** therefore draws from the **50 targets with ≥ 200 K562 cells that are carried by
the six-source union, both Orion donors and all three CD4 conditions** (`work/candidate_targets_ge200.csv`:
NABP2 664, MYRF 450, TRIM35 382, ZNF282 362, HMGXB3 342, ZNF821 336, TFCP2 330, DTNBP1 327, GLI3 324,
ZKSCAN1 323, RNF141 314, MED15 306, ZNF462 300, ZNF514 300, THAP7 299, IFNGR2 295, FOXO4 293, TWF2 284,
GPATCH8 284, ZNF507 279, SMARCC2 272, MED13 268, DNTTIP1 267, ZNF530 264, NAIF1 264, PHF23 261, ZHX3 261,
ANP32B 260, PCGF1 260, POR 260, EPC2 259, PPIF 256, RABEP2 254, ANKZF1 252, GIGYF2 251, DCP2 241, FOXN2 232,
BRPF1 230, E2F1 226, SNRK 224, FCGRT 224, ZNF251 219, UROS 217, PAN2 215, ZNF526 212, ZNF543 211, API5 206,
SMYD2 205, CYTH2 203, MBD4 200), at **200–400 cells each (cap 400)**, the file's native depth, 18,400
controls drawn uniformly over the 267 gemgroups (which is also how the perturbed cells are distributed, so
the depth match in §1 is automatic). At ≥ 150 cells the pool is 91.

Which union sources carry the 50 (perts files as above): `igvf_cardiomyocyte` 27, `dld1_lowmoi` 20,
`gse295214_h1` 16, `h1_2025` 4, `jiang_5line_noK562` 3, `myllia_a549` 3; 30 are carried by one source, 18 by
two, 2 by three or more. The mixture is close to the panel-wide one (87 / 67 / 64 / 20 / 9 / 14 of 300), so
the diagnostic probes the same source blend ARM6 shipped, with `igvf_cardiomyocyte` — the source R258 §6
names as the carrier of the union's `nmae` cost — the largest contributor on both.

### 2c. DE power: the scorer's Wilcoxon sets on these cells (`work/extract_and_de.py`)

**Method** (the scorer's, `docs/reference/vcc2026-metrics-BRIEF.md` §3; convention checked in the installed
scorer, `cell_eval2/de_compute.py:871-872`: scanpy `rank_genes_groups(method="wilcoxon", tie_correct=True)`
against the control group): all 272 panel targets present, up to 400 cells each (seed 0) against 18,400
control cells drawn uniformly (the same draw a harness would use), per-cell CPM to 1e6, genes with control
mean CPM > 5 (all 8,248 here), two-sided Mann-Whitney U with the normal approximation, tie-corrected (`tc`,
the scorer's convention) and uncorrected (`uc`, for sensitivity), Benjamini–Hochberg per target, α = 0.05,
target gene excluded; `gate` = significant genes with a finite non-zero LFC, i.e. `nmae`'s cohort rule
(≥ 10, BRIEF §7). Sign of the LFC uses the scorer's arithmetic-mean CPM ratio. No thinning (§1).

Comparators from real references at 400 cells / 20 k UMI: h1-test's median reference DE-set size 89
(R188, the number `POOL_N` is set to); the real scorer's `nmae` cohort on h1-test is **58 of 85**
perturbations (`experiments/221-cloud-scorer/cloud-output/ladder6-20260904/a15/run/agg_results.csv`,
`count` of `de_wilcoxon_lfc_nmae`), and on the challenge 209–261 of 300 (BRIEF §8).

**Result: the panel targets barely move K562's transcriptome.** 272 targets, 33.4 min, control median UMI 11,406
(`work/de_power_panel_targets.csv`, `work/summarise_de.log`):

| K562 cells | targets | `n_real` median (tie-corrected) | IQR | clear `nmae`'s ≥ 10 gate | uncorrected median |
|---|---:|---:|---|---:|---:|
| < 100 | 49 | 2 | [0, 10] | 13 (27 %) | 0 |
| 100–149 | 55 | 4 | [1, 9.5] | 14 (25 %) | 2 |
| 150–199 | 73 | 6 | [1, 18] | 28 (38 %) | 5 |
| 200–249 | 40 | 3.5 | [1.75, 18] | 15 (38 %) | 2 |
| 250–299 | 32 | 7.5 | [1, 26] | 13 (41 %) | 5.5 |
| 300–400 | 23 | 3 | [1, 6] | 2 (9 %) | 2 |
| **all 272** | 272 | **5** | — | **85 (31 %)** | — |

- 58 of 272 targets have **zero** significant genes; the fraction of significant calls that are negative is
  0.585 (median over targets) — a clean, depth-balanced regime (cf. 0.910 on the confounded Orion file,
  0.538 on h1-test, `docs/HARNESS-DEPTH-CONFOUND.md`), so the small sets are not a depth artefact.
- **The knockdowns themselves are real and detected**: for the 229 targets whose own gene is on the axis,
  the target gene's LFC is median **−3.05**, negative on 100 % and below −0.5 on 99 %. The CRISPRi worked;
  the downstream response in K562 is what is small. That is the panel's design (a deliberately
  non-common-essential set, `plolo.py:16-18`; panel ∩ K562-essential = 0), not a measurement problem.
- **More cells do not rescue it** (`work/de_power_curve.csv`, the three deepest targets subsampled at seed 1):

  | cells | 100 | 150 | 200 | 300 | 400 |
  |---|---:|---:|---:|---:|---:|
  | MYRF | 1 | 1 | 3 | 4 | 7 |
  | NABP2 | 1 | 1 | 1 | 3 | 9 |
  | PINK1 | 3 | 2 | 2 | 5 | 6 |

  At the mandate's 400 cells the two deepest candidates are still under the 10-gene gate.
- **On the 50 candidate targets** (`work/candidates_de_power.csv`): `n_real` median **6**, IQR [2, 45], 0–346;
  **17 / 50 clear the gate**, 15 reach 30, 4 reach h1-test's median of 89 (DCP2 346, ANP32B 134, DNTTIP1 113,
  FCGRT 103; then HMGXB3 85, PHF23 80, EPC2 77, MED13 65, PCGF1 64, FOXN2 63, PAN2 59, GIGYF2 54, ZNF462 50,
  E2F1 31, FOXO4 30, ZNF507 25, PPIF 24). 33 candidates would be dropped from `nmae` by the scorer's own
  cohort rule, among them the deepest (NABP2 400 cells → 9, MYRF 400 → 7, TRIM35 382 → 6, ZNF282 362 → 3).
  At ≥ 150 cells the fully-covered pool of 91 clears the gate on 34.
- Against the board's regime: h1-test's real cohort is 58/85 = 68 % with a median set of 89 genes at 400
  cells / 20 k UMI; the challenge's is 70–87 % (209–261 of 300). K562-GWPS on this panel is **31 %** with a
  median set of **5**, and 34 % (17/50) on the best-covered candidates — a different DE regime by an order of
  magnitude in set size and by half in cohort.

**Known-answer row — the same recipe on strong K562 knockdowns** (`work/ka_strong_targets.py`, same
18,400 control rows, same tie-corrected test, 3.75 GiB peak; `work/ka_strong_targets.csv`). If the panel's
small sets were a pipeline defect, these would be small too. They are not:

| target | cells | `n_real` | at 150 cells | target-gene LFC |
|---|---:|---:|---:|---:|
| RPS3 | 163 | **1,205** | 1,013 | −1.00 |
| HSPA5 | 299 | **425** | 72 | −0.78 |
| LMO2 | 273 | **342** | 104 | −0.67 |
| POLR2A | 231 | **204** | 164 | −0.97 |
| TAL1 | 400 | **136** | 26 | −2.67 |
| BRD4 | 239 | 20 | 6 | −0.66 |
| RPL3 | 400 | 18 | 5 | −0.37 (weak guide) |
| KLF1 | 157 | 14 | 11 | −2.54 |

Hundreds to a thousand significant genes where K562 biology responds, at the same depth and cell counts
the panel targets have; the panel's median of 5 is therefore the panel's biology in K562, and the
diagnostic's DE sets would be sparse for that reason alone.

## 3. The leak-audit surface, enumerated

Every production input that derives from Replogle K562-GWPS, with the file that proves it:

| # | input | where | K562 content | removable? |
|---|---|---|---|---|
| 1 | **`k562_gwps` donor** | `PRODUCTION_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")`, equal weights (`make_submission_pm_cd4.py:461-462, 1342-1352`; `scripts/def_pipeline.py:120`); frozen shas `:116-131` | one third of the donor mixture on the 269 targets it carries; the pseudobulk of the SAME cells the diagnostic would use as truth — a pure oracle leak on every one of them | yes: PLOLO's `k562_gwps` fold |
| 2 | **PerturbMatch** | `data/external/perturbmatch/meta.json`: "Replogle et al. (2022) K562 genome-wide CRISPRi Perturb-seq reanalysis by Yeung et al. (2026)", limma-voom LFC, 9,471 perts × 9,761 genes, 251 panel targets; blended 50:50 into the K562 slot on 248 targets (`PM_BLEND` `:606`, `EXPECTED_PM_REPLACEMENT_TARGETS = 248` `:1267`, `apply_perturbmatch_k562_replacement` `:1403-1525`) | the same raw K562 cells re-analysed: a second oracle leak on 248 targets; its term is multiplied by `k562_share` (`:1476-1495`), so it vanishes with the slot | yes, with #1 |
| 3 | **the whitening operator** `data/donors/whitener_union_pruned.npz` `e2f65ecf…` | `whitener_union_pruned.json`: donors `k562_gwps` (7,784 basis rows kept, 2,003 pruned), `orion_hct116` 10,807, `orion_hek293t` 10,803; `dropped_in_panel: []`, so the panel targets' own K562 rows ARE basis rows | second-order (covariance) exposure to the truth rows; diffuse but real | yes: 062's `held_k562_gwps` (`experiments/062-cd4-post-whiten/cloud-output/whiteners/held_k562_gwps.npz`, `13d988fc…`, 17,806 columns, Orion-only; its own json: basis rows `orion_hct116` 17,767 + `orion_hek293t` 17,977, no 12,000-row cap, no `basis_prune` key — **unpruned and uncapped**, `plolo.py:177-186`) |
| 4 | **essential prune list** `data/donors/k562_essential/perts.txt` | `basis_prune: essential`, 2,021 names (`whitener_union_pruned.json`) | a gene NAME list from the K562 essential screen, no response data; ∩ panel = 0 | not a leak of the truth; note the held operator predates pruning (#3) |
| 5 | `k562_essential` union-registry entry | `make_submission_pm_cd4.py:982-990`, "registered, never selected"; `VCC_UNION_SOURCES` unset selects the H1 pair (`:993-998`) | inert in every shipped build; 0/300 panel coverage | n/a |
| 6 | donor weights | `VCC_DONOR_WEIGHTS` unset = equal (`:549-582`); the rationale was measured with K562's sign accuracy on h1-test's DE pool (175) | development shaping, not an input | n/a |
| 7 | CD4 residual, H1-union residual, `V5_abc` basis | directions are CD4 / H1 / DLD-1 etc. (`cd4_deflation_basis/V5_abc.json` provenance: fitted on the staged CD4 residual); their solved fractions are sized against the whitened base's row norms (`solve_carrying_row_fraction` `:1897`; `apply_frozen_cd4_residual` `:2493`) | scalar dependence on the base only | re-solved on the K562-free base, as PLOLO already does |
| 8 | emission (`rc_callset`, spread, `--min-cpm 10` on the TARGET context's controls) | `HARNESS-DEPTH-CONFOUND.md` caveat 1; `build_h1_candidates.py:78-84` | none | n/a |
| 9 | stray artifacts | `data/donors/k562_gwps/svd_*.npy`, `transfer_w_*.npy` (Aug 21), `data/donors/194-debias-k562_gwps` | not referenced by the builder, driver or `src/` (grep) | n/a |
| 10 | **development shaping** | PLOLO's `k562_gwps` fold is one of the four folds behind every `pds` price (7/7 board signs, READ-FIRST §4); the identity's `nmae` pools all four folds (`E237/ADDENDUM-2-READING.md` §2, `six229.json`); 188, 175, 203, 216 all read K562 pseudobulk | K562 pseudobulk has selected every shipped knob — it is a development surface, not an untouched holdout | not removable; bounded by preregistering ONE comparison whose knobs are already frozen |

**"A model with every K562-derived input removed" already exists and is named:** PLOLO's `k562_gwps` fold
(`plolo.py:131-133`): sources `orion_hct116 + orion_hek293t`, operator 062 `held_k562_gwps`, **no PerturbMatch**
(`build_fold` blends PerturbMatch only when `k562_gwps` is a source, `plolo.py:555-560`), CD4 and union
residuals re-solved on that base; cached at `experiments/162-retrodicting-validator/work/fold-k562_gwps.npz`
(`pre` 269 × 18,533, `truth` 269 × 18,533, `wcols` 17,806, `targets` 269).

**That is a different model from production, plainly**: it drops one of three donors and the 248-target
PerturbMatch blend, and whitens with an Orion-only, unpruned, 17,806-column operator instead of the
three-donor pruned 17,942-column one. **So the diagnostic calibrates the INSTRUMENT — the identity's proxy
`nmae*` (top-89 pools, pseudobulk truth) against the real scorer's `de_wilcoxon_lfc_nmae` (Wilcoxon-BH sets on
400-cell-scale truth) on the SAME fold, targets and arm pair — not the shipped prediction.** That is the
right object for the mandate's question: the identity IS computed on this fold (`six229.py` pools four folds,
K562 included), so any DE penalty the proxy misses would be visible exactly there. What it cannot say is
what the shipped, K562-carrying prediction scores on K562 — that number is oracle-leaked by #1–#3 and must
never be read.

## 4. Cost

| step | what | wall clock | machine | evidence |
|---|---|---|---|---|
| download | none — file on disk, byte-exact | 0 | — | §1 |
| extraction | one streaming pass of the dense gzip `X`, keeping ≤ 400 cells × 272 targets + 18,400 uniform controls (seed 0; the row indices are saved in `work/extract_rows.npz`) | **1.7 min measured** | 2.6 GiB | `work/extract_and_de.log` |
| harness assembly | relabel `control → non-targeting`, `obs[target_gene]`, subset columns to the 7,679 official-axis genes (the prediction must carry the identical axis; 221 already emits on h1-test's 18,080-gene axis), **no thinning** (rate ≥ 1), 32–40 targets × 200–400 cells + 18,400 controls ≈ 30 k cells, csr gzip h5ad; then `check_harness_depth.py` (expected PASS, ratio ≈ 1.00, §1) | ≈ 10 min | < 3 GiB | `make_harness.py:120-135`; `check_harness_depth.py` |
| K562-free effect for TWO arms | PLOLO fold arms: **38 s / 4.96 GiB for 38 arms on the k562 fold** (R258 §9); the pair is two of them | minutes | 5 GiB | `experiments/258-*/RESULTS.md` §9 |
| rendering to counts | `profile_from_effect` + `emit_counts` on the harness's own controls and library sizes, exactly 221's path (`build_h1_candidates.py`: five 85-target × 400-cell candidates in **19:51** on the worker, ≈ 4 min each); NEW adapter from the fold objects (`pre`, `b_ctrl`, `gate`, residual fractions) to the harness axis, with two known-answer rows: the fold's `pds` for `h1union` must reproduce **0.8087859956721968** (`six229.json` `pds_per_fold/k562_gwps`) and the emitted pseudobulk must reproduce PLOLO's `deliver()` profile | **3–5 h of agent work**, ≈ 10 min compute | ≤ 6 GiB if streamed | `cloud-output/ladder5-20260904/prologue.log:179`; `plolo.py:262` |
| cloud scoring | h1-test reference points: `cell-eval2 run` **2:32–2:41 h** per candidate (85 perts, 52,400 cells, 18,080 genes; five logs), `baseline` 2.6 h, `prep-real-bundle` ≥ 4 h (R221 §7). DE work scales with perts × cells × genes: 40 × ~30 k × 7,679 / (85 × 52.4 k × 18,080) ≈ **0.11×** → `run` ≈ 20 min each, `baseline` ≈ 20 min, bundle ≈ 30–45 min | ≈ 2 h wall incl. launch, mostly waiting | r7i.2xlarge, on-demand for the bundle (R221: no spot on long steps) | `cloud-output/*/run.log` "Elapsed"; `cloud_score_worker.sh:140-147` |
| money | 221: $1.48 per six-member score, $1.37 of it the 2.5 h `run`; at 0.11× the workload, **≈ $2 for the bundle + two scores** | — | — | R221 §6 |
| reading | raw per-target `de_wilcoxon_lfc_nmae` of both arms from the run's per-perturbation table vs the identity's per-target `nmae*` on the same fold (`216/six216.py` pools, `wald89`; 38 s to re-run per fold) — compare RAW deltas, as 221 does (the scaled span `r − b` of a K562 harness is not the board's) | 1–2 h | — | `six229.py`; R221 §7 |

**Total: ≈ 1 agent-day (≈ 5–8 h of work, ≈ 3 h of waiting), ≈ $2, zero submission slots.** No step is a
major project. The one piece of new code is the fold→harness rendering adapter.

**Where the comparison lands.** The frozen pair is ARM6's variable — the six-source union at g = 0.35 on
top of the `h1union` base (`VCC_UNION_SOURCES` = the six, `VCC_H1_UNION_G` = 0.35 vs the pair at 0.28) —
because it is the ONLY arm class with a board reading of the identity: predicted −0.004783, realised
−0.00716 overall (**1.50×**, BR237 §2); its raw four-fold Δ`nmae*` is +0.01722 (`six229.json`
`stack6 g=0.35/d_nmae_vs_incumbent`). The diagnostic reads, on the same K562 fold and the same 32–40 targets:
(i) the identity's raw Δ`nmae*` per target, (ii) the scorer's raw Δ`de_wilcoxon_lfc_nmae` per target, and
asks whether (ii)/(i) is ≈ 1.5 (the board's ratio: the proxy is a conservative bound because the real sets
are larger/differently composed) or materially larger, and on which targets and genes the excess sits.
The **known-answer row** is the null pair `cd4svd5` vs `h1union` (`VCC_CD4_DEFLATE_K=5`), whose four DE
members were flat as priced on the board (`nmae` scaled +0.11047 vs +0.11074, BR249): the diagnostic must
read ≈ 0 on it, else the surface is mis-wired.

## 5. Verdict

**NO-GO.** The step that kills it is neither extraction (1.7 min, measured), nor the exclusion audit (the
K562-free model exists and is cached: PLOLO's fold), nor the depth ratio (1.0035, PASS without repair), nor
the money (≈ $2). It is **DE coverage — the biology of the surface**: the challenge panel is, by design,
non-essential in K562, and K562-GWPS answers its knockdowns with a median of **5** Wilcoxon-significant
genes per target; **31 % of the panel (85/272) clears `nmae`'s 10-gene gate, 34 % of the best-covered
candidates (17/50)**, and the deepest targets (400 cells) still read 6–9 genes. Cells cannot fix it (the
subsampling curve is flat), depth cannot fix it (11.3 k UMI is the file's ceiling, and depth is not the
limiting factor at 400 cells either), and selecting the responsive 17 makes it a 17-target cohort with
sets chosen on the truth — the board's `nmae` is read on 209–261 of 300 targets with sets of ~89 genes.
A Δ`nmae` of ARM6's size (raw +0.026 on the board, BR237) on 17 targets would sit at about one standard
error of the per-target scatter, and whatever the top-89 proxy misses on 5-gene sets is not what it misses
on 89-gene sets. **The diagnostic would measure a different DE regime, not the board's, on a model that is
not production's (§3), with a cohort too small to resolve the effect it is meant to bound.** That is a
well-supported negative: every number above is measured on the actual cells with the scorer's own recipe,
and the pipeline is validated in two directions (the target-gene knockdown is detected on 99 % of targets;
the strong-target known-answer row in §2c-KA).

What survives, for the record:

1. **The extraction and the leak audit are cheap and done** (`work/extract_rows.npz` reproduces the seed-0
   draw; §3 names every K562-derived input and the K562-free fold). If a K562 surface is ever wanted for a
   NON-DE question — `pds`/`mse` at count level, where PLOLO's pseudobulk fold already serves — it is a
   half-day, not a project.
2. **The mandate's recipe (400 cells / 20 k UMI) is unmeetable here on two counts**: 3 targets at 400 cells,
   and 0.56× the depth with thinning only able to go down.
3. **The real-DE-set question stays open and the only surface that answers it in the board's regime is the
   board itself** (BR237: the identity's leg at 150 %, sign right). If a local real-DE surface is still
   wanted, the candidates are H1-like references at Flex depth with dense responses — h1-test with its
   85-target oracle leak made explicit and priced (the `h1_2025` donor removed from the arm under test,
   which is exactly how PLOLO treats a held library), not a K562 surface whose DE sets are an order of
   magnitude sparser. That is a separate feasibility question, not started here.
4. If someone nonetheless wants the K562 number, the frozen comparison and its known-answer row are written
   in §4 and cost ≈ 1 agent-day and ≈ $2; the 17-target list is in `work/candidates_de_power.csv`. Read it as
   a 17-target anecdote, never as a calibration of the identity.

**GO / NO-GO: NO-GO** (coverage: K562-GWPS's real Wilcoxon sets on the panel are too sparse — median 5
genes, 31 % cohort — to stand in for the board's DE regime; not fixable by cells, depth or money).

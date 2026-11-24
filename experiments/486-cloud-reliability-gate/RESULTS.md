# 486 — Does a learned per-(target, gene) reliability gate transfer?

**Answer.** The learned gate does not transfer. On energy-weighted cosine — the pre-registered comparison, "beats" meaning the target-bootstrap 95% interval of the difference lies entirely above 0 — it is worse than the equal-weight mean on **4/4** held-out lines and worse than Testh's fixed formula on **3/4**. RPE1 is a tie against the fixed formula (+0.0001, interval [−0.0012, +0.0012]).

The fixed formula does transfer on the typical target and on discrimination. Unweighted cosine, normalised rank, and top-1 each beat the equal-weight mean on **4/4** lines. Energy-weighted cosine for the fixed formula beats equal-weight on **1/4** lines (K562, +0.0154). It loses on RPE1 (−0.0050) and Jurkat (−0.0062) and is tied on HepG2 (+0.0001, interval crosses 0).

Decision: on these four public lines there is transferable signal in the fixed reliability weight for rank and for the average target. Fitting the eight-feature sigmoid on the other lines does not raise energy-weighted cosine on a line it was not fit on. Do not treat Testh's learned-gate board move as something this specification reproduces.

| Pre-registered contrast (energy-weighted cosine) | Lines where the 95% CI is entirely above 0 |
|---|---:|
| learned − equal | **0 / 4** |
| learned − fixed | **0 / 4** (RPE1 tie; the other three lose) |
| fixed − equal | **1 / 4** (K562 only) |

## Question

Testh's board history says a fixed reliability weight raised pds 0.8054 → 0.8203, and a learned gate plus a per-target weight reached pds 0.8356 (mse 0.810 → 0.699). This experiment asks a narrower question on public data, with no production pipeline and no held-out-line leakage: **does a learned per-(target, gene) gate beat the equal-weight mean and Testh's fixed formula when the predicted cell line was not used to fit it, and on how many held-out lines?**

The rank below is a raw normalised rank on ~2,393 shared targets. It is not scaled board `pds`. The sign of the fixed-formula rank gain agrees with their disclosure. The magnitude is not a prediction of a board delta.

## Point estimates (CONFIRMED)

Leave-one-line-out. `ew` is energy-weighted cosine, `uw` unweighted cosine, `rank` mean normalised rank, `top1` fraction of targets whose own truth is strictly nearest. Source: `work/results/scores_by_line.csv`, from `uv run python evaluate.py` (log `work/evaluate.log`, `EVAL_EXIT:0`).

| Held out | n | equal ew | fixed ew | learned ew | equal uw | fixed uw | learned uw | equal rank | fixed rank | learned rank | equal top1 | fixed top1 | learned top1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| K562 | 2394 | 0.2170 | **0.2324** | 0.1903 | 0.1762 | **0.1891** | 0.1625 | 0.8447 | **0.8595** | 0.8253 | 0.0952 | **0.1216** | 0.0739 |
| RPE1 | 2393 | **0.2790** | 0.2740 | 0.2741 | 0.2279 | 0.2313 | **0.2326** | 0.8130 | **0.8269** | 0.8231 | 0.1408 | **0.1697** | 0.1592 |
| HepG2 | 2393 | 0.2715 | **0.2716** | 0.2551 | 0.2292 | **0.2330** | 0.2199 | 0.7775 | 0.7866 | **0.7926** | 0.0618 | **0.0731** | 0.0719 |
| Jurkat | 2393 | **0.1837** | 0.1775 | 0.1687 | 0.1589 | **0.1666** | 0.1536 | 0.8213 | **0.8311** | 0.8074 | 0.0832 | **0.1074** | 0.0694 |

Gene-permuted refit matches equal-weight to about 0.0001 on every line (K562 ew +0.00006, RPE1 −0.00010, HepG2 −0.00000, Jurkat −0.00002). Target shuffle rank is 0.500, 0.513, 0.506, 0.501 and top-1 is 0.0004 on every line (one target). Shuffle cosine is not zero (ew 0.039, 0.128, 0.106, 0.069) because responses share a component; it sits far below real transfer. **CONFIRMED:** the pipeline is not leaking a constant match, and a gate whose features are scrambled across genes finds nothing.

## Bootstrap (CONFIRMED)

4,000 paired target resamples, seed 486, percentile 95% interval. Full rows: `work/results/gains.csv`.

Energy-weighted cosine:

| Held out | learned − equal | learned − fixed | fixed − equal |
|---|---|---|---|
| K562 | −0.0267 [−0.0299, −0.0236] loses | −0.0421 [−0.0453, −0.0390] loses | **+0.0154 [+0.0142, +0.0167] beats** |
| RPE1 | −0.0050 [−0.0067, −0.0033] loses | +0.0001 [−0.0012, +0.0012] tie | −0.0050 [−0.0060, −0.0040] loses |
| HepG2 | −0.0164 [−0.0181, −0.0147] loses | −0.0165 [−0.0184, −0.0146] loses | +0.0001 [−0.0014, +0.0015] tie |
| Jurkat | −0.0150 [−0.0206, −0.0102] loses | −0.0088 [−0.0123, −0.0053] loses | −0.0062 [−0.0107, −0.0016] loses |

Fixed minus equal on the other three metrics. The interval is entirely above 0 on every line:

| Held out | uw | rank | top1 |
|---|---|---|---|
| K562 | +0.0129 [+0.0122, +0.0136] | +0.0148 [+0.0132, +0.0165] | +0.0263 [+0.0196, +0.0338] |
| RPE1 | +0.0033 [+0.0026, +0.0040] | +0.0139 [+0.0121, +0.0159] | +0.0288 [+0.0217, +0.0359] |
| HepG2 | +0.0038 [+0.0033, +0.0043] | +0.0091 [+0.0074, +0.0109] | +0.0113 [+0.0059, +0.0171] |
| Jurkat | +0.0077 [+0.0071, +0.0084] | +0.0097 [+0.0080, +0.0115] | +0.0242 [+0.0180, +0.0309] |

Learned minus equal is not uniformly negative off the pre-registered metric. Unweighted cosine beats equal-weight only on RPE1 (+0.0046). Rank and top-1 beat equal-weight on RPE1 (+0.0102, +0.0184) and HepG2 (+0.0150, +0.0100) and lose on K562 and Jurkat. Learned minus fixed on rank beats only HepG2 (+0.0059) and loses on the other three. Those rank gains do not change the pre-registered call: the fit was for energy-weighted cosine, and on that metric the learned gate loses.

## Why energy-weighted cosine and rank disagree (post-hoc, CONFIRMED as a description)

Not pre-registered. Computed after the scores, from `work/results/per_target.npz`, by `uv run python energy_quartiles.py`. Quartiles are of truth energy ‖y‖². The cell is the mean per-target cosine gain.

| Held out | fixed−equal Q1 (low) | Q4 (high) | learned−equal Q1 | Q4 | corr(fixed gain, energy) |
|---|---:|---:|---:|---:|---:|
| K562 | +0.0106 | **+0.0161** | +0.0024 | −0.0419 | +0.13 |
| RPE1 | +0.0121 | **−0.0065** | +0.0093 | −0.0092 | −0.30 |
| HepG2 | +0.0062 | +0.0014 | +0.0005 | −0.0226 | −0.18 |
| Jurkat | +0.0112 | +0.0032 | +0.0021 | −0.0133 | −0.25 |

The fixed gate lifts the low-energy quartile on every line. Energy-weighted cosine is dominated by the largest responses. On K562 those also improve, so `ew` rises. On RPE1 they get worse, so `ew` falls while unweighted cosine, rank, and top-1 rise. The learned gate's damage grows with energy on every line.

## Coefficients (CONFIRMED)

Nelder–Mead, three starts (zero plus two N(0, 0.5) draws), 500 evaluations each, 1,500 total. Every restart hit the evaluation cap (`success: false` in the fold JSON). That is the cap, not a crash. The three starts agree to about 0.001 on the training objective (HepG2's third start is 0.1820 against a selected 0.1832). Float32 training objective matches the float64 rescore to ~1e-8 on every fold, inside the 0.001 abort threshold. Source: `work/results/coefficients.csv`.

| Term | K562 held out | RPE1 held out | HepG2 held out | Jurkat held out |
|---|---:|---:|---:|---:|
| intercept | −0.004 | −1.350 | −1.708 | +0.930 |
| log_z2 | −0.455 | **+0.483** | **+0.448** | +0.472 |
| log_n_src | −2.771 | +0.332 | +0.321 | **−21.05** |
| log_ctrl | +3.182 | +0.860 | +0.613 | +2.663 |
| log_abs_mu | +0.263 | +0.093 | +0.310 | −0.200 |
| logit_sign_agree | −0.144 | **+0.580** | **+0.742** | −0.098 |
| log_cross_target_var | +2.051 | −0.080 | −0.339 | +1.736 |
| log_n_cells | −0.005 | +0.662 | +0.670 | +0.250 |
| log_source_mad | +0.004 | **−0.406** | **−0.486** | −0.209 |
| train ew (selected) | 0.2608 | 0.1918 | 0.1832 | 0.2659 |
| perm-refit train ew | 0.2391 | 0.1782 | 0.1683 | 0.2448 |
| outer learned ew | 0.1903 | 0.2741 | 0.2551 | 0.1687 |

When the training lines include K562 (RPE1 and HepG2 held out), the signs are a reliability gate: up-weight precise entries, sign agreement, more cells, and control expression; down-weight source MAD. Those two folds are also the smallest learned losses. Holding out K562, or holding out Jurkat, blows coefficients up (`log_n_src` −2.77 and −21, `log_ctrl` +3.18 and +2.66). A coefficient of −21 on `log_n_src` means the z-scored feature barely varies and the fit is amplifying a residual. The training objective still rises versus the gene-permuted refit on every fold (+0.022, +0.014, +0.015, +0.021), and the same coefficients then lose to equal-weight out of sample on energy-weighted cosine. That is overfitting the three training lines, worst where the linear predictor is unstable.

## Which features the fitted gate uses (CONFIRMED)

Each coefficient set to 0, no refit. `drop_vs_learned` is the fall in outer energy-weighted cosine. Positive means that feature was helping the learned prediction. The learned prediction is already worse than equal-weight, so a positive drop is not a transfer win. Source: `work/results/ablation.csv`.

| Feature | K562 | RPE1 | HepG2 | Jurkat | Mean drop |
|---|---:|---:|---:|---:|---:|
| log_ctrl | +0.0206 | +0.0138 | +0.0039 | +0.0027 | **+0.0102** |
| log_abs_mu | +0.0029 | +0.0005 | +0.0015 | −0.0024 | +0.0007 |
| log_n_cells | −0.0000 | +0.0016 | +0.0020 | −0.0016 | +0.0005 |
| log_z2 | −0.0107 | −0.0007 | +0.0006 | +0.0025 | −0.0021 |
| log_n_src | −0.0013 | +0.0000 | +0.0004 | −0.0078 | −0.0022 |
| log_source_mad | −0.0000 | −0.0068 | −0.0054 | −0.0009 | −0.0033 |
| logit_sign_agree | −0.0033 | −0.0101 | −0.0082 | −0.0007 | −0.0056 |
| log_cross_target_var | −0.0195 | −0.0016 | −0.0074 | +0.0006 | −0.0070 |

`log_ctrl` is the only feature whose removal hurts on all four lines. Zeroing z², sign agreement, source MAD, or cross-target variance improves the learned score on average. The fit is using control expression. It is not using the reliability ingredients, and the full predictor still loses to a constant weight.

## Sensitivity: n = cell count in the fixed formula (not a primary row)

`fixed_ncells` replaces `n_src` with the mean source cell count. It was locked as a sensitivity row before the scores. It beats equal-weight on energy-weighted cosine only for K562 (+0.0179, a bit larger than `n_src`'s +0.0154) and loses on RPE1 (−0.0118), HepG2 (−0.0062), and Jurkat (−0.0148), with intervals entirely below 0. Rank and top-1 still beat equal-weight on 4/4 lines. Unweighted cosine beats on K562 and Jurkat and ties on RPE1 and HepG2. Putting cell count into the fixed formula does not change the decision, and it is worse than `n_src` on energy-weighted cosine for the three lines where `n_src` already fails.

## Data (CONFIRMED: downloaded, checksummed, knockdown check passed)

Pseudobulk / bulk only. Files live in `work/data/` (gitignored). Checksums are SHA-256 of the bytes we read.

| File | URL | Bytes | SHA-256 |
|---|---|---:|---|
| `K562_essential_raw_bulk_01.h5ad` | https://ndownloader.figshare.com/files/35773070 | 79,766,954 | `80de95e54fcbca0e0537d569b43ec92fde6bd0482801505504baebe3118dcadf` |
| `K562_gwps_raw_bulk_01.h5ad` | https://ndownloader.figshare.com/files/35774443 | 374,587,922 | `7cec96b3b76169abbf6b6ab9d10bf00d71d942d89e63292351f745e130b154db` |
| `rpe1_raw_bulk_01.h5ad` | https://ndownloader.figshare.com/files/35775581 | 95,350,546 | `603c655f1cfa41d649baf3ae63fca224cc11f297e40d4ed59d390b1e8d2e2db2` |
| `GSE264667_hepg2_raw_singlecell_01.h5ad` | https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE264667&format=file&file=GSE264667%5Fhepg2%5Fraw%5Fsinglecell%5F01%2Eh5ad | 5,614,460,941 | `e1ad7c3c5a201c861a207a858aa7e59f5e6ac1955674c415f7de0d1dadadb52e` |
| `GSE264667_jurkat_raw_singlecell_01.h5ad` | https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE264667&format=file&file=GSE264667%5Fjurkat%5Fraw%5Fsinglecell%5F01%2Eh5ad | 9,366,490,264 | `ffbe15f2c8f7ffcfd7b0ba9e6937d4ebc2d03b0179fa8234648a59bcb82c04a3` |

Figshare article: doi [10.25452/figshare.plus.20029387](https://doi.org/10.25452/figshare.plus.20029387), the raw **bulk** pseudobulks (not the 8–66 GB single-cell files). Nadig et al. 2025 is GEO **GSE264667** (raw single-cell; no bulk file is deposited). X-Atlas/Orion was not used: the public release is hundreds of single-cell parquet batches (`Xaira-Therapeutics/X-Atlas-Orion`), not a pseudobulk, and does not fit this VM's memory alongside the required screens.

K562 genome-wide and K562 essential are the **same cell line**. They are inverse-variance merged into one K562 profile before any fold, so holding out K562 hides both screens. Lines in the test: **K562, RPE1, HepG2, Jurkat**.

### What the matrices are (CONFIRMED)

`log2((CPM + 1) / (NTC CPM + 1))` of the cell-count-weighted mean profile. Replogle bulk `X` is the per-cell mean UMI profile (one row per guide); guides are pooled by cell count. Nadig `X` is integer counts; cells are streamed in chunks of 8,192 and never densified as a full matrix.

Sampling variance of the log2 ratio is the delta method on per-cell means. Replogle: Poisson, plus across-guide variance / n_guides when a target has ≥2 guides. Nadig: empirical cell-level variance (Poisson only if n < 2), capped at 25.

Known-answer, measured on the built matrices: the targeted gene's own log2 fold change is a knockdown.

| Screen | Perturbations | Genes | NTC cells | Median cells/target | Median own-gene log2FC |
|---|---:|---:|---:|---:|---:|
| K562 essential | 2,057 | 8,561 | 12,425 | 121 | **−2.38** |
| K562 gwps | 9,866 | 8,246 | 89,151 | 178 | **−2.36** |
| K562 merged | 9,870 | 8,607 | 101,576 | 194 | **−2.38** |
| RPE1 | 2,393 | 8,748 | 13,303 | 72 | **−3.11** |
| HepG2 | 2,393 | 9,623 | 4,976 | 45 | **−2.84** |
| Jurkat | 2,393 | 8,881 | 12,013 | 83 | **−2.65** |

HepG2's NTC cell count is **4,976** and Jurkat's is **12,013**, the same counts recorded for these files in `experiments/350-same-line-priors-data/RESULTS.md`. Same for total Jurkat cells (262,956) and the perturbation counts (2,393). That is a file-identity check, not a score. Median cells/target here (HepG2 45, Jurkat 83) sit slightly below 350's 48 and 85.5 because this build keeps targets with fewer than 20 cells (min 2 and 1); 350 reported the range after a floor of 20.

Aligned evaluation: 6,647 genes in the intersection, 2,394 perturbations present in at least three lines (2,391 in all four). Rank drops the panel genes that sit on that axis (2,040 when K562 is held out, 2,039 otherwise), leaving 4,607 or 4,608 genes.

## Method (locked before the LOO numbers)

Code: `gate.py` (math), `build_replogle.py`, `build_nadig.py`, `evaluate.py`. Tests: `uv run python test_gate.py` (no download).

Prediction for a held-out line is `w * μ`, where `μ` is the equal-weight mean of the other lines' log2 fold changes and `w` is per (target, gene), in `[0.15, 1]`.

**z²** = `μ² / v_μ`, with `v_μ = (mean of source sampling variances) / n_src`. This is the squared z-statistic of the mean. Cross-source disagreement is **not** inside z².

**Fixed formula** (primary row): `w = 0.15 + 0.85 · z²/(z²+8) · n_src/(n_src+4)`. The `+4` prior is on the scale of a handful of lines, so `n` is `n_src`. A sensitivity row replaces `n_src` with the mean source cell count (`fixed_ncells`); it is not a primary row.

**Learned gate:** `w = 0.15 + 0.85 · sigmoid(b + β · z(features))`, eight features, z-scored on the task's own eligible entries:

1. `log_z2` — log(z² + ε)
2. `log_n_src` — log(number of sources with that entry)
3. `log_ctrl` — mean over sources of log(control CPM + 1)
4. `log_abs_mu` — log(|μ| + ε)
5. `logit_sign_agree` — logit of the fraction of sources matching sign(μ)
6. `log_cross_target_var` — log variance across targets of μ, cis entry removed
7. `log_n_cells` — log(mean source cell count for the target)
8. `log_source_mad` — log(MAD of the source log2 fold changes + ε)

Coefficients: Nelder–Mead, maximising the mean over inner leave-one-training-line-out folds of **energy-weighted cosine** (energy = ‖truth‖²). Starts: the zero vector (exactly the equal-weight cosine, because a constant `w` does not change cosine) plus 2 draws from N(0, 0.5), seed 486. The initial simplex step is 0.5 (SciPy's default step at the origin is 0.00025 and does not move the sigmoid). Each start is capped at 500 function evaluations. The objective is evaluated in float32; reported scores are float64, and the run aborts if those two disagree by more than 0.001 on the training objective. The held-out line's response is not an input to the fit (`test_fit_does_not_read_heldout_response`).

**Folds.** Leave one line out. Evaluate targets present in the held-out line and in at least two other lines. Inner folds use the same rule inside the training lines. Inner target counts were 2,391 (K562 held out) and 2,392 (the other folds).

**Cosine** drops the queried target's own gene. **Normalised rank** (pds-like) also drops every evaluated target's gene, then for each target reports `1 − (number of other targets with strictly higher cosine) / (n − 1)`. `top1` is the fraction with nothing strictly closer. Both cosines are the downstream profile, not the cis knockdown.

**Controls.**

- `perm_refit`: one permutation of the gene axis applied to the features in every training task and the test task; the gate is refit. Weights no longer describe the gene they multiply.
- `shuffle`: at test time, eligible source profiles (and their features) are reassigned to the wrong held-out targets. The fit is the real one.

**Uncertainty.** Target bootstrap, 4,000 draws, seed 486, percentile 95% interval. Draws are paired across methods. A method "beats" another on a line when the interval for the difference in energy-weighted cosine lies entirely above 0.

**Feature contributions.** Each fitted coefficient is set to 0 in turn (no refit) and the drop in outer energy-weighted cosine is recorded. That says what the fitted linear predictor uses. It is not a refit importance.

Definitions were not retuned after the scores. The only code changes after the matrices were built were memory: the fast objective no longer copies the feature tensor, source matrices are dropped after alignment, and glibc arenas are trimmed between fits. The first run was OOM-killed (exit 137, `work/evaluate_run1_oom.log`) during the HepG2 permutation refit, after K562 and RPE1 had finished. Those two folds' coefficients match the successful rerun to the printed digits (`train_ew` 0.2608 and 0.1918, same coefficient vectors). Wall time of the successful run: K562 2332 s, RPE1 3814 s, HepG2 4760 s, Jurkat 3605 s.

## Reproduce

```bash
cd experiments/486-cloud-reliability-gate
uv run python test_gate.py
uv run python build_replogle.py    # needs work/data/*.h5ad
uv run python build_nadig.py
uv run python evaluate.py          # writes work/results/
uv run python energy_quartiles.py  # post-hoc table from per_target.npz
uv run python format_summary.py
```

`evaluate.py` writes `work/results/summary.json`, `scores_by_line.csv`, `gains.csv`, `coefficients.csv`, `ablation.csv`, `fold_{line}.json`, and `per_target.npz`. Response matrices (`work/responses/*.npz`) and the h5ad downloads are gitignored. `per_target.npz` is tracked so the quartile table can be rebuilt without a four-hour refit.

## Coordinator's disposition (2026-09-27 06:50Z; brought to master from `cloud/486-reliability-gate`)
Accepted. A learned per-(target, gene) reliability gate fit on other lines does not transfer (0/4 held-out lines on
energy-weighted cosine), which agrees with 473's NO-GO for the fixed multiplier on the current (sharpened) bulk. Closed; not
re-proposed without a new mechanism. The fixed formula's rank/unweighted gains are noted, but sharpening at q 1.75 already
reweights the same genes (473).

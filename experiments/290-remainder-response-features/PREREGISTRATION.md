# 290 — PREREGISTRATION

**Frozen 2026-09-06, before any feature–remainder number exists.** Hashed into
`work/PREREGISTRATION.sha256`; every script aborts if the digest of this file does not match.

**Question.** 284 showed that after 21 build-time annotation features have predicted a gene's
cross-line conservation, the gene-level REMAINDER reproduces across two halves of the 28 sources
that share no source and no cell line (ρ +0.3277 vs a ceiling of +0.3816). What carries it?
This experiment asks the responses themselves: features computed on ONE half's sources, predicting
the OTHER half's remainder, so nothing on the target side is seen by the features.

**Hypothesis (R284 §6's described tail).** A response is conserved when it engages gene-specific
downstream programmes, and state-dependent when it engages the recipient's proliferation / stress
axis.

**Nothing here is buildable.** Analysis only: no build, no `vcc`, no `cell-eval2`, no download,
nothing under `submissions/` or `data/`, no file written outside this directory. The spent reserves
(`experiments/267-*/work/reserve/`, Broad-16 lines JVE127 / MG63 / SLR23 / UMRC3) are not read.

---

## 1. Frozen inputs (path + sha, never a description)

| object | path |
|---|---|
| the 28 banks' cross-line cosine matrix, gene axis, source keys and lines | `experiments/279-znf-state-dependence/work/held279.npz` (`Ccross`, `keys`, `lines`, `pairs_cross`, `c_held`, `axis`, `genes`) |
| the row reader used to rebuild them | `experiments/279-znf-state-dependence/build279.py` (`residual_rows`, `load_unit`, `pair_cos`), on `experiments/276-conservation-map/lib276.py` (`read_bank`, `common_axis_paths`, `load_source_path`, `cosine_rows`) and `experiments/279-znf-state-dependence/inv279.py` (`CANDIDATES`) |
| 284's primary-partition maps and remainders | `experiments/284-conservation-family-split/work/rem284.npz` (`C1, C2, r1, r2, m1, m2, ok, n1, n2, gene, fold, comp`) |
| 284's per-gene table | `experiments/284-conservation-family-split/work/pred284.csv` |
| 284's components and folds | `experiments/284-conservation-family-split/work/splits284.npz` (`comp`), `lib284.component_folds(comp, 5)` |
| the 21 annotation features and `r_bar` | `experiments/281-what-predicts-conservation/work/feat281.csv` (`set == "held"`, 17,106 genes) |
| hallmark gene sets | `data/external/broad16line/metadata/hallmark_symbols.gmt` |

All loading / metric / model / bootstrap functions are **copied verbatim** from `lib284.py` into
`lib290.py` (never imported across experiments: `rem284.py` checks its own preregistration digest).
Seed **20260906** throughout.

## 2. The partition (284's primary, declared before any 290 number)

Half 1 = **DLD-1, HCT116, HEK293T, HT-29, K562, MCF7, RPE1, THP-1, iPSC-Tian19, melanoma**
(10 lines, 14 sources, 81 cross-line pairs).
Half 2 = **A549, BxPC-3, Calu-3, H1, HAP1, HEK293, HepG2, Jurkat, TeloHAEC, cardiomyocyte,
iPSC-neuron** (11 lines, 14 sources, 87 cross-line pairs). The halves share no source and no line.

**Replication partitions, declared now, before any 290 number is read:** 284's **partition #2** and
**partition #3** from the same seeded draw (`rem284.py:123-131`), i.e. the next two entries of the
same list — no other partition may be substituted after the fact.

## 3. Known answers (run first; abort the experiment if any fails)

| row | required |
|---|---|
| **KA-1** 279's `c_held` over all 364 cross-line pairs, recomputed from the banks with the copied `residual_rows` / `load_unit` / `pair_cos` | max abs difference **0.0e+00** on all 17,106 genes vs `held279.npz['c_held']`; median +0.018605 |
| **KA-2** 284's primary-partition maps and remainders, rebuilt from `Ccross` + `feat281` + the component folds | `C1, C2, r1, r2, ok` reproduce `rem284.npz` at max abs difference **0.0e+00**; ceiling Spearman **+0.3816** and remainder Spearman **+0.3277** to 4 d.p., n = 2,848 |
| **KA-3** replication partitions #2 and #3 rebuilt the same way | ceiling / remainder reproduce `rem284.json` to 4 d.p. (#2 +0.4622 / +0.3389; #3 +0.5518 / +0.4043) |
| **KA-4** shuffled-target control (§7, control C1) | the observed out-of-fold ρ lies outside the shuffled null's spread |

A KA row that does not reproduce aborts the run (AGENTS.md:1375-1382).

## 4. The unit

**Genes with ≥ 1 cross-line pair in BOTH halves of the partition** — `ok` in `rem284.npz`,
**2,848 genes** for the primary partition. Every number is over this unit unless labelled otherwise.

## 5. Response-derived features — formulas, frozen

Let `S(h)` be the sources of half `h`. For source `s`, `blk_s` are the residual rows of the held
genes on the common 2,739-gene axis (`residual_rows`: each bank's rows centred on that bank's mean
over the held rows it carries); `nrm_s(g) = ||blk_s(g)||`; `u_s(g) = blk_s(g)/nrm_s(g)` is the
unit row that enters 279's cosine; `mu_s` is the mean raw row the centring subtracts.

Every per-source value is defined only where source `s` measures gene `g` with a non-zero row.
**Aggregation is the unweighted mean over the sources of the half that measure `g`**; a gene with no
such source gets NaN and the model's standardiser imputes the training median (`lib284.Standardiser`).

| column (suffix `_h1` / `_h2`) | formula |
|---|---|
| **(a)** `logmag` | `log( mean_s [ nrm_s(g) / median_g' nrm_s(g') ] )` — the scale-normalised residual row norm of 279's `mag`, averaged over the half's sources, logged |
| **(b)** `breadth` | `mean_s log1p( #{ j : |blk_s(g)[j]| > τ_s } )`, with **τ_s = the 99th percentile of \|blk_s\| over every (row, axis-gene) entry of that bank's held rows**. The compacted banks carry no control-variance matrix (only `b_ctrl`, a control mean), so the bank's own residual distribution is the noise scale; the rule fixes the expected count at 1 % of the axis (27.4 genes) per bank and makes the feature a *relative* DE breadth, comparable across banks with different effect units |
| **(c1)** `gen_cos` | `mean_s [ u_s(g) · mu_s/||mu_s|| ]` — signed alignment with the bank's generic response direction. (The bank's mean *residual* row is identically zero by construction — DEFINITION-276 §2 centres each bank — so the generic direction is the mean RAW row, the same object before centring.) |
| **(c2)** `gen_abscos` | `mean_s \| u_s(g) · mu_s/||mu_s|| \|` |
| **(c3)** `gen_pc5` | `mean_s Σ_{k≤5} ( u_s(g) · v_k^s )²` — fraction of the unit row's energy in the bank's top-5 right singular vectors, fitted on that bank's centred residual rows. Defined only for banks with **≥ 200 held rows** (else NaN for that source) |
| **(c4)** `gen_pc20` | the same with `k ≤ 20` |
| **(d)** `pe_<SET>`, 9 columns | `mean_s [ u_s(g) · w_H ]`, `w_H` = the unit-normalised indicator of hallmark set `H` restricted to the 2,739-gene axis. Frozen list: `E2F_TARGETS`, `G2M_CHECKPOINT`, `MYC_TARGETS_V1`, `P53_PATHWAY`, `UNFOLDED_PROTEIN_RESPONSE`, `APOPTOSIS`, `MTORC1_SIGNALING`, `OXIDATIVE_PHOSPHORYLATION`, `INTERFERON_GAMMA_RESPONSE` |
| **(d-summary)** `pe_prolif_abs` | `mean_s \| u_s(g) · w_P \|`, `w_P` = the unit-normalised sum of the `E2F_TARGETS`, `G2M_CHECKPOINT`, `MYC_TARGETS_V1` indicators on the axis |
| **(f)** `ontarget` | `mean_s u_s(g)[g]` where `g` is on the 2,739-gene axis, else NaN |
| **(f-ind)** `on_axis` | 1 if `g` is on the 2,739-gene axis, else 0 |

**(e) within-half consistency is NOT a feature.** The mean pairwise cosine of the gene's rows across
the half's sources IS that half's own conservation `C⁽¹⁾`; it is reported as the **strongest
baseline** ("conservation measured elsewhere predicts conservation here"), together with `r⁽¹⁾`
itself, which is 284's +0.3277 and the empirical ceiling for any half-1-derived predictor.

**Frozen model feature set (17 columns):** `logmag`, `breadth`, `gen_cos`, `gen_abscos`, `gen_pc5`,
`gen_pc20`, the 9 `pe_*`, `pe_prolif_abs`, `ontarget`, `on_axis` — i.e. (a), (b), (c), (d), (f).
Coverage (`n_src`, the number of the half's sources measuring the gene) is **not** in it; it is
reported as a separate baseline and as a labelled post-hoc addition.

## 6. Targets, models, splits

* **Primary target** `y = remainder_half2` (`r2` = `C⁽²⁾` minus the 21-feature ridge's out-of-fold
  prediction of `C⁽²⁾`), predicted from half-1 features. **Secondary target** `y = C⁽²⁾` raw.
  **Symmetric direction:** half-2 features → `remainder_half1` and `C⁽¹⁾`.
* **Models.** Ridge (**primary**; alpha by inner 5-fold CV over `lib284.ALPHAS`, 15 values) and a
  shallow `HistGradientBoostingRegressor` (**secondary**; `lib284.GBM_GRID`, 4 settings, inner
  5-fold CV). `lib284.fit_ridge` / `fit_gbm` verbatim.
* **Splits.** Out-of-fold prediction under 284's **grouped 5-fold over connected components of the
  HGNC-family graph** (`lib284.component_folds(comp, 5)`), so a test gene shares no family with its
  training rows.
* **Intervals.** Component block bootstrap (2,000 draws, block = 284's `component`) **and** 281's
  first-family block bootstrap (`lib284.block281`), both reported. Paired deltas from the same
  resample (`lib284.block_boot`, `paired_boot_rho`).

## 7. Baselines and controls

**Baselines** (same target, same unit, same folds):
B1 the annotation model — 284's S1 ridge on the 21 `lib284.FEATS`, refit to predict this target
out-of-fold; B2 the annotation model's prediction of `C⁽²⁾` as a single covariate; B3 reliability
`r_bar` alone; B4 magnitude alone (`logmag_h1`); B5 **`C⁽¹⁾` alone** — the within-half conservation,
the reference for "conservation measured elsewhere predicts conservation here"; B6 `r⁽¹⁾` alone —
284's own statistic, the empirical ceiling for a half-1-derived predictor; B7 the mean predictor;
B8 coverage alone (`log n_src_h1`).

**Controls:**
C1 **genes shuffled** — the training target permuted, 200 independent shuffles, the 281 KA-4′ form,
ridge only (the GBM is too slow for 200 refits; stated, not hidden).
C2 **same-half features** — half-2 features predicting `remainder_half2`: the in-sample leakage
reading, labelled as such, to show the size of the leak a careless design would have.
C3 **random programmes** — the 9 hallmark indicators replaced by random gene sets of matched size
drawn from the 2,739-gene axis, 20 draws; the rest of the feature block unchanged.

## 8. Frozen predictions and pass rules

Every rule is judged on the **primary target** (`remainder_half2` from half-1 features), the
**ridge**, out-of-fold, on the 2,848-gene unit, with the **component** block bootstrap.

* **P1 — response-derived features predict the held half's remainder.**
  PASS iff out-of-fold Spearman **≥ +0.10** AND the component-bootstrap 95 % interval **excludes 0**.
* **P2 — they beat the annotation model on the SAME target.**
  PASS iff the paired Δρ (response ridge − annotation ridge B1), 95 % interval from the same
  resample, is **entirely above 0**.
* **P3 — proliferation / stress engagement → LESS conserved.**
  Statistic: Spearman(`pe_prolif_abs_h1`, `remainder_half2`). PASS iff **< 0** and the interval
  excludes 0. Secondary readings reported beside it: a stress-axis summary
  (`|u·w_S|`, `w_S` = P53 + UPR + APOPTOSIS indicators, computed post-hoc from the frozen columns)
  and the 9 individual `pe_*` correlations.
* **P4 — genericness → LESS conserved.**
  Statistic: Spearman(`gen_pc5_h1`, `remainder_half2`). PASS iff **< 0** and the interval excludes 0.
  Secondary: `gen_pc20`, `gen_cos`, `gen_abscos`.

**Replication.** P1–P4 are re-read on partitions #2 and #3 and on the symmetric direction; those
are reported, not part of the frozen verdict.

**The experiment's one line reports P1–P4 with their measured numbers, whatever they are.** A
well-supported negative is a valid result; nothing here is re-tuned after a rule fails.

## 9. Post-hoc, labelled

Permutation importance on held genes; the 20 best- and worst-predicted genes; whether 284's two
tails (NMD / EJC machinery above prediction, proteostasis / condensin below) differ in (d) — which
programmes each tail's responses engage. Anything written after the preregistered numbers are read
is labelled POST-HOC and cannot change a verdict.

## 10. Multiplicity, stated in advance

Verdict cells: **4** (P1–P4). Everything else — 2 targets × 2 directions × 2 models × 2 bootstraps,
8 baselines, 3 controls, 3 partitions, per-feature correlations for 17 columns, the attribution and
the tails — is secondary, comparative or post-hoc and is labelled as such in `RESULTS.md`.

## 11. Limits accepted in advance

1. Both halves come from the same 28 banks on disk and the same 2,739-gene axis; the line-group
   partition removes shared sources and shared lines, not shared processing pipelines.
2. A bank's PCs and its generic direction are fitted on rows that include the gene's own row; with
   thousands of rows the self-contribution is negligible but it is not zero, and it is a
   *within-half* quantity either way, so it cannot see the target.
3. The component bootstrap is degenerate where the giant family component dominates the unit
   (284 §7.1); the first-family bootstrap is reported beside it.
4. The remainder is defined by 284's annotation model. "The features do not carry it" is a statement
   about those 21 features.

# 279 — PREREGISTRATION: are zinc-finger perturbations less conserved across cell lines?

Written 2026-09-06, experiment 279, **after** the metadata-only inventory (`inv279.py` →
`work/inventory279.json`, `work/inventory279_sources.csv`: `perts.txt` / `genes.txt` / the HGNC
gene-group table only — no `delta.npy` was opened, no response value of any kind was read) and
**before** any conservation number on a held gene was computed. Nothing below may be changed after
this file's sha256 is recorded in `work/PREREGISTRATION.sha256`; anything added later is labelled
POST-HOC in `RESULTS.md`.

## 0. The claim

From R276 §5, found **post hoc** on the 300 challenge panel targets: KRAB/C2H2 zinc-finger targets
(35/300 by the symbol rule `ZNF*` / `ZKSCAN*` / `ZSCAN*`) have a median `C_raw` of **+0.0044**
against **+0.0223** for the rest (Mann–Whitney p 1.1e-05), 23.3 % of the bottom decile against 6.7 %
of the top; their measured reliability is also lower (+0.041 vs +0.079, p 6.0e-03), and after the
reliability correction the gap reads +0.130 vs +0.234 at p 0.096 (n = 28).

> **H1 (one-sided).** Perturbations of C2H2 / KRAB zinc-finger genes produce responses whose
> direction is **less conserved across cell lines** than perturbations of other genes — i.e. their
> per-gene mean cross-line residual cosine (276's frozen statistic) is **lower**.

Mechanistic reason the claim is worth a test: many KRAB-ZNFs are lineage-restricted repressors of
transposable elements whose targets differ by cell type, so their knockdown response should be
state-dependent rather than shared.

## 1. Sets

* **Hypothesis-forming set (EXCLUDED from the test):** the 300 challenge panel targets
  (`data/challenge/pert_counts.csv`) — the map that produced the observation.
* **Held test set:** every **non-panel** perturbed gene that (a) appears in `perts.txt` of at least
  two included sources **from two different cell lines**, and (b) is one of the challenge's 18,533
  gene symbols (`data/challenge/gene_names.csv`). From the inventory: **17,106 genes**.
* **Included sources (S279), frozen:** every compact-donor bank in 276's scan
  (`experiments/276-conservation-map/work/inventory276_sources.csv`) that is single-line (no
  `AGGREGATE`), human, CRISPRi or CRISPR-KO (no CRISPRa), not a replicate-duplicate of another bank
  (`dld1_low1`, `dld1_low2` are out; their mean `dld1_lowmoi` is in), and carries **≥ 5 held-set
  genes** — the same admission rule as DEFINITION-276 §1 with "panel targets" replaced by "held
  genes". That is **28 sources / 24 distinct cell lines**: the 26 of 276 minus the five that carry
  0 held genes (`kolf2`, `feng_ipsc`, `cd4_rest`, `cd4_stim8hr`, `cd4_stim48hr`) plus the seven
  276 excluded only for carrying 0–1 **panel** targets (`k562_essential`, `rpe1`, `nadig_jurkat`,
  `nadig_hepg2`, `tian2019_ipsc`, `gasperini_k562`, `wessels_thp1`).
  `nadig_jurkat` is Jurkat = challenge context A, and is kept: it carries **0/300 panel targets**,
  so no panel response can leak through it; a variant without it is a declared secondary.

## 2. Group membership rule (decided before any response was read)

Source: the **HGNC** gene-group assignment, snapshot on disk at
`data/external/xcell/hgnc_protein_coding.txt` (HGNC complete protein-coding set, columns `symbol`
and `gene_group`; genenames.org gene-group taxonomy).

* **ZNF (test group)** = HGNC `gene_group` contains **"Zinc fingers C2H2-type"** OR
  **"KRAB domain containing"**. (763 and 373 genes genome-wide; 696 in the held set.)
* **NON-ZNF (comparison group)** = held gene with **no** gene_group whose name contains "Zinc
  finger" / starts with "Zinc". (16,060 in the held set.)
* **OTHER-ZF (excluded from both)** = any other zinc-finger family (CCCH, FYVE, matrin, CCHC,
  RANBP2, DHHC, MYND, ZZ, GATA, …). (350 in the held set.) They are neither the claim's family nor
  a clean control, so they are dropped from the primary contrast and reported separately.
* Declared secondaries: **KRAB-only** (`KRAB domain containing`), and **276's symbol rule**
  (`ZNF*` / `ZKSCAN*` / `ZSCAN*`) so the held numbers are directly comparable to the in-sample ones.

## 3. The statistic — 276's, unchanged

`experiments/276-conservation-map/lib276.py` is imported and used verbatim
(`read_bank`, `common_axis_paths`, `load_source_path`, `pair_list`, `cosine_rows`); only the target
list changes from the 300 panel targets to the 17,106 held genes.

* **Common gene axis G279** = the intersection over the 28 sources of their gene symbols ∩ the
  challenge's 18,533, in the challenge's own order (DEFINITION-276 §2). Measured in the inventory:
  **2,739 genes** (276's panel axis was 2,994 over 26 sources).
* **Residual**: each source's rows for the held genes it carries, restricted to G279, minus the mean
  over exactly those rows (each source centred on its own across-**held**-gene mean). No whitening,
  no scaling.
* **`C_held(g)`** = mean over all pairs of sources from **different cell lines** that both carry g of
  `cos(R_s[g], R_s'[g])`. `n_pairs(g)` recorded. Same-line pairs are excluded from the primary
  (they measure reproducibility, not conservation) and reported as `C_within(g)`.
* **`C_rel(g)`** = mean over pairs where both sources have a measured per-gene split-half reliability
  of `cos / sqrt(r_s · r_s')`, `r` Spearman–Brown-corrected to full depth (`2r½/(1+r½)`), pairs
  dropped where either `r ≤ 0.02` (DEFINITION-276 §3). Sources with halves for held genes: DLD-1
  (its two low-MOI replicate banks, no rebuild needed) and, if the time box allows, K562-GWPS,
  Orion HCT116 and Orion HEK293T (267's batch-parity accumulators, rebuilt with
  `truth267.write_half_bank`'s arithmetic into **this** experiment's `work/`).

## 4. The frozen prediction, effect size and pass rule

* **Predicted direction (one-sided):** `median(C_held | ZNF) − median(C_held | non-ZNF) < 0`.
* **Effect sizes reported:** the difference in medians (`Δmed`), the difference in means, and the
  **AUC** = P(C_ZNF < C_nonZNF) + ½ P(equal) — predicted **> 0.5**.
* **Unit of resampling: the perturbed gene.** Interval = two-sample bootstrap over genes, resampling
  ZNF and non-ZNF independently with replacement, B = 10,000, seed 20260906, percentile [2.5, 97.5].
  A one-sided Mann–Whitney p is reported alongside but is not the pass rule.
* **Confounds named in advance:** ZNF genes are (i) lowly expressed and (ii) weak responders, and
  both depress a cosine. Matching variables, both frozen here:
  - `expr(g)` = the **median over the 1,600 CCLE lines** of that gene's expression
    (`data/external/ccle_broad_2025/ccle_cache.npz`, the object 276 used; genes absent from CCLE are
    dropped from the matched analysis only);
  - `mag(g)` = the median over the included sources carrying g of
    `‖R_s[g]‖₂ / median_{g'} ‖R_s[g']‖₂` (per-source scale-normalised response norm on G279).
* **Matching rule:** deciles of `expr` × deciles of `mag`, both computed on the held set. Each ZNF
  gene is matched 1:1, **without replacement**, to a uniformly drawn non-ZNF gene from the **same
  (expr-decile, mag-decile) cell**; a ZNF gene whose cell is exhausted is dropped and counted.
  50 independent matching draws (seed 20260906); the matched effect is the mean over draws of the
  paired median difference, and its interval is a paired bootstrap over matched pairs (B = 2,000)
  with the matching re-drawn inside each resample.
* **PASS RULE — decided now.** **CONFIRMED** iff **both**
  1. the raw `Δmed` is negative and its 95 % bootstrap interval lies entirely below 0, **and**
  2. the matched `Δmed` is negative and its 95 % interval lies entirely below 0.

  Anything else is **NOT CONFIRMED**; an interval entirely above 0 is a refutation in the opposite
  direction and will be reported as such.

## 5. Controls, all frozen

* **(a) Size-matched random-family null.** 100 random draws of `n_ZNF` genes from the held non-ZNF
  pool (seed 20260906) → the null distribution and [2.5, 97.5] band of `Δmed`. Reported alongside:
  every HGNC gene group with ≥ 30 held genes ranked by its own `Δmed`, so the ZNF family's rank
  among **real** families is visible (a family-level null, stronger than the random one).
* **(b) The other large TF families.** **Homeoboxes** (union of every HGNC group whose name contains
  "homeobox") and **Basic helix-loop-helix proteins** — the same statistic and interval. This
  answers "is it ZNF-specific or is it any TF family?".
* **(c) Single-dataset replications** (is it a dataset mix?):
  - **(c1)** Orion HCT116 × Orion HEK293T **only** — one dataset, two lines;
  - **(c2)** Replogle K562-GWPS × Replogle RPE1 **only** — one dataset, two lines;
  - **(c3)** a **balanced core**: {`k562_gwps`, `orion_hct116`, `orion_hek293t`}, restricted to held
    genes carried by all three, so every gene has exactly the same three cross-line pairs and the
    "which sources measured this gene" confound is removed by construction.
* **(d) Shuffled-null of the map itself**: the target-shuffled pairing of DEFINITION-276 §4 (source
  j's row replaced by a random other held gene of j), **50 draws** (not 200 — the held map is 57×
  larger than the panel map), pooled mean and band, to show the held map is above chance at all.

## 6. In-sample reference (labelled, not a test)

The same three numbers (`Δmed`, AUC, matched `Δmed`) recomputed on the **300 panel targets** with
S279's sources and axis, and 276's published panel numbers quoted next to them.

## 7. Descriptor, declared POST-HOC in advance (item 5)

For the held ZNF genes: is low conservation **noise** (low split-half reliability) or **reproducible
line-specific direction** (true state dependence)? Split the ZNF genes with a measured reliability
into terciles of `r`; report `C_held` and `C_rel` per tercile against the non-ZNF genes matched on
the same reliability tercile. If ZNF conservation stays low in the **top** reliability tercile and
after the reliability correction, the deficit is state dependence; if it disappears, it is
measurement.

## 8. Known-answer rows (must reproduce before any held number is read)

* **KA-1**: re-running the frozen computation on **276's** 26 sources and 300 panel targets
  reproduces `work/conservation276.csv`'s `C_raw` column to < 1e-6 max abs difference, and
  `build276.json`'s 2,994-gene axis / 315 cross-line pairs / median +0.0198672 / share ≤ 0 = 16.0 %.
* **KA-2**: 276's post-hoc panel ZNF observation reproduces with the symbol rule: 35 targets, median
  +0.0044 vs +0.0223, Mann–Whitney p ≈ 1.1e-05.
* **KA-3**: the DLD-1 replicate-pair reliability median on the panel = **+0.086321** (rel276.json).
* **KA-4** (only if the halves are rebuilt): 271's per-target split-half reliability medians
  +0.04706 / +0.02756 / +0.02825 for K562 / HCT116 / HEK293T on the panel targets.

Any KA that does not reproduce aborts the experiment (AGENTS.md §1375-1382).

## 9. Scope

Measurement only. No build, no `vcc`, no `cell-eval2`, no model, no submission, nothing written
outside `experiments/279-znf-state-dependence/`, nothing under `data/` or `submissions/` touched.
`uv run python`, peak RSS under 8 GiB, seed 20260906 everywhere.

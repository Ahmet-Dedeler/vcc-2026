# 281 — What predicts a gene's cross-line conservation? Build-time gene features against 279's 17,106-gene map

**CONSERVATION IS PREDICTABLE FROM GENE FEATURES: YES by the frozen rule, PARTLY in substance** —
held-family Spearman **+0.320 [+0.256, +0.379]**, R² **+0.111 [+0.074, +0.146]** (block bootstrap over
the held families; ridge, 21 build-time features, 3,422 held genes in 481 held HGNC families); over
**reliability alone** by **−0.031 [−0.071, +0.012]** on its own and by **+0.056 [+0.040, +0.070]**
when added to it — and that margin falls to **+0.012 [−0.003, +0.026]** on the measurement-matched
balanced core, where a single measured reliability covariate reads ρ +0.520 against the features'
+0.382.

**The result that outlives the verdict is at the family level, not the gene level.** Out-of-fold
predicted family medians reproduce 279's whole family axis at Spearman **+0.811 [+0.641, +0.916]**
across 51 HGNC families (balanced core), and 279's zinc-finger deficit is **84 % predicted** from
build-time biology alone (observed Δmedian −0.01785, predicted −0.01506). The axis 279 named — core
machinery conserved, cell-identity families state-dependent — **is essentially the essentiality ×
expression axis**. Per gene, the same model explains 11 % of the variance (23 % on the core) and is
beaten on rank by one binary feature ("is this gene in Replogle's essential library", ρ +0.367).

**Nothing here is buildable and none is proposed.** Measurement only: no build, no `vcc`, no
`cell-eval2`, no model for the competition, nothing under `submissions/` or `data/`, no production
file touched. `uv run python` throughout. Peak RSS **1.22 GiB** (`feat281`, `ru_maxrss`
1,313,587,200 B), 0.37 GiB (`fit281`), 0.38 GiB (`core281`) — the ~6 GiB cap was never approached.

**Date 2026-09-06, 07:58–08:40 UTC (10:58–11:40 Ankara).** Deterministic, seed **20260906**.
`PREREGISTRATION.md` sha256 **`ca6b9314310ce944b47d8d7c015bae9a232c7f031f40b5bdab49e290938e967d`**,
frozen **2026-09-06 08:15:18 UTC** (`work/PREREGISTRATION.sha256`) — after the feature table existed
and its coverage was read, before a single feature–target association was computed. Everything
labelled POST-HOC was written after the preregistered numbers were read and cannot change the verdict.

One file was written outside this directory and removed: a `nohup` line for `rel281.py` ran from the
repo root and left a 176-byte `work/rel281.log` there holding only the interpreter's "can't open
file" error. It was deleted (`rm work/rel281.log`); the pre-existing root `work/` directory and its
other files were untouched. `find` over the repo for files newer than `PREREGISTRATION.md` now
returns, outside this directory, only other agents' work (`docs/`, `experiments/282-*`).

---

## 0. Known answers

| row | required | measured | verdict |
|---|---|---|---|
| **KA-1** 279's held map read off its deliverable | n 17,106; median `C_held` +0.0186; share ≤ 0 0.306; ZNF Δmedian −0.01047; n_ZNF 696 | **17,106**, **+0.018606**, **0.306**, **−0.01047**, **696** | PASS |
| **KA-2** my CCLE reader vs 276's own `ccle_breadth`/`ccle_median_tpm` on the 300 panel targets | max abs diff < 1e-6 | **0.00e+00** / **4.99e-07** | PASS |
| **KA-2b** my Replogle-essential flag vs 276's `essential` column | 300/300 | **300/300** | PASS |
| **KA-3** the DepMap artifact 269 downloaded | sha256 `d15514918130…`, 419,330,454 B, 1,150 models, 300/300 panel targets covered | **exact digest match**, **419,330,454 B**, **1,150**, **300/300** | PASS |
| **KA-4** shuffled-target control **as preregistered** | held Spearman's 95 % interval contains 0 | ridge **−0.052 [−0.084, −0.018]** | **FAIL — the row was mis-specified, see below** |
| **KA-4′** shuffled-target control, corrected (200 independent shuffles) | the observed ρ must sit outside the null's spread | null **−0.005 ± 0.098** (ridge, 200 draws), **−0.002 ± 0.058** (gbm, 40 draws); **0/200** and **0/40** draws reach the observed; z **+3.3** / **+5.5** | PASS |
| **KA-5** 279 §4(c3)'s balanced core, read off the column | ZNF median +0.01181, non-ZNF +0.02966, Δ −0.01785 on n 8,541 | **+0.011813 / +0.029663 / −0.017850**, n **8,541** | PASS |
| **regression** boot281 re-derives fit281's four held ρ | 0.0e+00 | **0.0e+00** on all four | PASS |
| **regression** fam281 re-derives core281's OOF ρ | +0.3398 | **+0.3398** | PASS |

**KA-4 as written was my error, recorded rather than papered over.** A gene bootstrap conditional on
*one* shuffle draw is not the null distribution: with the training labels shuffled the ridge
coefficients are a random projection of 21 correlated features, and the family-held test set arrives
in whole blocks, so a single draw lands at |ρ| ≈ 0.05 with a *narrow* gene-bootstrap interval around
it. The correct control resamples the shuffle: over 200 draws the null is centred at −0.005 with
**SD 0.098**, and 0/200 reach +0.320. Keeping α at the value chosen on the real data makes this null
*more* dispersed than a fully-refitted one would be (a re-chosen α would shrink toward the mean
predictor), so the control is conservative. **Its second consequence matters for every number below:
the gene bootstrap understates uncertainty under a family-held split, so every interval reported as
the result is a BLOCK (family) bootstrap** (`boot281.py`), which reads ±0.031 SD on the primary cell
against the gene bootstrap's ±0.016.

---

## 1. The feature table (item 2 of the mandate)

`feat281.py` → `work/feat281.csv` (17,406 rows = 17,106 held genes + the 300 panel targets, 43
columns) and `work/feat281.json` (provenance, coverage, KA rows). **21 frozen features, none derived
from a measured perturbation response.** Coverage is on the 17,106 held genes.

*(One correction against the preregistration: its §3 heading says "frozen, **22** columns". The
feature **set** is exactly the one its table specifies — 7 essentiality + 3 expression + 5 regulatory
+ 4 complex + 2 evolution — but that is **21**, not 22; the count in the frozen text was mine and was
wrong. Nothing downstream uses it; the fitted models use the 21 columns the table names.)*

| block | column | source (path / URL, terms) | coverage |
|---|---|---|---:|
| essentiality | `dm_mean_effect`, `dm_sd_effect`, `dm_frac_dependent` (< −0.5), `dm_frac_strong` (< −1.0), `dm_min_effect`, `dm_common_essential` (dependent in ≥ 90 % of models) | **DepMap 24Q2 Public** Chronos `CRISPRGeneEffect.csv`, 1,150 models × 18,443 genes, already on disk at `experiments/269-utility-router-feasibility/work/depmap/` (figshare 25880521; DepMap Public, CC BY 4.0) | 0.980 |
| essentiality | `replogle_essential` | Replogle 2022 essential library, `data/donors/k562_essential/perts.txt` (2,021 genes) — 276's convention | 1.000 |
| expression | `ccle_breadth` (share of lines with log2 TPM > 1), `ccle_median_tpm`, `ccle_sd` | CCLE / DepMap 2025 mRNA TPM, 1,600 lines × 19,209 genes, `data/external/ccle_broad_2025/ccle_cache.npz` — 276's convention (KA-2) | 0.987 |
| regulatory role | `log_collectri_regulon_size`, `collectri_indegree`, `is_tf_collectri` | CollecTRI via OmniPath, `data/external/regulon/collectri_omnipath.tsv` (64,516 edges, 1,201 TFs) | 1.000 |
| regulatory role | `log_dorothea_regulon_size`, `is_tf_dorothea` | DoRothEA A–C, `data/external/omnipath-dorothea/dorothea-abc.tsv` (32,629 edges, 429 TFs) | 1.000 |
| complexes | `humap_n_complexes`, `humap_max_complex_size` | **hu.MAP 2.0** (Drew et al. 2021), `humap2_complexes_20200809.txt`, 6,965 complexes — **fetched** to `work/`, 512,145 B, sha256 `a84c112cee83…`, `https://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt` (open download, no registration) | 1.000 |
| complexes | `reactome_n_complexes`, `reactome_max_complex_size` | **Reactome** `ComplexParticipantsPubMedIdentifiers_human.txt`, 16,169 human complexes, matched by the gene's first UniProt accession — **fetched** to `work/`, 3,690,987 B, sha256 `ad536e76c397…` (Reactome, CC0) | 1.000 |
| evolution | `log_n_paralogues` | **Ensembl BioMart** (useast mirror, current release, homologs page), 2,740,633 human paralogue pairs — **fetched** to `work/bm_paralog.tsv`, 70,310,896 B, sha256 `5e67b77cf031…` (Ensembl, no restriction on use) | 0.982 |
| evolution | `ortholog_depth` (0–5: mouse 19,267 / zebrafish 14,359 / fly 8,999 / worm 8,063 / yeast 3,806 human genes with a homolog) | Ensembl BioMart, same source, five one-column queries — `work/bm_orth_*.tsv`, 871,979 B total | 0.982 |

Total fetched: **~75 MB** into `work/` (85 MB on disk incl. the two test files), against the
mandate's 500 MB cap. `work/bm_fetch.sh` is the fetch script with the exact queries.

**Declared additions, never in the primary model.** `log_n_pairs` / `n_lines` (the nuisance
covariate: how many source pairs measured this gene) — excluded in advance because it is a property
of the map's measurement design, not of the gene; reported in §5. `fam::<family>` one-hot of the 64
HGNC families with ≥ 30 held genes — random split only. `r_bar` (the gene's measured split-half
reliability, 279 §6) — **the confounder**, reported as its own baseline.

**Two features the mandate named that do not exist here.** (i) 276's `shared_programme_share` is
computed from the response rows through the production whitener, i.e. from measured perturbation
responses — excluded by the mandate's own rule. (ii) A gene-age/phylostratum table: none was on disk
and none small and citable was reachable, so evolutionary depth is represented by the ortholog-depth
proxy above. DepMap's own portal is behind a Cloudflare challenge and CORUM's site is a SPA with no
scriptable download; both were replaced (269's on-disk DepMap artifact, hu.MAP + Reactome for
complexes).

---

## 2. Splits (item 1)

| split | construction | held |
|---|---|---:|
| **PRIMARY, family-held** | block = first HGNC `gene_group`, else `NOFAM:<cytoband>`; 2,425 blocks permuted with seed 20260906 and taken until ≥ 20 % of genes | **481 blocks, 3,422 genes (20.00 %)** |
| comparison, random-gene | seed 20260907 | 3,421 genes (20.00 %) |
| **second held set** | the 300 panel targets of 276 (disjoint from the held genes: overlap **0**) | 300 |
| descriptive | grouped 5-fold CV over the same blocks → out-of-fold predictions for all 17,106 | 17,106 |

Every preprocessing step (median imputation, standardisation, α and GBM hyperparameters by 5-fold
inner CV) is fitted inside the training rows only. Chosen on the primary cell: ridge α = 1,000; GBM
`learning_rate` 0.03, `max_iter` 200, `max_depth` 3.

---

## 3. The preregistered result (items 3–4)

**Family-held split, y = `C_held`, block (family) bootstrap, 2,000 draws.**

| arm | Spearman ρ | R² | AUC bottom decile |
|---|---|---|---|
| **ridge, 21 features** | **+0.3204 [+0.2558, +0.3791]** | **+0.1113 [+0.0744, +0.1457]** | **0.7030 [0.6645, 0.7366]** |
| **gbm (depth 3), 21 features** | +0.3160 [+0.2443, +0.3769] | +0.1151 [+0.0727, +0.1524] | 0.6755 [0.6351, 0.7101] |
| baseline: best single feature (`replogle_essential`, calibrated) | **+0.3669 [+0.3073, +0.4133]** | +0.0994 [+0.0621, +0.1288] | 0.5753 [0.5582, 0.5917] |
| baseline: **reliability alone** (`r_bar`) | **+0.3517 [+0.2899, +0.4092]** | **+0.1747 [+0.1099, +0.2295]** | **0.5065 [0.4685, 0.5447]** |
| features + reliability, ridge | +0.4075 [+0.3438, +0.4683] | +0.1955 [+0.1369, +0.2461] | 0.5710 [0.5341, 0.6088] |
| features + reliability, gbm | +0.4148 [+0.3481, +0.4771] | +0.2211 [+0.1443, +0.2871] | 0.5761 [0.5373, 0.6176] |
| baseline: mean predictor | n/a (constant) | −0.0007 [−0.0091, −0.0000] | 0.5000 |

**Paired margins** (same block resample):

| Δ Spearman | value | P(Δ ≤ 0) |
|---|---|---:|
| features − reliability alone | **−0.0312 [−0.0714, +0.0115]** | 0.929 |
| (features + reliability) − reliability alone, ridge | **+0.0556 [+0.0400, +0.0698]** | 0.000 |
| (features + reliability) − reliability alone, gbm | **+0.0618 [+0.0394, +0.0836]** | 0.000 |
| features − the single best feature | −0.0453 [−0.0852, −0.0064] | 0.990 |

**The frozen rule (PREREGISTRATION §6) fires YES**: the better of ridge/gbm has ρ ≥ 0.20 with its
interval entirely above 0, and R² > 0 with its interval entirely above 0.

Four things about that verdict, all measured:

1. **A single measured reliability covariate is at least as good as 21 biology features** (+0.352 vs
   +0.320; the paired Δ is −0.031 with P(Δ ≤ 0) 0.93). Features nevertheless carry information
   reliability does not: added to it they buy **+0.056 [+0.040, +0.070]** of ρ, and the partial
   association after removing `r_bar` on the rank scale is **+0.239** (raw +0.321, ρ(pred, `r_bar`)
   +0.306; `extra281.json`).
2. **The two arms find different things.** Reliability is *at chance* for the bottom decile
   (AUC 0.507 [0.469, 0.545]) while the features reach **0.703 [0.665, 0.737]**. Reliability predicts
   the *level* of conservation; the features find the **state-dependent tail**. This is the one place
   where the biology features clearly beat the confounder.
3. **One binary feature beats the whole model on rank** (`replogle_essential`, ρ +0.367, Δ +0.045
   over ridge) — but it cannot rank inside its two classes, so its bottom-decile AUC is only 0.575.
4. **The random-gene split reads *lower* than the family-held one** (ridge +0.2522 [+0.2181, +0.2833]
   gene bootstrap, R² +0.0763; gbm +0.2569): holding out whole families concentrates between-family
   variance in the test set, which is easier to rank. Adding the 64-family one-hot to the random
   split changes nothing (ridge +0.2577, gbm +0.2630) — **family identity carries no information the
   features do not already have.**

### 3.1 The second held set: the 300 panel targets

Trained on the family-held training genes, scored on 276's panel map (a different axis, 2,994 genes,
and 26 rather than 28 sources — a *transfer* held set, not an identically-distributed one). Block
bootstrap over its 184 family blocks.

| arm | Spearman ρ | R² | AUC bottom decile |
|---|---|---|---|
| ridge | +0.2517 [+0.0756, +0.3652] | +0.0567 [−0.0236, +0.1095] | 0.6041 [0.4934, 0.6912] |
| gbm | **+0.2879 [+0.1571, +0.3735]** | +0.0562 [−0.0633, +0.1418] | 0.6227 [0.5280, 0.7290] |
| the panel's own reliability (271 `ref_rel`) | +0.2983 [+0.1827, +0.4161] | — | — |
| Δρ features − reliability | −0.0597 [−0.2757, +0.1230] | | P(Δ ≤ 0) 0.690 |

The model **transfers**, at about 80–90 % of its held-family Spearman, with intervals roughly 3×
wider at n = 300. It does not beat the panel's own reliability there either.

*(The panel's reliability is the mean of 271's `ref_rel` on `k562_gwps` and `orion_hct116` only.
Programme rule 2 reserves the `orion_hek293t` fold's split-half B truths on the 300 panel targets;
that fold was not read. No entry in `docs/CONFIRMATION-LEDGER.md` is required.)*

### 3.2 The reliability-corrected target

`y = C_rel` (279 §6ii), family-held, block bootstrap:

| arm | ρ | R² | AUC10 | n |
|---|---|---|---|---:|
| ridge | +0.2117 [+0.1643, +0.2511] | +0.0227 [+0.0121, +0.0317] | 0.6724 [0.6388, 0.7087] | 2,098 |
| gbm | +0.1895 [+0.1308, +0.2488] | +0.0142 [−0.0028, +0.0301] | 0.6367 | 2,098 |
| reliability alone | +0.1930 [+0.1423, +0.2405] | +0.0181 | 0.5509 | 2,098 |
| Δρ features − reliability | +0.0188 [−0.0343, +0.0754] | | P(Δ ≤ 0) 0.248 | |
| **the same model on the panel's `C_rel`** | **+0.0033 [−0.1326, +0.1379]** | −0.0143 | 0.5053 | 227 |

Dividing out reliability costs a third of the rank signal (ρ +0.320 → +0.212) and four fifths of the
variance explained (R² +0.111 → +0.023) on the held map, and **erases it on the panel**. The panel
cell is 227 targets with a wide interval, but it is a clean negative at its own scope.

---

## 4. Attribution — which features carry it (item 4)

**Per-feature Spearman on the held-family test genes** (gene bootstrap, 1,000 draws), and
**permutation importance** of the GBM on the same genes (10 repeats, held split only):

| feature | held ρ | perm. ΔSpearman | perm. ΔR² |
|---|---:|---:|---:|
| `replogle_essential` | **+0.3669** [+0.3418, +0.3917] | +0.0501 | +0.0371 |
| *(`r_bar`, the confounder, for scale)* | *+0.3517* | — | — |
| `dm_frac_strong` | +0.3217 | +0.0000 | +0.0000 |
| `dm_frac_dependent` | +0.2849 | +0.0246 | +0.0274 |
| *(`log_n_pairs`, the nuisance covariate)* | *+0.2760* | — | — |
| `dm_common_essential` | +0.2741 | +0.0000 | +0.0000 |
| `dm_mean_effect` | −0.2705 | +0.0046 | +0.0054 |
| `dm_min_effect` | −0.2700 | +0.0043 | +0.0017 |
| `dm_sd_effect` | +0.2587 | −0.0034 | +0.0009 |
| `ccle_breadth` | +0.2487 | +0.0015 | +0.0011 |
| `ccle_median_tpm` | +0.2479 | **+0.0560** | +0.0134 |
| `humap_n_complexes` | +0.1740 | −0.0000 | +0.0008 |
| `humap_max_complex_size` | +0.1538 | +0.0012 | −0.0002 |
| `ortholog_depth` | +0.1279 | +0.0004 | +0.0004 |
| `reactome_max_complex_size` | +0.1189 | +0.0012 | +0.0008 |
| `log_n_paralogues` | **−0.1185** | +0.0036 | +0.0002 |
| `reactome_n_complexes` | +0.1182 | +0.0007 | +0.0018 |
| `ccle_sd` | −0.0652 | −0.0015 | −0.0011 |
| `log_collectri_regulon_size` | **−0.0168** [−0.0474, +0.0109] | −0.0000 | −0.0000 |
| `is_tf_collectri` | −0.0168 [−0.0430, +0.0128] | +0.0000 | +0.0000 |
| `log_dorothea_regulon_size` | −0.0104 [−0.0357, +0.0158] | +0.0000 | −0.0000 |
| `is_tf_dorothea` | −0.0104 [−0.0341, +0.0161] | +0.0000 | +0.0000 |
| `collectri_indegree` | +0.0015 [−0.0325, +0.0312] | +0.0010 | +0.0003 |

**Feature-block ablation**, family-held, ridge / gbm Spearman (`extra281.json`, and on the balanced
core in `rel281.json`):

| block | alone (raw) | alone (core) | dropped (raw) | dropped (core) |
|---|---|---|---|---|
| essentiality (7) | **+0.292 / +0.286** | **+0.332 / +0.345** | +0.252 / +0.268 | +0.317 / +0.328 |
| expression (3) | +0.246 / +0.255 | +0.295 / +0.292 | +0.296 / +0.289 | +0.346 / +0.341 |
| complexes (4) | +0.167 / +0.174 | +0.188 / +0.208 | +0.323 / +0.315 | +0.388 / +0.383 |
| evolution (2) | +0.149 / +0.166 | +0.178 / +0.199 | +0.321 / +0.314 | +0.386 / +0.393 |
| **regulon (5)** | **+0.000 / +0.015** | **+0.030 / +0.030** | +0.320 / +0.312 | +0.377 / +0.392 |
| all 21 | +0.320 / +0.316 | +0.382 / +0.390 | — | — |

**Three clean statements.** (i) *Essentiality and expression carry everything*: either block alone
reaches 77–90 % of the full model, and dropping both together is what costs. (ii) **Regulon size and
TF status carry nothing** — ρ +0.0002 as a five-feature block on the raw target, +0.030 on the core,
zero permutation importance; the mandate's regulon hypothesis reads at chance. (iii) Complex
membership and evolutionary depth are individually informative (ρ +0.15…+0.21) but **add nothing**
once essentiality and expression are in (dropping them moves ρ by ≤ 0.007).

---

## 5. The measurement-design control (POST-HOC, and the reason the verdict is "PARTLY")

`C_held` is built from however many cross-line pairs happen to measure a gene, and that number varies
by two orders of magnitude:

| pair count | n genes | sd of `C_held` | median `C_held` |
|---|---:|---:|---:|
| 1–2 | 5,420 | **0.1260** | +0.0161 |
| 3–4 | 2,634 | 0.0443 | +0.0068 |
| 5–8 | 3,117 | 0.0415 | +0.0148 |
| 9–16 | 3,497 | 0.0398 | +0.0128 |
| ≥ 17 | 2,438 | 0.0719 | +0.0702 |

Features that say "this gene is expressed and essential" also say "many sources measured it"
(ρ(prediction, `n_pairs`) on the test genes = **+0.659**). Read inside pair-count strata, the
primary model's held Spearman is:

| stratum | n (held-family test) | sd(C) | ρ | AUC10 |
|---|---:|---:|---:|---:|
| 1–2 pairs | 1,132 | 0.121 | **−0.029** | 0.496 |
| 3–8 pairs | 1,095 | 0.047 | +0.187 | 0.578 |
| 9–16 pairs | 659 | 0.043 | **+0.393** | 0.666 |
| ≥ 17 pairs | 536 | 0.076 | +0.268 | 0.535 |

**279's `c3_balanced_core` removes the confound by construction** — the conservation of every gene on
*exactly* the same three cross-line pairs {K562-GWPS × Orion HCT116, K562-GWPS × Orion HEK293T,
Orion HCT116 × Orion HEK293T}, 8,541 held genes (KA-5). Repeating the whole preregistered analysis
there (`core281.py`; 6,941 train / 1,600 test in 375 held blocks; block bootstrap):

| arm | ρ | R² | AUC bottom decile |
|---|---|---|---|
| features, ridge | **+0.3816 [+0.3032, +0.4514]** | **+0.2406 [+0.1765, +0.2876]** | 0.6107 [0.5554, 0.6702] |
| features, gbm | +0.3902 [+0.3108, +0.4596] | +0.2659 [+0.1906, +0.3225] | 0.6071 |
| **reliability alone** | **+0.5200 [+0.4469, +0.5855]** | **+0.5287 [+0.4558, +0.5803]** | 0.5479 [0.4915, 0.6053] |
| features + reliability, ridge | +0.5319 [+0.4608, +0.5954] | +0.5376 | 0.5699 |
| best single (`replogle_essential`) | +0.3345 | +0.1531 | 0.5660 |
| Δρ features − reliability | **−0.1404 [−0.1926, −0.0941]** | | P(Δ ≤ 0) 1.000 |
| Δρ (features + reliability) − reliability, ridge | **+0.0119 [−0.0029, +0.0257]** | | P(Δ ≤ 0) 0.046 |
| Δρ same, gbm | +0.0086 [−0.0108, +0.0270] | | P(Δ ≤ 0) 0.189 |

**Both directions of the control fire.** Predictability is *not* an artefact of the pair count — on
measurement-matched data it is **higher** (ρ +0.382, R² +0.241 vs +0.320 / +0.111). But once the
measurement is matched, **the gene features add essentially nothing to what the gene's measured
reliability already says**: +0.012 of ρ, interval touching zero. On the raw target the +0.056 margin
is largely the features supplying pair-count information that `r_bar` alone does not carry.

Out-of-fold over every gene (grouped 5-fold on the same blocks): raw target ρ **+0.2641 [+0.2491,
+0.2787]**, R² +0.0892, AUC10 0.6885 on 17,106 genes; balanced core ρ **+0.3398 [+0.3200, +0.3595]**,
R² +0.2272, AUC10 0.6019 on 8,541.

---

## 6. The family axis — the result that matters (POST-HOC)

Per gene the model is weak. Per **family** it reproduces 279 almost completely. Using the grouped
out-of-fold predictions on the balanced core, over the **51 HGNC families with ≥ 20 genes**:

> **Spearman(observed family median `C3`, out-of-fold predicted family median) = +0.8115
> [+0.6410, +0.9155]** (family bootstrap, 2,000 draws).

| family (sorted by observed median) | n | observed median | predicted median |
|---|---:|---:|---:|
| Mediator complex | 20 | +0.2965 | +0.1327 |
| DEAD-box helicases | 30 | +0.1109 | +0.1268 |
| RNA polymerase subunits | 25 | +0.0888 | +0.1445 |
| Ribosome biogenesis factors | 24 | +0.0872 | +0.1360 |
| RNA binding motif containing | 102 | +0.0789 | +0.0481 |
| Proteasome | 21 | +0.0723 | +0.1434 |
| AAA ATPases | 37 | +0.0621 | +0.0531 |
| … | | | |
| Forkhead boxes | 42 | +0.0174 | +0.0182 |
| PRD class homeoboxes | 44 | +0.0135 | +0.0155 |
| HOXL subclass homeoboxes | 50 | +0.0121 | +0.0195 |
| **Zinc fingers C2H2-type** | **654** | **+0.0117** | **+0.0202** |
| NKL subclass homeoboxes | 39 | +0.0080 | +0.0189 |
| PDZ domain containing | 30 | +0.0030 | +0.0322 |

**279's preregistered zinc-finger deficit, through this model**: observed Δmedian on the core
**−0.01785**, **predicted from build-time gene features −0.01506 — 84 % of it**. What 279 measured as
a family-specific biological property is, at the median, an expression-and-essentiality property:
zinc fingers are non-essential and lowly expressed, and the model knows their conservation deficit
without ever seeing a perturbation. The unexplained remainder at the median is −0.0028, smaller than
but consistent in sign with 279's expression/magnitude-matched −0.00635 [−0.01172, −0.00079] (a
different matching, on the full held map).

This does not overturn 279's §6 conclusion that the deficit is state dependence rather than
measurement — the model's predictors *are* biology, not noise, and 279's tercile ordering stands. It
relocates it: **the axis is not "zinc fingers are special", it is "how essential and how broadly
expressed the gene is", and a family label is a reading of where that family sits on it.**

### 6.1 The residual tail (item 5), with its artefact stated first

**On the raw target the residual tail is uninformative and I say so rather than reporting it as
biology**: all 20 genes furthest *above* prediction (`SYT13` +0.451, `MAGEC3` +0.509, `MT3` +0.541,
`HTR2A` +0.430, `KRT14`, `MTNR1A`, olfactory/neuronal genes) have **`n_pairs` = 1** and near-zero
CCLE expression — a single cross-line cosine drawn from a distribution with sd 0.126. They are the
tail of the map's noise, not conserved biology (`extra281.json`).

**On the balanced core the tail is real.** A second caveat applies to the *family* version of it and
is measured: a family predicted at the top of the ranking can only have a non-positive rank residual,
and ρ(predicted family median, residual family median) = **−0.744** — so the family residual table is
mostly that ceiling/floor, and §6's observed-vs-predicted comparison is the statistic to read. At the
**gene** level the tail is informative:

**20 genes far BELOW prediction** (core-essential, well-measured, and yet *not* conserved across
lines): `SRP54`, `VPS25`, `DYNLRB1`, `BANF1`, `LRR1`, `POLE2`, `NUS1`, `STX18`, `UBA52`, `SAP30BP`,
`CCNA2`, `CHEK1`, `CDC16`, `NXF1`, `GTF2B`, `CDCA5`, `RFC5`, `NDUFS5`, `GRWD1`, `NCAPD3`. Ten of the
twenty are cell-cycle / DNA-replication / chromosome-segregation genes (`CCNA2`, `CHEK1`, `CDC16`,
`CDCA5`, `POLE2`, `RFC5`, `NCAPD3`, `LRR1`, `GRWD1`, `BANF1`) — **the machinery whose knockdown
response is most obviously conditioned on the recipient's proliferative state.**

**20 genes far ABOVE prediction** (non-essential, lowly expressed, and yet conserved): `EHMT2`,
`KIF1A`, `PRSS21`, `CKLF`, `CFD`, `GABRE`, `TOMM7`, `PTK2B`, `FGF13`, `PAGE1`, `IRF6`, `PRAME`,
`FOXL1`, `NKX6-1`, `ZNF493`, `TSPAN32`, `RHOH`, `TNNI3`, `PRDM9`, `FLNC`. Their `r_bar` is 0.02–0.24,
so they are not the noise tail; two of them are zinc fingers, consistent with §6's finding that
zinc fingers sit *above* prediction once essentiality is accounted for.

---

## 7. Declared variants and controls

| cell | result |
|---|---|
| **+ the nuisance covariate** (`log_n_pairs`, `n_lines`), family-held | ridge ρ +0.3288 [+0.2988, +0.3571] (gene boot), R² +0.1114; gbm ρ **+0.3549** [+0.3262, +0.3835], R² +0.1209 — the GBM gains +0.039 of ρ from knowing how many sources measured the gene, the ridge +0.008. Excluding it in advance was the right call and cost little. |
| random split + 64-family one-hot | ridge +0.2577, gbm +0.2630 vs +0.2522 / +0.2569 without — family identity adds nothing. |
| shuffled training target (KA-4′) | 0/200 draws reach the observed ρ; null −0.005 ± 0.098. |
| mean predictor | R² −0.0007 [−0.0091, −0.0000] on every arm; the target's own scale is not being flattered. |
| gene bootstrap vs block bootstrap | ±0.016 vs ±0.031 SD on the primary cell — every interval reported as a result is the block one. |

---

## 8. Limitations, stated plainly

1. **Held ≠ new.** The target is 279's map, built from the same 28 banks; a shared analysis artefact
   in Orion or K562-GWPS would appear in both the target and (through expression/essentiality) the
   features. No new experiment was run.
2. **The features and the confounder are entangled by biology, not by leakage.** A gene that is
   essential and broadly expressed has a large response, which is *both* measured more reliably *and*
   plausibly more conserved. This experiment can separate the two statistically (§3, §5) but cannot
   decide which is causal; the balanced-core cell, where reliability wins by ρ +0.14, is the honest
   summary.
3. **The panel is a transfer set.** Different axis (2,994 vs 2,739 genes), different source set
   (26 vs 28), n = 300 — its intervals are 3× the held map's, and its `C_rel` cell is a null at
   n = 227.
4. **No gene-age table.** Evolutionary depth is a five-species ortholog-presence proxy, not a
   phylostratum. It reads ρ +0.128 alone and adds nothing; a real gene-age table could read
   differently.
5. **`replogle_essential` is a library-design label**, not a measurement: Replogle chose those 2,021
   genes as essential/expressed hits *before* the screen. That makes it build-time observable, but it
   also makes it partly a statement about what was worth measuring.
6. **Multiplicity.** Cells that decide the verdict: 2 (family-held ridge and gbm on `C_held`).
   Reported in total: ~70 cells across the primary, the random split, the family one-hot, the
   nuisance variant, `C_rel`, the panel, 24 per-feature Spearmans, 21 permutation importances, 11
   ablation cells × 2 targets, the balanced core, four pair-count strata, and the residual tables.
   All but the two are labelled secondary, descriptive or post-hoc.
7. **The 2,739-gene common axis** costs absolute cosine (279 §8.5); every number here is on that same
   axis for every gene, so the comparisons are fair and the absolute values are not the widest-axis ones.
8. **The reserve was not touched.** Programme rule 2's `orion_hek293t` half-B truths on the 300 panel
   targets were not read; the panel's reliability uses `k562_gwps` + `orion_hct116` only.

---

## 9. What this is worth

**Scientific track — three findings.**

* **F-A. The family axis of cross-line conservation is predictable from build-time gene biology**
  (Spearman +0.811 [+0.641, +0.916] between observed and out-of-fold predicted family medians over
  51 families), and it is essentially the **essentiality × expression** axis: essentiality alone
  reaches ρ +0.332 of the model's +0.382 on measurement-matched data, expression alone +0.295, while
  **regulon size and TF status read at chance (+0.030)** and complexes and evolutionary depth add
  nothing on top. 279's "core machinery → cell-identity families" ordering is recovered without ever
  looking at a perturbation.
* **F-B. Per gene, conservation is only weakly predictable, and a single measured reliability
  covariate does the job at least as well.** Held-family ρ +0.320 [+0.256, +0.379], R² +0.111 against
  reliability's +0.352 / +0.175; on the balanced core +0.382 against reliability's +0.520, and the
  features add +0.012 [−0.003, +0.026] to it. **The map is therefore not replaceable by a prior**:
  what a model could carry from gene annotation is the family-level ordering, not the per-gene value.
  The one thing the features do that reliability cannot is find the **state-dependent tail**
  (bottom-decile AUC 0.703 vs 0.507).
* **F-C. 279's zinc-finger deficit is 84 % predicted from essentiality and expression**
  (observed −0.01785 on the core, predicted −0.01506), and in the residual zinc fingers sit *above*
  prediction. The preregistered claim of 279 stands; its scope narrows further — the operative
  variable is where the gene sits on the essentiality/expression axis, and "zinc finger" is a label
  on that position.

**Competition track: nothing.** No arm, no dose, no rendering is proposed. The one action a result
like this could key — treating identifiably low-conservation rows differently — was already closed by
277 (−0.0043 [−0.0092, +0.0006] `pds`), and the per-gene R² here (0.11, 0.24 measurement-matched) is
far too small to key anything on. This experiment changes no production file.

**If it earns a stage 2**, the measurement that would move it is the one 276 §6, 279 §9 and 278 all
converged on independently: **depth**, not more features. Every arm here is bounded by the target's
own reliability — reliability alone explains 53 % of the variance of the balanced-core conservation,
and the gene features explain 24 %. A conservation map measured at 4× the cells per condition would
raise the ceiling for both, and only then is it worth asking whether biology adds to reliability.

---

## Files

| path | what |
|---|---|
| `PREREGISTRATION.md` + `work/PREREGISTRATION.sha256` | the frozen target, splits, the feature set, models, baselines, metrics and verdict rule (sha256 `ca6b9314…`, 2026-09-06 08:15:18 UTC) |
| `feat281.py` → `work/feat281.csv`, `work/feat281.json` | the feature table (17,406 × 43) with sources, coverage and KA-1/2/2b/3 |
| `work/bm_fetch.sh`, `work/humap2_complexes.txt`, `work/reactome_complexes_human.txt`, `work/bm_paralog.tsv`, `work/bm_orth_*.tsv` | the fetched public annotation tables and the exact BioMart queries (~75 MB) |
| `fit281.py` → `work/fit281.json`, `work/fit281.log` | the preregistered fits: primary, random split, family one-hot, reliability arms, nuisance variant, `C_rel`, the panel, per-feature Spearman, permutation importance |
| `boot281.py` → `work/boot281.json`, `work/boot281.log` | block (family) bootstrap intervals and the paired margins; regression checks against `fit281` at 0.0e+00 |
| `extra281.py` → `work/extra281.json`, `work/extra281.log` | KA-4 corrected, the calibrated single-feature baseline, the reliability partial, block ablation, grouped OOF, the raw-target residual tail |
| `core281.py` → `work/core281.json`, `work/core281.log` | the balanced-core repeat (KA-5), the pair-count stratification, the core OOF and residual tail |
| `fam281.py` → `work/fam281.json`, `work/fam281.log` | the family axis: observed vs predicted family medians, and the residual-table artefact |
| `rel281.py` → `work/rel281.json`, `work/rel281.log` | `C_rel` block bootstrap and the block ablation on the balanced core |
| `finalise281.py` | adds the balanced-core columns to the deliverable |
| **`work/pred281.csv`** | the deliverable: 17,406 rows × 19 columns — `gene`, `set`, `observed` (`C_held` / panel `C_raw`), `observed_rel`, `observed_core`, `block`, `gene_group1`, `split_family`, `split_random`, the four split-specific predictions, the grouped out-of-fold predictions on the raw target and the core, `rank_resid_oof_ridge`, `r_bar`, `n_pairs` |

Reproduce: `PYTHONDONTWRITEBYTECODE=1 uv run python {feat,fit,boot,extra,core,fam,rel,finalise}281.py`
in that order, from this directory (`work/bm_fetch.sh` first if the BioMart tables are absent).

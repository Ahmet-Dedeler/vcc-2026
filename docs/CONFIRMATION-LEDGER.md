# Confirmation ledger — reserved data that agents do not keep inspecting

Rule (`docs/RESEARCH-PROGRAMME-2026-09.md` §rules 2): reserved confirmation data is scored only for a preregistered claim,
once per claim, and every access is appended here with the claim it served. Reserves are declared by the experiment
that builds them and recorded below by the coordinator.

| reserve | declared by | path / sha | accesses (date, claim, result) |
|---|---|---|---|
| `orion_hek293t` half-B truths (split-half PLOLO; keys `batchB|truth`, `randB|truth`, `guideB|truth`, 294 × 17,806) | 267 | `experiments/267-replicate-validated-oracle/work/reserve/truthB-orion_hek293t.npz` sha256 `84d3879ae7f36c023e9101906ced42c90cfbe05b2fba883db143aa00684abfda` | 2026-09-05 (267 item 4, the single scoring that defined it: cross-validated oracle −0.011 A→B / −0.032 B→A on this fold); 2026-09-06 claim C1 (below) |
| four Broad-16 lines `JVE127`, `MG63`, `SLR23`, `UMRC3` (rule: sort by name, positions 4/8/12/16, declared before any residual) | 268 | `experiments/268-broad16-context-residuals/work/reserve.json` sha256 `1ce639b509d353708feded808914c63233bcee03513ace1c3f47cb93d071b9ce`; their single-cell files not downloaded, pseudobulk columns asserted unread | 2026-09-06 claim C2 (below), the single access: files fetched, built and scored once by 283 |
| **R3** — eight scPerturb h5ad files (Adamson 2016 10X010, Datlinger 2017, Datlinger 2021, Dixit 2016, Papalexi 2021 ECCITE RNA + arrayed RNA, Shifrut 2018, Tian 2019 day-7 neuron); no cell of any of them is in any bank, map row, donor table or fold | 292 | manifest `experiments/292-reserve-declaration/work/reserve292.json` sha256 `432cbda19ae445e3c566d70335505ec375e786c51cdba4146e3f79f385bad7b7` (per-file paths + sha256 + sizes); declaration `experiments/292-reserve-declaration/RESERVE.md` sha256 `f2f5ad3de8292b98e8fa02bdad3008f35e06aff7d7bd277b171f0a1c8b34c153` | declared 2026-09-06 15:55Z from metadata only; **2026-09-06 16:14–16:15Z claim C5 (below), the single access by experiment 295 — ALL EIGHT FILES SPENT** (C5a PASS ρ +0.5016 [+0.1525, +0.5659]; C5b PASS `gen_pc20` +0.3709 [+0.0635, +0.4148]) |
| **R4** — the JoungZhang 2023 hESC TF-ORF atlas, ONE file, ORF **overexpression** (1,145,823 cells, 3,369 ORF codes); ids resolved by 301 from metadata only, `X` never opened | 301 (candidate) → **DECLARED 2026-09-07 03:33Z by the coordinator, spent by experiment 302** | `data/external/scperturb/JoungZhang2023_atlas.h5ad` sha256 `8bd00eac883b469a1e82a2beca18f647e8023904920f233dbb8568e13503fc22`, 5,804,454,475 B; id map `experiments/301-metadata-acquisitions/work/joung_orf_map.csv` sha256 `a6b60664943abc28f3a6a7e9ef257a881adb097ca0bacc363a25f6ce29323211` | declared 2026-09-07 03:33Z; **2026-09-07 03:40–03:42Z claim C6 (below), the single access by experiment 302 — the file is SPENT** (C6a/C6b/C6c all FAIL by their rules and VOID by the prereg's own KA-E clause; median mirror cosine −0.0007 [−0.0047, +0.0016] against a shuffled null of +0.0008, bounded by an ORF split-half reliability of +0.025) |

### Claim C1 (preregistered 2026-09-06T05:36:16Z, before the access): 276's leave-the-held-line-out conservation predicts per-target transfer on the RESERVED HEK293T half-B truths

- Data: `experiments/267-replicate-validated-oracle/work/pertarget267-orion_hek293t.npz` key `mean|a23|batchB|pds` (the equal-weight mean's per-target `pds` scored on the reserved half-B batch truths — the single scoring 267 item 4 defined) vs `experiments/276-conservation-map/work/conservation276.csv` column `C_lolo_orion_hek293t` (HEK293T left out of the conservation score).
- Statistic: Spearman ρ over the targets present in both, paired target bootstrap 20,000 draws, seed 20260906.
- Pass rule (frozen now): ρ > 0 with the 95 % interval above 0. Fail: interval includes 0. One access; the result is appended below whatever it is.
- Result (2026-09-06 05:40Z, the single computed access; a first script errored on target alignment before computing anything): n = 293, Spearman ρ = **+0.361**, 95 % [+0.257, +0.458] → **PASS by the rule — RELABELLED 2026-09-06 13:40Z as split-half ROBUSTNESS, not an untouched confirmation: the half-B cells are part of the full HEK293T truth that 276 evaluated `C_lolo` against, and 267 item 4 had scored half B; freezing a test before reading a derived file does not restore independence from inspected parent data.** Comparators (not reserve): half-A ρ +0.334; 276's own deep-truth reading ρ +0.365. File: `experiments/276-conservation-map/work/confirm-C1-hek293t-halfB.json`.

### Claim C2 (preregistered 2026-09-06T08:46:12Z, before the access): F5's reproducible line-specific component on the four RESERVED Broad-16 lines

Full preregistration: `experiments/283-broad16-reserve-confirmation/PREREGISTRATION.md`, sha256
`65f74d562ee6d4a55d8b93e467338b1fd8622e7723d4f37c6f21cf9ce3808fc0` (hashed 08:47:21Z, before the first byte of the
reserved lines' single-cell data was fetched). One access, this claim only.

- Data: `JVE127`, `MG63`, `SLR23`, `UMRC3` from Figshare 33273600 `single_cell_data.zip`
  (`https://ndownloader.figshare.com/files/67643568`, 3,171,208,158 B), built by 268's `build_halves2.py` (seed 268,
  `EXCLUDE_ARM=1`, authors' cell assignment); everything fitted on 268's 12 working lines only (gene filter, centring,
  basal PCA, ridge, residual programme basis); 268's frozen 82 targets; paired target bootstrap 4,000 draws, seed 283.
- **C2a** (reproducibility): per reserved line, median over targets of the split-half cosine of the line-specific
  residual (response minus the 12-working-line shared response, cell halves as in 268). **PASS iff on ≥ 3 of 4 lines
  the median > 0.028 and the bootstrap 95 % lower bound > 0.028.** Calibration disclosed before the access: the same
  rule passes on only **8 of 12 working lines** (working per-fold medians span [0.0211, 0.0652]); secondaries frozen at
  the same time: whether each reserved median lies inside that working range, the fraction of targets above the line's
  own permutation null p95 (268: 0.262 working, 0.075 OR controls), the line's own OR-control median (268: 0.0183),
  and the guide-half version (268: 0.0355 vs 0.0461 cell halves).
- **C2b** (the shared response transfers): per reserved line, median over targets of cos(12-line shared response,
  reserved line's full response). **PASS iff on ≥ 3 of 4 lines the median lies inside 268's working-line
  leave-one-line-out per-fold range [0.0400, 0.0937]** (stated from `lolo268_K2.csv` before looking; overall working
  median 0.0658).
- **C2c** (basal state does not predict the residual): 268's model (c1 per-gene ridge, c2 programme-level, K = 2, J = 5)
  fitted on the 12 working lines only, predicting the reserved lines' residuals, vs baseline (a) = shared response
  alone. **PASS — i.e. 268's negative is confirmed — iff for both c1 and c2 the paired-target bootstrap 95 % interval
  of Δcos(full response) and of Δretrieval over (a) covers 0 or lies entirely below 0. FAIL if either interval for
  either arm lies entirely above 0.**
- Known-answer rows required before any reserved number is computed: the build re-run on working line `T3M5`
  reproduces 268's `halves/T3M5.npz` at 0.0e+00, and the 12-fold working-line analysis reproduces 268's
  `lolo268_K2.csv` at 0.0e+00. Abort on either.
- **Result (2026-09-06 08:48–08:58Z, the single access; four `*_GE_filtered_feature_bc_matrix` members fetched 08:48Z,
  built and scored once). Known-answer rows first, all before any reserved number: the build re-run on `T3M5` reproduced 268's
  `halves/T3M5.npz` at 0.0e+00; the 12-fold working-line run reproduced `lolo268_K2.csv` at 0.0e+00 over 45 columns ×
  984 rows (SHUFFLED arms included); each reserved line's pseudobulk reproduced the on-disk column 100/100. n = 82
  targets and G = 12,199 genes on all four lines, no target dropped.**
  - **C2a FAIL** (by the letter): 1 of 4 lines clears median > 0.028 with the bootstrap lower bound > 0.028 — `MG63`
    0.0502 [+0.0293, +0.0579] passes; `SLR23` 0.0563 [+0.0190, +0.0794] and `UMRC3` 0.0362 [+0.0273, +0.0495] fail on the
    interval, `JVE127` 0.0173 [+0.0131, +0.0265] on the level. Under the disclosed working-line rate 8/12, P(≤1 of 4) =
    0.111. **The secondaries all confirm the component**: paired against each line's own permutation null the residual
    reads **+0.0347 [+0.0229, +0.0468], P(≤0) = 0.0000**; the within-line OR-control null is the same on the reserve as
    on the working lines (0.0148 vs 0.0147) and all 4 lines exceed their own (1.8–3.9×); pooled median 0.0345
    [0.0265, 0.0451] vs working 0.0455 [0.0414, 0.0497] = 76 %; 3 of 4 medians inside the working range [0.0211, 0.0652];
    guide halves 0.45× the cell-half value vs 0.75× working. Post-hoc (labelled, not preregistered): the failures track
    **sequencing depth** — Spearman(UMI/cell, residual cosine) +0.650 (p = 0.0064) over 16 lines, the five shallowest
    lines all fail, and the reserve drew a depth-poor quarter (median 17,428 vs 30,944 UMI/cell).
  - **C2b PASS**: 3 of 4 medians inside the frozen [0.0400, 0.0937] — `MG63` 0.0595, `SLR23` 0.0543, `UMRC3` 0.0493;
    `JVE127` 0.0313 below it. Retrieval 0.684–0.768 vs the working lines' 0.700–0.819.
  - **C2c PASS** (268's negative confirmed on held data): c1 Δcos(full) **−0.0050 [−0.0077, −0.0022]**, Δretrieval
    −0.0033 [−0.0168, +0.0103]; c2 Δcos **+0.0000 [−0.0026, +0.0029]**, Δretrieval −0.0047 [−0.0178, +0.0086]. No
    interval above 0. Secondary (detectable, not usable): residual cosine c1 +0.0084, c2 +0.0101 vs +0.0003 for the
    wrong recipient; oracle line offset Δcos +0.1935 but Δretrieval −0.0345.
  - Files: `experiments/283-broad16-reserve-confirmation/RESULTS.md`, `work/c2_results.json`
    (sha256 `03de57ffb8511c7c82518a5edc056bbcf94f20f0d6d4c76af44d95cae1d4d0a7`), `work/lolo283_RESERVE.csv`
    (`4a368bbe2f5babd6219c62b7aef37e8514ddf0446573b188600fd5ef34104118`).
  - **Headline: C2a FAIL, C2b PASS, C2c PASS. The reserve is spent; these four lines are not to be used again.**

### Claims C3 and C4 (preregistered 2026-09-06T15:13:32Z, before any feature–outcome number on either surface): F13's response-geometry rule on two surfaces OUTSIDE the 28-bank conservation map

Full preregistration: `experiments/291-f13-independent-surfaces/PREREGISTRATION.md`, sha256
`3a4f5f8534469fba203268ab1a0d96afb2d9ce2e310ba898156ccb9edffa78da` (hashed 15:13:32Z; only counts and
column availability had been computed — no correlation, no fit, no outcome distribution).

**Neither surface is a declared reserve and neither is untouched data** (268 / 283 / 289 analysed Broad-16;
287 read the same GWCD4i supplementary table against ZNF labels, donor reproducibility and MASH DE counts).
What is frozen here is that neither has ever been compared with 290's response-geometry features or with any
model fitted on the 28-bank map. This is an out-of-map generalisation test, not an untouched confirmation.
The four reserved Broad-16 lines (`JVE127`, `MG63`, `SLR23`, `UMRC3`) are SPENT and were not read.

- Model (frozen, both claims): 290's 18 response-geometry columns aggregated over a stated bank subset, ridge
  fitted to predict 279's `C_held` on all 17,106 held genes (α by family-component grouped CV), then APPLIED to
  the surface's genes. No surface outcome enters any fit. Baseline `M_ann` = the same with 284's 21 annotation
  features.
- **C3 (Broad-16, 12 working lines, Cas9 knockouts, lines and dataset disjoint from all 28 banks).** Unit: the
  82 non-OR Broad-16 targets available on all 12 working lines and in the held map. Outcome: DEFINITION-276 §3
  on 268's full responses — mean over the 66 line pairs of the cosine of the line-centred residual rows
  (12,199-gene filter, `min_cpm` 5.0); second outcome 268's per-target `shared_full`. Statistic: Spearman,
  paired target bootstrap 2,000 draws, seed 20260906. Baselines: `M_ann`, `C_held` itself, magnitude.
  Null: 200 target shuffles. Depth (rule 9) reported, never in the model.
  **PASS iff** the `M_resp` interval is entirely above 0 AND ρ(`gen_pc20`, outcome) > 0 with its interval
  entirely above 0, on the mean-pairwise-cosine outcome.
- **C4 (GWCD4i, the authors' own published CD4-vs-K562 conservation).** Unit: the 1,773 genes of
  `K562_comparison.suppl_table.csv` with `logfc_pearson_r` that are in the held map, per condition; the frozen
  primary condition is **Stim48hr**. Features from the **23 K562-free banks** (excluding `gasperini_k562`,
  `jiang_K562`, `k562_essential`, `k562_gwps`, `perturbmatch`), so the K562 side of the outcome shares no data
  with the features; a labelled sensitivity arm includes them. Statistic: Spearman, component-block bootstrap
  over 284's components, 2,000 draws. KA before anything: the three `random_r*` columns read at chance
  (|median| < 0.01). Baselines: `M_ann`, `C_held`, `n_degs_MASH_K562`, the authors' `donor_correlation_mean`.
  **PASS iff all four hold in Stim48hr**: ρ ≥ +0.15 with the interval above 0; paired Δ over `M_ann` above 0;
  partial Spearman given `n_degs_MASH_K562` + `n_degs_MASH_Stim48hr` above 0; and `gen_pc20` AND `breadth`
  each positively associated with the outcome, intervals above 0.
- Known-answer rows required before any claim number: 290's primary cell rebuilt from the per-source store
  (+0.317385 at 0.0e+00); 279's `c_held` read from `held279.npz` (n 17,106, median +0.018605); 287's GWCD4i
  random-control row; 268's per-target decomposition recomputed from `work/halves/` reproducing
  `per_target_decomp_cpm5.csv` at 0.0e+00. Abort on any failure.
- Result: appended below whatever it is.

- **Result (2026-09-06 15:15–15:20Z, the single computation). Known-answer rows first, all before any claim number:** 290's primary
  cell rebuilt from `persrc290.npz` reads ρ **+0.31738512** (six decimals +0.317385, diff 1.2e-07) on 2,848
  genes with folds [2114, 201, 182, 162, 189]; 279's `c_held` read from `held279.npz` n **17,106**, median
  **+0.0186052**, max abs diff vs `znf279.C_held` 5.0e-07; the GWCD4i random controls at chance (medians
  **+0.00396 / +0.00365 / +0.00382**, real **+0.3030 / +0.3147 / +0.3117**); 268's per-target decomposition
  recomputed from `work/halves/` reproduces `per_target_decomp_cpm5.csv` at **5.5e-08 / 2.5e-08 / 8.2e-08**
  (float32 CSV) and `decomp268_cpm5.json`'s `cell_halves.pooled_full_depth` **0.090681 / 0.115379 / 0.793940**
  at 3.5e-08. (Note for later readers: `per_target_decomp.csv` and `decomp268.json` WITHOUT the `_cpm5`
  suffix are an eight-line partial run — the 12-line objects are the `_cpm5` ones.) Model fitted on the map:
  M_resp α 31.6, grouped-CV ρ **+0.3833** against `C_held`; M_ann α 3,162, **+0.2263**.
  - **C3 FAIL** (condition 1 only, narrowly). n = 82 non-OR targets, 66 line pairs, 12,199 genes;
    `C_broad16` median +0.02957, 0 targets ≤ 0. **M_resp ρ +0.2099 [−0.0130, +0.4263]** — interval includes 0
    → condition 1 FAILS; **`gen_pc20` +0.2952 [+0.0695, +0.5065]** → condition 2 PASSES. Baselines:
    **annotation ridge +0.2996 [+0.0826, +0.4894]** (paired M_resp − M_ann −0.0837 [−0.3570, +0.2062]);
    **`C_held` itself +0.2104 [−0.0043, +0.4144]** (paired +0.0021 [−0.1700, +0.1647]); magnitude +0.1518
    [−0.0685, +0.3639]; `breadth` +0.1857 [−0.0439, +0.3852]; out-of-fold variant +0.2152. Null: 200 target
    shuffles +0.0157 ± 0.1068, **10/200** reach the observed, z +1.8. Second outcome (268's `shared_full`,
    n 77): M_ann +0.3528 [+0.1347, +0.5474] vs M_resp +0.1118 [−0.1209, +0.3190]. Depth (rule 9, never in
    the model): ρ(conservation, median cells/KO) −0.081, but over the 66 LINE pairs ρ(pair cosine, pair
    geometric-mean UMI/KO) **+0.400** — depth acts on the line, not the target; the depth-partialled C3 is
    +0.2106. **Read as a low-power negative**: ±0.22 interval on 82 targets; the point estimate equals what
    the 28-bank map's own conservation achieves.
  - **C4 FAIL** (condition 2 only, in all three conditions). n = 1,773 genes, 408 component blocks.
    Stim48hr: **M_resp +0.2032 [+0.1296, +0.2517]** (cond 1 PASS, ≥ +0.15); **paired M_resp − M_ann −0.0641
    [−0.1236, +0.0112]**, P(Δ ≤ 0) 0.964 (cond 2 **FAIL** — annotation +0.2702 [+0.1922, +0.3111] is not
    beaten); **partial given `n_degs_MASH_K562` + `n_degs_MASH_Stim48hr` +0.0610 [+0.0296, +0.1563]**
    (cond 3 PASS); **`gen_pc20` +0.1858 [+0.1189, +0.2324]** and **`breadth` +0.1275 [+0.0160, +0.1555]**
    (cond 4 PASS). Null: 0/200 gene shuffles reach it, z +8.4. Other baselines: `C_held` +0.2739
    [+0.2104, +0.3319]; `n_degs_MASH_K562` +0.2959; the condition's own DE count **+0.5376**; the authors'
    `donor_correlation_mean` +0.1009 [−0.0615, +0.1819] (n 590). Labelled sensitivities: 28-bank features
    +0.2279; K562-free training target +0.1846; out-of-fold +0.1920. Rest and Stim8hr give the identical
    verdict pattern (1 PASS / 2 FAIL / 3 PASS / 4 PASS).
  - **What survives.** F13's inverted SIGN transfers to both surfaces as a single feature (`gen_pc20`
    +0.2952 and +0.1858, intervals above 0; surviving the authors' own magnitude at +0.0720). F13's HEADLINE
    — response geometry beats annotation by +0.3948 — does NOT: the ordering reverses on both surfaces and on
    268's shared-energy fraction. That +0.3948 was a property of 290's target (the remainder AFTER the
    annotation model), not of conservation. The 28-bank map itself predicts both surfaces as well as any model
    fitted on it (+0.2104 vs +0.2099; +0.2739 vs +0.2032). 290 §8.4's DNA-binder exception does not reproduce
    on the CD4 surface and inverts for annotated TFs (mean rank residual −0.0550).
  - Files: `experiments/291-f13-independent-surfaces/RESULTS.md`, `work/surfA291.json`, `work/surfB291.json`,
    `work/ka291.json`, `work/desc291.json`, `work/pred291_broad16.csv`, `work/pred291_gwcd4i.csv`.
  - **Headline: C3 FAIL (low-power, condition 1), C4 FAIL (condition 2). F13's sign is confirmed off-map on
    two surfaces; F13's "beats annotation" is withdrawn outside the map.** Neither surface was a reserve;
    both remain re-usable, and both have now been looked at with these features once.

### Reserve R3 (declared 2026-09-06T15:55Z by experiment 292, from METADATA ONLY — no claim yet)

Full declaration: `experiments/292-reserve-declaration/RESERVE.md`, sha256
`f2f5ad3de8292b98e8fa02bdad3008f35e06aff7d7bd277b171f0a1c8b34c153`. Inventory over every dataset on disk:
`work/inventory292.csv` (`2bb6ec7dc015e9bafb0b38e395317a6f533c1e48505a638836d8ec07e1a74e8a`, 104 data rows,
USED 62 / UNUSED 42). Digest list: `work/sha256.txt` (`315bca013d9ebb72ccb88921f3ad5cfd4c737481ee6cc8bc25a827e1fed2ccef`, 55 files).

**RESERVE R3 DECLARED: 8 datasets, 8 human genetic-perturbation sets, 174 perturbations in the held map
(pooled; 152 in the 7-file core), 1 new line; NOTHING COMPUTED BEYOND METADATA.**

**The files** (all under `data/external/scperturb/`; 1,575,574,586 B, 427,073 cells, 85,859 control cells):

| file | sha256 | bytes | line | modality | genes ≥ 50 cells (held map / panel) | cells/gene med, q10 | ctrl cells | guides/gene |
|---|---|---:|---|---|---|---|---:|---:|
| `AdamsonWeissman2016_GSM2406681_10X010.h5ad` | `e70fcd49808cab8d724de8d5a332940911206e1c8ef44cc7b568d048ed795c85` | 471,286,951 | K562 (repeat) | CRISPRi | 82 (80 / 2) | 583, 253 | 10,007 | 1 |
| `DatlingerBock2017.h5ad` | `2ee8b95ef1b2567ce5cbf83fc3a20bf4318f32dc007e1228f064fe428125f7da` | 39,085,166 | Jurkat (repeat) | Cas9 KO (CROP-seq) | 28 (28 / 0) | 150, 74 | 1,320 | 3 |
| `DatlingerBock2021.h5ad` | `9f8c35058d71442a904391dabf6574bbb24ae31ed00cef5d4e4ed498828f3592` | 33,600,780 | Jurkat (repeat) | Cas9 KO (scifi) | 20 (20 / 0) | 1,702, 1,062 | 4,497 | 2 |
| `DixitRegev2016.h5ad` | `a448af6aa250a6cdca472d26bcba5aec857381d8aeb03024b9358aebc2d6d23e` | 309,282,071 | K562 (repeat) | CRISPRi | 10 (8 / 2) | 2,435, 1,928 | 10,247 | 2.5 |
| `PapalexiSatija2021_eccite_RNA.h5ad` | `03e9602d124c261281936717d5795445805e6e704a52b0310bcde94225aa0929` | 147,215,278 | THP-1 (repeat) | Cas9 KO (ECCITE) | 23 (22 / 1) | 662, 294 | 2,386 | 4 |
| `PapalexiSatija2021_eccite_arrayed_RNA.h5ad` | `7f8b06d92955940c453b7bec5dec879d64c1b2e2667aa5b6f7bb4bbc6665d976` | 52,339,395 | THP-1 (repeat) | Cas9 KO, arrayed | 7 (7 / 0) | 726, 288 | 2,009 | 1 |
| `ShifrutMarson2018.h5ad` | `122c31ae44177f79b2d0441d8d10bf47028af0d6e536543e8c734048f8eecf29` | 253,843,410 | **primary human CD8⁺ T cells — the one NEW line** | Cas9 RNP KO | 18 (17 / 1) | 566, 404 | 39,813 | 2 (2 donors) |
| `TianKampmann2019_day7neuron.h5ad` (RESTRICTED) | `e5883ffef02c0da955b8d6f1dd3c92ee1dc4d73f43d6a1ed9562444097fca5df` | 268,921,535 | iPSC-neuron d7 (repeat) | CRISPRi | 24 (23 / 1) | 1,314, 299 | 15,580 | 1 |

`TianKampmann2019_day7neuron` is RESTRICTED because it is the **same study** as the map bank `tian2019_ipsc`
(Tian 2019): it adds a cell state, not an independent study, and any claim using it must say so.

**Access rule (frozen).** Outcomes — responses, pseudobulks, DE, cosines, conservation values, agreement
statistics, reliabilities, any feature derived from an expression value — may be computed on these files
**exactly once, for a preregistered claim written into this ledger before the first byte of `X` is read, by an
experiment the coordinator names in that preregistration.** Metadata reads (`obs`, `var`, `uns`, digests, cell
counts, label vocabularies) are free and unlimited. **No exploratory look**: no pilot, no "quick check", no
hyper-parameter chosen on the reserve, no second pass after seeing the first. The preregistration names the
statistic, the pass rule, the baselines, the controls and the known-answer rows; the result is appended here
whatever it is, and the file is then SPENT. Each file is spent independently — spending Adamson does not spend
Shifrut — and this ledger records which. Two riders: independence is a property of the CELLS, not of the moment
the test was frozen (R1's relabelling); and depth-condition the pass rule (programme rule 9 — cells per gene
span 150 to 2,435 across the reserve, q10 74 to 1,928).

**What R3 cannot confirm.** (i) It cannot measure cross-line conservation of its own: only 3 genes (`EGR1`,
`IRF1`, `SOCS1`) appear at ≥ 50 cells in two distinct reserve studies, 0 in three — its usable outcome is
per-gene agreement with the map's response measured in OTHER lines (291's C3/C4 shape). (ii) Four of the eight
files are Cas9 knockout, not the map's dominant CRISPRi: a claim must restrict (CRISPRi 115 genes ≥ 50 cells,
110 in the held map; Cas9 KO 69 / 67) or carry modality as a declared stratum. (iii) Per-file n is 7–82 → ONE
pooled claim with the file as the bootstrap block, never eight per-file claims (291's C3 failed condition 1 at
n = 82 with a ±0.22 interval). (iv) Panel coverage is 6 targets (`STT3A`, `SYVN1`, `E2F4`, `ELK1`, `IFNGR2`,
`STAT6`; 7 with `RFK`) — nothing competition-facing may be selected or tuned on it. (v) Only Shifrut carries a
biological-donor unit, so a donor reliability law is not computable. (vi) Comparability is measured, not assumed:
each file's `var` covers 2,636–2,739 of the map's 2,739 response-axis genes.

**Known-answer rows before any reserve number** (`work/scan292_ka.json`): 292's obs scan reproduces 278's
independent obs-only scan on `n_cells` and `n_labels` for **50/50** scPerturb files; 279's held gene list and
290's feature-store gene list are element-wise identical (17,106, same order), so the mandate's overlaps (a) and
(c) are one number; the 300 panel targets read from `conservation276.csv` equal `data/challenge/pert_counts.csv`;
and the sealed HCASMC asset re-digests to `b87452c2db1f508f5013bc15dff8088828ab931b35492d77ed6aea569fb5d9b9`,
matching 126's record in `docs/DEF-RUNBOOK.md`:296 (sha256 only — not decompressed, header not read).

**Not part of R3:** the sealed IGVF HCASMC asset (`data/external/hcasmc-igvf/IGVFFI0206LUDV.tsv.gz`) stays under
its own older rule, `docs/DEF-RUNBOOK.md` §5.-3 — one shot after the D/E/F freeze, competition-facing, 27 panel
targets. R3 does not authorise touching it and 292 did not.

### Claim C5 (preregistered 2026-09-06T16:09:10Z, before the first byte of any R3 file's `X` was read): F8's map score and F13's sign predicting transfer into the eight UNTOUCHED reserve datasets

Full preregistration: `experiments/295-reserve-r3-claim/PREREGISTRATION.md`, sha256
`a8a49c637c6c24e9503513d84576030c13b911b0dbc8bb89d455751ff78b8a7d` (hashed 16:09:10Z). **Experiment 295 is the
experiment named by R3's access rule to spend the reserve** (mandate 292 → 295). One access, all eight files, one
pass each; every file listed below is SPENT when the result is appended.

- **Files (all eight of R3, digests verified against `reserve292.json` before anything else — 8/8 MATCH):**
  `AdamsonWeissman2016_GSM2406681_10X010` `e70fcd49…`, `DatlingerBock2017` `2ee8b95e…`, `DatlingerBock2021`
  `9f8c3505…`, `DixitRegev2016` `a448af6a…`, `PapalexiSatija2021_eccite_RNA` `03e9602d…`,
  `PapalexiSatija2021_eccite_arrayed_RNA` `7f8b06d9…`, `ShifrutMarson2018` `122c31ae…`,
  `TianKampmann2019_day7neuron` `e5883ffe…` (RESTRICTED: same study as the map bank `tian2019_ipsc`).
- **Unit:** 174 distinct genes (≥ 50 cells, in 279's 17,106 held genes; 205 (file, gene) rows de-duplicated by
  taking the file with the most cells for a gene, a rule fixed before the access). Modality is a declared stratum
  (CRISPRi: Adamson, Dixit, Tian; Cas9 KO: Datlinger ×2, Papalexi ×2, Shifrut); Tian is flagged and the pooled
  result is given with and without it; cells per gene and UMIs per cell are declared covariates, reported next to
  every cell, **not** conditioned on (programme rule 9).
- **Outcome `A_g`:** cosine between the reserve file's residual response of `g` (`log1p(5e4·S/ΣS)` minus the file's
  own control profile, minus the mean over the file's own scored genes; the on-target column removed) and the
  28-bank map's CONSENSUS residual response for `g` over the banks whose LINE differs from the file's (Tian also
  excludes `tian_ipsc_neuron`), on 279's frozen 2,739-gene axis restricted to the file's `var`. Split-half
  reliability (8 draws) and cells/UMIs reported alongside.
- **C5a (F8 off-map on untouched data):** ρ(`C_held`, `A_g`) pooled. **PASS iff ρ ≥ +0.15 and the file-block
  bootstrap 95 % interval (2,000 draws, the FILE as the block) excludes 0.** Baselines reported, not in the rule:
  `M_ann` (284's out-of-fold annotation ridge), `mag`, the file's own reliability, the gene's DE count, `r_bar`,
  `C_held_excl`, and the partial Spearman given reliability and magnitude. Null: 200 gene shuffles.
- **C5b (F13's sign):** ρ(`gen_pc20`, `A_g`) with `gen_pc20` aggregated over the same non-excluded banks.
  **PASS iff the file-block 95 % interval lies entirely above 0.** `breadth` reported likewise.
- **Known-answer rows before any reserve number:** KA-1 all eight digests (PASS 8/8); KA-2 279's `c_held`
  n 17,106, median +0.0186052 (PASS); KA-3 291's C4 cell rebuilt from `pred291_gwcd4i.csv`, Stim48hr ρ +0.2739
  (PASS, with `gen_pc20` +0.1858 and `breadth` +0.1275); KA-4 the map-consensus construction reproduces
  `held279.npz` `Ccross` on 6,097 (pair, gene) cosines at **0.0e+00** (PASS); KA-5 computed library size ==
  `obs['ncounts']` for every cell; KA-6 a controls-vs-controls split reads ≈ 0 against the map consensus.
  Abort on KA-5 or KA-6 failure.
- **Multiplicity:** two verdict cells (C5a, C5b). Everything else is secondary and labelled.
- **Result (2026-09-06 16:14–16:15Z, the single access; 28 s, peak RSS 2.89 GiB). Known-answer rows first, all
  before any reserve number: KA-1 8/8 digests MATCH; KA-2 `c_held` n 17,106 median +0.01860516; KA-3 291's C4
  cell +0.2739 at n 1,773 (`gen_pc20` +0.1858, `breadth` +0.1275); KA-4 the map-consensus construction reproduces
  `Ccross` on 6,097 (pair, gene) cosines at 0.0e+00; KA-5 computed library size == `obs['ncounts']` at max abs
  diff 0 over all 427,073 cells of all 8 files; KA-6 the controls-vs-controls split reads |median cos| ≤ 0.014
  against the map consensus on every file. Every file's control count and every gene's cell count reproduced
  292's frozen grammar exactly (the first run ABORTED on that assertion — Adamson's 2,613 guide-unassigned cells
  were being mislabelled by a pandas code −1; fixed to 292's own serialisation and re-run from the start, no
  statistic having been computed). n = 174 pooled genes (205 (file, gene) rows), 8 file blocks.**
  - **C5a PASS**: ρ(`C_held`, `A_g`) = **+0.5016**, file-block 95 % **[+0.1525, +0.5659]**, gene bootstrap
    [+0.3709, +0.6083] — ρ ≥ +0.15 and the interval excludes 0. Null 200 gene shuffles +0.002 ± 0.072,
    **0/200** reach it, z **+6.9**. Baselines: `M_ann` +0.4211 [+0.1275, +0.5058] (paired Δ +0.0797 [−0.0747,
    +0.1767]); magnitude +0.3578 (paired Δ **+0.1483 [+0.0382, +0.2916]**, P(Δ≤0) 0.011); the file's own
    reliability +0.4672; the gene's DE count +0.3934; `r_bar` +0.4336; `C_held_excl` (every bank of the file's
    line dropped from both sides of each map pair) +0.4593; cells per gene −0.0760. **Partial given reliability
    and magnitude +0.3086 [+0.1389, +0.3964]**; adding log cells and log UMI/cell +0.3189 [+0.0638, +0.3989].
  - **C5b PASS**: ρ(`gen_pc20`, `A_g`) = **+0.3709**, file-block **[+0.0635, +0.4148]** (above 0);
    `breadth` +0.4365 [+0.1225, +0.5207]. `gen_pc20` does not beat annotation (paired Δ −0.0407 [−0.2244,
    +0.0942]) — as in C3/C4, only F13's SIGN travels off-map.
  - **Strata (declared):** pooled without Tian +0.5396 [+0.1628, +0.5749] / +0.3804 [+0.0017, +0.4235];
    CRISPRi (n 109) +0.5169 [−0.0435, +0.7381] / +0.3415 [−0.2619, +0.3567]; Cas9 KO (n 65) +0.2137 [−0.0766,
    +0.4260] / +0.1659 [−0.0437, +0.2971]. Same sign everywhere; neither stratum separated from 0 on 3 / 5
    blocks. Labelled secondaries: all 205 pairs +0.4928 / +0.3438; on-target column kept +0.5480 / +0.3754;
    disattenuated +0.4685 / +0.3520; **within-file ranks +0.3508 [+0.0773, +0.4852] / +0.2310 [+0.0441,
    +0.3069]** — both verdicts survive using only the ordering of genes inside one file.
  - **Depth (rule 9, reported, never conditioned):** median `A_g` +0.0381, 29.3 % ≤ 0; ρ(`A_g`, cells) −0.076,
    ρ(`A_g`, split-half r) +0.467, ρ(`A_g`, UMI/cell) +0.346, ρ(`A_g`, DE count) +0.393. Two files are shallow
    in ways the metadata could not show: `DatlingerBock2021` 460 UMI/cell (scifi), `TianKampmann2019_day7neuron`
    **15** UMI/cell (its labelled droplets reach into the background pool; `obs['ncounts']` confirms it).
  - **Descriptive (post-hoc, labelled):** best `HSPA5` +0.689, `XRN1` +0.639, `ASCC3` +0.590, `STAT1` +0.456
    (THP-1) — the proteostasis / translocon / spliceosome / proteasome end of 279's axis; worst `HYOU1` −0.251,
    `CBLB` −0.129 and `CD3D` −0.090 (the T-cell-identity genes in the new CD8⁺ context). **Papalexi pooled →
    arrayed** (7 genes): ρ(`A`) +0.643, reliability 0.262 → **0.500** and DE 1 → 1,563 in the arrayed screen with
    no gain in agreement with the map — removing guide calling raises reproducibility, not transferability.
    **Datlinger 2017 → 2021** (20 genes): ρ(`A`) +0.212, both sides near the noise floor (r 0.029 / 0.117).
    `IRF1` reads +0.049 (Dixit) / +0.177 (Papalexi pooled) / +0.167 (Papalexi arrayed).
  - Files: `experiments/295-reserve-r3-claim/RESULTS.md`, `work/reserve295.csv`, `work/stats295.json`,
    `work/access295.{csv,json,log}`, `work/map295.npz`, `work/paired295.json`.
  - **Headline: C5a PASS, C5b PASS. The 28-bank conservation map's per-gene score predicts transfer into eight
    datasets no object of the programme has ever touched (ρ +0.50; five new studies, four Cas9-knockout, one new
    cell type), and F13's inverted sign holds on a third surface — this one untouched. The map's advantage over
    the annotation model is NOT established (Δ +0.08, interval across 0), exactly as 291 found. Nothing
    competition-facing; nothing proposed for the model.**
  - **ALL EIGHT R3 FILES ARE NOW SPENT** — `AdamsonWeissman2016_GSM2406681_10X010`, `DatlingerBock2017`,
    `DatlingerBock2021`, `DixitRegev2016`, `PapalexiSatija2021_eccite_RNA`,
    `PapalexiSatija2021_eccite_arrayed_RNA`, `ShifrutMarson2018`, `TianKampmann2019_day7neuron`. They are not to
    be used for a claim again. The reserve is closed.

**Recoverable exclusion:** `JoungZhang2023_atlas.h5ad` (sha256
`8bd00eac883b469a1e82a2beca18f647e8023904920f233dbb8568e13503fc22`, 1,145,823 hESC cells, 3,369 TF ORFs) is the
largest untouched human genetic set on disk, but `obs.perturbation` is a bare numeric ORF id that no `obs`/`var`/
`uns` column resolves, and its modality is ORF overexpression. A ≤ 100 KB metadata download would resolve the
ids; that is a metadata read and would not spend a reserve. Nothing else about the file may be read until it is
declared.

### Reserve R4 — DECLARED 2026-09-07T03:33Z (candidate raised by 301 from metadata only; the JoungZhang 2023 hESC TF-ORF atlas, ids resolved)

**DECLARED, and spent by experiment 302 for claim C6 below.** 301's mandate was to remove the blocker recorded in
the "Recoverable exclusion" paragraph above, and it did; the coordinator declares the resulting candidate as
reserve **R4** and names **experiment 302** as the experiment that spends it, under the access rule below, for
the single preregistered claim C6. Independence (point (vii)) is hereby asserted by the coordinator on 292's
recorded evidence: `inventory292.csv` carries `JoungZhang2023_atlas` as `status` UNUSED / `reserve_tier`
`excluded-recoverable` ("no code reads it; only obs/var censuses"), and the file is absent from 279's
`work/inventory279_sources.csv` (the 28 map source banks). At declaration time nothing in the file had been read
beyond `obs` / `var` / `uns` — **`X` had never been opened**: no pseudobulk, no response, no DE, no reliability,
no agreement statistic.

**The file**

| file | sha256 | bytes | cells | line | modality | labels |
|---|---|---:|---:|---|---|---|
| `data/external/scperturb/JoungZhang2023_atlas.h5ad` | `8bd00eac883b469a1e82a2beca18f647e8023904920f233dbb8568e13503fc22` | 5,804,454,475 | 1,145,823 | H1/HUES hESCs (`obs.cell_line` = `hESCs`) | **ORF overexpression** (`obs.perturbation_type`), single perturbation per cell (`nperts` = 1) | 3,369 categories = 3,367 TF ORFs + 2 controls |

Cells: **1,018,010 perturbed**, 127,813 control (`GFP` 49,539 + `mCherry` 78,274). Two batches; `var` carries
`gene_symbol`, `ncells`, `ncounts`; `uns` is empty.

**The id resolution** (`experiments/301-metadata-acquisitions/RESULTS.md` §3; map
`experiments/301-metadata-acquisitions/work/joung_orf_map.csv`, meta `work/joung_map_meta.json`; sources with URL,
licence and sha256 in `work/sources301.md`). scPerturb's own curation notebook left `obs.perturbation` as the raw
category codes of the GEO file's `TF` column (`# TODO: Perturbation here is just a number. Find out what it
means!!!`). Codes were resolved by joining the atlas's SHARE-seq `cell_barcode` R-triples to GEO **GSE217460**'s
per-cell TF-barcode files `…_210322_TFmap_S01…S10.csv.gz` (public domain, NCBI GEO) and mapping the 24-nt ORF
barcode through **Joung et al. 2023 Table S1A** (`mmc2.xlsx`, CC BY 4.0, 3,550 ORFs). The sublibrary ↔ `P1` group
assignment is measured (six of ten separate at purity 1.000, none below 0.630, every runner-up ≤ 0.21);
**2,146 anchor codes are strictly increasing in TFORF number with 0 violations**, so the codes are the rank of
`TFORF####-SYMBOL` in sorted order; gaps were filled only where forced.

- **3,136 / 3,369 codes resolved** (2,146 anchored, 990 forced), 233 ambiguous (median 8 cells each, 5,917 cells).
- **99.48 % of cells covered** (1,139,906 / 1,145,823).
- **1,647 distinct TF symbols**; 1,069 with ≥ 50 cells, of which **993 are in 279's 17,106-gene held conservation
  map**. Per-ORF: 1,729 / 3,134 at ≥ 50 cells (median 65, q10 9).
- Confirmation the ordering is right, not assumed: Table S1A's last two entries are `TFORF3549` = GFP and
  `TFORF3550` = mCherry, and the map lands them on codes 3367 and 3368 — the file's two largest categories, i.e.
  exactly the two controls, **predicted by the rule rather than fitted to them**.

**Access rule (identical to R3's; FROZEN at declaration, 2026-09-07T03:33Z).** Outcomes — responses,
pseudobulks, DE, cosines, conservation values, agreement statistics, reliabilities, any feature derived from an
expression value — may be computed **exactly once, for a preregistered claim written into this ledger before the
first byte of `X` is read, by an experiment the coordinator names in that preregistration**. Metadata reads
(`obs`, `var`, `uns`, digests, cell counts, label vocabularies) are free and unlimited; 301's id resolution was
such a read and does not spend anything. No exploratory look, no pilot, no hyper-parameter chosen on it. Because
this is ONE file, it is spent as one unit — there is no per-file partial spend as R3 had.

**What it cannot confirm.**
(i) **Modality.** ORF **overexpression**, not CRISPRi knockdown. The sign of a response is not comparable to any
knockdown object in this programme — F8's map, 284's remainders, 290's geometry features and 291's surfaces are
all knockdown. A claim may use magnitude, geometry or rank structure; it may not compare signed directions to the
map without an explicit, stated re-orientation, and no such re-orientation has been validated here.
(ii) **Context.** One cell line (hESCs), one state, two batches, no biological donor unit — no line, state or
donor reliability law is computable, and it repeats an embryonic-stem context that 273 already flagged as
unrepresentative of the challenge's somatic references.
(iii) **Panel coverage is 0.** None of the 300 2026 panel targets is among the 1,647 resolved TF symbols
(0 / 300, and 0 among the 1,069 with ≥ 50 cells). Nothing competition-facing can be selected, tuned or checked on
it, and it cannot extend 293's construct-ceiling reading.
(iv) **Isoforms, not genes.** The library is TF **splice isoforms** — up to 19 ORFs per symbol, median 1 — so a
per-gene statistic must state its isoform-pooling rule before the access, not after.
(v) **Residual ambiguity.** 233 codes (5,917 cells, 0.52 %) are unresolved and must be dropped, not guessed;
a claim's n is the resolved set. Resolution can be raised by downloading the remaining 14 `TFmap` sublibraries
(~95 MB, same public sources) — a metadata read, still free.
(vi) **Depth.** Programme rule 9 applies: cells per resolved symbol run median 110, q10 12 — an order of
magnitude below R3's range at the low end, so any pass rule must be depth-conditioned.
(vii) **Independence is not yet asserted by the coordinator.** 301 re-read, but did not re-derive, 292's evidence:
`experiments/292-reserve-declaration/work/inventory292.csv` records `JoungZhang2023_atlas` with `status` **UNUSED**,
`reserve_tier` `excluded-recoverable` and evidence "no code reads it; only obs/var censuses", and the file is absent
from 279's `work/inventory279_sources.csv` (the 28 map source banks). The coordinator owns that assertion.

**Not part of this reserve:** `data/external/scperturb/JoungZhang2023_combinatorial.h5ad`
(sha256 `ce5188bd3f924f9c79fea6d0b2fb49e48a71d0cf6a2a5010ee9a9bd6e3289a01`, 834,073,300 B, 167,947 cells,
57 combinatorial TF labels already carried as symbols by scPerturb). It was never blocked, it is a combinatorial
design, and 301 did not open it beyond `obs`. R4 does not authorise touching it and 302 did not.

### Claim C6 (preregistered 2026-09-07T03:33:19Z, before the first byte of R4's `X` was read): in one cell line, does ORF overexpression reverse CRISPRi knockdown, and does conservation or knockdown geometry predict how strongly?

Full preregistration: `experiments/302-orf-vs-crispri-symmetry/PREREGISTRATION.md`, sha256
`69c2a92415b9e61149bdbee6eb3e61c61a8b23711954d237fcf8bd19d341b9e0` (hashed 03:33:19Z; only `obs` / `var` / `uns`
metadata, 301's id map and the on-disk `gse295214_h1` CRISPRi bank had been read). **Experiment 302 is the
experiment named by R4's access rule to spend the reserve.** One access, one streaming pass over `X`; the file is
SPENT when the result is appended.

- **Data.** R4 (`JoungZhang2023_atlas.h5ad`, `8bd00eac…`, 1,145,823 hESC cells × 37,528 genes, ORF
  overexpression) × the map's `gse295214_h1` CRISPRi bank (Pan et al. 2025, H1 hESCs, 2,354 perturbations ×
  18,533 genes, `experiments/115-gse295214-h1/work/compact-donor`). ORF ids from 301's `joung_orf_map.csv`
  (`a6b60664…`, 3,136 / 3,369 codes, 99.4836 % of cells).
- **Axis `A` = 15,710 genes**: H1 bank genes ∩ the bank's row-constant `measured_mask` (16,268; `delta` is
  exactly 0 elsewhere) ∩ the atlas's 37,528 distinct `var` symbols.
- **Unit.** 175 TF symbols with ≥ 50 pooled ORF cells (isoforms pooled by symbol; 419 codes, 212,329 cells) and
  an H1 bank row; **primary cohort `T100` = 148 with ≥ 100 cells**, secondary = the other 27.
- **Objects.** ORF profile `p_s = log1p(5e4·S_s/L_s)` (`L_s` = the total counts of the same cells over all
  37,528 columns); response `r_s = p_s − p_ctrl` (GFP + mCherry pooled); **residual `a_s = r_s − mean_{t∈T} r_t`,
  in which `p_ctrl` cancels exactly**; CRISPRi residual `b_s = d_s − mean_{t∈T} d_t` (279's `load_unit` path);
  the TF's own column deleted from both (144 / 148 have one on `A`); **`S_s = cos(a_s, b_s)`**. ORF split-half
  reliability from 16 seeded cell shards × 8 balanced 8/8 partitions (seed 20260907).
- **C6a (antisymmetry exists):** median `S_s` < 0 with the paired TF-bootstrap 95 % interval (20,000 draws)
  entirely below 0, AND #{`S_s` < −0.05} > #{`S_s` > +0.05}. **PASS iff both hold on `T100`.** Null: 200
  TF-shuffled pairings (expected median 0 by construction).
- **C6b (conservation predicts mirroring):** Spearman(`C_held`, −`S_s`) on the 142 `T100` TFs with a finite
  `C_held` in `znf279.csv`, paired TF bootstrap 10,000 draws. **PASS iff ρ > 0 and the interval is above 0.**
  Baselines reported, not in the rule: knockdown magnitude, ORF magnitude, ORF reliability, `r_bar`, and the
  partial Spearman given both magnitudes and the reliability.
- **C6c (F13's sign in a new modality):** Spearman(`gen_pc20` of the H1 knockdown row, from `persrc290.npz`,
  −`S_s`) on the same 142. **PASS iff ρ > 0 and the interval is above 0.**
- **Declared secondaries (labelled, no verdict):** |`S_s`|; the 50–99-cell cohort; per-isoform `S` for the 70
  TFs with ≥ 2 isoform codes at ≥ 50 cells; the positively aligned fraction and the top 20 each way; `S_s`
  disattenuated by `sqrt(rel_s)`; the map-wide `gen_pc20` and `breadth` as a second reading.
- **Known-answer rows.** KA-A control counts 49,539 / 78,274 (**PASS**); KA-B 249 / 233 / 148 reproducing the
  coordinator's metadata check (**PASS**); KA-C the H1 bank reproduces 279's `Ccross` on all 26 cross-line pairs
  involving `gse295214_h1`, 15,047 (pair, gene) cosines, max abs diff **0.0e+00** (**PASS**); KA-D the ORF map's
  sha256 (**PASS**); and inside the single pass, before any claim number: **KA-E (ABORT)** computed library size
  == `obs['ncounts']` for every cell; **KA-F** a controls-vs-controls split (GFP vs mCherry) reads median
  |cos| < 0.05 against the knockdown rows; **KA-G (ABORT)** the positive control — overexpression raises the
  TF's own transcript (median own-gene ORF residual > 0, ≥ 80 % positive) and CRISPRi lowers it (median < 0,
  ≥ 80 % negative); **KA-H** `X` holds non-negative integers. On an KA-E / KA-G failure the three verdicts are
  recorded VOID and the numbers reported as descriptive only.
- **Multiplicity:** three verdict cells. **Depth (rule 9):** cells per TF (median 474, q10 156), ORF UMI/cell
  (median 3,502), isoform count and reliability reported next to every cell; the rules are unconditioned, as
  declared. **Panel coverage 0/300 — nothing competition-facing is selected, tuned or checked.**
- **Result (2026-09-07 03:40–03:42Z, the single access: one streaming pass over `X`, 58 row blocks, 2.674e9
  nonzeros, 77 s, peak RSS 7.04 GiB; every statistic below computed from the sums that pass wrote,
  `work/sums302.npz` sha256 `01c2a69669ff17b694b349f1b50783b0caeec686a8781fea9e9b700f76d51671`).
  Known-answer rows first, all before any claim number: KA-0 atlas digest MATCH; KA-A GFP 49,539 / mCherry
  78,274 on codes 3367 / 3368; KA-B 249 / 233 / 148 to the digit; KA-C the H1 bank reproduces 279's `Ccross`
  on 26 pairs, 15,047 (pair, gene) cosines, at 0.0e+00; KA-D the ORF map digest; KA-F the GFP-vs-mCherry
  no-signal row median |cos| 0.0097 (max 0.0461, 0 % above 0.05) — PASS; KA-G the positive control PASS, the
  ORF residual's own-gene coordinate median +0.3368 with 95.1 % positive over 144 TFs and the CRISPRi
  residual's median −0.6998 with 100 % negative (an expression-level confirmation of 301's barcode-derived
  id map, of the cell→symbol assignment and of the axis); **KA-E FAIL and KA-H FAIL.**
  - **KA-E / KA-H, diagnosed from free metadata plus the collected sums — `X` was NOT reopened:** `X` holds
    half-integers (max fractional part 0.5) and its row sums do not equal `obs['ncounts']` (max abs diff
    9.96e3). The inconsistency is the FILE's: `sum(obs['ncounts'])` 4.18355e9 vs `sum(var['ncounts'])`
    4.04904e9 (3.2 % apart, `var['ncounts']` itself fractional), while `sum(obs['ngenes'])` =
    `sum(var['ncells'])` = 2,673,739,132 = exactly `X`'s nnz, so the structure read was the file's. The
    per-group ratio computed/`obs` runs 0.882–1.182 and is batch-keyed (0.940 on the 134 all-batch-0 groups,
    1.142 / 1.149 on the two all-batch-1 controls; Pearson with the batch-1 fraction +0.910): `obs['ncounts']`
    was recorded at a different processing stage, batch by batch. Nothing in the claim's arithmetic uses it —
    `log1p(5e4·S_s/L_s)` takes both terms from the same matrix. **But the preregistration's clause fires as
    written: the three verdicts are VOID, not PASS/FAIL, and the numbers below are descriptive.** The reserve
    is spent either way and cannot be reopened to repair the row.
  - **C6a — FAIL by its rule (VOID):** median `S_s` = **−0.0007**, paired TF bootstrap 95 % **[−0.0047,
    +0.0016]** (condition 1 fails); share < −0.05 **10.8 %** (16) vs > +0.05 **4.1 %** (6) (condition 2
    holds); 52.7 % negative. Null (200 TF-shuffled pairings) **+0.0008 ± 0.0017**, z **−0.83**, 19 % of draws
    ≤ observed. **Overexpression and knockdown responses are ORTHOGONAL, not opposed, at the cohort median.**
  - **C6b — FAIL by its rule (VOID):** Spearman(`C_held`, −`S_s`) = **−0.1143 [−0.2846, +0.0673]**, partial
    given both magnitudes and the ORF reliability **−0.1057** — the point estimate has the WRONG sign.
    Baselines (none predicts the sign): kd magnitude −0.013, `mag`(290) −0.038, ORF magnitude −0.047, ORF
    reliability +0.080, `r_bar` −0.046, cells −0.000, UMI/cell −0.103, isoform count −0.000, map-wide
    `gen_pc20` −0.036, `breadth` +0.014, batch-1 fraction +0.100.
  - **C6c — FAIL by its rule (VOID):** Spearman(`gen_pc20`, −`S_s`) = **+0.0583 [−0.1106, +0.2278]** (right
    sign, interval across 0), partial +0.0653.
  - **Why, measured:** the ORF side does not resolve a direction. Split-half reliability median **+0.0246**
    (Spearman–Brown full-sample **+0.048**; 5.4 % ≤ 0), and it obeys a single Spearman–Brown law over a 40×
    cell range — +0.011 at 154 cells, +0.014 at 321, +0.032 at 633, +0.096 at 1,636, **+0.245 at 6,222** —
    with **c = 9,562 cells [6,942, 11,126]** for ρ_full = 0.5 and ≈ 3,100 cells for the H1 Flex reference's
    0.243; the atlas's median resolved symbol has 110 cells. ρ(|S|, reliability) **+0.469**. With ρ_CRISPRi
    ASSUMED at 0.243 (273; not measured on this bank) the ceiling on any cosine is **0.108**, so the observed
    median implies a TRUE mirror cosine of **−0.006 [−0.044, +0.015]**; the best-measured stratum (`r` ≥ 0.20,
    n = 21) reads −0.0180 → **−0.055 [−0.196, +0.075]**. Reliability terciles drift negative (+0.0009 →
    −0.0009 → −0.0106) but no stratum's interval excludes 0 and the trend breaks at the top (`r` ≥ 0.30 and
    `cells` ≥ 2,000 both read positive). **A strong reversal is excluded; a weak one (≈ −0.05) is not.**
  - **Secondaries (declared).** |S|: median 0.0163, and **ρ(|S|, `gen_pc20`) = +0.3823 [+0.2261, +0.5192],
    partial given reliability, both magnitudes and log cells +0.3831** — F13's geometry feature predicts the
    SIZE of cross-modality overlap, not its sign (ρ(|S|, `C_held`) −0.110). 50–99-cell cohort (n 27) median
    +0.0024 [−0.0121, +0.0113]. Per-isoform: 219 rows / 70 TFs, median +0.0029; **two ORF splice-isoform
    constructs of one TF in the same cells agree at a median cosine of +0.021 over 365 pairs** — the same
    near-zero construct term 296 found for CRISPRi libraries, with the caveat that isoforms are different
    proteins (a lower bound). Disattenuated median −0.0070; 47.3 % positively aligned.
  - **Descriptive tails (post-hoc, labelled).** Most mirrored: `ELOF1` −0.232, `DAXX` −0.174, `MYC` −0.117,
    `ZFP90` −0.117, `NANOG` −0.108, `ZNF281` −0.096, `TGIF1` −0.085, `SUPT4H1` −0.076, `KLF16` −0.074,
    `POU5F1` −0.067, `NHP2` −0.065, `YAP1` −0.064, `KLF4` −0.064 — the hESC pluripotency network and core
    transcription machinery. Most aligned: `MAX` +0.134, `SATB1` +0.097, `SUB1` +0.092, `TP53` +0.057.
  - Files: `experiments/302-orf-vs-crispri-symmetry/RESULTS.md`, `work/symmetry302.csv`
    (sha256 `9c15824aed80d984b42a548ab5f8b46c4929adec1752038c6bf1b0bf043ace9f`), `work/stats302.json`,
    `work/sums302.npz`, `work/ka302_C.json`, `work/ka_e_diagnosis.json`, `work/sblaw302.json`,
    `work/strat302.json`, `work/desc302.json`, `work/absS302.json`, `work/isoform302.csv`.
  - **Headline: C6a / C6b / C6c all FAIL by their rules and are VOID by the preregistration's own KA-E
    clause. The measured content that stands: an ORF overexpression response and a CRISPRi knockdown response
    of the same TF in the same cell line are orthogonal at the median, not opposed (a strong mirror is
    excluded, a ≈ −0.05 one is not); conservation does not predict mirroring; F13's geometry predicts the
    size of cross-modality overlap and not its sign; and ORF overexpression obeys the same Spearman–Brown
    depth law as CRISPRi with c ≈ 9,600 cells, which is why a 1.15 M-cell atlas at 110 cells per TF cannot
    answer the question. Nothing competition-facing; nothing proposed for the model.**
  - **RESERVE R4 IS NOW SPENT** — `JoungZhang2023_atlas.h5ad` is not to be used for a claim again. The
    combinatorial companion file remains untouched and undeclared.

### C6 interpretation correction — 2026-09-07 coordinator review

The original result record above is retained, but its biological orthogonality,
latent-mirror bound, ORF depth-law, isoform-noise-ceiling and matching-specific
geometry interpretations are superseded by
`experiments/304-f15-audit/COORDINATOR-REVIEW-2026-09-07.md`.
The observed cosines reproduce from saved sums; their stronger interpretations
do not follow. **C6 stays VOID, R4 stays spent.** This audit is not a fresh claim
on the reserve and does not authorize another source read.

## Reserve R5 — DECLARED 2026-09-07 04:20Z (coordinator): the 510 `Nakatake2020_hESC_TF_screen` rows of `data/external/gse295214/pseudo_bulk_perturbation_database.h5ad`
- **Recovery note, 2026-09-07 (Codex coordinator).** The predecessor's interrupted
  write left an incomplete duplicate heading before this declaration. That
  fragment has been removed; the predecessor's declaration time is retained,
  not used as the time of any C7 preregistration. No C7 entry or access receipt
  existed at takeover, and no 303 worker had launched.
- **Proposed data and metadata, pending independent preflight.** The predecessor
  recorded 510 Nakatake rows in a 20,345 × 18,108 float32 database with `control`
  and `measure_mask` layers; metadata overlaps were 100 TFs with Pan2025, 58 with
  302's cohort and 13 with the database's Joung subset. These are cohort-design
  observations, not expression results. 303 must reproduce the counts and bind
  the precise row/column sets and source hash before release.
- **Independence and interpretation are NOT yet established.** Absence of a
  literal dataset name from a text search does not prove that no generic loader
  previously read the rows. Trace 115's `extract_compact_donor.py` and the map's
  inputs. Verify the original platform, cell-line identity, induction time and
  matrix transform. "hESC" alone does not mean the same line, and bulk profiling
  does not establish an unlimited biological reliability. The original draft's
  assertions to that effect are not carried forward as facts.
- **Access remains locked.** No Nakatake `X` or `control` values may be read
  until C7, its statistic, comparisons, multiplicity rule, known answers and
  exact source/code hashes are frozen here by the coordinator. One logged
  source access only; results are recorded regardless of outcome and the rows
  are then SPENT. Metadata, measurement masks and explicitly already-used
  Pan2025 rows are allowed for preflight. This does not authorize other unused
  datasets. R1–R4 remain spent; 302's saved sufficient statistics are descriptive
  context only, never a fresh independent confirmation.
- **Assigned to 303, PREACCESS ONLY.** Governing recovered mandate:
  `experiments/303-second-overexpression-library/HANDOFF-MANDATE.md`. A draft
  preregistration or a successful preflight is NOT authorization to access.
  The coordinator must review it and create a hash-bound release first.

## Claim C7 — frozen September 7, 2026, before R5 expression access

**Status at entry: FROZEN; no expression access yet.** The coordinator has
completed the pre-access review and independently rerun the loader, analysis
and access checks. The earlier PREACCESS-only worker authorization is not the
release; a separate coordinator release follows the commit of this record.

- Protocol:
  `experiments/303-second-overexpression-library/PREREGISTRATION.md`,
  SHA256 `e92006a21674982d351f1fb17c0ef831577bc08f63da33d831247f66fde3e342`.
- Exact source/code/data mapping:
  `experiments/303-second-overexpression-library/work/FROZEN-MANIFEST.json`,
  canonical manifest digest
  `7298dfe5258fc1b5bafe5576b7317178fa0d7f9233cfd9041625aa452de4bc48`.
  The digest uses sorted-key compact JSON as implemented by the access guard;
  it is not the pretty-printed file's byte digest.
- R5 source hash remains
  `0754d6206fe045744dfd9b21eb37c8e3435152c76b158b9810a84c64a029c0ac`.
  Exact rows 184–693, 510 targets; KD-overlap cohort 100; shared measured
  axis 15,309; primary axis 15,211 after excluding all 98 cohort own columns.
- **C7a:** observed median S < 0 and negative-tail excess > 0, with the
  maximum of the two directional permutation p values entering Holm.
  **C7b/C7c/C7d:** positive conservation→−S / geometry→−S /
  geometry→|S| association beyond the matching null. Four-cell Holm .05,
  unestimable cells retained at p=1.
- Null: 2,000 uniform whole-cohort OE permutations, KD/features fixed,
  seed 20260907. Intervals: 20,000 paired target resamples, descriptive
  and conditional on the frozen centered rows. No family-held, study-level
  or causal inference claimed.
- Scope: observed cross-study/cross-line alignment; exact OE input platform
  unresolved; no replicate reliability, disattenuation or equivalence bound.
  All58 and primary45 overlaps with spent 302 are descriptive only.
- Verification: seven loader known-answer rows pass (Pan delta/control
  bit-identical, 15,047 prior cosines exact); updated preflight **16/16**
  including a **40-test** guard suite; separate manifest-bound **40/40 guard**
  and **23/23 analysis** tests, including synthetic end-to-end and CPM→VOID.
  Receipts: `work/FROZEN-GUARD-TESTS.json`,
  `work/FROZEN-ANALYSIS-TESTS.json`, both bound to the canonical manifest.
- Execution: `scripts/freeze_c7.py issue` checks committed code/protocol and
  this ledger entry. Only the coordinator then invokes `access303.py` once.
  Marker is exclusive/fsynced before expression; interruption still spends
  R5. Afterwards `analyse303.py --analyse` reads only `nak303.npz`.
  Results and any validity failure will be recorded below, not hidden or
  repaired by a second reserve read.

This is not authorization to build or submit a prediction. R1–R4 stay spent;
no other unused row subset or dataset is included.

### C7 outcome — September 7, 2026, after the single access

**R5 SPENT at 06:01:00 UTC / 09:01 Ankara.** Freeze commit `45788c6`
precedes the access. Release SHA256
`81637e9c51d5e364164cda534fd704e7209cd609481dd10d6ef045db1611c78b`;
saved product `nak303.npz`
`40065fadf89afbca554cc3f88944d9269a05be087e0f92b2e7138738b90ad599`.
The 510 rows' CPM, mask, finite-value and shared-control checks pass.

**C7a PASS:** median observed cosine −.02633077, conditional target-bootstrap
95% [−.05203106, −.00948648]; 64/100 negative; 41 below −.05 versus
13 above +.05, tail excess .28 [.14, .41]. Both directional raw permutation
p values are 1/2001; four-cell Holm **.001999**.
**C7b/C7c/C7d FAIL:** rho −.12883 / +.07606 / +.08070, all adjusted p=1.
No feature or cutoff was selected after access.

`scripts/audit_c7.py` independently reproduces the full cosine matrix
(difference 0), permutations and four verdicts from the saved product,
without importing the readout kernel or reopening the source.
See `experiments/303-second-overexpression-library/COORDINATOR-RESULT.md`
and `work/COORDINATOR-VERIFICATION.json`. Uncentered/bank-wide sensitivities
stay negative, but are descriptive. Agreement with 302 shares its KD bank
and is not independent replication. C6 stays VOID.

This is conditional observed cross-study antisymmetry, not a universal inverse,
causal modality result, reliability estimate or cleared competition candidate.
All subsequent work may use the saved sufficient product only.

## Claim C8 — 307 (no reserve: this uses data already read by 305/301)

**Frozen 2026-09-07 16:00:01Z, before any per-construct candidate expression value was computed.**
`experiments/307-offtarget-enrichment/PREREGISTRATION.md`, sha256
`81f7cda4533bc522995d2168e05542c673edacf161f5368080f0223ac801d9e7`
(`work/PREREGISTRATION.sha256`).

**No reserve is declared or consumed.** 307 reuses 305's saved integer sufficient sums
(`experiments/305-t7-expression-test/work/sums305_{train,validation,test}.npz`) and 301's saved
sequence products. No h5ad is opened, nothing is re-streamed, no download, no new seed search.

**The claim.** 305's coordinator review asked for one aggregate candidate-enrichment test with the
construct as the experimental unit, against a joint null that preserves how the candidates were
selected (seed/PAM, promoter architecture, GC, mappability) **and** their expression, with the
candidate-specific exclusion 305's null omitted. C8 is that test.

- **Unit:** the gene-targeting construct. **Primary population: all 34** targeting constructs of the
  2025 library (13 failing + 21 working; 17 paired intended targets), each in its own split. The 31
  non-targeting constructs are a separately declared secondary population, treated as
  candidate-bearing observations selected by the same rule, **not** negative controls.
- **Verified before hashing (KA-2, PASS):** 301's rule reproduces its published counts — failing
  4/13, working 4/21, NTC 4/31 at L=12 — and the 8 candidate-bearing targeting constructs named in
  305's coordinator review are confirmed, carrying exactly 305's 11 primary pairs.
- **Candidate set:** 301's own rule — the hit set at `L*(c)`, the longest seed length in
  {12, 11, 10, 9} that yields a hit inside the CRISPRi active window after the 10 kb distal filter.
  All 34 targeting constructs have a non-empty set, so **no construct is dropped for lacking a
  strong candidate**; a construct with k = 0 contributes z = 0 to observed and null alike.
- **ONE primary statistic:** `S = (1/34) Σ_c sqrt(12·k_c)·(0.5 − r_c)`, where `r_c` is the mean
  percentile of the construct's candidate genes' control-referenced log1p-CP50k effects **within
  that construct's own matched non-candidate pools**. Positive = candidates more repressed.
- **Joint null:** candidate SETS are resampled, not genes independently. Each candidate is replaced
  by a gene matched on (a) `pi_off_L` — the measured probability that a protospacer drawn from the
  209,070-spacer hCRISPRi-v2 library hits that gene off-target at seed length L under 301's exact
  rule — plus active-window anchor count and mean TSS offset; (b) promoter GC, TSS count and 10 kb
  TSS density; (c) a mappability proxy (promoter-corpus 12-mer multiplicity) and a low-complexity
  proxy; (d) basal control expression decile. Hard strata: 10 expression deciles × 5 `pi_off`
  quintiles; soft match on the remaining seven covariates, M = 200. **Eligibility excludes every
  gene any 2025 construct intends to target (300 genes), the construct's own target, and every gene
  the construct's spacers hit at any seed length** — the exclusion 305's null omitted.
- **Pass rule (ONE verdict cell, α = 0.05 one-sided, no multiplicity correction needed):** PASS iff
  `p_joint ≤ 0.05` **and** the maximum leave-one-construct-out `p_joint` over the 34 deletions is
  `≤ 0.05`. **Attainable floor `1/(B+1) = 1.0e-5` with B = 100,000 draws** — unlike 305's design
  (Holm floor 11/158 ≈ 0.0696), this one can attain significance at its stated α, so a FAIL is
  informative. Interval: construct-level bootstrap over the 34 constructs, n = 34.
- **Seven known-answer rows** (KA-1 305's published effects from its sums; KA-2 301's counts, PASS;
  KA-3 NTC sums identical across splits; KA-4 a no-signal cross-construct row; KA-5 null
  calibration; KA-6 matching balance; KA-7 293's `kd_lfc`). A failure aborts the conclusions.

**Nothing here can prove off-target binding.** Expression is observational; the constructs differ in
on-target efficacy (`kd_lfc` −3.04 … +0.03), cell number (95 … 3,415), population composition and
handling. A positive would be corroboration of a sequence-linked association only; a negative means
"no aggregate expression corroboration under this test", not "off-targets do not occur". No model,
no submission and no production change follows from either outcome.

### C8 outcome — September 7, 2026, 16:13Z

**C8: UNRESOLVABLE-AT-THIS-RESOLUTION.** The pass rule's two conditions were both met —
`p_joint` = 9.9999e-06 (the null's floor) and maximum leave-one-construct-out p = 1.99998e-05 —
but **2 of the 7 preregistered known-answer rows failed**, and §9's frozen abort clause governs:
KA-4's no-signal row measured a residual bias of **+0.0683** against a `|mean| <= 0.05` rule (its
interval condition was met), and KA-6's matching balance missed on promoter-window count,
SMD **+0.2905** against a `<= 0.25` rule (one axis of nine; `pi_off` +0.001, basal expression
-0.027). Neither threshold was loosened after the fact. KA-1 (305's 23 published rows from its
saved sums, max |diff| **0.000e+00**), KA-2 (301's 4/13, 4/21, 4/31 and the 8 candidate-bearing
targeting constructs), KA-3 (the 31 NTC constructs' sums identical across all three splits, 38,176
shared control cells), KA-5 (null calibration, KS p 0.368) and KA-7 (293's `kd_lfc`, 2.22e-16) pass.

**Bounded estimate, reported instead of a verdict.** `S` = **+0.78567** standardized units **per
targeting construct**, construct-level bootstrap 95 % **[+0.38669, +1.16604]** (n = 34 constructs,
20,000 resamples); `Delta` = **-0.079012 log1p-CP50k [-0.136647, -0.030962]** per construct, about a
7.6 % reduction relative to matched genes; 83 candidate genes over 34 constructs. Null: 100,000
selection-matched candidate-set draws, mean +0.00005, sd 0.17173, 95 % [-0.334, +0.335];
attainable floor p **1.0e-5** (contrast 305, whose Holm floor 11/158 ~ 0.0696 made rejection
impossible for any effect).

**Two further nulls agree.** The declared second null (S5: a pseudo-construct is a random pair of
stringency-0 hCRISPRi-v2 protospacers rescanned through 301's exact rule, kept only at matching
`L*` and set size; matching-free scoring) gives W_obs +0.74470 vs -0.01033 [-0.360, +0.340],
p = 0.00149, all 34 constructs resolvable. The wrong-recipient construct permutation — immune to
the matched-pool construction — gives p = 2.0e-4 (0 of 5,000 derangements reach the observed;
level +0.07168 +- 0.00222). Post-hoc, repairing the KA-6 axis changes nothing (S +0.78813,
p 1.0e-5, SMD 0.0000).

**Specificity (post-hoc extension of declared secondary S3).** The same candidate genes scored in
the cells of each construct's **partner construct — same intended target, different protospacers** —
read **+0.04106 [-0.29505, +0.36933]** versus +0.78567 in their own construct's cells; own minus
partner **+0.74461 [+0.22730, +1.24514]**, paired Wilcoxon p 0.00459, n = 34 paired constructs.
This excludes intended-target downstream biology and generic candidate-gene properties and attaches
the association to the protospacer sequence — but not to a binding event.

**Cell-half stability (post-hoc).** Scored on 305's two deterministic cell halves, which are
referenced to disjoint sets of control constructs, S = +0.75776 [+0.38759, +1.09989] and
+0.69101 [+0.28470, +1.07550]; per construct the halves agree at Pearson +0.718 (p 1.8e-6), same
sign in 27/34. This is cell sampling, not biological replication.

**Declared secondaries.** T7's original contrast is absent at the aggregate level: failing (13)
z +0.7383 vs working (21) +0.8150, difference **-0.077 [-0.869, +0.699]**, Mann-Whitney p 0.595 —
off-targeting is still neither proved as the failing-construct mechanism nor excluded. The 31
non-targeting constructs, treated as candidate-bearing observations and not negative controls, give
S +0.2857 [-0.0691, +0.6358], p 0.0537 (`PRKCQ`, `DCTN2` and a third, `CDCP1`, are the nominal
members). `GRB10_P2 -> SNX8`: e -0.492396, matched-pool percentile 0.0000; partner `GRB10_P1`
+0.045325. 33 of 83 candidates sit below 305's CP50k 0.5 floor; above-floor only, S = +0.6907.
The effect is **not** concentrated at L = 12 (L* = 12, n 8: S +0.633 p 0.0366; L* = 11, n 16:
+1.014 p 3.0e-5; L* = 10, n 7: +0.496; L* = 9, n 3: +0.652) — the most awkward fact for a simple
seed-length-graded reading. The 2025 library is not unusually promiscuous (chi-square 3.13, p 0.371
against the hCRISPRi-v2 pseudo-construct expectation).

**Nothing here establishes off-target binding.** Expression is observational, the constructs differ
in on-target efficacy, cell number, population composition and handling, there is no biological
replicate of any construct, and all three splits share one control pool. The decisive next design
is a **seed-mutant control** (~24 re-cloned constructs, one 10x lane, order $3-6 k) and, for binding
itself, a dCas9 occupancy assay ($5-10 k); with the construct as the unit, more cells do not narrow
the interval, more constructs do. No model, no submission and no production change follows.
See `experiments/307-offtarget-enrichment/RESULTS.md` and `work/enrichment307.csv`.

---

## Claim C9 — 309 (calibration repair of C8; no reserve)

**Frozen 2026-09-07 16:27:26Z, before any repaired primary value was computed.**
`experiments/309-c9-calibration-repair/PREREGISTRATION.md`, sha256
`2fdf4531f3c319a052bd73eb3c76326954e2f75362c1a0540b3974288b8febca`
(`work/PREREGISTRATION.sha256`).

**No reserve is declared or consumed.** C9 reuses 305's saved integer sufficient sums and 301's/307's
saved sequence products. No h5ad is opened, nothing is re-streamed, no download, no new seed search.
No `vcc`, no `cell-eval2`, no build, nothing under `submissions/` or `data/`, no production change.

**What C9 is.** C8 measured `S = +0.786 [+0.387, +1.166]` (n = 34 constructs), joint-null p at its
floor 1.0e-5, and then **correctly aborted** to `UNRESOLVABLE-AT-THIS-RESOLUTION` because 2 of its 7
preregistered known-answer rows failed: **KA-4** (no-signal row) measured a residual
construct-independent bias of **+0.0683** against a `|mean| <= 0.050` rule (post-hoc at 5,000
derangements **+0.0717 +- 0.0022**, i.e. **9 % of the observed +0.786**), and **KA-6** (matching
balance) missed on one axis of nine, promoter-window count SMD **+0.2905** against `<= 0.250`.
**These are defects of the null/matching machinery, not of the hypothesis.** C9 repairs exactly those
two defects and re-tests the identical hypothesis at the identical α with the identical single
verdict cell. Repairing a demonstrably miscalibrated null and re-testing is legitimate; silently
relaxing a threshold is not, and 305's coordinator decision
(`experiments/305-t7-expression-test/COORDINATOR-DECISION-2026-09-07.md`) is binding precedent
against the latter. **No threshold in C8 was loosened, then or now.**

**This is a SECOND TEST OF THE SAME HYPOTHESIS ON THE SAME DATA.** Same 34 constructs, same cells,
same integer sums, same candidate genes. It is **not** an independent replication. Across C8 and C9
there are two verdict cells on one hypothesis and one data set; family-wise control at α = 0.05 over
the two attempts requires the smaller p <= 0.025, which the attainable floor 1.0e-5 satisfies
arithmetically — **but the real multiplicity is analytic**: C9's design was chosen after seeing C8's
estimate and diagnostics. Moreover the repaired point estimate is **largely predictable in advance**
from C8's own published post-hocs (P5, subtracting the wrong-recipient level: +0.71399
[+0.31501, +1.09436]; P2, repairing the KA-6 axis: S +0.78813) — so this preregistration is **not
blind to the outcome's magnitude**, and says so. Its value is the calibration, not surprise.
**C8's bounded estimate stands and is reported alongside C9 whatever C9 says; C9 does not replace or
retract it.**

**Repair 1 — KA-4 / residual bias. Variant (a): per-construct centring, chosen and justified from the
diagnostic.** For construct `c`, `b_c` = the mean of `c`'s own candidate set scored, through `c`'s own
matched pools, in the cells of every **other targeting construct of `c`'s split** (15 / 9 / 7 wrong
recipients); the repaired score is `z_c^rep = z_c − b_c` and `S = (1/34) Σ_c z_c^rep`. Chosen over
variant (b) (scalar subtraction) because: the pedestal is not a clean constant — C8's P3 puts the
per-construct wrong-recipient spread at [−0.295, +0.369], four to five times its own mean; the
deliverable here is a construct-level interval, which (b) leaves uncentred and (a) repairs; better
covariate balance does **not** remove the bias (C8's P2: +0.081), so it needs its own repair; and (a)
is exactly C8's partner-construct test generalised from one wrong recipient to all of them. (a)
reduces to (b) if the per-construct levels are homogeneous, so it is the weaker assumption.
**Direction, stated before running: this repair works AGAINST the effect** — it subtracts ≈ +0.072
from +0.786, ≈ **9 % down**, and makes passing harder on the point estimate; its effect on the null's
spread is not signed in advance and is reported either way. **Cost, stated before running:** under
(a) the repaired KA-4 row is zero in expectation **by construction**, so its content reduces to
verifying correct implementation and Monte-Carlo bracketing of zero — an unavoidable property of any
bias-subtraction repair, which is why the structurally different checks (KA-5, the real-spacer null,
the partner test) are kept.

**Repair 2 — KA-6 / matching.** Promoter-window count (the axis that missed) is added to the hard
strata: **10 basal deciles × 5 `pi_off_L` quintiles × 3 window classes = 150 cells**; soft matching
(M = 200 nearest on the seven covariates) and the widening thresholds are unchanged. The **full
nine-axis SMD table** for the repaired pools is reported and **every axis must meet `|SMD| <= 0.250`**.
Pool sizes, widenings and any construct whose pools cannot be built at the finer stratification are
reported by name with the resulting attainable resolution — **no construct is dropped silently**. One
new global exclusion is added because the centring requires recipient-independent pools: no gene that
is the `L*` candidate of any of the 65 constructs (134 genes, ≈ 0.9 % of the eligible universe) may
be a pool member. **Its direction is NOT conservative** (it could move the estimate up by ≈ 1.9 of
200 pool members), so it is reported as a declared sensitivity with the exclusion switched off.
**No variant may be selected, swapped or blended after any repaired primary value has been seen.**

**Unchanged from C8:** the hypothesis, α = **0.05 one-sided**, **ONE** verdict cell, the statistic
family `sqrt(12 k_c)(0.5 − r_c)`, the unit (targeting construct), the population (all 34), 301's `L*`
candidate rule, 305's response formula, the null family (matched candidate-set resampling),
B = 100,000, M = 200, the leave-one-construct-out condition, and the abort clause.

**Pass rule (ONE verdict cell, α = 0.05 one-sided):** PASS iff `p_joint <= 0.05` **and** the maximum
leave-one-construct-out `p_joint` over the 34 deletions is `<= 0.05`. Attainable floor
`1/(B+1) = 1.0e-5`. **`UNRESOLVABLE-AT-THIS-RESOLUTION` iff any of the NINE known-answer rows fails
its stated rule — C8's clause verbatim, not weakened because it has already fired once.**

**Nine known-answer rows.** C8's seven re-run (KA-1 305's 23 published rows; KA-2 301's counts **plus
a clause proving the inherited 209,070-spacer pool is the one CORRECTED for the Horlbeck `="SYMBOL"`
Excel-escaping bug**; KA-3 NTC sums identical across splits; KA-4 and KA-6 on the **unchanged**
machinery, governing by **reproduction** of C8's published values, not by the thresholds C8 aborted
on; KA-5 null calibration on both nulls; KA-7 293's `kd_lfc`) **plus KA-4R and KA-6R, the two repaired
rows, which must meet C8's ORIGINAL rules** (`|mean| <= 0.050` with the interval bracketing zero;
`|SMD| <= 0.250` on every one of nine axes).

**Declared secondaries, no verdict, kept from C8:** the real-spacer resampling null; the
wrong-recipient permutation (which under repair (a) is the same object as KA-4R and is reported as
such, not as an independent check); the partner-construct specificity test (uncentred, on the
repaired pools); the failing-vs-working contrast (301's original T7 question); the 31 non-targeting
constructs run identically and labelled **candidate-bearing observations, not negative controls**;
`S` by seed level; and four sensitivities (variant (b); centring over all 33; the global exclusion
off; each repair alone).

**What a positive would and would not establish.** It would **not** establish dCas9-KRAB off-target
binding: expression here is observational, and the partner-construct test narrows the explanation to
something carried by the protospacer sequence **without separating seed-mediated binding from a
construct artefact correlated with spacer identity**. A negative means "no aggregate expression
corroboration under this test", never "off-targets do not occur".

Deliverables: `experiments/309-c9-calibration-repair/RESULTS.md`, `work/enrichment309.csv`.

### C9 outcome — September 7, 2026, 16:45Z

**C9: PASS.** `S` = **+0.71293** standardized units **per targeting construct**, construct-level
bootstrap 95 % **[+0.32419, +1.07351]** (n = 34 constructs, 20,000 resamples); joint-null
**p = 1.99998e-05** (1 of 100,000 selection-matched joint draws reached the observed value;
attainable floor 1.0e-5); maximum leave-one-construct-out p **8.99991e-05** (deleting `GRB10_P2`);
`Delta` = **-0.078205 log1p-CP50k [-0.135970, -0.030017]**. Both preregistered pass conditions met.
**Nine known-answer rows, 9/9** — so the abort clause did not fire. **C8's estimate was
S = +0.786 [+0.387, +1.166] and it stands, unretracted, beside this one.**

**The repairs did what the preregistration said they would.** The KA-4 repair (per-construct
centring on the construct's own wrong-recipient distribution) works AGAINST the effect and moved the
estimate from **+0.78797** (uncentred on the repaired pools) to **+0.71293**, **-9.5 %**, as
predicted; it did not narrow the null (repaired null mean +0.00134, sd **0.17316** vs C8's +0.00005,
0.17173), so the p rose from the floor 1.0e-5 to 2.0e-5. The KA-6 repair moved it **+0.0023**
(+0.3 %). Decomposition: C8 as published +0.78567; KA-6 repair only +0.78797; KA-4 repair only
+0.70998; both +0.71293. **Variants (a) and (b) have the SAME point estimate by linearity**
(mean(z) - mean(b) = mean(z - b)), so the variant choice could not have been made for its effect on
the answer; they differ only in the interval — (a) [+0.32419, +1.07351], (b) [+0.18217, +1.23115].
The preregistered justification for (a) is confirmed by measurement: the per-construct
wrong-recipient level `b_c` has mean +0.07504 but **sd 0.42538**, range [-0.8093, +1.2951] — a
spread **5.7x** its own mean, so a scalar correction would have left it in place.

**Known-answer rows.** KA-1 max |diff| **0.000e+00**; KA-2(a) 4/13, 4/21, 4/31, 8 L=12 constructs,
83 pairs; **KA-2(b) the inherited 209,070-spacer pool is the CORRECTED one** — 0 Excel-escaped
labels, **0** own-target-in-own-hit-set violations at all four L, pool L=12 rate 0.160879 and
real-spacer rate 0.132353 both at |diff| 0.0e+00 vs C8's saved values, bug-equivalent rate **0.597**;
KA-3 identical NTC sums, 38,176 control cells; **KA-4 (unchanged machinery) +0.07168 (SE 0.00222)
reproducing C8's +0.0717** — it still fails C8's original <= 0.05 threshold, which is exactly the
defect repaired, and it governs here by reproduction; KA-5 KS p 0.1794 / 0.2456 on the repaired and
C8 nulls; **KA-6 (unchanged) reproduces C8's nine SMDs at 0.00e+00 including n_windows +0.2905**;
KA-7 2.220e-16; **KA-4R +0.00290 [-0.30574, +0.31191]** against C8's original rule; **KA-6R every
axis <= 0.25** (`log1p_n_windows` +0.2905 -> **0.0000**; worst axis now `map_mult` **+0.2460**,
inside the rule with 1.6 % to spare). C8's own primary was reproduced from C8's pools at
**|diff| 0.0e+00** (+0.78567, p 9.9999e-06).

**Attainable resolution at the finer stratification, no construct dropped.** 83 pools, mean size
180.4, minimum 25; **3 of 83** needed the declared widening and are named (`DPF2_P2`/`ITIH3` cell
n=18, `ZNF581_P2`/`RUNX2` n=9, `STAG2_P1`/`MAP3K7CL` n=10). All 34 constructs retain a full candidate
set. The one repair element whose direction FAVOURS the effect — the global L*-candidate exclusion —
is worth **+0.00034** (0.05 %), measured with it switched off.

**Secondaries (no verdict).** Real-spacer resampling null, all 34 constructs resolvable, 672 joint
draws: uncentred W +0.74441 p **0.00149**, centred W +0.67242 p **0.00149**, both at the floor.
Partner-construct specificity on the repaired pools: own +0.78797 vs partner **+0.04296
[-0.29799, +0.38107]**, own - partner +0.74500 [+0.22896, +1.25766], paired Wilcoxon p 0.00459,
n = 34 (all 34 partners are in the same split). 26/34 constructs positive, exact binomial p 0.00147.
Above the CP50k 0.5 floor only, S +0.65229. **The wrong-recipient permutation (p 2.0e-4 at its floor)
is under repair (a) the SAME OBJECT as KA-4R and is NOT an independent check** — a real cost of the
repair, which is why the sequence-resampling null carries more of the corroboration burden in C9.
**T7 is still unresolved**: failing (13) +0.8333 vs working (21) +0.6384, difference +0.1949
[-0.5584, +0.9434], Mann-Whitney p 0.7231 — the sign is opposite to C8's and the interval is seven
times the difference. **The effect is still not concentrated at L = 12** (L* = 11 n 16 +0.963
p 4.0e-5; L* = 12 n 8 +0.685 p 0.030; L* = 10 n 7 +0.197; L* = 9 n 3 +0.658). The 31 non-targeting
constructs, **candidate-bearing observations selected by the same rule and NOT negative controls**,
read S +0.43522 [+0.03554, +0.82811] p 0.00666 — 61 % of the targeting aggregate; on those
constructs intended-target biology cannot be the explanation at all, so this is corroborating, not a
control failure. Sens-2 (centring across all 33, cross-split) reads +0.78480: the preregistered
within-split centring is the more conservative choice and is the one that governs.

**Audit finding in C8's write-up (KA-0).** The 34 targeting constructs split **train 16 / test 10 /
validation 8**, not the "train 21 / validation 7 / test 6" stated in C8's `RESULTS.md` §1 and
`PREREGISTRATION.md` §2. C8's own `constructs307.csv` — the object every C8 computation read — has
16/10/8, so **no C8 number is affected**; the sentence is wrong.

**What this does and does not establish.** A PASS says genes selected by a 12-nt-seed + NGG promoter
match are, in aggregate, more repressed in the cells carrying that construct than matched genes the
rule did not select, **over and above the level those same genes read in the other constructs of the
split**. It does **NOT** establish dCas9-KRAB off-target binding: expression here is observational,
nothing measured is an occupancy assay, and the partner-construct test narrows the explanation to
something carried by the protospacer sequence **without separating seed-mediated binding from a
construct artefact correlated with spacer identity**. **C9 is a second test of the same hypothesis on
the same data — not an independent replication.** Its point estimate was predictable in advance from
C8's own post-hocs (P5 +0.71399, P2 +0.78813), the preregistration says so, and the honest summary is
the interval plus the fact that the sign and magnitude survive a calibration repair that works
against them. The decisive next design is unchanged: a **seed-mutant control** (~24 re-cloned
constructs, one 10x lane, order $3-6 k); for binding itself, a dCas9 occupancy assay ($5-10 k); with
the construct as the unit, more constructs narrow the interval and more cells do not.

No model, no submission, no production change. See
`experiments/309-c9-calibration-repair/RESULTS.md` and `work/enrichment309.csv`.

## Reserve R6 — DECLARED 2026-09-29T13:57Z by 534 (session "Public data scope"), before any draw

- **What:** four fresh Arc-like panels, S1 to S4.
  - Drawn by 238's `bundle238.py` rules: same gene axis, relabelled A/B/C controls, essential list subtracted.
  - One change: the pool also excludes every target of every `pert_counts.csv` on disk at draw time (7 distinct
    panels at declaration) and each earlier R6 draw. So R6 is disjoint from every tuned-on panel and internally
    disjoint.
  - Seeds: S1 534101 (efr rule), S2 534201 (efp rule), S3 534102 (efr rule), S4 534202 (efp rule).
  - The draws are written with their sha256 to `experiments/534-generalization-audit/work/partB/draws.json` by
    `draw534.py`.
- **Preregistration:** `experiments/534-generalization-audit/PREREGISTRATION-B.md`, sha256
  `74f108cf17cfd05972fac9517daa3449e83d90702b1108e68823aea941d5fd8c`.
- **Access rule:** two claims only.
  - C10a: the October bulk recipe's cumulative milestones M0 → M6 plus two controls, priced once now.
  - C10b: the final October recipe's exact bulk configuration, priced once, after the 10-05 A/B/C read.
  - No fitting, selection, dose search or debugging may read an R6 score. The known answers run on `efp`.

### Claim C10a (preregistered 2026-09-29T13:53Z, before the draw): does the October bulk recipe's gain hold on never-used Arc-like panels?

- Designs, statistics and the frozen reading rules: PREREGISTRATION-B §4–§7. The step verdicts are TRANSFERS / FRAGILE /
  FAILS on the R6 pooled paired-target bootstrap (4,000 draws, seed 534). Also recorded: retention vs A/B/C, the
  tuned-on check against `real`/`efp`/`efr`, and the controls (CD4 ×1.5, a permuted key) ≤ 0.
- The outcome is appended below after the single access, whatever it is.

### Outcome C10a (the single R6 access, 2026-09-29, by the 534 Part B worker): the October bulk recipe keeps 70 % of its A/B/C gain on never-used panels; every milestone TRANSFERS; the LODO key is the weak switch

*Appended 2026-09-29T17:41Z by the owner of 534 after checking the headline numbers against `results534B.json`. Numbers from
`experiments/534-generalization-audit/work/partB/results534B.json`; write-up in `experiments/534-generalization-audit/RESULTS-B.md`.*

- **Access.** R6 was drawn once by `draw534.py` at 2026-09-29T15:50:51Z, after the efp known answers had passed. The draws
  and their sha256 are in `work/partB/draws.json`:

  | panel | rule | seed | pool (essential-free) | sha256 prefix | K562 / HCT116 / HEK293T |
  |---|---|---|---|---|---|
  | S1 | efr | 534101 | 5,599 | 93237a3a | 300 / 300 / 300 |
  | S2 | efp | 534201 | 14,016 | 835d4b32 | 137 / 291 / 293 |
  | S3 | efr | 534102 | 5,184 | d9112014 | 300 / 300 / 300 |
  | S4 | efp | 534202 | 13,416 | 88a7802a | 126 / 293 / 293 |

  - Every panel has 300 targets, 0 essential, and 0 overlap with the 9 panel files on disk (7 distinct panels, 1,829
    excluded targets) or with each other.
  - All four panels were priced: 12 held-out cancer folds, 16 fresh operators. Each operator passed `_check_axis` on its
    own panel, and its axis differs from A/B/C's and from all 16 tuned-on (panel × donor set) axes.
  - The scores were read once, by `aggregate534.py`: 4,000 draws, seed 534. No fit, selection, dose or debugging step
    read an R6 score.
- **Known answers (efp, before the draw): all passed with |d| = 0.**
  - 438's stored rows reproduce exactly: M0–M3 and C-a, design and realised, both levelled twins, 3/3 folds. So do
    493's M4 row, 546o's M5 row and 522's `pctPROD` row for M6.
  - The fold caches are bit-identical, and so are the full-panel energies and level factors.
  - 72 prices from 438's own pricer match exactly; the H1 ×2 half on efp is +0.002545904997, as recorded.
- **Result.** 438's levelled route, k = 1, pooled over R6's 12 folds, 95 % paired target-bootstrap intervals. The
  A/B/C denominators are the record's values.

  | step | R6 pooled [95 % CI] (folds positive) | verdict | retention vs A/B/C |
  |---|---|---|---|
  | M0→M1, q 1.25 | +0.00652 [+0.00593, +0.00727] (12/12) | TRANSFERS | 0.55 |
  | M1→M2, H1 ×2 | +0.00099 [+0.00037, +0.00150] (8/12) | TRANSFERS | 0.42 |
  | M2→M3, q 1.5 | +0.00458 [+0.00403, +0.00526] (12/12) | TRANSFERS | 0.53 |
  | M3→M4, q 1.75 | +0.00282 [+0.00232, +0.00342] (12/12) | TRANSFERS | 0.61 |
  | M4→M5, LODO key | +0.00099 [+0.00033, +0.00181] (7/12) | TRANSFERS | 0.30 [0.10, 0.55] |
  | M5→M6, percentile cap | +0.00694 [+0.00517, +0.00828] (10/12) | TRANSFERS | 3.31 |

  - **Cumulative M0 → M6: +0.0229 [+0.0204, +0.0255] against +0.0327 on A/B/C, retention 0.70 [0.62, 0.78].** Against
    the exact designs priced on A/B/C's own panel (supplementary) it is 0.68 [0.61, 0.76].
  - By panel: S1 +0.0305, S2 +0.0250, S3 +0.0221, S4 +0.0138.
  - Calibrated board-unit projection: +0.0182 [+0.0157, +0.0206], against +0.0234 for the same designs on A/B/C.
- **Rules (PREREGISTRATION-B §7).**
  - No switch FAILS or is FRAGILE, so nothing is a switch-off candidate.
  - Tuned-on check: 1 of 3 October switches reads below half of its real/efp/efr mean. The LODO key reads +0.00099
    against half-mean 0.00147; q 1.75 and the percentile cap do not. So real/efp/efr are **not over-read** by the rule.
  - Controls:
    - C-a (CD4 ×1.5 on M4): −0.00087 [−0.00165, −0.00006]. 438's CD4 ×1 decision stands.
    - C-b (the LODO key permuted over targets): −0.0102 [−0.0116, −0.0085], 0/12 folds.
  - Surprise rule: no sign flip on M3–M6, and cumulative retention is above 0.3, so this is not a surprise.
- **For the coordinator (information only, nothing changes).**
  - The LODO key is the weakest October switch off A/B/C.
    - It is the one switch whose R6 value sits below half its tuned-on value.
    - It is negative on S4: −0.00060, 0/3, relative cosine −1.9 %.
  - The percentile cap carries most of R6's October value: +2.62 % relative cosine on A/B/C, +18.9 % on R6. It is
    ≈ 0 on S4, the panel with the weakest donor transfer (ret 0.21–0.30).
- **Deviation.** The owner set the gate floor to 27 GiB (`work/partB/GATE-OVERRIDE.json`, `swapguard534.sh`), because
  §9's 33 GiB free+inactive was unreachable (25.3–29.3 GiB) with the desktop apps open. No design, seed or reading rule
  moved.
- **C10b remains:** one more access, after the 10-05 read, for the final October recipe's exact bulk configuration.
  The four R6 operators per panel were deleted after pricing to free ~19 GB. They are deterministic re-fits (seed 0), and
  their sha256 values and exact commands are in `work/partB/state534.json` for C10b.

### Claim C10b — timing amendment, proposed 2026-09-30T21:23Z before any C10b access (owner: session "Explore contribution areas", experiment 551; agreed with the coordinator)

- **Proposal: C10b runs once, after the 10-04 A/B/C read, instead of after the 10-05 read** (PREREGISTRATION-B §3). The
  coordinator names the exact bulk configuration in a message after that read (~04:30Z); it is recorded verbatim here, with
  its timestamp, before the access. C10b runs 10-04 daytime.
- **Reason.** R6 prices the bulk channel only (pds, mse). The bulk switches still pending are all read by 10-04: `sharpen_pt`
  (10-02, 479's pair), `lodo_hi_pct` and the union source (10-03, 495 / 471c), the combination (10-04, 498⟨P⟩ / 497dr). The
  10-05 read (496a, `mde_sign`) is DE-only and cannot change what C10b prices. Running on 10-04 also keeps the heavy operator
  re-fits off 10-05, which is reserved for recipe_dr524, the October bundle, its inert known answer and a rehearsal.
- **Unchanged.** One access, one claim. No fitting, selection, dose search or debugging reads an R6 score. Known answers run on
  `efp` before any R6 number exists. The preregistration (`experiments/551-r6-c10b/PREREGISTRATION.md`: design tree, statistics,
  reading rules) is hashed and its sha256 recorded here before the access.
- **R6's state at this proposal.** One access so far (C10a: drawn 2026-09-29T15:50:51Z, priced and read once by
  `aggregate534.py`); nothing has read an R6 score since. The 16 fold operators are deleted; their sha256 values and exact
  commands are in `experiments/534-generalization-audit/work/partB/state534.json`. The fold caches stay in `…/partB/S1`–`S4`.
- **Scope.** A bulk switch adopted after this access (for example 459dr's H1 ×3 on 10-07, or 522's strength rider) is not
  covered by C10b. It goes to R7. If the read that fixes the last pending bulk switch is missing or late, C10b waits for it,
  and it does not run after 2026-10-06 12:00Z without a new amendment here.
- **Known answers (not an access), 2026-10-01 00:32Z:** 551's pricer reproduces every 534 / 438 / 493h / 546o / 522 row on
  `efp` with |d| = 0.0 on all three folds, including the per-target forms (`experiments/551-r6-c10b/RESULTS.md`).

- **F recorded 2026-10-04T05:40:10Z** (coordinator's message received 2026-10-04T05:39Z): `M6+P2g` (base M6, step P2g@M6). `experiments/551-r6-c10b/work/c10b-config.json` sha256 `113fa23d29e93b60d3ef5172e6932dcb1cf2f55e80f73cdba5e29d9ed33fc868`; PREREGISTRATION.md sha256 `f4d0b70bc006a19c3b99c71639947c1f0180eb43d51f2071ab98fe3c67377c6f`.
- **C10b outcome 2026-10-04T07:19:15Z** (`drive551.py access` 05:40:41–07:19:15Z, F = `M6+P2g`). **KA-REPRO exact:** all 16 R6 operators re-fitted to C10a's sha256 and axis; every fold cache bit-identical to C10a's; 39 C10a rows per fold and C10a's 8 pooled steps reproduced with |d| = 0.0. **P2g@M6 (per-target sharpening, the carried-only key tercile, on the percentile-cap base) TRANSFERS: +0.00208 [+0.00132, +0.00287] overall-equivalent, 11/12 folds > 0** (S1 +0.00030, S2 +0.00289, S3 +0.00277, S4 +0.00237). Calibrated board-unit projection +0.00093 [+0.00061, +0.00128]; the board read the same step at +0.0010 on cap 3 (479p2g) and on cap 6 (498p2g). Control (q permuted over targets) +0.00107 [+0.00019, +0.00202]; targeting (claim − control, paired) +0.00102 [+0.00010, +0.00190]. Cumulative M0 → F +0.02493 [+0.02207, +0.02797] (C10a's M0 → M6 +0.02285). Levelling sensitivity +0.00210. Consequence: nothing changes. **R6 is now spent (C10a + C10b) and closed.** `experiments/551-r6-c10b/work/results551.json` sha256 `edc4b9baba0e2565570c5286ea9368f7fbb5c8bb21efad64df5eca4d670b7995`, `tables551.md`.
## Reserve R7 — DECLARED 2026-09-30T21:43Z by 552 (session "Explore contribution areas", at the coordinator's request), before any draw or build

- **What:** R7a, four Arc-like bulk panels T1–T4 (534's draw rules; seeds 552101 / 552201 / 552102 / 552202; the exclusion
  now includes R6's four panels), and R7b, one cell-level DE surface in a cancer line, Sunshine et al. 2023 Calu-3 CRISPRi
  (`data/external/scperturb/SunshineHein2023.h5ad`, sha256 `4dd14fb68dc16057e1cf4a05d0bcd2efd26c1a9e4ef8c7d21afccaf2cfc84270`;
  its cells are SARS-CoV-2-infected, a stated caveat). Its prior touches (225's bulk donor sweep, 284/290's conservation maps,
  the same-line registry) are disclosed; no DE design has used it.
- **Declaration:** `experiments/552-r7-reserve/RESERVE.md`, sha256 `682395abaa0c98d9a15e1af4d3286d12fff1133bc0c9e79d19ebafb4848bcffe`.
- **Access rule:** one claim, once: the final October recipe (bulk on R7a, DE on R7b) against the recipe before the switches
  decided after C10a, after R6's last access (C10b) and after the recipe's DE switches are decided. The claim is
  preregistered and its sha256 recorded here before the first R7 number exists. No fitting, selection, dose search or
  debugging may read an R7 score. Void for R7b if an Oct 22 context matches Calu-3 before the access.
- **State:** declared; nothing drawn, built or read yet.
- **R7a drawn 2026-09-30T21:44Z** by `draw552.py` (5.8 s, after the declaration commit b4dc26f9): `work/draws.json` sha256
  `418bb6058736aa0d23cf3a5313bd8942390bee7257c3a9d0954cca21c39622bb`. Exclusion: 13 panel files (11 distinct panels, R6's four
  included), 3,029 targets. T1 efr 552101 pool 4,774, sha `2809543a…`; T2 efp 552201 pool 12,816, `184c0241…`; T3 efr 552102
  pool 4,376, `4bc68e21…`; T4 efp 552202 pool 12,216, `9e38b94f…`. Every panel 300 targets, 0 essential, 0 overlap with any
  panel file, R6 or each other; coverage K562 / HCT116 / HEK293T 300/300/300 (T1, T3), 114/290/293 (T2), 113/292/296 (T4).
  Controls byte-identical to 199's relabelled files, stored as clones. No operator, fold or score exists.
- **R7b:** not built yet (construction rules in RESERVE.md §3; metadata only until then).
- **Amendment 1 (2026-09-30 ~21:55Z, before any R7b build):** RESERVE.md §5, now sha256 `230c359116d76a9ab1e92eef2c702b805dce742dfd2c447bbbf73547ddcc05fe`. (1) The Sunshine file has
  no mock condition and no viral features (all 33,538 features human ENSG; metadata only), so the infected state stands.
  (2) Asymmetric rule: R7b can veto a DE switch only on a preregistered HARM beyond its own re-roll noise band; a null or
  small effect never counts against a switch; a positive is robustness evidence only. (3) Built this week outside the
  coordinator's windows; accessed only after the DE switches settle (not before the 10-07 reads).
- **Amendment 2 (2026-09-30 ~22:00Z, before the surface is built):** RESERVE.md §6, now sha256 `f3b8268f396852c619c1632f444cf7541f3a2c258b82503f1e36a5ac7d681205`. Control-only
  statistics, the target list, per-target cell counts and emitted predictions are prediction-side inputs (the October
  day computes the same from D/E/F controls), not reads. Any statistic involving a perturbed cell's expression, the
  cell-eval2 reference bundle and every score stay sealed until the access.
- **R7b surface built 2026-10-01T00:27Z** by `build_r7b.py` (construction only: counts copied, no statistic; 13 s,
  peak 6.7 GB): 107 targets, 27,396 cells (3,536 controls), cells per target 103–826 (median
  187); official genes mapped 18,333 / 18,533. `work/r7b/build_r7b.json` sha256 `6fd6e405ecd0689c2f95adad74ad8597a284147d2e138d2b79f84caadc5b3472`; truth
  `truth_r7b.h5ad` `bd35f450f24f988b40ba4f7a940ad74810a800e190f36c920419afb2446d1453`; controls `bundle/context_R.h5ad` `1506e016f891746c28729eba8ec1f40ed540876977c14b4390baa06e06f6cd35`.
  The truth file is opened again only by the access script.
- **Amendment 3 (2026-10-01 ~08:45Z, before any target cell's library size is read):** RESERVE.md §7, now sha256
  `b5193eb8c748dbfc70209f0519dba08666237c6a9f0233346b47d5952f65de99`. A human-library floor of 5,000 official-axis counts for controls and target cells (553: the production emitter
  cannot place R7b's bimodal libraries, q75/q25 7.24), BST2 excluded (emitter edge case), targets with < 30 floored truth
  cells dropped. Reads one library-size scalar per target cell, for the fixed threshold only. Rebuilt from the source into
  `work/r7b_floor5000/`; the unfloored build stays as built.
- **Floored R7b built 2026-10-01T08:38Z** (`build_r7b.py --floor 5000 --exclude-targets BST2 --min-cells-after 30`, from the
  source, 13 s, peak 6.9 GB): 106 targets (none dropped by the 30-cell rule), 18,759 cells, controls
  3,536 → 2,378 (equal to 553's `r7bd5000/context_R_5000.h5ad` cell for cell and count for
  count), floored truth cells per target 59–556 (median 127). `work/r7b_floor5000/build_r7b.json`
  sha256 `b17112aeaf2936a321673e0ed8fb38ef6e985d71d1a242202b8ae7c3f1d2c9f3`; truth `09435d8f924f95ab6a1bf89f20feda5d4e67394a05837393e23746b7e144bb83`. This floored surface is the one R7b's access reads.
- **R7 claim recorded 2026-10-14T23:41:44Z**: `experiments/552-r7-reserve/CLAIM-R7.md` sha256 `f134a63df5c2021cd9a03ff81e90311e52b81539d029a229750cc5022d0987bc` (F `M6+P2g`; DE block off, hook on, K default). Filled by `fill_claim552.py` from docs/october-decisions.json (lodo_hi_pct 10-03T10:02:43Z, sharpen_pt tercilec 10-05T23:19:45Z, mde_sign clean4:0.05 10-05T10:36:41Z, callset_block off 10-14T22:28:57Z, rc_callset_cpm off 10-06T18:30:57Z). Accessed early on the coordinator [400b6e]'s approval (10-14 ~23:42Z) by the science session [bf003b] (the reserve's owner session is gone). Out of the claim's grammar, so NOT confirmed by R7: de_uprn (630), de_dnmag (631, reads 10-15), PIE bulk (702).
- **R7a armed 2026-10-14T23:45:06Z**: `experiments/551-r6-c10b/work/r7a/claim.json` sha256 `4b160fb6ad64d6221bfb372f53516f8c7443be2688c5392bc129abddeceec369` (F `M6+P2g`, window 2026-10-14T23:45Z - 2026-10-16T12:00Z). Independent claim check by [17adab] (10-14 ~23:44Z): the claim encodes the decided recipe; the not-read statement for UPRN/DNMAG/PIE sits in the line above (written before any access), an accepted deviation from the draft's wording.

- 2026-10-16T04:36Z R7a ADMISSION-ONLY DEVIATION, recorded before any access ([d80a10], R7 owner): gate551.py pinned at blob da4f9878 (1e7c26b2: vm.swapusage total 0 = pass; this boot has no swap file). r7a551 allows exactly that blob; the 30 other instrument files are byte-identical to e7be21ca. A 04:33:14Z start had failed in plan() before access (no truth opened). Details: experiments/552-r7-reserve/RESULTS.md.

- 2026-10-16T04:38Z R7a WINDOW END EXTENDED TO 2026-10-16T18:00Z (coordinator [2d6f46]'s decision 10-16 ~04:40Z: the window protects the order of work, not the result; the reserve stays single-access). Recorded before the wrapper reads it. claim.json is untouched (its sha stays the one on the existing arming line); r7a551's window_end() reads this line.

- 2026-10-16T07:13Z R7a DRIVER PINNED AT BLOB a96203c241f5404eebf3b7afbabe751961eb4cb8 ([d80a10], R7 owner; recorded before any access resumes). The same admission-only deviation, carried into access()'s comparison with the record: the 06:34:06Z resume failed at "the instrument changed since the record" (the record holds gate551.py's efp-KA blob f36f5198) before touching any panel, and the record's driver blob (4b7f640d, at the arming) would have failed the next check too, because the two 10-16 changes above live in r7a551.py. access() now accepts gate551.py only at the pinned admission blob da4f9878 (with the record at its KA blob), and this driver only at its recorded blob or at the blob this line pins; every other instrument blob, every package version and every other external script must equal the record. No arithmetic changed. Proof before the resume: `experiments/552-r7-reserve/test_access552.py` runs access() to the truth boundary (panel() replaced by a sentinel): FAIL without this line (driver check), PASS with it.

- 2026-10-16T10:34Z R7a ADMISSION-ONLY DEVIATION 2, recorded before the access resumes ([d80a10], R7 owner; coordinator [2d6f46] OK 10:33Z, guard-only, same class as the first): T4:comp551 died at 10:26:42Z on 062's peak-RSS self-guard ("cd4_rest centering: peak RSS 8.430 GiB reached 8.000 GiB"; sens415:92 hard-codes 8.0 GiB; ~25 GiB free+inactive, so not memory pressure; R-NEW hit the same line this morning, lifted by readiness in def_lodo_key, 690j). T4's comp551 wrote nothing; nothing is aggregated or read. comp551 now runs through `experiments/552-r7-reserve/rss_lift552.py` (blob 884bd6e3), which replaces exactly one attribute of one module, `e162_062.enforce_memory_limit` (the module plolo registers from 062-cd4-post-whiten/run.py), by 690j's max(cap, 16 GiB) guard, in that process only. Why it is inert: 062/run.py lines 48-52 (`enforce_memory_limit`) are a pure check-and-raise, and `max_rss_gib` reaches only that function (calls at lines 189, 211, 224, 527, 548, 773; forwarded elsewhere; the CLI arg check at 1056); it never sizes chunks or changes any array. sens415.py, 062/run.py and every other pinned instrument file stay byte-identical. r7a551.py asserts the shim's blob, so the pin below covers both. Conditions (coordinator): the verdict stays the preregistered 12-fold statistic; nothing from T1-T3 is aggregated or opened before T4's price exists; if T4 fails again, nothing is aggregated and R7a is reported INCOMPLETE. Proof before the resume: the shim passes a cap below its process's peak that raises without it, on exactly 1 module; test_access552 FAILS without the pin below (driver check) and PASSES with it.

- 2026-10-16T10:34Z R7a DRIVER PINNED AT BLOB 423b8bfe9200c0631f9906b4a39f8715aa403361 ([d80a10]; supersedes a96203c2 for the resume after DEVIATION 2; access() takes the newest pin).
- **2026-10-16T11:58Z R7a VERDICT: TRANSFERS** (appended by the coordinator [2d6f46]; the access ran 10-16 09:14:35Z → 11:53:39Z, "END rc=0", one aggregate;
  R7 owner [d80a10]). PRIMARY F = M6+P2g vs M5, 12 folds on T1–T4: **+0.01056 [+0.00768, +0.01193], 12/12 > 0**, calibrated +0.01215;
  percentile cap +0.00839 [+0.00560, +0.00950] (12/12); P2g +0.00217 [+0.00144, +0.00300] (10/12), targeting over its permuted control
  +0.00349 [+0.00217, +0.00499]; retention vs R6 1.17 [0.85, 1.32]. Consequence by the claim: nothing changes. Caveat: T1–T4 are
  fresh-drawn panels, which over-read key steps ~×3 against Arc-drawn ones (FINDINGS §1, lab 569978b9), so this confirms sign and
  transfer to unseen targets, not a board-sized gain. Record: `experiments/551-r6-c10b/work/r7a/results_r7a.json`, `tables_r7a.md`,
  `experiments/552-r7-reserve/RESULTS.md` (4870e968). R7a is now SPENT. With R7b read on 10-14 (NO VETO, 550281d6), the whole R7 reserve is spent; the next sealed reserve is the lab's LZ.

## Reserve LZ — DECLARED 2026-10-15T23:06Z by the lab lead ([c61d4c], experiments/800-checker), appended by the coordinator ([2d6f46]) 2026-10-16
- Two Arc-like 300-target panels drawn by the lab's rule "arc800" (`experiments/800-checker/DESIGN.md` §3.1; `experiments/801-lab-folds/draw801.py`
  sha256 c4b8507e8f4c6a0f02069ec5766b2ee0b62452494a2df16b2a5246fc7f2a7ead) in commit e4cca4b2 (2026-10-15 23:03Z), after the six WORKING panels
  L1–L6, excluding every panel on disk at draw time (A/B/C, real, efp, efr, R6 S1–S4, R7a T1–T4) and each earlier lab draw:
  - **L7:** seed 800107, target-list sha256 97c30e7d11dddc0157219cd62a4312001553a7931b83bce8b3f087ebfddc9619, file
    `experiments/801-lab-folds/work/sealed/L7/pert_counts.csv` (sha256 33bb3256…, full value in SEALED.md).
  - **L8:** seed 800108, target-list sha256 32fd27d28425b12e102f7c15f16d01b95b10e32da9b37f9312d4df52c50694f1, file
    `experiments/801-lab-folds/work/sealed/L8/pert_counts.csv` (sha256 c8d9f1a3…, full value in SEALED.md).
- Held-out folds when opened: K562-GWPS, Orion HCT116, Orion HEK293T (PLOLO, 801's bulk engine) = 6 folds. Nothing has been built on them.
- **Access rule:** one access in total, by the lab lead, for one preregistered final candidate after its WORKING price exists; the claim file's sha
  goes into this ledger (through the coordinator) BEFORE the first operator fit; the result is appended here whatever it is. No fitting, selection,
  dose search or debugging may read an LZ score. Full text: `experiments/800-checker/SEALED.md`.
- **2026-10-20T09:20Z: LZ access CLAIM recorded (coordinator [3c7a7e]), before any LZ bundle, operator, dump or score:** `experiments/800-checker/LZ-CLAIM-811.md`,
  commit 1ddf4664, **sha256 3a8f725216075e793a9968b53f2ff0f28f845593ba49e48d3c659f009bd258c2** (verified against the committed, clean file at 09:19Z).
  - Candidate: 811's c0_l25 (frozen since 10-16). Folds L7/L8 × K562 / HCT116 / HEK293T = 6. Composite B0/B1/C0/C1; both increments preregistered
    (C1 − B1 swap ON, C0 − B0 swap OFF). check800 frozen, pooled, N = 1. PASS iff the deciding increment's pooled point > 0; FAIL ends 811's 6a dossier.
    Opened only after DEFo passes runbook §5 items 1–4.
  - **Coordinator's fixed mapping (SIXA-DECISION Amendment 4, 08ec9597, written before this claim existed):** "the branch production ships" means
    the CONSTRUCTION X2 ships. On a swap-OFF day: C0 − B0. On a swap-ON day: C1 − B1, unless Amendment 4 picked C0 from DEFo's C1 − C0, in which case
    C0 − B0. The construction is fixed from DEFo BEFORE LZ is opened, so this is still one read under a rule set in advance, with no new freedom (both
    increments are in the claim).
- **2026-10-22T07:49Z: LZ ACCESS RESULT (recorded by the coordinator [3c7a7e]); LZ is now SPENT (one access).** Claim 3a8f7252 (1ddf4664); access record
  experiments/800-checker/ACCESS-LZ-811.md (sha256 e967727d…), json experiments/811-lolo-model/work/run-comp811-LZ/agg811-811-LZ-6048cf1a.json (2e446425…), lab commit 734cca22.
  - Deciding increment (fixed from DEFo before LZ opened, SIXA Amendment 4b → X2 = C0): **C0 − B0 pooled over L7 + L8 × 3 held lines: +0.00999 [+0.00659, +0.01273], 6/6 folds,
    LT3 +0.00571 → PASS** (pooled point > 0). Reported, not decided: C1 − B1 +0.00645 [+0.00373, +0.00882], 6/6.
  - Every composite KA passed before any score was read (KA-705/707 on all 6 LZ folds |d| 0.0; KA-811; KA-L0 10/10 per run). Essential-screen counts L7 9, L8 20, so the arc800
    draws are expected to over-read (LZ ≈ 2.2× DEFo). Override: wrapper lzaccess811.sh e990d562…; no guard or engine edited. Nothing else was read from LZ.

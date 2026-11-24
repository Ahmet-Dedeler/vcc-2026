# 350 — Same-line priors: Jurkat bank measured, common axis is generic, HeLa/CAL33 have nothing

**Status: COMPLETE (descriptive inventory + light computation; no arm, no score claim).**

One-line verdict: the Nadig Jurkat bank is real and usable (2,393 targets, 0 on panel,
8,284 genes on the challenge axis), but its common response vector is a *generic* CRISPRi
axis — cosine 0.331 with HepG2's, 0.349 with K562's, −0.011 with our H1 design — so the
PLAN kill criterion **fired**: there is no line-specific same-line direction here to transfer.
Bonus finding: that generic axis is genuinely present in H1 *truth* (per-target cosine
0.108, positive on all 85 targets) and entirely absent from our whitened *design* (0.007) —
whitening removed it by construction. HeLa has no public knockdown-transcriptomic bank at
all; CAL33/HNSCC has none either. Jurkat is the only context with same-line data, and what
it offers is calibration-map targets, not a line-specific direction.

## Known-answer reproduction (own code, before any new measurement)

- Mean of `per-target-cap15.csv` rows with `metric == 'de_wilcoxon_lfc_nmae'` (58 rows)
  = **0.9508364630325662**, bit-identical (`==`) to `metrics-cap15.json` raw nmae and to
  343's recomputed value. **Match.**
- `effects.npz['cap15']` shape **(85, 18533)** as recorded. **Match.**

## What was measured (inputs, estimator, controls)

- **Inventory** (`work/compute350.py::load_headers`, h5py headers only): X shape/dtype,
  obs/var columns, `obs/__categories/gene` cardinality, control cell counts from
  `obs/gene` codes, gene-id type. Code path reused from `scripts/build_donor_table.py`
  (`label_col='gene'`, `gene_col='gene_name'`, `control='non-targeting'`). Byte-decode
  caveat found and handled: `var/__categories/gene_name` stores **bytes**.
- **Common vectors**: seeded-random 300 targets per bank (`default_rng(350)`), restricted
  to ≥20-cell targets (mirrors donor builder); X streamed in 8,192-row dense blocks,
  float64 group sums; CPM-of-means; **lfc = log2((CPM_pert+1e-9)/(CPM_ctrl+1e-9))** —
  exactly `exact_de323.de_from_counts`'s formula. Duplicate-symbol collapse (one dup per
  bank) mirrors `build_donor_table.py`. Common vector = mean over 300 targets; mapped to
  the 18,533 challenge axis. DE propensity per gene = fraction of 300 targets with
  |lfc|>0.5; Spearman vs context-A mean basal CPM (backed CSR, chunked).
- **Controls**: HepG2 as the non-A/B/C control line; H1 *truth* common vector from
  `reference.npz` sums (same log2 space, 85 targets) to separate "design misses axis"
  from "space-mismatch artefact"; K562-gwps donor-bank mean (`delta.npy`, log1p-space —
  caveated comparison); top-1% magnitude-drop robustness; library-depth check
  (controls vs perts mean/median UMI from `obs/UMI_count`, no X load).
- Outputs: `work/jurkat_common.npz`, `work/hepg2_common.npz` (challenge-axis `genes`,
  `vector` = mean log2-CPM-ratio, `mapped` mask, `propensity`, `targets`,
  `bank_genes`, `common_bank`), `work/inventory.json`, `work/compare.json`,
  `work/followup.json`.

## Results

### 1. Bank inventory (units: cells, genes, symbols)

| bank (file) | n_cells | X cols (raw) | unique symbols | obs target levels (excl NTC) | control cells | cells/target median [min, max] | UMI median ctrl / pert |
|---|---|---:|---:|---:|---:|---|---|
| Jurkat `nadig2025/jurkat_raw_singlecell.h5ad` (9.37 GB) | 262,956 | 8,882 | 8,881 | 2,393 | 12,013 | 85.5 [20, 2555] | 10,906 / 10,121 |
| HepG2 `nadig2025/hepg2_raw_singlecell.h5ad` (5.61 GB) | 145,473 | 9,624 | 9,623 | 2,393 | 4,976 | 48.0 [20, 536] | 16,020 / 15,630 |

- obs columns: `UMI_count, cell_barcode, gem_group, gene, gene_id, gene_transcript,
  mitopercent, sgID_AB, transcript, z_gemgroup_UMI` (+`__categories`); perturbation
  label column **`gene`** (categorical, 2,394 levels incl `non-targeting`); control
  label `non-targeting`. var columns include `gene_name` (symbol categories, bytes),
  `gene_id` (categorical Ensembl, 2,392 levels on obs side).
- Duplicate symbols collapsed: 1 per bank (same HSPA14-class pattern the donor builder notes).
- Depth check: controls only 6% (Jurkat) / 3% (HepG2) deeper than perturbed cells —
  **no depth confound** (contrast orion-hct116's 38%, `docs/HARNESS-DEPTH-CONFOUND.md`).
- Prebuilt donor tables exist: `data/donors/nadig_jurkat/` (2,270 perts ≥20 cells,
  log1p(TS=5e4) space) and `experiments/166-donor-acquisition/donors/nadig_hepg2/`
  (2,160 perts). Overlaps measured here, bytes-decoded:

| bank | panel∩targets (300) | var∩challenge-axis (18,533) | targets∩challenge-axis |
|---|---:|---:|---:|
| Jurkat | **0** | 8,284 / 8,881 | 2,176 / 2,393 |
| HepG2 | **0** | 9,024 / 9,623 | 2,176 / 2,393 |

Panel overlap 0 confirms 207 ("nadig_jurkat covers 0 of the panel") — every same-line
target is off-panel, i.e. calibration-map material only, never direct transfer.

### 2. Common vectors (units: log2 CPM-ratio; 300 random targets each, seed 350)

| quantity | Jurkat | HepG2 |
|---|---:|---:|
| entry-level fraction lfc<0 (300×G) | 0.495 | 0.525 |
| fraction of targets with mean lfc<0 | 0.713 | 0.893 |
| fraction of common-vector entries <0 | 0.618 | 0.714 |
| mapped genes on challenge axis | 8,284 | 9,024 |
| top-50 \|.\| share of vector energy | 0.280 | 0.325 |

Jurkat common extremes — most DOWN: CD40LG −3.40, TNFRSF4 −2.15, C9 −1.70, RORB −1.46,
PCOLCE2 −1.40, SCARB1 −1.34, AC010729.2 −1.21, BRMS1L −1.16, MNX1-AS1 −1.11, AC027601.6
−1.07, C6orf223 −1.04, LGALS1 −1.04, PIGB −1.03, SLC9A3R2 −0.97, MAN2C1 −0.96,
AC021683.3 −0.94, TXNRD2 −0.91, SLC38A5 −0.90, CNR2 −0.89, AL450326.1 −0.89, ACSL6
−0.89, TMEM143 −0.86, ZDHHC23 −0.86, CCNE1 −0.86, CDPF1 −0.83. Most UP: SNHG32 +0.37,
PBXIP1 +0.35, MALAT1 +0.33, PNRC1 +0.33, SAT2 +0.32, GAS5 +0.32, SLFN5 +0.32, H1F0
+0.31, TRAC +0.30, DDIT3 +0.30, BTG1 +0.29, SERINC1 +0.28, PPP1R15A +0.28, LIMD2
+0.27, BACH2 +0.27, SHISA2 +0.26, CD69 +0.26, LEF1-AS1 +0.26, ID2 +0.26, NBR2 +0.25,
GZMA +0.24, PNPLA8 +0.24, CDK19 +0.24, AL034397.3 +0.24, STK17B +0.24.
HepG2 common extremes — most DOWN: PDE4C −8.34, RRAD −5.59, AREG −4.02, LCN2 −3.93,
TUBA4B −3.77, EDN1 −3.36, CRYAB −3.26, UBD −3.22, SPINK1 −3.15, HEPN1 −3.14,
SERPINE2 −2.85, E2F2 −2.64, IGFBP2 −2.34, PLK4 −2.31, SPTLC3 −2.25, SPARC −2.16,
HIST1H1B −2.14, E2F8 −2.08, RAD54L −2.02, ADAMTS9 −1.99, IL18 −1.93, DDIAS −1.87,
SLC22A7 −1.87, ORC1 −1.86, C21orf58 −1.85. Most UP: CYGB +0.56, S100A6 +0.49,
CDKN1A +0.48, REEP2 +0.45, IGFBP1 +0.44, TP53INP2 +0.43, ATF3 +0.42, CAPN2 +0.42,
MOSPD1 +0.42, AC093866.1 +0.42, OPTN +0.41, GDF15 +0.41, NBPF1 +0.40, CDA +0.40,
MT2A +0.39, STOM +0.38, SMOX +0.38, SPINT1 +0.38, AC058791.1 +0.37, CPEB4 +0.37,
CHIC2 +0.36, SPOCK2 +0.36, TP53I3 +0.36, SERPINE1 +0.36, YPEL5 +0.35.
(Down tails reach −3 to −8 = 8–300×; per-gene control CPM was not re-streamed, so
low-count inflation of the extremes is undiagnosed — but §4's drop test shows no
conclusion rests on the tails.)

### 3. The three mandated cosines + same-space truth controls (units: cosine / sign agreement)

| pair (shared genes) | cosine | sign agreement (n) |
|---|---|---|
| **Jurkat common – HepG2 common** (7,200) | **0.331** | 0.658 (7,200) |
| **Jurkat common – H1 design (cap15 mean)** (8,262) | **−0.011** | 0.488 (8,262) |
| HepG2 common – H1 design (8,998) | −0.012 | 0.491 (8,998) |
| Jurkat common – K562-gwps bank mean¹ (6,908) | 0.349 | 0.669 (6,908) |
| HepG2 common – K562-gwps bank mean | 0.304 | — |
| H1 design – K562-gwps bank mean | −0.010 | — |
| H1 **truth** common – Jurkat common (8,172) | 0.134 | 0.564 (8,172) |
| H1 **truth** common – HepG2 common (8,903) | 0.125 | 0.569 (8,903) |
| H1 **truth** common – H1 design (17,383) | 0.004 | 0.517 (17,383) |

¹ Donor-bank means are log1p(TS·p) natural-log pseudobulk deltas vs our log2 CPM
ratios — direction comparison with stated space caveat, not a score.
H1 truth common (log2 space, norm ≈970, 64% down, diffuse: top-50 energy share 0.029
vs 0.28/0.33 for the banks).

- Per-target ceiling, same space, matched target order (verified identical 85-order):
  mean per-target cos(H1-truth-effect, JurkatCommon) = **0.1075 ± 0.0125, positive on
  85/85**; mean per-target cos(H1-truth-effect, our design) = **0.0068 ± 0.0127**
  (cross-space caveat on the design side: whitened log1p-deltas vs log2 truth).
- Robustness: dropping the top-1% |.| genes on either/both sides moves Jurkat–HepG2
  0.331 → 0.346/0.338/0.350 (broad-based, not tail-driven); truth–Jurkat 0.134 →
  0.103/0.065/0.056 (concentrated in large-|lfc| genes).

### 4. DE propensity: unusable as a key (units: fraction; Spearman)

- Spearman(propensity_Jurkat, context-A basal CPM) = **−0.563** (n=8,284, p≈0);
  HepG2 vs context-A = **−0.415** (n=9,024). High-basal genes are *less* often
  |lfc|>0.5 — the |lfc|>0.5 rate is dominated by low-count sampling noise, which
  crosses fixed thresholds most easily where counts are thinnest. A same-line
  "propensity" key from raw |lfc| rates would rank low-CPM genes first: backwards,
  and the same verdict as 344's null key by a different route.

### 5. What exists per context (web + local, accessed 2026-09-10)

| ctx | same-line perturbation transcriptomics | basal expression | what it can supply |
|---|---|---|---|
| **A = Jurkat** | **Nadig 2025** (GSE264667; local 9.37 GB h5ad measured above; GEO public): 2,393 targets, CRISPRi Perturb-seq. **GSE249595** (Song et al., Children's Nat'l; Jurkat E6 dCas9-KRAB, genome-wide CRISPRi CROP-seq, unstimulated + stimulated, 16 channels; 3.75 GB RAW tar NOT local — only guides fa + filelist in `data/external/gse249595/`; SuperSeries GSE247601). **Datlinger 2017 CROP-seq** (doi:10.1038/nmeth.4177; local `data/external/scperturb/DatlingerBock2017.h5ad`: 5,905 Jurkat cells, 32 TCR-panel targets). CRISPRai IL2 work (Nat Biotechnol 2025, doi:10.1038/s41587-024-02213-3, Jurkat + primary T). LINCS L1000: only 256 shRNA sigs, no knockdown bank. | CCLE `JURKAT_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE` in local `ccle_cache.npz` (axis overlaps challenge axis on 17,900/19,209) | 2,393 + genome-wide + 32 off-panel targets for a **calibration map** (341's use); common vector (generic, see §3); basal prior |
| **B = HeLa** | **Nothing broad.** No HeLa Perturb-seq/CROP-seq bank found (literature + local scPerturb catalog of 45 files). Nearest: CROPseq-multi (PMC10979941, 2024) — HeLa CRISPR-KO screens were *viability* (280 gene pairs); its HeLa scRNA-seq is a 5′-capture barcode benchmark, not a perturbation bank. Local cmap: **26,170 HELA siginfo rows are 100% compound/vehicle (+7 OE) — ZERO knockdown (trt_sh) signatures**; `genetic_sigs.parquet` has 0 HELA rows. JURKAT has only 256 trt_sh. | CCLE `HELA_CERVIX` local; DepMap fitness (Chronos) covers HeLa as usual | Basal prior + fitness only. No same-line direction, no calibration targets |
| **C = CAL33** | **Nothing.** CAL33 appears only in viability CRISPR (Brunello genome-wide KO under glutamine blockade, PMC12278720, 2025 — fitness, no transcriptomes) and non-perturbation omics (methylation/metabolites GSE255896, Fanconi-anemia KO clones). No Perturb-seq/CROP-seq in CAL33, FaDu, Detroit 562, Cal27, or SCC-25 found. Local cmap: 0 rows for all five; local DepMap extract (16 lines) covers none of A/B/C. | CCLE has `CAL33_…`, `FADU_…`, `DETROIT562_…`, `CAL27_…`, `SCC25_…` basal rows locally | Basal prior + fitness only |

License/size notes: GEO accessions public (GSE264667, GSE249595 RAW tar 3.75 GB —
over the 2 GB PLAN note threshold, **not downloaded**); scPerturb mirrors per-study
terms; LINCS L1000 level-5 data local under `data/external/cmap/` (trt_sh/trt_xpr
gctx + beta metadata); CCLE/DepMap local caches as above.

## Verdict against the PLAN kill criterion

**Kill criterion fired.** "Dead if Jurkat–HepG2 common cosine ≥ Jurkat–H1-design
cosine": 0.331 ≥ −0.011, by 0.34 with sign agreement 0.66 vs 0.49. The common
response axis is shared across Jurkat/HepG2/K562 (0.30–0.35) — a generic CRISPRi
programme, not Jurkat-specific information — so a **same-line common vector cannot
supply line-specific direction for context A**. The propensity route is independently
dead (§4: Spearman −0.56, noise-dominated, ranks backwards).
**What survives for the PI**: (i) 2,393 (+genome-wide +32) off-panel Jurkat targets
as calibration-map material — the only context where 341's construction is even
buildable; (ii) the generic shared axis is real in held-out-line truth (0.108/target,
85/85 positive) and our whitened design carries none of it (0.007) — restoring *some*
common mode is the only direction-relevant pointer here, but it is **not same-line
information** (any line's bank would do), it faces 341's retrieval guard and 317's
amplitude economics, and nmae reads per-target deviations where a shared additive
constant has no established price. Engineering-feasible; improves-no-member by default.

## Limits

- H1-panel (85 targets, one held-out line) development evidence, not an A/B/C forecast.
- 300-target subsamples (seed 350); full-bank commons not computed. Target lists saved
  in the npz files for exact reuse.
- Design-vs-truth comparisons cross effect spaces (whitened log1p-deltas vs log2
  ratios); the order-independent aggregate (0.004) and the matched-order per-target
  mean (0.007) agree, but a same-space design/truth cosine was not computed.
- Bank control-CPM gating was not re-streamed; extreme-tail (|lfc|>2) low-count
  inflation undiagnosed — mitigated by the top-1% drop test, which leaves the
  Jurkat–HepG2 cosine intact and only attenuates truth–Jurkat.
- GSE249595 paper/supplement (MOESM3 xlsx) content not inspected (no openpyxl locally);
  HeLa/CAL33 web search was 3 queries + catalog greps, not exhaustive.
- No retrieval/cosine-proxy gain is claimed as a score gain (rule 5).

## Exact commands to reproduce

```bash
# known answer (seconds)
uv run python -c "import pandas as pd,json,numpy as np; \
 d=pd.read_csv('experiments/323-current-family-validation/work/per-target-cap15.csv'); \
 n=d[d.metric=='de_wilcoxon_lfc_nmae']; m=json.load(open('experiments/323-current-family-validation/work/metrics-cap15.json')); \
 print(len(n),repr(n.value.mean()),repr(m['raw']['nmae']),n.value.mean()==m['raw']['nmae'])"
# main streaming compute (~6 min, <2 GB): headers + 300-target commons + K562/H1 comparisons + context-A basal + Spearman
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 uv run python experiments/350-same-line-priors-data/work/compute350.py
# follow-up, no X re-read except obs columns (~1 min): depth, H1-truth common, drop tests, per-target ceiling
uv run python experiments/350-same-line-priors-data/work/followup350.py
# load a common vector
uv run python -c "import numpy as np; z=np.load('experiments/350-same-line-priors-data/work/jurkat_common.npz'); print(sorted(z.files), z['vector'].shape)"
```

## Peak memory and wall time

- Wall: recon + probes ~40 min; `compute350.py` ~6 min (15 GB sequential streaming);
  `followup350.py` ~1 min. Total session ≈ 2 h of a 3 h budget.
- Peak RSS (analytic bound, not instrumented): one 8,192×9,624 float64 block ≈ 630 MB
  transient + 301×9,624 accumulators ≈ 23 MB + interpreter; context-A backed reads in
  2,000-row CSR chunks. **< 2 GiB at all times**; `memory_pressure` free stayed ≥70%.
  No process over 4 GiB was started; no scratch over 1 MB retained outside `work/`.

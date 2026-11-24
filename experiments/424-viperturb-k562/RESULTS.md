# 424: VIPerturb-seq K562 is a valid donor table, but it LOSES on PLOLO in every form: (i) −0.029, (ii) −0.022, (iii) the 30-target fill −0.0034 overall (k = 1, 3-fold cancer mean). No builder arm.

Status: **COMPLETE** (2026-09-23, 01:23–02:15Z). Mandate: `MANDATE.md` (PI, 01:35Z).
- Hygiene: no vcc, no cell-eval2, no builder edit, no git, nothing in `submissions/`. Writes went only to this directory
  and the new `data/donors/viperturb_k562/`. Every Python run used `uv run python`.
- No AWS: every bin was converted locally with a streaming RDS reader at ≤ 5.1 GiB RSS. Memory deviations are listed in §6.

## 0. The answer in ten lines

1. **Acquired and converted** all three genome-wide bins plus the manifest (10.2 GB, Zenodo 18460279). All four Zenodo md5s
   match. Local conversion used a new streaming XDR parser (`rds_stream.py`), 63–81 s per bin, peak ≤ 5.1 GiB. The parser
   matches R's own `readRDS` on binC with |d| = 0 over all 19,068 genes. Raw RDS files deleted afterwards.
2. **Built `data/donors/viperturb_k562/`** (17,489 targets × 19,068 genes, 7,949 `NO-TARGET` controls, median 48 cells per
   target) with the production `scripts/build_donor_table.build()` itself, unmodified. Only its two readers were stubbed.
   A rebuild with a different chunking is byte-identical.
3. **Coverage reproduces 413 exactly: 284/300 panel targets at ≥ 20 cells**. 254 are shared with K562-GWPS and 30 fill its
   gaps. Missing: 16 targets (PDK3 is in neither screen).
4. **Knockdown is real:** own-gene delta < 0 for 281/284 panel targets (median −0.52 in log-pseudobulk). K562-GWPS on the
   same targets: median −0.76. Per sgRNA, the panel median is log2 −1.25, and 27% of panel sgRNAs are weaker than −0.5.
5. **Agreement with K562-GWPS is positive but weak.**
   - 159's same-line statistic reads +0.0196 (ratio 7.4; permutation null −0.0006).
   - The known answer K562-essential × K562-GWPS reproduces 159 exactly: +0.2712.
   - On the 254 shared panel targets the centred cosine is +0.018 mean. **HCT116 agrees with K562-GWPS better (+0.026)
     than this second K562 screen does.**
6. **Why:** a VIPerturb panel row is about 8% reproducible signal.
   - Split-half r is +0.043, which gives ρ_V ≈ 0.083 by Spearman-Brown.
   - Only a third of that is shared between a target's own sgRNAs (guide split r +0.0125).
   - VIPerturb's agreement with K562-GWPS (+0.013 per half) sits at exactly that sgRNA-consistent level.
7. **PLOLO (420's instruments, operator and key fixed, every known answer |d| = 0):**
   - **(i) extra donor: −0.0291 / −0.0299** (k = 1 / 1.407, 3-fold cancer mean; −0.044 on the two folds where VIPerturb
     can act).
   - **(ii) averaged K562 row: −0.0222 / −0.0228.**
   - **(iii) 30-target fill: −0.0034 / −0.0042.**
   - Both Orion folds are negative for every arm. Bootstrap P(Δpds > 0) ≤ 0.17 on every Orion fold.
   - The shuffled-row control loses −0.033, barely more than (i): the damage is noise, not misplaced signal.
8. **Robust to the gene-axis reading.** Restricting VIPerturb to K562-GWPS's own gene axis, or to its expression floor, still
   loses in every arm. The best case is (iii) on K562-GWPS's axis at −0.0023.
   - Cell-count (inverse-variance) weights, about 0.2 on VIPerturb, lose −0.0076.
   - A 410-key refit including VIPerturb loses −0.027.
9. **The fill loses for a biological reason.**
   - The 30 panel targets absent from K562-GWPS are genes K562 barely expresses. Replogle's genome-wide screen targeted
     expressed genes.
   - Median 0.16 UMI per cell in VIPerturb's own controls, against 0.91 for the shared targets. EPHB2 has 0 counts;
     KIF21B and SPATA20 have 0.02–0.03.
   - A CRISPRi row for an unexpressed gene is a null profile plus noise. It made the filled targets worse on both Orion
     folds (their own board pds fell 0.741 → 0.705 and 0.737 → 0.684).
10. **Recommendation: no builder arm for (i), (ii) or (iii), for A/B/C or for October.** A K562 row cannot inform a target
    K562 does not express, and for expressed targets this screen is too shallow. The table stays on disk as a validated
    asset, not wired into the builder.

---

## 1. Acquisition record

| item | value |
|---|---|
| source | VIPerturb-seq genome-wide K562 CRISPRi (Bradu … Satija, bioRxiv 2026.02.12.705613), https://zenodo.org/records/18460279, CC BY 4.0 |
| downloaded | `genome_wide_manifest.txt` (399 kB) and `genome_wide_bin{A,B,C}.RDS` (3.56 / 3.78 / 2.89 GB = 10.2 GB), 01:23–01:35Z, about 8 MB/s (`work/dl/download.sh`, `download.log`) |
| integrity | md5 = Zenodo's for all four files. sha256 in `work/dl/sha256.json` and `provenance.json`: binA `eb88c050…`, binB `990da0c4…`, binC `0f24d283…` |
| format | gzip XDR v3 (written by R 4.5.2), Seurat v5, `Assay5` RNA with a single raw `counts` layer (dgCMatrix, all values integer). meta.data columns: orig.ident, nCount_RNA, nFeature_RNA, nCount_GDO, guide, gene, sample, nFeature_GDO |
| size in memory | `readRDS` of binC (the smallest bin) peaks at **16.7 GB**, so binB would be about 22 GB. Above the mandate's ~20 GB line, hence no R path |
| reader | `rds_stream.py`: a one-pass streaming parser of R serialization that spills the counts' `i` / `x` vectors to disk. `convert_bin.py`: per-label count sums in cell chunks |
| conversion | binA 1.68e9 nnz, 77 s, 4.1 GiB. binB 1.78e9 nnz, 81 s, 5.1 GiB. binC 1.36e9 nnz, 63 s, 3.4 GiB. Decompressed stream 16.6–21.7 GB per bin |
| reader audits (every bin) | per-cell colSums == meta `nCount_RNA` and nnz per column == `nFeature_RNA` for **every cell** (max \|d\| = 0). Σ group sums == Σ x |
| reader known answer | R `readRDS` on binC (`ka_r_binC.R`, 17.2 GB, run once under the memory gate). Group sums for NO-TARGET, PSKH1, RBM18, SIN3B, SND1 match with **\|d\| = 0 on all 19,068 genes** |
| AWS | **not used** ($0). `scripts/aws_reap.sh` (dry run, all regions, 02:07Z): one running instance, `peak-gateway` (t4g.small, 439 h, not launched by 424). 424 launched none. `--apply` was not run |
| disk | the RDS files (10.2 GB), per-bin sums (2.9 GB), split-half sums (2.2 GB), axis-variant tables (1.1 GB) and build chunks are **deleted**. Kept: the donor table (1.3 GB), `work/` (284 MB, mostly APFS clones of 420's fold caches and 3 × 42 MB panel-row caches). Rebuilding from scratch means re-downloading (about 20 min) and running `convert_bin.py` (about 4 min) |

## 2. The donor table `data/donors/viperturb_k562/`

- **Construction = K562-GWPS's, by the same code.**
  - `build_vip_donor.py` calls `scripts/build_donor_table.build()` (sha256 in `provenance.json`) and replaces only its
    h5ad reader and `group_sums` with stubs that serve the exact per-label count sums.
  - Everything else is production code: ≥ 20 cells, b = log1p(5e4 · P / ΣP), delta = b_p − b_ctrl, the donor's own gene
    axis, float32, sorted labels, meta.json.
  - To stay inside memory, the build ran on contiguous label chunks (rows are independent, since b_ctrl comes from the
    control alone).
- **Files:** `delta.npy` (17,489 × 19,068 float32), `b_ctrl.npy`, `perts.txt`, `genes.txt`, `meta.json` (same keys as
  K562-GWPS's). The new sidecar `provenance.json` carries URLs, md5 / sha256, the reader and its audits, cells per target
  for all 18,885 labels, each target's bin, the QC applied, and the builder sha.
- **meta:** control `NO-TARGET`, 7,949 cells; 1,396 labels dropped below 20 cells. Cells per target: min 20, median 48,
  max 141. Effect norm: median 11.9 (K562-GWPS 7.04, with a median of 179 cells).
- **Known answers:**
  - KA-bins: NO-TARGET has identical barcodes and identical summed counts in all three bins (used once). Every target is
    in exactly one bin. **Every label's cell count equals the manifest's `cell_count` (0 mismatches over 18,886 labels).**
  - KA-stub: all 17,490 builder groups served with matching cell counts.
  - KA-CHUNK: an independent 8-chunk rebuild (`build_axis_variants.py full`) is **byte-identical** (delta, b_ctrl, perts,
    genes, meta) to the 4-chunk table on disk.
- **QC applied, and why.**
  - No cell filter: the bins hold the authors' post-QC, single-guide-assigned cells.
  - No knockdown filter: K562-GWPS's table keeps every perturbation with ≥ 20 cells.
  - sgRNAs pooled per gene: K562-GWPS pools its dual-guide constructs.
  - Knockdown is measured and reported (§3.2 and `knockdown.csv`), not filtered on.
- **Gene axis.** The builder keeps the source's own axis: 19,068 Flex probes, 18,111 of them on the official axis.
  K562-GWPS's source file already carried only Replogle's 8,248 K562-expressed genes, so the axis rule is ambiguous. Both
  readings are priced (§4.3).

## 3. Validation (c)

### 3.1 Coverage (known answer: 413's 284)

**284 / 300 panel targets at ≥ 20 cells — PASS.** Panel cells: median 48, range 20–122.
- **30 fill targets** (absent from K562-GWPS): ABCD1, ANKRD52, CAPRIN2, DLAT, EPHB2, HEATR6, KIF21B, MAPK7, MLKL, NICN1,
  PARP3, PBLD, PHF11, PHF19, PSMB9, RAB11FIP5, SEMA4F, SLC44A1, SNN, SOCS5, SPATA20, TAF4, TBC1D19, TDRD7, TESK2, TMEM104,
  TMEM79, TTBK2, TXNDC16, ZC2HC1A.
- Shared with K562-GWPS: 254.
- **Missing** (< 20 cells): AKIRIN1, DHX35, EPC2, GPATCH8, HMGN3, LENG1, PAPOLA, PATL1, PDK3, POR, TFCP2, TMEM209, TRIM35,
  UNK, ZHX3, ZNF410. PDK3 is in neither screen.

### 3.2 Knockdown (each target's own gene goes down)

| set | own gene on the axis | down | < −0.1 | median own-gene delta |
|---|---:|---:|---:|---:|
| VIPerturb, panel (284) | 284 | **98.9 %** | 90.5 % | **−0.520** |
| VIPerturb, panel ∩ K562-GWPS (254) | 254 | 99.6 % | 93.7 % | −0.570 |
| K562-GWPS, same 254 (213 on its axis) | 213 | 100 % | 100 % | −0.761 |
| VIPerturb, all 17,489 targets | 16,128 | 74.7 % | 47.8 % | −0.069 (many targets are not expressed in K562) |
| K562-GWPS, all targets | 7,617 | 99.0 % | 96.9 % | −0.657 |

Per sgRNA (log2 own-target CPM in the guide's cells vs NO-TARGET, from per-cell counts):
- panel: 822 sgRNAs, median **−1.25**; 27.1 % weaker than −0.5;
- all 39,756 sgRNAs: median −0.39;
- 19 panel targets whose *best* sgRNA is weaker than −0.5: AGO1, BCL6, FCGRT, FZD2, GLI3, IFI27L2, MBD4, MLKL, MTHFSD,
  MTM1, PBLD, POLD4, PSMB9, SNN, TMCO6, TMEM104, ZNF25, ZNF251, ZNF462.

VIPerturb's knockdown is real but about 30 % shallower than K562-GWPS's on the same targets. Per-target values are in
`knockdown.csv`.

### 3.3 Effect agreement (159's line-pair statistic; known answer first)

`validate424.py`, using 159's `donorpairs.pair()` unchanged: same-target cosine, different-target cosine, their ratio (the
attenuation-invariant column), the centred ("resid") same-target cosine, and its permutation null.

| pair | n targets × genes | same | diff | **ratio** | resid same | resid null |
|---|---|---:|---:|---:|---:|---:|
| **KA: K562-essential × K562-GWPS** (159: 2003 × 6402, +0.2712 / +0.0329 / 8.26, resid +0.2857) | 2003 × 6402 | **+0.2712** | +0.0329 | 8.26 | +0.2857 | — |
| **VIPerturb × K562-GWPS** | 3000 × 5341 | **+0.0196** | +0.0026 | 7.39 | +0.0177 | −0.0006 |
| VIPerturb × K562-essential | 1443 × 6640 | +0.0518 | +0.0077 | 6.74 | +0.0477 | +0.0007 |
| VIPerturb × Orion HCT116 | 3000 × 15569 | +0.0051 | +0.0015 | 3.43 | +0.0036 | −0.0001 |
| VIPerturb × Orion HEK293T | 3000 × 15562 | +0.0046 | +0.0020 | 2.30 | +0.0027 | +0.0001 |
| VIPerturb × kolf2 | 258 × 17843 | −0.0131 | −0.0167 | 0.78 | +0.0041 | +0.0005 |
| (ladder) K562-GWPS × Orion HCT116 | 3000 × 5629 | +0.0322 | +0.0019 | 16.5 | +0.0323 | +0.0012 |
| (ladder) K562-GWPS × Orion HEK293T | 3000 × 5626 | +0.0289 | +0.0020 | 14.6 | +0.0281 | +0.0009 |
| (ladder) K562-essential × Orion HCT116 (159: identical) | 1800 × 6617 | +0.0999 | +0.0123 | 8.11 | +0.1032 | +0.0018 |

**Per-target distributions** (centred per-target cosine; `validate424.json` → `dist`):

| set | n | mean | median | > 0 | retrieval pds among the set |
|---|---:|---:|---:|---:|---:|
| VIPerturb × K562-GWPS, 254 shared panel targets (7,326 genes) | 254 | +0.0184 | +0.0118 | 72 % | 0.665 |
| **K562-GWPS × Orion HCT116, panel (cross-line!)** | 246 | **+0.0263** | +0.0191 | 74 % | **0.706** |
| VIPerturb × Orion HCT116, panel | 275 | +0.0078 | +0.0049 | 62 % | 0.603 |
| K562-essential × K562-GWPS (same line, strong targets) | 2003 | +0.2857 | +0.2418 | 99 % | 0.947 |

In the production E062 effect space on the panel (`diag_rows.py`, 245 targets carried by all four donors), the median
same-target cosines are:

| pair | median cosine |
|---|---:|
| K562-GWPS × VIPerturb | **+0.0097** |
| K562-GWPS × HCT116 | +0.0160 |
| K562-GWPS × HEK293T | +0.0164 |
| HCT116 × HEK293T | +0.0257 |
| HCT116 × VIPerturb | +0.0043 |
| HEK293T × VIPerturb | +0.0036 |

**Reading.**
- VIPerturb is clearly positive against K562. The ratio of 7.4 is near the same-line 8.3 and far above the null.
- Its signal is specific to K562: its ratios against the other lines are 2–3, against 15–17 for K562-GWPS.
- But the signal is attenuated: **per target, this second K562 screen agrees with K562-GWPS less than HCT116 does.**
- Agreement rises with effect size (`diag_depth.py`): mean cosine +0.009 for the weakest half of K562-GWPS effects, +0.063
  for the strongest 5 %.
- It does not rise with VIPerturb cell count (cell count is confounded with fitness).

### 3.4 Why so weak: split-half reliability (`convert_split.py`, `splithalf424.py`)

Random cell halves use disjoint NO-TARGET halves as controls. Guide halves compare the top sgRNA with the target's other
sgRNAs. The space is the table's, on K562-GWPS's gene axis minus the panel, centred.

| targets | n | cells (median) | split-half r | → ρ_V (Spearman-Brown) | guide-split r | VIPerturb half × K562-GWPS | full row × K562-GWPS |
|---|---:|---:|---:|---:|---:|---:|---:|
| panel | 271 | 49 | +0.0433 | **0.083** | **+0.0125** | +0.0133 | +0.0183 |
| strongest 300 K562-GWPS targets | 266 | 42 | +0.0752 | 0.140 | +0.0382 | +0.0472 | +0.0637 |
| random 700 shared | 665 | 49 | +0.0387 | 0.075 | +0.0105 | +0.0081 | +0.0111 |

A VIPerturb panel row is about 92 % noise. Two thirds of its reproducible part is not shared between the target's own
sgRNAs (guide-specific, or a guide without knockdown; 27 % of panel sgRNAs are weaker than log2 −0.5). What it shares with
K562-GWPS (+0.0133) is the same size as its sgRNA-consistent part (+0.0125). The chemistry is not wrong. The screen is
shallow: 48 cells, and knockdown about 30 % weaker than K562-GWPS's.

## 4. PLOLO pricing (d)

### 4.1 Instruments and known answers

- **Harness:** 420's, unchanged.
  - The PP composition; 410's key clip (0.3, 3) **held fixed** at each fold's value (on A/B/C, the pinned K3 artifact
    `d7f76a22`); 417a's sharpening q 1.25.
  - 414's ballast realisation at KMAX 16. Carriers are the arm's own unkeyed PP rows, as a builder that changes the base
    would make them.
  - The a23 read-back with board-emulated pds (λ 0.5).
  - 420's `Pricer` (route r417a: realised A/B/C render energy anchored at 417a's measured energy, the cliff-aware mse
    model, overall = (2.2138 Δpds_bd − 1.051 Δraw mse)/6), called unchanged.
- **Operator held fixed:** the held-line operator on each cancer fold, union_pruned `e2f65ecf` on kolf2 and A/B/C.
- **VIPerturb rows** are prepared exactly as every production donor: its own genome-wide column mean subtracted, scaled to
  EFFECT_SCALE / median centred norm, official axis.
- **Leakage:** VIPerturb is K562, so it is **excluded on the K562-GWPS fold**. Every arm equals the base there by
  construction. mean3 is 420's 3-fold cancer mean (with that zero); mean2 averages HCT116 and HEK293T only.
- **Known answers (all |d| = 0, before any arm number):**
  - KA-420: my base `PP|K3+sh1.25` equals comp420 (design and realised K16) on all four folds, and the fold K3 weights
    equal 418's W3.
  - KA-VIEW: the arm machinery with VIPerturb switched off equals the base `pre` bit-for-bit.
  - KA-ABC420: base A/B/C energies equal abc420.
  - KA-PP: constants vs the 410b meta match to 8.9e-6, 415's value.
  - KA-KEY: the 410 key recomputed on the three donors equals the artifact.
  - KA-P420: 420's pricer reproduces analyse420's `PP|K3+sh1.5` row.
  - KA-BASE: the base priced as an arm reads 0.
- **Controls:** the shuffled VIPerturb rows (`CTRL`) must lose, and do.

### 4.2 The table (primary: `data/donors/viperturb_k562`, full axis)

Δ vs 417a's bulk. Per-fold values are HCT116 / HEK293T, with kolf2 reported without a vote. The K562-GWPS fold is 0 by
construction.

| arm (targets touched per fold) | Δ rel agg cos | Δ pds_bd | Δ pds literal | Δ raw mse (k 1) | **Δoverall k = 1**: mean3 (mean2) | **k = 1.407**: mean3 (mean2) | kolf2 k 1 | A/B/C energy vs 417a | boot P(rel>0)/P(pds>0) hct · hek |
|---|---|---|---|---|---:|---:|---:|---|---|
| **(i)** extra donor, equal weight (275/279 rows) | −0.038 / −0.102 · kolf2 +0.098 | −0.037 / −0.086 · +0.016 | −0.033 / −0.081 | +0.111 / +0.128 | **−0.0291** (−0.0436) | **−0.0299** (−0.0448) | −0.0070 | 0.50 / 0.46 / 0.47 LO | 0.02/0.00 · 0.00/0.00 |
| **(ii)** one K562 row, 0.5/0.5 average + fill | −0.024 / −0.071 · +0.067 | −0.026 / −0.074 · +0.009 | −0.020 / −0.069 | +0.077 / +0.091 | **−0.0222** (−0.0333) | **−0.0228** (−0.0342) | −0.0053 | 0.63 / 0.60 / 0.61 LO | 0.04/0.00 · 0.02/0.00 |
| (ii−) average only, no fill | −0.021 / −0.050 · +0.049 | −0.024 / −0.071 · +0.007 | −0.018 / −0.063 | +0.068 / +0.077 | −0.0201 (−0.0302) | −0.0206 (−0.0309) | −0.0054 | 0.66 / 0.65 / 0.65 LO | 0.06/0.00 · 0.05/0.00 |
| (ii′) average weighted by cell counts (VIPerturb ≈ 0.2) + fill | +0.013 / +0.031 · +0.041 | −0.004 / −0.033 · +0.008 | −0.002 / −0.032 | +0.029 / +0.024 | −0.0076 (−0.0114) | −0.0073 (−0.0109) | −0.0003 | 0.82 / 0.80 / 0.80 LO | 0.93/0.17 · 0.89/0.00 |
| **(iii)** fill only (29/30 targets) | −0.026 / −0.074 · +0.002 | −0.003 / −0.004 · +0.002 | −0.004 / −0.009 | +0.013 / +0.030 | **−0.0034** (−0.0051) | **−0.0042** (−0.0064) | +0.0002 | 1.03 / 1.02 / 1.02 HI | 0.00/0.17 · 0.00/0.15 |
| (iii0) fill, K562 normrel median held (only the fill rows change) | −0.026 / −0.073 · +0.002 | −0.003 / −0.004 · +0.002 | −0.004 / −0.009 | +0.013 / +0.030 | −0.0034 (−0.0051) | −0.0042 (−0.0063) | +0.0002 | 1.03 / 1.02 / 1.02 HI | 0.00/0.17 · 0.00/0.15 |
| CONTROL: (i) with rows permuted | −0.090 / −0.122 · +0.043 | −0.047 / −0.090 · +0.003 | −0.041 / −0.085 | +0.134 / +0.143 | −0.0329 (−0.0494) | −0.0341 (−0.0511) | −0.0160 | 0.44 / 0.42 / 0.42 LO | 0.00/0.00 · 0.00/0.00 |
| (i) + 410 key refit on donors + VIPerturb | −0.070 / −0.119 · +0.081 | −0.034 / −0.082 · +0.018 | −0.033 / −0.082 | +0.100 / +0.114 | −0.0268 (−0.0402) | −0.0279 (−0.0419) | −0.0036 | 0.59 / 0.55 / 0.55 LO | 0.03/0.00 · 0.01/0.00 |
| base + that refit key only | −0.067 / −0.090 · −0.052 | +0.004 / +0.011 · +0.000 | +0.001 / −0.002 | +0.023 / +0.031 | −0.0013 (−0.0020) | −0.0026 (−0.0040) | −0.0028 | 1.01 / 1.01 / 1.02 | 0.05/0.90 · 0.00/0.99 |

Reading:
- **Every arm is negative on both folds where VIPerturb can act, at both transfer readings.** The pds loss alone decides it:
  (i)'s Δpds_bd −0.037 / −0.086 is −0.014 / −0.032 overall before any mse term.
  - The LO flags mean the mse number is an extrapolation: energies of 0.46–0.82 fall outside the 0.88–1.03 range the
    route was validated on.
  - The mse number only adds to the loss.
- **The shuffled control loses almost as much as the real rows** (−0.033 vs −0.029). VIPerturb's target-specific content
  buys back only +0.004 over pure noise. That is the §3.4 reliability made visible.
- **The render energy collapses.** Adding VIPerturb halves the A/B/C render energy of (i). The mean of four rows, one of
  which is mostly noise, carries 23 % less energy before whitening. Whitening then pushes VIPerturb's noise into directions
  the context gate drops.
- **(ii′) raises aggregate cosine** (+1.3 % / +3.1 %) but still costs pds (−0.004 / −0.033). The optimal VIPerturb weight
  is below the cell-count weight of about 0.2, i.e. near 0.
- **kolf2 rises for (i)/(ii) and for the shuffled control too.** That is 159's kolf2 blob (within-line cosine 0.24), not
  transfer, and why kolf2 has no vote.

### 4.3 The gene-axis reading does not rescue it

"Same gene axis rules" is priced both ways. Both variant tables were built with the same production code
(`build_axis_variants.py`) and priced end to end (`run_variant.sh`). All known answers are |d| = 0 again.

| arm, Δoverall k = 1 mean3 (mean2) | full probe axis (on disk) | K562-GWPS's gene axis (7,539) | K562-GWPS expression floor (9,209) |
|---|---:|---:|---:|
| (i) extra donor | −0.0291 (−0.0436) | −0.0129 (−0.0193) | −0.0261 (−0.0391) |
| (ii) averaged K562 row + fill | −0.0222 (−0.0333) | −0.0096 (−0.0144) | −0.0189 (−0.0284) |
| (ii−) averaged, no fill | −0.0201 (−0.0302) | −0.0084 (−0.0125) | −0.0169 (−0.0254) |
| (ii′) cell-count weights + fill | −0.0076 (−0.0114) | −0.0043 (−0.0064) | −0.0068 (−0.0101) |
| (iii) fill only | −0.0034 (−0.0051) | −0.0023 (−0.0034) | −0.0038 (−0.0058) |
| CONTROL shuffled | −0.0329 (−0.0494) | −0.0151 (−0.0226) | −0.0291 (−0.0437) |
| kolf2, (i) | −0.0070 | −0.0014 | −0.0065 |

K562-GWPS's own axis halves the damage, because it drops genes where VIPerturb is pure noise. Every arm stays negative on
both Orion folds, and on kolf2 too. 90 % of a VIPerturb row's energy already sits on genes above K562-GWPS's expression
floor, so the noise is on expressed genes: it is depth, not axis.

### 4.4 The 30-target fill, per target

Target-gene expression in VIPerturb's controls:
- the 30 fill targets: median **0.16 UMI/cell**; 12 of 30 below 0.1;
- the 254 shared targets: median 0.91;
- EPHB2 0.000, PSMB9 0.004, MLKL 0.008, SNN 0.015, KIF21B 0.020, SPATA20 0.032.

K562-GWPS lacks these targets because Replogle's genome-wide screen targeted K562-expressed genes. Of the 413 "thin five",
EPHB2, KIF21B and SPATA20 are essentially unexpressed in K562 (ABCD1 0.25, NICN1 0.20). Six of the 30 show no knockdown at
all: best sgRNA weaker than log2 −0.5, or no counts.

On PLOLO the fill **hurt the targets it fills** (`arrays424-*`):

| fold | fill targets | per-target cos, base → (iii) | improved | own board pds, base → (iii) | other targets' pds |
|---|---:|---|---:|---|---|
| HCT116 | 29 | +0.029 → +0.021 | 12/29 | 0.741 → 0.705 | unchanged (0.822 → 0.822) |
| HEK293T | 30 | +0.027 → +0.015 | 6/30 | 0.737 → 0.684 | unchanged (0.808 → 0.809) |
| kolf2 | 23 | +0.007 → +0.008 | 14/23 | 0.515 → 0.543 | unchanged |

### 4.5 What a refit would change

- **Key (measured).** Refitting 410's shared-energy key with VIPerturb as a fourth donor carries 30 more fold targets
  (259 → 289 on HCT116), but it reorders the weights (Spearman 0.50–0.54 on the Orion folds, 0.69 on kolf2). It costs −0.0013 on the base and does not rescue
  (i) (−0.0268 vs −0.0291).
- **CD4 / H1 / partial fractions** are re-solved inside every arm; it is the builder's rule. f_cd4 falls from 0.346 to
  0.298 in (i) on A/B/C, because the base median norm falls.
- **Operator (not measured; the mandate held it fixed).**
  - Refitting the whitening covariance with VIPerturb rows would add a mostly-noise K562 table (ρ_V ≈ 0.08) to C.
  - That damps the directions where VIPerturb is noisy, so it could recover part of the render-energy collapse.
  - It would also change every row's geometry and move the pinned operator (`e2f65ecf`), which the October runbook requires
    to reproduce.
  - It cannot create target-level signal VIPerturb does not carry: the shuffled control shows the loss is noise, and the
    sgRNA-consistent split-half r is 0.0125. Expected sign: still ≤ 0. Not worth a cloud refit.

## 5. Recommendation

**None of (i), (ii) or (iii) is worth a builder arm.** No builder change is proposed.
- (i) and (ii) inject a K562 row that is about 92 % noise into a three-line mean and lose 0.02–0.03 overall.
- (iii) fills 30 targets with rows for genes K562 does not express. It loses 0.003–0.005 and hurts exactly the filled
  targets.
- The loss holds at k = 1 and k = 1.407, on both Orion folds, on all three gene-axis readings, with cell-count weights and
  with a key refit.

If someone wants it anyway, the builder change would be:
- register `viperturb_k562` as a donor bank in `make_submission_pm_cd4.py`'s donor list (E062 preparation, pinned shas);
- then either add it to the normrel mean, for (i); or fold it into the K562 slot at 0.5 / 0.5 with PM50 scaled by 0.5, for
  (ii); or use it as a K562-slot fill behind a new knob such as `VCC_K562_FILL=viperturb`, for (iii).

**October D/E/F.**
- The fill is the only use that grows on a disjoint panel. The mechanism says it will not help.
- Targets outside K562-GWPS are mostly genes K562 does not express, and a K562 row cannot inform them. The fill targets on
  these folds had exactly one donor line, which is October's situation, and it still lost.
- Keep the table as a validated asset, not a donor.

**What VIPerturb might still be good for** (untested, low prior): a many-target *training* set for K562 response programmes,
where 17,489 noisy rows average out, and a per-sgRNA knockdown reference (`knockdown.csv`, `work/conv/*/cells.tsv.gz`). It
is not a per-target donor.

**For the record (FINDINGS §2 candidate):** "a second, shallower screen of a line we already hold (VIPerturb K562, 48
cells/target, ρ ≈ 0.08 per row) — as an extra donor, an averaged K562 row, or a gap fill — loses 0.003–0.03 on PLOLO. Gap
targets missing from a line's genome-wide screen are usually genes that line does not express."

## 6. Memory and rule deviations (said plainly)

The memory gate was checked before every heavy step (74–78 % free, free+inactive 32 GiB).
- **`build_vip_donor.py` peaked at 11.6 GB RSS for 7 s**, over the 8 GB local-Python cap.
  - Cause: 4-chunk float64 sums, their temporaries, and mmap'd per-bin sums.
  - The table was rebuilt with an 8-chunk, no-mmap builder at 3.95 GiB and is byte-identical (KA-CHUNK).
- **`validate424.py` peaked at 12.3 GB RSS.**
  - Cause: mostly file-backed mmap pages of the Orion donor tables (2.7 GB each) and VIPerturb's table, touched by 159's
    pair statistic.
  - Both runs were short and happened with 75 % free; nothing else was affected.
- **R known answer `ka_r_binC.R`: 17.2 GB** (plus a first attempt that stopped on a missing SeuratObject after reading, at
  16.7 GB).
  - This was R, not Python, within the mandate's per-bin ~20 GB line, run once under the gate.
- All pricing runs: ≤ 6.5 GiB. All bin conversions: ≤ 5.1 GiB.
- `scripts/aws_reap.sh` was run in its default dry-run mode only (read-only). The mandate asked for it only if AWS was used.

## 7. Files

| file | what |
|---|---|
| `rds_stream.py` | streaming reader for R XDR serialization (spills large vectors; ALTREP, S4 and pairlists handled; anything unknown raises) |
| `convert_bin.py`, `run_convert.sh`, `convert_bin.R`, `ka_r_binC.R` | bin → per-label sums (Python, streaming); the R route (abandoned on memory) and R's known answer |
| `build_vip_donor.py` | builds `data/donors/viperturb_k562` with the production `build_donor_table.build()` + provenance sidecar |
| `build_axis_variants.py`, `run_variant.sh` | chunk-invariance known answer; gene-axis variants and their pricing |
| `validate424.py`, `knockdown.csv`, `work/validate424.{json,log}` | §3.1–3.3 |
| `diag_rows.py`, `diag_depth.py` | E062-space row agreement and norms; agreement by effect size and cells |
| `convert_split.py`, `splithalf424.py`, `work/splithalf424.json`, `work/split/targets.json` | §3.4 |
| `lib424.py`, `comp424.py`, `abc424.py`, `analyse424.py` | §4: arms, PLOLO folds, A/B/C energies, pricing through 420's Pricer |
| `work/comp424*-<fold>.{json,log}`, `work/arrays424*-<fold>.npz`, `work/abc424*.{json,log}`, `work/analyse424*.{json,log,stdout}` | per-fold results (primary and `_gwpsaxis` / `_floor` variants) |
| `work/dl/` | manifest, `sha256.json`, download log (RDS deleted) |
| `work/conv/bin*/` | per-bin features, labels, cell counts, `cells.tsv.gz` (per-cell guide, library size, own-target count), `info.json` (audits) |
| `data/donors/viperturb_k562/` | `delta.npy` (17,489 × 19,068 f32), `b_ctrl.npy`, `perts.txt`, `genes.txt`, `meta.json`, `provenance.json` |

**Git note for the coordinator.** `data/donors/viperturb_k562/delta.npy` is 1.33 GB and is not gitignored (it shows as
untracked). Never stage it (AGENTS.md §10: nothing over 100 MB). Commit only `meta.json`, `provenance.json`, `perts.txt`,
`genes.txt` and `b_ctrl.npy` from that directory, if anything. Nothing in `experiments/424-viperturb-k562/` exceeds 45 MB.

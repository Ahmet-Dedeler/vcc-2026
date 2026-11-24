# 295 — claim C5 on reserve R3: does the conservation map predict transfer into eight untouched datasets?

**C5a (F8 on untouched data): PASS** (ρ **+0.5016** [+0.1525, +0.5659] file-block, [+0.3709, +0.6083] gene
bootstrap, n = 174; annotation **+0.4211**; magnitude **+0.3578**); **C5b (F13's sign): PASS** (`gen_pc20` ρ
**+0.3709** [+0.0635, +0.4148]; `breadth` **+0.4365** [+0.1225, +0.5207]); **CRISPRi-only +0.5169 [−0.0435,
+0.7381] / Cas9-only +0.2137 [−0.0766, +0.4260]**.

Preregistration `PREREGISTRATION.md` sha256 `a8a49c637c6c24e9503513d84576030c13b911b0dbc8bb89d455751ff78b8a7d`,
hashed and written into `docs/CONFIRMATION-LEDGER.md` as claim C5 at **2026-09-06T16:09:10Z**, before the first
byte of any reserve `X` was read. **One access, 2026-09-06 16:14–16:15Z, 28 s, peak RSS 2.89 GiB; all eight files
are now SPENT.** No build, no `vcc`, no `cell-eval2`, no download, nothing written outside this directory except
the two ledger entries.

---

## 1. Known-answer rows — all before any reserve number

| row | required | measured | verdict |
|---|---|---|---|
| **KA-1** | every R3 file's sha256 equals `reserve292.json` | **8/8 MATCH** | PASS |
| **KA-2** | 279's `c_held`: n 17,106, median +0.018605 | n **17,106**, median **+0.01860516** | PASS |
| **KA-3** | 291's C4 cell rebuilt from `pred291_gwcd4i.csv`: Stim48hr ρ(`C_held`) +0.2739 | **+0.2739** at n **1,773**; `gen_pc20` **+0.1858**, `breadth` **+0.1275** | PASS |
| **KA-4** | the map-consensus construction reproduces `held279.npz` `Ccross` | **6,097** (pair, reserve gene) cosines, max abs diff **0.0e+00** | PASS |
| **KA-5** | computed library size == `obs['ncounts']`, every cell of every file | max abs diff **0**, max relative diff **0**, on all 8 files (427,073 cells) | PASS |
| **KA-6** | a controls-vs-controls split reads no signal against the map consensus | per-file median cosine **+0.0094 / +0.0042 / +0.0048 / +0.0042 / −0.0137 / −0.0065 / +0.0002 / −0.0122** (‖·‖ < 0.014 everywhere) | PASS |

Two internal assertions also had to hold before anything was computed and did: every file's control-cell count and
every gene's cell count reproduced 292's frozen label grammar exactly (8 files, 205 gene groups, 427,073 cells).
The first run of the access **aborted** on that assertion — Adamson's 2,613 guide-unassigned cells (pandas code
−1) were being mislabelled as the last category. That is what the assertion is for; it was fixed to reproduce
292's serialisation (`""`, which its frozen `CTRL_PAT` counts as a control) and the run was repeated from the
start. No statistic had been computed at that point.

## 2. What was measured, per file

Response = `log1p(5e4·S/ΣS)` of the gene's pooled cells minus the file's own control profile, minus the mean over
the file's own scored genes, on 279's frozen 2,739-gene axis restricted to the file's `var`; the on-target column
removed. `A_g` = its cosine with the 28-bank map's consensus residual for that gene over the banks whose LINE is
not the file's (Tian also excludes `tian_ipsc_neuron`).

| file | modality | (file,gene) rows | primary rows | cells/gene (med) | UMI/cell (med) | split-half `r` | DE genes (med) | median `A_g` | share `A_g` ≤ 0 | median `C_held` |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `AdamsonWeissman2016_GSM2406681_10X010` | CRISPRi | 80 | 78 | 586 | 15,294 | +0.572 | 250 | **+0.1136** | 0.21 | +0.1334 |
| `DatlingerBock2017` | Cas9 KO | 28 | 8 | 152 | 6,688 | +0.026 | 1 | +0.0059 | 0.39 | +0.0127 |
| `DatlingerBock2021` | Cas9 KO | 20 | 19 | 1,702 | 460 | +0.117 | 45 | +0.0090 | 0.30 | +0.0165 |
| `DixitRegev2016` | CRISPRi | 8 | 8 | 2,434 | 13,877 | +0.077 | 27 | +0.0452 | 0.12 | +0.0385 |
| `PapalexiSatija2021_eccite_RNA` | Cas9 KO | 22 | 19 | 662 | 11,845 | +0.257 | 1 | +0.0089 | 0.50 | +0.0511 |
| `PapalexiSatija2021_eccite_arrayed_RNA` | Cas9 KO | 7 | 2 | 750 | 6,982 | +0.500 | 1,563 | +0.0024 | 0.43 | +0.0366 |
| `ShifrutMarson2018` (primary CD8⁺ T, the NEW context) | Cas9 KO | 17 | 17 | 568 | 5,358 | +0.072 | 1,416 | **−0.0035** | 0.59 | +0.0112 |
| `TianKampmann2019_day7neuron` (RESTRICTED) | CRISPRi | 23 | 23 | 1,015 | **15** | +0.060 | 13 | +0.0422 | 0.22 | +0.0458 |

Pooled: **median `A_g` +0.0381, 29.3 % ≤ 0**, n = 174 genes (205 (file, gene) rows de-duplicated by the
preregistered "most cells" rule). This is the same regime the programme measured elsewhere: 276/279's cross-line
median is +0.0186, 291's Broad-16 surface +0.0296.

**Depth, reported and never conditioned on (programme rule 9).** Two files are shallow in ways the metadata could
not show: `DatlingerBock2021` is scifi-RNA-seq at **460 UMI/cell**, and `TianKampmann2019_day7neuron`'s labelled
droplets carry a median of **15 UMI/cell** (the file holds all 182,790 barcodes, and guide labels reach deep into
the background pool; `obs['ncounts']` confirms it and KA-5 proves our sums are its sums). Across genes,
ρ(`A_g`, cells/gene) = **−0.076**, ρ(`A_g`, split-half `r`) = **+0.467**, ρ(`A_g`, UMI/cell) = **+0.346**,
ρ(`A_g`, DE count) = **+0.393** — the same ordering 280/285 found: what predicts agreement is the measured
reliability and the size of the response, not the cell count.

## 3. C5a — F8 off the map, on untouched data: **PASS**

`ρ(C_held, A_g)` over the 174 pooled genes = **+0.5016**; file-block bootstrap (2,000 draws, the 8 files as
blocks) **[+0.1525, +0.5659]**; gene bootstrap **[+0.3709, +0.6083]**. Rule: ρ ≥ +0.15 **and** the file-block
interval excludes 0 — both hold. Null: 200 gene shuffles read **+0.002 ± 0.072**, **0/200** reach the observed
value, **z +6.9**.

| baseline (same rows, reported, not in the rule) | ρ | file-block 95 % | paired Δ (`C_held` − it) |
|---|---:|---|---|
| `M_ann` — 284's out-of-fold annotation ridge | +0.4211 | [+0.1275, +0.5058] | +0.0797 [−0.0747, +0.1767], P(Δ≤0) 0.099 |
| `mag` — 279's response magnitude | +0.3578 | [+0.0087, +0.4452] | **+0.1483 [+0.0382, +0.2916], P(Δ≤0) 0.011** |
| the file's own split-half reliability | +0.4672 | [−0.0389, +0.6015] | +0.0721 [−0.0456, +0.2862] |
| the gene's DE count in the file | +0.3934 | [−0.0453, +0.6215] | +0.0671 [−0.2466, +0.3909] |
| `r_bar` — the map's own reliability | +0.4336 | [+0.0699, +0.5154] | +0.0611 [−0.0893, +0.2091] |
| `C_held_excl` — `C_held` recomputed with EVERY bank of the file's line dropped from both sides of each pair | +0.4593 | [+0.1663, +0.5125] | +0.0331 [−0.0030, +0.0612] |
| cells per gene | −0.0760 | [−0.2395, +0.2534] | — |

**Partial Spearman given the file's reliability and the map's magnitude: +0.3086 [+0.1389, +0.3964]**; adding log
cells and log UMI/cell: **+0.3189 [+0.0638, +0.3989]**. The association is not the reliability confound and not
the depth confound.

`C_held` is the best single predictor of transfer on this reserve, but only its margin over **magnitude** is
separated at 8 blocks (Δ +0.1483, P(Δ≤0) 0.011); the margin over the annotation model is +0.0797 with an interval
that includes 0 — consistent with 291's reading that on out-of-map surfaces the map's score and annotation are
close, and with 288's finding that magnitude is the live confound.

## 4. C5b — F13's SIGN off the map, on untouched data: **PASS**

`ρ(gen_pc20, A_g)` = **+0.3709**, file-block **[+0.0635, +0.4148]** — entirely above 0, so the rule passes.
`breadth` (reported, not in the rule) **+0.4365 [+0.1225, +0.5207]**. Both are positive: **broad, spectrally
central responses transfer better**, the inverted sign 290 found and 291 confirmed on two other surfaces
(+0.2952 Broad-16, +0.1858 GWCD4i). This is the third independent surface and the first that is untouched data.
`gen_pc20` does **not** beat the annotation model here either (paired Δ −0.0407 [−0.2244, +0.0942]) — as in 291,
F13's *headline* stays withdrawn off-map; only the sign travels.

## 5. Declared strata

| stratum | n | C5a ρ [file-block] | C5b ρ [file-block] |
|---|---:|---|---|
| **pooled (verdict, Tian included)** | 174 | **+0.5016 [+0.1525, +0.5659]** | **+0.3709 [+0.0635, +0.4148]** |
| pooled without Tian day-7 | 151 | +0.5396 [+0.1628, +0.5749] | +0.3804 [+0.0017, +0.4235] |
| CRISPRi (Adamson, Dixit, Tian) | 109 | +0.5169 [−0.0435, +0.7381] | +0.3415 [−0.2619, +0.3567] |
| Cas9 KO (Datlinger ×2, Papalexi ×2, Shifrut) | 65 | +0.2137 [−0.0766, +0.4260] | +0.1659 [−0.0437, +0.2971] |

Both strata point the same way and neither is separated from 0 on its own — 3 and 5 blocks give intervals of
±0.3 and ±0.25. The point estimates say the map, which is dominantly CRISPRi, predicts CRISPRi transfer about
2.4× better than knockout transfer; **that ordering is a point estimate on 8 blocks, not a measured difference.**
Dropping the same-study Tian file does not change either verdict.

## 6. Labelled secondaries (all preregistered, none a verdict)

| variant | C5a ρ [file-block] | C5b ρ [file-block] |
|---|---|---|
| all 205 (file, gene) rows | +0.4928 [+0.1968, +0.5697] | +0.3438 [+0.0617, +0.3953] |
| on-target column KEPT | +0.5480 [+0.1476, +0.6066] | +0.3754 [+0.0671, +0.4248] |
| disattenuated `A_g / √r` (164 rows with r > 0) | +0.4685 [+0.1765, +0.5385] | +0.3520 [+0.1190, +0.4018] |
| **within-file ranks** (kills every between-file offset) | **+0.3508 [+0.0773, +0.4852]** | **+0.2310 [+0.0441, +0.3069]** |
| `C_held_excl` as the predictor | +0.4593 [+0.1663, +0.5125] | — |

The within-file row matters: both verdicts survive when the only information used is the *ordering of genes inside
one file*, so neither is an artefact of files with better data also having genes the map happens to score high.

## 7. Descriptive, post-hoc, labelled

**Best-transferring** (all Adamson/K562 unless noted): `HSPA5` +0.689 (r 0.913), `XRN1` +0.639, `ASCC3` +0.590,
`TTI1` +0.527, `DHDDS` +0.490, `CHERP` +0.476, **`STAT1` +0.456 in Papalexi/THP-1** (r 0.855), `PSMD4` +0.441,
`HSPA9` +0.410, `SEC61G` +0.400. The head is proteostasis / ER-translocon / spliceosome / proteasome machinery —
279's "core machinery" end of the conservation axis, arriving in a study the map never saw.

**Worst:** `HYOU1` −0.251, `CBLB` −0.129 (Shifrut), `PTDSS1` −0.126, `UQCRQ` −0.118 (Tian), `CCND3` −0.113,
`CD3D` −0.090 (Shifrut), `SLC35B1` −0.087, `PDCD1LG2` −0.083 (Papalexi), `TARS` −0.083, `TIMM23` −0.077. The
T-cell-identity genes in the one NEW context sit at the bottom, which is what the axis predicts.

**The two same-study pairs R292 §2.3 named.**
*Papalexi pooled → arrayed* (7 genes, same lab, same paper, one gene per well vs a pooled screen with guide
calling): ρ(`A_pooled`, `A_arrayed`) **+0.643**, medians −0.020 vs +0.002, split-half reliability 0.262 vs
**0.500**. The arrayed measurement is twice as internally reproducible and calls 1,563 DE genes against 1 — a
per-well batch signature that reproduces within a gene without transferring to the map. **Removing the guide-calling
step raises reliability and DE yield sharply without raising agreement with the map: reproducibility and
transferability are different quantities**, the same dissociation 285's construct term and 267's oracle produced.
*Datlinger 2017 → 2021* (20 genes, same line, same TCR library, five years and a chemistry apart):
ρ(`A_2017`, `A_2021`) **+0.212**, medians +0.008 vs +0.009, reliability 0.029 vs 0.117. At 152 and 1,702 cells per
gene and 6,688 / 460 UMI/cell, both sides are near the noise floor; the pair is consistent with agreement but
cannot establish it.
*The three cross-study genes*: `IRF1` +0.049 (Dixit/K562, r 0.034), **+0.177 (Papalexi pooled, r 0.725), +0.167
(Papalexi arrayed, r 0.666)**; `EGR1` +0.042 / −0.002 / +0.019 (Dixit, Datlinger 2021, Datlinger 2017, all
r ≤ 0.12); `SOCS1` +0.033 (Shifrut) / +0.021 (Adamson).

## 8. Reading

The 28-bank conservation map is **not an artefact of the sources it was built from.** On eight datasets that no
bank, map row, donor table or fold has ever touched — five studies the programme had never opened, four of them a
different perturbation modality (Cas9 knockout, not the map's CRISPRi), one an entirely new cell type (primary
human CD8⁺ T cells) — a gene's map score predicts how well its response actually transfers at ρ +0.50, with the
association surviving reliability, magnitude, depth, and within-file ranking, and 0 of 200 shuffles reaching it.
F13's inverted sign travels to a third surface. Both were frozen before a single expression value was read.

Two honest limits. (i) The pooled interval is wide because 8 files is 8 blocks and one file (Adamson, 78 of the
174 primary genes) carries the pool; the file-block lower bound is +0.15, only just clear of the +0.15 level the
rule asks of the point estimate. (ii) The map's advantage over the *annotation* model is +0.08 with an interval
across 0 — this reserve confirms F8's substance, not any ordering of predictors, exactly as 291 found.

**Nothing here is competition-facing and nothing is proposed for the model.** R3 holds 6 of the 300 panel targets
and none of its lines is A/B/C; no threshold, feature or dose was selected on it. **All eight files are SPENT.**

## 9. Files

`PREREGISTRATION.md` (`a8a49c63…`, digest in `work/PREREGISTRATION.sha256`) ·
`work/meta295.py`, `work/meta295.json`, `work/reserve_genes295.csv` (metadata, pre-hash) ·
`work/buildmap295.py`, `work/map295.npz`, `work/ka_map295.json`, `work/buildmap295.log` (the map side, KA-4) ·
`work/access295.py`, `work/access295.csv`, `work/access295.json`, `work/access295.log` (**the single access**) ·
`work/stats295.py`, `work/stats295.json`, `work/stats295.log`, `work/paired295.py`, `work/paired295.json` ·
**`work/reserve295.csv`** — the per-gene table (file, gene, cells, guides, UMI/cell, `r_split`, DE count, `A`,
`A_onT`, `n_banks`, `C_held`, `C_held_excl`, `mag`, `logmag`, `gen_pc20`, `breadth`, `M_ann`, `r_bar`, modality,
line, `primary`).

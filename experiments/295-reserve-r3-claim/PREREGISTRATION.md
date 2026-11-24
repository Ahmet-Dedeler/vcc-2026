# PREREGISTRATION — claim C5 on reserve R3

**Experiment 295 (`experiments/295-reserve-r3-claim/`). Written and hashed BEFORE the first byte of `X` of any R3
file is read.** R3's access rule (`experiments/292-reserve-declaration/RESERVE.md` §3, mirrored in
`docs/CONFIRMATION-LEDGER.md`) names this experiment as the one that spends the reserve. **One access, all eight
files, one pass each; the result is appended to the ledger whatever it is; the eight files are then SPENT.**

The claim tests, on data no map, bank, donor table or fold has ever touched, the two live cross-line-conservation
findings of `docs/FINDINGS-2026-09.md` v1.9:

- **F8** — the 28-bank map's per-gene conservation score `C_held` is a reproducible property that predicts how well
  a gene's perturbation response transfers to a line the map did not measure it in (confirmed off-map by 291 on two
  previously-inspected surfaces: +0.274 on the GWCD4i authors' CD4-vs-K562 statistic, +0.210 on Broad-16).
- **F13's SIGN** — broad, spectrally central responses (`gen_pc20`) are MORE conserved (291: +0.186 and +0.295 on
  the same two surfaces). F13's *headline* ("geometry beats annotation") was withdrawn off-map by 291 and is **not**
  re-tested here.

---

## 0. What has been read before this file was hashed

**Metadata only, and nothing else.** (a) The R3 manifest `experiments/292-reserve-declaration/work/reserve292.json`
and the sha256 of each of the eight files (a digest of the bytes; explicitly free under the access rule) — all eight
**MATCH** the manifest. (b) 278's obs-only label counts and 292's frozen label grammars (`work/parse292.py`), giving
each file's perturbation labels, cells per gene, guides per gene and control labels. (c) Each file's `var` index
(gene symbols) and its overlap with the map's axis. (d) The HDF5 **structure** of each file: group and dataset
NAMES, the `encoding-type` / `encoding-version` / `shape` attributes of `X`, the length and dtype of `X/data`,
`X/indices`, `X/indptr`, and the (empty) `layers` / `obsm` group listings. No dataset holding an expression value
was sliced, indexed or summed. Structural attributes are declared here to be in the same class as `obs` / `var` /
`uns`: they are not derived from any expression value. Recorded in `work/meta295.json`.

What that metadata fixes, before any outcome exists:

| file | line | modality | genes ≥ 50 cells **and** in the held map | axis coverage | control cells | banks excluded from the consensus |
|---|---|---|---:|---:|---:|---|
| `AdamsonWeissman2016_GSM2406681_10X010` | K562 | CRISPRi | 80 | 2,732 / 2,739 | 10,007 | `gasperini_k562`, `jiang_K562`, `k562_essential`, `k562_gwps`, `perturbmatch` |
| `DatlingerBock2017` | Jurkat | Cas9 KO | 28 | 2,732 | 1,320 | `nadig_jurkat` |
| `DatlingerBock2021` | Jurkat | Cas9 KO | 20 | 2,636 | 4,497 | `nadig_jurkat` |
| `DixitRegev2016` | K562 | CRISPRi | 8 | 2,731 | 10,247 | the five K562 banks above |
| `PapalexiSatija2021_eccite_RNA` | THP-1 | Cas9 KO | 22 | 2,725 | 2,386 | `wessels_thp1` (carries 0 reserve genes anyway) |
| `PapalexiSatija2021_eccite_arrayed_RNA` | THP-1 | Cas9 KO | 7 | 2,723 | 2,009 | `wessels_thp1` |
| `ShifrutMarson2018` | primary CD8⁺ T (NEW) | Cas9 RNP KO | 17 | 2,739 | 39,813 | none — no bank is in this line |
| `TianKampmann2019_day7neuron` (RESTRICTED) | iPSC-neuron d7 | CRISPRi | 23 | 2,739 | 15,580 | `tian2019_ipsc` **and** `tian_ipsc_neuron` (same study / same iPSC line) |

**205 (file, gene) rows over 174 distinct genes**; 29 genes appear in more than one file (20 of them the
Datlinger 2017/2021 pair, 7 the Papalexi pooled/arrayed pair, 3 across studies: `EGR1`, `IRF1`, `SOCS1`).

The MAP side is also already built (`work/map295.npz`, `work/buildmap295.log`) — it reads only the 28 map banks,
never the reserve — and its construction check has already passed: see KA-4 below.

---

## 1. The outcome, per reserve row (computed exactly once)

Let `f` be a reserve file and `g` one of its genes (≥ 50 cells, in 279's 17,106 held genes; the frozen 292 grammar
resolves labels to genes and pools that gene's guides).

**Response space.** 279's frozen 2,739-gene axis (`held279.npz` `axis`), restricted to the genes `f`'s `var`
carries (2,636–2,739 per the table above).

**Reserve response.** With `S_G` the summed raw counts of a cell group `G` over the axis columns and `T_G` the
summed library size of those cells over ALL genes,
`prof(S_G) = log1p(5e4 · S_G / T_G)` (`lib285.prof`), `Δ_G = prof(S_G) − prof(S_ctrl)` with `ctrl` = the file's own
control cells (`lib285.delta`), and the residual `R_g = Δ_g − mean_{g' ∈ f} Δ_{g'}` over exactly the file's own
scored genes (`lib285.resid`, DEFINITION-276 §2 applied to the file). Note the control term cancels exactly in
`R_g`; it is still computed and reported, and it is used by the DE count and the no-signal control below.

**Map consensus.** For the same `g`, the mean over the map banks that carry it and whose LINE is not `f`'s line
(plus the two extra exclusions declared for Tian) of that bank's UNIT residual row, built by 279's own
`build279.load_unit` on the full 17,106-gene held list so the centring is the map's (`work/map295.npz`).

**`A_g` (the outcome).** The cosine between `R_g` and the map consensus for `g`, over the axis columns `f` carries,
**with the column of `g` itself removed from both sides** (the on-target column, which both sides drive by
construction). Labelled secondary, same access: `A_g^onT`, the identical cosine with that column kept.

**Split-half reliability `r_g`.** 8 draws (seed 20260906); each draw splits `g`'s cells at random into halves A and
B, forms each half's residual with that half's own across-gene mean, and takes the cosine of the two halves'
residual rows; `r_g` = the mean over the 8 draws. Reported next to every cell and used to disattenuate
(`A_g / sqrt(max(r_g, 0))`, reported only, never in a pass rule).

**Depth (programme rule 9), measured in the same access.** cells per gene, the median library size (UMIs/cell) of
`g`'s cells, and the file's own median UMIs/cell. **Reported next to every cell; the pass rules are unconditioned,
as declared.**

**DE count `d_g` (a magnitude proxy, not the scorer's statistic).** Welch t-test of `log1p(CPM)` per axis gene,
`g`'s cells vs the file's control cells, over axis genes with pseudobulk control CPM > 5, Benjamini–Hochberg at
q < 0.05; `d_g` = the number of genes called.

**The analysis unit.** **Primary: the 174 distinct genes.** A gene measured in more than one reserve file
contributes the row of the file with the MOST cells for it (ties → alphabetically first file). This rule is fixed
now, involves no outcome, and matches R292 §6's "one pooled claim over ~174 genes". Labelled secondary: all 205
(file, gene) rows.

---

## 2. The predictors and baselines (all fixed before the access)

| name | definition |
|---|---|
| `C_held` | 279's per-gene cross-line conservation, `held279.npz` `c_held` (mean over the 364 cross-line bank pairs) |
| `C_held_excl` | *labelled secondary*: the same mean restricted to pairs in which NEITHER bank is in the reserve file's line |
| `gen_pc20` | 290's per-source `pc20`, nanmean over the file's non-excluded banks (`lib291.Store.aggregate` arithmetic) |
| `breadth` | 290's per-source `breadth`, same aggregation (reported) |
| `logmag` | `log(nanmean(mag))` over the same banks (290) |
| `mag` | 279's per-gene response magnitude, `held279.npz` `mag` |
| `M_ann` | 284's out-of-fold annotation-ridge prediction of `C_held`, `pred284.csv` `pred_oof_raw` (21 annotation features, all-membership family folds) |
| `r_bar` | the map's own per-gene reliability, `pred284.csv` `r_bar` |
| `r_g`, `d_g`, cells, UMIs | measured in this access (§1) |

---

## 3. The two verdict cells

**Statistics.** Spearman ρ. Two bootstraps, 2,000 draws each, seed 20260906: (i) the **file-block** bootstrap —
resample the reserve FILES with replacement and take all rows of the drawn files (`lib291.boot_rho` with the file
as the group); (ii) a plain **gene** bootstrap over rows. Intervals are the 2.5/97.5 percentiles.

### C5a — F8 off-map, on untouched data
`ρ(C_held, A_g)` over the pooled primary rows.
**PASS iff ρ ≥ +0.15 AND the file-block 95 % interval excludes 0.** Otherwise FAIL.
Baselines, same rows, reported not in the rule: `M_ann`, `mag`, `r_g`, `d_g`, `r_bar`, `C_held_excl`; and the
partial Spearman `ρ(C_held, A_g | r_g, mag)` (`lib291.partial_spearman`).
**Null:** 200 shuffles of `A_g` against the gene labels; report the null mean, sd and the count reaching the
observed value.

### C5b — F13's SIGN off-map, on untouched data
`ρ(gen_pc20, A_g)` over the same rows.
**PASS iff the file-block 95 % interval lies entirely above 0.** `breadth` is reported the same way, not in the rule.

### Strata, declared now
Modality: **CRISPRi** = Adamson, Dixit, Tian day-7 (restricted); **Cas9 KO** = Datlinger ×2, Papalexi ×2, Shifrut.
Both verdict statistics are reported separately for each stratum and pooled. Tian day-7 is flagged as
same-study-as-a-map-bank and the pooled result is given **with and without it**. The verdict is the POOLED cell
including Tian (the map banks from that study are excluded from its consensus).

### Multiplicity
**Two verdict cells: C5a and C5b.** Everything else in this file and in `RESULTS.md` is secondary and labelled.

---

## 4. Known-answer rows — all before any reserve number is computed

| row | requirement | status when this file was hashed |
|---|---|---|
| **KA-1** | every R3 file's sha256 equals `reserve292.json` | **PASS**, 8/8 (`work/meta295.json`) |
| **KA-2** | 279's `c_held` re-read from `held279.npz`: n = 17,106, median +0.018605 | **PASS** (+0.01860516, n 17,106) |
| **KA-3** | 291's C4 cell rebuilt from its saved predictions `work/pred291_gwcd4i.csv`: Stim48hr ρ(`C_held`, `logfc_pearson_r`) = +0.2739 at n = 1,773, and `gen_pc20` +0.1858, `breadth` +0.1275 | **PASS**, exact to four decimals |
| **KA-4** | the map consensus construction reproduces `held279.npz` `Ccross`: the pairwise cosines of the stored per-bank unit rows equal `Ccross` on the reserve genes | **PASS** — 6,097 (pair, gene) cosines, max abs diff **0.0e+00** |
| **KA-5** | inside the access: the computed per-cell library size equals `obs['ncounts']` for every cell of every file (this is what proves `X` is the raw count matrix and that the row/column indexing is right) | to be measured |
| **KA-6** | inside the access: a controls-vs-controls split reads no signal — with the file's control cells split in half, `prof(cA) − prof(cB)` has a median cosine ≈ 0 against the map consensus rows of that file's genes | to be measured |

**Abort on any failure**: if KA-5 or KA-6 fails the run stops and the ledger records the abort, not a claim value.

---

## 5. Descriptive, post-hoc, labelled (not verdicts)

The best- and worst-transferring reserve genes; whether the Datlinger 2017 → 2021 same-line/same-library pair and
the Papalexi pooled → arrayed pair agree on `A_g` (the construct / guide-assignment readings R292 §2.3 named);
per-file medians; the disattenuated `A_g`; the (file, gene)-pair version of both verdicts.

## 6. What this claim does not do

No model, dose, threshold or feature is selected on R3, and nothing computed here is competition-facing (R3 holds 6
of the 300 panel targets and none of its lines is A/B/C). No build, no submission, no `vcc`, no `cell-eval2`, no
download. The two earlier reserves (HEK293T half-B; the four Broad-16 lines) are untouched, as is the sealed IGVF
HCASMC asset.

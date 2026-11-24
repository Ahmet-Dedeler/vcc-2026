# 225 — The space is swept. **DLD-1 is the only source that clears the gate alone, and sources SATURATE at ~1.8×.** The "they added datasets one at a time" reading of 215 is dead; a four-source stack is still worth **1.9× the gate** and should ride the next slot.

**Date:** 2026-09-04 · deterministic (seeds `20726`, `771103`) · **no `vcc` command, no `cell-eval2`,
no build, nothing under `submissions/`, no production file modified, nothing outside
`experiments/225-source-sweep/` written.** Peak RSS of a single process **5.7 GiB**; two runs
briefly overlapped at ~11.2 GiB combined (mandate cap 25 GB, and `cell-eval2` was verified absent
before every launch). Network: **13 files, 39.0 MB — 0.065% of INBOX-1's 60 GB budget**, every one
hashed and recorded in `work/PROVENANCE.json`; free space rose 175 → 188 GB over the run.
**3,045 scored rows** across seven scripts; work dir 612 MB, all of it reconstructible from bytes
already on disk.

> **All four coordinator INBOX messages arrived mid-run and all four are answered in §8.**
> INBOX-CONCENTRATION and INBOX-OCTOBER each required a new measurement; both were run, and
> **INBOX-OCTOBER changed the recommendation** (§7.2, §8.3).

---

## 0. The answer in twelve lines

1. **PHASE 0 PASSES, 20/20 gates.** 223's five headline rows reproduce through 225's own harness
   to the last reported digit: `dld1_lowmoi@0.35` **+0.00228**, `@0.28` **+0.00224**,
   `telohaec@0.28` **−0.00120**, `cmap_xpr@0.28` **−0.00077**, `kolf2@0.17` **−0.00174**; and the
   three baselines match 218's `resid218`/`compose218` rows at **Δ = 0.0e+00 exactly** on all three
   cancer folds. The known answer also reproduced **bit-identically across three separate
   processes**.
2. **The census is closed.** 83 assets scanned from their own bytes, plus GEO (192 human CRISPR
   single-cell series, 2025–2026), BioStudies/ArrayExpress (191 hits), CZI CELLxGENE (**388 public
   collections, ZERO CRISPR perturbation screens**), Zenodo (which 166 could not reach), Hugging
   Face, the IGVF portal (**98 released files across four content types**) and the Arc Virtual Cell
   Atlas. **Nothing panel-covering exists that we do not already hold.** X-Atlas/Pisces is still an
   empty repo (`usedStorage: 0`, unchanged since 2026-03-25, re-checked today).
3. **One genuinely new source was found, and it was on this disk all along.**
   `data/external/igvf-cardiomyocyte/IGVFFI0830FXFI.csv.gz` carries **87/300 panel targets** —
   more than DLD-1's 67, and **the largest partial-coverage source that had never been through the
   residual arm at all** (`kolf2` 273, `feng_ipsc` 182, `cmap_sh_cgs` 122 were all priced by 223).
   081 tested this asset in Aug 2026 through a *retrieval-gallery* adapter and closed it at
   −0.00073; **it had never been through the post-whitening residual arm**, which is exactly 208's
   rider rule — re-price a source against the CURRENT shipped stack, not the one it was measured on.
4. **It is real and it is far too small.** It passes **all three falsifiers** (F1 +0.718% vs shuffle
   −0.293%/−0.206%; F2 0.0015 pp; F3 2.30× its row-matched CD4 control) and prices at **+0.00038 =
   0.20× the shape-only gate**, with Δ`pds` positive on **1/3** folds. Its per-fold pattern gives it
   away: **+1.50%/+0.0023 on `orion_hek293t` and negative on the other two**. It is a
   single-recipient effect, not a source. Cause is visible in the table: it is a **significant-hits-
   only** deposit with a **median of 6 non-zero genes per panel row** against DLD-1's dense rows.
5. **Fifteen other census members and variants were priced and none clears the gate.** Best of the
   rest: `jiang_5line_noK562` **+0.00096 (0.50×)**, `sunshine_calu3` +0.00148 (0.78×) **which fails
   F1 — its own target-shuffled control scores higher**, `norman_k562_flip` +0.00031, `myllia_a549`
   +0.00027, `xucao_hek293` +0.00010. Everything else is inside noise or negative.
6. **A leak was found and repaired rather than argued away.** Jiang 2025's six lines **include
   K562**, which is the held line on the `k562_gwps` fold. Dropping that one column
   (`jiang_5line_noK562`) *raises* the price from +0.00079 to **+0.00096** — the contamination was
   costing, not flattering. Every stack number below uses the repaired table.
7. **PHASE 4, and it is the headline. Sources DO accumulate, but the curve is strongly concave and
   it saturates at ~1.8×.** Cumulative price over 222's arm: **+0.00228 (1 source) → +0.00315 →
   +0.00354 → +0.00376 (4 sources, 1.65×) → +0.00420 (all thirteen, 1.84×).** Each added source
   buys **40–45% of the previous one**; the fourth buys **10%** of what the first bought, and the
   remaining **nine together buy +0.00044 — 2% each**. A geometric series with that ratio has
   asymptote **+0.0040**, and thirteen sources already deliver **+0.00420. We are at the asymptote.**
8. **Take DLD-1 out and the whole ladder collapses.** Twelve sources without it price at
   **+0.00148 = 0.78× the gate (LOFO 0.26×)** and go `pds`-negative on a fold. Removing DLD-1 costs
   the four-member stack **1.02× of DLD-1's own value**; removing any *other* member costs at most
   **0.45×**. **This channel is one source times a 1.8× multiplier, not thirteen sources.**
9. **So the mandate's own decision rule fires on the negative branch, and not narrowly.** "If five
   sources give 5× DLD-1 the leaders' shape is explained… if they saturate at 1.5×, the hypothesis
   is dead and 224's reading wins." **Measured 1.65× at four, 1.84× at thirteen — and the whole
   channel is worth +4.6–5.5% of board cosine against `cqawesome`'s +106% climb. A factor of
   twenty.** Five public sources cannot be what those five steps were. **Hand the channel to 224.**
10. **But the stack is still the best arm on the table and it should ride.** `dld1_lowmoi ∪
    jiang_5line_noK562 ∪ igvf_cardiomyocyte ∪ myllia_a549` at `g ≈ 0.28` prices at **+0.00376 =
    1.98× the shape-only gate, +0.00353 = 1.86× on a leave-one-fold-out dose**, 3/3 folds positive
    on both cosine and `pds`, F1 and F3 passed at every dose — against DLD-1 alone at 1.20×/1.07×.
    **+0.00147 over 223's recommendation, for zero new data and zero new compute.**
11. **INBOX-CONCENTRATION, run on 223's DLD-1 as asked: its Δ`pds` gain IS concentrated** — top-3
    of ~290 targets carry **82%**, leave-top-3-out retains **25.2%** (cosine retains 54.5%). Better
    than the shipped `h1union`'s 14.6%, so a discount not a veto — **but the stack's cosine
    retention is flat at 54.5% → 56.1% → 56.6% as it grows from 1 to 4 to 13 sources**, which is
    what a distributed additive gain looks like and is a second reason to prefer the stack.
    **And a caveat for the whole project: four unrelated sources concentrate their Δ`pds` on the
    same targets (`DEAF1`, `ZNF324`, `ZNF22`, `IFNAR2`) — that is the instrument, not the sources.**
12. **INBOX-OCTOBER, and it inverts the ranking. `kolf2` (273/300 today) and `feng_ipsc` (182/300)
    project to EXACTLY ZERO on the October panel** — they are panel-*fitted* libraries.
    **`igvf_cardiomyocyte` collapses 87 → 36. DLD-1 goes UP, 67 → 73.** And **four screens written
    off as "0/300, training-only" — `nadig_hepg2`, `nadig_jurkat`, `rpe1`, `k562_essential` —
    become 30–33/300**, because today's panel excludes common-essentials by construction and
    October's does not. **They are already built and on disk and cannot be priced until the panel
    lands.** Three DEF-RUNBOOK known-answer rows reproduce exactly (150 / 291 / 293).

---

## 1. Known-answer rows — twenty gates, all pass

Rule 10. `phase0_225.py` → `work/phase0_225.json`. Every gate **aborts** the run on failure; the
script raises `SystemExit` if any is not PASS.

| # | gate | expected | measured | verdict |
|---|---|---|---|---|
| KA-1 | `f0` == 218 `resid218` `f0`, ×3 folds | Δ = 0 | `dcos 0.0e+00 dpds 0.0e+00` | **PASS ×3** |
| KA-1 | `cd4@0.34` == 218 `cd4_gn0.34`, ×3 folds | Δ = 0 | `0.0e+00` | **PASS ×3** |
| KA-1 | `base222` == 218 `cd40.34+h1PAIRCARRY0.28`, ×3 | Δ = 0 | `0.0e+00` | **PASS ×3** |
| KA-2 | 223 `dld1_lowmoi` @ `g = 0.35` | rel +1.405%, Δ`pds` +0.00256, price **+0.00228** | **+1.405% / +0.00256 / +0.00228** | **PASS** |
| KA-2 | 223 `dld1_lowmoi` @ `g = 0.28` | +1.403% / +0.00250 / +0.00224 | **identical** | **PASS** |
| KA-3 | 223 death `telohaec` @ 0.28 | −0.410% / −0.00145 / **−0.00120** | **identical** | **PASS** |
| KA-3 | 223 death `cmap_xpr` @ 0.28 | +1.048% / −0.00138 / **−0.00077** | **identical** | **PASS** |
| KA-3 | 223 death `kolf2` @ 0.17 | +0.544% / −0.00256 / **−0.00174** | **identical** | **PASS** |
| KA-4 | silent-zero detector live: `dld1` (22%) under the OLD median-over-all-rows rule | fraction exactly 0 | `0.00000`, ×3 folds | **PASS ×3** |
| KA-4 | `kolf2` (91%) under the same rule must **not** be void | fraction > 0 | > 0, ×3 folds | **PASS ×3** |

**KA-4 is the mandate's sizing trap, asserted rather than assumed.** Every one of the fifteen
sources priced here reads `f_per_unit = 0.0` under production's median-over-**every**-row rule —
`OLD 0.000` on all 1,134 sweep lines — so all of them would have scored **exactly zero and looked
dead** on the shipped instrument. The sweep uses `solve_fraction_carrying` throughout and asserts
`median_row_norm_carrying > 0` on every arm that claims a nonzero dose (F5, never fired).

**And it reproduces across processes.** `dld1_lowmoi @ 0.35` was scored independently in
`sweep225.json` and in `sweep225b.json` (different launch, different source list). The per-fold
`d_cos` and `d_pds` are **bit-identical to 12 decimal places on all three folds**.

Baselines, for anyone re-deriving these numbers:

| fold | `f0` cos / pds | `cd4@0.34` | `base222` (= 222's arm) |
|---|---|---|---|
| `k562_gwps` | 0.028213 / 0.781085 | 0.036352 / 0.795881 | **0.039436 / 0.805397** |
| `orion_hct116` | 0.023579 / 0.782055 | 0.029689 / 0.790140 | **0.032080 / 0.792213** |
| `orion_hek293t` | 0.020866 / 0.773842 | 0.024559 / 0.760550 | **0.026157 / 0.763623** |

---

## 2. PHASE 1 — the census

### 2.1 Method, and why it differs from 166's

166 swept the same universe against a **≥100-target bar**, because it was hunting a fifth
*panel-covering donor line*, and concluded — correctly, at that bar — that the universe is
exhausted. 223 then priced a **67-target** source and found it live. So the bar was wrong for the
residual channel, and **the sub-100 class had never been enumerated at all**. 225's census reports
every asset with **≥1** panel target, counted at the **≥20-cells-per-group** threshold the donor
builder actually uses — never at the deposit's advertised library size, which is always larger.

`census225.py` → `work/census225.json`. **83 assets scanned from their own bytes; 42 carry ≥1
panel target.**

### 2.2 The full census — every asset with ≥1 panel target

Panel counts are computed by 225 from the file, not quoted from a paper. "priced" means it went
through the residual arm in §3.

| source | lineage / line | modality | panel /300 | cells | status | 225 verdict |
|---|---|---|---:|---:|---|---|
| `ReplogleWeissman2022_K562_gwps` | myeloid, K562 | CRISPRi | **272** | 1.99M | **SHIPPED donor** | in the stack already |
| `orion_hek293t` | kidney, HEK293T | CRISPRi | **294** | — | **SHIPPED donor** | — |
| `orion_hct116` | colorectal, HCT116 | CRISPRi | **290** | — | **SHIPPED donor** | — |
| `kolf2` | iPSC | CRISPRi | **273** | — | priced by 223 | **−0.00174**, F3 fail |
| CD4+ T (GSE314342) | primary CD4 T | CRISPRi | ~284 | — | **SHIPPED residual** | the incumbent `cd4f034` |
| `feng_ipsc` | iPSC, 24–34 donors | CRISPRi | 182 | — | priced by 223 | **−0.00063**, F3 fail everywhere |
| `cmap_sh_cgs` | 9 cancer lines | shRNA (L1000) | 122 | — | priced by 223 | **−0.00234** |
| **`igvf_cardiomyocyte`** | **WTC11 cardiomyocyte** | **CRISPRi** | **87** | 1.05M | **PRICED HERE — new** | **+0.00038, 0.20× gate** |
| `cmap_xpr` | 9 cancer lines | CRISPR-KO (L1000) | 82 | — | priced by 223 | **−0.00077**, F3 fail |
| **`dld1_lowmoi`** | **colorectal, DLD-1** | **CRISPRi** | **67** | — | 223's live source | **+0.00228, 1.20× gate** |
| `telohaec` | aortic endothelial | CRISPRi | 65 | — | priced by 223 | **−0.00120** |
| `gse295214_h1` | hESC | CRISPRi | 64 | — | **in 222's union** | baseline, not a candidate |
| HUES8 ESC / HUES8 DE / WTC11 hepatocyte | 3 lineages | CRISPRi | **91** | 130k–1.05M | **NO EFFECT TABLE EXISTS** | §2.4 |
| HCASMC (IGVF) | coronary smooth muscle | CRISPRi | 27 | — | **SEALED by 126** | **not opened** — §2.5 |
| `altair_a549` | lung adeno, A549 | OPS | 20 | 1,001 | priced by 223 | **−0.00004**, F1 fail |
| `h1_2025` | hESC | CRISPRi | 20 | — | **in 222's union** | leaks 85/85 of h1-test |
| VCC2025 train / test / val | H1 hESC | CRISPRi | 13 / 8 / 4 | 221k/171k/99k | **in 222's union** | — |
| `myllia_a549` | lung adeno, A549 | CRISPRi | 14 | — | priced by 223 + here | **+0.00027, 0.14×** |
| **`jiang_6line` / `jiang_5line_noK562`** | 6 cancer lines | CRISPRi | **9** | — | **PRICED HERE — new** | **+0.00079 / +0.00096, 0.50×** |
| **`sunshine_calu3`** | lung, Calu-3 | CRISPRi | **6** | 90k | **PRICED HERE — new** | **+0.00148 but F1 FAIL** |
| **`tian_ipsc_neuron`** | iPSC-derived neuron | CRISPRi | **6** | 32k | **PRICED HERE — new** | **−0.00005**, F1 fail |
| `frangieh_melanoma` | patient melanoma | Cas9 KO | 5 | 218k | priced by 223 + here | **−0.00004** |
| **`xucao_hek293`** | kidney, HEK293 | CRISPRi | **5** | 98k | **PRICED HERE — new** | **+0.00010**, F1 fail |
| **`tian_ipsc_neuron_crispra`** (± flip) | iPSC neuron | **CRISPRa** | **4** | 21k | **PRICED HERE — new** | **+0.00024 / −0.00015** |
| **`norman_k562`** (± flip) | myeloid, K562 | **CRISPRa** | **3** | 111k | **PRICED HERE — new** | **−0.00032 / +0.00031**, F1 fail |
| Hon WTC11-benchmark (IGVF) | WTC11 iPSC | CRISPRi | **3** | — | **FETCHED HERE** | 60 targets only; not built |
| **`gasperini_k562`**, **`tian2019_ipsc`**, **`wessels_thp1`** | K562 / iPSC / THP-1 | CRISPRi | **1** each | 41k–276k | **PRICED HERE — new** | ≈ 0, all F1 fail |

**The zero-panel tail, enumerated so it never has to be looked at again.** Thirty-eight further
h5ads were scanned and carry **0/300** at ≥20 cells: Joung 2023 TF atlas (1.15M cells, 2,411
usable targets), sci-Plex 2/3/4, Replogle K562-essential and RPE1, Schiebinger, Gasperini
at-scale/highMOI, Lara-Astiaso ×3, McFarland, Zhao, Schraivogel TAP chr8/chr11, Aissa, Santinha,
Cui, Adamson ×3, Weinreb, Lotfollahi, Shifrut, Dixit, Datlinger 2017/2021, Chang, Liang, Papalexi
×4, Gehring, Xie. **166's structural explanation is confirmed on 38 independent files: the VCC
panel excludes common-essential genes by construction, so every essential-library screen reads
exactly 0 in any tissue.**

### 2.3 Registries swept, each closed with one number

| registry | swept today | result |
|---|---|---|
| **CZI CELLxGENE** | 388 public collections via the curation API | **ZERO CRISPR perturbation screens.** 20 collections match "perturb/crispr/screen" on text; every one is observational or a mouse-KO tissue atlas. The single chemical entry is sci-Plex, which we hold. |
| **GEO** | **192** human `Perturb-seq / CROP-seq / CRISPRi` GSEs, 2025–2026, pulled by esummary and ranked by sample count | Top four are GSE314342 (CD4, **have**), GSE295214 (**have**), GSE291147 (melanoma PerturbFate), GSE302335 (high-MOI methodology). **No genome-scale screen in a line we do not hold.** `GuEST-List` returns **0 hits** — Satija's Feb-2026 VIPerturb-seq genome-wide library has **no GEO deposition**. |
| **BioStudies / ArrayExpress** | 191 hits | PMC mirrors plus four E-MTAB records: E-MTAB-15032 (lncRNA CRISPRi in hESC), E-MTAB-14567 (alternative promoters, 166 already priced), E-MTAB-11999/12005 (organoid CROP-seq). **All non-coding or organoid; 0/300 by construction.** |
| **Zenodo** | newest 25 dataset records for `Perturb-seq`/`CRISPRi screen` — **the sweep 166 was WAF-blocked out of** | Three 2025–2026 records: the `perturbmatch` Replogle reanalysis (**have**), a CD4 CRE/promoter screen (non-coding), and a cardiac-fibroblast bulk+sc screen. **Nothing new.** 166's one gap is now closed. |
| **Hugging Face** | `perturb` + `crispr`, 120 datasets by last-modified | Only reprocessings of what we hold: `nicolas-lynn/vcc-perturb` (= VCC2025 H1), `theislab/chem-perturbridge` (= L1000, and chemical, so no gene targets), `ChanghaoKan/crispr-depmap` (fitness). **X-Atlas/Pisces re-checked today: `usedStorage: 0`, 3 files, `lastModified 2026-03-25`. Still empty.** |
| **IGVF Data Portal** | **98 released files** across `global differential expression per element`, `differential expression quantifications`, `element quantifications`, `guide RNA sequences` | **`differential expression quantifications` does not exist in the portal at all.** The *only* gene-level effect tables IGVF has released are **16 files** of `global differential expression per element`; 11 were fetched here and 4 are 126-sealed. §2.4. |
| **Arc Virtual Cell Atlas** | repo commits | Most recent commit **2026-08-28**, README-only. No new perturbation atlas since the 2026-01-12 scBaseCount release. |
| **DepMap 26Q1 / CCLE / Tahoe** | on-disk contents | DepMap is gene-effect *fitness*, CCLE is basal TPM, Tahoe's `pdde` parquet is keyed on `drug` — **none carries a gene-target × gene effect vector, so all three are 0/300 as residual sources by construction.** |

### 2.4 The 91/300 hepatocyte prize: it is not an access problem any more, it is an absence

223 named the HUES8 / definitive-endoderm / WTC11-hepatocyte arms (**91/300**, 131's audit) as "the
single largest unmeasured partial-coverage source… an access problem, not a science problem." 225
resolved it, and the answer is worse than 223 assumed.

- The IGVF portal has **exactly 16 released gene-level DE files**, enumerated above. **None of them
  is a Huangfu or Gersbach differentiation arm.** Those arms have released *guide libraries* (43
  Gersbach, 3 Huangfu) and *QC summaries*, and nothing downstream of them.
- The `IGVF/tf_perturb_seq` GitHub repo (1,057 blobs, walked in full) holds exactly **one**
  trans-DE table: `Hon_WTC11-benchmark`. **It was fetched and measured: 60 targets, 3 on panel.**
  The Huangfu HUES8 ESC/DE directories hold only energy-distance summaries and
  `target_by_target_matrix.csv` — a target×target distance matrix, **not a gene-level effect
  table**.
- So the 91/300 asset **has not been computed and released in any public form**. Certified Synapse
  access would get the raw H5MU (17–30 GB) and require computing the DE ourselves. **Given that the
  87/300 cardiomyocyte arm of the same library and the same consortium prices at 0.20× the gate,
  that is not worth 30 GB and a certification.** 131's NO-GO stands, now for a stronger reason.

### 2.5 The 126 seal was honoured

Experiment 126's `VALUE_ACCESS_BOUNDARY.v1.json` (`permanent: true`, `value_access_authorized:
false`) forbids decompression, header inspection, row counting and every value-level operation on
`IGVFFI0206LUDV`, pending an additive authorization 225 does not hold. **225 did not fetch, open,
decompress or parse that file, and did not fetch the three sibling HCASMC global-DE files
(`IGVFFI6376HTIF`, `IGVFFI8474JFZF`, `IGVFFI5886OWOF`) either** — taking them would have been a
reinterpretation of the boundary rather than a compliance with it. The 27/300 figure quoted in the
census is read from 126's own unsealed `VCC_TARGET_OVERLAP.v1.txt`. The exclusion is recorded in
`work/PROVENANCE.json`.

---

## 3. PHASE 2 — every fetchable member priced, marginal over 222's arm

`sweep225.py` → `work/sweep225.json` (1,134 rows) · `work/sweep225b.json` (the leak-repair run).
Frame: **`base222` = `cd4@0.34 + h1pair@0.28`**, so every number is the **marginal** value on top
of the arm 222 is shipping (208's rider rule; mandate falsifier 4). Sizing is the carrying-row
rule everywhere; CD4's own dose and sizing are never touched. Primary endpoint macro cosine
(±0.20% relative, 204); price is 218 §7.4's rate card under §8.4's corrected rule, **all six
members**, `k = 1`.

### 3.1 The table — best dose per source, three cancer folds

| source | panel | carry | `g*` | rel Δcos | Δ`pds` | Δ`q_w` | +cos | +pds | **price** | **× gate** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **`dld1_lowmoi`** *(known answer)* | 67 | 63 | **0.35** | **+1.405%** | **+0.00256** | +0.00016 | **3/3** | **3/3** | **+0.00228** | **1.20×** |
| `sunshine_calu3` | 6 | 6 | 0.45 | +0.264% | **+0.00190** | +0.00006 | 3/3 | **3/3** | +0.00148 | 0.78× |
| `jiang_5line_noK562` | 9 | 9 | 0.28 | +0.660% | +0.00109 | +0.00026 | 2/3 | 3/3 | +0.00096 | 0.50× |
| `jiang_6line` *(leaked)* | 9 | 9 | 0.28 | +0.628% | +0.00086 | +0.00026 | 2/3 | 3/3 | +0.00079 | 0.41× |
| **`igvf_cardiomyocyte`** | **87** | **86** | 0.10 | **+0.718%** | +0.00027 | +0.00000 | **3/3** | 1/3 | +0.00038 | 0.20× |
| `norman_k562_flip` | 3 | 3 | 0.45 | −0.072% | +0.00044 | −0.00002 | 1/3 | 3/3 | +0.00031 | 0.16× |
| `myllia_a549` | 14 | 12 | 0.28 | +1.308% | −0.00007 | +0.00008 | 3/3 | 2/3 | +0.00027 | 0.14× |
| `tian_ipsc_neuron_crispra` | 4 | 3 | 0.35 | +0.007% | +0.00033 | +0.00015 | 2/3 | 2/3 | +0.00024 | 0.12× |
| `gasperini_k562` | 1 | 1 | 0.45 | −0.037% | +0.00017 | +0.00002 | 1/3 | 1/3 | +0.00011 | 0.06× |
| `xucao_hek293` | 5 | 4 | 0.28 | +0.029% | +0.00013 | +0.00005 | 2/3 | 3/3 | +0.00010 | 0.05× |
| `wessels_thp1` | 1 | 1 | 0.28 | +0.015% | +0.00000 | +0.00001 | 2/3 | 0/3 | +0.00000 | 0.00× |
| `frangieh_melanoma` | 5 | 5 | 0.10 | +0.025% | −0.00006 | +0.00002 | 2/3 | 1/3 | −0.00004 | −0.02× |
| `tian_ipsc_neuron` | 6 | 5 | 0.10 | −0.034% | −0.00005 | −0.00002 | 0/3 | 1/3 | −0.00005 | −0.02× |
| `tian2019_ipsc` | 1 | 0–1 | 0.10 | −0.006% | −0.00009 | +0.00000 | 0/3 | 0/3 | −0.00007 | −0.04× |
| `tian_ipsc_neuron_crispra_flip` | 4 | 3 | 0.10 | −0.035% | −0.00019 | −0.00001 | 1/3 | 1/3 | −0.00015 | −0.08× |
| `norman_k562` *(as measured)* | 3 | 3 | 0.10 | −0.010% | −0.00042 | +0.00000 | 2/3 | 0/3 | −0.00032 | −0.17× |

**`dld1_lowmoi` is the only source above the gate alone, and by 2.4× over the next survivor.**

### 3.2 The CRISPRa sign test, which came out the way the biology says it should

Two of the census members are **activation**, not knockdown, screens — the opposite intervention to
what every recipient fold measures. Both were built as measured and as a sign-flipped twin, so
which sign helps (if either) is a measurement:

- **`tian_ipsc_neuron_crispra` is positive as measured (+0.00024) and its flip is negative
  (−0.00015).** So the CRISPRa response does *not* transfer as an inverted knockdown here.
- **`norman_k562` is the reverse**: negative as measured (−0.00032), positive flipped (+0.00031).

They disagree, both are inside 0.17× of the gate, and **`norman_k562_flip` fails F1** (its
target-shuffled control scores *higher*, +0.119% against −0.072%). **The CRISPRa class is dead in
both signs**, which retires it as a family and not one screen at a time.

### 3.3 The IGVF cardiomyocyte, in detail — the source this experiment was built to find

| fold | `g = 0.10` | `g = 0.17` | `g = 0.28` | `g = 0.35` | `g = 0.45` |
|---|---|---|---|---|---|
| `k562_gwps` | +0.440% / −0.00018 | +0.304% / −0.00084 | −0.246% / −0.00163 | −0.665% / −0.00219 | −1.299% / −0.00330 |
| `orion_hct116` | +0.604% / −0.00002 | +0.605% / −0.00046 | +0.216% / −0.00133 | −0.147% / −0.00212 | −0.729% / −0.00338 |
| **`orion_hek293t`** | **+1.110% / +0.00100** | **+1.428% / +0.00156** | **+1.504% / +0.00234** | **+1.398% / +0.00273** | **+1.157% / +0.00302** |

**It is one fold.** On `orion_hek293t` it looks like a genuine, well-dosed source with a clean
interior optimum; on the other two it is negative on `pds` at every dose. That is the signature of
a recipient-specific match, not of shared perturbation direction, and it is why the three-fold
price collapses to 0.20× the gate. **Two independent structural reasons it cannot do better:**

1. **The deposit is significant-hits-only.** Median **6** non-zero genes per panel row (mean 18.1,
   max 217) out of 8,236 columns — 1,575 measured coordinates over 87 rows. DLD-1's rows are dense.
   A residual row with six coordinates is almost all zero, and its effect-norm median on the panel
   is **2.65** against DLD-1's much larger rows.
2. **No independent library exists to replicate it with.** The other released file from the same
   analysis set, `IGVFFI5989UAVX`, is **97% identical** to the one on disk (26,712 of 27,483
   matched (target, gene) pairs agree exactly on `fc` and `num_cell`) — it is a re-release with
   extra annotation columns, not a second library. So the DLD-1-style two-library falsifier that
   made 223's case **cannot be run on this source at all**, and its single-fold pattern is the only
   evidence available. That is a reason for less confidence, not more.

---

## 4. PHASE 3 — the four mandatory falsifiers, every source

### F1 — target-shuffled control (218's derangement: geometry, per-row norms, support pattern and carrying count preserved exactly). Two seeds, `20726` and `771103`.

At each source's own price-argmax dose, three cancer folds, marginal over 222:

| source | `g*` | real | SHUF (seed 1) | SHUF2 (seed 2) | **F1** |
|---|---:|---:|---:|---:|---|
| **`dld1_lowmoi`** | 0.35 | **+1.405%** | −1.216% | −1.954% | **PASS** |
| **`myllia_a549`** | 0.28 | **+1.308%** | −0.113% | −0.391% | **PASS** |
| **`igvf_cardiomyocyte`** | 0.10 | **+0.718%** | −0.293% | −0.206% | **PASS** |
| `jiang_5line_noK562` | 0.28 | +0.660% | **+0.597%** | +0.338% | **weak pass — margin 0.06 pp, inside 204's ±0.20% resolution** |
| `tian_ipsc_neuron_crispra` | 0.35 | +0.007% | −0.447% | −0.498% | pass on sign, arm ≈ 0 |
| **`sunshine_calu3`** | 0.45 | +0.264% | **+0.286%** | **+0.421%** | **FAIL — its own shuffle scores higher** |
| `norman_k562_flip` | 0.45 | −0.072% | +0.119% | +0.119% | **FAIL** |
| `xucao_hek293` | 0.28 | +0.029% | +0.031% | −0.058% | **FAIL** |
| `wessels_thp1` | 0.28 | +0.015% | +0.015% | +0.015% | **FAIL — identical to its shuffle** |
| `tian_ipsc_neuron`, `tian2019_ipsc`, `gasperini_k562`, `norman_k562`, `tian_..._crispra_flip` | — | ≤ 0 | — | — | **FAIL** (arm negative) |

**F1 is what kills `sunshine_calu3`, and it is the most useful single result in this section.**
Ranked on the coordinator's own Δ`pds`-first filter (INBOX 1 §1), Calu-3 looked like the sweep's
second find: **Δ`pds` +0.00190, positive on 3/3 cancer folds at every one of five doses,
monotone in dose, price 0.78× the gate — from six panel rows.** It is not real. Permute which
target each of those six rows belongs to and the arm scores **better** (+0.286%/+0.421% against
+0.264%). Its Δ`pds` is a property of adding *any* six-row residual of that norm on those rows, not
of Calu-3's information. **Without a two-seed shuffle control this experiment would have reported a
false positive**, and it is exactly the failure mode 204 and 209 died of.

### F2 — non-carrying rows must not move (218's bar: |Δ| ≤ 0.018 pp)

Maximum excursion over the three cancer folds at each source's `g*`: **every source passes**, the
largest being `dld1_lowmoi` at **0.0125 pp** (`g = 0.35`) and `myllia_a549` at **0.0050 pp**;
`igvf_cardiomyocyte` reads **0.0015 pp**. Against carrying-row moves of **+5.1%** and **+4.2%**
that is two to three orders of magnitude apart. The arms add information on the rows they carry;
they do not reweight rows.

### F3 — the row-matched CD4 control (same rows, same sizing rule, same dose, CD4's information)

| source | `g*` | real | row-matched CD4 | margin | **F3** |
|---|---:|---:|---:|---:|---|
| **`myllia_a549`** | 0.28 | +1.308% | +0.259% | **5.04×** | **PASS** |
| **`dld1_lowmoi`** | 0.35 | +1.405% | +0.401% | **3.50×** | **PASS** |
| **`igvf_cardiomyocyte`** | 0.10 | +0.718% | +0.313% | **2.30×** | **PASS** |
| `jiang_5line_noK562` | 0.28 | +0.660% | **−0.248%** | real beats a *negative* control | **PASS** |
| `sunshine_calu3` | 0.45 | +0.264% | +0.172% | 1.54× | passes F3, but **F1 already failed** |
| `xucao_hek293` | 0.28 | +0.029% | +0.089% | **0.32×** | **FAIL** |
| `frangieh_melanoma` | 0.10 | +0.025% | +0.031% | 0.81× | **FAIL** |
| `norman_k562_flip`, `norman_k562`, `tian_ipsc_neuron`, `tian2019_ipsc`, `tian_..._flip` | — | ≤ 0 | — | — | **FAIL** |

> **One correction to 223's convention, stated because it changes a verdict.** 223 reported F3 as a
> *ratio* and every one of its controls happened to be positive. `jiang`'s row-matched CD4 control
> is **negative** (−0.248%), so the ratio flips sign (−2.5×) and reads as a catastrophic failure
> when the substantive comparison — *does the source beat the free alternative on those rows?* — is
> a clear pass. 225 scores F3 as a **difference**, and reports the ratio only for context.

### F4 — leak audit, every arm, stated explicitly

**No source is ever scored on the fold where it is the truth.** None of the fifteen is one of the
four held lines. Carrying rows per source per fold:

| source | `k562_gwps` | `orion_hct116` | `orion_hek293t` |
|---|---:|---:|---:|
| `igvf_cardiomyocyte` | 86/269 | 86/290 | 86/294 |
| `dld1_lowmoi` | 60/269 | 66/290 | 65/294 |
| `myllia_a549` | 13/269 | 12/290 | 13/294 |
| `jiang_6line` / `_5line_noK562` | 9/269 | 9/290 | 9/294 |
| `sunshine_calu3` | 6/269 | 6/290 | 6/294 |
| `tian_ipsc_neuron` | 6/269 | 5/290 | 6/294 |
| `xucao_hek293` | 4/269 | 5/290 | 5/294 |
| `frangieh_melanoma` | 5/269 | 5/290 | 5/294 |
| `norman_k562` (±flip) | 3/269 | 3/290 | 3/294 |
| `tian_ipsc_neuron_crispra` (±flip) | 3/269 | 4/290 | 4/294 |
| `gasperini_k562`, `wessels_thp1`, `tian2019_ipsc` | 1/269 | 1/290 | 0–1/294 |

**Every donor↔surface overlap that matters, stated:**

- **`h1_2025` covers 85/85 of `h1-test`'s targets** and is inside 222's shipped union, which is
  225's *baseline*. **`h1-test` is not used anywhere in 225** — every measurement is on 162's PLOLO
  challenge frame. This is the trap that cost experiment 157; it is flagged again because the
  baseline carries `h1_2025`.
- **Jiang 2025's six lines INCLUDE K562, which is the held line on the `k562_gwps` fold.** This is
  a real leak in the pooled consensus and 225 repaired it rather than noting it: `jiang_5line_noK562`
  drops that column. **The repair *raises* the price (+0.00079 → +0.00096)** — the contaminated
  column was adding noise, not signal — and every stack number in §5 uses the repaired table. The
  `k562_gwps` fold is where Jiang does *worst* (−0.28% at `g = 0.28`) either way, so no positive
  claim ever rested on the leak.
- **`xucao_hek293` is HEK293 and `orion_hek293t` is HEK293T** — near-identical lines. This would
  *flatter* Xu, and Xu fails F1 and F3 anyway, so no positive claim rests on it.
- **`norman_k562` and `gasperini_k562` are K562**, the held line on the `k562_gwps` fold. Same
  direction of bias, same outcome: both fail F1.
- **`igvf_cardiomyocyte`, `dld1_lowmoi`, `myllia_a549`, `sunshine_calu3`, `frangieh_melanoma`,
  `tian_*`, `wessels_thp1`** are none of the four held lines and share no line with any of them.
- **DLD-1 and HCT116 are both colorectal**, and 223 already showed the gain is not a lineage
  artefact (largest on the myeloid `k562_gwps` fold). 225 reproduces that ordering exactly:
  +1.646% k562 > +1.565% hct116 > +1.003% hek293t.

### F5 — non-vacuity, asserted in-loop

Every arm claiming a nonzero dose asserts `median_row_norm_carrying > 0` and the run aborts
otherwise. **Never fired.** Every arm records `median_row_norm_all_rows`,
`median_row_norm_carrying`, `carrying_rows`, `coverage_frac`, `fraction`, `per_unit`, the OLD-rule
fraction and a `VACUOUS` flag.

### F6 — two-library / multi-replicate replication

See §6.

---

## 5. PHASE 4 — the stack curve. **This is the experiment's answer.**

`stack225.py` → `work/stack225d.json` (headline) and `work/stack225.json` (the first pass, kept
for the record). Cumulative unions of the ranked survivors, each swept over six doses on all three
cancer folds, every one carrying its own F1 and F3 controls, members added in measured-price order.

**Two passes were run and only the second is quoted.** The first used `jiang_6line` (K562-leaked)
and included `sunshine_calu3`, which F1 then killed (§4). The second replaced Jiang with the
leak-repaired table and dropped Calu-3. The two agree on the finding and differ in the third
decimal: pass 1 read 1.84× at five members and 1.77× at thirteen; pass 2 reads 1.65× at four and
1.84× at thirteen. **The saturation is not an artefact of either construction.**

### 5.1 The curve

Members added in measured-price order, `jiang_5line_noK562` used throughout (the leak-repaired
table). `ALLSMALL` is the union of **every** census source with ≥1 panel target, including the ones
that are individually dead — the "many small sources" hypothesis given its very best shot, because
159's `1/k` law is about independent *estimates*, not about rows.

| stack | members | n | `g*` | carry | rel Δcos | Δ`pds` | +pds | **price** | **× DLD-1** | **× gate** | **LOFO** | LOFO × gate |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **CUM1** | `dld1` | 1 | 0.35 | 63 | +1.405% | +0.00256 | 3/3 | **+0.00228** | **1.00×** | 1.20× | +0.00204 | 1.07× |
| **CUM2** | + `jiang5` | 2 | 0.28 | 69 | +2.292% | +0.00345 | 3/3 | **+0.00315** | **1.38×** | 1.66× | +0.00315 | 1.66× |
| **CUM3** | + `igvf` | 3 | 0.28 | 126 | +2.679% | +0.00387 | 3/3 | **+0.00354** | **1.55×** | 1.86× | +0.00326 | 1.72× |
| **CUM4** | + `myllia` | 4 | 0.28 | 133 | +3.897% | +0.00376 | 3/3 | **+0.00376** | **1.65×** | 1.98× | +0.00353 | 1.86× |
| **ALLSMALL** | **all 13** | 13 | 0.35 | 147 | +3.953% | +0.00435 | 3/3 | **+0.00420** | **1.84×** | **2.21×** | **+0.00398** | **2.10×** |
| LOO −`jiang5` | `dld1+igvf+myllia` | 3 | 0.35 | 130 | +3.166% | +0.00263 | 3/3 | +0.00274 | 1.20× | 1.44× | +0.00263 | 1.38× |
| LOO −`igvf` | `dld1+jiang5+myllia` | 3 | 0.28 | 79 | +3.448% | +0.00335 | 3/3 | +0.00336 | 1.47× | 1.77× | +0.00336 | 1.77× |
| LOO −`myllia` | `dld1+jiang5+igvf` | 3 | 0.28 | 126 | +2.679% | +0.00387 | 3/3 | +0.00354 | 1.55× | 1.86× | +0.00326 | 1.72× |
| **LOO −`dld1`** | `jiang5+igvf+myllia` | 3 | 0.17 | 96 | +2.171% | +0.00122 | 3/3 | **+0.00143** | **0.63×** | **0.75×** | +0.00116 | 0.61× |
| **ALLSMALL −`dld1`** | **all 12 minus dld1** | 12 | 0.28 | 116 | +2.638% | +0.00117 | 2/3 | **+0.00148** | **0.65×** | **0.78×** | +0.00050 | 0.26× |

**Falsifiers on the stack itself.** Every row above beats its own target-shuffled control and its
own row-matched CD4 control at **every dose tested**. At `ALLSMALL`'s argmax: real **+3.953%**,
SHUF **−3.549%**, row-matched CD4 **+1.637%**. At `CUM1`'s: +1.405% / −1.216% / +0.401%. The stack
is real signal; it is the *amount* that is the finding.

### 5.2 The marginal source, in one column

| step | source added | marginal price | fraction of step 1 |
|---|---|---:|---:|
| 1 | `dld1_lowmoi` | **+0.00228** | 100% |
| 2 | `jiang_5line_noK562` | +0.00086 | 38% |
| 3 | `igvf_cardiomyocyte` | +0.00039 | 17% |
| 4 | `myllia_a549` | +0.00022 | 10% |
| 5–13 | **the other nine sources, together** | **+0.00044** | **19% total, ~2% each** |

Each added source buys roughly **40–45% of the previous one**. A geometric series with that ratio
has asymptote `0.00228 / (1 − 0.43) ≈ +0.0040` — **and thirteen sources already deliver +0.00420.
We are at the asymptote. There is no more of this channel to buy.**

### 5.3 Two independent proofs that it is DLD-1 and not "sources"

1. **Remove DLD-1 and the whole ladder collapses.** Twelve sources without it price at
   **+0.00148 (0.78× the gate, LOFO 0.26×)** and go `pds`-negative on one fold. Removing DLD-1 costs
   the stack **1.02×** of DLD-1's own standalone value; removing any *other* member costs at most
   **0.45×**. **The channel is one source plus a 1.8× multiplier, not thirteen sources.**
2. **Jiang's own six lines do not replicate each other** (§6.2). Six independently run screens of
   the same 218 targets in six different cancer lines — the purest available instance of "add
   another public dataset" — and **only one of the six is individually positive and survives F1**.

---

## 6. F6 — replication falsifiers

### 6.1 DLD-1's two independent libraries — 223's strongest evidence, reproduced through 225's harness

`sweep225c.py` run, marginal over 222, three cancer folds:

| table | `g` | rel Δcos | Δ`pds` | +pds | SHUF | row-matched CD4 | price | × gate |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `dld1_low1` alone | 0.28 | **+0.565%** | **+0.00079** | 2/3 | −0.798% | +0.464% | +0.00075 | 0.39× |
| `dld1_low2` alone | 0.28 | **+1.080%** | **+0.00267** | 3/3 | −0.890% | +0.464% | +0.00227 | 1.19× |
| **two-split mean** | 0.28 | **+1.403%** | +0.00250 | 3/3 | −0.751% | +0.464% | +0.00224 | 1.18× |
| **two-split mean** | 0.35 | **+1.405%** | +0.00256 | 3/3 | −1.216% | +0.401% | **+0.00228** | **1.20×** |

**Both libraries are independently positive on cosine and on `pds`, and on the primary endpoint the
mean beats both** (+1.403% against +1.080% and +0.565%) — the signature of averaging two noisy
estimates of one real signal. 223's F6 reproduces exactly. *(On price the mean and `low2` are
indistinguishable, +0.00224 against +0.00227; the mean's advantage is on cosine, which is the
primary endpoint, and its advantage over `low1` is 3×.)*

### 6.2 Jiang's six lines, run separately — **and this is the result that kills the hypothesis from the inside**

Jiang 2025 screens the **same 218 targets in six different cancer lines**. If "another public
dataset adds another independent estimate of the shared direction" were the mechanism behind the
leaders' climb, this is where it would show: six screens, one library, one lab, one protocol.

| line | lineage | rel Δcos @ 0.28 | Δ`pds` | +cos | +pds | SHUF | price | × gate |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `BXPC3` | pancreas | +0.689% | **+0.00071** | 2/3 | 2/3 | +0.225% | **+0.00068** | 0.36× |
| `MCF7` | breast | +0.384% | +0.00030 | 1/3 | 3/3 | **+1.065%** | +0.00032 | 0.17× — **F1 FAIL** |
| `A549` | lung adeno | +0.432% | −0.00035 | 2/3 | 1/3 | +0.344% | −0.00017 | −0.09× |
| `HT29` | colorectal | +0.463% | −0.00059 | 2/3 | 1/3 | +0.144% | −0.00033 | −0.17× |
| `K562` | myeloid | −0.167% | −0.00029 | 1/3 | 1/3 | −0.134% | −0.00024 | −0.13× |
| `HAP1` | near-haploid CML | −0.158% | −0.00091 | 1/3 | 1/3 | +0.001% | −0.00071 | −0.38× |
| **pooled 5-line consensus** | — | **+0.660%** | **+0.00109** | 2/3 | 3/3 | +0.597% | **+0.00096** | **0.50×** |

**One of six lines is positive and survives its shuffle. `MCF7` is beaten by its own target-shuffled
control by 2.8×.** The pooled consensus is worth more than any single line — which is the averaging
mechanism working — but it is worth **0.50× the gate**, and its own F1 margin (+0.660% real against
+0.597% shuffled) is **0.06 pp, inside 204's ±0.20% instrument resolution**. **`jiang` is therefore
NOT a credible standalone source**; it earns its place in the stack only as a marginal contributor,
and the honest three-member stack without it (`dld1 + igvf + myllia`) prices at **+0.00274, 1.44×
the gate**.

### 6.3 `igvf_cardiomyocyte` cannot be replication-tested, and that is a reason for less confidence

The only other released file from analysis set `IGVFDS6332VCTO` (`IGVFFI5989UAVX`, fetched and
compared here) agrees with the one on disk on **26,712 of 27,483** matched `(target, gene)` pairs —
same `fc`, same `num_cell`. It is a re-release with extra annotation columns, **not an independent
library**. No two-library split exists for this source anywhere in IGVF. Its single-fold pattern
(§3.3) is the only replication evidence available and it is negative.

---

## 7. Verdict — what this says about 215, and what should ride

### 7.1 The hypothesis this experiment was built to test is **dead**, and the number is not close

The mandate's own decision rule: *"if five sources give 5× DLD-1 the leaders' shape is explained…
if they saturate at 1.5×, the hypothesis is dead and 224's reading wins."*

**Measured: 1.65× at four sources, 1.84× at thirteen, asymptote ≈ 1.8×.** The saturating branch.

And the like-for-like comparison against 215 is a factor of twenty, **under every PLOLO→board
factor currently on the table** — which matters, because INBOX-BRANCH-A rules out converting with
a single scalar until 221 reconciles `cd4f034`'s **1.384** against `h1union`'s **1.000**:

| quantity | value |
|---|---:|
| `cqawesome`'s climb, `m0830b` → `m0903c`, 5 resolved steps | board cosine floor **0.06427 → 0.13211 = ×2.06, i.e. +106%** (215 §3.2) |
| **the ENTIRE exhaustively-swept public-source channel**, 13 sources, best stack | **+3.95% relative macro cosine, PLOLO, `cos` (per-target unweighted), folds `k562_gwps` + `orion_hct116` + `orion_hek293t`** |
| — at `h1union`'s measured factor **1.000** | +3.95% board |
| — at `cd4f034`'s measured factor **1.384** | +5.5% board |
| **ratio to their climb, across the whole factor range** | **1/27th to 1/19th** |

**The conclusion does not depend on which factor is right**, which is exactly why it is stated as a
range rather than a point. All prices elsewhere in this file are `k = 1` (no transfer factor
applied) and are labelled as such.

**Every public perturbation dataset on earth that shares a target with the 300-gene panel, unioned
together and dosed at its own optimum, moves the transfer cosine by about 5%. `cqawesome` moved
theirs by 106% in five steps.** Five public sources cannot be what those five steps were. The
"data generation is the bottleneck" argument may well be true for *them* — they have in-house
data — but it is not a route that is open to us, because the public half of it is now measured and
it is worth 5%.

**Recommendation to the coordinator: the direction channel's remaining budget should go to 224.**
215's incremental shape is not explained by public-source accumulation, and 225 has removed that
explanation with a number rather than an argument.

### 7.2 But the stack is the best arm currently on the table and it should ride the next slot

This is a negative result about a *hypothesis* and a positive result about an *arm*.

| rank | arm | dose | price (3 cancer folds) | LOFO (dose chosen out of sample) | × shape-only gate | falsifiers |
|---:|---|---|---:|---:|---:|---|
| **1** | **`dld1_lowmoi ∪ jiang_5line_noK562 ∪ igvf_cardiomyocyte ∪ myllia_a549`** | `g ≈ 0.28` | **+0.00376** | **+0.00353, 3/3 folds** | **1.86×** | F1 ✓ (2 seeds, all doses) · F2 ✓ · F3 ✓ · F4 ✓ · F5 ✓ |
| 2 | the union of **all 13** census sources | `g ≈ 0.35` | +0.00420 | +0.00398, 3/3 | 2.10× | same, but nine of its members are individually inside noise — **more moving parts for +0.00045** |
| 3 | `dld1_lowmoi ∪ igvf_cardiomyocyte ∪ myllia_a549` *(drop the F1-marginal Jiang)* | `g ≈ 0.35` | +0.00274 | +0.00263, 3/3 | 1.38× | the conservative choice |
| 4 | `dld1_lowmoi` alone — **223's recommendation** | `g ≈ 0.35` | +0.00228 | +0.00204, 3/3 | 1.07× | the incumbent proposal |
| — | anything **without** `dld1_lowmoi` | — | ≤ +0.00148 | ≤ +0.00116 | **≤ 0.78×** | **does not clear the gate** |

**Rank 1 is +0.00147 over 223's rank-1, for zero new data and zero new compute** — every table it
needs is on disk, three of the four were built in this experiment in under five minutes, and the
sizing rule it needs (`solve_carrying_row_fraction`) is already committed at
`scripts/make_submission_pm_cd4.py:1309`. **The only change is registering three more extra-residual
sources, which is 222's territory, not mine. CD4's dose and sizing are untouched throughout.**

**Two qualifications from the coordinator INBOXes, both measured in §8 and both load-bearing.**

- **Concentration (§8.2).** DLD-1's Δ`pds` gain retains only **25.2%** leave-top-3-out; its cosine
  gain retains **54.5%**. The *stack's* cosine retention is **56.1%** at four members and **56.6%**
  at thirteen — flat as the stack grows, which is the signature of a distributed, additive gain and
  is the reason to prefer the stack over the single source on more than headline price. **Rank on
  the cosine channel. Δ`pds` on PLOLO concentrates on the same handful of targets whatever the
  source (§8.2 item 3) and should not be the ranking key** — which also qualifies 223's and 218's
  Δ`pds`-first shortlists.
- **October (§8.3), and it changes what "redundant" means.** INBOX-OCTOBER's second-order point is
  correct and now measured: **the saturation curve being flat between four and thirteen sources is
  not a reason to drop the extra nine — it is October insurance, and it is nearly free.** The
  thirteen-source union carries **147** rows against the four-source union's **133** today, and it
  spreads the arm over nine more independent libraries at a cost of +0.00044 of headline price.
  Against a **disjoint** October panel, a stack leaning on one well-matched source degrades far
  worse than a broad one. **So rank 2 (all thirteen) is the better October arm even though rank 1
  is the cleaner September one**, and the choice between them is a risk decision, not a
  measurement one. My recommendation: **ship rank 1 now, and register the other nine before
  2026-10-22.**

**Classification: shape-only.** α, γ, λ and the emission law are untouched; the only renormalisation
is production's single global preserve scalar, present at every dose. It clears the shape-only gate
on every framing and clears **none** of the amplitude-class gate (+0.0089) on any framing —
`ALLSMALL`'s +0.00420 is 0.47× of it. As in 223, that classification is load-bearing and is not
argued both ways.

### 7.3 Things that are now closed and should not be re-opened

| claim | closed by | number |
|---|---|---|
| "There is an unswept public source with meaningful panel coverage" | §2.2, 83 assets scanned from their own bytes | The top of the sub-100 class is **87 (IGVF cardiomyocyte, priced 0.20×)** and **82 (`cmap_xpr`, 223 killed it)**. Below 67 the next is **20**. |
| "The 91/300 hepatocyte / HUES8 arms are an access problem" | §2.4, IGVF portal + the full 1,057-blob GitHub tree | **No gene-level effect table has been computed and released for them in any public form.** Not an access problem — an absence. |
| "X-Atlas/Pisces might have landed" | HF API, re-checked 2026-09-04 | `usedStorage: 0`, 3 files, `lastModified 2026-03-25`. |
| "CELLxGENE has perturbation screens we have not looked at" | 388 public collections | **Zero CRISPR perturbation screens.** |
| "Zenodo was never swept" (166's one acknowledged gap) | §2.3 | Swept. Three 2025–26 records, all held or non-coding. |
| "A CRISPRa atlas could serve as an inverted knockdown source" | §3.2, two screens × two signs | Both signs inside 0.17× of the gate; the two screens **disagree on which sign helps**; `norman_k562_flip` fails F1. Class retired. |
| "Adding more independent cell lines is the mechanism" | §6.2, Jiang's six lines | **1 of 6 individually positive and F1-surviving.** |
| "081's cardiomyocyte NO-GO might reverse under the corrected pricing" | §3.3 | It does not. It rises from −0.00073 to **+0.00038**, still **0.20×** the gate, one fold, **19% cosine retention** (§8.2) and **87 → 36 in October** (§8.3). Closed on four independent axes. |
| "a source's overlap with today's 300 is what it is worth" | §8.3, on 199's own two synthetic draws | **False, and by a lot.** `kolf2` 273 → **0**, `feng_ipsc` 182 → **0**, `igvf` 87 → 36, `k562_gwps` 269 → 150, while `dld1` 67 → **73** and four 0/300 screens → **30–33**. **Rank sources by breadth, not by today's overlap.** |

### 7.4 Honest limitations

- **The `kolf2` fold is excluded from every headline**, exactly as 218 and 223 excluded it
  (pluripotent recipient, z-score truth with a placeholder `b_ctrl`). All headline numbers are the
  three cancer folds.
- **The board transfer factor is not applied.** Every price is `k = 1`. §8.4's 1.17–1.39 was measured
  on arms that move ~93% of the panel's rows; these move 21–49% (63–147 carrying rows of 269–294).
  Applying it would raise `ALLSMALL` to +0.0049–0.0058 and would not change any verdict.
- **The stack's dose was chosen on the same three folds it is reported on**, so the LOFO column —
  dose selected on two folds, priced on the held third — is the honest number: **+0.00353 for rank 1,
  +0.00398 for rank 2, 3/3 folds positive for both.**
- **`igvf_cardiomyocyte` earns its place in the stack but not on its own.** It is 0.20× alone and a
  single-fold effect; in the union it adds +0.00039. If the coordinator wants the most defensible
  arm rather than the largest, rank 3 (`dld1 + igvf + myllia`, 1.38×) drops the F1-marginal member,
  and rank 4 (`dld1` alone, 1.07×) drops everything 225 added.
- **Δ`pds` is a fragile endpoint on this instrument and I used it as a co-primary anyway.** §8.2
  found four unrelated sources concentrating their Δ`pds` on the same handful of targets. Every
  Δ`pds` number in this file should be read with that discount, and the cosine column is the one I
  would defend. This does not change any verdict here — DLD-1 leads on both channels and the deaths
  fail on both — but it qualifies the *sizes*.
- **The four common-essential screens that come alive in October (§8.3 item 4) could not be priced**,
  because they carry 0 rows of today's panel and PLOLO can only score today's panel. Their October
  value is a coverage projection, not a measurement. **That is the single largest untested claim in
  this file** and it is untestable until 2026-10-22.
- **The HCASMC source (27/300) was not opened**, per 126's permanent boundary. If a future
  authorization is written, it is the one remaining unmeasured source in the census — and on the
  87-row cardiomyocyte's measured rate, 27 rows would be worth proportionally less.

---

## 8. The four coordinator INBOX messages, answered

All four arrived after this experiment launched. Two required new measurements; both were run.

### 8.1 INBOX-1 — disk budget: **used 0.065% of it**

| item | bytes | note |
|---|---:|---|
| authorised budget | **60 GB** | INBOX-1 rule 1 |
| **actually fetched from the network** | **39.0 MB** (13 files) | every one hashed in `work/PROVENANCE.json` |
| derived tables written (`work/tables/`) | 437 MB | 19 donor tables, all reconstructible from bytes already on disk |
| total work dir | 612 MB | |
| free space at start / at finish | 175 GB / **188 GB** | never approached the 100 GB abort threshold |

**No raw archive was fetched that needed deleting** — the census ranked smallest-per-panel-overlap
first (rule 3) and the answer turned out to be that **every decisive asset was already on disk**.
The largest single fetch was 15.0 MB. Rule 5's costed-recommendation case did arise once and is
answered in §2.4: **the 91/300 HUES8/hepatocyte arms would cost a Synapse certification plus a
17–30 GB H5MU plus computing the DE ourselves, and the 87/300 arm of the same library and
consortium measures at 0.20× the gate — so it is a named, sized NO-BUY, not an unfetched maybe.**

### 8.2 INBOX-CONCENTRATION — leave-top-k-out on both channels. **DLD-1's `pds` gain IS concentrated.**

`conc225.py` → `work/conc225.json`. 222's definition verbatim (per-target arrays out of
`score_effects`, leave-top-k-out by deleting those targets from the average), three cancer folds,
each channel ranked by **its own** per-target delta.

| arm | k | Δ`pds` full | Δ`pds` LOO | **`pds` retained** | top-k share | Δcos full | Δcos LOO | **cos retained** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **`dld1_lowmoi`** | 1 | +0.00256 | +0.00170 | **66.6%** | 34.7% | +0.000471 | +0.000370 | **78.6%** |
| **`dld1_lowmoi`** | **3** | +0.00256 | +0.00064 | **25.2%** | **81.9%** | +0.000471 | +0.000257 | **54.5%** |
| `CUM4` (4 sources) | 3 | +0.00376 | +0.00058 | 15.5% | 89.8% | +0.001192 | +0.000669 | **56.1%** |
| `ALLSMALL` (13) | 3 | +0.00435 | +0.00034 | 7.7% | 107.2% | +0.001188 | +0.000673 | **56.6%** |
| `igvf_cardiomyocyte` | 3 | +0.00027 | −0.00039 | **−144%** | −765% | +0.000219 | +0.000042 | **19.2%** |
| `jiang_5line_noK562` | 3 | +0.00109 | −0.00066 | **−61%** | 187% | +0.000188 | −0.000053 | **−28.5%** |
| *(222's `h1union`, for reference)* | 3 | — | — | *14.6%* | *99.8%* | — | — | *72.7%* |

**Three findings, and the third is the one that should change how this project reads Δ`pds`.**

1. **DLD-1's Δ`pds` is concentrated: the top three targets of ~290 carry 82% of it, and
   leave-top-3-out retains 25%.** It is *less* concentrated than 222's shipped `h1union`
   (14.6% retained) — which the board then scored at +0.0069 with all six members up — so on
   INBOX-BRANCH-A's rule this is **a discount, not a veto**. But it must be stated: 223's
   +0.00228 was never shown to be distributed and **it is not.**
2. **The cosine channel is the trustworthy one, and it is genuinely additive.** Retention is
   **54.5% / 56.1% / 56.6%** for one, four and thirteen sources — essentially constant as the
   stack grows. A concentrated gain would degrade as members were added; this does not. **Rank on
   cosine retention; treat Δ`pds` as the fragile channel.**
3. **🔴 The top-3 Δ`pds` targets are largely the SAME across unrelated sources.** `DEAF1` appears
   in the top-3 of `igvf_cardiomyocyte`, `CUM4` and `ALLSMALL`; `ZNF324` in `dld1`, `CUM4` and
   `ALLSMALL`; `ZNF22` and `IFNAR2` in `jiang` and the stacks. **These are four sources with
   different lineages, different labs and largely different carried rows, concentrating on the same
   handful of targets.** The INBOX names this case explicitly — *"if it concentrates on targets
   unrelated to it, suspect the instrument"* — and that is the reading here: **PLOLO's Δ`pds` has a
   few high-leverage targets that dominate any residual arm's measured gain, whatever the source.**
   This is a caveat on the *instrument*, and it applies to 223's and 218's Δ`pds`-ranked shortlists
   as much as to mine. **Handing it to 221 as an instrument question.**
4. **It sharpens two verdicts in §3 rather than changing them.** `igvf_cardiomyocyte` retains
   **19.2%** of its cosine gain and **−144%** of its `pds` gain; `jiang_5line_noK562` retains
   **−28.5%** and **−61%**. Both were already below the gate; both are now also **not distributed**.
   **Neither should ride on its own under any framing**, and their place in the stack is worth
   exactly the +0.00039 and +0.00086 the marginal column says and no more.

### 8.3 INBOX-OCTOBER — breadth, and **the ranking inverts. This is the most decision-relevant table here.**

`october225.py` → `work/october225.json`. Projected coverage under DEF-RUNBOOK §5's own two
synthetic draws, **reusing 199's bundles** (`def-bundle-r` = realistic, drawn from K562 ∩ Orion ∩
CD4; `def-bundle` = pessimistic, drawn from any donor's coverage) rather than inventing a sampling
model. Both draws are asserted disjoint from today's panel before anything is counted.

**Three known-answer rows first**, because this instrument is new: DEF-RUNBOOK §5 states
`k562_gwps` **269 → 150** on the pessimistic draw and Orion HCT116/HEK293T **→ 291 / 293**.
Measured here: **150**, **291**, **293**. Exact on all three.

| source | breadth (distinct targets) | today /300 | ∩ donor union | realistic /300 | **PESSIMISTIC /300** | **Δ vs today** |
|---|---:|---:|---:|---:|---:|---:|
| `orion_hek293t` | 17,977 | 294 | 17,977 | 300 | **293** | −1 |
| `orion_hct116` | 17,767 | 290 | 17,767 | 300 | **291** | +1 |
| `k562_gwps` | 9,787 | 269 | 9,787 | 300 | **150** | **−119** |
| `cmap_xpr` | 5,157 | 82 | 4,969 | 66 | **88** | **+6** |
| **`dld1_lowmoi`** | **4,715** | **67** | 4,506 | 67 | **73** | **+6** |
| `cmap_sh_cgs` | 4,342 | 122 | 4,103 | 88 | **66** | −56 |
| `telohaec` | 2,167 | 65 | 2,061 | 42 | **40** | −25 |
| **`igvf_cardiomyocyte`** | **1,888** | **87** | 1,818 | 60 | **36** | **−51** |
| `gse295214_h1` *(in 222's union)* | 2,354 | 64 | 2,347 | 57 | **34** | −30 |
| **`nadig_hepg2`** | 2,160 | **0** | 2,156 | 45 | **33** | **+33** |
| **`nadig_jurkat`** | 2,270 | **0** | 2,265 | 48 | **33** | **+33** |
| **`rpe1`** | 2,204 | **0** | 2,200 | 47 | **32** | **+32** |
| **`k562_essential`** | 2,021 | **0** | 2,017 | 34 | **30** | **+30** |
| `altair_a549` | 997 | 20 | 972 | 23 | 10 | −10 |
| `jiang_5line_noK562` | 218 | 9 | 216 | 3 | 4 | −5 |
| `myllia_a549` | 175 | 14 | 173 | 5 | 4 | −10 |
| `sunshine_calu3` | 197 | 6 | 171 | 7 | 3 | −3 |
| `h1_2025` *(in 222's union)* | 255 | 20 | 254 | 7 | **3** | −17 |
| **`kolf2`** | 273 | **273** | 273 | 0 | **0** | **−273** |
| **`feng_ipsc`** | 182 | **182** | 182 | 0 | **0** | **−182** |

**Four things fall out of this and every one of them is actionable.**

1. **`kolf2` and `feng_ipsc` go to EXACTLY ZERO.** They are panel-*fitted* libraries — 273 of 273
   and 182 of 182 targets on today's panel — and the October panel is disjoint by construction.
   They are the #3 and #5 sources by today's overlap and **#28 and #29 in October.** Anything that
   rests on them is worth nothing in the only submission that counts. (Both are already NO-GO on
   price, so nothing we plan to ship is affected — but the principle is now measured.)
2. **DLD-1 goes UP: 67 → 73.** It is a 4,715-target library, not a panel-fitted one. **The one
   source that clears the gate today is also the one that survives October**, which is a
   coincidence worth banking rather than assuming.
3. **`igvf_cardiomyocyte` collapses 87 → 36.** Its 87/300 is substantially an artefact of being a
   **TF** library against a **2.9×-TF-enriched** panel (166 §1). Combined with §8.2's 19% cosine
   retention, its case is now weak on three independent axes. **It stays in the stack for its
   +0.00039 and should not be relied on past October.**
4. **🔴 Four sources this project has written off as "0/300, training-only" become 30–33/300 in
   October** — `nadig_hepg2`, `nadig_jurkat`, `rpe1`, `k562_essential`. **This is exactly 166's
   own structural argument running in reverse**: today's panel *excludes common-essential genes by
   construction*, so every essential-library screen reads 0; the October panel is drawn from donor
   coverage **without that exclusion** (DEF-RUNBOOK §5: `panel ∩ essential list` **0 → 30–34**).
   Those four screens are **already built and on disk**, and each would carry about as many October
   rows as `igvf_cardiomyocyte` and half as many as DLD-1. **They cannot be priced today** — they
   carry 0 rows of today's panel, so PLOLO cannot see them — **and that is precisely why this needs
   to be recorded now rather than discovered on 2026-10-22.**

### 8.4 INBOX-BRANCH-A — baseline, and the units rule *(this file's §8.4, not AGENTS.md's)*

- **Baseline was already `h1union`.** Every number in this file is marginal over
  `base222 = cd4@0.34 + h1pair@0.28`, which is the shipped arm, and §1's KA-1 gate verifies it
  against 218's `compose218` row at Δ = 0 exactly. No re-work was needed.
- **Saturation was already the prioritised question** (Phase 4 got 1,188 of the 3,030 arms).
- **The units rule is obeyed**: §7.1's comparison is now given across the full measured factor
  range (1.000–1.384) rather than converted with one scalar, and every candidate is reported in
  **PLOLO macro cosine (`cos`, per-target unweighted) on the three named cancer folds** against
  **§8.4's +2.5–2.9% shape-only gate**. The `price` column is the project's standard `k = 1` rate
  card and is labelled uncalibrated throughout.
- **On "a source that is 70% distributed is still shippable":** DLD-1 is **54.5%** distributed on
  cosine and **25.2%** on `pds` (§8.2). That is a real discount against `h1union`'s 72.7%/14.6%
  mix, and it is the reason §7.2's rank 1 is recommended on the **stack's** constant-across-size
  cosine retention rather than on any single source's headline.

---

## 9. Files

| path | what |
|---|---|
| `harness225.py` | the source registry + `Fold`; imports 162/218/223's instrument, rewrites none of it |
| `phase0_225.py` → `work/phase0_225.json` | the 20 known-answer gates |
| `census225.py` → `work/census225.json` | 83 assets, panel overlap computed from each file |
| `build225.py` → `work/tables/`, `work/build225.json` | **19 donor tables built here** (1 IGVF, 8 Jiang, 8 scPerturb, 2 sign-flipped twins), 437 MB |
| `sweep225.py` → `work/sweep225{,b,c}.json` | 1,134 + 159 + 549 priced rows with all falsifier variants |
| `stack225.py` → `work/stack225{,d}.json` | 648 + 540 stack arms |
| `analyse225.py` → `work/analyse225.json` | §3–§5 tables · `work/stackfinal225.{txt,json}` · `work/lofo225.json` |
| `conc225.py` → `work/conc225.json` | **INBOX-CONCENTRATION**: leave-top-k-out on both channels, 5 arms × 3 folds |
| `october225.py` → `work/october225.json` | **INBOX-OCTOBER**: breadth + projected coverage on 199's two synthetic draws |
| `work/PROVENANCE.json` | 13 fetched files: URL, accession, bytes, sha256, source, retrieval date; and the 126-sealed exclusion |


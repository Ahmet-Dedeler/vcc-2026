# 298 — RESULTS: the four reviewer-facing checks, on saved outputs only

Run 2026-09-06 17:14–17:29Z. No build, no `vcc`, no `cell-eval2`, nothing written outside
`experiments/298-reviewer-checks/`. One raw-expression pass, the one the mandate authorises: a
streaming group-sum accumulation over the public K562-GWPS bank for item 3. **No R3 reserve file
was opened** — item 1 re-analyses 295's saved per-gene table only. Peak RSS 4.89 GiB (item 2),
3.39 GiB (item 3), both under the 6 GB cap. **No claim is amended here.**

**HEADLINE — C5 WITHOUT ADAMSON: C5a +0.2675 [+0.1045, +0.4294] file-block (HOLDS) / C5b +0.1945
[+0.0059, +0.2997] file-block (holds by the letter, FAILS on the gene bootstrap [−0.0261, +0.3797]);
HALLMARK-ONLY − RANDOM-MATCHED: +0.0590 [+0.0365, +0.2065]; STUDY vs CONSTRUCT in raw cosine:
+0.1680 vs +0.0106 (same 112 perturbations, full depth; difference +0.1574 [+0.1081, +0.2055]);
FAILING SPACERS WITH A SEED OFF-TARGET: not testable on disk (13 of 34 arms fail; neither the 2025
library's spacer sequences nor any promoter sequence is on this machine).**

Artifacts: `work/jackknife298.csv` · `work/jackknife298.json` · `work/ablation298.json` ·
`work/t8_298.json` · `work/t8_298.npz` · `work/seed298.csv` · `work/seed298.json` ·
`work/seed298_plan.md`. Scripts: `work/jack298.py`, `work/ablation298.py`, `work/t8_298.py`,
`work/seed298.py`. Logs: `work/*.log`, `work/*.nohup`.

---

## Item 1 — the single-file jackknife of claim C5

Source: `experiments/295-reserve-r3-claim/work/reserve295.csv` (174 primary genes, 8 files). The
statistics (`spearman`, `boot_rho` file-block, `boot_rho_rows` gene) are copied verbatim from
`stats295.py:36-64`; seed 20260906, 2,000 resamples, as 295.

**KA row.** C5a ρ(`C_held`, `A_g`) = **+0.501615** vs 295's reported +0.5016 (Δ 1.5e-05);
C5b ρ(`gen_pc20`, `A_g`) = **+0.370941** vs +0.3709 (Δ 4.1e-05). The residual is the reported
4-d.p. rounding. Pooled file-block intervals reproduce: C5a [+0.1521, +0.5656], C5b [+0.0627,
+0.4149] (295 reported [+0.153, +0.566] and [+0.064, +0.415]; the last digits differ because this
bootstrap draws from a fresh RNG stream, not 295's continued one). PASS.

**Composition — the concern is real.**

| file | n primary | share | line | modality | median `A_g` | median `C_held` |
|---|---:|---:|---|---|---:|---:|
| AdamsonWeissman2016_GSM2406681_10X010 | 78 | 44.8 % | K562 | CRISPRi | **+0.1316** | +0.1397 |
| TianKampmann2019_day7neuron | 23 | 13.2 % | iPSC-neuron | CRISPRi | +0.0422 | +0.0458 |
| DatlingerBock2021 | 19 | 10.9 % | Jurkat | Cas9-KO | +0.0099 | +0.0158 |
| PapalexiSatija2021_eccite_RNA | 19 | 10.9 % | THP-1 | Cas9-KO | +0.0224 | +0.0532 |
| ShifrutMarson2018 | 17 | 9.8 % | primary CD8 T | Cas9-KO | −0.0035 | +0.0112 |
| DatlingerBock2017 | 8 | 4.6 % | Jurkat | Cas9-KO | −0.0070 | +0.0091 |
| DixitRegev2016 | 8 | 4.6 % | K562 | CRISPRi | +0.0452 | +0.0385 |
| PapalexiSatija2021_eccite_arrayed_RNA | 2 | 1.1 % | THP-1 | Cas9-KO | +0.0061 | +0.0154 |

Adamson carries 78/174 genes and sits at 3–20× the other files' median agreement `A_g`.

**Leave-one-file-out** (`work/jackknife298.csv`; gene = plain gene bootstrap, file = file-block
bootstrap over the remaining files):

| dropped | n | C5a ρ | C5a gene 95 % | C5a file 95 % | C5b ρ | C5b gene 95 % | C5b file 95 % |
|---|---:|---:|---|---|---:|---|---|
| — (all 8) | 174 | +0.5016 | [+0.3707, +0.6129] | [+0.1521, +0.5656] | +0.3709 | [+0.2234, +0.4998] | [+0.0559, +0.4150] |
| **Adamson** | 96 | **+0.2675** | [+0.0611, +0.4396] | [+0.1045, +0.4294] | **+0.1945** | **[−0.0261, +0.3797]** | [+0.0059, +0.2997] |
| DatlingerBock2017 | 166 | +0.5015 | [+0.3659, +0.6262] | [+0.1578, +0.5650] | +0.3674 | [+0.2203, +0.4955] | [+0.0194, +0.4121] |
| DatlingerBock2021 | 155 | +0.5175 | [+0.3741, +0.6365] | [+0.1488, +0.5726] | +0.3828 | [+0.2236, +0.5134] | [+0.0619, +0.4219] |
| DixitRegev2016 | 166 | +0.4957 | [+0.3633, +0.6125] | [+0.1125, +0.5631] | +0.3933 | [+0.2505, +0.5168] | [+0.1236, +0.4233] |
| Papalexi eccite | 155 | +0.5176 | [+0.3700, +0.6358] | [+0.1079, +0.5752] | +0.3803 | [+0.2258, +0.5183] | [+0.0116, +0.4210] |
| Papalexi arrayed | 172 | +0.5084 | [+0.3697, +0.6211] | [+0.1808, +0.5702] | +0.3774 | [+0.2365, +0.5039] | [+0.1031, +0.4161] |
| Shifrut | 157 | +0.5060 | [+0.3722, +0.6248] | [+0.1422, +0.5636] | +0.3595 | [+0.2115, +0.4947] | [+0.0344, +0.3958] |
| Tian | 151 | +0.5396 | [+0.4042, +0.6573] | [+0.1484, +0.5752] | +0.3804 | [+0.2360, +0.5205] | **[−0.0048, +0.4242]** |
| **Adamson ALONE** | 78 | +0.5183 | [+0.3367, +0.6614] | — | +0.3363 | [+0.1413, +0.5007] | — |

**Verdict.** Adamson is the only file whose removal moves either statistic materially: every other
leave-one-out leaves C5a in +0.496…+0.540 and C5b in +0.360…+0.393. Dropping Adamson roughly halves
both — C5a +0.5016 → +0.2675, C5b +0.3709 → +0.1945.

- **C5a HOLDS without Adamson.** 295's frozen rule (ρ ≥ 0.15 and the lower bound > 0) is met on both
  bootstraps: +0.2675 with gene [+0.0611, +0.4396] and file-block [+0.1045, +0.4294]. Adamson alone
  reads +0.5183 — i.e. the claim is *stronger* in Adamson than out of it, but it is not carried by it.
- **C5b is MARGINAL without Adamson.** 295's rule (file-block lower bound > 0) is met by +0.0059 —
  it holds by the letter and by nothing else. The gene bootstrap on the same 96 genes includes zero
  ([−0.0261, +0.3797]), and dropping Tian instead already puts the pooled file-block bound below
  zero ([−0.0048, +0.4242]). **C5b should be reported as resting on the pooled 8-file set, with the
  file-block interval, and not as robust to leaving out either of the two CRISPRi files that carry
  it.** C5a should not.

*Caveat on the intervals.* With Adamson gone, the file-block bootstrap is narrower than the gene
bootstrap (7 blocks of 2–23 genes; the remaining files are homogeneous in `A_g`, so between-file
resampling moves the statistic less than within-file resampling does). Read the wider of the two.

---

## Item 2 — the FIT of hallmark identity alone, vs random matched gene sets

297 row 1c reads 290's control C3 as an ablation: replacing the nine hallmark projections inside
the FULL 18-column model moves +0.3174 → +0.2943, so "random matched sets recover 93 %". That is
the DROP, not the FIT, and it is not the number Wang et al. 2026's Hallmark contrast speaks to.
Refitted here against the same target (`remainder_half2`), the same component folds (2,848 genes,
437 components, fold sizes [2114, 201, 182, 162, 189]) and the same component block bootstrap
(2,000), with `lib290` loaded by path so the ridge, folds and bootstrap are literally 290's objects.

**KA row.** The full 18-column model refits to **+0.31738512466156** — 290's primary cell
+0.317385 to 1e-9 (asserted; the run aborts otherwise). The annotation-21 ridge refits to
−0.056642, 290's value. The matched random-set sizes reproduce ctrl290's exactly
(E2F 104, G2M 102, MYC 96, P53 44, UPR 43, APOPTOSIS 49, MTORC1 101, OXPHOS 79, IFNγ 38;
proliferation 228) and the draw sequence uses ctrl290's seed (`SEED + 91`).

| model | columns | OOF Spearman vs `remainder_half2` | component 95 % |
|---|---:|---:|---|
| **(i) HALLMARK-ONLY** (9 `pe_<HALLMARK>` + `pe_prolif_abs`) | 10 | **+0.1576** | [+0.1385, +0.2996] |
| **(ii) RANDOM-MATCHED** (same 10 columns, random gene sets, 20 draws) | 10 | **+0.0986** ± 0.0191 (range +0.0484 … +0.1270; **0/20 reach (i)**) | — |
| (iii) SPECTRAL-ONLY (`gen_pc5`, `gen_pc20`, `gen_cos`, `gen_abscos`, `breadth`, `logmag`) | 6 | +0.2424 | [+0.2156, +0.3236] |
| full frozen model (290's primary) | 18 | +0.3174 | [+0.2897, +0.3914] |
| annotation-21 ridge (290's comparator) | 21 | −0.0566 | [−0.1380, −0.0205] |

**The number 297 asked for**

| contrast | value | component 95 % | P(≤ 0) |
|---|---:|---|---:|
| **(i) − (ii), against the mean of the 20 draws** | **+0.0590** | **[+0.0365, +0.2065]** | 0.000 |
| (i) − (ii), against the BEST of the 20 draws | +0.0306 | [+0.0056, +0.1345] | 0.013 |
| (iii) − (i), spectral minus hallmark | +0.0847 | [−0.0301, +0.1024] | — |
| full − (iii), what the other 12 columns add over spectral | +0.0750 | [+0.0206, +0.1131] | — |

**Reading.** Hallmark *identity* carries real, non-generic signal on its own: fitted alone it reaches
+0.1576, it beats matched random gene sets by +0.0590 [+0.0365, +0.2065], and not one of the 20
random draws reaches it. **297's contradiction with Wang et al. 2026's Hallmark contrast dissolves
— it was an artefact of reading an ablation-inside-a-full-model as a fit.** What 290's C3 actually
established survives too, but it is a *redundancy* statement, not an absence: identity's fit is
largely spanned by response geometry — spectral columns alone reach +0.2424, above hallmark-only
(+0.0847, interval crossing zero, so the ordering is not resolved), and adding the hallmark block on
top of the full model buys only the +0.023 that 290 reported. Both sentences are true at once, and
F13's "the signal is response GEOMETRY, not programme identity" should be narrowed to "identity adds
little *once geometry is in the model*", which is not the same as "identity does not predict".

The ridge's chosen α tells the same story: hallmark-only picks α = 31.6–316 across folds (heavy
shrinkage, a weak diffuse signal), spectral-only picks α = 0.003–0.32 (a sharp one).

---

## Item 3 (T8) — the study-vs-construct contrast in BOTH statistics

F14 quotes a *reliability ratio* (across-study median cosine ÷ within-study median cosine = 0.79);
Replogle et al. 2022's same-line cross-experiment number is a *raw correlation* (≈ 0.50). Both are
reported here on the same K562 pair, plus the same two raw statistics for the construct pair.

**KA rows.** 285's K562 study row at n = 25 re-read: within **+0.0605**, across **+0.0479**, ratio
**0.7912** (want +0.0605 / +0.0479 / 0.791) — PASS. 285's construct row at n = 25: ρ_guide
**+0.0896** on **748** perturbations — PASS. The run aborts on either failure.

**(a) The reliability ratio (285's statistic, cell-matched sides).**

| pair | n | within (one source's cell halves) | across | ratio |
|---|---:|---:|---:|---:|
| STUDY, n = 25/side | 1,765 | +0.0605 [+0.0561, +0.0655] | +0.0479 [+0.0440, +0.0513] | **0.791** |
| STUDY, n = 50/side | 1,163 | +0.0942 | +0.0714 | 0.758 |
| CONSTRUCT, n = 25/side | 748 | +0.00632 | +0.00057 | **0.090** [−0.038, +0.232] |
| CONSTRUCT, n = 50/side | 653 | +0.01025 | +0.00063 | 0.062 [−0.040, +0.177] |

**(b) The raw statistics at FULL depth.** Study pair: 2,053 shared perturbations, common axis 8,202
genes (6,357 scored after removing the shared targets' own genes), median 170 (GWPS) / 121
(essential) cells per perturbation. Construct pair: 775 two-construct perturbations of K562-GWPS,
7,690 scored genes, median 212 (larger construct) / 111 (smaller) cells. `raw_cos` is the cosine of
the un-centred effect vector (`prof(S) − prof(ctrl)`), `raw_pearson` the Pearson of the same two
vectors across genes, `resid_cos` 285's panel-centred residual cosine.

| pair | n | raw cosine | raw Pearson | resid cosine | share ≤ 0 |
|---|---:|---:|---:|---:|---:|
| STUDY, all shared | 2,053 | +0.2150 [+0.1997, +0.2324] | +0.2064 [+0.1926, +0.2208] | +0.2384 | 0.046 |
| STUDY, 285's n = 50 cohort | 1,163 | +0.2284 [+0.2041, +0.2455] | +0.2172 [+0.1986, +0.2404] | +0.2645 | 0.044 |
| CONSTRUCT, all two-construct | 775 | +0.0057 [+0.0038, +0.0078] | +0.0045 [+0.0024, +0.0061] | +0.0033 | 0.395 |

**(c) The same 112 perturbations (study-shared ∩ two-construct) — where the inversion has to be read.**

| statistic | STUDY | CONSTRUCT | STUDY − CONSTRUCT (paired target bootstrap) | P(≤ 0) |
|---|---:|---:|---|---:|
| raw cosine | **+0.1680** [+0.1195, +0.2149] | **+0.0106** [+0.0035, +0.0176] | **+0.1574 [+0.1081, +0.2055]** | 0.0000 |
| raw Pearson | +0.1576 [+0.1099, +0.1987] | +0.0100 [+0.0009, +0.0178] | +0.1476 [+0.1002, +0.1899] | 0.0000 |
| resid cosine | +0.2349 [+0.2050, +0.2820] | +0.0050 [−0.0037, +0.0126] | +0.2298 [+0.1999, +0.2781] | 0.0000 |

**(d) By response-magnitude tertile** (geometric mean of the pair's residual norms) — the reviewer's
next question, since Replogle's ≈ 0.50 is quoted on strongly-responding perturbations while the rows
above use every shared perturbation and every scored gene:

| stratum | n | median magnitude | raw cosine | raw Pearson |
|---|---:|---:|---:|---:|
| STUDY T1 (weakest) | 685 | 5.82 | +0.0889 | +0.0852 |
| STUDY T2 | 684 | 8.17 | +0.2369 | +0.2288 |
| **STUDY T3 (strongest)** | 684 | 12.30 | **+0.4113** [+0.3833, +0.4310] | **+0.3798** [+0.3493, +0.4054] |
| CONSTRUCT T1 | 258 | 5.99 | +0.0064 | +0.0054 |
| CONSTRUCT T2 | 259 | 7.17 | +0.0052 | +0.0045 |
| **CONSTRUCT T3 (strongest)** | 258 | 9.05 | **+0.0057** [+0.0036, +0.0091] | **+0.0037** [+0.0013, +0.0070] |

**Reading.**
1. **The inversion is not an artefact of the ratio statistic.** In one common unit, at full depth,
   on the same 112 perturbations: a different study in the same line agrees at raw cosine +0.168,
   a second construct of the same gene at +0.011 — a 16× gap, +0.157 [+0.108, +0.206], P(≤ 0) = 0.
   F14's sentence survives the translation.
2. **Our raw study number (+0.21 pooled) and Replogle's ≈ 0.50 are the same statistic on different
   selections, and the gap is magnitude.** On the strongest-responding third the raw Pearson reads
   +0.380 and the raw cosine +0.411, next to +0.085 in the weakest third. The remaining distance to
   0.50 is what an expressed-gene / responding-gene filter and the deeper cells of that paper's own
   construction buy. The two numbers are reconcilable and should be quoted with their cohort.
3. **The construct term does not close with magnitude.** Study agreement rises 4.6× from the weakest
   to the strongest tertile (+0.089 → +0.411); construct agreement is flat (+0.0064 → +0.0057).
   Whatever separates two constructs of one gene is not a small-response problem — which is exactly
   what 296 found by a different route (failing constructs are reproducible but share no vector
   signature) and what 293 found (the term is knockdown failure).
4. 39.5 % of two-construct perturbations have a raw cosine ≤ 0 between their constructs, against
   4.6 % of shared perturbations between studies.

---

## Item 4 (T7) — the seed scan of failing spacers against promoters

**KA row.** 293's `native` n = 200 ρ_guide re-read: **+0.0745** [−0.0429, +0.3208] on 16 targets
(want +0.075 / 16) — PASS.

**Verdict: NOT TESTABLE ON DISK; item 4 stops at the plan, as the mandate directs.** Both halves of
the scan are missing and only one of them is downloadable under this mandate:

- **Spacer sequences — missing, not downloadable here.** The 2025 files' `obs/guide_id` holds
  protospacer *names*: `GENE_P1P2_A|GENE_P1P2_B`. **0 of 410 categories across the three splits are
  ACGT strings.** `obs` = `[_index, batch, guide_id, target_gene]`; `uns` is empty; `var` carries
  only gene ids and QC columns. The whole 2025 release on disk is three h5ads + `gene_names.csv` +
  three `pert_counts_*.csv` — no library table. The three sgRNA/library files that do exist on disk
  (`cd4_cell2026/.../sgrna_library_metadata.suppl_table.csv`, `gse249595/..._guides_reference.fa.gz`,
  `myllia-a549/GSE337804_sgRNA_lib.csv.gz`) all belong to other studies.
- **Promoter sequence — missing, would be downloadable.**
  `data/external/broad16line/metadata/gene_locations_biomart.csv` (86,402 rows) gives Ensembl id,
  chromosome, gene start/end and karyotype band: it fixes the TSS ± 500 bp *window* but carries no
  bases, and it lacks the strand column the TSS needs. No FASTA/2bit/genome file exists anywhere
  under `data/`. ~20,000 × 1,001 bp of promoter sequence is ≈ 6 MB gzipped, inside the 50 MB cap —
  but fetching it alone does not make the scan runnable.

**What was produced instead.** `work/seed298.csv` — the exact input list a future run needs: all
**34 arms** of the 17 two-construct targets, with construct id, cells, target-gene knockdown lfc and
the failing flag. `work/seed298_plan.md` — the full match logic (arms → spacers → TSS ± 500 windows
→ 9–12 bp PAM-proximal seed match with an NGG PAM, either orientation → the candidate gene's own
column in the failing construct's residual response → the null on the 21 non-failing arms), with
each step's data requirement marked present or missing.

**The failing population, which is measurable from saved tables:** 13 of 34 arms fail at
|lfc| < 1.0 — ARPC2_P2, CAST_P1, CAST_P2, DPF2_P2, GRB10_P2, HIRA_P1, IGF1R_P2,
KDM2B_ENST00000377071.4, KLF10_P1, MED15_P2, PCCB_P2, SMARCB1_P2, ZNF581_P2. One target (CAST) has
*both* constructs failing; five (DAXX, PRCP, SOX4, STAG2, TMSB4X) have neither. **13, not 34, is the
denominator any future "n with a candidate seed off-target" count must be read against**, and with
a published off-target rate of 0.8–4.7 % for a validated library (297), a 13-arm scan can produce an
existence proof but cannot estimate a rate.

---

## Provenance

| item | inputs (read-only) | statistics copied from |
|---|---|---|
| 1 | `295/work/reserve295.csv` | `295/work/stats295.py:36-64` (= `lib291.py:73-78, 198-207`) |
| 2 | `290/work/feat290.csv`, `pred290.csv`, `284/work/rem284.npz`, `splits284.npz`, `281/work/feat281.csv`, the 14 half-1 residual banks | `290/lib290.py` loaded by path (ridge, `component_folds`, `block_boot`); random-set construction `290/ctrl290.py:36-99`, same seed |
| 3 | `285/work/k562_{gwps,essential}_meta.npz`, `k562_{gwps,essential}_study_full.npy`, `k562_285.json`, `k562_guide.json`; one streaming pass over `data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad` (group sums only, cells never held) | `285/work/lib285.py` loaded by path (`prof`/`delta`/`resid`/`row_cos`/`boot_median`/`boot_median_diff`); cohort construction `k562_analyse.py:234-239`, `k562_guide.py:47-49` |
| 4 | `293/work/construct293.csv`, `ceiling293.json`; `vcc2025/*.h5ad` obs metadata only; `broad16line/metadata/gene_locations_biomart.csv` | — |

`construct293.csv` sha256 `4cca9aa6…` recorded in `work/seed298.json`. Nothing outside this experiment
directory was written or modified.

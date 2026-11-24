# 218 — Independence does **not** predict contribution. Alignment does. And the binding gate on donor acquisition was never independence — it is production's **median sizing rule**, which silently zeroes every source below 50% panel coverage.

**Date:** 2026-09-04 · deterministic (seeds `20260904`, `20726` for permutations) · **no `vcc`
command, no `cell-eval2`, no submission built, nothing under `submissions/` or `data/donors/`
touched, no production file modified, nothing posted anywhere, no account created.**
Everything written is inside `experiments/218-independent-donors/` except **98 MB** of new
downloads (`data/external/hspc_crispri500` 20 MB, `data/external/broad16line` 78 MB) — against
an 80 GB budget and a 15 GB reporting threshold, so no inventory report was owed.
Peak RSS **7.65 GiB** (`/usr/bin/time -l`, the 476-arm sweep). 8 scripts, **1,050 scored
configurations** across `resid218` (476), `compose218` (496) and `rowsplit218` (32).

---

## 0. The answer in ten lines

1. **The mandate's hypothesis is REFUTED by its own designed test.** Within each PLOLO fold,
   the Spearman correlation between a source's contribution and its **independence** from the
   production mixture is **−0.057** (per fold: +0.20, +0.86, −0.46, −0.82 — the sign flips).
   Against **matched-target transfer alignment** it is **+0.696** (+0.30, +0.89, +0.77, +0.83 —
   positive on 4 of 4).
2. **The decisive row needs no regression.** CD4, KOLF2 and Feng-iPSC score **identically 4/4**
   on every categorical independence axis (platform, lineage, transformation state, lab). Their
   contributions differ **6.1×**: **+11.56% / +1.90% / +5.08%** relative aggregate cosine at the
   matched shipped dose. The ordering is exactly their alignment ordering, and the *most*
   orthogonal source (Feng, `orth_frac` 0.975) is not the best.
3. **The right null is not zero, it is +4.61%** — 099's self-residual, i.e. whitening relaxation
   available for free from our own data. **Only CD4 clears it by a real margin (+6.96 pp).**
   Feng clears by +0.48 pp; **KOLF2 is 2.71 pp *below* it**, i.e. a maximally-independent iPSC
   line is worse than adding back our own effect.
4. **166's search was not wrong about the universe — but its bar was in the wrong place.**
   Under the new criterion I re-swept HuggingFace, GEO (167 human CRISPR scRNA series since
   2025-06), Europe PMC, BioStudies/ArrayExpress, CZ CELLxGENE (388 collections), the Arc
   atlas, scPerturb (64 h5ads re-censused with **166's own `census.py`, byte-identical copy**),
   and **figshare — which 166 never swept**. Three genuinely new items, all measured, all
   rejected: the **KOLF2.1J genome-scale raw atlas** (250 GB, CC-BY-4.0; its 46.7 GB
   panel-sized product carries **55/300** — measured remotely without downloading it), the
   **Broad 16-line** set (now posted; **0/300**, measured), and an **HSPC CRISPRi500** screen
   that is **mouse**. X-Atlas/Pisces is still an empty repository.
5. **The real gate is production's sizing rule, not availability.** `solve_cd4_fraction` sizes
   the residual by the **median row norm over every panel row**, so a source covering ≤50% of
   the panel has median 0 and the solve returns **exactly 0**. Two sources we have held for
   weeks are therefore **identically absent** from every build: Pan-2026 H1 (64/300) and
   `h1_2025` (20/300). This reproduces 200's silent zero and 207's h1-test finding on a third
   surface.
6. **Sized on their carrying rows instead, they are the largest direction arm this instrument
   has measured.** The union of the two (77/300 = 25.7% of the panel, only 7 targets shared)
   on top of the board incumbent `cd4f034`, at `g = 0.28`: **macro cosine +7.680% and `pds`
   +0.004887 across the three cancer-line folds**, and positive on cosine *and* `pds` on **all
   four folds at all three doses — 12 of 12**.
7. **Four falsifiers, all pass.** Target-shuffled control **strictly negative at every dose**
   (−0.60 / −1.42 / −3.16%); macro cosine moves *with* the energy-weighted one, so it is not a
   row-reweighting artefact; the **row-matched CD4 control** — same rows, same sizing, same
   dose, CD4's information instead — reproduces only **1/3.9** of the macro gain on the cancer
   folds and **costs** `pds` where the real arm gains it; and the rows the source does not carry
   move by **+0.00%** at every dose on every fold.
8. **Priced under §8.4's corrected rule, all six members, never carrying four at zero:**
   conservative (macro cosine, the three **cancer-line** folds only, uncalibrated `k = 1`)
   **+0.00919**, and **+0.01054 marginal** over the free row-matched-CD4 alternative — **5.5×**
   rule 13's shape-only gate. Leave-one-fold-out, dose chosen on the other three folds:
   **+0.02830**, positive on 4/4.
9. **`pds` goes UP, which no CD4 dose arm has ever done.** `cd4f034` paid −0.001345 scaled on
   `pds` for its cosine; this arm gains `pds` on every fold at every dose. That is the signature
   of adding a genuinely 4th independent estimate on the rows it covers, and it is consistent
   with 159's `1/k` law: at `g = 0.17` the carrying rows gain **+8.8% to +11.7%** cosine against
   the law's **+13.7%** prediction for `k` 3 → 4.
10. **Handed back, not built.** This needs a builder change (a coverage-aware sizing rule) that
    lives in `scripts/make_submission_pm_cd4.py` — outside my territory and outside this
    mandate. **CD4's own dose and sizing are untouched throughout** (216/217's territory).

---

## 1. Known-answer rows — nine, all pass

Rule 10. Every instrument here carries a configuration whose answer was already known, and
each script **aborts** on failure.

| # | gate | expected | measured | verdict |
|---|---|---|---|---|
| KA-1 | `cd4_gn0.1` on all 4 folds | 207's `f0.1` `cos_agg`/`pds` | **identical, Δ = 0.0e+00** | PASS |
| KA-2 | `f0` (residual off) on all 4 folds | 207's `f0` | **identical, Δ = 0.0e+00** | PASS |
| KA-3 | `kolfmask{0.05,0.1,0.2,0.4}` | 207's `kolf*` (its `alt_line_residual`) | **identical, Δ = 0.0e+00** | PASS |
| KA-4 | my multi-bank residual builder on the 3 CD4 banks, `gene_norm` | production `prepare_cd4_mean` | **max abs diff 0.0e+00**, votes identical, 4/4 folds | PASS |
| KA-5 | 159 §4 `k562_gwps × orion_hct116` | +0.0321 / +0.0020 / +0.0321 | +0.0323 / +0.0020 / +0.0322 | PASS |
| KA-6 | 159 §4 `k562_gwps × orion_hek293t` | +0.0289 / +0.0016 / +0.0286 | +0.0302 / +0.0017 / +0.0298 | PASS |
| KA-7 | 159 §4 `orion_hct116 × orion_hek293t` | +0.0306 / +0.0016 / +0.0300 | +0.0297 / +0.0019 / +0.0291 | PASS |
| KA-8 | `compose218` at `g = 0` | `resid218`'s `cd4_gn0.34` | identical to 1e−9 | PASS |
| KA-9 | non-vacuity, every arm | median residual row norm > 0 where claimed | see §7.1 | PASS |

KA-5/6/7 run **166's `independence.py` verbatim** through a byte-identical copy
(`shasum` matched before the run) inside 166's own statistical tolerance.

> ⚠️ **A trap that cost two runs and is worth recording.** 062's `enforce_memory_limit` reads
> `ru_maxrss`, a process **high-water mark** that `gc.collect()` cannot lower. 166's loader
> materialises three full donor tables (~6 GB), so the KA gate poisoned every later stage in the
> same process and the guard aborted at 14.5 GiB. **The gate now runs in its own subprocess.**
> Freeing the memory is not enough; the peak has to never happen in that process.

---

## 2. Phase 0a — independence, measured in the frame the residual is actually added in

`independence218.py` → `work/independence218.json` (36 rows). Per PLOLO fold, on the fold's
scored columns, for the residual matrix each source contributes:

| metric | meaning |
|---|---|
| `cos_pre` / `excess_pre` | per-target cosine to the base mixture — **redundancy** |
| `sub_10` | share of the source's energy inside the consensus's top-10 right-singular subspace (209's shared-subspace energy) |
| `orth_frac` | share of energy **outside** the consensus row space entirely |
| `cos_truth` / `excess_truth` | per-target cosine to the **held line's truth** — transfer alignment |

Four-fold means, external sources first:

| source | transformation | `orth_frac` | `sub_10` | `excess_pre` | `cos_truth` |
|---|---|---:|---:|---:|---:|
| **cd4** | NON-transformed **primary** human cells | 0.9389 | 0.0317 | +0.0277 | **+0.01812** |
| **feng_ipsc** | NON-transformed iPSC (24–34 donors) | 0.9753 | 0.0088 | +0.0130 | +0.01039 |
| **kolf2** | NON-transformed iPSC | 0.9681 | 0.0168 | +0.0072 | +0.00563 |
| gse295214_h1 | NON-transformed hESC | 0.9308 | 0.0339 | +0.0607 | +0.05466 † |
| telohaec | hTERT-immortalised primary | 0.9793 | 0.0078 | +0.0040 | +0.00319 |
| k562_gwps | production donor | 0.7667 | 0.0394 | +0.5720 | +0.01667 |
| orion_hct116 | production donor | 0.6387 | 0.1074 | +0.6545 | +0.02085 |
| orion_hek293t | production donor | 0.5967 | 0.1740 | +0.6586 | +0.02031 |
| perturbmatch | K562 reprocessing (same cells) | 0.8081 | 0.0527 | +0.5081 | +0.01580 ‡ |
| self_pre (consensus) | — | 0.3161 | 0.2607 | +0.9981 | +0.03005 |

† unmatched — §5 shows this number is mostly a target-set artefact and must not be quoted raw.
‡ `perturbmatch × k562_gwps` reads `cos_truth` **+0.639** on the `k562_gwps` fold — the same
cells reprocessed. That row is a **leak by construction**, it is excluded everywhere, and it is
recorded because it is the instrument demonstrating that it *can* detect a leak (19.1× the CD4
baseline on matched rows).

**The independence ordering is clean and it is the opposite of the useful ordering.** The three
external non-transformed sources are the *most* orthogonal to our mixture (0.94–0.98) and the
production donors the least (0.60–0.77) — exactly as the hypothesis says. It does not help.

---

## 3. Phase 0b — contribution, measured the way CD4 is measured

`resid218.py` → `work/resid218.json` (476 arms). Every source put through the **identical**
construction: an un-whitened post-whitening residual at a solved dose `f`, 162's PLOLO
instrument, production's operator/α/γ/λ/PerturbMatch blend/renderer/read-back all untouched.
Baseline `f0` = no residual. Primary endpoint `cos_agg` (204: instrument resolution ±0.20%
relative), `pds` secondary, `q_w` third.

At the **matched shipped dose `f = 0.10`**, 4-fold means (3 where the source is a consensus
member on one fold and cannot be used on its own):

| source | folds | rel Δ`cos_agg` | Δ`cos_agg` | Δ`pds` | Δ`q_w` | ×CD4 |
|---|---:|---:|---:|---:|---:|---:|
| **cd4** | 4 | **+11.562%** | +0.001995 | **+0.007005** | +0.000450 | 1.000 |
| feng_ipsc | 4 | +5.084% | +0.000631 | +0.002720 | +0.000112 | 0.440 |
| **self_pre (099's control)** | 4 | **+4.606%** | +0.000873 | −0.007788 | +0.000268 | 0.398 |
| orion_hct116 | 3 | +3.928% | +0.000704 | −0.005696 | +0.000218 | 0.340 |
| orion_hek293t | 3 | +3.699% | +0.000768 | −0.006830 | +0.000324 | 0.320 |
| k562_gwps | 3 | +3.298% | +0.000470 | +0.000143 | −0.000326 | 0.285 |
| perturbmatch | 3 | +2.770% | +0.000367 | −0.000406 | +0.000629 | 0.240 |
| kolf2 | 3 | +1.896% | +0.000423 | −0.000206 | +0.000576 | 0.164 |
| gse295214_h1 | 4 | **+0.000%** | 0.000000 | 0.000000 | 0.000000 | **VACUOUS** |
| telohaec | 4 | **+0.000%** | 0.000000 | 0.000000 | 0.000000 | **VACUOUS** |

CD4's `+11.562%` and the self-residual's `+4.606%` reproduce 207 §5.3's `+11.0%` and `+4.8%`
(207 quoted them at 4 significant figures off the same rows; the difference is rounding).

### 3.1 The right null is the self-residual, and almost nothing clears it

099 warned that adding back **our own** un-whitened effect recovers most of this channel for
free. It does: **+4.606%**. Measured against that floor rather than against zero:

| source | rel Δ`cos_agg` | excess over the free floor |
|---|---:|---:|
| **cd4** | +11.562% | **+6.956 pp** |
| feng_ipsc | +5.084% | +0.478 pp |
| orion_hct116 | +3.928% | −0.678 pp |
| orion_hek293t | +3.699% | −0.907 pp |
| k562_gwps | +3.298% | −1.308 pp |
| perturbmatch | +2.770% | −1.836 pp |
| **kolf2** | +1.896% | **−2.710 pp** |

**A maximally-independent, non-transformed, different-platform, different-lab iPSC line is
worse than using nothing new at all.** That is the cleanest single refutation of the mandate's
premise on this page. (It also explains 057's "KOLF2 loses after whitening" mechanistically
rather than empirically, and it holds under `gene_norm` *and* `raw` normalisation — both were
swept, both agree.)

---

## 4. Phase 0c — THE REGRESSION. This is the load-bearing result.

`regress218.py` → `work/regress218.json`. 27 non-vacuous, non-leaking (source, fold) points at
the matched dose, one endpoint, two rival predictor families fitted on the **same** points.

| hypothesis | predictor | n | r | r² | slope > 0 |
|---|---|---:|---:|---:|---:|
| **H1** | categorical novelty score (0–4) | 27 | +0.355 | 0.126 | 95.8% |
| **H1** | `orth_frac` (energy outside the consensus row space) | 27 | +0.176 | 0.031 | 84.2% |
| **H1** | `1 − sub_10` (outside the consensus top-10 PCs) | 27 | +0.146 | 0.021 | 82.3% |
| **H1** | `−excess_pre` (target-specific redundancy, negated) | 27 | +0.262 | 0.069 | 90.6% |
| **H0** | `cos_truth` (alignment to the recipient) | 27 | +0.244 | 0.060 | 91.6% |
| **H0** | `excess_truth` (target-**specific** transfer) | 27 | +0.246 | 0.060 | 92.0% |
| **H0, confound-controlled** | **matched-target alignment ratio vs CD4** (§5) | 23 | **+0.455** | **0.207** | **100.0%** |
| — | coverage | 27 | +0.113 | 0.013 | 69.2% |

The pooled fits are all weak, because pooling mixes two mechanisms (whitening relaxation for
consensus members, new information for external sources) and three folds of very different
truth quality. **The within-fold rank correlation removes both and it is unambiguous:**

| predictor | k562_gwps | orion_hct116 | orion_hek293t | kolf2 | mean |
|---|---:|---:|---:|---:|---:|
| **matched alignment** | +0.30 | +0.89 | +0.77 | +0.83 | **+0.696** |
| **independence (`orth_frac`)** | +0.20 | +0.86 | −0.46 | −0.82 | **−0.057** |

### 4.1 The decisive row — three sources, identical independence, 6.1× apart

| source | novelty | `orth_frac` | matched alignment | rel Δ`cos_agg` | ×CD4 |
|---|---:|---:|---:|---:|---:|
| **cd4** | 4/4 | 0.9389 | 1.000 (reference) | **+11.562%** | 1.00 |
| **feng_ipsc** | 4/4 | 0.9753 | 0.31–1.06 | +5.084% | 0.44 |
| **kolf2** | 4/4 | 0.9681 | 0.23–0.78 | +1.896% | 0.16 |

Every categorical independence axis the mandate names — third platform, third chemistry,
non-transformed rather than immortalised, different lab, different lineage — is **identical**
across these three. Their contributions differ by a factor of six, monotone in alignment. The
categorical score's r = +0.355 in the pooled fit is carried entirely by the *contrast between
external sources and consensus members*, and that contrast is confounded with alignment; inside
the external set, independence carries no information at all.

**So: independence is necessary bookkeeping, not the mechanism.** What CD4 has that KOLF2 does
not is 22M cells across 4 human donors and 3 stimulation states, which makes it an unusually
good estimator of the shared component — 207 §7.2 reached the same conclusion from a different
direction and this is the controlled measurement of it.

---

## 5. Phase 0d — the confound control, without which this document would have been wrong

`matched218.py` → `work/matched218.json`. 166 §4.2 measured that between the **same** two cell
lines the choice of target set moves the transfer cosine **×3.98** while the choice of line
moves it **×1.36**. A source carrying 61 of 290 targets therefore cannot have its `cos_truth`
compared to one carrying 282 of them.

**Raw, this experiment appeared to find that Pan's H1 hESC screen aligns 1.8–10.1× better with
our recipients than CD4 does.** On matched rows it collapses to 1.0–1.2× on the cancer lines
(paired-target bootstrap 95% CIs are in `work/matched218.json`; `n` = 58–63 per fold, 20 for
`h1_2025`):

| source | k562_gwps | orion_hct116 | orion_hek293t | kolf2 |
|---|---:|---:|---:|---:|
| gse295214_h1 **raw** (its rows vs CD4 on all of CD4's rows) | 1.77× | 2.13× | 2.31× | 10.09× |
| gse295214_h1 **matched** (CD4 on the same rows) | **1.018×** | **1.110×** | **1.221×** | **3.991×** |
| h1_2025 **raw** | 2.19× | 2.58× | 3.11× | 8.86× |
| h1_2025 **matched** | 1.806× | 1.156× | 1.842× | 3.061× |
| feng_ipsc matched | 0.315× | 0.310× | 0.465× | 1.055× |
| kolf2 matched | 0.277× | 0.226× | 0.247× | — |
| telohaec matched | 0.095× | 0.140× | 0.171× | 0.260× |
| **perturbmatch matched (LEAK)** | **19.114×** | 0.947× | 1.247× | 0.704× |

On the three **cancer-line** recipients — the class D/E/F will belong to (§2.1) — the H1 sources
are **1.0–1.8× CD4 on matched rows**, not 10×. The apparent advantage was the target set, and
166's control is what removed it. On the one **pluripotent** recipient they are 3–4×, which is a
state-match effect and is reported as such, not generalised.

---

## 6. Phase 1 — the search, under the new criterion

**Criterion:** human, transcriptome-wide perturbation→expression readout, maximally independent
of Replogle/Orion, **and usable panel coverage** — where §7 now shows "usable" means >50% of
the panel under production's median sizing rule, or a coverage-aware rule.

Registries swept: **HuggingFace** (`hf datasets ls`, 5 query families), **GEO** (E-utilities,
167 human CRISPR/CRISPRi/CRISPRa scRNA-seq series released 2025-06-01 → 2026-12-31, ranked by
sample count), **Europe PMC**, **BioStudies/ArrayExpress** (122 hits), **CZ CELLxGENE** (388
public collections, 20 perturbation-adjacent, none a genome-scale CRISPR screen), **Arc Virtual
Cell Atlas** (repo contents + 8 most recent commits — still only `scBaseCount`, `tahoe-100M`,
`virtual-cell-challenge`), **scPerturb** (64 on-disk h5ads re-censused), **figshare** (three
query families — **166 never swept this registry**, and it is where two of the three new items
were), **IGVF** (166's closure re-confirmed via the Stem Cell Reports benchmarking paper behind
it), **Zenodo** (**still WAF-blocked, HTTP 403 on the API with and without a browser UA** — 166's
one unturned stone remains unturned).

### 6.1 Coverage numbers, computed not eyeballed

`census166_copy.py` is **166's `census.py`, byte-identical** (`sha256`
`91291917d1f8374992d2c92d1b6aedc84622a3b5147d15b40b99b656b8e8ebae`, verified against the
original before running), run from this directory so its output lands in `work/census.json`.
64 files, 57 with a readable perturbation column. **Panel overlap outside K562-gwps's 272:
21, 13, 8, 6, 6 — nothing else on disk is within an order of magnitude of the bar.**

`sources218.py` does the same for the processed donor tables, which are not h5ads and so are
invisible to `census.py`:

| source | panel | fold coverage (k562 / hct116 / hek293t / kolf2) | median rule non-vacuous? |
|---|---:|---|---|
| orion_hek293t | 294 | 263 / 289 / 294 / 268 | yes |
| orion_hct116 | 290 | 260 / 290 / 289 / 264 | yes |
| **cd4** | 276 | 251 / 267 / 271 / 253 | yes |
| kolf2 | 273 | 249 / 264 / 268 / 273 | yes |
| k562_gwps | 269 | 269 / 260 / 263 / 249 | yes |
| perturbmatch | 251 | 248 / 242 / 245 / 231 | yes |
| feng_ipsc | 182 | 171 / 177 / 179 / 165 | yes (0.61) |
| **telohaec** | 65 | 58 / 62 / 63 / 59 | **NO — zeroed** |
| **gse295214_h1** | 64 | 61 / 62 / 63 / 58 | **NO — zeroed** |
| **h1_2025** | 20 | 20 / 20 / 20 / 20 | **NO — zeroed** |
| myllia_a549 | 14 | 13 / 12 / 13 / 10 | NO — zeroed |
| frangieh_melanoma | 5 | 5 / 5 / 5 / 4 | NO — zeroed |
| nadig_jurkat, nadig_hepg2, rpe1, k562_essential | 0 | 0 | n/a (166's disjointness) |

### 6.2 New items found, and the number that decided each

| item | organism / system | perturbations | **panel** | licence | size | verdict |
|---|---|---:|---:|---|---:|---|
| **KOLF2.1J genome-scale atlas, raw** — figshare+ `10.25452/figshare.plus.27261219`, article 27261219, published 2026-05-18 | human iPSC (KOLF2.1J), CRISPRi, non-transformed | **11,739** | ≤273 in the pan-genome file | **CC-BY-4.0** | **250.6 GB** (pan-genome 189.4, "strong" 46.7, +4 smaller) | **NO-GO.** The only in-budget product, `KOLF_Strong_Perturbations.h5ad` (46.7 GB), carries **1,656 targets and 55/300** — measured **remotely** over HTTP range reads without downloading a byte of `.X`. The panel-covering product is the 189.4 GB pan-genome file, which is CSC over cells so a panel-restricted extraction still reads ~all of it — over budget. And §3 measures the KOLF2 **line** at 0.16× CD4 and *below* the free relaxation floor, so a better-conditioned rebuild of the weakest source we hold is not where the value is. **166 did not have this article** (it held only the 30 MB webapp z-score slice); recording it closes the row. |
| **Broad 16-line Perturb-seq** — figshare `10.6084/m9.figshare.33273600`, posted 2026-08-18; preprint `10.64898/2026.08.24.746802` | 16 human cancer lines, 8 OncoTree lineages | **100** | **0/300, measured** | CC-BY-4.0 | 8.7 GB | **NO-GO on coverage.** 166 flagged this as "data not yet posted, realistically ~0". It is now posted and the overlap is **exactly 0**, computed from this repo's own experiment 111 target list (90 non-control targets, intersection with the panel empty). 166's estimate is confirmed as a measurement. Still the densest line×target matrix in public data and a real 161-class training asset. |
| **HSPC CRISPRi500** — figshare `10.6084/m9.figshare.31276981` + `31277494`, posted 2026-08-26 | **ex vivo haematopoietic stem/progenitor cultures** — exactly the "primary cell" class the mandate asks for | **520** | **0/300** | CC-BY-4.0 | 11.6 GB (21 MB DE table downloaded) | **NO-GO — wrong organism.** The DE table's response genes are `0610005C13Rik`, `0610009E02Rik`… and **0 of its 20 sampled targets are all-caps human symbols**. It is **mouse**. Cost of finding out: one 21 MB download. |
| **X-Atlas/Pisces** | HepG2, iPSC, Jurkat rest/activated, 16 contexts, 25.6M cells | — | — | CC-BY-NC-SA-4.0 | **0 bytes** | **Still empty**, re-checked 2026-09-04: 3 files, `used_storage: 0`, `last_modified` 2026-03-25. Unchanged in 5½ months. |
| E-MTAB-14004 "genome-wide genetic and environmental perturbations" | **Saccharomyces cerevisiae** | genome-wide | 0 | — | — | wrong organism |
| IGVF differentiating-hPSC Perturb-seq (Stem Cell Reports 2025-11, 1,996,260 cells) | hPSC → cardiomyocyte and other lineages | ≤24 gene targets per library | ≤24/300 | — | — | the benchmarking paper behind 166's IGVF closure; confirms it |
| GEO 2025-06 → 2026-09, 167 human series | largest: GSE314342 (CD4, held), GSE295214 (held), GSE271853, GSE291147 (143 genes), GSE306853 (5 genes) | ≤520 | — | — | — | **no genome-scale screen in a line we do not already hold** |
| `perturbai/wholebrain_crispr_atlas` (HF, 8M cells) | **mouse** brain, in vivo | genome-scale | 0 | CC-BY-4.0 | — | wrong organism |
| `xingjiepan/PerturbMulti` (HF) | imaging, 18 protein channels | — | — | CC-BY-4.0 | — | not a transcriptome readout |

**Phase 1 verdict: 166's conclusion survives the change of criterion.** Admitting primary cells,
iPSC-derived and differentiated systems, organoids, CRISPRa, base editing and non-10x
chemistries adds **no new usable source**. The one genuinely new panel-covering asset (the KOLF
raw atlas) is a better-conditioned version of the *worst* source we hold.

---

## 7. Phase 2 — the arm that fell out, and it is not an acquisition

### 7.1 The mechanism: the sizing rule is a coverage gate, and it is silent

`solve_cd4_fraction` (`scripts/make_submission_pm_cd4.py:1258`) sizes the residual by the ratio
of its **median row norm over every panel row** to the whitened base's. A source covering ≤50%
of the panel has median exactly 0, so `f_per_unit = 0` and the solved fraction is **0**. This is
200's silent zero, verified by 207 on h1-test, and reproduced here on a third surface: in
`resid218` the `telohaec` and `gse295214_h1` arms return `fraction 0.00000` and endpoints
**bit-identical to the baseline** at all six doses on all four folds.

Production's own code path *raises* on a fully-zero residual, so a CD4 build fails loudly — but a
**partial**-coverage source is not fully zero, it is 21% covered with a zero median, and nothing
anywhere catches that. **`resid218` reports `median_residual_row_norm`, `carrying_rows`,
`coverage_frac` and a `VACUOUS` flag on every one of its 476 arms** so the failure mode cannot
recur silently in this experiment.

Sizing on the **carrying rows** instead — `median(‖added_i‖ : ‖added_i‖ > 0)` — makes `f` mean
"on the rows this source covers, the added norm is `g`× the base's typical norm". For CD4
(93–98% coverage) the two definitions coincide; **CD4's sizing is never touched here.**

### 7.2 The arm: the H1 pair as a second residual on top of `cd4f034`

`compose218.py` → `work/compose218.json` (496 arms). Both residuals are solved against the
**same** denominator and added in **one** step followed by production's single global preserve
rescale. CD4 is pinned at the board incumbent's `VCC_CD4_TARGET_F = 0.34`. The second source is
the **union** of `gse295214_h1` (Pan 2026, 64/300) and `h1_2025` (VCC-2025 H1, 20/300) — they
share only 7 targets, so the union is **77/300 = 25.7%** of the panel, of which 67 are also
carried by CD4 and 10 are rows CD4 cannot reach at all.

Baseline is `cd4@0.34` alone. Per fold, relative gains:

| `g` | fold | rel Δ`cos_agg` | rel Δ**macro** cos | Δ`pds` |
|---:|---|---:|---:|---:|
| 0.17 | k562_gwps | +6.880% | +5.712% | +0.007509 |
| | orion_hct116 | +7.234% | +5.394% | +0.001881 |
| | orion_hek293t | +5.382% | +4.411% | +0.002693 |
| | kolf2 | +27.268% | +24.935% | +0.016371 |
| **0.28** | **k562_gwps** | **+10.912%** | **+8.481%** | **+0.009516** |
| | **orion_hct116** | **+11.507%** | **+8.053%** | **+0.002072** |
| | **orion_hek293t** | **+8.500%** | **+6.506%** | **+0.003072** |
| | kolf2 | +44.267% | +38.862% | +0.021628 |
| 0.45 | k562_gwps | +16.440% | +11.497% | +0.010542 |
| | orion_hct116 | +17.421% | +10.996% | +0.001165 |
| | orion_hek293t | +12.724% | +8.701% | +0.001242 |
| | kolf2 | +69.273% | +56.816% | +0.026566 |

**12 of 12 rows positive on all three endpoints.** The `kolf2` fold is 3–5× the others because
its recipient is pluripotent and both sources are pluripotent (§5's state-match effect); it is
**excluded from every headline number below** and reported only for completeness.

### 7.3 The four falsifiers

| # | test | result |
|---|---|---|
| 1 | **target-shuffled** second source (whole rows permuted by a derangement; geometry, per-row norms, support pattern and carrying count preserved exactly) | **strictly negative at every dose**: macro **−0.596% / −1.422% / −3.163%** at `g` = 0.17 / 0.28 / 0.45, and monotonically worse. No permutation tie — the failure mode that killed 204's arm |
| 2 | **macro cosine must agree with the energy-weighted one** (an arm touching 21% of rows could move `cos_agg` by inflating those rows) | agrees at every dose: **+7.680% macro against +10.307% agg** at `g` = 0.28 across the cancer folds. Both large, both positive |
| 3 | **row-matched CD4 control** — same rows, same sizing rule, same dose, CD4's information instead of the new source's; absorbs row selection and "more residual is better" entirely | reproduces only **1/3.9** of the macro gain on the cancer folds (**+1.970% vs +7.680%** at `g` = 0.28; **1/6.9** on the 4-fold macro) and **costs `pds`** (−0.003708) where the real arm **gains** it (+0.004887). Priced on the conservative endpoint the control is **negative, −0.00134** |
| 4 | **non-carrying rows must not move** | **+0.00%** at every dose on every fold — across 4 folds × 3 doses the largest excursion is **|Δ| = 0.018 pp** (`rowsplit218.py`). The carrying rows move **+8.8% to +23.2%** on the cancer folds |
| 5 | **leak-free `kolf2` fold** | included at full weight in the LOFO selection; excluded from the headline for the *opposite* reason (it is the most favourable fold, not the most leaky) |
| 6 | **non-vacuity** | every reported arm carries `median_second_row_norm > 0` and `frac_x > 0`, asserted in-loop; the sweep aborts otherwise |

A consistency check the arm was not designed to pass: at `g = 0.17` the **carrying rows** gain
**+8.8% to +11.7%** cosine. 159's `1/k` law predicts **+13.7%** for a 4th independent line
(`1/cos² = 66.5 + 1920/k`, `k` 3 → 4). **The arm behaves quantitatively like adding one more
independent donor on exactly the rows it covers**, which is what it is.

### 7.4 The price — all six members, §8.4's rule, never four at zero

Rate card, from `cd4f034`'s own board row (`210/BOARD-RESULT.md`) divided by its measured cosine
move (0.072280 → 0.090270, Δ = +0.017990). `mse`'s rate reproduces 197's independently derived
+0.00365 overall per +0.01 cosine to **0.5%**, which is the check that the card is not circular.

| member | `cd4f034` Δ scaled overall | per +0.01 absolute cosine |
|---|---:|---:|
| `mse` | +0.006537 | +0.003634 |
| `reach` | +0.004400 | +0.002446 |
| `nmae` | +0.001762 | +0.000980 |
| `fid` | +0.001524 | +0.000847 |
| `jac` | +0.000179 | +0.000100 |
| `pds` | −0.001345 | −0.000748 |
| **total** | **+0.013057** | **+0.007259** |

`pds` is also priced independently through 207's route (PLOLO raw Δ / 0.452 / 6 × 1.011, ×2.0
per §8.4) and the mixed total is what is reported.

| endpoint / frame | rel gain | Δ`pds` | price `k = 1` | price `k = 1.80` |
|---|---:|---:|---:|---:|
| `cos_agg`, 4 folds, `g = 0.28` | +18.797% | +0.009072 | +0.02035 | +0.03122 |
| `cos_agg`, 3 cancer folds, `g = 0.28` | +10.307% | +0.004887 | +0.01109 | +0.01705 |
| **macro cosine, 3 cancer folds, `g = 0.28` — HEADLINE** | **+7.680%** | **+0.004887** | **+0.00919** | **+0.01364** |
| macro, 3 cancer folds, `g = 0.17` | +5.172% | +0.004028 | +0.00674 | +0.00973 |
| macro, 3 cancer folds, `g = 0.45` | +10.398% | +0.004317 | +0.01073 | +0.01675 |
| **minus the free row-matched-CD4 alternative, `g = 0.28`** | | | **+0.01054** | |

**Gate (rule 13, shape-only, +0.0019):** α, γ, λ and the emission law are untouched; the only
renormalisation is the preserve scalar, which is present at every dose. The headline clears by
**4.8×**, the marginal-over-free-alternative by **5.5×**, and the amplitude-class gate
(+0.0089) by **1.03×** — so it clears under **either** classification, which is the point 207 §6.4
had to argue for and this arm does not. At `g = 0.45` the amplitude-class margin is 1.21×.

**Leave-one-fold-out** (207's falsifier 3 — the dose is chosen on the other three folds and
priced on the held one, selection and reporting frames never overlapping): **+0.02830**, with
every fold positive (`k562_gwps` +0.01974, `orion_hct116` +0.01346, `orion_hek293t` +0.01012,
`kolf2` +0.06988), all four selecting `g = 0.45`. The single components are `gse295214_h1`
+0.02151 and `h1_2025` +0.01207; the union is close to their sum, as mostly-disjoint row
coverage predicts.

### 7.5 What this reopens, explicitly

**Experiments 117 and 118 closed Pan-2026 H1 as a residual — on a `pds`-only price, before
§8.4.** 117's own held evaluation measured raw-PDS **+0.0023 macro on Jurkat** and +0.00006 on
HEK293T, all positive, and 118 recorded "source-only signal confirmed in aggregate, production
candidate NO-GO". §8.4 states that a direction arm priced that way is understated by **≈4.5×**
and that any NO-GO whose margin sits inside 4.5× must be re-read. **This is one of them, and it
was also being run through a sizing rule that set its dose to zero.**

---

## 8. What I did NOT measure, and the honest uncertainties

- **The dominant uncertainty is the PLOLO → board transfer for a 26%-coverage arm.** §8.4's
  1.80× factor and 207's 1.26× were both measured on arms that move ~93% of the panel rows.
  This one moves 25.7% of them, a lot. The **uncalibrated `k = 1` column is the headline**
  precisely because of this, and the macro (per-target, unweighted) endpoint is used rather than
  the energy-weighted one for the same reason — the board averages six metrics *per
  perturbation*, so an energy-weighted aggregate can over-credit a concentrated gain.
- **The rate card rests on one board point.** §8.4's own lesson is that a class constant
  estimated from one point is a **lower bound on magnitude** — which cuts both ways: the `pds`
  cost may be understated as easily as the gains.
- **No `cell-eval2` run, no build, no submission, no `fid`/`jac`/`reach` measured directly.**
  Rule 9: those three have `n_real` in a denominator or a union and are unreadable on every
  harness we hold. They are priced through the rate card and nowhere else.
- **`h1_2025` covers 85/85 of `h1-test`'s targets.** It leaks into *that harness* and into
  nothing here (it is neither a PLOLO source nor a PLOLO truth, and the challenge panel is not
  H1 data). **Any future screen of this arm must not use `h1-test`.** Flagging it because the
  same trap has already cost this project one experiment (157).
- **`gene_norm` vs `raw` was swept for every source and both agree in sign and ordering**; only
  `gene_norm` (production's CD4 mode) is tabulated above.
- **The `kolf2` fold is excluded from every headline**, though it is the most favourable. Its
  truth is a per-gene z-score matrix with a placeholder `b_ctrl` and the weakest absolutes of the
  four folds (`cos_agg` 0.0098 against 0.022–0.035).
- **The 189 GB KOLF pan-genome file was not downloaded**, so its exact panel coverage is inferred
  from the webapp panel (10,167 targets, 273 on the panel) rather than measured on the raw file.
  The 46.7 GB product's 55/300 **was** measured, remotely.
- **A coverage-aware sizing rule is a builder change and I did not make it.** Nothing under
  `scripts/`, `src/`, `data/donors/` or `submissions/` was touched.

---

## 9. Handoffs

1. **To whoever owns the builder:** `solve_cd4_fraction`'s median-over-all-rows sizing is a
   silent coverage gate. Two sources already on disk are identically absent from every build
   because of it. A `median over carrying rows` variant, applied to partial-coverage sources
   only and leaving CD4's own solve untouched, is what §7 measures. Production's fully-zero guard
   does not catch the partial case.
2. **To 216 / 217:** nothing here touches the CD4 dose or form. The second residual is priced
   **against `cd4@0.34`**, and 208's rule says a rider must be re-priced against the *current*
   incumbent — so if 217's `f = 0.80` lands, this arm needs re-pricing against that, not against
   0.34. Both `cd40.34+…` and `cd40.10+…` families are in `work/compose218.json` for exactly
   that reason.
3. **To 213 (transport law):** the carrying rows of a 4th well-aligned source gain +8.8% to
   +11.7% cosine at `g = 0.17`, against 159's `1/k` prediction of +13.7% for `k` 3 → 4. That is
   an independent, per-row confirmation of the law on data it was not fitted to.
4. **Donor acquisition is now closed under the right criterion.** Not because independence is
   satisfied, but because independence was never the binding variable: alignment is, coverage
   gates it, and no unheld source improves either. The remaining moves in this territory are
   (a) the sizing rule, measured here, and (b) re-extracting the CD4 atlas at genome scale,
   which is 216's.

---

## 10. Files

| script | what it produces |
|---|---|
| `sources218.py` | the candidate registry, per-fold panel coverage, and the median-rule vacuity prediction → `work/sources218.json` |
| `resid218.py` | every source as an un-whitened post-whitening residual, 2 normalisations × 6 doses × 4 folds + controls, 476 arms → `work/resid218.json` |
| `independence218.py` | the independence matrix in the residual frame + 166's KA gate in a subprocess → `work/independence218.json`, `work/ka218.json` |
| `matched218.py` | 166 §4.2's confound control: per-target cosines and matched-row alignment ratios → `work/matched218.json` |
| `regress218.py` | **the load-bearing regression**, the relaxation floor, within-fold Spearman → `work/regress218.json`, `work/contrib218.json` |
| `compose218.py` | the second-residual arms, the union arm, the carrying-row sizing rule, shuffled and row-matched-CD4 controls, 496 arms → `work/compose218.json` |
| `diag218.py` | macro-vs-aggregate, the row-matched control table, LOFO dose selection → `work/diag218.json` |
| `rowsplit218.py` | carrying vs non-carrying row decomposition → `work/rowsplit218.json` |
| `price218.py` | the six-member rate card and the gate → `work/price218_cd4{0.34,0.1}.json` |
| `census166_copy.py` | **166's `census.py`, byte-identical**, run from here → `work/census.json` |

```bash
cd experiments/218-independent-donors
uv run python sources218.py
uv run python resid218.py --arms all          # 476 arms, ~25 min, peak RSS 7.65 GiB
uv run python independence218.py              # KA gate runs in its own subprocess
uv run python matched218.py
uv run python regress218.py                   # the load-bearing result
uv run python compose218.py                   # 496 arms
uv run python diag218.py && uv run python rowsplit218.py && uv run python price218.py
uv run python census166_copy.py data/external
```

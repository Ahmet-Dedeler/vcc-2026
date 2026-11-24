# 293 — RESULTS: the guide-construct ceiling in the challenge's OWN regime (2025 H1 library, 10x Flex)

Preregistration `PREREGISTRATION.md`, sha256 `43465b6fbaeb5aff962f6c0794769b9e13e291a430673c2a01a5faa628de1b14`
(`work/PREREGISTRATION.sha256`), hashed **2026-09-06T16:02:51Z with only `obs` metadata read** — before any
expression value entered the experiment. Analysis only: no model, no build, no submission, nothing written
outside this directory. Peak RSS 7.61 GiB (train), whole run ≈ 30 min wall clock.

---

## 0. The headline

**CHALLENGE-LIBRARY CONSTRUCT CEILING (H1, Flex, 400 cells, 20 k UMI/cell): ρ_guide +0.075
[−0.010, +0.288] on 16 targets; between-construct cosine at 400/400 (SB) +0.046 [−0.006, +0.121] vs
same-construct +0.476 [+0.260, +0.717]; scorer `pds` (a) 1.000 → (b) 0.677; knockdown explains
essentially all of it (+0.687 [+0.197, +1.067] where both constructs reach |lfc| ≥ 1.0, −0.003
[−0.085, +0.181] where the weaker does not).**

Three readings behind that line, in order of how much they should change what the programme believes:

1. **F14 reproduces in Arc's own library, and at the low end of it.** Two constructs of one gene in the
   2025 H1 CRISPRi library agree at **ρ_guide 0.075** of what two halves of one construct's cells agree at
   — beside 285's 0.09 (K562-GWPS), 0.17 (Jurkat), 0.41 (HepG2), 0.02–0.06 (CD4), 0.76–0.81 (Broad-16
   Cas9). P-2 (the frozen prediction that a knockdown-validated dual-guide library would beat GWPS's 0.09,
   point guess 0.2–0.5) is **REFUTED**: Arc's library reads at GWPS's level, not above it.
2. **At the pooled median, two constructs of the same gene agree no better than two constructs of
   different genes.** The preregistered design had no cross-target null; one was added during the run and
   it is the most informative row in the file. (x) — construct 1 of target *t* against a construct of a
   *different* target — reads **+0.0167**, and **(b) − (x) = +0.0067 [−0.0333, +0.0825]**, an interval
   through zero at every rung and both depths. The scorer's rank statistic still separates them
   (`pds` (b) 0.677 vs (x) 0.539 vs (a) 1.000), because a rank can survive a cosine that a median cannot.
3. **The construct term here is knockdown failure, not construct-specific biology — and the targets that
   can measure it are the library's worst.** Split the 16 targets by the weaker construct's own on-target
   knockdown and ρ_guide goes from **−0.003 [−0.085, +0.181]** (11 targets, weaker construct |lfc| < 1.0)
   to **+0.687 [+0.197, +1.067]** (5 targets, both ≥ 1.0). But two-construct targets are *not* a random
   sample of the library: their 34 constructs knock down significantly worse than the 283 single-construct
   targets (38.2 % vs 16.6 % below |lfc| 1.0, Mann–Whitney **p = 0.0015**). **So 0.075 is a lower bound on
   the library-wide construct ceiling, and the knockdown-matched +0.687 is the estimate that applies to the
   ~83 % of targets whose single construct does work.** This is a genuine two-sided answer, and §6 states
   both.

---

## 1. Known-answer rows — 4/4 PASS, three of them at 0.0e+00 (`work/ka293.{py,json,log}`)

| row | required | measured | verdict |
|---|---|---|---|
| **KA-1** 273's lane-parity 200/200 residual median on `h1-test.h5ad` (85 targets, 48 lanes, disjoint control halves 9,208/9,192, 17,995 scored columns) | `0.24291977562969289` | `0.24291977562969289` (mean 0.330950 vs 273's 0.330950) | **PASS, |d| 0.0e+00** |
| **KA-2** 273's 5-target construct row from the same file — median, all five per-target resid values, all five `full` values, all five construct cell counts | median `0.013512255772155394`; DPF2 −0.031551, GRB10 +0.013512, MED15 −0.099169, PCCB +0.027146, SOX4 +0.353942; cells 37/363, 191/209, 150/250, 351/49, 338/62 | identical in every one of the 20 numbers | **PASS, max |d| 0.0e+00** |
| **KA-3** integer identity: lane half A + lane half B == the full row | 0 | counts |d| **0**, integer sums |d| **0.0e+00**, controls 0 / 0.0e+00 | **PASS** |
| **KA-4** the construct parser on `adata_Test.h5ad`'s obs == 273's five multi-construct targets; every h1-test barcode present in the source with a matching lane suffix | — | DPF2/GRB10/MED15/PCCB/SOX4 from the full split and from h1-test; barcodes found ✓; suffix == batch column ✓ | **PASS** |

`per_target` (the scorer's `pds_cosine`, == `cell-eval2` 0.16.0 to 5.6e-17 by 273's KA-2) is imported
verbatim from 218 through `lib285.py`; recorded sha256 in `work/ka293.json`.

### 1a. A fifth check that was not preregistered, and is the strongest one: the depth bridge

`h1-test.h5ad` is **depth-reduced** relative to the 2025 source files — median **19,965** UMI/cell against
**54,256** in `adata_Test.h5ad` (measured here). The 2026 challenge's own released control cells read
**20,340** (`data/challenge/context_A.h5ad`, 18,400 cells). So *h1-test is at the challenge's depth and the
2025 files are 2.7× deeper*, and "400 cells" only means something once the depth is stated.

Thinning every summed side by 19,965/54,256 = **0.3680** and re-running: the **test split's panel-wide
random-cell split at 200/200 over the same 85 targets reads +0.2457** against **273's lane-parity
+0.2429**. Different cells, a different split (random vs lane parity), a different centring pool, one
independent instrument — **1.2 % apart.** Everything below is therefore reported at both depths, with the
thinned pass labelled the challenge regime.

---

## 2. Inventory (metadata only; `work/inv293.py`, `work/inventory293.json`)

| split | cells | targets | constructs | NTC constructs | NTC cells | lanes | ≥ 2-construct targets |
|---|---:|---:|---:|---:|---:|---:|---:|
| train | 221,273 | 150 | 189 | 31 | 38,176 | 48 | **8** |
| validation | 98,927 | 50 | 85 | 31 | 38,176 | 48 | **4** |
| test | 170,846 | 100 | 136 | 31 | 38,176 | 48 | **5** |

**n = 17** (the mandate's expected 13 + validation's 4). The 38,176 non-targeting cells and their 31 NTC
constructs are the **same pool in all three splits**. Every construct spans all 48 Flex lanes except
`DPF2_P1` (42 lanes, 95 cells). Every target gene is on the 18,080-gene axis in every split. Parser check:
constructs = 31 NTC + targets + multi-extras in all three splits (189 = 31+150+8, 85 = 31+50+4,
136 = 31+100+5). 16 of the 17 pairs are two promoters (P1/P2); **KDM2B** is two transcript ids
(`ENST00000377069.4` / `ENST00000377071.4`). Cells per construct 95 … 3,415 (per-pair table in
`PREREGISTRATION.md` §1). Realised cohorts: **17** at n = 25 and 50, **16** at n = 100 and 200 (DPF2 out,
its P1 has 95 cells).

---

## 3. The matched ladder, pooled over the three splits

Each target is measured inside **its own split's** residual space (its own matched-n panel centring, its
own scored-column set), 8 draws averaged, paired target bootstrap 20,000 draws. (a) same construct's cell
halves; (b) construct 1 vs construct 2; (d) lane halves of one construct; (a2) the *other* construct's own
cell halves; (x) the cross-target null. Full table `work/construct293.csv`.

**Challenge regime (thinned to 19,965 UMI/cell):**

| n/side | k | (a) | (a2) | (b) | (d) | (x) null | (b) ≤ 0 | ρ_guide = med(b)/med(a) | `pds` (a) | (b) | (x) | SB(a)→400 | SB(b)→400 |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|
| 25 | 17 | +0.0596 | +0.0261 | +0.0053 | +0.0577 | +0.0032 | 0.41 | +0.089 [−0.078, +0.467] | 0.989 | 0.548 | 0.504 | +0.5035 | +0.0783 |
| 50 | 17 | +0.0986 | +0.0394 | +0.0085 | +0.0950 | +0.0016 | 0.18 | +0.087 [+0.020, +0.391] | 0.999 | 0.622 | 0.472 | +0.4668 | +0.0645 |
| 100 | 16 | +0.1901 | +0.0978 | +0.0261 | +0.1798 | +0.0218 | 0.31 | +0.137 [−0.007, +0.322] | 1.000 | 0.649 | 0.575 | +0.4842 | +0.0968 |
| **200** | **16** | **+0.3123** | **+0.1660** | **+0.0233** | **+0.3018** | **+0.0167** | **0.38** | **+0.075 [−0.010, +0.288]** | **1.000** | **0.677** | **0.539** | **+0.4759** | **+0.0456** |

**Native depth (54,256 UMI/cell), same 17 targets:**

| n/side | k | (a) | (a2) | (b) | (d) | (x) null | ρ_guide | `pds` (a) → (b) | SB(a)→400 |
|---:|---:|---:|---:|---:|---:|---:|---|---|---:|
| 25 | 17 | +0.1226 | +0.0697 | +0.0131 | +0.1155 | +0.0136 | +0.107 [−0.096, +0.553] | 0.994 → 0.561 | +0.6909 |
| 50 | 17 | +0.2144 | +0.1161 | +0.0193 | +0.1991 | +0.0184 | +0.090 [−0.010, +0.442] | 1.000 → 0.579 | +0.6858 |
| 100 | 16 | +0.3570 | +0.1941 | +0.0358 | +0.3475 | +0.0362 | +0.100 [+0.006, +0.332] | 1.000 → 0.669 | +0.6895 |
| 200 | 16 | +0.5192 | +0.3104 | +0.0387 | +0.4923 | +0.0378 | +0.075 [−0.043, +0.321] | 1.000 → 0.716 | +0.6835 |

Rows that hold at every rung and both depths:

* **ρ_guide is flat in the budget** — 0.075–0.137 across an 8× range of total UMIs and a 2.7× range of
  depth. That is the signature of a *ratio of a true correlation to a reliability*, not of a noise floor.
* **Lanes are free** (P-3 **CONFIRMED**): (a) − (d) = **+0.0105 [−0.0030, +0.0334]** at 200/200 challenge
  depth, ≤ +0.013 at every rung; (d)/(a) = 0.966. 285's "batch free" and 273's Spearman +0.994 hold on the
  2025 splits' own 48 lanes.
* **(b) is not separated from the cross-target null** at the median: (b) − (x) = +0.0020 / +0.0069 /
  +0.0043 / **+0.0067 [−0.0333, +0.0825]** at n = 25/50/100/200, every interval through zero; native depth
  −0.0006 / +0.0009 / −0.0004 / +0.0009. The scorer's rank statistic *does* separate them — `pds` (b)
  0.677 vs (x) 0.539 vs (a) 1.000 — so on the metric the board actually pays, swapping the construct costs
  **70 %** of the distance from a perfect same-construct match down to an unrelated gene, and keeps 30 %.
* **Sensitivity row passes**: (b) re-run with the fixed identity construct-1 → side A (rather than
  larger-construct → side A) reads +0.0258 against +0.0233 at 200/200, inside the draw noise.
* **Panel baselines** (the split's own reliability, all pool targets, challenge depth) —
  train +0.2558 (126 targets), validation +0.2637 (44), test +0.2457 (85) at 200/200; native
  +0.4405 / +0.4632 / +0.4288. The 16 two-construct targets' (a) = +0.3123 sits *above* the panel median,
  so they are not unusually weak responders.
* **Monte-Carlo stability**: three independent 8-draw realisations of the frozen design gave ρ_guide at
  200/200 of 0.060, 0.0747 and 0.0727 (a matched second seed reproduces (a) 0.3148 vs 0.3123, (b) 0.0229
  vs 0.0233, SB(a)→400 0.4788 vs 0.4759). The spread is ~⅟₃₀ of the bootstrap width.

---

## 4. Per target (challenge regime, n = 200/side; sorted by the disattenuated ratio)

`ρ_dis = (b)/√((a)·(a2))` corrects for the *second* construct's own reliability, so a construct with no
measurable response is no longer scored as a construct that disagrees. `|kd|` is the on-target knockdown
`log1p(50,000 p)` of the gene's own column at full depth — **excluded from every scored statistic**, so
nothing here is circular.

| target | split | cells big/small | (a) | (a2) | (b) | (d) | ρ_dis | `pds` (a)→(b) | kd big | kd small |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|
| SOX4 | test | 2785/558 | +0.610 | +0.469 | **+0.532** | +0.602 | **+0.995** | 1.000 → 1.000 | −2.094 | −1.189 |
| STAG2 | train | 911/662 | +0.665 | +0.807 | **+0.710** | +0.666 | **+0.968** | 1.000 → 1.000 | −1.149 | −1.795 |
| TMSB4X | train | 2499/2261 | +0.072 | +0.061 | +0.065 | +0.053 | **+0.968** | 0.842 → 0.745 | −3.037 | −2.737 |
| PRCP | train | 3415/916 | +0.150 | +0.133 | +0.103 | +0.130 | +0.730 | 1.000 → 0.987 | −1.569 | −1.539 |
| IGF1R | validation | 1614/442 | +0.356 | +0.213 | +0.171 | +0.323 | +0.621 | 1.000 → 0.999 | −1.086 | −0.662 |
| ZNF581 | validation | 1393/1109 | +0.233 | +0.074 | +0.055 | +0.221 | +0.416 | 1.000 → 0.560 | −1.748 | −0.169 |
| KLF10 | train | 1483/449 | +0.317 | +0.199 | +0.021 | +0.320 | +0.084 | 1.000 → 0.610 | −1.077 | −0.019 |
| SMARCB1 | validation | 1284/422 | +0.859 | +0.074 | +0.010 | +0.857 | +0.040 | 1.000 → 0.955 | −1.807 | −0.593 |
| GRB10 | test | 2136/1938 | +0.567 | +0.059 | −0.001 | +0.562 | −0.005 | 1.000 → 0.545 | −0.197 | −1.009 |
| ARPC2 | validation | 491/489 | +0.559 | +0.037 | −0.003 | +0.556 | −0.020 | 1.000 → 0.603 | −1.602 | −0.078 |
| KDM2B | train | 1148/994 | +0.308 | +0.568 | −0.013 | +0.284 | −0.030 | 1.000 → 0.745 | −1.355 | −0.772 |
| HIRA | train | 2888/519 | +0.268 | +0.882 | −0.050 | +0.255 | −0.103 | 1.000 → 0.482 | **+0.030** | −1.494 |
| CAST | train | 1358/894 | +0.550 | +0.051 | −0.026 | +0.538 | −0.157 | 1.000 → 0.269 | −0.237 | −0.482 |
| MED15 | test | 2178/1459 | +0.047 | +0.815 | −0.104 | +0.047 | −0.530 | 0.607 → 0.484 | −0.053 | −1.376 |
| DAXX | train | 982/350 | +0.130 | — | +0.026 | +0.110 | — | 1.000 → 0.775 | −1.184 | −1.672 |
| PCCB | test | 2358/393 | +0.055 | — | +0.049 | +0.044 | — | 0.729 → 0.478 | −1.595 | −0.187 |

(DPF2 qualifies only at n ≤ 50; DAXX_P2 and PCCB_P2 have < 400 cells so (a2) is undefined at n = 200.)

**The distribution is bimodal, not diffuse.** Median per-target ρ_dis **+0.062** [−0.030, +0.730] over the
14 measurable pairs: six pairs at 0.42–1.00 and six at ≤ 0. Nothing sits at the median.

**P-4 CONFIRMED in ordering, refuted in magnitude**: SOX4 is top-2 and MED15 is last, exactly 273's
ordering on the test five; the pooled 17-target median (+0.023 at 200/200 challenge depth) differs from
273's +0.0135 by less than its own bootstrap width, so the second half of P-4 does not hold.

**The two rows that matter most are HIRA and MED15.** In both, the *larger* construct (the one this
design uses as the reference) is the one that **does not knock the gene down** — HIRA_P1 `+0.030`,
MED15_P2 `−0.053` — and it nevertheless produces a **reproducible** response (a) = +0.268 and +0.047,
while the construct that *does* knock down is far more reproducible ((a2) +0.882 and +0.815) and points
somewhere else entirely ((b) −0.050 and −0.104). A construct that misses its target still delivers a
stable, repeatable transcriptome. That is what a 400-cell truth measures when the construct fails.

---

## 5. The mechanism, and the selection effect that bounds it (post-hoc, labelled)

### 5a. Knockdown of the weaker construct separates the two modes (`work/strat293.{py,json,log}`)

Challenge regime, n = 200/side:

| stratum | k | (a) | (a2) | (b) | ρ_guide | ρ_dis | `pds` (a)→(b) |
|---|---:|---:|---:|---:|---|---:|---|
| ALL | 16 | +0.312 | +0.166 | +0.023 | **+0.075 [−0.010, +0.288]** | +0.102 | 1.000 → 0.677 |
| both constructs \|lfc\| ≥ 1.0 | 5 | +0.150 | +0.301 | +0.103 | **+0.687 [+0.197, +1.067]** | +0.484 | 1.000 → 0.987 |
| weaker construct \|lfc\| < 1.0 | 11 | +0.317 | +0.137 | −0.001 | **−0.003 [−0.085, +0.181]** | −0.004 | 1.000 → 0.560 |
| weaker construct, upper half | 8 | +0.332 | +0.213 | +0.084 | +0.252 [+0.028, +0.836] | +0.314 | 1.000 → 0.971 |
| weaker construct, lower half | 8 | +0.293 | +0.074 | −0.002 | −0.006 [−0.218, +0.222] | −0.012 | 1.000 → 0.515 |
| both constructs *measurable* ((a2) ≥ 0.05) | 13 | +0.317 | +0.199 | +0.021 | +0.067 [−0.047, +0.365] | +0.084 | 1.000 → 0.745 |

Covariates across the 16 targets (Spearman, challenge regime, n = 200): the weaker construct's knockdown
**+0.626 (p = 0.009)** with the per-target ratio and **+0.644 (p = 0.007)** with (b) itself (signs stated
in *strength*, |lfc|; the raw lfc correlations are negative because knockdown is a negative number). The
mean of the two knockdowns **+0.706 (p = 0.002)**. Everything else is at chance: response magnitude
+0.253 (p = 0.345), |knockdown difference| +0.271 (p = 0.311), min cells per construct +0.112 (p = 0.680).
**P-6 CONFIRMED in sign and, unusually, exceeded in strength** — the frozen prediction said |ρ| < 0.6 and
"will not explain most of the spread"; it reads 0.63–0.71 and the stratification is the whole effect.

Note the row that rules out the lazy version of this story: restricting to pairs whose second construct is
merely **measurable** ((a2) ≥ 0.05, 13 of 16) does *not* rescue ρ_guide (+0.067). The failing constructs
are not silent — they are reproducible and wrong. It is knockdown, not signal-to-noise, that separates the
modes.

### 5b. The two-construct targets are the library's weakest — so 0.075 is a lower bound

On-target knockdown of every target in the three splits (its own gene's column, full depth):

| group | n | median lfc | \|lfc\| < 1.0 | \|lfc\| < 0.5 |
|---|---:|---:|---:|---:|
| single-construct targets | 283 | −1.515 | **16.6 %** | 3.2 % |
| the 34 individual two-construct constructs | 34 | −1.167 | **38.2 %** | 29.4 % |
| two-construct targets, cells pooled | 17 | −0.968 | 52.9 % | 23.5 % |

Mann–Whitney on |lfc|, the 34 constructs vs the 283 single-construct targets: **p = 0.0015**. A second
construct exists in this library disproportionately where constructs are weak — whichever way the causality
runs (a second promoter added because one did not suffice; two transcripts each needing their own vector).
**The construct term can only be measured where the library is worst.** Every ρ_guide in §3 is therefore a
*lower bound* on the library-wide value, and §5a's knockdown-matched +0.687 is the estimate that transfers
to the 83.4 % of targets whose one construct reaches |lfc| ≥ 1.0.

Across the whole 2025 panel (300 targets, pooled): median lfc −1.485, **18.7 % below |lfc| 1.0, 4.3 %
below 0.5, 0 % positive**. If the 2026 library's 300 single constructs fail at the same rate, roughly one
target in five carries a truth that is a reproducible response to something other than the intended
knockdown — and no model can predict those, because the response belongs to the vector, not the gene.

---

## 6. The ceiling, stated (`work/ceiling293.{py,json,log}`)

Variance model, stated so it can be checked: construct *k*'s true response is `c_k = g + u_k` (g the
gene-level mean over the construct distribution, u_k the construct-specific part, orthogonal in
expectation); the scored truth is a noisy 400-cell measurement `t_k = c_k + e`. With S = s²g + s²u + s²e:
(a) = (s²g+s²u)/S, (b) = s²g/S, ρ_guide = (b)/(a) = s²g/(s²g+s²u). A **perfect construct-level** prediction
scores `cos = √(a)`; a **perfect gene-level** prediction — the best any zero-shot entry can do, since no
2026 construct exists on disk and the model is never told which vector produced the truth — scores
`cos = √(b)`. (a) and (b) are Spearman–Brown extrapolated to the challenge's 400-cell budget (k = 400/n;
SB of a 200/200 split with k = 2 is one 400-cell measurement against another).

**Challenge regime, from n = 200/side, 16 targets:**

| quantity | value |
|---|---|
| reliability of a 400-cell truth against another 400-cell truth of the **same construct**, SB | **+0.476 [+0.260, +0.717]** |
| the same for **two different constructs** of the gene, SB | **+0.046 [−0.006, +0.121]** (via ρ_guide: +0.036 [−0.005, +0.137]) |
| ρ_guide | **+0.075 [−0.010, +0.288]** |
| fraction of a 400-cell truth's reproducible variance that is **construct-specific** (1 − ρ_guide) | **0.925** |
| ceiling for a **perfect construct-level** prediction, cos vs the scored truth | **+0.690 [+0.510, +0.847]** |
| ceiling for a **perfect gene-level (zero-shot)** prediction, cos vs the scored truth | **+0.213 [+0.000, +0.348]** |

**Knockdown-matched** (the regime 83.4 % of the library is in; ρ_guide +0.687 applied to the library-wide
SB(a)₄₀₀ = 0.476): perfect gene-level ceiling **cos ≈ +0.572** (interval on ρ_guide alone: +0.305 … +0.690).

So the honest bracket for what a zero-shot model can reach against a 400-cell H1 Flex truth is
**cos ≈ 0.21 … 0.57**, against **0.69** for a model that somehow knew the construct. At native 2025 depth
the same three numbers are 0.29 / — / 0.83 (the depth changes the ceiling; it does not change ρ_guide).

For scale, R273 §0.6 measured the incumbent's **per-target transfer cosine at ≈ 0.04 at the median** on
this reference, with the `pds` knee at 0.02. The construct ceiling is therefore **not** what is binding
today: even the pessimistic 0.21 is ~5× above where the pipeline sits. **The construct term is a real and
large ceiling, and it is not the current constraint.**

---

## 7. Against 285 (F14) and 273

| library / unit | ρ_guide | source |
|---|---|---|
| **2025 VCC H1, Arc's dual-guide library, Flex, challenge depth (this experiment)** | **+0.075 [−0.010, +0.288]**, 16 targets | §3 |
| the same, restricted to pairs where both constructs knock down | **+0.687 [+0.197, +1.067]**, 5 targets | §5a |
| K562-GWPS (748 two-construct perturbations) | 0.09 | 285 |
| Jurkat | 0.17 | 285 |
| HepG2 | 0.41 | 285 |
| CD4 (GWCD4i) | 0.02–0.06 | 285 |
| Broad-16, Cas9 (not CRISPRi) | 0.76–0.81 | 285 |

**F14 holds in the challenge's own library**, at the CRISPRi end of the range. Two amendments it earns:

* **F14's ρ_guide is confounded with library QC.** 285 read the construct as "the dominant irreducible term
  of a perturbation response". On this library the term is almost entirely *knockdown failure*: pairs where
  both constructs work read 0.687, pairs where one does not read −0.003. Whether 285's public-library
  values (GWPS 0.09, CD4 0.02–0.06) survive the same stratification is untested here and is the obvious
  next measurement — GWPS's single guides of unscreened efficiency are exactly the population where it
  would bite hardest, and Broad-16's 0.76–0.81 (Cas9, validated) is what the knockdown-matched stratum
  looks like.
* **273's n = 5 read was directionally right and quantitatively lucky.** 273's +0.0135 vs +0.243 becomes,
  on 17 targets and a matched design, +0.023 vs +0.312 at the same depth — the same story, and its five
  per-target signs reproduce (SOX4 top, MED15 bottom).

**285's "cells and UMIs are one axis" is independently reconfirmed here**: thinning the 2025 files by
0.368 reproduces 273's lane-parity h1-test reading to 1.2 % on the same 85 targets (§1a), and ρ_guide is
flat across a 2.7× depth change and an 8× cell change.

---

## 8. Frozen predictions, scored

| | prediction | outcome |
|---|---|---|
| **P-1** | ρ_guide at n = 100 below 0.5 | **CONFIRMED** — +0.137 [−0.007, +0.322] |
| **P-2** | ρ_guide above K562-GWPS's 0.09; point guess 0.2–0.5 | **REFUTED** — +0.075 (200/side) / +0.137 (100/side), at or below GWPS. The knockdown-matched stratum (+0.687) is *above* the guess, so the prediction was wrong about the population, not about validated constructs |
| **P-3** | (d) within ±0.03 of (a) at every n | **CONFIRMED** — max (a) − (d) = +0.0105 [−0.0030, +0.0334] |
| **P-4** | SOX4 top and MED15 bottom; pooled median differs from 273's +0.0135 by more than its bootstrap width | **HALF** — ordering confirmed; the pooled median (+0.023) is *within* the bootstrap width of 273's value |
| **P-5** | `pds` for (b) closer to chance than for (a), and the gap larger than the cosine gap implies | **CONFIRMED** — (a) 1.000 → (b) 0.677 against a cross-target null of 0.539: (b) keeps 30 % of the (x)→(a) distance while its *median cosine* keeps none |
| **P-6** | weaker construct's knockdown correlates with ρ_guide (Spearman > 0) but explains < 60 % of the spread | **CONFIRMED in sign, EXCEEDED in strength** — Spearman +0.626 (p = 0.009) and the stratification is essentially the whole effect |

---

## 9. Limits — read these before quoting anything above

1. **n = 17 targets, 16 at the rungs that matter, 5 in the knockdown-matched stratum.** Every interval is
   wide and the strata are small. The bimodality is visible in the raw per-target table (§4), not only in
   the summary — but 5 targets is 5 targets.
2. **The selection effect in §5b is the dominant caveat.** Two-construct targets are significantly weaker
   in knockdown than the rest of the library (p = 0.0015). The as-measured ρ_guide is a lower bound; the
   knockdown-matched value is an extrapolation from 5 pairs onto 283 targets.
3. **One line, one chemistry, one library.** H1 hESC, 10x Flex, Arc's 2025 CRISPRi panel. The 2026 panel
   shares 25/300 targets with it and is run in six *somatic* lines; 280 measured that primary CD4 T cells
   sit in a far lower reproducibility regime (+0.054 across donors). Nothing here says the 2026 contexts
   behave like H1 — §6's ceilings are H1-Flex-400 numbers.
4. **The (x) null is a cross-target pairing inside the same 17 targets**, so its cohort is small and its
   own interval is wide; "(b) is at the null" is a statement about medians with intervals through zero,
   not a proof of equality. The `pds` rows say the two are *not* identical.
5. **The variance model in §6 is a model.** It assumes g ⊥ u ⊥ e and exchangeable constructs. With two
   constructs per gene, s²g is estimated as their covariance; a systematic promoter effect shared by all
   P1 constructs would be absorbed into g and inflate the ceiling.
6. **Spearman–Brown from 200 → 400 is a 2× extrapolation** (the shortest one available) and 285 validated
   SB on real cell ladders to ±0.02 on this and three other sources; the 25 → 400 extrapolations in §3 are
   16× and are reported only for the flatness check.
7. **The knockdown is a bulk pseudobulk lfc on one gene column**, not a per-cell efficiency; a construct
   with 50 % of cells fully knocked down and 50 % untouched reads the same as uniform partial knockdown.
8. The (a2) column and the (x) null were **added during the run**, after the preregistered contrasts were
   computed and before any of them was interpreted; both are labelled as such. The three seed realisations
   in §3 quantify what that re-running cost in Monte-Carlo terms (ρ_guide 0.060 / 0.073 / 0.075).

---

## 10. What this means for the programme (no action proposed)

* **Nothing to build.** This is an instrument-and-ceiling measurement. No arm, no dose, no slot.
* **The construct ceiling does not currently bind.** The pessimistic zero-shot ceiling (cos 0.21) is ~5×
  the incumbent's median per-target transfer cosine (≈ 0.04, R273 §0.6). The 2026 panel sharing only
  25/300 targets with the 2025 library remains true and remains the reason no source on disk carries the
  2026 constructs — but the number that reading implied ("the construct term binds every entry") should be
  restated as **"the construct term is a ceiling at cos 0.21–0.57, well above the field."**
* **The finding that is worth its own line in FINDINGS: a CRISPRi construct that misses its target still
  produces a reproducible transcriptome.** HIRA_P1 (lfc +0.030) and MED15_P2 (lfc −0.053) have split-half
  reliabilities of +0.268 and +0.047 while their knocking-down partners read +0.882 and +0.815 and point
  elsewhere. At the 2025 panel's rate (18.7 % of targets below |lfc| 1.0), roughly one 2026 truth in five
  may be a reproducible response to the vector rather than to the gene. That is a property of the
  *benchmark*, not of any model, and it caps every entry equally.
* **The obvious follow-up, if a slot of attention is available:** apply §5a's stratification to 285's five
  public libraries. If GWPS's 0.09 is also mostly knockdown failure, F14's headline ("the guide construct
  is the dominant irreducible term") should be rewritten as "**guide efficacy** is the dominant irreducible
  term", which is a different — and more actionable — claim.

---

## Files

`PREREGISTRATION.md` (+ `work/PREREGISTRATION.sha256`) · `work/inv293.py`, `work/inventory293.json` ·
`work/lib293.py` (imports `experiments/285-reliability-benchmark/work/lib285.py` verbatim) ·
`work/ka293.{py,json,log}` · `work/run293.py`, `work/run293_<split>[_thin0.368][_s1].{json,npz,log}` ·
`work/agg293.py`, `work/agg293[_thin0.368][_s1].{json,log}` · `work/ceiling293.{py,json,log}` ·
`work/strat293.{py,json,log}` · **`work/construct293.csv`** (132 rows = 17 targets × 4 rungs × 2 depths,
31 columns: (a), (a2), (b), (d), (x), the fixed-identity sensitivity row, all four `pds`, the un-centred
cosines, per-construct cells and knockdown, response magnitude, the split's panel median and the target's
rank in it, ρ_dis).

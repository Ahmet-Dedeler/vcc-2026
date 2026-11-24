# 238 — The October asset re-priced on a panel that excludes essential genes

**Verdict, one line:** *the October four are carried by **common-essential targets** — **70.7–90.4%**
of their cosine gain and **10 of 10 negative `pds`** from the non-essential half, on the same panel
at matched dose — and on a panel that excludes essential genes as **Arc's one observed panel does**
the asset's **robust** value falls to **16%** (+0.000257 vs +0.001607 absolute aggregate cosine
after 236's own leave-top-3-out) with **`pds` changing sign** (+0.001904 → **−0.000978**). **236's
+0.0066…+0.0104 does not survive; the strategy claim is retracted.** The one draw that still reads
large (`efp`, +0.0077…+0.0104) is **20.6% robust** and **45–76% of it is a single gene, `RUVBL2`** —
the largest true effect on that panel and a core AAA+ ATPase that is essential in function and
merely absent from Replogle's list.*

**Date 2026-09-04.** Deterministic (seeds 7 / 199 / 1990 / 23323 / 233233 / 9090+n). **420 scored
arm×fold×panel rows plus 19 vacuity records** over **5 panels × 3 cancer folds + the `kolf2` fold**,
**7 operator re-fits**
(all axis-asserted), **3 CD4 re-extractions**, **109/109 known-answer arms**. No `vcc`, no
`cell-eval2`, no build, no submission, nothing written outside
`experiments/238-essential-free-panel/`. `scripts/make_submission_pm_cd4.py` and
`scripts/build_union_whitener.py` were never opened (235 owns them); the latter was *invoked* as the
operator fitter, which is how the runbook fits operators.

**Memory.** Gate checked before every heavy job and at every launch: `memory_pressure` 83–86% free,
`vm_stat` free+inactive 30–44 GiB, swap free 1.06–1.16 GiB. Jobs run **strictly serially**. Peak
`ru_maxrss` (`/usr/bin/time -l`): operator fits **11.96–12.33 GiB**, scoring **7.34 / 7.26 GiB**,
influence **6.07 GiB**, transferability **1.26 GiB**. The fits exceed the 10 GB guidance because
`ru_maxrss` counts the 4.2 GB basis memmap while resident; the anonymous working set is ~8 GiB, and
they were run at the runbook's own `--max-rss-gib 14`. `DEF-RUNBOOK` §2 records the same pair
(8.1 GiB at `eigh`, 12.9 GiB process high-water), as does 233 (8.62 / 12.20).

---

## 0. Phase 0 — the fact, established before anything was scored

### 0.1 The `0 / 34 / 30` row, reproduced, with its list provenance

**The list is `data/donors/k562_essential/perts.txt`** — `n = 2,021`,
sha256 `b1e5279534eb240f0f3487f0814d379e05b12799b21f6cb3818f5c7f55557e2d`. It is the perturbation
set of `data/external/scperturb/ReplogleWeissman2022_K562_essential.h5ad`, i.e. **Replogle 2022's
K562 essential-gene Perturb-seq screen**, built by `scripts/build_donor_table.py` (asserted at
runtime against that bank's own `meta.json`). It is the **same** list that
`build_union_whitener.py --basis-prune essential` uses to drop basis rows and that
`scripts/def_census.py` reports the panel overlap against, so this is a reproduction of the
published row and not a different statistic wearing its name.

| panel | n | ∩ essential | % | `DEF-RUNBOOK` §5 | verdict |
|---|---:|---:|---:|---:|---|
| `abc` — Arc's real 300 | 300 | **0** | 0.0 | 0 | **PASS** |
| `real` — 199 realistic draw | 300 | **34** | 11.3 | 34 | **PASS** |
| `pess` — 199 pessimistic draw | 300 | **30** | 10.0 | 30 | **PASS** |

Disjointness asserted, not assumed: `|abc ∩ real| = |abc ∩ pess| = 0`. Nine per-source coverage
rows also reproduce 225 §8.3 / `DEF-RUNBOOK` §5.-2 exactly (`k562_gwps` 269 → 150, Orion 291 / 293,
the four October sources 0 → 30–33).

### 0.2 236's headline reproduces — 109/109 known-answer arms at max |Δ| **0.0e+00**

`six238.py` is `236/six236.py` with the panel registry swapped; `bind238.py` rebinds
`panels233`'s attributes so `fold233` / `score233` are executed **verbatim**. Every arm for which
233 or 236 stored a value was re-run and asserted against it (`--ka` aborts on failure):

| source of the stored value | arms checked | max &#124;Δ`cos_agg`&#124; / &#124;Δ`pds`&#124; |
|---|---:|---|
| `236/work/six236.json` + `stack236.json`, `233/work/score233.json` + `oct233.json` | **109** | **0.0e+00** (tolerance 1e-12) |

236 §3.2's published table, re-derived here through `price238.py` (228's exact routes):

| panel | arm | rel `cos` | Δ`pds` | Δ`nmae` | `e` | composite (no reach) | 236 §3.2 |
|---|---|---:|---:|---:|---:|---:|---|
| pess | all four @0.34 | **+18.88%** | +0.00114 | +0.00297 | 0.180 | **+0.00664 … +0.00926** | identical |
| real | all four @0.34 | **+12.54%** | +0.00415 | +0.00030 | 0.433 | **+0.00637 … +0.00811** | identical |

**236 is internally sound. Nothing below is a correction to 236's arithmetic.**

### 0.3 The dependence, quantified before the new panel existed

**Essential fraction of each source's carrying rows, per fold, on all three existing panels**
(falsifier 4; the whole matrix is a known-answer row against 236 §1 — **28/28 exact**):

| panel | fold | `nadig_hepg2` | `nadig_jurkat` | `rpe1` | `k562_essential` | UNION |
|---|---|---|---|---|---|---|
| `abc` | all three | **0 VACUOUS** | **0 VACUOUS** | **0 VACUOUS** | **0 VACUOUS** | **0** |
| `real` | `orion_hct116` | 45 (**69%** ess) | 48 (**67%**) | 47 (**70%**) | 34 (**100%**) | 50 (**68%**) |
| `real` | `orion_hek293t` | 45 (69%) | 48 (67%) | 47 (70%) | 34 (100%) | 50 (68%) |
| `real` | `k562_gwps` | 45 (69%) | 48 (67%) | 47 (70%) | *(leak)* | 50 (68%) |
| `pess` | `orion_hct116` | 28 (**86%**) | 28 (**82%**) | 27 (**81%**) | 25 (**100%**) | 30 (**83%**) |
| `pess` | `orion_hek293t` | 29 (86%) | 29 (83%) | 28 (82%) | 26 (100%) | 31 (84%) |
| `pess` | `k562_gwps` | 32 (88%) | 32 (84%) | 31 (84%) | *(leak)* | 34 (85%) |

**67–88%, and 100% for `k562_essential`.** The mandate's "if it is ~100%, the conclusion follows
immediately" is nearly met; the remaining 12–33% is what Phase 2 had to measure rather than assume.

**And the version of the statement that needs no list at all.** All four are common-essential
libraries, so the union of their own perturbation lists — **2,389 targets, 84.6% of them on the
essential list** — is the entire set of targets the October asset can *ever* carry:

| panel | &#124;panel ∩ October-source union&#124; |
|---|---:|
| **`abc` — Arc's one observed panel** | **0 of 300** |
| `real` | 50 |
| `pess` | 35 |
| `efr` (new) | 9 |
| `efp` (new) | 9 |

**On the only real panel anyone has ever seen, the October asset carries nothing** — measured, with
no exclusion applied by me. That is 233's own KA-G, and it is the strongest single fact in this file.

---

## 1. Phase 1 — the honest panel

`bundle238.py` is `199/defbundle199.py` and `defbundle199r.py` with **exactly one line changed**:
the eligible pool has the common-essential list subtracted before the draw. Same seeds, same gene
axis, same manifest surgery, same control files.

| | `efp` (from 199's pessimistic pool) | `efr` (from 199's realistic pool) |
|---|---|---|
| pool | union of the three production donors' coverage | K562 ∩ Orion-HCT116 ∩ Orion-HEK293T ∩ CD4 atlas |
| | minus A/B/C, **minus the 2,021 essential genes** | minus A/B/C, **minus the 2,021 essential genes** |
| eligible before → after | 17,280 → **15,463** (removed 1,817) | 7,280 → **6,402** (removed 878) |
| seed | 199 (199's own) | 1990 (199r's own) |
| **&#124;panel&#124;** | **300** | **300** |
| **&#124;panel ∩ abc&#124;** | **0** ✅ | **0** ✅ |
| **&#124;panel ∩ essential&#124;** | **0** ✅ | **0** ✅ |
| ∩ 199's sibling draw | 3 | 21 |
| control `.h5ad` + `gene_names.csv` vs 199's | **4/4 byte-identical** | **4/4 byte-identical** |

Both properties are **asserted at build time**, not checked afterwards. The panels match Arc's one
observed panel on both.

**Coverage of the new panels** (the price of the exclusion, paid by everything, not just the
October sources):

| panel | k562_gwps | orion_hct116 | orion_hek293t | CD4 union |
|---|---:|---:|---:|---:|
| `abc` | 269 | 290 | 294 | 291 |
| `real` → `efr` | 300 → **300** | 300 → **300** | 300 → **300** | 300 → **300** |
| `pess` → `efp` | 150 → **146** | 291 → **293** | 293 → **296** | 174 → **186** |

**The production donors barely move.** Only the October four collapse.

**CD4 re-extraction** — known-answer row: the `abc` extraction reproduces the shipped cloud banks
**36/36 data files byte-identical** and 283/287/284, union 291, 10,113 genes. 047's gates fire on
`efp` (`only 186 panel targets map … require 290`) and are relaxed to 0 per `DEF-RUNBOOK` §5's
pre-decided response; they pass on `efr`.

**Operators** — seven fits through the production fitter with 233's flags, one per (panel, fold):

| job | n_cols | axis sha256 (16) | seconds | peak `ru_maxrss` |
|---|---:|---|---:|---:|
| `ka` (abc, 3 donors) | 17,942 | `b2cae0d381c3cd70` | 319.2 | 12.33 GiB — **MATCH `e2f65ecf…`** |
| `efr:k562_gwps` | 17,806 | `dd84dcbb6af42eca` | 390.1 | 12.17 GiB |
| `efr:orion_hct116` / `efr:orion_hek293t` | 17,942 | `5b2767cf036150f4` | 308.1 / 307.0 | 12.33 / 12.13 GiB |
| `efp:k562_gwps` | 17,807 | `a8751a3c01be4efd` | 303.7 | 11.96 GiB |
| `efp:orion_hct116` / `efp:orion_hek293t` | 17,942 | `b327501c5344e5ab` | 307.9 / 306.3 | 12.33 / 12.13 GiB |

> ⚠️ **`scripts/build_union_whitener.py` changed under this run** (235 owns it; mtime 12:09:41,
> sha `c36a7672…` → `9ea64b12…`). The `ka` fit above ran the **old** version; all six measurement
> operators ran the **new** one. So the known-answer row was **re-run against the new version**:
> it also reproduces `e2f65ecf…` and writes the documented A/B/C axis
> `edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f`
> (`work/whiten-ka2.log`, 13.68 GiB, deleted after verification). **Both builder versions are
> gated; all seven operators are verified.**

---

## 2. Phase 2 — the four sources on an essential-free panel

### 2.1 `k562_essential` is structurally VACUOUS, and the rest collapse to 9 rows

| panel | fold | `nadig_hepg2` | `nadig_jurkat` | `rpe1` | `k562_essential` | UNION |
|---|---|---:|---:|---:|---:|---:|
| `efr` | `orion_hct116` / `hek293t` | 9 | 9 | 9 | **0 VACUOUS** | **9** |
| `efr` | `k562_gwps` | 9 | 9 | 9 | *(leak)* | **9** |
| `efp` | `orion_hct116` | 9 | 7 | 8 | **0 VACUOUS** | **9** |
| `efp` | `orion_hek293t` | 9 | 7 | 8 | **0 VACUOUS** | **9** |
| `efp` | `k562_gwps` | 8 | 7 | 7 | *(leak)* | **8** |

**Carrying rows 30–31 → 9 and 50 → 9.** `k562_essential` reads **0 carrying rows on 4 of 4
non-leak folds** and the instrument refuses to price it — tautologically, since its own
perturbation list *is* the exclusion list, which is stated here rather than hidden. All three
sources carry **the same 9 rows**: the redundancy 236 measured at 71–72% is now ~100%.

### 2.2 The re-priced arms — null-corrected, CLEAN2, `g` = 0.34 (233's pre-registered dose)

| panel | arm | rows | rel `cos` | Δ`pds` | Δ`nmae` | `e` | **composite (no reach)** | with reach |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `pess` | all four | 30 | **+18.88%** | +0.00114 | +0.00297 | 0.18 | **+0.00664 … +0.00926** | +0.00994 … +0.01256 |
| `real` | all four | 50 | **+12.54%** | +0.00415 | +0.00030 | 0.43 | **+0.00637 … +0.00811** | +0.00902 … +0.01076 |
| **`efp`** | 3-stack (all that carry) | **9** | **+19.78%** | +0.00057 | −0.00054 | 0.70 | **+0.00770 … +0.01045** | +0.01097 … +0.01372 |
| **`efr`** | 3-stack (all that carry) | **9** | **+1.51%** | **−0.00138** | −0.00011 | 0.36 | **+0.00008 … +0.00029** | +0.00044 … +0.00065 |

Six-member board-unit split (null-corrected, CLEAN2):

| panel | `o_pds` | `o_mse` (k = 1 … 1.384) | `o_nmae` | `o_fid` | `o_jac` | [`o_reach`] |
|---|---:|---|---:|---:|---:|---:|
| `pess` | +0.00043 | +0.00683 … +0.00946 | −0.00083 | +0.00020 | +0.00001 | +0.00330 |
| `real` | +0.00155 | +0.00454 … +0.00628 | −0.00008 | +0.00029 | +0.00007 | +0.00265 |
| `efp` | +0.00021 | +0.00716 … +0.00991 | +0.00015 | +0.00019 | −0.00002 | +0.00327 |
| `efr` | **−0.00051** | +0.00055 … +0.00076 | +0.00003 | −0.00000 | +0.00002 | +0.00036 |

**The two essential-free draws disagree by 30×, and 236's own falsifier settles it.** Both are dose-monotone
(`efp` +2.92 / +5.84 / +11.66 / +19.78% at `g` = 0.05 / 0.10 / 0.20 / 0.34; `efr` +0.23 / +0.45 /
+0.90 / +1.51%), so the disagreement is not a dose artefact. **Relative** cosine is the misleading
column here — `efp`'s BASE `cos_agg` is 0.0174 / 0.0119 against `pess`'s 0.0212 / 0.0172, which
inflates any ratio. The absolute gain, with 236's leave-top-3-out applied:

| panel | rows | absolute Δ`cos_agg` | minus top 3 | **retained** | absolute Δ`pds` |
|---|---:|---:|---:|---:|---:|
| `pess` | 31 | +0.003614 | +0.001413 | **39.1%** | +0.000832 |
| `real` | 50 | +0.002592 | +0.001800 | **69.4%** | +0.002975 |
| **`efp`** | 9 | +0.002820 | **+0.000580** | **20.6%** | +0.000125 |
| **`efr`** | 9 | +0.000194 | **−0.000066** | **−33.8%** | **−0.002081** |

```
ROBUST (top-3-removed) absolute cosine gain
  199's draws (mean)          +0.001607
  essential-free draws (mean) +0.000257      =  16.0%
absolute d_pds
  199's draws (mean)          +0.001904
  essential-free draws (mean) -0.000978      =  SIGN FLIP
```

236 measured the same falsifier on its own panels and found the cosine gain **"never inverts, in
18 of 18 rows"** at ~55% mean retention. **On the essential-free panels it retains 20.6% and
−33.8%, and it inverts.** By 236's own robustness criterion neither draw carries an asset.

### 2.3 What carries the `efp` residual: three near-essential genes, and mostly one

`efp`'s nine carriers are **ANKRD17, ESYT1, INTS11, RUVBL2, SPEN, TADA1, TAF11, TBPL1, UTP11** —
Integrator, TFIID, TBP-like, SAGA, SSU-processome and R2TP components: core transcription and
ribosome-biogenesis machinery that is essential in function and simply absent from Replogle's
K562-essential library. `efr`'s nine are **ATP11B, ATP5PO, GSK3B, MEIS3, MMP17, MON1A, MYSM1, NF1,
RCCD1** — genuinely peripheral, and that panel reads ≈0.

| fold | cosine gain | top row | its share | ‖truth‖ rank of that row |
|---|---:|---|---:|---|
| `efp` / `orion_hct116` | +0.003472 | **RUVBL2** | **45%** | **4 of 293** (‖truth‖ 15.83 vs median 6.94) |
| `efp` / `orion_hek293t` | +0.002168 | **RUVBL2** | **76%** | **1 of 296** (‖truth‖ 22.26 vs median 6.20) |
| `efr` / `orion_hct116` | +0.000353 | ATP11B | 49% | — |
| `pess` / `orion_hct116` | +0.003785 | PRPF40A *(essential)* | 29% | — |
| `pess` / `orion_hek293t` | +0.003443 | RNF113A *(essential)* | 36% | — |

**On `pess` the top influences are the panel's essential genes** (PRPF40A, RNF113A, GTF3A, DIS3,
TMED2, PRPF18 — spliceosome, exosome, Pol III machinery). **On `efp` they are the near-essential
genes the exclusion list happened to miss, and the single largest true effect on the panel is one
of them.** That is the answer to "if the sources survive, say what carries them": they survive by
the same biology, sampled at the edge of the exclusion, in one or two rows out of 293.

### 2.4 The cleanest evidence — the within-panel decomposition, no panel change at all

Comparing `pess` with `efp` changes the panel *and* the targets. Restricting the **same** residual
on the **same** panel to its essential vs its non-essential carrying rows changes only which rows
carry it, and `solve_fraction_carrying` re-solves the dose on whichever rows survive, so the
**per-row magnitude is matched between the halves**.

| panel | arm | rows | rel `cos` | **share of the full arm** | Δ`pds` |
|---|---|---:|---:|---:|---:|
| `pess` | all four | 30 | +19.83% | 100% | +0.00168 |
| `pess` | └ **essential rows only** | 25 | **+15.83%** | **79.8%** | **+0.00199** |
| `pess` | └ **non-essential rows only** | 5 | +4.32% | 21.8% | **−0.00025** |
| `pess` | └ derangement (null) | 30 | +0.95% | 4.8% | +0.00054 |
| `real` | all four | 50 | +12.77% | 100% | +0.00306 |
| `real` | └ **essential rows only** | 34 | **+10.50%** | **82.2%** | **+0.00318** |
| `real` | └ **non-essential rows only** | 16 | +2.17% | 17.0% | **−0.00017** |
| `real` | └ derangement (null) | 50 | +0.23% | 1.8% | −0.00109 |

Over all ten (panel, source) cells measured this way — four singles and two stacks on each of the
two panels:

```
essential half:      70.7% - 90.4% of the cosine gain, from 22-34 rows;  d_pds POSITIVE  8/10
non-essential half:  10.1% - 28.9% of the cosine gain, from  4-16 rows;  d_pds NEGATIVE 10/10
```

**Ten of ten non-essential halves deliver negative `pds`.** `pds` is the member the board is
decided on.

### 2.5 The mechanism, measured model-free

233 §6's transferability ceiling — `cos(source_i, source_j)` on a target's own row, using neither
the truth nor the prediction, so it is a property of the panel and not of our machinery:

| panel | fold | ceiling, **essential** targets | ceiling, **non-essential** | ratio |
|---|---|---:|---:|---:|
| `real` | k562 / hct116 / hek293t | 0.07605 / 0.05219 / 0.05475 | 0.02767 / 0.01261 / 0.01468 | **2.75× / 4.14× / 3.73×** |
| `pess` | k562 / hct116 / hek293t | 0.12211 / 0.04634 / 0.05032 | 0.03998 / 0.01312 / 0.01618 | **3.05× / 3.53× / 3.11×** |

**A common-essential knockdown produces the same core stress / growth-arrest programme in every
cell line, so it transfers for free.** That is why the October sources look valuable on a panel
that contains such targets, and it is exactly the property Arc's panel excludes — for the assay
reason the mandate names: an essential knockdown kills the cell, so there is no response to score.

Panel-level means (mean ‖truth‖ shows no confound: 6.95–7.40 across all five panels):

| panel | `abc` | `real` | `pess` | `efr` | `efp` |
|---|---:|---:|---:|---:|---:|
| mean ceiling | 0.02767 | 0.02316 | **0.03158** | 0.02014 | 0.02425 |
| vs `abc` | 100% | 83.7% | **114.1%** | 72.8% | 87.6% |

`pess` is the **most** transferable of the four synthetic panels **because** it carries 30 essential
genes. Removing them drops it to 87.6% of Arc's. Arc's own panel remains above every random draw,
which is 233 §6's finding, reproduced here from a different direction.

---

## 3. The retention number on an essential-free panel — it gets WORSE

233's projection re-run through the same code. **233's own rows reproduce** (real 64.6% vs its
published 64.6%; pess 55.5% vs 55.6%):

| member (credit above the deranged null, ÷ `abc`) | `real` | `pess` | **`efr`** | **`efp`** |
|---|---:|---:|---:|---:|
| `cos_agg` | 93.4% | 83.7% | **37.7%** | **59.9%** |
| `pds@140` | 64.5% | 54.5% | 57.9% | 56.3% |
| `mse*` | 86.9% | 70.9% | **15.2%** | **38.4%** |
| `nmae` ⚠ | 43.8% | 49.6% | 40.7% | 58.8% |
| `fid` | 66.9% | 62.5% | 35.9% | 55.2% |
| `reach` | 61.9% | 53.1% | 32.4% | 41.1% |
| `jac` | 52.5% | 51.8% | 31.0% | 34.6% |
| **DIRECTION** (`pds@140`, `mse*`, `nmae`) | **65.0%** | **58.3%** | **37.9%** | **51.2%** |
| **ORDERING** (`fid`, `reach`, `jac`) | **60.4%** | **55.8%** | **33.1%** | **43.6%** |

Board-unit projection (233 §3.3's arithmetic; **a projection, not a board prediction** — it assumes
each board member's scaled credit scales like the local credit above the deranged null, and holds
`fid` fixed because our board `fid` is negative):

| | `pds` | `mse` | `nmae` | `fid` | `reach` | `jac` | **overall** | **of today** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| today (`h1union`, board) | 0.7454 | 0.1796 | 0.1107 | −0.0165 | 0.2229 | 0.0288 | **0.21183** | 100% |
| `real` | 0.4805 | 0.1561 | 0.0485 | −0.0165 | 0.1379 | 0.0151 | **0.13694** | **64.6%** |
| `pess` | 0.4065 | 0.1274 | 0.0549 | −0.0165 | 0.1183 | 0.0149 | **0.11759** | **55.5%** |
| **`efr`** | 0.4313 | 0.0273 | 0.0451 | −0.0165 | 0.0721 | 0.0089 | **0.09471** | **44.7%** |
| **`efp`** | 0.4198 | 0.0689 | 0.0652 | −0.0165 | 0.0917 | 0.0100 | **0.10652** | **50.3%** |

**233's 65% / 0.137 moves DOWN to 44.7–50.3% / 0.095–0.107 on an essential-free draw**, and the
member that moves is `mse` (86.9% / 70.9% → 15.2% / 38.4%), because `mse*` credit is `cos²` and the
aggregate-cosine retention falls to 37.7% / 59.9%. **Essential genes are the most transferable
stratum, so a panel without them is harder for the whole stack, not only for the October four.**
This cuts against the direction of my earlier framing and is reported as measured.

**But note the ceiling row in §2.5 before planning on it.** These are random draws; Arc *selects*
for cross-line reproducibility and gets 100% of its own ceiling without using essential genes at
all, while a random essential-free draw gets 72.8–87.6%. **233 §6's "lower bound of unknown slack"
applies with more force here, not less**, and the honest planning statement is that October's
retention sits somewhere between ~45% and no material loss.

---

## 4. `kolf2`, reported separately

| | on `abc` | on `real` / `pess` / `efr` / `efp` |
|---|---|---|
| `|panel ∩ coverage(kolf2)|` | **273** | **0 / 0 / 0 / 0** |
| fold evaluable targets | 273 | **0 — the fold does not exist** |
| BASE | `cos_agg` 0.01413, `pds` 0.63616, `nmae` 0.99464, `fid` 0.50618, `reach` 0.02674, `jac89` 0.00552 | — |
| deranged null | `cos_agg` 0.00026, `pds` 0.49524, `nmae` 1.00010, `fid` 0.49974, `reach` 0.01122, `jac89` 0.00532 | — |
| the four October sources | **all four VACUOUS (0 carrying rows)** | — |

`kolf2` is a **panel-fitted** library: **273 of its 273 targets are on today's 300** (225 §8.3).
On every disjoint panel that intersection is empty, so PLOLO's fourth fold has zero evaluable
targets and cannot be scored at all — **observed, not inferred**: `fold233.build('efp','kolf2')`
returns `n_targets 0, panel_coverage "0/300"`. The `abc` row is given above so the fold is on the
record; on `abc` it, too, prices all four October sources at **VACUOUS**.

---

## 5. Falsifiers

| # | falsifier | outcome |
|---|---|---|
| **1** | **operator re-fit per panel, axis asserted** | **7 fits**, every load verified against a recomputed sha256 of `union(donor coverage) − panel`. `n_cols` is **17,942 on four of five panels** and 17,806/17,807 on the two k562-held bases — i.e. it **cannot** identify the panel, which is the trap. Known-answer row `e2f65ecf…` **MATCH under both builder versions** (the shared file changed mid-run; §1). |
| **2** | **carrying-row sizing** | `per_unit` recorded both ways for every arm: the **all-rows** median returns **EXACTLY 0.0 in 56 of 56 arms** across all four synthetic panels (233 §5.1 measured 8; this extends it to every arm and to the essential-free panels), against a carrying-row `per_unit` of **0.97–1.46**. Every dose here is solved on carrying rows; the shipped rule would have delivered identically zero and the whole experiment would have read "no effect" for the wrong reason. |
| **3** | **`k562_essential` never on the `k562_gwps` fold** | enforced in code (`FORBIDDEN`), logged on every occurrence (`[panel/k562_gwps] … dropping ['k562_essential'] (same cell line as the truth)`). No source is ever used on the fold where it is the truth. |
| **4** | **essential fraction of carrying rows on every panel** | §0.3 and §2.1, per source per fold, **28/28 exact** against 236 §1's published carrying-row matrix. |
| **5** *(added)* | **two independent draws; if they disagree, say so** | They **disagree by 30×** on the composite and the file says so rather than quoting the favourable one. The disagreement is resolved by 236's own leave-top-3-out (20.6% vs −33.8% retention) and by naming the rows (§2.3) — not by averaging it away. |
| **6** *(added)* | **within-panel decomposition at matched per-row dose** | §2.4 — removes the "it is a different panel" confound entirely. 10/10 non-essential halves deliver negative `pds`. |
| **7** *(added)* | **positive control on the new panels** | `cd4_off` / `h1_off` / `resid_off` on every (panel, fold). On `efp`/`efr` removing the CD4 residual still costs **+17.74% / +12.33%** relative `cos_agg` and removing both residuals **+27.55% / +16.98%** (`abc` +22.38% / +40.28%; `pess` +10.11% / +24.99%). **The instrument resolves a real source on these panels perfectly well**, so the October four reading ≈0 is a property of **those sources**, not of the panel. (Side note, reported not chased: the H1 union's own contribution falls to **+3.38%** on `efr` against +10.31% on `abc` — 225's warning about `gse295214_h1` 64 → 34 and `h1_2025` 20 → 3, visible again.) |
| **8** | **a measured null, not an assumed zero** | every arm carries its own whole-row derangement at the same dose; every retention ratio is credit above the panel's and fold's own measured null (`pds` null 0.475–0.502, `fid` 0.499–0.500, `nmae` 1.000–1.020, `cos_agg` 0.0004–0.0020). |

---

## 6. What this says, and what it retracts

1. **RETRACTED: "wire the October four; they price at +0.0066…+0.0104."** That number is a
   property of the **30–50 common-essential targets 199's sampling model put in the panel**.
   Arc's one observed panel contains **zero targets from any of the four libraries**, and on a
   draw that matches Arc on that property the robust value is **16%** of 236's with `pds`
   changing sign. 236's arithmetic is correct; the panel it was given was not representative,
   and I wrote that mandate.
2. **`k562_essential` is dead on arrival** — its perturbation list *is* the essential list, so it
   carries exactly zero on any panel selected the way Arc selects. It should be removed from the
   October plan outright, not down-weighted.
3. **Registering the sources is still free and still correct.** Nothing here says un-register them:
   the registry entry costs nothing, the carrying-row sizing rule (falsifier 2) is needed by *any*
   low-coverage source, and if Arc's D/E/F panel does contain essential-adjacent targets the wiring
   is already there. **What must change is the expected value, which is ≈0, and the plan must not
   depend on it.** 236 §6.3 already said the composition buys +0.0006…+0.0015 over `rpe1` alone;
   the correction here is to the base, not to the increment.
4. **Do not re-run this measurement on a donor-coverage draw.** Any panel sampled from donor
   coverage inherits ~10% essential genes, because donor libraries screen them. `DEF-RUNBOOK` §5's
   two bundles should carry this caveat wherever the October sources are priced against them.
5. **The retention number moves the other way and that is the more consequential finding.**
   233's 55.6–64.6% becomes **44.7–50.3%** on an essential-free draw, driven by `mse` (`cos²`).
   Essential targets are **2.75–4.14× more cross-line transferable**, so they are the stratum our
   direction channel is best at — and Arc excludes them. Arc's own panel nevertheless out-transfers
   every random draw (§2.5), which is 233 §6's selection argument holding; **the planning range
   widens rather than narrowing.**
6. **What generalises beyond this experiment.** A source's value must be priced on a panel that
   matches the target panel's **composition**, not merely its **disjointness**. Disjointness was
   the property 199 designed for and it was not enough. The cheap general check is §0.3's
   list-free one: **intersect the candidate source's own perturbation list with the one real panel
   we have.** For these four it is 0 of 300, and that was on disk before any of this ran.

**No slot is proposed. Nothing was built. Nothing outside this directory was written.**

---

## 7. Files

`panels238.py` (five-panel registry, axis assertion, the `kolf2` fold) ·
`bind238.py` (rebinds 233's instrument onto it, so `fold233`/`score233` run verbatim) ·
`census238.py` (§0.1, §0.3, the October-source union) ·
`carry238.py` (per-fold carrying rows and their essential fraction; KA against 236 §1) ·
`bundle238.py` (§1, the two essential-free draws) ·
`cd4_238.py` (CD4 re-extraction, byte-identical KA) ·
`whiten238.py` (7 operator fits, gated and serial) ·
`six238.py` (all six members, every arm; the `--ka` gate; the within-panel split) ·
`price238.py` (228's board-unit routes, raw and null-corrected) ·
`retain238.py` (233's retention and board projection) ·
`strength238.py` (§2.5's confound check) ·
`influence238.py` (§2.3, which rows carry it) ·
`summary238.py` (the joined tables).

Artefacts in `work/`: `census238.json`, `carry238.json`, `bundle238.json`, `cd4_238.json`,
`whiten238.json`, `six238.json` (439 rows), `ka238.json` (109 rows), `price238.json`,
`retain238.json`, `strength238.json`, `influence238.json`, `absgain238.json`, `kolf2_238.json`,
plus every run log and the `/usr/bin/time -l` output. The seven operators (~9 GB) are kept so every
number is re-checkable; `ops/abc__production3.npz` is byte-identical to the shipped
`data/donors/whitener_union_pruned.npz` (`e2f65ecf…`) and exists only because the `kolf2` fold needs
the three-donor basis.

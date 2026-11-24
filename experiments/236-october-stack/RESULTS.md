# 236 — How the four October-only sources COMPOSE

**Verdict, one line:** *the October four compose to **+18.9% (pessimistic) / +12.5% (realistic)**
relative aggregate cosine above their own measured null, **saturating at 1.22–1.30× the best
single** — not 225's 1.84× — and the six-member composite is **+0.0066…+0.0104 (pess) /
+0.0057…+0.0081 (real)** before any `reach` credit, with **five of six members positive on both
panels once the null is subtracted** and `nmae` costing 1–12% of what `mse` earns. **They are
not four sources; they are one asset with four names** — **71–72% of their carrying-row slots
are redundant** (mean pairwise Jaccard **0.78–0.83**), the 4th member is worth **negative**, and
composing all of them buys only **+0.0006…+0.0015** over wiring `rpe1` alone, which is **below
the +0.0019 shape gate on its own**. **Wire them; do not plan around the 1.8× branch.**

**Date 2026-09-04.** Deterministic (seeds 7 / 2361 / 9090+n / 20726 / 23323 / 233233).
**1,044 scored arms** across 6 (panel, fold) cells, **18 axis-asserted operator loads**, 15
subsets × 6 doses × 2 panels, 42 six-member rows, 108 row-matched-control rows. **No `vcc`, no
`cell-eval2`, no build, no submission, nothing written outside
`experiments/236-october-stack/`.** No production file read or touched (235 owns
`scripts/make_submission_pm_cd4.py` and `scripts/build_union_whitener.py`; this experiment
never opened either). Peak resident **8.36 GiB** (`/usr/bin/time -l`: **6.82 / 8.19 / 8.36 GiB**
over 586 / 140 / 72 s for the three heavy jobs, run strictly serially — 13 minutes of compute in
total, because 233's six per-panel operators were reused rather than re-fitted). Memory
gate at every launch:
`memory_pressure` 88% free, `vm_stat` free+inactive 32–40 GiB, swap free 1,679 MiB.

**Nothing here is submittable.** The panels do not exist until 2026-10-22. This sizes an asset.

---

## 0. Phase 0 — known-answer rows. **Four gates, four passes, all at exactly zero.**

| # | gate | expected | measured | verdict |
|---|---|---|---|---|
| **KA-1** | every arm 233 also ran must reproduce it through this file | 157 rows | **157/157**, max &#124;Δ`cos_agg`&#124; **0.0e+00**, max &#124;Δ`pds`&#124; **0.0e+00** | **PASS** |
| **KA-2** | BASE — "the shipped stack as it will exist in October" — must be 233 §3.1's object on **all six** members | 36 values | **36/36 bit-identical** (`cos_agg`, `pds`, `nmae`, `fid`, `reach`, `jac89` × 6 panel-folds) | **PASS** |
| **KA-3** | the coverage census, recomputed from raw `perts.txt` + `pert_counts.csv` | 233's `census233.json` | **45/45 matched**: `k562_gwps` **269 → 150**, Orion **291 / 293**, `kolf2` 273 → 0, `feng_ipsc` 182 → 0, the four October sources **0 → 30–33** (pess) / **34–48** (real) | **PASS** |
| **KA-4** | panel disjointness, asserted not assumed | 0 | **&#124;abc ∩ real&#124; = 0**, **&#124;abc ∩ pess&#124; = 0** (&#124;real ∩ pess&#124; = 6, irrelevant) | **PASS** |

**233's published single-source table is reproduced to the last digit**, on the same
`DEF-RUNBOOK` §5 panels, the same per-panel re-fitted operators, the same per-panel re-extracted
CD4 banks, the same folds and the same read-back:

| panel | `rpe1` | `nadig_jurkat` | `nadig_hepg2` | `k562_essential` |
|---|---:|---:|---:|---:|
| pessimistic | **+17.66%** | **+17.27%** | **+14.82%** | **+11.71%** |
| realistic | **+10.49%** | **+6.44%** | **+5.51%** | **+8.85%** |

**One row 233 never ran, now filled in.** 233 priced `real` on two folds only; this run adds the
`real / k562_gwps` fold. It moves the realistic singles to **+11.94 / +9.59 / +5.31%** on the
full three-fold set (`k562_essential` stays 2-fold — see the leak rule).

---

## 1. 🔴 The mechanism, before the curve: **they carry the SAME ROWS**

This is the finding everything else follows from, and it is visible without any scoring.

| panel / fold | fold targets | `nadig_hepg2` | `nadig_jurkat` | `rpe1` | `k562_essential` | **UNION** |
|---|---:|---:|---:|---:|---:|---:|
| pess / `orion_hct116` | 291 | 28 | 28 | 27 | 25 | **30** |
| pess / `orion_hek293t` | 293 | 29 | 29 | 28 | 26 | **31** |
| pess / `k562_gwps` | 150 | 32 | 32 | 31 | *(leak)* | **34** |
| real / `orion_hct116` | 300 | 45 | 48 | 47 | 34 | **50** |
| real / `orion_hek293t` | 300 | 45 | 48 | 47 | 34 | **50** |
| real / `k562_gwps` | 300 | 45 | 48 | 47 | *(leak)* | **50** |

**108 row-slots collapse to 30 distinct rows** on the pessimistic draw; **174 collapse to 50** on
the realistic one — **71–72% redundancy** on the four-source folds, 64% on the three-source K562
folds. Mean pairwise Jaccard of the carrying-row sets:

| panel / fold | slots | union | redundancy | pairwise Jaccard (min–max, mean) |
|---|---:|---:|---:|---|
| pess / `orion_hct116` | 108 | 30 | **72.2%** | 0.733 – 0.931 (**0.821**) |
| pess / `orion_hek293t` | 112 | 31 | **72.3%** | 0.742 – 0.933 (**0.827**) |
| pess / `k562_gwps` | 95 | 34 | 64.2% | 0.853 – 0.939 (0.882) |
| real / `orion_hct116` | 174 | 50 | **71.3%** | 0.640 – 0.939 (**0.775**) |
| real / `orion_hek293t` | 174 | 50 | **71.3%** | 0.640 – 0.939 (**0.775**) |
| real / `k562_gwps` | 140 | 50 | 64.3% | 0.840 – 0.939 (0.892) |

The tightest pair is `nadig_hepg2 ∩ nadig_jurkat` — **27 of 28** (pess), **44 of 45** (real),
J = 0.93–0.94, same lab and same library — and it is exactly the pair the leave-one-out in §2.3
finds redundant.

**So this is not a coverage stack, it is four ESTIMATES of the same ~30–50 targets.** All four are
common-essential libraries (Nadig's two lines, Replogle's `k562_essential`, and `rpe1`, whose
panel-overlapping rows are the essential ones). The composition available is 159's `1/k`
averaging, not 225's coverage addition — and that is exactly what the curve below measures.

---

## 2. Phase 1 — THE STACK CURVE

Form: `compose218.union_rows` — the H1-union precedent already shipped in production (222) and
225's own stack form — as an un-whitened post-whitening residual on top of the shipped October
stack, **sized on carrying rows** (the all-rows median is `0.0000` for all four; the shipped rule
delivers identically zero — 233 §5.1). All **15** non-empty subsets were measured at **6** doses
on **both** panels, which yields the greedy order, all 24 orders and every leave-one-out from one
pass.

Fold set **CLEAN2** = `orion_hct116` + `orion_hek293t`, the two folds on which all 15 subsets are
defined, so the k-curve is apples-to-apples. Dose **g = 0.34**, 233's pre-registered dose (§2.4).

### 2.1 The curve, greedy order, both panels

| k | added | stack | **pess** rel `cos` | **sat ×** | Δ`pds` | shuf | **real** rel `cos` | **sat ×** | Δ`pds` | shuf |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | `rpe1` | `rpe1` | **+15.57%** | **1.00×** | −0.00095 | −0.84% | **+10.49%** | **1.00×** | +0.00271 | −0.67% |
| 2 | `k562_essential` | +`k562_essential` | **+19.25%** | **1.24×** | +0.00132 | −0.63% | **+12.19%** | **1.16×** | +0.00230 | −0.72% |
| 3 | `nadig_jurkat` | +`nadig_jurkat` | **+19.94%** | **1.28×** | +0.00236 | +0.10% | **+12.81%** | **1.22×** | +0.00190 | −0.01% |
| 4 | `nadig_hepg2` | all four | **+19.83%** | **1.27×** | +0.00168 | +0.95% | **+12.77%** | **1.22×** | +0.00306 | +0.23% |

The greedy order is the **same on both panels**. The best subset of each size is the greedy one at
every k on both panels.

**On all three folds** (the subsets without `k562_essential`, the only ones the K562 fold admits):

| k | stack | pess rel `cos` | sat × | real rel `cos` | sat × |
|---|---|---:|---:|---:|---:|
| 1 | `rpe1` | +17.66% | 1.00× | +11.94% | 1.00× |
| 2 | +`nadig_jurkat` | **+22.19%** | **1.26×** | **+14.33%** | **1.20×** |
| 3 | +`nadig_hepg2` | **+23.04%** | **1.30×** | +12.85% | 1.08× |

### 2.2 🔴 The saturation ratio, and the hypothesis it kills

```
SATURATION AT FOUR SOURCES:   1.27x  (pessimistic)     1.22x  (realistic)
                              1.30x / 1.08x at three, on the full three-fold set
225 on TODAY's panel:         1.65x at four,  1.84x at thirteen
```

**The mandate's decisive branch — "if the October four saturate like 225's sources, the stack is
~1.8× the best single (~+30% relative cosine, transformative)" — is FALSIFIED.** They saturate
**harder** than today's sources, at 1.22–1.30× rather than 1.65–1.84×, and §1 is why: 225's stack
added *rows*, this one adds *replicates of the same rows*.

**The marginal source, in one column** (pessimistic, CLEAN2, points of relative cosine):

| step | source added | marginal | fraction of step 1 |
|---|---|---:|---:|
| 1 | `rpe1` | **+15.57** | 100% |
| 2 | `k562_essential` | +3.68 | 24% |
| 3 | `nadig_jurkat` | +0.69 | 4% |
| 4 | `nadig_hepg2` | **−0.11** | **−1%** |

### 2.3 The fourth source is worth **negative**, and leave-one-out says which three matter

Leave-one-out of the 4-stack @ 0.34, CLEAN2 (cost in points of relative cosine):

| removed | pess | real |
|---|---:|---:|
| `rpe1` | **+2.05** | +1.53 |
| `nadig_jurkat` | +1.67 | +0.73 |
| `k562_essential` | +0.87 | **+1.71** |
| `nadig_hepg2` | **−0.11** | **−0.04** |

**Removing `nadig_hepg2` IMPROVES the stack on both panels.** It is redundant with
`nadig_jurkat` — same lab, same essential library, carrying rows overlapping 27/28 and 44/45.
**The asset is three sources.**

### 2.4 Order matters early and not at the end; and the dose is NOT identified

All 24 orders of the four, @0.34, CLEAN2, spread in points of relative cosine:

| after k sources | pess min / median / max | spread | real min / median / max | spread |
|---|---|---:|---|---:|
| 1 | +11.10 / +12.95 / +15.57 | **4.47** | +5.51 / +7.65 / +10.49 | **4.97** |
| 2 | +15.02 / +17.23 / +19.25 | 4.23 | +7.29 / +10.72 / +12.19 | 4.90 |
| 3 | +17.78 / +18.56 / +19.94 | 2.16 | +11.06 / +11.64 / +12.81 | 1.75 |
| 4 | +19.83 / +19.83 / +19.83 | **0.00** | +12.77 / +12.77 / +12.77 | **0.00** |

**Order costs up to 5 points if you stop early and exactly nothing if you take all four** (the
union is one object, so k=4 is order-free by construction). Picking wrong at k=1 costs 4.5–5.0
points — more than steps 2–4 combined buy.

⚠️ **The dose has no interior optimum on this grid.** The 4-stack is monotone through the last
dose tested on both panels (pess +3.01 → +6.00 → +11.86 → +19.83 → +25.86 → **+31.16%** at
g = 0.05…0.55; real +2.00 → … → **+19.57%**), so the argmax is the grid edge. **The headline is
quoted at 233's pre-registered g = 0.34 and NOT at the argmax** — choosing the dose on the
evaluation folds is precisely the free choice §8.1d faults 228 for. At the argmax the saturation
ratios are unchanged (**1.29× / 1.24×**), which is the useful part: *the composition ratio is
dose-invariant even though the level is not.*

---

## 3. Phase 2 — the six members, and `e = ΔM/ΔU`

228's arm survived read-back on cosine and `pds` and died at **−0.00702** on the composite because
its cosine was bought with MAGNITUDES (`e` = **0.072** against break-even 0.5). AGENTS.md §8.1d
therefore requires the question be answered before the arm is priced. 228's exact routes and
constants (`price228.py`) are reused unchanged.

### 3.1 🔴 The raw deltas disagreed between panels on `nmae`. **The derangement explains why.**

| CLEAN2 @0.34, 4-stack | pess | real |
|---|---:|---:|
| Δ`nmae` **raw** | **+0.00260** (worse) | **−0.00197** (better) |
| Δ`nmae` of its **own derangement** | −0.00037 | **−0.00227** |
| Δ`nmae` **above the null** | **+0.00297** | **+0.00030** |
| `e` raw | 0.229 | **3.211** |
| **`e` above the null** | **0.180** | **0.433** |

**The realistic panel's entire `nmae` credit is reproduced by a whole-row derangement of the same
residual at the same dose** — i.e. it is a dose artefact of adding residual mass, not a property of
the source. Once the null is subtracted (233 §3.2's own methodology: every number is credit above a
*measured* null) **the two panels agree in direction on all six members.** Falsifier 3 is
satisfied by correcting the estimator, not by asserting agreement.

### 3.2 The null-corrected six-member price

Arm minus its **own** derangement at the same dose, CLEAN2, g = 0.34. `k_cos` is carried as a
range (1.000 on `h1union`, 1.384 on `cd4f034`) because 219 showed it is not identified.

| arm | panel | rel `cos` | Δ`pds` | Δ`nmae` | Δ`fid` | **`e`** | **composite (no reach)** | with reach |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| **all four** | pess | **+18.88%** | +0.00114 | +0.00297 | +0.00036 | **0.180** | **+0.00664 … +0.00926** | +0.00994 … +0.01256 |
| **jurkat+rpe1+k562ess** | pess | **+19.84%** | +0.00254 | +0.00286 | +0.00049 | **0.206** | **+0.00763 … +0.01039** | +0.01142 … +0.01417 |
| `rpe1` alone | pess | +16.40% | +0.00309 | +0.00259 | +0.00041 | 0.190 | +0.00661 … +0.00889 | +0.01035 … +0.01263 |
| **all four** | real | **+12.54%** | +0.00415 | +0.00030 | +0.00053 | **0.433** | **+0.00637 … +0.00811** | +0.00902 … +0.01076 |
| **jurkat+rpe1+k562ess** | real | +12.82% | +0.00227 | +0.00045 | +0.00049 | 0.419 | +0.00568 … +0.00746 | +0.00795 … +0.00973 |
| `rpe1` alone | real | +11.15% | +0.00467 | +0.00059 | +0.00032 | 0.403 | +0.00580 … +0.00735 | +0.00772 … +0.00927 |

Board-unit member split for the pessimistic 4-stack: `pds` +0.00043, `mse` **+0.00683…+0.00946**,
`nmae` **−0.00083**, `fid` +0.00020, `jac` +0.00001, [`reach` +0.00330, rule-9 indicative].
Against the gates: **3.5–5.5× the +0.0019 shape gate**, **0.75–1.17× the +0.0089 amplitude gate**
without `reach`, and clear of both with it.

### 3.3 SIGNS or MAGNITUDES: **both, and the magnitude is 52–340× cheaper than 228's**

`e` is **0.18–0.43 on both panels — below the 0.5 break-even**, so this is not a clean sign arm
like `h1union`. But `e` alone is the wrong summary when the *quantity* of magnitude moved differs
by an order of magnitude. The decisive ratio is the `nmae` toll against the `mse` credit it buys:

| arm | rel `cos` | Δ`nmae` | `e` | **`nmae` toll ÷ `mse` credit** | composite |
|---|---:|---:|---:|---:|---:|
| **228's per-gene magnitude gain** | +4.55% | **+0.0373** | 0.072 | **4.6 … 6.3×** (the toll swamps the credit) | **−0.00702** |
| **236, pessimistic 4-stack** | +18.88% | +0.00297 | 0.180 | **0.088 … 0.121** | +0.0066 … +0.0093 |
| **236, realistic 4-stack** | +12.54% | +0.00030 | 0.433 | **0.013 … 0.018** | +0.0064 … +0.0081 |

(228's row is its published raw delta; both 236 rows are null-corrected, which is the
conservative direction — the raw pessimistic Δ`nmae` is +0.00260 and the raw realistic one is
**−0.00197**, i.e. the raw numbers flatter the arm on one panel and the correction removes that.)

**`nmae` charges 1–12% of the `mse` credit here, against 460–630% for 228.** Per point of relative
cosine the toll is **52× (pess) to 340× (real) smaller** than 228's.

**So the honest answer to §8.1d's question is: neither cleanly.** This is **not** a pure sign arm
— `e` is below break-even on both panels, so `nmae` does charge for it. It is also nothing like
228's magnitude arm: it moves **8.4× less** magnitude onto the reference's DE pool while buying
**3–4× more** cosine, so the toll is a rounding error against the credit instead of five times
it. **Five of the six members contribute positively on both panels** (`pds`, `mse`, `fid`, `jac`,
and indicatively `reach`); **`nmae` is the one that costs, at −0.00083 (pess) / −0.00008 (real)
board units against an `mse` leg of +0.0045…+0.0095.** No member turns the arm negative on either
panel, which is precisely what failed for 228 and for 229's `ARM6`.

---

## 4. Falsifiers — the five mandated, plus one this stack needed

| # | falsifier | outcome |
|---|---|---|
| **1** | **whole-row derangement** of each source/union at matched dose | run at **every** dose for **every** subset (**444** derangement arms). 4-stack @0.34: **+0.95% (pess) / +0.23% (real)** against **+19.83% / +12.77%** real — **21× and 56×**. Within ±1.1% at every dose on both panels. **And it earned its keep**: it caught the realistic panel's `nmae` credit as a dose artefact (§3.1) |
| **2** | `k562_essential` never on the fold holding `k562_gwps` | enforced in code (`FORBIDDEN`); every subset containing it is absent on that fold and reported as such, never silently averaged. No source is ever used on the fold where it is the truth |
| **3** | **both panels**; if they disagree in direction, conclude nothing | raw deltas disagreed on **`nmae` only**; the derangement identified the cause and the null-corrected deltas **agree in direction on all six members**. The CLEAN2 saturation ratio agrees to **0.05×** (1.27 vs 1.22) and the greedy order is **identical** |
| **4** | **leave-top-3-out retention**, `cos` and `pds` **separately** | **`cos` SURVIVES**: 30.6–73.2% of the gain retained and **still positive in 18/18 rows** (mean ≈ 55%). **`pds` DOES NOT**: after dropping the top 3 the `pds` gain is **negative in 13 of 18 rows**. **Quote the cosine channel; do not quote a per-fold `pds` gain from these sources** |
| **5** | **operator re-fit per panel, axis asserted** | **18/18** loads verified against a recomputed `sha256` of `union(donor coverage) − panel`; `n_cols` is **17,942 on every panel**, which is exactly why the production check cannot see it (233 §2: 98.8% mis-pairing with every production check passing, priced at −0.021 on `pds` alone) |
| **6** *(added)* | **row-matched control** — 218's, which 233 never ran for these sources | see §5 |

### 4.1 Falsifier 4, in full

Influence = that row's own leave-one-out effect on the **gain**; for `pds` the row also leaves the
**decoy set**, which is the honest form for a retrieval statistic. Retention is quoted as the
ratio of the context-mean gains (a mean of per-context ratios explodes where a denominator is
near zero, so it is not used). The `pds` column is given as the **absolute gain after the drop**,
because a ratio of two small signed numbers is unreadable.

| panel | fold | `cos` gain → after | **retained** | `pds` gain → **after** |
|---|---|---|---:|---|
| pess | `k562_gwps` | +0.01205 → +0.00776 | **64.4%** | +0.01317 → **+0.00568** |
| pess | `orion_hct116` | +0.00406 → +0.00179 | **44.1%** | +0.00154 → **−0.00097** |
| pess | `orion_hek293t` | +0.00352 → +0.00116 | **33.1%** | +0.00181 → **−0.00173** |
| real | `k562_gwps` | +0.00641 → +0.00463 | **72.2%** | −0.00116 → **−0.00281** |
| real | `orion_hct116` | +0.00380 → +0.00243 | **63.9%** | +0.00549 → **+0.00290** |
| real | `orion_hek293t` | +0.00219 → +0.00137 | **62.9%** | +0.00064 → **−0.00223** |

(4-stack rows — on the `k562_gwps` fold that is the three-source stack, `k562_essential` being
banned there; the singles and 3-stacks are in `work/six236.json` and follow the same pattern.)
**The top 3 of 150–300 targets carry ~45% of the cosine gain** — concentrated, as 25–50 carrying
rows must be, but the cosine gain **never inverts, in 18 of 18 rows**. The `pds` gain **does**,
in **13 of 18**. **The `pds` channel is not robust at fold level** and must not be quoted per
fold; the CLEAN2-mean Δ`pds` in §3.2 is positive on both panels and is the only form worth
carrying.

---

## 5. 🔴 The control that decides whether this is biology or bookkeeping

The four sources carry the *same* ~30–50 rows (§1), so the entire gain could be **more residual
mass on a favourable row subset** — a dose change we already control, dressed as a new source.
218 built the control that separates them: **CD4's own effect, restricted to the stack's own rows,
at the stack's own dose, with the same carrying-row sizing.** 233 never ran it for these sources.

**Realistic panel — the control is EXACTLY row-matched (50 rows vs 50):**

| dose | stack | **rowcd4** (same rows, CD4's effect) | randcd4 | Δ`pds` stack / rowcd4 |
|---:|---:|---:|---:|---|
| 0.05 | +2.00% | +0.72% | −0.10% | +0.00056 / −0.00027 |
| 0.20 | +7.76% | +2.65% | −0.69% | +0.00184 / −0.00140 |
| **0.34** | **+12.77%** | **+4.16%** | −1.58% | **+0.00306 / −0.00254** |
| 0.55 | +19.57% | +5.92% | −3.45% | +0.00442 / −0.00395 |

Pessimistic panel @0.34: stack **+19.83%**, rowcd4 **+2.54%**, randcd4 +0.44% (rowcd4 covers only
20 of the stack's 30 rows there, because CD4's own union collapses to 174/300 — so the pessimistic
control understates itself and the realistic row is the clean comparison).

**The stack is 3.1× (realistic, exactly row-matched) to 7.8× (pessimistic) its row-matched
control on cosine, and OPPOSITE IN SIGN on `pds`** — `rowcd4` is negative at every dose on both
panels while the stack is positive at every dose. **It is the biology, not the rows.**

---

## 6. What this says about October

1. **Wire them, in the union form, sized on carrying rows** — `rpe1` + `k562_essential` +
   `nadig_jurkat` are the three that pay; `nadig_hepg2` is optional (§2.3). The asset is real:
   it survives its own derangement by **21–56×**, survives a row-matched control by **3–8×**, and
   prices at **+0.0066…+0.0104 (pess) / +0.0057…+0.0081 (real)** on the six-member composite
   before any `reach` credit. Against 233 §3.3's own October projection of **0.1178 (pess) /
   0.1368 (real)**, that is **+5.6…+8.8% (pess) / +4.2…+5.9% (real) relative** — real, and not
   transformative.
2. **The 1.8× branch is dead; plan on 1.25×.** Anyone budgeting the October four as "four sources
   at 225's saturation curve" is over-budgeting by **~45%**. The reason is structural and will not
   change with a better dose or a better composition rule: **they carry the same rows.**
3. **The first source is 85–90% of the asset.** Composing all four buys **+0.0006…+0.0015** over
   `rpe1` alone — **below the +0.0019 shape gate by itself**, and one third of the ordering
   noise. **Do not treat "add the other three" as an arm.** Take them because they are free once
   the registry entry exists, not because they clear a gate.
4. **Pick `rpe1` first and it costs nothing to get the order wrong afterwards.** k=1 choice is
   worth 4.5–5.0 points of relative cosine — more than steps 2–4 combined.
5. **Drop `nadig_hepg2`, or keep it knowing it is worth −0.1 points.** It is `nadig_jurkat`'s
   twin (27/28 and 44/45 shared rows, same lab, same library).
6. **`nmae` does not close this arm**, which is the question §8.1d demands be asked first. `e` is
   0.18–0.43, below break-even, so it is not a pure sign arm — but the toll is **1–12% of the
   `mse` credit** against 228's 460–630%, because the *quantity* of magnitude moved is 8.4× smaller
   while the cosine bought is 3–4× larger.
7. **Quote the cosine channel, not a per-fold `pds` gain.** Falsifier 4 killed the latter.
8. **The dose is not identified here** and must not be fitted on these folds. 233's 0.34 is the
   honest operating point; the saturation ratio is dose-invariant even though the level is not.

**And the standing caveat that outranks all of the above: 233 §6.** These are **random** draws;
Arc **selects** its panel for cross-line reproducibility and "adequate representation in existing
public datasets". Every number here is a **lower bound of unknown slack** — and note the direction
of the bias for *this* question specifically: a selected panel will have **better** public
coverage, which raises these sources' carrying-row counts and could either raise the level (more
rows) or deepen the saturation (more shared rows). **Which of those it does is not established.**

---

## 7. Files

`stack236.py` (Phase 0 + the 15-subset × 6-dose curve, 894 arms, KA-1) ·
`census236.py` (KA-3, KA-4, independent of 233) · `kabase236.py` (KA-2, six members) ·
`six236.py` (Phase 2: six members, 188's `e = ΔM/ΔU`, falsifier 4) ·
`rowctl236.py` (§5, 218's row-matched control) · `price236.py` (228's exact routes, raw and
null-corrected) · `analyse236.py` (the curve, order sweep, LOO, overlap).
Artefacts in `work/`: `stack236.json`, `six236.json`, `rowctl236.json`, `census236.json`,
`kabase236.json`, `analyse236.json`, `price236.json`, and the three `/usr/bin/time -l` logs.
No operator was fitted — the six 233 per-panel operators this experiment needs (`pess` and
`real` × three cancer folds) were reused and axis-verified 18 times, which is why the whole
experiment cost 13 minutes of compute instead of an hour.

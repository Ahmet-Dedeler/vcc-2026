# 250 — Does the CD4 deflation TRANSFER to a disjoint panel, and where is its k knee there?

**Verdict, one line:** *the deflation's `pds` mechanism transfers to a disjoint 300-target panel
— **`pds` +0.0148 / +0.0099 (3/3 cancer folds, 9.9σ / 6.6σ) on the two full-coverage draws
(`real`, `efr`), controls clean** — but **its cosine cost is 4–14× larger off A/B/C** (−2.1 / −6.6 /
−5.2 / −3.5 % aggregate on `real` / `pess` / `efp` / `efr` against −0.47 % on A/B/C), so the day-of
re-fit at k = 5 clears the +0.0019 gate only where CD4 covers the whole panel (`real` +0.0048; `efr`
+0.0019, at the gate) and reads **≤ 0 on both low-coverage draws** (`pess` −0.0002, `efp` −0.0026);
**k = 5 sits BELOW the `pds` plateau on every disjoint panel** (k-max at 10–15; k = 5 delivers 79 / 61 /
47 / 52 % of it) though the six-member price is flat or falling past k = 5–10. **The projector is NOT the
same subspace across panels** (first principal angle 24–27°, fifth 88°; A/B/C's `V_5` captures 0.10–0.13
of a disjoint panel's centred energy against 0.20 for its own) — **and yet the PINNED A/B/C production
basis (`U_5` sha `f06e0262…`, the object 249 built and armed) beats the day-of re-fit on 4/4 disjoint
panels**, `pds` +0.0153 / +0.0037 / +0.0073 / +0.0107 at only −0.1 … −2.0 % cosine, price
+0.0061 / +0.0009 / +0.0014 / +0.0030, never negative; any foreign panel's basis does the same, so a
basis fitted on OTHER rows removes the shared T-cell hub and nothing else, while a basis fitted on the
rows it is applied to also removes panel-specific signal. **October should pin, not re-fit.***

**Date 2026-09-04.** Deterministic (246's seeds `20260904`(+offsets) / `20726`, 233's `233233`).
**No `vcc`, no `cell-eval2`, no build, no submission, no slot proposal, no production change; nothing
written outside `experiments/250-deflation-october/`** (every fold cache, operator and CD4 bank is read
from 233's / 238's / 246's / 249's directories; `plolo.WORK` is pinned to 250's `work/plolo/`; `git
status` and an mtime scan of every directory read confirm no writes). **Strictly serial**, memory gate
(`memory_pressure` ≥ 40 % free AND free+inactive ≥ 20 GiB AND swap free ≥ 400 MiB) logged before every
fold load (readings 74–77 % / 26–35 GiB / 918–1645 MiB); `/usr/bin/time -l` on all 36 heavy runs;
**peak RSS 4.12–5.11 GiB** (sweeps), 4.82–5.07 (bridges), 4.20–4.68 (cross), 1.03 (basis). 16 sweep
runs × 29–30 arms, 4 bridge runs × 22 arms, 15 cross runs × 11 arms, one basis run: **713 scored
arm×fold×panel rows**, ~17 min wall clock in total.

---

## 0. The answer in twelve lines

1. **Every known-answer row passes at 0.0e+00** (§1): 233's `BASE` and deranged `null` on all 15
   (panel, fold) cells through 233's own `build_final`; 250's `compose` reproduces that matrix bitwise;
   **40/40 bridge rows reproduce 246's json** (STACK, `B2_cd4_svd{3,5,10}`, the energy-matched random,
   the shuffled projector, the split-half projector, `C2_partial`, `E_cd4_svd5+partial`) on all four
   PLOLO folds; the A/B/C production residual equals 249's staged `cd4-off.npz` bitwise and its k = 5
   deflation reproduces the shipped audit including **`U_5` sha `f06e0262…`**; 233's `real`/`pess`
   retention rows and `analyse233.json`'s thinning rows re-derive at 0.0e+00.
2. **The two instruments differ only by the operator recipe** (§1.3): 246's `B2_cd4_svd5_resolved`
   (062 held operators) reads macro −0.809 % / agg −0.462 % / Δ`pds` +0.012263; through 233's
   essential-pruned re-fit operator (October's recipe) the same objects read **−0.486 % / −0.474 % /
   +0.012412**; 233's float64 consensus vs plolo's float32 storage moves nothing above 4e−11.
3. **Phase 1, the k ladder** (§2): on every disjoint panel Δ`pds` keeps rising to k = 10–15 while
   cosine keeps falling; **k = 5 is not inside the 10 % plateau anywhere** (`abc` 0.898, `real` 0.793,
   `pess` 0.607, `efp` 0.471, `efr` 0.520 of the k-max Δ`pds`). The price-optimal re-fit k is 5 on
   `abc`, 10–15 on `real`, 10 on `efr`, and no k pays on `pess` / `efp`.
4. **Phase 2, controls at k = 5** (§3): on `real` and `efr` every control 246 ran passes (energy-matched
   random −0.0047 … −0.0057 / −0.0020 … −0.0049, shuffled projector +0.0005 / −0.0010, split-half
   retains 100 % / 85–90 %, LT3O 62 % / 55 %, non-carrying rows exact to 1e−16); on `pess` a random
   29-dim subspace (seed 0) buys **+0.0058, more than the arm's +0.0045**, and LT3O is **−59 %**; on
   `efp` one fold is −0.0074 and LT3O **−70 %**. The two low-coverage draws are **DOES NOT**.
5. **Phase 3, identity** (§4): principal angles A/B/C ↔ disjoint 24–27°, 34–40°, 42–61°, 57–78°, 88–90°;
   cross-captured energy 0.103–0.133 vs own 0.198–0.232 (random-5 floor 0.0005). Not the same
   subspace — one to two shared directions plus panel-specific ones.
6. **Yet pinning wins** (§4.3–4.4): the pinned A/B/C basis costs 2–5 pp LESS cosine than the re-fit on
   every disjoint panel at equal or better `pds`, and pinning **any** other panel's basis does the same
   (`efp`: own re-fit +0.0035 at −5.2 % vs foreign +0.0073 … +0.0083 at −0.2 … −0.5 %). The extra
   energy a same-row fit removes (0.20 vs 0.10–0.14) is panel-specific signal, not hub.
7. **Phase 4** (§5): the `+partial` rider is +0.0010 / +0.0009 / +0.0002 / −0.0020 `pds` on the four
   disjoint panels (its cosine +0.5 … +1.1 %), shuffled ≈ 0 — not build-worthy in October either.
8. **`kolf2` apart** (§2.6): 0/273 of its targets sit on any disjoint panel (asserted), so it is
   evaluable on `abc` only: svd5 `pds` +0.00414 at −3.47 % (= 246 to 0.0e+00).
9. **Δq_w does not move off A/B/C** (+0.0003 / −0.0007 / −0.0001 / −0.0005 vs +0.0015): on a disjoint
   panel this is a pure retrieval arm with a cosine bill.
10. **Verdicts** (§6): `real` TRANSFERS, `efr` TRANSFERS (pays at k ≥ 7 re-fit, or pinned k = 5),
    `pess` DOES NOT, `efp` DOES NOT. Coverage separates them (CD4 carried 300 / 300 vs 174 / 186).
11. **The census falsifier is coverage, not the top-5 energy fraction** (§6.2): the fraction reads
    0.203 / 0.199 / 0.232 / 0.198 on the four panels — inside the A/B/C range on `pess`, which fails,
    and outside it on `efp`, which also fails. CD4 carried rows < ~250/300 is the reading that fired on
    both failures; DEF-RUNBOOK §5's pessimistic draw (174) fires it.
12. **Runbook row** (§7): unset until the 2026-09-05 board point; then ship on D/E/F as the **pinned
    A/B/C basis, k = 5** (needs a builder code path and a rehearsal before Oct 15); until that exists,
    re-fit k = 5 only at CD4 carried ≥ 290/300, else unset.

---

## 1. Known-answer rows — 0 failures

### 1.1 The instrument (233's, verbatim), on every panel

`bind250.py` rebinds `panels233` to a five-panel registry (233's artefacts for `abc`/`real`/`pess`,
238's for `efp`/`efr`) with every cache path under 250's `work/`; `fold233.consensus`,
`fold233.verify_op` (the axis assertion), `score233.cd4_prepared`, `score233.build_final`,
`compose218`'s sizing/merge and `E062.add_post_whiten_residual` are called as they are.

| row | reference | cells | result |
|---|---|---|---|
| `compose` (250) == `build_final` (233), same operator, same objects | bitwise (`np.array_equal`) | 16/16 (panel, fold) | **PASS**, 16/16 |
| `BASE233` (all-rows 0.34 + H1 union 0.28 carrying) | `score233.json` BASE (`abc`/`real`/`pess`); `six238.json` BASE (`efp`/`efr`) | 15 | **max\|Δ\| 0.0e+00** on `cos_agg`, `cos`, `pds`, `sgn_w`, 15/15 |
| `null233` (rows deranged, seed 233233) | the same files' `null` | 15 | **0.0e+00**, 15/15 |
| `abc`/`kolf2` | no 233 reference (233 dropped the fold); covered by the bridge | 1 | n/a — bridge row below |
| 233's retention ratios (credit above the null, `real`/`pess` ÷ `abc`), `cos_agg` and `pds` per fold | `retain233.json` | 12 | **0.0e+00** |
| `analyse233.json` thinning rows `orion_hct116|pess` (72.3715 / 55.1277), `orion_hek293t|pess` (73.2305 / 57.6296) | `analyse233.json` | 4 | **0.0e+00** |

### 1.2 The bridge — 250's copies of 246's constructions are bit-faithful (`bridge250.py`)

246's instrument exactly (062's held operators, 246's own fold caches, the shipped CD4 banks, today's
panel; production operator on `kolf2`), driven through `lib250`'s functions:

| arm (246 name) | k562_gwps | orion_hct116 | orion_hek293t | kolf2 |
|---|---|---|---|---|
| `STACK`, `B2_cd4_svd5`, `B2_cd4_svd5_resolved`, `B2_cd4_svd3_resolved`, `B2_cd4_svd10_resolved`, `B2_cd4_randE5`, `B2_cd4_svd5_SHUFDEF`, `B3_cd4_svd5_splithalf`, `C2_partial@0.34`, `E_cd4_svd5+partial@0.34` | 10/10 at **0.0e+00** | 10/10 at 0.0e+00 | 10/10 at 0.0e+00 | 10/10 at 0.0e+00 |

(`max|Δ|` over `cos_agg`, macro, `pds`, `a15_sgn_w`; tolerance 1e−9.) So the deflation, the
energy-matched random subspace, the shuffled projector, the split-half projector, the partial rows and
the merged composition are 246's to the bit before any of them touches a disjoint panel.

### 1.3 The A/B/C row the mandate asked for — and the difference between the two instruments

246's `B2_cd4_svd5_resolved` three-cancer-fold numbers through three constructions (`bridge250` +
`sweep250 --panel abc`):

| construction | base `cos_agg` k562 / hct116 / hek293t | rel agg | rel macro | Δ`pds` | Δq_w |
|---|---|---:|---:|---:|---:|
| **246 verbatim** (062 held operators, plolo float32 consensus) | 0.038339 / 0.031236 / 0.024323 | **−0.4617 %** | **−0.8087 %** | **+0.012263** | +0.00064 |
| the same objects through **233's re-fit A/B/C operator** (essential prune + 12,000-row cap = the recipe every disjoint panel's operator is fitted by) | 0.036983 / 0.029427 / 0.023895 | −0.4740 % | −0.4856 % | +0.012412 | +0.00152 |
| **233's instrument** (re-fit operator, `fold233` float64 consensus) — the `abc` sweep | same | −0.4740 % | −0.4856 % | +0.012412 | +0.00152 |

Per fold, 246 vs re-fit: k562 −1.63 / −1.78 % / +0.009987 vs −1.69 / −1.72 % / +0.009049; hct116
−1.95 / −1.91 % / +0.009446 vs −2.03 / −1.58 % / +0.010623; hek293t +2.20 / +1.27 % / +0.017355 vs
+2.29 / +1.84 % / +0.017564; `kolf2` identical (238's re-fit of the production operator reproduces
`e2f65ecf…`). **The two instruments differ only by the operator recipe**: base cosine −3.5 / −5.8 /
−1.8 % (233's KA-E), the deflation's Δ`pds` +0.00015 and its macro cosine cost halved (−0.81 → −0.49 %);
row-3 vs row-2 (float32 vs float64 consensus) differ by ≤ 4e−11 on every endpoint and 0 on `pds`.
Every number below is on the re-fit recipe, which is what October runs.

### 1.4 The production object (`basis250.py`)

| row | expected | measured | verdict |
|---|---|---|---|
| KA-P1 250's A/B/C production-shaped residual (300 rows, 233's re-extracted banks) vs 249's staged `cd4-off.npz` | bitwise | `mean_effect`, `gene_votes`, `target_condition_count`, target order all equal, max\|Δ\| 0.0 | **PASS** |
| KA-P2 `lib250.deflate(svd, 5)` on it vs 249's shipped `effects.cd4.deflation` | n_carried 291, energy removed 0.19866046506100588, top-5 0.1982359547102955, `U_5` sha `f06e0262…` | **identical, sha EQUAL** | **PASS** |

---

## 2. Phase 1 — the k ladder on the disjoint panels

Base = the panel's October stack: CD4 sized by 247's carrying-row rule at `g_cd4` 0.3443655278649641
(`def_pipeline.DEF_ENV`), H1 union at 0.28 carrying-row, operator = 233's / 238's per-panel fit; on
`abc` the production all-rows rule at f = 0.34 (the `BASE_carrying` row on `abc` reads −0.026 % /
−0.00001 `pds`: 247's numeric identity holds on the folds). Deflation is 246/249's construction
(per-gene centring over the carried rows, top-k left-singular projector, re-zeroed outside the
support), **resolved** (fraction re-solved on the deflated object, which is what `apply_frozen_cd4_
residual` does). Three cancer folds, `kolf2` apart. Price = 246 §4's convention: `cd4f034`'s rate card
on the macro cosine at k_cos = 1 plus PLOLO `pds` at transfer 0.7× / 1.0× / 2.0×; gate +0.0019
(shape-only class, 249 F5: amplitude −0.41 %).

Columns: rel aggregate cosine %, rel macro cosine %, Δ`pds`, Δq_w (`a15_sgn_w`), folds with `pds` > 0,
energy removed, solved CD4 fraction, price at 0.7× / 1.0× / 2.0×.

### 2.1 `abc` (today's panel, production sizing) — the control ladder

| k | agg % | macro % | Δ`pds` | Δq_w | +pds | E rm | frac | 0.7× | 1.0× | 2.0× |
|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|
| 2 | −0.023 | −0.582 | +0.00791 | +0.00107 | 3/3 | 0.118 | 0.387 | +0.0016 | +0.0025 | +0.0055 |
| 3 | −0.179 | −0.586 | +0.00936 | +0.00143 | 3/3 | 0.151 | 0.391 | +0.0020 | +0.0031 | +0.0066 |
| **5** | **−0.474** | **−0.486** | **+0.01241** | +0.00152 | 3/3 | 0.201 | 0.396 | **+0.0029** | **+0.0043** | **+0.0089** |
| 7 | −1.110 | −1.004 | +0.01126 | +0.00149 | 3/3 | 0.238 | 0.401 | +0.0022 | +0.0035 | +0.0077 |
| 10 | −1.594 | −1.439 | +0.01287 | +0.00067 | 3/3 | 0.283 | 0.414 | +0.0023 | +0.0038 | +0.0086 |
| 15 | −3.804 | −2.576 | +0.01382 | +0.00030 | 3/3 | 0.340 | 0.423 | +0.0018 | +0.0033 | +0.0084 |
| 5, fixed sizing | −1.254 | −0.964 | +0.01237 | +0.00146 | 3/3 | 0.201 | 0.374 | +0.0025 | +0.0039 | +0.0085 |
| 5, carrying @ `g_cd4` | −0.614 | −0.571 | +0.01236 | +0.00151 | 3/3 | 0.201 | 0.392 | +0.0028 | +0.0042 | +0.0088 |

Per fold at k = 5: k562 −1.69 % / +0.00905, hct116 −2.03 % / +0.01062, hek293t +2.29 % / +0.01756.
k-max Δ`pds` at k = 15; **k = 5 delivers 0.898 of it** — at the edge of the 10 % plateau (`[10, 15]`
by the strict rule); price-optimal k = 5; leave-one-fold-out k (chosen on two folds by price, read on
the third) +0.0035 (k = 5 / 5 / 10).

### 2.2 `real` — 199's realistic draw (CD4 300/300, K562 300, Orion 300; 34 essential targets)

| k | agg % | macro % | Δ`pds` | Δq_w | +pds | E rm | frac | 0.7× | 1.0× | 2.0× |
|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|
| 2 | −1.383 | −0.666 | +0.01495 | −0.00017 | 3/3 | 0.129 | 0.391 | +0.0034 | +0.0051 | +0.0107 |
| 3 | −1.629 | −0.775 | +0.01438 | +0.00054 | 3/3 | 0.159 | 0.394 | +0.0032 | +0.0048 | +0.0102 |
| **5** | **−2.106** | **−1.042** | **+0.01484** | +0.00029 | 3/3 | 0.203 | 0.398 | **+0.0031** | **+0.0048** | **+0.0103** |
| 7 | −2.457 | −1.429 | +0.01536 | +0.00045 | 3/3 | 0.238 | 0.401 | +0.0030 | +0.0047 | +0.0104 |
| 10 | −2.790 | −1.129 | +0.01758 | −0.00013 | 3/3 | 0.281 | 0.408 | +0.0038 | +0.0057 | +0.0123 |
| 15 | −3.569 | −1.513 | +0.01872 | −0.00026 | 3/3 | 0.336 | 0.414 | +0.0038 | +0.0059 | +0.0129 |
| 5, fixed | −2.525 | −1.340 | +0.01518 | +0.00012 | 3/3 | 0.203 | 0.370 | +0.0030 | +0.0047 | +0.0104 |
| 5, all-rows rule | −2.133 | −1.060 | +0.01488 | +0.00023 | 3/3 | 0.203 | 0.396 | +0.0031 | +0.0048 | +0.0103 |

Per fold at k = 5: **k562 −5.42 % / +0.00706**, hct116 −1.00 % / +0.02055, hek293t +0.11 % /
+0.01692 (the k562 fold alone prices **−0.0011** at 1.0×). k-max at 15; **k = 5 delivers 0.793** —
outside the plateau (`[10, 15]`, knee 10); price-optimal k = 15 (+0.0059; k = 10 +0.0057, k = 5
+0.0048 — a 0.001 spread, inside PLOLO's own `pds` σ in overall units); LOFO +0.0053.

### 2.3 `pess` — 199's pessimistic draw (CD4 174/300, K562 150, Orion 291/293; 30 essential targets)

| k | agg % | macro % | Δ`pds` | Δq_w | +pds | E rm | frac | 0.7× | 1.0× | 2.0× |
|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|
| 2 | −4.538 | −2.262 | +0.00476 | −0.00047 | 2/3 | 0.131 | 0.435 | −0.0004 | +0.0001 | +0.0019 |
| 3 | −5.222 | −2.278 | +0.00474 | −0.00009 | 2/3 | 0.161 | 0.441 | −0.0004 | +0.0001 | +0.0019 |
| **5** | **−6.602** | **−2.609** | **+0.00450** | −0.00074 | 3/3 | 0.212 | 0.446 | **−0.0007** | **−0.0002** | **+0.0015** |
| 7 | −7.442 | −3.184 | +0.00394 | −0.00021 | 2/3 | 0.253 | 0.450 | −0.0013 | −0.0008 | +0.0006 |
| 10 | −7.497 | −3.108 | +0.00473 | −0.00020 | 2/3 | 0.305 | 0.456 | −0.0010 | −0.0005 | +0.0013 |
| 15 | −7.892 | −3.130 | +0.00742 | −0.00022 | 2/3 | 0.373 | 0.473 | −0.0003 | +0.0005 | +0.0033 |
| 5, fixed | −6.574 | −2.523 | +0.00496 | −0.00060 | 3/3 | 0.212 | 0.415 | −0.0005 | +0.0000 | +0.0019 |
| 5, all-rows rule | −6.769 | −3.053 | +0.00335 | −0.00092 | 2/3 | 0.212 | 0.508 | −0.0013 | −0.0010 | +0.0003 |

Per fold at k = 5: k562 −7.33 % / +0.00725, hct116 −8.49 % / +0.00579, hek293t −3.98 % / +0.00047.
k-max at 15; k = 5 delivers 0.607; **no k prices above +0.0005 at 1.0×**; LOFO −0.0004. (The
all-rows-rule reference row `BASE233` sits at +0.60 % / −0.00166 `pds` relative to the carrying base:
247's over-dose, ×1.196, costs `pds` here.)

### 2.4 `efp` — 238's essential-free pessimistic draw (CD4 186/300, K562 146, Orion 293/296)

| k | agg % | macro % | Δ`pds` | Δq_w | +pds | E rm | frac | 0.7× | 1.0× | 2.0× |
|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|
| 2 | −2.709 | −3.738 | +0.00461 | +0.00014 | 2/3 | 0.145 | 0.432 | −0.0015 | −0.0010 | +0.0007 |
| 3 | −3.396 | −4.323 | +0.00506 | +0.00012 | 2/3 | 0.185 | 0.435 | −0.0018 | −0.0012 | +0.0006 |
| **5** | **−5.228** | **−5.475** | **+0.00354** | −0.00012 | 2/3 | 0.246 | 0.443 | **−0.0030** | **−0.0026** | **−0.0013** |
| 7 | −5.030 | −5.273 | +0.00383 | −0.00005 | 2/3 | 0.294 | 0.449 | −0.0028 | −0.0024 | −0.0010 |
| 10 | −5.686 | −5.264 | +0.00567 | +0.00001 | 2/3 | 0.349 | 0.460 | −0.0023 | −0.0017 | +0.0004 |
| 15 | −8.779 | −6.831 | +0.00750 | −0.00003 | 2/3 | 0.419 | 0.473 | −0.0030 | −0.0021 | +0.0007 |
| 5, fixed | −5.591 | −5.507 | +0.00406 | −0.00010 | 2/3 | 0.246 | 0.414 | −0.0029 | −0.0025 | −0.0010 |
| 5, all-rows rule | −4.780 | −5.566 | +0.00262 | −0.00009 | 2/3 | 0.246 | 0.486 | −0.0033 | −0.0031 | −0.0021 |

Per fold at k = 5: k562 −5.68 % / +0.00477, hct116 −2.56 % / +0.01325, **hek293t −7.44 % / −0.00741**.
k-max at 15; k = 5 delivers 0.471; **every k prices negative at 0.7× and 1.0×**; LOFO −0.0030.

### 2.5 `efr` — 238's essential-free realistic draw (CD4 300/300, K562 300, Orion 300)

| k | agg % | macro % | Δ`pds` | Δq_w | +pds | E rm | frac | 0.7× | 1.0× | 2.0× |
|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|
| 2 | −3.770 | −3.327 | +0.00550 | −0.00075 | 3/3 | 0.129 | 0.369 | −0.0010 | −0.0004 | +0.0017 |
| 3 | −3.376 | −2.435 | +0.00679 | −0.00074 | 3/3 | 0.155 | 0.372 | +0.0000 | +0.0008 | +0.0033 |
| **5** | **−3.523** | **−2.519** | **+0.00991** | −0.00052 | 3/3 | 0.199 | 0.379 | **+0.0008** | **+0.0019** | **+0.0056** |
| 7 | −3.177 | −2.245 | +0.01120 | −0.00077 | 3/3 | 0.234 | 0.386 | +0.0013 | +0.0026 | +0.0067 |
| **10** | −2.839 | −1.183 | **+0.01864** | +0.00023 | 3/3 | 0.270 | 0.391 | +0.0040 | **+0.0061** | +0.0130 |
| 15 | −3.579 | −1.649 | +0.01906 | +0.00051 | 3/3 | 0.315 | 0.400 | +0.0038 | +0.0059 | +0.0130 |
| 5, fixed | −3.699 | −2.593 | +0.01078 | −0.00038 | 3/3 | 0.199 | 0.352 | +0.0009 | +0.0022 | +0.0062 |
| 5, all-rows rule | −3.521 | −2.518 | +0.00990 | −0.00052 | 3/3 | 0.199 | 0.380 | +0.0008 | +0.0019 | +0.0056 |

Per fold at k = 5: k562 +0.42 % / +0.01677, hct116 −6.25 % / +0.00712, hek293t −4.74 % / +0.00585.
k-max at 15; k = 5 delivers 0.520 (plateau `[10, 15]`, knee 10); price-optimal k = 10 (+0.0061, 3.2×
the gate); at k = 5 the price is +0.0019, exactly the gate; LOFO +0.0056.

### 2.6 The knee, summarised — and `kolf2`

| panel | CD4 carried / 300 | k-max Δ`pds` (k) | Δ`pds`(5) / max | plateau (≥ 90 %) | price-best k (1.0×) | LOFO price |
|---|---:|---:|---:|---|---|---:|
| `abc` | 291 (fold 264–286) | +0.01382 (15) | **0.898** | 10, 15 | **5** (+0.0043) | +0.0035 |
| `real` | 300 | +0.01872 (15) | **0.793** | 10, 15 | 15 (+0.0059; 10 +0.0057; 5 +0.0048) | +0.0053 |
| `pess` | 174 | +0.00742 (15) | 0.607 | 15 | none pays | −0.0004 |
| `efp` | 186 | +0.00750 (15) | 0.471 | 15 | none pays | −0.0030 |
| `efr` | 300 | +0.01906 (15) | 0.520 | 10, 15 | **10** (+0.0061; 5 +0.0019) | +0.0056 |

**246's "k = 3–5 knee, saturating by 10" is an A/B/C property.** On every disjoint panel the `pds`
keeps rising to k = 10–15 and k = 5 is below the plateau; but the cosine bill also keeps rising, so in
price terms k = 5 … 10 are within ~0.001 of each other on `real` and k = 10 is clearly better on `efr`.
**`kolf2` apart:** 0 / 273 of its targets are on `real`, `pess`, `efp` or `efr` (asserted at runtime,
273 / 300 on `abc`), so it is evaluable on `abc` only: svd3 / 5 / 10 read `pds` +0.00206 / +0.00414 /
+0.00489 at −2.79 / −3.47 / −4.91 % aggregate (k = 5 equals 246 / 249 to 0.0e+00), the energy-matched
random −0.00246, split-half +0.00292, the pinned A/B/C basis +0.00370.

---

## 3. Phase 2 — the controls that made 246's claim, at k = 5 on each disjoint panel

Resolved sizing (the production construction) unless marked fixed; energy-matched random = 246's
bisection on the random-subspace dimension until it removes ≥ the energy svd-5 removes (two seeds);
shuffled projector = the top-5 fitted on a row-deranged copy and applied to the real rows; split-half =
basis fitted on a random half of the carried rows, applied to the other half, both ways (two seeds).

| panel | svd5 Δ`pds` (agg %) | randE5 / b (dims) | randE5 / b fixed | SHUFDEF | split-half k3 / **k5** / k10 (seed a, b) | LT3O `pds` retained (LT1O) | LT3O cosine | non-carrying rows: n, max\|Δ\| vs BASE × preserve ratio |
|---|---|---|---|---|---|---|---|---|
| `abc` | **+0.01241** (−0.47) | −0.00161 / −0.00661 (52 / 51) | −0.00040 / −0.00509 | +0.00009 | +0.0075, +0.0073 / **+0.0112, +0.0108** / +0.0127, +0.0121 | **56.6 %** (82.3 %) | loss broad (LT3O −0.00054 vs full −0.00023) | 5 / 8 / 8; ≤ 2.2e−16 |
| `real` | **+0.01484** (−2.11) | **−0.00573 / −0.00473** (59 / 58) | −0.00428 / −0.00289 | +0.00051 | +0.0134, +0.0130 / **+0.0148, +0.0148** / +0.0167, +0.0164 | **62.5 %** (85.5 %) | broad (−0.00071 vs −0.00040) | 0 (300/300 carried) |
| `pess` | +0.00450 (−6.60) | **+0.00581** / −0.00655 (29 / 29) | +0.00672 / −0.00525 | −0.00009 | +0.0027, +0.0030 / +0.0031, +0.0045 / +0.0054, +0.0042 | **−58.8 %** (26.4 %) | broad (−0.00094 vs −0.00044) | 27 / 122 / 123; ≤ 8.9e−16 |
| `efp` | +0.00354 (−5.23), 2/3 | +0.00154 / −0.00229 (35 / 36) | +0.00245 / −0.00118 | −0.00141 | +0.0033, +0.0050 / +0.0044, +0.0046 / +0.0054, +0.0062 | **−70.4 %** (29.9 %) | broad (−0.00131 vs −0.00083) | 11 / 113 / 113; ≤ 4.4e−16 |
| `efr` | **+0.00991** (−3.52) | **−0.00488 / −0.00203** (59 / 58) | −0.00206 / +0.00036 | −0.00104 | +0.0064, +0.0067 / **+0.0085, +0.0089** / +0.0125, +0.0125 | **55.0 %** (83.0 %) | broad (−0.00060 vs −0.00029) | 0 |

Readings.
- **`real` and `efr` reproduce 246's whole falsifier set**: random energy removal buys nothing
  (−0.002 … −0.006 at 3× the cosine cost), the shuffled projector is inert, the out-of-row basis keeps
  100 % / 85–90 % of the gain (as on A/B/C: 90 %), the gain is ~40 % top-3-concentrated (A/B/C 43 %),
  and rows CD4 does not carry move by exactly the preserve scalar (1e−16).
- **`pess`**: a random 29-dimensional subspace (seed a) removing the same 22 % of the energy buys
  **+0.0058**, more than the top-5 (+0.0045); seed b −0.0066. With 117–163 carrying rows the
  energy-matched control is ±0.006 and the arm sits inside it; the gain is carried by ≤ 3 targets (LT3O
  −59 %). **Not resolvable as a direction-selection effect** — a well-supported negative.
- **`efp`**: 2/3 folds, the hek293t fold −0.0074; the random control straddles the arm (+0.0015 /
  −0.0023 vs +0.0035); LT3O −70 %. Negative.
- The split-half retention stays high everywhere (even on `pess`/`efp`, 100 %+ of a small number),
  which is the same signature 246 read as "the subspace is a property of the source": it is — the
  problem on the low-coverage panels is not the basis but that the arm has nothing to buy there.

---

## 4. Phase 3 — the projector's cross-panel identity

### 4.1 Principal angles between production-shaped `V_5` bases (`basis250.py`)

`V_5` = gene-space top-5 right-singular directions of the centred carried rows of each panel's
production-shaped CD4 residual (all 300 panel rows from the panel's own re-extracted banks; 291 / 300 /
174 / 186 / 300 carried on `abc` / `real` / `pess` / `efp` / `efr`). The `abc` object is 249's shipped
one (§1.4).

| pair | θ₁ | θ₂ | θ₃ | θ₄ | θ₅ |
|---|---:|---:|---:|---:|---:|
| abc ↔ real | 26.9 | 39.1 | 46.3 | 56.6 | 88.6 |
| abc ↔ pess | 27.3 | 38.5 | 60.6 | 78.0 | 87.6 |
| abc ↔ efp | 24.2 | 39.6 | 50.7 | 77.7 | 89.8 |
| abc ↔ efr | 24.2 | 33.8 | 42.3 | 65.2 | 88.8 |
| real ↔ pess / efp / efr | 27.7 / 29.8 / 27.3 | 36.3 / 48.0 / 31.5 | 47.8 / 55.3 / 40.8 | 56.9 / 70.4 / 54.4 | 83.0 / 80.6 / 76.9 |
| pess ↔ efp / efr; efp ↔ efr | 28.8 / 25.6; 25.5 | 54.0 / 34.6; 42.3 | 60.3 / 55.6; 60.2 | 74.2 / 76.8; 73.9 | 85.7 / 89.3; 82.0 |

Fold-level bases vs the same panel's production basis: 0.0° on `real`/`efr` (full coverage: the fold's
rows are all 300), ≤ 9° on `abc` (25° for θ₅ on k562), up to 31° on `pess` k562 (117 rows) and 76° for
θ₅ on `efp` hct116 — so the fold-level re-fits stand for a day-of re-fit, and the top-5 is less stable
below ~180 rows.

### 4.2 Energy captured by a 5-dim basis (fraction of the panel's centred CD4 energy)

| panel | own top-5 | by A/B/C's `V_5` | ratio | by `real` / `pess` / `efp` / `efr` | random-5 floor | spectrum λ₁…λ₆ / Σ |
|---|---:|---:|---:|---|---:|---|
| `abc` | 0.1982 | 0.1982 | 1.000 | 0.1108 / 0.0956 / 0.1101 / 0.1163 | 0.0005 | 0.080 0.038 0.032 0.029 0.020 0.019 |
| `real` | 0.2025 | **0.1171** | 0.578 | — / 0.1153 / 0.1176 / 0.1391 | 0.0005 | 0.091 0.037 0.030 0.024 0.020 0.019 |
| `pess` | 0.1993 | **0.1028** | 0.516 | 0.1185 / — / 0.0979 / 0.1096 | 0.0005 | 0.069 0.055 0.028 0.024 0.022 0.020 |
| `efp` | 0.2316 | **0.1334** | 0.576 | 0.1416 / 0.1179 / — / 0.1258 | 0.0005 | 0.089 0.048 0.037 0.030 0.028 0.026 |
| `efr` | 0.1984 | **0.1249** | 0.630 | 0.1445 / 0.1088 / 0.1108 / — | 0.0005 | 0.089 0.040 0.026 0.022 0.022 0.020 |

(Cross-capture cells off the A/B/C column are from `basis250.json`'s `cross` block.) **They are not
nearly the same subspace.** One to two directions are shared at 24–39° and the fifth is orthogonal;
A/B/C's basis captures 52–63 % of what a panel's own top-5 captures — 200× the random floor, half the
in-sample value. The top-5 energy fraction itself is 0.198–0.232 on every panel, i.e. the "0.19–0.21"
census reading does not move with the panel.

### 4.3 The pinned A/B/C production basis on each panel (resolved; `pinABC{k}` = `R_c(I − V_kV_kᵀ)` after own-row centring)

| panel | k | Δ`pds` (folds, σ) | agg % | macro % | Δq_w | E rm | price 0.7× / 1.0× / 2.0× | re-fit at the same k: Δ`pds`, agg %, price 1.0× |
|---|---|---|---:|---:|---:|---:|---|---|
| `abc` | 5 | +0.01222 (3/3, 8.2) | −0.44 | −0.51 | +0.00151 | 0.200 | +0.0028 / **+0.0042** / +0.0087 | +0.01241, −0.47, +0.0043 |
| `real` | 3 | +0.01301 (3/3, 8.7) | −0.11 | +0.12 | −0.00013 | 0.107 | +0.0035 / +0.0049 / +0.0098 | +0.01438, −1.63, +0.0048 |
| `real` | **5** | **+0.01532** (3/3, 10.3) | **+0.06** | **+0.48** | +0.00000 | 0.120 | +0.0043 / **+0.0061** / +0.0118 | +0.01484, −2.11, +0.0048 |
| `real` | 10 | +0.01693 (3/3, 11.3) | −0.24 | +0.35 | −0.00022 | 0.154 | +0.0047 / +0.0066 / +0.0129 | +0.01758, −2.79, +0.0057 |
| `pess` | 3 | +0.00356 (2/3, 2.4) | −1.27 | −0.33 | −0.00025 | 0.094 | +0.0007 / +0.0011 / +0.0024 | +0.00474, −5.22, +0.0001 |
| `pess` | **5** | +0.00369 (3/3, 2.5) | −1.51 | −0.66 | −0.00002 | 0.106 | +0.0005 / **+0.0009** / +0.0023 | +0.00450, −6.60, −0.0002 |
| `pess` | 10 | +0.00472 (3/3, 3.2) | −1.93 | −1.09 | +0.00002 | 0.136 | +0.0004 / +0.0010 / +0.0027 | +0.00473, −7.50, −0.0005 |
| `efp` | 3 | +0.00522 (2/3, 3.5) | −0.92 | −2.52 | −0.00016 | 0.124 | −0.0005 / +0.0001 / +0.0021 | +0.00506, −3.40, −0.0012 |
| `efp` | **5** | **+0.00729** (3/3, 4.9) | −0.53 | −1.78 | −0.00008 | 0.140 | +0.0006 / **+0.0014** / +0.0041 | +0.00354, −5.23, −0.0026 |
| `efp` | 10 | +0.00962 (3/3, 6.4) | −1.08 | −1.94 | +0.00005 | 0.178 | +0.0011 / +0.0022 / +0.0058 | +0.00567, −5.69, −0.0017 |
| `efr` | 3 | +0.00842 (3/3, 5.6) | −2.22 | −1.88 | −0.00017 | 0.114 | +0.0008 / +0.0018 / +0.0049 | +0.00679, −3.38, +0.0008 |
| `efr` | **5** | **+0.01065** (3/3, 7.1) | −1.99 | −1.36 | −0.00033 | 0.128 | +0.0018 / **+0.0030** / +0.0070 | +0.00991, −3.52, +0.0019 |
| `efr` | 10 | +0.01227 (3/3, 8.2) | −2.45 | −1.83 | +0.00020 | 0.167 | +0.0019 / +0.0032 / +0.0078 | +0.01864, −2.84, +0.0061 |

Leave-top-3-out retention of the pinned k = 5 gain: `abc` 55.6 %, `real` **67.7 %**, `pess` −58 %,
`efp` **+32.9 %** (the re-fit: −70 %), `efr` 61.6 %. Non-carrying rows exact (1e−16) for the pinned
arm too. **On 4/4 disjoint panels the pinned basis prices above the day-of re-fit at k = 5**
(+0.0061 vs +0.0048; +0.0009 vs −0.0002; +0.0014 vs −0.0026; +0.0030 vs +0.0019) and is never
negative; its cosine cost is 2–5 pp smaller because it removes only 10–14 % of the energy — the part
that is shared with A/B/C — instead of 20–25 %.

### 4.4 Is A/B/C's basis special? Cross-panel pinning (`cross250.py`, 15 runs)

Each panel deflated with every other panel's production `V_k`, resolved, 3 cancer folds
(Δ`pds` / rel agg %; folds with `pds` > 0 in brackets):

| panel deflated ↓ / basis → | pin `abc` | pin `real` | pin `pess` | pin `efp` | pin `efr` | own re-fit |
|---|---|---|---|---|---|---|
| `abc`, k = 5 | +0.0122 / −0.4 (3) | +0.0112 / +1.0 (3) | +0.0102 / +0.6 (3) | +0.0098 / +0.0 (3) | +0.0119 / +1.0 (3) | +0.0124 / −0.5 (3) |
| `real`, k = 5 | +0.0153 / +0.1 (3) | *+0.0148 / −2.1 (3)* | +0.0161 / −0.8 (3) | +0.0144 / −0.6 (3) | +0.0163 / +0.0 (3) | +0.0148 / −2.1 (3) |
| `pess`, k = 5 | +0.0037 / −1.5 (3) | +0.0038 / −2.5 (2) | *+0.0038 / −6.6 (3)* | +0.0040 / −1.7 (2) | +0.0050 / −1.1 (2) | +0.0045 / −6.6 (3) |
| `efp`, k = 5 | +0.0073 / −0.5 (3) | +0.0083 / −0.3 (3) | +0.0073 / −0.5 (3) | *+0.0033 / −5.4 (2)* | +0.0075 / −0.2 (2) | +0.0035 / −5.2 (2) |
| `efr`, k = 5 | +0.0106 / −2.0 (3) | +0.0125 / −2.3 (3) | +0.0083 / −3.3 (3) | +0.0100 / −2.6 (3) | *+0.0099 / −3.5 (3)* | +0.0099 / −3.5 (3) |
| `abc`, k = 10 | +0.0129 / −1.6 | +0.0112 / +0.7 | +0.0115 / +0.6 | +0.0118 / +0.4 | +0.0140 / +1.4 | +0.0129 / −1.6 |
| `real`, k = 10 | +0.0169 / −0.2 | *+0.0176 / −2.8* | +0.0171 / −0.9 | +0.0153 / −0.5 | +0.0174 / −0.6 | +0.0176 / −2.8 |
| `pess`, k = 10 | +0.0047 / −1.9 | +0.0049 / −3.2 | *+0.0043 / −7.5* | +0.0039 / −1.9 | +0.0057 / −2.2 | +0.0047 / −7.5 |
| `efp`, k = 10 | +0.0096 / −1.1 | +0.0094 / −0.4 | +0.0093 / −0.7 | *+0.0061 / −5.3* | +0.0105 / −0.1 | +0.0057 / −5.7 |
| `efr`, k = 10 | +0.0123 / −2.5 | +0.0130 / −2.2 | +0.0127 / −2.7 | +0.0112 / −2.4 | *+0.0186 / −2.8* | +0.0186 / −2.8 |

(Italic = the panel's own production basis, which on `real`/`efr` is the fold re-fit itself.) Mean over
the four **foreign** bases vs the own re-fit, k = 5 — Δ`pds`, agg %, price 1.0× (worst foreign price):

| panel | own re-fit | foreign mean | foreign range of Δ`pds` |
|---|---|---|---|
| `abc` | +0.0124, −0.47 %, +0.0043 | +0.0108, **+0.68 %**, +0.0039 (+0.0030) | +0.0098 … +0.0119 |
| `real` | +0.0148, −2.11 %, +0.0048 | +0.0155, **−0.32 %**, +0.0059 (+0.0052) | +0.0144 … +0.0163 |
| `pess` | +0.0045, −6.60 %, −0.0002 | +0.0041, **−1.72 %**, +0.0009 (+0.0002) | +0.0037 … +0.0050 |
| `efp` | +0.0035, −5.23 %, −0.0026 | **+0.0076**, **−0.39 %**, +0.0017 (+0.0014) | +0.0073 … +0.0083 |
| `efr` | +0.0099, −3.52 %, +0.0019 | +0.0104, −2.53 %, +0.0024 (+0.0008) | +0.0083 … +0.0125 |

**A basis fitted on any other panel's rows is at least as good as the own-row fit on 5/5 panels and
costs 1–5 pp less cosine.** So the A/B/C basis is not special — the source is: the directions the
panels share (θ₁–θ₂ at 24–39°) are the T-cell programme hub 246 identified, and they carry the `pds`
gain; the directions a same-row fit adds on top (θ₃–θ₅, 42–90°) are panel-specific structure whose
removal costs cosine and buys nothing. This is 246's split-half result read across panels: an
out-of-row basis is the cheaper projector.

### 4.5 Which is safer for October, and why

**Pin the A/B/C production basis (k = 5, `U_5` sha `f06e0262…` → gene-space `V_5`) rather than
re-fitting on the D/E/F rows.** Four numbers: (i) price on the four disjoint panels +0.0061 / +0.0009 /
+0.0014 / +0.0030 vs the re-fit's +0.0048 / −0.0002 / −0.0026 / +0.0019 — the re-fit is negative on two
draws, the pinned basis on none; (ii) cosine cost −0.1 … −2.0 % vs −2.1 … −6.6 %; (iii) it is the object
whose board point is being bought on 2026-09-05 (`cd4svd5` is exactly this projector on A/B/C), so its
`pds` transfer factor will be a measured number before Oct 22 while a re-fit's would not; (iv) it removes
a day-of fit from the fail-closed path (no dependence on the census, on `n_carried`, or on the top-5
stability, which degrades below ~180 rows: θ₅ up to 76° between a fold re-fit and the production one on
`efp`). **What pinning needs and this experiment could not do (read-only mandate):** a builder code path
that applies a gene-space projector `R_c(I − V₅V₅ᵀ)` after centring on the D/E/F rows (the present
`deflate_cd4_residual` fits a row-space `U_k` on the rows it is applied to; its `basis` kwarg takes a
row-space matrix, not a gene-space one), a flag distinct from `VCC_CD4_DEFLATE_K` (rule 10: a pinned
basis and a re-fit are two variables), 249-style known-answer rows (on A/B/C the pinned projector must
reproduce the `cd4svd5` effect matrix `a28d5be2…` bit-for-bit, since on A/B/C's own rows the two
constructions coincide), and a rehearsal on 199's bundle. The support caveat: `V_5` has coordinates on
A/B/C's panel genes, which are zero on a D/E/F residual; the projection is re-zeroed outside the day's
support exactly as 246/249 do for partial rows.

---

## 5. Phase 4 — the `+partial` composition on the disjoint panels

CD4's one-/two-condition rows (production zeroes them), carrying-row sized at 0.34, composed with the
k = 5 re-fit deflation (resolved). Rows per panel: 13–15 (`abc`), 14 (`real`), 6–7 (`pess`), 3–4
(`efp`), 10 (`efr`).

| panel | `partial@0.34` alone: agg / macro / Δ`pds` (folds) | shuffled partial | `svd5+partial@0.34`: agg / macro / Δ`pds` (folds) | increment over `svd5`: Δ`pds` (folds), agg | LT3O of the composed | price 0.7× / 1.0× / 2.0× |
|---|---|---|---|---|---|---|
| `abc` | +1.20 % / +0.60 % / +0.00197 (3/3) | −0.23 % / −0.00009 | +0.73 % / +0.11 % / **+0.01438** (3/3) | +0.00197 (3/3), +1.21 % | 62.6 % | +0.0038 / +0.0054 / +0.0108 |
| `real` | +0.52 % / +0.52 % / +0.00097 (1/3) | −0.20 % / +0.00033 | −1.57 % / −0.51 % / +0.01581 (3/3) | +0.00097 (1/3), +0.54 % | 64.9 % | +0.0038 / +0.0055 / +0.0114 |
| `pess` | +1.06 % / +1.04 % / +0.00085 (2/3) | +0.18 % / +0.00004 | −5.53 % / −1.57 % / +0.00535 (3/3) | +0.00085 (2/3), +1.14 % | −33 % | +0.0003 / +0.0009 / +0.0029 |
| `efp` | +0.29 % / +0.27 % / +0.00023 (3/3) | −0.10 % / −0.00097 | −4.93 % / −5.20 % / +0.00374 (2/3) | +0.00020 (3/3), +0.32 % | −61 % | −0.0028 / −0.0024 / −0.0010 |
| `efr` | +0.21 % / −0.40 % / **−0.00199** (0/3) | +0.87 % / +0.00015 | −3.30 % / −2.91 % / +0.00792 (3/3) | −0.00200 (0/3), +0.17 % | 43 % | −0.0000 / +0.0008 / +0.0038 |

246's A/B/C reading (+0.0018 `pds`, +0.6 % macro, additive) reproduces on `abc` (+0.00197, additive to
0.0000). Off A/B/C the rider is +0.0010 / +0.0009 / +0.0002 / **−0.0020** with a small positive cosine —
the shuffled control sits at 0 — i.e. 6–14 rows of un-whitened CD4 on empty rows buy cosine, not
retrieval, and the sign of the retrieval flips on `efr`. **Not build-worthy as an October rider.**

---

## 6. Verdicts — one line per panel

Mechanism = `pds` > 0 on 3/3 folds, ≥ 3σ (162's LOO σ 0.001494 raw), both energy-matched random
seeds below half the arm, shuffled projector < +0.001, split-half retains ≥ 50 %, LT3O ≥ 25 %.
Pays = six-member price at 1.0× `pds` transfer above the +0.0019 shape-only gate.

| panel | verdict | k = 5 re-fit: Δ`pds` (folds, σ), agg %, price 0.7 / 1.0 / 2.0× | k to ship on Oct 22 | pinned A/B/C k = 5 price |
|---|---|---|---|---|
| `abc` (control) | TRANSFERS (246's own reading, re-fit recipe) | +0.01241 (3/3, 8.3σ), −0.47 %, +0.0029 / +0.0043 / +0.0089 | 5 (price-optimal; Δ`pds`(5)/max 0.90) | +0.0042 |
| **`real`** | **TRANSFERS** — every control passes | +0.01484 (3/3, 9.9σ), −2.11 %, +0.0031 / +0.0048 / +0.0103 | re-fit: any of 5 … 15 (price +0.0048 … +0.0059; `pds` plateau k ≥ 10, k = 5 at 0.79 of max); **pinned A/B/C k = 5: +0.0061** | **+0.0061** |
| **`efr`** | **TRANSFERS** — every control passes; pays at k ≥ 7 | +0.00991 (3/3, 6.6σ), −3.52 %, +0.0008 / +0.0019 / +0.0056 | re-fit: **10** (+0.0061; k = 5 at the gate, 0.52 of max); **pinned A/B/C k = 5: +0.0030** | +0.0030 |
| **`pess`** | **DOES NOT** — random control seed a +0.0058 > arm, LT3O −59 %, no k pays | +0.00450 (3/3, 3.0σ), −6.60 %, −0.0007 / −0.0002 / +0.0015 | unset (pinned A/B/C stays +0.0009, below the gate but ≥ 0) | +0.0009 |
| **`efp`** | **DOES NOT** — 2/3 folds, random control straddles the arm, LT3O −70 %, every k negative | +0.00354 (2/3, 2.4σ), −5.23 %, −0.0030 / −0.0026 / −0.0013 | unset (pinned A/B/C +0.0073 `pds` 3/3, price +0.0014, LT3O +33 %) | +0.0014 |

### 6.1 What separates the two outcomes

| | `real` | `efr` | `pess` | `efp` |
|---|---:|---:|---:|---:|
| CD4 carried rows / 300 (norm > 0) | 300 (286) | 300 (290) | **174 (167)** | **186 (182)** |
| K562 fold targets | 300 | 300 | **150** | **146** |
| essential targets on the panel | 34 | 0 | 30 | 0 |
| own top-5 energy fraction | 0.203 | 0.198 | 0.199 | 0.232 |
| energy captured by A/B/C's `V_5` | 0.117 | 0.125 | 0.103 | 0.133 |
| k = 5 re-fit price (1.0×) | +0.0048 | +0.0019 | −0.0002 | −0.0026 |

Coverage splits them; the essential-gene content does not (`efr` transfers with 0 essential targets,
`pess` fails with 30); the top-5 energy fraction and the cross-capture do not. This is the same
regime 233 §4.1 found for the CD4 residual itself (its `pds` contribution turns negative at 58 %
coverage): **where CD4 covers the panel, its hub is worth removing; where it covers half of it, the
residual has too little to give and its top directions are half noise** (fold-vs-production θ₅ 31–76°).
Two panels per outcome is what four draws allow; CD4 coverage and K562 coverage move together in
199's construction, so the two are not separated here.

### 6.2 The falsifier the D/E/F census would fire

**Not the top-5 energy fraction.** It reads 0.203 / 0.199 / 0.232 / 0.198 on `real` / `pess` / `efp` /
`efr` (day-of re-fit on the fold objects 0.199–0.265): `pess` sits inside the A/B/C range 0.19–0.21 and
fails, `efp` sits outside it and fails, `real`/`efr` sit inside and transfer. A census gate on it would
fire on the wrong panel. The readings that did separate the outcomes here, in the order the day can
read them:

1. **CD4 carried rows on the D/E/F panel** (`effects.cd4.fully_supported_targets` + partial rows):
   ≥ 290 / 300 on both transferring panels, ≤ 186 on both failing ones. DEF-RUNBOOK §5's pessimistic
   draw (CD4 union 174) is the failing regime. Propose: **re-fit deflation only if carried ≥ 250;
   below it, unset** (or pinned only, which never priced negative here).
2. **Energy captured by the pinned A/B/C `V_5`** (a 2 s computation on the day's residual): 0.103–0.133
   on the four disjoint panels, 0.198 in-sample, 0.0005 random. A reading < 0.05 means the hub is absent
   from the day's residual (nothing measured here covers that); > 0.19 means a panel that behaves like
   A/B/C's own rows.
3. `effects.cd4.deflation.energy_removed_frac` of a re-fit: 0.199–0.246 here (0.1987 shipped); a reading
   outside 0.18–0.27 is off the measured distribution.

---

## 7. Proposed `DEF-RUNBOOK` §3 decision-table row (text only — the runbook is not edited)

| flag / choice | SHIPS | validated by | if D/E/F look different |
|---|---|---|---|
| **`VCC_CD4_DEFLATE_K`** (249) | **Unset until the 2026-09-05 `cd4svd5` board point is read.** If it reads board raw `pds` up against the incumbent's 0.83572 (246 predicted +0.009 … +0.025 at 0.7 … 2.0×): ship the deflation on D/E/F **as the PINNED A/B/C basis, k = 5** (`U_5` sha `f06e0262…`, gene-space `V_5` from 249's staged residual, applied after centring on the D/E/F rows) — **never a day-of re-fit of the projector**. Pinning needs a builder path distinct from `VCC_CD4_DEFLATE_K` (rule 10), a 249-style byte-identity known answer on A/B/C (`a28d5be2…`) and a 199-bundle rehearsal, all before Oct 15; **until that exists the flag reads: `VCC_CD4_DEFLATE_K=5` (re-fit) only if the day's CD4 carried rows ≥ 290/300, otherwise unset** | 246 (A/B/C PLOLO `pds` +0.01226, 3/3, every falsifier); 249 (built, flag-off byte-inert); **250**: on four disjoint 300-target panels the day-of re-fit at k = 5 reads `pds` +0.0148 / +0.0045 / +0.0035 / +0.0099 (`real` / `pess` / `efp` / `efr`) at aggregate cosine −2.1 / −6.6 / −5.2 / −3.5 % (A/B/C −0.47 %), six-member price +0.0048 / −0.0002 / −0.0026 / +0.0019 — it clears +0.0019 only at full CD4 coverage; the pinned A/B/C basis reads `pds` +0.0153 / +0.0037 / +0.0073 / +0.0107 at −0.1 … −2.0 % cosine, price +0.0061 / +0.0009 / +0.0014 / +0.0030, never negative, better than the re-fit on 4/4; any foreign basis behaves the same (`efp` +0.0073 … +0.0083 vs own +0.0035); the `pds` knee off A/B/C is k = 10–15 but the price is flat past k = 5–10 | Census: **CD4 carried rows < 250/300 → the re-fit variant reads ≤ 0 (both failing draws: 174 / 186), ship unset**; the top-5 energy fraction (0.198–0.232 on all four, 0.198 shipped) does NOT discriminate and must not gate; record the fraction of the day's centred CD4 energy captured by the pinned `V_5` (0.10–0.14 measured; < 0.05 or > 0.19 is unmeasured territory) and `energy_removed_frac` (0.20–0.25 measured). Never change k on the day: no k pays on a failing panel and k = 5 … 10 are within 0.001 on a transferring one |

---

## 8. Honest limits

- **Random draws, not Arc-selected panels.** 199's / 238's panels are drawn from donor coverage;
  Arc's selection protocol produced a panel 1.3–2.0× more cross-line reproducible than a random draw
  (233 §6). October sits somewhere between `real`/`efr` (full coverage: transfers, pays) and
  `pess`/`efp` (half coverage: nothing to buy); the board's own CD4 coverage on Oct 22 decides which.
- **One instrument, one class of truth** — 162's PLOLO with donor-library pseudobulk truths on lines we
  already hold; the context half of the transition is unmeasured (233 §0). Board `pds` transfer for a
  gain has one calibration point (0.70×, `h1union`); prices are quoted at 0.7 / 1.0 / 2.0×. The
  2026-09-05 board point will fix this factor for the pinned object.
- **The cosine cost is fold-asymmetric everywhere** (`real` k562 −5.4 % while the Orion folds are
  ≈ 0; `efr` the reverse; `abc` hek293t +2.3 %). A reader weighting the worst fold prices the `real`
  re-fit at k = 5 at **−0.0011** and the pinned basis at +0.0014.
- **The pinned variant is not what 249 built.** 249's builder fits `U_k` on the rows it deflates; the
  pinned gene-space projector is a different code path that this read-only experiment did not write,
  test or rehearse. On A/B/C's own rows the two coincide (this is why the `cd4svd5` board point
  calibrates the pinned object), off A/B/C they do not (§4.3).
- **Two panels per outcome.** Coverage is the separating variable on these four, but CD4 coverage and
  K562 coverage co-vary in 199's draws; a full-coverage-CD4 / low-coverage-K562 panel was not available.
- **The `pess` random-control reading is two seeds** (+0.0058 / −0.0066): the honest statement is that
  the arm is inside the control's spread there, not that random removal "works".
- **`nmae` / `fid` / `reach` / `jac`** are priced only through 246 §4's rate card (cosine legs); 246's
  exact-`nmae` bound (|Δ| ≤ 0.0004) was not re-derived off A/B/C.
- The essential-free draws use 238's operators fitted through 235's changed fitter (238 §1) — both
  versions reproduced `e2f65ecf…`, and every operator loaded here passed the axis assertion.

---

## 9. Files

| file | what |
|---|---|
| `bind250.py` | the five-panel registry rebound onto `panels233` (233's / 238's artefacts read-only; every cache path under 250's `work/`); panel disjointness asserted at import |
| `gate250.py` | the memory gate (`--wait` polls; used before every heavy run) |
| `lib250.py` | staging through 233's verbatim modules, sizing (all-rows @ 0.34 vs carrying @ `g_cd4`), 246's `deflate` / energy-matched random / shuffled projector / split-half (byte copies), the gene-space basis (`gene_basis`, `deflate_gene_basis`, principal angles), the production A/B/C object (249's `cd4-off.npz`), `cd4_partial` per panel, 246's scorer (`resid218.score_effects`) |
| `sweep250.py` | one (panel, fold) per process: KA rows → `BASE` → k ladder → controls → pinned A/B/C → partial rows → non-carrying test; writes `work/sweep250-<panel>-<fold>.json`, `pertarget250-…npz`, `basis250-…npz` |
| `bridge250.py` | 246's instrument (062 operators) through 250's functions, 40 KA rows at 0.0e+00; then 233's re-fit A/B/C operator on the same objects → `work/bridge250-<fold>.json` |
| `basis250.py` | KA-P1/P2 (the shipped object and its deflation audit, bitwise), production-shaped `V_k` per panel, principal angles, cross-captured energy, spectra, fold-level angles → `work/basis250.json`, `basis250-production.npz` |
| `cross250.py` | every panel deflated with every other panel's production `V_5` / `V_10` → `work/cross250-<panel>-<fold>.json`, `work/cross250-summary.json` |
| `analyse250.py` | every table above → `work/analyse250.json`, `work/analyse250.log` |
| `run250.sh`, `runcross250.sh` | the serial drivers (gate + `/usr/bin/time -l`), logs under `work/logs/` |

```bash
cd experiments/250-deflation-october
./run250.sh          # 16 sweeps + 4 bridges + basis, ~11 min, peak RSS ≤ 5.11 GiB
./runcross250.sh     # 15 cross-pinned runs, ~6 min
uv run python analyse250.py | tee work/analyse250.log
```

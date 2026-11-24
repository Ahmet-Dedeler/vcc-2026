# 277 — Rendering the identifiable LOW-CONSERVATION rows differently: one-variable arms selected by 276's leave-the-held-line-out conservation map, priced on PLOLO and on the Flex-400 reference

## DOES NOT CLEAR — and the negative is clean.

**0 of 15 priced cells clears anything.** The headline arm the mandate names — **V1, the un-deflated CD4
residual on the bottom 30 % of rows by 276's `C_lolo`, deflated everywhere else** — reads Δ`pds`
**−0.00431 [−0.00917, +0.00055]** on the three cancer PLOLO folds, **1/3 folds positive**, `cos_agg`
**−1.10 %**, and prices **−0.00455 pessimistic / −0.00430 point** against a **+0.0019** gate: it is
2.4 gates on the *wrong side*. Its mirror (V2, deflate the class and un-deflate everything else) is worse
(−0.00825, 0/3); removing the residual on the class altogether (V3) is worse (−0.00707, 0/3); the boundary
arms (bottom 20 % / 40 %) are the same to within noise. **The one cell that is 3/3 positive is V4, 272's
eigenvalue floor applied only to the class** (+0.00063 [−0.00057, +0.00183], pessimistic **+0.00078** =
0.41× the gate) — and its own controls kill it: the identical floor built from a **shuffled** noise vector
reads **+0.00102** and from a **mirrored** one **+0.00078**, both ≥ the true one, so nothing in it is the
noise estimate.

**The row selector itself does not survive its proper null.** Against 24 size-matched random row sets per
fold running V1's identical construction, the conservation class sits at **z −0.88 / +0.83 / −0.98**
(8 % / 83 % / 21 % of draws at or below), **pooled z −0.53** — the map does not identify rows where
un-deflating the CD4 residual behaves differently from an arbitrary row set. The pipeline's own
donor–donor agreement reads pooled z −0.73 on the same test, i.e. no worse.

**The mechanism is measured, and it is the mandate's hypothesis run in reverse.** On the class's
*anti-aligned* rows (262's group e, base·truth ≤ 0) the un-deflated residual does help, **+0.03384** on 78
(target, fold) rows — but the same edit costs **−0.03858** on the class's 163 *marginal* rows (253 edited
rows in total across the three folds), and the marginal band is twice as large. That trade is not specific
to the class: random row sets show the same
signs (+0.056 / −0.016, +0.034 / −0.026, +0.020 / −0.021). Across the class the mean per-target cosine move
is **−0.00001**; 10 rows cross the 0.02 knee upward and 12 downward. **A different rendering of the CD4
residual does not move these rows across the knee — it reshuffles them across it at a small net loss**,
which is F4/262's finding ("the residual cannot repair a row whose whitened base points away from the held
line") holding for the *rendering* axis as well as the selector and dose axes.

**The competition-facing framing 276 §6.2 opened is therefore closed at this construction.** 276's
recommendation 1 (no dose from this map) already held; 277 adds: no *rendering* either, on the five shapes
tested (un-deflated / mirrored / removed / boundary-varied / operator-damped), on 3 PLOLO folds and the
Flex-400 reference, priced pessimistically.

**Date 2026-09-06, 05:40–06:20 UTC (08:40–09:20 Ankara).** Deterministic (seed 20260906). **Measurement
only: no build, no `vcc`, no `cell-eval2`, no submission, nothing under `submissions/` or `data/`, no
production file touched, no download.** `uv run python` throughout, one PLOLO fold per process, **peak RSS
4.90 GiB**, ~237 s of compute across 224 scored compositions.
`find data scripts src docs submissions experiments -newer work/DEFINITION-277.md` outside this directory
returns exactly one path, `data/leaderboard/2026-09-06T0614Z.json` — the board-polling cron's own snapshot,
not written by this experiment.

---

## 0. Known answers — 28 rows, all PASS; 25 of them at exactly 0.0e+00

| row | required | measured | verdict |
|---|---|---|---|
| **KA-1** V0 (`cd4svd5`) endpoints vs 257 `sweep257-<fold>.json` row `k5_s1` | < 1e-9 | **0.0e+00** on all 4 folds | PASS |
| **KA-1** V0 per-target `pds` / own cosine vs 262 `diag262-<fold>.npz` `cd4svd5|pds` / `|own` | 0.0 | **0.0e+00** on all 4 folds | PASS |
| **KA-1** `h1union` (the un-deflated stack) vs 257 `STACK` and 262 `h1union|pds` / `|own` | 0.0 | **0.0e+00** on all 4 folds | PASS |
| **KA-1b** 262's `cov|base_cos_truth` recomputed here | 0.0 | **0.0e+00** on all 4 folds | PASS |
| **KA-2** 273's `BASE` / `cd4h1:off` / `INCUMBENT` per-target arrays (a23, a15, effect level; 18 arrays) re-composed from 170's staged rows | 0.0 | **0.0e+00** | PASS |
| **KA-3** `mixed_whiten` with CLASS = ∅ == `plolo.Operator.apply` | 0.0 | **0.0e+00** on all 4 folds | PASS |
| **KA-4** BR249 retrodicted: `cd4svd5` vs `h1union` through the recipe at h1union's terms | +0.0043 ± 0.0008, and == 257's R1 to 1e-9 | `pds` **+0.00552**, `mse` **−0.00083**, DE-three −0.00035, `nmae` +0.00038 = **+0.00472**; \|Δ vs 257 R1\| **0.0e+00**; board realised **+0.0042634** | PASS |
| **KA-0 (floor277)** the spliced whitening path with no row edited == `render277`'s V0 | < 1e-12 | **0.0e+00** ×3 | PASS |
| **KA-1 (floor277)** the floor on ALL rows == 272's `floor:c1` per-fold Δ`pds` | < 1e-9 | **0.0e+00** ×3 | PASS |
| **KA (null277)** the seed-20260906 draw == `render277`'s `C_rand1` | < 1e-12 | 1.6e-16 / 7.3e-17 / 1.6e-16 (the two paths sum the same per-target array in a different order) | PASS |
| **F4 reproduction** on the Flex-400 reference: the CD4 residual's `cos_agg` gain over the base, un-deflated vs deflated | R273 §4: +30 % / +24 % | **+30.3 % / +24.1 %** | PASS |

`work/DEFINITION-277.md` (the class rule, the six constructions, the controls, the bands, the CLEARED
rule) was frozen at **2026-09-06 05:46 UTC**, sha256
`05a8939171874c24f22e5fd29c40ce1ad63c44544d951c0870feffc34052ddd1`, before `render277.py` composed
anything. Everything decided afterwards is labelled POST-HOC.

---

## 1. The class

**CLASS30(fold)** = the `floor(0.30 · n_eligible)` targets with the smallest `C_lolo_<fold>` in
`experiments/276-conservation-map/work/conservation276.csv` (sha `e5388e3f…`), ties by ascending symbol;
one target per cancer fold has no `C_lolo` and is never in a class. For the Flex-400 reference the score is
`h1_c_lolo` in `predict276.npz` (sha `a3a9f289…`) — **276 built it on h1-test's own 85 targets with the
identical construction and both H1 banks removed** (`predict276.py:263`, `lib276.HELD_LINE_SOURCES["H1"]`),
so it is already the H1-excluded score the mandate asks for; it is used as-is, not recomputed.

| fold | targets | with `C_lolo` | CLASS30 n | `C_lolo` cut | CD4-carried in class | CLASS20 / CLASS40 |
|---|---:|---:|---:|---:|---:|---:|
| `k562_gwps` | 269 | 268 | 80 | +0.00751 | 78 | 53 / 107 |
| `orion_hct116` | 290 | 289 | 86 | +0.00511 | 83 | 57 / 115 |
| `orion_hek293t` | 294 | 293 | 87 | +0.00635 | 84 | 58 / 117 |
| `kolf2` | 273 | 273 | 81 | +0.00798 | 79 | 54 / 109 |
| h1-test (Flex-400) | 85 | 85 | 25 | +0.01251 | 21 | — |

2–3 class rows per fold carry no complete CD4 row and are therefore un-editable by V1–V3; they are left in
the class (the mask is declared on the map, not on the residual) and contribute nothing of their own — the
only thing they see is the global fraction re-solve.

**How the fraction is re-solved.** V1's mixed residual `m1` takes the class rows from the raw CD4 object
and every other row from `lib246.deflate(cd4, "svd", 5)`; votes and counts are unchanged (`_as_res` copies
them). Its CD4 fraction is re-solved **once**, by production's own all-rows rule
(`plolo.solve_fraction(W, m1, 0.34)`) — a single scalar equal to 0.34 ÷ (median over every panel row of
‖m1‖ ÷ median ‖W‖). Because the un-deflated rows are longer the fraction falls: 0.395879 → 0.387298,
0.394269 → 0.385620, 0.397869 → 0.392070, a ratio of **0.978 / 0.978 / 0.985**. V3, which zeroes 84 rows,
goes the other way (**×1.129**: the all-rows median falls, so the surviving rows are up-dosed by 13 %) —
which is exactly why the fixed-fraction diagnostics exist. They separate the re-solve from the class edit
and change nothing that matters: V1 −0.00431 vs V1f −0.00435, V3 −0.00707 vs V3f −0.00684. **Delivered
amplitude** moves −0.06 % (V1), +0.44 % (V2), +0.14 % (V3), +0.53 % (V4); V2, V2f, V3f and V4 exceed 243's
0.2 % threshold and carry 197's amplitude leg in the price (V4: Δraw −0.00103).

---

## 2. The arm table (three cancer folds, vs `cd4svd5`; `kolf2` apart)

| arm | Δ`pds` | paired 95 % CI | folds > 0 | on edited rows | `cos_agg` | macro | Δ`nmae` | amp | LT3O | **PESS** | POINT | OPT | kolf2 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `V1_undefl` | −0.00431 | [−0.00917, +0.00055] | 1/3 | −0.01447 | −1.10 % | −0.31 % | +0.00011 | −0.06 % | +1.42 | **−0.00455** | −0.00430 | −0.00200 | −0.00066 |
| `V1f_undefl` | −0.00435 | [−0.00927, +0.00059] | 1/3 | −0.01466 | −0.86 % | −0.15 % | +0.00013 | −0.08 % | +1.42 | **−0.00402** | −0.00382 | −0.00188 | −0.00060 |
| `V2_mirror` | −0.00825 | [−0.01291, −0.00342] | 0/3 | −0.01171 (200 rows) | +1.56 % | +1.13 % | +0.00232 | +0.44 % | +1.19 | **−0.00365** | −0.00089 | +0.00057 | −0.00354 |
| `V2f_mirror` | −0.00860 | [−0.01336, −0.00371] | 0/3 | −0.01223 (200 rows) | +2.08 % | +1.44 % | +0.00260 | +0.43 % | +1.18 | **−0.00364** | −0.00002 | +0.00160 | −0.00347 |
| `V3_noresid` | −0.00707 | [−0.01209, −0.00189] | 0/3 | −0.02310 | +1.88 % | +0.12 % | +0.00041 | +0.14 % | +1.21 | **−0.00229** | +0.00071 | +0.00208 | −0.00428 |
| `V3f_noresid` | −0.00684 | [−0.01184, −0.00166] | 0/3 | −0.02309 | +0.13 % | −0.81 % | +0.00123 | +0.28 % | +1.22 | **−0.00342** | −0.00305 | −0.00209 | −0.00474 |
| **`V4_floor`** | **+0.00063** | [−0.00057, +0.00183] | **3/3** | +0.00184 | +1.52 % | −0.05 % | +0.00045 | +0.53 % | +0.30 | **+0.00078** | +0.00358 | +0.00390 | −0.00301 |
| `V5_20` (bottom 20 %) | −0.00392 | [−0.00834, +0.00054] | 1/3 | −0.01979 | −0.61 % | −0.32 % | +0.00011 | −0.14 % | +1.48 | **−0.00323** | −0.00309 | −0.00160 | +0.00059 |
| `V5_40` (bottom 40 %) | −0.00588 | [−0.01118, −0.00040] | 0/3 | −0.01480 | −1.41 % | −0.41 % | −0.00013 | +0.04 % | +1.31 | **−0.00594** | −0.00560 | −0.00260 | −0.00100 |
| `C_top30` (top 30 %) | −0.00311 | [−0.00588, −0.00043] | 0/3 | −0.01044 | +1.34 % | +0.60 % | +0.00110 | +0.07 % | +1.43 | **−0.00111** | +0.00114 | +0.00185 | −0.00107 |
| `C_agree30` (262's agreement) | −0.00464 | [−0.00901, −0.00058] | 0/3 | −0.01693 | −0.44 % | −0.26 % | +0.00148 | +0.15 % | +1.33 | **−0.00372** | −0.00342 | −0.00211 | +0.00394 |
| `C_agreeN` (agreement, size-matched) | −0.00387 | [−0.00833, +0.00035] | 0/3 | −0.01298 | −0.42 % | −0.23 % | +0.00146 | +0.13 % | +1.39 | **−0.00334** | −0.00305 | −0.00186 | +0.00410 |
| `C_rand1` | +0.00015 | [−0.00358, +0.00397] | 1/3 | +0.00040 | +0.72 % | +0.68 % | +0.00240 | +0.12 % | −13.5 | **−0.00055** | +0.00093 | +0.00109 | +0.00029 |
| `C_rand2` | −0.00349 | [−0.00813, +0.00109] | 0/3 | −0.01167 | −0.04 % | +0.14 % | +0.00015 | +0.17 % | +1.52 | **−0.00173** | −0.00170 | −0.00116 | −0.00329 |
| `C_rand3` | −0.00314 | [−0.00631, +0.00007] | 0/3 | −0.01076 | −0.30 % | −0.10 % | +0.00000 | −0.12 % | +1.41 | **−0.00211** | −0.00204 | −0.00115 | −0.00625 |

"on edited rows" is the mean Δ`pds` over the rows the arm actually changes (V2/V2f edit the *complement* of
the class, 200 rows). Leave-top-3-out is reported as the archive defines it (the fold-mean Δ after deleting
the three largest positive per-target Δ, divided by the full mean); on a **negative** arm the ratio exceeds
1 and carries no information — it is informative only for `V4_floor` (**+0.30**, i.e. three targets carry
70 % of its gain) and it is degenerate for `C_rand1`, whose full mean is ≈ 0.

**Pricing** (257's recipe, imported through `price272.price`, which itself known-answers against
`price257.price`): `pds` 2.2138 × k_pds × Δ`pds` / 6; `mse` through 197's identity at `cd4svd5`'s own
terms with 243's amplitude leg above 0.2 %; DE-three at a rate × Δcos_board; `nmae` = 188/216's exact
identity on the wald89 pool, the three folds concatenated. **Coherent** scenarios: per leg the end of its
band chosen by that leg's sign for that cell. V1–V3/V5/controls re-shape occupied CD4 rows — **the
deflation's own class** — so k_cos ∈ {1.36, 3.76} and k_pds ∈ {0.85, 1.22}; V4 changes the operator and
carries 272's band k_cos ∈ {0.94, 3.76}; DE-three ∈ {0.0425, 0.24} coherent with k_cos; the identity ×1.5
on a loss / ×0 on a gain (pessimistic), ×1.0 (point). POINT = the set that retrodicts BR249 (KA-4).

The six legs of V1 at the pessimistic set: `pds` **−0.00194**, `mse` −0.00152 (Δcos_board −0.00434),
`reach` −0.00075, `fid` −0.00026, `jac` −0.00003, `nmae` −0.00005 = **−0.00455**. Every member is priced;
none is unreadable.

**Cells clearing the pessimistic gate: 0. Cells clearing POINT: 1 (`V4_floor`, +0.00358) — and it fails
3/3-with-interval, retention and its own controls.**

### 2.1 The frozen CLEARED rule, condition by condition

| condition (DEFINITION-277 §5) | `V1_undefl` | `V4_floor` |
|---|---|---|
| pessimistic ≥ +0.0019 over `cd4svd5` | **NO** (−0.00455) | **NO** (+0.00078 = 0.41×) |
| `pds` > 0 on 3/3 cancer folds | **NO** (1/3) | yes (3/3) |
| paired interval above 0 | **NO** ([−0.00917, +0.00055]) | **NO** ([−0.00057, +0.00183]) |
| leave-top-3-out retention ≥ 50 % | n/a (negative arm) | **NO** (0.30) |
| every control below half the arm | n/a (negative arm); against the 24-draw null, pooled z **−0.53** | **NO** (shuffled +0.00102 and mirrored +0.00078 both ≥ +0.00063) |
| Flex-400 does not contradict the sign | **NO** (Flex +0.00336 vs PLOLO −0.00431) | not measured (V4 is an operator arm; the Flex read covers V1–V3) |

**DOES NOT CLEAR.**

---

## 3. What the class edit actually does

### 3.1 The anti-aligned gain is real and is outweighed by the marginal loss

262's groups, computed **on the edited rows only**, pooled over the three cancer folds:

| arm | group e (base·truth ≤ 0) | marginal (0, 0.05] | saturated > 0.05 |
|---|---:|---:|---:|
| `V1_undefl` | **+0.03384** (n 78) | **−0.03858** (n 163) | −0.00209 (n 12) |
| `V2_mirror` | −0.00220 (n 71) | −0.01731 (n 401) | +0.00006 (n 128) |
| `V3_noresid` | −0.05629 (n 78) | −0.00861 (n 163) | −0.00155 (n 12) |
| `C_rand1` | +0.05620 (n 48) | −0.01556 (n 167) | −0.00028 (n 38) |
| `C_rand2` | +0.03397 (n 45) | −0.02607 (n 173) | −0.00051 (n 35) |
| `C_rand3` | +0.02048 (n 49) | −0.02143 (n 176) | −0.00004 (n 28) |
| `C_top30` | −0.00268 (n 18) | −0.01732 (n 152) | +0.00008 (n 83) |
| `V4_floor` | +0.00498 (n 78) | +0.00054 (n 163) | −0.00052 (n 12) |

Reading. The class *is* enriched for the anti-aligned group — 78/253 edited rows (30.8 %) against
45–49/253 (18–19 %) for a random set, which is exactly the enrichment 276 §4.3 predicts. And un-deflating
*does* lift that group (+0.034). But (i) the random sets lift it just as much or more (+0.020 … +0.056), so
the lift is a property of the **group**, not of the map's selection; and (ii) the class's marginal rows —
twice as numerous — lose 1.8× more than a random set's marginal rows (−0.0386 vs −0.0156 … −0.0261), which
is the whole net loss. **V3 shows the same object from the other side**: deleting the residual on the class
costs −0.05629 on exactly those anti-aligned rows, i.e. the deflated CD4 residual is *already* worth
+0.056 there. There is nothing left for a different rendering to collect.

### 3.2 Nothing crosses the knee

| arm | crossings up | down | unchanged | mean cosine move | Δ`pds` up-crossers | Δ`pds` down-crossers |
|---|---:|---:|---:|---:|---:|---:|
| `V1_undefl` | 10 | 12 | 231 | −0.00001 | +0.116 | −0.107 |
| `V3_noresid` | 10 | 13 | 230 | −0.00088 | +0.081 | −0.164 |
| `C_rand1` | 6 | 7 | 240 | +0.00092 | +0.143 | −0.158 |
| `C_top30` | 3 | 10 | 240 | +0.00075 | +0.053 | −0.108 |
| `V4_floor` | 3 | 1 | 249 | +0.00030 | +0.020 | −0.029 |

The 0.02 knee behaves exactly as R273 §5 says it does — a crossing is worth ≈ ±0.11 `pds` for that target —
but the rendering produces as many downward crossings as upward ones, and the mean cosine move on the class
is **−0.00001**. The mandate's question ("can a different rendering move them across the 0.02 knee?") has a
measured answer: **it moves 22 of 253 class rows across the knee, 10 up and 12 down.**

### 3.3 The map *does* rank where un-deflating hurts — very weakly

Pooled by within-fold rank over 850 (target, fold) rows (276's convention),
ρ(`C_lolo`, V1's per-target Δ`pds`) = **+0.081 [+0.001, +0.157]**: low-conservation rows lose slightly more
from un-deflating. Per edited row, un-deflating the class costs **−0.01447** against **−0.01171** for the
complement (V2's own edited rows) — a 1.24× ratio. The pipeline's donor–donor agreement reads
ρ **+0.093 [+0.021, +0.167]** on the same axis, i.e. **at least as good as the map**. This is the same
verdict 276 §4.2/§4.6 reached on a different outcome.

---

## 4. The controls, and the selector's proper null

The frozen rule asks every control to sit below half the arm. With a **negative** arm that test is
degenerate, so the honest reading is the size-free one: does the map pick rows where V1's construction
behaves differently from an arbitrary row set of the same size? 24 random draws per fold, V1's identical
construction (`null277.py`):

| row set | k562 (z, pct of draws at or below) | hct116 | hek293t | pooled z |
|---|---|---|---|---:|
| `CLASS30` | −0.00469 (z −0.88, 8 %) | +0.00058 (z +0.83, 83 %) | −0.00883 (z −0.98, 21 %) | **−0.53** |
| `AGREE30` | −0.00509 (z −1.04, 4 %) | −0.00275 (z −0.27, 42 %) | −0.00608 (z −0.09, 42 %) | −0.73 |
| `AGREEN` (size-matched) | −0.00428 (z −0.71, 17 %) | −0.00198 (z −0.02, 42 %) | −0.00535 (z +0.14, 46 %) | −0.26 |
| `TOPCONS30` | −0.00171 (z +0.37, 67 %) | −0.00428 (z −0.77, 25 %) | −0.00335 (z +0.79, 79 %) | +0.20 |

Null mean ± sd of the fold-mean Δ`pds` over 24 draws: k562 −0.00258 ± 0.00240; hct116 −0.00193 ± 0.00304;
hek293t −0.00580 ± 0.00309. **The class is 0.53 σ below the random null pooled, with one fold on the other
side — not a selector for this effect**; and the mandate's stated requirement ("V1 must beat the
agreement-selected set to justify the map") reads the other way: the agreement set is *further* from the
null (−0.73) than the map's is.

The three seeds the frozen definition asked for, read as the definition asked (`C_rand1` +0.00015,
`C_rand2` −0.00349, `C_rand3` −0.00314; `C_agree30` −0.00464; `C_top30` −0.00311; mirror `V2` −0.00825),
say the same thing: **V1's −0.00431 sits inside the spread of its own controls.**

---

## 5. V4 — the one 3/3-positive cell, and why it is not the noise floor

V4 applies 272's A′ eigenvalue floor `max(λ, c·n)` at `c = 1` to the CLASS30 rows only and the incumbent
weights `λ^(−γ/2)` to every other row, with **one** global on-axis preserve; both residual fractions are
re-solved against the new base (263's `arm_parts`, which is what 272 did). The composition allows a per-row
operator choice cheaply — the expensive object is one matmul per weight vector, and any row set is then a
splice — so every cell below shares three matmuls per fold.

| cell | k562 | hct116 | hek293t | mean Δ`pds` | mean rel `cos_agg` |
|---|---:|---:|---:|---:|---:|
| **`floor:true:CLASS30` (= V4)** | +0.00008 | +0.00082 | +0.00098 | **+0.00063** | +1.52 % |
| `floor:shuffled:CLASS30` | +0.00099 | +0.00112 | +0.00096 | **+0.00102** | +1.06 % |
| `floor:mirrored:CLASS30` | +0.00044 | +0.00083 | +0.00107 | **+0.00078** | +1.10 % |
| `floor:true:TOPCONS30` | −0.00078 | +0.00019 | −0.00022 | −0.00027 | −0.50 % |
| `floor:true:RAND1 / 2 / 3` | +0.00043 / −0.00147 / −0.00129 | −0.00043 / +0.00062 / +0.00079 | −0.00045 / −0.00216 / −0.00128 | −0.00015 / −0.00100 / −0.00059 | ≈ +0.4 % |
| `floor:true:ALL` (= 272's `floor:c1`, KA at 0.0e+00) | −0.00449 | −0.00037 | −0.00169 | −0.00218 | +0.62 % |
| `gamma0.5:CLASS30` (POST-HOC) | +0.00024 | −0.00134 | −0.00216 | −0.00109 | +4.45 % |
| `gamma0.5:ALL` (POST-HOC) | −0.00578 | −0.00699 | −0.00811 | −0.00696 | +4.48 % |
| `gamma0.35:CLASS30` (POST-HOC) | −0.00135 | −0.00450 | −0.00650 | −0.00412 | +6.64 % |
| `gamma0.35:ALL` (POST-HOC) | −0.01608 | −0.02387 | −0.02298 | −0.02097 | +10.78 % |

Three readings.

1. **Restricting the floor to the class turns 272's global −0.00218 loss into a +0.00063 gain**, and the
   random row sets (mean −0.00058) and the top-conservation set (−0.00027) do not. That is the only place
   in this experiment where the map's row set does something the alternatives do not.
2. **But it is not the noise estimate.** Shuffling `n` across directions gives **+0.00102** and mirroring it
   gives **+0.00078** — both at or above the true vector. Whatever the class rows gain comes from capping
   the amplification of small eigenvalues *somewhere*, not from 272's measurement of where the noise is.
3. **And it is not "these rows want less whitening" either.** A lower γ on the class rows raises `cos_agg`
   (+4.5 % at γ 0.50, +6.6 % at γ 0.35) and **loses** `pds` (−0.00109, −0.00412) — 272 §0.10's mechanism
   holding on this subset: `pds` is bought by amplifying the noise-dominated directions, and damping them
   trades `pds` for aggregate cosine.

At +0.00078 pessimistic (0.41× the gate), an interval crossing zero, retention 0.30 and a kolf2 read of
−0.00301, V4 does not clear and is not a candidate. It is the one thread worth a sentence in the next
mandate, not a slot.

---

## 6. The Flex-400 re-read (h1-test, union OFF, 273's instrument, 85 targets)

Reference rows reproduce 273 exactly (KA-2 at 0.0e+00): BASE `cos_agg` 0.05727; the CD4 residual raises it
to 0.07464 un-deflated (**+30.3 %**) and 0.07107 deflated (**+24.1 %**) — R273/F4's +30 % / +24 %, the
observation that motivated V1. `pds` 0.871148 (base) / 0.882773 (un-deflated) / 0.883333 (incumbent).

| arm | Δ`pds` | paired 95 % CI | on class | off class | `cos_agg` | LT3O | PLOLO Δ`pds` |
|---|---:|---:|---:|---:|---:|---:|---:|
| `V1_undefl` | **+0.00336** | [−0.00448, +0.01232] | +0.01143 | +0.00000 | −0.29 % | −0.69 | −0.00431 |
| `V2_mirror` | −0.00350 | [−0.01232, +0.00449] | +0.00238 | −0.00595 | +5.31 % | +2.03 | −0.00825 |
| `V3_noresid` | −0.00238 | [−0.01120, +0.00491] | −0.01048 | +0.00099 | +0.43 % | +2.26 | −0.00707 |
| `V5_20` / `V5_40` | +0.00238 / +0.00238 | [−0.0049, +0.0108] / [−0.0059, +0.0116] | +0.00810 / +0.01143 | +0.00000 / −0.00139 | −0.27 % / −1.14 % | −1.16 / −1.34 | −0.00392 / −0.00588 |
| `C_top30` | −0.00392 | [−0.01092, +0.00084] | +0.00095 | −0.00595 | +6.06 % | +1.33 | −0.00311 |
| `C_agree30` | +0.00182 | [−0.00519, +0.00939] | +0.00571 | +0.00020 | −0.69 % | −1.36 | −0.00464 |
| `C_rand1 / 2 / 3` | +0.00126 / +0.00070 / −0.00028 | all crossing 0 | +0.0148 / −0.0081 / −0.0019 | −0.0044 / +0.0044 / +0.0004 | +0.04 / +2.51 / +0.19 % | −4.03 / −3.94 / +4.66 | +0.00015 / −0.00349 / −0.00314 |

**The Flex-400 reading is positive in sign and not established**, and it contradicts PLOLO. V1 reads
+0.00336 with an interval that crosses zero on both sides, `cos_agg` −0.29 %, and a leave-top-3-out of
−0.69 — i.e. deleting three of 85 targets flips the sign. Its whole move is 24 net rank-steps
(0.2857 total `pds` ÷ 1/84 per step) spread over 25 class rows, on a reference whose `pds` quantum is
1/84 = 0.0119. The random controls read +0.0013 / +0.0007 / −0.0003 and the agreement control +0.0018, so
V1 is at the top of its controls but inside their spread.

Under the frozen rule the Flex-400 condition is "does not contradict the sign" — here it *does* contradict
PLOLO's sign, but since PLOLO's own condition already fails by 2.4 gates, this is not what decides the
verdict. What it does establish is that **the observation that motivated V1 does not become a `pds` effect
even where it was measured**: on this reference the fully un-deflated CD4 residual carries a **+5.02 %**
`cos_agg` advantage over the deflated one (0.07464 vs 0.07107) and reads **−0.00056 `pds`** for it
(0.882773 vs 0.883333). That is R273 §0's finding again — every cosine effect in this pipeline is larger
than its `pds` consequence — now applied to the specific pair of renderings 276 §6.2 pointed at.

---

## 7. What this closes, and the one scientific residue

**Closed at this construction** (programme rule 6): *rendering the low-conservation class differently* — the
un-deflated CD4 residual on it, the mirror, its removal, the boundary at 20/30/40 %, and 272's eigenvalue
floor restricted to it — does not move `pds` on the three cancer PLOLO folds, is not distinguishable from
an arbitrary row set of the same size (pooled z −0.53 against a 24-draw null), and does not price above the
gate on any coherent scenario. Together with 262 (the selector axis), 263/264/265 (the arbiter axis),
271 (noise-aware combination) and 257 (the dose axis), **the CD4 residual's per-target axis is closed on all
four handles that were available: which rows, how much, from whom, and how rendered.**

276 §6.2's competition-facing framing is answered: **no.** 276's own recommendation (no dose, no arbiter,
extend on the depth axis rather than the context axis) stands unchanged and is now the only live branch.

**The residue that outlives the negative — and the correction it forces.** By conservation band, the k = 5
deflation's own per-target gain looks like a strong gradient:

| band | rows | mean `cd4svd5` `pds` | deflation Δ`pds` | share of its gain |
|---|---:|---:|---:|---:|
| bottom decile | 85 | +0.658 | +0.01904 | 15.6 % |
| 10–30 % | 170 | +0.696 | +0.01183 | 19.3 % |
| 30–50 % | 169 | +0.795 | +0.01577 | 25.6 % |
| 50–70 % | 171 | +0.831 | +0.00979 | 16.1 % |
| 70–90 % | 170 | +0.891 | +0.01359 | 22.2 % |
| top decile | 85 | +0.905 | +0.00145 | 1.2 % |

**That gradient is headroom, not biology.** Per target the map does not predict the deflation's gain at all:
ρ(`C_lolo`, deflation Δ`pds`) = **−0.021 [−0.092, +0.048]** on n = 850; controlling for the un-deflated
arm's own per-target `pds` — the room a row has left — it is **+0.059 [−0.011, +0.126]**; and
ρ(`C_lolo`, the un-deflated arm's own `pds`) is **+0.390 [+0.330, +0.448]**, which is 276's finding
(+0.387) and the whole of the apparent effect. The donor–donor agreement is equally uninformative here
(ρ −0.003; paired Δρ against the map −0.026 [−0.088, +0.039], P(≤ 0) 0.783). **Anyone reading a band table
of a per-target intervention on this pipeline must control for headroom; the deciles alone will manufacture
a 13× gradient out of a ceiling.** That is the reusable methodological point, and it applies to 262 §8's
band tables and to any future router that scores per-target actions.

**No preregistration draft.** Nothing cleared, so `PREREGISTRATION-DRAFT.md` was not written. For the
record, the production variable it would have needed is named and **not made**: a per-target boolean row
mask consumed by `deflate_cd4_residual(cd4, k, *, basis=None)` in `scripts/make_submission_pm_cd4.py:1993`
(which today deflates every carried row — `rows = np.flatnonzero(target_condition_count > 0)`, line 2021),
plus a frozen `C_lolo` table shipped as an artifact under `data/donors/` with a sha256 in its sibling json
and a builder-side assertion that the mask's target order is the day's panel order. Neither exists and
neither should be built on this evidence.

---

## 8. Files

* `RESULTS.md` (this file); `work/DEFINITION-277.md` (frozen, sha `05a89391…`) and `work/DEFINITION-277.sha256`.
* `lib277.py` — the class selector (`c_lolo`, `bottom_class`, `random_class`), `mix_rows` / `zero_rows`
  (a CD4 residual assembled from two prepared residuals by a row mask), `mixed_whiten` (per-row operator
  weights with production's single on-axis preserve), the groups and the leave-top-3-out.
* `render277.py` → `work/render277-<fold>.{npz,json,log}` — V0–V5 and the controls on the four PLOLO folds,
  with KA-1/KA-1b/KA-3.
* `flex277.py` → `work/flex277.{npz,json}` — the same arms on 273's Flex-400 instrument, with KA-2.
* `price277.py` → `work/price277.{json,log}` — the six-member price on three coherent scenarios, with the
  BR249 retrodiction (KA-4).
* `null277.py` → `work/null277-<fold>.{json,log}` — 24 size-matched random row sets per fold, V1's construction.
* `floor277.py` → `work/floor277-<fold>.{npz,json,log}` — V4 by row set and by noise vector, plus the
  POST-HOC γ arms; KA against 272's `floor:c1`.
* `analyse277.py` → `work/analyse277.{json,log}`; `tables277.py` → `work/tables277.md` (every table above,
  generated from the json so no number is typed twice).
* `work/plolo/` — byte copies of 262's fold caches, `shasum -c` 7/7 OK.

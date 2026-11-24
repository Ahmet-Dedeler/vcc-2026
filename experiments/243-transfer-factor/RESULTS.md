# 243 — The PLOLO→board cosine factor is **not a scalar**, and the reason is neither the fold set, nor the functional, nor amplitude. It is the **arm class**: the two CD4 dose steps agree at **1.36 ± 0.05** and the union step sits at **0.94–1.00**, on the measured (amplitude-corrected) board cosine. ARM6 is re-priced at **+0.0022 … +0.0038** on the direction legs; the board raw `mse` that decides the class model is pre-computed.

**Date:** 2026-09-04 · read-only everywhere except `experiments/243-transfer-factor/` · **no `vcc`,
no build, no `cell-eval2` run, no production change, no slot proposal.** Everything under
`uv run python`, strictly serial. **Peak RSS 0.77 GB** (`/usr/bin/time -l`, the ARM6 stream); the
three h5ad streams peaked at 0.72 / 0.68 / 0.77 GB and the inversion at 0.63 GB; every other script
< 0.1 GB. No PLOLO fold was loaded — every PLOLO number is read from a sibling's artifact and its path
is given. Memory gate before each stream (not required for a stream, recorded anyway): before the first
stream it was read in-session and printed (85% free, 35.5 GiB free+inactive; the redirect to a
file was skipped by a failed `cd` and is recorded here instead), before the second it is in
`work/memgate-stream2.log` (86% free, 34.4 GiB). Two functions of the `cell_eval2` *library*
(`_row_weights`, `_across_pert_budget`) are imported exactly as 197 and 211 did; the scorer CLI was
never invoked.

**Known-answer rows: 20 of 20 pass** (`ka243.py` → `work/ka243.log`), plus 24 artifact
cross-checks at 0.0e+00 (`matrix243.py`), the `b_ctrl` identity (096 vs 162, 0.0e+00), the
`cd4f034` stream reproducing 214's cached totals on **900/900 blocks**, and the `rc_callset` /
`cd4f034` identity terms reproducing 211's `bx211.json` and 197's a15 terms to 1e-9.

---

## 0. The answer in ten lines

1. **KA-1/2/3 reproduce.** 216's +17.98% and k_cos 1.384 (on its own rounded constants and
   transcribed raws), 222's +10.3065% / +7.6803%, 197's 0.079248 on `rc_hybrid` and the two
   floors 0.09027 / 0.09947 — all to the quoted digits. §1.
2. **Amplitude does not explain it.** With each arm's own delivered energy `A_c` and budget `B_c`
   measured on the shipped h5ads (three streams, ~13 s each), 197's identity gives the board
   cosine as a *measurement*: `rc_callset` 0.079249 → `cd4f034` **0.097850** → `h1union`
   **0.107332**, i.e. **+23.47%** and **+9.69%** against the floor's +24.88% / +10.31%. The
   amplitude/budget leg moved raw `mse` by −0.0015 on the dose step and +0.00005 on the union
   step; delivered amplitude changed −0.9% and +0.2%. The two board cosines are comparable
   objects. §3.
3. **No convention closes the gap.** Over 48 conventions (2 functionals × 2 tags × 2 fold sets ×
   2 aggregations × 3 board readings) **zero** admit one factor within ±15% of all three board
   points. The closest is +16%. §4.1.
4. **The fold set and the tag are not the explanation** (222 §4's hypotheses 1 and 2): dropping
   `kolf2` changes `h1union`'s factor 0.52–0.74 → 0.94–1.00 and `cd4f034`'s 1.31 → 1.41, still
   apart; the a15/a23 tag moves them by ≤ 0.06. §2, §4.1.
5. **The functional is decided by within-class consistency, not by taste.** On `cos_agg` the two
   CD4 dose points (a14→a15 and `cd4f034`, same knob) agree at **1.315 / 1.405** (±3.4%); on the
   macro cosine they read **1.41 / 1.81** (±13%). `cos_agg` — 197's own functional — is the one on
   which one knob gives one number. §4.1.
6. **The smallest model that retrodicts all three: two classes.** On {`cos_agg`, three cancer
   folds, mean of per-fold ratios, measured board cosine}: **k_dose = 1.36** (a23 tag; 1.37 a15)
   for a re-dose of rows CD4 already occupies (276/300), **k_newsource = 0.94** (a23; 1.00 a15)
   for a first source into 77 empty rows. The equivalent coverage-linear form is
   `k = 0.78 + 0.63·coverage`. Both are 2-parameter fits to 3 points; they differ only in what
   they predict next. §4.2.
7. **`pds` is not a second equation for the cosine.** 215's amplitude-free route reads the CD4
   dose step at **−7.2%** (ρ = −0.32) because `pds` turns over at `f ≈ 0.2` (216 §4.1), and the
   union step at +5.4% against a measured +9.7%. The artifact's own `A_c` is the right second
   equation, and it is cheap. §3.3.
8. **`k_pds` is class-dependent the same way** on the 3-cancer convention (1.15 / 0.95 / 0.70) and
   wildly so on 4 folds (1.19 / 2.02 / 0.38); the DE-four propagation per unit board cosine is
   0.17 / 0.42 / 0.24. §4.3.
9. **ARM6 re-priced with its OWN delivered terms** (streamed from `229-stack6-g035.h5ad`: amplitude
   **−1.24%** vs `h1union`, worth −0.00044 overall on its own): new-source class **+0.0024**
   (band +0.0022 … +0.0038), coverage model +0.0029, dose class +0.0044 — every reading ≥ +0.0015.
   229's six-member −0.00219 differs from these by exactly its `nmae` identity leg (−0.0048).
   **The board raw `mse` separates the models: ≥ 0.81135 → new-source class; 0.81002–0.81135 →
   coverage; ≤ 0.81002 → dose class; ≥ 0.81392 → below every model.** §5.
10. **The falsifier fires tomorrow, and the floor must not be used to read it**: at ARM6's
    amplitude the floor would report k = 0.64 for a true 0.94. Read the measured column of
    `work/arm6_243.json`. §5.2.

---

## 1. Known-answer rows (`ka243.py`, `work/ka243.log`) — 20/20

| # | row | got | want | source |
|---|---|---:|---:|---|
| KA-1a | PLOLO `cos_agg` f 0.10→0.34, 4 folds, ratio of fold-means | +17.9845% | +17.9845% | `216/work/board216.json` `rel_plolo`; fold rows `216/work/sweep216.json` |
| KA-1b/c | 4-fold mean `cos_agg` at f = 0.10 / 0.34 | 0.0200856 / 0.0236979 | 0.020086 / 0.023698 | 216 §4.4 |
| KA-1d–f | board floor at raw 0.87367 / 0.83635 → +24.8898% | 0.0722783 / 0.0902682 | `board216.json` `c_ship`/`c_inc`/`rel_board` | 216's constants `1.023611 / 2.074495` (line 47) |
| KA-1g | **k_cos** | **1.383957** | 1.383957 | 216 §4.4 "1.384" |
| KA-2a/b | 222 rel `cos_agg` / macro, 3 cancer folds, mean of per-fold ratios | **+10.3065% / +7.6803%** | +10.3065% / +7.6803% | `222/work/repro222.json` `summary` |
| KA-2c | 222 Δ`pds`, 3 cancer folds | +0.004887 | +0.004887 | same |
| KA-3a/b | 197's two constants from (Kc̄, ΣD, K) | 1.023611086 / 2.074494532 | `197/work/invert197.json` | exact |
| KA-3c/d | `rc_hybrid` raw 0.8738045055 → **measured** cosine | **0.079248508** | `197/work/rchybrid197.json` `rc_cos` | 0.0e+00 |
| KA-3e/f | floor at 0.83635 → 0.09027; at 0.81726 → 0.09947 | 0.0902683 / 0.0994705 | mandate | to the quoted digits |

Two things the rows record that the prose did not. (i) 216's k_cos was computed on **rounded**
constants and **transcribed** raws (`rc_callset` 0.87367 for a true 0.8738045, `cd4f034` 0.83635
for 0.8365398; 219 §6 already flagged the first); on full-precision constants and exact raws it is
1.3831. Immaterial, recorded. (ii) The mandate's KA-3 mixes two different inversions: `rc_hybrid`'s
0.079248 is the **full** identity with a15's own `A`, `B`, `claim`; `cd4f034`'s 0.09027 and
`h1union`'s 0.09947 are the one-sided **floor**. That distinction is the whole of §3.

**The exact board rows used everywhere below** (leaderboard snapshots, entry ids in
`lib243.py`): a14 raw `mse` 0.8827752792, `pds` 0.8327015979; `rc_callset` (= a15's profile row to
16 s.f., 202/197) 0.8738045055 / 0.8359866221; `cd4f034` 0.8365398403 / 0.8322779636; `h1union`
0.8172630319 / 0.8357227796.

---

## 2. Phase 1 — both arms on the same instrument, every way it can be defined (`matrix243.py`, `work/matrix243.json`)

Sources: a14→a15 from `experiments/207-cd4-dose/work/sweep207.json` (arms `a14`, `f0.1`);
`cd4f034` from `experiments/216-cd4-hard/work/sweep216.json` (`f0.1`, `f0.34`); `h1union` from
`experiments/218-independent-donors/work/compose218.json` (`cd40.34+none0`,
`cd40.34+h1PAIRCARRY0.28`, both tags) which 222's production-builder run
(`repro222.json`) reproduces at 0.0e+00. 24 cross-artifact checks (207 vs 216 rows, 218 base vs
216 f0.34, 222 vs 218) all at **0.0e+00**.

**Relative PLOLO gain (%)** — functional × tag × fold set × aggregation:

| step | `cos_agg` a23 4/pooled | 4/perfold | **3c/pooled** | **3c/perfold** | `cos_agg` a15 3c/perfold | macro a23 3c/perfold | macro a15 3c/perfold |
|---|---:|---:|---:|---:|---:|---:|---:|
| a14→a15 | 5.42 | 5.48 | 5.33 | **5.20** | 5.18 | 4.85 | 4.88 |
| `cd4f034` | **17.98** ← 216 | 18.86 | 17.19 | **16.70** | 16.57 | 12.94 | 13.05 |
| `h1union` | 13.96 | 18.80 | 10.47 | **10.31** ← 222 | 9.64 | **7.68** ← 222 | 7.28 |

Per fold (`cos_agg`, a23): a14→a15 6.18 / 5.57 / 3.85 / **6.32** (kolf2); `cd4f034` 20.18 / 18.16 /
11.77 / **25.34**; `h1union` 10.91 / 11.51 / 8.50 / **44.27**. `kolf2` inflates the union step by
4× more than it inflates the dose steps, which is why the fold set moves `h1union`'s factor most.

**Δ`pds`** (a23): a14→a15 +0.002761 (4) / +0.002866 (3c); `cd4f034` −0.001837 / −0.003916;
`h1union` +0.009072 / +0.004887.

**Board side, three readings** (§3 for the measured one):

| step | floor rel | **measured rel** | 215 `pds`-route rel |
|---|---:|---:|---:|
| a14→a15 | +6.37% | **+6.83%** | +6.83% (calibration point) |
| `cd4f034` | +24.88% | **+23.47%** | −7.22% |
| `h1union` | +10.31% | **+9.69%** | +5.43% |

**The factor, on the reference convention** (`cos_agg`, 3 cancer folds, per-fold mean):

| board reading | a14→a15 | `cd4f034` | `h1union` |
|---|---:|---:|---:|
| floor (222/216's) | 1.225 | 1.489 | 1.000 |
| **measured** (a23 tag) | **1.315** | **1.405** | **0.940** |
| measured (a15 tag) | 1.320 | 1.417 | 1.005 |

The full 48-row table is in `work/rule243.log`. 216's 1.384 and 222's 1.000 are both reproduced
and both sit on *different* conventions (4-fold pooled floor vs 3-cancer per-fold floor); put on the
same convention they become 1.49 vs 1.00 (floor) or 1.41 vs 0.94 (measured). **The discrepancy is
not a convention artefact; it survives every convention.**

---

## 3. Phase 2 — the board side is an inversion, and its slack is measured (`invert243.py`, `work/invert243.json`)

### 3.1 The decomposition

197 §1, per context: `MSE_c − 1 = (A_c + Kc_c − min(claim_c, B_c) − 2 X_c) / ΣD`, `raw = mean_c MSE_c`,
with `A_c = Σ_p w_p ‖δ̂_p‖²` (delivered energy, own target dropped), `B_c` the scorer's
`_across_pert_budget` on the delivered pseudobulk, `Kc_c` the reference's jackknife term, `claim_c`
capped at `cap_c = (Σ w_p) C_p`, `X_c` the alignment, `ΣD = 0.229141` the truth energy. Writing
`r_c = √(A_c/ΣD)` and `c = X̄ / (√ΣD · mean_c √A_c)` (216/197's functional):

```
raw = 1 + mean_c[ r_c² − (min(claim_c, B_c) − Kc_c)/ΣD ] − 2 c · mean_c √A_c / √ΣD
```

which is 197's `raw = r² − 2rc + 1` plus the budget/jackknife leg. The one-sided **floor**
`cos ≥ (1.023611 − raw)/2.074495` replaces the arm's own (`A`, `B`, `claim`) by the reference
constants, so **between two arms it attributes every change in raw to `c`**. Whether that is safe
depends on how much `A_c` and `B_c` moved — which the metas cannot say (the preserve scalar fixes a
*median row norm*, `A` is an *energy*, and the union loads 77 heavy rows) but the shipped h5ads can.

### 3.2 Measured on the shipped artifacts

Per-block panel totals streamed from `submissions/210-cd4f034.h5ad` and
`submissions/222-h1union-g028.h5ad` with 214's `block_totals` imported by value (`totals243.py`; the
`cd4f034` stream equals 214's cache on **900/900 blocks**, nnz 1,577,729,493 = the meta's; `h1union`
nnz 1,577,411,397 = its meta's). `rc_callset` from 214's cache. Terms by 211's `bx211.py`
convention, reproducing `bx211.json` and 197's a15 terms to 1e-9. `claim` at the cap (211
`chat211`: 0 rows under the cap on all three contexts for this emitter) — and it does not matter:
**all three contexts of `cd4f034` and `h1union` are on the `B` leg** (`B_c < cap_c`), so the result
is independent of `claim` unless `claim` fell below `B`, i.e. by > 1.2% in A, against an emitter
whose per-row `Ĉ` is 2.6× the cap.

| arm | raw `mse` | `A_A` / `A_B` / `A_C` | B/A | legs (A,B,C) | `X̄` | **measured cos** | floor | floor/meas | mean `r` |
|---|---:|---|---:|---|---:|---:|---:|---:|---:|
| `rc_callset` (= a15) | 0.8738045 | 0.23643 / 0.22643 / 0.23662 | 0.9909 | claim, B, B | 0.018317 | **0.079249** | 0.072214 | 0.911 | 1.0087 |
| `cd4f034` | 0.8365398 | 0.23291 / 0.22189 / 0.23205 | 0.9915 | B, B, B | 0.022411 | **0.097850** | 0.090177 | 0.922 | 0.9995 |
| `h1union` | 0.8172630 | 0.23363 / 0.22271 / 0.23290 | 0.9914 | B, B, B | 0.024625 | **0.107332** | 0.099469 | 0.927 | 1.0013 |

(a14 from 197's own stored terms: 0.074179; a15 from 197's terms 0.079248 vs `rc_callset` above
0.079249 — |d| 2.5e−7, the same object twice.)

### 3.3 Each step, split

| step | Δraw | Δraw from amplitude/budget leg | Δraw from alignment leg | delivered amplitude √(A_b/A_a)−1 | floor Δcos (rel) | **measured Δcos (rel)** | 215 `pds` route (rel, ρ) |
|---|---:|---:|---:|---:|---:|---:|---:|
| a14→a15 | −0.008971 | **+0.000920** | −0.009890 | −0.23% | +0.004324 (+6.37%) | **+0.005069 (+6.83%)** | +0.005069 (+6.83%, ρ 1.17) |
| `cd4f034` | −0.037265 | **−0.001532** | −0.035733 | −0.91% | +0.017963 (+24.88%) | **+0.018601 (+23.47%)** | −0.005723 (−7.22%, ρ −0.32) |
| `h1union` | −0.019277 | **+0.000053** | −0.019330 | +0.17% | +0.009292 (+10.31%) | **+0.009483 (+9.69%)** | +0.005315 (+5.43%, ρ 0.57) |

In floor-cosine units the amplitude/budget leg is worth **+0.00074 of the dose step's +0.01796
(4.1%)** and **−0.00003 of the union step's +0.00929 (0.3%)**. So the inversion's slack does differ
between a dose step and a new-source step — by an order of magnitude — but it is small on both, and
it runs the *wrong way* to reconcile the factors: correcting it takes the dose step's relative gain
*down* (24.9 → 23.5%) and the union step's *down* too (10.3 → 9.7%). **The two board cosines are
comparable objects, and the 1.4-vs-1.0 gap is not amplitude.**

**Sensitivity, per arm** (`d floor / d ln s` for a global rescale `s` of the delivered effect): a
**+1% amplitude** moves the floor cosine by **−1.73% of cos** on `rc_callset` (A on the `claim` leg,
0.12% above its cliff) but by **+0.88%** on `cd4f034` and **+0.73%** on `h1union` (all contexts on the
`B` leg, below the cliff). The metas' preserve scalar moved **−1.195%** (0.920359 → 0.909357) while
the measured delivered energy moved **+0.17%**: the median-norm preserve is not an amplitude
readout for a union that loads 77 rows, and the metas alone would have mis-signed this leg.

**The `pds` second equation.** 215's `Δcos = Δ(raw pds)/0.6479` reproduces a14→a15 exactly because
that is where 0.6479 was fitted. It reads the dose step at **−7.2%** (ρ = −0.32, 215's GEOMETRY
band) for an arm whose measured cosine rose 23.5% — because PLOLO `pds` peaks at `f ≈ 0.175` and
falls past it (216 §4.1) while the cosine keeps rising to `f ≈ 1.0`; and it reads the union step at
+5.4% against +9.7%. `pds` is a rank statistic with its own dose response; it is not a cosine
readout for either arm class. **The artifact's own `A_c` is the amplitude-free second equation, it
costs 13 s per arm, and every future direction arm should be inverted with it before its board
Δcos is quoted.**

**By-product, flagged for 211/187's owner, not acted on:** `cd4f034` and `h1union` sit *below* the
amplitude cliff in all three contexts (headroom `√(cap/B) − 1`: A **+0.45–0.60%**, B +4.2–4.4%, C
+7.4–7.6%; ARM6 +1.7 / +5.5 / +8.8%), whereas `rc_callset` was 0.12% *above* it in A. 211's
"claim at 100% of cap, headroom −0.124%" closure was measured on a15's profile. The A-context pot is
tiny (~+0.0001 overall); B/C read larger on the proxy cap but 197 §1.1 showed the proxy over-reads
C's true cliff by 1.75% and `abc-knee` paid for exactly that. Recorded, not re-opened.

---

## 4. Phase 3 — the rule (`rule243.py`, `work/rule243.{json,log}`)

### 4.1 One factor, any convention: **not identified**

For each of the 48 conventions the single factor `k*` minimising the worst |deviation| over the
three board points, and that worst deviation:

| convention (functional / tag / folds / agg / reading) | a14→a15 | `cd4f034` | `h1union` | `k*` | worst |
|---|---:|---:|---:|---:|---:|
| `cos_agg` a23 3c perfold **measured** (reference) | 1.315 | 1.405 | 0.940 | 1.149 | +22% |
| `cos_agg` a15 3c perfold measured | 1.320 | 1.417 | 1.005 | 1.193 | +19% |
| `cos_agg` a23 3c perfold floor (222's) | 1.225 | 1.489 | 1.000 | 1.220 | +22% |
| `cos_agg` a23 4 pooled floor (216's) | 1.175 | 1.383 | 0.738 | 1.010 | +37% |
| macro a23 3c perfold measured | 1.410 | 1.814 | 1.262 | 1.513 | +20% |
| macro a15 3c pooled measured (the closest of all 48) | 1.376 | 1.756 | 1.306 | 1.514 | **+16%** |
| any convention with the `pds` route | — | negative | — | n/a | n/a |

**0 of 48 pass ±15%.** The closest (macro, a15, pooled) passes only by spreading the two CD4
points to 1.38 and 1.76 — the same knob, 28% apart — which is a worse instrument, not a better rule.
On `cos_agg` the two CD4 points agree within **±3.4%** (a23) / ±3.6% (a15) on either fold set and
either aggregation; that is what fixes the functional. On the fold set: `kolf2` in the set moves
`h1union` from 0.94 to 0.52–0.74 and `cd4f034` from 1.41 to 1.25–1.31; out is the honest choice
for a pluripotent source scored on a pluripotent fold (218/222/228 all excluded it) and it is also
the set on which the dose points agree best.

### 4.2 The smallest model that retrodicts all three: two classes

On the reference convention — **`cos_agg` (a23 tag), three cancer folds, mean of per-fold ratios,
board cosine MEASURED with the arm's own `A_c`/`B_c`**:

| class | what the arm does | board points | **k_cos** |
|---|---|---|---:|
| **dose** | re-doses a residual on rows it already occupies (CD4, 276/300) | a14→a15 1.315, `cd4f034` 1.405 | **1.36** (±3.4%; a15 tag 1.37) |
| **new source** | first source into rows that were empty (H1 union, 77/300) | `h1union` 0.940 | **0.94** (a15 tag 1.00) |

Equivalent 2-parameter coverage form: `k = 0.778 + 0.632 · coverage` (a23; a15: `0.865 + 0.547·cov`)
— exact on the two coverage values by construction, so **the data cannot tell class from coverage;
only the next point can** (§5). This is 221 §9.1's hypothesis 3 and its coverage suggestion, now
with the amplitude confound removed. An alternative mechanism with the same ARM6 prediction — CD4
T-cell data transferring especially well into the T-lineage context A (JURKAT), so that the dose
class is really a *source-lineage* effect — cannot be separated from stack position by any board
point available before Oct 22, and is noted only because it would change what CD4 is worth on D/E/F.

### 4.3 The other legs are class-dependent the same way

| step | board Δ`pds` | `k_pds` (3c) | `k_pds` (4 folds) | DE-four Δoverall | per unit measured Δcos | `mse` Δoverall per unit Δcos |
|---|---:|---:|---:|---:|---:|---:|
| a14→a15 | +0.003285 | 1.146 | 1.190 | +0.000879 (219) | 0.173 | 0.310 |
| `cd4f034` | −0.003709 | 0.947 | **2.019** | +0.007864 | 0.423 | 0.351 |
| `h1union` | +0.003445 | 0.705 | 0.380 | +0.002267 | 0.239 | 0.357 |

Three things: the 3-cancer fold set also makes `k_pds` sane (216's 1.987 is a 4-fold number and
becomes 0.95 without `kolf2`); the `mse` exchange rate is the one constant in the table
(0.31–0.36 per unit cosine, 219's +0.3655); and the DE-four propagation for a *new source* (0.24)
sits near a14→a15's (0.17), not `cd4f034`'s (0.42) — the `rc_callset` call set does not by itself
double the DE conversion, as 219 §0.4 guessed; the dose to `f = 0.34` did.

### 4.4 The rule, stated

```
Quote a direction arm in PLOLO units as: rel cos_agg, a23 tag, three cancer folds, mean of
per-fold ratios.  Convert to the board with the CLASS factor, not a scalar:
    re-dose of occupied rows  (CD4 f-steps)      k_cos 1.36  (band 1.31-1.42)   k_pds 0.95-1.15
    first source into empty rows (a union)       k_cos 0.94  (band 0.94-1.00)   k_pds 0.70 (0.38-1.15)
Read the board with the MEASURED cosine (the arm's own A_c, B_c streamed off its h5ad, 13 s),
never the floor: the floor mis-reads an amplitude change of 1% as 0.7-1.7% of cosine.
mse leg: exact identity at the arm's own terms (d raw/d cos = -2.00 at h1union's amplitude);
DE-four: 0.17-0.42 x board d cos (0.24 for the union point); nmae is NOT priced by this rule.
```

---

## 5. Re-pricing ARM6 (237) and 228's per-gene arm

### 5.1 ARM6 — `229-stack6-g035`, `cd4@0.34 + union(six sources)@0.35`, coverage 176/300

PLOLO (229 `price229.json`, per fold vs `h1union`'s per-fold base in `compose218.json`): rel
`cos_agg` **+3.673%** on the three cancer folds (1.97 / 3.41 / 5.64; `kolf2` +3.85), rel macro
+2.865% (**= 229 §6.1's +2.865%**, known-answer), Δ`pds` +0.002687 (229: +0.00269).

**ARM6's own identity terms**, streamed from `submissions/229-stack6-g035.h5ad` (nnz 1,578,196,231
= its meta's; `work/totals-ARM6.npz`, `work/arm6_243.json`): `A` 0.22806 / 0.21734 / 0.22682, all
contexts on the `B` leg, **delivered amplitude −1.24% vs `h1union`**. On the `B` leg that costs raw
`mse` **+0.00250 at zero cosine change = −0.00044 overall** — a charge no PLOLO reading can see and
that the rule's price now carries.

| reading of the rule | k_cos | k_pds | DE-four rate | board Δcos | **predicted raw `mse`** | o_mse | o_pds | o_DE4 | **TOTAL** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **new-source class, point** | 0.94 | 0.70 | 0.24 | +0.003706 | **0.81244** | +0.00085 | +0.00070 | +0.00089 | **+0.00243** |
| new-source, a15-tag k | 1.00 | 0.70 | 0.24 | +0.003962 | 0.81193 | +0.00094 | +0.00070 | +0.00095 | +0.00258 |
| new-source, LOW band | 0.94 | 0.70 | 0.17 | +0.003706 | 0.81244 | +0.00085 | +0.00070 | +0.00064 | +0.00219 |
| new-source, HIGH band | 1.00 | 1.15 | 0.42 | +0.003962 | 0.81193 | +0.00094 | +0.00114 | +0.00168 | +0.00375 |
| coverage-linear (cov 0.587 → k 1.15) | 1.15 | 0.70 | 0.24 | +0.004528 | 0.81081 | +0.00113 | +0.00070 | +0.00108 | +0.00291 |
| dose class (union point was the outlier) | 1.36 | 0.70 | 0.42 | +0.005358 | 0.80917 | +0.00142 | +0.00070 | +0.00227 | +0.00438 |

Against: **229's six-member −0.00219**, whose `mse`+`pds` legs are +0.00031 + 0.00200 (229 §6.2) and
whose `nmae` leg is **−0.004783** from 216's identity instrument; and 237's pre-committed bins at
**±0.0015**. The rule's direction legs alone put ARM6 at **≥ +0.0022 under every model**, so:

- **Δ ≤ −0.0015** ⇒ the identity's `nmae` leg is real *and* at least −0.004 — the gate is
  validated at scale exactly as 237 says, and this rule's direction legs are confirmed as an
  upper bound on the arm.
- **Δ ≥ +0.0015** ⇒ the direction legs account for the whole result and the identity instrument
  mis-signed `nmae` for the third time in three.
- The two are **~0.0045 apart**, above 230's 0.0045 ordering-channel sd only just; the **`mse` leg
  is the clean part** (the board is deterministic to 1e−7 on it) and its prediction is
  **raw `mse` 0.8119–0.8124** under the class model.

### 5.2 The falsifier, pre-computed so tomorrow's raw `mse` reads directly

`work/arm6_243.json` converts any board raw `mse` into ARM6's **measured** cosine
(`X̄ = [mean_c(A_c + Kc_c − B_c) − (raw − 1)ΣD]/2`, `cos = X̄/(√ΣD · mean_c √A_c)` with ARM6's own
terms) and the realised factor `k = Δcos / (3.673% × 0.107332)`:

| board raw `mse` | measured cos | Δcos vs `h1union` | **k realised** | reading |
|---:|---:|---:|---:|---|
| ≤ 0.80768 | ≥ 0.11344 | ≥ +0.0061 | ≥ 1.55 | above every model |
| 0.80768 – 0.81002 | | | 1.25 – 1.55 | **dose class** — the union point was the outlier |
| 0.81002 – 0.81135 | | | 1.08 – 1.25 | **coverage-linear** |
| **0.81135 – 0.81392** | 0.1096 – 0.1113 | +0.0023 – +0.0040 | **0.75 – 1.08** | **new-source class** (point 0.81244) |
| ≥ 0.81392 | ≤ 0.10962 | ≤ +0.0023 | < 0.75 | below every model — the union factor itself was optimistic |

**Do not read this with the floor.** At the class prediction (raw 0.81244) the floor reports
k = **0.64** for a true 0.94: ARM6's −1.24% amplitude costs **0.30 units of k** in the floor
reading. The measured column is the instrument.

### 5.3 228's per-gene gain (`gain_reg@0.25`), coverage 300/300, occupied rows re-shaped

PLOLO (228 `analyse228.json` / `price228.json`): rel `cos_agg` **+4.328%** at production α
(+4.539% @adel), Δ`pds` +0.004704 (+0.004841), three cancer folds — 228 §3/§7's own numbers.

| class | k_cos | k_pds | board Δcos | raw `mse` | o_mse | o_pds | **mse + pds** |
|---|---:|---:|---:|---:|---:|---:|---:|
| dose / occupied rows, production α | 1.36 | 1.04 | +0.006314 | 0.80462 | +0.00222 | +0.00181 | **+0.00402** |
| dose / occupied rows, @adel | 1.36 | 1.04 | +0.006622 | 0.80400 | +0.00233 | +0.00186 | +0.00419 |
| new-source (lower bound) | 0.94 | 0.70 | +0.004368 | 0.80852 | +0.00153 | +0.00122 | +0.00276 |
| 228's own, k_cos 1.384 | — | — | — | — | +0.00227 | +0.00180 | +0.00408 |

The rule prices 228's direction legs at **+0.0028 … +0.0042** — it changes nothing about that arm,
because the rule's DE-four leg is a *sign-accuracy* propagation (§8.4) and 228 measured its arm
moving |t|-weighted sign accuracy by +0.00027 while moving `nmae`'s unmatched magnitude by +0.0373
raw (`e` = 0.072). **This rule does not price `nmae`, and `nmae` is what decides 228.** What ARM6's
board point does for 228: it tests the same identity instrument at a raw Δ`nmae` of +0.0172; 228's
arm sits at +0.0373, 2.2× further out. If ARM6 lands ≥ +0.0015, 228's `nmae` leg must be re-derived
before the arm can be read at all; if ≤ −0.0015, 228 stays dead a fortiori.

---

## 6. What this settles, and what it does not

**Settled.** (1) The 1.384 / 1.000 pair is a real, convention-independent, amplitude-independent
difference between two arm classes. (2) `cos_agg` on the three cancer folds is the PLOLO convention
on which one knob gives one factor. (3) The board cosine of any arm whose package is on disk is a
13-second measurement, not a floor, and the floor's amplitude slack (0.7–1.7% of cosine per 1% of
delivered amplitude) is large enough to flip a reading — it would have flipped ARM6's. (4) `pds` is
not a second equation for the cosine off the a14→a15 calibration point.

**Not settled, and the next board point decides it.** Whether the factor is a *class* constant or a
*coverage* function (ARM6 discriminates them at raw `mse` 0.81135 vs 0.81002–0.81135), and whether
CD4's 1.36 is stack position or T-lineage transfer into context A (nothing before Oct 22 can tell).

**Not attempted.** No `nmae` pricing — the mandate's rule is a cosine rule, and 237 exists because
the `nmae` instrument is the open question. No PLOLO fold load. No re-derivation of a number that
already existed.

---

## 7. Summary for `AGENTS.md` (≤ 10 lines, for the lead to paste)

```
243 (2026-09-04): the PLOLO->board cosine factor is NOT a scalar and the discrepancy is real.
- Both arms measured on the same instrument every way (2 functionals x 2 tags x 2 fold sets x 2
  aggregations x 3 board readings): 0 of 48 conventions admit one factor within +-15%.
- Amplitude is excluded: streaming each arm's own A_c/B_c off its h5ad (13 s) and inverting 197's
  identity gives MEASURED board cosines 0.0792 -> 0.0979 -> 0.1073 (+23.5%, +9.7%); the
  amplitude/budget leg is 4% of the dose step's dcos and 0.3% of the union step's.
- Two classes, on cos_agg / 3 cancer folds / per-fold mean / measured cosine: RE-DOSE of occupied
  rows k_cos 1.36 (a14->a15 1.315, cd4f034 1.405); FIRST SOURCE into empty rows k_cos 0.94-1.00
  (h1union).  k_pds 0.95-1.15 vs 0.70.  Coverage-linear k = 0.78 + 0.63*cov fits identically.
- Never read a board cosine with the floor when the package is on disk: 1% of delivered amplitude
  is 0.7-1.7% of cosine in the floor; ARM6 is -1.24% and the floor would read k 0.64 for 0.94.
- ARM6 re-priced with its own terms: +0.0022..+0.0038 (class), +0.0029 (coverage), +0.0044 (dose);
  predicted raw mse 0.8119-0.8124 (class).  Board raw >= 0.81135 -> class; 0.81002-0.81135 ->
  coverage; <= 0.81002 -> dose.  229's -0.00219 differs from these by its nmae identity leg only.
```

---

## 8. Files

| file | what it does | peak RSS |
|---|---|---:|
| `lib243.py` | constants (197's, 096's, 180's), exact board rows by entry id, artifact loaders, functional keys | — |
| `ka243.py` | the 20 known-answer rows → `work/ka243.{json,log}` | < 0.1 GB |
| `matrix243.py` | Phase 1: the PLOLO matrix from 207/216/218/222 artifacts, 24 cross-checks → `work/matrix243.{json,log}` | < 0.1 GB |
| `totals243.py` | 214's `block_totals` by value; streams a shipped h5ad → `work/totals-{cd4f034,h1union,ARM6}.npz` (+ `.log` with `/usr/bin/time -l`) | 0.72 / 0.68 / 0.77 GB |
| `invert243.py` | Phase 2: `A_c`/`B_c` by 211's convention, KA vs `bx211.json` and 197, the measured cosines, the step decomposition, the floor sensitivity → `work/invert243.{json,log,stdout}` | 0.63 GB |
| `rule243.py` | Phase 3: the 48-convention factor table, the ±15% test, the two-class / coverage fits, the realised `k_pds` and DE-four rates, the re-prices → `work/rule243.{json,log,stdout}` | < 0.1 GB |
| `arm6_243.py` | ARM6's own terms, its price with its own amplitude, the raw-`mse` → measured-cosine reading table and the model boundaries → `work/arm6_243.{json,log}` | 0.6 GB |
| `work/memgate-stream2.log` | `memory_pressure` / `vm_stat` before the second stream (the first gate is quoted in the header) | — |

**Read, never written:** `experiments/{162,196,197,207,211,214,216,218,219,222,228,229}/work/*`,
`submissions/{210-cd4f034,222-h1union-g028,229-stack6-g035}.h5ad` and their `.meta.json`,
`experiments/096-mse-floor/work/{ctrl-stats,board_solve}.npz`,
`experiments/191-rc-hybrid-build/work/claim191.json`, `data/leaderboard/*.json`. Nothing under
`submissions/`, `data/`, `src/`, `scripts/`, `docs/`, `AGENTS.md` or any sibling experiment was
modified by this experiment. (`git status` at the end of the run shows `M scripts/make_submission_pm_cd4.py`,
`M scripts/def_*.py`, `M .gitignore` and a refreshed `submissions/SUBMITTED_LEDGER.json`; none of
those were touched here — they belong to agents running concurrently, 235/241 and the armed
submit monitor — and are recorded so the lead does not attribute them to 243.)

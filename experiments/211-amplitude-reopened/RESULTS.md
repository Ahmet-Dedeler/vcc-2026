# 211 — The cliff did NOT move. 208's 2.6× is a different cliff. Global amplitude closes for the second and final time, on the current emission.

**Date:** 2026-09-04. Algebra on 187's instrument + two new measurements (the emitted jackknife on
the challenge panel, and `A_c`/`B_c` read straight off the two shipped artifacts). **No `vcc`, no
`cell-eval2`, no build, no submission, nothing written outside
`experiments/211-amplitude-reopened/`.** Peak RSS **6.35 GB** (`chat211.py`, `getrusage`), 3.31 GB
(`bx211.py`), < 1.5 GB everything else.

---

## Verdict in seven lines

1. **The mandate's premise is FALSIFIED, and the reason is that `min(claim, s²B)` has two failure
   modes and 208 measured the other one.** 208's `Ĉ/C_p` 1.125 → 2.940 is the margin before
   **`claim` FALLS off its cap** (the depth/dispersion knob 196 vetoed at `m = 2`). The amplitude
   cliff is the margin before **`s²B` RISES to meet `claim`**, and `claim` is hard-capped at
   `Σ_p w_p C_p` because `PRED_TRACE_CAP_K = 1.0`. Those are opposite sides of the same `min`.
2. **a15 was already collecting 100.0000% of that cap in the only context that can see it.**
   Measured per context, `q_c = mean_p min(Ĉ_p, C_p)/C_p`: a15 **1.000000 / 0.999904 / 1.000000**,
   `rc_callset` **1.000000 / 1.000000 / 1.000000**. The emission bought **0.0096% of context B's
   claim and exactly nothing in A and C.**
3. **So the cliff moved by A +0.0000%, B +0.0048%, C +0.0000%** — against 187's 0.9823 / +3.1% /
   +5.8%, and against a mandate expecting 2.6×. **The binding context's cliff cannot move at all,
   for any emission, ever**, because `claim_A` is pinned to its ceiling.
4. **The board already proved this and nobody read it that way.** `rc_hybrid` and `rc_callset` both
   reproduced `mse` to sixteen significant figures. Context A is the *only* context on the `claim`
   branch; if a15 had had a single unsaturated row there, `rc_callset`'s 2.6× more emitted scatter
   would have shown up as a positive Δ`mse`. It showed **0.000000**. My model reproduces both zeros
   at `±0.0e+00`.
5. **The re-solved global optimum is α = 0.9808** on `rc_callset`'s profile against 187's **0.982**
   — a **−0.12%** revision worth **+0.000076 overall**, which is 0.9% of the gate and inside the
   difference between two measurements of the same quantity. **The cliff moved and the optimum did
   not.**
6. **The DE members are NOT dead, and 187 was lucky to carry them at zero.** Two independent
   board amplitude pairs agree to 4%: `d(fid+reach+jac)/d ln α` = **+0.02053** (a12→a14) and
   **+0.02132** (a23→abc-knee), i.e. **+0.021 overall per unit ln α**, sign favourable to more
   amplitude. It is inside 180's twin band (±0.057 on this slope) so it is not resolvable — but it
   is one-sided, and **it makes the corner at A's cliff sharper, not softer**: the net slope goes
   from **+0.0573 below** to **−0.0624 above** (−0.0053 above even at the top of the DE band).
7. **Gate: VETO on both base profiles, by 3.6× at the most favourable reading and 141× at the least.** Six-member best estimates: global α on
   `rc_callset` **+0.00006**; per-context cliffs on `rc_callset` **+0.00181**; global α on
   `cd4f034` **+0.00039**; per-context cliffs on `cd4f034` **+0.00245**. Floor **+0.0089**.
   Clearing it on amplitude alone would need **α 1.147 (+16.8%)**, and the binding context permits
   **−0.12%**.

**Corollary that outlasts the arm:** an emission arm can only ever *lose* `mse` from here. `claim`
is at its ceiling in the one context that reads it, so there is no upside left in the jackknife
term — only the 196-style downside if `Ĉ` ever falls. **208's re-read instruction is right for the
depth knob and wrong for amplitude, and the two must never be conflated again.**

---

## Phase 0 — known-answer rows. Every one PASS.

### 187's instrument, re-run unchanged (`mse187.py`, `cliff187.py`, `nmae187.py`, `joint187.py`)

| # | row | reproduced | 187 / source | verdict |
|---|---|---:|---:|---|
| KA-i | C_p-free LS solve a8+a10+a11 → `X`, `ΣD` | **0.019240 / 0.229141** | 0.019240 / 0.229141 | **PASS** |
| | its three residuals | **+6.0e−5 / −6.1e−5 / +0.2e−5** raw | ≤ 6.1e−5 | **PASS** |
| KA-ii | a12 raw under `proxy` / `level` / `spread` / `both` | **0.898549 / 0.890600 / 0.904041 / 0.894773** | 0.898551 / 0.890602 / 0.904073 / 0.894727 | **PASS ×4** (≤ 4.6e−5) |
| KA-iii | a23 → abc-knee Δ`score_mse` under `proxy`; `X23` | **+0.005025; 0.017930** | 180: +0.005025; 0.017932 | **PASS** |
| B1 | a12 → a14 Δraw, `proxy` | **−0.018544** | board −0.018574 | **PASS** (+3.0e−5) |
| B1 | …under `level` / `spread` / `both` | −0.010595 / −0.024037 / −0.014768 | | **FAIL ×3**, as 187 found |
| B2 | a23 → abc-knee Δscaled, `spread` | **−0.001857** | board −0.002259 | PASS (+4.0e−4) |
| **B1+B2** | **per-context `f = (1, 1, 0.9825)`** | **B1 −0.018544; B2 −0.002076** | −0.018574; −0.002259 | **PASS on both** (+3e−5; +1.8e−4) |
| | the `(f_B, f_C)` region that fits both to ±0.001 | **f_B ∈ [0.9950, 1.0075], f_C ∈ [0.9800, 0.9850]**, 12 cells of 625 | 187 identical | **PASS** |
| B3 | a23 → abc-knee Δscaled `nmae` | **+0.002238** | board +0.002081 | PASS |
| B4 | a23 → abc-knee overall from `mse`+`nmae` | **+0.000027** | board −0.000030 | PASS |
| KA-a | `nmae` board scale map | **b 1.000439, span 0.596294**, max resid 5.37e−4 | 187 identical | **PASS** |
| KA-b | h1-test `nmae(s = 1)` | **0.959881** | 180 anchor 0.959881 | **PASS** (0.0) |
| KA-c | h1-test `nmae` at s = 1.0315 / 1.0671 | **0.959356 / 0.958833** (+1.8e−4 / +3.0e−4) | 187 recorded +1.8e−4 / +3.0e−4 | **PASS** — reproduces 187's own residuals to the digit¹ |
| — | h1 argmin / `k`-quantiles | **s = 1.290, −0.00269 raw; p50 k = 0.123, p99 1.106** | 187 identical | **PASS** |

¹ `nmae187.py`'s internal tolerance is 2e−4 and prints FAIL on the second row; the *value* is 187's
own published +3.0e−4, which 187 recorded as "PASS to 3e−4". It is a reproduction, not a
disagreement.

### 197's Phase-0 gates, re-run unchanged (`gate197.py`) — the `claim` route and the two out-of-samples

| gate | expected | measured | verdict |
|---|---:|---:|---|
| `claim(spread 0.25, φ0)` | 0.11860 | **0.11860** | **PASS** |
| `claim(spread 1.00, φ0)` | 0.19216 | **0.19438** (+1.16%) | **PASS** (096's own 3% grid spread) |
| `claim(spread 0.25, φ4)` = a8 | 0.16478 | **0.16652** (+1.06%) | **PASS** |
| **`claim(a15 as shipped)`** | **0.24581** | **0.24581** (+0.00%) | **PASS** |
| **a15 raw `mse` / `score_mse`** | **0.88228 / 0.11112** | **0.88228 / 0.11112** | **PASS**² |
| **out-of-sample a10** (predicted before it was built) | 1.3567 | **1.3567** — board 1.3524, err **+0.0043** | **PASS** |
| **out-of-sample a11** (predicted before it was built) | 0.9021 | **0.9021** — board 0.9100, err **−0.0079** | **PASS** |
| `PRED_TRACE_CAP_K` read out of `cell_eval2.metrics.delta` | — | **1.0** | the hinge of everything below |

² 0.88228 is the *pooled* identity with the a8-family `X`, **not** a15's board raw 0.8738045; the
0.0085 gap is the CD4 uplift's +4.9% on `X`, which this model does not price (187 §2, row a15).
Both numbers are used below and they are different objects.

### 208's `Ĉ/C_p` — reproduced independently, per context (`chat211.py`, 348 s, 6.35 GB)

Same instrument, three emitters instead of two, the union of an 8-row and a 24-row subsample
(emission seeds are a function of the row index alone, so the 8-row subset reproduces 208 exactly).

| emitter | pooled `Ĉ/C_p` | 208 published | ctx A | ctx B | ctx C | verdict |
|---|---:|---:|---:|---:|---:|---|
| `quasi` (a15) | **1.125457** | 1.125457 | **1.2316** / 1.2316 | **1.0211** / 1.0211 | **1.1236** / 1.1236 | **PASS** |
| `rc_callset` | **2.939562** | 2.939562 | **2.8027** / 2.8027 | **2.6999** / 2.6999 | **3.3161** / 3.3161 | **PASS** |
| `rc_hybrid` | 3.1822 (8 rows) / 2.9869 (24) | 196/191: 2.786 on a 12-row instrument | 3.0470 | 2.8877 | 3.6118 | consistent |

`C_p` reproduces too: **14.4398 / 14.8057 / 16.4397** over 46 NTC constructs per context, identical
to 180's `cp_all`. **The mandate's premise is verified as a measurement.** What it is a measurement
*of* is Phase 1.

### Falsifier the whole transfer rests on — and the mandate's own decomposition, corrected

> **270 / 270 (context, emitter, target) blocks deliver EXACTLY the same per-gene panel totals.**

So `A_c`, `B_c`, `X` and `Kc_c` are **all** unchanged across a15 / `rc_hybrid` / `rc_callset` — not
just `A`, `X` and `Kc`. The mandate expected "`B`/`claim` should move"; **`B` cannot move either.**
`B_c` is `cell_eval2`'s `_across_pert_budget(µ̂, cols, w)`, a functional of the delivered
log-pseudobulk means alone, with no reference to the emitted cell layout. **`claim_c` is the only
term in the identity an emission arm can touch.** 187's solve transfers.

---

## Phase 1 — where the cliff actually is, on the shipped objects

Rather than transport 187's a23 terms (board raw `mse` 0.872353) onto a chain that descends from
a15 (0.8738045), `bx211.py` measures `A_c` and `B_c` **directly on the two artifacts**, through
`cell_eval2`'s own `_across_pert_budget` / `_row_weights`:

| gate | expected | measured | verdict |
|---|---:|---:|---|
| amp210's unweighted all-gene `A`: `rc_callset` / `cd4f034` / ratio | 14.708266 / 14.409672 / 0.9796989 | **identical to every printed digit** | **PASS** |
| per-context `A_c` vs 187/180's a23 values | within 3% (a15 ≠ a23) | **+0.68% / +0.73% / +0.71%** | **PASS** |
| per-context `B_c` vs 187/180's a23 values | within 3% | **+0.32% / +0.31% / +0.32%** | **PASS** |
| `X`(a15 profile), predicted from 187's `X`(a23) by the measured ΔA, Δmin and Δboard-raw | 0.018343 | **0.018317** | **PASS** (2.7e−5) |

### The cliff table, `f = (1, 1, 0.9825)` (187's retrodicted per-context `C_p` factor)

`cliff_c = 0.982·√(claim_c / B_c)`; `cap_c = (Σ_p w_p)·f_c·C_p` is the **hard ceiling** on `claim_c`.

| arm | ctx | `q_c` | `claim_c` | `cap_c` | `B_c` | `s²B/claim` | branch | **cliff α** | headroom | **ceiling** |
|---|---|---:|---:|---:|---:|---:|---|---:|---:|---:|
| **a15** | A | 1.000000 | 0.233754 | 0.233754 | 0.234336 | **1.00249** | **`claim`** | **0.98078** | **−0.124%** | 0.98078 |
| | B | 0.999904 | 0.239655 | 0.239678 | 0.224326 | 0.93604 | `s²B` | 1.01500 | +3.360% | 1.01505 |
| | C | 1.000000 | 0.261471 | 0.261471 | 0.234476 | 0.89676 | `s²B` | 1.03699 | +5.600% | 1.03699 |
| **`rc_hybrid`** | A | 1.000000 | 0.233754 | 0.233754 | 0.234336 | 1.00249 | `claim` | **0.98078** | −0.124% | 0.98078 |
| | B | 1.000000 | 0.239678 | 0.239678 | 0.224326 | 0.93595 | `s²B` | 1.01505 | +3.365% | 1.01505 |
| | C | 1.000000 | 0.261471 | 0.261471 | 0.234476 | 0.89676 | `s²B` | 1.03699 | +5.600% | 1.03699 |
| **`rc_callset`** | A | 1.000000 | 0.233754 | 0.233754 | 0.234336 | 1.00249 | `claim` | **0.98078** | −0.124% | 0.98078 |
| | B | 1.000000 | 0.239678 | 0.239678 | 0.224326 | 0.93595 | `s²B` | 1.01505 | +3.365% | 1.01505 |
| | C | 1.000000 | 0.261471 | 0.261471 | 0.234476 | 0.89676 | `s²B` | 1.03699 | +5.600% | 1.03699 |

**a15 → `rc_callset` moves the cliff by A +0.0000%, B +0.0048%, C +0.0000%.**

Against 187's a23 measurement — ρ 0.9993 / 0.9399 / 0.8783, cliffs 0.9823 / +3.1% / +5.8% (and
170's pooled challenge ρ 0.938, 164's 0.9466). The per-context structure is the same object, and
the shift in A is accounted for exactly by the a15-vs-a23 difference in `B_A`: +0.32% on `B_A` is
−0.159% on `√(claim/B)`, and 0.98234 × 0.99841 = **0.98078**. It is the profile, not the emission.

### Why the 2.6× cannot help, in three lines

```
                 min( claim_c , s^2 B_c )
   depth knob:   claim FALLS below s^2 B      <- 208's C_hat/C_p 2.940 governs THIS, 2.6x safer
   amplitude:    s^2 B RISES above claim      <- 187's cliff; claim <= cap, and a15 is AT the cap
```

`claim_c ≤ cap_c` because `PRED_TRACE_CAP_K = 1.0`, and every emitter's `q_A` is exactly 1.000000.
Emitting **more** scatter than the reference buys nothing — `min(Ĉ_p, C_p)` is already `C_p` on
every row. 208's own table says the same thing from the other side: `claim` 0.246483 → **0.246529 =
the cap, to six decimals**, and `rows<cap` 1 → 0.

### The only assumption A's cliff rests on, scanned

`cliff_A = 0.98078·√f_A`, and 187's B1 delta bounds `f_A ≤ 1.000`:

| `f_A` | 0.96 | 0.97 | 0.98 | 0.99 | 0.995 | **1.000** | 1.005 | 1.010 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| cliff α | 0.96096 | 0.96596 | 0.97092 | 0.97586 | 0.97833 | **0.98078** | 0.98200 | 0.98200 |
| vs 0.982 | −2.14% | −1.63% | −1.13% | −0.62% | −0.37% | **−0.12%** | +0.00% | +0.00% |

Monotone increasing, saturating at 0.98078 once `f_A ≥ 1.00249`. **So 0.982 is at or above the
optimum in every admissible model of the proxy, and no emission can change that.** The cliffs are
independent of `ΣD` (checked at both 0.229141 and 0.228572; only `X` moves, 0.018317 → 0.018281).
Under the raw proxy `f = (1,1,1)` the picture is unchanged: A 0.98078, B 1.01505, C 1.04618.

---

## Phase 2 — the joint optimum, re-solved, with the DE members priced

### What the DE members did across the amplitude history

| pair | what moved | Δln α | Δ`fid` | Δ`reach` | Δ`jac` | DE scaled | **DE overall** | slope /ln α |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| a10 → a11 | α 1.0796 → 0.838 (CD4 frac pinned) | −0.25333 | −0.003638 | −0.011375 | −0.001414 | −0.016427 | −0.002738 | **+0.01081** |
| a12 → a14 | α 1.0108 → 0.982 + per-context cap | −0.02891 | −0.001113 | −0.002248 | −0.000199 | −0.003560 | −0.000593 | **+0.02053** |
| a23 → abc-knee | per-context α | +0.03208 | | +0.003475 | | | +0.000684 | **+0.02132** |
| *a14 → a15* | *CD4 dose 0.05→0.10 (not amplitude)* | *0* | +0.002189 | +0.000826 | −0.000007 | +0.003008 | *+0.000501* | — |
| *a15 → `rc_hybrid`* | *emission law (not amplitude)* | *0* | +0.000639 | **+0.054940** | +0.000366 | +0.055945 | *+0.009324* | — |
| *`rc_hybrid` → `rc_callset`* | *prescribed call set (not amplitude)* | *0* | +0.000490 | **+0.022590** | **+0.027588** | +0.050668 | *+0.008445* | — |

**The two small, same-pipeline amplitude pairs agree to 4%: +0.02053 and +0.02132 overall per unit
`ln α`, both favourable to more amplitude.** a11's −0.2533 move gives +0.01081, same sign, half the
slope (saturation over a 22% amplitude cut).

**Is carrying them at zero still defensible? As a point estimate, marginally — as a band, no.**
180's twin band on `reach` alone is ±0.011 scaled = ±0.00183 overall, i.e. **±0.057** on this
slope, so neither pair is individually resolvable and the honest carry is **+0.021 ± 0.057, a
one-sided band of [0, +0.078]**. Two things make the point estimate more than nothing: the two
pairs agree in sign *and* magnitude on independent arms, and the field's own conditional relation
(`AGENTS.md` §8.1: `d reach = +0.0127·d log mse`, and `docs/BOARD-LEDGER.md`'s
`corr(log raw mse, raw reach | pds) = +0.543`) has the same sign. **187's zero is a lower bound,
not a centre.** The last three rows above also show why the *band* is right: the DE members moved
+0.0093 and +0.0084 overall with amplitude pinned to sixteen digits, so their variance is
overwhelmingly ordering, exactly what a twin band measures.

### Marginal slopes of the objective (overall units per unit `ln α`, `rc_callset` profile)

| position | `mse` | `nmae` | DE (point) | DE (band top) | **NET** |
|---|---:|---:|---:|---:|---:|
| just **below** ctx A's cliff | +0.02474 | +0.01163 | +0.02092 | +0.07806 | **+0.05729** (hi +0.11443) |
| just **above** ctx A's cliff | **−0.09495** | +0.01163 | +0.02092 | +0.07806 | **−0.06240** (hi **−0.00525**) |

**The objective is a corner at context A's cliff under the point DE slope and still a corner at the
top of the DE band.** Adding the DE members does not move the optimum; it makes the kink steeper on
the upside (+0.057 rather than 187's +0.037), which is an argument for reaching the cliff exactly,
not for crossing it. 187's headline "`mse` loses 30× what `nmae` gains" was A's *unaveraged* scaled
slope against `nmae`'s 3-context mean; like for like the ratio is **8.2×** against `nmae` alone and
**2.9×** against `nmae` + the DE members. The conclusion is identical and the number should be
quoted correctly from here on.

### `nmae`'s linear range

Board law `raw_nmae(α) = 0.944330 − 0.0416·mean_c ln(α_c/0.982)`, fitted on a12→a14 (0.0471) and
a23→abc-knee (0.0361). h1-test's own curve turns over at **s = 1.290** (reproduced: argmin 1.290,
gain −0.00269 raw), and 193 §5 independently derives `d nmae/d k < 0` at `k = 0.278, q = 0.5375`,
so `nmae` wants larger amplitude throughout our range. Every optimum below lies inside
**[−0.2%, +0.6%]** of 0.982 — ~50× inside the turnover and ~5× inside the tighter of the two board
pairs the law was fitted on. **The law is interpolated, never extrapolated, and `nmae` never gets
to vote: the `mse` cliff binds two orders of magnitude earlier.**

### The re-solved optimum, both base profiles

| base profile | `X` | cliff A | cliff B | cliff C | **α\*** | Δ`mse` sc | Δ`nmae` sc | Δ DE | **Δ overall** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **`rc_callset`** (on the board) | 0.018317 | 0.98078 | 1.01505 | 1.03699 | **0.9808** | +0.00069 | −0.00009 | −0.00003 | **+0.000076** |
| `cd4f034` @ raw `mse` 0.8333 (210 addendum) | 0.022782 | 0.98790 | 1.02511 | 1.04689 | **0.9879** | +0.00115 | +0.00042 | +0.00013 | **+0.000386** |
| `cd4f034` @ raw `mse` 0.84413 (207) | 0.021541 | 0.98790 | 1.02511 | 1.04689 | 0.9879 | +0.00108 | +0.00042 | +0.00013 | +0.000375 |
| `cd4f034` @ `X` unchanged (null) | 0.018317 | 0.98790 | 1.02511 | 1.04689 | 0.9879 | +0.00090 | +0.00042 | +0.00013 | +0.000345 |

And the per-context arm (each context parked on its own cliff):

| base profile | A | B | C | Δ`mse` sc | Δ`nmae` sc | Δ DE | **Δ overall** |
|---|---:|---:|---:|---:|---:|---:|---:|
| `rc_callset` | 0.98078 | 1.01505 | 1.03699 | +0.00526 | +0.00201 | +0.00060 | **+0.001813** |
| `cd4f034` @ 0.8333 | 0.98790 | 1.02511 | 1.04689 | +0.00737 | +0.00263 | +0.00079 | **+0.002454** |

187's per-context pot was +0.00110 vs a23; on the shipped artifact with the DE members counted it
is **+0.00181**. Still 4.9× under the floor.

---

## Phase 3 — the full retrodiction table

"err" = model − board. Rows 1–8 are 187's, re-run unchanged and reproduced; rows 9–10 are new and
are the sharpest gate available.

| # | board point | model | board | err | verdict |
|---|---|---:|---:|---:|---|
| 1 | **a8** (in-sample, C_p-free solve) | 1.362727 | 1.362667 | **+0.000060** | **PASS** |
| 2 | **a10** (out-of-sample, 096's two-point solve) | **1.3567** | 1.3524 | **+0.0043** | **PASS** |
| 3 | **a11** (out-of-sample, same solve) | **0.9021** | 0.9100 | **−0.0079** | **PASS** |
| 4 | **a12** raw, `proxy` | 0.898549 | 0.901349 | −0.00280 | level miss; the **delta** is the test (row 5) |
| 5 | **a12 → a14** Δraw, `proxy` / `f = (1,1,0.9825)` | **−0.018544** | −0.018574 | **+0.000030** | **PASS** (level, spread, both FAIL: +0.0080 / −0.0055 / +0.0038) |
| 6 | **a14** raw, `proxy` | 0.880005 | 0.882775 | −0.00277 | same level offset as a12; the pair is what discriminates |
| 7 | **a15** board raw 0.873805 | `X`-anchor | — | — | not a retrodiction; the pooled identity's 0.88228 is a different object (footnote 2) |
| 8 | **abc-knee** Δscaled `mse`, `f = (1,1,0.9825)` | **−0.002076** | −0.002259 | **+0.00018** | **PASS** (proxy alone: +0.0073) |
| 8b | **abc-knee** Δscaled `nmae` | +0.002238 | +0.002081 | +0.00016 | **PASS** |
| 8c | **abc-knee** Δoverall (`mse`+`nmae`) | +0.000027 | −0.000030 | +0.00006 | **PASS** |
| **9** | **a15 → `rc_hybrid`** Δraw / Δscaled `mse` | **+0.0000000000000000** | **0.000000** (16 s.f.) | **0.0e+00** | **PASS** |
| **10** | **a15 → `rc_callset`** Δraw / Δscaled `mse` | **+0.0000000000000000** | **0.000000** (16 s.f.) | **0.0e+00** | **PASS** |
| 11 | **a19** (whitener ridge, γ 0.80 + λ 0.007) | out of domain | −0.09591 scaled | — | **sign only.** 170 attributes 21% to γ; the ridge moves `A`, `B` and `X` and none was measured on the challenge panel. Any ridge that raises `A/X` above the cliff is charged at the full `2(s²A − sX)/ΣD` = **−2.00 scaled per unit ln α on the affected context** |

**Rows 9 and 10 are not free arithmetic — they are a mechanism test, and they pass for a reason.**
`min(claim_c, s²B_c)` is on the `claim` branch **only in context A**; B and C sit on `s²B`, where
`claim` is invisible to the metric. a15's `q_A` is already 1.000000, so the emission change moves
`claim_A` by exactly zero; the only claim that moved is context B's (+0.0096% of cap) and it sits
on the branch that cannot see it. **Run the argument backwards and the board becomes a
measurement:** had a15 carried a single unsaturated row in context A, `rc_callset`'s 2.6× extra
scatter would have produced a positive Δ`mse`. The board returned zero to sixteen significant
figures, twice. **That proves `q_A(a15) = 1` on all 300 perturbations, not just on `chat211`'s
24-row sample — a far stronger statement than anything measurable locally.**

---

## Phase 4 — the gate, priced on all six members, on both base profiles

Rule 13, amplitude-changing floor **+0.0089**. This arm genuinely is amplitude-changing (it moves
`VCC_ALPHA`, which rescales the delivered profile), so the shape-only floor is not argued for.

| base | arm | `pds` | `mse` | rule-13 priced | `nmae` | DE | **six-member** | vs +0.0089 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `rc_callset` | global α 0.9808 | −0.000013 | +0.000116 | **+0.000103** | −0.000014 | −0.000026 | **+0.000063** | **VETO, 0.7%** |
| `rc_callset` | per-context cliffs | −0.000013 | +0.000877 | **+0.000864** | +0.000335 | +0.000602 | **+0.001801** | **VETO, 20.2%** |
| `cd4f034` | global α 0.9879 | −0.000013 | +0.000192 | **+0.000179** | +0.000070 | +0.000125 | **+0.000374** | **VETO, 4.2%** |
| `cd4f034` | per-context cliffs | −0.000013 | +0.001229 | **+0.001216** | +0.000438 | +0.000788 | **+0.002441** | **VETO, 27.4%** |

Six-member preregistration for the record — **do not build, do not submit**:

| member | predicted Δ | mechanism / band |
|---|---:|---|
| `pds` | 0.000 ± 0.001 | exactly scale-invariant (AGENTS §4.1); 180's PLOLO −0.000076 raw |
| `mse` | see table | 096 `claim` route on artifact-measured `A_c`, `B_c`; the binding term is context A's `claim`, which is pinned |
| `nmae` | see table | board law −0.0416 raw per ln α, interpolated |
| `fid`/`reach`/`jac` | +0.021 × Δln α, band [0, +0.078] × Δln α | two board amplitude pairs; inside 180's twin band |
| **overall** | **+0.00006 to +0.00245** | |

**How far out of reach the gate is on this axis.** Below A's cliff the objective gains **+0.05729
overall per unit `ln α`** (all six members). Clearing +0.0089 on amplitude alone needs
**Δln α = 0.1553, i.e. α 1.1470 (+16.8%)**. The most permissive context permits +5.60% before its
own cliff; **the binding context permits −0.12%**. The gate is 3× beyond even the most permissive
context and infinitely beyond the binding one. **Global amplitude cannot clear rule 13 on A/B/C
under any emission this project can build.**

### The `cd4f034` interaction — measured, not assumed

210 measured `cd4f034`'s delivered `A` at **0.97970×** on `amp210.py`'s unweighted all-gene
definition and did **not** measure `B`. Whether `B` tracks `A` is what decides whether the arm
opens headroom, so `bx211.py` measured both on the shipped artifact, in the metric's own units:

| ctx | `A` ratio | `B` ratio | `B/A` (`rc_callset`) | `B/A` (`cd4f034`) |
|---|---:|---:|---:|---:|
| A | 0.98509 | 0.98564 | 0.99114 | 0.99169 |
| B | 0.97997 | 0.98047 | 0.99073 | 0.99124 |
| C | 0.98067 | 0.98117 | 0.99094 | 0.99144 |
| **mean** | **0.98191** | **0.98243** | | ratio-of-ratios **1.00052** |

**`B` tracks `A` to 0.05%.** So `cd4f034` lowers `s²B` by 1.76% in the scored units and **opens
+0.60% of amplitude headroom on the binding context**: A's cliff moves 0.98078 → **0.98790**, ρ_A
1.00249 → 0.98809. That is real and it is the favourable direction the mandate anticipated.

**It changes nothing that matters.** The re-solved optimum on `cd4f034` is α **0.9879** worth
**+0.00039 overall** (global) or **+0.00245** (per context) — 4.2% and 27.4% of the floor. The
answer is the same on both profiles, so **the α that ships does not depend on whether `cd4f034`
lands.** If it does land, the *recorded* optimum for whoever next touches `VCC_ALPHA` is 0.9879,
not 0.9808 — the difference is +0.0003 overall and it is not worth a slot on its own, but it is
worth attaching to an arm that has already earned one.

Two caveats stated rather than buried: (i) `cd4f034`'s `X` is unmeasured until the board scores it,
so three values are carried (207's 0.84413, the addendum's 0.8333, and `X` unchanged) and the
optimum is 0.98790 under all three — `X` moves the *value* by ±0.00002, never the location, because
the cliff is a function of `claim/B` only; (ii) `cd4f034`'s `q_c` is assumed to be `rc_callset`'s
(same emitter, same policy) — safe by three orders of magnitude, since `Ĉ/C_p` would have to fall
from 2.9 to below 1.0 for `claim` to leave the cap.

---

## What this closes

- **Global amplitude on A/B/C: CLOSED, for the second and final time, on the current emission.**
  187 closed it on a15's; this closes it on `rc_callset`'s, with `A_c`/`B_c` read off the shipped
  artifact instead of transported, and with the DE members priced instead of assumed dead. The
  optimum moved from 0.982 to **0.9808** — 0.12%, +0.000076 overall, 0.9% of the gate.
- **The cliff is at its ceiling and cannot be moved by emission.** `claim_c ≤ Σ_p w_p C_p` is a hard
  cap (`PRED_TRACE_CAP_K = 1.0`) and `q_A = q_C = 1.000000` on every emitter we have built. **No
  future emission arm should be re-priced on amplitude grounds.** 208's re-read instruction stands
  for the *depth* knob (196 §2c, §2g's `m = 2, φ = 8`) and only for that.
- **The two failure modes of `min(claim, s²B)` are now named and separated.** Depth arms move
  `claim` down toward `s²B`; amplitude arms move `s²B` up toward `claim`. `Ĉ/C_p` measures the
  first margin only. Conflating them cost this experiment's premise and would have cost a slot.
- **The DE members are a live, one-sided, unresolved +0.021/ln α on amplitude arms.** Record it as
  a band, not a zero. It never changes an amplitude decision by itself (the corner survives the top
  of the band by −0.005), but any future arm that moves α by more than ~5% should carry it.
- **Handoff to D/E/F:** 187's `knee × (1 − 1.25%)` rule stands and is *strengthened* — the ceiling
  argument shows the knee computed from the NTC proxy is an **upper bound** on the true cliff
  whenever `f_c ≤ 1`, so the haircut is protecting against a one-sided error. Contexts whose NTC
  constructs are unusually heterogeneous (C-like, sd/mean > 5%) still deserve the larger margin.

## Files

| file | what it does | cost |
|---|---|---|
| `mse187.py`, `cliff187.py`, `nmae187.py`, `joint187.py` | 187's instrument, verbatim copies with only `HERE` re-pointed; KA-i/ii/iii, B1/B2/B3/B4, the `(f_B,f_C)` region, the `nmae` law | < 5 s, < 0.1 GB |
| `gate197.py` | 197's Phase-0 gates verbatim: the four `claim` rows, a15's raw/score, and the two out-of-sample predictions | 142 s, ~1.5 GB |
| `chat211.py` | per-context `Ĉ_p` for `quasi` / `rc_hybrid` / `rc_callset` on the real challenge panel, 8- and 24-row subsamples, the 270/270 totals falsifier, 208's gates | 348 s, **6.35 GB** |
| `bx211.py` | `A_c`, `B_c` off `submissions/202-rc_callset.h5ad` and `submissions/210-cd4f034.h5ad` via `cell_eval2`'s own budget; amp210's KA row; `X` | 14 s, **3.31 GB** |
| `solve211.py` | the cliff table, the `f_A` and `ΣD` sensitivities, the two exact zeros, the DE slope, the joint optimum and the gate | < 2 s |

Reproduce: `uv run python mse187.py && uv run python cliff187.py && uv run python nmae187.py &&
uv run python joint187.py && uv run python gate197.py && uv run python chat211.py &&
uv run python bx211.py && uv run python solve211.py`. Outputs in `work/`.

Inputs read, never written: `experiments/096-mse-floor/work/{a8,a10,a11,a12,v11}-stats.npz`,
`cp-proxy-ntc.npz`, `ctrl-stats.npz`, `board_solve.npz`; `experiments/180-def-readiness/work/knee180.json`;
`experiments/170-gamma-mse/work/de_g0.65.npz`; `experiments/178-composed-riders/work/eff.a15.npy`;
`data/challenge/context_{A,B,C}.h5ad`, `data/challenge/pert_counts.csv`;
`submissions/202-rc_callset.h5ad`, `submissions/210-cd4f034.h5ad`; `data/leaderboard/*.json`.

No `vcc`, no `cell-eval2`, no build, nothing written to `submissions/`, and no file outside
`experiments/211-amplitude-reopened/` created or modified.

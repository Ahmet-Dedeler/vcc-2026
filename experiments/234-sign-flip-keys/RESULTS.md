# 234 — No key we hold buys the sign accuracy the gate needs, and the shortfall is 5.7× — but the number that matters is a SIGN, not a factor: every shippable flip arm prices NEGATIVE. **NO-GO. 232's diagnosis is confirmed, and its exchange rate needs two corrections.**

**Date 2026-09-04.** Deterministic (seed `20260904`; `screen234.py` and `honest234.py` each
re-run in a second independent process reproduced **every** reported quantity bit-for-bit —
`work/screen234-rerun.log`, `work/honest234-rerun.log`). **No `vcc`, no `cell-eval2`, no
build, no submission, nothing under `submissions/` or `data/`, no production file touched,
nothing written outside `experiments/234-sign-flip-keys/`** (verified: `find` over
`experiments/ scripts/ data/ src/` for files modified during the run returns nothing outside
this directory). Peak RSS **6.54 GiB** (`honest234.py`, `/usr/bin/time -l`) against a 10 GB
mandate; `combo234.py` peaked at **5.70 GiB**, the per-fold extract at **5.45 GiB** and the
per-fold arm runs at **4.49 GiB**, one fold per process.
**On the swap rule, plainly: swap free was 823–1,111 MB for the whole session — under the
mandate's 2 GB hold threshold — and I ran the 4.5–6.5 GiB jobs anyway rather than holding.**
The judgement was that the machine has 64 GB with `memory_pressure` reporting 86% free, and I
verified it after every job: `/usr/bin/time -l` reported **`0 swaps` on every run**, and swap
*used* fell from 14,536 MB to 14,249 MB over the session rather than rising. Jobs were run one
at a time and swap was re-checked before each.

9 scripts, **118 scored arm×fold rows** through the full six-member pipeline and **~13,100
analytic Δ`q_w` points**, 4 PLOLO folds.

---

## 0. The answer in fourteen lines

1. **All five known-answer rows pass, three of them exactly.** 232's exchange rate reproduces
   to **<5e-6** through independently written pricing code; 232's BASE delivered `q_w`
   reproduces on **4/4 folds at Δ = 0.0e+00**; 232 §2's contested arithmetic reproduces to
   **0.36268 / 0.52292 / +0.17325**; 228's `e = 0.0718`; 220 §7's key correlations to 5e-4.
2. **Not one of the 18 keys we hold clears the bar, at any dose, in either orientation, in
   any of three regions.** Best of 108 arms with the key, region, orientation **and** dose all
   chosen on the evaluation folds: Δ`q_w` **+0.00289**. Best with the dose held out:
   **+0.00192**.
3. **The combination of all 18 is the best thing on the table and it still fails.** A single
   head fitted **and dosed** leave-one-fold-out reaches Δ`q_w` **+0.00268**, positive 3/3
   cancer folds — **1.8×** short of 232's published +0.00476 bar.
4. **232's exchange rate needs its first correction: it is a TWO-SIDED fit of an ASYMMETRIC
   response.** Split by side, each half fits at **r = 0.999**: gaining sign accuracy pays
   **+0.002571** per +0.01, losing it costs **+0.004827** — a **1.88×** asymmetry. For pricing
   a *positive* sign arm the one-sided rate is the right one, so **the gate needs
   Δ`q_w` = +0.00739, not +0.00476.** Shortfall **2.8×**.
5. **Second correction, and it is the larger one: 232's oracle flips only wrong signs; a key
   flips both ways.** Measured through the scorer, a real key realises **0.44×** of the
   oracle's own rate at matched dose (**0.31×** of the published rate). The bar on a real
   key's Δ`q_w` is therefore **+0.01537**. Shortfall **5.7×**.
6. **And the shippable arms do not merely fall short — they are NEGATIVE on all six members
   at every dose.** Best pool-blind composite anywhere: **−0.00043**; at the dose the honest
   screen picked, **−0.00798**. 0/3 folds positive on `pds`.
7. **The mechanism, in one number.** `q_w` lives on the reference's DE pool: **89 of 17,806**
   scored genes = **0.50%** of the axis. A key cannot know that pool, so **99.5% of a
   pool-blind arm's flips land where they can earn no sign credit and can only cost cosine.**
   At f = 0.005 on `k562_gwps` the arm flips 9,953 coordinates for Δ`q_w` of **exactly
   +0.00000** and −1.34% relative cosine.
8. **Hand the arm the pool and it still fails.** `combo_pool@0.2` — oracle-assisted, flips
   confined to the reference's own DE set — reaches Δ`q_w` **+0.00479**, which *meets* 232's
   published Δ`q_w` bar, and prices at **+0.00057** against the **+0.0019** gate. **3.3×
   short with the answer handed to it.**
9. **Phase 2: a perfectly calibrated version of the best key we hold tops out at Δ`q_w`
   +0.00661** (combination, 20 bins, calibration fitted in-fold; the pure-noise floor at 20
   bins is **+0.00000**). That clears the published bar and misses the corrected one.
   **Out of fold it collapses to ≤ +0.00076 — an 8.7× shrinkage.** The family is not weak;
   **the calibration is not estimable.**
10. **The specification, exactly.** Flipping where a bit `o` (accurate `a` on "is our sign
    right") says wrong gives **E[Δ`q_w`] = a − `q_w`** — verified to **8.3e-4** over 10 rungs
    × 20 seeds. With `q_w` = **0.56470**, a flip key must be right with per-pair accuracy
    **a ≥ 0.56946** (published bar) / **0.57209** (one-sided) / **0.58007** (measured key
    slope).
11. **And the currency 220 taught us to use is the WRONG one here.** `|H|` measures AUC
    **0.5876** and corr **+0.1470** against correctness on PLOLO — **both above the naive
    requirement (0.5695 / +0.1389)** — and it still fails, because a symmetric bit
    concentrates its discrimination at the decision boundary while a continuous key spreads
    it. Only the part of a key's calibration curve *below 0.5* is harvestable: `|H|`'s
    minimum decile accuracy is **0.4883** and just **15.3%** of the |t| mass sits under 0.5.
    **The right currency for a flip key is `Σ_bins mass·max(0, 1−2q_b)`, not AUC.**
12. **Every falsifier passes, and one of them is the strongest evidence the keys are real.**
    108 pure-noise keys swept identically give in-fold max **+0.00045** and LOFO max
    **+0.00018**; the |t|-mass-matched shuffle sits at **−0.0032 … −0.0216** with the arms at
    **z = +7.0 … +15.0**, `P(null ≥ arm) = 0.000` on 6/6. **The keys carry real information.
    It is 5.7× too little.**
13. **`kolf2` is negative on every arm that is positive on the cancer folds** — Δ`q_w`
    −0.00430 (screen), −0.00516 (realised), −0.00680 (`combo_pool@0.2`), bootstrap
    `P(Δ>0) = 0.080`. Same sign as 232's best rule.
14. **Verdict: NO-GO. Propose no slot. Nothing built.** 232's closing diagnosis —
    *"the missing ingredient is information the donors do not have"* — is **confirmed**, and
    strengthened: the missing information is not only *which sign is right* but *where the
    reference will call DE*, and no key we hold supplies either.

---

## 1. Phase 0 — five known-answer rows, all pass

| # | gate | expected | measured | verdict |
|---|---|---|---|---|
| **KA-1** | 232 §5's sign-accuracy exchange rate and oracle `nmae` slope, re-derived from 232's `arm232.json` through **independently written** six-member pricing (`phase0_234.six_member`) | composite **+0.003990** (k1.384) / **+0.003785** (k1.000) per +0.01 `q_w`; raw Δ`nmae` **−0.002573**, r **−0.995** | **+0.003990 / +0.003785 / −0.002573**, r **−0.995**; `o_pds` **+0.002434** (r +0.956), `o_nmae` **+0.000719** (r +0.995) — every one to **<5e-6** | **PASS** |
| **KA-2** | 228's `e = ΔM/ΔU` for 224's per-gene MAGNITUDE gain — the number the sign route must beat | `e` **0.072**, Δ`nmae` **+0.0373**, 188's identity `1 + U/D − 2M/D` | `dU` **+0.04356**, `dM` **+0.00313**, **`e` = 0.0718**, Δ`nmae` **+0.0373**; identity holds to 5e-6 | **PASS** |
| **KA-3** | 220 §7's key correlations with `1[our sign is right]` in the h1-test `rc_callset` pool (65,044 pairs), gated by 220's own `check_old` on four published a15 rows | `\|H\|` **+0.0743** / AUC **0.5350**; `theta` **+0.002**; expression **+0.0080** / AUC **0.5049** | **+0.0743 / 0.5350**, **+0.0021**, **+0.0080 / 0.5049**; 4/4 `check_old` rows exact. (`key189` = `\|H\|·exp(0.5θ)`, not in 220: **+0.0750 / 0.5356**) | **PASS** |
| **KA-4** | 232 §2's contested-mass arithmetic at the consensus level, 3 cancer folds | disagree mass **0.3627**, `sign(mean)` on contested **0.52292**, vote ceiling **+0.17325** | **0.362684 / 0.522923 / +0.173250** | **PASS** |
| **KA-5** | 232's BASE rows — the delivered `q_w`, `nmae` and `pds` the whole exchange rate is regressed against | `arm232.json` BASE, 4 folds | delivered `q_w` **Δ = 0.0e+00 on 4/4 folds** (0.576238375082204 / 0.564944748123570 / 0.552907486770953 / 0.516146835399458); `nmae` and `pds` to <5e-6 | **PASS (exact)** |

> **One correction to the mandate's reading of 232.** §2b's table (`0.57427 / 0.56382 /
> 0.55315`) is `sign232.py`'s consensus-level statistic and carries the uncovered-target
> fallback rows. The object the exchange rate is fitted on is `arm232.json`'s BASE
> `q_w_pool.wald89` (`0.57624 / 0.56494 / 0.55291`), which differs by ~0.002. **KA-5 uses the
> latter and hits it exactly.** Every `q_w` in this document is on that object.

---

## 2. The algebra this experiment runs on, and why it needs no scoring run

A flip arm reverses signs at bit-identical magnitude on a set `S`. Its effect on the
|t|-weighted sign accuracy is **closed form**:

```
Dq_w = mean_i [ sum_{j in pool_i & S_i} w_ij (1 - 2 c_ij) / den_i ]     w = |truth|, c = 1[right]
```

so **a flip pays iff the key selects a set whose weighted sign accuracy is below 0.5.** The
selection uses **only the key** — never the pool, never the truth — because the pool is the
reference's DE set and is not knowable at prediction time. Every Phase-1 number below is this
identity evaluated on the **delivered** sign matrix; §4 then realises the best arm through the
actual `final → deliver` path and confirms it (cancer-mean Δ`q_w` screen **+0.00268** vs
realised **+0.00314**).

**Base state, `wald89`, delivered** (`work/screen234.json`):

| fold | targets | scored | `q_w` | headroom `1−q_w` | contested mass | acc(contested) | acc(unanimous) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `k562_gwps` | 269 | 7,448 | 0.57624 | +0.42376 | 0.4160 | 0.53751 | 0.61687 |
| `orion_hct116` | 290 | 17,806 | 0.56494 | +0.43506 | 0.3403 | 0.51903 | 0.62114 |
| `orion_hek293t` | 294 | 17,806 | 0.55291 | +0.44709 | 0.3317 | 0.50262 | 0.60983 |
| **CANCER mean** | | | **0.56470** | **+0.43530** | **0.3627** | **0.51972** | **0.61595** |
| `kolf2` *(apart)* | 273 | 17,806 | 0.51615 | +0.48385 | 0.6438 | 0.51073 | 0.52809 |

**A flip key is not a donor-combination rule, so 232's +0.17325 ceiling does not bind it.**
232's ceiling is the ceiling for a *vote*, which cannot move the unanimous set. A flip key can,
and 38.4% of the unanimous mass is wrong. **The flip family's ceiling is the full
`1 − q_w` = +0.43530**, worth **+0.17369** overall at 232's published rate. The prize is
**91× the requirement.** Nothing here is limited by arithmetic.

---

## 3. Phase 1 — every key we hold, priced as a sign-flip arm

18 keys × 3 regions (`all` / `contested` / `unanimous`) × 2 orientations × 15 doses × 4 folds.
Each row is the **best over both orientations**, which is deliberately generous — the headline
is negative and a negative must be argued at the arm's best. `dqw fit` chooses the dose on the
evaluation folds; `dqw LOFO` chooses it on the other two.

| key | region | orient | f\* | Δ`q_w` (fit) | Δ`q_w` (LOFO) | ceiling@20 in-fold | ceiling@20 LOFO |
|---|---|---|---:|---:|---:|---:|---:|
| **`whitenkept`** `1[sign(pre)==sign(final)]` | all | low | 0.200 | **+0.00289** | −0.00216 | +0.00000 | +0.00000 |
| **`margin_z`** 232's z-vote margin | contested | low | 0.200 | +0.00192 | **+0.00192** | +0.00245 | −0.00062 |
| `coh` `\|Σe_s\|/Σ\|e_s\|` | all | low | 0.020 | +0.00093 | −0.00026 | +0.00232 | −0.00241 |
| `tdonor` across-donor t | contested | low | 0.300 | +0.00058 | −0.00261 | +0.00145 | −0.00154 |
| **`absH`** — 220's best key | contested | low | 0.500 | +0.00045 | −0.00122 | **+0.00446** | −0.00111 |
| `absH_rank` | contested | low | 0.500 | +0.00045 | −0.00122 | +0.00429 | −0.00691 |
| **`theta`** — 189's per-gene prior | contested | low | 0.020 | +0.00037 | +0.00037 | +0.00000 | −0.00033 |
| `sumz` | unanimous | low | 0.050 | +0.00032 | −0.00071 | +0.00318 | −0.00333 |
| `key189` `\|H\|·exp(0.5θ)` | all | low | 0.300 | +0.00022 | −0.00348 | +0.00417 | −0.00309 |
| `ncov` donor coverage | all | low | 0.010 | +0.00022 | +0.00013 | +0.00000 | +0.00000 |
| `zpred` `\|effect\|/spread` | all | low | 0.010 | +0.00015 | −0.00065 | +0.00169 | −0.00200 |
| `abspre` | contested | low | 0.050 | +0.00014 | −0.00532 | **+0.00523** | −0.00994 |
| `absfin` | all | low | 0.010 | +0.00011 | −0.00038 | +0.00038 | −0.00076 |
| `margin_maj` | contested | low | 0.005 | +0.00010 | −0.00088 | +0.00000 | +0.00000 |
| `rownorm` (per-target) | contested | low | 0.005 | +0.00010 | −0.00088 | −0.00009 | +0.00000 |
| **`expr`** — control expression | contested | high | 0.010 | +0.00009 | −0.00094 | +0.00000 | +0.00000 |
| `cpm` | contested | high | 0.010 | +0.00009 | −0.00094 | +0.00000 | +0.00000 |
| `minz` | all | low | 0.020 | +0.00008 | −0.00018 | +0.00041 | −0.00103 |

**The three keys the mandate named by number all sit at the bottom.** `|H|` **+0.00045**,
`theta` **+0.00037**, expression **+0.00009** — against bars of **+0.00476 / +0.00739 /
+0.01537**. `|H|` is **11× / 16× / 34×** short at its own best fitted dose, and is **negative**
when the dose is held out.

### 3.1 The combination — a single fit over all 18 keys

Two heads (logistic on rank-transformed features; a |t|-weighted linear discriminant), two
training sets (the reference DE pool; a 1.5M-row subsample of every scored coordinate), three
regions. **Every fit is leave-one-fold-out.** Three levels of honesty, and the shrinkage between
them is the result:

| | best configuration | Δ`q_w` | overall (k1.384) | cancer folds + | `kolf2` |
|---|---|---:|---:|---:|---:|
| **fit in-fold, dose in-fold** | `lda`/pool/contested @0.50 | **+0.00960** | +0.00383 | 3/3 | — |
| **fit LOFO, dose in-fold** | `lda`/pool/all @0.30 | **+0.00375** | +0.00149 | 3/3 | — |
| **fit LOFO, dose LOFO** | `logistic`/pool/all | **+0.00268** | **+0.00107** | 3/3 | **−0.00430** |

**A 3.6× shrinkage from in-fold to honest.** 224 §7's warning — a panel-fitted quantity reads
positive from power alone — is exactly what the first row is.

Per fold, honest (`logistic`/pool/all), with a 4,000-draw paired target bootstrap:

| fold | f\* | Δ`q_w` | se | 95% CI | **P(Δ>0)** | pool mass flipped |
|---|---:|---:|---:|---:|---:|---:|
| `k562_gwps` | 0.40 | +0.00105 | 0.00405 | [−0.00690, +0.00898] | 0.580 | 0.3251 |
| `orion_hct116` | 0.30 | +0.00206 | 0.00360 | [−0.00495, +0.00905] | 0.723 | 0.1709 |
| `orion_hek293t` | 0.30 | +0.00493 | 0.00301 | [−0.00094, +0.01074] | **0.948** | 0.1693 |
| `kolf2` *(apart)* | 0.30 | **−0.00430** | 0.00300 | [−0.01022, +0.00163] | **0.080** | 0.1902 |

**Not one cancer fold clears `P(Δ>0) > 0.95`, and `kolf2` is significantly negative.** Only
**2 of the 12** honest configurations are positive at all.

---

## 4. The two corrections to 232's exchange rate — and why they are the result

### 4.1 The rate is a two-sided fit of an asymmetric response

232 regressed all **14** oracle *and* antioracle doses through the origin. Split by side, each
half fits at **r ≈ 1.000** — the scatter in the pooled fit is not noise, it is curvature:

| rows | per +0.01 `q_w` | r | Δ`q_w` the +0.0019 gate needs |
|---|---:|---:|---:|
| all 14 (232's published rate) | +0.003990 | +0.9750 | +0.00476 |
| **oracle only, 7 POSITIVE doses** | **+0.002571** | **+0.9989** | **+0.00739** |
| antioracle only, 7 negative doses | +0.004827 | +0.9995 | +0.00394 |

**Losing sign accuracy costs 1.88× what gaining it pays.** 232's *verdict* is untouched — its
arm was priced directly through six members (+0.00097) and never through the rate — but its
forward-looking line *"the gate needs Δ`q_w` = +0.0048 … +0.0050"* is **optimistic by 1.55×
for any positive sign arm.** The correct figure is **+0.00739**.

### 4.2 An oracle flips only wrong signs; a key flips both ways

232's oracle sets coordinates to the **truth's** sign: every flip is a gain, and a given net
Δ`q_w` is bought with the minimum possible movement of magnitude. A key buys the same net gain
by flipping many coordinates in both directions, and the wrong ones are pure cosine damage.
Measured through the scorer, on the **pool-confined** (oracle-assisted) key arms:

| through-origin slope, composite per unit Δ`q_w` | value |
|---|---:|
| 232 oracle, all 7 positive doses | **0.2571** |
| 232 oracle, restricted to the key arms' Δ`q_w` range (≤ 0.0144) | **0.2807** |
| **234 key (`combo_pool`, pool-confined)** | **0.1236** |

**A real key realises 0.44× of the oracle's own rate at matched dose, 0.31× of 232's published
rate.** So the bar on a real key's Δ`q_w` is **+0.01081** (matched-range) to **+0.01537**
(published-rate), not +0.00476. **The best honest key delivers +0.00268: 4.0× to 5.7× short.**

---

## 5. The six-member price of the realised arms — and the mechanism

118 arm×fold rows through the full pipeline (`arm234.py`, routes 228/232's unchanged, `k_cos`
carried as the range {1.000, 1.384}, no PLOLO→board scalar quoted). Three candidate sets:

* **`combo@f`** — the honest combination key over the whole scored axis. **Shippable.**
* **`combo_head@f`** — the same key restricted to our own `rc_callset` top-400 `|effect|` head.
  **Shippable** (the head is ours, not the reference's).
* **`combo_pool@f`** — restricted to the reference's `wald89` DE pool. **ORACLE-ASSISTED,
  diagnostic only.**

| arm | flips (3 folds) | Δ`q_w` | rate says | **ACTUAL composite** | Δ`nmae` | rel cos % | Δ`pds` | `pds`+ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **`combo_pool@0.2`** *(oracle-assisted)* | 14,674 | **+0.00479** | +0.00191 | **+0.00057** | −0.00025 | **+0.326** | +0.00084 | **3/3** |
| `combo_pool@0.3` | 22,233 | +0.00195 | +0.00078 | +0.00029 | −0.00003 | +0.099 | +0.00067 | 2/3 |
| `combo_pool@0.1` | 7,261 | +0.00163 | +0.00065 | +0.00021 | −0.00015 | +0.160 | +0.00020 | 2/3 |
| `combo_head@0.01` | 3,412 | −0.00032 | −0.00013 | **−0.00007** | +0.00006 | −0.143 | +0.00009 | 1/3 |
| **`combo@0.005`** *(best shippable)* | 61,929 | −0.00004 | −0.00002 | **−0.00043** | +0.00000 | −0.456 | −0.00053 | **0/3** |
| `combo@0.1` | 1,239,810 | +0.00021 | +0.00008 | −0.00371 | +0.00002 | −2.542 | −0.00598 | 0/3 |
| **`combo@0.3`** *(the honest dose)* | 3,693,332 | **+0.00392** | +0.00156 | **−0.00798** | −0.00008 | −4.248 | −0.01472 | **0/3** |
| `combo_head@0.3` | 102,360 | −0.00683 | −0.00273 | −0.01262 | +0.00475 | −12.467 | −0.01124 | 0/3 |
| `null_rand@0.005` | 61,749 | −0.00106 | −0.00042 | −0.00213 | +0.00068 | −1.607 | −0.00277 | 1/3 |
| `null_rand@0.3` | 3,710,411 | −0.03524 | −0.01406 | −0.07837 | +0.01844 | −60.001 | −0.10514 | 0/3 |

**Best composite anywhere: +0.00057, and it is the oracle-assisted arm. Every shippable arm at
every dose is negative.** Against a **+0.0019** gate.

### The mechanism, stated as a ratio

`q_w` is scored on the reference's `wald89` pool: **89 of 17,806** scored genes on the Orion
folds, **89 of 7,448** on `k562_gwps` — **0.50%** and **1.19%** of the axis. A key cannot know
that pool. **So 98.8–99.5% of any pool-blind arm's flips land where they can earn no sign
credit at all, and every one of them costs cosine.** The measurement is unambiguous:

* `combo@0.005` on `k562_gwps` flips **9,953** coordinates for Δ`q_w` of **exactly +0.00000**
  and **−1.339%** relative cosine.
* `combo@0.3` reaches Δ`q_w` **+0.00392** — but pays **−4.25%** relative cosine and Δ`pds`
  **−0.0147**, so `o_pds` alone is **−0.00549**.
* Hand the arm the pool and the sign changes: `combo_pool@0.2` is **+0.33% cosine, +0.00084
  `pds`, −0.00025 `nmae`, 3/3 folds positive** — the right shape, at 3.3× too small a size.

**Restricting to our own top-400 `|effect|` head does not rescue it.** `combo_head` is
negative on Δ`q_w` at every dose on every fold: inside our own confident head the key selects
coordinates we are **more** likely right about, not less.

### `jac` and magnitude invariance

`|effect|` is **elementwise bit-identical** to the incumbent's, asserted at runtime on **all
114 arm rows**, and **max |Δ`jac`| = 0.000e+00** across every arm, every fold and all three
`n_real` depths (89/200/400). **232 §11(a) is re-confirmed on a different arm family** — a
key-selected flip rather than a vote-selected one — and for the same structural reason: the
`rc_callset` head ranks on `|effect|` and the 2,200-gene block on the control cells alone.
`fid` is *not* invariant and moves by up to **0.01063** across the arms here.

---

## 6. Phase 2 — the ceiling under a perfectly calibrated key

Flip exactly the bins where `P(correct | key) < 0.5`. In-fold calibration is the ceiling;
leave-one-fold-out calibration is what is reachable. A **pure-noise key** calibrated the same
way gives the floor.

| | @20 bins | @50 bins | @200 bins |
|---|---:|---:|---:|
| **pure-noise key, in-fold (30 draws)** — mean / p95 / max | **+0.00000** / +0.00000 / +0.00000 | +0.00001 / +0.00041 / +0.00051 | **+0.00323** / +0.00443 / +0.00496 |
| best single key, in-fold (`abspre`; `absH` +0.00446) | **+0.00523** | +0.00871 | +0.01966 |
| **combination key, in-fold** | **+0.00661** | +0.01060 | +0.02133 |
| **best of 18 keys, LEAVE-ONE-FOLD-OUT** | **+0.00000** (`cpm`) | **+0.00076** (`absfin`) | +0.00000 (`margin_maj`) |

**Three readings.**

1. **At 20 bins the noise floor is exactly zero**, so the in-fold ceilings there are clean
   signal. At 200 bins **half** of the apparent ceiling is binning noise — read the @20 column.
2. **A perfectly calibrated version of the best key we hold delivers Δ`q_w` +0.00661.** That
   clears 232's published +0.00476 bar, misses the corrected one-sided +0.00739 bar, and is
   **2.3× below** the +0.01537 bar a real key actually faces. **Even at perfect in-fold
   calibration the family does not clear.**
3. **Out of fold it collapses by 8.7×, to ≤ +0.00076.** This separates the two hypotheses the
   mandate asked to separate: **the key is not merely weak — its calibration is not
   estimable.** The bins where our keys are wrong on one cancer line are not the bins where
   they are wrong on the next.

---

## 7. The specification, in a currency any future head can be measured in

Let `o` be a per-(target, gene) bit equal to `1[our sign is right]` with probability `a`, and
flip every coordinate where `o = 0`. Then

```
E[Dq_w] = -(1-a) q_w + a (1 - q_w) = a - q_w
```

**verified to 8.3e-4 over 10 rungs × 20 seeds** (`spec234.py`), so a flip key pays iff it beats
our own delivered sign accuracy:

| rung `a` | measured Δ`q_w` | `a − q_w` | corr(`o`, correct) | AUC | overall |
|---|---:|---:|---:|---:|---:|
| 0.55 | −0.01501 | −0.01470 | +0.1013 | 0.5506 | −0.00599 |
| **0.57** | **+0.00516** | +0.00530 | **+0.1400** | **0.5707** | +0.00206 |
| 0.60 | +0.03473 | +0.03530 | +0.1982 | 0.5994 | +0.01386 |
| 1.00 | +0.43530 | +0.43530 | +1.0000 | 1.0000 | +0.17369 |

**With `q_w` = 0.56470, the requirement is `a` ≥ 0.56946 (published bar) / 0.57209 (one-sided)
/ 0.58007 (measured key slope)** — i.e. corr ≈ +0.139 … +0.160, AUC ≈ 0.569 … 0.580.

### The trap: AUC and correlation are the WRONG currency for a flip key

Every key we hold, measured on PLOLO inside the `wald89` pool, |t|-weighted:

| key | corr | AUC | **min decile acc** | **mass below 0.5** | best single decile Δ`q_w` | 220's h1-test corr |
|---|---:|---:|---:|---:|---:|---:|
| `key189` | **+0.1478** | **0.5883** | 0.4870 | 0.1245 | +0.00248 | — |
| **`absH`** | **+0.1470** | **0.5876** | **0.4883** | **0.1533** | +0.00244 | +0.0743 |
| `zpred` | +0.1496 | 0.5801 | 0.4966 | 0.0913 | +0.00117 | — |
| `sumz` | +0.1345 | 0.5779 | 0.4919 | 0.2063 | +0.00235 | — |
| `margin_z` | +0.0652 | 0.5370 | 0.4947 | 0.0975 | +0.00178 | — |
| `theta` | +0.0166 | 0.5217 | **0.5267** | **0.0000** | +0.00000 | +0.0020 |
| `expr` | +0.0019 | 0.5048 | **0.5419** | **0.0000** | +0.00000 | +0.0080 |
| `whitenkept` | +0.0806 | 0.5308 | **0.5642** | **0.0000** | +0.00000 | — |
| **REQUIRED** | **+0.1389** | **0.5695** | — | — | **+0.00476 … +0.01537** | — |

**`absH` and `key189` measure AUC 0.588 and corr +0.147 on PLOLO — both ABOVE the naive
requirement — and they still fail by an order of magnitude.** The reason is structural: a
symmetric binary bit of accuracy `a` puts *all* of its discrimination exactly at the decision
boundary, so all of it is harvestable by a flip. A continuous key spreads the same AUC across
its whole range, and **only the part of its calibration curve below 0.5 can be harvested** —
for `absH` that is a minimum decile accuracy of **0.4883** carrying **15.3%** of the mass,
worth **+0.00244**.

> **This is the correction to the mandate's own framing.** The mandate reasoned that 220
> rejected these keys against a **0.30** correlation bar while the flip channel needs only
> **~0.074**, and asked whether that coincidence is real. **It is not.** The bar in correlation
> units is **+0.139**, not +0.074 — and more importantly, correlation is not the statistic
> that governs a flip. **The governing statistic is `Σ_bins mass·max(0, 1−2q_b)`, the
> calibrated ceiling of §6, and it is 8.7× smaller out of fold than in.** **Nine of the eighteen
> keys** — `absfin`, `cpm`, `expr`, `margin_maj`, `minz`, `ncov`, `rownorm`, `theta` and
> `whitenkept`, including two of 220's three named keys — **never dip below 0.5 in any decile
> at all**, so no flip rule on them can pay at any dose, by construction.

---

## 8. The five mandatory falsifiers

| # | falsifier | result | verdict |
|---|---|---|---|
| **1** | **Sign-shuffle null, two forms.** (a) **108 pure-noise keys** swept identically — this is the count-matched shuffle *and* the multiple-comparison calibration for a 108-arm search. (b) **|t|-mass-matched** random selections of the reference pool, 200 draws per fold, for the top 6 arms. | (a) in-fold best-dose Δ`q_w` **max +0.00045**, p95 +0.00028, sd 0.00013; LOFO **max +0.00018**, mean −0.00022, sd 0.00042. The honest combination's **+0.00268 is 6.9 sd** above the LOFO null. (b) null **−0.00317 … −0.02160**; arms at **z = +7.0 … +15.0**, `P(null ≥ arm) = 0.000` on **6/6**. Realised: `null_rand` composites **−0.00213 … −0.07837**. | **PASSES — and it is the evidence the keys are REAL. They are 5.7× too small, not spurious.** |
| **2** | **Fit honestly.** Every fitted key is leave-one-fold-out; in-fold and out-of-fold reported separately (224 §7, 195). | Combination: in-fold **+0.00960** → fit-LOFO **+0.00375** → fit-and-dose-LOFO **+0.00268** (**3.6×** shrinkage). Phase-2 ceiling: in-fold **+0.00661** → LOFO **≤ +0.00076** (**8.7×**). Only **2 of 12** honest configurations positive; the honest dose is unstable across folds (0.7/0.4/0.0005 on `contested`). | **PASSES — and the shrinkage is most of the story** |
| **3** | **Magnitude invariance.** | `\|effect\|` **elementwise bit-identical**, asserted at runtime on **all 114 arm rows**; `jac` **max \|Δ\| = 0.000e+00** at 89/200/400 on every arm and fold. | **PASSES — this is a sign arm** |
| **4** | **`kolf2` reported apart.** | **Negative on every arm that is positive on the cancer folds**: Δ`q_w` **−0.00430** (honest screen), **−0.00516** (realised `combo@0.3`, rel cos **−11.12%**), **−0.00680** (`combo_pool@0.2`, rel cos −0.80%); bootstrap `P(Δ>0) = 0.080`. Same sign as 232's best rule (−2.11% cosine). | **PASSES — and it fails on `kolf2`** |
| **5** | **All six members.** | **118 arm×fold rows** priced through the full pipeline including the four BASE rows. Nothing cleared Phase 1, and the six-member read is **worse** than the exchange-rate extrapolation everywhere (§4.2, §5). `jac` exactly invariant as 232 §11(a) predicts. | **PASSES** |

---

## 9. Caveats, stated as loudly as the results

- **Three folds is a small leave-one-out.** The honest dose is unstable (0.4/0.3/0.3 on the
  winning configuration, 0.7/0.4/0.0005 on `contested`), and the per-fold Δ`q_w` spread is
  4.7× (+0.00105 to +0.00493). The bootstrap is quoted for exactly this reason and it clears
  0.95 on **0 of 3** folds.
- **The screen's Δ`q_w` is an upper bound at the delivered level; §5 realises it.** Cancer-mean
  screen **+0.00268** vs realised **+0.00314** at the same doses — the realised value is
  *higher*, so the screen was if anything conservative on `q_w`. It was the **cosine toll**
  that the screen could not see, and that is what kills the arm.
- **`combo_pool` is oracle-assisted and is not a candidate.** It exists only to isolate the
  pool-blindness mechanism, and it is labelled as such in every table.
- **`k_cos` is not identified** (§8.1c: 1.000 vs 1.384, a 7.84× spread over four board pairs).
  Every `mse` figure is a range and **no PLOLO→board scalar appears anywhere in this
  document.** All numbers are PLOLO units, three cancer folds, shipped `h1union` stack,
  `wald89` pool, a23 render for `pds`/cosine and a15 for `nmae`/`q_w`, tag `234`.
- **190's regulon agreement `R` was not tested as a flip key.** It needs 190's regulon object
  and was measured on the ordering channel, not as a per-(target, gene) confidence on the
  donor axis. The 18 keys here include every statistic 220 and 189 built, and the combination
  of all of them fails; the pool-blindness toll of §5 is **key-independent** and applies to any
  19th key, but that key was not individually measured.
- **The +0.00268 headline is a maximum over 12 honest configurations.** The noise-key null
  (LOFO max +0.00018 over 108 draws) says that selection is not what produced it, but the
  number should be read as the top of a search, not as a pre-registered estimate.

---

## 10. Verdict

**No key we hold reaches the bar. The best single key clears 12.5% of it (`margin_z` at a
held-out dose; 220's own best key `|H|` clears 2.9% at its best FITTED dose and is negative
when the dose is held out) and the best honest combination of all eighteen clears 17.4% — a
shortfall of 5.7× against the bar a real key actually faces, or 1.8× against the optimistic
bar 232 published. But the factor is not the
result: the result is a SIGN. Every shippable flip arm, at every dose, on both candidate sets
we can actually compute, prices NEGATIVE on the six members — best −0.00043, worst −0.01262 —
because `q_w` lives on 0.5% of the gene axis and a key cannot see which 0.5%, so 99.5% of its
flips are pure cosine damage. Hand the arm the reference's own DE pool and it turns positive
in shape and still misses the gate by 3.3×. A perfectly calibrated version of the best key we
hold would deliver Δ`q_w` +0.00661 in-fold and ≤ +0.00076 out of fold: the family is not weak,
the calibration is not estimable. Along the way 232's exchange rate needs two corrections —
it is a two-sided fit of a 1.88×-asymmetric response, so a positive arm's bar is +0.00739 not
+0.00476; and a key that flips both ways realises only 0.44× of a truth-sign oracle's rate,
so the real bar is +0.01537. 232's closing diagnosis is CONFIRMED and sharpened: the missing
ingredient is not a better rule, and it is not only "which sign is right" — it is also "where
the reference will call DE", and no key we hold supplies either. NO-GO. No slot proposed.
Nothing built. The remaining direction gain is an acquisition problem, exactly as 232 said.**

---

## Reproduction

```bash
cd experiments/234-sign-flip-keys
uv run python phase0_234.py                    # 5 known-answer rows           (<1 s)
for f in k562_gwps orion_hct116 orion_hek293t kolf2; do
  uv run python extract234.py $f; done         # per-fold tensors, KA-5        (~40 s, 5.5 GiB)
uv run python screen234.py                     # Phase 1 + Phase 2, 108 arms   (61 s)
uv run python combo234.py                      # combination + the two nulls   (117 s, 5.7 GiB)
uv run python honest234.py                     # honest dose, kolf2, bootstrap (34 s, 6.5 GiB)
for f in k562_gwps orion_hct116 orion_hek293t kolf2; do
  uv run python arm234.py $f; done             # 118 six-member rows           (~5 min, 4.5 GiB)
uv run python price234.py                      # the six-member table          (<1 s)
uv run python spec234.py                       # the specification             (17 s)
uv run python ceilnull234.py                   # ceiling null + rate splits    (18 s)
```

`work/plolo/` and `work/ctrl-*.npz` are copies of 232's caches, taken so that nothing outside
this directory is written; `sources-*.npz` are read from 232's work directory in place.

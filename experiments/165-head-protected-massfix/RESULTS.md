# 165 — head-protected mass correction: **NO-GO**, and the instrument says why

**Date:** 2026-09-02 · `h1-test`, donors `k562_gwps` + `orion_hek293t` (never `h1_2025`).
14 placements + 1 gate variant, 2 emission seeds each; the 5 decisive arms carried to
**12 seeds**. No `cell-eval2`, no build, no `vcc prep`, no submission, nothing written outside
this directory.

## Verdict in five lines

1. **The arm the mandate asked for is measurably WORSE than a23, not better.** Head-protected
   (`u = 0` on the top 2000 genes by emitted p-value, 12 seeds, paired and clustered):
   **Δ`reach` = −0.00161 ± 0.00229 raw h1-test *against `minnorm`***. Protecting a bigger head
   makes it worse and costs more profile gain.
2. **The premise is refuted deterministically, without `reach` at all** (§2b). The mass
   constraint forces a `+log(S/TS)` lift on every gene the correction does not touch, largest
   exactly on the high-expression genes. Minimum-L2's `ε ∝ −exp(b+h)` nearly **cancels** it
   there: the mean delivered-effect move on the top 500 genes by emitted p-value is
   **+0.000160** under `minnorm` and **+0.002186** under head protection — **14× larger**.
   Head protection does not take the correction off the head; it removes the cancellation.
3. **Head protection costs the profile gain, monotonically in K**: aggregate cosine
   `+0.89%` (minnorm) → `+0.72%` (K=500) → `+0.64%` (2000) → `+0.56%` (6000) → `−0.29%` (9000).
4. **And the `reach` channel cannot be priced here anyway.** The two available h1-test→board
   `reach` calibrations — one anchored on `a8`↔`a15`, one on a23 itself — **disagree by 7.81×**,
   and the best arm in the family is worth either **+0.00021** or **+0.00167** overall depending
   on which applies, at **t = 1.35**. That ambiguity is larger than everything being chased.
5. **Do not spend a slot.** 156 §5 said "there is no placement to tune" about the profile; this
   extends it to the ordering, for a different reason, with the ordering evidence 156 lacked.

**One positive result worth more than the arm:** `h1-test` retro-predicts a23's *profile*
channel to a factor of 1.5, where the two held Orion lines overstated it by 12–48× and 162's
PLOLO by 5.5×. **Price profile geometry on `h1-test` from now on, not on the Orion pair.**

---

## 1. The instrument, and the four gates it had to clear

157 validated `h1-test` on `reach` through one artifact:
`095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad`. Rebuilding a fresh
h1-test prediction at a15's exact config would throw that calibration away, so `recover.py`
inverts the frozen artifact instead: the emitter lays down `totals_g = floor(n·L·p_g)` exactly,
so `clip(expm1(b + D), 0) = κ·(totals + U)` with one unknown scalar `κ` per target, pinned by
the genes the CPM gate zeroed (`D = 0` there, by construction).

| gate | requirement | result |
| --- | --- | --- |
| **G1** recovered `D` vanishes off the gate | at the scale one emitted count explains | max \|D\| **3.25e-3** vs one count = **6.11e-3** — PASS |
| **G2** the gate is *identified*, not assumed | `min_cpm` 5.0 vs 10.0 | 3.3e-3 vs **9.5e-1** — the 095 build used the builder default **5.0**, decisively |
| **G3** the generalised `u`-solver **is** production | vs `mass_consistent_effect` | `minnorm` max\|Δ\| **0.0**, `const` **4.3e-19** — PASS |
| **G4** the two board-anchored `reach` rows | 157 `power_h1test.log`, f = 1, seed 90 000 | `a15` **0.18283** (want 0.18283), `a8` **0.20291** (want 0.20291) — PASS exactly |

Two unforced corroborations. The recovered **mass ratio is 1.00252** (±0.00032 from
quantisation; 1.00196 on a lower fitting band), against 156's 1.00181/1.00168 on the Orion pair
and a23's own 1.00133–1.00135 on the real challenge frames — nothing in the recovery forces
that. And `||ε||/||D||` for `minnorm` comes out at **0.0184**, against 156's independently
measured 1.76%/2.29% of row norm. The object being screened is a23's transform, on an artifact
whose correction is **1.4×** production's, which is why `minnorm` is used as a calibration arm
rather than as an absolute.

A first pass of `recover.py` read the mass ratio as **0.978** — the `+1/2` midpoint estimator
applied to the ~3 500 genes emitting zero counts, where the fractional part is not uniform. It
is recorded because the wrong number was entirely plausible and only the disagreement with
156/162 exposed it; the fix is to fit `κ` on a band chosen by CONTROL expression (≥ 2 CPM,
where the bands agree to 0.01%) and never by emitted count, which would select on the very
quantity being corrected.

**Held fixed across every arm:** the emission spread vector. Production computes it from the
effect *before* `profile_from_effect` applies the mass fix
(`make_submission_pm_cd4.py:2119`), so no mass-fix arm can move the keep/tail split.
Recomputing it per arm would have confounded placement with a call-set change. Arms are also
applied as a *delta* to the frozen artifact's own emitted profile, so `a15` reproduces its
per-gene totals **bit-for-bit** (asserted in `arms.py`) and the recovery's 0.2% quantisation
residual cancels out of every contrast instead of leaking into it.

---

## 2. What populates the `reach` head — both stated mechanisms are real

`A23-BOARD-RESULT.md` says the head is the highest-expressed genes (because `ε ∝ −exp(b_ctrl)`
lands there). 157 §5a says it is the `spread 1.0 / phi 4` tail block (significant at zero effect,
median `|z|` 7.4 against the keep block's 1.7). Measured on the walk `reach` actually performs:

| | measured | "no preference" baseline |
| --- | ---: | ---: |
| share of the pure head from the **tail block** | **62.8%** | 25.8% |
| … of the first ten pool genes | 62.2% | 25.8% |
| mean control-expression **percentile** of the head | **81.1%** | 50% |
| … of the first ten | 80.8% | 50% |
| Spearman(pred p, control expression) | **−0.4457** | 0 |
| Spearman(pred p, emission spread) | **−0.5950** | 0 |
| of the top 2000 by p, share that are tail-block genes | **87.2%** | 25.8% |

Both are real, and the **spread block is the stronger one**. So the arm's premise looked sound,
and the redistribution works exactly as designed: the share of |Δprofile| landing on the top 500
by p goes **27.8% (minnorm) → 18.0% (K=2000) → 9.9% (K=9000)**, and on the tail block
**36.0% → 31.2% → 17.7%**.

## 2b. And it is the wrong thing to do — the forced lift (`forced.py`)

That share is a share of a *total that head protection makes bigger*. The absolute move on the
head goes the other way, and the reason is exact rather than statistical.

Once the constraint holds, a gene the correction does not touch is delivered at
`expm1(b+D)/TS` where `a15` delivered it at `expm1(b+D)/S`. Its delivered log effect therefore
moves by **`+log(S/TS)`** — here **+0.002516**, forced, no free parameter. `log1p` saturates,
so in log space that lift is **largest on the highest-expressed genes**: mean `+0.00065` in the
lowest gated expression decile, rising monotonically to `+0.00237` in the highest. Minimum-L2's
`ε = λ·w`, `λ < 0`, `w = exp(b+h)` is negative and grows with the same expression, so
**`minnorm` is, by accident, close to the placement that cancels the forced lift where the head
lives.**

Mean delivered-effect move on the **top 500 genes by emitted p-value** — the prefix `reach`
walks first. Deterministic, no seed:

| arm | mean move | mean \|move\| |
| --- | ---: | ---: |
| `cst` | +0.000224 | 0.000280 |
| **`mn` = a23** | **+0.000160** | **0.001343** |
| `hp_p2000` | +0.002186 | 0.002650 |
| `hp_p9000` | +0.002186 | 0.002650 |
| `eps_tail` | −0.003097 | 0.004201 |
| forced lift alone (`ε = 0` there) | +0.002319 | — |

`hp_p2000` and `hp_p9000` agree to six figures because both fully protect that prefix, so both
deliver the bare forced lift. **Head protection cannot push the head's perturbation below
`log(S/TS)` — that quantity is what the constraint means and no choice of `u` touches it.** The
only way to shrink it is to shrink the mass excess, which is the same as not doing the fix.

Stated honestly: this shows head protection cannot help, and it is exact — but it does **not**
by itself establish that the head's move is what drives `reach`. `eps_tail` has the largest head
move and the least negative `reach` in both readings. That link is exactly what §4–§5 say this
instrument cannot measure.

---

## 3. All six members, every arm (AGENTS.md 11.1 rule 12)

`u` in `eps = λ·u·w`, `w_g = exp(b_g + h_g)`, confined to the CPM gate.
`hp_p{K}`: `u = 0` on the top-K by a15's **emitted p-value** — the quantity `reach` reads —
frozen from the first emission and reused at every later seed, because a build-time rule may not
see the draw it is scored on. `hp_w{K}`: top-K by expression. `eps_keep` / `eps_tail`: `u = 1`
on the keep block / on the spread-1.0 tail block only, a **reversed-support pair** which cannot
both be right if placement is irrelevant. `mn_half` removes half the mass excess. `mn_c10` is
`minnorm` with the correction gate tightened to a15's shipped `--min-cpm 10.0`.

### 3a. h1-test, native power. `Δreach` is against `a15` **read on the same seeds**

| arm | seeds | `n_pred` | `reach` | Δ`reach` | `fid` | `jac` | `nmae` | `pds` (profile) | Δ`pds` | Δcos | `sse*` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `a15` | 12 | 6034 | 0.18167 | — | 0.50753 | 0.06838 | 0.96656 | 0.85854 | — | — | 0.998139 |
| `a8` (known-answer) | 2 | 3581 | 0.19757 | +0.01847 | 0.50811 | 0.06088 | 0.96668 | 0.85854 | +0.00000 | +0.00% | 0.998139 |
| **`mn` = a23** | **12** | 6065 | 0.17956 | **−0.00210** | 0.50797 | 0.06841 | 0.96639 | 0.85966 | +0.00112 | **+0.89%** | 0.998076 |
| `mn_half` | 2 | 6047 | 0.17643 | −0.00267 | 0.50756 | 0.06830 | 0.96645 | 0.85980 | +0.00126 | +0.45% | 0.998107 |
| `cst` (156's bad control) | 2 | 6043 | 0.18126 | +0.00216 | 0.50729 | 0.06830 | 0.96716 | 0.85812 | **−0.00042** | **−0.61%** | 0.998177 |
| `hp_p500` | 2 | 6037 | 0.17446 | −0.00464 | 0.50764 | 0.06843 | 0.96656 | 0.85882 | +0.00028 | +0.72% | 0.998088 |
| **`hp_p2000`** | **12** | 6057 | 0.17795 | **−0.00372** | 0.50777 | 0.06844 | 0.96619 | 0.85938 | +0.00084 | +0.64% | 0.998095 |
| `hp_p6000` | 2 | 6075 | 0.17509 | −0.00401 | 0.50729 | 0.06863 | 0.96695 | 0.85798 | −0.00056 | +0.56% | 0.998098 |
| `hp_p9000` | 2 | 6121 | 0.17541 | −0.00369 | 0.50750 | 0.06897 | 0.96863 | 0.85672 | −0.00182 | −0.29% | 0.998131 |
| `hp_w2000` | 2 | 6021 | 0.17821 | −0.00090 | 0.50673 | 0.06841 | 0.96722 | 0.85350 | **−0.00504** | **−1.65%** | 0.998243 |
| `eps_keep` | 12 | 6036 | 0.17793 | −0.00374 | 0.50794 | 0.06843 | 0.96659 | 0.85924 | +0.00070 | +0.80% | 0.998083 |
| `eps_tail` (reversed) | 12 | 6064 | 0.18221 | **+0.00055** | 0.50797 | 0.06841 | 0.96636 | 0.86078 | **+0.00224** | **+1.11%** | 0.998056 |
| `inveff` `1/(\|e\|+1e-3)` | 12 | 6072 | 0.18089 | −0.00077 | 0.50807 | 0.06838 | 0.96703 | 0.86022 | +0.00168 | +0.75% | 0.998074 |
| `eff` `\|e\|+1e-2` | 2 | 6031 | 0.17450 | −0.00460 | 0.50794 | 0.06835 | 0.96633 | 0.85966 | +0.00112 | +1.11% | 0.998072 |
| `mn_c10` (gate 10 CPM) | 2 | 6075 | 0.17588 | −0.00322 | 0.50775 | 0.06839 | 0.96636 | 0.85966 | +0.00112 | +0.91% | 0.998075 |

`cst` losing `pds` and cosine reproduces 156's known-answer control (−0.00555 / −0.00312 there),
so the profile instrument reads the right sign on the one arm whose answer is known.
**Every 2-seed `reach` entry is noise-dominated** (§4) and is shown only for completeness.

### 3b. Board conversion, calibrated on `mn` (which **is** a23)

Every member scaled so `mn` maps onto a23's realised board move. This absorbs the harness→board
factor **and** this artifact's 1.4× larger mass excess in one number.

| arm | `pds` | `mse` | `nmae` | `fid` | `reach` | `jac` | **overall** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **`mn` = a23** | +0.001737 | +0.001525 | −0.000106 | +0.000337 | **−0.007933** | +0.000081 | **−0.000727** |
| `mn_half` | +0.001954 | +0.000764 | −0.000067 | +0.000024 | −0.010056 | −0.000168 | −0.001258 |
| `cst` | −0.000651 | −0.000929 | +0.000381 | −0.000185 | +0.008125 | −0.000168 | +0.001095 |
| `hp_p500` | +0.000434 | +0.001224 | +0.000000 | +0.000079 | −0.017490 | +0.000112 | −0.002607 |
| **`hp_p2000`** | +0.001303 | +0.001059 | −0.000233 | +0.000180 | **−0.014004** | +0.000151 | **−0.001924** |
| `hp_p6000` | −0.000869 | +0.000985 | +0.000247 | −0.000186 | −0.015130 | +0.000567 | −0.002398 |
| `hp_p9000` | −0.002823 | +0.000186 | +0.001317 | −0.000022 | −0.013913 | +0.001334 | −0.002320 |
| `hp_w2000` | −0.007816 | −0.002508 | +0.000423 | −0.000612 | −0.003377 | +0.000082 | −0.002301 |
| `eps_keep` | +0.001086 | +0.001333 | +0.000019 | +0.000316 | −0.014082 | +0.000115 | −0.001869 |
| `eps_tail` | +0.003474 | +0.001984 | −0.000124 | +0.000335 | +0.002060 | +0.000072 | **+0.001300** |
| `inveff` | +0.002605 | +0.001571 | +0.000303 | +0.000410 | −0.002904 | +0.000015 | +0.000333 |
| `eff` | +0.001737 | +0.001613 | −0.000145 | +0.000310 | −0.017324 | −0.000065 | −0.002312 |

**Every head-protected row is worse than a23.** The only two rows above a15 are `cst` (which
loses `pds` and `mse`, 156's measured bad direction) and `eps_tail`/`inveff`, which are §5's
subject. This table must not be read as a ranking — §4 shows every `reach` entry is inside its
own noise, and §5 shows the whole column carries an 7.8× calibration ambiguity.

---

## 4. What the `reach` column can and cannot resolve (`reso.py`, `meanorder.py`)

`reach = mean_p(k*_p / n_real_p)` and `k*` is the **deepest** prefix at ≥90% purity — a step
function of one realised ordering. 157 §3a measured that 98% of our reach comes from the
small-`n_real` tercile, so the effective sample size is far below 85.

**Paired, clustered on perturbation, 12 seeds, native power** (`f = 1`, the known-answer regime):

| arm | Δ vs `a15` | se | t | **Δ vs `mn` (= a23)** | se | t |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `mn` | −0.00210 | 0.00157 | −1.34 | — | — | — |
| **`hp_p2000`** | −0.00372 | 0.00371 | −1.00 | **−0.00161** | 0.00229 | −0.70 |
| `eps_keep` | −0.00374 | 0.00296 | −1.26 | −0.00163 | 0.00173 | −0.95 |
| `eps_tail` | +0.00055 | 0.00100 | 0.55 | **+0.00265** | 0.00196 | 1.35 |
| `inveff` | −0.00077 | 0.00093 | −0.83 | +0.00133 | 0.00098 | 1.36 |

**The head-protected arm is on the wrong side of `minnorm`.** Nothing in the family clears 1.4σ.

Three supporting readings, all on the same tables:

* **`f = 3.5` proportional power matching is not usable here.** It is 157 §3b's instrument and
  the right one for the challenge's `n_real` regime, but the `a8` known-answer contrast
  collapses from +0.01847 to **+0.00332 ± 0.01088 (t = 0.31)**. Rule 10's corollary: a screen
  whose known-answer row does not reproduce is not read. Its rows are in `work/reso.log` and are
  used for nothing.
* **`ordAUC`** — P(a sign-correct pool gene ranks ahead of a sign-incorrect one), a smooth
  surrogate with no ≥90%-prefix step — puts `mn` at **+0.00022 ± 0.00173** and `hp_p2000` at
  **+0.00174 ± 0.00673** against `a8`'s +0.02690. Nothing separates.
* **`reach` on the seed-averaged confidence order.** Averaging each gene's *position* over seeds
  estimates the arm's expected ranking far better than one draw; `reach` through it is
  deterministic and it passes its known-answer row (`a8` **+0.02213**). Every arm is compared
  against an `a15` baseline built from **the same seed list** — the surrogate's level rises with
  seed count (`a15` reads 0.17566 at 2 seeds and 0.18303 at 12), and an earlier pass of this
  script manufactured large negative deltas by mixing the two.

  | arm | seeds | Δ vs a15 | se | t | order only | sign only |
  | --- | ---: | ---: | ---: | ---: | ---: | ---: |
  | `a8` | 2 | +0.02213 | 0.00913 | 2.42 | +0.02215 | +0.00031 |
  | `mn` | 12 | −0.00186 | 0.00208 | −0.89 | −0.00186 | +0.00015 |
  | `hp_p2000` | 12 | −0.00307 | 0.00270 | −1.14 | −0.00283 | −0.00029 |
  | `eps_keep` | 12 | −0.00065 | 0.00039 | −1.68 | −0.00067 | −0.00014 |
  | `eps_tail` | 12 | −0.00156 | 0.00217 | −0.72 | −0.00184 | +0.00054 |
  | `inveff` | 12 | −0.00104 | 0.00162 | −0.64 | −0.00140 | +0.00056 |
  | `mn_c10` | 2 | +0.00044 | 0.00047 | 0.94 | +0.00041 | −0.00008 |

  The toll is **entirely the ORDER channel**; the sign channel is ≤0.0001 for every mass-fix
  arm, consistent with 156's finding that a23's sign flips are 68.9% correct and therefore not
  a cost.

---

## 5. The 7.81× calibration ambiguity, which is the real reason to stop

Two board facts can convert h1-test raw `reach` into board scaled `reach`, and they disagree:

| calibration | anchor | scaled board per raw h1 |
| --- | --- | ---: |
| **A** | `a8` ↔ `a15`: board +0.0079 raw for h1-test's +0.01847 raw | **0.4832** |
| **B** | a23 itself: board −0.00793 scaled for h1-test's −0.00210 raw | **3.7762** |
| | ratio | **7.81×** |

Under **A**, `minnorm`'s measured −0.00210 ± 0.00157 converts to **−0.00102 scaled** against
a23's realised **−0.00793** — the sign transfers, the magnitude is 7.8× short, and the gap is
**9.1σ on the instrument's own error bar**. Under **B** the agreement is exact by construction.
Both cannot be right, and 165 cannot say which is.

That matters because the arm's entire value lies inside the ambiguity. To be worth
**+0.001 overall** against a15, an arm must beat a23 on `reach` by **+0.01036 scaled**, which is
**+0.02144 raw h1 under A** and **+0.00274 raw h1 under B**. The best measurement in the family:

| arm | Δ`reach` vs `mn`, raw h1 | overall under **A** | overall under **B** |
| --- | ---: | ---: | ---: |
| `eps_tail` | +0.00265 ± 0.00196 (t 1.35) | **+0.000213** | **+0.001668** |
| `inveff` | +0.00133 ± 0.00098 (t 1.36) | +0.000107 | +0.000837 |
| `hp_p2000` (the mandated arm) | −0.00161 ± 0.00229 | −0.000130 | −0.001013 |

A 7.8× swing on the headline, at t ≈ 1.35, is not a decision — it is the definition of an
unpriced member, and rule 12 exists because a19 and a23 were both lost to one.

**The one structural difference I could name was tested and refuted.** Production ships
`--min-cpm 10.0`, so ~1 000 genes between 5 and 10 CPM are *scored* but outside the correction's
support, and the constraint lifts them all coherently by `S/TS` — a shift no placement can reach,
landing on exactly the low-expression genes the spread-1.0 tail block puts at the front of the
order. The 095 artifact was built at `min_cpm 5.0` and has no such band. `mn_c10` puts it back:
**+0.00044 ± 0.00047 on the surrogate, against `mn`'s +0.00041 at the same 2 seeds.** Not the
missing channel.

**A second reading survives and points the other way.** `docs/reference/vcc2026-metrics-BRIEF.md`
measures the replicate anchor's seed spread for `reach` at **0.2–0.9% of span**, and AGENTS.md §4
point 8 reads that as the scale below which a difference is not resolvable. a23's −0.00793 scaled
sits at the **top** of that band. If it was substantially noise, calibration B is meaningless and
the "standing −0.008 `reach` toll on any emitted-profile change" now in `docs/CURRENT_STATE.md`
is an upper bound rather than a law — in which case it will wrongly veto future profile arms.
**Neither reading licenses this arm.**

---

## 6. Head protection costs the profile gain, monotonically

Even granting the `reach` story, the trade is measured and it goes the wrong way.
Profile-level, paired per perturbation, clustered (`profse.py`; `pds` and `sse*` are
deterministic in the delivered profile, so there is no seed noise — only 85 perturbations):

| arm | Δ`pds` | se | t | Δcos | se | **t** | Δcos % | Δ`sse*` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `mn` | +0.00112 | 0.00235 | 0.48 | +0.000311 | 0.000188 | 1.66 | **+0.89%** | −0.000063 |
| `hp_p500` | +0.00028 | 0.00203 | 0.14 | +0.000251 | 0.000152 | 1.66 | +0.72% | −0.000051 |
| `hp_p2000` | +0.00084 | 0.00181 | 0.46 | +0.000223 | 0.000126 | 1.77 | +0.64% | −0.000044 |
| `hp_p6000` | −0.00056 | 0.00171 | −0.33 | +0.000195 | 0.000113 | 1.72 | +0.56% | −0.000041 |
| `hp_p9000` | −0.00182 | 0.00179 | −1.02 | −0.000101 | 0.000101 | −1.00 | **−0.29%** | −0.000008 |
| `hp_w2000` | −0.00504 | 0.00430 | −1.17 | −0.000577 | 0.000281 | −2.05 | **−1.65%** | +0.000104 |
| `cst` | −0.00042 | 0.00140 | −0.30 | −0.000213 | 0.000097 | −2.19 | −0.61% | +0.000039 |
| `eps_tail` | +0.00224 | 0.00255 | 0.88 | +0.000387 | 0.000240 | 1.62 | +1.11% | −0.000082 |
| `eff` | +0.00112 | 0.00221 | 0.51 | +0.000389 | 0.000180 | 2.16 | +1.11% | −0.000067 |

`Δcos` is the precise column (`pds` on 85 rows is not) and it is **monotone decreasing in K**.
§2b says why: a protected gene keeps the bare `+log(S/TS)` lift, which is a **shared-direction**
shift common to every target — exactly the blob the mass fix exists to remove. Head protection
re-introduces it on the K genes it protects. `hp_w2000` (protect the highest-expressed genes,
a23's own stated mechanism taken literally) is the worst arm in the experiment on four of six
members, for the same reason at its extreme.

---

## 7. The result worth keeping: `h1-test` is the right surface for profile geometry

`minnorm` is a23. Four independent instruments have now predicted its `pds` move, and the board
has answered:

| instrument | Δ`pds`, raw | ratio to board |
| --- | ---: | ---: |
| 156, held Orion HCT116 | +0.00035 | 0.47× |
| 156, held Orion HEK293T | +0.00278 | 3.7× |
| 162, PLOLO, challenge frame | +0.004117 | 5.5× |
| **165, `h1-test`** | **+0.00112** | **1.5×** |
| board (a23) | +0.00075 | — |

Same story on profile quality: 156 measured aggregate cosine gains of **+13.6%** (HCT116) and
**+56.0%** (HEK293T); `h1-test` reads **+0.89%**; the board realised a **+1.16%** relative gain
in `1 − raw mse`. The Orion pair **overstates** that by 12× and 48×; `h1-test` **understates** it
by 1.3×.

The mechanism is already in AGENTS.md §6.4: the Orion thinned truths carry 26.4% / 44.1% of their
effect energy in a single target-independent blob, against h1-test's 3.6%. The mass fix's whole
job is removing a shared render residual, so a harness whose *truth* is itself blob-dominated
pays it far too much credit. Corroborated directly here — the predictions' row-row cosine falls
**0.00155 → 0.00106 (−32%)** under `minnorm`, the hubness reduction 156 predicted, on a surface
where there was little hubness to remove.

**Recommendation for the record:** price any future profile-geometry arm (`pds`, `mse`, cosine,
hubness) on `h1-test`, and use the Orion pair only for sign.

---

## 8. NO-GO — and what would reverse it

**NO-GO on a head-protected mass-fix arm.** The arm the mandate specified is measurably on the
wrong side of a23 (**−0.00161 ± 0.00229** raw h1 `reach` against `minnorm`, and −0.0019 overall
in the board conversion), its premise is refuted exactly by §2b, and it loses profile gain
monotonically in K.

Not recommended either, and stated so it is not re-derived: `eps_tail` and `inveff` are the only
arms that look better than a23, at **t = 1.35 / 1.36** on `reach` and **t = 0.88 / 0.75** on
`pds`, and their headline value swings 7.8× with an unpinned calibration (§5). 156 §5's "do not
spend a slot on another member of this family" stands, now with the ordering evidence it lacked.

Three things would reverse this, in order of cost:

1. **A `reach` instrument with an error bar below 0.001 raw.** The blocker is that
   `k*_p / n_real_p` on 85 perturbations is a step function dominated by the small-`n_real`
   tercile. The route that does not need more emissions is a fully closed-form emitted p-order
   (`|F̄_g − ½| / sqrt(tc_g)`, Spearman +0.996 with the realised order, computable from the
   control CDF and the emitted count law) scored through `debench.Surface(key=...)` — that
   removes emission seed noise entirely instead of averaging it down. This experiment's
   seed-averaged order is a cheap approximation and reaches ±0.0021 on `mn` at 12 seeds.
2. **Pinning the h1-test→board `reach` calibration.** §5's factor of 7.81 is the single largest
   uncertainty in the whole write-up, and it is not reducible by more local compute — it needs a
   board point whose `reach` move is large enough to be unambiguous, which is a slot.
3. **A second `reach`-representative surface.** `h1-test` is the only one; repaired HCT116 misses
   a15's board `reach` by 2.6× and is not a substitute.

### What I did NOT measure

* **`fid` and `jac` levels are unreadable here** (rule 9: both carry `n_real`, 89 on h1-test
  against ~240 on the challenge). Across the whole mass-fix family their §3a columns move by
  ≤0.0014 (`fid`) and ≤0.0007 (`jac`) — `a8`'s larger `jac` move is an emission-law change, not
  a placement. Reported for completeness, converted only through the `mn` calibration, and not
  relied on.
* **`nmae` is measured but is not the board's `nmae`**: it is per-perturbation normalised, so
  aggregate-raw arithmetic on it is invalid (`debench.overall_delta` refuses it). The h1-test
  spread across the family is 0.96619–0.96863 (±0.13%), consistent with a23's realised −0.000106,
  and no placement moves it materially.
* **The challenge's own `reach` regime.** Everything is at `n_real` = 89. The one
  proportional-power instrument that reaches the challenge regime (`f = 3.5`) failed its
  known-answer row and was discarded.
* **Nine of the fifteen arms have only 2 emission seeds**, where the `reach` estimate is
  noise-dominated (clustered se 0.0014–0.0038) and same-seed paired deltas flip sign between the
  two draws (`hp_p2000` read +0.00038 at seed A and −0.00429 at seed B). Their §3a/§3b `reach`
  entries and everything derived from them are not evidence. The six arms that carry the
  decision — `a15`, `mn`, `hp_p2000`, `eps_keep`, `eps_tail`, `inveff` — have 12.
* **`min_cpm 10.0` beyond the correction gate.** `mn_c10` tightens the gate the correction may
  use; it does not rebuild the h1-test prediction at `min_cpm 10.0`, which would also change
  which genes carry effect at all. Different arm.
* **Anything downstream of the renderer** — no build, no `vcc prep`, no submission, and no file
  outside this directory modified.

## Files

`recover.py` (invert the frozen artifact → `work/effects.npz`) · `arms.py` (the generalised
`u`-solver, the emission, the DE tables) · `place.py` (what populates the head; writes the frozen
order) · `forced.py` (§2b, the irreducible `log(S/TS)` lift) · `price.py` (six members + board
conversion) · `reso.py` (paired clustered error bars, power matching, `ordAUC`) · `meanorder.py`
(seed-averaged order, channel split) · `profse.py` (profile-member error bars). Logs and `.npz`
under `work/`.

No `.diff` is included: the builder patch was conditional on GO. `mass_consistent_effect` in
`scripts/make_submission_pm_cd4.py` is unchanged.

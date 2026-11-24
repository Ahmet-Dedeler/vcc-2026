# 214 — the two cheap positives, composed on `cd4f034`

**Four results, in order of how much they change what we do next.**

1. **"PerturbMatch is a substitute for the CD4 dose" is REFUTED, and with the opposite sign.**
   Swept on PLOLO over seven doses × two PerturbMatch settings, PerturbMatch's `pds` value
   **RISES** with the dose — `+0.002856` raw at `f = 0.10`, `+0.003797` at `f = 0.34`,
   `+0.004310` at `f = 0.50` — and the interaction is **positive at every dose**
   (`+0.000960` raw at 0.34). They are weak **complements**, not substitutes. And the
   corollary the mandate was really after is a clean negative: **the `pds`-optimal dose is
   `f = 0.20` with PerturbMatch and `f = 0.20` without it.** 207's sweep was not contaminated
   by a redundant component and its optimum is not elsewhere.

2. **PerturbMatch never flipped. Two other things changed underneath it, and both are
   instrument facts rather than science.** Its effect-level `pds` value reproduces 058's
   `+0.00517` on a completely different instrument (**`+0.004785`** on PLOLO's 300-target
   panel against a `+0.00517` measured on a 238-target Orion gallery). What 200 read as
   `−0.001615` is (a) the same component read at the **delivered** level, which costs it 40%
   of its `pds`, and (b) the **first time its four DE members were priced on a
   board-validated surface**. 058, 064 and 069 measured `pds` alone; 061 did measure all six,
   on the **depth-confounded HCT116 harness** AGENTS §6.4 quarantines, and read `nmae`
   **+0.0037 scaled** where h1-test reads **−0.0151**. The sign flip is a harness change, not
   a model change.

3. **Arm A (drop PerturbMatch) is the largest unpriced number in this experiment and rule 13
   vetoes it anyway.** Six-member best estimate **+0.0088 (native) / +0.0097 (wald89)**,
   **+0.0054 / +0.0060** after the measured ×1.533 h1-test deflator, **+0.0046** if `reach`
   transfers at 164's 0.4277. But **every point of it is in `reach` and `nmae`** — rule 13's
   *unpriced* channels — while the one channel rule 13 does price, `pds`, is **negative**
   (`−0.001059`, signs `--0-`). Priced composite **−0.000830** against a **+0.0089** floor
   (the arm moves the delivered amplitude by **+4.28%**). **VETO. Do not spend a slot.**

4. **208 and 200 both under-priced their arms because they were measured on a renderer that
   has not shipped since a23.** `submissions/{a15,191-rc_hybrid,202-rc_callset,210-cd4f034}`
   all record `mass_fix` **absent**; only `a23-massfix` sets it, and a23 scored **−0.00073**
   and was abandoned. On the renderer that actually ships, Hamilton's `pds` is
   **+0.002069** at `f = 0.10` against 208's **+0.001124** — **1.84× larger** — and
   PerturbMatch's `pds` is **+0.001065** against 200's **+0.001885**. Re-price against the
   incumbent's *renderer*, not just its emitter.

**No arm ships.** Rule 13 vetoes all three. Nothing was built, nothing was submitted, no file
outside `experiments/214-cheap-positives/` was created or modified.

> **Status of the incumbent, read 2026-09-04 from the public leaderboard API:** our latest
> published entry is still **`rc_callset`, `31oB9yWy8CSNBDEkY4iW`, 0.1918638367**, submitted
> 2026-09-03T17:47Z. **`cd4f034` is built and prepped but NOT yet scored.** Every price below
> is against the *built* `cd4f034` package; the `jac` model is calibrated on `rc_callset`'s
> board point because that is the last scored one. Where `cd4f034` and `rc_callset` differ
> the difference is measured and stated.

---

## Phase 0 — known-answer rows. **Every gate PASS; the run aborts on any failure.**

### `pds` — PLOLO, 162's pre-committed 4-fold × A/B/C predictor (`report214.py`)

| gate | got | expected | |
|---|---:|---:|---|
| KA1 200's `SHIP` (= a23) per fold — `orion_hct116` | 0.791823 | 0.791823 | PASS (2.3e−07) |
| KA1 `orion_hek293t` / `k562_gwps` / `kolf2` | 0.777375 / 0.789121 / 0.610128 | same | PASS |
| KA1 4-fold mean | **0.742112** | 0.742112 | PASS (1.4e−07) |
| **KA1c the mandate's a15 / a23 pair — it is the `k562_gwps` fold** | a15 **0.795784**, a23 **0.789121** | 0.795784 / 0.789121 | **PASS** (9.1e−08 / 3.1e−07) |
| KA1b state-carry control (`noMassFix` ≡ `g_pm050_f010`) | worst \|Δ\| **0.00e+00** | 0 | PASS |
| KA2 mass fix (a15 − a23) | **−0.004014** scaled | −0.004014 | PASS (2.2e−07) |
| **KA3 PerturbMatch LOO (`noPM` − `SHIP`)** | **−0.011189** scaled | −0.011187 | **PASS** |
| KA4a union whitening LOO | −0.151888 | −0.151892 | PASS |
| KA4b CPM gate 10 vs 20 | −0.004061 | −0.004060 | PASS |
| KA4c CD4 residual `f = 0.10` LOO | −0.015496 | −0.015496 | PASS |
| **KA5 196/208's rounding rows, chalA, a23 renderer** | floor **0.739016**, ham **0.742224**, ∞ **0.742137** | identical | **PASS** (≤3.8e−07) |
| KA6 169's cached chalA rows at m = ∞ (4 folds) | 0.790712 / 0.790896 / 0.771621 / 0.615317 | identical | PASS |

KA3 and KA4a–c are 200's PerturbMatch row and three of its six load-bearing rows
(`+0.001885`, `+0.025593`, `+0.000684`, `+0.002611` overall), reproduced from 200's own
cached per-donor decomposition through this experiment's own composer.

### `mse` — 096's `claim` route on the real challenge panel (`mse214.py`)

| gate | got | expected | |
|---|---:|---:|---|
| `claim` flat 0.25/φ0 | 0.11860 | 0.11860 | PASS (−0.00%) |
| `claim` flat 1.00/φ0 | 0.19438 | 0.19216 | PASS (+1.16%) |
| `claim` flat 0.25/φ4 (a8's law) | 0.16652 | 0.16478 | PASS (+1.06%) |
| `claim` a15 as shipped | **0.24581** | 0.24581 | PASS (+0.00%) |
| A/B round trip at s = 1 | A 0.283614, B 0.281236 | `board_solve` | PASS (+0.000%) |
| a15 raw `mse` / `score_mse` | **0.88228 / 0.11112** | 0.88228 / 0.11112 | PASS |
| **200's own PerturbMatch geometry row** | **+0.000256** | +0.000256 | **PASS** |

### The challenge-side rounding tilt (`chal214.py`)

| profile | rule | mean rel | common log2 | **TILT log2** | rounding energy | `floor(x) == T` |
|---|---|---:|---:|---:|---:|---:|
| a15's own delivered totals | floor | −0.00229 | −0.000894 | **−0.01260** | 3.0905e−06 | **1.0000** |
| a15 | ham | +0.00007 | +0.000030 | **+0.00037** | 9.8522e−07 | 1.0000 |
| **`cd4f034`'s own delivered totals** | floor | −0.00229 | −0.000893 | **−0.01249** | 3.0891e−06 | **1.0000** |
| **`cd4f034`** | ham | +0.00007 | +0.000030 | **+0.00036** | 9.8435e−07 | 1.0000 |

196 expects **−0.01234 → +0.00039**; measured **−0.01260 → +0.00037**. **PASS.** *(One
convention had to be recovered before this gate would pass: the tilt is defined on the
**scorer's** DE universe, `mean per-cell CPM > 5` — `tilt208.py`:36 — not the builder's
pooled 10-CPM gate. Read on the builder's gate it halves to −0.00683, because the low decile
the tilt is defined on is exactly the 5–10 CPM band the builder gate removes.)*

**And the re-measure is the point:** the tilt on `cd4f034`'s *own* delivered profile is
**−0.01249**, within 0.9% of a15's, so 208's transfer of the h1-test `nmae` gain by the
measured challenge/h1-test tilt ratio survives the dose unchanged.

### `floor_tier`'s gate and 208's rider row (`emit214.py`, `six214.py`)

| gate | got | |
|---|---|---|
| `art_floor` ≡ `de202.rc_callset.npz`: `pred_p`, `pred_padj`, `pred_lfc`, `pred_z` | **bit-identical on 85 × 10,779, max \|diff\| 0** | PASS |
| …per-gene panel-total digest | **85/85** | PASS |
| …five members scored beside 202's own | `n_pred` 2634, `reach` 0.223341, `fid` 0.494362, `jac` 0.072975, `nmae` 0.966702 — **equal to 1e−12** | PASS |
| **`art_ham` ≡ 208's rider A** | `n_pred` **2629**, `reach` **0.22340**, `fid` **0.49429**, `jac` **0.07302**, `nmae` **0.96040** | **PASS, all five** |
| 208's marginal-to-`rc_callset` figure | DE half **+0.001705** here vs 208's +0.001702; with 208's own `pds` +0.001195 and `mse` +0.000002 → **+0.002902** vs 208's **+0.002899** | PASS |

### Rule 11 — path and sha256 for every operator (`inputs214.py`)

**14 operators on disk, 13 distinct sha256, exactly one matching the shipped digest.**

| path | sha256[16] | max eig | n | |
|---|---|---:|---:|---|
| **`data/donors/whitener_union_pruned.npz`** | **`e2f65ecf69b7e0b0`** | **5.894570** | **17942** | **SHIPPED** — matches `meta['operator']['sha256']` of a15, `rc_hybrid`, `rc_callset` and `cd4f034` |
| `data/donors/whitener_union.npz` | `5608e807230dea04` | 6.256702 | 17942 | the `--whitener-id` default trap (rule 10b) |
| `experiments/199-def-runbook/…_def199.npz` | `68dd20dc5354f401` | **5.938910** | 17942 | 208's near-collision: same name pattern, same shape, eigenvalue 0.75% away |
| `experiments/095-clean-harness/out/whitener_union.npz` | `5cf3899db9272ff2` | 8.036353 | 17942 | 200's near-collision pair |
| `experiments/095-clean-harness/out-h1test/whitener_union.npz` | `f80fe8d6ba59138b` | 8.080842 | 17745 | every h1-test rebuild here |
| `data/donors/whitener_union_full.npz` | `965d3cef36ac6adc` | 6.980217 | 17942 | |
| `experiments/095-…/out-essential`, `out-random`, `040-…`, `099-…`, `data/donors/whitener_full`, `140-…`×2 | 7 more | — | — | full table in `work/inputs214.json` |

`experiments/096-mse-floor/out-alpha/whitener_union.npz` is byte-identical to
`095-clean-harness/out/whitener_union.npz` (both `5cf3899db9272ff2`), which is why 14 files
carry 13 digests. **35 frozen inputs** are hashed in `work/inputs214.json`, including
`data/external/perturbmatch/delta.npy` = `efa20944f6f6485a…`, the digest 058 recorded.

---

## Phase 1 — what actually changed underneath PerturbMatch

### 1.1 The stack's history, and where PerturbMatch's `pds` value sits at each stage

The whole lineage shares **one effect matrix**. `effects.final_sha256` is `22c09f22c749bda8`
for **a23, `rc_hybrid` and `rc_callset` alike**, and `post_whiten_sha256` `b2426e15d76d` for
a23, `rc_hybrid`, `rc_callset` **and** `cd4f034`. So `rc_hybrid` and `rc_callset` are
*emission* arms and **cannot move PerturbMatch's `pds` price at all**; only the dose and the
renderer can.

| stage | dose | renderer | PerturbMatch's `pds` value, raw | overall |
|---|---|---|---:|---:|
| a15 / `rc_hybrid` / `rc_callset` | `f = 0.10` | shipped (no mass fix) | **+0.002856** | **+0.001065** |
| a23 — **200's `SHIP`, the row in the ablation table** | `f = 0.10` | a23's `minnorm` mass fix | **+0.005058** | **+0.001885** |
| **`cd4f034`** | `f = 0.34` | shipped (no mass fix) | **+0.003797** | **+0.001415** |
| 058/064's instrument (effect level, no renderer) | `f = 0` / `0.05` | — | **+0.004785** | +0.001784 |
| 058 as published (Orion holdouts, 238-target gallery) | — | — | **+0.00517** | — |

**PerturbMatch's `pds` has never been negative on any instrument, at any stage.** 058's
`+0.00517` reproduces as `+0.004785` on a different surface with a different truth, a
different gallery and a different donor set once the *level* is matched. What changed is:

* **the level.** 058, 061, 064 and 069 all scored the **effect**; production is scored after
  `profile_from_effect` (156). Rendering costs PerturbMatch **40% of its `pds`**
  (`+0.004785 → +0.002856`).
* **the renderer.** a23's mass fix recovers most of that (`+0.005058`), and 200 measured
  PerturbMatch on a23. **a23's renderer has not shipped since a23**, which scored −0.00073.
* **the accounting.** Nobody had ever priced PerturbMatch's four DE members. 200 did, and
  they are where the negative lives.

### 1.2 The substitution question — **refuted**

`pds` value of PerturbMatch = `pds(pm 0.50) − pds(pm 0.00)`, 4 PLOLO folds × challenge frames
A/B/C, shipped renderer, γ 0.65, α 0.982:

| CD4 `f` | `pds` pm .5 | `pds` pm 0 | **Δ raw** | fold sd | scaled | overall | /σ | signs | cos pm .5 | cos pm 0 | Δcos % |
|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| 0.00 | 0.732995 | 0.730158 | **+0.002837** | 0.003011 | +0.006275 | +0.001057 | 4.2 | `++0+` | 0.019776 | 0.019199 | +3.01 |
| 0.05 | 0.737448 | 0.734583 | +0.002865 | 0.003348 | +0.006337 | +0.001068 | 4.2 | `++0+` | 0.020871 | 0.020302 | +2.80 |
| 0.10 | 0.740297 | 0.737441 | +0.002856 | 0.003417 | +0.006319 | +0.001065 | 4.2 | `++0+` | 0.021842 | 0.021282 | +2.63 |
| **0.20** | **0.742229** | **0.738941** | +0.003288 | 0.003086 | +0.007272 | +0.001225 | 4.9 | `++0+` | 0.023408 | 0.022869 | +2.36 |
| 0.34 | 0.740255 | 0.736458 | **+0.003797** | 0.003045 | +0.008399 | +0.001415 | 5.6 | `++0+` | 0.024862 | 0.024350 | +2.10 |
| 0.50 | 0.735047 | 0.730737 | **+0.004310** | 0.003219 | +0.009535 | +0.001607 | 6.4 | `++0+` | 0.025723 | 0.025240 | +1.91 |
| 0.80 | 0.724089 | 0.720085 | +0.004004 | 0.003035 | +0.008857 | +0.001492 | 5.9 | `++0+` | 0.026059 | 0.025629 | +1.68 |

The mirror image — the dose's value at each PerturbMatch setting, and the interaction:

| CD4 `f` | Δ`pds` at pm 0.50 | Δ`pds` at pm 0.00 | interaction | overall |
|---:|---:|---:|---:|---:|
| 0.05 | +0.004453 | +0.004425 | +0.000028 | +0.000010 |
| 0.10 | +0.007302 | +0.007283 | +0.000020 | +0.000007 |
| **0.20** | **+0.009234** | **+0.008783** | +0.000451 | +0.000168 |
| 0.34 | +0.007260 | +0.006300 | **+0.000960** | +0.000358 |
| 0.50 | +0.002052 | +0.000579 | +0.001474 | +0.000549 |
| 0.80 | −0.008906 | −0.010073 | +0.001167 | +0.000435 |

**Three statements, all measured:**

1. **The interaction is positive at every dose and grows with it.** PerturbMatch and the CD4
   residual are weak **complements**. That is 064's finding (interaction `+0.000226` at
   `f = 0.05`) confirmed and extended: at `f = 0.05` this experiment reads `+0.000028`, the
   same order and the same sign.
2. **The `pds`-optimal dose is `f = 0.20` with PerturbMatch and `f = 0.20` without it.**
   Dropping PerturbMatch does **not** move the dose optimum, so 207's 25-point sweep was not
   run with a redundant component installed and its optimum is not somewhere else. *(207
   chose `f = 0.34` on the six-member price, not on `pds`; on `pds` alone both settings peak
   at 0.20 and `f = 0.34` costs `−0.00197` raw against that peak.)*
3. **The blend fraction 0.50 is at the optimum at both doses**, so 058's grid choice survives:

   | pm | `pds` at `f = 0.10` | `pds` at `f = 0.34` |
   |---:|---:|---:|
   | 0.00 | 0.737441 | 0.736458 |
   | 0.25 | 0.739968 | 0.739361 |
   | **0.50** | **0.740297** | **0.740255** |
   | 0.75 | 0.738957 | 0.739325 |
   | 1.00 | 0.736534 | 0.737059 |

### 1.3 The three candidate mechanisms, tested in the mandate's order

**The CD4 dose — ruled out.** §1.2: the sign is wrong. Raising the dose makes PerturbMatch
*more* valuable, not less.

**The operator — real, and it is the one that interacts.** PerturbMatch's `pds` value with the
whitener on and off:

| CD4 `f` | γ = 0.65 | γ = 0 | ratio | `pds` (γ 0, pm .5) | `pds` (γ 0, pm 0) |
|---:|---:|---:|---:|---:|---:|
| 0.00 | +0.002837 | +0.002650 | 0.93 | 0.669606 | 0.666956 |
| 0.10 | +0.002856 | +0.002207 | 0.77 | 0.676752 | 0.674544 |
| 0.34 | +0.003797 | +0.000598 | **0.16** | 0.687164 | 0.686565 |

Without the whitener, PerturbMatch is worth essentially **nothing** at the shipped dose. That
is consistent with 069's decision — PerturbMatch belongs on the *prediction* side and must
never enter covariance fitting — and it is a new fact: **the operator is what converts
PerturbMatch's orthogonal K562 processing into retrieval, and it needs a dose to work
against.** *(058's un-whitened control did read `+0.00523`, on its own 6,949-gene common axis
with no CD4 residual at all — the `f = 0.00` row here, `+0.002650`, is the comparable one.)*

**The emission / the renderer — the largest single factor.** PerturbMatch scored at the
effect level, at the delivered level, and at the delivered level through a23's mass fix:

| CD4 `f` | Δ intended (058/064's level) | Δ delivered (shipped) | Δ delivered (a23 mass fix) |
|---:|---:|---:|---:|
| 0.00 | **+0.004785** | +0.002837 | — |
| 0.10 | **+0.004785** | +0.002856 | **+0.005058** |
| 0.34 | +0.004170 | +0.003797 | +0.004482 |

**The renderer costs PerturbMatch 40% of its `pds`, and a23's mass fix gives it back.** 200
measured PerturbMatch through a23's mass fix; the stack has not shipped that flag since a23.

### 1.4 The DE half — where the negative actually lives, and how it re-prices

200's PerturbMatch row is `pds +0.001885`, `mse −0.000256`, `reach −0.001104`,
`fid +0.000075`, `jac +0.000299`, `nmae −0.002514`, sum **−0.001615** (sign convention:
positive = what removing the component would cost). **The `pds` and `mse` halves were never
negative; `reach` and `nmae` are the whole of it, and both are h1-test numbers.**

Re-measured on the **incumbent's emitter** (`rc_callset`'s prescribed call set), three matched
seeds, mass fix off:

| Δ h1-test raw, for **dropping** PerturbMatch | `reach` | sem | `fid` | `jac` | `nmae` | recall |
|---|---:|---:|---:|---:|---:|---:|
| **a15's `quasi` emitter (200)** | **+0.01491** | 0.00695 | — | — | — | — |
| **`rc_callset`/`cd4f034`'s emitter (here, `native`)** | **+0.02612** | **0.00009** | −0.00011 | −0.00065 | **−0.01148** | +0.00945 |
| …`wald89` pool (188's, the mandate's) | **+0.03139** | 0.00029 | −0.00011 | −0.00067 | −0.01059 | +0.00663 |

**PerturbMatch's `reach` cost roughly doubled against the current incumbent and became
resolvable**: the three-seed paired sem fell **77×**, from 0.00695 to **0.00009**, taking the
effect from **2.1σ** to far past any threshold. *(Both sems are three-seed and carry two
degrees of freedom, so the ratio is the trustworthy part, not the σ count.)* The mechanism is
that under `rc_callset` the delivered order **is** the `|H|` key — `rc_emit.py`:983,
delivered = prescribed ∪ unsuppressible — so `reach` stops reading emission noise and starts
reading the key directly. This is the exact inverse of 208's Hamilton finding on the same
emitter: the mechanism that *erased* Hamilton's `reach` half **amplifies** a key-changing
arm's.

**The mechanism, on the reference's own BH pool ordered by our `|H|` (`mech214.py`):**

| arm | P(1) | P(3) | P(10) | sign acc | ρ(\|H\|, ref \|lfc\|) | k = ours/theirs |
|---|---:|---:|---:|---:|---:|---:|
| incumbent geometry (PM 0.50) | 0.9310 | 0.8161 | 0.7138 | 0.5307 | **−0.0344** | 0.2552 |
| **A — drop PerturbMatch** | 0.8621 | 0.8161 | 0.7052 | **0.5379** | **−0.0001** | 0.2461 |
| artifact (095's geometry, PM 0) | 0.8448 | 0.8046 | 0.6931 | 0.5358 | −0.0089 | 0.3003 |

Every column is a mean over the 58 perturbations whose reference pool has ≥ 10 genes; ρ is the
per-perturbation Spearman of `|H|` against the reference's `|lfc|`, averaged, and both values
are small enough that only the *change* should be read. Dropping PerturbMatch makes the *head*
worse (P(1) −0.069) and the *body* better (whole-pool sign accuracy +0.0072) and removes a
small mean **anti**-correlation between `|H|` and the reference's `|lfc|` (ρ −0.0344 →
−0.0001). `reach` integrates the body, not the head, which is
`docs/BOARD-LEDGER.md` §5's caution ("three quarters of our `reach` comes from
between-perturbation heterogeneity") showing up as a sign disagreement between two statistics
of the same key.

The `nmae` half is **not** the aggregate magnitude ratio: `k` moves `0.2552 → 0.2461`, and at
`q = 0.537` the `nmae(k,q)` identity's slope below `k = 1` is `−0.074`, predicting
`Δnmae = +0.00068` — **17× too small and of the wrong sign** against the measured `−0.01148`.
That is 187's finding reproduced: the uniform-`k` identity is an artefact, so this is a
per-gene shape effect and has to be read off the measurement.

### 1.5 The instrument correction nobody had applied — **h1-test overstates PerturbMatch by ×1.533**

h1-test carries **two** production donors (`k562_gwps` + `orion_hek293t`; `h1_2025` is oracle
leakage, `orion_hct116` is not in 095's h1-test operator). PerturbMatch enters the **K562 slot
only**, so it owns **1/2** of a two-donor consensus and **1/3** of the shipped three-donor one.
Measured rather than argued (`disp214.py`):

| surface | donors | targets | PM-covered | ‖ΔE‖/‖E‖ | row median | row cos median |
|---|---:|---:|---:|---:|---:|---:|
| challenge (shipped) | 3 | 300 | 251 | **0.27059** | 0.24704 | 0.97081 |
| h1-test | 2 | 85 | 79 | **0.41486** | 0.36445 | 0.93314 |

*(`cd4f034`'s own meta records `perturbmatch_replacement_targets = 248` against the 251 E062's
`prepare_main_source` reports here; the three-target gap is between two both-supported tests
and is immaterial to the ratio.)*

**Ratio 1.533**, against the 1.500 the 1/2-vs-1/3 argument predicts. **Every h1-test reading of
PerturbMatch — 200's and this experiment's alike — overstates the shipped arm by that factor**,
and it is applied as a sensitivity throughout §2. `pds` needs no such correction: PLOLO runs on
the real 300-target panel with all three donors.

---

## Phase 2 — the composed price, marginal to `cd4f034`

### 2.1 `pds` — PLOLO, 4 folds, challenge frames A/B/C

Incumbent = `f = 0.34`, pm 0.50, mass fix off, `floor` rounding — `cd4f034` exactly.

| arm | rule | `pds` | Δ raw | fold sd | scaled | **overall** | /σ | signs | cos | Δcos % |
|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|
| incumbent `cd4f034` | floor | 0.735996 | +0.000000 | — | — | +0.000000 | 0.0 | `0000` | 0.024823 | +0.0000 |
| **A drop PerturbMatch** | floor | 0.733154 | **−0.002842** | 0.002458 | −0.006287 | **−0.001059** | −4.2 | `--0-` | 0.024314 | **−2.0512** |
| **B Hamilton rounding** | ham | 0.740111 | **+0.004115** | 0.002180 | +0.009102 | **+0.001534** | +6.1 | `++++` | 0.024849 | +0.1056 |
| **A∘B both** | ham | 0.736332 | +0.000336 | 0.002423 | +0.000742 | **+0.000125** | +0.5 | `+-+-` | 0.024340 | −1.9475 |
| *(reference)* no rounding at all | ∞ | 0.740255 | +0.004259 | 0.002167 | +0.009420 | +0.001587 | +6.3 | `++++` | 0.024862 | +0.1583 |

Per fold, **including the leak-free `kolf2`** (falsifier 3):

| arm | `orion_hct116` | `orion_hek293t` | `k562_gwps` | `kolf2` |
|---|---:|---:|---:|---:|
| incumbent `cd4f034` | 0.777958 | 0.749592 | 0.802447 | 0.613989 |
| A drop PerturbMatch | 0.776216 | 0.746662 | 0.802447 | 0.607291 |
| B Hamilton | 0.783908 | 0.754618 | 0.802840 | 0.619079 |
| A∘B | 0.781931 | 0.749387 | 0.802840 | 0.611169 |

`k562_gwps` is the `0` sign in `--0-`: with K562 held out, PerturbMatch has no slot to occupy,
so that fold is a structural tie, not a null result. **Arm B is `++++` on all four folds and
captures 97% of the `m = ∞` rounding ceiling.**

**Δcos is the direction test 204 sized.** Arm B moves the transfer cosine by **+0.106%
relative**, inside 204's ±0.20% resolution — **Hamilton is not a direction gain**, exactly as
208 found. Arm A moves it by **−2.05% relative**, ten times outside that floor: **dropping
PerturbMatch is a real direction loss.**

**Arm B re-priced, and 208's number corrected.** Hamilton's `pds` at each (dose, renderer):

| dose, renderer | frames | floor | ham | Δ raw | overall |
|---|---|---:|---:|---:|---:|
| `f 0.10`, a23's mass fix — **208's setting** | chalA | 0.739016 | 0.742224 | +0.003207 | **+0.001195** |
| `f 0.10`, a23's mass fix | A/B/C | 0.739018 | 0.742033 | +0.003015 | +0.001124 |
| `f 0.10`, **shipped renderer** | A/B/C | 0.734671 | 0.740222 | +0.005551 | **+0.002069** |
| `f 0.34`, a23's mass fix | A/B/C | 0.738653 | 0.740417 | +0.001764 | +0.000658 |
| **`f 0.34`, shipped renderer — `cd4f034`** | A/B/C | 0.735996 | 0.740111 | **+0.004115** | **+0.001534** |
| `f 0.34`, shipped renderer, PerturbMatch OFF | A/B/C | 0.733154 | 0.736332 | +0.003178 | +0.001184 |

**The growth is the renderer, not the dose.** a23's `minnorm` mass fix already absorbs about
half of `floor`'s systematic tilt, so measuring Hamilton through it under-states the arm that
ships by **1.84×**. The dose then takes some back (`+0.002069 → +0.001534`).

### 2.2 `mse` — 096's `claim` route, and the cosine channel separately

**Arm B — emission, `claim` on `cd4f034`'s own delivered profile with the incumbent's emitter:**

| rule | Ĉ/C_p | rows < cap | `claim` | branch | Δclaim | round E | Δprof | raw `mse` | score | **Δ overall** | nnz/blk |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| floor | 2.955 | **0 / 24** | 0.24653 | **SATURATED** | −0.00000 | 3.078e−06 | +0.00000 | 0.88228 | 0.11112 | +0.000000 | 1,733,290 |
| **ham** | **2.973** | **0 / 24** | **0.24653** | **SATURATED** | −0.00000 | **9.895e−07** | −0.00001 | 0.88227 | 0.11113 | **+0.000002** | 1,731,780 |

`claim` stays hard against its cap on both sides, with **zero** blocks under it, reproducing
208's `+0.000002` exactly. The PLOLO cross ratio at the incumbent's own configuration is
**floor 1.912 / ham 0.991** (per fold 2.610/1.523/1.853/1.661 and 1.192/0.921/0.938/0.911) —
196's mechanism intact: `floor`'s rounding error adds ~1.9× its own energy to the squared
error because it is *correlated* with the error already there; Hamilton's adds 0.99×, i.e. it
is uncorrelated noise.

**Arm A — geometry, `A`/`B`/`claim` as ratios with `rX = 1`:**

| pair | rA | rB | rClaim | raw `mse` | score | Δ overall |
|---|---:|---:|---:|---:|---:|---:|
| 200's own PerturbMatch row (a23 renderer, `f 0.10`) — **KA** | 1.03276 | 1.03447 | 1.00028 | 0.88083 | 0.11265 | **+0.000256** |
| PerturbMatch, shipped renderer, `f 0.10` | 1.03218 | 1.03402 | 1.00060 | 0.88069 | 0.11280 | +0.000280 |
| **ARM A, `f 0.34`** | **1.04283** | **1.04448** | 1.00014 | 0.88098 | 0.11249 | **+0.000229** |

`rX = 1` deliberately refuses to credit any direction change, so this is a **partial** price.
The direction half comes from PLOLO's cosine through 212's own measured exchange rate
`d(overall mse*)/d(cos) = +0.3655`: arm A's `Δcos = −0.000509` → **−0.000186**. **Net for arm
A: `+0.000229 − 0.000186 = +0.000043`** — a wash, with the two channels nearly cancelling. The
`A`-move is nearly free because **`B` tracks `A`** (1.04448 vs 1.04283), which is 211's finding
about `cd4f034` reappearing here.

### 2.3 The four DE members — h1-test, `rc_callset` emitter, three matched seeds

Every geometry arm is renormalised to the artifact's own median row norm (`h1parts200`'s
convention), so on this surface arm A is a **pure shape** comparison with no amplitude
confound. `n_pred` is **2600 on every arm and every seed** — the call set is prescribed.

| pool | arm | `reach` | sem | `fid` | `jac` | `nmae` | recall |
|---|---|---:|---:|---:|---:|---:|---:|
| native | incumbent | 0.19437 | — | 0.49382 | 0.07384 | 0.97298 | 0.48319 |
| native | **A** | **0.22042** | 0.00009 | 0.49371 | 0.07319 | **0.96151** | 0.49264 |
| native | **B** | 0.20647 | 0.00009 | 0.49404 | 0.07385 | 0.97195 | 0.48319 |
| native | **A∘B** | **0.22955** | 0.00038 | 0.49369 | 0.07319 | **0.96064** | 0.49263 |
| wald89 | incumbent | 0.18118 | — | 0.49382 | 0.07798 | 0.97577 | 0.55637 |
| wald89 | **A** | **0.21242** | 0.00029 | 0.49371 | 0.07731 | **0.96518** | 0.56300 |
| wald89 | **B** | 0.18024 | 0.00055 | 0.49404 | 0.07798 | 0.97609 | 0.55637 |
| wald89 | **A∘B** | 0.21002 | 0.00079 | 0.49369 | 0.07730 | 0.96542 | 0.56299 |

`matched` reproduces `native` to four decimals. The **artifact** family — the surface whose
`floor` row is bit-identical to `rc_callset` — reads `art_floor` 2634/0.22334/0.49436/0.07297/
0.96670 and `art_ham` 2629/**0.22340**/0.49429/0.07302/**0.96040**, i.e. 208's rider A exactly.

> **A disagreement this experiment is required to state.** Arm B's `reach` is **+0.01208 raw
> on the rebuild family** and **+0.00006 on the artifact family**, and its `nmae` is **−0.00103
> on the rebuild** against **−0.00630 on the artifact**. The rebuild's median row cosine to the
> shipped `H` is **0.813** (200 measured the same), so it is a *reconstruction* of the geometry
> and sits at a different point on `nmae`'s per-gene curve. **Arm B's price is taken from the
> ARTIFACT family** — that is the surface with a bit-identical known-answer row against the
> shipped package, and it is 208's own instrument. The rebuild family is used only for arm A,
> which cannot be expressed on an artifact at all, and for the composition's sign.

### 2.4 The six-member price (overall board units, marginal to `cd4f034`)

`reach` × **1.290** (`rc_hybrid`'s own measured h1-test → board factor for an order-changing
arm), with 164's **0.4277** as the sensitivity; `fid` as precision (challenge `n_pred` 2,599 ≫
`n_real` ≈ 222, so the yield cap is saturated); `jac` through 208's surviving S\* convention
with board `NP` pinned at **2,599** and `NR` = **227.3** calibrated on `rc_callset`'s own board
`jac`; `nmae` h1-test raw / span.

| member | instrument | **A** | **B** | **A∘B** |
|---|---|---:|---:|---:|
| `pds` | PLOLO 4-fold × 1.011 | **−0.001059** | **+0.001534** | +0.000125 |
| `mse` | 096 `claim`, A/B ratios | +0.000229 | +0.000002 | +0.000231 |
| `mse` | cosine channel (212's rate) | −0.000186 | +0.000010 | −0.000177 |
| `reach` | h1-test × 1.290 (`wald89`) | **+0.007621** | −0.000000 | +0.006939 |
| `fid` | h1-test precision | −0.000060 | +0.000123 | −0.000072 |
| `jac` | S\*, `NP` pinned 2,599 | +0.000265 | −0.000000 | +0.000265 |
| `nmae` | h1-test raw / span | **+0.002909** | −0.000086 | +0.002843 |
| **SUM, `wald89`** | | **+0.009719** | +0.001581 | **+0.010153** |
| **SUM, `native`** | | **+0.008791** | +0.004884 | **+0.012318** |
| SUM, `matched` | | +0.008781 | +0.004880 | +0.012307 |
| SUM, `reach` at 0.4277 (`wald89`) | sensitivity | +0.004625 | +0.001582 | +0.005515 |
| **SUM, DE deflated ×1/1.533** (`wald89`) | §1.5 | **+0.005986** | +0.001581 | **+0.006697** |
| SUM, DE deflated ×1/1.533 (`native`) | | +0.005380 | +0.004884 | +0.009258 |
| **SUM with arm B from the ARTIFACT row** | the calibrated surface | — | **+0.003249** | — |

**Composition, measured and not inferred (178's rule).** `A + B = +0.011300` (`wald89`) against
a **measured** `A∘B = +0.010153`: **sub-additive by −0.001147**. On `native`, `+0.013675`
against `+0.012318`, **−0.001356**. On the `pds` channel alone the sub-additivity is starker —
`A + B = +0.000474` against a measured `+0.000125`, `−0.000349` — because both arms move the
same delivered profile. **178's "gains compose sub-additively" holds here, unlike 208's pair,
and A∘B is worse than the naive sum in every pool.** It is not, however, worse than A alone.

---

## Phase 3 — the gate

**Arm classes, measured:**

| arm | digest test | delivered `A` ratio | median ‖row‖ ratio | class | toll | floor |
|---|---|---:|---:|---|---:|---:|
| **A** drop PerturbMatch | not an emission change | **1.04283** | 1.03099 | profile, **AMPLITUDE-CHANGING** | **0** | **+0.0089** |
| **B** Hamilton rounding | **0/900 identical**, max \|ΔT\| = 1, **UP on 6,807–7,318 genes, DOWN on 0** | 0.99513 | 0.99724 | **emission-law, profile-CHANGING** | **−0.002433** | +0.0019 |
| **A∘B** | 0/900 | 1.03775 \* | 1.02814 \* | both | −0.002433 | +0.0089 |

\* the A∘B amplitude row is the product of the two measured ratios, not an independent
measurement; the two are near-orthogonal (arm A moves the effect, arm B moves the integer
totals by ≤ 1 count each) but the product is an estimate and is labelled as one.

Arm A is a **donor-blend** change, i.e. rule 13's profile class, for which 169 falsified the
`reach` toll at `+0.0006`; it carries **no toll**. But it moves the **delivered** amplitude by
**+3.10% on the median row norm and +4.28% on `A`** (`head214.py`, five nested gene sets,
consistent on four of them; on the *top 2,000 by control CPM* it moves the other way, ×0.943,
which is itself a shape statement). That is **twice** the −2.03% that made 210 flag the
classification, and **4×** the 1.0% bound 210's mandate used, so **arm A is
amplitude-changing and its floor is +0.0089.** Per **211's standing rule α is NOT re-solved**:
amplitude is closed, the `A` move is already priced inside the claim route at `+0.000229`, and
re-solving α would make this a two-variable arm (rule 5). What it costs is stated rather than
recovered: arm A pushes the amplitude *toward* the `s²B`-vs-`claim` cliff, cutting the margin
from 5.64% to 1.16% in `s²`, i.e. the α headroom from +2.78% to +0.58%.

Arm B alters the per-gene panel totals themselves — 0/900 blocks identical, every move `+1` and
every move upward, exactly as 208 measured — so it **carries the −0.0146 toll and claims no
exemption**.

```
predicted = PLOLO d_pds/6 x 1.011  +  096 d_mse/6  -  0.0146/6 (emission-law, profile-changing)
```

| arm | priced `pds` | priced `mse` | toll | **composite** | floor | **verdict** | six-member best estimate |
|---|---:|---:|---:|---:|---:|---|---:|
| **A** drop PerturbMatch | **−0.001059** | +0.000229 | 0 | **−0.000830** | 0.0089 | **VETO** | +0.0060 … +0.0097 |
| **B** Hamilton rounding | +0.001534 | +0.000002 | −0.002433 | **−0.000898** | 0.0019 | **VETO** | +0.0032 (artifact row) |
| **A∘B** both | +0.000125 | +0.000231 | −0.002433 | **−0.002077** | 0.0089 | **VETO** | +0.0067 … +0.0123 |

Arm B also fails **without** the toll (`+0.001536` against `+0.0019`), which is 208's verdict
reproduced on a corrected, larger `pds`.

**All three arms are VETOED. No slot. Nothing is built.**

**And the honest statement of the tension, because it is the whole decision.** Arm A's
six-member estimate is the largest positive this experiment found, and **100% of it sits in
`reach` and `nmae`** — the two members AGENTS §11.1 rule 13 records as *not locally priceable*
(`reach` 2/6 signs, worse than chance; `nmae` 3/6) — while the one member that **is** priceable
6/6 on sign says the arm is **negative**. Rule 13's `+0.0089` floor **is** the realised RMS of
exactly those unpriced channels. Shipping arm A would be spending a slot on a number rule 13
was written to distrust, on a surface that overstates the component by a measured ×1.533, at a
`reach` transfer factor whose plausible range spans **0.39 to 1.76** and which was calibrated
on a *different arm class*. At 0.4277 the arm is `+0.0046`; deflated it is `+0.0060`; both are
below the floor. **The correct action is not to ship it and not to drop PerturbMatch.**

---

## Phase 4 — falsifiers

**F1 — the identity rows.** With both flags off this pipeline must reproduce `cd4f034`. A full
package rebuild is outside this mandate, so the identity is asserted on **every instrument**,
per-gene rather than on scalars:

| surface | identity | result |
|---|---|---|
| h1-test DE | `art_floor` ≡ `de202.rc_callset.npz` | **bit-identical, max \|diff\| 0 on 85 × 10,779 × 4 fields**; totals digest 85/85 |
| challenge profile | `floor(x) == T` on `cd4f034`'s own delivered totals | **1.0000 of sampled blocks**, all three contexts — the reconstruction of the intended profile re-delivers the incumbent's own integers |
| challenge artifacts | `cd4f034` vs `rc_callset` per-gene panel-total digest | **0/900 identical** (210 falsifier 3a: must be 0); first block `A\|ABCD1` `2cdbc2f0c36ba28b` vs `rc_callset`'s `bb6befc01caf7493`, the value 202/210 recorded |
| PLOLO | the incumbent's config reproduces 200's `SHIP`/a15 rows and 169's cached chalA rows | 11 gates, worst \|Δ\| **3.8e−07** |
| state-carry | two field-identical arm dicts | worst \|Δ`pds`\| **0.00e+00** |

**F2 — arm A's sign DOES depend on the instrument, and the two disagree.** Reported separately
as the mandate requires:

| instrument | what it reads | arm A |
|---|---|---:|
| **PLOLO** (300-target panel, 3 donors, challenge frames, delivered) | `pds` + transfer cosine | **−0.001059** overall, Δcos **−2.05%**, signs `--0-` |
| **h1-test** (85 targets, 2 donors, `rc_callset` emitter) | `reach`/`fid`/`jac`/`nmae` | **+0.0098 … +0.0107** overall |
| 200's own surface (h1-test at the `quasi` emitter, a23 renderer) | six members | +0.001615 |

**They disagree in sign, and the disagreement is structural, not noise**: PLOLO reads
*retrieval* (which PerturbMatch helps — the effect is 4.2–6.4× 162's leave-one-arm-out RMS at
every dose, with the same three-of-four fold signs throughout) and h1-test reads *sign ordering
and magnitude* (which it hurts — `+0.02612` raw `reach` against a three-seed sem of `0.00009`).
Both are measured with their own known-answer rows passing. **Per AGENTS §6.4 no absolute is quoted from the Orion
pair anywhere in this write-up**; 058's `+0.00517` appears once, as a published historical
value, and is compared only to a *relative* number.

**F3 — all four PLOLO folds, `kolf2` included**, are in §2.1. `kolf2` is the leak-free fold and
it is the largest single mover for arm A (0.613989 → 0.607291, −0.0067) and for arm B
(0.613989 → 0.619079, +0.0051).

**F4 — the composition is measured, not inferred**, and it is **sub-additive in every pool**
(§2.4). **A∘B is not worse than A alone** on any pool (`wald89` +0.010153 vs +0.009719;
`native` +0.012318 vs +0.008791), but it is worse than `A + B` by −0.0011 to −0.0014, and on
the `pds` channel alone the two arms eat 74% of each other (`A + B = +0.000474` against a
measured `+0.000125`, i.e. `−0.000349`).

**F5 — the format contract**, from the artifacts themselves:

| artifact | nnz | % of 4.75e9 | max counts/cell | % of 1e6 |
|---|---:|---:|---:|---:|
| `cd4f034` (built) | 1,577,729,493 | **33.2154** | 22,323 | 2.232 |
| `rc_callset` (scored) | 1,578,026,567 | 33.2216 | 22,323 | 2.232 |

Hamilton adds **6,319,146 counts** across the package (**+0.0878%** of the panel total
`S = 7,199,453,132 → 7,205,772,278`) and moves the emitted nnz by **−0.09%** on the challenge
blocks sampled here (208 measured **+0.04%** on a15's profile). Both caps clear by **3.0×** and
**45×** either way.

---

## What this changes for other experiments

* **200's ablation table was measured on a renderer that does not ship.** Its `SHIP` row is
  **a23** (`mass_fix = "minnorm"`); a15, `rc_hybrid`, `rc_callset` and `cd4f034` all record the
  flag **absent**. On the shipped renderer PerturbMatch's `pds` is `+0.001065`, not `+0.001885`
  — a **1.77× overstatement** — and the same correction applies in principle to every row of
  that table. It does not change any of 200's verdicts (the six load-bearing rows are 10–25×
  their gate) but it does change the sizes.
* **208's Hamilton price is 1.84× too small for the same reason.** `+0.002069` overall at
  `f = 0.10` on the shipped renderer against 208's `+0.001124`; `+0.001534` at `cd4f034`'s
  dose. It still does not clear `+0.0019`, so 208's verdict stands.
* **h1-test overstates any PerturbMatch reading by ×1.533, measured.** Two donors instead of
  three puts 1/2 of the consensus in the K562 slot instead of 1/3. This is a general property
  of the surface and applies to any future K562-side arm read on it.
* **`reach` under `rc_callset` is a near-noiseless read of the ordering key.** Three matched
  seeds give a `reach` sem of **0.00009** raw on an effect of 0.026, against 200's 0.00695 on
  the `quasi` emitter. **Any key-changing arm is now 77× cheaper to resolve than it was**, and
  that is probably worth more than either arm here: the ordering channel (§8.1) finally has a
  low-variance instrument.
* **207's CD4 dose sweep is clean.** The `pds`-optimal dose is `f = 0.20` with and without
  PerturbMatch, so the sweep was not run against a redundant component.
* **The CD4 residual is invisible on h1-test, re-proved not assumed.** The `f = 0.34` and
  `f = 0.10` h1-test effect matrices are **identical, max \|diff\| 0.000e+00**, because
  production's median-based sizing returns fraction 0 on a surface with 4 of 85 targets covered.
  Any h1-test result for `rc_callset` transfers to `cd4f034` unchanged — and no h1-test screen
  can ever price the dose.

## What could not be measured

* **Arm A's four DE members on the challenge.** They are h1-test numbers deflated by a measured
  ×1.533 donor-slot factor; the deflator corrects the *displacement*, not the *response*, and
  `reach`'s response to a key change is not linear in the displacement. A challenge-side read
  would need the truth we do not have.
* **`cd4f034`'s board members.** The package is built and prepped but unscored, so the `jac`
  model is calibrated on `rc_callset`'s board point and every "marginal to `cd4f034`" figure is
  marginal to a *predicted* incumbent.
* **Arm B on the rebuild geometry vs the artifact.** The two surfaces disagree by 200× on
  `reach` and 6× on `nmae`; only the artifact one has a bit-identical known-answer row, so it is
  the one used. A rebuild with a higher cosine to the shipped `H` than 0.813 would settle it.
* **Arm A against `rc_hybrid`.** Not measured; `rc_hybrid` is superseded and its call set is
  different from `rc_callset`'s, so its DE numbers would not transfer anyway.

## Handoff — the exact changes, and a conditional preregistration

**Neither arm is "one flag", and the mandate's premise that arm A is "one flag, no new code"
is wrong.** `PM_BLEND = 0.50` is a **hard-coded module constant** at
`scripts/make_submission_pm_cd4.py:572` with **no environment override** — the builder exposes
`VCC_ALPHA`, `VCC_ALPHA_REF`, `VCC_CD4_TARGET_F`, `VCC_CD4_FRACTION`, `VCC_DONOR_WEIGHTS`,
`VCC_MASS_FIX`, `VCC_RESAMPLE_INVERT`, `VCC_RESAMPLE_X`, `VCC_SIGN_REPAIR`, `VCC_WHITEN_GAMMA`
and `VCC_WHITEN_LAM`, and nothing for the blend. Both arms are source edits, both move
`hashes.implementation`, and arm A additionally touches the pinned-support gate at
`make_submission_pm_cd4.py`:3521 (`perturbmatch_replacement_targets` is frozen at **248** and
the builder **raises** if the observed count moves).

```python
# ARM A  scripts/make_submission_pm_cd4.py:572
- PM_BLEND = 0.50
+ PM_BLEND = float(os.environ.get("VCC_PM_BLEND", "0.50"))   # then VCC_PM_BLEND=0.0

# ARM B  src/vcc2026/donor_transfer.py:380   (mode="quasi" -- the ONE live call site of five;
#                                             277 is resample_mix's, 498/540/567 are off-path)
- totals = np.floor(n * L * profile_p).astype(np.int64)
+ totals = hamilton(n * L * profile_p)      # largest remainder; experiments/196/ham196.py
```

Arm B's edit is 208's, unchanged; `rc_hybrid_block` (`rc_emit.py`:875) routes through the same
call site, so one line moves a15, `rc_hybrid`, `rc_callset` and `cd4f034` together and
preserves the digest identity *between* them.

**Neither is recommended.** Rule 13 vetoes all three arms and the priced channel is negative
for arm A. What follows is a preregistration for a coordinator who chooses to override the
gate on arm A anyway — written now so the slot buys a falsifiable number rather than a hope.

**Preregistered per-member board prediction for ARM A, marginal to `cd4f034`** (which is
itself unscored; if `cd4f034` has not landed, subtract 207's predicted `+0.00437` first):

| member | prediction | basis |
|---|---:|---|
| `pds` | **−0.00106 ± 0.00055** overall (scaled −0.00636, raw ≈ −0.00284) | PLOLO 4-fold × 1.011, σ 0.000549, signs `--0-` |
| `mse` | **+0.00004 ± 0.00020**; raw `mse` within **0.0015** of 0.88228 | 096 claim `+0.000229` minus the cosine channel `−0.000186` |
| `reach` | **+0.0050 ± 0.0030** overall (raw ≈ +0.0056) | h1-test `+0.0261` × 1.290 ÷ 1.533; the ±spans 0.4277 to 1.76 |
| `nmae` | **+0.0019 ± 0.0010** overall (raw ≈ −0.0069) | h1-test `−0.01148` ÷ 1.533 |
| `fid` | **−0.00004 ± 0.0003** | h1-test precision, `n_pred` pinned at 2,599 |
| `jac` | **+0.0002 ± 0.0004**; challenge `n_pred` unchanged at **2,599** | the call set is prescribed and its size is asserted in code |
| **overall** | **+0.0060 ± 0.0035** | deflated `wald89` sum |

**Falsified if:** board `pds` falls by more than **0.006 raw** (2× the prediction — the arm
moves the profile, so `pds` is a real exposure); or raw `mse` moves by more than **0.005** in
either direction (the `A` move of +4.28% is the exposure and it should be near-free because
`B` tracks `A` to 0.16%); or challenge `n_pred` leaves **2,599 ± 1%**; or board `reach` moves
by less than **+0.002 scaled** (which would mean the ×1.290 factor does not apply to a
key-changing *geometry* arm and the whole DE half is an h1-test artefact).

**Pre-build assertions**, all checkable before a slot is spent:

| # | assertion | expected |
|---|---|---|
| 1 | rule-10 config diff vs `cd4f034` | **exactly one** key, `perturbmatch_blend` 0.5 → 0.0 |
| 2 | `meta['operator']['sha256']` | **`e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`** |
| 3 | `effects.whitening.after_all_median_norm` | **4.229295604454299**, drift < 1e−12 |
| 4 | `effects.cd4.fully_supported_targets` | **276** |
| 5 | `effects.post_whiten_sha256` | **must differ** — the arm changes the pre-whitening consensus |
| 6 | delivered median row norm vs `cd4f034` | **1.031 ± 0.005** (`head214.py`) — the arm is amplitude-changing and this is the number that says so |
| 7 | per-gene panel-total digest vs `cd4f034` | **0/900 identical** |
| 8 | challenge `n_pred` (`probe210.py`, 40 × 3) | **2,599 ± 1%** |
| 9 | package nnz | **1.578e9 ± 0.5%**, i.e. **33.2%** of the 4.75e9 cap |
| 10 | pinned-support gate | `perturbmatch_replacement_targets` must still read **248**, or the builder aborts |

## Reproduction

```
uv run python inputs214.py   # rule 11: 14 operators + 35 frozen inputs -> work/inputs214.json
uv run python pds214.py      # PLOLO grid, 36 arms x 4 folds x 3 frames -> work/pds214.json   (234 s)
uv run python report214.py   # the 11 pds known-answer gates + Phase 1  -> work/report214.json
uv run python chal214.py     # 900-block totals, tilt, class test, caps -> work/chal214.json  (30 s)
MAXT=8 uv run python mse214.py  # 096 claim route, both arms, 7 gates   -> work/mse214.json   (564 s)
uv run python emit214.py     # h1-test, 15 arms x 85 targets, rc_callset-> work/de214.*.npz  (2034 s)
uv run python six214.py      # 4 DE members, 3 pools, 2 KA rows         -> work/six214.json
uv run python head214.py     # call-set exposure + delivered amplitude  -> work/head214.json
uv run python disp214.py     # the h1-test / challenge deflator         -> work/disp214.json
uv run python mech214.py     # why PerturbMatch costs reach and nmae    -> work/mech214.json
uv run python price214.py    # the six-member price and rule 13's gate  -> work/price214.json
```

**Resources**, `/usr/bin/time -l` maximum resident set size: `mse214` **7.78 GiB** (heaviest),
`emit214` 6.50 GiB, `pds214` 6.28 GiB, `head214` 4.93 GiB, `disp214` 4.68 GiB, `chal214`
1.10 GiB — all under the 12 GB budget. No `vcc`, no `cell-eval2`, no submission build, nothing
written to `submissions/`, and no file outside `experiments/214-cheap-positives/` created or
modified. `src/` and `scripts/` were read only.

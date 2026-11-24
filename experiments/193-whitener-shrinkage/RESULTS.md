# 193 — The whitening operator's regulariser. **The ridge is optimal. The geometry is CLOSED.**

**Date:** 2026-09-03. Measurement only. No `vcc`, no `cell-eval2`, no build, no submission,
nothing written outside this directory. Peak RSS **9.2 GiB** (`structured193.py`), 6.3 GiB
(`screen193.py`), against a 14 GiB budget. Two `structured193.py` runs were deliberately NOT run
concurrently — 2 x 8.8 GiB would have breached it.

---

## Decision, in one line

**No covariance regulariser beats production's, on either side, and the premise that `C` is fit
on too little data is false — `pds` falls monotonically across the three union-axis operators as
the basis grows from 29,394 to 45,531 rows.** Fifty-five eigenvalue maps (54 non-identity) were
screened on PLOLO across four donor folds, plus a diagonal-target family that needs its own
eigendecomposition. Excluding one map that turns out to be the `gamma` axis in disguise (§2), the
best of them clears **+0.75 σ** of the instrument's own noise; the best composed six-member price
over all of them is **+0.00051** against a **+0.0019** floor; and every named estimator from the
literature — Ledoit–Wolf linear, Ledoit–Wolf analytical nonlinear, Marchenko–Pastur clipping, the
Ledoit–Péché RIE, and shrinkage toward a diagonal target — is negative. The one reading that came
back positive (`+1.24 σ`, one fold) inverted to `−0.0235` scaled when replicated on a second. Production's `f(ev) = ev^(-0.325)` with a
`1e-12` relative floor and **no ridge at all** sits at an interior optimum, verified on both
sides of the map-shape family.

---

## 0. Phase 0 — known-answer rows. **Six gates, five of them exact.**

Rule 11 in force throughout: every operator is loaded by path and asserted by sha256.

| # | gate | expected | measured | |
| --- | --- | ---: | ---: | :--- |
| 1 | map `prod(lam=0)` reproduces 162's stored `a15` rows, 4 folds × 3 challenge frames | — | worst \|Δ`pds`\| = **`0.00e+00`**, 12/12 bit-identical | PASS |
| 2 | `a15 → a23` mass fix, recomputed end to end through this experiment's code | `+0.004014` scaled (162/169) | **`+0.004014`** (raw `+0.001815`), \|diff\| `4.9e-07` | PASS |
| 3 | `union → union_pruned` at v11's config on `kolf2` — the operator-SWAP path | `+0.001360` raw (162 `a2-prune`) | **`+0.001360`**, \|diff\| `1.6e-07` | PASS |
| 4 | `rho = 0` structured arm: reconstruct `C` from the frozen operator, re-solve its eigensystem, re-whiten — on TWO folds | `0.611201` / `0.787706` (each frozen operator's own row) | **`0.611201`** / **`0.787706`**, \|d\| = `0.00e+00` | PASS |
| 5 | h1-test identity map reproduces 168d/173's base `pds` | `0.86134` | **`0.86134`**, \|d\| `4.5e-06` | PASS |
| 6 | `claim188.py`'s four `claim` gates, verbatim, before any `mse` arm | 0.11860 / 0.19216 / 0.16478 / 0.24581 | see §6 | — |

Gate 2 is the non-vacuous one for the renderer (it exercises `mass_consistent_effect`, which a
pure operator swap never touches); gate 3 is the non-vacuous one for the operator swap; gate 4
is the non-vacuous one for the eigen-solve; gate 5 for the h1-test transfer.

### Operator inventory used here — path, sha256, and what it was fitted on

| sha256 (first 16) | path | basis rows × cols | note |
| --- | --- | ---: | --- |
| `e2f65ecf69b7e0b0` | `data/donors/whitener_union_pruned.npz` | **29,394 × 17,942** | **ships (a15)**; 12,000/donor cap then 091's essential prune |
| `5608e807230dea04` | `data/donors/whitener_union.npz` | 33,787 × 17,942 | 12,000/donor cap, no prune. v11's operator |
| `965d3cef36ac6adc` | `data/donors/whitener_union_full.npz` | 45,531 × 17,942 | **no cap**: every donor row (9,787 + 17,767 + 17,977) |
| `13d988fcc3b129d2` | `062/…/held_k562_gwps.npz` | 35,744 × 17,806 | leakage-safe for the `k562_gwps` fold |
| `7cecbf7f36e40d22` | `062/…/held_orion_hct116.npz` | 27,764 × 17,942 | leakage-safe for the `orion_hct116` fold |
| `97e486e0d1a805e1` | `062/…/held_orion_hek293t.npz` | 27,554 × 17,942 | leakage-safe for the `orion_hek293t` fold |

**Correction to AGENTS.md §4.2, point 1.** It says the donor basis is `12,000 x 7,760`. No
operator on disk has that shape. The *first* whitener (`data/donors/whitener.npz`, the donor
INTERSECTION axis) was `18,000 × 7,312`; the **shipped** one is `29,394 × 17,942` on the donor
UNION axis minus the panel. The `12,000` is the per-donor row cap, not the basis size, and
`7,760` matches nothing. §4.2's mechanism is unaffected — the panel-basis-vs-donor-basis
argument stands — but the numbers should read `29,394 x 17,942` for production.

**Second correction, smaller but load-bearing for this mandate.** The mandate (and §4.2) describe
production as `P (C + lam I)^(-gamma/2)` with "`lam` a scalar ridge". **The shipped configuration
has `lam = 0.0`.** a15 whitens with `C^(-gamma/2)` and no ridge whatsoever; the only
regularisation present is `Whitener.fit`'s relative eigenvalue floor at `ev.max() * 1e-12`. The
one arm that ever shipped a non-zero ridge is a19 (`lam = 0.007`), which the board scored at
**−0.01352** and which was reverted. So "the regulariser is a ridge" was never quite true, and
the question this experiment had to answer is sharper than the mandate framed it: *should there
be a regulariser at all?*

---

## 1. Phase 1a — the spectrum. **There is no Marchenko–Pastur bulk to clip.**

Recomputed from each frozen operator's own eigenvalues (`work/spectrum-inventory.log`).

| operator | n | p | q = p/n | ev max | ev min | median | condition |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `union_pruned` (**ships**) | 29,394 | 17,942 | 0.6104 | 5.8946 | 5.895e−12 | 1.228e−03 | **1.000e+12** |
| `union` | 33,787 | 17,942 | 0.5310 | 6.2567 | 6.257e−12 | 1.378e−03 | **1.000e+12** |
| `union_full` | 45,531 | 17,942 | 0.3941 | 6.9802 | 6.980e−12 | 1.639e−03 | **1.000e+12** |
| `held_k562_gwps` | 35,744 | 17,806 | 0.4982 | 8.4259 | 8.426e−12 | 1.421e−03 | **1.000e+12** |
| `held_orion_hct116` | 27,764 | 17,942 | 0.6462 | 8.0363 | 8.036e−12 | 1.060e−03 | **1.000e+12** |
| `held_orion_hek293t` | 27,554 | 17,942 | 0.6512 | 3.6958 | 3.696e−12 | 1.201e−03 | **1.000e+12** |

**The condition number is exactly `1.000e+12` on all six**, which is `Whitener.fit`'s floor
`ev.max() * 1e-12`, not a property of the data. The floor is *binding*: 48 eigenvalues of the
shipped operator sit exactly on it, 144 are below `1e-9` relative and 2,247 (12.5 % of the
dimensions) below `1e-6` relative.

**The MP fit fails, and that is the result, not a failure of the fit.** Iterating
`sigma^2 = mean(ev below the edge)`, `edge = sigma^2 (1 + sqrt q)^2` collapses onto the floor on
every operator: it ends with 97–99 % of eigenvalues declared "spikes" and **0.00 % of the trace
inside the bulk**. There is no `sigma^2 (1 + sqrt q)^2` edge because the empirical spectrum is
not a spike-plus-bulk at all — it is a smooth decay over **twelve decades**. Quantiles of the
shipped operator:

```
   min     1%      5%      10%      25%      50%      75%      90%      99%      max
5.9e-12 9.9e-09 4.2e-07  2.9e-06  1.0e-04  1.2e-03  4.0e-03  7.9e-03  1.7e-02  5.89
```

Consequence: **eigenvalue clipping / MP denoising has nothing to clip.** `clipMP` in §2 is
numerically indistinguishable from the identity (+0.000002 scaled), because the fitted edge lands
at `3.0e-08`, below 95 % of the spectrum.

### 1b. What the whitener's gain actually multiplies

`f = ev^(-0.325)` is decreasing, so the operator's squared gain is concentrated where the
eigenvalues are smallest. **99.67 %** of `sum f^2` sits in the 12.5 % of directions below `1e-6`
relative, and the floored directions are amplified **505×** more than the median one. Read alone
that looks alarming: the operator is pouring almost all of its gain into a numerical null
space. §3 measures whether any
energy is actually there. (It is not, which is why the ridge has nothing to fix.)

---

## 2. Phase 1b — the `f(ev)` screen on PLOLO. **Nothing wins.**

Instrument: 162's PLOLO, a15 config (`alpha 0.982`, `cd4_f 0.10`, no mass fix), four donor
folds over the real 300-target panel, rendered into challenge frames A/B/C and scored against
the held library's own pseudobulk. Each fold uses its own leakage-safe operator; the arm changes
**only** `f(ev)`. `Operator.apply` restores the input's median row norm, so every arm here is
shape-only by construction — no amplitude ride-along.

PLOLO's own leave-one-arm-out sigma is `0.001494` scaled = **`0.000675` raw**; that is the
`x sigma` column. Full table: `work/phase1-combined.txt`. The head of it:

| map | hct116 | hek293t | k562 | kolf2 | mean Δ | scaled | +folds | ×σ |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: | ---: |
| `head^1.2` *(the γ axis — see below)* | +0.002211 | +0.004527 | +0.000328 | −0.000341 | +0.001681 | +0.003720 | 3/4 | **+2.49** |
| `cliprank12000` | +0.000262 | −0.000720 | +0.000703 | +0.001769 | +0.000504 | +0.001114 | 3/4 | +0.75 |
| `tail^1.2` | +0.000314 | −0.000341 | +0.000134 | +0.000507 | +0.000154 | +0.000340 | 3/4 | +0.23 |
| `ridge3e-06` | +0.000036 | +0.000066 | +0.000037 | −0.000009 | +0.000032 | +0.000072 | 3/4 | +0.05 |
| `clipMP` | +0.000008 | +0.000000 | +0.000005 | −0.000009 | +0.000001 | +0.000002 | 2/4 | +0.00 |
| **`prod(lam=0)` — production** | 0.787706 | 0.766498 | 0.795784 | 0.611201 | — | — | — | — |
| `ridge1e-05` | −0.000111 | +0.000058 | −0.000023 | −0.000184 | −0.000065 | −0.000144 | 1/4 | −0.10 |
| `tail^0.8` | −0.000167 | +0.000252 | −0.000305 | −0.000507 | −0.000182 | −0.000403 | 1/4 | −0.27 |
| `ridge0.0001` | −0.000449 | +0.000089 | −0.000203 | −0.000530 | −0.000273 | −0.000605 | 1/4 | −0.41 |
| `rie-LP-eta0.5` (RIE) | −0.000537 | +0.000015 | −0.002446 | −0.001647 | −0.001154 | −0.002552 | 1/4 | −1.71 |
| `cliprank2048` | +0.003536 | −0.000035 | −0.003042 | −0.005319 | −0.001215 | −0.002689 | 1/4 | −1.80 |
| `ridge0.007` *(a19's ridge)* | −0.006383 | −0.005653 | −0.004060 | −0.000754 | −0.004213 | −0.009320 | 0/4 | −6.24 |
| `ridge0.06` | −0.022901 | −0.019472 | −0.009872 | +0.000341 | −0.012976 | −0.028708 | 1/4 | −19.2 |
| `nls-LW2020` (analytical NLS) | −0.229082 | −0.221858 | −0.246384 | −0.107026 | −0.201088 | **−0.444884** | 0/4 | −298 |

Four things this table settles.

**(a) `lam = 0` is the ridge family's optimum.** Thirteen `lam` values from `1e-9` to `0.06`:
the largest reading anywhere on the axis is **+0.05 σ** (`lam = 3e-6`, +0.000072 scaled), the
axis is flat to ±0.05 σ below `1e-5`, and the loss grows **monotonically** above it —
−0.000144, −0.000359, −0.000605, −0.001128, −0.002807, −0.005686, −0.009320, −0.016790,
−0.028708 scaled at `lam` = 1e−5, 3e−5, 1e−4, 3e−4, 1e−3, 3e−3, 7e−3, 2e−2, 6e−2. `ridge0.007` reproduces a19's ridge and reads `−0.009320` scaled — the
right sign and the right order of magnitude for a19's board outcome, which is a free retrodiction
of an arm this screen never saw.

**(b) Ledoit–Wolf LINEAR shrinkage is inside that family, algebraically, and is priced at its own
optimum.** `d = (1-rho) ev + rho mu` equals `(1-rho)(ev + lam)` with `lam = rho mu / (1-rho)`, and
the median-row-norm renormalisation makes the leading constant invisible to every metric — so
*a lam sweep IS a rho sweep*. The optimal intensity is computable from the frozen spectrum
(RBLW and OAS, Chen et al. 2010, depend only on `tr C` and `tr C^2`):

| operator | mean ev | rho (RBLW = OAS to 5 dp) | equivalent ridge `lam` | PLOLO at that `lam` |
| --- | ---: | ---: | ---: | ---: |
| `union_pruned` (ships) | 3.31e−03 | **0.00338** | **1.12e−05** | −0.000144 scaled |
| `held_orion_hct116` | 3.40e−03 | 0.00201 | 6.83e−06 | ≈ 0 |
| `held_orion_hek293t` | 3.79e−03 | 0.00890 | 3.40e−05 | −0.000359 scaled |
| `held_k562_gwps` | 3.51e−03 | 0.00152 | 5.35e−06 | ≈ 0 |

**The statistically optimal linear shrinkage for this covariance is `rho = 0.0034`, and at that
intensity PLOLO reads −0.000144 scaled: zero, on the wrong side.** The whole linear family is
closed by measurement, not by argument.

**(c) The nonlinear estimators lose, and not because the port is broken.** `nls-LW2020` is a
direct implementation of Ledoit & Wolf (2020) and it destroys the arm (−0.44 scaled; on the
hct116 fold `pds` falls from 0.7877 to 0.5586, near chance). Before that is allowed to count as a result, `nlscheck193.py`
tests both estimators on synthetic data with a KNOWN population covariance, scoring PRIAL against
the oracle rotationally-invariant estimator `d_i = u_i' Sigma u_i`:

| population spectrum | span | analytical NLS | RIE (Ledoit–Péché) | MP clip | best ridge |
| --- | ---: | ---: | ---: | ---: | ---: |
| spike-plus-bulk, p=500 n=1000 | 2.7e+02 | **+97.24 %** | +96.53 % | +81.00 % | −1.83 % |
| spike-plus-bulk, p=1000 n=1600 | 5.3e+02 | **+98.66 %** | +98.33 % | +80.46 % | −1.43 % |
| power law over 6 decades, p=500 n=1000 | 3.0e+06 | **−5,972 %** | +82.12 % | −0.00 % | −1,788 % |
| power law over 6 decades, p=1000 n=1600 | 4.7e+06 | **−30,457 %** | +92.57 % | −0.00 % | −1,527 % |

**The port is correct** (97–99 % PRIAL on the spectrum shape the paper assumes, matching its
published range) **and it is mis-specified on the spectrum we actually have.** Ledoit–Wolf's
kernel bandwidth is proportional to the eigenvalue, `h * lambda_j` with `h = n^(-1/3)`; over
twelve decades that kernel is meaningless in the tail. This is a mechanism, not a shrug: it says
the analytical estimator is the *wrong tool for this matrix*, and it predicted the PLOLO collapse
before PLOLO was consulted.

**The RIE survives the specification test and still buys nothing.** `rie-LP-eta0.5` scores
**+82 % to +93 % PRIAL on the power-law spectrum** — it is a genuinely good estimator of this
covariance — and on PLOLO it reads **−0.002552 scaled, 1.7 σ negative.** That is the sharpest
statement in this experiment: *the estimator that is near-optimal in Frobenius loss on this exact
spectrum shape makes retrieval worse.* It is AGENTS.md §4.2 point 4 in a new form —
discriminability is not signal-to-noise, and it is not covariance-estimation accuracy either.

**(d) The optimum is INTERIOR, not a boundary.** Every regulariser in the literature lives on the
"compress the spectrum" side, so a sweep of them alone could only ever show that production sits
at an edge. `piecewise_power` bends the spectrum about its own median with a separate exponent
above and below, `theta = 1` being production exactly:

| direction | `theta` | scaled Δ`pds` | |
| --- | ---: | ---: | --- |
| tail compressed (all shrinkage estimators live here) | 0.80 / 0.90 / 0.95 | −0.000403 / −0.000224 / −0.000135 | worse |
| **production** | **1.00** | **0** | |
| tail stretched (anti-regularisation) | 1.05 / 1.10 / 1.20 | +0.000129 / +0.000100 / +0.000340 | ≤ 0.23 σ |

The tail axis is **flat to ±0.3 σ across a two-sided scan** and production sits in the middle of
it. There is no gradient to exploit in either direction.

### The one map with a gradient is the `gamma` axis in disguise

`head^theta` reads +2.49 σ at `theta = 1.2` and rises monotonically. It is not a new lever:
`m (ev/m)^theta` raised to `-gamma/2` is `ev^(-theta*gamma/2)` times a constant, so **`head^1.2`
is `gamma = 0.78` applied to the top half of the spectrum**, with `gamma = 0.65` kept below. The
mandate forbids re-sweeping `gamma`, 170 established it is a local optimum whose vertex is on the
shipped 0.65 (`mse` 0.8813 / **0.8738** / 0.8795 at 0.60 / 0.65 / 0.70) and that *every* gamma in
[0.35, 0.80] is net negative on all six members, and a19 shipped a partial version of this move
and lost −0.01352 on the board. §5 prices `head^1.2` anyway and it comes out **net negative**
(−0.002143) exactly as 170 predicts, on `nmae`. It is reported so the next agent knows the
gradient is real and already spoken for, not so anyone builds it.

---

## 3. Why every regulariser loses — the mechanism, per spectral band

`work/energy.log` decomposes the predicted effect matrix in the operator's own eigenbasis:
energy share before and after whitening, the mean between-row cosine inside the band (what
creates hubness), and the alignment with the held library's truth. Shipped operator, `kolf2`:

| band `ev/ev_max` | dims | E before | E after | mean row-row cos | truth alignment |
| --- | ---: | ---: | ---: | ---: | ---: |
| < 1e−9 | 144 | 0.00000 | 0.00018 | 0.0030 | 0.0022 |
| 1e−9 … 1e−6 | 2,103 | 0.00016 | 0.02496 | 0.0013 | 0.0007 |
| 1e−6 … 1e−4 | 4,671 | 0.03121 | 0.21390 | 0.0018 | 0.0026 |
| 1e−4 … 1e−3 | 8,157 | 0.35473 | 0.51873 | 0.0017 | 0.0049 |
| 1e−3 … 1e−2 | 2,845 | 0.43570 | 0.23596 | 0.0012 | 0.0158 |
| **1e−2 … 1** | **22** | **0.17819** | **0.00627** | **0.0300 / −0.0036** | **0.0341 / 0.0113** |

Two facts do all the work.

1. **The near-null tail holds no energy to regularise.** Below `1e-6` relative — the 12.5 % of
   directions carrying 99.67 % of the operator's squared gain — the prediction has
   **0.016 %** of its pre-whitening energy. A 505× relative amplification of ~zero is ~zero: it comes
   out at 2.5 % of the post-whitening energy. That is why `lam` is flat from `0` to `1e-5` and
   why the floor at `1e-12` is doing its job: *the ridge would be correcting a quantity that is
   not there.* (Mechanically: the covariance's null directions are directions in which no donor
   effect vector has any component, and the prediction is built from donor effect vectors.)
2. **Any regulariser strong enough to be felt reaches into the band that carries the answer.**
   75 % of the post-whitening energy is in `1e-4 … 1e-2`. A ridge only bites there once
   `lam ~ 1e-4`, and that is exactly where §2's table turns monotonically negative.

The table also quantifies §4.2's hubness story for the first time: **22 directions out of 17,942
carry 17.8 % of the raw energy at 20× the between-row cosine of everything else, and whitening
cuts them to 0.63 %.** Those same 22 directions carry the highest alignment with the truth — so
whitening is *deliberately* discarding the best-transferred part of the signal in order to remove
the part that makes every row look like every other row. That is why capping the head
(`wins0.25-0.9`) is the worst map in the whole screen at −0.116 scaled.

---

## 4. Phase 2 — the basis. **More rows does not help. Composition does.**

Leakage rule: all three union-axis operators contain rows from all three production donors, so
only the `kolf2` fold — a non-production donor that contributed no row to any of them — can
express an operator swap. That is 162's own constraint on its `a2-prune` arm, and gate 3 above
reproduces 162's number through this code path before anything new is measured.

`f(ev)` held at production's, a15 config, `kolf2`, challenge frame:

| operator | basis rows | sha | `pds` | Δ vs production | scaled |
| --- | ---: | --- | ---: | ---: | ---: |
| **`union_pruned` (ships)** | **29,394** | `e2f65ecf` | **0.611201** | — | — |
| `union` | 33,787 | `5608e807` | 0.610249 | **−0.000952** | −0.002105 |
| `union_full` | 45,531 | `965d3cef` | 0.609944 | **−0.001257** | −0.002781 |

Decomposed against the mandate's two questions:

* **"Orion carries 17,767 + 17,977 targets and the operator only uses 12,000 of each."** True,
  and using all of them is what `union_full` is. Holding composition fixed, lifting the cap adds
  **11,744 rows (+35 %) and costs −0.000305 raw `pds`** (`union → union_full`). That is 0.45 σ of
  PLOLO's noise — call it flat, certainly not a gain. **`C` is not fit on too little data.**
* **Composition still pays.** `union → union_pruned` removes 4,393 common-essential rows and is
  worth **+0.000952 raw** at a15's config (+0.001360 at v11's, gate 3). This is a clean
  independent confirmation of 091's `+0.002409` result on a completely different instrument
  (091: count-level HCT116 harness, `avg6`; here: PLOLO `pds` on the leak-free fold).

**So the operator's basis is size-saturated and composition-limited.** The lever that is left in
the basis is *which* perturbations, not *how many* — and that is 194's territory, not the
operator's. Recorded here so 194 does not re-derive it: at 29,394 rows the size axis is flat, and
the only measured way to move the basis is to remove rows that do not belong.

---

## 5. Phase 3 — the six-member price (rule 12), on the instrument each member is validated on

Rule 12 exists because of a19, which is *this exact family*: a whitener ridge changed the effect
vector's shape, `mse` read the profile directly, and a two-surface screen that never computed it
lost −0.01352. So no map here is disposed on `pds` alone.

**Instruments.** `pds` — PLOLO, the only board-validated member (6/6 signs, LOO RMS 0.001494
scaled), scaled/6 × 1.011. `nmae` — h1-test, 188 §2.4's exact identity
`nmae = 1 + (U − 2M)/D` on the **reference's own BH-significant pool**, measured median **89
genes/row**, which is the significance-selected pool 188 requires (a top-N-by-magnitude pool
inverts the sign). `reach` — h1-test × 157's `0.3930` board transfer factor. `fid` — h1-test
precision with the call set **held** at the incumbent's (rule 9: `n_real` is unreadable here).
`jac` — pinned to 0 by the held call set; stated, not claimed. `mse` — §6.

The h1-test arm is expressed by the operator-ratio identity `W_alt W_prod^-1 = V f_alt/f_prod V'`
applied to the recovered a15 effect, with a15's own operator (`e2f65ecf`, asserted). Two
approximations are carried and both inflate the arm if anything: `H` is recovered post-render,
and it contains the post-whitening CD4 residual (~10 % of each row) that the transform also
touches. A null result under them is the conservative reading.

| map | PLOLO Δ`pds` scaled | `pds` ovr | `nmae` ovr | `reach` ovr | `fid` ovr | **SUM(5)** | gate +0.0019 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| **`cliprank12000`** | +0.001114 | +0.000188 | +0.000178 | +0.000044 | +0.000100 | **+0.000510** | **fail (3.7×)** |
| `ridge1e-05` | −0.000144 | −0.000024 | −0.000011 | +0.000001 | +0.000047 | +0.000013 | fail |
| `clipMP` | +0.000002 | +0.000000 | −0.000000 | +0.000000 | +0.000000 | +0.000000 | fail |
| `tail^1.2` | +0.000340 | +0.000057 | +0.000045 | +0.000056 | −0.000159 | −0.000001 | fail |
| `floorrel1e-06` | −0.000019 | −0.000003 | +0.000000 | +0.000000 | −0.000003 | −0.000006 | fail |
| **`head^1.2`** *(the γ axis)* | **+0.003720** | **+0.000627** | **−0.002498** | −0.000017 | −0.000255 | **−0.002143** | fail |
| `tail^0.8` | −0.000403 | −0.000068 | −0.005447 | +0.000159 | +0.000121 | −0.005235 | fail |
| `ridge0.0001` | −0.000605 | −0.000102 | −0.005525 | +0.000058 | +0.000166 | −0.005403 | fail |
| `rie-LP-eta0.5` | −0.002552 | −0.000430 | −0.045883 | −0.000167 | +0.000625 | −0.045855 | fail |
| `ridge0.007` *(a19's ridge)* | −0.009320 | −0.001570 | −0.067952 | −0.000258 | +0.001043 | −0.068737 | fail |
| `nls-LW2020` | −0.444884 | −0.074963 | −0.077019 | −0.008276 | −0.006196 | −0.166454 | fail |

**Every map fails, and the best one fails by 3.7×.** Two readings worth keeping:

**`nmae` independently forbids regularisation, and the reason is exact.** On the reference's DE
pool the sign accuracy is `q = 0.5375` and the magnitude ratio is `k = 0.278`, so 183's landscape
identity `nmae = q|k−1| + (1−q)(k+1)` reduces at `k < 1` to `nmae = 1 + k(1 − 2q)` with
`1 − 2q = −0.075`. **`nmae` therefore improves with LARGER `k`** — bigger predictions on the
pool — and every regulariser shrinks `k` (0.2783 → 0.2774 → 0.2748 → 0.2548 as `lam` goes
0 → 1e−5 → 1e−4 → 0.007). This is a second, mechanistically independent reason `lam = 0` is
right, and it is the term that turns `head^1.2`'s +2.49 σ `pds` gain into a net loss.

**a19 is retrodicted twice by this screen without being fitted to.** `ridge0.007` reads −0.0093
scaled on `pds` and −0.068 overall on `nmae` — an arm this experiment reproduced from the map
family alone, whose board result (−0.01352 overall) was recorded a week earlier.

---

## 6. `mse` — 096's `claim` route on the challenge panel

`claim188.py`'s four known-answer gates were run verbatim first and all four pass:
`0.11860` / `0.19438` vs `0.19216` (+1.16 %) / `0.16652` vs `0.16478` (+1.06 %) / `0.24581`,
inside the instrument's own 3 % tolerance. The gated identity reads raw `mse` **0.87997**
against 164's ungated `0.88228`, and claim188 records that gating is worth exactly `-0.0025`
raw — so the baseline lands where it should.

**Rule 11 note on the object.** 096's cached object is v11's delivered effect and v11 whitened
with `data/donors/whitener_union.npz` (sha `5608e807`), NOT with a15's `union_pruned`. The
operator-ratio identity is only meaningful against the operator that did the whitening, so this
table is computed with `union` and is an extrapolation to a15's operator, which differs from it
by 091's prune alone.

| map | A | rA | B | claim | raw `mse` | Δ raw | score | **Δ overall** | branch |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| **`prod(lam=0)`** | 0.28280 | 1.00000 | 0.28106 | 0.24623 | 0.87997 | — | 0.11355 | — | `s^2 B` |
| `ridge1e-06` | 0.28283 | 1.00008 | 0.28109 | 0.24625 | 0.87996 | −0.00002 | 0.11357 | +0.000003 | `s^2 B` |
| `ridge1e-05` | 0.28276 | 0.99984 | 0.28102 | 0.24624 | 0.87995 | −0.00003 | 0.11358 | +0.000005 | `s^2 B` |
| `ridge0.0001` | 0.28237 | 0.99848 | 0.28064 | 0.24620 | 0.87994 | −0.00004 | 0.11359 | +0.000007 | `s^2 B` |
| `ridge0.007` *(a19's)* | 0.28383 | 1.00364 | 0.28201 | 0.24617 | 0.88026 | +0.00029 | 0.11325 | **−0.000050** | `s^2 B` |
| `floorrel1e-06` | 0.28284 | 1.00012 | 0.28110 | 0.24624 | 0.87995 | −0.00003 | 0.11357 | +0.000004 | `s^2 B` |
| `cliprank12000` | 0.28236 | 0.99843 | 0.28062 | 0.24615 | 0.87996 | −0.00001 | 0.11356 | +0.000002 | `s^2 B` |
| `clipMP` | 0.28280 | 1.00000 | 0.28106 | 0.24623 | 0.87997 | +0.00000 | 0.11355 | −0.000000 | `s^2 B` |
| `rie-LP-eta0.5` | 0.28497 | 1.00766 | 0.28324 | 0.24621 | 0.87992 | −0.00005 | 0.11360 | +0.000009 | `s^2 B` |
| `tail^0.8` | 0.28231 | 0.99826 | 0.28058 | 0.24621 | 0.87994 | −0.00004 | 0.11359 | +0.000007 | `s^2 B` |
| `tail^1.2` | 0.28310 | 1.00103 | 0.28134 | 0.24619 | 0.88004 | +0.00006 | 0.11348 | −0.000011 | `s^2 B` |
| `head^1.2` | 0.28096 | 0.99347 | 0.27928 | 0.24620 | 0.87975 | −0.00022 | 0.11378 | +0.000039 | `s^2 B` |

**`mse` is essentially blind to the eigenvalue map at fixed amplitude: every arm is inside
`|Δ overall| = 5e-05`, and every arm stays on the `s^2 B` branch — none of them goes near the
`claim` cliff.** The mechanism is that `Operator.apply` restores the median row norm, so the
map moves the effect's *direction* and barely its mean square (`rA` within 0.8 % on every arm,
0.36 % on a19's ridge).

**A flag for whoever revisits a19's decomposition.** 170 attributes **79 %** of a19's
`-0.01352` board `mse` loss to the ridge `lam` and 21 % to `gamma`. Measured directly here, the
ridge at fixed `gamma` and fixed median row norm is worth **−0.000050 overall** — 270× smaller
than that attribution implies. Three readings are possible and this experiment cannot
discriminate: (i) 096's object is v11's *delivered* effect, so a fresh build's amplitude
response is partly frozen into it; (ii) a19's loss is dominated by the `gamma 0.65 -> 0.80` move
crossing the `rho` cliff, which only an end-to-end rebuild expresses; (iii) 170's 79/21 split is
a model-based attribution rather than a direct measurement of the `lam` term. **It does not
change this verdict** — every map here is already negative or flat on `pds` and on `nmae`
independently — but the "79 % is the ridge" line should not be quoted as measured until someone
rebuilds a19 with `lam` as the only variable.

### The composed six-member table

| map | `pds` | `mse` | `nmae` | `reach` | `fid` | `jac` | **SUM** | gate +0.0019 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| **`cliprank12000`** | +0.000188 | +0.000002 | +0.000178 | +0.000044 | +0.000100 | 0 (held) | **+0.000512** | **fail, 3.7x** |
| `ridge1e-05` | −0.000024 | +0.000005 | −0.000011 | +0.000001 | +0.000047 | 0 | +0.000018 | fail |
| `clipMP` | +0.000000 | −0.000000 | −0.000000 | +0.000000 | +0.000000 | 0 | +0.000000 | fail |
| `ridge1e-06` | −0.000001 | +0.000003 | −0.000001 | +0.000000 | −0.000001 | 0 | +0.000000 | fail |
| `floorrel1e-06` | −0.000003 | +0.000004 | +0.000000 | +0.000000 | −0.000003 | 0 | −0.000002 | fail |
| `tail^1.2` | +0.000057 | −0.000011 | +0.000045 | +0.000056 | −0.000159 | 0 | −0.000012 | fail |
| **`head^1.2`** *(the γ axis)* | **+0.000627** | +0.000039 | **−0.002498** | −0.000017 | −0.000255 | 0 | **−0.002104** | fail |
| `tail^0.8` | −0.000068 | +0.000007 | −0.005447 | +0.000159 | +0.000121 | 0 | −0.005228 | fail |
| `ridge0.0001` | −0.000102 | +0.000007 | −0.005525 | +0.000058 | +0.000166 | 0 | −0.005397 | fail |
| `rie-LP-eta0.5` | −0.000430 | +0.000009 | −0.045883 | −0.000167 | +0.000625 | 0 | −0.045846 | fail |
| `ridge0.007` *(a19's)* | −0.001570 | −0.000050 | −0.067952 | −0.000258 | +0.001043 | 0 | −0.068787 | fail |
| `nls-LW2020` | −0.074963 | not run | −0.077019 | −0.008276 | −0.006196 | 0 | −0.166454 | fail |

**Nothing reaches a third of the floor. The single best map in fifty-five is `+0.000512`
against `+0.0019`.**

---

## 7. The structured (non-rotationally-invariant) target

`clip_rank` is a factor-model target and `ridge` is the identity target, but both keep
production's eigenvectors, so neither tests the mandate's fourth bullet. A **diagonal** target
does not: `C_rho = (1 - rho) C + rho diag(C)` has different eigenvectors and needs a fresh
eigendecomposition of the full 17,942 x 17,942 matrix. `diag(C)` is the per-gene variance of the
donor effect vectors — the "genes vary by different amounts but are uncorrelated" model — so
this arm asks directly how much of the operator's value is in the CORRELATIONS rather than in
the per-gene scales, which is the same structure that causes the hubness whitening exists to fix.

The `rho` ladder is run on `kolf2`, the only fold where the PRODUCTION operator is leakage-free;
the replication below uses `orion_hct116` with its own held operator. `rho = 0` is gate 4 on both:
reconstruct `C` from the frozen operator, re-solve its eigensystem from scratch, re-whiten, and
land back on the frozen operator's own row.

| `rho` | max eigenvalue | `pds` (challenge frame) | Δ vs `rho = 0` | scaled |
| ---: | ---: | ---: | ---: | ---: |
| **0 (gate 4)** | 5.8946 | **0.611201** — reproduces the frozen operator exactly | — | — |
| 0.05 | 5.6002 | 0.610783 | −0.000417 | −0.000924 |
| 0.25 | 4.4226 | 0.610483 | −0.000718 | −0.001589 |
| 0.75 | 1.4786 | 0.612040 | +0.000839 | +0.001857 |

On `kolf2` the response is **not monotone**: light shrinkage loses (`rho = 0.05`, `0.25`) but
heavy shrinkage — where the operator is most of the way to a per-gene variance normaliser — comes
back `+0.001857` scaled, `+1.24 σ` of PLOLO's leave-one-arm-out noise. That is one fold with no
fold-scatter error bar, against a per-fold scatter of 0.0002–0.006 raw for the maps in §2, i.e.
exactly the size of thing a single fold cannot resolve. **So it was replicated on a second fold
before being written down.**

| fold | operator | sha | `rho = 0` (gate 4) | `rho = 0.75` | Δ | scaled |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| `kolf2` | `whitener_union_pruned` | `e2f65ecf` | 0.611201 | 0.612040 | **+0.000839** | +0.001857 |
| `orion_hct116` | `held_orion_hct116` | `7cecbf7f` | 0.787706 | 0.777063 | **−0.010643** | **−0.023547** |
| | | | | **2-fold mean** | **−0.004902** | **−0.010845** |

**It does not replicate — it inverts, and the second fold's loss is 13× the first fold's gain.**
Gate 4 passes at `|d| = 0.00e+00` on both folds, so the two rows are measured through an
identical, verified code path. The `kolf2` `+1.24 σ` was fold noise, exactly as its own error bar
warned. The diagonal-target axis is negative on the mean and by a wide margin, and the honest
summary of it is: *light shrinkage toward the diagonal loses on the one fold that can see it,
and heavy shrinkage loses badly.*

That closes the last family in the mandate's list. Rotationally-invariant regularisers lose
(§2); the one non-rotationally-invariant structured target loses too, at every intensity that
survives replication. The correlation structure the diagonal target throws away is not noise the
operator should be protected from — it is what the operator is made of.

---

## 8. Verdict

**NO-GO on every arm. The ridge is already optimal, and the geometry is closed.**

`a15` stays exactly as it is: `WHITEN_GAMMA=0.65`, no ridge, operator
`data/donors/whitener_union_pruned.npz` sha
`e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`.

What is now closed, with numbers rather than argument:

1. **The ridge / linear-shrinkage family.** Thirteen `lam` values; flat to ±0.05 σ below
   `1e-5` and monotone decline above it; the statistically optimal Ledoit–Wolf intensity
   (`rho = 0.0034`, `lam = 1.1e-5`) measured at −0.000144 scaled. LW-linear is algebraically inside the ridge family under the operator's
   median-norm renormalisation, so this closes it entirely, not just at the points tested.
2. **The random-matrix family.** There is no MP bulk (0.00 % of the trace, spectrum spans twelve
   decades), so clipping is the identity. The analytical nonlinear estimator is *correctly
   implemented* (PRIAL +97 to +99 % on a spiked spectrum) and *mis-specified here* (−5,972 % on a
   power-law one); on PLOLO it costs −0.44 scaled.
3. **The RIE family — the strongest negative in the file.** Ledoit–Péché scores +82 to +93 %
   PRIAL on this spectrum's shape, i.e. it really is a better covariance estimate, and it reads
   −0.002552 scaled on `pds` and −0.0459 overall composed. **Improving the covariance estimate
   makes retrieval worse.**
4. **Both sides of the map-shape family.** The tail axis is flat to ±0.3 σ over a two-sided
   scan; production is interior, not at a boundary.
5. **The basis-size premise.** Lifting the 12,000/donor cap to every Orion row (+35 % rows) costs
   −0.000305 raw `pds`. Composition (091's prune) is worth +0.000952. Size is saturated.
6. **The structured (diagonal) target**, the only non-rotationally-invariant estimator in the
   mandate. Light shrinkage loses (−0.000924, −0.001589 scaled at `rho` 0.05, 0.25). Heavy
   shrinkage read `+1.24 σ` on one fold and **−0.023547 scaled on the second**, 2-fold mean
   −0.010845. It does not survive replication.

What is NOT closed, and belongs to someone else:

* **Basis composition.** The prune pays on two independent instruments. Which *other* rows do not
  belong is a donor-table question — **194's**, not the operator's.
* **`gamma`.** `head^theta` shows a real +2.49 σ `pds` gradient on the top half of the spectrum,
  and 170 + a19 + §5 all say it is paid for on `mse`/`nmae`. Nobody should reopen it without a
  magnitude map that survives `nmae`, which is 188's open item, not this one.

**Recommendation for AGENTS.md §4.2:** add a line saying the shipped operator is
`29,394 x 17,942` with `lam = 0`, that its regulariser is a `1e-12` relative eigenvalue floor,
and that 193 measured the ridge, linear-shrinkage, MP-clipping, RIE and nonlinear-shrinkage
families and found production optimal in all of them. This is the third independent closure of
the geometry (145 attenuation, 170 gamma, 193 regulariser) and it should stop the next agent from
spending a day on "the covariance estimator must be improvable".

---

## 9. Files

| file | what |
| --- | --- |
| `maps193.py` | the eigenvalue maps: ridge, LW/RBLW/OAS intensities, quantile and relative floors, winsorisation, rank clipping, MP clipping, Ledoit–Wolf 2020 analytical nonlinear shrinkage, Ledoit–Péché RIE, two-sided piecewise power |
| `screen193.py` | the PLOLO screen + known-answer gate 1 |
| `gate2_massfix.py` | known-answer gate 2 (a15 → a23 mass-fix price) |
| `nlscheck193.py` | the estimators' own known-answer row: PRIAL against the oracle RIE on synthetic spiked and power-law spectra |
| `phase2_basis.py` | the basis sweep + known-answer gate 3 |
| `energy193.py` | energy / hubness / truth-alignment by spectral band |
| `structured193.py` | the diagonal-target arm + known-answer gate 4 |
| `h1six193.py` | h1-test `pds` / `nmae` / `reach` / `fid` / `jac` + known-answer gate 5 |
| `mse193.py` | `mse` through 096's `claim` route + claim188's four gates |
| `work/phase1-combined.txt` | the full 55-map × 4-fold table |
| `work/*.json`, `work/*.log` | every number above, with the operator sha256 on each row |

## 10. What I did NOT measure, and why

* **A fresh eigendecomposition per basis size.** The basis sweep uses the three union-axis
  operators already on disk, so `union → union_full` is a clean cap-lift but
  `union_pruned → union` moves composition and size together. Gate 3 separates them at v11's
  config; a finer cap ladder (4,000 / 8,000 rows) was not built because the measured size
  gradient (−0.000305 over +35 % rows) is 0.45 σ and a further reduction cannot plausibly reach
  the +0.0014 raw needed to clear the gate.
* **`jac` as a real number.** The call set is held at the incumbent's, so `jac` is pinned to
  exactly 0 here. A map that changed the call set would need the emitter; none of these maps is
  close enough to the gate for that to change a decision.
* **Held-fold basis extensions.** 062's three held operators are already fitted on every donor
  row (no cap), so the "more data" question is only expressible on `kolf2`. Building capped
  variants of them would test the same gradient a third time.
* **`gamma`.** Out of mandate and closed by 170.
* **A fine `rho` ladder on the diagonal target, or all four folds at `rho = 0.75`.** Each point
  costs a full 17,942-dimension eigendecomposition (10–15 min wall on a contended machine). The
  budget was spent on cross-fold REPLICATION of the one non-negative reading rather than on
  resolution of the ladder, which is the right trade: the replication inverted the sign and a
  finer ladder would not have.

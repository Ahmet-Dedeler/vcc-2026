# 098 — `de_wilcoxon_lfc_nmae`, measured for the first time

**Headline: `nmae` is not a magnitude-calibration member. In our regime it is algebraically a
DIRECTION member, and its gap to the leader is the same direction gap `reach` and `fid` have.
Every calibration route is measured and closed; the only lever that moves it consistently on
both harnesses is the whitening exponent, which `experiment 097` priced at zero and which is
worth `+0.004 to +0.005 overall` — still not enough to pay its own `pds` bill on the worse
surface.**

Second headline, and it is the reusable one: **every `nmae` number this project has ever taken
was read on a reference gate 2–20x shallower than the board's, and the sign of the amplitude
response flips across that range.** At BH 0.05 both harnesses say "raise amplitude, `nmae`
improves". At challenge-like gate depth h1-test says the opposite — and the board, on its own
two amplitude moves, agrees with the deep gate.

---

## 0. What the metric actually is

From `docs/reference/vcc2026-metrics-FULL.md` §6, read exactly:

```
NMAE_p = sum_{g in S_p} |lfc_hat_{p,g} - lfc_{p,g}|  /  sum_{g in S_p} |lfc_{p,g}|
```

* `S_p` = the genes the **reference** calls significant (BH < 0.05), with finite reference lfc,
  `p`'s own target removed, `|S_p| >= 10`. **Chosen by the reference alone, identical for every
  submission.**
* A gate gene our table does not carry is read as a predicted lfc of **zero**, not masked.
* The denominator contains **only reference quantities**.

### (c) The normalisation is not a lever — this is settled, not estimated

`sum |lfc_ref|` is a functional of the withheld reference and of nothing we submit. The gate
membership, the gate-size cut and the denominator are all reference-side. **The only
submission-side quantity anywhere in this member is `lfc_hat`.** Route (c) is closed by reading,
with no experiment required.

### (d) The two things the spec makes obvious that nobody had used

**(d1) NMAE = 1 at zero prediction, exactly, and our prediction lives entirely in the linear
regime around that point.** Measured on the gate: `|lfc_hat|` is **0.138x (h1) / 0.168x (orion)**
of `|lfc_ref|` in the mean, and `|h| <= |r|` on **99.36% / 99.27%** of gate genes. On that side
`|h - r| = |r| - sign(r)*h` identically, so

```
NMAE  =  1  -  < sign(lfc_ref) , lfc_hat >_w ,        w_{p,g} = 1 / (P * sum_{S_p}|lfc_ref|)
```

Verified (`linear.py`): 0.97531 vs 0.97448 on h1 (residual `+0.00083`), 0.92747 vs 0.92387 on
orion (`+0.00359`).

**So the reference's fold-change MAGNITUDE enters only through the normaliser. What the metric
buys is signed predicted magnitude on genes the reference calls DE.** "Get the magnitudes right"
is the wrong objective for a submission at 0.14x amplitude; the right one is a directional inner
product. Everything below follows from that.

**(d2) The gate SIZE is the transfer variable, and it has never been controlled.** The challenge's
reference calls ~2,000 genes per row (`n_real`, back-solved from v14/v15's `fid` coverage:
1563/0.8950 = 1746 and 1911/0.9157 = 2087). Our harnesses call **89 (orion)** and **1,120 (h1)**
at BH 0.05. §4 shows the arm ranking inverts across that range.

---

## 1. Method, and what it is validated against

`nmae` is a deterministic functional of the predicted **profile**: the emitted cells' mean
per-cell CPM converges to `1e6 * prof`, so

```
lfc_hat[p,g] = log2( (1e6 * prof[p,g] + 1e-9) / (control mean CPM_g + 1e-9) )
prof[p]      = expm1(b_ctrl + eff[p]) / sum_g expm1(b_ctrl + eff[p])
```

`validate.py` drives production `emit_counts` (a12's tiered spread: top-8000 by `|effect|` at
0.25, tail 1.0, `phi` 4) on 40 rows per harness with two seeds and compares per row:

| harness | analytic | realised seed 0 | realised seed 1 | **twin gap** | realised − analytic | per-row corr |
|---|---:|---:|---:|---:|---:|---:|
| orion | 0.91764 | 0.91928 | 0.91930 | **0.00001** | +0.00165 | 0.99872 |
| h1-test | 0.97521 | 0.97551 | 0.97552 | **0.00001** | +0.00030 | 0.99919 |

Two consequences:

* The analytic evaluator is the metric, to `+0.0017 / +0.0003` in level and 0.999 per row. Every
  sweep below is therefore a real measurement, not a counterfactual.
* **The same-policy re-seeded twin gap for `nmae` is `1e-5`** — four orders below every number in
  this document. Unlike `reach` (twin gaps 0.021–0.028 raw, which retired that column entirely
  in 097), a per-arm twin on `nmae` is vacuous. The error bars here are **the two harnesses and
  the three gate depths**, and that is a stronger control than a seed. Twins are nonetheless
  carried through `gate6.py` and reported (they agree to `<= 0.0004`).

Full six-member emitted gate (`gate6.py`, 50 perturbations, twins) reproduces the analytic
`nmae` to `0.0005–0.0026` on every arm but one (`dewhiten CPM>200` on h1, 0.0054).

`pds` here is my own rank-cosine; it is **bit-identical** to `cell_eval2.metrics.discrimination`
(`0.865546` / `0.795339`, both harnesses, all digits).

---

## 2. Where the error is (measurement, `diag.py` / `diag2.py`)

| | h1-test | orion |
|---|---:|---:|
| scored rows / gate genes | 58 / 64,930 | 127 / 11,279 |
| mean `|lfc_ref|` on gate | 0.3639 | 0.4433 |
| mean `|lfc_hat|` on gate | 0.0502 | 0.0744 |
| **ratio** | **0.138** | **0.168** |
| sign agreement on gate | 0.5297 | 0.5801 |
| Pearson(h, r) | 0.1375 | 0.3624 |
| gate genes where we are WORSE than silence | **47.1%** | **42.0%** |
| gross credit from helped genes | +0.0837 | +0.1245 |
| damage from hurt genes | −0.0590 | −0.0520 |

**70% (h1) / 42% (orion) of our gross credit is destroyed by the genes we get wrong.** The whole
member is a signal-to-noise problem on per-gene direction.

Ceilings, for scale:

| what | h1 | orion |
|---|---:|---:|
| our prediction (a12) | 0.9753 | 0.9275 |
| best global scale on `lfc_hat` | 0.9625 | 0.9065 |
| best PER-ROW scale (oracle) | 0.9471 | 0.8832 |
| our sign at ORACLE magnitude | 0.8713 | 0.7801 |
| ORACLE sign at mean magnitude | **0.5686** | **0.4469** |

Board context: a12 raw `0.9448`, the leader `0.8923`, the half-depth replicate anchor `0.37–0.43`.
In credit terms (`1 − raw`) the leader has **1.95x** ours and the replicate **10.9x**. Closing the
`+0.0898` scaled gap needs raw `−0.0542`, i.e. **1.98x our current directional credit**.

---

## 3. Ranked routes, by measured Δ`nmae` on both harnesses

All shape arms are re-solved to a12's own `A` (what `--effect-scale` does), so shape and
amplitude are separate rows. Δ is raw `nmae`, **negative is better**. Twin gap `1e-5` throughout.

| # | route | buildable, no oracle? | Δ h1 @BH0.05 | Δ h1 @BH0.35 (|S|2600) | Δ orion @BH0.05 | Δ orion @BH0.60 (|S|2204) | verdict |
|---|---|---|---:|---:|---:|---:|---|
| **1** | **whitening `gamma` 0.65 → 0** | **yes** (config knob) | **−0.0109** | **−0.0149** | **−0.0185** | **−0.0124** | **the only consistent mover** |
| 1b | whitening `gamma` 0.65 → 0.45 | yes | −0.0095 | −0.0108 | −0.0076 | −0.0037 | same, gentler |
| 1c | de-whiten only CPM > 200 (1,173g) | yes | −0.0082 | −0.0105 | −0.0054 | −0.0015 | best `pds` ratio on h1 |
| 2 | global amplitude ×2 (no A hold) | yes | −0.0093 | **+0.0250** | −0.0211 | **+0.0112** | **sign flips with gate depth** |
| 3 | per-decile scale on `|lfc_hat|`, fitted on the OTHER harness | yes | — | — | — | — | **worse than a transferred flat scalar** (§5) |
| 4 | per-decile scale on control CPM, fitted on the OTHER harness | yes | — | — | — | — | same |
| 5 | sharpen `|lfc|^+0.25` at held A | yes | −0.0027 | +0.0005 | −0.0036 | −0.0028 | ≤0.004, sign-inconsistent |
| 6 | donor-agreement per-gene shrink (a13's key) at held A, `g`=0.25 | yes | −0.0023 | +0.0008 | −0.0069 | −0.0022 | sign-inconsistent; board-closed |
| 7 | donor-agreement AMPLIFY at held A | yes | +0.0085 | +0.0063 | +0.0309 | +0.0097 | **negative on both** |
| 8 | expression tilt `CPM^-0.5` at held A | yes | +0.0029 | +0.0047 | +0.0256 | +0.0079 | **negative on both** |
| 9 | per-perturbation reallocation (best buildable key) | yes | −0.0037 | −0.0013 | −0.0124 | −0.0006 | oracle ceiling only −0.016/−0.026 (§6) |
| 10 | remove the control-normalisation offset | yes | −0.0034 | −0.0039 | +0.0189 | +0.0224 | opposite signs, and costs `pds` 0.08–0.15 |
| 11 | constant directional bias on the gate | **structurally void** | ±0.0000 | ±0.0000 | ±0.0000 | ±0.0000 | §7 |

Only route 1 is negative on both harnesses at every gate depth, and monotone in the knob.

---

## 4. The gate-depth finding, and its board anchor

`gateladder.py` — same arms, gate relaxed by raising the BH threshold (a **proxy** for a more
powerful reference, labelled derived):

| harness | BH | mean `|S_p|` | raw `nmae` (a12) | best global lfc scale `c*` | gain at `c*` | weighted sign acc |
|---|---:|---:|---:|---:|---:|---:|
| h1 | 0.05 | 1,120 | 0.9791 | **1.90** | −0.0095 | 0.594 |
| h1 | 0.15 | 1,663 | 0.9840 | 1.42 | −0.0024 | 0.570 |
| h1 | 0.35 | 2,600 | 0.9914 | **0.90** | −0.0002 | 0.540 |
| h1 | 0.60 | 4,018 | 1.0040 | **0.48** | −0.0068 | 0.522 |
| orion | 0.05 | 90 | 0.9332 | **1.86** | −0.0197 | 0.686 |
| orion | 0.35 | 791 | 0.9678 | 1.78 | −0.0120 | 0.594 |
| orion | 0.60 | 2,204 | 0.9746 | 1.34 | −0.0029 | 0.568 |
| orion | 0.85 | 5,490 | 0.9948 | 0.64 | −0.0052 | 0.548 |

**The optimal amplitude falls from ~1.9 to below 1 as the gate deepens.** At the challenge's
~2,000 genes it is ~0.9–1.3, i.e. we are already there, and the whole magnitude-calibration
route is worth **`−0.002` to `−0.012` raw**, not the `−0.013 / −0.021` the BH-0.05 gate advertises.

**Board anchor (this is the part that makes it a result and not a curiosity).** The board has
measured two amplitude moves:

| pair | Δ ln amplitude | Δ`score_nmae` | implied d(raw nmae)/d(ln s) |
|---|---:|---:|---:|
| a10 → a11 (`s` 1.0 → 0.7762) | −0.2534 | +0.0026 | **+0.0062** |
| a11 → a12 (`s` 0.7762 → 0.9363) | +0.1874 | −0.0030 | **+0.0097** |

Harness readings of the same derivative:

| harness | BH 0.05 | BH 0.35 | BH 0.60 |
|---|---:|---:|---:|
| h1 | **−0.0179** | +0.0033 | +0.0289 |
| orion | **−0.0512** | −0.0266 | −0.0150 |

**At BH 0.05 both harnesses have the wrong sign. h1-test at BH 0.35–0.45 reproduces the board's
value. orion never does, at any depth.** So (i) every `nmae` reading in experiments 096 and 097
was taken on a cohort that ranks amplitude backwards, and (ii) orion's `nmae` derivative should
not be trusted at all — its 89-gene gate is 22x shallower than the board's.

Note this **reverses the premise this experiment was handed**: `mse` wants amplitude down and
`nmae` wants it down too, at challenge depth. They are aligned, not in tension. The exchange
rate settles it anyway: at a12's operating point one 1% amplitude cut costs `0.00144` scaled
`mse` and buys `0.00013` scaled `nmae` — **`mse` charges ~10x what `nmae` pays.** Amplitude is
closed for this member in both directions.

### The bracketing claim, checked on both arms as asked

| arm | h1-test | orion | **board** | bracketed? |
|---|---:|---:|---:|---|
| a12 | 0.9630 | 0.9189 | **0.9448** | **yes** |
| a13 | 0.9566 | 0.9538 | **1.0189** | **NO — the board is 0.065 outside both** |

**One arm of two. Do not lean on it.** Worse, on a13 the harnesses disagreed with the board in
*direction* on h1 and by 2.1x on orion. This document's deltas inherit that risk; nothing here
is safe at 8x.

And the gate ladder does **not** rescue a13 either: at challenge depth my reconstruction of its
key still reads `−0.0005 / −0.0010`. **Gate depth explains the amplitude sign flip; it does not
explain a13's board catastrophe, which remains unexplained.** Stated as an open defect.

---

## 5. Route (a) answered: the shape is a single global scalar, and per-gene calibration is a
## cell-line property

`diag2.py`, fitting an L1-optimal scale per stratum on one harness and applying it to the other:

| stratification | target | base | own-fit (in-sample) | **fitted on the other harness** | flat scalar from the other harness |
|---|---|---:|---:|---:|---:|
| `|lfc_hat|` decile | orion | 0.9275 | 0.8947 | **0.9129** | **0.9074** |
| `|lfc_hat|` decile | h1 | 0.9753 | 0.9603 | **0.9693** | **0.9672** |
| control-CPM decile | orion | 0.9275 | 0.8989 | **0.9107** | **0.9073** |
| control-CPM decile | h1 | 0.9753 | 0.9571 | 0.9618 | 0.9638 |

**In three of four cells the transferred ten-parameter calibration is WORSE than the transferred
one-parameter one.** By the rule in `RESEARCH-METHOD` §10 this is the `1.02x` shape — a per-gene
slope is a cell-line property. The flat scalar itself transfers essentially perfectly (h1's
`c* = 2.015` applied to orion gives 0.9074 against orion's own optimum 0.9065), which is what
makes the negative clean: the failure is the per-gene structure, not the fitting.

---

## 6. Route (b), and the one channel that is exactly `pds`-free

**Gate membership** is weakly predictable from control expression — median CPM in-gate vs
out-of-gate is **76.9 vs 45.3 (h1)** and **220.4 vs 39.8 (orion)** — i.e. the direction of the
competitor claim in `RESEARCH-METHOD` §10, on the *reference's* set rather than on ours. But it
is not usable: knowing a gene is more likely to be scored does not tell you which way it moves,
and a magnitude gain applied on an expression stratum is route 8, measured negative on both.

**The per-perturbation channel is exactly free on `pds`** (097 §9e: `pds_cosine` normalises each
row, so a per-row scale is bit-identical), which makes it the most valuable place to look. It is
also closed, with a ceiling (`perrow.py`, BH 0.35 / 0.60):

| harness | rows with `c_p < 0` | ORACLE per-row scale in [0,2] | best buildable key | its gain |
|---|---:|---:|---|---:|
| h1 @BH0.35 | 32.8% | −0.0163 | our mean `|lfc|` on gate | −0.0013 |
| h1 @BH0.60 | 56.9% | −0.0155 | (all keys positive = harmful) | +0.0009 |
| orion @BH0.35 | 1.6% | −0.0323 | donor cosine | −0.0015 |
| orion @BH0.60 | 2.4% | −0.0255 | donor cosine | −0.0008 |

`rho(donor cosine, c_p)` is **−0.166 / −0.061 / −0.074 on h1** and **+0.289 / +0.287 / +0.214 on
orion** — opposite signs, so the project's one transferring per-perturbation confidence score
does not transfer *for this member*. **Even a perfect per-row allocation is worth `+0.004` to
`+0.009` overall**, and the buildable part is `+0.0004`.

---

## 7. Route 11 is structurally void, which is worth recording

`orion`'s reference gate is 59.8% DOWN (unweighted) against h1's 52.5%, and applying a constant
`−0.05` to `lfc_hat` is worth **`−0.0238` on orion and `+0.0048` on h1**. That looks like a large
free lever and is not one, for two independent reasons:

1. It is a cell-line property with opposite signs on our two surfaces — the same trap that
   withdrew arm `a1-order-k2000`.
2. **It is not expressible.** Gated genes carry 99.4% of the profile's mass, so shifting all of
   them and renormalising cancels: the built arm measures `Δ = ±0.0000` at `delta = ±0.02` on
   both harnesses. A directional bias can only be bought on a *subset*, which requires predicting
   the subset, which is route (b).

---

## 8. The one live route, priced on all six members

`gate6.py` — production `emit_counts`, a12's tiered spread, `phi` 4, 50 perturbations, re-seeded
twins, amplitude held at a12's `A` so nothing here is an amplitude comparison.

| harness | arm | seed | n_pred | reach | fid | jac | nmae@0.05 | nmae@0.35 | nmae@0.60 | pds | kappa |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| h1 | a12 (g0.65) | 0 | 5854 | 0.0208 | 0.5099 | 0.0945 | 0.9650 | 1.0371 | 1.1077 | 0.8584 | 0.02074 |
| h1 | a12 (g0.65) | 1 | 5857 | 0.0174 | 0.5102 | 0.0947 | 0.9650 | 1.0371 | 1.1077 | | |
| h1 | **g0.45** | 0 | 5744 | 0.0265 | 0.5118 | 0.0951 | 0.9571 | 1.0272 | 1.0947 | **0.8611** | 0.02292 |
| h1 | **g0.45** | 1 | 5749 | 0.0209 | 0.5109 | 0.0952 | 0.9571 | 1.0272 | 1.0947 | | |
| h1 | g0.00 | 0 | 5508 | 0.0439 | 0.5147 | 0.0960 | 0.9567 | 1.0239 | 1.0864 | 0.8099 | 0.02944 |
| h1 | g0.00 | 1 | 5511 | 0.0433 | 0.5148 | 0.0960 | 0.9566 | 1.0238 | 1.0863 | | |
| h1 | de-whiten CPM>200 | 0 | 5815 | 0.0413 | 0.5118 | 0.0948 | 0.9635 | 1.0327 | 1.1021 | 0.8597 | 0.02252 |
| orion | a12 (g0.65) | 0 | 4326 | 0.0673 | 0.5578 | 0.0136 | 0.9324 | 0.9771 | 0.9901 | 0.7976 | 0.01737 |
| orion | a12 (g0.65) | 1 | 4324 | 0.0702 | 0.5584 | 0.0137 | 0.9320 | 0.9768 | 0.9899 | | |
| orion | **g0.45** | 0 | 4269 | 0.0857 | 0.5592 | 0.0138 | 0.9124 | 0.9672 | 0.9846 | 0.7779 | 0.01971 |
| orion | **g0.45** | 1 | 4271 | 0.0773 | 0.5599 | 0.0138 | 0.9121 | 0.9671 | 0.9845 | | |
| orion | g0.00 | 0 | 4192 | 0.0833 | 0.5651 | 0.0139 | 0.9002 | 0.9554 | 0.9726 | 0.6844 | 0.03385 |
| orion | g0.00 | 1 | 4198 | 0.0832 | 0.5657 | 0.0139 | 0.9000 | 0.9552 | 0.9724 | | |
| orion | de-whiten CPM>200 | 0 | 4332 | 0.0640 | 0.5580 | 0.0136 | 0.9240 | 0.9741 | 0.9879 | 0.7805 | 0.01917 |

`nmae` twins agree to `<= 0.0004`. `reach` twins differ by 0.003–0.016, straddling most between-arm
differences as usual — its column is reported and used for nothing but a sign.

### Cost of the best arm, member by member, in OVERALL points

`gamma` 0.65 → 0 (the full lever), amplitude held:

| member | h1 raw Δ | h1 overall | orion raw Δ | orion overall | how obtained |
|---|---:|---:|---:|---:|---|
| **`pds`** | **−0.0485** | **−0.0179** | **−0.1132** | **−0.0417** | measured, exact, `cell_eval2`-identical |
| `nmae` | −0.0132 | +0.0037 | −0.0175 | +0.0048 | measured at challenge-like depth |
| `mse` | (κ ×1.42) | +0.0115 | (κ ×1.95) | +0.0261 | **derived** — κ ratio onto the board's `kappa0` |
| `fid` (precision leg) | +0.0048 | +0.0027 | +0.0073 | +0.0041 | measured, emitted |
| `reach` | +0.0245 | +0.0046 | +0.0155 | +0.0029 | measured but **inside its twin floor on orion** |
| `jac` | +0.0015 | +0.0007 | +0.0003 | +0.0001 | sign only |
| **net** | | **+0.0052** | | **−0.0039** | |

**`pds` alone is 3.4x (h1) and 8.7x (orion) the `nmae` gain.** 097's closure of this route
survives — but its stated reason does not. 097 wrote *"`nmae` is NOT the binding member here — it
improves on orion and degrades on h1-test, disagreeing in sign, so it is read as zero."*
**Retracted here: `nmae` improves on BOTH surfaces at every gate depth, monotonically in `gamma`,
and is worth `+0.004 / +0.005` overall.** Two reasons 097 read it as zero, both convention:

* it measured at BH 0.05 only, the cohort §4 shows is the wrong one; and
* it re-solved every arm onto its own `rho` cliff, where harness amplitude collapses to
  **0.183x (h1) / 0.369x (orion)** of a12's `A` — because `claim/B` is 0.18/0.37 on the harnesses
  against **0.88 on the challenge**. `reconcile.py` shows that convention shrinks every `nmae`
  effect by ~4x on h1 and ~2x on orion. **The A-held convention is the one that resembles the
  board.**

---

## 9. Pre-registered prediction — the single best arm

**Arm `b1`: `WHITEN_GAMMA` 0.65 → 0.45, a12 otherwise unchanged, `--effect-scale` trimmed ~4% so
`A` is held.** Chosen over `gamma = 0` because the full lever is net-negative on orion by
`−0.0039` and `pds` is the member we lead board-wide.

Predicted, against a12's `0.16798741`:

| member | predicted raw | predicted Δ overall | basis |
|---|---:|---:|---|
| `pds` | 0.826 (Δ −0.008, range −0.020…+0.003) | **−0.0030** | h1 +0.0027 / orion −0.0197, board sits between |
| `nmae` | 0.937 (Δ −0.008) | **+0.0022** | h1 −0.0108 @BH0.35, orion −0.0056 @BH0.60 |
| `mse` | 0.883 (Δ −0.019) | **+0.0033** | κ ratio 1.10–1.14 onto `kappa0 = 0.0756`, derived |
| `fid` | +0.0017 | +0.0010 | measured precision leg, both surfaces |
| `reach` | — | 0 ± 0.005 | unreadable |
| `jac` | +0.0004 | +0.0002 | sign only |
| **overall** | **≈ 0.1717** | **+0.0037** (range −0.002 … +0.010) | |

**Falsifiers, in order of what they would retire:**

1. **raw `nmae` does not fall.** Retires the entire 098 reading: the gate-depth correction, the
   linearisation, and the claim that `nmae` responds to the whitener. Everything in this document
   downstream of §3 becomes void.
2. **raw `pds` falls by more than 0.020.** Then orion's `pds` sensitivity *understates* the board's
   and `WHITEN_GAMMA` must not be touched at any setting; 097's closure becomes unconditional.
3. **raw `mse` does not fall by at least 0.010.** Retires the κ-ratio transfer, which is the one
   *derived* number in the prediction and the largest positive term in it.
4. **overall below 0.16799.** The whitening route is closed for good with `nmae` correctly priced,
   and `nmae` has no remaining buildable lever at all.

**Recommendation, stated against my own arm: do not spend a slot on this.** Expected `+0.0037`
with a downside that lands on the one member we lead the board on, on a transfer path that was
wrong by 8x the last time it was used for `nmae` (a13). The measured value of 098 is the closure,
not the arm.

---

## 10. Measured / derived / assumed

**Measured** (executed against real data on both harnesses; twin gap `1e-5`):
analytic-vs-emitted agreement of the `nmae` evaluator; `|lfc_hat|/|lfc_ref|` = 0.138/0.168; the
linearisation residual; sign accuracy and credit decomposition by stratum; the L1-optimal global
and per-stratum scales and their cross-harness transfer; the gate-size ladder for `nmae`, `c*`
and sign accuracy; every arm's Δ`nmae` at three gate depths; `A`, `X`, `B` exactly; `pds`
exactly (`cell_eval2`-identical); emitted `n_pred`, `fid`, `jac`, `reach` with twins; the
structural nullity of a constant lfc bias.

**Derived** (a measured ratio carried into a frame it was not measured in):
the challenge's raw `mse` for any arm — the κ ratio is measured on a harness and mapped onto the
board's solved `kappa0 = 0.0756`, `sumD = 0.230010`, `claim = 0.24653`
(`A12-ARTIFACT-AUDIT.md`); the scaled-per-raw multipliers for `nmae` (−1.657) and `mse` (−1.052),
back-solved from the a12→a13 board pair; the challenge's `|S_p| ≈ 2,000` from v14/v15 `fid`
coverage; the BH-relaxation as a proxy for a deeper reference gate.

**Assumed** (stated, not tested):
that relaxing BH on our references produces a gene cohort resembling the challenge's deeper one
(it matches in size and in weighted sign accuracy; it is not the same experiment); that `claim`
stays pinned at the `#247` cap under a direction-only arm at held `A` (a12 sits at 99.75% of the
cap with 8.9% of headroom, so this is safe for small shape moves and *not* for amplitude moves);
that the board's `nmae` responds to the whitener the way both harnesses do — **untested, and the
one place a13 already proved a harness can be 8x wrong in the flattering direction.**

## 11. What this closes

Closed by measurement, on both harnesses: global magnitude calibration; per-gene and per-decile
magnitude calibration (does not transfer, worse than a scalar); expression-tilted amplitude;
donor-agreement per-gene reweighting for `nmae`; per-perturbation reallocation (oracle ceiling
`+0.004…+0.009` overall, buildable `+0.0004`); control-normalisation de-offsetting; a constant
directional bias (structurally inexpressible).

Left open, and now the only thing `nmae` can be bought with: **`<sign(lfc_ref), lfc_hat>` on the
reference's DE set — the per-gene direction.** It is the same quantity `reach` and `fid` need, and
`CURRENT_STATE.md` item 3 already names it as the only route to a durable first place. **`nmae`
should be folded into that programme and stopped being treated as a separate target.**

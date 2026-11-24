# 228 — 224's per-gene reliability gain SURVIVES the scorer read-back on cosine and `pds`, and is killed by `nmae`. **NO-GO at every dose, at fixed and at re-optimised α, under every pool, every `k_cos`, and every `jac` transfer.**

**Date 2026-09-04.** Deterministic (seeds `20260904`, `771103`). **No `vcc`, no `cell-eval2`, no
build, no submission, nothing under `submissions/` or `data/`, no production file touched,
nothing written outside `experiments/228-theta-gain-plolo/`.** Peak RSS **8.93 GiB** (one
process, `/usr/bin/time -l`), against a 12 GiB mandate. 11 scripts, **334 scored arm×fold
rows** (256 in the dose/α sweep, 46 amplitude-matched, 32 six-member). Preregistration committed at `d1fc985` **before the first fit**.

---

## 0. The answer in ten lines

1. **The gain is real after the read-back, and it is about half its effect-level size.** Wired as a
   per-gene multiplier on the production effect matrix and scored through `profile_from_effect` +
   the scorer read-back on 162's PLOLO folds, it delivers **+4.55% relative aggregate cosine** and
   **Δ`pds` +0.00484** on the three cancer folds at amplitude-matched α, `pds` **positive 3/3**.
   224's effect-level number was **+8.9%**. **156's warning fires at 1.9×, not at ∞.**
2. **None of it is amplitude.** At fixed α the arm is amplitude-neutral in the space `VCC_ALPHA`
   acts in but delivers **+7.38%** more amplitude in the space the scorer reads. Re-solving α to
   match the shipped stack's delivered norm (216's `@adel`, α 0.982 → **0.889**) leaves the cosine
   and `pds` gains **unchanged or slightly larger**. The direction gain is not amplitude.
3. **`nmae` is where it dies, and the mechanism is exact.** 188's identity, the pool 188
   validated: Δ`nmae` **+0.0373 raw** at the selected dose. Decomposed, `U/D` rises **+0.0436**
   while `M/D` rises **+0.0031** — efficiency `e = ΔM/ΔU = 0.072` against a break-even of **0.5**.
   **93% of the magnitude the gain moves onto the reference's DE pool is unmatched.**
4. **Why: the gain reshapes MAGNITUDES, not signs.** |t|-weighted sign accuracy moves
   **+0.00027**. `nmae`'s `M` is sign-matched `min(|pred|,|truth|)`, so it cannot collect a
   magnitude rebalance; `U` collects all of it. This is §4.1c's board-level
   `corr(cosine, nmae) = −0.394` reproduced **inside our own pipeline as a causal, dose-responsive
   relation**.
5. **Six-member price, three cancer folds, @adel, `k_cos = 1.384`:** `o_pds` **+0.00180**,
   `o_mse` **+0.00227**, `o_nmae` **−0.01043**, `o_fid` **−0.00009**, `o_jac` **−0.00058** →
   **−0.00702**. At the best dose of ten swept it is **−0.00128**. Gate **+0.0019**.
6. **It fails under every free choice, jointly.** Taking the un-validated `lfc89` pool, the larger
   `k_cos`, the absolute `jac` transfer, and the dose chosen on the evaluation folds themselves,
   the best composite anywhere is **+0.00177 — still below the gate.**
7. **All four nulls reproduce their 224 signs through this instrument.** gene-shuffle **−14.56%**,
   moment-match **−0.23%**, label-derangement **−16.28%** (−26.8% at `t = 1`), and the
   wrong-recipient basal rung reads **+1.71% against the right recipient's −6.19%** — 224's exact
   signature, wrong better than right.
8. **`fid` and `jac` are now locally priceable for this whole arm class, and that is reusable.**
   The shipped `rc_callset` emitter **prescribes** `n_pred = 400 + 2,200 = 2,600`, and the 2,200
   are a function of the control cells alone. So the yield term `min(1, n_pred/n_real)` is pinned
   at **1** and only the 400-gene head can move. Measured: Δ`fid` precision **−0.00017**,
   Δ`jac` **−3.3% relative**.
9. **Retention splits exactly as `h1union`'s did.** Cosine leave-top-3-out retains **68.9%**
   (`h1union` 72.7%) with bootstrap `P(Δ>0) ≥ 0.92` on 3/3; **`pds` retains −0.80** — the three
   best targets carry more than the whole Δ`pds` on two of three folds.
10. **The October property is real and it is not enough.** The gain is fitted on non-panel targets
    only; refitting with a disjoint basal-matched synthetic panel additionally withheld moves it by
    **r = +0.998 on 3/3 folds**, and a fit pool of only the three lines carrying **0/300** panel
    targets reproduces it at **r = +0.85–0.88**. Panel-independence is confirmed. The arm still
    does not clear.

---

## 1. Known-answer rows — five, all pass

| # | gate | expected | measured | verdict |
|---|---|---|---|---|
| **KA-1** | `t = 0`, CD4 @ 0.34 | 223 `focus223` `BASE:cd4` per fold | `Δcos = Δpds = ` **0.00e+00** on 4/4 folds | **PASS** |
| **KA-2** | `t = 0`, CD4 @ 0.34 + H1 union @ 0.28 | 223 `focus223` `BASE:base222` | **0.00e+00** on 4/4 folds | **PASS** |
| **KA-3** | the gain is 224's/189's object | 224: `r(gain, theta) = +0.81/+0.90`, `r(gain, cross-line reliability) = +0.938` | **+0.721** (189 `theta`) and **+0.902 / spearman +0.940** (cross-line reliability), on the 5,893 genes every fit line carries | **PASS** |
| **KA-4** | fit ∩ eval targets | 0 | 0, **asserted at runtime on every fold** | **PASS** |
| **KA-5** | `t = 0` per-gene panel-total digest | identical to the shipped stack's | identical (`t = 0` **is** the shipped stack, bit-for-bit) | **PASS** |

KA-1/KA-2 are exact zeros, not tolerances: the `t = 0` path reproduces 223's and 218's numbers to
the last bit before anything new is measured.

---

## 2. The gain, fitted on production's own objects (`gain228.py`)

Leave-one-donor-LINE-out inside the fit pool, **NON-PANEL targets only**, on the **official
18,533-gene axis**, through `E062.prepare_main_source` (production's loader, its per-line
`EFFECT_SCALE/raw_median` rescale, its mean-removal, its zeroing of the 300 panel coordinates),
`E062.build_baseline` (production's consensus) and the **fold's own** whitening operator at
`gamma = 0.65`, `lam = 0`. Nothing is remapped between fit and deploy.

| fold | fit lines | fit targets (non-panel) | eval targets (panel) | live genes | full-coverage genes |
|---|---:|---:|---:|---:|---:|
| `k562_gwps` | 5 | 1,767 | 269 | 17,949 | 6,287 |
| `orion_hct116` | 5 | 1,790 | 290 | 17,964 | 5,893 |
| `orion_hek293t` | 5 | 1,782 | 294 | 17,964 | 5,893 |
| `kolf2` | 6 | **1,758** | 273 | 17,964 | 5,893 |

`kolf2`'s 1,758 is 224 §5's own crosspanel number, reproduced from a different code path.

### 2.1 One estimator bug 224 could not have met, and one prior that is degenerate here

**The coverage bias.** On an intersection gene axis every line observes every gene, so 224's
`diag_fit` is unbiased. On the official axis `Y_ij = 0` for a gene the inner-held line does not
carry means **unobserved**, not "the response is zero", and summing `P_ij²` over those rows biases
`g_j` down by the gene's coverage fraction — a coverage weight wearing a reliability weight's
clothes. Uncorrected, the median gain was `0.0006` with 44% of genes negative. `accum_masked`
restricts each inner fold's contribution to the genes that fold's label line observes.

**The prior.** 224's `diag_fit` shrinks toward **1**. The natural scale of `g` is the transfer
cosine — ~0.20 on 224's raw surface, but **~0.02** in the whitened challenge-panel frame — so a
prior at 1 is 50× away from the data and the shrinkage grid collapses the per-gene spread in one
step. The consequence is measurable:

| estimator | `r` between gains fitted for DIFFERENT folds | verdict |
|---|---|---|
| `gain` — 224's verbatim shrink-toward-1 | **−0.205 to +0.940** (median +0.32) | selection is noise; on `orion_hct116` it picks `shrink = 30` and the multiplier is nearly flat |
| **`gain_reg`** — the same family shrunk toward the **pooled coefficient** | **+0.854 to +0.956** (spearman +0.67 to +0.88) | **reproduces 224's `r = 0.95` stability signature** |

**`gain_reg` is the arm.** `gain` is reported throughout as the literal 224 estimator; it
under-delivers (+1.32% LOFO against +4.33%) because its hyperparameter selection is unstable, not
because the object is different: where the two land in the same shrinkage regime they are the same
vector (`k562_gwps`, `r = +0.996`), and where they do not they are unrelated (`orion_hct116`,
`r = −0.022`, the fold where `shrink = 30` was selected and the multiplier went flat).

### 2.2 Mechanism — it is 189's `theta`, confirmed on this axis (`rel228.py`, `mech228.py`)

On the 5,893 genes every fit line carries, held `kolf2`, 6 lines, 1,758 non-panel targets:

| predictor (224's own definitions) | 224 §7 | **228** |
|---|---:|---:|
| per-gene **cross-line reliability** | **+0.938** (spearman +0.965) | **+0.902** (spearman **+0.940**) |
| `sd_line / sd_consensus` | −0.941 | −0.634 |
| 189's `theta` | +0.807 / +0.897 | **+0.721** |
| basal expression | +0.689 | −0.071 (full axis) / +0.64 vs the fold's `b_ctrl` |

On the **full** 18,533-gene axis `r(gain, theta)` is only **+0.10** Pearson (+0.48 Spearman) — the
11,000 partially-covered genes wash it out. That is a property of the axis, not of the object, and
it is why the full-coverage restriction is the honest comparison to 224.

---

## 3. The PLOLO price at production α (`arms228.py`, `analyse228.py`)

Marginal over the **shipped stack** (`base222` = three donors + PerturbMatch 50 + CD4 `f = 0.34`
+ H1 union `g = 0.28`, board **0.2118310**). Relative **aggregate** cosine `cos_agg`, three cancer
folds, `kolf2` apart.

| variant | `t` | k562 | hct116 | hek293t | **CANCER rel %** | **CANCER Δ`pds`** | kolf2 rel % |
|---|---:|---:|---:|---:|---:|---:|---:|
| **`gain_reg`** | 0.05 | +3.63 | +1.91 | +2.46 | **+2.667** | **+0.00316** | +2.50 |
| | 0.10 | +3.07 | +3.21 | +4.03 | **+3.437** | **+0.00405** | +3.92 |
| | **0.25** | +2.46 | +5.04 | +5.48 | **+4.328** | **+0.00470** | +4.91 |
| | 0.50 | +2.20 | +4.75 | +4.56 | +3.838 | +0.00351 | +5.02 |
| | 1.00 | +2.04 | +2.42 | +1.92 | +2.126 | +0.00010 | +5.10 |
| `gain` (224 verbatim) | 0.50 | +2.64 | +0.16 | +5.34 | +2.715 | +0.00320 | +5.16 |
| `gain_fullcov` | 0.50 | −1.49 | +0.07 | +5.48 | +1.358 | +0.00198 | +2.83 |

**Leave-one-FOLD-out dose selection** (218/207/223's rule — the dose is chosen on the other two
cancer folds, the price is read on the held one). For `gain_reg` the LOFO pick is **`t = 0.25` on
all three folds**, so the honest reading equals the plain one: **+4.328%, Δ`pds` +0.00470,
`pds` positive 3/3**.

`kolf2` is reported apart throughout and is **excluded from every headline** — it is pluripotent,
its `cos_agg` base is 0.0141 against the cancer folds' 0.024–0.038, and it inflated a coordinator
prediction on 2026-09-04 badly enough to force a band reissue.

---

## 4. THE `VCC_ALPHA` INTERACTION — the crux the mandate named, and it cuts the other way

The gain preserves the median row norm of the **effect matrix** over `wcols` — production's own
convention (`Operator.apply`, `add_post_whiten_residual`) and the space `VCC_ALPHA` acts in. It is
**not** neutral in the space the scorer reads, because `profile_from_effect` is nonlinear:

> **delivered median row norm +7.38%** on the three cancer folds at `t = 0.25`.

211 pins `claim` at 100.0000% of its cap in the binding context, so that excursion is amplitude we
do not have. 216's `@adel` policy re-solves α by a secant on `log(delivered median row norm)`
against the shipped stack's own value:

| fold | α\* | ratio | `cos_agg` base → arm | rel % | Δ`pds` |
|---|---:|---:|---|---:|---:|
| `k562_gwps` | **0.8133** | 1.00000 | 0.038339 → 0.039446 | +2.89 | +0.00123 |
| `orion_hct116` | **0.9303** | 1.00003 | 0.031236 → 0.032840 | +5.13 | +0.00226 |
| `orion_hek293t` | **0.9239** | 1.00003 | 0.024323 → 0.025684 | +5.59 | +0.01102 |
| **cancer mean** | **0.8892** | | | **+4.55** | **+0.00484** |
| `kolf2` (apart) | 0.8873 | 1.00002 | 0.014130 → 0.014820 | +4.89 | +0.00170 |

**How much of the gain is amplitude: none of it.** Matched: +4.55% / +0.00484. Unmatched:
+4.33% / +0.00470. The direction gain is very slightly **larger** once the amplitude excursion is
removed. Independently, the α sweep at fixed dose shows the gain's relative advantage *rising* as α
falls (+2.45% at `a = 0.85` → +1.83% at `a = 1.15`), which is the opposite of an amplitude proxy.

**Every number after this section is quoted at `@adel`**, because that is the configuration in
which 187's and 211's amplitude optimum is preserved.

---

## 5. `fid` and `jac`, priced — and a structural fact worth more than this arm

§8.4's lesson is that a direction arm priced `mse`-only is understated ≈2.51× because four members
are carried at zero. Here they are not carried at zero, and the reason they are *reachable* is a
property of the shipped emitter that had not been written down:

> **`rc_callset` PRESCRIBES the call set**: `head = 400` genes by `|effect|` inside the scorer's
> 5-CPM universe, then `2,200` more by the CONTROL's mean per-cell CPM, `n_pred = 2,600` exactly,
> enforced by `z_cap = 1.90` (`scripts/make_submission_pm_cd4.py:787-789`, `rc_emit.CallSet`).

Three consequences for **any** profile-class arm under the current emitter:

1. **`n_pred` cannot move.** It is a build constant. With the reference's `n_real` median 89,
   `min(1, n_pred/n_real)` is saturated at **1**, so `fid` **is** the sign precision of the call set
   and nothing else.
2. **2,200 of the 2,600 called genes cannot move** — that block is a function of the control cells
   alone (`ctrl228.py` measures the DE universe at 9,929 / 9,626 / 10,124 genes for A/B/C).
3. **Only the 400-gene head and the signs can move**, and both are measurable on PLOLO.

Measured, three cancer folds, @adel, `t = 0.25`:

| member | quantity | Δ | route |
|---|---|---:|---|
| `fid` | sign precision of the 2,600-gene call set | **−0.00017** | yield pinned at 1 → `o_fid` = **−0.00009** |
| `jac` | proxy Jaccard against the reference's top-89 by `|truth|` | **−0.000226** (**−3.3% relative**) | `o_jac` = **−0.00058** rescaled to the board's own raw `jac` 0.04096 |
| `reach` | ordering half: `|truth|`-weighted sign precision of the top-89 under key `|H|` | **−0.00408** | yield half carries `n_real`, unreadable (rule 9); composite quoted with and without |

`n_real` sensitivity: Δ`jac` is **−0.000226 / −0.000454 / −0.000789** at `n_real` = 89 / 200 / 400
— negative at every depth, and more negative the deeper the reference's DE set is taken.

---

## 6. `nmae` — the member that kills it, with the mechanism (`six228.py`, `adel228.py`)

188's exact identity, `nmae = mean_p (1 + (U_p − 2M_p)/D_p)`, `D` held at the incumbent's, on the
three cancer folds at `@adel`:

| pool | Δ`nmae` at `t = 0.02` | `0.05` | `0.10` | **`0.25`** | `0.50` | `1.00` |
|---|---:|---:|---:|---:|---:|---:|
| **`wald89`** (`|lfc|·√cpm`, **the construction 188 validated on h1-test**) | +0.0108 | +0.0175 | +0.0260 | **+0.0373** | +0.0469 | +0.0534 |
| `lfc89` | +0.0028 | +0.0044 | +0.0070 | +0.0076 | +0.0082 | +0.0083 |
| `h89` | +0.0110 | +0.0180 | +0.0215 | +0.0360 | +0.0437 | +0.0492 |

**Positive — worse — at every dose in every pool.** 188 warned the construction can invert the
sign; here it does not, it only changes the magnitude by 5×.

### The decomposition, which is the actual finding

At `t = 0.25`, `wald89`:

```
U/D   0.14504 -> 0.18861     dU = +0.04356    (predicted magnitude on the reference's DE pool)
M/D   0.08148 -> 0.08461     dM = +0.00313    (sign-matched min(|pred|,|truth|))
e = dM/dU = 0.072                             break-even is 0.5
```

**93% of the magnitude the gain moves onto the reference's DE pool is unmatched.** And the reason
is visible in a second number: the |t|-weighted sign accuracy `q_w` moves **+0.00027**. The gain
reshapes **magnitudes among genes whose signs it does not change** — which is exactly what raises
a cosine (a shape statistic) and exactly what `nmae` cannot pay for, since `M` is a *sign-matched
minimum*. The reliability weight up-weights high-reliability genes, which are high-expression
genes, which is precisely the population `|lfc|·√cpm` selects.

**187's `q_w` route is inapplicable to this arm class and must not be used on it.** It gives
`o_nmae = +0.000002`, i.e. "no effect", because its identity holds *at fixed per-gene magnitude*
and this arm is nothing but a per-gene magnitude change. Recorded so the next agent does not price
a magnitude-reshaping arm on `q_w`.

**Board context, and it points the same way.** §4.1c measured `corr(implied cosine, scaled nmae) =
−0.394` across the 55 teams above overall 0.15. This arm is that anti-correlation reproduced inside
our own pipeline as a *causal, dose-responsive* relation rather than a cross-sectional one.

**The honest caveat, stated as loudly.** This instrument mis-signed `nmae` on this arm family twice
— `cd4f034` (predicted −0.000152, realised **+0.001762**) and `h1union` (predicted −0.000108,
realised **+0.00588**). Both were predictions of order **1e−4**. This one is **+0.0373 raw**,
300× larger, dose-responsive, monotone across ten doses, and reproduced in all three pools. A bias
that flips a 1e−4 prediction does not overturn it.

---

## 7. The six-member price (`price228.py`)

Routes, per rule 12: `pds` = PLOLO 3 cancer folds `Δpds/0.452/6×1.011`; `mse` = 197's inversion
`Δraw = −2.074495·Δcos_board`, `/0.950162566/6`, with `Δcos_board = 0.099469 × rel × k_cos` and
**`k_cos` carried as a RANGE (1.000 on `h1union`, 1.384 on `cd4f034`, matched functional) because
it is not identified**; `nmae` = 188's exact identity `/0.5963/6`; `fid` = call-set sign precision
`/0.30/6` with the yield pinned; `jac` `/0.365/6`; `reach` ordering half, indicative only.

**Overall points, three cancer folds, `@adel`, `wald89`, `k_cos = 1.384`:**

| arm | α\* | rel cos % | `o_pds` | `o_mse` | `o_nmae` | `o_fid` | `o_jac` | **TOTAL** | +reach |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `gain_reg@0.02` | 0.939 | +1.89 | +0.00096 | +0.00095 | −0.00302 | −0.00005 | −0.00012 | **−0.00128** | −0.00113 |
| `gain_reg@0.05` | 0.916 | +2.81 | +0.00122 | +0.00141 | −0.00490 | −0.00006 | −0.00025 | **−0.00258** | −0.00258 |
| `gain_reg@0.10` | 0.904 | +3.60 | +0.00155 | +0.00180 | −0.00726 | −0.00008 | −0.00037 | **−0.00436** | −0.00441 |
| **`gain_reg@0.25`** | 0.889 | **+4.54** | **+0.00180** | **+0.00227** | **−0.01043** | **−0.00009** | **−0.00058** | **−0.00702** | −0.00736 |
| `gain_reg@0.50` | 0.879 | +4.09 | +0.00136 | +0.00205 | −0.01311 | −0.00003 | −0.00069 | −0.01042 | −0.01124 |
| `gain_reg@1.0` | 0.869 | +2.43 | +0.00007 | +0.00122 | −0.01493 | +0.00002 | −0.00081 | −0.01442 | −0.01548 |
| `gain@0.50` (224 verbatim) | 0.937 | +2.83 | +0.00120 | +0.00142 | −0.00959 | −0.00008 | −0.00045 | −0.00751 | −0.00789 |

**Gate: +0.0019** (rule 13, shape-only class; profile-class arm, so the `reach` toll is 0 per 169 —
tested in §9, not assumed).

**The ratio is what makes this robust.** `(o_pds + o_mse) / |o_nmae|` is **0.63 / 0.54 / 0.46 /
0.39 / 0.26 / 0.09** at `t = 0.02 … 1.0`. It is below 1 at the smallest dose measured and falls
monotonically. **There is no dose at which this arm is net positive.**

### The best case over every free choice, taken jointly

Un-validated `lfc89` pool **and** `k_cos = 1.384` **and** absolute `jac` transfer **and** the dose
chosen on the evaluation folds themselves:

| arm | best composite | with reach |
|---|---:|---:|
| `gain_reg@0.02` | +0.00101 | +0.00116 |
| `gain_reg@0.10` | +0.00125 | +0.00119 |
| **`gain_reg@0.25`** | **+0.00177** | +0.00143 |
| `gain_reg@1.0` | −0.00101 | −0.00207 |

**+0.00177 is the maximum reachable anywhere, against a +0.0019 gate.** The arm does not clear
even when every unresolved choice is settled in its favour.

---

## 8. The four falsifiers, re-run through THIS instrument — all reproduce

Same code path, one input changed. `@adel`, three cancer folds, `t = 0.25` unless noted.

| null | 224 | **228** | verdict |
|---|---:|---:|---|
| **gene-shuffle** (`g` permuted over genes) | −7.2% | **−14.56%** cos, Δ`pds` **−0.02701**, 0/3 | fails ✓ |
| **moment-matching** (`sd(Y)/sd(P)`, no target correspondence) | −3.5% | **−0.23%** cos, Δ`pds` **−0.00162**, 0/3 | fails ✓ |
| **label-derangement** | −25.7% | **−16.28%** at `t=0.25`, **−26.78%** at `t=1`, Δ`pds` −0.04156, 0/3 | fails ✓ |
| **wrong-recipient** (224's R4a 12-bin basal gain) | +0.41%, wrong > right | right recipient **−6.19%**; **wrong recipient +1.71%** | 224's signature ✓ |

The per-gene diagonal is recipient-independent by construction, so its own analogue of the
wrong-recipient test is the decoy/`pds` pair — and it passes: `pds` **rises** on 3/3 cancer folds
while the arm's cosine rises, i.e. the gain does not buy similarity by inflating every decoy.

**One thing to flag for whoever reads the null table.** `basal12_wrong` is *positive* on cosine and
`pds` (+1.81%, +0.00580) — it is a coarse 12-bin version of the same reliability weight, so it
carries some of the signal even when binned by the wrong line's basal profile. It fails the
six-member gate for exactly the same reason the arm does (`o_nmae` −0.00465, total **−0.00214**).
Its *relation to the right-recipient rung* — wrong beating right — is what reproduces 184/224.

---

## 9. Rule 13: the digest, and the toll (`six228.py`)

The per-gene panel-total digest is computed on the PLOLO frame's expected totals
(400 cells × the 20,109-count median library × the delivered profile) for every (context, target)
block:

| fold | blocks | `t = 0` | `t = 0.25` |
|---|---:|---|---|
| `k562_gwps` | 807 | identical (it **is** the shipped stack) | **differs on 807/807** |
| `orion_hct116` | 870 | identical | **differs on 870/870** |
| `orion_hek293t` | 882 | identical | **differs on 882/882** |
| `kolf2` | 819 | identical | **differs on 819/819** |

The eight aggregate digests were regenerated in a **second, independent process launch** and
reproduced bit-for-bit, which is this experiment's determinism check.

**Reading, against rule 13's single-valued table.** The digest differs, so the arm is **not**
profile-preserving — row 2 does not apply. But it is **not an emission-law change at all**: the
emitter, the spread policy, `keep_top`, the call-set sizes, `z_cap`, the dispersion and every seed
are untouched; only the effect matrix moves. That is **row 3, profile class** (operator, donors,
magnitude, CD4 dose), where 169 measured the toll at **+0.0006** and it is **0**. The composite in
§7 therefore charges **no `reach` toll**, and the test was run rather than assumed.

---

## 10. Retention and the paired bootstrap (`retain228.py`)

`gain_reg@0.25` at `@adel`, 4,000 paired-target draws.

| fold | Δ`cos_agg` | rel % | **cos retention, leave-top-3-out** | boot 95% CI | `P(Δcos>0)` | Δ`pds` | **`pds` retention** | `P(Δpds>0)` |
|---|---:|---:|---:|---|---:|---:|---:|---:|
| `k562_gwps` | +0.001113 | +2.90 | **0.407** | [−0.00042, +0.00272] | 0.923 | +0.00124 | −3.110 | 0.564 |
| `orion_hct116` | +0.001605 | +5.14 | **0.843** | [+0.00093, +0.00230] | **1.000** | +0.00226 | −0.011 | 0.740 |
| `orion_hek293t` | +0.001363 | +5.60 | **0.817** | [+0.00074, +0.00209] | **1.000** | +0.01102 | **0.716** | **1.000** |
| **cancer mean** | **+0.001360** | **+4.55** | **0.689** | | | **+0.00484** | **−0.801** | |
| `kolf2` (apart) | +0.000690 | +4.88 | 0.647 | [+0.00009, +0.00137] | 0.988 | +0.00170 | −0.952 | 0.610 |

**The cosine gain is distributed** — 68.9% survives dropping the three best targets, against
`h1union`'s 72.7% — **and the `pds` gain is not**: on `k562_gwps` and `orion_hct116` the top three
targets carry more than 100% of it, and its bootstrap `P(Δ>0)` is 0.56 and 0.74. This is
`h1union`'s exact split (100% of the cosine prediction realised, 70% of the `pds` one), and it
means the `pds` term in §7 — the largest positive term — is the **least** trustworthy one.

---

## 11. Marginal over the shipped stack, and 208's rider rule (`stack228.py`)

`gain_reg@0.25`, `a = 1`, three cancer folds, the same gain applied on three bases:

| base | `cos_agg` | + gain | rel % | Δ`pds` |
|---|---:|---:|---:|---:|
| bare whitened consensus | 0.021860 | 0.024202 | **+10.71** | +0.00548 |
| + CD4 @ 0.34 | 0.028332 | 0.029928 | **+5.63** | +0.00585 |
| **+ H1 union @ 0.28 (the shipped stack)** | 0.031299 | 0.032583 | **+4.10** | +0.00470 |

**The shipped stack already supplies ~62% of what the gain would have bought on the bare
consensus.** CD4 alone takes half of it. So even the direction half of this arm is substantially a
free rider on the two residual sources already shipped — 208's rule fires, and the headline number
in §3 is the correctly-discounted one.

---

## 12. §5 — October, and why the property is real but does not save the arm (`october228.py`)

The panel is replaced on 2026-10-22 with 300 targets disjoint from today's, and only the last entry
after that date is prize-eligible. **Every candidate in 225's source ladder is selected against
today's panel** — the source and above all the dose are chosen on the 269–294 panel targets a fold
can adjudicate. This gain is different, and the difference was tested, not asserted:

| test | result |
|---|---|
| the fit never sees a panel response | **by construction** — fit targets are the non-panel intersection of the fit pool, `fit ∩ eval = 0` asserted at runtime on 4/4 folds |
| refit with a **disjoint, basal-matched synthetic 300-gene panel** additionally withheld (1,722–1,743 fit targets instead of 1,767–1,790) | `r(gain_reg real, synthetic)` = **+0.99782 / +0.99819 / +0.99817** on the three cancer folds |
| refit on **only the three lines carrying 0/300 panel targets** (`rpe1`, `nadig_jurkat`, `nadig_hepg2`; 1,980 targets) — 224 §5's strictest control | `r` = **+0.885 / +0.854 / +0.878** |

**So which 300 genes are "the panel" does not enter the gain.** That is a genuinely unusual
property in this project's candidate set and it would have been worth real weight in the Oct-22
window — a gain fitted this way needs no re-selection when the panel changes, only a 30-second
refit against the new operator.

**Two honest limits.** (i) The **dose** `t` *is* selected on today's panel folds — a single scalar,
and LOFO picks the same `t = 0.25` for all three, but it is a panel dependence. (ii) The whitening
operator's gene axis is `union(donor coverage) − panel` and **must** be re-fitted on the day
(§8.3, `docs/DEF-RUNBOOK.md`); the gain must then be re-fitted against the new operator. That is
the runbook's job and is out of scope here; the code path costs ~30 s per fold.

**None of this changes the verdict.** A panel-independent arm that is worth −0.007 overall is worth
−0.007 overall in October too.

---

## 13. Leak audit — computed, not asserted

- **Fit pool excludes the PLOLO held line by cell-line GROUP**, not by table name
  (`K562` / `HCT116` / `HEK293T` / `KOLF2`), so no second library of the held line can enter.
- **`fit ∩ eval = 0`** asserted at runtime on every fold; **`fit ∩ panel = 0`** likewise.
- **`h1_2025` covers 85/85 of h1-test's targets and is inside 222's shipped union.** It is a
  production INPUT on every fold, is never a truth, and **is never in the fit pool**. `h1-test` is
  not used anywhere in this experiment.
- `cd4_rest`, `gse295214_h1` and `h1_2025` enter only as production residual inputs, identically in
  the baseline and in every arm, so no arm can gain from them differentially.
- **`debench.Surface.pds()` is never called and sizes nothing** (188 Phase 3).
- The operator per fold is the leave-the-held-line-out 062 operator for the three cancer folds and
  the shipped `whitener_union_pruned.npz` for `kolf2` (which contributed no row to it) — 162's own
  scheme, unmodified.
- Never touched: `cell-eval2`, any scorer, any operator on disk, `vcc`, `submissions/`, the CD4
  sizing rule, the H1 union sizing rule, or any file outside this directory. `git status` confirms
  the only paths this experiment changed are under `experiments/228-theta-gain-plolo/`.

---

## 14. Verdict, and what is handed on

> **NO-GO.** 224's per-gene reliability gain, wired as a per-gene multiplier on the production
> effect matrix and fitted leave-one-donor-line-out on non-panel targets, delivers a **real,
> distributed, amplitude-free +4.55% relative aggregate cosine and +0.00484 `pds`** on three cancer
> folds through the full read-back — and prices at **−0.0070 overall** (best dose **−0.0013**,
> absolute best case over every free choice **+0.0018**) against a **+0.0019** gate, because it
> moves per-gene magnitude onto exactly the genes `nmae` reads at an efficiency of **0.072**.
> **It does not survive the scorer read-back as an arm, though its direction half does.**

**What this closes.** 224 §10's specified successor is run. The gain is not shippable at any dose,
at fixed or re-optimised α, under any of the three `nmae` pools, either `k_cos`, or either `jac`
transfer. 224's own caution — *"nothing in this experiment justifies spending a slot"* — is
confirmed, and now with the six-member number that makes it a decision rather than a caveat.

**Three things that outlive the arm.**

1. **`nmae` is the binding constraint on any per-gene magnitude reshaping, and `e = ΔM/ΔU` is the
   statistic that decides it.** Break-even is 0.5; this arm reads 0.072. **Any future arm that
   changes per-gene magnitudes must report `e` before its cosine is believed** — and must not be
   priced on 187's `q_w` route, which is identically blind to it. This is the concrete,
   within-pipeline form of §4.1c's `corr(cosine, nmae) = −0.394`. **227 owns this axis and should
   have this number.**
2. **`fid` and `jac` are locally priceable for every profile-class arm under the current emitter**,
   because `rc_callset` prescribes `n_pred = 2,600` with 2,200 of its members fixed by the control
   cells alone. That retires 216's and rule 9's "unmeasurable" status for this arm class and costs
   one `argsort` per row. `ctrl228.py` + `six228.py:call_set` are the instrument.
3. **The `@adel` amplitude trap generalises.** Preserving the effect matrix's median row norm —
   production's own convention — is **not** amplitude neutrality in the space the scorer reads. A
   per-gene reweighting that preserves it delivered **+7.4%** more amplitude. Every future arm that
   reshapes the effect must be re-solved on the *delivered* norm before its members are read, or it
   will smuggle in an amplitude change 211 has closed.

**A fourth, smaller, worth writing down:** 224's `diag_fit` shrinks toward 1 while the natural scale
of the coefficient is the transfer cosine. On a surface where that cosine is ~0.2 the difference is
cosmetic; in the whitened challenge frame, where it is ~0.02, the prior is 50× off and the
hyperparameter selection becomes noise. Shrinking toward the pooled coefficient is the same family
with a prior on scale and it restores 224's own stability signature (`r = 0.95` across folds).

**Not built, not named.** No submission package, no model name, no slot recommendation beyond the
number. If the coordinator wants the direction half of this, the open question is whether a
**magnitude-neutral** form of the reliability weight exists — one that rebalances the *shape*
without concentrating predicted magnitude on the reference's DE pool. That is not obviously
constructible, because the gain *is* a magnitude.

---

## Files

| file | what |
|---|---|
| `PREREGISTRATION.md` | gate, wiring, four falsifiers, stop rule — committed `d1fc985` before the first fit |
| `gain228.py` | the LOLO per-gene gain on the official axis, both shrinkage families, all four nulls |
| `ctrl228.py` | the scorer's 5-CPM DE universe and control-CPM ranking per challenge context |
| `arms228.py` | the PLOLO dose grid + the α sweep, with the two exact known-answer rows |
| `six228.py` | `nmae` (188's identity, 3 pools), `reach` ordering half, the prescribed call set → `fid`/`jac`, the digest |
| `adel228.py` | the amplitude-matched α solve and all six members re-read there |
| `price228.py` | the six-member composite, `k_cos` as a range, best-case-over-everything |
| `retain228.py` | leave-top-3-out retention on cosine and `pds` separately + the paired bootstrap |
| `stack228.py` | 208's rider rule: the marginal as the stack grows |
| `october228.py` | the synthetic-panel refit and the 0/300-coverage fit pool |
| `mech228.py`, `rel228.py` | stability across folds, 189's `theta`, 224's cross-line reliability |
| `analyse228.py` | dose curves, LOFO dose selection, the α table |
| `work/*.json`, `work/*.log` | every number above |

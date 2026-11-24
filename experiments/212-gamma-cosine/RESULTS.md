# 212 — `WHITEN_GAMMA` re-swept with COSINE primary, at the re-solved α, on the current emission. **The cliff did not move and gamma stays closed.** Third and final closure.

**Date:** 2026-09-04 · no `vcc`, no `cell-eval2`, no build, no submission, nothing written
outside `experiments/212-gamma-cosine/`. Peak RSS 7.79 GiB (`plolo212.py`), all seven scripts
under 8 GiB against a 12 GB ceiling. Method and audit trail: [`METHOD.md`](METHOD.md).

---

## Verdict in eight lines

1. **The mandate's two numbers are not the same number, and the difference is one-sided.**
   Measured on 162's PLOLO, four folds, 300-target panel, challenge frames: the whitening
   operator's *entire* cost on the transfer cosine is **1.3575×** (γ 0.65 → 0). 197's gap is
   **1.655×**. **Abolishing the operator outright reaches cosine 0.10758 — still 1.219× short
   of the board's best direction — and costs −0.0687 raw `pds` = −0.02561 overall.** The gap
   is not inside the operator.
2. **170's ρ cliff did NOT move.** 208's "2.6× further away" is a different cliff. Measured on
   all three emissions at the same delivered profile (per-gene totals asserted identical 8/8),
   `claim` is **bit-equal in contexts A and C** and moves by a fraction of a percent in context B
   alone. Context A is the binding context. **170's and 178's algebra stands unchanged**, and
   experiment 211 reached the same verdict the same day from the shipped artifacts (§2.2).
3. **170 is replicated, and its confound is confirmed.** At α = 0.982 every γ in [0, 0.80] is
   net negative on `mse`, a V with its vertex on 0.65. Removing the confound by re-solving α
   at each γ turns the `mse` leg **positive everywhere below 0.65** — up to **+0.00987** at
   γ = 0.
4. **And `pds` eats it, to 3.8%.** The measured marginal rates along the gamma axis are
   **d(overall `mse`\*)/d(cos) = +0.3655** and **d(overall `pds`)/d(cos) = −0.3517**. Their
   sum is **+0.0138**, i.e. **they cancel to 3.8% of the `mse` leg.** Gamma is a knife-edge,
   and that is a property of the operator's power law against the board's own member scaling,
   not of the operating point.
5. **The two maximisers do differ, exactly as the mandate predicted — and it does not help.**
   Cosine is maximised at **γ = 0** (monotone, no interior optimum, all four folds).
   PLOLO `pds` is maximised at **γ = 0.80** (monotone the other way). The composite sits
   between them at **γ ≈ 0.55–0.65 on `rc_callset` (best +0.00013)** and **γ ≈ 0.525 on
   `cd4f034` (best +0.00042)**.
6. **Under the most generous combination available anywhere in this experiment** —
   `cd4f034`'s cosine leg, h1-test's `pds` (the surface that favours lower gamma), and the
   `nmae` credit with no α charge — the best arm is γ = 0.40 at **+0.00321**. It clears the
   shape-only floor and **fails rule 13's amplitude floor (+0.0089) by 2.8×**. The arm *is*
   amplitude-changing: α must fall 0.982 → **0.9254** at γ = 0.40.
7. **All four falsifiers pass.** The permuted-operator control is worse than the real operator
   at every γ > 0 on **both** channels and reduces to it at γ = 0 (the no-op check); all
   four PLOLO folds including leak-free `kolf2` agree in sign; the identity row is bit-exact
   with the operator sha asserted; and the two-surface criterion is re-run on four modern
   folds rather than the discredited Orion pair.
8. **NO-GO. No slot. Do not change `WHITEN_GAMMA`, and do not reopen it again.** This is the
   third closure (097 §9h, 170, this) and the first on the endpoint the gap is actually
   measured in, at each γ's own α, and on the emission that is on the board.

---

## 0. Phase 0 — known-answer rows. Twelve, all PASS.

| # | gate | expected | measured | verdict |
|---|---|---|---|---|
| 1 | `candidate_pre` sha256 | `0601336fdf1ec5c5…` | **identical** | **PASS** |
| 2 | a15 post-whitening audit, γ 0.65 / λ 0 — 8 scalars | a15's meta | worst rel **1.11e−15** | **PASS** |
| 3 | a19 post-whitening audit, γ 0.80 / λ 0.007 — 6 scalars | a19's meta | rel **0.00e+00** on all six | **PASS** |
| 4 | **the 20 shared effect arrays vs 170's cache** | bit-identity | **20 / 20 BIT-IDENTICAL** | **PASS** |
| 5 | `msekit` on the frozen 095 h1-test artifact | 2.0551018299736388 | **2.055101829973639**, rel 2.22e−16 | **PASS** |
| 6 | per-context `Kc` | 0.005084 / 0.005228 / 0.005919 | identical, rel ≤ 5.7e−5 | **PASS** |
| 7 | challenge-panel cap | 0.24653 | **0.24653** | **PASS** |
| 8 | challenge-panel `claim` at γ 0.65 | 164's 0.24581 / 170's 0.24610 | **0.24608** (−0.01% vs 170) | **PASS** |
| 9 | challenge-panel mean ρ at γ 0.65 | 170's 0.938 / 164's 0.9466 | **0.9384**, per-context **0.9926 / 0.9407 / 0.8821** (170: 0.9922 / 0.9402 / 0.8819) | **PASS** |
| 10 | PLOLO `equal` at γ 0.65 — 11 scalars on four folds | 204's published rows | rel **0.00e+00** on 11 of 11 | **PASS** |
| 11 | PLOLO `equal`'s `max_abs_dev_from_162_pre` | 204's own audit, 4 folds | exact on 4 of 4 | **PASS** |
| 12 | Δcos bootstrap half-width | 204's 4.3e−5 / 0.20% relative | **3.90e−5 / 0.178%** | **PASS** |

**Gate 4 is the one 170 could not run, and it is the strongest.** The whole 300 × 18,533 and
85 × 18,533 effect arrays at nine gammas, both surfaces, reproduced **bit for bit** from 170's
cache. That pins the operator, the CD4 solve, the gene axis and the BLAS path together in one
assertion. 170's own trap — `candidate_pre` reproduces only at `--chunk-rows 128`, because
PerturbMatch's centring mean accumulates per chunk — is carried in `gate212.py`'s abort message.

**Two more known-answer rows fell out of the sweep itself.** 178 §7 published PLOLO Δ`pds` for
gamma: **−0.00930** at γ 0.45 and **−0.02386** at γ 0.30. This experiment measures **−0.009297**
and **−0.023857**. Neither was targeted; both reproduce to four decimals.

### The one gate that does not reproduce exactly, and by how much

| quantity | 197's value | here | gap |
|---|---:|---:|---:|
| implied transfer cosine at γ 0.65 | 0.079249 | **0.077305** | **−2.45%** |

197 solved `X` against a15's board point using 191's per-context `claim` measured on the
shipped artifact and `Ā` = 0.233158; this experiment solves it against the same board point
using its own measured `claim` = 0.24608 and `Ā` = 0.231662. **−2.45% is inside 197's own
documented ±2% across-family systematic** (§1.2) and it is a level, not a ratio: every gamma
result below is a ratio to γ 0.65 and is unaffected. Reported rather than reconciled.

---

## 1. Rule 11 — path and sha256 for every operator, never a description

`optable212.py` resolves **17 operator files with 16 distinct digests** (`095/out` and
`096/out-alpha` are the same bytes). The five this experiment loads:

| role | path | sha256[16] | max eig | n |
|---|---|---|---:|---:|
| **production** | `data/donors/whitener_union_pruned.npz` | **`e2f65ecf69b7e0b0`** | 5.894570 | 17,942 |
| PLOLO fold `k562_gwps` | `experiments/062-…/whiteners/held_k562_gwps.npz` | `13d988fcc3b129d2` | 8.4259 | 17,806 |
| PLOLO fold `orion_hct116` | `…/held_orion_hct116.npz` | `7cecbf7f36e40d22` | 8.0363 | 17,942 |
| PLOLO fold `orion_hek293t` | `…/held_orion_hek293t.npz` | `97e486e0d1a805e1` | 3.6958 | 17,942 |
| PLOLO fold `kolf2` | `data/donors/whitener_union_pruned.npz` | `e2f65ecf69b7e0b0` | 5.894570 | 17,942 |

Production's basis is **29,394 × 17,942** (193's correction; §4.2's stated 12,000 × 7,760
matches no file on disk) and its condition number is 1e+12 as the eigenvalue floor is
5.895e−12. **The sha is asserted before the file is loaded** in all three scripts that load
one, because 170's warning stands: the identity check at production settings reproduces the
baseline by orthonormality and cannot detect a mismatched basis.

Two near-neighbours worth naming so nobody grabs them: `data/donors/whitener_union.npz`
(`5608e807…`, max eig 6.2567) and `experiments/199-def-runbook/work/whitener_union_pruned_def199.npz`
(`68dd20dc…`, max eig **5.9389** against production's **5.8946** — a 0.75% difference in a file
whose name differs by one suffix).

---

## 2. Phase 1 — the two cliffs, re-located on all three emissions. **THE CRUX.**

`mse`'s per-context numerator carries **two** saturations and they belong to different knobs:

```
num_c(s) = s²A_c + Kc_c − min( claim_c , s²B_c ) − 2 s X_c ,
claim_c  = (300/(G−1)) · mean_p min( Ĉ_p , K·C_p ) ,   K = PRED_TRACE_CAP_K = 1.0
```

* **(a) the ρ = 1 cliff** — `s²B_c` crossing `claim_c`. **This is the one gamma moves**, because
  gamma moves `A` and `B` (14% over [0.65, 0.35]). Its position is set by `claim_c`.
* **(b) the claim-saturation cliff** — `Ĉ_p` falling back through `K·C_p`. **This is the one the
  emission moves**, and it is the one 208 measured at 1.125 → 2.940.

### 2.1 Measured, three emissions, same delivered profile (`cliff212.py`)

Per-gene totals asserted identical **8/8 on every row** — the profile is pinned by construction
and proved, not claimed, so `A`, `B` and `X` cannot move and only `Ĉ_p` can.

| ctx | emitter | Ĉ/C_p (min) | rows under cap | **`claim_c`** | `s_cliff` = √(claim/B) at γ 0.65 |
|---|---|---:|---:|---:|---:|
| **A** | `quasi` (a15) | 1.1798 (1.1668) | 0 / 8 | **0.235434** | **1.0037** |
| **A** | `rc_hybrid` | 2.7148 (2.1393) | 0 / 8 | **0.235434** | **1.0037** |
| **A** | `rc_callset` (incumbent) | 2.5834 (2.0245) | 0 / 8 | **0.235434** | **1.0037** |
| B | `quasi` | **0.9974** (0.9884) | **5 / 8** | 0.238104 | 1.0318 |
| B | `rc_hybrid` | 2.5873 (2.0547) | 0 / 8 | **0.239126** | 1.0340 |
| B | `rc_callset` | 2.4514 (1.9288) | 0 / 8 | **0.239126** | 1.0340 |
| C | `quasi` | 1.0863 (1.0790) | 0 / 8 | **0.265026** | 1.0648 |
| C | `rc_hybrid` | 3.3044 (2.6241) | 0 / 8 | **0.265026** | 1.0648 |
| C | `rc_callset` | 3.0311 (2.3570) | 0 / 8 | **0.265026** | 1.0648 |

**208's by-product reproduced through an independent route** (170's `profile_from_effect`
staging rather than 196's `recon196` reconstruction): `quasi` **1.0878** against 208's 1.125
(−3.3%), `rc_callset` **2.6886** against 208's 2.940 (−8.6%). Same order, same ordering, same
conclusion. **The premise is verified.**

### 2.2 The answer, stated as the mandate asked it to be

> **208's "2.6× further away" refers to cliff (b), the claim-saturation cliff. Its effect on
> cliff (a), the ρ = 1 cliff, is a fraction of a percent in context B and EXACTLY ZERO — to the
> bit — in contexts A and C. Context A is the binding context: it has the smallest `s_cliff`
> (1.0037 against 1.0318 and 1.0648) and therefore sets the global α optimum. The cliff that
> 170 and 178 measured has not moved. Their algebra stands, and gamma stays closed on the old
> mechanism.**

The mechanism is arithmetic. In every emission `Ĉ_p ≥ K·C_p` on almost every row, so `claim_c`
is pinned at the cap and cannot move. The only context where a15's emission had rows *under*
the cap is B (5 of 8 here, Ĉ/C_p = 0.9974 — sitting exactly on the knee), and `rc_callset`
lifting it to 2.45 recovers them. **Contexts A and C are bit-identical across all three
emissions** — `claim_A` = 0.23543407548420192 and `claim_C` = 0.2650257970493955, the same
float three times.

> ✅ **Independently confirmed by experiment 211, published the same day — and its magnitude is
> the better one.** 211 asked the same question on the α axis and reached the same verdict —
> *"the cliff did NOT move; 208's 2.6× is a different cliff"* — measuring
> `q_c = mean_p min(Ĉ_p, C_p)/C_p` **on the two shipped artifacts**: a15
> **1.000000 / 0.999904 / 1.000000**, `rc_callset` **1.000000 / 1.000000 / 1.000000**, so the
> cliff moves **A +0.0000% / B +0.0048% / C +0.0000%**.
>
> **Where the two disagree, 211 is right.** This experiment reads context B's shift as ~+0.2% of
> the cliff position against 211's +0.0048%, a factor of 40, and the cause is on this side:
> 170's `chal.py` route, which this inherits, emits with `libs = median(control library)` for
> all 400 cells and reads 8 of 300 targets, so its `Ĉ_B` lands at 0.9974 of the cap where the
> real per-cell library distribution puts it at 0.999904. **Both routes agree that A and C are
> exactly pinned and that A binds** — which is the load-bearing part — and 211's number is
> measured on the bytes the board actually scored. **Use 211's.**
>
> The two also disagree by 0.5% on the re-solved global α: **0.98567** here (context A's cliff
> at `s` = 1.00374) against 211's **0.9808** on `rc_callset`'s profile, with 187's 0.9823 before
> either. **α is 211's territory and its number is the one to use.** All three agree the
> headroom is worth **~1e−4 overall**, and nothing in §3 depends on which is right: every gamma
> row is priced against the *same* baseline at *its* own α, so the level cancels out of every
> delta.

### 2.3 The ρ cliff across the whole gamma grid, a15 geometry (`chal212.py`)

| γ | `A` | `B` | `claim` | ρ (A / B / C) | `s_cliff` (A / B / C) | Ĉ/C_p (A / B / C) |
|---:|---:|---:|---:|:---:|:---:|:---:|
| 0.00 | 0.28726 | 0.28562 | 0.24570 | 1.000/1.000/1.000 | 0.904/0.919/0.958 | 1.164/0.990/1.077 |
| 0.30 | 0.26865 | 0.26716 | 0.24599 | 1.000/1.000/1.000 | 0.934/0.955/0.989 | 1.170/0.994/1.080 |
| 0.35 | 0.26470 | 0.26324 | 0.24603 | 1.000/1.000/1.000 | 0.940/0.962/0.997 | 1.170/0.994/1.081 |
| 0.40 | 0.26345 | 0.26199 | 0.24601 | 1.000/1.000/1.000 | 0.942/0.965/0.999 | 1.171/0.994/1.082 |
| 0.45 | 0.25672 | 0.25530 | 0.24612 | 1.000/1.000/0.977 | 0.955/0.978/1.012 | 1.171/0.995/1.082 |
| 0.50 | 0.25103 | 0.24965 | 0.24615 | 1.000/1.000/0.955 | 0.965/0.990/1.023 | 1.172/0.996/1.082 |
| 0.525 | 0.24743 | 0.24607 | 0.24614 | 1.000/1.000/0.942 | 0.972/0.997/1.030 | 1.172/0.996/1.082 |
| 0.55 | 0.24578 | 0.24442 | 0.24613 | 1.000/0.999/0.936 | 0.975/1.001/1.034 | 1.172/0.995/1.082 |
| 0.575 | 0.24405 | 0.24269 | 0.24608 | 1.000/0.992/0.929 | 0.978/1.004/1.038 | 1.173/0.995/1.082 |
| 0.60 | 0.24056 | 0.23923 | 0.24611 | 1.000/0.977/0.916 | 0.985/1.012/1.045 | 1.173/0.996/1.083 |
| 0.625 | 0.23635 | 0.23503 | 0.24608 | 1.000/0.960/0.900 | 0.994/1.021/1.054 | 1.175/0.995/1.083 |
| **0.65 (ships)** | **0.23166** | **0.23037** | **0.24608** | **0.993/0.941/0.882** | **1.004/1.031/1.065** | 1.173/0.995/1.083 |
| 0.70 | 0.22401 | 0.22276 | 0.24610 | 0.961/0.908/0.853 | 1.020/1.049/1.083 | 1.175/0.995/1.082 |
| 0.75 | 0.21831 | 0.21708 | 0.24611 | 0.938/0.884/0.831 | 1.033/1.064/1.097 | 1.175/0.995/1.083 |
| 0.80 | 0.20640 | 0.20523 | 0.24617 | 0.888/0.834/0.786 | 1.061/1.095/1.128 | 1.176/0.996/1.083 |

`Ĉ/C_p` is **flat in gamma to the third decimal** on all three contexts. Gamma cannot touch
cliff (b), which is the other half of §2.2's answer. This is the first time this column has been
run in **a15's geometry** (no mass fix) — the incumbent's own, since the production builder
*refuses* `VCC_MASS_FIX` together with `--emitter rc_hybrid|rc_callset`
(`scripts/make_submission_pm_cd4.py:2755`). 170's whole challenge column was a23's. The two
agree on ρ at γ 0.65 to 0.04%.

---

## 3. Phase 2 — the sweep, COSINE primary

### 3.1 The primary instrument: PLOLO, four folds, per-fold cosine (a15 render)

| γ | k562_gwps | orion_hct116 | orion_hek293t | **kolf2 (leak-free)** | **macro** | ratio | Δcos 95% CI (paired, 4,000 draws) |
|---:|---:|---:|---:|---:|---:|---:|:---|
| 0.00 | 0.035899 | 0.038140 | 0.035320 | 0.009240 | 0.029650 | **1.3575** | [+0.005440, +0.010295] |
| 0.30 | 0.034287 | 0.030776 | 0.027316 | 0.008589 | 0.025242 | 1.1557 | [+0.002584, +0.004270] |
| 0.35 | 0.033942 | 0.029855 | 0.026318 | 0.008471 | 0.024647 | 1.1284 | [+0.002155, +0.003493] |
| 0.40 | 0.033596 | 0.029017 | 0.025402 | 0.008358 | 0.024093 | 1.1031 | [+0.001749, +0.002786] |
| 0.45 | 0.033254 | 0.028251 | 0.024566 | 0.008250 | 0.023580 | 1.0796 | [+0.001367, +0.002133] |
| 0.50 | 0.032911 | 0.027550 | 0.023797 | 0.008137 | 0.023099 | 1.0575 | [+0.000996, +0.001533] |
| 0.525 | 0.032742 | 0.027221 | 0.023441 | 0.008085 | 0.022872 | 1.0472 | [+0.000822, +0.001252] |
| 0.55 | 0.032576 | 0.026911 | 0.023098 | 0.008030 | 0.022654 | 1.0372 | [+0.000650, +0.000982] |
| 0.575 | 0.032412 | 0.026607 | 0.022770 | 0.007974 | 0.022441 | 1.0274 | [+0.000482, +0.000721] |
| 0.60 | 0.032250 | 0.026313 | 0.022452 | 0.007923 | 0.022235 | 1.0180 | [+0.000317, +0.000472] |
| 0.625 | 0.032091 | 0.026031 | 0.022148 | 0.007871 | 0.022035 | 1.0088 | [+0.000156, +0.000232] |
| **0.65** | 0.031934 | 0.025759 | 0.021854 | 0.007821 | 0.021842 | 1.0000 | — |
| 0.70 | 0.031627 | 0.025251 | 0.021298 | 0.007722 | 0.021475 | 0.9832 | [−0.000437, −0.000299] |
| 0.75 | 0.031334 | 0.024755 | 0.020776 | 0.007621 | 0.021122 | 0.9670 | [−0.000852, −0.000590] |
| 0.80 | 0.031050 | 0.024293 | 0.020298 | 0.007523 | 0.020791 | 0.9519 | [−0.001238, −0.000867] |

**Every Δcos is resolved with P(> 0) = 1.000** (the bootstrap is paired, so the sign is
certain). The tightest arm's half-width is **3.90e−5 absolute = 0.178% relative**, reproducing
204's published 4.3e−5 / 0.20% — the instrument is the one 204 calibrated.

**The cosine has no interior maximum. It is monotone decreasing in gamma on all four folds.**

### 3.2 The two maximisers, side by side — they differ, and it does not help

| endpoint | instrument | maximiser | value at the maximiser |
|---|---|---|---|
| **transfer cosine** | PLOLO, 4 folds, a15 render | **γ = 0** (monotone, boundary) | 1.3575× the incumbent |
| **`pds`** | PLOLO, 4 folds, a23 render | **γ = 0.80** (monotone, boundary) | +0.00222 raw |
| `pds` | h1-test, `cell_eval2.discrimination_score` | γ = 0.575 (interior) | +0.00196 raw |
| `nmae` | h1-test, reference BH pool | γ = 0.30 (monotone) | −0.00412 raw (better) |
| **composite, `rc_callset`** | all of the above | **γ ≈ 0.55–0.65** | **+0.00013 overall** |
| **composite, `cd4f034`** | all of the above | **γ ≈ 0.525** | **+0.00042 overall** |

**The mandate's hypothesis is CONFIRMED in structure and REFUTED in consequence.** The cosine
optimum and the `pds` optimum are not the same point — they are at opposite ends of the range.
But the composite optimum barely moves, because the two members trade at almost exactly the
board's own exchange rate (§3.4).

**The h1-test / PLOLO `pds` disagreement is real and 178 §7 already adjudicated it:** the two
surfaces disagree in sign for gamma above 0.575. PLOLO is the one to believe (300 real panel
targets, ~290 decoys, deep pseudobulk truth, leakage-free operators). Both are carried below,
and the verdict is the same under either.

### 3.3 The sweep priced, at fixed α and at the RE-SOLVED α (`alpha212.py`)

`min(claim_c, s²B_c)` is evaluated **per context at every (γ, s)** — 170's own correction, that
linearising the `min` around the shipped point inverted its result. `s*` is found by a dense
8,001-point scan of the exact objective, not by assuming which kink binds.

| γ | board cos | **raw `mse` @ α 0.982** | Δovr | **α\*** | **raw `mse` @ α\*** | **Δovr\*** | Δovr `pds` | Δovr `nmae`\* | **TOTAL\*** |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.104940 | 0.969997 | −0.01688 | 0.8879 | 0.816982 | **+0.00987** | −0.02561 | — | **−0.01574** |
| 0.30 | 0.089339 | 0.929007 | −0.00969 | 0.9167 | 0.848659 | +0.00431 | −0.00890 | +0.00030 | −0.00428 |
| 0.35 | 0.087232 | 0.917590 | −0.00768 | 0.9234 | 0.852932 | +0.00356 | −0.00678 | +0.00031 | −0.00291 |
| 0.40 | 0.085274 | 0.916836 | −0.00755 | 0.9254 | 0.856916 | +0.00287 | −0.00493 | +0.00018 | −0.00189 |
| 0.45 | 0.083458 | 0.902116 | −0.00497 | 0.9373 | 0.860610 | +0.00222 | −0.00347 | +0.00015 | −0.00110 |
| 0.50 | 0.081754 | 0.890985 | −0.00302 | 0.9477 | 0.864086 | +0.00161 | −0.00228 | +0.00015 | −0.00052 |
| 0.525 | 0.080952 | 0.883432 | −0.00169 | 0.9545 | 0.865725 | +0.00132 | −0.00177 | +0.00014 | −0.00032 |
| 0.55 | 0.080178 | 0.881173 | −0.00130 | 0.9576 | 0.867315 | +0.00104 | −0.00133 | +0.00008 | −0.00021 |
| 0.575 | 0.079425 | 0.880841 | −0.00124 | 0.9608 | 0.868867 | +0.00077 | −0.00090 | −0.00002 | −0.00015 |
| 0.60 | 0.078695 | 0.878431 | −0.00081 | 0.9676 | 0.870371 | +0.00050 | −0.00055 | −0.00005 | −0.00009 |
| 0.625 | 0.077989 | 0.875126 | −0.00023 | 0.9761 | 0.871831 | +0.00025 | −0.00022 | −0.00003 | −0.00000 |
| **0.65 (ships)** | **0.077305** | **0.873790** | **0** | **0.9857** | **0.873250** | **0** | **0** | **0** | **0** |
| 0.70 | 0.076006 | 0.878789 | −0.00088 | 1.0018 | 0.875975 | −0.00048 | +0.00040 | −0.00018 | −0.00026 |
| 0.75 | 0.074756 | 0.883055 | −0.00163 | 1.0141 | 0.878636 | −0.00094 | +0.00068 | −0.00036 | −0.00062 |
| 0.80 | 0.073586 | 0.889049 | −0.00268 | 1.0420 | 0.881165 | −0.00139 | +0.00083 | −0.00050 | −0.00106 |

**Three readings.**

1. **170 is replicated.** The fixed-α column reproduces its V with the vertex on 0.65 and the
   same signs everywhere (170 measured −0.00390 at γ 0.50 and −0.00334 at 0.80; this reads
   −0.00302 and −0.00268 on a different geometry and a more generous cosine transfer). 170's
   result was correct.
2. **The confound was real.** Re-solving α turns the `mse` leg positive at every γ below 0.65,
   monotonically. **The mandate's diagnosis of 170's V as partly an α-mismatch artefact is
   CORRECT.** At γ 0.50 the arm moves from −0.00302 to **+0.00161** on `mse` alone.
3. **And it changes nothing, because `pds` moves by the same amount in the other direction.**
   Every TOTAL\* is ≤ 0. The maximum is γ 0.625–0.65 at **+0.00000 to −0.00000** — flat inside
   the instrument's own resolution across [0.55, 0.65], then negative in both directions.

`nmae`\* carries the α charge: 187's board-measured **d `nmae`_raw / d ln α = −0.0416** applied
to the re-solved α. The gamma credit (lower γ improves `nmae`) and the α charge (lower α hurts
it) **cancel to within ±0.0003 overall at every γ**, which is why the column is ~0 throughout.

### 3.4 The mechanism, and why gamma will stay closed

Measured on the two points straddling the incumbent (γ 0.60 → 0.70), in board-overall units per
unit of board cosine:

```
d(overall mse*)/d(cos)  =  +0.3655        <- exactly 8.2's published exchange rate
d(overall pds )/d(cos)  =  -0.3517        <- what whitening charges for it, on the gamma axis
                sum     =  +0.0138        <- 3.8% of the mse leg
```

**On the gamma axis the whitening operator converts bulk cosine into retrieval at 96.2% of the
rate the board pays for bulk cosine.** That is why `WHITEN_GAMMA` is a knife-edge, and it is a
statement about the operator's power law against the board's member scaling — **not** about the
current operating point, the current emission or the current dose. It is the quantitative form
of 197 §3's finding that we extract 26% more retrieval per unit of bulk cosine than the median
shrunk team, seen from the inside: the exponent that buys that extra retrieval is already set
where the marginal trade is fair.

The +3.8% residual tilts toward *lower* gamma. Integrated over the whole range it is worth
**+0.0004 at the very most**, which is what §3.3 and §4 measure directly.

### 3.5 Transfer-route sensitivity — the verdict under all four routes

`X(γ)` is transferred by the **cosine** ratio (the scale-invariant half) rather than by `X`'s own
ratio, which is where this differs from 170. Four routes, Δovr on `mse` at the re-solved α:

| γ | **plolo_a15** (primary) | plolo_a23 | h1 (`cos_deconf`) | x_ratio_170 (170's route) |
|---:|---:|---:|---:|---:|
| 0.35 | **+0.00356** | +0.00365 | +0.00410 | +0.00234 |
| 0.45 | **+0.00222** | +0.00226 | +0.00268 | +0.00098 |
| 0.50 | **+0.00161** | +0.00164 | +0.00199 | +0.00075 |
| 0.55 | **+0.00104** | +0.00106 | +0.00132 | +0.00032 |
| 0.65 | 0 | 0 | 0 | 0 |
| 0.80 | **−0.00139** | −0.00140 | −0.00196 | −0.00208 |

The routes span **2.3×** in magnitude (170's is the most conservative) and agree in sign and
ordering at every gamma. **The verdict is unchanged under all four**, because the `pds` leg
dominates in every one of them.

---

## 4. Phase 3 — six members, the gate, and the composition question

### 4.1 h1-test, all six, at the five gammas 170 did not run (`h1sweep212.py`)

Merged with 170's nine. `nmae` uses `debench.Surface.nmae`, whose gate is **the reference's own
BH < 0.05 significant set** with ≥ 10 genes per row — the significance-selected pool 188
identified as the challenge's construction (median `n_real` **90** here, against the mandate's
~89), not a top-N-by-magnitude pool. Donors: `k562_gwps` + `orion_hct116` + `orion_hek293t`;
`data/donors/h1_2025` covers 85/85 and is not used.

| γ | `mse` | `pds` | `nmae` | `reach` | `fid` | `jac` | `cos_deconf` | ratio | src |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|:--|
| 0.30 | 1.00620 | 0.86863 | **0.95576** | 0.16955 | **0.51083** | **0.06469** | 0.080342 | 1.1744 | 212 |
| 0.35 | 1.01012 | 0.87129 | 0.95603 | 0.18338 | 0.51027 | 0.06465 | 0.078531 | 1.1479 | 170 |
| 0.40 | 1.00449 | 0.87465 | 0.95659 | 0.16787 | 0.50941 | 0.06450 | 0.076792 | 1.1225 | 212 |
| 0.45 | 0.99939 | 0.87549 | 0.95723 | 0.18148 | 0.50830 | 0.06427 | 0.075020 | 1.0966 | 170 |
| 0.50 | 1.00571 | 0.87661 | 0.95769 | 0.17005 | 0.50890 | 0.06424 | 0.073313 | 1.0717 | 170 |
| 0.525 | 1.00744 | 0.87689 | 0.95805 | 0.18073 | 0.50911 | 0.06415 | 0.072484 | 1.0595 | 212 |
| 0.55 | 1.00520 | 0.87661 | 0.95839 | 0.17275 | 0.50878 | 0.06437 | 0.071661 | 1.0475 | 170 |
| 0.575 | 1.00576 | **0.87731** | 0.95889 | 0.16381 | 0.50800 | 0.06433 | 0.070826 | 1.0353 | 212 |
| 0.60 | 1.00853 | 0.87633 | 0.95930 | 0.16297 | 0.50889 | 0.06432 | 0.069998 | 1.0232 | 170 |
| 0.625 | 1.01266 | 0.87661 | 0.95959 | 0.17521 | 0.50869 | 0.06431 | 0.069205 | 1.0116 | 212 |
| **0.65** | 1.01597 | 0.87535 | 0.95988 | 0.17752 | 0.50816 | 0.06425 | 0.068411 | 1.0000 | 170 |
| 0.70 | 1.00304 | 0.87353 | 0.96122 | 0.17711 | 0.50772 | 0.06423 | 0.066854 | 0.9772 | 170 |
| 0.75 | 0.99604 | 0.87367 | 0.96236 | 0.18628 | 0.50777 | 0.06434 | 0.065230 | 0.9535 | 170 |
| 0.80 | 0.96449 | 0.87311 | 0.96400 | 0.17968 | 0.50664 | 0.06397 | 0.063683 | 0.9309 | 170 |

The five new rows interleave smoothly with 170's nine on `cos_deconf`, `nmae` and `A` — the two
tables are one surface. h1-test's own raw `mse` sits at 0.96–1.02, above the board's clamp point
0.98786, so `score_mse` is 0 there at every gamma: **h1-test cannot price `mse`**, which is why
§3.3 exists.

**The three DE members remain unreadable, exactly as 170 measured.** Its re-seeded twin at fixed
gamma moves `reach` by −0.0062 / +0.0054 / −0.0096 raw; the gamma deltas here run −0.0146 to
+0.0088, i.e. **inside or barely outside a two-seed floor estimate of the noise**. `fid` deltas
(+0.0027 … −0.0015) sit at the twin's ±0.0012; `jac` (±0.0005) is inside ±0.0003 twice over. All
three carry `n_real` (rule 9): h1-test's is 90 median / 766 mean.

### 4.2 Rule 13 — the class, the toll, and the floor

**Class: PROFILE (the operator's exponent).** Rule 13's table is single-valued on this: profile-class
arms — operator, donors, magnitude, CD4 dose — carry **no `reach` toll**; 169 falsified it,
measuring +0.0006. The −0.0146 toll belongs to emission-law arms that change the delivered
profile, and this is not one.

**Floor: +0.0089 (amplitude-changing), and it is the right one.** Lowering gamma raises `A` by
**14%** over [0.65, 0.35] on the challenge panel (§2.3), and the arm's whole value depends on
re-solving α — 0.982 → **0.9254** at γ 0.40, a 5.7% amplitude cut. There is no reading on which
this is shape-only. The mandate anticipated this and it is confirmed.

```
predicted = PLOLO d_pds/6 * 1.011  +  096 d_mse/6  -  0 (profile class)
SHIP only if predicted > +0.0089
```

**Best arm anywhere in this experiment: +0.00042. The gate is +0.0089. Short by 21×.**

### 4.3 The composition question — `rc_callset` vs `cd4f034` (`cd4f212.py`)

210's `cd4f034` (`VCC_CD4_TARGET_F` 0.10 → 0.34, built, submitting at 00:00 UTC) re-run through
the identical sweep. The CD4 residual is added **after** the operator and its fraction is
re-solved from a ratio of median row norms, so it interacts with gamma — weakly, but not at zero.

| γ | cos ratio (f = 0.34) | cos ratio (f = 0.10) | Δovr `mse`\* | Δovr `pds` | Δovr `nmae`\* | **TOTAL\*** | (f = 0.10) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.35 | 1.0997 | 1.1284 | +0.00321 | −0.00398 | +0.00031 | −0.00046 | −0.00291 |
| 0.45 | 1.0613 | 1.0796 | +0.00197 | −0.00191 | +0.00015 | +0.00022 | −0.00110 |
| 0.50 | 1.0441 | 1.0575 | +0.00142 | −0.00118 | +0.00015 | +0.00039 | −0.00052 |
| **0.525** | **1.0360** | 1.0472 | +0.00116 | −0.00087 | +0.00014 | **+0.00042** | −0.00032 |
| 0.55 | 1.0283 | 1.0372 | +0.00091 | −0.00060 | +0.00008 | +0.00039 | −0.00021 |
| 0.60 | 1.0136 | 1.0180 | +0.00044 | −0.00024 | −0.00005 | +0.00014 | −0.00009 |
| 0.65 | 1.0000 | 1.0000 | 0 | 0 | 0 | 0 | 0 |
| 0.80 | 0.9646 | 0.9519 | −0.00114 | −0.00018 | −0.00050 | −0.00182 | −0.00106 |

**Three findings, and one of them is a corroboration 210 will want.**

1. **The dose raises the absolute PLOLO cosine by +13.83%** at γ 0.65 (0.0218418 → 0.0248622,
   macro over four folds). 197 measured the f = 0.10 residual alone at **+6.8%** on the board;
   f = 0.34 giving +13.8% on an independent instrument is consistent and is the first
   *forward* measurement of the dose's direction gain. It implies a board cosine of **0.088**
   for `cd4f034`, up from 0.0773.
2. **Gamma's leverage on direction gets WEAKER as the dose rises** (ratio at γ 0.50: 1.0575 →
   1.0441) and its `pds` cost gets weaker too (−0.00613 → −0.00315 raw). Both legs shrink, the
   `pds` leg shrinks faster, and the optimum moves from γ ≈ 0.625 to **γ ≈ 0.525**, worth
   **+0.00042** instead of +0.00013.
3. **So the answer to "do they differ, and what does it depend on" is: yes, by 0.10 of gamma and
   by +0.0003 of overall, and it depends on the CD4 residual's share of the delivered
   direction.** Both are far inside the gate. `cd4f034` also moves context A's ρ cliff from
   `s` 1.0037 to **1.0141** (210 measured delivered `A` at 0.97970×), which is 1.03% of α
   headroom and belongs to 211.

*Convention note:* the `cd4f034` `mse` leg is priced through the exchange rate
d(overall `mse`\*)/d(cos) = +0.3655, which is **exact at the re-solved α** because the amplitude
cancels there. The challenge-panel `A`/`B`/`claim` at f = 0.34 were not re-measured (17 min);
given the answer is 21× under the gate, that was not worth the machine time and it is stated
rather than hidden.

### 4.4 The most generous combination available anywhere in this experiment

`cd4f034`'s cosine leg **+** h1-test's `pds` (the surface that favours lower gamma, and the one
178 §7 argued against) **+** the `nmae` credit with **no** α charge:

| γ | Δovr `mse`\* | Δovr `pds` (h1) | Δovr `nmae` | **TOTAL** |
|---:|---:|---:|---:|---:|
| 0.35 | +0.00321 | −0.00150 | +0.00106 | +0.00277 |
| **0.40** | **+0.00256** | **−0.00026** | **+0.00091** | **+0.00321** |
| 0.45 | +0.00197 | +0.00005 | +0.00073 | +0.00275 |
| 0.50 | +0.00142 | +0.00047 | +0.00060 | +0.00249 |
| 0.525 | +0.00116 | +0.00057 | +0.00051 | +0.00223 |
| 0.55 | +0.00091 | +0.00047 | +0.00041 | +0.00179 |

**Best: γ = 0.40 at +0.00321.** It clears the shape-only floor (+0.0019) and **fails the
amplitude floor (+0.0089) by 2.8×** — and the amplitude floor is the one that applies, because
this arm cuts α by 5.7%. Stated so that no future agent has to wonder whether a friendlier
instrument choice would have changed the answer. **It would not.**

---

## 5. Phase 4 — falsifiers

### F1 — 170's two-surface criterion, on four modern folds instead of the discredited pair

097's rule was "no γ is net positive on both h1-test and Orion", and 168d discredited the Orion
pair for large reshapings (3.8–10× overstatement with sign flips). This experiment does not use
the Orion pair; it uses **four PLOLO folds, each with the scored line held out and its own
leakage-free operator**, of which two happen to be Orion lines. Per-fold cosine ratios to γ 0.65:

| γ | k562_gwps | orion_hct116 | orion_hek293t | kolf2 |
|---:|---:|---:|---:|---:|
| 0.00 | 1.1242 | 1.4807 | 1.6162 | 1.1815 |
| 0.35 | 1.0629 | 1.1590 | 1.2043 | 1.0832 |
| 0.50 | 1.0306 | 1.0695 | 1.0889 | 1.0405 |
| 0.80 | 0.9723 | 0.9431 | 0.9288 | 0.9620 |

**Unanimous in sign on 4 of 4 folds at every gamma, with a 2× spread in magnitude.** `pds` is
likewise unanimous (monotone increasing in gamma on all four). **The two-surface criterion is
satisfied in the sense that matters — every surface agrees — and the arm is still negative,
because the disagreement was never what killed it.** No absolute is quoted from the Orion pair.

### F2 — the leak-free fold, reported

`kolf2` is not a production donor and contributed no row to any production operator, so that
fold uses the exact shipped operator (`e2f65ecf…`) with no leakage. It is the **second-least**
favourable to the "un-whiten for direction" story (1.1815× at γ = 0 against the macro's 1.3575×)
and it agrees in sign with the other three at every gamma. It is reported in §3.1 alongside the
other three, not folded away.

### F3 — the permuted-operator control (`falsify212.py`). **The one to watch, and it passes.**

Two controls, both with the **identical** eigenvalue spectrum, condition number and
amplification: `gene_perm` re-pairs every operator coordinate with a different gene (this is
AGENTS.md §8.3's Oct-22 trap made deliberate — the failure that passes every check the builder
runs), `spectrum_shuffle` keeps the pairing and permutes the eigenvalues among the eigenvectors.
Macro over four folds:

| γ | **real** cos / `pds` (a23 render) | **gene_perm** cos / `pds` | **spectrum_shuffle** cos / `pds` |
|---:|---:|---:|---:|
| 0.00 | 0.029967 / 0.673446 | **0.029967 / 0.673446** | **0.029967 / 0.673446** |
| 0.30 | 0.025376 / 0.718255 | 0.022026 / 0.668290 | 0.024064 / 0.659123 |
| 0.50 | 0.023167 / 0.735987 | 0.011692 / 0.650994 | 0.012728 / 0.648554 |
| 0.65 | 0.021885 / **0.742112** | 0.007454 / **0.627206** | 0.007584 / **0.631031** |
| 0.80 | 0.020821 / 0.744329 | 0.005594 / 0.609139 | 0.005370 / 0.610639 |

* **The no-op check passes.** At γ = 0 every eigenvalue is raised to the power 0, so the
  transform is `(M V) Vᵀ` — the identity by orthonormality — and a permutation cannot change
  it. `spectrum_shuffle` is **bit-equal** to `real` on all four folds and all six scalars.
  `gene_perm` is bit-equal on `pds` and `q_w` and agrees to **7.0e−9 relative** on the cosines
  and `rr`: permuting the columns changes the summation order of a pair of 17,942² GEMMs, which
  is the same last-bit effect 170 documented for the arm64/x86 `final_sha256`. A control that
  failed this check would not be a control.
* **Both controls lose on BOTH channels at every γ > 0**, and never come within reach of the
  real operator anywhere on the grid. **209's outcome — the random control as the maximiser —
  does not occur here.**
* **The mechanism check.** Under a meaningless operator `pds` falls monotonically as gamma
  rises (0.673 → 0.609); under the real one it *rises* (0.673 → 0.744). The `pds` gain that
  makes gamma a knife-edge is structural, not an artefact of the power map. If it were the
  latter, §3.4's cancellation would be a coincidence; it is not.
* At γ 0.65 `gene_perm` costs **−0.1149 raw `pds` = −0.0429 overall** and **−66% of the
  cosine**. That is the price of §8.3's silent failure, quantified, and it is a number the
  Oct 22 runbook can use.

### F4 — the identity row

γ = 0.65 reproduces the incumbent **exactly**: the `real` control's `a23_cos`, `a15_cos`,
`a23_pds`, `a15_pds` and `a15_sgn_w` are bit-equal to the primary sweep on all four folds, and
`gate212.py`'s effect arrays at γ 0.65 are bit-identical to 170's, which reproduced a15's own
audit to 1.1e−15 and a19's `final_sha256` exactly.

**But 170's warning is the operative one and it is honoured:** the numeric identity check at
production settings reproduces the baseline by orthonormality and **cannot** detect a mismatched
basis. So the gate is the **sha**, asserted before the file is loaded, in `gate212.py`,
`plolo212.py` and `falsify212.py`: `e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`.
`optable212.py` shows why that matters — 17 operator files, 16 distinct digests, one of which
(`199`'s D/E/F refit) differs from production in the second decimal of its max eigenvalue with a
filename one suffix away.

---

## 6. Verdict, in overall units

| arm | Δoverall, best available reading | rule 13 floor | verdict |
|---|---:|---:|---|
| γ 0.625, `rc_callset`, α re-solved | **−0.00000** | +0.0089 | **VETO** |
| γ 0.55, `rc_callset`, α re-solved | −0.00021 | +0.0089 | **VETO** |
| γ 0.525, `cd4f034`, α re-solved | **+0.00042** | +0.0089 | **VETO** — short by 21× |
| γ 0.40, most generous combination anywhere | +0.00321 | +0.0089 | **VETO** — short by 2.8× |
| γ 0.35, `rc_callset`, α 0.982 (170's frame) | −0.00768 | +0.0089 | **VETO** |
| γ 0.00 (no operator at all) | −0.01574 | +0.0089 | **VETO** |

**NO-GO. No slot. `WHITEN_GAMMA` stays at 0.65.**

**And do not reopen it.** The reason is §3.4 and it is structural: along the gamma axis the
operator trades bulk cosine for retrieval at **96.2%** of the board's own exchange rate between
`mse` and `pds`. That ratio is a property of the operator's power law and the board's member
scaling. It does not move when the emission changes (§2 measured that), it does not move when
the CD4 dose changes (§4.3 measured that — both legs shrank together), and it will not move when
α is re-solved (§3.3 measured that — the amplitude cancels exactly at the optimum). **A future
agent who finds a reason to re-open gamma needs a mechanism that breaks that 96.2%, not a new
endpoint or a new incumbent.**

---

## 7. What this corrects, closes and hands on

**Corrects.**

- **AGENTS.md §4.2's "the whitening operator does cost 1.5–2.1× on directional accuracy."**
  CONFIRMED in kind, and the number is measured here on the **transfer cosine** — the quantity
  197's inversion and §8.2's exchange rate actually use — as **1.3575× macro**, per fold
  **1.124 / 1.481 / 1.616 / 1.182**. The 1.5–2.1 range comes from 097's `kappa` on h1-test and
  the Orion pair (`docs/CURRENT_STATE.md` line 1433) and is met on 2 of 4 folds. The sentence
  should read **"1.12–1.62× per fold, 1.36× macro, on the transfer cosine."**
- **The mandate's own premise arithmetic.** "The cost of our own whitening and the size of our
  entire remaining gap are the same number" — they are not. **1.358× against 1.655×.** They
  differ by 22%, one-sided: **the whitening cost is strictly smaller than the gap, so the gap
  cannot be inside the operator even if the operator were free to remove.** This is the single
  most useful sentence in the document and it retires the whole line of attack.
- **The reading of 208's by-product that both this mandate and 211's carried** — that `Ĉ/C_p`
  1.125 → 2.940 creates "much more room before `s²B` reaches `claim`". It does not: the room is
  *set by* `claim`, and `claim` is the cap in both emissions. Two independent instruments
  measured the movement at +0.0048% (211, shipped artifacts) and ~+0.2% (here, reconstructed
  emission) in one context of three, and exactly zero in the other two. **Neither is 2.6×.**
- **170's transfer route.** It carried `X`'s h1-test ratio and let the challenge panel's own `A`
  move underneath it, which divides one surface's amplitude growth out by another's. The
  scale-invariant contraction is the **cosine** ratio. The two differ by up to **2.3×** (§3.5).
  170's verdict is unaffected; its `mse` magnitudes below γ 0.55 are conservative by ~2×.

**Closes.**

- `WHITEN_GAMMA`, for the third and final time: on the cosine endpoint, at each gamma's own α,
  on the current emission, on the current dose, with all four folds and two operator controls.
- The claim that 208's by-product reopens gamma. It reopens the **emission** cliff, which gamma
  cannot touch (`Ĉ/C_p` is flat in gamma to three decimals across all three contexts).

**Hands on.**

- **§3.4's cancellation is the reusable object.** `d(overall pds)/d(cos) = −0.3517` along the
  gamma axis, against the board's `+0.3655` for `mse`. Any future arm that buys direction by
  **flattening the geometry** pays that rate and is dead on arrival. **The only direction arms
  worth building are ones that raise the cosine without moving the operator's exponent** — which
  is 170's "norm-preserving, not norm-adding" criterion restated in the retrieval channel, and
  it is what 194's estimator work and the CD4 dose both are.
- **The CD4 dose's forward direction gain: +13.83% of PLOLO cosine at f = 0.34** (§4.3), the
  first forward measurement of it, corroborating 197's +6.8% board measurement at f = 0.10.
  If it transfers, `cd4f034` lands near board cosine **0.088** — and that is a *direction* gain
  of the kind §8.2 says is the whole remaining gap, arriving from the one lever that has ever
  produced one.
- **The Oct 22 trap, priced.** `gene_perm` — the operator pointed at the wrong genes, which
  passes every check the builder runs — costs **−0.0429 overall on `pds` alone** and 66% of the
  cosine (§F3). `docs/DEF-RUNBOOK.md` can quote that as the cost of getting the refit wrong.

---

## 8. What was NOT measured, and why

* **`cd4f034`'s challenge-panel `A` / `B` / `claim` at f = 0.34** — 17 minutes, and the `mse`
  leg is exact at the re-solved α through the verified exchange rate. The answer is 21× under
  the gate; re-measuring would not change it. Stated in §4.3, not hidden.
* **The h1-test six-member column is in a23's geometry** (`h1kit.builder()` forces
  `VCC_MASS_FIX=minnorm` and raises otherwise) while the incumbent ships without the mass fix.
  The control for it is that the PLOLO instrument carries **both** renders and they agree on
  every gamma ordering, and their cosine *ratios* differ by at most **0.87%** (at γ = 0, the
  extreme point; ≤ 0.3% everywhere inside [0.30, 0.80]).
* **`X_c` is taken uniform across contexts**, as 096, 170, 178 and 197 all do — one board
  equation, one unknown. A per-context `X` would move which context binds; context A binds by
  a 2.9% margin in `s_cliff`, so it would take a >2.9% per-context asymmetry in `X` to change
  it.
* **`reach`, `fid`, `jac` are reported as deltas only** and are unresolvable: 170's re-seeded
  twin moves `reach` by up to ±0.0018 overall at fixed gamma, larger than the entire effect.
  Rule 9's `n_real` caveat applies to all three.
* **`lam` was not swept.** 193/198 established production ships `WHITEN_LAM = 0.0`; there is no
  ridge, and the only direction available is upward, which a19 measured at −0.0959.
* **The operator's basis was not refit** (193 closed the eigenvalue map, 195 closed learned
  transforms) and **the donor tables were not touched** (194, 204, 209). Territory respected.
* **No submission was built, no slot spent, nothing written to `submissions/`, `cell-eval2`
  never run, `vcc` never run.** Every scored number comes from `cell_eval2`'s own primitives
  called on cached sufficient statistics, anchored to a real scorer run by gate 5.
* **One file outside this directory was touched and it is a bytecode cache.** Importing 204's
  `score204.py` refreshes `experiments/204-repression-weighted-mixing/__pycache__/score204.cpython-312.pyc`.
  No source, data or result file outside `experiments/212-gamma-cosine/` was written; 162's and
  204's fold caches were **cloned** into `work/plolo/` before `plolo.WORK` was redirected,
  precisely so that a cache miss could not write into a sibling's directory.

Peak RSS: `plolo212.py` **7.79 GiB**, `falsify212.py` 7.58, `cd4f212.py` 7.56, `h1sweep212.py`
6.61, `cliff212.py` 5.39, `gate212.py` 4.30, `chal212.py` 1.49 — all under the 12 GB ceiling,
and no two were run concurrently above 10 GB combined.

---

## 9. Files

| file | what it does |
|---|---|
| `gate212.py` | Phase 0: `candidate_pre` sha, the a15/a19 board-point audits, the 15-gamma operator sweep, and the 20/20 bit-identity gate against 170's cache → `work/gate212.json`, `work/eff/` |
| `plolo212.py` | **the primary instrument**: 162's PLOLO with gamma as an argument, 4 folds × 15 gammas, both renders, per-target cosines → `work/plolo212.json` |
| `boot212.py` | 204's 4,000-draw paired target-cluster bootstrap on Δcos and Δ`pds` → `work/boot212.json` |
| `chal212.py` | the challenge panel's `A`/`B`/`claim`/ρ/`s_cliff`/Ĉ per context per gamma, **in a15's geometry** → `work/chal212_a15.json` |
| `cliff212.py` | Phase 1: Ĉ/C_p and `claim` on three emissions at the same delivered profile, totals asserted 8/8 → `work/cliff212.json` |
| `alpha212.py` | the α re-solve, the four transfer routes, the marginal exchange rates, the six-member price, rule 13 → `work/alpha212.json` |
| `h1sweep212.py` | 170's six-member h1-test sweep at the five gammas 170 did not run → `work/h1sweep212.json` |
| `cd4f212.py` | the same sweep on `cd4f034`'s profile (`cd4_f` 0.34) → `work/cd4f212_summary.json` |
| `falsify212.py` | the permuted-operator and spectrum-shuffle controls, swept → `work/falsify212.json` |
| `optable212.py` | rule 11's operator table: 17 files, 16 digests → `work/optable212.json` |
| `METHOD.md` | the audit trail and the exact diff against 170's and 204's instruments |

**Inputs read, never written:** `experiments/170-gamma-mse/work/{pre,cd4,de_real_h1}.npz` and
`work/eff/`, `experiments/170-gamma-mse/{h1kit,h1sweep}.py`,
`experiments/096-mse-floor/work/{ctrl-stats,cp-proxy-ntc}.npz` and `{msekit,decompose}.py`,
`experiments/162-retrodicting-validator/plolo.py`,
`experiments/204-repression-weighted-mixing/{score204.py,work/plolo/*}` (cloned),
`experiments/153-de-bench/debench.py`, `scripts/make_submission_pm_cd4.py`,
`src/vcc2026/{donor_transfer,rc_emit}.py`, `data/donors/whitener_union_pruned.npz`,
`data/challenge/*`, `data/harness/h1-test.h5ad`. No production code was touched.

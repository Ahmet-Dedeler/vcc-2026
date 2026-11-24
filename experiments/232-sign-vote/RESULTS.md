# 232 — The majority sign is NOT better than the sign of the mean. The last untested construction in the consensus is closed. **NO-GO — and §8.1d is confirmed, causally, with an exchange rate attached.**

**Date 2026-09-04.** Deterministic (seed `20260904`; one arm re-run in a second independent
process reproduced **every** reported quantity and the panel-total digest **bit-for-bit**,
`work/det232.log`). **No `vcc`, no `cell-eval2`, no build, no submission, nothing under
`submissions/` or `data/`, no production file touched, nothing written outside
`experiments/232-sign-vote/`.** Peak RSS **8.82 GiB** (`phase0_232.py`, four folds in one
process, `/usr/bin/time -l`) against a 10 GB mandate; the other `/usr/bin/time -l`-measured runs
peaked at **3.14 GiB** (`prep232.py`) and **5.06 GiB** (`det232.py`, one full arm rebuild).
Every arm-scoring run was one process per fold for exactly this reason.
7 scripts, **282 scored arm×fold rows**, 4 PLOLO folds. Preregistration written before the first
Phase-1 number was read (`PREREGISTRATION.md`).

---

## 0. The answer in thirteen lines

1. **The disagreement rate is NOT tiny — it is a third of the scored DE pool.** |t|-weighted, on
   `wald89`, the covering donors disagree in sign on **36.27%** of the pool's mass on the three
   cancer folds (41.60 / 34.03 / 33.17%) and **64.38%** on `kolf2`. So the mandate's first
   terminal condition does not fire: there is plenty of contested mass.
2. **But `sign(mean)` on that contested mass is a coin flip: 0.52292.** So the whole prize
   available to *any* sign-combination rule — majority, weighted, or **oracle** — is
   `0.3627 × (1 − 0.52292) =` **+0.17325** of |t|-weighted sign accuracy. That ceiling is
   computed before any arm exists, and it is enormous.
3. **The unweighted majority is not a rule at production's PLOLO arity.** Each cancer fold has
   **two** covering donors after the hold-out, so **100% of contested pairs are 1-1 ties** —
   `mass_tie == mass_disagree` to five decimals on 3/3 folds — and the arm's mask is **empty at
   every dose**. `maj_strict`, which refuses to guess, scores **0.324 / 0.395 / 0.394** against
   the mean's 0.538 / 0.572 / 0.571.
4. **Where the majority IS a rule — arity 3 — it is WORSE than the mean.** On `kolf2`
   (all three production donors, leak-free truth) |t|-weighted sign accuracy is
   **0.50725 majority vs 0.50820 mean, −0.00095**. On the cancer folds with a third
   (non-production) voter added it is **0.56057 vs 0.56077, −0.00020**, positive on 1/3 folds.
   Scored end-to-end the majority arm is **−0.00375 overall**.
5. **The one implementable rule that beats the mean does so by 0.00145 — 0.84% of the ceiling.**
   A confidence-weighted vote using a per-(donor, gene) weight `|e_s| / sd_s(gene)` reaches
   **0.56150 vs 0.56005**, positive 3/3 cancer folds.
6. **Priced through all six members it is +0.00097 against a +0.0019 gate.** Best over every free
   choice (both `k_cos`, all three pools, both `jac` transfers, the dose picked on the evaluation
   folds themselves): **+0.00089 … +0.00097**. It fails by **2.0×**, and it is **negative on
   `kolf2`** (−2.11% cosine, Δ`pds` −0.00216, Δ`nmae` +0.00062).
7. **§8.1d's prediction is CONFIRMED, and this is the durable result.** A pure sign improvement at
   **bit-identical magnitude** moves `nmae` the *right* way: over a 14-dose two-sided oracle
   response, `Δnmae = −0.002573` per **+0.01** |t|-weighted sign accuracy, **r = −0.995**, and the
   efficiency `e = ΔM/ΔU` runs **29–40** against break-even 0.5 and 228's **0.072**.
8. **The head-to-head that settles it.** At a *matched* relative cosine gain of ≈ +4.4%,
   228's magnitude arm gave Δ`nmae` **+0.0373** and a composite of **−0.00702**; a pure sign
   route gives Δ`nmae` **−0.00889** and a composite of **+0.01170**. Same cosine, opposite
   `nmae` sign, a **0.0187** swing in overall points, decided entirely by *how* the cosine
   was bought.
9. **The exchange rate, new and reusable: +0.01 of |t|-weighted sign accuracy is worth
   +0.00399 overall** (k_cos 1.384; **+0.00379** at k_cos 1.000), r = 0.975 across 14 doses.
   So the gate needs **Δq_w = +0.0048 … +0.0050**. The vote delivers **+0.00090**.
10. **Which reframes the channel, and this is the number to carry forward.** A sign rule needs
    only **2.75% of the available disagreement-set headroom** to clear rule 13's gate. The best
    implementable one captures **0.52%**. The channel is not closed by arithmetic — it is
    **5.3× short on efficiency**, and the missing ingredient is information the donors do not have.
11. **Two structural facts fall out.** (a) Under the shipped `rc_callset` emitter a
    magnitude-preserving sign arm leaves `jac` **exactly invariant** — max |Δ`jac`| = **0.000e+00**
    over every arm, fold and depth — because the 400-gene head ranks on `|effect|` and the 2,200
    CPM block on the control cells alone. (b) **189's `theta` cannot serve as a vote weight at
    all**: it is per-GENE, so it factors out of `sign(Σ_s w(j)·sign(e_s))` exactly, at every arity.
12. **Verdict: NO-GO. Propose no slot. The last untested construction in the consensus is closed.**
13. **INBOX-1's second question, answered with the same instrument (§9): the whitening operator
    costs `q` −0.01370 (−2.43% relative), worth **−0.00547 overall** — 22% of the `pds` credit it
    buys, and never charged by 170 / 193 / 212 / 226, all of which scored the operator on
    direction metrics only.** It also explains **68%** of whitening's whole `nmae` cost, and the
    `k = 0.1286` it implies under 188's identity matches 187's independently measured 0.12.

---

## 1. Phase 0 — known-answer rows, three, all pass (`phase0_232.py`, `work/phase0_232.log`)

| # | gate | expected | measured | verdict |
|---|---|---|---|---|
| **KA-A** | 228's `nmae` decomposition and `e` for 224's per-gene gain (3 cancer folds, `@adel`, `wald89`) | `U/D` 0.14504→0.18861, `M/D` 0.08148→0.08461, `e = 0.072`, Δ`nmae` +0.0373 | `dU` **+0.04356**, `dM` **+0.00313**, **`e = 0.0718`**, Δ`nmae` **+0.0373**; 188's identity `1 + U/D − 2M/D` reproduces `nmae` to 5e-6 | **PASS** |
| **KA-B** | our own PLOLO baseline vs the shipped `h1union` stack | 223 `focus223` `cd4` and `base222` per fold | **Δcos = Δpds = 0.0e+00** on **8/8** rows (4 folds × 2 stacks), exact zeros, not tolerances | **PASS** |
| **KA-C** | 204: equal donor weighting maximises the transfer cosine | equal-weight `pds` 0.78912 / 0.79182 / 0.77738 / 0.61013; `share_4_1_1` −6.53% rel cosine | **0.789121 / 0.791823 / 0.777375 / 0.610126** (max Δ 2.8e-6); **7/7** alternatives below equal; `share_4_1_1` **−6.53%** | **PASS** |

KA-C reproduces 204's `Δcos` column **to 1e-6 on every arm** by an independent code path:

| arm | w (k562, hct116, hek293t) | 232 Δcos | 204 Δcos |
|---|---|---:|---:|
| `share_4_1_1` | (2.000, 0.500, 0.500) | **−0.001429** | −0.001429 |
| `share_1_2_2` | (0.600, 1.200, 1.200) | −0.000367 | −0.000367 |
| `minvar` | (0.618, 1.208, 1.174) | −0.000329 | −0.000328 |
| `a_self` | (1.246, 1.032, 0.722) | −0.000328 | −0.000327 |
| `k115` | (1.095, 0.952, 0.952) | −0.000045 | −0.000045 |
| `h115` | (0.952, 0.952, 1.095) | −0.000008 | −0.000008 |
| `k085` | (0.895, 1.053, 1.053) | −0.000002 | −0.000002 |

And a fourth, in `prep232.py`: **KA-0** — `E062.build_baseline` over the separated sources
reproduces 162's cached `pre` to **3.2e-8 … 7.4e-8 relative** on all four folds, so the objects
this experiment votes over are the production input and not a reconstruction of it.

---

## 2. Phase 1 — the sign-accuracy table, measured before anything was built (`sign232.py`)

|t|-weighted sign accuracy on the **`wald89`** pool (188's `|lfc|·√cpm`, top 89/row — the
construction 188 validated on h1-test), at the **consensus** level, i.e. on the object
`build_baseline` produces, before whitening. Ties resolved to `sign(mean)` (leave alone).

| fold | arity | **disagree mass** | of which TIES | `sign(mean)` | **majority** | **z-vote** | `maj_strict` | donor 0 | donor 1 | donor 2 | oracle | **ceiling** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `k562_gwps` | 2 | 0.4160 | **0.4160** | 0.53775 | **0.53775** | **0.53914** | 0.32356 | 0.53142 | 0.53175 | — | 1.0 | +0.20160 |
| `orion_hct116` | 2 | 0.3403 | **0.3403** | 0.57163 | **0.57163** | **0.57270** | 0.39543 | 0.55875 | 0.57240 | — | 1.0 | +0.16371 |
| `orion_hek293t` | 2 | 0.3317 | **0.3317** | 0.57077 | **0.57077** | **0.57265** | 0.39444 | 0.55886 | 0.56174 | — | 1.0 | +0.15444 |
| **CANCER mean** | 2 | **0.3627** | **0.3627** | **0.56005** | **0.56005** | **0.56150** | 0.37114 | 0.54968 | 0.55530 | — | 1.0 | **+0.17325** |
| `kolf2` *(apart)* | **3** | 0.6438 | 0.1545 | 0.50820 | **0.50725** | 0.50752 | 0.42452 | 0.50936 | 0.50033 | 0.50285 | 1.0 | +0.31520 |
| `k562_gwps` +kolf2 | 3 | 0.6471 | 0.0361 | 0.54798 | 0.54764 | 0.54956 | 0.52846 | 0.53983 | 0.53990 | 0.51064 | 1.0 | +0.29634 |
| `orion_hct116` +kolf2 | 3 | 0.6275 | 0.1323 | 0.56498 | **0.56644** | 0.56155 | 0.49777 | 0.55655 | 0.57260 | 0.51862 | 1.0 | +0.28856 |
| `orion_hek293t` +kolf2 | 3 | 0.6192 | 0.1437 | 0.56936 | 0.56764 | 0.56881 | 0.49215 | 0.55985 | 0.56613 | 0.52387 | 1.0 | +0.28073 |

Restricted to the **disagreement set** — the only place any rule can differ:

| fold | `sign(mean)` | majority | z-vote | donor 0 | donor 1 |
|---|---:|---:|---:|---:|---:|
| `k562_gwps` | 0.51543 | 0.51543 | **0.51880** | 0.50011 | 0.49989 |
| `orion_hct116` | 0.51895 | 0.51895 | **0.52253** | 0.47887 | 0.52113 |
| `orion_hek293t` | 0.53438 | 0.53438 | **0.54033** | 0.49024 | 0.50976 |
| **CANCER mean** | **0.52292** | **0.52292** | **0.52722** | 0.48974 | 0.51026 |
| `kolf2` (arity 3) | 0.51039 | **0.50903** | 0.50885 | 0.50982 | 0.49781 |

**Six readings, and together they are the result.**

1. **Contested mass is large and its sign content is a coin flip.** 36.27% of the pool's |t|
   weight is contested and `sign(mean)` gets **52.29%** of it right. On the **unanimous** set —
   43.49% of the mass, where every rule including the oracle returns the same answer — the
   incumbent scores **0.61574**. *All* of the recoverable deficit lives in the disagreement set,
   and that is where it is a coin flip.
2. **At production's PLOLO arity the majority does not exist.** Two donors ⇒ every contested pair
   is a tie: `mass_tie == mass_disagree` exactly on 3/3. This is not a small-sample accident, it
   is arithmetic, and it means **`sign(mean)` and the majority are the same estimator on the three
   cancer folds** — the mandate's premise (`+5, −1, −1`) needs three covering donors and PLOLO
   only ever has two on the folds that count.

   **The literal statistic the mandate asked for — the fraction of scored (target, gene) pairs
   where the mean and the majority actually DISAGREE — is therefore exactly `0.000%` on all three
   cancer folds at production arity.** Where the majority is a rule it is: **5.84%** on `kolf2`
   (295,551 of 5,059,509 pairs) and **5.52 / 5.54 / 15.34%** on the cancer folds with a third
   voter added. The confidence-weighted vote, which *is* well defined at arity 2, disagrees with
   the mean on **1.01 / 1.08 / 3.13%** of pairs — 269,198 coordinates across the three folds, of
   which 162,923 are actual sign changes once the whitened final's own sign is accounted for.
3. **At arity 3 the majority is measurable and it loses.** −0.00095 on `kolf2` (production's own
   three donors), −0.00020 averaged over the three cancer folds with a third voter added, positive
   on 1/3. Nothing here supports replacing the mean's sign with the majority's.
4. **The confidence-weighted vote is the only rule that beats the mean, by +0.00145, 3/3.** The
   weight is **per-(donor, gene)**: `z_s(i,j) = |e_s(i,j)| / sd_s(j)`, the donor's own response
   standardised by its own across-target spread for that gene. It is the `|t|`-analogue the
   mandate names — a t-like statistic for the donor's own effect, with the across-target sd as
   the only available standard error on a pseudobulk table.
5. **Why not 189's `theta`, and why not a cross-donor reliability.** A weighted sign vote is
   `sign(Σ_s w_s(i,j)·sign(e_s(i,j)))`. **Any weight that does not depend on `s` factors out of the
   sum exactly**, at every arity. `theta_g` is per-GENE, so it cancels identically and cannot be a
   vote weight at all. The same kills a symmetric leave-one-in cross-donor reliability at arity 2,
   where both donors necessarily receive the same number (`common232.rho_weights` documents it).
   **Per-donor coverage** is the third weight the mandate names: a per-donor *scalar* weight makes
   the arity-2 vote degenerate to "always trust donor k", and those rows are in the table —
   **0.48974 / 0.51026 on the disagreement set**, i.e. a coin flip, so that family is dead too.
6. **The ceiling is +0.17325 and nothing implementable comes near it.** The z-vote takes 0.84% of
   it at the consensus level and **0.52%** of it as delivered.

### 2b. Phase 1b — how much of a consensus sign even reaches the scorer

|t|-weighted agreement on the same pool, and the accuracy at each stage:

| fold | `sign(pre)` → whitened | → final | → delivered | acc `pre` | acc whitened | acc final | acc delivered |
|---|---:|---:|---:|---:|---:|---:|---:|
| `k562_gwps` | 0.8565 | 0.8189 | 0.7826 | 0.54813 | 0.54578 | 0.55701 | 0.57427 |
| `orion_hct116` | 0.8385 | 0.8104 | 0.8104 | 0.57163 | 0.55468 | 0.56396 | 0.56382 |
| `orion_hek293t` | 0.8547 | 0.8259 | 0.8245 | 0.57278 | 0.55097 | 0.55315 | 0.55315 |
| `kolf2` | 0.8634 | 0.8204 | 0.8204 | 0.50820 | 0.51012 | 0.51767 | 0.51744 |

*(`acc pre` here is `sign(pre)`, which is defined on the 0–5 fallback targets per fold where no
donor votes; §2's `sign(mean)` column returns 0 there and is scored as wrong. That is the whole
difference between, e.g., 0.54813 and 0.53775 on `k562_gwps`.)*

**The whitening operator and the two residual sources rewrite 17–22% of the consensus's signs by
|t| weight, and the renderer adds almost nothing on top.** So an upstream sign fix is diluted by a
factor ~0.8 before it is scored — which is why the arm in §3 is wired at the **final** effect
matrix rather than at the consensus, and why the delivered Δ`q_w` (+0.00090) is smaller than the
consensus-level one (+0.00145).

---

## 3. Phase 2 — the arm, and the proof that it is a SIGN arm (`arm232.py`)

Every arm is
```
final_arm            = final_prod                       everywhere except a mask M
final_arm[M]         = |final_prod[M]| * s_target[M]     s_target in {-1, +1}
```
`M` requires **≥ 2 covering donors**, a non-zero rule sign, **disagreement with `sign(mean)`**, and
a normalised vote margin ≥ `tau` (the dose). A coordinate whose target sign is 0 is left alone —
imposing a zero sign would zero a magnitude, which would make this a magnitude arm.

**Magnitude invariance is asserted elementwise at runtime** (`common232.impose_sign` raises if
`np.abs(out) != np.abs(final)` anywhere) and it held on **every one of the 282 rows**. Three
independent corroborations:

| check | result |
|---|---|
| `|effect|` elementwise vs incumbent | **bit-identical**, asserted on every arm |
| delivered median row norm ratio | max **0.11%** over all arms/folds — **no `@adel` α re-solve is needed**, unlike 228's +7.38% |
| `jac`, all depths 89/200/400 | **max |Δ| = 0.000e+00** — see §6 |

**And the mask touches exactly what it claims to.** The per-gene panel-total digest differs on
**777 / 807**, **777 / 870**, **777 / 882** blocks on the three cancer folds — `777 = 3 × 259`, and
**259 is exactly the number of targets with ≥ 2 covering donors** on each fold. On `kolf2` it is
**804 = 3 × 268**, and 268 is that fold's ≥2-coverage count (28 + 240). The arm provably moves
nothing on a singleton-covered target.

### Per-fold, the arm at its best dose (`vote_z@tau = 0`, the whole contested set)

| fold | sign flips | Δ`cos_agg` rel | Δ`pds` | Δ`nmae` (wald89) | Δ`q_w` |
|---|---:|---:|---:|---:|---:|
| `k562_gwps` | 94,667 | −0.018% | +0.00137 | +0.000085 | +0.00031 |
| `orion_hct116` | 35,240 | **+0.290%** | **+0.00264** | +0.000047 | +0.00088 |
| `orion_hek293t` | 33,016 | **+0.799%** | +0.00082 | **−0.000678** | +0.00149 |
| **CANCER mean** | 54,308 | **+0.357%** | **+0.00161** (3/3 positive) | **−0.000182** (1/3 improved) | **+0.00090** |
| `kolf2` *(apart)* | 110,152 | **−2.114%** | **−0.00216** | **+0.000617** | **−0.00168** |

The dose sweep is monotone away from `tau = 0`: every quantity shrinks toward zero as the
threshold rises (`vote_z@0.1` +0.095% rel cosine, `@0.3` −0.039%, `@0.7` −0.003%). **There is no
dose better than taking the whole contested set**, and the whole contested set is what §3 prices.

---

## 4. Phase 3 — the six-member price, led by `nmae` (`price232.py`)

Three cancer folds, shipped `h1union` stack, `wald89`, 228's routes unchanged. `k_cos` carried as
the **range** {1.000, 1.384} per §8.1c — it is not identified and no PLOLO→board scalar is quoted.

| arm | flips | Δ`q_w` | **Δ`nmae`** | `e = ΔM/ΔU` | rel cos % | Δ`pds` | pds+ |
|---|---:|---:|---:|---:|---:|---:|---:|
| **`vote_z@0`** (production, arity 2) | 162,923 | **+0.00090** | **−0.00018** | 6.80 | +0.357 | +0.00161 | **3/3** |
| `ctrl_mean@0` (CONTROL, same coordinates) | 106,274 | −0.00045 | +0.00006 | 14.8 | +0.185 | +0.00004 | 2/3 |
| `vote_z@0` (production + kolf2, arity 3) | 344,929 | −0.00028 | +0.00010 | −0.57 | −0.763 | +0.00050 | 2/3 |
| `ctrl_mean_maj@0` (control for the majority) | 538,281 | −0.00148 | +0.00056 | −27.6 | −0.910 | −0.00072 | 1/3 |
| **`vote_maj@0`** (MAJORITY, arity 3) | 825,393 | **−0.00258** | **+0.00135** | −14.2 | **−4.715** | −0.00205 | 2/3 |
| `null_rand@0` (falsifier 1) | 162,923 | −0.00178 | +0.00166 | 27.4 | −3.200 | −0.00219 | 0/3 |
| `null_wrong@0` (falsifier 2) | 124,312 | −0.00079 | +0.00071 | 8.18 | −2.543 | −0.00280 | 0/3 |
| `null_rand_maj@0` | 825,393 | −0.00744 | +0.00391 | 4.23 | −9.912 | −0.00971 | 0/3 |
| `null_wrong_maj@0` | 654,086 | −0.00821 | +0.00441 | 11.5 | −12.638 | −0.01322 | 0/3 |
| `vote_maj@*` (production, arity 2) | **0** | — | — | — | — | — | **EMPTY MASK at every dose** |

**Overall points. Gate +0.0019.**

| arm | `o_nmae` | `o_pds` | `o_mse` k1 | `o_mse` k1.384 | `o_fid` | `o_jac` | **TOTAL k1** | **TOTAL k1.384** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **`vote_z@0`** | **+0.00005** | +0.00060 | +0.00013 | +0.00018 | +0.00014 | **+0.00000** | **+0.00092** | **+0.00097** |
| `ctrl_mean@0` | −0.00002 | +0.00002 | +0.00007 | +0.00009 | +0.00006 | +0.00000 | +0.00013 | +0.00015 |
| `vote_z@0` (+kolf2) | −0.00003 | +0.00019 | −0.00028 | −0.00038 | +0.00000 | +0.00000 | −0.00011 | −0.00022 |
| `ctrl_mean_maj@0` | −0.00016 | −0.00027 | −0.00033 | −0.00046 | −0.00011 | +0.00000 | −0.00087 | −0.00099 |
| **`vote_maj@0`** | −0.00038 | −0.00076 | −0.00171 | −0.00236 | −0.00025 | +0.00000 | −0.00309 | **−0.00375** |
| `null_rand@0` | −0.00046 | −0.00082 | −0.00116 | −0.00160 | −0.00018 | +0.00000 | −0.00262 | −0.00306 |
| `null_wrong@0` | −0.00020 | −0.00104 | −0.00092 | −0.00127 | −0.00015 | +0.00000 | −0.00231 | −0.00266 |

`o_reach` (indicative, ordering half only) is **−0.000021** for the arm — negligible, and the arm
is profile-class so 169's toll is 0.

**Leading with `nmae`, as §8.1d requires.** The arm's Δ`nmae` is **−0.000182** — the **right**
sign, the first time in this programme that a cosine-positive arm has not been punished by `nmae`.
It is also **tiny and not consistent** (1/3 folds improved), and `o_nmae = +0.00005` is 5% of the
composite. §8.1d's wall is not what stops this arm; **size** is.

**The best case over every free choice, taken jointly** — both `k_cos` values, all three pools
(`wald89`/`lfc89`/`h89`), both `jac` transfers, and the dose chosen on the evaluation folds
themselves: **+0.00089 … +0.00097**. The maximum reachable anywhere is **+0.00097 against a
+0.0019 gate.** The arm does not clear even when every unresolved choice is settled in its favour.

---

## 5. §8.1d, tested and CONFIRMED — and the exchange rate that comes with it

The mandate asked for a falsifiable test of §8.1d's prediction that a sign arm should move `nmae`
the right way. The vote is too small to test it, so the same wiring was driven by an **oracle**
sign at seven doses and an **antioracle** sign at seven more — a controlled two-sided dose response
in |t|-weighted sign accuracy at **bit-identical magnitude**. This is the first isolated
measurement of the sign channel in the programme.

| dose | Δ`q_w` | **Δ`nmae`** | `e = ΔM/ΔU` | rel cos % | Δ`pds` | **TOTAL k1.384** |
|---|---:|---:|---:|---:|---:|---:|
| `oracle@1.0` | +0.41480 | **−0.08927** | **36.0** | +44.79 | +0.14188 | **+0.10414** |
| `oracle@0.5` | +0.20647 | −0.04444 | 38.0 | +22.37 | +0.08297 | +0.05649 |
| `oracle@0.25` | +0.10455 | −0.02237 | 32.2 | +11.21 | +0.04457 | +0.02945 |
| **`oracle@0.1`** | **+0.04129** | **−0.00889** | **40.5** | **+4.36** | +0.01784 | **+0.01170** |
| `oracle@0.05` | +0.02040 | −0.00451 | 34.1 | +2.23 | +0.00888 | +0.00587 |
| `oracle@0.02` | +0.01012 | −0.00215 | 29.3 | +1.06 | +0.00431 | +0.00284 |
| `antioracle@0.1` | −0.05386 | +0.01489 | 5.4 | −11.49 | −0.03173 | −0.02227 |
| `antioracle@1.0` | −0.53986 | +0.15187 | 7.6 | −84.36 | −0.45874 | −0.26107 |

**Three things, and they are the reason this experiment was worth running even though the arm fails.**

1. **`nmae` responds to signs with the right sign, monotonically, over 14 doses.**
   `Δnmae = −0.002573` per +0.01 `q_w`, **r = −0.995**. `e = ΔM/ΔU` sits at **29–40** against
   break-even **0.5** — i.e. the matched term `M` collects essentially all of it, which is exactly
   what §8.1d says a sign-matched `min(|pred|,|truth|)` must do when the signs are what improved.
2. **The head-to-head against 228, at matched cosine.** `oracle@0.1` buys **+4.36%** relative
   aggregate cosine; 228's `gain_reg@0.25` bought **+4.54%**. Same channel, same instrument, same
   folds:

   | | 228 `gain_reg@0.25` (MAGNITUDE) | 232 `oracle@0.1` (SIGN) |
   |---|---:|---:|
   | relative aggregate cosine | +4.54% | +4.36% |
   | Δ|t|-weighted sign accuracy | +0.00027 | **+0.04129** |
   | **Δ`nmae`** | **+0.0373** | **−0.00889** |
   | `e = ΔM/ΔU` | **0.072** | **40.5** |
   | `o_nmae` | −0.01043 | **+0.00249** |
   | **six-member composite** | **−0.00702** | **+0.01170** |

   **A 0.0187-point swing at identical cosine, decided only by whether the cosine was bought with
   signs or with magnitudes.** §8.1d is not an inference from the board any more; it is a
   measured, causal, two-sided contrast inside one instrument.
3. **The exchange rate.** Regressing each member on Δ`q_w` across the 14 oracle rows:

   | member | per **+0.01** `q_w` | r |
   |---|---:|---:|
   | `o_pds` | +0.002434 | 0.956 |
   | `o_mse` (k_cos 1.384) | +0.000740 | 0.982 |
   | **`o_nmae`** | **+0.000719** | **0.995** |
   | `o_fid` | +0.000097 | 1.000 |
   | `o_jac` | **+0.000000** | — (exactly pinned, §6) |
   | **composite (k 1.384)** | **+0.003990** | **0.975** |
   | composite (k 1.000) | +0.003785 | — |
   | raw Δ`nmae` | −0.002573 | −0.995 |
   | relative `cos_agg` | +1.477 % | 0.982 |

   **So rule 13's +0.0019 gate needs Δ`q_w` = +0.00476 (k 1.384) to +0.00502 (k 1.000).**

### What that says about the channel — the number to carry forward

| quantity | value |
|---|---:|
| ceiling on Δ`q_w` for **any** donor-sign-combination rule (consensus level) | **+0.17325** |
| the same, converted at the measured exchange rate | **≈ +0.069 overall** |
| Δ`q_w` needed to clear the +0.0019 gate | **+0.0048 … +0.0050** |
| **fraction of the ceiling a rule must capture to clear the gate** | **2.75%** |
| fraction captured by the best implementable rule (delivered) | **0.52%** |
| shortfall | **5.3×** |

**This is not a channel closed by arithmetic** — unlike §8.1c's ordering channel, where the floor
was 3.31× the entire prize. Here the prize is **36× the requirement** and the failure is
efficiency: `sign(mean)` is at 0.523 on the contested set and no function of the donors' own
numbers gets past 0.527. **The information needed to break a donor disagreement is not in the
donors.** Anything that would move this has to be a new *source*, not a new *rule* — which points
back at 218/225's acquisition programme, not at the mixing law.

---

## 6. Two structural facts, both reusable

**(a) `jac` is EXACTLY invariant under any magnitude-preserving sign arm.** Measured
**max |Δ`jac`| = 0.000e+00** across every arm, every fold and all three `n_real` depths
(89 / 200 / 400). The reason is a property of the shipped emitter, not of this arm: `rc_emit.CallSet`
picks its 400-gene head by `|effect|` inside the scorer's 5-CPM universe and its 2,200-gene block by
the **control cells'** mean per-cell CPM. Both are invariant to sign, so the entire 2,600-gene call
set is bit-identical and `n_pred` stays pinned at 2,600. `fid` is *not* invariant — it reads the
**sign** of the delivered lfc on that call set — and it moves by up to **0.01119** across the arms
here. **For the whole sign-arm class, `jac` may be carried at exactly zero and `fid` must not be.**

**(b) A per-gene weight cannot be a vote weight.** `sign(Σ_s w_s(i,j)·sign(e_s(i,j)))` is unchanged
by any weight independent of `s`. That rules out 189's `theta_g` outright, and it rules out any
symmetric cross-donor reliability at arity 2. Recorded so the next agent does not spend a run
weighting a vote by `theta`.

---

## 7. The five mandatory falsifiers

| # | falsifier | result | verdict |
|---|---|---|---|
| **1** | **Sign-shuffle null** — the same *number* of sign changes per target, drawn uniformly from coordinates with ≥ 2 covering donors, seed `20260904` | count-matched at 162,923 flips: Δ`q_w` **−0.00178**, rel cos **−3.20%**, Δ`pds` −0.00219, composite **−0.00306** against the arm's **+0.00097**. 0/3 folds positive on `pds` | **PASSES — the arm is not the flip count** |
| **2** | **Wrong-recipient null** — target `π(i)`'s mask and vote signs applied to target `i` under a derangement (184 died here: arm +0.0665, control +0.0667, null +0.0663) | composite **−0.00266**, rel cos **−2.54%**, Δ`q_w` −0.00079, 0/3 on `pds`. At arity 3 the majority's own wrong-recipient null is **−0.01334** against the majority's −0.00375 | **PASSES — no 184 signature** |
| **3** | **Leave-top-3-out retention, cosine and `nmae` separately** (222 §9: 14.6% vs 72.7%) | see table below | **the channels differ, and the arm fails significance on 2/3 cancer folds** |
| **4** | **Magnitude invariance** | `|effect|` **bit-identical**, asserted elementwise on all 282 rows; delivered median row-norm ratio within **0.11%**; `jac` **exactly** unchanged; digest differs on exactly `3 × (targets with ncov ≥ 2)` blocks | **PASSES — this is a sign arm, not a magnitude arm** |
| **5** | **Coverage honesty** | full census below; **singletons are left alone by construction** and the digest proves it | **PASSES** |

### Falsifier 3 — retention and the paired-target bootstrap (`retain232.py`, 4,000 draws)

| fold | Δ`cos_agg` | rel % | **ret. cos** | P(Δcos>0) | Δ`pds` | **ret. `pds`** | P(Δpds>0) | Δ`nmae` | **ret. `nmae`** | P(Δnmae<0) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `k562_gwps` | −0.000007 | −0.02 | n/m | 0.481 | +0.001373 | 0.443 | 0.895 | +0.000085 | n/m | 0.356 |
| `orion_hct116` | +0.000091 | +0.29 | **0.458** | 0.779 | +0.002645 | 0.357 | 0.931 | +0.000047 | n/m | 0.283 |
| `orion_hek293t` | +0.000194 | +0.80 | **0.709** | **0.978** | +0.000820 | 0.052 | 0.747 | −0.000678 | 0.465 | **0.982** |
| `kolf2` *(apart)* | −0.000298 | −2.11 | — | 0.000 | −0.002159 | — | 0.164 | +0.000617 | — | 0.001 |
| `ctrl_mean@0` — `orion_hct116` | +0.000141 | +0.45 | 0.533 | 0.946 | +0.000772 | −1.016 | 0.686 | −0.000056 | 0.011 | 0.918 |

**The channels do differ, exactly as 222 warned.** Cosine retention 0.458–0.709 (distributed);
`pds` retention **0.052–0.443** (concentrated — on `orion_hek293t` the three best targets carry
95% of the whole Δ`pds`, and `pds` is the arm's largest positive term). And the bootstrap is the
plainest statement of all: **only 1 of 3 cancer folds clears `P(Δcos>0) > 0.95`, none clears it on
`pds`, and `kolf2` is negative at `P(Δcos>0) = 0.000`.** `ctrl_mean` — the incumbent's own sign on
the identical coordinate set — beats the arm on cosine on one of the three folds (+0.45% vs +0.29%).
The arm wins on the fold-mean (+0.357% vs +0.185%, and Δ`q_w` +0.00090 vs −0.00045) but it is not
a robust win.

### Falsifier 5 — coverage census (|t|-weighted share of the `wald89` pool, production voters)

| fold | targets | **no covering donor** | ncov 1 (singleton) | ncov 2 | ncov 3 |
|---|---:|---:|---:|---:|---:|
| `k562_gwps` | 269 | 0.0549 *(5 targets)* | 0.0177 | **0.9274** | — |
| `orion_hct116` | 290 | 0.0006 *(0 targets)* | **0.2521** | 0.7473 | — |
| `orion_hek293t` | 294 | 0.0068 *(1 target)* | **0.2749** | 0.7183 | — |
| `kolf2` | 273 | 0.0042 | 0.0187 | 0.3104 | **0.6667** |

**Up to 27.5% of the scored pool's |t| weight sits on a gene only ONE donor carries.** Those pairs
have no vote at any arity and the arm leaves them untouched — confirmed by the digest, which
differs on exactly `3 × 259` blocks on each cancer fold (`259` = that fold's ≥2-coverage target
count) and `3 × 268` on `kolf2`. The 5 / 0 / 1 / 0 targets per fold that **no** donor carries take 062's
fallback row (`EFFECT_SCALE · mean_response / |mean_response|`), whose sign is a property of the
donors' mean response profile and not of any vote; those rows are excluded from the vote by
`ncov = 0` and from the `sign(pre) == sign(Σ_s e_s)` identity check, which then holds on
**100.0000%** of covered pairs on all four folds.

---

## 8. Caveats, stated as loudly as the results

- **PLOLO cannot express production's arity on the folds that count.** Production mixes three
  donors; each cancer fold has two. The majority-vs-mean question at arity 3 is answerable only on
  `kolf2` (which is pluripotent, reported apart, and excluded from every headline) and on an
  augmented `+kolf2` voter set that is *not* production. Both say the majority loses, and they
  agree, but neither is the shipped configuration. **This is the one place where a reader should
  discount the result**, and it is a limitation of the harness, not of the measurement.
- **`kolf2`'s `b_ctrl` is flat.** Every gene reads 55.2 CPM, so `√cpm` is constant and the
  `wald89` and `lfc89` pools are **identical** on that fold — a third independent reason its
  numbers are reported apart.
- **The exchange rate is fitted on oracle doses, which are large.** Its smallest rung
  (`oracle@0.01`, Δ`q_w` +0.00050) is the same order as the arm, and it lands at +0.00022
  composite against the rate's prediction of +0.00020 — so the extrapolation to the arm's scale
  is checked, not assumed.
- **`k_cos` is not identified** (§8.1c: 1.000 on `h1union`, 1.384 on `cd4f034`, a 7.84× spread
  across four board pairs). Every `mse` figure is quoted as a range and **no PLOLO→board scalar
  conversion appears anywhere in this document.** All numbers are PLOLO units, three cancer folds,
  shipped `h1union` stack, `wald89` pool, a23 render for `pds`/cosine and a15 for `nmae`/`q_w`,
  tag `232`.

---

## 9. INBOX-1 — answered: **the whitening operator costs 1.37 points of |t|-weighted sign accuracy, worth −0.00547 overall, and no operator experiment has ever charged it** (`gamma232.py`)

INBOX-1 item 1 (run the `z` variant once, report `M`/`U`/`e` and the six members, don't sweep) is
answered above and was in fact swept: **`vote_z@0` = +0.00097 against a +0.0019 gate**, Δ`nmae`
−0.00018, `e` 6.80, negative on `kolf2`, 1/3 folds significant. **Null. §§3–4 and §7 are the
one-line answer.**

Item 2 is the real one and it reproduces. Three cancer folds, shipped stack, **only `gamma`
varies** — `alpha`, `lam`, the operator file, CD4 @ 0.34, H1 union @ 0.28, the panel zeroing and
the render are all the incumbent's. `wald89`. `gamma = 0` is exactly "no whitening": the operator's
eigenbasis is **square (17,942 × 17,942)**, so at `gamma = 0` the map is `evecs @ evecs.T = I` plus
a preserve rescale of 1.0.

| gamma | `q` delivered | Δ`q` vs γ=0 | `q` whitened | `cos_agg` | `pds` | `nmae` | Δ`nmae` vs γ=0 |
|---|---:|---:|---:|---:|---:|---:|---:|
| **0** (no whitening) | **0.57828** | — | 0.56488 | **0.040855** | 0.719471 | **0.97688** | — |
| 0.35 | 0.57015 | −0.00813 | 0.55632 | 0.034275 | 0.773309 | 0.98041 | +0.00354 |
| 0.50 | 0.56722 | −0.01106 | 0.55373 | 0.032535 | 0.783215 | 0.98186 | +0.00498 |
| **0.65 (PRODUCTION)** | **0.56470** | **−0.01358** | **0.55048** | 0.031299 | **0.787077** | 0.98208 | **+0.00520** |
| 0.80 | 0.56380 | −0.01448 | 0.54817 | 0.030447 | 0.786205 | 0.98267 | +0.00579 |

**The operator's sign cost, isolated (consensus → whitened, γ = 0.65):**
`q` **0.56418 → 0.55048**, **Δ`q` = −0.01370 (−2.43% relative)**. At 232's measured exchange rate
that is **Δ`nmae` +0.00353** and **−0.00547 overall points**. The directly measured Δ`nmae`
across γ is **+0.00520**, so the sign channel explains **68%** of the whole `nmae` cost of
whitening — the rest is magnitude.

**An independent cross-check falls out.** 232's measured `dnmae/dq = −0.2573` implies
`k = 0.1286` under 188's identity `nmae = 1 + k(1 − 2q)`. **187 measured the |t|-weighted
median per-gene `k` at 0.12** by a completely different route. Two independent measurements,
0.1286 vs 0.12.

**But the cost must be quoted next to what the operator buys, or it is misleading.** On the same
three folds, production's γ against γ = 0: Δ`pds` **+0.06761** → `o_pds` **+0.02520**. The sign
cost is **−0.00547** — **22% of the `pds` credit**, real, and never previously charged. The
operator remains strongly net positive on this instrument.

> ⚠️ **One row here does not agree with `AGENTS.md` §4.2 and I am flagging it rather than
> resolving it.** On PLOLO the **un-whitened** consensus has a *higher* aggregate cosine
> (**0.040855 vs 0.031299, +30.6%**) and a much lower `pds` (0.71947 vs 0.78708). Hubness does not
> explain it — `rr` is 0.0031 at γ=0 against 0.0026 at γ=0.65 on `orion_hct116`, and the macro
> cosine moves the same way as the aggregate (0.0436 → 0.0321). §4.2 calls whitening the largest
> measured gain, and that was established **on the board through `mse`, i.e. through cosine**.
> **PLOLO's cosine and the board's mse-implied cosine disagree about γ**, and that is a
> harness-validity question, not something this experiment can settle. 170 measured the vertex
> end-to-end only inside [0.35, 0.80] and never at 0. **Do not read the γ=0 row as "whitening
> hurts cosine on the board."** Route it; I am not re-opening γ, per the inbox and per 170/a19.

---

## 10. Verdict

**The donors' majority sign is not better than the sign of their mean — it is a structural tie at
production's PLOLO arity and strictly worse (−0.00375 overall) at the only arity where it is a
rule. The best confidence-weighted vote prices at +0.00097 against a +0.0019 gate, is negative on
`kolf2`, and clears the paired bootstrap on 1 of 3 folds. §8.1d's wall is confirmed and it is NOT
what stops this arm: a pure sign move at bit-identical magnitude improves `nmae` at −0.002573 per
+0.01 sign accuracy with `e = 29–40`, and at matched cosine it beats 228's magnitude arm by 0.0187
overall points. The sign channel is worth ≈ +0.069 overall at oracle and a rule needs only 2.75%
of it to clear the gate — but every rule computable from the donors' own numbers captures 0.52%,
because `sign(mean)` on a contested pair is already a coin flip and nothing in the donors breaks
the tie. NO-GO. No slot proposed. Nothing built. The last untested construction in the consensus
mixing law is closed; the remaining direction gain is an acquisition problem, not a mixing one.**

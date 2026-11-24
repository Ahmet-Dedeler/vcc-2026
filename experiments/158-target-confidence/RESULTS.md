# 158 — the per-target scale is real, predictable, and not worth having

**Date:** 2026-09-02 · two held Orion lines, exact-space pseudobulk truth, deterministic
(no emitter, no seeds, no `cell-eval2`). Instrument: experiment 156's frozen a23 delivered
effect and `work/truth-{hct116,hek293t}.npz`.

## Verdict

**NO-GO.** A build-time key that predicts the mse-optimal per-target scale exists and is
solid — production's own `donor_cosine_key`, restricted to the scored-and-gated coordinates,
correlates with `s_p*` at Spearman **+0.395 / +0.465** and delivers **+16.2% / +29.6%**
aggregate transfer cosine (strict cross-surface) at exactly zero `pds` cost, with a strictly positive cluster
bootstrap and a permutation null at `p <= 0.0005` on both surfaces. It is not worth shipping
for two independent reasons, either of which alone is disqualifying:

1. **The prize is 10-30x smaller than the brief's prior.** The brief priced the oracle at
   +0.010 to +0.015 overall. Priced against AGENTS.md §4's own `mse` anchors and the
   ratio-of-sums definition, the **deployable rule** is worth **+0.0003 to +0.0029** overall
   and the **full oracle** +0.0014 to +0.0106, at the generous end of every assumption
   (§4). At the board's implied cosine the clamp bites and the member stays at exactly 0.
2. **`nmae` pays more than `mse` collects, and the ORACLE pays the most of all.** The
   recommended rule costs **−0.0023 / −0.0080** overall on `nmae` (γ=1 / γ=2), which is
   larger than the whole `mse` gain. And the *oracle* per-target scale — the one that
   maximises the `mse` gain by construction — costs **−0.0245** overall on `nmae`,
   **−0.0055** even re-levelled at its own `nmae` optimum. **The two members disagree about
   the per-target axis itself**, so no better key rescues this (§5).

The per-target magnitude axis therefore joins AGENTS.md §12's closed list, for the same
mechanism experiment 097 and `docs/CURRENT_STATE.md` already recorded: `nmae` is the binding
member on the entire shrinkage family. This is AGENTS.md 11.1 rule 2 in its exact original
form — the metric the parameter serves is the one that cannot move.

---

## 1. Reproduction of the brief's four baseline numbers — **PASS, all four exact**

`oracle.py`, production's own `mass_consistent_effect` called through
`156/massfix.py:mass_fix`, effects gated at a15's shipped `min_cpm = 10.0`, 132 targets and
18,401 scored genes (all 18,533 minus the surface's own target coordinates, matching
`pds_cosine`'s panel exclusion).

| quantity | HCT116 brief / measured | HEK293T brief / measured |
| --- | --- | --- |
| global-optimal single scale, `sse` | 0.999326 / **0.999326** | 0.999803 / **0.999803** |
| ORACLE per-target scale, `sse` | 0.998299 / **0.998299** | 0.999252 / **0.999252** |
| implied aggregate cosine | 0.02596→0.04125 (+59%) / **identical** | 0.01403→0.02735 (+95%) / **identical** |
| ORACLE MAGNITUDE ONLY | +0.000011 / **+0.000011** | +0.000005 / **+0.000005** |
| fraction with `s_p* <= 0` | 0.189 / **0.189** | 0.273 / **0.273** |
| `corr(s_p*, ‖d_p‖)` | +0.130 / **+0.130** | +0.134 / **+0.134** |

Two decompositions the brief did not carry, both from the same run:

- **ORACLE SIGN ONLY** (flip the rows whose `s_p* < 0`, one global scale) is worth
  **+0.000173 / +0.000211** of `sse`, i.e. 17% / 38% of the full oracle gain. So the lever
  is mostly *relative confidence among the correctly-signed rows*, not the anti-correlated
  fifth.
- **The positive-constrained oracle** `max(s_p*, 0)` reaches `sse` **0.998320 / 0.999305**
  against the unconstrained 0.998299 / 0.999252 — it recovers **98% / 91%** of the oracle
  gain. That matters because only a strictly positive scale is `pds`-free, so the reachable
  ceiling is essentially the whole oracle. All `capture` fractions below are quoted against
  this constrained ceiling, which is the honest denominator.

### `pds` invariance — verified numerically, not asserted

| scaling applied to the delivered rows | HCT116 `pds` | Δ | HEK293T `pds` | Δ |
| --- | ---: | ---: | ---: | ---: |
| none (baseline) | 0.8019315290 | — | 0.7544529262 | — |
| random lognormal(0, 1.5) per row | 0.8019315290 | **+0.00e+00** | 0.7544529262 | **+0.00e+00** |
| oracle, clipped to `>= 1e-6` | 0.8019315290 | **+0.00e+00** | 0.7544529262 | **+0.00e+00** |
| oracle **raw** (19% / 27% of rows `<= 0`) | 0.8467499422 | +4.48e-02 | 0.7568817950 | +2.43e-03 |

Exactly zero to the last printed digit for any strictly positive scale, on both surfaces, and
for every arm in §3. The last row is the brief's own warning made concrete: a scale that
reaches zero or flips sign is **not** a pure scaling and does move `pds` (here it would help,
but it needs the truth).

**One caveat that changes the deployment story.** `--shrink-pert-key` multiplies the effect
**before** the render, so what production would actually deliver is
`render(mass_fix(w·h))`, not `w·render(mass_fix(h))`. Those differ by 0.6% of the row norm,
and the pre-render path moves `pds` by **+4.05e-04 (HCT116) / −5.78e-04 (HEK293T)** —
sign-inconsistent across surfaces, so worst-surface −0.00058. The `sse` is unaffected to
3e-6. **The "exactly zero `pds`" property belongs to a post-render row scale, which is not
what the shipped flag does.** (`verify.py` check 4.)

---

## 2. Every key tried, and its association with `s_p*`

All keys are functions of donors, the recipient control profile, our own predicted effect, or
frozen external tables. The two surfaces have **disjoint** donor sets
(`104/work/<surface>/manifest.json`, `baseline.sources`): HCT116 ← `k562_gwps` +
`orion_hek293t`; HEK293T ← `k562_gwps` + `orion_hct116`. Nothing below reads the held
surface's perturbation responses.

Spearman correlation with `s_p*` and with the per-target cosine. **No selection happens in
this table** — both surfaces are shown for every key.

| key | ρ(s\*) HCT | ρ(s\*) HEK | ρ(cos) HCT | ρ(cos) HEK | sign transfers? |
| --- | ---: | ---: | ---: | ---: | :---: |
| **`donor_cos_live`** (donor cosine on scored∩gated genes) | **+0.395** | **+0.465** | +0.384 | +0.468 | yes |
| `donor_cos` (production's `donor_cosine_key`, full axis) | +0.375 | +0.411 | +0.364 | +0.406 | yes |
| `donor_cos_shared` (donor cosine on jointly-carried genes) | +0.375 | +0.412 | +0.364 | +0.406 | yes |
| `agree_share` (Σ`donor_agreement` / Σ mean\|donor\|) | +0.227 | +0.308 | +0.196 | +0.282 | yes |
| `agree_norm` (‖`donor_agreement`‖₂) | +0.246 | +0.250 | +0.311 | +0.299 | yes |
| `agree_frac` (fraction of genes with nonzero agreement) | +0.203 | +0.298 | +0.201 | +0.299 | yes |
| `basal_log_cpm` (target's own basal CPM in the recipient) | +0.161 | +0.210 | +0.172 | +0.202 | yes |
| `row_norm` (‖d_p‖) | +0.167 | +0.152 | +0.249 | +0.210 | yes |
| `donor_norm_geo` (geometric mean donor effect norm) | +0.158 | +0.150 | +0.222 | +0.200 | yes |
| `pm_support` (PerturbMatch carries the target) | +0.192 | +0.041 | +0.198 | +0.044 | weak |
| `donor_norm_ratio` (min/max donor norm) | +0.076 | +0.107 | +0.021 | +0.063 | yes |
| `kd_worst` / `kd_mean` (own-gene knockdown depth) | −0.064 / −0.120 | −0.098 / −0.087 | — | — | yes, tiny |
| `kd_k562_gwps` | −0.148 | −0.058 | — | — | weak |
| **`score_agreement_mean`** (exp 104 crossfit credit) | **+0.405** | **−0.441** | +0.443 | −0.436 | **NO — flips** |
| `score_agreement_wmean` | +0.392 | −0.435 | +0.433 | −0.427 | **NO — flips** |
| `score_full_mean` | +0.167 | −0.166 | +0.165 | −0.169 | **NO — flips** |
| `score_full_wmean` | +0.211 | −0.135 | +0.214 | −0.138 | **NO — flips** |

**`score_agreement_mean` is the reason this experiment was run both directions.** On either
surface *alone* it is the strongest key in the table — |ρ| = 0.405 / 0.441, beating
`donor_cos_live`. Its sign is **opposite** on the two surfaces, and used cross-surface it
**loses** 0.000512 / 0.000133 of `sse`, the worst two entries in §3. A single-surface screen
would have selected it and shipped a regression. Experiment 104 was already a NO-GO on its
own gates; this is a second, independent reading of the same overfit.

Not available and therefore untested: **per-donor cell counts for the target** (`meta.json`
carries only min/median/max, not a per-perturbation table) and **number of covering donors**
(degenerate — `vote_histogram` is `{"2": 132}` on both surfaces, so every target is carried
by exactly two donors and `--shrink-pert-missing` is a no-op here; on the real 300-target
panel 6 rows are not two-donor-carried and that policy is **unmeasured**). `row_nnz` is
degenerate on HEK293T (one distinct value) and was dropped.

### What each key achieves when used, hyperparameters fitted on the OTHER surface

`z_p = max(rank01(sign·key), floor)^gamma`, production's `pert_shrink_vector` form; one
global scale re-optimised afterwards (free, `pds`-neutral). `sign`, `gamma`, `floor` chosen
by argmax on the fit surface; `d sse` positive = `mse` improves.

| key | fit → eval | sign, γ, floor | eval `sse` | Δ`sse` | Δcos | fit → eval | Δ`sse` | Δcos |
| --- | --- | --- | ---: | ---: | ---: | --- | ---: | ---: |
| **`donor_cos_live`** | HEK→HCT | +1, 2.0, 0.3 | 0.999090 | **+0.000236** | **+16.2%** | HCT→HEK | **+0.000134** | **+29.6%** |
| `donor_cos` | HEK→HCT | +1, 2.0, 0.3 | 0.999102 | +0.000224 | +15.4% | HCT→HEK | +0.000098 | +22.4% |
| `donor_cos_shared` | HEK→HCT | +1, 2.0, 0.3 | 0.999104 | +0.000222 | +15.3% | HCT→HEK | +0.000099 | +22.6% |
| `agree_share` | HEK→HCT | +1, 2.0, 0.5 | 0.999191 | +0.000135 | +9.5% | HCT→HEK | +0.000094 | +21.7% |
| `agree_norm` | HEK→HCT | +1, 0.5, 0.0 | 0.999263 | +0.000063 | +4.6% | HCT→HEK | +0.000025 | +6.3% |
| `basal_log_cpm` | HEK→HCT | +1, 0.5, 0.3 | 0.999292 | +0.000034 | +2.5% | HCT→HEK | +0.000015 | +3.7% |
| `agree_frac` | HEK→HCT | +1, 1.0, 0.1 | 0.999322 | +0.000004 | +0.3% | HCT→HEK | +0.000038 | +9.3% |
| `row_norm` | HEK→HCT | +1, 0.25, 0.0 | 0.999303 | +0.000023 | +1.7% | HCT→HEK | +0.000009 | +2.3% |
| `donor_norm_geo` | HEK→HCT | +1, 0.25, 0.0 | 0.999304 | +0.000022 | +1.6% | HCT→HEK | +0.000009 | +2.2% |
| `pm_support` | HEK→HCT | +1, 0.25, 0.1 | 0.999312 | +0.000014 | +1.0% | HCT→HEK | −0.000000 | −0.1% |
| `donor_norm_ratio` | HEK→HCT | +1, 0.5, 0.2 | 0.999324 | +0.000002 | +0.2% | HCT→HEK | +0.000008 | +2.1% |
| `kd_k562_gwps` | HEK→HCT | −1, 0.25, 0.5 | 0.999315 | +0.000011 | +0.8% | HCT→HEK | −0.000000 | −0.1% |
| `kd_mean` | HEK→HCT | −1, 0.25, 0.1 | 0.999315 | +0.000011 | +0.8% | HCT→HEK | −0.000004 | −1.0% |
| `kd_worst` | HEK→HCT | −1, 0.25, 0.1 | 0.999327 | −0.000001 | −0.1% | HCT→HEK | +0.000001 | +0.2% |
| `score_full_wmean` | HEK→HCT | −1, 1.0, 0.5 | 0.999441 | **−0.000115** | −8.9% | HCT→HEK | −0.000032 | −8.5% |
| `score_full_mean` | HEK→HCT | −1, 2.0, 0.5 | 0.999581 | **−0.000255** | −21.2% | HCT→HEK | −0.000036 | −9.6% |
| `score_agreement_wmean` | HEK→HCT | −1, 2.0, 0.3 | 0.999829 | **−0.000503** | −49.7% | HCT→HEK | −0.000117 | −36.3% |
| `score_agreement_mean` | HEK→HCT | −1, 2.0, 0.3 | 0.999838 | **−0.000512** | −51.0% | HCT→HEK | −0.000133 | −42.9% |

Cross-surface hyperparameter transfer is essentially lossless for the donor-cosine family:
`donor_cos_live`'s own-surface best is +0.000263 / +0.000150 against the +0.000236 /
+0.000134 the other surface's argmax delivers, so the free parameters cost 10% of the gain,
not the gain.

---

## 3. The best deployable rule, and the discrete version

The grid is a **broad plateau, not a spike**, which is the property that makes a two-parameter
key on 132 targets believable. Δ`sse` × 1e6, every cell positive on both surfaces:

`donor_cos_live`

| γ | f=0.0 hct/hek | f=0.1 | f=0.2 | f=0.3 | f=0.5 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.25 | 91.5 / 42.0 | 89.9 / 42.6 | 84.4 / 39.5 | 75.2 / 34.3 | 48.9 / 20.9 |
| 0.50 | 142.7 / 73.3 | 146.5 / 75.6 | 146.5 / 72.9 | 137.7 / 65.7 | 95.9 / 41.7 |
| 1.00 | 192.0 / 114.6 | 197.4 / 117.0 | 209.5 / 117.8 | 215.6 / 113.1 | 177.5 / 80.8 |
| **2.00** | 207.8 / 146.6 | 209.1 / 147.1 | 217.9 / 148.6 | **236.3 / 150.0** | 263.2 / 133.9 |
| 4.00 | 147.0 / 137.4 | 147.0 / 137.4 | 147.7 / 137.5 | 152.0 / 138.1 | 198.9 / 142.2 |

`donor_cos` (production's exact key) has the same shape at 85-95% of the magnitude:
γ=1 f=0.3 reads 205.0 / 93.2; γ=2 f=0.3 reads 223.6 / 110.8.

### Arms, with the constrained-oracle capture fraction

| arm | surface | `sse` | Δ`sse` | Δcos | capture | Δ`pds` |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| **production default** `donor_cos`, γ 1.0, floor 0.3 — **0 DOF** | HCT116 | 0.999121 | +0.000205 | +14.2% | 0.204 | **+0.0e+00** |
| | HEK293T | 0.999710 | +0.000093 | +21.4% | 0.187 | **+0.0e+00** |
| `donor_cos_live`, γ 1.0, floor 0.3 (0 DOF, key swapped) | HCT116 | 0.999111 | +0.000216 | +14.9% | 0.214 | +0.0e+00 |
| | HEK293T | 0.999690 | +0.000113 | +25.5% | 0.227 | +0.0e+00 |
| **`donor_cos_live`, γ 2.0, floor 0.3** | HCT116 | 0.999090 | **+0.000236** | **+16.2%** | 0.235 | +0.0e+00 |
| | HEK293T | 0.999653 | +0.000150 † | +32.7% † | 0.301 | +0.0e+00 |
| `donor_cos_live`, γ 2.0, floor 0.5 | HCT116 | 0.999063 | +0.000263 ‡ | +17.9% ‡ | 0.261 | +0.0e+00 |
| | HEK293T | 0.999669 | **+0.000134** | **+29.6%** | 0.269 | +0.0e+00 |
| `donor_cos_live`, γ 4.0, floor 0.3 | HCT116 | 0.999174 | +0.000152 | +10.7% | 0.151 | +0.0e+00 |
| | HEK293T | 0.999665 | +0.000138 | +30.5% | 0.277 | +0.0e+00 |
| DISCRETE: worst 50% → 0.50 | HCT116 | 0.999211 | +0.000115 | +8.2% | 0.114 | +0.0e+00 |
| | HEK293T | 0.999733 | +0.000070 | +16.5% | 0.141 | +0.0e+00 |
| DISCRETE: worst 50% → 0.15 | HCT116 | 0.999253 | +0.000073 | +5.3% | 0.073 | +0.0e+00 |
| | HEK293T | 0.999711 | +0.000092 | +21.2% | 0.185 | +0.0e+00 |
| DISCRETE: worst 25% → 0.15 | HCT116 | 0.999310 | +0.000016 | +1.2% | 0.016 | +0.0e+00 |
| | HEK293T | 0.999773 | +0.000030 | +7.4% | 0.061 | +0.0e+00 |

† and ‡ are **own-surface** optima and are not leakage-clean: HCT116's argmax over the whole
grid is `floor 0.5` and HEK293T's is `floor 0.3`, so the strict cross-surface pair is
**+0.000236 (HCT, fitted on HEK) / +0.000134 (HEK, fitted on HCT)**, bolded. Both floors are
reported because the γ=2 row is flat in `floor` (§3 grid: 207.8→263.2 on HCT, 146.6→133.9 on
HEK) and the choice moves nothing in §6. The bootstrap below was run at `floor 0.3` on both
surfaces, so HEK's interval there sits at its own optimum and is the more generous of the two.

**The cheap discrete version is NOT the whole win — it is a third to a half of it.** Every
`(q, floor)` cell of the quantile rule is beaten by the continuous rank weight, and the best
discrete cell captures 7-19% of the ceiling against 24-30%. The brief's hypothesis is
refuted: the ordering among the *confident* targets carries as much of the lever as the
identification of the bad fifth does, which is the same thing the ORACLE SIGN ONLY
decomposition in §1 says (the anti-correlated fifth is only 17-38% of the gain).

Even the best rule collects only **24% / 30%** of the reachable ceiling.

### Bootstrap and permutation null

Target-cluster bootstrap: 20 Ward clusters of our own unit-normalised delivered rows (a
build-time grouping — targets predicted to point the same way share whatever systematic donor
error they carry, so resampling individual targets would understate the variance), 2000
paired resamples, the global scale re-optimised **inside** each resample so its own sampling
error is included. Permutation null: 2000 shuffles of the key across targets, weight form and
both hyperparameters held fixed.

| arm | surface | observed Δ`sse` | bootstrap 95% | perm `p` | null median | null p95 |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| `donor_cos_live` γ2 f0.3 | HCT116 | +0.000236 | [+0.000015, +0.000503] | **0.0005** | −0.000266 | −0.000130 |
| | HEK293T | +0.000150 | [+0.000044, +0.000269] | **0.0005** | −0.000077 | −0.000025 |
| `donor_cos` γ1 f0.3 (0 DOF) | HCT116 | +0.000205 | [+0.000048, +0.000417] | **0.0005** | −0.000105 | −0.000009 |
| | HEK293T | +0.000093 | [+0.000031, +0.000170] | **0.0005** | −0.000031 | +0.000006 |
| DISCRETE worst 50%→0.15 | HCT116 | +0.000073 | **[−0.000070, +0.000281]** | 0.0005 | −0.000238 | −0.000110 |
| | HEK293T | +0.000092 | [+0.000013, +0.000179] | 0.0005 | −0.000069 | −0.000021 |

`p = 0.0005` is the floor at 2000 permutations, i.e. no shuffle reached the observed value on
either surface for any arm. Note the null MEDIAN is **negative** — a random per-target
reweighting actively costs `sse`, because the global scale can no longer be optimal for every
row. So the observed positive values are not a free consequence of adding a per-row degree of
freedom. The discrete arm's HCT116 interval **includes zero**, which is the third independent
reason to prefer the continuous form.

**The statistical case is therefore strong and it is not the reason for the NO-GO.**

---

## 4. What it is worth on `mse` — and why the brief's prior is 10-30x too high

`expr_mse_unbiased_capped_norm` is `sum_p N_p / sum_p D_p` (metric BRIEF eq. 187), so the
per-target scale is separable and `s_p*` is its exact per-row argmin — that part of the brief
is right. The pricing is where it goes wrong. Writing `f` for the relative aggregate-cosine
gain, `c` for the board's aggregate cosine and `kappa = sum‖delta_p‖² / sum D_p >= 1` for the
denominator inflation of my uncorrected `sse` (`D_p` subtracts the truth's own sampling
noise; I do not), the reachable improvement at the `sse`-optimal global amplitude is

    d_raw_mse = -kappa * c^2 * ((1+f)^2 - 1),    d_overall = -d_raw_mse / 0.953 / 6

Δ overall, over the plausible range of both unknowns:

| scenario | `f` | c=.072, κ=1 | c=.072, κ=1.5 | c=.12, κ=1 | c=.12, κ=1.5 |
| --- | ---: | ---: | ---: | ---: | ---: |
| recommended rule, worst surface | +0.162 | +0.00032 | +0.00048 | +0.00088 | +0.00132 |
| recommended rule, best surface | +0.327 * | +0.00069 | +0.00104 | +0.00192 | +0.00288 |
| FULL ORACLE, worst surface | +0.589 | +0.00138 | +0.00207 | +0.00384 | +0.00576 |
| FULL ORACLE, best surface | +0.950 | +0.00254 | +0.00381 | +0.00705 | +0.01058 |

**The brief priced the oracle at +0.010 to +0.015 overall. The most generous single cell in
this table for the oracle is +0.0106, and that cell requires the best surface, the highest
plausible board cosine and a 50% denominator correction simultaneously. The deployable rule's
most generous cell is +0.0029 and its honest worst-surface cell is +0.0003.**

\* `f = +0.327` is the own-surface-optimal-floor variant (§3 †); the strict cross-surface
best-surface value is `f = +0.296`, which moves the rightmost cell from +0.00288 to +0.00263
and changes nothing.

Two further constraints on the upside, both from AGENTS.md:

- **§4.1b clamp.** `mse` is clamped to `[0, 1]` and its baseline anchor is 0.989. An arm that
  improves `sse` but leaves raw above 0.989 scores exactly 0. Recomputed at the board-implied
  `c = 0.072`, the recommended rule takes the *uncorrected* min-raw from 0.994816 to 0.992998
  (worst surface) — still above the baseline, hence **scaled 0**. The board's positive `mse`
  readings come from the jackknife/`rho` credit, which a per-target scale does not touch, so
  the arm's delta lands one-for-one on the raw value and is the entire effect.
- **`debench.overall_delta` refuses to convert `mse` and `nmae`** deltas, because both are
  per-perturbation normalised and are not functions of the reported aggregates (094 §1). The
  conversions in this section and in §5 are hand-computed from §4's anchors and are
  **estimates**, labelled as such.

**Flagged per AGENTS.md 11.1 rule 9.** `fid`, `jac` and `reach`'s normalisation all carry
`n_real`, which is 68-92 on the repaired HCT116 table against a challenge regime in none of
our harnesses' range. Every absolute `fid` / `jac` / `reach` / `n_real` number in this
experiment is unreadable; only orderings and same-surface ratios are used. The absolute
cosines (0.014-0.031) are also 3-5x below the board-implied 0.072 because the harness truth
carries 150 cells per target against the challenge's 400 (156's caveat) — which is exactly
why §4 is written in terms of the relative `f` and not the absolute `sse`.

---

## 5. All six members, and the finding that closes the lever

| member | effect of a positive per-target scale | how established |
| --- | --- | --- |
| `pds` | **EXACTLY 0.000000e+00**, both surfaces, every arm | measured, §1 (and ±5.8e-04 on the *pre-render* production path — flagged) |
| `mse` | +0.0003 to +0.0029 overall, likely clamped to 0 | measured Δ`sse`, extrapolated §4 |
| `reach` (keyed on \|effect\|) | **exactly invariant** — within-target \|effect\| ordering identical on 132/132 targets on both surfaces, and `reach` is a per-target statistic | measured, `final.py` |
| sign accuracy / head purity | 0.5417 → 0.5424 / 0.5421, i.e. inside noise | measured, `price.py` |
| **`nmae`** | **−0.0023 (γ=1) / −0.0080 (γ=2) overall — a LOSS** | measured, `price.py`, n=1 surface |
| `fid`, `jac` | unmeasured: both carry `n_real` (rule 9) and both need the emitter to recompute significance | stated, not laundered |

`reach` and within-target gene ordering are experiment 157's territory and this arm provably
does not move them, which is the clean part of the result.

### The `nmae` cost is structural, and the oracle pays the most

| arm | `nmae` | Δ`nmae` | Δ overall (est.) |
| --- | ---: | ---: | ---: |
| a23 baseline | 0.9479 | — | — |
| `donor_cos_live` γ1 f0.3, mean-preserved | 0.9561 | +0.0083 | **−0.00229** |
| `donor_cos_live` γ2 f0.3, mean-preserved | 0.9766 | +0.0287 | **−0.00796** |
| **ORACLE `max(s_p*,0)`, mean-preserved** | **1.0363** | **+0.0885** | **−0.02452** |

Two falsifiers, both of which fail to rescue it:

1. **Re-level every shape at its own `nmae`-optimal global amplitude** (production already has
   a free `VCC_ALPHA`, so holding the level at the incumbent's is not the fair comparison).
   The baseline improves too — its own optimum is `m* = 1.30`, `nmae` 0.9448 — and the arms
   land at **−0.00290 (γ=1) / −0.00554 (γ=2)**. The cost is slightly *worse*, not better.
2. **Is the lever intrinsically opposed, or is our key just weak?** Re-levelled at its own
   `nmae` optimum, the ORACLE per-target scale still costs **−0.00553**. The scale that
   maximises the `mse` gain by construction is `nmae`-pessimal. **The two members disagree
   about the per-target axis itself, so no better key fixes this.**

That is the same mechanism `docs/CURRENT_STATE.md` records for the whole a13 shrinkage family
and experiment 097 measured on 40 perturbations through the real rank-sum: `nmae` scores
absolute fold-change error over **the reference's** DE genes, which are not the genes our
donors agree about, so damping by donor confidence degrades it monotonically. This experiment
adds that **the degradation survives the oracle**, which the earlier work could not test
because it had no per-target truth to hand.

---

## 6. GO / NO-GO

**NO-GO.** Net estimate for the best deployable rule, worst surface, all measured members:

    mse   +0.00032 (optimistic +0.00288, most likely exactly 0 under the clamp)
    nmae  -0.00229 (gamma 1) / -0.00796 (gamma 2)
    pds    0.00000 exactly
    reach  0.00000 exactly
    fid, jac  unmeasured (rule 9 + emitter)
    ----------------------------------------
    net   negative at every gamma under every assumption in the section-4 table

The per-target magnitude axis is closed as an `mse` lever. **Do not reopen it without a new
mechanism that decouples `nmae` from the per-target amplitude** — for example a rule that
changes each row's *direction* rather than its scale, which is a different object and is not
`pds`-free.

**If it is ever revisited** (say `nmae`'s weight changes, or an emission change absorbs the
absolute-lfc cost), the configuration to run is the zero-DOF one, which spends no fitting
budget and is already implemented:

```bash
uv run python scripts/make_submission_pm_cd4.py \
  --shrink-pert-key donor-cos \
  --shrink-pert-gamma 1.0 \
  --shrink-pert-floor 0.3 \
  --shrink-pert-missing damp
  # (gamma/floor/missing are already the frozen defaults; the key flag alone reproduces it)
```

Δ`sse` +0.000205 / +0.000093, Δcos +14.2% / +21.4%, bootstrap [+0.000048, +0.000417] /
[+0.000031, +0.000170], perm `p <= 0.0005`, Δ`pds` exactly 0 as a post-render scale. Three
things to know before running it: (a) the flag scales the effect **pre-render**, which moves
`pds` by ±5.8e-04 rather than 0; (b) `--shrink-pert-missing` is a no-op on both harnesses
(132/132 two-donor) but binds on 6 of the real 300 targets and is **unmeasured**; (c) the
`donor_cos_live` variant that beats it by 10-60% is **not implemented** in the builder —
`donor_cosine_key` computes the cosine on the full panel axis, and restricting it to
`ctrl_cpm > 10` ∩ scored genes would be a one-line change and a new variable.

---

## 7. Known-answer checks — all four pass (`verify.py`)

Per AGENTS.md 11.1 rule 10's corollary, a screen with no configuration whose answer is
already known cannot detect the class of bug where the wrong numbers look as reasonable as
the right ones.

1. `rule.weight` vs production `pert_shrink_vector`, 200 randomised trials with forced ties
   and forced missing rows: **max\|diff\| 0.000e+00**.
2. Streaming per-gene donor mean vs `Donor.decompose()` on real `k562_gwps` data:
   **rel 2.0e-06** (production accumulates in float32).
3. `donor_cos` vs production `donor_cosine_key` on a synthetic two-donor `MultiDonorModel`
   with overlapping-but-unequal gene axes and a target only one donor carries:
   **max\|diff\| 1.1e-16**, 4/6 carried both ways.
4. The production pre-render path scored alongside the idealised post-render scale (§1).

Donor processing matches the build exactly: `MultiDonorModel(..., gene_norm=False, eb=False)`,
so the donor vector is `delta[p] − mean_over_all_perts(delta)` on the panel axis with zeros
off the donor's own axis. `_renorm` cancels out of every cosine and is applied explicitly to
the geometric-mean keys. Peak RSS 4.1 GB (Orion's 2.7 GB `delta.npy` is memory-mapped and the
mean is accumulated in row chunks).

## Files

`oracle.py` (§1, reproduction + `pds` invariance) · `keys.py` (§2, the 20 candidate keys,
streaming donor loader) · `rule.py` (§2-3, weight forms, cross-surface fitting, bootstrap,
permutation) · `final.py` (§3, grid stability, arms, all-six) · `price.py` (§4-5, board
pricing, `nmae` cost, both falsifiers) · `verify.py` (§7) · `work/*.json`.

No submission was built, no `cell-eval2` was run, no existing file was modified.

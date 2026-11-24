# 170 — `WHITEN_GAMMA` re-opened as an `mse` lever. It is a LOCAL OPTIMUM. Terminal NO-GO.

**Headline.** Lowering `WHITEN_GAMMA` does not improve `mse`; it makes it **worse**. Measured
end to end with the shipped operator, the challenge panel's raw `mse` has a **minimum at the
shipped `gamma = 0.65`** — `0.8813` at 0.60, **`0.8738` at 0.65**, `0.8795` at 0.70 — and it
degrades in both directions: `0.9261` at 0.35 and `0.8928` at 0.80. It is a **V with its vertex
on the shipped setting** — slope `−0.174` below it, `+0.127` above it. On the side a19 actually
moved, the measured slope is **4.8x smaller** than the `+0.6075` its board point implies.
**Every gamma in `[0.35, 0.80]` is net negative on all six members**, from `−0.0016` at 0.55 to
`−0.0071` at 0.35.

**The mechanism is the `rho` cliff, and it was not in anyone's model of this parameter.**
Lowering gamma raises the predicted effect's mean-square `A` (+14% from 0.65 to 0.35) and with
it the submission's own across-perturbation budget `B`. a15 sits essentially **on** the cliff
(`rho` = `0.99 / 0.94 / 0.88` per context, mean `0.938`; 164 measured `0.9466`). Raising `B`
past `claim` pushes all three contexts **above** it, where `min(claim, B)` saturates at `claim`
and the extra amplitude stops being refunded. 096 §(b) already proved `mse` is minimised
exactly at the cliff — this experiment shows **gamma moves the submission along that same
curve, and 0.65 is where it sits.** `WHITEN_GAMMA = 0.65` was fitted for `pds` in experiment
017; it turns out to also sit at the `mse` optimum, by coincidence.

**And a19 is re-attributed.** Its direction is now confirmed — raising gamma *does* hurt `mse`
— but its magnitude is not: gamma explains exactly **21%** of a19's `−0.0959`, measured at the same
gamma a19 shipped. **The other 79% is the ridge `lam`.** The `+0.6075` figure must never again be quoted as a gamma slope.

> ⚠️ **This corrects an intermediate version of this document.** A first pass translated the
> sweep with the `A − B` common-mode term alone, on the reasoning that `A − B` is flat in gamma
> (it is: `0.00117` → `0.00134` over the whole range). That is valid **only inside the
> `rho < 1` branch**, and the measurement shows lowering gamma leaves that branch. The
> difference is the whole result: `+0.0021` becomes `−0.0039` at gamma 0.50. The `min()` has to
> be evaluated per context at every gamma, not linearised around the shipped point.

Method and audit trail: [`METHOD.md`](METHOD.md).

---

## 1. Gates and controls

The sweep rebuilds production's own staging — `build_v9_pre_whiten_effects` →
`apply_perturbmatch_k562_replacement` → `apply_v9_whitening` → `apply_frozen_cd4_residual` →
`profile_from_effect` → `emit_counts` — by **calling the shipped functions**, and scores all six
members on **h1-test** through the **production `union_pruned` operator** with production's
exact three donors (`k562_gwps` + `orion_hct116` + `orion_hek293t`). `data/donors/h1_2025`
covers 85/85 of h1-test's targets and is pure oracle leakage; it is not used.

| # | known-answer gate (rule 10) | result |
| --- | --- | --- |
| 1 | `baseline_pre` + `candidate_pre` sha256 vs a15's shipped audit | **bit-identical** (`7678bab9…`, `0601336f…`) |
| 2a | a15's post-whitening audit, gamma 0.65 / lam 0 — 8 scalars | **PASS**, worst rel `1.1e-15` |
| 2b | a19's post-whitening audit, gamma 0.80 / lam 0.007 — 6 scalars | **PASS**, rel `0.00e+00` on all six |
| 3 | `msekit` on the frozen 095 h1-test artifact | **`2.055101829973639`** vs cell-eval2's own `agg_results.csv` **`2.0551018299736388`**, rel `2.2e-16` |
| 4 | 096's `mech_claim_grid2.py` re-run | `claim` **`0.11860`** at `spread 0.25 / phi 0`, cap **`0.24653`** — all 12 rows reproduce |
| 5 | challenge-panel `claim` at a15's policy and gamma | **`0.24610`** vs 164's **`0.24581`** (+0.12%) |
| 6 | challenge-panel `rho` at a15's policy and gamma | **`0.938`** vs 164's **`0.9466`** (−0.9%) |
| 7 | per-context `Kc` | **`0.005084 / 0.005228 / 0.005919`** — 096 §(e) exactly |
| 8 | `X_c` solved from a15's board point | **`0.017766`** vs 096's v11/v12 solve `0.018846` at `alpha 1.0796`, i.e. `0.01714` rescaled to a15's `0.982` before CD4 |

**Gate 1 initially FAILED, and the cause was real.** `candidate_pre` did not reproduce until
`--chunk-rows` was set to the builder's default **128**: PerturbMatch's centring mean is
accumulated one chunk at a time, so a different chunk size moves the last bits of `mean` and
every effect downstream. Exactly the silent mismatch rule 10 exists to catch.

**One thing is deliberately not asserted.** a15's *post*-whitening `final_sha256` does not
reproduce on this machine (a19's does). a15/a19 were built on Linux/x86_64 and this is arm64;
`whitener.apply` is a pair of 17,942 × 17,942 GEMMs and BLAS reduction order is
architecture-dependent. The evidence that this is last-bit noise and not a different
computation: all 14 scalars read off that object agree to ≤ `1.1e-15`, and the pre-whitening
object — which involves no large GEMM — is bit-identical.

### Controls, each with its reading when the mechanism is absent

| control | mass fix ON | mass fix OFF | published reference |
| --- | ---: | ---: | --- |
| render residual vs the intended effect | **0.0000%** | **2.67 – 3.45%** | 156: 2.83–3.35% on the challenge frame |
| emitted-pseudobulk rel L2 (profile norm) | 0.001583 | 0.001581 | a21: 0.0014 correct / 0.0155 mis-ordered repair |
| row-wise cosine, delivered effect at gamma 0.35 vs 0.80 | **median 0.962, min 0.896** | 0.962 / 0.897 | 1.000 ⇒ operator not applied |

The render control is the non-vacuous one: it is an identity **only because the mass fix makes
it one**, and switching the mass fix off moves it to 3%. The third bounds the whole lever —
sweeping gamma across `[0.35, 0.80]` rotates the delivered effect by a median of **16°**.

### The re-seeded twin: the emission noise floor on this surface

Same gamma, seed 0 vs seed 1, paired per target:

| gamma | Δ`mse` | Δ`pds` | Δ`nmae` | Δ`reach` | Δ`fid` | Δ`jac` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.45 | **0.00000** | **0.000000** | −0.00002 | **−0.00622** | +0.00124 | −0.00032 |
| 0.55 | **0.00000** | **0.000000** | −0.00008 | **+0.00542** | −0.00033 | −0.00034 |
| 0.65 | **0.00000** | **0.000000** | +0.00001 | **−0.00958** | +0.00016 | +0.00005 |

**`mse` and `pds` are bit-identical across seeds at all three gammas** — the quasi emitter's
exact-total repair pins the pseudobulk, so their gamma deltas carry no noise allowance at all.
`nmae` is seed-free to `1e-5`. **`reach` moves by `−0.0062`, `+0.0054`, `−0.0096` raw between
two seeds of the SAME arm — up to `±0.0018` overall**, larger than the entire gamma effect on
that member. `fid` carries `±0.0007`, `jac` `±0.0002`. **The three DE members cannot resolve
this lever**, and every number they contribute below should be read with that band.

---

## 2. The six-member table on h1-test

a15's flags plus `VCC_MASS_FIX=minnorm` (a23's shipped geometry, which is what the prescribed
build invocation produces). `mse` is scorer-exact through 096's `claim` route; `pds` is
`cell_eval2`'s own `discrimination_score` on the vcc2026 profile; `nmae` is exact;
`reach`/`fid`/`jac` carry `n_real` (rule 9) and are read as deltas only.

| gamma | `mse` | `pds` | `nmae` | `reach` | `fid` | `jac` | `A` | `X` | `B` | `claim` | `rho` | agg cos |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.35 | 1.01012 | 0.87129 | **0.95603** | 0.18338 | **0.51027** | **0.06465** | 0.08505 | **0.008960** | 0.08357 | 0.06701 | 1.000 | 0.06528 |
| 0.45 | 0.99939 | 0.87549 | 0.95723 | 0.18148 | 0.50830 | 0.06427 | 0.08233 | 0.008422 | 0.08089 | 0.06701 | 1.000 | 0.06236 |
| 0.50 | 1.00571 | **0.87661** | 0.95769 | 0.17005 | 0.50890 | 0.06424 | 0.08298 | 0.008263 | 0.08154 | 0.06701 | 1.000 | 0.06094 |
| 0.55 | 1.00520 | **0.87661** | 0.95839 | 0.17275 | 0.50878 | 0.06437 | 0.08248 | 0.008052 | 0.08104 | 0.06701 | 1.000 | 0.05957 |
| 0.60 | 1.00853 | 0.87633 | 0.95930 | 0.16297 | 0.50889 | 0.06432 | 0.08263 | 0.007872 | 0.08119 | 0.06701 | 1.000 | 0.05819 |
| **0.65 (ships)** | 1.01597 | 0.87535 | 0.95988 | 0.17752 | 0.50816 | 0.06425 | 0.08349 | 0.007734 | 0.08204 | 0.06701 | 1.000 | 0.05687 |
| 0.70 | 1.00304 | 0.87353 | 0.96122 | 0.17711 | 0.50772 | 0.06423 | 0.08093 | 0.007441 | 0.07951 | 0.06701 | 1.000 | 0.05558 |
| 0.75 | 0.99604 | 0.87367 | 0.96236 | 0.18628 | 0.50777 | 0.06434 | 0.07935 | 0.007189 | 0.07796 | 0.06701 | 1.000 | 0.05423 |
| 0.80 | 0.96449 | 0.87311 | 0.96400 | 0.17968 | 0.50664 | 0.06397 | 0.07367 | 0.006763 | 0.07236 | 0.06701 | 1.000 | 0.05294 |

Five things fall straight out of it.

1. **`claim` is `0.06701` at every gamma, to five decimals** — and on the challenge panel it is
   `0.24601 … 0.24619` across all seven gammas measured, likewise flat. The emission policy is
   unchanged and `min(Ĉ_p, C_p)` saturates, so **gamma cannot move the emission channel at
   all**. Derived, not assumed, as the mandate required.
2. **`X` is the clean channel.** `d(log X)/d(gamma) = −0.576 ± 0.035`, R² **0.974** over nine
   points. Alignment falls by more than half per unit gamma, monotonically, no scatter.
3. **`A − B`, the predictions' common mode, is flat** (`0.00148 … 0.00131` here, `0.00134 …
   0.00117` on the challenge panel). It is *not* what governs `mse` — see the correction at the
   top and §4.
4. **`A` moves a lot, and it is the term that matters.** On the challenge panel it runs
   `0.2648 → 0.2315` from gamma 0.35 to 0.65, a **14% fall**. On h1-test the same quantity is
   masked: `apply_v9_whitening` preserves the **median** row norm on the operator columns while
   `A` is a weighted **mean of squares**, and on 85 rows that wobbles (R² 0.62). The 300-row
   challenge panel resolves it cleanly (R² 0.997).
5. **`pds` has an interior maximum at gamma 0.50–0.55**, not at 0.65 and not monotone. 097 §9h
   found `pds` monotone increasing in gamma; it scored the **intended** effect, and 156
   established that the board sees the **delivered** one. On the delivered, mass-fixed object
   the curve turns over below the shipped setting. The gain is small (`+0.00126` raw) but it is
   **seed-free**.

⚠️ **h1-test cannot price the `mse` member.** Its own raw `mse` sits at `0.96–1.02`, above the
board's clamp point `0.98786`, so `score_mse` is clamped to 0 there at every gamma. That is why
§4 exists.

---

## 3. Measured slopes, against a19's implied `+0.61`

| quantity | slope in gamma | s.e. | R² |
| --- | ---: | ---: | ---: |
| `log X` (h1-test) | **−0.57628** | 0.03544 | **0.974** |
| raw `nmae` (h1-test) | **+0.01723** | 0.00132 | **0.961** |
| raw `pds` (h1-test) | −0.00065 | 0.00478 | 0.003 |
| raw `A` (challenge, 9 pts) | **−0.12787** | — | 0.993 |
| **raw `mse` (CHALLENGE, measured, §4) — below 0.65** | **−0.17437** | — | V-shaped, vertex at 0.65 |
| **raw `mse` (CHALLENGE, measured, §4) — above 0.65** | **+0.12687** | — | the side a19 moved |
| a19's board-implied `d(raw mse)/d(gamma)` | **+0.6075** | — | one point, `lam`-confounded |

**Cross-check against a different pipeline.** 097 §9h measured `kappa` on the *intended*,
un-mass-fixed effect with a13's shrinkage and got `kappa(0.45)/kappa(0.65) = 1.0703` on
h1-test. This experiment measures `X(0.45)/X(0.65) = 1.0889` on the *delivered*, mass-fixed
effect. Same sign, same size, from an independent reconstruction — the alignment channel
replicates. 097's orion ratio is `1.0408`, **41% smaller than h1-test's**; reported as a
labelled ordering cross-check only (168d: never size on the Orion pair), it says the `X`
channel may be up to 40% weaker on somatic contexts, which makes the verdict *more* negative.

---

## 4. The challenge panel: `A`, `B`, `claim` and `rho` measured at every gamma

`A`, `B` and `claim` are properties of our own submission and the challenge's own control
cells, so all three are computable with no access to the answers (096 §(a)). `sumD = 0.229141`
is 096's board-solved constant — a property of the data, identical for everyone. `X_c` is
solved from a15's **own board point** at gamma 0.65 (one equation, one unknown) and carried
across gamma by the ratio measured on h1-test. The identity is applied **per context**, with
`min(claim_c, B_c)` evaluated at each gamma:

```
   MSE_c - 1 = ( A_c + Kc_c - min(claim_c, B_c) - 2 X_c ) / sumD
```

| gamma | `A` | `B` | `claim` | mean `rho` | ctx at `rho`=1 | `X/X₀` | **raw `mse`** | **Δoverall from `mse`** |
| ---: | ---: | ---: | ---: | ---: | :-: | ---: | ---: | ---: |
| 0.35 | 0.26483 | 0.26350 | 0.24601 | 1.0000 | 3/3 | 1.1586 | **0.92610** | **−0.00918** |
| 0.45 | 0.25669 | 0.25539 | 0.24612 | 0.9924 | 2/3 | 1.0889 | **0.90965** | **−0.00629** |
| 0.50 | 0.25096 | 0.24968 | 0.24615 | 0.9852 | 2/3 | 1.0684 | **0.89602** | **−0.00390** |
| 0.55 | 0.24567 | 0.24441 | 0.24608 | 0.9784 | 1/3 | 1.0411 | 0.88530 | −0.00202 |
| 0.60 | 0.24043 | 0.23919 | 0.24619 | 0.9640 | 1/3 | 1.0179 | 0.88133 | −0.00132 |
| **0.65 (ships)** | 0.23152 | 0.23031 | 0.24610 | **0.9381** | 0/3 | 1.0000 | **0.87379** | **0** |
| 0.70 | 0.22386 | 0.22269 | 0.24610 | 0.9071 | 0/3 | 0.9621 | 0.87953 | −0.00101 |
| 0.75 | 0.21814 | 0.21699 | 0.24614 | 0.8837 | 0/3 | 0.9296 | 0.88450 | −0.00188 |
| 0.80 | 0.20624 | 0.20514 | 0.24613 | 0.8355 | 0/3 | 0.8744 | 0.89282 | −0.00334 |

**Raw `mse` has a minimum at the shipped gamma.** It is worse at 0.60 (`0.8813`) and worse at
0.70 (`0.8795`), and the curve rises monotonically away from 0.65 in both directions. The a23 anchor (raw `0.872517`) gives the same overall deltas to `±0.00005`,
so the choice of baseline does not matter.

### The mechanism

096 §(b) established that `mse(s)` **decreases up to the `rho` cliff `s_c = sqrt(claim/B)` and
increases after it — the minimum sits exactly on the cliff**, because below the cliff the
metric refunds the whole differential part of the amplitude through #348's budget, and above it
the refund saturates at `claim`.

Gamma is a second parameter that moves the submission along that same curve. Lowering it
flattens the operator less, leaves more energy in the dominant directions, and raises `A` (and
with it `B`) by 14% over `[0.65, 0.35]`. `B` crosses `claim` between gamma **0.60 and 0.65** in the
first context, between 0.50 and 0.55 in the second, and by 0.35 in all three. **Above the cliff the extra amplitude is pure cost**,
and it swamps the `+16%` of alignment that lowering gamma buys.

**So `WHITEN_GAMMA = 0.65` sits on the `mse` cliff.** It was fitted against `pds_cosine` in
experiment 017, long before the `mse` mechanism was understood, and lands at the `mse` optimum
by coincidence. Worth recording precisely because it means the parameter is **doubly** pinned,
and because any future change to `alpha`, the emission dispersion or `keep_top` that moves `B`
relative to `claim` will move this optimum with it.

### What this does to the hypothesis

The experiment was opened on the reading that gamma 0.50 would give raw `mse ≈ 0.783` — better
than the leader's `0.7767` — worth `+0.016` overall from `mse` alone.

**Measured: gamma 0.50 gives raw `mse = 0.896`, worth `−0.0039` overall.** The hypothesis is
**refuted in sign as well as in magnitude**. The naive back-extrapolation was wrong by `0.113`
of raw `mse`, and it pointed the wrong way.

### The re-attribution of a19

a19 moved gamma `0.65 → 0.80` **and** `lam 0 → 0.007`, and lost `0.09591` scaled `mse`. The
table **confirms a19's direction** — raising gamma past 0.65 does degrade `mse` — while
refuting its magnitude: the **measured** gamma-only cost at 0.80 is `0.01903` raw = `0.02003`
scaled, i.e. **21% of a19's `0.09591`**. Nothing was extrapolated — 0.80 is a swept point.

> **The remaining 79% of a19's `mse` collapse is the RIDGE, not gamma.** They are visibly
> different interventions in a19's own meta: `lam` cut operator amplification `25.74 → 22.03`,
> while gamma `0.65 → 0.80` *raises* it `25.74 → 47.29` (measured here). Flooring eigenvalues
> at `7e-3` when the operator's smallest are `5.9e-12` is a far larger reshaping than moving
> the exponent, and a19 confounded them.
>
> **Standing consequence: `+0.6075` is not a gamma slope and must not be quoted as one.** It is
> a (gamma, lam) composite dominated by `lam`. This does not open a new arm — the only
> direction available from the shipped `lam = 0` is upward, and upward is what a19 measured at
> `−0.0959` — but it does mean the eigenvalue floor is by far the most `mse`-sensitive knob on
> the operator, and any future operator work should price it first.

---

## 5. All six members priced (rule 12), and the composed gate (rule 13)

| gamma | `mse` | `pds` | `nmae` | `reach` | `fid` | `jac` | **TOTAL** | verdict |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0.35 | −0.00918 | −0.00150 | +0.00107 | +0.00110 | +0.00119 | +0.00019 | **−0.00713** | VETO |
| 0.45 | −0.00629 | +0.00005 | +0.00073 | +0.00075 | +0.00008 | +0.00001 | **−0.00466** | VETO |
| 0.50 | −0.00390 | +0.00046 | +0.00061 | −0.00141 | +0.00042 | −0.00000 | **−0.00382** | VETO |
| 0.55 | −0.00202 | +0.00046 | +0.00041 | −0.00090 | +0.00035 | +0.00006 | **−0.00163** | VETO |
| 0.60 | −0.00132 | +0.00036 | +0.00016 | −0.00274 | +0.00041 | +0.00003 | **−0.00309** | VETO |
| **0.65** | 0 | 0 | 0 | 0 | 0 | 0 | **0** | **ships** |
| 0.70 | −0.00101 | −0.00067 | −0.00037 | −0.00008 | −0.00024 | −0.00001 | **−0.00237** | VETO |
| 0.75 | −0.00188 | −0.00062 | −0.00069 | +0.00165 | −0.00022 | +0.00004 | **−0.00170** | VETO |
| 0.80 | −0.00334 | −0.00083 | −0.00114 | +0.00041 | −0.00085 | −0.00013 | **−0.00588** | VETO |

Rule 13's floors are **`+0.0019`** (shape-only) and **`+0.0089`** (operator/amplitude, a19's
class). **No arm is even positive**, so the floors are not the binding constraint — the sign is.

| member | how it was priced | confidence |
| --- | --- | --- |
| `mse` | 096's `claim` route: `A`, `B`, `claim`, `Kc` measured per context at every gamma; `sumD` board-solved; `X_c` solved from a15's own board point and carried by the h1-test ratio | good — three independent known-answer hits (`claim`, `rho`, `X_c`) at the shipped point |
| `pds` | `cell_eval2.discrimination_score`, exact on h1-test, **seed-free** | good; interior max at 0.50–0.55, worth `+0.00046` |
| `nmae` | exact, `debench.Surface.nmae` against h1-test's own reference DE table, seed-free | good; R² 0.96 in gamma, monotone |
| `reach` | realised emitted p-order on h1-test — the only surface that reproduces a board `reach` point | **poor**: the re-seeded twin moves it by `±0.0018` overall, larger than the effect |
| `fid` | realised call set, precision leg | fair; twin noise `±0.0007`; `n_real` regime caveat (rule 9) |
| `jac` | realised call set | flat; twin noise `±0.0002` |

**No member is unmeasured.** That is the rule-12 requirement a19 failed.

---

## 6. Verdict — TERMINAL NO-GO

**Do not change `WHITEN_GAMMA`. No slot.** The shipped `0.65` is a local optimum of the
composed six-member score *and* of raw `mse` on its own. All eight other settings swept in
`[0.35, 0.80]` are negative, by `−0.0016` to `−0.0071` overall. There is no operating point.

Do not reopen without a **materially different mechanism** — specifically, without first moving
`B` relative to `claim`. The `mse` optimum in gamma is set by where the submission sits on the
`rho` cliff, so it is not a property of the operator alone: change `VCC_ALPHA`, the emission
dispersion policy or `keep_top`, and this optimum moves with them. A joint `(alpha, gamma)`
screen is the only version of the question still open, and 164 already measured a15 at 99.7% of
the dispersion cap, which leaves very little room on the other axis.

### 097 was right, and this settles why

097 §9h closed this route on a **two-surface** rule ("every setting below 0.65 gains on h1-test
and loses on orion"), and 168d has since discredited the Orion pair for exactly this kind of
large reshaping — so the *rule* 097 used no longer holds. **Its conclusion is nevertheless
correct, and for a reason it did not have:** 097 priced `mse` analytically through
`MSE_min = 1 + (Kc − 2√(claim·W))/sumD`, which **assumes the arm sits ON the rho cliff**. 164
has since measured that a15 sits just below it, and this experiment measures that lowering
gamma pushes the arm *above* it. The assumption that made 097's `mse` estimate positive
(`+0.0108` overall at gamma 0 on h1-test) is precisely the one that fails. **097's
`+0.0022 / +0.0027 / +0.0023` h1-test net at gamma `0.15 / 0.30 / 0.45` should be read as
`−0.005` to `−0.007`.**

"Production's 0.65 is at or beside the optimum — do not change it" was right for a partly wrong
reason. It is now right for the right one.

### What would change this

The lever is `X`, and `d(log X)/d(gamma) = −0.576` is the cleanest measurement in this
experiment: lowering gamma really does buy `+16%` of alignment over `[0.65, 0.35]`. It converts
to negative `mse` only because it must be bought with amplitude the metric will not refund
above the cliff. **A transform that raised `X` without raising `A` would convert the full
`+16%`** — worth about `+0.0055` overall on `mse` at our current `X_c = 0.0178`, and **≥ 1.28x
that at the leader's implied cosine** (096 §(b)). That is the same standing conclusion as
157 / 158 / 159 / 160 / 164 — the binding constraint is the direction, not the levers on top of
it — but it now comes with a specific, testable shape for what a useful lever looks like:
**norm-preserving, not norm-adding.**

---

## 7. What was NOT measured, and why

* **`lam` (the ridge) was not swept.** The mandate was gamma alone and a19's confound is the
  reason. §4 shows `lam` carries **79%** of a19's `mse` move, so it is the larger object — but
  the only direction available from `lam = 0` is upward, and a19 measured that at `−0.0959`.
  Reopening it needs a mechanism, not a sweep.
* **The Orion pair was not used for sizing** (168d's standing rule). It enters only through
  097's already-published `kappa` ratio, labelled as an ordering cross-check.
* **`reach` is reported but not resolved.** Its between-arm spread on h1-test is smaller than
  the re-seed noise measured on the same surface (§1). Two seeds is a floor estimate, not a
  distribution.
* **`fid` and `jac` carry `n_real`** (rule 9): h1-test's is 90 median / 766 mean against the
  challenge's ~240 (settled by a21). Only their orderings transfer.
* **`X_c` on the challenge panel cannot be measured, only solved.** It is board-anchored at
  gamma 0.65 and carried across gamma by the h1-test `X` ratio. 097's orion number says that
  ratio could be 41% smaller on somatic contexts, which would make every arm *more* negative.
* **`A` and `B` on the challenge panel are computed on the INTENDED profile**, not the emitted
  one. The check that this is adequate is gate 6: the resulting `rho` at the shipped point is
  `0.938` against 164's emitted-artifact measurement of `0.9466`, a 0.9% agreement. Emitting
  all 900 challenge rows per gamma would remove the approximation and costs ~9 hours.
* **A no-mass-fix (a15-geometry) sweep was not run end to end.** The controls establish the
  render residual under both settings, and the `mse` translation anchored on a15's and a23's
  board points agrees to `±0.00005` overall, so the baseline choice does not move the verdict.
* **`fine.py` (the 0.025-resolution shape) was written but not run** — the 9-point emitted
  sweep and the 7-point challenge translation already resolve the shape, and the machine was
  saturated by two sibling agents.
* **No submission was built, no slot was spent, nothing was written to `submissions/`.** The
  submission hold is in force.
* **`cell-eval2` was never run.** Every scored number comes from its primitives called on
  cached sufficient statistics (096's `msekit`), anchored to a real `cell-eval2` run by gate 3.
* Peak RSS: `stage.py` 6.9 GiB, `gate.py` 4.2 GiB, `h1sweep.py` 6.6 GiB, `chal.py` 3.1 GiB —
  all under the 10 GB ceiling.

---

## 8. If it were shipped anyway

It should not be. Recorded only so the arm is resolvable. The prescribed invocation differs
from **a23** in one field (`whiten_gamma`); against **a15** it would differ in two
(`mass_fix` *and* `whiten_gamma`), which violates rule 5 and rule 10 — a15 does not ship
`VCC_MASS_FIX`.

```bash
VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.10 VCC_MASS_FIX=minnorm VCC_WHITEN_GAMMA=<g> \
uv run python scripts/make_submission_pm_cd4.py \
  --min-cpm 10.0 --overdispersion 4.0 --tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 \
  --whitener-id union_pruned \
  --whitener data/donors/whitener_union_pruned.npz \
  --whitener-metadata data/donors/whitener_union_pruned.json \
  --output submissions/<name>.h5ad
```

Config diff against a23's meta must be **exactly 1 of 46 fields**, and
`meta['operator']['sha256']` must equal `e2f65ecf…`. Preregistered prediction for the least-bad
arm, `gamma = 0.55`, against a23's `0.1733597079`:

| member | predicted Δ scaled | falsifier |
| --- | ---: | --- |
| `mse` | **−0.01211** | fires outside `[−0.025, 0.000]` |
| `pds` | +0.00276 | fires if `pds` falls |
| `nmae` | +0.00248 | fires if `nmae` rises |
| `reach` | −0.00539 | unreadable — `±0.0108` of twin noise |
| `fid` | +0.00209 | rule 9 |
| `jac` | +0.00033 | rule 9 |
| **overall** | **−0.00163** | fires above `+0.000` or below `−0.006` |

---

## 9. Files

| file | what it does |
| --- | --- |
| `stage.py` | production-faithful pre-whitening staging for the 300-target panel and h1-test's 85; the bit-exact `candidate_pre` gate |
| `gate.py` | applies the shipped `union_pruned` operator at each gamma + the CD4 residual; the two-board-point known-answer gate |
| `h1kit.py` | the h1-test surface: render, emit, DE, and the exact `pds` / `mse` routes |
| `h1sweep.py` | the six-member sweep, one gamma at a time, nothing frozen |
| `chal.py` | the challenge-panel `A` / `B` / `claim` / `rho` measurement and the `X_c` solve |
| `fine.py` | the gamma response at 0.025 resolution, without emission (written, not run) |
| `controls.py` | the three controls of §1, run with the mechanism present and absent |
| `report.py`, `decide.py` | the h1-test tables and rule 13's composed slot gate |
| `METHOD.md` | the audit trail, including two operational traps this experiment hit |

# 096 — why `expr_mse_unbiased_capped_norm` reads 1.3627, decomposed

> ## ⚠️ CORRECTED 2026-08-29 — `min(claim, B)` is a PER-CONTEXT quantity and this document averaged it
>
> Every `claim`, `ρ`, `s_c` and `MSE_min` below applies `min(claim, B)` to the **context means**.
> The scorer forms #348's budget inside each context's own scoring run, and `min` is concave, so
> the two agree only while all three contexts are on the same side of their own cliff. They are
> for v11, v12, a8 and a10; they are **not** for a11, a12 or a13 — the arms this document's own
> recipe puts on the cliff. Measured on the shipped a12 artifact the gap is **0.005141**, worth
> **+0.0218 raw**, and it is the whole of a12's prediction miss.
>
> Corrected constants, from three board points whose `claim` never touches the `C_p` proxy:
> **`X = 0.019240`, `ΣD = 0.229141`, deconfounded cosine `0.07556`** (this document: 0.018846 /
> 0.228572 / 0.0740; the a8+a11 re-anchoring in `A12-PREDICTION.md`: 0.018089 / 0.235508 /
> 0.0701 — that one was the contaminated fit). All three `C_p`-free board points then fit to
> `≤ 0.0001` raw.
>
> Also corrected here: intervention 1's `VCC_ALPHA = 1.012` is the **mean-form** argmin. The
> per-context argmin is **0.982**, worth a further `+0.0035` overall. And `pds`/`nmae` are NOT
> invariant to amplitude — two board points give slopes of `−0.013` and `−0.014` scaled per unit
> `VCC_ALPHA`.
>
> The identity (★) itself is unaffected; it is arithmetic. See **`A12-ARTIFACT-AUDIT.md`** for
> the measurements, the attribution ladder and what the decomposition may still be trusted for.

**Headline.** The excess is mostly forfeited sampling credit, not directional error — but not
entirely, and an earlier draft of this document overstated that (see the retraction in §(b)). Written in
the metric's own units, `MSE − 1 = (A − 2X + Kc − ρ·claim)/ΣD`, and on the real challenge panel
for the current entry the four terms are `+1.2379 − 0.1649 + 0.0237 − 0.7296 = +0.3671`
(observed `+0.3627`; the 0.0044 gap is the whole residual of the model). We over-pay `A` — our
predicted effect's squared norm is **1.24× the reference's entire corrected effect** — and we
**under-collect `claim`** — the sampling correction the scorer is willing to hand us is
`0.24653`, capped at the reference's own `C_p`, and we are claiming `0.16676`, i.e. **68 % of
what we are entitled to**. Both are fixable, and the second is nearly free.

**The recipe, with numbers:** raise the emitted cells' dispersion until `Ĉ_p ≥ C_p` on every
row (top-8,000 genes stay at a10's `spread 0.15, φ=16`, the tail goes to `spread 1.0, φ=4`) and
cut `VCC_ALPHA` from **1.0796 to 1.012** — a 6.3 % amplitude cut, not the 60 % the naive
least-squares reading implies. Predicted raw `mse` **0.878**, `score_mse` **+0.115**, overall
**+0.0192**, against a total DE cost of about **−0.0013** — `−0.0010` for the α cut and `−0.0003`
for the call-set churn, both measured rather than assumed. A donor-agreement shrinkage toward the
CONTROL on top is worth a further **+0.0102** (§5, which retracts an earlier panel-mean version).

Everything below is *measured*, *derived* or *assumed*, and each is labelled.

---

## 0. What was measured, and how it was checked

`experiments/096-mse-floor/msekit.py` calls `cell_eval2`'s own primitives —
`prep.bulk_lognorm_means`, `moments.jackknife_correction`,
`metrics.delta.mse_unbiased_capped`, `metrics.delta.distance_unbiased`, aggregated as
`ratio_of_sums` per `catalog.py:345` — on cached per-group sums. It is not a reimplementation.

**MEASURED, bit-exact.** On `data/harness/h1-test.h5ad` against the production prediction it
returns `2.0551018300` where the full `031-somatic` → `cell-eval2 run` pipeline returned
`2.0551018299736388` (`experiments/095-clean-harness/out-h1test/.../run/agg_results.csv`), with
the two legs matching to every printed digit. Every number in this document that is not a board
number comes through that path.

Config that matters: `bulk_target_sum = 5e4`, `control_source: real`, and
`filter_gene_min_cpm_cell = 5.0` is **DE-only** (`distances.py:98`, `de_compute.py:414`) — the
expression metric runs on all 18,533 genes minus each row's own target.

---

## (a) The decomposition

### The identity

For row `p`, on its own scored gene set (all genes bar its target), with `w_p = 1/G_p`:

```
N_p = ( ||b̂_p − b_p||²  − ρ·min(Ĉ_p, C_p) − C_p ) / G_p
D_p = ( ||b_p − b_ctrl||² − C_p − C_ctrl )        / G_p        MSE = Σ N_p / Σ D_p
```

Write `d_p = b_p − b_ctrl` (real effect) and `dĥ_p = b̂_p − b_ctrl` (predicted effect). Then
`||b̂_p − b_p||² = ||dĥ_p||² − 2⟨dĥ_p,d_p⟩ + ||d_p||²`, and defining

| term | definition | needs |
|---|---|---|
| `A` | `Σ_p w_p ‖dĥ_p‖²` | our file + the real control **only** |
| `X` | `Σ_p w_p ⟨dĥ_p, d_p⟩` | the withheld reference |
| `R` | `Σ_p w_p ‖d_p‖²` | the withheld reference |
| `K` | `Σ_p w_p C_p` | the reference (estimated here) |
| `Kc` | `Σ_p w_p C_ctrl` | the real control — **exact**, it ships in `data/challenge/` |
| `claim` | `Σ_p w_p min(Ĉ_p, C_p)` | our cells + `C_p` |
| `ρ` | `min(1, B/claim)`, `B = _across_pert_budget` | our file **only** |

gives `Σ N = A − 2X + R − ρ·claim − K` and `Σ D = R − K − Kc`, hence

```
        MSE − 1  =  ( A − 2X + Kc − ρ·claim ) / ΣD                             (★)
```

**MEASURED:** (★) reproduces `cell_eval2`'s own value to `8.9e-16` on h1-test and `2.1e-14` on
orion-thinned (`decompose.py`). It is an identity, not a fit.

### Which term is inflated — the two harnesses

| | h1-test (85 perts, 400 cells, 20k UMI) | orion-thinned (132, 150 cells, 12k UMI) |
|---|---:|---:|
| `MSE` | 2.0551 | 1.9405 |
| `A/ΣD` | **+1.4279** | **+1.5798** |
| `−2X/ΣD` | −0.1241 | −0.0818 |
| `+Kc/ΣD` | +0.0093 | +0.0187 |
| `−ρ·claim/ΣD` | −0.2581 | −0.5762 |
| `ρ` | 1.000 (budget/claim 5.46) | 1.000 (2.71) |
| aggregate cosine `X/√(A·R)` | 0.0432 | 0.0176 |

`A` is the inflated term on both surfaces. Note the *shape*: `A/ΣD ≈ 1.4–1.6` while
`2X/ΣD ≈ 0.1`. The prediction has roughly the right **magnitude** and almost no **direction**.

### The same decomposition on the REAL challenge panel

`data/challenge/context_{A,B,C}.h5ad` are the reference bundle's own control cells (18,400 =
46 ntc_id × 400 per `manifest.json`), and the metric uses the reference control as the origin on
both sides. So `b_ctrl`, `C_ctrl`, `A`, `B`, `Ĉ_p` are all **exactly** computable for a shipped
`.vcc` with no access to the answers. `mech_chal_stats.py` streams a 360,000-cell / 2.5e9-nnz
prediction group-by-group out of h5py in ~50 s.

`C_p` (the reference's own jackknife) is the one estimated quantity. **MEASURED proxy:** the
delete-1 jackknife of a single ntc construct's 400 cells — same cell count, same 20k median
depth, same construct-level structure as a reference perturbation group — reads
**14.544 / 14.772 / 16.372** on A / B / C (16 constructs each, sd 0.39 / 0.53 / 0.73). A pooled
46-construct 400-cell draw gives 14.547 / 14.691 / 16.814, so between-construct variance is
negligible and the proxy is stable to ~3 %.

| term, context-mean | v11 (`VCC_ALPHA` 1.0796, φ=0) | a8 = **Wildflower** (φ=4) | a10 (`s0.15 φ16`, built) |
|---|---:|---:|---:|
| `A` | 0.283614 | 0.282948 | 0.282948 |
| `B` (budget) | 0.281236 | 0.280503 | 0.280503 |
| `claim` | 0.118737 | 0.166762 | 0.169139 |
| `Kc` | 0.005410 | 0.005410 | 0.005410 |
| `B/claim` → `ρ` | 2.369 → **1.0** | 1.682 → **1.0** | 1.658 → **1.0** |
| `#247` cap fires on | 0 of 900 rows | 0 of 900 | 0 of 900 |

**So: (iii) is real and it is large.** `ρ = 1`; we forfeit nothing to the #348 budget. What we
forfeit is to the **#247 cap's ceiling from below**: `claim` may be as large as
`Σ_p w_p C_p = 0.24653`, and Wildflower claims `0.16676` — **our emitted cells are only 0.67× as
dispersed as the reference's** (`Ĉ_p` mean 10.30 against `C_p` 15.23), so we leave 0.0798 of
correction on the table. In MSE units that is **0.349**, which is most of our 0.363 excess.

### Where the two unknowns come from — a board-anchored solve, with a falsifier

`X` and `ΣD` are the only quantities in (★) we cannot compute. Two of our own board points
differ **only in `VCC_ALPHA`** (`v11 → v12` was one env variable; git `aaf9ff3`):

* v11, `VCC_ALPHA` 1.0796 → `s = 1.00000` → raw **1.5801** (`docs/STRATEGY.md` addendum 5)
* v12 "Alpha 0.40", `VCC_ALPHA` 0.40 → `s = 0.37050` → raw **0.9640** (leaderboard)

`A(s) = s²A`, `B(s) = s²B`, `X(s) = sX`, and `claim(s)` is flat — **MEASURED**: re-emitting v11's
profile at `s ∈ {0.15, 0.25, 0.37, 0.55, 0.85}` through the production `emit_counts` moves the
mean `Ĉ_p` by **0.01 %** (6.9024 → 6.9037), and the same emitter reproduces the shipped `Ĉ_p` to
1.0 %. Two equations, two unknowns:

```
   X = 0.018846        ΣD = 0.228572          (refits v11 to 1.5801 and v12 to 0.9640 exactly)
```

**FALSIFIER, and it holds.** a8/Wildflower was not used in the fit and differs from v11 in a
completely different way — `claim` 0.1187 → 0.1668 (φ=4) and `A` 0.28361 → 0.28295 (pruned
whitener). Predicted raw `mse` **1.3671**; observed **1.3627**; error **+0.0044 (+0.32 %)**.

**Second, independent check.** v15 "Half a step back" (raw 1.3236) has v11's profile exactly —
`tail-spread 1.0` / `keep-top 2000` changes only emission, and the exact-total repair leaves the
profile bit-identical. (★) then requires `claim_v15 = 0.17737`. **MEASURED** on the challenge
panel: uniform `spread 0.25` gives 0.1168, uniform `spread 1.0` gives 0.1922 — and v15's mixed
0.25/1.0 policy must land between them. It does, at 78 % of the way. Nothing was fitted to v15.

### The challenge panel, fully resolved

| quantity | value | `/ΣD` | status |
|---|---:|---:|---|
| `ΣD` — the reference's corrected effect | 0.228572 | 1.0000 | **solved** (2 board points) |
| `A` — our predicted amplitude, α=1.0796 | 0.283614 | +1.2408 | **measured** |
| `X` — alignment with the truth | 0.018846 | +0.0825 | **solved** |
| `Kc` — the control's correction | 0.005410 | +0.0237 | **measured, exact** |
| `claim` — φ=0 / φ=4 / **at the #247 cap** | 0.1187 / 0.1668 / **0.2465** | 0.519 / 0.730 / **1.079** | **measured / estimated cap** |
| `K = Σ w_p C_p` (proxy) | 0.246529 | 1.0785 | **estimated ±3 %** |
| `R = ΣD + K + Kc` (plug-in real distance) | 0.480511 | 2.1022 | derived — 47.6 % of it survives its own correction |
| **deconfounded cosine** `X/√(A·ΣD)` | **0.07402** | | derived |
| plug-in cosine `X/√(A·R)` | 0.05105 | | derived |
| unconstrained least-squares `s* = X/A` | 0.06645 | | derived |

Two readings worth stating plainly:

1. **Our effect direction is right about 7.4 % of the way**, aggregate-weighted. That is the
   number that bounds everything. It is not zero — a zero-signal prediction cannot beat the
   baseline at all — but it is small.
2. **Our amplitude is not absurd.** `√(A/ΣD) = 1.114`: our predicted effect is 11 % larger than
   the reference's *corrected* effect and 23 % *smaller* than its plug-in one. The prediction
   is not "twice optimal" in any physical sense; it is optimal-magnitude and near-orthogonal.

---

## (b) The optimal global α — and why it is **not** `X/A`

Substituting `dĥ_p → s·dĥ_p` into (★), with `claim` flat in `s` (measured above):

```
   MSE(s) − 1 = [ s²A − 2sX + Kc − min(s²B, claim) ] / ΣD
```

Two regimes, because `ρ = min(1, s²B/claim)` and `ρ·claim = min(claim, s²B)`:

* **`s ≥ s_c = √(claim/B)`** (`ρ=1`): slope `2(sA − X)/ΣD`, positive for all `s > X/A = 0.066`.
* **`s < s_c`** (`ρ<1`): the numerator becomes `s²(A − B) − 2sX + Kc`. **`A − B` is exactly the
  squared norm of our predictions' common mode** (MEASURED: 0.002378, i.e. **0.84 % of `A`**),
  so the slope `2(s(A−B) − X)/ΣD` is negative out to `s = X/(A−B) = 7.9`.

The function therefore **decreases up to `s_c` and increases after it**: the minimum sits
exactly on the ρ cliff.

> This is the substantive correction to the naive reading. Below the cliff the metric refunds
> the whole differential part of the amplitude through #348's budget, so shrinking further buys
> only the common mode (0.84 % of A) while giving up `2sX` of alignment credit and, at the
> cliff, the entire `claim`. **`X/A = 0.066` is not reachable; `s_c` is the answer.**

**MEASURED end-to-end.** `mech_alpha_h1.py` re-renders, re-emits through production
`emit_counts`, re-jackknifes and re-scores at each `s` on h1-test. Nothing frozen:

| `s` | 0.10 | 0.20 | 0.30 | 0.40 | **0.42** | 0.43 | 0.44 | 0.46 | 0.50 | 0.70 | 1.00 | 1.20 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `mse` | 1.0004 | 0.9886 | 0.9771 | 0.9659 | **0.9637** | 0.9651 | 0.9763 | 0.9996 | 1.0496 | 1.3680 | 2.0785 | 2.6758 |
| `ρ` | 0.056 | 0.220 | 0.492 | 0.874 | 0.963 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

Predicted cliff `s_c = √(0.039656/0.21551) = 0.4290`; measured argmin at `s = 0.42`, with `ρ`
reaching 1 between 0.42 and 0.43 and the curve turning up there. The V has a **kink**, not a
smooth minimum — the signature of the constraint rather than of the least-squares fit.

**On the challenge panel**, with the solved `X` and `ΣD`:

| emitter | `claim` | `s_c` | `VCC_ALPHA*` | raw `mse` | `score_mse` | overall |
|---|---:|---:|---:|---:|---:|---:|
| v11 `spread 0.25, φ=0` | 0.11874 | 0.6494 | 0.701 | 0.9210 | +0.0704 | +0.0117 |
| a8 `spread 0.25, φ=4` | 0.16676 | 0.7710 | 0.832 | 0.9029 | +0.0894 | +0.0149 |
| a10 `spread 0.15, φ=16` | 0.16914 | 0.7765 | 0.838 | 0.9021 | +0.0903 | +0.0150 |
| **at the #247 cap** | **0.24653** | **0.9375** | **1.012** | **0.8785** | **+0.1151** | **+0.0192** |

A useful closed form. At `s_c`, with `c = X/√(A·ΣD)`, `q = claim/ΣD` and `ε = (A−B)/B`, the
EXACT value is

```
        MSE_min  =  1 + Kc/ΣD  +  q·ε  −  2·c·√( q·(1+ε) )
```

`ε` is the predictions' common mode as a fraction of `B`, measured at **0.0087** for a8, so the
simpler form below — which drops `q·ε` and the `(1+ε)` factor — is accurate to 0.009 raw *for
families that preserve the common mode*. It is **not** safe for a family that shrinks toward the
panel mean, which is exactly what `_across_pert_budget` centres out; see the retraction in
intervention 5. Every number in intervention 1 was computed through the exact route
(`s²(A−B) + Kc − 2sX` with the `min(claim, s²B)` branch) and is unaffected. Dropping the two
terms:

```
        MSE_min  ≈  1 + Kc/ΣD  −  2·cos·√(claim/ΣD)
                 =  1.0237 − 2.077·cos          (at the cap, on this panel)
```

so the baseline `0.9879` needs `cos ≥ 0.0173`, and the leader's `0.8270` needs `cos ≥ 0.0947`
even at the cap. Our `0.0740` buys `0.870` at the cap — enough to clear the floor decisively,
not enough to reach the leader.

**RETRACTED — an earlier version of this section said the leader's directional edge was
`0.97×–1.28×` ours and "plausibly zero". That is wrong, and it is wrong in the direction that
mattered.** The table it rested on (contributed by the competitive-intel session, adopted here
without re-deriving it — the exact failure this document keeps documenting) inverted the role of
the cap. The correction:

`ΣD` and `Kc` are properties of the data and are the same for everyone; `claim` is a property of
the submission. Inverting the leader's 0.8270 requires `2·cos·√(claim/ΣD) = 1.02367 − 0.8270 =
0.19667`, so a leader collecting MORE credit needs LESS cosine:

| assumed leader emitter | `claim/ΣD` | coefficient | implied `cos` | vs our 0.0740 |
|---|---:|---:|---:|---:|
| φ=0, like our v11 | 0.51947 | 1.44149 | 0.13644 | 1.84× |
| φ=4, like our a8 | 0.72958 | 1.70831 | 0.11513 | 1.56× |
| **at the #247 cap — the maximum anyone can claim** | **1.07856** | **2.07707** | **0.09469** | **1.28×** |

**The cap is a ceiling on `claim`, so it yields the FLOOR on `cos`, not a mid-range estimate.**
The leader's aggregate cosine is **at least 0.0947, at least 1.28× ours**, and higher if their
emitter is anywhere short of the cap or their amplitude is off its optimum. There is no branch of
this table in which their direction is no better than ours; the earlier reading had the inequality
pointing the wrong way.

**The consequence is material and it changes the plan's ceiling.** Our own best possible raw, at
`cos = 0.0740` with `claim` at the cap, is **0.8699**. The leader is already at **0.8270**.
So **calibration alone does not reach them** — it reaches 4th on this member, not 1st. The honest
statement is: *most* of our deficit is forfeited noise credit, and closing it is cheap and worth
`+0.0192`; but there is also a real directional gap of at least 28 %, and it caps what
calibration can buy.

**Assumptions.** (1) `ΣD_c` and `X_c` are equal across the three contexts — the board publishes
only the pooled raw, so this cannot be checked; see (e). (2) `C_p` from the single-construct
control proxy, ±3 %. Sensitivity: a `C_p` 10 % lower moves `MSE_min` from 0.8784 to 0.8856,
i.e. the conclusion is insensitive. (3) `s` acts as an exact linear scalar on the log-space
effect. The pipeline's `VCC_ALPHA` does exactly that (`multi_donor.py:169`,
`make_submission_pm_cd4.py:970`) up to the renormaliser `Z ≈ 1`; the round-trip
`sums → profile → sums` at `s = 1` is exact to `2e-16` (`mech_alpha_chal.render`). One real
residual: a fresh build floors `n·L·profile` once more than this reconstruction, worth about
`+0.003` raw at `s = 0.4` on h1-test. It biases the answer *against* the intervention.

---

## (c) Per-gene / per-perturbation shrinkage — the oracle is huge, the buildable version is nil

Under the same constrained optimisation (maximise `Σ_g s_g V_g` subject to `Σ_g s_g² U_g = claim`,
with `U_g = Σ_p w_p dĥ_pg²`, `V_g = Σ_p w_p dĥ_pg d_pg`) the whole question collapses to one
scalar:

```
     MSE_min = 1 + ( Kc − 2√(claim · W) ) / ΣD
     W = X²/A                                   global scalar
     W = Σ_g V_g²/U_g                           per-gene
     W = Σ_p w_p ⟨dĥ_p,d_p⟩²/‖dĥ_p‖²            per-perturbation
```

The gain scales as **√W**, so a richer family is worth `√(W_family/W_global)`.

**MEASURED, oracle (fitted on the truth — an upper bound only):**

| family | h1-test `√(W/W_glob)` | orion `√(W/W_glob)` |
|---|---:|---:|
| global scalar | 1.000× | 1.000× |
| **per-gene** | **3.966×** | **6.548×** |
| per-perturbation | 2.013× | 1.969× |

**MEASURED, honest.** Two separate tests, both leak-free by rule 4:

*Cross-validated within a dataset* — fit `s_g` on half the perturbations, score on the other
half, ridge-shrunk toward the global coefficient, 30 random splits. Reported relative to the
same-protocol global baseline (λ→∞), which is the only fair comparator:

| λ | 0 | 1 | 3 | **10** | 30 | 100 | 1000 | →∞ (global) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| h1-test (42 fit rows) | 0.411 | 0.532 | 0.680 | **0.784** | 0.769 | 0.747 | 0.736 | 0.734 |
| orion (66 fit rows) | 1.087 | 1.152 | **1.230** | 1.215 | 0.997 | 0.832 | 0.754 | 0.745 |

so the best in-dataset per-gene gain is **1.068×** (h1-test) and **1.65×** (orion), rising with
the number of perturbations — which is encouraging for a 300-perturbation panel.

*Across datasets, which is the only construction a submission can actually ship* — fit `s_g` on
one cell line, apply it to the other:

| λ | 0 | 3 | 10 | 30 | **100** | 1000 | →∞ |
|---|---:|---:|---:|---:|---:|---:|---:|
| orion → h1-test | 0.143 | 0.287 | 0.579 | 0.936 | **1.023** | 1.004 | 1.000 |
| h1-test → orion | 0.122 | 0.355 | 0.726 | 0.974 | **1.009** | 1.002 | 1.000 |

**Verdict: no.** The optimal ridge is λ ≈ 100 — shrink the per-gene coefficients ~99 % of the
way to the global one — and the surviving gain is **+1 to +2 %** on `√W`, i.e. `≈ +0.003` raw
`mse`, which is inside the noise of everything else here. About **98 % of the 4–6.5× oracle
per-gene gain is unbuildable**: the per-gene coefficient is a cell-line property, not a model
property, and it does not transfer.

Per-perturbation is the same story from the other side. The oracle is 2.0×, but a
per-perturbation scalar has no out-of-sample form (the held-out unit *is* the perturbation), so
it must come from a *confidence score*. **MEASURED**, the buildable self-referential family
`s_p ∝ ‖dĥ_p‖^(−t)`:

| `t` | 0 | 0.25 | 0.50 | 0.75 | 1.00 | 1.25 |
|---|---:|---:|---:|---:|---:|---:|
| h1-test | 1.000× | 1.009× | **1.012×** | 1.009× | 1.002× | 0.993× |
| orion | 1.000× | 0.994× | 0.987× | 0.978× | 0.968× | 0.956× |

Nothing, and the two harnesses disagree in sign.

**The estimator I would actually ship** is therefore the global scalar, and it is fully
buildable *without any reference*, because the optimum is the **constraint**, not a fit:

```
   s* = √( claim / B ),      B     = cell_eval2.metrics.delta._across_pert_budget(our own bulks)
                             claim = Σ_p w_p min(Ĉ_p, Ĉ_p^ref)
```

`B` and `Ĉ_p` come from the submission alone; `Ĉ_p^ref` is the 400-cell single-construct
jackknife of the control cells the challenge ships us (measured above, sd 3 %). Nothing in `s*`
touches the answers. `mech_challenge.py` + `mech_alpha_chal.py` compute it in ~2 minutes from a
built `.vcc`.

The one shrinkage worth revisiting later is **per-gene fitted on a donor panel with ridge λ
selected by cross-cell-line CV** — the machinery is 20 lines and the measurement above says the
honest λ is large. Do not fit `s_g` against the harness reference and ship it; that is the
oracle leak of rule 4 and it is worth 0.14× out of sample.

---

## (d) Per-context clamping: **confirmed**, with a proof rather than an inference

**MEASURED.** Fit `score_mse = a·raw + c` on the six well-behaved positive-`score_mse` teams in
`data/leaderboard/lb_20260828T1325.json`:

```
   score_mse = clip( (0.987862879 − raw) / 0.950162566 , 0, 1 )
   max |residual| 2.1e-5 over Aginglab.com, Ibrahim Mansour, BioFarmer, IsoCell,
   black_handsome_egg, enough1  (raw 0.827 … 0.969)
```

`b = 0.98786` and span `0.95016` sit inside the spec's published ranges (`b` 0.986–0.992, span
"about 0.95"), so this is the scorer's own scale, recovered.

Three of the nine positive teams sit **above** that line and **none** sits below:

| team | raw | `score_mse` | affine map | residual |
|---|---:|---:|---:|---:|
| We527 | 1.3464 | 0.053198 | −0.377 | **+0.4306** |
| Immovable Myhal Rock | 0.9906 | 0.005867 | −0.0029 | +0.0087 |
| nanobot | 0.9918 | 0.001404 | −0.0041 | +0.0055 |

Zero of 385 teams have a negative `score_mse`. One-signed residuals are exactly what
`mean_c clip(z_c)` produces against `clip(mean_c z_c)`: `clip(z) − z = max(0, −z)` for every
`z ≤ 1`, so the residual is a sum of non-negative terms, always.

**And there is a proof, not just a signature.** We527 has two submissions:

```
   2026-08-27   raw 1.3576867034950901   score_mse 0.05319847662453766
   2026-08-28   raw 1.3464492032778235   score_mse 0.053198476624537694
```

**Two different raw values map to the same score to fifteen significant figures.** No global
map — monotone or not, affine or not — can be simultaneously (i) locally constant across
`[1.3464, 1.3577]`, (ii) positive there, and (iii) zero at 1.0040 (`athenen`) and at every one
of the other 375 zeroes. Per-context clamping does it trivially: We527's improvement lived
entirely inside the two contexts that were already clamped to 0, and the third context's raw was
byte-identical between the two builds. The board's `anchor_version`
`vcc2026-valA-r4+vcc2026-valB-r4+vcc2026-valC-r4` names three per-context anchor sets, and the
metric spec says contexts "are scored separately and averaged". **Confirmed.**

Inverting We527's residual: `3 × 0.430593 × 0.950163 = 1.2274` of raw is sitting above `b` in
the clamped contexts. Either one context at raw ≈ 0.836 with two at ≈ 1.60, or two at ≈ 0.912
with one at ≈ 2.21 — both reproduce their 0.0532 and the data cannot separate them.

### Strategic consequence

The member's contribution is `(1/3) Σ_c max(0, (b − u_c)/span)` — the upper clamp at 1 needs
`raw ≤ 0.0377` and is unreachable. That function is **convex** in the raw vector `(u_A,u_B,u_C)`.
Maximising a convex function over a convex feasible set puts the optimum at a **corner**:

* **If you cannot get all three contexts below `b`, concentrate.** Every unit of raw improvement
  spent in a context that stays above `b` is worth *exactly zero*, and the gradient above `b` is
  zero in both directions — there is no penalty for being terrible there. Excellent-in-one
  strictly dominates uniformly-mediocre.
* **If you can get all three below `b`, it does not matter how you spread it.** In that region
  the score is `(3b − Σu_c)/(3·span)`, a function of the **mean raw** alone, and concentration is
  exactly neutral until a context hits 0.0377.
* Corollary for any risky intervention: apply it **per context**, and never let a global
  parameter be set by a context that is hopeless. Also — a context already above `b` is a free
  place to spend amplitude on the DE members.

**For us specifically this does not bind** (see (e)): our predicted-side terms are near-identical
in all three contexts, and the intervention in §(b) moves all three together, from ≈1.34–1.40 to
≈0.87–0.89. We are in the "all three cross" regime, so the mean is the right objective. The
finding still matters as a *risk* statement — it means a partial fix is worth more than a linear
reading of the pooled raw suggests, because a context that crosses pays fully.

Board arithmetic, using the recovered map: one context to 0.90 = `+0.0051` overall; one context
to 0.83 = `+0.0092`; all three to 0.878 = `+0.0192`; all three to 0.83 = `+0.0277` (rank 1).

---

## (e) Our per-context `mse`

**Exactly measurable from what we hold** (`mech_challenge.py`, Wildflower = the current entry):

| ctx | `A_c` | `B_c` | `claim_c` | `Kc_c` | `C_ctrl` | `C_p(400)` | `K_c` | `num_c = A+Kc−claim−2X` | `MSE_c` * |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 0.286883 | 0.284448 | 0.162190 | 0.005084 | 0.3141 | 14.544 | 0.23543 | 0.092084 | 1.4029 |
| B | 0.274810 | 0.272386 | 0.164652 | 0.005228 | 0.3229 | 14.772 | 0.23913 | 0.077693 | 1.3399 |
| C | 0.287151 | 0.284676 | 0.173444 | 0.005919 | 0.3656 | 16.372 | 0.26503 | 0.081934 | 1.3585 |
| mean | 0.282948 | 0.280503 | 0.166762 | 0.005410 | | | 0.24653 | 0.083904 | 1.3671 (board 1.3627) |

\* **ASSUMED**, and this is the honest limit of the exercise: `MSE_c` needs `ΣD_c` and `X_c`, and
both live entirely on the withheld reference. The column above holds them at their fitted pooled
values. What *is* measured is that our predicted-side terms are nearly identical across contexts
(`A_c` spans 4.5 %), so **any large per-context spread in our `mse` must come from the reference
side, not from us.** Sensitivity: with `num_A = 0.0921` fixed, `ΣD_A ∈ [0.19, 0.27]` puts `MSE_A`
anywhere in `[1.34, 1.49]`.

One directional inference that *is* grounded in measurement: context C's sampling corrections are
12 % larger than A's (`C_p` 16.37 vs 14.54, `C_ctrl` 0.366 vs 0.314). Since `ΣD = R − K − Kc`, if
the three contexts' true effects `R_c` are comparable then `ΣD_C` is the smallest and **context C
is our worst**, by roughly `+0.13` raw. That is a hypothesis with a cheap falsifier: submit once
with a per-context α and read the pooled raw against the prediction.

Two things we cannot estimate at all: `X_c` (per-context alignment) and the split of the pooled
1.3627. If per-context numbers are needed, the only instrument is the board, and it takes one
submission with deliberately different α per context to resolve them.

---

## The interventions, ranked

| # | intervention | raw `mse` | Δoverall | cost | status |
|---|---|---:|---:|---|---|
| 1 | tail over-dispersion to the #247 cap + `VCC_ALPHA` 1.0796 → 1.012 | **0.878** | **+0.0192** | α cut −0.0010, call-set churn −0.0003 — **both priced** | measured on both harnesses |
| 2 | a10 as built at `VCC_ALPHA` 0.838 (one env var) | 0.902 | +0.0150 | reach −0.0036, nmae −0.0025 | measured, shippable today |
| 3 | overshoot the cap 1.3–1.5× as insurance | — | 0 | 0 | free |
| 4 | sub-gate-only over-dispersion (invisible to every DE member) | 0.898 | +0.0007 | 0 | measured, genuinely free |
| 5 | per-gene confidence shrinkage toward the CONTROL, donor-agreement α | **0.820** | **+0.0102** | must be realised | measured on both harnesses |
| 6 | raise the direction cosine | −0.021 raw per +0.01 | +0.0036 per +0.01 | the whole modelling problem | open |
| 7 | per-gene slope `s_g` fitted on a donor | ~1.36 | +0.0003 | — | **do not build** |

Baseline for every row: Wildflower, raw 1.3627, `score_mse` 0, overall 0.154361. `score_mse` uses
the recovered map; overall = `score_mse/6`. Reach costs use the board anchors
`b = 0.079022, r = 0.964145` (span 0.885).

### 1. Raise `claim` to the #247 cap, then cut α by 6 % — **+0.0192 overall**

`claim` is the only term in the numerator that emission moves, because the quasi emitter repairs
each gene's panel total exactly, so `A`, `B`, `X`, `pds` and `nmae` are all invariant to
`(spread, φ)` by construction (`donor_transfer.py:290-300`; verified as rule 5 requires — `nmae`
is printed as a guard in every exp095 realisation and does not move).

**MEASURED** `claim(spread, φ)` on the challenge panel (`mech_claim_grid2.py`,
`mech_split_policy.py`; cap = 0.24653):

| policy | `Ĉ_p` | `Ĉ/C_p` | `claim` | `VCC_ALPHA*` | raw `mse` | `score_mse` | overall |
|---|---:|---:|---:|---:|---:|---:|---:|
| `0.25, 0` (v11) | 7.33 | 0.48 | 0.1186 | 0.701 | 0.9210 | 0.0704 | +0.0117 |
| `0.15, 16` (a10, built) | 10.43 | 0.69 | 0.1688 | 0.838 | 0.9022 | 0.0902 | +0.0150 |
| `0.50, 4` | 21.62 | 1.42 | **0.2465** | 1.012 | **0.8785** | 0.1151 | **+0.0192** |
| **top-8000 `(0.15,16)` + tail `(1.0,4)`** | 17.68 | **1.16** | **0.2465** | **1.012** | **0.8785** | **0.1151** | **+0.0192** |
| top-5000 `(0.15,16)` + tail `(1.0,4)` | 33.31 | 2.19 | 0.2465 | 1.012 | 0.8785 | 0.1151 | +0.0192 |

Both halves are needed: raising `claim` to the cap at the shipped `VCC_ALPHA = 1.0796` reads
raw **1.0181** — still above the baseline — and cutting α alone tops out at 0.902 (row 2 of that
table). Together they read 0.8785.

The **split policy is the shippable one**: the machinery already exists
(`emission_spread_vector`'s `keep_top` / `keep_spread` / `tail_spread`, and `emit_counts` takes a
per-gene `overdispersion`), and top-8,000 reaches the cap at only 1.16× the reference's
dispersion, so nothing is over-cooked. Then set `VCC_ALPHA = 1.012`.

**Cost in `reach`.** Two components, both small, and the first is *negative*:

* The dispersion half **helps**. Board-measured: `a2-prune` (φ=0) reach raw 0.17015 →
  `a8` (φ=4) 0.18900, i.e. `+0.0213` scaled = **+0.0036 overall**. exp095's realisations agree
  (`s0.15_φ16` takes harness reach 0.0838 → 0.1565).
* The amplitude half costs. `VCC_ALPHA` 1.0796 → 1.012 is a **6.3 %** cut. Using the steepest
  slope our board supports (reach 0.1358 at α=0.40, 0.1890 at α=1.0796 → 0.0783 raw per unit α,
  an over-estimate because that pair is confounded), `Δreach ≤ −0.0054` raw = `−0.0061` scaled =
  **`−0.0010` overall**. `nmae` costs `−0.0007` overall on the same reasoning
  (`docs/STRATEGY.md`'s α ladder). `pds` is scale-invariant.

**Net = +0.0192 − 0.0010 − 0.0007 = +0.0175 overall**, i.e. `0.15436 → 0.1719`, rank ≈ 4–5 on
today's board. The `+0.0036` reach gain from φ is *not* counted — Wildflower already carries
φ=4, so it is banked in the baseline; exp095 prices a10's further emission changes separately.

**Risks — and the call-set cost is now MEASURED, not argued.** `mech_callability.py` runs the
real path on challenge context A: production `emit_counts`, then the scanpy-identical rank-sum
against the **real 18,400 control cells**, then BH at 0.05, over 15 perturbations.

| tail policy | `n_pred` | vs a10 | `Ĉ_p` | `Ĉ/C_p` | `claim` | fraction of a10's calls kept |
|---|---:|---:|---:|---:|---:|---:|
| a10 uniform `(0.15, φ16)` | 4772 | 1.000 | 10.18 | 0.70 | 0.16478 | 1.000 |
| **top-8000 by `|dĥ|`** | **5417** | **1.135** | 17.88 | 1.23 | **0.23543** | **0.913** |
| control CPM ≥ 40 | 5595 | 1.173 | 20.02 | 1.38 | 0.23543 | 0.858 |
| control CPM ≥ 20 | 5129 | 1.075 | 13.16 | 0.90 | 0.21298 | 0.865 |
| control CPM ≥ 10 | 4823 | 1.011 | 10.93 | 0.75 | 0.17691 | 0.900 |

Two things fall out, one reassuring and one that **retracts an assumption the recipe was resting
on**:

* **`n_pred` goes UP, not down** — +13.5 % under top-8000. `fid = k/max(n_pred, n_real)` is
  therefore in no danger from below; `n_real ≈ 4,757` stays covered. This is the opposite of the
  intuition that more dispersion means fewer calls, and it is consistent with
  `donor_transfer.py:261-265` — a near-deterministic emitter parks 400 identical values at the
  control's *mean*, above its median, and the rank test reads that as signal; restoring spread
  removes the artefact on those genes and lets genuinely-moved ones through.
* **The "spend the dispersion where nothing was callable" argument is FALSE on our submission,
  and I had assumed it.** MEASURED: our called genes sit at median control CPM **48.1**, our
  uncalled genes at **87.6** — expression predicts callability in the *opposite* direction here.
  Zhou Li's diagnosis (#49: callable genes at 107 CPM, uncallable at 38.7) does not reproduce on
  our emitter, because our callability is set by our own effect sizes and spread policy rather
  than by depth. Concretely, genes at control CPM ≥ 40 are only **62.7 %** of all our called
  genes, so a CPM-keyed tail over-disperses more than a third of them. The sub-gate route of §4
  is the only genuinely-invisible one, and §4 already measures it as too small (0.4 % of counts).

So the tail policy **churns 9–14 % of the call set** — top-8000-by-`|dĥ|` is the gentlest at
8.7 % churn and it still reaches the cap.

**And that churn is now PRICED, because `fid` factorises.** I had written this off as unreadable
under rule 1. That was too quick:

```
   F_p = k / max(n_pred, n_real) = (k/n_pred) · min(1, n_pred/n_real)
```

The **coverage** leg carries `n_real` and is indeed unreadable on our harnesses. The **precision**
leg `k/n_pred` carries no `n_real` at all — it needs only our own call set and the reference's
signs on those genes, both of which the harnesses have. And on the challenge `n_pred` 4772 already
exceeds `n_real ≈ 4757`, so coverage is **saturated at 1** and the tail policy raises `n_pred` to
5417, keeping it saturated. **Precision is therefore the only thing that can move `fid`, and the
whole risk reduces to one scalar: the directional precision of the calls that are NEW.**

**MEASURED** (`mech_precision.py`, both harnesses per rule 1, 25 perturbations each, same profile
emitted under both policies, exact DE path against the real controls, scorer's own-gene and
scored-universe exclusions applied):

| | h1-test | orion-thinned |
|---|---:|---:|
| `n_pred` a10 → tiered | 4198 → 6210 (1.479×) | 3598 → 4745 (1.319×) |
| precision `k/n_pred`, a10 → tiered | 0.5110 → 0.5062 | 0.5565 → 0.5456 |
| **precision of the NEW calls** | **0.4946** | **0.5131** |
| precision of the DROPPED calls | 0.4938 | 0.5132 |
| jaccard, a10 → tiered | 0.0751 → 0.0814 | 0.0121 → 0.0113 |

**The new calls land at chance, and the two surfaces bracket 0.500 from either side.** The
pessimistic branch — new calls actively anti-correlated with the reference — is refuted on both.
The dropped calls are at chance too, so the churn is chance-in, chance-out.

**The cost has a closed form.** When the churn is chance-in *and* chance-out, the above-chance
mass `E = k − 0.5·n_pred` is conserved — chance-level calls carry none of it, in either direction —
and only the denominator grows. So

```
   p_after − 0.5  =  (p0 − 0.5) · n_before / n_after
   Δscore_fid     =  −(1/0.29725) · (p0 − 0.5) · (1 − n_before/n_after)      [1/span = 3.3642]
```

Checked against the two harness measurements, which it was **not** fitted to:

| surface | `p0` | `n_pred` before → after | closed form | observed |
|---|---:|---:|---:|---:|
| h1-test | 0.5110 | 4198 → 6210 | 0.50744 | 0.5062 |
| orion-thinned | 0.5565 | 3598 → 4745 | 0.54284 | 0.5456 |

The two surfaces bracket it, h1-test 0.0012 below and orion 0.0028 above.

The formula states the transfer argument exactly: **the cost is proportional to the headroom you
have to lose.** Headroom above chance is +0.0110 on h1-test, +0.0565 on orion, and **+0.0041 on
the board** — we have the least to lose of the three, so the board pays least. Nothing in the cost
depends on the new calls being *good*; only on their not being actively anti-correlated, which is
now ruled out on both surfaces.

Transferring the *mechanism* rather than the magnitude (rule 1 done properly), with the board's own
constants — `p0 = 0.5041`, `n_pred` 4772 → 5417, 91.3 % kept, so 415 dropped and 1060 new:

| assumption about the new calls | precision | Δ`score_fid` | Δoverall |
|---|---:|---:|---:|
| at our own precision (optimistic) | 0.5041 | +0.0000 | +0.0000 |
| **at chance, drop leg accounted — the closed form** | **0.50361** | **−0.0016** | **−0.00027** |
| at chance, kept calls held at `p0` — *my first pass, wrong* | 0.50330 | −0.0027 | −0.00045 |
| 1.4 pp *below* chance (refuted) | 0.5013 | −0.0093 | −0.00155 |

**Correction to my own arithmetic, caught by the competitive-intel session.** My first pass held the
KEPT calls at `p0` while adding new calls at chance. That double-counts the dilution: I had also
measured the DROPPED calls at chance (0.4938 / 0.5132), and removing chance-level calls from an
above-chance set leaves the remainder slightly *enriched*, which partly offsets. The point estimate
is **−0.00027**, about 40 % cheaper than I quoted. It is inside the `±0.004` band either way and
changes no decision, but the number in the table is the corrected one.

The harness `Δprecision` (−0.005, −0.011) looks worse than the board's (−0.0008) for a reason
worth stating: **precision headroom above chance is +0.0110 on h1-test and +0.0565 on orion, but
only +0.0041 on the board.** We are 0.41 % better than a coin flip on the genes we call, so
diluting with chance-level calls costs almost nothing. `jac` moves +0.0062 / −0.0008 — opposite
signs, i.e. zero within noise, and its magnitude is genuinely unreadable here because its union
carries `n_real`.

**Verdict: the dispersion half costs about `−0.0003` overall, not the `−0.028` a pessimistic
reading allowed.** Combined with the α cut's `−0.0010`, intervention 1's total DE cost is
`≈ −0.0013` against a `+0.0192` gain. Still realise it through `realize.py` before shipping — but
it is no longer an open risk, it is a priced one.

Remaining risks. (i) A better key than `|dĥ|` is the production `rank_score`/`keep_key` path,
which is p-value ordered — MEASURED on the orion DE table, top-8,000-by-`|lfc|` holds only 77 %
of our called set, so the key is doing real work and the p-value order should do it better.
(ii) A very over-dispersed tail can produce empty rows, which `make_submission_pm_cd4.py:1285`
raises on.

### 2. Amplitude alone, no emission change — **+0.0150 overall**

Ship a10 (`s0.15 φ16`, already built and on disk) with `VCC_ALPHA = 0.838` instead of 1.0796.
Predicted raw `mse` 0.9021, `score_mse` +0.0903. Cost: a 22 % amplitude cut, so
`Δreach ≈ −0.0189` raw = `−0.0214` scaled = `−0.0036` overall, `nmae ≈ −0.0025` overall.
**Net ≈ +0.0089.** Strictly dominated by #1 but needs no new emitter code — one env variable.

### 3. Push `claim` past the cap as insurance — **+0.000, but it removes a failure mode**

`min(Ĉ_p, C_p)` saturates: once `Ĉ_p ≥ C_p` on every row there is no further gain and no
penalty. Since `C_p` is only *estimated* (±3 %), targeting `Ĉ_p ≈ 1.3–1.5 × Ĉ_p^ref` rather than
`1.0×` costs nothing in `mse` and removes the risk of under-claiming. Free.

### 4. Sub-gate-only over-dispersion — **+0.0007 overall, and it is genuinely free**

Genes below the DE gate (`filter_gene_min_cpm_cell = 5.0`) are invisible to all four DE members
and to `nmae`, but fully visible to `Ĉ_p`. **MEASURED**: 46–48 % of the 18,533 genes are under
5 CPM in the real controls — but they hold only **0.40–0.47 % of the counts**, so scattering all
of them at `spread 1.0, φ=16` lifts `Ĉ_p` from 10.43 only to 11.14 (`claim` 0.1688 → 0.1804,
`+0.0007` overall). Worth taking as a rider on #1, not worth building for.

### 5. Per-gene confidence shrinkage with a donor-agreement α — toward the CONTROL, **+0.0102**

**RETRACTED AND REPLACED. An earlier version of this section shrank toward the PANEL MEAN
(Miller et al., `μ_ID = α⊙μ_TD + (1−α)⊙μ_all`) and scored it at `+0.0044` to `+0.0097` overall.
That is wrong, and the sign is wrong: the panel-mean family is `−0.0030` against not shrinking at
all.** Caught by the exp096-shrinkage session while implementing it
(`mech_shrink_target.py`); re-derived independently here (`mech_shrink_exact.py`).

**The error.** I scored the family through `MSE_min = 1 + (Kc − 2√(claim·W))/ΣD`, which assumes
`A ≈ B`. `_across_pert_budget` is a **centred** sum of squares, so the panel mean *is* the common
mode it centres out: shrinking toward `d̄` collapses `B` while leaving `A`. Measured `(A−B)/A`:

| | unshrunk | toward panel mean, γ=1 | toward control, γ=1 |
|---|---:|---:|---:|
| h1-test | +0.0140 | **+0.0769** | +0.0149 |
| orion-thinned | +0.0119 | **+0.0530** | +0.0089 |

The dropped `q·ε` term is 4–6× larger in exactly the family the shortcut was used to score. Exact
minimisation over `s` (both ρ branches, `B'` recomputed by `_across_pert_budget`), as Δ raw `mse`
against the unshrunk prediction on its own optimum:

| family, γ=1, donor key | h1-test | orion-thinned |
|---|---:|---:|
| toward **panel mean** | −0.0049 | **+0.0152 (worse than not shrinking)** |
| toward **control** (`dĥ' = α⊙dĥ`) | **−0.0218** | **−0.0291** |

Toward-control wins on both harnesses at every γ tested. And `√(W/W_g)` still *rises* for the
panel-mean target (1.381 / 1.172) while the exact `mse` worsens — that divergence **is** the
correction, and it is why the shortcut had to be replaced rather than adjusted.

**A second, independent reason the panel mean was the wrong target.** The argument for it was that
the metric's 0-anchor is the mean perturbation response. But that anchor is the *reference's* mean
effect, and our panel mean does not estimate it: `cos(our d̄, reference mean effect)` = **+0.0162**
(h1-test), **−0.0664** (orion) — capture fraction `cos²` of 0.0 % and 0.4 %, and not even a shared
sign. The common mode is also a *free separable knob* (adding a constant to every row leaves `B`,
`claim` and `ρ` untouched), so folding it into an amplitude family can only cost.

**Corrected value.** Toward-control preserves `ε`, so the `√W` multiplier does transfer. Rule 1
takes the smaller of the two surfaces (1.415 h1, 1.551 orion):

| family, on top of intervention 1 | `ε` | `√W` mult | raw `mse` | `score_mse` | overall |
|---|---:|---:|---:|---:|---:|
| unshrunk | 0.0087 | 1.000 | 0.8785 | 0.1151 | +0.0192 |
| toward panel mean *(as published — retracted)* | 0.0530 | 1.172 | 0.8957 | 0.0970 | +0.0162 |
| **toward control, rule-1 figure** | 0.0149 | **1.415** | **0.8203** | **0.1763** | **+0.0294** |
| toward control, orion figure | 0.0089 | 1.551 | 0.7935 | 0.2046 | +0.0341 |

⚠️ **That `+0.0102` was itself still wrong, and the correction is mine.** It carried the
**blank-row artefact**: rows whose donor key is entirely zero have no defined rank order, and
`argsort` on an all-equal vector yields a full `0…G−1` ramp, so those rows were shrunk to nothing
and scored as if it helped (111/900 rows on the challenge panel). With `α = 1` on blank rows — "no
confidence information, do not shrink" — and the board solve re-anchored on the a11 result, the
honest figure is **`+0.0055` on top of intervention 1**, raw `mse` **0.8569**, which does **not**
pass the leader's 0.8270. Full working, the operating point and a pre-registered prediction:
`A13-OPERATING-POINT.md`.

**Buildability, measured on the real 300-target panel** (exp096-shrinkage,
`mech_shrink_coverage.py`): the agreement key is non-blank on **294/300** targets (6 carried by a
single donor), and a median of **1,988 genes per target** have a non-zero score — so a rank-2000
denominator lands almost exactly where the key stops carrying information, which is a useful
coincidence rather than a tuned choice. ⚠️ One build hazard: a target's own gene is its
top-agreement gene at median rank 0 and 282/294 sit inside rank 2000, but **12 do not** (max rank
10,064), and at γ=1 those would have their on-target knockdown shrunk ~150×. The on-target row is
excluded from all six metrics, so the scored consequence is nil, but a submission whose knockdowns
silently vanish on 12 targets should not ship without a guard.

### 6. Raise the direction cosine — **the only route past 0.85, and the expensive one**

`MSE_min ≈ 1.0237 − 2.077·cos` at the cap. Our `cos = 0.0740`; the leader's is **at least
0.0947** (§(b)). Each `+0.01` of cosine is worth `−0.0208` raw = `+0.0219` `score_mse` =
**`+0.0036` overall**, and it compounds with everything else. It is also the only lever that helps
`pds`, `fid` and `jac` at the same time.

**This is no longer optional if the target is first place.** Interventions 1–5 take us to
`0.878`, or `0.820` with the shrinkage; the leader is at `0.8270` today and the `#1` threshold is the
fastest-climbing rank on the board. Closing the last stretch needs `cos`, and `cos` is the whole
modelling problem — every experiment 001–095. No cheap version is identified here, which is
precisely why interventions 1–5 should ship first: they are what is available *now*.

### 7. A per-gene coefficient `s_g` fitted on a donor panel — **do not build**

Section (c): +1 to +2 % on `√W` out of sample, `≈ +0.003` raw `mse`, `≈ +0.0003` overall, against
a 4–6.5× oracle. The gap is not a tuning problem; a per-gene *regression slope* is a cell-line
property and does not transfer. Note this is a different family from #5, which does not fit a
slope at all — it interpolates between two profiles the model already produces, using an external
confidence score. That distinction is why one works and the other does not.

---

## What would falsify this

Pre-registered, in the order the board can test them:

1. **a10 as built** (`submissions/a10-nb-s015-phi16.vcc`, `VCC_ALPHA` 1.0796) scores raw
   `mse = 1.3567 ± 0.01`. Falsified outside `[1.33, 1.39]`.
2. **a10 rebuilt at `VCC_ALPHA = 0.838`** scores raw `mse = 0.902 ± 0.02` and `score_mse ≈ +0.090`.
   Falsified if `score_mse` stays 0, which would mean either `X` or `C_p` is badly wrong.
3. **Intervention 1** scores raw `mse = 0.878 ± 0.02`, `reach` within `−0.006` scaled of
   Wildflower's, `pds` and `nmae` unmoved beyond their seed noise. Falsified if `reach` falls by
   more than `0.02` scaled — that would mean the top-K key is leaking called genes into the
   over-dispersed tail.
4. **Intervention 5** raises `√W` by `≥ 1.17×`, i.e. raw `mse ≤ 0.854` when combined with #1.
   Falsified if a realised donor-agreement shrinkage moves `mse` by less than 0.01 raw, or if it
   costs more than 0.01 scaled on `pds` — it shrinks toward the panel mean, which `pds` reads.
5. **The tail policy's call-set churn is benign — now MEASURED, so this is a board prediction
   rather than a hope.** `n_pred` rises to ~5417 with coverage still saturated, the new calls sit
   at chance, and `score_fid` falls by `−0.0016 ± 0.004`. Falsified if a realised build moves
   `score_fid` by more than 0.01 scaled, or if `n_pred` lands below `n_real ≈ 4757`.
6. The identity (★) itself is not falsifiable — it is arithmetic — but the *solve* is: any new
   board point of ours must satisfy it to within `±0.01` raw given its measured `A`, `B`,
   `claim`, `Kc`. It already does for a8 (`+0.0044`) and, indirectly, for v15.

## Corrections to the record

* `docs/LOG.md` 2026-08-27 §1 says "`score_mse` and `score_nmae` are per-perturbation normalized
  then averaged — not linear in the reported aggregate. **Do not do arithmetic on them.**" For
  `score_mse` that is wrong in its diagnosis and its prescription: the map *is* affine in the
  reported aggregate, to `2.1e-5` over six teams spanning raw 0.827–0.969. The non-linearity is
  the **per-context clamp**, and it only bites on submissions whose contexts straddle `b`.
  Arithmetic on `score_mse` is fine below `b`.
* `docs/STRATEGY.md` addendum 5 §2 — "`mse` is a solved constraint, not a live one … we are
  already paying the maximum possible price for `mse`, which is to say none" — was right about
  the clamp and wrong about the conclusion. The clamp means further *degradation* is free; it
  never meant the member was unreachable. It is reachable for `+0.019` at a reach cost of
  `−0.001`, and the board's own history says so: Aginglab went `mse` 1.8219 → 0.9594 while
  `reach` went **up** 0.1672 → 0.2674 and `pds` was flat; Ibrahim Mansour went 0.9095 → 0.8703 →
  0.8348 with `reach` rising every time. **There is no measured mse/reach trade on this board.**
* `scripts/analyse_submission.py` inverts `mse` to a cosine via `cos = (r² + 1 − mse)/(2r)`.
  That form omits `Kc`, `claim` and the correction structure entirely — on this panel `claim/ΣD`
  is 0.52–1.08, so the omitted term is larger than the quantity being solved for. Its output is
  not usable; (★) replaces it.
* `docs/RESEARCH-METHOD.md` §9 lists "amplitude (as a global scalar, pre-φ)" as closed. The
  asymmetry note in that same section is exactly right: it was closed against the DE members
  while `mse` sat in the dead zone. Under `mse` it is open, and the optimum is not where the
  earlier work looked for it — it is at the `ρ` cliff, and its position depends on φ.

## Files

| file | what it does |
|---|---|
| `msekit.py` | scorer-exact `expr_mse_unbiased_capped_norm` from cached group sums (bit-exact vs the pipeline) |
| `decompose.py` | the (★) terms and the identity check |
| `mech_chal_stats.py` | streams a 30 GB prediction h5ad into per-group sums + jackknife in ~50 s |
| `mech_ctrl_stats.py` | `b_ctrl`, `C_ctrl` per challenge context |
| `mech_challenge.py` | all computable terms for a shipped `.vcc`, plus the `C_p` proxy |
| `mech_alpha_chal.py` | re-render at any α; exact round-trip at α=1 |
| `mech_alpha_h1.py` | end-to-end α curve on h1-test (re-emit, re-jackknife, re-score) |
| `mech_board_solve.py` | the two-point solve for `X`, `ΣD` and the a8 falsifier |
| `mech_shrinkage.py`, `mech_shrinkage2.py`, `mech_xgene.py` | the (c) families, CV and cross-cell-line |
| `mech_claim_grid.py`, `mech_claim_grid2.py`, `mech_split_policy.py`, `mech_subgate.py` | `claim(spread, φ)` and the split policies |
| `work/*.npz`, `work/*.log` | cached stats and raw run logs |

## A method note worth keeping

Two of this experiment's conclusions were reversed by asking whether a metric *factorises* before
declaring it unmeasurable.

* **RESEARCH-METHOD rule 1** says any metric with `n_real` in its denominator or its union is
  unreadable on our harnesses. True of `fid`'s value. **Not true of `fid`'s precision leg**, which
  is `k/n_pred` and carries no `n_real`. Writing `F_p = (k/n_pred)·min(1, n_pred/n_real)` split the
  member into an unreadable factor and a readable one, and the readable one was the entire risk —
  because on the challenge the unreadable factor is *saturated at 1* and stays there. Rule 1
  should be read as "the metric's **level** is unreadable", not "the metric is unreadable".

  **The negative case matters as much, or the rule gets over-applied.** Of the members rule 1
  names: `fid` factorises cleanly (precision × coverage). `reach` factorises *partly* — the
  purity-versus-normalised-depth curve is readable, which is the "ordering-driven part" rule 1
  already allows, while the `1/n_real` normalisation is not. **`jac` does not factorise at all**:
  `|R ∩ S| / |R ∪ S|` has `n_real` inside the union, added to our own set rather than multiplying
  it, and no rearrangement separates them — which is exactly why only its *sign* is quoted above.
  (`nmae` is not affected by rule 1 in the first place.) The rule earns its keep by telling these
  apart, not by rescuing all of them.
* The same move is what made the whole experiment work. `expr_mse_unbiased_capped_norm` looks like
  it needs the reference. Written as `(A − 2X + Kc − ρ·claim)/ΣD`, four of its six terms are
  computable from our own file plus the control cells we already hold, and only `X` and `ΣD`
  needed the board.

And one about competitors' published work: a rival's diagnosis is evidence about **their model**,
not about the task. Zhou Li's "callability is set by control expression" (107 CPM called vs 38.7
uncalled) is probably true of their pipeline and is the exact opposite of ours (48.1 called vs 87.6
uncalled), because our callability is set by our own effect sizes and spread policy. Test a
competitor's mechanism on our own artefact before building on it.

Finally, a rule that generalises past this member: **a per-gene regression slope is a CELL-LINE
property and transfers at 1.02×; a per-gene confidence score is a MODEL property and transfers at
1.17×.** That is why fitting `s_g = V_g/U_g` on a donor is worthless and why Miller's
`α = 1 − p_DEG` is not. When deciding whether a per-gene quantity is worth building, ask which of
the two it is.

## On the two harnesses bracketing the truth

RESEARCH-METHOD §1 says to run every arm on both harnesses because "their disagreement is your
error bar". This experiment tried to promote that to "they bracket the truth" and **the promotion
does not survive its own data.** Every cross-surface quantity here with a board value to compare:

| quantity | h1-test | orion-thinned | board | brackets? |
|---|---:|---:|---:|---|
| a10's hit rate into the reference's DE set | 0.0996 | 0.0123 | 0.0585 | **yes** |
| `expr_mse_unbiased_capped_norm` (the level) | 2.0551 | 1.9405 | 1.3627 | no — both high |
| plug-in cosine `X/√(A·R)` | 0.04316 | 0.01759 | 0.05105 | no — both low |
| deconfounded cosine `X/√(A·ΣD)` | 0.05192 | 0.03252 | 0.07402 | no — both low |

One for, three against. My earlier framing — "they bracket dimensionless quantities, not levels" —
is **retracted**: the cosine is dimensionless and does not bracket. (Retracted also: I first
compared h1's and orion's *plug-in* cosines against the board's *deconfounded* one. The
like-for-like pairs are in the table; the gap is smaller than I quoted but the sign is the same.)

**What actually separates them.** Both failures err in the *same direction on both surfaces* —
the harness thinks our prediction is worse than the board does, and higher `mse` and lower cosine
are the identical error. Two estimators sharing a bias cannot straddle the truth, however
dimensionless they are. So the criterion is not dimensionlessness:

> **A quantity can bracket only if the two surfaces' errors have opposite sign.** Before
> measuring, ask what systematic both harnesses share, and whether it pushes this quantity the
> same way on both. If it does, read the pair as a **bound**, not a bracket.

**What the shared systematic is — measured, not assumed.** The natural guess is reference sampling
noise: both harnesses have a noisier reference, and noise attenuates measured agreement.
`mech_noise_bias.py` tests it by thinning h1-test's reference and re-measuring:

| reference cells/pert | `R` | `K` | `ΣD` | plug-in cos | deconfounded cos |
|---:|---:|---:|---:|---:|---:|
| 400 | 0.2215 | 0.0670 | 0.1531 | 0.04316 | 0.05192 |
| 200 | 0.2873 | 0.1343 | 0.1516 | 0.03722 | 0.05124 |
| 100 | 0.4257 | 0.2691 | 0.1551 | 0.03158 | 0.05232 |
| 50 | 0.6908 | 0.5393 | 0.1501 | 0.02624 | 0.05630 |

**The plug-in cosine attenuates by 39 % over an 8× thinning; the deconfounded one does not move at
all** (±5 %, no trend), and `ΣD` is flat to 3 %. The metric's own correction does exactly what it
is designed to do. So reference noise explains the plug-in row of the first table and **not** the
deconfounded row — the tempting mechanism is right for one quantity and wrong for the other, which
is why it had to be measured.

**The shared systematic that does explain it is provenance, and it was mine to flag.** The harness
predictions come from `experiments/031-somatic/run.py` — donor transfer plus whitening, with no
PerturbMatch and no CD4. The board numbers come from `scripts/make_submission_pm_cd4.py`. **They
are different models**, and the production one is built for the challenge's contexts while the
harnesses are held-out cell lines. Both harnesses should therefore read a lower cosine, for a
modelling reason rather than a measurement one — and that is checkable *before* the run, by asking
whether the harness is executing the same model as the board.

**Registered as a hypothesis, not adopted as a prior** (n = 4, and RESEARCH-METHOD §2 has already
paid once for a pattern fitted to too few points). The test: the next quantity measured on both
surfaces gets its bracket-or-bound status **predicted in writing first**. Two are already due from
the realisation of intervention 1, and both are predicted **bound, pessimistic**:

* `Δ` precision from the tail policy. The shared systematic is precision headroom above chance —
  +0.0110 (h1), +0.0565 (orion), +0.0041 (board) — and the closed form says the cost is
  proportional to it. Both harnesses have more to lose, so both should read `Δscore_fid` **more
  negative** than the board's `−0.0016`. Falsified if either lands above it.
* The `n_pred` inflation ratio. Already observed at 1.479 (h1) and 1.319 (orion) against the
  challenge's measured 1.135 — both high, consistent with bound-not-bracket, and recorded here
  because it was measured before the rule was written rather than after.


## The same-artifact test: attempted, and it is ill-posed

The competitive-intel session's criterion for valid bracketing evidence is right in principle —
*only the identical prediction artifact scored on every surface measures the surfaces rather than
the models* — and it makes the table above `n = 0` rather than 1-for-3, because model and surface
vary together in every row. My "dimensionless quantities bracket" rule is not refuted; it is
**untested**, and so is the replacement.

The test looked nearly free. `orion-hct116-thinned`'s 132 targets are a **100 % subset** of the
challenge's 300, on a **bit-identical 18,533-gene axis** (h1-test is not usable: 5/85 targets and a
different axis). So the shipped `a8-nb-phi4.vcc` — the board's own Wildflower entry — can be
scored against the orion reference. `mech_same_artifact.py` does it:

| arm | `mse` | `A` | `X` | `ΣD` | deconfounded cos |
|---|---:|---:|---:|---:|---:|
| 031-somatic **on orion** | 1.9405 | 0.3145 | 0.00814 | 0.19906 | 0.0325 |
| a8 (shipped) ctx A **on orion** | **342.97** | 68.317 | 0.0895 | 0.19906 | 0.0243 |
| a8 (shipped) ctx B **on orion** | 343.20 | 68.463 | 0.1388 | 0.19906 | 0.0376 |
| a8 (shipped) ctx C **on orion** | 306.08 | 61.046 | 0.1227 | 0.19906 | 0.0352 |
| a8 (shipped) **on the board** | 1.3627 | 0.2830 | 0.0188 | 0.22857 | 0.0741 |

**`mse` 343 against 1.36. The test is not measuring the surface; it is measuring the cell line.**
A prediction artifact is *baseline-bound*: `b̂_p` is the challenge context's control profile plus an
effect, so scoring it against a different cell line's reference measures the control gap, which
swamps everything. Confirmed exactly — `A` is a SUM over rows, so the like-for-like comparison
multiplies the one-vector control gap by the 132 rows:

| ctx | `‖b_ctrl(challenge) − b_ctrl(orion)‖²/G` | × 132 | measured `A` | explained |
|---|---:|---:|---:|---:|
| A | 0.51682 | 68.221 | 68.317 | **99.9 %** |
| B | 0.51775 | 68.344 | 68.463 | 99.8 % |
| C | 0.46174 | 60.950 | 61.046 | 99.8 % |

The cell-line baseline is **548×** the per-row perturbation signal. `X` and the cosine are
corrupted by the same offset, so no column of that table is readable.

**So the criterion needs a second clause, and the second clause is what actually binds:**

> Cross-surface scoring of an identical artifact is meaningful only when the surfaces **share a
> baseline** — the same cell line and the same control cells. Otherwise the identical artifact
> measures the distance between cell lines, not between scorers.

No two of our three surfaces share a baseline, so **the bracketing question cannot be settled with
what we hold**, and RESEARCH-METHOD §1's own a8 example does not settle it either: what was run on
both harnesses was the same *emission policy* re-realised against each harness's own controls
(`realize.py` drives `emit_counts` on the harness), not the same artifact. Same policy, different
artifact — which is the confounded case, not the clean one.

Recorded as closed-until-new-evidence. Settling it would need one reference panel split into two
disjoint halves and both halves scored against the same submission — which is exactly what the
metric spec's own replicate anchor `r` does, and is the only construction here that holds the
baseline fixed while varying the surface.
### The trap in the construction that would reopen it — VERIFIED in the primary source

The split-half design already exists in the scorer, and it carries a documented failure mode that
would corrupt exactly this test. Read directly in
`.venv/lib/python3.12/site-packages/cell_eval2/anchor.py` rather than inferred:

* `_disjoint_halves` (imported from `ceiling.py`) and `_assert_disjoint_controls` are both called
  in `_score_one_split`, anchor.py:117-118 — the disjointness is *asserted*, not assumed.
* `control_source` is forced to `"pred"` in `_inner_config`, anchor.py:103, and the module
  docstring says it "is not optional".
* The measured size of getting it wrong, quoted verbatim: shared controls read **"0.5-2.3 %
  optimistic, consistently, on all 12 comparisons"** for `de_lfc_nmae`, and "`ceiling.py` measured
  the same mistake turning a 0.54 split-half reliability into 0.74."

And one line the summary of it did not carry, which makes the warning worse rather than better:

> "under this scheme an upward-biased anchor biases every submission's score **DOWNWARD**,
> uniformly, **with nothing in the output that looks wrong.**"

So the failure mode is not merely "it fabricates the agreement you are testing for" — it is silent
and uniform, which is the class of error this whole experiment has been finding. Anyone building
the split-half bracketing test must force `control_source="pred"` and keep the assertion.

**This also clears my own board solve, which depends on the anchor.** The scale I recovered —
`b = 0.987863`, span `0.950163`, so `r = 0.037700` — is fitted against the shipped anchors, and the
shipped anchors are built with disjoint controls and the assertion enforced. My recovered `r` also
lands inside the spec's independently published range of 0.028–0.045. Had the anchor carried the
shared-control bias, every scaled number in this document would have been uniformly wrong in the
same direction with nothing to show for it.

### One proposed edit to `docs/RESEARCH-METHOD.md`, NOT made here

This experiment was scoped to `experiments/096-mse-floor/`, so I have not edited the standing
rules. The change worth considering, for whoever owns that file:

> **§2 addendum — agreement between two agents is not a check.** Two agents reading the same spec
> and reaching the same conclusion have produced one reading held twice, not independent
> confirmation. In this experiment the first agreement reached was the wrong one: two independent
> analyses of the `mse` member both omitted the `Kc` and `ρ·claim` terms and agreed on a shrink of
> ~0.4×, where the answer is ~0.94×. What broke it was running `cell_eval2`'s own functions, not
> further argument. Convergence between readers should lower confidence relative to a single agent
> who has executed the code.

### Span sensitivity of this document's own numbers

The competitive-intel session's observation — a span bias is one multiplicative constant per
member, so within-member ordering is exactly preserved while cross-member arithmetic is not —
applies to the intervention table here, which weighs an `mse` gain against `reach`, `nmae` and
`fid` costs measured on four different spans. Tested rather than asserted, perturbing every span
by ±10 %:

| | Δ`mse` | Δ`reach` | Δ`nmae` | Δ`fid` | **net** |
|---|---:|---:|---:|---:|---:|
| spans 10 % smaller | +0.02131 | −0.00113 | −0.00078 | −0.00030 | **+0.01911** |
| **as published** | +0.01918 | −0.00101 | −0.00070 | −0.00027 | **+0.01720** |
| spans 10 % larger | +0.01744 | −0.00092 | −0.00064 | −0.00024 | **+0.01563** |

The `mse` gain is 12–19× the sum of the three costs throughout, so **the decision is invariant** —
no plausible anchor error makes intervention 1 not worth doing.

**The projected rank is not invariant, and that is worth stating plainly**, because the board is
dense exactly where we would land (0.17637 / 0.17054 / 0.17017 / 0.16890 — four teams inside
0.008):

| | overall | rank |
|---|---:|---:|
| net at published spans | 0.17156 | **3** |
| spans 10 % smaller | 0.17347 | 3 |
| spans 10 % larger | 0.16999 | **5** |
| net + Miller (orion end) | 0.17594 | 3 |

So the "rank ≈ 4–5" hedge earlier in this document is the honest form and the point estimate of 3
is not. This is the one class of claim here that inherits anchor risk: *whether to do it* does not,
*where it lands us* does.

### Board drift dominates anchor error, and the leader has already run this play

Two things checked against the ten snapshots in `data/leaderboard/`, both raised by the
competitive-intel session and both reproduced here rather than taken.

**1. The leader's entire last jump was this member.** Aginglab.com `GeroAI_v11 → v12`,
`0.18991 → 0.21122`, decomposed over the six published `score_*` fields:

| member | scaled before → after | contributes to overall |
|---|---:|---:|
| **mse** | +0.0410 → +0.1693 | **+0.02139** |
| nmae | +0.1784 → +0.1838 | +0.00090 |
| pds | +0.7082 → +0.7089 | +0.00012 |
| jac | −0.0044 → −0.0042 | +0.00003 |
| fid | +0.0034 → +0.0030 | −0.00008 |
| reach | +0.2128 → +0.2066 | −0.00104 |
| | **total** | **+0.02131** |

`+0.02139` of a `+0.02131` jump is `mse`; every other member nets `−0.00008` combined, and
**`reach` fell rather than rose**. ⚠️ Their `pds_cosine` is also flat across that jump
(0.819621 → 0.819937), but that does NOT establish that their direction was unchanged:
`pds_cosine` is a RANK statistic (§1 of the spec — the normalised rank of the correct match under
cosine distance), and a rank saturates while the bulk cosine `X/√(A·ΣD)` moves. The correct
inference about their direction comes from the cap argument in §(b), not from `pds`. That is this experiment's mechanism — buy the member, pay
nothing in `reach` — observed on someone else's build, from the board alone, with none of our
instrumentation involved. It is the strongest external check the analysis has.

**2. Time-to-submission dominates anchor error in the rank projection.** Measured drift of the
rank-3 threshold: `+0.00541`, `+0.00649`, `+0.00441` per day over the last three intervals
(decelerating, so the slowest is used). Against that, the ±10 % span spread moves our projected
overall by 0.00348 total — **one day of board drift is 1.27× the entire anchor uncertainty.**
Holding a fixed 0.17156 against a uniformly drifting board:

| days to submission | 0 | 1 | 2 | 3 |
|---|---:|---:|---:|---:|
| rank at the slow drift rate | 3 | 6 | 6 | 13 |
| rank at the prior interval's rate | 3 | 6 | 12 | 16 |

⚠️ Linear extrapolation of a noisy rate over a short window, and the `#1` figure in particular is
one team's single jump inside a 0.88-day interval. Treat the 3-day row as an upper bound on decay,
not a forecast. The direction is solid; the magnitude is not.

**Two consequences for the plan.** First, the right response to the span sensitivity above is not
to tighten the anchor estimate — it is to **submit sooner**, which is an argument for shipping the
α move immediately that is independent of the build-risk one. Second, the `#1` threshold is
climbing faster than any other rank, so intervention 1's `+0.0192` closes the gap to *today's*
leader and not necessarily to the one we will face; **intervention 5 (Miller, `+0.0044` to
`+0.0097`) is on the critical path rather than a follow-on**, because it is what still clears a
moving `#1`.

### The best-pds / worst-cosine tension, resolved by measurement

The competitive-intel session flagged a real tension: we hold the board's best `pds_cosine`
(0.8318 vs the leader's 0.8199) while being behind by **≥28 %** on aggregate cosine. There is no
contradiction — one is a normalised *rank* of the correct match, the other is bulk directional
agreement — but the combination predicts something specific: our accuracy should be **thick in
some perturbations and thin in others**, which is what lets a rank statistic stay high while the
bulk cosine stays low. `mech_concentration.py` tests it.

**Confirmed, and more extreme than expected:**

| | h1-test | orion-thinned |
|---|---:|---:|
| aggregate cosine | 0.0519 | 0.0325 |
| per-perturbation cosine, mean ± sd | +0.0351 ± 0.0434 | +0.0167 ± 0.0295 |
| p5 / p50 / p95 | −0.005 / +0.027 / +0.111 | −0.019 / +0.010 / +0.076 |
| perturbations with `cos ≤ 0` | 11.8 % | 28.8 % |
| **top quartile's share of total `X`** | **79.7 %** | **85.0 %** |

A quarter of the perturbations carry four-fifths of all the directional signal we have, and a
quarter to a ninth of them point the wrong way. That is exactly the shape a per-perturbation
confidence weighting exists to exploit.

**Which confidence scores actually see it — and only one does.** Spearman correlation with the
true per-perturbation cosine, and the `√W` multiplier from using each as a weight:

| score (all buildable, no reference) | ρ, h1 / orion | `√W` mult, h1 / orion |
|---|---:|---:|
| `‖dĥ_p‖` (self) | −0.001 / +0.253 | 0.595× / 0.983× |
| uniqueness `min_q (1−cos(dĥ_p, dĥ_q))` (self) | −0.337 / −0.102 | 0.698× / 0.915× |
| **donor agreement (k562_gwps × orion_hek293t)** | **+0.209 / +0.242** | **1.144× / 1.035×** |
| *oracle `cos_p`* (leak, ceiling only) | 1.000 / 1.000 | 1.838× / 1.676× |

Both self-referential scores are useless or actively harmful, and they disagree in sign between
surfaces. **The donor-agreement key is the only score positive on both**, at both the per-gene
level (§ intervention 5: 1.381× / 1.172×) and the per-perturbation level here. Per rule 1 take
the smaller: `1.035×` per-perturbation, against an oracle ceiling of `1.68×`, so the crude
min-max weighting used here captures only a fifth of what the family allows — the functional form
is untuned and that is where the remaining headroom sits.

**This unifies interventions 5 and 6.** They looked like separate lines of work — a calibration
trick and "the whole modelling problem". They are the same object seen twice: our directional
signal is concentrated, the concentration is partly visible through an external donor key, and
that key is the only thing in this experiment that has ever transferred. Raising the cosine and
weighting by confidence are two uses of one signal, and the donor key is the handle on both.

### Donor agreement, or external agreement? A second independent key says the latter

If the donor key is the only thing that transfers, the two readings imply different next steps:
donor *screens* carry the signal (tune that key harder), or *any external prior* does (add keys,
they may be complementary). A second external key that is **not a perturbation screen** separates
them. STRING network embeddings are already pinned in this repo —
`data/external/slim/gene_string_embeddings.v0.3.h5`, 21,688 symbols × 64 features, provenance in
`experiments/087-string-target-kernel/README.md`, which preregistered exactly this ablation and
**was never run**. Coverage here: 85/85 h1-test targets, 131/132 orion.

| per-perturbation score (all buildable) | ρ with true `cos_p`, h1 / orion | `√W` mult, h1 / orion | corr with donor key |
|---|---:|---:|---:|
| donor agreement (k562_gwps × orion_hek293t) | +0.209 / +0.242 | 1.144× / 1.035× | — |
| **kg_density** (STRING neighbourhood density) | **+0.190 / +0.187** | 0.979× / 1.029× | −0.053 / +0.304 |
| kg_consensus (STRING neighbour agreement) | +0.190 / **−0.058** | 0.982× / 0.972× | +0.161 / −0.140 |
| donor + kg_density, averaged | **+0.286** / +0.224 | 1.104× / 1.046× | — |

`kg_density` — how tightly a target gene sits in its STRING neighbourhood, a quantity with no
perturbation data in it at all — correlates with our true per-perturbation directional accuracy at
**ρ ≈ +0.19 on both surfaces**, nearly as strongly as the donor key, and is largely independent of
it (|corr| ≤ 0.30). That looked like two unrelated external priors seeing the same thing.

**⚠️ It partly is not, and the confound was caught before it was built on.** Dense STRING
neighbourhoods belong to well-studied hub genes, which are disproportionately essential, and
essential genes have large effects (TRADE: ~45 genes moved typically, >500 for essentials;
Raju reports per-perturbation coherence vs magnitude at Spearman 0.977). So `kg_density` may be
predicting the **true effect magnitude** rather than *our accuracy* — different objects for a
confidence weight. Partial Spearman with `cos_p`, controlling for the true `‖d_p‖`
(`mech_partial.py`):

| | h1-test | orion-thinned |
|---|---:|---:|
| true `‖d_p‖` vs `cos_p` (the confound's strength) | **+0.739** | +0.104 |
| `kg_density`: raw → **partial** | 0.190 → **−0.118** | 0.187 → **+0.178** |
| `donor agreement`: raw → **partial** | 0.209 → **+0.144** | 0.242 → **+0.229** |
| `‖dĥ_p‖` (self): raw → partial | −0.001 → +0.002 | 0.253 → +0.242 |

**The donor key survives partialling on both surfaces; `kg_density` does not.** On h1-test, where
magnitude predicts accuracy strongly (ρ = 0.739), `kg_density`'s signal is entirely magnitude and
its partial correlation goes *negative*. ⚠️ An earlier draft called that 0.739 "Raju's 0.977
reproduced". It is not: Raju measures *within*-perturbation cell-to-cell coherence against
magnitude, this measures *our model's agreement with truth* against magnitude. Related quantities,
a wide gap, and the direction corroborates while the number does not reproduce. Recorded because
"our number matches the literature" is the shape of agreement that has misled this experiment
twice. On orion, where the confound is
weak (ρ = 0.104), it survives. The two surfaces disagree precisely because the confound's strength
differs between them, which is rule 1 working as intended.

One more candidate for the combiner, noticed in the same table: `‖dĥ_p‖` — our OWN predicted
effect magnitude — partials to **+0.242 on orion**, essentially unchanged from its raw +0.253, so
on that surface it is a surviving accuracy signal and it is free. It reads +0.002 on h1-test, so
per rule 1 it is not usable alone; but it costs nothing to carry into the combiner's candidate set
rather than rediscovering it later.

**Honest reading:** the operative property is *external* agreement rather than donor screens
specifically — `kg_density` does carry something — but its apparent accuracy signal is at least
partly, and possibly entirely, a magnitude proxy, and these two surfaces cannot exclude "entirely".
The donor key is the only score here that predicts *where our predictions are trustworthy* after
effect size is removed. A magnitude prior is still useful — but it belongs in the model as an
effect-size prior, not in the confidence weight, and DepMap Chronos essentiality would be a more
direct instrument for it than a graph embedding.

**But the correlation does not yet convert into score.** Combining the two keys raises ρ on
h1-test (+0.286, better than either alone) and not on orion (+0.224, between them), and the `√W`
multiplier does not reliably improve (1.104× / 1.046× against the donor key's 1.144× / 1.035×).
Per rule 1 the transferable figure is the smaller of each pair: **1.035× → 1.046×**, a real but
marginal gain. `kg_consensus` flips sign between surfaces and should not be used.

The limiting factor is the **functional form, not the signal**: a min-max rescale averaged with
another min-max rescale is the crudest possible combiner, and the oracle ceiling for this family is
1.68×. Two independent keys that both correlate at ρ ≈ 0.2 and are mutually uncorrelated is the
raw material for a much better estimator than a mean of two normalised ranks.

**So the next experiment is specified rather than vague:** fit a per-perturbation confidence
weight with a tuned form, validated cross-cell-line the way §(c) validated the per-gene slope —
but build it on the donor key, and admit `kg_density` only after the magnitude confound is
excluded on a third surface. The keys exist and are pinned; §087's protocol already covers the
STRING half, and it should be run with `‖d_p‖` partialled out from the start.

One process note worth more than the result: the asset for this was pinned, SHA'd and
preregistered by experiment 087 four experiments ago and left unrun. The highest-yield measurement
available at the end of this experiment cost one script because someone had already done the hard
part and stopped.


### A fifth lesson, REGISTERED rather than promoted

Three instances in this experiment have the more elaborate construction losing to the plainer one:
the per-gene *slope* losing to the per-gene *confidence*, KG neighbour *consensus* losing to plain
neighbourhood *density*, and my CPM-keyed tail losing to plain top-K. It is tempting to call that
a lesson — "the more elaborate construction fails to transfer".

**It is not promoted here, and the reason is that it has the exact shape of the rule this
experiment already retracted:** n = 3, fitted post hoc, and cutting in a direction both readers
found satisfying. The bracketing rule was also n = 3, also post hoc, also satisfying, and turned
out to be n = 0 once the confound was found.

It is also not obviously the right factorisation of its own three cases. A competing reading
explains all three: what failed each time was the construction that **conditioned on
surface-specific structure** — neighbour consensus conditions on which neighbours co-respond *in
this dataset*, the CPM-keyed tail on *this context's* control depth, the per-gene slope is a
cell-line property by direct measurement (§(c)). The winners — donor agreement, plain top-K,
per-gene confidence — are all surface-independent. That reading covers the same evidence and
would make it *one* distinction found three times rather than two unrelated heuristics, since it
collapses into the external-versus-self-referential cut already established.

**Three instances cannot separate them, so the adjudication is pre-registered instead.** The
combiner experiment will compare functional forms of varying complexity, and the two readings come
apart there:

* *"elaborate fails"* predicts the **simplest** combiner wins, regardless of what its inputs are.
* *"surface-specific fails"* predicts complexity is **fine as long as every input is external**,
  so a tuned multi-key form should beat a simple one built on a self-referential input.

Write the prediction down before the run and see which holds. Until then neither is a rule.
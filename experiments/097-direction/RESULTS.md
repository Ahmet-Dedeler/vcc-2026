# 097 — what raises `X`? The routes to a higher directional cosine, ranked by measurement

**Headline.** The buildable route that raises the direction is a **per-perturbation confidence
weight built from the plain cosine between the donor screens' effect vectors for that
perturbation**, stacked on the per-gene donor-agreement damping. Measured through the scorer's
exact identity on both harnesses, it multiplies the deconfounded cosine by **1.48x** (h1-test
1.478, orion 1.610, transferable = the smaller), against a per-perturbation oracle ceiling of
**1.97x** and a per-gene oracle of **4.1x**. Mapped onto the challenge panel at the `#247` cap it
takes raw `mse` from a12's **0.8785 to 0.8047**, which is **past the leader's 0.8270** on this
member, **projected** worth **+0.0129 overall on top of a12** and **+0.0185 against the live a11
entry**. (Those are the **realised** multiplier put through the board map; the analytic multiplier
gives 0.8030 / +0.0132 / +0.0188. Realised is quoted throughout — §10 (vii).)

> ⚠️ **SUPERSEDED BY §§9e–9g, and in the direction that matters.** The multiplier above is
> `mse`-only. Once `pds` and the four DE members are priced, **this setting is net negative on
> orion (−0.0020)**, and the arm that survives on both surfaces is the conservative one already
> built: **a13, net +0.0215 (h1-test) / +0.0040 (orion)** on a gate of the actual artifact (§9g).
> The binding member is **`nmae`**, not `mse`, and it was unpriced everywhere in this thread —
> here, in the builder, and in `MSE-MECHANISM.md` — until §9f. ⚠️ §9f's rows labelled "a13" are
> **not** a13: this document's `gamma_gene` and the shipped builder's `--shrink-gamma` name
> different functions (§9g). Read §9g before acting on anything in §§1–9f.

> ⚠️ **Every "overall" number in this document is a projection, not a board measurement.** It goes
> through `score_mse = (0.987863 − raw)/0.950163`, an affine map recovered by 096 from six board
> teams, applied to a raw `mse` that is itself derived by transferring a harness-measured *ratio*
> onto 096's two-point board solve for `X` and `ΣD`. The measurement here is the multiplier on the
> direction; the overall points are arithmetic on top of it. Flagged by the `miller` session and
> worth repeating wherever these numbers are quoted.

**It is realised, not only computed.** Re-rendering the profile, driving the production
`emit_counts`, recomputing the jackknife on the emitted cells and calling the scorer returns
**1.473x (h1-test) / 1.580x (orion)** against the analytic 1.478 / 1.610 — agreement to 0.3 % and
1.9 % — with `claim` flat, no perturbation leaving the `#247` cap, and `(A−B)/B` falling rather
than rising. §6.

Five routes are dead and the numbers that kill them are here rather than in a footnote:

| route | verdict | the number |
|---|---|---|
| (c) shrink toward the grand mean of OUR OWN predicted effects | **dead** | `cos(our mean, the reference's mean) = +0.016 / −0.066`; capture is `cos²`, so 0.03 % / 0.4 % |
| (c′) any buildable estimator of the context's mean response | **dead** | best candidate (the donors' `m_full`, i.e. production's zeroed `BETA`) captures 1.1 % on h1-test and **inverts sign** on orion |
| (d) exploit that the on-target knockdown transfers 3–6x better | **dead by construction** | the scorer drops each row's own target coordinate from **every** member; counting it would raise `X` by **+61.6 % / +55.6 %**, and it is not counted |
| (a) donor re-weighting — global, `basal`, `precision`, or adding donors | **dead**, and it *validates* production's uniform weights | tuned global shares sit on opposite sides of 0.5 and each loses **−8 %** on the other surface; uniform is the crossing point and beats the best single donor by 1.17–1.28x |
| (e) the Ahlmann-Eltze cross-cell-line linear model, zero-shot | **dead** | its differential channel is **0.42x** the incumbent on orion, and with the intercept that carries its headline, **20/20 permuted maps beat the real fit** |

And the one route that was **large, positive and unpriced** is now **priced and closed**: dropping
the whitening operator (`gamma` 0.65 → 0) does read **1.495x / 2.136x** on the direction, the
biggest single-lever number here — and it is **net negative on both surfaces** (−0.0033 / −0.0080)
because `pds` charges 1.3x to 2.7x the directional gain for it. No `gamma` is positive on both.
**§9h — production's 0.65 should not be changed.**

And one correction to the standing frame, which changes how every arm in this project should be
priced: **`MSE_MECHANISM.md`'s master formula drops a term that the board charges 4.8x harder
than h1-test says it does.** See §1.

---

## 0. Measured / derived / assumed

**MEASURED** (scorer's own primitives, via `experiments/096-mse-floor/msekit.py` and
`decompose.py`, on `data/harness/h1-test.h5ad` and `data/harness/orion-hct116-thinned.h5ad`):
every `X`, `A`, `B`, `kappa`, `eps`, `claim`, `sumD` in this document; every cosine, every
Spearman and partial Spearman; the permutation nulls and bootstraps; the donor coverage of the
real 300-target panel; the top-K sign accuracies. `B` is `cell_eval2.metrics.delta.
_across_pert_budget` called directly, recomputed per arm, never assumed.

**MEASURED, end to end** (`realize_mse.py`): the realised `mse` for the recommended arms, with
the profile re-rendered, the production `emit_counts` driven, the delete-1 jackknife recomputed
on the emitted cells and `mse_unbiased_capped` called — so `claim`, `rho` and the `#247` cap all
move under the intervention. See §6.

**DERIVED**: the challenge-panel column of §5. It transfers only the *change* in `(kappa, eps)`
measured on the harnesses and applies it to the challenge's own levels, which come from a shipped
submission (`A`, `B`, `Kc`) plus MSE-MECHANISM's two-point board solve (`X`, `sumD`). Levels do
not transfer from a harness (096 measured harness cosine 0.052/0.033 against the board's 0.074);
ratios are the only thing quoted across surfaces.

**ASSUMED**: (i) that the *ratio* `kappa_arm/kappa_base` transfers from harness to board — this
is the load-bearing assumption of the whole document and its falsifier is in §8; (ii) that the
031-somatic harness model's response to a confidence weight resembles the production model's,
when production additionally carries PerturbMatch and CD4 — flagged as the largest risk in §7;
(iii) MSE-MECHANISM's solved `X = 0.018846` and `sumD = 0.228572`, which have now survived two
out-of-sample board predictions.

**NOT MEASURED, and it is the gap**: the four DE members under these arms. §7 gives a proxy that
points the right way and is not a substitute.

---

## 1. The objective, and a correction to the master formula

Everything is scored through the identity, written at its own optimum over the global amplitude:

```
   MSE(s) - 1 = [ s²A - 2sX + Kc - min(claim, s²B) ] / sumD          B = _across_pert_budget

   minimum on the rho cliff s_c = sqrt(claim/B):

   MSE_min = 1 + [ Kc + claim·eps - 2·sqrt(claim·sumD)·kappa ] / sumD

        kappa = X / sqrt(B · sumD)      the alignment, budget-referenced
        eps   = (A - B) / B             the common-mode share
```

`MSE-MECHANISM.md` uses `MSE_min ≈ 1 + Kc/sumD − 2·cos·sqrt(claim/sumD)` with
`cos = X/sqrt(A·sumD)`. That drops `claim·eps/sumD` entirely. For a **global scalar** the drop is
harmless — `eps` is invariant under `dh → s·dh` — and that is the family the formula was derived
for. For any family that moves the common mode relative to the differential it is not harmless,
and **the error is one-signed and optimistic**.

Size on our own current entry: `eps = 0.008716`, so the dropped term is
`0.24653 × 0.008716 / 0.228572 = 0.0154` of raw `mse`. MSE-MECHANISM's `1.0237 − 2.077·cos`
therefore reads **0.0154 low** at the cap for us and more for anything that inflates `eps`.

**And the penalty is priced by `claim/B`, which is far larger on the board than on either
harness** — measured: **0.183** (h1-test), **0.369** (orion), **0.879** (challenge at the `#247`
cap). a12's entire design is to raise `claim` to the cap, so *a12 makes the common-mode penalty
2.4x more expensive than orion says and 4.8x more expensive than h1-test says.* Any arm whose
harness `eps` rises should be re-priced before it ships.

This was found independently and simultaneously by the `miller` builder session from the other
side (the shrink-target measurement), and confirmed on both paths. Every arm below is reported
with its `eps`.

---

## 2. Where the signal is: the common/differential split

`common_mode.py`. With `w`-weighted means `dhbar`, `dbar`, the numerator splits **exactly and
orthogonally** into a common-mode leg and a differential leg:

| | h1-test | orion-thinned |
|---|---:|---:|
| reference's common-mode share of `R` | 3.75 % | 26.51 % |
| our common-mode share of `A` | 1.40 % | 1.19 % |
| common-mode share of `X` | +0.86 % | **−21.2 %** |
| `cos(dhbar, dbar)` | **+0.0162** | **−0.0664** |
| optimal scale on our own mean `c*` | +0.027 | **−0.460** |

Two consequences.

**A constant offset `u` added to every predicted effect is a separable knob.** It leaves `B`
(a centred sum of squares) untouched, leaves `claim` untouched, so it does not move the cliff or
the optimal amplitude. Its optimum is `u = dbar` and it is worth `Sw·||dbar||²` in the numerator,
independent of whatever the amplitude family is doing. This is the clean statement of route (c).

**Our own predicted mean is not that vector.** Capture is `cos²(v, dbar)`, so our own `dhbar`
captures 0.03 % (h1-test) and 0.4 % (orion), and the two surfaces do not agree on the sign of the
optimal scale. Route (c) as posed is dead on arrival.

---

## 3. Route (c): the ceiling is board-anchored, and it is small

The metric spec settles what the anchor is:
`docs/reference/vcc2026-metrics-FULL.md` line 67 — the baseline `b` **is** "the context's mean
perturbation response … assigned identically to every perturbation and then scored as an ordinary
submission". Its §2 adds that such a submission has identical group profiles, hence `B = 0` and
`rho = 0`, forfeiting the whole correction. Putting both into the identity with `A = X = Rc`:

```
   b = 1 + (Kc - Rc)/sumD   =>   Rc = Kc + (1-b)·sumD = 0.005410 + 0.012137·0.228572 = 0.008185
```

So **a perfect estimator of the challenge's mean perturbation response is worth 0.0358 raw `mse`,
= +0.0377 `score_mse`, = +0.0063 overall, and no more.** That is the entire oracle for route (c).

Cross-check against the harnesses, noise-corrected (the plug-in `Sw||dbar||²` contains the
reference's own sampling noise `K/P + Kc`, and orion has 150 reference cells per perturbation):

| | plug-in | noise floor | corrected | as raw `mse` |
|---|---:|---:|---:|---:|
| h1-test | 0.008314 | 0.002217 | 0.006097 | **+0.0398** |
| orion | 0.180331 | 0.007330 | 0.173001 | +0.8691 |
| **challenge (board-anchored)** | | | **0.008185** | **+0.0358** |

h1-test lands within 11 % of the board figure. **Orion's is 24x the board bound**, which no
property of this kind of data can produce — orion's common mode is a harness artefact (a
systematic offset between its control block and its perturbation block that survives the depth
repair), and orion should not be read on this quantity at all. Recorded because it is exactly the
"a quantity can bracket only if the two surfaces' errors have opposite sign" case: here they do
not even have the same order of magnitude.

### Every buildable candidate for that vector, and none of them works

`cmode_estimators.py`, `cos(v, dbar)` with the capture fraction `cos²`:

| candidate `v` (all buildable, none touches the reference) | h1-test | orion | capture h1 |
|---|---:|---:|---:|
| our own `dhbar` (Miller shrink-to-own-mean) | +0.0162 | −0.0664 | 0.03 % |
| **`m_full(k562_gwps + orion_hek293t)` — production's zeroed `BETA`** | **+0.1049** | **−0.3359** | **1.1 %** |
| `m` restricted to our own targets | +0.1149 | −0.0574 | 1.3 % |
| best single donor | +0.1795 | −0.5495 | 3.2 % |

**Sign inverts between surfaces on every one.** Under rule 1 that is not a signal. Production's
`BETA = 0.0` should stay 0. The donors' *panel mean* — the shared CRISPRi/stress response — is a
cell-line property and does not transfer, which is the same distinction §(c) of MSE-MECHANISM
drew for per-gene slopes.

**Independent confirmation from the other side.** The `miller` session measured the same family
through a different implementation, with `B` recomputed: shrink-toward-mean loses to
shrink-toward-zero on both harnesses at every gamma with both keys, **12 of 12**, and on h1-test
the mean target is worse than not shrinking at all. My exact-identity run reproduces it 4 of 4
(§5). The mechanism is the one in §1: shrinking toward the panel mean is shrinking toward exactly
what `_across_pert_budget` centres out, so it collapses `B` and inflates `eps` from ~0.012 to
0.022–0.031 — and `eps` is charged at `claim/B = 0.879` on the board.

---

## 4. Route (d): the on-target knockdown is worth exactly zero

`cross_donor2.py` measured that the target gene's own knockdown transfers across donors at
`rho = 0.647` while the downstream response transfers at `0.10–0.19`. The brief asked whether
that 3–6x can be exploited. It cannot, and the reason is structural rather than empirical.

**MEASURED** (`stack.py`): if the metric counted each row's own target coordinate, `X` would be
**+61.6 % larger on h1-test and +55.6 % larger on orion**. It does not count it —
`cell_eval2.metrics.delta._exclusion_cols` / `_drop_on_target` (issue #172) remove that
coordinate from both legs of the expression member, and the exclusion is panel-wide across every
other member too (`run._ONTARGET_EXCLUDING_FUNCS` covers `mse_unbiased_capped`,
`de_sig_jaccard` and `de_lfc_nmae`; `metrics.direction._exclude_on_target` covers fidelity and
reach; `metrics.discrimination` takes `exclude_target_gene=True` at `exclusion_scope="panel"`,
which drops all 300 panel genes from `pds` — verified in the primary source by the `miller`
session).

So the one part of the perturbation response that transfers cleanly between cell lines is the one
part the scorer deliberately refuses to pay for. **This is not a route; it is an explanation for
why our cosine is 0.074.** It also means any modelling effort spent on knockdown depth is worth
nothing on this competition's score, on every member.

---

## 5. Route (b): the per-perturbation weight — the finding

### The ceilings first

`perpert.py`, family `dh'_p = z_p·dh_p`, reported as the multiplier on `kappa` (equivalently on
the deconfounded cosine, equivalently on `sqrt(W)`):

| ceiling | h1-test | orion |
|---|---:|---:|
| ORACLE, graded, signed `z ∝ x_p/a_p` | **2.012x** | **1.967x** |
| ORACLE, graded, damping only `z ∝ max(0,x_p)/a_p` | 2.013x | 1.918x |
| ORACLE, hard gate on the wrong-way rows `z = 1[x_p>0]` | **1.070x** | **1.361x** |
| ORACLE, keep the top quartile by `x/a` | 1.588x | 1.650x |

**The brief's framing — "drop the wrong-way perturbations" — is the weak version of this route.**
The rows with `x_p ≤ 0` are 11.8 % / 28.8 % of the panel but carry only **−1.0 % / −17.5 %** of
`X`, so removing all of them with a *perfect* oracle buys 1.07x / 1.36x. Nearly all of the
per-perturbation oracle is in the **grading**, not the dropping. Confirmed from the other
direction by the hard top-fraction gates, none of which exceeds 1.17x at any fraction on either
surface and most of which are below 1.

### The key that works, and why it is new

`donor_cos_p` = the **plain cosine between the two donor screens' specific-effect vectors for
perturbation `p`**. 096's key (`donor096`) was the *sum* of the sign-gated geometric mean
`sqrt(|A|·|B|)·1[sign agree]`, which is a magnitude-weighted **total**; normalising it is the
entire difference. Measured, the gain above 1 roughly triples on h1-test (0.241 against 096's
0.087) and quadruples on orion (0.206 against 0.059).

Correlation with the true per-perturbation cosine, and the multiplier as a graded weight
`z_p = rank01(score)^gamma`:

| score (buildable, no reference) | rho, h1 / orion | **partial** rho, controlling true `‖d_p‖` | best multiplier h1 / orion |
|---|---:|---:|---:|
| `‖dh_p‖` (self) | −0.001 / +0.253 | +0.002 / +0.242 | 1.000 / 1.027 |
| uniqueness (self) | −0.337 / −0.102 | −0.132 / −0.106 | 1.000 / 1.000 |
| `donor096` — 096's key | +0.199 / +0.285 | +0.142 / +0.273 | 1.087 / 1.059 |
| `donor_signfrac` | +0.307 / +0.077 | +0.162 / +0.070 | 1.282 / 1.005 |
| **`donor_cos` — this experiment** | **+0.285 / +0.277** | **+0.091 / +0.263** | **1.241 / 1.206** |

`donor_signfrac` is the trap in that table: best on h1-test and nothing on orion. It is not
recommended and would have been the wrong pick from one surface.

The `gamma` curve is **broad and its peak agrees between the surfaces** — h1-test 1.208 / 1.230 /
1.239 / 1.241 and orion 1.166 / 1.190 / 1.203 / 1.206 at `gamma = 1 / 1.5 / 2 / 3`. The
recommendation is `gamma = 2`, which is interior on both.

⚠️ **The h1-test column of that table is two effects, and the next subsection separates them.**
h1-test has 4 rows with fewer than two donors, which get `donor_cos = 0` and are therefore damped;
orion has none. Over the two-donor rows alone the key reads **1.161 (h1) / 1.203 (orion)**, and the
rest of h1-test's 1.239 is the damping of the four evidence-free rows.

### Confound audit on the key

| check | h1-test | orion | reading |
|---|---:|---:|---|
| `rho(donor_cos, true ‖d_p‖)` | +0.301 | +0.211 | partly a size proxy — but it survives partialling on **both**, which `kg_density` did not |
| `rho(donor_cos, our own ‖dh_p‖)` | −0.123 | −0.012 | **not self-referential** — this is the 096 criterion that separates what transfers from what does not |
| `rho(donor_cos, donor096)` | +0.239 | +0.122 | largely a *different* signal from the one 096 measured |

### It is not selection noise

`perpert2.py`. Permutation null over the **whole procedure** — score shuffled, same grid, same
argmax, 1000 draws:

| score | harness | observed | null p50 | null p95 | p |
|---|---|---:|---:|---:|---:|
| `donor_cos` | h1-test | 1.266 | 1.000 | 1.053 | **0.001** |
| `donor_cos` | orion | 1.206 | 1.000 | 1.028 | **0.000** |
| `donor096` | h1-test | 1.087 | 1.000 | 1.050 | 0.019 |
| `donor_signfrac` | orion | 1.017 | 1.000 | 1.029 | 0.106 (null) |

Bootstrap over perturbations at a fixed `gamma = 3` (no in-sample selection): point 1.241 / 1.206,
p5 **0.832 / 0.996**, `P(mult > 1)` **0.81 / 0.95**. The h1-test lower tail is the honest weak
spot — 85 perturbations is not many — and it is the reason the two-surface agreement, not either
number alone, is the evidence.

### It stacks with the per-gene key

`exact.py`, everything through the exact identity with `B` recomputed:

| arm | `kappa` h1 | M h1 | `eps` h1 | `kappa` orion | M orion | `eps` orion |
|---|---:|---:|---:|---:|---:|---:|
| base (as shipped) | 0.05229 | 1.000 | +0.01424 | 0.03272 | 1.000 | +0.01208 |
| per-PERT `donor_cos^2` | 0.06474 | 1.238 | +0.01339 | 0.03931 | 1.202 | +0.01035 |
| per-GENE agreement^1 → **zero** (midrank, corrected) | — | **1.101** | +0.01317 | — | **1.525** | +0.00960 |
| per-GENE agreement^1 → zero (ordinal ranks, **retracted**) | 0.06148 | 1.176 | +0.01333 | 0.04822 | 1.474 | +0.00983 |
| per-GENE agreement^1 → **mean** | 0.06162 | 1.178 | **+0.03072** | 0.04349 | 1.329 | **+0.02245** |
| **STACK pert² × gene¹ (midrank, corrected)** | — | **1.467** | +0.01273 | — | **1.560** | +0.00877 |
| STACK pert² × gene² (midrank, corrected) | — | 1.540 | +0.01242 | — | 1.703 | +0.00841 |
| ORACLE per-perturbation | 0.10518 | 2.012 | +0.01257 | 0.06434 | 1.967 | +0.00931 |
| ORACLE per-gene slope (unbuildable, 096 §c) | 0.21428 | 4.098 | +0.08303 | 0.23675 | 7.236 | +0.23598 |

Two readings. **The two keys are largely independent**: 1.238 × 1.101 = 1.363 against a measured
1.467 on h1-test (super-multiplicative), and 1.202 × 1.525 = 1.833 against a measured 1.560 on
orion (sub-multiplicative) — so the stack is real on both but the pessimistic surface should be
quoted. **And every recommended arm moves `eps` in the SAFE direction** (it falls), which
is the opposite of the shrink-to-mean arm and is why the §1 correction does not eat the gain.

### On the challenge panel

`board.py`. Transfers only the *change* in `(M, eps)`, pessimistically (the worse surface on each),
onto the challenge's own levels `eps₀ = 0.008716`, `kappa₀ = 0.074428`, `claim = 0.24653` (the
`#247` cap = a12's emission policy), `Kc = 0.005410`, `sumD = 0.228572`.

> **This table is the ANALYTIC mapping** — it carries the analytic `(M, eps)` of the table above
> it, so that every row (including the arms that were never realised) comes from one source. The
> **realised** mapping is in §6 and is the one the headline and §8 quote: it gives the stack
> **0.8047 / +0.0129** against the analytic **0.8030 / +0.0132**. Do not mix rows between the two
> (§10 (vii)).


| arm | M used | Δ`eps` | raw `mse` | `score_mse` | Δ vs **a12** | Δ vs **a11 (live)** |
|---|---:|---:|---:|---:|---:|---:|
| a12 as built (global alpha only) | 1.000 | 0 | 0.8785 | 0.1151 | — | +0.0055 |
| per-pert `donor_cos²`, **no floor**, alone | 1.202 | −0.00085 | 0.8463 | 0.1490 | +0.0056 | +0.0112 |
| per-gene agreement^1 → zero alone (midrank) | 1.101 | −0.00107 | 0.8617 | 0.1328 | +0.0029 | +0.0085 |
| per-pert `donor_cos²`, floor 0.3, alone | 1.220 | −0.00080 | 0.8436 | 0.1518 | +0.0061 | +0.0116 |
| per-gene agreement^1 → **mean** | 1.178 | +0.01648 | 0.8687 | 0.1254 | +0.0017 | +0.0072 |
| STACK pert² × gene¹, no floor | 1.467 | −0.00151 | 0.8047 | 0.1928 | +0.0129 | +0.0185 |
| **STACK pert² floor 0.3 × gene¹ — RECOMMENDED** | **1.478** | −0.00146 | **0.8030** | 0.1946 | **+0.0132** | **+0.0188** |
| STACK pert³ floor 0.5 × gene¹ (grid edge) | 1.500 | −0.00150 | 0.7996 | 0.1982 | +0.0138 | +0.0194 |
| ORACLE per-perturbation | 1.967 | −0.00167 | 0.7272 | 0.2744 | +0.0265 | +0.0321 |

**The leader's `mse` (raw 0.82698) is beaten at `M ≥ 1.333`.** The recommended stack is at 1.478.

### The single-donor policy is worth more than the cosine grading, and it nearly went the wrong way

A target carried by fewer than two donors has no cross-donor cosine at all. There are 4 such rows
on h1-test (`AASS`, `CNN3`, `COL6A1`, `MED23`), **0 on orion**, and 6 of 300 on the challenge.
`singledonor.py` measures what they are:

> **they carry 13.1 % of the total predicted amplitude `A` and +0.5 % of the total `X`**, with
> mean per-perturbation cosine **−0.0012** against **+0.0369** for the other 81 rows.

Multiplier on `kappa` under three policies for those rows:

| policy | g=0.5 | g=1 | g=2 | g=3 |
|---|---:|---:|---:|---:|
| **(A) damp them** (`donor_cos = 0` ranks low) — h1-test | 1.149 | 1.207 | **1.239** | 1.241 |
| (B) keep them at `z = 1` — h1-test | 1.041 | 1.023 | **0.958** | 0.892 |
| (C) two-donor rows only, isolated — h1-test | 1.106 | 1.142 | **1.161** | 1.161 |
| all three policies — orion (132/132 have two donors) | 1.117 | 1.166 | 1.203 | 1.206 |

Row (C) is the clean isolation: **the graded cosine over rows that HAVE two donors is worth
1.161 / 1.203 on its own**, and damping the evidence-free rows supplies the rest of h1-test's
1.239. Policy (B) — which the production builder had adopted, on the principled ground that
penalising a target we know nothing about is a self-referential judgement — reads **0.958**, i.e.
worse than uniform, because it spends full amplitude on rows with no direction while damping rows
that have some. I endorsed (B) before measuring it; the measurement reverses that, and the
classification was wrong: **the availability of cross-donor evidence is a property of the donor
panels, not of our prediction, so it is exactly as external as the cosine.** Corroborated
independently — `n_donors` correlates with the true per-perturbation cosine at rho **+0.269**
(partial +0.239 after removing true effect size) on h1-test, the only surface where it varies.

⚠️ This is measured on **4 rows of 85 on one surface**, and orion cannot see it at all. That is
thin, and it is why the recommendation uses a FLOOR rather than the full damping.

### The key with THREE donors — the form production would actually build

The harness model is a two-donor LODO transfer, so every number above uses a two-donor key.
Production has three donor slots and would use the mean pairwise cosine. A third donor can be
added to the KEY without touching the prediction, which is the right directional test
(`threedonor.py`). h1-test may use `orion_hct116` (it is not that harness); orion may not.

| | rho with true `cos_p` | M at `gamma` = 1 / 2 / 3 (floor 0.3) |
|---|---:|---|
| h1-test, 2-donor (as measured throughout) | +0.285 | 1.206 / **1.251** / 1.247 |
| **h1-test, 3-donor pairwise (+`orion_hct116`)** | **+0.350** | 1.140 / **1.183** / 1.207 |
| orion, 2-donor — no leak-free third donor exists | +0.277 | 1.176 / **1.218** / 1.213 |

**The third donor makes the key a better predictor (rho +0.285 → +0.350) and a slightly weaker
weight (gain above 1 falls from 0.251 to 0.183).** Rank correlation and the multiplier are
different functionals, and only one of them is the score. Read this as: the production form works,
and the honest expectation for it is nearer **1.18–1.22** than the 1.25 the two-donor key shows —
so the stacked figure may be optimistic by roughly the same margin. It is measurable on only one
surface, which is why it is a caution rather than a revision.

The other donors cannot supply a third view: coverage of the harness targets is `rpe1` 21/85,
`nadig_jurkat` 21/85, `k562_essential` 18/85 on h1-test and **0/132, 0/132, 0/132** on orion
(`h1_2025` 12/132). That reproduces `CURRENT_STATE`'s "lineage-matched donors: closed on data"
from a different direction and closes the "add a donor" arm on coverage rather than on score.

### The floor is nearly free, and it buys both hazards

`floor.py` / `final_grid.py`, `z_p = max(rank01(donor_cos_p), f)^gamma`, best of a 33-cell grid by
the rule-1 minimum:

| arm | M h1 | M orion | **transferable** | `eps` h1 | `eps` orion |
|---|---:|---:|---:|---:|---:|
| pert g=3 floor 0.5 × gene¹ | 1.500 | 1.658 | 1.500 | +0.01274 | +0.00879 |
| pert g=3 floor 0.3 × gene¹ | 1.494 | 1.561 | 1.494 | +0.01259 | +0.00855 |
| **pert g=2 floor 0.3 × gene¹ — recommended** | **1.478** | **1.610** | **1.478** | +0.01278 | +0.00887 |
| pert g=2 floor 0.2 × gene¹ | 1.473 | 1.580 | 1.473 | +0.01274 | +0.00880 |
| pert g=2 floor 0.0 × gene¹ | 1.467 | 1.560 | 1.467 | +0.01273 | +0.00877 |
| gene¹ ALONE, no per-perturbation weight | 1.101 | 1.525 | 1.101 | +0.01317 | +0.00960 |

**The floor costs nothing on `mse` (1.467 → 1.478) and removes the entire empty-row hazard**, so it
should be taken.

### RETRACTED: "the per-perturbation half carries most of the stack"

An earlier draft of this section read the last row of that table and concluded that the
per-perturbation half carries most of the stack, "the reverse of what `MSE-MECHANISM.md`
intervention 5 would suggest." **That is an h1-test fact stated as a general one, and the `miller`
session caught it.** Which half carries the value *inverts* between the surfaces — realised
numbers, from §6:

| | per-pert alone | per-gene alone | **stack** |
|---|---:|---:|---:|
| h1-test | **1.246** | 1.098 | 1.473 |
| orion | 1.176 | **1.556** | 1.580 |

On h1-test the per-perturbation half's excess over 1 is more than twice the per-gene half's; on
orion it is the reverse and by a wider margin. Neither ordering is a property of the method.

**Two consequences, and the first strengthens the recommendation rather than weakening it.**

1. **Neither half alone is safe to ship.** Rule-1 minimum over the two surfaces, **realised**
   (the row above), with the analytic row alongside as the cross-check and **not** as a source to
   draw from:

   | source | per-gene alone | per-pert alone | **stack** | gap, stack − better half |
   |---|---:|---:|---:|---:|
   | **realised — quote this one** | **1.098** | **1.176** | **1.473** | **0.297** |
   | analytic (cross-check only) | 1.101 | 1.220 | 1.478 | 0.258 |

   The stack is the only unit whose value survives on *both* surfaces, and it is worth more than
   the sum of what either half can be trusted for — precisely because each surface's weak half is
   the other's strong one. Ship the pair or neither.
2. **The stack must not be read as decomposable.** If a later arm drops one half to buy back
   amplitude or call-set, the harness that appeared to justify dropping it may be the wrong one.
   The builder exposes the two as independent flags, both off by default, so that stays a decision
   rather than an assumption.

Note the fourth row against MSE-MECHANISM's intervention 5, which quoted **+0.0044 to +0.0097**
for the donor-agreement shrinkage: computed through the exact identity on the same family and the
same key, the shrink-**toward-mean** form it describes is worth **+0.0017**, and the
shrink-toward-**zero** form it never measured is worth +0.0029. That figure should be treated as
overstated for the family it names.

---

## 6. Realised, not just computed — and it confirms the analytic numbers

Every multiplier above is an analytic evaluation of the identity on a fixed `dh`. `realize_mse.py`
does the whole thing for real: for each arm it re-renders the profile from `b_ctrl + s·dh'`,
drives the **production `emit_counts`** (spread 0.5, phi 4 — a cap-reaching policy so `claim/B`
sits near the challenge's 0.879 rather than the shipped 0.18), recomputes
`jackknife_correction` on the emitted cells, and calls the scorer's own `mse_unbiased_capped`.
`claim`, `rho` and the `#247` cap therefore all move under the intervention instead of being
frozen. It probes at `s = 1` to read the arm's own `claim`, then evaluates its own cliff.

| harness | arm | realised `mse` | Δ vs base | realised `kappa` | **realised M** | analytic M | `eps` | `claim/B` |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| h1-test | base | 0.95371 | — | 0.05123 | 1.000 | 1.000 | +0.01844 | 1.177 |
| h1-test | per-pert `donor_cos²` floor 0.3 | 0.93858 | −0.01512 | 0.06384 | **1.246** | 1.261 | +0.01918 | 1.176 |
| h1-test | per-gene agreement^1 → zero | 0.94747 | −0.00623 | 0.05623 | **1.098** | 1.101 | +0.01804 | 1.178 |
| h1-test | **STACK** | **0.92400** | **−0.02970** | 0.07546 | **1.473** | 1.478 | +0.01781 | 1.180 |
| orion | base | 0.98290 | — | 0.02675 | 1.000 | 1.000 | +0.01684 | 1.186 |
| orion | per-pert `donor_cos²` floor 0.3 | 0.97580 | −0.00710 | 0.03147 | **1.176** | 1.220 | +0.01906 | 1.173 |
| orion | per-gene agreement^1 → zero | 0.95333 | −0.02957 | 0.04162 | **1.556** | 1.525 | +0.01586 | 1.175 |
| orion | **STACK** | **0.95214** | **−0.03076** | 0.04227 | **1.580** | 1.610 | +0.01579 | 1.181 |

**Realised against analytic: 1.473 vs 1.478 (h1-test) and 1.580 vs 1.610 (orion) — agreement to
0.3 % and 1.9 %.** (Two numerical coincidences worth flagging so nobody misreads
the tables. The *analytic* `floor = 0.2` arm in §5 also reads **1.473 / 1.580**; and the *analytic*
`no floor` arm in §5's challenge table also reads **0.8047 / +0.0129 / +0.0185**. Different arms,
same digits, both times. Given §10 (vii) that is worth stating rather than leaving to be
rediscovered — near-agreeing tables are exactly where a row gets taken from the wrong one.) The transferable realised figure is **1.473x**, against the analytic 1.478x. The
individual legs agree to 0.3–3.6 %. Nothing in the arm's value was an artefact of holding the
emission fixed.

**`claim` does not move under the shrinkage — the check all three routes in this thread needed.**
`claim/B` reads 1.173–1.186 across every arm on both harnesses, and on h1-test `claim` is
**0.067008 at every amplitude and for every arm, to six digits, sitting exactly on the `#247`
cap**. Every route here, `miller`'s, and `MSE-MECHANISM.md`'s all held `claim` at its unshrunk
value while the optimum `s* = sqrt(claim/B')` sits directly on it; a concordance between three
analyses that shared an assumption checks nothing (rule 10), and this is the run that could have
broken it.

⚠️ **A flat panel `claim` is necessary but not sufficient, and the `miller` session is right about
why.** `claim` sitting *on* the cap means `min(C_hat_p, K·C_p)` is saturated and set by the
**reference's** jackknife, so a flat total is equally consistent with "nothing moved" and with
"our margin is being eaten and has not run out yet". The failure mode is per-row: the
per-perturbation weight takes some rows to `z = 0.09`, and if any row's `C_hat_p` falls under
`K·C_p` that row leaves the cap and `s*` moves. The direction of the error matters — `claim` would FALL,
`s*` would fall with it, and the arm would be **over**-amplified, which costs `reach` and `nmae`
rather than `mse`, so it would not surface in any of the `mse` accounting in this thread.

**Measured (`capmargin.py`), and the caveat is discharged.** Real emission, real jackknife, each
arm at its own solved cliff:

| harness | arm | solved `s` | `claim` | rows AT the cap | margin `C_hat_p/(K·C_p)` min / p05 / median | low-`z` quartile min |
|---|---|---:|---:|---:|---|---:|
| h1-test | base | 0.5576 | 0.067008 | **85/85** | 1.223 / 1.541 / 1.65 | 1.561 |
| h1-test | STACK | 1.6956 | 0.067008 | **85/85** | **1.222** / 1.535 / 1.65 | 1.560 |
| orion | base | 0.8330 | 0.215616 | **0/132** | 0.173 / 0.388 / 0.47 | 0.173 |
| orion | STACK | 2.3748 | 0.215701 | **0/132** | **0.173** / 0.389 / 0.47 | 0.173 |

**No row changes side, and the whole margin distribution is unchanged to three decimals.**
`C_hat_p` is essentially invariant to the shrinkage row by row, not merely in the panel total.
And the two harnesses happen to bracket the two regimes: h1-test is **saturated** (85/85 at the
cap, so `claim` is set by the reference and cannot respond) while orion is **unsaturated** (0/132,
so `claim` is free to move) — and on orion, where it could have moved, it moved by **+0.04 %**.
That is the stronger of the two tests and it is the one that was available by accident.

**And the challenge sits in the safe regime by construction** — an argument the `miller` session
added and which I had not made. a12's entire design is to raise `claim` to the `#247` cap, i.e.
*into* the saturated h1-test case where rows are pinned and `claim` is set by the reference. So
the concern is smallest exactly where we ship, and orion — the unsaturated surface, where `claim`
was free to move and moved by 0.04 % — is the conservative bound rather than the relevant case.

**And the `miller` session's mechanism prediction beat mine.** I had argued the low-`z` rows were
where the margin was being spent, so the two risks were the same rows. They are not: on h1-test
the low-`z` quartile's minimum margin is **1.560** against the panel minimum of **1.222** — the
most-shrunk rows are among the *safest*. Their reasoning is why: `C_hat_p` is set by the emission
and by the profile's second moment, and every profile is a small perturbation of the same control
composition, so `C_hat_p` is nearly constant in `p` and the margin varies almost entirely through
the **reference's** `C_p`, which has nothing to do with our weights. My "the two risks are not
independent" claim is withdrawn.

### The board mapping, redone from the realised multipliers

§5's challenge table transfers the **analytic** `(M, eps)`. Redoing it from the realised ones, so
the quoted overall points and the measurement come from the same source:

| arm (realised `M`) | M used | Δ`eps` | raw `mse` | `score_mse` | Δ vs a12 |
|---|---:|---:|---:|---:|---:|
| a12 as built | 1.000 | 0 | 0.8785 | 0.1151 | — |
| per-pert `donor_cos²` floor 0.3 alone | 1.176 | **+0.00222** | 0.8537 | 0.1412 | +0.0044 |
| per-gene agreement^1 → zero alone | 1.098 | −0.00040 | 0.8629 | 0.1315 | +0.0027 |
| **STACK — RECOMMENDED** | **1.473** | −0.00063 | **0.8047** | 0.1928 | **+0.0129** |

Against the analytic mapping's +0.0061 / +0.0029 / +0.0132. The stack moves by 0.0003; the
per-perturbation half alone moves by 0.0017, because **its realised `eps` RISES** (+0.00222) where
the analytic one fell. That is the one place the two sources disagree qualitatively, and it is a
second reason not to ship that half alone: the `claim·eps` penalty of §1 works against it and only
the stack keeps `eps` below the base on both surfaces.

Also confirmed by the same run: `eps` **falls** under every recommended arm on both surfaces
(h1-test +0.01844 → +0.01781, orion +0.01684 → +0.01579), so the §1 correction does not eat the
gain — the opposite of the shrink-to-mean family, whose `eps` roughly doubles.

---

## 7. What this does NOT measure, and the two real risks

**Risk 1 — the DE members are unpriced.** Every number above is the expression member only. The
per-gene damping changes which genes clear the significance gate, and CURRENT_STATE item 6 is
explicit that a policy which cannot be produced by an actual emitter is not a result. One proxy,
which points the right way (`buildable.py`): sign accuracy of the top-K genes by `|dh|`, base →
per-gene^1,

| | K=50 | K=100 | K=300 | K=1000 |
|---|---:|---:|---:|---:|
| h1-test base | 0.5207 | 0.5225 | 0.5148 | 0.5174 |
| h1-test per-gene^1 (midrank) | 0.5224 | 0.5229 | **0.5206** | **0.5211** |
| orion base | 0.5520 | 0.5370 | 0.5190 | 0.5092 |
| orion per-gene^1 (midrank) | **0.5611** | **0.5452** | **0.5275** | **0.5197** |

Positive on both surfaces, at every K from 50 up. The per-perturbation weight leaves these
**bit-identical to the base**, which is also the mechanism check that the arm did what it claims:
a per-row scalar cannot reorder genes within a row (`de_proxy2.py`). But this is `|dh|` ordering, and exp095 measured that the realised p-value order
diverges from it badly — so this is suggestive, not priced. `realize.py` is still required.

**Risk 2 — the harness model is not the production model.** The predictions scored here come from
`experiments/031-somatic/run.py` (donor transfer + union whitening, two donors, no PerturbMatch,
no CD4). The board's come from `scripts/make_submission_pm_cd4.py`. `donor_cos` measures where the
*donor-transfer* prediction is trustworthy; PerturbMatch replaces half the K562 slot and CD4 adds
a residual, so the key may be weaker against the production prediction. It may also be *stronger*:
PerturbMatch is a third independent view of the same perturbation, and a four-way agreement key
(K562, PM, HCT116, HEK293T) is strictly more information than the two-way one measured here.
Untested either way.

**A build hazard, measured and then closed.** Holding the cliff after a per-perturbation weight
pushes the global amplitude up by `1/sqrt(mean z²)`: **1.70x at gamma = 1** and **2.15x at
gamma = 2** with no floor, so `VCC_ALPHA` 1.012 becomes ~1.7–2.2 and the bottom-ranked
perturbations emit essentially the control profile (15 targets with `z < 0.05` at gamma = 1, **66
at gamma = 2**, counted on the real panel by the `miller` session). Two things then break:
`make_submission_pm_cd4.py:1285` raises on an empty row, and — the larger one — `pds_cosine`'s own
docstring records that a **zero-norm predicted effect ties its whole row at distance 1.0 and scores
exactly 0.5** under the midrank tie policy. At `pds` raw 0.8333 and span 0.452 that is **0.00041
overall per affected target**, so 66 of them would cost `−0.027` — more than the entire `mse` gain.

Two mitigations, and together they close it. First, `pds` is **scale-invariant**: shrinking a
row's effect does not move its cosine distances at all, so the hazard is not "small z" but
"z small enough that the emitted counts quantise the direction away". a11's board result is direct
evidence that moderate rescaling is safe — `pds` **rose** `+0.0033` under a 22 % amplitude cut.
Second, the floor: with `z_p ≥ 0.3^2 = 0.09`, no row emits the control profile, and it costs
nothing (§5). The realisation in §6 is also the quantisation check — if quantisation bites, the
realised `M` for the per-perturbation arm will fall below the analytic 1.22.
`make_submission_pm_cd4.py:1285` raises on empty rows, and short of raising, a row with no calls
is an `n_pred` loss. **Floor the weight** (`z_p ← max(z_p, 0.2)` or similar); a floored variant is
included in the realisation run so the floor has a measured price rather than a guessed one.

### Buildability on the real panel — measured, not assumed

`buildable.py`, production's three donor slots against the 300 official targets:

| donor | targets carried | panel genes covered |
|---|---:|---:|
| `k562_gwps` | 269 / 300 | 7,679 / 18,533 |
| `orion_hct116` | 290 / 300 | 18,106 / 18,533 |
| `orion_hek293t` | 294 / 300 | 18,106 / 18,533 |

Donors per target: **259 have all three, 35 have two, 6 have one, 0 have none.** So the key
generalises as the mean pairwise cosine and is computable on **294/300**; 6 targets take the floor.
Its distribution on the real panel matches the harnesses (challenge median +0.0284, range
[−0.285, +0.439]; h1-test median +0.021, orion +0.026).

One build detail settled rather than left open: taking each pairwise cosine on the **union**
gene axis (which charges the sparse K562 donor for its zeros) versus on that **pair's shared**
axis gives a rank correlation of **+0.9956**, and `rho(key, n_donors)` is +0.029 / +0.045 either
way — so the key is not secretly counting donors and the axis choice does not matter.

---

## 8. Pre-registered prediction, with falsifiers

**The single best route: the per-perturbation `donor_cos` weight, stacked on the per-gene
donor-agreement damping toward ZERO.**

```
   donor_cos_p = mean pairwise cosine of the donors' specific-effect vectors for p,
                 = 0 when fewer than two donors carry p            (they are then damped)
   z_p         = max( rank01(donor_cos_p), 0.3 ) ** 2              (floor, gamma = 2)
   alpha_pg    = MIDRANK of the donor-agreement vector within row p, / G ; = 1 if that row has
                 fewer than two distinct values                    (shrink toward ZERO, not mean)
   dh'_pg      = z_p * alpha_pg * dh_pg ,  then re-optimise VCC_ALPHA onto the rho cliff
                 (compensation SOLVED from the realisation, not computed from a mean:
                  s rises 0.5576 -> 1.6956 on h1-test and 0.8330 -> 2.3748 on orion, i.e.
                  ~2.9-3.0x, so a12's VCC_ALPHA 1.012 -> ~3.0)
```

**Prediction, on top of a12 (`--keep-top 8000 --keep-spread 0.25 --tail-spread 1.0
--overdispersion 4.0`, `VCC_ALPHA` re-optimised, predicted raw `mse` 0.8785):**

| quantity | predicted | falsified if |
|---|---|---|
| raw `mse` | **0.805 ± 0.03** | outside [0.77, 0.86] |
| `score_mse` | +0.193 | below +0.150 |
| Δ overall vs a12 | **+0.0129** | below +0.006, or negative |
| Δ overall vs a11 (live 0.166115) | **+0.0185** → ≈ **0.1846** | the entry scores below a11 |
| `pds` | within ±0.005 | moves more than 0.01. The basis: `pds_cosine` is scale-invariant, so a graded per-row weight cannot move it except by quantising a row's emitted effect to zero, and the 0.3 floor prevents that; a11's board result corroborates (`pds` **rose** +0.0033 under a 22 % amplitude cut). A large fall would mean quantisation is biting and the floor must go up |
| `n_pred` | within ±15 % of a12's | falls below `n_real ≈ 4757`, which would put `fid`'s coverage leg back in play |

**The mechanism claim that would be falsified**: that `donor_cos` identifies where our prediction
is directionally right. The cheapest falsifier is local and does not cost a slot — build the arm,
score it on **both** harnesses through `realize_mse.py`, and check `kappa` rises by ≥1.35x on
each. If one surface delivers and the other does not, the key is surface-specific and should not
ship, exactly as `donor_signfrac` (1.282 / 1.005) should not.

**What would make me withdraw the recommendation before it is built**: the realised run in §6
showing `claim` moves materially under the shrinkage. Every route in this thread — mine, miller's
and MSE-MECHANISM's — holds `claim` at its unshrunk value while the optimum `s* = sqrt(claim/B)`
sits directly on it, and a concordance between three routes that all made the same assumption
checks nothing (rule 10). §6 is that check.

### Two settings I am NOT recommending, and why

The mse-only objective is **monotone to the edge of the grid** on both harnesses (`gamma = 4`,
`eta = 3` read 1.629 / 1.667 on the ordinal-rank grid of `validate.py` and were still climbing),
because concentrating all the amplitude on the highest-agreement (perturbation, gene) cells
maximises the objective by construction. I am not recommending the edge: it is the setting that
most aggressively shrinks the call set, and the board has already refuted call-set cuts twice
(v14, v15). `gamma = 2`, `eta = 1`, floor 0.3 is interior on both surfaces and is the conservative
choice; `gamma = 3, floor 0.5` reads 1.500 / 1.658 and is worth a further +0.0006 overall, which
is not worth the extra 30-odd near-control rows. The aggressive end is an upper bound, not a plan.

I am also not recommending `eta = 2` on the per-gene half (1.113 / 1.761 midrank). It is the arm
whose two surfaces disagree most, and orion — the surface that likes it — is the one whose
reference has 150 cells per perturbation and whose common mode is a 24x outlier against the board
anchor (§3). When a knob's payoff comes mostly from orion, discount it.

---

## 9. Ranked table of routes

Δcosine is the multiplier `M` on the deconfounded cosine, measured on both harnesses through the
exact identity. "Transferable" is the smaller of the pair (rule 1). Rows 1–3 were also realised
end to end (§6); where both exist the **realised** figure is the one quoted as transferable and
the analytic one is shown in parentheses as its cross-check — never the other way round, and never
mixed within a comparison (§10 (vii)). Rows 4 onward are analytic only.

| # | route | h1-test | orion | **transferable** | oracle ceiling | buildable? | Δ overall vs a12 |
|---|---|---:|---:|---:|---:|---|---:|
| 1 | **(b)+(096§5) STACK: per-pert `donor_cos²` floor 0.3 × per-gene agreement, to ZERO** | 1.478 *(realised 1.473)* | 1.610 *(realised 1.580)* | **1.473x** | 1.97x (pert) / 4.1x (gene) | **yes**, 294/300 coverage | **+0.0129** |
| 2 | (b) per-pert `donor_cos²` floor 0.3 alone | 1.261 *(realised 1.246)* | 1.220 *(1.176)* | **1.176x** | 1.97x | yes | +0.0044 |
| 3 | (096§5) per-gene agreement alone, toward ZERO | 1.101 *(realised 1.098)* | 1.525 *(1.556)* | **1.098x** | 4.1x | yes | +0.0027 |
| 4 | (096§5) per-gene agreement, toward MEAN — *as published* | 1.178 | 1.329 | 1.178x but `eps` **doubles** | — | yes, but loses to the zero target | +0.0017 |
| 5 | (b) 096's `donor096` key, per-perturbation | 1.087 | 1.059 | 1.059x | 1.97x | yes | +0.0018 |
| 6 | (b) drop the wrong-way perturbations (hard gate) | ≤1.165 | ≤1.078 | ≤1.078x | **1.07x / 1.36x oracle** | yes | ≤+0.002 |
| 7 | (c) shrink toward our own grand mean | — | — | **0.0** | +0.0063 overall total | yes, and worthless | ~0 |
| 8 | (c′) donors' mean effect as the offset (`BETA`) | +0.105 cos | **−0.336 cos** | **0.0** | +0.0063 overall total | yes, and sign-inverts | ~0 |
| 9 | (d) exploit the on-target knockdown | — | — | **0.0 by construction** | 0 | n/a | 0 |
| 10 | `donor_signfrac` (per-perturbation) | 1.282 | 1.005 | **1.005x** | 1.97x | yes, and does not transfer | ~0 |
| 11 | `‖dh_p‖`, uniqueness (self-referential keys) | ≤1.000 | ≤1.027 | ≤1.000x | 1.97x | yes, and null | 0 |
| 12 | (e′) drop the whitener (`gamma` 0.65 → 0) | 1.495 | 2.136 | 1.495x on `mse` | — | yes — **but GATED AND CLOSED, §9h** | **−0.0033 / −0.0080** |
| 13 | (a) `normrel^−1.5` per-perturbation donor weighting | 1.043 | 1.035 | 1.035x | — | yes, ~15 lines | +0.0009 |
| 14 | (a) global donor reweighting (`a3`) | argmax 0.60 | argmax 0.35 | **0.918x** | — | yes, and −8 % | negative |
| 15 | (a) `basal` weighting — the exp010 arm, re-tested under `X` | 1.007 | 0.990 | **0.990x** | — | yes, and a no-op | ~0 |
| 16 | (a) add `rpe1` / `nadig_jurkat` / `k562_essential` | — | — | **n/a** | — | **no** — 0/300 challenge coverage | 0 |
| 17 | (e) Ahlmann-Eltze linear, zero-shot between donors | 1.739 | **0.423** | **0.423x** | — | yes, and it loses by 58 % | negative |
| 18 | (e) Ahlmann-Eltze linear, supervised | 3.45 | **0.741** | 0.741x | — | **no** — the challenge ships no perturbation data | n/a |

Rows 12–18 are §§9b–9c. Row 12 is the largest single-lever number in the table and the one I am
**not** recommending on this evidence, because the whitener was fitted for `pds` and `pds` is the
member we lead; see §9c.

### Route (d), the remainder: two more per-gene keys with a mechanism

`extra_keys.py`, multiplier on `kappa`, per-gene weight `rank01(key)^g`:

| per-gene key | h1 ^1 | h1 ^2 | orion ^1 | orion ^2 | transferable |
|---|---:|---:|---:|---:|---:|
| control CPM rank (up-weight well-expressed genes) | 1.058 | 1.067 | 1.152 | 1.243 | **1.058–1.067** |
| control CPM rank, **inverted** | 0.671 | 0.492 | 0.429 | 0.159 | strongly negative |
| `>5 CPM` gate only | 1.006 | 1.007 | 1.055 | 1.060 | 1.006 |
| per-GENE donor cosine, panel-level (cheap) | 1.075 | 1.103 | 1.179 | 1.271 | **1.075–1.103** |
| per-GENE donor sign concordance, panel-level | 1.060 | 1.071 | 1.148 | 1.202 | 1.060–1.071 |

Both are positive on both surfaces and both are much weaker than the per-`(p,g)` agreement key
(1.101 / 1.525). The **panel-level** per-gene donor cosine is worth noting anyway because it is a
single length-`G` vector rather than a `300 × 18,533` matrix, so it is nearly free to build; it
captures roughly half of what the full matrix does. The control-CPM key is the weakest thing here
that is not zero, and its much larger gain on orion (whose reference has 150 cells against
h1-test's 400) is the signature of a low-`n` reference rather than of a modelling gain — treat
1.06 as the estimate, not 1.24.

### The per-gene key at THREE donors — RETRACTION, and the answer

`MultiDonorModel.donor_agreement` requires sign unanimity across **every** carrying donor. My
first reading was that this makes the production key a different, more aggressive estimator than
the one measured here, and that it should be replaced by a mean over donor PAIRS to hold the
graded fraction stable. **The first half is right and the second half is wrong**, both established
by measurement, and the `miller` session supplied the framing that makes it clean:

* With **two** donors, "mean over unordered pairs" has one pair and "unanimity across every
  carrying donor" IS that pair's sign test. The two are the same function, so **both harnesses
  measure both keys and cannot tell them apart.** My "you measured a different estimator" claim
  was wrong; the correct statement is that at three donors they are two different *continuations*
  of the same measured object.
* They continue in **opposite** directions, because unanimity is an INTERSECTION over pairs and
  pairwise is a UNION. Median graded fraction of the 18,533 panel genes, on the real 300 targets:
  **0.106 under unanimity, 0.756 under pairwise**, against the harness's 0.208. My claim that
  pairwise "keeps the graded fraction in band" is **retracted** — it overshoots as far as
  unanimity undershoots.

**h1-test can settle it, and does, because `orion_hct116` is a legal third donor there — which
makes the h1-test triple exactly production's donor set** (`k562_gwps` + `orion_hek293t` +
`orion_hct116`). `pairwise3.py`:

| key | graded fraction | `eta` = 0.5 | 1.0 | 2.0 |
|---|---:|---:|---:|---:|
| 2-donor (= unanimity = pairwise) | 0.208 | 1.070 | 1.101 | 1.113 |
| **3-donor UNANIMITY — production's frozen `donor_agreement`** | **0.109** | 1.162 | **1.278** | 1.378 |
| 3-donor PAIRWISE | 0.761 | 1.134 | 1.213 | 1.321 |

stacked with `z = max(rank01(donor_cos), 0.3)²` at `eta = 1`: **2-donor 1.478x, 3-donor unanimity
1.612x, 3-donor pairwise 1.504x.** The harness graded fractions (0.109, 0.761) reproduce the
challenge panel's (0.106, 0.756) to three decimals, so the two keys behave here as they do there.

**Three conclusions.** (i) Both three-donor forms **beat** the two-donor one, so the production key
is stronger than the estimator measured throughout this document, not weaker — my earlier caution
was wrong in sign. (ii) **Unanimity beats pairwise** and should be kept; `graded_fraction.py`
corroborates on both surfaces independently — holding the key fixed and dialling the graded
fraction from 0.21 down to 0.05 raises the multiplier on both (h1-test 1.101 → 1.138, orion
1.525 → 1.543 at `eta = 1`). Grading fewer genes is mildly better wherever it can be probed.
(iii) The structural objection stands but is outvoted, not refuted: under unanimity a two-donor
target still keeps 0.545 of its genes graded against a three-donor target's 0.106, so the
donor-count ordering is backwards — it just does not cost score on the one surface that can see it.

⚠️ **n = 1 surface, permanently.** orion's only legal third donor covers 12 of 132 targets, so
rule 1's two-surface test cannot be run on this comparison with anything we hold. What justifies
acting on one surface here is that the quantity is not harness-specific: it is production's own
three donor screens scored against a held-out reference.

### The per-gene key is not selection noise either

`genenull.py` — permutation shuffling the per-gene key **across genes within each row**, which is
the null that breaks the claim being made (`validate.py`'s null shuffled across rows instead and
tests something else):

| | observed | null p50 | null p99 | p |
|---|---:|---:|---:|---:|
| h1-test, `eta = 1` | 1.101 | 0.809 | 0.846 | **0.0000** |
| h1-test, `eta = 2` | 1.113 | 0.571 | 0.630 | **0.0000** |
| orion, `eta = 1` | 1.525 | 0.927 | 0.971 | **0.0000** |
| orion, `eta = 2` | 1.761 | 0.738 | 0.816 | **0.0000** |

Note the null sits **below 1**: a random per-gene damping destroys alignment, as it must. The
observed values are 30–95 % above the 99th percentile of that null on both surfaces. (Numbers
are the midrank key of §10; the ordinal-rank version read 1.176 / 1.474 at `eta = 1`.)


---

## 9b. Route (a): donor combination — measured, and it is closed except for one 3.5 % correction

Delegated and measured on both harnesses through the same exact identity, with the cached
whiteners re-used (`wcols` sha256 checked against the cached `.json` on both) and with the
prediction reconstructed **exactly** rather than through a surrogate — the quasi emitter's
exact-total repair means the only stochastic input is `L = int(median(lib_sizes))` per row, and
recovering those integers from the shipped prediction's own group sums reproduces the scorer's
`dh` to `1.8e-15` on 85/85 h1-test rows and 131/132 orion rows.

> ⚠️ **The analytic surrogate `dh = log1p(TS·profile) − b_ctrl` fails, and a per-perturbation
> cosine of 0.998 is not enough to save it.** It biases `sqrt(W)` by **+1.4 % (h1-test) and
> +19.2 % (orion)** — on orion, larger than every arm effect measured in this route. The reason is
> that the aggregate cosine is 0.04, so the signal lives in ~2 % of `dh`'s norm and a 0.2 %
> angular error is ten times the signal. Anyone tempted to approximate `dh` should read that twice.

| arm | h1-test | orion | **cross-surface** | verdict |
|---|---:|---:|---:|---|
| a1 uniform (production's `donor_weights`) | 1.000 | 1.000 | — | the baseline |
| a2 best single donor | 0.857 (k562) | 0.782 (hek293t) | — | **uniform beats it by 1.17–1.28x**, and *which* donor is better flips between surfaces |
| a3 global k562 share, tuned | argmax 0.60 | argmax 0.35 | **0.918 / 0.923** | **dead, −8 %** |
| a4 `basal` weighting | 1.007 | 0.990 | ~0.99 | **still closed**, and now with the mechanism: weights come out 0.531/0.469 and 0.515/0.485 — a near-no-op by construction |
| a5 `precision` (donor cell count) | — | — | 0.989 / 1.019 | **not buildable** — per-perturbation counts are discarded into min/median/max by `build_orion_donor.py:124`; its global proxy disagrees in sign |
| **a6 `normrel^−1.5`** | 1.043 | 1.035 | **1.035x** | **the only arm positive on both** |
| a6 hard per-perturbation donor selection | 0.662 | 0.591 | 0.591 | catastrophic — averaging *is* the mechanism |
| a6 `selfkd` (knockdown depth) | 0.981 | 1.027 | 0.981 | dead cross-surface |
| a7 add a donor | — | — | — | `rpe1`/`nadig_jurkat`/`k562_essential` cover **0/132** orion targets and **0/300** challenge targets. Dead on data. |

**`a3` is the cleanest negative in this experiment.** The two harnesses' tuned donor shares sit on
opposite sides of 0.5 (0.60 vs 0.35), the curve is steep, and each surface's argmax loses ~8 % on
the other — because orion is HCT116 and prefers its same-atlas relative HEK293T while h1-test
prefers Replogle K562. **Uniform is the crossing point.** This *validates* production's
`donor_weights() = {n: 1.0}` rather than challenging it, and it is exactly RESEARCH-METHOD §10's
"a cell line's response is not a model property".

`a6` is worth stating precisely because it is small and real: **down-weight the donor whose effect
for `p` is unusually large relative to its own panel median**, at exponent −1.5. At exponent −1
this is exactly *averaging the donors' unit-normalised effect vectors*, and the optimum sits just
past that. Splitting it: the direction-only leg reads 1.046 / 1.039 and the **magnitude-only leg
reads 1.001 / 0.999 — exactly nil**, which is 096 §(c)'s per-perturbation shrinkage result
reproduced from a third direction. Worth **+0.0008 to +0.0009 overall**; ~15 lines, and
`model.effect(..., weights=…)` already takes a per-target weight dict.

Jackknife over perturbations gives SE ≈ 0.020–0.029 on these multipliers, so `a6` is not
individually significant on either surface (z = +2.2 / +1.4); the evidence is that two surfaces
agree on sign *and* magnitude, 30–80σ outside the emitter noise floor (sd 0.0004 / 0.0012).

## 9c. Route (e): the Ahlmann-Eltze linear model — a clean negative, and one large surprise

Ahlmann-Eltze, Huber & Anders (*Nature Methods* 2025) report that nothing beat the training-set
mean on an L2 objective except a linear model whose perturbation embedding was pretrained on
another cell line. Measured here in both forms.

**(e1) The supervised ceiling does not transfer.** Fit `W` on half of a harness's perturbations
(donor-PCA embedding, ridge, control-centred intercept), score the held-out half: **3.50x on
h1-test, 11.20x on orion** on the raw `kappa`. Both collapse once the perturbation-**constant** is
removed: **3.45x / 0.741x**. The constant's cosine with the panel-average reference is **+0.871**,
so almost all of orion's apparent gain is the panel mean bought for free — a constant raises `X`
and leaves `B` untouched, so it inflates `kappa` at no cost. This arm is unshippable anyway: the
challenge contexts ship 18,400 control cells and **zero perturbation measurements**, so there is
nothing to fit `W` on.

**(e2) The shippable zero-shot form loses.** Fit `W` entirely between donors, apply to the
harness. Best arm reads 1.904 / 2.031 on raw `kappa` and **1.739 / 0.423 on the differential
channel** — a 58 % shortfall against the incumbent on the pessimistic surface. Without the
intercept it reads **−5.0 on orion**: sign-flipped, dead under rule 1. A permutation null
(shuffle the A↔B perturbation pairing, 20 seeds) is decisive: with the intercept, on orion,
**20/20 shuffled maps beat the real fit**. Cross-surface tuning is worse still — the h1-tuned
arm reads **−2.20** on orion, with 293 of 543 configurations flipping sign. **Verdict: no.** The
between-donor map explains 2.6 % of donor B's effect variance, and the paper's construction is
`e1`-shaped where our problem is `e2`-shaped.

**The surprise, and it is bigger than the route it came from: the WHITENER is costing direction.**
Sweeping `gamma` with everything else fixed:

| `gamma` | `kappa` h1 | `kappa` orion | M (emitted) h1 | M (emitted) orion |
|---:|---:|---:|---:|---:|
| **0.00 (no whitening)** | 0.07518 | 0.07958 | **1.495** | **2.136** |
| 0.45 | 0.05755 | 0.04499 | 1.133 | 1.206 |
| **0.65 (shipped)** | 0.05079 | 0.03738 | 1.000 | 1.000 |
| 1.00 | 0.03689 | 0.02526 | 0.732 | 0.679 |

Monotone on both surfaces. **`X` falls to 0.64x (h1) / 0.47x (orion) from `gamma` 0 → 0.65 while
`B` moves only 0.89x / 0.98x — the operator destroys alignment and does not buy budget in return.**
Challenge mapping at the pessimistic pair: raw `mse` **0.808–0.868** against the incumbent's 0.885,
i.e. **+0.003 (strict, common-mode removed) to +0.0136 overall**.

⚠️ **This is mse-only and it must not be acted on from this document.** `gamma = 0.65` was fitted
against `pds_cosine` in experiment 017 and `pds` is the board member we currently lead
(`score_pds` 0.7399, board maximum). exp017's measured `+0.022` macro `pds` for the operator is
worth roughly `+0.004` overall against an `mse` gain of `+0.003` to `+0.0136`, so the net is
genuinely uncertain and could be either sign. It is the largest single unmeasured question this
experiment surfaced, and it needs all six members on both harnesses before anyone touches
`WHITEN_GAMMA`.

---

## 9d. Two things about the EMITTER that bear on `X`, one verified and one small

**`X`, `A` and `B` are exactly invariant to `(spread, phi)` — verified, not cited.**
`MSE-MECHANISM.md` asserts this "by construction" (the quasi emitter repairs each gene's panel
total exactly) and a12's entire plan rests on it: a12 doubles `claim` through emission and assumes
the direction term is untouched. `emit_invariance.py` drives the production `emit_counts` on the
same profile:

| policy | `X` | `A` |
|---|---:|---:|
| `spread 0.25, phi 0` | 0.005072 | 0.322647 |
| `spread 0.50, phi 4` | **0.005072** | **0.322647** |

Bit-identical across a 2x spread change and `phi` 0 → 4. The assumption holds.

**The emitted DEPTH is not invariant, and it is a small real lever.** `emit_counts` uses
`L = int(median(lib_sizes))` and `totals = floor(n*L*profile)`, so the emitted depth decides which
low-expression genes are floored to zero — and a gene emitted as zero where the reference is also
near zero contributes to `X` for reasons that have nothing to do with biology. `depth_sweep.py`,
`kappa` relative to the shipped prediction:

| depth multiplier | 0.4 | 0.6 | 0.8 | 1.0 | 1.3 | 1.7 | 2.5 |
|---|---:|---:|---:|---:|---:|---:|---:|
| h1-test | 0.968 | 0.978 | 0.984 | 0.988 | 0.990 | 0.992 | **0.995** |
| orion | 0.529 | 0.684 | 0.761 | 0.814 | 0.855 | 0.889 | — |

Monotone increasing on both — **emit deeper** — but the challenge-like surface (400 cells, 20k UMI,
like h1-test) moves only **+0.4 % from 1.0x to 1.7x**, and 2.4 % over the whole 4x range. Orion's
much larger sensitivity is its lower depth (12k UMI), which puts far more genes near the floor.
Also measured: the *library-size distribution* matters more than its median — replacing sampled
real control libraries with a uniform library of the same mean moved `X` by **−21.5 % on orion and
−1.2 % on h1-test**. Production already samples real libraries, which is the good end of both. So:
a real mechanism, a **negligible lever on the challenge-like surface**, and a warning that any
reconstruction of a prediction that does not reproduce the library-size distribution will
mis-measure `X` by tens of percent. (That is how this was found: a reconstruction of mine did.)

One consequence for the record: the reference's *perturbation* cells are shallower than its
*control* cells on h1-test (median library 16,818 against 20,190, ratio 0.833) and slightly deeper
on orion (12,874 against 12,174, ratio 1.058). We emit at the control's median. Matching the
perturbation depth instead is a one-line change and, on the evidence above, worth well under 1 %.

---

## 9e. `pds` versus `gamma_gene` — the constraint the mse multiplier cannot see

Asked to confirm or refute independently the `mse-mechanism` result that raw `mse` is monotone in
`gamma_gene` while `pds` improves only to `gamma_gene <= 0.5` and turns negative past ~0.7, on the
strength of one surface (h1-test at three donors).

**Why a two-surface answer was affordable.** `pds_cosine` reads the pseudobulk in `bulk_lognorm`
space, and the quasi emitter's exact-total repair preserves each gene's panel total exactly, so
`pds` is invariant to `(spread, phi)` — asserted by 096, verified bit-exactly for `X`/`A`/`B` by
`emit_invariance.py`, and board-measured at `pds` Δ = 0.000000 under a8's `phi = 4`. At `s = 1` the
render round-trip is exact, so the predicted pseudobulk **is** `b_ctrl + dh'`. No emission needed;
`pds_gamma.py` calls the scorer's own `discrimination_score` on that.

| surface / key | `gg`=0 | 0.25 | 0.5 | 0.75 | 1.0 | 1.5 | 2.0 | 3.0 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| h1-test, 2-donor | 0.8671 | 0.8714 | **0.8720** | 0.8702 | 0.8674 | 0.8624 | 0.8580 | 0.8547 |
| **h1-test, 3-donor (production's key)** | 0.8671 | 0.8750 | 0.8797 | 0.8793 | **0.8815** | 0.8794 | 0.8756 | 0.8696 |
| **orion, 2-donor** | **0.7969** | 0.7950 | 0.7897 | 0.7846 | 0.7790 | 0.7727 | 0.7664 | 0.7612 |

**The specific claim is not confirmed; the concern behind it is, by a stronger route.**

* On h1-test with the **two**-donor key `pds` peaks at `gg = 0.5` and crosses back below baseline
  between 1.0 and 1.5 — at ≈**1.05**, not 0.7.
* On h1-test with the **three**-donor key — production's own — `pds` peaks at `gg = 1.0` and is
  **still above baseline at `gg = 3.0`**. It never turns negative anywhere in the swept range. So
  on the surface the original result was measured on, with the key production actually ships, the
  threshold does not exist.
* **On orion `pds` falls monotonically from `gg = 0`.** There is no positive `gamma_gene` that is
  free. This is the finding that matters, and it is the one a single surface could not produce.

So: `mse-mechanism`'s *threshold* is wrong and its *direction of concern* is right — and right for
a different reason than stated. Under rule 1 the transferable reading is orion's, and it says every
unit of `gamma_gene` is paid for in `pds`.

### `gamma_pert` is free for `pds` — exactly, and this settles half the disagreement

At every `gamma_gene` on every surface, `pds` is **bit-identical across `gamma_pert` = 0, 1 and 2**.
`pds_cosine` is scale-invariant per row and a per-perturbation weight is a per-row scalar, so it
cannot move the metric. **The `gamma_pert` half of the disagreement is therefore not constrained by
`pds` at all** — only `gamma_gene` is. At a13's operating point `gamma_pert = 2` beats
`gamma_pert = 1` on `kappa` on **both** surfaces (h1-test 0.06718 vs 0.06466, orion 0.04406 vs
0.04335), so `mse-mechanism`'s "`gamma_pert` 2.0 is worse than 1.0 at every `gamma_gene`" does not
reproduce here at `gg = 0.25` **on `mse`**. ⚠️ It does not follow that the knob is free, and I drew
that conclusion here and had to retract it: §9f prices `nmae` and `gamma_pert = 2` costs
**−0.0165 overall on orion**, taking the arm net negative. Free on one member is not free — §10
(viii). (It does reproduce on orion at `gg >= 1`, where `gp = 2` falls below
`gp = 1` — so their statement is right in the region they were working in and wrong at a13's.)

> ⚠️ **The floor is load-bearing for `pds`, not just for `n_pred` — but NOT for the rows I first
> named.** A near-zero-norm predicted effect ties its whole distance row at 1.0 and scores
> **exactly 0.5** under the midrank tie policy (`discrimination.py`'s own docstring). Unfloored at
> `gamma_pert = 2` the panel minimum `z` is **2.8e-06**, which is that row. Worth
> `(0.8333 − 0.5)/300` raw per affected target, and the span matters: `/0.452` then `/6` gives
> **−0.00041 overall each**, so six such rows are −0.0025.
>
> **My first version said those rows were the ones with fewer than two donors, and that is wrong**
> (the `miller` session measured it). `rank01` bottoms out at `0.5/n`, never 0, and on the
> challenge panel a `donor_cos` of 0.0 sits at roughly the **45th percentile** (median +0.0284
> over [−0.285, +0.439]), so `--shrink-pert-missing damp` puts the 6 uncarried targets at a **mid**
> rank — `z = 0.111` at `gamma = 2`, not near zero. The rows actually at risk are the
> **bottom-ranked CARRIED** targets, picked out by low cosine rather than by missing donors. On
> h1-test the two populations overlap, because its cosines are mostly positive so 0.0 does rank
> low there — which is exactly why the harness showed the effect through the missing rows. On the
> challenge they are disjoint. This is the same conditioning error as the low-`z`/low-margin one in
> §6, for the third time.
>
> **The distinction changes a decision**, which is why it is here: if the risk looks like it lives
> in 6 targets, `--shrink-pert-missing keep` looks like a fix. It is not — it moves those 6 to
> `z = 1` and leaves every genuinely low-cosine row at 2.8e-06. **The floor is the only fix.**

### The net, and why a13's conservative gamma is the right arm

Transferable (rule-1 minimum) at each config, against the unshrunk base:

| | `M` (mse) | Δ`pds` raw, h1 / orion | Δ overall, mse | Δ overall, `pds` (pessimistic) |
|---|---:|---:|---:|---:|
| **a13** (`gg` 0.25, `gp` 1.0) | 1.237 | +0.0079 / **−0.0019** | +0.0064 | **−0.0007** |
| 097's rec (`gg` 1.0, `gp` 2.0) | 1.513 | +0.0144 / **−0.0179** | +0.0139 | **−0.0066** |

Going from a13's `gamma_gene` to mine buys **+0.0075** on `mse` and gives back **−0.0059** on
`pds` on the surface that cannot be dismissed: **net +0.0016 pessimistic**, +0.0099 optimistic.
That is inside the uncertainty of everything else on the arm, and the four DE members were still
unpriced when it was computed. **a13 as built is the right arm to spend a slot on; 097's
`gamma_gene` is a follow-up, not a correction.** I also proposed `gamma_pert` 1.0 → 2.0 here as a
`pds`-free +0.0013. **§9f refutes that** — it is `pds`-free and `nmae`-expensive, and it takes the
arm negative on orion. Ship a13's `gamma_pert = 1.0`.


---

## 9f. GATE: the four DE members, priced — and they change the recommendation

`de_gate.py`, a13's exact configuration on both harnesses: production `emit_counts` with the
tiered spread vector (top-8,000 gated genes by |effect| at spread 0.25, tail at 1.0, `phi` 4),
scanpy-identical rank-sum against the harness's own controls, BH at 0.05, 40 perturbations, each
arm carrying a **re-seeded twin** as the noise floor. Each arm sits at its own solved cliff
amplitude (`s` = 0.428 / 0.808 / 1.069 / 1.394 on h1-test).

| harness | arm | `n_pred` | reach | reach twin | `fid` | `jac` | `nmae` |
|---|---|---:|---:|---:|---:|---:|---:|
| h1-test | a12 | 3835 | 0.1434 | 0.1227 | 0.5016 | 0.0725 | 0.9630 |
| h1-test | **a13** | 3744 | 0.1442 | 0.1228 | **0.5041** | 0.0745 | **0.9678** |
| h1-test | a13 + `gp`=2 | 3648 | 0.1101 | 0.0982 | 0.5022 | 0.0750 | 0.9768 |
| h1-test | 097's rec | 3537 | 0.1177 | 0.1220 | 0.5016 | 0.0746 | 0.9710 |
| orion | a12 | 3500 | 0.0552 | 0.0797 | 0.5706 | 0.0145 | 0.9189 |
| orion | **a13** | 3516 | 0.0743 | 0.0780 | **0.5742** | 0.0143 | **0.9369** |
| orion | a13 + `gp`=2 | 3456 | 0.0639 | 0.0526 | 0.5783 | 0.0145 | 0.9585 |
| orion | 097's rec | 3364 | 0.0661 | 0.0850 | 0.5817 | 0.0146 | 0.9625 |

**`reach` is unreadable at this sample size and no reach claim should be made from it.** The
same-policy re-seeded twins differ by **0.0207 (h1-test) and 0.0245 (orion)** raw, an order of
magnitude above every between-arm difference in the column. That is the emission-noise floor
exp095 warned about for `phi > 0`, and it dominates. Reported so nobody reads the reach column.

**`n_pred` is safe.** a13 lands within **2.4 % (h1-test) and +0.4 % (orion)** of a12. This is the
term that cost three board predictions (v14, v15, a8) and there is no coverage cliff here.

**`fid` rises, and it is the leg that transfers.** On both harnesses `n_pred` far exceeds
`n_real`, so `fid = k/n_pred` **is** the precision leg — the factor rule 1 says is readable.
a13 gains **+0.0025 / +0.0036** against a seed floor of 0.0001 / 0.0003.

**`nmae` is the cost, and nobody in this thread was pricing it.** a13 degrades it by
**+0.0048 / +0.0180** raw (higher is worse) against a seed floor of 0.0004 / 0.0003, and the
damage grows monotonically with shrinkage: `gp`=2 reads +0.0138 / +0.0396, 097's rec +0.0081 /
+0.0436. The mechanism is direct — per-gene damping toward zero shrinks the predicted fold
changes, and `nmae` is an absolute fold-change error scored over the *reference's* DE genes, which
are not the genes the donors agree on. `jac` disagrees in sign between surfaces and is read as
zero, as its non-factorising union requires.

### The net, member by member

Δ overall against a12, per surface and as the worst-member floor:

| arm | surface | mse | pds | fid | nmae | jac | **net** |
|---|---|---:|---:|---:|---:|---:|---:|
| **a13** | h1-test | +0.0064 | +0.0029 | +0.0013 | −0.0021 | +0.0009 | **+0.0095** |
| **a13** | orion | +0.0088 | −0.0007 | +0.0020 | −0.0075 | −0.0001 | **+0.0025** |
| **a13** | *worst-member floor* | | | | | | **−0.0006** |
| a13 + `gp`=2 | h1-test | +0.0077 | +0.0029 | +0.0006 | −0.0058 | +0.0013 | +0.0068 |
| a13 + `gp`=2 | orion | +0.0094 | −0.0007 | +0.0040 | **−0.0165** | +0.0000 | **−0.0037** |
| 097's rec | h1-test | +0.0139 | +0.0053 | +0.0005 | −0.0033 | +0.0010 | +0.0174 |
| 097's rec | orion | +0.0165 | −0.0066 | +0.0062 | **−0.0182** | +0.0000 | **−0.0020** |

⚠️ **The row labelled "a13" here is 097's power law at `gg = 0.25`, which is NOT a13 as built —
see §9g, which re-runs this gate on the builder's own form and supersedes the a13 rows below.**

**The arm tested here is the only one of these four positive on both surfaces.** Every more aggressive setting goes
negative on orion, and in both cases `nmae` is what does it. The worst-member floor for a13 is
−0.0006, i.e. approximately zero — that construction assumes every member simultaneously takes its
worse surface, which is conservative to the point of being unattainable, but it is the honest floor
and the margin on this arm is thin.

### Two retractions this forces on my own recommendations

**(1) "`gamma_pert` is free" is wrong — it is free on `pds` and it is not free on `nmae`.** §9e
established that `pds` is bit-identical across `gamma_pert`, and I proposed raising it 1.0 → 2.0 as
a free +0.0013. The DE gate refutes that: `gp`=2 costs **−0.0165 overall on orion through `nmae`**
and takes the arm net negative there. **Do not raise `gamma_pert`.** The error was reasoning about
one member's invariance and calling the knob free.

**(2) 097's own recommendation (`gg` 1.0, `gp` 2.0) does not survive the gate.** It is +0.0174 on
h1-test and **−0.0020 on orion**. §9e had already demoted it on `pds`; the DE gate demotes it
further, and for a different reason. Its mse advantage is real and its total is not.

**a13 as built is the arm to spend a slot on.** 097's `gamma_gene` is not a follow-up either,
unless `nmae` is addressed first — which is now the most valuable open question on this route,
because it is the term that eats the entire direction gain.


---

## 9g. RETRACTED: §9f's "a13" row was not a13 — and the real a13 is BETTER

**The two `gamma_gene` symbols in this project name different functions.** Written out:

```
   097, every table in this document:   alpha_g = ( ascending_midrank(agreement)_g / G ) ** gg
   the shipped builder, a13 as built:   alpha_g = exp( -gamma * descending_midrank(agreement)_g / 2000 )
```

A power law on a **fraction** against an exponential on a **count**. Different families: mine is
scale-free in `G`, the builder's carries an absolute 2,000-gene decay length. **No scalar converts
one to the other** — only settings that happen to agree on a chosen summary. Measured
(`scale_check.py`), with the *tie block* being the 87.2 % (h1-test, 3 donors) / 79.4 % (orion) of
each row whose agreement is exactly 0 and which therefore share one alpha:

| form | mean alpha | alpha on the tie block | M h1 / orion |
|---|---:|---:|---:|
| 097 `gg = 0.25` | 0.834 | 0.82 / 0.79 | 1.085 / 1.165 |
| 097 `gg = 1.0` | 0.500 | 0.45 / 0.40 | 1.278 / 1.525 |
| **builder `gamma = 0.25` — a13 as built** | **0.353** | 0.29 / 0.25 | **1.365 / 1.698** |
| builder `gamma = 1.0` — the old default | 0.077 | 0.0067 / 0.0037 | 1.515 / 2.102 |

By kappa: 097's `gg = 0.25` ≈ builder `gamma` **0.04–0.05**; 097's `gg = 1.0` ≈ builder **~0.18**;
builder `gamma = 1.0` ≈ 097's `gg` **2–3**.

**So §9f's row labelled "a13" was 097's power law at `gg = 0.25`, which is not a13.** a13 was built
with the literal flag `--shrink-gamma 0.25` on the shipped builder, so its alpha is the
*exponential* at 0.25 — mean alpha 0.353, not 0.834. ⚠️ Checked and stated because the
`miller` session raised the sharper possibility that this document's two pipelines disagreed with
*each other*: **they do not.** `pds_gamma.py` and `de_gate.py` both build `al` as
`rankdata(ag, "average")/G` and both raise it to `gg`. The multipliers and the DE gate are on one
consistent parameterisation. The defect was a **label**, not a split.

### The real a13, gated

`de_gate.py` re-run with `agreement_shrink_vector`'s own form (including its `<2 distinct values`
guard), same protocol — 40 perturbations, re-seeded twins, each arm at its own solved cliff:

| harness | arm | `n_pred` | `fid` | `jac` | `nmae` |
|---|---|---:|---:|---:|---:|
| h1-test | a12 | 3835 | 0.5016 | 0.0725 | 0.9630 |
| h1-test | builder `gamma`=0.05 | 3736 | 0.5034 | 0.0739 | 0.9672 |
| h1-test | **builder `gamma`=0.25 — a13** | 3539 | 0.4997 | 0.0736 | **0.9566** |
| h1-test | builder `gamma`=1.0 — old default | 3354 | 0.4962 | 0.0728 | 0.9709 |
| orion | a12 | 3500 | 0.5705 | 0.0146 | 0.9191 |
| orion | builder `gamma`=0.05 | 3509 | 0.5740 | 0.0145 | 0.9375 |
| orion | **builder `gamma`=0.25 — a13** | 3401 | 0.5806 | 0.0146 | **0.9538** |
| orion | builder `gamma`=1.0 — old default | 3208 | 0.5877 | 0.0151 | **1.0667** |

Δ overall against a12:

| arm | surface | mse | pds | fid | nmae | jac | **net** | `n_pred` |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| **a13 (`gamma`=0.25)** | h1-test | +0.0150 | +0.0044 | −0.0011 | **+0.0026** | +0.0006 | **+0.0215** | −7.7 % |
| **a13 (`gamma`=0.25)** | orion | +0.0214 | −0.0086 | +0.0057 | **−0.0145** | +0.0000 | **+0.0040** | −2.8 % |
| a13 | *worst-member floor* | | | | | | −0.0092 | |
| `gamma`=1.0 (old default) | h1-test | +0.0198 | −0.0017 | −0.0031 | −0.0034 | +0.0003 | +0.0120 | −12.5 % |
| `gamma`=1.0 (old default) | orion | +0.0313 | −0.0133 | +0.0097 | **−0.0615** | +0.0002 | **−0.0336** | −8.3 % |

**a13 as built is net positive on both surfaces (+0.0215 / +0.0040) and better than the arm §9f
mislabelled as it (+0.0095 / +0.0025).** The mislabelling was conservative, which is luck rather
than design.

**`gamma = 1.0` — the builder's default until tonight — is net −0.0336 on orion**, and its raw
`nmae` of **1.0667** is above the member's own baseline anchor (~1.0013): worse than predicting
nothing. The move to 0.25 was necessary, and the DE gate is what showed it.

**A consistency check worth recording.** 097's power-law `gg = 0.25` (mean alpha 0.834) nets
+0.0095 / +0.0025; the builder's `gamma = 0.05` (mean alpha 0.800) nets +0.0102 / +0.0025. **Two
different functional forms at nearly the same mean alpha land within 0.0007 of each other on both
surfaces** — so the family matters much less than the *amount* of shrinkage, which is the useful
generalisation and is why the mislabelling degraded to a mild arm rather than a wrong one.

⚠️ **`nmae` disagrees in SIGN between surfaces for a13** (+0.0026 on h1-test, −0.0145 on orion).
Under rule 1 its contribution is not reliably signed, which is why a13's worst-member floor
(−0.0092) is worse than the milder arm's (−0.0012) even though its net is better on both. a13 is
the better bet and it is not the safer one.


---

## 9h. GATE: the whitener, priced across all six members — CLOSED, it is load-bearing for `pds`

§9c reported the largest single-lever number in this experiment: dropping the whitening operator
(`gamma` 0.65 → 0) is worth **1.495 / 2.136** on directional accuracy, monotone on both surfaces.
I declined to recommend it because `WHITEN_GAMMA` was fitted for `pds_cosine` in experiment 017
and `pds` is the member we lead. `whiten_gate.py` now prices that trade.

**Protocol.** `experiments/031-somatic/run.py`'s build path reproduced exactly — the same
`MultiDonorModel(gene_norm=False, eb=False)`, the same `wcols` (rebuilt and **sha256-matched**
against each cached operator's meta: 17,745 h1-test / 17,942 orion), the **cached** whitener never
refitted, the same median renormalisation and global `alpha * effect_scale`. a13's shrinkage is
then applied on top in the builder's own form, and every arm is re-optimised onto its **own** rho
cliff before emission. Six members, 40 perturbations, **a re-seeded twin for every arm**.

> **Reconstruction check** (rule 3, and it is not clean on one surface): at `gamma = 0.65` the
> rebuilt `dh` matches the shipped prediction at per-row cosine **0.99922 (h1-test) / 0.99799
> (orion)**, but `kappa` reads **1.0137x / 1.1900x** the shipped value. That is §9b's lesson
> again — the aggregate cosine is ~2 % of `dh`'s norm, so a 0.2 % angular error is 10 % of the
> signal. **Absolute `kappa` levels from this reconstruction are not comparable to the shipped
> ones on orion.** Every number below is a RATIO within one reconstruction, where the same
> pipeline produced every arm, which is the only reading the check supports.

| harness | `gamma` | `kappa` | `pds` | `n_pred` | reach | reach twin | `fid` | `jac` | `nmae` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| h1-test | 0 | **0.10295** | **0.8383** | 3377 | 0.1551 | 0.1475 | 0.5043 | 0.0757 | 0.9657 |
| h1-test | 0.30 | 0.09253 | 0.8677 | 3451 | 0.1432 | 0.1507 | 0.5021 | 0.0748 | 0.9606 |
| h1-test | 0.45 | 0.08749 | 0.8731 | 3488 | 0.1599 | 0.1470 | 0.5011 | 0.0746 | 0.9587 |
| h1-test | **0.65 (ships)** | 0.08174 | 0.8767 | 3542 | 0.1226 | 0.1410 | 0.4994 | 0.0737 | 0.9570 |
| h1-test | 0.85 | 0.07686 | **0.8783** | 3606 | 0.1475 | 0.1186 | 0.4997 | 0.0729 | 0.9570 |
| orion | 0 | **0.07416** | **0.7058** | 3332 | 0.0944 | 0.0863 | 0.5849 | 0.0150 | 0.9414 |
| orion | 0.30 | 0.06629 | 0.7469 | 3366 | 0.0984 | 0.0972 | 0.5815 | 0.0146 | 0.9481 |
| orion | 0.45 | 0.06377 | 0.7621 | 3392 | 0.0884 | 0.1072 | 0.5818 | 0.0147 | 0.9503 |
| orion | **0.65 (ships)** | 0.06127 | 0.7724 | 3414 | 0.0752 | 0.0723 | 0.5804 | 0.0146 | 0.9530 |
| orion | 0.85 | 0.05931 | **0.7791** | 3446 | 0.0997 | 0.0712 | 0.5794 | 0.0145 | 0.9547 |

**`kappa` is monotone DECREASING in `gamma` and `pds` is monotone INCREASING, on both surfaces,
across the whole range.** The operator is exactly the trade experiment 017 bought it for, and both
halves are real.

### Δ overall against the arm that would ship (a13 at `gamma` = 0.65)

| `gamma` | surface | mse | **pds** | fid | nmae | jac | **net** | `n_pred` |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| **0** | h1-test | +0.0108 | **−0.0142** | +0.0027 | −0.0036 | +0.0009 | **−0.0033** | −4.7 % |
| **0** | orion | +0.0090 | **−0.0246** | +0.0025 | +0.0048 | +0.0002 | **−0.0080** | −2.4 % |
| 0.15 | h1-test / orion | +0.0083 / +0.0058 | −0.0069 / −0.0156 | | | | **+0.0022 / −0.0048** | |
| 0.30 | h1-test / orion | +0.0055 / +0.0035 | −0.0033 / −0.0094 | | | | **+0.0027 / −0.0032** | |
| 0.45 | h1-test / orion | +0.0029 / +0.0017 | −0.0013 / −0.0038 | | | | **+0.0023 / −0.0001** | |
| 0.85 | h1-test / orion | −0.0025 / −0.0014 | +0.0006 / +0.0025 | | | | −0.0021 / −0.0002 | |

**No `gamma` is net positive on both surfaces. The route is closed.**

**The exchange rate, stated plainly as asked.** Going 0.65 → 0:

* **h1-test: each +0.01 of `kappa` costs 0.0181 raw `pds`** — +0.0051 overall bought for −0.0067 spent, a ratio of **0.77**.
* **orion: each +0.01 of `kappa` costs 0.0517 raw `pds`** — +0.0070 bought for −0.0191 spent, a ratio of **0.37**.

`pds` alone more than eats the entire directional gain on both surfaces, before any DE member is
counted. Lowering `gamma` partway does not escape it: every setting below 0.65 is positive on
h1-test and negative on orion, with the crossing at `gamma ≈ 0.45` (−0.0001 on orion, inside
noise). Under rule 1 that is a null with no interior optimum — **production's 0.65 is already at
or beside the two-surface optimum and should not be changed.**

**`nmae` is NOT the binding member here** — it improves as `gamma` falls on orion (+0.0048) and
degrades on h1-test (−0.0036), disagreeing in sign, so it is read as zero. `fid` gains a small
amount at every lower `gamma` on both (+0.0025 to +0.0027 at `gamma`=0) and `jac` is flat. **The
whole result is `pds`.** `reach` is again unreadable: twin gaps of 0.0076–0.0289 straddle every
between-arm difference in the column, so it is reported with its twin and used for nothing.

### What this closes, and what it means

This was the largest open number in this experiment and it is now a measured negative. The
whitening operator is not a cost we are paying for `pds` out of ignorance — it is a trade that was
priced correctly in experiment 017 and remains correctly priced under the `mse` objective that
did not exist then. **RESEARCH-METHOD §9's asymmetry note ("a lever closed under one objective is
not closed under a new one") was the right reason to re-open it and the answer came back the same
way.** That is worth recording as a case where the asymmetry rule cost one measurement and
returned a null, which is what it is for.

The consequence for the project is the uncomfortable one the coordinator anticipated: **the
directional gain cannot be bought from the whitener, so the remaining direction has to come from
somewhere else entirely.** Of everything measured in this experiment, what survives is the
donor-agreement family (a13, §9g) at +0.0215 / +0.0040, and its ceiling is the per-perturbation
oracle at 1.97x — most of which is unbuildable. Nothing else here is both positive on two surfaces
and unexploited.


---

## 10. Corrections to my own numbers, made here rather than in a footnote

**(i) The per-gene key was ranked with ORDINAL ranks, and 79 % of every row is tied at zero.**
`np.argsort(np.argsort(v))` breaks ties by gene index. The agreement vector
`sqrt(|A|·|B|)·1[signs agree]` is exactly zero wherever either donor is silent — `k562_gwps`
covers only 7,679 of 18,533 panel genes — or the signs disagree, so a **median 79.2 % (h1-test)
and 79.4 % (orion) of each row is tied at zero**, and my first pass handed all of those genes an
arbitrary ramp ordered by position in the gene axis. Rows with fewer than two donors (4/85 on
h1-test, 0/132 on orion) got a pure index ramp. The same defect was found independently, and
first, in production's `agreement_shrink_vector` by the `miller` session.

Corrected to midranks, with `alpha = 1` on any row carrying fewer than two distinct key values
(`tiefix.py`):

| arm | ordinal (first pass) | **midrank (correct)** |
|---|---|---|
| per-gene ^1, h1 / orion | 1.176 / 1.474 | **1.101 / 1.525** |
| per-gene ^2, h1 / orion | 1.230 / 1.676 | **1.113 / 1.761** |
| STACK pert²×gene¹, h1 / orion | 1.461 / 1.521 | **1.467 / 1.560** |
| STACK pert²×gene², h1 / orion | 1.533 / 1.655 | **1.540 / 1.703** |
| per-pert `donor_cos²` (unaffected) | 1.238 / 1.202 | — |

**The recommendation survives (1.461 → 1.467); the per-gene-alone claim does not** — its
transferable value falls from 1.176x to 1.101x, because the h1-test half of that pair was what the
artefact was inflating. Every headline number above is the midrank one. Note what this says about
harness immunity: the blank-row half of the defect touched 4/85 and 0/132 rows here, so a harness
comparison could never have caught it — the same rule-10 point as the `claim` omission.

**(ii) The amplitude-compensation figure `1.75x` was the `gamma = 1` value under a `gamma = 2`
label.** `buildable.py` printed `mean(rank01²) = 0.327`, which is the compensation for
`z = rank01`, not for `z = rank01²`. Correct: **1.70x at gamma = 1, 2.15x at gamma = 2**. Caught by
the `miller` session cross-checking against its own build. The *arms* measured here were always
`gamma = 2`; only the diagnostic was mislabelled. But the hazard is 25 % larger than quoted, and
the population of near-control rows is **66 targets at gamma = 2 against 15 at gamma = 1** — which
is an argument for `gamma = 1` that I did not have when I wrote §8, and it is why the realisation
in §6 carries both.

**(iii) I endorsed the builder's single-donor policy before measuring it, and it is wrong.**
I wrote to the `miller` session that pinning a target with fewer than two donors at `z = 1` was
"the better choice", on their principled ground that penalising a target we have no evidence about
is a self-referential judgement. Then I measured it: those 4 h1-test rows carry **13.1 % of `A` and
+0.5 % of `X`**, and the policy takes the arm from **1.239x to 0.958x** — below uniform (§5). The
classification was the error: *the availability of cross-donor evidence is a property of the donor
panels, not of our prediction*, so it is exactly as external as the cosine itself, and it is
independently corroborated by `n_donors` (rho +0.269 with the true cosine, partial +0.239). The
recommendation now damps such targets, with a **floor** so they still emit — which is the part the
builder's objection was actually right about. Retracted to them in writing before this document
was finished.

**(iv) `MSE-MECHANISM.md` intervention 5's `+0.0044 to +0.0097` is overstated for the family it
names.** Through the exact identity, the shrink-toward-**mean** form that section describes is
worth **+0.0017** and the shrink-toward-**zero** form it never measured is worth **+0.0029**. The
`sqrt(W) = X²/A` shortcut it uses drops `A − B`, and the panel-mean target is the one family whose
`A − B` is what moves. Established on three independent implementations (mine, `miller`'s, and
`mech_miller.py`'s) which agree on the ordering 12/12 and 4/4 while disagreeing on levels by up to
0.027 raw `mse`.

**(v) "The per-perturbation half carries most of the stack" was an h1-test fact stated as a
general one.** Which half carries the value **inverts** between the surfaces — realised, per-pert
vs per-gene: **1.246 vs 1.098 on h1-test, 1.176 vs 1.556 on orion**. Caught by the `miller`
session. Retracted in full at §5, with the consequence that improves the recommendation: the
rule-1 minimum, **realised**, is 1.098 for the per-gene half alone and 1.176 for the
per-perturbation half alone against **1.473 for the stack**, so the stack is worth more than the
sum of what either half can be trusted for and **neither half alone is safe to ship**. This is the third time in this document
that a claim true on one surface was written as if it were true of the method; the other two are
(i) and (iii) above, and all three were caught by someone else re-running the slice.

**(vii) The "neither half alone is safe" sentence mixed two sources, and it sits where a later
arm will be justified against it.** I wrote the rule-1 minimums as "1.101 per-gene, 1.220 per-pert,
1.473 stack". The first two are **analytic** and the third is **realised**. All-realised they are
**1.098 / 1.176 / 1.473**; all-analytic, **1.101 / 1.220 / 1.478**. Caught by the `miller` session.
The conclusion is unaffected and the all-realised version is *stronger* — the gap between the stack
and the better half widens from 0.258 to **0.297** — but a mixed-source comparison is not a
comparison, and someone deciding whether to drop the per-perturbation half would have read 1.220
and thought it nearly worth shipping alone. §5 now quotes the realised row and labels it, with the
analytic row alongside as a cross-check rather than as a source.

This is the §10 pattern one level down. (i), (iii) and (v) were one-surface facts written as
properties of the method; this is **two tables that agree to 0.3–1.9 % — exactly the condition
under which nobody notices which one a number came from.** The countermeasure is the one already
adopted for the amplitude compensation in (ii): quote the form that cannot be got wrong, and say
which form it is. A tight realised-vs-analytic agreement is the reason both tables exist; it is not
a licence to take whichever number is to hand.

**The two failure modes have different triggers, which is what makes the habit cheap** (the
`miller` session's refinement, adopted verbatim). *Crossing* analytic with realised needs two
**concurrent** sources; *going stale* needs one source and a **revision**. The same grep catches
both, but the prompt to run it differs: the first is triggered by having two tables, the second by
having re-run anything. **"Re-ran a measurement? grep the old digits"** is the cheaper of the two
habits.

**The countermeasure was then tested, twice, and it worked both times.** After this entry was written the
`miller` session grepped its own builder for the same shape and found the superseded analytic
figures (`1.238 / 1.202` alone, `1.461 / 1.521` stacked — my pre-midrank numbers, superseded
twice) still quoted in **four places**. Neither of us found that by re-reading our own text; it
took writing the failure mode down and then searching for it. Both files now carry the realised
figures with the analytic named as cross-check-only. That is a second, independent instance of the
pattern in a different file, so it is a property of how two near-agreeing tables get used rather
than a quirk of this document — and the fix is mechanical (grep for the superseded digits), which
is the useful part. Re-running the grep after that turned up **two more** in the same file, of the
second kind: a flag's help still quoting a pre-midrank-fix measurement, and an `--effect-scale`
help text describing the amplitude cut as "roughly 10x" — a number *computed from a mean* sitting
beside the solved ~2.9–3.0x, which is (ii) again, in a file whose author had just quoted (ii)
approvingly. **Two files, six instances, none found by re-reading.** The sharper statement is that
the author of a number is no better placed than anyone else to spot it going stale, because the
stale value is one they personally believed — the grep works precisely because it does not care
who wrote it.

**(viii) "`gamma_pert` is free" — wrong, and wrong the same way three earlier entries were.**
§9e measured that `pds` is bit-identical across `gamma_pert`, and I proposed raising it 1.0 → 2.0
as a free gain. The DE gate then priced `nmae`: `gamma_pert = 2` costs **−0.0165 overall on
orion** and takes the arm net negative. The error was concluding a knob was free from *one*
member's invariance. That is (i), (iii) and (v)'s shape again — a true statement about one slice
generalised to the whole — this time across members rather than across surfaces.

**(ix) 097's headline recommendation does not survive its own DE gate.** `gamma_gene = 1.0,
`gamma_pert` = 2.0 is +0.0174 on h1-test and **−0.0020 on orion**. The `mse` multiplier that
drove §§1–9 is real and the total is not. The arm that survives both surfaces is the conservative
one the coordinator had already built. **The lesson is the one the whole document was set up to
avoid and did not: an objective that omits a member cannot rank arms that move it**, and `nmae`
was omitted from every calculation here, in `miller`'s builder, and in `MSE-MECHANISM.md`.

**And the sharper form of it, which the `miller` session put better than I had.** The builder's
`gamma = 1.0` default reads raw `nmae` **1.0667** against that member's own baseline anchor of
≈1.0013 — *worse than predicting nothing on that member*. It also carries the **highest
`sqrt(W)` anywhere in this experiment**: 1.515 / 2.102. So the `mse` surface did not merely fail
to see the problem — **it pointed hardest at exactly the worst arm in the family.** "mse-only
accounting is incomplete" understates it; the correct statement is that within this family the
`mse` objective is *anti-correlated* with the total at the aggressive end, so optimising it
without a gate actively selects the arm that loses.

**(x) Two different functions were both called `gamma_gene`, and one of them was mine.** This
document's tables use `alpha = (ascending_midrank/G)**gg`; the shipped builder uses
`alpha = exp(-gamma * descending_midrank / 2000)`. A power law on a fraction against an
exponential on a count — different families, no scalar conversion, and only one of them is
scale-free in `G`. §9f's rows labelled "a13" were mine at `gg = 0.25` (mean alpha 0.834); a13 as
built is the builder's at `gamma = 0.25` (mean alpha **0.353**). §9g re-gates the real artifact.

Two things make this the most dangerous entry in this list. It is the first defect here that would
have **reached a submission** — the coordinator would have built an arm from a gate that never
tested it. And it is invisible to every check the rest of §10 installed: both my pipelines were
internally consistent, the digits were not stale, and no two tables were crossed. **The failure was
a shared NAME across two codebases, which no grep of my own file could find.** The `miller` session
found it by asking for the formula written out. That is the countermeasure — when a number crosses
a codebase boundary, transfer the *definition*, not the value.

The parameterisation-independent quantity, which both sides should have been quoting: **mean
alpha**. My `gg = 0.25` (0.834) and the builder's `gamma = 0.05` (0.800) net within **0.0007** of
each other on both surfaces despite being different functional forms — the amount of shrinkage
transfers, the knob that produces it does not.

**(vi) The "overall" points are projections and are now labelled as such throughout.** They pass
through an affine map fitted on six board teams and a two-point board solve. The measured result
is the multiplier on the direction (analytic 1.478 / 1.610, realised 1.473 / 1.580); everything
downstream of it is arithmetic. Also flagged by the `miller` session.

**(xi) A divergence between the gated artifact and the builder, raised and then MEASURED AWAY.**
The `miller` session read `whiten_gate.py`'s per-perturbation line and found that it overrides an
uncarried target's rank *after* ranking (`where(ndon>=2, rank01, 0.0)` then floor), where the
builder's `damp` sets the SCORE to 0 and lets it rank with the panel. **The reading of my line is
correct.** The conclusion drawn from it is not, and only the measured panel settles it.

On the **real** 300-target panel (`work/challenge_donor_cos.npz` — the three production donors on
the actual targets, not a simulation), **17.0 % of targets have `donor_cos < 0`, so a 0.0 cosine
earns `rank01 = 0.1800`** — not the 0.43 a simulated distribution suggested. Since `0.18 < FLOOR`,
`max(0.18, 0.3)` and `max(0.0, 0.3)` are both **0.3**, and the two policies are **bit-identical on
all 300 rows**:

| `FLOOR` | rows differing | uncarried `z`: gate / builder |
|---:|---:|---|
| 0.00 | 6 | 0.0000 / 0.1800 |
| 0.10 | 6 | 0.1000 / 0.1800 |
| **0.18** | **0** | 0.1800 / 0.1800 |
| **0.30 — a13 ships here** | **0** | **0.3000 / 0.3000** |

**The two policies coincide for any floor ≥ 0.18, and a13's floor is 0.30.** No builder change is
needed and the gated artifact is the shipped one after all.

**The `donor_cos` distribution, labelled once so it stops drifting.** Reconciling the three rows
quoted across this thread turned up **two** axes of difference, not one — the statistic *and* the
population:

| source | statistic | population | values |
|---|---|---|---|
| `buildable.py`, §7, §9e | p05 / p25 / p50 / p75 / p95 | **294 carried** | −0.0428 / +0.0063 / +0.0284 / +0.0563 / +0.1147 |
| this document's range quotes | **min / max** | 294 carried | −0.2848 … +0.4392 |
| my own message of 20 minutes earlier | p05 … p95 | **300, uncarried 0-filled** | −0.0421 / +0.0057 / +0.0270 / +0.0547 / +0.1142 |

All three are correct and none is comparable to another without both labels. The `miller`
session's diagnosis — "yours are p05/p95 and mine are min/max" — is right on the statistic and
missed the population axis, which is the smaller of the two but is **mine**: I 0-filled the six
uncarried targets into a percentile in the very message that was arguing for labelling numbers.
**Canonical form for this document: p05…p95 over the 294 CARRIED targets, with min/max quoted
separately as the key's range.**

Worth keeping for the shape rather than the outcome. The divergence was real in the *code* and
absent in the *artifact*, because a third parameter — the floor — dominated both branches. Two
implementations can differ in a way that is invisible on every input they will ever see, and the
only way to know which case you are in is to run the actual population through both. **"These two
lines differ" and "these two artifacts differ" are separate claims and the second needs the data.**
It also puts a bound worth recording: the floor is doing more work than either policy, so the
`--shrink-pert-missing` choice is unmeasurable at `FLOOR = 0.3` and only becomes a live question
if the floor is ever lowered below 0.18.

### The one rule this list actually earned

Ten entries, and the `miller` session named what unifies the most important three. In **(vii)**
(analytic crossed with realised), **(x)** (two different functions both called `gamma_gene`), and
§9h's `kappa` caveat (a reconstruction reading 1.19x the shipped value on orion), **the fix was
never "get a better number" — it was "say which object the number is from".** Every one of those
three numbers was correct in its own frame and wrong the moment it crossed into another, and no
amount of extra precision would have helped, because precision was not what was missing.

The three failures also sit at increasing distance, which is why they needed three different
catches:

| distance | example | what finds it |
|---|---|---|
| two columns in one table | analytic vs realised, (vii) | grep the superseded digits |
| two codebases | `gamma_gene` vs `--shrink-gamma`, (x) | ask for the formula written out |
| two pipelines for one quantity | reconstructed vs shipped `kappa`, §9h | check the reconstruction against the artifact, and report the ratio even when the cosine looks perfect |
| two panels for one policy | rank of a 0.0 cosine, harness vs challenge, (xi) | state which panel a rank was computed on — and run the real one before believing a divergence |
| **code differs, artifact does not** | `damp` vs the gate: two policies, **zero** differing rows at `FLOOR` 0.30, (xi) | run the REAL population through both branches before believing a divergence |

⚠️ **The last row is not the same failure as the four above it, and the `miller` session was right
to insist it stay distinct.** Those four are all *a number carried into a frame where it did not
mean the same thing*. This one is *two definitions that genuinely differ, over an artifact that
does not care* — the code reading was correct and the conclusion drawn from it was not. The check
is the same shape (go to the data) but it protects against the **opposite** error: there, believing
two numbers agree when they do not; here, believing two implementations diverge when they cannot.

The last is the hardest, because a per-row cosine of **0.99922** reads as a faithful reconstruction
and is not one for a quantity that lives in 2 % of the vector's norm. **A validation statistic
chosen for the wrong scale will pass.** That is the same content as §9b's surrogate warning,
applied to my own pipeline rather than to someone else's — which is where it was hardest to apply
and where it nearly was not.

**So the standing rule, and it is cheap: when a number crosses a boundary, transfer the
definition, not the value.** Where the definition cannot travel — a reconstruction that cannot be
made bit-exact — transfer the *ratio within one frame* and say so, which is what §9h does and why
its verdict survives a 19 % level error.

## 11. Files

| file | what it does |
|---|---|
| `common_mode.py` | the orthogonal common/differential split of `X`, `A`, `R` |
| `cmode_estimators.py` | every buildable candidate for the reference's mean effect; route (c) |
| `perpert.py`, `perpert2.py` | the per-perturbation keys, gamma curves, permutation null, bootstrap |
| `stack.py` | do the two keys stack; the on-target check; the noise-corrected common mode |
| `validate.py` | the 2-D (gamma, eta) grid, cross-surface exponent choice, full-procedure null |
| `exact.py` | every arm through the exact identity with `B` recomputed — `kappa` and `eps` |
| `board.py` | the challenge-panel mapping with the board-anchored constants |
| `buildable.py` | donor coverage of the real 300; the DE-side top-K sign-accuracy proxy |
| `keycheck.py` | confound audit on `donor_cos`; the union-vs-pair-shared gene-axis build detail |
| `extra_keys.py` | the control-CPM and panel-level per-gene keys (route d remainder) |
| `genenull.py` | permutation null shuffling the per-gene key ACROSS genes |
| `tiefix.py` | the ordinal-vs-midrank correction of §10 |
| `realize_mse.py` | end-to-end realisation: re-render, re-emit, re-jackknife, re-score |
| `capmargin.py` | the per-row `#247` cap margin before and after the shrinkage |
| `singledonor.py` | isolates the 4 rows that carried h1-test's per-perturbation result |
| `floor.py`, `final_grid.py` | the (gamma, floor) grid and the settled operating point |
| `threedonor.py`, `pairwise3.py` | the per-perturbation and per-gene keys built with THREE donors |
| `graded_fraction.py`, `chal_alpha.py` | how the per-gene key's graded fraction moves with donors |
| `de_proxy2.py` | the DE-side top-K sign-accuracy proxy, corrected keys |
| `wkit.py`, `w_*.py` | route (a): donor combination (exact-emission reconstruction) |
| `lin_kit.py`, `lin_*.py` | route (e): the Ahlmann-Eltze linear model and the whitener sweep |
| `emit_invariance.py`, `depth_sweep.py` | is `X` invariant to the emission policy and to the emitted depth |
| `work/*.log`, `work/*.npz` | raw run logs and cached intermediates |

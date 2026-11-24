# 227 — `nmae` is NOT a separate axis from cosine. The board evidence for it is a clamp artefact, and Aginglab's edge is a LEVEL, not a direction.

**Date:** 2026-09-04. Measurement only. **No `vcc`, no `cell-eval2`, no build, no submission,
nothing written outside `experiments/227-nmae-axis/`.** Peak RSS **1.05 GB** (`plolo227.py`); every
other script under 0.2 GB. One job at a time; `cell-eval2` checked as not running before starting.

---

## Verdict in six lines

1. **`nmae` cannot be reached at fixed cosine.** Over four sign-fixing orderings — two of them
   *oracle* orderings chosen to be maximally cheap in cosine — the minimum cost of reaching
   Aginglab's `nmae` credit is **+5.0 % relative aggregate cosine** on h1-test and **+12.3 %** on
   PLOLO's three cancer folds. The mandate's gate was ±0.2 %. **Failed by 25–60×.**
2. **Candidate 1 (emission spread/dispersion) is closed, structurally and by measurement.** Over a
   16× spread range and a 16× dispersion range, `nmae` moves by at most **2.8e−3 raw**, and at
   production's own dispersion by **≤3.6e−4**. The quasi law's per-gene totals are
   `floor(N·L·profile)` — independent of `spread` — and every emitted cell carries the same library
   total, so the arithmetic mean of per-cell CPM **is** the design pseudobulk.
3. **Candidate 2 (the call set) is closed by the metric definition.** `nmae`'s gate is the
   *reference's* significant set, and our predicted lfc is a function of the design profile alone.
   `rc_hybrid` and `rc_callset` deliver identical per-gene totals count-for-count. **`rc_callset`'s
   +0.0084 overall cost `nmae` exactly zero.**
4. **The biggest channel in the whole decomposition is a per-row additive offset worth +0.0542 raw
   at oracle — larger than the entire gap to Aginglab — and it is not estimable.** It is a property
   of the reference's *gate selection* (corr **+0.975** with the gate's median true lfc, **−0.029**
   with the whole-row shift), and the best honest fit on our own observables has **R² = −0.66**.
5. **The mandate's board evidence does not survive.** 15 of the 55 teams above overall 0.15 report
   `score_mse` **exactly 0.00000** — their minimum raw `mse` is **1.93516**, far past the 0.987863
   clamp — so 197's inversion returns a cosine that was never measured. `corr(implied cosine,
   scaled nmae)` is **−0.394 with them and +0.022 without them.** The "catastrophic cosine, great
   `nmae`" cluster is not a direction phenomenon; it is teams who **never α-shrank**, sitting at
   k ≈ 1 where the identity gives `nmae → 2(1−q) ≈ 0.89`.
6. **Aginglab's six-tuple is not a direction story, and INBOX-1 is right about why.** On
   `pds_cosine` — reported **directly**, no inversion, no α assumption — the direction hypothesis
   predicts their raw `pds` should be **+0.019 (h1-test) to +0.041 (PLOLO) HIGHER** than ours. It is
   **0.0130 lower**. Of five constructions tested against the joint shape, four fail by an order of
   magnitude; the only survivor is a **per-perturbation response-LEVEL** prediction — the channel
   §4.3 sizes at +0.0542 raw oracle and shows we cannot estimate. Dialled to their `nmae` it
   overshoots their `pds` cost by **2.71×**; dialled to their `pds` cost it delivers **59 %** of
   their `nmae`. **That is the mechanism, with its residual stated.**

**Nothing here is an arm. No flag changes. `nmae` is handed back to the direction programme, and
one new competitive instrument (§6.4) and one board-reading rule (§6.2) go to 221.**

---

## 0. Phase 0 — known-answer rows. Six gates, all passing.

| # | row | expected | measured | verdict |
|---|---|---:|---:|---|
| KA-1 | board affine map retrodicts **our own** scaled `nmae` | +0.11074 | **+0.10994** (err −8.1e−4) | **PASS** |
| KA-1b | the same map on Aginglab | +0.18398 | +0.18430 (err +3.3e−4) | PASS |
| KA-1c | fitted span vs 188's board-fitted 0.6049 | 0.6049 | **0.603369** (−0.25 %) | PASS |
| KA-1d | fitted baseline vs spec range 1.0009–1.0017 | — | **1.000222** | PASS |
| KA-2 | 187 §1's gene-wise identity `1 − 2W_c_hi − (2q_lo−1)K_lo + K_hi` | exact | **≤ 2.2e−16** on 7 tables | **PASS** |
| KA-3 | 187 KA-b: h1-test raw `nmae` at s = 1 on 170's table | 0.959881 | **0.959881** | **PASS** |
| KA-4 | my `nmae` vs `debench.Surface.nmae()` | < 1e−12 | **0.0** on 7 tables | PASS |
| KA-5 | 187 §6: h1-test `nmae(s)` turnover / gain / slope | 1.29 / 0.0027 / −0.0179 | **1.290 / 0.00269 / −0.0181** | **PASS** |
| KA-6 | board facts: 55 teams >0.15, corr, the three named rows | −0.394; CellWM +0.19847@−1.648 | **−0.3938**; **+0.19847@−1.6475** | **PASS** |

The map is `scaled = (1.000222 − raw) / 0.603369`, fitted on the 522 entries with raw ∈ [0.75,
1.10]; median |residual| 1.3e−3, p95 4.7e−3. It is **not** globally affine — `score_nmae` is a mean
of three per-context scalings while the reported raw is a mean of three raws — which is why 188
warns against this arithmetic. Quoted with its residual, it retrodicts our row.

**KA-2 needed a correction to 187's prose.** 187 §1 writes the identity's two accuracies as
"k-weighted"; they are not the same weighting. The exact form is

```
nmae = 1 − 2·W_c_hi − (2·qK_lo − 1)·K_lo + K_hi
```

with `W_c_hi` the **|t|-weighted** mass of correct-sign genes *above* the kink and `qK_lo` the
**k-weighted** sign accuracy *below* it. Written that way it is exact to 2.2e−16 on every table;
written with one weighting for both it reads **0.017 high**. 187's numbers are unaffected (it used
the identity only for slopes), but anyone reusing the formula needs the distinction.

---

## 1. The structural fact everything below rests on

The spec (§6, and `docs/reference/vcc2026-metrics-FULL.md:342`) forms **both** sides' lfc against
the **same measured control**, so in the numerator that control cancels exactly:

```
lfc_hat − lfc  =  log2(m_hat_pert + ε) − log2(m_real_pert + ε)
```

**`nmae` is an L1 error on our predicted arithmetic-mean-CPM profile of the perturbed state, over
the reference's own significant gate, divided by the reference's |lfc|.** Three consequences, and
they dispose of two of the mandate's three candidates before any sweep:

- **It reads no call set.** The gate is reference-only, and our lfc is a function of the design
  profile alone. Nothing about which genes *we* call enters.
- **It reads no emitter.** Per-gene totals are the design's count-for-count and every emitted cell
  carries the same library total, so `mean_i(x_ig / T_i) = (Σ_i x_ig)/(N·L)` exactly.
  `src/vcc2026/rc_emit.py:38` and `scripts/make_submission_pm_cd4.py:3328` both say so in
  production's own words: *"`nmae` is pinned"*.
- **It has no compositional degree of freedom.** The scorer normalises each cell to 1e6, so
  `Σ_g m̂_g = Σ_g m_g = 1e6` and the mass-weighted mean lfc is pinned to ≈0 on both sides. Measured:
  the whole-row median lfc has sd **0.0152** on the truth against the gate's **0.2423**.

---

## 2. Candidate 1 — emission spread and dispersion. **CLOSED, measured.**

`emit227.py`, 25 h1-test rows, the same design emitted twelve times. The deterministic design
profile is the zero-emission reference row (`nmae` 0.948292900).

| arm | φ | Δ`nmae` vs design (seed 90 000 / 424 242) | cell totals |
|---|---:|---:|---|
| production spread vector | 4 | −2.8e−5 / −1.6e−4 | [18 251, 25 397] |
| uniform spread 0.25 | 4 | −2.2e−5 / +8.3e−5 | [19 068, 22 128] |
| uniform spread 1.0 | 4 | −1.6e−4 / −3.6e−4 | [17 334, 25 335] |
| uniform spread 4.0 | 4 | −1.6e−4 / −3.6e−4 | [17 334, 25 335] |
| uniform spread 1.0 | 1 | +1.4e−4 / −8.3e−5 | [18 358, 22 538] |
| uniform spread 1.0 | 16 | +2.5e−4 / **+2.8e−3** | [15 449, 34 685] |

**Max |Δ| over a 16× spread range = 3.6e−4.** The single 2.8e−3 excursion is at φ = 16, four times
production's dispersion, where the per-cell library totals swing 2.2×; that residue is the ≤1-count
library imbalance the emitter's weld leaves, not a Jensen effect. Independent confirmation from
data already on disk: 169's five matched seeds of the same arm spread **2.55e−4**.

**Against a 0.04487 raw gap, the entire emission channel is 1.3 % of it at φ = 16 and 0.8 % at
production's φ. 034's spread lever is real for `fid`/`reach` and is exactly zero here.**

---

## 3. Candidate 2 — the DE pool and the call set. **CLOSED by definition.**

The mandate asks what `rc_callset` (n_pred 4 480 → 2 599, +0.0084 overall) did to `nmae`. **Nothing,
and it could not have.** `rc_callset` differs from `rc_hybrid` in exactly one thing — `call_set_first`
in pass 2's |z| re-assignment — and both deliver pass 1's per-gene totals count for count
(`make_submission_pm_cd4.py:3326-3334`). Same profile ⇒ same lfc ⇒ same `nmae`. The +0.0084 was
bought on `fid`/`jac`/`reach` and charged nothing here.

188's pool finding is about *reconstructing the reference's gate in a harness*, not about our call
set, and it is used as specified below (`|lfc|·√cpm` at 89 genes/row on PLOLO).

---

## 4. Candidate 3 — the lfc estimator. The full decomposition.

`phase1.py`, h1-test, `de169.q1.00_mf` (a23-family production geometry): 85 rows, gate median
**434** genes, base `nmae` **0.966399**, `pds` 0.80182, `cos_agg` 0.06195, `q_w` 0.5679, K 0.3322,
mean |t| on the gate **0.409** log2.

### 4.1 There is no pocket of error to find

| stratum (control CPM) | genes | share of Σ|t| | share of Σ|h−t| | K | q_w | h ≡ 0 |
|---|---:|---:|---:|---:|---:|---:|
| 5–10 | 4 318 | 0.1340 | 0.1375 | 0.0407 | 0.5291 | 0.0000 |
| 10–20 | 6 510 | 0.1458 | 0.1346 | 0.0445 | 0.5883 | 0.0000 |
| 20–50 | 13 382 | 0.2359 | 0.2268 | 0.0674 | 0.5806 | 0.0000 |
| 50–150 | 21 884 | 0.2579 | 0.2496 | 0.0841 | 0.5356 | 0.0000 |
| >150 | 18 813 | 0.2258 | 0.2172 | 0.0954 | 0.6027 | 0.0000 |

Numerator share tracks denominator share to within 1 % in every stratum, and the same is true
across |t| quintiles. **We recover almost nothing, uniformly.** Two hypotheses die here:

- **No unmodelled band.** Production ships `MIN_CPM = 5.0`, matching the scorer's gate:
  `h` is exactly zero on **0.00 %** of 64 930 gate pairs. There is no "scored but not modelled"
  5–10 CPM band, and no ε-catastrophe (|h| > 5 log2 on 0.015 % of pairs).
- **No low-expression noise floor to buy back.** The 5–10 CPM band carries 13.4 % of the
  denominator and 13.7 % of the numerator — proportionate, not excess.

### 4.2 The ceiling ladder

| channel | raw `nmae` | gain vs shipped |
|---|---:|---:|
| predict zero everywhere | 1.000000 | −0.033601 |
| **as shipped** | **0.966399** | — |
| best global scalar (s = 0.80) | 0.963053 | +0.003346 |
| oracle per-CPM-stratum scalar | 0.957531 | +0.008868 |
| oracle per-row scalar | 0.936699 | +0.029700 |
| LOPO per-gene additive offset | 0.935183 | +0.031216 |
| **oracle per-row additive offset** | **0.912199** | **+0.054200** |
| our signs, oracle magnitude | 0.864211 | +0.102188 |
| oracle signs, our magnitude | 0.727268 | +0.239131 |

**The per-row additive offset is the only channel bigger than the gap to Aginglab.** Everything
else that is realisable-shaped is an order of magnitude too small.

### 4.3 What the offset is, and why it is unreachable

`phase2.py`. The L1-optimal per-row constant `δ_p` (mean −0.0161, sd **0.2019**, q05/q95
−0.213/+0.223) is:

| quantity | corr with `δ_p` | sd |
|---|---:|---:|
| median **true** lfc on the gate | **+0.9747** | 0.2423 |
| median true lfc over all 10 779 scored genes | **−0.0290** | 0.0152 |
| median **our** lfc on the gate | — | 0.0264 |

**`δ_p` is entirely a property of the reference's gate selection.** It is not a compositional shift
— the CPM normalisation forbids one, and the whole-row median is 16× tighter than the gate's. It is
the sign skew of *which genes the reference happened to call significant for that perturbation*.

It is not estimable from anything we hold:

| arm | Δ`nmae` | Δ`pds` | Δ`cos_agg` |
|---|---:|---:|---:|
| oracle offset (gate only, unrealisable) | **+0.05420** | +0.0398 | +8.47 % |
| oracle offset applied to all genes (the only shippable shape) | +0.05420 | **−0.0594** | **−6.35 %** |
| our own median shift × 1 | −0.00030 | −0.0008 | +0.33 % |
| our own median shift × 8 | −0.00467 | −0.0125 | +2.37 % |
| LOPO fit on our observables | **−0.04620** | −0.0213 | −10.70 % |

The LOPO fit has **corr(pred, oracle) = −0.012 and R² = −0.66** — worse than predicting the mean.
Shuffle nulls: −0.0869 (h1-test), **−0.1328** (PLOLO). The channel is real, broad (leave-top-3-out
retains 70.4 %), and orthogonal to everything we can observe. **And even at oracle, the shippable
shape costs 6.4 % of cosine, because we do not know the gate.**

---

## 5. Phase 1's gate — is `nmae` reachable at fixed cosine? **No.**

`phase3.py`. The frontier is **anchored at our own prediction**: to raise `q_w` we flip a random (or
adversarially chosen) subset of the genes we currently get *wrong*, leaving every magnitude and
every other gene untouched, so the f = 0 row reproduces our exact point. Aginglab's credit is
**1.6787×** ours.

| route to 1.6787× credit | `q_w` needed | Δ`cos_agg` | Δ`pds` | reachable? |
|---|---:|---:|---:|---|
| fix wrong signs, uniformly at random | 0.6095 | **+7.3 %** | +0.0250 | yes |
| **fix wrong signs, smallest \|h\| first (oracle, cheapest in cosine)** | 0.7026 | **+5.0 %** | +0.0189 | yes |
| fix wrong signs, largest \|t\|/\|h\| first (oracle) | 0.7392 | +10.5 % | +0.0311 | yes |
| fix wrong signs, largest \|h\| first | 0.5992 | +23.0 % | +0.0426 | yes |
| **magnitude shape only, at our pinned budget** | 0.5679 | −5.6 % | −0.0001 | **NO — saturates at 1.487×** |
| global rescale | 0.5679 | **0.0 %** | **0.0** | **NO — turns over below 1.10×** |

Two hard results:

- **The magnitude channel cannot get there at all.** Blending a full oracle on |t| into our
  magnitudes while holding our signs *and our L2 budget* saturates at **1.487×** at λ = 0.5 and then
  turns back. The budget is pinned by `mse` (187/211), so the only way past 1.49× is amplitude —
  and h1-test's rescale row reads `nmae` **1.2725 at 3×**, worse than predicting nothing.
- **The minimum cosine cost of the sign channel, over an oracle choice of which signs to fix, is
  +5.0 %.** The gate is 434 of 10 779 genes, but those genes carry enough energy that fixing them
  all raises `cos_agg` by **74 %**. `nmae` and cosine are not separable on the improvement side.

**PLOLO, three cancer folds** (`plolo227.py`, 188's validated `|lfc|·√cpm` pool at 89 genes/row;
base `nmae` 0.973521, `pds` 0.76203, `cos_agg` 0.022466, `q_w` 0.5585):

| arm | Δ`nmae` | `q_w` | Δ`pds` | Δ`cos_agg` |
|---|---:|---:|---:|---:|
| best global scalar | +0.010980 | 0.5585 | **+0.000000** | **+0.000 %** |
| ORACLE per-row offset | +0.085766 | 0.6922 | −0.0572 | +334.7 % |
| our own median shift ×1 | +0.000260 | 0.5597 | −0.0010 | +3.38 % |
| SHUFFLE oracle offset | **−0.132795** | 0.5918 | −0.1887 | +25.6 % |
| fix lo-|h| wrong signs, f = 0.10 | +0.001754 | 0.6005 | +0.0016 | +0.478 % |
| fix lo-|h| wrong signs, f = 0.20 | +0.006497 | 0.6423 | +0.0060 | +1.775 % |
| fix lo-|h| wrong signs, f = 0.40 | +0.025585 | 0.7276 | +0.0233 | +7.058 % |

**Exchange rate: +0.003660 raw `nmae` per 1 % relative `cos_agg`.** Aginglab's +0.04487 raw needs
**+12.3 % relative cosine**. h1-test's independent cross-check at the same credit ratio: **+5.0 %**.
Both are 25–60× the mandate's ±0.2 % gate.

**`kolf2`, reported separately as required** (base `nmae` 0.996222, `q_w` **0.5067** — barely above
chance, `cos_agg` 0.009560): it inflates every oracle wildly (the per-row offset reads **+0.2925**
there against +0.0858 on the cancer folds) and its global-scalar gain is +0.0001. It is not usable
for absolute `nmae` and is excluded from every number quoted above.

**One honest positive, and it is not new.** A global scalar is *exactly* free on both `pds` and
`cos_agg` (both are scale-invariant; measured to 0.000 000) and buys +0.010980 raw `nmae` on the
cancer folds. It is 187's closed amplitude channel — the cost is entirely in `mse` — and it is not
even coherent across folds: the optimal scalar is **1.22 / 0.80 / 0.50 / 0.88** on
k562 / hct116 / hek293t / kolf2, and the +0.011 mean is carried almost entirely by hek293t
(+0.0274 at s = 0.50), the fold AGENTS.md flags for a 44.1 % target-independent blob. **Do not read
it as a lever.**

### The one arm with an observable key, priced in full

Phase 2 measured `corr(row nmae, ‖pred‖) = +0.42` — rows we predict loudly score worse — so a
per-row amplitude keyed on that is the only honest arm the decomposition supports.

| endpoint | value |
|---|---|
| Δ`nmae` | **+0.004128** |
| Δ`pds` | **0.000000** (per-row scaling is exactly `pds`-invariant — internal check passed) |
| Δ`cos_agg` | **−7.85 %** |
| shuffle null | −0.000610 |
| leave-top-3-out retention | **−18.4 %** (top 3 rows carry 118.4 % of the gain) |

**Fails the cosine gate by 39× and fails the concentration falsifier outright.** Reported as a
failure, not a trade.

---

## 6. Phase 2 — retrodicting Aginglab.com

### 6.1 Amplitude is ruled out, three independent ways

| test | result |
|---|---|
| 187's board law `−0.0416 per ln α` | their Δraw needs **α × 2.94** |
| h1-test `nmae(s)` at 3× | **1.2725** — far worse than predicting nothing |
| 187 §3's `mse` charge at α = 1.10 alone | −0.1216 scaled; their `mse` edge is **−0.00132** |

They are at essentially our amplitude, and their `mse` says so.

### 6.2 The board evidence the mandate rests on is a clamp artefact

`mse` is clamped to [0, 1] and its baseline is 0.987863. **15 of the 55 teams above overall 0.15
report `score_mse` = 0.00000 exactly, and the smallest raw `mse` among them is 1.93516.** 197's
inversion applied to a clamped value returns a cosine that was never measured — which is where
CellWM's −1.648, Liaoyi Xu's −1.330 and AgingLab_BI's −2.640 come from.

| statistic | value |
|---|---:|
| `corr(implied cosine, scaled nmae)`, all 55 teams | **−0.3938** |
| the same, on the **40 teams where the inversion is defined** | **+0.0224** |
| mean scaled `nmae`, clamped 15 | +0.13878 |
| mean scaled `nmae`, unclamped 40 | +0.10926 |

**The entire negative correlation is the clamp.** The cluster is real as a phenomenon and its
mechanism is not direction: raw `mse` ≥ 1.94 means they never α-shrank, and at k ≈ 1 the identity
gives `nmae → 2(1−q) ≈ 0.89–0.90`, which is exactly what they score. **`mse`'s floor at 0 makes
that trade one-way cheap for them and expensive for us** — we would forfeit +0.1796 scaled to take
it. That is the "way to score `nmae` that costs direction", and it is a trade we already priced and
declined.

**Aginglab is not in that cluster.** Their `score_mse` is +0.17825, unclamped, so their inversion is
valid to the precision of the instrument — which is the whole question.

### 6.3 The falsification, done on `pds_cosine` — INBOX-1's correction, adopted

**INBOX-1 is right and it changes which evidence is load-bearing.** 197's inversion is a *one-sided
lower bound*, tight only for an α-optimised team (0.911 tight on our own row), so "their implied
cosine is below ours" cannot carry the argument. **`pds_cosine` is reported directly on the board.**
So the falsification is redone on it, using the Δ`pds` my own frontier produces on every arm
(`phase4.py`):

| | predicted Δ raw `pds_cosine` at their `nmae` edge |
|---|---:|
| h1-test, cheapest sign-fix ordering (oracle) | **+0.01885** |
| PLOLO, 3 cancer folds, +0.9099 raw `pds` per raw `nmae` | **+0.04083** |
| **observed** (theirs 0.82276, ours 0.83572) | **−0.01296** |

**The direction hypothesis is wrong by +0.032 (h1-test) to +0.054 (PLOLO) raw `pds_cosine`, and
`pds` has no jackknife term and no clamp, so nothing about the emitter can rescue it.** Their
`nmae` edge is not better direction. **Standing rule, and this project has now walked into the trap
twice (216 §4.4's functional mismatch, and this one): when comparing a rival's direction to ours,
use the board's `pds_cosine`, not the `mse` inversion.**

### 6.4 The instrument caveat, which stands on its own

`recon227.py`. Separately from the above — and this is the reusable by-product INBOX-1 asks for —
197's inversion is calibrated on **our** arms, and the scored numerator credits
`ρ·min(Ĉ_p, C_p)` — the jackknife of our **own emitted noise**, 187's `claim`. That term is a
property of the *emitter*, not the direction. From 187's cached a23 constants
(`claim` A 0.233754 / B 0.237923 / C 0.266129, mean 0.245935; `ΣD` 0.229141):

```
a 1 % relative change in the jackknife claim  →  Δraw mse = 0.010733
                                              →  read by 197 as 5.20 % of our cosine
```

- Their raw `mse` is **+0.00125** vs ours — 197 reads that as **−0.61 %** cosine, which is **0.12 %**
  of a claim difference.

**The board cannot resolve a team's cosine below a few percent unless their emission law matches
ours.** Combined with 197's own tightness of 0.911 on our row, the inversion is a good instrument
for ranking teams whose emission resembles ours and a poor one for small differences. **This does
not rescue the direction hypothesis — §6.3 kills it on `pds`, which the inversion never touches —
but every competitive read that turns on a few percent of implied cosine should carry this band.**

### 6.5 What DOES fit their six-tuple

`phase4.py`, h1-test, five constructions priced on all five deterministic members in scaled board
units. Target = Aginglab minus us: `nmae` **+0.07323**, `pds` **−0.03025**, `fid` +0.01970,
`jac` −0.03295, `reach` −0.02213.

| construction | Δ`nmae` | Δ`pds` | Δ`fid` | Δ`jac` | Δ`q_w` | verdict |
|---|---:|---:|---:|---:|---:|---|
| **TARGET** | **+0.07323** | **−0.03025** | +0.01970 | −0.03295 | ? | |
| concentrate on our top-1000 by \|h\| | −0.01142 | −0.04709 | −0.23056 | −0.13819 | −0.436 | **wrong sign on `nmae`** |
| same profile, widen the call set to 9 000 | **+0.00000** | +0.00000 | −0.00599 | −0.00743 | 0.000 | **exactly zero on `nmae`** |
| CPM tilt β = −0.15, row L2 fixed | +0.00509 | −0.05515 | 0 | 0 | 0.000 | signs right, **26× too little** |
| ORACLE: fix 25 % of gate signs, break 10 % off-gate | +0.10100 | **−0.96039** | 0 | 0 | +0.112 | **32× too much `pds`** |
| **per-row LEVEL offset × γ** (see below) | **dialable** | | 0 | 0 | 0.000 | **the only survivor** |

The level channel — §4.3's per-row offset, the one channel whose sign is right on every member —
dialled by γ and applied to the whole profile:

| γ | Δ`nmae` | Δ`pds` | Δ`cos_agg` |
|---:|---:|---:|---:|
| 0.20 | +0.02949 | −0.02231 | +0.37 % |
| **0.30** | **+0.04298** | **−0.03408** | +0.20 % |
| 0.40 | +0.05523 | −0.05019 | −0.20 % |
| **0.593 (hits their `nmae`)** | **+0.07323** | **−0.08195** | −1.5 % |
| 1.00 | +0.08983 | −0.13136 | −6.35 % |

**At their `pds` cost (γ ≈ 0.30) it delivers 59 % of their `nmae` edge; at their `nmae` it overshoots
their `pds` cost by 2.71×.** Nothing else tested is within an order of magnitude of the joint shape.

**Answer to the mandate's Phase 2.** Aginglab is not doing anything to the emitter, the call set, or
the amplitude. Their six-tuple is the signature of a model that predicts the **per-perturbation
response LEVEL** — how far, and in which direction, that perturbation moves its responders — which
is a *calibration* object, not a direction object. Our pipeline has no channel for it: we ship a
whitened, shrunk, mass-consistent **shape**, and §4.3 shows the level we do emit is uncorrelated
(R² = −0.66) with the truth's. The residual factor of 2.7 on `pds` is what separates a crude
uniform offset from a targeted one, and is the honest size of what is *not* retrodicted.

**Their tightness, as INBOX-1 asks:** their `mse` is unclamped and within 0.00125 raw of ours at a
`pds_cosine` 0.0130 *below* ours. A team at a lower cosine printing our `mse` is a team collecting
at least as much jackknife credit as we do — so **Aginglab is α-optimised, and 197's floor is
approximately tight on them.** It is *not* tight on the 15 clamped teams, where it is not defined at
all.

### 6.6 INBOX-2's question, answered: `q` and `pds` do **not** come apart

INBOX-2 asks whether `q` (|t|-weighted sign accuracy on the DE pool) can be raised at fixed `k`,
fixed cosine and fixed `pds`. **Measured: no, by a wide margin, in both directions.**

- Raising `q_w` by fixing gate signs and touching nothing else raises `pds` by **+0.095 scaled**
  (Δ`q_w` +0.111) — the two move together, hard.
- Forcing `pds` down while `q_w` rises requires breaking off-gate signs, and breaking just **10 %**
  of them costs **−0.960 scaled `pds`** — **32× the −0.030 the target allows.** The gate is 434 of
  10 779 genes and simply does not have enough leverage to be decoupled from the rest.

So the INBOX-2 hypothesis's "if it is false" branch is the one that fires — but its stated
alternative (the call set) is also false: **`nmae` reads no call set at all** (§3, and the
zero-to-five-decimals measurement in the table above). The real alternative is §6.5's level channel.
**224's gain vector was not on disk at the time of this run**, so it could not be tested directly as
a ready-made `q`-raiser; the general result above applies to any `q`-raiser, including it: **if 228
sees `q` rise, it will see `pds` rise with it, and the two arms are then redundant exactly as
INBOX-2 warns.** The one thing 228 should NOT expect is `q` up with `pds` flat.

---

## 7. Phase 3 — the mandatory falsifiers

1. **Cosine-pinned control.** Every reported `nmae` number in §4–§6 carries `pds` and 197's
   energy-weighted `cos_agg` from the same run. No exceptions. The two arms that gained `nmae` and
   are shippable in shape both lost cosine (−7.85 %, −6.35 %) and are reported as failures.
2. **Shuffle nulls.** Oracle offset −0.0869 (h1-test) / **−0.1328** (PLOLO); LOPO offset −0.0565;
   per-row amplitude −0.00061. Every fitted arm's null is negative and large, which is the correct
   signature for a row-specific quantity — and is why getting the row wrong is so expensive.
3. **Leave-top-3-out retention on `nmae`** (222 §9's protocol, per-target deletion):

   | arm | total | k = 1 | k = 3 | k = 5 |
   |---|---:|---:|---:|---:|
   | per-row amplitude (the honest arm) | +0.004128 | 44.5 % | **−18.4 %** | −49.0 % |
   | global scalar 0.80 | +0.003346 | 60.5 % | 19.3 % | −8.3 % |
   | oracle per-row offset | +0.054200 | 89.3 % | **70.4 %** | 56.1 % |

   **The channel that is broad is the one that cannot be estimated; the channels that can be
   estimated are three rows deep.** 222 §9 was right that the channels differ and must be measured:
   `nmae`'s realisable gains are *more* concentrated than `pds`'s, not less.
4. **October (panel replaced 2026-10-22, 300 disjoint targets).** The emission-side result is
   **panel-independent by construction, and it is a negative**: `nmae` reads a functional of the
   design profile that the emitter preserves exactly — `floor(N·L·profile)` per gene and an equal
   per-cell library — so no panel can make `spread` or `overdispersion` matter. The same holds for
   the call set. The offset channel is structurally panel-independent too (it is a property of how
   the reference selects its gate) and unpredictable on any panel. **The only `nmae` channel that
   *is* panel-dependent is the global scalar**, because its optimum sits against the per-context
   `mse` cliff — and 187 §7 already owns that with its `knee × (1 − 1.25 %)` D/E/F rule. **Nothing
   in this experiment adds October risk, and nothing in it reduces the value of the source ladder.**

---

## 8. What this closes, and what it hands on

- **`nmae` as a lever we can pull: closed.** It is a |t|-weighted L1 on the same predicted profile
  `pds` and `mse` read, restricted to 4 % of the genes. Its channels are sign accuracy on that gate
  (which drags `pds` and cosine with it in lock step, ±0.095 scaled per +0.11 `q_w`), amplitude
  (pinned by `mse`), and the response level (§4.3, unestimable from anything we hold). **There is no
  fourth thing, and none of the three is separable from cosine at the mandate's gate.**
- **But 187's dichotomy was still incomplete, and this is the finding to carry.** 187 called `nmae`
  "a direction problem, not an amplitude problem." **It is a *level* problem** — a per-perturbation
  scalar worth **+0.0542 raw at oracle**, 1.2× the whole gap to Aginglab, which is neither direction
  nor amplitude and which our pipeline has no channel for. We cannot estimate it (R² = −0.66) and
  §6.5 says Aginglab probably can. **That is the one genuinely open door this experiment found, and
  it belongs to whoever owns response-magnitude modelling — not to the emitter and not to `α`.**
- **The emitter is closed for `nmae`, permanently and by construction.** Do not sweep `spread`,
  `overdispersion`, `resample_x` or any call-set width against this member again. They are exactly
  zero, and the source says so.
- **Two board-reading rules, and they are the most reusable things here. Both go to 221.**
  1. **Compare a rival's direction with `pds_cosine`, never with the `mse` inversion.** 197's
     inversion is a one-sided bound whose slack is an *emitter* property: **a 1 % difference in the
     jackknife claim reads as 5.20 % of transfer cosine.** `pds_cosine` is reported directly, has no
     jackknife term and no clamp. This trap has now cost the project twice.
  2. **Filter `score_mse` > 0 before inverting.** 15 of the 55 top teams are clamped; any statistic
     over them mixes a measurement with a floor. The −0.394 that motivated this experiment is that
     mistake, and it flips to **+0.022** when corrected.
- **To 228 (INBOX-2's redundancy question), the cheap decisive number.** If the θ gain raises `q`,
  it will raise `pds` with it — measured slope **+0.86 scaled `pds` per unit `q_w`** on the gate, in
  both directions and over an oracle choice of which signs to fix. **`q` up with `pds` flat is not a
  thing that happens on this data**, so a `q` rise is not independent evidence of a second lever;
  price the two marginally, per 208's rider rule.
- **`nmae` will follow cosine; it will not lead it.** Nothing here changes the value of the source
  ladder, and nothing here competes with 224/225/226 for a slot.

## 9. Files

| file | what it does |
|---|---|
| `phase0.py` | the six known-answer gates; board affine map; the corrected identity → `work/phase0.json` |
| `phase1.py` | strata, zero audit, the ceiling ladder, cosine on every arm → `work/phase1.json` |
| `phase2.py` | what the per-row offset is; realisable estimators; shuffle nulls → `work/phase2.json` |
| `phase3.py` | the anchored (q_w, K) frontier, the separability ceilings, retention, the clamp → `work/phase3.json` |
| `phase4.py` | INBOX-1/2: the six-tuple fit, the `pds` falsification, the level channel, `q` vs `pds` → `work/phase4.json` |
| `emit227.py` | the emission falsifier: 12 emissions, 4 spreads × 3 dispersions × 2 seeds → `work/emit227.{json,log}` |
| `plolo227.py` | PLOLO replication, 3 cancer folds + `kolf2` separately, 188's validated pool → `work/plolo227.json` |
| `recon227.py` | 197's inversion sensitivity to the jackknife claim; the Aginglab reconciliation → `work/recon227.json` |

Inputs read, never written: `experiments/169-sharpening-six-member/de169.q1.00_mf*.npz`,
`de169.KA_a15_raw.npz`, `work/h1test-effect.npz`; `experiments/170-gamma-mse/work/de_g0.65.npz`;
`experiments/188-magnitude-calibration/work/{chalframes,pair-*}.npz`;
`experiments/187-joint-amplitude/work/mse187.json`;
`experiments/162-retrodicting-validator/work/fold-*.npz`;
`data/leaderboard/2026-09-04T0433Z.json`; `scripts/make_submission_pm_cd4.py` (imported, not
modified); `src/vcc2026/{donor_transfer,rc_emit}.py`.

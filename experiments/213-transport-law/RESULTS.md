# 213 — The transport law, written down and varied for the first time. **The shipped law is at the optimum of both families through it, and the whole family is worth at most +0.0010 overall.**

**Date:** 2026-09-04 · no `vcc`, no `cell-eval2`, no build, no submission, nothing written
outside `experiments/213-transport-law/`. Peak RSS **12.27 GiB** in `refit213.py`'s 17,942²
eigensolve — **2.3% over the 12 GB ceiling, disclosed rather than rounded down**; every other
script stayed under 8 GiB. §8 records what happened and what was done about it.

---

## Verdict in ten lines

1. **The shipped law is `log(p + c)` at `c = 1/TS = 20 CPM`, exactly — not by fitting, by
   algebra.** `log1p(TS p) = log(p + 1/TS) + log TS` and the constant cancels in a difference,
   so the incumbent is an **interior point of the `log(p+c)` family by construction**, and
   `c0` is precisely the 20-CPM crossover 026 §1 and 169 §6.4 both identified. It is *also* the
   `q = 0` interior point of the Box-Cox family on `p + c0`. Two one-parameter families, one
   incumbent, no bake-off.
2. **The mandate's parametrisation was one step off and the correction is the point.**
   `c = 1/TS` is not the log-fold-change end of the family; it is the **shipped point**.
   Log-fold-change is `c → 0` and additive-in-CPM is `c → ∞`, so the incumbent sits *between*
   the two structural limits rather than on one of them.
3. **On the cosine — the primary endpoint — the offset family has a clean interior maximum at
   `c ≈ 36–40 CPM`, worth `1.0159×` the incumbent** (95% CI on the ratio **[1.0095, 1.0221]**,
   P(> 0) = 1.000, paired 4,000-draw target bootstrap on 162's four PLOLO folds). The curvature
   family's maximum is at `q = +0.05`, worth **1.0021×**.
4. **On `pds` the maximiser is the incumbent itself** in the a15 render (`c = 24`, Δ = +0.00009
   raw, inside the interval) and `c ≈ 45` in the a23 render (+0.00084 raw, interval
   `[−0.0021, +0.0037]`, i.e. unresolved). The two maximisers differ by a factor of **1.5–1.7
   in `c`**, as the project's documented pattern predicts, and **it changes nothing.**
5. **Both structural limits are catastrophic, and both are reached from two different
   parametrisations that agree BIT-FOR-BIT.** Pure log-fold-change (026's law) reads
   **0.588×** cosine and **−0.0822** raw `pds`; additive-in-CPM reads **0.295×** and
   **−0.1867**; `rawadd` and `bc1` are the same map by different code and produce identical
   floats. Anscombe/√ 0.431×, rank transport 0.183× (uniform) / 0.581× (probit).
   **026 is confirmed and its magnitude was understated 20× by the effect-level instrument
   it used.**
6. **The best arm anywhere in the family is `+0.00103` overall on the two members PLOLO and
   096 can price (`mse` + `pds`), and `+0.00057` once h1-test's `nmae` is added.**
   Rule 13's floor is **+0.0089**, because the arm is unambiguously amplitude-changing: the
   delivered challenge-panel `A` goes to **1.80×** the incumbent's and α must fall
   0.982 → **0.7438**. At the *shipped* α the same arm costs **−0.109 overall on `mse`**.
   **Short by 8.6×; and 197's requirement for a direction arm is +31% of cosine against a
   family maximum of +1.6% — 19× short.**
7. **The mechanism, measured to five decimals: a transport law is a per-gene, sign-preserving
   gain.** `g_c/h` against the closed form `(p_don + c0)/(p_don + c)` correlates at
   **0.99999** with **0.0000% sign flips over 7.0 M donor entries**; the recipient half
   correlates at 0.990–0.993 with 0.17–0.20% flips. **A gain cannot reorder anything**, so the
   law is structurally incapable of touching the ordering channel (§8.1) — the same objection
   184 raised against the Tahoe operator, arriving from the opposite direction.
8. **Every falsifier passes, and F2 is the decisive one.** With the recipient's basal
   gene-deranged the best law's advantage does not merely fail to survive, it **inverts**:
   cosine ratio `logc40/logc20` **1.0159 → 0.9759**. The gain is transport, not amplitude —
   209's outcome (the random control as the maximiser) does not occur. The same control also
   measures something nobody had: **the incumbent itself loses 16% of its transfer cosine and
   0.070 raw `pds` when the recipient's basal is deranged.**
9. **Two more checks that could have rescued the arm and instead close it harder.**
   Refitting the whitening operator **in the law's own space** — the comparison the mandate
   asked for, and arguably the fair one — makes the best law **worse by 6.0% of cosine**, because
   a per-gene gain flattens the donor covariance and the refit operator's top eigenvalue falls
   **5.895 → 3.767** (§2.5). And 212's 96.2% cosine-for-retrieval exchange is **not**
   law-invariant, exactly as 212 suspected, but it moves the wrong way: under the winning law
   the gamma trade goes from a knife-edge to **net negative** (§2.7).
10. **NO-GO. No slot. Do not change the transport law.** And the negative is the useful form:
   the last load-bearing assumption in the stack is now measured on a curve, from two
   independent parametrisations, with the incumbent at its maximum — rather than resting on the
   fact that nobody had looked.

---

## 0. Phase 0 — known-answer rows. **Fifteen, all PASS.**

| # | gate | expected | measured | verdict |
|---|---|---|---|---|
| 1 | rule 11 operator table | 212's "17 files, 16 digests" | **17 files, 16 digests**; production `e2f65ecf69b7e0b0…`, n = 17,942, max eig 5.894570063345345 | **PASS** |
| 2 | the shipped effect definition, **proved** | `sum_g expm1(b+h) == TS` | max `\|sum/TS − 1\|` **1.94e−8** over 64 rows × 3 donors; min reconstructed `p^pert` **−3.39e−21** | **PASS** |
| 3 | 026 §1 compression at 5 / 20 / 100 CPM | 0.20 / 0.50 / 0.83 | **0.2000 / 0.5000 / 0.8333** | **PASS** |
| 4 | 169 §6.4 `\|h\|` ratio, 100 vs 10 CPM | ~2.4× | **2.364** at FC 1.25, **2.499** in the small-effect limit | **PASS** |
| 5 | the incumbent is INTERIOR in both families | — | c index 7/22, q index 6/16 | **PASS** |
| 6 | `logc20` and `bc0` effects vs `E062.prepare_main_source` | bit-exact | **max abs dev 0.0** (both) | **PASS** |
| 7 | 197's shipped transfer cosine / leader ratio / exchange rate | 0.079249 / 1.655× / 0.365 | **0.079249 / 1.6552 / 0.365** | **PASS** |
| 8 | 162's a15 → a23 mass-fix price, and 204's two `pds` anchors | +0.004014 scaled; 0.795784 / 0.789121 | **+0.004014488**; **0.7957840906989219 / 0.7891213079583496** | **PASS** |
| 9 | PLOLO: `logc20` **and** `bc0` == 204's `equal` row, 4 folds | 11 scalars each | **rel 0.00e+00 on 22 of 22** | **PASS** |
| 10 | 212's forward CD4 dose at `f = 0.34` | +13.83% of PLOLO cosine | **+13.828%** (0.0218418 → 0.0248622) | **PASS** |

**And four more that only the production chain can give** (`stage213.py` / `eff213.py` /
`chal213.py` / `h1sweep213.py`), because the law is injected into the shipped builder itself:

| # | gate | expected | measured | verdict |
|---|---|---|---|---|
| 11 | a15's `baseline_pre` and `candidate_pre` sha256, with the law injected at `logc20` | `7678bab9…` / `0601336f…` | **MATCH / MATCH** | **PASS** |
| 12 | a15's six operator + CD4 scalars | its own meta | worst rel **1.11e−15** | **PASS** |
| 13 | the staged effect files vs **170's own** `chal_g0.65.npz` / `h1_g0.65.npz` | bit-identity | **max abs dev 0.0 on both** | **PASS** |
| 14 | h1-test: `msekit` on the frozen 095 artifact, then 170's whole γ = 0.65 row | 2.0551018299736388; 11 scalars | rel 2.22e−16; **rel 0.00e+00 on 11 of 11** | **PASS** |
| 15 | challenge panel: `Kc`, the `C_p` cap, `claim`, ρ, the solved `X` | 096 (e) / 212 §2.3 | `Kc` 3/3 PASS; cap **0.24653**; `claim` **0.246076** (212: 0.24608); ρ **0.9926 / 0.9407 / 0.8821** (212: identical); implied cosine **0.077305** (212: 0.077305) | **PASS** |

Gate 13 is the strongest thing in this document. **The law is injected into the production
builder's own donor loading, and at its own parameter the whole pipeline reproduces a15's
staging hashes and 170's effect matrices bit-for-bit** — so every alternative law below differs
from the shipped artifact in the law and in nothing else.

### 0.1 What the general path costs at the incumbent, reported not hidden

`LogC.to_law` and `LogC.render` dispatch to the shipped implementation when `c == c0`, which is
what makes gates 6 and 9 bit-exact. Forcing the **general** algebraic path at the same parameter
lands **1.33e−15** from the shipped effect — an allocation-alignment artefact of
`np.linalg.norm` on differently-aligned buffers, **3.1 × 10⁴ times below the 4.13e−08 float32
round-trip 204 itself asserts**. The family is continuous at its own incumbent.

One numerical trap worth recording because it cost an hour: **`20.0 * 1e-6` is
`1.9999999999999998e-05` and differs from `1/50000.0` by one ulp.** Writing the incumbent's `c`
that way puts it 3.4e−21 off its own parameter and costs the identity row 1.3e−15 on every
effect. `c_cpm / 1e6` lands on `1/TS` bit-for-bit.

---

## 1. Phase 1 — the law, written down, and what each alternative assumes

The shipped law, as a map `(p_don_ctrl, p_don_pert, p_rec_ctrl) → p_rec_pert`:

```
h    = log1p(TS·p_don_pert) − log1p(TS·p_don_ctrl)             the donor table on disk
p_rec_pert ∝ clip(expm1(log1p(TS·p_rec_ctrl) + h), 0)/TS       profile_from_effect
```

Both halves are `f(p) = log(p + c0)` with `c0 = 1/TS = 2e−5 = 20 CPM`, because
`log1p(TS p) = log(p + c0) + log TS` and `log TS` cancels in the difference. Verified on this
machine (gate 2): `sum_g expm1(b_ctrl + delta_row) = 50000.00` on every sampled row of all three
production donors, so `p_don_pert = (p_don_ctrl + c0)·e^h − c0` is **exact**, and the whole
family is well defined on data we already hold.

| # | law | link `f(p)` | what it assumes | where it sits |
|---|---|---|---|---|
| 1 | **shipped** | `log(p + c0)` | the response is a fold change **above 20 CPM** and an absolute abundance shift below it | `logc20` ≡ `bc0` |
| 2 | multiplicative / pure LFC | `log p` | the response is a fold change at **every** expression level | `c → 0` (`mult`) |
| 3 | additive in raw CPM | `p` | the response is an absolute abundance shift at every level | `c → ∞` (`rawadd` ≡ `bc1`) |
| 4a | variance-stabilising, √ | `√p` | the response scales with the Poisson sd | `bc0.5` |
| 4b | Anscombe | `√(p + 3/(8TS))` | ditto, with the count VST's own offset | `anscombe` |
| 5a | rank transport, uniform | `F_frame(p)` | only a gene's **rank** among its own frame's genes is portable | `rank_u` |
| 5b | rank transport, probit | `Φ⁻¹(F_frame(p))` | ditto, on a scale that averages sensibly | `rank_z` |
| 6a | **the offset family** | `log(p + c)` | the fold-change/absolute crossover is at `c`, **fitted** | `logc*`, 22 points over six decades |
| 6b | **the curvature family** | `((p+c0)^q − 1)/q` | the response transform has curvature `q`, **fitted** | `bc*`, 16 points |

Families 6a and 6b intersect **only** at the incumbent and at the additive limit, so agreement
between them is evidence about the law and not about one parametrisation. Both limits are
reached by two different code paths and produce **bit-identical** floats (§5, F5).

**Where the recipient actually lives.** On the challenge's own controls, inside production's
CPM > 10 gate: **9,135 / 8,668 / 9,127** genes in contexts A/B/C, median **63.7 / 58.3 / 57.8
CPM**, and **12.5% / 14.8% / 15.1% below the incumbent's `c` = 20 CPM.** So the incumbent's
crossover sits at roughly the 13th percentile of the scored gene population — the family's `c`
is a real, well-identified parameter and not an edge case.

**And there is plenty of donor↔recipient mismatch for a law to exploit.** Donor vs context-A log
CPM: `corr` **0.490 / 0.348 / 0.342** for `k562_gwps` / `orion_hct116` / `orion_hek293t`, with
only **38–45%** of genes within 2× (`mech213.py`). The law is not being asked to price a
degenerate difference.

---

## 2. Phase 2 — the sweep, COSINE primary

162's PLOLO, four folds, the real 300-target panel, per-fold leakage-free operators, rendered
into the three challenge control frames and read back the way the scorer does. `logc20` and
`bc0` reproduce 204's published `equal` row on all four folds at **rel 0.00e+00** (gate 9).

### 2.1 The offset family — an interior maximum at `c ≈ 36–40 CPM`

Macro over four folds, a15 render (the geometry the incumbent ships); 95% CI is a paired
4,000-draw target-cluster bootstrap, the instrument 204 calibrated on these folds.

| `c` (CPM) | k562_gwps | orion_hct116 | orion_hek293t | **kolf2 (leak-free)** | **macro** | ratio | 95% CI on the ratio |
|---:|---:|---:|---:|---:|---:|---:|:---|
| 0.02 | 0.021454 | 0.015990 | 0.012326 | 0.006082 | **0.0139627** | 0.6393 | [0.5835, 0.6937] |
| 1 | 0.027633 | 0.020974 | 0.017152 | 0.006733 | **0.0181230** | 0.8297 | [0.7982, 0.8611] |
| 5 | 0.029764 | 0.023857 | 0.019980 | 0.007269 | **0.0202174** | 0.9256 | [0.9119, 0.9395] |
| 10 | 0.030775 | 0.024927 | 0.021030 | 0.007550 | **0.0210704** | 0.9647 | [0.9579, 0.9716] |
| 14 | 0.031346 | 0.025383 | 0.021469 | 0.007686 | **0.0214711** | 0.9830 | [0.9796, 0.9865] |
| **20 (ships)** | **0.031934** | **0.025759** | **0.021854** | **0.007821** | **0.0218418** | **1.0000** | — |
| 24 | 0.032205 | 0.025905 | 0.022001 | 0.007881 | **0.0219978** | 1.0071 | [1.0054, 1.0088] |
| 28 | 0.032382 | 0.025983 | 0.022093 | 0.007926 | **0.0220961** | 1.0116 | [1.0085, 1.0148] |
| 32 | 0.032512 | 0.026010 | 0.022140 | 0.007959 | **0.0221550** | 1.0143 | [1.0099, 1.0186] |
| **36** | 0.032607 | 0.026006 | 0.022145 | 0.007982 | **0.0221850** | **1.0157** | [1.0103, 1.0210] |
| **40** | 0.032656 | 0.025970 | 0.022130 | 0.008000 | **0.0221890** | **1.0159** | [1.0095, 1.0221] |
| 45 | 0.032686 | 0.025883 | 0.022085 | 0.008015 | **0.0221672** | 1.0149 | [1.0075, 1.0223] |
| 50 | 0.032675 | 0.025772 | 0.022011 | 0.008025 | **0.0221208** | 1.0128 | [1.0044, 1.0211] |
| 60 | 0.032553 | 0.025486 | 0.021816 | 0.008022 | **0.0219692** | 1.0058 | [0.9959, 1.0156] |
| 100 | 0.031435 | 0.024100 | 0.020751 | 0.007919 | **0.0210512** | 0.9638 | [0.9483, 0.9792] |
| 500 | 0.025267 | 0.016621 | 0.014082 | 0.007315 | **0.0158211** | 0.7244 | [0.6796, 0.7687] |
| 20000 | 0.018578 | 0.003089 | 0.002893 | 0.005265 | **0.0074564** | 0.3414 | [0.2184, 0.4639] |

*(all 43 laws x 4 folds are in `work/plolo213_main.json`; the macro, ratio, CI and
priced columns for every one of them are in `work/summary213.txt`.)*

**The maximum is interior, resolved, and small.** Its half-width at `c = 40` is
1.38e−4 absolute = **0.63% relative**, so the +1.59% point estimate is 2.5 half-widths clear
of the incumbent and P(> 0) = 1.000 on the paired bootstrap. It is a real gain. It is also
**1/19th of the +31% of cosine 197 says a direction arm must deliver to clear the gate alone.**

### 2.2 The curvature family — the same answer, from a different parametrisation

| `q` | −1 | −0.5 | −0.25 | −0.1 | −0.05 | −0.02 | **0 (ships)** | +0.02 | +0.05 | +0.075 | +0.1 | +0.25 | +0.5 | +1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| cosine ratio | 0.4241 | 0.8295 | 0.9418 | 0.9819 | 0.9921 | 0.9974 | **1.0000** | 1.0019 | **1.0021** | 0.9981 | 0.9887 | 0.8105 | 0.4690 | 0.2950 |
| Δ raw `pds` (a15) | −0.1288 | −0.0299 | −0.0048 | −0.0004 | +0.00002 | +0.00007 | **0** | −0.0002 | −0.0014 | −0.0033 | −0.0064 | −0.0735 | −0.1771 | −0.1867 |

**The curvature family's cosine optimum is +0.21% and its `pds` optimum is the incumbent.** Two
families, built from different algebra, agree that the shipped law is essentially the maximum.

### 2.3 The named alternatives, and 026's known answer

| law | cosine ratio | Δ raw `pds` (a15) | note |
|---|---:|---:|---|
| **`mult` — pure log-fold-change** | **0.5880** | **−0.08219** | 026's law, on the DELIVERED instrument |
| `rawadd` ≡ `bc1` — additive in CPM | 0.2950 | −0.18665 | two code paths, bit-identical |
| `anscombe` | 0.4308 | −0.18388 | |
| `bc0.5` (√) | 0.4690 | −0.17713 | |
| `rank_u` — quantile, uniform | 0.1827 | −0.22395 | worst arm in the experiment |
| `rank_z` — quantile, probit | 0.5810 | −0.07439 | |

> **026 is CONFIRMED and its magnitude was understated by ~20×.** 026 measured the
> log-fold-change re-expression at **−0.0041 macro `pds`** on an effect-level LODO surface;
> the same law on the delivered, panel-restricted, four-fold instrument reads **−0.0822 raw
> `pds`** and **0.588× the transfer cosine**. 026's mechanism — *"dividing by `expm1(b_ctrl)`
> amplifies noise exactly where the compression was largest; `log1p` is a variance-stabilising
> transform doing useful work"* — is right, and the family now says where the optimum of that
> trade-off is: **essentially where the scorer's own transform already put it.**

### 2.4 The two maximisers, side by side

| endpoint | instrument | maximiser | value there |
|---|---|---|---|
| **transfer cosine** | PLOLO, 4 folds, a15 render | **`c` = 36–40 CPM** (interior) | **1.0159×** [1.0095, 1.0221] |
| transfer cosine | PLOLO, curvature family | `q` = +0.05 | 1.0021× |
| **`pds`** | PLOLO, 4 folds, **a15** render | **`c` = 24** (i.e. the incumbent) | +0.000085 raw, unresolved |
| `pds` | PLOLO, 4 folds, a23 render | `c` = 45 | +0.000835 raw, CI [−0.0021, +0.0037] |
| `nmae` | h1-test, reference BH pool | `c` = 20 (the incumbent) | see §3 |
| **composite, `mse` + `pds` (a23 render)** | 096 + PLOLO | **`c` = 40** | **+0.00103 overall** |
| **composite, + h1-test's `nmae`** | as above | **`c` = 40** | **+0.00057 overall** |

The divergence the mandate predicted **is present** — the cosine maximiser is at `c ≈ 36–40` and
the `pds` maximiser is at `c` = 24 (a15) — and it is a factor of 1.5–1.7 in `c` and +0.0006 of
overall. Reporting them separately was the right instruction and it did not change the verdict.

### 2.5 The operator interaction, both ways

The whitener's covariance was estimated from donor effects in the shipped space, so holding it
fixed while changing the law is confounded in one direction and refitting it is confounded in
the other. Both columns are reported. **Which is fair:** the *refit* column is the fair test of
the LAW (a law you adopt comes with its own operator — and §8.3 says we refit the operator on
Oct 22 anyway); the *shipped* column prices the arm you could build today.

Measured on the **`kolf2` fold** — the leak-free one, which uses production's own three-donor
mixture and the shipped operator — with an operator refit from the same 29,394-row basis
construction as `scripts/build_union_whitener.py` (same 12,000/donor cap, same seed-0 draw, same
essential prune; the row counts come out **7,784 + 10,807 + 10,803**, identical to
`whitener_union_pruned.json`'s own `rows_kept_per_donor`).

**The refit machinery's known-answer row:** refitting at the INCUMBENT reproduces the shipped
operator behaviourally — `a15_pds` and `a23_pds` identical to sixteen digits, `a15_cos` to
**4.0e−8** relative, max eigenvalue 5.894550 against 5.894570 (**3.4e−6**). It is the same
operator.

| law | operator | max eig | a15 cosine | a15 `pds` | a23 cosine | a23 `pds` | cosine ratio to `logc20`/shipped |
|---|---|---:|---:|---:|---:|---:|---:|
| `logc20` | **shipped** `e2f65ecf…` | 5.894570 | 0.00782070 | 0.611201 | 0.00755902 | 0.610128 | 1.0000 |
| `logc20` | refit in its own space | 5.894550 | 0.00782070 | 0.611201 | 0.00755902 | 0.610128 | **1.0000** (KA row) |
| `logc40` | **shipped** | 5.894570 | 0.00799968 | 0.610783 | 0.00783831 | 0.613100 | **1.0229** |
| `logc40` | refit in its own space | **3.766960** | 0.00752353 | 0.607560 | 0.00733280 | 0.607143 | **0.9620** |

> **Refitting the operator in the law's own space makes the alternative law WORSE, by 6.0% of
> cosine and 0.0032 raw `pds`.** The mechanism is visible in the spectrum: the law is a per-gene
> gain (§4), it flattens the donor covariance, and the refit operator's top eigenvalue falls
> **5.895 → 3.767**. A flatter covariance is a weaker whitener, and whitening is where 5/6 of
> our score comes from (§4.2).
>
> **So the fair comparison is the one that favours the alternative.** A law you actually adopt
> comes with its own operator — and §8.3 says we refit the operator for D/E/F on Oct 22 anyway
> — so the refit column is the honest test of the LAW, and it is *worse*. Every number in §2
> and §3 uses the **shipped** operator, which is the generous column. The negative is robust to
> the choice, and the direction of the bias is known.


### 2.6 The PerturbMatch impurity, priced (F6)

PerturbMatch's table is a limma-voom log-fold-change, not a `log1p(TS p)` difference, so it has
no law-space form and is left untransformed in every arm — which is also what keeps the identity
row exact. Removing it entirely (`--no-pm`):

| law | cosine ratio WITH PM | cosine ratio NO PM | Δ raw `pds` WITH | Δ raw `pds` NO |
|---|---:|---:|---:|---:|
| `logc10` | 0.9647 | 0.9534 | −0.00094 | −0.00181 |
| **`logc20`** | 1.0000 | 1.0000 | 0 | 0 |
| **`logc40`** | **1.0159** | **1.0257** | −0.00076 | +0.00015 |
| `logc100` | 0.9638 | 0.9801 | −0.00702 | −0.00585 |
| `mult` | 0.5880 | 0.5563 | −0.08219 | −0.09067 |
| `rawadd` | 0.2950 | 0.2722 | −0.18665 | −0.18855 |

**The impurity is conservative**: without PerturbMatch the best law is *better* (1.0257 vs
1.0159), so keeping PM in the main sweep understates the alternative by ~1% of cosine and cannot
be what produces the negative.


### 2.7 212's 96.2% exchange, re-checked under the winning law (its own instruction)

212 §6 closed `WHITEN_GAMMA` on a structural argument — along the gamma axis the operator trades
bulk cosine for retrieval at 96.2% of the board's own rate, so gamma is a knife-edge — and
flagged that *"that ratio ... need not be law-invariant"* and should be re-checked under any law
that wins. `gamma213.py` runs the same 5-point gamma sweep on the same four folds under both
laws (γ 0.60 → 0.70, a15 render):

| law | board cosine at γ 0.65 | d(overall `pds`)/d(cos) | + the board's `mse` rate +0.3655 | what survives |
|---|---:|---:|---:|---:|
| `logc20` (ships) | 0.077305 | **−0.3069** | **+0.0586** | 16.0% of the `mse` leg |
| `logc40` | 0.078534 | **−0.4092** | **−0.0437** | **−11.9%** — the trade is net *negative* |

**The ratio is NOT law-invariant, and it moves the wrong way.** Under the family's cosine
maximiser the operator charges *more* retrieval per unit of bulk cosine than the board pays, so
gamma stops being a knife-edge and becomes strictly unpayable. Adopting the alternative law
would close `WHITEN_GAMMA` harder, not reopen it.

*(My incumbent-law rate, −0.3069, is 13% smaller in magnitude than 212's −0.3517 because 212
transferred the level through its own `cd4f034`/h1 routes and reported both renders while this
uses the a15 render and `chal213`'s own solved level. The sign, the ordering and the conclusion
— the two legs nearly cancel, tilting toward lower gamma — are the same, and gamma is 212's
territory: **its number is the one to use.**)*

---

## 3. Phase 3 — six members, and the gate

### 3.1 `mse` through 096's `claim` route, on the challenge panel, at the re-solved α

`chal213.py` measures `A`, `B`, `claim`, `Kc`, ρ and the ρ-cliff per context under each law's own
render, using `cell_eval2`'s own `_across_pert_budget` / `_row_weights` / `jackknife_correction`
and 096's control statistics. Its `logc20` row reproduces 212's exactly (§0, gate 15).

| law | `A` | `B` | `claim` | ρ (A/B/C) | `s_cliff` (A/B/C) | `rA` |
|---|---:|---:|---:|:---:|:---:|---:|
| **`logc20` (ships)** | **0.231662** | **0.230372** | **0.246076** | 0.9926/0.9407/0.8821 | 1.0037/1.0310/1.0648 | 1.000 |
| `logc40` | **0.417332** | 0.414919 | 0.245772 | 1.000/1.000/1.000 | 0.7575/0.7611/0.7892 | **1.8015** |

> 🔴 **The transport law is an amplitude lever in disguise, and this is the load-bearing number
> of the experiment.** A 1.6% direction gain arrives bundled with an **80% inflation of the
> delivered log-pseudobulk energy**, because the recipient half of the law multiplies the
> delivered effect by `(p_rec + c)/(p_rec + c0)` — a factor of up to `c/c0 = 2` on
> low-expression genes (§4). The pre-render amplitude barely moves
> (`after_all_median_norm` 4.2293 → 4.2592, +0.71%); the **delivered** amplitude nearly doubles.

Consequences, both exact through the identity:

| | `logc20` | `logc40` |
|---|---:|---:|
| board cosine (transferred by the scale-invariant PLOLO ratio) | 0.077305 | 0.078534 |
| **raw `mse` at the shipped α = 0.982** | 0.876396 | **1.498987** |
| Δ overall on `mse` at the shipped α | 0 | **−0.10921** |
| **re-solved α** | **0.9856** | **0.7438** |
| raw `mse` at the re-solved α | 0.873262 | 0.869112 |
| **Δ overall on `mse` at the re-solved α** | **0** | **+0.00073** |

The optimum sits exactly on context A's ρ cliff in both laws (096 (b)) — `s*` equals
`s_cliff(A)` to the grid's resolution in each — and the incumbent's own re-solve lands on
**α = 0.98560** against 212's **0.98567** on the same route (rel 7e−5). α is 211's territory and
its 0.9808 on `rc_callset`'s profile is the number that ships; nothing here depends on which,
because every row is priced against the same baseline at *its* own α.

The measured local exchange rate along the law axis is
**d(overall `mse`\*)/d(cos) = +0.5922** — *more* generous than the +0.3655 the gamma axis pays
(212 §3.4), because the law moves `A` and the cliff together. Every `mse` leg in
`work/summary213.txt` outside `logc20`/`logc40` is a linear extrapolation at that rate and is
labelled as such; the two exact points are the ones above. **The extrapolation is only
meaningful near the incumbent** — for `mult`, `rawadd` and the rank laws it is a linear
extension across a 40–80% cosine collapse and should not be read as a measurement; for those
arms the `pds` leg alone (−0.029 to −0.090 overall) is decisive several times over.

### 3.2 `nmae`, `reach`, `fid`, `jac` on h1-test, at each law's own α

`h1sweep213.py` — 170's own six-member pipeline, the law replacing gamma, the render replaced by
the law's. `logc20` reproduces 170's γ = 0.65 row on 11 scalars at rel 0.00e+00. Donors:
`k562_gwps` + `orion_hct116` + `orion_hek293t`; `data/donors/h1_2025` covers 85/85 and is NOT
used. `nmae`'s pool is the reference's own BH < 0.05 set, median `n_real` **90** (188's ~89).

| arm | α | `mse` | `pds` | `nmae` | `reach` | `fid` | `jac` | `n_pred` | `A` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `logc20`, seed 0 | 0.982 | 1.01597 | 0.87535 | 0.95988 | 0.17752 | 0.50816 | 0.06425 | 3803 | 0.08349 |
| `logc20`, seed 1 | 0.982 | 1.01597 | 0.87535 | 0.95989 | **0.16794** | 0.50832 | 0.06430 | 3803 | 0.08349 |
| `logc40`, α unchanged | 0.982 | 1.37540 | 0.87633 | 0.95952 | 0.21381 | 0.51055 | 0.06460 | 4446 | 0.14306 |
| **`logc40`, seed 0** | **0.7438** | 1.01030 | 0.87703 | 0.96152 | 0.19658 | 0.50737 | 0.06354 | 3754 | 0.08221 |
| **`logc40`, seed 1** | **0.7438** | 1.01030 | 0.87703 | 0.96153 | 0.19947 | 0.50777 | 0.06337 | 3761 | 0.08221 |
| `logc20` at `logc40`'s α (control) | 0.7438 | 0.93951 | 0.87633 | 0.96818 | 0.17226 | 0.50331 | 0.06294 | 3241 | 0.04804 |

**Read as deltas only** (rule 9: h1-test's `n_real` is 766 mean against the challenge's ~89
median, so `fid`, `jac` and `reach`'s normalisation do not transfer; and 162 measured `reach` at
**2 of 6 signs — worse than chance** — as a board predictor):

| member | Δ raw, `logc40` @ its own α − `logc20` @ its own α | span | Δ scaled | Δ overall |
|---|---:|---:|---:|---:|
| `nmae` | **+0.00164** (worse) | 0.60 | −0.00273 | **−0.00046** |
| `fid` | −0.00067 | 0.30 | −0.00223 | **−0.00037** |
| `jac` | −0.00082 | 0.365 | −0.00224 | **−0.00037** |
| `reach` | **+0.02529** | 0.895 | +0.02826 | **+0.00471** — *unpriceable* |

The two-seed spread on `reach` at fixed law is **0.0096 raw** for `logc20` and **0.0029** for
`logc40`, so +0.025 is ~2.5 seed-spreads on n = 2 per arm — suggestive, not resolved, and on
the one member 162 measured as anti-predictive. **It is reported and then set aside, and the
verdict is stated both with and without it.**

### 3.3 Rule 13 — the class, the toll, the floor, and the answer

**Class: PROFILE.** The law changes the delivered profile itself, not the emission law. Rule 13's
table is single-valued on this: profile-class arms carry **no `reach` toll** (169 falsified it at
+0.0006). The −0.0146 toll belongs to emission-law arms that change the delivered profile, and
this is not one.

**Floor: +0.0089 (amplitude-changing), and there is no reading on which it is not.** `rA` =
**1.80**, α must move **−24.3%**, and at the shipped α the arm loses 0.109 overall on `mse`
alone. 211's standing rule is respected: `claim` is at its ceiling in context A, so no arm here
is re-priced on amplitude grounds — the α re-solve here is a *cost* the arm must pay, not credit
it may claim.

```
predicted = PLOLO d_pds/6 * 1.011  +  096 d_mse/6  -  0 (profile class)
SHIP only if predicted > +0.0089
```

| arm | `mse`\* | `pds` (a23) | `nmae` | `fid` | `jac` | `reach` | **TOTAL** | floor | verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `logc40`, three priceable members | +0.00073 | +0.00030 | −0.00046 | — | — | — | **+0.00057** | +0.0089 | **VETO** |
| `logc40`, five members (h1 `reach` excluded) | +0.00073 | +0.00030 | −0.00046 | −0.00037 | −0.00037 | — | **−0.00017** | +0.0089 | **VETO** |
| `logc40`, all six, crediting h1's `reach` in full | +0.00073 | +0.00030 | −0.00046 | −0.00037 | −0.00037 | +0.00471 | **+0.00454** | +0.0089 | **VETO** |
| `logc36` (the cosine maximiser) | +0.00072 | +0.00028 | — | — | — | — | **+0.00100** | +0.0089 | **VETO** |
| `logc32` | +0.00066 | +0.00030 | — | — | — | — | +0.00095 | +0.0089 | **VETO** |
| `bc0.05` (curvature maximiser) | +0.00010 | +0.00017 | — | — | — | — | +0.00026 | +0.0089 | **VETO** |
| `mult` (026's law) | −0.01886 | −0.02902 | — | — | — | — | −0.04788 | +0.0089 | **VETO** |
| `rawadd` / `bc1` | −0.03227 | −0.08997 | — | — | — | — | −0.12224 | +0.0089 | **VETO** |

**The best arm in the whole family, under the single most generous reading available anywhere in
this experiment, is +0.00454 — and that reading credits, in full, a `reach` delta on the one
surface rule 9 forbids reading it on and the one member 162 measured at 2/6 signs. Under the
defensible reading it is +0.00057. The floor is +0.0089.**

---

## 4. The mechanism: a transport law is a per-gene, sign-preserving gain

To first order, with `dp = p_pert − p_ctrl`,

```
g_c        ≈  dp / (p_don + c)                             the donor half
dh         ≈  g_c · (p_rec + c) / (p_rec + c0)             the recipient half
so  dh_c / dh_c0  ≈  [(p_don + c0)/(p_don + c)] · [(p_rec + c)/(p_rec + c0)]        (*)
```

`mech213.py` measures both halves against (*) directly:

| half | statistic | measured |
|---|---|---|
| **donor** (`logc40`, `orion_hct116`, 400 rows) | corr(log measured ratio, log closed form) | **0.99999** |
| | median \|relative error\| | **0.0004** |
| | **sign flips** | **0.0000% of 7,013,224 live entries** |
| **recipient** (same effect matrix rendered two ways) | corr, contexts A/B/C | **0.9903 / 0.9922 / 0.9928** |
| | per-gene gain, median [p05, p95] | **1.24 [1.06, 1.59]** |
| | sign flips | 0.178% / 0.201% / 0.173% |
| | delivered row-norm ratio | ×1.263 / 1.278 / 1.275 |
| **both, plus the operator and CD4 stages** | row cosine to the shipped delivered effect | **0.9867 / 0.9850 / 0.9854** |
| | sign flips | 3.24% / 3.20% / 3.26% |
| | delivered row-norm ratio | ×1.319 / 1.344 / 1.336 |

**Three consequences, and they are the reusable part of this experiment.**

1. **A gain cannot reorder anything.** The law's own map flips **zero** signs on the donor side
   and 0.2% on the recipient side; the 3.2% that appear end-to-end are the operator and CD4
   stages reacting to a re-weighted input, not the law. So **the transport law is structurally
   incapable of buying the ORDERING channel** — the same objection 184 raised against the
   Tahoe cross-line operator (*"the winning family is a per-gene gain and a gain is
   sign-preserving"*), arriving here from the opposite direction and confirming it.
2. **The law can only ever pay for donor↔recipient basal MISMATCH.** (*) is exactly 1 wherever
   the two lines express a gene at the same level. There is plenty of mismatch to work with
   (log-CPM corr 0.34–0.49, only 38–45% of genes within 2×) — and the optimal exploitation of
   all of it is **+1.6% of cosine**.
3. **The delivered amplitude must move**, because the second bracket of (*) exceeds 1 for every
   `c > c0`. That is why an alternative law is an amplitude-changing arm under rule 13, and it
   is not avoidable by construction.

### 4.1 The question five experiments left open, answered

160, 161, 174, 175 and 104 all failed to make the recipient's basal profile useful **as a
feature**. None asked whether it was already entering **wrongly**. It is not — and it is
already entering **hard**:

* the shipped law's optimum in `c` is within **1.6% of cosine** of the family maximum, from
  **both** parametrisations;
* deranging the recipient's basal costs the shipped law **16% of its transfer cosine and 0.070
  raw `pds`** (§5, F2). The recipient's basal is load-bearing in the render already, and it is
  worth far more there than any of the five feature attempts recovered.

---

## 5. Phase 4 — falsifiers

### F1 — the incumbent is an interior point and reproduces the shipped numbers exactly

Interior in both families (index 7 of 22 in `c`, index 6 of 16 in `q`). At its own parameter it
reproduces: `E062.prepare_main_source` **bit-exactly**; 204's PLOLO `equal` row on 4 folds at
**rel 0.00e+00** (from `logc20` AND from `bc0`); a15's `baseline_pre` and `candidate_pre`
**sha256**; 170's two effect matrices **bit-for-bit**; 170's h1-test γ = 0.65 row on 11 scalars
at **rel 0.00e+00**; and 212's challenge-panel `claim`/ρ/`X`. **PASS.**

### F2 — the permutation control. **The one to watch, and it passes decisively.**

The recipient's control profile is gene-permuted before rendering, with the CPM gate left on the
TRUE basal so both arms carry effect on exactly the same genes. Three seeds, macro over four
folds, a15 render:

| law | true cosine | deranged cosine (3 seeds) | ratio | true `pds` | deranged `pds` |
|---|---:|---:|---:|---:|---:|
| **`logc20` (ships)** | 0.0218418 | 0.0183406 ± 0.0001380 | **0.840** | 0.740297 | 0.670467 |
| `logc5` | 0.0202174 | 0.0161236 ± 0.0000513 | 0.798 | 0.736693 | 0.704559 |
| **`logc40`** | 0.0221890 | 0.0178985 ± 0.0001614 | **0.807** | 0.739542 | 0.649604 |
| `logc100` | 0.0210512 | 0.0164158 ± 0.0001652 | 0.780 | 0.733272 | 0.632924 |
| `mult` | 0.0128440 | 0.0087921 ± 0.0002088 | 0.685 | 0.658110 | 0.639794 |
| `rawadd` | 0.0064444 | 0.0066007 ± 0.0001912 | **1.024** | 0.553647 | 0.583543 |

```
cosine ratio logc40 / logc20 :   TRUE 1.0159   DERANGED 0.9759
raw pds      logc40 − logc20 :   TRUE −0.000756  DERANGED −0.020863
```

**The gain does not merely fail to survive the derangement — it inverts.** So `logc40`'s +1.6%
is genuinely reading the recipient's own gene-level basal; it is transport, not amplitude.
**209's outcome — the random control as the maximiser over all its arms — does not occur here.**

Two by-products worth keeping:
* **the shipped law is not basal-free.** Deranging the recipient costs it 16% of cosine and
  0.070 raw `pds`, because `profile_from_effect` renormalises the profile and `expm1` is convex,
  so a wrong basal↔effect pairing changes the mass correction on every gene. §4.1.
* **`rawadd` is the one law that gets BETTER with a deranged recipient** (ratio 1.024). Additive
  transport reads the recipient's basal actively harmfully, which is the cleanest statement of
  why the `c → ∞` limit is the worst arm on the board.

**F2b — the donor-side derangement.** The donor's own basal is gene-permuted inside `to_law`.
This branch does not carry the PerturbMatch slot, so it must be read against the `--no-pm`
column (§2.6), where it reduces to a proper no-op at the incumbent:

| law | no-PM cosine | donor-deranged cosine | ratio |
|---|---:|---:|---:|
| `logc20` | 0.0212821 | **0.0212821** | **1.0000** — exact no-op, as it must be |
| `logc40` | 0.0218294 | 0.0197058 | 0.9027 |
| `logc100` | 0.0208582 | 0.0160742 | 0.7707 |
| `rawadd` | 0.0057924 | 0.0042942 | 0.7414 |
| `logc5` | 0.0202897 | degenerate (`nan`) | — |

`logc20`'s no-op is exact because at `c = c0` the law does not read the donor's basal at all;
that is the check that the control is a control. `logc5`'s donor derangement is degenerate — a
wrong small basal in the denominator produces unbounded effects, which is 026's mechanism
appearing as a failure mode rather than as a loss.

### F3 — the leak-free fold, reported

`kolf2` is not a production donor and contributed no row to any production operator, so that
fold uses the exact shipped operator (`e2f65ecf…`) with no leakage. It is carried in every table
above rather than folded away, and it agrees with the other three at every `c` in [24, 60]: per-fold cosine ratios at `c` = 40
are **1.0226 / 1.0082 / 1.0126 / 1.0229** (k562 / hct116 / hek293t / **kolf2**), so the
leak-free fold is one of the two *most* favourable to the alternative law and the negative is
not resting on it. The maximiser's location moves a little by fold — `c` = 45–50 on `k562_gwps`
and `kolf2`, `c` = 32–36 on the Orion pair — and every fold's maximum ratio is under 1.03. No absolute is quoted from the Orion pair.

### F4 — the `mse` identity, priced not asserted

The law changes the delivered profile, so `mse` cannot be claimed invariant, and it is not:
§3.1 prices it through 096's `claim` route on the challenge panel with `cell_eval2`'s own
primitives, at the shipped α (**−0.109 overall**) and at the re-solved α (**+0.00073**). The
delivered per-gene panel totals are **not** identical — `rA` = 1.80 — so `rc_hybrid`'s and
`rc_callset`'s digest-identity exemption does not apply and was never claimed.

### F5 — the two families' additive limits agree bit-for-bit

`rawadd` (`p' = p_rec + (p_don_pert − p_don_ctrl)`, written directly) and `bc1` (the same map
reached through the Box-Cox parametrisation, `((p+c0)^q − 1)/q` at `q = 1`) are computed by
different code and produce **identical floats** on all four folds and all six reported scalars:
cosine 0.0064444, `pds` 0.553647, a23 cosine −0.0018205, a23 `pds` 0.500779. A parametrisation
error in either family would have shown up here.

---

## 6. Verdict, in overall units

| arm | Δ overall, best available reading | rule 13 floor | verdict |
|---|---:|---:|---|
| `logc40` (cosine maximiser), 3 priceable members | **+0.00057** | +0.0089 | **VETO** — short by 16× |
| `logc40`, all six, crediting h1's unpriceable `reach` | **+0.00454** | +0.0089 | **VETO** — short by 2.0× |
| `logc36` | +0.00100 | +0.0089 | **VETO** |
| `bc0.05` (curvature maximiser) | +0.00026 | +0.0089 | **VETO** |
| `logc10` (toward log-fold-change) | −0.00231 | +0.0089 | **VETO** |
| `mult` — pure log-fold-change (026's law) | −0.04788 | +0.0089 | **VETO** |
| `rawadd` ≡ `bc1` — additive in CPM | −0.12224 | +0.0089 | **VETO** |
| `rank_u` — quantile transport | −0.12090 | +0.0089 | **VETO** |

**NO-GO. No slot. The transport law stays `log1p(TS·p)`.**

**And the reason is structural, not an operating point.** §4's identity says a transport law is
a per-gene gain `[(p_don+c0)/(p_don+c)]·[(p_rec+c)/(p_rec+c0)]`. Three things follow that do not
depend on `c`, on the emission, on the dose or on α:

1. it is **sign-preserving**, so it cannot touch the ordering channel at all;
2. it is **1 wherever donor and recipient agree**, so its whole budget is donor↔recipient basal
   mismatch — and the optimal use of all of that mismatch is measured here at **+1.6% of
   cosine** against the **+31%** 197 says a direction arm needs;
3. its second bracket **exceeds 1 for every `c > c0`**, so any move toward the additive end is
   automatically an amplitude-changing arm and automatically faces the +0.0089 floor.

**A future agent who wants to reopen the transport law needs a mechanism that is not a per-gene
function of the two basal profiles.**

---

## 7. What this corrects, closes and hands on

**Corrects.**

- **The mandate's own parametrisation of family (6).** `c = 1/TS` is not the log-fold-change end
  of `log(p + c)`; it is the **shipped point**. Log-fold-change is `c → 0`. The incumbent is
  therefore interior by algebra, which is a stronger statement than the mandate hoped for.
- **026's magnitude.** Its `lfc` re-expression loses **−0.0822 raw `pds`** and **41% of the
  transfer cosine** on the delivered instrument, not the −0.0041 it measured on the
  effect-level LODO surface. 156's *"every effect-level screen before this one scored an object
  the board never sees"* applies to 026 and the correction is a factor of 20.
- **The framing that the recipient's basal is unused.** It is heavily used — through the render,
  not as a feature — and deranging it costs the shipped law 16% of its transfer cosine and 0.070
  raw `pds`. 160/161/174/175/104 failed to add *more* basal information, not to add the first.

**Closes.**

- **The transport law**, on the cosine endpoint, at each law's own α, on the production chain,
  with two independent one-parameter families through the incumbent, two structural limits, four
  named alternatives, four PLOLO folds, and a permutation control that the winner fails.
- **The render**, as the same object: `profile_from_effect`'s `expm1(b + h)` is the `c = c0`
  member and every alternative render in the family loses or ties.

**Hands on.**

- **The identity (*) is the reusable object.** Any future arm that changes how a donor effect is
  expressed relative to expression level is a per-gene gain, is sign-preserving, and is bounded
  by donor↔recipient basal mismatch. That is a closed-form veto, not a screen.
- **A transport change is an amplitude arm.** `rA` = 1.80 for a 1.6% direction gain. If anyone
  proposes one again, the amplitude floor applies and the α re-solve is a cost.
- **`d(overall mse*)/d(cos) = +0.5922` along the law axis**, against +0.3655 along the gamma
  axis (212 §3.4). The law axis pays *more* per unit of cosine — it simply has almost no cosine
  to sell.
- **The staging harness.** `stage213.py` injects an arbitrary transformation of the donor tables
  into the *production* builder and reproduces a15's `candidate_pre` sha256 and 170's effect
  files bit-for-bit at the identity. Any future donor-side arm can be priced on all six members
  through it in about 12 minutes.

---

## 8. What was NOT measured, and why

* **`chal213`'s exact `A`/`B`/`claim` for laws other than `logc20` and `logc40`.** Each costs a
  full production staging (~6 min); the two measured points bracket the maximiser and every
  other `mse` leg in `work/summary213.txt` is a linear extrapolation at the **measured** local
  rate +0.5922 per unit cosine, labelled as such. The best arm is 8.6× under the gate, so a
  third point would not change the verdict.
* **h1-test's six members for laws other than `logc20` and `logc40`.** Same reason; and the
  three DE members are unreadable in absolute terms on that surface (rule 9) and unpriceable as
  board predictors (162: `reach` 2/6 signs).
* **The rank laws in a23 geometry.** Quantile transport has no analytic density, so the
  `minnorm` mass fix cannot be solved along its own gradient. They are reported in a15 geometry
  only, which is the geometry the incumbent ships.
* **PerturbMatch and the CD4 residual were not re-expressed.** Neither is a `log1p(TS p)`
  difference — PerturbMatch is a limma-voom LFC and the CD4 bank is a log2 fold change — so
  neither has a law-space form that is not an invention. Both are left exactly as production
  has them, which is also what keeps the identity row exact; F6 prices the PerturbMatch half
  and finds it **conservative**.
* **`WHITEN_GAMMA`, `VCC_ALPHA` and `VCC_CD4_TARGET_F` were not swept.** 212, 211 and 210 own
  them; the α re-solve here is a consequence of the law, not a free variable. Territory
  respected.
* **No submission was built, no slot spent, nothing written to `submissions/`, `cell-eval2`
  never run, `vcc` never run.** Every scored number comes from `cell_eval2`'s own primitives or
  from the production builder's own functions, anchored to a real scorer run by gate 14.
* **Files touched outside this directory: none.** 162's fold caches and 170's DE-reference cache
  were **cloned** into `work/` before `plolo.WORK` and `h1kit.WORK` were redirected, precisely
  so a cache miss could not write into a sibling's directory. Importing 204's `score204.py` and
  170's `stage`/`gate`/`h1kit`/`h1sweep` refreshes their `__pycache__/*.pyc`; no source, data or
  result file outside `experiments/213-transport-law/` was written.

**Peak RSS, and the one overrun.** `stage213.py` 7.12 GiB, `plolo213.py` ~7.5, `falsify213.py`
~6.5, `h1sweep213.py` ~4.5, `prep213.py` 4.07, `chal213.py` ~4.0, `gamma213.py` ~4.0,
`mech213.py` ~3.0, `eff213.py` ~5.5 — all comfortably inside the ceiling, and no two were run
concurrently above 10 GB combined.

`refit213.py` is the exception and it is recorded plainly. A first invocation asked for **two**
operator refits in one process and reached **15.53 GiB**: on macOS the pages of the 4.2 GB
float64 basis memmap stay resident after the file is unlinked, so the second law's basis build
started from the first law's high-water mark. **It was killed as soon as the reading was seen**,
and the run was repeated one law per process, where it peaks at **12.27 GiB** — still 2.3% over
the stated 12 GB ceiling, entirely inside `np.linalg.eigh` on a 17,942² float64 covariance
(2.57 GB in, 2.57 GB out, ~2× that in LAPACK workspace), which is irreducible at this basis
size. The three rows the first invocation had already **completed and written** before the
overrun — `logc20`/shipped, `logc20`/refit and `logc40`/shipped, all of which peaked at
10.71 GiB — are the ones reported; the overrun happened while assembling `logc40`'s basis, and
that row (`logc40`/refit) is the one the second, single-law invocation produced.

---

## 9. Files

| file | what it does |
|---|---|
| `law213.py` | the transport-law family: two one-parameter families through the incumbent, four named alternatives, two rank transports, the law-aware renderer and the law-aware mass fix |
| `prep213.py` | the donor-side sufficient statistics under a law (`E062.prepare_main_source`'s three streams, law-injected), cached per (donor, law) |
| `plolo213.py` | **the primary instrument**: 162's PLOLO with the law as an argument, 4 folds × 43 laws, both renders, per-target cosines → `work/plolo213_main.json` |
| `boot213.py` | 204's 4,000-draw paired target bootstrap on Δcos and Δ`pds` → `work/boot213_{a15,a23}.json` |
| `gate213.py` | Phase 0 rows 1–10: the operator table, the effect definition, 026/169's geometry, the interiority proof, 197/162/204's anchors, 212's CD4 dose |
| `stage213.py` | the law injected into the **production** builder's donor loading; a15's `baseline_pre`/`candidate_pre` sha256 gate |
| `eff213.py` | the shipped operator + CD4 residual at γ = 0.65, λ = 0; a15's six audit scalars; bit-identity against 170's effect files |
| `chal213.py` | the challenge panel's `A`/`B`/`claim`/ρ/`s_cliff` per law, the exact α re-solve, and the `mse` leg through 096 |
| `h1sweep213.py` | 170's six-member h1-test pipeline with the law replacing gamma and the law's own render |
| `mech213.py` | the mechanism: both halves of the law against the closed form, the sign-flip rate, the basal mismatch |
| `falsify213.py` | F2 and F2b — the deranged recipient and donor basal controls |
| `refit213.py` | the whitening operator refit in a law's own space, and the fold scored with both |
| `gamma213.py` | 212 §6's cosine-for-retrieval exchange, re-measured under each law on the gamma axis |
| `work/summary213.txt` | the consolidated priced table: all 43 laws, cosine, ratio, bootstrap CI, both renders' Δ`pds`, and the six-member legs |

**Inputs read, never written:** `data/donors/{k562_gwps,orion_hct116,orion_hek293t}`,
`data/donors/whitener_union_pruned.npz`, `data/external/{perturbmatch,kolf2,cd4t}`,
`data/challenge/*`, `data/harness/h1-test.h5ad`, `experiments/062-cd4-post-whiten/`,
`experiments/064-pm-cd4-composition/run.py`, `experiments/096-mse-floor/`,
`experiments/153-de-bench/debench.py`, `experiments/162-retrodicting-validator/`,
`experiments/170-gamma-mse/`, `experiments/204-repression-weighted-mixing/`,
`scripts/make_submission_pm_cd4.py`, `src/vcc2026/`. No production code was touched.

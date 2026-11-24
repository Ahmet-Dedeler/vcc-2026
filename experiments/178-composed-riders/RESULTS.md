# 178 — the sub-threshold riders, composed and priced as one production arm: **NO-GO. The rider family is closed.**

Date: **2026-09-02**, coordinator-run. Measurement only. No `cell-eval2`, no `vcc`, no count
build, nothing written into `submissions/`. Emission law held at a15's throughout
(`--tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 --overdispersion 4.0`), mass fix
`minnorm` as a23 shipped it, `lam` at 0. Peak RSS 9.81 GiB (the effects stage; ~5.8 GB of that
is memory-mapped donor `delta.npy` pages).

Two coordinator INBOX messages arrived mid-run and are answered here: **170's gamma verdict**
(`INBOX-URGENT-drop-gamma.md`) forced the `mse` route to be redone and added the gamma-free
composite `flip + w211`, which was built and priced; **175's fourth rider**
(`INBOX-fourth-rider.md`) is not priced, and §11 gives the quantitative reason it cannot rescue
the family.

---

## Verdict in eight lines

1. **Every rider and every composite is negative on all-six, by −0.005 to −0.030 overall,
   against rule 13's +0.0089 gate.** The best arm in the family is sign repair alone at
   **−0.0095**; the three-way composite is **−0.0298**.
2. **The `pds` costs compose SUB-additively (ratios 0.52–0.81) while the `mse` costs compose
   ADDITIVELY to SUPER-additively (0.955–1.226) and the DE gains compose sub-additively
   (0.29–0.84).** That is the scientific result, and it is the exact opposite of the asymmetry
   the experiment was opened on: **the gains discount when composed and the costs do not.**
3. **The binding mechanism is 170's ρ cliff, and it now has three more measurements.** a15 sits
   at ρ = 0.992 / 0.940 / 0.881 on contexts A/B/C. **Every rider except sign repair raises the
   effect's mean-square `A` by 11–28% and pushes ρ to 1.000 in all three contexts**, where
   `min(claim, B)` saturates and the extra amplitude is pure cost. Three composites clamp
   `score_mse` to its floor outright.
4. **Sign repair is the only norm-preserving rider, exactly as 170's criterion predicts**
   (`A` 0.23134 vs a15's 0.23152, ρ unchanged), and it is the only one with a **positive**
   `mse` — +0.0021, and +0.0022 under the other amplitude convention, so the sign is
   convention-free. It still dies, on `pds`.
5. **On the board-validated `pds` instrument every rider costs far more than the h1-test
   numbers that motivated it** — sign repair −0.0245 raw where 173 read −0.0042.
6. **The donor weighting dies on a mechanism: whitening consumes the sign gain that motivated
   it.** 175's +0.0043 for (2,1,1) over (1,1,1) is measured here at **−0.0012 after the
   operator**, and the arm's own `mse` cost is −0.0108.
7. **New capability, kept and proven safe.** `VCC_DONOR_WEIGHTS` and `VCC_SIGN_REPAIR` exist,
   default to the shipped behaviour, and the patched builder at a15's config reproduces the
   **incumbent a23 build's effects audit bit-for-bit, 9 fields of 9, including
   `effects.final_sha256`**. 46 unit tests pass (40 pre-existing + 6 new).
8. No slot spent. a23 unchanged.

---

## 1. What was built, and the rule-10 proof that the defaults are safe

Three switches in `scripts/make_submission_pm_cd4.py`:

| switch | default | what it does |
|---|---|---|
| `VCC_WHITEN_GAMMA` | `0.65` (already existed) | the operator's exponent |
| `VCC_DONOR_WEIGHTS` | unset → the literal `{name: 1.0}` | comma-separated weights in `PRODUCTION_DONORS` order, e.g. `2,1,1`; `MultiDonorModel.effect` renormalises over the donors that carry each target, so the values are relative |
| `VCC_SIGN_REPAIR` | unset → no-op | `cons`: on the operator's columns, wherever `sign(whitened) == -sign(pre-whitening consensus)` and both are non-zero, take the consensus's sign and keep the magnitude. Applied **after** `apply_v9_whitening`, **before** the CD4 residual |

Each emits its audit key only when it is on (like `mass_fix`), so an unflagged build's `meta` is
byte-identical to every prior build's rather than differing by a key.

**Rule 10, the effects stage reproduced with the patched builder** (`effects178.py`; the effects
stage is where all three riders live, and it is where the audit rule 10 diffs is written):

```
=== rule 10: the a15-config arm vs a23-massfix (incumbent, built HERE) ===
  effects.candidate_pre_sha256   0601336fdf1ec5c5  vs 0601336fdf1ec5c5  EXACT
  effects.final_sha256           22c09f22c749bda8  vs 22c09f22c749bda8  EXACT
  whitening.after_all_median_norm        4.229295604454299   vs 4.229295604454299   EXACT
  whitening.after_operator_median_norm  25.73792992728808    vs 25.73792992728808    EXACT
  whitening.before_operator_median_norm  4.122892101713494   vs 4.122892101713494   EXACT
  whitening.operator_preserve_scale      0.16018740098217019 vs 0.16018740098217019 EXACT
  cd4.fraction                           0.09472706196346513 vs 0.09472706196346513 EXACT
  cd4.whitened_base_median_row_norm      4.229295604454299   vs 4.229295604454299   EXACT
  operator.sha256                        e2f65ecf69b7e0b0    vs e2f65ecf69b7e0b0    EXACT
  EXACT-match differences: 0  (rule 10 clean)
```

Against **a15**'s shipped meta the same run reads **6 differences**, and they are not the patch.
a15 was built on AWS (`/work/exp089-.../input/whitener_union_pruned.npz`) and a23 on this
machine; **the two shipped metas already disagree with each other in exactly those six fields**
(`after_all_median_norm` 4.2292956044543 vs 4.229295604454299, `operator_preserve_scale` …024 vs
…019, `post_whiten_sha256`, `final_sha256`, `cd4.fraction`, `cd4.whitened_base_median_row_norm`)
**while agreeing exactly on `baseline_pre_sha256` and `candidate_pre_sha256`**, which localises
the difference to the operator matmul's summation order on different hardware, upstream of
anything this patch touches. 170 independently hit and documented the same arm64-vs-x86 effect.

The operator's **column set** was rebuilt from the model and sha-checked:
`b2cae0d381c3cd70…` = a15's `operator.columns_sha256`, 17,942 columns — the check rule 10's own
worked example (a19's silent `--whitener-id` default) says the frozen-input hash table cannot
make. Operator by path and sha (rule 11, never a description):
`data/donors/whitener_union_pruned.npz`, sha256
`e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123`, `basis_prune: essential`.
h1-test is H1 hESC and contributed no row to its basis; `data/donors/h1_2025` never enters.

**Per-arm difference count against the a15 arm's own written config** — exactly the intended
number in every case:

```
a15         0 []                                   effect matrix moved: False
g045 / g030 1 ['gamma']                            True
flip        1 ['sign_repair']                      True
w211        1 ['weights']                          True
flip_w211   2 ['weights', 'sign_repair']           True
g045_flip   2 ['gamma', 'sign_repair']             True
g045_w211   2 ['weights', 'gamma']                 True
all3_*      3 ['weights', 'gamma', 'sign_repair']  True
```

**Unit tests** (`uv run python scripts/test_make_submission_pm_cd4.py`): **46 pass**, peak RSS
0.183 GiB. Six are new and each asserts both directions:

- `donor_weights` OFF is the literal shipped dict; ON is the requested mix **keyed by position**
  (a `1,1,3` row asserts HEK293T got the 3, so a mis-assignment cannot pass); `(4,2,2)` lands
  bit-identically on `(2,1,1)`'s effect (the weights are relative); malformed, negative and
  all-zero values fail closed at import; the empty string is the OFF state.
- `sign_repair` OFF adds no meta key and leaves the effect matrix bit-identical; ON preserves
  `|effect|` **exactly**, touches no off-axis column, moves **exactly** the disagreeing entries,
  leaves planted zeros alone, is idempotent, and — the check that would fail if the mask were
  empty — is asserted **not** to be a no-op and to land on `sign(consensus)` wherever both are
  non-zero.

---

## 2. Instruments, and every known-answer row

| channel | instrument | known-answer row | result |
|---|---|---|---|
| `pds`, board-validated | PLOLO, challenge frame, 4 folds (`plolo178.py`) | 162's own `chal{A,B,C}/delivered` endpoints for `a15` and `a23`, 4 folds × 2 arms × 3 contexts | worst \|Δpds\| **0.00e+00** |
| — | the weighted mix at equal weights | must be `E062.build_baseline` **bit-for-bit** | max \|Δ\| **0.0**, all four folds |
| — | the k562-held fold | `w211` must equal `a23` **exactly** (no K562 to up-weight) | exact, 6 digits, on all three contexts |
| `mse`, real panel | 170's per-context identity (`chal178.py`) | a15 arm: `claim` **0.24610**, `A` **0.23152**, `B` **0.23031**, mean ρ **0.9381** | **+0.00%** on all four |
| — | — | `g045` arm: `A` **0.25669**, `B` **0.25539**, X ratio **1.0889** | +0.00 / +0.00 / **+0.11%** |
| — | — | `g045` raw `mse`, 170's independently measured **0.90965** | **0.90789**, −0.19% |
| `pds`/cos/sign, h1-test | delivered space, 168d's pipeline (`h1rebuild178.py`) | `signfix173.py` re-run: base `pds` **0.86134**, cos 0.03531, flip −0.00420, null −0.03389 | reproduced to the digit |
| donor machinery + weights | unscaled, un-whitened sign accuracy | 175's rows: K562 **0.5668/0.5670**, HCT116 **0.5366/0.5375**, HEK293T **0.5292/0.5308**, (1,1,1) **0.5697**, (2,1,1) **0.5740**, (3,1,1) **0.5742** | all **PASS** |
| whitening's sign anatomy | the rebuild's flip mask | 173: 13.9% of gated entries, 1.6% of energy, whitened sign right 0.4556 / consensus 0.5444 | **14.5%**, 2.31%, **0.4573 / 0.5427** |
| DE members, h1-test | 173's realisation (`realize178.py` + `six178.py`) | `KA_a15_raw` median `n_pred` **6034**, `reach` **0.18283**; 173's `base_mf` precision **0.50789**, nmae **0.96658**, reach **0.18432** | all **PASS**, exact |

**The one check that is not clean, reported rather than laundered.** The h1-test consensus is
**rebuilt** from the three production donors and the shipped `union_pruned` operator, because
gamma and the donor weights are not expressible as transforms of a frozen whitened matrix. Its
per-row cosine against the effect recovered from the 095 artifact is **median 0.672 (p05
0.542)** — 097 §9b's same-config reconstruction got 0.99922. These are different objects (3
donors vs 2, `union_pruned` vs the 095 h1-test operator, alpha 0.982 vs 1.08, no PerturbMatch or
CD4 residual on the rebuild), and its median `n_pred` is 4,768 against the recovery's 6,034.
**Every h1-test number below is a ratio inside one reconstruction, never an absolute across
two.**

---

## 3. The `mse` route was rebuilt mid-experiment, and the first version was wrong

`mse178.py` (superseded, kept for provenance) priced each arm at its own ρ cliff
`s* = √(claim/B)`, i.e. it silently re-tuned `VCC_ALPHA` per arm and thereby **linearised away
the `min(claim, B)` saturation**. That is the exact error 170's own first pass made and
corrected, and it flipped the sign of the `mse` column for every arm except sign repair
(`g045` read +0.0025; it is **−0.0062**).

`chal178.py` replaces it with 170's identity, applied **per context** at the **shipped**
amplitude, which is what the board reads:

```
MSE_c − 1 = ( A_c − 2 X_c + Kc_c − min(claim_c, B_c) ) / sumD
```

`A`, `B` and `claim` are properties of our own submission and are measured per arm per context
on the real 300-target panel through the builder's own `profile_from_effect` and
`emission_spread_vector`. `X` is solved from **a23's own board point** at the a15 configuration
and carried across arms by the ratio measured on h1-test, recovered from `endpoints`' own
outputs without a second pipeline (`x = (1−sse)·d/gstar`, `d` common to all arms).

**Two independent stagings agree to the fourth decimal.** 170 and this experiment share no code
beyond the builder and `cell_eval2`:

| quantity | 170 | 178 |
|---|---:|---:|
| a15 `claim` | 0.246100 | **0.246102** |
| a15 `A` / `B` | 0.23152 / 0.23031 | **0.231516 / 0.230315** |
| a15 mean ρ | 0.9381 | **0.938096** |
| γ 0.45 `A` / `B` | 0.25669 / 0.25539 | **0.256690 / 0.255387** |
| γ 0.45 h1-test `X/X₀` | 1.0889 | **1.0901** |
| γ 0.45 raw `mse` | 0.90965 | **0.90789** |

### The ρ cliff, measured for every rider

| arm | `A` | `B` | `claim` | ρ (A / B / C) | `X/X₁₅` | raw `mse` |
|---|---:|---:|---:|:---:|---:|---:|
| **a15/a23** | 0.23152 | 0.23031 | 0.24610 | **0.992 / 0.940 / 0.881** | 1.0000 | 0.87234 |
| **`flip`** | **0.23134** | **0.23015** | 0.24610 | **0.991 / 0.940 / 0.881** | 1.0762 | **0.86041** |
| `g045` | 0.25669 | 0.25539 | 0.24612 | 1.000 / 1.000 / 0.977 | 1.0901 | 0.90789 |
| `g030` | 0.26890 | 0.26755 | 0.24602 | 1.000 / 1.000 / 1.000 | 1.1832 | 0.93828 |
| `w211` | 0.26218 | 0.26036 | 0.24589 | 1.000 / 1.000 / 0.992 | 1.0487 | 0.93371 |
| `flip+w211` | 0.26226 | 0.26045 | 0.24585 | 1.000 / 1.000 / 0.992 | 1.1415 | 0.91957 |
| `g045+flip` | 0.25667 | 0.25537 | 0.24608 | 1.000 / 1.000 / 0.977 | 1.1331 | 0.90129 |
| `g045+w211` | 0.28451 | 0.28261 | 0.24573 | 1.000 / 1.000 / 1.000 | 1.1861 | **1.00721** |
| `all3` (γ 0.45) | 0.28480 | 0.28290 | 0.24575 | 1.000 / 1.000 / 1.000 | 1.2411 | **0.99979** |
| `all3` (γ 0.30) | 0.29649 | 0.29454 | 0.24564 | 1.000 / 1.000 / 1.000 | 1.3155 | **1.03963** |

**`claim` is flat at 0.2456–0.2461 across every arm** — the emission policy is unchanged and
`min(Ĉ_p, C_p)` saturates, so no rider can move the emission channel. Derived, not assumed.

**Sign repair is the only rider that does not move `A`** (−0.08%), keeps ρ exactly where a15 has
it, and converts its whole +7.6% alignment gain into `mse`. **Every other rider raises `A` by
11–28%, crosses the cliff, and pays for its alignment with amplitude the metric refuses to
refund.** The three arms marked in bold on `raw mse` sit **above the `score_mse` clamp point
0.98786**, so their scaled `mse` is exactly 0 and their delta (−0.02026) is a **censored lower
bound**, not a measurement — 162 flagged the same clamp for a11.

This is 170's design criterion — *"a transform that raises `X` without raising `A`"* — tested on
three riders it was not derived from, and it separates them cleanly.

---

## 4. The six-member table, every arm, in overall board units

Deltas against the a15/a23 incumbent. `pds` from PLOLO (calibrated ×1.011, /6); `mse` from
170's identity at the shipped amplitude; `fid` from h1-test precision (no `n_real`, rule 9),
twin-averaged; `nmae` and `jac` from h1-test; `reach` from h1-test at native BH, twin-averaged,
×0.393.

| arm | `pds` | `mse` | `fid` | `nmae` | `jac` | `reach` | **TOTAL** | **TOTAL\*** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| a15/a23 (base) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | −0.00243 |
| `g045` | −0.00347 | −0.00624 | +0.00074 | +0.00065 | +0.00041 | −0.00123 | **−0.00914** | −0.01034 |
| `g030` | −0.00889 | −0.01157 | +0.00151 | +0.00112 | +0.00022 | −0.00156 | **−0.01918** | −0.02005 |
| **`flip` (sign repair)** | **−0.00913** | **+0.00209** | +0.00140 | +0.00044 | +0.00016 | −0.00034 | **−0.00538** | −0.00748 |
| `w211` | +0.00013 | −0.01077 | +0.00051 | −0.00017 | +0.00003 | −0.00022 | **−0.01049** | −0.01271 |
| **`flip + w211`** (gamma-free) | −0.00631 | −0.00829 | +0.00220 | +0.00060 | +0.00011 | −0.00131 | **−0.01298** | −0.01411 |
| `g045 + flip` | −0.01021 | −0.00508 | +0.00178 | +0.00078 | −0.00006 | −0.00126 | **−0.01405** | −0.01522 |
| `g045 + w211` | −0.00213 | −0.02026† | +0.00138 | +0.00037 | +0.00031 | +0.00095 | **−0.01938** | −0.02276 |
| **`all3` (γ 0.45)** | −0.00707 | −0.02026† | +0.00218 | +0.00077 | +0.00014 | +0.00024 | **−0.02401** | −0.02668 |
| `all3` (γ 0.30) | −0.00929 | −0.02026† | +0.00211 | +0.00040 | +0.00002 | −0.00000 | **−0.02701** | −0.02944 |

† `score_mse` clamped to its 0 floor; the delta is the whole a15 credit lost and is a censored
lower bound. `TOTAL` uses the measured h1-test `reach` delta; `TOTAL*` replaces it with rule
13's −0.0146 scaled emitted-profile toll. **165's finding is binding — h1-test reproduces board
`reach` for emission-law changes but NOT for profile changes, and every arm here is a profile
change** — so neither `reach` column is trustworthy alone and both totals are given. The verdict
is identical under either.

**Resolution.** Matched re-seeded twins on five arms:

| pair | Δ precision | Δ nmae | Δ reach | Δ recall |
|---|---:|---:|---:|---:|
| `rb_a15` / twin | −0.00032 | +0.00003 | **−0.02777** | −0.00083 |
| `rb_g045` / twin | +0.00018 | −0.00021 | +0.01076 | +0.00519 |
| `rb_flip` / twin | +0.00017 | −0.00002 | +0.00096 | −0.01698 |
| `rb_w211` / twin | −0.00048 | +0.00007 | −0.00372 | −0.01414 |
| `rb_flip_w211` / twin | +0.00000 | +0.00005 | +0.00016 | — |
| `rb_all3_g045` / twin | +0.00089 | −0.00002 | +0.00526 | +0.00055 |

`fid` and `nmae` are resolved to ~1/5 and ~1/20 of the composite's own delta. **`reach` is not**
— the base's own twin gap (0.0278 raw) exceeds every between-arm move — and neither is `jac`,
whose recall twin gap runs to 0.017. Both are reported and used for nothing on their own. `mse`
and `pds` carry no emission noise at all (170: bit-identical across seeds, because the quasi
emitter's exact-total repair pins the pseudobulk).

---

## 5. The composition result — the gains discount, the costs do not

Composite delta ÷ sum of its riders' deltas, per member:

| composite | `pds` | `mse` | `fid` | `nmae` | `reach` |
|---|---:|---:|---:|---:|---:|
| **`flip + w211`** | **0.701** | **0.955** | 1.159 | 2.226 | 2.375 |
| `g045 + flip` | 0.811 | 1.226 | 0.833 | 0.718 | 0.804 |
| `g045 + w211` | 0.636 | 1.192† | 1.111 | 0.765 | −0.655 |
| `all3` (γ 0.45) | 0.567 | 1.359† | 0.825 | 0.838 | −0.135 |
| `all3` (γ 0.30) | 0.519 | 1.001† | 0.620 | 0.289 | 0.001 |

† the composite's `mse` is clamped, so these three ratios are lower bounds on the true
super-additivity, not measurements.

**Read across the two channels that are both large and resolved:**

- **`pds` — the dominant cost — composes SUB-additively, 0.52 to 0.81, on all five composites.**
- **`mse` — the other large channel, and also a cost for every arm but one — composes
  ADDITIVELY to SUPER-additively, 0.955 to 1.226 on the two unclamped composites.**
- **The DE gains — `fid` and `nmae`, the two resolved positives — compose sub-additively,
  0.29 to 0.84**, on every composite except `flip + w211`, whose gains are tiny in absolute
  terms (+0.0022 and +0.0006) and whose ratios sit on a small denominator.

So the pre-registered hope — *gains add, `pds` costs saturate* — is **refuted in the strongest
possible way: the costs that saturate are the ones that were already smaller, the cost that
does not saturate (`mse`) is the one the riders turned out to carry, and the gains are what
saturate.** Composing makes the exchange rate worse, not better.

**The mechanism, and why it settles the family.** Two distinct mechanisms, both measured:

1. **`pds` sub-additivity comes from overlapping populations.** All three riders act on the
   low-`|cons|` entries the whitener rotates hardest, and each removes work from the next: the
   flip fraction on the real 300-target panel falls **18.0% of on-axis entries at γ 0.65 →
   14.9% at γ 0.45 → 12.5% at γ 0.30 with the K562-heavy mix**. (The weighting does not share
   this: at γ 0.45 the flip fraction goes 14.9% → 15.6% when the mix changes, so `w211`'s
   `pds` sub-additivity with the other two is a geometry effect, not a population overlap.)
2. **`mse` super-additivity comes from the ρ cliff being a threshold, not a slope.** Each rider
   adds amplitude; once `B` crosses `claim` the refund saturates, so the second rider's
   amplitude is charged in full while the first's was partly refunded. That is why the two-way
   composites read ≥ 1.0 and the three-way ones clamp out entirely.

---

## 6. Rule 13's gate, with the arithmetic

```
predicted = PLOLO d_pds/6 × 1.011  +  d_mse/6  −  0.0146/6
```

Every arm reshapes the profile the renderer reads, so the `mse` term is in (170 measures it
directly on the real panel), the emitted-profile `reach` toll is in, and **the binding
threshold is the amplitude-class +0.0089**.

| arm | `pds` term | `mse` term | `reach` toll | **predicted** | gate | verdict |
|---|---:|---:|---:|---:|---:|---|
| `g045` | −0.003466 | −0.006235 | −0.002433 | **−0.012134** | +0.0089 | NO-GO |
| `g030` | −0.008893 | −0.011566 | −0.002433 | **−0.022893** | +0.0089 | NO-GO |
| **`flip`** | −0.009127 | **+0.002094** | −0.002433 | **−0.009467** | +0.0089 | NO-GO |
| `w211` | +0.000126 | −0.010765 | −0.002433 | **−0.013073** | +0.0089 | NO-GO |
| `flip + w211` | −0.006306 | −0.008285 | −0.002433 | **−0.017024** | +0.0089 | NO-GO |
| `g045 + flip` | −0.010213 | −0.005077 | −0.002433 | **−0.017724** | +0.0089 | NO-GO |
| `g045 + w211` | −0.002126 | −0.020264 | −0.002433 | **−0.024822** | +0.0089 | NO-GO |
| **`all3` (γ 0.45)** | −0.007073 | −0.020264 | −0.002433 | **−0.029770** | +0.0089 | NO-GO |
| `all3` (γ 0.30) | −0.009285 | −0.020264 | −0.002433 | **−0.031982** | +0.0089 | NO-GO |

**Sensitivity — re-tuning `VCC_ALPHA` onto each arm's own knee** (the coordinator's operational
corollary). This is the single most generous convention available: it hands every arm the
amplitude that maximises its own `mse` and charges it nothing for the change. It is *not* the
verdict, because `nmae`, `fid` and `reach` all move with amplitude and are not re-read there,
and because the required alphas run from 0.90 to 1.02 against a15's board-tuned 0.982.

| arm | ρ (A/B/C) | α\* | `mse` shipped | `mse` at knee | predicted at knee | verdict |
|---|:---:|---:|---:|---:|---:|---|
| `g045` | 1.000/1.000/0.977 | 0.9637 | −0.006235 | +0.001021 | −0.004879 | NO-GO |
| `flip` | 0.991/0.940/0.881 | **1.0151** | +0.002094 | +0.002174 | −0.009387 | NO-GO |
| `w211` | 1.000/1.000/0.992 | 0.9541 | −0.010765 | −0.000731 | −0.003038 | NO-GO |
| `flip + w211` | 1.000/1.000/0.992 | 0.9538 | −0.008285 | +0.001746 | −0.006993 | NO-GO |
| `all3` (γ 0.45) | 1.000/1.000/1.000 | 0.9150 | −0.020264 | +0.003101 | −0.006406 | NO-GO |
| `all3` (γ 0.30) | 1.000/1.000/1.000 | 0.8965 | −0.020264 | +0.004341 | −0.007378 | NO-GO |

Under the most generous reading anywhere in this experiment the least-bad arm is `g045 + w211`
at **−0.0029**. The gate is **+0.0089**. Nothing is within a factor of three of it.

---

## 7. The instrument disagreement, stated plainly

`pds` on h1-test and `pds` on PLOLO disagree **in sign** for gamma, for the donor weighting and
for `flip + w211`:

| arm | h1-test rebuild, Δ raw `pds` | PLOLO, Δ raw `pds` (4-fold, challenge frame) | 097 §9h, h1-test |
|---|---:|---:|---:|
| `g045` | **+0.00098** | **−0.00930** | **−0.0036** |
| `g030` | −0.00462 | −0.02386 | −0.0090 |
| `flip` | −0.00840 | −0.02448 | — (173: −0.0042) |
| `w211` | **+0.00532** | **+0.00034** | — |
| `flip + w211` | **+0.00546** | **−0.01692** | — |
| `all3` (γ 0.45) | +0.00672 | −0.01897 | — |

**PLOLO is the one to believe, for four reasons.** (a) It is the only `pds` instrument in this
project that has retrodicted board `pds` — 6/6 signs over six historical geometry arms, fitted
scale 0.899, leave-one-arm-out RMS 0.001494 scaled, exact permutation p = 0.0028 (162). (b) It
scores 269–294 decoys in the challenge frame; `pds` is a retrieval statistic and h1-test offers
85. (c) **097's own h1-test measurement agrees with PLOLO on gamma's sign and disagrees with
this rebuild**, so two of three readings say lowering gamma costs `pds`, and the outlier is the
reconstruction whose per-row cosine to the shipped artifact is 0.672. (d) 170 found that
`pds`-vs-gamma is **not** monotone on the *delivered*, mass-fixed object (interior maximum at
0.50–0.55) where 097 found it monotone on the *intended* one — so small `pds` differences on a
rebuilt h1-test object are exactly where this surface is least reliable. Directionally the two
surfaces do **not** conflict on the accuracy channel: every rider raises transfer cosine and
sign accuracy on both. They conflict only on the retrieval geometry.

Per-fold PLOLO scatter is real and reported: `all3` (γ 0.45) reads −0.0677 / −0.0475 / −0.0446 /
−0.0081 scaled across `orion_hct116` / `orion_hek293t` / `k562_gwps` / `kolf2` (sd 0.025).
**Every fold has the same sign**, which is what the verdict rests on. `flip + w211` reads
−0.0631 / −0.0369 / −0.0421 / −0.0076, same story.

**A structural control worth naming.** In the fold that holds out `k562_gwps` there is no K562
to up-weight, so `w211` must read **exactly** `a23` — it does, to six digits (0.790712 /
0.786717 / 0.789935 on A/B/C), and `g045_w211` reads exactly `g045`. That control would have
failed loudly if the weighting had been mis-plumbed, and unlike a row permutation it is not
invariant by construction: the other three folds all move.

---

## 8. What actually killed each rider

**`w211` — whitening consumes the sign gain, and then the cliff charges for the amplitude.**
175's motivation was +0.0043 pooled sign accuracy for (2,1,1) over (1,1,1) on the h1-test
reference DE pool. Measured here on the same pool, the same code path:

| | (1,1,1) | (2,1,1) | Δ |
|---|---:|---:|---:|
| **before** the operator (175's row, reproduced exactly) | 0.5697 | 0.5740 | **+0.0043** |
| **after** the operator (this experiment) | 0.5417 | 0.5405 | **−0.0012** |

The rider's entire premise is a pre-whitening quantity the operator removes. Downstream it is
`pds`-neutral (+0.00034 raw, fold sd 0.003) and worth +0.0005 on `fid` — but it raises `A` by
13%, pushes ρ to 1.000 in every context, and loses **−0.0108** on `mse`. This revises 019 in the
right direction without contradicting it: 019's losing configuration was 4:1, and (2,1,1) is not
4:1 — it is simply worth nothing and costs amplitude.

**`flip` (sign repair) — the right shape, the wrong `pds` price.** It is the only
norm-preserving rider and the only one with a positive `mse`, and 170's criterion picked it out
in advance. But PLOLO reads its `pds` at **−0.0245 raw** where 173's h1-test read −0.0042 — 6×
— on all four folds. The flip population is the same size on every surface (17–19% of on-axis
entries: 19.1 / 19.2 / 17.2 / 18.0% across the PLOLO folds, 18.0% on the real panel, 14.5% of
*gated* entries on h1-test against 173's 13.9%), so this is not a mask that misbehaved: partial
un-whitening costs far more retrieval on 290 decoys than on 85. **173's conclusion is confirmed
in mechanism and its price was under-estimated by 6×.**

**`g045` / `g030` — 170's answer, re-measured independently and agreeing to 0.2%.** Nothing new
is claimed here; §3's cross-staging is the contribution, and it makes 170's gamma closure a
two-pipeline result.

---

## 9. Rule 12 — all six members, with the instrument or the mechanism that pins each

| member | priced? | how, and the caveat |
|---|---|---|
| `pds` | **yes** | PLOLO, 4 folds, challenge frame, the only board-retrodicting instrument (162). Same sign on every fold for every arm. |
| `mse` | **yes** | 170's per-context identity on the real 300-target panel; `A`, `B`, `claim` and ρ measured per arm per context, `X` anchored on a23's board point and carried by the h1-test ratio. Reproduces 170 to +0.00% on six known-answer rows. Three composites clamp; their deltas are censored lower bounds. |
| `fid` | **yes** | h1-test precision (`mean_p k/n_pred`, carries no `n_real`, rule 9), twin-averaged. Composite delta +0.0022 against a twin spread of 0.0005. |
| `nmae` | **yes** | h1-test, twin-averaged. Composite delta +0.0008 against a twin spread of 0.00001. |
| `jac` | **partly** | h1-test overlap ratio rebuilt on the board's own `n_real`/`n_pred`. Recall's twin spread (0.017) exceeds every between-arm move, so it is reported as ≈0 and used for nothing. |
| `reach` | **NO** | 165: h1-test reproduces board `reach` for emission-law changes but **not** for profile changes, and the measured base twin gap here (0.0278 raw) is larger than the composite's own delta (0.0198). Handled by giving both the measured column and rule 13's −0.0146 toll; the verdict is unchanged either way. |

---

## 10. NO-GO, and what it closes

**The rider family is closed.** Not "sub-threshold and composable" — closed, on the
board-validated instrument, under both amplitude conventions, with and without gamma.

1. **`WHITEN_GAMMA` is now closed three times**: 097's two-surface screen, 170's
   production-staged `mse` sweep, and this experiment's independent re-staging (which
   reproduces 170's `A`, `B`, `claim`, ρ and raw `mse` to ≤0.2%).
2. **Sign repair is closed as a shippable object**, and its 173 estimate is corrected: the
   `pds` cost is 6× larger than the h1-test surface said. It survives as the *measurement*
   that 15.7% of the reference's DE genes are anti-informative on sign after whitening (0.4573
   here, 0.4556 in 173) — a fact about the operator, not a lever. It is also the only rider
   that satisfies 170's norm-preserving criterion, which is worth keeping as a design filter.
3. **Donor weighting is closed with a mechanism**: whitening consumes the pre-whitening sign
   gain, and the mixture change costs amplitude at the ρ cliff. Any future donor-mixture idea
   must be argued *after* the operator and must not raise `A`.
4. **"Compose the sub-threshold riders" is closed as a strategy for this family**, and the
   general lesson is sharper than the specific one: **the `pds` cost saturates when riders
   overlap, but the `mse` cost does not, because the ρ cliff is a threshold.** Two riders that
   each look affordable alone can be jointly unaffordable. Price composites on `A`, `B` and ρ
   *first*; if the composite crosses the cliff, nothing downstream matters.
5. **The synthesis in `docs/CURRENT_STATE.md` stands and is now quantitative.** Production sits
   at the `pds`-optimal point of the `pds`-versus-accuracy curve *and* at the `mse` optimum of
   the amplitude curve, and the two pins are independent. Only new information can move the
   predictor.

### Not built, and deliberately

For the record, the two composites that were priced would be built by exactly these invocations
— **and neither must be built**:

```
# the gamma-free composite the coordinator asked for: 2 env differences against a23
VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.1 VCC_MASS_FIX=minnorm \
VCC_DONOR_WEIGHTS=2,1,1 VCC_SIGN_REPAIR=cons \
uv run python scripts/make_submission_pm_cd4.py --whitener-id union_pruned \
  --whitener data/donors/whitener_union_pruned.npz \
  --whitener-metadata data/donors/whitener_union_pruned.json \
  --min-cpm 10.0 --overdispersion 4.0 --tail-spread 1.0 --keep-top 8000 --keep-spread 0.25 \
  --output <path>                                          # predicted overall  -0.0170

# the original three-rider composite: 3 env differences against a23
#   ... + VCC_WHITEN_GAMMA=0.45                            # predicted overall  -0.0298
```

`effects178.py` verified that the first differs from the a15/a23 arm's written config in
**exactly two** fields (`donor_weights`, `sign_repair`) and the second in **exactly three**
(plus `whiten_gamma`), and in nothing else.

---

## 11. What I did NOT measure

- **`reach` for any arm, credibly.** 165 forbids the h1-test reading for profile changes and no
  other `reach` instrument exists. Both readings are given; both leave the verdict unchanged.
- **The 095 recovery path for gamma and the weights.** Only the rebuild can express them, and
  the rebuild's per-row cosine to the recovery is 0.672. All h1-test numbers are ratios inside
  the rebuild; absolute levels there are not production's.
- **PerturbMatch and the CD4 residual on the h1-test rebuild.** Both omitted (the 095 artifact
  has neither). Both ARE present in the PLOLO and `mse` routes, which are the ones the verdict
  rests on.
- **`lam`.** Held at 0, per the mandate and per 170's re-attribution (the ridge did 79% of
  a19's `mse` damage, gamma 21%).
- **175's fourth rider — per-donor rank-1 deflation (`VCC_DONOR_DEFLATE`)** — handed over at
  19:45/20:20, after this experiment's arms were fixed. Not priced, and here is the
  quantitative reason it cannot rescue the family. Its own addendum 3 measures that the
  whitener at γ 0.65 **already annihilates the long-gene axis** (9.76% of the un-whitened
  consensus's energy, **0.11%** of the shipped `H`), and that projecting it out of `H` leaves
  sign accuracy unchanged at 0.5373. So its only live channel is as a better sign *source* for
  the `flip` rider: it lifts `cons`'s own accuracy 0.5449 → 0.5570. Even granting that the
  entire +0.012 carries through and **doubles** sign repair's DE gain (+0.0018 → +0.0036
  overall) and leaves its `mse` credit intact, `flip`'s predicted total moves from −0.0095 to
  −0.0077. It would have to make sign repair's `pds` cost (−0.0091, four folds, same sign on
  all four) **disappear and reverse** to reach +0.0089. Removing a direction the whitener has
  already flattened to 0.11% of the energy cannot do that. Worth measuring for its own sake as
  an operator-refit question; not worth measuring as a rescue for this family.
- **`zero` (173's cheaper variant of sign repair).** 173 measured it at a third of the `pds`
  cost for two-thirds of the cosine gain; at PLOLO's exchange rate that is ≈−0.0027 on `pds`
  against ≈+0.0014 on `mse`, still negative before the toll, so it was not run.
- **`mse178.py`'s numbers.** Superseded by `chal178.py` (§3) and kept only for provenance. Do
  not cite them.

---

## Files

`effects178.py` (rule-10 effects stage + the ten arms' challenge-panel effect matrices),
`h1rebuild178.py` (h1-test rebuild, delivered-space `pds`/cos/sign, the flip anatomy, the
random-flip null), `plolo178.py` (PLOLO with per-source folds so the weighting is expressible),
`realize178.py` + `six178.py` (the DE members), `chal178.py` (the `mse` route, 170's identity —
**this is the one to read**), `mse178.py` (superseded), `gate178.py` (rule 13), `table178.py`
(the six-member assembly). Logs and JSON in `work/`; DE tables `de178.*.npz`.

Builder changes: `scripts/make_submission_pm_cd4.py` (`VCC_DONOR_WEIGHTS`, `VCC_SIGN_REPAIR`,
`sign_repair_to_consensus`), `scripts/test_make_submission_pm_cd4.py` (six new tests, 46 pass).
No file outside those two and this directory was modified.

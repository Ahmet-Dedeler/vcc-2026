# 208 — the free riders, composed on `rc_callset`, and the `jac` convention adjudicated

**Three results, in order of how much they change what we do next.**

1. **The `jac` convention is settled, and NEITHER convention retrodicts both board points —
   but the variable that separates them is identified.** Fixed-NP gets the **sign wrong** on
   the `rc_hybrid → rc_callset` pair (−0.00810 predicted against **+0.01012 realised**);
   scaled-NP gets the sign wrong on `a15 → a21`. Inverting the model for each arm's *implied
   board recall* at its **own measured challenge `n_pred`** reconciles them: for the three
   **profile-preserving** arms it reproduces the h1-test recall to within **2.0% / 0.0% /
   7.1%**, and for `a21` — the one arm that moved the delivered profile — it is off by
   **1.77×**. **Scaled-NP is correct for profile-preserving arms; fixed-NP fits `a21` only by
   cancelling two errors.** Under the surviving convention 200's tier-vs-flat-1.00 trade
   re-prices from **+0.00558 for dropping the tier to −0.00259** — i.e. **keep the tier**.

2. **Rider B (flat 1.00) is dead twice over, and the second kill is structural.** Under
   `rc_callset` the call set is **PRESCRIBED**, so the spread policy cannot move `n_pred` by a
   single gene — measured identical on **85/85** h1-test targets. The entire ±0.0066 `jac`
   channel that 200's verdict turned on is **exactly zero** against the current incumbent, and
   what is left is a pure `reach` loss: **−0.00096 overall**.

3. **Rider A (Hamilton rounding) survives, at +0.0029, and it is still not a slot.** `pds`
   +0.001195 (PLOLO, reproducing 196 bit-identically), `mse` **0.00000 with `claim` saturated
   on all eight configurations**, DE **+0.00170** of which **100% is `nmae`**. Rule 13's priced
   composite is **+0.001197 against a +0.0019 floor → VETO under all three readings of the
   toll.** But **its `reach` half has evaporated against the new incumbent** (+0.00607 raw
   against `rc_hybrid`, **+0.00006** against `rc_callset`) while its `nmae` half transferred
   **unchanged to five decimals** — so it is *less* of a slot and *exactly as much* of a rider
   than 196 thought, and the reason is mechanical.

4. **The composition is additive, measured.** A + B = +0.001940; measured A∘B = **+0.001949**.
   Difference **+9e−6**. 178's "gains compose sub-additively" does not bite here — the two
   riders touch orthogonal objects (integer totals vs cell-to-cell placement).

**And a loud caveat this experiment is required to state: the Hamilton rider is NOT
pseudobulk-preserving.** Its per-gene panel totals differ from `rc_callset`'s on **900 of 900**
challenge blocks. The exemption that pinned four of six members for `rc_hybrid` and
`rc_callset` **does not apply**, and every member below is priced without it.

---

## Phase 0 — known-answer rows (rule 10). Every gate PASS; nothing was re-derived that could be reproduced.

| # | gate | expected | measured | |
|---|---|---|---|---|
| 1 | 196's **challenge** tilt, `floor` → `ham` | −0.01234 → +0.00039 log2 | **−0.01253 → +0.00030** (`tilt208`, 40 blocks × 3 contexts; per-context −0.01189 / −0.01291 / −0.01259) | PASS |
| 1b | 196's challenge rounding energy per rule | 2.837 / 3.126 / 3.254e−6 (floor), 9.288 / 9.736 / 10.55e−7 (ham) | **identical to every printed digit** | PASS |
| 2 | 196's h1-test tilt, `floor` → `ham` | −0.01352 → +0.00019 | **−0.01352 → +0.00019** | PASS |
| 3 | 196's PLOLO `pds`, `floor` | 0.739016 | **0.7390163577** (\|d\| **0.0e+00**) | PASS |
| 4 | 196's PLOLO `pds`, `ham` | 0.742224 | **0.7422236676** (\|d\| **0.0e+00**) | PASS |
| 5 | 169's cached PLOLO chalA rows at `m = ∞` | 0.790712 / 0.790896 / 0.771621 / 0.615317 | worst \|d\| **2.8e−07** | PASS |
| 6 | 196's `mse` `m=1` row (`quasi`, floor, tier) | Ĉ/C_p 1.12577 · claim 0.246529 · raw 0.882283 | **1.125457 / 0.246483 / 0.882283** | PASS |
| 7 | **200's flat-1.00 row** (`quasi`, floor, flat) | claim **0.24653**, raw `mse` **0.88228** bit-equal to a15's | **0.246529 / 0.882283** | PASS |
| 8 | `rc_callset` h1-test DE table, through THIS code path | bit-identical to `de202.rc_callset.npz` | **`pred_p`, `pred_padj`, `pred_lfc`, `pred_z` bit-identical on 85 × 10,779, max \|diff\| 0**; totals digest 85/85 | PASS |
| 9 | `rc_callset` h1-test four members | 2634 / 0.223341 / 0.494362 / 0.072975 / 0.966702 | identical to 1e−12 | PASS |
| 10 | `rc_hybrid` h1-test four members (191's) | 0.21681 / 6043 / 0.06827 / 0.96670 | **0.21681 / 6043 / 0.06827 / 0.96670** | PASS |
| 11 | **`rc_hybrid` challenge `n_pred` and order, off the artifact** | A 4,651 / B 4,480 / C 4,354, pooled **4,480**; ρ +0.9761/+0.9713/+0.9755, mean **+0.9743** | **identical, every digit** | PASS |
| 12 | **`rc_callset` challenge `n_pred` and order, off the artifact** | 2,599 / 2,599 / 2,599; ρ mean **+0.9993** | **2,599 pooled; ρ +0.9994/+0.9996/+0.9990, mean +0.9993** | PASS |
| 13 | package sha256 vs what the board scored | `e0b0677f…` / `2a894aa5…` | **both PASS** on 13,202,460,120 / 12,648,807,800 bytes | PASS |
| 14 | `rc_callset` artifact totals ≡ a15's recorded totals | 900/900 | **900/900** (first block `A\|ABCD1` digest `bb6befc01caf7493`, the value 191 recorded) | PASS |

### Rule 11 — path and sha256, never a description

**14** whitening operators are on disk under near-identical names. The shipped one, recorded in
`meta['operator']['sha256']` of **both** board entries, is:

| path | sha256[16] | max eig | n |
|---|---|---:|---:|
| **`data/donors/whitener_union_pruned.npz`** | **`e2f65ecf69b7e0b0`** | **5.894570** | **17942** |
| `data/donors/whitener_union.npz` | `5608e807230dea04` | 6.256702 | 17942 |
| `data/donors/whitener_union_full.npz` | `965d3cef36ac6adc` | 6.980217 | 17942 |
| `experiments/095-clean-harness/out/whitener_union.npz` | `5cf3899db9272ff2` | 8.036353 | 17942 |
| `experiments/096-mse-floor/out-alpha/whitener_union.npz` | `5cf3899db9272ff2` | 8.036353 | 17942 |
| `experiments/095-clean-harness/out-h1test/whitener_union.npz` | `f80fe8d6ba59138b` | 8.080842 | 17745 |
| `experiments/095-clean-harness/out-essential/whitener_union.npz` | `1fbd44b0872b60ba` | 7.734004 | 17942 |
| `experiments/095-clean-harness/out-random/whitener_union.npz` | `05a17699c8433e42` | 8.361587 | 17942 |
| **`experiments/199-def-runbook/work/whitener_union_pruned_def199.npz`** | `68dd20dc5354f401` | **5.938910** | 17942 |
| `experiments/040-union-gamma/out/whitener_union.npz` | `f950b268f35114d4` | 6.899901 | 17530 |
| `experiments/099-donors/work/whit/whitener_hek293t.npz` | `4502bb4f619b708a` | 3.294022 | 17942 |
| `data/donors/whitener_full.npz` | `dc6bd2cea9cea2f5` | 4.458631 | 7312 |
| `140-.../whitener-k562.npz`, `…-non_k562.npz` | `a68eb0d5…`, `269c41d1…` | 3.428861, 2.701425 | 4992 |

> ⚠️ **A NEW near-collision has appeared since 200 wrote its table.** 199's
> `whitener_union_pruned_def199.npz` has the shipped operator's **name pattern**, its **shape**
> (17,942) *and* a max eigenvalue of **5.938910** against the shipped **5.894570** — 0.75%
> apart, where 200's trap pair agreed to five significant figures on the eigenvalue but
> differed in name. This one differs in neither meaningfully. **Only the sha separates them.**

Nothing in this experiment loaded an operator: the effect matrix
(`experiments/178-composed-riders/work/eff.a15.npy`, `effects.final_sha256`
**`22c09f22c749bda8`**) was identified by sha against **both** packages' own meta and matched.
Twenty-four frozen inputs are hashed in `work/inputs208.json`. One incidental confirmation
worth recording: `experiments/202-callset-rider/work/de202.rc_hybrid.npz` is **byte-identical**
to `experiments/191-rc-hybrid-build/work/de191.rc_hybrid_seq.npz` (both `e6af1b95762aaf3f`).

---

## Phase 1 — the `jac` convention, adjudicated against the board

### 1.1 The two conventions, written down

Both share the mean-field Jaccard of a board perturbation,

```
J(R, NP, NR) = R·NR / (NP + NR − R·NR)                                   (M)
```

`R` = **recall** (the share of the reference's DE genes we call), transferred from h1-test;
`NP` = our call-set size on the **challenge**; `NR` = the challenge's `n_real`. They differ
only in `NP`:

```
F    fixed-NP            NP_arm = NP_base                       (169's, 200's incumbent)
S    scaled-NP, h1       NP_arm = NP_base · (n_pred^h1 ratio)
S*   scaled-NP, DIRECT   NP_arm = the arm's OWN MEASURED challenge n_pred
```

`NR` is **calibrated on the base point** — solved so `M` reproduces the base arm's own board raw
`jac` — so each convention gets its best shot rather than a borrowed constant.

### 1.2 The two board pairs, in opposite directions

| pair | challenge `NP` | board raw `jac` | h1-test recall | profile |
|---|---|---|---|---|
| `a15` → `a21` | 4,420 → **1,141** | 0.030080 → 0.023997 (**−0.006083**) | 0.669 → 0.271 | **MOVED** (board `nmae` −0.064214 scaled) |
| `rc_hybrid` → `rc_callset` | 4,480 → **2,599** | 0.03031 → 0.04043 (**+0.010120**) | 0.66929 → 0.49429 | **pinned** (`pds`/`mse` identical to 16 s.f.) |

### 1.3 The retrodictions — neither convention survives both

| pair | convention | NP used | pred `jac` | pred Δ | realised Δ | realised/pred | sign |
|---|---|---:|---:|---:|---:|---:|---|
| rc_hybrid→rc_callset (NR\*=206.0) | **F fixed-NP** | 4,480 | 0.02221 | **−0.00810** | +0.01012 | −1.25 | **WRONG** |
| | S scaled-NP (h1) | 1,953 | 0.04950 | +0.01919 | +0.01012 | 0.53 | OK |
| | **S\* scaled-NP (direct)** | **2,599** | **0.03766** | **+0.00735** | +0.01012 | **1.38** | **OK** |
| a15→a21 (NR\*=201.7) | **F fixed-NP** | 4,420 | 0.01197 | −0.01811 | −0.00608 | 0.34 | OK |
| | S scaled-NP (h1) | 1,401 | 0.03531 | +0.00523 | −0.00608 | −1.16 | **WRONG** |
| | S\* scaled-NP (direct) | 1,141 | 0.04244 | +0.01236 | −0.00608 | −0.49 | **WRONG** |

The `NP` exponent that hits the realised move exactly (`NP_arm = NP_base·(NP₁/NP₀)^θ`, θ=0 is F,
θ=1 is S\*) is **θ = +1.136** for the rc pair and **θ = +0.539** for a21.

### 1.4 The reconciliation — invert `M` for each arm's IMPLIED BOARD recall

At each arm's **own measured** challenge `NP`, with `NR` = 206.0 held for all four:

| arm | board `NP` | board raw `jac` | implied board recall | h1-test recall | implied/h1 | profile |
|---|---:|---:|---:|---:|---:|---|
| a15 | 4,420 | 0.03008 | 0.65586 | 0.66900 | **0.980** | pinned |
| rc_hybrid | 4,480 | 0.03031 | 0.66929 | 0.66929 | **1.000** | pinned |
| rc_callset | 2,599 | 0.04043 | 0.52920 | 0.49429 | **1.071** | pinned |
| **a21** | **1,141** | **0.02400** | **0.15326** | **0.27100** | **0.566** | **MOVED** |

**That is the finding.** With the call-set size taken at its measured value, the board's recall
tracks h1-test's to within 7% on every arm that held the delivered profile, and misses by
**1.77×** on the one arm that did not. `a21`'s recall did **not** transfer — its profile change
collapsed the challenge's recall far further than h1-test's — and fixed-NP "works" there only
because a wrong `NP` and a wrong recall push in opposite directions and partially cancel.

**Fixed-NP is not a convention; it is a false statement about a measured quantity.** We counted
`rc_callset`'s challenge call set at **2,599** on the built artifact against the real control
cells, before submission (reproduced here, Phase 0 gate 12). Asking the model to use 4,480 for
it is not a modelling choice.

**Verdict: use S\* — the arm's own measured challenge `n_pred` — for profile-preserving arms,
which is the class 200's trade belongs to.** It retrodicts `rc_callset` at **1.38×** (0.93× if
`NR` is taken at 192's mean band of 330 rather than the base-calibrated 206). Do **not** use it
for an arm that moves the delivered profile without first re-checking that its recall transfers.

### 1.5 The rule-8 sensitivity, which fails, and why that is informative

A per-perturbation form (`J = mean_p R_p·NR_p/(NP_p + NR_p − R_p·NR_p)` with h1-test's own
`n_real` scaled by `f`) gets the sign **wrong for both conventions at every `f` tested**
(F −0.0147…−0.0396, S\* −0.0022…−0.0149, against a realised +0.0101). The cause is measurable:
h1-test's `n_real` has **median 89 but mean 765**, with individual perturbations up to **5,484**,
and the two heaviest rows alone contribute −0.059 and −0.071 to the mean. 192 §4.2 puts the
board's `n_real` at median 31–45 / mean 270–383. **Scaling h1-test's distribution by a scalar
preserves a tail the board does not have** — a textbook instance of AGENTS §11.1 rule 9. The
mean-field form with a single calibrated `NR` is, for once, the *less* wrong instrument.

### 1.6 Re-pricing 200's trade

200's sign convention (positive = what removing the tier would cost):

| | `reach` | `fid` | `jac` | `nmae` | SUM |
|---|---:|---:|---:|---:|---:|
| tier 8000 − flat 1.00, **F fixed-NP** | −0.000344 | +0.001387 | −0.006590 | −0.000034 | **−0.005581** → dropping the tier *gains* +0.00558 |
| tier 8000 − flat 1.00, **S\* scaled-NP** | −0.000344 | +0.001387 | +0.001576 | −0.000034 | **+0.002585** → dropping the tier *costs* −0.00259 |

**Under the surviving convention, keep the tier.** 200's "cheapest live arm this experiment
found" is negative. And §2.2 below shows the question no longer even arises against the current
incumbent.

---

## Phase 2 — the composed price, marginal to `rc_callset`

Four arms, factorial, one process, seeds and everything else `rc_callset`'s bit for bit:
`rule ∈ {floor, ham} × spread ∈ {tier, flat}`. `floor_tier` **is** `rc_callset` (Phase 0 gate 8).

### 2.1 `pds` and the transfer cosine — PLOLO, 4 folds, challenge frame

| rule | k562_gwps | orion_hct116 | orion_hek293t | kolf2 | MEAN `pds` | Δ | overall | macro cos | Δcos |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| floor | 0.791724 | 0.785598 | 0.766792 | 0.611951 | 0.739016 | +0.000000 | +0.000000 | 0.021762 | +0.000000 |
| **ham** | 0.790948 | 0.790860 | 0.771285 | 0.615802 | **0.742224** | **+0.003207** | **+0.001182** | 0.021768 | **+0.000005** |
| ∞ | 0.790712 | 0.790896 | 0.771621 | 0.615317 | 0.742137 | +0.003120 | +0.001150 | 0.021784 | +0.000021 |

> **The cosine does not move, and that is the mechanism.** Δcos is **+0.02% relative**, an order
> of magnitude inside 204's **±0.20%** resolution, and the paired target bootstrap splits by
> fold — k562_gwps **−0.000285** [−0.000405, −0.000165], orion_hek293t **+0.000302**
> [−0.000265, +0.000852], the other two straddling zero. **Hamilton buys `pds` without buying
> direction.** It is the same signature as whitening (§4.2: cosine −32%, `pds` +8.1%): what
> `floor` costs is the **shared, systematic** part of the rounding error — a tilt common to all
> 300 predictions, which is a hubness term — not the information content. This also explains why
> raw `mse` moves by only −0.00001: 197's identity reads `mse` off the cosine, and the cosine is
> where nothing happened.

The spread policy is absent from this table **by identity, asserted not assumed**: `emit_counts`
computes `totals` from `(n, L, profile)` before `spread` is read, and `emit208.py` confirmed the
per-gene panel-total digest identical on **85/85** h1-test blocks for `floor_flat` vs
`floor_tier`. So B's `pds` and `mse` are exactly 0.

### 2.2 The four DE members on h1-test — and the structural kill of rider B

Pool `native` (the reference's own BH set and our own BH call set):

| arm | `n_pred` | recall | `reach` | `fid` | `jac` | `nmae` |
|---|---:|---:|---:|---:|---:|---:|
| `floor_tier` = **`rc_callset`** | 2634 | 0.49429 | 0.22334 | 0.49436 | 0.07297 | 0.96670 |
| `ham_tier` (**A**) | 2629 | 0.49423 | **0.22340** | 0.49429 | 0.07302 | **0.96040** |
| `floor_flat` (**B**) | **2634** | **0.49429** | **0.21937** | 0.49437 | **0.07297** | 0.96670 |
| `ham_flat` (**A∘B**) | 2629 | 0.49423 | 0.21948 | 0.49429 | 0.07302 | 0.96040 |

Pool `wald89` (188's significance-selected pool, `|lfc|·√cpm` at the reference's BH depth,
median `n_real` 89 — the construction the mandate specified, because a top-N-by-magnitude pool
inverts `nmae`'s sign): `reach` 0.22384 / **0.22389** / 0.21934 / 0.21939; `nmae` 0.96965 /
**0.96505** / 0.96965 / 0.96505. Pool `matched` (our call set truncated per row to
`floor_tier`'s own `n_pred`) reproduces `native` to four decimals.

> ### The kill: under `rc_callset` the spread policy cannot move the call set at all
>
> `floor_tier` and `floor_flat` deliver **identical `n_pred` on 85 of 85 targets**; so do
> `ham_tier` and `ham_flat`. `jac` is therefore **bit-identical** (0.07297 both) and recall is
> **bit-identical** (0.49429 both).
>
> **Mechanism.** `rc_hybrid_block`'s `spread` acts only on **pass 1**, whose job with a
> prescribed `call_set` is to supply the `|z|` **multiset** — not the set. The delivered set is
> `prescribed ∪ unsuppressible`, and `z_cap = 1.90` puts every complement gene's p-value at
> `2·sf(1.90) = 0.0574 > α`. Which genes clear BH is fixed by the prescription; the spread can
> only reorder values *inside* it. **200's entire tier-vs-flat verdict lived in `jac`, and
> against the current incumbent that channel is exactly zero.** What remains is a pure `reach`
> loss of **−0.00397 raw h1-test**, i.e. **−0.00096 overall.**

### 2.3 The six-member price (overall board units, marginal to `rc_callset`)

`reach` × **1.290** (`rc_hybrid`'s own measured h1-test → board factor for this arm class);
`fid` as precision (`n_pred` 2,599 ≫ the challenge's `n_real` ≈ 222, so the yield cap is
saturated and `fid` **is** precision); `jac` through Phase 1's surviving S\* convention with
board `NP` **pinned at 2,599 for every arm** (§2.2), so it can move only through recall, `NR`
calibrated on `rc_callset`'s own board point at **221.8**; `nmae` h1-test raw / span.

| member | instrument | **A** Hamilton | **B** flat 1.00 | **A∘B** |
|---|---|---:|---:|---:|
| `pds` | PLOLO 4-fold × 1.011 | **+0.001195** | 0.000000 *(identity)* | **+0.001195** |
| `mse` | 096 `claim`, challenge | +0.000002 | +0.000000 | +0.000002 |
| `reach` | h1-test × 1.290 | +0.000014 | **−0.000964** | −0.000936 |
| `fid` | h1-test precision | −0.000041 | +0.000005 | −0.000041 |
| `jac` | S\*, `NP` pinned 2,599 | −0.000002 | +0.000000 | −0.000002 |
| `nmae` | h1-test raw / span | **+0.001731** | +0.000000 | **+0.001731** |
| **SUM** | | **+0.002899** | **−0.000959** | **+0.001949** |
| SUM at `reach` ×0.4277 | sensitivity | +0.002890 | −0.000315 | +0.002575 |
| SUM, pool `wald89` | | +0.002432 | −0.001087 | +0.001339 |
| SUM, pool `matched` | | +0.002954 | −0.000959 | +0.002004 |

**Composition, measured rather than inferred (178's rule).** A + B = **+0.001940**; the measured
A∘B is **+0.001949**. The two riders are **additive to 9e−6** — they touch orthogonal objects
(integer gene totals vs the placement of those counts across cells), which is why the usual
sub-additivity does not appear. Same conclusion on `wald89` (−6e−6) and `matched` (+9e−6).

### 2.4 `mse` — the binding check, and it clears with room

a15 had only **5.6%** of margin above the knee. All eight configurations are **SATURATED with
zero blocks under the cap**:

| emitter | rule | spread | Ĉ/C_p | rows<cap | `claim` | branch | raw `mse` | Δ overall |
|---|---|---|---:|---:|---:|---|---:|---:|
| `quasi` (a15) | floor | tier | 1.125 | 1 / 24 | 0.24648 | SATURATED | 0.88228 | +0.00000 |
| `quasi` | floor | **flat** | 4.048 | 0 | 0.24653 | SATURATED | **0.88228** | +0.00000 |
| `quasi` | ham | tier | 1.130 | 0 | 0.24653 | SATURATED | 0.88227 | +0.00000 |
| **`rc_callset`** | **floor** | **tier** | **2.940** | **0** | **0.24653** | **SATURATED** | 0.88228 | +0.00000 |
| `rc_callset` | ham | tier | 2.956 | 0 | 0.24653 | SATURATED | 0.88227 | +0.00000 |
| `rc_callset` | floor | flat | 5.246 | 0 | 0.24653 | SATURATED | 0.88228 | +0.00000 |
| `rc_callset` | ham | flat | 5.250 | 0 | 0.24653 | SATURATED | 0.88227 | +0.00000 |

**A by-product worth keeping.** The incumbent's margin is no longer 5.6% — `rc_callset` sits at
**Ĉ/C_p = 2.940** against a15's **1.125**. The `mse` cliff that vetoed the depth knob (196 §2c)
is now **2.6× further away**, because `rc_callset`'s `z_cap` re-assignment leaves more emitted
cell-to-cell scatter. Any future emission arm that was blocked by that cliff should be re-priced.

### 2.5 What changed between 196 and here, and why

Same rider, same code, two different incumbents:

| Δ, h1-test raw, `native` pool | `reach` | `fid` | `jac` | `nmae` |
|---|---:|---:|---:|---:|
| 196: `rc_hybrid` floor → ham | **+0.00607** | +0.00082 | +0.00010 | **−0.00630** |
| 208: `rc_callset` floor → ham | **+0.00006** | −0.00007 | +0.00005 | **−0.00630** |

**`nmae` transfers unchanged to five decimals; `reach` loses 99% of its value.** The rounding
tilt perturbed the *p-value ordering* under `rc_hybrid` (ρ(−log p, |H|) = +0.9743 on the
challenge); `rc_callset` already drives that to **+0.9993**, so there is nothing left for the
rounding rule to repair. `nmae` reads the **delivered profile**, which no call-set policy
touches, so it is untouched. **196's informal +0.00337 DE upside against `rc_hybrid` was
double-counting with the arm that has since shipped; its formal +0.00120 `pds` price is exactly
unchanged.**

---

## Phase 3 — the gate

Rule 13 prices **only** `pds` (PLOLO) and `mse` (096's claim). Its **+0.0019** floor **is** the
realised RMS of the four DE members it does not price, so putting them on the left-hand side and
comparing against the same floor double-counts. Both numbers are given.

```
priced   = PLOLO d_pds/6 x 1.011  +  096 d_mse/6
         = +0.001182 x 1.011      +  +0.000002        =  +0.001197
```

Arm class: `mode="quasi"`, the spread vector, `overdispersion`, `DISPERSION`, `L`, every seed and
the whole emission law are **bit-identical**; the operator, the donors and `VCC_ALPHA` are
untouched. It is **shape-only** → floor **+0.0019**. It is *not* profile-preserving (Phase 4),
and AGENTS §11.1 rule 13's amendment is internally inconsistent about which toll that attracts
(its body charges arms that "move the delivered profile"; audit 181's parenthetical restricts the
−0.0146 to emission-law arms and cites 169's +0.0006 for profile arms). **All three readings are
carried and all three VETO:**

| arm | priced | toll | composite | gate | verdict | six-member best estimate |
|---|---:|---:|---:|---:|---|---:|
| **A** Hamilton | +0.001197 | 0 (196's reading) | **+0.001197** | 0.0019 | **VETO** | **+0.002899** |
| | | −0.001322 (a23's *realised* profile-arm `reach` cost) | −0.000125 | 0.0019 | VETO | |
| | | −0.002433 (emission-law toll, charged anyway) | −0.001236 | 0.0019 | VETO | |
| **B** flat 1.00 | +0.000000 | 0 | +0.000000 | 0.0019 | **VETO** | **−0.000959** |
| **A∘B** | +0.001197 | 0 | **+0.001197** | 0.0019 | **VETO** | **+0.001949** |

**Honest board conversion of the best estimate.** The +0.0029 for A is **60% `nmae`**, on the
member rule 13 records as **3/6 signs — indicative, not priced**. Its transfer has one anchor,
and it is measured rather than assumed: the tilt Hamilton removes is **−0.01253 log2** on the
challenge against **−0.01352** on h1-test, i.e. **0.927×** (`tilt208.py`); scaling the h1-test
`nmae` gain by that ratio gives **+0.001604** (native) / **+0.001172** (wald89) instead of
+0.001731 / +0.001264, so the honest range for A is **+0.0028 to +0.0029** (native) and
**+0.0023 to +0.0024** (wald89). The `jac` transfer factor
the mandate supplies (1.10×, or the 1.38× this experiment measured in §1.3) is immaterial here:
`jac` moves by **−2e−6** because the call set is pinned.

**Verdict: A is a rider, exactly as 196 said, and B is dead.** No slot.

---

## Phase 4 — falsifiers

### F1 — Hamilton is **NOT** pseudobulk-preserving. Say it loudly.

Streamed over all **900** blocks of `submissions/202-rc_callset.h5ad` (12,648,807,800 bytes,
sha `2a894aa5…`), per-gene panel totals recomputed from the artifact's own CSR:

```
GATE A   artifact totals ≡ a15's recorded totals                 900/900   PASS
F1       Hamilton totals ≡ rc_callset's                            0/900   <-- BY DESIGN
           max |ΔT| per block          min 1     median 1     max 1
           genes moved per block       min 6,798 median 6,955 max 7,322   (of 18,533)
           …inside the scorer's gate   min 4,384 median 4,704 max 4,819
           direction                   UP on 6,798–7,322, DOWN on 0
           panel total S               8,000,726 -> 8,007,622  (+0.0876%)
```

**Every move is +1 and every move is upward**, because `floor` was short by `E[frac] = 0.5` on
every gene and largest-remainder hands the missing ~6,955 counts back. So `pds`, `mse`, `nmae`
are **not** pinned for A or A∘B and none of them was priced by exemption: `pds` came from PLOLO
(§2.1), `mse` from the claim route on the real challenge panel (§2.4), `nmae` from h1-test
(§2.2). The **flat-1.00** rider *is* in the preserved class, and that was asserted per emitted
block (`emit208.py` digest 85/85; `mse208.py`'s `assert B.sum(0) == totals` on all 192 challenge
emissions), never claimed.

### F2 — the known-answer emission

`floor` at the shipped policy through this code path reproduces `rc_callset`'s h1-test DE table
**bit-identically** — `pred_p`, `pred_padj`, `pred_lfc`, `pred_z`, max \|diff\| **0** on all
85 × 10,779 entries, totals digest 85/85, and the four members identical to 1e−12. **This code
path is the one the board scored.**

### F3 — the two format caps

| arm | nnz / block | ratio | package nnz | % of 4.75e9 | max counts/cell | % of 1e6 |
|---|---:|---:|---:|---:|---:|---:|
| **`rc_callset` as built (artifact)** | — | 1.0000 | **1,578,026,567** | **33.22** | **22,323** | **2.23** |
| `rc_callset` floor/tier (re-emitted) | 1,735,422 | 1.0000 | 1,578,026,567 | 33.22 | 21,320 | 2.13 |
| **`ham`/tier (rider A)** | 1,736,078 | **1.0004** | 1,578,623,298 | **33.23** | 21,338 | 2.13 |
| `floor`/flat (rider B) | 1,423,799 | 0.8204 | 1,294,666,525 | 27.26 | 21,320 | 2.13 |
| `ham`/flat (A∘B) | 1,425,289 | 0.8213 | 1,296,021,350 | 27.28 | 21,338 | 2.13 |

Hamilton costs **+0.04% nnz** (+0.01 percentage points of the cap) and **+18 counts** on the
fullest cell. Both caps cleared by **3.0×** and **47×**. 196's finding that the 1e6 counts/cell
cap binds before the nnz cap is untouched — neither rider goes near either.

---

## Verdict

1. **The `jac` convention is adjudicated, with a mechanism.** Use the arm's **own measured
   challenge `n_pred`** (S\*) for arms that hold the delivered profile — it retrodicts
   `rc_callset`'s +0.01012 at 1.38× (0.93× at 192's `NR` band) where fixed-NP gets the **sign
   wrong**. Fixed-NP survives `a21` only because that arm's *recall* also failed to transfer,
   by 1.77×, and the two errors cancel. **Neither convention retrodicts both points, and the
   separating variable is whether the arm moved the profile.** Record this: it is the first
   time this project has had two board `jac` points in opposite directions to test with.

2. **200's tier-vs-flat-1.00 trade is closed, twice.** Under the surviving convention the tier
   is worth **+0.00259**, so keep it. And against the current incumbent the question is moot:
   `rc_callset` prescribes the call set, `n_pred` is identical on 85/85 targets, `jac` and recall
   are bit-identical, and flat 1.00 is a pure `reach` loss of **−0.00096**. **Do not ship it.
   192's handoff from 200 can be closed.**

3. **196's Hamilton rider is confirmed as a rider and re-confirmed as not a slot**, against a
   different incumbent, with all six members priced and no exemption claimed. **+0.0029** best
   estimate, **+0.001197** on rule 13's priced channels against a **+0.0019** floor.
   `claim` stays **saturated** with 0/24 blocks under the cap, nnz +0.04%, both format caps
   clear.

4. **Two mechanism results that outlast the arms.** (a) Hamilton's `pds` gain is **not** a
   direction gain — Δcos is **+0.02% relative**, inside the noise floor, so the +0.0032 raw `pds`
   is entirely the removal of a **shared systematic tilt**, the same hubness channel whitening
   works in. (b) A rider's value is **not** a property of the rider: Hamilton's `reach` half was
   worth +0.00607 raw against `rc_hybrid` and is worth **+0.00006** against `rc_callset`, because
   `rc_callset` already drove the delivered order to ρ **+0.9993**. **Riders priced against a
   superseded incumbent must be re-priced, not carried forward.**

5. **One by-product for the next emission arm.** The incumbent's `mse` margin is no longer 5.6%:
   `rc_callset` sits at Ĉ/C_p **2.940** against a15's **1.125**, so the saturation cliff is
   **2.6×** further away than when 196 measured it. Arms that the cliff vetoed — 196 §2g's
   `m=2, φ=8`, for one — deserve a re-read.

### Handoff — the exact change, and a preregistered rider prediction

Not a slot. **Attach to whichever arm next earns one.** One line, in one place:

```python
# src/vcc2026/donor_transfer.py:380   (mode="quasi" -- the ONE live call site of five;
#                                      line 277 is resample_mix's, 498/540/567 are off-path)
- totals = np.floor(n * L * profile_p).astype(np.int64)
+ totals = hamilton(n * L * profile_p)      # largest remainder; experiments/196/ham196.py
```

`rc_hybrid_block` (`rc_emit.py:875`) calls `emit_counts(..., mode="quasi")` for its pass-1
block, so this single line moves a15, `rc_hybrid` and `rc_callset` together and preserves the
digest identity **between** them — which is what pins `pds` and `mse` across the family.

Build unchanged otherwise (`experiments/202-callset-rider/build202.sh rc_callset`), then assert,
before spending anything:

| # | assertion | expected |
|---|---|---|
| 1 | challenge `n_pred` (`probe191.py`, 40 × 3) | **2,599 ± 1%**; the call set is prescribed, the rounding rule must not move it |
| 2 | package nnz | **1.5786e9 ± 0.5%**, i.e. **33.2%** of the cap |
| 3 | max counts/cell | **< 25,000** (2.5% of the 1e6 cap) |
| 4 | per-gene panel-total digest vs `rc_callset` | **0/900 identical** — and `max |ΔT| = 1`, all moves UP, on 900/900 |
| 5 | rule-10 config diff vs `rc_callset` | **exactly one key** |
| 6 | `meta['operator']['sha256']` | **`e2f65ecf69b7e0b0`** |
| 7 | `effects.whitening.after_all_median_norm` | **4.229295604454299**, drift < 1e−12 |

**Preregistered per-member board prediction, marginal to `rc_callset` (0.1918638367):**

| member | prediction | basis |
|---|---:|---|
| `pds` | **+0.00120** overall (scaled +0.00717, raw ≈ +0.00324) | PLOLO 4-fold ×1.011, σ 0.000549 |
| `mse` | **0.00000**, `claim` stays saturated; raw `mse` within **1e−4** of 0.87367 | 096 claim route, 0/24 blocks under cap |
| `nmae` | **+0.0016 ± 0.0003** overall (raw ≈ −0.0058) | h1-test −0.00630 raw × the measured 0.927 tilt ratio |
| `reach` | **0.0000 ± 0.0010** | h1-test +0.00006 raw; the order is already at ρ +0.9993 |
| `fid` | **−0.00004 ± 0.0003** | h1-test precision |
| `jac` | **0.0000**; challenge `n_pred` unchanged at 2,599 | the call set is prescribed |
| **overall** | **+0.0027 ± 0.0013** | |

**Falsified if:** board `pds` or `mse` moves by more than 1e−3 in the wrong direction (the
delivered totals moved, so unlike the last two arms they are *not* pinned and this is a real
exposure); or challenge `n_pred` leaves 2,599 ± 1%; or board `jac` moves by more than 0.002
scaled.

---

## Files and reproduction

```
uv run python inputs208.py    # rule 11: 14 operators + 24 frozen inputs   -> work/inputs208.json
uv run python check208.py     # 196's reconstruction ≡ 185's a15 totals
uv run python tilt208.py      # Phase 0 gate 1, the challenge-side tilt    -> work/tilt208.json
uv run python jac208.py       # Phase 1, the convention adjudication       -> work/jac208.json
uv run python emit208.py      # 4 arms x 85 targets on h1-test, 716 s      -> work/de208.*.npz
uv run python pds208.py       # PLOLO pds + macro cosine, 2 KA gates       -> work/pds208.json
MAXT=8 uv run python mse208.py# 8 configs, 096 claim route, 2 KA gates     -> work/mse208.json
uv run python six208.py       # 4 DE members, 3 pools, board conversion    -> work/six208.json
uv run python price208.py     # the composite and rule 13's gate           -> work/price208.json
uv run python probe191.py rc_hybrid=submissions/191-rc_hybrid.h5ad \
                          a15_local=submissions/202-rc_callset.h5ad --perts=40
uv run python falsify208.py   # 900-block digest + the format caps         -> work/falsify208.json
```

`probe191.py` is a **verbatim** copy of 202's copy of 191's tool (sha
`81535b8e4efebcc75a653750d7729d0f2d411668dd1788427044a3733d66f266`, identical), so its output
lands in this directory; its positional slots print as `rc_hybrid` / `a15_local` and here those
hold `rc_hybrid` / `rc_callset`. Its own "falsifier 3" FAIL line compares the two packages'
`n_pred` to ±2% and is meaningless when the two packages are deliberately different call sets.

**Resources.** Peak RSS, `/usr/bin/time -l`: `mse208` **6.22 GB** (the heaviest instrumented
step), `pds208` 1.53 GB, `falsify208` 0.89 GB, `tilt208` 0.44 GB. `emit208` was launched without
instrumentation; its footprint is dominated by 185's cached control block (10,779 × 18,400
float64 = 1.59 GB) plus four 85 × 10,779 float64 tables per field, and it ran to completion
alongside `mse208` with no memory pressure. **One overshoot, recorded rather than hidden:** the
verbatim `probe191.py` peaked at
**13.1 GiB** for 43 seconds, above the 12 GB budget — it materialises each context's 18,400 ×
18,533 control block densely, which is 191's design and not something this experiment changed.
It ran once, serially, with nothing else in flight.

No `vcc`, no `cell-eval2`, no submission build, nothing written to `submissions/`, and no file
outside `experiments/208-free-riders/` modified.

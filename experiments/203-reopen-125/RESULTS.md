# 203 — Experiment 125 reopened, adjudicated, and closed for the right reason

**Verdict: 198's diagnosis is CONFIRMED in full, and the reopening is a NEGATIVE.**

The rejection of 125 *was* a pure scale artefact — measured, not argued, and now measured twice
more. But the estimator underneath it does not pay. On the only `pds` instrument validated
against a board point it delivers **−0.002372 raw `pds`** (95% CI [−0.011331, +0.007291]) and
**−0.00147 transfer cosine** (95% CI [−0.002163, −0.000752]) against production. Its
composed rule-13 gate is **−0.001393** against a **+0.0089** floor, and its measured `nmae`
adds another **−0.002468**. The direction-only variant 125 itself proposed is a pure ordering
arm whose key loses to the shipped `|H|` key by **−0.00281** at 20% depth and **−0.01429** at
the reference's BH depth.

**"The rescale control explains everything and the estimator adds nothing" is the result.**
The family is now closed on the endpoints the competition actually reads.

There is one live hand-off, in §9.

---

## 0. What was run, and every known-answer row

| script | what it does | known-answer rows |
|---|---|---|
| `phase0.py` | re-derives 125's two macro deltas from 122's frozen truth with 125's own metric code | 3, all pass |
| `phase0b.py` | scale-free direction instruments in 125's own frame | 3, all pass |
| `prep203.py` | rebuilds the per-source production effect matrices for all four PLOLO folds | 4, all pass (≤7.4e−8) |
| `test_est203.py` | asserts the re-implemented estimator == 125's `run.py` | **231, all pass at tol 0.0** |
| `score203.py` | the a15 chain, four folds, eleven arms | 1 (a15 vs `162/work/arms.json`), passes |
| `boot203.py` | paired per-target bootstrap, 2,000 draws | 2, pass (2.2e−16; a15 self-interval exactly [0,0]) |
| `phase23.py` | `nmae` identity, ordering keys, mechanism | 3, all pass |
| `nmae203.py` | `nmae` at each arm's own re-solved amplitude | 2, both pass exactly |
| `reach203.py` | ordering proxies for `reach`/`fid`/`jac` | — (labelled non-priceable) |

No `vcc`, no `cell-eval2`, no submission build, no count emission, no write outside
`experiments/203-reopen-125/`. All donor streaming ran at `max_rss_gib=7.0`.

---

## 1. Phase 0 — the diagnosis, verified from first principles

### 1.1 The two macro deltas reproduce exactly

From `125/work/source-only-evaluation-v1/per-target-{credit,precision}-delta.npy`:

| | recomputed | recorded | error |
|---|---:|---:|---:|
| macro signed credit | −0.015049899010796 | −0.015049899010796 | **0.0e+00** |
| macro top-20% precision | +0.012174572820954 | +0.012174572820954 | **0.0e+00** |
| all three family credits | — | — | **0.0e+00** |
| all three family precisions | — | — | **0.0e+00** |

Stronger: re-running 125's `score_predictions` (`run.py:752-793`) against 122's frozen truth
(`122/work/source-only-phase1-v1/family-effects.npy`) reproduces both deltas to 8 decimals from
the saved float32 prediction arrays.

| | 122 (`equal_mean_agreement`) | 125 (posterior) | delta |
|---|---:|---:|---:|
| signed credit | **+0.01841455** | +0.00336465 | −0.01504990 |
| top-20% precision | **0.53740893** | 0.54958351 | +0.01217457 |

### 1.2 An independent third confirmation nobody has recorded before

The saved 122 baseline arrays are **post-transport**. Scoring them reproduces 122's own
`RESULTS.md` top-20% precisions **to nine decimal places in every family**:

| held family | this recomputation | `122/RESULTS.md` |
|---|---:|---:|
| Replogle K562 | 0.535544645 | 0.535544645 |
| Orion | 0.542422605 | 0.542422605 |
| CD4 | 0.534259555 | 0.534259555 |
| macro | 0.537408934918 | 0.537408935 |

…while its signed credit is **0.6145× 122's own +0.029965312**. The `abundance_exponential`
transport — a magnitude transformation applied identically to both arms — moved "signed credit"
by **−38.5%** and the precision endpoint by **exactly zero**. That is the whole homogeneity
story, visible inside 122's own published numbers.

### 1.3 Homogeneity — measured on the real truth and the real confidence

198 verified this against a random stand-in truth. Redone with 122's actual held effects and
125's actual confidence array:

| k | `credit(k·p)/credit(p)` | `precision20(k·p) − precision20(p)` |
|---:|---:|---:|
| 0.25 | 0.250000000000 | +0.000000000000e+00 |
| 0.50 | 0.500000000000 | +0.000000000000e+00 |
| 0.75 | 0.750000000000 | +0.000000000000e+00 |
| 2.00 | 2.000000000000 | +0.000000000000e+00 |
| 6.78 | 6.780000000000 | +0.000000000000e+00 |

max `|ratio − k|` = **8.9e−16**; max `|Δprecision|` = **0.0**.
**"Signed credit" is exactly degree-1 homogeneous and the precision endpoint is exactly
scale-invariant. The mandate is not void.**

### 1.4 The shrink, and how much of 125 is *not* a rescale

| statistic | value |
|---|---:|
| L1 ratio 125/122 | 0.178597 (**5.599×** shrink) |
| L2 ratio | 0.207457 (**4.820×**) |
| least-squares scalar | 0.147491 (**6.780×**) |
| per-target LS scalar, median | 0.149318, IQR [0.120870, 0.183100] |
| global cosine(125, 122) | 0.710950 |
| per-target cosine, median | 0.748658, IQR [0.702934, 0.785978] |
| coordinates whose sign flips | **35,702 / 231,888 = 15.40%** |

So 125 is *not* only a rescale — a seventh of its coordinates flip sign. That residual is what
§2 and §4 price.

### 1.5 The naming collision is real

`142/develop.py:1377-1385`:

```python
def signed_credit(prediction, truth):
    numerator = np.sum(np.sign(truth) * prediction, axis=1)
    denominator = np.sum(np.abs(prediction), axis=1)     # the PREDICTION's L1
```

against `125/run.py:774-776`, which divides by `np.sum(np.abs(truth))`. One name, two different
homogeneity degrees, and 125's NO-GO rested entirely on the non-invariant one. Confirmed.

---

## 2. The pure-rescale control — the decisive experiment

**In 125's own frame.** Take 122's prediction, multiply by one scalar, re-score:

| scalar | source | Δ signed credit | as % of the recorded −0.015050 | Δ precision |
|---:|---|---:|---:|---:|
| 0.147491 | least-squares | **−0.01569856** | **104.31%** | **+0.000e+00** |
| 0.178597 | L1 ratio | **−0.01512577** | **100.50%** | **+0.000e+00** |
| 0.207457 | L2 ratio | −0.01459432 | 96.97% | +0.000e+00 |
| 0.149318 | per-target median | −0.01566493 | 104.09% | +0.000e+00 |

Every scalar the shrink can reasonably be summarised by reproduces the recorded loss to within
±4%, and the scalar that reproduces it **exactly** is 0.182717 — inside the L1/L2 band.
**The recorded −0.01505 is the shrink, and nothing else.**

**In production units.** `shrink6.78` = a15's own effect divided by 6.7801, pushed through the
full production chain and scored on the board-validated instrument:

| endpoint | Δ vs a15 | 95% CI |
|---|---:|---|
| `pds` (PLOLO, challenge frame, 4-fold) | **+0.001073** | [−0.001199, +0.003441] |
| transfer cosine | **+0.000049** | [−0.000075, +0.000175] |
| `nmae` (188 identity, wald89 pool) | +0.01087 → **+0.00000 once alpha is re-solved** | — |
| `reach*` / `fid*` / `jac*` ordering proxies | −0.00009 / −0.00060 / −0.00005 | — |

A 6.78× amplitude change — the transformation that produced 125's −0.01505 and its NO-GO — is
**indistinguishable from zero on every endpoint the competition reads**, and its one real cost
(`nmae`) is undone **exactly** by re-solving alpha (`shrink6.78@alpha*`: Δ`nmae` = +0.00000,
delivered norm ratio 1.0000). 198's claim that `VCC_ALPHA` corrects the shrink for free is
now a measurement, not an inference.

---

## 3. A structural fact that matters more than the scale artefact

125 was fitted and scored on a 300 × 300 block: panel targets × **panel genes as responses**.

Production sets every one of those 300 response coordinates to **exact zero** —
`062/run.py:287` (`effects[:, panel_positions] = 0.0`) and `162/run_arms.py:58`. Experiment 127
measured what they are worth: **2.283% of reference DE genes and 2.284% of absolute reference
LFC mass**.

So the frozen 125 candidate cannot be scored on `pds` at all: it lives entirely inside
production's null space. Everything below therefore re-fits the *same estimator*, verified
bit-identical (§0, 231 rows at tol 0.0), on the axis production actually ships.

---

## 4. Direction, with the scale divided out

### 4.1 In 125's own frame

Four scale-free instruments, 125's frozen arrays, 122's frozen truth:

| arm | cosine | pds* | `q_w` | 142's signed credit |
|---|---:|---:|---:|---:|
| 122 `equal_mean_agreement` | 0.03811 | 0.61331 | 0.50086 | 0.031559 |
| 125 hierarchical posterior | 0.03558 | 0.60647 | 0.49911 | 0.033311 |
| **delta** | **−0.00253** | **−0.00684** | **−0.00174** | **+0.001752** |
| 122 × 0.147491 (control) | +0.0e+00 | +0.0e+00 | +0.0e+00 | +0.0e+00 |

Known answers: the rescale control moves all four by exactly 0.0; `credit125` scales by exactly
0.147491000000; a sign-randomised prediction reads `q_w` = 0.48758.

All three held families lose cosine and pds*; two of three lose `q_w`. The only scale-free
endpoint that improves is **142's** signed credit — which is the prediction's *own* magnitude
weighting, i.e. a statement that 125's biggest predictions are more often correctly signed. That
is a **confidence-ordering** property, not a direction property, and §4.3 confirms it directly.

### 4.2 In the production frame — the number that decides it

PLOLO, challenge frame (mean of contexts A/B/C, delivered), unweighted mean over all four folds
— 162's pre-committed instrument, which retrodicts 12 board points and scores `pds` 6/6 on sign.

| arm | shrink | `pds` | Δ`pds` | 95% CI | P(Δ>0) | cosine | Δcos | Δ`q_w` |
|---|---:|---:|---:|---|---:|---:|---:|---:|
| **a15 (production)** | 1.00 | 0.740297 | — | — | — | 0.02184 | — | — |
| eq122 (125's comparator) | 1.03 | 0.733404 | **−0.006893** | [−0.012251, −0.001358] | 0.007 | 0.02087 | −0.00097 | −0.00078 |
| post125 `eb_broad_prior` **(125's pick)** | 4.70 | 0.737925 | **−0.002372** | [−0.011331, +0.007291] | 0.309 | 0.02037 | **−0.00147** | −0.00049 |
| post125 `eb_conservative` | 6.60 | 0.742771 | +0.002474 | [−0.004369, +0.009006] | 0.765 | 0.02220 | +0.00036 | +0.00068 |
| post125 `eb_balanced` | 6.81 | 0.733662 | −0.006635 | — | — | 0.02058 | −0.00126 | −0.00017 |
| post125 `eb_low_floor` | 7.40 | 0.717925 | −0.022372 | — | — | 0.01777 | −0.00407 | −0.00201 |
| shrink6.78 (control) | 6.78 | 0.741371 | +0.001073 | [−0.001199, +0.003441] | 0.818 | 0.02189 | +0.00005 | −0.00100 |

Two readings, and they point the same way.

1. **125's selected setting loses.** `eb_broad_prior` — chosen in 14/15 of 125's outer folds —
   is −0.002372 `pds` and **significantly negative on cosine** (CI excludes zero). Its per-fold
   `pds` deltas are −0.00068, −0.01930, −0.00099, +0.01148: **1/4 folds positive.**
2. **125's own comparator is worse than production.** `eq122` sits −0.006893 below a15 with
   P(Δ>0) = 0.007. 125's entire experiment compared two arms that both lose to the shipped
   estimator. Relative to its own comparator the posterior *does* gain +0.004521 `pds` — but
   with 2/4 folds positive and a per-fold sd of 0.011, that is not a finding either.

**The amplitude re-solve, and what it cost.** The mandate required reporting the arm at its own
re-solved alpha. Matching a15's delivered median row norm takes `eb_broad_prior` to α = 4.61 and
**costs −0.004645 `pds`** (0.737925 → 0.733280); `eb_conservative` to α = 6.48, costing
−0.002549. The re-solve is a cost, not a correction: the renderer is nonlinear, so a smaller
effect stays nearer the linear regime and delivers a profile closer to the intended one — which
is exactly why `shrink6.78` *gains* +0.001073. **On `pds`, scoring the shrunk prediction at
a15's alpha was never the confound the mandate feared**; both alphas are reported and the
un-shrunk one is worse.

### 4.3 The direction-only variant 125 proposed (Phase 2)

In production terms this is a **pure ordering arm**: it changes only the key the DE call is
ranked by, so `pds`, `mse`, `nmae` and the cosine are bit-identical to the incumbent by
construction (§8.2 rule 2, board-verified twice). So the question is only whether the key is
better. Measured on the *identical* delivered a15 profile against the PLOLO folds' real truth:

| ranking key | top-20% precision | vs `|H|` | top-89 precision | vs `|H|` |
|---|---:|---:|---:|---:|
| **`|H|` (production, 185/191)** | **0.51696** | — | **0.58995** | — |
| 125 conservative confidence | 0.51495 | −0.00201 | 0.58391 | −0.00604 |
| 125 posterior confidence | 0.51415 | **−0.00281** | 0.57566 | **−0.01429** |
| 122 confidence | 0.51180 | −0.00516 | 0.57179 | −0.01816 |
| random | 0.48734 | −0.02962 | 0.48620 | −0.10374 |

125's confidence **does** beat 122's confidence (+0.00235 at 20% depth, +0.00387 at 89) — which
reproduces 125's +0.0122 in a completely different frame and is the honest vindication of its
result. **But both lose to the key production already ships.** 125's ordering gain was a gain
over a key we do not use.

Known answer: the `|H|` key's precision after a 6.78× shrink is 0.51699271 against 0.51695824
unshrunk — Δ = +3.45e−05, i.e. a positive rescale cannot reorder `|H|` (the residual is the
renderer's nonlinearity). The channels separate cleanly.

**Which channel this was measured in:** `q_w` and cosine, and separately the ordering key. They
are not conflated. §8.2's warning applies exactly — a confidence layer is not worth pursuing for
`reach`, and it did not raise cosine or `q_w`.

### 4.4 The 2×2 decomposition of 125's original +0.0122

Same instrument, 125's own frame, all four combinations of effect and confidence:

| effect | confidence | signed credit | top-20% precision | Δ precision |
|---|---|---:|---:|---:|
| 122 | 122 | +0.018415 | 0.537409 | — |
| 122 | **125** | +0.018415 | 0.547936 | **+0.010527 (86.5%)** |
| 125 | 122 | +0.003365 | 0.537652 | +0.000243 (2.0%) |
| 125 | 125 | +0.003365 | 0.549584 | +0.012175 |

**86.5% of 125's headline gain is the confidence array; 2.0% is the effect's signs.** And the
122-effect/125-confidence row costs **exactly zero** signed credit — the direction-only variant
was always free on 125's own rejection endpoint. It just isn't worth anything on the shipped
key (§4.3).

---

## 5. The mechanism — why the hierarchical layer does not help

The four frozen posterior settings span a 3.5× range in `pds`. What orders them?

| setting | median cos(posterior, plain equal mean) | Δ`pds` (4-fold) |
|---|---:|---:|
| `eb_low_floor` | 0.55027 | −0.022372 |
| `eb_balanced` | 0.67371 | −0.006635 |
| `eb_broad_prior` | 0.67988 | −0.002372 |
| `eb_conservative` | 0.78734 | +0.002474 |

**Pearson r = +0.9593** (n = 4). The ranking is monotone: the more of the plain average a
setting preserves, the better it scores. The variance floor is the active knob, and it works by
flattening the per-gene precision weights back toward uniform. **The hierarchical structure is
not adding information; it is subtracting it, and the settings that look best are the ones that
subtract least.** This is exactly §4.2's finding 3 — whitening "resurrects nothing else", and
EB shrinkage keeps its deficit — reproduced on a new estimator four experiments later.

---

## 6. The six-member price (rule 12) and the gate (rule 13)

**Arm class.** This arm changes the effect matrix, so it is **amplitude-changing and
profile-changing**: the +0.0089 floor applies, `mse` is **not** exempt (a19 lost −0.01352
exactly here), and rule 13's −0.0146 `reach` toll does **not** apply (169 falsified it for
profile arms at +0.0006; the `rc_hybrid` board point confirmed the toll belongs to
profile-*changing emission-law* arms).

| member | instrument | `eb_broad_prior` | `eb_conservative` | status |
|---|---|---:|---:|---|
| `pds` | PLOLO, challenge frame, 4-fold, ×1.011/6 | **−0.000884** | +0.000922 | priced (6/6 board signs) |
| `mse` | 096 `claim` route, via 197's affine inversion `raw = 1.023611 − 2.074495·cos` | **−0.000509** | +0.000125 | priced |
| `nmae` | 188 identity `1+(U−2M)/D`, lfc space, wald89 significance pool | **−0.002468** | −0.002749 | priced (third opinion, 162: 3/6 signs) |
| `reach` | ordering proxy only | +0.00685 raw* | +0.00890 raw* | **NOT priceable** — rule 9 |
| `fid` | ordering proxy only | −0.01441 raw* | +0.00171 raw* | **NOT priceable** — rule 9 |
| `jac` | ordering proxy only | +0.00347 raw* | +0.00471 raw* | **NOT priceable** — rule 9 |

`reach`/`fid`/`jac` are unpriced and I state plainly why: all three carry `n_real` in a
denominator or a union, `n_real` is 89 on this pool against ~2,000 on the challenge, 162 scores
`reach` at **2/6 correct signs — worse than chance** and `fid`/`jac` at none, and their call sets
require count emission this mandate forbids. The +0.0089 floor **is** the realised RMS of exactly
these unpriced channels for this arm class, which is how rule 13 charges the uncertainty.

**Composed gate:**

```
eb_broad_prior   :  -0.000884 (pds)  -0.000509 (mse)  =  -0.001393   vs +0.0089   VETO
eb_conservative  :  +0.000922 (pds)  +0.000125 (mse)  =  +0.001047   vs +0.0089   VETO
   with the measured nmae folded in:      -0.003861            -0.001702
```

**In 197's units.** The exchange rate is +0.00365 overall per +0.01 cosine, and a single
direction arm must deliver **+31% cosine** (= +0.0246 absolute on our 0.079249) to clear the gate
alone. Measured: `eb_broad_prior` **−0.00147** (−6.7%), `eb_conservative` **+0.00036** (+1.6%,
CI [−0.000056, +0.000791], includes zero). The best setting delivers **5% of what is required**,
and not significantly different from nothing.

**`nmae` at each arm's own re-solved alpha** (the amplitude-sensitive member; 187's ρ cliff):

| arm | α | delivered norm / a15 | `nmae*` | Δ | marginal refund `e` | Δ overall |
|---|---:|---:|---:|---:|---:|---:|
| a15 | 0.982 | 1.0000 | 0.98476 | — | — | — |
| shrink6.78 | 0.982 | 0.1475 | 0.99562 | +0.01087 | 0.545 | −0.003019 |
| **shrink6.78 @ α\*** | 6.658 | **1.0000** | **0.98476** | **+0.00000** | — | **−0.000000** |
| `eb_broad_prior` | 0.982 | 0.2308 | 0.99364 | +0.00889 | 0.541 | −0.002468 |
| `eb_broad_prior` @ α\* | 4.539 | 1.0545 | 1.02365 | +0.03890 | **−0.193** | −0.010804 |
| `eb_conservative` | 0.982 | 0.2080 | 0.99465 | +0.00990 | 0.544 | −0.002749 |
| `eb_conservative` @ α\* | 4.950 | 1.0011 | 0.99566 | +0.01091 | 4.855 | −0.003030 |

The sharpest number in the table is **`e = −0.193`**. 169 §5's threshold is not fitted: `nmae`
improves iff the marginal refund rate `e > 1/2`. At matched amplitude, every unit of L1 the
posterior spends relative to a15 lands on genes whose sign is **wrong**. The pure rescale, by
contrast, sits at `e = 0.545` and its cost vanishes entirely at α\*.

---

## 7. Falsifiers (Phase 4)

1. **Pure-rescale control — PASSES, and it explains everything.** §2: 100.5–104.3% of the
   recorded −0.01505 in 125's own frame, and +0.001073 `pds` / +0.000049 cosine / +0.00000
   `nmae`-at-α\* in production units. **The diagnosis is proven.** The residual — the 15.40% of
   coordinates that flip sign — is the estimator's real contribution, and §4 prices it at
   −0.002372 `pds` and −0.00147 cosine.
2. **125's frozen nulls.** Nothing is promoted, so no endpoint needs to beat them. Had I
   promoted `eb_conservative`'s `pds`, its own numbers already veto it: +0.002474 with a 95% CI
   of [−0.004369, +0.009006], P(Δ>0) = 0.765, 3/4 folds positive, and it was selected **post
   hoc** from a grid of four after being chosen in only 1/15 of 125's own inner folds.
3. **The precision endpoint's two standing failures.** 125's family-share gate (Orion supplied
   0.5667 of the positive gain against a 0.5 maximum) and its two complete-procedure nulls
   (p = 0.33, 0.13) still stand. **I did not clear them and I make no attempt to promote on
   precision.** §4.3 supersedes the question anyway: the endpoint they guard measures a ranking
   key that loses to the one we ship.

---

## 8. Verdict, in overall units

| claim under test | verdict |
|---|---|
| "Signed credit" is none of the six and is degree-1 homogeneous | **CONFIRMED** (8.9e−16) |
| The improved endpoint is exactly scale-invariant | **CONFIRMED** (0.0e+00) |
| The estimator shrinks 6.78× | **CONFIRMED** (LS 6.780, L1 5.599, L2 4.820) |
| 142 uses the same name for the scale-invariant version | **CONFIRMED** (`develop.py:1379`) |
| The recorded −0.01505 is the shrink and nothing else | **CONFIRMED** (100.5–104.3%) |
| The shrink is invisible to what the board reads | **CONFIRMED** (+0.001073 `pds`, +0.000049 cos) |
| **The +0.0122 was an unadjudicated DIRECTION gain** | **FALSE.** 86.5% of it is a confidence key, and that key loses to `|H|` by −0.00281 / −0.01429 |
| **The estimator raises cosine or `q_w`** | **FALSE.** −0.00147 cosine (CI excludes 0), −0.00049 `q_w` |

**Best arm in the family, priced: −0.0017 to −0.0039 overall against a +0.0089 gate.**
125's estimator was closed for the wrong reason in August. It is closed for the right reason now.

**What was actually salvaged.** Two things, both worth keeping:

- 125's estimator *is* better than 122's `equal_mean_agreement` (+0.004521 `pds`, 2/4 folds), and
  125's confidence *is* better than 122's confidence (+0.00235 / +0.00387). 125's original
  result was directionally correct about its own comparator. It was never compared to
  production, and production beats both.
- The mechanism (§5): performance across the four settings is explained at r = +0.96 by how much
  plain averaging each one preserves. **Empirical-Bayes shrinkage across panel targets is not an
  estimator improvement here.** That closes the "borrow strength across targets" route inside the
  panel, which §8.2 listed as one of the last two direction levers.

---

## 9. Hand-off (territory, §11.3)

One result is outside my territory and I am not acting on it.

**To 194 (donor tables).** `eq122` — per-row median-|effect| normalisation, average, rescale by
the median row scale — costs **−0.006893 `pds`** (CI [−0.012251, −0.001358], P(Δ>0) = 0.007) and
**−0.00097 cosine** against production's plain equal mean of raw effects. That is a real,
significant, *negative* result about a row-normalisation production does not use, measured on the
board-validated instrument. If 194 is considering per-row normalisation of the donor tables,
this is the price. It is also the reason 125's whole comparison was mis-anchored.

**To 193 (operator).** Nothing. The whitening operator was held at production's
`union_pruned` / `held_*` pair on every fold and never varied.

**Not handed off:** the `reach*`/`jac*` gains in §6 (+0.0069/+0.0035 for `eb_broad_prior`,
+0.0089/+0.0047 for `eb_conservative`). Per §8.2 we are already **ahead of the leader on `fid`
and `jac`**, and `reach`'s ordering channel is delivered at 95% efficiency by 185/191 against
`reach`'s 0.90 purity cap (189, 190). Rule 9 makes these numbers unreadable in board units
anyway. Chasing them would be spending the remaining time in the channel §8.2 measured at 0% of
the gap.

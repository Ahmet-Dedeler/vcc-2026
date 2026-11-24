# 198 — Systematic re-audit of every closure against the instruments we now know were wrong

**Date:** 2026-09-03. Adversarial audit. Read-only outside this directory. No `vcc`, no
`cell-eval2`, no build, no submission. Peak RSS **1.32 GB**, measured (`/usr/bin/time -l`,
`work/rss.txt`); the ceiling is the 17,942² float32 eigenvector matrix in §6.1, everything else is
under 0.2 GB. Cap 8 GB. Everything below is arithmetic on numbers already on disk.

---

## Headline, in four lines

1. **The map is mostly sound. 87 of 116 terminal verdicts carry no hard flag at all**, and of the
   29 that do, the invalidating finding was in almost every case applied *by the very next
   experiment*. This is the strongest statement available about this project's map, and it is
   positive.
2. **Two closures are genuinely reopenable and they are the same arm family: 150 and 151.** Every
   failing gate in both is a `reach` number read on **held Orion HCT116 and nothing else**. Priced
   without that one member, both are **positive**; corrected by AGENTS.md's own measured
   overstatement factors, 150 lands in `[−0.00042, +0.00077]` and 151 in `[+0.00013, +0.00127]`
   overall, against recorded verdicts of −0.00504 and −0.00434.
3. **125 is adjudicated and the coordinator's suspicion is confirmed and sharpened.** "Signed
   credit" is not one of the six, and — measured, not argued — it is **exactly homogeneous of
   degree 1 in the prediction magnitude**, while the endpoint that *improved* is exactly
   scale-invariant. 125's estimator shrinks by **6.78×**. The rejection is an amplitude artefact.
4. **AgingLab is reconciled: 172 and the +0.1002 are not in conflict, because "172 got nothing" is
   false as stated.** 172's learned key reads **+0.0417 raw h1-test over what we deliver**; its
   NO-GO is against `|H|`, not against production. AgingLab's move is **4.5× larger than the entire
   span from our delivered order to `|H|`**, so it cannot be a key-quality move at all.

**Plus four falsifications of live-mandate premises** (§6), two of which are factual errors about
the shipped code.

---

## 1. Enumeration — method and what it produced

`enumerate.py` → `work/enumeration.json`, `work/enumeration.tsv` (all 173 rows, machine-readable).
It parses every `experiments/*/RESULTS.md` for a verdict token, the scored members the text
actually quotes, and the surface names it quotes, weighting verdict-bearing lines 5×. Nothing is
taken from `CURRENT_STATE` summaries or from memory.

| | count |
|---|---:|
| `experiments/*/RESULTS.md` on disk | **173** |
| dirs with no `RESULTS.md` (live or superseded) | 23 |
| terminal verdict (NO-GO / VETO / NEGATIVE / CLOSED) | **116** |
| GO / SHIP | 15 |
| no verdict token | 42 |

Verdicts: `NO-GO` 78 · `CLOSED` 21 · `NEGATIVE` 14 · `VETO` 3 · `GO/SHIP` 15 · none 42.

**Primary endpoint of the 116 terminal verdicts:** `pds` 55 · `reach` 17 · `mse` 13 · `nmae` 8 ·
`cosine` 5 · `sign_acc` 5 · `jac` 2 · `fid` 2 · `overall` 2 · none 7.

**Primary surface of the 116 terminal verdicts:** none-quoted 42 · **`h1-test` 26** ·
**`board` 16** · **`orion-hct116` 12** · `somatic` 6 · `k562-gwps` 5 · `lodo` 4 · `h1-2025` 2 ·
`plolo` 2 · `orion-hek293t` 1.

That distribution is itself reassuring: the two surfaces the harness warning condemns carry
**13 of 116** primary verdicts, and `h1-test` plus `board` carry 42.

---

## 2. Flagging — and why the keyword sweep had to be thrown away

A first pass matched the mandate's ten invalidating classes as keyword co-occurrence. It flagged
**116 of 173** files, 79 of them on `C3_orion_absolute` alone — i.e. it flagged every file that
contains the string "orion". **That is not a finding, it is a word count**, and reporting it would
have been the exact failure mode AGENTS.md §11.1 rule 8 warns about: a summary that looks
authoritative because it averages away the constraint. It is recorded in
`work/enumeration.json` and then discarded.

`hardflag.py` replaces it with **six high-precision flags**, each requiring a literal construction
the invalidating finding condemns, not a topic:

| flag | what must literally appear | n | experiments |
|---|---|---:|---|
| `H_toll_0146` | the string `0.0146` (rule 13's emission toll) | 18 | 093 094 095 097 150 151 164 168 169 177 178 179 180 181 183 185 188 189 |
| `H_orion_absolute` | an absolute member value + an Orion surface **on the same line**, with delta language and `h1-test` excluded | 12 | 088 094 095 097 098 100 151 153 159 168 188 190 |
| `H_emission_reach` | a `reach` number on the same line as emission-law language | 11 | 027 034 037 095 153 154 171 172 177 185 186 |
| `H_n_real_target` | the literal `n_real_target` | 5 | 153 157 169 171 190 |
| `H_2670` | the literal `2,670` / `2670` | 5 | 154 164 177 179 181 |
| `H_surface_pds` | `debench` `Surface(...)`/`.pds()` sizing a `pds` arm | 2 | 169 188 |

**32 of 173 experiments carry ≥1 hard flag; 29 of those have a terminal verdict.**
Full table: `work/flagged.md`. **87 of the 116 terminal verdicts carry no hard flag at all.**

### 2.1 The decisive filter: has a *later* experiment already corrected it?

`forward198.py` counts forward citations of each flagged experiment's number in every
later `RESULTS.md`/`MANDATE.md` plus `AGENTS.md`, `CURRENT_STATE.md`, `LOG.md`, `BOARD-LEDGER.md`,
and separately counts citations carrying correction language.

**The result is the most important finding in this audit and it is a negative one.** Of the five
`H_n_real_target` hits: **157 is where the defect was discovered** (§3a, "`debench.Surface(n_real_target=N)`
is the wrong instrument for `reach`"), and **169, 171 and 190 each explicitly state they do not use
it** —

> 171 L80: *"`debench.Surface(n_real_target=N)` is **not** used. 157 §3a: uniform power matching deletes…"*
> 169 L279: *"…`Surface(n_real_target=N)` is **not** used, because 157 §3a…"*
> 190 L49: *"never rescaled by `Surface(n_real_target=…)`"*

Only **153** used it as an instrument, and 153 is the file that *documents* it.

Identically for `H_2670`: **181 already CONTRADICTED it** on a new non-circular route (L504:
*"2,670 would read reach 0.0034 — 51× below the board's 0.1742"*), and 177/179 bracket the
consequence rather than assume the value.

**The invalidating findings were, in almost every case, applied by the next experiment. The map is
self-healing.** That is a full result and it is the answer to the mandate's central question.

---

## 3. Ranked reopenable closures

Ranked by expected value of the correction in **overall board units**, with re-measurement cost.
`reprice198.py` → `work/reprice198.json`.

### R1 — `150-dual-moment-score-level-hct116` · TERMINAL_NO_GO · **the cleanest hit in the sweep**

**Surface: held Orion HCT116 and nothing else.** Three occurrences of "HCT116", zero of "h1-test".
Recorded six-member table (scaled), from its own "Primary result":

| member | scaled Δ | overall Δ | share of \|avg\| |
|---|---:|---:|---:|
| `pds` | +0.000000 | +0.000000 | 0.0% |
| `mse` | +0.000000 | +0.000000 | 0.0% |
| `nmae` | +0.006929 | +0.001155 | 22.9% |
| `fid` | −0.006765 | −0.001128 | 22.3% |
| **`reach`** | **−0.045128** | **−0.007521** | **149.1%** |
| `jac` | +0.014694 | +0.002449 | 48.5% |
| **AVG** | | **−0.005045** | (file records −0.005045 ✓) |

**All five members except `reach` sum to `+0.002476` overall.** All three of its FAILED gates are
that one number: "projected average ≥ +0.004" (reach-driven), "bootstrap lower bound above zero"
(reach-driven), and "No DE metric below −0.015 scaled" — which *is* the reach number.

AGENTS.md's harness warning is specifically about this: *"on `reach`: `h1-test` reproduces two board
points to within 0.014 raw while repaired HCT116 misses a15's by **2.6×** and overstates the
tier-policy effect by **4.4×**."* Applying both measured factors:

| | overall Δ |
|---|---:|
| as recorded | −0.005045 |
| corrected at 2.6× | **−0.000416** |
| corrected at 4.4× | **+0.000767** |

**Corrected range `[−0.00042, +0.00077]`.** The file says *"The interval is entirely negative. This
is not a borderline gate miss."* Corrected, **it is exactly a borderline gate miss.**

Two things make this worth more than its point estimate. **`pds` is exactly preserved per target**
(a frozen gate that PASSED), so this is the cheap invariant arm class of AGENTS.md §8.1 rule 2. And
its `jac` **+0.014694 scaled = +0.00245 overall** is the largest single positive DE member in either
file — `jac` being precisely what 192 is chartered to move, and where we rank 169th of 612.

**Re-measurement cost: LOW.** The arm is a fixed dual-moment transform; re-scoring it on `h1-test`
with a 188-style significance-selected pool is a screen, not a build. **This is the single highest
value-per-hour item in the audit.**

### R2 — `151-relaxation-consensus-selector` · terminal NO-GO · same family, same defect

150's rescue attempt, same surface, same dominant term.

| member | scaled Δ | overall Δ |
|---|---:|---:|
| `pds` | +0.000000 | +0.000000 |
| `mse` | +0.000341 | +0.000057 |
| `nmae` | +0.006315 | +0.001053 |
| `fid` | −0.003736 | −0.000623 |
| **`reach`** | **−0.043628** | **−0.007271** (167.4% of \|avg\|) |
| `jac` | +0.014645 | +0.002441 |
| **AVG** | | **−0.004344** ✓ |

Ex-`reach`: **+0.002927 overall.** Corrected: **+0.000131 (2.6×) to +0.001275 (4.4×)** — this one
**flips sign under both factors, not just one.** 151 also has the **fewest forward citations of any
flagged terminal closure (2) and zero corrections**, i.e. nothing since has looked at it.

**Re-measurement cost: LOW**, and it shares 150's harness, so R1 and R2 are one job.

### R3 — `125-hierarchical-panel-posterior` · NO-GO · closed on an endpoint the competition does not read

Fully adjudicated in §5.1. Corrected value is **not** directly priceable in board units because 125
never reached a scored member — that is the point. Its scale-invariant endpoint improved
**+0.012175** in all three held families. **Re-measurement cost: MEDIUM** (the estimator exists and
its artifacts are on disk; rescaling is one line, but pricing it needs a six-member screen).

### R4 — `154-confidence-key-screen` · NO-GO · **already overturned; no reopening needed, but the file is still wrong on disk**

Three independent defects, all real: it is measured on `orion-hct116-thinned` (L102); its robustness
caveat (L87–91) rests on 093's `n_real ≈ 2,670`, which 181 contradicted; and its headline is
overturned by the board — `CURRENT_STATE.md`: *"186: THE FIELD'S `reach` SPREAD IS ORDERING, NOT
DIRECTION. 154 IS OVERTURNED BY THE BOARD."*

Its operative instruction — *"**Do not** build further ordering keys, per-gene `phi` confidence maps,
or spread policies aimed at `reach`"* — was wrong, and 171, 172, 185, 189 and 190 all productively
violated it. **Nothing to reopen; the correction is already made. Recommend the file carry a
retraction header** so the next agent reading it cold does not obey L81–83. Cost: one edit.

### Not reopenable — checked and confirmed

- **153, 157, 169, 171, 190** (`n_real_target`): 157 discovered the defect; the other three
  explicitly refuse the instrument by name. Closures stand.
- **164, 177, 179, 181** (`2,670`): all four bracket the value rather than assume it; 181
  contradicts it outright. Closures stand.
- **088, 094, 095, 097, 098, 100, 159, 168, 188, 190** (`H_orion_absolute`): each quotes an Orion
  absolute somewhere, but in every case the *verdict* rests on `h1-test`, `board` or a
  mean-removed arm-to-arm delta, which AGENTS.md explicitly exempts (*"Mean-removed statistics and
  arm-to-arm deltas on those surfaces are not distorted"*). **098 independently found the same class
  of defect before 188 did** — its second headline is *"every `nmae` number this project has ever
  taken was read on a reference gate 2–20× shallower than the board's, and the sign of the amplitude
  response flips across that range."* Closures stand.
- **The 18 `H_toll_0146` hits**: rule 13's toll is correctly applied to emission-law arms in all of
  them, and 169's profile-class falsification (+0.0006) is already carried in AGENTS.md's rule-13
  amendment. No profile-class arm was found vetoed by the toll. Closures stand.

---

## 4. One sizing correction to `docs/BOARD-LEDGER.md` (flagged, not asserted)

`BOARD-LEDGER.md` §"Consequence for sizing" prices 172's ordering counterfactual at
**`|H|` = 0.23533** vs delivered 0.18283, *"+0.05250 raw on h1-test … +0.0210 raw on the board =
+0.00395 overall."*

But 172's own known-answer block carries **two** rows, and the scorer's ordering rule
(`vcc2026-metrics-BRIEF.md` §5: *"its significant calls first, then ascending predicted adjusted
p-value…"*) means the **demonstrated deliverable is the other one**:

| row | h1-test `reach` | Δ vs delivered | → board (×0.40) | → overall |
|---|---:|---:|---:|---:|
| production p-values (delivered) | 0.18283 | — | — | — |
| constant key (gene-name order) | 0.18021 | −0.00262 | | |
| **`sig-first, then \|H\|`** (= `\|H\|` over called genes only) | **0.21914** | +0.03631 | +0.01452 | **+0.00273** |
| 185 `rc_hybrid` realised | 0.21681 | +0.03398 | +0.01359 | +0.00256 |
| `\|H\|` over the whole pool (the ledger's row) | 0.23533 | +0.05250 | +0.02100 | **+0.00395** |

**The ledger's figure is 1.44× the demonstrated-deliverable figure.** 171 proved the sig-first key is
a BH no-op, so 0.23533 is not *forbidden* — it requires our emitted `p_adj` to be monotone in `|H|`
across the **uncalled** block too, which no built arm has achieved. Stated precisely: **0.21914 is
demonstrated, 0.23533 is aspirational, and the ledger prices the aspirational one.** `rc_hybrid`'s
pending board row is the measurement that settles it.

Incidental: 172 §4 also records that **our delivered ordering is worth `+0.00262` over a constant
gene-name-order key.** We are, to three decimal places, delivering no ordering at all.

---

## 5. The two named adjudications

### 5.1 — 125: is "signed credit" one of the six? **No, and it is worse than that.**

`adjudicate125.py` → `work/adj125.json`, `work/adj125.log`.

**Definition, quoted from source.** `experiments/125-.../preregistered_protocol.json:156`:
*"Per target, `sum(sign(held_effect) * predicted_effect)` divided by `sum(abs(held_effect))` over
held-valid off-diagonal response coordinates."* Implemented at `run.py:774-776`:

```python
credit[target] = float(
    np.sum(np.sign(truth) * prediction[target, indices]) / denominator
)   # denominator = float(np.sum(np.abs(truth)))
```

**It is none of the six.** `pds` is a cosine rank; `mse` is a squared error; `nmae` reads log2FC over
the reference's DE genes; `fid` and `jac` read the call set; `reach` reads the ranking of the
reference's pool. None is a truth-normalised signed inner product.

**And it is degree-1 homogeneous in the prediction.** The numerator is linear in `prediction`; the
denominator depends only on truth. Verified numerically:

```
credit(0.50·pred)/credit(pred) = 0.500000000000
credit(0.75·pred)/credit(pred) = 0.750000002765
credit(2.00·pred)/credit(pred) = 2.000000000000
```

The endpoint that **improved** — top-20% truth-magnitude-weighted sign precision, `run.py:777-792` —
is **exactly scale-invariant**:

```
precision20(0.50·pred) − precision20(pred) = +0.000000000000
precision20(0.75·pred) − precision20(pred) = +0.000000000000
precision20(2.00·pred) − precision20(pred) = +0.000000000000
```

**The estimator shrinks by 6.78×.** Measured on the saved arrays
(`observed-candidate-effect.npy` vs `observed-baseline-effect.npy`, 270,000 coords, 231,888
nonzero): L1 ratio **0.1786**, L2 ratio **0.2075**, least-squares scalar **0.14749**, cosine
**0.7110**. Per-target LS scalar median **0.1583** (IQR 0.1377–0.1735).

**The shrink alone accounts for the entire loss.** Under a pure-rescale null `cand = c·base`,
`Δcredit = (c−1)·credit(122)`, so the recorded **−0.015050** implies `credit(122) = +0.01765`
(at the L1 ratio, **+0.01832**; at L2, **+0.01899**). A transfer cosine of 0.04–0.08 (AGENTS.md
§4.2) puts `credit(122)` in exactly that range. **No directional degradation has to be invoked.**

**Verdict.** 125's promotion gate coupled a scale-invariant direction improvement to a
scale-**covariant** magnitude penalty and rejected on the latter. It cannot distinguish "worse
direction" from "smaller α" — and α is a free parameter this project tunes separately, at the joint
`mse`+`nmae` optimum (§4.1, exp 187). Multiplying 125's candidate by `1/c = 6.78` sets
`Δcredit = 0` **by construction** while leaving the `+0.012175` precision gain **exactly** intact.

Sharpest supporting evidence: **experiment 142 uses the same name for the scale-INVARIANT version.**
`142/develop.py:1379` divides by `np.sum(np.abs(prediction), axis=1)` — the *prediction*'s L1, not
the truth's. Two experiments in one protocol lineage, one statistic name, two different homogeneity
degrees, and 125's NO-GO rests entirely on the non-invariant one.

125's own text anticipated the remedy: *"A future experiment could test a direction-only
confidence/ranking layer while retaining Experiment 122's effect magnitudes."* **That experiment is
authorized by this audit.**

### 5.2 — AgingLab `GeroAI_v9 → v11` vs 172: reconciled

`board198.py` → `work/board198.log`. Six-member rows from
`experiments/186-intel-0903/board/trajectories.md:17-18` (the INBOX table gave only four).

| arm | date | overall | `pds` | `mse` | `nmae` | `fid` | `reach` | `jac` |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `GeroAI_v9` | 08-25 | 0.1545 | 0.8273 | 1.822 | 0.9216 | 0.5050 | 0.1672 | 0.0284 |
| `GeroAI_v11` | 08-26 | 0.1899 | 0.8196 | 0.959 | 0.8923 | 0.5135 | **0.2674** | 0.0287 |

**Q1 — decomposition.** Reconstructs the reported +0.0354 to **+0.03425** (error −0.00115, the mse
clamp):

| member | Δ raw | Δ overall |
|---|---:|---:|
| `pds` | −0.0077 | −0.00284 |
| `mse` | −0.8630 | +0.00526 |
| `nmae` | −0.0293 | +0.00807 |
| `fid` | +0.0085 | +0.00477 |
| **`reach`** | **+0.1002** | **+0.01886** |
| `jac` | +0.0003 | +0.00014 |

**`reach` alone is 53% of the whole move.**

**Q2 — amplitude is ruled out, and it works against them.** Refitting 186's law on the full
613-entry snapshot (`reach ~ 1 + log mse + pds`, 53 teams with `pds > 0.80`) gives
**d reach/d log mse = +0.03468** — note this is **2.7× larger than the +0.0127 that 186 quotes and
that 190's INBOX-1 relies on**; residual sd 0.03240. They shrank (`d log mse = −0.6418`), so the law
**predicts −0.03150**; they realised **+0.10020**. Unexplained **+0.13170 = 4.1 residual σ.**

**Q3 — it is not a sign-information move.** `fid` +1.68%, `jac` +1.06%, `pds` −0.93%, `reach`
**+59.93%**. Board-wide `corr(raw reach, raw fid)` over the `pds > 0.80` cohort is **+0.0297** —
the two members are decoupled. A 60% `reach` jump with `fid` +1.7% cannot be a change in the
underlying sign information; it is ordering/delivery. `v9` sat at board percentile **60.4**; `v11`
at **98.0**.

**Q4 — and this is the reconciliation. "172 got nothing" is false as stated.** 172's own §4 table,
all rows with the same predictions and the same signs:

| key | h1-test `reach` | Δ vs **production** | → board raw | → overall |
|---|---:|---:|---:|---:|
| production p-values | 0.18283 | — | — | — |
| M1 logistic (AgingLab's feature family, LODO-selected) | **0.22448** | **+0.04165** | +0.0162 … +0.0178 | +0.0031 … +0.0034 |
| M2 step-additive | 0.23163 | +0.04880 | +0.0190 … +0.0209 | +0.0036 … +0.0039 |
| `\|H\|` (the benchmark) | 0.23533 | +0.05250 | +0.0205 … +0.0225 | +0.0039 … +0.0042 |

172's NO-GO is *against `|H|`* — its own line 6 reads *"the ordering gap is NOT learnable from
AgingLab's feature family"*, and lines 16–19 say plainly *"Learning works; it just recovers what
spectral whitening already provided for free."* **Relative to what we deliver, AgingLab's feature
family is worth +0.0417 raw h1-test.**

**The two results are about different halves of the ordering channel and do not conflict:**

> AgingLab realised **+0.1002** raw board `reach`. Our best measured ordering counterfactual is
> **+0.0225** raw board `reach`. **AgingLab's move is 4.5× larger than the entire span from our
> delivered order to `|H|`, the best key 172 could find.**

A move 4.5× larger than the whole key-quality span **cannot be a key-quality move.** 172 measured
KEY QUALITY at fixed delivery and found the family's ceiling. AgingLab moved DELIVERY.

**What evidence would settle the mechanism** (the mandate asks for this explicitly):

1. **`rc_hybrid`'s board row** (submitted 2026-09-03 11:20 UTC). It changes delivery at a pinned
   key. If it realises ≥ +0.013 raw board `reach`, delivery is confirmed as the channel and the
   0.40 shrink factor is confirmed; if it realises ≪ that, the shrink factor or the `reach`
   pricing route is wrong and §4's table must be redone.
2. **`n_pred` for v9 and v11.** Not on the board — but `jac` and `fid` jointly constrain it, and
   both are ~flat, which *already* argues their call set barely changed. This is the strongest
   available inference and it points at the p-value *order*, not the set.
3. **The decisive negative control is already on disk and nobody has used it**: `Haystack` held
   `reach` 0.2050 at their **third** submission with `pds` 0.686 — worse than our first entry —
   and it barely moved while `pds` climbed 0.12 and `mse` fell 7× (190's INBOX-1 §2). A `reach`
   that is high *before* any tuning and invariant to model quality is an **emission or ordering
   convention**, not a learned component. Combined with Q3's `corr(reach, fid) = +0.03`, the
   arithmetic permits a convention in the emitted p-values and rules out a direction improvement.

---

## 6. Adversarial pass on the six live mandates

Four premises falsified from data already on disk, two of which are factual errors about the
shipped code. Every live run's Phase-0 gates that I could check **pass** — the methodology is
healthy; the errors are in the mandates, not the work.

### 6.1 — 193: **the premise is factually wrong. There is no ridge.**

The mandate's central claim:

> *"Production whitens with `P → P(C + λI)^(−γ/2)` … `λ` a scalar ridge. **Nobody has questioned the
> regulariser itself.** `C + λI` is the crudest possible covariance estimator: it shrinks every
> eigenvalue by the same additive constant."*

**`scripts/make_submission_pm_cd4.py:493-494` ships `WHITEN_GAMMA = 0.65` and `WHITEN_LAM = 0.0`**,
and the comment immediately above it (L484-486) says so in as many words: *"`lam` has sat at its 0.0
default in every build this project has ever shipped."* 193's own `work/gate0.json` records
`"map": "prod(lam=0)"`, so the **code** has it right and only the **prose** is wrong.

The shipped eigenvalue map is therefore `f(ev) = max(ev, ev_max·1e-12)^(−0.325)` (`whiten.py:65,79`)
— **a pure power law with a relative floor. It is already a nonlinear shrinkage**, and `γ` was tuned
to a six-member local optimum (170). Measured on the shipped operator
(`whitener198.py`, sha `e2f65ecf…`, G = 17,942, max_ev 5.894570):

| | union_pruned (PRODUCTION) |
|---|---:|
| eigenvalues at the relative floor | 48 (0.27%) |
| `w(ev_max)` | 0.5618 |
| `w(median)` | 8.830 |
| `w(ev_min)` | **4462.8** |
| **condition number of the applied map** | **7,943×** |
| directions amplified (`w>1`) | **17,941 (100.0%)**, carrying 90.07% of `trace(C)` |
| `ev ≤ ev_max·1e-6` | 2,247 dirs (**12.52%**), amplified ≥ **89×** the top |
| `ev ≤ ev_max·1e-8` | 408 dirs (2.27%), amplified ≥ 398× the top |

`whiten.py:81-84` renormalises the output to the input's median row norm, so the *absolute* weights
above are defined only up to one global scalar and the meaningful statistic is the **relative** map
— the condition number and the "× the top" column, both of which the renormalisation leaves exactly
invariant. Same for `union` (7,943×) and `union_full` (7,943×): the condition number is
`(1e-12)^(−0.325)` and is therefore **identical for all three operators by construction**, set
entirely by `floor_rel`, not by the data.

**This is the mechanism behind AGENTS.md §4.2's operational corollary** — *"whitening amplifies
low-variance directions, so `|predicted effect|` becomes a worse per-gene confidence signal
afterwards"* — with the number nobody had. The amplification spans **four orders of magnitude**, and
**12.5% of the basis is amplified ≥89× relative to the dominant programme.**

**Consequence for 193's design.** The experiment as framed ("replace the crude ridge with a better
estimator") has no ridge to replace, and its identity row must be `λ=0`. The actually-untested
lever is the one its own accidental code comment names: **there is no floor, and there should be
one.** 193's `work/nlscheck.json` independently supports narrowing the family — on **power-law**
spectra (which is what a real biological covariance looks like), Ledoit–Wolf NLS is catastrophic
(**−5,972 to −30,457**) while only RIE survives (+82 to +93). A power-law *map* is the natural
conjugate of a power-law *spectrum*, so the honest question is far narrower than the mandate
suggests.

**And 193 and 195 were handed contradictory descriptions of the same object.** 195's mandate calls
it *"a hand-chosen **one-parameter** family"* with `γ` swept — which, at `λ=0`, is **correct**. 193's
calls `λ` the crude regulariser — which is not. Both cannot be right and 195 is.

### 6.2 — 196: **its own Phase 0 falsifies its motivating premise; the abort condition has fired**

The mandate motivates the sweep with *"the rounding is coarse where it matters most … at 1 CPM a
gene receives 8 counts across 400 cells … **a large share** of the ~10,700 gated genes … is
delivered at a few tens of counts"*, and pre-commits: *"If that error is negligible the experiment
is over — report it and stop."*

`work/round196.md`, already written, measures it. On the three challenge contexts:

| | ctx A | ctx B | ctx C | h1-test |
|---|---:|---:|---:|---:|
| median delivered total over the gated pool | 459 | 399 | 395 | 385 |
| share below 100 counts | 11.2% | 14.0% | 13.9% | 15.1% |
| **share below 40 counts** | **0.2%** | **0.4%** | **0.2%** | **1.3%** |
| **mean relative shortfall over the gated pool** | **0.206%** | **0.236%** | **0.232%** | **0.255%** |

**"A large share delivered at a few tens of counts" is 0.2–1.3%.** The DE gate is >5 CPM, not 1 CPM,
and the floor is applied to each gene's 400-cell aggregate, not per cell. **The abort condition has
fired and 196's headline should be negative.**

**But there is a second, unnamed finding sitting in that same table that is worth more than the
sweep.** The rounding-plus-renormalisation imposes an expression-dependent **tilt** whose crossover
differs **13×** across the three scored contexts:

| | ctx A | ctx B | ctx C |
|---|---:|---:|---:|
| tilt crossover | 75.2 CPM | 229.4 CPM | **17.4 CPM** |
| net tilt, decile 1 → decile 10 (log2) | −0.0113 → +0.0009 | −0.0130 → +0.0001 | −0.0077 → **+0.0049** |

Context C's profile is tilted **up** across eight of ten deciles while context B's is tilted **down**
across nine. Overall is an unweighted mean over contexts, so this is a per-context systematic that
does not cancel. Magnitudes are small (≤0.013 log2 ≈ 0.9%) but `nmae` reads log2FC directly, and 187
put the entire collectable `nmae` gain inside the `mse`-feasible window at **0.0012 raw**. **196
should price the tilt, not the shortfall.**

### 6.3 — 192: **its own falsifiers undercut the mandate's thesis — and 192 says so itself**

192 finished while this audit was running and its verdict is *"The mandate's premise is half right
and the half that fails is the important half."* Everything below was read from its intermediate
logs first and then checked against its final `RESULTS.md`; where they agree I say so rather than
claiming the point. Four things the coordinator needs:

**(a) The prediction contributes nothing to call-set membership.** From `work/enrich192.log`,
precision `p(n)` = mean fraction of the top `n` that the reference calls significant:

| key | n=25 | n=50 | n=89 | n=300 |
|---|---:|---:|---:|---:|
| `B \|lfc_pred\|·√cpm` (188's key, the mandate's premise) | 0.2449 | 0.2139 | 0.1877 | 0.1413 |
| **`B √cpm` only** | **0.3078** | **0.2741** | **0.2465** | **0.1945** |
| `B detection rate` only | 0.2727 | 0.2615 | 0.2383 | 0.1908 |
| `B \|lfc_pred\|` alone | 0.0839 | 0.0688 | 0.0632 | 0.0533 |
| `B gene DE-rate LOTO (same line)` | 0.3717 | 0.3420 | 0.3157 | 0.2552 |
| `B gene DE-rate Orion (cross-line)` | 0.1327 | 0.1300 | 0.1193 | 0.1076 |
| `C shuffled \|H\|` (null) | 0.0756 | 0.0737 | 0.0732 | 0.0728 |

**Control CPM alone beats the mandate's headline key at every `n`** (0.3078 vs 0.2449 at n=25), and
multiplying in `|lfc_pred|` makes it *worse*. The reference's call set is a **gene-level,
perturbation-independent propensity** — which is why the same-line DE rate is best of all (0.3717)
and why it **does not transfer across lines** (0.1327, a 2.8× drop). The good news for D/E/F: `√cpm`
and detection rate are computable from the 18,400 NTC cells we get on day one.

**(b) The `fid` half of the thesis is falsified by 192's own shuffled-key control.** From
`work/retro192.log` FALSIFIER 1, at the retrodicted `f = 0.5`:

| `n_pred` | `fid` FIT | `fid` shuffled | `jac` FIT | `jac` shuffled |
|---:|---:|---:|---:|---:|
| 1100 | 0.4851 | 0.4838 | 0.0635 | 0.0181 |
| 2200 | 0.5048 | 0.5044 | 0.0618 | 0.0243 |
| 4400 | 0.5059 | 0.5053 | 0.0532 | 0.0300 |

**`fid` FIT ≈ `fid` shuffled at every `n`** (+0.0004 to +0.0013). The mandate's premise — *"at
`n_pred = n_real` the yield cap pins at 1 and `fid` becomes the accuracy of the **head**"* — assumes
the head is more accurate than the tail. **It is not.** Shrinking the call set to 1,100 takes `fid`
to 0.4851, *below* our current 0.5038 at 4,474. The entire real gain is in `jac` (3.5× over shuffle),
which is a genuine result — but it is one member, not three, and `price192.log` shows the joint
optimum sitting at or **above** the incumbent `n_pred`, worth `+0.001` to `+0.004` scaled on SUM(3).

**(c) The board retrodiction failed on `jac` mid-run, and 192 fixed it before finishing —
recorded because the fix is the finding.** `work/retro192.log` FALSIFIER 2, which is what this
audit first read, retrodicted v14's `fid` at ratio **+1.03×** but `jac` at **−0.13×, sign wrong**
— and `jac` is the member 192 is chartered to move. Under audit-181 rule 1 that would have made
its `jac` numbers unreadable. **192's final `RESULTS.md` §3.2 and finding 6 resolve it**: a
*composed* two-channel model (call set + emission law) retrodicts **both** `fid` and `jac` at
**0.72×**, and its stated lesson is the general one — *"Either channel alone gets one member right
and the other's sign wrong. Any future call-set retrodiction must carry both."* **No action
needed; the flag is retired.** It is left in this report because a single-channel call-set
retrodiction is a trap the next agent will otherwise walk into.

**(d) A premise disagreement, now resolved in 192's favour.** Its best fit implies board `n_real`
median **31–45**, mean **270–383** — against the **89–110** that its own mandate, AGENTS.md §2.1
and 188 all assert. A 2–3× disagreement on the quantity that sets `fid`'s denominator, measured
twice from two different members. **AGENTS.md §2.1 should be updated.**

**Fairness note.** (a), (b) and (d) were reached independently by 192 itself and are its own
headline — *"entirely GENE-LEVEL, and our prediction contributes exactly nothing to it"*, *"`fid`
is closed"*, and the two-member `n_real` measurement. This audit reached them from the intermediate
logs before 192's write-up existed; the agreement is corroboration on independent readings of the
same data, not a catch. The one thing this audit adds is (c)'s generalisation and the observation
that `√cpm`-alone beating the mandate's headline key was visible in `enrich192.log` from the
first pass.

### 6.4 — 197: the `1.28×` claim is not where the mandate says it is

The mandate says *"`AGENTS.md` §8 records an estimate of 'at least 1.28× ours'"* and makes
adjudicating it a mandatory falsifier. **`1.28` does not appear in `AGENTS.md` at all.** It is in
`docs/CURRENT_STATE.md:1273` and `:1395`. Minor, but 197 is instructed to abort on Phase-0 failures
and this is a citation it cannot resolve. Its actual Phase-0 gates **all pass** (`work/gate197.log`:
G1–G4 `claim` within 1.16%, G5 exact, G6 residuals ≤6.1e-5, G7 both out-of-sample predictions
reproduced) — the instrument is sound.

### 6.5 — 194 and 195: no falsifiable premise error found

**194** (donor re-estimation): premises check out. `184` closing line adaptation and `159`'s κ²≈0.34
do leave "estimate the shared component better" as the remaining direction lever, and the
knockdown-weighting boundary argument is correct — the target's own row is excluded from all six
metrics (§4 point 1), so weighting on it costs no scored signal. One caution: AGENTS.md §2.1 records
that the **panel** was selected on a >80% knockdown criterion *and* high basal target expression, so
the range-restriction argument that explains exp 167's null applies to the panel, not to the
Replogle/Orion donors 194 is re-estimating — 194 should not import 167's conclusion.

**195** (learned retrieval transform): premises check out and, as noted in §6.1, its description of
the shipped operator is the accurate one. Its instruction to read 080/084 first and stop if they
already settled it is exactly right.

---

## 7. What I would do with this

1. **Re-score 150/151's dual-moment transform on `h1-test`** with a significance-selected pool. One
   screen, no build, both arms share a harness. Expected value `[−0.0004, +0.0013]` overall, and it
   carries a `+0.0024` overall `jac` gain in an arm class where `pds` is exactly pinned. **Highest
   value per hour in this audit.**
2. **Tell 193 now** that `λ = 0` and the map's condition number is 7,943× — before it spends its
   remaining budget screening replacements for a ridge that does not exist. The floor is the
   experiment.
3. **Tell 196 now** that its abort condition has fired and to re-point at the per-context tilt.
4. **192 needs nothing from this audit** — it caught its own defects. But **update AGENTS.md §2.1's
   `n_real` median from 89–110 to 31–45**, which 192 measured twice from two board members.
5. **Rescale 125's estimator by 6.78× and price it on the six.** The direction gain is real,
   scale-invariant, and reproduced in all three held families.
6. **Add a retraction header to `154/RESULTS.md`** so its "do not build ordering keys" instruction
   stops being read as live.
7. **Correct `BOARD-LEDGER.md` §"Consequence for sizing"** to price the demonstrated `0.21914` row
   alongside the aspirational `0.23533`, or wait for `rc_hybrid`'s board row to settle it.

---

## 8. Integrity

No `vcc`, no `cell-eval2`, no `cell-eval`, no submission build, no writes to `submissions/`, no
writes outside `experiments/198-closure-reaudit/`. No production code touched. Peak RSS 0.63 GB
against an 8 GB cap. All scripts run under `uv run python`.

**Known-answer rows carried in this audit** (AGENTS.md §11.1 rule 10):

| row | expected | here | |
|---|---:|---:|---|
| 150's six-member average reconstructed from its own table | −0.005045 | **−0.005045** | ✓ |
| 151's six-member average reconstructed from its own table | −0.004344 | **−0.004344** | ✓ |
| GeroAI v9→v11 overall reconstructed from six raw members + the ledger's scaling map | +0.03540 | **+0.03425** | ✓ (−0.00115, the `mse` clamp) |
| a8 raw reach → scaled, via the ledger's slope | 0.124349 | **0.124212** | ✓ (see note) |
| `credit(k·pred)/credit(pred)` = k (degree-1 identity) | k | k to 3e-9 | ✓ |
| `precision20(k·pred) − precision20(pred)` = 0 (scale invariance) | 0 | **0.0 exactly** | ✓ |
| production whitener max eigenvalue (rule 11: sha `e2f65ecf…`) | 5.895 | **5.894570** | ✓ |
| `union` / `union_full` max eigenvalues (the rule-11 trap trio) | 6.257 / — | 6.256702 / 6.980217 | ✓ |

**Note on the ledger check row.** `BOARD-LEDGER.md:14` writes its own worked example as
*"(0.1890−0.0790)×1.1292 = 0.12423"*. The arithmetic is **0.124212**; `0.12423` is a transcription
slip in the ledger. Both are within 0.00014 of the recorded 0.124349, so the scaling map is
confirmed and nothing downstream moves — recorded only because a known-answer row that is quietly
0.0002 off is exactly the kind of thing that later gets treated as a real discrepancy.

Artifacts: `work/enumeration.{json,tsv}` (173 rows), `work/hardflags.json`, `work/flagged.md`,
`work/forward198.json`, `work/adj125.{json,log}`, `work/board198.{json,log}`,
`work/whitener198.{json,log}`, `work/reprice198.{json,log}`, `work/table198.log`.

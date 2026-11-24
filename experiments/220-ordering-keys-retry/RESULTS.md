# 220 — the ordering-key family is **CLOSED**, and the 77× tightening was **in the wrong axis**

Date **2026-09-04**. Measurement only. No `vcc`, no `cell-eval2`, no build, no `vcc prep`, no
submission, no `submissions/` write, no file outside `experiments/220-ordering-keys-retry/` created
or modified. **No donor table opened.** Peak RSS **3.63 GiB** (`emit220.py`; every other script is
under 0.6 GiB).

---

## Verdict in ten lines

1. **The 77× tightening is real and it reproduces exactly — and it is entirely in the EMISSION-SEED
   axis.** Both published sems reproduce from their own cached tables: `0.00695` on a15's `quasi`
   emitter, `0.00009` on `rc_callset`, **76×**. The **target-sampling** se of the *same arm* moves
   `0.01660 → 0.01531`, i.e. **1.08×, unchanged**. 189's `±0.00424` and 190's bootstrap CI are
   target-axis numbers. **The instrument that got 77× tighter was never the one limiting them.**
2. **Worse than "the tightening does not apply": the two instruments are the SAME OBJECT for a
   counterfactual key.** `reach(key)` reads only `sign(pred_lfc)` and the reference. The per-target
   `reach` vectors for `|H|` and for 189's key are **bit-identical** on `de164.a15.npz` and on
   `de202.rc_callset.npz` (`np.array_equal` True, max |diff| 0), even though 832 of 65,044 pool
   pairs carry a different predicted sign. **189 and 190 were already measured on the new
   instrument.**
3. **189's structural cap is right in its conclusion and WRONG in its arithmetic, and the correction
   is the most useful thing here.** `reach` is *not* a step function of front accuracy: the measured
   transfer curve gives **d `reach` / d acc@10 = +0.143**, and at `c = 0.30` mean front accuracy is
   **0.872 — below 0.90 — while `reach` is already +0.0379**. You do not have to cross 0.90 on
   average; you have to cross it on **more perturbations**. 189's "≤ +0.0005 raw" understates the
   family's response by **23×**.
4. **The family dies anyway, for a different and sharper reason: no available key moves the front.**
   189's key moves acc@10 by **+0.0103**; the gate needs **+0.0547**. Its own headline
   `+0.082` was a θ-tercile contrast *inside* `|H|` decile 1, not a deliverable front gain.
5. **The honest ceiling of the whole gene-level-key family is NEGATIVE.** Per-gene sign accuracy has
   split-half `ρ = +0.011`; used in fold it is worth **+0.0968** raw `reach`, used cross fold
   **−0.0114**. There is no transferable per-gene ordering signal on this surface at all.
6. **Both published effects are carried by one to three perturbations.** Drop the single
   largest-contributing target and 189's `+0.00408` becomes **+0.00166**; drop the best three and it
   is **−0.00191**. 190's `+0.00137` becomes **−0.00016** on dropping one target — **111% of its
   effect is a single perturbation** — and its leave-one-target-out λ estimate is **−0.00015**.
7. **The mandate's loophole is real in direction and closed by learnability.** 189's gain *is*
   concentrated where `reach` lives (160% of it inside the 31 targets carrying 94.4% of `reach`),
   and a perfect target gate would reach **+0.00843 raw = +0.00205 overall**, just over rule 13's
   `+0.0019`. But no stratifier computable from our own side correlates with it above `|r| = 0.10`,
   and the best split-half-honest gate is worth **+0.00044 raw = +0.00011 overall**, **17× below the
   gate**.
8. **190's one un-null-tested number is now null-tested and it dies.** Its `fid` row — precision
   `+0.00800` at `n = 50` on this instrument — reads **+0.00997 ± 0.00199 for a sign-shuffled
   network** (z = −0.99, P(null ≥ real) = 0.85).
9. **DELIVERED through the shipped emitter — the one thing `rc_callset` genuinely buys — both keys
   are worth about a thousandth of a `reach` point, and delivery attenuates them 4.5×.** 189:
   **+0.000896** (three-seed paired sem **0.000068**) against a counterfactual `+0.00408`; 190:
   **+0.001002**. On 188's significance-selected pool, whose `n_real` matches the board's, **both go
   negative** (−0.00778 and −0.00532, at every seed). Best six-member price **+0.000317** against
   rule 13's **+0.0019** — 6× short, with `pds`/`mse`/`nmae` pinned and the per-gene panel-total
   digest identical on **85/85** blocks for all 33 delivered arms.
10. **The mandate's decisive falsifier fires in all four cells.** Twelve delivered θ-shuffles and six
    delivered regulon-sign shuffles: `P(shuffle ≥ real)` = **0.167 / 0.833 / 0.750 / 1.000**. **190's
    real network is 1.2σ *worse* than its own sign-shuffle, beaten by 5 of 6 draws on one pool and
    6 of 6 on the other.** The mandate pre-committed: *"if a control still ties or beats the arm, the
    family is dead for good and that is the result."*

**VERDICT: the ordering-key family is closed, permanently.** Not by the 0.90 bar, which is
crossable, but by the absence of any per-(target, gene) confidence signal that beats `|H|`. The
closure now carries a **specification** instead of an argument (§7): a key clears rule 13's gate iff
it is a per-pair confidence with accuracy **≈ 0.639** on *"is our sign right"*, i.e. correlation
**≈ 0.30** with correctness inside the reference DE pool. `|H|` itself reads **0.074**; θ reads
**0.002**; expression **0.008**.

---

## 0. Phase 0 — every mandated known-answer row, on both instruments

`ka220.py`, `work/ka220.json`. **Every script in this directory calls `harness172.check` and aborts
on any failure** (rule 10). All rows passed on every run.

### A. the OLD instrument — 189/190's surface, `de164.a15.npz`, counterfactual keys

| row | here | expected |
|---|---:|---:|
| `reach(production order)` | **0.18283** | 0.18283 ✓ |
| `reach(\|H\|)` | **0.23533** | 0.23533 ✓ |
| `reach(sig-first, then \|H\|)` | **0.21914** | 0.21914 ✓ |
| `reach(ORACLE)` | **0.64553** | 0.64553 ✓ |
| **189's key `\|H\| × exp(+0.50 θ)`** | **0.23941**, Δ **+0.00408** | 0.23941 / +0.00408 ✓ |
| …its paired bootstrap over targets | se **0.00422**, P(better) **0.838** | 189: 0.00424 / 0.828 ✓ |
| **190's key `ct-prop-k25`, λ = −0.25** | **0.23670**, Δ **+0.00137** | 0.23670 / +0.00137 ✓ |
| …its 95% CI | **[−0.00137, +0.00516]** | 190: [−0.00140, +0.00520] ✓ |
| **H-oracle, plain 0/1 key** | **0.61862** | 172 §4.2 ✓ |
| **H-oracle, `\|H\|` breaking ties inside the block** | **0.62119** | the INBOX / 171 value ✓ |

*(Both H-oracle conventions reproduce to 5 decimals, so the 0.61862-vs-0.62119 discrepancy 189
recorded as "a tie convention" is exactly and only that. Gaps in this document are quoted against
0.61862 where 172's convention applies.)*

The 190 network rebuilds to `n_tfs = 593`, `topk = 25`, **42,022** non-zero pairs from
`collectri_rescued.csv` (`86c90b30f2cc75c1…`, the digest 190 recorded).

### B. the NEW instrument — `rc_callset`'s own emitted h1-test table

| row | here | expected |
|---|---:|---:|
| `reach(production)` — the **DELIVERED** order | **0.223341** | 214's 0.223341 ✓ |
| `fid` / `jac` / `nmae` / median `n_pred` | 0.494362 / 0.072975 / 0.966702 / 2,634 | 214 ✓ |
| `reach(\|H\|)` counterfactual | **0.23533** | identical to the old table |
| `reach(ORACLE)` | 0.64529 | 0.64553 on the old table |

### C. the two uncertainties — the mandate's premise, tested

| emitter | arm | three-seed mean | **SEED sem** | **TARGET bootstrap se** |
|---|---|---:|---:|---:|
| a15 `quasi` (200) | drop PerturbMatch | **+0.01491** | **0.00695** | **0.01660** |
| `rc_callset` (214) | drop PerturbMatch | **+0.02612** | **0.00009** | **0.01531** |
| | | | **76× tighter** | **1.08× — unchanged** |

Per-seed deltas: old `+0.028764 / +0.007031 / +0.008946`; new `+0.026052 / +0.026000 / +0.026297`.
Both published numbers reproduce exactly from their own cached three-seed tables. **On the new
emitter the target-sampling se is 167× the seed sem.**

**This is the whole mandate's premise, and it is falsified.** The mandate's arithmetic —
*"189's effect was 0.96σ on its instrument; on an instrument 77× tighter it would be ~80σ"* —
treats `±0.00424` as emission noise. It is not: 189's own `harness189.paired_bootstrap` resamples
**targets**, and 172 §4.1 established that axis. 85 perturbations are 85 perturbations under any
emitter.

---

## 1. What the new emitter *does* buy, and why it changes nothing here

Two of the three things an ordering key could gain from `rc_callset` turn out to be zero.

**(a) The counterfactual instrument is unchanged — bit-identically.** `Surface.reach(key)` reads the
reference's DE pool, the adjudicable mask, and `sign(pred_lfc)`. It never reads `pred_p`. The two
tables agree on `sign(pred_lfc)` for **98.72%** of pool pairs (832 of 65,044 differ), and the
differing pairs sit past `k*`, so the per-target `reach` vectors are **equal element-by-element**:

| key | old table | new table | `np.array_equal` | max \|diff\| |
|---|---:|---:|---|---:|
| `\|H\|` | 0.235334588 | 0.235334588 | **True** | 0 |
| 189 `λ = 0.50` | 0.239412972 | 0.239412972 | **True** | 0 |

**(b) What is genuinely new is DELIVERY.** Under a15's `quasi` emitter a key was a counterfactual;
under `rc_callset` the delivered order **is** the key (`rc_emit.py`:916 — the prescribed head is the
top 400 by `|key|`; `rc_emit.py`:922–926 — `|key|` then orders inside the block and inside its
complement; `rc_emit.py`:983 — delivered = prescribed ∪ unsuppressible). So a key can
be shipped and read back. `emit220.py` does exactly that, through the shipped `rc_hybrid_block`, and
its `absH` arm reproduces `de202.rc_callset.npz` **bit for bit** on `pred_p`, `pred_padj`,
`pred_lfc` and `pred_z` (85 × 10,779, max |diff| 0) with the per-gene panel-total digest identical
on 85/85 blocks.

**Delivery is lossy and the loss is measured:** delivered `reach` is **0.22334** against the
counterfactual `|H|`'s **0.23533**, i.e. the call-set block structure costs **−0.0120 raw** before
any key is changed.

**Delivery also moves the CALL SET, which no counterfactual screen could see.** The prescribed head
is selected on `|key|`, so a key change re-draws 400 genes per target:

| key | head-set Jaccard vs `\|H\|` |
|---|---:|
| `k189_l025` | 0.9205 |
| **`k189_l050`** | **0.8038** |
| `k189_l100` | 0.5285 |
| **`k190_km25`** | **0.9070** |
| `k190_p100` | 0.6166 |
| shuffled θ (12 draws) | 0.728 – 0.754 |

Delivering 189's key therefore replaces **20% of the call set**, which prices into `fid` and `jac`.
§5 reports what that costs.

---

## 2. 189 and 190, re-measured

### 2.1 The counterfactual, on both instruments — identical (§1a)

| key | reach | Δ vs `\|H\|` | target se | 95% CI | P(better) |
|---|---:|---:|---:|---|---:|
| `\|H\|` — the incumbent | 0.23533 | — | — | — | — |
| **189 `λ = 0.50`** | **0.23941** | **+0.00408** | 0.00422 | [−0.00378, +0.01272] | 0.838 |
| 189 `λ = 0.25` | 0.23691 | +0.00158 | — | — | — |
| **190 `ct-prop-k25`, λ = −0.25** | **0.23670** | **+0.00137** | 0.00169 | [−0.00137, +0.00516] | 0.784 |

### 2.2 Honest selection — the thing neither experiment did

189 selected λ on its evaluation surface and flagged it; 190 took a max over 13 λ × 8 networks.

| | in-sample max | **split-half held out** (400 splits) | **leave-one-target-out** |
|---|---:|---:|---:|
| **189** | +0.00408 at λ = +0.50 | **+0.00080 ± 0.00022** | +0.00408 (λ = 0.50 in 82/82 folds) |
| **190** | +0.00137 at λ = −0.25 | **−0.00200 ± 0.00020** | **−0.00015**, se 0.00077, P(>0) 0.463 |

The two cross-validations answer different questions and both matter. **LOO says 189's λ is stable**
— dropping one of 82 targets never moves the argmax — so its `+0.00408` is not a λ-selection
artefact. **Split-half says it does not survive halving the fitting set** (5.1× shrinkage), which is
the relevant regime for a 300-target panel we have never seen. **For 190 both agree and both are
negative**: its λ flips in the one fold that drops its single load-bearing target, and the held-out
value goes to zero or below.

### 2.3 The jackknife — how many perturbations carry the effect

| key | full | drop the single best target | drop the best three | targets that move at all |
|---|---:|---:|---:|---:|
| **189 `λ = 0.50`** | +0.00408 | **+0.00166** | **−0.00191** | 37 of 82 (22 up, 15 down) |
| **190 `λ = −0.25`** | +0.00137 | **−0.00016** | −0.00058 | 25 of 82 (14 up, 11 down) |

**60% of 189's effect and 111% of 190's is one perturbation.** The largest single-target
contributions for 189 are `+0.200 / +0.160 / +0.125` against losses of `−0.136 / −0.091 / −0.050`;
for 190 a single `+0.125` against `−0.054`. This is a target-count problem, and no emitter fixes it.

### 2.4 Permutation nulls at 400 seeds, with the selection inside the null

| key | real | null mean | null sd | z | P(null ≥ real) |
|---|---:|---:|---:|---:|---:|
| 189, θ permuted across the gene axis, λ fixed at 0.50 | +0.00408 | −0.00333 | 0.00379 | **+1.95** | **0.0225** |
| 189, with **λ re-selected inside the null** (120 draws) | +0.00408 | **+0.00133** | 0.00151 | — | **0.0667** |

189's key is **not separable from its own shuffle at α = 0.05 once the selection is inside the
null**, and a third of its raw value (`+0.00133` of `+0.00408`) is what any random gene-level tilt
buys through the max over eight λ. This is the same shape as 190's finding — where **eleven of
eighteen shuffle cells read higher than the real network** — and 189's own §4, where the call-set
gain at `n = 50` was `+0.0088` real against `+0.0091` shuffled.

### 2.5 190's `fid` row, the one thing it left un-null-tested

`fid220.py`. Precision at fixed `n`, macro over targets, on the `rc_callset` instrument, against 40
sign-shuffled networks with the best cell re-selected inside the null:

| `n` | `\|H\|` | 190 real | **null mean ± sd** | z | P(null ≥ real) |
|---:|---:|---:|---:|---:|---:|
| 10 | 0.5906 | +0.00000 | −0.00200 ± 0.00672 | +0.30 | 0.400 |
| 25 | 0.5412 | +0.00518 | +0.00616 ± 0.00397 | −0.25 | 0.625 |
| **50** | 0.5254 | **+0.00800** | **+0.00997 ± 0.00199** | **−0.99** | **0.850** |
| 100 | 0.5226 | +0.00071 | +0.00197 ± 0.00119 | −1.06 | 0.825 |
| 200 | 0.5181 | −0.00159 | −0.00036 ± 0.00072 | −1.71 | 0.950 |
| 400 | 0.5191 | −0.00003 | +0.00035 ± 0.00040 | −0.95 | 0.850 |

Best cell over the six sizes, re-selected inside the null: real **+0.00800**, null
**+0.01019 ± 0.00238**, P(null ≥ real) **0.850**. **190's gap is closed and it closes negative.**
(189's own key is *worse* than `|H|` at the front of this statistic: −0.0306 at `n = 10`.)

---

## 3. The structural claim — tested, and CORRECTED

`struct220.py`, `curve220.py`, `honest220.py`.

### 3.1 The prefix anatomy under `rc_callset` — the number 189 inherited from 172

| surface | reach | **median `k*`** | mean | p75 | p90 | `k* = 0` | targets carrying 94.4% |
|---|---:|---:|---:|---:|---:|---:|---:|
| a15, delivered (production order) | 0.18283 | **1** | 2.24 | 3 | 4 | 0.268 | 21 |
| a15, counterfactual `\|H\|` | 0.23533 | **3** | 7.06 | 6 | 10 | 0.122 | 31 |
| **`rc_callset`, delivered** | **0.22334** | **2** | 6.44 | 5 | 8 | 0.146 | **28** |
| `rc_callset`, counterfactual `\|H\|` | 0.23533 | **3** | 7.06 | 6 | 10 | 0.122 | 31 |
| `rc_callset`, ORACLE | 0.64529 | **67** | 470 | 640 | 1546 | 0.000 | 75 |

**172's `median k* = 3` reproduces exactly for `|H|` and it is 2 for the order `rc_callset` actually
delivers.** 172's "94.4% of `reach` in 32 of 85 targets" reproduces at **31 of 82 scored**.

### 3.2 The 0.90 bar is a per-perturbation ceiling function, and it is a COUNT

To have `k* ≥ d` a perturbation needs `ceil(0.9 d)` of its first `d` right — **5/5 at d = 5, 9/10 at
d = 10, 18/20 at d = 20**. Of the 82 scored targets, the number clearing it:

| key | d = 5 | d = 10 | d = 20 | d = 50 | median `k*` | reach |
|---|---:|---:|---:|---:|---:|---:|
| `\|H\|` (incumbent) | **24** | **12** | **5** | **2** | 3 | 0.23533 |
| delivered (production) | 25 | 8 | 5 | 2 | 2 | 0.22334 |
| 189 `λ = 0.50` | 21 | 9 | 3 | 1 | 3 | 0.23941 |
| 190 `λ = −0.25` | 23 | 14 | 4 | 2 | 2.5 | 0.23670 |
| **ORACLE** | **66** | **53** | **49** | **45** | 67 | 0.64529 |

**Neither key moves the count by more than two perturbations in either direction. The oracle moves
it by 40–60.** Mean front accuracy under `|H|` is 0.771 at d = 5, 0.693 at d = 10, 0.674 at d = 20;
the median shortfall — extra genes that would have to flip from wrong to right — is **1, 2 and 4**.

### 3.3 The transfer curve — and 189's arithmetic is wrong by 23×

An oracle-tilted key `|H| · exp(c · 1[correct])` interpolates continuously from `|H|` (c = 0) to the
oracle. Sweeping `c` measures the function 189's cap argument asserts (`rc_callset` instrument;
a15's is identical to 4 decimals):

| `c` | reach | Δ | acc@5 | acc@10 | acc@20 | median `k*` | P(`k*` ≥ 10) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | 0.23533 | +0.00000 | 0.7714 | 0.6931 | 0.6740 | 3.0 | 0.146 |
| 0.05 | 0.24069 | +0.00536 | 0.8086 | 0.7328 | 0.7140 | 3.0 | 0.183 |
| 0.10 | 0.24672 | +0.01138 | 0.8286 | 0.7672 | 0.7520 | 3.5 | 0.244 |
| 0.20 | 0.25507 | +0.01974 | 0.8800 | 0.8310 | 0.8170 | 5.5 | 0.415 |
| **0.30** | **0.27327** | **+0.03794** | 0.9057 | **0.8724** | 0.8630 | 10.0 | 0.524 |
| 0.50 | 0.31241 | +0.07708 | 0.9429 | 0.9086 | 0.9300 | 21.0 | 0.585 |
| 1.00 | 0.43568 | +0.20034 | 0.9657 | 0.9466 | 0.9750 | 44.5 | 0.610 |
| 6.00 | 0.64363 | +0.40829 | 0.9800 | 0.9638 | 0.9950 | 67.0 | 0.646 |

**Read the `c = 0.30` row.** Mean front accuracy at depth 10 is **0.872 — still below 0.90 — and
`reach` has already moved +0.0379.** `reach` is an average of a ceiling function over 82
perturbations, so raising accuracy moves the *distribution* across the bar long before the *mean*
crosses it.

**Local slope at the incumbent: `d reach / d acc@10 = +0.143`.** 189's own headline front gain
(+0.0820 inside `|H|` decile 1) would therefore be worth **+0.0117 raw `reach`** if it were a real
front gain — against 189's structural estimate of **≤ +0.0005**, i.e. **23× larger**.

**Be precise about what is wrong.** 189's *measurements* stand: reordering inside the top 2–10% of
`|H|` by θ really is worth `+0.00045 … +0.00141`, and this experiment reproduces that family
(§5.4's `k189_top10` / `k189_top25` rows deliver `+0.00062` and `+0.00155`). What fails is the
*mechanism* offered for them — *"lifting front accuracy 0.587 → 0.669 does not cross 0.90, so `k*`
does not move"* — and therefore the generalisation 190 §3.1 built on it: *"a prior confined to a
small, front-adjacent slice of pairs cannot move this metric no matter how good it is."* **It can.**
The `c = 0.05` row moves `reach` by `+0.0054` — seven times the whole 189 family — while lifting
acc@10 only from 0.693 to 0.733, nowhere near 0.90, and it touches nothing but the order. The bar
binds per perturbation, not on the average, so a small accuracy gain buys `reach` continuously.
**The real cap is the slope, +0.143, and the real constraint is §3.4.**

### 3.4 And the family dies anyway, because no key moves the front

| key | reach | acc@5 | **acc@10** | acc@20 | acc@50 |
|---|---:|---:|---:|---:|---:|
| `\|H\|` (incumbent) | 0.23533 | 0.7714 | **0.6931** | 0.6740 | 0.6298 |
| 189 `λ = 0.50` | 0.23941 | 0.7657 | **0.7034** (+0.0103) | 0.6750 | 0.6302 |
| 189 `λ = 0.25` | 0.23691 | 0.7743 | 0.6983 (+0.0052) | 0.6760 | 0.6353 |
| 190 `λ = −0.25` | 0.23670 | 0.7657 | **0.6862** (−0.0069) | 0.6630 | 0.6332 |

Rule 13's profile-preserving gate is **+0.0019 overall = +0.00783 raw h1-test `reach`** (×1.290,
÷span 1.12927, ÷6). On the measured slope that is **+0.0547 of acc@10**. 189's key delivers
**+0.0103** — 5.3× short — and 190's is **negative**.

**189's `+0.082` was never a front gain.** It is the accuracy gap between θ terciles *inside* `|H|`
decile 1 — a conditional contrast on a stratum, not a movement of the front once θ is allowed to
reorder that stratum against the rest of the pool.

**And there is a mechanism for why it cannot be, which nobody had measured.** θ is **2.4× flatter
exactly where the key reads**: over the 9,128 of 18,080 genes that enter *some* target's top-400
`|H|` head, `sd(θ) = 0.1397`, against **0.3399** outside and 0.2602 overall. A gene-level tilt whose
dynamic range collapses on the front cannot reorder the front. It shows up directly in the delivered
call set — the *real* θ at λ = 0.50 moves the head-set Jaccard to **0.804** while a **shuffled** θ
with the identical marginal distribution moves it to **0.728–0.754**, i.e. the real prior perturbs
the head *less* than its own permutation does (`loop220.py` §5).

### 3.5 The honest ceiling of the whole gene-level-key family: **negative**

189's class of key is `|H| × u_g` for a per-gene scalar `u`. Its best possible `u` is the gene's own
sign accuracy. Estimated on the same targets it is scored on, that is an oracle; estimated on a
disjoint half of the targets, it is the honest ceiling (12 random halves, tilt swept 1–8):

| surface | split-half ρ(per-gene accuracy) | IN-FOLD gene oracle | **CROSS-FOLD (honest)** |
|---|---:|---:|---:|
| a15 | **+0.0101** | +0.09683 | **−0.01047 ± 0.00160** |
| `rc_callset` | **+0.0108** | +0.09485 | **−0.01143 ± 0.00166** |

Cross-fold, by tilt: `1: −0.0105 · 2: −0.0197 · 4: −0.0361 · 8: −0.0510` — monotone downward, so
there is no tilt strength at which the family is positive.

**This is the decisive structural statement, and it is stronger than 189's own.** 189 measured
`split-half ρ = +0.0097` on per-gene accuracy and converted it into a `q_w`/`nmae` ceiling of
`+0.0019 overall`. Converted into `reach` — the member the family was actually proposed for — **the
same fact makes the ceiling negative**: a per-gene weight fitted on 41 targets and applied to 41
others *loses* `reach`. The in-fold number is `+0.0968`, **an inflation of ∞ in sign and ~9× in
magnitude over the `q_w` version's 22.3×.**

---

## 4. The loophole the mandate named — real in direction, closed by learnability

`loop220.py`. `reach` averages `k*/n_real` over perturbations and 94.4% of it lives in 31 of 82
targets, so a key that pays only where the front is nearly pure could pay even if the median is
untouched. Three questions, in order.

**(1) Is the gain where `reach` lives? YES.**

| key | Δ | of which inside the 31 top-`reach` targets |
|---|---:|---:|
| **189 `λ = 0.50`** | +0.00408 | **+0.00651 (160%)** |
| 189 `λ = 0.25` | +0.00158 | +0.00317 (201%) |
| 190 `λ = −0.25` | +0.00137 | +0.00000 (0%) |

189's gain is *more* than all of its own total inside those targets — the complement is net
negative. **The mandate's hypothesis is confirmed in direction.**

**(2) What is the ceiling of any gating rule?** Apply the key only where it helps (a target-level
oracle, and therefore an upper bound on every gate):

| key | plain | **gate = oracle** | → overall | vs the +0.0019 gate |
|---|---:|---:|---:|---|
| **189 `λ = 0.50`** | +0.00408 | **+0.00843** | **+0.00205** | **clears it — by 8%** |
| 189 `λ = 0.25` | +0.00158 | +0.00335 | +0.00081 | no |
| 190 `λ = −0.25` | +0.00137 | +0.00243 | +0.00059 | no |

**So the loophole is not closed by the structural cap.** A perfect per-target gate on 189's key
would clear rule 13. That is a real and previously untested finding and it is why this experiment
was worth running.

**(3) Is such a gate learnable from anything we can see? NO.** Every stratifier below is computed
from our own `H` and our own call set; none touches the reference. Correlation with the per-target
Δ`reach`:

| stratifier | 189 `λ 0.50` | 189 `λ 0.25` | 190 `λ −0.25` |
|---|---:|---:|---:|
| `\|H\|` row norm | −0.035 | −0.037 | −0.069 |
| `\|H\|` top-10 mean | −0.100 | −0.019 | −0.099 |
| `\|H\|` head concentration | −0.091 | +0.016 | −0.061 |
| `\|H\|` top10 / top100 | −0.054 | +0.032 | −0.036 |
| our own `n_pred` | −0.022 | −0.020 | −0.014 |
| θ of our own top 10 | −0.018 | +0.019 | −0.177 |
| `\|H\|` gap at rank 10 | −0.096 | −0.042 | −0.048 |

**Not one exceeds `|r| = 0.10` for 189.** Selecting the best (stratifier, direction, threshold)
gate on half the targets and scoring on the other half, 400 splits:

| key | in-sample best gate | **held out** | P(> 0) |
|---|---:|---:|---:|
| **189 `λ = 0.50`** | +0.00507 (`\|H\|` head concentration, low, top 75%) | **+0.00044 ± 0.00016** | 0.545 |
| 189 `λ = 0.25` | +0.00294 | −0.00029 ± 0.00009 | 0.340 |
| 190 `λ = −0.25` | +0.00208 | +0.00062 ± 0.00008 | 0.578 |

**+0.00044 raw = +0.00011 overall, 17× below the gate.**

**And the mechanism for why the loophole cannot be widened:**
`corr(per-target headroom, per-target incumbent reach) = **−0.900**`. The oracle's headroom is
`+0.40995` and only **13%** of it sits in the 31 targets carrying 94.4% of current `reach`; 94.4% of
the *headroom* is spread over **58** targets. **The perturbations where a key can help are, by
construction, not the ones already scoring** — so concentrating on the high-`reach` targets buys the
small slice of headroom that remains there.

---

## 5. Phase 2 — the six members, DELIVERED, and rule 13's gate

`emit220.py` delivers **14 keys as 20 arms** (3,850 s, peak RSS 3.53 GiB) through the shipped
`rc_hybrid_block` on h1-test, with three matched seeds for the three load-bearing arms; a second
run adds **12 more shuffle keys as 13 arms** (2,745 s, 3.63 GiB). `six220.py` prices them.

### 5.1 The class test — measured on all 85 blocks, not asserted

| test | result |
|---|---|
| arm `absH` ≡ `de202.rc_callset.npz` on `pred_p`, `pred_padj`, `pred_lfc`, `pred_z` | **bit-identical, 85 × 10,779, max \|diff\| 0** |
| …its per-gene panel-total digest | **85/85** |
| …scored beside 202's own table | `n_pred` 2634, `reach` 0.223341, `fid` 0.494362, `jac` 0.072975, `nmae` 0.966702 — **equal to 1e−12** |
| **every one of the 33 delivered arms vs `absH`, per-gene panel-total digest** | **85/85 identical — ALL PROFILE-PRESERVING** |
| measured `nmae` movement, every arm | **\|Δ\| ≤ 4 × 10⁻⁷** |

**So rule 13's class is `emission-law, profile-PRESERVING`: `reach` toll = 0, `pds` and `mse`
pinned at exactly 0.000000** (board-verified to sixteen significant figures on `rc_hybrid` and
again on `rc_callset`), and `nmae` measured at zero to seven decimals rather than assumed. Four of
six members are pinned by construction and the remaining two are `reach` and the `fid`/`jac` pair.

### 5.2 The delivered effect — and DELIVERY ATTENUATES IT 4.6×

Δ`reach` against `absH`, three matched seeds (`native` = the reference's own BH pool and our own
BH call set):

Each seed is paired against **its own** `absH` seed — 214's convention, and the estimand that
removes the shared emission stream:

| arm | seed 1 | seed 2 | seed 3 | **mean** | **seed sem** | counterfactual |
|---|---:|---:|---:|---:|---:|---:|
| `absH`, each seed against seed 1 (the incumbent's own spread) | +0.00000 | −0.00006 | +0.00006 | — | sd 0.00006 | — |
| **189 `λ = 0.50`** | +0.000766 | +0.000925 | +0.000996 | **+0.000896** | **0.000068** | +0.00408 |
| **190 `λ = −0.25`** | +0.001070 | +0.001074 | +0.000862 | **+0.001002** | **0.000070** | +0.00137 |

**189's counterfactual `+0.00408` delivers as `+0.00090` — an attenuation of 4.6×** — because
`rc_callset` promotes its entire 2,600-gene prescribed call set as one block before the key orders
anything, so the key only re-ranks *inside* that block. (The same mechanism costs the incumbent
`0.23533 → 0.22334` before any key changes.)

**And the two axes, on the delivered arm itself — this is the experiment in one row:**

| delivered arm | Δ`reach` | **seed sem** | **target-bootstrap se** | ratio | 95% CI over targets | P(better) |
|---|---:|---:|---:|---:|---|---:|
| **189 `λ = 0.50`** | +0.000766 | **0.000068** | **0.00392** | **58×** | [−0.00672, +0.00872] | **0.576** |
| **190 `λ = −0.25`** | +0.001070 | 0.000070 | 0.00183 | 26× | [−0.00216, +0.00502] | 0.706 |

214's instrument reproduces exactly and resolves the effect to **13σ on the seed axis**. On the axis
that decides whether it generalises to 300 unseen perturbations it is **0.20σ**, and `P(better)` is
**0.576** — a coin flip.

### 5.3 On the significance-selected pool the sign flips

`wald89` is 188's pool — the reference's top-`n_real` genes by `|lfc|·√cpm`, median `n_real` 89,
the regime closest to the board's (`n_real` ≈ 90–110 against h1-test's BH `n_real` = 222):

| arm | seed 1 | seed 2 | seed 3 | **mean** | seed sem |
|---|---:|---:|---:|---:|---:|
| **189 `λ = 0.50`** | −0.007400 | −0.008726 | −0.007213 | **−0.007780** | 0.000476 |
| **190 `λ = −0.25`** | −0.004798 | −0.006511 | −0.004649 | **−0.005319** | 0.000597 |

**Both keys LOSE `reach` on the pool whose `n_real` matches the board's, at every one of the
three seeds.** Rule 9's
warning applies to the level, not the delta, but a delta that changes sign between two pools of the
same reference is a fact about the key and it is reported rather than averaged away.

### 5.4 The dose ladder and the conditional forms, delivered (seed 1)

| arm | head Jaccard vs `\|H\|` | Δ`reach` native | Δ`reach` wald89 |
|---|---:|---:|---:|
| 189 `λ = 0.25` | 0.9205 | **−0.00219** | −0.00762 |
| **189 `λ = 0.50`** | 0.8038 | **+0.00077** | −0.00740 |
| 189 `λ = 1.00` | 0.5285 | **−0.00432** | −0.01359 |
| 189 top-10% of `\|H\|`, λ = 1.0 (189 §6.2's form) | 0.7188 | +0.00062 | −0.00649 |
| 189 top-25% of `\|H\|`, λ = 0.5 | 0.8340 | +0.00155 | −0.00677 |
| **190 `λ = −0.25`** | 0.9070 | **+0.00107** | −0.00480 |
| 190 `λ = +1.00` | 0.6166 | **−0.00913** | +0.00111 |

**The ladder is not monotone in either direction and it is not monotone in the same direction on the
two pools.** λ = 0.25 is negative while λ = 0.50 is positive and λ = 1.00 is negative again — the
signature of a quantity that is being read at the level of its own noise.

### 5.5 The six-member price, marginal to `rc_callset` / `cd4f034`

`reach` × **1.290** (185's board point for an order-changing arm); `fid` as precision; `jac` through
208's S\* convention with board `NP` pinned at 2,599 and `NR` = **221.8** solved so the mean-field
Jaccard reproduces `rc_callset`'s own board raw `jac` 0.04043; `nmae` raw / span.

| arm | pool | `pds` | `mse` | `reach` | `fid` | `jac` | `nmae` | **SUM** | gate | verdict |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| **189 `λ = 0.50`** | native | 0 | 0 | +0.000186 | −0.000204 | +0.000264 | −0.0000000 | **+0.000245** | +0.0019 | **VETO** |
| **189 `λ = 0.50`** | wald89 | 0 | 0 | −0.001797 | −0.000204 | +0.000214 | −0.0000000 | **−0.001786** | +0.0019 | **VETO** |
| **190 `λ = −0.25`** | native | 0 | 0 | +0.000260 | −0.000012 | +0.000069 | −0.0000000 | **+0.000317** | +0.0019 | **VETO** |
| **190 `λ = −0.25`** | wald89 | 0 | 0 | −0.001165 | −0.000012 | +0.000098 | −0.0000000 | **−0.001080** | +0.0019 | **VETO** |
| 189 `λ = 1.00` | native | 0 | 0 | −0.001050 | −0.000305 | −0.000239 | −0.0000000 | −0.001594 | +0.0019 | VETO |
| 190 `λ = +1.00` | native | 0 | 0 | −0.002218 | −0.000103 | −0.000181 | +0.0000000 | −0.002503 | +0.0019 | VETO |

**The best delivered arm in the whole experiment is `+0.000317` against a `+0.0019` gate — 6× short
— and it is 190's, whose own shuffled control beats it (§6).** At 164's stale 0.4277 transfer factor
every row is a third of what is shown. **Nothing ships. No slot is spent. Nothing was built.**

Note the `fid`/`jac` trade the counterfactual screens could not see: delivering 189's key moves
recall **+0.0068** (jac +0.000264 board) but costs precision (`fid` −0.000204 board), because 20% of
the 400-gene head is re-drawn. The two nearly cancel.

---

## 6. Phase 3 — the falsifiers

**F1 — the permutation control at the new precision. IT FIRES IN ALL FOUR CELLS, and the mandate
pre-committed to what that means:** *"run them, and if a control still ties or beats the arm, the
family is dead for good and that is the result."*

`null220.py`. Twelve θ-shuffles and six regulon-sign shuffles, each **delivered** through the shipped
`rc_callset` emitter. Every shuffle is a single seed, so every row is measured against `absH`
**seed 1** — which is why the real arms' means differ in the last digit from §5.2's paired values.
Known-answer row: the second emission run's `absH` is byte-identical to the first's on `pred_p`,
`pred_padj`, `pred_lfc`, `pred_z` and `digest`, and both reproduce `de202.rc_callset.npz` bit for
bit.

| pool | arm | **real** (3 seeds, sem) | **shuffled null** (n draws) | shuffles beating the arm | z | **P(shuffle ≥ real)** |
|---|---|---:|---:|---:|---:|---:|
| native | **189 `λ = 0.50`** | **+0.00090** (0.00008) | −0.00121 ± 0.00282 (12) | **2 of 12** | +0.75 | **0.167** |
| native | **190 `λ = −0.25`** | **+0.00100** (0.00004) | **+0.00141 ± 0.00033** (6) | **5 of 6** | **−1.24** | **0.833** |
| wald89 | 189 `λ = 0.50` | −0.00729 (0.00008) | −0.00357 ± 0.00380 (12) | **9 of 12** | −0.98 | **0.750** |
| wald89 | 190 `λ = −0.25` | −0.00483 (0.00005) | −0.00447 ± 0.00030 (6) | **6 of 6** | −1.22 | **1.000** |

**190's real network is systematically *worse* than its own sign-shuffle on both pools.** The
sign-shuffle preserves topology and head-set Jaccard (0.907 real against 0.907–0.909 shuffled), so
its null is tight — sd **0.00033** — and the real network sits **1.2σ below it**, beaten by 5 of 6
draws on `native` and 6 of 6 on `wald89`. That is 190's own 2026-09-03 finding, reproduced on the
DELIVERED instrument at ~30× tighter null precision instead of on a counterfactual one. **189 is not
separable from its own shuffle at α = 0.05 on either pool.**

*The delivered θ-shuffle null has sd **0.00282** — 41× the 0.000068 paired seed sem and of the same
order as the 0.00422 target-bootstrap se. That is the third independent demonstration that the axis
limiting this family is not the one the new emitter tightened.*

**F2 — the old-emitter comparison reproduces.** §0: `+0.00408 / se 0.00422 / P 0.838` against 189's
`+0.00408 / 0.00424 / 0.828`; `+0.00137 / CI [−0.00137, +0.00516]` against 190's `+0.00137 /
[−0.00140, +0.00520]`; both emitters' three-seed sems (`0.00695`, `0.00009`) reproduce exactly from
their own cached tables. The two instruments are measuring the same thing — in fact the *same
object* (§1a).

**F3 — the key is not selected on the surface it is scored on.** λ is chosen by split-half over
targets (400 splits) and by leave-one-target-out (§2.2); the target gate is chosen by split-half
(§4); the permutation nulls carry the selection inside them (§2.4). Every honest number is smaller
than its in-sample twin, and 190's changes sign.

**F4 — the arm class is proved by digest, not asserted.** 85/85 blocks identical on all 33 delivered
arms, `nmae` |Δ| ≤ 4e−7, `absH` bit-identical to the shipped `rc_callset` table (§5.1).

**F5 — no donor table was opened and no `n_real` was rescaled.** The effect matrix is 169's recovered
`H` for the frozen a15 artifact; `data/donors/h1_2025` (85/85 coverage, pure oracle leakage) is never
touched, and `Surface(n_real_target=…)` is never called anywhere in this directory
(`inputs220.py`).

---

## 7. What would reopen this — the specification, in a currency any future head can be measured in

`noisy220.py`. A per-(target, gene) confidence bit `o` that equals `1[our sign is right]` with
probability `a`, keyed as `|H| · exp(c · o)`, 20 seeds per cell:

| `a` | corr(`o`, truth) in pool | c = 0.2 | c = 0.5 | c = 1.0 | c = 3.0 |
|---:|---:|---:|---:|---:|---:|
| 0.50 | −0.000 | −0.00350 | −0.00687 | −0.01281 | −0.03070 |
| 0.55 | +0.099 | −0.00138 | −0.00254 | −0.00646 | −0.01968 |
| 0.60 | +0.200 | +0.00069 | +0.00230 | +0.00183 | −0.00860 |
| **0.65** | **+0.299** | +0.00270 | +0.00689 | **+0.00944** | +0.00520 |
| 0.70 | +0.399 | +0.00463 | +0.01217 | +0.01886 | +0.02220 |
| 0.80 | +0.599 | +0.00926 | +0.02680 | +0.05027 | +0.07682 |
| 0.90 | +0.800 | +0.01422 | +0.04793 | +0.12384 | +0.29988 |
| 1.00 | +1.000 | +0.01974 | +0.07708 | +0.20034 | +0.37916 |

**The gate (+0.00783 raw) is cleared at `a ≈ 0.639`.** In the two currencies any future confidence
head can be measured in — Pearson correlation with `1[our sign is right]` inside the reference DE
pool, and the Mann-Whitney AUC for the same event:

| signal | corr | **AUC** |
|---|---:|---:|
| **`\|H\|` — the incumbent key, the best ever measured here** | +0.0743 | **0.5350** |
| 189's θ | +0.0021 | 0.5107 |
| control expression | +0.0080 | 0.5049 |
| **REQUIRED to clear rule 13's `+0.0019`** | **+0.2992** | **0.6387** |

**A shippable ordering key needs a per-pair confidence at AUC ≈ 0.64 where `|H|` reads 0.535 — and
it has to be orthogonal to `|H|` on top of that.** Nothing measured anywhere in this project is
within an order of magnitude of the requirement. That is the closure, stated as a target rather than
as a verdict.

**And the size of the prize is bounded from the outside.** Board raw `reach`: ours **0.24983**,
`Aginglab.com` **0.2621**, field best (`luka`/PertMind) **0.2908**. Matching the best team in the
world on `reach` is worth `0.0410 × 1.12927 / 6 = ` **+0.0077 overall**; matching `Aginglab` is worth
**+0.0023**. 215 measured us at **+0.70σ** on the board's own `reach`-given-`pds` conditional against
`Aginglab`'s +1.00σ — that 0.30σ is **+0.0022 overall**. **`cd4f034`, one direction arm, delivered
+0.0131.** The ordering channel is not where the remaining points are, and this experiment is the
measurement that says so from the inside.

---

## 8. What this changes for other experiments

* **To 219 (the closure audit), which handed 189/190 here: §8.4's ×4.5 does not apply to this
  family, and re-pricing them would have been the wrong correction.** §8.4's two errors are (1) a
  direction gain propagates to five members that were carried at zero, and (2) PLOLO understates the
  board's cosine response. An ordering arm has **no members carried at zero** — `pds`, `mse` and
  `nmae` are pinned by construction and measured at zero here, and `reach`, `fid` and `jac` were the
  only things 189/190 ever priced — and it **does not move the transfer cosine at all**. The
  multiplier has nothing to multiply. **189 and 190 stay closed, and now they are closed on the
  measurement rather than on the margin.**
* **214's headline should be read narrowly.** *"`reach` under `rc_callset` is a near-noiseless read
  of the ordering key… any key-changing arm is now 77× cheaper to resolve"* is true of the
  **emission-seed** axis and false of the **target-sampling** axis, which is 167× larger on the same
  arm and 1.08× unchanged between the two emitters. The 77× is real and useful — it is what makes a
  `+0.0009` delivered effect resolvable at all — but it does not convert a target-limited screen into
  a resolved one. **h1-test has 85 perturbations under every emitter.**
* **189's mechanism paragraph should be corrected in place.** *"Lifting front accuracy 0.587 → 0.669
  does not cross 0.90, so `k*` does not move"* is false: at `c = 0.30` the measured front accuracy is
  **0.872** and `reach` has already moved **+0.0379**, because `reach` averages a ceiling function
  over 82 perturbations and the mean does not have to cross the bar for the distribution to. 189's
  *measurements* (`+0.00045` for the top-2% reordering, `≤ +0.0005` for that family of keys) are
  correct; the *generalisation* to "a prior confined to a front-adjacent slice cannot move this
  metric no matter how good it is", which **190 §3.1 adopted as an independent derivation**, is not.
  The correct cap is **d `reach` / d acc@10 = +0.143**, and the reason the family dies is that no
  key moves acc@10 by more than +0.010 against the +0.055 the gate needs.
* **190's `fid` gap is closed and closes negative** (§2.5). Its `+0.00555` at `n = 50` reads
  `+0.00800` on this instrument against a sign-shuffled null of `+0.00997 ± 0.00199`.
* **The `reach` channel's whole remaining prize is +0.0077 overall and it is bounded from outside.**
  Board raw `reach`: ours **0.24983**, `Aginglab.com` **0.2621**, field best (`luka`/PertMind)
  **0.2908**. Becoming the best `reach` team in the world is worth **+0.0077**; matching `Aginglab`
  is **+0.0023**; closing 215's measured **+0.70σ → +1.00σ** residual gap is **+0.0022**. `cd4f034`,
  one direction arm, delivered **+0.0131**. **§8.2 is right that everything left is direction, and
  this is the measurement that says so from inside the ordering channel.**
* **For anyone proposing a per-gene weight anywhere:** 189's split-half `ρ = +0.0097` on per-gene
  sign accuracy converts to a **negative** `reach` ceiling (−0.0114 cross-fold against +0.0949
  in-fold), not merely a small positive one. The in-fold statistic is not a weak version of the
  honest one; it has the wrong sign.

---

## 9. What I did NOT measure, and why

| not measured | why | can it move the verdict? |
|---|---|---|
| a **challenge-side** read of either key | needs a truth we do not have; the board is the only instrument and rule 7 forbids spending a slot on a `+0.0003` arm | no — the gate is 6× above the best delivered arm |
| 190's other seven networks, delivered | `ct-prop-k25` is the one 190 selected and the only one above its own λ = 0 row; the other four selected λ = 0, i.e. "do not use the prior" | no |
| a **line-specific** signed network, or a per-pair confidence conditionally better where `\|H\|` is weak | 190 §7.3 named these as the only reopening conditions and neither exists in any data we hold | §7 gives the numeric bar such a signal would have to clear |
| the **direction** use of either prior | 190 §6 priced it (third independent negative) and it is a direction arm, not an ordering arm — 213 owns the transport law | no |
| `k562_essential` or any other fitting line for θ | 189 §3.3 measured that line count does not convert (donor-side +0.0868 → +0.1047 buys +0.0014 raw `reach`) | no |
| the unseen D/E/F panel | untestable by construction; θ carries no target identity and 189's F3/F4 showed deleting the overlapping targets moves it by ρ 0.999 | — |

---

## 10. Reproduction

```
uv run python inputs220.py    # rule 11: path + sha256 for every frozen input   (5 s)
uv run python ka220.py        # Phase 0: both instruments, both sems            (10 s)
uv run python struct220.py    # prefix anatomy, honest selection, family ceiling, nulls (63 s)
uv run python curve220.py     # the front-accuracy -> reach transfer curve      (35 s)
uv run python honest220.py    # leave-one-target-out lambda, the bar as a count, jackknife (30 s)
uv run python loop220.py      # the loophole: concentration, gate ceiling, learnability (20 s)
uv run python fid220.py       # 190's un-null-tested fid row, null-tested       (60 s)
uv run python noisy220.py     # the confidence-signal specification (sec 7)     (25 s)
uv run python emit220.py      # DELIVER 20 arms through the shipped rc_callset  (3850 s)
ARMS=absH,shuf189_5,...,shuf190_6 OUTTAG=de220x uv run python emit220.py   # 13 more (2745 s)
uv run python six220.py       # the six members and rule 13's gate              (95 s)
uv run python null220.py      # F1: 12 + 6 DELIVERED permutation controls       (60 s)
```

`/usr/bin/time -l` maximum resident set size: `emit220` **3.53 / 3.63 GiB** (the two runs, the
heaviest by far), `ka220` 0.48 GiB, everything else <= 0.42 GiB — all far under the 10 GB budget.
Total wall clock ~ **1 h 55 min**, of which 1 h 50 min is the two emission runs.

### Artefacts

| file | what it produces |
|---|---|
| `common220.py` | both instruments, `reach_vec` / `kstar_vec` / `front_acc`, 189's θ, 190's `R`, the paired bootstrap |
| `inputs220.py` | rule 11 — path + sha256 for 15 frozen inputs → `work/inputs220.json` |
| `ka220.py` | **§0** — every known-answer row on both instruments, and both emitters' sems |
| `struct220.py` | **§2.2, §2.4, §3.1, §3.5** — prefix anatomy, honest λ, family ceiling, 400-seed nulls |
| `curve220.py` | **§1a, §3.3, §3.4** — the sign-identity check and the front-accuracy → `reach` curve |
| `honest220.py` | **§2.2, §2.3, §3.2** — leave-one-target-out λ, the bar as a count, the jackknife |
| `loop220.py` | **§4, §3.4** — the loophole, and θ's dynamic range on the front |
| `fid220.py` | **§2.5** — 190's `fid` row, null-tested |
| `noisy220.py` | **§7** — the noisy-oracle specification, in corr and AUC |
| `emit220.py` | **§5** — the keys delivered through the shipped `rc_hybrid_block`; the KA and the digest class test |
| `six220.py` | **§5.5** — the six members, both pools, three seeds, rule 13's gate |
| `null220.py` | **§6 F1** — the delivered permutation controls |

### Integrity

Every script calls `harness172.check` (or the `absH` ≡ `rc_callset` bit-identity) before it reports
anything and aborts on failure; all rows passed on every run. No sealed data opened, no donor table
opened, no `cell-eval2`, no build, no `vcc` invocation, no `submissions/` write, and no file outside
`experiments/220-ordering-keys-retry/` created or modified. `src/` and `scripts/` were read only.

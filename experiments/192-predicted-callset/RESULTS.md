# 192 — the call set. **188's key does not survive the substitution of our `lfc` for the reference's — and the call set is still worth +0.009, from a gene-level prior plus a protected `|H|` head.**

**Date:** 2026-09-03 · measurement only. No `cell-eval2`, no `vcc`, no build, no submission,
nothing written to `submissions/`, no file outside `experiments/192-predicted-callset/` modified.
Peak RSS **2.0 GB**, one process at a time.

> **This write-up uses the two constants the coordinator broadcast at 11:45 UTC** (`rc_hybrid`
> landed: `reach` +0.054940 scaled, overall +0.0093325). Rule 13's −0.0146 toll is **refuted for
> profile-preserving arms**, and the h1-test → board `reach` transfer factor for the designed-order
> class is **1.290**, not 0.4277. Every `reach` number below is given in **raw h1-test units** as
> well as converted, at **both** factors, so the conversion can be redone.

---

## Decision, in one line

**The mandate's central hypothesis is FALSIFIED and the mandate's objective is nevertheless MET.**
188's `|lfc|·√cpm` reconstructs the reference's DE set at 0.7219 — but that is an **oracle** number
built from the *reference's own* `lfc`, and substituting *our* `lfc` makes it **worse than the
control CPM alone** (fold 1.88 vs 2.05). Our predicted effect contributes **exactly zero** to
membership: a fitted logistic on (cpm, detection rate, cv²) and the same logistic with
`|lfc_pred|` **and** `|H|` added both read **2.37**.

What does work is a two-part call set that follows directly from the enrichment curve — a
**gene-level detectability prior** for membership (`jac`) with the **`|H|` head protected inside
it** for order (`reach`):

```
call set per (context, perturbation) =  top 400 by |H|   ∪   top 2,200 by control CPM
                                        (n_pred 2,600, against a15's 4,474)
```

**+0.0091 to +0.0111 overall**, marginal to `rc_hybrid`, positive at every `n_real` level tested
and at both `reach` transfer factors, moving **all three** DE members the right way at once
(`fid` +0.0010, `jac` +0.0129, `reach` +0.0112 raw h1-test at the retrodicted level). It preserves
the delivered per-gene panel totals by construction, so `pds`, `mse` and `nmae` are pinned.

**Recommendation: GO, as a rider on the next `rc_hybrid`-family build. It needs no slot of its own.**

Two by-products are worth as much as the arm:

1. **The board's own `n_real` distribution is now measured, twice, from two different members**
   (§4.2) — the quantity rule 9 says our harnesses cannot read and that has cost four board
   predictions. Median **31–45**, mean **270–383**; **0.35–0.5×** h1-test's.
2. **`n_real_p` is NOT predictable** (§3.3): LOTO Pearson **+0.24**, 6% of the variance, and on
   INBOX-1's own strength-bin criterion it calls a truly *negligible* perturbation *strong* **22%**
   of the time. Oracle per-perturbation sizing is worth **+0.036**; the buildable version is
   **−0.097**. **INBOX-1's arm is closed by its input, not by its payoff.**

---

## 0. Phase 0 — known-answer rows. **All pass, most exactly.** (`ka192.py`, `work/ka192.log`)

### 0.1 What "reproduces at 0.7219" actually is

188's own definition, read off `pool188.py:99-106` and re-implemented independently here:

```
for each perturbation p with n_real_p = |{g : padj_p,g < 0.05}| > 0:
    T_p(key) = the top n_real_p genes of row p by `key`
    o_p      = |T_p ∩ SIG_p| / n_real_p
overlap(key) = mean_p o_p            (unweighted over the 85 perturbations)
```

It is **not** a rank correlation and **not** an AUC. Because `|T_p| = |SIG_p| = n_real_p` exactly,
`o_p` is simultaneously the recall and the precision of the key at that row's own true size.
**Both the `lfc` and the `cpm` in 188's keys are the REFERENCE's**, so 0.7219 bounds the family; it
is not a buildable score. The whole of §1 is about the gap between that bound and what we can build.

| gate | expected | measured | verdict |
| --- | ---: | ---: | :-: |
| `overlap(|lfc_real|·√cpm)` | 0.7219 | **0.7219** | PASS |
| `overlap(|h_real|)` | 0.6406 | **0.6406** | PASS |
| `overlap(|lfc_real|)` | 0.3796 | **0.3796** | PASS |
| `overlap(cpm)` | 0.1512 | **0.1512** | PASS |
| median / mean `n_real` | 89 / 765 | **89 / 765.2** | PASS |
| a15 h1-test `reach` | 0.18283 | **0.18283** | PASS |
| a15 median `n_pred`, raw BH / debench | 6035 / 6034 | **6035 / 6034** | PASS |
| a15 `fid` | 0.50701 | **0.50701** | PASS |
| a15 `jac` | 0.06824 | **0.06824** | PASS |
| a15 `nmae` | 0.96659 | **0.96659** | PASS |
| `|H|` key `reach` | 0.23533 | **0.23533** | PASS |

**The random baseline for 188's statistic, which 188 did not report:** `mean_p n_real_p/G` =
**0.0736** analytically, 0.0737 measured over 5 shuffles. So `wald` is **9.8×** random and `cpm`
alone is **2.05×**. Every fold-enrichment below is against that denominator.

Donors: `k562_gwps` + `orion_hek293t` only, through 095's frozen h1-test artifact.
`data/donors/h1_2025` — the 85/85 same-line oracle — never enters.

---

## 1. Phase 1 — the enrichment curve. **This is the result the project has never had.**

`enrich192.py`, `work/enrich192.log`. 85 perturbations × 10,779 gated genes; 3 rows with
`n_real = 0` are dropped from every row-average and stated rather than hidden.

**Fold-enrichment `p(n) / (mean_p n_real_p / G)`:**

| key | n=25 | 50 | 89 | 150 | 300 | 600 | 1200 | 2400 | 4800 | **n = n_real_p** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **BUILDABLE** | | | | | | | | | | |
| `cpm` alone (challenge's own NTCs) | 4.18 | 3.73 | **3.35** | 3.01 | 2.64 | 2.30 | 1.95 | 1.64 | 1.33 | **2.05** |
| control detection rate alone | 3.71 | 3.55 | 3.24 | 2.93 | 2.59 | 2.23 | 1.92 | 1.61 | 1.32 | 1.99 |
| fit (cpm, nz, cv²) LOTO — **gene only** | 3.45 | 3.55 | 3.46 | 3.22 | 2.95 | 2.53 | 2.12 | 1.74 | 1.40 | **2.37** |
| fit (cpm, nz, cv²) **+ `|lfc_pred|` + `|H|`** | 3.45 | 3.54 | 3.47 | 3.21 | 2.95 | 2.54 | 2.12 | 1.75 | 1.40 | **2.37** |
| `|lfc_pred|·√cpm` (188's form, our lfc) | 3.33 | 2.91 | 2.55 | 2.26 | 1.92 | 1.64 | 1.39 | 1.22 | 1.08 | 1.88 |
| shift-AUC(`lfc_pred`) — exact power proxy | 2.87 | 2.58 | 2.31 | 2.09 | 1.84 | 1.64 | 1.45 | 1.30 | 1.15 | 1.81 |
| a15's own emitted p-order (**incumbent**) | 2.37 | 2.18 | 2.09 | 2.07 | 1.93 | 1.75 | 1.52 | 1.35 | 1.17 | **1.52** |
| `\|H\|` (185's key) | 2.27 | 1.88 | 1.59 | 1.40 | 1.22 | 1.12 | 1.05 | 1.02 | 0.99 | 1.55 |
| `\|lfc_pred\|` | 1.14 | 0.93 | 0.86 | 0.78 | 0.72 | 0.73 | 0.75 | 0.81 | 0.89 | **1.00** |
| gene DE-rate from Orion (cross-line) | 1.80 | 1.77 | 1.62 | 1.52 | 1.46 | 1.36 | 1.21 | 1.10 | 1.02 | 1.23 |
| *(not buildable)* gene DE-rate, same line, LOTO | 5.05 | 4.65 | 4.29 | 3.86 | 3.47 | 3.02 | 2.58 | 2.03 | 1.54 | *2.77* |
| **ORACLE** | | | | | | | | | | |
| `\|lfc_real\|·√cpm` | 9.11 | 8.39 | 7.87 | 7.38 | 6.64 | 5.74 | 4.72 | 3.55 | 2.18 | **9.81** |
| shift-AUC(`lfc_real`) | 9.01 | 8.37 | 7.80 | 7.30 | 6.49 | 5.53 | 4.43 | 3.18 | 1.96 | 9.07 |
| `\|h_real\|` | 8.86 | 8.16 | 7.52 | 6.94 | 6.13 | 5.28 | 4.38 | 3.35 | 2.13 | 8.71 |
| `\|lfc_real\|` | 7.11 | 6.33 | 5.72 | 5.20 | 4.52 | 3.87 | 3.24 | 2.63 | 1.96 | 5.16 |
| **CONTROL** | | | | | | | | | | |
| shuffled `\|H\|` within row | 1.03 | 1.00 | 0.99 | 0.99 | 0.99 | 1.00 | 1.00 | 0.99 | 1.00 | 1.04 |
| uniform random | 1.07 | 1.08 | 1.07 | 1.02 | 1.02 | 1.02 | 1.00 | 1.00 | 1.00 | 0.99 |

Precision `p(n)` and recall in the same units are in `work/enrich192.log` / `enrich192.json`.

### 1.1 Six readings, in order of importance

1. **Our prediction contributes nothing to membership.** The fitted key *with* `|lfc_pred|` and
   `|H|` is identical to the gene-only key at every `n` (2.37 vs 2.37 at `n_real_p`, 3.47 vs 3.46 at
   89). This is 188 §5.1's finding — `|p|` moves the true magnitude 9.5% over a 330× range — seen
   from the membership side instead of the magnitude side.
2. **`|lfc_pred|` alone is EXACTLY random** at `n = n_real_p` (fold 1.00) and **anti-informative
   below it** (0.72–0.93 at n = 50–600). Multiplying it into `√cpm` therefore *destroys*
   information: `|lfc_pred|·√cpm` **1.88** against `√cpm` alone **2.05**. **188's key does not
   survive the substitution of our lfc for the reference's.** That is the sentence this mandate was
   written to test, and the answer is no.
3. **An exact power proxy does not rescue it.** A shift-AUC key — the exact Mann-Whitney AUC of each
   gene's own empirical control CPM distribution multiplicatively shifted by the predicted fold
   change, which knows about zero-inflation and discreteness and approximates nothing — reads
   **1.81**, below `√cpm`. With the *reference's* lfc the same key reads **9.07**, so the
   instrument is sound; the input is not.
4. **The incumbent's call set is already a weak detectability set** (1.52) and `|H|` is worse (1.55
   at `n_real_p`, **0.99 at n = 4,800** — exactly random). 185's "`|H|` is no better than random for
   `jac`" is confirmed and localised: `|H|` is mildly enriched only in its first ~50 genes and is
   precisely uninformative over the bulk.
5. **The gene-level ceiling is 2.77 and it does not transfer.** A per-gene DE-rate estimated
   leave-one-target-out *on the same line* reads 2.77; estimated on Orion HCT116 and carried across
   it reads **1.23** — worse than `cpm` measured on the target line's own controls. §4.4 repeats the
   test for the fitted coefficients with the same answer.
6. **The oracle ceiling is 9.81 and `jac`'s replicate anchor is 0.3972.** Even a perfect membership
   predictor tops out where a repeat experiment agrees with itself less than half the time. The gap
   from 2.05 to 9.81 is per-gene directional and magnitude accuracy — §8.1's DIRECTION channel —
   not call-set mechanism.

---

## 2. Phase 2 — the joint frontier

`frontier192.py`, `work/frontier192.log`. `n_real` is scaled **proportionally** per perturbation
(`PowerSurface`), never `Surface(n_real_target=N)`. `reach` is reported through 191's `rc_hybrid`
delivery law — called genes first, `|H|` order inside the block and inside its complement.

### 2.1 `fid` is closed to the call set, oracle included

Over the whole grid — 7 keys × 10 sizes — the best `fid` beats a15's realised value by **+0.00088
raw at `f` = 1.0** (`|lfc_pred|·√cpm` at K = 6,400: 0.50789 vs 0.50701) and by **+0.00405 raw at
`f` = 0.5**, where the maximum belongs to the **ORACLE** `|lfc_real|·√cpm` at K = 3,200 (0.51106).
**Maximum available `fid` from the call set alone: +0.0022 overall, and at the level where it is
largest it belongs to a key we cannot build.** 185's reading 1 ("precision is 0.507–0.517 for every
key at every K") holds with the exact power key added and now with the oracle added.

The consequence for INBOX-5: its central row assumed `a(n_real) = 0.56` and priced `fid` at
**+0.0315**. Measured, `a(n_real)` is **0.502–0.508**. INBOX-5's *pessimistic* row (+0.0035) is the
correct one and it is the ceiling, not the floor. **`fid` is a sign-accuracy member; the call set
cannot reach it** — §5.3's +0.0010 raw is the whole of it, and it arrives as a side effect of
protecting the `|H|` head.

### 2.2 The sizing margin — real, one-sided, and it degenerates

The oracle margin sweep shows the asymmetry INBOX-5 §3 predicted, at `f` = 1.0, calling
`m · n_real_p` per perturbation by the fitted key:

| m | 0.5 | 1.0 | 1.5 | **2.0** | 3.0 | 5.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Δ(`fid`+`jac`) overall | −0.1290 | +0.0316 | +0.0351 | **+0.0361** | +0.0300 | +0.0247 |

**The optimum sits at m ≈ 2, the penalty below 1 is catastrophic and above it gentle.** That is the
oracle row. §3.3 shows the predicted row is worthless, so the buildable margin is not a margin
around a prediction — it is a global constant, and §5's `K` already contains it.

---

## 3. Phase 3 — the four falsifiers

### 3.1 Falsifier 1 — shuffled key at every `n`. **PASS, decisively.**

`work/enrich192.log`, `work/fineK192.log`:

| K | 1,400 | 1,800 | 2,200 | 2,600 | 3,000 | 3,400 | 3,800 | 4,474 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `cpm` key `jac` (f = 0.5) | 0.0559 | 0.0562 | 0.0551 | 0.0540 | 0.0526 | 0.0517 | 0.0506 | 0.0486 |
| **shuffled** `jac` | 0.0203 | 0.0226 | 0.0244 | 0.0257 | 0.0267 | 0.0280 | 0.0289 | 0.0301 |
| incumbent `jac` | 0.0376 | 0.0395 | 0.0407 | 0.0413 | 0.0415 | 0.0414 | 0.0413 | 0.0406 |

The shuffled control also lands at fold **0.99–1.04** at every `n` in §1's table, which is the
correctness check on the statistic itself.

### 3.2 Falsifier 2 — retrodicting v14. **2 of 3 members, both at 0.72×; `reach` inverts.**

`retro192.py`, `close192.py`, `cover192.py`, `assemble192.py`.

v14 = v11 + `--tail-spread 1.0 --keep-top 850`; board `fid` −0.03547, `jac` +0.00427,
`reach` +0.01373 scaled, at a realised challenge `n_pred` of **1,563** against v11's **5,961**
(164 §4). That is **two** variables — the call set's SIZE and its COMPOSITION relative to our own
`|lfc_pred|` ordering — and each of my first two models saw one and missed the other:

* keyed by a15's tiered p-order: `fid` −0.0366 (**1.03×**), `jac` −0.00057 (**sign wrong**);
* keyed by a8's flat p-order (v11's actual law): `jac` +0.00310 (**0.73×**), `fid` +0.0062
  (**sign wrong**).

**The composed model gets both**, each member read on the channel that carries it:

| member | predicted scaled | board scaled | ratio | sign | channel |
| --- | ---: | ---: | ---: | :-: | --- |
| `fid` | **−0.025373** | −0.035470 | **0.72** | OK | SIZE / coverage, **key-free** |
| `jac` | **+0.003095** | +0.004271 | **0.72** | OK | COMPOSITION, off our `|lfc_pred|` head |
| `reach` | −0.018132 | +0.013731 | −1.32 | **WRONG** | composition |

`fid`'s channel is key-free by construction: `fid = mean_p prec_p · min(1, n_pred/n_real_p)`, and
185 §4 measured precision flat at 0.507–0.517 for every key at every `K`, so **the entire shape of
`fid(n_pred)` is the coverage term** — a function of the `n_real` distribution and of nothing else.
`jac`'s rise is reproduced by moving the call set OFF our `|lfc_pred|` top-850, which is exactly
what §1.1 reading 2 predicts, since that head is *depleted* of the reference's DE genes.

**Both members land at 0.72× from one fitted parameter.** That two independent channels share a
scale factor is the strongest evidence here that the model is right and merely cold. I record it
and do **not** apply the 1/0.72 correction to the priced numbers.

`reach` inverts — the same member, on the same arm class, that audit 181 found this instrument
inverts on a21. It is reported as a check in §2 and §5, in raw h1-test units, at two conversion
factors.

### 3.3 Falsifier 3 — retrodicting a21, and the `n_real_p` predictor

a21: `n_pred` 4,420 → 1,141; board `fid` 0.5034 → 0.4945 (−0.029464 scaled), `jac` −0.016583. a21
also changed the **emission law**, so 185's own `a21_like` table supplies the second channel.

| member | call-set channel | emission-law residue (185's `a21_like`) | TOTAL predicted | board | ratio | sign |
| --- | ---: | ---: | ---: | ---: | ---: | :-: |
| `fid` | −0.02582 raw | +0.01119 raw | **−0.048621** | −0.029464 | 1.65× | OK |
| `jac` | −0.00401 raw | −0.00932 raw | **−0.036275** | −0.016583 | 2.19× | OK |

Signs right, sizes 1.7–2.2× hot. The **key-free** coverage model does better on a21 than the keyed
one, and a21 is out of sample for it: predicted raw `fid` **0.48824** against an observed
**0.49478** (−0.0065).

**The sizing predictor — INBOX-1's arm.** Leave-one-target-out linear fit of `log(n_real_p + 1)` on
`log‖H_p‖`, `log Σ|H_p|`, and the size of `H_p`'s own top-1% block:

```
Pearson  +0.2401 (p = 0.030)      Spearman +0.2269
residual sd 2.4682 in log space   ->  11.8x multiplicative error
constant-predictor residual sd 2.5102  ->  the model explains 6% of the variance
```

**On INBOX-1's own criterion** — separate the *negligible / subtle / strong* bins that Arc used to
build the panel — it also fails. True `n_real` terciles at 13 and 536, 82 scored rows:

| true \ predicted | negligible | subtle | strong |
| --- | ---: | ---: | ---: |
| negligible | 11 | 10 | **6** |
| subtle | 10 | 9 | 8 |
| strong | **6** | 8 | 14 |

Exact-bin accuracy **0.415** against a chance 0.333; `P(predict "strong" | truly "negligible")` =
**0.222** and `P(predict "negligible" | truly "strong")` = **0.214**. INBOX-1 asked for a predictor
that "never confuses negligible with strong"; this one does so **more than one time in five**.

Priced at `f` = 1.0 against a15's realised set, per-perturbation sizing by the fitted key:

| sizing | m=0.5 | m=1.0 | m=1.5 | m=2.0 | m=3.0 | m=5.0 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| **ORACLE** `n_real_p` | −0.1290 | +0.0316 | +0.0351 | **+0.0361** | +0.0300 | +0.0247 |
| **PREDICTED** `n̂_p` | −0.1494 | −0.1298 | −0.1114 | −0.0972 | −0.0794 | −0.0563 |

**The buildable version is 20–100× worse than a single global constant** (§5's +0.009), because an
11.8× error against a one-sided penalty drops a large stratum below `n_real_p`, where the loss is
linear. INBOX-1 is right that per-perturbation sizing is the structurally correct form and right
that it has never been on the board. **It should not go there: the idea is worth +0.036 and its
input does not exist.** Global `K` is not a Jensen mistake here — it is the Bayes action under a
predictor that explains 6% of the variance.

### 3.4 Falsifier 4 — leakage audit. **PASS; the effect is ≤ 0.00009 overall.**

5 of h1-test's 85 targets are in the challenge's 300-gene panel: `ANKZF1`, `MED15`, `MED25`,
`PLAGL2`, `ZNF32`. 292 of the 10,779 h1-test *genes* are panel targets, but genes enter only as
features, never as labels, and every metric already removes the row's own target.

Re-fitting the membership key with those five targets **deleted from every LOTO fit**
(`close192.py`):

| K | SUM(3), all 85 in the pool | SUM(3), 5 panel targets deleted | Δ |
| ---: | ---: | ---: | ---: |
| 1,600 | +0.00526 | +0.00535 | +0.00009 |
| 2,200 | +0.00728 | +0.00726 | −0.00002 |
| 2,800 | +0.00656 | +0.00652 | −0.00004 |
| 3,600 | +0.00548 | +0.00539 | −0.00009 |

**The recommended arm has no fitted parameters at all** — `cpm` and `|H|` are both objects we
already hold — so it is exempt by construction. The fitted logistic is reported for completeness
and is **not** in the shipped design (§4.4: it does not transfer across lines anyway).

---

## 4. What the board's own numbers say

### 4.1 The five-point board series is an instrument, and nobody had used it as one

164 §4 records five same-profile, emission-only arms with a **measured challenge `n_pred`** and a
**board raw `fid`**. Because `fid`'s shape is pure coverage (§3.2), those five points read out the
board's `n_real` distribution directly.

| arm | v14 | v15 | a8 | a15 | v11 |
| --- | ---: | ---: | ---: | ---: | ---: |
| challenge `n_pred` | 1,563 | 1,911 | 3,836 | 4,424 | 5,961 |
| board raw `fid` | 0.49126 | 0.49654 | 0.5021 | 0.5022 | 0.50178 |

### 4.2 Two independent estimates of the board's `n_real` level, and they agree

**(a) From `fid`'s coverage term** (`cover192.py`), fitting `fid(n) = P · mean_p min(1, n/n_real_p)`
with `n_real_p = round(f · n_real_p^{h1})`:

| f | 0.25 | 0.35 | **0.50** | 0.70 | 1.00 | 1.70 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| median `n_real` | 22 | 31 | **45** | 63 | 90 | 153 |
| RMS vs the five board points | 0.00432 | 0.00378 | **0.00126** | 0.00427 | 0.01230 | 0.02407 |

Best **f = 0.5**, fitted precision 0.50109. The flat-coverage null (`f ≤ 0.25`, coverage never
engages) sits at RMS 0.00432, so coverage explains ~71% of the board's own 0.011 `fid` range. A
**tail-capped** family (`n_real_p = min(cap, ·)`, motivated by a21's light-tail evidence) is
strictly worse — best cap 2,000 at RMS 0.00342 — and it misses a21 by −0.028 where the proportional
family misses by −0.0065.

**(b) From a15's own board `jac` 0.030080** at `n_pred` 4,474 (`price192.py`):

| f | 0.30 | 0.40 | 0.50 | 0.65 | 0.80 | 1.00 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| model `jac` / board `jac` | 0.91 | 1.14 | 1.35 | 1.64 | 1.90 | 2.20 |

Ratio 1.0 at **f ≈ 0.35**.

**The two agree to a factor of 1.4 in `f`**, from two different members and two different
mechanisms. Band **f ∈ [0.35, 0.50]**: median `n_real` **31–45**, mean **270–383**. Below audit
181's board-implied median of 90–110 and well below h1-test's 89/765. It is also **consistent with
INBOX-1's external fact** that the panel was deliberately binned *negligible / subtle / strong* by
DE-gene count: a designed low-`n_real` stratum is exactly what a small median with a long mean
looks like.

Caveat, stated plainly: a one-parameter fit to five points inside an assumed shape family, with an
advantage over the no-coverage null of RMS 0.00126 vs 0.00432. A measurement with a wide error bar,
not a constant. **Everything is therefore priced across the band, and at f = 0.65 as well.**

### 4.3 What that does to "the call-set channel is closed"

`docs/LOG.md` 2026-09-01 concluded from a21 that "`n_real ~240`, coverage did not engage at
`n_pred` 1,141, therefore the call-set channel is closed." The first half survives in spirit — the
board's `n_real` is small. The conclusion does not. Coverage **did** engage: at f = 0.5 the mean
coverage at `n_pred` 1,563 is 97.8% of its saturated value, which is exactly the −0.011 raw `fid`
the board shows between v14 and a15. And the conclusion conflated `fid` with `jac`. **`fid` is
closed** (§2.1 — no key reaches it, oracle included). **`jac` is not**, and its 1/span is 2.72.

### 4.4 The cross-line transfer, which is the buildability question

`orion192.py`. The 3-feature logistic fitted on **Orion HCT116-thinned**'s own reference DE calls
and its own control cells, then carried to h1-test using H1's own control statistics:

| | log cpm | log detection rate | log cv² |
| --- | ---: | ---: | ---: |
| fitted on h1-test (LOTO) | +1.673 | −0.091 | **+1.246** |
| fitted on Orion HCT116 | +0.892 | −0.279 | **−0.398** |

The cv² coefficient **changes sign**. Enrichment at `n = n_real_p`: same-line LOTO **2.37**,
Orion-transferred **2.01**, plain `cpm` **2.05**. **The fitted refinement does not transfer; plain
control CPM does, because it is measured on the challenge's own 18,400 NTC cells.** That is why the
shipped design uses `cpm` and no fitted parameter.

Absolute DE members are never quoted from the Orion surfaces (AGENTS.md harness warning: the
thinned pair carries a 26.4% / 44.1% target-independent technical blob). Only the transferred
*coefficients* are used, and they are scored on h1-test.

---

## 5. Phase 4 — the price

### 5.1 Three members are pinned, by measurement

| member | Δ overall | provenance |
| --- | ---: | --- |
| `pds` | **+0.000000** | 191 falsifier 1 — 900/900 delivered per-gene panel totals digest-identical; **board-verified 2026-09-03** by `rc_hybrid` to sixteen significant figures |
| `mse` | **+0.000000** | 185 §6 — `claim` at the reference cap 0.24653 for every rc arm; `rc_hybrid`'s board `mse` reproduced a15's to sixteen figures |
| `nmae` | **+0.000000** | 191 falsifier 2 — +0.00011 raw h1-test; `rc_hybrid`'s board `nmae` moved −0.00003 raw |

The call set is realised by prescribing which genes clear BH inside 191's emitter at unchanged
per-gene totals, so all three pins carry over unchanged.

### 5.2 The naive re-key: real, but it trades `reach` for `jac`

Re-keying the call set to control CPM (`work/final192.log` §B) buys `jac` and **loses `reach`, and
the loss is a KEY effect, not a size effect**: at fixed `K` = 4,474, changing only the key costs
−0.00797 raw h1-test `reach` at f = 0.5. That is 185 §4 reading 2 — "the order that maximises
`reach` and the set that maximises `jac` are different objects" — priced. At the board's own 1.290
factor the drag is −0.0015 to −0.0051 overall, comparable to the whole `jac` gain:

| f | K=2,200 | K=2,600 | K=3,400 | K=4,474 |
| --- | ---: | ---: | ---: | ---: |
| 0.35 | +0.00652 | +0.00579 | +0.00291 | +0.00100 |
| 0.40 | +0.00506 | +0.00475 | +0.00223 | +0.00217 |
| 0.50 | +0.00233 | +0.00307 | +0.00101 | +0.00084 |
| 0.65 | −0.00169 | −0.00012 | −0.00044 | −0.00081 |

(overall, at the 1.290 factor; the 0.4277 column is in `work/final192.log`.)

### 5.3 The protected design, which resolves the conflict. **This is the arm.**

`work/protect192.log`. Protect the `|H|` head **inside** the call set:

```
call = (top M by |H|)  ∪  (top K−M by control CPM)
```

The `|H|` head is the block 172 measured at sign accuracy **0.624** against a pool mean of 0.537,
so promoting it into the reach walk's front is exactly the thing rc_hybrid's ordering wants.

**Overall Δ at the 1.290 factor, marginal to `rc_hybrid` at `n_pred` 4,474:**

| K | M | f = 0.35 | f = 0.40 | f = 0.50 | f = 0.65 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2,200 | 0 | +0.00652 | +0.00506 | +0.00233 | −0.00169 |
| 2,200 | **400** | +0.01063 | +0.01004 | +0.00824 | +0.00239 |
| **2,600** | **400** | **+0.01109** | **+0.01011** | **+0.00910** | **+0.00437** |
| 2,600 | 800 | +0.00986 | +0.00972 | +0.00891 | +0.00416 |
| 3,400 | 800 | +0.00835 | +0.00794 | +0.00776 | +0.00527 |
| 3,400 | 1,600 | +0.00869 | +0.00826 | +0.00830 | +0.00544 |
| 4,474 | 1,600 | +0.00695 | +0.00660 | +0.00690 | +0.00463 |

**At the conservative 0.4277 factor**, `K` = 2,600 / `M` = 400 reads +0.00787 / +0.00766 / +0.00728
/ +0.00421 across the same four levels — **positive everywhere, at both factors.**

Members at the design point (`K` = 2,600, `M` = 400), h1-test raw, at f = 0.50 / 0.35:

| member | baseline (`rc_hybrid` @ 4,474) | arm | Δ raw |
| --- | ---: | ---: | ---: |
| `fid` | 0.50771 / 0.50771 | 0.50868 / 0.50902 | **+0.00097 / +0.00131** |
| `jac` | 0.04058 / 0.03082 | 0.05346 / 0.04306 | **+0.01288 / +0.01224** |
| `reach` | 0.26126 / 0.28910 | 0.27249 / 0.30892 | **+0.01123 / +0.01982** |

**All three move the right way at once**, which no arm in this family has done before, and `fid`
rises for the first time anywhere in this experiment.

### 5.4 The two-way toll

`rc_hybrid`'s board point refutes rule 13's −0.0146 toll **for profile-preserving arms**, and this
arm is profile-preserving by the same construction and the same digest proof. The coordinator's
amended rule also charges the toll when an arm "collapses the call set"; this arm cuts `n_pred`
4,474 → 2,600 (1.7×, against a21's 3.2× with recall 0.669 → 0.271), which is the **untested** half
of the rule. So, both ways, at `K` = 2,600 / `M` = 400:

| | f = 0.35 | f = 0.40 | f = 0.50 | f = 0.65 |
| --- | ---: | ---: | ---: | ---: |
| **no toll** (profile preserved; `rc_hybrid`'s finding) | +0.01109 | +0.01011 | +0.00910 | +0.00437 |
| **with the toll** (−0.0146 scaled `reach` = −0.00243) | +0.00866 | +0.00768 | +0.00667 | +0.00194 |
| ship gate, shape-only class | +0.0019 | +0.0019 | +0.0019 | +0.0019 |

**It clears the gate at every level, at both `reach` factors, with and without the toll.** The
single case that lands on the gate rather than above it is the worst corner — toll charged *and*
`f` = 0.65, i.e. an `n_real` level both of §4.2's estimates say is too high.

### 5.5 Verdict — GO, as a rider

**GO.** Not as its own slot: the call set is a parameter of an emission law that is already being
changed, so it should ride the next `rc_hybrid`-family build at zero marginal toll. Configuration
to preregister, one variable against `rc_hybrid`:

```
call set per (context, perturbation) = top 400 by |H|  ∪  top 2,200 by control CPM   (n_pred 2,600)
order                                = 191's rc_hybrid law, unchanged
everything else                      = rc_hybrid, unchanged
```

Pre-committed falsifiers for whoever builds it: (1) per-gene panel-total digests identical on
900/900 blocks — else `pds`/`nmae` are not pinned and the arm is void; (2) challenge `n_pred`
within ±3% of 2,600; (3) board `pds` and `mse` reproduce `rc_hybrid`'s to ≥ 12 significant figures;
(4) board `jac` **rises** — it is 88% of the predicted gain and the one member with no prior board
point in its favour.

---

## 6. What is now known that was not

1. **188's key does not survive the substitution of our `lfc` for the reference's.**
   `|lfc|·√cpm` at 0.7219 is an oracle; the buildable form reads **1.88** fold against **2.05** for
   the control CPM alone. The mandate's central hypothesis is measured and falsified.
2. **The membership signal available to us is entirely gene-level.** Adding `|lfc_pred|` and `|H|`
   to a gene-level logistic changes the enrichment by 0.00. `|lfc_pred|` alone is exactly random at
   `n_real_p` and *anti*-informative in its head.
3. **`fid` is closed to the call set, now with an oracle bound.** Over 7 keys × 10 sizes × 4
   `n_real` levels, including the reference's own `|lfc|·√cpm`, `fid` never exceeds a15's by more
   than **+0.0041 raw**. It is a sign-accuracy member; §8.1's DIRECTION channel owns it.
4. **`n_real_p` is not predictable at build time** — 6% of the variance, and it confuses
   *negligible* with *strong* 22% of the time on INBOX-1's own bin criterion. Per-perturbation
   sizing is worth **+0.036 with an oracle** and **−0.097** as built, and is dominated by a global
   constant. Closed, by its input.
5. **The board's own `n_real` distribution is measurable from the five-point `fid` series**, because
   `fid`'s shape is pure coverage. Two independent members put it at **0.35–0.5×** h1-test's,
   median **31–45**. This is the quantity rule 9 says our harnesses cannot read and that has cost
   four board predictions; it should now be *used* rather than assumed.
6. **v14's `fid` and `jac` are both retrodicted at 0.72×** by a composed model with one fitted
   parameter — but only once the arm is split into a SIZE channel and a COMPOSITION channel. Either
   alone gets one member right and the other's sign wrong. Any future call-set retrodiction must
   carry both.
7. **The set that maximises `jac` and the order that maximises `reach` conflict, and the conflict is
   repairable.** Re-keying the call set to detectability costs −0.008 raw `reach` at fixed size;
   protecting the top 400 by `|H|` inside the call set recovers it and adds **+0.011 raw** on top,
   because the reach walk promotes the called block and that block is now our highest-sign-accuracy
   genes. 185 §3.4's hybrid resolved this for the *order*; this resolves it for the *set*.
8. **`docs/LOG.md`'s "the call-set channel is closed" (2026-09-01) needs splitting.** `fid` is
   closed; `jac` is not; and coverage *did* engage on v14. The a21 reading generalised a `fid`
   observation onto a member `fid` shares nothing with.

---

## 7. Files

| file | what it produces |
| --- | --- |
| `ka192.py` | §0 — the Phase-0 gates, 188's statistic re-derived and re-implemented, the donor/panel overlap audit |
| `keys192.py` | the build-time feature block: `h ↔ lfc`, per-gene control detection rate / cv² / quantile grid, and the exact shift-AUC power table |
| `enrich192.py` | §1 — **the enrichment curve**: 16 keys × 9 sizes × 3 statistics, plus the LOTO logistic |
| `frontier192.py` | §2 — the (fid, jac, reach) frontier over key × size × `n_real` level, and the `n_real_p` predictor |
| `retro192.py` | §3.1, §3.3 — the keyed board-series fit, the a21 decomposition, the shuffled control at every size |
| `close192.py` | §3.2, §3.4 — the v14 composition test, the leakage re-fit, the priced optimum |
| `cover192.py` | §4.2(a) — the **key-free** fit of the board's `n_real` distribution from `fid`'s coverage term |
| `price192.py` | §4.2(b) — the second `n_real` estimate from a15's board `jac` |
| `orion192.py` | §4.4 — the cross-line coefficient transfer |
| `assemble192.py` | §3.2 — the composed v14 retrodiction, and the pre-broadcast ledger |
| `work/final192.log` | §5.2, §3.3 — the post-broadcast ledger at both `reach` factors, the SIZE-vs-KEY decomposition, the strength-bin test |
| `work/protect192.log` | §5.3 — **the protected design**, the (K, M) grid at four `n_real` levels |
| `work/fineK192.log` | §3.1, §5.2 — the fine `K` sweep and the shuffled control at each `K` |

Every table above is reproducible by re-running the script named beside it. `work/*.json` carries
the numbers machine-readably.

---

## 8. What I did NOT measure, and the limitations that matter

| not measured | why | can it move the verdict? |
| --- | --- | --- |
| the arm through an actual emission | priced as a prescribed call set on 185's frozen h1-test DE tables. 191 proved the emitter realises any prescribed set at bit-identical per-gene totals (900/900 digests; challenge `n_pred` +0.13%), and `rc_hybrid`'s board point confirms the profile identity — but no `de192.*` table was emitted here. | Not the sign. It would move `fid`/`jac` by the emitter's realisation error, measured at +0.13% on `n_pred`. **The builder must re-check falsifier 1.** |
| `mse` re-run through 096's `claim` on this arm | inherited from 185 §6 (`claim` at the reference cap for every rc arm — a bound, not an estimate) and now board-confirmed by `rc_hybrid` to sixteen figures. Re-running needs the challenge controls dense, 11.9 GB, at the RSS ceiling. | No. |
| the challenge's true `n_real` | needs the withheld reference. §4.2 measures it from the board at 0.35–0.5× h1-test's with a wide error bar. | **This is the binding uncertainty.** The design point is positive from f = 0.35 to 0.65; it was not tested above 0.65 in the protected form. |
| v14 emitted on h1-test | would need a fresh `--keep-top 850 --tail-spread 1.0` emission; the composition test used three constructed call sets at v14's own realised `n_pred` instead. | It would replace a constructed composition with a realised one and could move the 0.72× factor. |
| `M` finer than {0, 400, 800, 1,600} | four points already bracket a flat optimum (+0.0089 to +0.0111 across 400–1,600 at K ≥ 2,600). | No — the surface is flat, which is why 400 is safe to preregister. |
| `cell-eval2`, any build, any submission | hard constraints. Never run. | — |

**One of my own controls failed and is reported as a failure.** My first shift-AUC table used the
`log2` shift values as linear multipliers (`np.log1p(X*c)` with `c` negative), producing NaNs and an
oracle AUC key that scored **0.0891** — *below* the random baseline of 0.0736, for a key that must
be the best in the table. It was caught only because the oracle row was implausible, which is the
whole argument for carrying oracle rows: a *buildable* key at 0.0891 would have looked entirely
reasonable and I would have concluded the power proxy was useless. Fixed, the oracle shift-AUC
reads **0.6671** at `n_real_p`, fold 9.07.

**And one number in this write-up changed after it was first drafted.** The pre-broadcast ledger
(`assemble192.py`, §3.2's table) priced `reach` at the 0.4277 factor and charged the −0.0146 toll,
giving +0.0027 to +0.0056 for the naive re-key. The 11:45 UTC `rc_hybrid` result made both
constants wrong in opposite directions. Both versions are on disk; §5 is the current one.

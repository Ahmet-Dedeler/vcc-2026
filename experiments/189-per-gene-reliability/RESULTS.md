# 189 — the per-gene cross-line reliability prior: **REAL, ORTHOGONAL, AND NOT SHIPPABLE**

Date **2026-09-03**. Measurement only. No `vcc`, no `cell-eval2`, no build, no `vcc prep`, no
submission, no `submissions/` write, no file outside this directory modified. Peak RSS **7.5 GB**
(cap 10, four sibling agents).

---

## Verdict in seven lines

1. **The pilot reproduces, 18/18 rows, 16 of them to four decimals.** Three unstated conventions
   are now pinned (§1.1) and the reach known-answer rows reproduce exactly.
2. **The prior is real and it is enormous on the donor surface.** Leave-one-LINE-out over five
   cell lines: held-out decile spread **+0.1047** against a permutation null of **−0.0003 ±
   0.0025**, i.e. **z = +41.6**, and it beats control expression (+0.0809) on every one of the
   20 fold×pair rows but three.
3. **It is orthogonal to `|H|`** — `spearman(theta_g, gene-mean |H|) = +0.031` — and it **does**
   order sign accuracy where `reach` reads: inside the **first `|H|` decile** of h1-test's DE
   pool the top-theta tercile signs **0.6694 against 0.5873**, a gap of **+0.0820 ± 0.0403**
   (P(>0) = 0.977, z = +2.26 against its own permutation null).
4. **And it buys essentially nothing on any member.** Best ordering key `|H| × exp(+0.5 theta)`
   reads `reach` **0.23941** against `|H|`'s **0.23533** — **+0.00408 ± 0.00424, P(better) 0.828**,
   with lambda selected on the evaluation surface. Call-set precision: **the gain is entirely
   reproduced by the shuffled-theta control** (z = −0.04). `q_w`: **+0.00026**, worth **+0.000045
   overall**.
5. **The mechanism for (3)-versus-(4) is `reach`'s 90% bar.** The median pure prefix is 3 genes
   (172 §7.2). Lifting front accuracy 0.587 → 0.669 does not cross 0.90, so `k*` does not move;
   reordering only inside the top 2–10% of `|H|` is worth **≤ +0.0005** raw `reach`.
6. **The pilot's expression finding is an artefact of its gate, and 171 was right.** Replace the
   "above this line's own median |effect|" gate with a significance-like gate and
   `spearman(r_HE, expression)` goes **+0.5215 → −0.0696 → −0.3353**, expression's held-out
   decile spread goes **+0.0505 → −0.0186 → −0.0457**, while the prior's goes **+0.0402 →
   +0.1094 → +0.1422**. On h1-test, inside `|H|` decile 1, high expression signs **0.0724
   WORSE** (P(>0) = 0.019), and every expression tilt of `|H|` loses `reach` (best −0.01001).
   **185 should keep targeting `|H|`. The expression dominance of our p-order is a defect, not
   half signal.**
7. **Composed slot gate (rule 13): VETO.** +0.00077 overall from `reach` and +0.000045 from
   `nmae`, against a −0.00243 emission toll for delivering any key, = **−0.0016 net** against a
   **+0.0019** shape-only threshold.

**One thing outlives the verdict:** the honest per-GENE ceiling on h1-test is **+0.01065 `q_w`
= +0.0019 overall** (cross-fold oracle), and the in-fold oracle that a naive screen would report
is **22.3× inflated**. §6.3. Anyone proposing a per-gene weight for 188 should read that number
first.

---

## 0. Instruments, surfaces, and the known-answer gate

| | |
|---|---|
| deployment surface | **h1-test**, a15 emission — the only surface board-validated for `reach` on emission-law changes (157 §1) |
| DE table | `experiments/164-callset-artefact/de164.a15.npz` (85 × 10,779) |
| effect | `experiments/169-sharpening-six-member/work/h1test-effect.npz` (`H`, `b`, `cpm`, `gate`) |
| instrument | `experiments/172-sign-confidence-key/harness172.py`, **imported unmodified**, plus `reach_vec` (per-target, for the paired bootstrap) and `qw` |
| fitting surface | 5 cell LINES, 10 cross-line pairs, 10,725-gene working axis |
| peak RSS | 5.3 GB (pilot) / 6.9 GB (build) / 7.5 GB (tension). The first build ran at 13.8 GB and was rewritten gene-chunked; the chunked output is **bit-identical** on all 45 arrays. |

**Every script in this directory calls `harness189.check()` before it prints anything and aborts
if a row moves** (rule 10). All rows passed on every run:

| row | here | expected |
|---|---:|---:|
| `reach(production order)` | **0.18283** | 0.18283 ✓ |
| `reach(\|H\|)` | **0.23533** | 0.23533 ✓ |
| `reach(sig-first, then \|H\|)` | **0.21914** | 0.21914 ✓ |
| `reach(ORACLE)` | **0.64553** | 0.64553 ✓ |
| INBOX-4 §1 decile table, `\|H\|` and expression | both reproduce | ≤ 0.002 ✓ |
| `reach_vec().mean() == debench.reach()` | **0.18283**, 82 scored targets | exact to 1e-12 ✓ |

The **H-oracle** is **0.61862** on this instrument (172 §4.2); the mandate quotes the INBOX's
0.62119, which differs by a tie convention in the 0/1 key. Gaps below are quoted against 0.61862.

---

## 1. Phase 0 — the pilot reproduces

`pilot189.py`, `work/pilot189.json`. Deterministic; the whole table in one pass.

| row | here | pilot | delta |
|---|---:|---:|---:|
| common genes | **7,976** | 7,976 | 0 |
| three-way targets | **9,073** | 9,073 | 0 |
| genes with support | **7,916** | 7,916 | 0 |
| `r_HE` mean / sd | **0.5379 / 0.0467** | 0.5379 / 0.0467 | 0.0000 |
| spearman(`r_HE`, `r_HK`) | **+0.3833** | +0.3833 | 0.0000 |
| spearman(`r_HE`, `r_EK`) | **+0.3359** | +0.3359 | 0.0000 |
| spearman(`r_HK`, `r_EK`) | **+0.4763** | +0.4763 | 0.0000 |
| held-out decile of `r_HE` on `r_HK` | **0.5515 / 0.5113 / +0.0402** | same | 0.0000 |
| spearman(`r_HE`, `b_ctrl`) | **+0.5215** | +0.5215 | 0.0000 |
| spearman(`r_HE`, mean \|effect\|) | **+0.1125** | +0.1125 | 0.0000 |
| spearman(`r_HE`, across-target sd) | **+0.1578** | +0.1578 | 0.0000 |
| partialled, resid `r_HK` | **+0.2133** | +0.2133 | 0.0000 |
| partialled, resid `r_EK` | **+0.2251** | +0.2251 | 0.0000 |
| held-out decile spread, residual prior | **+0.0174** | +0.0173 | +0.0001 |
| held-out decile spread, expression alone | **+0.0505** | +0.0505 | 0.0000 |

### 1.1 Three conventions the pilot table did not state, now pinned

My first reading missed four rows. Rather than tune until they matched, each convention was
resolved by enumerating the alternatives and taking the one that lands on the published value.

| quantity | convention that reproduces | alternatives tried |
|---|---|---|
| `mean \|effect\|`, `across-target sd` | mean over the **HE pair only** (the pair whose agreement is being explained), centred blocks | 32, in `work/probe_confounder.py`; the three-line mean reads +0.0268 / +0.0748 |
| `b_ctrl` | **HCT116's alone** — the line common to `r_HE` and the held-out `r_HK` | 11, in `work/probe_bctrl.py`; the three-line mean reads +0.5238, the HE mean +0.5581 |
| decile means | **unweighted** over genes | support-weighted reads 0.5555 / 0.5126 / +0.0429 |

Worth recording for its own sake: **`mean |effect|` on the K562 block alone correlates −0.2389
with `r_HE`, and on the HEK293T block +0.2146.** The pilot's +0.1125 is an average across a sign
change, so "magnitude predicts reliability" is not a stable statement even inside the pilot.

### 1.2 One number the pilot did not report, and it matters

**259 of the 300 challenge panel targets sit inside the three-way pool.** 161's `core5` had zero
by construction. §7 (F4) is the control; it changes nothing, but the property is lost and is
recorded rather than assumed away.

---

## 2. Phase 1 — building the prior properly

### 2.1 The fitting inventory (`inventory189.py`)

Six donor tables, **five distinct cell LINES** — `k562_gwps` and `k562_essential` are two
libraries of one line (161 `common.py`; 159 measured a same-line pair at +0.2712 against
+0.03..+0.16 for genuine cross-line pairs), so a fold that splits them leaks. `h1_2025` is
**BANNED** in `common189.py` and is never opened.

| line | table | targets | genes on the h1-test axis | panel targets |
|---|---|---:|---:|---:|
| HCT116 | `orion_hct116` | 17,767 | 10,602 | 290 |
| HEK293T | `orion_hek293t` | 17,977 | 10,602 | 294 |
| K562 | `k562_gwps` | 9,787 | 7,101 | 269 |
| RPE1 | `rpe1` | 2,204 | 7,570 | 0 |
| JURKAT | `nadig_jurkat` | 2,270 | 7,560 | 0 |

Working gene axis = h1-test's DE axis ∩ genes carried by ≥2 tables = **10,725 of 10,779**.
Ten cross-line pairs, 1,986 to 17,698 shared targets each.

### 2.2 The estimator, and why raw pooling is wrong

**The ten pairs do not share a base rate**: pooled agreement runs **0.5221** (HEK293T|K562) to
**0.6147** (K562|JURKAT), and genes are carried by different pair SUBSETS (10,602 genes in the
Orion pair against 6,414 in K562|RPE1). Pooling raw counts therefore ranks a gene partly by
*which pairs happen to carry it*. Three estimators are carried through every table:

```
RAW    r_g = sum_p k_gp / sum_p n_gp                       the pilot's statistic
EB     beta-binomial posterior mean, (a,b) by marginal MLE  the mandate's shrinkage
HIER   logit(q_gp) = m_p + theta_g,  theta_g ~ N(0,tau^2)   the delivered prior
```

`tau` is fitted by the standard Laplace EB update: **tau = 0.3418** on all ten pairs.

### 2.3 What the shrinkage does: **nothing, and the arithmetic says why**

| | |
|---|---:|
| supported genes (`n_g` > 50) | 10,620 |
| `n_g` p10 / p50 / p90 | 4,686 / **13,141** / 13,345 |
| between-gene sd of raw `r_g` | **0.06671** |
| mean per-gene binomial SE `0.5/sqrt(n)` | **0.00640** |
| signal-to-noise | **10.4×** |
| fitted beta prior `a+b` | **35.9** → shrink weight at the median `n_g` = **0.99728** |
| `spearman(raw, EB)` | **+0.999882** |
| held-out decile spread: raw → EB → hierarchical | **+0.1040 → +0.1045 → +0.1047** |

**The mandate's instinct is right and the answer is that there is nothing to shrink.** With
~13,000 supported pairs per gene the binomial sampling error is 10× smaller than the real
between-gene spread, the marginal likelihood answers with a prior worth 36 pseudo-observations,
and the ORDER is unchanged to six decimals — so no endpoint in this experiment can move under
shrinkage, and none does. *(The counts are not independent Bernoulli draws — targets stay
correlated after centring — which is exactly why the marginal MLE lands on `a+b = 36` rather
than the method-of-moments 312: the excess dispersion is real between-gene signal.)*

### 2.4 Leave-one-LINE-out, and the scaling in fitting lines

`prior189.py`. Fit on pairs with neither endpoint = L, score on the pairs that involve L.

| held-out L | mean rho | prior d10 | raw d10 | EB d10 | **expression d10** |
|---|---:|---:|---:|---:|---:|
| HCT116 | +0.5160 | +0.0915 | +0.0904 | +0.0917 | +0.0730 |
| HEK293T | +0.4823 | +0.0775 | +0.0772 | +0.0782 | +0.0466 |
| K562 | +0.5128 | +0.0942 | +0.0935 | +0.0935 | +0.0800 |
| RPE1 | +0.4503 | +0.1275 | +0.1259 | +0.1260 | +0.0935 |
| JURKAT | +0.5287 | +0.1330 | +0.1329 | +0.1329 | +0.1116 |
| **MEAN** | **+0.4980** | **+0.1047** | +0.1040 | +0.1045 | **+0.0809** |

**Scaling with the number of fitting LINES** (all 5 held-out lines × all subsets of the
remaining lines of that size, `prior189.py` §3):

| fitting lines | fitting pairs | subsets | mean rho | mean d10 spread |
|---:|---:|---:|---:|---:|
| 2 | 1 | 30 | +0.4104 | +0.0868 |
| 3 | 3 | 20 | +0.4653 | +0.0979 |
| 4 | 6 | 5 | **+0.4980** | **+0.1047** |

Monotone. **Adding fitting lines helps on the donor surface — and §3.3 shows none of it reaches
h1-test.**

---

## 3. Phase 2 use 1 — the prior as an ORDERING KEY

`screen189.py`, `work/screen189.json`. Every key is `|H|` re-weighted by a gene-level scalar, so
it changes only the within-target order — the only thing `reach` lets us move (171 §0.1).

| key | **reach** | Δ vs `\|H\|` | SE | P(better) | top1% | top5% | macro1% | d1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| production order | 0.18283 | −0.05250 | 0.01540 | 0.000 | | | | |
| theta alone (gene-only) | 0.17450 | −0.06083 | 0.01995 | 0.000 | 0.522 | 0.540 | 0.591 | 0.551 |
| expression alone (gene-only) | 0.17589 | −0.05945 | 0.01476 | 0.000 | 0.555 | 0.538 | 0.579 | 0.540 |
| **`\|H\|` — THE INCUMBENT** | **0.23533** | — | — | — | **0.707** | 0.651 | **0.826** | 0.621 |
| `\|H\|` × raw `r_g` | 0.23781 | +0.00248 | 0.00284 | 0.804 | 0.716 | 0.655 | 0.831 | 0.619 |
| `\|H\|` × exp(+0.25 theta) | 0.23691 | +0.00158 | 0.00268 | 0.698 | 0.713 | 0.653 | 0.830 | 0.623 |
| **`\|H\|` × exp(+0.50 theta)** | **0.23941** | **+0.00408** | 0.00424 | **0.828** | 0.711 | 0.652 | 0.810 | 0.619 |
| `\|H\|` × exp(+1.00 theta) | 0.23154 | −0.00379 | 0.00543 | 0.238 | 0.639 | 0.632 | 0.766 | 0.614 |
| `\|H\|` × exp(−0.25 z_expr) | 0.22532 | −0.01001 | 0.00534 | 0.018 | 0.693 | 0.632 | 0.813 | 0.606 |
| `\|H\|` × exp(+0.25 z_expr) | 0.22484 | −0.01049 | 0.00614 | 0.037 | 0.680 | 0.632 | 0.789 | 0.606 |
| **ORACLE** | 0.64553 | +0.41019 | | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

Full 11 × 7 joint `(lambda, mu)` grid in `work/screen189.json`; its maximum is the same
0.23941, and **every cell with `mu != 0` is below the `mu = 0` column at the same lambda**.

**Three readings.**

* **The best row is +0.00408 with SE 0.00424 and lambda selected on h1-test over ten values.**
  172 §4.1 measured the paired-bootstrap SE of a reach difference at 0.005–0.016; this sits
  inside it, exactly as 172's own best learned key did (+0.00753 ± 0.00965).
* **Front precision is not improved.** Pooled top-1% 0.711 vs 0.707; macro top-1% **0.810 vs
  0.826, i.e. worse**. 172's calibration — *reach reads only the front, so precision at the top
  is the objective* — comes out negative here. Applied unconditionally, theta trades `|H|`
  information roughly one-for-one (§6.1 shows the conditional version does not).
* **`theta` alone and expression alone are both catastrophic** (0.1745, 0.1759 — below the
  production order's 0.18283). A gene-only key has no within-target structure beyond a fixed
  permutation; this is the same shape as 172 §6's target-level block.

### 3.3 More fitting lines does not help where it counts (`scaling189.py`)

| fitting set | pairs | tau | reach | Δ vs `\|H\|` | SE | P(better) | top5% |
|---|---:|---:|---:|---:|---:|---:|---:|
| 2 lines — **the pilot's own pair** (HCT+HEK) | 1 | 0.3970 | 0.23802 | +0.00269 | 0.00472 | 0.708 | 0.6457 |
| 2 lines (HCT+K562) | 1 | 0.0829 | 0.23777 | +0.00244 | 0.00170 | 0.942 | 0.6530 |
| 3 lines — the pilot's three-way | 3 | 0.3732 | 0.24000 | +0.00466 | 0.00427 | 0.862 | 0.6488 |
| 4 lines (+ RPE1) | 6 | 0.3538 | 0.23921 | +0.00387 | 0.00418 | 0.815 | 0.6515 |
| **5 lines — the delivered prior** | 10 | 0.3418 | 0.23941 | +0.00408 | 0.00424 | 0.828 | 0.6524 |

**Flat.** The donor-side +0.0868 → +0.1047 buys **+0.0014 raw reach**, one third of a standard
error. The apparent P = 0.942 on HCT+K562 is a `tau` artefact — that fit has `tau = 0.083`, a
4.8× weaker tilt — and `work/probe_matched_tilt.py` z-scores every prior to a common tilt
strength and collapses it: at matched `s = 0.15` the four fitting sets read +0.00452, +0.00271,
+0.00409, +0.00406 with SE 0.0040–0.0051.

---

## 4. Phase 2 use 2 — the prior as a CALL-SET FILTER

`members189.py`. I own **which** genes go in the call set at a given size; 185's
`INBOX-1-per-perturbation-sizing.md` owns **how many**, and is read but not duplicated. Because
`fid_p = k / max(n_pred_p, n_real_p)` reads sign precision directly and rule 9 forbids reading
absolute `fid` off this surface, the transferable quantity is **precision at fixed `n`**.

| key | n=10 | 25 | **50** | 100 | 200 | 400 | 800 | 1600 | 3200 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `\|H\|` (incumbent) | **0.5904** | 0.5431 | 0.5244 | 0.5229 | 0.5184 | 0.5191 | 0.5182 | 0.5141 | 0.5103 |
| `\|H\|` × exp(+0.50 theta) | 0.5623 | 0.5451 | 0.5287 | 0.5229 | 0.5175 | 0.5187 | 0.5188 | 0.5144 | 0.5102 |
| `\|H\|` × exp(+1.00 theta) | 0.5087 | 0.5073 | 0.5090 | 0.5143 | 0.5144 | 0.5143 | 0.5173 | 0.5151 | 0.5106 |
| **`\|H\|` × exp(+0.50 SHUFFLED theta)** | 0.5228 | 0.5366 | **0.5297** | 0.5201 | 0.5183 | 0.5167 | 0.5165 | 0.5139 | 0.5102 |
| `\|H\|` × exp(±0.50 z_expr) | 0.529/0.539 | 0.527/0.528 | 0.526/0.529 | | | | | | |

**The prior is worse than the incumbent at n = 10 and worse than its own SHUFFLE at n = 50.**
That is the whole result.

The one construction that looked alive — reordering by theta **inside** the top-X% of `|H|`, so
`|H|`'s stratification is preserved exactly — was priced properly in
`work/probe_stratified_precision.py` with a per-target bootstrap and a 30-seed permutation null:

| n_pred | `\|H\|` precision | stratified | Δ | boot SE | P(>0) | **null Δ** | null sd | **z** |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 25 | 0.5431 | 0.5435 | +0.0005 | 0.0117 | 0.518 | −0.0049 | 0.0089 | +0.60 |
| **50** | 0.5244 | 0.5331 | **+0.0088** | 0.0077 | 0.871 | **+0.0091** | 0.0069 | **−0.04** |
| 100 | 0.5229 | 0.5241 | +0.0012 | 0.0050 | 0.584 | +0.0021 | 0.0031 | −0.31 |
| 200 | 0.5184 | 0.5197 | +0.0012 | 0.0028 | 0.663 | +0.0016 | 0.0023 | −0.14 |
| 400 | 0.5191 | 0.5194 | +0.0003 | 0.0016 | 0.575 | +0.0007 | 0.0012 | −0.32 |

**The +0.0088 at n = 50 is +0.0091 for a shuffled prior.** Any perturbation of which genes fall
in the top 50 buys that; theta contributes nothing. **Use 2 is closed.**

*To 185:* the call-set family has one live lever left and it is yours — **sizing**, not
composition. Nothing in this experiment competes with INBOX-1.

---

## 5. Phase 2 use 3 — the prior as a WEIGHT on `q_w`

187 §5: with 99% of genes below the `k = 1` kink, `nmae = 1 − (2 q_w − 1) K`, so
**`d nmae = −2 K d q_w`** at fixed total predicted mass. Every weight below is renormalised so
`K` is bit-identical — this stays a direction question and hands the magnitude question to 188,
as the mandate requires. **No magnitude was changed.**

Baseline on this DE table: `q_w` **0.57516**, `K` **0.32140**, count accuracy **0.53435**.
*(187 quotes 0.621 / 0.199 on a23 geometry; this is a15's `de164` table, hence the offset. All
deltas below are within-table.)*

| weight `u` | `q_w` | Δ`q_w` | Δ`nmae` raw | Δ scaled | **Δ overall** |
|---|---:|---:|---:|---:|---:|
| exp(+0.25 theta) | 0.57541 | **+0.00026** | −0.00017 | +0.00027 | **+0.000045** |
| exp(+0.50 theta) | 0.57524 | +0.00008 | −0.00005 | +0.00009 | +0.000015 |
| exp(+1.00 theta) | 0.57174 | −0.00341 | +0.00219 | −0.00362 | −0.000604 |
| exp(+2.00 theta) | 0.54064 | −0.03452 | +0.02219 | −0.03668 | −0.006114 |
| raw `r_g` | 0.57606 | +0.00091 | −0.00058 | +0.00096 | +0.000161 |
| exp(+0.50 SHUFFLED theta) | 0.57358 | −0.00157 | | | |
| exp(±0.50 z_expr) | 0.5741 / 0.5760 | −0.0011 / +0.0009 | | | |

**Maximum achievable: +0.000045 overall.** Down-weighting low-`r_g` genes does not raise `q_w`.

---

## 6. Why a z = +41.6 signal converts to nothing — the mechanism

### 6.1 theta is orthogonal to `|H|` and it *does* order accuracy inside `|H|` strata

`work/probe_mechanism.py`, `work/probe_conditional.py`. Sign accuracy over h1-test's DE pool by
theta tercile, **within each `|H|` decile**:

| `\|H\|` decile | low-theta acc | high-theta acc | diff |
|---:|---:|---:|---:|
| **1 (the front of the reach walk)** | **0.5911** | **0.6623** | **+0.0712** |
| 2 | 0.5268 | 0.5656 | +0.0387 |
| 3–8 | 0.5026–0.5431 | 0.5249–0.5504 | +0.006 … +0.027 |
| 9 | 0.4988 | 0.4998 | +0.0009 |
| 10 | 0.5105 | 0.4863 | −0.0242 |
| **ALL** | 0.5261 | 0.5443 | **+0.0183** |

Decile 1 with error bars, per-target bootstrap and 30-seed permutation null:

| key | lo acc | hi acc | gap | boot SE | P(>0) | null mean | null sd | **z** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **theta (the prior)** | 0.5873 | 0.6694 | **+0.0820** | 0.0403 | **0.977** | +0.0062 | 0.0335 | **+2.26** |
| control expression | 0.6691 | 0.5967 | **−0.0724** | 0.0372 | 0.019 | −0.0057 | 0.0387 | −1.72 |
| gene-mean `\|H\|` itself | 0.6281 | 0.6392 | +0.0111 | 0.0414 | 0.602 | −0.0025 | 0.0283 | +0.48 |

Correlations confirming orthogonality: `spearman(theta, gene-mean |H|)` **+0.031**, inside the
pool **+0.044**; `spearman(theta, expression)` +0.256; `spearman(theta, reference |real_lfc|
marginal)` −0.199.

**So the prior is genuinely new information, it is strongest exactly where `reach` reads, and
`|H|` does not already contain it — the one thing 172 could not say about any of its keys.**

### 6.2 And it still cannot move `reach`, because the bar is 0.90

172 §7.2: the median pure prefix `k*` is **3 genes** and 94.4% of h1-test's `reach` lives in 32
targets. 157 §3b's identity is `oracle reach ≈ accuracy/0.9`. Lifting the front decile from
0.587 to 0.669 does not cross 0.90, so `k*` does not extend. Measured, with `|H|`'s
stratification preserved exactly (`work/probe_conditional.py`):

| arm | reach | Δ vs `\|H\|` | SE | P(better) |
|---|---:|---:|---:|---:|
| top 2% of `\|H\|` reordered by exp(+2.0 theta) | 0.23579 | **+0.00045** | 0.00029 | 0.944 |
| top 5% of `\|H\|` reordered by exp(+2.0 theta) | 0.23583 | +0.00049 | 0.00074 | 0.754 |
| top 10% of `\|H\|` reordered by exp(+1.0 theta) | 0.23580 | +0.00046 | 0.00095 | 0.706 |
| top 25% of `\|H\|` reordered by exp(+0.5 theta) | 0.23674 | +0.00141 | 0.00173 | 0.783 |

**A +0.082 front-accuracy gain, correctly applied, is worth +0.0005 raw `reach`.** That is the
sharpest available statement of why `reach` is model-limited rather than order-limited (154's
verdict), and it is now demonstrated with a key that *is* better at the front rather than one
that is not.

### 6.3 The honest per-GENE ceiling on h1-test, and a 22× trap

`genelevel189.py`. A naive screen computes a gene-level oracle on the same rows it scores. The
split-half is the falsifier: estimate per-gene accuracy on half the targets, read it on the other
half (12 seeds).

| | |
|---|---:|
| genes in the DE pool at all | 10,331 of 10,779 |
| targets per gene in the pool: median / p90 / max | **5** / 13 / 39 |
| split-half `rho(acc_A, acc_B)` | **+0.0097** |
| split-half decile spread A→B | **+0.0044** |
| `q_w` with an IN-FOLD gene oracle | 0.8144 |
| `q_w` with a CROSS-FOLD gene oracle | **0.5880** |
| `q_w` unweighted baseline | 0.5774 |
| **honest gene-level `q_w` headroom** | **+0.01065** → Δ`nmae` −0.00685 raw → **+0.00189 overall** |
| in-fold oracle headroom | +0.23707 — **inflated 22.3×** |

And the direct question: `spearman(theta_g, h1-test per-gene accuracy)` = **+0.0191** over
10,261 genes (+0.0830 restricting to the 1,986 genes in ≥10 pool targets); raw `r_g` +0.0244;
expression +0.0073.

**Two things follow.** (a) The whole per-gene weighting family on h1-test has a ceiling of
**+0.0019 overall even with the answer**, which is the shape-only slot threshold — so it is not a
place to spend a slot. (b) `members189.py`'s +0.0289/+0.0571 "oracle" rows are the inflated
in-fold statistic and **must not be quoted**; they are left in the table with this correction
attached.

---

## 7. Phase 3 — the falsifiers, all four run

### F1 — permutation null

| surface | real | null mean | null sd | verdict |
|---|---:|---:|---:|---|
| **donor LOLO** (20 fold×pair rows, 20 seeds each) | **+0.1047** | −0.0003 | 0.0025 | **PASS, z = +41.6** — 20/20 rows positive, per-row z +17.5 to +53.2 |
| **h1-test ordering key** (50 seeds, lam = 0.5) | **0.23941** | 0.23057 | 0.00389 | **PASS, z = +2.28**, P(null ≥ real) = 0.000 |
| h1-test, matched tilt, 4 fitting sets × 2 strengths | | | | z **+1.31 to +2.51**, all 8 positive |
| **h1-test call-set precision @50** | +0.0088 | **+0.0091** | 0.0069 | **FAIL, z = −0.04** |

The null is the strict one (§172 §5's lesson: a within-group shuffle leaks). theta is permuted
across the whole gene axis, which is the level the key reads.

### F2 — expression-only control

| where | prior | expression | who wins |
|---|---:|---:|---|
| donor LOLO, mean held-out decile spread | **+0.1047** | +0.0809 | prior, on 17 of 20 rows |
| h1-test, best tilt of `\|H\|` over the same grid | **+0.00408** | **−0.01001** | prior; **every** expression tilt loses |
| h1-test, DE-pool decile spread d1−d10 | `\|H\|` +0.1228 | **+0.0118** | expression is flat (known-answer row) |
| h1-test, outside the DE pool | `\|H\|` +0.0136 | +0.0030 | both flat |
| h1-test, inside `\|H\|` decile 1 | **+0.0820** | **−0.0724** | prior; expression has the **wrong sign** |

**Expression is not the headline.** It is a donor-surface artefact — see §8.

### F3 — the h1-test trap

`data/donors/h1_2025` (85/85 coverage) is BANNED in `common189.py` and never opened; the audit
prints the full donor↔surface overlap before every screen. The two Orion lines cover 84/85 and
K562 82/85 and **are** used, so the strict control is a prior refitted with those 85 targets
deleted from every donor pair (`build189.py --drop-h1test-targets`):

| variant | tau | rho vs main | reach@lam=0.5 | Δ vs `\|H\|` | prec@50 | `q_w` |
|---|---:|---:|---:|---:|---:|---:|
| MAIN — all targets | 0.3418 | +1.0000 | 0.23941 | +0.00408 | 0.5287 | 0.57524 |
| **F3 — h1-test's 85 targets deleted** | 0.3410 | **+0.9994** | **0.23941** | +0.00407 | 0.5280 | 0.57524 |
| **F4 — 300 panel targets deleted** | 0.3389 | **+0.9988** | **0.23943** | +0.00410 | 0.5282 | 0.57527 |

**No leak.** A gene-level statistic averaged over 9,000–17,700 targets cannot carry 85 of them.

### F4 — panel disjointness

**259 of the 300 challenge panel targets are in the fitting pool** (812 panel-target instances
across the ten pairs; the K562/RPE1/Jurkat pairs contribute zero by construction). 161's `core5`
had zero and this does not — so the property is **lost, and stated**. The control above shows it
is immaterial: the panel-free prior correlates +0.9988 with the main one and moves `reach` by
+0.00002. Pooled per-pair agreement is identical to 4 decimals with and without the panel
(e.g. HCT116|HEK293T 0.5242 → 0.5241).

---

## 8. The tension with 171, resolved — and it does redirect 185

The pilot's headline (2) was: *control expression is by itself the strongest per-gene
reliability signal, worth +0.0505 between deciles.* That reading is **wrong, and the mechanism
is its gate.**

The pilot gates each entry at *that line's own per-gene median |effect|* — an unconditioned
gate. Replace it with the significance-like gate every scored member actually applies (per
target, top-K by the Poisson-Wald proxy `|delta|·sqrt(cpm)`, 172 §1.1's validated DE-pool proxy,
median recall 0.828 of the real BH set) and the expression signal **inverts** while the prior
**strengthens** (`tension189.py`):

| gate | genes | rho(`r_HE`, expr) | **expr held-out d10** | **prior held-out d10** | mean `r` |
|---|---:|---:|---:|---:|---:|
| pilot: > own per-gene median \|effect\| | 7,916 | **+0.5215** | **+0.0505** | +0.0402 | 0.5379 |
| significance proxy: per-target top 5% | 1,103 | **−0.0696** | **−0.0186** | **+0.1094** | 0.6400 |
| significance proxy: per-target top 1% | 166 | **−0.3353** | **−0.0457** | **+0.1422** | 0.6852 |

**What the pilot measured is measurement precision, not transferability.** Two noisy measurements
of the same gene agree more often when the gene is highly expressed, and an unconditioned gate
lets that in. Every scored member conditions on significance first, and significance is itself
expression-driven — so the channel is already spent by the time any metric reads it.

Three independent confirmations on h1-test, none of which shares a line of code with the above:

* **Expression's DE-pool decile row is flat** — 0.540 … 0.528, d1−d10 = **+0.0118** — and it is
  one of the harness's own known-answer rows, so it was true before this experiment started.
* **Inside `|H|` decile 1, high expression signs 0.0724 WORSE** (P(>0) = 0.019, §6.1). Expression
  is *anti*-informative exactly at the front of the reach walk.
* **Every expression tilt of `|H|` loses `reach`, in both directions** (μ = −0.25: −0.01001;
  μ = +0.25: −0.01049), and the whole 11 × 7 joint grid is maximised at μ = 0.

### The message to 185

**171's diagnosis stands in full. The expression dominance of our emitted p-order
(`rho +0.485` against `−0.391` for `|H|`) is a defect and not half signal. Keep targeting `|H|`.**
The per-gene prior is a *different and orthogonal* object (`spearman(theta, gene-mean |H|)` =
+0.031, `spearman(theta, expression)` = +0.256) — it is not "expression in disguise" — but it is
worth +0.0005 raw `reach` when delivered correctly, so it does not change what 185 should build.

---

## 9. All six members (rule 12) and the composed slot gate (rule 13)

**This experiment builds nothing and emits nothing.** No prediction was rendered, no counts were
emitted, no call set was recomputed. Five of six are pinned by construction:

| member | Δ from this experiment | why |
|---|---|---|
| `pds` | **exactly 0** | reads the pseudobulk profile; none was produced |
| `mse` | **exactly 0** | reads the pseudobulk profile; no counts were emitted |
| `fid` | **exactly 0** | reads our own call set; none was recomputed |
| `jac` | **exactly 0** | reads our own call set; none was recomputed |
| `nmae` | **exactly 0** | our lfc is untouched; §5 is a counterfactual on the weights |
| `reach` | **counterfactual only** | the scorer reads our Wilcoxon p-values, not a key |

**The exemption is narrow and must not be inherited.** The moment a key is *delivered* the
emission law changes and rule 12 binds on the composition — that is 185's pricing, not this
document's.

**Pricing the best arm as if 185 delivered it perfectly:**

```
reach   +0.00408 raw x 1.1292 = +0.00461 scaled  ->  +0.00077 overall
nmae                                                 +0.000045 overall  (§5, best weight)
fid/jac  0        (§4: the gain is the shuffled control's)
pds/mse  0        (nothing built)
emission-law toll (rule 13, audit-181: EMISSION class)  -0.0146/6 = -0.00243 overall
                                                        ------------------------------
                                             predicted   -0.00161 overall
```

against a **+0.0019** shape-only threshold and **+0.0089** amplitude-changing. **VETO**, by an
order of magnitude, and rule 13 additionally records that `reach` is **not locally priceable at
all** (2/6 correct signs on our own board points) — so the +0.00077 is itself the least
trustworthy number in the ledger.

For scale: `|H|` closes **12.0%** of the H-oracle gap over the production order (172 §4.2); the
best key here closes **13.0%**. 172's best learned key closed 13.8%. **Three independent
families now land in the same 12–14% band, and the band is the answer.**

---

## 10. Artefacts

| file | what it produces |
|---|---|
| `common189.py` | the `r_g` definition, donor loaders, beta-binomial EB, rank partialling |
| `pilot189.py` | **§1** — the pilot, with an 18-row abort gate |
| `inventory189.py` | **§2.1** — line/pair inventory and the working gene axis |
| `build189.py` | **§2.2** — `(k_g, n_g)` for 10 cross-line pairs, gene-chunked; `--drop-h1test-targets` builds F3 |
| `prior189.py` | **§2.2–2.4** — hierarchical fit, LOLO folds, line-count scaling |
| `harness189.py` | the instrument: 172's harness + `reach_vec`, `qw`, `mass_ratio`, `audit`, `check` |
| `screen189.py` | **§3** — ordering keys, the joint grid, the first null |
| `members189.py` | **§4, §5** — call-set precision, `fid`/`jac`, `q_w` |
| `genelevel189.py` | **§6.3** — the split-half and the 22× oracle inflation |
| `tension189.py` | **§8** — the 171 tension, both surfaces |
| `falsify189.py` | **§7** — F1–F4 on one set of lines of code |
| `scaling189.py` | **§3.3, §2.3** — h1-test scaling and the shrinkage anatomy |
| `work/probe_*.py` | the convention hunts (§1.1) and the mechanism probes (§6.1, §6.2, §4) |

| artefact | sha256 |
|---|---|
| `work/prior189.npz` (the delivered prior: `theta`, `raw`, `eb`, `K`, `N`, 10,725 genes) | `faa1dc65d39424827c578b24aafdb08c7d00d17bc8f56dc38ce1e8ccbbb95deb` |
| `work/pairs189.npz` (per-pair counts, all + non-panel) | `b5f13664421ad4237c547fcce34a5b80f6ee93011be66761ba05ca86c279b5ee` |

---

## 11. What I did NOT measure, and why

| not measured | why | can it move? |
|---|---|---|
| a **delivered** arm's six members | building an emitter is 185's territory and a build is forbidden here | §9 prices the composition structurally |
| per-perturbation `n_pred` sizing | 185's `INBOX-1`; explicitly not duplicated | it is the live lever, and it is theirs |
| per-gene MAGNITUDE re-shaping | 188's territory; §5 hands over the achievable `q_w` and §6.3 the honest ceiling | ceiling is +0.0019 overall **with the answer** |
| `k562_essential` as a sixth fitting table | same LINE as `k562_gwps`; splitting them leaks (159) | it would add 2,003 same-line targets, not a line |
| a whitened analogue of theta | no leak-free per-fold whitening operator exists (172 §11) | §6.1 shows theta is already orthogonal to `\|H\|`, so the omission is not load-bearing |
| Tahoe / CD4 / KOLF2 lines | 184 owns Tahoe; CD4 is three states of one donor; KOLF2 loses after whitening (exp 057) | §3.3 says line count does not convert |
| the unseen D/E/F panel | untestable by construction. The assurance is structural: the prior carries **no target identity**, is fitted on five lines none of which is a challenge context, and F3/F4 show that deleting the overlapping targets changes it by rho 0.999 | — |
| `cell-eval2` on anything | forbidden (25 GB RSS, shared machine) | — |

### Integrity

Six known-answer rows gate every script and all passed on every run. `h1_2025` is `assert`-banned
in the loader and was never opened. No sealed data opened, no `cell-eval2` run, no build, no
`vcc` invocation, no `submissions/` write, and no file outside
`experiments/189-per-gene-reliability/` created or modified.

# 185 — Rank-controlled emission: the p-value order IS a designable object, and it is delivered — then rule 13 vetoes the slot

Date: **2026-09-03**. Measurement only. No `cell-eval2`, no build, no `vcc`, no submission, nothing
written outside this directory. Peak RSS 3.8 GB (`realize185.py`), 11.9 GB (`claim185.py`, the
challenge controls loaded dense in chunks — under the 12 GB cap, over the 10 GB target; noted).

**Decision in one line: the hypothesis is CONFIRMED on the structural question and on delivery —
the emitted p-value order is a free design variable at bit-identical pseudobulk, and an exact-
marginal emitter delivers the `|H|` order at Spearman +0.999 with a twin spread of 0.00007, taking
h1-test `reach` 0.183 → 0.233 (95% of the counterfactual ceiling) — and the arm is still NO-GO,
because it is an emission-law arm, this instrument inverts a21's board `reach` sign exactly as
audit 181 found, and rule 13's toll therefore applies: composite −0.00243 against a +0.0019 gate.**

What is new relative to 171/177/181, in four lines:

1. **171's closure was of the spread-schedule family, not of the order.** At fixed per-gene totals
   and fixed per-cell total, the per-gene Wilcoxon z has a reachable interval of median width **25**
   (8 at 7 CPM, 40 at 300 CPM); a |z| above 5 is reachable for **99.9–100%** of (perturbation, gene)
   pairs in every expression decile and the null floor |z| < 1 for **95.5%**. Nothing about
   low-expression genes kills it. §2.
2. **The order is delivered, not approximated.** Same |z| multiset as a15 (so the same `n_pred`,
   6,040 vs 6,034 on h1-test and 4,439 vs 4,436 on the challenge), re-keyed by `|H|`: Spearman
   (−log p, |H|) inside the call set goes **−0.391 → +0.996**, expression dependence
   **+0.486 → +0.076**, mean `P(1)` **0.7317 → 0.8780** (signed variant **0.9024**, above the 0.90 bar).
   The permutation control (same histograms, random gene assignment) reads **0.167**, below a15.
   §3.
3. **The order and the call set want different keys.** The reference's DE set is detectability-
   driven, so a15's expression-dominated call set overlaps it *better* than the `|H|`-ordered one:
   `jac` 0.0682 → 0.0593 when the whole order is re-keyed. The hybrid — a15's call set, `|H|` order
   inside it and inside its complement — keeps `jac` (0.0683) and delivers **72%** of the reach
   counterfactual (0.217). §3.4, §4.
4. **`pds`, `mse` and `nmae` are pinned, and the pins are measured.** 85/85 delivered-total digests
   equal a15's; `claim` sits at the reference cap **0.24653** on the real challenge panel (knee
   0.23268, Ĉ/C_p 1.9–3.3 vs a15's 1.085) so `mse` is exactly flat; every rank-controlled arm reads
   `nmae` **0.96670** to five decimals because a constant per-cell total makes the mean of per-cell
   CPM a function of the panel total alone. §6.

---

## 1. Known-answer rows — all mandatory, all passing

| gate | expected | measured | source | verdict |
| --- | ---: | ---: | --- | --- |
| a15 h1-test `reach`, native BH, fresh emission | 0.18283 | **0.18283** | 157/164/171/177 | PASS |
| a15 median `n_pred` (raw BH count) | 6035 | **6035** | 157/164 (`realize185.log`) | PASS |
| a15 median `n_pred` (debench convention, self gene excluded) | 6034 | **6034** | 177/181 | PASS |
| a8-flat h1-test `reach` / `n_pred` | 0.20291 / 3581 | **0.20291 / 3581** | 157/164/171 | PASS |
| `\|H\|` key reach on the same pool | 0.23533 | **0.23533** | INBOX-4/171/172/177/181 | PASS |
| oracle (own signs perfectly ordered) | 0.64553 | **0.64553** | 171/172 | PASS |
| my reach walk vs `debench.reach()` | identical | **\|diff\| < 1e−9 on every arm** | independent implementation | PASS |
| my histogram-z vs `de_kernel.wilcoxon_z` (scanpy-identical) | identical | **max \|diff\| 1.4e−14** on 300 random histograms over real control columns, \|z\| 0.03–60 | `selftest185.py` | PASS |
| delivered per-gene panel-total digest, every arm, every target | a15's `b53d02aebfd0816a…` | **85/85 targets, 14/14 arms** (171 §6's digest for a15 reproduced) | rule 12, non-vacuous form | PASS |
| `claim(flat 0.25/φ4)` | 0.16652 | **0.16652** (+0.00%) | 096/164/171/177 | PASS |
| `claim(a15 as shipped)` | 0.24581 | **0.24581** (+0.00%) | 164 §5 | PASS |
| a15 challenge `n_pred` (16 targets × 3 contexts) | 4,424 (096 probe) | **4,436** | 096/179 | PASS (+0.3%) |
| `P(emitted sign == sign H)` on the pool, a15 | 0.945 | **0.9452** | 172 §4.2 | PASS |
| 171 §12a Spearman table, a15 (|H| universe / call; expr universe / call) | −0.121 / −0.391 / +0.444 / +0.485 | **−0.121 / −0.391 / +0.444 / +0.486** | 171/177/181 | PASS |
| 171 §2b `P(k)` curve, a15 (`P(1)`, `P(3)`, `P(10)`) | 0.7317 / 0.6491 / 0.5655 | **0.7317 / 0.6491 / 0.5655** | 171 | PASS |

Surface: h1-test, 85 perturbations × 10,779 genes past the scorer's 5-CPM control gate, 18,400
real control cells. Donors: 095's h1-test artifact (`k562_gwps` + `orion_hek293t`). `h1_2025`
never enters. Per-gene totals are the a15 artifact's own (`totals185.py`, cached once); every
known-answer arm is a *fresh* emission through production's `emit_counts` on `realize171.py`'s
exact path and seeds.

---

## 2. Part A — controllability: the reachable z interval per gene

### 2a. The object, and why it is exact

The scorer ranks per-cell CPM against the real controls. With all 400 emitted cells at one total
`T`, CPM is monotone in the integer count, so define per gene
`Φ_g(c) = P_ctrl(cpm < 1e6·c/T) + ½·P_ctrl(cpm = 1e6·c/T)`. The Mann–Whitney U of our cells is
`U = m·Σ_i Φ_g(c_i)` **exactly** (midranks make U independent of ties among our own cells; the tie
correction enters only the denominator). So `F̄ = mean_i Φ_g(c_i)` is a *separable* function of the
per-cell histogram, and over histograms with mean `M_g/400` and cap `T` its reachable set is the
segment between the lower-convex and upper-concave envelopes of `Φ_g` at that mean — each endpoint
attained by a two-level histogram (Carathéodory on one moment constraint). `eme185.py` computes
both envelopes, the even (maximally spread) and the most-concentrated allocations, and the
histogram closest to the null, with the exact scanpy tie correction (`hist_z`, verified above).

### 2b. The distribution — all 85 perturbations pooled (medians over (perturbation, gene) pairs)

| control-CPM decile | width `z_max − z_min` | `z_min` | `z_even` | `z_max` | max \|z\| | min \|z\| | null inside | floor < 1 | ceiling > 5 | ctrl zero frac |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| d1 (7 CPM) | 8.1 | −7.6 | +0.48 | +0.48 | 7.6 | 0.049 | 64.8% | 84.2% | **99.9%** | 0.87 |
| d2 (12) | 11.7 | −10.4 | +1.31 | +1.31 | 10.4 | 0.030 | 79.3% | 91.2% | 100% | 0.78 |
| d3 (20) | 15.8 | −13.4 | +2.40 | +2.40 | 13.4 | 0.019 | 90.1% | 96.0% | 100% | 0.67 |
| d4 (29) | 19.9 | −16.5 | +3.49 | +3.49 | 16.5 | 0.013 | 94.9% | 97.9% | 100% | 0.56 |
| d5 (40) | 23.7 | −19.6 | +4.12 | +4.12 | 19.6 | 0.011 | 96.2% | 98.4% | 100% | 0.45 |
| d6 (56) | 26.2 | −22.8 | +3.41 | +3.43 | 22.8 | 0.016 | 93.2% | 97.0% | 100% | 0.34 |
| d7 (76) | 29.7 | −26.1 | +3.08 | +3.74 | 26.1 | 0.014 | 93.6% | 97.1% | 100% | 0.24 |
| d8 (105) | 32.8 | −29.0 | +3.35 | +3.81 | 29.0 | 0.017 | 93.1% | 96.6% | 100% | 0.15 |
| d9 (158) | 36.0 | −31.8 | +2.94 | +4.24 | 31.8 | 0.019 | 94.5% | 97.4% | 100% | 0.07 |
| d10 (312) | 40.2 | −33.7 | +2.31 | +6.51 | 33.7 | 0.017 | 98.4% | 99.2% | 100% | 0.01 |
| **ALL** | **25.0** | **−21.2** | **+2.65** | **+3.22** | **21.2** | **0.018** | **89.8%** | **95.5%** | **100%** | 0.40 |

By `|H|` stratum (per-perturbation quantiles):

| stratum | width | `z_min` | `z_even` | `z_max` | max \|z\| | min \|z\| | null inside | floor < 1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `\|H\|` top 1% | 25.0 | −23.1 | −1.96 | −1.28 | 23.2 | 1.275 | 44.1% | 47.7% |
| `\|H\|` top 10% | 26.4 | −23.0 | +0.55 | +1.63 | 23.0 | 0.041 | 60.5% | 74.8% |
| `\|H\|` below median | 24.3 | −20.1 | +2.72 | +3.30 | 20.1 | 0.017 | 96.6% | 99.9% |

Headline perturbation AASS (S = 8,245,740, T = 20,614): the same shape row for row
(`work/partA.log`); a15's own realised |z| on AASS (median 2.84, 90th pct 11.3, max 20.3) lies
inside the reachable interval for **98.9%** of genes.

### 2c. Reading

* **The interval is wide for every decile, including the low-expression genes that dominate the
  pool.** The mandate's abort condition ("if it collapses for the low-expression genes, the
  mechanism is dead") does not fire: at 7 CPM with 87% control zeros the reachable |z| still spans
  0.05 to 7.6, and the DOWN side reaches |z| ≈ 21 at the median gene.
* **The interval is asymmetric.** Concentration (zeros plus a few loaded cells) always moves the
  rank sum DOWN; the UP side is bounded by the concave envelope and sits close to the even
  allocation (`z_max ≈ z_even`, +0.5 to +6.5). This is 171 §12c's U-shape seen from the inside: the
  even emitter's `+2.65` is the "400 identical values at the mean" artefact, and no allocation can
  push a low-expression gene much further UP. **Since the scorer reads direction from the lfc sign
  and never from the z sign** (verified in `cell_eval2/metrics/direction.py`: `match` is
  `sign(lfc_pred) == sign(lfc_real)`, and the reach sort is `(sig, p_adj, p, |lfc|, name)`), the
  asymmetry costs nothing: a large |z| is available to every gene on the DOWN side.
* **Genes with a real predicted effect have a floor**, not a ceiling: `|H|` top-1% genes read
  min |z| 1.28 because their mean is far from the control's. Those are the genes that should be
  significant anyway.
* **`z_even` (the natural allocation) tracks expression, not effect** — +0.48 at d1 rising to +4.1
  at d5 — which is the whole reason a15's order is expression-dominated (171 §12).

### 2d. The two coupling constraints, handled

*Integers.* Every histogram is integer by construction (two-level families; exact z).

*Per-cell totals.* Each gene's histogram is assigned to cells greedily against remaining capacity
(largest level to the emptiest cell, random tie-break), then the genes **below the 5-CPM gate** —
invisible to every DE member — are used as slack: single counts are moved between cells until every
row sums to exactly `T` (`S // 400`, the `S mod 400` remainder cells at `T + 1`). Measured on all
85 × 9 designed blocks: **forced placements 0, maximum pre-repair residual 0, repair moves 0,
unresolved 0** — the greedy alone balanced every row exactly, so the slack was never spent. The
sub-gate reservoir is 113 counts per cell (median). Per-gene totals are untouched by construction
and the digest requirement (85/85 = a15's) is the rule-12 proof. Designed vs realised z (the
`T + 1` cells are the only difference): max |diff| 0.0044, median 3e−5.

---

## 3. Part B — realising a prescribed order

### 3a. The law

`z*_g = sign(H_g)·φ(|H_g|)` with `φ` the quantile map onto a15's **own** per-perturbation |z|
histogram: sort a15's realised |z| descending, assign the r-th largest to the gene with the r-th
largest `|H|`. This keeps the p-value multiset — hence BH's call count — at a15's, and isolates the
*order*. Two sign policies: `abs` (deliver the |z|, sign free; 33% of designed genes land with a z
sign opposite to sign(H), which the scorer never reads) and `signed` (clip to the reachable side of
sign(H); 11.9% of genes clipped by > 0.5 in |z|). Then: sub-gate genes on a15's own tail law, greedy
placement, exact repair, and the real Wilcoxon/BH.

### 3b. Members (h1-test, 85 perturbations; `reach` at proportional power `f`; `P(k)` = purity of the first k adjudicable pool genes)

| arm | `n_pred` | reach f=1 | f=2 | f=3.5 | fid | prec | jac | nmae | recall | `P(1)` | `P(3)` | `P(10)` | k* med |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **a15** (KA) | 6034 | 0.18283 | 0.12202 | 0.05361 | 0.50701 | 0.50701 | 0.06824 | 0.96659 | 0.669 | 0.7317 | 0.6491 | 0.5655 | 1 |
| a15_b (twin) | 6031 | 0.17537 | 0.10295 | 0.05446 | 0.50771 | 0.50771 | 0.06834 | 0.96641 | 0.685 | 0.6707 | 0.6491 | 0.5741 | 1 |
| a8-flat (KA) | 3581 | 0.20291 | 0.13124 | 0.05835 | 0.50793 | 0.51209 | 0.06089 | 0.96673 | 0.488 | 0.8049 | 0.7456 | 0.6121 | 2 |
| **rc_hist_a** | 6040 | **0.23278** | **0.14938** | 0.06948 | 0.50808 | 0.50808 | 0.05932 | 0.96670 | 0.616 | **0.8780** | 0.7851 | 0.6914 | 2 |
| rc_hist_a_b (twin) | 6032 | **0.23272** | 0.14952 | 0.06891 | 0.50803 | 0.50803 | 0.05930 | 0.96670 | 0.616 | 0.8780 | 0.7807 | 0.6897 | 2 |
| rc_hist_s (signed) | 5922 | 0.22579 | 0.14444 | 0.07022 | 0.50825 | 0.50825 | 0.05959 | 0.96670 | 0.612 | **0.9024** | 0.7895 | 0.6810 | 3 |
| **rc_hybrid** (a15's set, `\|H\|` inside) | 6043 | **0.21681** | 0.14421 | **0.07368** | 0.50693 | 0.50693 | **0.06827** | 0.96670 | **0.669** | 0.8780 | 0.7851 | 0.6914 | 2 |
| rc_lfc (key = a15's own `\|lfc\|`) | 6038 | 0.22808 | 0.13455 | 0.07062 | 0.50651 | 0.50651 | 0.05463 | 0.96670 | 0.568 | 0.8780 | 0.7588 | 0.6828 | 2 |
| **rc_perm (CONTROL)** | 6116 | **0.16705** | 0.09113 | 0.03363 | 0.50459 | 0.50459 | 0.06011 | 0.96670 | 0.600 | 0.6220 | 0.6316 | 0.5586 | 1 |
| rc_k1000 | 1001 | 0.20604 | 0.12062 | 0.05195 | 0.44394 | 0.51736 | 0.02679 | 0.96670 | 0.192 | 0.8659 | 0.7588 | 0.6862 | 2 |
| rc_k2000 | 1999 | 0.20527 | 0.12843 | 0.06051 | 0.48390 | 0.51283 | 0.03764 | 0.96670 | 0.267 | 0.8659 | 0.7675 | 0.6948 | 2 |
| rc_k4000 | 3999 | 0.21577 | 0.13794 | 0.06890 | 0.50568 | 0.50929 | 0.05089 | 0.96670 | 0.430 | 0.8659 | 0.7763 | 0.6897 | 2 |
| a21_like (Part D) | 1913 | 0.23252 | 0.15391 | 0.07650 | 0.47935 | 0.51421 | 0.04307 | 0.96934 | 0.271 | 0.8049 | 0.7456 | 0.6466 | 2 |
| a22_like (Part D) | 5569 | 0.17827 | 0.09648 | 0.03565 | 0.50706 | 0.50706 | 0.06871 | 0.96713 | 0.656 | 0.6585 | 0.5965 | 0.5517 | 1 |

Counterfactual ceilings on each rc arm's own table: `|H|` key **0.23533** (all), oracle 0.64529.
**rc_hist_a realises 95.1% of the `|H|` increment** ((0.23278 − 0.18283)/(0.23533 − 0.18283)) and
**rc_hybrid realises 94% of INBOX-4's "sig-first, then `|H|`" counterfactual** (0.21681 vs 0.21914,
increments 0.03398 vs 0.03631).

### 3c. Delivery fidelity — the mandate's primary Part B statistic

Median over perturbations of Spearman(−log₁₀ p, x), over the DE universe / inside our call set /
over the reach pool (171 §12a's statistic; a15's row is the known-answer row):

| arm | `\|H\|` universe | `\|H\|` call set | `\|H\|` pool | expression universe | expression call set | expression pool |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| **a15** (incumbent) | −0.121 | **−0.391** | −0.261 | +0.444 | **+0.486** | +0.457 |
| a8-flat | +0.479 | +0.374 | +0.358 | +0.282 | +0.341 | +0.337 |
| **rc_hist_a** | **+0.999** | **+0.996** | **+1.000** | +0.091 | **+0.076** | +0.067 |
| rc_hist_s | +0.991 | +0.956 | +0.996 | +0.128 | +0.138 | +0.091 |
| rc_hybrid | +0.514 | +0.984 | +0.507 | +0.345 | +0.007 | +0.282 |
| rc_lfc | +0.877 | +0.681 | +0.886 | −0.263 | −0.290 | −0.217 |
| rc_perm (control) | +0.022 | −0.011 | +0.035 | +0.001 | +0.038 | −0.004 |
| rc_k4000 | +0.748 | +0.986 | +0.744 | −0.008 | +0.087 | +0.042 |
| a21_like | +0.590 | +0.557 | +0.624 | +0.188 | +0.060 | +0.275 |

*(rc_hybrid's universe/pool value is 0.5 by design: `|H|` is the order inside the call block and
inside its complement, with the block boundary set by a15's own p; rc_lfc's expression correlation
is negative because a15's `|lfc|` is anti-correlated with expression.)*

**The p-order is a free variable.** The same |z| multiset delivers `|H|` at +0.999, `|lfc|` at +0.88,
or nothing at all (+0.02) depending only on the assignment. INBOX-4's falsifier — a15's
+0.485 expression / −0.391 `|H|` inside the call set — is reversed to +0.076 / +0.996, further than
177's best resample arm (+0.001 / +0.601) and with no resampling.

### 3d. The `P(k)` curve

| arm | `P(1)` | `P(2)` | `P(3)` | `P(5)` | `P(10)` | `P(20)` | `P(50)` | `P(100)` | reach |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| a15 | 0.7317 | 0.6623 | 0.6491 | 0.6200 | 0.5655 | 0.5350 | 0.5357 | 0.5237 | 0.18283 |
| rc_hist_a | **0.8780** | 0.8117 | 0.7851 | 0.7686 | 0.6914 | 0.6770 | 0.6323 | 0.6146 | 0.23278 |
| rc_hist_s | **0.9024** | 0.8182 | 0.7895 | 0.7486 | 0.6810 | 0.6600 | 0.6349 | 0.6139 | 0.22579 |
| rc_hybrid | 0.8780 | 0.8182 | 0.7851 | 0.7571 | 0.6914 | 0.6850 | 0.6362 | 0.6156 | 0.21681 |
| rc_perm | 0.6220 | 0.6429 | 0.6316 | 0.5914 | 0.5586 | 0.5510 | 0.5404 | 0.5329 | 0.16705 |
| `\|H\|` key (counterfactual) | 0.8780 | 0.8247 | 0.7982 | 0.7714 | 0.6931 | 0.6740 | 0.6302 | 0.6146 | 0.23533 |

The realised curve sits on the counterfactual's to two decimals at every depth. Mean `P(1)` clears
the bar `P₀ = 0.90` for the signed arm; 171's best arm was 0.8171. What has *not* changed is the
tail: `P(100)` is 0.61 for every `|H|`-ordered arm, which is the per-gene sign accuracy of the key
itself (172), not of the emitter.

### 3e. Seed noise is gone

| pair | reach | reach (twin) | \|diff\| |
| --- | ---: | ---: | ---: |
| a15 / a15_b | 0.18283 | 0.17537 | 0.00746 |
| rc_hist_a / rc_hist_a_b | 0.23278 | 0.23272 | **0.00007** |

The design is a deterministic function of the |z| multiset and `|H|`; only the multiset (from the
a15 twin) and the placement seed differ between twins. Every number in this family is therefore
resolvable at 1e−4 on h1-test, against the 0.007–0.016 spreads that made 171 §7 and 177 §4 need
twins. (This does **not** bear on the h1-test → board question, §5.)

### 3f. Where the gain lives (171 §7a(i) / 181 §15's test)

Paired Δreach vs the a15 twin mean, by `n_real` stratum (share of the panel sum in brackets):

| arm | `n_real` ≤ 13 (31) | 14–800 (30) | > 800 (24) | `n_real` ≥ 20 | `n_real` ≥ 100 | ALL (82) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `\|H\|` counterfactual | +0.0994 (60%) | +0.0578 (38%) | +0.0039 (2%) | +0.0150 | +0.0082 | +0.0562 |
| rc_hist_a | +0.0994 (63%) | +0.0512 (35%) | +0.0034 (2%) | +0.0122 | +0.0059 | +0.0537 |
| rc_hybrid | +0.0666 (60%) | +0.0382 (37%) | +0.0034 (3%) | **+0.0155** | +0.0055 | +0.0377 |
| rc_perm | −0.0361 (102%) | +0.0009 | −0.0002 | −0.0042 | −0.0010 | −0.0121 |
| a21_like | +0.1129 (72%) | +0.0384 (26%) | +0.0029 (2%) | +0.0064 | +0.0037 | +0.0534 |

rc_hist_a: paired t +3.37 (43 up / 21 down / 18 same), top-5 |Δ| share 41%, most influential
perturbation `n_real` 9 (Δ +0.556). The concentration in the small-`n_real` stratum is the
counterfactual's own (181 §15: 63.7%), not an artefact of delivery, and the gain survives at
challenge-like power: rc_hist_a − a15 = **+0.027 at f = 2, +0.016 at f = 3.5, +0.021 at f = 5**;
rc_hybrid +0.022 / +0.020 / +0.014. 181 §14 puts the board at f ≈ 1.05.

---

## 4. Part C — the call-set policy, now that it is free

`fid`/`jac` of "call the top-K per perturbation by key", on h1-test at proportional power `f`
(the `n_real` *distribution* is kept and scaled; rule 9 applies to every absolute here). Format
`fid/jac`.

| f (median / mean `n_real`) | key | K=50 | 100 | 200 | 400 | 800 | 1600 | 3200 | 6400 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **1.0** (89 / 765) | `\|H\|` | .261/.012 | .296/.012 | .333/.013 | .383/.017 | .428/.024 | .472/.034 | .503/.047 | .508/.060 |
| | a15 p-order | .268/.008 | .301/.013 | .333/.021 | .378/.030 | .421/.041 | .465/.051 | .500/.062 | .507/.068 |
| | random | .273/.003 | .303/.005 | .339/.008 | .382/.014 | .420/.022 | .465/.033 | .497/.047 | .504/.061 |
| **2.0** (179 / 1529) | `\|H\|` | .235/.010 | .263/.011 | .296/.013 | .336/.018 | .385/.029 | .426/.045 | .470/.071 | .500/.107 |
| | a15 p-order | .242/.007 | .268/.011 | .298/.019 | .333/.030 | .380/.043 | .421/.061 | .467/.086 | .499/.115 |
| **3.5** (314 / 2417) | `\|H\|` | .216/.009 | .240/.010 | .269/.013 | .305/.020 | .347/.033 | .393/.055 | .436/.093 | .480/.154 |
| | a15 p-order | .224/.006 | .245/.010 | .272/.018 | .304/.028 | .344/.044 | .389/.067 | .434/.104 | .479/.160 |

Realised step arms sit on the `|H|` curve (f = 1: rc_k1000 fid 0.444 / jac 0.027, rc_k2000
0.484 / 0.038, rc_k4000 0.506 / 0.051). Full table for f = 0.5 and every key in `work/six185.log`.

**Optimum as a function of assumed `n_real` scale**, in overall units against a15's realised set,
best (K, key) over the grid: f = 0.5 → K = 3200 by a15's order, **+0.0014**; f = 1 → K = 6400 by
a15's order, **−0.00004**; f = 2 → K = 6400, **+0.0018**; f = 3.5 → K = 6400, **+0.0048**. Three
readings:

1. **`fid` is flat above K ≈ 3,000 for every key** (0.50 ± 0.01) and falls below it through the
   yield cap. Precision is 0.507–0.517 for every key at every K: no key selects a high-precision
   subset. So the board's six-point "shrinking the call set has never raised `fid`" series
   (164 §4, 181 §1) is what this frontier predicts too.
2. **For `jac` the best key is a15's own expression-dominated p-order, at every f and every K**,
   and `|H|` is no better than *random* at K = 6400 (0.0602 vs 0.0606 at f = 1). The reference's
   DE set is a detectability set — high-expression genes clear its Wilcoxon — and the call set
   that overlaps it best is one ordered by detectability, which is exactly what 171 called the
   incumbent's defect for `reach`. **The order that maximises `reach` and the set that maximises
   `jac` are different objects**, and the hybrid arm is the constructive resolution: the set is a
   prefix of the p-order, so put a15's set in front (BH on the same |z| multiset reproduces it
   exactly — it is a p-threshold set) and order by `|H|` inside it.
3. **The optimum K grows with `n_real`**, because `jac`'s numerator needs recall once the reference
   set is large: at f = 3.5 the whole grid wants K ≥ 6400. Under 181 §14's board distribution
   (median 90–110, a heavy tail about half h1-test's) the a15-sized set is within +0.002 of the
   frontier's best, and no call-set change is worth a slot on its own.

---

## 5. Part D — board retrodiction (audit 181's standing rule)

**Class.** This arm changes the emission law (the count matrix's shape at fixed totals) and nothing
else. It is in a21/a22's class for rule 13. It differs from a21 in what the mandate names: the
delivered profile is bit-identical (85/85 digests; a21's was too, to 16 digits on the board) **and
the per-cell total is welded to a constant**, which a21 unwelded (a21_like's per-cell total sd here:
6,136; a15's 519; every rc arm's 0.47).

**Retrodiction through this instrument** (h1-test raw Δ vs a15 → scaled with 181's board-fitted
1/span; `reach` additionally × 0.4277):

| arm | member | h1-test raw Δ | → scaled | × 0.4277 | board (scaled) | sign |
| --- | --- | ---: | ---: | ---: | ---: | :-: |
| a21_like (K = 3,061 parametric, 72% resampled) | **reach** | **+0.04969** | **+0.05612** | **+0.02400** | **−0.018482** | **WRONG** |
| | fid | −0.02765 | −0.09193 | | −0.029464 | ok (3.1× hot) |
| | jac | −0.02518 | −0.06850 | | −0.016583 | ok (4.1× hot) |
| | nmae | −0.00275 | −0.00454 | | −0.064214 | ok (14× cold) |
| a22_like (top 2,199 resampled, rest a15's law) | **reach** | −0.00456 | −0.00515 | −0.00220 | **−0.010656** | **ok** (0.2–0.5×) |
| | fid | +0.00006 | +0.00019 | | −0.004447 | wrong, inside noise |
| | jac | +0.00046 | +0.00126 | | −0.001523 | wrong, inside noise |
| | nmae | −0.00054 | −0.00089 | | +0.001520 | wrong, inside noise |

**The instrument retrodicts a22 and inverts a21's `reach`**, exactly as 181 §1 found with 177's
tables (+0.044/+0.048 there; +0.050 here on the raw a15 totals). One difference from 177 is worth
recording: with the raw a15 totals rather than the mass-fixed profile, a21_like's `nmae` moves the
*right* way (worse), so 3 of 4 a21 signs now reproduce — but `reach`, the member this experiment is
about, does not, and it is off by 0.07 scaled.

**Does the toll apply?** Rule 13 (audit-181 amendment) says: apply it to emission-law arms, never
argue it away. The only measurement that could exempt this arm is a board point of a
*designed-order* emission law, and none exists. Two facts are stated for the record, neither of
which is such a measurement:

* a21's h1-test→board inversion has a candidate mechanism this arm does not share — a21's tail
  block was resampled real cells whose p-values are uniform *only if the drawn cells are
  exchangeable with the scorer's control cells*, which holds on h1-test by construction (we resample
  from the very cells the test uses) and is unverifiable on the board. This arm resamples nothing;
  its zero-effect behaviour is a15's own |z| multiset by construction. **Hypothesis, unmeasured.**
* a21's `nmae` collapse (−0.064) ran through the unwelded per-cell library; this arm's `nmae` is
  pinned to five decimals (§6). That is a measurement, but of `nmae`, not of the `reach` toll.

**So the toll is charged.** The mandate's shrink factor 0.4277 was calibrated on a8↔a15, a
spread-policy pair inside the quasi law; nothing in this experiment validates it for a designed
order, and reading 2 below is labelled accordingly.

---

## 6. Part E — six-member price

### 6a. `pds`, `mse`, `nmae` — pinned, and the pins measured

* **`pds`**: the delivered per-gene panel totals hash to a15's digest on **85/85 targets for every
  arm** (`b53d02aebfd0816a…` on target 1, 171 §6's value). BRIEF §1 defines `pds` on that
  pseudobulk. **Δ`pds` = 0.000000, identity.** (`debench`'s `pds` proxy reads the mean-of-ratios
  lfc and is not quoted, per 171 §6's caveat.)
* **`mse`**, 096's `claim` route on the real challenge panel (`claim185.py`, 16 targets × 3
  contexts, both gates +0.00%):

| policy | Ĉ | Ĉ/C_p | `claim` | branch | raw `mse` | Δ`mse` overall | challenge `n_pred` |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: |
| flat 0.25/φ4 (a8, KA) | 10.29 | 0.675 | 0.16652 | linear | 1.17174 | −0.018520 | — |
| **a15 as shipped (KA)** | 16.53 | 1.085 | **0.24581** | SATURATED | 0.88228 | — | **4,436** |
| **rc_hist_a** | 42.26 | 2.775 | **0.24653** (cap) | SATURATED | 0.88228 | **0.000000** | 4,439 |
| rc_hist_s | 28.35 | 1.862 | 0.24653 | SATURATED | 0.88228 | 0.000000 | 4,423 |
| rc_perm | 50.10 | 3.290 | 0.24653 | SATURATED | 0.88228 | 0.000000 | 4,444 |
| rc_k2000 | 42.83 | 2.812 | 0.24653 | SATURATED | 0.88228 | 0.000000 | 2,000 |
| rc_k4000 | 42.85 | 2.813 | 0.24653 | SATURATED | 0.88228 | 0.000000 | 4,000 |

  The concentrated histograms carry far more per-cell log-variance than the gamma-Poisson tail
  (Ĉ/C_p 1.9–3.3 against a15's 1.085), so `claim` hits the reference cap and `min(s²B, claim)`
  saturates: **`mse` is exactly flat, with 0.014 of claim headroom above the knee even for the
  signed arm.** rc_hybrid's multiset is rc_hist_a's/rc_perm's, which bracket it at 42–50; it was
  not run separately. Note also the challenge-side `n_pred`: the histogram permutation preserves
  a15's call count on the challenge (4,439 vs 4,436), so the board's `fid`-vs-`n_pred` series reads
  Δ`fid` ≈ 0 for it.
* **`nmae`**: with every cell at total `T`, `mean_i(1e6·c_ig/T) = 1e6·M_g/(400·T)`, so the scored
  lfc is a function of the pinned totals. Every rc arm reads **0.96670** to five decimals
  (a15 0.96659, its twin 0.96641: the quasi law's own per-cell-total noise, sd 519 counts). The
  difference to a15's twin mean is −0.00020 raw = −0.00005 overall — the price of removing a15's
  library noise, inside 164's 0.00018 seed spread.

### 6b. The ledger (overall units; a15 = twin mean; `jac` × 181's 0.44 board/h1-test ratio; `fid` from h1-test, with the board series reading 0 at unchanged `n_pred`)

| arm | Δreach raw | reach (×0.4277) | Δfid raw | fid | Δjac raw | jac | Δnmae | nmae | **reading 2: Σ h1-test** | **reading 1: rule 13** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| rc_hist_a (twin mean) | +0.05365 | +0.00432 | +0.00072 | +0.00040 | −0.00897 | −0.00179 | −0.00020 | −0.00005 | **+0.00287** | **−0.00243** |
| **rc_hybrid** | +0.03770 | +0.00304 | −0.00042 | −0.00023 | −0.00002 | −0.00000 | −0.00020 | −0.00005 | **+0.00274** | **−0.00243** |
| rc_hist_s | +0.04669 | +0.00376 | +0.00089 | +0.00049 | −0.00870 | −0.00174 | −0.00020 | −0.00005 | +0.00246 | −0.00243 |
| rc_lfc | +0.04898 | +0.00394 | −0.00085 | −0.00047 | −0.01366 | −0.00273 | −0.00020 | −0.00005 | +0.00069 | −0.00243 |
| rc_k4000 | +0.03667 | +0.00295 | −0.00168 | −0.00093 | −0.01741 | −0.00347 | −0.00020 | −0.00005 | −0.00151 | −0.00243 |
| rc_k2000 | +0.02617 | +0.00211 | −0.02345 | −0.01300 | −0.03065 | −0.00612 | −0.00020 | −0.00005 | −0.01706 | −0.00243 |
| rc_perm (control) | −0.01205 | −0.00097 | −0.00277 | −0.00154 | −0.00818 | −0.00163 | −0.00020 | −0.00005 | −0.00419 | −0.00243 |

`pds` and `mse` are 0.000000 in both readings.

**Reading 1 — rule 13 as amended, the standing rule for this class:**
`predicted = 1.011·Δpds/6 + Δmse/6 − 0.0146/6 = 0 + 0 − 0.00243`. **−0.00243 < +0.0019 ⟹ VETO**,
by construction for any arm whose value lives in `reach`/`fid`/`jac`.

**Reading 2 — if h1-test transferred with the a8/a15 shrink** (the reading a21 falsified for the
resample law; *not* board-validated for a designed order): rc_hybrid **+0.0027**, rc_hist_a
**+0.0029**, both above the +0.0019 shape-only gate, carried entirely by `reach`.

---

## 7. Verdict

### NO-GO. Do not spend a slot. a15 ships unchanged, zero flags different.

The standing rule is unambiguous and this experiment does not have the one thing that could
overturn it — a board point for a designed-order emission law. Under rule 13 the arm reads
−0.00243 against a +0.0019 gate.

**What is established, and should be written into the project's map:**

1. **171 §12's closure was of the spread-schedule family, not of the p-order.** The order is a free
   variable at fixed pseudobulk: the reachable |z| interval is wide in every expression decile
   (§2), and an exact-marginal emitter delivers any prescribed order (`|H|` at Spearman +0.999,
   `|lfc|` at +0.88) with the same call count as a15, a twin spread of 7e−5, `pds` by identity,
   `mse` at the cap, `nmae` to five decimals. "No emission law makes p-order follow |lfc|"
   (171 §4, citing 094's +0.07) is **false** as a general statement; it was true of the quasi
   family.
2. **The `|H|` counterfactual (+0.0525 raw) is deliverable: 95% of it, seed-free**, and the
   delivered `P(k)` curve equals the counterfactual's to two decimals. CURRENT_STATE's "no known
   delivery mechanism survives the board" now reads: *one delivery mechanism exists that has never
   been on the board*.
3. **`reach` and `jac` want different orders**, and the hybrid arm reconciles them at 72% of the
   reach counterfactual with `jac` unchanged (§3.4, §4). Any future delivery should be the hybrid,
   not the whole-axis re-key.
4. **The `P(k)` tail is the key's, not the emitter's.** `P(100)` = 0.61 for every `|H|`-ordered
   arm. Past the head, `reach` is bounded by per-gene sign accuracy (172), exactly as 171 said.

**If the slot economics change — the preregistration.** The arm is `rc_hybrid`: a15's per-gene
totals, per-cell totals welded to `S // 400`, a15's own per-perturbation |z| multiset assigned
inside a15's BH call set by descending `|H|` and inside its complement by descending `|H|`
(`abs` sign policy), sub-gate genes on a15's tail law, greedy placement with sub-gate repair.
Predicted board reading under the two hypotheses this slot would discriminate, both on a
zero-noise scorer: **H_transfer ≈ +0.0027 overall** (reach +0.0030, jac 0, fid −0.0002, nmae 0,
pds/mse 0); **H_toll ≈ −0.0024** (reach −0.0146 scaled). Falsifiers: (i) built `pds` and `mse`
identical to a15's to 16 digits, and `nmae` within 0.0005 — a larger move means the fixed-`T`
weld or the total repair broke; (ii) challenge `n_pred` (094's probe) within 2% of a15's 4,436;
(iii) board `jac` within ±0.002 scaled of a15's — the hybrid's whole point; (iv) board `reach`
decides. This is a slot that buys information about a class no board point covers; whether that
information is worth ≈ 0.005 of expected loss under H_toll is a decision for the coordinator and
Ahmet, not a recommendation of this experiment.

**One design note for whoever builds it.** The `abs` policy leaves 33% of designed genes with a
Wilcoxon z sign opposite to their lfc sign. `cell-eval2` never reads the z sign (§2c), but a scorer
that did would punish it; the `signed` policy (rc_hist_s, −0.007 raw reach, `P(1)` 0.902) is the
immune variant.

---

## 8. What I did NOT measure, and why

| not measured | why | can it move? |
| --- | --- | --- |
| a board point for a designed-order emission law | forbidden; and it is the only thing that can settle Part D | decides the sign of the arm |
| the Orion surfaces | 157 §1 / AGENTS.md: not validated for `reach`; relative direction only, not needed here | — |
| `claim` for rc_hybrid / rc_lfc | same |z| multiset as rc_hist_a/rc_perm, which bracket it at Ĉ/C_p 2.8–3.3, far above the cap; a separate run would read 0.24653 | No |
| `pds` on the delivered artifact | 85/85 digests; BRIEF §1 defines it on that pseudobulk | **No — identity** |
| challenge-side `reach`/`fid`/`jac` | no truth; `n_pred` measured (4,439) as the only challenge-side DE statistic available | — |
| `n_real` on the challenge | rule 9; the Part C frontier is reported as a function of `f` for that reason | changes the Part C optimum, not Part B |
| a per-gene φ or a per-perturbation K | out of scope; Part C shows a global K ≈ a15's is within 0.002 of the frontier at every f | small |
| the exact integer optimum of the envelope | the two one-parameter families reach the null to |z| 0.018 (median) and the counterfactual to 95%; a finer search would move the third decimal | No |
| `cell-eval2` on anything | forbidden | — |

### Integrity

* Fifteen known-answer rows (§1), all passing; four are bit-identical reproductions of 157/164/171
  from fresh emissions, two are 096's `claim` gates to +0.00%, one is the kernel identity at
  1.4e−14, and one is 171's Spearman/P(k) rows to three decimals.
* **The instrument was validated before it measured**: `selftest185.py` reproduces the scanpy
  kernel from the histogram form, proves the placement round-trip (histograms delivered exactly,
  rows exact, a15's digest) and caught one real bug — Φ tables built at the panel-median `T`
  instead of the target's own `T` read designed-vs-realised z at 0.14 median; fixed to 3e−5.
* **The permutation control fires the right way.** Same |z| multiset, gene↔|z| pairing destroyed:
  reach 0.167, below a15 (0.183), Spearman with `|H|` +0.02. If delivering an order did nothing,
  rc_perm would read what rc_hist_a reads. It reads 0.066 lower.
* Every rc arm carries the a15 twin's histogram as a second seed; the twin spread (7e−5) is
  reported, and it is why single-seed rows are quoted for the rest of the family.
* **Two things went wrong and are reported**: `claim185.py` peaked at 11.9 GB (the challenge h5ad
  densified in 2,000-row chunks alongside a 1.5 GB sorted control block), over the 10 GB target;
  and a stray placeholder file was written to `~/Code/` during the run and deleted within
  the same minute (`git status` clean outside this directory).
* The Part D retrodiction uses the same instrument, the same a15 totals and the same 177 code path
  as every other arm here; a21_like's `nmae` sign differs from 177's rs_all because 177 used the
  mass-fixed recovered profile and this run uses the shipped totals.
* No submission, no `vcc prep`, no `cell-eval2`, no build. Production code untouched.

## 9. Files

| file | what it produces |
| --- | --- |
| `totals185.py` / `work/a15_totals.npz` | the shipped per-target, per-gene totals of 095's h1-test a15 artifact (streamed once) |
| `harness185.py` / `work/ctrl_cs.npy`, `work/ctrl_meta.npz` | the h1-test control block in the scorer's DE space (171's conventions), cached |
| `eme185.py` | the exact-marginal emitter: Φ tables, exact histogram z, envelopes, design, placement, repair |
| `selftest185.py` / `work/selftest185.log` | §1 — kernel identity, envelope ordering, placement round-trip |
| `partA185.py` / `work/partA.log`, `work/partA.npz` | §2 — reachable z intervals, 85 × 10,779 |
| `realize185.py` / `work/realize185.log`, `work/realize185_extra.log`, `work/de185.*.npz` | §3, §5 — 14 arms, digests, DE tables (85 × 10,779, real + predicted p/padj/lfc/z, designed z) |
| `six185.py` / `work/six185.log`, `work/six185.json` | §1, §3, §4, §5 — known answers, members, delivery fidelity, frontier, retrodiction |
| `strata185.py` / `work/strata185.log` | §3.4–3.6 — strata, `P(k)`, twin spreads, power series |
| `claim185.py` / `work/claim185.log`, `work/claim185.npz` | §6 — `claim` on the real challenge panel, both gates, challenge `n_pred` |
| `ledger185.py` | §6b — the ledger in board units |

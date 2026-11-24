# 521 — Per-target DE call-set sizing, reopened by the strength key: NEGATIVE, now closed by its payoff where it matters (cancer lines) and not only by its input

2026-09-28, research worker (Claude Opus 5.5). Analysis only. No build, no `vcc`, no `cell-eval2`, no git. Nothing was
written outside this directory. Every process peaked ≤ 0.97 GB (`work/*.time`, `/usr/bin/time -l`), so no step needed the
33 GiB gate. The coordinator's `nreal521.py` is kept unchanged. Pricing surface: 435's production-emitter H1 stats (the
shipped DE construction: `rc_emit_k`, `CallSet(400, 2200, 1.9)`, clean4 q80 block prior, coarse MDE) on 3 emission seeds.
Predictor truth: H1 (323's reference) and 480's three cancer-fold gates.

## Answer

1. **The strength key is a real n_real predictor on H1, and much weaker on cancer lines.**
   - H1, leave-one-target-out: K562-GWPS alone r +0.60 (R² 0.35). The 4-source mean rank reaches r +0.55 … +0.61 (R² 0.29 … 0.36),
     against 192's +0.24 (R² 0.06). The residual is still ×8–9 in log space.
   - This does not depend on truth power: r is +0.54 … +0.62 against H1's own truth rebuilt at 250, 150 and 100 cells and at fold depth.
   - On the A/B/C panel targets in held-out cancer lines (480's gates), with every same-line source excluded, the key reaches only
     r +0.43 (K562), +0.25 … +0.27 (HCT116), +0.19 … +0.25 (HEK293T); R² 0.02 … 0.18.
2. **Ceiling against TODAY's shipped call set (H1, oracle n_real):** per-target best K +0.042 (f 1) … +0.084 (f 0.35).
   - Shrinking only the negligible targets (n_real < 20, known) is worth +0.006 … +0.010.
   - Calling nothing on the truly null targets is worth +0.010 at f 1: jac's empty-union rule scores those rows 1.
3. **The asymmetric risk is 8.5 : 1.** Shrinking one target to 10 calls gains +0.0004 overall on average if it is
   negligible, and loses −0.0034 if it is not.
   - A shrink therefore pays only if **89–97 %** of the shrunk targets are negligible.
   - H1's key reaches that purity in its extreme tail: 0 of its bottom 8 targets fall in the truth's top half.
   - The cancer folds do not: 7–46 % of their bottom 5–10 % fall in the truth's top half. Priced through H1's own utility
     curve, every small-K rank policy (≤ 25 calls) is **−0.0005 … −0.016** on every fold. The mildest shrinks (→ 400) range
     −0.007 … +0.0005.
4. **Policies with the real key on H1 (overall units, s0/s1/s2, seed-0 target CI):**
   - The LOTO Bayes action learns to do almost nothing at f 1 (−0.0000).
   - Simple rank policies run −0.005 … +0.004 in variant b, and every CI spans 0. Two variant-a rows (bottom 20–30 % → 400) exclude 0
     at f 0.35 (+0.006 [+0.0001, +0.013]) but lose at f 1 (−0.0025 / −0.0047).
   - The best of the 384-cell grid is "bottom 5 % by the threshold-averaged key → 10 calls": +0.0051 [+0.0004, +0.0118] at f 1,
     +0.0062 / +0.0067 at f 0.5 / 0.35. It clears the selection null (p < 0.005), but it is 4 targets (top-3 share 0.82–0.94)
     and a 10-call fid lottery: its expected value on H1's utility curve is +0.0003.
5. **Most of the positive rows at f 0.35 and in the compressed regime are the GLOBAL call-set size, not per-target information.**
   - Bayes(real key) − Bayes(permuted key) is −0.0015 … +0.0028 across regimes, seeds and variants.
   - Global K 1,200–1,600 prices ≈ +0.005 … +0.006 at f 0.35 and in the compressed regime (CIs exclude 0), about 0 at f 0.5, and −0.015 at f 1.
   - The sign at the board therefore depends on the board's n_real tail. 192's fid-coverage fit (f 0.5, proportional) says about 0;
     the board's jac level (≈ H1 at f 0.33) says positive.

**Verdict.** No per-target policy clears the gate once selection and transfer are priced: +0.003 at f 0.35–0.5 with a CI
excluding 0, not worse at f 1, on the lines the board uses. Per-target call-set sizing stays closed. 192 closed it by its
input (r 0.24). It now also fails on payoff: an r 0.55–0.60 key is not pure enough in its tail on cancer lines to carry an
8.5 : 1 asymmetric loss. **No build spec.** By-products (§8): the global call-set size is regime-sensitive under today's
emission; the board's weak end is thinner than proportional scaling assumes; "unexpressed target → null" is refuted.

---

## 1. Known answers (all pass; `ka521.py`, `check521.py`, `nreal521.py`)

| row | expected | measured |
|---|---|---|
| KA1 192's sizing table at f 1.0 with 192's own `frontier192` functions, FIT key, all 12 rows (ORACLE and PREDICTED × m 0.5 … 5) | 192 log, 5 decimals | all 12 PASS; **ORACLE m 2.0 +0.03608**, **PREDICTED m 2.0 −0.09717**; 192 predictor LOTO r +0.2401 |
| KA2 513's four H1 DE members from 435's stats | fid 0.5089291277653674, jac 0.0746985342012533, nmae 0.9427040517722616, reach 0.2692569775927277 | equal to ≤ 2e-16 |
| KA3 the coordinator's `nreal521.py` | K562 r +0.601, mean rank +0.547, AUC 0.746, bottom 20 %: 14/17 negligible | reproduced exactly |
| KA4 reconstructed call set: head = top 400 by \|key\|, block = top 2,200 by control CPM | prescribed 2,600 per target | 2,600 on every target. Realized calls 2,566 median: 399/400 head, 2,165/2,200 block, 7 unsuppressible/target (mean) |
| KA5 four members with the realized order, seeds 0/1/2 | 435 metrics json | exact on all 3 seeds |
| KA6 model order (calls in realized order, then \|key\| + jitter: the emitter's complement law) vs realized reach | realized | −0.0005 / +0.0033 / +0.0040 (2–4 targets differ). Every policy is priced model order vs model order |

## 2. The predictor (`pred521.py`, `thin521.py`, `expr521.py`)

Statistic: 518's `strength_for` (noise-normalised count of |z| > thr·noise genes over the bank's A/B/C-panel MAD scale),
thr 3/4/5. Combinations use the mean centred rank over the banks that carry the target. LOTO = leave-one-target-out linear
calibration of log1p n_real. Coverage = the fraction of targets inside the LOTO Gaussian predictive interval.

**H1 (85 targets, A/B/C power; h1_2025 never a source)**

| feature | n | r | ρ | LOTO sd (const) | ×err | LOTO R² | AUC(<20) | cover 50/80/90 |
|---|---:|---:|---:|---|---:|---:|---:|---|
| 192's design predictor (KA1) | 82 | +0.24 | +0.23 | 2.47 (2.51) | 11.8 | 0.06 | — | — |
| K562-GWPS thr4 | 82 | +0.601 | +0.597 | 2.12 (2.63) | 8.3 | 0.35 | 0.776 | 0.52/0.79/0.91 |
| Orion HCT116 thr4 | 84 | +0.382 | +0.383 | 2.49 | 12.0 | 0.13 | 0.666 | 0.37/0.77/0.93 |
| Orion HEK293T thr4 | 84 | +0.274 | +0.251 | 2.58 | 13.2 | 0.06 | 0.618 | 0.37/0.75/0.96 |
| VIPerturb K562 thr4 | 80 | +0.312 | +0.313 | 2.58 | 13.2 | 0.08 | 0.661 | 0.38/0.78/0.95 |
| S3 = K562 + HCT116 + HEK293T (nreal521) thr4 | 85 | +0.547 | +0.537 | 2.24 (2.66) | 9.4 | 0.29 | 0.746 | 0.41/0.74/0.91 |
| **S4 = S3 + VIPerturb thr4 (primary)** | 85 | **+0.551** | +0.545 | 2.23 | 9.3 | 0.29 | 0.762 | 0.46/0.73/0.94 |
| S4, thresholds 3/4/5 averaged | 85 | +0.588 | +0.579 | 2.17 | 8.7 | 0.34 | 0.775 | 0.46/0.72/0.92 |
| S4 with K562 ×2 | 85 | +0.606 | +0.599 | 2.13 | 8.4 | 0.36 | 0.783 | 0.45/0.73/0.92 |
| thr 3 / thr 5 (S4) | 85 | +0.511 / +0.445 | | | | 0.25 / 0.18 | 0.748 / 0.706 | |
| essential screens (RPE1, K562-ess, Jurkat, HepG2; 18–21 targets) | 18–21 | +0.05 … +0.72 | | | | | 0.44 … 0.88 | |
| S4 + essential screens | 85 | +0.536 | | | | 0.28 | 0.751 | |
| **non-essential H1 targets only (67; the panel is essential-free)**: S4 / S4 thr-avg / K562 | 67 / 67 / 64 | +0.540 / +0.558 / +0.510 | | 2.11 / 2.08 / 2.15 | 8.2 / 8.0 / 8.6 | 0.28 / 0.30 / 0.24 | 0.753 / 0.772 / 0.739 | |
| the target's own control CPM in H1 | 85 | −0.005 | | 2.69 | 14.8 | −0.03 | 0.502 | |

- **Power control:** the same keys against H1's truth rebuilt at 250 / 150 / 100 cells and at 150 cells + 0.6 thinning (the
  folds' depth; 480's gates) give S4 r +0.555 / +0.543 / +0.539 / +0.549 and K562 r +0.598 / +0.597 / +0.593 / +0.621.
  Predictability is power-invariant. The folds' weaker numbers below are a line/panel effect, not power.
- **Calibration:** the LOTO Gaussian predictive distribution is roughly calibrated (80 % intervals cover 0.72–0.79 on H1,
  0.80–0.82 on the folds) and homoscedastic (residual sd lower/upper half 2.2/2.3).

**Cancer folds (A/B/C panel targets in held-out lines; 150 cells, ~12k UMI; sources exclude the held line and same-line banks)**

| fold (n, median n_real) | best single source | production sources | all sources (+ KOLF2, CD4 robust, CD4 DESeq2 nDE) | cross-study only | KOLF2 + CD4 alone |
|---|---|---|---|---|---|
| K562 (267, 3; 64 nulls) | CD4 nDE r +0.33 | r +0.38, R² 0.14, AUC 0.68 | **r +0.43, R² 0.18, AUC 0.71** | = all | r +0.32 |
| Orion HCT116 (132, 68) | CD4 nDE +0.29 / K562 +0.27 | r +0.27, R² 0.06, AUC 0.66 | r +0.25, R² 0.05 | r +0.19, R² 0.02 | r +0.17 |
| Orion HEK293T (132, 141) | K562 +0.22 | r +0.25, R² 0.05, AUC 0.67 | r +0.19, R² 0.02 | r +0.18, R² 0.02 | r +0.07 |

The Orion thinned truths carry a target-independent technical blob (AGENTS §6.5). Even their unexpressed targets show
n_real 4–810, so their low end is noisy. The K562 fold is the cleanest cancer-line test: cross-study, no blob.

## 3. The ceiling against today's call set (`ceiling521.py`; H1, s0, variant b; board overall units)

| policy | f 1.0 | f 0.5 | f 0.35 | cmp |
|---|---:|---:|---:|---:|
| per-target oracle, best K from the grid | +0.042 [+0.027, +0.062] | +0.078 | +0.084 | +0.046 |
| oracle n_pred = 2 n_real (headcpm content) | +0.028 | +0.075 | +0.082 | — |
| negligible targets (n_real < 20, known) → 25 calls | +0.006 [−0.003, +0.014] | +0.009 | +0.010 | +0.007 |
| null targets (n_real = 0, known) → no calls (jac's empty union = 1) | +0.010 [−0.001, +0.027] | +0.038 | +0.038 | +0.010 |

Content when a target shrinks to K is 'headcpm', `CallSet(min(K, 400), max(K − 400, 0))`. That is 192 §5.3's own
construction: the |H| head is kept first, and block genes leave from the lowest CPM up. It is deliverable per target, and
the ceiling rows were best with it. 'keyord' and 'cpmord' are in `ceiling521.log`: dropping the head costs reach −0.005 …
−0.02.

## 4. Policies (`table521.py`, `policy521.py`, `summary521.py`)

Every per-target policy is priced exactly from a per-target row table. Each member is a mean of per-target rows, and a
target's row depends only on its own call set and lfc. Four regimes are used:
- f 1.0;
- f 0.5 / 0.35 (PowerSurface, proportional);
- **cmp**: n' = 2.47 n^0.745. It matches both board facts: the tail at f 0.425 (192 §4.2), and the board's nmae cohort
  209–261/300, i.e. 22 % of targets with a gate < 10. H1 has 32 % at f 1 and 42 % at f 0.35–0.5.

Two emitter variants are priced:
- **b**: dropped block genes keep their clean4 steering, so nmae is pinned. This needs a decoupled emitter.
- **a**: the production emitter with a smaller CallSet. Dropped block genes take the unsteered lfc (416 prior_off) plus the
  MDE add.

Base members (s0, model order):
- f 1: fid 0.50893, jac 0.07470, reach 0.26871, nmae 0.94270;
- f 0.5: 0.52481 / 0.05454 / 0.28020 / 0.94634;
- f 0.35: 0.52487 / 0.04402 / 0.31565 / 0.94702;
- cmp: 0.52487 / 0.04747 / 0.19166 / 0.92884.

**Variant b, value on seeds 0 / 1 / 2 [seed-0 95 % target CI], top-3 share of the gain (nan = no net gain)**

| policy (S4 key unless named) | f 1.0 | f 0.5 | f 0.35 | cmp |
|---|---|---|---|---|
| bottom 10 % → 25 calls | +0.0000/−0.0004/+0.0003 [−0.0053, +0.0055] | +0.0010/+0.0002/+0.0013 [−0.0042, +0.0063] | +0.0016/+0.0015/+0.0019 [−0.0039, +0.0072] | +0.0007/+0.0004/+0.0009 [−0.0046, +0.0062] |
| bottom 20 % → 400 | −0.0036/−0.0039/−0.0035 [−0.0111, +0.0016] | +0.0016/+0.0008/+0.0016 [−0.0036, +0.0067] t3 1.99 | +0.0037/+0.0038/+0.0038 [−0.0008, +0.0097] t3 1.12 | +0.0022/+0.0018/+0.0022 [−0.0014, +0.0061] |
| bottom 30 % → 400 | −0.0044/−0.0047/−0.0043 | +0.0016/+0.0008/+0.0017 | +0.0039/+0.0039/+0.0040 [−0.0007, +0.0098] | +0.0023/+0.0018/+0.0023 |
| two-level 10 % → 25, next 20 % → 400 | −0.0054/−0.0059/−0.0050 | +0.0007/−0.0002/+0.0011 | +0.0030/+0.0029/+0.0033 | +0.0012/+0.0007/+0.0014 |
| **Bayes action (LOTO), S4** | −0.0000/−0.0001/−0.0001 [−0.0010, +0.0008] | +0.0014/+0.0003/+0.0010 [−0.0019, +0.0056] | +0.0045/+0.0045/+0.0044 [+0.0009, +0.0095] (85/85 shrunk) | +0.0057/+0.0057/+0.0059 [+0.0025, +0.0091] (85/85 shrunk) |
| Bayes, S3 / K562 alone | −0.0000 / −0.0016 | −0.0008 / +0.0019 | +0.0045 / +0.0041 | +0.0049 / +0.0061 |
| grid best (384 cells): S4avg bottom 5 % → 10 | **+0.0051 [+0.0004, +0.0118]**, top3 0.93 | +0.0062 [+0.0008, +0.0132] | +0.0067 [+0.0012, +0.0140] | +0.0035 |

- Variant a (`summary521.log`) is within ±0.002 of b on the rank rows. Its Bayes rows are +0.0015 (f 1), −0.0017 (f 0.5),
  +0.0014 (f 0.35) and +0.0062 (cmp). Its nmae moves −0.005 … +0.002 raw (oracles included).
- Across the three emission seeds each row moves by ≤ 0.0008.
- **Degraded predictor** (Gaussian rank noise to the folds' r, 200 draws): at r 0.40 and at r 0.25 every rank policy loses
  at f 1 and f 0.5 (e.g. bottom 20 % → 400: −0.0056 / −0.0097 at f 1, P(≤ 0) 0.97 / 0.99). The Bayes action falls back to
  doing nothing (f 1: −0.0002 / 0.0000), or at f 0.35 and cmp to a global shrink (+0.003 … +0.006 regardless of r).

## 5. Separating the global size from per-target information (`global521.py`)

**Global K: every target at `CallSet(400, K − 400)`.** The table gives each value on seeds 0/1/2, the seed-0 CI, the
members (s0, raw), and the board columns `[reach ×1 | ×0 | ×−0.76]` with fid and jac at 0.72×.

| regime | K 1,200 | K 1,600 | K 2,000 | K 2,200 |
|---|---|---|---|---|
| f 1.0 (10 targets > 2,600) | −0.0277 | −0.0155 [−0.025, −0.007] | −0.0079 | −0.0050 [−0.008, −0.002] |
| f 0.5 (1 target > 2,600) | −0.0066 | −0.0003 [−0.007, +0.006] | +0.0007 [−0.002, +0.003] | +0.0006 |
| f 0.35 (0 targets > 2,600) | +0.0032 | **+0.0048 [+0.0009, +0.0100]**: fid +0.0036, jac +0.0029, reach +0.0077; board [+0.0039 \| +0.0024 \| +0.0013] | +0.0029 [+0.0019, +0.0041] | +0.0017 |
| cmp | **+0.0064 [+0.0031, +0.0098]**: fid +0.0059, jac +0.0047, reach +0.0053; board [+0.0049 \| +0.0039 \| +0.0031] | +0.0061 | +0.0035 | +0.0022 |

**The key's increment over the best global action, Bayes(real) − Bayes(permuted key):**

| variant | f 1.0 (S4 / S3 / S4avg) | f 0.5 | f 0.35 | cmp |
|---|---|---|---|---|
| b | −0.00004 / −0.00001 / −0.00003 | +0.0019 / −0.0005 / +0.0028 | +0.0017 / +0.0018 / +0.0009 | +0.0002 / −0.0005 / +0.0002 |
| a | +0.0015 / +0.0024 / +0.0024 | −0.0002 / −0.0015 / +0.0003 | −0.0002 / −0.0003 / +0.0005 | −0.0002 / −0.0001 / +0.0005 |

Rank policies were also layered on the best global K* (K* = 2,200 / 2,000 / 1,600 / 1,200 by regime):
- "S4 bottom 10 % → 10 calls" adds +0.0036 / +0.0057 / +0.0061 / −0.0000 by regime (variant b). Every CI spans 0 (e.g.
  f 1: [−0.0035, +0.0115]).
- Every larger or wider rank policy is ≤ +0.0017 or negative.

## 6. Selection and transfer (`select521.py`, `diag521.py`, `purity521.py`, `transfer521.py`, `transfer521b.py`)

- **Selection null.** Take the maximum over the 384-cell rank grid (4 keys × 8 q × 12 K) and compare it with the same
  maximum under 200 joint label permutations.
  - Real vs the null's 95th percentile: f 1 +0.0051 vs +0.0019; f 0.5 +0.0074 vs +0.0038; f 0.35 +0.0075 vs +0.0055
    (p ≤ 0.01). In cmp: +0.0047 vs +0.0048 (p 0.07–0.11).
  - So on H1 the key's extreme tail carries real information.
  - At f 0.35–0.5 the winning cells are "bottom 10 % → no calls". Their value comes from proportional scaling turning
    n_real 1 into 0 (8 nulls, 9 % of targets). The board's reach cohort (290–300) allows at most 3.3 %.
- **Who carries it** (`diag521.log`). The bottom 8 by strength are all negligible (n_real 1, 1, 8, 13, 6, 10, 8, 1). Yet at
  10 calls, fid becomes a 10-gene sign lottery: CNN3 (n_real 1) −0.0020, APOE −0.0006, TOR1A +0.0022, MFSD10 +0.0019. That
  is why K 10 vs K 25 swings +0.0036 → 0.0000. The three truly null targets rank 14th, 24th and 25th by strength, so the key
  cannot reach the null-oracle pot.
- **Break-even purity** (H1 per-target utilities, negligible = the truth's bottom 35 %):

  | regime, K | E[u \| negligible] | E[u \| not] | break-even purity |
  |---|---:|---:|---:|
  | f 1, K 10 | +0.00041 (sd 0.0010) | −0.00343 | 0.89 |
  | f 1, K 25 | +0.00010 | −0.00307 | 0.97 |
  | cmp, K 10 | +0.00021 | −0.00326 | 0.94 |
  | cmp, K 25 | +0.00017 | −0.00278 | 0.94 |
- **Tail purity**, as the share of the bottom 5 % / 10 % by strength in the truth's bottom 35 % | in its top half:

  | surface | bottom 5 % | bottom 10 % |
  |---|---|---|
  | H1 S4 | 0.75 \| 0.00 | 0.88 \| 0.00 |
  | K562 fold, all sources | 0.92 \| 0.08 | 0.85 \| 0.07 |
  | K562 fold, production sources | 0.54 \| 0.38 | 0.62 \| 0.27 |
  | HCT116 fold, all sources | 0.71 \| 0.14 | 0.62 \| 0.23 |
  | HEK293T fold, all sources | 0.57 \| 0.43 | 0.62 \| 0.31 |
- **Transfer price.** Each fold's selected targets are read at their truth percentiles on H1's LOTO utility-vs-percentile
  curve. H1's own q10 / K10 row comes out at +0.0010 on the curve vs +0.0036 realized, the gap being the 4-target lottery.

  | policy | K562 fold | HCT116 fold | HEK293T fold |
  |---|---|---|---|
  | f 1, bottom 5 % → 10 | −0.0005 … −0.0068 | −0.0021 … −0.0043 | −0.0048 … −0.0065 |
  | f 1, bottom 10 % → 10 | −0.0021 … −0.0096 | −0.0067 … −0.0159 | −0.0106 … −0.0117 |
  | f 1, → 400 | −0.0009 … −0.0043 | +0.0001 … −0.0046 | −0.0032 … −0.0071 |
  | cmp, → 400 | −0.0001 … −0.0017 | +0.0005 … −0.0013 | −0.0011 … −0.0029 |
- **Moderate shrinks** (bottom 20–50 % → 800–1,600 calls). The key does beat a random selection of the same size: at f 1,
  bottom 50 % → 1,600 gives −0.003 … −0.006 on the folds vs −0.0075 random. But the absolute value is ≤ 0 at f 1 and
  −0.002 … +0.0014 at f 0.5. At f 0.35 and in cmp, selected ≈ random (+0.0021 … +0.0033): that is the global channel.

## 7. Pins, deliverability, and what was not measured

- **pds is pinned by construction.** The per-gene totals are pass-1's and are independent of the call set;
  `rc_hybrid_block_k` asserts them bit-identical.
- **mse is pinned by the same argument that held for the zs rider.**
  - A call-set change only redesigns pass-2 histograms. Variant a also re-splits the dropped block genes across library
    classes.
  - The jackknife credit saturates at the truth's (jk_boost). Two board reads confirm that histogram-only changes leave mse
    bit-identical: rc_callset 09-03 (a 42 % call-set cut) and 446bj2z 09-28 (zs).
  - This is an analogy, not a measurement for this arm.
- **nmae is pinned in variant b.** It moves −0.005 … +0.002 raw in variant a.
- **Deliverability.** 'headcpm' shrinks are native `CallSet` shapes. A per-target size needs a builder patch, because
  `RC_CALLSET_HEAD` and `RC_CALLSET_CPM` are module constants. Variant b also needs the emitter's steering block decoupled
  from the call set; in `rc_emit_k` both come from `call_set.n_cpm`.
- **Not measured:**
  - a real emission of any shrink. Unsuppressible calls among dropped genes would blunt both gains and losses. It needs a
    ~7 GB harness run, so after 02:00Z; not run because no candidate needed it;
  - the q80 threshold and the kf mass balance re-solved over a smaller steered block in variant a (second order);
  - K > 2,600, which needs a model for newly added calls.
- **Units.** Everything here is H1 raw converted to board overall with 192's 1/span. The board columns in §5 apply
  `docs/H1-BOARD-CALIBRATION.md`'s call-set rule: fid and jac at 0.72×, and reach at ×1 (rc_callset), ×0, or ×−0.76 (v14).

## 8. By-products for the coordinator

1. **The global call-set size is regime-sensitive under today's emission.**
   - K 2,600 is right at f 1 and flat at f 0.5.
   - K 1,200–1,600 is +0.005 … +0.006 (H1) at f 0.35 and in the compressed surrogate, i.e. board +0.0013 … +0.0049.
   - The two board-based regime estimates disagree on exactly this: 192's fid-coverage fit (tail at f 0.5, tail-capped
     shapes rejected) against the jac level (today's board jac 0.0408 ≈ H1 at f 0.33).
   - Not recommended as a rider. As an information slot (DE-only, pds/mse pinned), a K 1,600 arm would read the board's
     tail directly through fid's coverage term. That also informs D/E/F sizing. It is a slot decision, not a spec.
2. **The board's weak end is thinner than proportional scaling assumes.**
   - The nmae cohorts 209–261/300 put 13–30 % of targets at a gate < 10, against H1's 32 % at f 1 and 42 % at f 0.35–0.5.
   - The reach cohort (290–300) caps nulls at 3.3 %, against 9 % at f ≤ 0.5.
   - Any lever that lives in small pools (reach, jac on weak targets) is over-priced by f 0.35–0.5 PowerSurfaces. Use f 1
     or the compressed map (`table521.cmp_counts`) beside them.
3. **"Unexpressed target → null" is refuted.**
   - Panel targets below 5 CPM in the line's controls are null 26 % of the time on the K562 fold (11/42), about the same as
     expressed targets (23 %). At 250 cells it is 0/11, and 0 on the Orion folds.
   - Calling nothing on them breaks even only at 39 %.
   - A/B/C have 1 / 5 / 6 such targets (`expr521.log`).

## 9. Files

- `ka521.py`: KA1–2. `check521.py`: KA4–6.
- `lib521.py`: state reconstruction and the member functions.
- `ceiling521.py`: §3.
- `pred521.py`, `thin521.py`, `expr521.py`: §2.
- `table521.py`: the per-target row tables for 3 seeds; `work/table521_s*.npz` is local, ~2 MB each.
- `policy521.py`: §4. `summary521.py`: the seed table. `global521.py`: §5.
- `select521.py`, `diag521.py`, `purity521.py`, `transfer521.py`, `transfer521b.py`: §6.
- `work/*.log` and `work/*.json` hold every number quoted. `work/*.time` holds the peak RSS of each run.
- Commands: `PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/<script>.py`. `table521.py` and
  `policy521.py` take the seed as their argument.

## Correction (2026-09-29, the 533 claims audit)
D10 (at f 0.35 only K 1,600 excludes 0: +0.0048 [+0.0009, +0.0100]; K 1,200 is +0.0032 [−0.0034, +0.0098]; in cmp both exclude 0). Source: `experiments/533-science-lead/RESULTS.md` §2.3.

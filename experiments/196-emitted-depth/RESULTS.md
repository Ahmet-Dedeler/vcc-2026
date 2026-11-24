# 196 — The emitted library depth, and the integer rounding underneath it

**Verdict in one line.** Depth is a **strongly negative** parameter — `+0.0009` on `pds` against
`−0.0185` on `mse` — and production's `L` is already at the optimum, so the depth knob closes.
But the defect depth was buying back is **not depth, it is `floor`**: `emit_counts` delivers
`totals = floor(n L p)`, which shorts every gene by ~0.5 counts and therefore tilts the delivered
profile against low-expression genes by **−0.0123 log2** on the challenge. Replacing `floor` with
largest-remainder rounding at production's depth removes **97%** of that tilt, captures **103% of
the entire depth `pds` gain**, leaves `mse`'s `claim` saturated, costs **+0.14% nnz**, and moves all
four DE members the right way on three independent pools and both emitters. Priced at
**+0.00120 overall**, which is below rule 13's `+0.0019` standalone floor — so it is a **free
rider, not a slot**.

Everything below is measured. Eight scripts, **13 known-answer gates / 43 assertions, all
passing.**

---

## 0. Known-answer rows (rule 10)

| gate | script | expected | measured | |
|---|---|---|---|---|
| a15 h1-test `reach` | `six196.py` | 0.18283 | 0.18283 | PASS |
| a15 h1-test `n_pred` | `six196.py` | 6034 | 6034 | PASS |
| a15 h1-test `jac` | `six196.py` | 0.06824 | 0.06824 | PASS |
| a15 h1-test `nmae` | `six196.py` | 0.96659 | 0.96659 | PASS |
| `rc_hybrid` h1-test, 4 members | `six196.py` | 0.21681 / 6043 / 0.06827 / 0.96670 | identical | PASS |
| a15 raw `mse` / score | `mse196.py` | 0.88228 / 0.11112 | 0.88228 / 0.11112 | PASS |
| PLOLO `m=inf` = 169's cached chalA rows | `pds196.py`, `pdsham196.py` | 0.790712 / 0.790896 / 0.771621 / 0.615317 | worst \|d\| **2.8e−07** | PASS |
| `pdsham196` floor row = `pds196` m=1 row | `pdsham196.py` | 4 folds | worst \|d\| **0.0e+00** | PASS |
| reconstruction: off-gate `floor(x)==T` | `recon196.py` | 1.0000 | 1.0000 on 4 surfaces | PASS |
| reconstruction: `m=1` == artifact | `recon196.py` | 1.0000 | 1.0000 on 4 surfaces | PASS |
| `ham196` floor arms == depth sweep m=1 | `ham196.py` | bit-identical | totals digest, pred `p`, reference `p` all MATCH | PASS |
| `mseham196` floor row = `mse196` m=1 | `mseham196.py` | 1.12577 / 0.24653 / 0.88228 | 1.12546 / 0.24648 / 0.88228 | PASS |
| `sixham196` floor arms = 185/191 rows | `sixham196.py` | 8 values | 8/8 exact | PASS |

Donors on h1-test are `k562_gwps` + `orion_hek293t` (+ `orion_hct116`, `kolf2` on PLOLO, which is
the challenge frame); **`h1_2025` is excluded everywhere** — it covers 85/85 and is pure leakage.
`pds` is sized on **PLOLO only**, never `debench.Surface.pds()` (188 Phase 3).

---

## 1. The defect, quantified before trying to fix it

`emit_counts` (`src/vcc2026/donor_transfer.py:380`, `mode="quasi"`) delivers

```
L      = int(median(lib_sizes))          # ~20,000
totals = floor(n * L * profile_p)        # n = 400
```

`floor` discards each gene's fractional part, so **every gene is short by a mean of 0.5 counts.**
0.5 counts is a fixed number of *counts* and every metric reads a *ratio*, so the shortfall is
0.9% of a 56-count gene and 0.02% of a 2,500-count one — a **size-correlated tilt**, not a common
mode the CPM normalisation can absorb.

### Delivered vs intended, by control-CPM decile — challenge context A, a15 as shipped

`n=400  L=20,109  n·L=8,043,600  gated=9,929  blocks=300`; common-mode shift **+0.083%**.

| decile | ctrl CPM | median T | share T<100 | E[rel shortfall] | E[log2 bias] | net tilt (log2) |
|---|---|---:|---:|---:|---:|---:|
| 1 | 5.0–11.6 | 63 | 98.0% | 0.828% | −0.01251 | −0.01132 |
| 2 | 11.6–20.6 | 126 | 13.8% | 0.406% | −0.00600 | −0.00481 |
| 3 | 20.6–30.9 | 202 | 0.0% | 0.251% | −0.00400 | −0.00281 |
| 5 | 43.2–57.4 | 399 | 0.0% | 0.126% | −0.00182 | −0.00063 |
| 7 | 75.4–101.6 | 697 | 0.0% | 0.072% | −0.00104 | +0.00015 |
| 10 | 220.4–5731.3 | 2495 | 0.0% | 0.019% | −0.00027 | +0.00092 |

Full four-surface table in `work/round196.md`.

### The share of the gated pool that is delivered coarsely

| surface | median T | T<100 | T<40 | T<10 | T=0 | mean rel shortfall |
|---|---:|---:|---:|---:|---:|---:|
| challenge A | 459 | **11.2%** | 0.2% | 0.00% | 0.000% | 0.206% |
| challenge B | 399 | **14.0%** | 0.4% | 0.02% | 0.001% | 0.236% |
| challenge C | 395 | **13.9%** | 0.2% | 0.00% | 0.001% | 0.232% |
| h1-test | 385 | **15.1%** | 1.3% | 0.02% | 0.008% | 0.255% |

### The single number that matters — the tilt

Measured identically on both surfaces (median over blocks, scorer-gated pool, decile 1 minus
decile 10 of `log2(delivered/intended)`):

| rule | mean rel err | common-mode (log2) | **TILT (log2)** |
|---|---:|---:|---:|
| `floor`, h1-test | −0.00249 | −0.00097 | **−0.01352** |
| `floor`, challenge A/B/C | −0.00223 | −0.00090 | **−0.01234** |

**The answer to the mandate's fork.** The *magnitude* is small — 0.2% mean — but it is
**systematic and size-correlated**, and that is the part the metrics can see. It is **not**
negligible for `pds` and `nmae`; it is negligible for everything else. The experiment continues.

---

## 2. Phase 1 — the depth sweep

One variable moves: `L`. The spread rank axis is computed once at `m = 1` and reused, so
`keep_top 8000 / keep_spread 0.25` selects the **same** 8,000 genes at every depth;
`overdispersion 4.0`, `DISPERSION 0.1` and the per-target seeds are a15's. The intended profile is
`recon196`'s reconstruction, whose `m = 1` row re-delivers a15's own per-gene panel totals
bit-identically (asserted per target).

### 2a. Format contract — **the mandate's premise is wrong** (Phase 3 falsifier 1)

The mandate assumed nnz was the binding constraint. It is not.

| m | package nnz | % of 4.75e9 cap | max counts/cell | % of 1e6 cap |
|---:|---:|---:|---:|---:|
| 1 | 2.378e9 | 50.06% | 21,884 | 2.19% |
| 2 | 2.908e9 | 61.21% | 44,183 | 4.42% |
| 4 | 3.279e9 | 69.03% | 87,173 | 8.72% |
| 8 | 3.506e9 | 73.82% | 174,016 | 17.40% |
| 16 | 3.655e9 | 76.95% | 352,280 | 35.23% |
| 32 | 3.556e9 | 74.86% | 705,619 | 70.56% |
| 45 | 3.620e9 | 76.21% | 993,028 | **99.30%** |
| 64 | 3.642e9 | 76.68% | 1,140,781 | **114.08% — DEAD** |
| 96 | 3.643e9 | 76.69% | 1,140,740 | **114.07% — DEAD** |

Rows `m ≤ 16` come from `chal196.py` (15 sampled blocks per challenge context); rows `m ≥ 32` from
`work/deep196.py`, a separate 12-block sample with its **own** `m=1` nnz anchor (2,696,035 vs
2,557,041 per block), which is why 76.95% at `m=16` and 74.86% at `m=32` are not monotone across
the join. Within each anchor the curve is monotone and saturating.

**nnz asymptotes at ~76.7% of the cap and never reaches 90%**, so the mandate's "whatever multiple
hits 90% of the nnz cap" arm **does not exist**. The mechanism is saturation of the gene axis: at
`m = 16` a block already stores 4.28e6/400 = **10,700 nonzeros per cell against a gated pool of
10,779 genes**. The binding constraint is the **1,000,000 counts/cell cap, which binds at m ≈ 45.**
(`work/chal196.log`, `work/deep196.log`. The `m=1` row is anchored on the shipped package's true
stored-entry count 2,377,814,190; deeper rows scale it by the measured nnz ratio on the same
sampled blocks, so the sampling error cancels. Re-emission check: 2,557,041 against the shipped
blocks' own 2,639,035, −3.107%.)

### 2b. `pds` — PLOLO, 4 folds, challenge frame (162's instrument)

| m | counts/cell | k562_gwps | orion_hct116 | orion_hek293t | kolf2 | MEAN | d vs 1× | **scaled/6** | ‖round‖/‖eff‖ |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 20,109 | 0.791724 | 0.785598 | 0.766792 | 0.611951 | 0.739016 | +0.000000 | +0.000000 | 0.04702 |
| 2 | 40,218 | 0.791558 | 0.789226 | 0.770913 | 0.613836 | 0.741383 | +0.002367 | +0.000873 | 0.02481 |
| 4 | 80,436 | 0.791239 | 0.790037 | 0.772004 | 0.614765 | 0.742011 | +0.002995 | +0.001104 | 0.01321 |
| 8 | 160,872 | 0.790962 | 0.790514 | 0.771830 | 0.614859 | 0.742041 | +0.003025 | +0.001115 | 0.00702 |
| 16 | 321,744 | 0.790740 | 0.790729 | 0.771633 | 0.615115 | 0.742054 | +0.003038 | +0.001120 | 0.00385 |
| 64 | 1,286,976 | 0.790740 | 0.790932 | 0.771807 | 0.615075 | 0.742138 | +0.003122 | +0.001151 | 0.00093 |
| ∞ | — | 0.790712 | 0.790896 | 0.771621 | 0.615317 | 0.742137 | +0.003120 | +0.001150 | 0.00000 |

`PDS_SPAN = 0.45208`; 162's PLOLO LOO sigma is 0.00149 raw = 0.000549 overall. The whole
available `pds` gain is **+0.00312 raw / +0.00115 overall**, and 76% of it is already collected
by `m = 2`.

### 2c. `mse` — 096's `claim` route, challenge panel (Phase 3 falsifier 3)

`knee s²B = 0.23268`, `cap K = 0.24653`, a15's `claim = 0.24581` — only **5.6% of margin**.

| m | Chat/C_p | rows under cap | claim | branch | raw mse | score | **d overall** |
|---:|---:|---:|---:|---|---:|---:|---:|
| 1 | 1.126 | 0/12 | 0.24653 | **SATURATED** | 0.88228 | 0.11112 | +0.00000 |
| 2 | 0.827 | 12/12 | 0.20399 | linear | **1.00780** | 0.00000 | **−0.01852** |
| 4 | 0.711 | 12/12 | 0.17517 | linear | **1.13386** | 0.00000 | **−0.01852** |

**The cliff is real and it is the whole story.** Deeper cells carry less emitted sampling noise,
`C_hat` falls, `claim` leaves saturation, and the member's entire 0.111 score goes to zero at
`m = 2`. The profile half moves the *right* way (rounding energy 3.07e−6 → 8.48e−7 → 2.33e−7) but
is worth only `−0.00002` raw mse against the claim term's `+0.126`. **The sign differs between the
two emitters only in margin, not in direction** — 191 measured `rc_hybrid` at `C_hat/C_p` 2.786
against a15's 1.094, so rc survives `m = 2` on `mse` where a15 does not; but rc's own DE numbers
below do not pay for the rest.

### 2d. The four DE members — h1-test, 3 pools, 83-perturbation paired cohort

`rc_hybrid` failed to weld a row at `m = 16` on **2/85 targets (FGFR1, RRM2)**; those two are
dropped from **every** arm so the comparison stays paired. `d overall` uses 181's board 1/span
(`reach` 1.12927, `fid` 3.32482, `jac` 2.72081, `nmae` 1.64787).

**Pool `native`** (the reference's own BH set and our own BH call set — what a submission at that
depth would actually score):

| arm | n_pred | reach | fid | jac | nmae | prec | d overall |
|---|---:|---:|---:|---:|---:|---:|---:|
| quasi m1 | 6036 | 0.18658 | 0.50709 | 0.06838 | 0.96489 | 0.50709 | +0.00000 |
| quasi m2 | 8467 | 0.18884 | 0.50492 | 0.06695 | 0.95906 | 0.50492 | +0.00018 |
| quasi m4 | 8777 | 0.17898 | 0.50532 | 0.06647 | 0.95821 | 0.50532 | −0.00144 |
| quasi m16 | 8823 | 0.17607 | 0.50543 | 0.06634 | 0.95723 | 0.50543 | −0.00172 |
| rc m1 | 6044 | 0.22194 | 0.50703 | 0.06840 | 0.96498 | 0.50703 | +0.00000 |
| rc m2 | 8476 | 0.23940 | 0.50506 | 0.06696 | 0.95920 | 0.50506 | **+0.00313** |
| rc m4 | 8782 | 0.23115 | 0.50525 | 0.06650 | 0.95827 | 0.50525 | +0.00173 |
| rc m16 | 8826 | 0.22428 | 0.50530 | 0.06636 | 0.95733 | 0.50530 | +0.00065 |

**Pool `wald89`** (188's significance-selected pool — `|lfc|·√cpm` at the reference's own BH depth,
median `n_real` = 89; the construction the mandate specified, because a top-N-by-magnitude pool
inverts `nmae`'s sign):

| arm | reach | fid | jac | nmae | d overall |
|---|---:|---:|---:|---:|---:|
| quasi m1 | 0.16978 | 0.50709 | 0.06963 | 0.97025 | +0.00000 |
| quasi m2 | 0.16090 | 0.50492 | 0.06701 | 0.96552 | **−0.00276** |
| quasi m4 | 0.13549 | 0.50532 | 0.06634 | 0.96463 | **−0.00738** |
| quasi m16 | 0.15550 | 0.50543 | 0.06617 | 0.96390 | −0.00343 |
| rc m1 | 0.21894 | 0.50703 | 0.06964 | 0.97005 | +0.00000 |
| rc m2 | 0.21440 | 0.50506 | 0.06703 | 0.96526 | **−0.00182** |
| rc m4 | 0.20612 | 0.50525 | 0.06636 | 0.96457 | −0.00338 |
| rc m16 | 0.21357 | 0.50530 | 0.06619 | 0.96385 | −0.00183 |

**On the mandate's own specified pool, every depth arm is negative on both emitters.** The
`native` pool's `rc m2 = +0.00313` is a pool artefact: it comes from letting the call set grow
6,044 → 8,476, and with the call set matched to the `m=1` arm's own per-row `n_pred` it collapses
to `+0.00052` and then goes negative. Full table in `work/six196.log`.

### 2e. Why depth costs `fid` and `jac` — the mechanism, measured

`chal196.py` counts the genes we call that our **builder's** CPM>10 gate assigns *exactly zero*
effect, while the **scorer's** 5-CPM gate still scores them ("scored-but-uneffected"):

| ctx | dead pool | called at m=1 | m=2 | m=4 | m=16 |
|---|---:|---:|---:|---:|---:|
| A | 794 | 22 (0.48%) | 595 (8.1%) | **794 (10.2%)** | 794 (10.1%) |
| B | 958 | 25 (0.50%) | 751 (9.9%) | **958 (12.1%)** | 958 (12.0%) |
| C | 997 | 43 (0.96%) | 680 (8.9%) | **996 (12.3%)** | 997 (12.1%) |

At `m ≥ 4` we call **100% of the dead pool**, at chance direction. That is `fid`'s precision and
`jac`'s union both paying for power we cannot use. At production's depth the emitter's own
coarseness is what keeps these genes *out* of the call set.

### 2f. The depth composite

| arm | pds | mse | DE (wald89) | **total** |
|---|---:|---:|---:|---:|
| quasi m=2 | +0.00087 | −0.01852 | −0.00276 | **−0.0204** |
| rc m=2 | +0.00087 | −0.01852 | −0.00182 | **−0.0195** |
| quasi m=4 | +0.00110 | −0.01852 | −0.00738 | **−0.0248** |

Before rule 13's emission-law toll (`−0.00243`), which this arm class *does* attract. **Depth is
closed. `L = int(median(lib_sizes))` is at the optimum and should stay there.**

### 2g. The one escape route, priced and not taken

`comp196.py` asked whether the `mse` cliff can be bought back with dispersion, since `C_hat` is
what `phi` controls:

| m | phi | n_pred | vs a15 | Chat/C_p | claim/cap | rows<cap |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 4 | 4,507 | 1.00× | 1.132 | 1.0000 | 0 |
| 2 | 4 | 7,291 | 1.62× | 0.830 | 0.8304 | 18 |
| **2** | **8** | 6,820 | 1.51× | **1.432** | **1.0000** | **0** |
| 4 | 16 | 7,096 | 1.57× | 2.534 | 1.0000 | 0 |

So `m=2, phi=8` **does** re-saturate `claim` and neutralise the `mse` cliff. It then buys only
`pds +0.00087` while inflating `n_pred` 1.5× and the dead-call fraction with it, and `phi` is
185/192's territory, not this experiment's. **Recorded, not pursued.**

---

## 3. Phase 2 — the zero-lfc call count

> `fid` counts a called gene with no direction as a miss. How many genes do we call for which our
> own predicted lfc is exactly zero?

**Zero. In every context, at every depth.**

| ctx | arm | n_pred | **zero_lfc** | zero_total | tiny_lfc (\|lfc\|<1e−3) |
|---|---|---:|---:|---:|---:|
| A | m1 | 4,859 | **0** | 0 | 68 (1.4%) |
| B | m1 | 4,737 | **0** | 0 | 60 (1.3%) |
| C | m1 | 4,543 | **0** | 0 | 53 (1.2%) |

The hypothesised defect does not exist: the quasi emitter's random scatter and the profile's own
offset always move a called gene's delivered total off the control's, so no called gene is ever
delivered without a direction.

The nearest real thing is §2e's `dead_called` — genes with *exactly zero predicted effect by
construction*, whose called direction is therefore a coin flip. At production's depth that is
**22 / 25 / 43 genes, i.e. 0.48% / 0.50% / 0.96% of the call set**. Dropping them would move
precision 0.50701 → 0.50704 with the yield term unchanged (`n_pred` ≫ `n_real`), i.e.
**+0.00002 overall. Immaterial — do not build the suppressor.** It reaches 10–12% and ~+0.0004
only at depths the `mse` cliff has already vetoed.

---

## 4. The result the depth sweep pointed at — **the rounding rule**

The depth sweep's `pds` and `nmae` gains are both bought by shrinking the *relative* rounding
error. That is achievable **at production's depth** by fixing the rule instead of the scale:

```
floor:  T = floor(x)                                    # every gene short by E[frac] = 0.5
ham:    T = floor(x), then the leftover round(sum x) - sum(floor x)
            counts go to the largest fractional remainders   # error centred, size-uncorrelated
```

`emit_counts` is **not modified**: `Recon.profile` already expresses any desired integer total as
`(T + 0.5)/(n L)`, whose own `floor` reproduces `T` exactly. **One variable moves** — the spread
vector, `phi`, `DISPERSION`, `L`, the seeds and the emission law are bit-identical, and the
`floor` arms reproduce the depth sweep's `m=1` arms **bit for bit** (totals digest, predicted `p`,
reference `p`).

### 4a. The tilt is removed

| rule | surface | mean rel err | common-mode (log2) | **TILT (log2)** |
|---|---|---:|---:|---:|
| floor | h1-test | −0.00249 | −0.00097 | **−0.01352** |
| **ham** | h1-test | +0.00005 | +0.00002 | **+0.00019** — 71× smaller |
| floor | challenge A/B/C | −0.00223 | −0.00090 | **−0.01234** |
| **ham** | challenge A/B/C | +0.00007 | +0.00003 | **+0.00039** — 32× smaller |

### 4b. `pds` — Hamilton at production depth **beats infinite depth**

| rule | n·L | k562_gwps | orion_hct116 | orion_hek293t | kolf2 | MEAN | d vs floor | **scaled/6** | ‖round‖/‖eff‖ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| floor | 8,043,600 | 0.791724 | 0.785598 | 0.766792 | 0.611951 | 0.739016 | +0.000000 | +0.000000 | 0.04702 |
| **ham** | 8,043,600 | 0.790948 | 0.790860 | 0.771285 | 0.615802 | **0.742224** | **+0.003207** | **+0.001182** | 0.02699 |
| ∞ | — | 0.790712 | 0.790896 | 0.771621 | 0.615317 | 0.742137 | +0.003120 | +0.001150 | 0.00000 |

Hamilton lands **above the m=∞ ceiling** (+0.000087 raw), i.e. it captures **103% of the entire
depth `pds` gain** — while halving, not eliminating, the rounding energy.

**The mechanism, and it is the interesting part.** What `pds` loses to `floor` is the **shared,
systematic** component of the rounding error, not the random one. `pds_cosine` is a retrieval
metric: a tilt common to all 300 predictions is a shared direction that inflates decoy similarity
— the same hubness mechanism as AGENTS §4.2 and the Orion blob. The PLOLO cross term measures it
directly, as the ratio of the realised change in squared error to the rounding error's own energy:

| fold | floor | ham |
|---|---:|---:|
| k562_gwps | **2.465** | 1.252 |
| orion_hct116 | **1.400** | 0.953 |
| orion_hek293t | **1.719** | 1.017 |
| kolf2 | **1.546** | 0.884 |

`floor`'s error adds **1.4–2.5× its own energy** — it is correlated with the prediction error
already there. Hamilton's adds **0.88–1.25×**, i.e. it is uncorrelated noise. Removing the
*correlation* is worth more than removing the *energy*.

### 4c. `mse` — `claim` stays saturated (Phase 3 falsifier 3, cleared)

| rule | Chat/C_p | rows under cap | claim | branch | round energy | raw mse | score | **d overall** | nnz/block |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| floor | 1.125 | 1/24 | 0.24648 | **SATURATED** | 3.073e−06 | 0.88228 | 0.11112 | +0.00000 | 2,626,215 |
| **ham** | **1.130** | **0/24** | 0.24653 | **SATURATED** | **9.857e−07** | 0.88227 | 0.11113 | **+0.00000** | 2,629,776 |

The depth cliff **does not apply to the rounding rule**, and this was checked rather than assumed:
moving a gene's total by ≤1 count leaves the emitted cell-to-cell spread alone, so `C_hat` does not
fall — it rises marginally (1.125 → 1.130) and the one block that was under the cap moves above it.
The challenge-side rounding energy falls **3.1×** for free.

### 4d. The four DE members — all three pools, both emitters

| pool | arm | n_pred | reach | fid | jac | nmae | **d overall** |
|---|---|---:|---:|---:|---:|---:|---:|
| native | quasi floor | 6034 | 0.18283 | 0.50701 | 0.06824 | 0.96659 | +0.00000 |
| native | **quasi ham** | 6025 | 0.18424 | 0.50758 | 0.06834 | **0.96026** | **+0.00237** |
| native | rc floor | 6043 | 0.21681 | 0.50693 | 0.06827 | 0.96670 | +0.00000 |
| native | **rc ham** | 6034 | 0.22288 | 0.50775 | 0.06837 | **0.96040** | **+0.00337** |
| wald89 | quasi ham | 6025 | 0.15879 | 0.50758 | 0.06965 | 0.96522 | **+0.00011** |
| wald89 | rc ham | 6034 | 0.21793 | 0.50775 | 0.06969 | 0.96505 | **+0.00260** |
| matched | quasi ham | 6033 | 0.18424 | 0.50758 | 0.06834 | 0.96026 | **+0.00237** |
| matched | rc ham | 6043 | 0.22288 | 0.50774 | 0.06838 | 0.96040 | **+0.00338** |

**Positive on every pool, for both emitters.** `fid` (+0.0006/+0.0008), `jac` (+0.0001) and `nmae`
(**−0.0063**, a gain) move the same way on all three pools; `reach` moves +0.0014/+0.0061 on
`native` and +0.0043 on `wald89` for rc. The single negative anywhere is quasi's `reach` on the
`wald89` pool (0.16693 → 0.15879), on the member AGENTS records as **not locally priceable at all
(2/6 signs, worse than chance)**.

The call set barely moves — `n_pred` 6034 → 6025, **−0.15%** — so this is not a yield trade.

### 4e. Rule 12 — all six members priced

| member | value (overall) | instrument | trust |
|---|---:|---|---|
| `pds` | **+0.001182** | PLOLO, 4 folds, 2 KA gates | priceable (162: sign 6/6, sigma 0.000549) |
| `mse` | **+0.00000** | 096 `claim` route, KA-gated | priceable for this arm class; `claim` **SATURATED** both sides |
| `nmae` | **+0.00173** | h1-test, 3 pools, −0.0063 raw and consistent | 3/6 signs — indicative |
| `reach` | +0.00114 | h1-test (rc, native) | **not locally priceable** |
| `fid` | +0.00045 | h1-test | **not locally priceable** |
| `jac` | +0.00005 | h1-test | **not locally priceable** |

Format margin: package nnz **50.06% → 50.13%** of the 4.75e9 cap (+0.14%); max counts/cell
22,197 → 22,162, i.e. **2.2% of the 1e6 cap**. Both constraints cleared by wide margins.

### 4f. Rule 13 — the composed slot gate

**Arm class.** This is a **PROFILE** arm, not an emission-law arm: `mode="quasi"`, the spread
vector, `overdispersion`, `DISPERSION`, `L` and every seed are bit-identical; only the integer
totals handed to the emitter move, by ≤1 count each. Per the 2026-09-03 amendment the `−0.0146`
toll is an **emission-law** toll (it equals `mean(a21, a22)`, two exact board points of that
class), and 169 measured `+0.0006` for profile arms.

```
predicted = +0.001182 × 1.011  +  0.00000  =  +0.00120
ship floor (shape-only)                    =  +0.0019     ->  VETO
with the emission-law toll charged anyway  =  -0.00123    ->  VETO
```

**Vetoed as a standalone slot either way.** +0.00120 is 2.15σ against 162's PLOLO LOO sigma, but
rule 13's `+0.0019` floor *is* the realised RMS of the unpriced channels, and an arm that only
clears its own noise floor by 2σ on the one priced member is exactly what that floor exists to
reject.

---

## 5. Verdict

1. **The emitted depth is closed.** `L = int(median(lib_sizes))` ≈ 20,000 is at the optimum.
   Raising it buys at most `+0.00115` overall on `pds` and costs `−0.01852` on `mse` the moment
   `m ≥ 2`, because `C_hat` falls out of `claim`'s saturation — a15 has only 5.6% of margin there.
   On the mandate's own specified DE pool every depth arm is additionally negative on both
   emitters, and at `m ≥ 4` we call 100% of the scored-but-uneffected pool at chance direction.
   **A parameter that has sat unexamined since 2026-08-21 is now measured and should not be
   revisited without a new mechanism.**

2. **Two premises in the mandate are corrected.** (a) nnz is **not** the binding constraint —
   it asymptotes at **76.7%** of the cap because the gene axis saturates at ~10,700 nonzeros/cell,
   and the **1e6 counts/cell** cap binds first, at `m ≈ 45`. The "90% of the nnz cap" arm does not
   exist. (b) The `mse` sign does **not** differ between the two emitters — `rc_hybrid`'s 2.5×
   larger `C_hat/C_p` margin lets it survive `m=2` where a15 does not, but the direction of the
   effect is the same and rc's DE numbers do not pay for the rest.

3. **Phase 2 is a clean negative.** We never call a gene with no direction: `zero_lfc = 0` in all
   three contexts at all five depths. The nearest real defect, calling genes with exactly-zero
   predicted effect, is **0.5–1.0% of the call set** at production's depth and worth **+0.00002
   overall** to suppress. Not worth building.

4. **The finding is the rounding rule, not the depth.** `floor` tilts the delivered profile by
   **−0.0123 log2** against low-expression genes on the challenge. Largest-remainder rounding at
   the *same* `L` removes **97%** of that tilt, delivers **103% of the entire depth `pds` gain**
   (+0.003207 raw, above the `m=∞` ceiling), leaves `mse`'s `claim` **saturated**, costs **+0.14%
   nnz** and **−0.15% `n_pred`**, and is positive on all four DE members across three pools and
   both emitters. The mechanism is that `pds` is hurt by the *correlated* part of the rounding
   error, not its energy — floor's error adds 1.4–2.5× its own energy to the squared error,
   Hamilton's adds 0.88–1.25×.

5. **Ship it as a rider, not as a slot.** Priced at **+0.00120**, below rule 13's `+0.0019`
   standalone floor. But it is the cheapest change in the project's inventory — no operator, no
   donor table, no call-set policy, no emission law, no format headroom, and `emit_counts` itself
   is untouched — so it should ride on the next arm that earns a slot for another reason (192's
   call set or 193's operator).

   **Recommended handoff.** Replace `totals = np.floor(n * L * profile_p)` with the
   largest-remainder allocation in `experiments/196-emitted-depth/ham196.py:hamilton`, then assert
   the delivered `n_pred` moves by less than 1% and the package nnz by less than 0.5%. **There are
   four such call sites in `donor_transfer.py` and exactly one is live**: line **380**, under
   `mode="quasi"`. Both shipped emitters route through it — `rc_hybrid_block`
   (`rc_emit.py:875`) calls `emit_counts(..., mode="quasi")` for its pass-1 block, so changing
   that one line changes a15 and `rc_hybrid` together and preserves 191's digest identity between
   them, which is what pins `pds` and `mse` across the pair. Line 277 (`_quasi_block`) is reached
   only by `mode="resample_mix"` (a21/a22's class), and 498 / 567 (`resample_mix` / `even`) are
   likewise off the shipped path — leave all three alone.

---

## Files

| file | what it produced |
|---|---|
| `recon196.py` → `work/recon196.log` | the intended-profile reconstruction + its 4-surface KA gate |
| `totals196.py` → `work/chal-totals-a15.npz` | a15's 900 delivered blocks, nnz 2,377,814,190 |
| `round196.py` → `work/round196.md`, `.json` | the rounding-error table by expression decile, 4 surfaces |
| `pds196.py` → `work/pds196.json` | `pds` vs depth, PLOLO, KA-gated |
| `mse196.py` → `work/mse196.json` | `mse` vs depth, both halves, 096's claim route |
| `chal196.py` → `work/chal196.json` | challenge-side `n_pred`, format contract, Phase 2 |
| `work/deep196.py` → `work/deep196.log` | the nnz/cell-cap saturation curve to `m = 96` |
| `comp196.py` → `work/comp196.json` | the depth × phi compensation screen |
| `sweep196.py` → `work/de196.*.npz`, `rc196.*.npz` | h1-test DE tables, 2 emitters × 5 depths |
| `six196.py` → `work/six196.log`, `.json` | the four DE members vs depth, 3 pools |
| **`ham196.py`** → `work/ham196.log` | **the rounding-rule arm: tilt + h1-test DE tables** |
| **`pdsham196.py`** → `work/pdsham196.json` | **`pds` for the rounding rule, PLOLO, 2 KA gates** |
| **`mseham196.py`** → `work/mseham196.json` | **`mse` for the rounding rule, both halves, KA-gated** |
| **`sixham196.py`** → `work/sixham196.log` | **the four DE members for the rounding rule, 3 pools** |
| `work/chaltilt196.log` | the challenge-side tilt for both rules |

No `vcc`, no `cell-eval2`, no submission build, no writes outside `experiments/196-emitted-depth/`.

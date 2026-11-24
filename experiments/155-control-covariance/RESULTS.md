# 155 — Control-cell covariance as a per-gene direction source

**Verdict: null.** The context's own control cells carry no usable per-gene knockdown
direction. The estimator as specified — every gene regressed on the target gene across
18,400 non-targeting controls, residualised on library size and the top 20 PCs — has a
paired transfer cosine of **−0.00091 ± 0.00089 (t(131) = −1.03)** against the incumbent's
**+0.01782 ± 0.00208 (t(131) = +8.59)**. Its split-half reproducibility on the controls
alone is **r = 0.028**. There is nothing there to transfer, and the failure is upstream of
the reference: the coefficients do not reproduce against *themselves*.

Everything that looks positive before that control is applied is the shared cell-state
programme (k = 0) or the reference pool's 59.7%-down marginal, and both are reproduced
exactly by dealing the coefficient rows to the **wrong** targets.

---

## 0. Frame, and reproduction of the target numbers

One fixed masking for every arm (`bench.py`): 153's `Surface('../095-clean-harness/
de_thinned.npz')`, intersected with 104's frozen HCT116 axis and the whitener's 17,942
canonical columns — **132 perturbations × 9,980 genes**, reference DE pool 10,998 gene-calls.
This is exactly 153 q9's gamma = 0.65 row.

| | sign acc | head purity | reach | pds |
| --- | ---: | ---: | ---: | ---: |
| target (brief) | 0.5366 | 0.6646 | 0.09040 | 0.7953 |
| **reproduced here** | **0.5366** | **0.6646** | **0.09040** | **0.7899** |

Sign accuracy, head purity and reach reproduce to the last digit. **`pds` does not: I get
0.7899, not 0.7953.** Checked four ways — full 10,261-gene gated axis (0.7899),
`exclude_panel=False` (0.7898), and 104's other two arms (`agreement` 0.7872, `full`
0.7993). None is 0.7953, so the brief's `pds` came from a configuration I cannot
reconstruct. Deltas below are against my own 0.7899.

`nmae` for the incumbent on this frame is 0.9690.

## 1. The estimator

`fit.py`. Controls only — the 18,400 `non-targeting` cells of
`orion-hct116-thinned.h5ad`; the 132 perturbation cohorts are opened only by the scorer,
through 095's cached table. 095's `norm_log` verbatim (CPM → 1e6, `log1p`). Covariates
`Z` = intercept, log library size, top *k* PCs. One residualisation and one cross-product:

```
R = X − Q(QᵀX),  Q = qr(Z);      beta[t,g] = ⟨R_t, R_g⟩ / ⟨R_t, R_t⟩
```

Candidate effect = **−beta** (knockdown lowers T). The convention is verified against the
data, not assumed: the reference's own-gene log2FC is negative for **132/132** targets,
median −1.855.

Nine of the 132 targets sit below the 5-CPM gate and are not columns of the scored axis;
they are carried as an appended block so every target still has its own regressor. Target
control CPM: median 39.8, quartile boundaries 17.9 / 39.2 / 78.5, min 0.04 (GLI3), max
655.4 (ANP32B). Median target zero-fraction 0.611.

## 2. Headline table — candidate vs incumbent

`final.py`. Both sign conventions scored explicitly. Candidates L2-matched to the
incumbent's per-row norm, so magnitude is never a hidden knob.

| arm | sign acc | head purity | reach | pds | nmae |
| --- | ---: | ---: | ---: | ---: | ---: |
| **INCUMBENT** 104 baseline, γ=0.65 | **0.5366** | **0.6646** | **0.09040** | **0.7899** | **0.9690** |
| null: all-down constant | 0.5974 | 0.5255 | 0.05872 | 0.5000 | 1.3723 |
| k=0 −beta *(the hypothesis)* | 0.4323 | 0.3907 | 0.11293 | 0.4669 | 1.0167 |
| k=0 +beta *(sign flipped)* | 0.5677 | 0.6093 | 0.17469 | 0.5331 | 0.9833 |
| k=5 −beta | 0.4835 | 0.4961 | 0.04263 | 0.5184 | 1.0004 |
| k=20 −beta | 0.4959 | 0.4892 | 0.04119 | 0.5245 | 1.0002 |
| k=20 +beta | 0.5041 | 0.5108 | 0.02682 | 0.4755 | 0.9998 |
| k=50 −beta | 0.4955 | 0.4788 | 0.03144 | 0.5187 | 1.0007 |

**The hypothesis as stated is refuted on its own endpoint.** With the specified 20 PCs
removed, −beta scores sign accuracy 0.4959 and head purity 0.4892 — chance, on both, to
three decimals. Nothing beats the incumbent on either.

`nmae` ≈ 1.000 for every k ≥ 5 arm is not informative: β magnitudes are ~10⁻² against
reference log2FCs of ~1, so the prediction reads as zero and `nmae` returns its
predict-nothing value.

Two rows need explaining before they are believed, and §4 and §5 do that:

* **the all-down constant beats the incumbent on sign accuracy** (0.5974 vs 0.5366);
* **k=0 +beta posts reach 0.1747, nearly twice the incumbent**, while being *worse* on
  head purity.

## 3. Blends

`final.py`, `(1−f)·incumbent + f·candidate`, candidate L2-matched per row.

| candidate | f | sign acc | head purity | reach | pds | Δacc | Δhead | Δreach | Δpds |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| k=0 −beta *(hypothesis)* | 0.02 | 0.5351 | 0.6655 | 0.08826 | 0.7831 | −0.0015 | +0.0009 | −0.00214 | −0.0068 |
| | 0.05 | 0.5341 | 0.6603 | 0.08809 | 0.7679 | −0.0025 | −0.0043 | −0.00231 | −0.0220 |
| | 0.10 | 0.5292 | 0.6525 | 0.08482 | 0.7394 | −0.0075 | −0.0121 | −0.00558 | −0.0505 |
| | 0.20 | 0.5195 | 0.6318 | 0.08174 | 0.6799 | −0.0172 | −0.0328 | −0.00866 | −0.1100 |
| k=20 −beta *(hypothesis)* | 0.02 | 0.5366 | 0.6672 | 0.09035 | 0.7890 | +0.0000 | +0.0026 | −0.00005 | −0.0009 |
| | 0.05 | 0.5361 | 0.6664 | 0.09125 | 0.7895 | −0.0005 | +0.0017 | +0.00085 | −0.0004 |
| | 0.10 | 0.5347 | 0.6664 | 0.09285 | 0.7890 | −0.0019 | +0.0017 | +0.00244 | −0.0009 |
| | 0.20 | 0.5363 | 0.6724 | 0.09255 | 0.7856 | −0.0004 | +0.0078 | +0.00214 | −0.0043 |
| k=0 +beta *(sign flipped)* | 0.02 | 0.5377 | 0.6664 | 0.09332 | 0.7927 | +0.0011 | +0.0017 | +0.00292 | +0.0028 |
| | 0.05 | 0.5392 | 0.6759 | 0.09522 | 0.7954 | +0.0025 | +0.0112 | +0.00481 | +0.0055 |
| | 0.10 | 0.5439 | 0.6776 | 0.09924 | 0.7814 | +0.0073 | +0.0130 | +0.00884 | −0.0085 |
| | 0.20 | 0.5534 | 0.6880 | 0.10122 | 0.7309 | +0.0167 | +0.0233 | +0.01082 | −0.0590 |

The hypothesis's own blends are flat-to-negative everywhere. The **sign-flipped** k=0 blend
is positive on all four endpoints at once at f = 0.05 — which is the shape this repo has
been burned by before, so §5 permutes it.

## 4. Why the raw numbers move at all — two nulls, both fatal

### 4.1 Nothing here beats the marginal

`strat.py`. The reference DE pool is **59.7% down-regulated**. A single constant all-down
vector therefore scores sign accuracy 0.5974 for free, above the incumbent's 0.5366 — the
incumbent itself predicts down for only 50.1% of pool genes.

The candidate predicts down for 44.6% (k=0) and 49.0% (k=20) of pool entries. Its aggregate
sign accuracy is a marginal-matching statistic and nothing more: no arm, in any stratum,
exceeds the stratum's own marginal (0.5833 for the best-measured 33 targets, 0.5974 overall).

The additive down-shift is nonetheless the one lever in this experiment that genuinely
flips signs, so it was priced (`blend.py` §3): `incumbent − c·(row sd)` gives, at
c = 0.10 / 0.50 / 2.00, Δreach +0.0047 / +0.0216 / +0.0350 against Δpds −0.0103 / −0.1312 /
−0.2431. Using 153's anchors that is **+0.00089 / +0.0041 / +0.0066 overall from `reach`
against −0.0038 / −0.0484 / −0.0896 from `pds`** — a losing trade at every point, by 4× at
its mildest. Same shape as the whitening trade in 153 §8, and closed for the same reason.

### 4.2 Target identity is not carrying any of it

Within-stratum row permutation, 200 seeds (`strat.py`) — rows dealt only among the
stratum's own targets, so every permuted row carries a comparable global offset:

| stratum | k | cand sign acc | permuted | z | cand head | permuted | z |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| top quartile CPM (33) | 0 | 0.5487 | 0.5883 ± 0.0336 | **−1.18** | 0.5451 | 0.5777 ± 0.0483 | −0.67 |
| top quartile CPM (33) | 20 | 0.5024 | 0.4929 ± 0.0101 | +0.94 | 0.4965 | 0.5243 ± 0.0278 | −1.00 |
| top half CPM (66) | 0 | 0.4817 | 0.5104 ± 0.0248 | −1.16 | 0.4707 | 0.5169 ± 0.0366 | −1.26 |
| top half CPM (66) | 20 | 0.4901 | 0.4876 ± 0.0076 | +0.33 | 0.5009 | 0.5112 ± 0.0219 | −0.47 |
| all 132 | 0 | 0.4323 | 0.4353 ± 0.0249 | −0.12 | 0.3907 | 0.4236 ± 0.0318 | −1.03 |
| all 132 | 20 | 0.4959 | 0.4920 ± 0.0052 | +0.76 | 0.4892 | 0.5078 ± 0.0145 | −1.28 |

Every z is in [−1.28, +0.94]. **Correct target assignment is worth nothing on any
endpoint, in any stratum, at any k.**

This also retracts an intermediate result of my own. A first pass permuted rows across all
132 targets while scoring only the top-33 stratum and reported z = **+2.88** for the k=0
candidate. That null is deflated: it deals low-CPM rows, which carry a much larger global
offset, into a high-CPM stratum, so the permuted arm is penalised for a reason unrelated to
target identity. Matched within-stratum, the same cell is **z = −1.18**. The +2.88 was an
artefact of the null, not a signal.

## 5. The sign-flipped blend, permuted — it is a gene-level effect

`plus.py`, `cosine.py`. The k=0 **+beta** arm, standalone and blended, against 50 and then
500 row-permutations of itself:

| | endpoint | real | permuted | z |
| --- | --- | ---: | ---: | ---: |
| standalone | sign acc | 0.5677 | 0.5735 ± 0.0254 | −0.23 |
| | head purity | 0.6093 | 0.5841 ± 0.0337 | +0.75 |
| | reach | 0.1747 | **0.1977** ± 0.0185 | −1.25 |
| | pds | 0.5331 | 0.5067 ± 0.0235 | +1.13 |
| blend f=0.05 (500 seeds) | sign acc | 0.5392 | 0.5404 ± 0.0020 | −0.61 |
| | head purity | 0.6759 | 0.6735 ± 0.0060 | +0.39 |
| | reach | 0.0952 | 0.0926 ± 0.0014 | +1.89 |
| | pds | 0.7954 | 0.7795 ± 0.0074 | +2.16 |
| blend f=0.10 (500 seeds) | reach | 0.0992 | 0.0941 ± 0.0023 | +2.21 |
| | pds | 0.7814 | 0.7590 ± 0.0121 | +1.85 |

Standalone, **the permuted arm outscores the real one on `reach`** (0.1977 vs 0.1747). The
entire two-fold reach advantage over the incumbent survives destroying target identity, so
it is a gene-level property of the k=0 coefficient matrix, not a per-target direction.

The blends leave wisps at z ≈ 1.9–2.2 on `reach` and `pds` — across 5 endpoints × 3
fractions of looking, and §6 shows the better-powered test puts them at zero.

## 6. The decisive measurement — paired transfer cosine

`cosine.py`. `pds` is a rank summary over 132 candidates and discards most of the
information. The underlying quantity is AGENTS.md §4.2's transfer cosine: does target *p*'s
coefficient row point at target *p*'s true effect? Paired per target — own cosine minus the
mean cosine of the other 131 predictions against that same truth, which removes every
per-target scale and magnitude effect — over the scored genes with all panel targets removed.

| arm | own | off-diag | paired | se | t(131) |
| --- | ---: | ---: | ---: | ---: | ---: |
| **INCUMBENT** 104 baseline | +0.01890 | +0.00108 | **+0.01782** | 0.00208 | **+8.59** |
| +beta k=0 | +0.03309 | +0.02626 | +0.00682 | 0.01200 | +0.57 |
| −beta k=0 *(the hypothesis)* | −0.03309 | −0.02626 | −0.00682 | 0.01200 | −0.57 |
| +beta k=20 | −0.00129 | −0.00038 | −0.00091 | 0.00089 | −1.03 |
| **−beta k=20 *(the hypothesis)*** | +0.00129 | +0.00038 | **+0.00091** | 0.00089 | **+1.03** |

Restricting to the best-measured targets, where the estimator has the most power, does not
rescue it — the standard error grows as fast as the estimate:

| arm | stratum | n | paired cos | se | t |
| --- | --- | ---: | ---: | ---: | ---: |
| +beta k=0 | top quartile CPM | 33 | +0.01768 | 0.03140 | +0.56 |
| +beta k=0 | top half CPM | 66 | +0.01335 | 0.02079 | +0.64 |
| +beta k=20 | top quartile CPM | 33 | −0.00162 | 0.00177 | −0.92 |
| +beta k=20 | top half CPM | 66 | −0.00223 | 0.00137 | −1.62 |

**The specified estimator (k=20) is tightly bounded at zero**: 95% CI [−0.0026, +0.0009],
which excludes the incumbent's +0.0178 by twenty standard errors. The k=0 variant is not
bounded — its 95% CI on the top quartile is [−0.044, +0.079] and contains both zero and the
incumbent's value — but that is a statement about noise, not evidence: **the point estimate
is the same size as the incumbent's while the standard error is 15× larger.** Per target,
this estimator is fifteen times noisier than what we already ship, so 132 targets cannot
resolve it even in principle.

## 7. Why — and it is not the reference's fault

`power.py`, controls only. Two diagnostics that never look at a perturbation:

**Coefficient-variance inflation.** Under the null, r ~ N(0, 1/√(N−k−2)) with N = 18,400.

| k PCs | sd(r) observed | sd(r) null | inflation | real variance fraction |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 0.02115 | 0.00737 | 8.23 | 0.879 |
| 5 | 0.00873 | 0.00737 | 1.40 | 0.287 |
| 10 | 0.00797 | 0.00737 | 1.17 | 0.143 |
| 20 | 0.00766 | 0.00738 | **1.08** | **0.072** |
| 50 | 0.00755 | 0.00738 | 1.05 | 0.045 |

**Split-half reproducibility** — fit on two disjoint halves of the control cells, correlate:

| k PCs | r(half1, half2) | after double-centering |
| ---: | ---: | ---: |
| 0 | 0.7792 | 0.7791 |
| 5 | 0.1583 | 0.1597 |
| 20 | **0.0276** | 0.0332 |
| 50 | 0.0130 | 0.0188 |

This is the mechanism, and it is a clean dichotomy:

* **k = 0**: 88% of coefficient variance is real and split-half r is 0.78 — a large,
  highly reproducible signal. It is shared cell-state programme membership (removing five
  PCs collapses it to 0.16), it survives double-centering, and §4.2/§5 show it carries zero
  target-specific direction. *Reproducible, and the wrong thing.*
* **k = 20**, the estimator the hypothesis specifies: inflation 1.08, split-half r = 0.028.
  *There is essentially nothing left to transfer.*

The power budget explains it. Median target expression is 39.8 CPM at a 12,176-UMI median
library — **~0.5 counts per cell**, with a median zero-fraction of 0.611. The regressor is
almost pure Poisson noise, and the resulting attenuation applies to the target and to every
downstream gene at once. 18,400 cells buy a correlation standard error of 0.0074; after
state removal the entire observed correlation distribution sits at 1.08× that.

## 8. One instrument caveat for 153

Worth flagging to whoever owns `debench.py`: on this surface **raw `reach` is not monotone
in per-gene direction quality.** The k=0 +beta arm posts reach 0.1747 at sign accuracy
0.5677 and head purity 0.6093; the incumbent posts 0.0904 at 0.5366 and 0.6646 — better
head purity, half the reach. Decomposing R_p = k*/n_real (`work/reachdiag.py`):

| arm | mean R | median k* | mean k* | frac k*=0 | R on small n_real | on large |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| incumbent | 0.0904 | 3.0 | 5.02 | 0.167 | 0.1147 | 0.0661 |
| +beta k=0 | 0.1747 | 1.0 | 14.77 | 0.356 | 0.1341 | 0.2153 |

The +beta arm is *worse* on the median perturbation (k* 1 vs 3, and 36% of perturbations
score zero) and wins on a heavy tail of long ≥90%-pure prefixes, concentrated on the
**large**-n_real perturbations. So `reach` is being paid by a target-independent consensus
direction that happens to be very accurate on a minority of perturbations. That is a real
mechanism for buying `reach` — and it is the mechanism the additive down-shift in §4.1
prices, at a `pds` cost 4× the gain.

## 9. What this does and does not close

**Closed.** Linear covariance of the target gene against the transcriptome across a single
context's control cells, at challenge-scale depth (18,400 cells, ~12k UMI), does not
estimate the knockdown direction. Not at any PC count, not with row/column/double centering,
not restricted to the best-measured targets, not blended at any fraction. The estimator does
not even reproduce against itself once cell state is removed.

**Not closed, and the specific things that would change the answer.** The failure is a power
failure, not a logical one, so it is worth naming the escapes rather than implying the idea
is dead:

1. **Depth.** The challenge contexts run ~20,000 median UMI against this harness's 12,176 —
   1.6× more counts, and the attenuation enters twice (regressor and response). Still far
   short of what §7 needs, but the number here is a lower bound, not the value.
2. **A denoised regressor.** T's measured expression is ~0.5 counts/cell. kNN or metacell
   smoothing of the *regressor only* is the standard fix and would raise its reliability by
   an order of magnitude. It was not tried here because smoothing on a PC-space neighbour
   graph and then residualising on the same PCs is circular, and resolving that honestly is
   its own experiment. This is the one refinement I would actually run next.
3. **More targets.** The §6 CI on the k=0 variant is [−0.044, +0.079] on 33 targets. The
   full 300-target panel would tighten it ~3×; that still would not reach the incumbent's
   resolution, but it would separate "zero" from "quarter of the incumbent".

Nothing here authorises a build, and nothing is promoted. `cell-eval2` was not run and no
submission was touched.

---

**Files.** `bench.py` (frame + endpoints), `fit.py` (estimator, controls only), `score.py`,
`probe.py` (centering sweep, stratification), `blend.py` (blends, down-shift pricing),
`power.py` (inflation, split-half), `strat.py` (matched permutation), `plus.py`
(sign-flipped permutation), `cosine.py` (paired transfer cosine), `final.py` (headline
table). Logs in `work/*.log`; coefficient matrices `work/coef_k{0,5,10,20,50}.npz`.
`fit.py` caches the 18,400 × 10,270 control matrix in the session scratchpad and rebuilds it
from the harness in ~6 s if it is gone.

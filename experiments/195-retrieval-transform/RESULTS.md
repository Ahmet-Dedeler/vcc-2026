# 195 — Learn the retrieval transform instead of choosing it

**Date:** 2026-09-03 · 3-fold LODO (K562-gwps / Orion-HCT116 / Orion-HEK293T), 259 VCC panel
targets per fold, 7,760-gene feat axis · **70,496 labelled cross-donor retrieval problems**

## Verdict

**NO-GO, and the hand-chosen family is at the optimum.**

`gamma = 0.65` was picked by sweeping the *panel* objective. Run the same family through the
**donor** objective — 70,496 labelled retrieval problems on a held-out perturbation split that
never sees the scored line — and it independently picks **the same gamma**. There is nothing
for a fit to find inside the family.

Outside the family, a low-rank-plus-diagonal correction fitted directly to a listwise ranking
loss **does learn something real and it does not cross a cell line**. Across 23 configurations
it gains up to **+0.0204** macro pds on held-out donor perturbations and loses on the held line
in **21 of 23**; the pre-committed arm (selected on donor dev, never on the panel) scores
**0.79776 → 0.79222, delta −0.00554, 95% CI [−0.01049, −0.00077], P(positive) = 0.0118**.

The mechanism is measured, not inferred: with the scored line untouched, an inner
leave-one-training-**LINE**-out puts the gain at **+0.00904 across unseen perturbations** and
**−0.00618 across an unseen cell line**. The correction dies at the line boundary, and it dies
by the same amount whether the unseen line is a donor line or the scored one. **The panel is
not special. The cell line is the boundary.**

This corroborates experiment 080's terminal NO-GO from an independent harness, with 117x more
training data, and adds the mechanism 080 could not isolate.

---

## 1. What 080 and 084 already settled (read first, as mandated)

**They are this idea, and 080 is terminal.** `experiments/080-listwise-metric-adapter` fitted

```
z' = normalize(z + beta * B U V^T B^T z)
```

— a rank-4/8/16 residual on the frozen post-whitening geometry, `B` a source-only randomized
right-singular basis — with a **listwise loss whose positive is the same perturbation in a
different donor and whose negatives are every other perturbation in that donor**. That is
exactly "fit the transform to the retrieval objective on labelled cross-donor problems", and
algebraically it is the same family this experiment's mandate specifies: the total operator is
`W_prod (I + low-rank)`, with `beta = 0` recovering the incumbent.

080's numbers:

| | value |
|---|---:|
| selected config | `listwise:r=8:beta=0.02:temp=0.03:reg=0.0001` |
| untouched **non-panel** macro PDS gain | **+0.023376** |
| complete-selection permutation-null p95 | +0.021217 |
| held panel, Orion HCT116 | **−0.0006884** |
| held panel, Orion HEK293T | +0.0006435 |
| KOLF2.1J confirmation | +0.0006868 |
| **primary macro (HCT116 + HEK293T)** | **−0.00002245** |
| cluster-bootstrap 95% CI | [−0.001863, +0.001803] |
| P(positive) | 0.4975 |
| `promotion_gate.passed` | **false** |

Result SHA-256 `9b12cdcb440081bd3910c72d00ba3e41f1a137b2024119c1c1a258c6ddb8ccb2`;
selection freeze `8a7335cf21ff6c576c5ce47eaf3a85d14d6e212a62e3680c8cddcf57204f8524`.
RESULTS.md is explicit: *"Do not integrate, retune, or run a production or count-level
confirmation for this adapter."*

Two things in 080's own numbers are worth reading again. Its in-domain gain **+0.023376**
cleared its permutation null p95 **+0.021217** by only **+0.0022** — the real fit beat a
deranged-label control by 10%. And its panel transfer was **−0.0000225**, i.e. exactly nothing.

**084 settled nothing scientific.** `experiments/084-listwise-runtime-study/RESULTS.md` reads
in full: *"Status: implementation and bounded benchmark pending."* It is a read-only FLOP/RSS
accounting of 080's optimizer and carries no measurement of the idea.

**What 080 left open, and why this experiment ran anyway.** 080 fitted on **600** non-panel
targets and could only report that an in-domain gain failed to transfer. It never asked (a)
whether the *family* is at the optimum of the direct objective — the cheapest question in the
whole area, and the one that decides whether "learn it" can help at all; (b) whether the
failure is data starvation, capacity, a query-distribution mismatch, or the cell-line boundary.
Sections 2–7 answer both on **70,496** training queries, and the answers close the area.

---

## 2. Phase 0 — known-answer rows (all pass)

The harness is `experiments/019-whiten-interact/interact.py` rebuilt on a compact cache
(`build_cache.py` → `harness.py`); 019 keeps two float32 copies of each donor's full
38,606-gene deviation matrix, which is ~16 GB for the two Orion lines and over this run's
ceiling. Nothing it computes leaves the 7,760-gene `feat` axis, so the rebuild streams each
`delta.npy` from a memmap and keeps only those columns. `_norm_scale`'s 25th-percentile floor
is still taken over the **full** axis, and the `default_rng(0)` draw order for the whitener
basis is replicated fold-by-fold, because both change the operator.

| row | AGENTS.md §4.2 / exp 019 | measured | delta |
|---|---:|---:|---:|
| 4:1 donors, `gamma = 0` (incumbent) | 0.7378 | 0.73779 | **−0.00001** |
| 4:1 donors, whitened | 0.7873 | 0.78726 | **−0.00004** |
| equal donors, `gamma = 0.55` | 0.7963 | 0.79629 | **−0.00002** |
| **equal donors, `gamma = 0.65`** | **0.7978** | **0.79776** | **−0.00004** |

Residuals are the float64-vs-float32 accumulation difference in the perturbation mean and the
`gene_norm` divisor. Macro cosine at `gamma = 0` reproduces §4.2's **0.0435** exactly.
`work/phase0.json`, `work/phase0.log`. Optimiser known-answer rows (reduced-coordinate cosine
vs full space `3.6e-16`; analytic vs finite-difference gradient `8.6e-9` relative; `M = 0`
bit-identical to the un-adapted operator; in-batch pds equal to the direct rank computation to
`1e-12`) in `test_fitlib.py`.

### Rule 11 — operator paths and sha256

The three operators this experiment fitted and used, `sha256` over the `.npz` bytes:

| path | max eigenvalue | sha256 |
|---|---:|---|
| `experiments/195-retrieval-transform/work/operators/op.k562_gwps.npz` | 1746.5413 | `b3ecd57dcb2d5ef2…` |
| `experiments/195-retrieval-transform/work/operators/op.orion_hct116.npz` | 1236.3177 | `13589045d9264db5…` |
| `experiments/195-retrieval-transform/work/operators/op.orion_hek293t.npz` | 645.7621 | `3149945b4ceab109…` |

Their sha256 over the raw `(evals, evecs)` arrays — the form quoted in `work/phase0.json` and
re-checked by `prep.py` before every use — are `70a86069208c7fb8…`, `5086ff76b9a62e58…`,
`94c2439fb0ed9f8b…`.

**These are NOT any of the operators on disk** and must never be swapped for one. They are
per-fold LODO operators fitted in-process on the `gene_norm`'d, scale-matched, 7,760-gene feat
basis, which is why their eigenvalues are three orders of magnitude larger than production's.
The nine pre-existing whiteners, none of which this experiment loaded, are tabulated with
paths, sha256 and max eigenvalues in `work/wrapup.json` → `operator_table` (`whitener_full`
4.4586, `data/donors/whitener_union` 6.2567, `whitener_union_full` 6.9802,
`whitener_union_pruned` 5.8946, `040/out/whitener_union` 6.8999, `095/out/whitener_union`
8.0364, `099` hek293t 3.2940, `140` k562 3.4289, `140` non-k562 2.7014).

Donor inputs are hashed too (`work/cache/manifest.json`): `k562_gwps/delta.npy`
`0c9422280fed7afb…`, `orion_hct116/delta.npy` `72c7a01dc2869bd3…`, `orion_hek293t/delta.npy`
`d95202d984d9322d…`, plus `rpe1` and `nadig_jurkat` added later for §6.

### One documentation inconsistency found

AGENTS.md §4.2's prose (line ~259) says whitening takes macro PDS *"0.7378 -> 0.8078 (+9.5%)"*.
§4.2's **own table** four lines earlier, exp 019's RESULTS.md, and this reproduction all say
**0.7978** (+8.1%). The cosine half of that sentence (0.0435 → 0.0295) reproduces at
**0.0435 → 0.0293**. The `0.8078` appears to be a typo for `0.7978`; nothing downstream in this
experiment depends on it, but it is the kind of number a later agent would try to hit.

---

## 3. Stage A — the family is already at the optimum of the direct objective

Before fitting anything, run the hand-chosen family through the objective the fit would use.
The donor objective is strictly held-line-safe: for fold `H`, one non-held donor's `gene_norm`'d
effect is the query and the *other* non-held donor's raw effect for the same perturbation is
the answer, over the perturbations both measure **minus every VCC target gene**. Gallery size
is 256, deliberately matched to the panel's 259. Dev is a 20% perturbation holdout.

Macro over the three folds (`work/stage_a.json`):

| gamma | panel pds (held LINE) | panel cosine | donor train pds | **donor dev pds** |
|---|---:|---:|---:|---:|
| 0 | 0.7210 | **0.0435** | 0.6804 | 0.6883 |
| 0.35 | 0.7815 | 0.0348 | 0.7269 | 0.7341 |
| 0.55 | 0.7963 | 0.0310 | 0.7405 | 0.7466 |
| **0.65** | **0.7978** | 0.0293 | **0.7427** | **0.7483** |
| 0.75 | 0.7971 | 0.0277 | 0.7422 | 0.7473 |
| 1.1 | 0.7849 | 0.0230 | 0.7283 | 0.7322 |
| 1.3 | 0.7749 | 0.0206 | 0.7165 | 0.7200 |

Refined (`work/wrapup.json` → `fine_gamma`):

| gamma | 0.55 | 0.58 | 0.60 | 0.62 | **0.65** | 0.68 | 0.70 | 0.75 | 0.80 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| panel pds | 0.79628 | 0.79724 | 0.79726 | 0.79735 | 0.79776 | **0.79784** | 0.79779 | 0.79713 | 0.79632 |
| donor dev pds | 0.74662 | 0.74745 | 0.74783 | 0.74811 | **0.74827** | 0.74823 | 0.74810 | 0.74731 | 0.74611 |
| panel cosine | 0.03095 | 0.03044 | 0.03010 | 0.02977 | 0.02928 | 0.02880 | 0.02849 | 0.02773 | 0.02699 |

**argmax(donor dev) = 0.65. argmax(panel) = 0.68, which beats 0.65 by 8e-5 — three orders of
magnitude inside the fold noise.** Two objectives that share no data agree on the family's one
parameter. `gamma = 0.65` is not a lucky borrowing from word-embedding post-processing; it is
where the retrieval objective puts it.

**Falsifier 3 passes on the family itself.** Cosine is maximised at `gamma = 0` over the full
grid — cosine would have rejected the whole +0.060 lever — and inside the fine window its
argmax is 0.45, still on the wrong side. Cosine is reported everywhere below and is never
allowed to select.

---

## 4. Stage B — the fit: learns on donors, loses on the line

`W(M) = W_prod (I + E_k M E_k^T)`, `M = U V^T + diag(a)`, `E_k` the top-k eigen-subspace of the
same donor covariance that defines `W_prod`. `V` is initialised at zero so **`M = 0` is
bit-identical to the incumbent** and the epoch-0 row is 0.79776 by construction, not by luck.
In-batch listwise softmax cross-entropy, Adam, `0.5 * reg * ||M||_F^2`. Every hyper-parameter
is selected on **donor dev**; the panel is scored once per configuration and selects nothing.

Because `E_k` is orthonormal the whole 7,760-dimensional problem collapses to `k` dimensions
plus two per-batch constants (`fitlib.py` docstring), which is what made 23 configurations, a
capacity ladder, a 12-seed permutation null and two line-generalisation studies affordable.

### Learning curve (fold K562-gwps, k=256, rank=8; `work/stage_b_main.json`)

| epoch | loss | train pds | dev pds | **panel pds** | panel cosine |
|---|---:|---:|---:|---:|---:|
| 0 | — | 0.6829 | 0.6808 | **0.80451** | 0.0302 |
| 1 | 5.2021 | 0.7069 | 0.7032 | 0.7901 | 0.0320 |
| 2 | 5.1455 | 0.7111 | 0.7057 | 0.7890 | 0.0325 |
| 5 | 5.1349 | 0.7119 | 0.7066 | 0.7874 | 0.0327 |
| 10 | 5.1332 | 0.7122 | 0.7065 | 0.7885 | 0.0329 |
| 20 | 5.1326 | 0.7151 | **0.7096** | 0.7920 | 0.0329 |
| 40 | 5.1327 | 0.7143 | 0.7082 | 0.7886 | 0.0325 |

The damage is done in **epoch 1** and never recovers. Note the third falsifier firing again in
miniature: the arm **raises** panel cosine 0.0302 → 0.0329 while **lowering** panel pds
0.8045 → 0.7886. A harness scoring cosine would have shipped this.

### Parameter count against sample size — overfitting is *not* the failure

27,648 batched fit queries in this fold; 70,496 labelled retrieval problems across all three
folds (fit + dev), each a 256-way ranking label.

| k | rank | params | queries/param | d train | **d dev** | **d panel** |
|---|---|---:|---:|---:|---:|---:|
| 64 | 0 (diag only) | 64 | 432 | +0.0007 | +0.0008 | −0.0009 |
| 64 | 8 | 1,088 | 25.4 | +0.0109 | +0.0100 | −0.0035 |
| 64 | 64 | 8,256 | 3.3 | +0.0138 | +0.0106 | −0.0050 |
| **256** | **8** | **4,352** | **6.4** | **+0.0134** | **+0.0113** | **−0.0055** |
| 256 | 32 | 16,640 | 1.7 | +0.0164 | +0.0112 | −0.0060 |
| 256 | 256 (full, no diag) | 131,072 | 0.21 | +0.0191 | +0.0107 | −0.0052 |
| 1024 | 16 | 33,792 | 0.8 | +0.0272 | +0.0102 | −0.0055 |

The 1,088-parameter model has 25 labelled queries per parameter, generalises cleanly to unseen
perturbations (+0.0100 dev), and still loses the line. Shrinking capacity by 2,000x (the
64-parameter diagonal) does not buy a positive; it buys a no-op. **Classical overfitting is not
what is happening** — the mandate's stated hypothesis is falsified, and something worse is.

### All 23 configurations

Adding a temperature / regularisation / learning-rate grid (`work/stage_b_tau.json`), the
23-config population:

- `d dev` ranges **+0.0000 to +0.0204**; `d panel` ranges **−0.0167 to +0.0016**
- `d panel` mean **−0.00485**, sd **0.00537**; **21 of 23 negative**
- **corr(d dev, d panel) = −0.227**; corr(d train, d panel) = −0.295

The two positive cells (+0.0016, +0.0010) are the most heavily regularised configs, i.e. the
ones closest to `M = 0`, and both are +1.2 sd above a population mean of −0.005 in 23 draws.
Selecting on donor dev — the only legitimate rule — picks `d dev = +0.0204`, `d panel = −0.0070`.
**Donor-side generalisation carries no information about held-line performance, and what
little it carries points the wrong way.**

### The pre-committed arm, with an interval

`k=256, rank=8, diag, tau=0.05, reg=1e-3, lr=0.01, gamma=0.65`, epoch chosen by max donor dev.
Paired bootstrap over panel perturbations, 20,000 replicates (`work/wrapup.json`):

| fold | base | arm | delta |
|---|---:|---:|---:|
| K562-gwps | 0.80451 | 0.79204 | −0.01247 |
| Orion-HCT116 | 0.80488 | 0.80074 | −0.00415 |
| Orion-HEK293T | 0.78389 | 0.78389 | +0.00000 (dev never improved; early stop kept the incumbent) |
| **macro** | **0.79776** | **0.79222** | **−0.00554** |

**95% CI [−0.01049, −0.00077], P(positive) = 0.0118.** The loss is significant, not noise. The
fold that loses most (K562-gwps, −0.0125) is the fold with the **most** training data — its
training pair is Orion↔Orion, 34,818 queries — and the fold that refuses to move is the one
where donor dev falls from epoch 1. More donor data buys a larger donor gain and a larger
held-line loss.

**Mandate abort criterion met: held-line pds does not beat the incumbent in any of 23
configurations.**

---

## 5. Falsifiers

### 5.1 Permuted labels — the real fit loses to doing nothing

Identical pipeline; only the query→answer map inside each batch is replaced by a derangement.
3 real seeds, 12 permuted seeds (`work/stage_c.json`):

| | real | permuted | permuted ≥ real |
|---|---:|---:|---:|
| **d dev** (donor) | **+0.01113** (sd 0.00022) | +0.00007 (sd 0.00012) | **0 / 12** |
| **d panel** (held line) | **−0.00614** (sd 0.00071) | +0.00022 (sd 0.00045) | **12 / 12** |

Both halves matter. On the donor objective the real fit beats every deranged twin by 50x its
own spread — **the learning signal is genuine**, this is not the Ahlmann-Eltze failure where 20
of 20 permutations won. On the held line the deranged twins — which learn nothing and leave
`M ≈ 0` — beat the real fit **12 out of 12**. A fit that has demonstrably learned the right
thing about donors is strictly worse than not fitting at all on the line that is scored.

### 5.2 The panel-basis artefact, reproduced on purpose

§4.2 point 1's signature had to be detectable before the Stage B negative could be trusted, so
it was reproduced as a positive control (`work/stage_d.json`).

**D2a — covariance fitted on the 259 panel predictions instead of the donor bank:**

| gamma | 0 | 0.35 | 0.65 | 0.9 | 1.3 | 2.0 | 3.0 |
|---|---:|---:|---:|---:|---:|---:|---:|
| macro panel pds | 0.7210 | 0.7520 | 0.7721 | 0.7838 | 0.7942 | **0.7977** | 0.5242 |

Monotone all the way to `gamma = 2.0` — three times past the donor basis's interior optimum,
with no peak — and then a numerical collapse when the floored near-zero eigenvalues of a
rank-259 basis get amplified. §4.2 point 1 confirmed.

**D2b vs D1 — capacity ladder, adapter fitted ON the panel (leakage) vs on donors (legitimate):**

| rank | panel-FIT (leakage) | donor-FIT panel | donor-FIT dev |
|---|---:|---:|---:|
| 0 | 0.8008 | 0.7821 | 0.7375 |
| 4 | 0.9158 | 0.7900 | 0.7539 |
| 16 | 0.9669 | 0.7906 | 0.7546 |
| 32 | 0.9826 | 0.7909 | **0.7556** |
| 64 | 0.9921 | 0.7878 | 0.7538 |
| 128 | 0.9961 | 0.7910 | 0.7534 |
| 256 | **0.9967** | 0.7880 | 0.7511 |

Fit on the thing you score and capacity buys the answer, monotonically, to 0.9967. The
legitimate donor fit does the opposite: an **interior optimum at rank 32** on donor dev and a
flat, always-below-baseline panel. **The harness can tell the two apart**, which is what makes
the negative worth believing.

### 5.3 Cosine is allowed to fall, and it does the reverse

Reported at every step above. Family sweep: cosine argmax `gamma = 0` (0.0435), pds argmax
0.65 — cosine rejects the project's largest gain. Learned arm: cosine **rises** 0.0302 → 0.0329
while pds **falls** 0.8045 → 0.7886. Cosine selected nothing anywhere in this experiment.

---

## 6. Stage E / F — where it actually breaks: the cell line

Stage B has one structural escape route. After a LODO fold holds out its scored line, **exactly
one line pair remains**, so nothing in the fit ever asks the correction to survive a change of
cell line. A pair-specific correction would look perfect on a perturbation holdout and still
die on the panel — precisely the observed pattern.

RPE1 and Nadig-Jurkat are never a LODO fold, so each fold can have **four** training lines and
six line pairs, and an inner **leave-one-training-LINE-out** becomes possible with the scored
line still untouched. (K562-essential is deliberately excluded: it is the same cell line as
`k562_gwps`, which *is* a scored fold.) Genes a line does not measure are stored as exact
zeros, so each pair's retrieval problem is posed on its shared sub-axis.

**Stage E — single-donor queries, 12 (fold, held-line) cells** (`work/stage_e.json`):

| rung | what is unseen | macro delta |
|---|---|---:|
| train | nothing | +0.01306 |
| dev | perturbations | **+0.00904** |
| **unseen LINE** | the gallery **cell line** | **−0.00618** |
| panel | gallery line **and** VCC panel targets | −0.00921 |

**Stage F — two-predictor queries, structurally identical to the panel query** (query = equal
weight mean of two donor lines, gallery = a third, unseen line), 12 cells
(`work/stage_f.json`):

| rung | macro delta |
|---|---:|
| train | +0.00748 |
| dev (perturbations) | +0.00082 |
| **unseen LINE** | **−0.00880** |
| panel | **−0.01273** (negative in **12 / 12** cells) |

This closes the last mundane explanation. Stages B/C/E train on single-donor queries while the
panel query is a two-donor average, and averaging is *not* neutral here — exp 019's whole
result is that whitening inverts the donor-weighting effect. Stage F builds the matched
problem, and the correction still fails to cross the line boundary, by more.

**The gain is real, it generalises across perturbations, and it stops dead at the cell line —
and it stops at a donor line by the same −0.006 to −0.009 that it stops at the scored line.**
The panel is not a special or adversarial target. The learned correction is line-pair-specific
structure, and with three lines total (§8.1: information-limited, κ² ≈ 0.34 shared across
lines, the `1/k` law needs ~10) there is not enough line diversity to fit a line-general one.
This is the same wall as exp 159's transfer ceiling and exp 166's "no further panel-covering
CRISPRi screen exists in public data", reached from a different direction.

---

## 7. Phase 2 — the six-member price is not owed, and here is the arithmetic

Rule 12 requires every member be priced or pinned. There is no arm to price: the pre-committed
configuration is **−0.00554 raw macro pds on LODO with P(positive) = 0.0118**, and no
configuration selected by any legitimate rule was positive. The composed gate (rule 13) is
therefore decided without spending PLOLO, `cell-eval2`, or a slot:

```
predicted = PLOLO d_pds/6 * 1.011 + 096 d_mse/6 (if amplitude/spread/phi moves)
          ≈ (-0.00554)/6 * 1.011                          = -0.00093
SHIP only if predicted > +0.0019   →  shortfall  -0.0028, before any unpriced channel
```

`mse` is **not** computed and cannot rescue this. Rule 12's exemption explicitly does **not**
cover anything touching the operator, and a19 lost **−0.01352** on `mse` from exactly this
class of change — a whitener ridge reshaping the effect vector. An unmeasured `mse` term here
can only make the composite worse or leave it unchanged. `nmae`, `reach`, `fid` and `jac` are
not locally priceable at all (rule 13, exp 162) and none of them can convert a negative `pds`
arm into a shippable one.

For calibration: the LODO macro `pds` gain a shippable arm of this class would need is
**+0.0113 raw**, which is 2.3x this screen's bootstrap CI half-width (0.00486) — **the screen has the
power to detect a qualifying arm.** It measured the opposite sign at P(positive) = 0.012.

Nothing was submitted, no counts were generated, `cell-eval2` was not run, and nothing outside
`experiments/195-retrieval-transform/` was written.

---

## 8. What this closes, and the one thing it does not

**Closed.**

1. *"`gamma = 0.65` was borrowed, not derived — fit it instead."* The direct retrieval
   objective, on data that never sees the scored line, puts the optimum at **0.65**. There is
   nothing to gain inside the family. (§3)
2. *"080 was data-starved at 600 targets."* 70,496 queries, 117x more, reproduce the same
   negative with a significant interval. (§4)
3. *"080's rank ≤ 16 was too small / its basis was wrong."* Capacity from 64 to 131,072
   parameters, three subspace sizes, four temperatures, three regularisations: 21 of 23
   negative, and `d panel` is *negatively* correlated with `d dev`. (§4)
4. *"The single-donor training query is the wrong problem."* The matched two-predictor problem
   fails harder. (§6)
5. *"Maybe it is overfitting and needs more regularisation."* At 25 queries per parameter with
   clean dev generalisation it is not overfitting; the most regularised configs are simply
   nearer the identity. (§4)

**Not closed, and out of this experiment's lane.** The **diagonal-only** arm — a free
eigenvalue map inside the existing family — was measured only as a byproduct and is
**193-whitener-shrinkage's territory**, not mine. For the record and to save 193 a run: on this
LODO harness a free diagonal on the top-64 eigen-directions fitted to the ranking objective is a
no-op (`d dev +0.0008`, `d panel −0.0009`, 64 parameters), and on the top-256 it is harmful
(`d panel −0.0157` at rank 0 in the D1 ladder). That is consistent with §3's finding that the
one-parameter spectral map is already at its optimum, but it is **not** a Marchenko–Pastur or
shrinkage result and should not be read as one. 192 (call set) and 194 (donor tables) are
untouched.

**The strategic reading.** `pds` is 5/6 of the score and its geometry lever is *finished*. The
transform is at the optimum of its own objective, and the residual structure a fit can find is
line-specific and does not transfer. Under §8.1's split this is a **DIRECTION**-channel result
and it says the same thing 159 and 166 said: the binding constraint is three cell lines, not
the shape of the operator applied to them. Operator work should stop; line diversity (184's
Tahoe transfer) and the ORDERING channel are where the remaining headroom is.

---

## 9. Reproduce

```bash
uv run python experiments/195-retrieval-transform/build_cache.py      #  7 s, 5.9 GB peak
uv run python experiments/195-retrieval-transform/phase0.py           # 116 s, 6.6 GB  MUST PASS
uv run python experiments/195-retrieval-transform/test_fitlib.py      #  gradient / algebra checks
for f in k562_gwps orion_hct116 orion_hek293t; do
  uv run python experiments/195-retrieval-transform/prep.py --fold $f
done                                                                  # 11 s each, 7.3 GB peak
uv run python experiments/195-retrieval-transform/stage_a.py          # 13 s
uv run python experiments/195-retrieval-transform/stage_b.py --grid main --epochs 40 \
  --out experiments/195-retrieval-transform/work/stage_b_main.json    # ~2 min
uv run python experiments/195-retrieval-transform/stage_b.py --grid tau --epochs 25 \
  --out experiments/195-retrieval-transform/work/stage_b_tau.json     # ~7 min
uv run python experiments/195-retrieval-transform/stage_c.py --seeds 12 --epochs 25   # ~9 min
uv run python experiments/195-retrieval-transform/stage_d.py --epochs 25              # ~5 min
uv run python experiments/195-retrieval-transform/extend_cache.py
for f in k562_gwps orion_hct116 orion_hek293t; do
  uv run python experiments/195-retrieval-transform/prep_multi.py --fold $f
done                                                                  # ~3 min total
uv run python experiments/195-retrieval-transform/stage_e.py --k 256 --rank 8 --epochs 25
uv run python experiments/195-retrieval-transform/stage_f.py          # 37 s
uv run python experiments/195-retrieval-transform/wrapup.py           # 37 s, 8.3 GB peak
```

Peak RSS across the whole experiment **8.27 GB** (against a 12 GB ceiling), total wall clock
**~30 min**. Logs and JSON in `work/`. The projected-query arrays under `work/prep_multi/`
(7.7 GB, derived, rebuilt by the `prep_multi.py` loop above) were deleted after the run; the
operators under `work/operators/` are kept because their sha256 are the rule-11 record.

# 318 — Is PER-TARGET magnitude redistribution a real lever, separate from the uniform amplitude axis?

**Date:** 2026-09-08 · preregistration `PREREGISTRATION.md`, **two digests both recorded in
`work/PREREGISTRATION.sha256`**: v1 `cfefadcb…` at 01:01:34Z, v2 `76aac610…` at 01:03Z. The
single edit between them is to §2 only — every arm renormalised to one whole-matrix effect
norm, because at the first spelling the rule's anti-correlation with the row norm would have
made rule-vs-shuffle a 1.4 % amplitude contrast instead of an assignment contrast. The
hypothesis, the frozen predictions and the decision rule are unchanged, and **both digests
predate the first worker launch at 01:04:22Z, every scored build, and every score.** · instrument: the real `cell-eval2` 0.16.0 (`--preset vcc2026`, DE backend `scanpy`)
on `data/harness/h1-test.h5ad`, six disposable EC2 workers · **no submission built, nothing
under `submissions/` written, no `vcc` command run, no local `cell-eval2`, local peak RSS
4.18 GiB.**

> **318 (per-target magnitude redistribution): NO-GO — and NOT SEPARABLE FROM ITS SHUFFLED
> CONTROL.** The replayed profile, isolated at fixed direction and fixed whole-matrix norm,
> buys **+0.000093 overall — 1.0 % of the +0.0089 gate** — at **z = −0.06** against three
> entropy-matched shuffles of its own multipliers (best shuffle +0.000543 = 6 % of the gate).
> **Six times the measured redistribution reads −0.001101**: the axis is negative and steepens
> with dose. And the premise does not hold — **half of the +0.00075 is `pds`, which per-target
> magnitude cannot move**, and **the "exactly shape-only" parity arm redistributes delivered
> magnitude MORE (sd(log) 0.0485) than the "distortion" it is supposed to isolate (0.0329)**,
> because the norm match is enforced one stage before an operator whose per-row gain spans
> 0.591–1.286.

---

## 1. Step 1 — the distortion, measured where it acted

### 1.0 The mandate's numbers are FIXTURE numbers. Here are the real ones.

The mandate opens on "delivered per-target row norms 0.806–2.155, whole matrix +6.0 %, where
the parity arm holds every row at 1.0000000". Those come from `experiments/316-parity-arm/
work/test316.json`, which measured them on **313's synthetic fixture** — 20 targets × 60 genes,
a third of the columns "panel", against production's 1.62 %.

`profile318.py` measured the same quantity on the real thing: the **delivered log-effect** of
every (context, target) block of the two 12.6 GB predictions that produced the board points,
against the challenge's own control profiles, on all 300 panel targets × 3 contexts.

    p(t,c) = the 400 emitted count vectors of target t in context c, summed
    effect = log1p(50000 * p / p.sum())  -  b_ctrl(c)          # what every DE member reads
    m(t,c) = || effect_shipped || / || effect_parity ||

| | fixture (316) | **real panel (318)** |
|---|---|---|
| `m` median | — | **1.011353** |
| `m` range | 0.806 – 2.155 | **0.9321 – 1.2349** |
| sd(log m) | — | **0.032885** |
| whole-matrix norm ratio | **1.0600** | **1.0085497** |

**The distortion is 6× smaller in spread and 7× smaller in mean than the mandate's brief
states.** Everything below is priced on the measured profile.
(`work/profile318.json`, 900 blocks, peak RSS 1.58 GiB, 13.6 s.)

### 1.1 The profile is STRUCTURED — three independent ways

| test | result |
|---|---|
| reproducibility across the 3 contexts (same effect matrix, different controls and library draws) | Pearson(log m) **0.9992 / 0.9992 / 0.9995**, ICC ≈ **1.000** |
| coherence as a DOSE family: `VCC_NORMREL_EXP` −1.5 vs −3.0, delivered log ratio vs the incumbent | Spearman **+0.9931**, Pearson +0.9805 |
| predictability from build-time observables (300 targets, context means) | see below |

So the answer to *"is the profile structured or effectively arbitrary"* is **structured, and
not marginally**: it is a deterministic, dose-coherent function of the target.

**What predicts `log m(t)`** (Spearman over 300 panel targets, `work/predict318.log`):

| predictor | Spearman | Pearson |
|---|---:|---:|
| the parity arm's OWN delivered magnitude ratio | **−0.5789** | −0.4624 |
| \|log m\| vs the parity row's rotation from the incumbent (1 − cos) | **+0.5717** | +0.6396 |
| the incumbent's own delivered row norm | **−0.3044** | −0.2376 |
| the parity row's rotation from the incumbent (signed) | −0.1246 | +0.0690 |

Read together with §1.2 these say one thing: **`m(t)` is the amount by which the operator's
per-row gain changes when the row rotates.** Where the row moved more, its delivered magnitude
moved more (|log m| vs rotation +0.57); and it is largely the parity arm's own redistribution
being partly undone (−0.58). The row-norm relation (−0.30: a smaller predicted response is
amplified) is the one that is directly implementable at build time, and it is what §2's rule
arm uses.

*A note on the two row-norm numbers.* Against the PARITY arm's own delivered row norm the
Spearman reads −0.414; against the INCUMBENT's it reads **−0.304**. The first is partly
circular (the parity norm is the denominator of `m`), so −0.304 is the honest figure and the
one quoted here. The preregistration cites −0.414; the sign, the monotone direction and the
rule built from it are unaffected, and the rule was frozen before any score existed.

**The single-donor control — 6 targets the arm provably cannot touch.** `census.json` records
6 of the 300 panel targets carried by exactly ONE production donor (`C5orf22`, `CDK8`,
`GLYR1`, `IL10RB`, `SPATA20`, `UNK`); on those, a per-target donor weighting is the identity
whatever the exponent. Their delivered magnitude still moves, and it moves by a **constant**:

| | log m (shipped/parity) | log r (parity/incumbent) |
|---|---|---|
| the 6 single-donor targets | +0.01065 … +0.01153 | −0.03656 … −0.03811 |
| all 300, mean \|·\| | 0.02536 | — |
| the other 294, mean \|·\| | 0.02565 | — |

So the profile splits into a **global** part — a fixed +1.1 % (shipped vs parity) and −3.7 %
(parity vs incumbent), the whitening's global preserve scale and the emitter renormalisation
responding to the *whole* matrix — and a **per-target** part carried by the 294 targets the
weighting can act on. Neither part was chosen by anyone; both are downstream of a construction
decision taken two stages earlier.

### 1.2 🔴 The confound cannot be broken by the parity arm, because the parity arm is not magnitude-neutral either

316's KA-3 certifies the parity arm "exactly shape-only": max |delivered row norm / incumbent
− 1| = **4.44e-16**. That statement is true, and it is made on the **pre-whiten** effect
matrix — one stage before the operator. Measured on what the scorer reads
(`work/family318.log`, four board-scored 12.6 GB predictions on disk):

| arm | board overall | delivered ratio vs incumbent: median | range | **sd(log)** | whole-matrix | median row cos vs incumbent |
|---|---:|---:|---|---:|---:|---:|
| `cd4svd5` incumbent | 0.2160944 | 1 | — | 0 | 1 | 1 |
| **parity (316)** | 0.2169733 | 0.972819 | 0.9234 – 1.3067 | **0.048539** | 0.9965743 | 0.98708 |
| **shipped (312, exp −1.5)** | 0.2177189 | 0.990479 | 0.9132 – 1.3328 | **0.044288** | 1.0050947 | 0.98762 |
| exp −3.0 (312b) | 0.2164511 | 0.980194 | 0.8549 – 1.4612 | 0.070191 | 1.0033888 | 0.95570 |

**The parity arm redistributes delivered magnitude MORE than the shipped arm does
(sd 0.0485 vs 0.0443), and both redistribute more than they differ from each other (0.0329).**
The mandate's premise — "the same weighting, made exactly shape-only" versus "an accidental
per-target norm distortion" — does not survive contact with the delivered matrices. There is no
magnitude-neutral member of this family on the board.

**Why, measured** (`mechanism318.py`, the production operator, no donor model, KA row at
0.000000):

| rotation applied at an IDENTICAL pre-whiten norm | delivered norm ratio: median | range | sd(log) |
|---|---:|---|---:|
| **KA — cos 1.000000 (identity)** | **1.000000** | **1.000000 – 1.000000** | **0.000000** |
| cos 0.999 (2.6°) | 0.984705 | 0.7888 – 1.5030 | **0.118196** |
| cos 0.98702 (the parity arm's own median) | 0.997418 | 0.7021 – 1.7840 | **0.159239** |
| cos 0.95 | 0.987083 | 0.7003 – 1.6203 | 0.158178 |

The operator's per-row gain runs **0.591 – 1.286, sd(log) 0.132** across the 300 challenge
rows. It is strongly anisotropic, so **a pre-whiten norm match is not a delivered norm match**:
a 2.6° rotation at an exactly matched pre-whiten norm already redistributes delivered magnitude
3.6× more than the entire shipped-vs-parity difference. (The rotation here is into a *random*
orthogonal direction and is therefore an upper bound; the real rotations lie inside the donor
subspace, which is why the observed 0.0485 is smaller. The direction of the argument does not
depend on that.)

### 1.3 What the +0.00075 is actually made of

Raw and scaled members of the two board points, from their captured records
(`data/board-entries/`), scaled through 317 §1.4's fitted map:

| member | raw Δ (shipped − parity) | scaled Δ | share of the total |
|---|---:|---:|---:|
| `pds` | +0.001014 | **+0.002246** | 50 % |
| `mse` | −0.003788 | **+0.003988** | 89 % |
| `nmae` | +0.001855 | **−0.003075** | −69 % |
| `fid` | +0.000268 | +0.000901 | 20 % |
| `reach` | +0.000629 | +0.000710 | 16 % |
| `jac` | −0.000005 | −0.000015 | 0 % |
| **Σ/6** | | **+0.000793** | (board: **+0.0007456**) |

**Half of it is `pds`, and `pds` is the one member a per-target magnitude cannot pay.** `pds`
is a per-target rank statistic; a positive per-target scalar cannot reorder a target's own
genes. 317 measured this on the real scorer: a **uniform ×1.20** returned
`Δ raw pds = +0.000000`, and the whole 0.70 → 2.00 span moved it by 3.5e-3 (count flooring
only). §3's `ptamp6` tests it directly at six times the measured redistribution.

---

## 2. The isolated arm — construction

Because §1.2 shows there is no magnitude-neutral member of the board family, the isolation has
to be built rather than borrowed. `build318.py` is `experiments/317-operating-point/
build317.py` with one structural change: `effect_scale` — the delivered log-effect amplitude,
applied in `profile_from_effect`, the exact point 317's uniform amplitude knob acts — becomes a
**vector over targets**. Effect matrix, library draw, emission seed, spread policy, call set
and emitter are byte-identical to `base` in every arm, so the ONLY variable is how magnitude is
distributed across targets. Direction is held exactly fixed: a positive per-target scalar
cannot rotate a row.

`q(t)` is the measured 300-value `log m` profile transferred to h1-test's 85 targets by matched
quantile of the target's own effect row norm, descending — the observable §1.1 found the
profile monotone in, taken at its **strong form** (Spearman = −1 by construction, against the
measured −0.304). That choice can only help the hypothesis. Every arm is then renormalised to
ONE whole-matrix effect norm, **1.0085497 — the real shipped arm's own measured value** — so
`ptunif` carries the whole-matrix component alone and every other contrast is pure
redistribution. A common scalar preserves sd(log m) exactly.

| arm | `m(t)` | median | range | sd(log m) | whole-matrix (effect) |
|---|---|---:|---|---:|---:|
| `base` | 1 | 1.000000 | — | 0 | 1.0000000 |
| `ptunif` | flat | 1.008550 | — | 0 | 1.0085497 |
| `ptrule` | `exp(q)` | 1.021119 | 0.965451 – 1.144537 | 0.030964 | 1.0085497 |
| `ptshuf1` | permuted, seed 1 | 1.007275 | 0.952362 – 1.129020 | 0.030964 | 1.0085497 |
| `ptshuf2` | permuted, seed 2 | 1.004342 | 0.949589 – 1.125732 | 0.030964 | 1.0085497 |
| `ptshuf3` | permuted, seed 3 | 1.004740 | 0.949965 – 1.126178 | 0.030964 | 1.0085497 |
| `ptamp6` | `exp(6q)` | 1.057226 | 0.755254 – 2.096471 | 0.185783 | 1.0085497 |

### 2.1 Gates — every one of them, on every worker

| gate | claim | measured |
|---|---|---|
| 1 | a15 and a19's operator anchors (221's) | **PASS**, every scalar ≤ 1e-11 relative |
| 2 | the h1 effect is 170's artifact | `‖d‖F/‖E‖F` **6.462e-15**, `1−min row cos` 2.22e-16 |
| 3 | the CD4 fraction is linear in `f` | rel **2.22e-16** |
| 4 | a per-target magnitude arm moves EVERY block | **0/85** digest-identical to `base`, 6/6 arms |
| **5** | **`base` IS 221's `rc_callset`** | sha256 **`048b9286030ff0124514134b55d935fdb79a4d002cc04d590d83676b32553c77`** on **6/6 workers** |
| **6** | the shuffles are ENTROPY-MATCHED: same multiset up to one scalar | dev **≤ 2.0e-16**, sd(log) diff ≤ 1.4e-17 |
| **7** | the emitter DELIVERS the redistribution it was asked for | see below |

Gate 7 matters because `effect_scale` acts on the log-effect *before* `profile_from_effect`
renormalises, so the delivered ratio need not equal `m(t)`. Measured on the emitted blocks
against the same control profile the scorer reads:

| arm | asked sd(log m) | **delivered** sd(log) | rank agreement | delivered whole-matrix |
|---|---:|---:|---:|---:|
| `ptrule` | 0.030964 | **0.030728** | **0.9987** | 1.0097237 |
| `ptshuf1` | 0.030964 | 0.030732 | 0.9996 | 1.0094181 |
| `ptshuf2` | 0.030964 | 0.030776 | 0.9997 | 1.0085944 |
| `ptshuf3` | 0.030964 | 0.030711 | 0.9996 | 1.0086917 |
| `ptamp6` | 0.185783 | **0.184742** | **0.9998** | 1.0157383 |

The emitter delivers 99.2–99.4 % of the requested spread at rank agreement ≥ 0.9987. One
residual to carry: at dose 6 the log1p renormalisation pushes the *delivered* whole-matrix
norm to 1.01574 against 1.00855 for the others, so `ptamp6` carries ≈ +0.7 % more delivered
amplitude than its uniform twin. `ptunif` measures the price of exactly that kind of uniform
excess, so it can be subtracted rather than assumed away (§4).

---

## 3. The six members, on the real scorer

### 3.1 KA-1 — the baseline is 221's recorded row

`base` is not re-scored. `cell-eval2 run` reads only (prediction, reference, version, preset);
all four are pinned identically — prediction sha256 `048b9286…` (gate 5, 6/6 workers),
reference sha256 **`06c6baea8c62ef3d38cacdbf58e104908eb93f0ef3b25aca56b86d28240f9692`**
(printed by 317's KA-1 worker, by 317's `amp120` worker and by every 318 worker), `cell-eval2`
0.16.0, `--preset vcc2026`. 317 §1.2 proved the row reproduces bit for bit from the same object
on a fresh worker. The reused row against 221's record: **max |diff| 4.76e-13** (the residual is
the 12-digit truncation of the recorded constants). **KA-1 PASS.**

### 3.2 RAW six members, h1-test, 85 targets

| arm | `pds` | `mse` | `nmae` | `fid` | `reach` | `jac` |
|---|---:|---:|---:|---:|---:|---:|
| **base** | 0.873810 | 1.015978 | 0.960046 | 0.495058 | 0.226704 | 0.074204 |
| `ptunif` | 0.874230 | 1.024439 | 0.959843 | 0.495168 | 0.226800 | 0.074203 |
| `ptrule` | 0.874790 | 1.026352 | 0.959459 | 0.495106 | 0.226670 | 0.074204 |
| `ptshuf1` | 0.874370 | 1.025232 | 0.959939 | 0.495130 | 0.226713 | 0.074201 |
| `ptshuf2` | 0.874790 | 1.024514 | 0.959736 | 0.495141 | 0.228507 | 0.074199 |
| `ptshuf3` | 0.874370 | 1.025040 | 0.960179 | 0.495092 | 0.226702 | 0.074199 |
| `ptamp6` | 0.873389 | 1.036929 | 0.958286 | 0.495047 | 0.226866 | 0.074213 |

Deltas vs `base`, and the same deltas priced in board overall points through 317 §1.4's fitted
map at transfer factor 1.0 (the optimistic reading):

| arm | d`pds` | d`mse` | d`nmae` | d`fid` | d`reach` | d`jac` | **d overall** |
|---|---:|---:|---:|---:|---:|---:|---:|
| `ptunif` | +0.000420 | +0.008461 | −0.000204 | +0.000110 | +0.000095 | −0.000000 | **−0.001194** |
| `ptrule` | +0.000980 | +0.010374 | −0.000587 | +0.000048 | −0.000034 | +0.000000 | **−0.001276** |
| `ptshuf1` | +0.000560 | +0.009254 | −0.000108 | +0.000072 | +0.000009 | −0.000003 | **−0.001346** |
| `ptshuf2` | +0.000980 | +0.008536 | −0.000310 | +0.000083 | +0.001803 | −0.000005 | **−0.000667** |
| `ptshuf3` | +0.000560 | +0.009062 | +0.000133 | +0.000034 | −0.000002 | −0.000005 | **−0.001403** |
| `ptamp6` | −0.000420 | +0.020951 | −0.001760 | −0.000011 | +0.000161 | +0.000009 | **−0.003316** |

**Every arm in the design is negative, and `mse` is the whole story.** (`mse` is a per-panel
aggregate; raw up is worse.)

Paired per-target intervals (unit: **one perturbation on h1-test**, paired against `base` —
same effect matrix, same library draw, one multiplier vector moved; 10,000-resample percentile
bootstrap, seed 0; `*` = 95 % excludes zero; n = 85 `pds`/`fid`/`jac`, 58 `nmae`, 82 `reach`;
`expr_mse_unbiased_capped_norm` has no per-target value and gets no interval):

| arm | `pds` | `nmae` | `reach` |
|---|---|---|---|
| `ptunif` | +0.00042 [−0.00084, +0.00168] | **−0.00020 [−0.00034, −0.00006] \*** | +0.00010 [−0.00023, +0.00050] |
| `ptrule` | +0.00098 [−0.00056, +0.00294] | **−0.00059 [−0.00092, −0.00023] \*** | −0.00003 [−0.00025, +0.00016] |
| `ptshuf1` | +0.00056 [−0.00042, +0.00168] | −0.00011 [−0.00051, +0.00028] | +0.00001 [−0.00025, +0.00028] |
| `ptshuf2` | +0.00098 [−0.00042, +0.00252] | −0.00031 [−0.00076, +0.00009] | +0.00180 [−0.00012, +0.00542] |
| `ptshuf3` | +0.00056 [−0.00070, +0.00182] | +0.00013 [−0.00041, +0.00061] | −0.00000 [−0.00013, +0.00014] |
| `ptamp6` | −0.00042 [−0.00266, +0.00140] | **−0.00176 [−0.00312, −0.00040] \*** | +0.00016 [−0.00009, +0.00041] |

### 3.3 The uniform control does almost all of the work

`ptunif` — a **flat** ×1.008550, zero redistribution, the whole-matrix component alone, which
is 317's axis — costs **−0.001194 overall by itself**. That is **94 % of `ptrule`'s −0.001276**.
Its raw `mse` cost is +0.008461 against the +0.0068 the preregistration predicted from 317's
own amplitude curve (1.24× the prediction, same sign, same order): 317's uniform axis is
reproduced independently at a dose 24× smaller than its smallest rung.

Subtracting it — for each arm, `corrected = Δ(arm) − k·Δ(ptunif)` with
`k = (delivered whole-matrix norm − 1)/(ptunif's excess)`, so the residual is what
**redistribution alone** bought:

| arm | delivered whole excess | k | d`pds` | d`mse` | d`nmae` | **d overall, redistribution only** |
|---|---:|---:|---:|---:|---:|---:|
| `ptunif` | +0.008480 | 1.000 | 0 | 0 | 0 | **0** (by construction) |
| `ptrule` | +0.009724 | 1.147 | +0.000499 | +0.000673 | −0.000354 | **+0.000093** |
| `ptshuf1` | +0.009418 | 1.111 | +0.000094 | −0.000143 | +0.000118 | **−0.000021** |
| `ptshuf2` | +0.008594 | 1.013 | +0.000555 | −0.000039 | −0.000104 | **+0.000543** |
| `ptshuf3` | +0.008692 | 1.025 | +0.000130 | +0.000390 | +0.000341 | **−0.000180** |
| `ptamp6` | +0.015738 | 1.856 | −0.001200 | +0.005248 | −0.001383 | **−0.001101** |

(The uniform correction is linear in the excess. The uniform `mse` response is mildly convex —
measured slope 0.998 raw `mse` per unit `s` at +0.85 %, against 317's local 0.8 — so the
linear subtraction slightly **under**-removes the uniform part at `ptamp6`'s larger excess.
That bias makes `ptamp6`'s residual look worse than it is; it cannot manufacture the negative
sign at dose 1, where `k` is within 15 % of 1.)

### 3.4 The separability test, as declared

Declared in the preregistration before any score: *`ptrule` is separable only if
|Δ(`ptrule`) − mean Δ(shuffles)| > 2 × sd(shuffles) on at least one member.*

| member | `ptrule` | shuffle mean | shuffle sd | **z** |
|---|---:|---:|---:|---:|
| **overall (priced)** | −0.001276 | −0.001139 | 0.000410 | **−0.33** |
| overall, uniform-corrected | +0.000093 | +0.000114 | 0.000380 | **−0.06** |
| `pds` | +0.000980 | +0.000700 | 0.000243 | +1.15 |
| `mse` | +0.010374 | +0.008951 | 0.000372 | **+3.83 (worse)** |
| `nmae` | −0.000587 | −0.000095 | 0.000222 | −2.22 (better) |
| `fid` | +0.000048 | +0.000063 | 0.000026 | −0.59 |
| `reach` | −0.000034 | +0.000603 | 0.001039 | −0.61 |
| `jac` | +1.1e-08 | −4.1e-06 | 1.2e-06 | +3.50 (span 5e-6, meaningless) |

**On the composite — the only thing that decides — `ptrule` is at z = −0.33 as built and
z = −0.06 uniform-corrected. It is indistinguishable from a random assignment of its own
multipliers, and it is on the wrong side of the shuffle mean.** The two members where it does
separate are `mse` at z = +3.83, where the *rule* is the **worse** assignment, and `jac`, whose
whole span is 5e-6.

The best pure-redistribution number in the entire design is **`ptshuf2` at +0.000543** — an
arbitrary permutation. That is 29 % of the +0.0019 shape-only gate and 6 % of the +0.0089
amplitude gate, and it sits 1.1 shuffle-sd from its own siblings.

### 3.5 `pds` — the preregistered structural prediction, and the honest reading

Predicted before scoring: |Δ raw `pds`| < 0.0005 at dose 1 and < 0.004 at dose 6.

| | `ptunif` | `ptrule` | `ptshuf1` | `ptshuf2` | `ptshuf3` | `ptamp6` | 317's re-seeded twin |
|---|---:|---:|---:|---:|---:|---:|---:|
| Δ raw `pds` | +0.00042 | +0.00098 | +0.00056 | +0.00098 | +0.00056 | **−0.00042** | **+0.00070 [−0.00056, +0.00210]** |

**The literal 0.0005 threshold is missed at dose 1 — and the substance is confirmed anyway,
three ways.** (i) Every dose-1 value sits inside 317's re-seeded twin's own paired interval,
i.e. inside this surface's emission-noise floor for `pds`; every paired interval above includes
zero. (ii) The values are **quantised and repeat exactly across independent arms** (+0.000980
for both `ptrule` and `ptshuf2`; +0.000560 for both `ptshuf1` and `ptshuf3`): these are a
handful of discrete rank flips in an 85-target rank statistic, not a systematic response.
(iii) Decisively, **six times the redistribution moves `pds` LESS and in the opposite
direction** (−0.00042). A dose that reverses sign as it grows is not a response.

So the mandate's arithmetic can be closed: the board's `pds` +0.002246 scaled — **half of the
+0.00075 attributed to "the distortion"** — cannot be bought by per-target magnitude. It is the
direction difference between applying the weighting inside and outside the norm match.

---

## 4. The price

### 4.1 Against the gate

An amplitude-changing arm must clear **+0.0089** (AGENTS.md §8.4, §3). In raw units of a single
member, at transfer factor 1.0 (317 §4.1): `mse` −0.05073, `nmae` −0.03222, `pds` +0.02411.

| | overall | vs +0.0089 | vs +0.0019 (shape-only) |
|---|---:|---:|---:|
| `ptrule`, the isolated arm as built | **−0.001276** | — (negative) | — |
| `ptrule`, redistribution only | **+0.000093** | **1.0 %** | 4.9 % |
| best of the 3 shuffled controls, redistribution only | +0.000543 | 6.1 % | 29 % |
| **`ptamp6` — six times the measured dose**, redistribution only | **−0.001101** | — (negative) | — |
| `ptunif` — 317's uniform axis at the same whole-matrix norm | −0.001194 | — (negative) | — |

**The power clause fires.** At **6× the redistribution actually present in the shipped arm**
(delivered sd(log) 0.1847 against 0.0307), the composite is **−0.001101** — further from the
gate, not nearer. The axis has a measured sign and it is negative; this is not an
underpowered null.

### 4.2 Is the +0.00075 reachable deliberately?

**No, and the +0.00075 was not a magnitude effect in the first place.** Four measured reasons:

1. **Half of it is `pds`** (+0.002246 of the +0.000793 scaled sum × 6), and `pds` does not
   respond to per-target magnitude on this scorer (§3.5).
2. **The magnitude-sensitive part of the axis is negative when isolated.** `mse` is 89 % of the
   board delta's positive side; here, redistribution costs `mse` (+0.000673 raw at dose 1 for
   the rule, +0.005248 at dose 6) rather than buying it, and the *uniform* component that
   accompanied it on the board costs −0.001194 overall by itself.
3. **The profile carries no assignment information.** Replaying it beats a random permutation
   of its own multipliers by z = −0.06.
4. **Even the profile is only partly reachable.** Out of fold (5 folds, 300 panel targets), a
   build-time rule on the three observables of §1.1 predicts `log m(t)` at **R² 0.333,
   Spearman +0.645** (`work/rule318.log`) — it could reproduce about 58 % of the sd. `ptrule`
   was deliberately built at the **strong form** (Spearman −1 with the row norm) — more
   aggressive than any achievable rule — and still did not separate.

**What the +0.00075 most likely is.** It is the board delta between two constructions that
differ in *direction* (median delivered row cosine 0.99838, min 0.9485) as well as magnitude,
and it is **0.17 × the 0.004502 ordering-noise sd AGENTS.md §3 records between digest-identical
arms**, 0.39 × the shape-only gate and 0.084 × the amplitude gate. Nothing in this experiment
identifies a lever inside it.

### 4.3 What this does NOT say

- It does not locate the board's amplitude optimum. h1-test cannot (317 §2.2): it is scored
  with the H1 union OFF (oracle-leaked, 85/85) and its raw `mse` at the shipped point is 1.016
  against the board's 0.81671. **No absolute `mse` level here is read as a board opportunity**,
  and `ptunif`'s −0.001194 is quoted as *"a uniform +0.85 % costs about a thousandth on this
  surface"*, not as a board price.
- It does not price a redistribution arm for the board: the h1-test → board transfer factor for
  this arm class is unidentified. Factor 1.0 is the optimistic reading and is the one used; any
  factor below 1 shrinks a result that is already negative or at 1 % of the gate.
- It does not say the shipped arm should be withdrawn. It says the +0.00075 between it and the
  parity arm is not a per-target magnitude lever, and that neither arm is magnitude-neutral.
- Every interval is on **85 h1-test targets**, and `reach` is the noisiest member on this
  surface by an order of magnitude (317 §1.5).

---

## 5. Cost and cleanliness

| | |
|---|---|
| instances | 6 × `r7i.2xlarge` on-demand, `us-east-1`, one arm each |
| wall time | 2.680 – 2.869 h each |
| **total spend** | **$9.12** compute + EBS (cap **$12**); ≈ $1.52 per six-member score |
| termination | **6 of 6 terminated**, verified by `describe-instances` (`work/cost318.txt`); the monitor terminated each worker the moment its `agg_results.csv` landed, before the ~4 h `prep-real-bundle` stage |
| local | peak RSS **4.18 GiB** (`mechanism318.py`; the profile scripts 1.7 GiB); **no local `cell-eval2`**, ever; `uv run python` throughout |
| writes | only inside `experiments/318-per-target-magnitude/` |
| submission | **none.** No `vcc` command run, nothing under `submissions/` written, `scripts/make_submission_pm_cd4.py` and `scripts/def_pipeline.py` untouched (the two 12.6 GB builds and `submissions/249-cd4svd5-f034.h5ad` were read only) |

One pre-existing `t4g.small` (`i-0a4311287d708e2ac`, running since 2026-09-04T19:00Z) is **not
ours** and was left alone, exactly as 317 recorded.

**Artifacts.** `profile318.py` (the delivered profile), `mechanism318.py` (the operator-stage
attribution), `build318.py` (the arms; 317's gates 1–5 kept, 6 and 7 added), `read318.py` (the
pricing and separability), `launch318.sh`, `monitor318.sh`, `wait318.sh`,
`aws/prologue-*.sh`; every number above is in `work/*.json` / `work/*.log` / `work/*.txt` and
the six `agg_results.csv` / `results.csv` in `cloud-output/`.

---

## 6. Verdict

**318 (per-target magnitude redistribution): NO-GO — and NOT SEPARABLE FROM ITS SHUFFLED
CONTROL.** Isolated on the real six-member scorer, holding direction and the whole-matrix norm
exactly fixed, the replayed profile buys **+0.000093 overall (1.0 % of the +0.0089 gate)**,
which is **z = −0.06 against three entropy-matched shuffles of its own multipliers** (best
shuffle +0.000543, itself only 6 % of the gate); **six times the measured redistribution is
−0.001101**, i.e. the axis is negative and gets more negative with dose; and the uniform
component that came with it costs **−0.001194** on its own. The +0.00075 the shipped arm bought
over the parity arm is not reachable this way, because **half of it is `pds`, which per-target
magnitude cannot move** (6× the dose moves `pds` by −0.00042, inside 317's re-seed floor and
opposite in sign), and because **the parity arm is not magnitude-neutral either** — its own
delivered redistribution (sd(log) 0.0485) is *larger* than the difference the mandate calls
"the distortion" (0.0329), since the norm match is enforced one stage before an operator whose
per-row gain spans 0.591–1.286.

**Spend: $9.12 of a $12 cap. Nothing competition-facing was selected, tuned or built.**


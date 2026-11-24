# 029/030 — Magnitude sharpening: a large `pds` lever that `nmae` may pay for

**Date:** 2026-08-21 · LODO, 259 VCC-panel perturbations, on top of whitening (γ = 0.65,
equal donors). `pds` is exact; `nmae` is the experiment 018 proxy.

## What was being tested, and what actually happened

The hypothesis was that `nmae` rewards restraint — a gene predicted as zero contributes
exactly `|true|` (ratio 1), a gene predicted confidently wrong contributes more — so shrinking
low-confidence genes toward zero should buy `nmae` and cost `pds`.

**It did the opposite.** Weighting each predicted effect by the cross-donor confidence raised
to a power moved macro PDS 0.8017 → 0.8405 and made `nmae` *worse*.

## Experiment 030: it is not consensus, it is sharpening

The confidence signal is the geometric mean of the donors' magnitudes and the prediction is
roughly their arithmetic mean, so `P · c^λ` is close to `sign(P)·|P|^(1+λ)` — a power
transform of the prediction's own magnitudes, using no cross-donor information at all. Testing
the two separately:

| variant | q = 1.5 | q = 2 | q = 2.5 |
|---|---|---|---|
| **`self` — `sign(P)·|P|^q`, no cross-donor info** | **0.8293** | **0.8399** | **0.8415** |
| `geo` — weight by cross-donor geometric magnitude | 0.8236 | 0.8332 | 0.8387 |
| `mask` — `geo` plus a sign-agreement gate | 0.8037 | 0.8212 | 0.8329 |

`self` matches or beats `geo` everywhere, and the sign-agreement gate is worst. **The
mechanism is the shape of the magnitude distribution, not consensus.** `gene_norm` and
whitening between them leave the predicted magnitudes flatter than the truth's, and a
one-parameter exponent fixes it — with no extra information required.

`pds` gain at q = 2: **+0.038 macro PDS**, on top of whitening's +0.070.

## The cost, and the honest caveat

`nmae` at our actual operating magnitude (α = 0.92 of the true effect size, inverted from v3's
board `mse`):

| q | pds (scaled) | nmae @ α=0.92 (scaled) | sum |
|---|---|---|---|
| 1.0 | 0.667 | 0.036 | 0.703 |
| **1.5** | 0.728 | −0.014 | **0.714** |
| 2.0 | 0.751 | −0.076 | 0.675 |
| 2.5 | 0.755 | −0.158 | 0.597 |

At a fixed emitted magnitude, `nmae` gives back roughly what `pds` gains, and the best net is
a marginal +0.011 scaled (≈ +0.002 overall) at q = 1.5.

**But `pds` is exactly α-invariant and `nmae` is not**, so the two can in principle be
decoupled: sharpen for `pds`, then re-tune α for `nmae`. At each variant's own best α,
q = 2 beats q = 1 by +0.058 scaled ≈ +0.010 overall.

**And this is where the proxy runs out.** It puts the `nmae`-optimal α at ~0.75; the
count-level run (exp 023, the real metric) measures the optimum at **α = 1.08**, with 0.7
*worse*. They disagree, so the decoupling argument rests on a number the proxy gets wrong.

**Not shipped, and not resolvable with the tools that produced it.** It needs the somatic
count-level harness (`data/harness/orion-hct116.h5ad`, 132 real VCC-panel targets in HCT116),
which is exactly why that was built tonight. The `pds` gain is exact and large; the `nmae` cost
is estimated by an instrument known to be miscalibrated on the axis that decides it.

---

# 030b/c — what the sharpening is *not*, and the principle that survives

Two follow-ups, both of which killed an explanation I had written down as fact.

## It is not distribution matching

The claim was that `gene_norm` and whitening leave the predicted magnitude distribution
flatter than the truth's. That is measurable, and it is barely true:

| stage | sd(log \|effect\|) | top-decile mass share |
|---|---|---|
| donor mean, pre-whitening | 1.1143 | 0.2642 |
| whitened, gamma = 0.85 | 1.1377 | 0.2633 |
| **truth (held-out line)** | **1.1955** | **0.2720** |

A power transform multiplies `sd(log|e|)` by exactly `q`, so the **matching** exponent is
`1.1955 / 1.1377 = 1.05`. The **PDS-optimal** exponent is 2.0–2.5. At q = 2 the prediction is
far *spikier* than the truth and PDS improves anyway. Whatever sharpening is doing, it is not
making the prediction resemble the target's shape.

## It is not selection either

If the transform works by suppressing the unreliable tail of small per-gene entries, its hard
counterpart should do at least as well — and `Jessica A` sits on the board with
`K562 top128 transfer v1` and the only raw `mse` below 1, so the hard version is a real
strategy someone is running.

| variant | macro PDS |
|---|---|
| full (no transform) | 0.8011 |
| top-128 | 0.8155 |
| top-1,000 | 0.8021 |
| top-2,000 | 0.7912 |
| top-4,000 | 0.7990 |
| **q = 2.5 (smooth)** | **0.8409** |
| top-1,000 then q = 1.5 | 0.8225 |

Hard truncation at every k underperforms the smooth transform, and combining them is worse
than the smooth transform alone.

## What the mechanism actually is

Multiplicative shrinkage whose strength *decreases* with the observed magnitude:
`|e| -> |e| * |e|^(q-1)`. If per-gene noise is roughly additive and comparable across genes,
then small entries have the worst signal-to-noise, and shrinking them harder is the standard
answer. The transform keeps every gene's contribution — weak but real — rather than deleting
it, which is why it beats a cutoff.

## The principle, which now has two independent instances

**Smooth reweighting beats hard truncation, in every basis this project has tried.**

- Across perturbations: SVD truncation loses monotonically (exp 004); smooth spectral
  flattening gains +0.070 (exp 017).
- Within a perturbation: top-k selection gains little or loses (above); smooth power
  sharpening gains +0.038.

It also dissolves an apparent contradiction. Experiment 004 concluded "the low-variance tail
carries signal, do not truncate", and this concludes "the small-magnitude tail is noise,
suppress it". Both are right: a low-variance *SVD direction* is a coherent pattern shared
across perturbations, while a small *per-gene entry* is one number in one prediction. Different
tails, different objects.

Whether the `pds` gain survives its `nmae` cost is a count-level question, and the q = 1.5 and
q = 2.0 arms of experiment 031 are running on the somatic harness to answer it.

---

# 030d/e — the two transforms compose, and sharpening replicates

## Joint (gamma, q), 3-line LODO on the VCC panel

| gamma \ q | 1.0 | 1.5 | 2.0 | 2.5 |
|---|---|---|---|---|
| 0.00 | 0.7210 | 0.7578 | 0.7826 | 0.7932 |
| 0.65 | 0.8017 | 0.8293 | 0.8399 | **0.8415** |
| 0.85 | 0.8011 | 0.8287 | 0.8385 | 0.8409 |
| 1.05 | 0.7954 | 0.8234 | 0.8341 | 0.8363 |

Whitening alone +0.080, sharpening alone +0.072, **both +0.120** — ~79% additive. They are
reweighting different things: directions *across* perturbations versus magnitudes *within*
one. The optimal gamma drifting 0.85 → 0.65 as q rises is the small overlap showing.

## Replication on the independent 4-line LODO

Three predictors per fold, and only 20 of the 226 shared perturbations on the VCC panel, so
this is a different gene population as well as a different fold structure.

| config | k562 | hct116 | hek293t | h1 | macro |
|---|---|---|---|---|---|
| gamma 0, q 1 | 0.8096 | 0.7566 | 0.7351 | 0.8484 | 0.7874 |
| gamma 0.65, q 1 | 0.8536 | 0.8079 | 0.8063 | 0.8781 | 0.8365 |
| **gamma 0.65, q 2** | 0.8686 | 0.8404 | 0.8570 | 0.8916 | **0.8644** |

**+0.028 from sharpening on top of whitening**, consistent across all four folds, with the
peak at q = 2 here against 2.5 on the 3-line harness — a broad plateau either way.

## Order of operations: sharpen AFTER whitening

The two do not commute and there is no first-principles argument for either order.

| | q = 1.5 | q = 2 | q = 2.5 |
|---|---|---|---|
| whiten, then sharpen | **0.8565** | **0.8644** | **0.8641** |
| sharpen, then whiten | 0.8554 | 0.8586 | 0.8508 |

"After" wins at every exponent and the gap widens with q. Sharpening first distorts the
magnitudes the covariance is estimated against, so the operator is fitted to a space the
prediction no longer lives in.

---

# 030f — the third basis works, and it is redundant with the second

There was one obvious basis left after whitening (directions across perturbations) and row
sharpening (magnitudes within a perturbation): **a gene's profile across the panel**. The same
argument applies — for gene g the 300 values say how much it responds to each perturbation, and
the large ones are the better-estimated ones. Doing both axes at once is Sinkhorn-style
bi-normalisation, which is standard in retrieval.

| row q \ col r | 1.0 | 1.25 | 1.5 | 2.0 |
|---|---|---|---|---|
| 1.0 | 0.8017 | 0.8173 | 0.8283 | **0.8395** |
| 1.5 | 0.8293 | 0.8380 | 0.8419 | **0.8433** |
| 2.0 | 0.8399 | 0.8429 | 0.8432 | 0.8400 |

Sinkhorn, alternating 1.35 on both axes: 0.8430 after two rounds, **0.8150 after four**.

**Column sharpening alone reaches 0.8395; row sharpening alone reaches 0.8399. Doing both
gets 0.8433.** Either axis gets you to essentially the same place, and combining them adds
+0.0034 — inside the zone where this harness stops resolving.

That is a boundary worth having. Whitening and sharpening were ~79% additive because they
reweight genuinely different structure. Row and column sharpening are ~10% additive because
they are two routes to removing *the same* redundancy: the unreliable small-magnitude entries,
approached once by row and once by column.

**The magnitude-reweighting family is exhausted at ~0.843.** There is no fourth axis to try
here, and iterating the two we have makes it worse. Any further gain has to come from
somewhere structurally different — which, on current evidence, means more donor lines rather
than more transforms.

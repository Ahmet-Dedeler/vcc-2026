# 095 — the DE questions, re-asked on a harness that is not depth-confounded

**Date:** 2026-08-27 · harness `data/harness/orion-hct116-thinned.h5ad`
(`--sampling greedy --thin-controls`, depth ratio `1.0004`, zero-fraction gap `−0.0015`) ·
prediction rebuilt by `experiments/031-somatic/run.py` at the production grid
`g=0.65 q=1.0 a=1.08 spread=0.25`, union axis, LODO from `k562_gwps` + `orion_hek293t` ·
scored on the anchors fitted from 780 live board entries.

## What this file concludes, in six lines

1. **Every emission counterfactual in this project was wrong.** Reordering a DE table by
   fiat said `+0.00583`; driving the real emitter says `−0.00021`. `realize.py` (§8) is the
   gate that catches this, and it has a `5e-4` noise floor.
2. **`spread` cannot steer the confidence ordering.** Nine policies, §8b — everything that
   preserves `n_pred` lands at or below zero. It is a strong lever on *membership* and a weak
   one on *order*.
3. **Gamma-Poisson over-dispersion can.** §10: φ=4 at production spread takes reach
   `0.084 → 0.161` with `nmae` invariant. That is the only emission change measured here
   that moves the ordering.
4. **Its optimum is set by the scorer, not by realism.** §13–14: a per-gene φ matched to each
   gene's own control dispersion *loses* to a single global φ.
5. **No harness we hold can settle a call-set question** (§5): `n_real` is 68 here, 322 on the
   broken harness, 222 on h1-test — all pin `fid`'s coverage term at 1, so all three say
   "truncate", and the board has refuted that twice.
6. **The remaining `reach` headroom is a model problem.** §16: over-dispersion captures 14%
   of the ordering ceiling; the other `+0.0875` needs better per-gene *direction*.

## Why

Everything experiments 092, 093 and 094 measured about direction, `fid`, `reach`, `jac` or
`n_real` was read off `data/harness/orion-hct116.h5ad`, whose control cells are 38% deeper
than its perturbation cells. `docs/HARNESS-DEPTH-CONFOUND.md` has the diagnosis; this
re-runs the measurements.

## Results

| quantity | broken harness | **repaired** | h1-test | board (v15) |
| --- | ---: | ---: | ---: | ---: |
| reference DE genes negative | 0.908 | **0.609** | 0.500 | — |
| `n_real` median | 322 | **68** | 222 | — |
| our `n_pred` median | ~5,340 | 4,020 | — | ~6,000 |
| sign accuracy on the reference's DE set | 0.542 | **0.574** | — | — |
| raw `fid` | 0.531 | 0.541 | — | **0.4967** |
| raw `jac` | 0.034 | 0.012 | — | 0.0321 |

### 1. `|lfc|` ordering beats the production ordering, and the result replicates

`reach` walks the reference's DE genes in **our** confidence order. The pool and each gene's
correctness are fixed by the model; only the order is free.

| ordering | raw reach | scaled | Δoverall vs production |
| --- | ---: | ---: | ---: |
| random permutation | 0.05151 | −0.0311 | −0.00480 |
| **production** (`sig` first, then `p_adj`) | 0.07702 | −0.0023 | — |
| **by `\|predicted lfc\|`** | **0.10082** | +0.0246 | **+0.00448** |
| oracle (correct signs first) | 0.62498 | +0.6168 | +0.10318 |

The broken harness measured the same comparison at **+0.00497**. Agreement to `0.0005`
across a harness repair that moved the reference's direction balance from 0.908 to 0.609
means **this one is not an artefact of the confound**. Production ordering also beats random
by `+0.0048`, so our emitted p-value order does carry real information — it is simply worse
than the model's own magnitude ordering.

### 2. What this does and does not rescue from experiment 094

094 found `Spearman(our p_adj, sign correct) = +0.50` on the broken harness and built arm
`a1-order-k2000` around promoting the top 2,000 genes by `|effect|`. The `+0.50` was
retracted — it was our up-first ordering multiplied by that screen's 90.8%-down reference.

**But the part of a1's rationale that survives is the part measured here**: an `|lfc|`-keyed
order is worth about `+0.0045` over the emitted-p order, on both harnesses. a1's mechanism
was a two-level approximation of exactly that. The reason given for withdrawing it — that a
directionally mixed head cannot clear the 0.9 purity bar — is **not** the operative
objection, because purity is computed on the **lfc sign**, which `spread` cannot move.

### 3. `--keep-top 2000` is the optimum, and it beats the full ordering

`--keep-spread` below `--tail-spread` lowers the top-K genes' p-values, moving them toward
the front of the order the scorer walks. Modelling that promotion directly (`pick_k.py`):

| policy | raw reach | scaled | Δoverall |
| --- | ---: | ---: | ---: |
| production (K=0) | 0.07702 | −0.0023 | — |
| promote top-250 | 0.08720 | +0.0092 | +0.00192 |
| promote top-500 | 0.09027 | +0.0127 | +0.00249 |
| promote top-1000 | 0.09774 | +0.0211 | +0.00390 |
| **promote top-2000** | **0.10799** | **+0.0327** | **+0.00583** |
| promote top-3000 | 0.10717 | +0.0318 | +0.00568 |
| promote top-4000 | 0.09223 | +0.0149 | +0.00287 |
| promote top-6000 | 0.07980 | +0.0009 | +0.00052 |
| *full `\|lfc\|` order* | 0.10082 | +0.0246 | +0.00448 |

A broad peak at 2,000–3,000, and **the two-level promotion beats the full `|lfc|` ordering**
— because it keeps the production p-order *within* the promoted block, which is better than
`|lfc|` there. So the right policy is not "order everything by magnitude", it is "promote
the block the model is confident about, and leave the emitted order alone inside it".

The K=2000 that arm `a1-order-k2000` was built with on 2026-08-27 was chosen before any of
this was measured, on the broken harness, and lands on the peak.

**The counterfactual models promotion only, and the real mechanism does two things.**
`--keep-spread 0.0` lowers the top-K's p-values, which promotes the ones already called —
but only **37.9%** of the top-2000 by `|predicted lfc|` are already significant (min 22.4%
over perturbations), so up to ~1,240 per perturbation can become **new** calls, a **31%**
growth in `n_pred`. The table above credits none of that. Whether the extra calls help or
hurt is genuinely open: they are the highest-`|effect|` genes we hold, so their directional
precision should be above our 0.497 average — but `fid`'s coverage term is already
saturated on the board, so any dilution costs. **This is why the arm is a board measurement
and not a shipped default.**

### 4. Three cautions on the numbers above

* **`n_real` collapses to 68** on the repaired harness, against 322 broken and 222 on
  h1-test. Thinning the controls costs detection power, so `k*` is a handful of genes and
  the reach numbers rest on a thin base. The *ranking* of orderings is stable across both
  harnesses; the magnitudes are not.
* **`fid` still does not transfer.** 0.541 here against **0.4967** on the board. With
  `n_real` at 68 and `n_pred` at 4,020 the coverage term is pinned at 1 and `fid` is pure
  precision, which is not the regime the challenge is in. Do not tune a call-set size here.
* **The residual direction skew is 0.609, not 0.500.** Thinning fixes depth, not the
  unmatched batch effects between a control pool from batches 1–20 and perturbation cells
  from all 109.

### 5. The `fid` curve is still measured in the wrong regime — and the board disagrees with it

`fid_curve.py`, call set = our top-*m* genes by predicted adjusted p-value:

| m | raw `fid` | raw `jac` | Δ(fid+jac) overall |
| ---: | ---: | ---: | ---: |
| 250 | 0.6101 | 0.0221 | **+0.04341** |
| 1,000 | 0.5895 | 0.0175 | +0.02978 |
| 2,000 | 0.5694 | 0.0146 | +0.01718 |
| **4,020 (production)** | 0.5408 | 0.0120 | — |
| 6,000 | 0.5284 | 0.0103 | −0.00775 |
| 10,000 | 0.5168 | 0.0084 | −0.01511 |

Monotone: shrink the call set, `fid` and `jac` both rise. That is the same shape experiment
092 got on the **broken** harness, and **the board refuted it twice** — v14 and v15 both
truncated and both lost.

The reason is regime, not confounding. Here `n_real` is 68, so `min(1, m/n_real)` is pinned
at 1 for every m in the table and `fid` degenerates to pure precision, which of course
improves as you keep only your best calls. On the challenge `n_real` is large enough that
the coverage term bites. **No harness we hold — broken at 322, repaired at 68, h1-test at
222 — is in the challenge's regime**, so none of them can settle a call-set question.

The only measurements in the right regime are our own three board points, and they run the
other way:

| | `n_pred` | raw `fid` |
| --- | ---: | ---: |
| v14 | ~2,317 | 0.49126 |
| v15 | ~2,370 | 0.49654 |
| **v11** | **~5,340** | **0.50178** |

Monotone **increasing** in `n_pred` — but reading that as "more calls is better" is reading
a gradient past the point where its mechanism stops. Inverting it properly, using the
repaired harness's `fid(m)` as a precision curve (valid there because coverage is pinned at
1) shifted to the board's level:

| entry | `n_pred` | raw `fid` | implied precision | **implied coverage** |
| --- | ---: | ---: | ---: | ---: |
| v14 | 2,317 | 0.49126 | 0.5327 | 0.9221 |
| v15 | 2,370 | 0.49654 | 0.5319 | 0.9335 |
| **v11** | **5,340** | **0.50178** | 0.5018 | **1.0000** |

**Coverage is already saturated at v11.** The v14 → v11 rise was coverage *recovery*, and it
is spent. Beyond v11, `fid` is pure precision and declining:

| projected `n_pred` | precision = `fid` | Δoverall vs v11 |
| ---: | ---: | ---: |
| 6,000 | 0.4989 | −0.00164 |
| **7,860 (a3's +31%)** | **0.4930** | **−0.00494** |
| 9,000 | 0.4902 | −0.00649 |

So a3 is **not** a free ordering test: its promotion is worth about `+0.0058` and its
call-set growth about `−0.0049`, which nets to roughly zero. The fix is to hold `n_pred`
fixed — raise the tail enough to offset the genes the promotion adds — so the ordering
channel is measured on its own.

### 6. The best buildable ordering key, and why we cannot spend it yet

`reach`'s oracle is worth `+0.103`. Which *computable* key gets closest? Tested as full
orderings of the reference's DE pool, with the two donors the LODO prediction actually used
(both carry all 132 targets):

| key | raw reach | Δoverall |
| --- | ---: | ---: |
| production (emitted p-order) | 0.07702 | — |
| `\|predicted lfc\|` | 0.10082 | +0.00448 |
| geometric magnitude `√(\|A\|·\|B\|)` alone | 0.09172 | +0.00277 |
| `min(\|A\|,\|B\|)` | 0.08832 | +0.00213 |
| cross-donor **sign agreement** first, then `\|lfc\|` | 0.11050 | +0.00630 |
| **`√(\|A\|·\|B\|)` × sign agreement** | **0.11302** | **+0.00678** |
| oracle | 0.62498 | +0.10318 |

**Experiment 013's key is confirmed, five days later, against the real metric** — it beats
`|lfc|` by 51% and is the best buildable ordering we have found.

**But it does not survive the mechanism we can currently ship.** `--keep-top K
--keep-spread 0` is a *two-level block promotion*, not a full ordering, and it keeps the
emitted p-order **inside** the promoted block:

| K | promote by `\|effect\|` | promote by agreement |
| ---: | ---: | ---: |
| 500 | +0.00249 | +0.00358 |
| 1,000 | +0.00390 | +0.00355 |
| **2,000** | **+0.00583** | +0.00341 |
| 3,000 | +0.00568 | +0.00114 |

The ranking inverts. `|effect|`'s top-2,000 is a set whose *internal* emitted order is
already informative; agreement selects a different, smaller-effect set whose internal order
is worse, and the block mechanism cannot fix that because it does not reorder within the
block.

**So the agreement key needs a graded, monotone spread map** — `spread_g = f(rank_g)` across
the whole axis — rather than the two levels `--keep-spread` provides. That is the next
concrete build, and it is worth roughly `+0.001` over today's arm on this harness, with more
available because a graded map can approach the full ordering rather than approximating it.

Today's arm therefore ships `|effect|` at K=2000, which is the peak of what is buildable.

### 7. How harness reach deltas scale to the board — one data point, use with care

v14 is the only artefact where the same policy was measured on both surfaces:

| | harness raw reach | board raw reach |
| --- | ---: | ---: |
| production / v11 | 0.014652 | 0.1711 |
| `--tail-spread 1.0 --keep-top 850` / v14 | 0.017459 | 0.1833 |
| **Δ** | **+0.0028** | **+0.0122** |

The board delta is **4.4×** the harness delta. If that amplification held, today's modelled
`+0.031` raw would imply a board move well past the observed board maximum, which it
obviously cannot — so the relationship is not a constant factor, and the honest range for
arm a5 is roughly **+0.002 to +0.013 overall, most likely ~+0.005**.

Recording it because the direction matters: the harness *understates* reach movement, which
is the opposite of how this repo has historically corrected proxy numbers (`docs/V9-RESULT.md`
has h1-test *overstating* the DE metrics by 4–5×). Do not apply one correction factor to
both surfaces.

### 8. REALISED reach kills the counterfactual — the spread channel cannot steer the order

Everything above is a **counterfactual**: it reorders a fixed DE table by fiat. The
mechanisms we can ship set a per-gene `spread`, and the realised p-value order is whatever
the Wilcoxon then produces. `realize.py` closes that gap — it drives the **unmodified
production `emit_counts`** with each policy's spread vector, runs the exact DE path, and
scores what the scorer would actually see. 60 perturbations, paired RNG streams, and a
**same-policy second seed as a noise control**:

| policy | realised reach | `n_pred` med | Δoverall |
| --- | ---: | ---: | ---: |
| production | 0.07475 | 3,990 | — |
| **production, different seed** | 0.07449 | 3,982 | **−0.00005** ← noise floor |
| block `\|effect\|` K=2000 (arm a5) | 0.07362 | 4,245 | **−0.00021** |
| graded `\|effect\|` (0.00 → 0.50) | 0.07063 | 3,918 | −0.00077 |
| graded **agreement** (arm a6) | 0.06449 | 3,911 | −0.00193 |

**Emission-seed noise is 5e-5, so these are precise.** And the counterfactual's `+0.00583`
for the K=2000 block promotion realises as **−0.00021** — not merely smaller, but the wrong
side of zero.

**Why.** The counterfactual treated "promoted" as a hard primary sort key. The real
mechanism only *lowers p-values*, and predicted p-values span hundreds of orders of
magnitude while `spread` moves them by a few. The promoted genes were mostly significant and
already near the front; the ones that were not do not travel far enough to matter. **The
spread channel changes which genes are called far more readily than it changes their
order.**

That is the same trap experiment 093 fell into, one level up: 093 checked that its
counterfactual was realisable and still mis-predicted the board, and here a counterfactual
that looked robust across two harnesses does not survive contact with the emitter at all.

### 8a. DISCLOSURE: the block-promotion arms had an oracle leak, and it strengthens the result

`realize.py`'s block policies (`k500` … `k6000`, and `graded_effect`) selected the promoted
set with

```python
rl  = log2(expm1(Br[:, keep]).mean(0) / ref_lin)   # Br = the REAL perturbation cells
eff[keep] = np.abs(rl)                              # "|lfc| stands in for |effect|"
```

`rl` is the **reference's** log fold change, not our prediction's. So those arms ranked genes
by an **oracle magnitude we could never compute at build time**. That is a real error in the
measurement and it is disclosed here rather than quietly fixed.

**It biases those arms optimistically, and they still failed.** Every block policy that
preserved `n_pred` landed at or below zero *while being handed the oracle's own magnitude
ordering*. So the conclusion — `spread` cannot steer the confidence order — is **stronger**
than §8/§8b state, not weaker: even a key we cannot build does not rescue the mechanism.

Two things are unaffected:

* **the over-dispersion arms**, which take a scalar `φ` and rank nothing (a7, a8, and every
  φ row in §10, §12 and §15 are clean);
* **the `graded_agree` arm**, which used cross-donor deltas — legitimately buildable.

And one comparison must be read with care: §8's headline contrast (`+0.00583` counterfactual
against `−0.00021` realised) pits a counterfactual ranked on **our predicted** `|lfc|`
against a realisation ranked on the **oracle's**. The realisation had the easier key and
still lost, so the direction of the contrast holds — but they are not the same policy, and
anyone re-running this should fix the key before quoting the gap as a like-for-like number.

### 8b. The K sweep: nothing escapes, and the exceptions prove the rule

| policy | realised reach | `n_pred` med | Δoverall |
| --- | ---: | ---: | ---: |
| production | 0.07475 | 3,990 | — |
| production, different seed | 0.07449 | 3,982 | −0.00005 |
| k500 | 0.07420 | 4,023 | −0.00010 |
| k1000 | 0.07226 | 4,074 | −0.00047 |
| **k2000 (arm a5)** | 0.07362 | 4,245 | **−0.00021** |
| k4000 | 0.07122 | 4,562 | −0.00067 |
| k6000 | 0.07139 | 4,922 | −0.00063 |
| k2000, tail 0.40 | 0.07519 | **3,378** | +0.00008 |
| k2000, tail 0.60 | 0.07698 | **2,499** | +0.00042 |
| graded 0.0 → 1.0 | 0.07602 | **2,594** | +0.00024 |

Every policy that leaves `n_pred` near production is at or below zero. **The only three that
are positive all shrink the call set by 15–40%** — and that trade is not hypothetical, the
board measured it twice: v14 at `n_pred` 1,563 and v15 at 1,911 lost `0.0029` and `0.0007`
against v11, because `fid`'s coverage term is saturated at v11 and gives back far more than
`reach` gains. `k2000_tail060`'s `+0.00042` of reach against a 40% call-set cut is the same
losing trade at a smaller scale.

**So the emission ordering channel is closed.** Not "unpromising" — measured, with a 5e-5
noise floor, across nine policies spanning both mechanisms and the full useful range of K.

### 9. What this costs Addendum 7

`STRATEGY.md` Addendum 7 says the confidence ordering is "a free parameter, orthogonal to
the predicted expression profile". The first half stands — the ordering is what `reach`
reads, and moving it cannot touch `pds`, `mse` or `nmae`. **The second half needs
qualifying: free is not the same as steerable.** With the emitter we have, `spread` is a
weak lever on order and a strong lever on membership, so any policy that reorders enough to
matter also shrinks the call set enough to lose `fid` — which is precisely how v14 and v15
lost, and why every arm in the table above lands at or below zero.

Steering the order therefore needs something `spread` cannot do: a change to the emitted
*distribution shape* (the gamma-Poisson extension exp093 §7 probed and never validated), or
a change to the profile itself, which is no longer free because it moves `pds` and `nmae`.

### 10. Over-dispersion IS the lever, and there is an operating point that keeps `n_pred`

`spread` caps at Poisson. Gamma-Poisson scatter (`Var = mu + phi*mu^2`, same exact-total
repair) does not, and it is the only emission change measured here that moves the order.
Realised through `emit_counts(..., overdispersion=phi)`, 40 perturbations, `production_seed2`
as the noise control (`−0.00046`):

| policy | `n_pred` | reach | `fid` | `jac` | `nmae` | Δreach | Δfid | Δjac | **ΔTOTAL** |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| production | 3,934 | 0.0838 | 0.5465 | 0.0133 | 0.9295 | — | — | — | — |
| spread 0.05, φ=20 | 5,112 | 0.0853 | 0.5351 | 0.0115 | 0.9295 | +0.0003 | −0.0064 | −0.0008 | −0.0069 |
| spread 0.10, φ=8 | 4,644 | 0.1111 | 0.5408 | 0.0117 | 0.9294 | +0.0051 | −0.0032 | −0.0007 | +0.0012 |
| spread 0.25, φ=0.5 | 3,740 | 0.0891 | 0.5511 | 0.0131 | 0.9295 | +0.0010 | +0.0025 | −0.0001 | +0.0035 |
| **spread 0.25, φ=2** | **3,218** | **0.1378** | 0.5644 | 0.0137 | 0.9293 | **+0.0102** | +0.0100 | +0.0002 | **+0.0204** |
| spread 0.25, φ=8 | 2,513 | 0.1086 | 0.5787 | 0.0196 | 0.9289 | +0.0047 | +0.0180 | +0.0029 | +0.0256 |
| spread 0.50, φ=0.5 | 1,692 | 0.1425 | 0.5917 | 0.0182 | 0.9293 | +0.0111 | +0.0253 | +0.0022 | +0.0386 |
| spread 1.0, φ=0.30 | 692 | 0.1961 | 0.5762 | 0.0288 | 0.9290 | +0.0212 | +0.0167 | +0.0070 | +0.0448 |

**`nmae` is invariant to four decimals across every row**, which is the exact-total repair
holding exactly as it must.

**Read the `fid` column with the coverage caveat.** `n_real` is 68 on this harness, so
`min(1, n_pred/n_real)` is pinned at 1 and every call-set cut looks free. On the challenge
`n_real` is far larger and coverage is saturated only at v11's ~5,961, so the big-φ rows —
`spread 1.0, φ=0.30` at `n_pred` 692, i.e. ~1,050 on the challenge — would take the same
`fid` beating v14 took. **Their `+0.045` is not real.**

**`spread 0.25, φ=2` is the operating point** precisely because it is the one that buys
reach without paying that: `+64%` reach (`0.0838 → 0.1378`) for an **18%** call-set cut, to
`n_pred` ~4,880 on the challenge, where coverage is still ~0.99. Board-anchored, that is
roughly `+0.010` reach and `+0.006` fid, so **≈ +0.016 overall** — and unlike every
counterfactual in this file, it was measured through the emitter that will actually produce
it.

### 10b. Credit: the 2026-08-26 session had this idea and did not execute it

Searching the prior session transcripts turns up, from `Memory Management experiment review`
(2026-08-26 01:06, immediately after v14/v15 were submitted):

> "…genes we do not [call] are emitted with realistic overdispersion, so they fall down the
> p-value order and out of our called set."

and, from `CCLE matching and per-perturbation confidence` (2026-08-22):

> "…the negative-binomial overdispersion used to turn a mean profile into 400 count vectors.
> **It is a real hyperparameter…**"

So the mechanism was sketched twice and never measured — the same fate as experiment 093
§7's `phi` probe, which ran on the confounded harness and was filed as "not shipped". What
was missing was not the idea but an instrument that could tell a real emission effect from a
counterfactual, which is what `realize.py` and its seed control are.

Recording it because it changes what to conclude: this is a **confirmed and quantified**
hypothesis, not a new one, and the fact that it survived three independent arrivals raises
rather than lowers confidence in it.

### 11. The scoring path is scorer-exact — verified against cell-eval2 itself

`cell-eval2` scored the same prediction on the same repaired harness independently of every
number above:

```
[1/1] g0.65_q1.0_a1.08_spread0.25: avg6=0.1368 pds=0.7953(0.653)
      nmae=0.9309(0.117) reach=0.0770(0.006) fid=0.5408 jac=0.0120
```

| metric | official `cell-eval2` | this experiment's own pipeline |
| --- | ---: | ---: |
| `reach` | 0.0770 | **0.07702** |
| `fid` | 0.5408 | **0.5408** |
| `jac` | 0.0120 | **0.0120** |

Agreement to the precision the scorer prints. So `de_tables.py`, `analyse.py` and — the one
that matters — `realize.py` are all reading the same quantity the competition reads, and
the realised-emission numbers in §8, §8b and §10 are scorer-exact rather than a
reimplementation that happens to correlate.

### 12. The φ optimum, and the seed variance over-dispersion introduces

Sweeping φ at production spread (0.25), realised, 40 perturbations:

| φ | `n_pred` | reach | Δreach | Δfid | Δjac | ΔTOTAL |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 (production) | 3,934 | 0.0838 | — | — | — | — |
| 0.5 | 3,740 | 0.0891 | +0.0010 | +0.0025 | −0.0001 | +0.0035 |
| 1.0 | 3,528 | 0.1113 | +0.0052 | +0.0050 | −0.0002 | +0.0101 |
| 1.5 | 3,357 | 0.1288 | +0.0085 | +0.0079 | −0.0001 | +0.0163 |
| 2.0 | 3,218 | 0.1378 | +0.0102 | +0.0100 | +0.0002 | +0.0204 |
| 3.0 | 3,018 | 0.1604 | +0.0144 | +0.0132 | +0.0007 | +0.0283 |
| **4.0** | **2,863** | **0.1675** | **+0.0158** | +0.0151 | +0.0011 | **+0.0320** |
| *4.0, second seed* | 2,880 | *0.1535* | *+0.0131* | +0.0156 | +0.0011 | *+0.0298* |
| 5.0 | 2,690 | 0.1418 | +0.0109 | +0.0165 | +0.0017 | +0.0291 |
| 6.0 | 2,636 | 0.1233 | +0.0075 | +0.0175 | +0.0022 | +0.0271 |
| 8.0 | 2,513 | 0.1086 | +0.0047 | +0.0180 | +0.0029 | +0.0256 |

**Reach peaks at φ=4 and falls away on both sides.** Past it the scatter is wide enough
that genes start being called significant in the *wrong* direction, which is the same
non-monotonicity exp093 §7 saw at φ=1.0 on the broken harness.

**Over-dispersion introduces real emission variance, and it must be controlled for.** Two
φ=4 draws give reach `0.1675` and `0.1535` — a spread of `0.014` raw, against production's
`0.0015` across its own two seeds, i.e. **an order of magnitude more seed noise**. Any future
φ result quoted from a single seed is unreliable at the `±0.002` overall level. The honest
φ=4 figure is the mean of the two, `0.1605`.

**And the `n_pred` transfer is now calibrated.** Probing the built φ=2 artefact against the
real challenge controls gives `n_pred` **4,523** (A 4,685 / B 4,878 / C 4,006) against v11's
5,961: a **24%** cut where this harness predicted **18%**. So the harness understates the
call-set cost of over-dispersion by about a third, and any φ chosen here should be checked
against that before a slot is spent.

### 13. Per-gene φ loses to a single global φ — and the first attempt had a bug worth recording

The principled form of over-dispersion is a **per-gene** φ fitted so each gene's emitted CV
matches its own control cells'. `emit_counts` now accepts a per-gene vector. It is worse:

| policy | `n_pred` | reach | Δreach | Δfid | Δjac | ΔTOTAL |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| production | 3,934 | 0.0838 | — | — | — | — |
| **global φ=4** | 2,863 | **0.1675** | **+0.0158** | +0.0151 | +0.0011 | **+0.0320** |
| per-gene, cap 8 | 2,456 | 0.1311 | +0.0089 | +0.0178 | +0.0026 | +0.0293 |
| per-gene, cap 20 | 1,898 | 0.1183 | +0.0065 | +0.0147 | +0.0042 | +0.0254 |
| per-gene, cap 50 | 1,760 | 0.1253 | +0.0078 | +0.0120 | +0.0048 | +0.0246 |

Global φ=4 wins on reach **and** preserves more of the call set. Matching the controls
gene by gene spends dispersion where the control happens to be noisy, which is not where the
reach ordering needs it — the control CV is dominated by expression level and dropout, not
by how confident our prediction is.

**The first attempt was a bug, and the shape of it generalises.** I wrote

    phi_g = CV_ctrl^2 - 1/mu_s

which is the formula for a gene whose *entire* total is scattered. In the quasi emitter only
`spread` of it is: the per-cell count is `(even part) + NB(mu_s, phi)` with
`mu_s = spread * T/n`, but the count's **mean** is the full `T/n`. So

    Var = mu_s + phi*mu_s^2,   CV^2 = Var / (T/n)^2,
    phi_g = (CV_ctrl^2 * (T/n)^2 - mu_s) / mu_s^2.

Dividing by the scattered mean instead of the full mean made the target look
already-exceeded and assigned **φ = 0 to 97.9% of genes** — the arm ran, produced `n_pred`
3,948 against production's 3,934, and read as a clean null. It would have been written up as
"per-gene φ does nothing" had the φ distribution not been checked directly.

**A null result whose mechanism was never verified is not a result.** The check that caught
it was two lines: print the assigned φ's percentiles and the fraction at zero.

### 14. Realism is not the objective

The control cells' own dispersion, over the 10,261 gated genes:

| | p10 | p25 | **p50** | p75 | p90 |
| --- | ---: | ---: | ---: | ---: | ---: |
| control CV | 0.914 | 1.191 | **1.591** | 2.210 | 2.966 |
| implied φ (`CV²`) | 0.836 | 1.418 | **2.532** | 4.884 | 8.798 |

The dispersion-matching φ at the median gene is **2.53**; the metric optimum measured today
is **4.0**. Close enough that "we are under-dispersed" is a fair description of the
direction, but **the optimum is set by the scorer, not by biology** — and §13 shows that
matching biology gene-by-gene actively loses. Over-dispersion is a lever that reorders
p-values, and it should be tuned as one.

### 15. The (spread, φ) frontier — the map for the next arm

`spread` sets how much of each gene's total is scattered at all; `φ` sets how heavy that
scatter is. They trade against each other, and the binding constraint is `n_pred`, because
`fid`'s coverage term is saturated at v11 and charges for any large call-set cut:

| policy | `n_pred` (Δ vs production) | reach *(2-seed mean where available)* | harness ΔTOTAL |
| --- | ---: | ---: | ---: |
| production | 3,934 | 0.0838 | — |
| spread 0.15, φ=6 | 4,008 *(+2%)* | 0.1194 | +0.0082 |
| **spread 0.20, φ=4** | **3,494 *(−11%)*** | 0.1451 | +0.0189 |
| spread 0.25, φ=3 | 3,014 *(−23%)* | 0.1552 | +0.0272 |
| **spread 0.25, φ=4 — arm a8** | **2,872 *(−27%)*** | **0.1605** | **+0.0309** |
| spread 0.30, φ=3 | 2,427 *(−38%)* | 0.1429 | +0.0328 |
| spread 0.35, φ=2 | 2,170 *(−45%)* | 0.1631 | +0.0394 |

Two things to read off it.

**The harness ΔTOTAL column is monotone in the call-set cut, and that is the coverage
artefact, not a real ranking.** `n_real` is 68 here so every cut looks free on `fid`; on the
challenge the cut is charged. Rank these by `reach` and `n_pred`, never by the harness total.

**`spread 0.20, φ=4` is the coverage-safe point** — it keeps 89% of the call set (which
projects to `n_pred` ~5,070 and coverage ~0.995 on the challenge) and still takes reach from
`0.084` to `0.145`. a8 sits mid-frontier at a 27% cut, which the probe put at `n_pred` 3,836
and coverage ~0.975 on the challenge.

If a8's board result shows the `fid` coverage term biting harder than projected, **`spread
0.20, φ=4` is the retreat**, and it costs only about a quarter of the reach gain.

### 16. How much of the ordering ceiling over-dispersion actually buys: 14%

The oracle ordering — correct signs first — is invariant to the emitter, because it depends
only on which genes we get the *sign* right on, and `lfc` is untouched by `spread` or `φ`.
So it is a fixed ceiling to measure against:

| ordering | raw reach | scaled | Δoverall vs production |
| --- | ---: | ---: | ---: |
| production | 0.08378 | +0.0054 | — |
| **φ=4 (arm a8)** | **0.16052** | +0.0921 | **+0.01445** |
| oracle | 0.62498 | +0.6168 | +0.10191 |

**Over-dispersion captures 14% of the ordering headroom**, leaving **`+0.0875` overall**
between φ=4 and a perfect confidence signal.

That is the honest framing of what shipped today: a large step on a channel that is mostly
still empty, and the remainder is **not** an emission problem. `lfc` signs are fixed by the
model, so closing the rest means predicting per-gene direction better — which is where
`nmae` (`+0.0140` to the leader) and `fid` (`+0.0064`) live too, and where every team ahead
of us is winning. The emitter can now *express* a good ordering; it cannot invent one.

### 17. WHY it works: the head gets genuinely purer, at every depth

`reach` can rise for two very different reasons — the ordering genuinely improves, or a
lucky run of correct signs lands at the front. Only the first transfers. `purity.py` measures
the profile directly: sign purity of the first *k* adjudicable genes of the reference's DE
pool, in our confidence order, 30 perturbations, paired seeds.

| depth into the reference's DE pool | production | **φ=4** | Δ |
| --- | ---: | ---: | ---: |
| first 1 | 0.8846 | **1.0000** | +0.1154 |
| first 2 | 0.8462 | 0.9615 | +0.1154 |
| first 5 | 0.7385 | 0.8385 | +0.1000 |
| first 10 | 0.6962 | 0.8269 | +0.1308 |
| first 20 | 0.6192 | **0.7750** | **+0.1558** |
| first 50 | 0.5663 | 0.7000 | +0.1337 |
| first 100 | 0.5391 | 0.6345 | +0.0955 |

**Purity improves at every depth, by +0.10 to +0.16.** This is not a reshuffle that happened
to help — the emitter is now ranking by something that genuinely predicts whether we have
the sign right, and it holds from the first gene to the hundredth.

That is the strongest argument that the gain transfers. `reach`'s 0.9 pass mark makes the
metric itself brittle — a single wrong call near the front truncates `k*` — but the
*underlying* quantity it samples has moved a long way, and it moved in the same direction
everywhere rather than at one depth.

The mechanism is the one §9 predicted: Poisson noise is too weak to reorder p-values that
span hundreds of orders of magnitude, so production's ordering is dominated by each gene's
control-distribution skewness. Widening the emitted distribution lets a gene's actual
predicted effect compete with that artefact.

### 18. The mechanism, confirmed directly

§9's claim was that `spread` cannot reorder because predicted p-values span hundreds of
orders of magnitude while Poisson noise moves them by a few, so the order is dominated by
each gene's control-distribution skewness rather than by our model. If over-dispersion works
for that reason, it should make the emitted p-order track our own `|predicted effect|` more
closely. It does, monotonically:

| φ | Spearman(p_adj, `\|predicted effect\|`) |
| ---: | ---: |
| 0 (production) | **0.060** |
| 2 | 0.170 |
| 4 | 0.223 |
| 8 | 0.275 |

Production's `0.060` matches the `+0.015`–`+0.07` measured independently on the real
challenge artefact and on the HCT116 tables: **our confidence order has been very nearly
independent of our own model's magnitudes**, which is the defect in one number.

**And the nuance that matters:** φ=8 has the *best* alignment (`0.275`) and *worse* reach
(`0.109` against φ=4's `0.168`). Alignment is necessary, not sufficient — past φ=4 the extra
scatter degrades the signal faster than the improved ordering recovers it, which is the same
non-monotonicity §12 saw in reach and exp093 §7 saw at φ=1.0. **The optimum is a balance
between letting the model's signal through and not drowning it**, and it sits at φ≈4 on this
harness.

### 19. Left running: the experiment 091 replication on the repaired harness

`run_prune_arm.sh` re-runs exp091's `essential` and `random` prune arms against
`orion-hct116-thinned.h5ad`. The essential arm was still in `cell-eval2` after 65 minutes
when this session ended; the random arm had not started.

**Its value is now limited, and worth stating so nobody over-reads it.** The board measured
the pruned basis exactly at `+0.0001257` (`experiments/091-basis-prune/BOARD-RESULT.md`), so
the replication can no longer tell us anything about the operator. What it *can* tell us is
how well the **repaired harness predicts the board** — the broken one said `+0.000876` and
over-stated by 7×. If the repaired harness lands near `+0.0002`, the repair improved its
predictiveness and that is a useful calibration for every future arm. If it lands at
`+0.0009` again, the depth fix did not help this class of measurement.

Results land in `out-essential/results.json` and `out-random/results.json`.

### 20. Selective over-dispersion fails, and the failure explains the mechanism

a8 lost `jac` to its call-set cut, so the obvious repair is to over-disperse only the genes
we do *not* trust: `φ = 0` on the top-K by predicted effect, `φ = 6` on the rest. Keys are
the **buildable** `|predicted effect|` (§8a), 30 perturbations:

| policy | `n_pred` | reach | Δreach | Δfid | Δjac | ΔTOTAL |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| production | 4,058 | 0.0944 | — | — | — | — |
| **uniform φ=4** | 2,942 | **0.1708** | **+0.0144** | +0.0151 | +0.0009 | **+0.0304** |
| keeptight K=1000 | 2,694 | 0.1278 | +0.0063 | +0.0179 | +0.0020 | +0.0261 |
| keeptight K=2000 | 2,708 | 0.1192 | +0.0047 | +0.0173 | +0.0018 | +0.0237 |
| keeptight K=3000 | 2,736 | 0.1224 | +0.0053 | +0.0180 | +0.0017 | +0.0250 |

**It costs a third of the reach gain at the same call-set size** — the `n_pred` values are
within 8% of uniform φ=4, so this is not a coverage effect. Exempting the trusted genes is
simply worse.

**And that explains what over-dispersion actually does.** It is *not* "push the bad genes
out of the call set" — if it were, protecting the good ones would help. It is "let every
gene's predicted effect compete with its control-distribution skewness", which §18 measured
directly as `Spearman(p_adj, |effect|)` rising `0.060 → 0.223`. Exempting the top-K leaves
exactly the genes that matter most still ordered by skewness, and those are the ones at the
head of the pool where `reach`'s 0.9 purity bar is decided.

So the `jac` cost is **not separable** from the reach gain by this route. Recovering it needs
a smaller call-set cut overall — the `spread 0.20, φ=4` point on the §15 frontier — not a
selective exemption.

## Files

| file | what |
| --- | --- |
| `de_tables.py` | scorer-exact real+predicted DE tables for any `--real` / `--pred` pair |
| `analyse.py` | reference balance, `n_real`, sign accuracy, and reach under four orderings |
| `grid.json` | the single production config |
| `run_prune_arm.sh` | queues the experiment 091 prune replication on this harness |

---

## 21. `spread` and `phi` are separable knobs — and a8's operating point is not on the frontier

Everything before this section tuned `spread` and `phi` as if they were one dial. They are not.
Sweeping each independently on the repaired harness (40 perturbations, `realize.py`, the
unmodified production emitter) separates their jobs cleanly:

- **`spread` sets `n_pred`.** It is the fraction of each gene's panel total that scatters at
  all, so it controls how many genes clear significance. At `phi=4`, `n_pred` runs
  5,414 → 2,863 as `spread` goes 0.025 → 0.25, monotonically.
- **`phi` sets how much ordering information that scatter carries**, at whatever `n_pred`
  `spread` has chosen.

### The frontier

`reach` at matched `n_pred`, all at 40 perturbations on `orion-hct116-thinned`:

| spread | phi | `n_pred` | `reach` |
| ---: | ---: | ---: | ---: |
| 0.25 | 0 *(production)* | 3,934 | 0.0838 |
| 0.05 | 4 | 5,238 | 0.0784 |
| 0.05 | 16 | 5,136 | 0.0884 |
| 0.10 | 4 | 4,800 | 0.0919 |
| 0.10 | 8 | 4,644 | 0.1111 |
| 0.10 | 16 | 4,403 | 0.1222 |
| 0.10 | 32 | 4,172 | 0.1209 |
| 0.15 | 4 | 4,161 | 0.1157 |
| **0.15** | **8** | **3,906** | **0.1339** |
| **0.15** | **16** | **3,553** | **0.1565** |
| 0.15 | 32 | 3,232 | 0.1374 |
| 0.175 | 16 | 3,120 | 0.1358 |
| 0.20 | 4 | 3,481 | 0.1473 |
| 0.25 | 4 *(**a8, shipped**)* | 2,863 | 0.1675 |

**Two results follow.**

**(a) The optimal `phi` rises as `spread` falls.** At `spread` 0.25 the peak is `phi=4`
(§ earlier: 0.1675, and `phi` 5 and 6 fall to 0.1418 and 0.1233). At `spread` 0.15 the peak
is `phi=16`. At `spread` 0.10 it is `phi` 16–32. Tuning `phi` at fixed `spread=0.25`, which
is all we did before a8, could never find this.

**(b) At matched `n_pred`, low-spread/high-`phi` dominates high-spread/low-`phi`.**
`s0.15_phi16` reaches 0.1565 at `n_pred` 3,553; `s0.20_phi4` reaches 0.1473 at 3,481.
**Better on both axes.** The arm built earlier today as tomorrow's slot 1 —
`a9 = spread 0.20, phi 4` — is therefore Pareto-dominated, and was rebuilt as
`a10 = spread 0.15, phi 16`.

**The sharpest single point:** at `spread 0.15, phi 8` the emitter holds `n_pred` at 3,906
against production's 3,934 — a 0.7% difference, i.e. **no call-set cut at all** — and lifts
`reach` from 0.0838 to 0.1339, **+60% for free on the coverage term**. Whatever the board
charges for cutting the call set, that arm does not pay it.

### Why this matters for the board specifically

a8's board result cost `-0.0034` on `jac` and delivered only `+0.0064` of a predicted
`+0.0390` on `fid`, both traced in `A8-BOARD-RESULT.md` to its `n_pred` cut
(5,961 → 3,836). Every arm in the table above reaches a given `reach` at a **higher**
`n_pred` than the `spread`-only route does, so each one pays less of that charge.

### Two honest limits on this section

**Seed noise at high `phi` is an order of magnitude worse than at production.** Two seeds of
`s0.15_phi16` give `reach` 0.1565 and 0.1688 — a spread of **0.0123**, against production's
~0.0015 noise floor. The `s0.15_phi16` vs `s0.20_phi4` margin is ~1.2x that spread, so
**result (b) is directionally supported but not tightly bounded**; the `n_pred` half of the
domination (3,553/3,544 vs 3,481) is solid because counts are near-noiseless.

**This rests on ONE harness, which is exactly the failure mode `A8-BOARD-RESULT.md` warned
about.** The h1-test cross-check could not be run: `realize.py` assumes the 18,533-gene
production panel and h1-test carries 18,080, and only 3 perturbations are available to it
there. h1-test is the harness that predicted *zero* for `phi=4` and was closest to the board
in absolute `reach`, so its silence here is a genuine gap, not a formality.
**Fixing `realize.py`'s gene-space alignment is the highest-value next piece of tooling.**
Until then a10 is a one-harness recommendation and should be read as such.

## 22. Both harnesses now run — and they overturn a8's operating point

Section 21 rested on one harness and named the h1-test cross-check as the gap. That gap is
now closed. Two things were wrong in `realize.py`, not one:

1. **Gene axis.** It assumed the 18,533-gene production panel; h1-test carries 18,080. Fixed
   by mapping the prediction onto the harness axis once at load and dropping from the gate
   any gene the prediction does not carry — emitting zeros for those would have fed a fake
   all-zero column into the rank test and scored it as a real gene.
2. **The prediction file, which was the bigger error.** It was pointed at the *challenge*
   prediction, whose 300 targets overlap h1-test's labels in only **3 perturbations** — so
   the earlier "h1-test only has 3 perturbations" reading was an artefact of my invocation,
   not a property of the harness. With the matched `h1-test-baseline-pred.h5ad` it runs the
   full 40.

### h1-test, 40 perturbations, matched prediction

| policy | `n_pred` | `reach` | `fid` | `jac` | ΔTOTAL |
| --- | ---: | ---: | ---: | ---: | ---: |
| production | 3,794 | 0.0794 | 0.5232 | 0.0682 | +0.00000 |
| production_seed2 *(noise floor)* | 3,800 | 0.0807 | 0.5236 | 0.0682 | +0.00047 |
| **`s0.15_phi16` (a10)** | 2,789 | 0.0872 | 0.5248 | **0.0682** | **+0.00236** |
| `s0.20_phi4` (a9) | 2,502 | 0.0849 | 0.5236 | 0.0637 | −0.00078 |
| `s0.25_phi4` (**a8, shipped**) | 2,054 | 0.0897 | 0.5141 | 0.0673 | **−0.00358** |

**a10 is the only arm that scores positive, at 5x the noise floor. a8's shipped point scores
negative.**

### The harnesses disagree about a8, and that disagreement is readable

| ranking | 1st | 2nd | 3rd |
| --- | --- | --- | --- |
| orion-hct116-thinned (`n_real` 92) | **a8** +0.0320 | a10 +0.0223 | a9 +0.0192 |
| h1-test (`n_real` 322) | **a10** +0.0024 | a9 −0.0008 | a8 −0.0036 |

They invert on a8 — and the direction of the inversion is exactly what `A8-BOARD-RESULT.md`
predicts. Orion's `n_real` is 92, so its `fid` and `jac` reward any call-set cut
mechanically; a8 cuts `n_pred` hardest (2,863 on orion, 2,054 on h1-test), so it collects the
most of that artefact and orion ranks it first. h1-test's `n_real` is 3.5x larger and the
artefact is correspondingly weaker. **The board's `n_real` is ~2,000 — an order of magnitude
beyond h1-test — so the trend across `n_real` 92 → 322 → 2,000 runs against a8, not for it.**

**The board already confirmed this trend once.** a8's two disappointments were `jac`
(−0.0034, wrong sign) and `fid` (+0.0064 of a predicted +0.0390) — both charges for cutting
`n_pred`. h1-test reproduces both: a8 is the only arm there that loses `fid` (0.5232 →
0.5141), and it loses `jac` too. **a10 holds `jac` at 0.0682, flat to four decimals against
production — the exact member a8 gave away on the board — while still taking `reach`
(+0.0078 raw) and `fid` (+0.0016).**

### Verdict

**a10 (`--tail-spread 0.15 --overdispersion 16.0`) replaces a9 as the next submission**, and
on this evidence it is a candidate to beat a8 rather than merely to probe around it. The
case does not rest on either harness alone: orion supplies the frontier geometry (§21, where
a10 dominates a9 on both axes), h1-test supplies the `n_real`-corrected ranking, and the
board's own a8 post-mortem supplies the extrapolation between them.

**What would falsify it:** if a10 scores below a8's `0.1543607`, then `reach` magnitude
dominates the `n_pred` charge on the challenge, orion's ranking was right after all, and the
correct move is back toward higher `spread` — test `s0.30_phi4` next, not lower `spread`.

### a10's board-side `n_pred`, measured

Probed reference-free against the real `data/challenge/context_*.h5ad` controls
(`probe_de.py`, 12 perturbations per context), so this is the challenge's own gene space and
depth, not a harness proxy:

| context | `n_pred` median | IQR | call set % UP |
| --- | ---: | --- | ---: |
| A | 4,882 | [4,387, 5,430] | 62.5% |
| B | 5,020 | [4,901, 5,285] | 62.6% |
| C | 4,370 | [4,047, 4,688] | 66.9% |

**Mean ~4,757, against a8's 3,836 and v11's 5,961.** a10 recovers ~55% of the call set a8
gave up, which is the whole point of the arm: the board charged a8 `-0.0034` on `jac` and
capped its `fid` gain at `+0.0064` for that cut.

The confidence ordering is also intact — `frac UP` runs 0.785 → 0.195 monotonically across
confidence deciles in every context, so the extra calls a10 makes are not noise being
readmitted at the bottom of the ranking.

**Regression note.** Section 21's orion numbers were measured before the `realize.py` gene-axis
patch. Re-run after it, they reproduce bit-exactly — production `n_pred` 3,934 / `reach`
0.08378, `s0.15_phi16` 3,553 / 0.15654 — so the patch is behaviour-preserving whenever the
harness and prediction panels already match, and only changes the mixed-panel path it was
written for. `realize_regress.log` holds the run.

# 094 — `reach` is a ranking metric, and the ordering channel is not steerable

> **READ THE ADDENDA FIRST.** This file's original headline — "our ranking is inverted",
> `Spearman(p_adj, correct) = +0.50` — is **RETRACTED**. It was our up-first ordering
> multiplied by one screen's depth-confounded direction distribution, and h1-test puts that
> distribution at 0.538 rather than 0.908. §3 and §4 must be read narrowly; the addendum
> beginning "the asymmetry does NOT transfer" says exactly which parts survive.
>
> The conclusion the file actually ends at, established in experiment 095 by driving the
> **real emitter** rather than reordering tables: the confidence ordering is *free* — moving
> it cannot touch `pds`, `mse` or `nmae` — but it is **not steerable with `spread`**. Nine
> policies, `5e-4` noise floor, everything that preserves `n_pred` at or below zero. What
> does work is **gamma-Poisson over-dispersion**, and that is what arms a7 and a8 ship.

**Date:** 2026-08-27 · tables `experiments/092-call-set-truncation/de_tables.npz`
(132 HCT116 perturbations × 10,258 gated genes, real and predicted p / p_adj / lfc) ·
scaling anchors fitted from 780 live leaderboard entries.

## Verdict *(as originally written — see the banner above; the +0.50 is retracted)*

Our confidence ordering — the thing `reach` walks — is **anti-correlated with sign
correctness**, `Spearman(our p_adj, correct) = +0.5009`. The scorer starts at the most
significant end, which is our *least* reliable end. This is an emission artefact, not a
modelling gap, and it is worth more than any other single lever on the board.

## 1. Where the score actually is

Fitting `scaled = (raw − b)/(r − b)` on 780 board entries reproduces four members to
`R² > 0.9999`, and every anchor lands inside the published spec range:

| member | `b` | `r` | span | raw Δ for **+0.01 overall** | our raw |
| --- | ---: | ---: | ---: | ---: | ---: |
| `pds` | 0.50020 | 0.95228 | 0.452 | +0.0271 | **0.8317 (board max)** |
| `fid` | 0.51239 | 0.80964 | **0.297** | **+0.0178** | 0.4967 |
| `reach` | 0.07902 | 0.96414 | 0.885 | +0.0531 | 0.1815 |
| `jac` | 0.03038 | 0.39727 | 0.367 | +0.0220 | 0.0321 |

`score_avg` is the plain mean of the six scaled members (max |diff| 3.3e-16 over all 780),
so **Δoverall = Δscaled/6** exactly.

Gap to the leader (Aginglab.com, `0.18991`), in overall points:

| member | reach | nmae | fid | mse | jac | pds |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| they gain | **+0.0162** | +0.0142 | +0.0095 | +0.0068 | −0.0015 | −0.0047 |

## 2. Amplitude is not the lever — this is now settled

Regressing each raw metric on `log(raw mse)` across the top 60 teams:

| metric | slope /log-unit | R² | **our residual** | z |
| --- | ---: | ---: | ---: | ---: |
| reach | +0.0191 | **0.097** | **−0.0056** | −0.15 |
| fid | +0.0094 | 0.165 | +0.0025 | +0.18 |
| pds | −0.0218 | 0.164 | +0.0499 | +1.54 |

Within teams, over 288 consecutive submission pairs, **`Δscore_avg` is uncorrelated with
`Δlog mse` (r = −0.07, p = 0.24)**. Doubling amplitude buys +0.013 reach ≈ +0.0025 overall.

Two board facts finish the argument:

* **We527** ships `STATE baseline` at raw mse `1.3577` — within 2.6% of our `1.3236` — and
  scores reach `0.2529` against our `0.1815`.
* **Aginglab** went `GeroAI_v9 → v11` by *halving* amplitude (mse `1.822 → 0.959`) and
  *gaining* `+0.100` reach.

So Addendum 5's "alpha's marginal direction is up" is **retired**. The v13 α=2.0 package was
built on 2026-08-25 (`s3://…/output/exp089-production-pm-cd4/20260825T193634Z/`) and never
submitted; on this evidence it should stay unsubmitted.

## 3. The defect

`reach` = mean over perturbations of `k*/n_real`. The pool is the **reference's** significant
set; it is ordered by **our** confidence (our significant calls first, then ascending
predicted `p_adj`, then `p`, then `|lfc|`); `k*` is the **deepest** prefix whose sign purity
clears `0.9`. The pool and each gene's sign-correctness are fixed by the model. **The only
free variable is the order.**

On the HCT116 tables, ordering counterfactuals (call set held at production):

| ordering | raw reach | Δoverall |
| --- | ---: | ---: |
| production (`sig` first, then `p_adj`) | 0.01465 | — |
| random permutation | 0.01194 | −0.00051 |
| by `\|lfc\|` | 0.04107 | +0.00497 |
| our significant calls **last** | 0.05125 | +0.00689 |
| oracle (correct signs first) | **0.59987** | **+0.11020** |

**Putting our own significance flag last beats putting it first.** The reason:

| | sign accuracy on the reference's DE set |
| --- | ---: |
| genes we **do** call significant | **0.3446** |
| genes we **do not** call significant | **0.8259** |

paired `t = −54.0`, `p = 5.6e-90`; the inversion holds in **99.2%** of perturbations
individually. Within-perturbation deciles by our `p_adj` run
`0.371, 0.267, 0.229, 0.270, 0.402, 0.597, 0.752, 0.816, 0.857, 0.902` — monotone in the
wrong direction over the last six.

## 4. The mechanism

The reference's DE genes are **90.8% down-regulated**. Our predictions are balanced (46.8%
up). Conditional accuracy: predicting **down** is right **92.2%**, predicting **up** is right
**10.8%**. And our most-confident decile is **83.5% up-predictions** against 46.8% overall.

That skew is the artefact `src/vcc2026/donor_transfer.py:261` already described and nobody
connected to scoring: a minimum-variance emitter places 400 cells near the predicted **mean**,
which for right-skewed counts sits above the **median**, so the rank-sum against 18,400 real
controls reports "significantly up" at zero true effect. Control expression does not explain
it — the partial `Spearman(p_adj rank, correct | ctrl CPM)` is still `+0.4521`.

The closed form confirms it. With `n = 400`, `m = 18,400`:

```
z_g = sqrt(12·n·m/(N+1)) · (F̄_g − ½) / sqrt(tc_g),    F̄_g = mean of the control CDF
                                                            over our emitted count law
```

`|F̄_g − ½|/sqrt(tc_g)` predicts our own p-order at **Spearman +0.996** — better than a
re-run of the emitter with a different seed reproduces it (+0.994). It is computable **before
a single cell is emitted**. And `|predicted effect|` — the key `--keep-top` ranks on —
correlates with it at **+0.071**. Our ranking is very nearly independent of our own model's
signal and is instead a ranking of each gene's control-distribution skewness.

Two corrections fall out: no gene saturates (0 of 10,258; `pred_p == 0` never occurs), so the
`|lfc|` tiebreaker is dead code, and `experiments/093-spread-truncation/RESULTS.md` §2's claim
that the quantile key "saturates at 0/1 for thousands of genes" is wrong.

## 5. Why the channel is free

`emit_counts(mode="quasi")` repairs each gene's panel total **exactly**, so the pseudobulk —
and therefore `pds`, `expr_mse` and `nmae` — cannot move. Board-verified in v14: `pds` moved
by **exactly 0.000000** and `nmae` by 1e-5 across a full-range spread change. Only `fid`,
`jac` and `reach` are exposed.

## 6. What v14 and v15 actually measured

| | tail | keep | keep spread | n_pred | fid | reach | jac | overall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| **v11** | — | — | 0.25 | ~5340 | **−0.0357** | 0.1041 | 0.0031 | **0.150183** |
| v15 | 1.0 | 2000 | 0.25 | ~2370 | −0.0533 | 0.1155 | 0.0048 | 0.149436 |
| v14 | 1.0 | 850 | 0.25 | ~2317 | −0.0712 | 0.1178 | 0.0073 | 0.147273 |

Both arms moved `reach` and `jac` the right way and gave it all back on `fid`. v14 and v15
sit at nearly the same `n_pred` (+2%) yet differ by `+0.0053` raw `fid`, so that pair is a
*composition* effect, not a size effect. Inverting the size step v15 → v11 with the harness
precision curve puts coverage at **≈0.93 at n_pred 2370 and ≈1.0 at 5340**: `fid`'s coverage
term is **already saturated at v11**, so growing the call set further cannot pay, and
shrinking it costs coverage faster than precision recovers.

Both arms raised the tail and left the kept genes at `SPREAD`, which only ever *shrinks*
`n_pred`. That was the wrong end of the channel to move.

## 7. The arm

`--keep-spread` (new) decouples the two ends. Arm **a1-order-k2000** is

```
--tail-spread 0.25 --keep-top 2000 --keep-spread 0.0
```

The tail is bit-identical to v11, so the call set barely moves and `fid`'s coverage term is
left alone; the 2,000 genes per target carrying the largest predicted effect emit with **zero**
Poisson scatter, which lowers their p-values and lifts them up the confidence order that
`reach` walks. Their sign accuracy is `0.677` against `0.542` pooled, so promoting them raises
the purity of the head — **without** betting on the 90.8% down-asymmetry, which is not yet
confirmed to transfer off this harness.

## 8. Caveats, stated plainly

* **The HCT116 harness is reach-blind.** Raw reach `0.0147` there against `0.1815` on the
  board; `k*` is **1** on the median perturbation and exactly 0 on 31.8% of them. The
  *direction* of every ordering result above is trustworthy; the *magnitudes* are not.
  Reach work belongs on `h1-test` or the board.
* The 90.8% down-asymmetry is measured on one reference. Arm a1 deliberately does not depend
  on it.
* §6's coverage inversion borrows the precision curve from the harness. It is an inference,
  not a measurement.

## Files

| file | what |
| --- | --- |
| `r1_anatomy.py` | ordering counterfactuals incl. oracle and anti-oracle |
| `r2_signal.py` | sign accuracy by call status and by within-perturbation decile |
| `r3_mechanism.py` | global-flip, direction, expression-level and per-perturbation checks |
| `r4_realizable.py` | fid/jac/reach for realizable (call-set, ordering) pairs on board anchors |
| `aws/fastbuild.sh` | 18-minute build path reusing the cached worker tree |

---

# Addendum, same day — the asymmetry does NOT transfer, and §3–4 must be read narrowly

Sections 3 and 4 above were measured entirely on the Orion HCT116 harness. Two
measurements taken afterwards change how they should be read. Both are recorded here
rather than edited away, because the reasoning that produced the wrong design is the
part worth keeping.

## A. h1-test says the reference's DE set is BALANCED

`data/harness/h1-test.h5ad`, 40 perturbations with ≥50 cells, 12,000 control cells,
scorer-exact gate and Wilcoxon path, `lfc = log2((mean CPM_pred + 1e-9)/(mean CPM_ctrl + 1e-9))`:

| dataset | fraction of reference-significant genes with **negative** lfc | n_real median |
| --- | ---: | ---: |
| Orion HCT116 | **0.908** | 304 |
| **h1-test, padj < 0.05** | **0.538** (median 0.523, IQR [0.459, 0.613]) | 253 |
| h1-test, padj < 0.01 | 0.529 | — |
| h1-test, top 100 by padj | 0.495 | — |

exp088 §5 established h1-test as the **reach-representative** harness (raw reach
`0.23–0.24` against the board's `0.169`; HCT116 is `0.015`). On that harness the
reference is a coin flip on direction.

**So HCT116's 90.8% is a property of that one screen, not of the assay.** The
`+0.5009` "inverted ranking" in §3 is the product of two facts multiplied together —
our ordering is up-first, and *that* reference is overwhelmingly down. On a balanced
reference an up-first ordering is not systematically wrong, it is merely uninformative.

**Everything in §3–4 that depends on the down-asymmetry is therefore retired**:
the `sig_last` result, the "predicting down is right 92.2%" split, and the whole
down-first ordering family. What survives is narrower and still true: the ordering
channel is free, `|predicted effect|` is nearly orthogonal to the order the scorer
actually walks (Spearman +0.07), and the order is analytically pre-computable.

## B. On the real challenge artefact, `spread` is a SIGN knob, not a confidence knob

Auditing `submissions/v11.vcc` directly against the real `data/challenge/context_*.h5ad`
controls — no ground truth needed, every quantity is a property of our own prediction:

* **`n_pred` ≈ 6,000** per perturbation (A 6,261 / B 6,066 / C 5,676 median), which is
  **61% of the entire gated universe** and ~12% above the 5,340 the HCT116 harness implied.
* Our confidence ranking is essentially a **signed** lfc ranking: within-perturbation
  decile 1 is **85–89% up**-predictions, decile 10 is **5–11%**, monotone throughout, and
  `Spearman(ascending p-rank, descending signed lfc) = +0.69`.
* `Spearman(p_adj rank, |lfc| rank) = +0.015`. The two orderings are unrelated — which is
  why `--keep-top`, which ranks on `|effect|`, could never have steered the order.
* In **linear** CPM space **100.0%** of gated genes have mean > median, and the control's
  own mean sits at the **58.3rd percentile**. Emitting 400 cells at the predicted mean
  starts every gene 8.3 percentile points above the median — a free "up" at zero effect.
  (In log1p space only 37–42% are right-skewed; the log framing is a red herring.)
* Our emission is **under-dispersed ~2.3×**: per-gene CV 0.38–0.40 across our 400 cells
  against 0.84–0.92 across the real controls.
* **75–81% of our significant calls have `|lfc| < 0.10`**, and ~35% are internally
  contradictory — the rank-sum earns its p-value in the up direction while the reported
  lfc says down.

Uniform spread sweep on one real challenge perturbation:

| spread | n_sig | down-conversion | sign agreement (z vs lfc) |
| ---: | ---: | ---: | ---: |
| 0.00 | 7,993 | 0.68 | 0.55 |
| **0.25 (v11)** | 6,515 | 0.44 | 0.65 |
| 0.50 | 4,086 | 0.13 | 0.83 |
| 1.00 | 1,270 | 0.03 | 0.90 |

Low spread concentrates the cells at the mean, which is above the control median, so the
gene **reads up**; high spread drops the emitted median at fixed panel total, so it
**reads down** — and suppresses significance far faster than it flips sign.

**This killed arm a1 before it was submitted.** a1 gave `spread = 0.0` to the top 2,000
genes by `|effect|`, a set that is ~61% down / 39% up, so it would have driven a
directionally mixed block to the head of the ordering — a head that cannot clear the 0.9
purity bar at any plausible accuracy. The build was stopped and never submitted.

## C. What this leaves

The sober reading: **our `reach` of 0.1815 is close to what our model's per-gene
directional confidence actually supports, and there is no free re-ordering trick.** The
teams at 0.25–0.27 have a genuinely better confidence signal — and the board's own
strongest pattern is that every entry naming **STATE** lands at reach 0.21–0.27 across
raw mse spanning 1.36 to 1,774, i.e. amplitude-invariant and model-driven.

Retired today, with the measurement that retired each:

| lever | verdict | evidence |
| --- | --- | --- |
| amplitude / alpha | **dead** | Δscore ⊥ Δlog mse over 288 board pairs (p = 0.24); exp088 α 1.2 → −0.0037, 2.0 → −0.0610 |
| `min_cpm` gate | **at its peak** | exp088b: 8 → +0.0028, **10 → +0.00455**, 15 → −0.0056, 30 → −0.0129 |
| call-set truncation | **at its peak** | v11 (0.150183) > v15 (0.149436) > v14 (0.147273) |
| down-first ordering | **invalid** | h1-test reference is 53.8% down, not 90.8% |
| `\|effect\|`-keyed ordering | **cannot steer** | Spearman(p_adj rank, \|lfc\| rank) = +0.015 |

Every scalar knob on this pipeline is now at or adjacent to its measured optimum. The
remaining headroom is in the **operator** and in a genuinely better per-gene directional
signal — not in emission parameters.

That makes **experiment 091** the one measured, positive, unshipped result we hold, and it
is shipped today: `data/donors/whitener_union_pruned.npz`, sha256 `e2f65ecf…`, 4,393
common-essential rows dropped from the basis, `dropped_in_panel: []`.


---

# Control result: the fast build path is bit-identical to v11

`aws/fastbuild.sh` with no flags produced a package hashing to

```
24a1282b1ca69e17a5b1dfa3c4279817dedff725b0f885fdf785e4b28e5c22ef
```

which is **byte-for-byte `submissions/v11.vcc`**, built on 2026-08-24 by the exp089
ceremony runner on a different day. Three things follow, and all three matter:

1. **The 18-minute path is the 3h20m path.** The ceremony runner's extra ~3 hours were
   entirely re-downloading the donor tree, the controls bundle and the Python environment
   around a build whose own `/usr/bin/time -v` reads 16:55 (12:08 on a warm cache). Nothing
   it verifies is lost: the whitener sha256 is still asserted inside the builder, the test
   suite still runs fail-closed before every build, and each arm still writes its own
   `build-report.json`.
2. **The pipeline is exactly deterministic** across machines and days — every RNG stream,
   the library-size draw, and the emission all reproduce.
3. **This session's code changes are provably inert by default.** `--keep-spread` and the
   `WHITENER_REGISTRY` / `--whitener-id` rework changed the builder, and an unflagged build
   still emits the same bytes. That is the guarantee `test_emission_spread_is_the_shipped_
   scalar_unless_a_tail_is_asked_for` asserts at unit level, now confirmed end to end on a
   3.9 GB artefact.

It also means the final slot's fallback costs nothing: `submissions/v11.vcc` on disk *is*
the rebuilt artefact, so restoring our best-known submission needs no build and no download.


---

# Arm a3 — the ordering channel, resurrected on evidence

Arm `a1-order-k2000` was withdrawn mid-build on the grounds that a directionally mixed head
cannot clear the 0.9 purity bar. **That objection was wrong**, and experiment 095 says so:
purity is computed on the **lfc sign**, which `spread` provably cannot move, so promoting a
mixed block is not self-defeating. The premise that *was* wrong — the `+0.50` inverted
ranking — stays retracted.

What replaced it is a measurement rather than an argument. On the depth-repaired harness:

| | Δoverall |
| --- | ---: |
| ordering the reference's DE pool by `\|predicted lfc\|` | +0.00448 |
| **promoting the top 2,000 by `\|effect\|` ahead of the rest** | **+0.00583** |
| the same on the broken harness | +0.00497 |

The `|lfc|` result agrees to `0.0005` across a harness repair that moved the reference's
direction balance from 0.908 to 0.609, so it is not a confound artefact. And block promotion
**beats** a full `|lfc|` sort, because it keeps the emitted p-order *inside* the promoted
block — the right policy is "promote what the model is confident about, then leave the
emitted order alone", not "sort everything by magnitude".

**Arm a3** is therefore `a1` unchanged: v11 base, shipped operator,
`--tail-spread 0.25 --keep-top 2000 --keep-spread 0.0`. The tail is bit-identical to v11, so
this is a single-variable measurement against a known board point (`0.1501827`).

It is deliberately **not** built on `a2-prune`. Stacking two untested changes would make
neither attributable, and the board is noiseless for a fixed artefact, so two clean
single-variable measurements are worth more than one combined guess.

## What a3 actually tests, and why the harnesses cannot answer it

`--keep-spread 0.0` does two things, and the promotion curve credits only the first:

1. it promotes the top-2,000 up the confidence order — worth `+0.00583` above;
2. only **37.9%** of them were already significant, so it can add ~1,240 calls per
   perturbation, growing `n_pred` by up to **31%**.

Every harness we hold scores (2) at about `−0.005`, because `n_real` is 322 / 68 / 222 and
`fid`'s coverage term is pinned at 1, so shrinking always looks good. **The board says the
opposite** — raw `fid` `0.49126 → 0.49654 → 0.50178` as `n_pred` goes `2317 → 2370 → 5340`.

a3 is the first measurement of that gradient **above** v11, in the only regime that counts.
Both of its effects are free on `pds`, `mse` and `nmae` by construction, so whatever it
returns is a clean reading of the `fid`/`jac`/`reach` surface.


## a3, measured on the real challenge controls before spending a slot

`probe_de.py` recomputes the scorer's Wilcoxon against `data/challenge/context_*.h5ad`.
No ground truth involved — `n_pred` is a property of our own prediction and the real
controls. Median over 12 perturbations per context:

| arm | ctx A | ctx B | ctx C | mean |
| --- | ---: | ---: | ---: | ---: |
| `v11-control` | 6,098 | 6,192 | 5,592 | 5,961 |
| **`a3-order-k2000`** | 6,264 | 6,486 | 5,828 | **6,193** |

**+3.9%, not the +31% projected from the harness.** The harness projection assumed only
37.9% of the top-2,000 were already significant, which is true there — but on the challenge
we already call ~6,000 of ~9,100 gated genes, so nearly all of the promoted block was inside
the call set already and the promotion mostly re-ranks rather than adds.

That collapses a3's `fid` exposure from about `−0.005` to about `−0.0004`, and makes it
what it was meant to be: **a near-pure measurement of the ordering channel.** a4, the
`n_pred`-neutral variant staged in `aws/launch_a4.sh`, is not needed and was not built.

Also note the challenge `n_pred` of ~5,961 for v11 against the ~5,340 the HCT116 harness
implied — one more place that harness understated the real geometry.

*(The probe's own `lfc` used `expm1(mean(log1p(CPM)))`, a geometric mean, where cell_eval2
uses the arithmetic mean of per-cell CPM. That biases its reported up/down split and it has
been corrected in the script; `n_pred` comes from the rank test and is unaffected, so the
table above stands.)*


## The board's own x-axis was wrong too

Probing the v14 and v15 h5ads still on the worker gives their **real** challenge `n_pred`,
which is not what the HCT116 harness implied:

| artifact | harness estimate | **real challenge** | board raw `fid` |
| --- | ---: | ---: | ---: |
| v14 (`keep 850`) | ~2,317 | **1,563** | 0.49126 |
| v15 (`keep 2000`) | ~2,370 | **1,911** | 0.49654 |
| v11 | ~5,340 | **5,961** | 0.50178 |
| **a3** | — | **6,193** | *pending* |

Every earlier analysis of the call-set channel — including experiment 093's BOARD-RESULT and
this file's own first pass — placed those points on a fabricated axis. v14 and v15 are 33%
and 19% *lower* than believed, and v11 is 12% higher.

Re-inverting with the real numbers, using the repaired harness's `fid(m)` as a precision
shape (valid there because its coverage is pinned at 1), level-shifted by `0.9493`:

| entry | real `n_pred` | raw `fid` | implied precision | implied coverage |
| --- | ---: | ---: | ---: | ---: |
| v14 | 1,563 | 0.49126 | 0.5489 | **0.8950** |
| v15 | 1,911 | 0.49654 | 0.5422 | **0.9157** |
| v11 | 5,961 | 0.50178 | 0.5018 | **1.0000** |

A 22% call-set increase (v14 → v15) bought `+0.017` coverage; a 3.1× increase (v15 → v11)
bought `+0.088`. **So the whole rise across our three board points is coverage, not
precision**, and it is spent by v11. Projected past it:

| `n_pred` | precision = `fid` | Δoverall vs v11 |
| ---: | ---: | ---: |
| **6,193 (a3)** | **0.5010** | **−0.00044** |
| 7,000 | 0.4984 | −0.00188 |
| 8,000 | 0.4952 | −0.00366 |

**a3's `fid` exposure is −0.0004 against a modelled reach gain of +0.0058.** It also
delivers the first board reading of the `fid` gradient on the far side of coverage
saturation, which no harness we hold can produce.


## Arm a6 — the graded map on the agreement key, built to be probed

Experiment 095 §6 found the best buildable ordering key (`√(|A|·|B|)` gated on cross-donor
sign unanimity, `+0.00678` as a full ordering) and found that the *block* mechanism cannot
spend it. a6 builds the mechanism that can:

```
--whitener-id union_pruned --tail-spread 0.25 \
--spread-map graded --keep-key agreement --spread-lo 0.00 --spread-hi 0.50
```

`spread_lo`/`spread_hi` are chosen so the **mean** spread over the gated axis is `0.25`,
which is v11's uniform value — a linear ramp from 0.00 to 0.50 averages 0.25. The intent is
to hold `n_pred` near production while letting the realised p-value order follow the
agreement key rather than the emitted-noise order.

**It is built to be probed, not assumed.** Two things about it are genuinely unmeasured:

1. **`n_pred` under a graded map.** The block map was probed at `+3.9%` on the real
   challenge controls; nothing analogous exists for a ramp, and `fid`'s coverage term is
   saturated at v11 so a large move costs.
2. **Mechanism fidelity.** The `+0.00678` was a counterfactual that assumed the p-order
   follows the key *exactly*. A graded spread ramp only approximates that — the realised
   order is a compromise between the key and the underlying effect sizes. The block map's
   `+0.00583` is the more faithfully realised of the two.

So a6 is only submitted if its probe puts `n_pred` close to v11's ~5,961. Otherwise a5 —
whose mechanism is validated end to end — takes the slot. **The measurement slot 2 is
actually buying is "what is the ordering channel worth on the board", and both arms deliver
that**; a5 delivers it with less variance.

One known blemish: a6 was already building when the single-donor fallback landed, so its 6
targets carried by only one donor rank by gene index rather than by `|effect|`. Worth about
`−0.0002` overall — inside what a6 is measuring.

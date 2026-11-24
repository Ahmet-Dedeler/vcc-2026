# 317 — Can OUR pipeline trade `mse` for `nmae` + `reach`? The operating point as the variable.

**Date:** 2026-09-07/08 · preregistration `PREREGISTRATION.md`, committed `98b2d9a` **before any
build** · instrument: the real `cell-eval2` 0.16.0 (`--preset vcc2026`, DE backend `scanpy`) on
`data/harness/h1-test.h5ad`, eleven disposable EC2 workers, **$17.18 total, $1.56 per
six-member score, all eleven terminated, nothing left running** · **no submission built, nothing
under `submissions/` touched, no `vcc` command run, no local `cell-eval2`, local peak RSS < 1 GB.**

> **317 (operating-point trade `mse` → `nmae`+`reach`): DOES NOT EXIST as a competition lever.**
> **Peak composite +0.000131 overall on the emission axis (7 % of the +0.0019 shape-only gate,
> 1.4× a re-seed) and +0.002374 on the amplitude axis even with `mse` counted FREE (27 % of the
> +0.0089 gate an amplitude change must clear), against 0.2177 shipped. It pays `mse` at
> −0.2183…−1.5210 raw and buys `nmae` −0.0029 raw and `reach` +0.0076 raw, both with 95 % paired
> intervals through zero — an exchange of about 35 : 1 against us.** §4.1c is re-established on
> stronger evidence than it had.

---

## Gate 1 — the instrument

### 1.1 What can score six members, and what it costs

| instrument | six members? | cost / evaluation | status here |
|---|---|---|---|
| **cloud `cell-eval2`** (`scripts/cloud_score.sh`, r7i.2xlarge on-demand) | **yes, RAW** | **$1.56 measured (221: $1.48), ≈ 2.8 h wall, zero local RSS** | **the instrument used** |
| local `cell-eval2` | yes | free, ~25 GB RSS | **forbidden by this agent's standing rules** (one instance, shared machine). Never run. |
| h1-test / `instrument273.py` | `pds` only (5.6e-17) | free | cannot answer a six-member question |
| PLOLO | `pds` of geometry arms | free | wrong channel |
| the board | all six | one of two daily slots | not spent — 317 is a measurement, not a candidate |

**The SCALED composite is not available on h1-test**, and rebuilding it was not worth the budget.
`cell-eval2 score` needs a *real bundle* (`prep-real-bundle`, five disjoint half-splits, ≈ 10× one
`run`). 221 tried twice and both failed —
`s3://…/output/cloud-score/bundle-20260904/STATUS` still reads `FAILED 1 (completed=0)`, and
`ka2-20260904` burned 157 min on the `--baseline`-is-a-directory bug. **So the deliverable is the
RAW six-member row per operating point, converted to board overall points through the board's own
scaling map (§1.4).**

### 1.2 KA-1 — the instrument reproduces a recorded six-member row, bit for bit

221's `rc_callset.h5ad` was re-scored from the same S3 object on a fresh worker:

| member | 221 recorded | 317 re-score | \|diff\| |
|---|---:|---:|---:|
| `pds` | 0.873809523810 | 0.873809523810 | 0.00e+00 |
| `mse` | 1.015977751768 | 1.015977751768 | 0.00e+00 |
| `nmae` | 0.960046422899 | 0.960046422899 | **4.44e-16** |
| `fid` | 0.495058277089 | 0.495058277089 | 0.00e+00 |
| `reach` | 0.226704287761 | 0.226704287761 | 0.00e+00 |
| `jac` | 0.074203737608 | 0.074203737608 | 0.00e+00 |

**KA-1 PASS**, max \|diff\| 4.44e-16 (`work/sweep317.txt`). Five members bit-exact.

### 1.3 KA-2/3 — every arm is one knob off a byte-identical rebuild of a board-anchored arm

Each worker rebuilds `base` from `experiments/170-gamma-mse/work/pre.npz` + the shipped operator
and gates its sha256 against 221's. **10/10 workers reported PASS:**

```
[GATE 5] base sha256   048b9286030ff0124514134b55d935fdb79a4d002cc04d590d83676b32553c77
         221 rc_callset 048b9286030ff0124514134b55d935fdb79a4d002cc04d590d83676b32553c77  -> PASS
```

221's gates 1–3 were kept verbatim and passed on every worker: the a15 and a19 operator anchors
to ≤ 1e-11 relative, 170's h1 effect at `||d||F/||E||F 6.462e-15`, `1−min_row_cos 2.22e-16`, and
the CD4 f-linearity at `rel 2.22e-16`. `base` (= `rc_callset`) is itself a **board point**,
overall 0.1918638.

### 1.4 KA-5 — the board scaling map, and its out-of-sample retrodiction

`scaled = (raw − b)/(r − b)`, fitted per member on the 2026-09-07T2153Z snapshot (780 entries),
restricted to each member's **unclamped core** (690 of 780 rows sit exactly on `score_mse = 0`
and 36 on the `nmae` floor of −6; fitting through them reads `mse`'s span as −5.24 instead of
−0.95):

| member | d(scaled)/d(raw) | `b` | span `r−b` | R² | n | AGENTS.md §4 spec |
|---|---:|---:|---:|---:|---:|---|
| `pds` | +2.214250 | 0.500265 | +0.451620 | 0.999966 | 780 | b 0.500, span 0.43–0.48 ✓ |
| `mse` | −1.052616 | 0.987865 | −0.950014 | **0.999996** | 59 | b 0.986–0.992, span 0.95–0.96 ✓ |
| `nmae` | −1.657458 | 1.000193 | −0.603333 | 0.997781 | 689 | b 1.0009–1.0017, span 0.57–0.63 ✓ |
| `fid` | +3.362662 | 0.512426 | +0.297383 | 0.999927 | 780 | b 0.505–0.522, span 0.28–0.33 ✓ |
| `reach` | +1.128866 | 0.079013 | +0.885845 | 0.999968 | 780 | b 0.047–0.097, span 0.86–0.93 ✓ |
| `jac` | +2.725356 | 0.030380 | +0.366924 | 0.999847 | 780 | b 0.021–0.037, span 0.34–0.39 ✓ |

All six land inside the published spec ranges. `score_avg` is the plain mean of the six `score_*`
columns (max \|diff\| **3.3e-16** over 778 rows), so `Δoverall = Σ Δscaled / 6` exactly.

Applied **out of sample** to five of our own scored board rows (`a11`, `a14`, `a15`, `abc-knee`,
`cd4svd5`; raw and scaled both recorded in `docs/BOARD-LEDGER.md` / `docs/READ-FIRST.md` §1) it
retrodicts the board's published scaled columns with worst-member errors

```
pds 2.5e-03   mse 6.8e-05   nmae 1.7e-03   fid 5.4e-04   reach 2.0e-04   jac 2.1e-04
```

(`work/ka5-map-retrodiction.txt`). Those are level offsets that cancel in a delta; the slopes,
which are what price a delta, are fitted at R² ≥ 0.9978.

**One documented limit, because it changes a number.** The `mse` slope is valid on raw
`mse` ∈ [0.82, 0.94], where our arms live. Above ≈ 0.95 the per-context clamps flatten it, so
extrapolating it to the leader's 0.99153 over-charges their `mse` loss (−0.0307 overall against
the board's published −0.01793). **The leader's `mse` row must be read off the board, not off
this map** — §3 does exactly that.

### 1.5 Noise floor, and the experimental unit

`a15` vs `a15_seed1` — 221's re-seeded twin, identical construction, emission seed only, through
this same instrument. **This is ONE seed pair, not an sd**; the mandate's "034 sd 0.0004 /
0.0012" is from a different harness and a different emitter generation and is not the floor here.

| | `pds` | `nmae` | `fid` | `reach` | `jac` |
|---|---:|---:|---:|---:|---:|
| raw Δ (aggregate) | +0.000700 | +0.000048 | +0.000777 | −0.003319 | +0.000104 |
| paired per-target, mean [95 %] | +0.00070 [−0.00056, +0.00210] | +0.00005 [−0.00014, +0.00022] | +0.00078 [−0.00034, +0.00185] | **−0.00332 [−0.01585, +0.00841]** | +0.00010 [−0.00019, +0.00040] |
| n targets | 85 | 58 | 85 | 82 | 85 |

Priced through §1.4 the twin is **+0.000096 overall**. Note `reach`'s ±0.016 paired interval —
**`reach` is the noisiest member on this surface by an order of magnitude**, and no `reach`
result below that width should be read as a signal.

**Experimental unit for every interval below: one perturbation (target gene) on h1-test, PAIRED
between the arm and `base`** — same effect matrix, same library draw, one knob moved — n ≤ 85,
10,000-resample percentile bootstrap, seed 0. `nmae` reads n = 58: the scorer omits 5 targets for
an empty reference gate and 28 for fewer than 10 gated genes, and that set is reference-side only,
so it is identical for every arm. `reach` reads n = 82 (3 targets have no confident reference DE
genes). **`expr_mse_unbiased_capped_norm` is a single panel-level aggregate with no per-target
value, so it gets no interval — only the aggregate delta.**

---

## Gate 0 — the structural result, which decides half the question before any score

`emit_counts(mode="quasi")` computes each gene's panel total as `floor(n·L·profile_g)` and then
**repairs the scattered part back to that exact integer**. Measured on the real emitter
(`invariance317.py`, 18,080 genes × 400 cells, `work/gate0-invariance.txt`):

```
overdispersion 0.0 / 1.0 / 2.0 / 8.0 / 16.0   totals identical: True   sum|d| 0
keep_spread 0.25 -> 0.10                      totals identical: True   sum|d| 0
tail_spread 1.0  -> 0.70                      totals identical: True   sum|d| 0
keep_top    8000 -> 16000                     totals identical: True   sum|d| 0
emission seed 7 -> 11 (libraries fixed)       totals identical: True   sum|d| 0
NEGATIVE CONTROL: profile changed             totals identical: False  sum|d| 1,022,505
```

Confirmed in the sweep's own builds (gate 4 on every worker) and then **on the real scorer**:

```
od1 / od8 / ks010 / ts070 / kt16000        vs base :  85/85 blocks digest-identical  PASS
amp070 / amp085 / amp120 / amp145 / amp200 vs base :   0/85                          PASS
```

**Every emission arm returned `Δpds = +0.000000` and `Δmse = +0.000000` from `cell-eval2`, and the
paired per-target `pds` interval is [+0.000000, +0.000000] on 85/85 targets.** So **no emission
knob can pay the `mse` side of the trade.** Of the six knobs the mandate names, exactly one —
amplitude — moves `mse` at all.

The one route by which an emission knob *could* touch `mse` is the jackknife `claim` term
(AGENTS.md §8.1c: *"`mse` is not invariant — the jackknife term reads our emitted noise"*). That
channel is **already at its cap**: 211 measured `claim` at 1.000000 / 0.999904 / 1.000000 of cap
on A/B/C. On h1-test it is unreadable in either direction because the surface's raw `mse` is
**1.016**, past the board's `b` = 0.9879, where `score_mse` clamps to 0.

**The amplitude knob is what it says it is.** `--effect-scale s` is applied in
`profile_from_effect` on the delivered log-effect — exactly the quantity §4.1's `δ̂ = α·d`
algebra is written on. Its equivalence to `VCC_ALPHA → 0.982·s` was checked on the production
functions (`amp_equiv317.py`): `apply_v9_whitening`'s preserve factor is a ratio of norms of the
same input, so it is **exactly scale-invariant** and its output is exactly `c ×` the output at
`E` (max rel 1.4e-15 over c ∈ {0.70, 0.85, 1.20, 1.45, 2.00}); the CD4 residual added at a solved
*fraction* is equivariant to **≤ 9.2e-4 relative**, and on h1-test the CD4 bank covers 4/85
targets. So the two knobs agree to ~1e-3 on one leg touching 4 targets.

---

## Gate 2 — the sweep. All six members at every operating point.

**RAW six members on h1-test**, one variable per row, `base` = 221's `rc_callset` rebuilt
byte-identically:

| arm | variable | `pds` | `mse` | `nmae` | `fid` | `reach` | `jac` |
|---|---|---:|---:|---:|---:|---:|---:|
| **base** | shipped | 0.873810 | 1.015978 | 0.960046 | 0.495058 | 0.226704 | 0.074204 |
| `od1` | `--overdispersion` 4.0 → 1.0 | 0.873810 | 1.015978 | 0.960046 | 0.495058 | 0.227219 | 0.074188 |
| `od8` | `--overdispersion` 4.0 → 8.0 | 0.873810 | 1.015978 | 0.960047 | 0.495054 | 0.225087 | 0.074204 |
| `ks010` | `--keep-spread` 0.25 → 0.10 | 0.873810 | 1.015978 | 0.960047 | 0.495058 | 0.226826 | 0.074204 |
| `ts070` | `--tail-spread` 1.0 → 0.70 | 0.873810 | 1.015978 | 0.960046 | 0.495059 | 0.227507 | 0.074158 |
| `kt16000` | `--keep-top` 8000 → 16000 | 0.873810 | 1.015978 | 0.960046 | 0.488690 | 0.229049 | 0.064934 |
| `amp070` | amplitude ×0.70 | 0.876050 | 0.944540 | 0.970778 | 0.494617 | 0.228306 | 0.074205 |
| `amp085` | amplitude ×0.85 | 0.875630 | **0.930874** | 0.964700 | 0.494908 | 0.228406 | 0.074211 |
| `amp120` | amplitude ×1.20 | 0.873810 | 1.234248 | **0.957129** | 0.495451 | 0.228400 | 0.074191 |
| `amp145` | amplitude ×1.45 | 0.873249 | 1.567668 | 0.957843 | 0.496073 | **0.234329** | 0.074134 |
| `amp200` | amplitude ×2.00 | 0.872549 | 2.537005 | 0.996304 | **0.498828** | 0.228292 | 0.074224 |

**Deltas vs `base`, priced in board overall points** (transfer factor 1.0 = the optimistic
reading; the spread-policy arm class's measured h1-test → board `reach` factor is 0.39–0.4277,
i.e. h1-test *overstates* by 2.3–2.6×, so the second column is the honest one for emission arms):

| arm | class | d`pds` | d`mse` | d`nmae` | d`fid` | d`reach` | d`jac` | **d overall** | ×0.4277 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `od1` | emission | +0.000000 | +0.000000 | −0.000000 | −0.000000 | +0.000514 | −0.000015 | **+0.000090** | +0.000038 |
| `od8` | emission | +0.000000 | +0.000000 | +0.000000 | −0.000005 | −0.001617 | +0.000000 | −0.000307 | −0.000131 |
| `ks010` | emission | +0.000000 | +0.000000 | +0.000000 | +0.000000 | +0.000121 | +0.000000 | +0.000023 | +0.000010 |
| `ts070` | emission | +0.000000 | +0.000000 | −0.000000 | +0.000001 | +0.000803 | −0.000046 | **+0.000131** | +0.000056 |
| `kt16000` | emission | +0.000000 | +0.000000 | −0.000000 | −0.006368 | +0.002345 | −0.009270 | **−0.007338** | −0.003139 |
| `amp070` | amplitude | +0.002241 | −0.071438 | +0.010732 | −0.000441 | +0.001601 | +0.000001 | +0.010450 | — |
| `amp085` | amplitude | +0.001821 | −0.085104 | +0.004654 | −0.000150 | +0.001702 | +0.000007 | **+0.014556** | — |
| `amp120` | amplitude | +0.000000 | +0.218270 | −0.002918 | +0.000393 | +0.001696 | −0.000013 | −0.036953 | — |
| `amp145` | amplitude | −0.000560 | +0.551690 | −0.002203 | +0.001015 | +0.007625 | −0.000070 | −0.094413 | — |
| `amp200` | amplitude | −0.001261 | +1.521028 | +0.036257 | +0.003769 | +0.001588 | +0.000021 | −0.274903 | — |

**Paired per-target intervals** (unit: one perturbation; `*` = 95 % interval excludes zero). Full
table in `work/sweep317.txt`; the ones that matter:

| arm | `nmae` | `fid` | `reach` | `jac` |
|---|---|---|---|---|
| `od1` | −0.000000 [0, 0] | −0.000000 [−1.4e-5, +1.3e-5] | +0.000514 [−0.000665, +0.001716] | −0.000015 [−4.6e-5, 0] |
| `ts070` | −0.000000 [0, 0] | +0.000001 [−0.000114, +0.000104] | +0.000803 [−0.000012, +0.001903] | −0.000046 [−0.000144, +0.000007] |
| `kt16000` | −0.000000 [0, 0] | −0.006368 [−0.014557, +0.000441] | +0.002345 [−0.000314, +0.006422] | **−0.009270 [−0.013499, −0.005692] \*** |
| `amp120` | −0.002918 [−0.005784, +0.000012] | +0.000393 [−0.000077, +0.000860] | +0.001696 [−0.000126, +0.005229] | −0.000013 [−4.6e-5, +1.3e-5] |
| `amp145` | −0.002203 [−0.008772, +0.004275] | **+0.001015 [+0.000097, +0.001985] \*** | +0.007625 [−0.002289, +0.025284] | −0.000070 [−0.000198, +0.000049] |
| `amp200` | **+0.036257 [+0.013584, +0.061536] \*** | **+0.003769 [+0.001469, +0.006641] \*** | +0.001588 [−0.000353, +0.005213] | +0.000021 [−0.000494, +0.000613] |

**`nmae` is EXACTLY pinned on every emission arm** (`−0.000000 [0, 0]` on 58 targets) — 227's
structural closure, reproduced as a paired interval on the real scorer.

### 2.1 The amplitude curves, and an independent confirmation of 187

Six points on `s` (`ampcurve317.py`, `work/ampcurve317.txt`):

| member | s=0.70 | 0.85 | **1.00** | 1.20 | 1.45 | 2.00 | turning point |
|---|---:|---:|---:|---:|---:|---:|---|
| `mse` | 0.944540 | **0.930874** | 1.015978 | 1.234248 | 1.567668 | 2.537005 | minimum on the grid at **s ≈ 0.85** |
| `nmae` | 0.970778 | 0.964700 | 0.960046 | **0.957129** | 0.957843 | 0.996304 | parabola vertex **s = 1.2038** (max resid 2.1e-3) |
| `reach` | 0.228306 | 0.228406 | 0.226704 | 0.228400 | **0.234329** | 0.228292 | vertex s = 1.484 (max resid 3.6e-3 — inside its own noise) |
| `fid` | 0.494617 | 0.494908 | 0.495058 | 0.495451 | 0.496073 | **0.498828** | monotone up |
| `jac` | 0.074205 | 0.074211 | 0.074204 | 0.074191 | 0.074134 | 0.074224 | flat (range 9e-5) |
| `pds` | 0.876050 | 0.875630 | 0.873810 | 0.873810 | 0.873249 | 0.872549 | flat (range 3.5e-3 over a 2.9× amplitude span) |

- **187's `nmae` turnover is confirmed on this surface, independently.** 187 (board-side algebra,
  2026-09-03) predicted a uniform-`k` increase is safe on `nmae` "up to ≈ +20 %, turnover at
  +29 %". **Measured here: the minimum is at s = 1.20 (+20 %) and `nmae` is already worse at
  s = 1.45.** The total collectable `nmae` gain over the whole axis is **0.00292 raw**, against
  187's "0.0012 raw inside the `mse`-feasible window" — same order, ours larger only because
  this reading ignores the `mse` constraint 187 respected.
- **`pds` is not bit-exactly scale-invariant through a count emitter.** §4.1's algebra is exact on
  the continuous profile; the emitter floors counts, so a 2.9× amplitude span moves `pds` by
  3.5e-3. Small, but it is not zero, and the amplitude arms are correctly classified as
  amplitude-changing rather than shape-only.

### 2.2 Instrument limitation, stated plainly

**h1-test's own `mse` optimum is at s ≈ 0.85, not at the shipped point.** That is *not* a board
opportunity and `amp085`'s +0.014556 must not be read as one. h1-test must be scored with the H1
union OFF (it is oracle-leaked, 85/85 targets — AGENTS.md §6.4), so the prediction on that
surface is materially weaker than the board's; its raw `mse` at the shipped point is **1.016**
against the board's **0.81671**, and a weaker direction wants a smaller amplitude. **h1-test
cannot locate the board's amplitude optimum**; 187 and 211 locate it, at 0.982, where we ship.
What h1-test *does* measure, and what the hypothesis is actually about, is the **DE members'
response to amplitude** — and `amp085` moves the wrong way for the hypothesis anyway: it **buys**
`mse` and **loses** `nmae` (+0.004654 raw, interval excludes zero).

---

## Gate 3 — the answer

### 3.1 The emission axis: `mse` and `pds` pinned, and nothing left over

| | overall | vs the +0.0019 shape-only gate |
|---|---:|---|
| best single emission point (`ts070`) | **+0.000131** | **7 %** |
| five-point emission ORACLE (each member's best, unreachable) | **+0.000442** | **23 %** |
| re-seeded twin, same pricing | +0.000096 | 5 % |
| at the class's measured h1→board `reach` factor 0.4277 | +0.000056 | 3 % |

The best emission point is **1.4× a re-seed** and its only moving member (`reach` +0.000803) has
a paired interval [−0.000012, +0.001903] that includes zero. **The emission axis is empty.**

`kt16000` — growing the call set — is the worst arm in the sweep at **−0.007338 overall**,
carried by `jac` −0.009270 (interval excludes zero) and `fid` −0.006368, for `reach` +0.002345
(interval includes zero). **094's "growing the call set further cannot pay" is confirmed on the
real scorer, six members at once, for the first time.**

### 3.2 The amplitude axis: the trade exists, at 35 : 1 against us

For `s` ∈ [1.00, 1.20] the members do move in opposite directions — `mse` +0.218270 raw against
`nmae` −0.002918 and `reach` +0.001696. **That is the trade, and it is priced:**

| | overall |
|---|---:|
| what `amp120` buys (`nmae` + `reach` + `fid` + `jac` + `pds`) | **+0.001340** |
| what `amp120` pays (`mse`) | **−0.038293** |
| ratio | **28.6 : 1 against** |
| `amp145`: buys +0.002374, pays −0.096787 | **40.8 : 1 against** |
| `amp200`: buys −0.008060 (`nmae` has turned), pays −0.266843 | catastrophic |

**And the ceiling, computed with `mse` counted FREE** — the most generous possible reading, the
one that grants the leader's premise that a team at the `mse` clamp pays nothing more
(`ceiling317.py`):

| | overall | vs the +0.0089 amplitude gate |
|---|---:|---|
| best single amplitude point with `mse` free (`amp145`) | **+0.002374** | **27 %** |
| five-point amplitude ORACLE with `mse` free (each member's best from a different `s`) | **+0.005190** | **58 %** |
| `nmae`'s entire contribution to that oracle | +0.000806 | 9 % |
| `reach`'s entire contribution | +0.001435 | 16 % |

**Even if `mse` were free, the whole amplitude axis cannot clear the gate an amplitude-changing
arm must clear.** And `mse` is not free: our live entry's `score_mse` is **0.18013**, so the
entire pot on the `mse` side is **+0.03002 overall** — every point of which the measured exchange
rate spends at a loss.

### 3.3 How far this is from the leader's shape

Not on the same axis at all. **Measured on the board itself** (`board_trade317.py`; 41 snapshots
→ **2,754 distinct submission ids**, 793 teams, **1,961 consecutive same-team pairs**):

**Tier A — `mse` readable and moved** (both ends raw `mse` ∈ [0.70, 1.05], |Δ| > 0.002), n = 285:

| a +1.0 raw `mse` GIVEN UP brings | median | mean [95 % over pairs] |
|---|---:|---|
| Δ raw `nmae` (lower is better → + is worse) | **+0.0812** | +0.0272 [−0.1884, +0.2120] |
| Δ raw `reach` | **−0.1912** | −0.2938 [−0.5704, −0.0155] |

- gave up `mse` (n = 112): Δ`score_avg` mean **−0.01475** [−0.02723, −0.00379], **33.9 % improved**
- bought `mse` (n = 173): Δ`score_avg` mean **+0.01538** [+0.00796, +0.02374], **83.2 % improved**

**Tier B — direction-stable** (|Δ raw `pds`| < 0.002 *and* |Δ raw `mse`| > 0.01), n = 75 pairs /
39 teams — `pds` is scale-invariant, so holding it nearly fixed removes most model rewrites:

| | median | mean [95 %] |
|---|---:|---|
| Δ raw `nmae` / Δ raw `mse` | **+0.0555** | −0.0184 [−0.1986, +0.1121] |
| Δ raw `reach` / Δ raw `mse` | **−0.0029** | −0.0529 [−0.1947, +0.0940] |
| **overall points per +1.0 raw `mse` given up** | **−0.17583** | −0.04064 [−0.18785, +0.12963] |

- gave up `mse` (n = 27): median **−0.00374**, **25.9 % improved**
- bought `mse` (n = 48): median **+0.00787**, **85.4 % improved**

**Tier C — exactly-pure calibration** (Δ raw `pds` = 0 to 1e-9): **n = 7 pairs, only 2 giving up
`mse`.** Underpowered; same signs (Δ`nmae`/Δ`mse` +0.0570, Δ`reach`/Δ`mse` −0.1389) and nothing
more should be read into it.

**And the leader is an outlier inside a bad band, not a team at a clever operating point.** 182 of
778 rows sit in their raw-`mse` band [0.95, 1.05]:

| score | raw `mse` | raw `nmae` | raw `reach` | raw `pds` | team |
|---:|---:|---:|---:|---:|---|
| **0.24018** | 0.9915 | **0.8717** | **0.3711** | 0.8683 | **ViRtUaL cElL builders** |
| 0.17602 | 1.0234 | 0.9150 | 0.2198 | 0.8140 | Zhenxing Fang |
| 0.16525 | 0.9822 | 0.9185 | 0.2105 | 0.8131 | xiaoyan chen |
| 0.16450 | 0.9993 | 0.9350 | 0.2190 | 0.8288 | emira |
| 0.15560 | 0.9540 | 0.9343 | 0.1857 | 0.8044 | Huiyan Qian |

The band's **median `score_avg` is −0.25606.** The leader beats the next-best row *in their own
`mse` band* by **+0.043 raw `nmae` and +0.151 raw `reach`**. Nobody reached that shape by
choosing an operating point. **They have DE-level accuracy that nobody else at that `mse` has;
that is a model, not a calibration.** Reading their row as a reachable operating point would be
§4.1c's retracted error again, in a new coordinate.

---

## 4. Is 094's saturation a property of the METRIC or of our CONSTRUCTION?

**Both — and they are different terms, which is why the mandate's suspicion is right.**

094's sentence is about **`fid`**, not `reach`. Its §6: *"Inverting the size step v15 → v11 with
the harness precision curve puts coverage at ≈ 0.93 at `n_pred` 2370 and ≈ 1.0 at 5340: **`fid`'s
coverage term is already saturated at v11**, so growing the call set further cannot pay."*

- **The saturated thing is a hard cap IN THE METRIC.** `fid` carries the yield factor
  `min(1, n_pred/n_real)`. Once our called set is at least the reference's, that factor is exactly
  1 and cannot rise, while a larger set lowers precision and enlarges `jac`'s union. **Our
  construction sits on the far side of that cap** (094: `n_pred` ≈ 5,961 on the challenge panel
  vs ≈ 5,340 on the harness), so it binds for us. **317's `kt16000` confirms it on the real
  scorer:** doubling the low-spread set to 16,000 genes gives `fid` −0.006368 and `jac` −0.009270
  for `reach` +0.002345 — **−0.007338 overall**.
- **`reach` is NOT saturated, in the metric or for us.** Fitted baseline `b = 0.079013`,
  replicate anchor `r = b + span = 0.964858` (AGENTS.md §4 quotes 0.958–0.978 — the fit lands
  inside). Across 778 board rows raw `reach` runs **0.04551 → 0.37115**, median 0.15157. **Ours,
  0.27428, is 22.0 % of the way from baseline to replicate; the leader's 0.37115 is 33.0 %.**
  Neither is near the ceiling.
- **What IS closed for `reach` is the ORDERING sub-channel, by arithmetic, not the level.**
  231 / §8.1c: the remaining prize for becoming the best `reach` team on Earth is **+0.002689
  overall** — 0.6 sd of the 0.004502 ordering noise — and its transfer factor is unidentified
  (four board pairs give 0.390 / 1.2927 / 1.4346 / 3.0563). Between teams, amplitude explains
  almost none of the field's `reach`: raw `reach` on `log(raw mse)` gives **R² 0.189 over the top
  60, 0.010 over the top 200, 0.004 over all 778** (`work/reach_amp317.txt`).

**So: 094's saturation is a metric cap on `fid`'s yield term that our construction has already
reached. It is not a statement about `reach`, and `reach` has ~3× its current value still on the
table in metric terms — reachable only through better signs on the reference's DE genes. That is
the direction programme, not the operating point.**

### 4.1 The bars this had to clear

Raw units of a single member, at transfer factor 1.0:

| gate | `pds` | `mse` | `nmae` | `fid` | `reach` | `jac` |
|---|---:|---:|---:|---:|---:|---:|
| shape-only **+0.0019** | +0.00515 | −0.01083 (−1.3 %) | −0.00688 (−0.7 %) | +0.00339 | +0.01010 (**+3.7 %**) | +0.00418 |
| amplitude-changing **+0.0089** | +0.02411 | −0.05073 (−6.2 %) | −0.03222 (−3.4 %) | +0.01588 | +0.04730 (**+17.2 %**) | +0.01959 |

---

## 5. What this re-establishes, and what it does not

**Re-established, on stronger evidence than it had.** §4.1c's *"there is ONE route, not two:
transfer cosine"* survives the live board's counter-example. The new evidence is of three kinds
the 2026-09-04 argument did not have: (i) **all six members measured on the real scorer at ten
operating points of our own pipeline**, each one knob off a byte-identical rebuild of a board
point; (ii) a **structural proof, gated in the build and confirmed by the scorer**, that the
whole emission axis leaves `pds` and `mse` at exactly zero, so only amplitude can make the trade
at all; (iii) the **board's own 1,961 same-team pairs**, which price surrendering `mse` at
−0.176 overall points per raw unit and show 74 % of the teams that tried it got worse.

**Consistent with the prior record, not in tension with it.** 187 (the `mse`+`nmae` joint
optimum at 0.982, `nmae` turnover at +29 %) is independently confirmed here at +20 %. 211 (`claim`
at cap). 227 (`nmae` exactly pinned on emission and call-set changes — reproduced as
`−0.000000 [0, 0]` on 58 targets). 231/§8.1c (the ordering channel closed by arithmetic).
034's `reach` +14–18 % is **not** contradicted: that was `spread` 0 → 0.35 on the *even-emitter*
generation; the shipped point already carries `spread` 0.25 / `--overdispersion 4.0`, and the
residual moves from here are +0.0005 to +0.0008 raw.

**What 317 does NOT establish.** It does not say the leader's row is unreachable — it says it is
not reachable *by moving our operating point*, and that their band is where weak whole-profile
predictions live, so their advantage is a DE-level modelling advantage. It does not locate the
board's amplitude optimum (h1-test cannot; §2.2). It does not price a `reach` arm for the board:
the h1-test → board `reach` factor for this arm class is 0.39–0.4277 from two routes, and 221
found no factor at all for two other classes. And every interval here is on **85 h1-test targets
with the H1 union off**, a surface whose absolute `mse` sits 0.20 raw above the board's.

**No slot is implicated. Nothing here is competition-facing, selected, or tuned.**

---

## 6. Cost and cleanliness

| | |
|---|---|
| instances | 11 × `r7i.2xlarge` on-demand, `us-east-1` |
| wall time | 164.1–180.3 min each |
| **total spend** | **$17.18** (cap was $25) |
| **per six-member score** | **$1.56** (221 measured $1.48) |
| termination | **11 of 11 terminated**, verified by `describe-instances` (`cost317.py`); the collector terminated each worker the moment its `agg_results.csv` landed in S3, before the ~4 h `prep-real-bundle` stage |
| local | peak RSS < 1 GB; **no local `cell-eval2`**; `uv run python` throughout |
| writes | only inside `experiments/317-operating-point/` |

One pre-existing `t4g.small` (`i-0a4311287d708e2ac`, running since 2026-09-04T19:00Z, ≈ $0.017/h)
is **not ours** and was left alone; flagging it for `scripts/aws_reap.sh`.

**Artifacts.** `build317.py` (the sweep builder, 221's four gates kept + gates 4/5 re-purposed),
`invariance317.py`, `amp_equiv317.py`, `price317.py`, `ka_map317.py`, `board_trade317.py`,
`reach_amp317.py`, `sweep317.py`, `paired317.py`, `ampcurve317.py`, `ceiling317.py`,
`cost317.py`, `monitor317.sh`, `aws/prologue-*.sh`; every number above is in `work/*.txt` /
`work/sweep317.json`, and the eleven `agg_results.csv` / `results.csv` in `cloud-output/`.

---

**317 (operating-point trade `mse` → `nmae`+`reach`): DOES NOT EXIST (peak composite +0.000131
overall on the emission axis and +0.002374 on the amplitude axis with `mse` counted FREE, against
0.2177 shipped — 7 % of the shape-only gate and 27 % of the amplitude gate; it pays `mse`
+0.2183…+1.5210 raw and buys `nmae` −0.0029 raw and `reach` +0.0076 raw, both with paired
intervals through zero, an exchange of 28–41 : 1 against us).**

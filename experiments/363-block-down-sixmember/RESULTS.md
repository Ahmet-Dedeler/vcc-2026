# 363 — Block-down is dead: fid loses (−0.0036, resolved) and pds pays 5× the bar (−0.027 overall)

Status: **NEGATIVE** (kill criterion fails on both legs: pds price exceeds the bar 5-fold;
fid moves the wrong way with its entire CI below zero).

## Known-answer reproduction (own code, `run363.py`)

`cap15` (unmodified `effects.npz::cap15`) vs
`experiments/323-current-family-validation/work/metrics-cap15.json`:

| check | result |
|---|---|
| emitted-count digest sha256 | **equal** (`547cf0ab…bc3a529f`) |
| panel-total digest sha256 | **equal** |
| raw pds / mse / fid / reach / jac | **exact** (0.0) |
| raw nmae | diff **1.1e-16** (1 ulp; per-target rows bit-exact, summation/serialization noise; gated ≤5e-16, recorded) |
| mse legs (num/den per target) | **exact** |
| stats arrays (targets/sums/jk/p/q/lfc) | **exact** (`array_equal`, NaN-aware) |
| per-target rows (510 rows) | **exact** |

Match. Changed arms are interpreted. Design-stage checks: reconstructed block =
top-400 |design| head + top-2200 `ref_lin` non-head per target (361's recipe bit-for-bit):
**2200/2200 called on all 85 targets** (pooled frac 1.0000); 2 universe genes absent
from the effects axis zeroed (`TBCE-1`, `HSPA14-1`); |.| preserved elementwise on every
changed arm; flipped fractions of carrying block entries — down 0.5006, up 0.4980;
realized block sign — down 0.9986 negative (rest exact zeros), up 1.0000 non-negative,
random 0.4979 negative (seed 20260910).

## What was measured

Four arms through 323's exact emitter + scorer from this directory (`run363.py`;
fresh `make_rngs(0)` per arm so library picks and noise draws are shared across arms):
`cap15`, `block_down` (design[block] = −|design|), `block_random` (iid 50/50 sign,
seed 20260910), `block_up` (+|design|). Call set, magnitudes, head untouched by
construction. Contrasts: paired 10,000-target bootstrap (seed 20260909), mse as
whole-panel ratio with claim cap recomputed per resample (identity ≤1e-10 passes on
all arms); H1-local scaling with 323's anchors. Plus per-target fid Δ vs own block
down-tilt, delivered-vs-design sign tracking by |design| decile, and fid contribution
by head/block/other region (`report363.py`).

## Results — six raw members per arm (means over 85 targets)

| arm | pds | mse | nmae | fid | reach | jac | H1-local overall |
|---|---|---:|---:|---:|---:|---:|---:|
| cap15 | 0.889636 | 0.974946 | 0.950836 | 0.497555 | 0.254347 | 0.074760 | 0.280953 |
| block_down | 0.816527 | 0.969708 | 0.954766 | 0.493920 | 0.260293 | 0.074767 | 0.252951 |
| block_random | 0.897619 | 0.989464 | 0.958852 | 0.491399 | 0.210425 | 0.074759 | 0.271837 |
| block_up | 0.856303 | 0.987784 | 0.960709 | 0.485312 | 0.201861 | 0.074704 | 0.252851 |

Paired Δ vs cap15 (point [95% CI]; bold = entire interval resolved):

| contrast | pds | mse | nmae | fid | reach | jac | overall |
|---|---:|---:|---:|---:|---:|---:|---:|
| down−cap15 | **−0.073109** [−0.1106,−0.0395] | −0.005238 [−0.0154,+0.0251] | +0.003929 [−0.0010,+0.0093] | **−0.003635** [−0.0073,−0.0000] | +0.005945 [−0.0245,+0.0363] | +0.000007 | **−0.028002** [−0.0469,−0.0122] |
| random−cap15 | +0.007983 [−0.0035,+0.0204] | **+0.014518** [+0.0097,+0.0203] | **+0.008016** [+0.0045,+0.0121] | **−0.006156** [−0.0095,−0.0029] | **−0.043922** [−0.0846,−0.0099] | −0.000001 | **−0.009116** [−0.0207,−0.0003] |
| up−cap15 | **−0.033333** [−0.0657,−0.0017] | **+0.012837** [+0.0011,+0.0499] | **+0.009872** [+0.0054,+0.0146] | **−0.012243** [−0.0168,−0.0079] | **−0.052486** [−0.1094,−0.0003] | −0.000056 | **−0.028102** [−0.0471,−0.0128] |

(mse/nmae positive = worse. jac never moves: call sets ~identical, jaccard vs cap15
0.99990 / 1.00000 / 0.99693; n_pred mean 2599.5 / 2599.3 / 2599.5 / 2593.3.)

## Verdict against the PLAN.md kill criterion

**Not board-worthy — dead on the pds leg alone, and the fid leg points the wrong way.**

- pds price of block_down: raw −0.073109 → H1-scaled −0.073109/0.44366947/6 =
  **−0.02746 overall**, against the bar of 0.0055 (half the projected A-side fid
  gain). **Exceeds the bar 5.0×.** mse shows no price (−0.005238, n.s.) — the entire
  cost is pds.
- fid: −0.003635 with CI entirely negative, against 361's exact-arithmetic prediction
  of +0.0081. The tilt gain does not survive emission (mechanism below).
- Random-vs-coherent test (PLAN's discriminator): block_random costs pds **nothing**
  (+0.0080, n.s.) while coherent down costs −0.073 and coherent up −0.033 (both
  resolved). **The price is the common component, not re-randomisation** — 155's
  uniform-shift lesson holds for sign-coherent blocks.
- Tilt extrapolation for the PI: per-target fid Δ vs own block down-tilt, slope
  +0.00185 per 0.01 tilt (r = +0.37, n = 85; tilt mean 0.5178, range 0.444–0.640).
  At Jurkat's 0.540 (+0.0222 tilt) the projected A-side fid Δ ≈ −0.0036 + 0.0041 ≈
  **+0.0005 raw** — i.e. ~zero even where the tilt is 4× H1's — while the pds price
  is tilt-independent common-mode. No context salvages this.

## Why 361's +0.0081 does not survive emission (measured mechanism, two halves)

361 flipped *delivered* lfc signs with q/p fixed. Emission re-runs DE from counts,
and two things break:

1. **Transmission is magnitude-gated; small effects are sticky.** Delivered-vs-arm-design
   sign agreement by |design| decile: cap15 and block_random transmit identically
   (0.62 → 1.00 across deciles — the emitter is a faithful, sign-blind transmitter),
   but block_down agrees at only 0.055/0.137/0.287 in the bottom three deciles: with
   shared picks/noise draws across arms, small-effect delivered signs are pinned by the
   noise realization, not the design sign. Forcing signs moves delivery only in the top
   half of the block.
2. **Collateral damage on untouched genes exceeds the tilt gain.** The whole fid loss sits
   on the block (head/other contributions identical to 4 decimals). On *flipped*
   (positive→down) genes delivered-vs-truth accuracy moves 0.4924 → 0.4953 (+0.003 on
   93.6k genes — the tilt gain is real but tiny, since transmission is partial). On
   *kept-negative* genes, whose design is bit-identical, accuracy FALLS 0.5245 → 0.5114
   (−0.013 on 93.1k genes; delivered signs agree across arms only 72%; non-block genes
   83%). The coherent down-shift of 2,200 high-CPM genes re-routes library mass (CPM
   coupling) and scrambles delivery of everything else. Net block contribution
   0.41751 → 0.41399 = −0.0035.
3. **Profile vs call dissociation.** pds/mse read continuous profiles (the coherent shift
   moves all 85 profiles → retrieval collapse) while fid/reach/nmae read thresholded
   calls (call set jaccard ~1.0). One intervention, opposite readings — the same reason
   361's delivered-lfc counterfactual could not see the pds cost.

Mirror note: block_up is worse everywhere than block_down (fid −0.0122 vs −0.0036),
asymmetric as the 51.8% down-tilt demands — the tilt is real, it is just not collectable
this way.

## Limits

H1 development panel (85 targets, one held-out line); not an A/B/C forecast — the Jurkat
sentence above is a linear extrapolation, not a measurement. Bootstrap CIs are
target-resampling intervals on H1, not transfer intervals. nmae cohort 58 (27 ineligible),
reach cohort 82 (3 NaN) — per-member valid-cohort means (report342 convention). One
ulp (1.1e-16) nmae known-answer tolerance, everything else bit-exact. No claim about which
mechanism the leader uses.

## Exact commands to reproduce

```
cd ~/Code/virtual-cell-challenge-2026
memory_pressure | tail -1   # proceed only if free >= 50%
uv run python experiments/363-block-down-sixmember/run363.py --arms cap15            # ~149 s, known-answer gate
uv run python experiments/363-block-down-sixmember/run363.py --arms block_down,block_random,block_up  # ~423 s
uv run python experiments/363-block-down-sixmember/report363.py                      # bootstrap + tilt + tracking
```

Outputs in `experiments/363-block-down-sixmember/work/`: `design.json` (flip/block
diagnostics), `blocks.npz` (partition), `stats-<arm>.npz`, `per-target-<arm>.csv`,
`metrics-<arm>.json`, `known-answer.json`, `results.json`, `per-target-363.csv`
(85 rows: tilt, n_pred/n_real, per-arm fid and Δ), `execution.json`.

## Peak memory and wall time

Peak RSS 5.44 GB (ru_maxrss bytes; cap15 arm 149 s; three changed arms 423 s total;
report363.py seconds-scale on saved arrays);
sequential arms, one process, NUMBA/OMP/POLARS/OPENBLAS threads 4. System free 72–75%
throughout; 12 GiB cap never approached. No downloads, no reserve reads, no submissions,
no files written outside `experiments/363-block-down-sixmember/`.

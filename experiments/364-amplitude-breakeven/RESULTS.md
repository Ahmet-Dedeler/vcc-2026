# 364 — Uniform amplitude never reopens: no break-even at any sign accuracy, including oracle

Status line: **COMPLETE** (kill criterion met — the route is dead, and the oracle
endpoint proves it is structurally dead, not just dead at current accuracy).

## Known-answer reproduction (own code, `breakeven364.py`)

| quantity | reproduced | recorded (`metrics-cap15.json` / `per-target-cap15.csv`) | match |
|---|---|---:|---|
| nmae mean (58 eligible) | 0.95083646 | 0.95083646 | max per-target diff 4.4e-16 |
| fid mean (85) | 0.49755452 | 0.49755452 | max per-target diff 0.0 |
| reach mean (82 scored, 3 NaN) | 0.25434726 | 0.25434726 | max diff 0.0, NaN set identical |
| eligible targets | 58 | 58 (343) | yes |

Conventions: gate = {ref q<0.05} ∩ finite-nonzero ref lfc ∩ not-own-target;
eligible = gate ≥ 10; reach ordering sig-first/q/p/|lfc|/gene (361); means are
nanmeans over the valid cohort per member (report342 convention).

## What was measured (inputs, estimator, controls)

Inputs (read-only): 323 `stats-cap15.npz`, `reference.npz`,
`per-target-cap15.csv`, `h1-anchors.json`; 317 `sweep317.json`.
Delivered gate accuracy: **pooled |l|-weighted 0.5619, mean-over-targets 0.5880**
(== 343), magnitude ratio 0.2065, 64,931 eligible gate pairs, 29,345 wrong-sign.

Estimator: for each pooled-accuracy target p ∈ {0.60,…,0.75}, 20 seeded draws —
shuffle all wrong-sign gate pairs, flip the prefix (sign only; |lfc|, q, p
untouched) until pooled accuracy ≥ p (achieved p = target +2e-5 every draw).
At each (draw, p, a ∈ {1, 1.2, 1.45, 2}): L = a·s_new·|lh| on gate genes,
off-gate unchanged; exact recompute of nmae/fid/reach. pds/mse from 317's
recorded raw deltas (table below), assumed p-independent. H1-local scaled Δ via
`h1-anchors.json` spans; jac carried at 0 (q untouched ⇒ call sets fixed);
overall(p,a) = mean of six scaled deltas vs delivered (p0, a=1).

Controls: (i) known-answer gate above — PASS or the run stops; (ii) fid/reach
invariance across a verified at **0.0e+00** (fid provably: sign(a·x)=sign(x),
a>0; reach empirically: the |lfc| tie-break never fires on continuous q/p);
(iii) delivered-sign amplitude leg vs 317's independent measurement —
×1.2: ours −0.00319 raw vs 317's −0.00292; ×1.45: −0.00224 vs −0.00220;
×2.0: +0.0150 vs +0.0363. Same sign throughout, first two within 9% across
different surfaces/generations; (iv) oracle anchor nmae(a=1) = 0.80339 ==
343's "all signs right" 0.803.

## 317 numbers used (recorded, `sweep317.json` arms, base = rc_callset rebuild)

| arm | Δpds raw | Δmse raw |
|---|---:|---:|
| amp120 (×1.2) | +1.1e-16 (≈0) | **+0.21827016** |
| amp145 (×1.45) | −0.00056022 | **+0.55168990** |
| amp200 (×2.0) | −0.00126050 | **+1.52102760** |

## Results table 1 — raw nmae(p, a), mean over 58 eligible targets (± sd over 20 draws)

| p (pooled / mean-over-targets achieved) | ×1 | ×1.2 | ×1.45 | ×2.0 |
|---|---:|---:|---:|---:|
| delivered (0.5619 / 0.5880) | 0.95084 | 0.94765 | 0.94860 | 0.96581 |
| 0.60 (0.6000 / 0.6238±0.0023) | 0.93810±0.00085 | 0.93254±0.00101 | 0.93070±0.00121 | 0.94250±0.00157 |
| 0.62 (0.6200 / 0.6421±0.0030) | 0.93138±0.00113 | 0.92458±0.00136 | 0.92127±0.00163 | 0.93026±0.00210 |
| 0.65 (0.6500 / 0.6706±0.0037) | 0.92119±0.00116 | 0.91249±0.00140 | 0.90694±0.00169 | 0.91156±0.00226 |
| 0.70 (0.7000 / 0.7172±0.0035) | 0.90465±0.00123 | 0.89286±0.00148 | 0.88365±0.00177 | 0.88115±0.00228 |
| 0.75 (0.7500 / 0.7644±0.0054) | 0.88779±0.00176 | 0.87286±0.00210 | 0.85996±0.00252 | 0.85028±0.00333 |
| **oracle 1.00** (deterministic) | **0.80339** | 0.77275 | 0.74137 | 0.69595 |

fid(p): 0.50123 / 0.50319 / 0.50613 / 0.51100 / 0.51590 / oracle 0.54039.
reach(p): 0.25828 / 0.26059 / 0.26444 / 0.27640 / 0.29410 / oracle 1.000
(reach=1.0 at oracle is arithmetic: every prefix 100%-pure ⇒ deepest ≥90%
prefix is the whole pool — vacuous, reported for completeness, not priced as a forecast).

## Results table 2 — overall(p, a), H1-local points vs delivered (six-member mean)

| p | ×1 (sign gain alone) | ×1.2 | ×1.45 | ×2.0 |
|---|---:|---:|---:|---:|
| 0.60 | +0.00517 | −0.03105 | −0.08849 | −0.26008 |
| 0.62 | +0.00795 | −0.02791 | −0.08494 | −0.25571 |
| 0.65 | +0.01223 | −0.02309 | −0.07947 | −0.24899 |
| 0.70 | +0.02030 | −0.01413 | −0.06947 | −0.23695 |
| 0.75 | +0.02958 | −0.00395 | −0.05822 | −0.22364 |
| oracle | +0.19634 | +0.16727 | — | — |

Member legs of the marginal ×1.2 step (×1.2 vs ×1 at same p), scaled overall points:

| p | nmae buys | mse pays | fid/reach/pds/jac | net |
|---|---:|---:|---|---:|
| 0.60 | +0.00957 | −0.22689 | 0 | −0.03622 |
| 0.65 | +0.01499 | −0.22689 | 0 | −0.03532 |
| 0.75 | +0.02572 | −0.22689 | 0 | −0.03353 |
| oracle | +0.05275 | −0.22689 (assumed) | 0 | −0.02902 |

## The break-even line

- **a*(p)** (amplitude the sign gain at p can carry while staying ≥ 0 overall,
  interpolated): 1.029 / 1.044 / 1.069 / 1.118 / **1.176** at p = 0.60…0.75.
  Even at p=0.75 only ×1.18 rides free — and that free ride is the sign gain
  paying for the amplitude, not amplitude paying for itself.
- **p*(a) does not exist for any a>1 at any p ≤ 0.75** — no zero crossing.
  Linear extrapolation of the last segment: p*(×1.2) ≈ 0.77, p*(×1.45) ≈ 1.01
  (impossible), p*(×2.0) ≈ 1.6 (impossible).
- **Oracle kills the extrapolation**: at perfect signs the marginal ×1.2 leg is
  still **−0.0290** (nmae +0.0528 scaled vs mse −0.2269). For it to break even,
  the true mse cost would need to shrink **4.3×** (0.218 → 0.051 raw) through
  better alignment — against 317's board-validated 28–41:1 exchange the other way.
- **mse-free ceiling** (most generous reading, 317 §3.2 style): marginal ×1.2 at
  oracle with mse completely free = +0.0528/6 = **+0.00879 vs the +0.0089
  amplitude gate — fails by a hair with every advantage granted**
  (perfect signs AND free mse). Nothing weaker can pass.
- **Pay per +0.01 pooled accuracy** (OLS over the 5 p points): +0.00162 overall
  at ×1 (+0.00180 / +0.00202 / +0.00243 at ×1.2/1.45/2.0 — amplitude steepens the
  slope but never pays its own intercept).

## Verdict against the PLAN.md kill criterion

**Killed, at both readings.** Literal (no p ≤ 0.75 makes any a>1 net-positive):
met — ×1.2 is −0.0039 even at p=0.75. Functional (only non-uniform amplitude
could matter): the oracle shows even p=1.0 leaves the uniform ×1.2 leg at
−0.029, and the mse-free-ceiling + perfect-signs case still misses the gate.
**Uniform amplitude is structurally closed, not accuracy-gated.** Mechanism:
mse taxes the whole 18,533-gene profile while nmae collects only on the
reference gate (magnitude ratio there is 0.21 even delivered) — a whole-profile
tax for a gate-local prize. This is 317/227/347's verdict with the accuracy
dimension mapped out and the endpoint measured instead of extrapolated.

## What the running sign levers would need (plain statement)

At ×1 (no amplitude), each +0.01 of pooled gate accuracy is worth ≈ +0.0016
H1-local overall (nmae ≈ 70%, fid ≈ 15%, reach ≈ 15% of the gain at p=0.65).
346's same_ridge β=1 (+0.035 design-stage mean-over-targets, 0.600→0.635) would,
if it transferred 1:1 to delivered pooled accuracy (it will not — design≠
delivered, and it spends the recipient's own Pan-H1 labels, unavailable for
A/B/C), land near pooled ≈ 0.60: **≈ +0.005 overall at ×1, while ×1.2 on top
still nets ≈ −0.036 marginal. Take any sign gain at ×1; amplitude stays off
regardless.** 362/363 should be judged on their sign gain at ×1 only — no
compound with amplitude exists to be collected.

## Assumptions (as mandated, all conservative toward amplitude except the last)

1. pds/mse deltas from 317's h1-test surface applied as raw deltas on the 323
   surface; independence of p assumed. Direction of error: flipping wrong→right
   signs strictly reduces ‖pred−truth‖² on the gate, so the true mse cost at
   high p is SMALLER than assumed — the assumed leg overcharges amplitude at
   high p. Bounded above: even zero mse cost fails the gate at oracle (+0.00879
   vs +0.0089), so this assumption cannot reverse the verdict.
2. pds improvement from better signs NOT modeled (sign flips change the profile
   direction) — understates the sign gain at ×1; overall(p,a) is a lower bound.
   Does not help amplitude's marginal leg (pds is scale-invariant; 317 measures
   |Δpds| ≤ 0.0013 over the whole 2× span).
3. 317's mse deltas embed that surface's jackknife claim; the 323 surface's
   claim differs. No usable mse-vs-accuracy identity found in 227/188 for an
   H1-local adjustment (227's 0.010733-per-1% is claim-sensitivity, not
   accuracy-dependence) — hence the bounding argument in (1) instead.
4. Random-subset flips are an optimistic sign model (a real key flips both ways;
   234: real keys realise 0.44× the oracle rate) — true p from a lever costs
   more than shown here.

## Limits

H1 development panel (58 eligible targets, one held-out line); not an A/B/C
forecast. H1-local overall points (323 anchors), not board points. Sign flips
are expressive counterfactuals, not an implementable key. The p grid is pooled
accuracy; mean-over-targets runs ≈ +0.024 above it on this panel (table 1) —
quote which one when comparing with 346 (mean-over-targets, design-stage).

## Exact commands to reproduce

```
cd ~/Code/virtual-cell-challenge-2026
memory_pressure | tail -1   # proceed only if free >= 50%
uv run python experiments/364-amplitude-breakeven/breakeven364.py  # ~90 s, writes work/breakeven.json
uv run python experiments/364-amplitude-breakeven/oracle364.py     # ~5 s, writes work/oracle364.json
uv run python experiments/364-amplitude-breakeven/price364.py      # ~1 s, writes work/priced364.json
```

## Peak memory and wall time

breakeven364.py 87 s, oracle364.py 2 s, price364.py < 1 s; peak RSS ≈ 195 MB
(all arrays ≤ 85×10,779 float64; reach loop dominates). No process exceeded
4 GiB; no scratch files kept. work/ holds `breakeven.json`, `oracle364.json`,
`priced364.json` (+ scripts and PLAN.md alongside).

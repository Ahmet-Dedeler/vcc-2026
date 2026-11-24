# 300 — The one live operator thread, tested as a MECHANISM, not swept as a knob: does the whitener's rank-retrieval gain come from the SPREAD of its training rows over many directions? Three preregistered arms at the incumbent's exact row count — the reversed-reliability tilt inside the shape-only band, the essential-gene prune's interaction with it, and a direct spread manipulation — with mirror controls, leave-top-3-out retention, the coherent-pessimistic price, and the Flex reference WITH an interval. ≤ 3.5 hours. No build, no submission.

You are a research subagent. Work in `experiments/300-spread-tilt/`. **Edit nothing outside it.** No build, no `vcc`,
no `cell-eval2`, nothing under `submissions/`, `data/harness/` or `data/`. `uv run python`; at most two fold
processes in flight; peak RSS under ~11 GiB per process; `memory_pressure` ≥ 40 % free before each. **Time box:
3.5 hours.** `docs/READ-FIRST.md` §3 binds you (one variable per arm; every arm priced on the coherent-pessimistic
rule; known-answer rows first; the gate; the ordering-noise SD 0.0045 is three times 299's effect — say so wherever
an interval is quoted). Reserves are spent and not needed. This is NOT a knob sweep around `cd4svd5`'s parameters:
every arm below tests one stated mechanism at the incumbent's row count, and the doses are fixed before any number.

## Why (R299 §5–§8; R272 §0.10; R286 §6)
299 found that fitting the whitener on reliability-SELECTED rows loses (below every control; a row-count effect below
full rank, null at full rank), and — post hoc, at the incumbent's exact row count and donor proportions — that a
rank-preserving tilt of the training rows AWAY from the reproducible ones gains `pds` on 3/3 folds (+0.00094
[+0.00017, +0.00173] shape-only; +0.00139 at a steeper, amplitude-changing dose), passes its mirror control, and
reproduces in sign on 273's Flex instrument (+0.0014). It priced 0.74 × the shape-only gate with retention 47.8 %:
two narrow misses. The reading offered — the metric is paid by the donor rows' SPREAD over many directions, and
reproducible rows concentrate the covariance — is a mechanism, and it makes three testable predictions: (P1) the
gain is monotone in the tilt dose inside the shape-only band and saturates before the amplitude leaves it; (P2) it
interacts with 062's essential-gene prune (which already removes the MOST reproducible rows and is itself worth
+0.0021 on Flex): tilting on top of the prune gains less than tilting on the un-pruned banks, or nothing; (P3) a
DIRECT spread manipulation — weighting rows by the inverse of their alignment with the bank's leading directions
(low-leverage rows up-weighted), with reliability untouched — gains at least as much as the reliability tilt, and
the two are collinear (their per-target Δ`pds` correlate). If P3 holds, the reliability axis was a proxy; if the
tilt gains and P3 does not, reliability itself matters. **Nothing here is a candidate unless it clears the frozen
gate; the point is to learn the mechanism with the same machinery that would price a candidate.**

## READ FIRST
`experiments/299-reliable-training-rows/{RESULTS.md, PREREGISTRATION.md, rel299.py, run299.py, price299.py,
mech299.py, flex299.py}` and its `work/rel/<line>.npz` (per-row reliabilities — reuse; do not recompute) and
`work/sweep299-<fold>.json` (the `wgt:*` rows to reproduce as KA); 272's `op272.py` (`fit_operator` with row
weights; the noise spectrum); 263's `lib263.py`; 261's `price261.py`; 265's `boot265.py`; 273's `lib273.py` /
`instrument273.py`; `scripts/build_union_whitener.py` (062's prune + 12,000-row cap: reproduce the production
operator's row set as closely as its seed allows and state the residual difference).

## Protocol (write `PREREGISTRATION.md`, hash into `work/PREREGISTRATION.sha256`, before any new fold number)
1. **Known answers:** 299's `wgt:mild`, `wgt:mildrev`, `wgt:steeprev` per-fold Δ`pds` reproduced at ≤ 1e-9 through
   this experiment's plumbing; `cd4svd5` == 263 at 0.0e+00; the BR249 retrodiction +0.00472 at 0.0e+00; 273's
   INCUMBENT on the Flex frame at 0.0e+00.
2. **Arms (frozen doses, all at the incumbent's row count and donor proportions; every panel row kept):**
   (A) the reversed tilt at five doses w = 1 + d·(0.5 − u), d ∈ {0.5, 1.0, 1.5, 2.0, 2.5} (d = 1 is 299's `mildrev`,
   d = 1.8 its `steeprev`), read per fold with delivered amplitude, so the shape-only band (±0.5 %) is located on
   the dose axis; (B) the same five doses applied ON TOP of the essential-gene prune (062's row set), and the prune
   alone as its own KA against the production operator; (C) the direct spread manipulation: w = 1 + d·(0.5 − v)
   where v is the row's rank quantile of its energy share in the bank's top-50 PCs (fitted on the un-weighted rows),
   at d ∈ {1.0, 2.0}; (D) controls: the mirror of every A / C dose (w reflected about 1), and a random-weight
   control with the same weight distribution assigned to rows at random (2 seeds). Read every arm on the three
   cancer folds: per-target `pds`, own cosine, `cos_agg`, amplitude, 262's groups, knee crossings, the `nmae`
   identity leg, paired target bootstrap (20,000), 3/3 rule, leave-top-3-out retention; price with 261's recipe.
3. **Gate (frozen):** pessimistic ≥ +0.0019 shape-only / +0.0089 otherwise; 3/3 with the pooled interval above 0;
   retention ≥ 50 %; the arm above its mirror AND its random-weight control by more than half its Δ`pds`. State
   P1–P3's outcomes with the rules: P1 = the A-dose curve monotone and the largest shape-only dose's Δ`pds` above
   d = 1's with the paired interval above 0; P2 = (B at the best A dose) − (A at that dose) with its interval; P3 =
   Spearman(per-target Δ`pds` of C, of A at matched dose) and C's own gate reading.
4. **The Flex reference (273's union-OFF instrument) WITH an interval:** the best shape-only A dose and the best C
   dose vs their matched refit, paired target bootstrap over the 85 targets (the refit is one eigendecomposition per
   arm; budget two).
5. **Mechanism row:** 272's statistic (share of directions with noise > eigenvalue; whitened energy on them) for the
   best A dose and the best C dose, next to the incumbent's 76.8 % / 86.8 % and 299's `sel` 99.9 %.

## Deliverable
`PREREGISTRATION.md` (hashed), `RESULTS.md`: KA table; the A / B / C × fold table with controls, price, gate per
arm; P1–P3; the Flex readings with intervals; the mechanism row; `work/pertarget300-<fold>.npz`; one line:
**SPREAD MECHANISM: P1 … / P2 … / P3 …; best shape-only arm … Δ`pds` … [lo, hi], 3/3 …, retention …, controls …,
priced … vs +0.0019 → CLEARS / DOES NOT CLEAR; FLEX … [lo, hi]**. If an arm clears every condition, say so plainly
and stop — the coordinator decides about a build; do not build. Then stop.

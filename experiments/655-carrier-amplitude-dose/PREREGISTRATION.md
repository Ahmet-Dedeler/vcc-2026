# 655 — DE carrier amplitude dose curve under production's zero repair: ×1.5 and ×2.0 (PREREGISTRATION)

Council session [7772dc]; the coordinator's range 654–655 (low priority). Written and committed BEFORE any 655 run, after
654's seed 0 on the production-faithful construction (Amendment A1). Class: **scorer mechanics** (the DE carriers' design
amplitude; the bulk is re-levelled in production, ballast-pinned here).

## Why
654 seed 0, with `--zero-repair -5`: carriers ×1.25 vs ×1, overall-eq:
- f 1 +0.0013, f 0.5 +0.0015, f 0.35 +0.0011, cmp +0.0012;
- nmae −0.0042 raw at f 0.35; fid / jac / reach unchanged;
- head DOWN delivered |lfc| +22.5 %.

631 put the lfc-space optimum for called DOWN genes at m 2. So the dose curve's peak may lie above ×1.25.

## Condition to run
- 655 runs only if 654's production-faithful ×1.25 − ×1 is > 0 at f 0.35 on at least 2 of 3 seeds.
- If not, 655 is withdrawn unrun, and this file records that.

## Arms (one variable: `eff_de_native` × s, bulk = `ps_bulk(eff)` unchanged; all with `--zero-repair -5`)
- `x150` (s 1.5) and `x200` (s 2.0), seeds 0, 1, 2. They are new design branches in `run655_496.py`, a copy of 654's
  `run654_496.py` that adds exactly these two branches (diff recorded).
- References (no new runs): 654's base (= base_zr, identical by A1's KA) and 654's x125_zr, seeds 0–2.

## Known answer
`run655.py cpm x125 0 … --zero-repair -5` == 654's `stats-x125_s0_bps_zr-5_prodem_…_dnk250zs.npz` bit for bit (q, p, lfc,
sums). This exercises the copy's scaling path end to end.

## Measures
Exactly 654's `analyze654.py --zr` measures, run on 655's arms against 654's base:
- members at 521's four regimes, overall-eq;
- 631's regions (share and delivered |lfc|);
- calls per target;
- zero-total counts (`zerocount654.py`).

## Predictions (registered now)
- **D1:** the curve is concave with an interior peak. ×1.5 − ×1 overall-eq at f 0.35 is > ×1.25 − ×1 (the peak is above
  1.25). ×2.0 − ×1 is < ×1.5 − ×1 (past the peak).
- **D2:** fid, jac and reach stay within ±0.001 / ±0.0005 / ±0.01 raw of base through ×1.5. Calls per target rise
  monotonically: 2,571 → 2,597 at ×1.25.

## Bins (overall-eq vs ×1 at f 0.35 AND f 0.5, every seed the same sign; the best dose s* of {1.25, 1.5, 2.0})
- **s* reads ≥ +0.0015 at both regimes:** a candidate October knob through `VCC_ALPHA` (α × s*). Its re-level absorbs the
  bulk. It needs, in order:
  - a builder KA showing the α path moves the carriers exactly as `eff_de` × s does here;
  - a held-out read on a cancer surface with zero repair (615's Orion rows or 614's lever rows, once trusted);
  - the A/B/C board read with its own prereg.
  Not October before those.
- **+0.0005 … +0.0015:** a weak candidate. Record it; the same steps apply only if 615 / 614 agree.
- **< +0.0005:** amplitude is a small lever. Record it and close the family (654 + 655).

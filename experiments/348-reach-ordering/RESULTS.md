# 348-reach-ordering — ordering lever CLOSED (best design-time key +0.0264 < +0.03 bar)

Status: **NEGATIVE** (clean negative with controls; kill criterion fired).

## Known-answer reproduction

Independent polars walk of BRIEF §5
(pool = ref-sig genes minus own target gene; rank = sig-first, p_adj↑, p↑,
|lfc|↓, gene↑; adjudicable-only depth; k\* = deepest k with P(k) ≥ 0.9 else 0;
R = k\*/N_conf; NaN iff pool empty) vs
`work/per-target-cap15.csv`, plus the installed scorer
(`prepare_de` + `de_direction_reach`, universe='adjudicated',
corrected=False) through the exact `evaluate323.score()` frame construction:

| check | result |
|---|---|
| scorer vs recorded, max abs diff (85 targets) | **0.0** |
| my walk vs recorded, max abs diff | **0.0** |
| targets off by > 1e-9 | **0 / 85** |
| panel mean (nanmean, 82 scored) | **0.25434726355678233** (exact) |
| NaN targets (empty pool, identical in all three) | PCCB, TRAFD1, TSKU |
| truth cross-check (parquet vs reference.npz) | max\|dq\| = 0, max\|dlfc\| = 0 |
| unadjudicable genes inside any pool | 0 |

`work/known-answer.json`. Proceeded only after this passed.

## What was measured

Same delivered lfc SIGNS throughout; only the rank key over the
reference-significant pool changes. Design matrix from `effects.npz::cap15`
mapped by gene name (universe coverage 10,777/10,779; the 2 missing genes get
|design| = 0, i.e. rank last — negligible). Common vector = mean design
effect over the 85 targets. Random = 100 seeds (`default_rng(3480000+seed)`,
full-axis permutation, pool-masked walk). Panel mean = nanmean over the same
82 scored targets in every scheme (NaN set verified identical). Paired-target
bootstrap (10,000 resamples) for CIs in follow-up.

## Results (units: raw reach, panel mean over 82 H1 targets)

| scheme | mean | Δ vs current | paired 95% CI of Δ | pre-reg? |
|---|---|---:|---:|---|
| O1 current (sig, p_adj, p, \|lfc\|, gene) | 0.25435 | — | — | yes |
| O2 random (floor) | 0.19869 ± 0.01201 (range 0.171–0.227) | −0.0557 | — | yes |
| O3 \|delivered lfc\| desc | 0.26815 | +0.0138 | [−0.0092, +0.0377] | yes |
| O4 \|design effect\| desc | 0.28075 | **+0.0264** | [+0.0059, +0.0514] | yes |
| O5 agree-with-common first, then \|design\| | 0.25244 | −0.0019 | [−0.0275, +0.0240] | yes |
| O6 oracle correct-sign first (ceiling) | 0.67478 | +0.4204 | [+0.3696, +0.4686] | yes |
| O7a predicted-down first, then \|lfc\| | 0.27195 | +0.0176 | [−0.0133, +0.0496] | yes |
| O7b disagree-first (reverse of O5) | 0.23073 | −0.0236 | [−0.0558, +0.0053] | yes |
| O7c sig-first then \|lfc\| (no p-refinement) | 0.24713 | −0.0072 | [−0.0216, +0.0060] | yes |
| O7d sig-first then \|design\| — POST-HOC | 0.25486 | +0.0005 | — | **no** |
| O7e down-first then \|design\| — POST-HOC | 0.27883 | +0.0245 | — | **no** |

Per-target sign of the O4−O1 delta: positive on 35%, negative on 15%, zero on
the rest. Anatomy: the gain lives in small-pool targets —
mean Δ = +0.0674 for n_conf ≤ 20 (n = 32) vs +0.0002 for n_conf > 20 (n = 50);
corr(n_conf, Δ) = −0.16. Largest single-target moves are ±1 pool step on
pools of 6–9 genes (ZNF32 +0.167 at n_conf = 6; COASY −0.222 at n_conf = 9).

### P(k) profile, current ordering (mean over targets with depth ≥ k)

| k | n targets | mean P(k) | frac targets with P(k) ≥ 0.9 |
|---:|---:|---:|---:|
| 5 | 70 | 0.806 | 0.457 |
| 10 | 58 | 0.762 | 0.328 |
| 20 | 50 | 0.711 | 0.180 |
| 50 | 47 | 0.673 | 0.128 |
| 100 | 41 | 0.647 | 0.073 |
| 200 | 38 | 0.625 | 0.026 |

Purity is already below 0.9 at k = 5. k\* = 0 on 7.3% of scored targets
(6/82: AASS, LBR, SET, SOX4, TOR1A, TSFM) — every one is a first-position
miss with **no** recovery at any deeper prefix. k\* < 10 on 74.4%.

### Mechanism notes (read off the pattern, not new runs)

1. Pure size orderings beat emitted-significance orderings: O4 > O3 > O1, and
   bolting the significant-first block back onto |design| destroys the whole
   gain (O7d − O4 = −0.0259, back to 0.2549 ≈ O1). The emitter's p-value
   prefix is a *worse* ranker than raw effect size; the within-block
   p-refinement helps only slightly (O1 − O7c = +0.0072).
2. The common-vector agreement key is dead: O5 ≈ 0 (−0.0019) while its exact
   reverse O7b is −0.0236 — the key carries a weak signal in the *agree*
   direction, worth nothing for depth.
3. The 343-motivated down-first split adds nothing once size is used
   (O7e − O4 = −0.0019); with |lfc| it is +0.0176 but its CI crosses zero.
4. Oracle +0.42 vs best design-time +0.026: reach is bound by per-gene sign
   information (343's verdict), not by ranking — consistent with 232/234.

## Verdict against the PLAN.md kill criterion

**FIRED → the ordering lever is closed.** Best pre-registered design-time
ordering (O4, +0.0264) misses the ≥ +0.03 bar; the two post-hoc combos do not
exceed it either (best +0.0245, and selected post-hoc so discounted anyway).
O4's CI (+0.0059…+0.0514) overlaps the bar, but the bar is on the point
estimate by design, and the gain concentrates in n_conf ≤ 20 targets where
one pool step is 0.05–0.17 — the coarsest, least generalizable stratum.

The conditional emitter analysis is therefore **not triggered**; for the
record, the analytic expectation (no new run): a pure within-call-set
re-ranking would leave `fid` and `jac` exactly invariant (both read only the
significant set + signs); O4-style re-ranking crosses the significance
boundary (O7d proves the block matters), so it would change the set and its
fid/jac side effects would need a full 323 scoring run to price. Not worth
one: the reach leg itself fails the gate.

## Limits

- H1 development panel (85 targets), not A/B/C; absolute reach values do not
  transfer, only the ordering-class conclusion.
- O4's edge is small-pool-driven; the A/B/C panel has a designed low-n_real
  stratum (§2.1), so the direction could plausibly differ there — but that
  cuts both ways and does not clear any gate.
- Post-hoc O7d/O7e were chosen after seeing O4; reported for mechanism only.
- Random floor (0.199) uses full-axis permutations; pool-restricted shuffles
  would read marginally higher — immaterial to the verdict.

## Exact commands to reproduce

```bash
uv run python experiments/348-reach-ordering/reach348.py known    # bit-exact KA check
uv run python experiments/348-reach-ordering/reach348.py order    # Table 1 + per-target CSV
uv run python experiments/348-reach-ordering/reach348.py profile  # P(k), k* diagnostics
uv run python experiments/348-reach-ordering/followup348.py       # bootstrap CIs + anatomy
uv run python experiments/348-reach-ordering/explor348.py         # post-hoc O7d/O7e
```

Outputs: `work/known-answer.json`, `work/orderings.json`,
`work/per-target-orderings.csv`, `work/profile.json`, `work/delta-ci.json`.
Inputs read (never written): `experiments/323-current-family-validation/work/`
`{stats-cap15,reference,reference_de.parquet,per-target-cap15,metrics-cap15,effects}.npz/csv/json/parquet`.

## Peak memory and wall time

Wall ≈ 4 s (known) + ≈ 3 min (order, dominated by 100-seed random walks) +
≈ 7 s (profile) + ≈ 20 s (bootstrap/anatomy); total < 5 min. Peak RSS <
1 GiB (largest object 85 × 10,779 float64 ≈ 7 MB; long frame ≈ 0.9 M rows
polars, transient). `memory_pressure` free was 75% at start; no guard needed.

# 342 — Hamilton integer-total rounding, remeasured on the current family

Registration, September 9, 2026. Initial implementation by Claude (no
execution); reviewed and corrected by Codex before freeze. Codex runs and verifies.
Freeze before any
generation. No model fitting, new source, reserve read, network, production
default change, package or submission. This is a bounded remeasurement of an
already-tested, still-unsubmitted rider, not a new method.

## Question

The production emitter delivers per-gene panel totals `floor(n * L * profile_p)`
(`src/vcc2026/donor_transfer.py`, `emit_counts`, `mode="quasi"`). Every gene is
short by its own fractional part, a fixed count that reads as a size-correlated
tilt in ratio metrics. Does largest-remainder (Hamilton) rounding of those totals
improve the CURRENT family (323's cap15 composition, rc_callset emission,
H1 union OFF) on actual emitted counts and all six members, beyond a control
that restores the same integer mass at random coordinates?

## Prior evidence this does not repeat

- 196 (`ham196.py`): on a15/rc_hybrid, Hamilton removed the h1-test size tilt
  (log2 −0.01352 → +0.00019) and bought NMAE and reach.
- 208 (`emit208.py`, RESULTS): on rc_callset, raw h1-test NMAE 0.96670 → 0.96040,
  reach +0.00006, fidelity −0.00007, Jaccard +0.00005, PLOLO PDS +0.0032 raw,
  raw MSE −0.00001; totals differ from the incumbent on 900/900 challenge blocks,
  every move +1 and upward, ~6,955 of 18,533 genes per block, panel mass +0.088%.
  Verdict then: a rider, vetoed only by rule 13's toll.
- 230: the toll was mis-keyed; Hamilton's delivered profile is 1.91× closer to the
  intended one (relL1 8.78e-4 → 4.60e-4). Left live, uncleared, never submitted.

Those measurements used older compositions and the emitter's own DE arithmetic.
None used 323's exact evaluator, the cap15 effects, the H1-excluding operator,
the current CD4/PM composition, or a random-allocation control.

## Fixed construction (identical to 323 cap15)

Effects: `experiments/323-current-family-validation/work/effects.npz` entry
`cap15`, mapped by symbol to the H1 reference axis exactly as `evaluate323`
does; hash checked against 323's `compose.json`. The SAME effects feed all
three arms. Alpha .982, CD4 f .34 / k 5, PM production default, whitening
gamma .65 / lambda 0 on the H1-excluding operator, H1 union OFF, mass fix unset.

Emission: production `rc_hybrid_block` with prediction/library-pick seed 0,
emission seed offset 1,000,003 (one stream across the 85 targets, as in
production), jitter seed 1000 + i, 400 cells, min CPM 10, tail spread 1,
keep 8000 / spread .25 block map, overdispersion 4, dispersion .1,
call set head 400 / CPM 2200 / z cap 1.9, `call_set_first=True`.
Controls, libraries, CPM gates and rank-control inputs come from 323's frozen
`reference.npz`; truth never enters construction.

## Three arms, one variable

Let `x = n * L * profile_p` at the production rounding site, `base = floor(x)`,
`k = round(sum(x)) - sum(base)` (Python round, half to even).

| arm | rule |
|---|---|
| `floor` | `np.floor(x).astype(np.int64)`, production's expression verbatim |
| `hamilton` | `k` increments to the largest fractional parts; ties to the lower gene index (`argsort(-frac, kind="stable")`) |
| `random_up` | `k` increments uniformly without replacement among coordinates with fractional part > 0, drawn from `np.random.default_rng(3420909 + target_index)`; target_index is the position in 323's sorted reference target order |

Every coordinate of every arm differs from `base` by 0 or 1. The rounding choice
never consumes the emission generator: the seamed statement hands the rule ONE
argument, the intended vector, and the control's generator is constructed
inside the per-block adapter. `random_up` restores the same integer mass as
Hamilton at random coordinates, so the pair separates common-mode mass
restoration from remainder-targeted allocation.

## The seam (rounding342.py, seam342.py)

`rc_hybrid_block` calls `emit_counts(mode="quasi")` twice per block: pass 1 on
the full profile (its totals drive the |z| multiset, the per-cell weld, the
designed histograms and the delivered-total assertion) and pass 2 on the
profile with scorer-gated genes zeroed (the sub-gate tail block). The adapter
`BlockRounder` decides ONCE on the first call and serves later calls only if
they are coordinate-wise restrictions of the first, returning the restricted
decision; anything else raises. Floor is restriction-consistent by itself and
computes production's expression on every call.

`seam342.Seam` parses production `donor_transfer.py`, replaces exactly one
statement inside the `mode == 'quasi'` branch
(`totals = np.floor(n * L * profile_p).astype(np.int64)` →
`totals = __vcc342_round_totals__(n * L * profile_p)`), compiles the result
into a fresh namespace that borrows production's `_quasi_block`,
`_repair_totals` and `_stochastic_round`, and binds one rounder per block into a
fresh function object. No global NumPy or module monkeypatch; the production
module object is untouched; the other three `totals = np.floor(...)` sites stay.
The attestation (production source hash, the single line difference, 4 → 3
untouched floor sites, seamed-function hash) is hash-bound at freeze and
re-derived at verify. Counts are never added after emission: the new totals go
through the real pass-1 / rank / placement / weld machinery, which asserts
delivered totals equal the rule's totals per block.

## Integrity requirements before any changed arm is read

1. Pure tests (`test342.py`): known fractional totals, zero/integer entries,
   stable ties, mass conservation, 0/1 differences, random-control mass and
   dedicated generator, unchanged floor path through `emit_counts` and through
   `rc_hybrid_block` (bit-identical blocks, identical emission-stream state),
   Hamilton totals flowing through the machinery, injection-site specificity
   (other modes bit-identical with a trapping hook, unbound seam refuses).
   `verify` runs them.
2. Re-emit all 85 cap15 blocks under `floor` through the seam and require, against
   323's recorded cap15: the whole emitted-count digest, the panel-total digest,
   per-gene totals, every saved statistic (`sums`, `jk`, `p`, `q`, `lfc`), every
   per-target row, both MSE legs and all six raw metrics EXACTLY (zero
   difference). 323's recorded values: PDS .8896358543417369, MSE
   .9749462418923823, NMAE .9508364630325662, fidelity .4975545182292535, reach
   .25434726355678233, Jaccard .07475983298131939. Any inequality stops the run
   with the failure preserved in `work/known-answer.json`; the requirement is
   not weakened to a tolerance.
3. Per block and arm: rounder called exactly twice, delivered and scored totals
   equal the rule's totals, mass matches the rule, 0/1 differences, increments only at
   fractional coordinates, moved count equals `k` (0 for floor).
   Floor deliberately retains `sum(floor(x))`, NOT `round(sum(x))`. Only the
   two changed arms restore rounded mass. Codex caught and corrected the
   draft's false integrity failure for the intended floor deficit before
   freeze or response measurement, with a new pure regression test.
4. Bootstrap identities per arm: full-panel MSE and official `_across_pert_budget`
   reproduced within 1e-10, per-target means within 1e-12, no non-finite
   resample, positive denominators.

## Measurements

Scoring is `exact_de323` plus cell-eval2's own PDS, jackknife MSE and DE
metrics, exactly `evaluate323.block_stats` / `evaluate323.score`; never the
emitter's older DE arithmetic. For each arm save `stats-{arm}.npz`
(sufficient statistics), `per-target-{arm}.csv`, `metrics-{arm}.json`
(six raw metrics, digests, MSE legs) and `rounding-{arm}.json` (per-target `k`,
moved genes inside/outside the scorer gate, panel total, call-set delivered
size, and for `random_up` the Jaccard of its chosen set with Hamilton's).

Report all six raw metrics per arm and the H1-local scaled overall (323's
`h1-anchors.json`, MSE clamped to [0,1], NMAE floored at −6). Raw MSE and its
jackknife credit are recomputed for every arm and in every resample, never
assumed pinned. Paired 10,000-target resamples, seed 20260909, identical across
arms and members, whole-panel MSE ratio and claim cap recomputed per resample.
Contrasts: hamilton − floor, hamilton − random_up, random_up − floor, each with
overall and per-member points and 95% percentile intervals. No board conversion
factor and no A/B/C forecast.

Rounding changes the pseudobulk totals of every block. No competition member,
including PDS and MSE, is assumed invariant; the report records how many blocks
keep floor's panel-total digest (expected 0 of 85 for both changed arms).

## Prospective gates (LOCAL_GAIN_REQUIRES_SEPARATE_ABC_REVIEW iff all pass)

1. All integrity requirements above pass.
2. Hamilton local H1 overall point ≥ +0.001 versus floor.
3. Paired overall 95% lower bound > 0 versus floor AND versus random_up.
4. No raw MSE or NMAE contrast (hamilton versus floor or random_up) whose entire
   95% interval is harmful (lower bound > 0 for these lower-is-better members).
5. No raw fidelity, reach or Jaccard contrast (same two comparators) whose entire
   95% interval is harmful (upper bound < 0).

A pass means only that a local H1 gain exists and requires a separate A/B/C
review; it is not a cleared package, slot, or production change. Any failure is
recorded as NO_LOCAL_ADVANCEMENT and preserved; nothing is retuned, no arm is
added, and no gate is relaxed. PDS contrasts are reported with the same harm
flag but are not gated by this registration.

Expected directions from 208, stated so they can be falsified: raw NMAE falls by
roughly 0.006 under Hamilton, raw PDS rises by a few thousandths, fidelity,
reach and Jaccard move within noise. If `random_up` matches Hamilton, the gain is
common-mode mass restoration; if Hamilton beats it with a positive lower bound,
the remainder targeting itself carries the gain.

## Resources and stopping

One process, under `scripts/run_with_memory_guard.py --max-rss-gib 16` and the
repository launch gate; 40-minute in-process budget (2,400 s, checked before each
arm, target block and report); the run must not start after 19:35 UTC on
September 9 without a new coordinator decision passed as `--decision`; handoff
20:58 UTC. Arms run in the order floor, hamilton, random_up, and the floor
known answer must pass before the second arm starts. Every output stays in
`experiments/342-hamilton-current-family/work/`; binary statistics are
git-ignored, compact receipts are retained.

## Not selected, not implemented

The empirical-null call-set proposal and block-sign repair were not selected and
are absent from this experiment.

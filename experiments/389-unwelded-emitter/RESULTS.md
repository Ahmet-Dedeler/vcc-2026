# 389-unwelded-emitter — unwelding the DE channel from the bulk channel: mechanism works, no significant gain at feasible doses

- Status: COMPLETE (S0, S1, S2 + controls done; S3 skipped — 388 fired its kill criterion, no qualifying sign estimator).
- One-line verdict: **NEGATIVE with a working mechanism.** K-heterogeneous libraries decouple the arithmetic-mean (DE) channel from per-gene totals exactly as designed — pds/mse bit-identical in all six K-class arms, steering exact to 1e-15 — but steering the call set at the largest class-feasible dose (≈1.2×, not the nominal 1.48×) moves H1-local overall by +0.00004 with CI [−0.00106, +0.00089]: null. Nominal 1.48×/2.0× doses are jointly class-infeasible on 81/85 targets (per-gene targets fine; the class-level pigeonhole binds). The specified shuffled control collapses to the unsteered state (deranged joint demand infeasible 85/85) and is reported as such.

## Known-answer reproduction

- S0 PASS (`work/known-answer-cap15.json`): `run389.py cap15` reproduces
  `experiments/323-current-family-validation/work/metrics-cap15.json` — emitted-count
  digest equal, panel-total digest equal, mse legs exact; absolute diffs
  fid/jac/mse/pds/reach 0.0, nmae 1.1e-16 (threaded-scoring summation order, well
  within 1e-6 tolerance).
- K=1 PASS (`work/known-answer-k1.json`): `rc_hybrid_block_k` with K=1 reproduces
  cap15 with all six absolute diffs exactly 0.0 and equal emitted-count digest —
  the generalised function is bit-identical to production through its K=1 branch
  (verbatim prefix + literal `EME.design`/`EME.place`, identical RNG draw order).
- Kernel proof (`work/kernel389.json`): `hist_z_k == wilcoxon_z` to 8.9e-16 on real
  control rows (3 synthetic multi-class panels × 24 genes). The per-class
  generalisation (own Phi table per class at v=log1p(1e6·c/L_k); pairs grouped by
  exact CPM value) is scorer-exact.

## What was measured (inputs, estimator, controls)

- H1 development panel (85 targets × 400 cells), 385's exact cap15 chain; only the
  block function varies (`run389.py`). Effects/reference/builder read-only from 323;
  `CallSet(400, 2200, 1.9)`; picks stream identical across arms.
- Classes: 3 equal-fraction classes at h1-test control library q25/50/75, measured
  from `data/harness/h1-test.h5ad` non-targeting cells (n=18,400, mean 21,253,
  CV 0.351; q = 16079/20191/25331), scaled per block to sum S (e.g. target 0:
  S=8,245,700 → levels 16151/20281/25444, counts 134/133/133).
- Arms: `cap15` (production), `k1` (K=1, must match), `k3plain` (S1: K=3, no
  steering), `k3s125/k3s148/k3s200` (S2: call-set steered to s× design LFC in log2
  space, delivered at per-target feasible dose s*≤s), `k3shuf` (specified shuffled
  control: resolved targets deranged, seed 20260913).
- Statistics: paired 10,000-target bootstrap, seed 20260913, H1-local scaling with
  323 anchors, MSE ratio + claim cap recomputed per resample (`report389.py` →
  `work/results389.json`; staged `results389-s1.json`, `results389-s2.json`).
- Invariants asserted per block: (1) per-gene totals bit-identical to production;
  (2) per-cell totals equal class levels; (3) pooled |z| re-measured exactly;
  (4) steered genes' delivered arithmetic mean within half a one-count move
  (clips flagged, exempted). All held on all 85×6 K-class blocks (0 unresolved,
  max steering err 9.6e-16).

## Results

### Six raw members per arm (H1 panel)

| arm | pds | mse | nmae | fid | reach | jac | H1-local overall |
|---|---:|---:|---:|---:|---:|---:|---:|
| cap15 | 0.889636 | 0.974946 | 0.950836 | 0.497555 | 0.254347 | 0.074760 | 0.280953 |
| k1 | same | same | same | same | same | same | 0.280953 |
| k3plain (S1) | same | same | 0.950814 | 0.497611 | 0.256617 | 0.074765 | 0.281415 |
| k3s125 (deliv ~1.20) | same | same | 0.949111 | 0.497562 | 0.253062 | 0.074765 | 0.281202 |
| k3s148 (deliv ~1.20) | same | same | 0.949286 | 0.497557 | 0.252240 | 0.074765 | 0.280990 |
| k3s200 (deliv ~1.20) | same | same | 0.949486 | 0.497562 | 0.252265 | 0.074765 | 0.280939 |
| k3shuf (collapsed, see Controls) | same | same | 0.950814 | 0.497611 | 0.256617 | 0.074765 | 0.281415 |

pds/mse bit-identical in every K-class arm (per-gene sums digest equal cap15's;
mse legs exact — the jackknife is placement-blind at 16 digits). "same" above is
exact, not rounded.

### Paired contrasts vs cap15 (point [95% CI], H1-local overall + raw members)

| contrast | overall | nmae | fid | reach | jac | mse/pds |
|---|---:|---:|---:|---:|---:|---|
| k1 − cap15 | 0 [0, 0] | 0 | 0 | 0 | 0 | 0 |
| k3plain − cap15 (S1) | +0.00046 [−0.00025, +0.00144] | −0.00002 [−0.00009, +0.00005] | +0.00006 [−0.00005, +0.00016] | +0.00227 [−0.00135, +0.00731] | +0.00001 [−0.00000, +0.00001] | 0 |
| k3s125 − cap15 | +0.00025 [−0.00086, +0.00114] | −0.00173 [−0.00370, +0.00024] | +0.00001 [−0.00007, +0.00008] | −0.00129 [−0.00618, +0.00212] | +0.00001 | 0 |
| k3s148 − cap15 | +0.00004 [−0.00106, +0.00089] | −0.00155 [−0.00368, +0.00054] | +0.00000 [−0.00008, +0.00008] | −0.00211 [−0.00660, +0.00042] | +0.00001 | 0 |
| k3s200 − cap15 | −0.00001 [−0.00114, +0.00087] | −0.00135 [−0.00358, +0.00085] | +0.00001 [−0.00007, +0.00008] | −0.00208 [−0.00657, +0.00045] | +0.00001 | 0 |
| k3shuf − cap15 | identical to k3plain (bit-identical statistics; nmae 1 ulp scoring noise) | | | | | |

### S1 lattice cost (kill rule: loss ≥ 0.005 overall, CI excl. zero)

Not triggered: S1 overall is +0.00046 [−0.00025, +0.00144] — a full order of
magnitude above the kill floor, CI including zero. Heterogeneity alone costs
nothing resolved: |z| drift median 0.018, p90 0.13, max 1.37 (916k genes);
delivered call set Jaccard vs incumbent median 1.000, min 0.998 (79/85 targets
identical); delivered size median 2600 both arms. Placement works harder
(forced ~3.6k vs 0, repair moves ~3.6k vs 0) but resolves everything
(unresolved 0/85; post-repair cell totals exact by assert).

### Reachable-LFC ceiling (kill rule: |ΔLFC| < 0.1 log2 for typical call-set genes)

Not triggered: log2(Amax/Amin) = 0.656 for EVERY gated gene with M>0 (n=916,214;
call-set subset n=220,999 identical) — every gene fits in any single class
(n_k·L_k ≈ 2.2M ≫ max M), so the range is exactly log2(Lmax/Lmin). Predicted clips
at optimistic caps (221k call-set gene-targets): 1.25× → 51 (0.02%); 1.48× → 297
(0.13%); 2.0× → 4,262 (1.9%).

### S2 dose reality: the class-level pigeonhole binds (the experiment's main mechanism finding)

Per-gene 1.48× targets are individually feasible (ceiling 0.656 log2), but JOINTLY
the steered demand exceeds class mass: on target ACTN1 the 1.48× split demands
+16k net counts into class0 beyond its true capacity (sub-gate mass occupies the
rest; even fully drained, gated_0 ≤ mass_0 is hard — measured gated_0 − mass_0 =
+1246 with only 1292 raise vs 1308 lower genes: up-gene mass dominates). Global
repair then destroys per-class balance chasing global extremes and stalls on
all-gated over-full cells (diagnosed from residuals: class0 uniformly −35, no
sub-gate counts to swap). Fix: per-target feasible-dose bisection
(`feasible_dose`, margin 6000 counts/class for fragmentation slack) delivering
s* ≤ nominal, exact at s*. Delivered dose: k3s125 med 1.197 (26/85 at nominal);
k3s148 med 1.197, range 1.032–1.480 (4/85 at nominal); k3s200 med 1.195, max 1.707
(0/85 at nominal). **All three nominal doses collapse to the same ≈1.2× delivered
dose — there is no dose curve to read above ~1.25× at these classes.** Backstop
`rebalance_within_class` (gated within-class moves, M_k-preserving) never fired
(0 moves all panels; margin sufficed).

### S2 verdict

At delivered ≈1.2× on the 2,600-gene call set: nmae −0.0014…−0.0017 (all CIs
include zero; best k3s125 [−0.00370, +0.00024]), reach −0.0013…−0.0021 (all CIs
include zero), fid flat to ±0.00008 (steering scales magnitudes; signs come from
the design, unchanged), jac +0.000005 (call set pinned by the |z|-multiset
re-key), overall +0.00025…−0.00001 (all CIs include zero). Tension (inference, not
a claim): amplifying designed magnitudes helps magnitude-weighted nmae
directionally while the changed realised-|z| order degrades reach directionally;
nothing resolved. Consistent with 343 (signs are the gap, amplitude is not) and
with 388 (below).

### Controls

- Shuffled control v1 (`k3shuf`, as specified: resolved targets deranged across
  genes, seed 20260913): the deranged joint demand is class-infeasible at margin
  on 85/85 targets → every block fell back to the mass-proportional split
  (dose 1.0, no steering asserts) → statistics bit-identical to `k3plain`
  (sums/jk/p/q/lfc arrays exact; both digests equal; five members JSON-equal,
  nmae 1 ulp threaded-scoring noise). It "does not gain", trivially — but it
  tests the S1 lattice, not deranged steering. Honestly reported as a collapsed
  control, not a passed one.
- Shuffled control v2 (deranged design-LFC multipliers — same marginal dose, no
  gene-correct information): sampled 3 targets, all infeasible (dose 1.0
  fallback). Spot measurement of deranged s=1 demand vs hard class mass:
  +6k…+40k over on some class per target (targets 0/1/5/40; unmasked demand).
  Deranged steering is structurally unplaceable with exact cell totals at these
  classes — panel not run. This bounds the control the assignment asked for:
  only near-true-direction joint demand fits the pigeonhole.
- S3 skipped: 388 (`experiments/388-de-channel-prize/RESULTS.md`, COMPLETE)
  fired its kill criterion on both clauses — no estimator beats E0's sign
  accuracy by ≥0.03 (best E1 +0.012); the decoupled channel's priced prize is
  amplitude-only. No qualifying sign estimator exists to steer by.

## Verdict against the PLAN.md kill criteria

- S1 loss rule (≥0.005 H1-local overall loss, CI excl. zero): NOT triggered
  (+0.00046 [−0.00025, +0.00144]). Correctly proceeded.
- Ceiling rule (reachable |ΔLFC| < 0.1): NOT triggered (0.656 everywhere).
- Implicit S2 bar (resolved overall gain at feasible dose): NOT MET — overall
  null at all three nominal doses (delivered ≈1.2×). The channel is mechanically
  decoupled and priced; it buys nothing significant without better signs.

## Limits

- H1-panel numbers are development evidence, not A/B/C forecasts; no board
  conversion applied.
- The feasible-dose boundary (≈1.2×) is specific to K=3 equal-count classes at
  h1-test q25/50/75 levels; wider classes (e.g. 5/95% quantiles) move more
  arithmetic mean per count and would bind later — not tested (would break
  S1/S2 lattice comparability; flagged for the PI, not claimed).
- Range/clip precomputation (`range389.py`) used optimistic caps and 323's
  recorded totals; delivered clips/doses are per-block asserted instead.
- `place_k` gained the gated-rebalance backstop + sub-first ordering AFTER S1;
  both are gated no-ops on the S1 path (backstop fires only when unresolved>0,
  which S1 never had; unsteered split/repair code byte-untouched) — S1 numbers
  stand without re-run. K=1 bit-identity re-proven after the refactor (target-0
  twin + full-panel digest equality predates the refactor; the refactor touched
  only K>1 branches and factored shared code — verified by re-running the
  target-0 twin post-refactor: identical).
- H1-local overall uses 323's anchors (same as 342/385).

## Exact commands to reproduce

- `uv run python experiments/389-unwelded-emitter/run389.py <cap15|k1|k3plain|k3s125|k3s148|k3s200|k3shuf>` — one arm at a time (~140–185 s, ~5 GiB; memory-gate `memory_pressure | tail -1` free ≥ 40% first)
- `uv run python experiments/389-unwelded-emitter/report389.py k1 k3plain k3s125 k3s148 k3s200 k3shuf` → `work/results389.json`
- `uv run python experiments/389-unwelded-emitter/range389.py` → `work/range389.json`
- `uv run python experiments/389-unwelded-emitter/test389_smoke.py` (kernel proof + target-0 twin + K=3 drift; ~50 s)
- All with `NUMBA_CACHE_DIR=experiments/389-unwelded-emitter/work/numba-cache`, 4 threads; env set inside the scripts.

## Peak memory and wall time

- Wall per arm (from receipts): cap15 494 s, k1 302 s (both under CPU contention
  with the 387 production build + first numba compile), k3plain 141 s, k3s125
  168 s, k3s148 183 s, k3s200 184 s, k3shuf 159 s. Bisection adds ~1.5 s/target
  on steered arms. Bootstrap report ~30 s. Smoke 48 s.
- Memory: 12 GiB RSS cap respected throughout; system-wide free never below 61%
  at any gate check (production build ~22 GiB co-resident). No per-process RSS
  instrumentation was added; no OOM, no spill, no scratch left behind
  (numba-cache 2.4 MiB retained for rebuilds).

# 362 — PLAN: six-member price of a same-line common vector (oracle ceiling + transferable comparison)

## Hypothesis
345 measured at the DESIGN stage that adding the same-line LOO common vector
`Cmean` (leave-one-target-out mean of the other 84 H1 targets' reference lfc) to the
cap15 design raises |lfc|-weighted sign accuracy on the reference gate 0.600 → 0.656
at α=1 (scale s = pooled-median |D| / pooled-median |Cmean| ≈ 2.7549), exact NMAE
0.9485 → 0.9338. No donor-line common vector carries any of it (best transferable
+0.0026, grid-max, vs +0.02 gate; PI re-derivation: K562/Jurkat/HepG2 commons score
0.527–0.534 weighted vs design's 0.600).

For context A (= Jurkat) a same-line common vector IS buildable from the Nadig Jurkat
bank (`experiments/350-same-line-priors-data/work/jurkat_common.npz`), but 350 killed
it as a direction (generic CRISPRi axis, cosine 0.331 with HepG2's, −0.011 with our H1
design), and 155 showed a uniform additive down-shift bought reach at 4× its pds cost.
Before anyone builds for A, this experiment puts the H1-oracle construction through
the FULL six-member emission scorer: does the +0.056 sign gain survive emission, and
what does pds pay? This is an ORACLE (uses other targets' H1 truth) — a development
ceiling, not a candidate. A clean negative closes the route for every line.

## Exact inputs (all read-only; nothing outside this directory is written)
- `experiments/345-common-program-structure/work/stage345.npz`: `Cmean` (85×10,779,
  per-target LOO mean of other targets' gate-masked reference lfc on the universe),
  `Dfull` (cap15 design on the universe), `genes` (universe, 10,779; first genes
  NOC2L/ISG15/AGRN…), `targets` (85, same order as effects), `elig` (58 eligible target
  indices, int array), `gate` (85×10,779 bool), `R` (reference lfc on gate else 0).
  Verified: `elig` dtype int64, shape (58,); targets identical to effects targets.
- `experiments/323-current-family-validation/work/effects.npz`: `cap15` (85×18,533 on
  `genes` TSPAN6/TNMD/DPM1… = 319 composition) — the design the harness emits. Universe
  genes mapped to this axis BY NAME; the 2 universe genes absent from the effects axis
  (`TBCE-1`, `HSPA14-1`, per 345) have no slot and are dropped (their Dfull entries are
  0; sign 0 = miss, exactly the 343/345 convention).
- `experiments/350-same-line-priors-data/work/jurkat_common.npz`: `genes` (18,533
  challenge axis), `vector` (mean log2-CPM-ratio over 300 seeded-random Jurkat targets),
  `mapped` mask (8,284 mapped). Mapped to the universe by gene name; 0 where unmapped.
- 323 harness (imported read-only, never written): `run323.py` (`builder`, `map_axis`),
  `evaluate323.py` (`load_reference`, `block_stats`, `score`, emission constants),
  `report323.py` (`mse_bootstrap`, `scaled`, `ci`), `work/reference.npz`,
  `work/reference_de.parquet`, `work/metrics-cap15.json` (known answer),
  `work/per-target-cap15.csv` (known answer), `work/stats-cap15.npz` (known answer),
  `work/h1-anchors.json` (H1-local scaling anchors).
- Template for own-directory reuse: `experiments/342-hamilton-current-family/run342.py`
  (`emit_arm`: replicates `evaluate323.emit_and_score` bit-exactly but writes into the
  experiment's own `work/`) + `report342.py` (paired 10,000-target bootstrap, seed
  20260909, MSE ratio recomputed per resample).

## Estimator
1. **Scales (recomputed, stated):** s_oracle = pooled-median |Dfull| over the eligible
   gate / pooled-median |Cmean| over the eligible gate (345's construction;
   expect ≈2.7549). s_generic = same rule with the Jurkat universe-mapped vector in
   place of Cmean (same J for all targets; pooled over eligible gate). The shuffled arm
   reuses s_oracle (a gene permutation preserves the pooled magnitude multiset exactly).
2. **Design-stage weighted sign accuracy (343's functional)** for every arm on the
   modified universe rows: per eligible target, gate from `reference_de.parquet`
   (q<.05 & finite & nonzero lfc & gene≠target), wacc = Σ|l|·[sign(pred)==sign(l)]/Σ|l|,
   mean over the 58. Predicts the emission read; checks the +0.056 oracle gain is
   actually installed before any scorer time is spent.
3. **Modified full-axis designs** (universe genes only; all other effects-axis genes
   bit-identical to cap15):
   - `cap15` — unmodified (known answer).
   - `oracle_a05` — D[t,g] + 0.5·s_oracle·Cmean[t,g].
   - `oracle_a1` — D[t,g] + 1.0·s_oracle·Cmean[t,g].
   - `oracle_a1_normpreserve` — oracle_a1 row rescaled on the universe genes to the
     original row's norm there (n0=‖Dfull[t]‖, n1=‖modified[t]‖, ×n0/n1; separates sign
     information from added energy).
   - `generic_a1` — D[t,g] + 1.0·s_generic·J[g] with J = Jurkat common vector on the
     universe (0 where unmapped).
   - `shuffled_a1` — oracle_a1 with one GLOBAL gene permutation π of the 10,779 universe
     genes applied identically to every Cmean row (seed `default_rng(362)`), same
     s_oracle: added energy without information.
4. **Emission + scoring:** own-directory replica of `evaluate323.emit_and_score`
   (production builder `R.builder()`, H1 union OFF; cells 400, min_cpm 10, tail_spread
   1.0, keep_top 8000, keep_spread 0.25, block map [0, 0.60], dispersion 0.1,
   overdispersion 4.0, CallSet(400,2200,1.9) call_set_first, jitter 1000+i, seed 0 —
   the exact 342 constants), effects mapped to the reference axis by `R.map_axis`,
   `E.block_stats` + `E.score` per arm. ~142 s / ~5 GiB per arm, sequential, one process.
5. **Reporting:** all six raw members per arm, H1-local scaled overall with 323's anchors
   (`report323.scaled`: mse clipped [0,1], nmae floored −6; overall = mean of six
   scaled), paired 10,000-target bootstrap (seed 20260909) Δ vs cap15 for each member
   and the overall (report342 pattern, MSE panel-ratio recomputed per resample).

## Controls / falsifiers
- **Known answer (gate, must pass before any changed arm is interpreted):** unmodified
  cap15 must reproduce `metrics-cap15.json` EXACTLY — emitted-count digest
  `547cf0ab…`, panel digest `f699b3d5…`, all six raw members
  (pds 0.8896358543417369, mse 0.9749462418923823, nmae 0.9508364630325662,
  fid 0.4975545182292535, reach 0.25434726355678233, jac 0.07475983298131939),
  mse_legs, every stats array, every per-target row. If it does not match, stop and
  report that instead of any arm result.
- **shuffled_a1** (energy without information): any member movement shared by
  oracle_a1 and shuffled_a1 is an amplitude/energy effect, not sign information.
- **oracle_a1_normpreserve** (information without energy): separates the two.
- **generic_a1** (transferable axis): PI's re-derivation predicts ≈ no sign gain
  (0.534); prices what a cross-line generic axis does to the six members.
- Design-stage wacc must reproduce 345's grid (α=0: 0.6002, α=0.5: 0.6527,
  α=1: 0.6564) before emission; otherwise the construction is wrong, not the idea.

## Kill criterion (from the mandate; the verdict section must answer it verbatim)
The route is worth building for A only if `oracle_a1` (or `_normpreserve`) improves
the H1-local overall with a paired interval above zero AND pds's raw loss is smaller
than the summed gains on nmae+fid+reach in scaled units. If pds pays more than the DE
members gain even for the oracle, the route is dead for every line (record it so).
Also state whether `generic_a1` moves anything — if it improves mse or pds without
sign gain, say so explicitly.

## Resources / limits
Memory 71% free at plan time (`memory_pressure`); cap 12 GiB RSS, threads
`NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4`, sequential arms, no full-matrix densification
(emission blocks only), scratch deleted afterwards. No reserve paths, no paid APIs, no
cloud, no downloads. H1-panel numbers are development evidence, not A/B/C forecasts;
no retrieval/cosine/proxy gain will be called a score gain.

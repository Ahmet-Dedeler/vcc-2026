# 194 — Re-estimate the donor effect tables. They were built on 2026-08-22 and never revisited, and everything sits on them.

You are a research subagent. Work ONLY in `experiments/194-donor-reestimation/` (new donor tables
go in `data/donors/194-*/` — that is the ONE exception).
No `vcc`, no `cell-eval2`, no submission build, no changes to any existing file in `data/donors/`.
Peak RSS ≤ 14 GB — **stream, never materialise a cell-level matrix.** `uv run python`.
This is a **long compute job** and that is fine; use the wall clock.

## READ FIRST
`AGENTS.md` §5 (data), §8.1, §9 (64 GB is a real constraint), §11.1 rules 1–13. Then:
- `scripts/build_donor_table.py` and `scripts/build_orion_donor.py` — **the estimators under test.**
- `experiments/166-donor-acquisition/RESULTS.md` — no further panel-covering line exists.
- `experiments/159-transfer-ceiling/RESULTS.md` — κ² ≈ 0.34, we sit at ~30% of the infinite-donor
  context-independent ceiling, and the `1/k` law.
- `experiments/184-tahoe-transfer/RESULTS.md` — the line-specific residual is **not predictable**,
  closed in the strong form on 50 cell lines. So the shared component is all there is, and
  estimating it better is the only route left on the direction channel.
- `experiments/010-lodo`, `016-modz`, `039-coverage-weight`, `015-pooling`, `022-graveyard`.

## The gap

Every number in this project rests on three files built on 2026-08-22 by a **plain pseudobulk
difference**: sum the counts of a target's cells, sum the control's, normalise each to TS, take
`log1p` of both, subtract. `min_cells = 20`. For K562-gwps there is **no guide-level QC and no
knockdown weighting at all**; Orion filters on the dataset's own `pass_guide_filter`. Median cells
per target: K562 **179**, HCT116 **154**, HEK293T **203**; minimum **20** in all three.

That estimator has never been compared against an alternative, and it is upstream of `pds`, `mse`,
`nmae` and every sign. 159 says we are at ~30% of the context-independent ceiling and 184 says the
line-specific 70% is unreachable — **so the only remaining direction lever is estimating the shared
component better from the cells we already have.**

Candidates, in rough order of expected value. Each is a standard move in this literature and none
has been tried here:
1. **Weight or filter cells by observed on-target knockdown.** A guide that did not repress its
   target contributes a control cell to the perturbed pool and dilutes the effect toward zero. The
   knockdown is directly observable as the target gene's own expression drop. **Note the boundary:
   the target's own row is excluded from every scored metric, so using it to weight cells costs
   nothing and throws away no scored signal.** The 2026 challenge itself used a >80% knockdown
   criterion (§1).
2. **Guide-level aggregation before target-level.** Average per guide, then across guides, rather
   than pooling all cells — this down-weights guides with many cells and bad knockdown, and is what
   Replogle's own pipeline does.
3. **Robust / shrunken location** instead of the pooled sum: trimmed mean, median-of-guides, or
   MODZ over guides (016 tested MODZ at the *consensus* level, not within a donor).
4. **Batch / plate structure.** The Orion long-gene technical axis (`AGENTS.md` harness warning) is
   a library artefact that loads on the donor tables at 9.8% of the un-whitened consensus energy.
   Removing plate effects *within* a donor before differencing is a different operation from
   deflating the axis afterwards.
5. **Cell-quality filtering** (UMI depth, mitochondrial fraction, doublets) and `min_cells`.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
Rebuild **one** existing donor table with the current estimator and require it to reproduce the
shipped `delta.npy` **bit-identically** (or, if the build is not deterministic, to ≤1e-6 on every
recorded summary in its `meta.json`). If you cannot reproduce the incumbent you cannot measure an
improvement against it. Record every source file's sha256.

### Phase 1 — the estimators
Build each candidate as a variant table. Keep the gene axis, the target list and the output format
byte-compatible with the incumbent so downstream code is unchanged. Report cells retained, cells
reweighted, and the change in each table's own effect-norm distribution.

### Phase 2 — score them where it counts
The endpoint is **cross-line transfer**, not within-line fit — a better within-line estimate that
does not transfer is worthless here. Score on:
- **PLOLO** `pds` (162's instrument, the only one validated against a board `pds` point, four donor
  folds, challenge frame). **`debench.Surface.pds()` must never size a `pds` arm** (188 Phase 3).
- **Per-gene sign accuracy `q`** on held lines, and the **|t|-weighted `q_w`** (187 §5: `nmae`
  moves only through `q_w`).
- 159's **κ²** and the effective-donor count, recomputed on the new tables — if a better estimator
  raises κ², that is the cleanest possible statement of the gain.

Hold out the scored line from everything, including the whitening basis. `data/donors/h1_2025`
covers 85/85 of h1-test's targets and must never enter an h1-test screen.

### Phase 3 — falsifiers
1. A **permuted-weight control** for the knockdown weighting: shuffle the weights across cells
   within a target. The real weighting must beat it.
2. **Does it survive whitening?** 022 found that under whitening, EB shrinkage, inverse-variance
   pooling, low-rank truncation and sign-only all kept their exact deficits, and only donor
   weighting changed rank. Score every arm **through the production operator**, not before it.
3. **Rule 12**: an estimator change moves the effect's shape, so price `mse` (096's `claim` route)
   and `nmae` (h1-test, **significance-selected pool** per 188) as well as `pds`.

## Deliverable
`RESULTS.md`: the reproduction row, each estimator's transfer scores against the incumbent, κ²
before and after, the falsifiers, and a verdict. If the plain pseudobulk difference is already the
best of these, say so — that is a genuine and useful closure of the last direction lever.

## Territory
192 owns the call set. 193 owns the whitening operator. 195 audits. You own the **donor tables**.

# 380 — PLAN

## Record check (done first, per the assignment)

Read `experiments/156-*/RESULTS.md`, `175-*/RESULTS.md`, `204-*/RESULTS.md`, `209-*/RESULTS.md`,
`246-*/RESULTS.md` §0/§1, `249-*/RESULTS.md`, and `345-*/RESULTS.md`. Also inspected
`experiments/323-current-family-validation/work/compose.json` directly to see what `cap15`
actually contains today.

**What each experiment deflated, on which object/surface, and its verdict:**

| exp | object deflated | surface | verdict |
|---|---|---|---|
| 156 | the **render's** mass-fix correction geometry / a pre-render "centroid" (flat cross-target mean) on two held Orion donor panels | exact-space pseudobulk truth (HCT116, HEK293T), deterministic, no emitter | centroid removal alone is a **no-op** (+0.00000 / +0.00098 pds); the real fix (`massfix-minnorm`) is a different mechanism (mass-constraint direction), not a shared-subspace deflation |
| 175 | per-**donor** rank-1 PCA (each of K562/HCT116/HEK293T's own top eigenvector), pre-consensus | H1-test sign-accuracy pool (donor tables only, no render) | local sign-accuracy **gain** (+0.0073 consensus) — but this is a proxy (sign accuracy on a pool, not the 6-member harness) |
| 204 | the **donor table's** target-independent mean ("blob"), checked as a nuisance in a Phase-1 attenuation estimator; not really a shared-subspace deflation | production frame, algebraic + Phase-1 regression, no rescoring | blob is **<1% of table energy**, moving it changes nothing; the real finding was that donor attenuation weighting itself doesn't move `pds` (K562 share has partial R²=0.978, attenuation direction ≈ 0.05) |
| 209 | the Orion-**donor**-pair's shared excess subspace (rank 2 = the "long-gene axis"), pre-whitening | 162's PLOLO, actual production a23 render (not a cosine-only toy), 4 folds | **NO-GO**: every real deflation arm is cosine-negative (−1.2%…−9.2%), a **permutation control is the maximizer**; production's whitener already removes ~91% of this subspace; "~0.1% of delivered energy left to take, and taking it costs more than it's worth" |
| 246 | (a) 099/218's separate H1-union "self-residual" object's own shared subspace, tested on top of the shipped stack; (b) the shipped **CD4 residual's own** top-k shared subspace (one un-merged source, before it is added into the final design) | 162's PLOLO, actual production render, 4 folds | (a) **NO-GO** — cosine gain vanishes once composed with the stack ("the stack already carries the un-whitened relaxation"); (b) **BUILD-WORTHY** — CD4-alone deflation (k=5) buys `pds` +0.012 raw for −1.2% aggregate cosine, priced +0.0023…+0.0082 overall |
| 249 | built (b) above into production: `VCC_CD4_DEFLATE_K=5` in `scripts/make_submission_pm_cd4.py`, applied to the CD4 residual **before** it is sized/merged into the final design | full build pipeline, byte-level KA against 246 | **BUILT AND SHIPPED.** Confirmed directly: `experiments/323-current-family-validation/work/compose.json` → `arms.cap15.cd4.deflation = {k: 5, n_carried: 80, energy_removed_frac: 0.3206, top_k_energy_frac_of_centred: 0.3202}`. So the CD4 component's own top-5 deflation is **already inside `cap15`** — CLOSED, not to be re-proposed. |
| 345 | **measurement only, no deflation performed.** SVD energy shares of the reference (`R`) and of `cap15`'s design **restricted to the reference's gate support** (`Dsup`, 259 targets × 4,951 genes) | descriptive, H1 | `Dsup` rank-1 alone carries **14.01%** of `Dsup`'s own energy (rank-2 26.22%, rank-3 37.74%, rank-5 48.70%) — i.e. **even after CD4's own deflation is baked in**, the fully-merged, whitened, preserved `cap15` design still carries a large shared component on the support that scores. 345 only tried *adding* an external common vector (donor means, LOO truth) to the design; it never tried *removing* the design's own top-k directions and rescoring. |

**Conclusion of the record check: NOT closed on the object this mandate asks about.** Every prior
deflation acted on a donor table, the render's correction term, or one un-merged residual source
(CD4 alone) *before* it enters the final design. Nobody has taken the SVD of the actual delivered
`effects.npz['cap15']` array (85 × 18,533, all sources merged, CD4's own deflation already
included, whitened, preserved — i.e. exactly what production emits) and deflated *that* object's
own top-k target-shared directions, then rescored with the actual 6-member counts-based harness.
345 shows there is something there to remove (14% rank-1 on gate support) that survived every
upstream fix. **Proceeding, per the assignment's own instruction for this case.**

## Known-answer reproduction (before any new measurement, rule 3)

Reuse 323's frozen `work/effects.npz['cap15']` and `work/reference.npz` (read-only). Reimplement
the emission+scoring loop from `evaluate323.emit_and_score` / `342/run342.py:emit_arm` **inside
this experiment's own script**, writing only to `experiments/380-.../work/`, and run it once on
`cap15` unmodified (k=0, identity — the `floor`-equivalent). This must reproduce
`experiments/323-current-family-validation/work/metrics-cap15.json`'s `raw` six-member dict,
`panel_total_digest_sha256`, and `emitted_count_digest_sha256` **exactly**. If it does not
reproduce, stop and report that instead of proceeding to new arms.

## Hypothesis

Deflating the top-k target-shared singular directions of the FINAL `cap15` design (all sources,
post-preserve — the literal object production emits) raises `pds` on the 6-member H1 harness,
because 345 shows a real, sizeable (14% at k=1 on gate support) shared component survives there,
and the general mechanism established by 156/175/204/209/246 is that shared/common directions in
an object that reaches the scorer un-deflated cost `pds` (hubness) while often carrying little
of the per-target DE signal. This is a genuinely new object relative to every closed experiment
above (see table).

## Construction (`deflate380.py`)

Inputs: `experiments/323-current-family-validation/work/effects.npz` (`cap15`, `genes`,
`targets`), `experiments/323-current-family-validation/work/reference.npz` (via
`evaluate323.load_reference()`), `experiments/323-current-family-validation/work/h1-anchors.json`
(323/342's frozen H1-local anchors, read-only, used only to compute the "overall" scaled score
exactly as 342 does), `experiments/323-current-family-validation/report323.py`
(`mse_bootstrap`, `ci`, `scaled` — imported, not copied).

For a target×gene matrix `X` (85 × 18,533, the **native** effects.npz axis — SVD is computed on
this object exactly as named in the mandate, before `run323.map_axis` reindexes it onto the
reference's own (smaller) gene axis for emission):

1. **Centering choice, decided at k=1 by actually running both through the 6-member harness**
   (cheap: ~150-250s/arm, well inside budget) rather than guessing from a cheap proxy:
   - `uncentered`: `Xc = X`.
   - `centered`: `Xc = X - mean_over_targets(X)` (per-gene mean across the 85 rows).
2. Economy SVD `U, S, Vt = svd(Xc, full_matrices=False)` (rank ≤ 85, trivial cost).
3. Remove top-k: `Xc_def = Xc - U[:, :k] @ diag(S[:k]) @ Vt[:k]`.
4. If centered, add the mean back (so only the *above-mean* shared pattern is removed, not the
   flat mean itself — reported separately since 156/204 already show the flat mean is <1% of
   energy and is not the object in question here).
5. **Norm-restore per row**: rescale each target row of the deflated matrix back to that row's
   *original* (pre-deflation) L2 norm on the native axis, so amplitude is not a confound (this is
   the mandate's explicit instruction and matches 246/249's `deflate()` convention).
6. `energy_removed_frac = sum(S[:k]**2) / sum(S**2)` (on whichever centering convention is in
   force), reported per k.

Controls, same recipe otherwise:
- `random_k3`: replace the real top-3 `U[:, :3]` with a random orthonormal 3-dim basis in the
  85-row space (QR of a seeded Gaussian 85×3 matrix — the exact convention used in
  `experiments/246-direction-hypotheses/lib246.py`'s `deflate(..., mode="rand")`). A literal
  3-dimensional random subspace of an ≤84-rank row space cannot be made to carry the same energy
  as the true top-3 by chance (expected ≈ 3/84 ≈ 3.6% vs the true top-3's actual share) — this is
  noted explicitly rather than silently forcing a match; both energies are reported side by side.
- `shuffled_k3`: apply ONE fixed random permutation to the **gene axis** (columns) of `Xc`, take
  the top-3 right-singular vectors of the *permuted* matrix (same singular values as the real
  arm's spectrum would have on a relabelled genome), map them back to the true gene coordinates
  via the inverse permutation, and remove *that* (gene-identity-scrambled but shape-legitimate)
  basis from the true (unpermuted) `Xc`. This isolates whether the true top-3's specific gene
  identity (not just "some same-shaped pattern") is what carries the effect.

## Emission + scoring (`run380.py`)

Byte-for-byte the same call sequence as `evaluate323.emit_and_score` / `342/run342.py:emit_arm`
(`b.profile_from_effect`, `b.emission_spread_vector`, `b.emission_rank_axis`,
`rc_emit.rc_hybrid_block` with production's default `emit_counts`, `E.block_stats`, `E.score`),
imported from `run323`/`evaluate323`, **never edited**, with every write redirected to this
experiment's own `work/` (never `experiments/323-.../work/`). One arm per process invocation;
`memory_pressure | tail -1` checked before each (gate: proceed only if free ≥ 50%).

Arms, in order: `cap15` (KA, k=0) → `deflate_k1_centered`, `deflate_k1_uncentered` (decide which
to keep) → `deflate_k2`, `deflate_k3`, `deflate_k5` (chosen centering) → `random_k3`,
`shuffled_k3` (chosen centering). 8 harness runs total.

## Statistics (`report380.py`)

Paired 10,000-target bootstrap, seed 20260909 (323/342's own seed and resampling scheme,
reused via `report323.resample`-equivalent code inlined here since `report323.main()` is a
full-file script, not an importable function — only its `mse_bootstrap`/`ci`/`scaled` helpers are
imported). "H1-local overall" = `scaled(raw, anchors.baseline, anchors.replicate).mean(axis=1)`
exactly as 342 does, using 323's frozen `h1-anchors.json` unchanged. Contrasts: each arm minus
`cap15` (the actual production baseline, not `uniform`/`plain15`), on overall and all six raw
members, with 95% percentile CIs.

## Kill criterion (from the assignment's "Reading" section)

Board-worthy only if, for some k: (a) `pds` rises with a 95% CI excluding 0, **and** (b) the
overall H1-local score rises with a 95% CI excluding 0, **and** (c) `random_k3` does NOT
similarly resolve a positive `pds`/overall move. Otherwise: NEGATIVE. Report the `mse`/DE member
price explicitly regardless of outcome (rule 5), for the PI to weigh against 249's board-transfer
ratios.

## Memory / time

Native SVDs are 85×18,533 — trivial (<1s, <1 GiB). The only ~5 GiB step is the emission+scoring
loop per arm (323/342's own measured cost: ~142 s, ~5 GiB). 8 arms sequential, memory-gated
before each. Budget 3 h; expected wall time under 45 minutes including overhead.

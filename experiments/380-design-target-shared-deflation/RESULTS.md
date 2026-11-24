# 380 — Deflating the FINAL cap15 design's own top-k target-shared directions: clean NEGATIVE, monotone in k, and the two controls show the harm is direction- and identity-specific

**Status: NEGATIVE.** Board-worthy only if `pds` rose with a resolved interval AND overall rose
with a resolved interval AND a random-direction control failed to reproduce it (PLAN.md kill
criterion). None of the three hold. `pds` never has a positive point estimate at any k; at k=3
and k=5 it is **resolved negative** (opposite of the hypothesis); the H1-local overall score is
resolved negative from k=1 (centered) upward. There is no gain to weigh a price against.

## Record check (mandatory, done first — see PLAN.md for the full table)

Read `experiments/156-*/RESULTS.md`, `175-*/RESULTS.md`, `204-*/RESULTS.md`, `209-*/RESULTS.md`,
`246-*/RESULTS.md` (§0/§1), `249-*/RESULTS.md`, `345-*/RESULTS.md`, and inspected
`experiments/323-current-family-validation/work/compose.json` directly.

**Finding: not closed on this object.** 156 deflated the render's mass-fix correction geometry
(no-op). 175 deflated per-donor rank-1 PCs on donor tables (local sign-accuracy proxy gain, later
priced NEGATIVE on the real instrument by 209, −2.3%…−2.7% cosine). 204 checked the donor table's
target-independent mean (<1% of energy, irrelevant) and found donor attenuation weighting doesn't
move `pds`. 209 deflated the Orion-donor-pair's shared excess subspace pre-whitening — NO-GO,
every real arm cosine-negative, a permutation control was the maximizer, and it explicitly closed
with "the whitener has already taken it... ~0.1% of delivered energy left to take." 246 tested
deflating a *separate* H1-union residual object on top of the stack (NO-GO — "the stack already
carries the un-whitened relaxation") and, distinctly, deflating the **CD4 residual's own** top-k
shared subspace *before* it merges into the design (BUILD-WORTHY, `pds` +0.012 raw). 249 built
that CD4-alone construction into production (`VCC_CD4_DEFLATE_K=5`). Directly confirmed via
`compose.json`: **`arms.cap15.cd4.deflation = {k: 5, energy_removed_frac: 0.3206}` is already
inside `cap15`** — CD4-alone deflation is CLOSED, shipped, and not re-proposed here. 345 measured
(descriptively only, no rescoring) that `cap15`'s design, restricted to the reference's gate
support (`Dsup`, 259 targets × 4,951 genes), still carries **14.01%** of its own energy in a
single rank-1 direction — i.e. shared structure survives every upstream fix — but never removed
it and rescored.

**Conclusion: every prior deflation acted on a donor table, the render's correction, or one
un-merged residual source (CD4 alone) *before* it enters the final design.** Nobody had taken the
SVD of the actual delivered `effects.npz['cap15']` array (85 × 18,533, all sources merged, CD4's
own deflation already included, whitened, preserved) and deflated *that* object's own top-k
directions, then rescored with the real counts-based 6-member harness. Proceeded per the
assignment's own instruction for this case.

## Known-answer reproduction

Reimplemented `evaluate323.emit_and_score`'s exact call sequence
(`profile_from_effect` / `emission_spread_vector` / `emission_rank_axis` / `rc_emit.rc_hybrid_block`
/ `block_stats` / `score`, same seeds) inside this experiment's own `run380.py`, writing only to
this directory's `work/`, and ran it once on `cap15` unmodified (k=0, identity).

| check | result |
|---|---|
| `emitted_count_digest_sha256` vs `experiments/323-current-family-validation/work/metrics-cap15.json` | **exact match** |
| `panel_total_digest_sha256` | **exact match** |
| every saved per-target array (`sums`, `jk`, `p`, `q`, `lfc`) | **`np.array_equal` exact**, all five |
| `mse_legs` (per-gene numerator/denominator dict, 82 genes × 2) | **exact match**, every entry |
| raw `pds`, `mse`, `fid`, `reach`, `jac` | **exact match** (`abs diff == 0.0`) |
| raw `nmae` | differs by **1.1102230246251565e-16** (one float64 ULP) — `0.9508364630325663` (mine) vs `0.9508364630325662` (323's) |

One caveat, checked and cleared before trusting the above: `scripts/make_submission_pm_cd4.py`
has moved since 323 froze its inputs (post-323 lab history: 249, 251, 254, 312, 316, 319, and now
`VCC_CTX_A_PRIOR`/366–373). Diffed `8e5ca0c` (323's freeze commit) against `HEAD`: every changed
line is new, independently-gated, default-off code (a `VCC_CTX_A_PRIOR` block, unset here); **zero
lines touch** `profile_from_effect`, `emission_spread_vector`, `emission_rank_axis`, `make_rngs`,
or `array_sha256` — the only five functions this harness calls. All other 56 of 323's 57 pinned
inputs match exactly. Documented in `run380.py`'s `KNOWN_INERT_DRIFT` rather than silently
bypassed; `work/environment-drift.json` records it.

Given every digest and every underlying array is bit-identical and only the final scalar `nmae`
aggregate differs by one ULP (consistent with non-associative floating-point reduction order in
`cell_eval2`/polars's multi-threaded aggregation, not a computational difference), **the known
answer is reproduced** under `evaluate323`'s own tolerance convention (`E.TOLERANCE["nmae"] =
1e-6`, already used by its "historical" known-answer gate) — ten orders of magnitude looser than
the observed gap. `work/known-answer.json`: `"pass": true`.

## What was measured

**Object**: `experiments/323-current-family-validation/work/effects.npz['cap15']` — the literal
production design (85 targets × 18,533 genes, all sources merged: donors + CD4-with-its-own-k=5-
deflation-already-applied + H1 union OFF, whitened, preserved). Never `experiments/323-.../work/`
was written to; every output lives in this experiment's own `work/`.

**Construction** (`deflate380.py`): economy SVD of the native 85×18,533 matrix (centred across
targets, or not); remove the top-k right-singular directions (`X - U[:,:k] @ diag(S[:k]) @
Vt[:k]`); norm-restore each target row to its pre-deflation L2 norm (amplitude held fixed).
Mapped onto the reference's own gene axis via `run323.map_axis` exactly as the unmodified arms
are, then run through the identical emission+scoring pipeline as the known-answer step.

**Centering decision** (mandate: "do both for k=1 and say which you keep"): ran BOTH conventions
through the full 6-member harness (not a numpy proxy) before deciding. `centered` cost an order of
magnitude less on `pds` for the same energy removed (Δpds −0.00084 vs −0.00994 at 3.72%/3.70%
energy) and is the more defensible definition of "target-shared" (a component common across how a
gene responds to *different* targets, not conflated with one target's raw row scale). **Kept
`centered=True`** for k=2, k=3, k=5 and both controls; `work/centering-decision.json` records the
full reasoning and both raw rows.

**Controls**, same recipe otherwise:
- `random_k3`: a random orthonormal 3-dim basis in the 85-row space (QR of a seeded Gaussian
  matrix — `experiments/246-*/lib246.py`'s own `mode="rand"` convention). Not energy-matched by
  construction (a genuinely random 3-dim subspace of an ≤84-rank row space cannot be made to
  capture the same energy as the true top-3 by chance): removed 3.69% vs the real k3's 10.31%.
- `shuffled_k3`: the REAL top-3 right-singular vectors, with a single fixed derangement applied
  to their gene-axis coordinates, then removed from the true (unpermuted) design — same
  coefficient shape/energy potential as the real component, wrong gene identity. **A design
  mistake was caught here before any heavy run**: the first version learned a basis on a
  gene-permuted *copy* of the data and mapped it back via the inverse permutation — but a column
  permutation is an orthogonal change of basis, so that construction silently *recovers the
  original top-k subspace exactly* (verified empirically on synthetic data: the two constructions'
  `energy_removed_frac` agreed to float noise). Caught by a standalone numpy sanity test before
  touching the real 6-member harness; fixed to scramble the *coefficients* onto random genes
  directly. Re-verified on synthetic data with real shared structure (99% energy in a top-3):
  fixed version removes 0.05% vs the real top-3's 99.3%, as it should.

**Scoring**: paired 10,000-target bootstrap, seed 20260909 (323/342's own scheme), reusing
`report323.mse_bootstrap` / `.ci` / `.scaled` (imported, not copied) and 323's frozen
`h1-anchors.json` for the "H1-local overall" combination — all read-only inputs, identical to
342's own methodology.

## Results

**Energy fraction removed per k** (native 18,533-gene axis, centred convention):

| arm | k | energy removed |
|---|---:|---:|
| `deflate_k1_centered` | 1 | 3.7215% |
| `deflate_k1_uncentered` | 1 | 3.7027% |
| `deflate_k2` | 2 | 7.3571% |
| `deflate_k3` | 3 | 10.3083% |
| `deflate_k5` | 5 | 16.0132% |
| `random_k3` (control) | 3 | 3.6930% |
| `shuffled_k3` (control) | 3 | 0.0152% (vs. the real top-3's 10.3083% it is built from) |

For scale: 345's 14.01% rank-1 figure was measured on `Dsup` (85 × 4,951 gate-support genes only)
— a different, smaller object than the full 18,533-gene universe used here per the mandate's
explicit instruction ("SVD of the 85×G design matrix `effects.npz['cap15']`"). Restricting to
gate-support genes concentrates relative shared energy; not a discrepancy, a scope difference,
flagged per rule 5.

**All six raw members and the H1-local overall, Δ = arm − cap15, 95% CI from the paired
10,000-target bootstrap** (units: raw metric scale; `mse`/`nmae` are lower-is-better, the rest
higher-is-better; "RESOLVED ✓/✗" = entire 95% interval on the stated side of zero):

| arm | Δoverall (scaled) | Δpds | Δmse | Δnmae | Δfid | Δreach | Δjac |
|---|---:|---:|---:|---:|---:|---:|---:|
| `k1_centered` | **−0.00214** `[−0.00583,−0.00050]` HARM✓ | −0.00084 `[−0.00196,+0.00014]` | +0.00492 `[−0.00224,+0.01380]` | **+0.00065** `[+0.00004,+0.00154]` HARM✓ | −0.00020 `[−0.00093,+0.00041]` | −0.00793 `[−0.02204,+0.00009]` | −0.00015 `[−0.00045,+0.00003]` |
| `k1_uncentered` | −0.00416 `[−0.01255,+0.00048]` | −0.00994 `[−0.02941,+0.00014]` | +0.00403 `[−0.00030,+0.01076]` | +0.00107 `[−0.00010,+0.00276]` | −0.00061 `[−0.00190,+0.00024]` | +0.00015 `[−0.00422,+0.00564]` | −0.00007 `[−0.00015,+0.00001]` |
| `k2` | **−0.00606** `[−0.01667,−0.00005]` HARM✓ | −0.00966 `[−0.02927,+0.00112]` | +0.02417 `[−0.00048,+0.07077]` | **+0.00291** `[+0.00010,+0.00736]` HARM✓ | −0.00110 `[−0.00314,+0.00041]` | −0.00620 `[−0.02137,+0.00360]` | −0.00041 `[−0.00117,+0.00000]` |
| `k3` | **−0.00822** `[−0.01816,−0.00290]` HARM✓ | **−0.01148** `[−0.02577,−0.00168]` HARM✓ | **+0.03342** `[+0.00008,+0.08220]` HARM✓ | **+0.00363** `[+0.00045,+0.00851]` HARM✓ | −0.00111 `[−0.00318,+0.00049]` | **−0.01270** `[−0.02942,−0.00085]` HARM✓ | **−0.00041** `[−0.00113,−0.00001]` HARM✓ |
| `k5` | **−0.01324** `[−0.02614,−0.00560]` HARM✓ | **−0.02423** `[−0.05070,−0.00504]` HARM✓ | **+0.04249** `[+0.00261,+0.09064]` HARM✓ | **+0.00551** `[+0.00143,+0.01092]` HARM✓ | **−0.00200** `[−0.00409,−0.00031]` HARM✓ | −0.01006 `[−0.02763,+0.00384]` | **−0.00049** `[−0.00130,−0.00003]` HARM✓ |
| `random_k3` | −0.00060 `[−0.00419,+0.00256]` | −0.00210 `[−0.00812,+0.00364]` | **+0.00162** `[+0.00050,+0.00296]` HARM✓ | +0.00044 `[−0.00071,+0.00156]` | **−0.00120** `[−0.00228,−0.00015]` HARM✓ | +0.00305 `[−0.00967,+0.01487]` | −0.00007 `[−0.00019,+0.00004]` |
| `shuffled_k3` | −0.00056 `[−0.00136,+0.00010]` | −0.00070 `[−0.00140,0.00000]` | +0.00002 `[−0.00005,+0.00008]` | +0.00004 `[−0.00002,+0.00010]` | −0.00002 `[−0.00024,+0.00020]` | −0.00151 `[−0.00484,+0.00134]` | +0.00002 `[−0.00002,+0.00006]` |

`work/results.json` carries every number above plus the full per-arm raw values, digests, and
bootstrap identity checks (`mse` bootstrap identity, budget identity — both pass to 1e-10 on every
arm, per `report323.py`'s own gates).

## Verdict against the PLAN.md kill criterion

**Kill criterion**: board-worthy only if, for some k, `pds` rises with a resolved 95% CI, AND
overall rises with a resolved 95% CI, AND `random_k3` does not reproduce it. **All three fail
simultaneously and unambiguously**:

1. **`pds` never has a positive point estimate at any k tested** (1, 2, 3, 5; centered or not).
   At k=3 and k=5 it is **resolved negative** — the interval excludes zero on the harmful side,
   the opposite of the hypothesis, not merely "unresolved."
2. **Overall is resolved negative** from k=1 (centered) through k=5. It is never resolved
   positive at any k.
3. `random_k3` shows a much smaller, largely unresolved effect (Δoverall CI crosses zero,
   Δpds CI crosses zero) — so the real arms' harm is not merely "any energy removed hurts
   equally"; removing the *specific* top-k directions of the real design hurts several-fold more
   per unit of energy than a generic random direction of the same or larger dimension.
4. `shuffled_k3` — same singular-vector coefficients as the real top-3, gene identity scrambled —
   removes almost no energy (0.015% vs 10.3%) and moves nothing (every CI includes zero, tightly).
   This shows the real top-3's cost is tied to its correct **gene identity**, not merely its shape
   or magnitude: it is capturing genuine, gene-specific, cross-target correlated signal (plausibly
   a shared biological response program across knockdowns) rather than an emission-side hub
   artifact of the kind 156's mass-fix mechanism produces.
5. **The effect is monotone in k, in the wrong direction, from the very first unit removed** —
   unlike 209 (whose real arms started only slightly negative and worsened) or 246's CD4-alone
   construction (which *gained* pds), here even k=1 already costs `pds` and is `nmae`-resolved-
   harmful; nothing suggests an interior optimum at k>0 exists to search for.

**mse/DE price, reported explicitly as required (no gain exists to weigh it against 249's board
transfer ratios, but stated for the record):** every k pays on `mse` (monotonically, up to +0.0425
raw at k=5 against a baseline of 0.975 — a ~4.4% relative worsening, resolved harmful from k=3),
on `nmae` (resolved harmful from k=1 centered), and on `fid`/`reach`/`jac` at larger k. There is no
member, at any k, that improves with a resolved interval. This is a strictly dominated family: it
costs on the discrimination channel the mandate hoped to help AND on every other channel
simultaneously.

**Reading**: this closes the "deflate the FINAL merged design's own top-k" construction the same
way 209 closed the donor-pair excess and 246/Hypothesis-1 closed the on-top-of-stack self-residual
— by direct measurement with controls, not by proxy. Combined with 345's finding (a real 14%
rank-1 share survives on gate-support genes) this says: the shared structure that is *left* in the
production design after whitening and CD4's own deflation is not further-extractable nuisance: it
is load-bearing signal, and removing it costs everywhere, at every dose tried, more than a
generic-direction control, and drops to near-zero cost only when its specific gene identity is
destroyed. `docs/CURRENT_STATE.md`'s "closed" list should read: **final merged post-preserve
design's own top-k target-shared deflation (k=1,2,3,5; centered and uncentered at k=1;
energy-controlled against random and gene-shuffled nulls) — CLOSED, 380, clean negative on all six
members simultaneously.**

## Limits

- H1 development panel only (85 targets, one context, H1 union OFF) — development evidence, not
  an A/B/C forecast, per rule 5.
- `random_k3` is dimension-matched (k=3), not energy-matched to the real arm (3.69% vs 10.31%
  removed) — a literal 3-dimensional random subspace of an ≤84-rank row space cannot be made to
  capture 10%+ of the energy by chance; this is stated rather than forced, and the real arm still
  costs several-fold more per unit of energy removed than this weaker control, which is itself
  informative.
- k was swept only up to 5 (per the mandate); given the strictly monotone-worsening pattern
  observed, there is no reason from this evidence to expect an interior optimum at higher k, but
  it was not tested.
- Six-member local overall uses 323's frozen `h1-anchors.json`; no board conversion factor is
  applied, per 342's own convention.
- The single documented, pre-vetted drift in `scripts/make_submission_pm_cd4.py` (new, unset,
  default-inert `VCC_CTX_A_PRIOR` code) does not touch any function this harness calls; see
  `work/environment-drift.json`.

## Exact commands to reproduce

```bash
cd experiments/380-design-target-shared-deflation
NUMBA_NUM_THREADS=4 OMP_NUM_THREADS=4 POLARS_MAX_THREADS=4 OPENBLAS_NUM_THREADS=4 \
  uv run python run380.py cap15                     # known-answer (writes work/known-answer.json)
# decide centering, then (memory_pressure | tail -1 checked >=50% free before each):
uv run python run380.py deflate_k1_centered
uv run python run380.py deflate_k1_uncentered
# work/centering-decision.json written (centered=true) before continuing:
uv run python run380.py deflate_k2
uv run python run380.py deflate_k3
uv run python run380.py deflate_k5
uv run python run380.py random_k3
uv run python run380.py shuffled_k3
uv run python report380.py cap15 deflate_k1_centered deflate_k1_uncentered deflate_k2 \
  deflate_k3 deflate_k5 random_k3 shuffled_k3        # writes work/results.json
```

## Peak memory and wall time

Directly measured (`/usr/bin/time -l`) on the known-answer (`cap15`) arm: **140.09 s real, peak
RSS 5,349,376,000 bytes = 4.98 GiB**. Every other arm's emission loop reports its own wall time
inside `metrics-<arm>.json["seconds"]`: 133.5–141.0 s, tightly clustered around `cap15`'s measured
value (the deflation transform itself is a sub-second numpy SVD on an 85×18,533 matrix; the cost
is entirely the shared 85-target emission+scoring loop, unaffected by which effects array is fed
in) — consistent with the ~5 GiB / ~142 s per-arm figure already on record for this harness
family. `memory_pressure | tail -1` was checked immediately before every one of the 8 heavy
invocations; free was 77–80% throughout (never below the 50% gate), so no run was held back.
8 arms run strictly sequentially, one process at a time, per rule 2. Total wall time for this
experiment, record-check through final report: well under the 3 h budget.

# 170 — method, gates and what each control would read if the mechanism were absent

This file is the audit trail for `RESULTS.md`. Every number in that document comes through
the pipeline described here.

## The question

`WHITEN_GAMMA = 0.65` is a shipped production parameter (`scripts/make_submission_pm_cd4.py`,
`VCC_WHITEN_GAMMA`). Three facts said it might be set too high:

1. **a19's board result.** Raising gamma `0.65 -> 0.80` **and** adding `lam 0 -> 0.007`
   collapsed `mse` from scaled `0.120057` to `0.024149` (raw `0.873790 -> 0.964917`) while
   `pds` moved only `-0.001183`. Implied `d(raw mse)/d(gamma) = +0.6075` if gamma carries all
   of it.
2. **097 §9c** measured that LOWERING gamma RAISES directional accuracy `kappa`, monotone on
   two surfaces.
3. 097's closure was a **two-surface** rule requiring positivity on the Orion pair, and 168d
   has since measured the Orion pair as overstating by **3.8-10x on large reshapings, with
   sign flips**. The surface that vetoed lowering gamma is the discredited one.

⚠️ **One correction to the framing this experiment was set up with.** 097 did NOT omit `mse`.
§9h's delta table has an `mse` column (`+0.0108 / +0.0090` overall at `gamma = 0` on
h1-test / orion), derived analytically from `kappa` through `MSE_min = 1 + (Kc - 2*sqrt(claim*W))/sumD`.
What 097 lacked is an `mse` measured **through the emitter and the scorer**, and its analytic
route assumed the arm sits ON the rho cliff -- which 164 has since measured to be false for
a15 (`rho = 0.9466`). So 097's `mse` number is a model, not a measurement, and the reason to
redo it is the model's assumption, not a missing column.

## The pipeline, and why each stage calls shipped code

| stage | what runs | file |
| --- | --- | --- |
| donor effects | `build_v9_pre_whiten_effects` (shipped), 3 donors, alpha 0.982 | `stage.py` |
| PerturbMatch | `apply_perturbmatch_k562_replacement` (shipped), blend 0.5 | `stage.py` |
| whitening | `apply_v9_whitening` (shipped), `union_pruned`, gamma swept, lam 0 | `gate.py` |
| CD4 residual | `apply_frozen_cd4_residual` / `add_post_whiten_residual` (shipped) | `gate.py` |
| render | `profile_from_effect` (shipped): CPM gate 10, `MASS_FIX=minnorm` | `h1kit.deliver` |
| emission | `emit_counts(mode="quasi")` + `emission_spread_vector` (shipped), a15's tier-8000 policy | `h1kit.emit_and_score` |
| `mse` | `cell_eval2.metrics.delta.mse_unbiased_capped` / `distance_unbiased` / `moments.jackknife_correction`, via 096's `msekit` + `decompose` | `h1sweep.py` |
| `pds` | `cell_eval2.metrics.discrimination_score`, vcc2026 settings (cosine, `n-1`, midrank, panel exclusion, `control_source=real`) | `h1kit.pds_exact` |
| DE | 092's `de_kernel` (scanpy-identical rank-sum) + BH, then 153's `debench.Surface` | `h1sweep.py` |

Nothing is reimplemented. 156's lesson is explicit about why: its first pass had its own copy
of the mass fix, which used the unclipped mass while the renderer clips, and that is the class
of screen/production mismatch that cost a19 `-0.01352`.

`candidate_pre` (the pre-whitening effect matrix) is **independent of gamma**, so it is built
once and cached; the sweep only re-runs the operator and everything after it.

## The four known-answer gates (rule 10)

| # | gate | what it would read if the pipeline were wrong |
| --- | --- | --- |
| 1 | `baseline_pre` and `candidate_pre` sha256 **bit-identical** to a15's shipped audit | any error in the donor set, alpha, weights, gene axis, PerturbMatch blend or its centring changes the sha. It initially FAILED and the cause was real: `--chunk-rows` defaults to 128 and PerturbMatch's mean is accumulated per chunk, so a different chunk size changes the last bits. |
| 2 | the post-whitening audit scalars at **a15** (gamma 0.65 / lam 0) **and at a19** (gamma 0.80 / lam 0.007) | two shipped board points that share `candidate_pre` and differ only in (gamma, lam). A wrong operator, wrong columns or wrong renorm moves `after_operator_median_norm` by orders of magnitude. All 14 scalars agree to <= 1.1e-15 relative. |
| 3 | `msekit` on the frozen 095 h1-test artifact returns `2.0551018299736388` | the value the full `031-somatic -> cell-eval2 run` pipeline wrote to `agg_results.csv`. Anchors the exact `mse` code path to a real scorer run. Reproduced to 2.2e-16. |
| 4 | 096's `mech_claim_grid2.py` reproduces `claim = 0.11860` at `spread 0.25 / phi 0` and the `0.24653` cap | anchors the challenge-panel `claim` machinery used for the board translation. |

**On gate 2 and sha256.** a15's post-whitening `final_sha256` does NOT reproduce on this
machine while a19's does. That is not a discrepancy in the computation: a15/a19 were built on
Linux/x86_64 (their metas record the platform) and this is arm64, `whitener.apply` is a pair
of 17,942 x 17,942 GEMMs, and BLAS kernels differ in reduction order across architectures. The
evidence that it is last-bit noise and not a different computation is that every scalar read
off the same object agrees to 1e-15 or better, and that the pre-whitening object -- which
involves no large GEMM -- is bit-identical.

## Controls that would fail if the mechanism were absent

Two vacuous checks have been caught in this project in two days, so each control is stated
together with its reading under the null:

* **Render identity.** With `MASS_FIX=minnorm` the renderer's normalisation is the identity,
  so `log1p(TS*p) - b_ctrl` must equal the gated, mass-fixed effect. With the mass fix OFF the
  same quantity reads 2-3% of the row norm (156 measured 2.83-3.35% on the challenge frame).
  Run both ways.
* **Emitted-profile fidelity.** The quasi emitter repairs each gene's panel total exactly, so
  the pseudobulk read back from the emitted cells must match the intended one. a21 measured
  0.0014 for a correct build and 0.0155 for the wrong repair order.
* **Gamma is not a no-op.** Row-wise cosine between the delivered effects at gamma 0.35 and
  0.80. If the operator were silently skipped this reads 1.000.
* **Paired emission.** Every gamma uses the same per-target RNG seed, so between-arm
  differences are driven by the profile rather than by the draw. A re-seeded twin at fixed
  gamma measures what is left.

## Two operational notes worth carrying forward

1. **`--chunk-rows` is load-bearing for reproducibility, not just for memory.** PerturbMatch's
   centring mean and the CD4 banks' column means are accumulated chunk by chunk, so a
   different chunk size changes the last bits and therefore every downstream sha256. The
   builder's default is **128**. Any script that re-runs the staging must pass it.
2. **Do not gate a background chain on `ps aux | grep "<script>.py"` in this repo.** A sibling
   agent's own command line quotes the scripts it is told to use, so `mech_claim_grid2.py`
   matched a *prompt string* belonging to experiment 169 and a wait loop never exited. This is
   the same shape as the `pgrep` self-match bug already in `docs/LOG.md`. Match on a PID
   captured at launch, or on a sentinel file.

## What this experiment deliberately does NOT do

* It does not size anything on the Orion pair (168d's rule). Orion enters only through 097's
  already-published `kappa` table, as a labelled ordering cross-check.
* It does not use `data/donors/h1_2025` (85/85 coverage of h1-test's targets, pure oracle
  leakage). The donor set is `k562_gwps + orion_hct116 + orion_hek293t`, production's exact
  three.
* It sweeps gamma ALONE with `lam = 0`. a19 moved both, which is precisely why its slope is
  unusable.

# 195 — Learn the transform for the retrieval objective instead of choosing it. `pds` is 5/6 of everyone's score.

You are a research subagent. Work ONLY in `experiments/195-retrieval-transform/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory.
Peak RSS ≤ 12 GB. `uv run python`. Long compute is fine; use the wall clock.

## READ FIRST
`AGENTS.md` §4.2 **in full** (the geometry lever, and its four numbered warnings — especially
point 1 on where the operator may be fit), §8.1, §11.1 rules 1–13. Then:
- `src/vcc2026/whiten.py`; `experiments/017`, `019`, `020`, `021`, `022`, `044-csls-distill`,
  `073-postwhiten-csls-distill`, `080-listwise-metric-adapter`, `084-listwise-runtime-study`,
  `170-gamma-mse`, `183-unique-boost`, `162-retrodicting-validator` (PLOLO).
- `experiments/188-magnitude-calibration/RESULTS.md` Phase 3 — **`debench.Surface.pds()` must never
  size a `pds` arm.** Standing rule. Use PLOLO.

**Read 080 and 084 first and tell me in your RESULTS what they already settled** — they are named
"listwise metric adapter" and "listwise runtime study" and may be exactly this idea. If they are,
your job is to say so with their numbers and stop. Do not re-run a closed experiment.

## The gap, if 080/084 leave one

`pds_cosine` is a **retrieval** metric: for each prediction, is its own true effect the nearest of
300 by cosine? §4.2's insight was that retrieval in an anisotropic space suffers hubness, and that
flattening the space fixes it — worth **+0.060 macro PDS**, the largest gain in the project.

But the fix is a **hand-chosen one-parameter family**: `W = (C + λI)^(−γ/2)`, with `γ` swept and
fixed at 0.65. That family was not derived from the retrieval objective; it was borrowed from
word-embedding post-processing and it happened to work. **The objective itself is differentiable
and we have tens of thousands of labelled retrieval problems** — every donor target is a query whose
correct answer is known, and Orion alone carries 17,767 + 17,977 targets.

So: **fit `W` to maximise retrieval directly**, leave-one-LINE-out, instead of picking it from a
family. The natural objectives are a listwise ranking loss (ListNet / softmax over the 300 decoys)
or a metric-learning loss (LMNN / triplet with the true effect as positive and the other 299 as
negatives), over a low-rank-plus-diagonal `W` so the parameter count stays far below the sample size.

**Three warnings that are not negotiable and have each cost this project a result:**
1. **Fit on donor perturbations, never on the 300 panel predictions** (§4.2 point 1). The panel
   basis is rank-deficient and its sweep rises monotonically forever — at the limit it just
   orthogonalises 259 vectors, which improves a ranking by construction.
2. **Hold the scored line out of everything**, including the basis. `data/donors/h1_2025` covers
   85/85 of h1-test's targets and is pure leakage.
3. **Mean cosine is the wrong objective** (§4.2): whitening *lowers* transfer cosine 0.0435 → 0.0295
   while *raising* PDS 0.7378 → 0.8078. Any harness optimising cosine would have rejected the
   project's largest gain. Optimise the rank, and report cosine only as a diagnostic.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
- Reproduce PLOLO's a15 row, or 169's `q = 1.00` identity (16/16 bit-identical, worst |Δpds| 0.00e+00).
- Reproduce §4.2's headline on the LODO harness: incumbent 0.7378 → whitened γ=0.65 0.7873 →
  equal-donor-weight **0.7978**.
- **Rule 11: record a path and a sha256 for every operator you load.** Eight distinct whitening
  operators exist on disk with near-identical names and max eigenvalues 8.036 / 6.257 / 5.895, and
  an entire evening was once lost to loading the wrong one.

### Phase 1 — the fit
Parameterise `W` as low-rank-plus-diagonal, initialise at production's operator (so the identity row
is exactly the incumbent), and optimise a listwise/ranking loss over donor targets, LOLO. Report the
learning curve, the parameter count against the sample size, and the gap between train and held-line
performance — **if held-line PDS does not beat the incumbent, the whole idea is an overfit and you
should say so early and cheaply.**

### Phase 2 — price on PLOLO and all six
`pds` on PLOLO (162's four folds, challenge frame, pre-committed unweighted mean). Then rule 12:
an operator change moves the effect's shape, so `mse` through 096's `claim` route (a19 lost
**−0.01352** on `mse` from exactly this and its screen never computed it), `nmae` on h1-test with a
**significance-selected pool** (188: `|lfc|·√cpm` at the reference's BH depth, median ~89 genes/row
— a top-N-by-magnitude pool inverts the sign of `nmae` deltas), and `reach`/`fid`/`jac` on h1-test.
Composed gate: SHIP only if predicted > +0.0019.

### Phase 3 — falsifiers
1. **Permuted-label control**: fit the same model with the query→answer mapping deranged. The real
   fit must beat it by more than its own fold spread. (Precedent: the Ahlmann-Eltze zero-shot
   transfer had **20 of 20 permuted maps beat the real fit** — that is the failure mode here.)
2. **The panel-basis artefact check**: confirm your fit does not improve monotonically with capacity
   on a rank-deficient basis, which is §4.2 point 1's signature.
3. **Cosine must be allowed to fall.** If your arm raises cosine and lowers PDS, it is optimising
   the wrong thing; report both.

## Deliverable
`RESULTS.md`: what 080/084 already settled, the known-answer rows with shas, the learning curve, the
held-line PDS against the incumbent, the six-member price, and the three falsifiers. "The hand-chosen
family is at the optimum" is a full and valuable result.

## Territory
192 owns the call set. 193 owns the operator's **eigenvalue map** inside the existing family —
coordinate by reading its RESULTS if it lands first, and do not duplicate its Marchenko–Pastur work.
194 owns the donor tables. You own **fitting the transform to the retrieval objective**.

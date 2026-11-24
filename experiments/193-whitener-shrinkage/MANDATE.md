# 193 — The whitening operator's regulariser is a ridge. It should not be. `pds`, where we are 3rd of 612.

You are a research subagent. Work ONLY in `experiments/193-whitener-shrinkage/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory.
Peak RSS ≤ 14 GB. `uv run python`.

## READ FIRST
`AGENTS.md` §4.2 (the geometry lever — this is the single largest measured gain in the project),
§8.1, §11.1 rules 1–13. Then:
- `src/vcc2026/whiten.py`, `scripts/build_union_whitener.py`.
- `experiments/017`, `019`, `020`, `021`, `022` (the original whitening result and its graveyard),
  `091-basis-prune`, `145-spectral-trace-bins`, `170-gamma-mse` (gamma is a six-member local
  optimum — **do not re-sweep gamma**), `183-unique-boost`, `162-retrodicting-validator` (PLOLO).
- `experiments/188-magnitude-calibration/RESULTS.md` §Phase 3 — **`debench.Surface.pds()` must
  never size a `pds` arm; use PLOLO.** This is now a standing rule.

## The gap

Production whitens with `P → P (C + λI)^(−γ/2)`, `C = BᵀB/n` over donor effect vectors, `γ = 0.65`,
`λ` a scalar ridge. That is worth **+0.060 macro PDS** — the largest gain this project has recorded.
Every subsequent attempt to improve the geometry has moved `γ` or `λ`, or added a term.

**Nobody has questioned the regulariser itself.** `C + λI` is the crudest possible covariance
estimator: it shrinks every eigenvalue by the same additive constant. The whole point of whitening
here is to flatten an anisotropic spectrum for a **retrieval** problem, and the sample covariance's
eigenvalues are biased in a known, direction-dependent way — large ones up, small ones down — which
is exactly the error a uniform ridge cannot correct. The statistics literature has better estimators
for precisely this, and they are cheap:

- **Ledoit–Wolf linear shrinkage** and, better, **analytical NONLINEAR shrinkage** (Ledoit & Wolf
  2020), which shrinks each eigenvalue by its own optimal amount.
- **Eigenvalue clipping / Marchenko–Pastur denoising** (random-matrix theory): eigenvalues inside
  the MP bulk are noise and should be replaced by their mean, not merely damped.
- **Rotationally-invariant estimators (RIE)**, the optimal-in-Frobenius family.
- Shrinking toward a **structured** target (diagonal, block, or a factor model) rather than `I`.

All of these are drop-in replacements for `(C + λI)^(−γ/2)`: they produce a different eigenvalue
map `f(ev)` in the same `W = V f(ev) Vᵀ` form the codebase already uses.

**Second, independent question: `C` is fit on too little data.** §4.2 says the donor basis is
12,000 × 7,760. Orion alone carries **17,767 (HCT116) + 17,977 (HEK293T)** targets and K562-gwps
9,787. Report what the shipped operator was actually fit on and whether using more of it helps.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
- Reproduce PLOLO's a15 row exactly (162's four donor folds × challenge frames). 169 §2 records
  `q = 1.00` reproducing 162's a23 rows at **worst |Δpds| = 0.00e+00, 16/16 bit-identical**, and the
  a15 → a23 mass-fix price **+0.004014 scaled**. Reproduce at least one of those.
- The operator inventory in `RESULTS.md` §13 lists **eight** distinct whitening operators on disk
  with near-identical names and max eigenvalues 8.036 / 6.257 / 5.895. **Rule 11: record a path and
  a sha256, never a description, and never infer an operator from a filename.** An entire evening
  was once lost to this. Assert the sha of every operator you load.

### Phase 1 — the eigenvalue map
Recompute the donor covariance's spectrum. Report the eigenvalue distribution against the
Marchenko–Pastur prediction for the actual `(n, p)` — how much of the spectrum is bulk (noise) and
how much is signal. Then implement the alternative maps above as `f(ev)` and screen them on
**PLOLO** (`pds`, challenge frame, 162's pre-committed 4-fold unweighted mean), with the ridge at
production's setting as the identity row.

### Phase 2 — the basis
Refit `C` on progressively more donor effect vectors, and on pruned/unpruned variants, holding
`f(ev)` at production's. Report `pds` against the incumbent. **The leakage rule is absolute:** fit
the operator on donor perturbations, never on the 300 panel predictions (§4.2 point 1 — the panel
basis's sweep rises monotonically forever, which is an artefact), and when scoring a held line,
hold it out of the basis (`RESULTS.md` §13 records which operators are leakage-safe).

### Phase 3 — price all six (rule 12)
An operator change moves the effect vector's shape, so **`mse` is not exempt** — a19 lost −0.01352
on `mse` from exactly this and its screen never computed it. Price `mse` through 096's `claim`
route on the challenge panel, `nmae` on h1-test emitted **with a significance-selected pool**
(188: build the pool with `|lfc|·√cpm` at the reference's BH depth, median ~89 genes/row — a
top-N-by-magnitude pool inverts the sign of `nmae` deltas), and `reach`/`fid`/`jac` on h1-test.
Composed gate: SHIP only if predicted > +0.0019 (shape-only).

## Deliverable
`RESULTS.md`: known-answer rows with sha256s, the MP/bulk decomposition, the `f(ev)` screen on
PLOLO, the basis-size sweep, the six-member price, and a verdict. If the ridge is already optimal,
say so — that closes the geometry for good and is worth knowing.

## Territory
192 owns the call set. 194 owns the donor tables. 195 audits. You own the **operator**.

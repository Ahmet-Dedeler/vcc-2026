# 204 — The three donors' achieved repression differs 3.3×, and we weight them equally. Derive the weights instead of fitting them.

You are a research subagent. Work ONLY in `experiments/204-repression-weighted-mixing/`.
No `vcc`, no `cell-eval2`, no submission build, no writes outside your directory.
**Peak RSS ≤ 10 GB** — three other agents and one submission build are running on this 64 GB box.
`uv run python`.

## READ FIRST
`AGENTS.md` §2.1, §4.2 **in full** (equal donor weighting is one of its four numbered findings),
**§8.2** (the entire remaining gap is direction, sized at +0.00365 overall per +0.01 of cosine),
§11.1 rules 1–13. Then:
1. `experiments/194-donor-reestimation/RESULTS.md` — **the measurement this mandate exists for, and
   its explicit hand-off.** Read its Phase 0 known-answer rows and reuse its streaming code.
2. `experiments/019-whiten-interact/`, `010-lodo/`, `039-coverage-weight/`, `022-graveyard/` — where
   equal weighting was established and where donor re-weighting was closed.
3. `experiments/197-cosine-inversion/RESULTS.md` — the exchange rate and the identifiability bound.
4. `experiments/162-retrodicting-validator/` — PLOLO, the only `pds` instrument validated against a
   board point. **`debench.Surface.pds()` must never size a `pds` arm** (188 Phase 3).

## The measurement, and why it reopens a closed question

194 measured, for the first time, each production donor's **achieved repression** — the target
gene's own residual expression after CRISPRi:

| donor | median residual | implied knockdown | share of targets clearing the challenge's own >80% criterion |
|---|---:|---:|---:|
| `k562_gwps` | **0.155** | ~85% | **47.6%** |
| `orion_hct116` | **0.300** | ~70% | — |
| `orion_hek293t` | **0.510** | ~49% | **4.4%** |

**A 3.3× spread, and we average the three with equal weights.** `AGENTS.md` §4.2 established equal
weighting (exp 019) and 022 re-confirmed it under whitening — but **all of that predates this
measurement.** 194's own closing line: *"the equal-weight donor consensus was established before
anyone measured that difference — that belongs to whoever owns donor mixing."* You own it.

**Why attenuation is not a scale nuisance here.** Write donor `i`'s effect for a target as
`d_i = a_i·s + ε_i`, with `s` the shared cross-line component, `a_i` the donor's attenuation, and
`ε_i` line-specific plus sampling. `pds_cosine` is scale-invariant **per row**, so a *global*
rescale of one donor is invisible — but the consensus **averages across donors before the row is
normalised**, so a donor contributing at half amplitude contributes half the signal and its full
noise. The minimum-variance estimator of `s` weights by `a_i/σ_i²`, not uniformly. With
`a ∝ (0.845, 0.700, 0.490)` that is roughly **1.7× more weight on K562 than on HEK293T**, before any
noise term.

**And this is a mechanism, not a fit** — which matters more than the size. `AGENTS.md` §4.2's
closure of donor re-weighting is explicitly about *tuned* shares: *"tuned shares land on opposite
sides of 0.5 and each loses 8% on the other surface — which validates production's uniform
weights."* A weight **derived from each donor's own measured repression** has no free parameter to
overfit and transfers to D/E/F unchanged. That is a different object from a tuned share.

## Work

### Phase 0 — known-answer rows (mandatory, abort on failure)
- Reproduce 194's repression table from the donor cell-level data independently. State exactly how
  residual expression is computed (which cells, which normalisation, self-gene only) — 194's number
  is your target and your definition must match it or you are measuring something else.
- Reproduce PLOLO's a15 row (162's four folds, challenge frame, pre-committed unweighted mean) and
  the a15 → a23 mass-fix price **+0.004014** scaled. 194 reports its scoring path reproduces both to
  **0.00e+00**; reuse it.
- **Rule 11: record a path and a sha256 for every operator you load.** The shipped operator is
  **29,394 × 17,942**, sha **`e2f65ecf`** (193's correction to §4.2's stated 12,000 × 7,760, which
  matches no file on disk), and nine near-identically-named operators exist.

### Phase 1 — estimate `a_i` properly, per donor
Median residual is a summary. Estimate each donor's attenuation as a **scalar per donor** in a way
that is comparable across platforms: the natural route is the shared-target regression slope
`d_i ≈ a_i·(consensus of the others)` over the targets a donor shares with the rest, which cancels
the unknown `s`. Report `a_i` from at least two independent routes (self-gene residual, and the
shared-target slope) and say whether they agree. **They may not**, because the Orion and Replogle
platforms differ in depth and chemistry — and if they disagree, that disagreement is the result and
must not be averaged away.

**One trap, from 194 and from the standing harness warning:** on Orion the knockdown signal
*inverts* (correlation −0.232 / −0.159 against the target's measured residual) and knockdown
weighting **inflates the target-independent blob from 1.00% → 10.45%** on HEK293T. That blob is the
long-gene / nuclear-fraction technical axis. So any `a_i` you estimate from Orion self-gene
residuals may be measuring the platform, not the biology. **Deflate the blob before estimating, and
report `a_i` with and without.**

### Phase 2 — score the derived weights
Score `w_i ∝ a_i` and `w_i ∝ a_i/σ_i²` on **PLOLO** `pds`, plus the **implied transfer cosine**
(197's route, so the answer lands in §8.2's units) and **|t|-weighted `q_w`** (187 §5). Include as
rows: equal weights (the incumbent, your identity row), 194's per-target de-attenuation
(+0.00068 ± 0.00288, so a null), and the tuned shares 019/022 closed.

Also test the **de-attenuate-then-average** form — divide each donor's effect by its own `a_i`
before the equal-weight mean — which is algebraically different from re-weighting and is the form
the model above actually implies.

### Phase 3 — falsifiers, all mandatory
1. **Permuted-attenuation control**: shuffle `a_i` across donors and re-score. 194's own knockdown
   falsifier read a *higher* p(Δ>0) than the real arm, which is exactly the trap here.
2. **Hold the truth fixed.** 194 found that rebuilding the *truth* with the same estimator swings
   `kd_post` from −0.00156 to **+0.00433** — a 0.006 swing from the design choice alone. Score
   against the incumbent truth, always.
3. **Both surfaces, and the leak-free fold.** Report all four PLOLO folds including `kolf2`, which
   162 flags as the only leak-free one, and never quote an absolute from the Orion pair.
4. **Rule 12, six members.** A mixing change moves the effect's shape, so `mse` is **not** exempt
   (a19 lost −0.01352 exactly here) — price it through 096's `claim` route; `nmae` on h1-test with a
   **significance-selected pool** (188). Donors on h1-test: `k562_gwps` + `orion_hek293t` only.
5. Gate: +0.0019 shape-only. Rule 13's `reach` toll does **not** apply to a profile-class arm
   (169 falsified it; the `rc_hybrid` board point confirmed the toll belongs to profile-*changing*
   emission arms).

## Deliverable
`RESULTS.md`: the reproduced repression table, `a_i` from two routes with the blob controlled, the
scored weightings in `pds` / cosine / `q_w`, the six-member price, the falsifiers, and a verdict in
overall units. **"Equal weighting is optimal even given the 3.3× attenuation spread" is a full
result** — and given §4.2 lists equal weighting as one of its four load-bearing findings, it is worth
knowing on a mechanism rather than on a fit.

## Territory
196 emitted depth · 200 shipped-stack ablation · 202 the call-set build. 193 closed the operator's
eigenvalue map and hands you **basis composition** as the live lever there — note it, but do not
pursue it in this run. **You own donor mixing and attenuation.**

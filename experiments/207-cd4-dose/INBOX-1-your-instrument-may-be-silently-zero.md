# INBOX 1 for 207 — h1-test cannot see CD4's DE members at all, and there is a sharper cosine instrument

Coordinator (lead), 2026-09-03 18:35 UTC. Two findings from sibling agents that landed after your
mandate was written. **Verify both yourself; they are claims, not instructions.**

## 1. 🔴 On h1-test, the CD4 residual is IDENTICALLY ZERO — silently

From `experiments/200-production-ablation/RESULTS.md`, its one flagged finding for the builder:

> *"CD4's DE members are unmeasurable on h1-test because production's `solve_cd4_fraction` sizes the
> residual by the **median** row norm, and with 4/85 coverage that median is 0 → `frac = 0.0` →
> `reb_noCD4` is byte-identical to `reb_ship`."*

**CD4 covers only 4 of h1-test's 85 targets.** The shipped sizing rule takes a *median* over the
panel's row norms, so with 81/85 rows at zero the median is zero, the solved fraction is zero, and
**the arm you are sweeping does nothing on that surface** — while every gate still runs and every
number still looks reasonable. This is exactly the failure class `AGENTS.md` §11.1's rule-10
corollary exists for: *"the wrong numbers look exactly as reasonable as the right ones."*

**What to do about it:**
- **Add a non-vacuity assertion to every h1-test row**: the delivered effect at fraction `f` must
  differ from the `f = 0` artifact by a stated norm. If it does not, the row is void — do not report
  it, and say so.
- Your Phase 3 asks for `nmae`/`reach`/`fid`/`jac` on h1-test. **Those members may be unmeasurable
  for this arm on that surface.** If so, say it plainly and price them as unmeasured with the
  mechanism, rather than reporting zeros as evidence of no effect. Rule 12 explicitly permits
  "unmeasured, and here is why it cannot move" — it does not permit reporting a silent zero.
- On the **challenge** the coverage is 276/300, so the median is not zero and the residual is live.
  That asymmetry is the whole point: **the surface that can see this arm is the one with no truth.**
  Your `pds`/cosine instruments (PLOLO, 197's route) run on the challenge frame and are unaffected.
- 200 recommends a guard in the builder. Do not add one — that is a production change and outside
  your mandate. Record it and hand it back.

## 2. A sharper instrument for exactly your primary endpoint

From `experiments/204-repression-weighted-mixing/RESULTS.md`, which spent its run measuring
direction arms on PLOLO and reports as a by-product:

> **Δcos on PLOLO is measurable to ±0.20% relative (≈ ±0.00006 overall)** — sharper than Δ`pds` for
> direction arms.

Your mandate makes cosine the primary endpoint. **Use PLOLO's Δcos as the discriminating statistic
and report its interval**, rather than reading direction through Δ`pds`, which 197 showed is the
member least sensitive to it (whitening buys us 26% more retrieval per unit of cosine than the
field, which is precisely why 062's `pds`-selected dose may be wrong).

At 197's exchange rate of **+0.00365 overall per +0.01 of cosine**, a ±0.20% relative interval on a
0.0792 baseline is ±0.00016 cosine ≈ **±0.00006 overall**. That resolves the entire hypothesis: 206's
H1 range is +0.004 to +0.0133 overall, i.e. **+0.011 to +0.036 of cosine — 70× to 225× the
instrument's noise floor.** If the effect is real you will see it unambiguously; if you measure zero
at that precision, the hypothesis is dead and that is a clean, powerful negative.

## 3. Context you should have, since 204 answered the adjacent question

204 tested repression-derived donor weighting — the other lineage/attenuation idea — and it is a
**decisive NO-GO**: equal weighting is the maximum of the transfer cosine over all 29 alternatives,
its permutation falsifier **ties** the best real arm, and **the negative has power** (the model's own
inputs predicted +8.5% to +16.8% cosine; measured **−1.50%, CI [−2.38%, −0.67%]**).

That does **not** bear on your hypothesis — 204 re-weighted the three *existing* donors, while you
are dosing a *fourth, lineage-matched* one, and 206 §5's argument is precisely that a lineage-matched
donor contributes to the **line-specific** component that 159's `1/k` law treats as noise. But it
does mean the burden of proof is high on this family, and that a permutation control tying your best
arm is the outcome to watch for. 204 also found `a_i` is **not identified** (four independent
estimates of HEK293T's relative attenuation span 0.504→1.254), so do not build anything on an
attenuation scalar.

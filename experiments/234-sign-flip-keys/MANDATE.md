# 234 — Can ANY key we already hold buy the +0.005 of sign accuracy the gate needs? 232 built the instrument; use it.

You are a research subagent. Work in `experiments/234-sign-flip-keys/`. **No `vcc`, no submission
build, no production change. No `cell-eval2` while another agent runs one — check first.** Peak RSS
≤ 10 GB; **check `memory_pressure` and `sysctl vm.swapusage` before any job over ~4 GB and hold if
swap free < 2 GB** — a job that exceeds memory here is killed **without a traceback** (219, 226).
`uv run python`. Three other agents are live (221, 229, 233).

## What 232 just established, and why this experiment is now cheap and sharp

**Read `experiments/232-sign-vote/RESULTS.md` in full first.** It closed sign *voting* and, in doing
so, handed us a calibrated channel:

- **§8.1d confirmed causally.** A **pure sign improvement at bit-identical magnitude** moves `nmae`
  the right way: **Δ`nmae` = −0.002573 per +0.01 |t|-weighted sign accuracy, r = −0.995**, with
  efficiency `e = ΔM/ΔU` of **29–40** against break-even 0.5 (228's magnitude arm: **0.072**).
- **The head-to-head.** At a matched **+4.4%** relative cosine gain: 228's magnitude route gave
  Δ`nmae` **+0.0373** and composite **−0.00702**; the pure sign route gives Δ`nmae` **−0.00889** and
  composite **+0.01170**. **A 0.0187 swing decided entirely by HOW the cosine was bought.**
- **The exchange rate (new, reusable): +0.01 of |t|-weighted sign accuracy = +0.00399 overall**
  (`k_cos` 1.384; +0.00379 at 1.000), r = 0.975 over 14 doses. **So the gate needs
  Δq_w = +0.0048 … +0.0050.**
- **The headroom is enormous and untouched.** Contested mass is **36.27%** of the |t|-weighted pool
  and `sign(mean)` there is a coin flip at **0.52292**, so the ceiling for *any* sign rule is
  **+0.17325** of sign accuracy. **A rule needs 2.75% of that headroom to clear.** The best
  implementable donor-only rule captures **0.52%** — **5.3× short**.
- **A magnitude-preserving sign arm leaves `jac` EXACTLY invariant** (max |Δ`jac`| = 0.000e+00 over
  every arm, fold and depth) — so these arms are cheap to price.

**232's closing diagnosis: "the missing ingredient is information the donors do not have."** **Test
that claim.** It is the single most consequential open statement in the project.

## The question, stated as arithmetic you can check on day one

To clear, you need **Δq_w ≈ +0.005** overall. All of it must come from the **contested** subset
(36.27% of mass, currently **0.52292** accurate), so you need to lift accuracy *there* by
**+0.005 / 0.3627 ≈ +0.0138**, i.e. from **0.523 to ≈ 0.537**.

**That is a very low bar, and it is much lower than the bar the same keys were rejected against.**
220 closed the ordering-key family with a **specification**: a key clears rule 13's *ordering* gate
iff it is a per-pair confidence with accuracy **≈ 0.639** on "is our sign right" (correlation
**≈ 0.30** with correctness). It then measured: **`|H|` = 0.074**, `theta` = 0.002, expression =
0.008. **Those keys were rejected against 0.30. Here the requirement is ~0.074.** Whether that
coincidence is real is exactly what you must determine — **do not assume it.**

## Work

### Phase 0 — known-answer rows (mandatory)
Reproduce, through your own code: **232's exchange rate** (+0.00399 overall per +0.01 `q_w`) and its
oracle response slope (Δ`nmae` −0.002573 per +0.01, r −0.995); **228's `e = 0.072`** for 224's gain;
**220's measured key correlations** (`|H|` 0.074, `theta` 0.002, expression 0.008) with correctness
inside the reference DE pool. **If these do not reproduce, stop.**

### Phase 1 — price every key we hold as a SIGN-FLIP arm
The construction is **flip the sign where the key says we are likely wrong**, magnitudes
**bit-identical**. Sweep the flip fraction as the dose. Keys to test, all already on disk:
`|H|` (220's best), 189's `theta`, expression, `|effect|`/spread ratios, donor disagreement
magnitude, and any per-(target,gene) statistic 220 or 189 built. **Also test their combination** —
a single logistic/rank fit on the contested subset, fitted **leave-one-fold-out**.

**Report, for each, the realised Δq_w and the overall via 232's exchange rate.** The bar is
**+0.0048…+0.0050**.

### Phase 2 — the ceiling under a realistic key
Report the **oracle-restricted-to-this-key-family** number: what a perfectly calibrated version of
the best key would deliver. That separates "the key is weak" from "the family is weak".

## Falsifiers, all mandatory
1. **Sign-shuffle null** — flip the same number of signs at random, matched on |t| mass.
2. **Fit honestly** — any fitted key is fitted **leave-one-fold-out**; report in-fold and out-of-fold
   separately. 224 §7 showed a panel-fitted quantity reads ~0 from power alone; 195 cleared on held
   perturbations and failed at the line boundary.
3. **Magnitude invariance** — prove `|effect|` is bit-identical. If it is not, §8.1d says the arm is
   dead and you have measured the wrong thing.
4. **`kolf2` reported separately** — 232's best rule was **negative on `kolf2`** (−2.11% cosine),
   and `kolf2` is pluripotent while A/B/C are all cancer lines.
5. **All six members** for anything that clears Phase 1. `jac` should be **exactly invariant** by
   232 §11(a) — assert it, and if it is not, your arm is not magnitude-preserving.

## Do not
- **Do not re-open sign VOTING** — 232 closed it, including the arity finding that PLOLO's cancer
  folds leave two donors so every contested pair is a 1-1 tie.
- **Do not use 189's `theta` as a vote weight** — 232 §11(b): it is per-GENE and factors out of
  `sign(Σ_s w(j)·sign(e_s))` **exactly, at every arity**. As a *flip* key it is still testable.
- **Do not convert PLOLO to board with a scalar** beyond 232's own exchange rate, and say which
  `k_cos` you used (the factor is not identified: 1.384 vs 1.000 on two arms).

## Deliverable
`RESULTS.md`: the Phase-0 rows, a table of every key with its Δq_w and overall, the Phase-2 ceiling,
the falsifiers, and a one-line verdict. **"No key we hold reaches +0.005, and the shortfall is X×"
is a complete and valuable result** — it would confirm 232's diagnosis and tell us the remaining
work is *acquiring sign information*, not *combining what we have*. **Propose no slot; build
nothing.**

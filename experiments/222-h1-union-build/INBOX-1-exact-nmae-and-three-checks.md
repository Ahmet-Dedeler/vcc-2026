# INBOX 1 for 222 — price the arm on 216's EXACT `nmae`, not on a rate card, plus three cheap checks

Coordinator (lead), 2026-09-04 02:30 UTC. All four items use artifacts you have already staged.
None of them blocks the build; do them while the emission runs or immediately after.

## 1. 🔴 Run the arm through 216's exact instrument before I submit

`h1union` is currently priced only through 218's **per-cosine rate card**, and the preregistration's
`nmae` prediction is the words *"up or flat"*. **Twice today a proxy got a sign wrong that the exact
identity got right:**
- 216 flipped `f = 0.40` from the proxy's **+0.00025** to the exact **−0.00029**, which is what
  closed the dose;
- my own `cd4f080` extrapolation failed because I scaled `nmae` as a positive member past its
  turning point — the identity in `AGENTS.md` §8.4 says it is optimal at `k = 1` with slope
  **−0.074 below and +1 above**.

**Do not let me ship a third arm on a proxy.** `experiments/216-cd4-hard/six216.py` implements 188's
exact `nmae` identity on the four PLOLO folds together with `pds` and 096/197 `mse`. You have
`work/effects-on028.npy`, `work/effects-off.npy` and rebuilt folds in `work/plolo/`. Run the
`g = 0.28` arm through it and report the exact six-member price. If it disagrees with 218's rate
card, **the exact number wins** and I want it before the slot is spent.

## 2. Concentration over the 77 carried targets — pre-committed reading
220 found both published ordering effects were carried by **one to three perturbations**. 218 never
reported a leave-top-k-targets-out on the union's Δcos or Δ`pds`; `compose218.json` stores per-target
arrays for the **baselines only**. Recompute from your on/off effects through `plolo`.

**The reading is pre-committed in `PREREGISTRATION-ADDENDUM.md` §5:** if the top three targets carry
more than half of Δ`pds` on the cancer folds, a positive board result is **weaker** evidence for
223's stack, because the board would be measuring a handful of perturbations rather than a source.

## 3. Measure the CD4 row-norm ratio — one line, and it is the largest silent October risk
I verified `solve_cd4_fraction` (`scripts/make_submission_pm_cd4.py:1267`) and the docstring is
explicit: the medians are *"taken over every row — including rows that carry no residual."* So below
50% coverage it **zeroes** (218's finding), and **between 50% and ~90% it silently OVER-doses**,
because the zeros pull `median_all_rows` below the carrying-row median and
`solved = target_f / f_per_unit` rises by exactly `median_carrying / median_all_rows`.

On A/B/C this is ≈ a no-op: CD4 carries 276/300, so the all-rows median sits at the ~46th percentile
of carrying rows. **In `docs/DEF-RUNBOOK.md` §5's pessimistic draw CD4's union coverage is 174/300**,
where the all-rows median falls to the **~14th percentile** — and 216 prices `f = 0.55` at −0.0022
and `f = 0.80` at −0.0064 against the argmax. **The board is dark when this would fire and production
raises only on `f_per_unit <= 0`, so it is silent.**

**While CD4 is staged, report `median(‖r_i‖) / quantile_0.14(‖r_i‖)` over its carrying rows.**
If it is ~1.05 the risk is retired; if it is ~1.6 it is the largest October item on the books. **Do
not change any sizing rule for CD4** — that is out of scope for this slot. Just measure it.

## 4. Leak warning your mandate does not carry
**`data/donors/h1_2025` covers 85/85 of h1-test's targets** and it is *inside this arm's union*.
**Neither `h1union` nor any of its components may be priced or validated on h1-test** — that is
oracle leakage, not a weak instrument. PLOLO's challenge frame is the instrument. 223's mandate says
this; yours does not, so it is said here. It also means 221's forthcoming six-member local→board
table applies to this arm **only** through the PLOLO frame.

## 5. One reconciliation for your write-up
218 §7.2 says **10** of the 77 union rows are rows CD4 cannot reach; your staged frame says
`targets_outside_cd4_support = 1`. Those use different definitions (complete-mean3 versus any-vote).
**One sentence naming both definitions** prevents the next auditor re-deriving it.

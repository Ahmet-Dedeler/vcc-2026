# How to do science on this project

Standing instructions for anyone — human or agent — running an experiment here. Every rule
below was paid for with a wasted submission slot or a wrong conclusion. Read this before
proposing anything.

## 1. The board is the only ground truth FOR THE CHALLENGE SCORE. Everything else is a proxy with a known bias — for that score.

> Scoped 2026-09-05 (`docs/RESEARCH-PROGRAMME-2026-09.md`): a public experiment is real evidence about its OWN population
> (the scientific track); whether a result transfers to the challenge is a separate, measured question, never an assumption.

We hold three measurement surfaces and they disagree systematically:

| surface | what it is | its known lie |
| --- | --- | --- |
| counterfactual / analytic | recompute a metric from stored tables | ignores that emission is stochastic; has been wrong by sign |
| harness (`realize.py`) | drive the real emitter, real DE path, real scorer | `n_real` is 92 (orion) or 322 (h1-test) against the challenge's ~2,000 |
| the board | the actual scorer on the actual reference | 2 submissions/day, ~1 hour turnaround |

**Rule: any metric with `n_real` in its denominator or its union is unreadable on our
harnesses.** That is `fid`, `jac`, and `reach`'s normalisation — three of the four DE members.
Only their *ordering-driven* parts transfer. This artefact has cost three predictions
(v14's `fid`, v15's `fid`, a8's `jac`); it will cost a fourth if you forget it.

**Rule: a metric's LEVEL can be unreadable while a FACTOR of it transfers.** The rule above
says `fid`, `jac` and `reach` are unreadable on our harnesses. Read that as *unreadable in
level*, not unmeasurable. Before declaring a member off-limits, ask whether it factorises into
an `n_real`-carrying leg and one that does not:

- `fid = (k/n_pred) * min(1, n_pred/n_real)` — the coverage leg carries `n_real` and is
  unreadable; **the precision leg `k/n_pred` carries none and transfers.** On the challenge
  `n_pred` already exceeds `n_real`, so coverage is saturated and our `fid` of 0.5041 *is* our
  directional precision — we are 0.41% better than a coin flip on the genes we call.
- `expr_mse_unbiased_capped_norm`, written as `(A - 2X + Kc - rho*claim) / sumD`, has four of
  six terms computable from our own submission plus control cells we already hold. Only `X`
  and `sumD` needed the board. That factorisation is what made experiment 096 possible.
- `jac`'s union carries `n_real` and does **not** factorise — which is why only its sign
  should ever be quoted from a harness.

Two of the four DE members factorise and two do not. The rule earns its keep by telling them
apart, not by rescuing all of them.

**Rule: run every emission arm on BOTH harnesses.** Their disagreement is your error bar, and
the *direction* of the disagreement is usually diagnostic — orion's small `n_real` mechanically rewards call-set
cuts, so when orion likes an arm that h1-test dislikes, suspect the cut.

**Correction (2026-08-28): "same policy" is not "same artifact."** This section used to claim
the two harnesses "bracketed the truth for a8." That example does not establish what it was
cited for. `experiments/095-clean-harness/realize.py` consumes a *profile* and drives
`emit_counts` against each harness's own control cells — it never scores the shipped `.vcc`.
So a8 on the harnesses and a8 on the board are the same emission **policy**, not the same
artifact, run through different models against different baselines.

The distinction is not academic. Scoring the shipped `a8-nb-phi4.vcc` directly against
`orion-hct116-thinned` returns `mse = 343` against the board's `1.36`, and the control-profile
gap alone explains **99.8-99.9%** of it: **a prediction artifact is baseline-bound**, so
scoring one against a different cell line measures the distance between cell lines, which here
is **548x** the per-row perturbation signal. No two of our three surfaces share a baseline, so
whether the harnesses genuinely bracket cannot currently be settled either way. Treat
cross-surface agreement as suggestive, and say which of policy or artifact was held fixed.

## 2. Predict before you measure, in writing, with a falsifier.

Every arm gets a pre-registered prediction file: the expected sign and magnitude of each
member, the mechanism you believe produces it, and **the observation that would prove you
wrong**. Write it before the run finishes. `experiments/094-reach-ordering/A7-PREDICTION.md`
is the template.

This is not ceremony. Predicted `+0.0213`, delivered `+0.0041` — the 5x miss is the most
informative number we have about how our harnesses transfer, and it only exists because the
prediction was written down first.

**A single-member objective can point hardest at the worst arm.** This is stronger than
"optimising one member is incomplete", and it is the finding that most changes how to work here.

Tuning the shrinkage family on `expr_mse` alone, every measurement any of three analyses ran
named `gamma = 1.0` the best setting in the family — its `sqrt(W)` of **1.515 / 2.102** is the
highest number anywhere in experiment 096. Scored on all six members, that same arm reads
**−0.0336 on orion**, and its raw `nmae` of **1.0667** sits *above* that member's own baseline
anchor of ~1.0013: **worse on `nmae` than predicting nothing at all.**

So the mse surface did not merely fail to see the problem. It ranked the arms in an order that
was **anti-correlated** with the objective over the range that mattered. A partial objective is
not a noisy version of the full one; it can be an actively misleading one, and the more
confidently it discriminates, the more it misleads.

**Practical rule: no knob in a family ships on a single-member measurement.** Every default in
the shrinkage implementation was set that way; one of the five happened to get a full gate and
was wrong by enough to flip a harness's sign. The remaining four carry explicit
`UNVALIDATED DEFAULT` blocks naming the member each was chosen against — the honest state, and
better than a file that reads as though someone had checked.

**A validation statistic chosen at the wrong scale will pass.** The direction study's
reconstruction of our own prediction matched the shipped artifact at a per-row cosine of
**0.99922 / 0.99799** — which reads as faithful and is not, for the quantity being measured.
The aggregate directional cosine lives in about **2% of the effect vector's norm**, so a 0.2%
angular error is **10% of the signal**, and the reconstructed `kappa` came out at 1.19x the
shipped value on one surface.

Before trusting a check, ask what fraction of the thing you are validating the statistic
actually constrains. A near-perfect agreement on the whole vector says nothing about a quantity
that lives in its residual.

**When a number crosses a codebase boundary, transfer the definition, not the value** — and
where the definition cannot travel, **transfer the ratio within one frame and say so**. That
second clause is what let a verdict survive a 19% level error: the decisive member was measured
directly and everything else was a within-reconstruction ratio. Three failures this session sat
at increasing distance and needed three different catches:

| distance | example | what finds it |
| --- | --- | --- |
| two columns of one table | analytic vs realised multipliers | grep the superseded digits |
| two codebases | `gamma = 0.25` as an exponential and as a power law | ask for the formula written out |
| two pipelines, one quantity | reconstructed vs shipped `kappa` | check against the artifact, and report the ratio even when the cosine looks perfect |
| two panels, one policy | the harness vs challenge rank of a zero cosine | state which population a rank was computed over — **and run the real one before believing a divergence** |

The fourth row earned its last clause the hard way. Two independent readings of the same source
concluded that a gated configuration and a built artifact diverged on six rows. Both reasoned
correctly from the code and from a correct median. Running the actual 300-target panel through
both branches found **zero** rows differing — a floor dominated both, which no amount of reading
would have revealed. **"These two lines differ" and "these two artifacts differ" are separate
claims, and only the second one matters.**

Note this row runs **opposite** to the three above it. Those are all "a number was carried into
a frame where it did not mean the same thing." This one is "two definitions genuinely differ and
the artifact does not care." Same check — go to the real object — protecting against the
opposite conclusion.

**The bound it produced is worth more than the divergence that turned out not to exist:**
`--shrink-pert-missing` is *unobservable* at `--shrink-pert-floor >= 0.18`, because every branch
floors to the same weight. Below 0.18 the policy choice goes live and is unmeasured on every
surface we hold. **So lowering that floor is two changes, not one.**

In all three the fix was **"say which object the number is from"**, never "get a better
number". Each was correct in its own frame and wrong on crossing into another; precision was
never what was missing. The
most dangerous defect of the 2026-08-28 session was the only one that would have reached a
submission, and every guard in this document missed it. Two codebases each had a knob called
`gamma = 0.25`: one an exponential in an absolute rank count, the other a power law in a
normalised rank. Both implementations were internally consistent, no digits were stale, no
tables were crossed — the defect was a shared *name*, which no grep of either file could find.

The tell that would have caught it: quote the **effect**, not the setting. Mean `alpha` — the
actual amount of shrinkage — agreed to within 0.0007 across both functional forms, while the
`gamma` producing it differed by 5x. If two parties quote a knob value, they have agreed on
nothing; if they quote the quantity the knob controls, a mismatch is visible immediately.

**A surface can be structurally incapable of exhibiting the defect you are asking it about.**
This is the sharpest version of the rule below and the one that actually bit. A rank-tie bug in
the shrinkage key annihilated genes by their position in the gene array, on **111 of 900 rows**
of the challenge panel. Three analyses agreed 12/12 on a conclusion computed from that vector,
and when the bug was fixed **every number moved by less than 0.001** — which read as
reassurance and was not. The artefact touched **0 of 132 orion rows and 4 of 85 h1-test rows**.
The harnesses had not survived the fix; they had never been able to show the bug at all.

Before treating a harness result as evidence about a defect, ask whether the harness *contains*
the condition that produces it — not whether it is noisy, whether it is **applicable**.

**Agreement between two analyses is not a check if they could not have disagreed.** In
experiment 096 two independent derivations — one from the metric spec, one pre-registered from
the algebra — converged on the same model of the expression-error member and the same
recommendation to shrink amplitude to ~0.4x. The answer is ~0.94x. Both had dropped the same
two correction terms, each individually larger than the quantity being solved for. The
convergence read as confirmation and was the opposite: **one reading held twice.** Two readers
of the same document should *lower* your confidence relative to one agent who has executed the
code. What broke it was running `cell_eval2`'s own functions and finding a term we had both
assumed was binding was not.

## 3. A null result whose mechanism was never verified is not a result.

We once "closed" per-gene over-dispersion on a clean null. The null was a bug: wrong algebra
assigned φ=0 to 97.9% of genes. **Before believing a null, print the intermediate quantity and
confirm the intervention actually happened.** Distribution, not summary statistic.

## 4. Guard against measuring yourself.

Two real incidents:
- **A harness confound.** `orion-hct116` had control cells 38% deeper than perturbation cells
  in 132/132 perturbations. A rank test reads depth as signal. It silently invalidated seven
  experiments. Every harness now passes `scripts/check_harness_depth.py` before use.
- **An oracle leak.** A ranking key used the *reference's* own lfc — unbuildable at submission
  time. Ask of every key: **could I compute this without the answer?**

## 5. Verify the intervention is bit-exact where it should be, and only where it should be.

The exact-total repair means `pds`/`nmae` must not move under an emission change. When a8
returned `pds` Δ = 0.000000, that was the guard firing correctly. If a quantity you believe is
invariant moves, stop and find out why before reading anything else in the run.

## 6. Understand the metric's shape before optimising it.

`expr_mse_unbiased_capped_norm` is clamped to `[0, 1]`. We spent weeks reading "mse routes are
null" when the truth was that we were 37% inside a **dead zone** where every improvement
scores exactly zero. A clamped metric is a cliff, not a gradient; a plateau in your measurement
may be the clamp, not the mechanism.

Read `docs/reference/vcc2026-metrics-FULL.md` for the member you are moving. Specifically ask:
is it a mean over perturbations or a single panel-wide ratio? Is it clamped? Is it scale-
sensitive or rank-based? Is it normalised by something that differs between harness and board?

## 7. Report faithfully, including against yourself.

Retract in the title, not in a footnote. When you find your own error — the oracle leak, the
inverted-ranking claim, the φ=3-vs-φ=4 seed mismatch — write it up as prominently as the
finding it corrects. Record the seed noise floor alongside every effect, and never compare a
single seed against a mean.

## 8. Cost discipline

AWS workers are `r7i.8xlarge` and stopped by default. Start, build, fetch, **stop**. A full
build is ~20 min via `experiments/094-reach-ordering/aws/fastbuild.sh`. Never leave one
running; never chain a stop into a script that another arm is still using.

**Never `git add -A` while a subagent is editing.** On 2026-08-28 a docs-only commit of mine
used it while an agent was mid-edit on `scripts/make_submission_pm_cd4.py`, capturing a state
where a new argument was referenced in `run()` but not yet defined in the parser. `master`
failed the test suite and was pushed that way. Stage explicit paths, and **run the test suite
before pushing whenever the working tree contains anyone else's work in progress.**

If it happens, the repair that does not destroy the agent's tree is to rewrite the index entry
alone:

```
BLOB=$(git show <good-rev>:<path> | git hash-object -w --stdin)
git update-index --cacheinfo 100644,$BLOB,<path>
git commit -m "..."
```

`git checkout <rev> -- <path>` would overwrite the working file and throw away whatever the
agent had not yet finished.

## 9. What is already closed — do not re-open without new evidence

amplitude (as a global scalar, pre-φ), `min_cpm`, global `spread`, call-set truncation,
`basal_beta`, basis-prune tuning, spread-based ordering, per-gene φ, selective
over-dispersion. Each was closed by measurement; the write-ups say which surface and why.

**But note the asymmetry:** "closed" means *it did not move the members we were then
measuring*. Amplitude was closed against the DE members while `mse` was in its dead zone and
therefore invisible. A lever closed under one objective is not closed under a new one.


## 10. Two failure modes seen repeatedly in this project

**Agreement between two derivations is not a check if they share an omission.** In experiment
096 two independent analyses both modelled the expression-error member as pure geometry
(`raw = r^2 - 2rc + 1`), both reproduced the board's ordering, and both recommended shrinking
amplitude to ~0.4x. Both had dropped the same two terms — `Kc` and `rho*claim` — which are
individually **larger than the quantity being solved for**. The agreement read as
confirmation and was not. Before trusting a concordance, check that the two routes could have
disagreed.

**A competitor's published diagnosis is evidence about their model, not about the task.** A
public writeup reported that DE callability is set by control expression (median CPM 107 for
called genes vs 38.7 for uncalled). Measured on our own submission the ordering is
*reversed* — 48.1 called vs 87.6 uncalled — because our callability is set by our own effect
sizes and spread policy. Their finding is almost certainly true of their pipeline. It is not a
property of the data, and building on it would have sent the tail policy in the wrong
direction.

**External signals transfer; self-referential ones do not.** This is the most reusable thing
experiment 096 produced, and it *predicts* failures rather than only explaining them. Our
directional signal is extremely concentrated — the top quartile of perturbations carries
**79.7% / 85.0%** of all of it across the two harnesses, and **11.8% / 28.8%** point the wrong
way. A confidence score that could identify which quartile is which would be worth a great
deal. Measured correlation with the true per-perturbation cosine, h1-test / orion:

| score | rho with truth | usable? |
| --- | ---: | --- |
| `‖dhat_p‖` — our own predicted magnitude | −0.001 / +0.253 | no, and it flips sign |
| prediction uniqueness — our own | −0.337 / −0.102 | actively harmful |
| **`donor_agreement` — cross-donor** | **+0.209 / +0.242** | **the only one positive on both** |

Both scores computed **from our own prediction** fail and disagree across surfaces; the one
computed from **independent donors agreeing with each other** works on both. The same key is
the only per-gene shrinkage that transferred. Prefer a signal whose evidence comes from
outside the model being scored.

**A per-gene slope is a cell-line property; a per-gene confidence score is a model property.**
Measured transfer across our two harnesses: fitting a per-gene shrinkage coefficient on one
cell line and applying it to the other gains `1.02x` — i.e. nothing, shrink 99% of the way
back to a global scalar. Keying the same shrinkage on `donor_agreement` instead gains
`1.17-1.38x` on both surfaces. Use this to decide what is worth building: quantities that
describe *the model's confidence* transfer across contexts; quantities that describe *a cell
line's response* do not.

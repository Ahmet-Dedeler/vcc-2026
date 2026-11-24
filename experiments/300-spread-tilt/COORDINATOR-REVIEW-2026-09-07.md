# Coordinator review of 300 — 2026-09-07

This is a wording/scope audit of existing results, not a new experiment.
The original results and preregistration remain intact. **The submission
decision is unchanged: no tested arm cleared the frozen gate.**

## What the evidence supports

`RESULTS.md:91–110` reports five feasible positive tilt doses (one post-hoc)
and two direct-spread doses. The primary tilt at d=1 passes its three controls
but fails price (+0.00140 vs +0.00190) and retention (47.8% vs 50%). Larger
measured tilt doses fail their applicable amplitude-changing gate, and their
retention is lower. This is sufficient for the recorded NO-GO.

## Corrections to carry into summaries

1. **A measured grid is not a structural impossibility proof.**
   `RESULTS.md:163–166` itself records non-monotone delivered amplitude
   (+0.469%, +0.233%, +0.528%, +0.867%, +0.751%). Five points do not establish
   that every intermediate dose must fail, and the threshold change at 0.5%
   is a decision rule, not a biological law. The supported conclusion is:
   **none of the tested dose constructions clears these frozen criteria**.
   This does not recommend searching finer doses or weakening the criteria.

2. **Do not compare raw PDS directly with an overall-score standard deviation.**
   `RESULTS.md:295–296` compares +0.00150 raw PDS with 0.004502 overall.
   These are different units. The full six-member forecast for that dose is
   +0.00296 overall (`RESULTS.md:98`), with its own transfer assumptions.
   Moreover, `AGENTS.md:1436–1438` describes unexplained variation across
   different emission arms with identical panel totals, not repeat-run noise
   from submitting an identical artifact. `READ-FIRST.md` records repeat
   scoring as deterministic. Call this unmodelled ordering-channel variation,
   not intrinsic scorer randomness; do not convert it into a sampling sigma
   for this operator experiment.

3. **The completed controls belong to the primary dose.**
   All three controls passed for d=1 (`RESULTS.md:124–133`). The missing
   mirror ladder and spread mirrors are declared at `RESULTS.md:332–335`.
   Do not summarize this as every dose having passed all controls. Those
   omissions do not rescue arms that already fail other mandatory criteria.

4. **Partial collinearity is not a localization experiment.**
   A top-50 leverage weighting correlating with the reliability weighting
   while producing a smaller effect does not alone prove the operative signal
   lies outside the top 50 directions. The weighting and operator response
   differ nonlinearly. The outside-top-50 explanation is a hypothesis consistent
   with these results, not a uniquely identified mechanism.

5. **Unresolved differences stay unresolved.**
   The prune interaction and the Flex contrasts include zero. Marginal
   confidence intervals on separate effects do not establish an ordering
   between them; a paired contrast is required. The reported point ordering
   may be described as such.

## Replacement summary

The tested reversed-reliability tilt shows a small local positive effect, with
the primary dose exceeding its mirror and random controls. No tested dose
meets all precommitted submission criteria. Direct top-50 leverage weighting
captures part of the observed pattern but not its full magnitude. The prune
interaction and transfer to the Flex reference remain unresolved. These
results close the tested submission candidates, not all transformations of
these data or all future modelling approaches.

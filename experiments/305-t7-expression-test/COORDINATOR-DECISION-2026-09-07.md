# Coordinator interpretation rule — before reviewing candidate expression outcomes

The worker identified the empirical-null resolution limit in its preregistration
and disclosed it. That is useful. However, `PREREGISTRATION.md` §7 then changes
the primary "corroborated" conjunction to an **unadjusted p ≤ 0.01**, relegates
Holm to a reported number and adds an across-split null. The takeover mandate
required multiplicity control over all primary candidate pairs and prohibited
substituting another split for an insufficient same-split null.

**Coordinator decision: retain the original evidence standard.** A nominal
outlier may be reported descriptively, but it is not an accepted multiplicity-
controlled primary corroboration. Do not silently weaken the standard because
the original design lacks resolution. This decision is recorded from metadata
and the preregistration, before the coordinator has reviewed any candidate-gene
expression outcome.

## Resolution calculation

There are 11 primary construct–gene pairs. Same-split targeting-construct counts
are 158 (train), 54 (validation) and 105 (test), including the tested construct.
After its exclusion the minimum +1-corrected rank p is respectively 1/158,
1/54 and 1/105. With the full 11-test family, even the smallest attainable
first Holm value is 11/158 ≈ 0.0696, above 0.05. A gene failing a measurement
filter is not evidence against its null; retain all 11 predeclared tests in the
formal family (untestable hypotheses cannot be claimed rejected).

The validation pair cannot even attain the worker's nominal p ≤ 0.01:
its floor is 1/54 ≈ 0.0185. Such a miss is a resolution limitation, not
evidence of absent repression.

## Reporting

- Preserve the worker's preregistration, timings and all numbers unchanged.
  Disclose its pre-expression deviation from the mandate.
- Report direct effects, both cell halves, control effects and paired-construct
  specificity. These observations can be informative without a formal rejection.
- Label the worker's raw-rank conjunction **NOMINAL / DESCRIPTIVE**, not an
  accepted primary discovery. Report the original multiplicity-controlled
  criterion separately, with its attainable-resolution limitation.
- Any across-split null is secondary/descriptive and carries between-split
  comparability assumptions. It cannot replace the frozen same-split primary.
- No statistical outcome here proves off-target binding or excludes off-target
  mechanisms; the sequence scan and expression data are observational.

## Experimental-unit correction

The 11 primary pairs belong to **8** candidate-bearing targeting constructs
(CAST_P1, GRB10_P2, KDM2B_ENST00000377071.4, KLF10_P1, DPF2_P1,
IGF1R_P1, MED15_P1, TMSB4X_P2), not 9 as §1/§8 state.
The original screen had 17 paired intended targets / 34 targeting constructs;
only 8 targets bear the primary candidates. State which population each
summary or resampling calculation actually uses.

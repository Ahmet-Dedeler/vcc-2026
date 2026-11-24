# Coordinator review of 305 — promising observations, no accepted primary discovery

305 completed successfully. The pre-outcome coordinator decision in
`COORDINATOR-DECISION-2026-09-07.md` remains binding; the worker's replacement
of familywise control with nominal p ≤ .01 is not adopted.

Independent reconstruction from integer sufficient sums, without importing
305's analysis helpers: `scripts/audit_round5.py`, output
`work/round5-coordinator-review.json`. All 23 rows' pooled effects, half effects,
fold changes and same-split ranks reproduce to maximum absolute error
2.22e-16. This checks arithmetic, not independent biology.

## Accepted observations

| pair | log1p-CP50k effect | CPM fold change | cell halves | full-11 Holm p |
|---|---:|---:|---|---:|
| GRB10_P2 → SNX8 | -.492396 | .45070 | -.51923 / -.46554 | .10476 |
| DPF2_P1 → CEP120 | -.632763 | .34455 | -.65016 / -.60395 | .10476 |

The partners GRB10_P1 and DPF2_P2 have effects +.0453 and -.0005 on these
candidate genes. This supports construct-specific association. It does not
exclude differences in on-target efficacy, population composition or handling,
and is not proof of seed-mediated binding.

The 11 primary pairs span 8 intended targets/constructs; 3 pairs fail the
measurement floor and stay in the formal family with p=1. **Zero of 11 pass.**
The primary rank test lacks enough resolution to reject at the required FWER;
absence of rejection is not evidence of absence of repression.

NTC sums are identical across splits. The 12 comparator evaluations are
**4 unique construct–gene pairs**, not 12 independent controls. Two of these
show nominal strong repression (DCTN2, PRKCQ). Because their sequence matches
were selected by the same rule, “non-targeting” is a design label: they are
candidate-bearing observations, not necessarily negative controls for this
specific mechanism.

## Post-hoc results: useful diagnostics, not repaired confirmation

The expression-matched base-rate and lane-matched analyses were introduced
after outcomes. Keep them, but do not use their Poisson-binomial
P=1.8e-9 as a confirmatory rescue:

- Several candidate genes share a construct, and all comparisons share control
  data. Independent Bernoulli trials are not established.
- Expression matching alone does not reproduce selection by seed/PAM,
  promoter architecture, GC content, TSS density or mappability.
- Zero observed rates in finite comparison pools are uncertain, not known
  zero probabilities.
- The script's all-gene null does not implement the candidate-specific
  exclusion of constructs targeting each null gene. It is not literally the
  same test for every arbitrary gene.
- The claim that these controls “can only weaken” is false: the report uses
  them to introduce a very small new P value and a positive headline.
- No change under lane reweighting excludes only that measured reweighting
  explanation, not all batch/handling confounding.

The small exceptional sets for SNX8 and the NTC candidates are **worth a
properly designed follow-up**. CEP120 is less compelling because 120 genes
under the same 95-cell DPF2_P1 construct meet the nominal rule.

## Next scientific step, not a production proposal

Before examining new outcomes, specify one aggregate candidate-enrichment
test with construct/target as the unit; preserve the number and structure of
candidate matches in a sequence/promoter- and expression-matched joint null.
Use new guide constructs or an independent screen where possible. A
sequence-mutant or occupancy assay would test the mechanism directly.

Existing halves can audit implementation and estimate cell-sampling error;
they cannot become new biological replicates. The existing data do not
establish that off-targeting explains failing-vs-working differences, and do
not exclude off-targeting generally. No NTC is removed, no model is changed,
and no submission is justified by this experiment.

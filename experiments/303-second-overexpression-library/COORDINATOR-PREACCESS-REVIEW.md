# 303 pre-access coordinator review — September 7, 2026

The preflight worker finished successfully without a reserve outcome.
Its PREACCESS and DRAFT files remain the historical record. This review
specifies corrections before the final preregistration and any release.

## Corrections

1. **58 is not an error.** All 175 scored rows of `symmetry302.csv` have S_g.
   R5 overlaps 58 of them and 45 of the 148 primary rows. Keep both
   descriptive comparisons; do not silently narrow the original cohort.
2. **Distinguish file and content hashes.** The cohort symbol-list hash
   fca12254… is not the hash of the multi-column COHORT.tsv. The final
   release must bind the actual file and its row identities.
3. **Same columns in the primary null.** Remove all cohort target-transcript
   columns (98 present) from both response matrices, for every pairing:
   15,309 measured columns become 15,211 scored columns. The draft's rule
   removes one column for a matched pair and two for most shuffled pairs.
   Retain that draft rule as a labelled sensitivity, not primary inference.
4. **Uniform permutations for the primary reference distribution.** Sample
   from the full target-permutation group, fixed points allowed, rather than
   only derangements. Hold KD rows/features fixed and permute OE matching.
   Any formal interpretation remains conditional on target-label
   exchangeability; it is not a mechanistic causal null.
5. **Direction conditions.** C7a requires observed median S < 0 and negative
   tail excess > 0 as well as the conjunction permutation/Holm rule. The
   three association cells require a positive observed correlation in their
   declared direction. “Less positive than the null” is not antisymmetry.
6. **Full family and missingness.** Four cells always enter Holm; an
   unestimable cell gets p=1 and cannot pass. Report bootstrap intervals as
   conditional descriptive target-resampling intervals, not study replication.
7. **Complete code before access.** The worker produced a statistic kernel
   but no complete end-to-end analysis. Complete and test it before freezing;
   do not invent analyses after observing R5.
8. **Guard must verify current inputs.** A cached source.sha256 or cached
   metadata digest is not live verification. Bind source, current row
   identities, cohort/axis, KD bank, feature inputs, known-answer outputs,
   final preregistration and all analysis/access code. An invocation of the
   public run_access function must not bypass marker/release checks.

## Scientific scope

The comparison is cross-study and treated as cross-line (Nakatake SEES3
versus Pan H1), not a matched-assay repeat. Exact platform provenance of the
producer's input table is unresolved. The processed file has one profile per
target and a shared Emerald control; it supplies no replicate structure.
No bulk-perfect-reliability assumption, disattenuation or latent-equivalence
claim is allowed. A failure of the transform/measurement validity gates
voids inference without being relabelled a biological negative.

The test remains worth doing as a bounded observed-response comparison.
It is not a competition candidate and cannot settle the general question
of gain-versus-loss-of-function symmetry. No family-adjusted generalization
claim will be made from the conditional target bootstrap.

## Ownership and release

The two code review workers have disjoint ownership (access guard/tests;
analysis kernel/driver/tests). Neither may open R5. The coordinator alone
freezes PREREGISTRATION.md, commits it and the tested code, enters C7 in the
ledger, then issues one hash-bound RELEASE.json. Until then R5 stays locked.

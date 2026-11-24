# 305 — Expression evidence for the already-frozen T7 seed candidates

The user asked Codex to continue the interrupted Round 5 on 2026-09-07.
This is the missing third mandate: finish the EXPRESSION check that 301
explicitly could not perform from 293's scalar-only saved output.

## Scope

Own only `experiments/305-t7-expression-test/`. Read READ-FIRST, AGENTS, the
research programme, 301 RESULTS and seed_contrast.py, 298's seed298_plan.md
and seed298.csv, and 293's RESULTS / PREREGISTRATION / work/lib293.py,
run293.py and agg293.py. Use `apply_patch` for authored files and `uv run python`.
Do not edit the ledger, shared docs, old experiment code, data, or submissions.
No model building, vcc, cell-eval2, publication, commits, new agents, downloads,
or new seed search. R1–R5 are off limits. The 2025 train/validation/test files
already used by 293 are permitted read-only; they are not a new reserve.

Time box 2 hours, peak RSS <6 GiB, stream the three files SERIALly. Check the
repo's memory gate before each stream. Avoid loading a whole single-cell matrix
or duplicating the old 400-cell grid. Prefer reusable sufficient sums, and
do not invoke old scripts that write outside your directory.

303 is validating a second overexpression library without reading its reserve;
304 is auditing F15 and novelty. Do not read or write their directories.

## Question and frozen scope

301's distal 12-nt, PAM-anchored, active-window candidate set is fixed. It is
4/13 failing constructs and 4/21 working constructs, NOT evidence of equivalence
or absence of off-targeting (Fisher P=0.679 is inconclusive). Determine whether
those specific candidate genes are unusually repressed in the construct's
existing expression data, and whether that repression repeats across disjoint
cell/control halves. This is observational corroboration, NOT causal proof of
off-target binding and NOT a wet-lab design task.

Before reading any NEW candidate expression, write and hash a preregistration
with the exact candidate-pair table, source file identities, split/construct
mapping, response formula, seeds, nulls and pass rules. The identities and
failing labels are already known; explicitly call this an extension on used
data, not an untouched confirmation. Do not change seed length, window,
distance filter, candidate gene, failing definition or effect threshold after
seeing expression.

## Implementation and required checks

1. Reconstruct the primary candidate table from 301's saved hits and its exact
   filters; reproduce its 4/13, 4/21 and 4/31 background counts, handling duplicate
   spacers/hits as one construct-gene pair. Freeze all evaluated pairs, not just
   failing-arm hits. Retain unresolved/unmeasured candidates as NOT TESTABLE.
2. Inspect 293's saved artifacts first. Reuse sufficient sums if available;
   otherwise stream exact construct and same-split non-targeting control sums
   from the original 2025 files. Keep the full library total for normalization
   even if storing only the columns required for candidates/positive controls.
   Save reproducible per-construct, per-half sufficient statistics and hashes.
3. Freeze a common response formula from lib293/lib285. Report the direct
   control-referenced effect for candidate repression, not only a residual:
   being below the panel mean is not the same as going below controls.
   Report residual empirical ranks as a separate statistic.
4. Known answers BEFORE candidate conclusions: 34 construct identities and the
   13/21 labels; gene order and missingness; a synthetic known fold change;
   chunked sums vs a direct small non-reserve slice; a same-NTC-halves null;
   selected own-target effects through 293's ORIGINAL cell selection/regime
   reproduced at a declared dtype-aware tolerance. Do not require a full-depth
   re-extraction to equal a differently thinned historical number. Document
   exact comparable configurations. Abort conclusions on a failed check.
5. Two disjoint deterministic cell halves per construct, with DISJOINT
   non-targeting controls from the same split/batch where feasible. Report
   cell and UMI counts, batch/lane balance and whether any control is shared.
   No cell bootstrap labelled as biological replication. Do not invent
   biological replicates when the file lacks them.
6. Primary per-pair evidence: negative direct effect on both independent
   halves AND an unusually negative pooled effect relative to that candidate
   gene under OTHER constructs in the same split. Exclude the candidate's
   own-targeting constructs and the tested construct from this empirical null.
   State the empirical rank p-value with +1 correction and Holm correction
   over ALL primary candidate pairs. Predeclare a minimum null size and
   mark insufficient cases NOT TESTABLE; do not substitute another split.
   Report effect magnitude even if a statistical criterion passes.
7. A specificity check: compare the SAME candidate gene under the other
   construct of the SAME intended target when present. This helps distinguish
   intended-target downstream biology from a construct-specific association.
   It still does not prove an off-target mechanism. Non-targeting guide
   candidates and control-vs-control effects are negative comparators.
8. Report failing-versus-working summaries, but cluster uncertainty at the
   intended target (17 paired targets), not 34 independent constructs. The
   experiment is underpowered for a population off-target rate. Absence of
   a supported candidate is "no expression corroboration under this test,"
   never "off-targets do not explain failing constructs."

## Deliverables

`PREREGISTRATION.md` + hash; `work/input_manifest.json`;
`work/candidates305.csv`; saved sufficient sums/metadata;
`work/expression305.csv` with effects, both halves, empirical ranks, Holm p,
paired-construct comparator, cell counts and testability;
`work/results305.json`; tests and scripts; `RESULTS.md`.

Give one line: primary pairs evaluated / untestable, number meeting all
preregistered corroboration conditions, counts by failing/working/NTC group,
and the strongest limitation. State all deviations, whether fixed before or
after outcomes, and the peak RSS/runtime. Then stop.

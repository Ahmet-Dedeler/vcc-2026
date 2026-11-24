# Scoped Claude implementation request — do not run or freeze

Ahmet authorizes Claude Code assistance through his existing Max subscription.
Codex coordinates and verifies. You own ONLY files under
`experiments/342-hamilton-current-family/`. Read other repository files as
needed, but never edit them. No agents, commands, model calls, network,
reserve access, response-data inspection, commits, numerical runs, package,
production-default changes, or submissions. Return an implementation summary.

Implement a bounded **current-family remeasurement** of the still-unsubmitted
Hamilton integer-total rounding rider. This is not a new mathematical method:
196/208 tested it on older compositions, and 230 left it live but uncleared.
Read their actual results/code and 323's evaluator before writing.

Primary sources:

- `experiments/196-emitted-depth/ham196.py`
- `experiments/208-free-riders/emit208.py`, `RESULTS.md`
- `experiments/230-toll-magnitude/RESULTS.md`
- `src/vcc2026/donor_transfer.py`, `src/vcc2026/rc_emit.py`
- `experiments/323-current-family-validation/{run323.py,evaluate323.py,exact_de323.py}`
- `experiments/337-availability-consensus/` for a recent all-six-member
  count runner/report and independent dense-budget audit.

Design constraints, to encode in a prospective PREREGISTRATION:

1. Same cap15 effects, controls, seeds, profile, spread tier, call set,
   whitening, CD4 and PM as 323. H1 union remains OFF. No model fitting.
2. Three arms: exact original floor baseline; Hamilton largest remainder;
   uniform-random up-rounding control. For each block let x=n*L*profile_p
   at the SAME production rounding site, base=floor(x), residual count
   round(sum(x))-sum(base). Hamilton assigns that many increments to the
   largest fractional parts, stable gene-index tie breaking.
   The control assigns the SAME number of increments uniformly without
   replacement among coordinates with fractional part >0, using a separate
   generator seeded 3420909+target_index. Never consume emission RNG for
   this rounding choice. Every coordinate differs from floor by 0 or 1.
3. Modify only this allocation site through an experiment-local, auditable
   adapter. No global NumPy monkeypatch. Do not simply add counts after
   final emission: the new totals must go through the real pass-1/rank/
   placement machinery. Keep every other emission operation unchanged.
   A carefully checked local copy or explicit AST seam is acceptable if
   its exact intended differences are attested and floor is bit-identical.
4. Re-emit all 85 cap15 baseline blocks and require the whole emitted-count
   digest, per-gene totals and all six metrics to match 323 EXACTLY before
   interpreting changed arms. Save sufficient statistics and per-target
   rows. Use exact_de323 scoring, not the emitter's older arithmetic.
5. Report all six raw metrics and H1-local scaled overall; raw MSE and its
   jackknife credit must be recomputed, never assumed pinned. Use the same
   paired 10,000-target bootstrap with a fixed seed, including recomputation
   of the whole-panel MSE budget per resample. No board conversion factor.
6. Prospective local-advancement gates: all integrity checks pass; Hamilton
   overall point >=+.001 versus baseline; paired overall lower bound >0
   versus baseline AND random-up control; no raw MSE or NMAE contrast whose
   entire 95% interval is harmful; and no fidelity/reach/Jaccard contrast
   whose entire 95% interval is harmful. A pass means only
   LOCAL_GAIN_REQUIRES_SEPARATE_ABC_REVIEW, not a cleared package or slot.
7. Preserve any failure and do not tune. Explicitly state rounding changes
   pseudobulk totals and no competition member is assumed invariant.
8. Unit tests must cover known fractional totals, zero/integer entries,
   stable ties, mass conservation, 0/1 differences, random-control mass/
   independent RNG, unchanged floor path and injection-site specificity.

Deliver code, pure tests and draft registration in your owned directory.
Provide `freeze`, `verify`, `run` entry points, but DO NOT invoke them.
Freeze must hash-bind code, registration, effects/reference/evaluator and
production source dependencies before generation. It must refuse a dirty/
uncommitted registered implementation when running. Do not let your worker
read numerical .npz/.npy response artifacts while implementing.

Numerical execution is Codex's responsibility: one process at a time,
16 GiB tree RSS ceiling, repository launch gate, 40-minute run budget.
It must not start after 19:35 UTC without a new coordinator decision.
Handoff is 20:58 UTC; return your implementation promptly.

Do not implement the empirical-null call-set proposal or block-sign repair:
those were NOT selected. If a requested seam cannot preserve the baseline,
state the blocker instead of weakening its exact known-answer requirement.

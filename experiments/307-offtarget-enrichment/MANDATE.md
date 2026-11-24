# 307 — The repression leads, tested properly: ONE aggregate candidate-enrichment statistic with the CONSTRUCT as the experimental unit, against a joint null that preserves how the candidates were selected (seed/PAM, promoter architecture, GC, mappability) AND their expression

You are a research subagent. Work in `experiments/307-offtarget-enrichment/`. **Edit nothing outside it** except
appending your preregistration and its result to `docs/CONFIRMATION-LEDGER.md` under a heading you create,
`## Claim C8 — 307 (no reserve: this uses data already read by 305/301)`. No `vcc`, no `cell-eval2`, no build,
nothing under `submissions/` or `data/`, no download. `uv run python`. Peak RSS under ~6 GB. **Time box: 3.5 hours.**
Siblings (306 and others) run in parallel; do not read or write their directories.

## Why
305 asked whether the seed-match candidates frozen by 301 are repressed in the 2025 expression data. Its own
coordinator review (`experiments/305-t7-expression-test/COORDINATOR-REVIEW-2026-09-07.md`, and the pre-outcome
`COORDINATOR-DECISION-2026-09-07.md`, both binding) concluded: **zero of the 11 preregistered primary pairs passes
Holm**, the same-split rank test cannot attain the required familywise error at this resolution (the smallest
attainable first Holm value is 11/158 ≈ 0.0696), and the post-hoc Poisson-binomial P = 1.8e-9 is **not** an accepted
rescue because (i) candidate genes share constructs and all comparisons share control data, so independent Bernoulli
trials are not established; (ii) expression matching alone does not reproduce selection by seed/PAM, promoter
architecture, GC content, TSS density or mappability; (iii) zero observed rates in finite pools are uncertain, not
known zero; (iv) the all-gene null did not implement the candidate-specific exclusion of constructs targeting each
null gene. Two observations survive as **worth a properly designed follow-up**: `GRB10_P2 → SNX8` (partner
`GRB10_P1` effect +0.045) and the two non-targeting comparators (`DCTN2`, `PRKCQ`); `DPF2_P1 → CEP120` is weaker
because 120 genes under that same 95-cell construct meet the nominal rule.

**Your job is the aggregate test the review specified, designed and frozen BEFORE you look at any outcome.** It is
one statistic with the construct as the unit, not 11 per-pair tests. A well-supported negative is the expected
result and is a valid result.

## READ FIRST
`experiments/305-t7-expression-test/` — `COORDINATOR-DECISION-2026-09-07.md` and `COORDINATOR-REVIEW-2026-09-07.md`
FIRST, then `PREREGISTRATION.md`, `RESULTS.md`, and its `work/` (`lib305.py`, `sums305.py`, `candidates305.py`,
`candidates305.csv`, `expression305.csv`, `results305.json`, the saved integer sum `.npz` files, `input_manifest.json`).
`experiments/301-metadata-acquisitions/RESULTS.md` and its `work/` (`vcc2025_spacers.csv`, the promoter fetches, the
seed-scan code) — the candidate selection rule you must reproduce and preserve. `docs/RESEARCH-PROGRAMME-2026-09.md`
rules 1, 5, 6, 9. `docs/FINDINGS-2026-09.md` F14 and the T7 wording.

**Reuse 305's saved integer sufficient sums wherever they answer the question.** Re-streaming the expression data is
allowed only if a statistic genuinely needs a quantity the sums do not contain; if you re-stream, say exactly why,
and carry a known-answer row proving your stream reproduces 305's sums (it reproduced them to 2.22e-16).

## The design (write `PREREGISTRATION.md`, hash it into `work/PREREGISTRATION.sha256`, append it to the ledger as C8, BEFORE computing any candidate outcome)
- **Unit:** the targeting construct (and, separately declared, the non-targeting construct). State the exact
  population: 305's review says the 11 primary pairs belong to **8** candidate-bearing targeting constructs
  (CAST_P1, GRB10_P2, KDM2B_ENST00000377071.4, KLF10_P1, DPF2_P1, IGF1R_P1, MED15_P1, TMSB4X_P2) out of 34 targeting
  constructs / 17 paired intended targets. Verify that count yourself and report what you find.
- **Statistic (ONE primary):** an aggregate enrichment of repression among each construct's candidate set relative to
  the same construct's matched non-candidate genes, pooled over constructs with the construct as the resampled unit.
  Define it precisely (e.g. a per-construct standardised mean candidate effect, or a per-construct rank statistic,
  pooled by a construct-level bootstrap or sign test). It must be computable for **every** construct, so that a
  construct with no strong candidate contributes evidence rather than being dropped.
- **The joint null — the heart of this mandate.** Resample candidate sets, not genes independently. For each
  construct, draw pseudo-candidate sets of the SAME SIZE whose members match the real candidates jointly on:
  (a) the sequence rule that selected them (seed length/PAM class and where in the promoter window the match falls —
  reproduce 301's rule; if an exact re-run of the scan is too expensive, state the approximation and bound its
  effect), (b) promoter GC and TSS density, (c) mappability or a stated proxy, (d) basal expression decile in the
  control cells. **Implement the candidate-specific exclusion 305's null omitted:** a null gene drawn for construct
  c must exclude any gene that some construct in the library intends to target, and must exclude c's own intended
  target. Report how many draws are feasible and the resulting resolution.
- **Pass rule:** state it before you compute it, with the number of verdict cells, the multiplicity correction over
  them, and the attainable-resolution floor (as 305's review computed: say the smallest p your null can produce).
  **If the design cannot attain significance at your stated α, say so IN THE PREREGISTRATION and report the test as
  a bounded estimate with an interval instead of a verdict.** Do not weaken the standard after seeing outcomes —
  that is the specific failure the coordinator rejected in 305.
- **Declared secondaries (reported, no verdict):** the failing-vs-working construct contrast (301's original T7
  question); the two non-targeting comparators treated as what they are — candidate-bearing observations selected by
  the same rule, not negative controls; the `GRB10_P2 → SNX8` pair with its partner construct; per-construct cell
  counts and the measurement floor.
- **Known-answer rows (before any outcome):** 305's pooled effects reproduce from its saved sums; the candidate
  selection reproduces 301's counts; a no-signal row (a construct against another construct's candidate set) sits at
  the null; the NTC sums are identical across splits (305 found this — confirm and say what it implies).

## Rules of evidence
- Name the experimental unit in every interval and every p (programme rule 1).
- A construct that fails a measurement floor stays in the family; an untestable hypothesis is not a rejected one.
- Nothing here can prove off-target binding: expression is observational and the constructs differ in on-target
  efficacy, population composition and handling. Say what a positive would and would not establish, and say it in
  the abstract of your RESULTS.md, not only in a limitations paragraph.
- If the honest answer is "the existing data cannot resolve this", that is the deliverable — state what design
  would (new constructs, an independent screen, a sequence-mutant or occupancy assay) and what it would cost.

## Deliverable
`PREREGISTRATION.md` (hashed, in the ledger as C8 before outcomes), `RESULTS.md` (the KA rows; the population and
unit; the null's construction and its attainable resolution; the ONE primary statistic with its interval and
verdict-or-bound; the secondaries, labelled; the limits), `work/enrichment307.csv` (per construct: candidate count,
matched-null draws, the statistic, its null distribution summary, cells, floor), and one line:
**C8 (aggregate candidate-repression enrichment, construct as unit): PASS / FAIL / UNRESOLVABLE-AT-THIS-RESOLUTION
(statistic … [lo, hi]; null …; attainable floor p …).** No model, no submission, no production change. Then stop.

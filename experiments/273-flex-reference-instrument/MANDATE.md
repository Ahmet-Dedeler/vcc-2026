# 273 — A RELIABLE-REFERENCE instrument: h1-test (10x Flex, exactly 400 cells per target, split-half reliability 0.48) with the H1 union OFF so it is not leaked — validated against the cloud scorer — then re-read the noise-limited verdicts (272's operator damping, 271's shrinkage, 249's deflation) on it. ≤ 5 hours. Measurement only. No build, no submission.

You are a research subagent. Work in `experiments/273-flex-reference-instrument/`. **Edit nothing outside it.** No build,
no `vcc`, no `cell-eval2` RUN (the cloud scorer's existing outputs are your known answers; do not launch it), nothing
under `submissions/`, `data/harness/` or `data/`. `uv run python`. Peak RSS under ~10 GB, one heavy load at a time
(`data/harness/h1-test.h5ad` is 0.68 GB, `h1-full.h5ad` 1.55 GB; `h1-test-baseline-pred.h5ad` 2.7 GB — backed reads).
**Time box: write up after ~5 hours.** `docs/RESEARCH-PROGRAMME-2026-09.md` binds you.

## Why (R267 §3, R271 §5, R272 §0.10)
Every local verdict of the programme was read on 162's PLOLO, whose truths are 150–200-cell Orion/K562 references with
split-half residual reliability ≈ 0.03–0.05 at half depth (0.05–0.09 at full depth) — a third of Orion targets have no
measurable direction. 272 found that `pds` on those truths is bought by amplifying the directions where donor rows
are mostly sampling noise, and that removing the noise removes the signal. **Two things are unknown and decide how
to read the last week: (i) is that true on a reference reliable enough to resolve a target's direction, and (ii) what
is the achievable per-target cosine / `pds` regime when the reference is Flex at 400 cells — the challenge's own
count and chemistry?** `data/harness/h1-test.h5ad` is exactly that (159 §1.3: resid reliability **0.481 measured at
400 cells**; `scripts/make_harness.py`: 2026-faithful downsampling), scored by the real `cell-eval2` on the cloud in 221
(nine of ten statistics bit-exact). Its trap is donor leakage: `h1_2025` covers 85/85 of its targets (AGENTS.md "TWO
TRAPS"), so **every arm here runs with the H1 union OFF** (`VCC_H1_UNION_G` unset / the pair absent) — the base is the
three production donors + CD4 (deflated k = 5, re-solved) only; that configuration is not leaked. State plainly that
H1 is embryonic and the challenge lines are somatic (159/010): the instrument reads the reliability REGIME, not the
biology of A/B/C.

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; `experiments/159-transfer-ceiling/` (`reliability.py`, `shards.py`: how h1-test is
loaded, split, and its reliability computed — reproduce 0.481 / 0.3300 as KA-1); `experiments/221-cloud-scorer/RESULTS.md`
(which predictions were scored on h1-test by the real scorer, their `pds` / `mse` rows and the prediction files —
KA-2 must reproduce one of those `pds` numbers locally to the digit with the scorer's own per-target rank statistic;
if none of 221's scored predictions is union-free, use the ones that predate the union (`rc_callset`, `cd4f034` era)
and say so); `scripts/make_harness.py`; `experiments/262-where-pds-is-lost/diag262.py` (the per-target `pds` statistic,
262 KA-3) and `experiments/246-direction-hypotheses/lib246.py` (`score_effects`); `experiments/272-noise-corrected-operator/op272.py`
and `RESULTS.md` §0 (constructions A′ c ∈ {0.25, 1, 4} and the cross-fitted operator); `experiments/271-*/est271.py`
(EB row shrinkage); `experiments/257-*/lib257.py` (the deflation on/off); `experiments/162-retrodicting-validator/plolo.py`
(how a prediction for a held context is composed from the donors + CD4 in the residual space — you need the same
composition for the H1 context: its control profile is in `h1-test.h5ad`); AGENTS.md §"TWO TRAPS", §8.1b.

## Stage 1 (this mandate)
1. **The instrument.** Compose the union-OFF prediction for the H1 context (production operator `union_pruned` —
   note it was fitted on donors that include no H1 rows — donors + CD4 deflated k = 5 re-solved, no union), render it
   as the fold instrument renders (`profile_from_effect`, minnorm), and score per-target `pds` and cosine against
   h1-test's truth with the scorer's statistic. **KA-1:** 159's reliability rows reproduce. **KA-2:** for a prediction
   221 scored on the cloud, the local per-target `pds` mean equals the scorer's `pds` (state the max |Δ|); if the
   local statistic cannot be made to agree, stop and report why — the instrument is not validated.
2. **The reference's reliability by unit** (extend 267's table to a Flex reference): split-half residual cosine at
   200 / 200 by batch (h1-test has batches? say what units exist), by guide if present, random cells; the share of
   targets ≤ 0; Spearman–Brown to 400; compare with 267's Orion / K562 numbers side by side.
3. **The conversion curve on a reliable reference** (262 §4.2's knee): per-target `pds` vs observed cosine on h1-test;
   where the knee sits (PLOLO: 0.02 observed ≈ 0.07–0.09 signal cosine, R271 §5); the fraction of targets in the
   below-chance / marginal / saturated regimes.
4. **Re-read three noise-limited verdicts on this reference (union OFF everywhere; each vs its own union-OFF
   incumbent):** (a) 272's A′ operator damping at c ∈ {0.25, 1, 4} and, if 272's cross-fitted `orion_hct116` operators
   are reusable for the H1 composition, the cross-fitted operator; (b) 271's EB row shrinkage; (c) 249's CD4 deflation
   on vs off (the one arm with a board point: it bought +0.0149 raw `pds` on A/B/C — its sign here is the check that
   the instrument sees what the board saw). Report per-target Δ`pds` and Δcosine with paired target bootstraps (85
   targets — say so; the intervals will be wide), and the same by 262-style groups defined ON THIS reference.
5. **The reading:** for each of (a)–(c), does the PLOLO verdict hold, reverse, or blur on a reference with ρ 0.48?
   Specifically: does damping the noise-dominated directions still lose `pds` when the truth can resolve directions?

## Deliverable
`RESULTS.md`: the instrument and its two KA rows (or the honest failure); the reliability table beside 267's; the
conversion curve; the three re-reads with bootstraps; one line each: **(a) HOLDS / REVERSES / BLURS, (b) …, (c) …**,
and a one-paragraph statement of what the Flex-400 regime implies for the board's `pds` ceiling. No price (85 targets
of an embryonic line do not price a board arm), no preregistration, no model name, no submission. Then stop.

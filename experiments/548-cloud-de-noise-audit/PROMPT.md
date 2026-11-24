You are an adversarial auditor for the Arc Virtual Cell Challenge 2026 repo (this checkout, branch master). Read AGENTS.md
first (§4 scoring, §11 method rules, especially rule 13 and §11.1c.4), then docs/FINDINGS.md §0–§2 and §8.1b/§8.1c (search
for "0.004502" and "digest-identical"), docs/FINAL-PHASE-FEEDBACK-PLAN.md §1.4, experiments/547-de-noise-check/ (RESULTS.md,
pairs547.py, pairs547.json) and experiments/537-arm-pipeline/RESULTS.md (§0 A and §2, a worker's draft; it may change).
This VM has NO gitignored data; data/board-entries/*.json and submissions/*.meta.json ARE tracked and are your main inputs.

Experiment 548. Write ONLY under experiments/548-cloud-de-noise-audit/. Commit to your branch only. Never run `vcc`, never
contact anyone, never touch GitHub issues/PRs, never edit files outside your directory.

## The question
AGENTS rule 13 says the DE ("ordering") channel carries sd ≈ 0.0045 overall between digest-identical arms, and the
final-phase plan says a DE arm needs > ~0.009 overall (or two seeds) before adoption. A worker (537) and the coordinator
(547) now claim that this number is a spread of different EMISSION DESIGNS measured in experiment 230, not the re-roll noise
of one design, and that arms which keep the DE design move the DE members by only ~0.0002–0.0004 overall on the board.
If that is right, every DE bin and the final-phase adoption rule are ~10–20× too strict. If it is wrong, relaxing them would
let noise pass as signal. Find out which, adversarially: try to break the claim.

## Tasks
1. **Reproduce 547** from the tracked files (its pairs, rms and quadrature total). Challenge its pair selection: is every pair
   truly DE-identical (compare every DE-relevant config key in the metas, not only the ones 547 names; the builder writes the
   emitter sha, the DE sign prior, MDE, zero repair, DE magnitude, call-set keys, seeds)? Is `sameline` handled correctly?
   Are there more pairs 547 missed, or pairs it should have excluded? Report the rms with and without each questionable pair.
2. **Trace the 0.0045.** Find experiment 230's record (grep the repo: `experiments/230-*`, FINDINGS §8.1b "digest-identical",
   the History sections). Recompute it from whatever 230 recorded if the inputs are tracked; otherwise quote exactly what it
   says it measured (which arms, which emission laws, which member) and state whether those arms shared an emission design.
3. **The re-roll noise of one design.** Where can a same-design, different-seed difference be measured? (a) board pairs
   whose metas differ only in an emission seed (search every meta for seed-like keys); (b) H1 seed triplets in tracked
   experiment outputs (528, 517, 435, 480, 490: look for per-seed six-member tables in RESULTS.md or tracked json/csv) and
   how they scale from H1's target count to the board's 900 target-rows (300 targets × 3 contexts). Give the estimate with
   its assumptions; say plainly if the tracked data cannot answer it.
4. **Verdict.** For three cases give the sd to use (overall units) with the evidence: (i) bulk-only arms' DE members;
   (ii) DE-knob arms that share every emission draw with their base (prior-only arms, 468); (iii) DE-knob arms that re-draw
   the emission. Then propose exact replacement wording for AGENTS rule 13 (one sentence) and FINAL-PHASE-FEEDBACK-PLAN
   §1.4, or say why the current wording should stay.

## Deliver
`experiments/548-cloud-de-noise-audit/` with your scripts (runnable with `uv run python …` from the repo root on tracked data
only), their outputs (json/csv), and `RESULTS.md`: the reproduction, the challenges and their effect, the 230 trace, the
re-roll estimate, the verdict table, the proposed wording, and what you could not determine. Commit to branch
`cloud/548-de-noise-audit`. Be exact; cite file:line for every number you take from the record.

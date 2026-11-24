# 605 — adversarial cloud review of the council's three board-bound patches (Cursor cloud agent, Grok 4.7 xhigh)

You are reviewing code in the GitHub repo Ahmet-Dedeler/virtual-cell-challenge-2026 (branch master). Write ONLY to a new
branch `cloud/605-council-review`, in `experiments/605-cloud-council-review/`. Do not edit any other path. Do not run `vcc`.
The VM has NO gitignored data (no challenge h5ad, no donors, no harness), so this is a code-and-logic review. Where a
claim needs data, say what you would check and why.

## Context (read first)
- `AGENTS.md`: §4 scoring, §6.2 slot rules, §11.1/11.1c method rules.
- `docs/reference/vcc2026-metrics-BRIEF.md`: the scorer's six members, exact definitions.
- `experiments/600-breakthrough-council/SYNTHESIS.md`: what the council concluded today and why.

## What to review (each is an UNAPPLIED patch the coordinator plans to apply on 2026-10-06 with known-answer rebuilds)
1. **630 rank-neutral UP pass** (`experiments/630-rank-neutral-up/`: `src_rank_up.py`, `patch_builder630.py`,
   `test630.py`, `test_rank_up630.py`, `PREREGISTRATION-DRAFT.md`, `RESULTS.md`). Claim: after zs, it raises one carrier
   per uncalled UP gene just above the control max, paid by ballast donations, with per-gene totals and per-cell
   libraries exact, |ΔAUC| ≤ 1e-4 per gene, and 0 call flips (a guard), so pds / mse / fid / jac are exactly unchanged
   and only nmae moves.
   Hunt for:
   - any path where a per-cell library or a per-gene total changes;
   - any path where the scorer's BH call set can change (the guard's coverage);
   - AUC accounting with ties (the scorer uses Wilcoxon with the normal approximation on per-cell CPM);
   - the interaction with zs and the ballast;
   - knob parsing;
   - the unset path being byte-identical.
2. **612 key-gated hook** (`experiments/612-key-gated-hook/`, `patch612.py`, `test612.py`). Claim: a source-agnostic
   `VCC_DE_MDE_SIGN_KEYGATE=<q>` that, per target, keeps the key's sign where |key| is in the target's top half and
   the hook source's sign elsewhere, and is a no-op when unset. Check the "top half of nonzero |key|" definition, ties,
   targets with no nonzero key, and the composition order 598 → 612 → 630.
3. **604 identity-prior plumbing** (`experiments/604-identity-prior/`: `def_identity_prior.py`, `patch604.py`,
   `test604.py`, `build604.py`, `make_refs604.py`, `RESULTS.md`). Claim: `identity_prior()` computes
   `score = z(P_C4 − 0.5) + λ·z(I)` with Flex-only references (one per identity call, never the context's own), and
   quantile-maps it onto CLEAN4's P. λ 0 reproduces CLEAN4 exactly (asserted). `patch604` adds
   `--programme-source identity:<lam>` to the POST-BUNDLE `scripts/def_sameline.py` and leaves every other path byte for
   byte.
   Hunt for:
   - quantile-map edge cases (ties, NaN, identity-only interpolation outside the covered score range);
   - the guard's de-duplication logic;
   - z-scores over different gene sets per context;
   - the psign convention (−1 = down, P > 0.5);
   - whether the day-of path on a relabelled bundle really equals the board form;
   - anything in def_sameline's matched-context handling (rows / signs / w) that the identity branch could break.

## Deliverable
`experiments/605-cloud-council-review/REVIEW.md` on your branch. For each finding give:
- severity (BLOCKER / MAJOR / MINOR / NOTE);
- file:line;
- the concrete failure scenario (inputs → wrong output);
- a suggested fix (code if short).

Also list what you verified as correct, with the reasoning. Be adversarial: a NULL review is fine if nothing is wrong,
but every claimed property above should be either confirmed with a line-level argument or refuted.
Commit and push your branch. Do not open a PR.

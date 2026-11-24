# 607 — second adversarial cloud review: the 10-06 apply-day chain at a PINNED commit (Cursor cloud agent, Grok 4.7 xhigh)

You review code in the GitHub repo Ahmet-Dedeler/virtual-cell-challenge-2026 at the commit named in the launch request
(`source.ref`). **Before anything else, run `git rev-parse HEAD` and write it at the top of REVIEW.md.** If it is not the
pinned commit, stop and say so.
- Write ONLY to the branch `cloud/607-review2`, in `experiments/607-cloud-review2/`. Do not edit any other path. Do not
  run `vcc`.
- The VM has NO gitignored data: no challenge h5ad, no donors, no harness. Small tracked npz files ARE in git (e.g.
  `experiments/604-identity-prior/work/*.npz`, `experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz`).
  Use them.
- Where a claim needs data you don't have, say exactly what you would check and why.

A first review (605, `experiments/605-cloud-council-review/REVIEW.md` + `RESULTS.md`) read a STALE commit. Read both:
its findings and the triage say what was already found. **Do not re-report a 605 finding unless it is still present
at the pinned commit; then say "still present" and give the line.**

## Context (read first)
- `AGENTS.md`: §4 scoring, §6.2 slot rules, §11.1 method rules.
- `docs/reference/vcc2026-metrics-BRIEF.md`: the scorer's six members.
- `experiments/606-applyday-composition/compose606.py` + `RESULTS.md`: the apply order and the composed shas.

## What to review (UNAPPLIED patches; the coordinator applies them on 2026-10-06 in this order)
641 → 598 → 612 → 612dd → 630 → 642 → 604 → 644 (606's chain; `test598.canon_tree` builds the post-bundle tree).

1. **630 rank-neutral UP pass**, after its 605 fix (`experiments/630-rank-neutral-up/src_rank_up.py`,
   `patch_builder630.py`, `test630.py`, `test_rank_up630.py`, `RESULTS.md`).
   - Verify the donor-gain fix: no gene can be accepted with a lower mean or a flipped lfc sign. Construct
     adversarial blocks (carrier library above every ballast library; mixed donors; ties at the control max) and RUN
     them through `up_pass2` / `apply` (numpy only, no data needed).
   - Per-cell libraries and per-gene totals stay exact.
   - The guard: what happens when calls still flip after `GUARD_ROUNDS`?
   - The knob grammar (`VCC_DE_UPRN`), and that the unset builder path is byte-identical.
   - The fix is commit 69972aa1. Its owner says production exposure is NONE structurally: ballast rows hold ~995k counts and carriers ≤ 25k on 85/85 H1 blocks. Target these adversarial blocks: carriers whose library is ≥ some ballast library (mixed eligibility); the `cross=True` budget loop with ineligible donors interleaved; quantum / rounds with a recipient reused across rounds; `own_col`; blocks with K = 0 ballast rows. Also check the fail-SAFE guard: after 5 rounds the target keeps its zs block with `ur_guard_failed = 1`.
2. **612 / 612dd key gate** (`experiments/612-key-gated-hook/patch612.py`, `patch612dd.py`, `test612.py`). The 605
   tie-at-cut finding: still present? Check `VCC_DE_MDE_SIGN_KEYGATE` parsing, the no-op when unset, and the def_day
   pass-through.
3. **604 identity prior v3** (`experiments/604-identity-prior/def_identity_prior.py`, `patch604.py`, `test604.py`,
   `build604.py`, `RESULTS.md` top section). Claims:
   - side = sign(zc + λ·zi) with zc = (P − ½)/sd centred at ½, on CLEAN4-covered genes only;
   - confidence = CLEAN4's own |P − ½| distribution assigned by the rank of |score| (stable, ties by own |P − ½|);
   - identity-only genes keep CLEAN4's psign with P NaN;
   - λ 0 == CLEAN4 bit for bit.
   RECOMPUTE the v3 artifacts' properties from the tracked npz files (`ident604v3_flex_lam{0,0.5,1}_ABC.npz`,
   `flex_refs604.npz`, CLEAN4's npz) and confirm or refute each claim. Also check the patch's anchor, the
   matched / unmatched context handling in `def_sameline.py`, and that `identity@<lam>` parses.
4. **644 same-line skip** (`experiments/644-runbook-amendments/patch_skip644.py`, `test_skip644.py`). A skipped context
   must take the unmatched path for every source, including 604's identity source.
5. **The composed chain.** Re-run `uv run python experiments/606-applyday-composition/compose606.py` if `test598.canon_tree`
   works without data. If it does, report the final shas and compare them with 606's RESULTS (latest table). If it does
   not, say which input it needs.

## Deliverable
`experiments/607-cloud-review2/REVIEW.md` on your branch. For each finding give:
- severity (BLOCKER / MAJOR / MINOR / NOTE);
- file:line;
- the concrete failure scenario (inputs → wrong output; prefer an executed counterexample);
- a short fix.

Also list every claimed property you CONFIRMED, with the line-level argument or the executed check. Be adversarial. A
NULL review is fine. Commit and push your branch. Do not open a PR.

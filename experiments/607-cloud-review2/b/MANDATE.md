# 607b — delta review: what changed in the 10-06 chain AFTER 607's pinned commit (Cursor cloud agent, Grok 4.7 xhigh)

You review code in the GitHub repo Ahmet-Dedeler/virtual-cell-challenge-2026 at the commit named in the launch request.
**First run `git rev-parse HEAD` and write it at the top of REVIEW.md. If it is not the pinned commit, stop.**
- Write ONLY to the branch `cloud/607b-delta`, in `experiments/607-cloud-review2/b/`. Do not edit any other path. Do not
  run `vcc`.
- The VM has no gitignored data. Small tracked npz files are in git; use them.
- A parallel review (607, branch `cloud/607-review2`) covers the chain at commit 7dbb1007. Do not repeat it. Review only
  the DELTA below.

## The delta (all UNAPPLIED; the coordinator applies the chain on 2026-10-06; freeze 10-05 18:00Z)
1. **630's module, d8c39c3f → 7e6b3e52** (`experiments/630-rank-neutral-up/src_rank_up.py`; commit c30d64ce).
   `apply()` wraps the pass in try/except AssertionError and returns the INPUT block with `ur_invariant_failed = 1`.
   - Is it bit-identical to d8c39c3f when nothing fires? Diff the two versions (`git show 69972aa1:…` vs HEAD).
   - Do all three paths (normal, guard-failed, invariant-failed) return the same 19 `DIAG_KEYS`?
   - Can the except branch swallow a non-invariant error that should abort (a real bug masked as "fail-safe")?
   - Does a returned input block still satisfy what the builder expects next (dtype, shape, libraries)?
   Run `experiments/630-rank-neutral-up/test_neggain630.py` and `test_rank_up630.py` if they need no data.
2. **630dd** (`experiments/630-rank-neutral-up/patch630dd.py`, `test630dd.py`): the def_day.sh pass-through
   `DEF_DAY_DE_UPRN` → `VCC_DE_UPRN`. Six anchors (DOC / V / E / C / F / A) on the ORIGINAL DE_MAG lines.
   - C checks the meta (variant + module 7e6b3e52) and prints a NON-FATAL per-context report from the final rc_audit
     dump, ending with `630dd: uprn FAILSAFE_TOTAL <n>`.
   - Claims to check: it commutes with 612dd byte for byte; unset = no behaviour change; the validation regex equals
     the builder's; the report can never abort def_day (`set -uo pipefail`, no `-e`); the meta check reads the right
     keys (`config.de_uprn.variant`, `.module_sha256`).
   Run `uv run python experiments/630-rank-neutral-up/test630dd.py` (it builds scratch trees with `test598.canon_tree`;
   if that needs data you don't have, say what).
3. **604's octile flip** (`experiments/604-identity-prior/def_identity_prior.py` `flip_vector` / `identity_prior(...,
   flip_octiles)`, `build604.py --flip`, `test604flip.py`; commit 0b87b2ba).
   - Claim: P_C4 → 1 − P_C4 on the covered genes of the masked octiles (625's `octile` array, sha-pinned, matched by
     name), before scoring; the pinned psign is negated where P ≠ ½; None / 00000000 are bit-identical; λ 0 + a mask
     returns exactly the flipped prior.
   - Recompute from the tracked npz files: the octile partition of CLEAN4's covered genes; F1–F5 of test604flip.
   - Check that the provenance can't make two different masks produce the same artifact sha, and that a mask cannot
     touch an uncovered gene.

## Deliverable
`experiments/607-cloud-review2/b/REVIEW.md` on your branch:
- each finding with severity (BLOCKER / MAJOR / MINOR / NOTE), file:line, the concrete failure scenario (prefer an
  executed counterexample) and a short fix;
- every claim you CONFIRMED, with the executed check.

A NULL review is fine. Commit and push your branch. Do not open a PR.

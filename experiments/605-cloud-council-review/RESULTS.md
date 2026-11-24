# 605 — adversarial cloud review of 630 / 612 / 604: triage against the CURRENT code ([4d0207], 2026-10-02 13:20Z)

The Cursor agent (`bc-f6800f81…`, grok-4.7-xhigh) finished at ~13:05Z. Its review (`REVIEW.md`, branch
`cloud/605-council-review` 49ac0cc0) read commit **2788ccd9 (11:54Z)**. That is stale: 604's directory, 630's builder
patch and module (`patch_builder630.py`, `src_rank_up.py`), and 612 v2 (`VCC_DE_MDE_SIGN_KEYGATE`) all landed after it.
So its "not in this commit" verdicts are artifacts of the snapshot. Its code findings were re-checked here against
master as of 13:15Z.

| finding | applies to current code? | owner | action |
|---|---|---|---|
| **630 BLOCKER: donor gain ≤ 0 clamped to 1e-300.** If the carrier's library ≥ every ballast library, `cw[a] − min(cw[bal]) ≤ 0`. The ask becomes unbounded, every in-budget ballast count moves, and the mean (and the lfc sign) can FALL while totals, libraries and the U cap all pass | **YES**: `src_rank_up.py:190` has the same clamp. `ZS.fill` also takes from donors whose library is smaller than the carrier's, which is a negative gain per count. The revert check (`du` bounds + "carrier above cmax") never checks the mean or the sign. v3's `quantum` caps the ask per round, so the damage is bounded, not zero | d4cd5f | Fix before the 10-06 apply: donors only where `cw[bal] < cw[a]` (i.e. library LARGER than the carrier's; the review's prose says "smaller library", but its code `cw[bal] < cw[a]` is right). Revert the gene if its mean falls or its lfc sign changes. Then measure how often the condition fired on the real A/B/C blocks and on H1 |
| 630 MAJOR: the call guard is fail-open after 5 rounds | **Fixed in the builder module**: `apply()` raises `AssertionError` after `GUARD_ROUNDS` = 5. Note that a raise aborts the whole build. On D/E/F day, reverting that target to the pre-pass block (X0) and recording it is safer than an abort | d4cd5f | consider revert-instead-of-raise |
| 630 MAJOR: mse is not structurally invariant (the jackknife is over cells) | YES as a statement. It was measured equal on H1 (`jk_delta` ≈ +3.3, mse repr-identical), but nothing pins it | d4cd5f | prereg wording: "mse measured unchanged on H1", not "exactly 0 by construction"; check per-target Ĉ vs C on the A/B/C capacity blocks |
| 630 NOTE: reach is not pinned by 0 call flips; \|ΔAUC\| does not pin the scorer's z (tie term) | YES (measurement, not structure) | d4cd5f | the board read covers it |
| 612 MAJOR: `k < np.quantile(k[nz], q)` keeps every gene tied AT the cut on the strong side, so a point mass at the cut gives 0 weak genes | YES: `patch612.py:76–77` (v2 uses the same helper) | 21dd4c | Low risk with continuous float keys. Real if keys tie (capped or rounded design rows). A stable rank split makes q exact. The knob is unset by default on 10-06 (KGID needs a WIN), so this is not urgent |
| 612 MINOR: malformed `clean4kg:` parse error; test612 T4b aborts without the harness npz | v1 knob form only / test hygiene | 21dd4c | optional |
| 612 NOTE: the builder's weak half is pinned to CLEAN4, so 604 cannot ride the hook | **Fixed by 612 v2**: KEYGATE is source-agnostic (`file:` reads 604's pdown) | — | none |
| 612 × 630 do not commute: 612's new UP genes become 630 targets and can be sign-flipped by the 630 BLOCKER | YES, through the BLOCKER | d4cd5f | fixed by the 630 fix |
| 604 (reviewed absent). Its "safe rule" for a port: Flex-only, de-dup by identity call, never the own call, no identity-only genes, a quantile map, NaN → CLEAN4, λ 0 exact | **v3 satisfies it** (`604/RESULTS.md` v3; KAs 39/39). The builder reads `programme_<c>` directly (`make_submission_pm_cd4.py:5687–5689`), so def_sameline's bank `programme()` rule (NaN → −1) is never on 604's path. One quantified deviation: 185 covered genes have no I (H1's pinned profile has only its > 5 CPM genes). Their SIDE equals CLEAN4's, and their confidence moves by mean \|ΔP\| 0.0004–0.003. That matters only for a hook τ > 0, and τ is 0 | 4d0207 | none |

**Bottom line for the coordinator.** The review's "do not apply" verdict came mostly from the stale snapshot. The ONE
finding that blocks 630 as written is the donor-gain clamp. It needs a donor filter plus a mean/sign revert in
`src_rank_up.py`, then a re-run of test630 and the KAs before 10-06. 612 is safe to apply with the knob unset. 604 is
unaffected.

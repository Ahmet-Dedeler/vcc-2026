# 607 — second adversarial cloud review (Grok 4.7 xhigh, pinned 7dbb1007): triage ([4d0207], 2026-10-02 ~13:55Z)

The agent (`bc-6e3c6282…`, branch `cloud/607-review2` 7e5cfea6) checked out exactly the pinned commit, built the full
10-06 chain, and reproduced all twelve 606 pins (13:35Z table). Verdict: **the 10-06 chain holds. No BLOCKER, no MAJOR.**
Three MINORs and several NOTEs. Everything claimed was CONFIRMED by executed checks:
- 630 (fixed module): 0 mean drops and 0 sign flips over 160 random mixed-library blocks; the fail-safe guard; own_col;
  the knob grammar.
- 612: tie fill on 16,160 draws; the no-op when unset.
- 604 v3: artifacts recomputed from the claim by an independent assignment.
- 644: a skipped context takes the unmatched path, including 604's identity source.

| finding | sev | owner | status |
|---|---|---|---|
| 604: a pinned reference sharing a bundle context's NAME overwrote that context in the pool dict (a rehearsal reusing ids A/B/C with other calls) | MINOR | 4d0207 | **FIXED** (pinned refs are renamed `pin:<name>` on collision). test604flip F6 shows the pre-fix module dropping "B" at similarity 1.0 on the same input, and the fixed module keeping `pin:B` (−0.39). The committed A/B/C artifacts are byte-identical (F1) |
| 644: `DEF_DAY_SAMELINE_SKIP=E, F` (with a space) passes bash, then Python refuses " F"; it fails closed | MINOR | [537757] | relay: `strip()` in both |
| 630: `VCC_DE_UPRN=b_k100_q20r0p10` (rounds 0) is accepted and moves nothing | MINOR | d4cd5f | relay: reject rounds < 1 (and quantum < 1) in parse_spec and the builder regex; 630dd's validation must follow the same grammar |
| 604: `identity:` / `identity@` with an empty dose raised ValueError (traceback); a negative dose was accepted | NOTE | 4d0207 | **FIXED** in patch604: SystemExit with the source string, and a dose < 0 refused. test604 T9 added |
| 604: a score of exactly 0 copied the donor's side | NOTE | 4d0207 | **FIXED**: P = ½ (psign keeps the gene's own CLEAN4 sign). It never occurs on the artifacts |
| 612: the `run612.py` pricing harness still keeps a tie at the cut on the strong side (the patch is fixed) | NOTE | 21dd4c | relay: prices and builds diverge only on tied targets; count them |
| 630: the `rnup630.py` harness still fails open; its v1 still clamps the gain (non-shipping) | NOTE | d4cd5f | relay |
| 630: only the carrier nearest the control max is tried; if it is deeper than every ballast row, the gene is skipped (a missed raise, never a wrong one) | NOTE | d4cd5f | relay (optional: try the next candidate) |
| 630: "production exposure none structurally" overstates it. The filter is right for any block; it is idle only when every ballast library exceeds every carrier library, which is an empirical H1 claim | NOTE | d4cd5f | relay (wording) |
| 644: the skip log says "context unchanged" while a programme source then writes that context's programme | NOTE | [537757] | relay: log "unmatched (644)" |
| mse is not pinned by fixed totals (jackknife over cells) | NOTE | d4cd5f | already recorded |

**Not checkable without data, now on the coordinator's 10-06 KA list:**
- an all-unset rebuild against the pre-chain builder: same recipe, same prediction hash (the unset path is
  byte-identical);
- H1's per-target `min(ballast library) > max(carrier library)` on the seed-0 dumps and the A/B/C capacity blocks;
- the count of priced 612 targets that tie at the cut.

**Pins after the 604 fixes (`606/work/compose606_607fix.log`, ALL PASS; with 604f appended also ALL PASS):**
- moved: def_sameline `f147f8bbe114` (was d964aa03), def_identity_prior `7aeb66b44f28` (was 60e98488);
- unchanged: def_day `c2704f3d`, builder `2c6cfe7d`, rank_up `7e6b3e52`, def_identity 7073406f, helpers fb715e38,
  def_pipeline 76314af6, october_recipe 9dbcd914, def_preflight 17f0836f, def_postbuild 557c2b52, gate524 b16e60c0.

607b (the delta review: 630's try/except module 7e6b3e52, 630dd, 604's flip; pinned 7ed5e52c) is running.

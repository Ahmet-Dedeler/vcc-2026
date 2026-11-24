# 819 judge run: 813-public-data (TeloHAEC genome-wide fold line), RESULTS

- **Judged:** `PROPOSAL-frozen.md`, sha256 `46a7a212dfb3b89844de1a56702fdf2caa92fff5cc487aeb364f52e92d6907ab` (re-hashed).
  There was no earlier VERDICT-*.md, so this is a first submission.
- **Verdict:** **REVISE**. The full text is in `VERDICT.md`, sha256 `75b2c11c4c96ec8be193e571bdd76b4ce1aae71ae059b61c81452aa047aaf4c0`.
- **Why not KILL:** use (a), TeloHAEC as a held-out truth line, is not a closed construction. The record's TeloHAEC
  closures (105 direct source, 225 stack donor, 135 atlas) all used the 2,167-target pilot as a source or prediction target.
- **Why not GO:** four substantive gaps.
  1. The kill gate KA-SG uses an OLS slope that errors-in-variables attenuation sends to ≈ 0 for an unbiased estimator
     (≤ ~10 single-guide cells per target).
  2. The P2 falsifier can only fail on 702, and 702's rows need a PIE TeloHAEC inference that the plan does not contain.
  3. `--fold-ops prod` has no operator on any L panel: the L dumps were built with `--op none` and held fits only.
  4. There is no decision role for a non-cancer fold, and no leakage rule for TeloHAEC rows that are already in 811's bank list.

  Smaller items: 105 is uncited, P3's levels are confounded by donor count, and the data-terms check is not yet a hard gate.
- **Compute:** lookups only, plus one sub-second JSON read. Cloud: none launched by this run. At finish, `labcloud.sh status`
  showed two instances, both from other lanes (`lab801-abcp-e`, `de803-b9`); this run did not touch them. Lab total logged
  $37.88.
- **Known answers:** n/a for a judge run (no screen). Each cited number was checked against its source file:
  - 271: +0.047 / +0.028 / +0.028;
  - 615: 0.060 / 0.049, bar 0.10;
  - 811: Orion pair +0.011…+0.015;
  - 225: telohaec@0.28 −0.00120;
  - coverage JSON: 300 / 298–300 in library, ≥ 100 on 298 / 295–299.
- **Left:** 813 resubmits once with R1–R7. Conditions C1–C4 carry into a GO.

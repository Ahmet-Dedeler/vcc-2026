# VERDICT: 813 PROPOSAL (TeloHAEC genome-wide CRISPRi as a fourth fold line + an 811 bank)

Judge run, 2026-10-16 ~09:10Z (clock: `date -u` 2026-10-16T09:10:39Z at writing). Judged the FROZEN copy
`experiments/819-judge/813-public-data-20261016T090443Z/PROPOSAL-frozen.md`, sha256
`46a7a212dfb3b89844de1a56702fdf2caa92fff5cc487aeb364f52e92d6907ab` (re-hashed here; it matches the mandate's prefix). There is
no earlier `VERDICT-*.md` in `experiments/813-public-data/`, so this is a first submission. Nothing was computed beyond
lookups and one sub-second JSON read (`work/meta/telohaec_rep1_coverage_v2.json`).

## Verdict: **REVISE**

The idea is not closed, and it is cheap: a $3 kill gate for a new instrument. It is the right next data move for the lab.
Its stage-0 gates are mostly sound. But four things are wrong or missing:
1. one kill-gate statistic (KA-SG) fails by construction on a correct estimator;
2. the falsifier of the instrument (P2) rests on a known answer (702) whose rows cannot be built on this line as planned,
   and the other two rows are near-certain positives;
3. stage 1's "production operator" does not exist for any L panel;
4. the proposal does not say how the folds enter a decision, or how TeloHAEC rows are kept out of the candidates priced
   on them.

The resubmission can be short. Items R1–R7 below are the whole list.

## The JUDGE.md checks

### 1. Is the mechanism closed already? **No, for use (a). For use (b), one closure is missing from §7.**
- **(a) TeloHAEC as a held-out truth line is new.** No RESULTS.md on disk uses TeloHAEC as a fold truth for bulk arms. Every
  earlier use of the line was as a source or prediction target, on the 2,167-target pilot:
  - 105: a direct source;
  - 223/225: a stack donor, `telohaec@0.28` −0.00120 (`experiments/225-source-sweep/RESULTS.md` l.22, l.159);
  - 122/130/135: the atlas confirmation, terminal NO-GO (`experiments/135-telohaec-atlas-confirmation/RESULTS.md` l.3, l.40).

  The FINDINGS §2 closure "Direction via a fifth donor line: sources saturate (218/225/229/248)"
  (`docs/FINDINGS.md` l.885) is about adding a donor to the production stack. Use (a) adds no donor, so the closure does not
  apply.
- **(b) is a vote in 811's consensus.** 811 runs that construction itself; it is a live lead, not a closure (`docs/LAB-LOG.md`
  l.268–274). Adding the genome-wide rows is pre-counted as a separate 811 try (`experiments/811-lolo-model/PROPOSAL-2.md`
  l.24–26, l.157).
- **Missing from §7: 105.** `experiments/105-telohaec-residual/RESULTS.md` l.14–17 and l.31–33:
  - The TeloHAEC rows as a direct source for HCT116/HEK293T read pds 0.558 / 0.536 against K562's 0.699 / 0.662
    (−0.141 / −0.127). The equal blend loses −0.018 on both.
  - It closes "direct TeloHAEC transfer… reconsidered only with a new mechanism that explains why a model can recover
    useful information despite the large, consistent direct-transfer deficit."

  This is not a closure of (a). But it is the record's direct evidence that TeloHAEC sits far from every cancer donor. It bears
  on P3 (levels) and on whether a fold this far away ranks arms the way the board does. For (b) it is the relevant closure.
  811's pilot ablation already reads TeloHAEC slightly negative (≤ −0.0007, `experiments/811-lolo-model/RESULTS.md` l.294).
  → R6.

### 2. Does the prediction have a falsifier and bins written before any number? **Partly. One gate is broken; the trust falsifier is weak.**
- **The stage-0 gates are concrete and preregistered.** G0-KD, G0-POWER and G0-BLOB have numeric bars (frozen §4, l.128–136).
  Their reference numbers check out against the record:
  - 271's +0.047 / +0.028 / +0.028, and 6 / 32 / 34 % at or below zero (`experiments/271*/RESULTS.md` l.1, l.9);
  - 615's bar 0.10 and its measured 0.060 / 0.049 (`experiments/615-orion-hek-harness/PREREGISTRATION-blob.md` l.48;
    `RESULTS.md` l.33–34).

  KA-SIM's three pass conditions are explicit (l.70–74), and the estimator choice is by simulation only (l.74). Good.
- **KA-SG fails by construction (frozen l.76–82).** The gate requires the OLS slope of A's rows on the single-guide rows to
  lie in [0.8, 1.25], over the 100 stage-0 targets. Two numbers settle it, both from the proposal itself:
  - **Single-guide rows are almost pure noise.** Single-guide cells are ≤ ~15 % of 1,335,923 cells (l.77–78;
    `RESULTS.md` l.61). That is at most ~2.0e5 cells over 19,276 targets, i.e. ≤ ~10 cells per target.
  - **A's own rows are ~90 % noise too.** Its predicted split-half residual cosine is +0.03…+0.08 (l.149).

  An OLS slope of one noisy row on another is attenuated toward 0 by the regressor's reliability (errors-in-variables). Here
  the reliability is far below 0.1, so the expected slope for an UNBIASED estimator is ≈ 0, not 1. As written, the gate would
  STOP a correct estimator and record the false mechanism "the multi-guide rows are biased" (l.82).
  - The second KA-SG condition (the cross-target mean of A − single ≈ 0) is consistent. But it is blind to μ̄, because both
    estimators carry −μ̄. It tests only k-dependent technical effects, and should say so.
  - The proposal also does not say whether the single-guide cells are removed from A's carriers before the comparison.
    Shared cells would put shared noise into both sides.

  → R1.
- **The P2 falsifier (frozen l.152–160) discriminates almost nothing.**
  - **495 cap 6 and 498p2g are near-certain positives on any fold.** They read SHIP on 18/18 BULK-W and 3/3 BULK-ABC folds
    (`docs/LAB-LOG.md` l.127–129).
  - **479sp2 is excluded**, because no lab surface prices its class (LAB-LOG l.206–212).
  - **That leaves 702 as the only row that can fail, and it cannot be built as planned.**
    - On the folds, 702's arm blends the HELD line's own PIE rows (`experiments/804-known-answers/PROTOCOL.md` l.14;
      `arm702.py` l.6–11).
    - For K562 those came from a dedicated PIE inference job with the `replogle_xdataset` checkpoint
      (`experiments/804-known-answers/pie804/pie804.sh` l.2–3, l.29–31).
    - No PIE TeloHAEC rows exist. Stages 0–2 contain no PIE inference: no checkpoint, no TeloHAEC control mean or context
      text in PIE's format, no cost (frozen l.122–126).
    - "PIE never saw TeloHAEC" (l.156) is asserted without a citation. 700's dataset list (tahoe, jiang, arc_vcc_25,
      orion; `experiments/700-pie-transfer/RESULTS.md` l.49) is consistent with it, but the claim should cite that.
  - So as written, TRUSTED (3/3 on ABC/telohaec) cannot be reached, and DIAGNOSTIC (2/3) is reached by default.

  → R2.
- P3 (frozen l.162–164) is a useful blob sentinel, but its comparison is confounded (see check 3).

### 3. Is it measurable on the checker's surfaces? **The bulk members are. The DE members are correctly pinned. The fold build is under-specified.**
- **Members:** bulk only (pds literal + pdsl, and mse through r). DE stays with DE-H1 / DE-CD4 (frozen l.178–179), which is
  `experiments/800-checker/DESIGN.md` §4 l.65–72 and `PHASEB.md` l.6–16. That is acceptable for an instrument: it prices no
  arm of its own.
- **The operator does not exist (frozen l.125: "801 dumps for L1–L6/telohaec… with the production operator
  (`--fold-ops prod`)").** In 801, `prod` means `ctx.full_op`, the panel's 3-donor operator
  (`experiments/801-lab-folds/lib801.py` l.334–346). Every L-panel dump was built with `--op none` and per-held-line fits only
  (`experiments/801-lab-folds/cloud801.sh` l.102–109; `work/L1/ops/` holds only `held_*.json`). So no full-donor operator
  exists for any L panel.
  - Stage 1 therefore needs six new full 3-donor operator fits, one per L panel. They are not in the work column, the time
    or the cost.
  - Fits are not bit-reproducible across machines or thread counts (`experiments/801-lab-folds/RESULTS.md` l.30). Each needs
    a recorded sha and a reproducibility rule.

  Only ABC/telohaec has its operator: production's, as for the ABC/kolf2 fold (801 RESULTS l.326). → R3.
- **The engine has no `telohaec` entry.** `FOLD_SOURCES` (`lib801.py` l.79–84) and the bank registry need one, presumably
  MAIN, as for kolf2. The proposal names "the engine's KAs" (l.125). It must also say that the change leaves every existing
  fold's stored rows bit-identical, and how that is checked. → R3.
- **P3's level comparison is confounded by donor count.** A TeloHAEC fold priced with MAIN has three cancer donors. Every
  cancer fold has two (`lib801.py` l.79–82). "Below every cancer fold's" (l.162) therefore mixes cross-lineage difficulty
  with a 3-vs-2 donor stack. → R4.

### 4. Does it use any seen data as evidence of gain? **No gain is claimed. One leakage rule is missing for when the folds are used.**
- The folds price nothing of their own. The A/B/C board reads enter only as known-answer calibration, the protocol's use
  (DESIGN §7 l.86–94). No sealed reserve is touched: L1–L6 and ABC only. LZ is not named (frozen l.12), and that is correct.
- **Missing rule:** 811's day bank list already contains the TeloHAEC pilot as a voting line (`experiments/811-lolo-model/PROPOSAL-2.md`
  l.57). Use (b) would add the genome-wide rows. On a TeloHAEC-held fold, any candidate that carries TeloHAEC rows (pilot or
  genome-wide) must drop them, as PLOLO does. Otherwise the fold is a seen surface for that candidate, which DESIGN §7 l.93–94
  (702 rule 2) allows to bound a cost but never to supply a gain. The proposal runs (a) and (b) together and does not say
  this. → R5.
- Stage 1 builds bank rows for all 19,276 targets. That includes the targets of the sealed LZ panels
  (`docs/CONFIRMATION-LEDGER.md`, "Lab reserve LZ"). A bank is not a dump or a score, so DESIGN §10 l.121–123 is not broken.
  But nobody in the lane may build or read a TeloHAEC fold on an LZ panel. → condition C2.

### 5. Is the try budget honest? **Yes.**
- 0 configurations priced on BULK-W by this lane (frozen l.15–16, l.189–190).
- The estimator choice (2) is made on simulation only (l.187). Controls, gate and stage-0 targets are fixed (l.188). Nothing
  is re-tuned after a failed gate (l.190).
- The 811 configuration is counted in 811's family (l.165), as 811 PROPOSAL-2 l.157 already provides.
- The wear entries go to 800's ledger as `surface=telohaec` (l.191–192).
- Under R1, a redefined KA-SG must be fixed before any real row is read, and calibrated on KA-SIM's twin. That keeps the
  estimator count at 2.

### 6. Is the expected gain worth its cost? **Yes, once the plan is complete.**
- It is an instrument, not an arm, so there is no board gain to set against the class noise.
- What it buys:
  - fold lines beyond the one-lab Orion pair, whose mutual partnership is the strongest in 811's ablation (+0.011…+0.015,
    `experiments/811-lolo-model/RESULTS.md` l.291);
  - near-full coverage on every panel (verified here: in library 300 / 298–300; ≥ 100 rep-1 incidences on 298 / 295–299 per
    panel, `work/meta/telohaec_rep1_coverage_v2.json`);
  - a same-day fourth fold line on D/E/F's real panel (frozen l.139–141).
- At ≤ $3 for the kill gate and ≤ $20 in all, the price is small against the lab's $800.
- The cost table must grow by the six operator fits (R3) and by the PIE inference, or by its removal (R2).

## What is missing (REVISE: fix each item, resubmit once)

- **R1 — Re-specify KA-SG so that an unbiased estimator passes it.**
  - Make the single-guide cells disjoint from A's carriers and background for that test.
  - Replace the OLS slope with a statistic that noise does not attenuate. For example, a reliability-corrected slope
    ⟨A, S⟩ / ⟨A_h1, A_h2⟩-type ratios, or an errors-in-variables (Deming/total least squares) slope with the noise variances
    from split halves.
  - Run the exact statistic on KA-SIM's simulated data first (the real assignment, the single-guide subset) to show its pass
    rate under an unbiased truth, and set the [lo, hi] bar from that, before any real row is read.
  - State what the cross-target-mean condition can and cannot see (it is blind to μ̄).
- **R2 — Make P2 able to fail.** Either:
  - add a PIE inference stage for TeloHAEC: checkpoint (and why it is unseen, with the citation), the TeloHAEC control mean
    and context text in PIE's input format, instance, cost, and the sha of the output rows; or
  - drop 702 from ABC/telohaec, restate P2 on the rows that can be built, and add at least one board-pointed
    discriminating negative the engine can express on a TeloHAEC fold (DESIGN §9 l.105–115 lists the candidates; name which
    are expressible). Or say plainly that the fold's trust then rests on the two positives alone, and lower its bin names
    accordingly (no "TRUSTED" on two near-certain positives).
- **R3 — The fold build, completely.**
  - For each of L1–L6/telohaec, name the operator: a new full 3-donor fit per panel with the 801 day recipe
    (`--basis-prune essential`, the 12k cap). Give the fit job, time, cost, sha recording, and the reproducibility rule.
  - For ABC/telohaec, the production operator's sha.
  - The `lib801` / registry change for the new fold line (sources = MAIN), and a KA that every existing fold's stored rows
    stay bit-identical after it.
- **R4 — Remove P3's donor-count confound,** or drop P3 as a blob sentinel. For example, also price the incumbent on each
  TeloHAEC fold with each 2-donor cancer subset, so the level comparison against the cancer folds is like for like, and state
  the P3 bins on that.
- **R5 — The decision role and the leakage rule, preregistered.**
  - **Leakage:** on a TeloHAEC-held fold, every candidate's TeloHAEC rows (pilot and genome-wide) are dropped, and check800
    refuses a gain read otherwise (a declared-seen surface).
  - **Decision role:** say how the TeloHAEC folds enter a check800 decision: reported apart like kolf2, a veto, or pooled
    into BULK-W. Say what a candidate that gains on the cancer folds and loses on TeloHAEC gets. All six challenge lines are
    cancer lines (`AGENTS.md` l.23) and TeloHAEC is not, so "a stricter transfer test" (frozen l.23–25) is not automatically
    "a more D/E/F-like test". The rule must not let a non-cancer line overturn a cancer-fold read without a stated reason.
- **R6 — Novelty: cite 105 and answer it.**
  - For (a): why a fold this far from every donor (105: TeloHAEC rows are 0.13–0.14 pds below K562 as a source for the
    Orion lines) still orders arms like the board. P2/R2 is the test; say so.
  - For (b): name the mechanism that differs from 105's direct-source construction, or leave (b) wholly to 811's own judged
    try and delete P4 here.
- **R7 — Data terms as a hard gate.** "The consortium's external-use terms are re-read at stage 0, before any fetch"
  (frozen l.208–209) becomes a STOP gate with a recorded file (the terms text, sha256) and the reading in RESULTS.md. If the
  terms forbid this use, no byte is fetched.

## Conditions that will carry into a GO (stated now so the resubmission can meet them)
- **C1:** stage 0 runs alone. Stage 1 is launched only after every stage-0 gate passes and the lab lead has read the
  stage-0 RESULTS.
- **C2:** no TeloHAEC fold, dump or score on an LZ panel. Bank rows for LZ targets may exist on S3, but they are not read.
- **C3:** the estimator (A or B) and the KA-SG statistic are frozen by sha before the first real row is read. A failed gate
  is recorded as a negative with its mechanism, with no re-tune. The DLD-1 fallback (frozen §9) is a new proposal, not an
  automatic branch.
- **C4:** cloud via `labcloud.sh` only, on-demand, Project=lab, ≤ $50 per launch, the lane cap as resubmitted. Nothing
  > 200 MB to the Mac, and no Mac transfer at all in a WINDOW.

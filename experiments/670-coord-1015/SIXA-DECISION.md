# Which replacement final entry (6a) to ship after R0: the decision rule

Written 2026-10-19 09:31Z by the coordinator [223d39], BEFORE any D/E/F data exists and before C1's board read (10-20 s1).
How a candidate is BUILT and GATED is in `experiments/690-def-readiness-1015/OCT22-CHECKLIST.md` §C 6a (a)–(f); this file fixes only
WHICH candidate goes into the single replacement. The board is dark for final entries, so this is decided entirely on evidence
fixed in advance. Only the last final entry counts; R0 is the default and stays final unless a candidate clears everything below.

## Inputs (all produced without the coordinator's judgement)
1. **708's day-of decision** on D/E/F's own folds: `experiments/708-def-panel-repricing/work/decision708-d708def.txt`, written by
   wait708 (pid 66950) after the release and R0, with R705=CLEARED, R707=WEAK, R710=POSITIVE. It says whether the **swap** (707's
   hook, AMENDMENT 2's relaxed rule) turns ON and whether the **stack weights** (H1 / CD4) move.
2. **A/B/C composition reads already on the board:** R-ABC2 (R0's recipe) COMPOSES (D −0.00006); **R-ABC2H** (R0 + swap)
   **COMPOSES** (D −0.00012, entry ndt8swt7). There is no A/B/C composition read for a weights move (708 has its own guards, 710 = no block).
3. **C1** (R0 + swap + 811 c0_l25) on the board, 10-20 s1: AS PREDICTED / ABOVE / HARM against 0.2757814 (lab PREREG-C1-ABC, 333dddc8).
4. **The lab's 811 dossier** (docs/LAB-DEF-RUNBOOK.md Amendments 1–3): DEFo 2.5 % bound > 0 and DEF point > 0 for the increment over the
   branch that ships (C1 − B1 if the swap is ON, C0 − B0 if OFF), no Arc-drawn panel HARM (C1's board read counts), K / DE-H1 clean, and
   one LZ sign check, routed through the coordinator into CONFIRMATION-LEDGER.

## Candidates and eligibility
- **X0 = R0.** Always eligible; the default.
- **X1 = R0 + 708's decided changes** (the swap if ON, plus any weights move 708 decides). Eligible iff 708's decision moves something
  AND, for the swap, R-ABC2H did not read BELOW (it COMPOSES). Built exactly as 697 (R0's recipe + DEF_DAY_DUMP_HOOK=dumphook707.sh)
  for a swap-only X1, or through `october_recipe.py decide …` with decision708 as evidence for a weights move (checklist 6a (a)).
- **X2 = X1 + 811** (or R0 + 811 if the swap is OFF). Eligible iff X1's conditions hold for its base, AND C1 is not HARM, AND the lab's
  dossier passes its preregistered rule on the branch that ships, AND the LZ sign check agrees.

## Choice
1. Build in order of simplicity: X1 first, then X2 if eligible (each ~4–6 h on def_day; there is time until 11-04).
2. Ship the most complete candidate that passed EVERY checklist 6a gate (a)–(f): X2 if it passed, else X1 if it passed, else nothing
   (R0 stays final). A failed gate on X2 does not block X1: the "one attempt" is one SUBMISSION, not one build.
   **No fix-and-retry** (agreed with readiness, OCT22-CHECKLIST 6a (d)): each candidate is built at most once. A build interrupted
   before its gates (reboot, OOM, crash) may be rebuilt with unchanged inputs. A gate failure shown to be a judge false-FAIL (as 697 / 698)
   is re-gated on the finished build, not rebuilt. A real gate failure retires that candidate.
3. Never ship a candidate whose increment has only a board read and no D/E/F-fold read, or only a fold read and no A/B/C board read:
   the X1 swap has both (R-ABC2H + 708), and X2's 811 has both (C1 + DEFo/DEF).
4. Timing: decide within 24 h of 708's decision file and the lab's dossier; submit by 10-27 if possible; never start a replacement upload
   after 11-04 12:00Z. The submit goes through arm1015 with a screened, claimed name and a name_guard beside it (6a (f)).
5. Record the decision, its inputs (file shas, entry ids) and the try count in docs/PLAN, docs/NEXT-SLOT-DECISION and this file, before
   arming the submit.

## Amendment 1 (2026-10-20 02:35Z by the clock; before any D/E/F data, after C1's read and before any C0 number)
- **Input 3 is filled:** C1 (RMxJtslBpY2GkwIMHOiB) read 0.2777649 vs pred 0.2757814, D +0.0019835 → **ABOVE, not HARM**. On the swap-ON
  branch, X2's A/B/C condition holds.
- **Gap found in rule 3:** if 708 turns the swap OFF, X2 = R0 + 811 = C0. Its increment C0 − B0 has DEFo/DEF fold reads but no A/B/C board
  read of its own (C1 carries the swap). By rule 3 that X2 would not be eligible. C1's read is not stretched to cover it.
- **Fix:** a board read of **C0** (R-ABC2's recipe + 811's c0_l25 on the unswapped stack dump, def_day / A/B/C) in one of the three slots
  left (10-20 s2, 10-21 s1, 10-21 s2). The lab writes the hook, its KAs and the prediction before the build; readiness builds it as 699.
  HARM if D ≤ −0.0005 against that prediction, or if C0 < R-ABC2 (0.2739849) − 0.0005 = 0.2734849. HARM closes X2 on the swap-OFF branch
  only. Not HARM fills input 3 for that branch.
- **If C0 is not read before the release** (no READY, a failed build, a lost slot), X2 on the swap-OFF branch stays ineligible. That
  branch then ships X1 (if 708 moves the weights) or R0. Nothing else changes: the swap-ON branch is decided by C1 as before.
- Hard stop: no def_day build for C0 starts after 10-21 12:00Z, so the machine is clean for watch691's R0 on 10-22.

## Amendment 2 (2026-10-20 03:39Z by the clock; before any D/E/F data and any W1 number): the H1 weights guard (673 W1)
- W1 = R0's exact recipe + DEF_DAY_H1_STACK=x2.5 (CD4 ×1), on def_day / A/B/C, 10-21 s1 (673 PREREGISTRATION 5c4f4bc5; reader `read695.py --w1`).
- **If W1 reads BELOW** (D ≤ −0.0005 vs pred 0.2744673): no 708 move may raise H1 above ×2 in the final. Mechanism: readiness re-arms wait708 with
  R710=HARM before the release, so decide708 applies the block itself. That is safe because wait708's log is nobody's ready-log. At decision time the
  coordinator also checks that decision708 shows no H1 > ×2 (backstop). COMPOSES / ABOVE changes nothing.
- **Coverage:** a Part A move to H1 ×2.5 alone has W1 (R0's recipe) + 710 (ALLSW). A move to H1 ×2.5 + CD4 ×1.5 (the 10-16 rehearsal's pick) has
  W1 for the H1 part and 457dr (+0.0031, October path) for the CD4 part, with no joint read. Rule 3 counts component reads for a weights move. Any CD4
  ×1.5 move still needs AMENDMENT 1's fold confirmation (708), and 438's October-like panels ruled CD4 ×1.5 out, so that is where it is guarded.

## Amendment 3 (2026-10-20 04:32Z by the clock; before any D/E/F data and any W1/W2 number): the CD4 weights guard (674 W2)
- W2 = W1 + DEF_DAY_CD4_STACK=x1.5 (the 708 rehearsal's pick), 10-21 s2, read against W1 (674 PREREGISTRATION 12562add + Amendment 1).
- **If W2 reads BELOW:** X1 may not carry a CD4 dose above ×1, whatever decision708 says. The coordinator applies this when choosing the 6a
  candidate; decide708 is not changed. If decision708 moves CD4 above ×1, X1 is built with CD4 ×1 and only its other changes. If that leaves
  nothing, X1 is not eligible.
- COMPOSES / ABOVE: no change (CD4 ×1.5 is A/B/C-positive by 438; only a BELOW informs the final).

## Input record (2026-10-20 07:31Z): C0 read
- **C0 (4rq3qHIuv8bHAxpTkUeu) 0.2780624 vs pred 0.2763696, D +0.0016928 → ABOVE, not HARM.** Amendment 1's condition holds. The swap-OFF X2 (R0 + 811)
  is eligible on its A/B/C read; it still needs the DEFo / DEF rule on C0 − B0 and the LZ sign check.
- Information only, not a rule change: on A/B/C, C1 − C0 = −0.0002975. The branch is still chosen by 708 on D/E/F's folds.

## Amendment 4 (2026-10-20T07:39Z by the clock; before any D/E/F data): on a swap-ON day, the swap must earn its place on top of 811
- **Trigger** (readiness flagged it after C0's read): on A/B/C the bulk members, which are exact, give C0 − R-ABC2 = +0.0038118 and C1 − R-ABC2 = +0.0035022.
  Once 811 is in, the swap adds −0.00031. ABCp's information read says the opposite: C1 − C0 = (C1 − B1) + (B1 − B0) − (C0 − B0) = +0.00153 + 0.00169 −
  0.00269 = +0.00053. 708's swap rule prices the swap ALONE over R0 (B1 − B0), never on top of 811, so nothing fixed in advance answers this for D/E/F.
- **Rule** (applies only when 708 turns the swap ON and X2 is eligible on the swap-ON branch, i.e. C1 − B1 passed the lab's bars):
  1. The lab's dossier reports **C1 − C0 on DEFo** (pooled point and per fold), from the same composite arms (C1:C0 is one more aggtag811 pair).
  2. If that point is **> 0 on ≥ 2/3 folds**, X2 = C1 (R0 + swap + 811), as before.
  3. Otherwise, if C0 passes its own bars on DEFo/DEF (C0 − B0: pooled 2.5 % bound > 0, DEF point > 0), **X2 = C0 (R0 + 811, no swap)**. The swap then
     ships in no candidate except X1 = R0 + swap, the fallback if X2 fails its gates.
  4. Otherwise X2 = C1 (the base rule; this case needs C0 − B0 weak while C1 − B1 passes, which a positive swap makes rare).
- **Why this is not board-mining:** the A/B/C pair only raises the question. The choice is made on D/E/F's own folds, with bars fixed now. Both
  constructions already have an A/B/C board read (C1 ABOVE, C0 ABOVE), so rule 3 holds for either. Ties go to the candidate with fewer changes.
- **Build count:** still one X2 build (C1 or C0, chosen before building), so "one attempt" and "no fix-and-retry" are unchanged. The LZ sign check
  reads on the increment of the construction that ships (C1 − B1 or C0 − B0).
- Tries: this adds one comparison (C1 vs C0) to the swap-ON branch; the swap-OFF branch is unchanged (X2 = C0 on its own bars).

## Amendment 5 (2026-10-20T09:31Z by the clock; before any D/E/F data): X2 only on the base 811 was priced on
- **Gap:** the lab prices 811 on B0/B1 = R0's bulk, i.e. weights h2c1 and 705's θ rider ON. If 708 moves the weights or VETOES the rider, "X2 = X1 + 811"
  would ship 811 on an unpriced base. Neither the folds nor the board measured that combination.
- **Rule:** X2 is always the PRICED construction: **R0 + 811 (C0), or R0 + swap + 811 (C1)**, chosen by the swap branch and Amendment 4. It never takes 708's
  weights move or rider veto. X1 = R0 + all of 708's decided changes, as before.
- **When nested** (708 moves neither the weights nor the rider), nothing changes: X2 ⊇ X1, so ship X2 if it passes, else X1, else R0.
- **When not nested** (a weights move or a rider veto): ship X2 if it passes every 6a gate and its 811 rule (DEFo, DEF, LZ), else X1, else R0. Reason, fixed
  now: 811 carries the most verification of any candidate (two A/B/C board reads ABOVE, DEFo + DEF + LZ). A weights move has fold reads plus component
  board reads only. Still one submission, each candidate built at most once.
- The mechanical form of the whole rule (Amendments 1–5) is `experiments/670-coord-1015/decide6a.py`. It prints what to build, in order, with the env
  changes over R0's recipe, and the eligibility reasons.

## Input record (2026-10-20T11:34Z): W1 re-gated, W2 dropped
- W1 (673) false-FAILed gate673 on bulk_sharpen_pt's data-derived fields (an H1 dose moves the key ratio's normalisation). It was re-gated on the finished build
  (regate673.py, READY w1 1337cc9b) and is armed for 10-21 s1. Amendment 2 is unchanged (w1route.sh applies a BELOW).
- W2 (674) was not built. **Amendment 3 never applies** (no W2 read): a 708 CD4 move is guarded only by 708's own AMENDMENT 1 fold confirmation. decide6a.py
  reads W2 as NOT READ, which blocks nothing.

## Amendment 4b (2026-10-20T15:20Z by the clock; before any D/E/F data): the swap must CLEAR a bound on top of 811, not just lean positive
- **Trigger** (the lab, 11ef50f7 follow-up): on both Arc-drawn panels with board truth, the fold instrument reads C1 − C0 slightly positive (LABC +0.00023 on 2/3
  folds; ABCp ≈ +0.0005 on 2/3 by difference), while the board's exact bulk members read −0.0003. Amendment 4's "point > 0 on ≥ 2/3 folds" would keep the
  swap where A/B/C truth says it adds nothing. The instrument leans positive on this comparison.
- **Rule (replaces Amendment 4 item 2):** X2 = C1 only if **C1:C0 on DEFo has ci95[0] > 0** (literal-pds, the same bootstrap as the runbook's bars). Otherwise X2 = C0 if
  C0 − B0 passes its own bars, else C1 (item 4 unchanged). On both calibration panels this picks C0, the board's answer.
- Why it is allowed now: it is calibrated on A/B/C instruments against A/B/C board truth, before any D/E/F number exists, and it can only REMOVE a change
  (burden of proof on the extra lever). decide6a.py is updated to match (test_decide6a).

## Input record (2026-10-21 00:55Z): W1 read
- **W1 (WULHAm2fytPtmAb3ryVN) 0.2746084 vs pred 0.2744674, D +0.0001410 → COMPOSES.** Amendment 2 does not fire: a 708 H1 move above ×2 stays eligible for X1
  (decide6a reads W1 from submit-1021s1.log). All A/B/C inputs are now fixed: R-ABC2H COMPOSES, C1 ABOVE, C0 ABOVE, W1 COMPOSES, W2 not read.

## Input record and X1 ruling (2026-10-22T06:11Z, after R0 and day708; before any X1/X2 build)
- **R0** = final entry k2OKr0zPvJq48BCtdUx8 (published 06:05Z, members hidden; package 5c82ea9e, identity fa832988).
- **708's decision** (decision708-d708def.txt, sha256 ff8243e2d7af9e96…, DECISIONS WRITTEN 06:09:50Z): WEIGHTS → h3c1.5 (H1 ×3, CD4 ×1.5; confirmation under the production key
  +0.00177, floor-kept +0.00172 → MOVE), RIDER 705 ON (r +0.00032, 3/3), SWAP707 ON (s|r +0.00182, 3/3, lt3 +0.00075; floor-kept +0.00169). 710 POSITIVE → no block.
- **decide6a.py** (all board bins as recorded: R-ABC2H COMPOSES, C1 ABOVE, C0 ABOVE, W1 COMPOSES, W2 NOT READ): swap branch ON; NOT nested (weights move);
  **X1 eligible = R0 + DEF_DAY_H1_STACK=x3 + DEF_DAY_CD4_STACK=x1.5 + DEF_DAY_DUMP_HOOK=experiments/707-strong-tier-direction/dumphook707.sh**. X2 waits for the lab's
  DEFo/DEF dossier; by Amendment 5 it is the priced C1 or C0 alone, without the weights move. Build order: X1 now. Ship order: X2 (if eligible and passing), else X1, else R0.
- Coverage: H1 ×3 above W1's ×2.5 is covered by Amendment 2's rule (W1 COMPOSES → no block on H1 > ×2; 708 AMENDMENT 2's 710 condition is likewise "H1 above ×2").
  CD4 ×1.5: Amendment 3 inactive (W2 not read). 708's D/E/F held-out confirmation is the preregistered evidence, and it overrules 438's October-like history.
- **6a (d) ruling on R0's known day trips, fixed before X1 builds:**
  1. DEF_DAY_CD4_THRESHOLDS=relaxed is part of every D/E/F build's day recipe from the start (the panel maps 287/300). It is not a retry.
  2. Build X1 WITHOUT DEF_DAY_ALLOW_LOW_DOSE. If postbuild's ONLY failing row is 554's and every lambda is ≥ 0.85, re-gate the FINISHED build with the override, exactly as R0
     (554's preregistered 0.85–0.9 rule: a judge waiver on a finished build, not a fix-and-retry). Any lambda < 0.85, or any other real FAIL, retires X1.
  3. If def_day ends with no safe build, X1 retires (R0 stands; no plan B for a replacement).
  4. The config diff vs R0's meta may contain exactly the weights' keys (bulk_load, bulk_load_sha256, bulk_ctx_scale, and bulk_sharpen_pt with data-derived fields only:
     the W1 lesson) plus the hook export. r0 stays on disk as the reference.

## Input record (2026-10-22T07:01Z): the lab's 811 D/E/F dossier (ad212d84) → X2 = C0, pending LZ
- All KAs pass on DEFo and DEF. K / DE-H1 clean = YES (C0 max |dK| vs inc 4.46 %; C1 6.09 % on HCT116, covered by the DE-H1 read and both A/B/C board reads).
- **decide6a** on the 8 agg811-811-composite jsons (DEFo + DEF; pairs C1:B1, C0:B0, C1:C0), --lab-clean YES, --lz NONE, verified by the coordinator:
  - C1:B1 DEFo +0.00299 [+0.00054, +0.00566] 3/3, LT3 +0.00142; DEF +0.00222 → the swap-ON branch passes §5 items 1–2.
  - C1:C0 DEFo +0.00033, **ci95[0] −0.00071**, 3/3 folds → the swap does NOT clear its bound on top of 811 (Amendment 4b). Under the old 4 ("> 0 on ≥ 2/3") this would have
    kept C1; 4b, set before any D/E/F number, picks C0.
  - C0:B0 DEFo **+0.00446 [+0.00175, +0.00733] 3/3, LT3 +0.00291; DEF +0.00363** → C0 passes its own bars.
  - **X2 = C0 = R0's day recipe + DEF_DAY_DUMP_HOOK=experiments/811-lolo-model/dumphook811c0.sh** (no swap, no weights move; Amendment 5's priced construction).
    lz_increment C0 − B0. Eligible once LZ reads PASS (pooled point > 0).
- Plan: X2 builds right after X1's def_day ends (one builder), with the same 6a (d) ruling (relaxed CD4 from the start; 554 re-gate only if it is the only FAIL and every
  lambda ≥ 0.85). Ship order unchanged: X2 if eligible and passing every 6a gate, else X1, else R0. One submission.

## X2 bank gap and the conditions for GO-X2 (2026-10-22T07:05Z; before any X2 build)
- Readiness's read-only preflight of dumphook811c0.sh on r0's stack dump (ada1baa2) with the D/E/F bundle REFUSED: "136 targets were never extracted into 811's banks"
  (log experiments/690-def-readiness-1015/work/def-day/hook811c0_on_r0/hook.log). The pinned banks manifest (37e34335…) covers panels ABC + L1–L6 only. The lab's DEFo/DEF
  composite runs built their own rows on EC2, so the priced numbers exist but the Mac hook cannot reproduce them yet.
- **Ruling:** X2 is still Amendment 5's priced construction **iff the hook builds exactly the rows that were priced.** Banks are data inputs (target rows from public donors), not
  design choices. The method, λ, bank list and centring stay as frozen. Conditions, ALL required before the coordinator touches experiments/670-coord-1015/work/GO-X2:
  1. the lab extends 811's banks to D/E/F's 300 targets with the SAME pinned extraction code (ext811, unchanged), commits, and re-pins BANKS_MANIFEST_SHA in hook811.py
     (the only code change; reported with its diff);
  2. a known answer: the hook's C0 bulk rows on D/E/F's real panel == the DEFo kit's C0 construction for the same targets (|d| = 0, as 017cc763 on ABCp), and the hook
     still reproduces its old A/B/C known answers (λ 0 → input content; λ 0.25 on 2dfe327d → 96204332);
  3. readiness's preflight of the re-pinned hook on r0's stack dump succeeds (its output = X2's content reference);
  4. **LZ reads PASS** on C0 − B0.
  Any one failing closes X2; X1 then ships if it passes its gates, else R0 stays. No other change to 811 is allowed.
- 2026-10-22T07:25Z: **the lab's bank plan accepted** (ETA ≈ 08:45Z, EC2). (a) 801's kit811ext.py with KIT811_PANELS=DEF, unchanged: the code that built the priced DEFo/DEF rows; old rows kept byte for
  byte, D/E/F rows appended. (b) One allowed deviation: a ~30-line format adapter defbanks811.py. kit811ext's manifest has no "out" key, so the adapter copies the npz files BYTE
  FOR BYTE into work/banks and writes ext811's manifest schema with the DEF panel added. Condition: it lists each copied npz's sha256 == the extension's output. hook811.py
  changes only BANKS_MANIFEST_SHA. Both diffs go to the coordinator. (c) KAs: kasame811c0 (hook == composite811 C0, |d| 0) on DEFo ×3, DEF ×3, ABCp ×3 and L1; composite811 on DEFo
  with the new banks == the stored priced DEFo arrays (|d| 0); dumphook811c0 λ 0.25 on 2dfe327d → 96204332 and λ 0 → input content. All of this satisfies GO-X2 conditions 1–2.
- 2026-10-22T07:49Z: **GO-X2 condition 4 holds: LZ PASS** on C0 − B0 (+0.00999 [+0.00659, +0.01273], 6/6; CONFIRMATION-LEDGER). The lab's runbook §5 items 1–5 all hold for X2 = C0. Still
  needed for GO-X2: 1 (re-pinned banks + diffs + npz shas), 2 (the KAs) and 3 (readiness's preflight of the re-pinned hook on r0's stack dump).
- 2026-10-22T08:10Z: **GO-X2 conditions 1 and 2 (EC2) hold** (lab eb79f7f5, verified by the coordinator): hook811.py's ONLY change is the BANKS_MANIFEST_SHA line 37e34335 → 4cc6d82a (hook now
  c7cff26d); the manifest on disk hashes to 4cc6d82a with panels ABC, L1–L6, DEF; dumphook811c0.sh unchanged (297c3832); committed, clean; old rows byte-identical; measured
  masks unchanged. KAs: kasame811c0 (hook == composite811 C0) PASS on DEFo ×3, DEF ×3, ABCp ×3, L1; composite811 with the new banks reproduces the priced DEFo/DEF arrays
  (402/402); on EC2, old and new banks give byte-identical A/B/C output (ba7a6477, ≤ 1.3e-15 from the Mac's 96204332 = platform float). **Still open:** the Mac file-sha KA
  (→ 96204332) and readiness's preflight on r0's stack dump. Both were asked for at 2026-10-22T08:10Z.
- 2026-10-22T08:11Z: **the adapter condition holds, verified independently**: recomputing sha256 from disk, all 18 banks in work/banks == kit811ext's out_sha256 == the list's copied column, and
  all 18 in work/banks.pre-def-37e34335 == the pinned 37e34335 manifest == kit811ext's merge base (npz-sha-list.txt, a91cbe93; defbanks811.py c52d9b35 asserts the same).
  Still open for GO-X2: the Mac file-sha KA (→ 96204332) and readiness's preflight.
- 2026-10-22T08:13:31Z: **GO-X2 WRITTEN** (experiments/670-coord-1015/work/GO-X2). All four conditions hold: (1) banks re-pinned (hook811.py c7cff26d, manifest 4cc6d82a committed at 77c50953,
  adapter 18/18 verified); (2) KAs incl. the Mac file-sha KA (λ 0.25 on 2dfe327d → 96204332 byte for byte; λ 0 → input; lab c1d35afa); (3) readiness's preflight PASS 08:10:40Z
  (hook on r0's stack dump → 3b4f0c05; carriers / perts / genes / env equal); (4) LZ PASS. X2 builds after X1's def_day in readiness's chain_x1x2 (its own preflight + 6a gate).

## X1 result and the gate6a judge fixes (2026-10-22T10:14Z)
- **X1 READY** 3528fb51 (submissions/x1def.vcc; ready log experiments/692-def-r0/work/x1g-ready/x1def.ready.log): CONVERGED at lvl2 (factors D 1.0835 / E 0.9557 / F 1.0681, rho
  0.9783 / 0.9715 / 0.9703; lvl0 deleted on E 1.173, lvl1 on F 1.0745). Postbuild failed ONLY 554's row, D prior_lambda 0.8975 (R0's value), so the (d)(2) re-gate ran (KA reproduced 42
  rows; the override flipped only 554; zero-gene 0; prep rc 0; novelty OK). **gate6a X1 PASS** (gate6a_x1.json; judge 360f3125, commit 777b6db6).
- **Ruling on readiness's three judge fixes, made before the PASS** (none touches the build): (1) relative paths resolved: a crash fix, accepted. (2) bulk_load compared by
  sha, not file name: stricter, accepted. (3) line (e): the pinned builder embeds its work dir and each cd4_compact meta.json sha, and those metas record a MEASURED peak_rss_gib.
  New rule: identical after normalising those; the metas may differ only in peak_rss_gib; every other implementation entry byte-equal. This is a judge false-FAIL on provenance and
  measurement (as 698's float floor), and the R0-vs-R0 negative control still fails every candidate row. Accepted. **Condition: X2 is judged by the same committed gate6a
  (360f3125), with no further edits.**
- **X1 is HELD as the fallback, not submitted.** Ship order X2 > X1 > R0, one submission. X2 launched 10:11:55Z (the chain's preflight reproduced 3b4f0c05); ETA ≈ 15:00–16:30Z.
  If X2 fails a real gate, or never reaches READY, X1 is submitted.

## THE 6a DECISION (2026-10-22T17:49Z): ship X2 = C0
- **X2 READY e39d2ba6** (submissions/x2def.vcc; x2g-ready): CONVERGED at lvl3 (D 1.0757 / E 0.8924 / F 1.0150; rho 0.9767 / 0.9757 / 0.9781). Postbuild failed ONLY 554's row (D prior_lambda
  0.8994 ≥ 0.85), so it was re-gated per (d)(2). **gate6a X2 PASS** with the same judge (gate6a_x1.py 360f3125, unchanged since 777b6db6). Diff vs R0 is exactly {bulk_ctx_scale,
  bulk_load_sha256}; stack dump == R0's by content; the used dump == the 811 hook on R0's stack dump by content (3b4f0c05); same tree. The gate's NOTE on bulk_sharpen_pt is ~1e-15
  float noise (edges / key_ratio_c), with behavioural fields identical (counts [100, 102, 98], moves, carried-key sha): accepted under the W1 rule.
- **decide6a** with every input (board bins; decision708; the 8 DEFo/DEF jsons; LZ PASS; lab clean YES; X1 PASS, X2 PASS): **SHIP = X2.**
- Armed 2026-10-22T17:49Z as slot R1: arm1015 R1 x2def (ready x2g-ready/x2def.ready.log, run x2g-ready/x2def.run.log), name from the queue ("LACK OF COURAGE IN THE MODERN SOFT WORLD IS A
  VERY BIG OPPORTUNITY."; the reserved final name is kept for any Nov 5 re-upload), prior = prior_relay stub of R0, cap 11-04 12:00Z, name_guard R1. **X1 is retired, never
  submitted.** One replacement submission. Tries for the 6a switch: candidates built 2 (X1, X2), submitted 1.

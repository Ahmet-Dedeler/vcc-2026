# Lab runbook for the D/E/F release (10-22) — and 811's 6a preregistration

Written 2026-10-17 22:33Z by the clock, BEFORE any D/E/F data exists. The lab lead runs it. Production owns `data/challenge-def/` (the lab
only reads it), the slots and the recipe. Every step below is EC2 or seconds-long local arithmetic. Each step has its gate,
and a failed known answer stops the step.

## 0. Gates (all must hold before step 1)
- `bash experiments/800-checker/window800.sh` prints **CLEAR**. The release window runs from release detection until R0's
  "SUBMITTED ok", or 10-21 20:00Z → 10-22 14:00Z, whichever is earlier. Inside it: no lab local job, no Mac transfer > 1 GB,
  and at most 2 lab agents.
- `data/challenge-def/pert_counts.csv`, `gene_names.csv` and `context_*.h5ad` exist (production's download). Record their sha256
  in LAB-LOG before reading anything.
- Same-panel branch: if D/E/F's targets equal A/B/C's, P-A / P-B are VOID and steps 2–4 still run. The D/E/F folds are then
  LABC with D/E/F's frames.

## 1. The free predictions (seconds, local)
`uv run python experiments/814-panel-rule/predict814.py data/challenge-def/pert_counts.csv` (sha256 `3c7e2658a32e2a8d…`, frozen).
- P-A: essential-screen targets; support ≤ 3, fail ≥ 8.
- P-B: targets in prop10 deciles 9–10; support ≥ 130, fail ≤ 95.
- Log both verdicts. A P-A fail means Arc drew D/E/F differently, so the essential-screen reading of the panel bias must be
  revisited before any v2-based number is used.

## 2. D/E/F-panel folds (EC2, ~16 min each, ~$0.5 each; the two run in parallel)
```
bash experiments/801-lab-folds/kit/kit801.sh DEFo data/challenge-def/pert_counts.csv data/challenge-def/gene_names.csv --frames data/challenge-def
bash experiments/801-lab-folds/kit/kit801.sh DEF  data/challenge-def/pert_counts.csv data/challenge-def/gene_names.csv
bash experiments/801-lab-folds/kit/kit801-collect.sh DEFo ; bash experiments/801-lab-folds/kit/kit801-collect.sh DEF
```
- **DEFo** (D/E/F's own control frames, ~0.66 GB Mac upload) is the decision surface. **DEF** (A/B/C's relabelled frames) is
  reported.
- KA: the engine's fold KAs (per-donor mean == LODO key; PT0 = 0) on every fold.

## 3. Arms on the D/E/F folds (EC2)
`bash experiments/801-lab-folds/kit/make-jobs.sh DEFo` (and DEF), then the two `labcloud.sh launch` lines from
`make-jobs.sh`'s header: 810's sgu1 / 707, and 811's c0_l25 / c1_l25 / src_l25. Their own KAs gate them: 810's recomposed
parts == the dumps; 811 reproduces its L1 known answer. ka804's arm set is already produced by the kit.

## 4. The panel read (seconds, local or EC2)
`uv run python experiments/800-checker/panelrank800.py --fresh A --arc DEFo --arc DEF` (sha256 `d5f9e27c87c98679…`; known answer
passed 2026-10-17: `--arc LABCk` equals LABC on 15/15 cells, max |d| 1.7e-17, 0 rank mismatches). Then
`overread800.py --fresh A --abc …` against DEFo. Question: does D/E/F's panel rank like A/B/C's among the A-only draws
(Arc's selection), or inside the spread (the bias was one draw's idiosyncrasy)? Report it; it decides nothing by itself.

## 5. 811's 6a candidate: the preregistered read (written now)
- **Candidate (frozen, unchanged since 10-16 10:49Z):** c0_l25 = the many-line consensus, day bank list, λ 0.25, centring off.
  `m811.py` sha256 `e00192e644828a0d…`, read through the kit's 811 job (`run-kit811.sh`, unchanged arithmetic).
- **Surface and N:** DEFo, the three cancer held lines, vs inc = pct+P2g, literal pds (pdsl reported). This is the family's
  first read on this surface, so N = 1 there and there is no best-of-N shift. Its selection happened on other surfaces (W,
  Arc-drawn A/B/C), whose wear is counted there.
- **Standing evidence already on record:** ABCp +0.00264 [−0.00003, +0.00561] 3/3; LABC +0.00230 3/3 (both WEAK after their
  N); A-only +0.00638 18/18; DE-H1 through 803's hook: no DE cost (nmae −0.0004 … −0.0010 raw).
- **6a DOSSIER only if ALL hold:**
  1. DEFo pooled 2.5 % bound > 0 (target bootstrap, 4,000 draws, seed 800), with all 3 folds > 0 or the leave-top-3-out mean
     > 0;
  2. DEF (A/B/C frames) point > 0;
  3. no Arc-drawn panel reads it HARM (ABCp, LABC: points > 0, already true);
  4. |dK_rel| on DEFo ≤ 5 % on every fold, or the DE-H1 read stands (it does: no cost);
  5. the sealed LZ confirmation (L7, L8, one access, through the coordinator per `docs/CONFIRMATION-LEDGER.md`): point > 0,
     with each LZ panel's essential-screen count reported (~17 each; arc800-drawn, so expect an over-read).
- **Projection reported, not decided on:** board ≈ DEFo ÷ ρ, with ρ ∈ [1.44, 1.72] (Arc-drawn key / direction, ABCp KAs).
- **Outcomes:**
  - ALL hold → one DOSSIER to the coordinator, who decides the 6a replacement (one attempt, none after 11-04 12:00Z).
  - Point > 0 with bound ≤ 0 → WEAK, no dossier.
  - Point ≤ 0 → the family is closed for the final.
- Nothing is re-tuned after a D/E/F number exists. Any other arm read on DEF / DEFo is report-only and counted in the wear
  ledger.

## 6. Costs and timing
- Steps 2–4 ≈ $3 on EC2 in total.
- Earliest start: R0 submitted (or 10-22 14:00Z), then ~1 h to the panel read and ~1.5 h to 811's read.
- Rough budget check: lab spend $39 of $800.

## Amendment 1 (2026-10-17 22:35Z by the clock, before any D/E/F data; coordinator's correction of the same day)
**Why.** §C 6a allows ONE replacement attempt. The likely attempt is R0 (which carries 705's θ 0.25 rider) PLUS 707's strong-tier
swap, which 708's day path predicts ON. So 811 must be priced as an **increment on top of what would ship**, not against R0 alone.

**New definitions** (on every fold; keys, sharpening and render as each arm defines them):
- **B0** = inc (pct+P2g) + 705's θ 0.25 rider. This is R0's bulk.
- **B1** = B0 + 707's swap (β 0.25, top z_src tertile toward the sources' unit mean, norm-preserving). This is R0 + swap.
- **C1** = B1 with 811's c0_l25 consensus blend (λ 0.25, day bank list, centring off) applied to B1's composition rows. The
  order is: swap first, then 811. Keys are unchanged.
- **C0** = B0 with the same 811 blend: the no-swap branch.

**Decision rule** (replaces §5's items 1–2 and keeps 3–5):
- The **increment C1 − B1 on DEFo** must have a pooled 2.5 % bound > 0 (all 3 folds > 0 or leave-top-3-out > 0).
- **C1 − B1 on DEF** must have a point > 0.
- Also reported: the pair C1 − B0 (swap + 811 together) and C0 − B0 (811 without the swap).
- The dossier names which to ship and why:
  - if production's R0 day path ends with the swap ON, the rule reads on C1 − B1;
  - if the swap is OFF, it reads on C0 − B0, with the same bars.
- N = 1 per branch on DEFo. Both branches are preregistered here, so reading both is not a selection.
- **LZ:** the coordinator confirmed (same message) that the lab keeps the sealed LZ reserve for this one sign check. The claim
  file's sha256 goes into `docs/CONFIRMATION-LEDGER.md` through the coordinator BEFORE the first LZ operator fit. On LZ the sign
  check is read on the same increment as the branch that decides.

**Implementation.** 811 builds the composite arms (B0, B1, C0, C1) as a new script that calls the existing arithmetic unchanged:
ka804's 705 and 707 constructions, m811's consensus, engine801's price. Known answers:
- B0 and B1 reproduce ka804's 705 and 707 arms on existing folds;
- C0 at λ 0 equals B0, and C1 at λ 0 equals B1, bit-exact;
- C0 with the θ rider off equals m811's c0_l25.
Deadline 10-21 18:00Z, before the release window. An information read of C1 − B1 and C0 − B0 on ABCp and LABC (counted in
wear, decides nothing) gives the day its expected size.

## Amendment 2 (2026-10-17 23:14Z): the exact 811 commands (composite arms built and KA'd by 811, `composite811.py` `6ec25788…`)
After kit801 + collect for each tag (DEFo first, then DEF):
```
bash experiments/811-lolo-model/kit/make-jobs811c.sh DEFo
bash experiments/800-checker/cloud/labcloud.sh launch --name lab811c-kit-defo --jobs experiments/811-lolo-model/kit/jobs811c-DEFo.txt \
  --type r7i.4xlarge --minutes 45 --parallel 1 $(sed 's/^/--sync /' experiments/811-lolo-model/kit/sync811c-DEFo.txt) \
  --collect experiments/811-lolo-model/work/run-comp811kit-DEFo
```
- Aggregate with `aggtag811.py` (pairs C1:B1, C1:B0, C0:B0); the exact lines and json fields are in `811/RESULTS.md` §000000.
- The kit job re-runs composite811's KAs on L1, ABCp/k562 and the tag's own folds before it reads anything.
- **Information read** (before the day, counted in wear, decides nothing): C1 − B1 = ABCp **+0.00153 [−0.00096, +0.00428] 2/3**,
  LABC +0.00142 2/3. C0 − B0 = ABCp +0.00269 [+0.00005, +0.00564] 3/3. B1 − B0 (the swap alone) = ABCp +0.00169 3/3. On top of
  the swap, 811's increment is ~55 % of its no-swap size: the swap already moves the strong targets toward the sources' mean.

## Amendment 3 (2026-10-18 08:32Z by the clock, before any D/E/F data and before any C1 number)
**The coordinator's offer:** a board read of C1 (R-ABC2H's recipe = R0 + 707's swap, plus 811's c0_l25) on A/B/C, in 10-20 s2 or
10-21 s1, built by readiness on def_day from the lab's UNAPPLIED `dumphook811.sh`. It runs only if the lab reports READY by
10-19 12:00Z (hook KAs clean).
- **Added to §5's rule (item 3, "no Arc-drawn panel reads HARM"):** C1's board read is HARM if D ≤ −0.0005 against its
  preregistered prediction, or if C1 reads below R-ABC2H's board read by ≥ 0.0005. Either one closes 811's 6a route for the
  final, whatever DEFo reads.
- **Prediction** (written once R-ABC2H's board read exists, before C1 is built): pred = R-ABC2H's read + 811's increment
  deflated, in board overall units. The increment is ABCp C1 − B1 +0.00153 [−0.00096, +0.00428] (LABC +0.00142) ÷ ρ ∈ [1.44,
  1.72] (Arc-drawn ABCp, key / direction), central ≈ +0.0009. Bins ±0.0005 around pred.
- A non-HARM read does not qualify 811 by itself. The 6a decision stays the DEFo read of §5 / Amendment 1.

## Amendment 4 (2026-10-20T08:20Z by the clock, by the coordinator; before any D/E/F data): report C1 − C0 on DEFo
- **Why:** on A/B/C the board's exact bulk members give C1 − C0 = −0.00031 (swap alone +0.00091, 811 alone +0.00408, both +0.00378). Once 811 is in,
  the swap adds nothing there. ABCp's information read implies the opposite sign (+0.00053). SIXA-DECISION Amendment 4 therefore decides a
  swap-ON day's X2 on D/E/F's own folds.
- **What the dossier adds:** the pair **C1:C0** in aggtag811 on DEFo (pooled point + per-fold signs) and on DEF (point), beside C1:B1, C1:B0 and C0:B0.
  Same composite arms, no new build.
- **How it is used (SIXA Amendment 4, tightened by 4b on the lab's calibration note):** if 708 turns the swap ON and C1 − B1 passes §5's bars,
  X2 = C1 only if C1:C0 on DEFo has literal-pds **ci95[0] > 0**.
  Otherwise X2 = C0, if C0 − B0 passes the same bars; else C1. The LZ sign check reads the increment of the construction that ships.
- **Commands (lab, 2026-10-20T14:22Z by the clock):** 811/RESULTS.md §000000's aggtag811 loop now runs C1:B1 C1:B0 C0:B0 **C1:C0** (`--parallel 4`).
  Format test PASS (811 §00000000): the keys `C1|C0|DEFo|pds` / `…|DEF|pds` and `per_fold[].overall` parse in decide6a.py. Pass all eight jsons to `--agg`.
- **Lab caution (an information read, decides nothing):** on both Arc-drawn panels the folds read C1 − C0 as weakly positive (LABC +0.00023, 2/3 folds;
  ABCp ≈ +0.0005–0.0006 by difference, 2/3), but A/B/C's exact bulk members read −0.0003. For this pair, then, "> 0 on 2/3 folds" kept the swap
  on both panels where the board says it adds nothing. The rule is the coordinator's to keep or tighten before D/E/F data.

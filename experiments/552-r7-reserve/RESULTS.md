# 552 — Reserve R7 (running record)

Owner: session "Explore contribution areas". Declaration: `RESERVE.md` (+ amendments 1–2). Ledger: `docs/CONFIRMATION-LEDGER.md`,
"Reserve R7". Purpose: one sealed confirmation of the FINAL October recipe (bulk on R7a, DE on R7b), after R6's last access
(C10b, 551) and after the DE switches settle (not before the 10-07 reads), on the coordinator's message.

## Status (2026-09-30 22:05Z)
- **Declared** 21:43Z (b4dc26f9), before any draw. Amendment 1 (aae29b72): no mock / uninfected cells can be identified in the
  Sunshine file (no condition column, no viral features), so the infected state stands; asymmetric veto-only-on-HARM rule.
  Amendment 2 (b9fca152): control-only, prediction-side inputs are not reads.
- **R7a drawn** 21:44Z (85afb153): T1–T4, 300 targets each, 0 essential, disjoint from R6 and the 11 panels on disk;
  `work/draws.json` sha256 `418bb605…`. Controls are clones of 199's relabelled files. No operator, fold or score exists.
- **R7b surface:** `build_r7b.py` (construction only) is queued behind 551's efp KA and gate551 (`wait552.sh`); the machine's
  free+inactive sits near 20 GiB tonight because desktop apps hold ~12 GB in the compressor.
- **R7b instrument:** worker 553 (local Opus 5.5) generalises the real-emitter DE harness to a surface; KA on h1-test; the
  access script is written and tested on h1-test, never run on R7b. Deadline 10-06 12:00Z.

## Status update (2026-10-01 08:40Z)
- R7b surface (unfloored) built 00:27Z: 107 targets, 27,396 cells (ledger).
- **553 found R7b as built unfit for the production emitter** (bimodal control libraries, q75/q25 7.24; 107/107 placements
  fail; steering dose 0 under any class rule). **Amendment 3** (70affa19): floor 5,000 official-axis counts for controls
  and target cells, BST2 excluded, < 30 cells dropped. **Floored surface built 08:38Z** (f402a8c4): 106 targets, 18,759
  cells, 2,378 controls identical to 553's floored prediction side, 59–556 truth cells per target (median 127).
- 553's harness: all H1 known answers bit for bit; access tested end to end on h1-test; R7b-shaped rehearsal band sd
  0.0005–0.0015 overall-equivalent. Access: after the 10-07 reads, on the coordinator's message, with board-retrodiction rows.

## Candidates considered for R7b (metadata only; RESERVE.md §3 has the reasons)
| dataset | line | modality | depth (median UMI) | targets (≥100 cells) | verdict |
|---|---|---|---:|---:|---|
| SunshineHein2023 | Calu-3 (lung adenocarcinoma) | CRISPRi | 17,626 | 140 (128 carried by ≥ 2 donors; 22 essential) | **chosen** (SARS-CoV-2-infected) |
| XuCao2023 | HEK293 | CRISPRi | 2,155 | 169 | too shallow; HEK293 ≈ a donor line |
| FrangiehIzar2021 | melanoma | Cas9 KO, co-culture | 10,988 | 231 | knockout under immune co-culture; immune targets |
| Myllia / Altair A549, Jiang 5 lines | A549 etc. | CRISPR | — | — | used as banks; A549 may be an October line |
| RPE1 (Replogle) | RPE1 | CRISPRi | — | — | used by 430/516/517/521 |
| McFarland, Chang, Zhao, Gehring, Lotfollahi | various | drugs | — | — | not genetic |
| Wessels 2023 | THP-1 | Cas13 | 8,782 | 5 | too few targets |

## Requirement for R7b's claim preregistration (noted 2026-10-01 ~00:45Z, from worker 537's 10-01 read)
- **Board calibration point:** 446bj2zcp (the call-propensity block) read NULL on 10-01: −0.00079 overall vs 446bj2z; raw fid
  −0.00358, jac +0.00399, reach −0.00335. Against H1's real-emitter price the board gave jac 0.48× its predicted gain,
  fid 1.95× its predicted cost, reach −0.38×; 528's public-fold capture screen (K562 / HEK293T / HCT116) had predicted
  near-full transfer (0.87× H1's capture gain). Source: `experiments/537-arm-pipeline/RESULTS.md` §0 "10-01 READ", §4.
- **So R7b's claim must carry board-retrodiction rows** (AGENTS rule 3, "validate the validator"): inside the single access,
  before the claim's own rows are read, R7b prices the knobs the board has already read (the block: `lcap-kept` vs `cpm`; zs;
  and the hook if 496a has read by then) and reports each member's ratio to H1's price next to the board's (block: jac
  ~0.5×, fid cost ~2×). The preregistration fixes how far R7b may sit from the board's pattern before its claim rows are
  labelled "instrument not board-like" (reported, never used to veto). 553's access script already prices these MODEs.

## 553 final (2026-10-02 01:10Z): the R7b instrument is ready (`553/RESULTS.md` v6)
- Access: `access553.py r7bd5000 --owner-go <claim> …` runs the whole read in one command (~25 min, < 4 GB); its dry run
  verified the floored surface's 5 outputs by sha256, the target list and the controls against the prediction side, and
  the 10 control halves (10-02 00:54Z). The truth has not been opened.
- **The band to preregister is `band_total`** (emission seeds + split-half truths + control re-draws: the control sample
  moves an R7b-shaped reading as much as everything else together). On the h1r7 rehearsal `band_sd` 0.0004–0.0013 vs
  `band_total` 0.0006–0.0018 vs a between-draw spread 0.0004–0.0031.
- **Retrodiction must be against the SHAPE-matched rehearsal, not H1-full:** R7b's shape alone (H1 cut to R7b's cell and
  control counts) already gives the call-propensity block a board-like pattern (jac 0.45×, fid 1.4×, reach −0.05× of
  H1-full; board 0.48×, 1.95×, −0.38×). `--retro work/access553_rehearsal5` reports each member over the rehearsal mean.
- Shape-matched expectations: CP +0.0014, CPH +0.0015, K 2,000 −0.0016 overall-equivalent; the hook ~0.
- **Pricing other DE arms on a power-reduced H1:** `553/RESULTS.md`, last section (commands, cost, and why it is a regime
  instrument, not a board-shaped one: A/B/C's board has H1-full's 400 cells x 18,400 controls).
- **R7a's access machinery is ready (10-02, before C10b):** `551/r7a551.py` (plan / ka-cd4 / record / access) and
  `551/aggregate_r7a551.py` run 551's instrument unchanged on T1–T4 (its 28 imported files byte-identical to `e7be21ca`, the
  efp KA). `plan` passes (bundles == draws.json, no T output, not armed); KA-CD4 (efp) PASS; the aggregation's smoke test on
  efp's known-answer outputs runs (not a reading). The claim is drafted in `CLAIM-R7-DRAFT.md` (rules fixed; F and the DE
  arms filled after the DE reads); nothing is recorded in the ledger yet.
- **555 (cloud adversarial review, Grok 4.7 xhigh) merged and answered (d9cbc6f0, 90ec5d3d):** 3 BLOCKERs (pin missed
  `ntclevel417` — already fixed in 11fcda0b; `--owner-go` never read the ledger; unbanded pairs dropped silently by band553)
  and 9 SHOULD-FIX items, every one fixed with its own test or adopted into the claim draft (K switches report-only on R7b;
  board-like = CP raw jac > 0 and raw fid < +0.001; tails stated). `555/RESPONSE.md` has the table.


---

## R7 ACCESS, 2026-10-14/15 (science session [bf003b]; the owner session is gone; coordinator [400b6e] approved early access)

**Order of events (exact):**
1. `fill_claim552.py --bulk-config M6+P2g --block off --hook on --k default` wrote `CLAIM-R7.md`, sha256 f134a63d…, at
   2026-10-14T23:41:44Z.
2. Ledger line in `docs/CONFIRMATION-LEDGER.md` `## Reserve R7`, committed cbf443c3 at 23:42:13Z and pushed. It also states
   the switches R7 cannot read: de_uprn (630), de_dnmag (631), PIE bulk (702).
3. R7b dry run (truth not opened), then the R7b access started 23:42:3xZ (`access553.py … --owner-go CLAIM-R7.md`).
4. The coordinator's request for an independent check by [17adab] arrived after step 3. [17adab] then checked the frozen
   claim post hoc: "the claim encodes the decided recipe correctly". The not-read statement sits in the ledger line, not in the
   frozen file, which is accepted as a recorded deviation (the draft's sentence was written for identity-type sources).
5. R7a recorded 23:45:09Z (`claim.json` 4b160fb6, window 10-14T23:45Z – 10-16T12:00Z; ledger line committed 1ed7c731) and
   armed behind a wrapper that starts `r7a551.py access` only when no def_pipeline build is in its build stage
   (`work/r7a_access_1015.sh`).

**R7b reading** (`553/work/access553_r7bd5000-xBST2/access553.md`; 106 targets, BST2 excluded, seeds 0–2; 218 s, peak 3.4 GB):

| effect (overall-equivalent at R7b's anchors) | value | seed sd | band sd | effect / band sd | target bootstrap 95 % | raw nmae / fid / reach / jac |
|---|---:|---:|---:|---:|---|---|
| **hook pair `cpm+c4h` vs `cpm`** (veto-eligible) | **−0.00230** | 0.00050 | 0.00092 | −2.5 | [−0.00358, −0.00118] | +0.0051 / +0.0000 / −0.0008 / −0.0000 |
| CP `lcap-kept` vs `cpm` (board-likeness) | +0.00156 | 0.00052 | 0.00169 | — | [−0.00375, +0.00570] | −0.0018 / **+0.0077** / −0.0100 / **+0.0018** |

- **Board-likeness (the claim's rule, applied first): FAILS.** CP raw jac +0.0018 > 0 passes, but CP raw fid +0.0077 is not
  < +0.001. The board read CP with fid −0.0036 (10-01). By the claim, every R7b row is "instrument not board-like" and **none
  can veto.**
- **Verdict for the CLEAN4 hook: NO VETO** (not board-like). The 10-control-half step that feeds band553's band_total did not
  run: `control half r7bd5000c0A failed: 552's floored control bundle differs from the prediction side's controls`.
  access553's floored-surface check compares each half's controls with the full floored bundle, so it fails on every half by
  construction (a harness defect, never exercised by the H1 rehearsal). The verdict does not depend on it. The
  `--ctrl-splits` band (cells + controls, 0.00126; effect/band −1.8) is reported, and per the claim it is not the band of any
  rule.
- **What it says, report-only:** on Sunshine Calu-3 (3′, another lab, infected cells) the CLEAN4 hook moves nmae the WRONG
  way, raw +0.0051, while the board read it right (496a: −0.00158). The shape-matched rehearsal predicted ~0. With the
  instrument failing its board-likeness test, this is a caveat for D/E/F, not a decision.
- **R7a false start, 2026-10-15 02:32:57–02:34:22Z (stopped by [bf003b]).** build702 exited NOT READY at 02:32:32Z, which opened
  R7a's wrapper (its conditions are "no build702 / dn671 process alive, no def_pipeline build"). What the access touched before
  the stop, from `work/r7a/state.json`, `gate.log` and `logs/`:
  - plan checks: bundle shas vs draws.json and the instrument pin;
  - T1's CD4 banks (`T1:cd4` done 02:33:05Z: public CD4 T-cell donor banks for T1's targets);
  - the start of T1's first operator fit, `fit:held_k562_gwps` (gate OPEN at 02:33:05Z). Its log shows only the
    NON-held-out donor tables loading (orion_hct116 / orion_hek293t onto the union axis). It was killed before completing
    and is not marked done.

  No held-out truth (K562 GWPS rows of T1's panel) was loaded, and no full551 / comp551 / scoring phase ran. **No R7a truth
  was opened.** The resumed run is the single access: it skips `T1:cd4` and redoes the fit from scratch. The wrapper was
  re-armed behind build702's rebuild (pid 30635). build702.sh's one process runs r1 → ntclevel → r2, so its builder-free
  ntclevel gap is covered by "no build702.sh alive".
- **R7a:** pending, behind the wrapper; results appended here when they land.

## R7a, 10-16: a failed start (no truth opened) and an admission-only deviation, recorded before the access (2026-10-16T04:36Z, [d80a10])
- **04:33:14Z failed start.** sci1016's step 3 ran r7a551.py access. It stopped in `plan()` (access.log line 13), BEFORE
  any access: "instrument files differ from the efp-KA'd code at e7be21ca: gate551.py". The coordinator's swap fix
  (1e7c26b2) had changed gate551.py, which is inside R7a's pinned import closure. **No truth was opened** and no T step
  ran past T1's earlier state (state.json unchanged). The reserve remains unspent; the access is still single.
- **The obstacle:** this boot has no swap file (vm.swapusage total = 0.00M), so the KA'd gate rule "swap free ≥ 400 MiB"
  can never pass. R7a would wait at its gate indefinitely; its window ends 10-16 12:00Z.
- **The deviation (R7 owner's decision; coordinator's option (c)):**
  - gate551.py is restored to 1e7c26b2's blob da4f9878e6a02b89d0ee746c6c820f711cf0ad7a. Its only change vs the KA'd blob
    f36f5198 is two lines in swap_free_mib(): "total = 0.00M -> inf".
  - r7a551.py's integrity check allows exactly that file at exactly that blob (ALLOWED_ADMISSION) and records it in the
    plan as `allowed_admission_diffs`.
  - The gate decides WHEN a step runs, never what it computes. All 30 other instrument files stay byte-identical to
    e7be21ca, so the arithmetic is the efp-KA'd code. Dry-run `plan()` PASS at 04:35:45Z.
- **Re-run:** behind build705 and before readiness's run690, per the coordinator's order (`work/r7a_after705.sh`).

## R7a, 10-16 06:34Z: a second failed resume (no panel touched) and the driver pin (2026-10-16T07:14Z, [d80a10])
- **What failed:** `r7a_after705.sh` started the wrapper at 06:34:03Z, right after build705. plan() passed (the admission pin
  above), then access() stopped at its next check, `p0["instrument"]["blobs"] == claim["plan"]["instrument"]["blobs"]`.
  The record (claim.json, 10-14) holds gate551.py at its efp-KA blob f36f5198, so the pinned admission blob differs there.
  The check after it would have failed too: the record pins this driver (r7a551.py) at 4b7f640d, its blob at the arming.
  The two 10-16 changes (the admission pin and the ledger window extension) live in r7a551.py. Both checks sit before
  arm_gate() and before any panel step: **nothing was opened**, and state.json is unchanged (T1 up to full551, from 10-15).
- **My miss:** I proved the 04:36Z deviation with a plan() dry run, which does not exercise access()'s comparison with the
  record. The new test exercises the whole pre-access path.
- **The fix (same deviation, no arithmetic change):** access() now accepts gate551.py only at plan()'s pinned admission blob
  da4f9878, and only when the record holds the KA blob. It accepts r7a551.py only at its recorded blob or at the blob pinned
  by the ledger's newest "R7a DRIVER PINNED AT BLOB" line in the R7 section (a96203c2, written 07:13Z before any resume).
  Every other instrument blob, every package version and every other external script (the whitener, extract_cd4,
  aggregate_r7a551) must still equal the record.
- **Proof before the resume:** `test_access552.py` runs access() with panel() replaced by a sentinel. Without the pin line it
  FAILS at the driver check (the instrument check passed). With it, it PASSES and reaches the truth boundary (first panel
  T1, F M6+P2g); the window reads OPEN (to 18:00Z by the extension line).
- **Re-run:** the coordinator places it in the queue. R-NEW (readiness's run690) comes first.

## R7a, 10-16 09:14Z: the access resumed (after a waiter self-match), [d80a10]
- **The resume waiter (`work/r7a_after690.sh`) would never have fired.** Its trigger required `! alive after690.sh`, and its
  own file name `r7a_after690.sh` contains that string, so it matched itself. This is AGENTS 6.2's file-name trap (my
  miss). run690 ended at 09:12:08Z.
- At 09:14:32Z R7a's wrapper (`r7a_access_1015b.sh`) was started directly. build707 and readiness's chain690j both hold on
  the wrapper by name, so nothing raced it. The waiter was stopped at 09:14:51Z.
- **09:14:35Z R7a ACCESS:** every pre-access check passed, including the driver pin a96203c2 and the window extension to
  18:00Z. T1's earlier steps were skipped as done. T1:comp551 passed 551's gate on thin margins: pressure free 55 % against
  a 55 % floor, free+inactive 27.3 GiB, swap free 478 MiB against 400 (total 2 GiB).

## R7a, 10-16 10:26Z: T4's comp551 hit 062's peak-RSS guard; DEVIATION 2 and the resume ([d80a10])
- **Done:** T1–T3 complete, plus T4's CD4 banks, four fits, axes and full551. **T4:comp551 died at 10:26:42Z** on "cd4_rest
  centering: peak RSS 8.430 GiB reached 8.000 GiB". sens415:92 hard-codes 8.0 GiB for 062's enforce_memory_limit, which
  counts the process's whole peak. Memory was not short: free+inactive was ~25 GiB. R-NEW died on the same line this
  morning. T4's comp551 wrote nothing. Nothing is aggregated, and no T1–T3 price has been opened.
- **DEVIATION 2** (coordinator OK 10:33Z, guard-only, the same class as the first):
  - comp551 runs through `rss_lift552.py` (blob 884bd6e3), which replaces exactly one attribute of one module,
    `e162_062.enforce_memory_limit`, by readiness's 690j guard max(cap, 16 GiB).
  - 062/run.py lines 48–52 are a pure check-and-raise. `max_rss_gib` reaches only that function (calls at 189, 211, 224, 527,
    548, 773), so the deviation changes only whether the process raises.
  - sens415.py, 062/run.py and every other pinned file are unchanged.
  - r7a551.py pins the shim's blob. Ledger lines (10:34Z): the deviation, and "R7a DRIVER PINNED AT BLOB 423b8bfe".
- **Proof before the resume:**
  - The shim passes a 0.01 GiB cap that raises without it, and lifts exactly 1 module (`e162_062`).
  - test_access552 FAILS without the new pin (driver check) and PASSES with it (truth boundary reached, window OPEN).
- **Conditions (coordinator):**
  - The verdict is the preregistered 12-fold statistic.
  - If T4 fails again, nothing is aggregated and R7a is reported INCOMPLETE.
  - T4 runs after build707 exits.
- **Resume wrapper:** `work/r7a_access_1016c.sh`, whose name matches gate602's yield pattern. It appends "END rc=" to
  access.log for the coordinator's GO and the lab's HOLD.

## R7a VERDICT: TRANSFERS (access complete 2026-10-16T11:53:39Z, [d80a10])
The preregistered statistic (CLAIM-R7.md §2), quoted from `experiments/551-r6-c10b/work/r7a/logs/aggregate_r7a551.log`
(`results_r7a.json` sha256 142fee03…, `tables_r7a.md` 4a7a4351…):

**PRIMARY (F = M6+P2g vs M5, the post-C10a switches):** +0.01056 [+0.00768, +0.01193], 12/12 folds > 0: **TRANSFERS**;
calibrated +0.01215.

| step | adds | T1 | T2 | T3 | T4 | pooled [95 % CI] (n+) | calibrated |
|---|---|---:|---:|---:|---:|---|---|
| M5→M6 | percentile cap P 96.647 (495 / 522) | +0.00946 (3/3) | +0.00428 (3/3) | +0.00750 (3/3) | +0.01233 (3/3) | +0.00839 [+0.00560, +0.00950] (12/12) | +0.01119 |
| P2g@M6 | per-target sharpening P2g on the percentile-cap base | +0.00248 (3/3) | +0.00276 (3/3) | +0.00180 (2/3) | +0.00163 (2/3) | +0.00217 [+0.00144, +0.00300] (10/12) | +0.00096 |
| C-P@M6 | control: P2g's q permuted over targets | −0.00308 (0/3) | +0.00192 (3/3) | −0.00159 (1/3) | −0.00255 (1/3) | −0.00133 [−0.00286, +0.00005] (5/12) | −0.00162 |

- Targeting (P2g@M6 − its control, paired): +0.00349 [+0.00217, +0.00499].
- Primary retention vs R6: 1.17 [0.85, 1.32] (R6 +0.00902).
- Cumulative M0 → F: +0.02488 [+0.02206, +0.02758]; retention vs R6 1.00 [0.88, 1.11].
- **Consequence (the claim's rule): nothing changes.** The prior was ≈ +0.0069 plus P2g's C10b value. It held.
- **Magnitude, not sign:** T1–T4 were drawn by our rules, not by Arc. On such fresh-drawn panels key-class steps read ×3.08
  [1.91, 6.76] their size on an Arc-drawn panel (the lab's 10-16 panel report; shape ×1.19). Most of the primary is the
  key-class percentile cap. So R7a confirms that the decided bulk switches transfer to unseen targets with the right sign;
  it is not a board-sized gain. The A/B/C board reads of the same steps were 495 +0.00163 and 479p2g +0.00104.
- **Process:** two admission-only deviations, both recorded before use, neither touching the arithmetic:
  - gate551's swap rule;
  - 062's peak-RSS guard, lifted for comp551 only.
  Three failed starts opened nothing: plan-stage, the record comparison, and the T4 guard. The window was extended once, to
  18:00Z, by the coordinator. The access was single.

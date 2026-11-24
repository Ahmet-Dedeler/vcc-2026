# The Lab — running log (experiments 800–899)

Lead: session "VCC lab: verifier-first research (800-899)" [c61d4c], started 2026-10-15 ~22:25Z (Ahmet's local date 10-16).
Charter: `docs/LAB-2026-10-16.md`. This log is the lab's record: plan, decisions, every cloud launch and its cost, and the
checker's known-answer report. Experiment records live in `experiments/8NN-*/`. Newest notes at the bottom of each section.

## Plan (Phase A, the checker)
- **800-checker**: the contract (`experiments/800-checker/DESIGN.md`), the one command `check800.py`, the EC2 runner
  (`cloud/labcloud.sh`), statistics (paired target bootstrap, pooled folds), calibration (over-read ratio per class),
  wear ledger and best-of-N.
- **801**: fresh Arc-like panels by the lab's rule "arc800" (strength mix and coverage matched to A/B/C, measured on the
  GWCD4i pilot, a line never used as a fold); 6 WORKING panels × 3 cancer folds; 2 SEALED panels drawn and recorded only;
  the bulk fold engine (551/708's arithmetic, vendored, known answers |d| = 0) run on EC2.
- **803**: the DE engine: the production emitter on h1-test and h1-trainval (618) on EC2, candidate = emitter flags, six
  members per seed, known answers vs 618/435.
- **804**: the known-answer retrodiction: 495, 498p2g, 479sp2, 702, zs, 496a hook, DN (+ extras), hit rate and
  over-read ratio, written below.
- **805 (A.2)**: a fresh Flex DE surface from GWCD4i CD4 T cells (unworn, the challenge's chemistry).

## Dossier deadlines (coordinator, 10-16 09:32Z; plan `docs/PLAN-2026-10-15.md`, its 09:33Z note)
- **10-19 s2:** dossier at the coordinator by **10-18 18:00Z**. **10-20 s1/s2:** by **10-19 18:00Z**. 10-21: last A/B/C slot, spare
  unless R-ABC3 needs it.
- Path: dossier → the coordinator re-runs its KA → the production queue builds it (never the lab) → board read → the October recipe
  only with the lab's held-out read as well (11.1c.2). **The recipe freezes 10-20 23:59Z.**
- After the freeze, a lab candidate can still ship as a **6a replacement final entry**: one attempt, on the checker's read alone,
  **none after 11-04 12:00Z**. A candidate not ready by 10-19 18:00Z takes that route. One well-priced dossier beats two rushed.
- **Lab queue against these dates (lead, 10-16 09:33Z):**
  1. 811's one frozen reopening try with 813's TeloHAEC genome-wide rows (needs 813 stage 0 → 1a bank) and 803's DE hook for
     rule 6. On time for 10-18 18:00Z only if every gate passes first time; otherwise 6a.
  2. Grade 707's board read (10-17 s2 / 10-18 s1) against both instruments' committed predictions: 062 vs ABCp decides which
     A/B/C-panel surface the lab trusts.
  3. After the release window closes (≥ 10-22 14:00Z): **D/E/F-panel folds**: D/E/F's own 300 targets on the public lines, with
     the arc800 rules. This is the most October-like instrument the lab can have, and the 6a decision would rest on it.
  4. Before any 6a dossier: one confirmation on the sealed LZ reserve (L7, L8), recorded through the coordinator.

## Notes
- 2026-10-15 22:26Z: started. Read the charter, AGENTS, FINDINGS §0–2, PLAN-2026-10-15, 522/702/706/708/438, 534/551 (R6
  machinery), 614–618 (DE surfaces). Facts that shape the design:
  - Public lines other than H1 carry almost no DE power at Arc-like depth (Orion clean truth: median 1–2 DE calls per
    target; K562: 5; H1: 89). So DE pricing lives on H1 (and, new, GWCD4i Flex) and bulk pricing on PLOLO folds.
  - The bulk folds measure direction (relative aggregate cosine r, literal pds) and convert to board units; 702 failed on
    that conversion (point estimate, seen surface, 41 % of the change off the K562 axis). 479sp2 failed because its board
    loss was nmae through ballast, a member the bulk folds do not see.
  - A production build holds ~28 GB on the Mac (make_submission_pm_cd4): the lab runs nothing heavy locally.
- 2026-10-15 22:40Z: design written (`experiments/800-checker/DESIGN.md`); runner written; smoke test next.

- 2026-10-15 22:45Z: runner smoke test PASS (r7i.xlarge, 1.5 min, $0.005: code tarball, S3 mirror, venv, parallel jobs,
  collect, terminate verified).
- 2026-10-15 22:48Z: workers launched (Opus 5.5, high effort, detached; `experiments/800-checker/lab_subagent.sh`):
  801 (fresh arc800 panels + vendored bulk engine, $40), 803 (DE engine on h1-test / h1-trainval, $60); 22:53Z: 805
  (GWCD4i Flex DE surface, power gate first, $80). GWCD4i's cells sit in us-east-1 (free transfer to our instances).
- 2026-10-15 22:52Z: `stats800.py` KA PASS: it reproduces 522's exact arm statistics on 13 folds (52 arms, |d| ≤ 1e-13)
  and lands 522's bootstrap intervals within sampling noise. First reading off 522's own arrays (not the KA): on A/B/C's
  panel folds the bulk engine's mse side tracks 495 (fold mse −0.011 raw vs board −0.0120) but its literal pds is the
  wrong sign (fold +0.0016 vs board −0.0013): key changes move pds only through emission non-linearities the folds
  emulate imperfectly. Calibration will therefore be per member, not per arm.
- 2026-10-15 22:53Z (commit fe4dba9d): 804 PROTOCOL frozen (hit rule, surfaces, units) before any checker number; amendment 1 (22:55Z, 55d22384):
  bulk rows decided on the combined bulk + DE-side-effect price where expressible. Old instruments' record on the same
  rows: **8/13 hits; all 5 misses were false SHIPs** (479sp2, 702, K20, 458d, 471c). The checker must learn to say no.

- 2026-10-15 22:56Z: read the charter's **Amendment 1** (committed 22:28Z by [730b99]'s review, after this session first
  read the charter). Compliance actions, 22:57–23:02Z: (1) the launcher now tags instance + volume `Project=lab` and refuses
  a launch cap > $50 or a UTC day > $150; (2) all three workers stopped (their argv carried a production gate pattern via
  my launcher's system prompt; nothing of theirs was running on EC2); the launcher's prompt no longer names any gate
  script; (3) 801 resumed by session id at 22:59Z with a note on the amendment; 803 and 805 stay paused through the
  23:30–01:45Z window (≤ 2 lab agents: the lead + 801) and resume at 01:47Z by `experiments/800-checker/resume/later.sh`
  (detached); (4) commits from now add explicit files only (earlier commits added small lab-only files by directory;
  checked: nothing outside lab paths, nothing large). Timestamps in this log, DESIGN and PROTOCOL corrected to commit times.

- 2026-10-15 23:10Z: coordinator [2d6f46] clarified Amendment 1's waiter rule: watch691 / r0ready692 sleep until 10-21 22:00Z
  and a sleeping waiter does not count; the rule bites in the reserved window 10-21 20:00Z → 10-22 14:00Z (extended to "until
  R0 is submitted" if R0 slips). Outside the windows, local jobs ≤ 8 GB through gate602 are allowed (never beside def_day /
  run_arm_generic / run_ka_527 / vcc prep|submit); sub-second < 100 MB aggregations any time outside the windows. The lab
  keeps heavy work on EC2.

- 2026-10-16 00:40Z: **deviation recorded:** 801 ran its engine KA and closure discovery locally at 23:05–23:15Z (peak
  0.6–2.4 GiB, 6–16 s each, outside the 23:30 window) WITHOUT gate602, while production's DNv2 build was running. My
  mandate had allowed "light (< 2 GB) local jobs"; the charter says local jobs go through gate602 and never beside a
  production build. No harm observed (memory 79 % free), but the rule was broken by my wording. Fixed for every launch from
  now: the worker prompt allows only sub-second < 100 MB aggregations locally; anything else through gate602 or on EC2.

- 2026-10-16 00:45Z: coordinator amendment: the reserved window also opens the moment D/E/F's release is detected
  (`data/challenge-def/manifest.json`, or "RELEASE detected" in 670's guard_release.log) and runs until R0 is submitted
  ("SUBMITTED ok" in 670's submit-R0.log), besides 10-21 20:00Z → 10-22 14:00Z. Implemented as
  `experiments/800-checker/window800.sh` (prints WINDOW / CLEAR, also the daily 23:30–01:45Z); the launcher refuses a Mac
  stage > 1 GB in a window, the worker launcher refuses a second worker in a window, and the worker prompt tells every
  worker to run it before any local job.

- 2026-10-16 02:49Z: **801 first pass DONE** ($10.55 of $40): 18 fresh L-fold dumps (L1–L6 × K562 / HCT116 / HEK293T,
  held-out operators fitted on EC2 and kept on S3), the A/B/C folds (4) and efp KA folds (3) through one engine
  (`engine801.price(F, arm, base)` → per-target x, a, d, pds, pdsl, K, change energy). Every stored known answer reproduces:
  |d| = 0.0 on the Mac (96 rows), ≤ 1.6e-15 on dumps rebuilt on EC2 (66 rows). Findings worth keeping:
  - **dump708.py (production's 708 day-of tool) crashes after writing its npz** (json `default=float` on a Path in
    `spec["fold_ops"]`): reported to the coordinator with the one-token fix; not touched by the lab.
  - Operator fits are not bit-reproducible across machines / thread counts (three sha256 for one efp fit); prices move ≤ 1.4e-11.
  - **Panels matched to A/B/C on coverage, CD4 strength mix and basal expression still have much heavier percentile-cap
    tails** (fold cap 7.8–19.9 vs A/B/C 5.2–7.0), so key-class gains read 3–10× larger on L folds than on A/B/C's own folds.
    What makes A/B/C's key tail light is not the strength mix; that is an open question with a direct October consequence
    (D/E/F is Arc-selected like A/B/C). Key-class prices on BULK-W must go through the over-read ratio.
- 2026-10-16 02:55Z: ka804 launched on EC2 (r7i.8xlarge, 22 folds, cap $5.68); 801 resumed for the 710 compositions.

- 2026-10-16 03:05Z: **805 stage 0, CD4 Rest: preregistered verdict FAIL (blob), power PASS.** 100 random non-essential
  knockdown-significant targets, 400 cells from all 4 donors, 18,400 pooled controls: n_real median 34.5, 94 % ≥ 10, 79 % ≥ 20
  (bar 20 / 60 %); knockdown strong (own-gene remaining median 0.11). But 615's split-half blob share is 0.120 > 0.10 (h1-test
  0.049). Depth is 7.5k UMI on the official axis (0.38× Arc's 20k; no thinning possible). A prereg amendment (before any gate
  number) requires every donor to fill its quota: targets carried by some donors only read n_real 121 vs 35, because donor
  genes are called DE against a single pooled control. Stim48hr runs next by the prereg's rule; ~$1.3 so far.

- 2026-10-16 03:25Z: disk fell to 61 GiB (production builds). The runner's `collect` had pulled 805's S3-side intermediates
  (9.1 GB) and 801's dump copies (1.7 GB) to the Mac; both deleted after checking S3 holds identical objects (sizes equal;
  801's own work/dumps keep one local copy). `collect` now fetches only STATUS, worker.log, jobs/ and declared repo/ paths.
  Disk back to 70 GiB. Lab local footprint ~8 GB (cap 30).

- 2026-10-16 03:30Z: **805 stops at stage 0 (preregistered FAIL on the blob bar in both conditions; power passes in both).**
  Stim48hr: n_real median 36, 88 % ≥ 10; blob split-half 0.1227 (bar 0.10; h1-test 0.049); depth 10.9k UMI (0.55× Arc).
  Rest: 34.5 / 94 %, blob 0.1201, 7.6k UMI. So GWCD4i carries Arc-like DE power at 400 cells but a shared across-target
  direction ~2.5× h1-test's; it is not a clean DE truth by the rule written before the data. A descriptive diagnostic (blob
  vs depth / long-gene axes) is running; the verdict stands. DE pricing stays on DE-H1 (h1-test, h1-trainval).

- 2026-10-16 03:35Z: coordinator: disk dipped to 53 GiB at 03:30Z; lab asked to stay ≤ 15 GB locally (D/E/F + R0 need ~40 GB).
  A second 16 GB of 805 intermediates had come down through a wait still running the old frozen runner; deleted after the
  same S3 check. Lab now 8.1 GB apparent (~4.5 GB real: 801's bundles are APFS clones). Rule from now: fold dumps and
  intermediates live on S3; pull only what a local aggregation needs, delete it after, log it. Disk 80 GiB free.

- 2026-10-16 03:45Z: 805 final ($1.90; committed bed1d14d). Blob diagnostic: not Orion's long-gene axis (cos +0.05 / +0.015
  vs h1-test +0.037); partly a library-depth axis (16–21 % of the shared direction's energy vs h1-test 1.4 %; |cos(PC1,
  depth)| 0.57–0.62 vs 0.06). With depth projected out the split-half share is still 0.108 / 0.108 (> 0.10). A future attempt
  needs a new prereg: all-donor quota, a depth-matched control pool, Stim48hr (0.55× Arc depth). Power is not the obstacle.

- 2026-10-16 04:00Z: 801's follow-up DONE ($0.83): the 710 compositions on all 22 folds (L folds h2c1 / h2.5c1 bit-identical
  KA; A/B/C folds h2c1.5 = 518's fST to 3.4e-14). Local disk pass: 801's 18 L dumps and 22 design files deleted after their
  S3 mirror copies were size-checked (the 4 efp KA dumps, not mirrored, stay); collected repo/ copies of finished runs pruned
  (2.2 GB; they remain at output/lab/<run>/repo/ on S3; 803's b1 copy restored because 803 is live). Lab local 6.1 GB
  apparent; disk 108 GiB free. 710 is being priced on EC2 (ka804-710).

## Known-answer report, part 1: the bulk rows (2026-10-16 03:20Z; `experiments/804-known-answers/`)
Checker = 801's engine on 18 fresh L folds (BULK-W) and A/B/C's own 3 cancer folds (BULK-ABC), frozen units (literal pds +
mse through 518's identity), paired target bootstrap 4,000. Arrays: run 20261016T030032Z-ka804-b ($0.30); grading
20261016T031052Z-grade804 ($0.01). Board deltas from docs/BOARD-LEDGER.md (`board_pairs804.py`).

| row | board | BULK-W checker [95 %] | BULK-ABC checker [95 %] | grade (W / ABC) | old instrument |
|---|---:|---|---|---|---|
| 495 cap 3 → 6 | +0.00163 | +0.00776 [+0.0067, +0.0085] SHIP, 18/18 | +0.00239 [+0.0009, +0.0039] SHIP, 3/3 | HIT / HIT | HIT |
| 479p2g | +0.00104 | +0.00335 [+0.0027, +0.0041] SHIP | +0.00364 [+0.0024, +0.0050] SHIP | HIT / HIT | HIT |
| 498p2g | +0.00101 | +0.00299 [+0.0024, +0.0037] SHIP | +0.00347 [+0.0022, +0.0048] SHIP | HIT / HIT | HIT |
| 479sp2 | −0.00071 | +0.00379 [+0.0030, +0.0045] SHIP | +0.00363 [+0.0024, +0.0050] SHIP | **MISS / MISS** | MISS |
| 702 PIE blend (K562 folds) | −0.00255 | −0.00048 [−0.0023, +0.0015] NULL, 1/6 | −0.00236 [−0.0054, +0.0006] HARM (1 fold) | **HIT / HIT (HARM)** | MISS |

- **Bulk hit rate 4/5 on both surface sets (old instruments 3/5 on the same rows).** The new catch is 702: run through
  production's own fold arithmetic (composition → key → sharpening → render → ballast), PIE's K562 rows LOSE direction
  (r −1.96 % on A/B/C's K562 fold, −0.39 % on the fresh K562 folds), where the old read (PIE vs an Orion donor mean, ew-cosine)
  said +3.7 %. On A/B/C's K562 fold the checker's −0.00236 sits on the board's −0.00255 (ρ 0.93).
- **The miss is 479sp2**, and it is a class failure: the folds price strength-tercile sharpening like key-tercile sharpening
  (+0.0036 … +0.0038 both), while the board split them (+0.0010 vs −0.0007; S_P2's loss was pds −0.0011 raw and nmae +0.0018
  raw through ballast). Neither literal nor board-emulated (λ 0.5) fold pds sees it (both +0.007 … +0.008 raw vs board
  −0.0011). Per-target sharpening that re-assigns q by a strength statistic is NOT priceable on the bulk folds; it needs the
  emission (DE-H1 with the six members from emitted cells, requested from 803).
- **Over-read ratios (checker / board, overall):** BULK-ABC: key 1.46, shape 3.44–3.51, external 0.93. BULK-W: key 4.75,
  shape 2.96–3.23, external 0.19. The fresh panels' heavy key tails (801) over-read key changes 4.75×. By member on BULK-ABC:
  mse tracks the board (ρ_mse 0.89 key, 1.37 / 1.45 shape, 1.65 external); literal pds is over-read ~5× for sharpening and
  wrong-signed for the key change (board-emulated pds λ 0.5 gets the key change's sign right: −0.00098 vs board −0.00131;
  a diagnostic found after the read, not yet a validated unit).
- **Calibration file:** `experiments/800-checker/ka800.json` (`calib804.py`). Use: price on BULK-ABC (the board's panel) for
  the board read and on BULK-W for D/E/F transfer, divide by the class ρ, and treat sharpening-by-statistic changes as
  unpriced on bulk folds.
- **Prospective (committed 1b1f4a3c, 03:15Z, before any of these reads):** 705 rider → CLEARED, board ≈ +0.0001; 707
  strong-tier swap → NULL (≈ +0.0002: mse better, pds worse on A/B/C's panel). **710** (H1 ×2.5, committed c5185d56,
  04:03Z) → NULL on A/B/C (−0.00005 [−0.0005, +0.0004], P(≤ 0) 0.57: mse slightly better, pds slightly worse); fresh folds
  WEAK +0.0003. The checker does not support raising H1 above ×2 on A/B/C-like panels.

## Known-answer report, part 2: the DE rows and the verdict (2026-10-16 07:10Z)
DE-H1 = 803's engine (the production emitter on h1-test 85 targets and h1-trainval 126, 3 seeds; per-target members from the
emitted cells; KA: its cell-eval2 cross-check on one run matches per target to 4e-15, and every 618 / 619 / 612 pair delta
reproduces to ≤ 1.4e-8 raw, UPRN to 2e-4 = re-roll). Frozen protocol: fid + jac + nmae at f 0.35, paired bootstrap 4,000,
both surfaces equal weight. Graded on EC2 (20261016T070036Z-gradede804, $0.01).

| row | board | checker f 0.35 [95 %] | grade | ρ | f 0.61 (post-hoc call-set rule) | old instrument |
|---|---:|---|---|---:|---|---|
| zs | +0.00055 | +0.00051 [+0.00035, +0.00069] SHIP | HIT | 0.93 | +0.00056 | HIT |
| CLEAN4 hook (496a) | +0.00056 | +0.00053 [+0.00029, +0.00081] SHIP (h1tv board form +0.00040) | HIT | 0.95 (0.72) | +0.00042 | HIT |
| DN d_m200 | +0.00037 | +0.00109 [+0.00065, +0.00157] SHIP | HIT | 2.9 | +0.00053 | HIT |
| UPRN (630arm) | +0.00093 | +0.00279 [+0.00197, +0.00366] SHIP | HIT | 3.0 | +0.00281 | HIT |
| idall3 vs plain hook | −0.00050 | +0.00029 [+0.00003, +0.00056] SHIP | **MISS** | — | +0.00036 | MISS |
| K 2,000 | −0.00072 | +0.00071 [+0.00033, +0.00106] SHIP | **MISS** | — | **−0.00059 HARM (ρ 0.82)** | MISS |
| K 1,600 | −0.00211 | +0.00234 [+0.00118, +0.00333] SHIP | **MISS** | — | **−0.00188 HARM (ρ 0.89)** | HIT |

### The verdict (PROTOCOL.md §6)
- **Core rows at the frozen settings: all five positives ship (495, 498p2g, zs, hook, DN) and 702 does not; 479sp2 ships.**
  The pass criterion ("does not ship 702 or 479sp2, ships ≥ 4 of 5") therefore FAILS on one row, and by the protocol the
  class it cannot price is named: **per-target sharpening re-assigned by a statistic** (479sp2), priced only with a range
  spanning zero until DE-H1 can carry bulk arms (803 follow-up).
- **All graded rows, headline surfaces, frozen settings: 9 / 12 hits** (bulk 4/5 on BULK-W, DE 5/7 on DE-H1). The old
  instruments on the same 12 rows: **8 / 12**, and every one of their misses was a false SHIP (479sp2, 702, K20, idall3).
  The lab's new catch is 702 (NULL on fresh K562 folds, HARM on A/B/C's K562 fold at ρ 0.93); it loses K16, which I-662 had.
- **Call-set size is a regime question, as 662 / 619's amendment said:** at f 0.61 both K rows are HARM at board size
  (ρ 0.82 / 0.89); at f 0.35 both are mis-signed. Since f 0.61 was chosen after the f 0.35 read, it enters as a CLASS RULE
  for future call-set candidates (read at f 0.61, with this calibration), not as a validated hit. With it: 11 / 12 minus
  idall3 and 479sp2 = 10 / 12.
- **idall3 stays unexplained**: every instrument (h1-test, 618, I-662, now this engine) puts the identity hook ahead of the plain
  hook; the board put it 0.0005 behind. Identity-type sign priors are a class this checker over-reads; price them with a
  range spanning zero.
- **Calibration (ka800.json):** DE magnitude passes ρ 0.93–3.0 (median 2.9: add-based passes over-read ~3×, zs does not);
  DE sign hook ρ 0.95; call-set at f 0.61 ρ 0.86. Bulk: see part 1.
- **What this means for new candidates:** a candidate gets a SHIP from the checker only on a lower bound > 0 after
  best-of-N, on the surface set that carries its class; its board projection divides by its class ρ; classes with a sign
  failure (sharpening-by-statistic, identity-type sign priors) are reported with a range including zero. The prospective
  rows (705 / 707 / 710) are the next, out-of-sample test.

- 2026-10-16 07:20Z: coordinator accepted Phase A; Reserve LZ appended to docs/CONFIRMATION-LEDGER.md. Two conditions for
  Phase B, applied: the f 0.61 call-set rule is labelled post hoc and non-decisive until a later board read confirms it
  (PHASEB.md, ka800.json `post_hoc`); and why pds tracks poorly while mse tracks well is written down before any lane leans
  on pds (PHASEB.md: scale-invariance + the ~300-way ranking; key changes are invisible to pds by construction; direction
  changes move pds by the neighbours they cross, and fold truths are spaced unlike the board's).

- 2026-10-16 07:30Z: for 708 (science [d80a10], cc coordinator): class ratios in 522's units (0.0710·r + 0.1245·Δpds_lit,
  708's score708 unit) on BULK-ABC: key 1.25 (fresh heavy-tail panels 4.48), sharpening 1.63–1.78, external rows 0.68
  (fresh 0.14). 708's arms: r (rider) → key, divide by 1.25 unless D/E/F's fold cap hi ≫ 7; s (707 swap) → direction,
  provisional 0.68 (702's operation), uncalibrated until 707 reads; A:<design> → composition, no board point, range 0.68–1.25
  until 710 reads. 707 reads +0.00107 in 522 units vs +0.00044 in lab units on BULK-ABC (negative pds, positive mse), so its
  board read also tests 522's 0.34 pds discount. Informational only.

- 2026-10-16 07:55Z: **803 follow-up (sharpening on the real H1 emitter, $1.21): DE-H1 does NOT reproduce 479sp2 either.**
  P2g − uniform has the board's signs (pds +0.0027…+0.0033, mse −0.0023…−0.0051 raw, nmae ~0). S_P2 − P2g: pds +0.0043 /
  −0.0001 (board −0.0026), nmae +0.00014 / +0.00017 (board +0.0017, a tenth), mse disagrees between surfaces (the harness has
  no re-level). K moves only +0.45 / −0.20. 803's reading: the gap is the key and panel composition (21 H1 targets sit at the
  key floor and fill P2g's low tier), not emission. **Verdict on the instrument: sharpening re-assigned by a statistic is not
  priceable on any current lab surface** (bulk or emission); proposals in that class get a range spanning zero and need a
  board read. Untested next step (not scheduled): S_P2 − P2g on DE-H1 with a key built to A/B/C's floor and cap shares.

- 2026-10-16 08:05Z: **caveat on BULK-ABC's calibration (712, 9bf8b6c9):** A/B/C's PLOLO fold operators are 062's cloud whiteners
  (basis = every donor row, K562 essentials in, no per-donor cap), while L folds, real/efp/efr, production e2f65ecf and the 10-22
  day path use `--basis-prune essential` + the 12k cap. So "A/B/C light key tail" and BULK-ABC's key ρ (1.46 lab / 1.25 in 522
  units) may be partly an operator effect; the raw delta-row-norm tail is the same on every panel (712). 801 is refitting A/B/C's
  held-out operators with the day recipe (ABCp folds) and re-running the KA arms; until then the key ratios are flagged as
  possibly instrument-specific.

- 2026-10-16 08:12Z: **Phase B, first proposals.**
  - 810 (gated un-whitening of the donor base on the strong tier): judge **REVISE** (VERDICT-1.md): add the closed global
    construction as a same-run known-answer comparator, fix mechanism numbers, honest tries, a decision rule, provenance.
  - 812 (DE emission): **NO-GO on all four hypotheses, no compute asked** (NTC-guide control bias +0.0002, real but 0–5
    genes per context; sign-confidence block: jac pays for fid; cis neighbours already in 425; head size 400 optimal). Its
    reading: the field's DE gap is ≈ +0.04–0.05 of sign accuracy on the truth's responding genes (per-target model skill),
    i.e. 811's territory; the emission side is at its optimum for the sign information we have. Process note: its four
    screens ran (~$0.40 EC2, 1–13 s gated local jobs) before a proposal existed, reported as exploratory with try counts.
  - 811 (many-line consensus, per-line centring + reliability weights, nested selection L1–L3 → confirm L4–L6 + BULK-ABC;
    8 configs): sent to the judge. Its baseline: more public lines carry more transferable direction on fresh folds (17/18),
    but against the incumbent it is a W-only effect (A/B/C NULL at λ 0.10, HARM at 0.50: mse up, literal pds pays), with the
    extra lines bringing 2–4× the truth's across-target shared energy.

- 2026-10-16 08:28Z: **810 PROPOSAL-2 → GO** (judge, with binding conditions G1–G4): one priced arm (sgu1: base→pre on the top
  z_src tertile, β 1, H1/CD4 kept, row norms kept) charged N = 6, against inc, 707 and the CLOSED global form sgu1all as a
  same-run known answer (STOP if the checker reads sgu1all clearly positive: then it cannot see the γ knife-edge). G1: no
  dossier on 062-operator BULK-ABC alone (re-price on ABCp first). G2: HOLD-707 projections divide by max(ρ, 1). G3: missing KA
  key = failed KA. G4: N = 6 in candidate.json, 063 on the closed list. Budget ≤ $5. 811 → REVISE (drop agreement weights,
  one bank list on every surface, H1 out on DE-H1, decision rule). 812 → building skill812 (DE sign-skill scorer).

- 2026-10-16 08:42Z: **811 PROPOSAL-2 → GO** (judge, G1–G5): many-line consensus (the D/E/F day's bank list, no KOLF2), per-line
  centring × λ {0.10, 0.25} = 4 arms; nested (select on L1–L3, confirm once on L4–L6 + BULK-ABC); decision set ABCp; falsifier F1
  restated on the many-line increment; DE-ineligible arms (|dK| > 5 %) cannot be selected; ≤ $3. 811 itself predicts NULL/WEAK on
  A/B/C (no dossier expected): the run prices the consensus family for October either way.

- 2026-10-16 08:35Z: **812: skill812 built** (`experiments/812-de-signs/skill812.py`): a DE sign-skill scorer for any H1-free
  source on h1-test + h1-trainval, per class (a0, a1, head a1, u1, by expression bin) vs CLEAN4 and vs the incumbent's
  delivered sign, paired bootstrap, priced at H1's per-class rates. KA: CLEAN4 a0 0.5171 / 0.5225 (700's signsplit row);
  the production emitter's realized signs reproduce rates812 to 5e-5. Information runs: the bulk delta's own sign beats
  CLEAN4 nowhere usable; the 12-line public mean beats it only on u1 (+0.0013 [+0.00007, +0.0025] overall-eq at H1 rates),
  which does not survive best-of-3, and loses on a0 (−0.0025). 812's own prior: NO-GO (516 closed the same object: an H1-only
  gain). Its PROPOSAL-u1hook (a cheap held-out stage 0 before any DE-H1 price) is with the judge.

- 2026-10-16 08:45Z: lead check of 811's PROPOSAL-2-ADDENDUM (sha256 2a3c8e95…, m811.py e00192e6…): G1–G5 covered as the judge
  wrote them (ABCp detected by fold id with inc = pct+P2g there and only the engine's own KAs; F1 restated against src_lX; DE
  eligibility before selection; N 9 / 10 branch; the align811 correction; wear rows to work/wear811.jsonl). Accepted.

- 2026-10-16 08:55Z (corrected from 09:10Z; commit time): **810 priced (GO run, $1.23 lane total): no dossier.** sgu1 (gated base un-whitening): BULK-ABC −0.00200
  [−0.0053, +0.0011] (pds eats 158 % of the mse gain), BULK-W +0.0089 [+0.0054, +0.0130] 18/18; it loses to 707 paired on ABC
  (−0.0024). decide810 → DE-H1 FIRST (K moves > 5 % on 12 W folds); past that it would CLOSE (C1, C4). Lead decision
  (`810/LEAD-DECISION-1.md`): DE-H1 not funded; family CLOSED for the A/B/C-panel route; one ≤ $0.5 information read on ABCp
  (day-recipe operators) because the candidate undoes the whitener and BULK-ABC uses 062's. **Checker validation, new known
  answer:** the closed global un-whitening (sgu1all) reads −0.0099 on BULK-ABC with pds eating 232 % of mse, as 212/097 found:
  the checker sees the γ knife-edge. A KA tolerance ruling (relative 1e-12 on unnormalised sums; the failure was 3.6e-15
  relative) is in `810/RULING-KA.md`. check800 now prices BULK-ABC, BULK-ABCp and BULK-W separately (the pooled number is
  labelled "NOT a decision number"); 810's ledger row (+0.0071 pooled) was produced before that fix.

- 2026-10-16 08:59Z (corrected from 09:20Z; commit time): **811 priced (GO run, lane $0.92): no dossier, a lead.** Many-line consensus (the D/E/F day's bank list),
  λ 0.25, centring off (c0_l25): W-tune +0.0085 [+0.0062, +0.0111] 9/9 (N 10); held-out W-conf **+0.0068 [+0.0042, +0.0096] 9/9,
  +0.0023 after e(7)**; 062-BULK-ABC +0.0009 [−0.0016, +0.0036] 2/3 (reported only), with the many-line increment over sources-only
  +0.0031 [+0.0012, +0.0049] 3/3. F1 does not fire: centring removes the shared energy (0.022 → 0.0005) yet keeps 89 % of the
  increment on fresh folds, so it is target-specific direction (class transfer). No dossier because ballast K moves +5.86 % on one
  W-conf fold and DE-H1 cannot yet price external rows (rule 6). Pending: the ABCp read (801), DE-H1's new bulk_transform hook
  (803, now commissioned), λ 0.10 on ABCp, and the same model with TeloHAEC / KOLF2 full rows (813). 811's 56 wear rows merged.

- 2026-10-16 09:04Z (corrected from 09:06Z; commit time): **812 PROPOSAL-u1hook → KILL** (judge): the closed 516 construction (same-target public rows signing
  uncalled genes through 496's hook; 516 priced the rider at ≤ +0.001 and its plain-mean form lost on the one valid held-out line,
  K562 −0.0023). Lane 812 is PARKED with two products: four closed DE-emission hypotheses (NTC control bias, sign-confidence block,
  cis front, head size) and the skill812 instrument, which 811's H1-free predictions use as their DE objective. 813's TeloHAEC
  proposal is with the judge.

- 2026-10-16 09:16Z: **813 PROPOSAL (TeloHAEC genome-wide CRISPRi as a fourth fold line) → REVISE** (judge,
  `813/VERDICT-1.md`): "the right next data move for the lab", a $3 kill gate. Fix list R1–R7: KA-SG as written fails an
  unbiased estimator (noise-attenuated OLS; needs a reliability-corrected or TLS slope, bar set on KA-SIM first); P2 cannot
  fail (702's PIE rows are not buildable on TeloHAEC, the other rows are near-certain positives); stage 1 had no operator for
  the L panels (fresh 3-donor day-recipe fits); P3 donor-count confound; decision role and leakage unstated; 105 uncited;
  data terms not a STOP gate. Lead calls sent with the resume (`800-checker/resume/813-verdict1.md`): TeloHAEC folds are
  reported apart, never pooled into BULK-W and never a veto alone (all six challenge lines are cancer lines); a cancer-fold
  gain with a TeloHAEC CI upper bound < 0 carries a written flag; P4 / the 811 bank use is deleted here and stays 811's own
  pre-counted try; registry changes come as an UNAPPLIED patch. 813 resumed to write PROPOSAL-2.md; no download before a GO.

- 2026-10-16 09:19Z: **801 ABCp registered** (A/B/C folds with DAY-recipe held-out operators: `build_union_whitener.py
  --basis-prune essential`, 12k cap; operator sha256 k562 `48cb88b6…`, hct116 `a0e27947…`, hek293t `5d982462…`; 062 KA PASS ×3,
  worst |d| 8.3e-14; five runs, $1.32). Answer to 712's question: the A/B/C-vs-L percentile-cap gap is about **one third operator
  recipe, two thirds panel** (cap hi 6.15/6.52/7.03 → 6.18/7.84/8.39 vs L medians 10.4/14.8/11.6). **KAs regraded on ABCp**
  (`804/grade804p.py`, 2.6 s, 40 MB): the same grades as 062 on all six board-pointed rows (5 HIT incl. 702 HARM −0.0020 vs board
  −0.0025; 479sp2 still MISS); key ρ 1.44 (062 1.46), shape ρ 3.05, external 0.79. **Prospective, committed before the board
  reads (7a38fb29):** 705 NULL on both instruments; **707 discriminates them**: ABCp +0.00178 [+0.0005, +0.0032] SHIP 3/3
  (r +3.17 %) → board +0.0006 … +0.0023, vs 062's NULL (+0.0001 … +0.0005). A board read above +0.0005 favours ABCp. 810 (info
  read) and 811 (decision read, N = 10 branch) relaunched on ABCp at 09:18Z.

- 2026-10-16 09:23Z: **Coordinator HOLD on lab local jobs** (message 09:21Z): none, gate602 or not, until R7a's access wrapper writes
  "END rc=" (expected 10:45–11:00Z; R7a's gate opened on thin margins). Implemented in `800-checker/window800.sh` (rule 4):
  while `800-checker/HOLD` exists and its log (`551-r6-c10b/work/r7a/access.log`, the only 551/552 file the lab reads) lacks the
  pattern, it prints "WINDOW hold (local jobs only)" but exits 0, so the agent count (windowmgr, lab_subagent) and EC2 are
  unaffected. It lapses by itself when the pattern appears. The four new sessions the coordinator saw are the lab's
  headless workers (judge-813, 813, 810, 811). 813 froze PROPOSAL-2.md; it is with the judge (`819/813-public-data-20261016T092345Z`).
  EC2: `ka804-710p` (row 710 on ABCp, cap $0.9) running.

- 2026-10-16 09:28Z: **710 on ABCp: KA FAIL, no number read** (run `ka804-710p`, $0.09). 062's A/B/C design `fin_h2c1.5` differs
  from the ABCp fold's fin by 1.20 / 3.75 / 4.03 (the compositions live in the operator's whitened space), so the frozen KA
  (|d| < 1e-12) stopped all three folds, as it should. 801 resumed (follow-up 3, ≤ $2, by 10-17 20:00Z) to re-cut the H1 ×2 /
  ×2.5 compositions on ABCp with the day operators; `price710p.py` then prefers `ABCp__<held>.npz`. check800 now projects every
  block with its own surface's ρ (0468bef4); ABCp's ρ entries are written by calib804 after the hold lapses.

- 2026-10-16 09:29Z: **813 PROPOSAL-2 → GO, stage 0 only (≤ $5)** (judge, `813/VERDICT-2.md`): every VERDICT-1 item fixed
  (KA-SG now a split-half unattenuated ratio with simulated bounds; 702 dropped from P2 with the PIE training-set citation and
  479sp2 reported only; bins renamed, no "TRUSTED"; TeloHAEC folds reported apart, never a veto; data terms PASS, sha
  `fc17237b…`). Binding G1–G6: 0b only after the lead commits freeze813.json's sha; P2 gets a third bin UNRESOLVED; operator-fit
  guard trips are infra reruns; explicit L-fold dump wiring; KA-SG fallback stated; bank never on the Mac. Lane cap stage 0 ≤ $5,
  stage 1 ≤ $6 (released by the lead after stage-0 RESULTS). 813 resumed for stage 0a.

- 2026-10-16 09:31Z: **811 ABCp decision read → WEAK, no dossier** (N = 10 branch; $0.22, lane $1.14 of $3). c0_l25 on BULK-ABCp
  (vs the October incumbent pct+P2g) **+0.00264 [−0.00003, +0.00561] 3/3**, after e(10) +0.00042 [−0.00224, +0.00339]: §5 rule 1
  fails on its own (pdsl +0.00331, bound after e(10) −0.00114), and rule 6 (DE-H1) still blocks independently. ABCp reads ≈ +0.0017
  above 062's read (+0.00092), as 707 does (+0.00178 vs +0.00044): both operator and incumbent differ, so this run does not
  separate them. Family NOT closed (point > 0); reopening frozen to one counted try with the 813 lines. λ 0.10 not read (it
  would be an A/B/C-chosen dose; 811 advises against). The lead keeps the 811 construction as the lab's strongest transfer
  lead on fresh lines (W-conf +0.0068, +0.0023 after e(7)); 803's hook answers rule 6.

- 2026-10-16 09:36Z: **810 ABCp information read** ($0.09; price810 / decide810 shas unchanged; 42 KA keys ≤ 2.8e-13, incl. an
  independent reproduction of 801's ka-abcp 707 and inc arrays). On BULK-ABCp (vs cap6+P2g(abc)): **707 +0.00178 [+0.0005, +0.0032]
  3/3**, sgu1 **+0.00012** [−0.0032, +0.0035] (062: −0.0020; pds share 1.58 → 0.97), sgu1all −0.0049. decide810 as frozen: DE-H1
  FIRST, then CLOSE on C4 (sgu1 − 707 −0.0017). **Answer: about one fifth of the W/ABC split (0.0109 → 0.0088) was 062's whitener;
  the rest is the panel.** The day operator lifts 707 as much as sgu1, so the order on A/B/C's panel is unchanged. No dossier; the
  family stays closed for the A/B/C route. Wear: sgu1 N 7, 707 N 2 on BULK-ABC-type reads.
- **Open question this raises (lead, for the 6a route):** W folds read direction changes strongly positive (sgu1 +0.0089 18/18,
  811 c0_l25 +0.0068 9/9) while A/B/C's panel reads them near zero, on the SAME truth lines and now the same operator recipe. Is
  A/B/C's panel inside the spread of the six arc800 panel draws (then W, averaging panels, is the better predictor of a new draw
  like D/E/F), or systematically outside it (then arc800's draw rule misses how Arc builds panels, and A/B/C-panel folds are the
  better proxy)? The per-fold records already on disk answer it per held line (rank of ABCp/<held> among L1–L6/<held>); the lead
  runs that aggregation once the coordinator's hold lapses.

- 2026-10-16 10:00Z: **813 stage 0a: KA-SIM FAIL as run → lead ruling** (`813/RULING-0a.md`). Run `s0a813`, $0.094, metadata
  only (180.6 MB of 66.3 GB; X never opened). Estimator A chosen (median cos 0.765 vs B 0.745). KA-SIM (ii) 0.200 > 0.10 on the
  injected guide-load design. 813 found two code departures from its frozen design (guide-count strata capped at 12; a twin at
  4.7× fewer incidences) and that KA-SG's statistics are log-bias dominated on single-guide rows. Ruling: ONE corrected 0a
  (freeze813-v2, ≤ $0.15, diff limited to the preregistered strata and twin, bars unchanged, a v2 failure is final); KA-SG →
  G5's "uninformative" fallback, no new statistic; v1 kept as a negative: **coarsening guide count leaks a load-dependent offset
  (0.20 of a median effect)**, which binds later estimators here. 710 on ABCp: KA PASS on 801's re-cut designs (run
  ka804-710p-b, $0.10); its grade waits for the hold.

- 2026-10-16 10:33Z: **813 corrected 0a (v2) FAILS (ii) → TeloHAEC is NOT a lab fold line under the frozen bars** (run
  `s0a813v2`, $0.129; lane total $0.223; no expression byte fetched; freeze813-v2 sha `ed11ef33…`). v2 as preregistered: (i) A vs
  matched-incidence twin 0.971 PASS; (ii) offset 0.084 PASS uninjected, **0.143 FAIL** under an injected guide-load effect (0.5 ×
  median effect per guide); (iii) PASS; estimator A. Negatives recorded (binding on this lane): coarsened k leaks 0.200; uncapped
  stratum matching still leaks 0.143 (likely sparse high-k strata where carriers lack matched non-carriers; untested); first-order
  deconvolution (B) fixes neither. A re-proposal would need a high-k rule fixed in advance, with v2's simulation as its KA.
  **Lead decision: lane 813 PARKED.** DLD-1 (#2) is a cancer line but thin on every panel (65–82 targets), multiply perturbed (the
  same estimator problem) and 116's split-half cosine 0.070 warns of marginal power; KOLF2.1J cells (#3) are pluripotent and a
  ~190 GB pull. Consequence for 811: its frozen reopening (with 813's lines) has no TeloHAEC bank; KOLF2 full rows would be the
  only route, and the lead decides that after the panel-rank analysis says which instrument the 6a route rests on.

- 2026-10-16 10:41Z: **Instrument finding: A/B/C's panel is NOT exchangeable with the arc800 panel draws** (`800/panelrank800.py`
  on EC2, run `agg-abcp-1016b`, $0.005; table `800-checker/work/panelrank800.txt`). For each held line, ABCp/<held> ranked among
  L1–L6/<held> on the same arms: **bottom of the seven in 14/21 (arm, held) cells (3.0 expected), mean percentile 0.19**. By class:
  **key** (495 cap6−cap3, 705) bottom in 6/6 (z −1.8 … −5.5); **direction** (707, sgu1, 811 c0_l25) bottom in 8/9; **shape**
  (479p2g, 479sp2) inside the spread (ranks 2–6). ABC062 shows the same pattern. This explains the W folds' key over-read (ρ 4.75 vs
  1.46 on A/B/C) and the recurring W-wins / A/B/C-null split of direction candidates (522/706/707/810/811): the fresh panels
  reward key and direction changes more than Arc's panel does. **Confound:** ABC-panel folds use A/B/C's fold rules (CD4 ×1.5, its
  controls and banks), L folds D/E/F rules. 801 follow-up 4 builds LABC (A/B/C's 300 targets on the L-fold recipe) and a per-fold
  truth-strength profile to separate panel from rules (≤ $6, by 10-17 12:00Z). **Until it reports, the lab treats BULK-W reads of
  key or direction changes as optimistic for any panel Arc draws (D/E/F included).** The calibration now holds ABCp's own ρ
  (key 1.44, shape 3.05, external 0.79). 710 on ABCp: −0.00004 [−0.0005, +0.0004] NULL, committed as a prospective row before its
  10-18 board read.

- 2026-10-16 10:46Z: **Preliminary per-class over-read factors** (`800/overread800.py`, run `overread-1016`, $0.011; two-level
  bootstrap over the six panel draws and over targets, 4,000 draws): fresh panels / A/B/C-panel folds, pooled over arms × held
  lines. **vs ABCp:** key **3.78 [2.15, 10.9]** (495, 705; 6 cells), direction **4.72 [1.88, 27.5]** (707, sgu1, 811 c0_l25; 9
  cells), shape **1.16 [0.80, 1.74]** (479p2g, 479sp2). vs ABC062: key 3.39 [2.00, 8.27], direction unbounded (A/B/C sum CI spans
  0), shape 0.98 [0.69, 1.39]. Consistent with the board-calibrated ratios (key ρ_W/ρ_ABC 4.75/1.46 = 3.25; shape 3.10/3.48 =
  0.89). Coordinator's four asks for the 10-20 freeze (panel vs rules with a number; the day-of instrument on an Arc-drawn panel;
  per-class factors with CI; sign disagreements) are answered in one instrument report once 801's LABC lands (due 10-17 12:00Z).
  Item 4 from the record so far: no DECIDED switch's fresh-panel price disagrees in sign with its A/B/C board read (lodo,
  lodo_hi_pct, sharpen_pt all agree); 705 θ (off) and 707 await their reads.

- 2026-10-16 10:49Z: **803 bulk_transform hook delivered; 811 c0_l25 is DE-neutral on DE-H1** (one launch, r7i.8xlarge 1.54 h, $3.30;
  identity transform reproduces E0's g1pre/zs on both surfaces × 3 seeds bit-equal on every stats array, nmae ≤ 4e-15 = cell-eval2's
  last-bit noise). 811's banks were re-extracted for the 211 H1 targets with 811's own reader (source shas match; H1 group held
  out; 2–10 lines per target). c0_l25 vs inc: fid / jac / reach ≤ 4e-5 or in seed noise; **nmae improves** (f 1: −0.00044 ± 0.00018
  h1-test, −0.00045 ± 0.00005 h1-trainval) although ballast K moves +3 % / +9–10 %; bulk on the H1 line: mse better, pds −0.0094
  (h1-test). **Rule 6 is answered for 811: no DE cost.** 811 stays WEAK on its decision set (ABCp rule 1), so no dossier.
  810 sgu1 on DE-H1 not run (needs the whitened base part and a 1.3 GB operator; ~2 h; not funded: the family is closed for the
  A/B/C route).
- **Lead plan for 811 (the 6a route).** c0_l25 (λ 0.25, centring off, day bank list) is frozen as the lab's first 6a candidate,
  to be read ONCE on D/E/F-panel folds after the release window (an unworn surface, N = 1 there), with DE-H1 already neutral.
  The sealed LZ folds (L7, L8) were drawn by arc800's rule and so carry the panel bias found today; they stay sealed, but a 6a
  confirmation should rest on Arc-drawn panels (D/E/F's own, LABC) rather than on LZ alone.

- 2026-10-16 10:51Z: **814-panel-rule opened** (instrument lane; ≤ $5 before the judge): which target properties separate A/B/C's
  panel from arc800 draws, and a draw rule v2 under which A/B/C's panel is typical (bar written in the mandate: mean percentile ≥
  0.35, bottom ≤ 6/21 on new v2 folds). Already ruled out by the lead (file comparison): **essential genes**, since A/B/C's panel
  holds 0 of the 2,021 K562-essential perts, like every L panel. Lead hypothesis: arc800 bins strength on GWCD4i (CD4) DE counts
  while Arc binned on an H1 pilot; the H1-2025 proxy is the first feature 814 tests. 801 keeps LABC and the per-fold profile.
  Phase B watch is now `800-checker/phaseb_watch.sh` (finds every 8NN worker itself).

- 2026-10-16 11:14Z: **Panel instrument report sent to the coordinator** (`docs/LAB-PANEL-REPORT-2026-10-16.md`, 569978b9). It is the PANEL:
  LABC (A/B/C's targets on the fresh-fold recipe, 801, $0.17) ranks bottom in 14/21 cells (mean percentile 0.19), like ABCp;
  A/B/C's fold rules add only +0.4–0.5 pp of r. Over-read (fresh ÷ Arc-drawn, like-for-like): key **3.08 [1.91, 6.76]**,
  direction **6.5 (upper bound open; lower bound 2.0)**, shape **1.19 [0.81, 1.81]**. LABC = the def_day recipe on an Arc-drawn
  panel, so day-of folds on D/E/F's own panel should not carry the over-read. No DECIDED switch disagrees in sign with its board
  read. 810 LABC: sgu1 −0.00015, 707 +0.00115 3/3; 811 LABC: c0_l25 +0.0023 3/3 (WEAK), many-line increment +0.0036 3/3 (holds),
  sources-only blend negative. **Consequence for the lab:** BULK-W prices key and direction candidates 3–6× too high for any
  Arc-drawn panel; from now on key/direction candidates are decided on ABCp + LABC (and D/E/F-panel folds after the release),
  with BULK-W reported only. Shape candidates may still be priced on BULK-W (factor ≈ 1).

- 2026-10-16 11:18Z: **Science gets LABC (option a), under registration** (coordinator's relay, 10-16). Science prices dose partners on
  LABC with score708 (707 β 0.125 / 0.5, 705 θ 0.5, H1/CD4 design points). The dumps are already on the Mac in dump708
  format, read-only. Protocol: `800-checker/wear/LABC-PROTOCOL.md`: register configs + scorer sha + dump shas before reading;
  the lab logs declared rows and family N (Arc-drawn set {062, ABCp, LABC}: 707 starts at 3, 705 at 3, H1/CD4 at 2); winners are
  quoted at the best-of-N bound. Caveat sent: LABC is one Arc-drawn panel, so a dose chosen on it is A/B/C-panel-chosen.

- 2026-10-16 11:19Z: Heartbeat: no lab instance running; 814 alive; Data volume 93 GiB free. Lab local footprint was 12 GB (du, clones
  counted twice) against the 15 GB cap: removed the `out/repo` copies of four 801 runs under `800-checker/cloud/runs/`
  (abcp-d, abcp-e, labc: byte-identical to the repo files; designs-abcp2: the design npz, held on S3 in the input mirror and
  the run output). Status, logs and job records kept. Footprint now 8.6 GB.

- 2026-10-16 11:41Z: **814 finds the panel property: essential-screen targets.** A/B/C's panel holds **0** targets of the RPE1 ∪ Jurkat
  essential-gene screens (`data/donors/rpe1`, `nadig_jurkat`); an arc800 draw is expected to hold 16.9 ± 3.7 (A/B/C percentile
  0.000; L1–L6 hold 14–24). arc800 removed only K562's essential list. Dropping those ~17 targets from the L folds moves ABCp's
  rank from 0.194 (14/21 bottom) to **0.357 (5/21)**; direction 0.10 → 0.39. Same-count controls do less (strongest truths
  0.276, strongest predictions 0.241, random 0.207), so it is WHICH targets, not their strength. Second axis: source-bank
  responsiveness (A/B/C median at pct 1.000). Mechanism offered: Arc's cells are sorted ~14 days after the guide with ≥ 80 %
  knockdown, so fitness genes do not survive. Cost $0.041. **PROPOSAL (draw rule v2 "arc814")** to the judge
  (`819/814-panel-rule-20261016T114114Z`): v2 = arc800 minus the essential screens, with responsiveness-weighted draws; six
  independent panels (draw shas pre-recorded); validation on 18 new folds ≤ $5 against the lead's bar (mean pct ≥ 0.35, bottom ≤
  6/21). Free prediction for 10-22: D/E/F's panel holds ≤ 3 essential-screen targets. **Consequence now:** the fresh panels'
  direction over-read is mostly an essential-gene effect; candidates whose fresh-panel gain sits on those targets will not carry
  to an Arc-drawn panel.

- 2026-10-16 11:48Z: **814 PROPOSAL → REVISE (judge, `814/VERDICT-1.md`; text only), and a lead error on record.** The essential-screen
  cut is NOT a new finding. FINDINGS §8.1g / LOG l.200–212 already record that Arc's panel avoids the K562, RPE1, HepG2 and Jurkat
  essential screens, and that synthetic panels must satisfy that criterion; four earlier experiments walked past it. **The lab's
  arc800 rule (DESIGN §3.1, written by the lead) applied only the K562 list, and the lead's 814 mandate repeated the narrowing
  ("essential genes are NOT the difference"). Both were wrong.** 814 recovered the rule from the data. Fixed now: DESIGN §3.1
  amendment (every lab draw removes RPE1 ∪ Jurkat [∪ HepG2] essential-screen targets; L1–L6 and LZ violate it). Part B
  (responsiveness / KD / coverage weights) has a precedent in 522 §4b (key tails stayed heavy) and is fitted to A/B/C, so its
  validation is in-sample; the judge asks for: a novelty section, B labelled A/B/C-fitted with a second D/E/F falsifier (prop10
  decile distribution, bins now), the same-panel branch, the pool pinned by sha, and the decision role (per-class factors
  re-measured; a class below ~0.25 keeps a factor > 1). 814 resumed to write PROPOSAL-2.md.

- 2026-10-16 11:49Z: **Preregistered: A-only panels (the corrected default rule) by 801, follow-up 5** (≤ $5, by 10-17 18:00Z). Six
  independent arc800 draws with the RPE1 ∪ Jurkat essential screens also removed (seeds 801301–801306), 18 folds on the L recipe,
  then ka804 / 810 / 811 arms. **Bar, before any number:** ABCp mean percentile ≥ 0.35 and bottom ≤ 6/21 among A1–A6 (LABC
  beside). Per-class factors re-measured; a class below 0.25 keeps a factor > 1. 814's part B (weights) stays with the judge.

- 2026-10-16 11:54Z: **Coordinator's local-job HOLD released**: R7a's access.log wrote "END rc=0" at 11:53:39Z. The HOLD file is removed;
  window800 prints CLEAR. Lab local jobs resume under the usual rules (gate602 for anything beyond sub-second aggregations,
  never beside a production build). No lab local job ran during the hold (09:21–11:53Z): every aggregation went to EC2.

- 2026-10-16 11:58Z: **713 (science) LABC dose partners read and logged: no partner qualifies** (registered rule; 522 units).
  707 swap: β 0.125 − 0.25 −0.00053 [−0.00105, −0.00005] 0/3, best-of-5 bound −0.00135; β 0.5 − 0.25 −0.00085 [−0.00183, +0.00015]
  1/3, bound −0.00243. 705 θ: 0.5 − 0.25 −0.00013 [−0.00042, +0.00018] 1/3, best-of-4 bound −0.00057. References vs score708's
  incumbent: 707 β 0.25 **+0.00159 [+0.00061, +0.00261] 3/3**; 705 θ 0.25 +0.00005 [−0.00024, +0.00036] 2/3. Reading: on LABC the
  swap has an interior optimum near β 0.25 and the rider is null, consistent with the fresh-panel over-read. One panel. The lab
  checked the bound arithmetic (CI low − e(N)·sd); KA0/KA1 are science's report.

- 2026-10-16 12:00Z: **814 PROPOSAL-2 → REVISE on part B (judge, `814/VERDICT-2.md`); lead decision: B parked, no third round.** All
  VERDICT-1 items fixed. B (responsiveness / KD / coverage weights) has no falsifier of its own beyond A, names no decision it could
  change before the D/E/F-panel folds supersede it (10-22), and its A-only variant duplicated 801 follow-up 5. **P-B stays as the
  lab's frozen 10-22 test** (D/E/F's 300 targets over prop10's frozen deciles; ~$0.02; bins in 814 PROPOSAL-2 l.123–129), read
  beside P-A (D/E/F holds ≤ 3 essential-screen targets). B's folds return only as a new proposal, if P-B supports AND a decision is
  named that D/E/F-panel folds cannot serve. **DESIGN §3.1 amendment 2:** the essential-screen cut is the full screen library
  (`814/work/flags/esslib.txt`, 2,394, `24fffbfe…`), not the 2,351 donor-bank lists. 801's A1–A6 were already drawn with the 2,351
  cut (11:52Z); they are not redrawn, and their count of the 7 differing genes is reported. Lane 814: done for now (spend $0.04).

- 2026-10-16 12:02Z: 713 part 2 (science, H1 dose: A:h2.5c1, A:h3c1 vs h2c1 on LABC) registered (ee791df7) and logged (family "H1/CD4
  design" N 4). Needs LABC dumps with designs: **801 follow-up 6 queued** (`800-checker/resume/801-followup6.md`; LABCd, KA =
  non-design arrays identical; by 10-18 06:00Z), to start when follow-up 5 exits.

- 2026-10-16 12:45Z: **801 follow-up 5 done: six A-only panels, 18 folds** ($2.68; draws `draws801a.json` d4565ad2…, recorded before any
  fold; A/B/C's exact strata; 0 K562-essential, 0 donor-list essential-screen, 0 overlap with L7/L8 or A/B/C; fresh day-recipe
  fits; per-donor mean == LODO key ≤ 2.9e-14). **Recorded deviation** (DESIGN amendment 2): full-library (esslib) hits per panel
  A1 1, A2 0, A3 1, A4 2, A5 2, A6 0 (CDK7 ×3, POT1, CEBPG, CMTR2; 6 of 1,800 targets), against 14–24 per L panel. First look at the
  key arm (495, r): ABCp is still below all six A folds on K562 and HEK293T and inside the spread on HCT116. The preregistered read
  waits for the direction arms: 810 (sgu1, 707) and 811 (c0_l25) relaunched on the A folds; panelrank800 / overread800 take
  `--fresh A`. 801 follow-up 6 (LABCd for science) started.

- 2026-10-16 12:56Z: **801 follow-up 6 done: LABCd** (LABC dumps + designs h2.5c1, h3c1; run lab801-labcd, $0.13; KA: every non-design
  array bit-identical to LABC, fin == h2c1). Shas ba5433a2… / f48b21de… / bc524dc5…, sent to science for 713 part 2.

- 2026-10-16 13:00Z: **810 on the A-only folds: the direction over-read REMAINS** (price810 frozen; 288 KA keys ≤ 4.2e-15; $0.34). sgu1
  **+0.00869 [+0.0044, +0.0130] 17/18**, 707 **+0.00472 [+0.0032, +0.0062] 17/18**, against ABCp +0.00012 / +0.00178 and LABC
  −0.00015 / +0.00115: the A-only panels read like L (+0.00894 / +0.00485), not like A/B/C. So dropping the essential-screen
  targets from L folds (814's preview, 0.10 → 0.39 on direction) is not the same as drawing new panels without them: the
  replacements, drawn from the same CD4-strength strata, carry the gains too. The essential-screen cut is still Arc's rule and
  stays in DESIGN, but it does NOT explain the direction over-read by itself. The formal 21-cell read waits for 811's arm.
  sgu1's A-only gain sits in a few targets per fold (leave-top-3-out −0.0004; 707 +0.0017).

- 2026-10-16 13:12Z: **713 part 2 (science, H1 dose on LABCd) logged: no H1 candidate.** vs h2c1 (522 units): H1 ×2.5 −0.00012
  [−0.00048, +0.00028] 2/3 (best-of-4 bound −0.00068); H1 ×3 −0.00064 [−0.00133, +0.00014] 1/3 (bound −0.00172; HARM on LABC);
  paired ×2.5 → ×3 −0.00052 [−0.00085, −0.00014] 0/3. On an Arc-drawn panel, raising H1 above ×2 at CD4 ×1 does not pay, which
  agrees with 710's NULL predictions on 062 / ABCp. Bounds re-derived by the lab.

- 2026-10-16 13:19Z: Heartbeat: Data volume 72 → 73 GiB free (production builds are consuming disk); lab footprint was 15 GB by du (at
  the cap). Removed the `out/repo` copies of runs lab801-aonly and lab801-labcd (byte-identical to the repo) and aonly811. The
  latter held 811's A-panel CD4 banks (`811/work/cd4a`, 168 files), which live only on S3 at
  `s3://<your-s3-bucket>/output/lab/20261016T125901Z-aonly811/repo/experiments/811-lolo-model/work/cd4a/` (copy
  server-side to the input mirror if an EC2 job needs them). Footprint now 11 GB.

- 2026-10-16 13:21Z: **A-only panels: preregistered read FAILS** (`800/work/panelrank800-A.txt`, `overread800-*-freshA.txt`; run
  aonly-final-1016). Among A1–A6 (essential screens removed), ABCp mean percentile **0.26, bottom 9/21** (bar ≥ 0.35 and ≤ 6/21);
  LABC 0.25 / 9; ABC062 0.25 / 11. Better than among L1–L6 (0.19 / 14), not enough. Over-read factors, A-only fresh ÷ LABC: **key
  2.53 [1.53, 5.43]** (was 3.08 on L), **direction 6.0 (upper bound open)** (was 6.5), **shape 0.99 [0.68, 1.51]**; vs ABCp: key 3.09
  [1.73, 8.54], direction 4.36 [1.63, 24.5], shape 0.97. **Reading:** the essential-screen cut (Arc's stated rule, kept in DESIGN)
  removes part of the key over-read and almost none of the direction over-read. Something else about Arc's panel, or the
  idiosyncrasy of one draw, makes key and direction changes pay less on it. Only D/E/F's panel (the second Arc draw, 10-22) can
  separate the two; 814's P-A / P-B predictions are frozen for it. **Instrument policy unchanged:** key and direction candidates
  are decided on Arc-drawn panels (ABCp, LABC; D/E/F's own after the release); fresh panels (L or A) are reported with the
  factors above; shape may use fresh panels. Draw-rule family: try 1 (A-only) FAIL; a second rule needs the judge.

- 2026-10-16 13:46Z: **D/E/F-panel fold kit ready** (801 follow-up 7, $0.51). `bash experiments/801-lab-folds/kit/kit801.sh DEF
  data/challenge-def/pert_counts.csv data/challenge-def/gene_names.csv` builds, for the three cancer held lines, CD4 banks, fresh
  day-recipe operators, dump801 folds and ka804's L arm set on EC2 (~16 min per panel). `kit/make-jobs.sh <TAG>` writes the 810 /
  811 job files. **Dry-run KA (TAG LABCk on A/B/C's files):** operators bit-identical to ABCp/LABC's; dumps byte-identical to
  LABC; ka804 arrays ≤ 2.8e-13 (pricing thread count 2 vs 3). The 810/811 kit jobs are syntax-checked only; their own KAs gate
  them on the day. Note: lab folds use A/B/C's relabelled control frames; a D/E/F-controls variant is not built.
  **Day plan (after R0 is submitted and window800 is CLEAR):** run the kit on D/E/F's panel → 810 + 811 arms → panelrank / overread
  with DEF beside LABC → 814's P-A / P-B → read 811 c0_l25 once (N = 1 on DEF).

- 2026-10-16 13:51Z: **811 PROPOSAL-3 (sign skill on responding genes): no new construction, one optional $0.10 diagnostic** → judge. Every
  per-target public-line sign source is already closed (516, 621, 525, 634; 812's 12-line mean KILL). The many-line gain is
  target discrimination (pds), not DE sign (803 §8: Δfid ≤ 4e-5). Only new information (a public line near D/E/F's lineages at
  depth) could change that, and none is on disk. Diagnostic D3: skill812 on c0's day-list consensus U (u1 / a1 / a0 vs the
  incumbent, bins written; family 811-signs N 1); below 812's bar it closes the sign-skill objective for the lane. c0_l25 stays
  the frozen 6a candidate.

- 2026-10-16 13:56Z: **811 sign skill on responding genes (PROPOSAL-3 / D3): KILLED** (judge, 13:54Z; `811/VERDICT-3.md`). Proposed: no
  construction; D3 = skill812 on c0_l25's day-list consensus U. Closed because U on H1's targets is 812's `pubmean` (KILLed as
  H-U1PUB; a0 loss −0.0023 … −0.0025 already fails 812's bar), 516 already read a public same-target consensus essential-free
  on H1 (+0.0118 [−0.0007, +0.0226], did not carry to K562), and the above-bar branch would be a new bank list for a killed hook
  source. Mechanism: the many-line increment is target discrimination (pds), c0 does not touch the DE carriers (803 §8.1; DE
  members ≤ 4e-5). **Reopens only with** a pinned DE truth on a held-out non-H1 line (e.g. a K562-GWPS cell fold, 805's DE-CD4 with
  CD4 out of U, a validated KOLF2.1J cell surface) on which a centred source clears 812's bar, or a public line near D/E/F's
  lineages found after the release. **c0_l25 unchanged** (frozen 6a candidate). No budget granted.

- 2026-10-16 13:56Z: **Note on 811's reopening path (lead, after a grep):** the "K562-GWPS cell fold" DE truth it names was already built and
  closed by production as **614** (K562 real-emitter DE harness: trust gate FAILED 10-02, re-run FAILED 10-03; FINDINGS l.230–233).
  The truth has too little power (median n_real 5; 8,248 genes; jac level ~0.006 vs the board's ~0.04), and it misread G2 / G3 by
  9–22×. So a K562 DE surface is not an open route. The remaining reopeners are 805's DE-CD4 (stage-0 blob FAIL here) or a new
  public line after 10-22. No lab lane is opened on it now.

- 2026-10-16 18:18Z: Heartbeat: Data volume fell to 62 GiB free (production builds). The lab removed its local copies of 21 fold dumps whose
  S3 mirror copies match by size: A1–A6 × 3 and LABCk × 3 (`801/work/dumps/`, source `s3://…/input/lab/repo/experiments/
  801-lab-folds/work/dumps/`). Now 64 GiB free. Kept locally: LABC and LABCd (science reads them), abcp / abc062, and the three
  efp dumps (no S3 copy). Lab pricing runs on EC2, so nothing local needs the removed dumps.

- 2026-10-16 18:19Z: Coordinator: the disk dip is production's and expected. R-NEW's def_day holds its level h5ads (lvl0 + lvl2 23.9 GB, lvl3
  ~12 GB building), so the low point is **~52 GiB at ~19:30Z**, below the 60 GiB floor. run690 deletes them on rc 0, so it recovers
  to ~100 GiB before build710 / R-ABC2. **Lab action: none, beyond staying ≤ 11 GB and starting no local job** until it recovers.

- 2026-10-17 01:18Z: **Prospective row 705 graded** (board lstr705, 10-17 s1: overall −0.0002874; bulk primary Δ(pds+mse) −0.00024;
  pds −0.000308 raw, mse +0.000546 raw (worse); DE re-roll only; CLEARED). Both instruments called it right (not HARM). On the
  range: **062's +0.0001 [−0.0001, +0.0004] misses (board below); ABCp's ≈0 [−0.0003, +0.0002] hits**; ABCp also had pds's sign and
  size (−0.00039 vs −0.00031 raw), while 062 had it wrong. Both had mse's sign wrong (the board arm was re-levelled, r3, which the
  folds don't emulate). Declared non-discriminating beforehand, so this is weak evidence for ABCp, not a test. **707 (10-17 s2 /
  10-18 s1) is the discriminating read.** File: `804/work/prospective-graded.json`.

- 2026-10-17 01:32Z: **Prospective row 707 graded: the board sides with ABCp.** Board sw707 vs ALLSW (entry AYbAiHRj…): overall
  +0.0010330; **bulk primary Δ(pds+mse) +0.000842** (raw pds -0.000520, mse -0.005899 = better); DE +0.00019 (re-roll). In like
  units (the lab's rows are board overall units, not 522 units):
  - **062:** committed NULL on 707's bins (≈ +0.0002, range −0.0006 … +0.0007) → **missed**: the read is +0.00034 above the NULL bin's
    top. Cause, from the members: 062 over-read the pds loss ~5× (−0.0027 vs −0.00052 raw) against mse ~1.4×, so its net cancelled.
  - **ABCp:** committed ABOVE the NULL bin, +0.0006 … +0.0023 → **inside** (implied ρ 2.1; checker +0.00178).
  - LABC (reported before the read in the panel report, not a prospective row): +0.00115, implied ρ 1.4.
  **One read counts against 062's call; it validates neither instrument.** Lab consequence: ABCp (day-recipe operators) stays the
  A/B/C-panel decision surface beside LABC; 062's folds are reported only. With 705 (weak, same side), the record now reads 2/2 for
  ABCp's ranges and 0/2 for 062's. 710's row is still prospective (10-18 s1).

- 2026-10-17 02:21Z: **707 added to the checker's known answers; the direction class is calibrated** (run calib-707-1017, $0.01; board pair
  707-vs-ALLSW, overall +0.0010330, in `804/work/board_pairs.json`; previous grade / ka800 files kept as `*-1016.json`).
  bulk-direction ρ = checker ÷ board overall: **W 4.69, ABCp 1.72, ABC062 0.43** (062's row is now a MISS: checker +0.00044 WEAK
  vs board +0.00103). Key ρ for comparison: W 4.75, ABCp 1.44, ABC062 1.46. So ABCp reads key and direction alike (1.4–1.7×),
  fresh panels over-read both ~4.7×, and the fresh ÷ ABCp direction ratio 2.7 agrees with the panel report (LABC lower bound 2.0;
  ABCp 4.72 [1.88, 27.5]). n = 1 per class, so a range, not a law. ka800.json now projects direction-class candidates with these.

- 2026-10-17 08:20Z: Heartbeat: Data volume at **57 GiB free (below the 60 GiB floor)**, from production builds. Lab action: uploaded 801's
  efp / efprefit dumps to the S3 input mirror (they had no copy) and removed 13 local fold dumps that S3 holds at equal size:
  efp × 4, LABCd × 3, abcp × 3, abc062 × 3. Now 58 GiB. Kept locally: LABC (the dump science reads). 801's dumps dir is now 275 MB;
  the lab prices on EC2. Coordinator told.

- 2026-10-17 08:25Z: Coordinator: the dip is planned (rider690's level files, deleted at its end), and production is moving 7 superseded
  packages to S3 (+~19 GB by ~09:45Z). **The uplink is shared until ~09:45Z: lab Mac-side transfers stay < 1 GB until then**
  (EC2 ↔ S3 unaffected). The lab's ≤ 9 GB footprint is fine.

- 2026-10-17 10:25Z: **814's 10-22 reader frozen**: `experiments/814-panel-rule/predict814.py` sha256 `3c7e2658…` (inputs pinned by sha;
  constants at full precision; feature table for 20,111 genes, `feat6_allgenes.csv.gz` 4d6a3320…, EC2 $0.014). Bins: **P-A**
  (D/E/F targets in the essential-screen library) support ≤ 3, partial 4–7, fail ≥ 8; **P-B** (targets in prop10 deciles 9–10 of
  300) support ≥ 130, partial 96–129, fail ≤ 95; same-panel → VOID. KA on A/B/C: P-A 0, P-B 164 (VOID, as required). Reader
  sanity: L1 24 / 87 (both FAIL), efp 9 / 48 (both FAIL). On 10-22: `uv run python experiments/814-panel-rule/predict814.py
  data/challenge-def/pert_counts.csv`, once.

- 2026-10-17 10:27Z: **Hazard fixed: recycled worker pids.** The watch reported "805 alive". In fact 805's stored pid (9077) had been
  recycled by an unrelated `sleep 30`; 805 has not run since 10-16. windowmgr TERMs whole process groups of "extra" workers in a
  reserved window, so a bare `kill -0` liveness test could have hit a stranger's group, production's included. Fixed: every lab
  liveness test now requires the pid to be a live process **named `claude`** (windowmgr stopped, patched and relaunched detached,
  pid 9627, before tonight's 23:30Z window; lab_subagent.sh's "already running" check; phaseb_watch.sh). The stale pid files of
  the exited workers (803, 805, 810–814) were renamed `claude-run.pid.exited`.

- 2026-10-17 10:37Z: **801 follow-up 8: the kit takes `--frames <bundle>`** ($0.50). KA: `--frames data/challenge` reproduces LABC (dumps
  byte-identical; fresh fits bit-identical; ka804 ≤ 2.8e-13). On D/E/F's own frames only the control baselines, expression gates
  and β scale move in the dumps; compositions, key, truth and scored genes are shared, so the two fold variants differ only in how
  the same effect is written into and read from each context's counts. **10-22 commands, both variants:**
  `kit801.sh DEF data/challenge-def/pert_counts.csv data/challenge-def/gene_names.csv` (A/B/C frames) and
  `kit801.sh DEFo … --frames data/challenge-def` (D/E/F's own frames; ~0.66 GB Mac upload, under the 1 GB limit).

- 2026-10-17 22:33Z: **10-22 runbook + 811's 6a preregistration committed before any D/E/F data** (`docs/LAB-DEF-RUNBOOK.md`). Steps:
  gates → 814's frozen predictions → kit DEFo (D/E/F's own frames, decision surface) + DEF → 810/811 arms → panelrank800 `--arc
  DEFo --arc DEF` → 811's c0_l25 read on DEFo (N = 1; 6a dossier only if the DEFo 2.5 % bound > 0, DEF point > 0, no Arc-drawn
  HARM, the K / DE-H1 condition and the sealed LZ sign check all hold). panelrank800 gained `--arc TAG`; its known answer passed
  (LABCk == LABC on 15/15 cells, max |d| 1.7e-17, 0 rank mismatches; through gate602, seconds).

- 2026-10-17 22:35Z: **Runbook Amendment 1 (coordinator's correction): 811 is priced as an increment on top of what would ship.** 6a
  allows one attempt, and the likely one is R0 (with 705's θ) + 707's swap. New arms: B0 = inc + θ, B1 = B0 + swap, C0 = B0 + 811,
  C1 = B1 + 811 (swap first). Rule: DEFo C1 − B1 2.5 % bound > 0 and DEF point > 0 if the swap ships ON; C0 − B0 with the same bars
  if it ships OFF; pairs reported. LZ (one access) is confirmed as the lab's, via the coordinator's CONFIRMATION-LEDGER, with the
  claim sha before the first operator fit. 811 resumed to build composite811 (KAs vs ka804 and m811; information read on ABCp / LABC;
  kit job copy; ≤ $2; by 10-21 18:00Z).

- 2026-10-17 23:14Z: **811 composite arms built** ($0.60; `composite811.py` 6ec25788…; m811 / ka804 / engine801 unchanged; KAs: B0 == ka804 705,
  B1 (rider off) == ka804 707, ≤ 2.8e-13; λ 0 identities bit-exact; C0 (rider off) == stored c0_l25 0.0; kit dry run on LABCk ==
  direct LABC read, 15/15 arm files |d| 0). Information read on Arc-drawn A/B/C panels: **C1 − B1 (811 on top of R0 + swap)
  ABCp +0.00153 [−0.00096, +0.00428] 2/3, LABC +0.00142 2/3**, i.e. it would NOT pass Amendment 1's bar there. C0 − B0 (no
  swap) ABCp +0.00269 [+0.00005, +0.00564] 3/3. The swap absorbs ~45 % of 811's gain. DEFo decides (N = 1); exact commands in
  runbook Amendment 2. Wear merged (8 rows).

- 2026-10-18 01:18Z: 710's board read moved: production's plan (`NEXT-SLOT-DECISION.md` l.2624) runs **h1s710 in 10-19 s1** (10-18 s1 is
  R-ABC2, s2 ku38). The lab's 710 rows (062 NULL −0.00005; ABCp NULL −0.00004, committed 10-16) stay prospective. Before grading,
  check that h1s710 is the construction the rows priced (H1 union ×2.5 vs ×2 at A/B/C's CD4 ×1.5, key and sharpening unchanged).

- 2026-10-18 04:18Z: **The f 0.61 call-set rule is REFUTED by the board** (production's ku38, 10-18: reach rises monotonically with K 2,600 →
  3,400 → 4,200, +0 / +0.0054 / +0.0094 raw; both H1 dials, f 0.87 and f 0.61, got the direction wrong). The lab had labelled f 0.61
  post hoc and non-decisive (coordinator, 10-16); it is now closed. **Consequence:** DE-H1's de-callset class cannot price
  call-set size, and any lab candidate that moves ballast K reports its call-set effect as UNMEASURED. ka800.json's `post_hoc`
  note had been dropped by the 10-17 calib804 regeneration; it is restored with this record, and calib804 now carries `post_hoc`
  across regenerations.

- 2026-10-18 08:32Z: **Coordinator offers a board read of C1 on A/B/C** (R-ABC2H + 811's c0_l25; 10-20 s2 or 10-21 s1), needing by 10-19
  12:00Z an UNAPPLIED dump hook (811 then 707 order: swap first) with KA-off == dumphook707 byte for byte, plus a preregistered read.
  Accepted conditionally: 811 is building `dumphook811.sh` (deadline 10-19 10:00Z; NOT READY if unclean). Runbook Amendment 3: a
  HARM board read of C1 (D ≤ −0.0005 vs pred, or ≥ 0.0005 below R-ABC2H) closes 811's 6a route. Pred ≈ R-ABC2H + 0.0009
  (ABCp C1 − B1 ÷ ρ 1.44–1.72), written once R-ABC2H reads.

- 2026-10-18 08:59Z: **C1 hook READY** (811's `dumphook811.sh` 6cccd69a…, `hook811.py` c5a5614a…, committed 503725de; UNAPPLIED; EC2 $0.28).
  KAs: EC2 KA-off byte-identical to dumphook707 on the same machine; KA-same-construction |d| 0.0 on ABCp × 3 + L1 (so the
  dump-space hook IS composite811's C1, and the gap cannot flip the HARM call); 690g's validator PASS. **The lead re-ran readiness's
  exact KA-off on the Mac** (gate602, ~17 s, ≤ 1.7 GB, window CLEAR): HOOK811_LAMBDA=0 on 697's dump `2dfe327d` → **`af27efbd`**
  byte for byte (= dumphook707's production output). At λ 0.25 the same input → `83b44060…` (the C1 dump readiness should
  reproduce); no bulk_ctx. The hook reads 811's git-ignored banks on the Mac (sha-checked); they must stay on disk until C1 is built.

- 2026-10-18 20:14Z: **Reboot recovery (Mac restarted ~20:10Z).** Lab state after the restart:
  - windowmgr relaunched detached (pid 3455), with the patched liveness check;
  - the heartbeat cron survived with this session;
  - no lab instances or workers were running;
  - the C1 hook files are unchanged and committed (dumphook811.sh 6cccd69a…, hook811.py c5a5614a…), and 811's banks are intact
    (399 MB, 19 files).
  Production's own recovery was already underway: armready-1020s1 "re-armed after the reboot by the coordinator", and c1811's
  package was built READY at 11:56Z, before the reboot. Production's waiters and watchers stay production's to re-arm
  (AGENTS §6.2).

- 2026-10-19 01:31Z: **Prospective row 710 graded: both instruments' ranges contain the read** (board h1s710 UW4LsxTU…, vs ALLSW: overall
  +0.0004824, **bulk +0.000380** (mse +0.00225 scaled, pds ~0), DE +0.000103 re-roll; its own bin POSITIVE). The construction matches the lab's
  row (H1 ×2.5 vs ×2 at A/B/C's CD4 ×1.5, key and sharpening unchanged). 062 (−0.00005) and ABCp (−0.00004) both predicted NULL in −0.0005 …
  +0.0004: the read sits at the range top, so both under-read the point by ~0.0004 and both called not-HARM correctly. Non-discriminating,
  as declared. **Final tally of the committed rows:** ABCp's ranges contain 3/3 (705, 707, 710), 062's 1/3. **A pattern across 707 and 710:**
  the literal-pds leg of the A/B/C-panel folds predicts pds LOSSES the board does not show (707: board −0.00052 raw vs ABCp −0.0013, 062
  −0.0027; 710: board ~0 vs −0.0008 / −0.0007), while the mse leg has the right sign. So the checker's bulk read is pessimistic on pds
  for composition / direction changes. That matters for C1's read (10-20 s1), where pds is ~half of 811's mse gain.

- 2026-10-19 02:18Z: **C1's prediction filled in** (`811/PREREG-C1-ABC.md`, 333dddc8): R-ABC2H read **0.2748936** (entry ndt8swt7…, 01:49Z;
  production's pred 0.2750179). **pred = 0.2757814**; AS PREDICTED in [0.2752814, 0.2762814]; **HARM if C1 ≤ 0.2752814 or C1 <
  0.2743936.** Correction: commit 333dddc8's MESSAGE says 0.2757819 / 0.2752819, typed before the script printed. The FILE holds the
  computed values above, and they are the binding ones (no amend, per Amendment 1). C1 reads in 10-20 s1.

- 2026-10-20 02:37Z: **C1 read ABOVE on the board: 0.2777649** (pred 0.2757814, D +0.0019835; +0.0028713 over R-ABC2H, nearly all mse: Δmse +0.01773
  scaled, Δpds −0.00059). Not HARM, so 811's 6a route stays open (swap-ON branch). The board gave 811's increment **1.9× ABCp's read**
  (+0.00153; implied ρ 0.533 vs the 1.72 deflation used), the third case of the folds' pessimism on composition / direction changes
  (they over-read pds losses). C1 is a NEW BEST entry (0.2777649 > sw707 0.2771181). The coordinator offers a board read of **C0**
  (811 on R-ABC2's UNSWAPPED stack, the swap-OFF branch of SIXA rule 3), READY by 10-20 14:00Z for 10-20 s2, else 10-21 s1/s2.
  Lab is building it now (hook + KAs + PREREG-C0-ABC before any C0 number).

- 2026-10-20 02:48Z: **C0 hook READY** (swap-OFF branch board read, offered for 10-20 s2):
  - `811/dumphook811c0.sh` = hook811.py (unchanged, c5a5614a…) straight on the stack dump, no swap707; UNAPPLIED.
  - **PREREG-C0-ABC committed BEFORE any C0 number** (e603fe73): ρ = class median 1.1280 (707 1.7232, C1's 811 increment 0.5329);
    pred **0.2763696**; HARM if C0 ≤ 0.2758696 or < 0.2734849.
  - Mac KAs (gate602): HOOK811_LAMBDA=0 on 2dfe327d → output **byte-identical to the input** (2dfe327d…); λ 0.25 → C0 dump
    **96204332…**; 690g validator PASS (mean cosine 0.9825).
  - EC2 KA-same-construction (kasame811c0, $0.13): hook == composite811's C0 rows, |d| 0.0, same changed targets, KA-scale PASS,
    on ABCp × 3 + L1.

- 2026-10-20 08:18Z: **C0 read ABOVE on the board: 0.2780624, a new best** (pred 0.2763696, D +0.0016928; **+0.0040775 over R-ABC2**, nearly all
  mse: Δmse +0.02314, Δpds −0.00027; DE +0.00027 noise). Implied ρ 0.66 (C1: 0.53): both 811 board reads exceed the folds by
  1.5–1.9×. **The interaction:** swap alone +0.00091, 811 alone +0.00408, both +0.00378, so with 811 in, the swap adds nothing
  (C1 − C0 −0.0003). Calibration record, bulk-direction on ABCp, now three board points: 707 1.72, C1's 811 increment 0.53, C0's
  0.66 (class median 0.66). The lab's direction-class deflation should no longer be 1.72; ka800 is updated on the next calibration
  run. 811's 6a candidate is the strongest A/B/C board gain the team has seen in weeks; its 6a gate stays the DEFo read after the
  release (runbook §5 / Amendments 1 and 3).

- **2026-10-20T14:22Z — C1:C0 wired for SIXA Amendment 4 (coordinator OUTBOX 07:40Z / 09:34Z).** 811's exact 10-22 aggtag811 loop gains C1:C0
  (RESULTS §000000, `--parallel 4`). Format test on the LABCk dry run's arrays relabelled DEFo/DEF (scratchpad, agg811 unchanged, gate602, 12 MB):
  decide6a parses `C1|C0|DEFo|pds` and its per-fold values (811 §00000000). By-product, an information read: C1 − C0 on LABC +0.00023, 2/3 folds,
  vs −0.0003 on the board. So the Amendment 4 rule would keep the swap on both Arc-drawn panels; flagged to the coordinator (runbook Amendment 4
  caution + `experiments/800-checker/OUTBOX.md`). Wear row in 811's wear811.jsonl and 800's ledger (n_tries 12, decides nothing).
  The LZ mapping (C0 − B0 when DEFo's C1 − C0 fails on a swap-ON day) is accepted as recorded in CONFIRMATION-LEDGER e9b4bfd6; no addendum.

- **2026-10-20T15:20Z — SIXA Amendment 4b adopted (coordinator 1e9a9695, on the lab's 14:22Z caution).** A swap-ON day keeps the swap on top of 811 only if
  C1:C0's DEFo ci95[0] > 0; else C0 if C0 − B0 passes its bars, else C1. Verified in decide6a.py (`c1_over_c0`), the runbook's Amendment 4 and
  the coordinator's 17/17 tests. On the format-test jsons it returns False on LABC (ci95[0] −0.00052): C0. 811 RESULTS' reading lines updated.

- **2026-10-22T00:11Z — lab lead awake for the D/E/F release (an inspector session's wake-up, sent with Ahmet's OK).** The lead had been down since
  10-20 15:21Z (reboot). The release was detected 10-22 00:02Z and watch691 unpacked data/challenge-def at 00:04Z: D/E/F = MOLM13 / U87MG / SKLMS1
  by production's match, NEW-PANEL, 129 targets shared with A/B/C, gene axis identical. window800 says WINDOW: daily until 01:45Z, then
  release→R0. R0 has not started on its own (CD4 gate); readiness is launching it by hand, ETA ~4–6 h. Re-armed: windowmgr (pid 18097), the
  hourly heartbeat cron, and a background waiter that wakes the lead on CLEAR or on R0 arm trouble. Frozen shas re-checked, all match:
  predict814 3c7e2658, panelrank800 d5f9e27c, m811 e00192e6, composite811 6ec25788, agg811 d34d1e26, aggtag811 af68906a, kit801 1ad7545c,
  kit801-collect fcd52585, make-jobs 5a4da962, make-jobs811c c5b40246, run-kit811c 442e867a. Nothing has been read from data/challenge-def by the
  lab. LAB-DEF-RUNBOOK starts at CLEAR, with step 0's sha record. Note: the runbook's §0 says "whichever is earlier", but window800 rule 2 holds the
  lab until R0's SUBMITTED ok with no time cap. window800 is the gate.

- **2026-10-22T05:19Z — R0 stopped on postbuild; the lab stays held.** At 05:17:56Z def_day's postbuild failed 1 of 42 checks: 554's dose guard, D's median
  prior_lambda 0.8975 < 0.9 (E and F pass). arm1015 R0 logged NOT SUBMITTING / END. The builder was restored to e99bbd4a. This is production's,
  and the lab touches nothing. A heads-up went to the coordinator, and a note is in the lab OUTBOX. window800 stays WINDOW until R0's SUBMITTED ok;
  the lead's waiter was re-armed (wakes on CLEAR or on a new submit-R0.log line, cap 10-23 06:00Z).

- **2026-10-22T06:24Z — LAB-DEF-RUNBOOK started (R0 SUBMITTED 06:05:02Z, final entry k2OKr0zPvJq48BCtdUx8; coordinator's go at ~06:10Z).**
  - **Window.** window800 rule 3 still prints WINDOW until 14:00Z: as coded, the fixed window is not shortened by R0. The coordinator, who owns the
    window, said start. Every step from here is EC2 or a Mac transfer < 1 GB, which the window allows anyway. Local Python goes only through gate602
    (production's X1 is building).
  - **§0 shas (data/challenge-def):** context_D 97dfff32570748677fef6ee732955f71af5e7596a3747a11aac36903206a0d13,
    context_E 1942a0b6e00331da4503a9f67101a002564675bc775c2b8534a9628941436e93, context_F 5d835a036045fdf3849b5c3f8633282498aba75038ef77c999e596815a2ee38c,
    gene_names 25bfa66715e186bebabce7ac788bbcea47e2bf59ca70be1f8f3a06f2f0e47201 (== A/B/C's), manifest 7fd5a2ea490c40b6d156dd189db2c8cbddedabab13dec1c27b69089a04b34e28,
    pert_counts 6a472047d58ddb2a59b2345f962dabaf4120276c381a5e97ec19ea03f640b118 (shas pasted from `shasum` output).
  - **§1 free predictions (predict814 3c7e2658, frozen 10-16):** 129 targets shared with A/B/C (NEW-PANEL).
    - **P-A SUPPORT:** 0 essential-library targets (support ≤ 3). Arc drew D/E/F by the same essential-screen rule, so the A-only reading of the
      panel bias stands.
    - **P-B SUPPORT:** 143 / 300 in prop10 deciles 9–10 (support ≥ 130). Deciles [2, 2, 19, 22, 19, 28, 39, 24, 58, 85].
    - json `experiments/814-panel-rule/work/predict/predict814-6a472047d58ddb2a.json`.
  - **§2 kits launched:** DEFo 20261022T061527Z-lab801-kit-defo (own frames), DEF 20261022T061958Z-lab801-kit-def.
    - The first DEFo call died silently before launch: kit801's frames loop runs `aws s3 ls` on files that are not yet in the mirror, and under
      `set -e` + pipefail that exits.
    - Worked around without editing kit801: staged the four frames files with `labcloud.sh stage` (0.64 GB, sizes match) and moved the partial dir to
      `work/DEFo.aborted-20261022T0615Z`. A kit801 fix (`|| true` on that `aws s3 ls`) is a follow-up for after the day.

- **2026-10-22T07:01Z — 811's D/E/F read (runbook §5): items 1–4 PASS on the swap-ON branch; Amendment 4b → X2 = C0.** Every KA passes.
  - DEFo C1 − B1: +0.00299 [+0.00054, +0.00566], 3/3 folds. DEF: +0.00222.
  - C1 − C0 on DEFo: +0.00033, bound -0.00071 ≤ 0, so the swap does not clear it.
  - C0 − B0 on DEFo: +0.00446 [+0.00175, +0.00733], 3/3 folds. DEF: +0.00363.
  - K / DE-H1: clean. C0's max |dK| is 4.46 %. C1's is 6.09 % on HCT116, but the DE-H1 read stands and the A/B/C board reads were ABOVE.
  - decide6a (read-only): X2 = C0, LZ increment C0 − B0. Details are in 811 RESULTS §000000000; 8 wear rows (N = 1). The LZ access is next.

- **2026-10-22T07:21Z — §3 report-only arms and §4 panel read.**
  - 811's kit (m811: inc / c0_l25 / c1_l25 / src_l25) on DEFo and DEF: rc 0. Its L1 KAs and the L1 CD4 compact all hold at worst |d| 0.
    Arrays are in `811/work/run-m811kit-{DEFo,DEF}`.
  - 810's kit (sgu1 / 707 parts) FAILED at `kit810.sh prep` on both tags: "…/bundle/context_A.h5ad: No such file".
    - kit810's prep copies the D/E/F controls but lacks cloud801 prep's manifest step for kit bundles (two csv files only), so it was never kit-ready
      for a new panel. The LABCk dry run exercised 801/811 only.
    - Cost $0.02. 810 is closed for the A/B/C route and report-only here, so it is skipped for the day. The fix is a follow-up.
  - panelrank800 `--fresh A --arc DEFo --arc DEF` is queued in gate602 (ticket 496) behind production's X1 build (no local job beside a builder).
    It decides nothing and will run when the builder is done.
  - LZ kits running: LZ7 20261022T070218Z, LZ8 20261022T070608Z.

- **2026-10-22T07:48Z — LZ sign check PASS (one access, now SPENT).** **C0 − B0 (deciding): +0.00999 [+0.00659, +0.01273], 6/6 folds > 0, LT3 +0.00571** (LZ7 K562 / HCT116 / HEK293T, LZ8 K562 / HCT116 / HEK293T: +0.00929 / +0.00912 / +0.01651 / +0.00611 / +0.00811 / +0.01080; pdsl +0.01044). C1 − B1 (reported): +0.00645 [+0.00373, +0.00882], 6/6. Every KA passed first. Essential-screen counts: L7 9, L8 20.
  Runbook §5 items 1–5 all hold for X2 = C0. The 6a dossier and the ledger line go to the coordinator. The X2 banks fix (defbanks811, coordinator
  items 1–2) is running on EC2.

- **2026-10-22T08:09Z — X2 banks fix done on EC2 and installed on the Mac.** hook811.py was re-pinned to banks manifest 4cc6d82a, its only change.
  - KA-same-C0 passes on 10 folds; the priced DEFo/DEF arrays are reproduced exactly (402/402).
  - On EC2, A/B/C old vs new is byte-identical. The Mac reference differs from EC2 by ≤ 1.3e-15, which is platform arithmetic.
  - The Mac file-sha KA is queued in gate602 behind X1. Details: 811 RESULTS §00000000000.

- **2026-10-22T08:11Z — Mac file-sha KA PASS** (beside X1, coordinator's OK). Re-pinned hook + new banks, λ 0.25 on 2dfe327d → 96204332 byte for byte;
  λ 0 → input content. The lab's X2 items are complete; only readiness's preflight remains before GO-X2.

- **2026-10-22T08:13Z — GO-X2 written by the coordinator (c3f24c91).** The dossier was verified and LZ is recorded in CONFIRMATION-LEDGER (79868244).
  The coordinator committed work/banks/manifest.json at 4cc6d82a (77c50953). X2 = C0 builds after X1 (~11:30Z → ~16:30Z).
  **FREEZE:** work/banks and hook811.py stay untouched until X2 is submitted or retired (`experiments/811-lolo-model/FREEZE-X2.md`).

- **2026-10-22T10:10Z — §4 panel read (panelrank800 d5f9e27c `--fresh A --arc DEFo --arc DEF`; gate602 10:10Z, 4 s, rc 0): D/E/F's panel ranks like
  A/B/C's.** Reported only; it decides nothing.
  - Ranked among the six A-only lab draws on the same held line, the mean percentile is 0.26 on DEFo and 0.29 on DEF (8/18 cells at the bottom vs
    2.6 expected). That matches ABCp's 0.26 (9/21), LABC's 0.25 and ABC062's 0.25; exchangeable would be 0.50.
  - **Reading:** the over-read of lab draws against Arc panels is Arc's selection, not one draw's quirk. It holds even after the essential-screen cut,
    and D/E/F, drawn by Arc independently, sits where A/B/C sits. The Arc-drawn calibration (ABCp ρ) is the right one for D/E/F. The lab's fresh
    draws (L, A, LZ) over-read key and direction changes, so the size of LZ's +0.00999 is not a projection.
  - **One exception:** 811's c0_l25 on DEFo K562 ranks 7/7, above every A-only draw. On D/E/F 811 reads as strong as on the lab draws, the only
    arm/line with that pattern.
  - sgu1 has no D/E/F row (810's kit prep failed, report-only).
  - json `experiments/800-checker/work/panelrank800-A-DEFo-DEF.json`, log `work/panelrank800-def.log`.

- **2026-10-22T17:50Z — 811 SHIPS (X2 = C0, R0 + c0_l25) as the replacement final entry** (coordinator; 6a decision 63803087, READY e39d2ba6, decide6a
  SHIP = X2; X1 retired unsubmitted). The lab's dossier (DEFo/DEF read, LZ PASS, bank fix + KAs) carried it.
  - **Freeze extended to Nov 5:** 811's banks, hooks and records stay unchanged (`experiments/811-lolo-model/FREEZE-X2.md`).
  - **For the write-up:** 811 is the first lab lane that went from a fresh-fold lead (10-16) through Arc-drawn reads (ABCp, LABC), two A/B/C board
    reads (C1, C0 ABOVE) and the D/E/F read to the sealed LZ, all with frozen code.

- **2026-10-22T19:17Z — X2 = C0 (R0 + 811) is SUBMITTED as the final entry:** 18:20:14Z, entry S5RxOrm4ZbqtnPRIX8NX (production's submit-R1.log;
  members hidden in the final phase). The 811 freeze holds to Nov 5.

## Phase B plan (starts only after the checker passes its known answers)
Four idea lanes with disjoint territory, one judge. Every idea: PROPOSAL.md (mechanism, class, prediction with bins, what
would falsify it, compute) → judge VERDICT.md (GO / REVISE / KILL) → only GO gets cloud → priced by check800 on BULK-W /
DE-H1 / DE-CD4 → DOSSIER.md → one message to the coordinator.
- **B1 direction on strong targets** (706/707's lead): what makes the strong tier's direction wrong on held-out lines
  (stack dilution by H1/CD4 on strong targets, donor disagreement, gene-wise noise), priced per tier.
- **B2 a new cross-dataset model, leave-one-line-out**: every public Perturb-seq line we hold (K562 GWPS/ess, RPE1,
  Nadig HepG2/Jurkat, Orion HCT116/HEK293T, VCC25 H1, GWCD4i, kolf2 panel rows, Jiang lines), trained without the
  held-out line, handed to the checker as `predictions` rows; first baseline: consensus × context gene-sensitivity.
- **B3 the DE channel**: sign accuracy on called genes (fid 0.005 vs the field's 0.04–0.09), priced on DE-H1 and DE-CD4.
- **B4 data**: any public line or screen we lack (Pisces, 2026 releases), wired first as checker folds, then as donors.
- **Judge**: a separate headless run per proposal (`experiments/8NN/VERDICT.md`), reads the proposal against FINDINGS §2
  (closed families) and asks for the mechanism and the prediction before compute.

## Cloud launches and cost (every lab instance; budget ≤ $800 through Nov 5)
| UTC | run | type / state | cap / hours | cost | instance |
|---|---|---|---|---|---|
| 2026-10-15T22:44:21Z | 20261015T224325Z-smoke | r7i.xlarge on-demand $0.2646000000/h | cap 25 min | launch, cap $0.16 | i-02da357d313328f83 |
| 2026-10-15T22:45:52Z | 20261015T224325Z-smoke | terminated | 0.016 h | actual $0.005 | i-02da357d313328f83 |
| 2026-10-15T23:02:14Z | 20261015T230018Z-pie804 | r7i.4xlarge on-demand $1.0584000000/h | cap 100 min | launch, cap $2.0 | i-0e1ebc069f3dba50f |
| 2026-10-15T23:21:38Z | 20261015T231934Z-lab801-folds-a | r7i.8xlarge on-demand $2.1168000000/h | cap 120 min | launch, cap $4.63 | i-0df17e19d15552634 |
| 2026-10-15T23:35:19Z | 20261015T230018Z-pie804 | terminated | 0.552 h | actual $0.596 | i-0e1ebc069f3dba50f |
| 2026-10-15T23:58:16Z | 20261015T235614Z-lab801-folds-b | r7i.12xlarge on-demand $3.1752000000/h | cap 90 min | launch, cap $5.33 | i-08c504da3297d280b |
| 2026-10-15T23:58:26Z | 20261015T231934Z-lab801-folds-a | terminated | 0.614 h | actual $1.312 | i-0df17e19d15552634 |
| 2026-10-16T00:38:03Z | 20261015T235614Z-lab801-folds-b | terminated | 0.663 h | actual $2.12 | i-08c504da3297d280b |
| 2026-10-16T00:40:31Z | 20261016T003830Z-lab801-folds-c | r7i.16xlarge on-demand $4.2336000000/h | cap 75 min | launch, cap $6.03 | i-0e492095745f5c8ca |
| 2026-10-16T00:51:09Z | 20261016T003830Z-lab801-folds-c | terminated | 0.177 h | actual $0.755 | i-0e492095745f5c8ca |
| 2026-10-16T00:53:31Z | 20261016T005128Z-lab801-folds-d | r7i.16xlarge on-demand $4.2336000000/h | cap 75 min | launch, cap $6.03 | i-061726412e0521c2c |
| 2026-10-16T01:02:51Z | 20261016T005128Z-lab801-folds-d | terminated | 0.156 h | actual $0.662 | i-061726412e0521c2c |
| 2026-10-16T01:05:08Z | 20261016T010308Z-lab801-folds-e | r7i.16xlarge on-demand $4.2336000000/h | cap 75 min | launch, cap $6.03 | i-0050124561eafad51 |
| 2026-10-16T01:14:51Z | 20261016T010308Z-lab801-folds-e | terminated | 0.162 h | actual $0.691 | i-0050124561eafad51 |
| 2026-10-16T01:16:56Z | 20261016T011500Z-lab801-folds-f | r7i.16xlarge on-demand $4.2336000000/h | cap 75 min | launch, cap $6.03 | i-0df43b81da232a7ea |
| 2026-10-16T01:49:51Z | 20261016T014750Z-probe805 | r7i.xlarge on-demand $0.2646000000/h | cap 30 min | launch, cap $0.18 | i-026ea3d261ed51235 |
| 2026-10-16T01:55:11Z | 20261016T011500Z-lab801-folds-f | terminated | 0.638 h | actual $2.714 | i-0df43b81da232a7ea |
| 2026-10-16T01:57:27Z | 20261016T015337Z-de803-b1 | r7i.8xlarge on-demand $2.1168000000/h | cap 330 min | launch, cap $12.12 | i-010bbfb72dce11919 |
| 2026-10-16T01:58:28Z | 20261016T015504Z-lab801-folds-g | r7i.16xlarge on-demand $4.2336000000/h | cap 75 min | launch, cap $6.03 | i-06ad0c1ccd752d531 |
| 2026-10-16T02:11:58Z | 20261016T021043Z-s0-805 | i4i.4xlarge on-demand $1.3730000000/h | cap 240 min | launch, cap $5.75 | i-025d7e3de84e1a0f5 |
| 2026-10-16T02:16:10Z | 20261016T021043Z-s0-805 | terminated | 0.07 h | actual $0.098 | i-025d7e3de84e1a0f5 |
| 2026-10-16T02:18:15Z | 20261016T021703Z-s0b-805 | i4i.4xlarge on-demand $1.3730000000/h | cap 240 min | launch, cap $5.75 | i-039d13fe07a5b132d |
| 2026-10-16T02:19:02Z | 20261016T014750Z-probe805 | terminated | 0.475 h | actual $0.136 | i-026ea3d261ed51235 |
| 2026-10-16T02:22:43Z | 20261016T021703Z-s0b-805 | terminated | 0.074 h | actual $0.104 | i-039d13fe07a5b132d |
| 2026-10-16T02:23:55Z | 20261016T022334Z-s0c-805 | i4i.4xlarge on-demand $1.3730000000/h | cap 240 min | launch, cap $5.75 | i-05caae220adc7dc69 |
| 2026-10-16T02:26:23Z | 20261016T015504Z-lab801-folds-g | terminated | 0.465 h | actual $1.981 | i-06ad0c1ccd752d531 |
| 2026-10-16T02:28:46Z | 20261016T022732Z-lab801-ka-h | r7i.4xlarge on-demand $1.0584000000/h | cap 60 min | launch, cap $1.26 | i-04ab16d5ce928b498 |
| 2026-10-16T02:34:18Z | 20261016T022334Z-s0c-805 | terminated | 0.173 h | actual $0.242 | i-05caae220adc7dc69 |
| 2026-10-16T02:36:42Z | 20261016T023623Z-s0d-805 | i4i.4xlarge on-demand $1.3730000000/h | cap 240 min | launch, cap $5.75 | i-077c50bdf8fa06eb0 |
| 2026-10-16T02:46:00Z | 20261016T022732Z-lab801-ka-h | terminated | 0.288 h | actual $0.311 | i-04ab16d5ce928b498 |
| 2026-10-16T02:49:56Z | 20261016T024934Z-s0e-805 | i4i.4xlarge on-demand $1.3730000000/h | cap 240 min | launch, cap $5.75 | i-0c473ccdcd8e00b91 |
| 2026-10-16T02:50:13Z | 20261016T024954Z-de803-b2 | r7i.2xlarge on-demand $0.5292000000/h | cap 180 min | launch, cap $1.75 | i-09fe1bf9283ea24a5 |
| 2026-10-16T02:54:46Z | 20261016T025326Z-ka804-a | r7i.8xlarge on-demand $2.1168000000/h | cap 150 min | launch, cap $5.68 | i-09b9b49d014102d34 |
| 2026-10-16T02:56:00Z | 20261016T023623Z-s0d-805 | terminated | 0.322 h | actual $0.449 | i-077c50bdf8fa06eb0 |
| 2026-10-16T03:00:03Z | 20261016T025326Z-ka804-a | terminated | 0.088 h | actual $0.189 | i-09b9b49d014102d34 |
| 2026-10-16T03:01:49Z | 20261016T030032Z-ka804-b | r7i.8xlarge on-demand $2.1168000000/h | cap 150 min | launch, cap $5.68 | i-08634c3c8b7425ed5 |
| 2026-10-16T03:10:13Z | 20261016T030032Z-ka804-b | terminated | 0.14 h | actual $0.3 | i-08634c3c8b7425ed5 |
| 2026-10-16T03:11:34Z | 20261016T031052Z-grade804 | r7i.xlarge on-demand $0.2646000000/h | cap 25 min | launch, cap $0.16 | i-01f0cdc049827fd04 |
| 2026-10-16T03:13:38Z | 20261016T015337Z-de803-b1 | terminated | 1.27 h | actual $2.716 | i-010bbfb72dce11919 |
| 2026-10-16T03:13:44Z | 20261016T031052Z-grade804 | terminated | 0.036 h | actual $0.01 | i-01f0cdc049827fd04 |
| 2026-10-16T03:18:19Z | 20261016T031543Z-lab801-designs-i | r7i.16xlarge on-demand $4.2336000000/h | cap 60 min | launch, cap $4.96 | i-00494abd1ced97a50 |
| 2026-10-16T03:19:21Z | 20261016T031751Z-de803-b3 | r7i.4xlarge on-demand $1.0584000000/h | cap 150 min | launch, cap $2.88 | i-03cdeff3574f83d83 |
| 2026-10-16T03:19:40Z | 20261016T031543Z-lab801-designs-i | terminated | 0.004 h | actual $0.017 | i-00494abd1ced97a50 |
| 2026-10-16T03:21:17Z | 20261016T031921Z-de803-b4 | r7i.8xlarge on-demand $2.1168000000/h | cap 150 min | launch, cap $5.7 | i-090cc1ac507a4706b |
| 2026-10-16T03:22:24Z | 20261016T031833Z-lab801-designs-j | r7i.16xlarge on-demand $4.2336000000/h | cap 60 min | launch, cap $4.96 | i-0bd99b88c1fa2d473 |
| 2026-10-16T03:26:17Z | 20261016T024934Z-s0e-805 | terminated | 0.606 h | actual $0.845 | i-0c473ccdcd8e00b91 |
| 2026-10-16T03:26:24Z | 20261016T032411Z-blob805 | r7i.xlarge on-demand $0.2646000000/h | cap 45 min | launch, cap $0.25 | i-04157b33662137b98 |
| 2026-10-16T03:28:34Z | 20261016T032411Z-blob805 | terminated | 0.036 h | actual $0.01 | i-04157b33662137b98 |
| 2026-10-16T03:31:13Z | 20261016T033002Z-blob805b | r7i.xlarge on-demand $0.2646000000/h | cap 45 min | launch, cap $0.25 | i-0abe52d9fcb61215f |
| 2026-10-16T03:32:34Z | 20261016T031833Z-lab801-designs-j | terminated | 0.17 h | actual $0.722 | i-0bd99b88c1fa2d473 |
| 2026-10-16T03:34:24Z | 20261016T033002Z-blob805b | terminated | 0.053 h | actual $0.015 | i-0abe52d9fcb61215f |
| 2026-10-16T03:42:24Z | 20261016T034108Z-lab801-designs-k | r7i.2xlarge on-demand $0.5292000000/h | cap 50 min | launch, cap $0.55 | i-061e2703cd2e20add |
| 2026-10-16T03:52:10Z | 20261016T034108Z-lab801-designs-k | terminated | 0.163 h | actual $0.09 | i-061e2703cd2e20add |
| 2026-10-16T03:56:19Z | 20261016T035452Z-ka804-710 | r7i.8xlarge on-demand $2.1168000000/h | cap 90 min | launch, cap $3.55 | i-0378703457482683b |
| 2026-10-16T04:02:40Z | 20261016T035452Z-ka804-710 | terminated | 0.106 h | actual $0.226 | i-0378703457482683b |
| 2026-10-16T04:16:24Z | 20261016T031921Z-de803-b4 | terminated | 0.919 h | actual $1.966 | i-090cc1ac507a4706b |
| 2026-10-16T04:30:14Z | 20261016T031751Z-de803-b3 | terminated | 1.182 h | actual $1.276 | i-03cdeff3574f83d83 |
| 2026-10-16T04:51:40Z | 20261016T045110Z-de803-b5 | r7i.8xlarge on-demand $2.1168000000/h | cap 200 min | launch, cap $7.49 | i-0761d4a0022f0e8cf |
| 2026-10-16T05:40:03Z | 20261016T024954Z-de803-b2 | terminated | 2.831 h | actual $1.56 | i-09fe1bf9283ea24a5 |
| 2026-10-16T06:03:45Z | 20261016T045110Z-de803-b5 | terminated | 1.201 h | actual $2.569 | i-0761d4a0022f0e8cf |
| 2026-10-16T06:08:35Z | 20261016T060813Z-de803-b6 | r7i.2xlarge on-demand $0.5292000000/h | cap 120 min | launch, cap $1.19 | i-0004e3121ca9b8cbe |
| 2026-10-16T06:09:23Z | 20261016T060908Z-de803-b7 | r7i.large on-demand $0.1323000000/h | cap 40 min | launch, cap $0.13 | i-019c808e53fc8937d |
| 2026-10-16T06:11:31Z | 20261016T060908Z-de803-b7 | terminated | 0.036 h | actual $0.006 | i-019c808e53fc8937d |
| 2026-10-16T06:50:48Z | 20261016T060813Z-de803-b6 | terminated | 0.704 h | actual $0.388 | i-0004e3121ca9b8cbe |
| 2026-10-16T07:00:54Z | 20261016T070036Z-gradede804 | r7i.xlarge on-demand $0.2646000000/h | cap 25 min | launch, cap $0.16 | i-0940b2b9400b39580 |
| 2026-10-16T07:03:02Z | 20261016T070036Z-gradede804 | terminated | 0.036 h | actual $0.01 | i-0940b2b9400b39580 |
| 2026-10-16T07:09:16Z | 20261016T070853Z-de803-b8 | r7i.8xlarge on-demand $2.1168000000/h | cap 150 min | launch, cap $5.7 | i-0df0fd3974d522595 |
| 2026-10-16T07:16:52Z | 20261016T071527Z-s810-decomp-a | r7i.8xlarge on-demand $2.1168000000/h | cap 60 min | launch, cap $2.5 | i-060735625d76fefef |
| 2026-10-16T07:25:15Z | 20261016T071527Z-s810-decomp-a | terminated | 0.14 h | actual $0.299 | i-060735625d76fefef |
| 2026-10-16T07:26:05Z | 20261016T072436Z-s810-decomp-b | r7i.8xlarge on-demand $2.1168000000/h | cap 60 min | launch, cap $2.5 | i-0cf3d211250bc3bb7 |
| 2026-10-16T07:26:13Z | 20261016T072449Z-base811 | r7i.8xlarge on-demand $2.1168000000/h | cap 45 min | launch, cap $1.96 | i-0e3ada75154eecdad |
| 2026-10-16T07:35:09Z | 20261016T072449Z-base811 | terminated | 0.149 h | actual $0.319 | i-0e3ada75154eecdad |
| 2026-10-16T07:35:31Z | 20261016T072436Z-s810-decomp-b | terminated | 0.157 h | actual $0.336 | i-0cf3d211250bc3bb7 |
| 2026-10-16T07:37:52Z | 20261016T073629Z-s810-dbg | r7i.xlarge on-demand $0.2646000000/h | cap 30 min | launch, cap $0.19 | i-0190d04d300f2f5fc |
| 2026-10-16T07:43:15Z | 20261016T070853Z-de803-b8 | terminated | 0.566 h | actual $1.211 | i-0df0fd3974d522595 |
| 2026-10-16T07:44:12Z | 20261016T073629Z-s810-dbg | terminated | 0.106 h | actual $0.03 | i-0190d04d300f2f5fc |
| 2026-10-16T07:44:26Z | 20261016T074407Z-de812-b1 | r7i.xlarge on-demand $0.2646000000/h | cap 40 min | launch, cap $0.23 | i-050d309cd0f6a429c |
| 2026-10-16T07:45:44Z | 20261016T074418Z-s810-decomp-c | r7i.2xlarge on-demand $0.5292000000/h | cap 40 min | launch, cap $0.46 | i-0bba8a1a9fb11e39f |
| 2026-10-16T07:47:11Z | 20261016T074655Z-agg811 | c7i.4xlarge on-demand $0.7140000000/h | cap 30 min | launch, cap $0.48 | i-04f68cf9d0e0f1898 |
| 2026-10-16T07:47:37Z | 20261016T074717Z-de812-b2 | r7i.xlarge on-demand $0.2646000000/h | cap 40 min | launch, cap $0.23 | i-0d19ce7f181b89762 |
| 2026-10-16T07:50:07Z | 20261016T074655Z-agg811 | terminated | 0.049 h | actual $0.036 | i-04f68cf9d0e0f1898 |
| 2026-10-16T07:50:46Z | 20261016T075028Z-de812-b3 | r7i.xlarge on-demand $0.2646000000/h | cap 40 min | launch, cap $0.23 | i-00de348c9c377ccdd |
| 2026-10-16T07:52:04Z | 20261016T075142Z-align811 | r7i.xlarge on-demand $0.2646000000/h | cap 40 min | launch, cap $0.23 | i-06d920b868388adb8 |
| 2026-10-16T07:53:20Z | 20261016T074418Z-s810-decomp-c | terminated | 0.122 h | actual $0.067 | i-0bba8a1a9fb11e39f |
| 2026-10-16T07:59:16Z | 20261016T075850Z-de812-b4 | r7i.xlarge on-demand $0.2646000000/h | cap 40 min | launch, cap $0.23 | i-0493977dd5d4f5fcb |
| 2026-10-16T08:00:11Z | 20261016T075142Z-align811 | terminated | 0.135 h | actual $0.039 | i-06d920b868388adb8 |
| 2026-10-16T08:00:59Z | 20261016T080038Z-align811b | r7i.xlarge on-demand $0.2646000000/h | cap 40 min | launch, cap $0.23 | i-0bb91c5e1d6ec9059 |
| 2026-10-16T08:03:18Z | 20261016T080249Z-de812-b5 | r7i.xlarge on-demand $0.2646000000/h | cap 40 min | launch, cap $0.23 | i-0404ba6404a17569e |
| 2026-10-16T08:05:36Z | 20261016T080512Z-align811c | r7i.xlarge on-demand $0.2646000000/h | cap 40 min | launch, cap $0.23 | i-072327032f3859e0d |
| 2026-10-16T08:06:36Z | 20261016T080623Z-de812-b6 | r7i.xlarge on-demand $0.2646000000/h | cap 30 min | launch, cap $0.18 | i-060a5691f9e384700 |
| 2026-10-16T08:08:41Z | 20261016T080512Z-align811c | terminated | 0.052 h | actual $0.015 | i-072327032f3859e0d |
| 2026-10-16T08:08:54Z | 20261016T080732Z-lab801-abcp-a | r7i.8xlarge on-demand $2.1168000000/h | cap 75 min | launch, cap $3.03 | i-0f1fc2bc50dee6d9a |
| 2026-10-16T08:08:57Z | 20261016T074407Z-de812-b1 | terminated | 0.409 h | actual $0.117 | i-050d309cd0f6a429c |
| 2026-10-16T08:09:02Z | 20261016T074717Z-de812-b2 | terminated | 0.357 h | actual $0.102 | i-0d19ce7f181b89762 |
| 2026-10-16T08:09:06Z | 20261016T075028Z-de812-b3 | terminated | 0.306 h | actual $0.088 | i-00de348c9c377ccdd |
| 2026-10-16T08:09:11Z | 20261016T075850Z-de812-b4 | terminated | 0.165 h | actual $0.047 | i-0493977dd5d4f5fcb |
| 2026-10-16T08:09:16Z | 20261016T080249Z-de812-b5 | terminated | 0.1 h | actual $0.029 | i-0404ba6404a17569e |
| 2026-10-16T08:09:20Z | 20261016T080623Z-de812-b6 | terminated | 0.046 h | actual $0.013 | i-060a5691f9e384700 |
| 2026-10-16T08:09:29Z | 20261016T080038Z-align811b | terminated | 0.142 h | actual $0.041 | i-0bb91c5e1d6ec9059 |
| 2026-10-16T08:21:46Z | 20261016T081954Z-s810-price-a | r7i.8xlarge on-demand $2.1168000000/h | cap 30 min | launch, cap $1.43 | i-031e65f97f8d744ae |
| 2026-10-16T08:26:09Z | 20261016T082526Z-de812-b7 | r7i.xlarge on-demand $0.2646000000/h | cap 60 min | launch, cap $0.32 | i-017f64619e4984740 |
| 2026-10-16T08:27:38Z | 20261016T081954Z-s810-price-a | terminated | 0.089 h | actual $0.191 | i-031e65f97f8d744ae |
| 2026-10-16T08:28:26Z | 20261016T080732Z-lab801-abcp-a | terminated | 0.326 h | actual $0.697 | i-0f1fc2bc50dee6d9a |
| 2026-10-16T08:31:16Z | 20261016T082526Z-de812-b7 | terminated | 0.086 h | actual $0.025 | i-017f64619e4984740 |
| 2026-10-16T08:31:19Z | 20261016T082954Z-s810-price-b | r7i.8xlarge on-demand $2.1168000000/h | cap 30 min | launch, cap $1.43 | i-01afbaef915fe32b9 |
| 2026-10-16T08:32:22Z | 20261016T083058Z-lab801-abcp-b | r7i.4xlarge on-demand $1.0584000000/h | cap 50 min | launch, cap $1.08 | i-0559348f80ed2f209 |
| 2026-10-16T08:34:55Z | 20261016T083336Z-m811s1 | r7i.8xlarge on-demand $2.1168000000/h | cap 40 min | launch, cap $1.78 | i-0934b97ac8f189b29 |
| 2026-10-16T08:38:55Z | 20261016T082954Z-s810-price-b | terminated | 0.127 h | actual $0.272 | i-01afbaef915fe32b9 |
| 2026-10-16T08:40:37Z | 20261016T083336Z-m811s1 | terminated | 0.095 h | actual $0.204 | i-0934b97ac8f189b29 |
| 2026-10-16T08:41:14Z | 20261016T084059Z-agg811s1 | c7i.4xlarge on-demand $0.7140000000/h | cap 30 min | launch, cap $0.48 | i-03ee1c5719501e433 |
| 2026-10-16T08:42:22Z | 20261016T083058Z-lab801-abcp-b | terminated | 0.167 h | actual $0.18 | i-0559348f80ed2f209 |
| 2026-10-16T08:44:49Z | 20261016T084059Z-agg811s1 | terminated | 0.06 h | actual $0.044 | i-03ee1c5719501e433 |
| 2026-10-16T08:45:04Z | 20261016T084311Z-s810-price-c | r7i.xlarge on-demand $0.2646000000/h | cap 30 min | launch, cap $0.19 | i-01bbb36269554227c |
| 2026-10-16T08:46:28Z | 20261016T084511Z-m811s2 | r7i.8xlarge on-demand $2.1168000000/h | cap 40 min | launch, cap $1.78 | i-00c03dfca44ce7256 |
| 2026-10-16T08:47:11Z | 20261016T084540Z-lab801-abcp-c | r7i.4xlarge on-demand $1.0584000000/h | cap 50 min | launch, cap $1.08 | i-067c328d7031cbad6 |
| 2026-10-16T08:51:47Z | 20261016T084511Z-m811s2 | terminated | 0.089 h | actual $0.19 | i-00c03dfca44ce7256 |
| 2026-10-16T08:51:50Z | 20261016T084311Z-s810-price-c | terminated | 0.104 h | actual $0.03 | i-01bbb36269554227c |
| 2026-10-16T08:52:21Z | 20261016T085206Z-agg811s2 | c7i.4xlarge on-demand $0.7140000000/h | cap 30 min | launch, cap $0.48 | i-03a6cffbcbc4a9be5 |
| 2026-10-16T08:54:28Z | 20261016T084540Z-lab801-abcp-c | terminated | 0.122 h | actual $0.132 | i-067c328d7031cbad6 |
| 2026-10-16T08:55:09Z | 20261016T085206Z-agg811s2 | terminated | 0.047 h | actual $0.035 | i-03a6cffbcbc4a9be5 |
| 2026-10-16T08:56:16Z | 20261016T085444Z-lab801-abcp-d | r7i.4xlarge on-demand $1.0584000000/h | cap 50 min | launch, cap $1.08 | i-0969bf13f82bdfe7b |
| 2026-10-16T09:04:42Z | 20261016T085444Z-lab801-abcp-d | terminated | 0.141 h | actual $0.152 | i-0969bf13f82bdfe7b |
| 2026-10-16T09:06:50Z | 20261016T090515Z-lab801-abcp-e | r7i.4xlarge on-demand $1.0584000000/h | cap 50 min | launch, cap $1.08 | i-00834cd516591e3c2 |
| 2026-10-16T09:10:48Z | 20261016T091016Z-de803-b9 | r7i.8xlarge on-demand $2.1168000000/h | cap 180 min | launch, cap $6.77 | i-06b460fe3905bdf52 |
| 2026-10-16T09:15:38Z | 20261016T090515Z-lab801-abcp-e | terminated | 0.147 h | actual $0.159 | i-00834cd516591e3c2 |
| 2026-10-16T09:20:16Z | 20261016T091854Z-m811s2b | r7i.8xlarge on-demand $2.1168000000/h | cap 30 min | launch, cap $1.43 | i-0f2eb49fac987caed |
| 2026-10-16T09:22:25Z | 20261016T092043Z-s810-abcp-a | r7i.2xlarge on-demand $0.5292000000/h | cap 40 min | launch, cap $0.46 | i-04a9f0f20e3c029f1 |
| 2026-10-16T09:22:33Z | 20261016T092059Z-ka804-710p | r7i.4xlarge on-demand $1.0584000000/h | cap 40 min | launch, cap $0.9 | i-0ead74686dd9029b5 |
| 2026-10-16T09:25:59Z | 20261016T091854Z-m811s2b | terminated | 0.095 h | actual $0.204 | i-0f2eb49fac987caed |
| 2026-10-16T09:26:26Z | 20261016T092610Z-agg811s2b | c7i.2xlarge on-demand $0.3570000000/h | cap 20 min | launch, cap $0.18 | i-0bf70c1a0694b2dcf |
| 2026-10-16T09:27:39Z | 20261016T092059Z-ka804-710p | terminated | 0.085 h | actual $0.092 | i-0ead74686dd9029b5 |
| 2026-10-16T09:29:06Z | 20261016T092610Z-agg811s2b | terminated | 0.045 h | actual $0.017 | i-0bf70c1a0694b2dcf |
| 2026-10-16T09:31:02Z | 20261016T092938Z-lab801-designs-abcp | r7i.2xlarge on-demand $0.5292000000/h | cap 45 min | launch, cap $0.51 | i-041c7ed8e8551168a |
| 2026-10-16T09:31:48Z | 20261016T092043Z-s810-abcp-a | terminated | 0.157 h | actual $0.086 | i-04a9f0f20e3c029f1 |
| 2026-10-16T09:32:04Z | 20261016T092938Z-lab801-designs-abcp | terminated | 0.003 h | actual $0.002 | i-041c7ed8e8551168a |
| 2026-10-16T09:32:39Z | 20261016T093113Z-lab801-designs-abcp2 | r7i.2xlarge on-demand $0.5292000000/h | cap 45 min | launch, cap $0.51 | i-0c36d513d28b87e53 |
| 2026-10-16T09:38:04Z | 20261016T093745Z-s0a813 | r7i.xlarge on-demand $0.2646000000/h | cap 60 min | launch, cap $0.32 | i-0938b5cbdae0e0b44 |
| 2026-10-16T09:41:20Z | 20261016T093113Z-lab801-designs-abcp2 | terminated | 0.145 h | actual $0.08 | i-0c36d513d28b87e53 |
| 2026-10-16T09:43:58Z | 20261016T094216Z-ka804-710p-b | r7i.4xlarge on-demand $1.0584000000/h | cap 40 min | launch, cap $0.9 | i-01eade716f716c908 |
| 2026-10-16T09:49:15Z | 20261016T094216Z-ka804-710p-b | terminated | 0.088 h | actual $0.095 | i-01eade716f716c908 |
| 2026-10-16T09:57:43Z | 20261016T093745Z-s0a813 | terminated | 0.328 h | actual $0.094 | i-0938b5cbdae0e0b44 |
| 2026-10-16T10:03:00Z | 20261016T100242Z-s0a813v2 | r7i.xlarge on-demand $0.2646000000/h | cap 40 min | launch, cap $0.23 | i-0828e9366fc6f360f |
| 2026-10-16T10:29:55Z | 20261016T100242Z-s0a813v2 | terminated | 0.449 h | actual $0.129 | i-0828e9366fc6f360f |
| 2026-10-16T10:34:41Z | 20261016T103421Z-agg-abcp-1016 | r7i.xlarge on-demand $0.2646000000/h | cap 30 min | launch, cap $0.18 | i-00bd336626faabce9 |
| 2026-10-16T10:36:22Z | 20261016T103421Z-agg-abcp-1016 | terminated | 0.028 h | actual $0.008 | i-00bd336626faabce9 |
| 2026-10-16T10:36:58Z | 20261016T103636Z-agg-abcp-1016b | r7i.xlarge on-demand $0.2646000000/h | cap 30 min | launch, cap $0.18 | i-006b0ee4694cfd336 |
| 2026-10-16T10:38:39Z | 20261016T103636Z-agg-abcp-1016b | terminated | 0.019 h | actual $0.005 | i-006b0ee4694cfd336 |
| 2026-10-16T10:43:26Z | 20261016T091016Z-de803-b9 | terminated | 1.544 h | actual $3.302 | i-06b460fe3905bdf52 |
| 2026-10-16T10:43:45Z | 20261016T104318Z-overread-1016 | r7i.xlarge on-demand $0.2646000000/h | cap 40 min | launch, cap $0.23 | i-0aca6bc80716358ef |
| 2026-10-16T10:44:16Z | 20261016T104243Z-lab801-labc | r7i.4xlarge on-demand $1.0584000000/h | cap 60 min | launch, cap $1.26 | i-002f7c1a07f111cd9 |
| 2026-10-16T10:46:00Z | 20261016T104318Z-overread-1016 | terminated | 0.038 h | actual $0.011 | i-0aca6bc80716358ef |
| 2026-10-16T10:53:42Z | 20261016T104243Z-lab801-labc | terminated | 0.158 h | actual $0.17 | i-002f7c1a07f111cd9 |
| 2026-10-16T10:56:24Z | 20261016T105604Z-labc-pre-1016 | r7i.xlarge on-demand $0.2646000000/h | cap 40 min | launch, cap $0.23 | i-0c3d796c60215dd27 |
| 2026-10-16T10:58:36Z | 20261016T105604Z-labc-pre-1016 | terminated | 0.037 h | actual $0.011 | i-0c3d796c60215dd27 |
| 2026-10-16T10:59:02Z | 20261016T105743Z-s810-labc-a | r7i.2xlarge on-demand $0.5292000000/h | cap 40 min | launch, cap $0.46 | i-0de961b254024f15d |
| 2026-10-16T10:59:12Z | 20261016T105755Z-labc811 | r7i.8xlarge on-demand $2.1168000000/h | cap 15 min | launch, cap $0.89 | i-0103a298561dc85aa |
| 2026-10-16T11:04:10Z | 20261016T105755Z-labc811 | terminated | 0.07 h | actual $0.149 | i-0103a298561dc85aa |
| 2026-10-16T11:04:43Z | 20261016T110424Z-agg811labc | c7i.2xlarge on-demand $0.3570000000/h | cap 20 min | launch, cap $0.18 | i-07ec9f6ea3501c23f |
| 2026-10-16T11:06:22Z | 20261016T110536Z-feat814a | r7i.xlarge on-demand $0.2646000000/h | cap 40 min | launch, cap $0.23 | i-0896f19f559fbcda9 |
| 2026-10-16T11:06:54Z | 20261016T110424Z-agg811labc | terminated | 0.036 h | actual $0.014 | i-07ec9f6ea3501c23f |
| 2026-10-16T11:07:42Z | 20261016T105743Z-s810-labc-a | terminated | 0.145 h | actual $0.08 | i-0de961b254024f15d |
| 2026-10-16T11:09:37Z | 20261016T110536Z-feat814a | terminated | 0.054 h | actual $0.016 | i-0896f19f559fbcda9 |
| 2026-10-16T11:11:01Z | 20261016T111040Z-labc-final-1016 | r7i.xlarge on-demand $0.2646000000/h | cap 40 min | launch, cap $0.23 | i-006c8a044ab19ba37 |
| 2026-10-16T11:13:13Z | 20261016T111040Z-labc-final-1016 | terminated | 0.037 h | actual $0.011 | i-006c8a044ab19ba37 |
| 2026-10-16T11:18:55Z | 20261016T111838Z-sim814a | r7i.large on-demand $0.1323000000/h | cap 25 min | launch, cap $0.08 | i-0c2e41c33eb130b2d |
| 2026-10-16T11:21:07Z | 20261016T111838Z-sim814a | terminated | 0.037 h | actual $0.006 | i-0c2e41c33eb130b2d |
| 2026-10-16T11:22:53Z | 20261016T112237Z-sim814b | r7i.large on-demand $0.1323000000/h | cap 25 min | launch, cap $0.08 | i-0dc5390642bd425ea |
| 2026-10-16T11:24:37Z | 20261016T112237Z-sim814b | terminated | 0.02 h | actual $0.003 | i-0dc5390642bd425ea |
| 2026-10-16T11:27:24Z | 20261016T112704Z-sim814c | r7i.large on-demand $0.1323000000/h | cap 25 min | launch, cap $0.08 | i-07b587796a4b8bb54 |
| 2026-10-16T11:29:39Z | 20261016T112704Z-sim814c | terminated | 0.038 h | actual $0.006 | i-07b587796a4b8bb54 |
| 2026-10-16T11:30:48Z | 20261016T113030Z-sim814d | r7i.large on-demand $0.1323000000/h | cap 25 min | launch, cap $0.08 | i-035f4f24b8a6b0eea |
| 2026-10-16T11:33:03Z | 20261016T113030Z-sim814d | terminated | 0.038 h | actual $0.006 | i-035f4f24b8a6b0eea |
| 2026-10-16T11:34:13Z | 20261016T113355Z-sim814e | r7i.large on-demand $0.1323000000/h | cap 25 min | launch, cap $0.08 | i-07e94e2b07b37e35e |
| 2026-10-16T11:35:48Z | 20261016T113355Z-sim814e | terminated | 0.027 h | actual $0.004 | i-07e94e2b07b37e35e |
| 2026-10-16T11:54:50Z | 20261016T115240Z-lab801-draw-a | r7i.xlarge on-demand $0.2646000000/h | cap 30 min | launch, cap $0.19 | i-09c5527011e9efd75 |
| 2026-10-16T12:01:53Z | 20261016T115240Z-lab801-draw-a | terminated | 0.117 h | actual $0.034 | i-09c5527011e9efd75 |
| 2026-10-16T12:04:39Z | 20261016T120234Z-lab801-aonly | r7i.16xlarge on-demand $4.2336000000/h | cap 75 min | launch, cap $6.04 | i-0e278dea983945514 |
| 2026-10-16T12:41:56Z | 20261016T120234Z-lab801-aonly | terminated | 0.622 h | actual $2.647 | i-0e278dea983945514 |
| 2026-10-16T12:47:14Z | 20261016T124527Z-lab801-labcd | r7i.4xlarge on-demand $1.0584000000/h | cap 40 min | launch, cap $0.9 | i-05ba8095861de06ca |
| 2026-10-16T12:48:08Z | 20261016T124635Z-s810-aonly-a | r7i.8xlarge on-demand $2.1168000000/h | cap 45 min | launch, cap $1.97 | i-057be5be9f86f0002 |
| 2026-10-16T12:54:40Z | 20261016T124527Z-lab801-labcd | terminated | 0.124 h | actual $0.134 | i-05ba8095861de06ca |
| 2026-10-16T12:57:34Z | 20261016T124635Z-s810-aonly-a | terminated | 0.157 h | actual $0.337 | i-057be5be9f86f0002 |
| 2026-10-16T13:00:50Z | 20261016T125901Z-aonly811 | r7i.8xlarge on-demand $2.1168000000/h | cap 40 min | launch, cap $1.79 | i-0d7440dca74170bfb |
| 2026-10-16T13:08:49Z | 20261016T125901Z-aonly811 | terminated | 0.133 h | actual $0.285 | i-0d7440dca74170bfb |
| 2026-10-16T13:09:35Z | 20261016T130915Z-agg811aonly | c7i.4xlarge on-demand $0.7140000000/h | cap 20 min | launch, cap $0.36 | i-019429be8f8c38d83 |
| 2026-10-16T13:12:17Z | 20261016T130915Z-agg811aonly | terminated | 0.045 h | actual $0.033 | i-019429be8f8c38d83 |
| 2026-10-16T13:12:56Z | 20261016T131231Z-agg811aonly2 | c7i.4xlarge on-demand $0.7140000000/h | cap 20 min | launch, cap $0.36 | i-08350b2dea0bc64de |
| 2026-10-16T13:15:56Z | 20261016T131231Z-agg811aonly2 | terminated | 0.05 h | actual $0.037 | i-08350b2dea0bc64de |
| 2026-10-16T13:17:53Z | 20261016T131726Z-aonly-final-1016 | r7i.xlarge on-demand $0.2646000000/h | cap 40 min | launch, cap $0.23 | i-016a1dae494579593 |
| 2026-10-16T13:21:07Z | 20261016T131726Z-aonly-final-1016 | terminated | 0.054 h | actual $0.016 | i-016a1dae494579593 |
| 2026-10-16T13:29:52Z | 20261016T132753Z-lab801-kit-labck | r7i.8xlarge on-demand $2.1168000000/h | cap 60 min | launch, cap $2.51 | i-07287b1bd11f24220 |
| 2026-10-16T13:44:08Z | 20261016T132753Z-lab801-kit-labck | terminated | 0.238 h | actual $0.509 | i-07287b1bd11f24220 |
| 2026-10-17T02:18:41Z | 20261017T021825Z-calib-707-1017 | r7i.xlarge on-demand $0.2646000000/h | cap 30 min | launch, cap $0.18 | i-0d38c8033ca5753c8 |
| 2026-10-17T02:20:50Z | 20261017T021825Z-calib-707-1017 | terminated | 0.036 h | actual $0.01 | i-0d38c8033ca5753c8 |
| 2026-10-17T10:21:13Z | 20261017T102050Z-tab814p | r7i.xlarge on-demand $0.2646000000/h | cap 30 min | launch, cap $0.18 | i-01f845f12a6e742a3 |
| 2026-10-17T10:21:58Z | 20261017T102000Z-lab801-kit-labcf | r7i.8xlarge on-demand $2.1168000000/h | cap 60 min | launch, cap $2.51 | i-01f9450c0c0bfa474 |
| 2026-10-17T10:24:44Z | 20261017T102050Z-tab814p | terminated | 0.05 h | actual $0.014 | i-01f845f12a6e742a3 |
| 2026-10-17T10:35:53Z | 20261017T102000Z-lab801-kit-labcf | terminated | 0.232 h | actual $0.497 | i-01f9450c0c0bfa474 |
| 2026-10-17T22:41:57Z | 20261017T224004Z-comp811 | r7i.8xlarge on-demand $2.1168000000/h | cap 25 min | launch, cap $1.25 | i-0d904f47e22ae50ac |
| 2026-10-17T22:51:29Z | 20261017T224004Z-comp811 | terminated | 0.159 h | actual $0.34 | i-0d904f47e22ae50ac |
| 2026-10-17T22:52:28Z | 20261017T225152Z-agg811comp | c7i.2xlarge on-demand $0.3570000000/h | cap 20 min | launch, cap $0.18 | i-05f3c864e28363cfd |
| 2026-10-17T22:58:24Z | 20261017T225259Z-lab811c-kit-labck | r7i.4xlarge on-demand $1.0584000000/h | cap 45 min | launch, cap $0.99 | i-0a458b59e1ce90278 |
| 2026-10-17T22:59:38Z | 20261017T225152Z-agg811comp | terminated | 0.12 h | actual $0.045 | i-05f3c864e28363cfd |
| 2026-10-17T23:09:26Z | 20261017T225259Z-lab811c-kit-labck | terminated | 0.184 h | actual $0.199 | i-0a458b59e1ce90278 |
| 2026-10-17T23:10:29Z | 20261017T230951Z-aggtag811-labck | c7i.2xlarge on-demand $0.3570000000/h | cap 20 min | launch, cap $0.18 | i-091e4762d865fbf0f |
| 2026-10-17T23:12:49Z | 20261017T230951Z-aggtag811-labck | terminated | 0.039 h | actual $0.015 | i-091e4762d865fbf0f |
| 2026-10-18T08:38:33Z | 20261018T083702Z-kahook811 | r7i.4xlarge on-demand $1.0584000000/h | cap 40 min | launch, cap $0.9 | i-0bf9186f38a982426 |
| 2026-10-18T08:45:54Z | 20261018T083702Z-kahook811 | terminated | 0.123 h | actual $0.133 | i-0bf9186f38a982426 |
| 2026-10-18T08:48:19Z | 20261018T084656Z-kahook811b | r7i.4xlarge on-demand $1.0584000000/h | cap 40 min | launch, cap $0.9 | i-0c564c9ab510829bc |
| 2026-10-18T08:56:15Z | 20261018T084656Z-kahook811b | terminated | 0.132 h | actual $0.143 | i-0c564c9ab510829bc |
| 2026-10-20T02:41:00Z | 20261020T023935Z-kahook811c0 | r7i.4xlarge on-demand $1.0584000000/h | cap 40 min | launch, cap $0.89 | i-015ea5195c4fa9245 |
| 2026-10-20T02:48:23Z | 20261020T023935Z-kahook811c0 | terminated | 0.123 h | actual $0.133 | i-015ea5195c4fa9245 |
| 2026-10-20T08:28:55Z | 20261020T082833Z-calib-c0c1-1020 | r7i.xlarge on-demand $0.2646000000/h | cap 20 min | launch, cap $0.13 | i-0f2b72d9181c5023c |
| 2026-10-20T08:30:05Z | 20261020T082833Z-calib-c0c1-1020 | terminated | 0.02 h | actual $0.006 | i-0f2b72d9181c5023c |
| 2026-10-22T06:19:55Z | 20261022T061527Z-lab801-kit-defo | r7i.8xlarge on-demand $2.1168000000/h | cap 60 min | launch, cap $2.51 | i-05954be0a164ef0fd |
| 2026-10-22T06:24:01Z | 20261022T061958Z-lab801-kit-def | r7i.8xlarge on-demand $2.1168000000/h | cap 60 min | launch, cap $2.51 | i-084077e5ea2634dcd |
| 2026-10-22T06:33:20Z | 20261022T061527Z-lab801-kit-defo | terminated | 0.224 h | actual $0.479 | i-05954be0a164ef0fd |
| 2026-10-22T06:38:05Z | 20261022T061958Z-lab801-kit-def | terminated | 0.235 h | actual $0.502 | i-084077e5ea2634dcd |
| 2026-10-22T06:38:22Z | 20261022T061527Z-lab801-kit-defo | terminated | 0.308 h | actual $0.659 | i-05954be0a164ef0fd |
| 2026-10-22T06:39:06Z | 20261022T061958Z-lab801-kit-def | terminated | 0.251 h | actual $0.538 | i-084077e5ea2634dcd |
| 2026-10-22T06:43:11Z | 20261022T064006Z-lab811c-kit-defo | r7i.4xlarge on-demand $1.0584000000/h | cap 45 min | launch, cap $0.99 | i-0e6eff2935290d08d |
| 2026-10-22T06:46:09Z | 20261022T064311Z-lab811c-kit-def | r7i.4xlarge on-demand $1.0584000000/h | cap 45 min | launch, cap $0.99 | i-0a5958ae064ae0a47 |
| 2026-10-22T06:48:51Z | 20261022T064617Z-lab810-kit-defo | r7i.4xlarge on-demand $1.0584000000/h | cap 60 min | launch, cap $1.26 | i-077ea1bc90813528a |
| 2026-10-22T06:51:48Z | 20261022T064851Z-lab811-kit-defo | r7i.4xlarge on-demand $1.0584000000/h | cap 60 min | launch, cap $1.26 | i-06e90d543fa2d2813 |
| 2026-10-22T06:54:23Z | 20261022T065148Z-lab810-kit-def | r7i.4xlarge on-demand $1.0584000000/h | cap 60 min | launch, cap $1.26 | i-0e772d2e90756fd8d |
| 2026-10-22T06:56:41Z | 20261022T064006Z-lab811c-kit-defo | terminated | 0.225 h | actual $0.243 | i-0e6eff2935290d08d |
| 2026-10-22T06:56:52Z | 20261022T064311Z-lab811c-kit-def | terminated | 0.179 h | actual $0.193 | i-0a5958ae064ae0a47 |
| 2026-10-22T06:57:19Z | 20261022T065423Z-lab811-kit-def | r7i.4xlarge on-demand $1.0584000000/h | cap 60 min | launch, cap $1.26 | i-09d650382985d1d5f |
| 2026-10-22T06:57:52Z | 20261022T065734Z-aggtag811-defo | c7i.2xlarge on-demand $0.3570000000/h | cap 20 min | launch, cap $0.18 | i-06e0fa5f9a3af276c |
| 2026-10-22T06:58:16Z | 20261022T065756Z-aggtag811-def | c7i.2xlarge on-demand $0.3570000000/h | cap 20 min | launch, cap $0.18 | i-02954ccf7a09a454d |
| 2026-10-22T06:59:31Z | 20261022T065734Z-aggtag811-defo | terminated | 0.028 h | actual $0.01 | i-06e0fa5f9a3af276c |
| 2026-10-22T06:59:40Z | 20261022T065756Z-aggtag811-def | terminated | 0.024 h | actual $0.009 | i-02954ccf7a09a454d |
| 2026-10-22T07:06:05Z | 20261022T070218Z-lab801-kit-lz7 | r7i.8xlarge on-demand $2.1168000000/h | cap 60 min | launch, cap $2.51 | i-0cadbc0e1c15aa382 |
| 2026-10-22T07:10:00Z | 20261022T070608Z-lab801-kit-lz8 | r7i.8xlarge on-demand $2.1168000000/h | cap 60 min | launch, cap $2.51 | i-0ad1f2314bed9d6ea |
| 2026-10-22T07:10:35Z | 20261022T064617Z-lab810-kit-defo | terminated | 0.363 h | actual $0.392 | i-077ea1bc90813528a |
| 2026-10-22T07:10:44Z | 20261022T064851Z-lab811-kit-defo | terminated | 0.316 h | actual $0.341 | i-06e90d543fa2d2813 |
| 2026-10-22T07:10:52Z | 20261022T065148Z-lab810-kit-def | terminated | 0.275 h | actual $0.297 | i-0e772d2e90756fd8d |
| 2026-10-22T07:11:02Z | 20261022T065423Z-lab811-kit-def | terminated | 0.229 h | actual $0.247 | i-09d650382985d1d5f |
| 2026-10-22T07:19:55Z | 20261022T070218Z-lab801-kit-lz7 | terminated | 0.231 h | actual $0.493 | i-0cadbc0e1c15aa382 |
| 2026-10-22T07:24:39Z | 20261022T070608Z-lab801-kit-lz8 | terminated | 0.244 h | actual $0.523 | i-0ad1f2314bed9d6ea |
| 2026-10-22T07:25:03Z | 20261022T070218Z-lab801-kit-lz7 | terminated | 0.316 h | actual $0.676 | i-0cadbc0e1c15aa382 |
| 2026-10-22T07:25:43Z | 20261022T070608Z-lab801-kit-lz8 | terminated | 0.262 h | actual $0.561 | i-0ad1f2314bed9d6ea |
| 2026-10-22T07:29:28Z | 20261022T072630Z-lab811c-kit-lz7 | r7i.4xlarge on-demand $1.0584000000/h | cap 45 min | launch, cap $0.99 | i-0ea791651d3654c64 |
| 2026-10-22T07:32:23Z | 20261022T072928Z-lab811c-kit-lz8 | r7i.4xlarge on-demand $1.0584000000/h | cap 45 min | launch, cap $0.99 | i-0dc5ec7d280a877c7 |
| 2026-10-22T07:45:35Z | 20261022T074226Z-defbanks811 | r7i.4xlarge on-demand $1.0584000000/h | cap 45 min | launch, cap $0.99 | i-0cb33f8ee60f7b83a |
| 2026-10-22T07:45:50Z | 20261022T072630Z-lab811c-kit-lz7 | terminated | 0.273 h | actual $0.295 | i-0ea791651d3654c64 |
| 2026-10-22T07:45:59Z | 20261022T072928Z-lab811c-kit-lz8 | terminated | 0.227 h | actual $0.245 | i-0dc5ec7d280a877c7 |
| 2026-10-22T07:46:38Z | 20261022T074617Z-agglz811 | c7i.2xlarge on-demand $0.3570000000/h | cap 20 min | launch, cap $0.18 | i-04d6123952fcfa713 |
| 2026-10-22T07:47:55Z | 20261022T074617Z-agglz811 | terminated | 0.022 h | actual $0.008 | i-04d6123952fcfa713 |
| 2026-10-22T07:55:27Z | 20261022T074226Z-defbanks811 | terminated | 0.165 h | actual $0.178 | i-0cb33f8ee60f7b83a |
| 2026-10-22T07:58:31Z | 20261022T075804Z-kaabc811 | r7i.xlarge on-demand $0.2646000000/h | cap 20 min | launch, cap $0.13 | i-098a31a47b0a15e1b |
| 2026-10-22T08:02:44Z | 20261022T075850Z-kaabc811b | r7i.xlarge on-demand $0.2646000000/h | cap 20 min | launch, cap $0.14 | i-0e504e29af84158df |
| 2026-10-22T08:08:26Z | 20261022T075850Z-kaabc811b | terminated | 0.095 h | actual $0.027 | i-0e504e29af84158df |

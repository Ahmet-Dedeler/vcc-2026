# 534 — Generalization audit: does our climb solve the task, or fit A/B/C?

Session "Public data scope [2e83f3]", 2026-09-29, at Ahmet's request ("we don't want something random that increases but
doesn't actually solve the challenge"). Part A below is read from the written record only: `docs/BOARD-LEDGER.md`,
`data/board-entries/*.json`, and each arm's BOARD-RESULT / PREREGISTRATION / RESULTS. Nothing was recomputed. Part B
(the sealed mock final, `PREREGISTRATION-B.md`) is running as a worker and reports in `RESULTS-B.md`.

## 0. Answer

1. **Most of the climb is transfer, and it was tested on unseen data before the board read it.** Over the +0.0567 from
   `cd4svd5` (0.2160944, 09-05) to 446bj2z (0.2728260, 09-28):
   - **≈ 68–70 % is transfer-class** (+0.0386 … +0.0397): better predicted responses. Each step was priced on
     held-out lines (PLOLO cancer folds, kolf2 or H1) first. It landed on the board inside or above its preregistered
     band. The bulk switches were then re-priced on October-like panels.
   - **≈ 27 % is scorer mechanics** (+0.0151): the MDE, the zero-gene repair plus the re-level, and zs. These carry
     to D/E/F because the scorer does not change, and the one fitted part (the re-level) is re-derived from D/E/F's own
     controls on the day. They do not make the model better at biology.
   - **≈ 4–5 % was A/B/C-specific** (+0.002 … +0.003): the CD4 ×1.5 half of 417e. It won on A/B/C, lost on all three
     October-like panels, and October already ships CD4 ×1 (438, 457dr).
2. **The guard works in both directions.**
   - Twelve of the 30 scored tests since 09-04 read negative or flat and were never adopted (§3). In five of them, a
     local held-out instrument had predicted a gain and the board said no.
   - In the other direction, one A/B/C board gain (CD4 ×1.5) was caught by the October-like panels and excluded.
3. **Gaps found, and what was done about each** (§5):
   - No sealed test data was left: all four earlier reserves are spent. This session declared R6 and runs Part B on it.
   - The October decision record cites board reads, not the held-out reads. §4's table is handed to the coordinator
     for `docs/october-decisions.json`.
   - Three small bulk steps (normrel, cap, partial rows; +0.0037 together) have weak support off the board.
   - The DE channel has no cancer-line held-out surface. Its October switches rest on H1 plus the A/B/C board.
   - `AGENTS.md` §11.1c now makes the three-class label, the second reading and the sealed confirmation standing rules.

## 1. What counts as evidence here

- **Transfer** means a better prediction of responses.
  - Its evidence is a gain on data the model never saw: a PLOLO fold (a donor line held out and predicted from the
    others), the kolf2 fold, or h1-test without H1 donors.
  - The A/B/C board is itself held-out evidence of this kind (three real cancer lines no donor comes from). A step
    that was predicted locally and then read inside its band on the board has passed two independent tests.
- **Scorer mechanics** carry to D/E/F only as long as they are not fitted to A/B/C.
- **A/B/C-specific** means it helped A/B/C and failed off it.

## 2. The climb, step by step (every step is a scored board entry)

Δ is overall board units vs the previous best. "Before the board" is the held-out read the preregistration used.
"October-like" is 438's / 522's levelled k = 1 re-pricing on `real` / `efp` / `efr`, three cancer folds each.

| step | Δ overall | change | class | read on unseen data before the board | board vs that read | October-like panels | October recipe |
|---|---:|---|---|---|---|---|---|
| 312 | +0.0016 | per-target donor weights (normrel −1.5) | transfer | PLOLO Δpds +0.0026 [0.0000, +0.0053], 3/3 cancer folds; 097 replication | in band | not re-priced as a step (inside the PP base) | yes (base) |
| 319 | +0.0009 | cap form of the weights | transfer, small, **board-derived** | none: the hypothesis came from the board's own dose curve (312/312b). Afterwards, 323's H1 retrodiction got the sign right with an interval including 0 | small | not re-priced as a step | yes (base) |
| 319p | +0.0012 | CD4 partial-condition rows | transfer, small | 323 first; +0.0013 on a second base (254) | composes | not re-priced as a step | yes |
| 319p2 | +0.0003 | partial dose 0.5 | flat | — | flat bin | — | yes (dose) |
| 392 | +0.0025 | cross-line sign prior on the CPM block (CLEAN4) | transfer (signs) | H1; leave-one-bank-out (398) | fid top bin; nmae/reach did not follow H1 | n/a (DE) | yes |
| 395 | +0.0013 | keep the strongest 20 % of calls' own signs | transfer (sign delivery) | H1 optimum | in band | n/a (DE) | yes |
| 410b | +0.0126 | per-target magnitude key (raw mse −0.080), made free by ballast cells | transfer, delivered by mechanics | PLOLO, fold-held operators: aggregate cosine +0.27 … +0.30 on each held cancer fold (CIs above 0), kolf2 +0.50 | transfers; the welded form (410) read −0.0044 | the key permuted over targets loses on every panel: A/B/C −0.0223, real −0.0235, efp −0.0149, efr −0.0048 | yes (LODO form since 446bj2) |
| 417a | +0.0070 | bulk sharpening q 1.25 | transfer (shape) | the kolf2 fold predicted pds +0.0063 (board +0.0064); cancer folds over-read 2.5× | k_sh 0.42 / 0.65 | later q steps positive 3/3 on every panel (438 §3.4) | yes |
| 417e | +0.0051 | stack: H1 ×2 + CD4 ×1.5 | H1 half transfer; CD4 half A/B/C-specific | PLOLO +0.0037 (k = 1) | k ≈ 1.3 | H1 ×2: +0.0051 real (essential-inflated) / +0.0026 efp / +0.0003 efr. **CD4 ×1.5: −0.0016 / −0.0011 / −0.0024** | H1 ×2 yes; **CD4 ×1** |
| 433a | +0.0051 | coarse MDE on called genes | scorer mechanics | H1 nmae −0.0096 | 1.76× H1 | n/a | yes |
| 433h | +0.0095 | zero-gene repair (24 zeroed gate genes) + re-level to the mse cliff | mechanics: bug fix + calibration | repair predicted (read 3× its top bin); re-level exactly its price | as predicted | the re-level is re-derived from D/E/F's own controls; 453dr validated the October path on A/B/C inside its band | yes |
| 433k | +0.0046 | q 1.5 | transfer (shape) | 438 calibrated +0.0044 | hit | +0.0041 … +0.0049 (k = 1), 3/3 folds on each panel | yes |
| 446b | +0.0023 | q 1.75 | transfer (shape) | +0.0021 | hit | a further +0.0020 … +0.0030 (k = 1), 3/3 | yes (q = 1.75) |
| 446bj2 | +0.0022 | LODO per-target key | transfer | 438: +0.0027 on A/B/C | mse past its band | +0.0031 efp / +0.0032 real / +0.0017 efr | yes |
| 446bj2z | +0.0005 | zs exact-rank pass on carriers | scorer mechanics | H1 nmae −0.0013 … −0.0029 raw | hit (−0.0022; other five members bit-identical) | n/a | yes |

The steps sum exactly to 0.2728260 − 0.2160944 = +0.0567316.

**Attribution caveats.**
- 410b's gain is the key's mse gain (+0.014 overall-equivalent), with the ballast preventing the nmae cost that the
  welded 410 paid. Both are needed. It is booked as transfer because the information is the key's, and the permuted-key
  control shows that information is real on every panel.
- 417e's split uses 438's A/B/C prices (H1 +0.0024, CD4 +0.0017 at k = 1) and 457dr's board read of the CD4 half
  (+0.0031). Hence the ranges in §0.

**Before 09-05.** The climb from a15 (0.1741, 08-30) to `cd4svd5` came from transfer-class geometry and donor changes:
whitening, α, the CD4 residual, the H1 union (+0.0069, 12/12 PLOLO fold × dose cells positive, all six members rose),
and the CD4 deflation. Each was PLOLO-priced first (FINDINGS H.2, §8.x). Not re-audited here.

## 3. Tested and rejected: the guard at work

| entry | Δ vs its base | what it was | why it was not adopted |
|---|---:|---|---|
| ARM6 | −0.0036 | six union sources | gate test; sources saturate |
| 312b | −0.0013 | normrel −3.0 | dose curve turned over |
| 316 | −0.0008 | normrel parity (PM blend) | removed a distortion that was buying mse |
| 319d3 | −0.0003 | cap −3.0 | flat; dose family closed |
| 319pA | −0.0013 | Jurkat-calibrated rows for context A | H1 mispredicted the DE gains; pds paid |
| 387 | −0.0015 | consensus sign repair through the bulk | pds cost 2.3× H1's |
| 399 | −0.0016 | neighbourhood tilt | validated on H1 and leave-one-line-out, did not carry to the panel |
| 410 | −0.0044 | the key, welded | nmae cost; fixed by ballast (410b) |
| 417c | −0.0017 | unwhitened carriers | H1 mis-signed nmae |
| 446bja | −0.0011 | Jurkat programme as A's block prior | predicted fid +0.0015 … +0.0040; read +0.0002 |
| 446br2 | +0.0004 | zero-repair floor −2 | partial bin; −5 stays |
| 446bd2 | +0.0004 | DLD-1 as a second union residual | flat; stays out |

Two further entries (453dr, 457dr) validated the October path.
- In 5 of the 12 rejections (319pA, 399, 417c, 446bja, 446bd2), a local held-out instrument had predicted a gain and
  the board said no.
- In 387, the local read had itself predicted the loss, and the slot tested it. 312b and 319d3 were
  preregistered dose-curve points that mapped the far side of a peak.

**Our local instruments alone over-admit.** The pairing of a local held-out read with a preregistered A/B/C read is what
keeps false gains out.

## 4. The October switches: board read and held-out read (for `docs/october-decisions.json`)

| switch | value / status | A/B/C board read | held-out read | verdict under AGENTS 11.1c |
|---|---|---|---|---|
| `q` | 1.75, DECIDED | 433k +0.0046 (q 1.5), 446b +0.0023 (q 1.75), both in their bands | 438 §3.4: q 1.5 step +0.0041 … +0.0049 and q 1.75 step +0.0020 … +0.0030 (k = 1), 3/3 folds on real/efp/efr (A/B/C +0.0085 / +0.0046); the kolf2 fold predicted 417a's board pds | transfer; retention ≈ 0.5 |
| `lodo` | 1, DECIDED | 446bj2 +0.0022 | 438 §3.5: +0.0031 efp, +0.0032 real, +0.0017 efr (A/B/C +0.0027); the permuted-key control is negative everywhere | transfer; retention ≈ 1 |
| `lodo_hi_pct` | PENDING (10-03) | 495 = cap 6 vs 446bj2's cap 3 on A/B/C; PLOLO predicted +0.0019 … +0.0029 (522: the percentile cap equals cap 6 on A/B/C, both +2.62 % direction vs cap 3, 3/3) | 522, vs cap 3: +19.4 % efp, +11.6 % efr, +14.1 % real (3/3 each); the percentile cap vs cap 6, a further +5.95 % / +6.16 % / +5.41 % (CIs above 0) | two readings. A flat 10-03 read would be a failed prediction by the same instrument that prices the October value, so it counts against the switch and keeps cap 3 as preregistered (coordinator's correction, 2026-09-29; an earlier draft of this row said the A/B/C read could only show harm, which was wrong) |
| `sharpen_pt` | PENDING (10-02) | the 479 pair pending | PLOLO on the LODO stack: P2g +0.0035, S_P2 +0.0038 (3/3 cancer folds each); kolf2 +0.0020 / +0.0032; 9/9 October-like folds (493, 540c) | held-out positive; the board decides the pair |
| `cd4_stack` | ×1, DECIDED | 457dr: the CD4 ×1.5 half is +0.0031 on A/B/C | 438: CD4 ×1.5 −0.0011 efp, −0.0024 efr, −0.0016 real | **A/B/C-specific: the model case.** ×1 ships |
| `de_mag` | zs, DECIDED | 446bj2z nmae −0.0022 raw, other five bit-identical | H1: −0.0013 … −0.0029 predicted | scorer mechanics; A/B/C capacity is 0.70× H1's |
| `programme_source` | PENDING (09-30) | 458d vs 446b pending | H1 only, and H1 is the programme's own line (462e/f), so it says nothing about other lines. Public sign sources read ≈ 0 on every cancer fold (516), and lines share ~1/3 of a programme (462) | **the A/B/C read is the only cross-line test.** Treat a small positive as board-only in 11.1c's sense |
| `prior_gate` | PENDING (09-30) | 458g vs 458d pending | H1 only | as `programme_source` |
| `mde_sign` | PENDING (10-05) | 496a vs 446bj2 pending | H1, 3 seeds: nmae −0.00146 raw (CLEAN4 hook). 517's Arc hook is 12× that on H1, but it is H1's own programme | CLEAN4 is cross-line, so the rider is transfer-class. The Arc-hook form inherits `programme_source`'s caveat |
| `sameline_w` / `sameline_flags` | 0.5 / noprog, DECIDED | rows untestable on A/B/C; the programme failed for A (446bja −0.0011) | 450: a same-line bank from another lab on the H1 truth, +0.053 … +0.079 at 35 % coverage; a K562 replication | held-out only; applies only if a D/E/F line matches a registry bank (520: unlikely) |
| `max_levels`, `cd4_thresholds` | operational | — | 448/452 rehearsals; 481 F1 | n/a |

## 5. Gaps, and what was done

1. **No sealed data was left.** All four reserves (R3 to R5, Broad-16's four lines) are spent. Every surface in use has
   informed many choices (h1-test; the PLOLO folds; `real`/`efp`/`efr` through 438, 522 and 540–546). **Done:** R6 is
   declared in `docs/CONFIRMATION-LEDGER.md`: four fresh Arc-like panels, disjoint from every panel on disk. Part B
   prices the October bulk milestones on them once now (C10a) and once for the final recipe after 10-05 (C10b).
2. **The October record cites board reads.** `docs/october-decisions.json` stores each switch's board evidence; the
   held-out reads live in RESULTS files. **Done:** §4 is handed to the coordinator, who agreed to add a held-out line
   per switch. AGENTS 11.1c rule 2 makes it standing.
3. **Three small bulk steps have weak off-board support.** normrel −1.5 (+0.0016), the cap (+0.0009) and the partial
   rows (+0.0012):
   - normrel had a 3/3 PLOLO read before the board; the cap was chosen from the board's own dose curve; H1 got all
     their signs right, with intervals that include 0 (323);
   - they were never re-priced as separate steps on October-like panels. They sit inside every build's base.
   - Low stakes (+0.0037 together). **Part C measured it** (`RESULTS-C.md` §1: 438's levelled route, k = 1, pooled over
     real/efp/efr's 9 folds):
     - normrel HOLDS: +0.00112 [+0.00041, +0.00182], 8/9 folds.
     - The cap FAILS: −0.00011 [−0.00049, +0.00028], 3/9 folds.
     - The CD4 partial rows FAIL: −0.00026 [−0.00083, +0.00026], 4/9 folds.
     - Both FAILs are ≈ 0 off A/B/C, and both help on A/B/C's own panel: +0.00028, and +0.00094 [+0.00017, +0.00173].
       That `abc` row scores A/B/C's targets against the donor lines' truth. So the A/B/C-specificity is a panel
       effect, not an A/B/C-line effect.
     - Owner's recommendation to the coordinator: keep both for October. Across the four panels, keeping has the
       lower maximum regret (partial: 0.0003 vs 0.0009; cap: 0.0001 vs 0.0003). Both effects are below dayof438's 0.0010
       move threshold, and a `def_day` change carries engineering risk.
     - A default-off knob would make the switch readable on D/E/F feedback (FINAL-PHASE-FEEDBACK-PLAN) if members
       return.
4. **The DE channel has no cancer-line held-out surface.**
   - The pending DE switches (`programme_source`, `prior_gate`, `mde_sign`) rest on H1 plus the A/B/C board, and the
     programme sources are H1's own programme.
   - This is also where the leaders beat us (reach, nmae).
   - **Open, for the coordinator:** a cancer-line DE check would need cell-level truth from a line that is neither a
     donor nor a spent reserve, which is not on disk today.
   - Until then, rule 2 treats a DE switch whose only cross-line read is A/B/C as board-only.
5. **The scorer-mechanics share (~27 %) is D/E/F-dependent.** The MDE's value scales with sign accuracy, the re-level
   with depth. Both are re-derived or re-read on the day (the re-level from D/E/F's controls; the ρ* re-read in
   FINAL-PHASE-FEEDBACK-PLAN day 2 if members return). Nothing to do before Oct 22.

## 6. What this does not show

- **It does not estimate the D/E/F score.** Retention of the transfer share on new lines is measured only for the bulk
  switches, and only on panels drawn from the donor lines' own coverage. Part B adds a sealed version of that
  measurement. A line-level estimate needs data we do not have.
- **Every number here is quoted from the record.** 410b's and 417e's splits are approximate (§2). No instrument was
  re-run; the known-answer discipline for new numbers is Part B's.

# Lab instrument report: the fresh-panel over-read (written 2026-10-16 11:14Z by the clock)

For the coordinator's four questions (message of 10-16 ~10:45Z) before the 10-20 23:59Z freeze. **An instrument report, not a
dossier.** Every number is in board units (overall) unless marked r (relative aggregate cosine). Files:
`experiments/800-checker/work/panelrank800.txt`, `overread800-{LABC,ABCp,ABC}.txt`, `experiments/801-lab-folds/RESULTS.md` §14,
`experiments/810-strong-direction/RESULTS.md` §11, `experiments/811-lolo-model/RESULTS.md` §0000.

## 1. Panel or rules? **The panel.**
- **LABC** = A/B/C's 300 targets built with the exact fresh-fold recipe: dump801 (a copy of dump708, the def_day fold builder),
  D/E/F rules, day-recipe operators, the same held truth lines (K562 GWPS, Orion HCT116, Orion HEK293T), and inc = pct+P2g.
  801 built it once, for $0.17.
- With every fold on the SAME recipe, A/B/C's panel ranks at the **bottom of its held line's seven folds in 14/21 (arm, held)
  cells, where 3.0 are expected; mean percentile 0.19 (0.50 if exchangeable)**. This is identical to ABCp, the A/B/C rules
  with the day operators.
- Key arm 495 (cap6 − cap3), as r:

  | held line | LABC | ABCp | L1–L6 |
  |---|---|---|---|
  | K562 | +2.0 % | +1.5 % | +8.4 … +14.0 % |
  | HCT116 | +3.7 % | +3.4 % | +5.0 … +9.2 % |
  | HEK293T | +4.2 % | +3.7 % | +7.4 … +13.6 % |

  A/B/C's fold rules account for only +0.4–0.5 points of r.
- **Not the cause:**
  - essential genes: 0/300 on A/B/C, and 0 on every L panel;
  - the held-line truth-strength profile: identical truth, inside the L spread, weak and inconsistent correlation with the gain.
- 801's inference: the difference sits on the source side, where A/B/C's LODO key has a light tail. 814 (opened today) is
  finding the target property and a draw rule v2.

## 2. The day-of folds on an Arc-drawn panel: **they behave like A/B/C, not like the fresh draws.**
- LABC IS that proxy: the def_day fold recipe applied to the one Arc-drawn panel we hold. On it, key and direction gains
  shrink as on A/B/C's own folds (bottom 14/21, as above).
- So, on the assumption that Arc drew D/E/F's panel the way it drew A/B/C's, 708's day-of folds on D/E/F's own panel should
  NOT carry the fresh-panel over-read. The stack weights decided there are immune in that sense.
- The over-read belongs to panels drawn by a rule other than Arc's (the lab's arc800, and plausibly production's
  October-like panels).

## 3. Per-class over-read factor (fresh arc800 panels ÷ Arc-drawn panel)
Two-level bootstrap over the six panel draws and over targets; 4,000 draws; pooled over arms × held lines as a ratio of sums.

| class | arms | vs LABC (like-for-like) | vs ABCp | lower bound (fresh CI low ÷ A/B/C CI high) |
|---|---|---|---|---|
| **key** | 495 cap6−cap3, 705 θ | **3.08 [1.91, 6.76]** | 3.78 [2.15, 10.9] | 1.6 |
| **direction** | 707, sgu1, 811 c0_l25 | **6.5, upper bound open** (the A/B/C-panel sum's CI reaches 0: [−0.0033, +0.0243]) | 4.72 [1.88, 27.5] | **2.0** |
| **shape** | 479p2g, 479sp2 | **1.19 [0.81, 1.81]**: no deflation | 1.16 [0.80, 1.74] | 0.7 |

- Board cross-check: the checker's calibrated key ratio W ÷ A/B/C = 4.75 ÷ 1.46 = 3.25, and shape 3.10 ÷ 3.48 = 0.89.
- **Applying them to the October-like prices.** These are the lab's factors for its arc800 panels. Production's efp, efr, real
  and R6 panels were drawn by other rules, but R6's percentile cap reading 3.31× its A/B/C value says the same bias is likely.
  - **lodo_hi_pct** (+5.4 … +6.2 % r off-board): key ÷ 3.08 gives **≈ +1.8–2.0 %** (range +0.8 … +3.2 %). It agrees in sign
    with its A/B/C board read: 495, cap 6 vs cap 3, +0.00163.
  - **lodo** (+0.0017 … +0.0032): key ÷ 3.08 gives ≈ +0.0006 … +0.0010 (range +0.0003 … +0.0017). It has its own board read,
    +0.0022 (446bj2), which stands.
  - **705 θ** (+0.7 … +1.9 % r): key ÷ 3.08 gives ≈ +0.2 … +0.6 %. On the Arc-drawn panels 705 itself reads ≈ 0:
    LABC +0.00050 / +0.00007 / −0.00016 by held line; ABCp −0.00005 [−0.0004, +0.0003].
  - **707**: direction ÷ ≥ 2.0. On the Arc-drawn panels it reads LABC +0.00115 [−0.0001, +0.0025] 3/3 and ABCp +0.00178
    [+0.0005, +0.0032] 3/3, against +0.0048 on the fresh panels.

## 4. Sign disagreements (fresh-panel price vs A/B/C board read): **none among the DECIDED switches.**
- lodo: board +0.0022, held-out +. Agrees.
- lodo_hi_pct: board 495 +0.00163, held-out +. Agrees.
- sharpen_pt: board 479p2g +0.00104, held-out +. Agrees.
- The one known disagreement, 479sp2 (S_P2: fresh + and Arc-drawn +, board −0.00071), was not adopted. Its miss is not a panel
  effect; every lab surface misses it.
- Pending: 705 θ (OFF; 10-17 s1). Fresh +, Arc-drawn ≈ 0, so a board read ≤ 0 would be a disagreement, but it changes nothing:
  the switch is off.
- Pending: 707 (10-17 s2 / 10-18 s1). Arc-drawn panels are positive, so the same sign is predicted.

## Tries and cost
- panelrank800 and overread800: one version each, no tuning.
- Comparators: ABC062, ABCp, LABC.
- Arms: 7, all read from arrays that existed before this question, plus one LABC read per arm.
- Wear: 810 sgu1 N 8, 707 N 3, sgu1all N 3; 811 c0_l25 N 11. Rows merged into `800-checker/wear/ledger.jsonl`.
- Lab cost of the whole panel study ≈ $0.6, all on EC2, with no local job during the coordinator's hold.

## Side result (811)
- On LABC, 811's **many-line increment over sources-only holds: +0.0036 [+0.0018, +0.0055] 3/3**, only 1.7× smaller than on the
  fresh panels.
- The sources-only blend turns negative on A/B/C's panel (LABC −0.0014, 062 −0.0022).
- 811's full arm stays WEAK (LABC +0.0023 3/3; after e(11) +0.00006). Still no dossier.

## Addendum (2026-10-16 13:21Z by the clock): what the essential-screen cut does and does not explain
- **814 (11:41Z):** A/B/C's panel holds 0 targets of the RPE1 ∪ Jurkat essential screens, against 16.9 ± 3.7 per arc800 draw.
  That is Arc's stated criterion (FINDINGS §8.1g), which the lab's arc800 rule broke. Dropping those targets from L1–L6 moves
  A/B/C's mean percentile from 0.19 to 0.36.
- **The preregistered test FAILS on new draws.** 801 drew six panels with the screens removed (A1–A6, 18 folds, same recipe).
  Among them A/B/C reaches only **0.26, bottom 9/21** (bar ≥ 0.35 and ≤ 6/21).
  - sgu1 +0.0087 and 707 +0.0047 on the new panels (17/18 each), as on L.
  - The 811 consensus arm +0.0064, 18/18.
- **Factors, A-only fresh ÷ LABC:**
  - key **2.53 [1.53, 5.43]**;
  - direction **≈ 6.0, upper bound open**;
  - shape **0.99 [0.68, 1.51]**.
- **So the cut explains part of the key over-read and little of the direction over-read.** Re-reading a panel with the
  essential-screen targets dropped is NOT a sufficient deflation. Use the LABC factors and decide key / direction on Arc-drawn
  panels.
- **Open question:** whether the rest is a property of Arc's selection or the idiosyncrasy of one draw. D/E/F's panel (the
  second Arc draw) is the evidence. 814's predictions for it are frozen: P-A (≤ 3 essential-screen targets) and P-B (prop10
  deciles).

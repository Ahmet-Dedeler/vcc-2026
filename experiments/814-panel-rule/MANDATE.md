# 814 — Why is Arc's panel not like the lab's panel draws? A draw rule v2 (instrument lane)

Lab lead's mandate, 2026-10-16 10:51Z. Read `docs/LAB-2026-10-16.md` (charter + Amendment 1) and
`experiments/800-checker/PHASEB.md` first. Then read the LAB-LOG notes of 10:41Z and 10:46Z, and
`experiments/800-checker/work/panelrank800.txt` and `overread800-ABCp.txt`.

## Objective
The lab's fresh folds (L1–L6, 801) were drawn by the "arc800" rule (`experiments/800-checker/DESIGN.md` §3.1,
`experiments/801-lab-folds/draw801.py`). On the same held truth lines, the same arms and the same operator recipe, A/B/C's real
panel sits at the bottom of the six draws for key and direction changes:
- 14/21 cells at the bottom; mean percentile 0.19;
- over-read factors 3.8× (key) and 4.7× (direction);
- shape changes are exchangeable (1.16×).
D/E/F's panel will be drawn by Arc the way A/B/C's was, so a fold set that over-reads key and direction changes misprices every
candidate of those classes for the final.

Your job:
1. Find which target-level properties separate A/B/C's 300 targets from arc800 draws.
2. Propose a draw rule v2 under which A/B/C's panel is typical.

**This lane is instrument work only. It prices no candidate.**

## Already known (do not redo)
- Essential genes are NOT the difference. A/B/C's panel holds 0 of the 2,021 K562-essential perts
  (`data/donors/k562_essential/perts.txt`), and so does every L panel (lead, 10:51Z).
- arc800 matches donor-coverage pattern and a strength mix binned by the **GWCD4i Rest downstream-DE count**, and applies a
  basal-expression floor. Arc binned by an **H1 pilot screen** (AGENTS §2.1; FINDINGS H.1, l.911).
- Separately, 801 (follow-up 4) is building LABC (A/B/C's panel on the L-fold recipe) and a per-FOLD truth-strength profile
  (median d, top-30 share). Do not build folds; LABC and the profile are 801's.

## Steps
1. **Target features.** For A/B/C's 300, for the arc800 pool, and for each of L1–L6, compute what the data on disk allow, each
   with its source path and sha256:
   - truth effect strength per held line (K562 GWPS, Orion HCT116, Orion HEK293T): ‖δ‖ and the number of DE genes where
     available;
   - VCC-2025 H1 effect strength where the target is covered (`data/donors/h1_2025`; any H1 DE count you can derive). This is
     the closest public proxy to Arc's H1 pilot bins;
   - GWCD4i DE count; basal expression in A/B/C's controls and in each held line; knockdown efficiency where recorded;
   - donor coverage count; the source-agreement statistic z_src if a table exists (707 / 522 used it).
2. **Where A/B/C sits.** Compare A/B/C's panel against 1,000 arc800-rule redraws from the same pool (draw801's code, new seeds
   ≥ 814000, never the LZ seeds 800107/800108), giving A/B/C's percentile on each feature. Name the 1–3 features where it is
   extreme.
   - Read the per-target values of the five candidate arms on folds already on disk (the arrays panelrank800 reads). Check
     whether those features explain which targets carry the key and direction gains on L folds but not on A/B/C's panel.
3. **PROPOSAL.md** for draw rule v2: the matched features and bins, with a simulation showing A/B/C's panel is typical under v2
   (feature percentiles 0.1–0.9).
   - Plan the validation: new v2 draws (new seeds), day-recipe fits, dumps by 801, then panelrank800 against them. The bar,
     written now: A/B/C's mean percentile ≥ 0.35 and bottom-count ≤ 6/21.
   - Give cost, class (instrument) and try count. The judge reviews it before any fold is built.

## Rules
- Write only in `experiments/814-panel-rule/`. No git. No vcc. Nothing in data/.
- Local compute only when `bash experiments/800-checker/window800.sh` prints `CLEAR`. Right now it prints "WINDOW hold (local
  jobs only)": the coordinator's hold runs until R7a ends or 18:15Z. Until then, compute on EC2 via
  `experiments/800-checker/cloud/labcloud.sh` (on-demand, Project=lab). Lane budget ≤ $5 before the judge. File reads are
  fine any time.
- Never use or read the sealed reserves (LZ L7/L8 targets; anything in 551/552) or production's R6/R7 panels beyond their
  target lists' exclusion in the pool, which draw801 already applies.
- Report in `RESULTS.md` (features table, A/B/C percentiles, what explains the gains). Then write PROPOSAL.md and stop.

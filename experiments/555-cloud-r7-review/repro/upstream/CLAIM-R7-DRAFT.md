# R7 claim — DRAFT (not the claim: nothing here is recorded in the ledger, no R7 number exists)

Drafted 2026-10-02 ~01:40Z by the reserve's owner (session "Explore contribution areas"), before C10b (10-04) and before the
DE reads of 10-05 … 10-07. The rules below are meant to be frozen as written; only the configurations marked **FILL** are
taken from `docs/october-decisions.json` when the DE switches settle. Then this file is copied to `CLAIM-R7.md`, its sha256
goes into the ledger's R7 section, and only after that does any R7 step run (RESERVE.md §4).

## 1. What R7 confirms

One claim, once. **The final October recipe against the recipe before the switches decided after C10a**: its bulk
configuration on R7a (four Arc-like panels, fresh in targets), its DE configuration on R7b (Calu-3, floored). R7 confirms; it
does not choose. No fitting, selection, dose search or debugging may read an R7 score, and nothing in this claim is
re-run with other settings after a number exists.

## 2. R7a — the bulk channel

- **F = FILL** (a `designs551.CONFIGS` key: `M5`, `M6`, `M5+P2g`, `M6+P2g`, `M5+S_P2`, `M6+S_P2`), from `lodo_hi_pct` (10-03)
  and `sharpen_pt` (10-02: `tercile:2.5,2.0,1.75`, upgrade to `tercilec` queued as 598). designs551's P2g is the carried-only
  form (cut points from the carried targets' key weights, non-carried MID, every floor target LOW: 493h `P2g_carried`), the
  form October ships once 598 is applied. Read side by side on 10-02 (598's patch, the builder's `tercile` tie rules,
  lib551.q_p2g): the same carried-only cut points, floor guard, cap-tie and non-carried rules; the numerical tie (tiers on
  efp) belongs to 598's known answers. If October keeps the plain `tercile`, this claim prices the carried-only form and
  says so.
- **Instrument:** 551's code, unchanged: the 28 repo files its processes import are byte-identical to their blobs at
  `e7be21ca`, where the efp known answers passed (|d| 0.0). `r7a551.py plan` proves this; `record` pins it, and `access`
  refuses to run if it changed. New code is limited to the driver (`r7a551.py`) and the aggregation
  (`aggregate_r7a551.py`), both pinned at record. KA-CD4 (efp) passed 10-02: the CD4 step reproduces 238's efp banks (38
  files byte-identical, 4 json records equal but for their `peak_rss_gib`).
- **Procedure (per panel, T1 … T4 in order):** CD4 banks (gates on, then relaxed, as 534) → spec.json → four operator fits,
  fresh (534's commands without an expected sha; each sha256 and axis recorded) → axis check (each axis is its own panel's;
  not A/B/C's, not a tuned-on panel's, not R6's) → `full551 T --config F` → `comp551 T <3 folds> --level --config F` →
  operators deleted (sha256 and bytes kept). Then `aggregate_r7a551.py`.
- **PRIMARY:** F against M5, the sum of F's post-C10a steps (M5 → M6 if F's base is M6; F's per-target step on its base, if
  any). If F = M5, the primary is the cumulative M0 → M5. 438's levelled route at k = 1; the mean over 12 folds; paired
  target bootstrap, 4,000 draws, seed 552. **Verdict:** TRANSFERS if D > 0 and the 95 % lower bound > 0; FRAGILE if D > 0
  and lb ≤ 0; FAILS if D ≤ 0.
- **Reported, not judged:** each post-C10a step alone; the per-target step's control (q permuted over targets, seed 493) and
  the targeting value (step − control, paired); the cumulative M0 → F with its retention against R6 (C10a's pooled
  milestones, C10b's step) and against A/B/C (534's record; F's step on A/B/C's own panel if 551's abc run exists); the
  calibrated board-unit projection (designs551.K_CAL).
- **Prior prediction (fixed now):** M5 → M6 on R6 (C10a) **+0.00694** [+0.00517, +0.00828], 10/12 folds (A/B/C record
  +0.0021); P2g's step on R6 is C10b's reading (10-04), with the PLOLO prior +0.0035 (534 §4, 3/3 cancer folds). R7a uses the
  same donors on fresh targets, so the expected retention against R6 is ~1. If F = M6 + P2g, F vs M5 is predicted at about
  +0.0069 + C10b's P2g value, TRANSFERS. FAILS would contradict R6 and the PLOLO panels.
- **Consequence:** TRANSFERS: nothing changes. FRAGILE: flagged to the coordinator. FAILS: the post-C10a switches go to the
  coordinator as candidates to switch off (the coordinator decides; R7 does not).
- **Cost and window:** 16 fits (~5 min, ≤ 14 GiB each) plus full551/comp551 (≤ 8 GB), ~3–4 h serially, at C10b's floors
  (22 GiB free+inactive + pressure 55 % + swap 400 MiB), only inside the coordinator's recorded window.

## 3. R7b — the DE channel

- **Surface:** `r7bd5000` (552's floored truth, amendment 3: 106 targets, 2,378 controls; BST2 excluded).
- **Arms (FILL from the decided DE switches, in 553's MODE grammar):** the base `cpm` (the recipe before the post-C10a DE
  switches); the final recipe's MODE; for every decided DE switch, the final recipe with that switch off (so each switch is
  priced as a switch-on-recipe pair); and the retrodiction arm `lcap-kept` (the call-propensity block alone, already read
  on the board). Map: CP = `lcap-kept`, hook = `+c4h`, CPH = `lcap-kept+c4h`, K 2,000 = `-n1600`. Any arm not yet on disk
  is emitted before the access (`run553.py r7bd5000 <MODE> base <seed> … --exclude-targets BST2`, seeds 0–2, and seed 0 on
  each of the 10 control halves; emitting reads no truth).
- **Statistic:** per arm (vs `cpm`) and per pair, the overall-equivalent effect at R7b's own anchors, seeds 0–2, with
  `band_total` = sqrt(band_sd² + sigma_c²) (seeds + split-half truths + the control re-draw from the 10 halves, 553 §"the
  band").
- **Board-likeness first (552's requirement, rule 3):** CP's raw member effects on R7b next to the shape-matched rehearsal's
  mean (`--retro work/access553_rehearsal5`). The board read CP (10-01) with jac + and fid −. R7b is **board-like** if its CP
  reading has raw jac > 0 and raw fid < 0. Otherwise every claim row is reported as "instrument not board-like" and none of
  them can veto.
- **HARM veto (amendment 1, asymmetric):** a decided DE switch is a HARM candidate only if its switch-on-recipe effect is
  **≤ −2.5 × band_total** of that pair AND R7b is board-like. The 2.5 allows for the between-draw spread, which on the
  rehearsal ran ~1.3× band_total (draws553: 0.83–3.0, 90 %), and keeps the false-veto rate per switch ≈ 1–3 %. A HARM goes
  to the coordinator as a veto candidate. A null or small effect never counts against a switch; a positive one is
  robustness evidence only.
- **Prior prediction (fixed now; the five-draw means of the shape-matched rehearsal, `553/work/draws553.json`, overall-
  equivalent):** CP +0.0014 (between-draw sd 0.0019), CPH +0.0015 (0.0019), hook ~0 (0.0004), K 2,000 −0.0016 (0.0016),
  K on CPH −0.0017 (0.0023), CPHK20 −0.0001 (0.0031). Expected: no HARM bin for CP, CPH or the hook. A K-size harm reads at
  about −1 to −2 band, so a HARM there needs Calu-3 to show it more strongly than H1 did.
- **Command (one access, ~25 min, < 4 GB, the truth opened only here):**
  `access553.py r7bd5000 --owner-go CLAIM-R7.md --exclude-targets BST2 --arms <FILL> --seeds 0,1,2 --pairs <FILL>
  --ctrl-splits 5 --retro work/access553_rehearsal5 --halves 0,1,2,3,4` (after its `--dry-run`).
- **Void** for R7b if an Oct 22 context matches Calu-3 before the access (the ledger's rule).

## 4. Order of operations

1. Fill §2's F and §3's arms and pairs from `docs/october-decisions.json` (quote the decided values and their timestamps).
   Freeze the file as `CLAIM-R7.md`; record its sha256 in the ledger's R7 section and commit, before any R7 step.
2. R7a: `r7a551.py plan`, then `r7a551.py record --claim experiments/552-r7-reserve/CLAIM-R7.md --config F --window-start …
   --window-end … --received … --message-file …`, then `r7a551.py access` inside the window.
3. R7b: emit any missing arm (no truth), `--dry-run`, then the access with `--owner-go`.
4. Append both outcomes to the ledger; report to the coordinator (numbers, verdicts, and the class of each change).

## 5. Count of tries (AGENTS §11.1c.4)

R7a: one access with F fixed before it. R7b: one access with its arms and pairs fixed before it. No R7 number has been
read by anyone before this claim is recorded.

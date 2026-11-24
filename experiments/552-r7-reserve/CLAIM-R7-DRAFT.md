# R7 claim — DRAFT (not the claim: nothing here is recorded in the ledger, no R7 number exists)

Drafted 2026-10-02 ~01:40Z by the reserve's owner (session "Explore contribution areas"), before C10b (10-04) and before the
DE reads of 10-05 … 10-07. Revised 02:30Z with every finding of 555's adversarial review (Grok 4.7 xhigh;
`experiments/555-cloud-r7-review/RESULTS.md`, dispositions in `RESPONSE.md` there). The rules below are meant to be frozen as
written; only the configurations marked **FILL** are taken from `docs/october-decisions.json` when the DE switches settle.
Then this file is copied to `CLAIM-R7.md`, its sha256 goes into the ledger's `## Reserve R7` section, and only after that
does any R7 step run (RESERVE.md §4).

## 1. What R7 confirms

One claim, once. **The final October recipe against the recipe before the switches decided after C10a**: its bulk
configuration on R7a (four Arc-like panels, fresh in targets), its DE configuration on R7b (Calu-3, floored). R7 confirms; it
does not choose. No fitting, selection, dose search or debugging may read an R7 score, and nothing in this claim is
re-run with other settings after a number exists. Only the readings named here count: a later re-report (new `--pairs`,
`--retro`, `--nb`, panels, or a different band) is exploratory, is labelled so, and can never confirm or veto anything.

## 2. R7a — the bulk channel

R7A_F = FILL

- **F** is the line above: one `designs551.CONFIGS` key (`M5`, `M6`, `M5+P2g`, `M6+P2g`, `M5+S_P2`, `M6+S_P2`), from
  `lodo_hi_pct` (10-03) and `sharpen_pt` (`tercilec:2.5,2.0,1.75` once 598 lands on 10-06). `r7a551.py record` reads that
  line and refuses a `--config` that differs from it, and refuses `FILL`.
- **The instrument stays at `e7be21ca` and F is one of the six keys.** designs551's P2g is the carried-only form (cut
  points from the carried targets' key weights, non-carried MID, every floor target LOW: 493h `P2g_carried`), which is
  598's `tercilec` (same cut points, floor guard, cap-tie and non-carried rules, read side by side on 10-02; 598's A/B/C
  known answer moves exactly 493h's 4 targets). designs551 has no plain-`tercile` design. If October ships anything F's key
  cannot express, that needs a new pin and a new efp known answer BEFORE `record`. It is never a note added after the score.
- **Instrument pin:** 551's code, unchanged. The 30 repo files its processes can import are byte-identical to their blobs
  at `e7be21ca`, where the efp known answers passed (|d| 0.0). The pin is the import closure over every import statement,
  including function-level ones (lib551's `ntclevel417`, which the T-panel claims use). The package versions are recorded,
  and the environment is the efp KA's (uv.lock and the venv unchanged since 08-21). `plan` proves all of this, `record` pins
  it, and `access` refuses on any change. New code is limited to the driver (`r7a551.py`) and the aggregation
  (`aggregate_r7a551.py`), both pinned at record. KA-CD4 (efp) passed 10-02: the CD4 step reproduces 238's efp banks (38
  files byte-identical, 4 json records equal but for their `peak_rss_gib`).
- **If F = M5** (no post-C10a bulk switch), R7a is **not accessed**: `r7a551.py access` says so and stops before any T
  panel step, and the reserve stays sealed for a later claim (C10b's rule for C10a designs).
- **Procedure (per panel, T1 … T4 in order):** CD4 banks (gates on, then relaxed, as 534) → spec.json → four operator fits,
  fresh (534's commands without an expected sha; each sha256 and axis recorded) → axis check (each axis is its own panel's;
  not A/B/C's, not a tuned-on panel's, not R6's) → `full551 T --config F` → `comp551 T <3 folds> --level --config F` (each
  fold's json must carry `config == F`) → operators deleted (sha256 and bytes kept; resumable). Then
  `aggregate_r7a551.py`, once (it refuses a second write, and any `--nb` or panel list other than the preregistered ones).
  Heavy steps wait for 551's gate at C10b's floors (22 GiB free+inactive, also with any `GATE-OVERRIDE.json` present;
  pressure 55 %; swap 400 MiB), inside the coordinator's window only.
- **PRIMARY:** F against M5, the sum of F's post-C10a steps (M5 → M6 if F's base is M6; F's per-target step on its base,
  if any). 438's levelled route at k = 1; the mean over 12 folds; paired target bootstrap, 4,000 draws, seed 552 (the plain
  and the calibrated rows share the draws). **Verdict:** TRANSFERS if D > 0 and the 95 % lower bound > 0; FRAGILE if D > 0
  and lb ≤ 0; FAILS if D ≤ 0.
- **Reported, not judged:** each post-C10a step alone; the per-target step's control (q permuted over targets, seed 493) and
  the targeting value (step − control, paired); the cumulative M0 → F; its retention against R6 (C10a's pooled milestones,
  and C10b's step only if C10b's recorded F equals this F; otherwise no R6 retention is printed); retention against A/B/C,
  whose M3 → M6 denominators are record534's text values from other routes and are labelled so by the aggregator; the
  calibrated board-unit projection (designs551.K_CAL).
- **Prior prediction (fixed now):** M5 → M6 on R6 (C10a) **+0.00694** [+0.00517, +0.00828], 10/12 folds (A/B/C record
  +0.0021); P2g's step on R6 is C10b's reading (10-04), with the PLOLO prior +0.0035 (534 §4, 3/3 cancer folds). R7a uses the
  same donors on fresh targets, so the expected retention against R6 is ~1. If F = M6 + P2g, F vs M5 is predicted at about
  +0.0069 + C10b's P2g value, TRANSFERS. FAILS would contradict R6 and the PLOLO panels.
- **Consequence:** TRANSFERS: nothing changes. FRAGILE: flagged to the coordinator. FAILS: the post-C10a switches go to the
  coordinator as candidates to switch off (the coordinator decides; R7 does not).
- **Cost and window:** 16 fits (~5 min, ≤ 14 GiB each) plus full551/comp551 (≤ 8 GB), ~3–4 h serially; window 10-07
  06:00–12:00Z (tentatively agreed with the coordinator 10-02), R7b first.

## 3. R7b — the DE channel

- **Surface:** `r7bd5000` (552's floored truth, amendment 3: 106 targets, 2,378 controls; BST2 excluded).
- **Arms (FILL, in 553's MODE grammar):** the base `cpm` (the recipe before the post-C10a DE switches, zs included); the
  final recipe's MODE; for every decided DE switch, the final recipe with that switch off (each switch is then a
  switch-on-recipe pair); and the retrodiction arm `lcap-kept`. Map: CP = `lcap-kept`, hook = `+c4h`, CPH =
  `lcap-kept+c4h`, K 2,000 = `-n1600`, K 1,600 = `-n1200`. Every arm the 10-07 tree can produce is emitted already
  (coordinator, 10-02).
  - ARMS = FILL
  - PAIRS = FILL (every switch-on-recipe pair)
  - HALF_ARMS / HALF_PAIRS = FILL: the veto-eligible arms and pairs (below). They must be on the 10 control halves (the six
    core arms are); the access lists every other effect as "no band_total" in plan.json before the truth opens.
- **The pairs follow from the decided switches by this rule (fixed now, filled mechanically by `fill_claim552.py`).**
  B = the call-propensity block, H = the CLEAN4 hook, K = the call-set size (default, 2,000 or 1,600).
  MODE(B, H, K) = (`lcap-kept` if B else `cpm`) + (`-n1600` for K 2,000, `-n1200` for K 1,600) + (`+c4h` if H).
  The final recipe's arm is MODE(B, H, K). The block and the hook are judged at the default call-set size, on
  core = MODE(B, H, default): for each of B and H that is on, the pair core vs core-without-it. These are the veto-eligible
  pairs, and their arms are among the six on the control halves. If K is not the default, K's pair is final vs core,
  report-only. When B and H are both on, the hook's sign-shuffle control (`lcap-kept+c4h-s` vs `lcap-kept+c4h`) is
  reported. `lcap-kept` (CP vs cpm, the board-likeness reading) is always an arm. The tree's recipes and their pairs:

  | decided | final arm | veto-eligible pairs (on the halves) | report-only pairs |
  |---|---|---|---|
  | B | `lcap-kept` | `lcap-kept:cpm` | — |
  | H | `cpm+c4h` | `cpm+c4h:cpm` | — |
  | B + H | `lcap-kept+c4h` | `lcap-kept+c4h:cpm+c4h` (block), `lcap-kept+c4h:lcap-kept` (hook) | `lcap-kept+c4h-s:lcap-kept+c4h` |
  | B + H + K 2,000 | `lcap-kept-n1600+c4h` | the two of B + H, at the default size | `lcap-kept-n1600+c4h:lcap-kept+c4h` (K), the sign control |
  | K only (2,000 or 1,600), or nothing decided | — | — | R7b is **not accessed**: no veto-eligible switch, so the reading could not confirm or veto anything (the reserve stays sealed) |
- **Statistic:** per arm (vs `cpm`) and per pair, the overall-equivalent effect at R7b's own anchors, seeds 0–2, with
  `band_total` = sqrt(band_sd² + sigma_c²) from band553 (seeds + split-half truths + the control re-draw of the 10 halves;
  553 "the band"). `band_sd_ctrl` (the `--ctrl-splits` band) is reported and is not the band of any rule.
- **Board-likeness first (552's requirement, AGENTS rule 3):** R7b is **board-like** if its CP reading (`lcap-kept` vs `cpm`,
  the raw member effect means, never the retro ratio) has raw jac > 0 and raw fid < +0.001. The board read CP (10-01) with
  jac + and fid −. The +0.001 allowance is one rehearsal fid sd: all five rehearsal draws pass it, while a strict fid < 0
  failed one (h1r7s2, +0.00006; 555 finding 9). `--retro work/access553_rehearsal5` prints each member against the
  rehearsal's five-draw mean (CP raw jac +0.0041, fid −0.0018). If R7b is not board-like, every row is reported as
  "instrument not board-like" and none can veto.
- **HARM veto (amendment 1, asymmetric), for the veto-eligible switches only: the call-propensity block, the CLEAN4 hook,
  and their combination.** A decided switch of these is a HARM candidate only if its switch-on-recipe effect is **≤ −2.5 ×
  that pair's band_total**, the pair is on the control halves, AND R7b is board-like. A HARM goes to the coordinator as a
  veto candidate. A null or small effect never counts against a switch, and a positive one is robustness evidence only.
  - The false-veto rate per switch is 0.6 % if the surface's noise equals band_total. The rehearsal's between-draw spread
    is wider (between-draw sd / band_total on h1r7: CP 1.57, CPH 1.44, hook 0.42, hook-on-CP 0.78). At that spread the null
    tail is CP 5.5 %, CPH 4.1 %, and ~9.5 % for the three together if independent (555 finding 8, `calibrate_harm.py`).
- **Call-set size switches (K 2,000, K 1,600, K on any recipe) are report-only: no veto.** On one R7b-shaped surface K's
  between-draw spread is 2.5× its band_total, and its rehearsal prior is already negative (K on CPH −0.00165, between-draw
  sd 0.0023). H1-shaped data price K below zero while the board decides it. A zero-centred line would fire ~40 % of the time
  by reproducing that known price. A prior-centred line (−0.00395) has power ~0.16 against a K-size harm. One Calu-3 reading
  cannot resolve a K-sized effect, so the claim says so and does not pretend the 2.5 line will (555 finding 8).
- **Identity-type sign sources (604 / 622 / 639: the identity-axis sign prior and its variants) are REPORT-ONLY on R7b,
  never a veto** (the coordinator, 2026-10-02). R7b's truth (Sunshine Calu-3) is 3′ data, while the rider's references
  are Flex (A/B/C). On 3′-only references the identity score reads −0.05 … −0.09 on HEK293T (622/639), so a Calu-3
  reading would mix chemistry into the signal. The fill rule's tree has no such switch, and 553 has no arm for it. If one
  is decided before the access, the claim says in writing, before it is frozen, that R7b does not read it.
- **Prior prediction (fixed now; the five-draw means of the shape-matched rehearsal, `553/work/draws553.json`, overall-
  equivalent):** CP +0.0014 (between-draw sd 0.0019), CPH +0.0015 (0.0019), hook ~0 (0.0004), K 2,000 −0.0016 (0.0016),
  K on CPH −0.0017 (0.0023), CPHK20 −0.0001 (0.0031). Expected: no HARM for CP, CPH or the hook; R7b board-like.
- **Instrument pin:** 553's code as committed in 34ad97c9, plus the 555 fixes (commit 90ec5d3d, 2026-10-02): the sealed access
  arms only if the claim's sha256 is in the ledger's `## Reserve R7` section; malformed `--pairs` are refused before any
  phase; effects without a band are listed in plan.json before the truth opens and recorded by band553 with
  `band_total = None` (no longer dropped silently). The fix left every output of the synthetic end-to-end access unchanged
  (|d| ≤ 2.2e-16, cell-eval2's nmae jitter), and `test553b.py` passes. The blobs are written into the frozen claim and
  checked before the access.
- **Command (one access, ~25 min, < 4 GB, the truth opened only here):**
  `access553.py r7bd5000 --owner-go experiments/552-r7-reserve/CLAIM-R7.md --exclude-targets BST2 --arms $ARMS --seeds 0,1,2
  --pairs $PAIRS --half-arms $HALF_ARMS --half-pairs $HALF_PAIRS --ctrl-splits 5 --retro work/access553_rehearsal5
  --halves 0,1,2,3,4` (after its `--dry-run` with the same arms, pairs and halves).
- **Void** for R7b if an Oct 22 context matches Calu-3 before the access (the ledger's rule).

## 4. Order of operations

1. Fill §2's `R7A_F` line and §3's ARMS / PAIRS / HALF_ARMS / HALF_PAIRS and blobs from `docs/october-decisions.json` and the
   coordinator's message (quote the decided values and their timestamps). Freeze the file as `CLAIM-R7.md`; record its
   sha256 in the ledger's `## Reserve R7` section and commit, before any R7 step.
2. R7a: `r7a551.py plan`, then `r7a551.py record …`; append the "R7a armed" line it prints to the ledger's R7 section
   (`access` requires exactly one); then `r7a551.py access`.
3. R7b: blobs as pinned, the claim's sha in the R7 section (the access also checks it), `--dry-run`, then the access.
4. Append both outcomes to the ledger; report to the coordinator (numbers, verdicts, and the class of each change).

## 5. Count of tries (AGENTS §11.1c.4)

R7a: one access with F fixed before it. R7b: one access with its arms and pairs fixed before it. No R7 number has been
read by anyone before this claim is recorded.

## 6. The 10-07 runbook (window 06:00–12:00Z)

```bash
cd ~/Code/virtual-cell-challenge-2026
# 0. ~04:30Z: the coordinator's message gives F and the DE arm list. Fill, then freeze under the final name.
cp experiments/552-r7-reserve/CLAIM-R7-DRAFT.md experiments/552-r7-reserve/CLAIM-R7.md   # then edit the FILLs only
shasum -a 256 experiments/552-r7-reserve/CLAIM-R7.md      # -> a line in the ledger's '## Reserve R7' section
git add experiments/552-r7-reserve/CLAIM-R7.md docs/CONFIRMATION-LEDGER.md && git commit -m "R7 claim recorded" && git push
# 1. R7a: pin and arm (light); append the printed 'R7a armed' line to the ledger's R7 section and commit
uv run python experiments/551-r6-c10b/r7a551.py plan
uv run python experiments/551-r6-c10b/r7a551.py record --claim experiments/552-r7-reserve/CLAIM-R7.md --config <F> \
    --window-start 2026-10-07T06:00Z --window-end 2026-10-07T12:00Z --received <UTC> --message-file <message.txt>
# 2. R7b first at 06:00Z: blobs as pinned in section 3, the dry run, then the access
uv run python experiments/553-r7b-harness/access553.py r7bd5000 --dry-run --exclude-targets BST2 --arms $ARMS \
    --seeds 0,1,2 --pairs $PAIRS --half-arms $HALF_ARMS --half-pairs $HALF_PAIRS --halves 0,1,2,3,4
/usr/bin/time -l uv run python experiments/553-r7b-harness/access553.py r7bd5000 \
    --owner-go experiments/552-r7-reserve/CLAIM-R7.md --exclude-targets BST2 --arms $ARMS --seeds 0,1,2 \
    --pairs $PAIRS --half-arms $HALF_ARMS --half-pairs $HALF_PAIRS --ctrl-splits 5 --retro work/access553_rehearsal5 \
    --halves 0,1,2,3,4
# 3. R7a after R7b's process exits (~3-4 h): resumable; done steps are skipped on a re-run
nohup uv run python experiments/551-r6-c10b/r7a551.py access > experiments/551-r6-c10b/work/r7a/access.log 2>&1 &
# 4. Outcomes into the ledger's R7 section (both readings, verdicts, the board-likeness and HARM lines); report.
```

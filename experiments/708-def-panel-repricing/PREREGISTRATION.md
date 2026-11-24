# 708: post-release re-pricing of panel-dependent switches on D/E/F's REAL panel (preregistration)

Written 2026-10-15 ~21:45Z by the science session [d80a10], on the coordinator [2d6f46]'s YES (21:40Z). Written before
D/E/F exists and before any 705 / 707 board read.

## Why
Several October switches depend on the PANEL: 438 §4.4 showed `abcdef` ≈ `abc` within 0.00014, so the October-like draws
differ from A/B/C because their panels differ, not because of the day-of rules. 706 found the same for the strong-tier
direction swap. D/E/F's panel is drawn by Arc the same way A/B/C's was. After release, held-out PUBLIC lines (K562-GWPS,
Orion HCT116, Orion HEK293T as PLOLO folds) can price these switches on D/E/F's own 300 targets.
- Class: transfer, re-fitted on the day (AGENTS 11.1c.1). Nothing is fitted to D/E/F's truth (none exists), and there are
  no board reads in the dark.
- Only the LAST final entry counts, and the final runs until 11-05. So a move can ship as a later final entry with its
  own build and gates, after the automated R0 entry (10-22), which stays the decided recipe.

## Part A: the stack weights (H1 × CD4), `dayof438.py` verbatim
- Command (10-22, after R0 is submitted; ≤ 8.8 GB, ~30 min, plus ~15 min for the held-out fold operators):
  `uv run python experiments/438-october-stack-repricing/dayof438.py --bundle data/challenge-def --op <day's whitener>
  --cd4-root <day's cd4_compact> --tag def1022 --fold-ops held --held-ops <dir> --board-q 1.25`
  - The held-out fold operators are fitted on the day: `build_union_whitener.py --donors <the other two>
    --basis-prune essential`, 13 GB each, serial.
  - If that is impossible, `--fold-ops prod` (438 §6: composition steps within 0.00013).
- **The rule, verbatim from 438 §5 (recommend438.py):** the incumbent is H1 ×2, CD4 ×1. Candidates are H1 ×{1, 1.5, 2, 2.5,
  3} × CD4 ×{0.75 … 2}. A candidate moves only if it beats the incumbent on the LEVELLED k = 1 price by ≥ +0.0010 on the
  mean of the three cancer folds, on 3/3 folds, and still after leave-top-3-out. The most robust mover wins. q is the
  board's; q itself is not re-decided here.
- **Actionability:**
  - CD4 ×1.5 exists on the day path (`DEF_DAY_CD4_STACK=x1.5`).
  - An H1 dose other than ×2, or a CD4 dose other than ×1 / ×1.5, needs a dump-subshell pass-through from readiness
    with an inert known answer.
  - If that is not ready for the final build, ship the best SUPPORTED candidate that passes the rule, else the incumbent.

## Part B: the late switches (705's rider, 707's swap), on the same folds
- `dump708.py` = 522's `dump522.py` arithmetic on the D/E/F panel's folds: fold LODO-V1 key, percentile cap P 96.647,
  q 1.75 row sharpening, render, ballast KMAX 16; the 3 cancer folds; held-out operators as in Part A.
- `score708.py` prices, each against the incumbent percentile-cap key: (r) 522's θ 0.25 floor-preserving rider; (s)
  707's top-z-tertile β 0.25 swap; (rs) both. Price = 522's board units (0.0710 · r · ret + 0.1245 · Δpds_lit), k 1,
  ret 1.
- **Rules, by each switch's A/B/C board read:**
  - **705 CLEARED:** the rider ships as its prereg says. 708 can only VETO it: veto iff its D/E/F-panel price is ≤ −0.0005
    on the mean with ≥ 2/3 folds negative.
  - **707 TRANSFERS:** ships with dumphook707 and a rehearsal, unless 708 vetoes it (same veto rule).
  - **707 WEAK / NULL:** ships only if (s) passes 438's rule verbatim on the D/E/F panel (≥ +0.0010 mean, 3/3 folds,
    positive leave-top-3-out). This is 707 AMENDMENT 1.
  - **705 or 707 HARM:** off; 708 never reopens a HARM.
  - If both ship, (rs) must also pass the veto rule (an interaction check).

## Rehearsal before 10-21 (condition: or 708 is OFF)
- **Part A:** dayof438.py on the relabelled A/B/C bundle with 438's `abcdef` inputs: the A/B/C production operator
  (data/donors/whitener_union_pruned.npz), A/B/C's CD4 banks (data/external/cd4t/cloud-compact/cd4_compact), and the held
  fold operators of 062 (experiments/062-cd4-post-whiten/cloud-output/whiteners), with --board-q 1.25. Its table and
  recommendation must reproduce 438's stored `abcdef` outputs exactly: every comp438 row, and the RECOMMENDATION line
  under the final rule, obtained by `--skip-runs` on 438's stored abcdef outputs before the rehearsal runs.
- **Part B:** dump708.py on the same inputs must reproduce 438's stored `abcdef` "h2c1|V1_5_3+sh1.5" row on every fold
  (|d| < 1e-12, dump522's own known-answer check), and its per-donor terms must equal lodo_key (1e-12).
- **When:** step 7 after the science chain (sci1016.sh) ends, through gate602 (≤ 9 GB), or any earlier builder-free gap
  ≥ 45 min.

## Tries
- Part A: one frozen rule over 438's grid, which is the rule's own search; leave-top-3-out guards it.
- Part B: three frozen configurations (r, s, rs), no re-selection on D/E/F.
- Nothing here is tuned after D/E/F exists.

## What it does not do
- It never changes the R0 entry.
- It never reads D/E/F truth.
- It does not re-decide q, the key, the DE channel or any switch outside Parts A and B.

## AMENDMENT 1 (2026-10-15T21:39Z, before any measurement): Part A's winner is confirmed under the production key
- dayof438 prices the weights with the 410 key (K3) at scalar q. October ships the LODO percentile key (P 96.647) with
  tercilec sharpening. The K3 table cannot see a key × weight interaction.
- **Confirmation:** a Part A mover ships only if, in Part B's dumps (fold LODO-V1 key, percentile cap, q 1.75 render), the
  winner's composition (`FC.final("h<H>c<C>")`) priced against the incumbent `h2c1` has a mean ≥ 0 over the three cancer
  folds (522 units, k 1, ret 1). Otherwise the weights stay H1 ×2 / CD4 ×1.
- This can only block a move, never create one. The rehearsal reports it for abcdef's winner (H1 ×2.5, CD4 ×1.5) as
  information; on A/B/C the board already carries CD4 ×1.5 under the LODO key (ALLSW's dump is
  dump_h1x2_cd4x15_lodo05hi6).

## AMENDMENT 2 (2026-10-15T21:58Z, before any D/E/F data): the risk posture for the late switches (coordinator [2d6f46], 10-15 ~22:00Z)
Decision quoted from the coordinator's message: "Part A weights keep 438's rule verbatim (they will have board points).
707's swap and the rider θ get a RELAXED rule because their A/B/C reads predict ~0 by construction while October-like
draws say +0.001…+0.004 … Reason: we are behind on A/B/C and the payoff is convex, so positive-EV transfer levers with real
held-out support beat the lowest-regret default; negative-EV stays out."
- **Part A (weights):** unchanged: 438's rule verbatim, plus AMENDMENT 1's confirmation. 710 adds one condition: if 710's
  A/B/C board read of H1 ×2.5 is HARM (≤ −0.0003 overall), no 708 move may raise H1 above ×2.
- **Rider (r) and 707 swap (s), one rule each, replacing Part B's earlier routes:** adopt on D/E/F iff all four hold:
  - (a) the switch's A/B/C board read is not HARM (705: not HARM; 707: TRANSFERS, WEAK or NULL);
  - (b) the mean over the 3 D/E/F-panel folds is ≥ +0.0005 (522 units, k 1, ret 1);
  - (c) ≥ 2/3 folds are positive;
  - (d) the leave-top-3-out mean is positive.
- **Interaction:** if both are adopted, (rs) must pass (b)–(d) too; otherwise keep the one with the larger mean.
- **Supersedes:** Part B's veto-only route for 705 CLEARED and 707 TRANSFERS, and 707 AMENDMENT 1's 438-rule route for
  WEAK / NULL.
- **Note (science, for the record):** the rider's expected D/E/F-panel price is about +0.0005 (522: +0.0004 … +0.0009),
  so (b) sits at its centre and the rider may miss it by a hair. That follows from the decision as written.
- Implemented in decide708.py. The decision record links OCT22-CHECKLIST §C step 6a (245b8b1c), the gate for any
  replacement entry.

## AMENDMENT 3 (2026-10-15T22:02Z, before any D/E/F data): the rider's route (coordinator [2d6f46], 10-15 ~22:05Z)
The coordinator's clarification, quoted: "On (b): I meant the D/E/F-panel threshold for 707's swap. The RIDER's route is
different. It already has 9/9 October-like folds with all CI lower ends > 0, so if 705 reads CLEARED on A/B/C,
lodo_strength θ 0.25 is DECIDED ON before the 10-20 freeze and ships in R0. 708 then only VETOES it, if its D/E/F-panel
fold mean is ≤ −0.0005 (a HARM on the real panel); it does not need to re-earn +0.0005."
- **Rider:** 705 CLEARED → ON in R0. 708 vetoes it iff the D/E/F-panel mean of (r) is ≤ −0.0005. A veto ships as a
  replacement entry with the rider OFF, through OCT22-CHECKLIST 6a. 705 HARM → off, and 708 does not reopen it.
- **707 swap:** AMENDMENT 2's relaxed rule, unchanged: not HARM on A/B/C, plus D/E/F-panel mean ≥ +0.0005, ≥ 2/3 folds
  positive and leave-top-3-out > 0.
- **The swap is judged on top of what ships (science's consequence of the above).** When the rider is ON (and not vetoed),
  the swap's rule is applied to its MARGINAL price over the rider: arm (s|r) = (rs) priced against (r), with its own
  leave-top-3-out. When the rider is OFF, it is applied to (s) against the incumbent. This replaces AMENDMENT 2's
  "(rs) must pass too", and it never asks the rider to re-earn anything.
- Implemented in score708.py (arm `s|r`) and decide708.py.

## AMENDMENT 4 (2026-10-16T07:35Z, before any D/E/F data): Arc's basal floor in the held line, and board-sized columns (coordinator [2d6f46] approved, 10-16 ~07:33Z)
**Why.** Arc's panel rule (AGENTS 2.1) puts only targets with high basal expression in a screen. In a PLOLO fold the
held-out public line plays the context being predicted, so the Arc-like fold scores only targets expressed in THAT line.
`lib438.build_fold438` keeps every panel target the held bank perturbed. Measured on 438's fold caches (`floor708.py`,
c450c150; 801's A/B/C threshold 13.4122 CPM reproduced): targets below the floor in the held line show no measurable
knockdown there (0 % reach an own-gene log-effect < −0.5, against 24–76 % above it). On Arc's panels they are also less
predictable from the donors. They are 13–24 % of each Orion fold on A/B/C / real and about half of efp's; GWPS does not
measure 16–29 % of targets at all.
**What it changes on 10-22 (Part B and Part A's confirmation):**
- **Floor F**, from D/E/F's own controls by 801's rule: the D/E/F panel's 5th percentile of pseudobulk control CPM
  (1e6 × gene total / all counts per context file, averaged over D/E/F). No A/B/C number enters; on 180's relabelled A/B/C
  bundle the rule must give 13.4122 CPM (the rehearsal's known answer).
- **Floor-kept targets, per fold:** the held bank's control CPM on the official axis (1e6 · expm1(b_ctrl) / sum, as
  `build_fold438`'s gate) ≥ F. A target the held bank does not measure counts as below.
- **Every arm is priced twice** by `score708.py`, on all targets (the preregistered primary, unchanged) and on the
  floor-kept targets. That covers r, s, rs, s|r and A:<design>, each with its own leave-top-3-out inside its version.
- **Robustness rule:** a 708 ACTION needs its rule to pass in BOTH versions. The actions are the swap turning ON, the rider
  VETO, and Part A's AMENDMENT 1 confirmation (A:<design> mean ≥ 0). Otherwise there is no action: the swap stays OFF, the
  rider stays as R0 shipped it, and the weights stay h2c1. The rule can only block an action. On real / efp / efr it never
  binds:
  - all → floor-kept prices: swap 0.00478 → 0.00466, 0.00525 → 0.00609, 0.00805 → 0.00527;
  - rider: no veto in either version.
- If the floor step fails or its record is missing, `decide708.py` blocks every action and says why.
**Board-sized columns, for information only.** `score708.py` prints each fold's percentile cap `hi`, the tail statistic
(A/B/C 5.2–7.0, the lab's L panels 7.8–19.9, 199's real 9.3–16.0). It also prints each arm's price divided by the lab's
SAME-PANEL class ratio ρ (`experiments/800-checker/ka800.json`, `bulk-<class>|ABC`: key 1.46, shape 3.48, external 0.93).
The same-panel ratio applies because 708 prices on D/E/F's own panel, so D/E/F's tail is already in the number. BULK-W's
×4.75 compares L-panel prices with A/B/C board reads, so it mixes instrument bias with a panel effect. Until
`work/class708.json` assigns each arm a class (the coordinator asked the lab), all three alternatives are printed. The
preregistered thresholds do not change.
**Tries:** +1, the robustness variant. It is a fixed rule, not a search: one floor (801's rule, applied to D/E/F) and one
combination (AND).
- **AMENDMENT 4 note (07:36Z, information column only):** the ρ values quoted above are the lab's own units (ρ_lab). In 522
  units (score708's) the lab gives (10-16 ~07:35Z, grade804.json): r key/weight 1.25 (n = 1; 4.48 on its heavy-tail fresh
  panels, so ~4.5 if D/E/F's fold cap hi ≫ 7); s and s|r direction 0.68, provisional until 707's board read; A:<design>
  composition 0.68–1.25, no board point until 710. `work/class708.json` holds them, and score708 prints price / ρ_522.

## AMENDMENT 5 (2026-10-16T12:51Z, before any D/E/F data): the rehearsal's Part A known answer compares like-labelled stages (coordinator [2d6f46] OK, 10-16 ~13:00Z)
**What failed, and why it is not drift.** rehearse708 ran dayof438 on 180's relabelled A/B/C bundle at 12:24Z. It compared the
printout with 438's stored abcdef outputs and failed with 41 of 48 lines differing. The two runs take different branches of
the per-context claim rule, by design (ntclevel417's docstring; price438 and full438 line 7):
- A/B/C-labelled runs use 187's board-fitted claims (f = 1.00 / 1.00 / 0.9825; 187's identity ratios 1.0000 / 0.9927 /
  1.0000) with the fitted EV margin sd 0.002.
- Any D/E/F-labelled bundle, the relabelled rehearsal bundle included, uses f = 1 claims from its own NTC cells with sd 0.01.
  That is the 10-22 path.
- The inputs are identical: the three control matrices of 180's bundle and data/challenge match byte for byte (indices,
  data, indptr, sums), and so do the NTC ids. cp46 recomputed today on data/challenge equals the rehearsal's to 8 digits
  (B 14.80573414, C 16.43965014).
- Measured decomposition (`kaA708.py` on the graded output):
  - claims ref / got = **1.00000 / 0.99266 / 0.98249**, i.e. 187's fitted f;
  - **all unlevelled cells are identical on 17/17 design rows**;
  - every rho = ref rho × the claim ratio, to rounding (max deviation 0.0010 in the 3-decimal table,
    0.000117 in the header; the PP row reads 1.000 / 0.99287 / 0.98263);
  - the retention line is identical;
  - **the RECOMMENDATION is identical** (H1 ×2.5, CD4 ×1.5, q 1.25);
  - only levelled cells differ, by design.
- Reported, not compared: the q 1.5 RULE CHOICE (H1 ×2.5 vs ref H1 ×3, both at CD4 ×1.5). The rule reads levelled prices,
  and q 1.5 is not the board's q.

**The corrected Part A known answer (`experiments/708-def-panel-repricing/kaA708.py`):**
- (i) every claim-independent cell is identical: the design table's unlevelled columns on every row, and the retention line;
- (ii) every rho equals the reference rho × claim_ref / claim_got per context, within 0.002 (table) / 0.0004 (header)
  relative;
- (iii) the RECOMMENDATION is identical.
- A negative control fails: one perturbed unlevelled cell gives FAIL.
- **Graded now: PART A KA PASS** on the 12:24Z output `work/dayof438-r708abcdef.txt` (sha256 7f054544…, git blob 39aca2a1;
  commit 20cb2bc8), result in `work/kaA708-r708abcdef.json`.

**Procedure from here:**
- rehearse708 grades Part A with kaA708 (reusing that graded output), then runs Part B (dump708's 438 known answer, |d| < 1e-12)
  and Part C (score708), unchanged.
- rehearse708b then runs day708 end to end, AMENDMENT 4's floor step included. It starts only while no build or R-NEW chain
  is alive and before 01:30Z, or else after R-ABC2's chain (rabc2_695) has exited, so it cannot run into R-ABC2's 10-17 build.
- **decide708 now requires the rehearsal on file:** kaA708's json PASS, rehearse708's REHEARSAL PASS and rehearse708b's
  completed day path. Otherwise every 708 action is blocked: 708 is OFF and R0 ships as decided.
- Tries: unchanged. No choice was re-made; the comparison was re-specified to like-labelled stages.

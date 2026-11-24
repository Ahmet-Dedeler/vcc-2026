# 708 — D/E/F-panel re-pricing: running record

The preregistration and its amendments are in `PREREGISTRATION.md`. The day path is `day708.sh`; the rehearsals are
`rehearse708.sh` (Part A known answer + Part B dumps) and `rehearse708b.sh` (day708 end to end on 180's relabelled A/B/C).
Neither rehearsal has run yet: both are queued behind the science chain.

## 2026-10-16 — the basal floor in the held line (AMENDMENT 4) and the key tail

**Question (the coordinator, relaying the lab's 801 diagnostics):** efp / efr carry 156 / 68 targets below A/B/C's basal
floor, which Arc's selection excludes. Does that explain why October-like panels price the strong-tier swap (706/707) high
while A/B/C's panel does not? And should 708's fold truth be restricted on 10-22?

**Which floor.** In a PLOLO fold the held-out PUBLIC line (K562-GWPS, Orion HCT116, Orion HEK293T) is the context being
predicted. So the Arc-like condition is that the target is expressed in that line, not in A/B/C. The floor follows 801's
rule: the panel's 5th percentile of pseudobulk control CPM, which `floor708.py` reproduces on A/B/C as 13.4122 CPM to
|d| 0. A fold target is below it when its held-line control CPM, on the official axis (build_fold438's own b_ctrl
arithmetic), is lower, or when the held bank does not measure the gene at all.

**What below-floor targets are** (`floor708.py`, 438's fold caches + the banks' delta; `work/floor708.json`):

| panel / fold | n | unmeasured | below 13.41 | own-gene kd median (share < −0.5): below / above | cos(pre, truth) median: below / above |
|---|---|---|---|---|---|
| abc / K562 | 269 | 43 | 47 | −0.30 (0.00) / −0.76 (0.76) | 0.015 / 0.033 |
| abc / HCT116 | 290 | 0 | 56 | −0.14 (0.00) / −0.56 (0.59) | 0.018 / 0.033 |
| abc / HEK293T | 294 | 0 | 38 | −0.10 (0.00) / −0.39 (0.28) | 0.019 / 0.033 |
| real / HCT116 | 300 | 0 | 71 | −0.03 (0.00) / −0.51 (0.52) | 0.003 / 0.026 |
| efp / HCT116 | 293 | 0 | 156 | −0.00 (0.00) / −0.51 (0.50) | 0.018 / 0.023 |
| efr / HEK293T | 300 | 0 | 69 | −0.03 (0.00) / −0.36 (0.26) | 0.022 / 0.017 |

Below-floor targets show no measurable knockdown in the held line. None reaches an own-gene log-effect of −0.5; above the
floor 26–76 % do. Their truth norms are as large as the rest, so their fold truth is mostly noise or a line-specific
response. On Arc-like panels they are also less predictable from the donors. GWPS's expression axis does not contain
16–29 % of panel targets. Full table: `work/floor708.json`.

**Does it move 708's prices?** `floorprice708.py` runs score708's exact arithmetic on 522's fold dumps, once on all
targets and once on floor-kept targets. The all-target column equals score708's self-test v2 to every digit. Results
(`work/floorprice708.txt`):

| panel | arm | all targets | floor-kept | verdict change |
|---|---|---|---|---|
| real | r | +0.00051 (3/3) | +0.00042 (3/3) | none |
| real | s | +0.00478 (3/3) | +0.00466 (3/3) | none |
| real | s\|r | +0.00476 | +0.00465 | none |
| efp | r | +0.00099 (3/3) | +0.00055 (2/3) | none (no veto either way) |
| efp | s | +0.00525 | +0.00609 | none |
| efr | r | +0.00129 (3/3) | +0.00104 (2/3) | 438's move rule only (not a 708 rule) |
| efr | s | +0.00805 | +0.00527 | none |

**Reading.** The floor changes the swap's price by −35 % to +16 % and the rider's by at most −0.0005, and it never flips a
708 decision. So the floor does not explain the October-like versus A/B/C gap.

**What differs is the key tail.** The percentile cap `hi` (P 96.647 of key / median over carried targets) is:
- 199's `real` draw: 9.3 / 16.0 / 10.3. This is our own October-like draw, not an Arc panel.
- efp: 8.8 / 17.2 / 9.6.
- efr: 8.8 / 16.8 / 10.3.
- The lab's L panels, matched to A/B/C on coverage, strength and A/B/C basal expression: 7.8–19.9.
- A/B/C alone: 5.2–7.0.

Restricting to floor-kept targets leaves the caps about where they are (real 8.9 / 17.6 / 8.6). Below-floor targets hold
only 6–26 % of the key mass and are not concentrated in the top decile. A/B/C's light tail is the odd one out. Whether
D/E/F shares it is measured on the day: score708 prints each fold's `hi`.

**Decision (AMENDMENT 4, coordinator-approved 10-16 ~07:33Z):**
- On 10-22 the floor F comes from D/E/F's own controls by 801's rule.
- Every arm is priced on all targets and on floor-kept targets.
- A 708 action (swap ON, rider veto, Part A confirmation) must pass in both versions.
- Board-sized columns (price / same-panel ρ in 522 units, the lab's `class708.json`) are printed for information only.
- Tries +1 (the robustness variant).

**Implementation and tests** (58b52f85; `work/test_amend4_708.txt`):
- floor708 day mode; score708 `--floor-json`; decide708's both-version rule (a missing floor blocks every action); day708's
  floor step.
- KA1: score708's new primary json equals the HEAD code's, exactly.
- KA2: the floor column equals floorprice708, |d| 0.
- KA3: equals self-test v2.
- On 180's relabelled A/B/C bundle F = 13.412203349132618, |d| 0; day708 asserts it in rehearse708b.
- The day708 block runs under bash 3.2, rc 0.
- test_decide708: 6 cases PASS.

**On the over-read (the coordinator's item 2).** BULK-W's ×4.75 is the checker's price on fresh L panels divided by the
board's read on A/B/C. It mixes instrument bias with the heavy-tail panel effect. 708 prices on D/E/F's own panel, so its
tail is already in the number, and the right correction is the same-panel ratio. In 522 units the lab gives:
- r (key / weight): 1.25; 4.48 if D/E/F's tail is heavy.
- s (direction): 0.68, provisional until 707's board read.
- A:<design> (composition): 0.68–1.25.

Information only; the preregistered thresholds stand.

## 2026-10-16/17 — the rehearsals: Part A KA PASS (AMENDMENT 5), Part B / C, and the day path end to end (recorded 2026-10-16T21:33Z)
- **rehearse708:** REHEARSAL PASS. Part A was graded by kaA708, like-labelled (AMENDMENT 5; the graded output is commit 20cb2bc8).
  Part B reproduced dump708's 438 known answer on 3/3 folds, and Part C (score708) ran.
- **rehearse708b:** day708 ran end to end on 180's relabelled A/B/C bundle. It used the 10-22 procedure: 453's day operator
  66610be7 plus three held-out fold operators fitted fresh with the day recipe (`--basis-prune essential`, 12k cap).
  **day708 rc 0, DECISIONS WRITTEN** (21:06:18Z). AMENDMENT 4's floor step reproduced 801's 13.412203349132618 CPM (|d| 0).

### (a) Which step changed the choice: Part A's own rule, on the day path's operators. Not AMENDMENT 1, not the floor, not q.
- rehearse708b compares the two runs' Part A RECOMMENDATION lines. Those come from dayof438's 438 rule, before Part B, so
  AMENDMENT 1's production-key confirmation and the floor step had no mover to judge.
- The day path's Part A printed "→ design incumbent (h2c1)", and decide708 wrote "WEIGHTS stay H1 x2 / CD4 x1 (no mover passed
  438's rule)". q is 1.25 in both runs.
- The difference is the operators. The rehearsal's Part A used the A/B/C production operator e2f65ecf and 062's unpruned held
  fold operators, so that it could reproduce 438's stored computation (its known answer). The day path uses the day's operator
  and fresh day-recipe fold operators, which is what runs on 10-22.
- With those operators the levelled K3 surface at q 1.25 is lower, and its high-dose corner loses the HEK293T fold:

| levelled k=1, q 1.25 | rehearsal (062 ops) | day path (day-recipe ops) |
|---|---|---|
| incumbent H1 ×2 / CD4 ×1 | +0.00239 | +0.00215 |
| A/B/C weights H1 ×2 / CD4 ×1.5 | +0.00408 | +0.00352 |
| H1 ×2.5 / CD4 ×1.5, per fold (K562 / HCT116 / HEK293T) | +0.0071 / +0.0054 / +0.0010 | +0.0069 / +0.0046 / +0.0000 |
| H1 ×3 / CD4 ×1.5, per fold | +0.0077 / +0.0058 / +0.0004 | +0.0075 / +0.0050 / −0.0006 (2/3) |
| argmax | H1 ×3 / CD4 ×1.75 +0.00476 | H1 ×3 / CD4 ×1.5 +0.00393 |
| movers (paired ≥ +0.0010, 3/3, leave-top-3-out > 0) | several → H1 ×2.5 / CD4 ×1.5 | **none** → incumbent |
| fold retention ret_f | 1.001 / 0.999 / 1.000 | 0.986 / 0.973 / 0.989 |

- This matches 713's production-key read on LABC, the same day recipe built independently by the lab. There, H1 above ×2 does
  not pay at CD4 ×1: ×2.5 −0.00012, ×3 −0.00064.

### (b) "DAYPATH different choice" is an expected outcome of the day path, not a defect
- rehearse708b's own header made the comparison informational: "a different choice is reported, not failed". Its pass
  condition is day708 rc 0 with decisions written. That is what decide708's rehearsal-on-file check requires (AMENDMENT 5), and
  that is the right condition.
- The rehearsal's Part A and the day path use different operators by design: one tests the code path against 438's stored
  outputs, the other is the 10-22 procedure.
- The day path is also validated end to end. Its Part B prices on A/B/C's panel equal the lab's independently built LABC
  (EC2, ABCp's fits) to 5 digits: the rider (θ 0.25) +0.00005 (2/3) against LABC's +0.00005, and the swap (β 0.25) +0.00159
  (3/3) against 713's +0.00159. The operator fits are not bit-reproducible across machines, yet the prices agree.

### (c) What it predicts for 10-22, if D/E/F prices like A/B/C (an Arc-drawn panel)
- **Weights:** Part A makes no move, so R0's decided weights (H1 ×2, CD4 ×1) are the final weights.
- **Rider:** r ≈ +0.00005 (floor-kept +0.00007), far from the −0.0005 veto. The rider stays as R0 shipped it (ON if 705 reads
  CLEARED).
- **Swap:** s|r +0.00157 (3/3; leave-top-3-out +0.00097), floor-kept +0.00165 (3/3; +0.00113). Both clear AMENDMENT 2's relaxed
  rule (mean ≥ +0.0005, ≥ 2/3 folds, leave-top-3-out > 0) in both versions.
  - **So if 707's A/B/C board read is not HARM, an A/B/C-like D/E/F panel turns the swap ON:** a replacement entry through
    OCT22-CHECKLIST §C 6a. It is the one 708 action this rehearsal predicts. The rehearsal's own FINAL line shows swap OFF only
    because it ran with 707=NONE.
  - Board-sized (same-panel ρ 0.68, provisional): ≈ +0.0023.
- **Tail:** the day path's fold caps on A/B/C's panel are 6.11 / 7.82 / 8.44 (062 operators: 5.2–7.0; the lab's L panels:
  7.8–19.9). D/E/F's own caps will print on the day.

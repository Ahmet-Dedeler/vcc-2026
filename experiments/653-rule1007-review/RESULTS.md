# 653 — adversarial read of the 10-07 / 10-08 rule (537's `rule1007_draft537.py` v3, sha 8367c89d) and its §5.1 tree

> **FROZEN (coordinator, c5cdd3bc, 10-02 17:14Z):** `rule1007_v4_653.py` IS the 10-07 / 10-08 slot rule, sha256
> `e00af6c7d126d4768dfc9f8c5d1541f44d80d98bba40f10fefe0383fd18128ce`. chain1007 asserts that sha. Never edit the file.
> It runs without `--arc-hook`, and with `--scored=` on 10-08. R4 is adopted as a CPH prereg amendment.

2026-10-02 17:05–17:13Z (commit time), council session [7772dc], at the coordinator's request, before the rule freezes (due before 10-06
00:01Z). Read-only toward 537: its file is unchanged (sha checked before and after).

**Files (all here):**
- `repro653.py`: every finding below, reproduced against v3 itself;
- `rule1007_v4_653.py`: a suggested v4 with selftests;
- `rule1007_v3_to_v4.diff`;
- `equiv653.py`: v3 vs v4 over every bin state.

## 0. Verdict
**The 10-07 lists are correct for every reachable input.** Today's inputs give:
- 10-07 s1: **CPH**;
- 10-07 s2: **K20** if k16 reads ≥ MID, else **459dr**.

The reachable tree is small. CP is already read (NULL + rescue), so every CP-POS/PARTIAL branch is dead, and MZ can never
enter a 10-07 list (`cp_on` is True through the rescue flag).

**Two findings matter before the freeze:**
- **R1:** the 10-08 list can silently drop 459dr.
- **R6:** the rule ignores 524dr, whose LAG bin says "before any other October work".

**Three are cheap hardening:**
- R2: a missing capture crashes the rule.
- R3: it has no bulk-leak falsifier.
- R4: CPH "SHIP" ignores reach and nmae.

**One is a test gap:** R5, the selftest that would catch R1 accepts either answer.

R7–R11 are scheduling and record items.

## 1. Findings
| # | finding | evidence (repro653.py unless noted) | fix (v4 implements R1–R4, R8, R11) |
|---|---|---|---|
| **R1** | **The 10-08 list excludes the heads of the RECOMPUTED 10-07 lists, not the arms actually scored.** On 10-08 an operator who passes K20's id (it WAS scored, on 10-07 s2) makes the recomputed s2 head 459dr. d8 then drops 459dr although it was never read. | Expected branch (CP NULL + rescue, 496a FAILS, k16 POS), CPH read NO: v3 gives `['ALLSW']` with K20's id passed and `['459dr', 'ALLSW']` with `none`. 128 sweep states mis-list (equiv653). The docstring says the 4th argument means "scored on 10-06", but nothing stops the natural misuse on 10-08. | v4 takes `--scored=<tags scored on 10-07>`. Without it, it excludes only CPH / K20 when their ids are given and never a recomputed head; the chain's novelty gate refuses any package that was in fact scored. With `--scored` = the actual picks, **v4 equals v3 in all 384 states** (equiv653) |
| **R6** | **524dr (10-06 s1) is listed as an input (DR) in §5.1 but the rule never reads it.** 524's prereg: on LAG, "before any other October work, find the lagging member". A switch that did not reach the emitter on the October path comes out of the recipe. | rule1007_draft537.py has no DR argument; 524's PREREGISTRATION bins | Read 524dr against 453dr and Σ (524's reader) BEFORE the 10-07 waiter fires. On LAG or a falsifier, 10-07 s1 becomes the diagnosis (or the slot is held), because every further production-path read is worth less than an October path that drops a switch. Not coded in v4 (Σ needs 524's reader) |
| R2 | A missing capture (`data/board-entries/<id>.json`) raises FileNotFoundError, so a waiter that calls the rule dies with no pick. | R2 in repro653 | v4 reads it as UNKNOWN, prints a WARNING and reports `MISSING=` in the summary |
| R3 | No falsifier. 537 §5.2 says "pds and mse must equal 446bj2z's to ≤ 1e-6 raw in every DE arm (falsifier: a bulk leak — stop)". A k16 that moved pds still bins POS. | R3 in repro653 | v4 bins it LEAK (treated as off) for CP / H / K / CPH. 496a's base is 446bj2 (its package name has no zs), and the pin is checked against that base |
| R4 | CPH "SHIP" = Δ(fid + jac) ≥ +0.0012 only. With the hook on, nmae and reach move too, so an arm whose reach or nmae loss makes the overall delta ≤ 0 still ships the block into October. | R4: fid + jac +0.0013 with overall −0.0005 → SHIP | v4 adds an overall floor (Δ score_avg vs 446bj2z > 0). Adopt it as a CPH prereg amendment before 10-07's read |
| R5 | The selftest's last 10-08 assertion accepts `['ALLSW']` OR `['459dr']`, so it passes with R1's bug. | v3 line 162 | v4 asserts `['459dr', 'ALLSW']` under both argument conventions, and adds R1–R4 / R8 regressions |
| R7 | **ALLSW and S1 are not built.** Every other named package exists: CPH `.vcc` 10-01 01:02Z; K16 / K20 / MZ 09-29; 459dr 09-26; 496a 09-28. In the expected branch the 10-08 list is `[459dr, ALLSW]` (or `[ALLSW]` under R1), so 10-08 s2 has nothing READY unless ALLSW is built by ~10-07 18:00Z. | `ls submissions/` | Append the coordinator's post-tree queue (13:08Z, re-ordered 14:53Z: 630's dose pair, 524dr2, factorial / all-switches, the hook-form arm, 625) as d8's tail, each gated on READY, so a slot is not left empty. 630's pair needs the 10-06 apply + KA first |
| R8 | MZ / CPMZ (the Arc hook) appear when 496a ≥ PARTIAL. 652 F1 recommends dropping the Arc hook from the October path: circular held-out read, an attenuation fitted to 458d's A/B/C read, and 629. Also, 496's prereg says PARTIAL is "not shipped without a second read", yet here PARTIAL triggers a different arm. | 652 F1 | v4 lists MZ / CPMZ only with `--arc-hook`. If kept, preregister that no A/B/C read switches October to the Arc hook |
| R9 | k16 (10-06 s2) must be captured before the 10-07 waiter evaluates the rule. A slow 10-06 upload (09-29 crawled for hours) leaves K UNKNOWN. | — | v3 / v4 both fall back safely (no K20 on 10-07; K20 re-enters 10-08 through `k_off`). Capture k16 immediately after it publishes (`capture_board_entry.py`) and check the file before the waiter fires |
| R10 | 459dr is a 09-26 package. It reads H1 ×3 vs ×2 on 453's October path, before LODO, zs and tercile, so it calibrates 438's H1 axis on an older path. | 459 PREREGISTRATION | Still a valid pair vs 453dr (its question is the instrument's H1 axis). If it reads ≥ +0.0015 (wire dayof438's H1 choice), confirm on the current path inside 524dr2 before October uses it |
| R11 | Scaled deltas are anchor-dependent. All 38 captured entries are on anchors `valA/B/C-r4`, but a mid-stream re-anchor would corrupt the scaled bins (CP / K / CPH) silently. | equiv check of `anchor_version` | v4 refuses on an anchor mismatch |

## 2. Reachability (what the tree can still do)
- **10-07.**
  - CP is fixed: NULL + rescue, so `cp_on` is True and `k_on` = (k16 ≥ MID and K20 unscored).
  - Reachable outcomes: s1 = CPH (if READY), s2 = K20 or 459dr.
  - MZ is unreachable on 10-07. So are the selftest's cpPOS / cpNUL branches (correct, but dead).
- **10-08.**
  - `block_on` = CPH SHIP. Its re-priced chance is low: −0.0001 … +0.0009 against the +0.0012 bar.
  - Expected: CPH NO → `[K20 if k16 ≥ MID and still unread] + [MZ if H on] + [459dr, ALLSW]`.
  - CPHK20 and CPMZ are reachable only on SHIP, and need builds on 10-07 (§5.3: not before CPH reads SHIP).
- **C (495, 10-03)** only adds S1 to the tail, and S1 is unbuilt.

## 3. What v4 is and is not
- **Each change is one of R1–R4, R8 or R11, each tested.** Selftest: v3's 7 branches + 10-08 checks + regressions,
  PASS. On today's inputs v4 prints exactly v3's lists.
- **equiv653: with `--arc-hook` and `--scored` = the actual 10-07 picks, v4 == v3 in all 384 states.** v4 changes nothing
  when v3 is used as documented. It only removes R1's foot-gun and adds the falsifiers.
- **What it does not do:** 524dr's bin (R6) and the post-tree tail (R7). Both are the coordinator's chain logic.
- **It is a suggestion:** 537's file is untouched, and the coordinator decides what to freeze.

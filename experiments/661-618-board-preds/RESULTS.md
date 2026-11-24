# 661 — 618-registered board predictions for the queued DE arms (coordinator's 661 items)

[21dd4c], 2026-10-03. Predictions committed before each board read:
- PRED-661-625.md (e8edc6cb);
- PRED-661-630.md, PRED-661-dnmag.md and PRED-661-idall3.md (00ad7165).

## (4) 625's 8 single-octile prior-only flips on 618 (`flip618.py`, `score661oct.py`, `work/score661oct.json`; 21:36Z)
- **KA661 PASS:** `flip618 h1 none` == 612's c4kg s0 stats, bit for bit (pre-zs and zs). The copy reproduces flip625's
  base.
- **Runs:** 618 base (c4kg on the trainval surface, zero repair −5) plus octiles 0–7, 3 seeds each.
- **Bulk:** pinned. Sums differ from the base in at most 4 cells of 1 target, by 1 count (relative 1e-9). That is zero
  repair's floor, as in 618's hook runs.

| octile | Δnmae f 0.35 (seed sd) | Δfid f 0.35 | Δcalls / target | Δnmae f 1 | Δfid f 1 | 625 rule f 0.35 | f 1 |
|---|---|---|---|---|---|---|---|
| 0 | +0.00436 (1e-5) | −0.00336 | +4.7 | +0.00643 | −0.00328 | fail (degrades) | fail |
| 1 | +0.00298 (4e-5) | −0.00217 | +2.3 | +0.00314 | −0.00223 | fail | fail |
| 2 | +0.00246 (2e-5) | −0.00135 | +1.2 | +0.00192 | −0.00150 | fail | fail |
| 3 | +0.00060 (4e-5) | −0.00016 | −0.1 | +0.00150 | −0.00023 | fail | fail |
| 4 | **−0.00332** (5e-5) | **−0.00278** | **+156** | −0.00243 | +0.00279 | fail (fid < −0.002) | **pass** |
| 5 | **−0.00236** (2e-5) | **−0.00643** | **+353** | −0.00065 | +0.00190 | fail (fid) | **pass** |
| 6 | +0.00089 (6e-5) | −0.01225 | +1186 | +0.00168 | +0.00185 | fail | fail |
| 7 | +0.00394 (7e-5) | −0.01787 | +1238 | +0.00512 | −0.00636 | fail | fail |

The bar at f 0.35 is max(0.25 × median |x| = 0.00068, 2σ). The base has 2,569 calls per target.

- **Measurability (coordinator's test): every octile moves nmae by ≫ 2 per-seed sd** (effects 6e-4 … 6e-3 vs sd
  1e-5 … 7e-5).
- **618's prediction for 625's rule:**
  - at f 0.35, the board-like regime, **no octile passes**: 0–3 and 6–7 degrade, and 4–5 lower nmae but fail the fid
    no-harm check;
  - at f 1, octiles **{4, 5}** pass;
  - both 4 and 5 sit on large call-set expansions (+156 / +353 calls per target): **mechanical, not sign quality**.
- **Against the registration (PRED-661-625):**
  - measurable ✓;
  - {0, 1, 2, 6, 7} degrade ✓;
  - octile 3 near the bar ✓ (it degrades, +0.0006 < bar);
  - octile 5 lowers nmae with dcalls ≥ +100 ✓, and octile 4 behaves the same way (registered as "possibly 4");
  - the registered passing set {5}(+3/4) matches at f 1 ({4, 5}); at f 0.35 the fid check empties it.

## What it means for 625's 9 board slots
- 618 is an H1-line surface whose programme agrees with CLEAN4, the same as h1-test. On such a line 625's probes read NULL
  at board-like regimes, or "pass" octiles 4 / 5 for a mechanical reason (more calls) at f 1. Neither is a CLEAN4
  sign error.
- **Clean octiles vs confounded octiles:**
  - Octiles 0–3 flip with no call-set change (|Δcalls| ≤ 5 per target). Their nmae / fid read sign quality.
  - Octiles 4–7 flip DOWN priors to UP and expand the realised call set by 156–1,238 calls per target. Their reads are
    dominated by call-set mechanics.
- **So a 625 read can only find a real A/B/C disagreement with CLEAN4 in octiles 0–3.** A pass in 4–7 on the board would
  most likely repeat this mechanism, and a D/E/F mask built from it would carry call expansion, not better signs.
- **Suggestion (the coordinator's call):**
  - if 625 is spent, restrict it to octiles 0–3 plus one combination entry: 5 slots instead of 9, freeing 4;
  - add the rule "a passing octile also needs |Δcalls| ≤ 10 per target in its board rc diagnostics".
  - If A/B/C lines are as CLEAN4-consistent as H1, 618 says even those 5 read NULL. The expected value is then ~0
    against the +0.002–0.003 if a real disagreement exists.

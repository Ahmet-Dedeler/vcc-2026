# 711: the call-set size curve's next point, K 4,200 (CPM 3,800) on 446bj2z, and the October K rule (preregistration)

Written 2026-10-16 ~05:10Z by the science session [d80a10], on the coordinator [2d6f46]'s decision (10-16 ~05:05Z). At this
time neither 711's build nor the POST HOC H1 size curve (701 run701b / reach701) exists.

## Why
K-UP (701, ku30 = 446bj2z + CPM 3,000, K 3,400; entry HAbelEdgYYziqwfjTgdI, 0.2740045) read NULL on its preregistered
primary Δ(fid + jac) (+0.000123 overall), with full delivery (3,339–3,373 of 3,400 calls per target). But reach moved
+0.005398 raw (+0.00108 overall) against flat reach at K 1,600 (+0.00035) and K 2,000 (+0.00060).
- Call-set size is a DE knob that leaves pass 1 identical (537 §2.5). Its paired noise is 0.3–0.6× an independent re-roll
  (reach σ ≈ 0.0004–0.0008 raw), so the jump is ~7–13σ, probably a real design effect.
- October needs the size curve, not a re-draw: if reach keeps paying above 3,400, the October K should move. This read is
  the next dose point, on the same base.

## The arm (one variable vs 446bj2z, entry `be90smZUDcIV7j5MgQ03`, 0.2728260)
- **446bj2zku38** = 446bj2z's exact env + `VCC_RC_CALLSET_CPM=3800`, so every target's call set is CallSet(400, 3800):
  K 4,200. The knob's validated range is [100, 5000].
- Built by run_arm_generic.sh against 446bj2z's meta, exactly as 701 / 523's K16 / K20.
- Expected config diff: `rc_callset_cpm` (+ `de_sign_emitter_sha256`, the emitter pin, which has no output effect).
- Pass 1 is untouched, so pds and mse must equal 446bj2z's to ≤ 1e-6 raw.

## Class and noise
**Scorer mechanics.** No biology changes. It carries to D/E/F because the final panel is drawn by the same design and
read by the same scorer. DE-knob paired noise as above.

## Readings (all vs 446bj2z, and vs ku30 beside them)
- **Primary: Δreach raw.** Also reported: Δ(fid + jac) in overall units (701's statistic), Δnmae raw and Δoverall.
- **Prediction (mechanism: the added ~800 block genes become called, so more of the reference's DE genes enter reach's walk
  early):** if the mechanism is monotone, Δreach ≥ ku30's +0.0054. fid / jac stay near flat (size above n_real neither adds
  coverage nor dilutes much), and jac may slip further (ku30: −0.00088 raw).
- **Falsifiers (STOP and diagnose before reading):**
  - pds or mse moves by more than 1e-6 raw (a bulk leak);
  - a config diff beyond the expected keys;
  - nmae moves by more than ±0.002 raw;
  - realised calls (the build's rc audit, via 701's deliver701.py) grow by less than half the prescribed +800 over
    ku30. That reads as "capacity-bound", not as "size does not pay".

## The October K rule (the coordinator's, fixed now, before 711's read and before the H1 curve)
October K = the size with the largest board Δoverall among {2,600 (446bj2z, Δ 0), 3,400 (ku30, +0.0011785), 4,200 (711)},
**only if** both hold:
- **(i)** the POST HOC H1 runs show reach rising with K at paired-seed significance. `reach701.py`'s rule: at f 0.87, all
  3 per-seed slopes of reach on CPM (sizes 2200 / 2800 / 3000 / 3400 / 3800) are > 0, the mean slope is > 2 × its seed
  standard error, and the mean paired Δreach at CPM 3800 is ≥ the one at CPM 3000;
- **(ii)** that size's Δoverall exceeds the next-smaller size's by ≥ +0.0005.

Otherwise K stays 2,600 (`rc_callset_cpm` off). If adopted, it is a decision flip of `rc_callset_cpm` through the
existing DEF_DAY_RC_CALLSET_CPM pass-through, plus one rehearsal.

## Tries
- Board, on the call-set size: K 1,600, 2,000, 2,600 and 3,400 so far (4 sizes); 5 after this read.
- H1: 537's earlier 2 sizes × 3 seeds, plus run701b's 4 sizes × 3 seeds (POST HOC: labelled as such; the rule above was
  fixed before they ran).
- The rule picks among 3 board sizes with a +0.0005 step margin and an H1 mechanism gate, so best-of-3 luck alone cannot
  pass it.

## Amendment 1 to the October K rule (the coordinator [2d6f46], 2026-10-16T11:34Z, BEFORE 711's entry exists)
Science [d80a10] found an edge the text above handles badly. K 4,200 can lead but beat K 3,400 by less than +0.0005 while K 3,400 has
already beaten K 2,600 by +0.0011785. As written, the rule then keeps K at 2,600 and discards a size that cleared its own step. That was
not the intent: (ii) exists so that a size counts only when it is distinguishable from the next-smaller one, not to punish a curve that
flattens.
**Amended (ii), the step-down:** start at the size with the largest board Δoverall. If its Δoverall exceeds the next-smaller size's by
≥ +0.0005, adopt it. Otherwise step down to the next-smaller size and apply the same test against ITS next-smaller size. Repeat until a
size passes, or until K 2,600 (which needs no test). Condition (i) is unchanged and applies to any adopted size above 2,600.
- Under the amendment, the only case that changes is "4,200 leads by < +0.0005 over 3,400", which now adopts K 3,400 (if (i) holds).
- Try count unchanged (3 board sizes). The margin and the H1 gate still stop best-of-3 luck. It costs no new board read.
- read711.py implements the amended rule; its dry-runs are in experiments/670-coord-1015/READ-1017.md.

## Amendment 2: the H1 dial test (the coordinator [2d6f46], 2026-10-16T12:37Z, BEFORE ku38's package or read exists)
**Facts at this time.** Condition (i) at the preregistered dial f 0.87 FAILS (reach701.json, 62553c07): H1 reach falls with K (slope
−0.0104 per 1,000 CPM, 3/3 seeds), with Δreach vs K 2,600 of −0.016…−0.022 at every larger size. That dial already fails to retrodict the
one board point that exists: at K 3,400 (CPM 3000) it says −0.0177 [corrected, see below], and the board's ku30 read +0.0054 raw. The dial came from FINDINGS 0b ("the board
sits near H1's f 1 regime for SIZE changes"). 662's calibration puts the board's call-set level at f ≈ 0.61. Science's exploratory reads
at f 0.61 (seen after the result, so they decide nothing alone): K 3,400 +0.0046 (ku30's board +0.0054), K 4,200 −0.0016.
**The test.** ku38's board read is unseen and discriminates the two dials:
- f 0.61 predicts board Δreach (ku38 vs 446bj2z) ≈ −0.0016: K 4,200 falls BELOW ku30.
- f 0.87's shape and the monotone mechanism (this prereg's own prediction) predict ≥ +0.0054.
- Threshold, fixed now: **+0.0019 raw** (about the midpoint; ≥ 4 paired-noise σ from each prediction, σ ≈ 0.0004–0.0008).
**Consequence for October K, replacing condition (i) ONLY if the board confirms f 0.61:**
- ku38 Δreach ≤ +0.0019 → **f 0.61 is validated out of sample** (it retrodicts ku30 and predicts ku38). Condition (i′) replaces (i):
  at the candidate from amended (ii), H1's paired Δreach vs K 2,600 at f 0.61 is > 0 on 3/3 seeds and its mean > 2 × se. Computed
  with reach701's exact arithmetic, regime 0.61 instead of 0.87, into experiments/701-kup-callset/work/reach701_f061.json, which
  science writes BEFORE ku38's read. If (i′) fails at the candidate, step down as in Amendment 1.
- ku38 Δreach > +0.0019 → f 0.61 is refuted, condition (i) at f 0.87 stands, and K stays 2,600.
- Either way, an adopted K enters October only after R-ABC3 rehearses it on def_day (gate695 allows K-only) before the 10-20 23:59Z freeze.
**Tries.** Board sizes 5 after this read. Dials considered 2 (0.87 preregistered, 0.61 from 662, which existed before this prereg).
The dial is chosen by an unseen board read, not by fit, so a best-of-2 dial choice is not rewarded unless the out-of-sample prediction holds.
- **Correction (2026-10-16T12:43Z, science [d80a10] caught it; factual only, no rule or threshold changes):** Amendment 2 first quoted f 0.87's
  K 3,400 value as −0.0217. That number is CPM 3400 (K 3,800); at K 3,400 (CPM 3000) f 0.87 gives −0.0177. The conclusion stands: f 0.87
  fails to retrodict ku30 (+0.0054). reach701_f061.json (6da80c1f, written before ku38's read): K 3,400 +0.00422 / +0.00316 / +0.00646,
  mean +0.00462, se 0.00097 → PASS; K 4,200 −0.00123 / −0.00172 / −0.00193, mean −0.00163 → FAIL.

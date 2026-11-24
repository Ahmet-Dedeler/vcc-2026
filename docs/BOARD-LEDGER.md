# Board ledger — every scored arm of ours, all six members

## 2026-10-22 — X2 = R0 + 811 (C0), the REPLACEMENT FINAL ENTRY: S5RxOrm4ZbqtnPRIX8NX, published, members hidden. This is the entry that counts unless a later one replaces it.
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| X2 = R0's day recipe (+ relaxed CD4) + DEF_DAY_DUMP_HOOK=dumphook811c0.sh (811 c0_l25, no swap, no weights move), re-gated with DEF_DAY_ALLOW_LOW_DOSE=1 | Oct 22 18:20 (upload 17:49–18:08) | hidden | hidden | hidden | hidden | hidden | hidden | hidden | S5RxOrm4ZbqtnPRIX8NX |
- Build: CONVERGED at lvl3 (D 1.0757 / E 0.8924 / F 1.0150; rho 0.9767 / 0.9757 / 0.9781). Package e39d2ba6, prediction identity 94ae91b6. gate6a PASS (judge 360f3125).
  554 re-gate (D prior_lambda 0.8994). Chosen by SIXA-DECISION + Amendments 1–5 via decide6a. Evidence for 811's increment: A/B/C board C0 − R-ABC2 +0.0041
  (ABOVE); D/E/F folds C0 − B0 +0.00446 [+0.00175, +0.00733] 3/3; sealed reserve LZ +0.00999 [+0.00659, +0.01273] 6/6. X1 (weights + swap) built and passed but was not
  submitted (one replacement; X2 ranks first). Name from the queue; the reserved final name is unused.

## 2026-10-22 — R0, the FIRST FINAL ENTRY (D/E/F): k2OKr0zPvJq48BCtdUx8, published, members hidden until the Nov 19 announcement
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| R0 = the frozen October recipe (decisions 04be2af3) on D/E/F (MOLM13 / U87MG / SKLMS1; NEW-PANEL) + DEF_DAY_CD4_THRESHOLDS=relaxed (287/300 map), re-gated with DEF_DAY_ALLOW_LOW_DOSE=1 | Oct 22 06:05 (upload 05:38–05:52) | hidden | hidden | hidden | hidden | hidden | hidden | hidden | k2OKr0zPvJq48BCtdUx8 |
- Build: def_day converged at lvl3, factors D 1.0737 / E 0.9063 / F 1.0294, rho 0.9762 / 0.9764 / 0.9749. lvl0 (E 1.5555) and lvl1 (D 1.0151, F 1.0800) were past the cliff and
  deleted. Package 5c82ea9e, prediction identity fa832988. Postbuild 42/42 with 554's dose override (D prior_lambda 0.8975 vs floor 0.9, a preregistered 0.85–0.9 case);
  zero-gene 0; FAILSAFE 0. is_final true, partition test. Only the LAST final entry before Nov 5 23:59Z counts. A 6a replacement (SIXA-DECISION) may follow.

## 2026-10-21 — W1 0.2746084 (D +0.0001410 vs pred 0.2744674): H1 ×2.5 on R0's exact recipe COMPOSES; no block on a 708 H1 move
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| W1 = R-ABC2's recipe + DEF_DAY_H1_STACK=x2.5 on def_day / A/B/C (673; re-gated by the coordinator, regate673.py), vs pred | Oct 21 00:53 (upload 00:01–~00:25) | 0.2746084 | 0.80956 | 0.38416 | 0.18453 | 0.00551 | 0.23542 | 0.02847 | WULHAm2fytPtmAb3ryVN |
| Δ vs R-ABC2 (nnHQfcZH) | | +0.0006235 | +0.00037 | +0.00161 | +0.00028 | +0.00009 | +0.00139 | +0.00001 | |
- Reader `read695.py --w1` (rc 0): pred = R-ABC2 0.2739849 + 710's +0.0004824 = 0.2744674; **D = +0.0001410 → COMPOSES**. **No block:** a 708 move of
  H1 above ×2 stays possible (SIXA Amendment 2 does not fire). w1route.sh read it at 00:53:43Z and left wait708 (R710=POSITIVE) alone, as designed.
- Bulk (Δpds + Δmse)/6 +0.000330 (on ALLSW's base, 710 gave +0.000380: the H1 step composes on R0's recipe). DE +0.000294 (raw reach +0.0014): re-roll class.
- Class: transfer (composition weights), guard only. Tries for the H1 dose on the board: ×1, ×2 (many), ×2.5 on ALLSW (710), ×2.5 on R0's recipe (this).
  W2 (CD4 ×1.5) was not built, so SIXA Amendment 3 never applies.

## 2026-10-20 — C0 0.2780624 (D +0.0016928 vs pred 0.2763696): R0 + 811's c0_l25 without the swap reads ABOVE; the swap-OFF X2 has its A/B/C read
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| C0 = R-ABC2's recipe with DEF_DAY_DUMP_HOOK=dumphook811c0.sh (811 λ 0.25 on the unswapped dump) on def_day / A/B/C (699), vs pred | Oct 20 07:29 (upload 07:04–~07:20) | 0.2780624 | 0.80892 | 0.40569 | 0.18525 | 0.00543 | 0.23464 | 0.02845 | 4rq3qHIuv8bHAxpTkUeu |
| Δ vs R-ABC2 (nnHQfcZH) | | +0.0040775 | −0.00027 | +0.02314 | +0.00099 | +0.00001 | +0.00060 | −0.00001 | |
- Reader `read695.py --c0` (rc 0): pred = R-ABC2 0.2739849 + ABCp's C0 − B0 +0.00269 ÷ ρ 1.1280 = 0.2763696; **D = +0.0016928 → ABOVE**. Not HARM
  (> 0.2758696; > 0.2734849). **Consequence (SIXA Amendment 1):** X2 on the swap-OFF branch has its A/B/C read. Both X2 branches now have one.
- Bulk (Δpds + Δmse)/6 = +0.003812, nearly all mse (+0.0231 scaled). DE members +0.000266 (raw reach +0.0005): about 1 re-roll sd, noise class.
- **The interaction (information, 4 board points on def_day):** swap alone (R-ABC2H − R-ABC2) +0.0009087; 811 alone (C0 − R-ABC2) **+0.0040775**;
  both (C1 − R-ABC2) +0.0037800. **C1 − C0 = −0.0002975**: once 811 is in, the swap adds nothing on A/B/C (−1.2 re-roll sd; pds −0.0027,
  mse +0.0008). 811 on top of the swap = 70 % of its no-swap size on the board (ABCp said ~55 %). Both move the strong targets toward the
  sources' mean, so they overlap. The 6a rule is unchanged: whether the swap ships is 708's call on D/E/F's own folds, and it is not
  re-decided from A/B/C reads (11.1c.3). The DEFo dossier's C1 − B0 vs C0 − B0 pair is where the same question gets a D/E/F answer.
- **Instrument (second 811 point):** ABCp C0 − B0 +0.00269 → board +0.00408 (1.52×); C1 − B1 +0.00153 → board +0.00287 (1.88×). Both 811
  increments are under-read by ABCp, not over-read. Class-median ρ 1.128 still predicted 59 % of it.
- Class: transfer (composition rows). Tries for 811 on the board: 2 (C1, C0). **Team's highest board overall** (prev. C1 0.2777649).

## 2026-10-20 — C1 0.2777649 (D +0.0019835 vs pred 0.2757814): R0 + swap + 811's c0_l25 reads ABOVE, so 811's 6a route stays open
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| C1 = R-ABC2H's recipe with DEF_DAY_DUMP_HOOK=dumphook811.sh (swap β 0.25, then 811 λ 0.25) on def_day / A/B/C (698), vs pred | Oct 20 00:35 (upload 00:01–~00:30) | 0.2777649 | 0.80625 | 0.40650 | 0.18572 | 0.00555 | 0.23409 | 0.02848 | RMxJtslBpY2GkwIMHOiB |
| Δ vs R-ABC2H (ndt8swt7) | | +0.0028713 | −0.00059 | +0.01773 | +0.00012 | −0.00001 | −0.00002 | +0.00000 | |
- Reader `read695.py --c1` (rc 0): pred = R-ABC2H 0.2748936 + ABCp's C1 − B1 +0.00153 ÷ ρ 1.7232 = 0.2757814; **D = +0.0019835 → ABOVE**.
  Not HARM on either line (> 0.2752814; > R-ABC2H − 0.0005 = 0.2743936). **Consequence:** input 3 of SIXA-DECISION is filled (C1 not HARM);
  X2 still needs the lab's DEFo/DEF dossier on the branch that ships and the LZ sign check. Readiness recorded the read in 698/RESULTS.
- The whole gain is bulk: (Δpds + Δmse)/6 = +0.0028564, nearly all mse (+0.0177 scaled; pds −0.0006). DE members +0.000015 (re-roll class:
  the hook changes only the stack dump, so pass 1 is redrawn with the same law). Raw reach −0.000018.
- **Instrument finding (one read):** ABCp said +0.00153 for C1 − B1 and the board gave +0.00287, 1.9× ABCp and 3.2× the deflated
  prediction. The lab's ρ 1.72 (fresh panels over-read direction) does not hold for this increment; on 707 and 710 the folds also
  predicted pds losses the board did not show. One read: it counts against the deflation here, it does not re-calibrate ABCp.
- Class: transfer (composition rows). Tries for 811 on the board: 1 (this one). Team's highest board overall (prev. sw707 0.2771181).

## 2026-10-19 — R-ABC2H 0.2748936 (D −0.0001244 vs pred 0.2750179): R0 + 707's swap COMPOSES on def_day, so the 6a route stands
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| R-ABC2H = R-ABC2's recipe + DEF_DAY_DUMP_HOOK=dumphook707.sh on def_day / A/B/C (697), vs pred | Oct 19 01:49 (upload 01:21–~01:45) | 0.2748936 | 0.80684 | 0.38878 | 0.18560 | 0.00556 | 0.23411 | 0.02848 | ndt8swt7ILYn3qMJYZ9a |
- Reader `read695.py --rabc2h`: pred = R-ABC2 0.2739849 + 707's +0.0010330 = 0.2750179; **D = −0.0001244 → COMPOSES** (±0.0005). The swap on
  top of R0's exact recipe adds +0.0009087 on A/B/C (vs +0.0010330 on ALLSW's base). **Consequence:** the 6a route (R0 + swap via 708's day-of
  rule) has its composition read and its def_day set-value rehearsal on A/B/C (and rider690-H on the new panel). It ships only if 708 turns the swap
  ON from D/E/F's own folds. C1's prediction is now filled in by the lab (333dddc8): 0.2757814; HARM if C1 ≤ 0.2752814 or < 0.2743936.

## 2026-10-19 — h1s710 0.2765675 (+0.0004824 vs ALLSW): the H1 dose step reads POSITIVE, so 708 gets no block on H1 > ×2
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| h1s710 = ALLSW with the stack re-dumped at VCC_H1_UNION_G=0.70 (re-levelled), vs ALLSW | Oct 19 01:20 (upload 00:19–~01:00; name renamed after a content-filter refusal) | 0.2765675 | 0.81530 | 0.39814 | 0.17619 | 0.00562 | 0.23572 | 0.02844 | UW4LsxTUkeLumDCBotfU |
- Reader `experiments/710-h1-dose-board/read710.py`. Primary Δoverall **+0.0004824** ≥ +0.0003 → **POSITIVE: no block** (R710=POSITIVE for decide708);
  the instrument's H1 step reads as calibrated on Arc's panel. Bulk (pds + mse)/6 +0.000380 (mse +0.00225 scaled; pds ~0); DE members
  +0.000103 (re-roll). Prediction +0.0002 (−0.0005 … +0.0008): inside, above the point. **October:** nothing changes in R0 (H1 ×2 DECIDED).
  708's day-of rule may still move H1 on D/E/F's own folds (its A/B/C day path recommended no move). 710 cannot block it, and its own
  production-key confirmation still applies. Class: transfer (composition weights). Tries: board H1 points ×1, ×2, and this one.

## 2026-10-18 — ku38 0.2736785 (+0.0008526 vs 446bj2z): the dial test refutes f 0.61; October K stays 2,600
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| ku38 = 446bj2z + VCC_RC_CALLSET_CPM=3800 (K 4,200), under OVERRIDE-554, vs 446bj2z | Oct 18 03:15 (upload 01:40–~03:00) | 0.2736785 | 0.81523 | 0.38042 | 0.17366 | 0.00714 | 0.24386 | 0.02177 | aiUKL3IOpD4miuD5LpEu |
- Reader `experiments/711-kup-dose/read711.py`. Falsifiers clean (Δpds 0, Δmse −2.2e-8, |Δnmae| 0.00044). **Primary Δreach raw +0.009424** (ku30 +0.005398):
  above Amendment 2's +0.0019, so **the f 0.61 dial is REFUTED**, condition (i) at f 0.87 stands, and reach701 says it fails. Δoverall: K 3,400
  +0.0011785, K 4,200 +0.0008526. The best size is 3,400, its condition fails, and the step-down ends at 2,600 → **`rc_callset_cpm` stays off**.
  R-ABC3 is not built. Raw DE vs 446bj2z: reach +0.0094, fid +0.0005, jac −0.0025, nmae +0.0004.
- **Finding:** the board's reach rises monotonically with K (2,600 → 3,400 → 4,200: 0 / +0.0054 / +0.0094 raw), as the mechanism predicted.
  BOTH H1 dials got the direction wrong (f 0.87 says it falls everywhere; f 0.61 said 4,200 falls back), so the H1 harness cannot price
  call-set size. The overall optimum on the board is ~3,400 (jac pays above it), but it has no held-out confirmation, so per §11.1c it
  stays out of October. Tries on size: 5 board points (1,600 / 2,000 / 2,600 / 3,400 / 4,200).

## 2026-10-18 — R-ABC2 0.2739849 (D −0.0000596 vs pred 0.2740445): the October recipe + θ 0.25 COMPOSES on def_day
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| R-ABC2 = R-ABC + DEF_DAY_LODO_STRENGTH=0.25 on def_day / A/B/C (695), vs pred | Oct 18 01:40 (upload 00:01–~01:27) | 0.2739849 | 0.80919 | 0.38255 | 0.18426 | 0.00542 | 0.23403 | 0.02846 | nnHQfcZHsT3as9EafA3i |
- Reader `experiments/670-coord-1015/read695.py --s705 o7ObvpGu`: pred = R-ABC 0.2743318 + 705's Δ −0.0002874 = 0.2740445; **D = −0.0000596
  → COMPOSES** (±0.0005). θ 0.25 stays DECIDED for R0; nothing reverts. Readiness's split (695 RESULTS d336db8f): Δ vs R-ABC −0.0003469; bulk
  (pds + mse)/6 −0.000322 vs 705's −0.000240 on ALLSW (exact members: a small real interaction of −0.00008 with R-ABC's base); DE −0.000025
  (re-roll). R0's recipe is now board-read on def_day as built. Frozen from this read: 697 R-ABC2H pred 0.2750179 (4ca54b40); 696 R-ABC3
  pred 0.2751634 if K is adopted (1e995bf1). The displayed row is now 0.2740 (rank 83); the best entry stays sw707 0.2771181.

## 2026-10-17 — sw707 0.2771181 (+0.0010330 vs ALLSW, NEW BEST): the strong-tier swap reads WEAK on its primary (Δ(pds+mse) +0.00084)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| sw707 = ALLSW + 707's strong-tier swap β 0.25 (VCC_BULK_LOAD …sw025, re-levelled), vs ALLSW | Oct 17 01:22 (upload 00:29–01:21) | 0.2771181 | 0.81415 | 0.40206 | 0.17618 | 0.00566 | 0.23620 | 0.02846 | AYbAiHRjuHdil2zCjlsk |
- Reader `experiments/707-strong-tier-direction/read707.py`. Primary Δ(pds + mse) in overall units = **+0.000842** (scaled pds −0.001123, mse +0.006174;
  raw pds −0.000472, mse −0.005871): bin **WEAK** (+0.0005 … +0.0010) → **not adopted into R0** (below the bar for a late, code-bearing switch).
  DE members +0.000191 overall (re-roll; falsifier clean). The prediction was −0.0004 … +0.0001 (best-of-N corrected): the read beat it.
  **The swap is still live for October through 708's day-of route** (AMENDMENT 2's relaxed rule: not HARM on A/B/C, D/E/F-panel folds decide;
  the 10-16 day-path rehearsal predicts s|r +0.00157, 3/3, so it would turn ON as a 6a replacement entry if D/E/F prices like A/B/C).
  Class: transfer (direction). Instruments: 062 folds called NULL, the lab's Arc-drawn ABCp / LABC called positive (7a38fb29). The lab grades them.
  Tries: 1 board arm (β 0.25; 713 found β 0.125 / 0.5 worse on LABC).

## 2026-10-17 — lstr705 0.2757977 (−0.0002874 vs ALLSW): the θ 0.25 LODO strength rider is CLEARED by its SAFETY read (primary −0.00024 > −0.0005)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| lstr705 = ALLSW + 522's θ 0.25 rider (VCC_BULK_LOAD …s025, re-levelled r3), vs ALLSW | Oct 17 00:28 (upload 00:01–00:28) | 0.2757977 | 0.81441 | 0.39530 | 0.17588 | 0.00568 | 0.23510 | 0.02841 | o7ObvpGu05zvmCwsFtiC |
| ALLSW (base) | Oct 14 | 0.2760851 | 0.81527 | 0.39588 | 0.17599 | 0.00558 | 0.23536 | 0.02843 | 4WNmsNfdQhYx6bfib0dR |
- Reader `experiments/705-lodo-strength-safety/read705.py`. Primary Δ(pds + mse) in overall units = **−0.000240** (scaled pds −0.000855, mse −0.000583;
  raw pds −0.000308, mse +0.000546): bin **CLEARED** (only ≤ −0.0005 vetoes). Prediction +0.0001 (−0.0003 … +0.0004): inside, at the low end.
  DE members −0.0000476 overall (re-roll only; falsifier |DE| > 0.0015 clean). **Consequence:** `lodo_strength` θ 0.25 ships for October via
  DEF_DAY_LODO_STRENGTH=0.25 (530, env-only). It enters R-ABC2 (10-18 s1) and rider690. Class: transfer (a key change). Its October value:
  fresh panels +0.7 … +1.9 % deflated by the key over-read (÷ ~2.5–3.1) ≈ +0.2 … +0.6 %; LABC (Arc-drawn) +0.00005; this A/B/C board −0.00029
  overall. So ≈ 0, a near-free option, as the lab predicted. Tries: 1 board arm on θ (522 chose θ among 4 on its folds).

## 2026-10-16 — K-UP 0.2740045 (+0.0011785 vs 446bj2z): call-set size above 2,200 reads NULL on its primary (Δ(fid+jac) +0.00012); reach +0.0054 raw unregistered
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| ku30 = 446bj2z + VCC_RC_CALLSET_CPM=3000 (K 3,400), vs 446bj2z | Oct 16 04:59 (sub. 04:34) | 0.2740045 | 0.81523 | 0.38042 | 0.17440 | 0.00857 | 0.23933 | 0.02608 | HAbelEdgYYziqwfjTgdI |
- Primary Δ(fid + jac) = +0.000123 overall (fid +0.000946, jac −0.000881 raw): bin **NULL** → `rc_callset_cpm` stays off (K 2,600). Falsifiers clean (pds 0,
  mse −5.9e-8 raw; fid up / jac down as a larger set should; nmae +6.6e-5 raw). **Not preregistered:** reach +0.005398 raw (+0.00108 overall), ~4.1× the
  same-law re-roll rms (0.0013 raw); K16 / K20 (size DOWN) left reach flat, so the effect is asymmetric. Adoption needs a preregistered replication
  (a pass-1 re-draw at CPM 3,000 with reach as primary), proposed for 10-18 s2. Dose curve now K 1,600 / 2,000 / 2,600 / 3,400.

## 2026-10-16 — R-ABC 0.2743318 (+0.0012429 vs 524dr2): the resolved October recipe COMPOSES on def_day (95 % of the predicted +0.0013033)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| R-ABC (694) = the resolved October recipe on def_day on A/B/C (524dr2 + UPRN + DNMAG v2), vs 524dr2 | Oct 16 04:21 (sub. 03:55) | 0.2743318 | 0.81077 | 0.38291 | 0.18412 | 0.00546 | 0.23430 | 0.02844 | PooEbFClgTLMQoXiOCQP |
- pds, mse, fid, jac raw Δ exactly 0; nmae −0.004458 raw (scaled +0.00747); reach −1.1e-5 raw. Bin COMPOSES [+0.0010, +0.0016]. UPRN + DNMAG add on the
  October path (tercilec + CLEAN4 hook) as they did on the 446bj2z path (+0.00093 + +0.00037). gate694's config-diff row was overridden (path-only
  diffs; gate v2 then passed the same package). Rank 74 (displayed row). Record: experiments/694-rabc/ (OVERRIDE-gate694.md, RESULTS).

## 2026-10-15 — 702 (PIE bulk blend) 0.2735387 (−0.0025464 vs ALLSW): HARM, both bulk members worse; October `dump_hook` off
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 702 = ALLSW + VCC_BULK_LOAD → PIE-blended cap-6 dump (w 0.25, seen targets only, norm-preserving), re-levelled, vs ALLSW | Oct 15 04:35 (sub. 04:01) | 0.2735387 | 0.80621 | 0.39077 | 0.17590 | 0.00562 | 0.23431 | 0.02843 | U7BlIZ53wTlLbaZTTET2 |
- Primary Δ(pds + mse) −0.0023631 overall (pds −0.00419 raw, mse +0.00488 raw, i.e. worse); DE members −0.00018. Predicted +0.0008 … +0.0017
  (Δmse an improvement of 0.0013 … 0.0031). Bin **HARM**: not adopted. The held-out instruments (K562 GWPS seen targets, h1-test seen) over-read by
  ~+0.004 overall and got mse's SIGN wrong. A PIE-citing team's gain (0.2524 → 0.2993) came on a weaker base. Reader: experiments/702-pie-bulk-arm/read702.py.

## 2026-10-15 — DN 0.2741292 (+0.0003736 vs 630arm): the DOWN magnitude pass reads PARTIAL (Δnmae raw −0.001646); October adopts d_m200, flagged
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| DN = 630arm + VCC_DE_DNMAG=d_m200 (631), vs 630arm | Oct 15 01:49 (sub. 01:12, upload 20 min) | 0.2741292 | 0.81523 | 0.38042 | 0.18234 | 0.00550 | 0.23287 | 0.02842 | wzaVLFTmlAMW3Mfw0dD8 |
- Five members bit-identical (reach +9.6e-8 raw); falsifiers clean; build fail-safes 0. Prediction −0.0018 … −0.0048 raw: read at 0.55× the centre.
  Record: `experiments/671-dnmag-arm/BOARD-RESULT-DN.md`.

## 2026-10-14 — ALLSW 0.2760851 NEW BEST (+0.0011661 vs 498p2g; R1 PASS, R2 ADDITIVE: the DE switches add exactly on the cap6+P2g bulk)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| ALLSW = 498p2g + zs + CLEAN4 hook 0.05 (no re-level), vs 498p2g | Oct 14 23:04 (sub. 22:16, upload 35 min) | 0.2760851 | 0.81527 | 0.39588 | 0.17599 | 0.00558 | 0.23536 | 0.02843 | 4WNmsNfdQhYx6bfib0dR |
- I = +0.0000629 overall (ADDITIVE); R1 v2 +0.0003202 ≤ 0.0006 PASS; prediction 0.2760235 (read +0.0000616). F1 fired on Δpds raw −3.716e-6
  = exactly one rank flip (1/(3·300·299)) from the 28 documented ±1 rounding swaps: diagnosed benign. Rank 59.
  Record: `experiments/658-prereg-1007-builds/BOARD-RESULT-ALLSW.md`.

## 2026-10-07/08 (recorded 10-14 22:30Z, after a 7-day idle gap) — CPH 0.2734755 (+0.0006496 vs 446bj2z; callset_block does NOT ship); 459dr 0.2673086; 630arm 0.2737556 (+0.0009297 vs 446bj2z: UPRN PARTIAL); 524dr2 0.2730889 (+0.0006835 vs 524dr)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 446bj2zcph = 446bj2zcp + CLEAN4 hook (call-propensity block + hook), vs 446bj2zcp / 446bj2z | Oct 7 (sub. 00:01) | 0.2734755 | 0.81523 | 0.38042 | 0.17760 | −0.00225 | 0.23059 | 0.03926 | FRGEXV8WiqxcX5Es6SyI |
| 459dr = October path with the H1 ×3 axis, q 1.75, vs 453dr | Oct 7 (sub. 11:24) | 0.2673086 | 0.81089 | 0.35349 | 0.17225 | 0.00543 | 0.23340 | 0.02839 | 8F3f1Ay3vwrtFdMbjgTU |
| 630arm = 446bj2z + VCC_DE_UPRN (rank-neutral UP magnitude, K 1.0), vs 446bj2z | Oct 8 (sub. 00:01) | 0.2737556 | 0.81523 | 0.38042 | 0.18010 | 0.00550 | 0.23287 | 0.02842 | QJRJK6OC6Zw1kjwnZqid |
| 524dr2 = the resolved October recipe (524dr + tercilec + mde_sign), vs 524dr | Oct 8 (sub. 00:39) | 0.2730889 | 0.81077 | 0.38291 | 0.17665 | 0.00546 | 0.23432 | 0.02844 | dVHLws6KdMxudq0Yn5ap |
- **CPH:** primary Δ(fid + nmae) vs 446bj2zcp **+0.00120** overall (fid +0.00124, nmae −0.00199 raw; pds/mse bit-identical; jac +0.000003 raw):
  POSITIVE on its own bin (≥ +0.0004), so the hook repairs the block. But the prereg's fired branch ("the read vs 446bj2z decides
  whether the block ships at all") uses 535's bin: Δ(fid + jac) vs 446bj2z = **+0.00052 < +0.0012**, so **the block does not ship**:
  `callset_block` → off (to be decided in october-decisions.json). vs 446bj2z: fid −0.00234, jac +0.00399, reach −0.00200, nmae −0.00189 raw.
- **630arm:** Δnmae raw **−0.00305**, the other five members bit-identical: **PARTIAL** (−0.0040 … −0.0015). U075 was never run, so
  the dose rule cannot read; per the bin "October adopts (mechanism-only, nothing fitted), flagged". H1's centre −0.005 … −0.011 over-read.
- **524dr2 vs 524dr:** +0.00068 (nmae −0.00173, reach +0.00127 raw, bulk ~0): tercilec + mde_sign compose on def_day as expected.
- **459dr:** +0.00095 vs 453dr (mse −0.0039, reach +0.0019, pds −0.0005 raw) but −0.0051 vs 524dr (mse +0.028 raw): the H1 ×3 axis
  is not October material next to the decided recipe.
- **Context (10-14 22:13Z live board):** these are ±0.001 reads while the field moved by 0.03–0.05. 61 teams now sit above our best
  (0.2749190; 30 on 10-08), 19 above 0.30, Illumina AI 0.4729. The pack beats us on every member but jac (fid 0.04–0.09 scaled vs our
  0.005). See NEXT-SLOT-DECISION 10-14.

## 2026-10-06 — K16 0.2707192 (−0.0021068 vs 446bj2z): K 1,600 HURTS (fid −0.0041 raw); October K stays 2,200; 524dr 0.2724054 (+0.0060 vs 453dr = Σ + 0.0006): October recipe PASS
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 446bj2zk16 = 446bj2z + VCC_RC_CALLSET_CPM=1200 (K 1,600), vs 446bj2z | Oct 6 18:15 (sub. 05:31; scoring queued ~12.5 h) | 0.2707192 | 0.81523 | 0.38042 | 0.17544 | −0.00785 | 0.23334 | 0.02773 | Imj3EXViNco8rT4eLLlE |
| 524dr = the October recipe on def_day (LODO pct, tercile, zs; hook off), vs 453dr | Oct 6 23:04 (2nd slot; 1st entry cancelled stuck; reboot-resumed) | 0.2724054 | 0.81065 | 0.38287 | 0.17400 | 0.00553 | 0.23293 | 0.02845 | KTzIDR0Occ2jYzklkARE |
- Primary Δ(fid + jac) −0.00234 overall → bin ≤ −0.0005 (keep 2,200). With K20's −0.00086, the dose rule adopts nothing, so
  **rc_callset_cpm DECIDED off**. I-662 (Δfid −0.0116 predicted) FAILS; the read −0.0041 puts the board near f 0.87 for size changes.
  The 524dr slot-1 entry VqBgraiBxlaK0WP5efOF was cancelled while still queued (5.5 h, no machine; no slot spent).
  Record: `experiments/523-callset-size/BOARD-RESULT-446bj2zk16.md`.

## 2026-10-05 — 496a 0.2728318 (+0.0005552 vs 446bj2): the plain CLEAN4 uncalled-gene hook TRANSFERS (Δnmae −0.00158 raw); idall3 0.2746732 (−0.0002458 vs 498p2g): the identity hook is BELOW the plain hook
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 496a = 446bj2 + VCC_DE_MDE_SIGN=clean4:0.05 (the uncalled-gene sign hook), vs 446bj2 | Oct 5 01:03 (sub. 00:01, Wi-Fi; scoring ~40 min) | 0.2728318 | 0.81523 | 0.38042 | 0.17364 | 0.00542 | 0.23386 | 0.02842 | HoQY4xrzYer5y232oh85 |
| idall3 = 498p2g + the ungated v3 identity hook (604 artifact λ 0.5, no KEYGATE), vs 498p2g | Oct 5 01:37 (sub. 01:04 after 496a's capture) | 0.2746732 | 0.81528 | 0.39588 | 0.16943 | 0.00562 | 0.23340 | 0.02842 | 4R0FMEtcLtvZMmg2X0G5 |
- nmae −0.00158 raw (+0.00040 overall), reach +0.0009 raw (unpriced, +0.00016), fid/jac pinned, pds bit-identical, mse −4.6e-8
  raw (DE-only arms move the jackknife mse by O(1e-8); the falsifier is 1e-6). BIN TRANSFERS → October `mde_sign` = clean4:0.05,
  decided after recipe_dr524 resolves (658). 618's I-662 price (−0.0014) holds; h1-test's (−0.0022) over-reads 1.4×. CPH's T 0.797
  → CPH primary +0.0008 (P8), below SHIP. Record: `experiments/496-uncalled-sign-hook/BOARD-RESULT-496a.md`.
- idall3: nmae +0.00027 raw (worse than no hook), reach −0.0007, fid/jac pinned. Primary Δ(fid + nmae) −0.000111 overall, i.e.
  −0.000501 vs Δ_496a → BELOW (not HARM). October keeps the plain hook. h1-test (WIN +0.0007 … +0.0010), 618 (TIE +0.00007)
  and I-662 (+0.0002) all over-predicted; 618 was the closest. Record: `experiments/663-idall3-arm/BOARD-RESULT-idall3.md`.

## 2026-10-04 — 498p2g 0.2749190 NEW BEST (+0.0010083 vs 495; +0.0026424 vs 446bj2: cap 6 + P2g compose additively); K20 0.2721097 (−0.0007163 vs 446bj2z: fid falls)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 498p2g = 495 (cap-6 LODO key) + P2g per-target sharpening, re-levelled, vs 495 | Oct 4 05:37 (sub. 00:01, network drop 00:59–05:21) | 0.2749190 | 0.81528 | 0.39588 | — | — | — | — | JWSghc0o7MN5D6o8Iadc |
| 446bj2zk20 = 446bj2z + VCC_RC_CALLSET_CPM=1600 (K 2,000), vs 446bj2z | Oct 4 16:07 (hotspot, 143 resume tries) | 0.2721097 | 0.81523 | 0.38042 | 0.17467 | −0.00008 | 0.23358 | 0.02884 | kJfe2dP9LEJawDEf0pMm |
- K20: Δ(fid + jac) −0.00086 overall (fid −0.0017 raw), bulk exactly 0. The 610 card predicted +0.0044, wrong in SIGN: for
  call-set size the board is f 1-like. October K stays 2,200 pending K16. Record: `experiments/523-callset-size/BOARD-RESULT-446bj2zk20.md`.
- Increment vs 495: pds +0.0013, mse −0.0027 raw (bulk +0.00095 overall), DE flat. P2g HOLDS (rule fixed 05:25Z before the
  score), so C10b's F = M6+P2g. Record: `experiments/498-combined-bulk/BOARD-RESULT-498p2g.md`.

## 2026-10-03 — 495 0.2739107 NEW BEST (+0.0016340 vs 446bj2: the LODO key at cap 6 TRANSFERS); 471c 0.2696638 (−0.0003790 vs 446b: HARM bin, bulk flat)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 495 = 446bj2 with the LODO key at cap 6 (lodo6, hi6 load, re-levelled), vs 446bj2 | Oct 3 09:59 (sub. 00:01, Mac asleep 01:45–09:59) | 0.2739107 | 0.81239 | 0.39305 | 0.17018 | 0.00569 | 0.23369 | 0.02847 | 9Hxef2F6Clf4G2hIkNPq |
| 471c = 446b with the h1_2025_full union swap, re-levelled, vs 446b | Oct 3 13:19 (paused 10:53–12:08, hotspot) | 0.2696638 | 0.81467 | 0.36545 | 0.17155 | 0.00552 | 0.23240 | 0.02839 | yzNoMVih0pKqoDjdqrXm |
- mse −0.0120 raw (gain), pds −0.0013 raw (small cost); DE pinned and inside 537's 3σ noise bands. BIN TRANSFERS, so
  **October `lodo_hi_pct` DECIDED 96.64697522769474** (the percentile form; == cap 6 on A/B/C). The PLOLO instrument predicted
  +0.0019 … +0.0029: same sign, ~15–40 % over-predicted. Record: `experiments/495-lodo-cap6/BOARD-RESULT-495.md`.
- 471c: bulk ~0 (mse −0.0002, pds ~0 raw), reach −0.0024 raw inside 537's re-roll band. Bin HARM (≤ −0.0002); keep the default
  union. Record: `experiments/471-h1full-union/BOARD-RESULT-471c.md`.
- Delivery: the upload died at 01:37Z (DNS) while the Mac slept on battery with the lid closed; the 536 resume block
  finished it on wake.

## 2026-10-02 — 479p2g 0.2733131 NEW BEST (+0.0010364 vs 446bj2: per-target key-tercile sharpening TRANSFERS); 479sp2 0.2715625 (−0.0007141: strength terciles FAIL)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 479p2g = 446bj2 + per-target sharpening tercile:2.5,2.0,1.75 (key-ratio terciles), vs 446bj2 | Oct 2 00:27 | 0.2733131 | 0.81870 | 0.38436 | 0.17102 | 0.00551 | 0.23190 | 0.02839 | ioNMsGwrG2EvYuoB5Zu7 |
| 479sp2 = 446bj2 + stercile:2.5,2.0,1.75 (strength terciles), vs 446bj2 | Oct 2 00:55 | 0.2715625 | 0.81283 | 0.38073 | 0.16822 | 0.00550 | 0.23373 | 0.02837 | NPiJoPH47VYhiWQHhoOG |
- P2g: pds +0.0016, mse −0.0037 raw, DE in noise: TRANSFERS; S_P2: pds −0.0010, nmae +0.0018 raw (low-tercile ballast): FAILS.
  **October `sharpen_pt` DECIDED tercile:2.5,2.0,1.75**; carried-only cut points (493h) queued as an upgrade. Record:
  `experiments/479-pertarget-sharpen/BOARD-RESULT-479-pair.md`.

## 2026-10-01 — 460dr 0.2692306 (+0.0028694 vs 453dr: the October path's day-of LODO key transfers); 446bj2zcp 0.2720318 (−0.0007942 vs 446bj2z: the call-propensity block alone is NULL)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 460dr = 453dr's October path + `DEF_DAY_LODO=1` (day-of LODO key), vs 453dr | Oct 1 00:27 | 0.2692306 | 0.81433 | 0.36313 | 0.17145 | 0.00547 | 0.23262 | 0.02838 | LrkhsIXWMjEsvRXSGToa |
| 446bj2zcp = 446bj2z + 528's call-propensity DE block (`kept` steering), vs 446bj2z | Oct 1 01:01 | 0.2720318 | 0.81523 | 0.38042 | 0.17454 | −0.00640 | 0.22914 | 0.03926 | K1Q0y2OFYjTNRqh9jKfk |
- 460dr vs 453dr: mse −0.0131 raw (0.90× the A/B/C key's board effect), pds +0.0010 raw (5.0×; the A/B/C key's pds Δ was only
  +0.0002), DE inside noise (inside every 3σ band of 537's preregistered noise check). Bin "same sign, outside the band:
  keep it ON"; `DEF_DAY_LODO=1` stands. Record: `experiments/460-lodo-dayof/BOARD-RESULT-460dr.md`.
- 446bj2zcp vs 446bj2z: pds identical, mse −2.3e-8 raw (DE-only on the board to that precision), jac +0.0040, fid −0.0036, reach −0.0034 raw; primary
  Δ(fid + jac) −0.00018 overall: NULL (predicted +0.0020 … +0.0025: H1 over-read jac 2.3× and missed the fid loss). The
  swapped-in genes are called but unsteered. 537's rescue condition met: 446bj2zcph (block + CLEAN4 hook, 549, READY) is the
  10-07 candidate. Record: `experiments/535-callprop-arm/BOARD-RESULT-446bj2zcp.md`.

## 2026-09-30 — 458d 0.2615288 (−0.0085 vs 446b: Arc's H1 programme HARMS as the DE block prior); 458g 0.2650949 (the gate recovers +0.0036)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 458d = 446b with Arc's H1 programme (all public H1 cells) as every context's DE block prior | Sep 30 06:01 | 0.2615288 | 0.81472 | 0.36520 | 0.16887 | −0.01267 | 0.20478 | 0.02827 | b7UI3NRy1rdVNRRyZwmU |
| 458g = 458d with the design-strength gate q33 (weak third back to CLEAN4) | Sep 30 06:35 | 0.2650949 | 0.81472 | 0.36520 | 0.17017 | −0.00800 | 0.22012 | 0.02836 | JL8llNLHaTyDQzm798TV |
- 458d vs 446b: fid −0.0055 raw (HARM), reach −0.0267, nmae +0.0015, pds/mse bit-identical, no jac leak. 458g vs 458d +0.0036
  (the gate helps) but 458g vs 446b −0.0049. **October: no programme** (`programme_source off`, `prior_gate off`).
  Slot 1's upload died when the Mac slept on battery (02:39Z) and published after its retry (06:01Z). Record:
  `experiments/458-arc-lab-programme/BOARD-RESULT-458d-458g.md`.

## 2026-09-29 — 457dr 0.2695096 (+0.0031484 vs 453dr: the CD4 half transfers on the October path; October keeps ×1); 446bd2 0.2704629 (+0.0004200 vs 446b: DLD-1 flat, stays out)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 457dr = 453dr's October path with the stack at CD4 ×1.5, vs 453dr | Sep 29 01:50 | 0.2695096 | 0.81556 | 0.36445 | 0.17129 | 0.00534 | 0.23205 | 0.02837 | 1FfirZsGsAZvaXYTd39J |
| 446bd2 = 446b + DLD-1 as the bulk's second union residual (g 0.2), vs 446b | Sep 29 10:52 | 0.2704629 | 0.81963 | 0.36570 | 0.17111 | 0.00550 | 0.23246 | 0.02838 | R5uE5Rve3R9nWXGbeSLM |
- Bin ≥ +0.0010 (band +0.0005 … +0.0035, predicted +0.0017 … +0.0022): all bulk (mse −0.0143 raw, pds +0.0017), DE pinned.
  `cd4_stack` stays ×1 (re-open threshold +0.0040). With 453dr, the October path at production's CD4 dose = 446b −0.0005
  (bulk +0.0000): the path's whole A/B/C bulk cost was the CD4 half. Record: `experiments/457-cd4-half-def/BOARD-RESULT-457dr.md`.
- 446bd2 vs 446b: pds +0.0021 raw (top of its band), mse −0.0005 raw (below the −0.001 … −0.004 prediction), reach −0.0023
  raw (realisation noise): overall +0.0004, bin FLAT, DLD-1 stays out of the bulk and the October stack. Record:
  `experiments/445-dld1-bulk-union/BOARD-RESULT-446bd2.md`.
- The server was slow: slot 1 published 1 h 49 min after upload start; slot 2's upload (446bd2, from 03:08Z) crawled
  (148 MB of 2.77 GB after 3.6 h, ~0.7 MB/s when re-checked).

## 2026-09-28 — 446bj2z NEW BEST 0.2728260 (+0.0005493 vs 446bj2: zs transfers, pure nmae); 453dr 0.2663613 (the October path, VALIDATED)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 453dr = the October path (`def_day.sh`) on the relabelled A/B/C bundle, vs 446b | Sep 28 12:37 | 0.2663613 | 0.81201 | 0.34939 | 0.17167 | 0.00546 | 0.23127 | 0.02836 | q3VuTImh4FRGKXH2waCp |
| **446bj2z = 446bj2 + `VCC_DE_MAG=dn_k250_zs`** | Sep 28 13:05 | **0.2728260** | 0.81523 | 0.38042 | **0.17452** | 0.00550 | 0.23287 | 0.02842 | be90smZUDcIV7j5MgQ03 |
- Both slots fired 12 h late: the 09-27 ~20:23Z reboot killed the 09-28 waiters; re-armed 12:11Z with the same arguments.
- 453dr vs 446b: bulk −0.00309 (preregistered band −0.0045 … −0.0005: **VALIDATED**), DE −0.00059 (pinned), no hard
  failure. The October pipeline costs exactly its priced amount on A/B/C (mostly CD4 ×1 vs ×1.5, on mse). Record:
  `experiments/453-abc-dress-rehearsal/BOARD-RESULT-453dr.md`.
- 446bj2z vs 446bj2: nmae −0.002179 raw (predicted −0.0013 … −0.0029), pds/mse/fid/reach/jac **bit-identical**: zs
  TRANSFERS, October `de_mag = dn_k250_zs`. Record: `experiments/456-zs-riders/BOARD-RESULT-446bj2z.md`.

## 2026-09-27 — 446bj2 NEW BEST 0.2722766 (+0.0022338 vs 446b): the LODO key transfers
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| **446bj2 = 446b with the stack under the LODO key (lo 0.5, hi 3), re-levelled** | Sep 27 00:27 | **0.2722766** | 0.81523 | **0.38042** | 0.17122 | 0.00550 | 0.23287 | 0.02842 | YvVru3DAF89Jiyr34olR |
- mse −0.0145 raw (past the predicted −0.004 … −0.012); pds +0.0002; DE members pinned. key_ok, so October runs
  `DEF_DAY_LODO=1`. Record: `experiments/446-q175-level/BOARD-RESULT-446bj2.md`.
| 446br2 = 446b with the repair floor −2 | Sep 27 01:04 | 0.2704006 | 0.81474 | 0.36513 | 0.17356 | 0.00557 | 0.23503 | 0.02837 | wE7usCaR6TAlPp16TR0V |
- 446br2: nmae −0.0014 raw (partial bin), calls and bulk pinned, overall +0.0004 < +0.0006. **−5 stays.** Record:
  `experiments/446-q175-level/BOARD-RESULT-446br2.md`.

## 2026-09-26 — 446b NEW BEST 0.2700428 (+0.0022834 vs 433k); then 446bja 0.2689767 (−0.0010661 vs 446b: the same-line programme does not transfer)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| **446b = 433k, q 1.5 → 1.75 at matched ρ** | Sep 26 00:37 | **0.2700428** | **0.81472** | **0.36520** | 0.17132 | 0.00557 | 0.23507 | 0.02837 | UDgjOKL9gI1udAexq6ue |
| 446bja = 446b + A's block prior from Nadig Jurkat (`clean4q80ja`) | Sep 26 01:03 | 0.2689767 | 0.81472 | 0.36520 | 0.17273 | 0.00624 | 0.22660 | 0.02837 | AsHteaDAXCWQo57fLYti |
- Raw pds +0.0032 and mse −0.0076, both in their transfer bins. 433k's record predicted +0.0021 for this step; the board
  gave +0.0023.
- nmae +0.0013 raw (−0.0004 overall) sits just outside the "pinned" bin. It is the ballast: median K 4 → 6 for the
  sharper render, so the carriers are re-drawn on fewer cells. q stops at 1.75.
- Record: `experiments/446-q175-level/BOARD-RESULT-446b.md`.
- 446bja: pds/mse exactly pinned, fid +0.0002 raw (predicted +0.0015 … +0.0040: no transfer), nmae −0.0010 (neutral), reach
  −0.0079 (the small-gate cost H1 showed). October keeps `noprog`. Record: `experiments/450-same-line-programme/BOARD-RESULT-446bja.md`.

## 2026-09-25 — 433h then 433k — two NEW BESTS: 0.2631120 (+0.0094597), then 0.2677594 (+0.0046475)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| **433h = 433a + zero repair −5 + 417h's re-level** | Sep 25 00:35 | **0.2631120** | 0.79726 | 0.33835 | **0.17503** | 0.00551 | 0.23409 | 0.02844 | wMoDAi1mca08FA2Bfqkd |
| **433k = 433h with q 1.5 + 417k's re-level** | Sep 25 01:02 | **0.2677594** | **0.80760** | **0.35717** | 0.17358 | 0.00553 | 0.23428 | 0.02838 | T76JWD3tOajk5IYvzwUW |
- 433h: the zero repair gave nmae −0.0247 raw, 3× its preregistered top bin: a handful of −35 cells in A/B/C gates had
  dominated whole targets. The re-level gave mse −0.0134 raw.
- 433k: pds +0.0046 and mse −0.0179 raw with the DE members pinned. The levelled q step transferred at 438's calibrated
  rate (+0.0044 predicted, +0.0046 read).
- Records: `experiments/433-coarse-mde-arm/BOARD-RESULT-433h.md`, `BOARD-RESULT-433k.md`.

## 2026-09-25 — 433k (433h + q 1.5 + its own re-level) — NEW BEST 0.2677594 (+0.0046 vs 433h)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| **433k = 433h, q 1.25 → 1.5 + 417k re-level** | Sep 25 01:02 | **0.2677594** | **0.80760** | **0.35717** | 0.17358 | 0.00553 | 0.23428 | 0.02838 | T76JWD3tOajk5IYvzwUW |
- Raw pds +0.0046 and mse −0.0179; DE flat. The q step transfers at k_sh ≈ 0.6.
- Record: `experiments/433-coarse-mde-arm/BOARD-RESULT-433k.md`.

## 2026-09-25 — 433h (433a + zero repair + re-level) — NEW BEST 0.2631120 (+0.0095 vs 433a)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| **433h = 433a + `VCC_DE_ZERO_REPAIR=-5` + re-level** | Sep 25 00:34 | **0.2631120** | 0.79726 | **0.33835** | **0.17503** | 0.00551 | 0.23409 | 0.02844 | wMoDAi1mca08FA2Bfqkd |
- Raw nmae −0.0247 (the zero-gene repair, 3× its predicted maximum); raw mse −0.0134 (the re-level, exactly its
  price); DE calls flat.
- Record: `experiments/433-coarse-mde-arm/BOARD-RESULT-433h.md`.

## 2026-09-24 — 433a (coarse MDE on 417e) — NEW BEST 0.2536523 (+0.0051 vs 417e)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| **433a = 417e + `VCC_DE_MDE=0.1,0.1,0.3,coarse`** | Sep 24 01:27 | **0.2536523** | 0.79613 | 0.32419 | **0.13311** | 0.00539 | 0.23467 | 0.02843 | FCbBcvGxF1HOFPxPISUC |
- Raw nmae −0.0168 (1.76× H1). pds and mse are bit-identical to 417e (16 digits). reach +0.0024. First entry above 0.25.
- Record: `experiments/433-coarse-mde-arm/BOARD-RESULT-433a.md`.

## 2026-09-24 — 417e (bulk stack H1 ×2 + CD4 ×1.5 via LOAD) — NEW BEST 0.2485853 (+0.0051 vs 417a)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| **417e = 417a + `VCC_BULK_LOAD` stack** | Sep 24 00:48 | **0.2485853** | 0.79613 | **0.32419** | 0.10532 | 0.00537 | 0.23202 | 0.02849 | wPP8U70kj4mOIvtA83CO |
- Raw mse −0.0262 and pds +0.0014 vs 417a; the DE members are pinned (|Δ| ≤ 0.0005).
- The stack transferred at k ≈ 1.3 of 420's PLOLO price (+0.0037 at k = 1): composition steps transfer at about full
  rate, shape steps (sharpening) at about half.
- Record: `experiments/417-bulk-sharpen/BOARD-RESULT-417e.md`.

## 2026-09-23 — 417a (bulk-only sharpening q 1.25 through ballast) — NEW BEST 0.2435334 (+0.0070 vs 410b)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| **417a = 410b + `VCC_BULK_SHARPEN_Q=1.25`** | Sep 23 00:27 | **0.2435334** | 0.79325 | **0.29657** | 0.10507 | 0.00523 | 0.23259 | 0.02848 | T88GsCjp53Bmk8xQcJMS |
| 417c = 417a + `VCC_CARRIER_GAMMA=0` | Sep 23 00:54 | 0.2418746 | 0.79391 | 0.29653 | **0.08769** | 0.00923 | 0.23532 | 0.02856 | 0xe3j7aqD2bypXKkCTBA |
- Raw pds +0.0064 and raw mse −0.0249 vs 410b. The DE members are pinned (|Δ| ≤ 1e-4 raw except reach +0.0011).
- Implied sharpening transfer is k_sh 0.42 (pds) / 0.65 (mse), consistent.
- Record: `experiments/417-bulk-sharpen/BOARD-RESULT-417a.md`.
- 417c's unwhitened carriers gained fid (+0.0012 raw) and reach (+0.0024) but lost nmae (+0.0089 raw), −0.0017 overall.
  H1 had predicted nmae −0.0068, so it mis-signed nmae for a carrier-design change for the second time (after 392).
  Record: `BOARD-RESULT-417c.md`.

## 2026-09-22 — 410 (welded per-target bulk key) and 410b (ballast) — 410b is the NEW BEST (+0.0126 vs 395)
| arm | UTC | overall | pds | mse | nmae | fid | reach | jac | entry |
|---|---|---:|---:|---:|---:|---:|---:|---:|---|
| 410 wshare_lo0.3_hi3 (welded) | Sep 22 | 0.2195152 | 0.77817 | **0.27123** | 0.01676 | 0.00263 | 0.21959 | 0.02871 | EETok0iqIbzQFg4CTYMQ |
| **410b = 410 + ballast (unwelded)** | Sep 22 | **0.2365722** | 0.77873 | **0.27040** | 0.10525 | 0.00520 | 0.23131 | 0.02854 | 4WgVi6FIE4iVw9PmRUP3 |
Raw mse 0.730213 (−0.0797 vs 395): the per-target key transfers; the weld cost nmae +0.054 raw. `experiments/410-pertscale-revisit/BOARD-RESULT-410.md`.

## September 7–8, 2026 — four previously omitted entries

Backfilled September 9 from the complete records under `data/board-entries/`.
These are **published scaled members**, not values inferred through a fitted
raw-to-scaled map. Scores displayed below are rounded; full precision and the
exact sum-of-six check are in `docs/review-2026-09-09/state.json`.

| arm | UTC date | overall | pds | mse | nmae | fid | reach | jac |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 312 plain −1.5 | Sep 7 | 0.217718885 | 0.790389028 | 0.180131402 | 0.102766825 | −0.015950255 | 0.220258776 | 0.028717532 |
| 312b plain −3.0 | Sep 7 | 0.216451143 | 0.792408225 | 0.176341809 | 0.099804062 | −0.016385172 | 0.217906494 | 0.028631440 |
| 316 PM parity | Sep 8 | 0.216973287 | 0.788136704 | 0.176144233 | 0.106153676 | −0.016860247 | 0.219531421 | 0.028733937 |
| **319 cap −1.5** | Sep 8 | **0.218659739** | 0.788313923 | 0.181012812 | 0.107514506 | −0.015759119 | 0.222191666 | 0.028684649 |

IDs, in row order: `Q519I2pLWRCW7PAEx5ul`, `1n5LFk6PSEHq15UHENIE`,
`2K5ppxkowMSj2ugoZqYb`, `xuDelEik7zGWMXdNu0R4`.
The novelty ledger now includes all four retained package and content hashes.
Their original captures had null package hashes: the new hashes are clearly
labelled retrospective local verification, not retrospective server attestation.

**319 is the current best recorded artifact.** Its gain is +0.000940855 over
312 and +0.002565337 over cd4svd5. It is not the largest historical improvement:
the CD4 dose change alone gained +0.0130557. Rank and score are different;
the September 9 check was rank 5 of 821 with this same best score.

**Historical sections below retain their original dates and conclusions.**
They do not certify that a missing package was unsubmitted, that a quoted
incumbent is still current, or that a fitted scaling relation is exact.

Built 2026-09-03 from the `*BOARD-RESULT.md` files in `experiments/`. **This is the only place
the full six-member trajectory exists in one table.** Every pricing model, retrodiction and gate
should be fit and checked against it.

Two facts that make it trustworthy:
- **The board scorer is deterministic.** A byte-identical resubmission of `a15-cd4f010.vcc`
  (entry `5QfX14axLAG8FF2rSDJe`) reproduced every member to float representation (`docs/LOG.md`,
  2026-08-31). There is no run-to-run scorer variance to subtract; every delta below is signal.
  Re-seeded twins of the *same config* are a different matter and do differ.
- The scaling map is exact against these rows. Fitted on 612 live board entries:
  `scaled = slope·(raw − b)` with `pds` +2.2138, `mse` −1.051, `nmae` −1.6532, `fid` +3.3645,
  `reach` +1.1292, `jac` +2.7263. Check: a8 raw reach 0.1890 → (0.1890−0.0790)×1.1292 = 0.12423
  against a recorded 0.124349. Anchors have not moved since August.

## Scaled members (these average to Overall)

| arm | date | overall | `pds` | `mse` | `nmae` | `fid` | `reach` | `jac` |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| a2-prune | 08-26 | 0.1503084 | 0.736561 | 0.000000 | 0.094339 | −0.035087 | 0.102951 | +0.003086 |
| **a8** (`--overdispersion 4.0`) | 08-27 | 0.1543607 | 0.736561 | 0.000000 | 0.094301 | −0.028694 | **0.124349** | −0.000353 |
| a10 | 08-28 | 0.1541841 | 0.736561 | 0.000000 | 0.094338 | −0.027897 | 0.123595 | −0.001493 |
| a11 (`VCC_ALPHA` 1.0796→0.838) | 08-28 | 0.1661150 | 0.739889 | 0.082041 | 0.096982 | −0.031535 | 0.112220 | −0.002907 |
| a12 (tiered spread) | 08-29 | 0.1679874 | 0.737885 | 0.091163 | 0.093974 | −0.031473 | 0.116790 | −0.000413 |
| a13 (donor-agreement shrinkage) | 08-29 | *−0.0203 vs a12* | | | *raw 1.0185* | | 0.097121 | |
| a14 (per-context cap) | 08-30 | 0.1704242 | 0.738614 | 0.110615 | 0.091972 | −0.032586 | 0.114542 | −0.000612 |
| **a15** (CD4 residual) | 08-30 | **0.1740865** | 0.745871 | 0.120057 | 0.094239 | −0.030397 | 0.115368 | −0.000619 |
| a19 (whitener ridge) | 08-31 | 0.1605690 | | *−0.01352 on `mse`* | | | | |
| a21 (`resample_x 0.14`) | 09-01 | | | | −0.064 | −0.029 | **−0.0185** | −0.017 |
| a23 (mass fix) | 09-02 | 0.1733565 | | | | | 0.107435 | |
| abc-knee (per-context α) | 09-03 | 0.1740437 | 0.747726 | 0.119322 | 0.096214 | −0.029477 | 0.110910 | −0.000433 |

## Raw members, where recorded

| arm | `pds` | `mse` | `nmae` | `fid` | `reach` | `jac` |
|---|---:|---:|---:|---:|---:|---:|
| a2-prune | 0.831773 | — | — | 0.5021 | 0.1701 | — |
| **a8** | 0.831773 | — | — | 0.5041 | **0.1890** | — |
| a11 | 0.833307 | 0.909989 | — | — | 0.178243 | — |
| a12 | — | 0.901349 | — | — | — | — |
| a14 | 0.832702 | 0.882775 | 0.945764 | 0.502881 | 0.180615 | 0.030081 |
| a15 | 0.835987 | 0.873805 | 0.944251 | 0.503536 | 0.181173 | 0.030080 |
| a23 | — | 0.872517 | — | — | — | — |
| abc-knee | 0.836786 | 0.874491 | 0.943171 | 0.503819 | 0.177303 | 0.030145 |

## What the ledger says that no single board result did

**1. `reach` peaked at a8 and has never recovered.** Raw `0.1890` on 2026-08-27; today
`0.177303`. Everything since — the α move (a11, −0.0114 scaled), the shrinkage arm (a13,
−0.0197), the mass fix (a23, −0.0079), the per-context α (abc-knee, −0.0045) — took a bite out of
it, and only a12's tiered spread (+0.0046) and a15 (+0.0008) gave any back. Cumulative
**−0.0134 scaled = −0.0022 overall since a8.**

Each trade was individually correct — the same period bought `mse` from raw 1.3524 to 0.8745,
worth **+0.1193 scaled = +0.0199 overall** — so this is not a list of mistakes. It is a record of
a member that has been paying for the others for six weeks and was never separately defended.

**2. Shrinkage does NOT explain our `reach` deficit, and the board proves it.** Among the 31
teams with raw `pds > 0.80` **and raw `mse < 1.0`** — i.e. shrunk at least as hard as we are —
the mean raw `reach` is **0.2106**, against our **0.1773**. Fourteen of them are listed in
`docs/CURRENT_STATE.md`; `Aginglab.com` holds raw `mse` 0.8270 with reach **0.2621**, and
`cqawesome` is *better* shrunk at 0.7604 with **0.2592**.

There is a genuine amplitude→`reach` coupling on the board — partial `corr(log raw mse, raw reach
| raw pds)` = **+0.543** over the 51 teams with `pds > 0.80`, and the unshrunk `mse > 3` bin
averages reach 0.2662 — but it does not reach down to us. **We are below the shrunk cohort's own
mean, so the deficit is ours, not the regime's.**

| raw `mse` bin | n | mean `reach` | mean `fid` | mean `nmae` | mean `pds` | mean overall |
|---|---:|---:|---:|---:|---:|---:|
| 0.0–1.0 | 31 | 0.2106 | 0.4984 | 0.9311 | 0.8205 | 0.17632 |
| 1.0–1.5 | 6 | 0.1984 | 0.5057 | 0.9455 | 0.8122 | 0.14898 |
| 1.5–3.0 | 5 | 0.1863 | 0.5055 | 0.9447 | 0.8160 | 0.14783 |
| > 3.0 | 9 | 0.2662 | 0.5106 | 0.9401 | 0.8142 | 0.16576 |
| **us** | | **0.1773** | 0.5038 | 0.9432 | **0.8368** | 0.17404 |

**3. The mechanism is on record and it is ours alone.** 171 §12c measured that `|z|` is U-shaped
in the emitted spread with its minimum at the shipped **0.25**, and that a15's tiered policy
(`--keep-top 8000 --keep-spread 0.25 --tail-spread 1.0`) therefore **parks the signal block at
the minimum-|z| point while the noise tail sits at the high-|z| end**: 61–79% of the head
positions of the reach walk are TAIL-block genes at sign accuracy 0.5212, against 0.6629 for the
KEEP block. a8, our `reach` high-water mark, used a **flat** 0.25 spread and no tier.

Matching the shrunk cohort's mean is **+0.0333 raw = +0.0063 overall**; matching Aginglab is
**+0.0848 raw = +0.0160 overall**.

## 3b. The 2026-09-03/04 entries — three arms in two UTC days, and what each one proved

| entry | arm | date | overall | Δ vs prev | the one variable |
|---|---|---|---:|---:|---|
| `nQUJZhZb3MW2uFJHBX5e` | `rc_hybrid` | 09-03 | 0.1834190095 | **+0.0093325** | the emitted p-ORDER, at bit-identical pseudobulk |
| `31oB9yWy8CSNBDEkY4iW` | `rc_callset` | 09-03 | **0.1918638367** | **+0.0084448** | the CALL SET (`top 400 by |H| ∪ top 2,200 by control CPM`, `n_pred` 4,480 → 2,599) |
| `QSe3JlNUYZ2BzvlcsL14` | `cd4f034` | 09-04 | **0.2049195830** | **+0.0130557** (pred +0.0044 / +0.0057; 210 `BOARD-RESULT.md`) | `VCC_CD4_TARGET_F` 0.10 → **0.34** — the first DIRECTION arm to clear the gate |

### Scaled members

| member | a15 | `rc_hybrid` | `rc_callset` |
|---|---:|---:|---:|
| `pds` | 0.745871 | 0.7458708290005253 | **0.7458708290005253** |
| `mse` | 0.120057 | 0.12005661519665856 | **0.12005661519665856** |
| `nmae` | 0.094239 | 0.09429017536041247 | 0.09429042916917090 |
| `fid` | −0.030397 | −0.029758 | −0.029268 |
| `reach` | 0.115368 | **0.170308** | **0.192898** |
| `jac` | −0.000619 | −0.000253 | **+0.027335** |

**`pds` and `mse` are identical across all three to SIXTEEN significant figures.** That is the
single most useful fact in this file. It was a construction argument on 09-02; after two submissions
it is a **repeated board measurement**, and it holds across an order change *and* a 79% call-set
change (jaccard 0.2072, `n_pred` −42%). It is what licenses pricing any emission-side arm on the DE
members alone.

### What the two ordering arms bought, in world rank

| member | 09-03 morning | after `rc_callset` |
|---|---:|---:|
| `jac` | 172 / 640 | **1 / 646** |
| `reach` | **192 / 640** | 29 / 646 |
| `pds` | 8 / 640 | 7 / 646 |
| `mse` | 26 / 640 | 26 / 646 |
| `nmae` | 80 / 640 | 83 / 646 |
| `fid` | 190 / 640 | 194 / 646 |
| **overall** | **24 / 640** | **10 / 646** |

And on the field's own conditional expectation, our `reach` residual went **−1.92σ → +0.70σ**
(215 §by-products). The ordering hole that 186 identified on 09-03 is closed.

### Two constants these points fixed, both of which had been wrong

> **Correction 2026-09-04 (221, on the real scorer):** the h1-test→board `reach` factor for the order-arm class is **0.769×** (h1-test OVERSTATES), not 1.290× as item 2 below says; the board half reproduces to four decimals, the whole disagreement is in the local half (`experiments/221-cloud-scorer/RESULTS.md`). Item 2 is kept as written for the record.

1. **Rule 13's `−0.0146` scaled `reach` toll does not apply to profile-PRESERVING emission arms.**
   `rc_hybrid` was charged it and realised **+0.054940**. The toll belongs to arms that move the
   delivered profile; both its calibration points (a21, a22) did. See the rule-13 table in
   `AGENTS.md` §11.1.
2. **The h1-test → board `reach` factor is not a constant.** For designed-order arms it is **1.290**
   — h1-test *understates* the board by 3× — against the **0.4277** that every prediction in this
   project had passed through (and the 0.390 the a8/a15 pair implied, §4 below). Both of those came
   from spread-policy arms. **Do not apply one factor across arm classes.**

---

## 4. A second, independent read on the h1-test → board `reach` shrink factor

164 §6 measured the factor at **0.4277** by one route. The a8/a15 pair gives another, though a
confounded one:

| | h1-test | board |
|---|---:|---:|
| a8-flat | 0.20291 | 0.1890 |
| a15 | 0.18283 | 0.181173 |
| Δ (a8 − a15) | **+0.02008** | **+0.007827** |
| implied shrink | | **0.390** |

**0.390 against 164's 0.4277 — agreement to 9%.** Treat it as corroboration, not as a measurement:
a8 and a15 differ in far more than the spread policy (α, the CD4 residual, the per-context cap), so
this pair is not the single-variable comparison the factor deserves. It is still the only second
opinion on a constant that every `reach` prediction in this project passes through, and both routes
say a h1-test `reach` delta arrives on the board at **roughly 40% of its size**.

**Consequence for sizing.** 172's ordering counterfactual — our own `|H|` key reading **0.23533**
against our delivered **0.18283**, same predictions and same signs — is **+0.05250 raw on h1-test**.
At a 0.40 shrink that is **+0.0210 raw on the board = +0.00395 overall**, and it lands our `reach`
at ~0.198, still below the shrunk cohort's 0.2106. So delivering `|H|` is necessary and **not
sufficient**; a better key (189, 190) has to compose with it.

## 5. A caution for anyone modelling `reach`: a homogeneous model is off by 4–5×

Coordinator, 2026-09-03. Fit a per-gene accuracy profile `a(k) = 0.537 + (a₀ − 0.537)·e^{−(k−1)/τ}`
to the *measured* cumulative purity curves — a15's `P(1)/P(3)/P(10)` = 0.7317 / 0.6491 / 0.5655 and
`rc_hybrid`'s 0.8780 / 0.7851 / 0.6914 — then simulate `k* = max{k : P(k) ≥ 0.9}` over `n_real` = 89
and average. It gives **reach 0.027 (a15) and 0.042 (rc_hybrid)** against measured **0.18283 and
0.21681**.

**So roughly three quarters of our `reach` does not come from the average head-accuracy profile at
all — it comes from between-perturbation heterogeneity**, i.e. from perturbations that are either
unusually accurate or have a small `n_real`. That is the same fact 172 reported from the other side
(*"94.4% of h1-test's reach lives in 32 of 85 targets"*) and the same one behind the harness trap in
`AGENTS.md` (*"98% of our reach comes from the small-`n_real` tercile"*), and `AGENTS.md` §2.1 now
explains why such a stratum exists at all: **the panel was deliberately binned by DE-gene count and
sampled across bins, so a low-`n_real` stratum is guaranteed by construction.**

**Two consequences.**
1. **Never model or price `reach` in aggregate.** Any statement of the form "head accuracy `p`
   implies reach `r`" is wrong by a factor of several. It must be computed per perturbation and
   averaged, over the real `n_real` distribution — and pools must be scaled **proportionally**,
   never flattened.
2. **The marginal value of average head accuracy is smaller than it looks, and the value sits in
   the stratum that already scores.** A `reach` arm should report its effect on the
   reach-carrying subpopulation separately from its effect on the mean.

Recorded as a caution rather than a finding: the simulation's purpose was to price a head-accuracy
target, it failed to reproduce the measurement, and the failure is the useful part.

---

## Arm: `h1union` (g=0.28) — entry `iUziqE1gewfcPKzYUDaB`, submitted 2026-09-04 03:52 UTC, published ~04:16Z: **0.2118310, rank 3, Δ +0.0069114 vs `cd4f034`** (222 `BOARD-RESULT.md`; full record `data/board-entries/iUziqE1gewfcPKzYUDaB.json`). THE INCUMBENT pending the 249 point.

- package `submissions/222-h1union-g028.vcc`, sha256
  `f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e`, 2,824,785,920 bytes
  (h5ad 12,648,807,800 bytes)
- model name `The Nucleus Declined To Comment`; **description empty** (§14)
- upload **md5-verified** by the server at 2,824,785,920 bytes; scoring job
  `vcc-metrics-prod-iuziqe1gewfcpkzyudab`
- second of the two 2026-09-04 UTC slots; the first was `cd4f034` (0.2049196), and `cd4f080`
  was cancelled unspent
- built by 222 behind `VCC_H1_UNION_G`; default path proven **byte-identical** to `cd4f034`
  (`7036e89f…`); CD4's solved fraction did not move by one bit, so the arm is **single-variable**

**Pre-registered reading (3 cancer folds, `kolf2` excluded as in 218):** relative MACRO cosine
**+7.6803%**, Δ`pds` **+0.004887**, relative aggregate cosine +10.3065%; 12 of 12 fold×dose cells
positive; target-shuffled control strictly negative at every dose (−0.5965/−1.4217/−3.1633%);
non-carrying-row excursion 0.01778 pp. Every arm reproduces 218 to **0.00e+00**.

**My board estimate, on 216's realised transfer factors** (k_cos 1.384, k_pds 1.987, measured on
`cd4f034`'s own board point): **≈ +0.0035 overall → ≈ 0.2084, rank 3.** This is **below** the
preregistration's +0.0055–0.0105 band, which was built on the `kolf2`-inflated 4-fold Δ`pds`
(+0.009166) rather than the 3-cancer +0.004887. **Read the band as reissued.**

**Two caveats 222 raised, both honest, neither a stop condition:**
1. Exact `nmae` is **slightly negative** (−0.000108 overall), inside falsifier 4's ±0.002 bound so
   it does not fire, but contradicting the preregistration's "up or flat" sign. The same instrument
   mis-signs its own family's board point, so treat it as a magnitude bound, not a sign.
2. **Δ`pds` is concentrated**: top-3 targets per cancer fold carry **99.8%** of the mean, and
   leave-top-3-out retains only **14.6%**. **The cosine gain is not concentrated — leave-top-3-out
   retains 72.7%.** The addendum's §5 pre-commitment fires: a positive board result is *weaker*
   evidence for 223's stack than it would otherwise be.

**Why the slot was spent anyway.** Cosine is the channel that matters (§4.1c) and its gain is
distributed; the arm is single-variable; and the board point directly calibrates **225's
source-stacking program**, which is running now. The concentration split makes the result *more*
informative, not less: if the board gain tracks the distributed cosine prediction rather than the
concentrated `pds` one, the mechanism is real transfer rather than three lucky targets.


---

## Housekeeping — reclaimed build artifacts (2026-09-04)

Deleted, all deterministically reproducible from committed, gated build scripts, with their
verification results recorded in the owning experiment's `RESULTS.md`:

| file | GB | why safe |
|---|---:|---|
| `submissions/217-ka-cd4f034-f034.h5ad` | 12.1 | **full sha256 identical** to `210-cd4f034.h5ad`, which is kept |
| `submissions/210-ka-rc_callset-f010.h5ad` | 12.1 | **full sha256 identical** to `202-rc_callset.h5ad`, which is kept |
| `submissions/hy-kt8000.h5ad` | 18.8 | **never submitted** (179: package built, 6/7 gates pass, NO slot spent); the arm priced ≈ **−0.0006** and is superseded by four shipped arms. `.meta.json` config record retained. |

**Reclaimed 2026-09-04 06:35 (77 GB), all build INTERMEDIATES whose submitted `.vcc` survives:**
`abc-knee.h5ad` (17.7 GB, `abc-knee.vcc` kept), `a23-massfix.h5ad` (17.7, `a23-massfix.vcc` kept),
`191-a15_local.h5ad` (17.7, `a15-cd4f010.vcc` kept), `201-rc_hist_a.h5ad` (11.7 — 201 carries
`SUPERSEDED.md` and was **never submitted**), `222-ka-cd4f034-f034.h5ad` (11.8 — 222 finished; a
**verified byte-duplicate** of `210-cd4f034.h5ad`, which is kept). Every `.meta.json` config record
was retained, and **all 34 historical `.vcc` artifacts are untouched.**

**Still deletable if space is needed**: `191-rc_hybrid.h5ad` (`191-rc_hybrid.vcc` kept), and 229's
two known-answer rebuilds **once 229 has finished** — its `229-ka-h1union-g028.vcc` was
independently verified by the coordinator to hash to `f13ce40b…`, byte-identical to the shipped
`222-h1union-g028.vcc` at 2,824,785,920 bytes, which is the proof that 229's registry change is
inert.

**Never delete**: the current incumbent's package and `.vcc` (`222-h1union-g028.*`), the previous
incumbent (`210-cd4f034.h5ad`), or any artifact an agent is still running against.

## Arm: `ARM6` (stack6, six union sources, g=0.35) — entry `ydAemVdWmsTx9xhDWcou`, submitted 2026-09-05 00:02:21Z, published 00:28:02Z: **0.2081906728, rank 4, Δ −0.0036403 vs `h1union`. A GATE TEST, not an arm; the incumbent is unchanged.**

- package `submissions/229-stack6-g035.vcc`, sha256
  `7a3c2d3e414f8d9dfba0bd1c42f2b93994dedfda73cf207cea9b7aabdc78291c`; model name `I'll just keep moving forward.`;
  **description empty** (§14); upload 2.6 GiB in 12m43s; job `vcc-metrics-prod-ydaemvdwmstx9xhdwcou`
- **full record** `data/board-entries/ydAemVdWmsTx9xhDWcou.json` (the API stops resolving it once the 249 entry
  publishes); reading `experiments/237-gate-validation/BOARD-RESULT.md`, transcript `reading-ydAemVdWmsTx9xhDWcou.txt`
- first of the two 2026-09-05 UTC slots; the second is `cd4svd5` (249), fired 00:30:21Z after this record was captured

| member | raw | scaled | Δ scaled vs `h1union` | overall contribution |
|---|---:|---:|---:|---:|
| `pds` | 0.8423782980 | 0.7602722152 | +0.01490 | **+0.00248** |
| `mse` | 0.8122683282 | 0.1848194290 | +0.00525 | **+0.00088** |
| `nmae` | 0.9601566451 | 0.0677845938 | **−0.04296** | **−0.00716** |
| `fid` | 0.5075662030 | −0.0161824516 | +0.00029 | +0.00005 |
| `reach` | 0.2769826076 | 0.2235221007 | +0.00058 | +0.00010 |
| `jac` | 0.0409876541 | 0.0289281499 | +0.00009 | +0.00002 |

**What it proved.** The 216/188 exact `nmae` identity — mis-signed on both earlier, small-magnitude checks — is
**validated at scale**: predicted −0.004783 overall, realised **−0.00716**, same sign, conservative by a third. 243's
class rule is confirmed on a **fourth** point: raw `mse` 0.81227 against the point prediction 0.81244, realised
**k_cos 0.962** in the new-source bin (0.94–1.00), not the coverage bin. Direction legs excluding `nmae` +0.00352,
inside 243's +0.0022 … +0.0038. The six-source union's `pds`+`mse` gain (+0.00336, raw `pds` +0.00666 = 1.9× the
`h1union` step) is real and its `nmae` cost is the entire loss — the open hypothesis is an `nmae`-neutral rendering
of the same directions, to be priced on the now-validated identity (237 `BOARD-RESULT.md` §3).

## Arm: `cd4svd5` (`VCC_CD4_DEFLATE_K=5`) — entry `Xv0ayvuoHpsKYWv2iu42`, submitted 2026-09-05 00:30:23Z, published 00:54:57Z: **0.2160944022, rank 3 of 688, Δ +0.0042634 vs `h1union`. THE INCUMBENT.**

- package `submissions/249-cd4svd5-f034.vcc`, sha256
  `6e8b41c6244b32dbe2a9765c7d2f2293b2032c49e59b455982b95bc6bb074533`; model name `You can do it too.`;
  **description empty** (§14); second of the two 2026-09-05 UTC slots, after `ARM6`
- **full record** `data/board-entries/Xv0ayvuoHpsKYWv2iu42.json`; reading `experiments/249-cd4-deflate/BOARD-RESULT.md`
- one variable: the CD4 residual's top-5 target-shared directions deflated before sizing (249; amplitude −0.409 %)

| member | raw | scaled | Δ scaled vs `h1union` | overall contribution | P249 predicted |
|---|---:|---:|---:|---:|---:|
| `pds` | 0.8506428837 | 0.7780868669 | **+0.03271** | **+0.00545** | +0.0032 / +0.0046 / +0.0091 |
| `mse` | 0.8218061907 | 0.1747722967 | −0.00480 | **−0.00080** | −0.0003 |
| `nmae` | 0.9339736196 | 0.1104721876 | −0.00027 | −0.00005 | −0.0001 … +0.0003 |
| `fid` | 0.5073946332 | −0.0168327489 | −0.00036 | −0.00006 | −0.0001 |
| `reach` | 0.2750615691 | 0.2213925724 | −0.00155 | −0.00026 | −0.0002 |
| `jac` | 0.0408976017 | 0.0286752386 | −0.00016 | −0.00003 | ≈ 0 |

**What it proved.** A `pds`-first re-shape of the rows CD4 occupies transfers to the board at **k_pds 1.22** (raw
+0.01492 vs PLOLO +0.01226) — five of six members inside their priced bands. Its cosine cost transfers at **k_cos 3.76
on `cos_agg`** (board −1.735 % vs PLOLO −0.4617 %): a third class beside re-dose (1.36) and new-source (0.94), and
the reason the arm netted +0.0043 instead of its +0.0086 upside. Prereg falsifier 2 fired (raw `mse` +0.00454 >
+0.0030); falsifiers 1, 3 did not; 4, 5 pass off the package. October: the pinned `V5_abc` basis SHIPS (RB §3.1
trigger fired); k is never raised on it.

## Arm: `319d3` (`VCC_NORMREL_EXP=-3.0 VCC_NORMREL_MODE=cap`) — entry `5LOjb2VLiIBh0JG6aAgg`, submitted 2026-09-09 22:36:26Z, published by 23:02Z: **0.2183911409, rank 13, Δ −0.000269 vs 319 (cap −1.5). BIN 2, FLAT — the cap dose family is CLOSED at −1.5.**

- package `submissions/319d3-cap30.vcc`, sha256 `6ddbf76d2c7b914fae09f606f5da44262f25e415870552f211d606fa08c33d74`; model name `a moth in the lint trap`; description empty; first of the two 2026-09-09 UTC slots (the second expired unused by the frozen reading)
- **full record** `data/board-entries/5LOjb2VLiIBh0JG6aAgg.json`; reading `experiments/319-cap-variant/BOARD-RESULT-d3.md`
- one variable: `normrel_exp` −1.5 → −3.0 on the cap form

| member | scaled | Δ scaled vs 319 |
|---|---:|---:|
| `pds` | 0.793795 | **+0.005481** |
| `mse` | 0.181321 | +0.000308 |
| `nmae` | 0.101252 | **−0.006263** |
| `fid` | −0.015686 | +0.000073 |
| `reach` | 0.221041 | −0.001151 |
| `jac` | 0.028624 | −0.000060 |

**What it proved.** The 2×2 plain/cap × −1.5/−3.0 is complete on the board: plain −0.00127, cap −0.00027, interaction +0.0010. The cap form neutralises plain's MSE/reach loss with dose but keeps the PDS-for-NMAE trade (pds +0.0055, nmae −0.0063). No further exponent earns a slot on either form. 326's local raw-MSE harm (+0.015 on H1) did not transfer: board `mse` flat.

## Arm: `254 cd4svd5+partial` (`VCC_CD4_DEFLATE_K=5 VCC_CD4_PARTIAL_G=0.34`) — entry `PVHrJHweHD5gvOf43ftL`, submitted 2026-09-09 23:10:09Z, published 23:42Z: **0.2174269300, rank 13, Δ +0.0013325 vs `cd4svd5` (its own base). BIN 1 — the partial rows TRANSFER, via `mse` (raw −0.0045).**

- package `submissions/254-cd4svd5-partial-f034.vcc`, sha256 `eea5ed1802a2337dc3d69b38bec4e7962835ba15b194c650015049b7ff035bf0`; model name `the intern watered the plastic plant`; description empty; second 2026-09-09 UTC slot, re-decided at 23:10Z from the unsubmitted-package inventory
- **full record** `data/board-entries/PVHrJHweHD5gvOf43ftL.json`; reading `experiments/254-partial-rows/BOARD-RESULT.md`
- one variable on `cd4svd5`: `cd4_partial_g` absent → 0.34. Contributions: mse +0.00079, reach +0.00027, nmae +0.00015, fid +0.00014, pds −0.00001, jac 0.
- **What it proved.** CD4's 15 partial-condition rows carry transferable expression information. The drafted PDS route did not fire; MSE did. Next: `VCC_CD4_PARTIAL_G=0.34` as a one-variable rider on 319 (cap −1.5), priced on 323 first.

## Arm: `319p` (319 + `VCC_CD4_PARTIAL_G=0.34`) — entry `SvzrTOQJAAaNjawzSyrC`, submitted 2026-09-10 00:16:37Z, published 00:46Z: **0.2198189778, Δ +0.0011593 vs 319. BIN 1 — composes. NEW BEST.**

- package `submissions/319p-cap15-partial.vcc`, sha256 `e8bd9cb5796f2bcd4971b413020338a7a8eb55a64ca08824562d5f5dff29ab4e`; model name `the stapler has tenure`; description empty; first 2026-09-10 UTC slot
- **full record** `data/board-entries/SvzrTOQJAAaNjawzSyrC.json`; reading `experiments/319-cap-variant/BOARD-RESULT-p.md`
- one variable on 319: `cd4_partial_g` absent → 0.34. mse +0.00425 scaled (+0.00071), reach +0.0017, nmae +0.0006, fid +0.0006, pds −0.0002, jac 0.
- **What it proved.** The partial-condition CD4 rows are additive across bases (+0.00133 on cd4svd5, +0.00116 on cap). The composition is now cd4svd5 + cap normrel −1.5 + partial rows @0.34.

## Arm: `319pA` (319p + `VCC_CTX_A_PRIOR` blend_rows β=0.35, Jurkat-calibrated rows for 259 A targets) — entry `sEa0LAsWFnXOoXSPJZ59`, submitted 2026-09-10 02:14Z, published 02:40Z: **0.2185396201, Δ −0.0012794 vs 319p. BIN 3 — pds pays; closed at this construction.**

- package `submissions/319pA-cap15-partial-ctxA035.vcc`, sha256 `fe325f45abfea30666fc63f81ad1739b695b6cc1e6a85e83c8b267aa972d48d9`; model name `the fern unionized`; second 2026-09-10 UTC slot; builder patch c47d2343 (inert when unset)
- **full record** `data/board-entries/sEa0LAsWFnXOoXSPJZ59.json`; reading `experiments/375-ctx-a-arm/BOARD-RESULT.md`
- one variable on 319p: context A rows blended toward the Nadig-Jurkat-calibrated prediction. pds −0.0090 scaled (−0.0015 overall), mse +0.0026, fid +0.0011, nmae 0, reach −0.0023, jac 0.
- **What it proved.** Same-line calibration at 26.7% axis coverage costs pds on the real context at the rate the H1 harness priced and buys almost nothing on the DE members there. The source-frame retrieval proxy (374) mispredicted the pds cost; H1 mispredicted the DE gains. 319p stays best.

## Arm: `387 sign repair` (`VCC_SIGN_REPAIR=cons` on 319p2) — entry `LRbQ4eOEGLjW7yK82mVa`, submitted 2026-09-13 07:06:50Z, published 07:32Z: **0.2186049, rank 19, Δ −0.00151 vs 319p2. Lower bin — pds pays.**
- package `submissions/387-cap15-partial50-signcons.vcc`, sha256 `986f7665488a54e3e542e71f3e779c26901fd5002e69ec0f7e49454904f7061c`; model name `the lint roller has opinions`; first 2026-09-13 UTC slot
- **full record** `data/board-entries/LRbQ4eOEGLjW7yK82mVa.json`; reading `experiments/387-sign-repair-arm/BOARD-RESULT.md`
- one variable on 319p2. Contributions: pds −0.00375, mse +0.00086, reach +0.00067, fid +0.00052, nmae +0.00019, jac 0.
- **What it proved.** Consensus sign repair through the bulk buys DE members inside noise (+0.0014) for a resolved pds cost 2.3× H1's; sign information must be delivered outside the bulk (392).

## Arm: `319p2` (319p with `VCC_CD4_PARTIAL_G=0.5`) — entry `lqWNXX5YaYArdglu6Stn`, submitted 2026-09-11 19:10Z, published 19:41Z: **0.2201159, rank 15, Δ +0.000297 vs 319p. Bin 2 (flat) — partial-rows dose family closed; best-scored package until 392.**
- package `submissions/319p2-cap15-partial50.vcc`, sha256 `a3f770aa9f51d56ca44da39e54bf2f4079ec1586a8ab33bcaf389a71e7a4d58e`; model name per `experiments/319-cap-variant/BOARD-RESULT-p2.md`; first 2026-09-11 UTC slot.
- **full record** `data/board-entries/lqWNXX5YaYArdglu6Stn.json`; reading `experiments/319-cap-variant/BOARD-RESULT-p2.md`.
- one variable on 319p: `cd4_partial_g` 0.34 → 0.5. Raw: pds 0.854738 / mse 0.809898 / nmae 0.934784 / fid 0.507933 / reach 0.277798 / jac 0.040903.

## Arm: `392` (`VCC_DE_SIGN_PRIOR=clean4` on 319p2: CLEAN4 per-gene sign tilt on the called CPM block through the decoupled DE channel) — entry `TJTOj9qrbiobJBpqLSa3`, submitted 2026-09-13 07:46:36Z, published 08:10Z: **0.2226221, rank 17 at publication, Δ +0.0025062 vs 319p2. fid TOP bin; overall bin "partial transfer, incumbent". NEW BEST.**
- package `submissions/392-cap15-partial50-designclean4.vcc`, sha256 `20b987a4173157f11421b7dd25f1cb115cd8c97c03d9a4989a6b12dfda5de4ab`; model name `the stapler went on sabbatical`; second 2026-09-13 UTC slot; builder 260e03aa…, emitter `src/vcc2026/rc_emit_k.py` a12cf8f9…, prior 7e532f84….
- **full record** `data/board-entries/TJTOj9qrbiobJBpqLSa3.json`; reading `experiments/392-de-sign-prior-arm/BOARD-RESULT.md`; registration `PREREGISTRATION.md` (commit 175d8454).
- one variable on 319p2 (five `de_sign_*` config keys). Raw: pds 0.854738 (exact) / mse 0.809898 (exact) / nmae 0.938837 (+0.0041 worse) / fid 0.516157 (+0.0082) / reach 0.273878 (−0.0039) / jac 0.040898. Contributions (overall points): fid +0.00437, nmae −0.00116, reach −0.00070, jac 0, pds 0, mse 0.
- **What it proved.** A target-agnostic cross-line sign tilt carries A/B/C information for the highly expressed block, and the unwelded per-cell library channel delivers it with the bulk untouched (pds/mse identical to the printed digit). The members that read the strongest block genes (nmae, reach) did not follow the H1 forecast → 395 (q80 protection). Read four days late: nobody was at the machine 09-13 08:11Z → 09-16 05:00Z; four slots expired.

## Arm: `395` (`VCC_DE_SIGN_PRIOR=clean4q80` on 392: the strongest 20 % of block calls keep their delivered sign) — entry `RzcGA26xZEigSuXgUp5t`, submitted 2026-09-16 05:50:22Z, published 06:21Z: **0.2239363, rank 16 at publication, Δ +0.0013142 vs 392. Bin "+0.0003 … +0.0015", incumbent (Δfid raw −0.0019 ≥ −0.002). NEW BEST.**
- package `submissions/395-cap15-partial50-designclean4q80.vcc`, sha256 `19435f099a82026c40c107db7cd309e97e53342d5fba013b46038463afa18ddb`; model name `the thermostat has a union rep`; first 2026-09-16 UTC slot; builder 05e49916…, emitter 6ff0b866…, prior 7e532f84….
- **full record** `data/board-entries/RzcGA26xZEigSuXgUp5t.json`; reading `experiments/395-protected-signs-arm/BOARD-RESULT.md`; registration `PREREGISTRATION.md` (commit d86fd5b3).
- one variable on 392 (`de_sign_protect_q` 0.8; the knob value and the emitter re-pin are the other two config keys). Raw: pds 0.854738 (exact) / mse 0.809898 (exact) / nmae 0.937101 (−0.0017, better) / fid 0.514228 (−0.0019) / reach 0.283909 (+0.0100) / jac 0.040891. Contributions: reach +0.00188, nmae +0.00047, fid −0.00102, jac 0, pds 0, mse 0.
- **What it proved.** Protecting the strongest block calls from the marginal tilt's flips buys reach on A/B/C (the front of the order) at a third of the fid gain; the q50 dose is authorised.

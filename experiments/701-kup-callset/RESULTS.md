# 701 — K-UP: the call-set size K 3,400 (ku30) on the board, and the H1 size curve behind the October K rule

- **Board (10-16 s2), recorded in docs/BOARD-LEDGER.md:** ku30 = 446bj2z + CPM 3,000 (K 3,400), entry HAbelEdgYYziqwfjTgdI,
  0.2740045, against 446bj2z's 0.2728260 (+0.0011785). The preregistered primary, Δ(fid + jac), read NULL (+0.000123 overall).
  reach moved +0.005398 raw. Delivery was full (3,339–3,373 of 3,400 calls per target; deliver701.py).
- **Follow-ups:** 711 (K 4,200 on the board, 10-18 s2) and the October K rule, written in 711's PREREGISTRATION with the
  coordinator's Amendment 1 (step-down). The rule needs condition (i), the H1 size curve below (run701b / run701c +
  reach701.py).

## The October K rule's condition (i): FAILS (H1 size curve, POST HOC as preregistered; recorded 2026-10-16T12:34Z, [d80a10])
- **reach701.py (the rule's statistic, f 0.87):** the per-seed slopes of reach on CPM are −0.01169, −0.01106 and −0.00851 per
  1,000 CPM. The mean is −0.01042 (se 0.00097). Mean Δreach against CPM 2200 (K 2,600): 2800 −0.01741, 3000 −0.01767, 3400
  −0.02166, 3800 −0.01635. **Condition (i) FAILS** (f 1.0 also fails: mean slope −0.00715, se 0.00033). By 711's rule
  (Amendment 1 included), **October K stays 2,600 (rc_callset_cpm off) whatever 711's board read shows.** The rule needs
  both the board and the H1 mechanism.
- **EXPLORATORY, run701d (K 5,000 = CPM 4600):** the step 3800 → 4600 is −0.00440 (se 0.00190) at f 0.87 and −0.00416 at
  f 1.0. The coordinator's fallback (1), a K 5,000 slot, is dead.
- **EXPLORATORY, regimes outside the rule (after seeing the preregistered result, so it cannot decide anything):** mean Δreach
  against CPM 2200 at f 0.61 is 2800 +0.00041, 3000 +0.00462, 3400 −0.00209, 3800 −0.00163, 4600 −0.00436. At f 0.35 it is
  +0.00383, +0.00607, +0.00632, −0.00112, −0.00099. So at the lower regimes H1 REPRODUCES the board's K 3,400 gain (ku30:
  board reach +0.0054 raw; H1 f 0.61 +0.0046) and predicts no further gain at K 4,200 (CPM 3800). At f 0.87 / 1.0 every
  larger size loses.
- **The conflict to resolve before anyone leans on it:** 711's preregistration took f 0.87 from FINDINGS 0b ("the board
  sits near H1's f 1 regime for SIZE changes"). 662 puts the board's call-set level at f ≈ 0.61, and the lab's calibration
  has DE-H1 tracking call-set changes at f 0.61 (ρ 0.86, no sign failures) but not at f 0.35. If 0.61 is the right dial
  for size, the board / H1 disagreement disappears, and K 3,400 (not 4,200) is where the reach gain sits. That claim needs
  its own preregistered test: 711's board read is one, H1 at f 0.61 predicting K 4,200 below K 3,400. Choosing the regime
  after seeing this curve is not allowed (AGENTS 11.1c.4).

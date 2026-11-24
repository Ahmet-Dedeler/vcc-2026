# 477 — PREREGISTRATION (third session, 2026-09-26 23:35Z, at the coordinator's request; before any number)

480's donor-consensus rider (A_pg → A_pg · 2^(k · cons_pg), cons = the panel donors' same-target consensus log2 ratio,
pc 1 CPM), priced ANALYTICALLY through 490's exact emitter reconstruction (`recon_h1.py`) on 435's three coarse-MDE H1
seeds. The ballast is carried as delivered, and calls / p / reach order are held. `price477.py`. This is the cheap first
step. A real price needs an emitter hook and an H1 run (class capacity, call leakage).

- Arms: mixed k 0.1 / 0.2 / 0.3, magnitude-only k 0.1 / 0.2. k 0 is the known answer (Δ exactly 0).
- Prediction (from 480's direct-lfc price, −0.0030 raw nmae at k 0.2 on h1-test): mixed k 0.2 nmae −0.002 … −0.003 raw.
  fid / reach move slightly through the sign flips the mixed form makes. jac is 0 (calls held). Overall-equiv
  (475's slopes) +0.0006 … +0.0010.
- GO for an emitter hook + H1 run (not a slot):
  - mixed k 0.2 overall-equiv > 0 on 3/3 seeds, with the pooled bootstrap lower bound > 0 on at least 2;
  - nmae ≤ −0.0015 raw;
  - the non-essential subset's overall-equiv > 0.
  Otherwise the rider is closed at this construction.

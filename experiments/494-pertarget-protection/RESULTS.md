# 494 — Per-target q80 protection by design-strength tercile: null on H1 (≤ +0.0003 overall-equiv). Closed.

2026-09-27 00:35Z, third session. Exploratory, analytic through 490's exact reconstruction.
- 435's three coarse-MDE H1 seeds, calls / p / order held, cell-eval2's own DE functions.
- Script `price494.py`, output `work/price494.{json,log}`, 34 s.
- The protection quantile is set per target by the design row's strength tercile (weak / mid / strong = 29 / 28 / 28
  targets). kf and the MDE balance are re-derived per target.
- Known answers: the reconstruction's diag known answers, the harness tables, and the base arm == the delivered tables
  exactly (Δ 0 on every member).

| arm (protection weak / mid / strong) | fid | reach | nmae | overall-equiv (3 seeds) |
|---|---|---|---|---|
| W (1.0 / 0.8 / 0.8) | +0.0004 … +0.0005 | −0.0006 … −0.0013 | −0.0004 … −0.0005 | +0.0003 / +0.0002 / +0.0001 (CI ±0.004) |
| W7 (1.0 / 0.8 / 0.7) | +0.0005 … +0.0006 | −0.0006 … −0.0013 | −0.0003 … −0.0004 | +0.0004 / +0.0002 / +0.0001 |
| S7 (0.8 / 0.8 / 0.7) | +0.0001 | 0 | +0.0001 … +0.0002 | +0.00002 |
| R (control: 0.7 / 0.8 / 1.0) | −0.0005 | −0.0007 … −0.0008 | +0.0004 … +0.0005 | −0.0005 (3/3) |

- Removing protection on the weak third trades a little fid and nmae for reach, all within ±0.0013 raw.
- The reverse loses, so the direction is right, but the size is inside H1's ordering noise, and H1 transfers ordering at
  ~0.4×.
- **Closed:** q80 is near-optimal per target as well as globally (395's finding extends).

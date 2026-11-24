# 400 — Uncalled-gene steering from neighbourhood confidence is null end to end on H1; the shipped 399 construction reads fid +0.0185 / reach +0.0174 / nmae −0.0006

Status: **COMPLETE** (τ ∈ {0.2, 0.1, 0.05}) (PI-run, 2026-09-16
05:41–06:05Z; `run400.py` = 394's harness with `k_emit400.py` = 394's emitter copy + an `extra_mask` OR'ed into the steered set;
385's exact call chain; cap15 known answer reproduced exactly first, `work/known-answer-cap15.json`).

## Arms (K = 3 library classes, q80 protection, per-target prior = 399's construction on the H1 targets: six-bank votes of the 25
nearest GWPS perturbations, shrink 400, family-marginal default; 82/85 targets with a row)
| arm | pds | mse | nmae (lower better) | fid | reach | jac |
|---|---|---|---|---|---|---|
| cap15 (production) | 0.88964 | 0.97495 | 0.95084 | 0.49755 | 0.25435 | 0.07476 |
| **k3knn (= 399)** | same | same | **0.95023** | **0.51602** | **0.27176** | 0.07472 |
| k3knn_u20 (+ uncalled genes at conf ≥ 0.2) | same | same | 0.95064 | 0.51594 | 0.27023 | 0.07485 |
| k3knn_u10 (conf ≥ 0.1) | same | same | 0.95017 | 0.51647 | 0.26713 | 0.07502 |
| k3knn_u05 (conf ≥ 0.05) | same | same | 0.94914 | 0.51901 | 0.26292 | 0.07513 |
pds/mse identical to cap15 in every arm (bulk untouched). k3knn vs cap15: fid +0.0185, reach +0.0174, nmae −0.0006, jac −0.00004
— the strongest DE-side H1 row of any construction (391's clean4: fid +0.0122, reach +0.0175, nmae +0.0024). Dose: λ median 1.0,
2 % of targets below 1; shrink k median 0.57; 880 flips per target (of 2,200 block calls).

## Uncalled steering (vs k3knn)
u20: nmae +0.0004, fid −0.0001, reach −0.0015, jac +0.0001 (443 extra genes steered per target at full dose); u10: nmae −0.0001,
fid +0.0005, reach −0.0046, jac +0.0003; u05: nmae −0.0011, fid +0.0030, reach −0.0088, jac +0.0004 — the more uncalled genes are
steered, the more reach falls while fid/nmae creep up; the board's reach slope (+0.188/raw) makes the trade negative. The analytic forecast (397c: nmae −0.0015 resolved at τ = 0.2) does not survive delivery:
the extra genes enter the sign-preserving mass balance and the class split, and the walk's order shifts. **Closed at this
construction.**

## Limits
H1 development panel; one line; reach noise ±0.03–0.04. The 399 board reading is the A/B/C instrument.

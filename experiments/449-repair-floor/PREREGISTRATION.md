# 449 — the zero-repair floor −5 → −2 (446br2). Preregistration, 2026-09-25 01:35Z
Written before any 09-26/27 entry exists and before 446br2 is built.

## Why
433h's repair (−35 → −5 on 24 zero-total gated cells) gained nmae −0.0247 on the board, 3× the preregistered bound that
assumed every zero sat in a typical gate. That only adds up if the zeros sat in WEAK targets' small gates (Σ|l| ≈ 5–20),
where one cell's error carries several units of that target's nmae. The 446b audit (`work/audit-446b-extremes.txt`, from
the restored package) has 37 gated cells below −4.5 or above +3:
- the repaired cells sit at −5.0 ± 0.3;
- most are TRANS genes that recur across unrelated targets (MTFP1 under ZNF410 and ZBTB48, SOD1 under ZDHHC4 and
  ANKRD52, H3F3B under RAB5A and UNK, UBXN10 under MAPK12 and CDK8, CNFN under AKT2 in A and B);
- most sit in low-key targets (w 0.3–1.5).
They are K562-only knockdowns: strong non-cis K562 knockdowns replicate in H1 about 12–20 % of the time, and when they
do, the H1 truth median is about −0.5. A gated cell predicted at −5 against a truth of −1 … −2 still costs 3–4 units. At a
floor of −2 it costs 0–1. Only the 2 cis cells (ANKZF1→ATG9A, PPIP5K2→GIN1) have H1-like truths near −4.
`VCC_DE_ZERO_REPAIR` already takes the floor as its value, so no code changes.

## The arm
446br2 = 446b with `VCC_DE_ZERO_REPAIR=-2` (config diff exactly `de_zero_repair`).
- The repair pays each added count from the same cell's largest sub-gate gene, after the ballast.
- The ballast rows, every library and every other tested gene are unchanged (asserted by the module).
- pds and mse move only by the ~24 blocks' tiny per-gene total shift.

## Bins (446br2 − 446b, raw)
- **nmae:** Δ ≤ −0.002 means the shallower floor pays; the trans zeros' truth is moderate. **Prediction: −0.002 … −0.008.**
  |Δ| < 0.001 means flat: keep −5. Δ > +0.001 means the gated truths are deep: keep −5, and consider a cis-only deep
  floor.
- **fid, jac, reach:** |Δ| ≤ 0.001 / 0.001 / 0.01. A repaired cell at 25 % of control (−2) is still called down. Anything
  outside means calls moved: diagnose.
- **pds, mse:** |Δ| ≤ 0.0005 raw.
- **Overall:** prediction +0.0006 … +0.0026. Ship −2 in the October stack if nmae lands in its first bin.

## 446br2 verified 03:03Z (still before any 09-26/27 entry)
Package `submissions/446br2-…-mde010coarse-zr2.vcc`, sha256 `9bbc8d2a…`, model name "the pipette rack learned origami".
- Config diff vs 446b is exactly `de_zero_repair`.
- 0 zero-total gated genes; novelty rc 0.
- Render energy vs 410b 0.9601 / 1.0280 / 1.0783 (446b: 0.9603 / 1.0282 / 1.0785). The ~0.02 % bulk move is the
  repair's per-gene total shift, as the pds/mse bin allows.
- Extremes audit (`work/audit-446br2-extremes.txt`): **13** gated cells beyond −4.5 / +3, against 446b's 37. The
  repaired cells now sit near −2. The remaining 13 are delivered knockdowns, not repairs.

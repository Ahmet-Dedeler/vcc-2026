# 391 — CPM-block signs from the 4-bank cross-line prior, end to end on H1: fid +0.0122 resolved, pds/mse exactly pinned, the deranged-prior control loses

Status: **COMPLETE, POSITIVE on fid** (PI-run, 2026-09-13 07:03–07:14Z). Board arm: 392.

## Known answer
`run391.py cap15` (389's harness, 385's chain) reproduced 323's `metrics-cap15.json`: emitted-count digest, panel-total
digest, mse legs and all six raw members exact (the driver aborts on any mismatch).

## Arms (K = 3 library classes at h1 control quartiles; emitter = snapshot of 389's `k_emit389.py` sha 3a375dfa… plus a
"prior" target mode; call set CallSet(400, 2200, 1.9) unchanged)
- `k3c5b`: CLEAN4 prior sign (K562-GWPS, K562-essential, Nadig Jurkat, Nadig HepG2; uncovered → down) on the CPM-block
  genes at the delivered magnitude; down-signed magnitudes × one per-target factor k so the block's arithmetic mass is
  unchanged (390 C5b); dose λ ∈ [0, 1] bisected for class capacity like 389's `feasible_dose`.
- `k3c5bshuf`: the same prior deranged across gated genes (seed 20260913), same balance and dose logic.
- `k3down`: every block gene down, same balance and dose logic (390 R2).

## Results (raw Δ vs cap15, paired 10,000-target bootstrap, report391.py → work/results391.json)
| arm | pds | mse | fid | reach | nmae (lower better) | jac | H1-local overall |
|---|---|---|---|---|---|---|---|
| **k3c5b** | **0** | **0** | **+0.0122 [+0.0067, +0.0182]** | +0.0175 [−0.019, +0.055] | +0.0024 [−0.0011, +0.0066] | −0.00001 | +0.0052 [−0.0029, +0.0132] |
| k3c5bshuf | 0 | 0 | −0.0019 [−0.0063, +0.0024] | −0.0267 [−0.072, +0.012] | **+0.0046 [+0.0022, +0.0071] harm** | +0.00001 | −0.0069 [−0.016, +0.001] |
| k3down | 0 | 0 | +0.0006 [−0.0018, +0.0033] | +0.0009 [−0.019, +0.018] | **+0.0026 [+0.0010, +0.0042] harm** | +0.00001 | −0.0005 [−0.0043, +0.0030] |
End to end reproduces 390's analytic C5b (fid +0.0116, reach +0.018, nmae +0.0024) almost exactly.

## Diagnostics (per target, q10 / median / q90)
| arm | dose λ | targets at λ = 1 | targets λ ≤ 0.5 (no flips possible) | shrink k | block genes flipped |
|---|---|---|---|---|---|
| k3c5b | 0.95 / 1.00 / 1.00 | 87 % | 1 % | 0.26 / 0.32 / 0.37 | 1,076 / 1,113 / 1,146 |
| k3c5bshuf | 0.94 / 1.00 / 1.00 | 87 % | 1 % | 0.30 / 0.35 / 0.39 | 1,076 / 1,105 / 1,135 |
| k3down | 0.30 / 0.52 / 1.00 | 22 % | **45 %** | 0 | 1,053 / 1,093 / 1,121 |
**All-down fails on class capacity, not on the idea:** sending every block gene down asks all block mass to move toward the
high-library class, so the feasible dose stalls near the flip threshold (λ = 0.5) in half the targets. The prior moves
genes both ways and runs at full dose; the deranged control runs at the same dose, so it is a fair targeting control.

## A/B/C audit of the 319p2 package (20 sampled targets per context, `work/abc_sign_audit.json`)
Block delivered down fraction A 0.500 / B 0.496 / C 0.526; the prior would flip 53.1 / 51.7 / 52.7 %; prior coverage of
the block 87 / 79 / 77 %; block mass-weighted delivered LFC ≈ 0 (+0.003 / +0.004 / −0.002). Same shape as H1.

## Limits
H1 development panel (85 targets, one line; the line where 371 found the prior weakest). Reach is unresolved. The prior
hurts a HEK293T-like line in 371. Board transfer is what 392 reads.

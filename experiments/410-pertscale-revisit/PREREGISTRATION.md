# 410 — per-target bulk magnitude from the whitened shared-energy key (welded arm), preregistered 2026-09-22 before any build

## One variable
`VCC_PERT_SCALE=wshare_lo0.3_hi3` on **395's configuration** (entry `RzcGA26xZEigSuXgUp5t`, 0.2239363, the incumbent).
Every final effect row p is multiplied by `w_p * c`:

- `key_p = mean_{i<j} sym<W d_i, d_j> / |W mean_i d_i|^2` over the donors carrying p among k562_gwps / orion_hct116 /
  orion_hek293t (raw bank rows, W = union_pruned e2f65ecf, gamma 0.65, no preserve constant) — the cross-fitted plug-in
  estimate of the mse-optimal per-target scale `<P_p,T_p>/|P_p|^2`, a held donor standing in for the truth;
- `w_p = clip(key_p / median, 0.3, 3)`; 6 targets carried by < 2 donors get 1.0;
- `c` preserves the whole-matrix effect norm over non-panel columns (the arm moves the per-target DISTRIBUTION of
  magnitude only, not the global amplitude).

Artifact `work/pert_scale_410_lo0.3_hi3.npz`, sha256 `d7f76a22f62ad6991d7220ce8a65218db38b2d71dcfe3c3781ec6a8f70409e0e`
(`build_weights410.py 0.3 3`). On the panel: 294/300 carried, 15.6 % of keys <= 0, 27.9 % of targets at the floor,
14.3 % at the cap; top keys TAF4, SYVN1, SMARCA5, PQBP1, UPF3B, ACLY, SIN3B.

**Expected config diff vs 395's meta, exactly three keys:** `pert_scale` (absent -> `wshare_lo0.3_hi3`),
`pert_scale_sha256` (absent -> `d7f76a22…`), `de_sign_emitter_sha256` (`6ff0b866…` -> `dbbb9121…`: the 09-16 guard for
the no-feasible-dose path, which 395's build never reached — output-identical for 395's configuration).

## Evidence before the build (PLOLO, deep pseudobulk truth, fold-held operators; `check410.py`)
| held fold | rho(key, s*) | aggregate cosine gain at lo 0.3 hi 3 | 95 % CI (target bootstrap) |
|---|---:|---:|---|
| k562_gwps | +0.576 | +0.299 | [+0.209, +0.394] |
| orion_hct116 | +0.601 | +0.267 | [+0.169, +0.364] |
| orion_hek293t | +0.627 | +0.287 | [+0.184, +0.385] |
| kolf2 (production operator) | +0.314 | +0.500 | [+0.312, +0.701] |

158's weaker key (rank of donor cosine) read +16 %/+30 % on the thinned Orion surfaces. pds is invariant to a strictly
positive row scale (158 measured the pre-render path at ±0.0006).

## Why this is worth a slot now (158 was closed for reasons that no longer hold)
158 was closed because (1) mse was clamped at the implied cosine (stale since 09-04: we score +0.187 scaled) and
(2) the welded nmae paid more than a clamped mse collected. This arm is still WELDED: the DE channel reads the scaled
rows. It measures the one thing no local surface can: **does the per-target key transfer to A/B/C mse**, and what the
welded DE members charge for it. If mse transfers, the next arm unwelds (ballast cells: a few 1e6-count cells set the
pseudobulk while carrying 1/400 of the per-cell mean and rank weight).

## Predictions (395's raw members: pds 0.854738, mse 0.809898, nmae 0.937101, fid 0.514228, reach 0.283909, jac 0.040891)
- **pds**: within ±0.002 raw of 0.854738 (falsifier F1: a larger move means the render is not a row scale).
- **mse (primary)**: raw Δ −0.03 … −0.07 (aggregate cosine +0.27…+0.30 at our 0.107 via 197's slope −2.0745 per unit
  cosine, discounted for an unknown PLOLO→board factor for this arm class).
- **nmae**: worse (welded; 158: −0.0023 … −0.0080 overall). fid / reach / jac: sign unknown (call counts move with row
  magnitude through the emitter's pass-1 |z| multiset; the floor targets lose calls, the cap targets gain them).

## Bins (primary = raw mse vs 0.809898)
| Δ raw mse | reading | next |
|---|---|---|
| ≤ −0.030 | the key transfers to A/B/C | build the unwelded (ballast) version; dose the cap upward |
| −0.030 … −0.010 | partial transfer | unwelded version at this dose; examine per-context if exposed |
| > −0.010 | no transfer on A/B/C | close the per-target key; do not build ballast for it |
Overall vs 0.2239363 is reported, not decisive: the welded DE cost is part of what is being measured.
Standing: the board shows the latest entry; a lower displayed row is not corrected. Never on rank.

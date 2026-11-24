# 040 — Union-axis gamma calibration on the faithful H1 test harness

> **Correction, 2026-08-22:** the H1 comparator quoted below is quarantined.
> The old runner excluded the 2026 target panel before applying the dense
> whitening transform, rather than excluding the actual H1 target panel.
> Eighty of 85 H1 self-target coordinates could therefore mix into otherwise
> scored coordinates. No H1 truth entered the covariance fit, and the direction
> of the bias is unknown, but the validation was not production-faithful.
>
> The four files currently under `out/runs/` are **H1** predictions
> (`52,400 x 18,080`), despite `run.log` saying “HCT116 excluded.” That message
> was hard-coded in the runner and was misleading. A cloud scorer rejected an
> attempted HCT116 comparison immediately because the real file has 18,533
> genes. These files inherit the H1 target-mask mismatch and must not be used.
> The interrupted local logs and rejected cloud attempt are not results.

`v9` changes the whitening axis from the 7,221-gene donor intersection to the
17,530-gene donor union on the H1 test harness. Experiment 038 selected
`gamma=0.65` from LODO, but the production-sized union operator has only one
count-level measurement.

This sweep holds every other `v9` choice fixed and tests the broad neighborhood
around 0.65 on `h1-test`: 85 perturbations, 400 perturbed cells per target, and
18,400 controls. The existing `gamma=0.65` result is the comparator:

`avg6=0.1895, pds=0.8655, nmae=0.9561, reach=0.2373, fid=0.5061, jac=0.0638`.

No challenge submission should change from `v9` unless a candidate beats that
result by more than scorer/seed noise and remains consistent with LODO.

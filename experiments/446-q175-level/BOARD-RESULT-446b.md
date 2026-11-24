# 446b BOARD RESULT — 0.2700428, NEW BEST (+0.0022834 vs 433k); the levelled q 1.5 → 1.75 step transfers at the calibrated rate

Entry `UDgjOKL9gI1udAexq6ue`, submitted 2026-09-26 00:01Z by `arm0926_s1.sh` ("the autoclave learned to knit",
description empty). Upload 22 min; scored and captured 00:37Z, rank 7. Package `446b-…-sharpen175-stackload-ctxlevel-
mde010coarse-zr5.vcc` sha256 `e171cb39…`. Reference 433k `T76JWD3tOajk5IYvzwUW` 0.2677594. Bins: `PREREGISTRATION.md`.
Reader: `read0926.py`.

| member | Δ raw vs 433k | Δ scaled | Δ overall |
|---|---:|---:|---:|
| pds | **+0.003170** | +0.007121 | +0.001187 |
| mse | **−0.007642** | +0.008027 | +0.001338 |
| nmae | +0.001346 | −0.002267 | −0.000378 |
| fid | +0.000012 | +0.000041 | +0.000007 |
| reach | +0.000700 | +0.000786 | +0.000131 |
| jac | −0.000003 | −0.000008 | −0.000001 |
| **overall** | | | **+0.0022834** |

## Against the bins
- **pds +0.0032 raw: "transfers" (≥ +0.002). mse −0.0076 raw: "transfers" (≤ −0.005).** Overall +0.0023, inside the
  predicted +0.001 … +0.004. 433k's record priced exactly this step at +0.0021 (438's grid at 417a's calibrated k).
- **nmae +0.0013 raw is just outside the "pinned" bin (≤ 0.001). Diagnosed (00:45Z) from the two builds' output
  audits: it is the ballast, not a carrier leak.**
  - The sharper render needs more ballast mass: median ballast cells per perturbation rose from 4 to 6
    (`rc_callset/ballast_K`), and the median ballast ratio fell from 0.73 to 0.58.
  - With two fewer carrier cells, each perturbation's carriers are a fresh draw. `n_pred_pass1` median 4391 → 4374,
    moves 3952 → 3911.
  - The deep ballast cells enter the per-cell arithmetic means with K/400 of the weight. So a sharper bulk reaches
    the mean LFC a little.
  - The same sign showed on the previous q step (433h → 433k nmae +0.0008 raw). It is a small structural cost of q,
    not a defect: −0.0004 overall against +0.0025 in the bulk.
- fid, jac and reach are inside their bins.

## Consequences
- **446b is the incumbent. q 1.75 is the bulk default for A/B/C and October** (overall bin: ≥ +0.001 with pds ≥ 0).
  DEF-RUNBOOK's "THE DAY IN ONE PAGE" already passes 1.75. `def_pipeline.py`'s SHIPPED_ENV q follows, with a new row-1
  known answer, after chain0926's builds.
- The q curve stops here. 438 prices q 2 at no calibrated gain over 1.75 (A/B/C +0.0003; `real`/`efp` −0.0005 / −0.0002),
  and the ballast's nmae cost grows with q.
- The 09-27 rule is in **Branch 1** (446b − 433k +0.0023 ≥ +0.001 and pds +0.0032 ≥ −0.0005). 1a or 1b follows from
  446bja's read.

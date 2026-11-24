# 415 — a production leave-one-donor-out per-target key (mandate, 2026-09-22 19:35Z)

**Current state (supersedes docs/READ-FIRST.md):** read `experiments/410-pertscale-revisit/RESULTS.md`, `BOARD-RESULT-410.md`,
`BOARD-RESULT-410b.md`, `key410.py`, `check410.py`, `build_weights410.py`. 410's key
`key_p = mean_{i<j} ½(<W d_i,d_j> + <W d_j,d_i>) / |W mean d|²` over the three production donors (raw rows, W = the
fold/production operator without its preserve constant) estimates the mse-optimal per-target scale of the whitened donor
BASE. It transferred on the board (raw mse −0.080 at clip(0.3,3)). But the shipped bulk row also carries a CD4 residual
(~35 % of the row norm) and an H1-union residual (77 targets). Adding CD4/KOLF2 as extra PAIRS in the key destroyed it
(key2_410.py) — because the fold prediction had no CD4 in it.

**Your territory (others: 413 scans for data, 414 re-prices bulk arms — do not do their work):** build
`key_p = mean_j <P^{−j}_p, d_{j,p}> / |P_p|²`, where P is the production-composition bulk row (whitened base + CD4
residual + H1 union residual, as experiments/228-theta-gain-plolo/six228.py `build_final` / arms228.py builds it on PLOLO
folds) and P^{−j} is that row rebuilt WITHOUT donor j (d_j = donor j's raw row), i.e. an unbiased estimate of <P,T>
covering the CD4/H1 components. Validate on PLOLO (162's folds; three cancer folds; kolf2 separately) against 410's key
evaluated on the SAME production-composition fold predictions: Spearman with the oracle s*_p = <P_p,T_p>/|P_p|², relative
aggregate transfer cosine at clip(0.3,3) and clip(0.3,∞) (energy-preserving), 1,000× target bootstrap CI.
**Known answer first:** reproduce 410's key on the folds (rho +0.576/+0.601/+0.627; gains in check410.json) with its
original fold predictions before switching to the production composition.
**If and only if** your key beats 410's on all three cancer folds, write the panel artifact in 410's exact format
(`targets, w, key, carried, n_donors` npz + json sidecar with sha256; A/B/C panel, production operator e2f65ecf, the
Oct-22 axis guard from build_weights410.py) for clip(0.3,3) and clip(0.3,∞) — do NOT touch the builder.
**Deliverable:** `experiments/415-lodo-key/RESULTS.md` with the comparison table and a one-line recommendation.
≤ 8 GB RSS (check `memory_pressure` first; builds of ~20 GB run tonight); no builder edits, no vcc, no cell-eval2; write
only in experiments/415-lodo-key/.

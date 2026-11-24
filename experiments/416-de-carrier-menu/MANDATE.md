# 416 — rank DE-carrier designs under ballast, with seed-averaged H1 six-member scores (mandate, 2026-09-22 20:10Z)

**Current state:** read `docs/FINDINGS.md` §0–§1, `experiments/410-pertscale-revisit/RESULTS.md` (§5, §7),
`BOARD-RESULT-410b.md`. Since 410b the bulk (pds, mse) and the DE channel (fid, reach, jac, nmae) are separately
designable: ballast cells pin the pseudobulk while the 400−K carrier cells alone set the DE table. Our DE members are now
the gap (board nmae 0.936 vs a pack at 0.88–0.91 with similar pds/mse). Two H1 facts to respect: (1) reach moves ±0.013
raw between K-class emission REALISATIONS (drop1 control), so single-seed H1 reach deltas under ~0.015 are noise;
(2) nmae/fid/jac realisation noise is small (≤ 3e-4).

**Your territory (others: 413 scans for new data, 414 re-prices BULK arms — do not do their work; do not touch the bulk):**
using the H1 six-member harness `experiments/410-pertscale-revisit/run410b.py` (ballast, ratio mode, jackknife boost;
its `base395` known answer is `work/h1/metrics-base395.json`), score CARRIER variants with the bulk held at the incumbent
(eff_bulk = the shipped cap15 effects), each at THREE emission seeds (vary the RNG seed passed to `b.make_rngs`; keep the
bulk identical) so realisation noise averages out. Copy the harness into your own directory; never edit 410's files.
Carrier designs (one variable each vs the base carriers):
  (a) base carriers, 3 seeds — the noise floor;
  (b) whitening strength of the carrier design: gamma 0.0 (unwhitened), 0.3, 0.65 (base), 0.9 — norm-matched rows
      (the pre-whitening H1 design is in `experiments/323-current-family-validation/work/stage.npz`, operator
      `whitener_h1_pruned.npz` there; see `run323.py`);
  (c) 228's per-gene gain on the carriers (gain_reg, t 0.25; `experiments/228-theta-gain-plolo/work/gain-kolf2.npz`);
  (d) the DE sign prior OFF (K-class carriers without `lfc_target`) and variants of its protect quantile (q 0.7, 0.9);
  (e) any other carrier idea you can justify from FINDINGS §1 or the 38x–40x experiments, at most two.
Report per design: mean ± sd over seeds of all six members (pds/mse must stay at the base's values — that is the
ballast falsifier), and the predicted overall at board slopes (pds +2.2138, mse −1.051, nmae −1.6532, fid +3.3645,
reach +1.1292, jac +2.7263 per raw unit; divide by 6). Rank designs by the DE-member sum; flag any whose gain exceeds
2 sd of the seed noise.

**Deliverable:** `experiments/416-de-carrier-menu/RESULTS.md` with the table, the noise floor, and ONE recommended
carrier change for a board arm (or "none clears noise"). H1 arms run strictly serially (~5 min, ≤ 12 GB each); check
`memory_pressure` (≥ 40 % free and free+inactive ≥ 20 GiB) before each; a ~20 GB production build may be running.
No builder edits, no vcc, no cell-eval2; write only in experiments/416-de-carrier-menu/.

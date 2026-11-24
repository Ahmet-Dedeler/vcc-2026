# 414 — re-price the BULK arms that died on nmae, now that ballast unwelds the DE channel (mandate, 2026-09-22 19:35Z)

**Current state (supersedes docs/READ-FIRST.md):** read `experiments/410-pertscale-revisit/BOARD-RESULT-410b.md`,
`BOARD-RESULT-410.md`, `RESULTS.md` and `PREREGISTRATION-410c.md` / `-410d.md` there. 410b proved on the board that ballast
cells unweld the DE channel: pds/mse read the BULK (pseudobulk), fid/reach/jac/nmae read the carriers, and the carriers
stay exactly as 395 builds them. So any BULK change now costs nothing on the DE members. Many bulk arms were closed only
because "nmae pays" (the weld): 158 (per-target scale), 224/228 (per-gene reliability gain), 169 (sharpening), 187
(amplitude / cliff), 318. Tomorrow's board arms are 410c (410b with the UNCAPPED per-target key, `wshare_lo0.3_inf`) and
410d (410c + 228's gain_reg at t 0.25 on the final bulk row).

**Your territory (others: 413 scans for new data, 415 builds a leave-one-donor-out key — do not do their work):**
price bulk-only designs on PLOLO (experiments/162-retrodicting-validator/plolo.py, cached folds in its work/) with the
PRODUCTION composition (whitened base + CD4 residual + H1 union, as 228 built it: experiments/228-theta-gain-plolo/
arms228.py `price` / six228.py `build_final`) and the render read-back (plolo.deliver / 228's path), three cancer folds
(k562_gwps, orion_hct116, orion_hek293t), kolf2 reported separately:
  (a) 228's gain (gain-<fold>.npz, gain_reg) at t ∈ {0.1, 0.25, 0.5, 1.0}, applied (i) to the whitened base BEFORE the
      residuals (228's way) and (ii) to the FINAL row (410d's way), each on top of 410's per-target key (key410.py,
      clip(0.3,∞), energy-preserving — see experiments/410-pertscale-revisit/check410.py);
  (b) 169's sharpening and any other magnitude reshape closed on nmae (find them in 158/169/187/224/318 RESULTS.md);
  (c) amplitude: with the key-weighted bulk, where is the mse-optimal global / per-context scale? Use 187's per-context
      cliff model (experiments/187-joint-amplitude/work/cliff187.json; MSE_c(s)−1 = [s²A_c + K_c − min(claim_c, s²B_c)
      − 2sX]/ΣD) with X raised by the key's measured aggregate-cosine gain. Does the optimum move above context A's
      cliff (α 0.9823)? Per-context recommendations with predicted Δ raw mse.
Report pds (rendered), relative aggregate cosine, predicted Δ raw mse via 197's slope (raw mse = 1.0236 − 2.0745·cos;
the board's measured cosine for 395 ≈ 0.107 and 410b's raw mse is 0.730994), and predicted Δ overall (scaled slopes:
pds +2.2138, mse −1.051 per raw unit; divide by 6).

**Known answers first:** 410's key gains (check410.json: clip(0.3,3) +0.299/+0.267/+0.287) and 228's gain_reg t 0.25 rows
(experiments/228-theta-gain-plolo/work/arms228-*.json) must reproduce before any new number.

**Deliverable:** `experiments/414-bulk-reprice/RESULTS.md`: ranked table of bulk designs with the numbers above, and ONE
recommended next bulk arm as a single variable on 410c. ≤ 8 GB RSS (check `memory_pressure` before heavy steps: builds of
~20 GB run on this Mac tonight); no builder edits, no vcc, no cell-eval2; write only in experiments/414-bulk-reprice/.

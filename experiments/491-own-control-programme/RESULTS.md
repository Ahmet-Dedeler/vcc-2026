# 491 — The recipient's own control cells as the DE block-programme source: real but small on H1. No arm.

2026-09-26 22:40Z, third session (analysis only; no build, no `vcc`, no `cell-eval2`). Scripts: `frr491.py` (491a),
`pcs491.py` (491b). Outputs: `work/frr491.{json,log}`, `work/pcs491.{json,log}`.

## Question
160 found a line's own control subspace holds 0.71–0.89 of its MEAN perturbation response (Orion, closed then only for
pds). The DE channel's block prior is exactly a mean-response sign programme, today imported from other lines (CLEAN4)
or, on 09-30, from Arc's H1 bank (458d). If the recipient's controls carry the programme, every context, including
D/E/F, would get a same-line, same-lab programme for free. 441 had tested single control PCs (≤ 0.18 correlation, PC2
gate-block 0.557) and the identity axis (ps6: fid-like +0.0014); it had not tested the subspace capture on H1 nor
transporting CLEAN4 through the recipient's covariance (projection onto the top-k control PCs; the reweighting /
fluctuation-response form Cov_ctrl(y_g/mean_g, h) with h the cells' score on CLEAN4).

## Known answers (all reproduce)
- CLEAN4 fid-like 0.5244 and LOTO oracle 0.5544 (462f); gate-block 0.5799 and oracle 0.6455 (441 §2).
- Three targets' lfc recomputed from their own cells against my control CPM equal the reference table to 1.1e-5.

## Results (h1-test, 85 targets; fid-like = sign agreement with every target's lfc on the 2,200 block genes)
| quantity | value |
|---|---|
| share of H1's mean-response energy inside the top-k control PCs, all genes / block genes | k 20: 0.20 / 0.08; k 50: 0.23 / 0.085; k 200: 0.25 / 0.09 (random k/G: 0.005 at k 50) |
| best own-control programme, fid-like (all controls) | `proj_k2 + CLEAN4` 0.5300 vs CLEAN4 0.5244: **+0.0055** [+0.0012, +0.0095] target bootstrap, P(≤0) 0.008 |
| same, programme from control half A, truth vs half B | 0.5254 vs 0.5230: **+0.0024** |
| gate-block (significant block genes) | +0.0019 [−0.010, +0.013] |
| reweighting / FRR forms (`frr_lin`, `frr_z`) | correlate 0.27 / 0.22 with the mean response (CLEAN4 0.18) but fid-like ≤ +0.003 |
| by n_real tercile (fid-like) | weak +0.007, mid +0.004, strong +0.006; the weak tercile's gate genes are only 163 in total, so its gate-block 0.43 → 0.58 is ~10 genes |

What the axis is (491b): H1's control PC1 is the G2/M cell cycle (UBE2C, AURKA, … vs histones); PC2 is S-phase/E2F
(DTL, CCNE1, CDC25A, UHRF1, MCM10, E2F1) against a hypoxia/glycolysis stress set (DDIT4, SLC16A3, EGLN1, PDK1, P4HA1,
ACSS2). H1's mean response lies along −PC2 (block r −0.28): perturbation lowers S-phase/E2F genes and raises the
hypoxia/glycolysis set. H1's PC1 matches A/B/C's PC1 (|cos| 0.81–0.85); H1's PC2 matches A's PC2 at 0.45 and is spread
over several PCs in B and C (≤ 0.28 each).

## Reading
- In the challenge's own chemistry (Flex), the control subspace holds only ~9 % of the block-gene mean response at k 50,
  not 160's 71–89 % (Orion, 10x 3′, where a nuclear-fraction blob is a shared technical axis). The own-control route has
  a low ceiling on H1.
- The usable part is one biological axis (E2F/S-phase vs hypoxia-glycolysis stress). It is worth about +0.002–0.006
  fid-like on H1, the same class as 441 ps6's identity axis (+0.0014). At the recorded H1 → board fid transfer
  (0.07–0.67×, `docs/H1-BOARD-CALIBRATION.md`) that is ≤ +0.002 overall, with the k choice (2) fitted on H1 and a
  different PC order in B and C.
- **No arm.** If 458d reads positive on 09-30, a same-lab programme will dominate this. If a future arm needs a
  recipient-side programme for a line with no bank, `proj_k20` (fid-like +0.002, less k-fragile) is the form to start
  from, and the halves check is mandatory.

## Memory note
`frr491.py` peaked at 5.7 GB RSS for ~20 s, over the 4 GB cap agreed with the peer session for tonight; `pcs491.py`
peaked at 2.9 GB. Later runs keep one float32 matrix at a time.

## Addendum 491c (22:50Z) — combined recipient-side programme: CLEAN4 + own-control transport + identity axis
`combo491.py`, `work/combo491.{json,log}`. z-scored on the block and summed; identity = 441 ps6's I (H1 controls' log mean
CPM minus the A/B/C contexts' mean, same chemistry; its block correlation with H1's mean response is +0.28 for −I).
Target bootstrap 4,000 draws (Δ vs CLEAN4, [2.5 %, 50 %, 97.5 %]):

| programme | fid-like Δ (all controls) | fid-like Δ (halves) | gate-block Δ |
|---|---|---|---|
| + proj_k2 | +0.0054 [+0.0011, +0.0096] | +0.0024 [−0.0016, +0.0063] | +0.0016 [−0.011, +0.012] |
| + identity (λ 2) | +0.0034 [−0.0004, +0.0069] | +0.0036 [−0.0001, +0.0070] | +0.0134 [−0.007, +0.030] |
| **+ proj_k2 + identity (λ 1)** | **+0.0064 [+0.0025, +0.0101]** | **+0.0040 [+0.0001, +0.0076]** | **+0.0131 [−0.003, +0.026]** |
| + proj_k20 + identity (λ 1) | +0.0018 | +0.0009 | +0.0052 |

- The two recipient-side axes add. The combination is the best recipient-side programme measured on H1: fid-like +0.004
  … +0.006, gate-block +0.013.
- At the recorded H1 → board transfer for sign priors (fid 0.07–0.67×), that is ≤ +0.002 overall. It stays a rider-class
  lever, and k 2 is fitted on H1.
- **October use, if ever:** for a context with no same-line bank, and only if 458d's 09-30 read says CLEAN4 remains the
  block prior. The identity axis for D reads D's controls against E, F and A/B/C. Not an arm now.
- Memory: 4.38 GB peak for 8 s (the centred copy in `proj`), again just over the 4 GB cap agreed for tonight.

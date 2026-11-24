# 441 — Where the DE channel's sign gap lives on H1, and what can close it (2026-09-24, coordinator, analytic, no build)

Surface: h1-test truth (`323/work/reference.npz`, 85 targets; gate = q < 0.05, finite lfc, own gene out, targets with gate
≥ 10 → 58), our H1 carriers (`426/work/h1/stats-base_s1_bps_mde0.1cap0.1min0.3coarse.npz`, i.e. the shipped coarse MDE;
`391/work/stats-{cap15,k3c5b}.npz` for design vs CLEAN4-prior carriers). Accuracies pooled over (target, gene) unless
stated. Motivation: the board leaders pinned their bulk around 09-17 and have since moved only DE members (cqawesome
nmae 0.869 → 0.829, fid 0.508 → 0.535, reach 0.319 → 0.334 with pds/mse frozen to 4 decimals; nailong's 09-23 step
nmae −0.025). Our nmae is 0.920.

## 1. Our nmae is sign-limited, not magnitude-limited
| true \|lfc\| bin | share of gate \|l\| mass | our sign acc | our mean \|lfc\| delivered |
|---|---:|---:|---:|
| 0–0.25 | 21 % | 0.564 | 0.10 |
| 0.25–0.5 | 33 % | 0.555 | 0.14 |
| 0.5–1 | 29 % | 0.549 | 0.13 |
| 1–2 | 12 % | 0.604 | 0.14 |
| 2–4 | 4 % | 0.703 | 0.27 |
Delivered magnitude is flat (~0.13) across bins, which is L1-optimal for a ≈ 0.55 (the MDE rule). More magnitude needs
better signs.

## 2. The sign gap is on the genes we do not call, and the context's own mean response would close it
Gate sign accuracy (unweighted / |l|-weighted):
| subset | ours (design) | ours (+CLEAN4 prior) | CLEAN4 alone | **H1 leave-one-target-out mean response (oracle)** |
|---|---|---|---|---|
| called (25k gate genes) | 0.567 / 0.626 | 0.598 / **0.649** | 0.578 / 0.571 | 0.639 / 0.640 |
| uncalled (40k) | 0.536 / 0.540 | 0.536 / 0.541 | 0.568 / 0.571 | **0.616 / 0.630** |
| block (2,200 top-CPM) | 0.553 | 0.587 | 0.580 | **0.645** |
| all gate | 0.548 / 0.562 | 0.560 / 0.568 | 0.573 | **0.625 / 0.632** |
Per-target means: ours 0.560 / 0.580 weighted; the oracle mean response 0.612 / 0.618. The leaders' inferred gate sign
accuracy (434f: ~0.62–0.66) is exactly this oracle's level. On the genes we call our signs already beat it (0.649).

## 3. Donor routes to the context's mean response do not get there
| predictor of the context mean | gate-block | gate-nonblock | all block genes (fid-like, vs ref lfc sign) |
|---|---:|---:|---:|
| H1 LOO mean (oracle) | 0.645 | 0.614 | 0.570 |
| CLEAN4 (4 cancer banks, all perturbations) | **0.580** | **0.567** | **0.535** |
| K562-GWPS all-perturbation mean | 0.570 | 0.539 | 0.531 |
| K562-GWPS mean over the H1 panel's own targets (LOO) | 0.562 | 0.526 | 0.528 |
| control-cell PCs (real H1 NTC cells, log1p CPM, best of top 10, orientation oracle) | 0.557 (PC2) | — | — |
- Restricting a donor's mean to the panel's targets makes it worse, not better: the context mean is context-specific,
  not panel-specific.
- The control cells' top PCs correlate ≤ 0.18 with the mean response (PC1 +0.06). 160's "the control basis captures
  the mean response" does not hold on H1 as a sign source.
- Delivering CLEAN4 signs on uncalled genes through the MDE is 426's closed route (class over-fill; analytic prize only
  −0.003 raw nmae).

## Verdict
The DE gap to the leaders is the context's own shared response programme on uncalled genes (+0.05–0.08 sign accuracy),
and nothing we hold estimates it for an unseen line better than CLEAN4. Probing it on A/B/C through the board would not
transfer to D/E/F (their mean response is unobservable). No arm. If a future source gives cancer-line mean responses
closer to the recipient (e.g. many-line CRISPRi atlases), this table is the instrument to price it: the prize is the gap
between the CLEAN4 row and the oracle row.

## Postscript — is there a new many-line source? (2026-09-24 21:45Z, web + API checks)
- **X-Cell / X-Atlas-Pisces (Xaira): still nothing downloadable.** HF model `Xaira-Therapeutics/X-Cell` lastModified
  2026-03-16, 272,925 bytes (README + PNG, no weights); GitHub last commit 2026-03-16; HF dataset
  `X-Atlas-Pisces` lastModified 2026-03-25 (README/LICENSE only). 205's verdict stands; `scripts/pisces_watch.sh` keeps
  polling (Pisces would add genome-scale Jurkat = same-line for context A).
- **The 100M-cell genome-wide CRISPRi/CRISPRa atlas "across dozens of human cell types" (2026, ScienceDirect P588) was
  generated on Illumina's own PIP-seq kit.** It is Illumina's data, not public. This is the most likely source of the
  board leader "Illumina AI" (model "PerturbationAI") and its all-member lead (pds 0.888, mse 0.556, nmae 0.810,
  reach 0.430): many-line same-chemistry-class data is exactly the oracle row of section 3.
- Newly published but already on disk or pluripotent-only: the KOLF2.1J genome-scale iPSC atlas (Nat Biotech 2026),
  the 34-line iPSC trans map, VIPerturb-seq (424: negative as a K562 donor).

## Postscript 2 — a line-specific signal exists but is weak (22:15Z)
- **"Identity genes go down."** Genes expressed relatively higher in the recipient than in the three donor lines
  (rel = log pseudobulk at TS 5e4, H1 controls minus the mean of K562/HCT116/HEK293T `b_ctrl`) tend to fall in the
  recipient's mean response: fraction down 0.47 in the bottom decile of rel → 0.68 in the top decile; Spearman −0.128
  over 10,713 genes. Absolute expression carries nothing (Spearman −0.001).
- As a sign source it is weak: −rel alone gives gate sign accuracy 0.545 (all genes, no coverage gaps), CLEAN4 0.573 on
  its 57 % coverage. Combined (z-scored CLEAN4 + λ·(−rel), λ 0.25–1): gate 0.564, uncalled 0.553, block fid-like
  0.5250–0.5256 vs CLEAN4's 0.5244. No arm; the table stays as the instrument for any future recipient-side feature.

## Postscript 3 — choosing the CPM block by prior confidence: fid up a little, jac down more (20:20Z)
Emulated on H1 (head = top 400 by |delivered lfc|, block = 2,200 by a key, CLEAN4 flips on the unprotected 80 % of the
block, design signs elsewhere; the emulation reads fid 0.5079 / jac 0.0728 against the delivered 0.5097 / 0.0747):
| block key | fid | precision | jac |
|---|---:|---:|---:|
| control CPM (shipped) | 0.5079 | 0.5232 | 0.0728 |
| CPM rank among CLEAN4 conf ≥ 0.02 | 0.5101 | 0.5252 | 0.0565 |
| CPM rank among conf ≥ 0.05 | 0.5133 | 0.5284 | 0.0457 |
| log CPM + 20·conf | 0.5118 | 0.5268 | 0.0623 |
| log CPM + 40·conf | 0.5133 | 0.5284 | 0.0547 |
fid +0.002…+0.005 raw (+0.007…+0.018 scaled) against jac −0.005…−0.027 (−0.014…−0.07 scaled): net negative. The block
is chosen by control CPM because jac needs power; that choice stands.

## Postscript 4 — withholding the MDE from targets whose signs are worse than chance (21:45Z)
On H1 (shipped coarse MDE vs no MDE, seed 0), 8 of 58 gated targets have |l|-weighted gate sign accuracy < 0.5 and the
MDE hurts 14 of them (Spearman(acc, ΔMDE) −0.62). The ORACLE gate (MDE withheld where acc < 0.5) moves the mean ΔMDE from
−0.00966 to −0.01072 raw: −0.0011, about +0.0003 overall at the board's 1.76×. Any real predictor gets less. Not pursued.

## Postscript 5 — control-cell cell-cycle and cell-size axes are not the mean response either (21:55Z)
From H1's own NTC cells: a G1-shift axis (mean log CPM of the 40 % least-cycling cells by Tirosh S/G2M scores, 42 + 52
markers found, minus the mean of all cells) and a small-cell axis (40 % lowest libraries minus all). They correlate
+0.029 and +0.078 with H1's mean response and give gate sign accuracy 0.514 and 0.515 (64,940 gate pairs). The
"perturbations push cells toward G1 / smaller" hypothesis does not supply the recipient's shared programme on H1.

## Postscript 6 — the identity axis, measured in the SAME chemistry, is the best recipient-side signal (22:30Z)
The July 2026 preprint "Perturbation response decomposition…" (xinyizhanglab; repo already in `data/external/
perturbation-decomposition`, fig1_dsup) reports that the cell-line template (mean response) lies partly along the
control IDENTITY axis, I_c = mean ctrl_c − mean ctrl_others, on four same-chemistry essential screens. Postscript 2's
identity test used cross-chemistry donors (10x 3′ K562/HCT116/HEK293T vs Flex H1). Redone within the challenge's own
chemistry (H1 Flex controls vs the A/B/C Flex controls, log mean CPM):
- corr(H1 mean response, I_H1) = **−0.236** (Spearman −0.244; cos² 5.6 %), against −0.128 cross-chemistry.
- As a sign predictor, −I_H1 gives gate sign accuracy 0.556–0.568 on ALL genes. CLEAN4 gives 0.573 on its 57 % coverage.
- Combined, z(CLEAN4) + λ z(−I) with λ 1–2: gate-block **0.580 → 0.592**, gate 0.573 (now full coverage), block-all
  (fid-like) 0.5244 → 0.5258.
- Reading:
  - a real, recipient-side, zero-shot-computable signal ("identity genes go down under perturbation");
  - worth ≈ +0.01 sign accuracy on the significant block genes (reach, nmae) and ≈ nothing on fid;
  - for A/B/C it would read ctrl_c − mean(other challenge contexts), for D/E/F ctrl_D − mean(E, F, plus A/B/C);
  - not built: expected board value ~+0.001. Delivery on uncalled genes is still 426's closed route.

## Postscript 7 — after 433h's board read: strong donor knockdowns are real in the recipient, but already delivered (00:50Z 09-25)
433h's zero repair gained nmae −0.0247 on the board. Most of 428's 24 zero-total genes (22 trans, "donor-row
knockdowns") must have sat in A/B/C reference gates, which is far above the ~1 % base rate. On H1:
- **Strong non-cis donor knockdowns replicate at 20–26 %** against a 3.7 % base rate: K562-GWPS rows 0.202 (n = 367),
  Orion HCT116 0.263 (n = 114), HEK293T 0.25 (n = 12). The strongest (< −1 log units) replicate every time (n = 7).
- The K562 rows (which may share protospacers with Arc's library via the common CRISPRi parent collection) do NOT
  replicate better than Orion's (a different library). No evidence for shared-guide off-targets; this is conserved
  biology.
- **Already delivered:** among gated strong-knockdown genes our delivered lfc median is −0.34 against a true −0.47, with
  sign 0.76. Deepening them (to −1, −1.5, or the donor minimum in log2) buys at most −0.0007 raw nmae on H1.
- The repair's board gain came from removing catastrophic −35 cells, not from magnitude in general.

# 473 — PREREGISTRATION: 410d2 (228's per-gene reliability gain m, render-gated, clipped) re-priced on 446b's bulk

Written 2026-09-26 22:40Z by the third session at the coordinator's request, BEFORE any number on the 446b base exists.
Analysis only: no build, no builder edit, no `vcc`, no `cell-eval2`. The coordinator keeps the decision and any build.

## Question
414 priced 228's per-gene gain m, amplitude-matched and bounded, at **+0.0025 (kolf2-like) … +0.0046 (cancer)** overall
on 410c's bulk (`price410d.json` `matched_overall`), then shelved it. The bulk has since moved: stack H1 ×2 + CD4 ×1.5
(`VCC_BULK_LOAD` of `dump_h1x2_cd4x15.npz`), sharpening q 1.75, per-context levelling. What is m worth as an INCREMENT
over 446b (0.2700428), and is it render-safe?

## Inputs (sha256)
- 446b's bulk: `experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz` `bb19f0e5…` (the bulk after the K3 key; 446b's
  meta `bulk_load_sha256` is the same).
- 446b's carriers (PP composition, unkeyed): `experiments/417-bulk-sharpen/work/dump_ka417a_g065.npz` `aed59b0b…`
  (`carriers`; 446b's own composition env is PP: cd4_target_f 0.34, h1_union_g 0.28).
- m: `experiments/410-pertscale-revisit/work/bulk_gene_gain_g228reg_t025.npz` `da937cfa…` (the builder's pinned
  `g228reg_t025` artifact; on PLOLO folds the fold-specific 228 gain at t 0.25, as 414 did).
- 446b's measured render energies: `experiments/446-q175-level/work/ntclevel-446b.json` `71e2e6d2…` (A, B, claim per
  context); factors A 1.1004, B 1.1517, C 1.1720; q 1.75; K3; ballast KMAX 16.

## Arms (one variable each vs the base `stack|K3+sh1.75`, 446b's bulk)
- **A1 (primary):** m clipped to [1/2, 2], multiplied into the keyed bulk BEFORE sharpening, one constant restoring the
  energy on the rendered (gated) non-panel genes summed over the three contexts (the builder's `…rg` rule), then
  sharpening q 1.75 (row renorm) and the per-context factor.
- A2 (reference): unclipped m, same placement and constant (410d2 as it would be built on the load path).
- A3 (reference): m clipped to [1/3, 3].
- A1post (reference): the clipped m applied AFTER sharpening with each row's gated norm restored (a different builder
  placement; sharpening at q makes a pre-sharpening m act like m^q on the sharpened shape).
- C1 (control): A1 with m permuted over genes. Must lose.

## Instruments
1. **Render safety on the real A/B/C panel** (`abc473.py`, light): 414/420's realised route (plolo `profile_from_effect`
   minnorm, 414's ballast emulation at KMAX 16 with the carriers' design profiles) on the dump. Per context: render-space
   energy ratio arm/base (A_c, B_c over non-panel genes), max per-target energy ratio, the largest single-gene mass share
   of any block, SYVN1's ratio. **Known answers:** the base's realised A_c/B_c equal 446b's measured `ntclevel-446b.json`
   values within 1 % (the realised route's validated accuracy is 0.4–0.7 %, 420 §2), and the emulated ballast K
   distribution matches 446b's meta (median 6, max 16).
2. **Direction on PLOLO folds** (`arms473.py`, one fold per process, ~6.5 GB, run after 04:30Z or behind a vm_stat
   free+inactive ≥ 33 GiB gate): lib420's `stack|K3` composition, the fold's 228 gain. **Known answers first:**
   `stack|K3+sh1.75` designed and realised rows == `experiments/420-bulk-stack/work/comp420-<fold>.json` exactly; 414's
   `Kinf+gf@0.25` == `arms414-<fold>.json` (my gain code). Per-target arrays saved for the bootstrap.
3. **Price** (`analyse473.py`): 420's cliff-aware route re-anchored at 446b's measured energies; each arm and the base
   re-levelled to its own EV-optimal per-context factor (187's cliff model, `level_opt`), overall =
   (2.2138 Δpds_board-emulated − 1.051 Δraw mse)/6, k ∈ {1, 1.407}; cancer folds and kolf2 reported separately;
   target bootstrap (4,000 draws) on Δcos_agg and Δpds per fold.
4. Also re-run 414's `price410d.py` pricing on its saved inputs: `matched_overall` must reproduce +0.0045474 (cancer)
   and +0.0025158 (kolf2) exactly (the coordinator's KA of 414's bounded-arm row). 414's emulation rows on 410c's bulk
   cannot be re-run: the 410c artifact (`the-kettle-applied-for-tenure.h5ad`) and its `.vcc` were swept; the real-panel
   known answer is item 1's, on 446b.

## Predicted range (before any 446b-base number)
- Render safety: clipped m keeps context energy within ±3 % (414: ×0.994 / ×1.002 / ×0.995 on 410c), but q 1.75 makes
  the pre-sharpening m act like m^1.75, so the max per-target ratio may approach 1.5. P(item 1 fails) ≈ 0.3.
- Direction: sharpening already concentrates energy on large-|h| genes, which overlap high-m (reliable, expressed)
  genes, so it absorbs part of 228's gain. **Predicted Δoverall vs 446b: +0.000 … +0.003 (cancer-like, k 1, re-levelled);
  kolf2-like −0.001 … +0.002.**

## GO rule (all four)
1. A/B/C render safety: energy ratio in [0.97, 1.03] in every context AND max per-target ratio ≤ 1.5 AND the largest
   single-gene mass share ≤ 2× the base's.
2. PLOLO cancer folds: Δoverall > 0 on ≥ 2 of 3 folds (k 1), and the pooled target-bootstrap 95 % lower bound of the
   relative Δcos_agg > 0.
3. Predicted board Δoverall ≥ +0.0015 (cancer-like, k 1, re-levelled) AND ≥ 0 kolf2-like.
4. C1 (permuted m) negative on ≥ 2 of 3 cancer folds.
All four → GO: a builder patch (m on the loaded bulk, clipped [1/2, 2], render-gated) and a build are the coordinator's.
Any failure → NO-GO, recorded with the failing item and its mechanism. A3/A2/A1post are reported, not selected.

# 410 — per-target bulk magnitude from a whitened shared-energy key: real on every surface, and the DE channel pays for the weld

Status: **PARTIAL** (PLOLO + H1 complete; board arm 410 building 2026-09-22 17:30Z; ballast arm 410b in preparation).

## Known-answer reproduction
- H1 harness (`run410h.py base395`, 400's k3c5bq80 construction via `run400.py`, read-only): reproduces
  `experiments/394-protected-block-signs/work/metrics-k3c5bq80.json` — emitted-count digest `5d4b3b82…` and panel-total
  digest `f699b3d5…` equal; all six members |Δ| ≤ 1.1e-16 (`work/h1/known-answer-base395.json`).
- PLOLO folds are 162's cached `fold-*.npz` with their registered operators (sha-checked by `plolo.Operator`).

## 1. The key (key410.py)
`key_p = mean_{i<j} ½(<W d_i, d_j> + <W d_j, d_i>) / |W mean_i d_i|²` over the donors that form the prediction (raw bank
rows, W without the on-axis preserve constant). It is the cross-fitted plug-in estimate of the mse-optimal per-target
scale `s*_p = <P_p, T_p>/|P_p|²`: a held donor stands in for the truth, so each donor's own noise drops out of the numerator.

## 2. PLOLO (deep pseudobulk truth, fold-held operators) — `size410.py`, `robust410.py`, `check410.py`
| fold | base agg cos | rho(key, s*) | gain clip(0.3,3) [95% CI] | gain uncapped | oracle |
|---|---:|---:|---|---:|---:|
| k562_gwps | 0.02702 | +0.576 | +0.299 [+0.209, +0.394] | +0.386 | +0.715 |
| orion_hct116 | 0.02121 | +0.601 | +0.267 [+0.169, +0.364] | +0.341 | +0.672 |
| orion_hek293t | 0.01868 | +0.627 | +0.287 [+0.184, +0.385] | +0.350 | +0.707 |
| kolf2 (production op.) | 0.00686 | +0.314 | +0.500 [+0.312, +0.701] | +0.936 | +1.685 |
Relative aggregate-cosine gains are unbiased under truth noise (the key never reads the truth; noise only rescales |T|).
158's rank-of-donor-cosine key reads +0.08…+0.29 on the same folds (`size410.json`, key `cos`).

**Negative side results (both recorded so nobody re-runs them):**
- *Target expression in the recipient line* (`expr410.py`): targets unexpressed in the held line (K562: 43 targets with
  CPM < 1) carry the same transferable signal (cos +0.023 vs +0.027…+0.034; same truth norms). An expression gate on top of
  the key costs −0.01…−0.03 of relative cosine on every fold. Consistent with 167. Closed.
- *Extra lines in the key* (`key2_410.py`): adding KOLF2 or CD4 (3 conditions) to the pairwise shared energy DESTROYS it
  (rho 0.58 → 0.09…0.51; gains go negative uncapped). The key must be built from exactly the donors that form the
  prediction — it estimates that prediction's projection on the truth.

## 3. H1, six members end to end, welded (`run410h.py`, `cmp410h.py`; Δ vs base395)
| arm | pds | mse | nmae | fid | reach | jac | board-slope overall |
|---|---:|---:|---:|---:|---:|---:|---:|
| ps_lo03_hi3 (key, energy-preserving) | +0.00462 | **−0.04858** | +0.02377 | −0.00138 | −0.01581 | −0.00054 | −0.00033 |
| ps_shuf (same weights, permuted) | +0.00126 | +0.03016 | +0.02257 | −0.00316 | −0.03033 | −0.00121 | −0.01906 |
- The key's information is worth **−0.079 raw mse** against the permuted control on H1, a line no PLOLO fold saw.
- **nmae pays the same with or without the key** (+0.024 vs +0.023): the cost is the energy-preserving shrink of the median
  target's fold changes (final scale median 0.59), not the assignment. That is 158's wall, and it is a WELD cost.
- Welded, the mse gain (≈ +0.0085 overall at board slopes) is eaten by nmae (−0.0066) and reach (−0.0030).

## 4. Consequence: unweld (410b, ballast cells)
The scorer's pseudobulk sums raw counts; per-cell CPM gives every cell 1/400 weight in the DE table. K cells of ~999k counts
set a block's per-gene totals while barely moving arithmetic means or ranks. `run410b.py` (H1) and `patch_ballast410.py`
(builder) implement: carriers = 395's DE design on 400−K cells, ballast = alloc(S, scaled profile) − carrier totals.
Jackknife margin on H1: our C_hat/C = 1.71…2.39 (min 1.706), so the carriers' leave-one-out share may fall to ~59% before
any mse credit is lost; ballast cells whose composition differs from the block total add their own variance.

## 5. Ballast on H1 (`run410b.py`; Δ vs base395; controls share the carriers' realisation)
| arm | pds | mse | nmae | fid | reach | jac | note |
|---|---:|---:|---:|---:|---:|---:|---|
| bal_ident (design-pinned, w=1, K=1) | −0.00588 | −0.00442 | +0.00000 | +0.00007 | −0.01253 | −0.00001 | pins bulk to the bare design |
| bal_ident_ratio (ratio mode, w=1) | −0.00056 | +0.00050 | +0.00002 | +0.00005 | −0.01253 | −0.00001 | bulk restored |
| drop1 (same 399 carriers, NO ballast) | +0.00028 | +0.00005 | +0.00002 | +0.00007 | −0.01408 | −0.00001 | realisation control |
| bal_lo03_hi3_ratio (key bulk, no jk boost) | +0.00560 | +0.00969 | +0.00026 | +0.00012 | −0.00801 | +0.00006 | jk credit lost |
| dropK_lo03_hi3 (same carriers, no ballast) | +0.00014 | +0.00019 | +0.00018 | −0.00004 | −0.00937 | −0.00001 | realisation control |
| de_med_lo05_hi18_ratio (CARRIERS key-scaled, bulk base) | −0.00224 | +0.03825 | +0.01300 | +0.00099 | −0.00247 | +0.00009 | DE-side key scaling hurts nmae |
| bal_lo03_hi3_ratio_sink_jk (target-gene mass sink, 60/85 used) | −0.00686 | −0.05107 | +0.00026 | +0.00015 | −0.00575 | +0.00006 | exact scaling of the read-back HURTS H1 pds |
| **bal_lo03_hi3_ratio_jk (SHIPPED 410b construction)** | **+0.00560** | **−0.04959** | +0.00026 | +0.00015 | −0.00575 | +0.00006 | C_hat/C min 1.50; board-slope overall **+0.0097** |

Readings:
- **Ballast unwelds the DE channel.** Against the same-carrier controls, nmae/fid/jac move ≤ 3e-4 and reach within its
  realisation noise. The reach deltas vs base395 (−0.008…−0.014) are the K-class emitter's realisation noise
  (drop1 = same 399 carriers without any ballast reads −0.014): **H1 reach differences below ~0.015 are not findings for
  these emitters** (the 0.0033 floor in AGENTS §8.1h was measured on the old quasi emitter).
- **Design-pinned ballast removes a real pattern.** The incumbent emitter delivers 0.5–100-count genes −0.002…−0.0035
  log units below the design (`emitbias410.py`), a shared component worth +0.006 pds on H1; ratio mode keeps it.
  (It is not the reference's Jensen bias: that is only −0.0003…−0.0013, `jensen410.py`.)
- **The jackknife credit is the binding constraint.** With K up to 16 ballast cells the carriers' leave-one-out share
  falls by (S_c/S)²: K=16 blocks read C_hat/C 0.19–0.23, 20/85 targets below 1, and the lost credit (+0.0097 raw mse)
  swamped the direction gain the welded arm realised (−0.049). Fix: pair-wise anti-correlated ballast composition noise
  (pair sums exact → pseudobulk and pair contribution to arithmetic means unchanged) sized to the carriers' own jackknife.
- **The key is a BULK key, not a DE key.** Scaling the carriers per target by the key (median 1, clip 0.5–1.8) with the
  bulk held fixed costs nmae +0.013.

## 6. Stacking with 228's per-GENE gain (PLOLO, effect level, `stack410.py`; relative aggregate cosine vs the whitened base)
| fold | gene gain (gain_reg, t 0.25) | key clip(0.3,3) | key uncapped | gain + key clip | gain + key uncapped |
|---|---:|---:|---:|---:|---:|
| k562_gwps | +0.120 | +0.299 | +0.386 | +0.455 | +0.516 |
| orion_hct116 | +0.150 | +0.267 | +0.341 | +0.458 | +0.531 |
| orion_hek293t | +0.142 | +0.287 | +0.350 | +0.459 | +0.524 |
The two stack ~multiplicatively. 228 measured the gene gain's read-back retention at ~0.5 (+8.9 % → +4.55 %), and its
pds gain is concentrated in three targets; its nmae cost (−0.0104 overall) was a weld cost that ballast removes.

**Mass sink (H1):** delivering exactly w x the carriers' read-back effect (target gene absorbs the mass) moves H1 pds to
−0.0069 vs ratio mode's +0.0056; welded scaling read +0.0046 on H1 but −0.0038 raw on the BOARD. H1 and the board disagree
in sign on pds for this arm class, so neither the sink nor H1's pds reading is actionable; parked.

## 7. More H1 DE-channel negatives (ballast pins the bulk, only the carriers move)
| carriers | Δ nmae | Δ fid | Δ reach | note |
|---|---:|---:|---:|---|
| unwhitened donor mean (norm-matched) | +0.0028 | +0.0011 | −0.0221 | whitened carriers stay (388 agrees) |
| ×1.25 | +0.0025 | +0.0001 | −0.0105 | |
| ×1.5 | +0.0123 | +0.0009 | −0.0137 | |
| ×2.0 | +0.0287 | +0.0021 | −0.0107 | DE amplitude is not the lever |

## 8. A better key: 415's production leave-one-donor-out key (worker, Opus 5.5 xhigh; KA rows re-checked)
`key_p = mean_j <P^{-j}_p, d_j> / |P_p|^2` on the SHIPPED composition (whitened base + CD4 + H1 residuals) beats 410's
key on all three cancer folds (+1.2…+5.4 % relative on top of 410's weighting at clip 3, +4…+10 % uncapped); the whole
gain is the residual term <R, d_j>. Artifacts `experiments/415-lodo-key/work/pert_scale_415pp_lo0.3_hi3.npz` (5f1ae21f)
and `_hi1e+09` (5a12079f); rho with 410's key 0.943 on the panel (verified). Queued as the next key swap after 410c/410d.

(board results appended below as they land)

## §9. 410d artifact, measured (2026-09-22 21:05Z): the withdrawal was right, and by far more than 414 predicted

414's `cliffcheck414.py` on the built 410d artifact (`build_d/the-lighthouse-keeper-went-remote.h5ad`) gives
render-space energy vs 410b of **A 1.208 / B 1.192 / C 2.159**. Context C more than doubled. That puts rho at
1.25 / 1.15 / 1.98 against 410b's 1.03 / 0.96 / 0.92, and the amplitude-only raw mse change vs 410b is **+0.500**:
about −0.09 overall, before any direction gain.

The effect-space energy match hid the render-space blow-up in one context. The builder's ballast diagnostics agree:
clip mass up to 1,016,304 counts per block, trim up to 952,304, and jk_ratio max 3.11.

**Lesson for 410d2 (render-gated, summed over contexts):** a single normalisation constant over the three contexts cannot
stop one context absorbing the energy. Any gene-gain arm must match render-space energy PER CONTEXT and pass
`cliffcheck414.py` before upload. 410d2 is shelved until it does.

**Provenance note (added 21:25Z, from 416's known-answer audit).** The `bal_ident_ratio` H1 row in §2's table (line ~55)
was produced BEFORE the ballast-cell cap fix: its ballast cells held 999,000 counts, with no 995,000 trim. The committed
`run410b.py` trims to 995,000 (`ballast_trim` 4000 on every target) and does not reproduce that row bit for bit.
- The rows 416 re-derived with the committed code are the harness's known answers: `base_s0 --bulk ps` equals
  `bal_lo03_hi3_ratio_jk` bit for bit, as does `x125`/`pre`.
- The ratio-mode conclusion (the bulk is restored at w = 1) is unaffected.

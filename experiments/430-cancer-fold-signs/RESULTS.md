# 430: "unwhitened donor-mean signs beat our whitened composition" does not hold for A/B/C. On H1 it comes from the essential-gene donors; on the literal PLOLO folds it comes from same-study siblings. Cross-study, on cancer lines, the donor mean is +0.001 (CI spans 0) and a blend +0.002. No carrier-sign slot

Status: **COMPLETE** (2026-09-23, worker, Opus 5.5).
- Analysis only: PLOLO folds through 418's composition machinery, plus an analytic H1 bridge on stored files.
- No vcc, no cell-eval2, no builder or `src/` edit, no H1 emission, no git.
- Nothing was written outside `experiments/430-cancer-fold-signs/` (checked with `find -newer` over every imported experiment dir and `data/donors`).
- Peak RSS: 5.98–7.15 GB per fold process, 3.3 GB for the H1 bridge, all serial. `memory_pressure` was 70–75 % free before every run, and free+inactive was 28.8 GiB at the start.

## 0. The answer

1. **On the literal PLOLO cancer folds the answer is yes.**
   - Built without L, the raw equal-weight donor mean (M3L) beats the production composition (PP) on the H1-validated gate proxy by **+0.0111 [+0.0068, +0.0154], 3/3 folds**.
   - The mandate's blends 0.5·PP + 0.5·w·M3L gain **+0.0065 / +0.0076 / +0.0084** (w 1 / 1.5 / 2), 3/3.
   - γ 0 (PPg0) gains **+0.0111 [+0.0074, +0.0146]**.
2. **The gain rides on the same-study sibling.**
   - On the Orion folds, the donor pool contains the other Orion line. Per fold, M3L reads +0.0134 (hct116) and +0.0193 (hek293t), but **+0.0006 on the K562 fold**, whose two sources are both cross-study.
   - Built only from donors of a different study than L, which is the board's situation (no donor shares A/B/C's Arc 10x Flex study):

     | arm | Δ vs whitened composition | per fold |
     |---|---|---|
     | raw donor mean | **+0.0012 [−0.0007, +0.0031]** | +0.0006 / −0.0004 / +0.0035 |
     | blend w 1.5 | **+0.0018 [+0.0006, +0.0029]** | |
     | γ 0 | **+0.0042 [+0.0026, +0.0058]**, 3/3 | |

   - kolf2 (non-cancer, all sources cross-study) reads every unwhitened candidate **below** PP: M3L −0.013 … −0.029, blends −0.005 … −0.015.
3. **H1's +0.017 (427e) cannot transfer to A/B/C.**
   - It comes from K562-essential, RPE1 and Nadig-Jurkat. They cover 18–21 of the 85 H1 targets, and those targets hold **60 % of H1's real-gate cells**. The same three donors cover **0 of the 300 A/B/C targets** (the panel has no essential genes).
   - With only the donors that cover the panel (K562-GWPS and the two Orion lines), the same H1 comparison reads **−0.0017** (0.5589 vs our delivered 0.5606).
   - On H1's 38 targets that no essential-screen donor covers, the donor mean scores 0.5454 against our delivered 0.5486.
   - Re-running 427e's own nmae arithmetic with the panel-available donors on those 38 targets:

     | | 427e as published (six donors, all targets) | panel-available donors, 38 targets |
     |---|---:|---:|
     | plain MDE | −0.0098 | −0.0111 |
     | blend + MDE | −0.0245 | −0.0125 |
     | **blend's increment over the plain MDE** | **−0.0147** | **−0.0014** |

4. **Relative to H1's +0.017, in the same proxy units** (z1000 gate, raw truth):
   - Cancer cross-study, donor mean vs whitened composition: **+0.0012**, about **0.07×** the H1 headline.
   - Against the delivered-sign emulation (PPprior), the literal fold numbers are +0.0055 (donor mean) and +0.0020 (blend w 1.5).
   - Under the same proxy on H1, the six-donor mean vs delivered is +0.0058; the panel-available mean vs delivered is −0.0007.
5. **Board prediction for a carrier-sign arm (donor mean or blend into the carriers): about 0 overall, sign +/−. Not worth a slot.**
   - Sign members: the blend buys 0.43× of γ 0's cross-study sign gain. 417c's γ 0 bought fid +0.0012 and reach +0.0024 raw on the board (+0.0011 overall), so the blend's share is about **+0.0005 overall**, a tenth of the ordering-channel sd (0.0045).
   - nmae is unpriced and carries 417c's downside (−0.0029 overall). H1 and this screen's L1 nmae proxy both mis-sign 417c's nmae (below), so neither can price it.
   - The panel-transferable part of 427e's idea is the MDE itself (426a), already armed.

## 1. Known answers (every script aborts on failure; all passed before any candidate number was read)

| row | expected | measured |
|---|---|---|
| KA-PP: `compose_arm`(defaults) vs `sens415.compose_pp` | bit-identical | max \|d\| = 0, s identical, 4/4 folds |
| KA-PP: effect-level aggregate cosine of PP with PLOLO truth (live scored rows) vs `sens415.json` `base_agg` | 1e-12 | **\|d\| = 0**: k562 0.039933927090, hct116 0.030822698203, hek293t 0.025167292995, kolf2 0.015285134638 |
| KA-DLV: PP delivered through the a23 read-back, per-target P·T and P·P in frames A/B/C (the per-target cosine with L's row) vs 418's `PP\|none__x/__a` and **420's `PP\|K3__x/__a`** | 1e-12 relative | **\|d\| = 0**, 4/4 folds. Delivered per-target cosine medians (A/B/C): k562 0.0243/0.0239/0.0235; hct116 0.0188/0.0208/0.0195; hek293t 0.0176/0.0184/0.0174; kolf2 0.0052/0.0043/0.0055 |
| KA-TRUTH: my raw δ_L − column mean vs PLOLO's centred truth, scored columns | float32 | max \|d\| ≤ 1.2e-7 (truth max 1.2–3.4) |
| H1: 427e's gate sign accuracy, delivered / donor average (covered) | 0.5606 / 0.5776 | **0.5606 / 0.5776**, coverage 0.995 |
| H1: 427e's nmae table (donor × 1 / × 1.5 alone, blend × 1.5 alone / + MDE, ours + MDE) | −0.0007 / −0.0090 / −0.0120 / −0.0245 / −0.0098 | −0.00069 / −0.00901 / −0.01197 / −0.02447 / −0.00982 |
| **Board retrodiction (the γ 0 family has a board point, 417c)** | board fid at 0.47× H1's prediction, nmae worse (+0.0089 raw) | **Signs: pass for the cross-study view.** Cross-study γ 0 is +0.0042 against H1's +0.0083 under the same gate, a ratio of **0.51**; literal PLOLO is 1.34×. **nmae proxy: FAIL.** The L1 proxy says γ 0 improves nmae (−0.0037 cross-study, −0.0070 literal), but the board said worse. Read nothing about board nmae from this screen. |

## 2. Instruments

**Folds.** PLOLO's own folds (`plolo.FOLDS`):
- Held-out L ∈ {k562_gwps, orion_hct116, orion_hek293t}, with kolf2 as the non-cancer reference.
- Each fold uses its held operator (kolf2: production) and its sources: the other two production donors (all three for kolf2).
- Targets are the panel targets L covers.
- Counted targets are those covered by the fold's sources and by the cross-study donor set: 264 / 260 / 263 (kolf2 273).

**Truth.** L's own donor δ row (`data/donors/<L>/delta.npy`), pert minus control, in log1p-CP50k units.
- Primary is RAW, which is what a DE test against controls reads.
- Secondary is CENTRED on L's all-perturbation mean (PLOLO's truth).
- kolf2's table is a z-score matrix whose baseline is not control-relative, so the CLEAN4 sign scores 0.23 there. For kolf2, read the centred truth only.

**Gate universe.** PLOLO's scored columns (L measures the gene; panel columns out, which also removes the target itself), and L's control CPM > 5.

**Gate proxies.** The mandate's top-K |δ_L| (K 50 / 200 / 1000) and |δ_L| ≥ τ (0.1 / 0.2), plus **z200 / z1000**: top-K by |δ_L|·√(control CPM), a significance proxy.
- Validated on H1 (`h1bridge430.py`). H1 is the one surface with both the real reference DE gate and a pseudobulk δ from the same cells.
- Rank agreement of 24 sign candidates between each proxy and the real gate (Spearman):

  | proxy | top50 | top200 | top1000 | τ 0.1 | τ 0.2 | z200 | **z1000** |
  |---|---:|---:|---:|---:|---:|---:|---:|
  | Spearman | 0.22 | 0.52 | 0.73 | 0.74 | 0.66 | 0.82 | **0.94** |

- The |δ| gates miss the real gate's many small-|δ|, highly expressed genes, and so mis-rank per-gene sign priors. **z1000 is the primary gate below.**
- The top-K and z gates have a fixed size per target, so their pooled and per-target-mean accuracies coincide. The τ gates report both.

**Metrics.**
- Sign accuracy: a zero prediction scores ½.
- L1 nmae proxy: mean over targets of Σ|k·l̂ − δ| / Σ|δ| on the gate, at the candidate's own L1-optimal global k (exact: the weighted median of δ/l̂ with weights |l̂|/S_p).
- Paired bootstrap over targets, 2,000 draws. The cancer mean resamples each fold independently and averages the three fold Δs.

**Candidates** (all without L):

| candidate | construction |
|---|---|
| **PP** | The 418/420 production composition: normrel cap −1.5 + PM50 base, whitened γ 0.65, + CD4 (k5) + H1 union + CD4 partial, one preserve, panel zeroed. |
| **PPprior** | PP with the CLEAN4 (LOLO) sign on L's top-2,200 CPM block, strongest 20 % protected. An emulation of the delivered clean4q80 signs; on H1 the same emulation of the carriers reads 0.5655 on the real gate vs our delivered 0.5606 (carriers alone 0.5523). |
| **PPg0** | γ 0, with residual doses re-solved and rows norm-matched to PP (416's g000). |
| Wonly / pre | Whitened / unwhitened donor base, no residual. |
| **M3L** | Mandate (3): equal-weight raw δ mean of `data/donors/*` minus L. K562-essential, RPE1 and Nadig-Jurkat cover 0 of 300 panel targets, so M3L is the raw mean of PP's own two sources. |
| M3S | Minus L's whole study. |
| M3Lc / M3Sc | Centred versions of M3L / M3S. |
| **C4L** / C4S | Mandate (4): CLEAN4 per-gene sign with L's bank / study removed. |
| **BLw** | Mandate (5): 0.5·PP + 0.5·w·M3L, rows scaled to PP's norm first. |
| PPxs, PPxsg0, BXw | Orion folds only: PP, γ 0 and the blend built from the K562-GWPS source alone (cross-study). |

## 3. The table: cancer-fold mean, raw truth, Δ vs PP (95 % CI), per fold k562 / hct116 / hek293t

### Primary gate z1000 (H1-validated)

| candidate | sign acc | Δ [CI] | folds + | per fold | L1 nmae proxy | Δ nmae [CI] | k per fold |
|---|---:|---|---:|---|---:|---|---|
| PP (production) | 0.5208 | 0 | — | — | 0.9854 | 0 | 0.81 / 0.62 / 0.51 |
| PPprior (≈ delivered signs) | 0.5264 | +0.0056 [+0.0030, +0.0082] | 3/3 | +0.0077 / +0.0055 / +0.0036 | 0.9853 | −0.0001 | 0.82 / 0.62 / 0.50 |
| **PPg0 (γ 0)** | 0.5318 | **+0.0111 [+0.0074, +0.0146]** | 3/3 | +0.0038 / +0.0125 / +0.0169 | 0.9783 | −0.0070 [−0.0092, −0.0048] | 0.84 / 0.71 / 0.60 |
| Wonly (whitened base, no residual) | 0.5164 | −0.0044 [−0.0056, −0.0033] | 0/3 | −0.0063 / −0.0045 / −0.0024 | 0.9901 | +0.0047 | |
| pre (unwhitened base, no residual) | 0.5293 | +0.0085 [+0.0044, +0.0125] | 2/3 | −0.0022 / +0.0102 / +0.0175 | 0.9823 | −0.0031 | |
| **M3L (donor mean, literal)** | 0.5319 | **+0.0111 [+0.0068, +0.0154]** | 3/3 | **+0.0006** / +0.0134 / +0.0193 | 0.9799 | −0.0055 [−0.0086, −0.0023] | 0.70 / 0.71 / 0.57 |
| M3Lc (centred) | 0.5297 | +0.0089 [+0.0046, +0.0133] | 2/3 | −0.0024 / +0.0107 / +0.0185 | 0.9816 | −0.0037 | |
| **BLw1** | 0.5273 | +0.0065 [+0.0043, +0.0089] | 3/3 | +0.0017 / +0.0077 / +0.0101 | 0.9808 | −0.0046 | 0.89 / 0.83 / 0.68 |
| **BLw1.5** | 0.5284 | +0.0076 [+0.0050, +0.0104] | 3/3 | +0.0018 / +0.0092 / +0.0119 | 0.9802 | −0.0052 | 0.71 / 0.68 / 0.56 |
| **BLw2** | 0.5292 | +0.0084 [+0.0055, +0.0115] | 3/3 | +0.0021 / +0.0100 / +0.0132 | 0.9798 | −0.0055 | 0.59 / 0.58 / 0.48 |
| M3S (donor mean, other studies only) | 0.5193 | −0.0015 [−0.0034, +0.0005] | 1/3 | +0.0006 / −0.0033 / −0.0018 | 0.9889 | +0.0035 | |
| **C4L (CLEAN4 sign)** | 0.5235 | +0.0027 [−0.0016, +0.0070] | 2/3 | +0.0042 / −0.0010 / +0.0049 | 0.9896 | +0.0042 | 0.04 (±1 units) |

### Every mandate gate, Δ pooled accuracy vs PP (cancer mean; per-target-mean Δ in parentheses for τ)

| candidate | top50 | top200 | top1000 | τ 0.1 | τ 0.2 | z200 | z1000 |
|---|---|---|---|---|---|---|---|
| PP acc | 0.5513 | 0.5311 | 0.5180 | 0.5146 (pt 0.5195) | 0.5224 (pt 0.5518) | 0.5369 | 0.5208 |
| PPprior | +0.0079 | +0.0056 | +0.0038 | +0.0032 (+0.0034) | +0.0074 (+0.0097) | +0.0162 | +0.0056 |
| PPg0 | +0.0100 | +0.0091 | +0.0082 | +0.0092 (+0.0076) | +0.0122 (+0.0113) | +0.0142 | +0.0111 |
| M3L | +0.0081 | +0.0071 | +0.0072 | +0.0083 (+0.0061) | +0.0128 (+0.0080) | +0.0156 | +0.0111 |
| BLw1 | +0.0078 | +0.0068 | +0.0050 | +0.0053 (+0.0044) | +0.0081 (+0.0059) | +0.0095 | +0.0065 |
| BLw1.5 | +0.0089 | +0.0076 | +0.0055 | +0.0061 (+0.0049) | +0.0094 (+0.0082) | +0.0107 | +0.0076 |
| BLw2 | +0.0088 | +0.0078 | +0.0060 | +0.0066 (+0.0052) | +0.0105 (+0.0082) | +0.0113 | +0.0084 |
| C4L | −0.0202 | −0.0014 | +0.0052 | +0.0102 (+0.0034) | +0.0170 (−0.0203) | +0.0034 | +0.0027 |
| L1 nmae proxy Δ, BLw1.5 | −0.0082 | −0.0059 | −0.0037 | −0.0040 | −0.0097 | −0.0078 | −0.0052 |

Every M3L and BLw row carries the K562 fold at +0.0003 … +0.0024 and the Orion folds at 5–20× that. The full tables (both truths, every candidate, per-target CIs, k) are in `work/table430.md` and `work/analyse430.json`.

### Cross-study only (the board-like reading): K562 fold as is; Orion folds with the K562-only base. Δ vs the whitened composition

| gate | γ 0 | raw donor mean | blend w 1 | blend w 1.5 | blend w 2 |
|---|---|---|---|---|---|
| **z1000** | **+0.0042 [+0.0026, +0.0058]** (+0.0038 / +0.0025 / +0.0063) | **+0.0012 [−0.0007, +0.0031]** (+0.0006 / −0.0004 / +0.0035) | +0.0016 [+0.0006, +0.0026] | **+0.0018 [+0.0006, +0.0029]** | +0.0018 [+0.0005, +0.0031] |
| z200 | +0.0059 [+0.0032, +0.0085] | +0.0018 [−0.0015, +0.0050] | +0.0034 | +0.0039 | +0.0038 |
| top200 | +0.0043 [+0.0022, +0.0066] | −0.0012 [−0.0041, +0.0017] | +0.0018 | +0.0017 | +0.0018 |
| τ 0.1 | +0.0034 [+0.0021, +0.0046] | +0.0003 [−0.0012, +0.0019] | +0.0010 | +0.0012 | +0.0013 |
| z1000, centred truth | +0.0054 | −0.0002 | +0.0009 | +0.0010 | +0.0009 |

### kolf2 (non-cancer, centred truth, all sources cross-study): Δ vs PP

| gate | PPg0 | M3L | BLw1 / 1.5 / 2 | M3Lc |
|---|---|---|---|---|
| top50 | −0.0042 [−0.0114, +0.0032] | −0.0286 | −0.0116 / −0.0134 / −0.0148 | −0.0192 |
| top200 | −0.0061 [−0.0107, −0.0016] | −0.0239 | −0.0098 / −0.0109 / −0.0128 | −0.0172 |
| top1000 | −0.0026 [−0.0054, +0.0001] | −0.0128 | −0.0045 / −0.0053 / −0.0063 | −0.0095 |

## 4. H1 bridge (`h1bridge430.py`; stored files only, no emission)

Real gate = 323's `reference_de.parquet` (427e's rule; 58 targets). All values are sign accuracy.

| candidate | real gate, pooled | real, per-target mean | z1000 proxy | Δ vs CAR on the real gate [CI] |
|---|---:|---:|---:|---|
| D, our delivered LFC | 0.5606 | 0.5589 | 0.5352 | |
| CAR, shipped carriers (γ 0.65 + CD4) | 0.5523 | 0.5716 | 0.5281 | 0 |
| CARprior, the PPprior emulation on CAR | 0.5655 | 0.5700 | 0.5378 | |
| CARg0, γ 0 carriers (417c) | 0.5654 | 0.5822 | 0.5364 | +0.0131 [+0.0024, +0.0241] |
| **M6, six-donor mean (427e)** | **0.5772** | 0.5899 | 0.5410 | +0.0249 [+0.0096, +0.0391] |
| **M3P, the three panel-available donors** | **0.5589** | 0.5813 | 0.5345 | +0.0066 [−0.0052, +0.0184] |
| blend 0.5·CAR + 0.5·1.5·M3P | 0.5602 | 0.5804 | 0.5343 | +0.0079 [+0.0006, +0.0172] |
| 427e blend 0.5·D + 0.5·1.5·M6 | 0.5761 | | 0.5397 | |
| same with M3P | 0.5626 | | 0.5348 | |

**Where M6's advantage lives.** The 21 targets covered by K562-essential, RPE1 or Jurkat (20 with a real gate) hold 60.2 % of the real-gate cells.

| stratum | D | CAR | CARg0 | M6 | M3P |
|---|---:|---:|---:|---:|---:|
| 20 essential-covered targets | 0.5685 | 0.5566 | 0.5761 | **0.5982** | 0.5677 |
| 38 other targets (panel-like) | 0.5486 | 0.5459 | 0.5490 | 0.5454 | 0.5454 |

## 5. Mechanism

- **Whitening costs signs; the CD4/H1 residual pays some back.**
  - On the literal folds (z1000), whitening the base loses 0.013 (Wonly −0.0044 vs pre +0.0085).
  - The residual adds +0.0044 on the whitened base.
  - Skipping whitening with the residual kept (γ 0) gives +0.011 literal and +0.004 cross-study.
  - This is the sign half of 417c's board read: fid and reach rose. It is not new, and its nmae cost is already measured on the board.
- **The raw donor mean adds nothing over γ 0.**
  - Cross-study it is +0.001 against γ 0's +0.004.
  - Literally they tie at +0.011, both carried by the Orion sibling line, which shares L's study-level structure.
  - Whitening, fitted on the donors' shared signal covariance, damps exactly that shared direction. That is why un-whitening looks strong when a sibling is in the pool.
- **Centring is a small effect.** Raw beats centred by +0.002 on raw truth; the order reverses on centred truth.
- **The CLEAN4 per-gene sign alone** reaches PP's accuracy on the significance gate (+0.0027, CI spans 0) and beats it on τ 0.1 pooled (+0.010), but loses per target on τ 0.2 (−0.020). The delivered prior (PPprior) is +0.0056, 3/3.

## 6. What this predicts for a carrier-sign arm on the board

- **"Unwhitened donor mean into the carriers" (alone or blended):** Δoverall ≈ **+0.0005 from fid/reach**, with nmae unknown and a downside the size of 417c's (up to −0.003).
  - The fid/reach figure scales 417c's measured +0.0011 by the cross-study sign ratio, blend/γ 0 = 0.0018/0.0042.
  - **Expected sign: ≈ 0, possibly negative. Do not spend a slot.**
- **427e's "blend + MDE":** its panel-transferable increment over 426a's MDE is ≈ −0.0014 nmae on H1's panel-like stratum (analytic, 427e arithmetic), not −0.0147. Delivering it is a carrier-design change, the class where H1 mis-signs nmae. **Shelve it.** 426a carries the value.
- **If a carrier sign lever is ever reopened,** γ 0 is the one with cross-study support (+0.004, 3/3), and it has already been read on the board: signs up, nmae down, net −0.0017. The open question is magnitude structure, not signs, which agrees with 422.
- **Instrument lessons for future sign screens:**
  - Price cross-study (drop same-study siblings), not literal PLOLO. The literal fold over-read 417c's sign transfer 1.34× vs the board's 0.47×; cross-study read 0.51×.
  - Use a significance-weighted gate (|δ|·√expression), not top-|δ|.
  - Never price nmae from a pseudobulk L1 proxy for carrier changes. It mis-signed 417c, as H1 did.

## Files

| file | what it does |
|---|---|
| `donors430.py` | Donor panel rows, column means and control profiles on the official axis → `work/donor430-*.npz`. |
| `fold430.py <fold>` | Known answers, candidates and per-target metrics → `work/fold430-<fold>.{json,npz,log}`. |
| `analyse430.py` | Bootstrap, cancer mean and cross-study view → `work/analyse430.json`, `work/table430.md`. |
| `h1bridge430.py` | Proxy validation, 427e known answers, essential-target split, 427e's nmae table with M3P → `work/h1bridge430.{json,log}`. |

Reproduce:
```
uv run python experiments/430-cancer-fold-signs/donors430.py
for f in k562_gwps orion_hct116 orion_hek293t kolf2; do uv run python experiments/430-cancer-fold-signs/fold430.py $f; done
uv run python experiments/430-cancer-fold-signs/analyse430.py
uv run python experiments/430-cancer-fold-signs/h1bridge430.py
```
Each fold takes about 15 s and peaks at ≤ 7.2 GB.

# 702 — PIE as a bulk-direction source, bulk-only, seen targets (preregistration)

Written 2026-10-14 ~23:10Z by the science session [bf003b]. At this time no PIE prediction has been read: the cloud run
`20261014T225548Z` is still downloading assets. Neither the w nor this arm's package exists, and no 10-15 entry exists.
Slot: 10-15 s2 if READY and preregistered by 18:00Z, else 10-16 s1 (coordinator [400b6e]'s rule).

## The arm (one variable vs ALLSW, `658allsw`, the October recipe on the A/B/C path)
- **702pie** = ALLSW's exact env, with `VCC_BULK_LOAD` swapped from the cap-6 LODO dump (d58750b8…) to
  `make_pie_dump700.py --base <cap-6 dump> --seen work/pie_train_perts.json --w <w>` (no `--carriers`). Then re-levelled as 498's arms
  were: r1 at ALLSW's factors, `ntclevel417`, r2 at factors × EV, with at most one r3, and rho within 0.5 % of 0.9793 in
  A, B and C. The script is `build702.sh`.
- The blend, per target with PIE evidence ("seen"), on non-panel genes: u = base/|base|, v = PIE/|PIE|,
  new = |base| · normalise((1 − w)·u + w·v). Seen targets have a non-empty mean-of-contexts PIE effect row, and the seen
  count is recorded in the dump env. Every row keeps its non-panel norm, so ||bulk||/||carriers|| (the P2g tercile variable),
  the per-target magnitudes (the LODO key) and the energy split are unchanged. The carriers (DE design) are not touched.
- PIE row = mean over A/B/C of PIE's lfc_pred mapped into the scorer's effect space with each context's control pseudobulk
  (log1p(5e4 p_ctrl 2^lfc) − log1p(5e4 p_ctrl), renormalised). A shared row, because VCC_BULK_LOAD holds one matrix. A
  per-context version (672) is a later arm if PIE's A/B/C rows differ materially (median cross-context cosine < 0.9).
- **"Seen" rule (transfers to D/E/F):** target ∈ the perturbations of PIE's pinned train split (PIE_splits@396ab956
  `replogle_xdataset/train.json`, any dataset; 10,967 perturbations). It is a function of (panel, split) only, so it is
  computable on the day for the D/E/F panel. On A/B/C, 154 / 300 targets are seen.
- Expected config diff vs ALLSW: `bulk_load`, `bulk_load_sha256`, `bulk_ctx_scale` (+ `de_sign_emitter_sha256` if the pin moved).

## Leakage audit (done before any number; `audit700.py` → `work/audit700.json`)
PIE xdataset's training and validation rows, and its pooled evidence, come only from tahoe / jiang / arc_vcc_25 (train 146 +
val 44 perturbations) / orion. Replogle has weight 0 and no split keys. 0 of h1-test's 85 targets are in VCC25 train/val.
PIE has seen H1's 2025 control pool as an input, and none of h1-test's responses. The K562 GWPS truth (Replogle) is outside
training and evidence. The A/B/C lines (Jurkat, HeLa, CAL-33) are new contexts for PIE (Jurkat appears only as an
evaluation-only Replogle context). Orion HCT116 / HEK293T are in PIE's evidence, so no Orion-based instrument is used here.

## Class
**Transfer** (a better prediction of the responses). AGENTS 11.1c.2 requires a held-out read besides the board read.
- Held-out instruments, both PIE-clean: (1) K562 GWPS rows of the 269 covered panel targets (`eval700_k562.py`; baseline =
  the Orion HCT116/HEK293T donor mean, the plainest donor transfer); (2) h1-test (85 targets; `eval700_h1.py`; baseline =
  435's production H1 key, built without H1 donors).
- Both are read on the SEEN subset as well as on all targets, so the mask's value carries to a different seen share on D/E/F.

## Choosing w (fixed now, before any PIE number)
Grid w ∈ {0.25, 0.5, 0.75} in the blend above (w = 1 would be PIE alone). w* = the grid value with the highest
energy-weighted cosine to the truth on K562's seen targets (the cleaner surface and a cancer line). h1-test's seen subset must
agree in sign: blend(w*) − base > 0 in energy-weighted cosine. If h1-test disagrees, or no w gains on K562 seen, there is no
arm: 702 is recorded as NULL on held-out data and is not built. Tries: 3 values × 2 surfaces, all reported.

## Primary statistic and bins (board, vs ALLSW's own entry)
pds and mse are bit-reproducible, so the bulk read is exact. A bulk change that moves ballast K re-draws pass 1, so the DE
members carry the re-roll sd ≈ 0.00025 overall (537 §2.5).
- **Primary: Δ(pds + mse) in overall units** = (scaled Δpds + scaled Δmse) / 6.
- **TRANSFERS, ≥ +0.0015:** adopt for October, via the same dump route on D/E/F: PIE on D/E/F's controls + the D/E/F seen
  rule, blended into the day-of LODO dump at w*. It needs a def_day pass-through and a rehearsal by 10-21 ([17adab]).
- **WEAK, +0.0005 … +0.0015:** adopt only if both held-out reads gained at w* and the DE members did not fall by more than 2σ
  (0.0005 overall); otherwise do not adopt.
- **NULL, −0.0005 … +0.0005:** not adopted. PIE stays a candidate for the DE side (703: SAMELINE noprog / pt) only.
- **HARM, ≤ −0.0005:** not adopted; record the held-out instrument's over-read.
- Falsifiers (stop and diagnose before reading): DE members move by more than 6σ (0.0015 overall), when a bulk-only arm can
  move them only through the re-roll; or the config diff is not exactly the expected keys.

## Predictions
Filled in by an amendment after the held-out reads and BEFORE the build: the w*, the K562 and h1-test energy-weighted cosine
gains (all targets and seen), and the board prediction. The board prediction uses 518's identity (+1 % energy-weighted
direction X ≈ −0.0041 raw mse ≈ +0.0007 overall) and 417's shape-arm transfer factor k ≈ 0.4–0.65 for pds.

## AMENDMENT 1 (2026-10-14T23:19Z): held-out reads, w*, prediction (written before the dump, the build or any 702 entry)
Reads from the cloud predictions (run 20261014T225548Z; `eval700_k562.py`, `eval700_h1.py` → 700/work/eval700_*.json).
Energy-weighted cosine to the truth (ew), and mini-pds (rank of the right target among all evaluated truths):

| surface (targets) | base | PIE alone | blend seen w 0.25 | w 0.5 | w 0.75 |
|---|---|---|---|---|---|
| K562, seen (135) ew | 0.0442 (Orion donor mean) | 0.0265 | **0.0470** | 0.0449 | 0.0365 |
| K562, all (269) ew / mini-pds | 0.0377 / 0.7370 | 0.0223 / 0.5813 | 0.0391 / 0.7357 | 0.0380 / 0.7209 | 0.0338 / 0.6996 |
| h1-test, seen (44) ew | 0.0971 (435's production key) | 0.1049 | **0.1153** | 0.1255 | 0.1201 |
| h1-test, all (85) ew / mini-pds | 0.0961 / 0.8847 | 0.0882 / 0.6468 | 0.1046 / 0.8810 | 0.1093 / 0.8619 | 0.1068 / 0.8346 |

- **Rule result: w\* = 0.25** (K562 seen maximum). h1-test seen agrees in sign (+0.018). So the arm is GO at w 0.25. Tries: 3 w
  values × 2 surfaces, as frozen.
- Seen share: K562 135/269, h1-test 44/85, A/B/C 154/300.
- Cost seen in the instruments: discrimination falls a little (mini-pds −0.0013 K562, −0.0037 h1-test). PIE alone
  discriminates much worse than our direction (0.58 / 0.65 vs 0.74 / 0.88): its rows share more response across targets.
- **Board prediction (vs ALLSW):** Δmse from 518's identity (+1 % X ≈ +0.0007 overall) at the instruments' X gain (+3.7 % K562,
  +8.8 % h1-test, all targets) with 417's shape-arm factor k ≈ 0.5: **+0.0013 … +0.0031**. Δpds from the mini-pds cost (raw
  −0.0013 … −0.0037, span 0.45): **−0.0005 … −0.0014**. **Primary Δ(pds + mse): +0.0008 … +0.0017**, i.e. the WEAK/TRANSFERS
  boundary. DE members: re-roll only (|Δ| ≲ 0.0005 overall).
- Exploratory, NOT this arm (try 2, found after these reads): centring PIE's rows across targets keeps h1-test's ew gain
  (0.1031) without the pds cost (mini-pds +0.0019 H1, +0.0017 K562), but loses K562's ew gain (0.0377 = base). It is a
  follow-up arm (702c) only after this read.
- 703 (PIE as block sign prior) is WITHDRAWN on its held-out read: on h1-test our call set's non-significant genes (721's
  a0) get 0.5115 from PIE vs 0.5171 from CLEAN4 (700/work/signsplit700.json).

## AMENDMENT 2 (2026-10-14T23:31Z): construction fix found at dump time, before the build or any 702 entry
- **Bug:** PIE's lfc head is supervised only on tested DE genes (≥ 5 CPM). On genes below 5 CPM in the controls its lfc_pred
  is unconstrained. On A/B/C: |lfc| median 8.4, p95 26, on 8,344 of 18,262 genes; on the ≥ 5 CPM genes the median is 0.024 and
  the max < 5. Through p_ctrl·2^lfc these values blew the A/B/C effect rows up: median |row| 20 (vs our 1.8), and the cosine
  of our rows to PIE's was 0.02.
- **Fix (one rule, the scorer's own DE filter):** PIE's lfc is zeroed on genes below 5 CPM in the context's control pseudobulk
  (`make_pie_dump700.py` MIN_CPM, `eval700_h1.py`, `eval700_k562.py`, `explore_center700.py`).
- **Held-out reads re-run with the fix:** K562 is identical (Replogle's genes are all ≥ 5 CPM); h1-test moves in the 4th
  decimal (seen ew 0.1153 at w 0.25, all ew 0.1046 / mini-pds 0.8812). **w\* stays 0.25**; the prediction stands.
  Files: 700/work/eval700_{k562,h1}_v2.json.
- **Arm input (fixed):** `experiments/702-pie-bulk-arm/work/dump_pie702_w025_seen.npz`, sha256 2bc63f5878f786212861498172e540cef35f377b35dde9fe1d0f1e60a42e16dd.
  It is the cap-6 dump d58750b8 with bulk rows blended on the 154 seen targets; carriers, perts, genes and env are byte-copied
  (--keep-env). Max |non-panel norm ratio − 1| is 1e-14; median cos(new, base) is 0.959; median cos(base, PIE) is 0.063.
- PIE's A/B/C effect rows differ by context (median cross-context cosine A-B 0.55, A-C 0.43, B-C 0.45). By this prereg's
  note, a per-context arm (672's VCC_BULK_LOAD_CTX) is the follow-up candidate if this shared arm reads TRANSFERS or WEAK.

## AMENDMENT 3 (2026-10-15T02:35Z): expected config diff gains `bulk_sharpen_pt` (float noise in its audit), before the rebuild
- r1 (02:32Z) stopped at run_arm_generic's exact-diff check: got `bulk_load, bulk_load_sha256, bulk_sharpen_pt`.
- The only moved field of `bulk_sharpen_pt` is the audit's tercile `edges`: [0.22146409043439688, 0.6508407286972743] vs ALLSW's
  [0.22146409043439647, 0.650840728697275], relative 2e-15 / 1e-15. counts (100/100/100), floor_to_low 97, tie_to_mid 0 and q
  are identical.
- Why: the builder reads each target's tercile from the LOADED dump (v = ||bulk|| / ||carriers|| on non-panel genes). The
  blend restores every row's norm to 1e-14 relative, not bit for bit, so the quantile edges move in the last digits.
- **Per-target tiers: 0 of 300 targets change tercile.** The builder's rule was replicated on both dumps. The nearest v to
  either edge is 6e-4 relative away, 10 orders of magnitude above the noise. So the sharpening every target receives is
  ALLSW's. On D/E/F, def_day computes the terciles from the day's dump in the same way, so PIE's norm-preserving blend
  leaves them unchanged there too.
- **EXPECT** is now `bulk_load, bulk_load_sha256, bulk_sharpen_pt` (r1) and the same + `bulk_ctx_scale` (r2). The arm itself
  is unchanged, with the same dump 2bc63f58. The r1 artifact of 02:32Z (never packaged) is moved aside and its h5ad deleted.
- R7a's access had started at 02:32:57Z, because build702's exit opened its wrapper. It was stopped at 02:34Z while T1's first
  fit was still waiting at 551's gate (CD4 banks done, no fit run), and is re-armed behind this rebuild (resumable).

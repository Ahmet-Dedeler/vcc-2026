# 525 — Public DE-channel designs, 2026-09-29

**Answer.** No public design that was opened and tied to a board row produces the mid-pack DE members that beat ours. The cluster around `cellpred-0.9.25` is a shared private pipeline: the string is absent from every repo opened here, from PyPI, and from the Hugging Face model index. The design many teams actually publish and reproduce is AtlasShift (kaipengm2 @ `d24ce4f`). Its bit-identical board copies sit at raw fid 0.5076, nmae 0.9183, reach 0.1493, jac 0.0300. That is worse than our best entry on all four DE members (fid 0.5144, nmae 0.8949, reach 0.2852, jac 0.0408). chill and cqawesome need a sign-accuracy gap of the size FINDINGS §1 (434f, 480) already prices; nothing public supplies that generator.

No scorer was run. No challenge matrix was downloaded. `vcc` was not run.

Board snapshot: `work/leaderboard.json`, 1,211 entries, newest `submission_date` 2026-09-29T00:33:33Z. Raw members below are that file. Descriptions are free text the teams wrote; they are evidence of what was claimed, and only of mechanism where the code was opened.

---

## 1. The cellpred cluster is real, and it is not public

Our row (team Ahmet Dedeler, model `the fume hood took up yoga`, id `be90smZUDcIV7j5MgQ03`, 58 submissions, 2026-09-28T12:38Z): raw pds 0.867428, mse 0.626408, nmae 0.894912, fid 0.514360, reach 0.285247, jac 0.040803. Overall 0.272826. Matches the brief.

Ava Brown, model `cellpred-0.9.25`, id `DoDQu9GyOJqU1Xq3AHgI`, 28 submissions, 2026-09-28T04:29Z, empty description: pds 0.859673, mse 0.761375, nmae 0.880679, fid 0.527837, reach 0.294092, jac 0.033311, overall 0.255540.

Four lookalikes agree with Ava at seed scale and disagree in the fourth decimal. Deltas versus Ava:

| team | model | n | Δfid | Δjac | Δnmae | Δreach | Δpds |
|---|---|---:|---:|---:|---:|---:|---:|
| SQ W | `sq-0929a` | 43 | +2.5e-4 | +2.9e-5 | +9.9e-4 | −6.1e-4 | −3.1e-3 |
| HY H (org Insight) | `hh-v36` | 28 | −3.0e-6 | +1.3e-6 | +5.2e-3 | −2.6e-3 | +6.2e-4 |
| Anna Xie (ゆい) | `M31` | 28 | +1.9e-5 | +4.4e-5 | +1.2e-3 | +1.6e-3 | −6.3e-3 |
| Vivian Xu | `run_2026_28` | 28 | −8.6e-5 | +1.3e-5 | +4.8e-3 | −5.3e-3 | −6.7e-4 |

Shared submission counts (28 for Ava, HY H, Anna Xie, Vivian Xu) and empty descriptions. Not one file: no six-vector is bit-identical, and no jac value is shared at full precision.

Yunfan Jin uses the same model string `cellpred-0.9.25` and is a different configuration (6 submissions, 2026-09-28T16:46Z): versus Ava, fid −0.0128 (0.5150), jac −0.00205 (0.03126), nmae −0.0091 (0.8716), reach −0.0265 (0.2676), mse −0.0248 (0.7366). Same software name, different knobs.

Versus us, Ava is fid +0.0135, reach +0.0088, nmae −0.0142, jac −0.0075, with a weaker bulk (mse 0.761 vs 0.626, pds 0.860 vs 0.867). That pattern is “more calls, acceptable signs, slightly better magnitudes on the reference gate.” Flattening extra genes at ±1% cannot move nmae, because nmae is L1 of log2FC on the reference’s significant genes (metrics brief §5). Call-set expansion alone moves fid up and jac down; it does not buy Ava’s nmae.

**Search for the code.**

- `gh search repos` for `virtual cell challenge 2026`, `vcc2026`, and `virtualcellchallenge` (this session): no repository whose name or description contains `cellpred`.
- `gh search code cellpred-0.9.25`: HTTP 429, twice, after a short retry. The code-search index was not read. Repo trees and blob fetches (below) were.
- PyPI `https://pypi.org/pypi/cellpred/json`: HTTP 404.
- Hugging Face models `?search=cellpred` and `?search=vcc2026`: empty lists.
- Grep of the opened sources (AtlasShift, Matin, mmyen, ferry, danielchen26, biona, forrest, ewissink, Cai review, Rajeeva, heoneyzi, zihanye): no `cellpred`.

Historical “CellPred” papers (a 2010 tissue-composition tool and a 2020 cellular-data paper) are unrelated and were not treated as this model. **Conclusion: `cellpred-0.9.25` is private.** The 0.25 cluster is a shared generator plus small knobs, not a public tutorial.

chill (`Rain All Night`, id `QuwHD9vVMscMbviOiWgw`, 24 submissions, empty description): nmae 0.835077, fid 0.530326, reach 0.307637, jac 0.036774. cqawesome (`m0948pxav`, member Chenxi Qiu, 77 submissions, empty description): nmae 0.827207, fid 0.537906, reach 0.336551, jac 0.038393. FINDINGS §1 (434f) maps those nmae values to about 0.65–0.66 gate sign accuracy under an iid constant-magnitude model, an upper bound if magnitudes follow confidence. 480: each +0.01 of gate-wide sign accuracy is about −0.006 to −0.008 raw nmae, so chill’s 0.835 from our H1 0.9427 needs about +0.08. Public emitters do not publish that sign source.

---

## 2. The public starter is AtlasShift, and it loses the DE channel

Repo [kaipengm2/Virtual-Cell-Challenge-2026](https://github.com/kaipengm2/Virtual-Cell-Challenge-2026), HEAD `d24ce4fdae0cd7cba3ba29546cd8737094eae98a` (2026-09-23, “update readme”). Local copies `work/src/kaipeng/{predict,model,pack,prepare}.py` match the GitHub blob sizes (predict.py 8,885 bytes, blob `81607242…`; model.py 18,688 bytes, blob `b0b048e7…`).

Board author: Kaipeng Zheng, Shanghai Innovation Institute, model `AtlasShift-NoiseMatchedK3`, overall 0.180419 (2026-09-23). Description links the repo. Raw: pds 0.823935, mse 0.894165, nmae 0.901046, fid 0.508885, reach 0.180105, jac 0.030690. README claims public rank 82 and score 0.1545618019 on 2026-09-10. That older snapshot is what the reproductions hit; the author’s current row has moved.

### 2.1 Emission, quoted

`predict.py` `control_template` (lines 130–160) draws `cells * pool_k` real control cells without replacement (`pool_k = 4`, seed `20260910 + context`), sorts them by library size, averages composition within each pool of 4, and sets each predicted cell’s depth to the rounded mean library of its pool. The context mean composition and the depth-weighted bulk composition are computed from all controls.

`predict()` (lines 191–215) builds two desired profiles and calls `dual_moment_counts`:

- `log2fc` space, amplitude 0.6, clip 3, against the control mean composition.
- `bulk_delta` space, amplitude 0.3, clip 3, against the bulk composition.
- Sources `K562_GWPS`, `HCT116`, `HEK293T`, `H1_2025` with weights `[2, 1, 1, 2]`, common-response subtract 1 (`predict.py` 22–28, 181–185).
- CD4 family weight 0.5 (`add_cd4_family`, line 198).
- Promoter-neighbor prior: genes with TSS distance ≤ 5000 bp are capped at a remaining fraction that ramps from 0.15 at ≤500 bp to 1 at 5 kb (`apply_promoter_prior`, `model.py` 240–268; call at `predict.py` 204–206 with fraction 0.15).

`dual_moment_counts` (`model.py` 332–444) fits cell compositions so the arithmetic mean matches the CPM probability and the depth-weighted mean matches the bulk probability (100 iterations, tolerance 2e-4), stochastic-rounds residuals onto the integer row depths, then `repair_bulk_totals` (lines 273–329) moves counts between genes while preserving each cell’s depth. There is no per-target call-set size, no BH threshold, and no gene that is forced significant. Significance is whatever Wilcoxon does to this reduced-variance template. Signs and magnitudes are the fused public log2FC / bulk-delta, clipped, plus the promoter cap.

Fork [mjh72849-cpu/Virtual-Cell-Challenge-2026](https://github.com/mjh72849-cpu/Virtual-Cell-Challenge-2026) HEAD `8543e658d0499e3e158d8481e191f126f70b73be` (2026-09-15, “Update readme”): `model.py` and `predict.py` are the same GitHub blobs as kaipeng (`b0b048e7…`, `81607242…`). Only the README differs.

### 2.2 Board copies

| team | model | overall | fid | nmae | reach | jac | evidence |
|---|---|---:|---:|---:|---:|---:|---|
| Becker Lab | `atlasshift_repro` | 0.154560 | 0.507605 | 0.918270 | 0.149301 | 0.029995 | description: “Reproduction of kaipengm2/Virtual-Cell-Challenge-2026 @ d24ce4f (AtlasShift), unchanged” |
| Yunke Zhang | `statbase-rerun` | 0.154560 | same six raw members, bit-identical to Becker | | | | empty description; vectors equal |
| cccccc | `Version 5-1` | 0.154553 | 0.507601 | 0.918259 | 0.149265 | 0.029993 | description cites source commit `8543e658…` (the mjh fork). Within 1e-5–1e-8 of Becker, not bit-identical |
| gugugaga | `AtlasShift-Public-Repro` | 0.152198 | 0.507264 | 0.918818 | 0.150915 | 0.029992 | description: five-source reproduction, “exact … official integer count emitter” |
| Faye Jiang | `AtlasShift-5source-reproduction` | 0.155706 | 0.507426 | 0.918592 | 0.156412 | 0.030002 | empty description; name only |
| xihun luo | `xihun2` | 0.124736 | 0.505583 | 0.939802 | 0.129332 | 0.030056 | description cites kaipeng commit `8543e658` |
| Wuhan Zhang | `AtlasShift-pool1-20260928` | 0.062707 | 0.317703 | 0.919945 | 0.230652 | 0.028682 | empty description. Versus Becker: fid −0.1899, reach +0.0814, nmae +0.0017, mse +0.0048, jac −0.0013, pds −0.0010 |
| RenormAI | `AtlasShift-DualMoment-Flip-K448-v1` | 0.158845 | 0.507709 | 0.918155 | 0.172390 | 0.029697 | description: “AtlasShift parity with K448 dual-moment rank-based count shaping”. pds identical to Becker |

Wuhan’s model name is the one public variance ablation: dropping the pool-4 template (the name says `pool1`; the submission’s code was not opened) collapses fid by 0.19 and raises reach by 0.08. Pooling is what creates AtlasShift’s calls. Even the pooled row loses to us on fid, nmae, reach, and jac. Later private AtlasShift-named rows (Dalang Yu `AtlasShift C57 C51 Restored` nmae 0.8967 / reach 0.2256; Jessy Liu `AtlasShift v54` nmae 0.8879 / reach 0.2129) improve nmae without reaching our reach, and their code is not the published emitter.

This dual-moment construction is the family FINDINGS history already prices as 147–152 (“moved the order and lost reach”). Do not re-propose it.

---

## 3. Other opened emitters

Blob sizes of the local copies match the GitHub contents API at the SHAs below.

### 3.1 biona001 — calibrated NB, code for the NB not in the repo

[biona001/Virtual-Cell-Challenge-2026](https://github.com/biona001/Virtual-Cell-Challenge-2026) HEAD `cf8dd807540b6477baedf635e60cb79148048591` (2026-09-09). README score table through 2026-09-09. Best public row in that table: b11 overall 0.1318, raw fid 0.498, nmae 0.961, reach 0.220, jac 0.026. Never fid 0.528.

What the repo does contain: `scripts/build_ntc_bootstrap.py` (blob 15,092 bytes). It samples controls with replacement and copies counts unchanged (module docstring, lines 1–7; `rng.integers` at line 178). README maps that null to entry `CGbdzCOTvHYYGIYWQQLA`, overall −0.2987, raw fid 0.003, jac 0.000.

What is prose only, in `results/submissions.md` (README: “Full method notes are kept off-repo”):

- b01A: NB that double-counted library variance → ~436 spurious down calls with ~5 genes injected. Raw fid 0.463.
- b02: depth-calibrated NB (dispersion conditional on library size; ~1 spurious call per 400 null cells) plus full shrunken K562 |z|>1 (~2,500 genes/target). Raw fid 0.472, nmae 0.991. Entry `xxOfLsSjxpyvy2pvSvNT`.
- b05: three sources, ~7,100 genes injected, scorer calls ~560–900. Overall 0.1113, raw fid 0.490.
- b09: |z|>3 collapsed fid via coverage (raw fid 0.376).
- b06: global amplitude ×0.4 worse on every metric.

Board team for the GitHub owner was not matched beyond those entry ids, which the API no longer returns as current rows. The generator that produced b02–b11 was not opened.

### 3.2 ferry-ai — gamma-Poisson control model, measured null

[ferry-ai/virtual-cell-challenge-2026](https://github.com/ferry-ai/virtual-cell-challenge-2026) HEAD `517c49e587204cca748eb00dc49333227616c2ca` (2026-09-28). `src/vcc2026/generator.py` blob 11,784 bytes.

`ControlModel.sample` (lines 208–227): rates from a GMM over PCA state plus log library and kernel regression of nearest controls; genes with phi>0 get a gamma multiplier; then Poisson. Fold change multiplies the rate. Explicitly not a copy of control cells. Dispersion (`_fit_dispersion`, lines 160–205) matches control variance, and for mid-sparsity genes bisects phi so the zero fraction matches.

`docs` checkpoint 0020 (`work/src/ferry/0020-generatore.md` §3.3), context A, 1,500 controls, 4 pseudo-targets × 100 cells, zero effect: real cells 0 significant genes; trial-01 Poisson-on-mean 93 significant/target, 89% up; ControlModel GMM 0.2; KDE 0.0.

Checkpoint 0021: submitted t02, entry `49Gvtu504clN1mIu8T2V`, overall −0.0928, raw fid 0.2534, jac 0.0157, reach 0.1382, nmae 0.9728. First submission `PNn227rxP3bVByS37W41` overall 0.0459, raw fid 0.4580. `--max-calls` (0021 §6) keeps the K genes with largest |log fold change| and sets the rest exactly to zero. Their two-point anchors on panel `vcc2026-val-1`: fid base 0.5123 / replicate 0.8101; nmae 1.0012 / 0.3932; reach 0.0790 / 0.9687; jac 0.0304 / 0.3953; pds 0.5021 / 0.9495. Those are their solve, not a replacement for the spec ranges. They once misread fid as requiring k ≤ n_conf and corrected it in the same note. Current board team for these entry ids: not resolved (the API returns only the current entry).

### 3.3 danielchen26 — two-point columns, shift knob measured flat

[danielchen26/vcc2026-stage1-lab](https://github.com/danielchen26/vcc2026-stage1-lab) HEAD `2cf61dad5a7a744e5943524f1e6054aa194a64e7` (2026-09-23). `src/vcclab/decoder.py` blob 9,531 bytes, saved as `work/src/daniel/decoder.py`.

`design_cells` (lines 72–183): bootstrap 400 real control CPM rows without replacement; optional `lfc_all` sets the mean on every gate gene (mass-centered); for each gene in `r_set`, a two-point column (0, s) whose nonzero fraction is bisected over 24 steps so the mean control quantile hits `0.5 ± shift`; row-renormalise to 1e6; Hamilton or stochastic rounding so counts equal CPM. `force_mean=True` pins the empirical column mean and, by his reading of cell-eval2, causes the scorer to discount the jackknife (he quotes rho 0.484). Docstring lines 126–130: the per-gene shift knob is flat. Scheduling shift by source Wald z (7 slopes, mean shift fixed at 0.10) moved reach only +0.0001 to +0.0009, because the head already hits and purity 0.90 dies at ~0.52 per-gene sign transferability.

E02 (`experiments/E02-…/RESULT.md`, saved): SELENOT, context A, mean locked at 78.18 CPM (lfc +0.322). At fixed mean, raising concentration t past ~1.30 makes the Wilcoxon displacement cross zero and then go negative while lfc stays positive (t=3.00: d −0.131, p 5e-20). Rank sign and mean sign are separable in this two-point family. Measured on one gene.

He claims E[|R|]≈288 in the docstring. FINDINGS 513’s A/B/C median n_real is 31–45. 288 is his claim, possibly a mean of a heavy tail, and was not remeasured here. E08 call-set K sweep is marked 待做 in his notes. `lfc_all` is the same idea as our 478/496 uncalled-gene sign hook (496: nmae −0.00146 raw, weak). Do not propose it again.

Board mapping of this GitHub account: unverified. DanielChow (CUHK, overall 0.1418), Daniel T. (0.0770), and Daniel Murray (0.1867) are different display names and were not linked to this repo by any description.

### 3.4 Matin Gerami — resample controls, multiply by the fold

[MMatinGerami/zero-shot-perturbation](https://github.com/MMatinGerami/zero-shot-perturbation) HEAD `ca4be8145c13571c1ed26e54c973a6d3ba6ccb92` (2026-09-28). `src/zsp/emit.py` blob 1,038 bytes:

```python
def emit_cells(...):
    ratio = np.divide(pred_mean, control_mean, ..., where=control_mean > 0)
    idx = rng.integers(0, control_cells.shape[0], n_cells)
    X = control_cells[idx].toarray() * ratio[None, :]
    # stochastic round
```

With-replacement resample of real controls, per-gene multiplicative fold, stochastic round. No call-set control.

Board: Matin Gerami, `weighted-transfer-6src-v3`, overall 0.0787, raw fid 0.4585, nmae 0.9447, reach 0.2166, jac 0.0250, mse 3.86, pds 0.688. README local cell-eval2 numbers are HepG2/Jurkat holdouts, not the board. The README states that damping magnitude helps mse and hurts every DE metric. That matches biona’s b06 and does not approach the cellpred cluster.

### 3.5 mmyen — identical cells, and a control-resample linreg

[mmyen/virtual-cell-challenge-2026](https://github.com/mmyen/virtual-cell-challenge-2026) HEAD `769d304abe9c9f823ccba761414a76337c9e691a` (2026-09-27).

- `models/baseline/model.py` (3,330 bytes), `predict_profiles`: one rounded mean profile, `vstack` of that same row `cells_per_pert` times (lines 68–73 of the saved copy). 400 identical cells.
- `models/pc_linreg/model.py` (14,292 bytes), `predict`: resample real controls without replacement, multiply by `2**effect` or `1+effect`, round; the target gene is forced to `log2(1 − 0.9)` in log2 space (lines 273–284).

Board: Meilin Yen, `pc linreg model (log 2 fold)`, overall −0.1895, raw fid 0.1944, jac 0.0042, reach 0.0731, nmae 1.0029, pds 0.4974. Identical-cell or near-null emission. Ilyes Baali’s blog (below) reports the same failure mode at about −0.81 for 400 identical cells.

### 3.6 Forrest Sheldon and ewissink — null copies of controls

[forrestsheldon/virtual-cell-challenge-2026](https://github.com/forrestsheldon/virtual-cell-challenge-2026) HEAD `0b4d97488fbd6aeacee8b324adc071f8fe31489f` (2026-09-28). `scripts/baselines/resample_controls.py`: 400 controls balanced over 46 NTC guides (32×9 + 14×8), seed 0, no effect (`sample_rows`, lines 23–33). Report `control-resampling-seed0.md`: entry `Td8NrsH2aoBlrcqxbeDc`, overall −0.303693, scaled fid −1.721. Current board row `Baselines II - a=1+global`: overall 0.0854, raw fid 0.3934, nmae 0.9123, reach 0.2457, jac 0.0320, mse 1.573, pds 0.7574. The transfer builder was not line-audited; the null was.

[ewissink/vcc2026](https://github.com/ewissink/vcc2026) HEAD `0792d8d48b7da3562714353f4340efcb2e6e68ff`. `pipeline/baseline_random_controls.py` (9,350 bytes, matches local): draw 400 real controls, unmodified, then `np.rint`. Board team not mapped. Docstring cites the identical-cell failure as motivation for submitting real cells.

### 3.7 Ilyes Baali’s blog — control-cell template, generator withheld

Fetched 2026-09-29: https://ilyesbaali.me/blog/virtual-cell-challenge/ (HTTP 200, 167,010 bytes, `work/src/ilyes/blog.html`), page date September 1, 2026.

He discloses: predict only the change and apply it to control cells drawn from the target context. 400 identical cells scored about −0.81. Combining six screens and “taking care over how the simulated cells are generated” moved him from rank 167 to 58 then 45. He stops short of the better generator because the competition is running. His scaling formulas are the simplified ones Cai’s review already flags as not `(u−b)/(r−b)` (`CaiZongyuan` `participant-evidence-review.md`, 2026-09-14, opened). Self-reported ranks. No current board team was matched to this post. Cross-type effect cosine ~0.02 and the A≈Jurkat / B≈HeLa / C≈CAL33 guess are his claims; the identity guess matches our 041/356 and is still his inference.

### 3.8 Description-only rows (no code opened)

These are board text. They are not verified mechanisms.

| team | model | overall | fid | nmae | reach | jac | what the description says |
|---|---|---:|---:|---:|---:|---:|---|
| William | `transfer-protocol-v15b` | 0.2323 | 0.5092 | 0.8973 | 0.2601 | 0.0253 | pooled 4-cell columns for the top 20% of genes; DE amplitude 1.2; bulk amplitude A/B/C 0.78/0.78/0.81; sources K562 GWPS, X-Atlas HCT116+HEK293T, CD4. `gh search` for the recipe repo was not completed (code search 429). |
| NiuLai (UConn Health) | `Rank28 Recipe d1a1p1…` | 0.1890 | 0.5105 | 0.9070 | 0.2330 | 0.0258 | “Reconstruction inspired by the public transfer-protocol recipe, not an exact reproduction.” Same 4-cell pooling on the top 20%. |
| Rafael | `H044 … DE-magnitude emitter` | 0.2193 | 0.5126 | 0.8804 | 0.2078 | 0.0314 | per-cell DE magnitude emitter, fold-change-optimal quantile, centred sign, scatter guard; 17+ donors including KOLF2.1J. Closest public-text nmae to Ava, with reach 0.208 far below ours and below Ava. |
| OmicsMan | `R76_FlatCalls` | 0.1536 | 0.5022 | 0.9342 | 0.1661 | 0.0354 | 1,200–2,100 extra mid-expression genes flattened at their own pseudobulk mean ±1% with the source sign; claims pseudobulk unchanged and call count 2,000–3,000. nmae got worse. ±1% cannot explain Ava’s nmae. |
| Bill M (Monad Labs) | `Confident transfer…` | 0.1587 | 0.5111 | 0.9317 | 0.2542 | 0.0303 | deterministic emitter: mass in each gene’s top 5% of cells. mse 5.22. |
| LibraryPrep+FACS | `SynthCell V0.02R01` | 0.1503 | 0.5124 | 0.9389 | 0.2181 | 0.0353 | pooled real-control sampling k=4, magnitude ×2, mean-preserving gamma CV=0.5, dependent rounding. mse 7.07. |
| vcc (AdamWang) | `agree30_b2_v1` | 0.1800 | 0.5111 | 0.9385 | 0.2028 | 0.0311 | pool_k 4, 400 cells. Self-reported H1 probe: reach 0.0553→0.0823, nmae 1.1198→1.0906. |

William’s full recipe, as scored, is worse than us on fid, nmae, reach, and jac. The isolated “variance only on the top 20% of genes” piece is the only part worth a harness test, and Daniel’s shift result plus FINDINGS 513 say reach is sign-limited. It is not one of the three tests below.

---

## 4. Named repos with no DE emitter

| repo | HEAD / pushed | what was opened | board |
|---|---|---|---|
| [CaiZongyuan/virtual-cell-2026](https://github.com/CaiZongyuan/virtual-cell-2026) | HEAD `eb808d9073b326eaad92b5098797fb59cd25ef14`, pushed 2026-09-28 | `participant-evidence-review.md`, `REPOS.md`. Index of upstream clones and a warning that Ilyes’s formulas are simplified. `effect_calibration/predict.py` was listed and not downloaded. Not a frontier emitter. | BioCai `selective_pm_full` overall 0.1494 exists; not linked by evidence |
| [Rajeeva-IITM/vcc2026](https://github.com/Rajeeva-IITM/vcc2026) | pushed 2026-09-19 | README: library-size scale to 5e4, median of the 1,000 lowest-expressing controls, `np.rint`. Placeholder, no scored model | unmapped |
| [heoneyzi/Virtual-Cell-Challenge-2026](https://github.com/heoneyzi/Virtual-Cell-Challenge-2026) | pushed 2026-09-12, archived | CV brief: Jiheon Kang, 6-member team, frozen raw-count stage, PDS 0.5794 on a Jiang24 IFNG/BxPC3 proxy. Explicitly not a leaderboard score. No emitter | unmapped |
| [zihanye-med-sys/virtual-cell-challenge-2026](https://github.com/zihanye-med-sys/virtual-cell-challenge-2026) | pushed 2026-09-27 | README claims the validation controls via Git LFS. Tree pointers were 134 bytes. Not downloaded (challenge matrices stay off this VM). `VC.py` is 0 bytes | unmapped |
| [rikku1983/arc_virtual_cell_challenge_2026](https://github.com/rikku1983/arc_virtual_cell_challenge_2026) | HEAD `dfa584ff4d6d8d3e2012dd0d3e15c6915acb8080`, 2026-09-26, repo size 0 | empty README | none |
| [Liuqi3372/vcc2026](https://github.com/Liuqi3372/vcc2026) | pushed 2026-09-22 | AtlasShift-like layout. `model.py` size 18,627 vs kaipeng 18,688. Not byte-diffed | unmapped |
| [danielchen26/vcc2026-dossier](https://github.com/danielchen26/vcc2026-dossier) | pushed 2026-08-27 | scorer-algebra workbook, older than stage1-lab. Not an emitter | see §3.3 |
| [0x5916/vcc2026-detection](https://github.com/0x5916/vcc2026-detection) and `vcc2026-training` | 2026-09-15 | tree is ONNX, speech, image capture. Name collision. Source not read line by line | none |
| [forrestsheldon/vcc2026-h1-benchmark](https://github.com/forrestsheldon/vcc2026-h1-benchmark) | pushed 2026-09-28 | not opened this session beyond Cai’s citation of commit `d28dd049`. A local H1 harness, not an emission recipe | Forrest’s board row is §3.6 |
| [saberhq/sidechain](https://github.com/saberhq/sidechain) | pushed 2026-09-23 | not opened. 306 already recorded no disclosed score. Board team Sidechain `SER-7abefn` overall 0.1131, description names an emission dial at lambda 0.5 | description only |
| [algaebrown/virtualcellchallenge](https://github.com/algaebrown/virtualcellchallenge) | pushed 2026-08-29 | not opened | unmapped |
| [dung-h/VirtualCellChallenge](https://github.com/dung-h/VirtualCellChallenge) | 2026-02-28 | search description: scGPT baseline. 2025-era. Not opened | — |
| HF `barthazian/vcc2026-broad-development-20260908` | README fetched, 1,789 bytes | 256 files, 5.1 GB of checkpoints and prediction arrays. Card: campaign paused, not an official score, no emission algorithm. Not downloaded | — |
| HF `yraj/vcc_2026_dataset`, `bendidiihab/vcc-2026-data` | API metadata only | yraj README returned HTTP 401 in the earlier fetch. No algorithm read | — |

Prior sweep `experiments/306-recon-oct/RESULTS.md` (2026-09-07): no public repo near the frontier; best described entry then was rank 42 at 0.175. That is still true on 2026-09-29 for anything whose code was opened. The new public material since 306 is AtlasShift’s full emitter (and its reproductions), ferry’s ControlModel, danielchen26’s decoder, Matin’s emit, and mmyen’s two baselines. All of them lose the DE comparison above.

---

## 5. Search limits (unverified, not silent)

- GitHub code search (`cell-eval2`, `vcc prep`, `pert_counts.csv`, `context_A.h5ad`, `cellpred-0.9.25`, `dual_moment_counts`): HTTP 429 this session. Inventory used repo search plus the named trees.
- EuropePMC quoted query `"Virtual Cell Challenge" AND FIRST_PDATE:[2026-08-20 TO 2026-09-29]`: HTTP 503. An unquoted query returned HTTP 200 with hitCount 16,582 (parsed as loose terms); the first hit is the organiser paper, Cell 2026, doi 10.1016/j.cell.2026.08.004. No competitor DE-emission preprint was opened. arXiv API returned HTTP 301 and zero parsed entries.
- Zhihu, CSDN, WeChat, Bilibili, and Japanese pages were not fetched into `work/`. 306 and the earlier pass of this experiment found announcement reposts, not emitters. Treat those surfaces as unchecked.
- Kaggle: no notebook was opened.
- `MaurOrtizz/AtlasShift` appeared in an earlier search and was not opened. A LinkedIn “AtlasShift” job name also exists. Do not conflate them with kaipengm2.

---

## 6. Three h1-test designs we are not already running

Ours, from FINDINGS §1, so these tests have a baseline to differ from: unwelded per-cell libraries; ballast cells so bulk magnitude moves without an nmae cost; a ~2,200-gene CPM block with the CLEAN4 cross-line sign prior and q80 protection; coarse MDE `VCC_DE_MDE=0.1,0.1,0.3,coarse`; `dn_k250_zs`; zero-gene repair `VCC_DE_ZERO_REPAIR=-5`. Sign sources are at the public-data limit (516 negative). The nmae gap is signs on weak and uncalled genes (480: oracle signs at our magnitudes recover −0.254 raw on H1, 2.6× the oracle-magnitude pot). Local H1 nmae in 480 is 0.9427 raw. That is not the board nmae 0.8949. Keep the units apart.

**Known-answer row for all three.** Record the current production emission’s h1-test six raw members and the artifact hash in the same run, before the switch. This experiment did not rescore and does not invent that digest. Splice-off (the new knob at its identity setting) must match that hash. Direction members of digest-identical arms reproduce to 16 figures on the board. Ordering-channel noise is sd ≈ 0.0045 overall (FINDINGS instruments). On H1, reach deltas under ~0.015 raw are noise for K-class emitters.

### Test 1 — minimal-U call-set expansion

**Change.** On the production emitter, add one pre-ranked tranche of genes that are currently uncalled. For each gene in the tranche, binary-search the Wilcoxon U (Daniel’s `design_cells` bisection, `decoder.py` 157–174 at `2cf61dad`) until the gene just clears BH 0.05 in the prior’s sign, at the smallest |log2FC| that does so. Restore column totals of genes that are already scored, and restore the pseudobulk (ballast or a depth-preserving column repair of the `repair_bulk_totals` kind). Tranche size 0 is the known-answer.

**What should move.** fid raw +0.005 to +0.015 and jac raw −0.003 to −0.008 if the tranche is mostly extra calls, which is the fid-up / jac-down half of Ava’s gap (fid +0.0135, jac −0.0075). nmae stays put if reference-gate |log2FC| is pinned: nmae does not read genes the reference did not call, and a just-significant |log2FC| on a reference-gate gene we currently miss is a different, larger edit. reach −0.01 to +0.005 and seed-noisy, because a new significant gene ranks before a non-significant reference gene. 521 already priced a global K of 1,200–1,600 at about +0.005 overall-equivalent on H1 at f 0.35, and about 0 to +0.005 on the disputed board regime (f~0.5 from the fid fit → ~0; f~0.33 from the jac level → +0.0013 to +0.0049). That is smaller than chill or cqawesome and only part of Ava’s fid.

**Why this is not a closed construction.** 192, 385, and 521 closed per-target K taken from a strength key, and closed shrinking a call set below n_real. 426 and 435 closed raising the MDE dose, because the classes over-fill and calls leak. Minimal-U does not raise |log2FC| on genes already at the MDE. 400 identical cells is the jac-destroying extreme (Ilyes’s −0.81, mmyen’s board jac 0.004), not a sized tranche. OmicsMan’s ±1% flatten is the same slogan with a magnitude too small to clear BH on mid-expression genes; this test uses the smallest U that does clear, and it will show whether that U is large enough to tax nmae.

### Test 2 — rank/mean sign concordance on our emitter

**Change.** On h1-test, count genes where the Wilcoxon sign and the arithmetic-CPM log2FC sign disagree under the incumbent emission. If the rate is 0, the arm is a no-op and must hash-match the known-answer. If it is not, reallocate counts inside those genes at fixed gene total and fixed arithmetic mean until the two signs agree, and stay past α so the call set (jac) does not change.

**What should move.** fid and reach, by about the fraction of calls whose rank sign is currently flipped relative to the truth. A 2% flip is the order of Ava’s +0.014 fid only if those flips are wrong versus truth; the first number the test owes is the flip rate, not a score. nmae, pds, and mse stay pinned because means and bulk totals are fixed.

**Why this is not a closed construction.** 426, 435, and 484 move means and the shrink. 422 changes the intended sign (sign-only carriers). E02 shows the flip exists in Daniel’s two-point family at high concentration (SELENOT, t=3, d −0.131 at lfc +0.32), which resembles an extreme-class coarse MDE, and it was measured on one gene in his family. Daniel’s shift dial (reach +0.0001 to +0.0009) is a different knob: it reorders |d| when signs already agree. That measurement stands. Do not rerun the shift dial.

### Test 3 — null-column dispersion match

**Change.** Genes outside the intended call set are emitted as resampled control columns, or as ferry’s zero-matched gamma-Poisson (`generator.py` 178–227 at `517c49e5`). Intended DE columns and the pseudobulk totals stay the production emitter. An effect=0 arm is part of the run.

**Known-answer inside the test.** Effect 0 must produce ~0 significant genes per target. Ferry measured that: real cells 0, ControlModel 0.0–0.2, Poisson-on-one-mean 93 with 89% up (0020). biona’s prose says a calibrated NB gives ~1 spurious call per 400 null cells, and an uncalibrated NB gave ~436. Experiment 153 (`experiments/153-de-bench/RESULTS.md` §18), at the then-shipped median |effect| 0.015, measured production spurious calls 19% versus 0% for real-cell multiplicative resampling, with `C_hat` 0.656→1.000. That 19% is pre-MDE and pre-coarse. Remeasure it. Do not quote it as the current emitter. Splice-off matches the incumbent digest.

**What should move.** jac up if the production emitter still leaks calls on null columns. fid up if those leaks had the wrong sign, down if they were accidental true positives. This moves jac the opposite way from Ava (Ava’s jac is lower than ours). It does not explain the cluster. It prices whether our current calls on non-design genes are real.

**Why this is not a closed construction.** 344 was an empirical-null sizing key (NTC constructs show zero DE), not a generator zero-rate match. 317 closed a global dispersion phi at 7% of the gate. 153 resampled every gene at pre-MDE magnitudes and found no single phi that is both `C_hat`=1 and low spurious calls. A split that leaves the call columns on the production emitter and replaces only the null columns was not that experiment.

### Explicitly not proposed

- AtlasShift dual-moment matching of CPM and bulk (147–152; public copies lose all four DE members; Wuhan’s pool1 name shows the variance knob moves fid and reach inside a band that is still below us).
- A global phi (317).
- Identical-cell emission (Ilyes ≈ −0.81, self-reported; Meilin Yen board fid 0.194, jac 0.004).
- `lfc_all` / signs on uncalled genes as a new arm (496, about +0.0002 overall-equivalent).
- Daniel’s confidence-to-p shift (he measured it flat).
- William’s transfer-protocol as a package (board DE worse than ours on every member).

## Correction (2026-09-29, the 533 claims audit)
D5 (danielchen26's decoder has no mapped board row: its board members are unmeasured; its own shift knob measured flat) and D7 (the +0.08 gate sign accuracy is from the BOARD nmae gap 446b 0.8960 − chill 0.8351 at 480's rate; from H1's 0.9427 it would be ~+0.13 … +0.18). Source: `experiments/533-science-lead/RESULTS.md` §2.3.

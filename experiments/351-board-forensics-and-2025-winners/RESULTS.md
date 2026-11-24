# 351 — Board forensics: rivals gain nmae+mse with jac untouched; 2025 winners say statistics beat pure ML

Status: **COMPLETE** (Parts A+B+C; no harness runs per mandate — analysis of local snapshots + web only).

## Known-answer reproduction (own code, before any new measurement)

| check | expected (record) | reproduced | match |
|---|---|---:|---|
| 319 six scaled members + overall (`data/board-entries/xuDelEik7zGWMXdNu0R4.json` vs same id in snapshots) | 0.21865973944346428, pds 0.7883139225, mse 0.1810128120, nmae 0.1075145060, fid −0.0157591186, reach 0.2221916662, jac 0.0286846490 | identical | **0.00e+00** on all seven |
| 323 `raw.nmae` (`metrics-cap15.json`, cf. 343) | 0.9508364630325662 | 0.9508364630325662 | exact |
| 323 `raw.reach` | 0.25434726355678233 | …823 (last digit) | 1e-16 (float repr) |
| leader−us gap (`docs/SPRINT-HANDOFF-2026-09-09.md`: 0.0443827472; nmae +0.02282, reach +0.01868, fid +0.00860, pds +0.00490, mse −0.00886, jac −0.00177) | as listed | all six + sum | **≤1e-8** |
| `mean(six scaled) == score_avg` on 319 | — | diff 2.78e-17 | pass |

Dedupe guard: 4,767 distinct ids across 69 snapshots; all 840 ids of the temporally-latest
snapshot (`2026-09-09T2137Z.json`) present; sort is by parsed snapshot datetime, not filename
(`lb_*` vs `2026-*` prefixes would otherwise put Aug-30 last).

## What was measured (inputs, estimator, controls)

- Inputs: `data/leaderboard/*.json` (69 snapshots, 2026-08-21 → 2026-09-09T2137Z), stdlib `json`
  only, serial load. Own entries: `data/board-entries/*.json`.
- Dedupe by entry `id`; per-team best = max `score_avg`; top-30 from latest snapshot only.
- Consecutive pairs per team ordered by `submission_date`; Δ on scaled members; Pearson +
  Spearman (rank-robustness); strata: all pairs (3,823), improvers Δoverall>0 (2,416), top-30
  teams (347), unclamped-mse pairs with both `score_mse`>0.001 (543; top-30 subset 205).
- Outputs: `work/entries.csv` (4,767 rows), `work/top30.csv`, `work/improvements.csv`,
  `work/leader_trajectory.csv`, `work/correlations.json`, `work/jac_above.csv` (header-only).

## Results — Part A

### A1. Top 30 (latest snapshot, 840 entries, 2026-09-09T2137Z). Full table: `work/top30.csv`.

Top-30 medians (scaled): pds **0.760485**, mse **0.185225**, nmae **0.147903**,
fid **−0.008921**, reach **0.20546**, jac **0.004073**, overall 0.215542.
We ("Ahmet Dedeler", 0.21866) are **rank 13**, profile **HLLLHH** vs median
(HIGH pds/reach/jac; LOW mse/nmae/fid).

### A2. The 12 teams above us beat us on nmae+mse — and on almost nothing else

Overall-point contributions (Δmember/6), leader row reproduces the handoff exactly:

| rank | team | Δoverall | pds | mse | nmae | fid | reach | jac |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | ViRtUaL cElL builders | +0.0444 | +0.0049 | −0.0089 | +0.0228 | +0.0086 | +0.0187 | −0.0018 |
| 2 | cqawesome | +0.0299 | +0.0071 | +0.0182 | +0.0095 | −0.0028 | +0.0024 | −0.0045 |
| 3 | nailong king | +0.0112 | −0.0024 | +0.0114 | +0.0051 | +0.0011 | +0.0000 | −0.0041 |
| 4 | brocade | +0.0088 | +0.0020 | +0.0082 | +0.0043 | −0.0003 | −0.0012 | −0.0041 |
| 5 | Dongmei Cui | +0.0083 | −0.0022 | +0.0018 | +0.0083 | +0.0051 | −0.0015 | −0.0033 |
| 6 | TOBE | +0.0081 | +0.0003 | +0.0107 | +0.0013 | +0.0055 | −0.0071 | −0.0027 |
| 7 | Aginglab.com | +0.0079 | −0.0077 | +0.0067 | +0.0115 | +0.0037 | −0.0008 | −0.0054 |
| 8 | ZHIRUN HAN | +0.0052 | −0.0021 | +0.0133 | +0.0022 | −0.0047 | +0.0000 | −0.0036 |
| 9 | トルイ | +0.0046 | +0.0001 | +0.0131 | −0.0018 | −0.0025 | +0.0001 | −0.0044 |
| 10 | zhongni ye | +0.0041 | −0.0024 | +0.0042 | +0.0047 | +0.0014 | +0.0003 | −0.0041 |
| 11 | tom yoy | +0.0026 | −0.0028 | −0.0026 | +0.0110 | +0.0011 | +0.0006 | −0.0047 |
| 12 | Hengcui Ge | +0.0004 | −0.0046 | −0.0017 | +0.0105 | +0.0051 | −0.0059 | −0.0029 |

Count of the 12 beating us per member: **mse 10/12, nmae 11/12**, pds 5/12, fid 8/12,
reach 1/12 materially (leader +0.0187; next-best +0.0024), **jac 0/12**.
Seven of twelve trail us on pds. The gap is a levels problem on nmae+mse, not a ranking problem.

### A3. No amplitude trade on the board — nmae gains arrive with mse flat/up

| stratum | corr(Δnmae,Δmse) Pearson / Spearman | corr(Δreach,Δpds) | corr(Δnmae,Δreach) | n |
|---|---|---|---|---:|
| all consecutive pairs | **+0.010 / +0.052** | +0.582 / +0.488 | −0.149 / +0.073 | 3,823 |
| improvers (Δoverall>0) | −0.043 / −0.017 | +0.518 / +0.427 | −0.216 / −0.023 | 2,416 |
| top-30 teams | +0.086 / +0.032 | +0.301 / +0.177 | +0.116 / +0.294 | 347 |

Caveat: 84% of all entries sit at the mse clamp (0.0), so Δmse=0 trivially flatters the
flat/up share. Honest stratum (both endpoints `score_mse`>0.001): corr **−0.061** (n=543);
pairs with Δnmae>0 keep/gain mse in **70%** (205/293); top-30 unclamped: corr **+0.187**,
**65%** (75/116). 317's 28–41:1 amplitude trade is OURS (our operating point), not the
field's: rivals routinely gain nmae without paying mse. reach and pds move together;
nmae vs reach Pearson is outlier-driven (Spearman ≈ 0 to +0.29).

### A4. Leader trajectory (11 entries, `work/leader_trajectory.csv`)

test1 (broken, −0.966) → test2 (0.138) → test4 (0.182) → seven `PromoterAI_precomputed`
entries: flat ~0.17, then Sep-6 wild excursions (fid −0.77/−0.39/−0.86 with jac up to
+0.055 at overall ≤0.114), then the winning break **Sep-7 17:06Z: +0.1267 to 0.2402,
carried by a fid repair (−0.855→−0.002, +0.1422 overall-points)** plus reach +0.0056,
nmae +0.0043, pds +0.0023, with mse −0.0198 and jac −0.0079 paying; **Sep-8: +0.0219 to
0.2621 via mse +0.0072, fid +0.0066, nmae +0.0052**; Sep-9: +0.0010 flat. Current:
pds 0.8177, mse 0.1279, nmae 0.2445, fid +0.0358, reach 0.3343, jac 0.0181.
Lesson: the leader won with a fid repair followed by an mse+nmae step — both members we
trail on — while never touching our jac lead.

### A5. Jaccard: nobody is above us

Zero teams hold a best entry with `score_jac` above ours (0.0286846490); `work/jac_above.csv`
is header-only. The only entries ever above 0.029 are the leader's own Sep-6/7 excursions
(jac +0.0554/+0.0478) at overall 0.114/0.051 with fid ≤ −0.77 — a high-jac/low-fid regime,
not a competitive one. Our jac lead is intact; no defense needed.

### A6. Where big steps come from (competitive zone, `prev score_avg ≥ 0.15`, 187 steps ≥ +0.005)

Mean overall-point carry: **mse +0.0069, pds +0.0033, reach +0.0016, nmae +0.0015**,
fid +0.0005, jac +0.0002. At prev ≥ 0.20 (17 steps): pds +0.0049, mse +0.0030 lead; biggest
single elite nmae step is TOBE's +0.1019 scaled (+0.0170 overall-points, 2026-09-08).

## Results — Part B (web, accessed 2026-09-10; quotes ≤ 15 words)

2025 task in one line: few-shot H1 prediction (150 training + 50 validation + 100 test
perturbations), three metrics — PDS (L1 rank), DES (DE-set overlap), MAE — with MAE
"effectively excluded from competitive optimization" (wrap-up §4) because almost all models
scored below baseline on it. 2025 page frozen leaderboard confirms: 3rd-place TransPert won
with mae **4.219** against 1st/2nd at 1.026/1.069
(`https://virtualcellchallenge.org/2025`).

1. **1st — BM_xTVC (BioMap Research; Yucheng Guo, Qirong Yang), model xTrimoSCPerturb.**
   Class: hybrid — "new single-cell pretrained model using an improved scFoundation
   architecture, then integrated it with protein model embeddings" (wrap-up §5). Training:
   "selected a subset of the public perturbation datasets listed on the challenge website".
   Cells/counts: "trained on pseudo-bulk RNA-seq data, repeating outputs multiple times to
   satisfy DES evaluation requirements". DE handling: the key line — "Purely AI-based
   approaches did not consistently outperform statistical baselines" on DES/MAE, so they
   "incorporated training-set DEG frequency and mean expression levels as explicit features";
   fused loss "PDS + DES + small-weight MAE". Note: the only public code repo under the
   winner's name (`github.com/sxiweiLi/BM_xTVC_Solution_VCC`) is a third-party
   ESM-2→expression MLP reimplementation, not the official solution; official details are
   "forthcoming academic paper" only.
2. **2nd — XLearning Lab (Xi Peng, Sichuan University), model X.** Class: explicitly
   metric-driven conditional generation. "PDS carried roughly twice the weight of DES in
   final scoring"; "shifted from noisy single-cell to pseudo-bulk data representation";
   "Fully Connected Network using aggregated control expression, ESM-2 protein embeddings
   for perturbed genes, and UMI count indicators"; "Residual learning predicted perturbation
   deltas"; training "exclusively used publicly available Perturb-seq datasets, including a
   small H1 dataset available as part of PerturbAtlas". Group's public org
   (`github.com/XLearning-SCU`) hosts PerturbFormer, a "virtual cell foundation model".
3. **3rd — Outlier (Chicago/Dartmouth/HKU), model TransPert.** Class: "statistical framework
   for cross-cell-line perturbation prediction using only summary-level data, namely
   pseudo-bulk perturbation profiles and the differential testing results from the Wilcoxon
   test". Mechanism: "combines gene-level perturbation summaries from multiple reference
   cell lines using similarity-aware aggregation, then refines these estimates by predicting
   which genes are likely to be significantly affected in the target cell line"; "global
   linear scaling optimizes the predictions for PDS", scaling via CV. Method paper:
   arXiv:2511.16954 (ICLR 2026 spotlight) — "PDS is highly sensitive to the choice of
   similarity or distance measure and to the scale"; L1/L2 PDS "behave very differently
   from cosine-based measures, even after norm matching"; scaling makes rankings "behave
   increasingly like cosine similarity".
   Generalist Prize — Cleopatra (Altos Labs, `go-with-the-flow`): "flow matching generative
   model", "customized U-Net architecture", "pretrained on approximately 7 million
   high-quality single cells" incl. internal screens, then fine-tuned on h1ESC.
4. **2026 cycle:** no public write-ups or forum posts found (Discord not accessible; search
   2026-09-10). Live intelligence is the board itself: leader `ViRtUaL cElL builders` is
   org **Illumina** (members incl. Eyal Ben David, Kyle Farh's group area) with model name
   **`PromoterAI_precomputed`** — matching Illumina's public PromoterAI ("deep learning
   model for predicting the impact of promoter variants on gene expression", Science 2025;
   `github.com/Illumina/PromoterAI`; precomputed scores for "all human promoter SNVs").
   Whether literal or misdirection is unverifiable; Illumina's July-2026 press release also
   documents a proprietary **Billion Cell Atlas** ("350 million cells sequenced", "200
   disease-relevant cell lines"). Related 2026 literature (not VCC entries): OCOO-T
   (arXiv:2606.12838, minimalist flow-matching Transformer, Tahoe100M/Replogle/PBMC SOTA
   claims), SCALE (arXiv:2603.17380, conditional-transport foundation model), VCBench
   (arXiv:2604.27646: "naive dataset aggregation can even reduce performance"; metrics give
   "substantially different model rankings").

## Results — Part C (≤3 proposals; nearest priors checked by grep of `experiments/*/RESULTS.md`)

**P1 — Recipient-conditioned DE-gate predictor (moves jac/fid/nmae-gate).**
TransPert's transferable piece is "predicting which genes are likely to be significantly
affected in the target cell line" — a supervised gate model, not a response model. Ours is a
fixed rule (top-400 |H| ∪ top-2,200 CPM); 343 shows called genes are sign-right at 0.638 vs
0.552 uncalled, and only 34.6% of gate mass sits in our set. Nearest priors: 192 (call set
worth +0.009 from a gene-level prior + protected |H| head — PRECEDENT this direction pays),
344 (empirical-null key dead — different source: supervised cross-line features, not null
structure), 341 (own-line labels judged by retrieval — different endpoint: gate membership
judged by jac/fid, never retrieval). Different: per-(target,gene) classifier trained on the
six source lines (features: donor DE frequency à la BM_xTVC's winning trick, recipient basal
CPM, donor agreement), applied to the held-out line. Cheap falsifier on 323: swap the call
set at matched n_pred, require Δ(jac+fid) > 0 on paired bootstrap; kill if ≤ 0.

**P2 — Precision-first call-set shrink via internal sign agreement (moves fid).**
Our fid (−0.0158) trails the top-30 median (−0.0089); 8/12 teams above us beat us there while
0/12 beat us on jac — we can spend call-set size for precision. 343 hands us the key our own
outputs already contain: designed-vs-delivered sign agreement splits accuracy 0.649/0.558
(0.638/0.552 delivered). Nearest priors: rc_callset (size cut 4,480→2,599 paid via jac —
different mechanism: agreement-gated membership, not rank truncation), 344 (null key —
different key), 342 (Hamilton rounding — different knob). Different: require designed≡delivered
sign agreement for call-set membership at fixed n_pred (backfill by |H|). Cheap falsifier on
323: Δfid at fixed n_pred; kill if the precision gain does not survive the coverage term
`min(1, n_pred/n_real)`.

**P3 — Promoter-regulatory prior for the sign deficit (moves nmae/reach-sign half).**
343's deficit is 0.588 |lfc|-weighted sign accuracy with up-genes at 0.516; every tested sign
key is closed pool-blind (234) — but no experiment has used *regulatory-sequence* information,
and the 2026 leader's org+model-name pair points at exactly that (PromoterAI, Illumina).
Nearest priors: 333 (GenePT vectors quarantined for collisions — different modality and
failure mode), 334/336 (description text for response prediction — different source), 155
(control covariance direction — different source). Different: per-gene regulatory features
computable locally with no downloads (CpG/promoter conservation/GC from the reference genome
annotation already vendored or in `data/external`; or PromoterAI-style TSS-window scores if
licensed) tested as a STRATIFIER of sign errors first — correlate with per-gene sign
correctness on 323's saved `stats-cap15.npz` + `reference_de.parquet` (pure arithmetic, no
harness run). Falsifier: |corr| ≈ 0 → dead before any model is built. "Nothing new" was
considered and rejected: sequence-derived direction priors are genuinely untested here.

## Verdict against the PLAN.md kill criteria

- A: kill = dedupe/known-answer failure — all four known answers reproduce (three exact,
  one to 1e-16); dedupe guards pass. Not killed.
- B: kill = no public 2025 write-ups — found: official wrap-up with all three methods +
  TransPert paper + repos. Not killed.
- C: kill = everything maps onto CLOSED — P1/P2/P3 each differ in source or endpoint from
  066/155/192/227/232/234/317/318/333/334/336/341/342/344 (stated per proposal). Not killed.

## Limits

- Board forensics is validation-panel (A/B/C) evidence; 343/H1 numbers are development
  evidence. Correlations are observational (shared-team methods, clamp artefacts, junk-entry
  jumps) — never a causal lever or score forecast. Improvement pairs double-count teams with
  many entries; elite-zone strata (prev ≥ 0.15/0.20) are the honest read.
- 2025 lessons transfer partially: 2025 was few-shot single-context H1 with L1-PDS; 2026 is
  zero-shot multi-context with cosine-PDS + six members. TransPert's scaling trick is
  2026-irrelevant (cosine is scale-invariant); its gate-prediction step is the live piece.
- Leader-reading is circumstantial: model names may be misdirection (ours are, per §14).

## Exact commands to reproduce

```bash
memory_pressure | tail -1   # proceed only if free ≥ 50% (was 75%)
python3 experiments/351-board-forensics-and-2025-winners/work/ka_check.py
python3 experiments/351-board-forensics-and-2025-winners/work/build_boards.py
python3 experiments/351-board-forensics-and-2025-winners/work/followup.py
```

## Peak memory and wall time

- `build_boards.py` (heaviest): peak footprint ≈ 270 MiB (`/usr/bin/time -l`), wall ≈ 15 s,
  serial stdlib JSON over 69 snapshots (largest ~1 MB). Total Part A compute < 2 min.
- Part B: web search/fetch only, 2026-09-10. No downloads, no harness/builder/scorer runs,
  no files written outside `experiments/351-board-forensics-and-2025-winners/`.

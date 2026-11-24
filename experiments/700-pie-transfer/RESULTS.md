# 700 — Arc's PIE as a source for our predictions: setup, leakage audit, held-out reads

Science session [bf003b], from 2026-10-14 22:20Z. Status: **in progress.** The sections below fill in as reads land.

## Why
- The field moved while we were idle (07-14 Oct): we fell from #24 to #66 at 0.2731, and 19 teams sit above 0.30.
- Between 10-08 and 10-14 the median top-25 team gained, in scaled units, pds +0.039, mse +0.042, nmae +0.014, fid +0.007,
  reach +0.009 ([400b6e]'s decomposition). That is mostly bulk.
- Arc released **PIE** on 10-02/10-05 (bioRxiv 10.64898/2026.10.02.756297; github ArcInstitute/pie; HF arcinstitute/pie).
  It is a pretrained Perceiver IO that outputs, per (context, perturbation, gene), p_de, lfc_pred (log2 FC of mean CP10K)
  and delta_p.
  - Its training labels are cell-eval2's own DE recipe: Wilcoxon, a 5 CPM filter, BH per perturbation.
  - Inputs: knowledge sources (ESM2, NCBI/gene/perturbation/context text, STRING, DepMap), the context's control mean, and
    per-gene response evidence pooled from training (same perturbation in other contexts; same context, other
    perturbations).
  - The xdataset checkpoint has had 483 HF downloads since 10-05.
  - Team "123 123" cites PIE in its board description and went 0.2524 → 0.2993: pds +0.070, mse +0.100, nmae +0.066 scaled.
  - Board names containing "PIE" with naive emission sit ≤ 0.2467. So PIE's value for us is as a SOURCE inside our stack,
    not as a standalone model.
- Licence: the model is non-commercial, but Outputs are carved out (AUP §1); cite Verma et al. 2026 in the write-up
  ([020377], 680).

## Setup (reproducible)
- **Code:** arc-pie in its own Python 3.12 venv (`experiments/700-pie-transfer/.venv`, `uv venv --python 3.12`;
  `uv pip install arc-pie`).
- **Weights:** `arcinstitute/PIE_replogle_xdataset` `best_auprc.ckpt` (966 MB, sha256 e60a1304…), trained on Tahoe 0.68,
  Orion 0.30, Jiang 0.01, VCC25-H1 0.01.
- **Assets:** pinned HF commits resolved into `data/external/pie/data` by `work/fetch_assets700.py` (~75 GB: ncbi_text
  43 GB, Tahoe preprocessed 22 GB, Orion 6 GB, the rest < 3 GB). Downloaded 10-14 22:26–23:30Z. The download was paused once
  (22:30–22:51Z) for the coordinator's ALLSW upload: big downloads starve uploads, so pause them in announced upload windows.
- **New contexts:**
  - `work/make_ctrl_h5ad700.py` writes log1p(CP10k) control h5ads; `pie prep controls_only=true` runs on the challenge gene
    axis ∩ PIE gene_text (18,262 / 18,533 genes).
  - A → `jurkat` (PIE's stored Replogle Jurkat text), B → `CVCL_0030` (HeLa), C → `CVCL_1108` (CAL-33).
  - The B and C texts were rendered with PIE's own `build_context_entry` from Cellosaurus 56.0 (PIE's pinned release; the
    Jurkat render reproduces the stored text byte for byte). They were embedded with text-embedding-3-large through the
    Hack Club AI proxy: re-embedding the stored Jurkat text gives cos 0.99997 against PIE's vector.
  - The extended source is at `data/external/pie/local_sources/context_text`.
- **D/E/F:** `ccle_context_text700.py` precomputes context texts for every CCLE line (1,749 of 1,981 DepMap IDs mapped to
  Cellosaurus; records PIE's age parser rejects, e.g. ">50Y", are skipped and listed in its log), so the day needs no
  embedding API.
- **The evidence cache does not fit this Mac.** Arc's build allocates 9 float64/int64 (perts × genes) accumulators plus 2
  per-dataset int64 cubes: ~11.3k × ~25k → ~55 GB peak. It was built on one disposable EC2 r7i.8xlarge (`cloud700.sh` /
  `cloud700_worker.sh`, run 20261014T225548Z, instance i-0234538365f437e85, terminate-on-shutdown, dead-man 300 min) with
  Arc's unmodified code, together with PIE inference on three query sets. The artifact is copied to S3 for local reuse on
  D/E/F day.

## Leakage audit (before any number; `audit700.py` → `work/audit700.json`)
- Datasets with train/val rows, and the evidence datasets: tahoe, jiang, arc_vcc_25 (train 146 + val 44 perts), orion.
  Replogle has weight 0 and no split keys.
- h1-test: 0 of its 85 targets are in VCC25 train/val. PIE has seen H1's 2025 control pool as an input (ARC_H1 control
  means), and none of h1-test's responses. Same-context evidence is keyed by (dataset, context), and our h1test dir is its
  own dataset, so none is served.
- K562 GWPS truth (Replogle): outside training and evidence; clean.
- Coverage: PIE has same-perturbation evidence ("seen" = in the pinned train split's perturbations, 10,967) for **154 / 300**
  A/B/C panel targets (148 via Orion, which PIE filtered to ≥ 30 % knockdown; our own Orion donor tables cover ~290). Every
  metric is reported on all targets and on the seen subset.

## Cloud run 20261014T225548Z (r7i.8xlarge, 247 GiB, us-east-1)
- Asset download on EC2: ~2 min for ~75 GB. Evidence build: 6 min, **peak RSS 67.5 GB** (`infer_k562.log`, /usr/bin/time -v).
  It cannot be built on this 64 GB Mac.
- Inference on CPU: K562 269 rows, 3 min; h1test 85 rows, 2 min 50 s, peak 15.7 GB with the evidence cached; abc 900 rows
  (pending).
- Inference abc (900 rows × 18,262 genes): ~16 min on CPU (3.4 cores used).
- Artifacts: `s3://<your-s3-bucket>/output/700-pie/20261014T225548Z/` (out/pred_*.parquet, evidence/ 22.8 GB),
  copied locally to data/external/pie/cache/evidence/87a01957… (aws s3 sync, 10-14 23:33Z, outside the coordinator's upload
  windows).
- **Cost:** instance i-0234538365f437e85 ran 22:56–23:27Z (0.52 h × $2.117 = $1.10); EBS negligible; egress of ~23 GB to the
  Mac ≈ $2.1; S3 storage ≈ $0.53/month for the evidence. ≈ $3.2 one-off, on AWS credits. Termination verified
  (`describe-instances`: terminated); `scripts/aws_reap.sh`: 0 flagged (only the long-running peak-gateway t4g.small).

## Reads (held-out, PIE-clean)
**Bulk direction** (`eval700_k562.py`, `eval700_h1.py`). ew = energy-weighted cosine to the truth's pseudobulk effect, panel
genes removed; mini-pds = rank of the right target among all evaluated truths.

| | K562 all ew / mini-pds | K562 seen ew | h1-test all ew / mini-pds | h1-test seen ew |
|---|---|---|---|---|
| ours (K562: Orion donor mean; H1: 435's production key) | 0.0377 / 0.7370 | 0.0442 | 0.0961 / 0.8847 | 0.0971 |
| PIE alone (lfc → effect space) | 0.0223 / 0.5813 | 0.0265 | 0.0882 / 0.6468 | 0.1049 |
| blend, seen only, w 0.25 (702's arm) | 0.0391 / 0.7357 | **0.0470** | 0.1046 / 0.8810 | 0.1153 |
| blend, seen only, w 0.5 | 0.0380 / 0.7209 | 0.0449 | 0.1093 / 0.8619 | 0.1255 |
| exploratory: centred PIE rows, seen, w 0.25 | 0.0377 / 0.7387 | — | 0.1031 / 0.8866 | — |

PIE alone discriminates targets much worse than our direction: its rows share more response across targets. In a light
blend it adds a little: mse-relevant direction goes up and pds-relevant discrimination goes slightly down. Centring removes
the pds cost and also K562's ew gain.

**A/B/C predictions, and a construction bug found at dump time** (702 AMENDMENT 2): PIE's lfc head is unsupervised below 5
CPM. On A/B/C its |lfc| there is median 8.4, p95 26 (8,344 of 18,262 genes), and that blew the effect rows up. PIE's lfc is
now zeroed below the scorer's 5-CPM filter everywhere. The held-out reads above are unchanged by the fix (K562 identical;
h1-test 4th decimal). PIE's A/B/C effect rows differ by context (median cross-context cosine A-B 0.55, A-C 0.43, B-C 0.45)
and are nearly orthogonal to ours (median cosine 0.063), yet they help in a light blend on both held-out surfaces.

**DE signs: PIE loses to what we ship** (`signsplit700.py`, h1-test, the production emitter's call set = head | CPM block,
adjudicable genes):

| sign source | called & truth-sig (n 25,061) | called & NOT sig = 721's a0 (n 195,858) | head, truth-sig |
|---|---|---|---|
| our emitted lfc | **0.6007** | 0.5154 | **0.7222** |
| CLEAN4 prior | 0.5778 | **0.5171** | 0.5753 |
| our design key | 0.5667 | 0.5076 | 0.7222 |
| PIE lfc | 0.5594 | 0.5115 | 0.5725 (seen 0.648) |
| PIE where p_de ≥ 0.5, else CLEAN4 | 0.5656 | 0.5140 | 0.5692 |

- PIE's real DE strength is RANKING: AUPRC of truth DE genes 0.175 vs our |key| 0.148 and our emitted p 0.142 (random 0.074).
  But its signs are no better than ours.
- n_real: PIE Σp_de Spearman 0.43 vs 521's strength key 0.55–0.60, so per-target sizing stays closed.
- **703 (PIE as block sign prior) WITHDRAWN** on this read; 721's rule needed PIE a0 ≥ 0.54.

**a0 sources tried** (`a0sources700.py` → `work/a0sources700.json`; h1-test, called & not significant, n 195,858). Nothing
beats CLEAN4. This closes the question on public/derivable sources, together with 651 (library coupling) and 610
(CPM coupling, expression range):

| source | a0 | coverage |
|---|---|---|
| CLEAN4 (today) | **0.5171** | 0.82 |
| our emitted lfc | 0.5154 | 1.00 |
| our design/bulk key | 0.5076 | 1.00 |
| panel-mean of our key (shared programme) | 0.5013 | 1.00 |
| PIE lfc | 0.5115 | 1.00 |
| library-size coupling, ORACLE per-target shift (sign β_g·D_t) | 0.5008 | 1.00 |
| logistic, 5-fold target cross-fit, all of the above + log CPM | 0.5170 | 1.00 |
| same without the oracle term | 0.5172 | 1.00 |

The per-target library shift is tiny on h1-test (median D −1.0 %, q10 −6.5 %, q90 +3.0 %). [020377]'s bound from 721: even an
oracle context-mean direction scores 0.552 on these genes. So the pack's ~0.53 is per-target skill on small effects, not a
shared programme or a normalisation artefact.


## Where it stands (2026-10-14 ~23:55Z) and what runs next
- **702 (PIE bulk, w 0.25, seen targets, 5-CPM fix; dump 2bc63f58):** build702.sh is armed (r1 → ntclevel → r2) behind DN's
  build and KAu. Its READY line feeds the coordinator's 10-15 s2 chooser (`READY pie702*` before 18:00Z, else K-UP).
  Reader: `702/read702.py` (KA: ALLSW vs itself = 0). Follow-ups 702c (centred) / 702x (per-context, 672) are preregistered
  as a choice rule in `702/PREREGISTRATION-FOLLOWUPS.md`.
- **703 withdrawn** (PIE signs lose to CLEAN4 on the held-out a0). **704** (α ×1.25) is preregistered and unarmed; it is
  built only if DN reads NULL/HARM. **701** (K-UP) is the s2 fallback; 721's analysis says low value.
- **R7 sealed reserve** (coordinator-approved early access):
  - R7b read: hook pair −0.00230 overall-eq, but R7b is not board-like by the claim's own test (CP fid +0.0077), so no veto;
    the hook stays.
  - R7a is armed behind 702's and DN's chains.
  - Record: `552/RESULTS.md`.
- **D/E/F day path for PIE (insurance until 702 reads):** `dumphook700.sh` for [17adab]'s DEF_DAY_DUMP_HOOK (690g). The
  evidence (22.8 GB) and the CCLE context texts (1,786) are local, so the day needs no network. A mock-day rehearsal on A/B/C
  plus the context-specificity test (`ctxspec700.py`, Replogle 4 lines) are queued behind R7a (`queue_gap700.sh`).
- 721 (news lane's board-as-instrument) replicated independently (`replicate721.py`): +0.0123 vs their +0.0124.

## PIE line CLOSED (2026-10-15T04:50Z, [bf003b])
- **Bulk direction: board HARM.** 702 (w 0.25, seen targets) read 0.2735387, −0.0025464 vs ALLSW (primary −0.0023631: pds
  −0.0042 raw, mse +0.0049 raw worse). The held-out surfaces over-read it. Post-mortem with bootstrap CIs and the gene-axis
  split: `702/RESULTS.md`. 702c and 702x are not built: their preregistered rule needed pds and mse to move in opposite
  directions, and they did not. `dump_hook` is decided **off** for October (99df0055).
- **Signs:** lose to CLEAN4 on the held-out a0 (0.5115 vs 0.5171; 703 withdrawn).
- **Target strength:** Σp_de Spearman 0.42 vs our strength key's 0.55–0.60 (h1-test, seen context). Per-target sizing stays
  closed.
- **DE-gene ranking:** the one property where PIE beats us on held-out (AUPRC 0.175 vs 0.148; unseen-target subset 0.151
  vs 0.137 / 0.146). It is measured only on h1-test, a seen context, and it would need a call-set use (which genes to
  call), never a bulk use. It is not pursued for October: the call set's membership is CallSet(400, K), and 528/CPH closed
  propensity-block changes on the board.
- The queued mock-day and ctxspec jobs were stopped before they started (`work/queue_gap700.nohup`). The hook code
  (`dumphook700.sh`), the CCLE context texts and the cloud predictions stay on disk for the write-up.
- Why PIE helped teams such as "123 123" and not us: their base direction was weaker than PIE's, while ours (the donor
  stack plus LODO and sharpening) is better than PIE's on A/B/C. 702's blend moved our rows toward a worse direction.
  The field's PIE gains are evidence about their bases, not about PIE's ceiling.

- **Open idea from the board (2026-10-15T21:31Z, news lane [730b99]):** team "Chen Zhang" (0.0679 → 0.2791 on 10-15) describes
  "PIE p_de rank on stage-2 padding candidates (single factor; PDS/MSE bit-identical to base)". That uses PIE to choose
  which genes PAD the call set, not to sign them. It matches the one property where PIE beats us on held-out (DE-gene
  ranking, AUPRC 0.175 vs 0.148). For us the analogue is the CallSet's CPM block (85 % of our calls, ranked by control CPM),
  re-ranked by PIE p_de within the gated universe.
  - The closest closed family is 528/CPH's call-propensity block (lcap key; board +0.0006, `callset_block` does not ship).
    A PIE-ranked block is a different key, but the same class: a DE prior on membership.
  - Its held-out price needs a non-seen-context DE truth: K562 GWPS cells (8,248 genes, dense h5ad) with per-target
    Wilcoxon calls. h1-test is seen-context for PIE, so it may only bound a cost (702 rule 2).
  - Not pursued before the 10-21 freeze: it needs the day-of PIE inference path (cloud evidence + context text) on top
    of a membership change. Recorded for the write-up and for any post-freeze work.

## PIE local data: what can be removed safely, and how to re-fetch it (2026-10-16T03:33Z, [d80a10])
The coordinator asked (10-16 03:3xZ, disk at 53 GiB) to delete PIE's local data and cache, a deliberate exception to
AGENTS §10's "keep data/external/" for re-downloadable model caches. At 03:32Z the disk showed **79 GiB free**, above the
60 GiB floor. Nothing has been deleted: that exception is Ahmet's call. This record makes the deletion safe and reversible
whenever he makes it.
- **Removable (~89 GB):**
  - `data/external/pie/data` (68 GB): PIE's Hugging Face datasets, under `data/hf/datasets/arcinstitute/<name>/<revision>/`:
    - pie_sources cb1aaa4e7655605bdc70a9bd77bbd62016b8c7d7 (41 GB)
    - pie_tahoe100m 597d6ded9f9a148331388484ffb193489121c702 (20 GB)
    - pie_x_atlas_orion 0950c8aa5ce5b33dc9d6f0c5eade3ae445ca0ffb (5.8 GB)
    - pie_replogle_nadig_essential 20c9faef76fc96fdc809871340bd93499413dfe7 (727 MB)
    - pie_jiang 99242356dfcf58e049d867703c5743df8d1790bf (210 MB)
    - pie_arc_vcc_25 e7dda6065316959958b87583ddf0a15f6e69c37d (47 MB)
  - `data/external/pie/cache` (21 GB): the evidence cache 87a01957365041383ac2cbae3a9cf44be6f48e5dbc69ca045f14fa270441d151.
- **Re-fetch:**
  - Datasets: `hf download arcinstitute/<name> --repo-type dataset --revision <revision> --local-dir …`, or
    `work/fetch_assets700.py` (the layout PIE expects).
  - Evidence: `aws s3 sync s3://<your-s3-bucket>/output/700-pie/20261014T225548Z/evidence/87a01957365041383ac2cbae3a9cf44be6f48e5dbc69ca045f14fa270441d151/ data/external/pie/cache/evidence/87a01957365041383ac2cbae3a9cf44be6f48e5dbc69ca045f14fa270441d151/ --region us-east-1`
    (64 objects, 22.85 GB, checked 10-16 03:32Z).
- **Keep (small, used by the closed-line write-up and any re-run):**
  - `runs/` (922 MB, incl. replogle_xdataset/best_auprc.ckpt);
  - `abc_ctrl/` (2.4 GB), `h1_ctrl/` (903 MB), `local_sources/` (23 MB: the CCLE context texts), `splits/`, `prep/`;
  - the cloud predictions in `work/cloud/`.
- **What it costs:** none of 702 / 705 / 707 / 708 / 710 or the October recipe reads these directories (PIE is closed;
  `dump_hook` is off). A PIE re-run needs the re-fetch first (~1 h download).

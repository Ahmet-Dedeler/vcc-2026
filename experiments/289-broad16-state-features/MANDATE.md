# 289 — Does the line-specific response component on Broad-16 become predictable when the STATE features are the source paper's own (p53 status, HPV-E6 status, lineage, proliferation programme, copy-number where available) instead of two basal expression PCs? 268's leave-one-line-out rerun on the 12 working lines, preregistered, plus the paper's two mechanistic predictions tested on our data. Analysis only, ≤ 4 hours. No build, no submission.

You are a research subagent. Work in `experiments/289-broad16-state-features/`. **Edit nothing outside it.** No build,
no `vcc`, no `cell-eval2`, nothing under `submissions/` or `data/`; small public tables (≤ 1 GB total, cited, sha256
recorded) may be fetched into `work/`. `uv run python`. Peak RSS under ~8 GB. **Time box: 4 hours.** Three sibling
experiments (284, 285, 286) run in parallel; do not read or write their directories. `docs/RESEARCH-PROGRAMME-2026-09.md`
binds you (preregister before any held number; name the unit; the three controls; rule 9: depth-condition pass rules).
**The four reserved Broad-16 lines (JVE127, MG63, SLR23, UMRC3) are SPENT (claim C2) and must not be read** — `load268.py`
asserts this; keep the assertion. Every number is on the 12 working lines under leave-one-line-out.

## Why (287's audit, contradiction C-1; `experiments/287-novelty-audit/NOVELTY-AUDIT.md` rows 9–10, §3)
268 (F5) found a reproducible, high-rank line-specific response component and reported it "not basal-state-predictable"
because a linear ridge on TWO basal expression PCs gained nothing on the full response under LOLO (predicted-residual
cosine +0.018 vs +0.001 wrong recipient; full-response Δcos −0.007 per gene, +0.000 programme-level; R268 §4). The
source paper (Maffa et al., *Dissecting context-dependent cancer vulnerabilities using Perturb-seq*, bioRxiv
10.64898/2026.08.24.746802, 2026-08-24; the dataset's own authors) predicts line-specific responses from cell STATE in
the mechanistic sense: p53 activation after UBE3A / MDM2 knockout only in the HPV-E6 cervical lines (SKGII, C4I), and
IER3IP1 dependence from a hemizygous 18q / 11q co-deletion. The two statements can both hold (transcriptome-wide linear
model on expression PCs vs per-gene predictions from genotype and viral status). F5's headline has been narrowed to "not
predictable from basal EXPRESSION by a linear model at 11 training lines" (FINDINGS v1.2). **This experiment measures
whether the paper's kind of state feature predicts the residual where expression PCs did not — on our data, our
statistic, our controls.**

## READ FIRST
268's `RESULTS.md` §2 (the decomposition, the split-half ceiling), §4 (the LOLO protocol, arms (a)–(nn), the WRONG-
recipient and SHUFFLED-label controls, the ceiling, the paired target bootstrap; table at lines 185–195) and its code
`work/load268.py` (halves, `build_responses`, `control_profiles`, META), `work/lolo268.py` (the fold loop: refit
everything inside the 11 training lines; K basal PCs; ridge λ = 1; J = 5 residual programmes; `retrieval`), `work/
lolosum268.py`; 283's `RESULTS.md` (C2c: 268's basal model on the reserved lines — do NOT re-run it) and `work/
lib283.py`; 287's audit rows 9–10, 24 and §3 C-1 / T4; `docs/FINDINGS-2026-09.md` F5 (v1.2). Data: 268's
`work/halves/<line>.npz` (12 working lines; cell halves A / B, guide halves v1 / v2, full), `data/external/broad16line/
metadata/cell_line_metadata.csv` (line, ArxSpan id, batch, OncotreeLineage), `hallmark_symbols.gmt` (MSigDB hallmark
sets incl. P53_PATHWAY, E2F_TARGETS, G2M_CHECKPOINT), `experiments/269-utility-router-feasibility/work/depmap/
CRISPRGeneEffect.csv` + `Model.csv` (DepMap 24Q2 Chronos gene effects for 1,150 models — the 16 lines' ArxSpan ids
are in the metadata), `data/external/ccle_broad_2025/` (expression TPM; `data_clinical_sample.txt` with Oncotree
lineage and TMB), `data/external/depmap26q1/extracted/processed/` (the paper's own per-KO pseudobulks / LFC / FDR
matrices, `cell_quality_covariates.csv` = UMIs per line).

## The frozen state features (write them in `PREREGISTRATION.md`, hashed, before any LOLO number)
Build-time observable for each of the 12 working lines, none derived from the perturbation responses:
1. **p53 functional status** — proxied from DepMap gene effects on disk: a line with a strong MDM2 dependency and a
   non-negative TP53 effect is p53-wild-type; a line with no MDM2 dependency and/or a TP53 knockout benefit is p53-
   deficient (state the exact rule and thresholds before looking at responses; report the per-line calls; where the
   paper or DepMap's mutation annotation names TP53 status directly and is reachable, record agreement).
2. **HPV-E6 status** — the paper's HPV-E6 cervical lines (SKGII, C4I) = positive; every other line negative unless the
   paper states otherwise (cite the sentence).
3. **Lineage** — OncotreeLineage from the metadata (11 lineages over 12 lines: use it only as the "nearest same-
   lineage line" arm and as a shuffled control, never as a one-hot in a ridge with 11 training rows).
4. **Proliferation and stress programmes of the BASAL state** — hallmark E2F_TARGETS + G2M_CHECKPOINT and
   P53_PATHWAY scores of each line's control (OR-KO) profile (mean log-CPM over the set minus the all-gene mean), fitted
   inside the training lines as 268 fitted its PCs.
5. **Copy number, arm-level** (the paper's 18q / 11q co-deletion feature) — ONLY if a public DepMap arm-level or
   segment CN table for these ArxSpan ids is reachable at ≤ 1 GB (cite; sha256); otherwise declare it absent and run
   without it. Do not spend more than 30 minutes on acquisition.
6. **Sequencing depth per line** (UMIs per cell from `cell_quality_covariates.csv` / 268's inventory) — a NUISANCE
   covariate, reported alongside every result (rule 9), never a predictor.

## Protocol
1. **Known answer first:** re-run 268's LOLO arm (c1) with K = 2 basal PCs from `lolo268.py`'s code path and reproduce
   R268 §4's row at 0.0e+00 (predicted-residual cosine +0.0181 [+0.0134, +0.0229]; full-response Δcos −0.0071; retrieval
   +0.0003). Abort if it does not reproduce.
2. **State-feature LOLO (frozen arms):** the same fold loop, the same ridge (λ = 1) and programme-level (J = 5) models,
   with the basal PCs replaced by (S1) the frozen state vector {p53 status, HPV-E6, E2F/G2M score, P53_PATHWAY score
   (+ arm CN if present)}, (S2) each state feature alone, (S3) state vector + the 2 basal PCs; the same three
   readouts (cosine of the predicted residual to the held residual; gain on the FULL response, per gene and
   programme-level; within-line retrieval) with the same controls — WRONG recipient (every other training line's
   features), SHUFFLED labels (features permuted across training lines, 5 permutations), the split-half ceiling —
   paired target bootstrap over targets, 12 folds. Pass rule (frozen): an arm "predicts the residual from state" if
   its predicted-residual cosine beats BOTH controls with the paired interval above 0 AND its full-response Δcos is
   > 0 with the interval above 0 on the programme-level readout.
3. **The paper's two mechanistic predictions, as preregistered per-target tests on the 12 working lines** (state each
   as a claim with a rule before computing): (P1) after UBE3A and MDM2 knockout, the P53_PATHWAY hallmark score of the
   line-specific residual is higher in the HPV-E6 lines (SKGII, C4I) than in the p53-wild-type non-HPV lines and than
   in the p53-deficient lines — report the per-line values, the ordering, and a permutation over line labels; (P2) if
   IER3IP1 and the 18q / 11q features exist on our lines, the paper's co-deletion lines show the stronger IER3IP1
   knockout response (fitness / transcriptome), else state that P2 is untestable here and why. Say whether each
   reproduces in OUR residual object, not only in the paper's DE tables.
4. **Per-target gating (the mechanism-specific model class), frozen:** for each target, LOLO-predict the SIGN of the
   held line's residual projection on the P53_PATHWAY and E2F/G2M directions from the p53 / HPV / proliferation
   features (a 1–2-parameter logistic or sign rule per target); report the held-line accuracy against the shuffled
   control and the fraction of targets whose gating beats chance — the reading that says whether "state predicts
   response" is true target-by-target even where the transcriptome-wide ridge gains nothing.
5. **Depth as a covariate (rule 9):** every fold's readout next to the held line's UMIs per cell; Spearman of the
   per-line gain with depth; state whether the state-feature gain, if any, is concentrated in the deep lines.
6. **Descriptive (post-hoc, labelled):** the 20 targets whose line-specific residual is best predicted from state, and
   which feature carries each.

## Deliverable
`PREREGISTRATION.md` (hashed into `work/PREREGISTRATION.sha256`), `RESULTS.md`: the KA row; the per-line state-feature
table with the rule that produced each call; the LOLO table (S1 / S2 / S3 next to 268's (c1) / (c2) / (nn), with both
controls and the ceiling); P1 / P2; the per-target gating result; the depth covariate; the descriptive tail;
`work/state289.csv` (per line: features, depth) and `work/lolo289.csv` (per fold × arm); one line: **STATE FEATURES
PREDICT THE LINE-SPECIFIC RESIDUAL: YES / NO (S1 predicted-residual cosine … vs wrong-recipient … / shuffled …;
full-response Δcos …); PAPER'S P1 REPRODUCES ON OUR DATA: YES / NO; P2: YES / NO / UNTESTABLE; PER-TARGET GATING beats
chance on … % of targets**. No model use for the competition is proposed. Then stop.

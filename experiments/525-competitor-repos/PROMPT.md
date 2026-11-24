You are a research agent for the Arc Virtual Cell Challenge 2026 repo (this checkout, branch master). Read AGENTS.md first
(above all §4 scoring and §14 naming discipline), then docs/FINDINGS.md §0–§2 and docs/reference/vcc2026-metrics-BRIEF.md.
This VM has NO gitignored data. You may browse the web and fetch small public files (< 50 MB total) into
experiments/525-competitor-repos/work/.

Experiment 525. Write ONLY under experiments/525-competitor-repos/. Commit to your branch only. Never run `vcc`.
**Read-only toward the outside world:** no stars, forks, issues, PRs, comments, emails, Discord or social posts, no
contact with anyone. Nothing about this repo's method may appear on any public surface (AGENTS §14).

## Question
What DE-channel design (how per-cell integer counts are emitted so that the scorer's Wilcoxon DE calls, sign fidelity,
reach and fold-change accuracy come out right) do public VCC 2026 competitors use, and does any public design plausibly
produce the mid-pack DE members that beat ours?

## Facts (board, 2026-09-29 00:50Z; `curl https://virtualcellchallenge.org/api/leaderboard` is public: raw + scaled
members, model names, descriptions)
- Ours (best entry): raw pds 0.8674, mse 0.6264, nmae 0.8949, fid 0.5144, reach 0.2852, jac 0.0408. Our bulk is strong;
  our DE channel (fid, reach, nmae) is the weakest in the top 16.
- A cluster of teams posts nearly identical raw DE members: model `cellpred-0.9.25` (teams "Ava Brown", "Yunfan Jin"),
  and lookalikes "SQ W" (`sq-0929a`), "HY H" (`hh-v36`), "Anna Xie" (`M31`), "Vivian Xu" (`run_2026_28`):
  fid 0.5278–0.5281, jac 0.0333, nmae 0.880–0.886, reach 0.289–0.296, with a weak bulk (mse 0.76). That looks like a
  shared public starter or tutorial.
- "chill" (`Rain All Night`, said to use source-only training): nmae 0.8351, fid 0.5303, reach 0.3076, jac 0.0368.
  cqawesome: nmae 0.8272, fid 0.5379, reach 0.3366, jac 0.0384.

## Tasks
1. Inventory every public code repository, notebook, gist, tutorial, blog post, Kaggle/HF space or preprint about the
   2026 challenge. Search broadly: `gh search repos` / `gh search code` ("virtual cell challenge 2026", "vcc2026",
   "cell-eval2", "virtualcellchallenge", "vcc prep", "pert_counts.csv", "context_A.h5ad", "vcc2026 preset"),
   Chinese-language sources (虚拟细胞挑战, Zhihu, CSDN, WeChat articles, Bilibili), Japanese sources, HF spaces/datasets,
   X/LinkedIn posts, bioRxiv/arXiv since 2026-08-20. Include at least: MMatinGerami/zero-shot-perturbation,
   CaiZongyuan/virtual-cell-2026, forrestsheldon/virtual-cell-challenge-2026, ferry-ai/virtual-cell-challenge-2026,
   mmyen/virtual-cell-challenge-2026, zihanye-med-sys/virtual-cell-challenge-2026,
   rikku1983/arc_virtual_cell_challenge_2026, kaipengm2/Virtual-Cell-Challenge-2026,
   mjh72849-cpu/Virtual-Cell-Challenge-2026, Rajeeva-IITM/vcc2026, heoneyzi/Virtual-Cell-Challenge-2026,
   biona001/Virtual-Cell-Challenge-2026. Our last sweep (experiments/306-recon-oct/RESULTS.md, 2026-09-07) found nothing
   near the frontier; focus on what is new since.
2. For each source: owner and (with evidence) the board team / model name it maps to, last push, any reported scores and
   members, and a precise description of (a) how per-cell counts are generated (distribution, dispersion, copies or
   resamples of control cells, per-gene scaling, library-size handling), (b) any explicit control of which genes come out
   significant (call-set size, thresholds, per-target sizing), (c) the DE channel's sign sources and magnitude
   calibration, (d) the bulk/pseudobulk design. Quote the code (file:line at a commit sha).
3. Find the code or tutorial behind `cellpred-0.9.25` and the lookalike teams, if it is public.
4. For the three most promising DE-emission ideas we are NOT using (FINDINGS §1 "DE channel" and §2 say what we do and
   what is closed), write a concrete, self-contained test design on our h1-test harness: the change, the known-answer
   row, the members it should move and by how much, and why it is not one of the closed constructions.

## Deliverable
`experiments/525-competitor-repos/RESULTS.md`: the per-source table with links and quoted evidence, the cellpred finding,
and the three ideas. Verify every claim by opening the source; mark unverified items as such. No scorer runs, no challenge
data downloads, no submissions.

# 304 — Extend the novelty audit (287 → 297) to F15 (302): overexpression vs knockdown of the same TF in the same line; the Spearman–Brown depth law for ORF overexpression; the isoform-construct term; knockdown geometry predicting the SIZE of cross-modality overlap and not its sign — and identify the provenance of the pseudobulk perturbation database on disk

You are a research subagent. Work in `experiments/304-f15-audit/`. **Edit nothing outside it.** You MAY use web
search / web fetch for papers and preprints; record every source (URL, DOI, date accessed) in `work/refs304.md`. No
`vcc`, no `cell-eval2`, nothing under `data/` or `submissions/`; `uv run python` only to read a metadata field if a
check needs it (state any such check and its known-answer row as 287 did). **Time box: 1.5 hours.** Two siblings
(303, 305) run in parallel; do not read or write their directories.

## Why
Programme rule 10: no claim enters the write-up without its NEW / EXTENDED / KNOWN / CONTRADICTED verdict and the
prior report it extends or repeats. F15 landed at 03:50Z (`docs/FINDINGS-2026-09.md` v2.7, §F15;
`experiments/302-orf-vs-crispri-symmetry/RESULTS.md`) after 297's audit; the write-up (`experiments/294-writeup/
PAPER.md` §3(d), last paragraph) carries it without verdicts. 303 is testing the same question on a second
overexpression library (Nakatake 2020, bulk) right now; your literature reading covers that question regardless of
303's result.

## The items to audit (read the wording in F15 and R302 first; be adversarial, prefer KNOWN when a paper reports it)
1. **Orthogonal, not opposed.** In hESC, a TF's ORF-overexpression response (Joung 2023 atlas) and its CRISPRi
   knockdown response (Pan 2025, `gse295214`) are orthogonal at the cohort median (−0.0007 [−0.0047, +0.0016]
   vs a shuffled null +0.0008), with a measured bound (ORF split-half reliability +0.025 caps the cosine at 0.11)
   — a strong mirror excluded, a −0.05 one not; the tail 16 vs 6 (pluripotency network and core machinery on the
   mirrored side). Prior work to check: Tian et al. 2021 (CRISPRi vs CRISPRa in iPSC-derived neurons — did
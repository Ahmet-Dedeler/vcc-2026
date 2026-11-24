# 275 — The scientific track's first deliverable: a findings note, numbers-first, on what the 2026-09-05/06 programme established about perturbation-response transfer. Writing and re-verification only. ≤ 3 hours. No new measurement, no build, no submission.

You are a research subagent. Work in `experiments/275-findings-writeup/`. **Edit nothing outside it** except the one
deliverable file `docs/FINDINGS-2026-09.md` (create it). No build, no `vcc`, no `cell-eval2`, no new fold runs. `uv run
python` only to re-read saved JSON/CSV numbers. Peak RSS under ~2 GB. **Time box: 3 hours.**

## Why
`docs/RESEARCH-PROGRAMME-2026-09.md` keeps two tracks; the scientific one needs its results written as findings that
stand on their own — with the experimental unit named, the number, its interval, its limit — not as a chain of
negative competition verdicts. Seven experiments produced findings that outlive their verdicts. Write them up so a
reader who never saw the repo can judge them, and so the coordinator can decide what is worth confirming on the
reserved data (`docs/CONFIRMATION-LEDGER.md`) and what is worth a paper.

## READ FIRST (and re-verify every number you quote against the file it came from; cite the file and section)
`docs/RESEARCH-PROGRAMME-2026-09.md`; `docs/READ-FIRST.md` v4.2 §4–§6; the RESULTS.md of **267** (reference reliability;
the replicate-validated oracle), **271** (donor-row reliability; noise-aware estimation; the 400-cell statement),
**272** (the operator's in-sample noise learning; the spectrum; retrieval needs the noise-dominated directions),
**273** (the exact instrument; Flex-400 reliability; the guide construct; the conversion curve; the regime),
**262** (per-target `pds` structure; the conversion curve; the deficit decomposition), **263/265** (the oracle ceiling and
its permutation test), **268** (Broad-16 variance decomposition), **266/260** (K562 DE coverage); BR237 and BR249 for the
board points the instruments were calibrated on; `experiments/159-transfer-ceiling/RESULTS.md` §1 (the earlier
reliability laws these confirm or correct).

## The note (`docs/FINDINGS-2026-09.md`, ≤ 2,500 words plus tables)
1. **Setting** (one paragraph): the task (predict 300 perturbation responses in held cell lines from donor lines), the
   instruments, the two tracks, what "board-validated" means here.
2. **Findings**, each as: claim → the measurement (unit, n, statistic, interval) → the limit → what it changes. Cover:
   (F1) Perturb-seq response reproducibility at 150–200 cells on a non-essential panel (267 §3, 271 §1: medians, share
   ≤ 0, batch vs random halves, Spearman–Brown; 273: Flex-400 +0.243, 0/85 ≤ 0; the guide construct +0.014).
   (F2) The donor-selection oracle is reference noise (267: −0.018 cross-validated vs +0.056 in-sample; 265: p = 0.011
   on one fold, at chance on three); what survives (the reliable tercile).
   (F3) The whitening operator learns its training noise (271 §1.2, 272 §1: 77–79 % noise-dominated directions
   carrying 87–92 % of whitened energy; cross-fit +0.02 vs −0.14/−0.24) and why retrieval needs those directions
   (signal 5–10 % of row energy spread over them).
   (F4) Per-target `pds` is a threshold in the transfer cosine at ≈ 0.02 — a property of the rank statistic (273 §3),
   with the regime shares on PLOLO and on Flex-400; the deficit decomposition (262 §5).
   (F5) Context dependence on Broad-16: a reproducible line-specific component (11.5 %, 4× null) that is high-rank and
   not basal-state-predictable at 11 training lines (268).
   (F6) DE-set sparsity on K562-GWPS for a non-essential panel (260/266: median 5–6 significant genes, unchanged by
   matched controls) — what it means for local DE metrics.
   (F7) The board calibrations: the `nmae` identity's sign and magnitude (BR237), the PLOLO→board class factors and
   their sensitivity (BR249, 257–259), the "8 σ" correction (265).
3. **What is and is not established** — a table: finding / unit resampled / n / independent replication? / could be
   confirmed on the reserved data (267's HEK293T half-B, 268's four lines) with which preregistered test.
4. **Implications** for the two tracks in ≤ 10 lines each, no speculation beyond the numbers.
5. **Methods appendix**: the exact statistic (273 §1), the split units, the bootstrap and permutation conventions, the
   pessimistic-pricing rule, the leakage rules used.

Style: numbers first, every one sourced (`file:section`); intervals as [lo, hi]; no adjectives about significance the
numbers do not support; state every limit (H1 embryonic; three cancer lines; one instrument; 85 targets).

## Deliverable
`docs/FINDINGS-2026-09.md` (the note) and `RESULTS.md` in this directory (a checklist of every number verified against
its source, and any discrepancy you found between a RESULTS.md headline and its own tables — report those, do not
silently fix them). No model name, no submission. Then stop.

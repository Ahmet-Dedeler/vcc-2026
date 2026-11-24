# 306 — Reconnaissance, 2026-09-07: what data and what published method exists TODAY that was not available when this repo's data was frozen (2026-08-30), and what is measurably true about the methods now scoring 0.22–0.24 on the board

You are a research subagent. Work in `experiments/306-recon-oct/`. **Edit nothing outside it.** You MAY use web
search and web fetch freely; record every source (URL, DOI/accession, date accessed, and the exact sentence you are
relying on) in `work/sources306.md`. **No download into `data/`**; anything you fetch goes under
`experiments/306-recon-oct/work/` and must be under 200 MB total. No `vcc`, no `cell-eval2`, no build, nothing under
`submissions/`. `uv run python` only. **Time box: 3 hours.** Siblings run in parallel; do not read or write their
directories.

## Why
The competition track is closed at every construction around the incumbent (`docs/READ-FIRST.md` §5;
`docs/RESEARCH-PROGRAMME-2026-09.md`, "State of the competition track after Round 4"). The programme's own conclusion
is that the next gain needs **new data** or **a new modelling idea with a stated mechanism**, not another pass over
the same constructions. The last reconnaissance (270, `experiments/270-data-reconnaissance/RESULTS.md`) was
2026-09-05 and the data on disk was frozen 2026-08-30. `cd4svd5` is 0.2160944, rank 3/762; the leader `cqawesome` is
0.238932, a gap of +0.0228. **Your job is to find what changed and what is knowable — not to propose a model.**

## READ FIRST
`docs/READ-FIRST.md` (all of it; §5 closed families, §6 open items), `experiments/270-data-reconnaissance/RESULTS.md`
(§4.3 on X-Cell/Pisces, §5's recommendations, and its full on-disk inventory table — **do not repeat what it already
settled**), `AGENTS.md` §5 (data) and §3 (timeline), `docs/FINDINGS-2026-09.md` §3 (what is and is not established).

## Part A — new data since 2026-08-30 (the primary deliverable)
For each item: is it released, what exactly does it contain (lines, targets, cells/depth, perturbation type,
chemistry), does it touch the 300-target panel or the three held lines, what is the licence, what is the download
size, and **what would it buy against a named open question** (F5 needs more LINES; F14 needs a SECOND CONSTRUCT per
target; 293 needs the challenge's own library in a public screen; 268's stage 2 needs contexts).
1. **X-Atlas/Pisces** (`huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Pisces`) and **X-Cell** — 270 found
   `usedStorage 0` and "Coming Soon" on 2026-09-05. Re-check the dataset repo, the model repo, the GitHub org and any
   preprint. State the `usedStorage`, file list and last-modified date you actually observe.
2. **Arc Institute** — any release, blog post, dataset drop, model or challenge update since 2026-08-30 (Virtual Cell
   Challenge pages, `arcinstitute.org`, their GitHub org, the STATE model repos and weights).
3. **GEO / ArrayExpress / scPerturb / GPSAdb / Zenodo / figshare** — new Perturb-seq or CRISPRi/a screens deposited
   since 2026-08-01, especially: any screen in **K562, HCT116, HEK293T, Jurkat, HepG2, H1/hESC, primary CD4** or any
   line with ≥ 50 of the 300 panel targets; any screen with **two or more constructs per target**; any screen using
   the **2025 VCC library** (its spacers are in `experiments/301-metadata-acquisitions/work/vcc2025_spacers.csv` —
   you may read that file).
4. **Preprints since 2026-08-01** (bioRxiv/arXiv) on perturbation-response transfer, virtual cell models, or the
   challenge itself. Note any that report a number on this challenge's metric.

## Part B — the board's frontier, from public sources only
The challenge's own leaderboard, forum/Discourse/Discord if public, GitHub repos that reference the Virtual Cell
Challenge or `cell-eval`, and any write-up by a named team. Questions: what methods are teams reporting; has anyone
published a score near 0.24 with a method description; is there any public statement from Arc about the
October 22 test release or the scoring. **Public sources only — no accounts, no logins, no scraping behind auth, and
do not contact anyone.** If a source is a forum post, quote it and attribute it; treat every claim in it as a claim,
not a fact, and say so.

## Part C — one falsifiable check on something you find
Pick the single most consequential factual claim from Part A or B that can be checked against data **already on
disk** (inventory in 270 §4 and `data/external/`), and check it. State it as a known-answer row: what you expected,
what you measured, and whether it passed. Example shape: "dataset X claims N of the 300 panel targets" → count them
against `data/challenge/`'s panel list. Keep it under 30 minutes and under 4 GB RSS.

## Rules of evidence
- A dataset is **RELEASED** only if you observed downloadable files and their sizes; "announced", "in press" and
  "coming soon" are separate states. Say which you observed and paste the evidence.
- Distinguish **what a paper claims** from **what its deposited data contains**.
- Where a source is paywalled or unreachable, say so; do not infer contents from a title.
- Negative results are the expected outcome for most items and are valuable. Do not inflate.

## Deliverable
`RESULTS.md` with: (1) a one-line verdict; (2) the Part A table (item | state | contents | panel/line overlap |
licence | size | what it would buy | evidence URL); (3) Part B's findings with quotes and attributions; (4) Part C's
known-answer row; (5) **a ranked shortlist of at most three acquisitions worth a download**, each with the open
question it closes, its size, and what it would cost to use. `work/sources306.md` with every source. No model is
proposed, nothing is downloaded into `data/`, no submission is discussed. Then stop.

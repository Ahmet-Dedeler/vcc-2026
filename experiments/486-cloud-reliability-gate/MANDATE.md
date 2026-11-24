# 486 — Does a learned per-(target, gene) reliability gate transfer to held-out cell lines? (Cursor Cloud, Grok 4.7 xhigh)

Your experiment directory is `experiments/486-cloud-reliability-gate/`. This is a self-contained modelling study on public
data you download to your VM; it does not touch the production pipeline.

## Why

Our model predicts a new cell line's response to knocking down gene t by combining other lines' measured responses to t
(the "donor transfer"). The board's #3 team (Testh) disclosed a lever we lack: weighting each (target, gene) entry of the
transferred response by a learned reliability gate. On their board history a fixed-formula weight raised pds
0.8054 → 0.8203, and a learned gate plus a per-target weight took it to 0.8356 (mse 0.810 → 0.699). Their disclosed
gate: `w = 0.15 + 0.85 · sigmoid(linear in 8 truth-free features)`: log z² (source mean² / variance), log n_src, log
control expression, log |μ|, logit cross-source sign agreement, log cross-target variance of the gene, log n_cells, log
source MAD; coefficients fitted by Nelder–Mead on leave-one-context-out folds to maximise energy-weighted cosine. Their
fixed formula: `w = 0.15 + 0.85 · z²/(z²+8) · n/(n+4)`. The lead's team will separately re-price a hand-built version
inside the real pipeline (experiment 473). Your job is the clean public-data test: is there transferable signal in a
learned gate, and how much?

## Data (download to the VM; pseudobulk / "bulk" summaries are enough)

- Replogle et al. 2022 (figshare+ doi 10.25452/figshare.plus.20029387): K562 genome-wide, K562 essential, RPE1 essential
  (their "bulk" / normalised pseudobulk h5ad files, and per-perturbation cell counts).
- Nadig et al. 2025: HepG2 and Jurkat essential-wide CRISPRi Perturb-seq (find the public release: GEO or figshare).
- Optional if tractable on your disk: X-Atlas/Orion (Hugging Face `xaira-therapeutics/x-atlas-orion`) HCT116 / HEK293T.
Record exact URLs, file sizes and checksums in RESULTS.md.

## Method

1. Build per-line response matrices on a common gene axis: for each perturbation, the log fold change of its pseudobulk
   against that line's non-targeting control, plus per-entry variance / cell-count estimates where the data allow.
2. Folds: leave one line out. For the held-out line, predict each shared target's response from the other lines only.
   Evaluate on targets present in the held-out line and at least two others.
3. Rows (every table): (a) equal-weight mean of the other lines (the known-answer baseline); (b) Testh's fixed formula;
   (c) the learned gate fitted on the TRAINING folds only (nested: never fit on the held-out line); (d) controls: the
   learned gate with features permuted across genes, and with targets shuffled.
4. Metrics: mean cosine between predicted and true response (energy-weighted and unweighted), and a pds-like
   discrimination: for each target, the normalised rank of its own true response among all held-out targets by cosine
   distance, after removing the panel targets' own genes.
5. Report gains with target-bootstrap 95 % intervals (resample targets, 4,000 draws) per held-out line, the fitted
   coefficients, and which features carry the gain. Say plainly whether the learned gate beats the fixed formula and the
   equal-weight mean out of sample, and on how many held-out lines.

## Ground rules for every cloud teammate

- You run on a Cursor Cloud VM with this repo and its `uv` environment (`uv run python ...`; `cell-eval2` is installed).
  The gitignored data (`data/challenge`, `data/donors`, `data/harness`, `data/external`, `submissions/`) is NOT on your
  VM; do not assume it exists. Production runs on the lead's macOS machine (bash 3.2, BSD grep, APFS); you are on Linux,
  so reason explicitly about any shell or tool behaviour that differs (AGENTS.md §6.2 lists the traps).
- Read `AGENTS.md` and `docs/FINDINGS.md` §0–§2 first. The record is large: grep `docs/` and `experiments/*/RESULTS.md`
  before claiming something is new.
- Write ONLY inside your experiment directory. Commit only that directory to your branch. Do not open a PR, do not merge,
  do not edit production files (`scripts/`, `src/`, `docs/`, other experiments); put proposed patches as separate files
  in your directory.
- Nothing public: no posting, issues, PRs, accounts, emails or messages. Never run `vcc`.
- Evidence standard: every claim carries a reproduction (a test, a script, a command and its output) that the lead can
  re-run. Label each finding CONFIRMED (reproduced) or SUSPECTED (reasoned). A well-supported negative is a result.
- Commit `RESULTS.md` early and update it as you go; finish with a complete `RESULTS.md` whose first lines give the
  answer. Budget: up to ~8 hours.

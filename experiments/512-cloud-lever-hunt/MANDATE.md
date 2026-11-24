# 512 — fresh-model lever hunt: what could close the DE-channel gap that the record has NOT tried? (cloud agent, code and record only)

Repo `Ahmet-Dedeler/virtual-cell-challenge-2026` (private), branch `master`. No gitignored data on this VM (no challenge
controls, donors or harness surfaces); work from the repo's text, code, tracked board snapshots and public literature.
Never run `vcc`, never touch credentials, never post anything anywhere outside this repo (AGENTS.md §14: nothing about the
method may become public). Commit ONLY to `cloud/512-lever-hunt`, ONLY inside `experiments/512-cloud-lever-hunt/`. No PR.

## Why

We are 7th of ~1,200 teams on the validation board (best entry 0.2723 overall; Illumina AI 0.357, cqawesome 0.321,
Testh 0.306). The final is scored on a new panel (contexts D/E/F, released Oct 22; only the last entry before Nov 5
counts), so only levers that transfer to a new panel and new cell lines matter. About 80 % of our gap to the top two sits
in the four DE members:

| team | nmae raw (lower better) | fid raw | reach raw | pds | mse raw |
|---|---:|---:|---:|---:|---:|
| Illumina AI | 0.814 | 0.547 | 0.428 | — | — |
| cqawesome | 0.827 | 0.538 | 0.337 | — | — |
| chill | 0.835 | 0.530 | 0.308 | — | — |
| the "pack" (5 near-identical accounts) | 0.89–0.91 | 0.528 | 0.29 | — | — |
| **us** | **0.896** | **0.514** | **0.285** | — | — |

(Members for every team and every day are in `data/leaderboard/*.json`, tracked.)

The record's current reading (read these before proposing anything):
- `docs/FINDINGS.md` §0–§2: current state, conclusions, and **§2, the closed list** (do not re-propose without a new
  mechanism).
- `docs/reference/vcc2026-metrics-BRIEF.md`: the exact metric definitions.
- `experiments/480-nmae-instrument/RESULTS.md`: the nmae gap is sign information on weak responses and uncalled genes
  (oracle signs at our magnitudes recover −0.254 on H1; oracle magnitudes −0.0985; our sign accuracy on uncalled gate genes
  0.54); no class-based magnitude rule beats the shipped MDE.
- `experiments/192-predicted-callset/RESULTS.md`, `experiments/385-per-target-callset-size/RESULTS.md`: fid is closed to
  the call set (precision 0.507–0.517 for every key at every size; n_real not predictable).
- `docs/BOARD-INTEL-2026-09-26.md`, `docs/H1-BOARD-CALIBRATION.md`, `experiments/434-*/RESULTS.md`,
  `experiments/478-*/RESULTS.md`, `experiments/496-uncalled-sign-hook/RESULTS.md`, `experiments/462*`, `experiments/467*`,
  `experiments/470*` (programme priors), `experiments/423*`, `experiments/426*`, `experiments/429*` (DE magnitude delivery).
- The decisive pending test: 09-30's board read of Arc's own H1 programme (per-gene P_down from Arc's public 2025 H1
  screen) as every context's DE block prior (`experiments/458-arc-lab-programme`, `experiments/468-prior-gate`).

## What to do

1. **Board forensics (tracked data).** From `data/leaderboard/*.json` (dedupe entries by `id`), reconstruct the raw and
   scaled member trajectories of the top 10 distinct teams. List every step where some members stayed bit-identical and
   others moved (like cqawesome's nmae-only steps, Testh's 09-23 DE-only jump, chill's first entry at nmae 0.835 with a
   weak bulk). For each, use the metric definitions to bound which constructions could produce it (e.g. an nmae move with
   pds, mse, fid, jac identical means arithmetic means moved with the call set and every call's sign unchanged). Say which of
   those constructions our record has tried, citing the experiment.
2. **Structural levers of the scorer.** From the BRIEF alone, list properties of the scoring that reward a construction
   independent of biological information (we already exploit: ballast cells that set the pseudobulk at 1/400 weight in
   means and ranks; unwelded per-cell libraries; exact-U mass moves; the MDE add; the zero-gene repair). Find any we have
   not used, with the member it moves and why it is legitimate (a fitted-model prediction, not a lookup of real effects).
3. **Literature and public methods, 2025–2026.** Methods for predicting perturbation response direction in unseen cell
   lines (State, scGPT, GEARS, linear baselines, programme or module models, VCC 2025 winners' write-ups, anything on
   bioRxiv/arXiv). Extract only what is new relative to the record: an information source, a model class, or a training
   target we have not used.
4. **Proposals.** Up to six, ranked by expected transfer to a NEW panel and new lines. Each: mechanism tied to a metric
   definition; the information source; which members move and which way; the closest item in FINDINGS §2 and exactly
   why this is different; the existing script that would price it on h1-test (e.g. `experiments/480-nmae-instrument/*`,
   `experiments/496-uncalled-sign-hook/run496.py`, `experiments/385-per-target-callset-size/run385.py`) and what to change;
   a prior on the size in overall units with a reason. A proposal that is a renamed closed item is worse than no proposal.

## Return

`experiments/512-cloud-lever-hunt/RESULTS.md`: (1) the forensics table and step list with bounds; (2) the structural
list; (3) the literature findings with links; (4) the ranked proposals. Put any analysis scripts next to it. Mark every
claim VERIFIED (you computed or read it) or REPORTED (secondary). Short sentences, numbers first.

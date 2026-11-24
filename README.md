# Virtual Cell Challenge 2026

Solo entry, research archive, and experimental codebase for the [Arc Institute Virtual Cell Challenge 2026](https://virtualcellchallenge.org) (Aug-Nov 2026): zero-shot prediction of single-cell transcriptional responses to CRISPRi knockdowns across unseen human cell lines.

I finished **Rank 59 out of 800+ final competing teams** (5,400+ registrants, top ~7%). Everything here was built on a single MacBook (Apple M5 Max, 64 GB RAM) plus disposable CPU EC2 instances for offline `cell-eval2` scoring, with **zero GPUs and zero neural networks**.

I am open-sourcing the full research record (code, 650+ experiment write-ups, preregistrations, negative results, submission metadata, and 313 validation-phase leaderboard snapshots) because I am not planning to work in this area in the near term, and the negative results, metric mechanics, and validation tooling are worth more out in the open than sitting in a private archive.

## Final Result

Counting final submission: `X2 = R0 + 811` (entry ID `S5RxOrm4ZbqtnPRIX8NX`, submitted 2026-10-22 on test contexts D = MOLM13, E = U87MG, F = SKLMS1).

| Metric | `cell-eval2` key | Scaled Score (Test D/E/F) |
|---|---|---:|
| **Overall** | unweighted mean (6 metrics x 3 contexts) | **0.2441** |
| Perturbation discrimination (`pds`) | `pds_cosine` | `0.7636` |
| Expression accuracy (`mse`) | `expr_mse_unbiased_capped_norm` | `0.2726` |
| DE direction reach (`reach`) | `de_wilcoxon_direction_reach_raw` | `0.1879` |
| DE log-FC accuracy (`nmae`) | `de_wilcoxon_lfc_nmae` | `0.1494` |
| DE direction fidelity (`fid`) | `de_wilcoxon_direction_fidelity_yield_raw` | `0.0565` |
| DE significance overlap (`jac`) | `de_wilcoxon_sig_jaccard` | `0.0348` |

Peak validation-phase score on contexts A/B/C (Jurkat, HeLa, CAL33): `0.2781` (`C0`). For context on the final test board, the top three winning teams were #1 Illumina BioInsight AI Lab (`CellularAI`, ~0.440), #2 FlowState (~0.335), and #3 Xaira (`X-Cell Transfer`, ~0.326).

## How the Model (`X2`) Works

There are no neural network weights or checkpoints. The final model is a deterministic, linear-algebraic donor-transfer pipeline fitted on public Perturb-seq screens and calibrated per context from the challenge's own non-targeting control cells. Full technical specification: [`docs/FINALIST-DESCRIPTION-DRAFT.md`](docs/FINALIST-DESCRIPTION-DRAFT.md) and [`docs/october-decisions.json`](docs/october-decisions.json).

1. **Public donor stack (`src/vcc2026/donor_transfer.py`, `src/vcc2026/multi_donor.py`)**: Per-target pseudobulk response vectors across public CRISPRi screens (Replogle K562/RPE1, Nadig HepG2/Jurkat, X-Atlas/Orion HCT116/HEK293T, Pan 2025 H1, VCC 2025 H1, primary CD4+ T cells, DLD-1, Jiang 2025, Sunshine 2023, Wessels, Myllia), mapped onto the 18,533-gene axis.
2. **Shrinkage whitening (`src/vcc2026/whiten.py`)**: Whitening operator fitted on donor signal covariance (shrinkage gamma = 0.65) to transfer cross-line perturbation structure, combined with an H1 residual union (weight 2x) and a low-rank deflated primary CD4+ T-cell residual (weight 1x).
3. **Per-target magnitude key (`scripts/def_lodo_key.py`)**: Leave-one-donor-out shared-energy key measuring cross-donor agreement per target, capped at the 96.6th percentile, plus a 0.25 strength rider.
4. **Tercile sharpening**: Row-norm-preserving magnitude sharpening `sign(h) * |h|^q` with `q = 2.5 / 2.0 / 1.75` across the strong, middle, and weak key terciles.
5. **18-screen many-line consensus blend (`experiments/811-lolo-model/`)**: Norm-preserving blend (`lambda = 0.25`) toward an equal-weight cross-line consensus direction extracted across 18 public screen banks. Validated on D/E/F leave-one-line-out folds (+0.0045, 95% CI [+0.0018, +0.0073]) and a sealed holdout reserve (+0.0100, 6/6 folds) before release-day submission.
6. **Per-context noise-calibrated levelling (`scripts/def_day.sh`)**: Iteratively levels each unseen context's effect scale against a noise target computed strictly from that context's own non-targeting controls (to stay on the safe side of `mse`'s asymmetric cliff). On D/E/F this converged in 3 levels (factors D: 1.0757, E: 0.8924, F: 1.0150).
7. **Differential-expression layer (`src/vcc2026/rank_up.py`, `src/vcc2026/dn_mag2.py`)**: Control detectability + predicted-effect call-set prior, cross-line sign consensus (`clean4:0.05`) with top-20% protection, minimum-detectable-effect floor, rank-neutral UP adjustment, and DOWN magnitude pass.
8. **Decoupled carrier/ballast count renderer + zero-gene repair (`src/vcc2026/rc_emit_k.py`, `src/vcc2026/zero_repair.py`)**: Renders 400 sparse integer-count cells per perturbation. A small set of deep ballast cells (within the 1,000,000 count/cell cap) carries the predicted pseudobulk profile, while the remaining carrier cells follow the context's control library quantiles and carry the DE ranking and call-set design. Gated genes whose rendered total would collapse to zero receive a floor repair.

## Key Findings from 650+ Experiments

The full record lives in [`docs/FINDINGS.md`](docs/FINDINGS.md) and each `experiments/<NNN>/RESULTS.md`. Four takeaways that held up across the entire competition:

- **Bulk and DE are separate functionals of the same count matrix.** In `cell-eval2`, `pds` and `mse` read log-pseudobulk of summed counts, while `fid`, `reach`, `jac`, and `nmae` read per-cell arithmetic CPM and Wilcoxon ranks. Designing those two channels separately through a carrier/ballast renderer was one of the largest levers in the entire project.
- **Donor diversity beat model complexity.** On zero-shot unseen cell lines, adding independent Perturb-seq lines and whitening cross-line agreement beat every fine-tuned foundation model or deep architecture we tested (including Arc's State baseline and PIE). Same-line shared response programmes proved lab-specific and anti-informative when transferred across lines.
- **Offline harnesses can silently lie.** Depth-confounded harnesses (see [`docs/HARNESS-DEPTH-CONFOUND.md`](docs/HARNESS-DEPTH-CONFOUND.md)), H1 call-set power mismatches, and un-repaired zero-total genes (which cost ~31 error per collapse when a local harness skips production's zero-repair) caused multiple false positives early on until we enforced production-faithful scoring and known-answer retrodiction checks.
- **Discipline compounds faster than cleverness.** Preregistering pass/fail bins before reading board scores, holding out sealed validation reserves (`R6`, `R7`, `LZ`) that could only be opened once, requiring byte-for-byte known-answer checks on production code changes, and rehearsing the entire Oct 22 release-day pipeline (`docs/DEF-RUNBOOK.md`) prevented overfitting to the 2-slot/day validation board.

## Repository Layout

| Path | Contents |
|---|---|
| `src/vcc2026/` | Core library: donor transfer, whitening, carrier/ballast emitter, UP/DOWN magnitude passes, zero-gene repair |
| `scripts/` | Production submission builder (`make_submission_pm_cd4.py`), format writer (`build_submission.py`), release-day automation (`def_day.sh`), EC2 `cell-eval2` runner (`cloud_score.sh`), dataset downloaders, and leaderboard archive tool (`reconstruct_leaderboard.py`) |
| `experiments/` | 650+ numbered experiments with their code (`.py`/`.sh`), `RESULTS.md`, `PREREGISTRATION*.md`, `BOARD-RESULT*.md`, figures, and summary tables (including `811-lolo-model/`, the final consensus blend) |
| `docs/` | Master findings (`FINDINGS.md`), method description (`FINALIST-DESCRIPTION-DRAFT.md`), frozen October switches (`october-decisions.json`), release-day runbook (`DEF-RUNBOOK.md`), confirmation and board ledgers, and vendored `cell-eval2` metric specs (`docs/reference/`) |
| `data/leaderboard/` | Lossless archive of **313 public validation-phase leaderboard snapshots** (Aug 21 - Oct 22, 2026; 24,189 unique entry versions across 2,077 teams; compressed from 415 MB to 7.3 MB) |
| `data/board-entries/` | Captured official API records for this team's validation and final entries |
| `submissions/` | `SUBMITTED_LEDGER.json` and `*.h5ad.meta.json` build manifests and hashes for every submission |
| `analysis/` | Standalone diagnostic scripts and notes |
| `AGENTS.md` | Working rules and methodology used to run the lab across coding-agent sessions |

## Historical Leaderboard Dataset

`data/leaderboard/` preserves 313 snapshots of the public VCC 2026 validation leaderboard captured between 2026-08-21 and 2026-10-22 (before the leaderboard went dark for the final phase). Because raw snapshots repeat the same entries across timestamps, they are stored losslessly in two gzipped JSONL files (`entries.jsonl.gz` + `snapshots.jsonl.gz`, 7.3 MB total).

```bash
# Summary statistics
python scripts/reconstruct_leaderboard.py --summary

# Inspect any team's validation trajectory
python scripts/reconstruct_leaderboard.py --team "Ahmet Dedeler"

# Unpack all 313 original JSON files
python scripts/reconstruct_leaderboard.py --unpack /tmp/vcc2026-leaderboard/
```

See [`data/leaderboard/README.md`](data/leaderboard/README.md) for the exact schema.

## Reproducing / Upstream Data

Raw single-cell matrices (`.h5ad`), derived donor tables (`data/donors/`), and `.vcc` submission packages are not stored in this repository. To rebuild from scratch:

```bash
uv sync
vcc datasets download controls               # official non-targeting controls
./scripts/download_2025_h1.sh                # 2025 VCC H1 dataset
```

Links and citations for the remaining public Perturb-seq screens (Replogle, Nadig, X-Atlas/Orion, Pan 2025, primary CD4+ T cells, Jiang 2025, DLD-1) are documented in [`docs/FINDINGS.md`](docs/FINDINGS.md) and [`docs/FINALIST-DESCRIPTION-DRAFT.md`](docs/FINALIST-DESCRIPTION-DRAFT.md).

## Open Science & Beneficial Use

This repository is released openly to advance transparent, reproducible science in computational biology and virtual cell modeling. Meaningful progress in perturbation biology depends on shared baselines, honest accounting of negative results, open evaluation mechanics, and reproducible experimental records, not only winning headlines.

By publishing this codebase, experimental log, and historical benchmark data, the intent is to support beneficial scientific research: understanding gene regulation, improving therapeutic target discovery, strengthening single-cell benchmarks, and helping researchers and students build better methods faster. Anyone building on this work is strongly encouraged to apply it responsibly toward public-interest science and human health, and to refrain from any harmful biological application or misuse.

## License

- **Code** (`src/`, `scripts/`, `analysis/`, `experiments/**/*.py`, `experiments/**/*.sh`): licensed under the **Apache License, Version 2.0** ([`LICENSE`](LICENSE)). Apache 2.0 is used deliberately so academic labs, independent researchers, and biotech teams can freely inspect, reuse, and build on the software with clear copyright and patent terms.
- **Documentation, research notes, figures, and leaderboard archive** (`README.md`, `docs/`, `experiments/**/*.md`, `data/leaderboard/`, `submissions/`): licensed under **Creative Commons Attribution 4.0 International (CC BY 4.0)**.
- See [`NOTICE`](NOTICE) for attribution details and upstream dataset notes.

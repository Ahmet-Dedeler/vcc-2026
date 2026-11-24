# DRAFT §3(j) — What a transfer model's DE signs miss: the recipient's own shared response programme

*Draft of 2026-09-24 (coordinator). MERGED into PAPER.md as §3(j) (F18) on 2026-09-25; PAPER.md is the live text. Sources: experiments 440, 441 (RESULTS.md, all
postscripts). No provenance verdict yet: the novelty audit (287/297 style) must run before this becomes a claim row.
Likely KNOWN in part: Replogle 2022 describes a common transcriptional response shared by many perturbations.
Candidate NEW parts: its line-specificity at the resolution of per-gene DE signs, and the negative results on
recovering it for an unseen line.*

## The measurement
On the challenge's H1 reference (h1-test: 85 targets, 400 cells each, the scorer's own Wilcoxon table), take each
target's significant genes (the set the fold-change member scores). Compare the sign of the reference lfc with the sign
of four predictors:
1. our cross-line transfer model;
2. a leave-one-target-out mean of the OTHER targets' reference responses in the same line (the recipient's shared
   programme; an oracle, because it reads the recipient's perturbed cells);
3. a cross-line consensus of four other lines' mean responses (K562 genome-wide, K562 essential, Jurkat, HepG2);
4. features of the recipient's unperturbed cells.

| predictor | significant genes we call | significant genes we do not call (≈60 %) | all |
|---|---:|---:|---:|
| transfer model | **0.649** (|l|-weighted) | 0.540 | 0.568 |
| recipient's shared programme (oracle) | 0.640 | **0.630** | **0.632** |
| cross-line consensus of mean responses | 0.571 | 0.571 | 0.573 |

**Result 1.** Where the transfer model is confident enough to call a gene, its signs match or beat the recipient's own
shared programme. On the majority of significant genes it does not call, it is near chance (0.54). The shared
programme gets 0.63 there. The DE accuracy gap between a zero-shot transfer model and an oracle sits almost entirely in
that programme.

## The programme is line-specific, and the recipient's controls do not reveal it
- Other lines' mean responses reach only 0.57–0.58. Restricting a donor line's mean to the recipient's own panel targets
  makes it worse (K562: 0.570 → 0.562). The programme is a property of the line, not of the panel.
- The recipient's unperturbed cells do not recover it:
  - top principal components of the control cells: best 0.557, |corr| with the mean response ≤ 0.18;
  - a G1-shift axis (least-cycling control cells minus all): 0.514;
  - a small-cell (low-library) axis: 0.515.
- One real recipient-side signal: the control identity axis. Genes expressed relatively higher in the recipient's
  controls than in other lines' controls tend to go down across its perturbations. The size depends on measuring it
  within one chemistry:
  - cross-chemistry (Flex H1 vs 10x-3′ donors): Spearman −0.13, sign 0.545;
  - within Flex (H1 vs the challenge's A/B/C controls): corr −0.236, cos² 5.6 %, sign 0.556–0.568 on every gene;
  - on top of the cross-line consensus it lifts significant block genes 0.580 → 0.592.
  This is the mechanism the July 2026 decomposition preprint reports (template along the identity axis; xinyizhanglab,
  fig. 1 supp.). Here it is measured on a non-essential Flex panel, where it is real but explains only ~6 % of the
  template.

## …but a screen of the SAME line from another lab does carry it (added 2026-09-25, experiment 450)
Pan et al. 2025 (Whitehead; GSE295214) ran a genome-scale CRISPRi Perturb-seq screen in H1 hESCs. It differs from the
2025 VCC H1 truth in lab, library and chemistry (10x 5′ vs Flex). Its programme (per-gene fraction of perturbations
going down, h1-test's targets excluded), on genes both priors cover:

| predictor | gate-block | all gate | fid-like block |
|---|---:|---:|---:|
| cross-line consensus (CLEAN4) | 0.581 | 0.574 | 0.524 |
| **same line, other lab (Pan)** | **0.625** | **0.594** | **0.539** |
| recipient's own programme (oracle) | 0.653 | 0.638 | 0.555 |

- The other-lab screen closes 62 % of the gap on the block genes.
- Per target: +0.057 [95 % CI +0.024, +0.091], better on 37 / 49 targets.
- Its essential-gene perturbations alone carry the programme (0.627). So it is a property of the line, stable across
  lab, chemistry and target set, not of any particular perturbation set.
- The same screen's own rows for the 30 h1-test targets it covers go further:
  - gate signs 0.764 vs 0.551 for the cross-line transfer model (26 / 26 targets);
  - target-centred direction 0.260 vs 0.069 (30 / 30).
- Replication on a cancer line: with VIPerturb-seq K562 (Satija lab) as the truth, Replogle's K562 screen beats the
  best cross-line bank 2.4× in direction on essential targets.

**Through the full prediction pipeline** (the challenge's own scorer on the H1 surface, the production emitter), blending
the other lab's rows into our cross-line prediction for the 30 covered targets (35 % of the panel) raises the score by
**+0.08** (board units): mse −0.17, nmae −0.07, reach +0.07, pds +0.011, fid +0.03. That is larger than any single
change in our record. Two negatives bound the claim:
- a screen of a different line of the same tissue (DLD-1 for HCT116) does not beat cross-line banks;
- a same-line L1000 shRNA signature set (CMap HepG2 for Nadig HepG2) is worse than other lines' L1000 signatures.
The benefit needs the same line and a Perturb-seq-grade measurement.

**Result 2.** The line-specific shared programme is unobservable from the recipient's unperturbed cells, but it
transfers from any other perturbation screen of the same line. So a "virtual cell" for an unseen line is bounded by
whether anyone has perturbed that line, not by who did it or how.

## A benchmark fact recovered from public submissions
The fold-change set member (Jaccard) scores 1.0 when both sides call nothing. About fifty public control-copy entries,
which call nothing on any target, all read exactly 0 on the validation board. So every one of the 300 targets has at
least one significant non-target gene in each of the three validation lines. Arc's selection ("a strong set of
perturbations and responses") holds at the level of each (target, line) pair.

## Why it matters
- Zero-shot prediction in an unseen line is bounded less by target-specific effects, which transfer models already call
  well, than by the line's generic response to being perturbed. No public data we hold estimates that for a line never
  perturbed.
- A line-specific shared programme is the cheapest thing a new screen could measure: tens of random perturbations would
  do. It is also the most valuable single covariate for transfer.
- The leaderboard leader appears to train on a 100M-cell, many-line proprietary atlas (Illumina PIP-seq). That is
  exactly the resource that supplies the oracle row above.

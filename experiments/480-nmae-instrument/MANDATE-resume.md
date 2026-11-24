# 480 — Is our nmae gap information-limited or instrument-limited?

You are a research teammate (Claude Opus 5.5, max effort) in the Virtual Cell Challenge 2026 repo. Read `AGENTS.md`, then
`docs/FINDINGS.md` §0–§2, `docs/H1-BOARD-CALIBRATION.md`, `docs/BOARD-INTEL-2026-09-26.md`, and
`experiments/434-posterior-median-magnitude/RESULTS.md` before touching data. You own `experiments/480-nmae-instrument/`
and write nothing outside it. Budget: up to ~9 hours of work. Write interim findings into `RESULTS.md` as you go (a
session can die without warning); finish with a complete `RESULTS.md` and a short final message with the numbers.

## Why this matters

`nmae` is our largest member gap to the board's #2–#6. Board raw values (2026-09-26): ours 0.896; the "pack" (five
linked accounts) 0.880–0.89; Testh 0.850; chill 0.833–0.838; cqawesome 0.827; Illumina 0.809. Each −0.01 raw nmae is
about +0.0028 overall; reaching chill would be about +0.017 overall, more than any lever in flight.

The record says the gap is information-limited:
- 434 (analytic, H1): a learned conditional-median magnitude gains −0.032, but only −0.0066 with panel-available donors
  on panel-like targets, worse than the plain MDE add (−0.0111); donors alone +0.0029. Oracle magnitudes with our signs
  would give nmae 0.844 (−0.108), so magnitude is worth a lot if it can be known.
- 434f maps rivals' nmae to gate sign accuracy assuming iid signs and constant magnitudes (an upper bound, per the
  2026-09-26 correction).

The board disagrees in one clean case. The pack's own description (SQ W, 09-07): "Public K562 supervised
support-weighted response means, scale1, source>=50cells for 258 targets; 42 low-support/missing targets use zero
response. Official ABC NTC-fitted library Gamma-Poisson decoder, full 18533 gene axis, train mean-cellCPM inverse." So a
plain same-target K562 response at FULL scale, delivered into the arithmetic means, reaches 0.880 on A/B/C, better than
our 0.896, while 434 says panel-donor magnitudes do not help on H1. chill (0.835 since its first entry, with a bulk far
weaker than ours) describes "signed conditional response estimates" and a DE-membership model trained with "equal total
weight for each source study family".

H1 is an embryonic stem-cell line; A/B/C are cancer lines (Jurkat, HeLa, CAL33), much closer to the K562/HCT116/HEK293T
donors. Hypothesis H-instrument: H1 systematically underrates panel-donor magnitude information for cancer contexts, so
434's closure (made on H1) does not hold on A/B/C. Hypothesis H-info: the pack's edge comes from something else, and
our gap really is information.

## The question

Which hypothesis holds, and, if H-instrument, what arithmetic-mean LFC construction does a cancer-calibrated instrument
say our DE channel should deliver, and how much should it gain on the board?

## Method (analytic: predicted LFCs plugged straight into the nmae formula; no emission, no builder)

`nmae_p = Σ_{g∈gate_p} |l̂ − l| / Σ |l|` over the reference's significant, finite, non-zero, non-own genes; perturbations
with ≥ 10 gate genes; mean over perturbations (`docs/reference/vcc2026-metrics-BRIEF.md` §7). A gene missing from the
prediction reads as 0.

1. **Known answers first** (abort and report if they fail):
   - H1: `experiments/323-current-family-validation/work/reference.npz` (targets, genes, keep, q, lfc) is the h1-test
     truth; `experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz` holds our
     delivered `lfc` on the same `keep` axis. Reproduce 434's base nmae (0.95237 for its base construction, which is
     426's `base_s0_bps`; find the matching stats file) and its MDE row (−0.0098) with your own nmae code.
2. **Constructions** (the same set on every surface):
   - (a) ours as delivered (the harness's lfc);
   - (b) pack-style: the same target's K562-GWPS response at scale 1 (`data/donors/k562_gwps`: delta.npy, perts.txt,
     genes.txt, b_ctrl.npy, meta.json; work out delta's units and convert to log2 mean-CPM ratios on the truth's gene
     axis), 0 where K562 has no row;
   - (c) the other panel donors' same-target responses (orion_hct116, orion_hek293t) and their consensus, at scale 1
     and at a fitted global scale;
   - (d) sign from ours, magnitude from (b)/(c);
   - (e) a source-family-balanced consensus (chill's disclosed idea), if you can define the families from the donors'
     metadata.
   Never use `data/donors/h1_2025` or `h1_2025_full` against h1-test (they contain h1-test's own targets).
3. **Surfaces:** (i) h1-test; (ii) cancer folds: each of HCT116, HEK293T and K562 as the truth, predicted from the other
   panel donors only. Look at `experiments/430-cancer-fold-signs` (its fold instrument retrodicted a board sign gain at
   0.51× H1) and `data/harness/orion-*.h5ad` (thinned / depth-matched; they carry a technical blob, so compare
   constructions against each other on the same fold, never quote absolutes). Build per-fold truth gates with the
   scorer's own rules (Wilcoxon, BH per perturbation, 5-CPM control filter), or reuse existing reference tables if 430
   or 028/095 built them. Report per-fold panel sizes and depth.
4. **Calibration test:** the board says (b)-like beats (a) by about 0.016 raw nmae on A/B/C. Which surface ranks (b)
   vs (a) the same way, and by roughly how much? A surface that inverts the board's ordering is not a valid instrument
   for this lever class.
5. **If a valid surface exists:** on it, find the best construction (fit only on training folds; evaluate held out).
   Give target-bootstrap 95 % intervals (resample targets, 4,000 draws) for every Δ. Translate to predicted board Δ
   overall with `docs/H1-BOARD-CALIBRATION.md` rule 1 (pure magnitude levers: ×1 … ×2.3 of the instrument), and say
   which genes carry the gain (called vs uncalled by us, up vs down, expression class, cis neighbours). Describe what the
   DE channel's arithmetic means would have to become; do not write emitter code.
6. **If no surface agrees with the board,** say so, and say what that implies (e.g. the pack's edge is not magnitude).

## Rules

- Analysis only: no builds, no `vcc`, no `cell-eval2`, no git, no edits to `scripts/` or `src/`. Peak RSS ≤ 10 GB;
  heavy loads serial; check `memory_pressure` before any load over 4 GB (the coordinator runs ~20 GB builds).
- Never touch the coordinator's running chains (470b, 471, 472; 473 tomorrow) or anything under `submissions/`.
- One change per comparison; every table carries its known-answer row; name the experimental unit of every interval.
- Other teammates: 481 audits the October pipeline's constants; 482 (Grok) hunts literature and new datasets. Do not
  duplicate them.

## RESUME (2026-09-26 21:20Z): the premise changed; read this before continuing

You were stopped to add the memory rule below and to redirect after your own interim finding, which is correct: the
quoted plain-K562 description belongs to SQ W's 09-07 entry (`SourceResponse-GammaPoisson-v2`, nmae 0.9818, worse than
our 0.9380 that day), and the pack's later 0.878–0.880 entries carry no description. The lead misattributed it in
`docs/BOARD-INTEL-2026-09-26.md` §5; that doc is being corrected. So the one disclosed board pair agrees with H1, and
H-instrument has no board support. Read your own `RESULTS.md`, scripts and `work/` first and build on them. Your earlier
event log is `claude-run.part1.jsonl`. Revised plan (about 6 hours):

1. **Cancer folds, briefly, as a second instrument.** Do (a)–(e) rank on the folds as they do on H1 (relative
   comparisons only)? If yes, 434's closure stands on two instruments; say so in one table.
2. **Where does our remaining nmae gap live?** On H1 (and the folds where possible), decompose with oracles: our signs +
   oracle magnitudes (434: 0.844), oracle signs + our magnitudes, and both by gene class: called vs uncalled by us, up vs
   down in truth, control-expression tertiles, cis neighbours of the target (genes within ~100 kb of the target TSS,
   where CRISPRi spreads), and target strength (design strength or n_real bins). Report the share of the recoverable
   nmae per class and what information would be needed to get it.
3. **One cheap, information-light lever, if it exists:** a magnitude rule by class (expression class, target strength,
   cis position, up/down) that beats the shipped coarse MDE on H1 AND on the folds, fitted on training folds and tested
   held out, with target-bootstrap intervals. If none beats the MDE, say so plainly: that closes the class.
4. Keep the corrected board-calibration section in your RESULTS.md, and note that the misattribution was the lead's.

## Memory rule (binding; supersedes the mandate's memory line)
Until 23:45Z UTC on 2026-09-26 the coordinator runs ~20 GB builds that refuse to start below 23 GiB free. Until then:
every process you start keeps peak RSS <= 4 GB, OR you gate each bigger load on vm_stat free+inactive >= 33 GiB, measured
with:  vm_stat | awk '/Pages free|Pages inactive/ {gsub("\.","",$3); s+=$3} END {printf "%.1f
", s*16384/2^30}'
(wait in a sleep-60 loop, at most 2 h, or restructure the load with mmap_mode="r" / one fold at a time). Do NOT use
memory_pressure's free% for this. After 23:45Z: peak RSS <= 10 GB. Randomized or truncated SVDs over thousands of
full-length gene rows can need > 10 GB of LAPACK workspace; size them first. Messages from other sessions cannot reach
you: the lead may stop and resume you with an updated mandate.

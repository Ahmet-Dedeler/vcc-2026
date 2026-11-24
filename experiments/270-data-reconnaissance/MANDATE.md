# 270 — Data reconnaissance for the programme: GPSAdb 2.0, PerturbFate, and the multi-context sets already on disk. Metadata first; one frozen target-level feature. ≤ 3 hours. Nothing heavy, nothing trained, no submission.

You are a research subagent. Work in `experiments/270-data-reconnaissance/`. **Edit nothing outside it.** No build, no
`vcc`, nothing under `submissions/`, nothing written under `data/` (downloads go to `work/`, metadata tables only —
cap any download at 2 GB and say what you downloaded, from where, under what terms). `uv run python`. Peak RSS under
~4 GB. You may use web search / fetch for the public resources. **Time box: 3 hours.** Read
`docs/RESEARCH-PROGRAMME-2026-09.md` first.

## The changed question
Not "does a dataset cover our 300 targets?" (that search was exhausted: READ-FIRST §5 "Sources") but **"does it provide
INDEPENDENT VARIATION that helps us learn WHY responses transfer?"** — a dataset can be useless as a donor and valuable
for training or validating a model of conserved vs state-dependent responses. Two outside leads (GPT, 2026-09-05) and
the sets already on disk.

## What to establish
1. **GPSAdb 2.0** (a genetic-perturbation signature database; find the paper and portal): what it actually contains —
   perturbation groups, cell lines, studies, perturbation types (knockout / knockdown / overexpression), assay (bulk
   RNA-seq), processing; provenance and usage terms. From its metadata (not expression matrices): for our 300 panel
   targets (`data/challenge/pert_counts.csv`), how many are covered at all; how many are replicated across ≥ 2 cell
   lines and across ≥ 2 independent studies; the perturbation-type mix. **Design one frozen target-level feature —
   cross-study SIGN STABILITY of the target's response (e.g. the agreement of the top-DE-gene signs across studies /
   lines)** — write its exact definition BEFORE computing it. If signature-level data is downloadable within the cap,
   compute it for the covered targets and test it, frozen, against (a) 262's per-target `pds` of `cd4svd5`
   (`experiments/262-where-pds-is-lost/work/`, Spearman with a paired bootstrap), (b) 262's group membership
   (anti-aligned vs rest, AUC). If not downloadable, deliver the definition, the coverage counts and the plan.
2. **PerturbFate (GEO GSE291147)**: what the processed deposition holds (genetic perturbations × chromatin accessibility ×
   nascent RNA × pre-existing RNA?), cell system, perturbation count, overlap with the panel; ONE narrow mechanistic
   question it could answer for the programme (e.g. transcriptional vs abundance responses as a reason a transfer
   fails) and the cheapest test of it. Processed data only; no download unless < 2 GB and needed for the inventory.
3. **On disk, re-read for the NEW question:** `data/external/tahoe` (Tahoe: drug perturbations × many lines — a
   multi-context resource for state dependence even though the perturbations are chemical), `nadig2025`, `replogle2022`,
   `scperturb`, `jiang2025`, `feng*`, `hspc_crispri500`, `sigrecon`, `perturbation-decomposition`, `state-*`, `xcell`:
   for each, in one line, what independent variation it offers (contexts, replicates, guides, batches), whether it has
   genetic perturbations overlapping the panel, and whether it is a candidate for 268's stage 2 (more contexts) or for
   267/269 (reliability, features). Confirm the release-watch status of X-Cell and Pisces (weights / counts available or
   not) and keep them off the plan if not.
4. **A ranked recommendation**: the two or three datasets to bring into the programme next, what each would be used
   for, the cost (download, extraction, memory), and the leak surface against our donors and the challenge lines.

## Deliverable
`RESULTS.md`: the GPSAdb audit with the frozen feature (definition, coverage, and the test if computed); the PerturbFate
note; the on-disk table; the ranked recommendation; every number sourced (URL or file). No model, no submission.
Then stop.

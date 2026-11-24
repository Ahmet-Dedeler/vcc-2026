# 279 — A preregistered biological claim from the conservation map, tested on HELD targets: are KRAB-ZNF / zinc-finger perturbations less conserved across cell lines than other perturbations? ≤ 4 hours. Measurement only. No build, no submission.

You are a research subagent. Work in `experiments/279-znf-state-dependence/`. **Edit nothing outside it.** No build, no
`vcc`, no `cell-eval2`, nothing under `submissions/` or `data/`. `uv run python`. One raw-cell extraction at a time,
peak RSS under ~8 GB. **Time box: 4 hours.** `docs/RESEARCH-PROGRAMME-2026-09.md` binds you: PREREGISTER before you
look, held targets never seen when the hypothesis was formed, refit inside the split, paired bootstraps, the unit named.

## Why (R276 §5, §6.3b; F8)
276's conservation map found, POST HOC, that the zinc-finger family (ZNF*) sits in the low-conservation tail with a
raw effect at p ≈ 1e-5 and named it as a claim to preregister and test on held targets: **perturbations of KRAB-ZNF /
C2H2 zinc-finger transcription factors produce responses that are less conserved across cell lines (more
state-dependent) than perturbations of other genes.** The claim is biologically plausible (many KRAB-ZNFs are
lineage-restricted repressors of transposable elements with cell-type-specific targets) and it is testable on targets
the panel never contained: the same datasets on disk carry hundreds of non-panel ZNF and non-ZNF perturbations.

## READ FIRST
`docs/RESEARCH-PROGRAMME-2026-09.md`; 276's `RESULTS.md` §2 (the frozen statistic — REUSE `build276.py`'s conservation
computation exactly; do not redefine it), §5 (the post-hoc ZNF observation: which targets, which effect size, how it
was found — the hypothesis-forming set), `work/DEFINITION-276.md`, `work/conservation276.csv`; 267's half-banks and 271's
`halfrows271.py` (for reliability correction where halves exist); the datasets with many non-panel genetic
perturbations across ≥ 2 lines on disk (Replogle K562 GWPS + RPE1 essential/GWPS, Nadig HepG2/Jurkat, Broad-16's 90,
Jiang's 218 — check which carry ZNF genes in ≥ 2 lines; `scripts/build_donor_table.py` for the response object).

## Protocol (write `PREREGISTRATION.md` with a sha256 BEFORE computing any held-target conservation)
1. **Hypothesis-forming set (exclude):** every panel target (300) — the map that produced the observation. **Held test
   set:** all NON-panel perturbed genes measured in ≥ 2 lines on disk, split into ZNF (HGNC zinc-finger families: C2H2
   ZNF*, KRAB-domain; take the family list from `data/external/omnipath-dorothea` / `regulon` / a gene-family table on
   disk or HGNC's public family page — cite it; decide the membership rule before looking at any response) and non-ZNF.
2. **Prediction, frozen:** the ZNF held targets' conservation (276's statistic, raw and, where halves exist,
   reliability-corrected) is lower than the non-ZNF held targets' — one-sided; effect size = the difference in medians and
   the AUC; pass rule = the paired-bootstrap 95 % interval of the median difference excludes 0 in the predicted direction
   AND the effect survives matching on expression level and response magnitude (ZNFs are lowly expressed and their
   responses are weak: state the confounds and the matching rule in advance — e.g. within-decile of CCLE median
   expression and of response norm).
3. **Controls:** (a) a random gene family of the same size drawn from the held set (100 draws: the null distribution of
   the median difference); (b) the OTHER large TF family on disk (e.g. homeobox or bHLH) — is the effect ZNF-specific or
   "any TF family"; (c) the same test within a single dataset (Replogle K562 vs RPE1) to show it is not a dataset mix.
4. **Then, and only then, compute.** Report the held-set sizes, the raw and matched effect, the AUC, the null band, the
   two controls; and the same numbers recomputed on the panel targets (the forming set) for reference, clearly labelled
   as in-sample.
5. **Mechanistic descriptor (post-hoc, labelled):** for the ZNF held targets, is low conservation associated with low
   response reliability (noise) or with reproducible line-specific direction (true state dependence)? Use the split-half
   reliability where halves exist to split the two.

## Deliverable
`PREREGISTRATION.md` (hashed, timestamped, before any held number), `RESULTS.md`: the held-set inventory, the frozen
prediction and pass rule, the result with intervals, the null band and the two controls, the in-sample reference, the
noise-vs-state-dependence descriptor, and one line: **CONFIRMED / NOT CONFIRMED on held targets (median difference …,
AUC …)**. Write the per-gene table to `work/znf279.csv`. No model, no submission. Then stop.

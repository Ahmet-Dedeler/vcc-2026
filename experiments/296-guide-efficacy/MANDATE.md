# 296 — Is F14's "construct term" guide EFFICACY? The knockdown stratification 293 found in the challenge's library, applied to 285's public libraries (K562-GWPS 748 two-construct perturbations; Nadig Jurkat 124; HepG2 70; Broad-16 Cas9 97 × 12 lines; the CD4 bank's guide pairs), plus the "failing construct still gives a reproducible transcriptome" phenomenon and each library's construct-failure rate. Preregistered, ≤ 4 hours. Analysis only. No build, no submission.

You are a research subagent. Work in `experiments/296-guide-efficacy/`. **Edit nothing outside it.** No build, no
`vcc`, no `cell-eval2`, nothing under `submissions/`, `data/harness/` or `data/`; no download. `uv run python`,
streaming / row-blocked reads; peak RSS under ~10 GB, one raw pass at a time, `memory_pressure` ≥ 40 % free before
any pass ≥ 8 GiB. **Time box: 4 hours.** One sibling (294, the write-up) runs in parallel; do not read or write its
directory. Both old reserves and R3 are SPENT; none is needed here.

## Why (R293 §5, §7, §10; F14 in `docs/FINDINGS-2026-09.md` v2.0)
285 (F14) measured the guide construct as the dominant irreducible term of a perturbation response: two constructs of
one gene agree at ρ_guide 0.09 (K562-GWPS, 748 pairs), 0.17 (Jurkat), 0.41 (HepG2), 0.02–0.06 (CD4), 0.76–0.81
(Broad-16, Cas9) of the same-construct cell split. 293 then measured Arc's own 2025 library and found the same low
number (0.075) — but split by the WEAKER construct's on-target knockdown it is bimodal: **+0.687 [+0.197, +1.067]**
where both constructs reach |lfc| ≥ 1.0, **−0.003 [−0.085, +0.181]** where the weaker does not; the knockdown of the
weaker construct predicts the ratio at Spearman +0.63; two-construct targets are the library's weakest; and a construct
that misses its target still produces a REPRODUCIBLE transcriptome (split-half +0.27 and +0.05 for two failing
constructs whose working partners read +0.88 and +0.82 and point elsewhere). 285 §6b's GWPS terciles reached only
ρ ≈ 0.24 at the strongest knockdown, with a tercile rule, at n = 25, and with knockdown predicting the cell split as
strongly as the guide split. **If the public-library construct term is also mostly knockdown failure, F14's headline
becomes "guide efficacy is the dominant irreducible term" — a different and more actionable claim (a validated
construct transfers; a failing one is a reproducible response to the vector). If it is not, the construct term is
construct-specific biology and F14 stands as written.**

## READ FIRST
285's `RESULTS.md` §5–§6b, §11 and `work/lib285.py` (the response object, the matched-n ladder, the two readouts, the ρ
construction), `work/k562_guide.{py,npz,json}` (the 748 pairs, both readouts), `work/k562_kd.{py,npz,json}` (the
per-construct on-target lfc and the tercile stratification — the starting point), `work/nadig_285.py`,
`work/broad_285.py`, `work/cd4_285.py`, `work/csplit_285.py`; 293's `RESULTS.md` §3–§6 and `work/{lib293.py, run293.py,
strat293.py, ceiling293.py, construct293.csv}` (the exact stratification rule, the (a2) column, the cross-target null
(x), the disattenuated ρ_dis, the variance model). Data: `data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad`
(guide identities; 775 two-construct perturbations), `data/external/nadig2025/*_raw_singlecell.h5ad` (`sgID_AB`),
268's `work/halves/<line>.npz` (Broad-16 `v1` / `v2` guide halves), 280's bank (`bank280_X.npy`, `obs280.npz`: guides
per (target, donor, condition)); the on-target column of each construct's own gene where the gene is on the axis.

## Protocol (write `PREREGISTRATION.md`, hashed into `work/PREREGISTRATION.sha256`, before any stratified number)
1. **Known answers:** 285's K562 guide row (748 pairs, n = 25: (b) +0.0006, (a) +0.0063, `pds` 0.633 → 0.509) and its
   §6b tercile table at 0.0e+00 from its saved arrays; 293's ALL / both-≥1.0 / weaker-<1.0 rows (+0.075 / +0.687 /
   −0.003) re-read from `construct293.csv`; the union identity of halves.
2. **Frozen rule, every library:** the weaker construct's on-target knockdown |lfc| (its own gene's column, full
   depth, excluded from every scored statistic) with 293's cut at 1.0 AND a library-relative cut (the library's
   median single-construct |lfc|); strata: both constructs ≥ cut / weaker < cut; at the largest matched n each
   library supports (K562 n = 25 and 50; Jurkat / HepG2 25; Broad-16 at its UMI budgets; CD4 at its budgets).
   Report per stratum: (a) same-construct cell split, (a2) the other construct's own split, (b) between constructs,
   (x) the cross-target null (constructs of different genes), ρ_guide = (b)/(a), ρ_dis = (b)/√((a)(a2)), the scorer's
   `pds` (a) / (b) / (x); paired target bootstrap. **Frozen predictions:** (P1) in every CRISPRi library ρ_guide in the
   both-≥1.0 stratum exceeds the weaker-<1.0 stratum with the paired interval above 0; (P2) the both-≥1.0 ρ_guide is
   ≥ 0.4 in K562-GWPS (293's +0.687 minus its width); (P3) Broad-16 (Cas9, validated guides) shows no such split
   (its constructs mostly work); (P4) knockdown of the weaker construct predicts ρ_dis at Spearman > +0.3 in each
   CRISPRi library. State each outcome.
3. **The failing-construct phenomenon:** among constructs with |lfc| < 0.5 (no knockdown), the same-construct
   split-half reliability vs that of working constructs (≥ 1.0), per library; the cosine of a failing construct's
   response to its working partner's; to the library's mean response (is it a generic stress / vector signature?);
   to other failing constructs' responses (a shared vector response would show as agreement ACROSS failing
   constructs of different genes — the (x) null restricted to failing pairs). Frozen prediction (P5): failing
   constructs' split-half reliability > 0 with the interval above 0 (reproducible), and their cross-construct
   agreement exceeds the working constructs' (x) null (a shared vector signature) — state both outcomes.
4. **Library failure rates:** the share of single-construct targets with |lfc| < 1.0 and < 0.5 per library (K562-GWPS
   panel and non-panel separately; Nadig; Broad-16; CD4 from `guide_kd_efficiency.suppl_table.csv` if on disk under
   `data/external/cd4_cell2026/suppl_tables_2026-08-24/`), next to 293's 18.7 % for the 2025 panel — the "one truth in
   five" estimate's range across libraries; and, for the 2026 challenge, the fraction of the 300 panel targets whose
   K562-GWPS construct fails (an indicative number only; the 2026 constructs are different).
5. **Descriptive (post-hoc, labelled):** whether promoter choice or transcript-level constructs (293's KDM2B) explain
   any of the split; the 20 most and least efficient constructs per library.

## Deliverable
`PREREGISTRATION.md` (hashed), `RESULTS.md`: KA rows; the stratified table per library (both cuts); P1–P5 with their
rules; the failing-construct phenomenon; the failure-rate table; the descriptive tail; `work/efficacy296.csv` (per
(library, target, construct): cells, |lfc|, (a), (a2), (b), (x), ρ_guide, ρ_dis, `pds`); one line: **GUIDE EFFICACY
EXPLAINS THE CONSTRUCT TERM: YES / PARTLY / NO (K562 both-≥1.0 ρ_guide … vs weaker-<1.0 …; Jurkat …; HepG2 …;
Broad-16 …; CD4 …); FAILING CONSTRUCTS REPRODUCIBLE: YES / NO (split-half … ; shared vector signature: YES / NO);
FAILURE RATES: …**. No model, no build, no submission. Then stop.

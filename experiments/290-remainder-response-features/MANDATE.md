# 290 — What carries the reproducible per-gene conservation remainder that annotation cannot? RESPONSE-DERIVED features measured on one half of the source set, predicting the other half's remainder — the organising question asked of the responses themselves. Preregistered, ≤ 4 hours. Analysis only. No build, no submission.

You are a research subagent. Work in `experiments/290-remainder-response-features/`. **Edit nothing outside it.** No
build, no `vcc`, no `cell-eval2`, nothing under `submissions/` or `data/`; no download (everything is on disk). `uv run
python`. Peak RSS under ~8 GB (the banks are read row-blocked / memmapped as 279 reads them). **Time box: 4 hours.**
Three sibling experiments (285, 286, 289) run in parallel; do not read or write their directories.
`docs/RESEARCH-PROGRAMME-2026-09.md` binds you (preregister before any held number; name the unit; strong baselines;
rule 10: state what is new against 284 and 281). Both reserves are SPENT: never read `experiments/267-*/work/reserve/`
or the Broad-16 lines JVE127 / MG63 / SLR23 / UMRC3.

## Why (R284 §5–§8, FINDINGS v1.3 F12)
284 showed that after 21 build-time annotation features have predicted a gene's cross-line conservation, the gene-level
REMAINDER reproduces across two halves of the 28 sources that share no source and no cell line (ρ +0.328 [+0.241,
+0.373] vs a ceiling of +0.382; 20/20 partitions; +0.277 after measured reliability and response magnitude are also
removed; +0.322 on genes whose own knockdown cannot contribute). So per-gene conservation biology exists beyond the
family / essentiality × expression axis, and annotation does not carry it. Its described tail (mRNA-surveillance and
exon-junction machinery conserved MORE than essentiality predicts; proteostasis and chromosome-segregation machinery
LESS) suggests a hypothesis: **a response is conserved when it engages gene-specific downstream programmes, and
state-dependent when it engages the recipient's proliferation / stress axis.** The clean test uses the responses
themselves as features — measured on ONE half of the sources — to predict the OTHER half's remainder, so nothing on the
target side is seen by the features.

## READ FIRST
284's `RESULTS.md` §5–§6, `rem284.py` (the admissible partitions; the primary partition #1: half 1 = DLD-1, HCT116,
HEK293T, HT-29, K562, MCF7, RPE1, THP-1, iPSC-Tian19, melanoma; half 2 = A549, BxPC-3, Calu-3, H1, HAP1, HEK293, HepG2,
Jurkat, TeloHAEC, cardiomyocyte, iPSC-neuron), `work/rem284.npz` (keys `C1, C2, r1, r2, m1, m2, ok, n1, n2, gene, fold,
comp`: the two halves' maps, remainders, magnitudes, pair counts, components) and `work/pred284.csv` (`C_half1`,
`C_half2`, `remainder_half1`, `remainder_half2`, `component`, `gene_group_all`, `r_bar`, the S1 annotation predictions);
279's `build279.py` (`residual_rows`, `load_unit`, `pair_cos`: how each source's residual rows are read on the common
2,739-gene axis; `work/held279.npz` `keys` / `lines` / `pairs_cross` / `Ccross`); 276's `work/inventory276_sources.csv`
(the 28 sources: path, dataset, line, assay); 268's `RESULTS.md` §2 (the shared / line-specific / noise decomposition —
the programme-engagement idea) and `data/external/broad16line/metadata/hallmark_symbols.gmt` (hallmark sets:
E2F_TARGETS, G2M_CHECKPOINT, P53_PATHWAY, UNFOLDED_PROTEIN_RESPONSE, MTORC1, MYC_TARGETS, INTERFERON, …); 281's feature
table `experiments/281-…/work/feat281.csv` (the annotation baseline features).

## Protocol (write `PREREGISTRATION.md`, hashed into `work/PREREGISTRATION.sha256`, BEFORE any feature–remainder number)
1. **Known answers:** rebuild 284's primary partition maps C⁽¹⁾ / C⁽²⁾ and remainders from `rem284.npz` / `pred284.csv`
   and reproduce the ceiling +0.3816 and the remainder ρ +0.3277 at 0.0e+00 (same genes, same statistic); reproduce
   279's `c_held` over all pairs at 0.0e+00 from the banks (279's KA construction) to prove the row reader.
2. **Response-derived features, computed on HALF 1's sources only** (per gene, aggregated over the half-1 sources that
   measure it; every feature has a stated formula and a stated aggregation):
   (a) **response magnitude** — mean residual row norm; (b) **breadth** — number of genes with |residual| above a fixed
   quantile of the bank's control-noise (a DE-count proxy on the residual rows; state the rule); (c) **genericness** —
   cosine of the gene's residual row to the bank's mean residual row over all its perturbations, and the fraction of
   the row's energy inside the bank's top-5 / top-20 PCs (fitted on that bank's rows); (d) **programme engagement** —
   the row's projection (signed, normalised) onto each hallmark set's indicator direction for a frozen list of ~10
   programmes (proliferation: E2F_TARGETS, G2M_CHECKPOINT, MYC_TARGETS_V1; stress: P53_PATHWAY, UNFOLDED_PROTEIN_
   RESPONSE, APOPTOSIS; metabolic: MTORC1_SIGNALING, OXIDATIVE_PHOSPHORYLATION; immune: INTERFERON_GAMMA_RESPONSE;
   plus the |projection| onto the proliferation axis as one summary); (e) **within-half consistency** — the mean
   pairwise cosine of the gene's rows across half-1 sources (this is half 1's own conservation C⁽¹⁾; it is the
   STRONGEST baseline: conservation in one half predicting the other half's remainder is the ceiling-like reference
   and must be reported as such, not as a feature); (f) **on-target knockdown** where the gene is on the axis (else
   missing, with an indicator). Then the SAME features from HALF 2's sources, for the symmetric direction.
3. **Targets:** y = half 2's REMAINDER (`remainder_half2`, i.e. C⁽²⁾ minus the annotation model's out-of-fold
   prediction) — primary; y = C⁽²⁾ raw — secondary; and the symmetric direction (half-2 features → half-1 remainder).
   Unit: genes with ≥ 1 pair in both halves (2,848 in the primary partition); component block bootstrap (284's
   `component`) AND 281's first-family block bootstrap, both reported; a second and third admissible partition from
   284's list as replication (declare which before looking).
4. **Models and splits:** ridge and a shallow GBM, hyperparameters by inner CV inside the training genes; grouped
   5-fold CV over components (as 284) for the out-of-fold prediction; per-feature Spearman with intervals.
   **Baselines:** the annotation model (284's S1 ridge, refit on the same genes; and its prediction as a covariate);
   reliability `r_bar`; magnitude alone; C⁽¹⁾ alone (the within-half conservation, the reference ceiling for
   "conservation measured elsewhere predicts conservation here"); mean predictor. **Controls:** genes shuffled
   (200 shuffles, the 281 KA-4′ form); features computed on the SAME half as the target (the in-sample leakage
   reading, labelled, to show the size of the leak a careless design would have); features with the programme
   projections replaced by random gene sets of matched size (20 draws).
5. **Frozen predictions and pass rules:** (P1) response-derived features (a)–(d), (f) predict the held half's
   remainder at out-of-fold Spearman with the component interval excluding 0 and ≥ 0.10; (P2) they beat the annotation
   model on the SAME target (paired Δρ interval above 0); (P3) the proliferation / stress engagement summary is
   NEGATIVELY associated with the remainder (Spearman < 0, interval excluding 0) — responses on the recipient's
   proliferation / stress axis are the state-dependent ones; (P4) genericness is negatively associated with the
   remainder. State each rule's outcome; the experiment's one line reports P1–P4.
6. **Descriptive (post-hoc, labelled):** attribution (permutation importance on held genes); the 20 genes best and
   worst predicted; whether 284's two tails (NMD / EJC machinery above, proteostasis / condensin below) are explained
   by (d) — which programmes each tail's responses engage.

## Deliverable
`PREREGISTRATION.md` (hashed), `RESULTS.md`: the KA table; the feature table with formulas and coverage; the
primary and symmetric results with both bootstraps and the replication partitions; the baselines and the three
controls; P1–P4 with their rules; attribution and the tails; `work/feat290.csv` (per gene: the half-1 and half-2
features), `work/pred290.csv` (per gene: targets, out-of-fold predictions, split); one line: **RESPONSE-DERIVED
FEATURES PREDICT THE HELD-HALF REMAINDER: YES / NO (ρ … [lo, hi]; annotation …; C⁽¹⁾ alone …); BEAT ANNOTATION:
YES / NO; PROLIFERATION / STRESS ENGAGEMENT → LESS CONSERVED: YES / NO (ρ …); GENERICNESS → LESS CONSERVED: YES / NO**.
No model use for the competition is proposed. Then stop.

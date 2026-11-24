# 332 — pre-outcome depth-analysis addendum

September 9, 2026, approximately 13:06 UTC. Main freeze **87adc5db**.
The main fits have started; **no new biological prediction accuracy,
arm contrast or terminal result has been inspected**. Known-answer baseline
identities and engineering timings only have been read.

The initial protocol omitted research-programme rule 9's explicit
depth-conditioned pass check. This addendum adds that check before looking
at outcomes, without changing the frozen code, data, model, fit or original
gates. Report it as a pre-outcome amendment, not part of the original freeze.

Fixed metadata: `data/external/depmap26q1/extracted/processed/cell_quality_covariates.csv`,
SHA-256 `1cf07cb227aa18821a6c08942bbcaaeead337a568e49081fa8c644cb27f7b6bf`.
Use `total_umis_per_ko`, not the misleadingly named `avg_umis_per_ko`
(the latter is on the per-cell depth scale). Metadata-only selection:

- Lower-depth six: C4I, SNU8, KMRC20, U343, T3M5, HCT15.
- Higher-depth six: ONE58, RPMI7951, KYSE140, HKGZCC, SKGII, KYSE450.

For each stratum, report primary-minus-source-mean changes in all five
source-stage readouts and primary absolute-error improvement over every
required comparator in the original protocol. Use its identical paired
target and crossed line/target bootstrap, with six lines in each stratum.
Report Spearman correlation of per-line primary-minus-source-only
absolute-error improvement with log total UMIs/KO as descriptive (12 lines).

**Additional necessary advancement gate:** absolute-error improvement over
each of the six original required comparators must be nonnegative at the
point estimate in **both** depth strata. This does not relax any original
gate or select a depth cutoff from outcomes. If the original test fails,
this addendum cannot rescue it. These depth groups do not separate assay
batch from biology and do not create additional experimental replication.

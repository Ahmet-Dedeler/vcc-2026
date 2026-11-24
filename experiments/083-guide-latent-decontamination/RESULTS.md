# 083 — Guide-latent decontamination results

Date: **2026-08-23**

Amended: **2026-08-30**

## Decision

**STOPPED AT STRUCTURAL PREFLIGHT. Do not launch heavy compute.**

The local Replogle K562 GWPS H5AD contains only **8 VCC-panel targets with at
least two distinct guide-pair observations**, and only **6** of those enter
the exact experiment-064 K562/PerturbMatch source slot, versus the
predeclared experiment-078 minimum of **100**. All eight are identifiable in
the pinned seed metadata, so missing sequence annotation is not the limiting
factor. The experimental unit itself lacks the required panel replication.

This is not a negative biological result. No guide pseudobulk cache, model fit,
held-donor effect evaluation, count generation, or submission was produced. It
is a design-feasibility rejection under the frozen hypothesis-3 promotion
contract.

On 2026-08-30, experiment 109 made one guarded cache-generation attempt using
this extractor. It streamed the complete source matrix, then failed the
unchanged `2e-5` compact-reconstruction gate before publishing a cache or
starting any fit. The cause was an implementation mismatch: guide effects
correctly excluded guides with fewer than 20 cells, but target reconstruction
incorrectly pooled only those retained guides, unlike the compact donor
builder, which pools all cells for targets with at least 20 total cells. The
extractor was repaired before any fit or held-outcome access by separating
those two estimands.

## Exact local evidence

| quantity | observed |
|---|---:|
| Cells | 1,989,578 |
| Measured genes | 8,248 |
| Guide-pair categories | 11,187 |
| Target labels | 9,867 |
| Perturbation targets with >=2 retained guide-pair constructs | 775 |
| VCC targets represented in pooled K562 compact | 272 / 300 |
| VCC targets with a retained >=20-cell guide-pair row | 269 / 300 |
| VCC targets with >=2 guide-pair constructs | **8 / 300** |
| VCC targets with >=2 sequence-identifiable constructs | **8 / 300** |
| Exact arm-D-evaluable sequence-identifiable multi-guide targets | **6 / 300** |
| Required for promotion | **100 / 300** |

The eight targets are `CAT`, `MED15`, `MYD88`, `MYRF`, `NABP2`, `PDCD4`,
`ZNF514`, and `ZNF768`; `ZNF514` and `ZNF768` lack PerturbMatch rows and
therefore keep the exact fallback. The experiment-068 seed-evidence targets (`PQBP1`,
`SHPRH`, `FZD2`, `GLYR1`) each have one K562 guide-pair construct and are
therefore preserved exactly by the specified fallback.

See `source-audit.json` for the machine-readable evidence and provenance.

## Implementation status

The reusable experiment package is complete:

- read-only H5AD/seed-metadata preflight;
- streaming extractor with >=20-cell guide effects and independent all-cell
  compact-target reconstruction;
- provenance-checked sequence/seed feature builder;
- leakage-safe fit/tune/test latent model;
- exact identity and unsupported-target fallbacks;
- pre-panel freeze boundary;
- independent experiment-064 arm-D reconstruction/evaluator;
- source-cluster bootstrap, permutation nulls, shrink control, and promotion
  gates;
- fail-closed AWS launcher/runner;
- sub-1-GiB local synthetic tests.

The implementation intentionally does not consume any experiment-080 cache or
evaluator artifact.

## What would reopen the experiment

Only a source that supplies at least 100 VCC-panel targets with two or more
independent guide-level expression observations can pass the first gate.
Paired protospacers within one observed construct do not count as two
expression replicates. If such a source appears, rerun `inspect_source.py`;
the guarded pipeline can then proceed without changing the preregistered
model, nulls, or promotion thresholds.

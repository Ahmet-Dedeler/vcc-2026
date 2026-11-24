# 074 — SIGrecon compact-signature audit

Date: 2026-08-23

## Decision

**Do not add SIGrecon to the production donor stack.**

The July 2026 Zenodo release is compact and valid, but both objects cover the
same 1,297 common-essential perturbations and have **0/300 overlap** with the
VCC 2026 panel:

| object | perturbations | VCC-panel overlap |
|---|---:|---:|
| `perturbseq.k562` | 1,297 | 0 |
| `perturbseq.rpe1` | 1,297 | 0 |

Each perturbation contains ranked gene-symbol arrays (`up` and `up_full`),
not numeric effect magnitudes. It therefore cannot provide direct donor rows,
and it is weaker than the existing numeric K562/RPE1 tables for training
generic cross-line structure.

The two checksum-verified source files and their provenance are recorded in
`data/external/sigrecon/meta.json`. No count generation, cloud job, scorer, or
challenge submission was run.

## Resource note

The Python `rdata` compatibility parser successfully decoded both complete R
objects. Peak local RSS was 5.55 GiB, with system free-memory percentage still
reported at 90%. No further full-object parsing is planned.

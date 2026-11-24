# Experiment 116 — DLD-1 low-MOI source artifact

Status: **source-only build complete and verified; formal held-context promotion
gate did not pass**

Date: **2026-08-30**

## Acquisition

Only the two pre-allowlisted official GEO LFC matrices were acquired. Each
object was downloaded twice into distinct partial files, and acceptance
required exact frozen HTTP metadata, exact byte count, complete gzip
validation, equal compressed and uncompressed SHA-256 digests, and a
byte-for-byte comparison. Both transfers began at byte zero in this run.

| Split | Bytes | Accepted SHA-256 | Uncompressed bytes | Uncompressed SHA-256 |
| --- | ---: | --- | ---: | --- |
| Low1 | 155,201,822 | `aaa70afe2dfefd22b72fe6294b87887bc6ee928590e419f459bcdf9616fba517` | 352,165,262 | `22218107cb60b71f5639fa6c371e69ee405bbc0d3a4e4f8087c096870eb9b313` |
| Low2 | 155,251,513 | `4e9bd64d8a3223dcb50d9c2a5b844698099a81dc2f86c38e4c0e270b4c87c390` | 352,617,876 | `6a2bf290d614d0740ba8d145618552479f63345fbfb7db61e13a39a4ffa03587` |

The 26.5 GB raw archive and every other processed object remain absent.

## Canonical artifact

The streaming build produced:

- 18,533 canonical response coordinates;
- 6,451 original author columns preserved separately in Low1 and Low2;
- 5,656 exactly mapped non-NTC promoter labels;
- 4,715 canonical target symbols;
- 759 targets with multiple promoter columns;
- 2,827 source response rows, of which 2,587 map exactly to the canonical
  response axis and 240 remain unsupported;
- 300/300 panel response coordinates set to exact zero and
  `measured=false` in every transfer array.

`author_lfc.h5` preserves finite canonical author values in float64 with an
exact measured mask, including NTC and unsupported-target columns. The
transfer arrays remain split-specific and use the frozen coordinatewise
finite median for duplicate promoters. A full second pass over both source
CSVs confirmed exact author values, exact masks, exact aggregation, untouched
fill rows, and exact panel exclusion for all 2,587 mapped response rows per
split.

## Source-only reliability

Across all 4,715 canonical target rows:

| Diagnostic | Result |
| --- | ---: |
| Common measured non-panel response coordinates | 2,536 |
| Median Low1/Low2 cosine | 0.069630 |
| Mean Low1/Low2 cosine | 0.092976 |
| Cosine p10 / p90 | 0.021754 / 0.175576 |
| Positive cosine rows | 4,608 / 4,715 |
| Median centered Pearson | 0.069773 |
| Median top-100 sign agreement | 0.610 |
| Mean top-100 sign agreement | 0.629 |
| Rows with cosine >=0.1 and top-100 sign agreement >=0.6 | 1,338 |

The frozen broad-coverage reliability mask has zero passing rows because the
author matrices contain only 2,536 common non-panel canonical response genes,
13.91% of the full 18,233-coordinate non-panel axis, below the preregistered
80% requirement. This is a coverage failure, not a missing-value discrepancy:
the common measured count is exactly 2,536 for every target.

Promoter variants are weakly concordant with each other:

| Diagnostic | Low1 | Low2 |
| --- | ---: | ---: |
| Finite promoter-pair comparisons | 1,202 | 1,202 |
| Median pairwise cosine | 0.008179 | 0.007284 |
| Mean pairwise cosine | 0.008681 | 0.008599 |

Two PECAM1 pairs have no common measured coordinates. This argues against
post-hoc selection of a favored promoter and supports retaining all promoter
columns and the deterministic frozen aggregation.

## VCC-overlap split retrieval

Exactly 67 VCC targets are represented. With all 300 VCC panel response
coordinates excluded:

| Diagnostic | Result |
| --- | ---: |
| Low1 query → Low2 gallery PDS | **0.943464** |
| Low2 query → Low1 gallery PDS | **0.950701** |
| Macro PDS | **0.947083** |
| Bootstrap 95% interval | **[0.915197, 0.971732]** |
| Permutation null p95 | **0.560380** |
| Permutation p-value | **0.00009999** |
| Perfect bidirectional macro rows | 38 / 67 |

The weakest target is `ZNF251` at macro PDS 0.2121; the next weakest is
`ALG12` at 0.6288. The remaining retrieval signal is exceptionally strong
despite modest absolute split cosine, matching the competition's
rank-discrimination objective better than a magnitude-similarity objective.

## Decision

The exact frozen promotion gate is **NO-GO**: seven of nine checks pass, but
median target cosine (0.0696 versus required 0.2) and broad-axis reliable-row
fraction (0 versus required 0.5) fail.

Therefore Experiment 116 does **not** authorize a held-context run. The source
is nevertheless scientifically promising for a future, separately frozen
method that is explicitly designed around sparse-response retrieval: its
67-target split PDS is far above chance and highly significant. Such a new
experiment would need its own source-only candidate definition and held
protocol; this result cannot be used to revise Experiment 116 after the fact.

## Key immutable hashes

```text
PREVALUE_FREEZE.v1.json
9d3d6dc2489f648c9a0586f2a4c27e99a6344ed1fd915e8d23f0ecab74c5fc7b

ACQUISITION_RECEIPT.v1.json
8115283b263381ed283d6ffc6e72f40a2de6882bc7917c623dae7bd26787f454

MATRIX_VALUE_ACCESS_LOCK.v1.json
14ce20c2f249d1fd6ad53344e39240b79929751c117ad24ae25eecd115e524c5

artifact-manifest.json
c592e2ec51cb56e7df085224b3ff75d8735ce87280683fd0967135472e22be79

result.json
6e90fd909059f19e8140235ba62bd458b482802813c7679c5bdca4b33d3bfc28

author_lfc.h5
03d1bdec0556a8a143093fdda270defb3bc6756d9c4c07db554d41cbb38210e5

transfer_lfc_low1.npy
9eb958b76556d6f25aa94c6989b965b0ee895eba69a6cdb73b3e25dfc1ea7f6c

transfer_lfc_low2.npy
c72b095bf12b3a20cb7d4a6d8bffaf1e3e47b2ca5d1e656779f1699876368533

transfer_measured_low1.npy
8c7019eac573af5b5ac05865c427b8f2e7f4be3447821ebb09c20c97f0f7da64

transfer_measured_low2.npy
8c7019eac573af5b5ac05865c427b8f2e7f4be3447821ebb09c20c97f0f7da64
```

The complete per-file ledger is
`work/compact-source-v1/artifact-manifest.json`.

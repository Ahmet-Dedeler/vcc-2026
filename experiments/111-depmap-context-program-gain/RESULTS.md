# Experiment 111 — source-only result

Status: **NO-GO**

Date: **2026-08-30**

## Bottom line

The preregistered DepMap context-program operator produced a small positive
nested source-test gain, but it did not clear the frozen promotion gates:

| Measure | Result | Gate | Pass |
| --- | ---: | ---: | :---: |
| Macro raw PDS delta | `+0.00490196` | `>= +0.010` | no |
| Nonnegative outer batches | `5 / 6` | `>= 5 / 6` | yes |
| Worst outer-batch delta | approximately `0` | `>= -0.003` | yes |
| Target-bootstrap lower 95% | `-0.00095316` | `>= 0` | no |
| Context-label null p95 | `+0.00680828` | observed above p95 | no |
| Target-pairing null p95 | `+0.00322259` | observed above p95 | yes |
| Shallow-effect half | `+0.00944081` | `>= 0` | yes |
| Deep-effect half | `+0.00036311` | `>= 0` | yes |

The effect beats the target-pairing null but not the context-label null. The
largest outer gain is confined to `apollo31` (`+0.01960784`); the other five
batches range from numerical zero to `+0.00326797`. This does not establish a
robust control-conditioned context effect.

The final tune-only selector preferred
`rr4_cr4_l1_b1_cap5` (response rank 4, four control PCs, ridge 1, beta 1,
five-degree cap), with tune delta `+0.01034858`. Because the untouched nested
test and null gates failed, this arm was **not** frozen as a candidate.

## Frozen source inventory

- Six authenticated allowlisted CSV members only.
- Selective compressed payload: `300,011,660` bytes.
- Published source directory: `data/external/depmap26q1/extracted`.
- Pseudobulk matrix: 34,648 response rows and 1,593 context-target columns.
- SCEPTRE LFC/FDR matrices: 22,966 response rows and 1,440 columns.
- Modeling axis: 15,211 globally complete source genes on the challenge gene
  axis after removing all 300 VCC targets and all 100 Broad targets.
- 90 non-control source targets; 87 are complete across all 16 contexts.
- Six experimental batches, 16 cell lines, 11 lineages.

## Outer results

| Batch | Raw PDS delta | Selected arm |
| --- | ---: | --- |
| `apollo31` | `+0.01960784` | `rr4_cr4_l1_b1_cap5` |
| `apollo32` | approximately `0` | `rr8_cr4_l1_b1_cap5` |
| `apollo33` | `+0.00326797` | `rr4_cr4_l1_b1_cap3` |
| `apollo34` | `+0.00326797` | `rr4_cr2_l1_b1_cap5` |
| `apollo35` | `+0.00326797` | `rr8_cr4_l1_b1_cap5` |
| `apollo38` | approximately `0` | `rr4_cr4_l1_b1_cap5` |

The deterministic wrong-context diagnostic gives only `+0.00068083`, while
the real mapping gives `+0.00490196`. That difference is suggestive, but the
100 complete context-label permutations reach p95 `+0.00680828`
(`p = 0.1683`), so it is not reliable enough to spend held truth.

## Acquisition recovery

The first selective attempt failed closed with:

```text
ERROR: crispr_table.csv row 1442 gene_effect is empty
```

The exact source member contains 70 such rows, all olfactory controls, with
only `gene_effect` and `dependency_strength_class` blank and
`is_dependent=False`. The recovery amendment pins the exact member SHA-256 and
the exact 70-pair digest; no analysis-target missingness was allowed.

The second attempt atomically published the six validated files, then failed
while writing the receipt because the old helper admitted only
experiment-local roots:

```text
ERROR: path escapes ~/Code/virtual-cell-challenge-2026/experiments/111-depmap-context-program-gain: ~/Code/virtual-cell-challenge-2026/data/external/depmap26q1/extracted
```

The receipt was recovered from the already-published exact files after full
local hash/schema validation and a fresh remote central-directory audit. The
selected payload was not downloaded a third time.

## Resource use

- First selective attempt: 74.16 s wall, 133,070,848 bytes maximum RSS.
- Second selective attempt: 120.83 s wall, 133,021,696 bytes maximum RSS.
- Receipt recovery and full local schema audit: 35.44 s wall, 118,292,480
  bytes maximum RSS.
- Source-only nested fit plus 100 context-label and 100 target-pairing complete
  selections: 1,009.15 s wall, 1,948,368,896 bytes maximum RSS; 21 unique
  response-basis fits.

## Immutable outputs

`work/source-only-v1/`:

| File | SHA-256 |
| --- | --- |
| `artifact-manifest.json` | `a02bd4967ac5eda9ec7cd171d3c332482fdfe973b7375ae0701c303da2056648` |
| `result.json` | `3b5855b4a69042af4958711f4ed4475cea9009550fb3a2d2259a4aa7ff44aaaa` |
| `context-null-values.npy` | `09c6e79feaea858a3897345156ae82f68aa6f4e4c2a6c07a2bbac1ea816f9382` |
| `target-null-values.npy` | `80fbf5223f6458388bf5bc998f972ff11b39216ace564500caae570ddfcf2876` |
| `target-bootstrap-values.npy` | `94a9f35e1b6da7ceb4eccb94d633aab769edce87fda02f4f57b0f412e7c81521` |
| `genes.txt` | `c699b0ee0d23cf091f12e64511c41b6d651e3c8acec2e1f5d555705e4d82b05e` |
| `targets.txt` | `6e952a48ca1e8f7f7c5326f41b789d394fe23f311374e83e2b3edf023b4ce1b7` |

The artifact-manifest verification passes for all six payload files.
`candidate.npz` and `candidate.json` are absent by design.

## Decision

Close this exact context-program family as a source-only **NO-GO**. Do not open
H1 or Orion truth for it, do not fit a post-outcome variant, and do not create
a held-evaluation proposal.

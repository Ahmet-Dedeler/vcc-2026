# Experiment 119 results

Status: **source acquisition authenticated and compact source verified**

No held perturbation truth, existing outcome directory, evaluation result,
count-generation path, or submission path was accessed.

## Authenticated source

- file: `data/external/altair-a549/altair-scrna.h5ad`
- bytes: `281474102`
- SHA-256:
  `4632972e6c2871603f4ca9108489787d7761aae05dc6a4bf3488f1ed78a63c9f`
- official rows × genes: `1001 × 9999`
- targeting rows: `1000`
- pooled non-targeting row: `1`
- represented cells: `606075`

## Verified semantics

`.X` is the official OPS `log2_fold_change` matrix. The frozen `ALG12` API
response agrees with the H5AD row to the API's five-decimal serialization
precision, and `.X / lfcSE` reproduces the source `stat` layer.

The compact artifact preserves:

- `.X` as `log2_fold_change.npy`;
- `layers/neg_log10_fdr` as `neg_log10_fdr.npy`;
- finite-LFC support as `measured_mask.npy`;
- finite-FDR support as `fdr_mask.npy`;
- every unsupported coordinate as NaN.

No normalization, rescaling, imputation, alias resolution, or symbol merging
is performed.

## Coverage

- source response genes: `9999`
- exact canonical symbol matches: `7843`
- source symbols not on the VCC axis: `2156`
- canonical genes absent from the source: `10690`
- source-supported targeting rows: `997/1000`
- exact VCC target-label overlap: `21/300`
- VCC targets with measured source effects: `20/300`

The all-NaN targeting rows are `COQ10B`, `HSD17B12`, and `NEK1`.
`Nontargeting` is also all-NaN and is intentionally excluded from the compact
target axis.

## Compact artifact

Location:

```text
experiments/119-altair-a549/work/compact-source/
```

- shape: `1000 × 18533`
- total artifact bytes: `185979662`
- finite canonical LFC values: `7819471`
- finite canonical FDR values: `6088888`
- artifact-manifest SHA-256:
  `4b5415b8f73db391d6a073b63066a719d1abb69c0ef9383bbf0ffc5687e818b0`
- compact receipt SHA-256:
  `6f26ed7994a07ba8bb9ff6c1763fe475843af8a8801464cd10d714abf5aaae62`

Spot checks against the authenticated source were bit-exact for finite
float32 values, including `ALG12/CPLX2` and `ZNF496/EEPD1`. The `COQ10B` row
remains entirely NaN with a false measured mask.

## Verification

```text
12 unit/integration tests passed
ruff: all checks passed
source audit: verified
compact artifact: verified
git diff --check: clean
```

This artifact is a source representation only. Its log2 fold changes are not
a VCC donor-space `delta.npy` and have not been tested for predictive value.

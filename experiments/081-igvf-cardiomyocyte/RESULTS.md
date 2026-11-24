# 081 — IGVF cardiomyocyte sparse-direction transfer

Date: **2026-08-23**

## Decision

**NO-GO. Do not add this IGVF cardiomyocyte adapter to the production
predictor.**

The frozen non-panel selection chose signed-square-root log2FC, global L2
scaling, and fraction `0.20`. That candidate looked positive on the untouched
non-panel test, but its target-cluster interval still crossed zero. More
decisively, the one-time complete VCC-panel evaluation was negative in all
three held-donor folds.

| stage | K562 delta | HCT116 delta | HEK293T delta | macro delta |
|---|---:|---:|---:|---:|
| non-panel tune | +0.003917 | +0.008588 | +0.006105 | **+0.006203** |
| untouched non-panel test | +0.004193 | +0.000165 | +0.005392 | **+0.003250** |
| one-time VCC-panel final | −0.000584 | −0.000898 | −0.000718 | **−0.000733** |

The panel rejection is not a gallery-coverage artifact:

- all 259 targets shared by K562, HCT116, and HEK293T remained in the final
  retrieval gallery;
- 81 targets had at least one usable IGVF response edge after removing every
  VCC target-gene coordinate;
- all 178 unsupported targets were exact incumbent fallbacks, with measured
  PDS delta exactly zero in every fold;
- among the 81 modified targets, mean macro delta was `−0.002345`.

## Frozen selection

The experiment used 1,514 non-panel targets with at least one official-axis
IGVF response and complete K562/HCT116/HEK293T coverage:

| split | targets | use |
|---|---:|---|
| fit | 869 | log2FC clipping and source-side scale constants only |
| tune | 314 | representation, scale mode, and fraction selection |
| test | 331 | untouched non-panel pre-gate |

The tune rule maximized worst-fold PDS delta, then macro delta, then preferred
the smaller fraction. Identity was mandatory. It selected:

```text
representation = signed_sqrt
scale_mode      = global_l2
fraction        = 0.20
```

All three tune folds were positive. Nearby raw-log2FC and sign candidates were
also positive, so this was not a one-cell grid accident. The transfer did not
generalize to the VCC panel.

## Uncertainty and nulls

### Untouched non-panel test

- Macro delta: `+0.003250`
- Target-cluster bootstrap 95% CI:
  `[-0.001052, +0.009113]`
- Positive bootstrap probability: `0.9288`
- Exact-sparsity sign-null p95: `+0.001112`
- Matched-label null p95: `+0.002571`

The point estimate exceeded both null p95 values and the `+0.003` effect gate,
but the cluster-bootstrap lower bound was negative. The non-panel pre-gate
therefore failed before considering the panel.

### One-time VCC-panel final

- Macro delta: `−0.000733`
- Target-cluster bootstrap 95% CI:
  `[-0.004303, +0.002609]`
- Positive bootstrap probability: `0.33865`
- Exact-sparsity sign-null p95: `+0.000085`
- Matched-label null p95: `+0.000231`

The observed panel result was below both null thresholds and negative in every
held fold. It failed every scientific promotion criterion other than format,
gallery completeness, and fallback preservation.

## Interpretation

The compact DE artifact contains a real transferable signal on the broad
non-panel TF catalogue, but it is not stable on this VCC target panel.
The likely reasons are:

1. **lineage specificity:** cardiomyocyte TF responses contain strong
   programmes that transfer among broad TF perturbations but not reliably to
   the VCC panel;
2. **extreme sparsity:** the filtered table has a median of only a few reported
   genes per perturbation, so a small number of high-impact targets dominate
   the non-panel mean;
3. **panel composition shift:** positive non-panel outliers such as `ISL1` and
   `GATA4` have no panel analogue strong enough to offset negative panel
   targets such as `PIN1`, `PLAGL2`, and `MYRF`.

Do not rescue this result by choosing a different fraction or target subset
from the panel outcomes. That would violate the frozen one-time panel gate.
A future IGVF revisit needs materially richer released directions—an
unfiltered full-gene effect table or another compact lineage—not panel-driven
retuning of this sparse artifact.

## Data and provenance

Only the two authorized small files were downloaded:

| file | bytes | official MD5 | observed status |
|---|---:|---|---|
| `IGVFFI0830FXFI.csv.gz` | 1,887,277 | `39097b9eb95bb1e43e245480d242fea9` | match |
| `IGVFFI8270UPKB.csv.gz` | 348,103 | `766b1fd06fe7d641104d2569a735a63f` | match |

Schema audit:

- filtered DE: 27,578 rows, 1,889 target labels, 8,236 normalized response
  symbols, and 120 duplicate target-response pairs aggregated by
  cell-count-weighted mean;
- guide manifest: 14,151 tab-delimited rows and 2,062 gene symbols;
- raw shared coverage: 1,525 non-panel targets and 85 panel targets;
- usable official-axis coverage: 1,514 non-panel targets, 84 panel targets
  before coordinate exclusion, and 81 panel targets after exclusion.

The exact portal/API provenance and observed SHA-256 hashes are recorded in
`PROVENANCE.json`.

## Leakage and geometry audit

- Every held fold used the other two main donors as predictors.
- Each frozen gamma-0.65 union operator excludes the full 300-gene VCC panel
  and was fitted without its held donor.
- All 1,514 non-panel proxy target coordinates were additionally zeroed before
  transformation and excluded from proxy retrieval, preventing on-target
  knockdown from selecting the adapter.
- Fit targets set robust clipping and scale constants only.
- Tune targets selected the representation/scale/fraction.
- Test targets were opened after selection.
- Target-specific held VCC-panel rows were extracted only after the non-panel
  test completed.
- The IGVF residual was added after incumbent whitening and only on released
  target-response edges, preserving every missing target and gene.
- Bootstrap clusters came from incumbent source-prediction similarity, not
  held-truth outcomes.
- No count emission, `cell-eval2`, cloud work, or challenge submission ran.

## Resource and validation audit

- Peak measured process-tree RSS: **0.69 GiB**
- Internal ceiling: `0.90 GiB`
- External process-tree ceiling: `0.98 GiB`
- System free memory remained approximately 90%.
- Eight synthetic tests passed.
- Ruff passed with no findings.
- Aggregate scores were independently recomputed from `per_target_pds.csv`.
- Missing-target deltas were independently verified as exactly zero.
- Result-embedded hashes for `run.py`, `test_run.py`, and `PROVENANCE.json`
  match the current files.

Authoritative artifacts:

- `results.json`
- `tune_grid.csv`
- `per_target_pds.csv`
- `PROVENANCE.json`
- `run.py`
- `test_run.py`

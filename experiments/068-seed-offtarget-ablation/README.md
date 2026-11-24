# 068 — Seed-off-target source ablation

## Question

Does removing only the public-source effect vector implicated by a reported
seed-driven off-target improve effect-level PDS over frozen experiment-064 arm
`D`?

The four predeclared VCC-panel targets are:

| target | reported off-target | implicated source |
|---|---|---|
| `PQBP1` | `TIMM17B` | Replogle K562 and Xaira/Orion HCT116 |
| `SHPRH` | `CENPP` | Xaira/Orion HCT116 |
| `FZD2` | `P2RX7` | Xaira/Orion HCT116 |
| `GLYR1` | `UBN1` | Zhu/Marson CD4 |

The source claims are verified from the compressed candidate tables in the
public `AustinHartman/perturb_seed_reproduce` repository at commit
`198fe3d84263002e9cfa17e1a2178d4e8dd0a669`, not copied from prose.

## Frozen intervention

- Comparator: experiment-064 arm `D`.
- K562 claim: remove the entire 50:50 incumbent/PerturbMatch K562 slot because
  both versions derive from the same Replogle assay.
- HCT116 claim: remove only the Orion HCT116 source row.
- CD4 claim: remove the complete post-whitening CD4 mean3 residual. The public
  table measures the Zhu CD4 guide in the 8-hour condition; treating the
  sequence-level seed match as applying to the CD4 source family is the one
  conservative extrapolation tested here.
- Do not zero the reported off-target gene or inject any literature effect.
- Keep the frozen gamma-0.65 operators, source scaling, 50:50 PerturbMatch
  blend, CD4 fraction 0.05, panel-coordinate exclusion, and experiment-064
  arm-B preserve scalar.

If all ordinary somatic source rows are removed for a target, retain the
existing experiment-064 missing-target fallback: the normalized mean source
response. This changes only source-row availability, not the operator or any
other target.

## Held-truth cases

Primary cases use the exact two experiment-064 donor-held-out geometries:

- held HCT116: K562/PerturbMatch + HEK293T sources;
- held HEK293T: K562/PerturbMatch + HCT116 sources.

Only cases where the implicated source is actually present are estimands.
`GLYR1` has no HCT116 or HEK293T truth row, so it cannot be evaluated in the
native 064 folds.

Secondary cases score the same two frozen 064 prediction geometries against
independent held KOLF truth. KOLF is never a prediction source or operator
source. This extension is reported separately and is used only to cover
`GLYR1` and replicate the other rules; it cannot by itself authorize shipping.

Each query ranks against every VCC target available in that held-truth bank:
290 HCT116, 294 HEK293T, or 273 KOLF targets. All 300 panel-gene coordinates
remain excluded from transformations and retrieval.

## Matched null

For each reported target, controls are sampled deterministically from other
VCC targets with:

- the same held-truth cases;
- the same source-ablation mask in each case;
- the exact same K562, PerturbMatch, companion-donor, CD4-condition, and
  held-truth availability pattern.

Up to 64 controls are sampled with seed `6802026`; smaller exact-matching
classes are exhausted. The target-level null applies the reported target's
source-removal rule to each matched control and averages its paired PDS delta
over the same held cases. The overall null draws one matched control per
reported target for 20,000 deterministic replicates.

No held-truth effect, arm score, target identity embedding, or challenge result
is used to select controls or define a rule.

## Decision rule

Ship only if all of the following hold:

1. the native experiment-064 cases have positive mean PDS delta;
2. at least three of four target-cluster deltas are positive and none is
   negative;
3. the four-target mean exceeds the matched-null 95th percentile;
4. the low-memory reconstruction reproduces available experiment-064 arm-D
   ranks exactly.

Otherwise reject the hard-coded source-ablation rule.

## Safety

- One Python process, one BLAS thread, no subprocesses.
- Internal stop threshold: 0.70 GiB RSS, below the hard 0.75 GiB cap.
- Stored `.npz` eigenvectors are streamed directly in row blocks; they are
  never loaded or memory-mapped.
- Public source tables are streamed and decompressed in memory-bounded mode.
  Each compressed object is checked against a 50,000,000-byte download cap.
- No `cell-eval2`, cloud launch, count emission, submission, or write outside
  this experiment directory.

# 058 — PerturbMatch K562-GWPS reanalysis results

Date: 2026-08-22

## Verdict

**Do not replace the incumbent K562 donor with PerturbMatch. Keep a 50:50
K562-processing ensemble as a promising follow-up, but do not change v9 from
this experiment alone.**

PerturbMatch is slightly worse by itself and as a full replacement after the
gamma-0.65 source whitening used by the current model. Its errors are
complementary, however: averaging the incumbent and PerturbMatch K562 vectors
50:50 improves both independent Orion holdouts, both before and after
whitening.

The production-shaped primary comparison is:

| raw effects, source-whitened gamma 0.65 | macro PDS | delta | HCT116 | HEK293T |
|---|---:|---:|---:|---:|
| incumbent K562 + Orion baseline | 0.794525 | — | 0.804099 | 0.784952 |
| PerturbMatch replacement + Orion | 0.793958 | -0.000567 | 0.801635 | 0.786282 |
| 10% PerturbMatch in K562 slot | 0.796227 | +0.001702 | 0.805819 | 0.786636 |
| 25% PerturbMatch in K562 slot | 0.797920 | +0.003395 | 0.808141 | 0.787700 |
| **50% PerturbMatch in K562 slot** | **0.799693** | **+0.005168** | **0.808974** | **0.790412** |
| 75% PerturbMatch in K562 slot | 0.798514 | +0.003989 | 0.806829 | 0.790200 |
| 90% PerturbMatch in K562 slot | 0.796538 | +0.002012 | 0.804453 | 0.788622 |

Both held lines independently select the 50% blend from the predeclared grid.
The paired target-cluster bootstrap interval for its macro delta is
`[-0.000763, +0.011311]`, so the two-fold gain is encouraging but not enough
to call production-safe.

## Direct K562 direction

The same complementarity exists before adding the second source donor:

| raw effects, source-whitened gamma 0.65 | macro PDS | delta vs incumbent K562 |
|---|---:|---:|
| incumbent K562 direct | 0.762100 | — |
| PerturbMatch direct | 0.758084 | -0.004016 |
| 10% PerturbMatch blend | 0.765991 | +0.003891 |
| 25% PerturbMatch blend | 0.769413 | +0.007313 |
| **50% PerturbMatch blend** | **0.771540** | **+0.009440** |
| 75% PerturbMatch blend | 0.767321 | +0.005221 |
| 90% PerturbMatch blend | 0.761710 | -0.000390 |

The two K562 processing estimates have mean target-wise cosine `0.620396`
(median `0.618822`, p10 `0.570612`, p90 `0.682372`) on the evaluated axis.
PerturbMatch is therefore not merely a rescaled copy. It is weaker alone but
adds useful orthogonal direction when averaged with the incumbent estimate.

## Unwhitened control

The result is not created by the whitener:

| raw effects, no whitening | macro PDS | delta |
|---|---:|---:|
| incumbent K562 + Orion baseline | 0.702939 | — |
| PerturbMatch replacement + Orion | 0.705226 | +0.002287 |
| 10% PerturbMatch in K562 slot | 0.704260 | +0.001321 |
| 25% PerturbMatch in K562 slot | 0.706441 | +0.003501 |
| **50% PerturbMatch in K562 slot** | **0.708169** | **+0.005230** |
| 75% PerturbMatch in K562 slot | 0.707620 | +0.004680 |
| 90% PerturbMatch in K562 slot | 0.706618 | +0.003679 |

The direct unwhitened comparison is incumbent `0.720021` versus PerturbMatch
`0.717991` (`-0.002030`). A 50:50 direct K562-processing ensemble reaches
`0.734390` (`+0.014369`) and improves both held folds.

Gene normalization is weaker: the unwhitened full 50:50 ensemble gains only
`+0.002181`, while full replacement loses `-0.002154`. Raw effects remain the
primary representation.

## Protocol

- Valid holdouts: Orion HCT116 and Orion HEK293T.
- K562 is deliberately omitted as a holdout because both K562 inputs derive
  from the same Replogle experiment.
- Retrieval gallery: all **238** VCC targets shared by PerturbMatch, incumbent
  K562, HCT116, and HEK293T.
- Feature axis: **6,949** official genes shared by all four banks.
- All **300** challenge target-gene coordinates are excluded from scoring and
  whitening, not merely the 238 query targets.
- Every query is ranked against the complete 238-target gallery.
- Source centers and covariance bases exclude all 300 panel perturbation rows.
- Held-line truth remains raw and uncentered.
- The whitening operator is candidate-independent and uses only incumbent
  K562 plus the non-held Orion source. PerturbMatch cannot tune the geometry it
  is being judged under.
- Basis: up to 6,000 non-panel perturbations per source, deterministic seed.
- No count emission, `cell-eval2`, challenge prediction, or submission.

## Acquisition and exact format

The authoritative Genentech Zenodo artifact is
`lfc_matrix_ALL.csv.gz`, 748,423,850 bytes:

- source MD5: `0017076a5cc7b9979febf3bd19bd2157`;
- source SHA-256:
  `bca0ce46f573d316500ee6006e9f063d2f42a9745040b3e37554d3eec44b8305`;
- license: CC BY 4.0;
- source shape: 9,761 gene rows x 9,471 `<target>_K5` columns;
- VCC coverage: 251/300 targets.

The streaming extractor writes a 9,471 x 9,761 Fortran-order float32
`delta.npy` without retaining or materializing the decompressed CSV.

- local `delta.npy` SHA-256:
  `efa20944f6f6485a7cc7f484c7ab70bd2ba491a6e07b231f306e104b2e1e4cf2`;
- compact archive SHA-256:
  `2fe98cb106851f48bbfcccd98ec909488c88eabf80caf84f1508419046d9117b`.

The separate processed SCE is described as top-3,000 HVGs; the public `ALL`
LFC CSV itself has 9,761 gene rows. The extractor initially treated the prose
as the matrix shape, failed closed at row 3,001, then was corrected to derive
and strictly verify the complete 9,761-row stream.

## Memory and verification

- Extraction peak RSS: **0.809 GiB**.
- LODO plus two 6,949-dimensional eigendecompositions peak RSS:
  **2.920 GiB**.
- Local synthetic tests: **4 passed**, 157 MiB maximum RSS.
- Cloud evaluator exit status: `0`.
- Result SHA-256:
  `1037fce3033230f4baacf195015dceddee8e5e48fc43eb6edbe8f10615fe9f1b`.
- Full outputs: `results.json`, `results.csv`, `per_target_pds.csv`,
  `cloud-extract.log`, and `cloud-evaluate.log`.

## Decision

PerturbMatch should not replace K562, and it should not be added to v9 yet.
The single next worthwhile test is a **50:50 incumbent/PerturbMatch K562 slot
on the exact production union axis**, with the operator and all other model
choices frozen. That test needs another dedicated heavy run and is
intentionally outside experiment 058.

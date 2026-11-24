# 065 — PerturbMatch + post-whitening CD4 HCT116 count result

Date: **2026-08-23**

## Decision

**Positive enough to build an exploratory production-shaped successor to v9,
but not direct submission authority.**

The frozen joint candidate improved the clean 132-target Orion HCT116
count-level comparator:

| metric | baseline raw | candidate raw | scaled change |
| --- | ---: | ---: | ---: |
| PDS cosine | `0.7978834143` | `0.8025098311` | **`+0.0102241256`** |
| MSE normalized | `2.0537529303` | `2.0631046951` | `0` |
| direction fidelity | `0.5281190914` | `0.5269334706` | `-0.0039520693` |
| direction reach | `0.0150248423` | `0.0145488827` | `-0.0005312049` |
| DE Jaccard | `0.0337623585` | `0.0345378913` | `+0.0020960347` |
| LFC NMAE | `0.9570910100` | `0.9545580735` | `+0.0042124339` |
| **six-metric average** | **`0.1216402396`** | **`0.1236484596`** | **`+0.0020082200`** |

The PDS gain is the main positive signal. NMAE and Jaccard also improve, while
fidelity and reach regress slightly. MSE remains above the official baseline
and earns no scaled credit.

Relative to PerturbMatch alone from experiment 061, the joint candidate raises
raw PDS from `0.7981725653` to `0.8025098311` and the six-metric proxy from
`0.1217308225` to `0.1236484596`.

## Interpretation

This is the first production-shaped count-level evidence that the frozen
experiment-064 PerturbMatch plus post-whitening CD4 direction survives beyond
effect-level PDS. The gain is modest and measured on one held cell line, so it
authorizes only construction and inspection of the exact challenge-context
package.

Experiment 065 is not mathematically identical to v9 production:

- it uses alpha `1.08` rather than `1.0796`;
- it uses a baseline-derived operator/global rescaling rather than v9's exact
  candidate-specific preservation semantics; and
- its control selection and count emission share one RNG stream, whereas v9
  uses separate deterministic streams.

Experiment 071 therefore treats this result as an exploratory packaging gate,
not as permission to submit.

## Execution and preservation

- Run ID: `20260823T040540Z`
- Worker: `i-0343d04439564b709`, on-demand `r7i.2xlarge`
- `cell-eval2`: `0.16.0`, preset `vcc2026`
- Build peak process-tree RSS: `8.8168 GiB`
- Scorer peak process-tree RSS: `15.6500 GiB`
- Exit code: `0`
- Candidate H5AD SHA-256:
  `d960cbc2ee2482c3894e8ded5ddc7aa96af65b039526ff1ef7a2f68294b3fde2`
- Durable cloud prefix:
  `s3://<your-s3-bucket>/output/exp065-pm-cd4-hct-count/20260823T040540Z/`
- Compact local mirror:
  `cloud-output/20260823T040540Z/`

All downloaded compact artifacts present in `artifact-sha256.txt` were
checksum-verified. The 2.28 GB candidate H5AD was deliberately not downloaded;
its independent checksum file and manifest entry are preserved locally.

No challenge submission was made.

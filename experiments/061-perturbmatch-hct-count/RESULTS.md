# 061 — PerturbMatch HCT116 count-level result

Date: 2026-08-23

## Decision

**Pass the non-negative prerequisite, but do not ship PerturbMatch alone.**

The frozen 50:50 incumbent/PerturbMatch K562-slot ensemble produced a very
small positive six-metric change against the exact clean gamma-0.65 HCT116
baseline:

| metric | baseline raw | candidate raw | favorable raw change | scaled change |
|---|---:|---:|---:|---:|
| PDS cosine | 0.7978834143 | 0.7981725653 | +0.0002891511 | +0.0006390079 |
| MSE normalized | 2.0537529303 | 2.0678142616 | -0.0140613314 | 0 |
| direction fidelity | 0.5281190914 | 0.5265800842 | -0.0015390072 | -0.0051300238 |
| direction reach | 0.0150248423 | 0.0143803864 | -0.0006444559 | -0.0007192588 |
| DE Jaccard | 0.0337623585 | 0.0345113612 | +0.0007490027 | +0.0020243315 |
| LFC NMAE | 0.9570910100 | 0.9548484973 | +0.0022425127 | +0.0037294408 |
| **six-metric average** | **0.1216402396** | **0.1217308225** | — | **+0.0000905829** |

The count-level result is directionally consistent with experiment 058, but
far smaller than its effect-level PDS gain. Improvements in PDS, Jaccard, and
NMAE are nearly cancelled by losses in fidelity and reach. This is below a
practical shipping threshold and should not replace v9 by itself.

The predeclared prerequisite for experiment 065 was merely a non-negative
six-metric average. That gate passes, so the frozen PerturbMatch plus
post-whitening CD4 composition advances to its one count-level confirmation.

## Execution and evidence

- Worker: `i-0343d04439564b709`, on-demand `r7i.2xlarge`.
- Candidate: 2,282,477,819 bytes.
- Candidate SHA-256:
  `61990b2bbef737ec7f1fdd028735adcac0c29c922aa2357a193aac8bae4daf71`.
- Candidate aggregate SHA-256:
  `f469c2c445cc86eb77cfccce141a4639f23fee2bbcb0eb9000162bb9eed544de`.
- Baseline aggregate SHA-256:
  `f1011cd12786d2196949daa6ce3bf328d52c0afc9d46966e8d3e882f6200cccf`.
- Result SHA-256:
  `636ebb12382effdf782e3e12310517c54d6f3170d087cc6c49f521a4e9912334`.
- Exit code: `0`.
- Durable S3 prefix:
  `s3://<your-s3-bucket>/output/exp061-perturbmatch-hct-count/attempt2/`.
- Compact result, aggregate, logs, manifest, and exit code are preserved under
  `cloud-output/`.
- No challenge submission was made.


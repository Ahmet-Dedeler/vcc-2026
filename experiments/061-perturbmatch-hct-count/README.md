# 061 — PerturbMatch K562-slot ensemble at HCT116 count level

Status: candidate built successfully; isolated AWS worker is waiting for the
clean gamma-0.65 baseline before serial `cell-eval2` scoring.

This experiment is the single production-shaped follow-up authorized from
experiment 058. It compares the existing clean HCT116 gamma-0.65 baseline with
one candidate:

- keep Orion HEK293T unchanged;
- replace the incumbent K562 slot with a 50:50 ensemble of independently
  normalized raw effects from K562-GWPS and PerturbMatch;
- if PerturbMatch lacks a target, retain the incumbent K562 slot exactly;
- if PerturbMatch lacks a gene, retain the incumbent K562 coordinate exactly;
- load the already fitted, candidate-independent union whitener from experiment
  040 rather than fitting another operator;
- require that all 300 VCC target coordinates are absent from that operator;
- freeze gamma 0.65, alpha 1.08, CPM gate 5, quasi exact-total emission,
  spread 0.25, and seed 0;
- score with `cell-eval2==0.16.0`, preset `vcc2026`, serially on AWS.

The baseline is not rebuilt. The comparison waits for the existing clean
gamma-0.65 HCT116 `agg_results.csv` and fails closed if it never appears.
No challenge submission is created or sent.

## Current run

- worker: `i-0343d04439564b709`, on-demand `r7i.2xlarge`, `us-east-1a`;
- output prefix:
  `s3://<your-s3-bucket>/output/exp061-perturbmatch-hct-count/attempt2/`;
- candidate size: `2,282,477,819` bytes;
- candidate SHA-256:
  `61990b2bbef737ec7f1fdd028735adcac0c29c922aa2357a193aac8bae4daf71`;
- build peak RSS: `9.948 GiB`;
- target-coordinate exclusion: all 300 verified absent;
- modified HCT targets: 123/132; nine PerturbMatch-missing targets were
  preserved bit-for-bit from the incumbent.

The first cloud attempt is quarantined: two runners reached the same HDF5
output path and one failed on the file lock. The retry adds a process-wide
`flock`, removes partial local output before building, and writes to the
separate `attempt2` S3 prefix. Only the `attempt2` artifact above is valid.

# Experiment 136 — DepMap residual margin screen

Status: **terminal source-development NO-GO**

The rejected Experiment 132 Ridge response was re-expressed as a residual
around the exact same-target source mean and evaluated only on Experiment
132's already-consumed 613-target development split.

## Result

- Every nonzero global dose from `0.0025` through `1.0` reduced macro PDS.
- At dose `0.0025`, macro PDS already fell by `-0.0002590`; every scored
  context was negative.
- Directional signed credit increased at small doses, reaching `+0.010110`
  at dose `0.20`, but macro PDS simultaneously fell by `-0.021078`.
- A deployment-available rule that chose the largest dose preserving observed
  source-context PDS rank still lost `-0.020302` macro PDS.
- A stricter nearest-decoy-margin rule selected identity for all but one
  target/context and produced effectively zero gain.
- Even an invalid oracle choosing one dose per target from the held development
  outcomes achieved only `+0.001926` macro PDS, below the project-wide
  `+0.003` gate, and its K562-family delta remained negative
  (`-0.000528`).

The residual contains broad directional signal, but the retrieval cost begins
immediately and its oracle ceiling is too small. Do not open Experiment 132's
untouched 410-target source-test split for this proposal family. A future
margin-constrained experiment needs a different residual proposal, not a more
elaborate router around this Ridge output.

No untouched test outcome, held challenge outcome, `cell-eval2`, counts,
challenge artifact, submission, commit, or push was used.

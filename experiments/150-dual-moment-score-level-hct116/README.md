# Experiment 150 — score-level dual-moment arbitration on held HCT116

Status: **completed; `TERMINAL_NO_GO`**

Date: **2026-08-31**

## New hypothesis

Experiments 148 and 149 used the same fixed dual-moment transformation and
both improved projected overall score with exact PDS preservation:

- H1-test: `+0.007213` over all 85 targets;
- held HEK293T: `+0.012919`, bootstrap 95% interval
  `[+0.009955, +0.015857]`.

The metric carrying the gain differed. H1 improved mainly through NMAE/reach,
while HEK293T improved mainly through fidelity/Jaccard. Experiment 149 is still
a formal no-go under its original mechanism-specific gates.

This new experiment asks a different, preregistered question: does the same
fixed transformation improve the total six-metric score in a third held
context, regardless of which DE member supplies the gain?

## Frozen construction

The implementation and all parameters are copied unchanged from Experiment
149. The held context is Orion HCT116-thinned. Its own perturbation donor is
excluded; predictions use K562-GWPS plus Orion-HEK293T, fixed PerturbMatch,
the essential-pruned held-HCT operator, and the solved a15 CD4 residual.

Only the harness's contiguous 18,400 non-targeting rows may be read during
construction. HCT116 perturbation expression remains sealed until scoring.

## Frozen score-level gate

- exact per-target PDS equality;
- projected six-member average delta at least `+0.004`;
- target-bootstrap 95% lower bound above zero;
- at least two of NMAE, reach, fidelity, and Jaccard have nonnegative scaled
  deltas;
- no one of those four DE metrics has scaled delta below `-0.015`;
- scaled MSE does not decrease.

Raw MSE remains a diagnostic. When both arms are above the official baseline,
the competition clamps both to scaled zero and raw movement cannot affect the
leaderboard objective.

HCT116-thinned retains some atlas/batch limitations, so a pass is supporting
evidence for a production build, not direct authorization to submit.

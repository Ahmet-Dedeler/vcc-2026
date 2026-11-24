# Experiment 110 implementation scaffold

Date: 2026-08-30

This scaffold implements the frozen `README.md` and
`preregistered_protocol.json` without reading source expression, STATE output,
or held truth. It does not train a checkpoint, construct a real candidate, or
authorize evaluation.

## Files

- `arbitration.py` contains pure NumPy/scikit-learn transformations, fitting,
  LOCO selection, target-block bootstrap, and permutation-null utilities.
- `artifact_contract.py` defines exact NPZ and manifest contracts for the
  development selection freeze and truth-free held prediction freeze.
- `test_arbitration.py` and `test_artifact_contract.py` use synthetic arrays
  only.

The frozen preregistration files were not changed.

## Numerical contract

The selector feature matrix has exactly nine columns in this order:

1. `a14_abs_lfc`
2. `state_abs_lfc`
3. `signed_log_magnitude_ratio`
4. `absolute_log_magnitude_ratio`
5. `a14_cross_donor_agreement`
6. `control_log_cpm`
7. `target_level_a14_state_cosine`
8. `state_control_bank_sign_agreement`
9. `state_control_bank_variance`

The builder accepts no gene, target, or context identifiers. Opaque target
blocks are carried separately for weighting, bootstrap, and within-target
permutation. Shared blocks remain shared when contexts are pooled, so every
perturbation receives one total unit of fit weight.

STATE must be supplied as four aligned control-bank predictions. The aggregate
is their mean, computed only after sign-agreement and variance features are
formed.

The fixed grid is the exact 144-member Cartesian product in the protocol.
Each LOCO rotation fits normalization only on its two training contexts. A
development selection is valid only when all three rotations have one unique,
identical winner.

Sparse replacement changes only threshold-selected disagreement coordinates.
The row cap scales that sparse delta by bisection until its PDS-space cosine
distance from a14 is within the selected limit. Unsupported coordinates remain
bit-identical to a14. The writer rejects a nonzero incumbent value on any
forbidden coordinate and requires the final candidate to be exactly zero
there.

The null path performs exactly 100 deterministic within-target label
permutations. Every repetition reruns all LOCO fits and all 144 candidates,
records the complete rotation score tensor, selects a macro winner, and refits
that null selector on all development contexts.

## Artifact contract

`development-context-freeze` stores one source context and independently
reconstructs its permitted features, STATE/reference labels, disagreement
coordinates, direction mask, and target blocks from the bound a14, four-bank
STATE, control, agreement, and public-reference arrays. Its manifest requires
the context-matched held-out STATE checkpoint claim and binds that checkpoint
manifest.

`development-selection-freeze` stores:

- the exact feature names, context order, and 144-row grid;
- all observed LOCO scores and the common winning index;
- frozen train-only mean, scale, coefficient, and intercept;
- all 100 null LOCO score tensors, selected configurations, label hashes, and
  refit model parameters.

`held-prediction-freeze` stores no truth, reference LFC, significance label, or
score. It stores:

- gene and target axes;
- a14, all four STATE banks, and their verified mean;
- support, forbidden, and geometry masks;
- disagreement-pair coordinates and permitted features;
- observed probabilities, sparse replacement mask, cap scales, and candidate;
- compact null pair selections and cap scales.

All manifests use exact key sets, bind the unchanged preregistration and
implementation hashes, bind every input path and SHA-256, reject path escape,
and refuse overwrite. A selection freeze must bind all three schema-valid
development contexts. A held freeze must bind the exact selection manifest.
`load_held_with_selection()` recomputes observed probabilities and all 100 null
masks from the frozen model parameters before accepting the pair.

## Lightweight verification

```bash
uv run python experiments/110-state-disagreement-arbitration/test_arbitration.py
uv run python experiments/110-state-disagreement-arbitration/test_artifact_contract.py
ruff check experiments/110-state-disagreement-arbitration
ruff format --check experiments/110-state-disagreement-arbitration
```

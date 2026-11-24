# 338 — fixed-rank moment freedom is too narrow: no integer construction

September 9, 2026. **DOES_NOT_EARN_INTEGER_CONSTRUCTION.**
Frozen at **31fad5a9**. Engineering relaxation only; no response accuracy
was evaluated and no real perturbed reference values were loaded.

## Frozen gate result

All 128 target/gene pairs solved both variable-depth, rank-constrained
extrema. All original-depth known answers and rank instruments passed.
However, the required fraction with at least 5% symmetric arithmetic-CPM
room was not met in three of four targets:

| target | genes with >=5% room in both directions | fraction | required >=25% |
|---|---:|---:|---|
| HAX1 | 4/32 | 12.5% | fail |
| SBDS | 6/32 | 18.75% | fail |
| PLAGL2 | 9/32 | 28.125% | pass |
| ZNF850 | 6/32 | 18.75% | fail |

The median symmetric relative room is **0.002679 (0.268%)** with fixed
rank sum, versus **0.200003 (20.000%)** without that constraint. The
constrained range is -0.080078 to +0.200024. A negative value means the
old arithmetic mean is outside the new feasible interval: changing the
row depths can force a one-directional shift even while preserving U.

These are continuous histogram bounds under the fixed .25 depth contrast,
baseline zero counts and declared finite count support. They are not a
universal ceiling on dual-moment methods. The gate stays failed; no support,
depth, threshold or target selection will be retuned from these outcomes.

## Verification

- Eight pure tests pass, including analytic extrema, saved histogram
  constraints, fractional-input rejection and constraint nesting.
- All 85 baseline per-gene sums and the whole emitted-count digest exactly
  reproduce experiment 323 before the LPs run.
- All **768 extrema** (original depths, variable depths with U, and variable
  depths without U) have successful solver status.
- Independent audit reconstructs all saved histogram equalities and
  primal/dual certificates without calling the construction's Phi, bounds,
  or LP solver.
- Maximum saved-histogram gene-total residual: **3.57e-9 counts**.
  Maximum scaled equality residual: **7.34e-9**.
  Maximum absolute scaled primal-dual gap: **2.61e-11**.
- Independent grouped-control U versus pooled ranks agrees within
  **1.40e-9**; Python-integer tie-corrected z versus the scorer kernel agrees
  within **1.78e-15**.
- All frozen gates and the terminal decision are independently reproduced.

The audit verifies the recorded relaxation and selected baseline blocks,
not a second full-baseline re-emission. It does not establish joint
integer row feasibility, unchanged tie-corrected z/BH calls, improved
biological signs, or a competition gain. All six competition members
remain unpriced for any hypothetical new emitter.

Run **16:43:18.753–16:45:31.280 UTC**, **132.53 seconds**.
Driver peak RSS 5,178,785,792 bytes; sampled tree peak **4.71 GiB**.
Run and independent audit exited successfully. The short audit finished
between memory samples; its 0.01 GiB guard sample is not its peak estimate.

- Inputs SHA:
  `45af8c049d6d71880c1a638c847683f5173b28bec1b28c9fa07805469cd9b966`
- Results SHA:
  `738a4bd183b91e6f0677c71651211ca85d99a86d010bbcef239395e44b5cc10c`
- Audit SHA:
  `358d22d2443e81c7aee5fd8b63595d6518f0ab56debb326000ee9b8e0f99c739`

No candidate, new prediction, production change, package or submission.

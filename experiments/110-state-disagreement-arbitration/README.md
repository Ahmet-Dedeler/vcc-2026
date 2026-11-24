# 110 — STATE/a14 disagreement arbitration

Status: **NO-GO before cross-fitted checkpoint training**

Date: **2026-08-30**

## Question

Can STATE improve the downstream sign/ranking metrics as a sparse specialist
without giving back a14's leading PDS geometry?

The public STATE baseline is complementary to a14: it loses PDS but gains MSE,
NMAE, fidelity, reach, and Jaccard. Experiment 106 tests only global
composition/effect/direction blends. This experiment instead learns whether
STATE or a14 is more likely to have the correct sign on a particular
target-by-gene pair, and changes only high-confidence disagreements.

## Candidate

For a14 LFC `a[p,g]` and row-norm-matched STATE LFC `s[p,g]`, fit

```text
q[p,g] = P(sign(s[p,g]) == sign(reference[p,g]) | prediction-only features)
```

only on coordinates where `sign(a) != sign(s)`. The candidate is

```text
h[p,g] = a[p,g] + lambda * I(q[p,g] >= tau) * (s[p,g] - a[p,g])
```

followed by the frozen row-level angular cap and exact incumbent fallback.

The selector may use only:

- `abs(a)` and `abs(s)`;
- signed and absolute a14/STATE magnitude ratios;
- a14 cross-donor agreement;
- control log-CPM;
- target-level a14/STATE cosine;
- STATE sign agreement and variance across four independently aggregated
  control-bank predictions.

Gene identity, target identity, literature values, STRING, regulons, and held
context labels are forbidden features.

## Why this is distinct

- Experiment 101 had no independent alternative sign source.
- Experiment 104 could rescale a14 but prohibited sign flips and new
  coordinates.
- Experiment 108 applied a lossy regulon projection directly and did not learn
  when to trust STATE over a14.
- Experiment 106 applies one global blend weight to every supported
  coordinate.

## Frozen development protocol

1. Use only the fixed step-20,000 `last.ckpt` from experiment 106. Never use
   `best.ckpt`.
2. Run STATE inference separately for the four frozen control banks; do not
   tune the checkpoint, control partition, or inference seed.
3. Build a14 analogues and strictly context-held STATE predictions for the
   already-local RPE1, Jurkat, and HepG2 public support contexts. The main
   experiment-106 checkpoint already holds out HepG2. Train two additional
   fixed-step checkpoints from the same pinned source/config with RPE1 and
   Jurkat respectively removed from training. Predictions from a checkpoint
   that trained on the evaluated context are forbidden.
4. Train only on a14/STATE sign-disagreement pairs whose reference LFC is
   defined and significant. Give every perturbation equal total weight.
5. Use weighted ridge logistic regression with no gene or target identifiers.
6. Select one common configuration by leave-one-context-out rotation over
   RPE1, Jurkat, and HepG2. Every rotation must prefer the same configuration.
7. Freeze model parameters, feature normalization, threshold, blend fraction,
   angular cap, axes, candidate arrays, and null predictions before opening
   the experiment-110 H1 or HEK293T evaluation.

The candidate grid is fixed in `preregistered_protocol.json`.

## Held evaluation

- H1-test: target-held confirmation.
- HEK293T: context-held confirmation.
- HCT116: PDS-only safety diagnostic; its absolute direction metrics remain
  invalid because of the documented depth/batch confounding.

Every gallery is complete for the frozen shared target axis. PDS excludes all
300 VCC panel coordinates. Directional credit excludes each row's own target
coordinate. Truth is raw and uncentered.

The H1/HEK293T data have been inspected by earlier hypotheses, so these are
hypothesis-specific confirmations rather than pristine untouched datasets.
The protocol is frozen before experiment 106 reveals STATE's performance on
those surfaces.

The two cross-fitted development checkpoints are expensive but necessary:
without them, RPE1/Jurkat labels would leak through the already-trained STATE
predictor and make the selector's apparent accuracy meaningless.

## Promotion gates

All must pass:

- oracle disagreement routing headroom at least `+0.015` directional credit
  on both H1 and HEK293T;
- the same configuration wins every development-context rotation;
- learned macro directional-credit delta at least `+0.006`, nonnegative on
  both primary surfaces, with target-bootstrap lower 95% above zero;
- top-20% weighted sign-precision gain at least `+0.02` macro and nonnegative
  on both primary surfaces;
- complete-gallery raw PDS delta no lower than `-0.001` on either primary
  surface and nonnegative in the macro;
- observed directional-credit gain exceeds all 100 frozen within-target label
  permutations;
- unsupported targets/genes are exactly a14 and all forbidden coordinates are
  zero.

Only after all effect-level gates pass may one sequential count-level run be
performed on H1 and HEK293T. It must improve NMAE on both, keep fidelity
precision and reach ordering nonnegative, and change total calls by at most
5%.

Failure means **NO-GO**. This experiment authorizes no challenge build,
submission, commit, or push.

## Result

Experiment 106 and the H1 domain check in experiment 113 closed this
`state_sm` arbitration family before the two expensive leave-context-out
checkpoints were trained.

The row-norm-matched disagreement oracle was `+0.14504` directional credit on
H1 but only `+0.00686` on HEK293T. The preregistered rule required at least
`+0.015` on **both** primary surfaces, so the family fails even with perfect
knowledge of which disagreements to accept.

An independent code audit also found that the scaffold used the
reference-significant mask both to fit and to define candidate application
pairs. Significance is unavailable at deployment, so that implementation
would have understated collateral changes on nonsignificant coordinates.
This did not affect the NO-GO decision: the oracle already failed under the
more favorable restricted population.

No RPE1-out or Jurkat-out checkpoint was trained, no count-level run was
performed, and no challenge artifact was built or submitted.

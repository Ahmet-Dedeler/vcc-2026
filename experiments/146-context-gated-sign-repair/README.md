# Experiment 146 — context-gated sign repair

Status: **revised implementation scaffold; source prepare v4 is not authorized**.

The first authorized prepare attempt on August 31, 2026 failed closed before
creating a package because six pseudobulk-complete response genes were absent
from the authenticated FDR axis. The temporary directory and lock were removed.
The revised scaffold freezes the deterministic FDR-axis intersection described
below; the v1 authorization is retained only as an audit record and cannot bind
the revised code.

The v2 attempt also failed closed before package creation: the 4,974 count and
six exclusions passed, but the preregistered diagnostic hash used a
noncanonical separator. V3 uses the runner's own `sha256_lines()` definition
with real newline separators. The v2 authorization is likewise retained only
as an immutable audit record.

The v3 attempt passed the corrected axis checks, then failed closed because the
published FDR matrix has 558,367 nonfinite entries over the 1,289 available
context-target blocks, including 38 blocks with no finite FDR values. V4
freezes those exact missingness counts and the all-missing block-set hash.
Nonfinite FDR is unlabelled and never imputed; unavailable pseudobulk pairs
retain all-NaN FDR rows and are never copied from the wider FDR table.

## Question

Can a nonlinear model use production-available donor evidence and a destination
control profile to identify downstream coordinates whose selected-slice
a15-core sign is wrong, then repair only a very small, energy-bounded subset
without giving back the existing retrieval advantage?

The intervention is deliberately narrow. It never learns a replacement effect
vector, changes an effect magnitude, introduces a new nonzero coordinate, or
fits an operator on the prediction panel. For a selected coordinate:

```text
d'[p,g] = -d[p,g]
```

Every unselected coordinate is bit-identical to the incumbent. Since a sign
flip preserves `abs(d[p,g])`, the row norm is preserved before reconstruction.
The selected coordinates must also satisfy a per-row squared-energy budget.

This is distinct from:

- Experiment 100, which learned confidence but did not flip signs;
- Experiment 101, whose sparse flip model used the earlier fixed
  STRING/context feature family and failed to identify a >60% wrong-sign tail;
- Experiment 104, which was a sign-preserving linear-credit adapter;
- Experiment 111, which fitted a low-rank linear context-program operator; and
- donor weighting, which changes whole donor vectors rather than routing signs
  per target and response gene.

## Actual source interface

The authenticated DepMap 26Q1 extraction contains:

- 16 cell-line contexts in six experimental batches;
- 90 experimental perturbation targets before isolation filtering;
- exactly 81 permitted evaluation targets after excluding the nine targets
  overlapping sealed 410 rows:
  `COX5B, METAP1, MTPAP, NRBP1, POLR2D, RAB18, SLC39A9, TRAIP, UBE2D3`;
- ten olfactory-control perturbations used to construct context control
  profiles;
- 22,966 response rows in the published LFC/FDR matrices and 34,647 response
  rows in the pseudobulk matrix; and
- zero perturbation-target overlap with the 300 VCC targets.

The 90-target source axis, the exact nine-target exclusion, and the resulting
81-target axis are SHA-256 frozen. The parser and authenticated schema are
reused from Experiment 111. Perturbation and control pseudobulks are
reconstructed at `target_sum=50,000` before `log1p` differencing. The response
axis is Experiment 142's frozen 4,992-gene axis intersected first with genes
complete across every permitted DepMap pseudobulk column and then with the
authenticated FDR axis, after excluding all 300 VCC targets. The final axis has
4,974 genes with SHA-256
`fefbe82822b0cdbcc56f49231520e5f52b8bcbd7f400faa5aafdb0602fd8aa79`.
The six pseudobulk-complete genes absent from the FDR table are frozen as
`FOXN3, PIGH, CLN6, PIF1, FBXO33, SPTSSA`.

Metric masks are separate and mirror the competition:

- Source-panel PDS excludes the union of all 81 permitted DepMap target
  coordinates and all 300 VCC target coordinates.
- Signed-credit labels, intervention eligibility, and evaluation exclude all
  300 VCC coordinates plus only the current row's own perturbation target.

The implementation asserts that these masks are not accidentally reused for
one another.

The source baseline is an a15-core, selected-slice reconstruction on the 81
permitted DepMap targets:

- equal K562-GWPS, Orion-HCT116, and Orion-HEK293T production donors;
- `alpha=0.982`;
- a gamma-0.65 signal-covariance whitener fitted from only the selected donor
  slices; and
- exclusion of all 81 evaluation-target coordinates and all 300 VCC target
  coordinates before fitting or applying the whitener.

PerturbMatch and CD4 are deliberately absent from this cheap source stage. A
source pass can only advance to a separate production-integration audit.

Complete donor `delta.npy` containers are never hashed or loaded. The prepare
authorization binds each donor's `genes.txt`, `perts.txt`, and `meta.json`,
plus the frozen Experiment 140 source ledger. After the source-access lock is
held, each `delta.npy` is opened read-only with `mmap_mode="r"` and only the
81-target by response-gene slice is copied. The nine sealed-overlap target
rows are excluded before every centering, scale, covariance, and support
statistic.

Both Orion gene axes contain duplicate symbols. Their full source shape is
preserved, and symbol lookup resolves to the last occurrence, matching the
production donor loader. Five permitted targets (`AURKB`, `KIF11`, `MIPEP`,
`MYC`, and `TYMS`) have K562-GWPS support only. They remain valid baseline and
feature rows; the intervention still requires at least two donors and cannot
alter them.

No challenge control `.h5ad` is opened. Only `gene_names.csv`,
`pert_counts.csv`, and `manifest.json` may be read from `data/challenge`.

## Strict two-axis cross-fitting

Context-only LODO is insufficient, and contexts sharing an experimental batch
are not independent. Experiment 146 therefore uses exactly 30 outer cells:
`six leave-one-experimental-batch-out folds × five whole-target folds`.

The frozen 81 targets are deterministically sorted by
`(sha256(lowercase target), target)` and assigned round-robin to five folds
with sizes `17, 16, 16, 16, 16`.

For outer experimental batch `b` and outer target fold `f`:

1. All outcomes from every context in batch `b` are inaccessible.
2. Every occurrence of every target in fold `f` is absent from fitting,
   outcome-dependent weighting or normalization, calibration, and candidate
   selection.
3. On the other four target folds, inner whole-target cross-fitting holds out
   one fold at a time, fits on the remaining three, and predicts the inner
   target fold across contexts outside batch `b`.
4. Those four prediction blocks are concatenated before any arm is scored
   against the complete 64/65-target inner gallery. Because the authenticated
   source has seven missing context-target pairs, candidate selection uses only
   non-held contexts with all 64/65 inner targets observed; incomplete contexts
   may contribute their available rows to fitting but never a shortened
   selection gallery. Mini-gallery scoring is forbidden.
5. There is no hyperparameter grid. The single fixed full-feature HGB arm is
   enabled only by its complete-gallery inner result; otherwise the fold uses
   exact identity.
6. The model is refit on contexts outside batch `b` and targets outside `f`.
7. It predicts fold `f` in every context of held batch `b` exactly once.

After all five target folds are complete, the cross-fitted rows form one
prediction per context. PDS uses every observed target in that held context
(79–81 targets because of the seven authenticated missing pairs) after
removing all 81 source-panel target coordinates. Thus every reported target
and batch were absent everywhere their outcomes could affect fitting or
intervention selection.

Control PCA is target-independent and uses controls only. It is nevertheless
fit only on contexts outside the held experimental batch. Donor effects for
held targets remain available because they are exactly the information
available at production time.

## Fixed estimator and features

The primary estimator is
`HistGradientBoostingClassifier(max_depth=3, learning_rate=0.05,
max_iter=100, l2_regularization=1, early_stopping=False)`.
Each fit uses at most 500,000 deterministically sampled pair rows, with equal
total weight per context-target block and inverse-class balancing.
Fits are sequential; the runner rejects overlapping fits.

Features are fixed before source evaluation:

1. incumbent sign, log magnitude, and within-target magnitude rank;
2. real-donor mean, median, MAD, median absolute effect, positive fraction,
   support count, and agreement with the incumbent;
3. control-similarity-weighted donor effect and sign vote;
4. nearest-control donor effect;
5. destination and donor-weighted basal expression for the response gene;
6. their difference, mean donor similarity, and maximum donor similarity;
7. target-level incumbent norm and donor pairwise cosine agreement; and
8. eight source-control PCs plus eight destination-minus-donor PC
   coordinates.

The label is whether the incumbent and source truth signs disagree among genes
with authenticated source FDR at most `0.05`. Rows with zero/nonfinite
incumbent or truth are unavailable. Each label block excludes the row's own
perturbation target coordinate and all 300 VCC coordinates. The marginal
sample weight is `abs(incumbent)` normalized to equal total mass per
context-target block, matching continuous signed-credit sensitivity.

The primary arm and five matched controls use identical outer and inner splits:

- full-feature HGB, the only promotable arm;
- matched zero-tail control using the full-HGB arm's exact selected mask and
  budgets;
- deterministic donor-agreement scoring;
- full-feature weighted logistic regression;
- donor-only weighted logistic regression; and
- same-capacity donor-only HGB (`max_iter=100`).

## Single fixed intervention

There is exactly one stage-1 intervention policy:

```text
score quantile                 = 0.995
total squared-energy cap       = 0.0003
PDS-axis squared-energy cap    = 0.0003
selected coordinate count cap = 64
signed-L1-mass cap             = 0.005
```

Eligible coordinates require at least two production donors. Selection is
descending score with coordinate index as the deterministic tie break. Each
cap is enforced while coordinates are added. There is no 3×3 grid and no
outcome-driven threshold calibration.

## Source metrics and falsifiers

- Complete available-gallery raw cosine PDS on each held context (79–81
  observed targets) after removing all 81 permitted DepMap target coordinates
  and all 300 VCC coordinates.
- Continuous signed credit:

  ```text
  sum(sign(truth) * prediction) / sum(abs(prediction))
  ```

  computed at frozen source FDR `0.05` after removing the row's own target and
  all 300 VCC coordinates.

- `abs(incumbent)`-weighted wrong-sign precision among flipped coordinates.
- Five deterministic held-target-block summaries.
- Twenty-thousand experimental-batch bootstrap draws.
- Twenty-thousand whole-target-block bootstrap draws, with every target first
  aggregated using the same equal-batch estimand as the principal score.
- Nine whole-target refit permutation nulls. Within each outer cell, targets
  are stratified by their exact availability pattern across non-held contexts.
  Each null deranges complete outcome blocks only within identical-pattern
  classes, leaves singleton classes fixed and records them, refits the same
  fixed full HGB, and evaluates the unchanged held batch × held target-fold
  cell.
- Frozen FDR sensitivity at `0.01`, `0.05`, and `0.10`.

## Source promotion gates

All are mandatory:

1. Batch-bootstrap 95% lower bound for macro signed-credit gain is above zero.
2. Target-block-bootstrap 95% lower bound is above zero.
3. At least 12 of 16 held contexts improve in signed credit.
4. Every held target fold has nonnegative macro signed-credit gain.
5. Weighted flip precision is at least 65% in at least 12 contexts.
6. Macro raw PDS delta is at least `-0.002`, with no context below `-0.010`.
7. Full HGB beats zero-tail and same-capacity donor-only HGB with positive
   batch-bootstrap lower bounds; it must also beat the remaining matched
   controls.
8. The observed gain exceeds the p95 whole-target refit-permutation null.
9. Macro gain is positive at each frozen FDR threshold.
10. No context supplies more than 30% of positive signed-credit gain.
11. All source/sealed exclusions, batch and target isolation, norm identity,
    total-energy, PDS-energy, selected-count, signed-L1, resource, and manifest
    checks pass exactly.

A failed gate closes this mechanism. Thresholds, folds, features, and seeds
must not be changed after outcomes are observed.

## Immutable stages

No authorization file is part of this scaffold.

When separately authorized, source work has two commands:

```bash
uv run python experiments/146-context-gated-sign-repair/run.py prepare \
  --confirm I_CONFIRM_EXPERIMENT_146_SOURCE_PREPARE_V4

uv run python experiments/146-context-gated-sign-repair/run.py evaluate \
  --confirm I_CONFIRM_EXPERIMENT_146_SOURCE_EVALUATE_V1
```

`prepare` writes one atomic, no-overwrite source package under
`work/source-prepared-v1/`. `evaluate` reads only that validated package and
writes one atomic, no-overwrite result under `work/source-evaluation-v1/`.
Every file is size- and SHA-256-bound by an artifact manifest.
Both stage authorizations bind `scripts/run_with_memory_guard.py`; evaluation
loads and validates the prepared arrays only after the process-tree monitor is
active.

The commands require separate, externally created authorization documents:

```text
SOURCE_ONLY_PREPARE_AUTHORIZATION.v4.json
SOURCE_ONLY_EVALUATE_AUTHORIZATION.v1.json
```

The v4 prepare authorization and the evaluation authorization are intentionally
absent. The retained v1/v2/v3 prepare authorizations bind earlier code and are
invalid for this revision. There is no development, sealed-outcome,
count-generation, scorer, build, or submission command in `run.py`.

## Lightweight checks

These checks are synthetic/static only:

```bash
uv run python -m unittest \
  experiments/146-context-gated-sign-repair/test_run.py

uv run ruff check \
  experiments/146-context-gated-sign-repair/run.py \
  experiments/146-context-gated-sign-repair/test_run.py

uv run ruff format --check \
  experiments/146-context-gated-sign-repair/run.py \
  experiments/146-context-gated-sign-repair/test_run.py

python3 -m json.tool \
  experiments/146-context-gated-sign-repair/preregistered_protocol.json \
  >/dev/null

uv run python experiments/146-context-gated-sign-repair/run.py audit
```

Both authorized stages must run under
`scripts/run_with_memory_guard.py` with an 18 GiB process-tree RSS ceiling and
a 25% system-free-memory floor. The monitor report is persisted inside the
immutable stage and validated before downstream use. Fits run one at a time.

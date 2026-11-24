# 051 — CD4 breadth-adaptive blend

## Goal

Test whether primary-CD4 effect breadth can decide **how much** CD4 response to
blend into the existing somatic-donor consensus for each perturbation.

The experiment compares:

- fixed CD4 fractions;
- smooth target-adaptive fractions derived only from compact-CD4 effect
  breadth and QC;
- the selected adaptive rule against a mean-matched constant control; and
- the selected adaptive rule against deterministic shuffled-breadth controls.

It is an effect-level PDS experiment. It does not run `cell-eval2`, emit
single-cell counts, alter production code, build a submission, or submit
anything.

## Leak-free split

The split and order are hard-coded:

1. Sweep hyperparameters on the HCT116 and HEK293T held-donor folds.
2. Select the best fixed rule and best adaptive rule using only those two
   folds.
3. Freeze both rules and the overall winner.
4. Evaluate K562 once as the lineage-analogue holdout, together with only the
   predeclared matched-constant and shuffled-breadth controls.

K562 PDS never appears in the development rows or selection functions. K562
may act as an input donor while HCT116 or HEK293T is held out, as in ordinary
LODO; it is never used as a target during hyperparameter selection.

The code never reads `K562_comparison.suppl_table.csv`, leaderboard results, or
challenge outcomes. Adaptive weights use only each compact CD4 condition's
`obs.csv`.

## Adaptive rule

For each target, the confidence inputs are:

- `log1p(n_total_de_genes)`, falling back to `n_downstream` or
  `n_up_genes + n_down_genes`; and
- a QC composite from available on-target significance, inverted low
  expression/off-target flags, guide count, and target-cell count.

For Stim8hr, the rule uses Stim8hr metadata. For the three-condition CD4 mean,
it averages metadata across Rest, Stim8hr, and Stim48hr.

Breadth and QC are robustly standardized. The unnormalized confidence is

```text
sigmoid(slope * (breadth_z + qc_strength * qc_z))
```

It is then monotonically rescaled and clipped so its mean is exactly the
candidate's fixed fraction. Therefore every adaptive candidate has a
mean-matched constant control. The shuffled control permutes breadth across
target labels while leaving the QC labels and the full weight-generating rule
unchanged.

Predeclared grid:

- fixed mean fraction: `0.10, 0.20, 1/3, 0.50, 0.67`;
- adaptive slope: `0.5, 1.0, 2.0`;
- QC strength: `0.0, 0.5, 1.0`;
- CD4 source: Stim8hr or the equal three-condition mean;
- prediction mode: centered raw or centered mean-absolute gene-normalized.

Within `0.001` development macro PDS, selection prefers the configuration with
the stronger worst fold and then the simpler raw/no-whitening/Stim8hr rule.

## Gene and perturbation axes

- Perturbations are the VCC-panel intersection across all three somatic banks
  and all three compact CD4 conditions.
- The master feature axis is the official VCC gene order represented by any
  bank, after removing **all 300 VCC target genes before centering,
  normalization, blending, or whitening**.
- Each fold is scored only on genes measured by its held donor.
- Missing source genes are zero-filled on that held-donor score axis.
- Truth is always the held donor's centered raw effect.

All source `delta.npy` files and all prepared fold matrices are memory mapped.
Prepared matrices live in a temporary directory and are deleted after the
run.

## Required guarded run

Do not invoke the real-data CLI directly. It refuses to run unless
`--guarded` is present. The canonical command is:

```bash
uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib 4 --min-free-percent 15 -- \
  uv run python experiments/051-cd4-breadth-blend/run.py \
    --guarded \
    --cd4-root /path/to/cd4_compact
```

The evaluator also has an internal `3.75 GiB` stop. BLAS is limited to two
threads. No `cell-eval2` process is started.

Outputs, when real compact data exists:

```text
experiments/051-cd4-breadth-blend/results.json
experiments/051-cd4-breadth-blend/development_rows.csv
experiments/051-cd4-breadth-blend/k562_holdout_rows.csv
```

## Optional production-union whitening

Exact 17k-gene union-whitener fitting is not safe to perform locally under the
4 GiB cap. The evaluator can apply cloud-fitted, fold-specific operators via
`--whitener-dir`, but it rejects an operator unless its provenance proves:

- the held donor was not a fit source;
- challenge outcomes were not used;
- the source list is exactly the two non-held somatic donors;
- all 300 panel targets were excluded; and
- the operator gene list exactly equals the non-held donors' production union
  on that fold's score axis.

Expected layout:

```text
whitener-root/
  k562_gwps/
    meta.json
    genes.txt
    evals.npy
    evecs.npy
  orion_hct116/
    ...
  orion_hek293t/
    ...
```

`meta.json` contract:

```json
{
  "schema": "vcc2026-exp051-fold-union-whitener-v1",
  "held_out_donor": "orion_hct116",
  "fit_uses_held_truth": false,
  "fit_uses_challenge_outcomes": false,
  "axis_kind": "production_union",
  "source_donors": ["k562_gwps", "orion_hek293t"],
  "excluded_targets_sha256": "<sha256 of panel symbols joined by newline>",
  "operator_genes_sha256": "<sha256 of genes.txt joined by newline>"
}
```

The default run omits whitening. This keeps the critical CD4 blend test
available immediately when compact effects arrive without pretending that an
all-donor production operator is a valid LODO operator.

## Tests

Synthetic tests write only to temporary directories under this experiment:

```bash
uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib 4 --min-free-percent 15 -- \
  uv run python experiments/051-cd4-breadth-blend/test_run.py
```

They verify monotone mean-matched adaptive weights, the HCT116+HEK293T/K562
split, all-panel-coordinate exclusion, shuffled controls, memory-mapped
inputs, rejection of held-truth whitening, and refusal of unguarded CLI runs.

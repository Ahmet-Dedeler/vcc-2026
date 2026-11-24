# 118 — PDS-null uncertainty-weighted cross-panel residual

Status: **source-only first stage**. This directory can read only the three
hash-pinned CD4 compact conditions, their SE/FDR/QC artifacts, and the canonical
VCC gene and target axes. It cannot read Pan values, held harness truth, Nadig
or KOLF outcomes, challenge truth, counts, submissions, or `cell-eval2`.

## Question

Can a small cross-condition CD4 residual improve downstream gene direction
while being mathematically invisible to VCC PDS?

VCC PDS removes all 300 panel-gene coordinates before cosine retrieval. This
experiment therefore places a residual only on off-diagonal panel coordinates,
sets the perturbed target coordinate to exact zero, and sets every non-panel
coordinate to exact zero. Every row is centered over its mutable supported
panel coordinates and repaired in one deterministic coordinate so its
correctly rounded compensated float64 panel sum is exactly zero.

This first stage does not combine the residual with an incumbent model and
does not evaluate any held surface. It asks only whether the residual transfers
directionally across the three source CD4 conditions.

## Frozen source data

The only biological values allowed are:

- `cd4_rest`
- `cd4_stim8hr`
- `cd4_stim48hr`

For each condition the runner reads `log_fc.npy`, `lfc_se.npy`,
`adj_p_value.npy`, `valid_log_fc.npy`, `base_mean.npy`, `obs.csv`,
`genes.txt`, and `perts.txt`. It also reads the canonical
`data/challenge/gene_names.csv` and `data/challenge/pert_counts.csv` axes.
Every path and SHA-256 is pinned in `SOURCE_ONLY_AUTHORIZATION.json`.

A condition-level target passes QC only when:

- on-target knockdown is significant;
- no neighboring-gene knockdown is flagged;
- no distal off-target is flagged;
- target expression is not flagged low;
- at least two guides are present.

A source estimate is usable when target QC passes, `valid_log_fc` is true,
LFC/SE/FDR are finite, and SE is positive. Held-condition directional truth
additionally requires FDR < 0.05 and nonzero LFC.

## Candidates

The exact grid is four estimators times three doses, plus identity:

```text
estimator ∈ {unweighted, inverse_variance, r_mu, r2_mu}
eta       ∈ {0.25, 0.50, 1.00}
identity  = exact zero
```

For valid training-condition estimates `l_c` with standard errors `s_c`:

```text
mu = sum(l_c / s_c^2) / sum(1 / s_c^2)
z  = abs(mu) * sqrt(sum(1 / s_c^2))
r  = (2 Phi(z) - 1)
     * abs(sum((1 / s_c^2) * sign(l_c))) / sum(1 / s_c^2)
```

The estimators are the unweighted mean, `mu`, `r * mu`, and `r^2 * mu`;
each is multiplied by `eta`. A target-response pair needs at least two valid
training conditions. Only off-diagonal panel responses are mutable.

## Source-only selection

There are three leave-one-CD4-condition-out folds. Whole target rows stay
together. For each held condition and target:

1. use only held-condition off-diagonal genes with valid, FDR-significant,
   nonzero truth;
2. compute signed linear credit
   `sign(truth) @ prediction / sum(abs(truth))`;
3. rank those genes by `abs(prediction)` and compute truth-magnitude-weighted
   sign precision in the top 20%.

Metrics are averaged equally over evaluable targets in each condition.
A nonidentity candidate is eligible only if every fold has credit >= 0 and
top-20% precision >= 0.5. Eligible candidates are ordered by:

1. larger worst-fold signed linear credit;
2. larger worst-fold top-20% sign precision;
3. smaller `eta`;
4. simpler estimator (`unweighted`, `inverse_variance`, `r_mu`, `r2_mu`).

Exact ties prefer identity.

## Frozen nulls

Before any effect LFC/SE/FDR array is loaded, `build` freezes:

- 100 whole-target label derangements, shared across conditions;
- 100 response-label derangements within ten deterministic pooled-CD4
  expression bins.

Expression bins are formed by sorting source-observed panel genes by the
median of `log1p(base_mean)` over all rows and conditions, breaking ties by
canonical panel order, then splitting the sorted genes into ten balanced
bins. Panel genes absent from the compact source remain fixed.

Each null repeats the complete 13-arm leave-one-condition-out selection while
leaving held-condition truth labels fixed. A selected nonidentity candidate
survives the source-only gates only if its observed worst-fold credit is
strictly above the 95th percentile of both complete-selection null
distributions.

## Commands

All commands are intentionally local and source-only:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -B \
  experiments/118-pds-null-panel-residual/run.py audit

PYTHONDONTWRITEBYTECODE=1 python3 -B -m unittest \
  experiments/118-pds-null-panel-residual/test_run.py

PYTHONDONTWRITEBYTECODE=1 python3 -B \
  experiments/118-pds-null-panel-residual/run.py build

PYTHONDONTWRITEBYTECODE=1 python3 -B \
  experiments/118-pds-null-panel-residual/run.py verify
```

`build` refuses to overwrite `work/source-freeze-v1/`. The output contains
the 13 final all-condition panel candidates, the fail-closed selected full
residual, both frozen permutation banks, fold and null results, and a
hash-bound freeze manifest. It never generates cells, counts, or a submission.

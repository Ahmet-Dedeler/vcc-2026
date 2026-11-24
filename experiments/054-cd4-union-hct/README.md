# 054 — CD4 mean3 under the clean HCT union whitener

## Goal

On the real 132-target Orion-HCT116 held-out panel, compare:

1. the current K562-gwps + Orion-HEK293T equal-source baseline after the
   experiment-040 clean-HCT union whitener; and
2. the predeclared fixed `0.50` blend with the primary-CD4
   Rest/Stim8hr/Stim48hr mean, blended before that same operator.

This is a frozen comparison, not an HCT hyperparameter search. The primary
candidate is:

```text
gene_norm + cd4_mean3 + fraction 0.50
+ preserve_baseline missing-gene handling
+ production_self whitening renormalization
```

The harness also reports the raw-CD4 version, one required structural control
for missing CD4 genes, and a diagnostic fixed-scalar whitening policy.

## What is and is not measured

The reference side is computed directly from the real HCT H5AD:

```text
log1p(50,000 * perturbation_pseudobulk)
  - log1p(50,000 * HCT_control_pseudobulk)
```

The prediction side remains an effect matrix. PDS excludes all 300 challenge
target coordinates and uses the full HCT target set as the retrieval gallery.
The covered-target and all-three-condition summaries subset query rows only;
they do not make an easier, smaller retrieval gallery.

The experiment:

- does not run `cell-eval2`;
- does not emit single-cell counts;
- does not fit or modify a whitener;
- does not build or submit a challenge artifact; and
- does not select a configuration on HCT.

## Exact baseline and whitener contract

The source baseline reproduces experiment 031's relevant operations:

- center each donor over all of its perturbations;
- measure each donor's median centered-effect norm on every mapped official
  gene, including duplicate Orion columns as the existing code did;
- map duplicate symbols by the last source occurrence, matching the existing
  repeated-index assignment;
- normalize each available donor target effect to the fixed `6.817` source
  scale; and
- average the available K562 and HEK293T source votes equally.

The only accepted operator is:

```text
experiments/040-union-gamma/out_hct_clean/whitener_union.npz
experiments/040-union-gamma/out_hct_clean/whitener_union.json
```

`run.py` independently derives the union `wcols` from K562 + HEK293T on the
official 18,533-gene HCT axis, removes the full 300-target challenge panel, and
requires the adjacent JSON to equal the experiment-040 metadata byte-for-byte
as parsed:

```json
{
  "axis": "union",
  "donors": ["k562_gwps", "orion_hek293t"],
  "n_cols": 17942,
  "wcols_sha256": "b2cae0d381c3cd7065a3de33f1a6da5514577f2b3187b6dc84786ee789e27e4d",
  "excluded_targets": ["<the sorted 300-gene challenge panel>"]
}
```

The current checkout's independently derived hashes are:

```text
challenge panel in CSV order:
c9c4c9a69f76afd4507e9a7619fe9925c19ff34e5878c7854493683b23bea5ca

clean-HCT wcols int64 bytes:
b2cae0d381c3cd7065a3de33f1a6da5514577f2b3187b6dc84786ee789e27e4d
```

If the clean artifact or metadata is missing or differs, the harness stops. It
never falls back to the quarantined
`experiments/040-union-gamma/out/whitener_union.npz` and contains no covariance
fit path.

The 1.23 GB `evecs.npy` member is streamed out of the NPZ into scratch and
memory-mapped. Float32 eigenvector blocks are promoted to float64 in bounded
chunks, matching `Whitener.load` numerically without allocating a full
2.46 GB float64 copy.

Two whitening-renormalization policies are evaluated:

- `production_self` (**primary**): reproduce
  `scripts/make_submission.py --whiten-renorm preserve`; each candidate
  computes and applies its own
  `before_operator_median_norm / after_operator_median_norm` scalar.
- `baseline_fixed` (diagnostic): reuse the baseline's scalar for every
  candidate, isolating the linear blend/operator effect from the panel-level
  renormalization change.

The dense operator is identical in both policies. A CD4-absent target is
identical to baseline before whitening under both. It remains exactly
identical after whitening only under `baseline_fixed`; under
`production_self`, the candidate-specific panel scalar can move that row. The
induced maximum absolute difference is reported rather than treated as an
error.

## CD4 construction

For each of `cd4_rest`, `cd4_stim8hr`, and `cd4_stim48hr`:

1. map compact genes to the official HCT axis;
2. remove all 300 challenge target coordinates;
3. center every gene on all perturbations available in that CD4 condition;
4. optionally divide by the source-only mean absolute centered effect, with
   the established 25th-percentile floor;
5. restore the normalized source's median effect norm; and
6. scale the condition to the same fixed `6.817` source norm as the current
   donor baseline.

For an HCT target present in fewer than three CD4 conditions, `mean3` means the
mean of the available condition votes. A target absent from all three
conditions is always copied exactly from baseline.

The missing-gene policies are intentionally separate:

- `preserve_baseline`: blend only where at least one available CD4 condition
  measured the gene; leave baseline unchanged elsewhere.
- `zero_vote`: map unavailable CD4 genes to zero and apply the global 50%
  convex blend, halving baseline effects at those coordinates.

PDS is reported for:

- all HCT targets;
- targets covered by at least one CD4 condition; and
- targets shared by all three CD4 conditions.

The result grid is fixed at eight comparisons:

```text
2 CD4 modes
x 2 missing-gene policies
x 2 whitening-renormalization policies
```

No value in this grid is selected on HCT. The predeclared primary remains
`gene_norm / preserve_baseline / production_self`.

## Exact cloud paths and command

Assumptions on the clean-HCT worker:

```text
/work/repo
/work/repo/data/harness/orion-hct116.h5ad
/work/repo/data/harness/orion-hct116.h5ad.meta.json
/work/repo/data/donors/k562_gwps/
/work/repo/data/donors/orion_hek293t/
/work/repo/data/external/cd4t/cloud-compact/cd4_compact/
/work/repo/experiments/040-union-gamma/out_hct_clean/whitener_union.npz
/work/repo/experiments/040-union-gamma/out_hct_clean/whitener_union.json
```

The scratch filesystem needs at least 1.3 GB free for the extracted
eigenvector member. From `/work/repo`, run exactly one process:

```bash
mkdir -p /work/vcc054-scratch

uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib 8 --min-free-percent 10 --poll-seconds 2 -- \
  uv run python experiments/054-cd4-union-hct/run.py \
    --guarded \
    --max-rss-gib 7.5 \
    --scratch-parent /work/vcc054-scratch \
    --output experiments/054-cd4-union-hct/results.json
```

Expected outputs:

```text
experiments/054-cd4-union-hct/results.json
experiments/054-cd4-union-hct/summary.csv
experiments/054-cd4-union-hct/per_target_pds.csv
```

If the worker does not already contain the downloaded compact CD4 directory,
copy the preserved 187.8 MiB compact archive there and extract it so the final
path is exactly:

```text
/work/repo/data/external/cd4t/cloud-compact/cd4_compact/
```

Do not start this run until experiment 040 has finished writing both
`whitener_union.npz` and `whitener_union.json`.

## Synthetic validation

The tests touch no real data and are intended to stay below 1 GiB:

```bash
uv run python scripts/run_with_memory_guard.py \
  --max-rss-gib 1 --min-free-percent 15 --poll-seconds 1 -- \
  uv run python experiments/054-cd4-union-hct/test_run.py
```

They verify:

- the frozen 0.50 protocol and primary candidate;
- the three requested target subsets;
- exact absent-target preservation before whitening under both policies and
  after whitening under `baseline_fixed`;
- reporting, without rejection, of the absent-target change induced by
  production's candidate-specific preserve scalar;
- the structural difference between zero-vote and preserve-baseline
  missing-gene behavior;
- strict rejection of whitener metadata drift;
- direct numerical agreement with experiment 031's `MultiDonorModel`,
  including duplicate-source handling;
- memory-mapped NPZ extraction; and
- refusal to run the real-data CLI without the external-guard acknowledgement.

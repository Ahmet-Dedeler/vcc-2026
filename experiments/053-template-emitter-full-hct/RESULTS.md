# 053 — Full-HCT streaming template emitter

Date: 2026-08-22

## Full-HCT decision

**Reject fixed-margin template emission. Keep v9's quasi emitter.**

The production-shaped HCT116 run preserved PDS exactly, as required by the
column-total invariant, but lost `0.070689` on the reference-scaled six-metric
average:

| metric | quasi raw | template raw | favorable raw change | scaled change |
|---|---:|---:|---:|---:|
| PDS cosine | 0.7978834143 | 0.7978834143 | 0 | 0 |
| MSE normalized | 2.0537529303 | 1.7966571507 | +0.2570957796 | 0 |
| direction fidelity | 0.5281190914 | 0.3998806920 | -0.1282383994 | -0.4274613315 |
| direction reach | 0.0150248423 | 0.0347673651 | +0.0197425228 | +0.0220340657 |
| DE Jaccard | 0.0337623585 | 0.0392606798 | +0.0054983213 | +0.0148603279 |
| LFC NMAE | 0.9570910100 | 0.9772758242 | -0.0201848141 | -0.0335686249 |
| **six-metric average** | **0.1216402396** | **0.0509509792** | — | **-0.0706892605** |

The realistic control-derived support pattern improves Jaccard and reach, but
its severe fidelity loss dominates the score. The lower raw MSE remains above
the official baseline and therefore earns no scaled credit. This closes the
template-emission path for the current model.

The successful run exited `0`. Compact outputs are preserved under
`cloud-output/g0.65_q1.0_a1.08_spread0.25/`; the template aggregate SHA-256 is
`d029c8aec95c24ff2de2b357875d8b5afc85eef0fecaa1d64f71b3e469229d72`.
No challenge submission was made.

## Outcome

A cloud-ready transformer now takes one complete, valid experiment-040 HCT116 prediction and
re-emits every 400-cell perturbation block with the fixed-margin control-template method.

The implementation:

- reads the prediction and real HCT116 files in backed mode;
- reads only one perturbation or control chunk at a time;
- keeps only one 400 x gene-axis dense working block;
- preserves every perturbation's integer per-gene column totals exactly;
- asserts scorer pseudobulk equality before writing and after reopening the completed file;
- copies control rows unchanged;
- copies all non-`X` H5AD objects and root attributes from the source, preserving `obs`, `var`,
  ordering, categories, indices, and auxiliary metadata;
- streams CSR `data` and `indices` into growable HDF5 datasets;
- writes to a temporary sibling and atomically renames only after success;
- rejects split/non-contiguous perturbation blocks, mismatched HCT116 gene or target axes,
  non-integer counts, incorrect cell counts, empty cells, per-cell totals above 1,000,000, and
  global stored-entry counts above 4.75 billion.

No production script was changed. No `cell-eval2` command was run locally and no submission was
built or sent.

## Input contract

The transformer intentionally fails closed unless:

- prediction and real HCT116 files have identical gene names and order;
- prediction and real files contain the same perturbation set;
- every non-control prediction label occupies one contiguous 400-row block;
- a control block exists in the prediction;
- real HCT116 contains at least 400 control cells.

The experiment-040 H1 artifacts currently present locally do not satisfy the HCT116 gene/target
contract and are not usable as full experiment-053 inputs. The cloud runner expects the clean
HCT116 outputs under:

```text
experiments/040-union-gamma/out_hct_clean/runs/<gamma-tag>/pred.h5ad
```

## Streaming design

For each perturbation:

1. Read its 400 source rows as sparse CSR.
2. Validate raw integer counts and calculate exact integer column totals.
3. Rescale one fixed, library-quantile-stratified set of 400 real HCT116 controls to the block's
   grand total.
4. Allocate each gene across those templates according to its real control support and abundance.
5. Repair row totals using within-gene transfers, which cannot change column totals.
6. Assert exact column-total and scorer-pseudobulk equality.
7. Append the emitted CSR block to the output HDF5 file and release block memory.

Control rows are copied in 256-row sparse chunks. A second backed pass reopens the final file and
checks:

- `obs` equality;
- `var` equality;
- exact control counts;
- exact per-target column totals;
- bitwise-equal scorer bulk-log vectors.

Each target's column-total SHA-256 digest is recorded in the transform report.

## Comparison with quasi

| Property | Existing quasi, spread 0.25 | Fixed-margin template |
|---|---|---|
| Per-gene total | Exact | Exact and copied from quasi |
| PDS plug-in profile | Fixed | Bitwise identical to quasi |
| MSE plug-in profile | Fixed | Bitwise identical to quasi |
| Library-size shape | Nearly constant | HCT116-control shape, globally rescaled |
| Gene zero/support pattern | Synthetic even plus random scatter | Derived from real HCT116 controls |
| DE rank behavior | Strong over-calling in experiment 049 | Better Jaccard/sign/top-500 proxy in 3/3 tested targets |
| LFC NMAE proxy | Better in experiment 049 | Worse in 3/3 tested targets |
| Compute | Fast vectorized allocation | Slower gene-wise allocation plus row repair |
| Peak working matrix | One 400 x genes count block | One 400 x genes count block plus support mask |

Because the pseudobulk is exact, PDS cannot change. MSE's plug-in term also cannot change, but its
jackknife term can change with the cell layout. The intended gain is entirely through more
realistic ranks, zeros, and variance for the DE metrics.

## Validation performed

### Synthetic end-to-end test

A 16-cell, six-gene fixture with two four-cell perturbations and eight controls was transformed,
written, reopened, and checked.

- Output `obs`: exact
- Output `var`: exact
- Auxiliary `uns`: preserved
- Control count rows: exact
- Both perturbation column totals: exact
- Both scorer pseudobulk vectors: exact
- Dry-run mode: exercised emission without creating an output H5AD
- Score collector: paired higher/lower-is-better deltas verified
- Tests: **3 passed**
- Guarded test peak: **0.14 GiB**

### Real-size dry run

Because no clean full HCT116 prediction was present locally yet, the existing 52,400 x 18,080
H1-layout prediction and matching H1 test reference were used only to exercise the real-size
memory path. This is not a biological result and does not rehabilitate the quarantined H1 model
validation.

- Controls streamed: 18,400
- Templates: 400
- Exercised block: AASS, 400 x 18,080
- Exact pseudobulk: yes
- Script peak RSS: **0.804 GiB**
- Outer process-tree peak RSS: **0.81 GiB**
- Hard guard: **4.0 GiB**
- Output H5AD written: no

Full transformation remains block-serial, so target count affects runtime and output size rather
than peak block memory. The internal limit defaults to 3.8 GiB and the supplied cloud runner wraps
the process with a 4.0 GiB process-tree guard.

## Cloud scoring runner

`run_cloud_after_gamma.sh` expects experiment 040's clean HCT116 gamma outputs to already exist on
the worker. Set `GAMMA_TAG` to process exactly one predeclared gamma directory; when it is unset,
the script processes every gamma directory it finds. For each selected directory it:

1. transforms the quasi prediction under the 4 GiB guard;
2. reuses experiment 040's existing quasi `agg_results.csv` by default;
3. runs template `cell-eval2` scoring serially;
4. writes a paired raw-metric comparison;
5. optionally uploads compact reports and aggregate CSVs when `S3_OUTPUT_URI` is set.

To rescore quasi on the same worker instead of reusing its existing result, set
`SCORE_QUASI=1`. No scorer processes are parallelized.

Example:

```bash
REPO=/work/repo \
GAMMA_ROOT=/work/repo/experiments/040-union-gamma/out_hct_clean \
GAMMA_TAG=g0.65_q1.0_a1.08_spread0.25 \
S3_OUTPUT_URI=s3://BUCKET/output-template-hct \
bash experiments/053-template-emitter-full-hct/run_cloud_after_gamma.sh
```

The script deliberately does not shut down or terminate the worker; lifecycle management remains
with the calling EC2 job.

## Risks and decision boundary

1. **No full HCT116 result yet.** The transformer is tested, but its leaderboard-relevant metric
   effect remains unmeasured until the clean gamma sweep produces a valid HCT116 prediction.
2. **Experiment 049 was a cross-layout proxy.** Its 3/3 DE improvements are encouraging, but its
   H1-prediction/HCT116-truth mismatch prevents treating the magnitude as reliable.
3. **LFC tradeoff.** Template emission worsened the LFC NMAE proxy on all three experiment-049
   targets. An official six-metric comparison is mandatory before adoption.
4. **MSE jackknife can move.** Exact pseudobulk freezes the plug-in profile, not the prediction-side
   sampling correction.
5. **Shared templates.** Reusing one stratified 400-control set for every perturbation may impose
   repeated cell structure. Target-stable seeds change allocation but not the source templates.
6. **Support repair.** Exact row margins require moving counts away from some original template
   support. Per-target retained, removed, and introduced nonzero counts are recorded.
7. **Runtime and disk.** Template allocation is materially slower than quasi, and output sparsity
   can differ. The writer enforces the stored-entry cap but sufficient disk space must be checked
   on the worker.
8. **Harness-only controls.** The transformer preserves the HCT116 scorer file's control block.
   A challenge submission has a different contract and must not be produced directly from this
   harness artifact.

Decision: run the paired cloud scoring after the clean HCT116 gamma sweep, then adopt template
emission only if the reference-scaled gain in Jaccard/fidelity/reach exceeds losses in NMAE and
any MSE jackknife movement. PDS should be exactly unchanged.

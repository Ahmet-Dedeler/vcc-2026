# 114 — DepMap 24Q4 co-essentiality for weak incumbent rows

Status: **source acquisition and deterministic feature construction are
complete; deterministic role-input builders are implemented, but the real role
bundles have not been built or approved. Candidate freeze, held-truth access,
fitting, and evaluation remain locked**.

The official source file and the amended `17,916 × 1,178` float32 feature
matrix were built on 2026-08-30 without opening held perturbation truth.
The feature NPZ SHA-256 is
`7f42b585f788db95cdc345d44b44926bd76ce0fc5ae94c544fd2f6929b9ddc51`;
its logical matrix hash is
`2444fadb3feec310b38823dd962694aea487ea9556e776c7fa545de7ebfe648e`.
All 17,916 rows are finite and unit-normalized, and ten deterministic source
rows reproduced the pinned float32 transform bit-for-bit. The guarded
source-only suite passes 38/38 tests. A separate review is still required
before authorizing source-candidate freeze.

## Question

Can a target-similarity kernel built from DepMap 24Q4 CRISPR fitness profiles
improve only the a14 rows where the incumbent has weak source support?

This is deliberately narrower than experiments 082 and 087. Those experiments
fit all-target ESM2/STRING residual corrections and failed their non-panel
pre-gates. Experiment 114 does not learn a global residual. It replaces only a
frozen source-only subset with a leave-target-out neighbor direction and
retains the incumbent row norm.

Experiment 111 also uses DepMap, but for a different hypothesis: context
program gains from a repeated perturbation-by-context tensor. Experiment 114
uses the older DepMap 24Q4 viability matrix only as a gene-similarity
representation.

## Frozen eligible rows

The exact hash-pinned source target lists report 269/300 K562-GWPS coverage,
290/300 Orion-HCT116 coverage, 294/300 Orion-HEK293T coverage, and no target
lacking all three sources.

`submissions/a15-cd4f010.vcc.meta.json` is not treated as the a14 incumbent.
It is retained only as an audit cross-check of the pre-CD4
`neither_supported_preserved` PerturbMatch support field shared by a14 and
a15. Eligibility is derived directly from the panel, K562, PerturbMatch, and
Orion target lists, all hash-pinned in `SOURCE_MANIFEST.json`.

The candidate may touch only:

- 28 exact a14 `neither_supported_preserved` rows, where neither K562-GWPS nor
  PerturbMatch supports the target;
- 13 additional rows carried by fewer than all three real donors.

Their disjoint lists and ordered hashes are frozen in
`preregistered_protocol.json`; the union contains 41 targets. Every other row
must be bit-identical to the incumbent.

## Exact representation

The implementation follows `pphouse/virtual_cell` commit
`79cd46024dcca8e6afd96bb0bca6bfcbbe77978c`:

1. read DepMap Public 24Q4 `CRISPRGeneEffect.csv`;
2. parse `SYMBOL (ENTREZ)` columns to `SYMBOL`;
3. transpose to genes × 1,178 cell lines;
4. replace every nonfinite value with `NaN`;
5. subtract each gene's finite-value mean;
6. replace missing centered entries with exact float32 zero;
7. L2-normalize each gene row in float32;
8. drop genes with fewer than two finite values; a constant row survives as
   the exact zero feature vector, matching the upstream implementation.

The raw co-essentiality profile is used. Its top-50 neighborhood transform was
inferior upstream and is not an arm here.

The blend is the upstream published weight `1.0`: concatenate row-normalized
STRING and DepMap blocks and normalize again. Where both blocks exist, each
contributes half of the resulting cosine similarity.

The exact DepMap artifact is pinned to Figshare article `27993248`, file
`51064667`, 428,678,699 bytes, MD5
`6edf7ade09b9b34199210b559d4745d3`. It may be downloaded only through the
allowlisted Figshare file ID into `work/source-features/`, with exact byte
count and MD5 verification before atomic publication. The acquisition receipt
also binds the observed SHA-256, protocol, source manifest, and implementation.

## Source-candidate freeze boundary

`SOURCE_CANDIDATE_FREEZE_BOUNDARY.json` is the second-stage trust root. It
authenticates the already-built DepMap artifact and receipt, the STRING source,
the canonical 300-target and 18,533-gene axes, and the exact 41-row eligibility
hashes. It deliberately contains no approved roles and keeps its own
`source_candidate_freeze_authorized` flag false.

A future independently reviewed role entry must pin one directory under
`work/roles/<role>/`, one exact role manifest, and one exact output directory.
It must also declare:

- the held-surface source ID, which must appear in the forbidden-source list;
- the complete allowed source-ID list;
- the exact incumbent endpoint identity, such as an a14-shaped endpoint;
- the exact source-effect space, including rebase, whitening, masking, and
  scaling order.

The role manifest must contain only two artifact groups, `incumbent` and
`source_bank`. Each group contains hash-pinned `artifact` and `receipt`
records. Both receipts must identify a hash-pinned builder and every source
input, assert `source_only=true`, `contains_held_truth=false`, and
`sealed_before_held_truth=true`, and bind every logical array hash.

The incumbent artifact schema is:

```
targets, genes, effects, incumbent_fallback,
source_carrier_count, expected_source_count
```

Its target and gene axes must exactly equal the challenge files. Its fallback
flags and carrier counts must reproduce the hash-pinned K562, PerturbMatch, and
Orion target-list derivation and therefore the exact preregistered 41 rows.

The source-bank artifact schema is:

```
targets, genes, effects, covered, source_ids
```

Its gene axis must be canonical, its source IDs must exactly equal the role
allowlist, the held-surface source must be absent, and all 300 panel response
coordinates must be literal zero and marked uncovered.

The freeze command rejects an arbitrary or self-authenticated input manifest.
The role manifest path and hash must first be added to the boundary by a
separate reviewed amendment. Evaluation remains independently locked after a
successful source-only freeze.

## Deterministic role-input builders

`ROLE_SOURCE_MANIFEST.json` defines the only two authorizing roles:

- `h1`: sources are `k562_gwps`, `orion_hct116`, and `orion_hek293t`; the
  forbidden held source is `h1_2025`.
- `orion_hek293t`: sources are `k562_gwps` and `orion_hct116`; all HEK293T
  effect and control arrays are forbidden. Its `perts.txt` may be read only to
  reproduce the already-preregistered global 41-row eligibility mask.

There is deliberately no `orion_hct116` role because that harness is
depth-confounded and is not an authorizing surface.

The manifest binds the exact a14 settings: alpha 0.982, equal role-donor
weights, `gene_norm=false`, empirical Bayes disabled, beta 0, PerturbMatch
K562-slot blend 0.5, CD4 fraction 0.05 with target solve disabled, whitening
gamma 0.65 and lambda 0, common-essential pruning after the 12,000-row cap,
and NumPy seed 0 consumed once in role-donor order. It also binds the
experiment-107 role operator geometry and all source/implementation hashes.

`build_role_inputs.py` refuses to run a real build while experiment 107's
strict artifact builder is active. Once that build finishes, it authenticates
the corresponding experiment-107 operator without opening its held H5AD,
rebuilds the exact role-specific a14 incumbent on the canonical
`300 × 18,533` axes, and writes:

```
work/roles/<role>/incumbent.npz
work/roles/<role>/incumbent-receipt.json
work/roles/<role>/source-bank.npz
work/roles/<role>/source-bank-receipt.json
work/roles/<role>/manifest.json
```

The source bank is finite by construction. Before reading donor effect values,
the builder computes, from the authenticated STRING and DepMap features, the
union of the exact positive top-25 leave-query-target-out neighbors needed by
all 41 eligible queries in all three non-control arms. Only source rows for
that closed target set are retained. Rows are ordered by the frozen role
source order and then lexicographically by target.

Each retained donor row uses the exact a14 donor decomposition and
role-normalization arithmetic, is masked to the role operator coordinates,
then transformed by the role-specific experiment-107 operator at gamma 0.65.
The bank-global median norm multiplier is omitted because it is one positive
scalar shared by every source row and cancels exactly when candidate rows are
restored to their incumbent norm. PerturbMatch and CD4 remain incumbent-only;
they are not silently relabeled as real donors in the source bank. All
canonical panel and held-surface response coordinates are literal zero and
uncovered.

The NPZ writer fixes member order, timestamps, permissions, dtype, and C-order,
and every artifact is published through a non-overwriting staging directory.
Receipts bind the builder, role-source contract, upstream operator, every
consumed source file, canonical axes, coverage, logical arrays, transforms,
and held-source policy.

## Arms

The four experimental arms are:

1. `current_fallback` — exact incumbent comparator;
2. `string_neighbor`;
3. `depmap_neighbor`;
4. `string_depmap_blend_w1` — preregistered primary.

`identity` is an additional mandatory control and must be array-identical to
`current_fallback`.

For every query, every source signature bearing that same target label is
removed before neighbor selection. Similarities are positive cosine; exactly
the first 25 positive neighbors are retained after descending stable sort, so
ties resolve by the frozen sorted source-target order; weights are then raised
to power 3. Source rows are averaged per response coordinate only where that
coordinate was measured. All 300 VCC target-gene coordinates are removed
before aggregation and remain exactly zero. A valid candidate is scaled back
to the incumbent row's allowed-coordinate L2 norm. Confidence is the absolute
value of the final stored effect array after dtype conversion.

Missing features, no positive neighbors, zero/nonfinite candidate rows, and
zero-norm incumbent rows all return the incumbent row bit-for-bit.

## Evaluation contract

No evaluation is currently authorized. Before any held perturbation truth may
be opened, the runner must freeze:

- all real arm matrices;
- source-blind confidence arrays;
- all 100 target-label permutation-null candidates;
- exact target/gene/eligibility axes;
- input, implementation, and output hashes.

PDS is raw and uncentered, uses the complete predeclared gallery, midranks
ties, and excludes all 300 VCC target response coordinates.

Directional precision selects the top 20% of response genes by absolute
frozen predicted effect before truth access. A selected gene is correct only
when the held reference marks it directional and its sign matches. Selection
is therefore source-blind.

The 100 nulls permute query-target feature labels within the two frozen
eligibility strata and repeat primary-arm construction.

Promotion requires all of:

- nonnegative raw-PDS delta on every clean held surface;
- macro raw-PDS delta at least `+0.003`;
- nonnegative directional-precision delta on every clean directional surface
  and strictly positive macro improvement;
- observed macro PDS gain above the target-label permutation p95;
- exact identity/current-fallback equality;
- exact incumbent fallback for all ineligible or unsupported rows;
- unchanged target coverage and every leakage/provenance invariant.

## Current commands

Read-only contract audit:

```bash
uv run python experiments/114-depmap-coessentiality/run.py manifest
```

The only authorized large download:

```bash
uv run python experiments/114-depmap-coessentiality/run.py \
  download-depmap --confirm-large-download
```

Build the deterministic source feature artifact after the verified download:

```bash
uv run python experiments/114-depmap-coessentiality/run.py build-depmap
```

Both commands use fixed paths under
`experiments/114-depmap-coessentiality/work/source-features/`, refuse existing
outputs and symlinks, and never open held perturbation outcomes.

Synthetic tests:

```bash
uv run python -m pytest -q \
  experiments/114-depmap-coessentiality/test_run.py \
  experiments/114-depmap-coessentiality/test_build_role_inputs.py
```

Lightweight role readiness check:

```bash
uv run python experiments/114-depmap-coessentiality/build_role_inputs.py check
```

After experiment 107 has finished successfully, build the real role bundles
strictly serially:

```bash
uv run python experiments/114-depmap-coessentiality/build_role_inputs.py \
  build-role --role h1 --confirm-source-only-real-build

uv run python experiments/114-depmap-coessentiality/build_role_inputs.py \
  build-role --role orion_hek293t --confirm-source-only-real-build
```

Do not run the second command until the first exits. These commands build only
source artifacts; they do not freeze candidates or access held outcomes.

`freeze` and `score` remain fail-closed. `freeze` additionally requires a
boundary-approved role manifest, but the approved-role map is currently empty.
A separately reviewed amendment must pin each real role bundle and set both
authorization flags before any source candidate is frozen. A second review is
required before held-truth scoring.

No count generation, `cell-eval2`, challenge build, submission, commit, push,
or write outside this directory is part of experiment 114.

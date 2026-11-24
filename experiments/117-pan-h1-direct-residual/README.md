# 117 — Pan2025 H1 direct residual

Status: **one-shot held evaluation complete — NO-GO.** The immutable result is
in `work/held-evaluation-v1/`; no counts or submission were produced.

The first authorized build exposed one source-only serialization edge case:
`memoryview.cast("B")` rejects legitimate arrays with a zero-length dimension.
Logical-array hashing now preserves the existing `exp107-array-v1` prefix,
dtype, and shape header and contributes no payload bytes when `nbytes == 0`.
No role was published by the failed attempt, and the output path remains absent.

Date: **2026-08-30**

## Question

Can the public GSE295214/Pan2025 H1 CRISPRi effects improve the exact a14
direction as a narrow same-target residual, especially for the 57 VCC-2026
panel rows absent from the older VCC-2025 H1 donor?

The source contains 64 of 300 VCC targets after author QC. Seven already occur
in the VCC-2025 H1 donor and 57 are new:

- all-supported arm: all 64 Pan-covered panel rows;
- primary novel-only arm: the 57 rows absent from VCC-2025 H1;
- every other target and every unsupported response coordinate: exact
  incumbent fallback.

## Frozen scientific design

The candidate uses Experiment 107's source-only exact-a14 role bundles. For
each held role, the relevant incumbent and its serialized leave-surface-out
whitener are authenticated before arrays are loaded. Pan is never inserted
into the whitening basis.

For a Pan-supported target:

1. align the compact Pan effect and measured mask to the canonical 18,533-gene
   axis;
2. remove the union of all 300 VCC panel coordinates and the role's complete
   held-target axis;
3. intersect with the role operator columns and the row-specific Pan measured
   mask;
4. zero-pad the supported raw source row into the operator axis, apply the
   exact serialized a14 whitener (`gamma=0.65`, `lambda=0`), then mask again;
5. match the transformed source norm to the incumbent norm on that exact
   supported coordinate set;
6. form and persist the actual normalized residual
   `scaled_pan - destination_incumbent`;
7. mark a row effective only if its support, source norm, incumbent norm, and
   realized residual norm are finite and nonzero;
8. blend only effective supported coordinates:

   `candidate = incumbent + f * (scaled_pan - incumbent)`.

The primary fraction is fixed at `f=0.10`. `f=0.05` and `f=0.20` are frozen
dose diagnostics, not alternative promotion candidates. The identity arm is
bit-exact `f=0`.

No missing target, unmeasured source coordinate, operator-missing coordinate,
excluded panel coordinate, held-gallery coordinate, nonfinite row, zero-norm
row, or zero-residual row may be imputed. It remains bit-identical to the
incumbent. The frozen arrays persist separate all-supported and novel-only
effective-row masks.

Directional ranking uses absolute **actual normalized residual**, not the
transformed Pan source vector. It is deterministic: descending magnitude with
canonical gene index as the tie-break, and fallback rows have only `-1`
rank entries.

## Validation roles

- `h1`: same-context replication diagnostic only. It cannot authorize
  promotion. Its 30 Pan-overlap rows are all already present in VCC-2025 H1,
  so the primary novel-only arm is intentionally exact identity here.
- `orion_hek293t`: mandatory clean cross-context gate.
- `nadig.gallery_1` through `nadig.gallery_6`: one mandatory clean Jurkat
  cross-context surface, aggregated with equal gallery weight. This is a
  reused development benchmark, not an independent confirmation, and cannot
  authorize promotion without HEK293T.

All roles use the exact a14 held-role geometry frozen by Experiment 107:
leave-surface-out donors, serialized operators, `alpha=0.982`, equal donor
weights, PerturbMatch fraction `0.5`, CD4 fraction `0.05`, and exclusion of
the VCC panel plus the complete held target gallery.

Every role's exact target-symbol overlap and novel-target overlap are enforced
against `SOURCE_INPUT_MANIFEST.json`; global 64/57 panel coverage is not
accepted as a substitute for a role-level check.

Promotion aggregation is also frozen:

- each surface raw-PDS delta is the complete-gallery arithmetic mean;
- the six Nadig galleries receive equal weight in one Jurkat macro;
- HEK293T and that Jurkat macro then receive equal 1/2 context weight;
- H1 is excluded from every promotion statistic.

For each target-label null, the frozen destination-to-source map moves the
source direction and its support together. Norm matching is recomputed against
the destination incumbent; observed residuals are never shuffled directly.
Only effective rows are eligible for permutation.

The 20,000-draw whole-cluster bootstrap derives each role/draw PCG64 seed from
the first eight little-endian bytes of SHA-256 over the frozen namespace,
role, and zero-based draw separated by NUL bytes. It samples the role's sorted
unique clusters with replacement, reconstructs complete-gallery surface
means, applies the same equal-context aggregation, and uses NumPy's `linear`
(Hyndman-Fan type 7) 0.025 quantile as the lower bound.

## Promotion boundary

The primary comparison is `pan_novel_f100` versus exact identity. Promotion
requires every condition in `preregistered_protocol.json`, including:

- nonnegative raw-PDS and directional-credit deltas on HEK293T;
- nonnegative raw-PDS delta in every Nadig gallery and positive Jurkat macro
  directional credit;
- macro cross-context raw-PDS gain of at least `+0.003`;
- gain above all frozen target-label permutation nulls;
- a positive target-cluster bootstrap lower bound;
- exact fallback and exclusion invariants.

H1 can provide a mechanism/reproducibility diagnostic for `pan_all_f100`, but
no H1 result can replace either mandatory cross-context gate.

## Current authorization

`SOURCE_CANDIDATE_FREEZE_BOUNDARY.json` authorizes only the source-candidate
build. It pins:

- the compact receipt SHA-256;
- the artifact-manifest SHA-256;
- the compact `delta.npy` and `measured_mask.npy` bindings;
- the source manifest, acquisition amendment, download receipt, extractor, and
  original public pseudobulk hash/size chain;
- this builder's final SHA-256;
- explicit approval for every role.

The builder authenticates those bindings before reading the compact source.
It refuses to overwrite an existing output directory. The authorization does
not permit held-truth access, evaluation, `cell-eval2`, count generation,
challenge packaging, or submission.

```bash
uv run --with pytest pytest \
  experiments/117-pan-h1-direct-residual/test_build_source_candidates.py
uv run ruff check experiments/117-pan-h1-direct-residual
```

Safe preflight:

```bash
uv run python experiments/115-gse295214-h1/extract_compact_donor.py \
  --verify-only \
  --output-dir experiments/115-gse295214-h1/work/compact-donor
uv run python \
  experiments/117-pan-h1-direct-residual/build_source_candidates.py \
  verify-scaffold
```

Authorized source-only build command, to be run only after the final
independent review:

```bash
uv run python \
  experiments/117-pan-h1-direct-residual/build_source_candidates.py \
  build \
  --output-root \
  experiments/117-pan-h1-direct-residual/work/source-candidates
```

This experiment never runs `cell-eval2`, emits cells/counts, builds a challenge
package, or submits.

## Held result

The repaired revision-2 evaluator passed independent review before opening any
held value. It then evaluated only the predeclared `pan_novel_f100` arm against
exact identity on HEK293T and the six Nadig Jurkat galleries. H1 used only the
predeclared `pan_all_f100` diagnostic arm and could not authorize promotion.

The candidate was consistently positive but too small:

- cross-context raw-PDS delta: `+0.0011948386`;
- HEK293T raw-PDS delta: `+0.0000578302`;
- Jurkat macro raw-PDS delta: `+0.0023318469`;
- target-cluster bootstrap lower 95%: `+0.0008382064`;
- largest frozen target-label null: `+0.0004315281`;
- HEK293T directional-precision delta: exactly `0`;
- H1 same-context diagnostic raw-PDS delta: `+0.0194677871`.

The frozen `+0.003` cross-context minimum and strictly positive HEK293T
directional-credit gate both failed. Do not tune the fraction or reopen this
held outcome. The exact Pan direct-residual candidate is closed.

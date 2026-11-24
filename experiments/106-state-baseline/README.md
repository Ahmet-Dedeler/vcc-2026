# 106 — official STATE VCC baseline reproduction

Status: **NO-GO for the reproduced `state_sm` step-20,000 checkpoint**

Date: **2026-08-30**

## Question

Can the official Arc Institute STATE VCC recipe provide an independent
prediction channel that improves the current a14 submission when validated
outside the challenge leaderboard?

This is not a claim that a model name on the public board identifies an exact
method. The board observation only motivates reproducing Arc's public recipe.

## Frozen upstream sources

- STATE repository:
  `https://github.com/ArcInstitute/state.git`
- Pinned `main` commit:
  `9bbfe78a434a55205e4de834e1ea99f85f7a3add`
- Official VCC notebook:
  Google Drive id `1QKOtYP7bMpdgDJEipDxaJqOchv7oQ-_l`
- Downloaded notebook SHA-256:
  `0b3888b9a36e6fbfc056b9e5585d825aa5a97a92f34d3bc2ff69cba064f21422`
- Official support bundle:
  `https://storage.googleapis.com/vcc_data_prod/datasets/state/competition_support_set.zip`
- Expected compressed bytes:
  `8716992349`
- Recorded HTTP ETag:
  `"-CODBsJzT95ADEAI="`
- Verified archive SHA-256:
  `f1d5fd56a2eee240e5585fac52bb75b0b830a7be6a5f6de1bf773bff4f507f7c`
- Verified GCS CRC32C:
  `MBP7uA==`
- Verified uncompressed bytes:
  `48,691,267,232`

## First-stage reproduction

Follow the official notebook exactly before changing architecture:

```text
data.kwargs.toml_config_path=competition_support_set/starter.toml
data.kwargs.num_workers=8
data.kwargs.batch_col=batch_var
data.kwargs.pert_col=target_gene
data.kwargs.cell_type_key=cell_type
data.kwargs.control_pert=non-targeting
data.kwargs.perturbation_features_file=competition_support_set/ESM2_pert_features.pt
training.max_steps=40000
training.ckpt_every_n_steps=20000
model=state
```

The notebook performs challenge inference from the step-20,000 checkpoint.
Both step 20,000 and step 40,000 may be retained, but the first reproduction
must not select between them using challenge leaderboard outcomes.

## Reproduction audit

Verified on **2026-08-30**:

- The local clone is detached at the pinned commit
  `9bbfe78a434a55205e4de834e1ea99f85f7a3add`.
- The notebook on disk matches the frozen SHA-256 above.
- The isolated environment uses CPython 3.11.15. Its generated `uv.lock` has
  SHA-256 `ce7ee30ac4dbcfbb879500747ce04d0e310c87fafbff6c09eff36d2064d588a1`.
- The resolved environment currently contains PyTorch 2.13.0 and Transformers
  5.16.1. STATE itself did not ship a lockfile at the pinned commit, so this
  dependency resolution is part of the experiment provenance.
- All 37 upstream tests pass in the isolated environment.
- PyTorch reports that MPS is built and available on this Mac, and a basic
  matrix operation succeeds on MPS.

The notebook is not internally immutable. Its current source says
`model=state`, but the saved successful execution output is a 162-million
parameter model with `cell_set_len=128`, `hidden_dim=672`, and four Llama
layers. That output matches the repository's `model=state_sm` profile, while
the pinned commit's current `model=state` profile uses `cell_set_len=512`,
`hidden_dim=768`, and eight Llama layers. The notebook output also records a
6.69 GB support download, whereas the current server object is expected to be
8,716,992,349 bytes. Therefore these are treated as two distinct candidates:

1. **current-source arm:** the notebook command exactly as written at the
   pinned commit (`model=state`);
2. **archived-output arm:** `model=state_sm`, matching the architecture that
   actually produced the notebook's recorded training output.

The first bounded timing smoke test will compare these profiles without
opening challenge outcomes. Full training will use only a profile that fits
memory and has a defensible provenance label.

Upstream training selects CUDA when available and otherwise forces CPU, even
when Apple MPS is available. `mps-runtime.patch` adds an opt-in
`STATE_USE_MPS=1` device path for training and inference. CUDA and default CPU
behavior remain unchanged. This is a runtime-placement patch only; it does not
alter data, architecture, loss, optimizer, or predictions by design.

The live leaderboard snapshot used to prioritize this work contained 447
teams. `GZ YAI` was rank 2 at `0.1860290893` with the public model name
`STATE baseline`, versus a14 at rank 10 and `0.1704242443`. STATE was lower
than a14 on PDS but higher on MSE, NMAE, fidelity, reach, and Jaccard. That
metric complementarity is why the validation plan includes a frozen
STATE/a14 blend rather than treating the official baseline only as a
replacement candidate.

Archive verification and extraction are implemented by
`verify_and_extract.py`. It checks the exact HTTP content length, the GCS
CRC32C value, computes SHA-256, rejects path traversal and symbolic-link
members, checks 20 percent disk headroom against the archive's declared
uncompressed size, extracts through a temporary directory, and records a
member-level manifest. Its safety checks have a separate synthetic test.

The current archive has two packaging inconsistencies relative to the saved
notebook execution:

1. it is wrapped in a top-level `competition_support_set/` directory, so the
   extractor unwraps that single matching directory instead of creating a
   duplicate nested path;
2. all six training matrices contain log1p-scale values, but only
   `competition_train.h5` carries the standard AnnData `uns/log1p` marker.
   `cell-load` 0.10.4 rejects the mixed metadata before training.

The verified extraction remains untouched. `prepare_support_view.py` creates
an APFS copy-on-write view and adds only the missing empty `uns/log1p`
metadata mappings to the five donor files. `prepare_local_config.py` verifies
the frozen `starter.toml` hash and rewrites only its Colab data root to that
working view. The expression matrices, split declarations, and source archive
are not modified.

`run_smoke.sh` performs a guarded one- or two-step MPS training smoke test
against the verified support bundle. It disables implicit CPU fallback so an
unsupported MPS operation fails visibly, disables external W&B logging, uses
one process with a 40 GiB tree limit, and accepts only the two predeclared
profiles `state` and `state_sm`.

`run_training.sh` is the resumable full runner. Its primary invocation is
`bash experiments/106-state-baseline/run_training.sh state 20000`; the current
training code writes `last.ckpt` and `final.ckpt`, so step 20,000 is run as an
explicit first stage rather than assuming the notebook's stale
`step=20000.ckpt` filename still exists. A later invocation with 40,000 steps
resumes from `last.ckpt`.

`benchmark_architecture.py` is the data-independent preflight. It instantiates
the selected profile at the notebook's exact 18,080 expression dimensions and
5,120 perturbation-feature dimensions, then runs a real forward pass, energy
loss, backward pass, and optimizer step on MPS.

Both profiles passed that preflight with CPU fallback disabled:

| profile | parameters | cell set | Metal allocated | one warm optimization step |
| --- | ---: | ---: | ---: | ---: |
| `state_sm` | 162,626,484 | 128 | 3.05 GiB driver | 0.146 s |
| `state` | 213,551,124 | 512 | 4.18 GiB driver | 0.232 s |

These are synthetic steady-state timings, not full-loader throughput
estimates. They establish that both exact architectures fit comfortably and
that their forward, energy-loss, backward, and optimizer paths execute
natively on MPS without fallback.

Both profiles also passed real-data training:

| profile | batch size | steps | peak process-tree RSS | observed rate |
| --- | ---: | ---: | ---: | ---: |
| `state` | 1 | 2 | 3.63 GiB | 1.08 step/s |
| `state_sm` | 1 | 2 | 1.36 GiB | 2.24 step/s |
| `state` | 16 | 2 | 5.72 GiB | 0.12 step/s |
| `state_sm` | 16 | 10 | 4.61 GiB | 0.50 step/s including startup |

The sustained 20,000-step `state_sm` run is faster after warm-up
(approximately 2.1–2.3 steps/s in its CSV telemetry) and is the first full
arm because it matches the notebook's saved successful architecture. The
larger current-source `state` profile remains a separate later arm; its
two-step timing projects to roughly 46 hours on this machine.

The support TOML has no explicit validation split. `cell-load` substitutes the
held-out HepG2 test loader as validation, so `best.ckpt` is not an admissible
selection artifact. Use the preregistered fixed 20,000-step `last.ckpt`,
record its hash and global step, and call the HepG2 result transductive because
HepG2 controls remain in training.

The current checkpoint callback ignores `training.ckpt_every_n_steps` and
saves on `training.val_freq`; `infer` defaults to `final.ckpt` while `predict`
defaults to `last.ckpt`. Every evaluation and inference command must therefore
pass an explicit immutable checkpoint path.

Upstream inference is not a submission writer. It expects the 18,080-gene
STATE order, densifies its input, and emits floating-point log-expression.
The 2026 challenge uses 18,533 genes and raw sparse integer counts. Challenge
inference therefore requires a separate name-aligned, streaming adapter:
18,533 controls -> the 18,080 STATE input order -> predicted STATE residual ->
18,533 challenge order, with zero residual on the 456 challenge-only genes,
followed by the existing count emitter. Never use `state tx infer --all-perts`;
the saved ESM map has 19,790 keys, not the 300-target challenge panel.

## Held-context implementation

The fixed-step evaluation path is implemented locally:

- `freeze_checkpoint.py` refuses to read a live training checkpoint, requires
  global step 20,000 and the exact `state_sm` model contract, then copies and
  hashes the checkpoint, config, hparams, axes, dependency lock, saved
  perturbation map, and a clean snapshot of the patched STATE source into an
  immutable directory. Inference imports that frozen source snapshot rather
  than the mutable checkout.
- `prepare_held_controls.py` builds one deterministic 18,400-cell control-only
  reservoir per Orion context. Selection observes only rows labelled as
  controls and never uses perturbation expression or perturbation depth.
- `freeze_a14_endpoints.py` exports the exact held-context a14 analogue
  (`alpha=0.982`, `CD4_FRACTION=0.05`) without reading held perturbation
  expression, using the same independently frozen control-only input as STATE.
- `state_effect_adapter.py` accepts control rows only. It aligns genes by exact
  symbol, applies `log1p` to raw counts without CP10K normalization, uses one
  frozen bank of four 128-cell sets per context for every target, rejects
  missing/zero perturbation embeddings, and runs the model directly instead of
  the dense upstream `tx infer` command.
- STATE outputs are interpreted as final log-expression, converted with
  `expm1`, and aggregated before 50k log-pseudobulk construction. The shared
  gene block is mass-matched to the control block; challenge-only genes retain
  control expression and therefore exact zero STATE effect.
- `score_state_blends.py` opens held truth only after inference is sealed. It
  compares the exact incumbent, STATE with incumbent fallback on unsupported
  genes, and fixed composition/effect/direction blend families at weights
  `0.10,0.20,0.35,0.50,0.65,0.80`. Every candidate is converted through one
  final 50k-normalized profile with no post-normalization coordinate overwrite,
  and PDS uses `cell-eval2`'s cosine-distance convention.
- `run_held_screen.sh` executes checkpoint freeze, exact a14 export, STATE
  inference, and scoring sequentially under memory guards. It refuses to start
  while the STATE trainer is active. Both contexts are scored on their exact
  common target intersection. The sorted HCT116 and HEK293T axes are pinned at
  132 targets each, and their exact 120-target intersection is revalidated
  during both allowlist construction and scoring.

The current synthetic/local contract suite contains 51 tests. It covers the
18,077/456/3 axis split, exact checkpoint shape/step requirements, deterministic
control-only reservoir, control-bank reuse, target-homogeneous 128-cell
sequences, streamed count aggregation, frozen runtime and incumbent provenance,
unsupported-gene fallback, panel masking, 50k mass closure, official zero-norm
PDS behavior, blend invariants, and truth-free inference signatures. It also
regresses the PyArrow column-projection API, bytecode-excluded source snapshot
hashing, exact 18,400-control count, held-operator path/hash binding, and common
gallery count/order hashes.

## Result

The fixed 20,000-step run completed and was frozen successfully. On the
120-target HCT116/HEK293T common gallery, pure STATE was effectively random by
PDS:

| candidate | HCT116 PDS | HEK293T PDS | macro delta vs a14 |
| --- | ---: | ---: | ---: |
| a14-shaped incumbent | 0.80917 | 0.79454 | 0 |
| STATE fallback | 0.50252 | 0.50028 | **-0.30046** |

Every preregistered global composition, effect, and direction blend lost PDS.
Even the smallest 10% direction blend lost `-0.01513` on HCT116 and `-0.03634`
on HEK293T; HEK293T directional credit also fell `-0.00603`.

Experiment 113 then tested the same immutable checkpoint on the in-domain H1
held-target panel. Pure STATE again scored at chance (`PDS 0.50448`) versus
the frozen incumbent's `0.86345`. This rules out treating the Orion result as
only an out-of-context transfer failure.

The current `state_sm` checkpoint is therefore not a production candidate.
The larger current-source `state` architecture remains a scientifically
distinct untested arm, but the public leaderboard model string `STATE
baseline` is not evidence that this exact checkpoint or architecture produced
that score.

Full artifacts and metrics are in
`work/runs/20260830T064800Z/a14/summary.json`; the H1 diagnostic is in
`../113-state-h1-domain-check/`.

## Decision contract

1. Verify the support archive size and SHA-256 before extraction.
2. Inspect the archive schema and `starter.toml`; do not silently rewrite
   dataset partitions.
3. Install STATE at the pinned commit in an isolated environment.
4. Run a small CPU/MPS compatibility smoke test before committing to a full
   training run.
5. The official challenge prediction is an artifact, not validation truth.
   Do not tune an ensemble weight on the challenge board.
6. The Orion screen is a panel-target diagnostic on an abundance-selected
   partial gallery. It may falsify a blend family, but it cannot select or
   authorize a production blend weight.
7. Exclude all 300 VCC target-gene coordinates from PDS geometry.
8. Use one guarded heavy process at a time. No parallel `cell-eval2` runs.
9. No challenge submission is authorized by acquisition or training alone.

The reported `plugin_nmse_proxy` is an effect-level squared-error diagnostic.
It is not the official count-level, own-target-masked, jackknife-corrected
`expr_mse_unbiased_capped_norm`.

## Scope

- Writes for the reproduction belong under `experiments/106-state-baseline/`
  and `data/external/state-vcc/`.
- No current production builder, submission artifact, or prior experiment is
  modified during acquisition and smoke testing.

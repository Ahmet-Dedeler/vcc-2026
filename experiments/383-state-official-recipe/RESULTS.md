# 383 — STATE official recipe: pricing the working-checkpoint ensemble leg

**Status: IN PROGRESS** (download running; known-answer passed; training pending)

## 1. Recipe + provenance (deliverable 1, exact official 2026 recipe)

The 2026 challenge-site baseline path and the `ArcInstitute/state` VCC path converge
on ONE recipe (the on-disk notebook is the 2025 H1 few-shot version of it; the train
command is identical — only the inference template differs: 2025 `competition_val_template.h5ad`
vs 2026 the three new contexts' control cells + 300-gene target list):

- Challenge site: https://virtualcellchallenge.org/ (2026: zero-shot prediction across new
  cellular contexts; participants receive non-targeting control profiles + 300 gene IDs per
  line; scored by cell-eval `--preset vcc2026`; accessed 2026-09-11). The Cell 2026 paper
  (S0092-8674(26)00931-1) confirms the design: six CRISPRi lines, three released as
  validation inputs. No separate 2026 STATE notebook exists publicly; the recipe below is
  the current official one.
- Model docs: https://github.com/ArcInstitute/state, README "Train STATE for the Virtual
  Cell Challenge" → Colab https://drive.google.com/drive/folders/1QKOtYP7bMpdgDJEipDxaJqOchv7oQ-_l
  (same Drive ID as 106 froze on 2026-08-30; README unchanged at pinned commit; accessed
  2026-09-11). Evaluation framework https://github.com/ArcInstitute/cell-eval.
- Repo pin (read-only, local): `data/external/state-vcc/state` @
  `9bbfe78a434a55205e4de834e1ea99f85f7a3add` (2026-07-23). Working tree carries 106's
  uncommitted MPS runtime patch (`mps-runtime.patch`: accelerator `mps` in `_train.py`,
  `map_location=cpu`+`.to(device)` in `_infer.py`/`_predict.py`, gated on `STATE_USE_MPS=1`).
  Venv `data/external/state-vcc/state/.venv`: CPython 3.11, torch 2.13.0, MPS built+available
  (verified this session).
- Notebook on disk `data/external/state-vcc/official_vcc_notebook.ipynb`
  (sha256 `0b3888b9...92f34d3bc2ff69cba064f21422` per 106 PROVENANCE.json), cells 9/10/14/17:
  - Data: `https://storage.googleapis.com/vcc_data_prod/datasets/state/competition_support_set.zip`
    (8.7 GB; uncompressed 48.7 GB; 12 zip members: 6 training `.h5`
    `competition_train/k562_gwps/rpe1/jurkat/k562/hepg2` + `competition_val_template.h5ad`,
    `competition_test_template.h5ad`, `starter.toml`, `gene_names.csv`, `ESM2_pert_features.pt`).
    Expected sha256 `f1d5fd56a2eee240e5585fac52bb75b0b830a7be6a5f6de1bf773bff4f507f7c`
    (106 `bundle-manifest.json`). Support bundle DELETED from this machine after 106
    (only notebook + repo remained); re-downloading into `experiments/383-state-official-recipe/data/`
    (PLAN.md notes the >2 GB download).
  - Train (cell 14, exact): `state tx train
    data.kwargs.toml_config_path=competition_support_set/starter.toml data.kwargs.num_workers=8
    data.kwargs.batch_col=batch_var data.kwargs.pert_col=target_gene
    data.kwargs.cell_type_key=cell_type data.kwargs.control_pert=non-targeting
    data.kwargs.perturbation_features_file=competition_support_set/ESM2_pert_features.pt
    training.max_steps=40000 training.ckpt_every_n_steps=20000 model=state
    output_dir=competition name=first_run` (wandb tags only).
    NOTE: notebook source says `model=state`, but its saved output (162 M params,
    cell_set_len=128, hidden 672, 4 llama layers) matches the repo's `model=state_sm`
    profile; the pinned commit's `model=state` = cell_set_len 512, hidden 768
    (106 README + smoke-state/config.yaml re-verified this session). 106/113 trained
    `state_sm` to step 20,000 → held-line pds 0.50 (no identity). The full `state` profile
    was only smoke-configured, never trained. THIS experiment trains `model=state`.
  - Inference (cell 17, exact): `state tx infer --model-dir competition/first_run
    --checkpoint competition/first_run/checkpoints/step=20000.ckpt
    --adata <context control cells> --pert-col target_gene`, then
    `cell-eval prep -i prediction.h5ad -g gene_names.csv`.
  - Starter split (106 `runs/starter.local.toml`): all `replogle_h1` cell types train
    except `"replogle_h1.hepg2" = "test"` (whole-line zeroshot) — the identity-check split.
- Board context (snapshot `data/leaderboard/2026-09-11T1910Z.json`, 893 entries, read this
  session): best STATE-flavoured entry `(STATE+c2431c)` (GZ YAI) 0.1912, pds_cosine 0.8172,
  nmae scaled 0.1394 — matches the assignment's 0.17–0.19 / 0.82 / 0.137 (it is a blend,
  not the pure recipe). Pure-recipe-shaped attempts cluster at pds ≈ 0.50 with negative
  overalls: `state` (cl9) pds 0.5161/avg −0.1874, `Sylvalso-STATE40K` pds 0.5091/avg −0.0626,
  `VCC26_STATE_B0` pds 0.4959, `STATE-5ctx-HVG-delta-40k` pds 0.4984 — the public
  corroboration of 106/113/182's no-identity finding. A working identity-bearing STATE
  checkpoint has therefore never been priced; that is what this experiment trains.

## 2. Known-answer reproduction

`known_answer383.py`: 323's `load_reference()` (hash-checked) + `score()` on saved
`work/stats-cap15.npz`, outputs to `work/known-answer.json`. Result this session:

| member | recomputed | recorded (`metrics-cap15.json` raw) | |diff| |
|---|---|---|---|
| pds | 0.8896358543417369 | 0.8896358543417369 | 0.0 |
| mse | 0.9749462418923823 | 0.9749462418923823 | 0.0 |
| nmae | 0.9508364630325662 | 0.9508364630325662 | 0.0 |
| fid | 0.4975545182292535 | 0.4975545182292535 | 0.0 |
| reach | 0.25434726355678233 | 0.25434726355678233 | 0.0 |
| jac | 0.07475983298131939 | 0.07475983298131939 | 0.0 |

PASS (all ≤ 1e-12; exact 0.0 — same code path on saved stats, as expected).
Caveat: 323's `E.verify()` manifest gate no longer passes because the builder evolved
(committed `c47d2343` VCC_CTX_A_PRIOR chain); working tree clean, so this is committed
history, not corruption. `load_reference()` hash checks on the statistics that enter
`score()` pass. Environment: repo `.venv` python3.12, cell_eval2 as installed in the repo venv
(`evaluation-inputs.json` scorer-version check skipped with `verify()`; the DE-metric
versions are those of the current checkout).

## 3. Throughput probe (feasibility measurement)

`train383.sh 40 probe-state` (exact recipe config, `model=state`, batch default 16,
num_workers 8, MPS via `STATE_USE_MPS=1`): 40 steps in 223 s wall → **0.18 steps/s**
(lightning reports 0.24 it/s steady-state; wall includes ~40 s startup). Peak RSS
23.59 GiB (cap 40; system free dipped to 27–43% during the run). Checkpoint 2.5 GB.
Projection: recipe inference checkpoint step-20000 ≈ 31 h; full 40k ≈ 62 h. 5,000 steps
≈ 7.7 h — **slower than the 2 h / 5,000-step stop-rule bar by ~4×**. Full recipe
infeasible on this Mac within the night. Largest feasible run in a ~4.5 h window:
≈ 2,900 steps = **14.5% of the way to the recipe's step-20000 inference checkpoint**
(7.25% of the 40k-step recipe). For reference, 106's `state_sm` ran at 1.25 it/s;
the full `state` profile (cell_set_len 512, hidden 768) is ~7× slower per step.
Launched: `CKPT_EVERY=1000 train383.sh 2900 official-state`. CORRECTION logged 22:15Z:
at the pinned commit `get_checkpoint_callbacks()` (`src/state/tx/utils/__init__.py:130`)
ignores its `_ckpt_every_n_steps` argument — the notebook's `training.ckpt_every_n_steps`
is a DEAD knob in this code version; checkpoints (`best.ckpt` on val_loss, `last.ckpt`)
save only on validation ticks, `every_n_train_steps=val_freq=2000`. So intermediate
trajectory = `last.ckpt`@step-2000 + `final.ckpt`@step-2900 only (no dynamics impact;
save-only difference from the recipe). First checkpoint due ≈ 22:55Z. Identity check +
H1 scoring to follow on the final checkpoint.

## 4. Identity-check design + truth provenance (gate: pds >= 0.7; stop <= 0.55)

Recipe-native held set: `competition_val_template.h5ad` = 98,927 H1 cells (dense float32
log1p, max 7.2): 50 named perturbations (161–2,925 cells each) + 38,176 non-targeting
controls. Truth `work/truth_val.npz`: chunked per-pert mean(X) − mean(controls), 50×18,080,
delta norms 0.96–4.0 (median 1.69), mean|delta| 0.0094 — subtle H1 effects, as expected.
Inference template `work/templates/identity_controls.h5ad`: 1,000 val-template controls
(seed 383) + `identity.tsv` (50 perts × 100 cells via infer's TSV-pad path, which clones
controls and relabels). Metric `identity383.py`: cosine-rank retrieval of each TRUE delta
vs all 50 PREDICTED deltas (pds direction, midrank ties). Training-space note: all support
`.h5` matrices are log1p (train nonzero-data max 6.7/mean 1.5; hepg2 dense max 7.1), model
config `log1p_from_raw_counts=false`, so H1/A/B/C templates from raw counts are log1p'd at
build (`build_templates383.py`); the val-template identity leg needs no transform by me,
which makes it the clean pipeline control. Pert-map check on the probe run:
19,792 ESM2 5120-dim entries covering **300/300** panel targets — the recipe has an
inference path to every challenge target (unlike 182's one-hot checkpoints). Gene spaces:
support 18,080 = h1-test var order (asserted at build); challenge 18,533 overlaps by
18,077 (456 challenge genes get zero delta at A/B/C conversion).
...[truncated 615 chars]
## 5. CLOSURE by the PI, 2026-09-13 06:00Z — training killed by the 2026-09-12 00:22Z reboot; step-2000 checkpoint has NO identity

The `official-state` run (model=state, 2,900-step target) died in the machine reboot at 2026-09-12 00:22Z with
`last.ckpt`/`best.ckpt` saved at the step-2000 validation tick (23:09Z). Inference on that checkpoint
(`work/infer-identity-step2000.log`, MPS, TSV-pad path, `ANNDATA_ALLOW_WRITE_NULLABLE_STRINGS=1` needed for the
output write) and `identity383.py` give:

| quantity | value |
|---|---:|
| pds_identity (cosine-rank retrieval, 50 H1 val perts) | **0.5049** |
| mean own-target cosine | −0.0089 |
| verdict (gate ≥ 0.7 PASS / ≤ 0.55 STOP) | **STOP-NO-IDENTITY** |

Record: `work/identity_step2000.json`. This is 5 % of the recipe's 40k steps, so it does not price the full recipe;
but the public board does: every pure-STATE-shaped entry on 2026-09-13 sits at pds ≈ 0.50 on A/B/C (`state` cl9
0.5161, Sylvalso-STATE40K 0.5091, VCC26_STATE_B0 0.4959, STATE-5ctx-HVG-delta-40k 0.4984), i.e. the official recipe
has no zero-shot identity on a NEW line for anyone, and the PLAN's blend hypothesis was an in-domain (H1) check that
cannot transfer to A/B/C. **No cloud training is bought; the STATE family is closed as an independent source.**
The 98 GB support set and 10 GB of checkpoints under `data/` and `runs/` are gitignored and can be deleted.

## Storage (2026-10-01T13:14Z, coordinator, with Ahmet's OK conditional on an independent review)
Deleted by exact path (`git clean -xdf -- …`): `data/competition_support_set/`, `data/competition_support_set-log1p/`,
`runs/` (all gitignored, untracked). An independent read-only review (Opus agent) agreed: the kept zip
`data/competition_support_set.zip` passes `unzip -tq`, its sha256 matches 106's `bundle-manifest.json`, and all 11 files'
CRC32 match the zip, so the support set is rebuildable byte for byte; the log1p copy is `prepare_support_view.py`'s clone
plus an empty `uns/log1p` group; the checkpoints' one result (identity 0.5049) is in tracked `work/identity_step2000.json`.
The two support folders were APFS clones of each other, so the space freed was 55.4 GiB (44.2 → 99.6 GiB free), not the
~109 GiB `du` showed. Kept: the zip and every tracked file.

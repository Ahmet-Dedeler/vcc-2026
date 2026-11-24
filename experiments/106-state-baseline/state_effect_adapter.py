#!/usr/bin/env python3
"""Leakage-safe STATE inference adapter for effect-level VCC evaluation.

The adapter reads only non-targeting rows from the held-context H5AD, plus a
frozen incumbent effect bundle and immutable checkpoint manifest. Perturbation
truth is deliberately absent from the API and its expression rows are never
indexed. The result is a frozen scorer-space effect matrix that can be scored
later by a separate process.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import re
import shutil
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
STATE_GENE_PATH = (
    ROOT / "data/external/state-vcc/competition_support_set-log1p/gene_names.csv"
)
PANEL_PATH = ROOT / "data/challenge/pert_counts.csv"

TARGET_SUM = 50_000.0
CONTROL_LABEL = "non-targeting"
BASE_SEED = 20_260_830
N_CONTROL_SETS = 4
SET_LEN = 128
EXPECTED_STEP = 20_000
EXPECTED_PROFILE = "state_sm"
EXPECTED_SOURCE_COMMIT = "9bbfe78a434a55205e4de834e1ea99f85f7a3add"
EXPECTED_INPUT_DIM = 18_080
EXPECTED_PERT_DIM = 5_120
EXPECTED_HIDDEN_DIM = 672
EXPECTED_CONTROL_ROWS = 18_400
MAX_LOG_EXPRESSION = 14.0
HELD_OPERATOR_PATHS = {
    "hct116": HERE / "operators/held_orion_hct116_pruned.npz",
    "hek293t": HERE / "operators/held_orion_hek293t_pruned.npz",
}


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_strings(values: Sequence[str]) -> str:
    return hashlib.sha256(("\n".join(map(str, values)) + "\n").encode()).hexdigest()


def tree_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    files = sorted(
        item
        for item in path.rglob("*")
        if item.is_file()
        and "__pycache__" not in item.parts
        and item.suffix not in {".pyc", ".pyo"}
    )
    for item in files:
        relative = item.relative_to(path).as_posix().encode()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        with item.open("rb") as handle:
            while block := handle.read(8 << 20):
                digest.update(block)
    return digest.hexdigest()


def runtime_environment() -> dict[str, str]:
    packages = sorted(
        f"{dist.metadata.get('Name', 'UNKNOWN')}=={dist.version}"
        for dist in importlib.metadata.distributions()
    )
    return {
        "python_executable": str(Path(sys.executable).resolve()),
        "python_version": platform.python_version(),
        "torch_version": importlib.metadata.version("torch"),
        "transformers_version": importlib.metadata.version("transformers"),
        "packages_sha256": hashlib.sha256(
            ("\n".join(packages) + "\n").encode()
        ).hexdigest(),
    }


def stable_seed(*parts: object) -> int:
    payload = "\0".join(map(str, parts)).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def read_axis(path: Path, preferred_header: str | None = None) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows:
        raise ValueError(f"empty axis file: {path}")
    first = rows[0][0].strip()
    headers = {"gene_name", "target_gene", "gene", "symbol"}
    if preferred_header is not None and first == preferred_header:
        rows = rows[1:]
    elif first in headers:
        rows = rows[1:]
    values = tuple(row[0].strip() for row in rows if row and row[0].strip())
    if not values:
        raise ValueError(f"axis contains no values: {path}")
    if len(values) != len(set(values)):
        raise ValueError(f"axis contains duplicate symbols: {path}")
    return values


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def atomic_npz(path: Path, **values: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with tmp.open("wb") as handle:
        np.savez_compressed(handle, **values)
    tmp.replace(path)


@dataclass(frozen=True)
class CheckpointSpec:
    path: Path
    sha256: str
    global_step: int
    profile: str
    source_commit: str
    source_patch_sha256: str
    source_tree: str
    source_tree_sha256: str
    config_sha256: str
    hparams_sha256: str
    var_dims_sha256: str
    pert_map_sha256: str
    batch_map_sha256: str
    cell_type_map_sha256: str
    uv_lock_sha256: str
    support_archive_sha256: str
    state_gene_axis_sha256: str
    starter_toml_sha256: str
    runtime_packages_sha256: str
    runtime: Mapping[str, object]
    training_recipe: Mapping[str, object]
    support_files: Mapping[str, object]

    @classmethod
    def load(cls, path: Path) -> "CheckpointSpec":
        raw = json.loads(path.read_text())
        checkpoint = (path.parent / raw["checkpoint"]).resolve()
        spec = cls(
            path=checkpoint,
            sha256=str(raw["checkpoint_sha256"]),
            global_step=int(raw["global_step"]),
            profile=str(raw["profile"]),
            source_commit=str(raw["source_commit"]),
            source_patch_sha256=str(raw["source_patch_sha256"]),
            source_tree=str(raw["source_tree"]),
            source_tree_sha256=str(raw["source_tree_sha256"]),
            config_sha256=str(raw["config_sha256"]),
            hparams_sha256=str(raw["hparams_sha256"]),
            var_dims_sha256=str(raw["var_dims_sha256"]),
            pert_map_sha256=str(raw["pert_map_sha256"]),
            batch_map_sha256=str(raw["batch_map_sha256"]),
            cell_type_map_sha256=str(raw["cell_type_map_sha256"]),
            uv_lock_sha256=str(raw["uv_lock_sha256"]),
            support_archive_sha256=str(raw["support_archive_sha256"]),
            state_gene_axis_sha256=str(raw["state_gene_axis_sha256"]),
            starter_toml_sha256=str(raw["starter_toml_sha256"]),
            runtime_packages_sha256=str(raw["runtime_packages_sha256"]),
            runtime=dict(raw["runtime"]),
            training_recipe=dict(raw["training_recipe"]),
            support_files=dict(raw["support_files"]),
        )
        spec.validate()
        return spec

    def validate(self) -> None:
        if self.global_step != EXPECTED_STEP:
            raise ValueError(
                f"checkpoint step {self.global_step} != fixed step {EXPECTED_STEP}"
            )
        if self.profile != EXPECTED_PROFILE:
            raise ValueError(f"checkpoint profile {self.profile!r} is not state_sm")
        if self.source_commit != EXPECTED_SOURCE_COMMIT:
            raise ValueError("STATE source commit changed")
        if self.training_recipe.get("max_steps") != EXPECTED_STEP:
            raise ValueError("training recipe does not bind the fixed step")
        if self.training_recipe.get("train_seed") != 42:
            raise ValueError("training recipe does not bind train_seed=42")
        if not self.support_files:
            raise ValueError("checkpoint manifest lacks bound support files")
        if not self.path.is_file():
            raise FileNotFoundError(self.path)
        observed = sha256_file(self.path)
        if observed != self.sha256:
            raise ValueError(f"checkpoint hash changed: {observed} != {self.sha256}")


@dataclass(frozen=True)
class IncumbentSpec:
    manifest_path: Path
    frozen_path: Path
    surface: str
    frozen_sha256: str
    controls_h5ad_sha256: str
    controls_manifest_sha256: str
    held_operator: Path
    held_operator_sha256: str
    held_operator_metadata_sha256: str
    target_order_sha256: str
    gene_axis_sha256: str

    @classmethod
    def load(
        cls,
        manifest_path: Path,
        frozen_path: Path,
        controls_h5ad: Path,
        expected_surface: str,
    ) -> "IncumbentSpec":
        raw = json.loads(manifest_path.read_text())
        expected = {
            "experiment": "106-state-baseline",
            "stage": "held-incumbent-freeze",
            "surface": expected_surface,
            "sealed_before_perturbation_truth": True,
            "endpoint": "a14-shaped",
            "alpha": 0.982,
            "cd4_fraction": 0.05,
            "basis_prune": "essential",
        }
        for key, value in expected.items():
            if raw.get(key) != value:
                raise ValueError(
                    f"incumbent manifest {key}={raw.get(key)!r}, expected {value!r}"
                )
        if frozen_path.resolve() != (manifest_path.parent / "frozen.npz").resolve():
            raise ValueError("incumbent NPZ is not adjacent to its manifest")
        operator_path = (ROOT / str(raw["held_operator"])).resolve()
        expected_operator = HELD_OPERATOR_PATHS[expected_surface].resolve()
        if operator_path != expected_operator:
            raise ValueError("incumbent manifest names the wrong held operator")
        operator_metadata = operator_path.with_suffix(".json")
        spec = cls(
            manifest_path=manifest_path.resolve(),
            frozen_path=frozen_path.resolve(),
            surface=str(raw["surface"]),
            frozen_sha256=str(raw["frozen_sha256"]),
            controls_h5ad_sha256=str(raw["controls_h5ad_sha256"]),
            controls_manifest_sha256=str(raw["controls_manifest_sha256"]),
            held_operator=operator_path,
            held_operator_sha256=str(raw["held_operator_sha256"]),
            held_operator_metadata_sha256=str(raw["held_operator_metadata_sha256"]),
            target_order_sha256=str(raw["target_order_sha256"]),
            gene_axis_sha256=str(raw["gene_axis_sha256"]),
        )
        if sha256_file(spec.frozen_path) != spec.frozen_sha256:
            raise ValueError("incumbent NPZ hash differs from its manifest")
        if sha256_file(controls_h5ad) != spec.controls_h5ad_sha256:
            raise ValueError("control-only H5AD differs from incumbent manifest")
        controls_manifest = controls_h5ad.with_suffix(
            controls_h5ad.suffix + ".manifest.json"
        )
        if sha256_file(controls_manifest) != spec.controls_manifest_sha256:
            raise ValueError("control-only manifest differs from incumbent manifest")
        if sha256_file(spec.held_operator) != spec.held_operator_sha256:
            raise ValueError("held operator differs from incumbent manifest")
        if (
            not operator_metadata.is_file()
            or sha256_file(operator_metadata) != spec.held_operator_metadata_sha256
        ):
            raise ValueError("held operator metadata differs from incumbent manifest")
        control_metadata = json.loads(controls_manifest.read_text())
        if control_metadata.get("stage") != "held-control-only-input":
            raise ValueError("control input manifest has the wrong stage")
        if control_metadata.get("surface") != expected_surface:
            raise ValueError("control input manifest has the wrong surface")
        if control_metadata.get("contains_perturbation_rows") is not False:
            raise ValueError("control input manifest permits perturbation rows")
        if control_metadata.get("control_rows") != EXPECTED_CONTROL_ROWS:
            raise ValueError("control input manifest has the wrong row count")
        if control_metadata.get("expression_read_policy") != (
            "only predicate-filtered control rows are materialized by experiment code"
        ):
            raise ValueError("control input manifest has the wrong expression policy")
        expected_labels = {
            "hct116": {"Non-Targeting": 165_777},
            "hek293t": {"Non-Targeting": 218_838},
        }[expected_surface]
        if control_metadata.get("accepted_control_label_counts") != expected_labels:
            raise ValueError("control input manifest has the wrong accepted labels")
        return spec


def validate_runtime_provenance(
    checkpoint: CheckpointSpec,
    run_dir: Path,
    state_genes: Sequence[str],
) -> None:
    run_dir = run_dir.resolve()
    if run_dir != checkpoint.path.parent.resolve():
        raise ValueError("runtime directory differs from checkpoint manifest directory")
    frozen = {
        "config.yaml": checkpoint.config_sha256,
        "hparams.yaml": checkpoint.hparams_sha256,
        "var_dims.pkl": checkpoint.var_dims_sha256,
        "pert_onehot_map.pt": checkpoint.pert_map_sha256,
        "batch_onehot_map.torch": checkpoint.batch_map_sha256,
        "cell_type_onehot_map.torch": checkpoint.cell_type_map_sha256,
        "uv.lock": checkpoint.uv_lock_sha256,
        "starter.local.toml": checkpoint.starter_toml_sha256,
        "runtime-packages.txt": checkpoint.runtime_packages_sha256,
    }
    for name, expected in frozen.items():
        path = run_dir / name
        if not path.is_file():
            raise FileNotFoundError(path)
        observed = sha256_file(path)
        if observed != expected:
            raise ValueError(f"frozen runtime artifact changed: {name}")
    source_root = (run_dir / checkpoint.source_tree).resolve()
    if run_dir not in source_root.parents:
        raise ValueError("frozen STATE source escapes the checkpoint directory")
    if not source_root.is_dir():
        raise FileNotFoundError(source_root)
    if tree_sha256(source_root) != checkpoint.source_tree_sha256:
        raise ValueError("frozen STATE source tree changed")
    if sha256_strings(state_genes) != checkpoint.state_gene_axis_sha256:
        raise ValueError("STATE gene axis changed")
    observed_runtime = runtime_environment()
    for key, expected in checkpoint.runtime.items():
        if observed_runtime.get(key) != expected:
            raise ValueError(
                f"runtime {key}={observed_runtime.get(key)!r}, expected {expected!r}"
            )
    if observed_runtime["packages_sha256"] != checkpoint.runtime_packages_sha256:
        raise ValueError("installed package set differs from frozen runtime")


@dataclass(frozen=True)
class GeneAlignment:
    source_genes: tuple[str, ...]
    state_genes: tuple[str, ...]
    source_cols: np.ndarray
    state_cols: np.ndarray
    source_only: tuple[str, ...]
    state_only: tuple[str, ...]

    @classmethod
    def build(
        cls, source_genes: Sequence[str], state_genes: Sequence[str]
    ) -> "GeneAlignment":
        source = tuple(map(str, source_genes))
        state = tuple(map(str, state_genes))
        if len(source) != len(set(source)):
            raise ValueError("source gene axis contains duplicates")
        if len(state) != len(set(state)):
            raise ValueError("STATE gene axis contains duplicates")
        source_index = {gene: index for index, gene in enumerate(source)}
        source_cols = []
        state_cols = []
        for state_col, gene in enumerate(state):
            source_col = source_index.get(gene)
            if source_col is not None:
                source_cols.append(source_col)
                state_cols.append(state_col)
        source_set = set(source)
        state_set = set(state)
        return cls(
            source_genes=source,
            state_genes=state,
            source_cols=np.asarray(source_cols, dtype=np.int64),
            state_cols=np.asarray(state_cols, dtype=np.int64),
            source_only=tuple(gene for gene in source if gene not in state_set),
            state_only=tuple(gene for gene in state if gene not in source_set),
        )

    def align_and_log_controls(self, raw: object) -> np.ndarray:
        import scipy.sparse as sp

        if raw.shape[1] != len(self.source_genes):
            raise ValueError("control matrix does not match source gene axis")
        shared = raw[:, self.source_cols]
        if sp.issparse(shared):
            shared = shared.toarray()
        shared = np.asarray(shared, dtype=np.float32)
        if not np.isfinite(shared).all() or np.any(shared < 0):
            raise ValueError("controls must be finite and nonnegative")
        out = np.zeros((shared.shape[0], len(self.state_genes)), dtype=np.float32)
        out[:, self.state_cols] = shared
        np.log1p(out, out=out)
        return out

    def project_mean_to_source(
        self,
        state_mean: np.ndarray,
        control_mean: np.ndarray,
    ) -> np.ndarray:
        state_mean = np.asarray(state_mean, dtype=np.float64)
        control_mean = np.asarray(control_mean, dtype=np.float64)
        if state_mean.shape[-1] != len(self.state_genes):
            raise ValueError("STATE mean does not match STATE gene axis")
        if control_mean.shape != (len(self.source_genes),):
            raise ValueError("control mean does not match source gene axis")
        out = np.broadcast_to(
            control_mean, state_mean.shape[:-1] + control_mean.shape
        ).copy()
        predicted_shared = state_mean[..., self.state_cols]
        control_shared = control_mean[self.source_cols]
        predicted_mass = predicted_shared.sum(axis=-1, keepdims=True)
        control_mass = float(control_shared.sum())
        if not np.isfinite(predicted_mass).all() or np.any(predicted_mass <= 0):
            raise ValueError("predicted shared-gene mass must be positive and finite")
        predicted_shared = predicted_shared * (control_mass / predicted_mass)
        out[..., self.source_cols] = predicted_shared
        return out


def sample_control_sets(
    control_rows: np.ndarray,
    context: str,
    *,
    base_seed: int = BASE_SEED,
    n_sets: int = N_CONTROL_SETS,
    set_len: int = SET_LEN,
) -> np.ndarray:
    rows = np.asarray(control_rows, dtype=np.int64)
    if rows.ndim != 1 or len(rows) == 0:
        raise ValueError("control row set is empty or malformed")
    if n_sets <= 0 or set_len <= 0:
        raise ValueError("n_sets and set_len must be positive")
    rng = np.random.default_rng(stable_seed(base_seed, context))
    return rng.choice(rows, size=(n_sets, set_len), replace=True)


def build_perturbation_batch(
    embedding: np.ndarray, n_sets: int, set_len: int = SET_LEN
) -> np.ndarray:
    embedding = np.asarray(embedding, dtype=np.float32)
    if embedding.shape != (EXPECTED_PERT_DIM,):
        raise ValueError(
            f"perturbation embedding shape {embedding.shape} != ({EXPECTED_PERT_DIM},)"
        )
    if not np.isfinite(embedding).all() or float(np.linalg.norm(embedding)) <= 0:
        raise ValueError("perturbation embedding must be finite and nonzero")
    return (
        np.broadcast_to(embedding, (n_sets, set_len, EXPECTED_PERT_DIM))
        .reshape(n_sets * set_len, EXPECTED_PERT_DIM)
        .copy()
    )


def streamed_mean_counts(prediction_blocks: Sequence[np.ndarray]) -> np.ndarray:
    total: np.ndarray | None = None
    n_rows = 0
    for block in prediction_blocks:
        values = np.asarray(block, dtype=np.float64)
        if values.ndim != 2:
            raise ValueError("prediction blocks must be two-dimensional")
        if not np.isfinite(values).all():
            raise ValueError("STATE emitted non-finite log-expression")
        clipped = np.clip(values, 0.0, MAX_LOG_EXPRESSION)
        counts = np.expm1(clipped)
        current = counts.sum(axis=0, dtype=np.float64)
        total = current if total is None else total + current
        n_rows += len(values)
    if total is None or n_rows == 0:
        raise ValueError("no prediction rows")
    return total / n_rows


def bulk_profile(mean_counts: np.ndarray, target_sum: float = TARGET_SUM) -> np.ndarray:
    values = np.asarray(mean_counts, dtype=np.float64)
    if values.ndim not in (1, 2):
        raise ValueError("mean counts must be one- or two-dimensional")
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("mean counts must be finite and nonnegative")
    total = values.sum(axis=-1, keepdims=True)
    if np.any(total <= 0):
        raise ValueError("mean counts have zero total")
    return np.log1p(target_sum * values / total)


def scorer_effect(predicted_mean: np.ndarray, control_mean: np.ndarray) -> np.ndarray:
    return bulk_profile(predicted_mean) - bulk_profile(control_mean)


def apply_coordinate_mask(
    effect: np.ndarray, genes: Sequence[str], excluded: set[str]
) -> np.ndarray:
    out = np.asarray(effect, dtype=np.float64).copy()
    mask = np.asarray([str(gene) in excluded for gene in genes], dtype=bool)
    out[..., mask] = 0.0
    return out


def load_target_embeddings(
    targets: Sequence[str], pert_map: Mapping[object, object]
) -> np.ndarray:
    lookup = {str(key): value for key, value in pert_map.items()}
    rows = []
    for target in targets:
        if str(target) not in lookup:
            raise KeyError(f"target missing from saved perturbation map: {target}")
        value = lookup[str(target)]
        if hasattr(value, "detach"):
            value = value.detach().cpu().numpy()
        array = np.asarray(value, dtype=np.float32)
        if array.shape != (EXPECTED_PERT_DIM,):
            raise ValueError(f"{target}: embedding shape {array.shape}")
        if not np.isfinite(array).all() or float(np.linalg.norm(array)) <= 0:
            raise ValueError(f"{target}: embedding must be finite and nonzero")
        rows.append(array)
    return np.stack(rows)


def validate_model_contract(model: object, state_genes: Sequence[str]) -> None:
    expected = {
        "input_dim": EXPECTED_INPUT_DIM,
        "output_dim": EXPECTED_INPUT_DIM,
        "pert_dim": EXPECTED_PERT_DIM,
        "cell_sentence_len": SET_LEN,
        "hidden_dim": EXPECTED_HIDDEN_DIM,
        "output_space": "all",
        "predict_residual": True,
    }
    for name, value in expected.items():
        observed = getattr(model, name, None)
        if observed != value:
            raise ValueError(f"model contract {name}={observed!r}, expected {value!r}")
    if getattr(model, "batch_encoder", None) is not None:
        raise ValueError("state_sm checkpoint unexpectedly uses a batch encoder")
    if getattr(model, "gene_decoder", None) is not None:
        raise ValueError("state_sm checkpoint unexpectedly uses a gene decoder")
    hparams = getattr(model, "hparams", {})
    checkpoint_genes = tuple(map(str, hparams.get("gene_names", ())))
    if checkpoint_genes != tuple(map(str, state_genes)):
        raise ValueError("checkpoint gene_names differ from frozen STATE axis")


def sum_control_counts(
    x: object, rows: np.ndarray, n_genes: int, chunk_rows: int = 256
) -> np.ndarray:
    import scipy.sparse as sp

    total = np.zeros(n_genes, dtype=np.float64)
    for start in range(0, len(rows), chunk_rows):
        block = x[rows[start : start + chunk_rows], :]
        if sp.issparse(block):
            total += np.asarray(block.sum(axis=0)).ravel()
        else:
            total += np.asarray(block, dtype=np.float64).sum(axis=0)
    if not np.isfinite(total).all() or np.any(total < 0) or total.sum() <= 0:
        raise ValueError("invalid aggregate control counts")
    return total


def load_runtime(
    checkpoint: CheckpointSpec,
    run_dir: Path,
    state_genes: Sequence[str],
    device_name: str,
) -> tuple[object, object, object]:
    os.environ.setdefault("STATE_USE_MPS", "1")
    os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "0")
    frozen_source_root = (run_dir / checkpoint.source_tree).parent
    sys.path.insert(0, str(frozen_source_root))

    import torch
    from state.tx.models.state_transition import StateTransitionPerturbationModel

    if device_name == "mps":
        if not torch.backends.mps.is_available():
            raise RuntimeError("MPS requested but unavailable")
        device = torch.device("mps")
    elif device_name == "cpu":
        device = torch.device("cpu")
    else:
        raise ValueError("device must be mps or cpu")

    pert_map_path = run_dir / "pert_onehot_map.pt"
    if sha256_file(pert_map_path) != checkpoint.pert_map_sha256:
        raise ValueError("saved perturbation map hash changed")
    pert_map = torch.load(pert_map_path, map_location="cpu", weights_only=False)
    model = StateTransitionPerturbationModel.load_from_checkpoint(
        checkpoint.path,
        map_location="cpu",
        weights_only=False,
    )
    validate_model_contract(model, state_genes)
    model.to(device).eval()
    return model, pert_map, torch


def infer_effects(
    *,
    controls_h5ad: Path,
    incumbent_npz: Path,
    incumbent_manifest: Path,
    context: str,
    checkpoint_manifest: Path,
    run_dir: Path,
    output_dir: Path,
    evaluation_run_id: str,
    n_sets: int,
    base_seed: int,
    device_name: str,
    extra_mask_paths: Sequence[Path],
) -> None:
    import anndata as ad

    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {output_dir}")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,127}", evaluation_run_id):
        raise ValueError("evaluation run id is malformed")
    checkpoint = CheckpointSpec.load(checkpoint_manifest)
    expected_surface = context.lower()
    incumbent_spec = IncumbentSpec.load(
        incumbent_manifest,
        incumbent_npz,
        controls_h5ad,
        expected_surface,
    )
    state_genes = read_axis(STATE_GENE_PATH)
    if len(state_genes) != EXPECTED_INPUT_DIM:
        raise ValueError("frozen STATE axis length changed")
    validate_runtime_provenance(checkpoint, run_dir, state_genes)

    with np.load(incumbent_npz, allow_pickle=False) as incumbent:
        source_genes = tuple(incumbent["genes"].astype(str))
        targets = tuple(incumbent["targets"].astype(str))
        incumbent_effect = incumbent["baseline"].astype(np.float32)
        incumbent_bctrl = incumbent["bctrl"].astype(np.float32)
    if sha256_strings(targets) != incumbent_spec.target_order_sha256:
        raise ValueError("incumbent target order differs from its manifest")
    if sha256_strings(source_genes) != incumbent_spec.gene_axis_sha256:
        raise ValueError("incumbent gene axis differs from its manifest")
    if incumbent_effect.shape != (len(targets), len(source_genes)):
        raise ValueError("incumbent effect shape is malformed")
    if incumbent_bctrl.shape != (len(source_genes),):
        raise ValueError("incumbent control profile shape is malformed")
    for name, values in {
        "incumbent_effect": incumbent_effect,
        "incumbent_bctrl": incumbent_bctrl,
    }.items():
        if not np.isfinite(values).all():
            raise ValueError(f"{name} contains non-finite values")

    alignment = GeneAlignment.build(source_genes, state_genes)
    handle = ad.read_h5ad(controls_h5ad, backed="r")
    try:
        observed_genes = tuple(map(str, handle.var_names))
        if observed_genes != source_genes:
            raise ValueError("control and incumbent gene axes differ")
        if "target_gene" not in handle.obs:
            raise KeyError("controls file lacks obs['target_gene']")
        labels = np.asarray(handle.obs["target_gene"].astype(str))
        control_rows = np.flatnonzero(labels == CONTROL_LABEL)
        noncontrol_rows = int(np.count_nonzero(labels != CONTROL_LABEL))
        if len(control_rows) == 0:
            raise ValueError("no non-targeting controls")
        if len(control_rows) != EXPECTED_CONTROL_ROWS:
            raise ValueError(
                f"STATE inference requires exactly {EXPECTED_CONTROL_ROWS} controls; "
                f"found {len(control_rows)}"
            )
        if noncontrol_rows:
            raise ValueError(
                "STATE inference requires a control-only H5AD; "
                f"found {noncontrol_rows} non-control rows"
            )
        selected = sample_control_sets(
            control_rows,
            context,
            base_seed=base_seed,
            n_sets=n_sets,
            set_len=SET_LEN,
        )
        raw_selected = handle.X[selected.reshape(-1), :]
        control_sets = alignment.align_and_log_controls(raw_selected).reshape(
            n_sets, SET_LEN, len(state_genes)
        )
        control_sum = sum_control_counts(handle.X, control_rows, len(source_genes))
        control_mean = control_sum / len(control_rows)
    finally:
        handle.file.close()

    model, pert_map, torch = load_runtime(checkpoint, run_dir, state_genes, device_name)
    computed_bctrl = bulk_profile(control_mean).astype(np.float32)
    bctrl_max_abs_diff = float(
        np.max(np.abs(computed_bctrl.astype(np.float64) - incumbent_bctrl))
    )
    if bctrl_max_abs_diff > 2e-5:
        raise ValueError(
            "incumbent and control-file baselines differ: "
            f"max_abs_diff={bctrl_max_abs_diff:.6g}"
        )
    ctrl_cpm = 1e6 * control_mean / control_mean.sum()
    embeddings = load_target_embeddings(targets, pert_map)
    device = next(model.parameters()).device
    state_effects = np.zeros((len(targets), len(source_genes)), dtype=np.float32)
    predicted_means = np.zeros_like(state_effects)
    source_only_set = set(alignment.source_only)
    source_only_cols = np.asarray(
        [gene in source_only_set for gene in source_genes], dtype=bool
    )
    output_min = np.inf
    output_max = -np.inf
    clipped_low_values = 0
    clipped_high_values = 0

    controls_flat = control_sets.reshape(n_sets * SET_LEN, len(state_genes))
    control_tensor = torch.from_numpy(controls_flat).to(device)
    with torch.inference_mode():
        for row, target in enumerate(targets):
            pert = build_perturbation_batch(embeddings[row], n_sets)
            batch = {
                "ctrl_cell_emb": control_tensor,
                "pert_emb": torch.from_numpy(pert).to(device),
            }
            pred = model(batch, padded=True)
            pred_np = pred.detach().cpu().numpy().astype(np.float32, copy=False)
            if pred_np.shape != (n_sets * SET_LEN, len(state_genes)):
                raise ValueError(f"{target}: unexpected output shape {pred_np.shape}")
            if not np.isfinite(pred_np).all():
                raise ValueError(f"{target}: STATE output contains non-finite values")
            output_min = min(output_min, float(pred_np.min()))
            output_max = max(output_max, float(pred_np.max()))
            clipped_low_values += int(np.count_nonzero(pred_np < 0.0))
            clipped_high_values += int(np.count_nonzero(pred_np > MAX_LOG_EXPRESSION))
            state_mean = streamed_mean_counts([pred_np])
            full_mean = alignment.project_mean_to_source(state_mean, control_mean)
            predicted_means[row] = full_mean.astype(np.float32)
            effect = scorer_effect(full_mean, control_mean)
            effect[source_only_cols] = 0.0
            state_effects[row] = effect.astype(np.float32)

    excluded = set(read_axis(PANEL_PATH, preferred_header="target_gene"))
    for path in extra_mask_paths:
        excluded.update(read_axis(path))
    state_effects = apply_coordinate_mask(state_effects, source_genes, excluded).astype(
        np.float32
    )

    if np.any(state_effects[:, source_only_cols] != 0):
        raise AssertionError("source-only genes acquired nonzero STATE effect")
    masked_cols = np.asarray([gene in excluded for gene in source_genes], dtype=bool)
    if np.any(state_effects[:, masked_cols] != 0):
        raise AssertionError("masked coordinates are not exactly zero")
    for name, values in {
        "control_mean": control_mean,
        "ctrl_cpm": ctrl_cpm,
        "state_effects": state_effects,
        "predicted_means": predicted_means,
    }.items():
        if not np.isfinite(values).all():
            raise ValueError(f"{name} contains non-finite values")

    final_output = output_dir
    final_output.parent.mkdir(parents=True, exist_ok=True)
    output_dir = Path(
        tempfile.mkdtemp(prefix=f".{final_output.name}.", dir=final_output.parent)
    )
    frozen_path = output_dir / "frozen.npz"
    atomic_npz(
        frozen_path,
        genes=np.asarray(source_genes, dtype="U"),
        targets=np.asarray(targets, dtype="U"),
        bctrl=incumbent_bctrl,
        control_mean=control_mean.astype(np.float64),
        ctrl_cpm=ctrl_cpm.astype(np.float64),
        state_supported=np.asarray(
            [gene not in source_only_set for gene in source_genes], dtype=bool
        ),
        sampled_control_rows=selected.astype(np.int64),
        incumbent=incumbent_effect,
        state=state_effects,
        state_predicted_mean=predicted_means,
    )
    manifest = {
        "experiment": "106-state-baseline",
        "evaluation_run_id": evaluation_run_id,
        "stage": "effect-freeze",
        "sealed_before_perturbation_truth": True,
        "inference_api_has_no_truth_argument": True,
        "expression_rows_read": "non-targeting controls only",
        "noncontrol_rows_present_but_not_indexed": noncontrol_rows,
        "context": context,
        "surface": expected_surface,
        "controls_h5ad": str(controls_h5ad.relative_to(ROOT)),
        "controls_h5ad_sha256": sha256_file(controls_h5ad),
        "controls_manifest_sha256": sha256_file(
            controls_h5ad.with_suffix(controls_h5ad.suffix + ".manifest.json")
        ),
        "incumbent_npz": str(incumbent_npz.relative_to(ROOT)),
        "incumbent_npz_sha256": sha256_file(incumbent_npz),
        "incumbent_manifest": str(incumbent_manifest.relative_to(ROOT)),
        "incumbent_manifest_sha256": sha256_file(incumbent_manifest),
        "held_operator": str(incumbent_spec.held_operator.relative_to(ROOT)),
        "held_operator_sha256": incumbent_spec.held_operator_sha256,
        "held_operator_metadata_sha256": (incumbent_spec.held_operator_metadata_sha256),
        "checkpoint_manifest": str(checkpoint_manifest.relative_to(ROOT)),
        "checkpoint_manifest_sha256": sha256_file(checkpoint_manifest),
        "checkpoint_sha256": checkpoint.sha256,
        "global_step": checkpoint.global_step,
        "profile": checkpoint.profile,
        "targets": len(targets),
        "target_order_sha256": sha256_strings(targets),
        "source_genes": len(source_genes),
        "state_genes": len(state_genes),
        "shared_genes": int(len(alignment.source_cols)),
        "source_only_genes": len(alignment.source_only),
        "state_only_genes": list(alignment.state_only),
        "source_axis_sha256": sha256_strings(source_genes),
        "state_axis_sha256": sha256_strings(state_genes),
        "base_seed": base_seed,
        "control_sets": n_sets,
        "set_len": SET_LEN,
        "sampled_control_rows_sha256": hashlib.sha256(selected.tobytes()).hexdigest(),
        "normalization": "log1p raw counts input; expm1 then 50k log-pseudobulk",
        "incumbent_bctrl_max_abs_diff": bctrl_max_abs_diff,
        "shared_mass_matching": True,
        "source_only_policy": "preserve control expression; exact zero STATE effect",
        "output_bounds": {
            "required_min": 0.0,
            "required_max": MAX_LOG_EXPRESSION,
            "observed_min": output_min,
            "observed_max": output_max,
            "clip_policy": "np.clip matching pinned STATE tx infer",
            "clipped_low_values": clipped_low_values,
            "clipped_high_values": clipped_high_values,
            "clip_applied": bool(clipped_low_values or clipped_high_values),
        },
        "masked_symbols": sorted(excluded),
        "masked_symbols_sha256": sha256_strings(sorted(excluded)),
        "frozen_sha256": sha256_file(frozen_path),
        "implementation_sha256": sha256_file(Path(__file__)),
    }
    try:
        atomic_json(output_dir / "manifest.json", manifest)
        os.replace(output_dir, final_output)
    except BaseException:
        shutil.rmtree(output_dir, ignore_errors=True)
        raise
    print(final_output / "manifest.json")


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--controls-h5ad", required=True, type=Path)
    p.add_argument("--incumbent-npz", required=True, type=Path)
    p.add_argument("--incumbent-manifest", required=True, type=Path)
    p.add_argument("--context", required=True)
    p.add_argument("--checkpoint-manifest", required=True, type=Path)
    p.add_argument(
        "--run-dir",
        type=Path,
        default=HERE / "runs/official-state_sm",
    )
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--evaluation-run-id", required=True)
    p.add_argument("--control-sets", type=int, default=N_CONTROL_SETS)
    p.add_argument("--base-seed", type=int, default=BASE_SEED)
    p.add_argument("--device", choices=("mps", "cpu"), default="mps")
    p.add_argument("--extra-mask", action="append", default=[], type=Path)
    return p


def main() -> None:
    args = parser().parse_args()
    infer_effects(
        controls_h5ad=args.controls_h5ad.resolve(),
        incumbent_npz=args.incumbent_npz.resolve(),
        incumbent_manifest=args.incumbent_manifest.resolve(),
        context=args.context,
        checkpoint_manifest=args.checkpoint_manifest.resolve(),
        run_dir=args.run_dir.resolve(),
        output_dir=args.output_dir.resolve(),
        evaluation_run_id=args.evaluation_run_id,
        n_sets=args.control_sets,
        base_seed=args.base_seed,
        device_name=args.device,
        extra_mask_paths=tuple(path.resolve() for path in args.extra_mask),
    )


if __name__ == "__main__":
    main()

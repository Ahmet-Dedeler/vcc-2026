#!/usr/bin/env python3
"""Freeze the completed fixed-step STATE checkpoint and its inference inputs."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import pickle
import platform
import shutil
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
STATE_ROOT = ROOT / "data/external/state-vcc/state"
DEFAULT_RUN = HERE / "runs/official-state_sm"
DEFAULT_OUTPUT = HERE / "frozen/state_sm-step20000"
EXPECTED_STEP = 20_000
EXPECTED_PROFILE = "state_sm"
EXPECTED_COMMIT = "9bbfe78a434a55205e4de834e1ea99f85f7a3add"
SUPPORT_ROOT = ROOT / "data/external/state-vcc/competition_support_set-log1p"
SUPPORT_FILES = (
    "competition_train.h5",
    "k562_gwps.h5",
    "rpe1.h5",
    "jurkat.h5",
    "k562.h5",
    "hepg2.h5",
    "gene_names.csv",
    "ESM2_pert_features.pt",
)


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def atomic_copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = destination.with_name(f".{destination.name}.{os.getpid()}.tmp")
    shutil.copyfile(source, tmp)
    os.replace(tmp, destination)


def atomic_json(path: Path, value: object) -> None:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def state_training_processes(run_dir: Path) -> list[str]:
    result = subprocess.run(
        ["ps", "-axo", "pid=,command="],
        check=True,
        capture_output=True,
        text=True,
    )
    resolved = run_dir.resolve()
    needle = str(resolved)
    output_dir_arg = f"output_dir={resolved.parent}"
    name_arg = f"name={resolved.name}"
    return [
        line.strip()
        for line in result.stdout.splitlines()
        if ("state tx train" in line or "run_with_memory_guard.py" in line)
        and (needle in line or (output_dir_arg in line and name_arg in line))
    ]


def git_head(path: Path) -> str:
    return subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def git_diff_sha256(path: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(path), "diff", "--binary"],
        check=True,
        capture_output=True,
    )
    return hashlib.sha256(result.stdout).hexdigest()


def git_untracked_source_files(path: Path) -> tuple[str, ...]:
    result = subprocess.run(
        [
            "git",
            "-C",
            str(path),
            "status",
            "--porcelain=v1",
            "--untracked-files=all",
            "--",
            "src/state",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    return tuple(
        line[3:] for line in result.stdout.splitlines() if line.startswith("?? ")
    )


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


def runtime_environment() -> dict[str, object]:
    packages = sorted(
        f"{dist.metadata.get('Name', 'UNKNOWN')}=={dist.version}"
        for dist in importlib.metadata.distributions()
    )
    payload = ("\n".join(packages) + "\n").encode()
    return {
        "python_executable": str(Path(sys.executable).resolve()),
        "python_version": platform.python_version(),
        "torch_version": importlib.metadata.version("torch"),
        "transformers_version": importlib.metadata.version("transformers"),
        "packages": packages,
        "packages_sha256": hashlib.sha256(payload).hexdigest(),
    }


def hash_bound_file(path: Path) -> dict[str, object]:
    before = path.stat()
    digest = sha256_file(path)
    after = path.stat()
    if (
        before.st_size != after.st_size
        or before.st_mtime_ns != after.st_mtime_ns
        or before.st_ino != after.st_ino
    ):
        raise RuntimeError(f"file changed while hashing: {path}")
    return {
        "path": str(path.relative_to(ROOT)),
        "bytes": before.st_size,
        "sha256": digest,
    }


def validate_training_recipe(
    config_path: Path, starter_path: Path
) -> dict[str, object]:
    import yaml

    config = yaml.safe_load(config_path.read_text())
    expected_config = {
        ("data", "name"): "PerturbationDataModule",
        ("data", "kwargs", "toml_config_path"): str(starter_path.resolve()),
        ("data", "kwargs", "num_workers"): 8,
        ("data", "kwargs", "pin_memory"): False,
        ("data", "kwargs", "batch_col"): "batch_var",
        ("data", "kwargs", "pert_col"): "target_gene",
        ("data", "kwargs", "cell_type_key"): "cell_type",
        ("data", "kwargs", "control_pert"): "non-targeting",
        ("data", "kwargs", "perturbation_features_file"): str(
            (SUPPORT_ROOT / "ESM2_pert_features.pt").resolve()
        ),
        ("training", "batch_size"): 16,
        ("training", "lr"): 0.0001,
        ("training", "max_steps"): EXPECTED_STEP,
        ("training", "train_seed"): 42,
        ("training", "val_freq"): 2000,
        ("training", "ckpt_every_n_steps"): EXPECTED_STEP,
        ("name",): "official-state_sm",
    }
    for keys, expected in expected_config.items():
        observed: object = config
        for key in keys:
            if not isinstance(observed, dict) or key not in observed:
                raise ValueError(f"training config lacks {'.'.join(keys)}")
            observed = observed[key]
        if observed != expected:
            raise ValueError(
                f"training config {'.'.join(keys)}={observed!r}, expected {expected!r}"
            )

    starter = tomllib.loads(starter_path.read_text())
    expected_dataset = str(
        SUPPORT_ROOT / "{competition_train,k562_gwps,rpe1,jurkat,k562,hepg2}.h5"
    )
    expected_starter = {
        ("datasets", "replogle_h1"): expected_dataset,
        ("training", "replogle_h1"): "train",
        ("zeroshot", "replogle_h1.hepg2"): "test",
    }
    for keys, expected in expected_starter.items():
        observed: object = starter
        for key in keys:
            if not isinstance(observed, dict) or key not in observed:
                raise ValueError(f"starter TOML lacks {'.'.join(keys)}")
            observed = observed[key]
        if observed != expected:
            raise ValueError(
                f"starter TOML {'.'.join(keys)}={observed!r}, expected {expected!r}"
            )
    return {
        "config_sha256": sha256_file(config_path),
        "starter_sha256": sha256_file(starter_path),
        "max_steps": EXPECTED_STEP,
        "train_seed": 42,
        "batch_size": 16,
        "learning_rate": 0.0001,
        "zeroshot_context": "hepg2",
    }


def validate_existing_freeze(output_dir: Path) -> None:
    manifest_path = output_dir / "checkpoint-manifest.json"
    if not manifest_path.is_file():
        raise FileExistsError(f"incomplete frozen checkpoint directory: {output_dir}")
    manifest = json.loads(manifest_path.read_text())
    expected = {
        "global_step": EXPECTED_STEP,
        "profile": EXPECTED_PROFILE,
        "source_commit": EXPECTED_COMMIT,
        "immutable": True,
    }
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise ValueError(
                f"existing checkpoint manifest {key}={manifest.get(key)!r}, "
                f"expected {value!r}"
            )
    file_hashes = {
        "checkpoint.ckpt": "checkpoint_sha256",
        "config.yaml": "config_sha256",
        "hparams.yaml": "hparams_sha256",
        "var_dims.pkl": "var_dims_sha256",
        "pert_onehot_map.pt": "pert_map_sha256",
        "batch_onehot_map.torch": "batch_map_sha256",
        "cell_type_onehot_map.torch": "cell_type_map_sha256",
        "uv.lock": "uv_lock_sha256",
        "starter.local.toml": "starter_toml_sha256",
        "runtime-packages.txt": "runtime_packages_sha256",
    }
    for filename, key in file_hashes.items():
        path = output_dir / filename
        if not path.is_file() or sha256_file(path) != manifest.get(key):
            raise ValueError(f"existing frozen artifact changed: {filename}")
    source_root = output_dir / str(manifest["source_tree"])
    if tree_sha256(source_root) != manifest.get("source_tree_sha256"):
        raise ValueError("existing frozen STATE source tree changed")


def read_checkpoint_metadata(path: Path) -> dict[str, object]:
    import torch

    checkpoint = torch.load(
        path,
        map_location="cpu",
        weights_only=False,
        mmap=True,
    )
    try:
        hparams = dict(checkpoint.get("hyper_parameters", {}))
        genes = tuple(map(str, hparams.get("gene_names", ())))
        return {
            "global_step": int(checkpoint.get("global_step", -1)),
            "input_dim": int(hparams.get("input_dim", -1)),
            "output_dim": int(hparams.get("output_dim", -1)),
            "pert_dim": int(hparams.get("pert_dim", -1)),
            "hidden_dim": int(hparams.get("hidden_dim", -1)),
            "cell_set_len": int(hparams.get("cell_set_len", -1)),
            "output_space": str(hparams.get("output_space", "")),
            "predict_residual": bool(hparams.get("predict_residual", False)),
            "gene_decoder_bool": bool(hparams.get("gene_decoder_bool", True)),
            "batch_encoder": bool(hparams.get("batch_encoder", False)),
            "gene_names": genes,
        }
    finally:
        del checkpoint


def read_csv_axis(path: Path) -> tuple[str, ...]:
    rows = path.read_text().splitlines()
    if rows and rows[0] in {"gene_name", "gene", "symbol"}:
        rows = rows[1:]
    values = tuple(row.strip().split(",")[0] for row in rows if row.strip())
    if len(values) != len(set(values)):
        raise ValueError(f"duplicate gene names in {path}")
    return values


def validate_checkpoint_metadata(
    metadata: dict[str, object], expected_genes: tuple[str, ...]
) -> None:
    expected = {
        "global_step": EXPECTED_STEP,
        "input_dim": 18_080,
        "output_dim": 18_080,
        "pert_dim": 5_120,
        "hidden_dim": 672,
        "cell_set_len": 128,
        "output_space": "all",
        "predict_residual": True,
        "gene_decoder_bool": False,
        "batch_encoder": False,
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(
                f"checkpoint {key}={metadata.get(key)!r}, expected {value!r}"
            )
    if metadata["gene_names"] != expected_genes:
        raise ValueError("checkpoint gene axis differs from frozen support axis")


def freeze(run_dir: Path, output_dir: Path) -> None:
    if output_dir.exists():
        validate_existing_freeze(output_dir)
        print(output_dir / "checkpoint-manifest.json")
        return
    active = state_training_processes(run_dir)
    if active:
        raise RuntimeError(
            "STATE training is still active; checkpoint is mutable:\n"
            + "\n".join(active)
        )
    if git_head(STATE_ROOT) != EXPECTED_COMMIT:
        raise ValueError("STATE source checkout moved from the pinned commit")
    expected_patch_sha256 = sha256_file(HERE / "mps-runtime.patch")
    observed_patch_sha256 = git_diff_sha256(STATE_ROOT)
    if observed_patch_sha256 != expected_patch_sha256:
        raise ValueError(
            "STATE source diff no longer matches the frozen MPS runtime patch"
        )
    untracked_source = git_untracked_source_files(STATE_ROOT)
    if untracked_source:
        raise ValueError(
            "STATE source contains untracked files: " + ", ".join(untracked_source)
        )

    starter_path = run_dir.parent / "starter.local.toml"
    sources = {
        "checkpoint.ckpt": run_dir / "checkpoints/last.ckpt",
        "config.yaml": run_dir / "config.yaml",
        "hparams.yaml": run_dir / "version_0/hparams.yaml",
        "var_dims.pkl": run_dir / "var_dims.pkl",
        "pert_onehot_map.pt": run_dir / "pert_onehot_map.pt",
        "batch_onehot_map.torch": run_dir / "batch_onehot_map.torch",
        "cell_type_onehot_map.torch": run_dir / "cell_type_onehot_map.torch",
        "uv.lock": STATE_ROOT / "uv.lock",
        "starter.local.toml": starter_path,
    }
    for source in sources.values():
        if not source.is_file():
            raise FileNotFoundError(source)

    training_recipe = validate_training_recipe(sources["config.yaml"], starter_path)
    support_files = {
        name: hash_bound_file(SUPPORT_ROOT / name) for name in SUPPORT_FILES
    }
    support_view_manifest = hash_bound_file(HERE / "support-view-manifest.json")
    bundle_manifest = hash_bound_file(HERE / "bundle-manifest.json")
    runtime = runtime_environment()

    final_output = output_dir
    final_output.parent.mkdir(parents=True, exist_ok=True)
    output_dir = Path(
        tempfile.mkdtemp(prefix=f".{final_output.name}.", dir=final_output.parent)
    )
    source_tree_before = tree_sha256(STATE_ROOT / "src/state")
    source_hashes_before = {name: sha256_file(path) for name, path in sources.items()}
    try:
        for name, source in sources.items():
            atomic_copy(source, output_dir / name)
        frozen_source = output_dir / "source/state"
        shutil.copytree(
            STATE_ROOT / "src/state",
            frozen_source,
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo"),
        )
        source_hashes_after = {
            name: sha256_file(path) for name, path in sources.items()
        }
        if source_hashes_after != source_hashes_before:
            raise RuntimeError("a source artifact changed while it was being frozen")
        if tree_sha256(STATE_ROOT / "src/state") != source_tree_before:
            raise RuntimeError("STATE source changed while it was being frozen")
        if tree_sha256(frozen_source) != source_tree_before:
            raise RuntimeError("frozen STATE source differs from its source")
        frozen_hashes = {name: sha256_file(output_dir / name) for name in sources}
        if frozen_hashes != source_hashes_before:
            raise RuntimeError("frozen artifact hashes differ from their sources")
        packages_path = output_dir / "runtime-packages.txt"
        packages_path.write_text("\n".join(runtime["packages"]) + "\n")
        if sha256_file(packages_path) != runtime["packages_sha256"]:
            raise RuntimeError("runtime package snapshot hash mismatch")

        metadata = read_checkpoint_metadata(output_dir / "checkpoint.ckpt")
        state_genes = read_csv_axis(
            ROOT
            / "data/external/state-vcc/competition_support_set-log1p/gene_names.csv"
        )
        validate_checkpoint_metadata(metadata, state_genes)
        with (output_dir / "var_dims.pkl").open("rb") as handle:
            var_dims = pickle.load(handle)
        if tuple(map(str, var_dims.get("gene_names", ()))) != state_genes:
            raise ValueError("var_dims.pkl gene axis differs from checkpoint axis")
        import yaml

        hparams = yaml.safe_load((output_dir / "hparams.yaml").read_text())
        if tuple(map(str, hparams.get("gene_names", ()))) != state_genes:
            raise ValueError("hparams.yaml gene axis differs from checkpoint axis")
        hparam_contract = {
            "input_dim": 18_080,
            "output_dim": 18_080,
            "pert_dim": 5_120,
            "hidden_dim": 672,
            "cell_set_len": 128,
            "output_space": "all",
            "predict_residual": True,
            "gene_decoder_bool": False,
            "batch_encoder": False,
        }
        for key, expected in hparam_contract.items():
            if hparams.get(key) != expected:
                raise ValueError(
                    f"hparams.yaml {key}={hparams.get(key)!r}, expected {expected!r}"
                )

        provenance = json.loads((HERE / "PROVENANCE.json").read_text())
        manifest = {
            "checkpoint": "checkpoint.ckpt",
            "checkpoint_sha256": frozen_hashes["checkpoint.ckpt"],
            "global_step": metadata["global_step"],
            "profile": EXPECTED_PROFILE,
            "source_commit": EXPECTED_COMMIT,
            "source_patch_sha256": observed_patch_sha256,
            "source_tree": "source/state",
            "source_tree_sha256": source_tree_before,
            "config_sha256": frozen_hashes["config.yaml"],
            "hparams_sha256": frozen_hashes["hparams.yaml"],
            "var_dims_sha256": frozen_hashes["var_dims.pkl"],
            "pert_map_sha256": frozen_hashes["pert_onehot_map.pt"],
            "batch_map_sha256": frozen_hashes["batch_onehot_map.torch"],
            "cell_type_map_sha256": frozen_hashes["cell_type_onehot_map.torch"],
            "uv_lock_sha256": frozen_hashes["uv.lock"],
            "starter_toml_sha256": frozen_hashes["starter.local.toml"],
            "training_recipe": training_recipe,
            "support_files": support_files,
            "support_view_manifest": support_view_manifest,
            "bundle_manifest": bundle_manifest,
            "runtime": {
                key: value for key, value in runtime.items() if key != "packages"
            },
            "runtime_packages_sha256": runtime["packages_sha256"],
            "support_archive_sha256": provenance["support_archive"]["sha256"],
            "state_gene_axis_sha256": __import__("hashlib")
            .sha256(("\n".join(state_genes) + "\n").encode())
            .hexdigest(),
            "model_contract": {
                key: value for key, value in metadata.items() if key != "gene_names"
            },
            "source_run": str(run_dir.relative_to(ROOT)),
            "immutable": True,
            "runtime_import": "frozen source snapshot",
            "selection_rule": "fixed step 20000; best.ckpt forbidden",
        }
        atomic_json(output_dir / "checkpoint-manifest.json", manifest)
        os.replace(output_dir, final_output)
    except BaseException:
        shutil.rmtree(output_dir, ignore_errors=True)
        raise
    print(final_output / "checkpoint-manifest.json")


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run-dir", type=Path, default=DEFAULT_RUN)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return p


def main() -> None:
    args = parser().parse_args()
    freeze(args.run_dir.resolve(), args.output_dir.resolve())


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Freeze Experiment 121 null and bootstrap plans before source-value access."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
AUTHORIZATION_PATH = HERE / "AUTHORIZATION_MANIFEST.json"
COMMON_TARGETS_PATH = HERE / "common_targets.txt"
DEFAULT_OUTPUT = HERE / "work" / "preregistered-freeze-v1"


class FreezeError(RuntimeError):
    """A preregistration input, plan, or immutability invariant failed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FreezeError(message)


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def json_bytes(payload: dict[str, Any]) -> bytes:
    return (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"JSON root is not an object: {path}")
    return value


def read_axis(path: Path) -> tuple[str, ...]:
    require(path.read_bytes().endswith(b"\n"), f"axis lacks final newline: {path}")
    values = tuple(path.read_text(encoding="utf-8").splitlines())
    require(all(value and value == value.strip() for value in values), "invalid axis")
    require(len(values) == len(set(values)), f"duplicate labels in {path}")
    return values


def validate_bindings(authorization: dict[str, Any]) -> None:
    sources = authorization["source_bindings"]
    for source in sources.values():
        if "path" in source:
            records = (source,)
        else:
            records = source.values()
        for record in records:
            path = ROOT / record["path"]
            require(path.is_file() and not path.is_symlink(), f"missing input: {path}")
            if path.suffix != ".npy":
                require(
                    sha256_file(path) == record["sha256"],
                    f"text binding changed: {path}",
                )


def validate_target_axis(protocol: dict[str, Any]) -> tuple[str, ...]:
    targets = read_axis(COMMON_TARGETS_PATH)
    altair = read_axis(
        ROOT / "experiments/119-altair-a549/work/compact-source/perts.txt"
    )
    myllia = read_axis(
        ROOT / "experiments/120-myllia-a549/work/compact-donor/perts.txt"
    )
    exact_intersection = tuple(sorted(set(altair) & set(myllia)))
    expected = tuple(protocol["target_axis"]["labels"])
    require(targets == expected, "common_targets.txt differs from protocol")
    require(targets == exact_intersection, "target axis is not the exact intersection")
    require(len(targets) == 20, "target count changed")
    return targets


def generate_derangements(
    count: int,
    target_count: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.Generator(np.random.PCG64(seed))
    identity = np.arange(target_count)
    result = np.empty((count, target_count), dtype=np.int16)
    for replicate in range(count):
        while True:
            permutation = rng.permutation(target_count)
            if np.all(permutation != identity):
                result[replicate] = permutation
                break
    return result


def generate_bootstrap(
    count: int,
    target_count: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.Generator(np.random.PCG64(seed))
    return rng.integers(
        0,
        target_count,
        size=(count, target_count),
        dtype=np.int16,
    )


def validate_plans(
    derangements: np.ndarray,
    bootstrap: np.ndarray,
    target_count: int,
) -> None:
    identity = np.arange(target_count)
    require(derangements.ndim == 2, "derangements must be two-dimensional")
    require(bootstrap.ndim == 2, "bootstrap must be two-dimensional")
    require(
        derangements.shape[1] == target_count,
        "derangement target count changed",
    )
    require(bootstrap.shape[1] == target_count, "bootstrap target count changed")
    require(
        np.all(np.sort(derangements, axis=1) == identity),
        "a null row is not a permutation",
    )
    require(
        np.all(derangements != identity[None, :]),
        "a null row is not a derangement",
    )
    require(
        np.all((bootstrap >= 0) & (bootstrap < target_count)),
        "bootstrap index out of range",
    )


def safe_output(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    require(
        resolved != HERE and resolved.is_relative_to(HERE),
        f"output must be below {HERE}",
    )
    return resolved


def code_bindings() -> dict[str, str]:
    names = (
        "AUTHORIZATION_MANIFEST.json",
        "common_targets.txt",
        "evaluate.py",
        "freeze_plans.py",
        "preregistered_protocol.json",
        "test_evaluate.py",
    )
    return {name: sha256_file(HERE / name) for name in names}


def build(output: Path = DEFAULT_OUTPUT) -> Path:
    output = safe_output(output)
    require(not output.exists(), f"refusing to overwrite frozen output: {output}")
    protocol = read_json(PROTOCOL_PATH)
    authorization = read_json(AUTHORIZATION_PATH)
    require(
        protocol["status"]
        == "frozen_before_either_source_expression_array_was_opened",
        "protocol status changed",
    )
    require(
        authorization["created_before_source_expression_value_access"] is True,
        "authorization timing declaration changed",
    )
    validate_bindings(authorization)
    targets = validate_target_axis(protocol)

    null_spec = protocol["target_label_derangement_null"]
    bootstrap_spec = protocol["whole_target_bootstrap"]
    derangements = generate_derangements(
        int(null_spec["replicates"]),
        len(targets),
        int(null_spec["seed"]),
    )
    bootstrap = generate_bootstrap(
        int(bootstrap_spec["replicates"]),
        len(targets),
        int(bootstrap_spec["seed"]),
    )
    validate_plans(derangements, bootstrap, len(targets))

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent)
    )
    try:
        np.save(temporary / "target-derangements.npy", derangements)
        np.save(temporary / "whole-target-bootstrap.npy", bootstrap)
        plan_files = {
            name: {
                "sha256": sha256_file(temporary / name),
                "size_bytes": (temporary / name).stat().st_size,
            }
            for name in (
                "target-derangements.npy",
                "whole-target-bootstrap.npy",
            )
        }
        manifest = {
            "schema_version": 1,
            "experiment": protocol["experiment"],
            "status": "plans_frozen_before_source_expression_value_access",
            "target_count": len(targets),
            "target_order": list(targets),
            "numpy_version": np.__version__,
            "plans": plan_files,
            "code_and_protocol": code_bindings(),
            "source_bindings": authorization["source_bindings"],
            "safety": {
                "source_expression_arrays_opened": False,
                "held_truth_opened": False,
                "existing_outcomes_opened": False,
                "cell_eval2": False,
                "counts_generated": False,
                "submission_paths_opened": False,
            },
        }
        (temporary / "freeze-manifest.json").write_bytes(json_bytes(manifest))
        os.replace(temporary, output)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return output


def verify(output: Path = DEFAULT_OUTPUT) -> None:
    output = safe_output(output)
    require(output.is_dir() and not output.is_symlink(), "freeze output missing")
    manifest = read_json(output / "freeze-manifest.json")
    require(
        manifest["status"] == "plans_frozen_before_source_expression_value_access",
        "freeze status changed",
    )
    require(manifest["code_and_protocol"] == code_bindings(), "frozen code drift")
    authorization = read_json(AUTHORIZATION_PATH)
    validate_bindings(authorization)
    targets = validate_target_axis(read_json(PROTOCOL_PATH))
    require(tuple(manifest["target_order"]) == targets, "frozen target order changed")

    for name, record in manifest["plans"].items():
        path = output / name
        require(path.stat().st_size == record["size_bytes"], f"{name} size changed")
        require(sha256_file(path) == record["sha256"], f"{name} hash changed")

    protocol = read_json(PROTOCOL_PATH)
    expected_derangements = generate_derangements(
        int(protocol["target_label_derangement_null"]["replicates"]),
        len(targets),
        int(protocol["target_label_derangement_null"]["seed"]),
    )
    expected_bootstrap = generate_bootstrap(
        int(protocol["whole_target_bootstrap"]["replicates"]),
        len(targets),
        int(protocol["whole_target_bootstrap"]["seed"]),
    )
    actual_derangements = np.load(
        output / "target-derangements.npy",
        allow_pickle=False,
    )
    actual_bootstrap = np.load(
        output / "whole-target-bootstrap.npy",
        allow_pickle=False,
    )
    require(
        np.array_equal(actual_derangements, expected_derangements),
        "derangement content is not deterministic",
    )
    require(
        np.array_equal(actual_bootstrap, expected_bootstrap),
        "bootstrap content is not deterministic",
    )
    validate_plans(actual_derangements, actual_bootstrap, len(targets))


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("build", "verify"))
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    if args.command == "build":
        output = build(args.output)
        print(output)
    else:
        verify(args.output)
        print("preregistered freeze verified")


if __name__ == "__main__":
    main()

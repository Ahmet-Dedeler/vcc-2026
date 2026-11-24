#!/usr/bin/env python3
"""Freeze Experiment 122 null assignments before effect values are loaded."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL = HERE / "preregistered_protocol.json"
INVENTORY = HERE / "SOURCE_INVENTORY.json"
AUTHORIZATION = HERE / "SOURCE_ONLY_AUTHORIZATION.json"
PREREGISTRATION = HERE / "PREREGISTRATION.sha256"
DEFAULT_OUTPUT = HERE / "work" / "preregistered-freeze-v1"

EXPECTED_PREREGISTRATION = {
    "preregistered_protocol.json": (
        "446950dbf9ffe5f139df60a353849f993ef66a35a352feb232dab7c576756236"
    ),
    "SOURCE_INVENTORY.json": (
        "cb407c3ced26c4055f076147826216bb46b80ac8111704913de93d23b614c82b"
    ),
    "SOURCE_ONLY_AUTHORIZATION.json": (
        "25bc4e404280ecf47dc64ec69e85c8770991e28f953c2e93613b61c4b61e060d"
    ),
}


class FreezeError(RuntimeError):
    """A preregistration or immutable-output invariant failed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FreezeError(message)


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"JSON root must be an object: {path}")
    return value


def json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def read_panel() -> tuple[str, ...]:
    path = ROOT / "data/challenge/pert_counts.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ["target_gene"], "panel CSV header changed")
        panel = tuple(str(row["target_gene"]) for row in reader)
    require(len(panel) == 300, "panel must contain exactly 300 targets")
    require(len(panel) == len(set(panel)), "panel contains duplicate targets")
    return panel


def verify_preregistration() -> None:
    observed_lines = PREREGISTRATION.read_text(encoding="utf-8").splitlines()
    expected_lines = [
        f"{digest}  {name}" for name, digest in EXPECTED_PREREGISTRATION.items()
    ]
    require(observed_lines == expected_lines, "PREREGISTRATION.sha256 changed")
    for name, expected in EXPECTED_PREREGISTRATION.items():
        require(sha256_file(HERE / name) == expected, f"{name} hash changed")
    authorization = read_json(AUTHORIZATION)
    require(authorization["source_only_authorized"] is True, "source-only not authorized")
    require(
        authorization["protocol_sha256"]
        == EXPECTED_PREREGISTRATION["preregistered_protocol.json"],
        "authorization protocol binding changed",
    )
    require(
        authorization["source_inventory_sha256"]
        == EXPECTED_PREREGISTRATION["SOURCE_INVENTORY.json"],
        "authorization inventory binding changed",
    )


def derangements(draws: int, size: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    identity = np.arange(size, dtype=np.int16)
    output = np.empty((draws, size), dtype=np.int16)
    for draw in range(draws):
        for _ in range(10000):
            candidate = rng.permutation(size).astype(np.int16)
            if np.all(candidate != identity):
                output[draw] = candidate
                break
        else:
            raise FreezeError(f"could not draw derangement {draw}")
    return output


def validate_bank(bank: np.ndarray, draws: int, size: int) -> None:
    require(bank.shape == (draws, size), "derangement bank shape changed")
    require(np.issubdtype(bank.dtype, np.integer), "derangements are not integer")
    identity = np.arange(size)
    require(np.all(bank != identity[None, :]), "derangement has a fixed point")
    require(
        all(np.array_equal(np.sort(row), identity) for row in bank),
        "derangement is not bijective",
    )


def safe_output(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    require(resolved != HERE and resolved.is_relative_to(HERE), "unsafe output path")
    return resolved


def build(output: Path = DEFAULT_OUTPUT) -> Path:
    verify_preregistration()
    output = safe_output(output)
    require(not output.exists(), f"refusing to overwrite {output}")
    protocol = read_json(PROTOCOL)
    panel = read_panel()
    nulls = protocol["complete_selection_nulls"]
    draws = int(nulls["draws_per_family"])
    target = derangements(draws, len(panel), int(nulls["target_label"]["seed"]))
    response = derangements(
        draws,
        len(panel),
        int(nulls["response_label"]["seed"]),
    )
    validate_bank(target, draws, len(panel))
    validate_bank(response, draws, len(panel))

    output.parent.mkdir(parents=True, exist_ok=True)
    stage = output.parent / f".{output.name}.staging"
    require(not stage.exists(), f"staging directory exists: {stage}")
    stage.mkdir()
    try:
        with (stage / "target-label-derangements.npy").open("xb") as handle:
            np.save(handle, target, allow_pickle=False)
        with (stage / "response-label-derangements.npy").open("xb") as handle:
            np.save(handle, response, allow_pickle=False)
        (stage / "panel-targets.txt").write_text(
            "".join(f"{target_name}\n" for target_name in panel),
            encoding="utf-8",
        )
        files = {}
        for path in sorted(stage.iterdir()):
            files[path.name] = {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
        manifest = {
            "schema_version": 1,
            "experiment": "122-panel-response-atlas",
            "phase": "phase_1_source_only",
            "created_before_effect_matrix_loading": True,
            "draws_per_null_family": draws,
            "panel_targets": len(panel),
            "preregistration": EXPECTED_PREREGISTRATION,
            "files": files,
        }
        (stage / "freeze-manifest.json").write_bytes(json_bytes(manifest))
        os.sync()
        stage.rename(output)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return output


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    verify_preregistration()
    output = safe_output(output)
    manifest = read_json(output / "freeze-manifest.json")
    require(manifest["preregistration"] == EXPECTED_PREREGISTRATION, "freeze drift")
    for name, record in manifest["files"].items():
        path = output / name
        require(path.is_file() and not path.is_symlink(), f"missing frozen file: {name}")
        require(path.stat().st_size == record["size_bytes"], f"size changed: {name}")
        require(sha256_file(path) == record["sha256"], f"hash changed: {name}")
    target = np.load(output / "target-label-derangements.npy", allow_pickle=False)
    response = np.load(output / "response-label-derangements.npy", allow_pickle=False)
    validate_bank(target, int(manifest["draws_per_null_family"]), 300)
    validate_bank(response, int(manifest["draws_per_null_family"]), 300)
    require(
        tuple((output / "panel-targets.txt").read_text(encoding="utf-8").splitlines())
        == read_panel(),
        "frozen panel axis changed",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    result = verify(args.output) if args.verify else {"output": str(build(args.output))}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

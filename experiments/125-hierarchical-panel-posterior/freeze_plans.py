#!/usr/bin/env python3
"""Freeze Experiment 125 folds, nulls, and bootstrap plans before effect loading."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
from pathlib import Path
from typing import Any, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL = HERE / "preregistered_protocol.json"
AUTHORIZATION = HERE / "SOURCE_ONLY_AUTHORIZATION.json"
PREREGISTRATION = HERE / "PREREGISTRATION.sha256"
DEFAULT_OUTPUT = HERE / "work" / "preregistered-freeze-v1"


class FreezeError(RuntimeError):
    """A preregistration or plan invariant failed."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FreezeError(message)


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"JSON root must be an object: {path}")
    return value


def read_panel() -> tuple[str, ...]:
    source = ROOT / "data/challenge/pert_counts.csv"
    with source.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == ["target_gene"], "panel header changed")
        panel = tuple(str(row["target_gene"]) for row in reader)
    require(len(panel) == 300, "panel size changed")
    require(len(panel) == len(set(panel)), "panel contains duplicates")
    return panel


def stable_target_folds(panel: Sequence[str], n_folds: int) -> np.ndarray:
    order = sorted(
        range(len(panel)),
        key=lambda index: (
            hashlib.sha256(
                f"experiment-125-target-fold|{panel[index]}".encode()
            ).hexdigest(),
            index,
        ),
    )
    folds = np.empty(len(panel), dtype=np.int8)
    for rank, index in enumerate(order):
        folds[index] = rank % n_folds
    return folds


def random_derangement(rng: np.random.Generator, size: int) -> np.ndarray:
    identity = np.arange(size, dtype=np.int16)
    for _ in range(10000):
        candidate = rng.permutation(size).astype(np.int16)
        if np.all(candidate != identity):
            return candidate
    raise FreezeError(f"failed to draw derangement of size {size}")


def target_derangements(draws: int, size: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return np.stack([random_derangement(rng, size) for _ in range(draws)])


def response_derangements(
    draws: int,
    expression_bins: np.ndarray,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    size = len(expression_bins)
    output = np.tile(np.arange(size, dtype=np.int16), (draws, 1))
    for draw in range(draws):
        for bin_index in sorted(set(map(int, expression_bins)) - {-1}):
            members = np.flatnonzero(expression_bins == bin_index)
            require(members.size >= 2, "expression bin too small")
            local = random_derangement(rng, len(members))
            output[draw, members] = members[local]
    return output


def validate_target_bank(bank: np.ndarray, draws: int, size: int) -> None:
    identity = np.arange(size)
    require(bank.shape == (draws, size), "target bank shape changed")
    require(np.issubdtype(bank.dtype, np.integer), "target bank is not integer")
    require(np.all(bank != identity[None, :]), "target derangement fixed point")
    require(
        all(np.array_equal(np.sort(row), identity) for row in bank),
        "target derangement is not bijective",
    )


def validate_response_bank(
    bank: np.ndarray,
    bins: np.ndarray,
    draws: int,
) -> None:
    size = len(bins)
    identity = np.arange(size)
    require(bank.shape == (draws, size), "response bank shape changed")
    require(np.issubdtype(bank.dtype, np.integer), "response bank is not integer")
    require(
        all(np.array_equal(np.sort(row), identity) for row in bank),
        "response derangement is not bijective",
    )
    for row in bank:
        absent = bins < 0
        require(np.array_equal(row[absent], identity[absent]), "absent response moved")
        for bin_index in sorted(set(map(int, bins)) - {-1}):
            members = np.flatnonzero(bins == bin_index)
            require(np.all(row[members] != members), "response fixed point")
            require(
                set(map(int, row[members])) == set(map(int, members)),
                "response crossed expression bin",
            )


def validate_bootstrap(bank: np.ndarray, draws: int, size: int) -> None:
    require(bank.shape == (draws, size), "bootstrap shape changed")
    require(np.issubdtype(bank.dtype, np.integer), "bootstrap bank is not integer")
    require(np.all((bank >= 0) & (bank < size)), "bootstrap index out of range")


def verify_preregistration() -> dict[str, Any]:
    lines = PREREGISTRATION.read_text(encoding="utf-8").splitlines()
    expected: dict[str, str] = {}
    for line in lines:
        digest, name = line.split("  ", 1)
        expected[name] = digest
    require(
        set(expected) == {"preregistered_protocol.json", "SOURCE_ONLY_AUTHORIZATION.json"},
        "preregistration file list changed",
    )
    for name, digest in expected.items():
        require(sha256_file(HERE / name) == digest, f"{name} hash changed")
    authorization = read_json(AUTHORIZATION)
    require(
        authorization["protocol_sha256"]
        == expected["preregistered_protocol.json"],
        "authorization protocol binding changed",
    )
    for relative, digest in authorization["allowed_input_files"].items():
        path = (ROOT / relative).resolve()
        require(path.is_relative_to(ROOT), f"input escapes repository: {relative}")
        require(path.is_file() and not path.is_symlink(), f"missing input: {relative}")
        require(sha256_file(path) == digest, f"input hash changed: {relative}")
    return authorization


def safe_output(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    require(resolved != HERE and resolved.is_relative_to(HERE), "unsafe output path")
    return resolved


def write_npy(path: Path, array: np.ndarray) -> None:
    with path.open("xb") as handle:
        np.save(handle, np.ascontiguousarray(array), allow_pickle=False)


def build(output: Path = DEFAULT_OUTPUT) -> Path:
    authorization = verify_preregistration()
    protocol = read_json(PROTOCOL)
    output = safe_output(output)
    require(not output.exists(), f"refusing to overwrite {output}")
    panel = read_panel()
    exp122_panel = tuple(
        (
            ROOT
            / protocol["inputs"]["panel_targets"]
        ).read_text(encoding="utf-8").splitlines()
    )
    require(panel == exp122_panel, "Experiment 122 panel axis drift")
    bins = np.load(
        ROOT / protocol["inputs"]["response_expression_bins"],
        allow_pickle=False,
    )
    require(bins.shape == (300,), "response expression-bin shape changed")
    n_folds = int(protocol["nested_validation"]["target_folds"])
    folds = stable_target_folds(panel, n_folds)
    draws = int(protocol["nulls"]["draws_each"])
    target_bank = target_derangements(
        draws,
        len(panel),
        int(protocol["nulls"]["target_label"]["seed"]),
    )
    response_bank = response_derangements(
        draws,
        bins,
        int(protocol["nulls"]["response_label_within_expression_bin"]["seed"]),
    )
    bootstrap_draws = int(protocol["bootstrap"]["replicates"])
    bootstrap_rng = np.random.default_rng(int(protocol["bootstrap"]["seed"]))
    bootstrap = bootstrap_rng.integers(
        0,
        len(panel),
        size=(bootstrap_draws, len(panel)),
        dtype=np.int16,
    )
    validate_target_bank(target_bank, draws, len(panel))
    validate_response_bank(response_bank, bins, draws)
    validate_bootstrap(bootstrap, bootstrap_draws, len(panel))
    require(
        np.array_equal(np.bincount(folds, minlength=n_folds), np.full(n_folds, 60)),
        "target folds are not balanced",
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    stage = output.parent / f".{output.name}.staging"
    require(not stage.exists(), f"staging directory exists: {stage}")
    stage.mkdir()
    try:
        write_npy(stage / "target-folds.npy", folds)
        write_npy(stage / "target-label-derangements.npy", target_bank)
        write_npy(stage / "response-label-derangements.npy", response_bank)
        write_npy(stage / "response-expression-bins.npy", bins)
        write_npy(stage / "target-cluster-bootstrap.npy", bootstrap)
        (stage / "panel-targets.txt").write_text(
            "".join(f"{target}\n" for target in panel),
            encoding="utf-8",
        )
        files = {
            path.name: {
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for path in sorted(stage.iterdir())
        }
        manifest = {
            "schema_version": 1,
            "experiment": "125-hierarchical-panel-posterior",
            "phase": "source_only_nested_evaluation",
            "created_before_source_effect_loading": True,
            "protocol_sha256": sha256_file(PROTOCOL),
            "authorization_sha256": sha256_file(AUTHORIZATION),
            "allowed_input_count": len(authorization["allowed_input_files"]),
            "target_folds": n_folds,
            "null_draws_each": draws,
            "bootstrap_replicates": bootstrap_draws,
            "files": files,
        }
        (stage / "freeze-manifest.json").write_bytes(json_bytes(manifest))
        os.sync()
        stage.rename(output)
        for path in output.iterdir():
            path.chmod(0o444)
        output.chmod(0o555)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return output


def verify(output: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    verify_preregistration()
    output = safe_output(output)
    manifest = read_json(output / "freeze-manifest.json")
    for name, record in manifest["files"].items():
        path = output / name
        require(path.is_file() and not path.is_symlink(), f"missing frozen file: {name}")
        require(path.stat().st_size == record["size_bytes"], f"size changed: {name}")
        require(sha256_file(path) == record["sha256"], f"hash changed: {name}")
    folds = np.load(output / "target-folds.npy", allow_pickle=False)
    target = np.load(output / "target-label-derangements.npy", allow_pickle=False)
    response = np.load(output / "response-label-derangements.npy", allow_pickle=False)
    bins = np.load(output / "response-expression-bins.npy", allow_pickle=False)
    bootstrap = np.load(output / "target-cluster-bootstrap.npy", allow_pickle=False)
    validate_target_bank(target, int(manifest["null_draws_each"]), 300)
    validate_response_bank(response, bins, int(manifest["null_draws_each"]))
    validate_bootstrap(bootstrap, int(manifest["bootstrap_replicates"]), 300)
    require(folds.shape == (300,), "target fold shape changed")
    require(
        tuple((output / "panel-targets.txt").read_text(encoding="utf-8").splitlines())
        == read_panel(),
        "panel axis changed",
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

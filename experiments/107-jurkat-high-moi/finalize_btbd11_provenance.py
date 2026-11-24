#!/usr/bin/env python3
"""Finalize Experiment 107 role manifests after the BTBD11 validation repair."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
VALIDATOR_PATH = HERE / "strict_nadig_validation.py"
CHAIN_PATH = HERE / "BTBD11_VALIDATION_PROVENANCE_AMENDMENT.json"
ARTIFACT_ROOT = HERE / "work" / "strict-validation"
ARCHIVE_ROOT = HERE / "work" / "strict-validation-pre-final-provenance"
RECEIPT_PATH = HERE / "work" / "btbd11-final-provenance-reseal-receipt.json"


def _load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


STRICT = _load_module("exp107_final_provenance_validator", VALIDATOR_PATH)


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: dict[str, Any], mode: int = 0o444) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(temporary_name, mode)
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def verify_chain() -> dict[str, Any]:
    chain = STRICT.load_json(CHAIN_PATH)
    if chain.get("schema_version") != 1:
        raise RuntimeError("unsupported BTBD11 provenance amendment schema")
    if chain.get("experiment") != "107-jurkat-high-moi":
        raise RuntimeError("BTBD11 provenance amendment belongs to another experiment")
    for name, binding in chain["bindings"].items():
        path = ROOT / binding["path"]
        observed = STRICT.file_binding(path, root=ROOT)
        if observed != binding:
            raise RuntimeError(f"BTBD11 provenance binding changed: {name}")
    return chain


def _verify_npz_and_arrays(manifest: dict[str, Any]) -> None:
    npz_path = ROOT / manifest["npz"]["path"]
    if STRICT.file_binding(npz_path, root=ROOT) != manifest["npz"]:
        raise RuntimeError(f"{manifest['role']}: NPZ binding changed")
    with np.load(npz_path, allow_pickle=False) as archive:
        observed = {
            name: STRICT.array_binding(np.asarray(archive[name]))
            for name in archive.files
        }
    if observed != manifest["arrays"]:
        raise RuntimeError(f"{manifest['role']}: logical arrays changed")


def finalize() -> dict[str, Any]:
    if RECEIPT_PATH.exists() or ARCHIVE_ROOT.exists():
        raise RuntimeError("final provenance archive or receipt already exists")
    verify_chain()
    contract = STRICT.validate_contract()
    manifests = {
        role: STRICT.expected_role_bundle_paths(role)[1]
        for role in STRICT.REQUIRED_ARTIFACT_ROLES
    }
    if not all(path.is_file() and not path.is_symlink() for path in manifests.values()):
        raise RuntimeError("one or more required role manifests are missing or symlinked")

    pending = dict(manifests)
    completed: set[str] = set()
    records: list[dict[str, Any]] = []
    while pending:
        ready = sorted(
            role
            for role, path in pending.items()
            if set(STRICT.load_json(path)["parent_roles"]) <= completed
        )
        if not ready:
            raise RuntimeError("role-manifest dependency graph is cyclic")
        for role in ready:
            manifest_path = pending.pop(role)
            payload = STRICT.load_json(manifest_path)
            _verify_npz_and_arrays(payload)
            expected_paths = tuple(STRICT.expected_transitive_paths(role, contract))
            expected_relative = {
                path.resolve().relative_to(ROOT.resolve()).as_posix()
                for path in expected_paths
            }
            old_relative = {item["path"] for item in payload["transitive_inputs"]}
            allowed_new = {
                path.resolve().relative_to(ROOT.resolve()).as_posix()
                for path in STRICT.REPAIR_PROVENANCE_PATHS
            }
            if not expected_relative - old_relative <= allowed_new:
                raise RuntimeError(f"{role}: unexpected transitive-input additions")
            for item in payload["transitive_inputs"]:
                path = ROOT / item["path"]
                observed = STRICT.file_binding(path, root=ROOT)
                if observed == item:
                    continue
                if path.resolve() in {
                    STRICT.ROLE_BUILDER_PATH.resolve(),
                    VALIDATOR_PATH.resolve(),
                }:
                    continue
                if (
                    path.resolve().is_relative_to(ARTIFACT_ROOT.resolve())
                    and path.suffix == ".json"
                ):
                    continue
                raise RuntimeError(
                    f"{role}: unrelated transitive input changed: {item['path']}"
                )
            archive_path = ARCHIVE_ROOT / manifest_path.relative_to(ARTIFACT_ROOT)
            archive_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(manifest_path, archive_path)
            archive_path.chmod(0o444)
            replacement = json.loads(json.dumps(payload))
            replacement["transitive_inputs"] = [
                STRICT.file_binding(path, root=ROOT)
                for path in sorted(
                    {path.resolve() for path in expected_paths},
                    key=str,
                )
            ]
            old_binding = STRICT.file_binding(manifest_path, root=ROOT)
            atomic_json(manifest_path, replacement)
            records.append(
                {
                    "role": role,
                    "npz": payload["npz"],
                    "old_manifest": old_binding,
                    "new_manifest": STRICT.file_binding(manifest_path, root=ROOT),
                }
            )
            completed.add(role)

    binding = STRICT.bind_artifact_manifest(
        STRICT.ARTIFACT_MANIFEST_PATH,
        contract=contract,
    )
    receipt = {
        "schema_version": 1,
        "experiment": "107-jurkat-high-moi",
        "status": "btbd11_validation_and_provenance_repair_complete",
        "provenance_amendment": STRICT.file_binding(CHAIN_PATH, root=ROOT),
        "role_count": len(records),
        "roles": records,
        "artifact_binding_sha256": binding["binding_sha256"],
        "held_outcomes_opened": False,
        "npz_bytes_changed": False,
    }
    atomic_json(RECEIPT_PATH, receipt)
    return receipt


if __name__ == "__main__":
    print(json.dumps(finalize(), indent=2, sort_keys=True))

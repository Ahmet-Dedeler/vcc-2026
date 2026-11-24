#!/usr/bin/env python3
"""Reseal unchanged Experiment 107 artifacts after the BTBD11 code repair.

The failed build had already produced 23 valid source-only role NPZs.  The
repair changes only the previously unreachable integration branch for a score
coordinate absent from the incumbent operator.  This tool verifies every
existing artifact and non-builder input, preserves the original manifests,
and updates only the strict-artifact-builder file binding.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import tempfile
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ARTIFACT_ROOT = HERE / "work" / "strict-validation"
ARCHIVE_ROOT = HERE / "work" / "strict-validation-pre-btbd11-repair"
RECEIPT_PATH = HERE / "work" / "btbd11-reseal-receipt.json"
PARENT_RECEIPT_PATH = (
    HERE / "work" / "btbd11-parent-binding-reconciliation-receipt.json"
)
AMENDMENT_PATH = HERE / "BTBD11_OPERATOR_SUPPORT_REPAIR.json"
BUILDER_PATH = HERE / "strict_artifact_builder.py"
VALIDATOR_PATH = HERE / "strict_nadig_validation.py"
CONTRACT_PATH = HERE / "strict_validation_contract.json"


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def array_binding(values: np.ndarray) -> dict[str, Any]:
    array = np.ascontiguousarray(np.asarray(values))
    digest = hashlib.sha256()
    digest.update(b"exp107-array-v1\0")
    digest.update(array.dtype.str.encode("ascii"))
    digest.update(b"\0")
    digest.update(json.dumps(list(array.shape), separators=(",", ":")).encode("ascii"))
    digest.update(b"\0")
    digest.update(memoryview(array).cast("B"))
    return {
        "dtype": array.dtype.str,
        "shape": list(array.shape),
        "c_order": bool(array.flags.c_contiguous),
        "sha256_array": digest.hexdigest(),
    }


def file_binding(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    if not resolved.is_relative_to(ROOT.resolve()):
        raise RuntimeError(f"path escapes repository: {path}")
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f"path is missing or symlinked: {path}")
    return {
        "path": resolved.relative_to(ROOT.resolve()).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def atomic_json(path: Path, value: dict[str, Any], mode: int = 0o444) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def load_amendment() -> dict[str, Any]:
    amendment = json.loads(AMENDMENT_PATH.read_text(encoding="utf-8"))
    if amendment.get("schema_version") != 1:
        raise RuntimeError("unsupported repair amendment schema")
    if amendment.get("experiment") != "107-jurkat-high-moi":
        raise RuntimeError("repair amendment belongs to another experiment")
    bindings = amendment["bindings"]
    current = {
        "strict_validation_contract": file_binding(CONTRACT_PATH),
        "strict_nadig_validation": file_binding(VALIDATOR_PATH),
        "repaired_strict_artifact_builder": file_binding(BUILDER_PATH),
    }
    for name, observed in current.items():
        expected = bindings[name]
        if observed["path"] != expected["path"] or observed["sha256"] != expected["sha256"]:
            raise RuntimeError(f"repair binding changed: {name}")
        if "bytes" in expected and observed["bytes"] != expected["bytes"]:
            raise RuntimeError(f"repair byte count changed: {name}")
    return amendment


def verify_and_reseal(*, apply: bool) -> dict[str, Any]:
    amendment = load_amendment()
    old_builder = amendment["bindings"]["old_strict_artifact_builder"]
    new_builder = amendment["bindings"]["repaired_strict_artifact_builder"]
    manifests = sorted(ARTIFACT_ROOT.rglob("*.json"))
    expected_count = amendment["existing_artifact_reuse"]["expected_existing_role_count"]
    if len(manifests) != expected_count:
        raise RuntimeError(
            f"expected {expected_count} existing role manifests, found {len(manifests)}"
        )

    records: list[dict[str, Any]] = []
    for manifest_path in manifests:
        if manifest_path.is_symlink() or not manifest_path.is_file():
            raise RuntimeError(f"unsafe role manifest: {manifest_path}")
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
        role = payload.get("role")
        if not isinstance(role, str) or not role:
            raise RuntimeError(f"missing role in {manifest_path}")
        npz_binding = payload.get("npz")
        if not isinstance(npz_binding, dict):
            raise RuntimeError(f"{role}: missing NPZ binding")
        npz_path = ROOT / str(npz_binding["path"])
        if file_binding(npz_path) != npz_binding:
            raise RuntimeError(f"{role}: NPZ binding changed")
        with np.load(npz_path, allow_pickle=False) as archive:
            observed_arrays = {
                name: array_binding(np.asarray(archive[name])) for name in archive.files
            }
        if observed_arrays != payload.get("arrays"):
            raise RuntimeError(f"{role}: logical array binding changed")

        builder_records = [
            item
            for item in payload.get("transitive_inputs", [])
            if item.get("path") == old_builder["path"]
        ]
        if len(builder_records) != 1 or builder_records[0] != old_builder:
            raise RuntimeError(f"{role}: old builder binding is not exact")
        for item in payload["transitive_inputs"]:
            if item["path"] == old_builder["path"]:
                continue
            if file_binding(ROOT / item["path"]) != item:
                raise RuntimeError(f"{role}: non-builder transitive input changed")

        old_manifest_binding = file_binding(manifest_path)
        replacement = json.loads(json.dumps(payload))
        replacement["transitive_inputs"] = [
            new_builder if item["path"] == old_builder["path"] else item
            for item in replacement["transitive_inputs"]
        ]
        records.append(
            {
                "role": role,
                "npz": npz_binding,
                "old_manifest": old_manifest_binding,
                "new_manifest_sha256": hashlib.sha256(
                    (
                        json.dumps(replacement, indent=2, sort_keys=True) + "\n"
                    ).encode("utf-8")
                ).hexdigest(),
                "path": manifest_path.relative_to(ROOT).as_posix(),
                "replacement": replacement,
            }
        )

    if apply:
        if RECEIPT_PATH.exists() or ARCHIVE_ROOT.exists():
            raise RuntimeError("repair archive or receipt already exists")
        for record in records:
            manifest_path = ROOT / record["path"]
            archive_path = ARCHIVE_ROOT / manifest_path.relative_to(ARTIFACT_ROOT)
            archive_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(manifest_path, archive_path)
            archive_path.chmod(0o444)
            atomic_json(manifest_path, record["replacement"])
        receipt = {
            "schema_version": 1,
            "experiment": "107-jurkat-high-moi",
            "status": "existing_source_only_role_manifests_resealed",
            "repair_amendment": file_binding(AMENDMENT_PATH),
            "old_builder": old_builder,
            "new_builder": new_builder,
            "role_count": len(records),
            "roles": [
                {key: value for key, value in record.items() if key != "replacement"}
                for record in records
            ],
            "held_outcomes_opened": False,
            "npz_bytes_changed": False,
        }
        atomic_json(RECEIPT_PATH, receipt)
    return {
        "status": "applied" if apply else "validated",
        "role_count": len(records),
        "roles": [record["role"] for record in records],
    }


def reconcile_parent_bindings(*, apply: bool) -> dict[str, Any]:
    """Refresh only role-manifest bindings that changed during the first reseal."""

    load_amendment()
    manifests = {
        json.loads(path.read_text(encoding="utf-8"))["role"]: path
        for path in ARTIFACT_ROOT.rglob("*.json")
    }
    if len(manifests) != 23:
        raise RuntimeError(f"expected 23 live role manifests, found {len(manifests)}")
    pending = dict(manifests)
    completed: set[str] = set()
    records: list[dict[str, Any]] = []
    while pending:
        ready = sorted(
            role
            for role, path in pending.items()
            if set(json.loads(path.read_text(encoding="utf-8"))["parent_roles"])
            <= (completed | (set(manifests) - set(pending)))
        )
        if not ready:
            raise RuntimeError("role-manifest dependency graph is cyclic")
        for role in ready:
            manifest_path = pending.pop(role)
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            parent_roles = set(payload["parent_roles"])
            replacements: list[dict[str, Any]] = []
            changed_paths: list[str] = []
            for item in payload["transitive_inputs"]:
                source_path = ROOT / item["path"]
                observed = file_binding(source_path)
                if observed == item:
                    replacements.append(item)
                    continue
                relative = source_path.resolve().relative_to(ARTIFACT_ROOT.resolve())
                if relative.suffix != ".json":
                    raise RuntimeError(
                        f"{role}: non-parent transitive input changed: {item['path']}"
                    )
                parent_role = None
                for candidate in parent_roles:
                    candidate_paths = [
                        path.resolve()
                        for path in ARTIFACT_ROOT.rglob("*.json")
                        if json.loads(path.read_text(encoding="utf-8")).get("role")
                        == candidate
                    ]
                    if source_path.resolve() in candidate_paths:
                        parent_role = candidate
                        break
                if parent_role is None:
                    raise RuntimeError(
                        f"{role}: changed JSON is not a declared parent: {item['path']}"
                    )
                archived = ARCHIVE_ROOT / relative
                if not archived.is_file() or archived.is_symlink():
                    raise RuntimeError(f"{role}: archived parent manifest is missing")
                if (
                    archived.stat().st_size != item["bytes"]
                    or sha256_file(archived) != item["sha256"]
                ):
                    raise RuntimeError(
                        f"{role}: old parent binding does not match preserved manifest"
                    )
                replacements.append(observed)
                changed_paths.append(item["path"])
            replacement = json.loads(json.dumps(payload))
            replacement["transitive_inputs"] = replacements
            old_binding = file_binding(manifest_path)
            if apply and changed_paths:
                atomic_json(manifest_path, replacement)
            records.append(
                {
                    "role": role,
                    "old_manifest": old_binding,
                    "new_manifest": file_binding(manifest_path)
                    if apply
                    else {
                        "path": old_binding["path"],
                        "bytes": len(
                            (
                                json.dumps(replacement, indent=2, sort_keys=True) + "\n"
                            ).encode("utf-8")
                        ),
                        "sha256": hashlib.sha256(
                            (
                                json.dumps(replacement, indent=2, sort_keys=True) + "\n"
                            ).encode("utf-8")
                        ).hexdigest(),
                    },
                    "changed_parent_paths": changed_paths,
                }
            )
            completed.add(role)
    if apply:
        if PARENT_RECEIPT_PATH.exists():
            raise RuntimeError("parent-binding reconciliation receipt already exists")
        atomic_json(
            PARENT_RECEIPT_PATH,
            {
                "schema_version": 1,
                "experiment": "107-jurkat-high-moi",
                "status": "parent_role_manifest_bindings_reconciled",
                "repair_amendment": file_binding(AMENDMENT_PATH),
                "role_count": len(records),
                "roles": records,
                "held_outcomes_opened": False,
                "npz_bytes_changed": False,
            },
        )
    return {
        "status": "applied" if apply else "validated",
        "role_count": len(records),
        "changed_roles": [
            record["role"] for record in records if record["changed_parent_paths"]
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--apply",
        action="store_true",
        help="preserve old manifests and atomically reseal live manifests",
    )
    parser.add_argument(
        "--reconcile-parent-bindings",
        action="store_true",
        help="refresh only parent JSON bindings after the initial reseal",
    )
    args = parser.parse_args()
    if args.reconcile_parent_bindings:
        result = reconcile_parent_bindings(apply=args.apply)
    else:
        result = verify_and_reseal(apply=args.apply)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

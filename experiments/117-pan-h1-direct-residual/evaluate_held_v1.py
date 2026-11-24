#!/usr/bin/env python3
"""One-shot held evaluation for Experiment 117.

Importing this module and running ``authenticate`` are truth-safe. The fixed
truth files are not statted, hashed, or opened until ``run`` has atomically
created and fsynced the permanent outcome-access receipt.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
import zipfile
from collections.abc import Iterable, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXPERIMENT = "117-pan-h1-direct-residual"
PROTOCOL_PATH = HERE / "HELD_EVALUATION_PROTOCOL.v1.json"
AUTHORIZATION_PATH = HERE / "HELD_EVALUATION_AUTHORIZATION.v1.json"
ENTRYPOINT_PATH = Path(__file__).resolve()
ACCESS_RECEIPT_PATH = HERE / "work/held-evaluation-v1-access-receipt.json"
RESULT_DIRECTORY = HERE / "work/held-evaluation-v1"
EXPECTED_PROTOCOL_SHA256 = (
    "32a4433cf66eb131a3dad10b2da215b48c209c2af1bda388cbf467f2e74cfd9c"
)
SOURCE_ROOT = HERE / "work/source-candidates"
SOURCE_TOP_PATH = SOURCE_ROOT / "source-candidate-inputs.json"
PRIMARY_ARM = "pan_novel_f100"
H1_DIAGNOSTIC_ARM = "pan_all_f100"
IDENTITY_ARM = "identity_f000"
H1_HELD_ARMS = (IDENTITY_ARM, H1_DIAGNOSTIC_ARM)
AUTHORIZING_HELD_ARMS = (IDENTITY_ARM, PRIMARY_ARM)
PRIMARY_FRACTION = 0.10
DIRECTIONAL_FRACTION = 0.20
NULL_REPLICATES = 99
BOOTSTRAP_REPLICATES = 20_000
BOOTSTRAP_NAMESPACE = "experiment-117-pan-direct-cluster-bootstrap-v1"
NADIG_ROLES = tuple(f"nadig.gallery_{index}" for index in range(1, 7))
CROSS_CONTEXT_ROLES = ("orion_hek293t", *NADIG_ROLES)
ALL_ROLES = ("h1", *CROSS_CONTEXT_ROLES)
RESULT_ARRAYS_NAME = "evaluation-arrays.npz"
RESULT_SUMMARY_NAME = "evaluation-summary.json"
RESULT_MANIFEST_NAME = "result-manifest.json"
TS = 50_000.0
REQUIRED_ARRAYS = {
    "targets",
    "genes",
    "score_indices",
    "incumbent_effects",
    "source_effects",
    "coordinate_supported",
    "target_supported",
    "target_novel_vs_vcc2025_h1",
    "normalized_residual",
    "target_effective_all",
    "target_effective_novel",
    "identity_f000",
    "pan_novel_f050",
    "pan_novel_f100",
    "pan_novel_f200",
    "pan_all_f050",
    "pan_all_f100",
    "pan_all_f200",
    "directional_rank_indices",
    "directional_rank_indices_all",
    "target_cluster_labels",
    "target_label_null_permutations",
}


class EvaluationError(RuntimeError):
    """A frozen authorization, input, metric, or output invariant failed."""


@dataclass(frozen=True)
class SourceBundle:
    role: str
    manifest_path: Path
    npz_path: Path
    manifest: dict[str, Any]
    arrays: dict[str, np.ndarray]


@dataclass(frozen=True)
class ReferenceDE:
    genes: np.ndarray
    targets: np.ndarray
    lfc: np.ndarray
    significant: tuple[np.ndarray, ...]
    gate_indices: np.ndarray


@dataclass(frozen=True)
class AccessLease:
    path: Path
    descriptor: int
    sha256: str


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise EvaluationError(message)


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON constant: {value}")


def load_json(path: Path) -> dict[str, Any]:
    try:
        with Path(path).open(encoding="utf-8") as handle:
            value = json.load(
                handle,
                object_pairs_hook=_reject_duplicate_keys,
                parse_constant=_reject_constant,
            )
    except (OSError, ValueError, TypeError) as error:
        raise EvaluationError(f"cannot read valid JSON: {path}") from error
    _require(isinstance(value, dict), f"JSON root must be an object: {path}")
    return value


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256(
        "".join(f"{value}\n" for value in values).encode("utf-8")
    ).hexdigest()


def sha256_array(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values))
    _require(not array.dtype.hasobject, "object arrays cannot be bound")
    digest = hashlib.sha256()
    digest.update(b"exp107-array-v1\0")
    digest.update(array.dtype.str.encode("ascii"))
    digest.update(b"\0")
    digest.update(json.dumps(list(array.shape), separators=(",", ":")).encode("ascii"))
    digest.update(b"\0")
    if array.nbytes:
        digest.update(memoryview(array).cast("B"))
    return digest.hexdigest()


def array_binding(values: np.ndarray) -> dict[str, Any]:
    array = np.ascontiguousarray(np.asarray(values))
    return {
        "dtype": array.dtype.str,
        "shape": list(array.shape),
        "c_order": True,
        "sha256_array": sha256_array(array),
    }


def sha256_json_payload(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _canonical_relative(relative: str, label: str) -> Path:
    _require(
        isinstance(relative, str)
        and relative
        and relative == Path(relative).as_posix()
        and not Path(relative).is_absolute()
        and ".." not in Path(relative).parts,
        f"{label} path is not canonical repository-relative POSIX",
    )
    return Path(relative)


def _reject_symlink_components(path: Path, label: str) -> None:
    resolved_root = ROOT.resolve()
    candidate = Path(os.path.abspath(path))
    _require(
        candidate.is_relative_to(resolved_root),
        f"{label} escapes repository",
    )
    relative = candidate.relative_to(resolved_root)
    current = resolved_root
    for part in relative.parts:
        current = current / part
        _require(not current.is_symlink(), f"{label} contains a symlink")
    _require(
        candidate.resolve().is_relative_to(resolved_root),
        f"{label} resolves outside repository",
    )


def _source_file(relative: str, label: str) -> Path:
    path = ROOT / _canonical_relative(relative, label)
    _reject_symlink_components(path, label)
    _require(path.is_file(), f"{label} is missing")
    return path.resolve()


def _verify_source_binding(binding: Mapping[str, Any], label: str) -> Path:
    _require(
        isinstance(binding, Mapping)
        and set(binding) == {"path", "bytes", "sha256"},
        f"{label} binding schema changed",
    )
    path = _source_file(str(binding["path"]), label)
    _require(path.stat().st_size == int(binding["bytes"]), f"{label} size changed")
    _require(sha256_file(path) == str(binding["sha256"]), f"{label} hash changed")
    return path


def _validate_truth_declaration(binding: Mapping[str, Any], label: str) -> None:
    """Validate declaration syntax only; never touch the declared truth path."""
    _require(
        isinstance(binding, Mapping)
        and set(binding) == {"path", "bytes", "sha256"},
        f"{label} truth binding schema changed",
    )
    _canonical_relative(str(binding["path"]), label)
    _require(int(binding["bytes"]) > 0, f"{label} truth size is invalid")
    digest = str(binding["sha256"])
    _require(
        len(digest) == 64 and all(char in "0123456789abcdef" for char in digest),
        f"{label} truth hash is invalid",
    )


def _verify_truth_binding_after_receipt(
    binding: Mapping[str, Any],
    label: str,
    lease: AccessLease,
) -> Path:
    _require(lease.path == ACCESS_RECEIPT_PATH.resolve(), "invalid access lease")
    _require(
        ACCESS_RECEIPT_PATH.is_file()
        and sha256_file(ACCESS_RECEIPT_PATH) == lease.sha256,
        "outcome-access receipt changed before truth verification",
    )
    return _verify_source_binding(binding, label)


def _decode_strings(values: np.ndarray) -> tuple[str, ...]:
    return tuple(
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in np.asarray(values).tolist()
    )


def _validate_array(
    archive: Mapping[str, np.ndarray],
    declared: Mapping[str, Any],
    name: str,
) -> np.ndarray:
    _require(name in archive and name in declared, f"source array missing: {name}")
    value = np.asarray(archive[name])
    _require(array_binding(value) == declared[name], f"source array drifted: {name}")
    return value


def _role_directory(role: str) -> Path:
    return SOURCE_ROOT / role.replace(".", "/")


def load_source_bundle(
    role: str,
    role_binding: Mapping[str, Any],
) -> SourceBundle:
    _require(role in ALL_ROLES, f"unexpected source role: {role}")
    manifest_path = _verify_source_binding(
        role_binding["manifest"], f"{role} source manifest"
    )
    npz_path = _verify_source_binding(role_binding["npz"], f"{role} source NPZ")
    expected_directory = _role_directory(role).resolve()
    _require(manifest_path.parent == expected_directory, f"{role} directory changed")
    _require(npz_path.parent == expected_directory, f"{role} NPZ directory changed")
    manifest = load_json(manifest_path)
    _require(
        manifest.get("experiment") == EXPERIMENT
        and manifest.get("role") == role
        and manifest.get("status")
        == "frozen_before_experiment_117_held_outcome_access",
        f"{role} source manifest identity changed",
    )
    _require(
        manifest.get("npz", {}).get("sha256") == role_binding["npz"]["sha256"]
        and manifest.get("npz", {}).get("bytes") == role_binding["npz"]["bytes"],
        f"{role} source manifest-to-NPZ binding changed",
    )
    _require(
        manifest.get("safety")
        == {
            "held_truth_opened": False,
            "cell_eval2_run": False,
            "counts_generated": False,
            "submission_built": False,
        },
        f"{role} source safety declaration changed",
    )
    with np.load(npz_path, allow_pickle=False) as archive:
        _require(set(archive.files) == REQUIRED_ARRAYS, f"{role} array inventory changed")
        arrays = {
            name: np.asarray(_validate_array(archive, manifest["arrays"], name)).copy()
            for name in sorted(REQUIRED_ARRAYS)
        }
    return SourceBundle(
        role=role,
        manifest_path=manifest_path,
        npz_path=npz_path,
        manifest=manifest,
        arrays=arrays,
    )


def _assert_source_invariants(bundle: SourceBundle, protocol: Mapping[str, Any]) -> None:
    arrays = bundle.arrays
    targets = _decode_strings(arrays["targets"])
    genes = _decode_strings(arrays["genes"])
    _require(len(targets) == len(set(targets)), f"{bundle.role} target duplicates")
    _require(len(genes) == len(set(genes)), f"{bundle.role} gene duplicates")
    base = np.asarray(arrays["incumbent_effects"])
    identity = np.asarray(arrays["identity_f000"])
    _require(
        base.shape == (len(targets), len(genes)),
        f"{bundle.role} incumbent shape changed",
    )
    _require(
        np.array_equal(base.view(np.uint8), identity.view(np.uint8)),
        f"{bundle.role} identity is not bit-exact",
    )
    support = np.asarray(arrays["coordinate_supported"], dtype=bool)
    effective_all = np.asarray(arrays["target_effective_all"], dtype=bool)
    effective_novel = np.asarray(arrays["target_effective_novel"], dtype=bool)
    _require(
        support.shape == base.shape
        and effective_all.shape == effective_novel.shape == (len(targets),),
        f"{bundle.role} support shapes changed",
    )
    _require(
        not np.any(effective_novel & ~effective_all),
        f"{bundle.role} novel mask escapes all-supported mask",
    )
    for population, effective in (
        ("novel", effective_novel),
        ("all", effective_all),
    ):
        allowed = support & effective[:, None]
        fallback = ~allowed
        for fraction in (50, 100, 200):
            candidate = np.asarray(arrays[f"pan_{population}_f{fraction:03d}"])
            _require(candidate.shape == base.shape, f"{bundle.role} arm shape changed")
            _require(
                np.array_equal(
                    np.ascontiguousarray(candidate[fallback]).view(np.uint8),
                    np.ascontiguousarray(base[fallback]).view(np.uint8),
                ),
                f"{bundle.role} frozen fallback changed",
            )
    mappings = np.asarray(arrays["target_label_null_permutations"])
    _require(
        mappings.shape == (NULL_REPLICATES, len(targets))
        and mappings.dtype == np.dtype("<i4"),
        f"{bundle.role} null assignments changed",
    )
    labels = np.asarray(arrays["target_cluster_labels"])
    _require(
        labels.shape == (len(targets),) and np.all(labels >= 0),
        f"{bundle.role} cluster labels changed",
    )
    required_arms = set(protocol["source_candidate_contract"]["required_arms"])
    _require(required_arms.issubset(arrays), f"{bundle.role} frozen arms are incomplete")


def _validate_nadig_mapping(
    role: str,
    bundle: SourceBundle,
    plan: Mapping[str, Any],
    declared: Mapping[str, Any],
) -> None:
    number = int(role.rsplit("_", 1)[1])
    gallery = plan["galleries"][number - 1]
    _require(gallery["gallery"] == number, f"{role} plan gallery changed")
    targets = _decode_strings(bundle.arrays["targets"])
    genes = np.asarray(_decode_strings(bundle.arrays["genes"]))
    score_indices = np.asarray(bundle.arrays["score_indices"], dtype=np.int64)
    score_genes = tuple(genes[score_indices].tolist())
    _require(tuple(map(str, gallery["targets"])) == targets, f"{role} targets changed")
    _require(
        tuple(map(str, gallery["score_genes"])) == score_genes,
        f"{role} score genes changed",
    )
    observed = {
        "targets_ordered_sha256": sha256_lines(targets),
        "score_genes_ordered_sha256": sha256_lines(score_genes),
        "target_indices_in_truth": array_binding(
            np.asarray(gallery["target_indices_in_nadig"], dtype=np.int64)
        ),
        "score_indices_in_truth": array_binding(
            np.asarray(gallery["score_indices_in_nadig"], dtype=np.int64)
        ),
        "score_indices_in_official": array_binding(
            np.asarray(gallery["score_indices_in_official"], dtype=np.int64)
        ),
    }
    _require(observed == dict(declared), f"{role} frozen role mapping changed")


def _validate_h5ad_mapping(
    role: str,
    bundle: SourceBundle,
    declared: Mapping[str, Any],
) -> None:
    targets = _decode_strings(bundle.arrays["targets"])
    genes = np.asarray(_decode_strings(bundle.arrays["genes"]))
    score_indices = np.asarray(bundle.arrays["score_indices"], dtype=np.int64)
    observed = {
        "target_mapping": "exact_symbol_in_frozen_source_target_order",
        "pds_gene_mapping": "exact_symbol_from_frozen_score_indices",
        "directional_gene_mapping": (
            "official_order_exact_symbol_intersection_with_truth_axis_no_alias_merge"
        ),
        "targets_ordered_sha256": sha256_lines(targets),
        "score_genes_ordered_sha256": sha256_lines(
            genes[score_indices].astype(str)
        ),
        "score_indices_in_official": array_binding(score_indices),
    }
    _require(observed == dict(declared), f"{role} frozen role mapping changed")


def held_arms_for_role(role: str) -> tuple[str, str]:
    if role == "h1":
        return H1_HELD_ARMS
    _require(role in CROSS_CONTEXT_ROLES, f"unknown held-evaluation role: {role}")
    return AUTHORIZING_HELD_ARMS


def _held_output_contract() -> dict[str, Any]:
    return {
        "h1_scored_and_persisted_arms": list(H1_HELD_ARMS),
        "authorizing_surface_scored_and_persisted_arms": list(
            AUTHORIZING_HELD_ARMS
        ),
        "authorizing_surfaces": list(CROSS_CONTEXT_ROLES),
        "forbidden_arms_must_never_be_read_scored_or_persisted": True,
        "nadig_score_columns_must_match_frozen_bundle_order_before_slicing": True,
    }


def _protocol_held_outcome_exposure() -> dict[str, Any]:
    return {
        "h1_scored_and_persisted_arms": list(H1_HELD_ARMS),
        "authorizing_surface_scored_and_persisted_arms": list(
            AUTHORIZING_HELD_ARMS
        ),
        "authorizing_surfaces": list(CROSS_CONTEXT_ROLES),
        "forbidden_held_outcomes": [
            "pan_novel_f050",
            "pan_novel_f200",
            "pan_all_f050",
            "pan_all_f200",
            "pan_all_f100_on_authorizing_surfaces",
            "pan_novel_f100_on_h1",
        ],
        "forbidden_arms_must_never_be_read_scored_or_persisted": True,
        "nadig_score_columns_must_exactly_equal_frozen_bundle_score_genes_in_order_before_truth_slicing": True,
    }


def _score_held_arms(
    bundle: SourceBundle,
    truth_score: np.ndarray,
) -> dict[str, np.ndarray]:
    score = np.asarray(bundle.arrays["score_indices"], dtype=np.int64)
    return {
        arm: pds_per_target(
            np.asarray(bundle.arrays[arm], dtype=np.float64)[:, score],
            truth_score,
        )
        for arm in held_arms_for_role(bundle.role)
    }


def _nadig_score_mapping(
    bundle: SourceBundle,
    gallery: Mapping[str, Any],
    truth_genes: Sequence[str],
) -> tuple[np.ndarray, np.ndarray]:
    score_truth = np.asarray(gallery["score_indices_in_nadig"], dtype=np.int64)
    score_official = np.asarray(bundle.arrays["score_indices"], dtype=np.int64)
    genes = np.asarray(tuple(map(str, truth_genes)))
    official = np.asarray(_decode_strings(bundle.arrays["genes"]))
    _require(
        score_truth.ndim == score_official.ndim == 1
        and len(score_truth) == len(score_official)
        and np.all((score_truth >= 0) & (score_truth < len(genes)))
        and np.all((score_official >= 0) & (score_official < len(official))),
        f"{bundle.role} Nadig score mapping is invalid",
    )
    _require(
        np.array_equal(genes[score_truth], official[score_official]),
        f"{bundle.role} Nadig score gene order changed",
    )
    return score_truth, score_official


def authenticate_pretruth() -> dict[str, Any]:
    """Authenticate all source-side state without touching any truth path."""
    _require(
        sha256_file(PROTOCOL_PATH) == EXPECTED_PROTOCOL_SHA256,
        "held-evaluation protocol changed",
    )
    protocol = load_json(PROTOCOL_PATH)
    authorization = load_json(AUTHORIZATION_PATH)
    _require(
        protocol.get("schema_version") == 1
        and protocol.get("freeze_revision") == 2
        and protocol.get("experiment") == EXPERIMENT
        and protocol.get("status")
        == "frozen_before_experiment_117_held_truth_access_after_independent_audit_repair",
        "held-evaluation protocol identity changed",
    )
    _require(
        authorization.get("schema_version") == 1
        and authorization.get("freeze_revision") == 2
        and authorization.get("experiment") == EXPERIMENT
        and authorization.get("status")
        == "held_evaluation_v1_repaired_and_refrozen_pending_independent_review",
        "held-evaluation authorization identity changed",
    )
    _require(
        protocol.get("held_outcome_exposure")
        == _protocol_held_outcome_exposure(),
        "held-outcome exposure protocol changed",
    )
    _require(
        authorization.get("held_output_contract") == _held_output_contract(),
        "held-output authorization contract changed",
    )
    _require(
        authorization.get("authorization")
        == {
            "source_candidate_read": True,
            "held_truth_access_after_receipt": True,
            "held_evaluation_after_receipt": True,
            "model_refit": False,
            "arm_or_fraction_selection": False,
            "cell_eval2": False,
            "count_generation": False,
            "challenge_build": False,
            "submission": False,
            "commit": False,
            "push": False,
        },
        "held-evaluation authorization map changed",
    )
    _require(
        authorization.get("execution_contract")
        == {
            "access_receipt": ACCESS_RECEIPT_PATH.relative_to(ROOT).as_posix(),
            "result_directory": RESULT_DIRECTORY.relative_to(ROOT).as_posix(),
            "receipt_must_be_exclusive_fsynced_and_read_only_before_truth_hash_or_load": True,
            "receipt_survives_failure_and_prevents_replay": True,
            "caller_supplied_paths": False,
            "result_overwrite": False,
            "result_manifest_immutable": True,
            "required_null_replicates": NULL_REPLICATES,
            "required_cluster_bootstrap_replicates": BOOTSTRAP_REPLICATES,
            "required_cross_context_roles": list(CROSS_CONTEXT_ROLES),
            "h1_diagnostic_cannot_promote": True,
        },
        "held-evaluation execution contract changed",
    )
    _require(
        authorization["protocol"]
        == {
            "path": PROTOCOL_PATH.relative_to(ROOT).as_posix(),
            "bytes": PROTOCOL_PATH.stat().st_size,
            "sha256": EXPECTED_PROTOCOL_SHA256,
        },
        "authorization does not bind the exact held-evaluation protocol",
    )
    evaluator_binding = authorization["evaluator"]
    _require(
        evaluator_binding["path"] == ENTRYPOINT_PATH.relative_to(ROOT).as_posix()
        and evaluator_binding["bytes"] == ENTRYPOINT_PATH.stat().st_size
        and evaluator_binding["sha256"] == sha256_file(ENTRYPOINT_PATH),
        "authorization does not bind this evaluator",
    )
    predecessor_paths = {
        name: _verify_source_binding(binding, f"predecessor {name}")
        for name, binding in authorization["predecessors"].items()
    }
    top_binding = authorization["source_candidates"]["top_manifest"]
    top_path = _verify_source_binding(top_binding, "source-candidate top manifest")
    _require(top_path == SOURCE_TOP_PATH.resolve(), "source-candidate root changed")
    top = load_json(top_path)
    _require(
        top.get("status") == "source_candidates_frozen_before_held_outcome_access"
        and top.get("held_truth_opened") is False,
        "source-candidate top manifest is not pretruth-frozen",
    )
    _require(set(top["roles"]) == set(ALL_ROLES), "top-manifest role inventory changed")
    _require(
        tuple(authorization["source_candidates"]["roles"]) == ALL_ROLES,
        "authorization role order changed",
    )
    _require(
        set(authorization["role_bindings"]) == set(ALL_ROLES)
        and set(authorization["role_mappings"]) == set(ALL_ROLES),
        "authorization role inventory changed",
    )
    _require(
        set(authorization["truth_bindings"])
        == {"h1", "orion_hek293t", "nadig_jurkat"},
        "authorization truth inventory changed",
    )
    plan = load_json(predecessor_paths["experiment_107_pretruth_plan"])
    _require(
        plan.get("status") == "pretruth_complete"
        and plan.get("nadig_outcome_opened") is False,
        "Experiment 107 pretruth plan changed",
    )
    bundles: dict[str, SourceBundle] = {}
    for role in ALL_ROLES:
        declared = authorization["role_bindings"][role]
        _require(
            top["roles"][role]
            == {
                "manifest": declared["manifest"]["path"],
                "sha256": declared["manifest"]["sha256"],
            },
            f"{role} top-manifest binding changed",
        )
        bundle = load_source_bundle(role, declared)
        _assert_source_invariants(bundle, protocol)
        if role.startswith("nadig."):
            _validate_nadig_mapping(
                role,
                bundle,
                plan,
                authorization["role_mappings"][role],
            )
        else:
            _validate_h5ad_mapping(
                role,
                bundle,
                authorization["role_mappings"][role],
            )
        bundles[role] = bundle
    for label, binding in authorization["truth_bindings"].items():
        if label == "nadig_jurkat":
            for part, item in binding["files"].items():
                _validate_truth_declaration(item, f"{label} {part}")
        else:
            _validate_truth_declaration(binding["harness"], f"{label} harness")
            _validate_truth_declaration(
                binding["reference_de"], f"{label} reference DE"
            )
    _require(
        tuple(authorization["truth_bindings"]["h1"]["scored_arms"])
        == H1_HELD_ARMS,
        "H1 held-arm exposure changed",
    )
    _require(
        tuple(authorization["truth_bindings"]["orion_hek293t"]["scored_arms"])
        == AUTHORIZING_HELD_ARMS
        and tuple(
            authorization["truth_bindings"]["nadig_jurkat"]["scored_arms"]
        )
        == AUTHORIZING_HELD_ARMS,
        "authorizing held-arm exposure changed",
    )
    _require(
        not ACCESS_RECEIPT_PATH.exists() and not ACCESS_RECEIPT_PATH.is_symlink(),
        "held-truth access was already claimed",
    )
    _require(
        not RESULT_DIRECTORY.exists() and not RESULT_DIRECTORY.is_symlink(),
        "held-evaluation result already exists",
    )
    return {
        "protocol": protocol,
        "authorization": authorization,
        "authorization_sha256": sha256_file(AUTHORIZATION_PATH),
        "predecessor_paths": predecessor_paths,
        "plan": plan,
        "top": top,
        "bundles": bundles,
    }


def _receipt_payload(authenticated: Mapping[str, Any]) -> dict[str, Any]:
    authorization = authenticated["authorization"]
    return {
        "schema_version": 1,
        "experiment": EXPERIMENT,
        "status": "held_truth_access_committed_before_any_truth_hash_or_load",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "authorization": {
            "path": AUTHORIZATION_PATH.relative_to(ROOT).as_posix(),
            "bytes": AUTHORIZATION_PATH.stat().st_size,
            "sha256": authenticated["authorization_sha256"],
        },
        "protocol": dict(authorization["protocol"]),
        "evaluator": dict(authorization["evaluator"]),
        "source_candidate_top_manifest": dict(
            authorization["source_candidates"]["top_manifest"]
        ),
        "truth_bindings": authorization["truth_bindings"],
        "truth_hash_or_values_opened_before_receipt": False,
        "cell_eval2_run": False,
        "counts_generated": False,
        "challenge_built": False,
        "submission_performed": False,
    }


@contextmanager
def exclusive_access_receipt(
    authenticated: Mapping[str, Any],
) -> Iterable[AccessLease]:
    """Create one permanent fsynced receipt and hold its exclusive lock."""
    path = ACCESS_RECEIPT_PATH
    path.parent.mkdir(parents=True, exist_ok=True)
    _require(not RESULT_DIRECTORY.exists(), "result directory already exists")
    payload = (
        json.dumps(
            _receipt_payload(authenticated),
            indent=2,
            sort_keys=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    try:
        descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o400)
    except FileExistsError as error:
        raise EvaluationError("held-truth access was already claimed") from error
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        written = 0
        while written < len(payload):
            written += os.write(descriptor, payload[written:])
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
        parent = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(parent)
        finally:
            os.close(parent)
        receipt_hash = hashlib.sha256(payload).hexdigest()
        _require(sha256_file(path) == receipt_hash, "access receipt write changed")
        yield AccessLease(path=path.resolve(), descriptor=descriptor, sha256=receipt_hash)
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)


def pds_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    predicted = np.asarray(prediction, dtype=np.float64)
    observed = np.asarray(truth, dtype=np.float64)
    _require(
        predicted.ndim == 2
        and predicted.shape == observed.shape
        and predicted.shape[0] >= 2,
        "PDS arrays do not align",
    )
    _require(
        np.all(np.isfinite(predicted)) and np.all(np.isfinite(observed)),
        "PDS arrays contain nonfinite values",
    )
    result = np.full(len(predicted), 0.5, dtype=np.float64)
    truth_norm = np.linalg.norm(observed, axis=1)
    truth_unit = np.divide(
        observed,
        truth_norm[:, None],
        out=np.zeros_like(observed),
        where=truth_norm[:, None] > 0,
    )
    rows = np.flatnonzero(np.linalg.norm(predicted, axis=1) > 0)
    if rows.size:
        unit = predicted[rows] / np.linalg.norm(predicted[rows], axis=1)[:, None]
        similarity = unit @ truth_unit.T
        own = similarity[np.arange(len(rows)), rows][:, None]
        rank = (similarity > own).sum(axis=1).astype(np.float64)
        rank += 0.5 * ((similarity == own).sum(axis=1) - 1)
        result[rows] = 1.0 - rank / (len(predicted) - 1)
    return result


def aggregate_cross_context(surface_values: Mapping[str, float]) -> dict[str, float]:
    _require(
        set(surface_values) == set(CROSS_CONTEXT_ROLES),
        "HEK293T and all six Nadig galleries are mandatory",
    )
    values = {role: float(surface_values[role]) for role in CROSS_CONTEXT_ROLES}
    _require(
        all(math.isfinite(value) for value in values.values()),
        "cross-context value is nonfinite",
    )
    jurkat = math.fsum(values[role] for role in NADIG_ROLES) / 6.0
    macro = math.fsum((values["orion_hek293t"], jurkat)) / 2.0
    return {
        "orion_hek293t": values["orion_hek293t"],
        "nadig_jurkat": jurkat,
        "cross_context": macro,
    }


def deterministic_seed(namespace: str, role: str, draw: int) -> int:
    payload = f"{namespace}\0{role}\0{draw}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def cluster_bootstrap(
    per_target_delta: Mapping[str, np.ndarray],
    cluster_labels: Mapping[str, np.ndarray],
    *,
    repetitions: int = BOOTSTRAP_REPLICATES,
) -> dict[str, Any]:
    _require(
        set(per_target_delta) == set(CROSS_CONTEXT_ROLES)
        and set(cluster_labels) == set(CROSS_CONTEXT_ROLES),
        "bootstrap requires HEK293T and all six Nadig galleries",
    )
    _require(0 < repetitions <= BOOTSTRAP_REPLICATES, "invalid bootstrap count")
    draws = np.empty(repetitions, dtype=np.float64)
    for draw in range(repetitions):
        surface: dict[str, float] = {}
        for role in CROSS_CONTEXT_ROLES:
            delta = np.asarray(per_target_delta[role], dtype=np.float64)
            labels = np.asarray(cluster_labels[role], dtype=np.int64)
            _require(
                delta.shape == labels.shape and np.all(labels >= 0),
                f"{role} bootstrap inputs do not align",
            )
            unique = np.unique(labels)
            rng = np.random.Generator(
                np.random.PCG64(deterministic_seed(BOOTSTRAP_NAMESPACE, role, draw))
            )
            sampled = rng.choice(unique, size=len(unique), replace=True)
            counts = np.bincount(sampled, minlength=int(unique[-1]) + 1)
            weights = counts[labels]
            denominator = int(weights.sum())
            _require(denominator > 0, f"{role} bootstrap draw is empty")
            surface[role] = float(np.dot(delta, weights) / denominator)
        draws[draw] = aggregate_cross_context(surface)["cross_context"]
    return {
        "repetitions": repetitions,
        "lower_95": float(np.quantile(draws, 0.025, method="linear")),
        "upper_95": float(np.quantile(draws, 0.975, method="linear")),
        "draws": draws,
    }


def normalized_residual_from_mapping(
    incumbent: np.ndarray,
    source: np.ndarray,
    support: np.ndarray,
    novel: np.ndarray,
    mapping: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    base = np.asarray(incumbent, dtype=np.float64)
    added = np.asarray(source, dtype=np.float64)
    covered = np.asarray(support, dtype=bool)
    novel_rows = np.asarray(novel, dtype=bool)
    assignment = np.asarray(mapping, dtype=np.int64)
    _require(
        base.shape == added.shape == covered.shape
        and assignment.shape == novel_rows.shape == (len(base),),
        "null reconstruction arrays do not align",
    )
    _require(
        np.all((assignment >= -1) & (assignment < len(base))),
        "null assignment is out of bounds",
    )
    residual = np.zeros_like(base)
    effective = np.zeros(len(base), dtype=bool)
    effective_support = np.zeros_like(covered)
    for destination in np.flatnonzero((assignment >= 0) & novel_rows):
        source_row = int(assignment[destination])
        _require(
            novel_rows[source_row],
            "frozen null assignment crosses the novel stratum",
        )
        mask = covered[source_row]
        if not np.any(mask):
            continue
        source_values = added[source_row, mask]
        incumbent_values = np.asarray(base[destination, mask], dtype=np.float64)
        if not np.all(np.isfinite(source_values)) or not np.all(
            np.isfinite(incumbent_values)
        ):
            continue
        source_norm = float(np.linalg.norm(source_values))
        incumbent_norm = float(np.linalg.norm(incumbent_values))
        if (
            not math.isfinite(source_norm)
            or source_norm <= 0
            or not math.isfinite(incumbent_norm)
            or incumbent_norm <= 0
        ):
            continue
        row = source_values * (incumbent_norm / source_norm) - incumbent_values
        if not np.all(np.isfinite(row)) or float(np.linalg.norm(row)) <= 0:
            continue
        residual[destination, mask] = row
        effective[destination] = True
        effective_support[destination, mask] = True
    return residual, effective, effective_support


def null_candidate_score(
    bundle: SourceBundle,
    mapping: np.ndarray,
) -> np.ndarray:
    arrays = bundle.arrays
    base = np.asarray(arrays["identity_f000"], dtype=np.float64)
    residual, effective, effective_support = normalized_residual_from_mapping(
        base,
        arrays["source_effects"],
        arrays["coordinate_supported"],
        arrays["target_novel_vs_vcc2025_h1"],
        mapping,
    )
    score = np.asarray(arrays["score_indices"], dtype=np.int64)
    candidate = base[:, score].copy()
    selected = effective_support[:, score] & effective[:, None]
    candidate[selected] += PRIMARY_FRACTION * residual[:, score][selected]
    return candidate


def directional_precision_delta(
    bundle: SourceBundle,
    reference_sign: Callable[[int, np.ndarray], np.ndarray],
    *,
    arm: str,
    effective_mask: str,
    ranking_name: str,
) -> dict[str, np.ndarray | float]:
    arrays = bundle.arrays
    base = np.asarray(arrays["identity_f000"], dtype=np.float64)
    candidate = np.asarray(arrays[arm], dtype=np.float64)
    effective = np.asarray(arrays[effective_mask], dtype=bool)
    ranking = np.asarray(arrays[ranking_name], dtype=np.int64)
    support = np.asarray(arrays["coordinate_supported"], dtype=bool)
    incumbent_values: list[float] = []
    candidate_values: list[float] = []
    rows: list[int] = []
    for row in np.flatnonzero(effective):
        count = int(support[row].sum())
        _require(count > 0, f"{bundle.role} effective row has no support")
        selected_count = max(1, int(math.ceil(DIRECTIONAL_FRACTION * count)))
        selected = ranking[row]
        selected = selected[selected >= 0][:selected_count]
        _require(
            len(selected) == selected_count
            and np.all(support[row, selected]),
            f"{bundle.role} frozen directional ranking is incomplete",
        )
        truth_sign = np.asarray(reference_sign(int(row), selected), dtype=np.int8)
        _require(
            truth_sign.shape == selected.shape
            and np.all(np.isin(truth_sign, (-1, 0, 1))),
            f"{bundle.role} reference signs are invalid",
        )
        candidate_correct = (truth_sign != 0) & (
            np.sign(candidate[row, selected]).astype(np.int8) == truth_sign
        )
        incumbent_correct = (truth_sign != 0) & (
            np.sign(base[row, selected]).astype(np.int8) == truth_sign
        )
        candidate_values.append(float(np.mean(candidate_correct)))
        incumbent_values.append(float(np.mean(incumbent_correct)))
        rows.append(int(row))
    _require(rows, f"{bundle.role} has no directional rows")
    incumbent_array = np.asarray(incumbent_values, dtype=np.float64)
    candidate_array = np.asarray(candidate_values, dtype=np.float64)
    delta = candidate_array - incumbent_array
    return {
        "rows": np.asarray(rows, dtype=np.int64),
        "incumbent": incumbent_array,
        "candidate": candidate_array,
        "delta": delta,
        "surface_delta": float(np.mean(delta)),
    }


def _sum_rows(matrix: Any, rows: np.ndarray, genes: int) -> np.ndarray:
    result = np.zeros(genes, dtype=np.float64)
    ordered = np.sort(np.asarray(rows, dtype=np.int64))
    _require(ordered.size > 0, "held outcome group has no cells")
    for start in range(0, len(ordered), 512):
        result += np.asarray(
            matrix[ordered[start : start + 512]].sum(axis=0),
            dtype=np.float64,
        ).ravel()
    _require(np.all(np.isfinite(result)) and result.sum() > 0, "invalid held sum")
    return result


def _load_h5ad_truth(path: Path, targets: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    import anndata as ad

    handle = ad.read_h5ad(path, backed="r")
    try:
        genes = np.asarray(handle.var_names).astype(str)
        labels = np.asarray(handle.obs["target_gene"].astype(str))
        expected = tuple(map(str, targets))
        _require(
            tuple(sorted(set(labels) - {"non-targeting"})) == expected,
            f"{path.name} target axis changed",
        )
        control_sum = _sum_rows(
            handle.X,
            np.flatnonzero(labels == "non-targeting"),
            len(genes),
        )
        control = np.log1p(TS * control_sum / control_sum.sum())
        effects = np.empty((len(expected), len(genes)), dtype=np.float64)
        for row, target in enumerate(expected):
            total = _sum_rows(handle.X, np.flatnonzero(labels == target), len(genes))
            effects[row] = np.log1p(TS * total / total.sum()) - control
    finally:
        handle.file.close()
    _require(np.all(np.isfinite(effects)), "held effects are nonfinite")
    return genes, effects


def _load_reference_de(path: Path) -> ReferenceDE:
    with np.load(path, allow_pickle=False) as archive:
        genes = np.asarray(archive["genes"]).astype(str)
        gate = np.asarray(archive["gate"], dtype=bool)
        targets = np.asarray(archive["perts"]).astype(str)
        lfc = np.asarray(archive["lfc"], dtype=np.float64)
        pointers = np.asarray(archive["sig_ptr"], dtype=np.int64)
        indices = np.asarray(archive["sig_idx"], dtype=np.int64)
    gate_indices = np.flatnonzero(gate)
    _require(
        lfc.shape == (len(targets), len(gate_indices))
        and pointers.shape == (len(targets) + 1,)
        and pointers[0] == 0
        and pointers[-1] == len(indices),
        "reference-DE arrays are malformed",
    )
    significant = tuple(
        indices[pointers[row] : pointers[row + 1]] for row in range(len(targets))
    )
    return ReferenceDE(
        genes=genes,
        targets=targets,
        lfc=lfc,
        significant=significant,
        gate_indices=gate_indices,
    )


def _h5ad_reference_signer(
    bundle: SourceBundle,
    held_genes: np.ndarray,
    reference: ReferenceDE,
) -> Callable[[int, np.ndarray], np.ndarray]:
    _require(
        np.array_equal(reference.genes, held_genes),
        f"{bundle.role} reference-DE axis differs from held axis",
    )
    target_lookup = {target: index for index, target in enumerate(reference.targets)}
    held_lookup = {gene: index for index, gene in enumerate(held_genes)}
    gate_lookup = {
        int(full): position for position, full in enumerate(reference.gate_indices)
    }
    targets = _decode_strings(bundle.arrays["targets"])
    genes = np.asarray(_decode_strings(bundle.arrays["genes"]))

    def signs(row: int, selected: np.ndarray) -> np.ndarray:
        target = targets[row]
        _require(target in target_lookup, f"reference DE lacks target {target}")
        ref_row = target_lookup[target]
        significant = set(map(int, reference.significant[ref_row]))
        result = np.zeros(len(selected), dtype=np.int8)
        for offset, official_index in enumerate(map(int, selected)):
            held_index = held_lookup.get(str(genes[official_index]))
            if held_index is None or held_index not in gate_lookup:
                continue
            gate_index = gate_lookup[held_index]
            if gate_index not in significant:
                continue
            value = float(reference.lfc[ref_row, gate_index])
            if math.isfinite(value) and value != 0:
                result[offset] = 1 if value > 0 else -1
        return result

    return signs


def _evaluate_h5ad_surface(
    bundle: SourceBundle,
    truth: Mapping[str, Any],
    lease: AccessLease,
) -> dict[str, Any]:
    harness = _verify_truth_binding_after_receipt(
        truth["harness"], f"{bundle.role} harness", lease
    )
    reference_path = _verify_truth_binding_after_receipt(
        truth["reference_de"], f"{bundle.role} reference DE", lease
    )
    targets = _decode_strings(bundle.arrays["targets"])
    held_genes, held_effects = _load_h5ad_truth(harness, targets)
    axes = truth["axes"]
    _require(
        len(held_genes) == axes["truth_gene_count"]
        and sha256_lines(held_genes) == axes["truth_gene_ordered_sha256"],
        f"{bundle.role} held gene axis changed",
    )
    _require(
        sha256_lines(targets) == axes["target_ordered_sha256"],
        f"{bundle.role} target axis changed",
    )
    official_genes = np.asarray(_decode_strings(bundle.arrays["genes"]))
    held_lookup = {gene: index for index, gene in enumerate(held_genes)}
    score = np.asarray(bundle.arrays["score_indices"], dtype=np.int64)
    score_genes = official_genes[score]
    _require(
        sha256_lines(score_genes.astype(str))
        == axes["score_gene_ordered_sha256"],
        f"{bundle.role} score axis changed",
    )
    _require(
        all(gene in held_lookup for gene in score_genes),
        f"{bundle.role} truth lacks a score gene",
    )
    truth_score = held_effects[
        :, np.asarray([held_lookup[gene] for gene in score_genes], dtype=np.int64)
    ]
    _require(
        tuple(truth["scored_arms"]) == held_arms_for_role(bundle.role),
        f"{bundle.role} held-arm exposure changed",
    )
    arm_scores = _score_held_arms(bundle, truth_score)
    reference = _load_reference_de(reference_path)
    directional = directional_precision_delta(
        bundle,
        _h5ad_reference_signer(bundle, held_genes, reference),
        arm=truth["directional_arm"],
        effective_mask=truth["effective_mask"],
        ranking_name=truth["directional_ranking"],
    )
    comparator = arm_scores[IDENTITY_ARM]
    deltas = {
        arm: arm_scores[arm] - comparator
        for arm in held_arms_for_role(bundle.role)
        if arm != IDENTITY_ARM
    }
    return {
        "pds": arm_scores,
        "pds_delta": deltas,
        "directional": directional,
        "held_genes": held_genes,
        "truth_score": truth_score,
    }


def _read_axis_after_receipt(path: Path) -> tuple[str, ...]:
    values = tuple(path.read_text(encoding="utf-8").splitlines())
    _require(values and len(values) == len(set(values)), f"{path.name} axis invalid")
    return values


def _evaluate_nadig(
    bundles: Mapping[str, SourceBundle],
    truth: Mapping[str, Any],
    plan: Mapping[str, Any],
    lease: AccessLease,
) -> dict[str, dict[str, Any]]:
    _require(
        tuple(truth["scored_arms"]) == AUTHORIZING_HELD_ARMS,
        "Nadig held-arm exposure changed",
    )
    paths = {
        name: _verify_truth_binding_after_receipt(
            binding, f"Nadig {name}", lease
        )
        for name, binding in truth["files"].items()
    }
    genes = _read_axis_after_receipt(paths["genes"])
    targets = _read_axis_after_receipt(paths["targets"])
    axes = truth["axes"]
    _require(
        len(genes) == axes["gene_count"]
        and sha256_lines(genes) == axes["gene_ordered_sha256"],
        "Nadig gene axis changed",
    )
    _require(
        len(targets) == axes["target_count"]
        and sha256_lines(targets) == axes["target_ordered_sha256"],
        "Nadig target axis changed",
    )
    score_mappings: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for role in NADIG_ROLES:
        gallery = plan["galleries"][int(role.rsplit("_", 1)[1]) - 1]
        score_mappings[role] = _nadig_score_mapping(
            bundles[role],
            gallery,
            genes,
        )
    values = np.load(paths["delta"], mmap_mode="r", allow_pickle=False)
    _require(
        values.shape == (len(targets), len(genes))
        and np.issubdtype(values.dtype, np.floating),
        "Nadig truth matrix changed",
    )
    gene_lookup = {gene: index for index, gene in enumerate(genes)}
    output: dict[str, dict[str, Any]] = {}
    for role in NADIG_ROLES:
        bundle = bundles[role]
        gallery = plan["galleries"][int(role.rsplit("_", 1)[1]) - 1]
        rows = np.asarray(gallery["target_indices_in_nadig"], dtype=np.int64)
        score_truth, score_official = score_mappings[role]
        _require(
            tuple(targets[index] for index in rows)
            == _decode_strings(bundle.arrays["targets"]),
            f"{role} target mapping changed",
        )
        truth_score = np.asarray(values[np.ix_(rows, score_truth)], dtype=np.float64)
        arm_scores = _score_held_arms(bundle, truth_score)
        official_genes = np.asarray(_decode_strings(bundle.arrays["genes"]))

        def signs(
            row: int,
            selected: np.ndarray,
            *,
            truth_rows: np.ndarray = rows,
            canonical: np.ndarray = official_genes,
        ) -> np.ndarray:
            result = np.zeros(len(selected), dtype=np.int8)
            for offset, official_index in enumerate(map(int, selected)):
                truth_column = gene_lookup.get(str(canonical[official_index]))
                if truth_column is None:
                    continue
                value = float(values[int(truth_rows[row]), truth_column])
                if math.isfinite(value) and value != 0:
                    result[offset] = 1 if value > 0 else -1
            return result

        directional = directional_precision_delta(
            bundle,
            signs,
            arm=truth["directional_arm"],
            effective_mask=truth["effective_mask"],
            ranking_name=truth["directional_ranking"],
        )
        identity = arm_scores[IDENTITY_ARM]
        output[role] = {
            "pds": arm_scores,
            "pds_delta": {
                arm: arm_scores[arm] - identity
                for arm in AUTHORIZING_HELD_ARMS
                if arm != IDENTITY_ARM
            },
            "directional": directional,
            "truth_rows": rows,
            "truth_score": truth_score,
        }
    return output


def _evaluate_nulls(
    bundles: Mapping[str, SourceBundle],
    truth_score: Mapping[str, np.ndarray],
    observed_delta: Mapping[str, np.ndarray],
) -> np.ndarray:
    gains = np.empty(NULL_REPLICATES, dtype=np.float64)
    for draw in range(NULL_REPLICATES):
        surface: dict[str, float] = {}
        for role in CROSS_CONTEXT_ROLES:
            bundle = bundles[role]
            mapping = np.asarray(
                bundle.arrays["target_label_null_permutations"][draw],
                dtype=np.int64,
            )
            candidate = null_candidate_score(bundle, mapping)
            null_pds = pds_per_target(candidate, truth_score[role])
            identity_pds = pds_per_target(
                np.asarray(bundle.arrays["identity_f000"], dtype=np.float64)[
                    :, np.asarray(bundle.arrays["score_indices"], dtype=np.int64)
                ],
                truth_score[role],
            )
            surface[role] = float(np.mean(null_pds - identity_pds))
        gains[draw] = aggregate_cross_context(surface)["cross_context"]
    _require(
        observed_delta.keys() == truth_score.keys(),
        "observed/null role sets differ",
    )
    return gains


def promotion_report(
    pds_surface: Mapping[str, float],
    directional_surface: Mapping[str, float],
    null_gains: np.ndarray,
    bootstrap: Mapping[str, Any],
    invariants_pass: bool,
) -> dict[str, Any]:
    pds_aggregate = aggregate_cross_context(pds_surface)
    directional_aggregate = aggregate_cross_context(directional_surface)
    null = np.asarray(null_gains, dtype=np.float64)
    _require(
        null.shape == (NULL_REPLICATES,) and np.all(np.isfinite(null)),
        "promotion requires exactly 99 finite null gains",
    )
    draws = np.asarray(bootstrap["draws"], dtype=np.float64)
    _require(
        bootstrap["repetitions"] == BOOTSTRAP_REPLICATES
        and draws.shape == (BOOTSTRAP_REPLICATES,)
        and bootstrap["lower_95"]
        == float(np.quantile(draws, 0.025, method="linear")),
        "promotion bootstrap is incomplete",
    )
    checks = {
        "hek293t_raw_pds_nonnegative": pds_surface["orion_hek293t"] >= 0.0,
        "hek293t_directional_precision_positive": (
            directional_surface["orion_hek293t"] > 0.0
        ),
        "all_six_nadig_raw_pds_nonnegative": all(
            pds_surface[role] >= 0.0 for role in NADIG_ROLES
        ),
        "nadig_macro_directional_precision_positive": (
            directional_aggregate["nadig_jurkat"] > 0.0
        ),
        "cross_context_macro_raw_pds_at_least_0_003": (
            pds_aggregate["cross_context"] >= 0.003
        ),
        "bootstrap_lower_95_positive": bootstrap["lower_95"] > 0.0,
        "observed_gain_exceeds_all_99_nulls": (
            pds_aggregate["cross_context"] > float(np.max(null))
        ),
        "all_source_and_fallback_invariants": bool(invariants_pass),
    }
    return {
        "passes": bool(all(checks.values())),
        "checks": checks,
        "raw_pds": pds_aggregate,
        "directional_precision": directional_aggregate,
        "null_maximum": float(np.max(null)),
        "bootstrap_lower_95": float(bootstrap["lower_95"]),
        "h1_used_for_promotion": False,
    }


def _json_safe(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float):
        _require(math.isfinite(value), "result JSON contains nonfinite value")
    return value


def _write_json_exclusive(path: Path, payload: Mapping[str, Any]) -> None:
    encoded = (
        json.dumps(_json_safe(payload), indent=2, sort_keys=True, allow_nan=False)
        + "\n"
    ).encode("utf-8")
    descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o400)
    try:
        written = 0
        while written < len(encoded):
            written += os.write(descriptor, encoded[written:])
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
    finally:
        os.close(descriptor)


def _write_npz_exclusive(path: Path, arrays: Mapping[str, np.ndarray]) -> None:
    _require(not path.exists(), f"refusing to overwrite {path}")
    with zipfile.ZipFile(
        path,
        mode="x",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=6,
        allowZip64=True,
    ) as archive:
        for name in sorted(arrays):
            value = np.ascontiguousarray(np.asarray(arrays[name]))
            _require(not value.dtype.hasobject, f"object result array: {name}")
            info = zipfile.ZipInfo(f"{name}.npy", (1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = 0o400 << 16
            with archive.open(info, "w", force_zip64=True) as handle:
                np.lib.format.write_array(handle, value, allow_pickle=False)
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    os.chmod(path, 0o444)


def _output_binding(actual: Path, final: Path) -> dict[str, Any]:
    return {
        "path": final.relative_to(ROOT).as_posix(),
        "bytes": actual.stat().st_size,
        "sha256": sha256_file(actual),
    }


def write_immutable_results(
    summary: Mapping[str, Any],
    arrays: Mapping[str, np.ndarray],
    inputs: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    destination = RESULT_DIRECTORY.resolve()
    _require(
        destination.parent == (HERE / "work").resolve(),
        "result directory changed",
    )
    _require(
        not destination.exists() and not destination.is_symlink(),
        "immutable result directory already exists",
    )
    try:
        destination.mkdir(mode=0o700)
    except FileExistsError as error:
        raise EvaluationError("immutable result directory already exists") from error
    try:
        arrays_path = destination / RESULT_ARRAYS_NAME
        summary_path = destination / RESULT_SUMMARY_NAME
        manifest_path = destination / RESULT_MANIFEST_NAME
        _write_npz_exclusive(arrays_path, arrays)
        _write_json_exclusive(summary_path, summary)
        manifest = {
            "schema_version": 1,
            "experiment": EXPERIMENT,
            "status": "immutable_held_evaluation_v1_complete",
            "outputs": {
                "arrays": _output_binding(
                    arrays_path, destination / RESULT_ARRAYS_NAME
                ),
                "summary": _output_binding(
                    summary_path, destination / RESULT_SUMMARY_NAME
                ),
            },
            "inputs": list(inputs),
        }
        manifest["binding_sha256"] = sha256_json_payload(manifest)
        _write_json_exclusive(manifest_path, manifest)
        directory = os.open(destination, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
        os.chmod(destination, 0o555)
        parent = os.open(destination.parent, os.O_RDONLY)
        try:
            os.fsync(parent)
        finally:
            os.close(parent)
        return manifest
    except BaseException:
        if destination.exists():
            for path in destination.iterdir():
                try:
                    path.chmod(0o400)
                except OSError:
                    pass
            try:
                destination.chmod(0o500)
            except OSError:
                pass
        raise


def _file_binding(path: Path) -> dict[str, Any]:
    resolved = Path(path).resolve()
    _require(resolved.is_relative_to(ROOT.resolve()), "input escapes repository")
    return {
        "path": resolved.relative_to(ROOT.resolve()).as_posix(),
        "bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def add_role_result_arrays(
    output: dict[str, np.ndarray],
    role: str,
    result: Mapping[str, Any],
) -> None:
    allowed = held_arms_for_role(role)
    pds = result["pds"]
    _require(
        isinstance(pds, Mapping) and tuple(pds) == allowed,
        f"{role} result exposes a forbidden held arm",
    )
    safe = role.replace(".", "_")
    for arm in allowed:
        output[f"{safe}_{arm}_pds"] = np.asarray(pds[arm], dtype=np.float64)
    directional = result["directional"]
    output[f"{safe}_directional_rows"] = np.asarray(
        directional["rows"], dtype=np.int64
    )
    output[f"{safe}_incumbent_directional_precision"] = np.asarray(
        directional["incumbent"], dtype=np.float64
    )
    output[f"{safe}_candidate_directional_precision"] = np.asarray(
        directional["candidate"], dtype=np.float64
    )
    output[f"{safe}_directional_precision_delta"] = np.asarray(
        directional["delta"], dtype=np.float64
    )


def run_real_evaluation() -> dict[str, Any]:
    """Create the one-shot receipt, then evaluate every mandatory frozen role."""
    authenticated = authenticate_pretruth()
    with exclusive_access_receipt(authenticated) as lease:
        auth = authenticated["authorization"]
        h1 = _evaluate_h5ad_surface(
            authenticated["bundles"]["h1"],
            auth["truth_bindings"]["h1"],
            lease,
        )
        hek = _evaluate_h5ad_surface(
            authenticated["bundles"]["orion_hek293t"],
            auth["truth_bindings"]["orion_hek293t"],
            lease,
        )
        nadig = _evaluate_nadig(
            authenticated["bundles"],
            auth["truth_bindings"]["nadig_jurkat"],
            authenticated["plan"],
            lease,
        )

        pds_surface = {
            "orion_hek293t": float(np.mean(hek["pds_delta"][PRIMARY_ARM])),
            **{
                role: float(np.mean(nadig[role]["pds_delta"][PRIMARY_ARM]))
                for role in NADIG_ROLES
            },
        }
        directional_surface = {
            "orion_hek293t": float(hek["directional"]["surface_delta"]),
            **{
                role: float(nadig[role]["directional"]["surface_delta"])
                for role in NADIG_ROLES
            },
        }
        per_target_delta = {
            "orion_hek293t": np.asarray(hek["pds_delta"][PRIMARY_ARM]),
            **{
                role: np.asarray(nadig[role]["pds_delta"][PRIMARY_ARM])
                for role in NADIG_ROLES
            },
        }
        clusters = {
            role: np.asarray(
                authenticated["bundles"][role].arrays["target_cluster_labels"],
                dtype=np.int32,
            )
            for role in CROSS_CONTEXT_ROLES
        }
        bootstrap = cluster_bootstrap(per_target_delta, clusters)

        truth_score = {
            "orion_hek293t": np.asarray(hek["truth_score"], dtype=np.float64),
            **{
                role: np.asarray(nadig[role]["truth_score"], dtype=np.float64)
                for role in NADIG_ROLES
            },
        }
        null_gains = _evaluate_nulls(
            authenticated["bundles"],
            truth_score,
            per_target_delta,
        )
        promotion = promotion_report(
            pds_surface,
            directional_surface,
            null_gains,
            bootstrap,
            invariants_pass=True,
        )

        summary = {
            "schema_version": 1,
            "experiment": EXPERIMENT,
            "status": "held_evaluation_v1_complete",
            "created_at_utc": datetime.now(UTC).isoformat(),
            "primary_comparison": "pan_novel_f100_minus_identity_f000",
            "h1_diagnostic": {
                "promotion_authorizing": False,
                "arm": H1_DIAGNOSTIC_ARM,
                "raw_pds_delta": float(
                    np.mean(h1["pds_delta"][H1_DIAGNOSTIC_ARM])
                ),
                "directional_precision_delta": float(
                    h1["directional"]["surface_delta"]
                ),
            },
            "surface_raw_pds_delta": pds_surface,
            "surface_directional_precision_delta": directional_surface,
            "cross_context_raw_pds": aggregate_cross_context(pds_surface),
            "cross_context_directional_precision": aggregate_cross_context(
                directional_surface
            ),
            "target_label_null": {
                "replicates": NULL_REPLICATES,
                "maximum": float(np.max(null_gains)),
            },
            "target_cluster_bootstrap": {
                "replicates": BOOTSTRAP_REPLICATES,
                "lower_95": bootstrap["lower_95"],
                "upper_95": bootstrap["upper_95"],
            },
            "promotion": promotion,
            "access_receipt": _file_binding(ACCESS_RECEIPT_PATH),
            "cell_eval2_run": False,
            "counts_generated": False,
            "challenge_built": False,
            "submission_performed": False,
        }
        arrays: dict[str, np.ndarray] = {
            "target_label_null_cross_context_gains": null_gains,
            "target_cluster_bootstrap_cross_context_draws": bootstrap["draws"],
        }
        for role, result in {
            "h1": h1,
            "orion_hek293t": hek,
            **nadig,
        }.items():
            add_role_result_arrays(arrays, role, result)

        input_paths = [
            AUTHORIZATION_PATH,
            PROTOCOL_PATH,
            ENTRYPOINT_PATH,
            ACCESS_RECEIPT_PATH,
            SOURCE_TOP_PATH,
            *authenticated["predecessor_paths"].values(),
        ]
        for bundle in authenticated["bundles"].values():
            input_paths.extend((bundle.manifest_path, bundle.npz_path))
        for truth_group in auth["truth_bindings"].values():
            if "files" in truth_group:
                input_paths.extend(
                    ROOT / binding["path"]
                    for binding in truth_group["files"].values()
                )
            else:
                input_paths.extend(
                    ROOT / truth_group[name]["path"]
                    for name in ("harness", "reference_de")
                )
        unique: list[Path] = []
        seen: set[Path] = set()
        for path in map(Path, input_paths):
            resolved = path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                unique.append(resolved)
        input_bindings = [_file_binding(path) for path in unique]
        for role, bundle in authenticated["bundles"].items():
            declared = auth["role_bindings"][role]
            _require(
                _file_binding(bundle.manifest_path) == declared["manifest"]
                and _file_binding(bundle.npz_path) == declared["npz"],
                f"{role} source candidates changed during evaluation",
            )
        manifest = write_immutable_results(summary, arrays, input_bindings)
        return {
            "status": summary["status"],
            "passes": promotion["passes"],
            "result_directory": RESULT_DIRECTORY.relative_to(ROOT).as_posix(),
            "result_binding_sha256": manifest["binding_sha256"],
        }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    commands.add_parser(
        "authenticate",
        help="verify source-only authorization without touching truth paths",
    )
    run = commands.add_parser("run", help="perform the one authorized held evaluation")
    run.add_argument(
        "--confirm-open-held-truth",
        action="store_true",
        help="required acknowledgement that the fixed held outcomes will be opened",
    )
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.command == "authenticate":
        authenticated = authenticate_pretruth()
        print(
            json.dumps(
                {
                    "experiment": EXPERIMENT,
                    "status": authenticated["authorization"]["status"],
                    "roles": list(authenticated["bundles"]),
                    "truth_paths_touched": False,
                    "access_receipt_absent": True,
                    "result_directory_absent": True,
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    _require(
        args.confirm_open_held_truth,
        "run requires --confirm-open-held-truth",
    )
    print(json.dumps(run_real_evaluation(), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

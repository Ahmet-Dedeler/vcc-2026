#!/usr/bin/env python3
"""Strict reused-development-benchmark validation for Experiment 107.

Nadig Jurkat was already used elsewhere in this repository and is therefore
not an independent confirmation set.  This module still enforces a hard
experiment-local access boundary: ``pretruth`` hashes but never parses the
Nadig outcome, and ``outcome-check`` validates a fully bound plan and writes a
non-replayable access receipt before NumPy opens ``delta.npy``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from pathlib import PurePosixPath
from typing import Any, Mapping, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CONTRACT_PATH = HERE / "strict_validation_contract.json"
SOURCE_DIRECTORY = HERE / "work" / "frozen-source-model"
ARTIFACT_MANIFEST_PATH = HERE / "work" / "strict-validation-inputs.json"
PRETRUTH_PLAN_PATH = HERE / "work" / "strict-validation-plan.json"
OUTCOME_ACCESS_RECEIPT_PATH = HERE / "work" / "nadig-outcome-access-receipt.json"
NADIG_DIRECTORY = ROOT / "data" / "donors" / "nadig_jurkat"
CHALLENGE_DIRECTORY = ROOT / "data" / "challenge"
DONOR_ROOT = ROOT / "data" / "donors"

EXPECTED_CONTRACT_SHA256 = (
    "8a6b874dddcfb2d4f172dacb6ed608c321b6039b82df15789364d727e1e572ee"
)
SOURCE_MODEL_KEYS = frozenset(
    {
        "exported_effects",
        "coefficients_binary_units",
        "projected_coefficients",
        "nuisance_response_slopes",
        "nuisance_target_slopes",
        "response_basis",
        "nullspace_basis",
        "evidence_support",
        "zero_variance",
        "nuisance_aliased",
        "fitted",
        "projection_diagonal",
        "transferable",
        "target_scales",
        "singular_values",
        "retained_singular_directions",
        "target_symbols",
        "response_symbols",
        "response_vcc_indices",
        "response_exportable",
        "selected_lambda",
        "lambda_grid",
        "held_cell_counts",
        "fold_mean_gains",
        "mean_held_cell_gains",
        "fold_selected_effects",
        "fold_zero_variance",
        "fold_nuisance_aliased",
        "fold_projection_diagonal",
        "fold_transferable",
        "fold_support_training_channels",
        "fold_support_distinct_called_guides",
        "fold_support_guides_with_three_untreated_cells",
        "fold_support_untreated_carrier_cells",
        "fold_support_training_channel_count",
        "fold_support_supported",
        "final_support_training_channels",
        "final_support_distinct_called_guides",
        "final_support_guides_with_three_untreated_cells",
        "final_support_untreated_carrier_cells",
        "final_support_training_channel_count",
        "final_support_supported",
        "fold_unmodeled_indptr",
        "fold_unmodeled_target_burden",
        "fold_stability_pairs",
        "fold_effect_cosine",
        "fold_effect_sign_agreement",
        "projection_nullspace_annihilation_error",
        "projection_idempotence_error",
        "projection_gauge_shift_invariance_error",
    }
)
DONOR_ORDER = ("k562_gwps", "orion_hct116", "orion_hek293t")
AUTHORIZING_HELD_SURFACES = ("h1", "orion_hek293t")
REQUIRED_ARTIFACT_ROLES = (
    *(f"donor.{name}" for name in DONOR_ORDER),
    *(f"incumbent.nadig.gallery_{number}" for number in range(1, 7)),
    *(f"operator.nadig.gallery_{number}" for number in range(1, 7)),
    *(f"integration.nadig.gallery_{number}" for number in range(1, 7)),
    "incumbent.h1",
    "incumbent.orion_hek293t",
    "operator.h1",
    "operator.orion_hek293t",
    "integration.h1",
    "integration.orion_hek293t",
)
PERMUTATION_REPLICATES = 99
BOOTSTRAP_REPLICATES = 20_000
FRACTIONS = (0.0, 0.05, 0.1, 0.2)
ROLE_BUNDLE_SCHEMA_VERSION = 1
ROLE_BUNDLE_STATUS = "frozen_before_experiment_107_nadig_outcome_access"
OFFICIAL_GENE_COUNT = 18_533
ROLE_BUNDLE_ROOT = HERE / "work" / "strict-validation"
ROLE_BUILDER_PATH = HERE / "strict_artifact_builder.py"
PARENT_PROTOCOL_PATH = HERE / "preregistered_protocol.json"
REPAIR_PROVENANCE_PATHS = (
    HERE / "BTBD11_OPERATOR_SUPPORT_REPAIR.json",
    HERE / "BTBD11_VALIDATION_PROVENANCE_AMENDMENT.json",
    HERE / "reseal_btbd11_repair.py",
    HERE / "finalize_btbd11_provenance.py",
    HERE / "work/btbd11-reseal-receipt.json",
    HERE / "work/btbd11-parent-binding-reconciliation-receipt.json",
)
BASE_ROLE_PROVENANCE_PATHS = (
    CONTRACT_PATH,
    PARENT_PROTOCOL_PATH,
    ROLE_BUILDER_PATH,
    Path(__file__).resolve(),
    *REPAIR_PROVENANCE_PATHS,
)
SOURCE_MODEL_PROVENANCE_PATHS = (
    SOURCE_DIRECTORY / "source-model.npz",
    SOURCE_DIRECTORY / "manifest.json",
    HERE / "model_input.py",
    HERE / "modeling.py",
    HERE / "fit_source.py",
    HERE / "source_artifacts.py",
)
PRODUCTION_IMPLEMENTATION_PATHS = (
    ROOT / "scripts/make_submission_pm_cd4.py",
    ROOT / "scripts/build_submission.py",
    ROOT / "scripts/make_submission.py",
    ROOT / "src/vcc2026/donor_transfer.py",
    ROOT / "src/vcc2026/multi_donor.py",
    ROOT / "src/vcc2026/whiten.py",
    ROOT / "experiments/061-perturbmatch-hct-count/build_candidate.py",
    ROOT / "experiments/062-cd4-post-whiten/run.py",
)
OPERATOR_IMPLEMENTATION_PATHS = (
    ROOT / "scripts/build_union_whitener.py",
    ROOT / "src/vcc2026/whiten.py",
)


class StrictValidationError(RuntimeError):
    """A frozen validation boundary or statistic contract was violated."""


@dataclass(frozen=True)
class FrozenSourceModel:
    directory: Path
    manifest: Mapping[str, Any]
    manifest_sha256: str
    model_sha256: str
    exported_effects: np.ndarray
    target_symbols: np.ndarray
    response_symbols: np.ndarray
    response_vcc_indices: np.ndarray
    response_exportable: np.ndarray
    transferable: np.ndarray


@dataclass(frozen=True)
class GalleryAxes:
    number: int
    targets: tuple[str, ...]
    target_indices_in_nadig: np.ndarray
    target_indices_in_source: np.ndarray
    source_transferable: np.ndarray
    score_genes: tuple[str, ...]
    score_indices_in_nadig: np.ndarray
    score_indices_in_source: np.ndarray
    score_indices_in_official: np.ndarray
    whitener_columns: np.ndarray
    basis_rows: Mapping[str, np.ndarray]


@dataclass(frozen=True)
class NestedSelection:
    fractions_by_held_gallery: tuple[float, ...]
    baseline_gallery_means: np.ndarray
    selected_gallery_means: np.ndarray
    held_gallery_gains: np.ndarray
    macro_baseline: float
    macro_selected: float
    macro_gain: float


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _reject_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON constant: {value}")


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(block_bytes), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_strings(values: Sequence[str]) -> str:
    payload = "".join(f"{value}\n" for value in values).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_int64(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values, dtype=np.int64))
    return hashlib.sha256(array.tobytes(order="C")).hexdigest()


def sha256_array(values: np.ndarray) -> str:
    """Hash an array with its dtype and shape so reorders cannot alias."""

    array = np.asarray(values)
    _require(array.dtype.hasobject is False, "object arrays cannot be bound")
    contiguous = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(b"exp107-array-v1\0")
    digest.update(contiguous.dtype.str.encode("ascii"))
    digest.update(b"\0")
    digest.update(
        json.dumps(list(contiguous.shape), separators=(",", ":")).encode("ascii")
    )
    digest.update(b"\0")
    digest.update(contiguous.tobytes(order="C"))
    return digest.hexdigest()


def array_binding(values: np.ndarray) -> dict[str, Any]:
    """Bind an array's exact logical representation, independent of NPZ encoding."""

    array = np.asarray(values)
    _require(array.dtype.hasobject is False, "object arrays cannot be bound")
    contiguous = np.ascontiguousarray(array)
    return {
        "dtype": contiguous.dtype.str,
        "shape": list(contiguous.shape),
        "c_order": True,
        "sha256_array": sha256_array(contiguous),
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


def _is_sha256(value: object) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _confined_regular_file(root: Path, relative: str, label: str) -> Path:
    pure = PurePosixPath(relative)
    _require(
        not pure.is_absolute()
        and relative == pure.as_posix()
        and pure.parts
        and ".." not in pure.parts,
        f"{label} path is not canonical repo-relative POSIX",
    )
    root = Path(root).resolve()
    candidate = root.joinpath(*pure.parts)
    cursor = root
    for part in pure.parts:
        cursor = cursor / part
        _require(not cursor.is_symlink(), f"{label} path traverses a symlink")
    try:
        candidate.resolve().relative_to(root)
    except ValueError as error:
        raise StrictValidationError(f"{label} escapes the binding root") from error
    _require(candidate.is_file(), f"{label} is missing or not a regular file")
    return candidate


def file_binding(path: Path, *, root: Path = ROOT) -> dict[str, Any]:
    root = Path(root).resolve()
    path = Path(path)
    try:
        relative = path.resolve().relative_to(root).as_posix()
    except ValueError as error:
        raise StrictValidationError(
            f"artifact is outside binding root: {path}"
        ) from error
    regular = _confined_regular_file(root, relative, "artifact")
    return {
        "path": relative,
        "bytes": regular.stat().st_size,
        "sha256": sha256_file(regular),
    }


def _verify_file_binding(
    expected: Mapping[str, Any],
    *,
    root: Path,
    label: str,
) -> dict[str, Any]:
    _require(
        set(expected) == {"path", "bytes", "sha256"},
        f"{label} file binding schema changed",
    )
    relative = expected.get("path")
    _require(isinstance(relative, str) and relative, f"{label} path is invalid")
    _require(
        isinstance(expected.get("bytes"), int) and int(expected["bytes"]) >= 0,
        f"{label} byte count is invalid",
    )
    _require(_is_sha256(expected.get("sha256")), f"{label} SHA-256 is invalid")
    path = _confined_regular_file(root, relative, label)
    observed = file_binding(path, root=root)
    _require(observed == dict(expected), f"{label} changed")
    return observed


def _role_kind(role: str) -> str:
    kind = role.split(".", 1)[0]
    _require(
        kind in {"donor", "operator", "incumbent", "integration"},
        f"unknown artifact role kind: {role}",
    )
    return kind


def _surface_contract(role: str, contract: Mapping[str, Any]) -> Mapping[str, Any]:
    parts = role.split(".")
    if len(parts) >= 3 and parts[1] == "nadig":
        _require(
            len(parts) == 3 and parts[2].startswith("gallery_"),
            f"invalid Nadig role: {role}",
        )
        try:
            number = int(parts[2].removeprefix("gallery_"))
        except ValueError as error:
            raise StrictValidationError(
                f"invalid Nadig gallery role: {role}"
            ) from error
        galleries = contract["gallery_construction"]["galleries"]
        _require(1 <= number <= len(galleries), f"gallery role is out of range: {role}")
        frozen = galleries[number - 1]
        _require(frozen["gallery"] == number, f"gallery contract order changed: {role}")
        return {
            **frozen,
            "surface": f"nadig.gallery_{number}",
            "target_count": contract["gallery_construction"]["targets_per_gallery"],
            "target_axis_sha256": frozen["target_axis_sha256_newline_terminated"],
            "score_axis_sha256": frozen["score_axis_sha256_newline_terminated"],
            "operator_donors": contract["incumbent"]["donor_order"],
        }
    _require(len(parts) == 2, f"invalid held-surface role: {role}")
    surface = parts[1]
    held = contract["held_context_gate"]["surfaces"]
    _require(surface in held, f"unknown held-surface role: {role}")
    return held[surface]


def _expected_parent_roles(role: str, contract: Mapping[str, Any]) -> tuple[str, ...]:
    kind = _role_kind(role)
    if kind == "donor":
        return ()
    surface = role.split(".", 1)[1]
    frozen = _surface_contract(role, contract)
    donors = tuple(f"donor.{name}" for name in frozen["operator_donors"])
    if kind == "operator":
        return donors
    if kind == "incumbent":
        return (*donors, f"operator.{surface}")
    return (f"incumbent.{surface}", f"operator.{surface}")


def expected_role_bundle_paths(
    role: str,
    *,
    bundle_root: Path = ROLE_BUNDLE_ROOT,
) -> tuple[Path, Path]:
    kind, suffix = role.split(".", 1)
    directory = Path(bundle_root) / suffix.replace(".", "/")
    return directory / f"{kind}.npz", directory / f"{kind}.json"


def _surface_provenance_paths(
    role: str,
    contract: Mapping[str, Any],
) -> tuple[Path, ...]:
    if ".nadig.gallery_" in role:
        return (
            NADIG_DIRECTORY / "perts.txt",
            NADIG_DIRECTORY / "genes.txt",
        )
    surface = _surface_contract(role, contract)
    return (ROOT / str(surface["harness_path"]),)


def expected_transitive_paths(
    role: str,
    contract: Mapping[str, Any],
) -> tuple[Path, ...]:
    """Return the exact production transitive-input set for one role bundle."""

    kind = _role_kind(role)
    if kind == "donor":
        donor = role.split(".", 1)[1]
        return (
            *BASE_ROLE_PROVENANCE_PATHS,
            *(
                DONOR_ROOT / donor / filename
                for filename in (
                    "delta.npy",
                    "b_ctrl.npy",
                    "genes.txt",
                    "perts.txt",
                    "meta.json",
                )
            ),
        )
    frozen = _surface_contract(role, contract)
    donor_order = tuple(frozen["operator_donors"])
    parent_bundle_paths = tuple(
        path
        for parent in _expected_parent_roles(role, contract)
        for path in expected_role_bundle_paths(parent)
    )
    surface_paths = _surface_provenance_paths(role, contract)
    if kind == "operator":
        donor_inputs = tuple(
            DONOR_ROOT / donor / filename
            for donor in donor_order
            for filename in ("delta.npy", "genes.txt", "perts.txt")
        )
        return (
            *BASE_ROLE_PROVENANCE_PATHS,
            *SOURCE_MODEL_PROVENANCE_PATHS,
            CHALLENGE_DIRECTORY / "gene_names.csv",
            CHALLENGE_DIRECTORY / "pert_counts.csv",
            DONOR_ROOT / "k562_essential/perts.txt",
            *OPERATOR_IMPLEMENTATION_PATHS,
            *donor_inputs,
            *parent_bundle_paths,
            *surface_paths,
        )
    if kind == "incumbent":
        donor_inputs = tuple(
            DONOR_ROOT / donor / filename
            for donor in donor_order
            for filename in (
                "delta.npy",
                "b_ctrl.npy",
                "genes.txt",
                "perts.txt",
                "meta.json",
            )
        )
        perturbmatch = tuple(
            ROOT / "data/external/perturbmatch" / filename
            for filename in ("delta.npy", "genes.txt", "perts.txt", "meta.json")
        )
        cd4 = tuple(
            ROOT / "data/external/cd4t/cloud-compact/cd4_compact" / condition / filename
            for condition in ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
            for filename in (
                "delta.npy",
                "b_ctrl.npy",
                "genes.txt",
                "perts.txt",
                "meta.json",
            )
        )
        return (
            *BASE_ROLE_PROVENANCE_PATHS,
            *SOURCE_MODEL_PROVENANCE_PATHS,
            CHALLENGE_DIRECTORY / "manifest.json",
            CHALLENGE_DIRECTORY / "gene_names.csv",
            CHALLENGE_DIRECTORY / "pert_counts.csv",
            *PRODUCTION_IMPLEMENTATION_PATHS,
            *donor_inputs,
            *perturbmatch,
            *cd4,
            *parent_bundle_paths,
            *surface_paths,
        )
    return (
        *BASE_ROLE_PROVENANCE_PATHS,
        *SOURCE_MODEL_PROVENANCE_PATHS,
        *parent_bundle_paths,
        *surface_paths,
    )


def _require_array(
    arrays: Mapping[str, np.ndarray],
    bindings: Mapping[str, Any],
    name: str,
    *,
    dtype: str | None = None,
    dtype_kind: str | None = None,
    shape: tuple[int, ...] | None = None,
) -> np.ndarray:
    _require(name in arrays, f"role bundle array is missing: {name}")
    value = np.asarray(arrays[name])
    _require(name in bindings, f"role bundle array binding is missing: {name}")
    _require(
        set(bindings[name]) == {"dtype", "shape", "c_order", "sha256_array"},
        f"role bundle array binding schema changed: {name}",
    )
    _require(
        bindings[name] == array_binding(value), f"role bundle array changed: {name}"
    )
    _require(value.flags.c_contiguous, f"role bundle array is not C-order: {name}")
    if dtype is not None:
        _require(value.dtype.str == dtype, f"role bundle array dtype changed: {name}")
    if dtype_kind is not None:
        _require(
            value.dtype.kind == dtype_kind, f"role bundle array kind changed: {name}"
        )
    if shape is not None:
        _require(value.shape == shape, f"role bundle array shape changed: {name}")
    return value


def validate_role_bundle(
    role: str,
    paths: Sequence[str],
    contract: Mapping[str, Any],
    *,
    root: Path = ROOT,
) -> dict[str, Any]:
    """Validate one role-owned NPZ+JSON bundle and all transitive input hashes."""

    _require(
        len(paths) == 2 and paths[0].endswith(".npz") and paths[1].endswith(".json"),
        f"artifact role {role} must name exactly one NPZ followed by one JSON manifest",
    )
    root = Path(root).resolve()
    npz_path = _confined_regular_file(root, paths[0], f"artifact role {role} NPZ")
    manifest_path = _confined_regular_file(
        root, paths[1], f"artifact role {role} manifest"
    )
    if root == ROOT.resolve():
        expected_paths = tuple(
            path.relative_to(ROOT).as_posix()
            for path in expected_role_bundle_paths(role)
        )
        _require(
            tuple(paths) == expected_paths,
            f"artifact role {role} bundle path changed",
        )
    manifest = load_json(manifest_path)
    _require(
        set(manifest)
        == {
            "schema_version",
            "experiment",
            "status",
            "role",
            "kind",
            "npz",
            "arrays",
            "parent_roles",
            "transitive_inputs",
            "parameters",
        },
        f"artifact role {role} manifest schema changed",
    )
    _require_exact(manifest, "schema_version", ROLE_BUNDLE_SCHEMA_VERSION)
    _require_exact(manifest, "experiment", "107-jurkat-high-moi")
    _require_exact(manifest, "status", ROLE_BUNDLE_STATUS)
    _require_exact(manifest, "role", role)
    kind = _role_kind(role)
    _require_exact(manifest, "kind", kind)
    _require(
        manifest.get("npz") == file_binding(npz_path, root=root),
        f"artifact role {role} NPZ binding changed",
    )
    parents = manifest.get("parent_roles")
    _require(
        parents == list(_expected_parent_roles(role, contract)),
        f"artifact role {role} parent roles changed",
    )
    transitive = manifest.get("transitive_inputs")
    _require(
        isinstance(transitive, list) and transitive,
        f"artifact role {role} has no transitive input bindings",
    )
    transitive_paths: set[str] = set()
    for index, binding in enumerate(transitive):
        _require(
            isinstance(binding, dict),
            f"artifact role {role} transitive binding is invalid",
        )
        observed = _verify_file_binding(
            binding,
            root=root,
            label=f"artifact role {role} transitive input {index}",
        )
        _require(
            observed["path"] not in transitive_paths,
            f"artifact role {role} repeats a transitive input",
        )
        transitive_paths.add(observed["path"])
    if root == ROOT.resolve():
        expected_transitive = {
            path.relative_to(ROOT).as_posix()
            for path in expected_transitive_paths(role, contract)
        }
        _require(
            transitive_paths == expected_transitive,
            f"artifact role {role} transitive input set changed",
        )
    _require(
        isinstance(manifest.get("parameters"), dict),
        f"artifact role {role} parameters are invalid",
    )
    bindings = manifest.get("arrays")
    _require(
        isinstance(bindings, dict), f"artifact role {role} array bindings are missing"
    )
    try:
        with np.load(npz_path, allow_pickle=False) as archive:
            arrays = {name: np.asarray(archive[name]) for name in archive.files}
    except (OSError, ValueError, KeyError) as error:
        raise StrictValidationError(
            f"cannot load artifact role {role}: {error}"
        ) from error
    _require(set(arrays) == set(bindings), f"artifact role {role} array schema changed")

    if kind == "donor":
        expected_keys = {"genes", "perts", "delta_shape", "b_ctrl_shape"}
        _require(set(arrays) == expected_keys, f"donor role {role} array keys changed")
        genes = _require_array(arrays, bindings, "genes", dtype_kind="U")
        perts = _require_array(arrays, bindings, "perts", dtype_kind="U")
        delta_shape = _require_array(
            arrays, bindings, "delta_shape", dtype="<i8", shape=(2,)
        )
        b_ctrl_shape = _require_array(
            arrays, bindings, "b_ctrl_shape", dtype="<i8", shape=(1,)
        )
        _require(
            genes.ndim == perts.ndim == 1, f"donor role {role} axes are not vectors"
        )
        _require(
            tuple(map(int, delta_shape)) == (len(perts), len(genes))
            and int(b_ctrl_shape[0]) == len(genes),
            f"donor role {role} dimensions do not align",
        )
    else:
        surface = _surface_contract(role, contract)
        targets = int(surface["target_count"])
        score_genes = int(surface["score_gene_count"])
        if kind == "operator":
            columns = int(surface["whitener_column_count"])
            donors = tuple(surface["operator_donors"])
            expected_keys = {"evals", "evecs", "columns"} | {
                f"basis_{name}" for name in donors
            }
            _require(
                set(arrays) == expected_keys, f"operator role {role} array keys changed"
            )
            evals = _require_array(
                arrays, bindings, "evals", dtype="<f8", shape=(columns,)
            )
            evecs = _require_array(
                arrays, bindings, "evecs", dtype="<f4", shape=(columns, columns)
            )
            column_values = _require_array(
                arrays, bindings, "columns", dtype="<i8", shape=(columns,)
            )
            _require(
                sha256_int64(column_values) == surface["whitener_columns_int64_sha256"],
                f"operator role {role} column axis changed",
            )
            _require(
                np.all(np.isfinite(evals)) and np.all(evals > 0),
                f"operator role {role} eigenvalues are invalid",
            )
            _require(
                np.all(np.isfinite(evecs)),
                f"operator role {role} eigenvectors are non-finite",
            )
            _require(
                np.array_equal(column_values, np.unique(column_values))
                and np.all(column_values >= 0)
                and np.all(column_values < OFFICIAL_GENE_COUNT),
                f"operator role {role} columns are not sorted unique indices",
            )
            probe = np.linspace(0, columns - 1, num=min(32, columns), dtype=np.int64)
            gram = evecs[:, probe].astype(np.float64).T @ evecs[:, probe].astype(
                np.float64
            )
            _require(
                np.allclose(gram, np.eye(len(probe)), rtol=5e-5, atol=5e-5),
                f"operator role {role} sampled eigenvectors are not orthonormal",
            )
            expected_basis = surface["whitener_basis_rows_by_donor"]
            for donor in donors:
                values = _require_array(
                    arrays,
                    bindings,
                    f"basis_{donor}",
                    dtype="<i8",
                    shape=(int(expected_basis[donor]["count"]),),
                )
                _require(
                    sha256_int64(values) == expected_basis[donor]["int64_sha256"],
                    f"operator role {role} basis rows changed for {donor}",
                )
                _require(
                    np.array_equal(values, np.unique(values)) and np.all(values >= 0),
                    f"operator role {role} basis rows are not sorted unique for {donor}",
                )
        elif kind == "incumbent":
            expected_keys = {
                "targets",
                "genes",
                "score_genes",
                "score_indices",
                "effects",
            }
            _require(
                set(arrays) == expected_keys,
                f"incumbent role {role} array keys changed",
            )
            target_axis = _require_array(
                arrays, bindings, "targets", dtype_kind="U", shape=(targets,)
            )
            _require(
                sha256_strings(target_axis.astype(str))
                == surface["target_axis_sha256"],
                f"incumbent role {role} target axis changed",
            )
            genes = _require_array(
                arrays,
                bindings,
                "genes",
                dtype_kind="U",
                shape=(OFFICIAL_GENE_COUNT,),
            )
            score_axis = _require_array(
                arrays,
                bindings,
                "score_genes",
                dtype_kind="U",
                shape=(score_genes,),
            )
            _require(
                sha256_strings(score_axis.astype(str)) == surface["score_axis_sha256"],
                f"incumbent role {role} score axis changed",
            )
            score_indices = _require_array(
                arrays,
                bindings,
                "score_indices",
                dtype="<i8",
                shape=(score_genes,),
            )
            _require(
                len(np.unique(score_indices)) == len(score_indices)
                and np.all((0 <= score_indices) & (score_indices < len(genes))),
                f"incumbent role {role} score indices are invalid",
            )
            effects = _require_array(
                arrays,
                bindings,
                "effects",
                dtype="<f8",
                shape=(targets, OFFICIAL_GENE_COUNT),
            )
            _require(
                len(set(genes.astype(str))) == len(genes)
                and len(set(target_axis.astype(str))) == len(target_axis)
                and len(set(score_axis.astype(str))) == len(score_axis),
                f"incumbent role {role} contains duplicate axis symbols",
            )
            _require(
                np.array_equal(
                    genes[score_indices].astype(str), score_axis.astype(str)
                ),
                f"incumbent role {role} score indices do not map its score axis",
            )
            _require(
                np.all(np.isfinite(effects)),
                f"incumbent role {role} effects are non-finite",
            )
            bundle_contract = contract.get("role_bundle_contract", {})
            expected_gene_hash = bundle_contract.get(
                "official_gene_axis_sha256_newline_terminated"
            )
            if expected_gene_hash is not None:
                _require(
                    sha256_strings(genes.astype(str)) == expected_gene_hash,
                    f"incumbent role {role} official gene axis changed",
                )
        else:
            expected_keys = {
                "targets",
                "genes",
                "score_genes",
                "score_indices",
                "incumbent_effects",
                "source_effects",
                "transferable",
                "coordinate_supported",
                *{f"candidate_f{int(round(value * 1000)):03d}" for value in FRACTIONS},
                *{f"effective_f{int(round(value * 1000)):03d}" for value in FRACTIONS},
            }
            _require(
                set(arrays) == expected_keys,
                f"integration role {role} array keys changed",
            )
            target_axis = _require_array(
                arrays, bindings, "targets", dtype_kind="U", shape=(targets,)
            )
            _require(
                sha256_strings(target_axis.astype(str))
                == surface["target_axis_sha256"],
                f"integration role {role} target axis changed",
            )
            genes = _require_array(
                arrays,
                bindings,
                "genes",
                dtype_kind="U",
                shape=(OFFICIAL_GENE_COUNT,),
            )
            score_axis = _require_array(
                arrays,
                bindings,
                "score_genes",
                dtype_kind="U",
                shape=(score_genes,),
            )
            _require(
                sha256_strings(score_axis.astype(str)) == surface["score_axis_sha256"],
                f"integration role {role} score axis changed",
            )
            score_indices = _require_array(
                arrays,
                bindings,
                "score_indices",
                dtype="<i8",
                shape=(score_genes,),
            )
            _require(
                len(np.unique(score_indices)) == len(score_indices)
                and np.all((0 <= score_indices) & (score_indices < len(genes))),
                f"integration role {role} score indices are invalid",
            )
            incumbent = _require_array(
                arrays,
                bindings,
                "incumbent_effects",
                dtype="<f8",
                shape=(targets, OFFICIAL_GENE_COUNT),
            )
            source_effects = _require_array(
                arrays,
                bindings,
                "source_effects",
                dtype="<f8",
                shape=(targets, OFFICIAL_GENE_COUNT),
            )
            transferable = _require_array(
                arrays,
                bindings,
                "transferable",
                dtype="|b1",
                shape=(targets,),
            )
            coordinate_supported = _require_array(
                arrays,
                bindings,
                "coordinate_supported",
                dtype="|b1",
                shape=(OFFICIAL_GENE_COUNT,),
            )
            _require(
                len(set(genes.astype(str))) == len(genes)
                and len(set(target_axis.astype(str))) == len(target_axis)
                and len(set(score_axis.astype(str))) == len(score_axis),
                f"integration role {role} contains duplicate axis symbols",
            )
            _require(
                np.array_equal(
                    genes[score_indices].astype(str), score_axis.astype(str)
                ),
                f"integration role {role} score indices do not map its score axis",
            )
            expected_coordinates = np.zeros(OFFICIAL_GENE_COUNT, dtype=bool)
            expected_coordinates[score_indices] = True
            _require(
                np.array_equal(coordinate_supported, expected_coordinates),
                f"integration role {role} coordinate support changed",
            )
            _require(
                np.all(np.isfinite(incumbent)) and np.all(np.isfinite(source_effects)),
                f"integration role {role} effects are non-finite",
            )
            _require(
                not np.any(
                    source_effects[
                        ~(transferable[:, None] & coordinate_supported[None, :])
                    ]
                ),
                f"integration role {role} source changed an unsupported entry",
            )
            _require(
                int(transferable.sum()) == int(surface["transferable_target_count"]),
                f"integration role {role} transferable target count changed",
            )
            bundle_contract = contract.get("role_bundle_contract", {})
            expected_gene_hash = bundle_contract.get(
                "official_gene_axis_sha256_newline_terminated"
            )
            if expected_gene_hash is not None:
                _require(
                    sha256_strings(genes.astype(str)) == expected_gene_hash,
                    f"integration role {role} official gene axis changed",
                )
            for fraction in FRACTIONS:
                suffix = f"f{int(round(fraction * 1000)):03d}"
                candidate = _require_array(
                    arrays,
                    bindings,
                    f"candidate_{suffix}",
                    dtype="<f8",
                    shape=(targets, OFFICIAL_GENE_COUNT),
                )
                effective = _require_array(
                    arrays,
                    bindings,
                    f"effective_{suffix}",
                    dtype="|b1",
                    shape=(targets,),
                )
                assert_exact_fallback(
                    incumbent,
                    candidate,
                    effective if fraction else np.zeros_like(transferable),
                    coordinate_supported,
                )
                if fraction == 0.0:
                    _require(
                        np.array_equal(
                            candidate.view(np.uint8), incumbent.view(np.uint8)
                        )
                        and not np.any(effective),
                        f"integration role {role} fraction-zero branch changed",
                    )
                else:
                    _require(
                        np.all(~effective | transferable),
                        f"integration role {role} effective support exceeds source support",
                    )
    return {
        "role": role,
        "kind": kind,
        "npz": file_binding(npz_path, root=root),
        "manifest": file_binding(manifest_path, root=root),
        "arrays": bindings,
        "parent_roles": parents,
        "transitive_inputs": transitive,
    }


def bind_artifact_manifest(
    manifest_path: Path,
    *,
    required_roles: Sequence[str] = REQUIRED_ARTIFACT_ROLES,
    root: Path = ROOT,
    contract: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Bind every exact file named by the frozen validation-input manifest."""

    root = Path(root).resolve()
    manifest_path = Path(manifest_path)
    manifest_binding = file_binding(manifest_path, root=root)
    manifest = load_json(manifest_path)
    _require_exact(manifest, "schema_version", 1)
    _require_exact(manifest, "experiment", "107-jurkat-high-moi")
    _require_exact(
        manifest, "status", "frozen_before_experiment_107_nadig_outcome_access"
    )
    artifacts = manifest.get("artifacts")
    _require(isinstance(artifacts, dict), "artifact manifest artifacts are missing")
    roles = tuple(required_roles)
    _require(
        set(artifacts) == set(roles),
        "artifact manifest role set is incomplete or changed",
    )
    bound_roles: dict[str, dict[str, Any]] = {}
    seen_paths: set[str] = set()
    frozen_contract = validate_contract() if contract is None else contract
    for role in roles:
        paths = artifacts[role]
        _require(
            isinstance(paths, list)
            and all(isinstance(value, str) and value for value in paths),
            f"artifact role {role} must contain a path list",
        )
        _require(len(paths) == len(set(paths)), f"artifact role {role} repeats a path")
        entries: list[dict[str, Any]] = []
        for relative in paths:
            _require(
                relative not in seen_paths,
                f"artifact path is assigned twice: {relative}",
            )
            seen_paths.add(relative)
            path = _confined_regular_file(root, relative, f"artifact role {role}")
            entries.append(file_binding(path, root=root))
        bundle = validate_role_bundle(role, paths, frozen_contract, root=root)
        bound_roles[role] = {"files": entries, "bundle": bundle}
    for role, value in bound_roles.items():
        for parent in value["bundle"]["parent_roles"]:
            _require(
                parent in bound_roles,
                f"artifact role {role} parent is missing: {parent}",
            )
            parent_paths = {
                bound_roles[parent]["bundle"]["npz"]["path"],
                bound_roles[parent]["bundle"]["manifest"]["path"],
            }
            transitive_paths = {
                entry["path"] for entry in value["bundle"]["transitive_inputs"]
            }
            _require(
                parent_paths.issubset(transitive_paths),
                f"artifact role {role} does not bind parent bundle {parent}",
            )
        if value["bundle"]["kind"] == "integration":
            surface = role.split(".", 1)[1]
            incumbent_role = f"incumbent.{surface}"
            incumbent_hash = bound_roles[incumbent_role]["bundle"]["arrays"]["effects"][
                "sha256_array"
            ]
            _require(
                value["bundle"]["arrays"]["incumbent_effects"]["sha256_array"]
                == incumbent_hash,
                f"artifact role {role} does not embed the exact incumbent",
            )
            operator_role = f"operator.{surface}"
            integration_path = root / artifacts[role][0]
            operator_path = root / artifacts[operator_role][0]
            integration_manifest = load_json(root / artifacts[role][1])
            with (
                np.load(integration_path, allow_pickle=False) as integration,
                np.load(operator_path, allow_pickle=False) as operator,
            ):
                assert_operator_missing_score_fallback(
                    genes=np.asarray(integration["genes"]),
                    score_indices=np.asarray(
                        integration["score_indices"], dtype=np.int64
                    ),
                    operator_columns=np.asarray(
                        operator["columns"], dtype=np.int64
                    ),
                    incumbent=np.asarray(
                        integration["incumbent_effects"], dtype=np.float64
                    ),
                    source=np.asarray(
                        integration["source_effects"], dtype=np.float64
                    ),
                    candidates={
                        float(fraction): np.asarray(
                            integration[
                                f"candidate_f{int(round(fraction * 1000)):03d}"
                            ],
                            dtype=np.float64,
                        )
                        for fraction in FRACTIONS
                    },
                    parameters=integration_manifest["parameters"],
                    label=f"artifact role {role}",
                )
    payload = {
        "manifest": manifest_binding,
        "roles": bound_roles,
    }
    payload["binding_sha256"] = sha256_json_payload(payload)
    return payload


def verify_artifact_bindings(
    expected: Mapping[str, Any],
    manifest_path: Path,
    *,
    required_roles: Sequence[str] = REQUIRED_ARTIFACT_ROLES,
    root: Path = ROOT,
    contract: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    observed = bind_artifact_manifest(
        manifest_path,
        required_roles=required_roles,
        root=root,
        contract=contract,
    )
    _require(observed == expected, "bound artifact files changed or were reordered")
    return observed


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(
            Path(path).read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
            parse_constant=_reject_constant,
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        raise StrictValidationError(
            f"cannot read strict JSON {path}: {error}"
        ) from error
    if not isinstance(value, dict):
        raise StrictValidationError(f"{path} must contain one JSON object")
    return value


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise StrictValidationError(message)


def _require_exact(mapping: Mapping[str, Any], key: str, expected: Any) -> None:
    _require(mapping.get(key) == expected, f"{key} changed from its frozen value")


def validate_contract(path: Path = CONTRACT_PATH) -> dict[str, Any]:
    """Validate the byte-frozen contract and its parent protocol."""

    path = Path(path)
    _require(path.is_file() and not path.is_symlink(), "contract is missing or unsafe")
    _require(
        sha256_file(path) == EXPECTED_CONTRACT_SHA256,
        "strict validation contract bytes changed",
    )
    contract = load_json(path)
    _require_exact(contract, "schema_version", 3)
    _require_exact(contract, "experiment", "107-jurkat-high-moi")
    _require_exact(
        contract,
        "status",
        "repaired_before_experiment_107_nadig_outcome_access",
    )
    _require_exact(
        contract,
        "benchmark_role",
        "reused_development_benchmark_not_independent_confirmation",
    )
    parent = contract.get("parent_protocol")
    _require(isinstance(parent, dict), "parent_protocol is missing")
    parent_path = ROOT / str(parent.get("path", ""))
    _require(parent_path.is_file(), "parent protocol is missing")
    _require(
        sha256_file(parent_path) == parent.get("sha256"),
        "parent protocol hash changed",
    )
    held = contract.get("held_truth")
    _require(isinstance(held, dict), "held_truth is missing")
    _require_exact(held, "effect_file", "delta.npy")
    _require_exact(held, "target_axis_file", "perts.txt")
    _require_exact(held, "gene_axis_file", "genes.txt")
    _require_exact(held, "effect_values_forbidden_before_source_model_freeze", True)
    _require_exact(held, "control_profile_forbidden_for_selection", True)
    _require_exact(held, "metadata_effect_norm_forbidden_for_selection", True)
    _require_exact(
        contract["source_integration"],
        "blend_grid",
        list(FRACTIONS),
    )
    bundle_contract = contract.get("role_bundle_contract")
    _require(isinstance(bundle_contract, dict), "role_bundle_contract is missing")
    _require_exact(bundle_contract, "schema_version", ROLE_BUNDLE_SCHEMA_VERSION)
    _require_exact(bundle_contract, "official_gene_count", OFFICIAL_GENE_COUNT)
    _require(
        _is_sha256(bundle_contract.get("official_gene_axis_sha256_newline_terminated")),
        "official gene-axis hash is invalid",
    )
    _require_exact(bundle_contract, "effect_dtype", "<f8")
    _require_exact(bundle_contract, "eigenvalue_dtype", "<f8")
    _require_exact(bundle_contract, "eigenvector_dtype", "<f4")
    _require_exact(bundle_contract, "index_dtype", "<i8")
    _require_exact(bundle_contract, "mask_dtype", "|b1")
    _require_exact(bundle_contract, "all_arrays_c_order", True)
    _require_exact(
        contract["target_label_permutation_null"],
        "replicates",
        PERMUTATION_REPLICATES,
    )
    _require_exact(
        contract["cluster_bootstrap"],
        "replicates",
        BOOTSTRAP_REPLICATES,
    )
    authorization = contract.get("authorization_bindings")
    _require(isinstance(authorization, dict), "authorization_bindings is missing")
    _require_exact(
        authorization,
        "validation_implementation",
        "experiments/107-jurkat-high-moi/strict_nadig_validation.py",
    )
    _require_exact(
        authorization,
        "artifact_manifest",
        ARTIFACT_MANIFEST_PATH.relative_to(ROOT).as_posix(),
    )
    _require_exact(authorization, "artifact_manifest_schema_version", 1)
    _require_exact(
        authorization,
        "artifact_manifest_status",
        ROLE_BUNDLE_STATUS,
    )
    _require_exact(
        authorization,
        "role_value",
        "exactly_two_unique_role_owned_files_in_npz_then_json_order",
    )
    _require_exact(
        authorization,
        "role_bundle_schema_version",
        ROLE_BUNDLE_SCHEMA_VERSION,
    )
    _require_exact(authorization, "role_bundle_status", ROLE_BUNDLE_STATUS)
    _require_exact(
        authorization,
        "required_artifact_roles",
        list(REQUIRED_ARTIFACT_ROLES),
    )
    _require_exact(
        authorization,
        "nadig_artifacts_hashed_before_array_open",
        ["perts.txt", "genes.txt", "delta.npy"],
    )
    _require_exact(
        authorization,
        "outcome_access_receipt_policy",
        "exclusive_create_before_numpy_load_and_reject_any_replay",
    )
    _require_exact(
        authorization,
        "outcome_access_receipt",
        OUTCOME_ACCESS_RECEIPT_PATH.relative_to(ROOT).as_posix(),
    )
    _require_exact(
        authorization,
        "pretruth_plan",
        PRETRUTH_PLAN_PATH.relative_to(ROOT).as_posix(),
    )
    held_gate = contract.get("held_context_gate")
    _require(isinstance(held_gate, dict), "held_context_gate is missing")
    _require_exact(
        held_gate,
        "authorizing_directional_surfaces",
        list(AUTHORIZING_HELD_SURFACES),
    )
    _require_exact(held_gate, "diagnostic_directional_surface", None)
    surfaces = held_gate.get("surfaces")
    _require(
        isinstance(surfaces, dict) and tuple(surfaces) == AUTHORIZING_HELD_SURFACES,
        "held-surface contract changed",
    )
    for surface in AUTHORIZING_HELD_SURFACES:
        frozen = surfaces[surface]
        _require(
            isinstance(frozen, dict)
            and frozen.get("surface") == surface
            and isinstance(frozen.get("target_count"), int)
            and int(frozen["target_count"]) > 0
            and _is_sha256(frozen.get("target_axis_sha256"))
            and isinstance(frozen.get("score_gene_count"), int)
            and int(frozen["score_gene_count"]) > 0
            and _is_sha256(frozen.get("score_axis_sha256"))
            and isinstance(frozen.get("operator_donors"), list)
            and frozen["operator_donors"]
            and set(frozen["operator_donors"]).issubset(DONOR_ORDER)
            and isinstance(frozen.get("whitener_column_count"), int)
            and int(frozen["whitener_column_count"]) > 0
            and _is_sha256(frozen.get("whitener_columns_int64_sha256"))
            and isinstance(frozen.get("whitener_basis_rows_by_donor"), dict)
            and set(frozen["whitener_basis_rows_by_donor"])
            == set(frozen["operator_donors"]),
            f"held-surface contract is incomplete: {surface}",
        )
    _require_exact(
        contract["promotion_gates"],
        "nadig_alone_cannot_authorize_promotion",
        True,
    )
    return contract


def _read_axis(
    path: Path,
    label: str,
    *,
    require_unique: bool = True,
) -> tuple[str, ...]:
    path = Path(path)
    _require(path.name in {"perts.txt", "genes.txt"}, f"unsafe {label} filename")
    _require(path.is_file() and not path.is_symlink(), f"{label} is missing or unsafe")
    try:
        values = tuple(path.read_text(encoding="utf-8").splitlines())
    except (OSError, UnicodeError) as error:
        raise StrictValidationError(f"cannot read {label}: {error}") from error
    _require(values and all(values), f"{label} contains an empty symbol")
    if require_unique:
        _require(len(values) == len(set(values)), f"{label} contains duplicate symbols")
    return values


def _read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    try:
        with Path(path).open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            _require(reader.fieldnames == [column], f"{path} header changed")
            values = tuple(str(row[column]) for row in reader)
    except (OSError, UnicodeError, csv.Error) as error:
        raise StrictValidationError(f"cannot read {path}: {error}") from error
    _require(values and all(values), f"{path} contains an empty value")
    _require(len(values) == len(set(values)), f"{path} contains duplicate values")
    return values


def _validate_source_arrays(
    arrays: Mapping[str, np.ndarray],
    manifest: Mapping[str, Any],
    *,
    expected_targets: int,
    expected_responses: int,
    expected_exportable: int,
) -> None:
    _require(set(arrays) == SOURCE_MODEL_KEYS, "source model keys changed")
    targets = np.asarray(arrays["target_symbols"]).astype(str)
    responses = np.asarray(arrays["response_symbols"]).astype(str)
    effects = np.asarray(arrays["exported_effects"])
    transferable = np.asarray(arrays["transferable"])
    exportable = np.asarray(arrays["response_exportable"])
    vcc_indices = np.asarray(arrays["response_vcc_indices"])
    _require(targets.shape == (expected_targets,), "source target axis changed")
    _require(responses.shape == (expected_responses,), "source response axis changed")
    _require(
        effects.shape == (expected_targets, expected_responses),
        "source effect shape changed",
    )
    _require(
        transferable.shape == (expected_targets,), "transferable mask shape changed"
    )
    _require(transferable.dtype == np.bool_, "transferable mask is not boolean")
    _require(exportable.shape == (expected_responses,), "exportable mask shape changed")
    _require(exportable.dtype == np.bool_, "exportable mask is not boolean")
    _require(vcc_indices.shape == (expected_responses,), "response index shape changed")
    _require(
        np.issubdtype(vcc_indices.dtype, np.integer),
        "response indices are not integers",
    )
    _require(
        np.all(vcc_indices >= -1), "response indices contain invalid negative values"
    )
    _require(
        np.all(vcc_indices[exportable] >= 0), "exportable responses lack VCC indices"
    )
    _require(
        len(set(targets.tolist())) == len(targets), "source targets are not unique"
    )
    _require(
        len(set(responses.tolist())) == len(responses),
        "source responses are not unique",
    )
    _require(
        np.all(targets != "") and np.all(responses != ""),
        "source axes contain empty symbols",
    )
    _require(np.all(np.isfinite(effects)), "source effects contain non-finite values")
    _require(
        np.count_nonzero(exportable) == expected_exportable,
        "exportable response count changed",
    )
    _require(
        np.array_equal(effects[~transferable], np.zeros_like(effects[~transferable])),
        "unsupported source rows are not exact zero fallback",
    )
    _require(
        manifest.get("targets") == expected_targets, "manifest target count changed"
    )
    _require(
        manifest.get("response_genes") == expected_responses,
        "manifest response count changed",
    )
    _require(
        manifest.get("exportable_response_genes") == expected_exportable,
        "manifest exportable response count changed",
    )
    final_support = manifest.get("final_support")
    _require(isinstance(final_support, dict), "source final_support is missing")
    _require(
        final_support.get("transferable_targets")
        == int(np.count_nonzero(transferable)),
        "source transferable count disagrees with manifest",
    )
    selected_lambda = float(np.asarray(arrays["selected_lambda"]).reshape(()))
    _require(math.isfinite(selected_lambda), "selected source lambda is non-finite")
    _require(
        manifest.get("selected_lambda") == selected_lambda,
        "selected source lambda changed",
    )
    lambda_grid = np.asarray(arrays["lambda_grid"], dtype=np.float64)
    _require(
        np.array_equal(lambda_grid, np.asarray([1.0, 3.0, 10.0, 30.0, 100.0])),
        "source lambda grid changed",
    )
    _require(
        manifest.get("lambda_grid") == lambda_grid.tolist(),
        "manifest lambda grid changed",
    )
    held_folds = manifest.get("held_channel_folds")
    _require(
        isinstance(held_folds, list)
        and len(held_folds) == 4
        and all(isinstance(fold, list) and len(fold) == 4 for fold in held_folds),
        "source held-channel fold structure changed",
    )
    fold_count = len(held_folds)
    lambda_count = lambda_grid.size
    for name in (
        "evidence_support",
        "zero_variance",
        "nuisance_aliased",
        "fitted",
        "projection_diagonal",
        "target_scales",
    ):
        _require(
            np.asarray(arrays[name]).shape == (expected_targets,),
            f"{name} shape changed",
        )
    for name in ("coefficients_binary_units", "projected_coefficients"):
        _require(
            np.asarray(arrays[name]).shape == effects.shape,
            f"{name} shape changed",
        )
        _require(np.all(np.isfinite(arrays[name])), f"{name} is non-finite")
    for name in ("evidence_support", "zero_variance", "nuisance_aliased", "fitted"):
        _require(np.asarray(arrays[name]).dtype == np.bool_, f"{name} is not boolean")
    held_cell_counts = np.asarray(arrays["held_cell_counts"])
    fold_mean_gains = np.asarray(arrays["fold_mean_gains"])
    mean_held_cell_gains = np.asarray(arrays["mean_held_cell_gains"])
    _require(
        held_cell_counts.shape == (fold_count,)
        and np.issubdtype(held_cell_counts.dtype, np.integer)
        and np.all(held_cell_counts > 0),
        "held cell counts changed",
    )
    _require(
        fold_mean_gains.shape == (fold_count, lambda_count)
        and np.all(np.isfinite(fold_mean_gains)),
        "fold mean gains changed",
    )
    _require(
        mean_held_cell_gains.shape == (lambda_count,)
        and np.all(np.isfinite(mean_held_cell_gains)),
        "mean held-cell gains changed",
    )
    _require(
        np.asarray(arrays["fold_selected_effects"]).shape
        == (fold_count, expected_targets, expected_responses),
        "fold selected effects shape changed",
    )
    _require(
        np.all(np.isfinite(arrays["fold_selected_effects"])),
        "fold selected effects are non-finite",
    )
    for name in (
        "fold_zero_variance",
        "fold_nuisance_aliased",
        "fold_transferable",
        "fold_support_supported",
    ):
        value = np.asarray(arrays[name])
        _require(
            value.shape == (fold_count, expected_targets) and value.dtype == np.bool_,
            f"{name} shape or dtype changed",
        )
    expected_training_channels = np.asarray(
        [
            [
                channel
                for channel in range(1, 17)
                if channel not in set(map(int, held_fold))
            ]
            for held_fold in held_folds
        ],
        dtype=np.int64,
    )
    fold_training_channels = np.asarray(arrays["fold_support_training_channels"])
    _require(
        np.issubdtype(fold_training_channels.dtype, np.integer)
        and np.array_equal(fold_training_channels, expected_training_channels),
        "fold support training-channel axes changed",
    )
    for name in (
        "fold_support_distinct_called_guides",
        "fold_support_guides_with_three_untreated_cells",
        "fold_support_untreated_carrier_cells",
        "fold_support_training_channel_count",
    ):
        value = np.asarray(arrays[name])
        _require(
            value.shape == (fold_count, expected_targets)
            and np.issubdtype(value.dtype, np.integer)
            and np.all(value >= 0),
            f"{name} shape, dtype, or values changed",
        )
    final_training_channels = np.asarray(arrays["final_support_training_channels"])
    _require(
        np.issubdtype(final_training_channels.dtype, np.integer)
        and np.array_equal(final_training_channels, np.arange(1, 17)),
        "final support training-channel axis changed",
    )
    for name in (
        "final_support_distinct_called_guides",
        "final_support_guides_with_three_untreated_cells",
        "final_support_untreated_carrier_cells",
        "final_support_training_channel_count",
    ):
        value = np.asarray(arrays[name])
        _require(
            value.shape == (expected_targets,)
            and np.issubdtype(value.dtype, np.integer)
            and np.all(value >= 0),
            f"{name} shape, dtype, or values changed",
        )
    final_supported = np.asarray(arrays["final_support_supported"])
    _require(
        final_supported.shape == (expected_targets,)
        and final_supported.dtype == np.bool_,
        "final support mask changed",
    )
    _require(
        np.array_equal(final_supported, np.asarray(arrays["evidence_support"])),
        "final support mask disagrees with estimator evidence support",
    )
    _require(
        np.all(
            ~np.asarray(arrays["fold_transferable"])
            | np.asarray(arrays["fold_support_supported"])
        )
        and np.all(~transferable | final_supported),
        "transferable rows exceed evidence-derived support",
    )
    fold_projection = np.asarray(arrays["fold_projection_diagonal"])
    _require(
        fold_projection.shape == (fold_count, expected_targets)
        and np.all(np.isfinite(fold_projection)),
        "fold projection diagonal changed",
    )
    unmodeled_indptr = np.asarray(arrays["fold_unmodeled_indptr"])
    unmodeled_values = np.asarray(arrays["fold_unmodeled_target_burden"])
    _require(
        unmodeled_indptr.shape == (fold_count + 1,)
        and np.issubdtype(unmodeled_indptr.dtype, np.integer)
        and unmodeled_indptr[0] == 0
        and np.all(np.diff(unmodeled_indptr) >= 0)
        and unmodeled_indptr[-1] == unmodeled_values.size,
        "packed fold unmodeled-target burden changed",
    )
    _require(
        unmodeled_values.ndim == 1
        and np.all(np.isfinite(unmodeled_values))
        and np.all(unmodeled_values >= 0),
        "fold unmodeled-target burden is invalid",
    )
    stability_pairs = np.asarray(arrays["fold_stability_pairs"])
    expected_pairs = np.asarray(
        [
            (left, right)
            for left in range(fold_count)
            for right in range(left + 1, fold_count)
        ],
        dtype=np.int64,
    )
    _require(
        np.issubdtype(stability_pairs.dtype, np.integer)
        and np.array_equal(stability_pairs, expected_pairs),
        "fold stability pairs changed",
    )
    for name in ("fold_effect_cosine", "fold_effect_sign_agreement"):
        value = np.asarray(arrays[name], dtype=np.float64)
        _require(
            value.shape == (len(expected_pairs), expected_targets)
            and not np.any(np.isinf(value)),
            f"{name} shape or values changed",
        )
    projection_report = manifest.get("projection_postconditions")
    _require(
        isinstance(projection_report, dict)
        and projection_report.get("validated_after_serialization") is True,
        "serialized projection postconditions are missing",
    )
    projection_errors = {
        "projection_nullspace_annihilation_error": "nullspace_annihilation_error",
        "projection_idempotence_error": "projection_idempotence_error",
        "projection_gauge_shift_invariance_error": "gauge_shift_invariance_error",
    }
    tolerance = float(projection_report.get("tolerance", math.nan))
    _require(
        math.isfinite(tolerance) and tolerance >= 0.0, "projection tolerance is invalid"
    )
    for array_name, report_name in projection_errors.items():
        value = float(np.asarray(arrays[array_name]).reshape(()))
        _require(
            math.isfinite(value)
            and value <= tolerance
            and projection_report.get(report_name) == value,
            f"{array_name} changed or exceeds tolerance",
        )
    solver = manifest.get("final_solver")
    _require(isinstance(solver, dict), "source final_solver is missing")
    _require(
        solver.get("all_components_converged") is True, "source solver did not converge"
    )
    _require(
        isinstance(solver.get("maximum_relative_residual"), (int, float))
        and math.isfinite(float(solver["maximum_relative_residual"]))
        and float(solver["maximum_relative_residual"])
        <= float(solver.get("relative_tolerance", -1)),
        "source solver residual exceeds tolerance",
    )


def _load_frozen_source_model(
    directory: Path,
    *,
    expected_targets: int,
    expected_responses: int,
    expected_exportable: int,
    verify_current_implementations: bool,
) -> FrozenSourceModel:
    directory = Path(directory)
    _require(
        directory.is_dir() and not directory.is_symlink(),
        "source directory is missing or unsafe",
    )
    manifest_path = directory / "manifest.json"
    _require(
        manifest_path.is_file() and not manifest_path.is_symlink(),
        "source manifest is missing or unsafe",
    )
    manifest = load_json(manifest_path)
    expected_manifest = {
        "status": "frozen_before_held_outcome_access",
        "experiment": "107-jurkat-high-moi",
        "primary_condition": "untreated",
        "held_truth_opened": False,
        "nadig_expression_or_scores_opened": False,
        "challenge_outcomes_opened": False,
        "model_file": "source-model.npz",
    }
    for key, value in expected_manifest.items():
        _require_exact(manifest, key, value)
    if verify_current_implementations:
        implementation_files = {
            "model_input_implementation_sha256": HERE / "model_input.py",
            "modeling_implementation_sha256": HERE / "modeling.py",
            "fit_implementation_sha256": HERE / "fit_source.py",
            "source_artifacts_implementation_sha256": HERE / "source_artifacts.py",
        }
        for key, implementation_path in implementation_files.items():
            _require(
                manifest.get(key) == sha256_file(implementation_path),
                f"source implementation changed: {implementation_path.name}",
            )
    model_path = directory / "source-model.npz"
    _require(
        model_path.is_file() and not model_path.is_symlink(),
        "source model is missing or unsafe",
    )
    model_sha = sha256_file(model_path)
    _require(manifest.get("model_sha256") == model_sha, "source model hash changed")
    try:
        with np.load(model_path, allow_pickle=False) as archive:
            arrays = {name: archive[name] for name in archive.files}
    except (OSError, ValueError, KeyError) as error:
        raise StrictValidationError(
            f"cannot load frozen source model: {error}"
        ) from error
    _validate_source_arrays(
        arrays,
        manifest,
        expected_targets=expected_targets,
        expected_responses=expected_responses,
        expected_exportable=expected_exportable,
    )
    return FrozenSourceModel(
        directory=directory.resolve(),
        manifest=manifest,
        manifest_sha256=sha256_file(manifest_path),
        model_sha256=model_sha,
        exported_effects=np.asarray(arrays["exported_effects"], dtype=np.float64),
        target_symbols=np.asarray(arrays["target_symbols"]).astype(str),
        response_symbols=np.asarray(arrays["response_symbols"]).astype(str),
        response_vcc_indices=np.asarray(arrays["response_vcc_indices"], dtype=np.int64),
        response_exportable=np.asarray(arrays["response_exportable"], dtype=bool),
        transferable=np.asarray(arrays["transferable"], dtype=bool),
    )


def load_frozen_source_model(directory: Path = SOURCE_DIRECTORY) -> FrozenSourceModel:
    """Load the production frozen source artifact, failing closed."""

    return _load_frozen_source_model(
        directory,
        expected_targets=18_465,
        expected_responses=374,
        expected_exportable=349,
        verify_current_implementations=True,
    )


def construct_gallery_targets(
    source_targets: Sequence[str],
    nadig_targets: Sequence[str],
    vcc_targets: Sequence[str],
    *,
    take_first: int = 1_800,
    gallery_size: int = 300,
) -> tuple[tuple[str, ...], ...]:
    """Construct SHA-sorted disjoint galleries without reading outcome values."""

    eligible = set(source_targets) & set(nadig_targets)
    eligible.difference_update(vcc_targets)
    ordered = sorted(
        eligible,
        key=lambda symbol: hashlib.sha256(f"exp107|{symbol}".encode("utf-8")).digest(),
    )
    _require(len(ordered) >= take_first, "not enough eligible targets")
    selected = tuple(ordered[:take_first])
    _require(take_first % gallery_size == 0, "gallery partition is not exact")
    return tuple(
        selected[start : start + gallery_size]
        for start in range(0, take_first, gallery_size)
    )


def score_axis(
    response_symbols: Sequence[str],
    response_exportable: Sequence[bool],
    nadig_genes: Sequence[str],
    vcc_targets: Sequence[str],
    gallery_targets: Sequence[str],
) -> tuple[str, ...]:
    responses = np.asarray(response_symbols).astype(str)
    exportable = np.asarray(response_exportable, dtype=bool)
    _require(
        responses.ndim == exportable.ndim == 1 and responses.shape == exportable.shape,
        "source response and exportable axes do not align",
    )
    allowed = set(nadig_genes) - set(vcc_targets) - set(gallery_targets)
    return tuple(sorted(set(responses[exportable]) & allowed))


def whitener_columns(
    official_genes: Sequence[str],
    donor_genes: Mapping[str, Sequence[str]],
    vcc_targets: Sequence[str],
    gallery_targets: Sequence[str],
    *,
    donor_order: Sequence[str] = DONOR_ORDER,
) -> np.ndarray:
    official_index = {symbol: index for index, symbol in enumerate(official_genes)}
    _require(
        len(official_index) == len(official_genes), "official gene axis has duplicates"
    )
    union: set[str] = set()
    for name in donor_order:
        _require(name in donor_genes, f"missing donor gene axis: {name}")
        union.update(set(donor_genes[name]) & set(official_index))
    union.difference_update(vcc_targets)
    union.difference_update(gallery_targets)
    return np.asarray(
        sorted(official_index[symbol] for symbol in union), dtype=np.int64
    )


def select_whitener_basis_rows(
    donor_perts: Mapping[str, Sequence[str]],
    essential_targets: set[str],
    gallery_targets: set[str],
    *,
    row_cap: int = 12_000,
    donor_order: Sequence[str] = DONOR_ORDER,
) -> dict[str, np.ndarray]:
    """Mirror production row sampling, then prune essential and held-gallery rows."""

    rng = np.random.default_rng(0)
    result: dict[str, np.ndarray] = {}
    for name in donor_order:
        _require(name in donor_perts, f"missing donor perturbation axis: {name}")
        perts = tuple(donor_perts[name])
        rows = (
            np.sort(rng.choice(len(perts), row_cap, replace=False))
            if len(perts) > row_cap
            else np.arange(len(perts), dtype=np.int64)
        )
        # Preserve the production common-essential duplicate policy: only the last
        # row for a duplicated perturbation receives its symbol.
        last = {symbol: index for index, symbol in enumerate(perts)}
        row_names: list[str | None] = [None] * len(perts)
        for symbol, index in last.items():
            row_names[index] = symbol
        keep = np.asarray(
            [
                row_names[index] not in essential_targets
                and perts[index] not in gallery_targets
                for index in rows
            ],
            dtype=bool,
        )
        result[name] = rows[keep]
    return result


def reconstruct_gallery_axes(
    source: FrozenSourceModel,
    *,
    nadig_targets: Sequence[str],
    nadig_genes: Sequence[str],
    official_genes: Sequence[str],
    vcc_targets: Sequence[str],
    donor_genes: Mapping[str, Sequence[str]],
    donor_perts: Mapping[str, Sequence[str]],
    essential_targets: set[str],
    contract: Mapping[str, Any],
) -> tuple[GalleryAxes, ...]:
    """Reconstruct and verify every frozen gallery, score axis, and whitener plan."""

    galleries = construct_gallery_targets(
        source.target_symbols.tolist(),
        nadig_targets,
        vcc_targets,
    )
    construction = contract["gallery_construction"]
    combined = tuple(symbol for gallery in galleries for symbol in gallery)
    _require(
        len(set(source.target_symbols.tolist()) & set(nadig_targets) - set(vcc_targets))
        == construction["eligible_target_count"],
        "eligible target count changed",
    )
    _require(
        sha256_strings(combined)
        == construction["combined_order_sha256_newline_terminated"],
        "combined gallery order changed",
    )
    ambiguous = contract["known_coverage"]["ambiguous_source_target"]
    _require(
        combined[int(ambiguous["combined_zero_based_index"])] == ambiguous["symbol"],
        "ambiguous source target position changed",
    )
    source_target_index = {
        symbol: index for index, symbol in enumerate(source.target_symbols)
    }
    source_response_index = {
        symbol: index for index, symbol in enumerate(source.response_symbols)
    }
    nadig_target_index = {symbol: index for index, symbol in enumerate(nadig_targets)}
    nadig_gene_index = {symbol: index for index, symbol in enumerate(nadig_genes)}
    official_index = {symbol: index for index, symbol in enumerate(official_genes)}

    plans: list[GalleryAxes] = []
    for number, targets in enumerate(galleries, start=1):
        frozen = construction["galleries"][number - 1]
        _require(
            sha256_strings(targets) == frozen["target_axis_sha256_newline_terminated"],
            f"gallery {number} target axis changed",
        )
        scores = score_axis(
            source.response_symbols.tolist(),
            source.response_exportable,
            nadig_genes,
            vcc_targets,
            targets,
        )
        _require(
            len(scores) == frozen["score_gene_count"],
            f"gallery {number} score count changed",
        )
        _require(
            sha256_strings(scores) == frozen["score_axis_sha256_newline_terminated"],
            f"gallery {number} score axis changed",
        )
        columns = whitener_columns(
            official_genes,
            donor_genes,
            vcc_targets,
            targets,
        )
        _require(
            len(columns) == frozen["whitener_column_count"],
            f"gallery {number} whitener column count changed",
        )
        _require(
            sha256_int64(columns) == frozen["whitener_columns_int64_sha256"],
            f"gallery {number} whitener columns changed",
        )
        rows = select_whitener_basis_rows(
            donor_perts,
            essential_targets,
            set(targets),
        )
        _require(
            sum(len(value) for value in rows.values()) == frozen["whitener_basis_rows"],
            f"gallery {number} whitener basis row count changed",
        )
        target_source_indices = np.asarray(
            [source_target_index[symbol] for symbol in targets], dtype=np.int64
        )
        transferable = source.transferable[target_source_indices].copy()
        if ambiguous["symbol"] in targets:
            transferable[targets.index(ambiguous["symbol"])] = False
        plans.append(
            GalleryAxes(
                number=number,
                targets=targets,
                target_indices_in_nadig=np.asarray(
                    [nadig_target_index[symbol] for symbol in targets], dtype=np.int64
                ),
                target_indices_in_source=target_source_indices,
                source_transferable=transferable,
                score_genes=scores,
                score_indices_in_nadig=np.asarray(
                    [nadig_gene_index[symbol] for symbol in scores], dtype=np.int64
                ),
                score_indices_in_source=np.asarray(
                    [source_response_index[symbol] for symbol in scores], dtype=np.int64
                ),
                score_indices_in_official=np.asarray(
                    [official_index[symbol] for symbol in scores], dtype=np.int64
                ),
                whitener_columns=columns,
                basis_rows=rows,
            )
        )
    return tuple(plans)


def pds_per_target(
    prediction: np.ndarray,
    truth: np.ndarray,
    supported: np.ndarray | None = None,
) -> np.ndarray:
    """Official midrank PDS, with unsupported prediction rows fixed at 0.5."""

    predicted = np.asarray(prediction, dtype=np.float64)
    observed = np.asarray(truth, dtype=np.float64)
    _require(
        predicted.ndim == 2 and predicted.shape == observed.shape,
        "PDS arrays do not align",
    )
    _require(predicted.shape[0] >= 2, "PDS needs at least two targets")
    _require(
        np.all(np.isfinite(predicted)) and np.all(np.isfinite(observed)),
        "PDS arrays are non-finite",
    )
    mask = (
        np.ones(predicted.shape[0], dtype=bool)
        if supported is None
        else np.asarray(supported, dtype=bool)
    )
    _require(mask.shape == (predicted.shape[0],), "PDS support mask does not align")
    result = np.full(predicted.shape[0], 0.5, dtype=np.float64)
    truth_norm = np.linalg.norm(observed, axis=1)
    truth_unit = np.divide(
        observed,
        truth_norm[:, None],
        out=np.zeros_like(observed),
        where=truth_norm[:, None] > 0,
    )
    rows = np.flatnonzero(mask)
    pred_norm = np.linalg.norm(predicted[rows], axis=1)
    nonzero = pred_norm > 0
    rows = rows[nonzero]
    if rows.size:
        pred_unit = predicted[rows] / np.linalg.norm(predicted[rows], axis=1)[:, None]
        similarity = pred_unit @ truth_unit.T
        own = similarity[np.arange(len(rows)), rows][:, None]
        rank = (similarity > own).sum(axis=1).astype(np.float64)
        rank += 0.5 * ((similarity == own).sum(axis=1) - 1)
        result[rows] = 1.0 - rank / (predicted.shape[0] - 1)
    return result


def _fraction_key(value: float) -> float:
    number = float(value)
    _require(math.isfinite(number), "blend fraction is non-finite")
    return number


def nested_blend_selection(
    pds_by_fraction: Mapping[float, Sequence[np.ndarray]],
    fractions: Sequence[float] = FRACTIONS,
) -> NestedSelection:
    """Perform six leave-one-gallery-out selections with the frozen tie-break."""

    grid = tuple(sorted(_fraction_key(value) for value in fractions))
    _require(grid and grid[0] == 0.0, "fraction grid must include zero")
    _require(set(map(float, pds_by_fraction)) == set(grid), "PDS fraction grid changed")
    gallery_count = len(pds_by_fraction[grid[0]])
    _require(gallery_count == 6, "strict nested selection requires six galleries")
    means = np.empty((len(grid), gallery_count), dtype=np.float64)
    for fraction_index, fraction in enumerate(grid):
        galleries = pds_by_fraction[fraction]
        _require(len(galleries) == gallery_count, "fraction gallery counts disagree")
        for gallery_index, values in enumerate(galleries):
            array = np.asarray(values, dtype=np.float64)
            _require(
                array.ndim == 1 and array.size > 0,
                "gallery PDS is empty or not one-dimensional",
            )
            _require(np.all(np.isfinite(array)), "gallery PDS is non-finite")
            means[fraction_index, gallery_index] = array.mean()
    selected: list[float] = []
    held_means = np.empty(gallery_count, dtype=np.float64)
    for held in range(gallery_count):
        training = np.delete(means, held, axis=1).mean(axis=1)
        best = np.max(training)
        chosen_index = int(np.flatnonzero(training == best)[0])
        selected.append(grid[chosen_index])
        held_means[held] = means[chosen_index, held]
    baseline = means[grid.index(0.0)]
    gains = held_means - baseline
    return NestedSelection(
        fractions_by_held_gallery=tuple(selected),
        baseline_gallery_means=baseline.copy(),
        selected_gallery_means=held_means,
        held_gallery_gains=gains,
        macro_baseline=float(math.fsum(baseline.tolist()) / gallery_count),
        macro_selected=float(math.fsum(held_means.tolist()) / gallery_count),
        macro_gain=float(math.fsum(gains.tolist()) / gallery_count),
    )


def pds_fraction_grid(
    incumbents: Sequence[np.ndarray],
    sources: Sequence[np.ndarray],
    truths: Sequence[np.ndarray],
    transferable: Sequence[np.ndarray],
    *,
    fractions: Sequence[float] = FRACTIONS,
) -> dict[float, tuple[np.ndarray, ...]]:
    """Build complete-gallery PDS arrays for the frozen blend grid."""

    _require(
        len(incumbents) == len(sources) == len(truths) == len(transferable) == 6,
        "strict PDS grid requires six galleries",
    )
    result: dict[float, tuple[np.ndarray, ...]] = {}
    for fraction in sorted(map(float, fractions)):
        gallery_scores: list[np.ndarray] = []
        for base, source, truth, support in zip(
            incumbents, sources, truths, transferable, strict=True
        ):
            candidate, _ = blend_with_source(base, source, fraction, support)
            gallery_scores.append(pds_per_target(candidate, truth))
        result[fraction] = tuple(gallery_scores)
    return result


def final_fraction(
    pds_by_fraction: Mapping[float, Sequence[np.ndarray]],
    fractions: Sequence[float] = FRACTIONS,
) -> float:
    grid = tuple(sorted(map(float, fractions)))
    objectives = []
    for fraction in grid:
        galleries = pds_by_fraction[fraction]
        objectives.append(
            float(
                np.mean(
                    [
                        np.asarray(values, dtype=np.float64).mean()
                        for values in galleries
                    ]
                )
            )
        )
    best = max(objectives)
    return grid[objectives.index(best)]


def permute_transferable_rows(
    source: np.ndarray,
    transferable: np.ndarray,
    *,
    permutation_zero_based: int,
    gallery_one_based: int,
) -> np.ndarray:
    """Permute labels among transferable source rows with the frozen SeedSequence."""

    _require(
        0 <= permutation_zero_based < PERMUTATION_REPLICATES,
        "permutation index is outside 0..98",
    )
    _require(1 <= gallery_one_based <= 6, "gallery number is outside 1..6")
    values = np.asarray(source)
    mask = np.asarray(transferable, dtype=bool)
    _require(
        values.ndim == 2 and mask.shape == (values.shape[0],),
        "permutation arrays do not align",
    )
    output = np.zeros_like(values)
    rows = np.flatnonzero(mask)
    rng = np.random.Generator(
        np.random.PCG64(
            np.random.SeedSequence(
                [1_072_026, 107, permutation_zero_based, gallery_one_based]
            )
        )
    )
    output[rows] = values[rng.permutation(rows)]
    return output


def deterministic_label_permutations(
    source: np.ndarray,
    transferable: np.ndarray,
    *,
    gallery_one_based: int,
) -> tuple[np.ndarray, ...]:
    """Return exactly the 99 frozen within-gallery label permutations."""

    return tuple(
        permute_transferable_rows(
            source,
            transferable,
            permutation_zero_based=replicate,
            gallery_one_based=gallery_one_based,
        )
        for replicate in range(PERMUTATION_REPLICATES)
    )


def blend_with_source(
    incumbent: np.ndarray,
    source: np.ndarray,
    fraction: float,
    transferable: np.ndarray,
    coordinate_supported: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Norm-match and blend only supported cells, preserving exact fallbacks."""

    base = np.asarray(incumbent)
    added = np.asarray(source, dtype=np.float64)
    _require(base.ndim == 2 and added.shape == base.shape, "blend arrays do not align")
    _require(np.all(np.isfinite(base)), "incumbent contains non-finite values")
    rows = np.asarray(transferable, dtype=bool)
    columns = (
        np.ones(base.shape[1], dtype=bool)
        if coordinate_supported is None
        else np.asarray(coordinate_supported, dtype=bool)
    )
    _require(rows.shape == (base.shape[0],), "blend target support does not align")
    _require(
        columns.shape == (base.shape[1],), "blend coordinate support does not align"
    )
    fraction = float(fraction)
    _require(fraction in FRACTIONS, "blend fraction is outside the frozen grid")
    candidate = base.copy()
    effective = rows.copy()
    if fraction == 0.0:
        return candidate, np.zeros_like(rows)
    for row in np.flatnonzero(rows):
        source_row = added[row, columns]
        incumbent_row = np.asarray(base[row, columns], dtype=np.float64)
        source_norm = float(np.linalg.norm(source_row))
        incumbent_norm = float(np.linalg.norm(incumbent_row))
        if (
            not math.isfinite(source_norm)
            or source_norm <= 0
            or not math.isfinite(incumbent_norm)
        ):
            effective[row] = False
            continue
        scaled = source_row * (incumbent_norm / source_norm)
        candidate[row, columns] = (
            (1.0 - fraction) * incumbent_row + fraction * scaled
        ).astype(base.dtype, copy=False)
    assert_exact_fallback(base, candidate, effective, columns)
    return candidate, effective


def assert_operator_missing_score_fallback(
    *,
    genes: np.ndarray,
    score_indices: np.ndarray,
    operator_columns: np.ndarray,
    incumbent: np.ndarray,
    source: np.ndarray,
    candidates: Mapping[float, np.ndarray],
    parameters: Mapping[str, Any],
    label: str,
) -> tuple[str, ...]:
    """Prove exact fallback for scored genes absent from the donor operator."""

    gene_axis = np.asarray(genes).astype(str)
    scores = np.asarray(score_indices, dtype=np.int64)
    columns = np.asarray(operator_columns, dtype=np.int64)
    base = np.asarray(incumbent, dtype=np.float64)
    added = np.asarray(source, dtype=np.float64)
    _require(
        gene_axis.ndim == scores.ndim == columns.ndim == 1,
        f"{label} operator-fallback axes are invalid",
    )
    _require(
        base.ndim == 2 and added.shape == base.shape and base.shape[1] == len(gene_axis),
        f"{label} operator-fallback arrays do not align",
    )
    missing = scores[~np.isin(scores, columns, assume_unique=False)]
    missing_genes = tuple(gene_axis[missing].tolist())
    declared = parameters.get("operator_missing_score_genes")
    if missing_genes:
        _require(
            declared == list(missing_genes),
            f"{label} operator-missing score-gene declaration changed",
        )
        _require(
            parameters.get("operator_missing_score_policy")
            == "retain_in_score_axis_require_zero_incumbent_and_source",
            f"{label} operator-missing score policy changed",
        )
        incumbent_missing = np.ascontiguousarray(base[:, missing])
        source_missing = np.ascontiguousarray(added[:, missing])
        _require(
            np.array_equal(
                incumbent_missing.view(np.uint8),
                np.zeros_like(incumbent_missing).view(np.uint8),
            ),
            f"{label} operator-missing incumbent coordinate is not positive zero",
        )
        _require(
            np.array_equal(
                source_missing.view(np.uint8),
                np.zeros_like(source_missing).view(np.uint8),
            ),
            f"{label} operator-missing source coordinate is not positive zero",
        )
        for fraction, candidate in candidates.items():
            values = np.asarray(candidate, dtype=np.float64)
            _require(
                values.shape == base.shape,
                f"{label} candidate shape changed for fraction {fraction}",
            )
            _require(
                np.array_equal(
                    np.ascontiguousarray(values[:, missing]).view(np.uint8),
                    incumbent_missing.view(np.uint8),
                ),
                f"{label} operator-missing candidate changed for fraction {fraction}",
            )
    else:
        _require(
            declared in (None, []),
            f"{label} declares nonexistent operator-missing score genes",
        )
    return missing_genes


def permutation_null_distribution(
    incumbents: Sequence[np.ndarray],
    sources: Sequence[np.ndarray],
    truths: Sequence[np.ndarray],
    transferable: Sequence[np.ndarray],
    *,
    repetitions: int = PERMUTATION_REPLICATES,
    fractions: Sequence[float] = FRACTIONS,
) -> np.ndarray:
    """Repeat complete nested selection for all 99 frozen label permutations."""

    _require(
        repetitions == PERMUTATION_REPLICATES,
        "strict permutation null requires exactly 99 replicates",
    )
    _require(
        len(incumbents) == len(sources) == len(truths) == len(transferable) == 6,
        "strict permutation null requires six galleries",
    )
    gains = np.empty(repetitions, dtype=np.float64)
    for replicate in range(repetitions):
        permuted = [
            permute_transferable_rows(
                source,
                support,
                permutation_zero_based=replicate,
                gallery_one_based=gallery + 1,
            )
            for gallery, (source, support) in enumerate(
                zip(sources, transferable, strict=True)
            )
        ]
        grid = pds_fraction_grid(
            incumbents,
            permuted,
            truths,
            transferable,
            fractions=fractions,
        )
        gains[replicate] = nested_blend_selection(grid, fractions).macro_gain
    return gains


def assert_exact_fallback(
    incumbent: np.ndarray,
    candidate: np.ndarray,
    transferable: np.ndarray,
    coordinate_supported: np.ndarray | None = None,
) -> None:
    base = np.asarray(incumbent)
    changed = np.asarray(candidate)
    _require(
        base.shape == changed.shape and base.dtype == changed.dtype,
        "fallback arrays changed shape or dtype",
    )
    rows = np.asarray(transferable, dtype=bool)
    columns = (
        np.ones(base.shape[1], dtype=bool)
        if coordinate_supported is None
        else np.asarray(coordinate_supported, dtype=bool)
    )
    _require(
        rows.shape == (base.shape[0],) and columns.shape == (base.shape[1],),
        "fallback masks do not align",
    )
    fallback = ~(rows[:, None] & columns[None, :])
    left = np.ascontiguousarray(base[fallback]).view(np.uint8)
    right = np.ascontiguousarray(changed[fallback]).view(np.uint8)
    _require(np.array_equal(left, right), "unsupported target or coordinate changed")


def verify_exact_fallback_evidence(
    evidence: Mapping[str, Any],
) -> dict[str, Any]:
    """Verify six galleries directly or from a hash-pinned NPZ evidence artifact."""

    if set(evidence) == {"artifact_path", "artifact_sha256"}:
        artifact_path = Path(str(evidence["artifact_path"]))
        expected_sha256 = evidence["artifact_sha256"]
        _require(_is_sha256(expected_sha256), "fallback artifact hash is invalid")
        _require(
            artifact_path.is_file() and not artifact_path.is_symlink(),
            "fallback evidence artifact is missing or unsafe",
        )
        observed_sha256 = sha256_file(artifact_path)
        _require(
            observed_sha256 == expected_sha256,
            "fallback evidence artifact hash changed",
        )
        expected_keys = {
            f"gallery_{gallery}_{name}"
            for gallery in range(1, 7)
            for name in (
                "incumbent",
                "candidate",
                "transferable",
                "coordinate_supported",
            )
        }
        try:
            with np.load(artifact_path, allow_pickle=False) as archive:
                _require(
                    set(archive.files) == expected_keys,
                    "fallback evidence artifact schema changed",
                )
                direct = {
                    "incumbents": tuple(
                        archive[f"gallery_{gallery}_incumbent"]
                        for gallery in range(1, 7)
                    ),
                    "candidates": tuple(
                        archive[f"gallery_{gallery}_candidate"]
                        for gallery in range(1, 7)
                    ),
                    "transferable": tuple(
                        archive[f"gallery_{gallery}_transferable"]
                        for gallery in range(1, 7)
                    ),
                    "coordinate_supported": tuple(
                        archive[f"gallery_{gallery}_coordinate_supported"]
                        for gallery in range(1, 7)
                    ),
                }
                result = verify_exact_fallback_evidence(direct)
        except (OSError, ValueError, KeyError) as error:
            raise StrictValidationError(
                f"cannot verify fallback evidence artifact: {error}"
            ) from error
        return {
            **result,
            "artifact_sha256": observed_sha256,
        }

    _require(
        set(evidence)
        == {"incumbents", "candidates", "transferable", "coordinate_supported"},
        "fallback evidence fields changed",
    )
    incumbents = evidence["incumbents"]
    candidates = evidence["candidates"]
    transferables = evidence["transferable"]
    coordinate_support = evidence["coordinate_supported"]
    _require(
        all(
            isinstance(values, Sequence)
            for values in (incumbents, candidates, transferables, coordinate_support)
        )
        and len(incumbents)
        == len(candidates)
        == len(transferables)
        == len(coordinate_support)
        == 6,
        "fallback evidence requires six galleries",
    )
    galleries: list[dict[str, Any]] = []
    for gallery, (incumbent, candidate, rows, columns) in enumerate(
        zip(
            incumbents,
            candidates,
            transferables,
            coordinate_support,
            strict=True,
        ),
        start=1,
    ):
        base = np.asarray(incumbent)
        changed = np.asarray(candidate)
        row_mask = np.asarray(rows)
        column_mask = np.asarray(columns)
        _require(
            base.ndim == 2
            and changed.ndim == 2
            and np.all(np.isfinite(base))
            and np.all(np.isfinite(changed)),
            "fallback arrays must be finite matrices",
        )
        _require(row_mask.dtype == np.bool_, "fallback target mask is not boolean")
        _require(
            column_mask.dtype == np.bool_,
            "fallback coordinate mask is not boolean",
        )
        assert_exact_fallback(base, changed, row_mask, column_mask)
        galleries.append(
            {
                "gallery": gallery,
                "incumbent_sha256": sha256_array(base),
                "candidate_sha256": sha256_array(changed),
                "transferable_sha256": sha256_array(row_mask),
                "coordinate_supported_sha256": sha256_array(column_mask),
            }
        )
    payload = {
        "schema_version": 1,
        "gallery_count": 6,
        "galleries": galleries,
    }
    payload["binding_sha256"] = sha256_json_payload(payload)
    return payload


def apply_supported_principal_whitener(
    source_supported: np.ndarray,
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
    supported_operator_positions: np.ndarray,
    *,
    gamma: float = 0.65,
    lam: float = 0.0,
) -> np.ndarray:
    """Apply the supported-by-supported principal block of a frozen whitener."""

    source = np.asarray(source_supported, dtype=np.float64)
    values = np.asarray(eigenvalues, dtype=np.float64)
    vectors = np.asarray(eigenvectors, dtype=np.float64)
    positions = np.asarray(supported_operator_positions, dtype=np.int64)
    _require(
        vectors.shape == (values.size, values.size),
        "whitener eigensystem does not align",
    )
    _require(
        source.ndim == 2 and source.shape[1] == positions.size,
        "source/operator support does not align",
    )
    _require(
        np.all((0 <= positions) & (positions < values.size)),
        "operator positions are out of range",
    )
    _require(
        np.all(np.isfinite(values)) and np.all(values > 0),
        "whitener eigenvalues are invalid",
    )
    weights = (values + float(lam)) ** (-float(gamma) / 2.0)
    selected_vectors = vectors[positions, :]
    principal = (selected_vectors * weights) @ selected_vectors.T
    return source @ principal


def reciprocal_knn_clusters(
    incumbent: np.ndarray,
    *,
    gallery_one_based: int,
    k: int = 5,
    minimum_cosine: float = 0.80,
) -> np.ndarray:
    values = np.asarray(incumbent, dtype=np.float64)
    _require(values.ndim == 2 and values.shape[0] >= 2, "cluster input is invalid")
    _require(np.all(np.isfinite(values)), "cluster input is non-finite")
    norms = np.linalg.norm(values, axis=1)
    unit = np.divide(
        values, norms[:, None], out=np.zeros_like(values), where=norms[:, None] > 0
    )
    similarity = unit @ unit.T
    np.fill_diagonal(similarity, -np.inf)
    count = min(max(1, int(k)), len(values) - 1)
    neighbors = np.argsort(-similarity, axis=1, kind="stable")[:, :count]
    directed = np.zeros((len(values), len(values)), dtype=bool)
    directed[np.arange(len(values))[:, None], neighbors] = True
    reciprocal = directed & directed.T & (similarity >= minimum_cosine)
    parent = np.arange(len(values), dtype=np.int64)

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = int(parent[index])
        return index

    for left, right in zip(*np.where(np.triu(reciprocal, k=1)), strict=True):
        left_root, right_root = find(int(left)), find(int(right))
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)
    roots = [find(index) for index in range(len(values))]
    labels = {root: number for number, root in enumerate(sorted(set(roots)))}
    return np.asarray(
        [f"g{gallery_one_based}:c{labels[root]}" for root in roots],
        dtype="U32",
    )


def _cluster_resample_weights(
    labels: np.ndarray,
    *,
    replicate_zero_based: int,
    gallery_one_based: int,
) -> np.ndarray:
    labels = np.asarray(labels).astype(str)
    unique = np.unique(labels)
    _require(unique.size > 0, "bootstrap gallery has no clusters")
    rng = np.random.Generator(
        np.random.PCG64(
            np.random.SeedSequence(
                [1_072_026, 107, 20_000, replicate_zero_based, gallery_one_based]
            )
        )
    )
    sampled = rng.choice(unique, size=unique.size, replace=True)
    multiplicity = {label: int(np.count_nonzero(sampled == label)) for label in unique}
    return np.asarray([multiplicity[label] for label in labels], dtype=np.int64)


def gallery_stratified_cluster_bootstrap(
    pds_by_fraction: Mapping[float, Sequence[np.ndarray]],
    cluster_labels: Sequence[np.ndarray],
    *,
    repetitions: int = BOOTSTRAP_REPLICATES,
    fractions: Sequence[float] = FRACTIONS,
) -> dict[str, Any]:
    """Repeat nested selection inside deterministic whole-cluster resamples."""

    grid = tuple(sorted(map(float, fractions)))
    gallery_count = len(cluster_labels)
    _require(gallery_count == 6, "strict bootstrap requires six galleries")
    _require(
        0 < repetitions <= BOOTSTRAP_REPLICATES, "invalid bootstrap repetition count"
    )
    arrays: dict[float, tuple[np.ndarray, ...]] = {}
    for fraction in grid:
        _require(fraction in pds_by_fraction, "bootstrap fraction is missing")
        values = tuple(
            np.asarray(item, dtype=np.float64) for item in pds_by_fraction[fraction]
        )
        _require(len(values) == gallery_count, "bootstrap gallery count changed")
        arrays[fraction] = values
    labels = tuple(np.asarray(item).astype(str) for item in cluster_labels)
    for gallery in range(gallery_count):
        size = arrays[grid[0]][gallery].size
        _require(
            size > 0 and labels[gallery].shape == (size,),
            "bootstrap labels do not align",
        )
        _require(
            all(label.startswith(f"g{gallery + 1}:") for label in labels[gallery]),
            "bootstrap cluster identities are not gallery-prefixed",
        )
        _require(
            all(arrays[fraction][gallery].shape == (size,) for fraction in grid),
            "bootstrap PDS axes do not align",
        )
    draws = np.empty(repetitions, dtype=np.float64)
    for replicate in range(repetitions):
        means = np.empty((len(grid), gallery_count), dtype=np.float64)
        for gallery in range(gallery_count):
            weights = _cluster_resample_weights(
                labels[gallery],
                replicate_zero_based=replicate,
                gallery_one_based=gallery + 1,
            )
            denominator = int(weights.sum())
            _require(denominator > 0, "bootstrap draw has no targets")
            for fraction_index, fraction in enumerate(grid):
                means[fraction_index, gallery] = float(
                    np.dot(arrays[fraction][gallery], weights) / denominator
                )
        selected_held = np.empty(gallery_count, dtype=np.float64)
        for held in range(gallery_count):
            training = np.delete(means, held, axis=1).mean(axis=1)
            chosen = int(np.flatnonzero(training == np.max(training))[0])
            selected_held[held] = means[chosen, held]
        draws[replicate] = float((selected_held - means[0]).mean())
    observed = nested_blend_selection(pds_by_fraction, grid)
    return {
        "repetitions": repetitions,
        "point_estimate": observed.macro_gain,
        "lower_95": float(np.quantile(draws, 0.025, method="linear")),
        "upper_95": float(np.quantile(draws, 0.975, method="linear")),
        "draws": draws,
    }


def gain_concentration(
    per_target_deltas: Sequence[np.ndarray],
    cluster_labels: Sequence[np.ndarray],
    observed_nested_macro_gain: float,
) -> dict[str, Any]:
    _require(
        len(per_target_deltas) == len(cluster_labels) == 6,
        "concentration requires six galleries",
    )
    gain = float(observed_nested_macro_gain)
    if not math.isfinite(gain) or gain <= 0:
        return {
            "passes": False,
            "fraction": math.inf,
            "maximum_positive_cluster_sum": None,
        }
    contributions: dict[str, float] = {}
    target_count = 0
    for gallery, (values, labels) in enumerate(
        zip(per_target_deltas, cluster_labels, strict=True),
        start=1,
    ):
        delta = np.asarray(values, dtype=np.float64)
        cluster = np.asarray(labels).astype(str)
        _require(
            delta.shape == cluster.shape and delta.ndim == 1,
            "concentration arrays do not align",
        )
        _require(np.all(np.isfinite(delta)), "concentration deltas are non-finite")
        _require(
            all(label.startswith(f"g{gallery}:") for label in cluster),
            "cluster identities are not gallery-prefixed",
        )
        target_count += delta.size
        for label in np.unique(cluster):
            contributions[label] = contributions.get(label, 0.0) + float(
                delta[cluster == label].sum()
            )
    _require(target_count == 1_800, "concentration denominator requires 1,800 targets")
    maximum = max((value for value in contributions.values() if value > 0), default=0.0)
    fraction = maximum / (1_800.0 * gain)
    return {
        "passes": bool(fraction <= 0.5),
        "fraction": float(fraction),
        "maximum_positive_cluster_sum": float(maximum),
    }


def _validate_nested_selection(nested: NestedSelection) -> None:
    fractions = tuple(map(float, nested.fractions_by_held_gallery))
    baseline = np.asarray(nested.baseline_gallery_means, dtype=np.float64)
    selected = np.asarray(nested.selected_gallery_means, dtype=np.float64)
    gains = np.asarray(nested.held_gallery_gains, dtype=np.float64)
    _require(
        len(fractions) == 6 and all(value in FRACTIONS for value in fractions),
        "nested selection must contain six frozen fractions",
    )
    _require(
        baseline.shape == selected.shape == gains.shape == (6,),
        "nested selection must contain six gallery values",
    )
    _require(
        np.all(np.isfinite(baseline))
        and np.all(np.isfinite(selected))
        and np.all(np.isfinite(gains)),
        "nested selection contains non-finite values",
    )
    _require(
        np.array_equal(gains, selected - baseline),
        "nested gallery gains disagree with selected minus baseline",
    )
    expected_baseline = float(math.fsum(baseline.tolist()) / 6)
    expected_selected = float(math.fsum(selected.tolist()) / 6)
    expected_gain = float(math.fsum(gains.tolist()) / 6)
    _require(
        all(
            math.isfinite(value)
            for value in (
                nested.macro_baseline,
                nested.macro_selected,
                nested.macro_gain,
            )
        )
        and math.isclose(
            nested.macro_baseline,
            expected_baseline,
            rel_tol=1e-12,
            abs_tol=1e-15,
        )
        and math.isclose(
            nested.macro_selected,
            expected_selected,
            rel_tol=1e-12,
            abs_tol=1e-15,
        )
        and math.isclose(
            nested.macro_gain,
            expected_gain,
            rel_tol=1e-12,
            abs_tol=1e-15,
        )
        and math.isclose(
            nested.macro_gain,
            nested.macro_selected - nested.macro_baseline,
            rel_tol=1e-12,
            abs_tol=1e-15,
        ),
        "nested macro summaries are internally inconsistent",
    )


def _validate_held_context_metrics(
    metrics: Mapping[str, Any],
    held_gate: Mapping[str, Any],
) -> tuple[dict[str, bool], dict[str, dict[str, float]]]:
    required = tuple(held_gate["authorizing_directional_surfaces"])
    _require(
        required == AUTHORIZING_HELD_SURFACES,
        "authorizing held-surface set changed",
    )
    _require(
        isinstance(metrics, Mapping) and set(metrics) == set(required),
        "H1 and HEK293T held-context gates are both required",
    )
    pds_minimum = float(held_gate["pds_delta_minimum_each_surface"])
    directional_minimum = float(
        held_gate["directional_credit_delta_minimum_each_authorizing_surface"]
    )
    checks: dict[str, bool] = {}
    normalized: dict[str, dict[str, float]] = {}
    for surface in required:
        values = metrics[surface]
        _require(
            isinstance(values, Mapping)
            and set(values) == {"pds_delta", "directional_credit_delta"},
            f"{surface} held-context metrics are incomplete",
        )
        pds_delta = float(values["pds_delta"])
        directional_delta = float(values["directional_credit_delta"])
        _require(
            math.isfinite(pds_delta) and math.isfinite(directional_delta),
            f"{surface} held-context metrics are non-finite",
        )
        checks[f"{surface}_pds_delta_at_least_{pds_minimum}"] = pds_delta >= pds_minimum
        checks[f"{surface}_directional_credit_delta_at_least_{directional_minimum}"] = (
            directional_delta >= directional_minimum
        )
        normalized[surface] = {
            "pds_delta": pds_delta,
            "directional_credit_delta": directional_delta,
        }
    return checks, normalized


def promotion_gate_report(
    pure_source_gallery_pds: Sequence[float],
    nested: NestedSelection,
    bootstrap: Mapping[str, Any],
    permutation_null_gains: np.ndarray,
    concentration: Mapping[str, Any],
    held_context_metrics: Mapping[str, Any],
    *,
    fallback_evidence: Mapping[str, Any],
    contract: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Apply Nadig development gates plus H1/HEK293T promotion gates."""

    frozen = validate_contract() if contract is None else dict(contract)
    gates = frozen["promotion_gates"]
    held_gate = frozen["held_context_gate"]
    _validate_nested_selection(nested)
    pure = np.asarray(pure_source_gallery_pds, dtype=np.float64)
    null = np.asarray(permutation_null_gains, dtype=np.float64)
    _require(
        pure.shape == (6,) and np.all(np.isfinite(pure)),
        "pure-source PDS must contain six finite values",
    )
    _require(
        null.shape == (PERMUTATION_REPLICATES,) and np.all(np.isfinite(null)),
        "permutation null must contain exactly 99 finite gains",
    )
    lower = bootstrap.get("lower_95")
    _require(
        isinstance(lower, (int, float)) and math.isfinite(float(lower)),
        "bootstrap lower bound is invalid",
    )
    _require(
        bootstrap.get("repetitions") == BOOTSTRAP_REPLICATES,
        "promotion requires exactly 20,000 bootstrap replicates",
    )
    draws = np.asarray(bootstrap.get("draws"), dtype=np.float64)
    _require(
        draws.shape == (BOOTSTRAP_REPLICATES,) and np.all(np.isfinite(draws)),
        "promotion bootstrap draws must contain exactly 20,000 finite values",
    )
    expected_lower = float(np.quantile(draws, 0.025, method="linear"))
    _require(
        float(lower) == expected_lower,
        "bootstrap lower bound disagrees with its bound draws",
    )
    _require(
        isinstance(bootstrap.get("point_estimate"), (int, float))
        and float(bootstrap["point_estimate"]) == nested.macro_gain,
        "bootstrap point estimate disagrees with nested gain",
    )
    concentration_fraction = concentration.get("fraction")
    _require(
        isinstance(concentration_fraction, (int, float))
        and math.isfinite(float(concentration_fraction))
        and float(concentration_fraction) >= 0.0,
        "concentration fraction is invalid",
    )
    fallback = verify_exact_fallback_evidence(fallback_evidence)
    held_checks, normalized_held = _validate_held_context_metrics(
        held_context_metrics,
        held_gate,
    )
    pure_macro = float(math.fsum(pure.tolist()) / pure.size)
    checks = {
        "pure_source_each_gallery_minimum": bool(
            np.all(pure >= float(gates["pure_source_pds_each_gallery_minimum"]))
        ),
        "pure_source_macro_minimum": bool(
            pure_macro >= float(gates["pure_source_macro_pds_minimum"])
        ),
        "nested_gain_each_gallery_minimum": bool(
            np.all(
                nested.held_gallery_gains
                >= float(gates["nested_blend_pds_delta_each_gallery_minimum"])
            )
        ),
        "nested_macro_gain_minimum": bool(
            nested.macro_gain >= float(gates["nested_blend_macro_pds_delta_minimum"])
        ),
        "bootstrap_lower_95_strictly_above_0": bool(float(lower) > 0.0),
        "observed_gain_strictly_exceeds_permutation_maximum": bool(
            nested.macro_gain > float(np.max(null))
        ),
        "maximum_cluster_fraction_at_most_0_5": bool(
            math.isfinite(float(concentration_fraction))
            and float(concentration_fraction)
            <= float(gates["maximum_gain_fraction_from_one_target_cluster"])
        ),
        "exact_fallback": fallback["gallery_count"] == 6,
        **held_checks,
    }
    return {
        "passes": bool(all(checks.values())),
        "checks": checks,
        "pure_source_macro_pds": pure_macro,
        "nested_macro_gain": nested.macro_gain,
        "permutation_maximum_gain": float(np.max(null)),
        "bootstrap_lower_95": float(lower),
        "concentration_fraction": float(concentration_fraction),
        "fallback_evidence": fallback,
        "held_context_metrics": normalized_held,
        "benchmark_role": frozen["benchmark_role"],
    }


def _plan_payload(
    source: FrozenSourceModel,
    galleries: Sequence[GalleryAxes],
    *,
    nadig_targets: Sequence[str],
    nadig_genes: Sequence[str],
    nadig_artifacts: Mapping[str, Any],
    artifact_bindings: Mapping[str, Any],
) -> dict[str, Any]:
    validation_implementation = file_binding(Path(__file__), root=ROOT)
    return {
        "status": "pretruth_complete",
        "experiment": "107-jurkat-high-moi",
        "benchmark_role": "reused_development_benchmark_not_independent_confirmation",
        "contract_sha256": EXPECTED_CONTRACT_SHA256,
        "validation_implementation": validation_implementation,
        "source_manifest_sha256": source.manifest_sha256,
        "source_model_sha256": source.model_sha256,
        "source_directory": source.directory.relative_to(ROOT.resolve()).as_posix(),
        "artifact_bindings": artifact_bindings,
        "nadig_outcome_opened": False,
        "nadig_target_count": len(nadig_targets),
        "nadig_target_axis_sha256": sha256_strings(nadig_targets),
        "nadig_gene_count": len(nadig_genes),
        "nadig_gene_axis_sha256": sha256_strings(nadig_genes),
        "nadig_artifacts": nadig_artifacts,
        "galleries": [
            {
                "gallery": item.number,
                "targets": list(item.targets),
                "target_indices_in_nadig": item.target_indices_in_nadig.tolist(),
                "target_indices_in_source": item.target_indices_in_source.tolist(),
                "source_transferable": item.source_transferable.tolist(),
                "score_genes": list(item.score_genes),
                "score_indices_in_nadig": item.score_indices_in_nadig.tolist(),
                "score_indices_in_source": item.score_indices_in_source.tolist(),
                "score_indices_in_official": item.score_indices_in_official.tolist(),
                "whitener_columns": item.whitener_columns.tolist(),
                "basis_rows": {
                    name: item.basis_rows[name].tolist() for name in DONOR_ORDER
                },
            }
            for item in galleries
        ],
    }


def write_immutable_json(path: Path, payload: Mapping[str, Any]) -> None:
    """Create one JSON record exactly once and make it read-only."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    _require(not path.is_symlink(), f"immutable JSON path is a symlink: {path}")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    try:
        descriptor = os.open(path, flags, 0o444)
    except FileExistsError as error:
        raise StrictValidationError(
            f"immutable JSON already exists; replay refused: {path}"
        ) from error
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(path, 0o444)
        directory_descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    except BaseException:
        path.unlink(missing_ok=True)
        raise


def _require_production_path(observed: Path, expected: Path, label: str) -> None:
    _require(
        Path(observed).resolve() == Path(expected).resolve(),
        f"{label} is not the frozen production path",
    )


def bind_nadig_artifacts(
    outcome_directory: Path = NADIG_DIRECTORY,
) -> dict[str, Any]:
    """Hash the Nadig axes and opaque delta bytes without loading the array."""

    directory = Path(outcome_directory)
    _require_production_path(directory, NADIG_DIRECTORY, "Nadig directory")
    artifacts = {
        "target_axis": file_binding(directory / "perts.txt", root=ROOT),
        "gene_axis": file_binding(directory / "genes.txt", root=ROOT),
        "delta": file_binding(directory / "delta.npy", root=ROOT),
    }
    payload = {
        "array_values_opened": False,
        "artifacts": artifacts,
    }
    payload["binding_sha256"] = sha256_json_payload(payload)
    return payload


def verify_nadig_artifacts(
    expected: Mapping[str, Any],
    outcome_directory: Path = NADIG_DIRECTORY,
) -> dict[str, Any]:
    observed = bind_nadig_artifacts(outcome_directory)
    _require(
        observed == expected,
        "Nadig target, gene, or delta artifact changed after pretruth freeze",
    )
    return observed


def _prepare_pretruth_plan_from_source(
    source: FrozenSourceModel,
    contract: Mapping[str, Any],
    artifact_manifest_path: Path,
) -> dict[str, Any]:
    nadig_targets = _read_axis(NADIG_DIRECTORY / "perts.txt", "Nadig target axis")
    nadig_genes = _read_axis(NADIG_DIRECTORY / "genes.txt", "Nadig gene axis")
    official_genes = _read_csv_column(
        CHALLENGE_DIRECTORY / "gene_names.csv", "gene_name"
    )
    vcc_targets = _read_csv_column(
        CHALLENGE_DIRECTORY / "pert_counts.csv", "target_gene"
    )
    donor_genes = {
        name: _read_axis(
            DONOR_ROOT / name / "genes.txt",
            f"{name} gene axis",
            require_unique=False,
        )
        for name in DONOR_ORDER
    }
    donor_perts = {
        name: _read_axis(
            DONOR_ROOT / name / "perts.txt",
            f"{name} perturbation axis",
            require_unique=False,
        )
        for name in DONOR_ORDER
    }
    essential = set(
        _read_axis(DONOR_ROOT / "k562_essential" / "perts.txt", "essential target axis")
    )
    artifact_bindings = bind_artifact_manifest(
        artifact_manifest_path,
        required_roles=contract["authorization_bindings"]["required_artifact_roles"],
        contract=contract,
    )
    nadig_artifacts = bind_nadig_artifacts(NADIG_DIRECTORY)
    galleries = reconstruct_gallery_axes(
        source,
        nadig_targets=nadig_targets,
        nadig_genes=nadig_genes,
        official_genes=official_genes,
        vcc_targets=vcc_targets,
        donor_genes=donor_genes,
        donor_perts=donor_perts,
        essential_targets=essential,
        contract=contract,
    )
    return _plan_payload(
        source,
        galleries,
        nadig_targets=nadig_targets,
        nadig_genes=nadig_genes,
        nadig_artifacts=nadig_artifacts,
        artifact_bindings=artifact_bindings,
    )


def prepare_pretruth_plan(
    source_directory: Path = SOURCE_DIRECTORY,
    output_path: Path = PRETRUTH_PLAN_PATH,
    artifact_manifest_path: Path | None = None,
) -> dict[str, Any]:
    """Freeze the exact source/input binding without parsing Nadig outcomes."""

    artifact_manifest = (
        ARTIFACT_MANIFEST_PATH
        if artifact_manifest_path is None
        else Path(artifact_manifest_path)
    )
    _require_production_path(source_directory, SOURCE_DIRECTORY, "source directory")
    _require_production_path(
        artifact_manifest,
        ARTIFACT_MANIFEST_PATH,
        "artifact manifest",
    )
    _require_production_path(output_path, PRETRUTH_PLAN_PATH, "pretruth plan")
    contract = validate_contract()
    source = load_frozen_source_model(source_directory)
    payload = _prepare_pretruth_plan_from_source(
        source,
        contract,
        artifact_manifest,
    )
    write_immutable_json(output_path, payload)
    return payload


def authorize_outcome_phase(
    source_directory: Path,
    plan_path: Path,
    artifact_manifest_path: Path | None = None,
) -> tuple[FrozenSourceModel, dict[str, Any]]:
    """Validate every pretruth prerequisite before any outcome path is touched."""

    artifact_manifest = (
        ARTIFACT_MANIFEST_PATH
        if artifact_manifest_path is None
        else Path(artifact_manifest_path)
    )
    _require_production_path(source_directory, SOURCE_DIRECTORY, "source directory")
    _require_production_path(plan_path, PRETRUTH_PLAN_PATH, "pretruth plan")
    _require(
        Path(plan_path).is_file() and not Path(plan_path).is_symlink(),
        "pretruth plan is missing or unsafe",
    )
    _require_production_path(
        artifact_manifest,
        ARTIFACT_MANIFEST_PATH,
        "artifact manifest",
    )
    contract = validate_contract()
    source = load_frozen_source_model(source_directory)
    plan = load_json(plan_path)
    expected = _prepare_pretruth_plan_from_source(
        source,
        contract,
        artifact_manifest,
    )
    _require(
        plan == expected, "pretruth plan differs from exact metadata reconstruction"
    )
    verify_artifact_bindings(
        plan["artifact_bindings"],
        artifact_manifest,
        required_roles=contract["authorization_bindings"]["required_artifact_roles"],
        contract=contract,
    )
    return source, plan


def load_nadig_outcome_after_authorization(
    source_directory: Path,
    plan_path: Path,
    outcome_directory: Path = NADIG_DIRECTORY,
    *,
    artifact_manifest_path: Path | None = None,
    receipt_path: Path | None = None,
) -> tuple[FrozenSourceModel, dict[str, Any], np.ndarray]:
    """Record one authorized access, then load the reused benchmark effects."""

    receipt = (
        OUTCOME_ACCESS_RECEIPT_PATH if receipt_path is None else Path(receipt_path)
    )
    _require_production_path(outcome_directory, NADIG_DIRECTORY, "Nadig directory")
    _require_production_path(
        receipt,
        OUTCOME_ACCESS_RECEIPT_PATH,
        "outcome-access receipt",
    )
    _require(
        not receipt.exists() and not receipt.is_symlink(),
        "outcome-access receipt already exists; replay refused",
    )
    source, plan = authorize_outcome_phase(
        source_directory,
        plan_path,
        artifact_manifest_path,
    )
    observed_nadig = verify_nadig_artifacts(
        plan["nadig_artifacts"],
        outcome_directory,
    )
    outcome_path = Path(outcome_directory) / "delta.npy"
    receipt_payload = {
        "schema_version": 1,
        "status": "nadig_outcome_access_committed",
        "experiment": "107-jurkat-high-moi",
        "benchmark_role": "reused_development_benchmark_not_independent_confirmation",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "contract_sha256": EXPECTED_CONTRACT_SHA256,
        "validation_implementation": plan["validation_implementation"],
        "source_manifest_sha256": source.manifest_sha256,
        "source_model_sha256": source.model_sha256,
        "pretruth_plan": file_binding(plan_path, root=ROOT),
        "artifact_bindings_sha256": plan["artifact_bindings"]["binding_sha256"],
        "nadig_artifacts": observed_nadig,
    }
    write_immutable_json(receipt, receipt_payload)
    try:
        values = np.load(outcome_path, allow_pickle=False)
    except (OSError, ValueError) as error:
        raise StrictValidationError(f"cannot load Nadig outcome: {error}") from error
    _require(
        values.shape == (int(plan["nadig_target_count"]), int(plan["nadig_gene_count"]))
        and np.all(np.isfinite(values)),
        "Nadig outcome matrix is invalid",
    )
    return source, plan, values


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    pretruth = commands.add_parser("pretruth")
    pretruth.add_argument("--source-directory", type=Path, default=SOURCE_DIRECTORY)
    pretruth.add_argument(
        "--artifact-manifest",
        type=Path,
        default=ARTIFACT_MANIFEST_PATH,
    )
    pretruth.add_argument("--output", type=Path, default=PRETRUTH_PLAN_PATH)
    outcome = commands.add_parser("outcome-check")
    outcome.add_argument("--source-directory", type=Path, default=SOURCE_DIRECTORY)
    outcome.add_argument(
        "--artifact-manifest",
        type=Path,
        default=ARTIFACT_MANIFEST_PATH,
    )
    outcome.add_argument("--plan", type=Path, required=True)
    outcome.add_argument("--outcome-directory", type=Path, default=NADIG_DIRECTORY)
    outcome.add_argument("--receipt", type=Path, default=OUTCOME_ACCESS_RECEIPT_PATH)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        if args.command == "pretruth":
            report = prepare_pretruth_plan(
                args.source_directory,
                args.output,
                args.artifact_manifest,
            )
            print(
                json.dumps(
                    {"status": report["status"], "output": str(args.output)},
                    sort_keys=True,
                )
            )
        else:
            _, plan, outcome = load_nadig_outcome_after_authorization(
                args.source_directory,
                args.plan,
                args.outcome_directory,
                artifact_manifest_path=args.artifact_manifest,
                receipt_path=args.receipt,
            )
            print(
                json.dumps(
                    {
                        "status": "authorized_outcome_loaded",
                        "shape": list(outcome.shape),
                        "galleries": len(plan["galleries"]),
                        "receipt": str(args.receipt),
                    },
                    sort_keys=True,
                )
            )
    except (StrictValidationError, OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

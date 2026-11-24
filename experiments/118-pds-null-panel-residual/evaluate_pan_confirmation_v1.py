#!/usr/bin/env python3
"""One-shot Pan confirmation evaluator for Experiment 118.

Importing this module and running ``authenticate`` never calls ``numpy.load``
on a Pan artifact. ``run`` creates and fsyncs a permanent exclusive access
lock before entering the single Pan-value-loading phase.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import zipfile
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXPERIMENT = "118-pds-null-panel-residual"
PROTOCOL_PATH = HERE / "PAN_CONFIRMATION_PROTOCOL.v1.json"
AUTHORIZATION_PATH = HERE / "PAN_CONFIRMATION_AUTHORIZATION.v1.json"
ENTRYPOINT_PATH = Path(__file__).resolve()
ACCESS_LOCK_PATH = HERE / "work" / "pan-confirmation-v1-access.lock"
RESULT_DIRECTORY = HERE / "work" / "pan-confirmation-v1"
RESULT_SUMMARY = "summary.json"
RESULT_ARRAYS = "resamples.npz"
RESULT_MANIFEST = "result-manifest.json"
CONFIRMATION_TOKEN = "I_CONFIRM_EXPERIMENT_118_PAN_CONFIRMATION_V1"
SELECTED_CANDIDATE = "unweighted_eta1p00"
EXPECTED_TARGETS = 64
BOOTSTRAP_DRAWS = 20_000
BOOTSTRAP_SEED = 1181152026
NULL_DRAWS = 100


class ConfirmationError(RuntimeError):
    """Raised when authentication, evaluation, or immutability fails."""


@dataclass(frozen=True)
class AccessLock:
    path: Path
    sha256: str


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ConfirmationError(message)


def reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_nonfinite_json(value: str) -> None:
    raise ValueError(f"nonfinite JSON constant: {value}")


def load_json(path: Path) -> dict[str, Any]:
    try:
        with path.open(encoding="utf-8") as handle:
            value = json.load(
                handle,
                object_pairs_hook=reject_duplicate_keys,
                parse_constant=reject_nonfinite_json,
            )
    except (OSError, ValueError, TypeError) as error:
        raise ConfirmationError(f"cannot read valid JSON: {path}") from error
    require(isinstance(value, dict), f"JSON root is not an object: {path}")
    return value


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_array(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values))
    require(not array.dtype.hasobject, "object arrays cannot be bound")
    digest = hashlib.sha256()
    digest.update(b"experiment-118-pan-confirmation-array-v1\0")
    digest.update(array.dtype.str.encode("ascii"))
    digest.update(b"\0")
    digest.update(json.dumps(list(array.shape), separators=(",", ":")).encode())
    digest.update(b"\0")
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def canonical_relative(relative: str, label: str) -> Path:
    path = Path(relative)
    require(
        isinstance(relative, str)
        and bool(relative)
        and relative == path.as_posix()
        and not path.is_absolute()
        and ".." not in path.parts,
        f"{label}: noncanonical repository-relative path",
    )
    return path


def source_path(relative: str, label: str) -> Path:
    path = ROOT / canonical_relative(relative, label)
    absolute = Path(os.path.abspath(path))
    require(absolute.is_relative_to(ROOT.resolve()), f"{label}: path escapes repo")
    current = ROOT.resolve()
    for component in absolute.relative_to(ROOT.resolve()).parts:
        current = current / component
        require(not current.is_symlink(), f"{label}: symlink component forbidden")
    require(absolute.is_file(), f"{label}: missing file")
    return absolute


def verify_binding(binding: Mapping[str, Any], label: str) -> Path:
    require(
        isinstance(binding, Mapping)
        and set(binding) == {"path", "bytes", "sha256"},
        f"{label}: binding schema changed",
    )
    path = source_path(str(binding["path"]), label)
    require(path.stat().st_size == int(binding["bytes"]), f"{label}: size changed")
    require(sha256_file(path) == str(binding["sha256"]), f"{label}: hash changed")
    return path


def read_csv_axis(path: Path, column: str) -> list[str]:
    with path.open(newline="") as handle:
        values = [str(row[column]) for row in csv.DictReader(handle)]
    require(values and len(values) == len(set(values)), f"{path}: invalid axis")
    return values


def read_text_axis(path: Path) -> list[str]:
    values = path.read_text().splitlines()
    require(values and len(values) == len(set(values)), f"{path}: invalid axis")
    return values


def validate_protocol() -> dict[str, Any]:
    protocol = load_json(PROTOCOL_PATH)
    require(protocol["candidate"]["name"] == SELECTED_CANDIDATE, "candidate changed")
    require(protocol["candidate"]["reselection"] is False, "reselection enabled")
    require(
        protocol["candidate"]["selection_json_sha256"]
        == "3d849a0fb1227b3e968861fb891d954513a73de12c80130a2b20059f880028d4",
        "selection hash changed",
    )
    require(protocol["intersection"]["expected_target_count"] == EXPECTED_TARGETS, "target count changed")
    require(protocol["bootstrap"]["draws"] == BOOTSTRAP_DRAWS, "bootstrap count changed")
    require(protocol["bootstrap"]["seed"] == BOOTSTRAP_SEED, "bootstrap seed changed")
    require(protocol["nulls"]["target_label"]["draws"] == NULL_DRAWS, "target null count changed")
    require(
        protocol["nulls"]["response_label_within_expression_bin"]["draws"]
        == NULL_DRAWS,
        "response null count changed",
    )
    require(
        protocol["linear_mean_application"]["candidate_residual_exponentiation"]
        is False,
        "candidate residual exponentiation enabled",
    )
    forbidden = ("model_fit", "candidate_reselection", "counts", "cell_eval2", "challenge_build", "submission", "commit", "push")
    require(
        all(protocol["authorization"][key] is False for key in forbidden),
        "forbidden action authorized",
    )
    return protocol


def verify_exact_hash_map(bindings: Mapping[str, str], label: str) -> None:
    require(isinstance(bindings, Mapping) and bool(bindings), f"{label}: empty")
    for relative, expected in sorted(bindings.items()):
        path = source_path(str(relative), f"{label}:{relative}")
        require(sha256_file(path) == str(expected), f"{label}: hash changed: {relative}")


def authenticate_pre_pan() -> dict[str, Any]:
    """Authenticate all inputs without loading any Pan NumPy value array."""
    protocol = validate_protocol()
    authorization = load_json(AUTHORIZATION_PATH)
    require(authorization["status"] == "go_for_independent_review", "authorization status changed")
    require(authorization["scope"] == "experiment_118_pan_confirmation_v1", "authorization scope changed")
    require(authorization["protocol_sha256"] == sha256_file(PROTOCOL_PATH), "protocol hash changed")
    require(
        authorization["implementation_files"][
            "experiments/118-pds-null-panel-residual/evaluate_pan_confirmation_v1.py"
        ]
        == sha256_file(ENTRYPOINT_PATH),
        "evaluator hash changed",
    )
    verify_exact_hash_map(authorization["implementation_files"], "implementation")
    verify_exact_hash_map(authorization["experiment_118_source_inputs"], "Experiment 118 source input")
    verify_exact_hash_map(authorization["experiment_118_frozen_outputs"], "Experiment 118 frozen output")

    source_auth_path = source_path(
        "experiments/118-pds-null-panel-residual/SOURCE_ONLY_AUTHORIZATION.json",
        "Experiment 118 source authorization",
    )
    source_auth = load_json(source_auth_path)
    require(
        source_auth["allowed_input_files"]
        == authorization["experiment_118_source_inputs"],
        "source authorization input map differs",
    )
    freeze_path = source_path(
        "experiments/118-pds-null-panel-residual/work/source-freeze-v1/freeze-manifest.json",
        "Experiment 118 freeze manifest",
    )
    freeze = load_json(freeze_path)
    require(freeze["source_files"] == authorization["experiment_118_source_inputs"], "freeze source map differs")
    for name, binding in freeze["artifacts"].items():
        relative = f"experiments/118-pds-null-panel-residual/work/source-freeze-v1/{name}"
        require(
            authorization["experiment_118_frozen_outputs"].get(relative)
            == binding["sha256"],
            f"freeze output binding differs: {name}",
        )
    selection_path = source_path(
        "experiments/118-pds-null-panel-residual/work/source-freeze-v1/selection.json",
        "Experiment 118 selection",
    )
    require(
        sha256_file(selection_path)
        == "3d849a0fb1227b3e968861fb891d954513a73de12c80130a2b20059f880028d4",
        "selection hash changed",
    )
    selection = load_json(selection_path)
    require(
        selection["observed_selected_before_null_gates"] == SELECTED_CANDIDATE
        and selection["authorized_selected_after_all_source_gates"]
        == SELECTED_CANDIDATE
        and selection["nonidentity_survives_source_only_gates"] is True,
        "frozen selected candidate changed",
    )

    pan_paths = {
        name: verify_binding(binding, f"Pan compact {name}")
        for name, binding in authorization["experiment_115_compact"].items()
    }
    artifact_manifest = load_json(pan_paths["artifact_manifest"])
    compact_receipt = load_json(pan_paths["compact_receipt"])
    require(
        compact_receipt["artifact_manifest"]["sha256"]
        == authorization["experiment_115_compact"]["artifact_manifest"]["sha256"],
        "compact receipt manifest hash differs",
    )
    for key, filename in (
        ("delta", "delta.npy"),
        ("measured_mask", "measured_mask.npy"),
        ("b_ctrl", "b_ctrl.npy"),
        ("genes", "genes.txt"),
        ("perts", "perts.txt"),
    ):
        require(
            artifact_manifest["files"][filename]["sha256"]
            == authorization["experiment_115_compact"][key]["sha256"],
            f"Pan artifact manifest differs for {filename}",
        )

    canonical_genes = read_csv_axis(
        source_path("data/challenge/gene_names.csv", "canonical genes"),
        "gene_name",
    )
    canonical_targets = read_csv_axis(
        source_path("data/challenge/pert_counts.csv", "canonical targets"),
        "target_gene",
    )
    pan_genes = read_text_axis(pan_paths["genes"])
    pan_targets = read_text_axis(pan_paths["perts"])
    frozen_vcc_targets = read_text_axis(pan_paths["vcc_targets"])
    require(tuple(pan_genes) == tuple(canonical_genes), "Pan gene axis is not canonical")
    intersection = [target for target in canonical_targets if target in set(pan_targets)]
    require(len(intersection) == EXPECTED_TARGETS, "Pan exact-symbol target intersection changed")
    require(
        tuple(frozen_vcc_targets) == tuple(intersection),
        "Pan frozen VCC-target intersection changed",
    )
    require(
        not os.path.lexists(ACCESS_LOCK_PATH),
        "permanent Pan access lock already exists",
    )
    require(
        not os.path.lexists(RESULT_DIRECTORY),
        "immutable Pan result directory already exists",
    )
    return {
        "protocol": protocol,
        "authorization": authorization,
        "authorization_sha256": sha256_file(AUTHORIZATION_PATH),
        "pan_paths": pan_paths,
        "canonical_genes": canonical_genes,
        "canonical_targets": canonical_targets,
        "pan_targets": pan_targets,
        "exact_target_intersection": intersection,
    }


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float):
        require(math.isfinite(value), "result contains nonfinite float")
    return value


def write_json_exclusive(path: Path, payload: Mapping[str, Any]) -> None:
    encoded = (
        json.dumps(json_safe(payload), indent=2, sort_keys=True, allow_nan=False)
        + "\n"
    ).encode()
    descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o400)
    try:
        offset = 0
        while offset < len(encoded):
            offset += os.write(descriptor, encoded[offset:])
        os.fsync(descriptor)
        os.fchmod(descriptor, 0o444)
    finally:
        os.close(descriptor)


def create_permanent_access_lock(authenticated: Mapping[str, Any]) -> AccessLock:
    ACCESS_LOCK_PATH.parent.mkdir(parents=True, exist_ok=True)
    require(not RESULT_DIRECTORY.exists(), "result directory appeared before lock")
    payload = {
        "schema_version": 1,
        "experiment": EXPERIMENT,
        "stage": "pan_confirmation_v1",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "pid": os.getpid(),
        "authorization_sha256": authenticated["authorization_sha256"],
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "entrypoint_sha256": sha256_file(ENTRYPOINT_PATH),
        "pan_value_access_occurs_only_after_this_fsynced_lock": True,
        "permanent_even_on_failure": True,
        "single_attempt": True,
    }
    write_json_exclusive(ACCESS_LOCK_PATH, payload)
    parent = os.open(ACCESS_LOCK_PATH.parent, os.O_RDONLY)
    try:
        os.fsync(parent)
    finally:
        os.close(parent)
    return AccessLock(ACCESS_LOCK_PATH, sha256_file(ACCESS_LOCK_PATH))


def verify_lock(lock: AccessLock) -> None:
    require(lock.path == ACCESS_LOCK_PATH, "wrong access lock path")
    require(lock.path.is_file(), "access lock missing")
    require(sha256_file(lock.path) == lock.sha256, "access lock changed")


def load_pan_values_once(
    authenticated: Mapping[str, Any], lock: AccessLock
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The evaluator's sole Pan NumPy-loading boundary."""
    verify_lock(lock)
    paths = authenticated["pan_paths"]
    # Re-hash immediately after the permanent lock to close the preflight/load
    # gap. Hashing is opaque-byte authentication; these are the only np.load
    # calls for Pan values.
    for key in ("delta", "measured_mask", "b_ctrl"):
        binding = authenticated["authorization"]["experiment_115_compact"][key]
        require(sha256_file(paths[key]) == binding["sha256"], f"{key} changed after lock")
    delta = np.load(paths["delta"], mmap_mode="r", allow_pickle=False)
    measured = np.load(paths["measured_mask"], mmap_mode="r", allow_pickle=False)
    b_ctrl = np.load(paths["b_ctrl"], mmap_mode="r", allow_pickle=False)
    return delta, measured, b_ctrl


def exact_zero_sum(values: np.ndarray) -> bool:
    return math.fsum(map(float, values)) == 0.0


def repair_to_sum(
    values: np.ndarray, mutable: np.ndarray, target_sum: float
) -> None:
    indices = np.flatnonzero(mutable)
    require(indices.size > 0, "no mutable coordinates for sum repair")
    order = [int(indices[-1])]
    order.extend(
        int(index)
        for index in indices[np.argsort(np.abs(values[indices]), kind="stable")]
        if int(index) != order[0]
    )
    for index in order:
        for _ in range(4):
            residual = math.fsum(map(float, values)) - target_sum
            if residual == 0.0:
                return
            before = float(values[index])
            values[index] -= residual
            if float(values[index]) == before:
                break
    raise ConfirmationError("could not enforce exact compensated sum")


def repair_zero_sum(values: np.ndarray, mutable: np.ndarray) -> None:
    repair_to_sum(values, mutable, 0.0)


def intersected_panel_residual(
    selected_panel: np.ndarray,
    destination_target: int,
    measured_panel: np.ndarray,
    *,
    source_target: int | None = None,
    response_mapping: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    n = selected_panel.shape[0]
    require(selected_panel.shape == (n, n), "selected panel residual must be square")
    require(measured_panel.shape == (n,), "measured panel mask shape changed")
    source = destination_target if source_target is None else int(source_target)
    require(0 <= source < n, "source target out of bounds")
    mapping = (
        np.arange(n, dtype=np.int64)
        if response_mapping is None
        else np.asarray(response_mapping, dtype=np.int64)
    )
    require(mapping.shape == (n,), "response mapping shape changed")
    row = np.asarray(selected_panel[source, mapping], dtype=np.float64).copy()
    mutable = np.asarray(measured_panel, dtype=bool).copy()
    mutable[destination_target] = False
    mutable &= np.isfinite(row)
    row[~mutable] = 0.0
    if np.any(mutable):
        row[mutable] -= float(np.mean(row[mutable]))
        repair_zero_sum(row, mutable)
    require(row[destination_target] == 0.0, "target self residual changed")
    require(np.all(row[~mutable] == 0.0), "non-intersection residual changed")
    require(exact_zero_sum(row), "intersected residual is not exact zero-sum")
    return row, mutable


def validate_frozen_permutations(
    target_permutations: np.ndarray,
    response_permutations: np.ndarray,
    expression_bins: np.ndarray,
) -> None:
    identity = np.arange(300)
    require(
        target_permutations.shape == response_permutations.shape == (100, 300),
        "frozen permutation shape changed",
    )
    require(expression_bins.shape == (300,), "frozen expression-bin shape changed")
    for row in target_permutations:
        require(np.array_equal(np.sort(row), identity), "target null is not bijective")
        require(np.all(row != identity), "target null has a fixed point")
    for row in response_permutations:
        require(np.array_equal(np.sort(row), identity), "response null is not bijective")
        absent = expression_bins < 0
        require(np.all(row[absent] == identity[absent]), "source-absent response moved")
        for bin_index in sorted(set(map(int, expression_bins)) - {-1}):
            members = np.flatnonzero(expression_bins == bin_index)
            require(
                set(map(int, row[members])) == set(map(int, members)),
                "response null crossed an expression bin",
            )
            require(np.all(row[members] != members), "response null has a fixed point")


def profile_mass(values: np.ndarray, panel_mask: np.ndarray) -> float:
    return math.fsum(
        (
            math.fsum(map(float, values[panel_mask])),
            math.fsum(map(float, values[~panel_mask])),
        )
    )


def apply_additive_linear_mean(
    control_linear: np.ndarray,
    panel_indices: np.ndarray,
    panel_residual: np.ndarray,
    target_index: int,
    mutable_panel: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[str, bool]]:
    control = np.asarray(control_linear, dtype=np.float64)
    candidate = control.copy()
    require(panel_residual.shape == panel_indices.shape, "panel residual shape changed")
    mutable = (
        np.asarray(panel_residual != 0, dtype=bool)
        if mutable_panel is None
        else np.asarray(mutable_panel, dtype=bool).copy()
    )
    require(mutable.shape == panel_residual.shape, "mutable panel mask changed")
    mutable[target_index] = False
    panel_candidate = control[panel_indices].copy()
    panel_candidate += panel_residual
    repair_to_sum(
        panel_candidate,
        mutable,
        math.fsum(map(float, control[panel_indices])),
    )
    candidate[panel_indices] = panel_candidate
    applied_residual = panel_candidate - control[panel_indices]
    panel_mask = np.zeros(len(control), dtype=bool)
    panel_mask[panel_indices] = True
    control_mass = profile_mass(control, panel_mask)
    candidate_mass = profile_mass(candidate, panel_mask)
    invariants = {
        "nonpanel_linear_mean_identity": bool(np.array_equal(candidate[~panel_mask], control[~panel_mask])),
        "target_diagonal_linear_mean_identity": bool(
            candidate[panel_indices[target_index]] == control[panel_indices[target_index]]
        ),
        "panel_mass_conservation": bool(
            math.fsum(map(float, candidate[panel_mask]))
            == math.fsum(map(float, control[panel_mask]))
        ),
        "total_mass_conservation": bool(candidate_mass == control_mass),
        "candidate_nonnegative": bool(np.all(candidate >= 0)),
    }
    if control_mass > 0 and candidate_mass > 0:
        normalized_control = control / control_mass
        normalized_candidate = candidate / candidate_mass
        invariants["pds_visible_normalized_linear_identity"] = bool(
            np.array_equal(
                normalized_candidate[~panel_mask], normalized_control[~panel_mask]
            )
        )
        invariants["pds_visible_normalized_log1p_identity"] = bool(
            np.array_equal(
                np.log1p(50_000.0 * normalized_candidate[~panel_mask]),
                np.log1p(50_000.0 * normalized_control[~panel_mask]),
            )
        )
    else:
        invariants["pds_visible_normalized_linear_identity"] = False
        invariants["pds_visible_normalized_log1p_identity"] = False
    invariants["all_pass"] = all(invariants.values())
    return candidate, applied_residual, invariants


def application_invariant_summary(
    selected_panel: np.ndarray,
    control_linear: np.ndarray,
    panel_indices: np.ndarray,
    measured_panel: np.ndarray,
    destination_targets: np.ndarray,
) -> dict[str, bool]:
    require(
        measured_panel.shape == (len(destination_targets), len(panel_indices)),
        "application invariant input shape changed",
    )
    invariant_rows: list[dict[str, bool]] = []
    for row, target in enumerate(destination_targets):
        residual, mutable = intersected_panel_residual(
            selected_panel,
            int(target),
            measured_panel[row],
        )
        _, _, invariants = apply_additive_linear_mean(
            control_linear,
            panel_indices,
            residual,
            int(target),
            mutable,
        )
        invariant_rows.append(invariants)
    require(invariant_rows, "no application-invariant rows")
    names = sorted(invariant_rows[0])
    summary = {
        name: all(row[name] for row in invariant_rows)
        for name in names
    }
    summary["all_pass"] = all(summary.values())
    return summary


def target_metrics(
    residual: np.ndarray, truth: np.ndarray, eligible: np.ndarray
) -> tuple[float, float, int]:
    indices = np.flatnonzero(eligible & np.isfinite(truth) & (truth != 0))
    require(indices.size > 0, "target has no evaluable response coordinate")
    y = truth[indices]
    p = residual[indices]
    denominator = float(np.abs(y).sum())
    require(denominator > 0 and math.isfinite(denominator), "invalid truth magnitude")
    credit = float(np.sign(y) @ p) / denominator
    selected_count = max(1, int(math.ceil(0.20 * len(indices))))
    order = np.argsort(-np.abs(p), kind="stable")[:selected_count]
    weights = np.abs(y[order])
    correct = (np.sign(p[order]) == np.sign(y[order])).astype(np.float64)
    precision = float(correct @ weights) / float(weights.sum())
    return credit, precision, len(indices)


def evaluate_mapping(
    selected_panel: np.ndarray,
    control_panel: np.ndarray,
    truth_panel: np.ndarray,
    measured_panel: np.ndarray,
    destination_targets: np.ndarray,
    *,
    target_mapping: np.ndarray | None = None,
    response_mapping: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    credits = np.empty(len(destination_targets), dtype=np.float64)
    precisions = np.empty(len(destination_targets), dtype=np.float64)
    eligible_counts = np.empty(len(destination_targets), dtype=np.int32)
    for row, target_index in enumerate(destination_targets):
        source_target = (
            int(target_index)
            if target_mapping is None
            else int(target_mapping[int(target_index)])
        )
        residual, mutable = intersected_panel_residual(
            selected_panel,
            int(target_index),
            measured_panel[row],
            source_target=source_target,
            response_mapping=response_mapping,
        )
        panel_candidate = np.asarray(control_panel, dtype=np.float64).copy()
        panel_candidate += residual
        repair_to_sum(
            panel_candidate,
            mutable,
            math.fsum(map(float, control_panel)),
        )
        applied_residual = panel_candidate - control_panel
        credit, precision, count = target_metrics(
            applied_residual, truth_panel[row], mutable
        )
        credits[row] = credit
        precisions[row] = precision
        eligible_counts[row] = count
    return credits, precisions, eligible_counts


def bootstrap_targets(
    credits: np.ndarray, precisions: np.ndarray
) -> tuple[np.ndarray, np.ndarray, dict[str, float], dict[str, float]]:
    require(credits.shape == precisions.shape and credits.ndim == 1, "bootstrap input mismatch")
    rng = np.random.Generator(np.random.PCG64(BOOTSTRAP_SEED))
    credit_draws = np.empty(BOOTSTRAP_DRAWS, dtype=np.float64)
    precision_draws = np.empty(BOOTSTRAP_DRAWS, dtype=np.float64)
    for draw in range(BOOTSTRAP_DRAWS):
        indices = rng.integers(0, len(credits), size=len(credits))
        credit_draws[draw] = float(np.mean(credits[indices]))
        precision_draws[draw] = float(np.mean(precisions[indices]))
    return (
        credit_draws,
        precision_draws,
        {
            "lower_95": float(np.quantile(credit_draws, 0.025, method="linear")),
            "upper_95": float(np.quantile(credit_draws, 0.975, method="linear")),
        },
        {
            "lower_95": float(np.quantile(precision_draws, 0.025, method="linear")),
            "upper_95": float(np.quantile(precision_draws, 0.975, method="linear")),
        },
    )


def null_metrics(
    selected_panel: np.ndarray,
    control_panel: np.ndarray,
    truth_panel: np.ndarray,
    measured_panel: np.ndarray,
    destination_targets: np.ndarray,
    target_permutations: np.ndarray,
    response_permutations: np.ndarray,
) -> dict[str, np.ndarray]:
    require(target_permutations.shape == (NULL_DRAWS, 300), "target null bank changed")
    require(response_permutations.shape == (NULL_DRAWS, 300), "response null bank changed")
    output = {
        "target_credit": np.empty(NULL_DRAWS, dtype=np.float64),
        "target_precision": np.empty(NULL_DRAWS, dtype=np.float64),
        "response_credit": np.empty(NULL_DRAWS, dtype=np.float64),
        "response_precision": np.empty(NULL_DRAWS, dtype=np.float64),
    }
    for draw in range(NULL_DRAWS):
        credit, precision, _ = evaluate_mapping(
            selected_panel,
            control_panel,
            truth_panel,
            measured_panel,
            destination_targets,
            target_mapping=target_permutations[draw],
        )
        output["target_credit"][draw] = float(np.mean(credit))
        output["target_precision"][draw] = float(np.mean(precision))
        credit, precision, _ = evaluate_mapping(
            selected_panel,
            control_panel,
            truth_panel,
            measured_panel,
            destination_targets,
            response_mapping=response_permutations[draw],
        )
        output["response_credit"][draw] = float(np.mean(credit))
        output["response_precision"][draw] = float(np.mean(precision))
    return output


def promotion_report(
    macro_credit: float,
    macro_precision: float,
    credit_interval: Mapping[str, float],
    precision_interval: Mapping[str, float],
    nulls: Mapping[str, np.ndarray],
    target_count: int,
    all_invariants: bool,
) -> dict[str, Any]:
    thresholds = {
        "exact_target_count": target_count == EXPECTED_TARGETS,
        "all_targets_evaluable": target_count == EXPECTED_TARGETS,
        "all_application_invariants": all_invariants,
        "macro_credit_positive": macro_credit > 0.0,
        "credit_bootstrap_lower_positive": float(credit_interval["lower_95"]) > 0.0,
        "macro_precision_above_0p55": macro_precision > 0.55,
        "precision_bootstrap_lower_above_0p50": float(precision_interval["lower_95"]) > 0.50,
        "credit_above_target_null_p95": macro_credit
        > float(np.percentile(nulls["target_credit"], 95, method="linear")),
        "credit_above_response_null_p95": macro_credit
        > float(np.percentile(nulls["response_credit"], 95, method="linear")),
        "precision_above_target_null_p95": macro_precision
        > float(np.percentile(nulls["target_precision"], 95, method="linear")),
        "precision_above_response_null_p95": macro_precision
        > float(np.percentile(nulls["response_precision"], 95, method="linear")),
    }
    return {
        "decision": "GO" if all(thresholds.values()) else "NO_GO",
        "all_gates_pass": all(thresholds.values()),
        "gates": thresholds,
        "thresholds_frozen_before_pan_load": True,
    }


def write_npz_exclusive(path: Path, arrays: Mapping[str, np.ndarray]) -> None:
    require(not path.exists(), f"refusing to overwrite {path}")
    with zipfile.ZipFile(
        path,
        mode="x",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=6,
        allowZip64=True,
    ) as archive:
        for name in sorted(arrays):
            value = np.ascontiguousarray(np.asarray(arrays[name]))
            require(not value.dtype.hasobject, f"object result array: {name}")
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


def output_binding(path: Path) -> dict[str, Any]:
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def write_immutable_results(
    summary: Mapping[str, Any],
    arrays: Mapping[str, np.ndarray],
    authenticated: Mapping[str, Any],
    lock: AccessLock,
) -> dict[str, Any]:
    require(not RESULT_DIRECTORY.exists(), "result directory already exists")
    RESULT_DIRECTORY.mkdir(mode=0o700)
    try:
        summary_path = RESULT_DIRECTORY / RESULT_SUMMARY
        arrays_path = RESULT_DIRECTORY / RESULT_ARRAYS
        manifest_path = RESULT_DIRECTORY / RESULT_MANIFEST
        write_npz_exclusive(arrays_path, arrays)
        write_json_exclusive(summary_path, summary)
        manifest = {
            "schema_version": 1,
            "experiment": EXPERIMENT,
            "stage": "pan_confirmation_v1",
            "status": "immutable_one_shot_result",
            "access_lock": {
                "path": lock.path.relative_to(ROOT).as_posix(),
                "sha256": lock.sha256,
            },
            "authorization_sha256": authenticated["authorization_sha256"],
            "protocol_sha256": sha256_file(PROTOCOL_PATH),
            "entrypoint_sha256": sha256_file(ENTRYPOINT_PATH),
            "outputs": {
                "summary": output_binding(summary_path),
                "arrays": output_binding(arrays_path),
            },
        }
        write_json_exclusive(manifest_path, manifest)
        directory = os.open(RESULT_DIRECTORY, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
        os.chmod(RESULT_DIRECTORY, 0o555)
        parent = os.open(RESULT_DIRECTORY.parent, os.O_RDONLY)
        try:
            os.fsync(parent)
        finally:
            os.close(parent)
        return manifest
    except BaseException:
        if RESULT_DIRECTORY.exists():
            for child in RESULT_DIRECTORY.iterdir():
                try:
                    child.chmod(0o400)
                except OSError:
                    pass
            try:
                RESULT_DIRECTORY.chmod(0o500)
            except OSError:
                pass
        raise


def run_one_shot(confirm: str) -> dict[str, Any]:
    require(confirm == CONFIRMATION_TOKEN, "exact --confirm-pan-access token required")
    authenticated = authenticate_pre_pan()
    lock = create_permanent_access_lock(authenticated)
    try:
        delta, measured_mask, b_ctrl = load_pan_values_once(authenticated, lock)
        require(delta.shape == (2354, 18533) and delta.dtype == np.dtype("<f4"), "Pan delta schema changed")
        require(
            measured_mask.shape == delta.shape
            and measured_mask.dtype == np.dtype("bool"),
            "Pan measured-mask schema changed",
        )
        require(b_ctrl.shape == (18533,) and b_ctrl.dtype == np.dtype("<f4"), "Pan control schema changed")

        frozen_root = HERE / "work" / "source-freeze-v1"
        selected_full = np.load(
            frozen_root / "selected_residual_full.npy", allow_pickle=False
        )
        target_permutations = np.load(
            frozen_root / "target_label_permutations.npy", allow_pickle=False
        )
        response_permutations = np.load(
            frozen_root / "response_label_permutations.npy", allow_pickle=False
        )
        expression_bins = np.load(
            frozen_root / "response_expression_bins.npy", allow_pickle=False
        )
        require(selected_full.shape == (300, 18533), "selected residual shape changed")
        require(target_permutations.shape == (100, 300), "target null shape changed")
        require(response_permutations.shape == (100, 300), "response null shape changed")
        require(expression_bins.shape == (300,), "expression-bin shape changed")
        validate_frozen_permutations(
            target_permutations, response_permutations, expression_bins
        )

        genes = authenticated["canonical_genes"]
        targets = authenticated["canonical_targets"]
        pan_targets = authenticated["pan_targets"]
        gene_index = {gene: index for index, gene in enumerate(genes)}
        target_index = {target: index for index, target in enumerate(targets)}
        pan_row = {target: index for index, target in enumerate(pan_targets)}
        intersection = authenticated["exact_target_intersection"]
        destination_targets = np.asarray(
            [target_index[target] for target in intersection], dtype=np.int64
        )
        pan_rows = np.asarray([pan_row[target] for target in intersection], dtype=np.int64)
        panel_indices = np.asarray([gene_index[target] for target in targets], dtype=np.int64)
        panel_mask = np.zeros(len(genes), dtype=bool)
        panel_mask[panel_indices] = True
        require(np.all(selected_full[:, ~panel_mask] == 0), "selected residual changed PDS-visible coordinates")
        require(np.all(selected_full[np.arange(300), panel_indices] == 0), "selected residual diagonal changed")
        selected_panel = np.asarray(selected_full[:, panel_indices], dtype=np.float64)

        control_linear = np.expm1(np.asarray(b_ctrl, dtype=np.float64))
        require(np.all(np.isfinite(control_linear)) and np.all(control_linear >= 0), "invalid Pan control")
        pan_delta_panel = np.asarray(delta[np.ix_(pan_rows, panel_indices)], dtype=np.float64)
        truth_panel = np.expm1(
            np.asarray(b_ctrl[panel_indices], dtype=np.float64)[None, :]
            + pan_delta_panel
        ) - control_linear[panel_indices][None, :]
        measured_panel = np.asarray(
            measured_mask[np.ix_(pan_rows, panel_indices)], dtype=bool
        )
        require(np.all(np.isfinite(truth_panel)), "Pan linear truth is nonfinite")

        observed_credit, observed_precision, eligible_counts = evaluate_mapping(
            selected_panel,
            control_linear[panel_indices],
            truth_panel,
            measured_panel,
            destination_targets,
        )
        invariant_summary = application_invariant_summary(
            selected_panel,
            control_linear,
            panel_indices,
            measured_panel,
            destination_targets,
        )

        credit_draws, precision_draws, credit_interval, precision_interval = (
            bootstrap_targets(observed_credit, observed_precision)
        )
        nulls = null_metrics(
            selected_panel,
            control_linear[panel_indices],
            truth_panel,
            measured_panel,
            destination_targets,
            target_permutations,
            response_permutations,
        )
        macro_credit = float(np.mean(observed_credit))
        macro_precision = float(np.mean(observed_precision))
        promotion = promotion_report(
            macro_credit,
            macro_precision,
            credit_interval,
            precision_interval,
            nulls,
            len(observed_credit),
            bool(invariant_summary["all_pass"]),
        )
        summary = {
            "schema_version": 1,
            "experiment": EXPERIMENT,
            "stage": "pan_confirmation_v1",
            "status": "one_shot_pan_confirmation_complete",
            "created_at_utc": datetime.now(UTC).isoformat(),
            "candidate": SELECTED_CANDIDATE,
            "candidate_reselected": False,
            "application": "direct_additive_linear_mean_no_residual_exponentiation",
            "target_count": len(intersection),
            "targets": intersection,
            "eligible_response_count_per_target": eligible_counts.tolist(),
            "macro_target_equal_signed_credit": macro_credit,
            "macro_target_equal_top20_weighted_sign_precision": macro_precision,
            "bootstrap": {
                "draws": BOOTSTRAP_DRAWS,
                "seed": BOOTSTRAP_SEED,
                "signed_credit": credit_interval,
                "top20_weighted_sign_precision": precision_interval,
            },
            "nulls": {
                "draws_each": NULL_DRAWS,
                "target_label": {
                    "signed_credit_p95": float(np.percentile(nulls["target_credit"], 95, method="linear")),
                    "top20_precision_p95": float(np.percentile(nulls["target_precision"], 95, method="linear")),
                },
                "response_label_within_expression_bin": {
                    "signed_credit_p95": float(np.percentile(nulls["response_credit"], 95, method="linear")),
                    "top20_precision_p95": float(np.percentile(nulls["response_precision"], 95, method="linear")),
                },
            },
            "application_invariants": invariant_summary,
            "promotion": promotion,
            "access_lock": {
                "path": lock.path.relative_to(ROOT).as_posix(),
                "sha256": lock.sha256,
            },
            "pan_value_load_phases": 1,
            "counts_generated": False,
            "cell_eval2_run": False,
            "challenge_built": False,
            "submission_performed": False,
        }
        arrays = {
            "target_indices": destination_targets,
            "per_target_signed_credit": observed_credit,
            "per_target_top20_weighted_sign_precision": observed_precision,
            "eligible_response_counts": eligible_counts,
            "bootstrap_signed_credit": credit_draws,
            "bootstrap_top20_weighted_sign_precision": precision_draws,
            **nulls,
        }
        manifest = write_immutable_results(summary, arrays, authenticated, lock)
        return {
            "summary": summary,
            "result_manifest_sha256": sha256_file(
                RESULT_DIRECTORY / RESULT_MANIFEST
            ),
            "manifest": manifest,
        }
    except BaseException as error:
        if not RESULT_DIRECTORY.exists():
            failure = {
                "schema_version": 1,
                "experiment": EXPERIMENT,
                "stage": "pan_confirmation_v1",
                "status": "one_shot_pan_confirmation_failed_after_permanent_lock",
                "created_at_utc": datetime.now(UTC).isoformat(),
                "error_type": type(error).__name__,
                "error": str(error),
                "access_lock": {
                    "path": lock.path.relative_to(ROOT).as_posix(),
                    "sha256": lock.sha256,
                },
                "rerun_authorized": False,
            }
            RESULT_DIRECTORY.mkdir(mode=0o700)
            write_json_exclusive(RESULT_DIRECTORY / "failure.json", failure)
            os.chmod(RESULT_DIRECTORY, 0o555)
        raise


def authenticate_command() -> None:
    authenticated = authenticate_pre_pan()
    print(
        json.dumps(
            {
                "status": "GO-for-independent-review",
                "pan_values_loaded": False,
                "authorization_sha256": authenticated["authorization_sha256"],
                "protocol_sha256": sha256_file(PROTOCOL_PATH),
                "entrypoint_sha256": sha256_file(ENTRYPOINT_PATH),
                "exact_symbol_target_intersection_count": len(
                    authenticated["exact_target_intersection"]
                ),
                "access_lock_absent": True,
                "result_directory_absent": True,
            },
            indent=2,
            sort_keys=True,
        )
    )


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("authenticate", "run"))
    parser.add_argument("--confirm-pan-access")
    args = parser.parse_args(argv)
    if args.command == "authenticate":
        require(
            args.confirm_pan_access is None,
            "authenticate does not accept a Pan access token",
        )
        authenticate_command()
    else:
        result = run_one_shot(str(args.confirm_pan_access))
        print(json.dumps(json_safe(result["summary"]), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Experiment 141: cross-context calibration of same-target Pan residuals.

The 613-target development split is already consumed. The 410-target split is
deliberately unreachable from this program.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
import resource
import sys
import tempfile
from pathlib import Path
from typing import Any, Mapping, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL = HERE / "preregistered_protocol.json"
RUNNER = HERE / "run.py"
TESTS = HERE / "test_run.py"
PLAN_DIR = HERE / "work/preregistered-plan-v1"
SOURCE_MANIFEST = HERE / "SOURCE_INPUT_MANIFEST.v1.json"
AUTHORIZATION = HERE / "PREOUTCOME_AUTHORIZATION.v1.json"
ACCESS_LOCK = HERE / "work/development-outcome-access.v1.lock"
RESULT_DIR = HERE / "work/development-result-v1"
DEVELOPMENT_FREEZE = HERE / "DEVELOPMENT_RESULT_FREEZE.v1.json"
INFERENCE_DIR = HERE / "work/conditional-inference-v1"
INFERENCE_FREEZE = HERE / "INFERENCE_RESULT_FREEZE.v1.json"

EXP132 = ROOT / "experiments/132-response-decomposition-ridge"
EXP140 = ROOT / "experiments/140-batch-matched-controls"
PARENT132_RUN = EXP132 / "run.py"
PARENT140_RUN = EXP140 / "run.py"
PARENT140_PLAN = EXP140 / "work/preregistered-plan-v2"
PARENT140_IDENTITY = EXP140 / "work/identity-reproduction-v5"
DEV_TARGETS = EXP132 / "work/preregistered-freeze-v2/selection_targets.txt"
SEALED_TEST_TARGETS = EXP132 / "work/preregistered-freeze-v2/test_targets.txt"
RESPONSE_GENES = EXP132 / "work/preregistered-freeze-v2/response_genes.txt"
PAN_ROOT = EXP132.parent / "115-gse295214-h1/work/compact-donor"

EXPECTED_SHA256 = {
    PARENT132_RUN: "e78e079fc3b4a715e8b412b0a48bfe44a2f3f6cfcadb8a65a789351c53877778",
    PARENT140_RUN: "3bea7cb91e7b5d86cb41d2921b156de164e79d9c3dc0837c20f416c4d3fc1a93",
    PARENT140_PLAN / "plan-arrays.npz": "bdd17742d8691ff5174631c11c00c7d3dfb23f51c9e225ca672549608d466d5b",
    PARENT140_PLAN / "plan-manifest.json": "2f91872b1318ab7f0ee608600f5acfaaad0397f508e833707e28c2f379e6007c",
    PARENT140_IDENTITY / "artifact-manifest.json": "1967a66f463ba78b02c674c09274c9c399fb68b405ff9b0deb045ec4772a5649",
    EXP140 / "SOURCE_HASH_MANIFEST.v5.json": "70c860e308bdf433fbd8e6812e95c5e24b090173f2cc3ac4e6d689ab9195a149",
    PAN_ROOT / "artifact-manifest.json": "773895e0dcff2b6ae4aa9db22f98cd1376315a62b1a7cd193ef559330cf25cd2",
    PAN_ROOT / "compact-donor-receipt.json": "2d24b8fcf86cf6fd6c5fc0a02d1da5525393692e9f67d70f865c72189f875348",
    PAN_ROOT / "delta.npy": "d0fa6555f9a62ce503d2ad55c8bfda954322657b52dcd48a24e29c330e718c7a",
    PAN_ROOT / "measured_mask.npy": "da468865f861be0714037c33357b93b9f747351f7478866e6e0dd6dff9d1ebcb",
    PAN_ROOT / "perts.txt": "7060076272d411d4ee00b898cc26eedcb1eaaac5f5e94e2c60f3b45978320747",
    PAN_ROOT / "genes.txt": "c2b8dec58c02810d64018164155b655c4e7c7be451a3ae579965db193af33723",
    DEV_TARGETS: "74a070a538455ffc2a528dae3094d21e047905463779b2a5206d681f5ca6ebb5",
    SEALED_TEST_TARGETS: "bd37bcacb4eea90da11902bbb8f4736068a7ff67c001739e3407de2b4bf73fc6",
    RESPONSE_GENES: "24743754eead5b8ddd81634bd12feac1df431eb2b9581197755b75c4c90949c8",
}

FAMILIES = ("k562", "rpe1", "jurkat")
FAMILY_DATASETS = {
    "k562": ("k562_gwps", "k562_essential"),
    "rpe1": ("rpe1",),
    "jurkat": ("nadig_jurkat",),
}
DATASET_FAMILY = {
    dataset: family
    for family, datasets in FAMILY_DATASETS.items()
    for dataset in datasets
}
SCORED_DATASETS = tuple(DATASET_FAMILY)
CALIBRATION_TEACHERS = {
    "k562": ("rpe1", "jurkat"),
    "rpe1": ("jurkat",),
    "jurkat": ("rpe1",),
}
CANDIDATES = {
    "identity": None,
    "scalar_control": 0.0,
    "diag_rho025": 0.25,
    "diag_rho050": 0.50,
    "diag_rho075": 0.75,
}
PRIMARY_CANDIDATE = "diag_rho050"
DEV_COUNT = 613
RESPONSE_COUNT = 4992
FOLDS = 5
NULL_DRAWS = 49
BOOTSTRAP_DRAWS = 20_000
TARGET_NULL_SEED = 141_202_610
COORD_NULL_SEED = 141_202_611
EPS = 1e-12
MIN_GAIN = 0.003
MIN_DIAG_OVER_SCALAR = 0.001
MAX_CONCENTRATION = 0.5
GIB = 1024**3
DEVELOPMENT_CONFIRM = "I_CONFIRM_EXPERIMENT_141_DEVELOPMENT_OUTCOME"


class ContractError(RuntimeError):
    """Raised when a frozen experiment boundary is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def import_path(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha256_file(path: Path, block: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(block):
            digest.update(chunk)
    return digest.hexdigest()


def read_lines(path: Path) -> tuple[str, ...]:
    values = tuple(
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    )
    require(bool(values), f"{path}: empty")
    require(len(values) == len(set(values)), f"{path}: duplicate values")
    return values


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".npz"
    )
    os.close(descriptor)
    try:
        with Path(temporary).open("wb") as handle:
            np.savez_compressed(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if sys.platform == "darwin" else value / 1024**2


def memory_check(label: str, limit_gib: float) -> None:
    peak = peak_rss_gib()
    print(f"[rss] {label}: {peak:.3f} GiB", flush=True)
    if limit_gib > 0 and peak >= limit_gib:
        raise MemoryError(
            f"{label}: peak RSS {peak:.3f} GiB reached {limit_gib:.3f} GiB"
        )


def derangements(
    rng: np.random.Generator, draws: int, size: int
) -> np.ndarray:
    result = np.empty((draws, size), dtype=np.int16)
    base = np.arange(size, dtype=np.int16)
    for draw in range(draws):
        for _ in range(10_000):
            candidate = rng.permutation(base)
            if np.all(candidate != base):
                result[draw] = candidate
                break
        else:
            raise ContractError("failed to create derangement")
    return result


def create_plan() -> None:
    require(not PLAN_DIR.exists(), "plan already exists")
    validate_static_hashes()
    targets = read_lines(DEV_TARGETS)
    sealed = read_lines(SEALED_TEST_TARGETS)
    genes = read_lines(RESPONSE_GENES)
    require(len(targets) == DEV_COUNT, "development target count drifted")
    require(len(sealed) == 410, "sealed target count drifted")
    require(set(targets).isdisjoint(sealed), "development/test overlap")
    require(len(genes) == RESPONSE_COUNT, "response count drifted")
    pan_perts = read_lines(PAN_ROOT / "perts.txt")
    gallery = set(targets) | set(sealed)
    pan_center_rows = np.asarray(
        [index for index, target in enumerate(pan_perts) if target not in gallery],
        dtype=np.int16,
    )
    require(len(pan_center_rows) > DEV_COUNT, "Pan centering pool unexpectedly small")
    with np.load(PARENT140_PLAN / "plan-arrays.npz", allow_pickle=False) as data:
        folds = np.asarray(data["development_folds"], dtype=np.int8)
        blocks = np.asarray(data["development_blocks"], dtype=np.int16)
        bootstrap = np.asarray(data["bootstrap_block_draws"], dtype=np.uint16)
    require(folds.shape == (DEV_COUNT,), "fold shape drifted")
    require(set(folds.tolist()) == set(range(FOLDS)), "fold labels drifted")
    require(blocks.shape == (DEV_COUNT,), "block shape drifted")
    require(bootstrap.shape[0] == BOOTSTRAP_DRAWS, "bootstrap draws drifted")
    target_nulls = derangements(
        np.random.default_rng(TARGET_NULL_SEED), NULL_DRAWS, DEV_COUNT
    )
    coordinate_nulls = derangements(
        np.random.default_rng(COORD_NULL_SEED), NULL_DRAWS, RESPONSE_COUNT
    )
    PLAN_DIR.mkdir(parents=True)
    atomic_npz(
        PLAN_DIR / "plan-arrays.npz",
        development_folds=folds,
        development_blocks=blocks,
        bootstrap_block_draws=bootstrap,
        target_label_nulls=target_nulls,
        coordinate_nulls=coordinate_nulls,
        pan_center_rows=pan_center_rows,
    )
    manifest = {
        "schema_version": 1,
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "created_before_experiment_141_development_access": True,
        "development_targets": len(targets),
        "sealed_targets": len(sealed),
        "response_genes": len(genes),
        "fold_sizes": [int(np.sum(folds == fold)) for fold in range(FOLDS)],
        "target_blocks": int(blocks.max()) + 1,
        "bootstrap_draws": len(bootstrap),
        "target_label_null_draws": len(target_nulls),
        "coordinate_null_draws": len(coordinate_nulls),
        "pan_center_rows": len(pan_center_rows),
        "pan_center_row_index_sha256": hashlib.sha256(
            pan_center_rows.tobytes()
        ).hexdigest(),
        "seeds": {
            "target_label": TARGET_NULL_SEED,
            "coordinate": COORD_NULL_SEED,
        },
        "files": {
            "plan-arrays.npz": {
                "bytes": (PLAN_DIR / "plan-arrays.npz").stat().st_size,
                "sha256": sha256_file(PLAN_DIR / "plan-arrays.npz"),
            }
        },
    }
    atomic_json(PLAN_DIR / "plan-manifest.json", manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))


def validate_static_hashes() -> None:
    for path, expected in EXPECTED_SHA256.items():
        require(path.is_file(), f"missing input: {path}")
        require(sha256_file(path) == expected, f"input hash changed: {path}")


def validate_plan() -> dict[str, np.ndarray]:
    require(PLAN_DIR.is_dir(), "plan missing")
    manifest = read_json(PLAN_DIR / "plan-manifest.json")
    require(manifest["development_targets"] == DEV_COUNT, "development count")
    require(manifest["sealed_targets"] == 410, "sealed count")
    require(manifest["response_genes"] == RESPONSE_COUNT, "response count")
    expected = manifest["files"]["plan-arrays.npz"]
    path = PLAN_DIR / "plan-arrays.npz"
    require(path.stat().st_size == expected["bytes"], "plan bytes changed")
    require(sha256_file(path) == expected["sha256"], "plan hash changed")
    with np.load(path, allow_pickle=False) as data:
        arrays = {name: np.asarray(data[name]) for name in data.files}
    require(arrays["development_folds"].shape == (DEV_COUNT,), "fold shape")
    require(arrays["development_blocks"].shape == (DEV_COUNT,), "block shape")
    require(
        arrays["target_label_nulls"].shape == (NULL_DRAWS, DEV_COUNT),
        "target null shape",
    )
    require(
        arrays["coordinate_nulls"].shape == (NULL_DRAWS, RESPONSE_COUNT),
        "coordinate null shape",
    )
    require(arrays["pan_center_rows"].ndim == 1, "Pan center row shape")
    require(len(arrays["pan_center_rows"]) > DEV_COUNT, "Pan center row count")
    require(
        hashlib.sha256(arrays["pan_center_rows"].tobytes()).hexdigest()
        == manifest["pan_center_row_index_sha256"],
        "Pan center row hash changed",
    )
    return arrays


def freeze() -> None:
    require(not SOURCE_MANIFEST.exists(), "source manifest exists")
    require(not AUTHORIZATION.exists(), "development already authorized")
    validate_static_hashes()
    validate_plan()
    sources = {
        str(path.relative_to(ROOT)): {
            "bytes": path.stat().st_size,
            "sha256": expected,
        }
        for path, expected in EXPECTED_SHA256.items()
    }
    implementation = {
        str(path.relative_to(ROOT)): {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in (HERE / "README.md", PROTOCOL, RUNNER, TESTS)
    }
    manifest = {
        "schema_version": 1,
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "created_before_experiment_141_development_access": True,
        "sources": sources,
        "implementation": implementation,
        "plan_manifest_sha256": sha256_file(PLAN_DIR / "plan-manifest.json"),
        "sealed_410_effect_access_authorized": False,
    }
    atomic_json(SOURCE_MANIFEST, manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))


def authorize_development(confirm: str | None) -> None:
    require(confirm == DEVELOPMENT_CONFIRM, "explicit development confirmation missing")
    require(SOURCE_MANIFEST.is_file(), "source manifest missing")
    require(not AUTHORIZATION.exists(), "authorization exists")
    require(
        not ACCESS_LOCK.exists()
        and not RESULT_DIR.exists()
        and not DEVELOPMENT_FREEZE.exists()
        and not INFERENCE_DIR.exists()
        and not INFERENCE_FREEZE.exists(),
        "development attempted",
    )
    validate_static_hashes()
    validate_plan()
    manifest = read_json(SOURCE_MANIFEST)
    require(
        manifest["plan_manifest_sha256"]
        == sha256_file(PLAN_DIR / "plan-manifest.json"),
        "plan changed after freeze",
    )
    for relative, expected in manifest["implementation"].items():
        path = ROOT / relative
        require(path.stat().st_size == expected["bytes"], f"bytes changed: {relative}")
        require(sha256_file(path) == expected["sha256"], f"hash changed: {relative}")
    authorization = {
        "schema_version": 1,
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "explicit_post_review_confirmation": confirm,
        "development_613_effect_access_authorized": True,
        "conditional_nulls_authorized_only_after_raw_gate": True,
        "sealed_410_effect_access_authorized": False,
        "cell_eval2_authorized": False,
        "count_generation_authorized": False,
        "challenge_build_or_submission_authorized": False,
        "production_integration_authorized": False,
        "source_manifest_sha256": sha256_file(SOURCE_MANIFEST),
        "protocol_sha256": sha256_file(PROTOCOL),
        "plan_manifest_sha256": sha256_file(PLAN_DIR / "plan-manifest.json"),
    }
    atomic_json(AUTHORIZATION, authorization)
    print(json.dumps(authorization, indent=2, sort_keys=True))


def validate_authorization() -> None:
    validate_static_hashes()
    validate_plan()
    require(SOURCE_MANIFEST.is_file(), "source manifest missing")
    require(AUTHORIZATION.is_file(), "authorization missing")
    manifest = read_json(SOURCE_MANIFEST)
    for relative, expected in manifest["sources"].items():
        path = ROOT / relative
        require(path.stat().st_size == expected["bytes"], f"bytes changed: {relative}")
        require(sha256_file(path) == expected["sha256"], f"hash changed: {relative}")
    for relative, expected in manifest["implementation"].items():
        path = ROOT / relative
        require(path.stat().st_size == expected["bytes"], f"bytes changed: {relative}")
        require(sha256_file(path) == expected["sha256"], f"hash changed: {relative}")
    auth = read_json(AUTHORIZATION)
    require(auth["development_613_effect_access_authorized"] is True, "development denied")
    require(auth["sealed_410_effect_access_authorized"] is False, "sealed test enabled")
    require(auth["cell_eval2_authorized"] is False, "cell-eval2 enabled")
    require(
        auth["challenge_build_or_submission_authorized"] is False,
        "challenge action enabled",
    )
    require(
        auth["source_manifest_sha256"] == sha256_file(SOURCE_MANIFEST),
        "source manifest changed",
    )
    require(auth["protocol_sha256"] == sha256_file(PROTOCOL), "protocol changed")
    require(
        auth["plan_manifest_sha256"]
        == sha256_file(PLAN_DIR / "plan-manifest.json"),
        "plan changed after authorization",
    )


def validate_directory_manifest(directory: Path) -> None:
    manifest_path = directory / "artifact-manifest.json"
    manifest = read_json(manifest_path)
    for relative, expected in manifest["files"].items():
        path = directory / relative
        expected_bytes = expected.get("bytes", expected.get("size_bytes"))
        require(expected_bytes is not None, f"manifest size missing: {path}")
        require(path.is_file(), f"manifest member missing: {path}")
        require(path.stat().st_size == expected_bytes, f"member bytes changed: {path}")
        require(sha256_file(path) == expected["sha256"], f"member hash changed: {path}")


def load_modules() -> tuple[Any, Any]:
    return (
        import_path("exp141_parent132", PARENT132_RUN),
        import_path("exp141_parent140", PARENT140_RUN),
    )


def load_baselines() -> dict[str, np.ndarray]:
    baselines = {
        "k562": np.asarray(
            np.load(PARENT140_IDENTITY / "baseline-k562.npy"), dtype=np.float64
        ),
        "rpe1": np.asarray(
            np.load(PARENT140_IDENTITY / "baseline-non_k562.npy"),
            dtype=np.float64,
        ),
        "jurkat": np.asarray(
            np.load(PARENT140_IDENTITY / "baseline-non_k562.npy"),
            dtype=np.float64,
        ),
    }
    for family, values in baselines.items():
        require(values.shape == (DEV_COUNT, RESPONSE_COUNT), f"{family}: baseline shape")
        require(np.all(np.isfinite(values)), f"{family}: baseline nonfinite")
    return baselines


def apply_operator(matrix: np.ndarray, path: Path) -> np.ndarray:
    with np.load(path, allow_pickle=False) as data:
        evals = np.asarray(data["evals"], dtype=np.float64)
        evecs = np.asarray(data["evecs"], dtype=np.float64)
    require(evals.shape == (RESPONSE_COUNT,), "operator eigenvalue shape")
    require(evecs.shape == (RESPONSE_COUNT, RESPONSE_COUNT), "operator eigenvector shape")
    weights = evals ** (-0.65 / 2.0)
    values = np.asarray(matrix, dtype=np.float64)
    return ((values @ evecs) * weights) @ evecs.T


def load_pan_geometry(
    parent132: Any,
    baselines: Mapping[str, np.ndarray],
    plan_arrays: Mapping[str, np.ndarray],
    limit_gib: float,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    targets = read_lines(DEV_TARGETS)
    responses = read_lines(RESPONSE_GENES)
    bank = parent132.Bank.from_root("pan_h1", PAN_ROOT)
    selected = np.asarray(
        parent132.load_effect_rows(bank, targets, responses), dtype=np.float64
    )
    row_index = {target: row for row, target in enumerate(bank.perts)}
    column_index = {gene: column for column, gene in enumerate(bank.genes)}
    rows = np.asarray([row_index[target] for target in targets], dtype=np.int64)
    columns = np.asarray([column_index[gene] for gene in responses], dtype=np.int64)
    measured = np.load(PAN_ROOT / "measured_mask.npy", mmap_mode="r")
    require(
        np.all(np.asarray(measured[np.ix_(rows, columns)], dtype=bool)),
        "Pan selected support is incomplete",
    )
    delta = np.load(PAN_ROOT / "delta.npy", mmap_mode="r")
    center_rows = np.asarray(plan_arrays["pan_center_rows"], dtype=np.int64)
    require(np.all((center_rows >= 0) & (center_rows < len(bank.perts))), "Pan center row bounds")
    gallery_rows = set(rows.tolist())
    sealed_names = set(read_lines(SEALED_TEST_TARGETS))
    sealed_rows = {
        row_index[target] for target in sealed_names
    }
    require(gallery_rows.isdisjoint(set(center_rows.tolist())), "development row in Pan center")
    require(sealed_rows.isdisjoint(set(center_rows.tolist())), "sealed row in Pan center")
    mean = np.zeros(RESPONSE_COUNT, dtype=np.float64)
    chunk = 256
    for lo in range(0, len(center_rows), chunk):
        selected_rows = center_rows[lo : lo + chunk]
        mean += np.asarray(delta[selected_rows][:, columns], dtype=np.float64).sum(axis=0)
    mean /= len(center_rows)
    centered = selected - mean
    memory_check("Pan centered", limit_gib)

    transformed_by_family = {}
    for key, families in (
        ("k562", ("k562",)),
        ("non_k562", ("rpe1", "jurkat")),
    ):
        transformed = apply_operator(
            centered, PARENT140_IDENTITY / f"whitener-{key}.npz"
        )
        for family in families:
            transformed_by_family[family] = transformed.copy()
        memory_check(f"Pan transformed {key}", limit_gib)
    matched = norm_match_pan(transformed_by_family, baselines)
    return matched, transformed_by_family


def norm_match_pan(
    transformed: Mapping[str, np.ndarray],
    baselines: Mapping[str, np.ndarray],
    *,
    row_permutation: np.ndarray | None = None,
    coordinate_permutation: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    result = {}
    for family in FAMILIES:
        values = np.asarray(transformed[family], dtype=np.float64)
        if row_permutation is not None:
            permutation = np.asarray(row_permutation, dtype=np.int64)
            require(permutation.shape == (DEV_COUNT,), "Pan row permutation shape")
            values = values[permutation]
        if coordinate_permutation is not None:
            permutation = np.asarray(coordinate_permutation, dtype=np.int64)
            require(
                permutation.shape == (RESPONSE_COUNT,),
                "Pan coordinate permutation shape",
            )
            values = values[:, permutation]
        base = baselines[family]
        base_norm = np.linalg.norm(base, axis=1)
        pan_norm = np.linalg.norm(values, axis=1)
        valid = (base_norm > EPS) & (pan_norm > EPS)
        require(np.all(valid), f"{family}: degenerate Pan or baseline row")
        scaled = values * (base_norm / pan_norm)[:, None]
        require(
            np.max(np.abs(np.linalg.norm(scaled, axis=1) - base_norm)) < 1e-8,
            f"{family}: Pan norm matching failed",
        )
        result[family] = scaled
    return result


def truth_by_family(truth: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    return {
        "k562": 0.5
        * (
            np.asarray(truth["k562_gwps"], dtype=np.float64)
            + np.asarray(truth["k562_essential"], dtype=np.float64)
        ),
        "rpe1": np.asarray(truth["rpe1"], dtype=np.float64),
        "jurkat": np.asarray(truth["nadig_jurkat"], dtype=np.float64),
    }


def fit_calibration(
    destination_family: str,
    baselines: Mapping[str, np.ndarray],
    pan: Mapping[str, np.ndarray],
    family_truth: Mapping[str, np.ndarray],
    training_families: Sequence[str],
    training_indices: np.ndarray,
    rho: float,
) -> np.ndarray:
    indices = np.asarray(training_indices, dtype=np.int64)
    numerator = np.zeros(RESPONSE_COUNT, dtype=np.float64)
    denominator = np.zeros(RESPONSE_COUNT, dtype=np.float64)
    q = pan[destination_family][indices] - baselines[destination_family][indices]
    for family in training_families:
        residual = (
            family_truth[family][indices]
            - baselines[destination_family][indices]
        )
        numerator += np.sum(q * residual, axis=0)
        denominator += np.sum(q * q, axis=0)
    total_denominator = float(denominator.sum())
    pooled = (
        float(numerator.sum()) / total_denominator
        if total_denominator > EPS
        else 0.0
    )
    diagonal = np.full(RESPONSE_COUNT, pooled, dtype=np.float64)
    np.divide(
        numerator,
        denominator,
        out=diagonal,
        where=denominator > EPS,
    )
    coefficient = (1.0 - rho) * pooled + rho * diagonal
    return np.clip(coefficient, 0.0, 1.0)


def candidate_prediction(
    name: str,
    held_family: str,
    training_indices: np.ndarray,
    query_indices: np.ndarray,
    baselines: Mapping[str, np.ndarray],
    pan: Mapping[str, np.ndarray],
    family_truth: Mapping[str, np.ndarray],
) -> np.ndarray:
    base = baselines[held_family]
    if name == "identity":
        return base
    rho = CANDIDATES[name]
    require(rho is not None, f"{name}: missing rho")
    teachers = CALIBRATION_TEACHERS[held_family]
    coefficient = fit_calibration(
        held_family,
        baselines,
        pan,
        family_truth,
        teachers,
        training_indices,
        rho,
    )
    output = base.copy()
    query = np.asarray(query_indices, dtype=np.int64)
    q = pan[held_family][query] - base[query]
    output[query] = base[query] + q * coefficient
    return output


def score_queries(
    predictions: Mapping[str, np.ndarray],
    baselines: Mapping[str, np.ndarray],
    truth: Mapping[str, np.ndarray],
    query: np.ndarray,
    gallery: np.ndarray,
    parent140: Any,
) -> dict[str, dict[str, np.ndarray]]:
    result = {}
    for dataset in SCORED_DATASETS:
        family = DATASET_FAMILY[dataset]
        result[dataset] = {
            "pds": parent140.pds_subset(
                predictions[family], truth[dataset], query, gallery
            )
            - parent140.pds_subset(
                baselines[family], truth[dataset], query, gallery
            ),
            "signed_credit": parent140.signed_credit(
                predictions[family][query], truth[dataset][query]
            )
            - parent140.signed_credit(
                baselines[family][query], truth[dataset][query]
            ),
        }
    return result


def empty_scores(size: int) -> dict[str, dict[str, np.ndarray]]:
    return {
        dataset: {
            metric: np.empty(size, dtype=np.float64)
            for metric in ("pds", "signed_credit")
        }
        for dataset in SCORED_DATASETS
    }


def aggregate_scores(
    scores: Mapping[str, Mapping[str, np.ndarray]]
) -> dict[str, Any]:
    datasets = {
        dataset: {
            metric: float(np.mean(scores[dataset][metric]))
            for metric in ("pds", "signed_credit")
        }
        for dataset in SCORED_DATASETS
    }
    families = {"pds": {}, "signed_credit": {}}
    per_target = {}
    macro = {}
    for metric in ("pds", "signed_credit"):
        family_values = {
            family: np.mean(
                np.stack(
                    [scores[dataset][metric] for dataset in FAMILY_DATASETS[family]]
                ),
                axis=0,
            )
            for family in FAMILIES
        }
        families[metric] = {
            family: float(values.mean())
            for family, values in family_values.items()
        }
        per_target[metric] = np.mean(
            np.stack([family_values[family] for family in FAMILIES]), axis=0
        )
        macro[metric] = float(per_target[metric].mean())
    return {
        "datasets": datasets,
        "families": families,
        "per_target": per_target,
        "macro": macro,
    }


def evaluate_procedure(
    folds: np.ndarray,
    baselines: Mapping[str, np.ndarray],
    pan: Mapping[str, np.ndarray],
    truth: Mapping[str, np.ndarray],
    parent140: Any,
    candidate_names: Sequence[str] = tuple(CANDIDATES),
) -> dict[str, Any]:
    family_truth = truth_by_family(truth)
    all_indices = np.arange(DEV_COUNT, dtype=np.int64)
    scores_by_candidate = {
        name: empty_scores(DEV_COUNT)
        for name in candidate_names
    }
    fold_records = []
    for fold in range(FOLDS):
        query = all_indices[folds == fold]
        training = all_indices[folds != fold]
        for name in candidate_names:
            predictions = {
                family: candidate_prediction(
                    name,
                    family,
                    training,
                    query,
                    baselines,
                    pan,
                    family_truth,
                )
                for family in FAMILIES
            }
            score = score_queries(
                predictions,
                baselines,
                truth,
                query,
                all_indices,
                parent140,
            )
            for dataset in SCORED_DATASETS:
                for metric in ("pds", "signed_credit"):
                    scores_by_candidate[name][dataset][metric][query] = score[dataset][metric]
        fold_records.append(
            {
                "fold": fold,
                "targets": len(query),
                "training_targets": len(training),
                "candidates": list(candidate_names),
            }
        )
        print(f"[outer fold {fold}] scored {len(candidate_names)} candidates", flush=True)
    summaries = {
        name: aggregate_scores(scores)
        for name, scores in scores_by_candidate.items()
    }
    return {
        "fold_records": fold_records,
        "scores": scores_by_candidate,
        "summaries": summaries,
    }


def raw_gate(result: Mapping[str, Any]) -> dict[str, bool]:
    primary = result["summaries"][PRIMARY_CANDIDATE]
    scalar = result["summaries"]["scalar_control"]
    return {
        "primary_pds_families_nonnegative": all(
            value >= -EPS for value in primary["families"]["pds"].values()
        ),
        "primary_signed_credit_families_nonnegative": all(
            value >= -EPS
            for value in primary["families"]["signed_credit"].values()
        ),
        "primary_macro_pds_at_least_0.003": primary["macro"]["pds"] >= MIN_GAIN,
        "primary_macro_signed_credit_at_least_0.003": (
            primary["macro"]["signed_credit"] >= MIN_GAIN
        ),
        "primary_beats_scalar_by_0.001_pds": (
            primary["macro"]["pds"] - scalar["macro"]["pds"]
            >= MIN_DIAG_OVER_SCALAR
        ),
    }


def output_manifest(directory: Path) -> dict[str, Any]:
    return {
        "files": {
            path.name: {
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in sorted(directory.iterdir())
            if path.is_file() and path.name != "artifact-manifest.json"
        }
    }


def load_observed_scores() -> dict[str, dict[str, np.ndarray]]:
    with np.load(RESULT_DIR / "per-target-observed.npz", allow_pickle=False) as data:
        primary = empty_scores(DEV_COUNT)
        for dataset in SCORED_DATASETS:
            for metric in ("pds", "signed_credit"):
                primary[dataset][metric] = np.asarray(
                    data[f"crossfit__{PRIMARY_CANDIDATE}__{dataset}__{metric}"],
                    dtype=np.float64,
                )
    return primary


def validate_development_result() -> dict[str, Any]:
    require(RESULT_DIR.is_dir(), "development result missing")
    require(DEVELOPMENT_FREEZE.is_file(), "development result freeze missing")
    freeze_value = read_json(DEVELOPMENT_FREEZE)
    require(
        freeze_value["authorization_sha256"] == sha256_file(AUTHORIZATION),
        "development authorization binding changed",
    )
    require(
        freeze_value["artifact_manifest_sha256"]
        == sha256_file(RESULT_DIR / "artifact-manifest.json"),
        "development manifest binding changed",
    )
    validate_directory_manifest(RESULT_DIR)
    require(
        freeze_value["result_sha256"] == sha256_file(RESULT_DIR / "result.json"),
        "development result binding changed",
    )
    return read_json(RESULT_DIR / "result.json")


def block_bootstrap(
    values: np.ndarray, labels: np.ndarray, draws: np.ndarray
) -> np.ndarray:
    labels = np.asarray(labels, dtype=np.int64)
    samples = np.asarray(draws, dtype=np.int64)
    block_count = int(labels.max()) + 1
    sums = np.bincount(
        labels,
        weights=np.asarray(values, dtype=np.float64),
        minlength=block_count,
    )
    sizes = np.bincount(labels, minlength=block_count).astype(np.float64)
    return sums[samples].sum(axis=1) / sizes[samples].sum(axis=1)


def positive_share(values: Sequence[float]) -> float:
    positive = np.maximum(np.asarray(values, dtype=np.float64), 0.0)
    return 1.0 if positive.sum() <= EPS else float(positive.max() / positive.sum())


def concentration_report(
    scores: Mapping[str, Mapping[str, np.ndarray]],
    blocks: np.ndarray,
) -> dict[str, Any]:
    summary = aggregate_scores(scores)
    result = {}
    labels = np.asarray(blocks, dtype=np.int64)
    for metric in ("pds", "signed_credit"):
        block_values = np.bincount(
            labels,
            weights=np.asarray(summary["per_target"][metric], dtype=np.float64),
        )
        result[metric] = {
            "family": positive_share(
                tuple(summary["families"][metric][family] for family in FAMILIES)
            ),
            "target_block": positive_share(tuple(block_values.tolist())),
        }
    return result


def run_inference(limit_gib: float) -> None:
    validate_authorization()
    require(not INFERENCE_DIR.exists(), "conditional inference exists")
    require(not INFERENCE_FREEZE.exists(), "conditional inference freeze exists")
    observed_result = validate_development_result()
    require(observed_result["raw_gate_passes"] is True, "raw development gate failed")
    require(
        read_json(DEVELOPMENT_FREEZE)["conditional_inference_authorized"] is True,
        "conditional inference not authorized",
    )
    arrays = validate_plan()
    parent132, parent140 = load_modules()
    validate_directory_manifest(PARENT140_IDENTITY)
    validate_directory_manifest(PAN_ROOT)
    parent140.validate_identity()
    parent140.validate_source_hashes(parent140.development_source_paths())
    baselines = load_baselines()
    _, transformed_pan = load_pan_geometry(
        parent132, baselines, arrays, limit_gib
    )
    truth = parent140.load_development_truth(limit_gib)

    null_values = {
        kind: {
            metric: np.empty(NULL_DRAWS, dtype=np.float64)
            for metric in ("pds", "signed_credit")
        }
        for kind in ("target_label", "coordinate")
    }
    for kind in ("target_label", "coordinate"):
        for draw in range(NULL_DRAWS):
            if kind == "target_label":
                pan = norm_match_pan(
                    transformed_pan,
                    baselines,
                    row_permutation=arrays["target_label_nulls"][draw],
                )
            else:
                pan = norm_match_pan(
                    transformed_pan,
                    baselines,
                    coordinate_permutation=arrays["coordinate_nulls"][draw],
                )
            null_result = evaluate_procedure(
                arrays["development_folds"],
                baselines,
                pan,
                truth,
                parent140,
                candidate_names=(PRIMARY_CANDIDATE,),
            )
            observed_null = null_result["summaries"][PRIMARY_CANDIDATE]["macro"]
            for metric in ("pds", "signed_credit"):
                null_values[kind][metric][draw] = observed_null[metric]
            if draw % 5 == 4:
                print(f"[{kind}] {draw + 1}/{NULL_DRAWS}", flush=True)
                memory_check(f"{kind} null", limit_gib)

    primary_scores = load_observed_scores()
    observed_summary = aggregate_scores(primary_scores)
    bootstrap = {}
    bootstrap_values = {}
    for metric in ("pds", "signed_credit"):
        samples = block_bootstrap(
            observed_summary["per_target"][metric],
            arrays["development_blocks"],
            arrays["bootstrap_block_draws"],
        )
        bootstrap_values[metric] = samples
        bootstrap[metric] = {
            "draws": len(samples),
            "lower_95": float(np.quantile(samples, 0.025, method="linear")),
            "upper_95": float(np.quantile(samples, 0.975, method="linear")),
        }
    concentration = concentration_report(
        primary_scores, arrays["development_blocks"]
    )

    null_report = {
        kind: {
            metric: {
                "draws": len(null_values[kind][metric]),
                "p95": float(
                    np.quantile(
                        null_values[kind][metric],
                        0.95,
                        method="linear",
                    )
                ),
                "add_one_p": float(
                    (
                        1
                        + np.sum(
                            null_values[kind][metric]
                            >= observed_summary["macro"][metric]
                        )
                    )
                    / (NULL_DRAWS + 1)
                ),
            }
            for metric in ("pds", "signed_credit")
        }
        for kind in ("target_label", "coordinate")
    }
    checks = {}
    for metric in ("pds", "signed_credit"):
        checks[f"{metric}_bootstrap_lower_positive"] = (
            bootstrap[metric]["lower_95"] > 0.0
        )
        for kind in ("target_label", "coordinate"):
            checks[f"{metric}_above_{kind}_null"] = (
                observed_summary["macro"][metric]
                > null_report[kind][metric]["p95"]
            )
        checks[f"{metric}_family_concentration"] = (
            concentration[metric]["family"] <= MAX_CONCENTRATION
        )
        checks[f"{metric}_target_block_concentration"] = (
            concentration[metric]["target_block"] <= MAX_CONCENTRATION
        )
    serializable = {
        "schema_version": 1,
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "stage": "conditional_inference_after_raw_gate",
        "primary_candidate": PRIMARY_CANDIDATE,
        "observed": {
            key: value
            for key, value in observed_summary.items()
            if key != "per_target"
        },
        "bootstrap": bootstrap,
        "nulls": null_report,
        "positive_gain_concentration": concentration,
        "checks": checks,
        "passes": all(checks.values()),
        "decision": (
            "PASS_FOR_EXACT_A15_DEVELOPMENT_FOLLOWUP_ONLY"
            if all(checks.values())
            else "NO-GO"
        ),
        "untouched_410_test_accessed": False,
        "production_integration_authorized": False,
        "peak_rss_gib": peak_rss_gib(),
    }
    temporary = INFERENCE_DIR.with_name(f".{INFERENCE_DIR.name}.tmp")
    require(not temporary.exists(), "temporary inference directory exists")
    temporary.mkdir(parents=True)
    atomic_json(temporary / "result.json", serializable)
    saved = {}
    for kind, metrics in null_values.items():
        for metric, values in metrics.items():
            saved[f"null__{kind}__{metric}"] = values
    for metric, values in bootstrap_values.items():
        saved[f"bootstrap__{metric}"] = values
    atomic_npz(temporary / "nulls-and-bootstrap.npz", **saved)
    atomic_json(temporary / "artifact-manifest.json", output_manifest(temporary))
    os.replace(temporary, INFERENCE_DIR)
    atomic_json(
        INFERENCE_FREEZE,
        {
            "schema_version": 1,
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "development_result_freeze_sha256": sha256_file(DEVELOPMENT_FREEZE),
            "artifact_manifest_sha256": sha256_file(
                INFERENCE_DIR / "artifact-manifest.json"
            ),
            "result_sha256": sha256_file(INFERENCE_DIR / "result.json"),
            "sealed_410_test_authorized": False,
            "production_integration_authorized": False,
        },
    )
    print(json.dumps(serializable, indent=2, sort_keys=True))


def develop(limit_gib: float) -> None:
    validate_authorization()
    require(not ACCESS_LOCK.exists(), "development already accessed")
    require(not RESULT_DIR.exists(), "development result exists")
    require(not DEVELOPMENT_FREEZE.exists(), "development result freeze exists")
    require(not INFERENCE_DIR.exists(), "conditional inference already exists")
    require(not INFERENCE_FREEZE.exists(), "conditional inference freeze exists")
    arrays = validate_plan()
    parent132, parent140 = load_modules()
    validate_directory_manifest(PARENT140_IDENTITY)
    validate_directory_manifest(PAN_ROOT)
    parent140.validate_identity()
    parent140.validate_source_hashes(parent140.development_source_paths())
    atomic_json(
        ACCESS_LOCK,
        {
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "development_613_effects_accessed": True,
            "sealed_410_effects_accessed": False,
            "authorization_sha256": sha256_file(AUTHORIZATION),
        },
    )
    baselines = load_baselines()
    pan, _ = load_pan_geometry(parent132, baselines, arrays, limit_gib)
    truth = parent140.load_development_truth(limit_gib)
    result = evaluate_procedure(
        arrays["development_folds"],
        baselines,
        pan,
        truth,
        parent140,
    )
    checks = raw_gate(result)
    serializable = {
        "schema_version": 1,
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "stage": "source_only_613_target_development",
        "candidate_grid": list(CANDIDATES),
        "primary_candidate": PRIMARY_CANDIDATE,
        "cross_fitted_folds": result["fold_records"],
        "cross_fitted_candidates": {
            name: {
                key: value
                for key, value in summary.items()
                if key != "per_target"
            }
            for name, summary in result["summaries"].items()
        },
        "raw_gate": checks,
        "raw_gate_passes": all(checks.values()),
        "decision": (
            "PASS_TO_CONDITIONAL_INFERENCE_ONLY"
            if all(checks.values())
            else "NO-GO"
        ),
        "untouched_410_test_accessed": False,
        "cell_eval2_run": False,
        "counts_or_submission_built": False,
        "peak_rss_gib": peak_rss_gib(),
    }
    temporary = RESULT_DIR.with_name(f".{RESULT_DIR.name}.tmp")
    require(not temporary.exists(), "temporary result directory exists")
    temporary.mkdir(parents=True)
    atomic_json(temporary / "result.json", serializable)
    saved = {}
    for name, scores in result["scores"].items():
        for dataset in SCORED_DATASETS:
            for metric in ("pds", "signed_credit"):
                saved[f"crossfit__{name}__{dataset}__{metric}"] = scores[dataset][metric]
    atomic_npz(temporary / "per-target-observed.npz", **saved)
    atomic_json(temporary / "artifact-manifest.json", output_manifest(temporary))
    os.replace(temporary, RESULT_DIR)
    atomic_json(
        DEVELOPMENT_FREEZE,
        {
            "schema_version": 1,
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "authorization_sha256": sha256_file(AUTHORIZATION),
            "artifact_manifest_sha256": sha256_file(
                RESULT_DIR / "artifact-manifest.json"
            ),
            "result_sha256": sha256_file(RESULT_DIR / "result.json"),
            "conditional_inference_authorized": serializable[
                "raw_gate_passes"
            ],
            "sealed_410_test_authorized": False,
        },
    )
    print(json.dumps(serializable, indent=2, sort_keys=True))


def audit() -> None:
    payload: dict[str, Any] = {
        "plan_exists": PLAN_DIR.exists(),
        "source_manifest_exists": SOURCE_MANIFEST.exists(),
        "authorization_exists": AUTHORIZATION.exists(),
        "development_accessed": ACCESS_LOCK.exists(),
        "development_result_exists": RESULT_DIR.exists(),
        "development_result_freeze_exists": DEVELOPMENT_FREEZE.exists(),
        "conditional_inference_exists": INFERENCE_DIR.exists(),
        "conditional_inference_freeze_exists": INFERENCE_FREEZE.exists(),
        "sealed_410_code_path_exists": False,
    }
    if AUTHORIZATION.exists():
        validate_authorization()
    if RESULT_DIR.exists() or DEVELOPMENT_FREEZE.exists():
        payload["development"] = validate_development_result()
    if INFERENCE_DIR.exists() or INFERENCE_FREEZE.exists():
        require(INFERENCE_DIR.is_dir(), "inference result directory missing")
        require(INFERENCE_FREEZE.is_file(), "inference freeze missing")
        freeze_value = read_json(INFERENCE_FREEZE)
        require(
            freeze_value["development_result_freeze_sha256"]
            == sha256_file(DEVELOPMENT_FREEZE),
            "inference development binding changed",
        )
        require(
            freeze_value["artifact_manifest_sha256"]
            == sha256_file(INFERENCE_DIR / "artifact-manifest.json"),
            "inference manifest binding changed",
        )
        require(
            freeze_value["result_sha256"]
            == sha256_file(INFERENCE_DIR / "result.json"),
            "inference result binding changed",
        )
        validate_directory_manifest(INFERENCE_DIR)
        payload["inference"] = read_json(INFERENCE_DIR / "result.json")
    print(json.dumps(payload, indent=2, sort_keys=True))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "command",
        choices=(
            "plan",
            "freeze",
            "authorize-development",
            "develop",
            "inference",
            "audit",
        ),
    )
    parser.add_argument("--max-rss-gib", type=float, default=16.0)
    parser.add_argument("--confirm")
    args = parser.parse_args()
    if args.command == "plan":
        create_plan()
    elif args.command == "freeze":
        freeze()
    elif args.command == "authorize-development":
        authorize_development(args.confirm)
    elif args.command == "develop":
        develop(args.max_rss_gib)
    elif args.command == "inference":
        run_inference(args.max_rss_gib)
    else:
        audit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

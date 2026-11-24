#!/usr/bin/env python3
"""Experiment 138: DLD1 replicated sparse-sign gating.

Stages:

* metadata: axes/manifests only; no effect arrays.
* plan: metadata and authenticated DepMap features only; no effect arrays.
* freeze: hash and bind every input before Experiment 138 outcome access.
* develop: read only the inherited 613 development rows.
* test: read the inherited 410 test rows only after every development gate.

All writes are confined to this experiment directory.
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
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL = HERE / "preregistered_protocol.json"
README = HERE / "README.md"
TESTS = HERE / "test_run.py"
PLAN_DIR = HERE / "work" / "preregistered-plan-v1"
AUTHORIZATION = HERE / "PREVALUE_AUTHORIZATION.v1.json"
DEVELOPMENT_LOCK = HERE / "work" / "development-access.v1.lock"
DEVELOPMENT_DIR = HERE / "work" / "development-result-v1"
CANDIDATE_FREEZE = HERE / "DEVELOPMENT_CANDIDATE_FREEZE.v1.json"
TEST_AUTHORIZATION = HERE / "TEST_ACCESS_AUTHORIZATION.v1.json"
TEST_LOCK = HERE / "work" / "untouched-test-access.v1.lock"
TEST_DIR = HERE / "work" / "untouched-test-result-v1"

PARENT_DIR = ROOT / "experiments/132-response-decomposition-ridge"
PARENT_RUN = PARENT_DIR / "run.py"
PARENT_AUTHORIZATION = PARENT_DIR / "SOURCE_ONLY_AUTHORIZATION.v2.json"
PARENT_PLAN_DIR = PARENT_DIR / "work" / "preregistered-freeze-v2"
PARENT_PLAN_MANIFEST = PARENT_PLAN_DIR / "plan-manifest.json"
PARENT_ALL_TARGETS = PARENT_PLAN_DIR / "all_targets.txt"
PARENT_SELECTION_TARGETS = PARENT_PLAN_DIR / "selection_targets.txt"
PARENT_TEST_TARGETS = PARENT_PLAN_DIR / "test_targets.txt"
PARENT_RESPONSE_GENES = PARENT_PLAN_DIR / "response_genes.txt"
PARENT_RANDOMNESS = PARENT_PLAN_DIR / "frozen-randomness.npz"

DLD1_DIR = ROOT / "experiments/116-dld1-lowmoi/work/compact-source-v1"
DLD1_MANIFEST = DLD1_DIR / "artifact-manifest.json"
DLD1_TARGETS = DLD1_DIR / "target_symbols.txt"
DLD1_GENES = DLD1_DIR / "canonical_genes.txt"
DLD1_LOW1 = DLD1_DIR / "transfer_lfc_low1.npy"
DLD1_LOW2 = DLD1_DIR / "transfer_lfc_low2.npy"
DLD1_MASK1 = DLD1_DIR / "transfer_measured_low1.npy"
DLD1_MASK2 = DLD1_DIR / "transfer_measured_low2.npy"

EXPECTED_PARENT_RUN_SHA256 = (
    "e78e079fc3b4a715e8b412b0a48bfe44a2f3f6cfcadb8a65a789351c53877778"
)
EXPECTED_PARENT_AUTH_SHA256 = (
    "9bca6e620359e44b5440d0b8588412f40f274304ae67d9f55fd7e7d03341e197"
)
EXPECTED_PARENT_PLAN_SHA256 = (
    "170d55ed5f812eb85f5503af4bb9d7dec8d17dce6e58fe966d85bf0eec5185d1"
)
EXPECTED_ALL_TARGETS_FILE_SHA256 = (
    "f5a9d65cc1d01765dc5c61649fb876c9e2ddf4f316791b28439324129fea4d51"
)
EXPECTED_SELECTION_FILE_SHA256 = (
    "74a070a538455ffc2a528dae3094d21e047905463779b2a5206d681f5ca6ebb5"
)
EXPECTED_TEST_FILE_SHA256 = (
    "bd37bcacb4eea90da11902bbb8f4736068a7ff67c001739e3407de2b4bf73fc6"
)
EXPECTED_RESPONSE_FILE_SHA256 = (
    "24743754eead5b8ddd81634bd12feac1df431eb2b9581197755b75c4c90949c8"
)
EXPECTED_PANEL_FILE_SHA256 = (
    "f57edd7b912ebd718efc7ee9d0f334772513e7cc418d133ce525470e373b3276"
)
EXPECTED_DLD1_MANIFEST_SHA256 = (
    "c592e2ec51cb56e7df085224b3ff75d8735ce87280683fd0967135472e22be79"
)
EXPECTED_DLD1_TARGETS_SHA256 = (
    "db4847da3fd2312d66bc7e4808dcea597ed56af3fb2dee4bb91d8f25a801867d"
)
EXPECTED_DLD1_GENES_SHA256 = (
    "c2b8dec58c02810d64018164155b655c4e7c7be451a3ae579965db193af33723"
)

DEVELOPMENT_COUNT = 613
TEST_COUNT = 410
RESPONSE_COUNT = 4_992
DEVELOPMENT_SUPPORTED = 140
TEST_SUPPORTED = 94
FOLD_COUNT = 5
NULL_DRAWS = 49
BOOTSTRAP_DRAWS = 20_000
TARGET_NULL_SEED = 138_202_610
RESPONSE_NULL_SEED = 138_202_611
BOOTSTRAP_SEED = 138_202_612
MIN_MACRO_GAIN = 0.003
MAX_CONCENTRATION = 0.5
EPS = 1e-12
GIB = 1024**3
DEVELOP_CONFIRMATION = "I_CONFIRM_EXPERIMENT_138_DEVELOPMENT_ONLY"
TEST_AUTHORIZATION_CONFIRMATION = (
    "I_REVIEWED_EXPERIMENT_138_DEVELOPMENT_AND_AUTHORIZE_TEST"
)
TEST_CONFIRMATION = "I_CONFIRM_EXPERIMENT_138_UNTOUCHED_TEST"


class ContractError(RuntimeError):
    """Raised when an immutable experiment boundary is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def load_parent_module() -> Any:
    spec = importlib.util.spec_from_file_location("experiment132_parent", PARENT_RUN)
    require(spec is not None and spec.loader is not None, "cannot load Experiment 132")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


PARENT = load_parent_module()


@dataclass(frozen=True)
class Arm:
    name: str
    k: int | None
    beta: float

    @property
    def complexity(self) -> tuple[int, int, float]:
        if self.k is None:
            return (0, 0, 0.0)
        return (1, self.k, self.beta)


ARMS = (
    Arm("identity", None, 0.0),
    Arm("k050_beta0p025", 50, 0.025),
    Arm("k050_beta0p050", 50, 0.05),
    Arm("k100_beta0p025", 100, 0.025),
    Arm("k100_beta0p050", 100, 0.05),
    Arm("k200_beta0p025", 200, 0.025),
    Arm("k200_beta0p050", 200, 0.05),
)
ARM_BY_NAME = {arm.name: arm for arm in ARMS}


@dataclass(frozen=True)
class Proposal:
    low1: np.ndarray
    low2: np.ndarray
    measured1: np.ndarray
    measured2: np.ndarray
    supported: np.ndarray


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if os.uname().sysname == "Darwin" else value / 1024**2


def memory_checkpoint(stage: str, limit_gib: float) -> None:
    peak = peak_rss_gib()
    print(f"[rss] {stage}: {peak:.3f} GiB", flush=True)
    if limit_gib > 0 and peak >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak:.3f} GiB reached {limit_gib:.3f} GiB"
        )


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode("utf-8")).hexdigest()


def read_lines(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    )
    require(bool(values), f"{path}: empty axis")
    if unique:
        require(len(values) == len(set(values)), f"{path}: duplicate axis")
    return values


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def atomic_write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def atomic_savez(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".npz"
    )
    os.close(descriptor)
    try:
        np.savez_compressed(temporary, **arrays)
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def output_manifest(
    directory: Path, *, exclude: set[str] | None = None
) -> dict[str, Any]:
    excluded = exclude or set()
    return {
        "files": {
            path.name: {
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in sorted(directory.iterdir())
            if path.is_file() and path.name not in excluded
        }
    }


def validate_static_metadata() -> dict[str, Any]:
    expected_hashes = {
        PARENT_RUN: EXPECTED_PARENT_RUN_SHA256,
        PARENT_AUTHORIZATION: EXPECTED_PARENT_AUTH_SHA256,
        PARENT_PLAN_MANIFEST: EXPECTED_PARENT_PLAN_SHA256,
        PARENT_ALL_TARGETS: EXPECTED_ALL_TARGETS_FILE_SHA256,
        PARENT_SELECTION_TARGETS: EXPECTED_SELECTION_FILE_SHA256,
        PARENT_TEST_TARGETS: EXPECTED_TEST_FILE_SHA256,
        PARENT_RESPONSE_GENES: EXPECTED_RESPONSE_FILE_SHA256,
        PARENT.PANEL_PATH: EXPECTED_PANEL_FILE_SHA256,
        DLD1_MANIFEST: EXPECTED_DLD1_MANIFEST_SHA256,
        DLD1_TARGETS: EXPECTED_DLD1_TARGETS_SHA256,
        DLD1_GENES: EXPECTED_DLD1_GENES_SHA256,
    }
    for path, expected in expected_hashes.items():
        require(path.is_file(), f"missing bound input: {path}")
        require(sha256_file(path) == expected, f"bound input hash changed: {path}")

    parent_manifest = read_json(PARENT_PLAN_MANIFEST)
    require(
        parent_manifest["selection_count"] == DEVELOPMENT_COUNT,
        "development count drifted",
    )
    require(parent_manifest["test_count"] == TEST_COUNT, "test count drifted")
    require(
        parent_manifest["response_gene_count"] == RESPONSE_COUNT,
        "response count drifted",
    )

    development = read_lines(PARENT_SELECTION_TARGETS)
    test = read_lines(PARENT_TEST_TARGETS)
    all_targets = read_lines(PARENT_ALL_TARGETS)
    responses = read_lines(PARENT_RESPONSE_GENES)
    panel = set(PARENT.read_csv_column(PARENT.PANEL_PATH, "target_gene"))
    dld1_targets = set(read_lines(DLD1_TARGETS))
    require(len(development) == DEVELOPMENT_COUNT, "development axis length changed")
    require(len(test) == TEST_COUNT, "test axis length changed")
    require(len(responses) == RESPONSE_COUNT, "response axis length changed")
    require(set(development).isdisjoint(test), "development/test targets overlap")
    require(
        set(development) | set(test) == set(all_targets),
        "development/test union differs from the frozen gallery target set",
    )
    require(
        set(responses).isdisjoint(panel),
        "response axis contains a VCC panel target coordinate",
    )
    require(
        set(responses).isdisjoint(all_targets),
        "response axis contains a gallery target coordinate",
    )
    development_active = np.asarray(
        [index for index, target in enumerate(development) if target in dld1_targets],
        dtype=np.int16,
    )
    test_active = np.asarray(
        [index for index, target in enumerate(test) if target in dld1_targets],
        dtype=np.int16,
    )
    require(
        len(development_active) == DEVELOPMENT_SUPPORTED,
        f"development DLD1 support changed: {len(development_active)}",
    )
    require(
        len(test_active) == TEST_SUPPORTED,
        f"test DLD1 support changed: {len(test_active)}",
    )
    return {
        "development_targets": development,
        "test_targets": test,
        "response_genes": responses,
        "development_active_indices": development_active,
        "test_active_indices": test_active,
        "development_target_logical_sha256": sha256_lines(development),
        "test_target_logical_sha256": sha256_lines(test),
        "response_gene_logical_sha256": sha256_lines(responses),
    }


def metadata_report() -> dict[str, Any]:
    metadata = validate_static_metadata()
    return {
        "experiment": "138-dld1-replicate-sign",
        "effect_arrays_opened": False,
        "development_targets": len(metadata["development_targets"]),
        "development_dld1_supported": len(metadata["development_active_indices"]),
        "development_exact_fallback": DEVELOPMENT_COUNT
        - len(metadata["development_active_indices"]),
        "test_targets": len(metadata["test_targets"]),
        "test_dld1_supported": len(metadata["test_active_indices"]),
        "test_exact_fallback": TEST_COUNT - len(metadata["test_active_indices"]),
        "response_genes": len(metadata["response_genes"]),
        "parent_plan_manifest_sha256": sha256_file(PARENT_PLAN_MANIFEST),
        "dld1_manifest_sha256": sha256_file(DLD1_MANIFEST),
    }


def contiguous_folds(size: int, folds: int = FOLD_COUNT) -> np.ndarray:
    require(size >= folds, "too few targets for folds")
    labels = np.empty(size, dtype=np.int8)
    for fold, indices in enumerate(np.array_split(np.arange(size), folds)):
        labels[indices] = fold
    return labels


def derangements(
    rng: np.random.Generator, *, draws: int, size: int, dtype: np.dtype[Any]
) -> np.ndarray:
    require(size > 1, "derangement requires at least two entries")
    base = np.arange(size)
    result = np.empty((draws, size), dtype=dtype)
    for draw in range(draws):
        for _ in range(100_000):
            candidate = rng.permutation(base)
            if np.all(candidate != base):
                result[draw] = candidate
                break
        else:
            raise ContractError(f"failed to produce derangement draw {draw}")
    return result


def bootstrap_draws(block_labels: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    count = int(np.max(block_labels)) + 1
    require(count > 1, "bootstrap block count too small")
    return rng.integers(
        0,
        count,
        size=(BOOTSTRAP_DRAWS, count),
        dtype=np.uint16,
    )


def create_plan() -> None:
    require(not PLAN_DIR.exists(), f"{PLAN_DIR} exists; refusing overwrite")
    require(not AUTHORIZATION.exists(), "authorization already exists")
    require(
        not DEVELOPMENT_DIR.exists()
        and not CANDIDATE_FREEZE.exists()
        and not TEST_AUTHORIZATION.exists()
        and not TEST_DIR.exists(),
        "later-stage artifact already exists",
    )
    metadata = validate_static_metadata()
    development = metadata["development_targets"]
    test = metadata["test_targets"]

    development_features = PARENT.feature_rows(development)
    test_features = PARENT.feature_rows(test)
    development_blocks = PARENT.reciprocal_feature_blocks(development_features)
    test_blocks = PARENT.reciprocal_feature_blocks(test_features)
    with np.load(PARENT_RANDOMNESS, allow_pickle=False) as parent_randomness:
        require(
            np.array_equal(test_blocks, parent_randomness["test_block_labels"]),
            "test DepMap block construction differs from Experiment 132",
        )

    development_folds = contiguous_folds(DEVELOPMENT_COUNT)
    require(
        [int(np.sum(development_folds == fold)) for fold in range(FOLD_COUNT)]
        == [123, 123, 123, 122, 122],
        "development fold sizes changed",
    )

    target_rng = np.random.default_rng(TARGET_NULL_SEED)
    response_rng = np.random.default_rng(RESPONSE_NULL_SEED)
    bootstrap_rng = np.random.default_rng(BOOTSTRAP_SEED)
    development_target_nulls = derangements(
        target_rng,
        draws=NULL_DRAWS,
        size=DEVELOPMENT_SUPPORTED,
        dtype=np.int16,
    )
    test_target_nulls = derangements(
        target_rng,
        draws=NULL_DRAWS,
        size=TEST_SUPPORTED,
        dtype=np.int16,
    )
    response_nulls = derangements(
        response_rng,
        draws=NULL_DRAWS,
        size=RESPONSE_COUNT,
        dtype=np.int16,
    )

    PLAN_DIR.mkdir(parents=True)
    atomic_savez(
        PLAN_DIR / "frozen-plan-arrays.npz",
        development_fold_labels=development_folds,
        development_active_indices=metadata["development_active_indices"],
        test_active_indices=metadata["test_active_indices"],
        development_target_nulls=development_target_nulls,
        test_target_nulls=test_target_nulls,
        response_coordinate_nulls=response_nulls,
        development_block_labels=development_blocks,
        test_block_labels=test_blocks,
        development_bootstrap_draws=bootstrap_draws(development_blocks, bootstrap_rng),
        test_bootstrap_draws=bootstrap_draws(test_blocks, bootstrap_rng),
    )
    plan = {
        "schema_version": 1,
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "created_before_experiment_138_effect_array_access": True,
        "parent_selection_targets_file_sha256": sha256_file(PARENT_SELECTION_TARGETS),
        "parent_test_targets_file_sha256": sha256_file(PARENT_TEST_TARGETS),
        "parent_response_genes_file_sha256": sha256_file(PARENT_RESPONSE_GENES),
        "development_count": DEVELOPMENT_COUNT,
        "development_supported_count": DEVELOPMENT_SUPPORTED,
        "development_folds": [
            int(np.sum(development_folds == fold)) for fold in range(FOLD_COUNT)
        ],
        "test_count": TEST_COUNT,
        "test_supported_count": TEST_SUPPORTED,
        "response_gene_count": RESPONSE_COUNT,
        "development_block_count": int(development_blocks.max()) + 1,
        "test_block_count": int(test_blocks.max()) + 1,
        "null_draws_each": NULL_DRAWS,
        "bootstrap_draws_each": BOOTSTRAP_DRAWS,
        "candidate_grid": [arm.name for arm in ARMS],
        "arrays": {
            "frozen-plan-arrays.npz": {
                "bytes": (PLAN_DIR / "frozen-plan-arrays.npz").stat().st_size,
                "sha256": sha256_file(PLAN_DIR / "frozen-plan-arrays.npz"),
            }
        },
    }
    atomic_write_json(PLAN_DIR / "plan-manifest.json", plan)
    print(json.dumps(plan, indent=2, sort_keys=True))


def validate_plan() -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    require(PLAN_DIR.is_dir(), "metadata plan missing")
    plan = read_json(PLAN_DIR / "plan-manifest.json")
    require(
        plan["created_before_experiment_138_effect_array_access"] is True, "unsafe plan"
    )
    require(plan["development_count"] == DEVELOPMENT_COUNT, "development count changed")
    require(plan["test_count"] == TEST_COUNT, "test count changed")
    require(plan["response_gene_count"] == RESPONSE_COUNT, "response count changed")
    require(
        plan["candidate_grid"] == [arm.name for arm in ARMS], "candidate grid changed"
    )
    array_path = PLAN_DIR / "frozen-plan-arrays.npz"
    expected = plan["arrays"]["frozen-plan-arrays.npz"]
    require(
        array_path.stat().st_size == expected["bytes"], "plan array byte count changed"
    )
    require(sha256_file(array_path) == expected["sha256"], "plan array hash changed")
    with np.load(array_path, allow_pickle=False) as data:
        arrays = {key: np.asarray(data[key]) for key in data.files}
    require(
        arrays["development_fold_labels"].shape == (DEVELOPMENT_COUNT,),
        "development fold shape changed",
    )
    require(
        arrays["development_target_nulls"].shape == (NULL_DRAWS, DEVELOPMENT_SUPPORTED),
        "development target-null shape changed",
    )
    require(
        arrays["test_target_nulls"].shape == (NULL_DRAWS, TEST_SUPPORTED),
        "test target-null shape changed",
    )
    require(
        arrays["response_coordinate_nulls"].shape == (NULL_DRAWS, RESPONSE_COUNT),
        "response-null shape changed",
    )
    return plan, arrays


def used_input_paths() -> tuple[Path, ...]:
    paths = [
        PARENT_RUN,
        PARENT_AUTHORIZATION,
        PARENT_PLAN_MANIFEST,
        PARENT_ALL_TARGETS,
        PARENT_SELECTION_TARGETS,
        PARENT_TEST_TARGETS,
        PARENT_RESPONSE_GENES,
        PARENT_RANDOMNESS,
        PARENT.PANEL_PATH,
        PARENT.DEPMAP_PATH,
        PARENT.DEPMAP_RECEIPT,
        DLD1_MANIFEST,
        DLD1_TARGETS,
        DLD1_GENES,
        DLD1_LOW1,
        DLD1_LOW2,
        DLD1_MASK1,
        DLD1_MASK2,
    ]
    for root in PARENT.SOURCE_ROOTS.values():
        paths.extend((root / "perts.txt", root / "genes.txt", root / "delta.npy"))
    paths.extend(
        (
            PARENT.SOURCE_ROOTS["pan_h1"] / "measured_mask.npy",
            PARENT.SOURCE_ROOTS["pan_h1"] / "artifact-manifest.json",
        )
    )
    return tuple(dict.fromkeys(paths))


def implementation_paths() -> tuple[Path, ...]:
    return (README, PROTOCOL, HERE / "run.py", TESTS)


def freeze() -> None:
    require(not AUTHORIZATION.exists(), "authorization already exists")
    require(
        not DEVELOPMENT_LOCK.exists()
        and not DEVELOPMENT_DIR.exists()
        and not CANDIDATE_FREEZE.exists()
        and not TEST_AUTHORIZATION.exists()
        and not TEST_LOCK.exists()
        and not TEST_DIR.exists(),
        "effect-stage artifact already exists",
    )
    validate_static_metadata()
    plan, _ = validate_plan()
    bindings: dict[str, Any] = {}
    for path in used_input_paths():
        require(path.is_file(), f"missing input: {path}")
        relative = str(path.relative_to(ROOT))
        print(f"[hash] {relative}", flush=True)
        bindings[relative] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    implementation = {
        str(path.relative_to(ROOT)): sha256_file(path)
        for path in implementation_paths()
    }
    authorization = {
        "schema_version": 1,
        "experiment": "138-dld1-replicate-sign",
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "created_before_development_or_test_outcome_access": True,
        "development_access_authorized": True,
        "test_access_requires_passing_candidate_freeze": True,
        "cell_eval2_authorized": False,
        "counts_or_submission_authorized": False,
        "commit_or_push_authorized": False,
        "protocol_sha256": sha256_file(PROTOCOL),
        "plan_manifest_sha256": sha256_file(PLAN_DIR / "plan-manifest.json"),
        "parent_authorization_sha256": sha256_file(PARENT_AUTHORIZATION),
        "dld1_manifest_sha256": sha256_file(DLD1_MANIFEST),
        "implementation": implementation,
        "inputs": bindings,
        "plan": {
            "development_count": plan["development_count"],
            "test_count": plan["test_count"],
            "response_gene_count": plan["response_gene_count"],
            "candidate_grid": plan["candidate_grid"],
            "null_draws_each": plan["null_draws_each"],
            "bootstrap_draws_each": plan["bootstrap_draws_each"],
        },
    }
    atomic_write_json(AUTHORIZATION, authorization)
    print(json.dumps(authorization, indent=2, sort_keys=True))


def validate_authorization() -> dict[str, Any]:
    require(AUTHORIZATION.is_file(), "PREVALUE_AUTHORIZATION.v1.json missing")
    value = read_json(AUTHORIZATION)
    require(
        value["development_access_authorized"] is True, "development access disabled"
    )
    require(value["cell_eval2_authorized"] is False, "cell-eval2 unexpectedly enabled")
    require(
        value["counts_or_submission_authorized"] is False, "counts unexpectedly enabled"
    )
    require(value["protocol_sha256"] == sha256_file(PROTOCOL), "protocol changed")
    require(
        value["plan_manifest_sha256"] == sha256_file(PLAN_DIR / "plan-manifest.json"),
        "plan changed",
    )
    require(
        value["parent_authorization_sha256"] == sha256_file(PARENT_AUTHORIZATION),
        "parent authorization changed",
    )
    require(
        value["dld1_manifest_sha256"] == sha256_file(DLD1_MANIFEST),
        "DLD1 manifest changed",
    )
    for relative, expected in value["implementation"].items():
        require(
            sha256_file(ROOT / relative) == expected,
            f"implementation changed: {relative}",
        )
    for relative, expected in value["inputs"].items():
        path = ROOT / relative
        require(path.is_file(), f"input missing: {relative}")
        require(
            path.stat().st_size == expected["bytes"], f"input size changed: {relative}"
        )
        require(
            sha256_file(path) == expected["sha256"], f"input hash changed: {relative}"
        )
    validate_static_metadata()
    validate_plan()
    return value


def load_dld1_proposal(
    targets: Sequence[str], response_genes: Sequence[str]
) -> Proposal:
    dld1_targets = read_lines(DLD1_TARGETS)
    dld1_genes = read_lines(DLD1_GENES)
    target_index = {target: index for index, target in enumerate(dld1_targets)}
    gene_index = {gene: index for index, gene in enumerate(dld1_genes)}
    require(
        all(gene in gene_index for gene in response_genes),
        "Experiment 132 response gene missing from DLD1 canonical axis",
    )
    supported = np.asarray([target in target_index for target in targets], dtype=bool)
    supported_positions = np.flatnonzero(supported)
    source_target_columns = np.asarray(
        [target_index[targets[index]] for index in supported_positions], dtype=np.int64
    )
    source_gene_rows = np.asarray(
        [gene_index[gene] for gene in response_genes], dtype=np.int64
    )
    shape = (len(dld1_genes), len(dld1_targets))
    output_shape = (len(targets), len(response_genes))
    values: dict[str, np.ndarray] = {}
    for name, path, dtype in (
        ("low1", DLD1_LOW1, np.float32),
        ("low2", DLD1_LOW2, np.float32),
        ("measured1", DLD1_MASK1, np.bool_),
        ("measured2", DLD1_MASK2, np.bool_),
    ):
        matrix = np.load(path, mmap_mode="r")
        require(matrix.shape == shape, f"{path.name}: shape changed")
        aligned = np.zeros(output_shape, dtype=dtype)
        aligned[supported_positions] = np.asarray(
            matrix[np.ix_(source_gene_rows, source_target_columns)].T,
            dtype=dtype,
        )
        values[name] = aligned
    require(np.all(np.isfinite(values["low1"])), "DLD1 Low1 contains nonfinite values")
    require(np.all(np.isfinite(values["low2"])), "DLD1 Low2 contains nonfinite values")
    require(
        np.all(values["low1"][~values["measured1"]] == 0.0),
        "DLD1 Low1 unmeasured coordinates are not zero",
    )
    require(
        np.all(values["low2"][~values["measured2"]] == 0.0),
        "DLD1 Low2 unmeasured coordinates are not zero",
    )
    return Proposal(
        low1=values["low1"],
        low2=values["low2"],
        measured1=values["measured1"],
        measured2=values["measured2"],
        supported=supported,
    )


def permute_target_labels(
    proposal: Proposal, active_indices: np.ndarray, permutation: np.ndarray
) -> Proposal:
    active = np.asarray(active_indices, dtype=np.int64)
    perm = np.asarray(permutation, dtype=np.int64)
    require(perm.shape == active.shape, "target permutation shape changed")
    require(
        np.array_equal(np.sort(perm), np.arange(len(active))),
        "invalid target permutation",
    )
    source = active[perm]

    def remap(values: np.ndarray) -> np.ndarray:
        output = np.asarray(values).copy()
        output[active] = values[source]
        return output

    return Proposal(
        low1=remap(proposal.low1),
        low2=remap(proposal.low2),
        measured1=remap(proposal.measured1),
        measured2=remap(proposal.measured2),
        supported=proposal.supported.copy(),
    )


def permute_response_coordinates(
    proposal: Proposal, permutation: np.ndarray
) -> Proposal:
    perm = np.asarray(permutation, dtype=np.int64)
    require(
        perm.shape == (proposal.low1.shape[1],),
        "response permutation shape changed",
    )
    require(
        np.array_equal(np.sort(perm), np.arange(len(perm))),
        "invalid response permutation",
    )
    return Proposal(
        low1=proposal.low1[:, perm],
        low2=proposal.low2[:, perm],
        measured1=proposal.measured1[:, perm],
        measured2=proposal.measured2[:, perm],
        supported=proposal.supported.copy(),
    )


def stable_top_k_mask(
    magnitude: np.ndarray, eligible: np.ndarray, k: int
) -> np.ndarray:
    require(magnitude.shape == eligible.shape, "top-k input shape mismatch")
    result = np.zeros_like(eligible, dtype=bool)
    canonical_index = np.arange(magnitude.shape[1], dtype=np.int64)
    for row in range(magnitude.shape[0]):
        indices = canonical_index[eligible[row]]
        if len(indices) == 0:
            continue
        order = np.lexsort((indices, -magnitude[row, indices]))
        chosen = indices[order[: min(k, len(indices))]]
        result[row, chosen] = True
    return result


def actionable_coordinates(proposal: Proposal, k: int) -> tuple[np.ndarray, np.ndarray]:
    common = (
        proposal.measured1
        & proposal.measured2
        & np.isfinite(proposal.low1)
        & np.isfinite(proposal.low2)
    )
    sign1 = np.sign(proposal.low1)
    sign2 = np.sign(proposal.low2)
    same_nonzero = (sign1 == sign2) & (sign1 != 0.0)
    top1 = stable_top_k_mask(np.abs(proposal.low1), common, k)
    top2 = stable_top_k_mask(np.abs(proposal.low2), common, k)
    actionable = common & same_nonzero & top1 & top2
    actionable &= proposal.supported[:, None]
    consensus = np.where(actionable, sign1, 0.0).astype(np.int8)
    return actionable, consensus


def transform_incumbent(
    incumbent: np.ndarray, proposal: Proposal, arm: Arm
) -> tuple[np.ndarray, dict[str, Any]]:
    base = np.asarray(incumbent, dtype=np.float64)
    require(base.shape == proposal.low1.shape, "incumbent/proposal shape mismatch")
    if arm.k is None:
        return base.copy(), {
            "arm": arm.name,
            "actionable_coordinates": 0,
            "effective_targets": 0,
            "fallback_exact": True,
            "sign_preserved": True,
            "norm_preserved": True,
        }

    actionable, consensus = actionable_coordinates(proposal, arm.k)
    output = base.copy()
    match = np.sign(base) == consensus
    output[actionable & match] *= np.exp(arm.beta)
    output[actionable & ~match] *= np.exp(-arm.beta)
    before = np.linalg.norm(base, axis=1)
    after = np.linalg.norm(output, axis=1)
    scalable = proposal.supported & (before > EPS) & (after > EPS)
    output[scalable] *= (before[scalable] / after[scalable])[:, None]
    output[~proposal.supported] = base[~proposal.supported]

    require(np.all(np.isfinite(output)), f"{arm.name}: nonfinite candidate")
    require(
        np.array_equal(output[~proposal.supported], base[~proposal.supported]),
        f"{arm.name}: unsupported fallback changed",
    )
    require(
        np.array_equal(np.signbit(output), np.signbit(base)),
        f"{arm.name}: sign bit changed",
    )
    require(
        np.array_equal(output == 0.0, base == 0.0),
        f"{arm.name}: zero pattern changed",
    )
    require(
        np.allclose(
            np.linalg.norm(output[proposal.supported], axis=1),
            before[proposal.supported],
            rtol=1e-12,
            atol=1e-12,
        ),
        f"{arm.name}: row norm changed",
    )
    return output, {
        "arm": arm.name,
        "actionable_coordinates": int(np.count_nonzero(actionable)),
        "effective_targets": int(np.count_nonzero(np.any(actionable, axis=1))),
        "supported_targets": int(np.count_nonzero(proposal.supported)),
        "fallback_targets": int(np.count_nonzero(~proposal.supported)),
        "fallback_exact": True,
        "sign_preserved": True,
        "zero_pattern_preserved": True,
        "norm_preserved": True,
    }


def pds_queries(
    prediction: np.ndarray, truth: np.ndarray, query_indices: np.ndarray
) -> np.ndarray:
    require(prediction.shape == truth.shape, "PDS shape mismatch")
    similarity = cosine_similarity(prediction, truth)
    gallery = np.arange(len(truth), dtype=np.int64)
    return pds_from_similarity(similarity, query_indices, gallery)


def cosine_similarity(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    require(prediction.shape == truth.shape, "PDS shape mismatch")
    require(len(prediction) > 1, "complete gallery too small")
    return PARENT.row_unit(prediction) @ PARENT.row_unit(truth).T


def pds_from_similarity(
    similarity: np.ndarray,
    query_indices: np.ndarray,
    gallery_indices: np.ndarray,
) -> np.ndarray:
    require(similarity.ndim == 2, "PDS similarity must be a matrix")
    require(similarity.shape[0] == similarity.shape[1], "PDS gallery must be square")
    queries = np.asarray(query_indices, dtype=np.int64)
    gallery = np.asarray(gallery_indices, dtype=np.int64)
    require(len(gallery) > 1, "PDS gallery too small")
    gallery_position = {int(value): index for index, value in enumerate(gallery)}
    require(
        all(int(query) in gallery_position for query in queries),
        "every PDS query must be present in its gallery",
    )
    selected = similarity[np.ix_(queries, gallery)]
    own_columns = np.asarray(
        [gallery_position[int(query)] for query in queries], dtype=np.int64
    )
    own = selected[np.arange(len(queries)), own_columns]
    greater = np.sum(selected > own[:, None], axis=1)
    equal = np.sum(selected == own[:, None], axis=1)
    midrank = 1.0 + greater + 0.5 * (equal - 1)
    return 1.0 - (midrank - 1.0) / (len(gallery) - 1.0)


def signed_credit_queries(
    prediction: np.ndarray, truth: np.ndarray, query_indices: np.ndarray
) -> np.ndarray:
    queries = np.asarray(query_indices, dtype=np.int64)
    return PARENT.signed_credit_per_target(prediction[queries], truth[queries])


def score_arm(
    arm: Arm,
    identity_by_group: Mapping[str, np.ndarray],
    truth_by_dataset: Mapping[str, np.ndarray],
    proposal: Proposal,
    identity_similarity: Mapping[str, np.ndarray],
) -> tuple[
    dict[str, dict[str, np.ndarray]],
    dict[str, Any],
    dict[str, np.ndarray],
]:
    count = proposal.low1.shape[0]
    queries = np.arange(count, dtype=np.int64)
    score: dict[str, dict[str, np.ndarray]] = {}
    invariants: dict[str, Any] = {}
    candidates: dict[str, np.ndarray] = {}
    candidate_similarity: dict[str, np.ndarray] = {}
    for group in PARENT.HELD_GROUPS:
        candidates[group], invariants[group] = transform_incumbent(
            identity_by_group[group], proposal, arm
        )
    for dataset in PARENT.SCORED_DATASETS:
        group = PARENT.DATASET_GROUP[dataset]
        if arm.k is None:
            candidate_similarity[dataset] = identity_similarity[dataset]
        else:
            candidate_similarity[dataset] = cosine_similarity(
                candidates[group], truth_by_dataset[dataset]
            )
        candidate_pds = pds_from_similarity(
            candidate_similarity[dataset], queries, queries
        )
        identity_pds = pds_from_similarity(
            identity_similarity[dataset], queries, queries
        )
        candidate_credit = signed_credit_queries(
            candidates[group], truth_by_dataset[dataset], queries
        )
        identity_credit = signed_credit_queries(
            identity_by_group[group], truth_by_dataset[dataset], queries
        )
        score[dataset] = {
            "pds": candidate_pds - identity_pds,
            "signed_credit": candidate_credit - identity_credit,
        }
    return score, invariants, candidate_similarity


def score_all_arms(
    identity_by_group: Mapping[str, np.ndarray],
    truth_by_dataset: Mapping[str, np.ndarray],
    proposal: Proposal,
) -> tuple[
    dict[str, dict[str, dict[str, np.ndarray]]],
    dict[str, Any],
    dict[str, dict[str, np.ndarray]],
    dict[str, np.ndarray],
]:
    scores: dict[str, dict[str, dict[str, np.ndarray]]] = {}
    invariants: dict[str, Any] = {}
    candidate_similarity: dict[str, dict[str, np.ndarray]] = {}
    identity_similarity = {
        dataset: cosine_similarity(
            identity_by_group[PARENT.DATASET_GROUP[dataset]],
            truth_by_dataset[dataset],
        )
        for dataset in PARENT.SCORED_DATASETS
    }
    for arm in ARMS:
        (
            scores[arm.name],
            invariants[arm.name],
            candidate_similarity[arm.name],
        ) = score_arm(
            arm,
            identity_by_group,
            truth_by_dataset,
            proposal,
            identity_similarity,
        )
    return scores, invariants, candidate_similarity, identity_similarity


def aggregate_score(
    score: Mapping[str, Mapping[str, np.ndarray]], indices: np.ndarray
) -> dict[str, Any]:
    selected = np.asarray(indices, dtype=np.int64)
    datasets: dict[str, Any] = {}
    for dataset in PARENT.SCORED_DATASETS:
        datasets[dataset] = {
            metric: float(np.mean(score[dataset][metric][selected]))
            for metric in ("pds", "signed_credit")
        }
    macro: dict[str, float] = {}
    per_target: dict[str, np.ndarray] = {}
    for metric in ("pds", "signed_credit"):
        k562 = 0.5 * (
            score["k562_gwps"][metric][selected]
            + score["k562_essential"][metric][selected]
        )
        values = np.mean(
            np.stack(
                [
                    k562,
                    score["rpe1"][metric][selected],
                    score["nadig_jurkat"][metric][selected],
                ]
            ),
            axis=0,
        )
        per_target[metric] = values
        macro[metric] = float(np.mean(values))
    return {"datasets": datasets, "macro": macro, "per_target": per_target}


def eligible_for_selection(summary: Mapping[str, Any]) -> bool:
    return all(
        summary["datasets"][dataset][metric] >= -EPS
        for dataset in PARENT.SCORED_DATASETS
        for metric in ("pds", "signed_credit")
    )


def select_arm(
    scores: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
    indices: np.ndarray,
) -> tuple[Arm, dict[str, Any]]:
    summaries = {arm.name: aggregate_score(scores[arm.name], indices) for arm in ARMS}
    eligible = [arm for arm in ARMS if eligible_for_selection(summaries[arm.name])]
    require(
        any(arm.name == "identity" for arm in eligible), "identity became ineligible"
    )
    best = ARM_BY_NAME["identity"]
    for arm in eligible:
        current = summaries[arm.name]["macro"]
        incumbent = summaries[best.name]["macro"]
        if current["pds"] > incumbent["pds"] + EPS:
            best = arm
        elif abs(current["pds"] - incumbent["pds"]) <= EPS:
            if current["signed_credit"] > incumbent["signed_credit"] + EPS:
                best = arm
            elif (
                abs(current["signed_credit"] - incumbent["signed_credit"]) <= EPS
                and arm.complexity < best.complexity
            ):
                best = arm
    return best, summaries


def nested_selection(
    scores: Mapping[str, Mapping[str, Mapping[str, np.ndarray]]],
    fold_labels: np.ndarray,
    candidate_similarity: Mapping[str, Mapping[str, np.ndarray]],
    identity_similarity: Mapping[str, np.ndarray],
) -> dict[str, Any]:
    count = len(fold_labels)
    oof = {
        dataset: {
            metric: np.empty(count, dtype=np.float64)
            for metric in ("pds", "signed_credit")
        }
        for dataset in PARENT.SCORED_DATASETS
    }
    folds: list[dict[str, Any]] = []
    all_indices = np.arange(count, dtype=np.int64)
    for fold in range(FOLD_COUNT):
        held = all_indices[fold_labels == fold]
        inner = all_indices[fold_labels != fold]
        inner_scores: dict[str, dict[str, dict[str, np.ndarray]]] = {}
        inner_local = np.arange(len(inner), dtype=np.int64)
        for arm in ARMS:
            inner_scores[arm.name] = {}
            for dataset in PARENT.SCORED_DATASETS:
                inner_scores[arm.name][dataset] = {
                    "pds": pds_from_similarity(
                        candidate_similarity[arm.name][dataset],
                        inner,
                        inner,
                    )
                    - pds_from_similarity(
                        identity_similarity[dataset],
                        inner,
                        inner,
                    ),
                    "signed_credit": scores[arm.name][dataset]["signed_credit"][inner],
                }
        selected, inner_summaries = select_arm(inner_scores, inner_local)
        held_summary = aggregate_score(scores[selected.name], held)
        for dataset in PARENT.SCORED_DATASETS:
            for metric in ("pds", "signed_credit"):
                oof[dataset][metric][held] = scores[selected.name][dataset][metric][
                    held
                ]
        folds.append(
            {
                "fold": fold,
                "held_count": len(held),
                "selected_arm": selected.name,
                "inner_macro": inner_summaries[selected.name]["macro"],
                "held": {
                    "datasets": held_summary["datasets"],
                    "macro": held_summary["macro"],
                },
            }
        )
    summary = aggregate_score(oof, all_indices)
    final_arm, all_summaries = select_arm(scores, all_indices)
    return {
        "folds": folds,
        "oof_score": oof,
        "oof_summary": summary,
        "final_arm": final_arm,
        "all_development_summaries": all_summaries,
    }


def positive_share(values: Sequence[float]) -> float:
    positive = np.maximum(np.asarray(values, dtype=np.float64), 0.0)
    total = float(positive.sum())
    if total <= EPS:
        return 1.0
    return float(positive.max() / total)


def concentration_report(
    score: Mapping[str, Mapping[str, np.ndarray]],
    block_labels: np.ndarray,
) -> dict[str, Any]:
    report: dict[str, Any] = {}
    for metric in ("pds", "signed_credit"):
        family_values = [
            0.5
            * (
                float(np.mean(score["k562_gwps"][metric]))
                + float(np.mean(score["k562_essential"][metric]))
            ),
            float(np.mean(score["rpe1"][metric])),
            float(np.mean(score["nadig_jurkat"][metric])),
        ]
        macro_per_target = np.mean(
            np.stack(
                [
                    0.5
                    * (score["k562_gwps"][metric] + score["k562_essential"][metric]),
                    score["rpe1"][metric],
                    score["nadig_jurkat"][metric],
                ]
            ),
            axis=0,
        )
        block_sums = np.bincount(
            np.asarray(block_labels, dtype=np.int64),
            weights=macro_per_target,
            minlength=int(np.max(block_labels)) + 1,
        )
        report[metric] = {
            "family_positive_gain_share_max": positive_share(family_values),
            "target_block_positive_gain_share_max": positive_share(block_sums),
        }
    return report


def bootstrap_report(
    score: Mapping[str, Mapping[str, np.ndarray]],
    block_labels: np.ndarray,
    draws: np.ndarray,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    report: dict[str, Any] = {}
    arrays: dict[str, np.ndarray] = {}
    for metric in ("pds", "signed_credit"):
        macro_per_target = np.mean(
            np.stack(
                [
                    0.5
                    * (score["k562_gwps"][metric] + score["k562_essential"][metric]),
                    score["rpe1"][metric],
                    score["nadig_jurkat"][metric],
                ]
            ),
            axis=0,
        )
        values = PARENT.block_bootstrap(
            macro_per_target,
            np.asarray(block_labels, dtype=np.int64),
            np.asarray(draws, dtype=np.int64),
        )
        arrays[metric] = values
        report[metric] = {
            "draws": len(values),
            "lower_95": float(np.quantile(values, 0.025, method="linear")),
            "upper_95": float(np.quantile(values, 0.975, method="linear")),
        }
    return report, arrays


def null_distribution_development(
    base_proposal: Proposal,
    active_indices: np.ndarray,
    permutations: np.ndarray,
    kind: str,
    identity: Mapping[str, np.ndarray],
    truth: Mapping[str, np.ndarray],
    fold_labels: np.ndarray,
    limit_gib: float,
) -> dict[str, np.ndarray]:
    output = {
        "pds": np.empty(len(permutations), dtype=np.float64),
        "signed_credit": np.empty(len(permutations), dtype=np.float64),
    }
    for draw, permutation in enumerate(permutations):
        proposal = (
            permute_target_labels(base_proposal, active_indices, permutation)
            if kind == "target"
            else permute_response_coordinates(base_proposal, permutation)
        )
        scores, _, candidate_similarity, identity_similarity = score_all_arms(
            identity, truth, proposal
        )
        nested = nested_selection(
            scores, fold_labels, candidate_similarity, identity_similarity
        )
        for metric in output:
            output[metric][draw] = nested["oof_summary"]["macro"][metric]
        if draw % 5 == 4:
            memory_checkpoint(f"{kind} null {draw + 1}/{len(permutations)}", limit_gib)
    return output


def null_distribution_test(
    arm: Arm,
    base_proposal: Proposal,
    active_indices: np.ndarray,
    permutations: np.ndarray,
    kind: str,
    identity: Mapping[str, np.ndarray],
    truth: Mapping[str, np.ndarray],
    limit_gib: float,
) -> dict[str, np.ndarray]:
    output = {
        "pds": np.empty(len(permutations), dtype=np.float64),
        "signed_credit": np.empty(len(permutations), dtype=np.float64),
    }
    indices = np.arange(base_proposal.low1.shape[0], dtype=np.int64)
    for draw, permutation in enumerate(permutations):
        proposal = (
            permute_target_labels(base_proposal, active_indices, permutation)
            if kind == "target"
            else permute_response_coordinates(base_proposal, permutation)
        )
        identity_similarity = {
            dataset: cosine_similarity(
                identity[PARENT.DATASET_GROUP[dataset]], truth[dataset]
            )
            for dataset in PARENT.SCORED_DATASETS
        }
        score, _, _ = score_arm(arm, identity, truth, proposal, identity_similarity)
        summary = aggregate_score(score, indices)
        for metric in output:
            output[metric][draw] = summary["macro"][metric]
        if draw % 10 == 9:
            memory_checkpoint(
                f"test {kind} null {draw + 1}/{len(permutations)}", limit_gib
            )
    return output


def null_report(
    observed: Mapping[str, float],
    target_nulls: Mapping[str, np.ndarray],
    response_nulls: Mapping[str, np.ndarray],
) -> dict[str, Any]:
    report: dict[str, Any] = {}
    for kind, values_by_metric in (
        ("target_label", target_nulls),
        ("response_coordinate", response_nulls),
    ):
        report[kind] = {}
        for metric, values in values_by_metric.items():
            report[kind][metric] = {
                "draws": len(values),
                "p95": float(np.quantile(values, 0.95, method="linear")),
                "add_one_p": float(
                    (1 + np.sum(values >= observed[metric])) / (len(values) + 1)
                ),
            }
    return report


def fold_checks(nested: Mapping[str, Any]) -> dict[str, bool]:
    checks: dict[str, bool] = {}
    for fold in nested["folds"]:
        fold_id = fold["fold"]
        for dataset in PARENT.SCORED_DATASETS:
            for metric in ("pds", "signed_credit"):
                checks[f"fold_{fold_id}_{dataset}_{metric}_nonnegative"] = (
                    fold["held"]["datasets"][dataset][metric] >= -EPS
                )
    return checks


def stage_checks(
    *,
    nested: Mapping[str, Any] | None,
    fixed_summary: Mapping[str, Any],
    bootstrap: Mapping[str, Any],
    nulls: Mapping[str, Any],
    concentration: Mapping[str, Any],
    invariants_pass: bool,
) -> dict[str, bool]:
    checks: dict[str, bool] = {}
    observed_summary = nested["oof_summary"] if nested is not None else fixed_summary
    if nested is not None:
        checks["selected_final_arm_nonidentity"] = (
            nested["final_arm"].name != "identity"
        )
        checks.update(fold_checks(nested))
        checks["nested_macro_pds_at_least_0.003"] = (
            observed_summary["macro"]["pds"] >= MIN_MACRO_GAIN
        )
        checks["nested_macro_signed_credit_at_least_0.003"] = (
            observed_summary["macro"]["signed_credit"] >= MIN_MACRO_GAIN
        )
    for dataset in PARENT.SCORED_DATASETS:
        checks[f"fixed_{dataset}_pds_nonnegative"] = (
            fixed_summary["datasets"][dataset]["pds"] >= -EPS
        )
        checks[f"fixed_{dataset}_signed_credit_nonnegative"] = (
            fixed_summary["datasets"][dataset]["signed_credit"] >= -EPS
        )
    checks["fixed_macro_pds_at_least_0.003"] = (
        fixed_summary["macro"]["pds"] >= MIN_MACRO_GAIN
    )
    checks["fixed_macro_signed_credit_at_least_0.003"] = (
        fixed_summary["macro"]["signed_credit"] >= MIN_MACRO_GAIN
    )
    checks["pds_bootstrap_lower_positive"] = bootstrap["pds"]["lower_95"] > 0.0
    checks["signed_credit_bootstrap_lower_positive"] = (
        bootstrap["signed_credit"]["lower_95"] > 0.0
    )
    for metric in ("pds", "signed_credit"):
        checks[f"{metric}_above_target_null_p95"] = (
            observed_summary["macro"][metric] > nulls["target_label"][metric]["p95"]
        )
        checks[f"{metric}_above_response_null_p95"] = (
            observed_summary["macro"][metric]
            > nulls["response_coordinate"][metric]["p95"]
        )
        checks[f"{metric}_family_concentration_at_most_0.5"] = (
            concentration[metric]["family_positive_gain_share_max"] <= MAX_CONCENTRATION
        )
        checks[f"{metric}_block_concentration_at_most_0.5"] = (
            concentration[metric]["target_block_positive_gain_share_max"]
            <= MAX_CONCENTRATION
        )
    checks["all_candidate_invariants"] = invariants_pass
    return checks


def load_stage(
    targets: Sequence[str], response_genes: Sequence[str], limit_gib: float
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], Proposal]:
    effects = PARENT.load_split_effects(targets, response_genes, limit_gib)
    contexts = PARENT.context_effects(effects)
    identity = PARENT.source_means_by_held(contexts)
    proposal = load_dld1_proposal(targets, response_genes)
    memory_checkpoint("loaded Experiment 132 outcomes and DLD1 proposal", limit_gib)
    return effects, identity, proposal


def run_development(confirm: str | None, limit_gib: float) -> None:
    require(confirm == DEVELOP_CONFIRMATION, "development confirmation string missing")
    validate_authorization()
    require(not DEVELOPMENT_LOCK.exists(), "development access already attempted")
    require(not DEVELOPMENT_DIR.exists(), "development result exists")
    require(not CANDIDATE_FREEZE.exists(), "candidate freeze exists")
    require(
        not TEST_LOCK.exists() and not TEST_DIR.exists(), "test stage already exists"
    )
    metadata = validate_static_metadata()
    _, plan_arrays = validate_plan()
    DEVELOPMENT_LOCK.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_json(
        DEVELOPMENT_LOCK,
        {
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "authorization_sha256": sha256_file(AUTHORIZATION),
            "targets": DEVELOPMENT_COUNT,
            "test_rows_accessed": False,
            "irreversible_development_access": True,
        },
    )
    DEVELOPMENT_DIR.mkdir(parents=True)
    targets = metadata["development_targets"]
    response_genes = metadata["response_genes"]
    truth, identity, proposal = load_stage(targets, response_genes, limit_gib)
    require(
        np.array_equal(
            np.flatnonzero(proposal.supported),
            plan_arrays["development_active_indices"],
        ),
        "development DLD1 support changed after plan",
    )
    (
        scores,
        invariants,
        candidate_similarity,
        identity_similarity,
    ) = score_all_arms(identity, truth, proposal)
    nested = nested_selection(
        scores,
        plan_arrays["development_fold_labels"],
        candidate_similarity,
        identity_similarity,
    )
    final_arm = nested["final_arm"]
    all_indices = np.arange(DEVELOPMENT_COUNT, dtype=np.int64)
    fixed_summary = aggregate_score(scores[final_arm.name], all_indices)
    bootstrap, bootstrap_arrays = bootstrap_report(
        nested["oof_score"],
        plan_arrays["development_block_labels"],
        plan_arrays["development_bootstrap_draws"],
    )
    concentration = concentration_report(
        nested["oof_score"], plan_arrays["development_block_labels"]
    )
    target_null_arrays = null_distribution_development(
        proposal,
        plan_arrays["development_active_indices"],
        plan_arrays["development_target_nulls"],
        "target",
        identity,
        truth,
        plan_arrays["development_fold_labels"],
        limit_gib,
    )
    response_null_arrays = null_distribution_development(
        proposal,
        plan_arrays["development_active_indices"],
        plan_arrays["response_coordinate_nulls"],
        "response",
        identity,
        truth,
        plan_arrays["development_fold_labels"],
        limit_gib,
    )
    nulls = null_report(
        nested["oof_summary"]["macro"], target_null_arrays, response_null_arrays
    )
    invariants_pass = all(
        all(bool(value) for key, value in group.items() if isinstance(value, bool))
        for arm_report in invariants.values()
        for group in arm_report.values()
    )
    checks = stage_checks(
        nested=nested,
        fixed_summary=fixed_summary,
        bootstrap=bootstrap,
        nulls=nulls,
        concentration=concentration,
        invariants_pass=invariants_pass,
    )
    result = {
        "schema_version": 1,
        "stage": "development",
        "final_arm": final_arm.name,
        "nested_folds": nested["folds"],
        "nested_oof": {
            "datasets": nested["oof_summary"]["datasets"],
            "macro": nested["oof_summary"]["macro"],
        },
        "final_arm_all_development": {
            "datasets": fixed_summary["datasets"],
            "macro": fixed_summary["macro"],
        },
        "bootstrap": bootstrap,
        "nulls": nulls,
        "concentration": concentration,
        "candidate_invariants": invariants,
        "checks": checks,
        "passes": all(checks.values()),
        "test_rows_accessed": False,
        "production_authorized": False,
        "peak_rss_gib": peak_rss_gib(),
    }
    atomic_write_json(DEVELOPMENT_DIR / "result.json", result)
    oof_arrays = {
        f"oof__{dataset}__{metric}": nested["oof_score"][dataset][metric]
        for dataset in PARENT.SCORED_DATASETS
        for metric in ("pds", "signed_credit")
    }
    atomic_savez(
        DEVELOPMENT_DIR / "statistics.npz",
        **oof_arrays,
        target_null__pds=target_null_arrays["pds"],
        target_null__signed_credit=target_null_arrays["signed_credit"],
        response_null__pds=response_null_arrays["pds"],
        response_null__signed_credit=response_null_arrays["signed_credit"],
        bootstrap__pds=bootstrap_arrays["pds"],
        bootstrap__signed_credit=bootstrap_arrays["signed_credit"],
    )
    atomic_write_json(
        DEVELOPMENT_DIR / "artifact-manifest.json",
        output_manifest(DEVELOPMENT_DIR, exclude={"artifact-manifest.json"}),
    )
    freeze_value = {
        "schema_version": 1,
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "authorization_sha256": sha256_file(AUTHORIZATION),
        "development_access_lock_sha256": sha256_file(DEVELOPMENT_LOCK),
        "development_result_sha256": sha256_file(DEVELOPMENT_DIR / "result.json"),
        "development_manifest_sha256": sha256_file(
            DEVELOPMENT_DIR / "artifact-manifest.json"
        ),
        "selected_final_arm": final_arm.name,
        "development_passes": result["passes"],
        "untouched_test_access_authorized": result["passes"],
        "test_target_count": TEST_COUNT,
        "test_rows_accessed": False,
        "production_authorized": False,
    }
    atomic_write_json(CANDIDATE_FREEZE, freeze_value)
    print(json.dumps(result, indent=2, sort_keys=True))


def validate_candidate_freeze() -> dict[str, Any]:
    require(CANDIDATE_FREEZE.is_file(), "development candidate freeze missing")
    value = read_json(CANDIDATE_FREEZE)
    require(
        value["authorization_sha256"] == sha256_file(AUTHORIZATION),
        "authorization changed after development",
    )
    require(
        value["development_access_lock_sha256"] == sha256_file(DEVELOPMENT_LOCK),
        "development access lock changed",
    )
    require(
        value["development_result_sha256"]
        == sha256_file(DEVELOPMENT_DIR / "result.json"),
        "development result changed",
    )
    require(
        value["development_manifest_sha256"]
        == sha256_file(DEVELOPMENT_DIR / "artifact-manifest.json"),
        "development manifest changed",
    )
    require(value["selected_final_arm"] in ARM_BY_NAME, "unknown frozen arm")
    result = read_json(DEVELOPMENT_DIR / "result.json")
    validate_candidate_freeze_semantics(value, result)
    return value


def validate_candidate_freeze_semantics(
    freeze_value: Mapping[str, Any], development_result: Mapping[str, Any]
) -> None:
    require(
        freeze_value["selected_final_arm"] == development_result["final_arm"],
        "candidate freeze arm differs from development result",
    )
    require(
        freeze_value["development_passes"] is development_result["passes"],
        "candidate freeze pass flag differs from development result",
    )
    require(
        freeze_value["untouched_test_access_authorized"]
        is development_result["passes"],
        "candidate freeze test-access flag differs from development result",
    )


def authorize_test(confirm: str | None) -> None:
    require(
        confirm == TEST_AUTHORIZATION_CONFIRMATION,
        "post-development test-authorization confirmation string missing",
    )
    validate_authorization()
    candidate_freeze = validate_candidate_freeze()
    require(candidate_freeze["development_passes"] is True, "development did not pass")
    require(not TEST_AUTHORIZATION.exists(), "test authorization already exists")
    require(
        not TEST_LOCK.exists() and not TEST_DIR.exists(), "test stage already exists"
    )
    value = {
        "schema_version": 1,
        "experiment": "138-dld1-replicate-sign",
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "created_after_independent_development_review": True,
        "untouched_test_access_authorized": True,
        "authorization_sha256": sha256_file(AUTHORIZATION),
        "candidate_freeze_sha256": sha256_file(CANDIDATE_FREEZE),
        "development_result_sha256": sha256_file(DEVELOPMENT_DIR / "result.json"),
        "selected_final_arm": candidate_freeze["selected_final_arm"],
        "production_authorized": False,
    }
    atomic_write_json(TEST_AUTHORIZATION, value)
    print(json.dumps(value, indent=2, sort_keys=True))


def validate_test_authorization() -> dict[str, Any]:
    require(TEST_AUTHORIZATION.is_file(), "post-development test authorization missing")
    value = read_json(TEST_AUTHORIZATION)
    require(
        value["created_after_independent_development_review"] is True,
        "test authorization is not review-gated",
    )
    require(
        value["untouched_test_access_authorized"] is True,
        "untouched-test access disabled",
    )
    require(
        value["authorization_sha256"] == sha256_file(AUTHORIZATION),
        "prevalue authorization changed before test",
    )
    require(
        value["candidate_freeze_sha256"] == sha256_file(CANDIDATE_FREEZE),
        "candidate freeze changed after test authorization",
    )
    require(
        value["development_result_sha256"]
        == sha256_file(DEVELOPMENT_DIR / "result.json"),
        "development result changed after test authorization",
    )
    candidate_freeze = validate_candidate_freeze()
    require(
        value["selected_final_arm"] == candidate_freeze["selected_final_arm"],
        "test authorization arm differs from candidate freeze",
    )
    return value


def run_test(confirm: str | None, limit_gib: float) -> None:
    require(confirm == TEST_CONFIRMATION, "untouched-test confirmation string missing")
    validate_authorization()
    candidate_freeze = validate_candidate_freeze()
    validate_test_authorization()
    require(candidate_freeze["development_passes"] is True, "development did not pass")
    require(
        candidate_freeze["untouched_test_access_authorized"] is True,
        "untouched-test access denied",
    )
    require(
        candidate_freeze["selected_final_arm"] != "identity", "identity was selected"
    )
    require(not TEST_LOCK.exists(), "untouched-test access already attempted")
    require(not TEST_DIR.exists(), "untouched-test result exists")
    metadata = validate_static_metadata()
    _, plan_arrays = validate_plan()
    TEST_LOCK.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_json(
        TEST_LOCK,
        {
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "authorization_sha256": sha256_file(AUTHORIZATION),
            "candidate_freeze_sha256": sha256_file(CANDIDATE_FREEZE),
            "targets": TEST_COUNT,
            "irreversible_untouched_test_access": True,
        },
    )
    TEST_DIR.mkdir(parents=True)
    targets = metadata["test_targets"]
    response_genes = metadata["response_genes"]
    truth, identity, proposal = load_stage(targets, response_genes, limit_gib)
    require(
        np.array_equal(
            np.flatnonzero(proposal.supported), plan_arrays["test_active_indices"]
        ),
        "test DLD1 support changed after plan",
    )
    arm = ARM_BY_NAME[candidate_freeze["selected_final_arm"]]
    identity_similarity = {
        dataset: cosine_similarity(
            identity[PARENT.DATASET_GROUP[dataset]], truth[dataset]
        )
        for dataset in PARENT.SCORED_DATASETS
    }
    score, invariants, _ = score_arm(
        arm, identity, truth, proposal, identity_similarity
    )
    indices = np.arange(TEST_COUNT, dtype=np.int64)
    fixed_summary = aggregate_score(score, indices)
    bootstrap, bootstrap_arrays = bootstrap_report(
        score,
        plan_arrays["test_block_labels"],
        plan_arrays["test_bootstrap_draws"],
    )
    concentration = concentration_report(score, plan_arrays["test_block_labels"])
    target_null_arrays = null_distribution_test(
        arm,
        proposal,
        plan_arrays["test_active_indices"],
        plan_arrays["test_target_nulls"],
        "target",
        identity,
        truth,
        limit_gib,
    )
    response_null_arrays = null_distribution_test(
        arm,
        proposal,
        plan_arrays["test_active_indices"],
        plan_arrays["response_coordinate_nulls"],
        "response",
        identity,
        truth,
        limit_gib,
    )
    nulls = null_report(
        fixed_summary["macro"], target_null_arrays, response_null_arrays
    )
    invariants_pass = all(
        bool(value)
        for group in invariants.values()
        for key, value in group.items()
        if isinstance(value, bool)
    )
    checks = stage_checks(
        nested=None,
        fixed_summary=fixed_summary,
        bootstrap=bootstrap,
        nulls=nulls,
        concentration=concentration,
        invariants_pass=invariants_pass,
    )
    result = {
        "schema_version": 1,
        "stage": "untouched_source_test",
        "candidate": arm.name,
        "observed": {
            "datasets": fixed_summary["datasets"],
            "macro": fixed_summary["macro"],
        },
        "bootstrap": bootstrap,
        "nulls": nulls,
        "concentration": concentration,
        "candidate_invariants": invariants,
        "checks": checks,
        "passes": all(checks.values()),
        "decision": "PASS_FOR_SEPARATE_REVIEW_ONLY"
        if all(checks.values())
        else "NO-GO",
        "production_authorized": False,
        "peak_rss_gib": peak_rss_gib(),
    }
    atomic_write_json(TEST_DIR / "result.json", result)
    per_target = {
        f"{dataset}__{metric}_delta": score[dataset][metric]
        for dataset in PARENT.SCORED_DATASETS
        for metric in ("pds", "signed_credit")
    }
    atomic_savez(
        TEST_DIR / "statistics.npz",
        **per_target,
        target_null__pds=target_null_arrays["pds"],
        target_null__signed_credit=target_null_arrays["signed_credit"],
        response_null__pds=response_null_arrays["pds"],
        response_null__signed_credit=response_null_arrays["signed_credit"],
        bootstrap__pds=bootstrap_arrays["pds"],
        bootstrap__signed_credit=bootstrap_arrays["signed_credit"],
    )
    atomic_write_json(
        TEST_DIR / "artifact-manifest.json",
        output_manifest(TEST_DIR, exclude={"artifact-manifest.json"}),
    )
    print(json.dumps(result, indent=2, sort_keys=True))


def audit() -> None:
    payload: dict[str, Any] = {
        "metadata": metadata_report(),
        "plan_exists": PLAN_DIR.exists(),
        "authorization_exists": AUTHORIZATION.exists(),
        "development_accessed": DEVELOPMENT_LOCK.exists(),
        "development_result_exists": DEVELOPMENT_DIR.exists(),
        "candidate_freeze_exists": CANDIDATE_FREEZE.exists(),
        "test_authorization_exists": TEST_AUTHORIZATION.exists(),
        "untouched_test_accessed": TEST_LOCK.exists(),
        "untouched_test_result_exists": TEST_DIR.exists(),
    }
    if PLAN_DIR.exists():
        payload["plan"] = validate_plan()[0]
    if AUTHORIZATION.exists():
        payload["authorization_sha256"] = sha256_file(AUTHORIZATION)
    if CANDIDATE_FREEZE.exists():
        payload["candidate_freeze"] = validate_candidate_freeze()
    print(json.dumps(payload, indent=2, sort_keys=True))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "command",
        choices=(
            "metadata",
            "plan",
            "freeze",
            "develop",
            "authorize-test",
            "test",
            "audit",
        ),
    )
    parser.add_argument("--confirm")
    parser.add_argument("--max-rss-gib", type=float, default=24.0)
    args = parser.parse_args()
    if args.command == "metadata":
        print(json.dumps(metadata_report(), indent=2, sort_keys=True))
    elif args.command == "plan":
        create_plan()
    elif args.command == "freeze":
        freeze()
    elif args.command == "develop":
        run_development(args.confirm, args.max_rss_gib)
    elif args.command == "authorize-test":
        authorize_test(args.confirm)
    elif args.command == "test":
        run_test(args.confirm, args.max_rss_gib)
    else:
        audit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

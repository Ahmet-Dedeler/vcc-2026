#!/usr/bin/env python3
"""Experiment 139: route existing source vectors using Jurkat landmarks."""

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

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
README = HERE / "README.md"
PROTOCOL = HERE / "preregistered_protocol.json"
TESTS = HERE / "test_run.py"
PLAN_DIR = HERE / "work" / "preregistered-plan-v1"
AUTHORIZATION = HERE / "PREVALUE_AUTHORIZATION.v1.json"
DEV_LOCK = HERE / "work" / "development-access.v1.lock"
DEV_DIR = HERE / "work" / "development-result-v1"
CANDIDATE_FREEZE = HERE / "DEVELOPMENT_CANDIDATE_FREEZE.v1.json"
TEST_AUTHORIZATION = HERE / "TEST_ACCESS_AUTHORIZATION.v1.json"
TEST_LOCK = HERE / "work" / "untouched-test-access.v1.lock"
TEST_DIR = HERE / "work" / "untouched-test-result-v1"

EXP132 = ROOT / "experiments/132-response-decomposition-ridge"
PARENT_RUN = EXP132 / "run.py"
PARENT_AUTH = EXP132 / "SOURCE_ONLY_AUTHORIZATION.v2.json"
PARENT_PLAN = EXP132 / "work/preregistered-freeze-v2"
ALL_TARGETS = PARENT_PLAN / "all_targets.txt"
DEV_TARGETS = PARENT_PLAN / "selection_targets.txt"
TEST_TARGETS = PARENT_PLAN / "test_targets.txt"
RESPONSE_GENES = PARENT_PLAN / "response_genes.txt"
PARENT_RANDOMNESS = PARENT_PLAN / "frozen-randomness.npz"
PANEL_PATH = ROOT / "data/challenge/pert_counts.csv"
SOURCE_DIR = ROOT / "experiments/107-jurkat-high-moi/work/frozen-source-model"
SOURCE_MODEL = SOURCE_DIR / "source-model.npz"
SOURCE_MANIFEST = SOURCE_DIR / "manifest.json"

EXPECTED = {
    PARENT_RUN: "e78e079fc3b4a715e8b412b0a48bfe44a2f3f6cfcadb8a65a789351c53877778",
    PARENT_AUTH: "9bca6e620359e44b5440d0b8588412f40f274304ae67d9f55fd7e7d03341e197",
    PARENT_PLAN
    / "plan-manifest.json": "170d55ed5f812eb85f5503af4bb9d7dec8d17dce6e58fe966d85bf0eec5185d1",
    ALL_TARGETS: "f5a9d65cc1d01765dc5c61649fb876c9e2ddf4f316791b28439324129fea4d51",
    DEV_TARGETS: "74a070a538455ffc2a528dae3094d21e047905463779b2a5206d681f5ca6ebb5",
    TEST_TARGETS: "bd37bcacb4eea90da11902bbb8f4736068a7ff67c001739e3407de2b4bf73fc6",
    RESPONSE_GENES: "24743754eead5b8ddd81634bd12feac1df431eb2b9581197755b75c4c90949c8",
    PANEL_PATH: "f57edd7b912ebd718efc7ee9d0f334772513e7cc418d133ce525470e373b3276",
    SOURCE_MODEL: "03759a7b2cafa86d2643a58669235aca7fc38a33599aa403c67439f8979a7e26",
    SOURCE_MANIFEST: "ae7cf5bf374e05d18080e48c7318a6405d96b77edd2e5e1bfe513ca81bff4a68",
}
DEV_COUNT, TEST_COUNT, RESPONSE_COUNT = 613, 410, 4992
DEV_SUPPORTED, TEST_SUPPORTED, LANDMARK_COUNT = 609, 409, 138
FOLDS, NULL_DRAWS, BOOTSTRAP_DRAWS = 5, 49, 20_000
TARGET_SEED, LANDMARK_SEED, BOOTSTRAP_SEED = 139_202_610, 139_202_611, 139_202_612
MIN_NONZERO_LANDMARKS = 16
MIN_GAIN, MAX_CONCENTRATION, EPS = 0.003, 0.5, 1e-12
GIB = 1024**3
DEV_CONFIRM = "I_CONFIRM_EXPERIMENT_139_DEVELOPMENT_ONLY"
TEST_AUTH_CONFIRM = "I_REVIEWED_EXPERIMENT_139_DEVELOPMENT_AND_AUTHORIZE_TEST"
TEST_CONFIRM = "I_CONFIRM_EXPERIMENT_139_UNTOUCHED_TEST"


class ContractError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def import_parent() -> Any:
    spec = importlib.util.spec_from_file_location("experiment132_for_139", PARENT_RUN)
    require(spec is not None and spec.loader is not None, "cannot load Experiment 132")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


PARENT = import_parent()


@dataclass(frozen=True)
class Arm:
    name: str
    rho: float
    temperature: float

    @property
    def identity(self) -> bool:
        return self.rho == 0.0

    @property
    def complexity(self) -> tuple[float, float]:
        return (self.rho, -self.temperature)


ARMS = (
    Arm("identity", 0.0, 1.0),
    Arm("rho025_temp025", 0.25, 0.25),
    Arm("rho025_temp010", 0.25, 0.10),
    Arm("rho050_temp025", 0.50, 0.25),
    Arm("rho050_temp010", 0.50, 0.10),
)
ARM_BY_NAME = {arm.name: arm for arm in ARMS}


@dataclass(frozen=True)
class LandmarkEvidence:
    values: np.ndarray
    supported: np.ndarray


def sha256_file(path: Path, block: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while piece := handle.read(block):
            digest.update(piece)
    return digest.hexdigest()


def read_lines(path: Path) -> tuple[str, ...]:
    values = tuple(x.strip() for x in path.read_text().splitlines() if x.strip())
    require(values and len(values) == len(set(values)), f"invalid axis: {path}")
    return values


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    try:
        with os.fdopen(fd, "w") as handle:
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
    fd, temporary = tempfile.mkstemp(dir=path.parent, suffix=".npz")
    os.close(fd)
    try:
        np.savez_compressed(temporary, **arrays)
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def peak_rss() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if os.uname().sysname == "Darwin" else value / 1024**2


def memory_check(label: str, limit: float) -> None:
    value = peak_rss()
    print(f"[rss] {label}: {value:.3f} GiB", flush=True)
    if limit > 0 and value >= limit:
        raise MemoryError(f"{label}: {value:.3f} GiB >= {limit:.3f} GiB")


def source_metadata() -> dict[str, Any]:
    for path, expected in EXPECTED.items():
        require(path.is_file(), f"missing input: {path}")
        require(sha256_file(path) == expected, f"input hash changed: {path}")
    development, test, responses = (
        read_lines(DEV_TARGETS),
        read_lines(TEST_TARGETS),
        read_lines(RESPONSE_GENES),
    )
    all_targets = read_lines(ALL_TARGETS)
    panel = set(PARENT.read_csv_column(PANEL_PATH, "target_gene"))
    with np.load(SOURCE_MODEL, allow_pickle=False) as archive:
        targets = np.asarray(archive["target_symbols"]).astype(str)
        transferable = np.asarray(archive["transferable"], dtype=bool)
        source_responses = np.asarray(archive["response_symbols"]).astype(str)
        exportable = np.asarray(archive["response_exportable"], dtype=bool)
    target_index = {name: i for i, name in enumerate(targets)}
    response_index = {name: i for i, name in enumerate(source_responses)}
    response_set = set(responses)
    landmarks = tuple(
        name
        for name, usable in zip(source_responses, exportable, strict=True)
        if usable and name in response_set
    )
    landmark_parent_indices = np.asarray(
        [responses.index(x) for x in landmarks], dtype=np.int16
    )
    landmark_source_indices = np.asarray(
        [response_index[x] for x in landmarks], dtype=np.int16
    )

    def support(axis: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
        source_indices = np.asarray(
            [target_index.get(x, -1) for x in axis], dtype=np.int32
        )
        supported = np.asarray(
            [i >= 0 and transferable[i] for i in source_indices], dtype=bool
        )
        return source_indices, supported

    dev_source, dev_support = support(development)
    test_source, test_support = support(test)
    require(
        len(development) == DEV_COUNT and len(test) == TEST_COUNT, "target count drift"
    )
    require(len(responses) == RESPONSE_COUNT, "response count drift")
    require(set(development).isdisjoint(test), "development/test targets overlap")
    require(
        set(development) | set(test) == set(all_targets),
        "development/test union differs from frozen gallery targets",
    )
    require(
        set(responses).isdisjoint(panel),
        "response axis contains a VCC panel target coordinate",
    )
    require(
        set(responses).isdisjoint(all_targets),
        "response axis contains a gallery target coordinate",
    )
    require(len(landmarks) == LANDMARK_COUNT, "landmark count drift")
    require(int(dev_support.sum()) == DEV_SUPPORTED, "development support drift")
    require(int(test_support.sum()) == TEST_SUPPORTED, "test support drift")
    return {
        "development": development,
        "test": test,
        "responses": responses,
        "landmarks": landmarks,
        "landmark_parent_indices": landmark_parent_indices,
        "landmark_source_indices": landmark_source_indices,
        "development_source_indices": dev_source,
        "development_supported": dev_support,
        "test_source_indices": test_source,
        "test_supported": test_support,
    }


def metadata_report() -> dict[str, Any]:
    data = source_metadata()
    return {
        "effect_values_opened": False,
        "development_targets": len(data["development"]),
        "development_transferable": int(data["development_supported"].sum()),
        "development_fallback": int((~data["development_supported"]).sum()),
        "test_targets": len(data["test"]),
        "test_transferable": int(data["test_supported"].sum()),
        "test_fallback": int((~data["test_supported"]).sum()),
        "response_genes": len(data["responses"]),
        "shared_exportable_landmarks": len(data["landmarks"]),
    }


def contiguous_folds(size: int) -> np.ndarray:
    labels = np.empty(size, dtype=np.int8)
    for fold, indices in enumerate(np.array_split(np.arange(size), FOLDS)):
        labels[indices] = fold
    return labels


def derangements(rng: np.random.Generator, draws: int, size: int) -> np.ndarray:
    base = np.arange(size)
    output = np.empty((draws, size), dtype=np.int16)
    for draw in range(draws):
        for _ in range(100_000):
            candidate = rng.permutation(base)
            if np.all(candidate != base):
                output[draw] = candidate
                break
        else:
            raise ContractError("could not create derangement")
    return output


def make_bootstrap(labels: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    count = int(labels.max()) + 1
    return rng.integers(0, count, size=(BOOTSTRAP_DRAWS, count), dtype=np.uint16)


def create_plan() -> None:
    require(not PLAN_DIR.exists(), "plan exists")
    require(not AUTHORIZATION.exists(), "authorization exists")
    data = source_metadata()
    dev_features = PARENT.feature_rows(data["development"])
    test_features = PARENT.feature_rows(data["test"])
    dev_blocks = PARENT.reciprocal_feature_blocks(dev_features)
    test_blocks = PARENT.reciprocal_feature_blocks(test_features)
    with np.load(PARENT_RANDOMNESS, allow_pickle=False) as parent:
        require(
            np.array_equal(test_blocks, parent["test_block_labels"]),
            "test blocks differ from Experiment 132",
        )
    folds = contiguous_folds(DEV_COUNT)
    require(
        [int(np.sum(folds == x)) for x in range(FOLDS)] == [123, 123, 123, 122, 122],
        "fold sizes changed",
    )
    target_rng = np.random.default_rng(TARGET_SEED)
    landmark_rng = np.random.default_rng(LANDMARK_SEED)
    bootstrap_rng = np.random.default_rng(BOOTSTRAP_SEED)
    PLAN_DIR.mkdir(parents=True)
    atomic_npz(
        PLAN_DIR / "plan-arrays.npz",
        development_folds=folds,
        development_blocks=dev_blocks,
        test_blocks=test_blocks,
        development_bootstrap=make_bootstrap(dev_blocks, bootstrap_rng),
        test_bootstrap=make_bootstrap(test_blocks, bootstrap_rng),
        development_source_indices=data["development_source_indices"],
        development_supported=data["development_supported"],
        test_source_indices=data["test_source_indices"],
        test_supported=data["test_supported"],
        landmark_parent_indices=data["landmark_parent_indices"],
        landmark_source_indices=data["landmark_source_indices"],
        development_target_nulls=derangements(target_rng, NULL_DRAWS, DEV_SUPPORTED),
        test_target_nulls=derangements(target_rng, NULL_DRAWS, TEST_SUPPORTED),
        landmark_nulls=derangements(landmark_rng, NULL_DRAWS, LANDMARK_COUNT),
    )
    arrays = PLAN_DIR / "plan-arrays.npz"
    manifest = {
        "schema_version": 1,
        "created_before_source_effect_or_outcome_access": True,
        "development_targets": DEV_COUNT,
        "development_transferable": DEV_SUPPORTED,
        "test_targets": TEST_COUNT,
        "test_transferable": TEST_SUPPORTED,
        "response_genes": RESPONSE_COUNT,
        "landmarks": LANDMARK_COUNT,
        "fold_sizes": [int(np.sum(folds == x)) for x in range(FOLDS)],
        "development_blocks": int(dev_blocks.max()) + 1,
        "test_blocks": int(test_blocks.max()) + 1,
        "candidate_grid": [x.name for x in ARMS],
        "null_draws_each": NULL_DRAWS,
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "arrays": {"bytes": arrays.stat().st_size, "sha256": sha256_file(arrays)},
    }
    atomic_json(PLAN_DIR / "plan-manifest.json", manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))


def validate_plan() -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    manifest = read_json(PLAN_DIR / "plan-manifest.json")
    require(manifest["candidate_grid"] == [x.name for x in ARMS], "grid changed")
    path = PLAN_DIR / "plan-arrays.npz"
    require(path.stat().st_size == manifest["arrays"]["bytes"], "plan size changed")
    require(sha256_file(path) == manifest["arrays"]["sha256"], "plan hash changed")
    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: np.asarray(archive[name]) for name in archive.files}
    require(arrays["development_folds"].shape == (DEV_COUNT,), "fold shape changed")
    require(
        arrays["development_target_nulls"].shape == (NULL_DRAWS, DEV_SUPPORTED),
        "development target-null shape changed",
    )
    require(
        arrays["test_target_nulls"].shape == (NULL_DRAWS, TEST_SUPPORTED),
        "test target-null shape changed",
    )
    require(
        arrays["landmark_nulls"].shape == (NULL_DRAWS, LANDMARK_COUNT),
        "landmark-null shape changed",
    )
    return manifest, arrays


def used_inputs() -> tuple[Path, ...]:
    paths = [
        PARENT_RUN,
        PARENT_AUTH,
        PARENT_PLAN / "plan-manifest.json",
        ALL_TARGETS,
        DEV_TARGETS,
        TEST_TARGETS,
        RESPONSE_GENES,
        PARENT_RANDOMNESS,
        PANEL_PATH,
        PARENT.DEPMAP_PATH,
        PARENT.DEPMAP_RECEIPT,
        SOURCE_MODEL,
        SOURCE_MANIFEST,
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


def freeze() -> None:
    require(not AUTHORIZATION.exists(), "authorization exists")
    require(
        not DEV_LOCK.exists()
        and not DEV_DIR.exists()
        and not CANDIDATE_FREEZE.exists()
        and not TEST_AUTHORIZATION.exists()
        and not TEST_LOCK.exists()
        and not TEST_DIR.exists(),
        "outcome-stage artifact exists",
    )
    source_metadata()
    plan, _ = validate_plan()
    inputs = {}
    for path in used_inputs():
        print(f"[hash] {path.relative_to(ROOT)}", flush=True)
        inputs[str(path.relative_to(ROOT))] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    implementation = {
        str(path.relative_to(ROOT)): sha256_file(path)
        for path in (README, PROTOCOL, HERE / "run.py", TESTS)
    }
    atomic_json(
        AUTHORIZATION,
        {
            "schema_version": 1,
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "created_before_development_or_test_outcome_access": True,
            "development_authorized": True,
            "test_requires_passing_development_freeze": True,
            "protocol_sha256": sha256_file(PROTOCOL),
            "plan_sha256": sha256_file(PLAN_DIR / "plan-manifest.json"),
            "implementation": implementation,
            "inputs": inputs,
            "plan": plan,
            "cell_eval2_counts_submission_authorized": False,
        },
    )


def validate_authorization() -> None:
    require(AUTHORIZATION.is_file(), "authorization missing")
    auth = read_json(AUTHORIZATION)
    require(auth["development_authorized"] is True, "development denied")
    require(auth["protocol_sha256"] == sha256_file(PROTOCOL), "protocol changed")
    require(
        auth["plan_sha256"] == sha256_file(PLAN_DIR / "plan-manifest.json"),
        "plan changed",
    )
    for relative, expected in auth["implementation"].items():
        require(sha256_file(ROOT / relative) == expected, f"code changed: {relative}")
    for relative, expected in auth["inputs"].items():
        path = ROOT / relative
        require(path.stat().st_size == expected["bytes"], f"size changed: {relative}")
        require(sha256_file(path) == expected["sha256"], f"hash changed: {relative}")
    source_metadata()
    validate_plan()


def load_landmarks(
    targets: Sequence[str],
    source_indices: np.ndarray,
    supported: np.ndarray,
    plan: Mapping[str, np.ndarray],
) -> LandmarkEvidence:
    with np.load(SOURCE_MODEL, allow_pickle=False) as archive:
        effects = np.asarray(archive["exported_effects"])
        transferable = np.asarray(archive["transferable"], dtype=bool)
        output = np.zeros((len(targets), LANDMARK_COUNT), dtype=np.float64)
        rows = np.flatnonzero(supported)
        indices = np.asarray(source_indices[rows], dtype=np.int64)
        require(np.all(transferable[indices]), "planned target lost transferability")
        output[rows] = effects[
            np.ix_(indices, np.asarray(plan["landmark_source_indices"], dtype=np.int64))
        ]
    require(np.all(np.isfinite(output)), "nonfinite landmark evidence")
    return LandmarkEvidence(output, np.asarray(supported, dtype=bool))


def target_null(
    evidence: LandmarkEvidence, permutation: np.ndarray
) -> LandmarkEvidence:
    active = np.flatnonzero(evidence.supported)
    perm = np.asarray(permutation, dtype=np.int64)
    require(perm.shape == active.shape, "target null shape changed")
    output = evidence.values.copy()
    output[active] = evidence.values[active[perm]]
    return LandmarkEvidence(output, evidence.supported.copy())


def landmark_null(
    evidence: LandmarkEvidence, permutation: np.ndarray
) -> LandmarkEvidence:
    perm = np.asarray(permutation, dtype=np.int64)
    require(perm.shape == (evidence.values.shape[1],), "landmark null shape changed")
    return LandmarkEvidence(evidence.values[:, perm], evidence.supported.copy())


def row_cosine(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    numerator = np.sum(left * right, axis=-1)
    denominator = np.linalg.norm(left, axis=-1) * np.linalg.norm(right, axis=-1)
    return np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator, dtype=np.float64),
        where=denominator > EPS,
    )


def routing_weights(similarity: np.ndarray, arm: Arm) -> np.ndarray:
    require(similarity.ndim == 2 and similarity.shape[1] == 5, "similarity shape")
    if arm.identity:
        return np.full(similarity.shape, 0.2, dtype=np.float64)
    logits = similarity / arm.temperature
    logits -= logits.max(axis=1, keepdims=True)
    softmax = np.exp(logits)
    softmax /= softmax.sum(axis=1, keepdims=True)
    weights = (1.0 - arm.rho) / 5.0 + arm.rho * softmax
    require(np.all(np.isfinite(weights)), "nonfinite weights")
    require(np.all(weights >= 0.0), "negative weights")
    require(
        np.allclose(weights.sum(axis=1), 1.0, rtol=0.0, atol=1e-14),
        "weights not simplex",
    )
    return weights


def routed_candidate(
    source_stack: np.ndarray,
    landmark_indices: np.ndarray,
    evidence: LandmarkEvidence,
    arm: Arm,
) -> tuple[np.ndarray, dict[str, Any]]:
    sources = np.asarray(source_stack, dtype=np.float64)
    require(sources.shape[0] == 5, "exactly five sources required")
    identity = np.mean(sources, axis=0)
    if arm.identity:
        return identity.copy(), {
            "fallback_exact": True,
            "simplex": True,
            "convex_hull": True,
        }
    source_landmarks = sources[:, :, np.asarray(landmark_indices, dtype=np.int64)]
    source_landmarks = np.moveaxis(source_landmarks, 0, 1)
    similarity = row_cosine(source_landmarks, evidence.values[:, None, :])
    nonzero = np.count_nonzero(evidence.values, axis=1)
    source_norms = np.linalg.norm(source_landmarks, axis=2)
    usable = (
        evidence.supported
        & (nonzero >= MIN_NONZERO_LANDMARKS)
        & np.all(source_norms > EPS, axis=1)
        & (np.linalg.norm(evidence.values, axis=1) > EPS)
    )
    weights = routing_weights(similarity, arm)
    candidate = np.einsum("ns,sng->ng", weights, sources, optimize=True)
    candidate[~usable] = identity[~usable]
    require(np.array_equal(candidate[~usable], identity[~usable]), "fallback changed")
    require(np.all(np.isfinite(candidate)), "candidate nonfinite")
    lower, upper = sources.min(axis=0), sources.max(axis=0)
    require(np.all(candidate >= lower - 1e-12), "candidate below convex hull")
    require(np.all(candidate <= upper + 1e-12), "candidate above convex hull")
    return candidate, {
        "fallback_exact": True,
        "simplex": True,
        "convex_hull": True,
        "supported_targets": int(evidence.supported.sum()),
        "usable_targets": int(usable.sum()),
        "fallback_targets": int((~usable).sum()),
    }


def build_sources(effects: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    contexts = PARENT.context_effects(effects)
    return {
        held: np.stack(
            [contexts[name] for name in PARENT.CONTEXTS if name != held], axis=0
        )
        for held in PARENT.HELD_GROUPS
    }


def cosine_similarity(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    require(prediction.shape == truth.shape, "PDS shape mismatch")
    require(len(prediction) > 1, "PDS gallery too small")
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
    positions = {int(value): index for index, value in enumerate(gallery)}
    require(
        all(int(query) in positions for query in queries),
        "every PDS query must be present in its gallery",
    )
    selected = similarity[np.ix_(queries, gallery)]
    own_columns = np.asarray(
        [positions[int(query)] for query in queries], dtype=np.int64
    )
    own = selected[np.arange(len(queries)), own_columns]
    greater = np.sum(selected > own[:, None], axis=1)
    equal = np.sum(selected == own[:, None], axis=1)
    midrank = 1.0 + greater + 0.5 * (equal - 1)
    return 1.0 - (midrank - 1.0) / (len(gallery) - 1.0)


def score_arm(
    arm: Arm,
    source_stacks: Mapping[str, np.ndarray],
    truth: Mapping[str, np.ndarray],
    evidence: LandmarkEvidence,
    landmark_indices: np.ndarray,
    identity_similarity: Mapping[str, np.ndarray],
) -> tuple[
    dict[str, dict[str, np.ndarray]],
    dict[str, Any],
    dict[str, np.ndarray],
]:
    candidates, invariants = {}, {}
    for group in PARENT.HELD_GROUPS:
        candidates[group], invariants[group] = routed_candidate(
            source_stacks[group], landmark_indices, evidence, arm
        )
    output, candidate_similarity = {}, {}
    indices = np.arange(len(evidence.values), dtype=np.int64)
    for dataset in PARENT.SCORED_DATASETS:
        group = PARENT.DATASET_GROUP[dataset]
        identity = np.mean(source_stacks[group], axis=0)
        if arm.identity:
            candidate_similarity[dataset] = identity_similarity[dataset]
        else:
            candidate_similarity[dataset] = cosine_similarity(
                candidates[group], truth[dataset]
            )
        output[dataset] = {
            "pds": pds_from_similarity(candidate_similarity[dataset], indices, indices)
            - pds_from_similarity(identity_similarity[dataset], indices, indices),
            "signed_credit": PARENT.signed_credit_per_target(
                candidates[group], truth[dataset]
            )
            - PARENT.signed_credit_per_target(identity, truth[dataset]),
        }
    return output, invariants, candidate_similarity


def all_scores(
    source_stacks: Mapping[str, np.ndarray],
    truth: Mapping[str, np.ndarray],
    evidence: LandmarkEvidence,
    landmark_indices: np.ndarray,
) -> tuple[
    dict[str, Any],
    dict[str, Any],
    dict[str, dict[str, np.ndarray]],
    dict[str, np.ndarray],
]:
    scores, invariants, candidate_similarity = {}, {}, {}
    identity_similarity = {
        dataset: cosine_similarity(
            np.mean(source_stacks[PARENT.DATASET_GROUP[dataset]], axis=0),
            truth[dataset],
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
            source_stacks,
            truth,
            evidence,
            landmark_indices,
            identity_similarity,
        )
    return scores, invariants, candidate_similarity, identity_similarity


def aggregate(score: Mapping[str, Any], indices: np.ndarray) -> dict[str, Any]:
    idx = np.asarray(indices, dtype=np.int64)
    datasets = {
        dataset: {
            metric: float(np.mean(score[dataset][metric][idx]))
            for metric in ("pds", "signed_credit")
        }
        for dataset in PARENT.SCORED_DATASETS
    }
    per_target, macro = {}, {}
    for metric in ("pds", "signed_credit"):
        values = np.mean(
            np.stack(
                [
                    0.5
                    * (
                        score["k562_gwps"][metric][idx]
                        + score["k562_essential"][metric][idx]
                    ),
                    score["rpe1"][metric][idx],
                    score["nadig_jurkat"][metric][idx],
                ]
            ),
            axis=0,
        )
        per_target[metric], macro[metric] = values, float(values.mean())
    return {"datasets": datasets, "per_target": per_target, "macro": macro}


def select_arm(
    scores: Mapping[str, Any], indices: np.ndarray
) -> tuple[Arm, dict[str, Any]]:
    summaries = {arm.name: aggregate(scores[arm.name], indices) for arm in ARMS}
    eligible = [
        arm
        for arm in ARMS
        if all(
            summaries[arm.name]["datasets"][dataset][metric] >= -EPS
            for dataset in PARENT.SCORED_DATASETS
            for metric in ("pds", "signed_credit")
        )
    ]
    best = ARM_BY_NAME["identity"]
    for arm in eligible:
        current, old = summaries[arm.name]["macro"], summaries[best.name]["macro"]
        key = (current["pds"], current["signed_credit"])
        old_key = (old["pds"], old["signed_credit"])
        if key[0] > old_key[0] + EPS or (
            abs(key[0] - old_key[0]) <= EPS
            and (
                key[1] > old_key[1] + EPS
                or (
                    abs(key[1] - old_key[1]) <= EPS and arm.complexity < best.complexity
                )
            )
        ):
            best = arm
    return best, summaries


def nested(
    scores: Mapping[str, Any],
    folds: np.ndarray,
    candidate_similarity: Mapping[str, Mapping[str, np.ndarray]],
    identity_similarity: Mapping[str, np.ndarray],
) -> dict[str, Any]:
    count = len(folds)
    oof = {
        dataset: {metric: np.empty(count) for metric in ("pds", "signed_credit")}
        for dataset in PARENT.SCORED_DATASETS
    }
    records = []
    all_indices = np.arange(count)
    for fold in range(FOLDS):
        held, inner = all_indices[folds == fold], all_indices[folds != fold]
        inner_scores = {}
        inner_local = np.arange(len(inner), dtype=np.int64)
        for arm in ARMS:
            inner_scores[arm.name] = {}
            for dataset in PARENT.SCORED_DATASETS:
                inner_scores[arm.name][dataset] = {
                    "pds": pds_from_similarity(
                        candidate_similarity[arm.name][dataset], inner, inner
                    )
                    - pds_from_similarity(identity_similarity[dataset], inner, inner),
                    "signed_credit": scores[arm.name][dataset]["signed_credit"][inner],
                }
        arm, summaries = select_arm(inner_scores, inner_local)
        held_summary = aggregate(scores[arm.name], held)
        for dataset in PARENT.SCORED_DATASETS:
            for metric in ("pds", "signed_credit"):
                oof[dataset][metric][held] = scores[arm.name][dataset][metric][held]
        records.append(
            {
                "fold": fold,
                "selected_arm": arm.name,
                "inner_macro": summaries[arm.name]["macro"],
                "held_datasets": held_summary["datasets"],
                "held_macro": held_summary["macro"],
            }
        )
    final, summaries = select_arm(scores, all_indices)
    return {
        "folds": records,
        "oof": oof,
        "oof_summary": aggregate(oof, all_indices),
        "final_arm": final,
        "final_summary": summaries[final.name],
    }


def positive_share(values: Iterable[float]) -> float:
    values = np.maximum(np.asarray(tuple(values), dtype=np.float64), 0.0)
    return 1.0 if values.sum() <= EPS else float(values.max() / values.sum())


def concentration(score: Mapping[str, Any], blocks: np.ndarray) -> dict[str, Any]:
    result = {}
    for metric in ("pds", "signed_credit"):
        families = [
            0.5
            * (
                score["k562_gwps"][metric].mean()
                + score["k562_essential"][metric].mean()
            ),
            score["rpe1"][metric].mean(),
            score["nadig_jurkat"][metric].mean(),
        ]
        macro = np.mean(
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
        block_values = np.bincount(blocks.astype(int), weights=macro)
        result[metric] = {
            "family": positive_share(families),
            "target_block": positive_share(block_values),
        }
    return result


def bootstrap(
    score: Mapping[str, Any], blocks: np.ndarray, draws: np.ndarray
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    report, arrays = {}, {}
    for metric in ("pds", "signed_credit"):
        macro = np.mean(
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
        values = PARENT.block_bootstrap(macro, blocks.astype(int), draws.astype(int))
        arrays[metric] = values
        report[metric] = {
            "lower_95": float(np.quantile(values, 0.025, method="linear")),
            "upper_95": float(np.quantile(values, 0.975, method="linear")),
        }
    return report, arrays


def run_nulls(
    kind: str,
    permutations: np.ndarray,
    evidence: LandmarkEvidence,
    source_stacks: Mapping[str, np.ndarray],
    truth: Mapping[str, np.ndarray],
    landmarks: np.ndarray,
    folds: np.ndarray | None,
    fixed_arm: Arm | None,
    limit: float,
) -> dict[str, np.ndarray]:
    result = {
        metric: np.empty(len(permutations)) for metric in ("pds", "signed_credit")
    }
    for draw, permutation in enumerate(permutations):
        changed = (
            target_null(evidence, permutation)
            if kind == "target"
            else landmark_null(evidence, permutation)
        )
        scores, _, candidate_similarity, identity_similarity = all_scores(
            source_stacks, truth, changed, landmarks
        )
        summary = (
            nested(scores, folds, candidate_similarity, identity_similarity)[
                "oof_summary"
            ]
            if fixed_arm is None
            else aggregate(scores[fixed_arm.name], np.arange(len(evidence.values)))
        )
        for metric in result:
            result[metric][draw] = summary["macro"][metric]
        if draw % 7 == 6:
            memory_check(f"{kind} null {draw + 1}", limit)
    return result


def null_report(
    observed: Mapping[str, float],
    target_values: Mapping[str, np.ndarray],
    landmark_values: Mapping[str, np.ndarray],
) -> dict[str, Any]:
    return {
        name: {
            metric: {
                "p95": float(np.quantile(values[metric], 0.95, method="linear")),
                "add_one_p": float(
                    (1 + np.sum(values[metric] >= observed[metric]))
                    / (len(values[metric]) + 1)
                ),
            }
            for metric in ("pds", "signed_credit")
        }
        for name, values in (
            ("target_label", target_values),
            ("landmark_coordinate", landmark_values),
        )
    }


def checks(
    nested_result: Mapping[str, Any] | None,
    fixed: Mapping[str, Any],
    boot: Mapping[str, Any],
    nulls: Mapping[str, Any],
    conc: Mapping[str, Any],
    invariant: bool,
) -> dict[str, bool]:
    observed = nested_result["oof_summary"] if nested_result else fixed
    output = {}
    if nested_result:
        output["nonidentity"] = not nested_result["final_arm"].identity
        for fold in nested_result["folds"]:
            for dataset in PARENT.SCORED_DATASETS:
                for metric in ("pds", "signed_credit"):
                    output[f"fold{fold['fold']}_{dataset}_{metric}"] = (
                        fold["held_datasets"][dataset][metric] >= -EPS
                    )
    for dataset in PARENT.SCORED_DATASETS:
        for metric in ("pds", "signed_credit"):
            output[f"{dataset}_{metric}"] = fixed["datasets"][dataset][metric] >= -EPS
    for metric in ("pds", "signed_credit"):
        output[f"macro_{metric}"] = observed["macro"][metric] >= MIN_GAIN
        output[f"fixed_macro_{metric}"] = fixed["macro"][metric] >= MIN_GAIN
        output[f"bootstrap_{metric}"] = boot[metric]["lower_95"] > 0.0
        output[f"target_null_{metric}"] = (
            observed["macro"][metric] > nulls["target_label"][metric]["p95"]
        )
        output[f"landmark_null_{metric}"] = (
            observed["macro"][metric] > nulls["landmark_coordinate"][metric]["p95"]
        )
        output[f"family_concentration_{metric}"] = (
            conc[metric]["family"] <= MAX_CONCENTRATION
        )
        output[f"block_concentration_{metric}"] = (
            conc[metric]["target_block"] <= MAX_CONCENTRATION
        )
    output["invariants"] = invariant
    return output


def load_stage(
    targets: Sequence[str], limit: float
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    effects = PARENT.load_split_effects(targets, read_lines(RESPONSE_GENES), limit)
    return effects, build_sources(effects)


def manifest(directory: Path) -> dict[str, Any]:
    return {
        "files": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in sorted(directory.iterdir())
            if path.is_file() and path.name != "artifact-manifest.json"
        }
    }


def develop(confirm: str | None, limit: float) -> None:
    require(confirm == DEV_CONFIRM, "development confirmation missing")
    validate_authorization()
    require(
        not DEV_LOCK.exists()
        and not DEV_DIR.exists()
        and not CANDIDATE_FREEZE.exists(),
        "development already attempted",
    )
    data, (_, plan) = source_metadata(), validate_plan()
    atomic_json(
        DEV_LOCK,
        {"created_utc": dt.datetime.now(dt.UTC).isoformat(), "test_accessed": False},
    )
    DEV_DIR.mkdir(parents=True)
    truth, source_stacks = load_stage(data["development"], limit)
    evidence = load_landmarks(
        data["development"],
        plan["development_source_indices"],
        plan["development_supported"],
        plan,
    )
    scores, invariants, candidate_similarity, identity_similarity = all_scores(
        source_stacks, truth, evidence, plan["landmark_parent_indices"]
    )
    selected = nested(
        scores,
        plan["development_folds"],
        candidate_similarity,
        identity_similarity,
    )
    fixed = selected["final_summary"]
    boot, boot_arrays = bootstrap(
        selected["oof"], plan["development_blocks"], plan["development_bootstrap"]
    )
    conc = concentration(selected["oof"], plan["development_blocks"])
    target_values = run_nulls(
        "target",
        plan["development_target_nulls"],
        evidence,
        source_stacks,
        truth,
        plan["landmark_parent_indices"],
        plan["development_folds"],
        None,
        limit,
    )
    landmark_values = run_nulls(
        "landmark",
        plan["landmark_nulls"],
        evidence,
        source_stacks,
        truth,
        plan["landmark_parent_indices"],
        plan["development_folds"],
        None,
        limit,
    )
    nulls = null_report(
        selected["oof_summary"]["macro"], target_values, landmark_values
    )
    invariant = all(
        v
        for arm in invariants.values()
        for group in arm.values()
        for v in group.values()
        if isinstance(v, bool)
    )
    gate = checks(selected, fixed, boot, nulls, conc, invariant)
    result = {
        "stage": "development",
        "final_arm": selected["final_arm"].name,
        "folds": selected["folds"],
        "nested_oof": {
            "datasets": selected["oof_summary"]["datasets"],
            "macro": selected["oof_summary"]["macro"],
        },
        "final_arm_all_development": {
            "datasets": fixed["datasets"],
            "macro": fixed["macro"],
        },
        "bootstrap": boot,
        "nulls": nulls,
        "concentration": conc,
        "invariants": invariants,
        "checks": gate,
        "passes": all(gate.values()),
        "test_accessed": False,
    }
    atomic_json(DEV_DIR / "result.json", result)
    atomic_npz(
        DEV_DIR / "statistics.npz",
        target_null_pds=target_values["pds"],
        target_null_credit=target_values["signed_credit"],
        landmark_null_pds=landmark_values["pds"],
        landmark_null_credit=landmark_values["signed_credit"],
        bootstrap_pds=boot_arrays["pds"],
        bootstrap_credit=boot_arrays["signed_credit"],
    )
    atomic_json(DEV_DIR / "artifact-manifest.json", manifest(DEV_DIR))
    atomic_json(
        CANDIDATE_FREEZE,
        {
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "authorization_sha256": sha256_file(AUTHORIZATION),
            "development_result_sha256": sha256_file(DEV_DIR / "result.json"),
            "development_manifest_sha256": sha256_file(
                DEV_DIR / "artifact-manifest.json"
            ),
            "final_arm": selected["final_arm"].name,
            "development_passes": result["passes"],
            "test_authorized": result["passes"],
        },
    )
    print(json.dumps(result, indent=2, sort_keys=True))


def validate_candidate() -> dict[str, Any]:
    value = read_json(CANDIDATE_FREEZE)
    require(
        value["authorization_sha256"] == sha256_file(AUTHORIZATION),
        "authorization changed",
    )
    require(
        value["development_result_sha256"] == sha256_file(DEV_DIR / "result.json"),
        "development changed",
    )
    require(
        value["development_manifest_sha256"]
        == sha256_file(DEV_DIR / "artifact-manifest.json"),
        "manifest changed",
    )
    result = read_json(DEV_DIR / "result.json")
    validate_candidate_semantics(value, result)
    return value


def validate_candidate_semantics(
    freeze_value: Mapping[str, Any], development_result: Mapping[str, Any]
) -> None:
    require(
        freeze_value["final_arm"] == development_result["final_arm"],
        "candidate freeze arm differs from development result",
    )
    require(
        freeze_value["development_passes"] is development_result["passes"],
        "candidate freeze pass flag differs from development result",
    )
    require(
        freeze_value["test_authorized"] is development_result["passes"],
        "candidate freeze test flag differs from development result",
    )


def authorize_test(confirm: str | None) -> None:
    require(confirm == TEST_AUTH_CONFIRM, "test-authorization confirmation missing")
    validate_authorization()
    frozen = validate_candidate()
    require(frozen["development_passes"] is True, "development did not pass")
    require(not TEST_AUTHORIZATION.exists(), "test authorization exists")
    require(not TEST_LOCK.exists() and not TEST_DIR.exists(), "test already attempted")
    value = {
        "schema_version": 1,
        "experiment": "139-jurkat-donor-simplex",
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "created_after_independent_development_review": True,
        "test_authorized": True,
        "authorization_sha256": sha256_file(AUTHORIZATION),
        "candidate_freeze_sha256": sha256_file(CANDIDATE_FREEZE),
        "development_result_sha256": sha256_file(DEV_DIR / "result.json"),
        "final_arm": frozen["final_arm"],
        "production_authorized": False,
    }
    atomic_json(TEST_AUTHORIZATION, value)
    print(json.dumps(value, indent=2, sort_keys=True))


def validate_test_authorization() -> dict[str, Any]:
    require(TEST_AUTHORIZATION.is_file(), "post-development test authorization missing")
    value = read_json(TEST_AUTHORIZATION)
    require(
        value["created_after_independent_development_review"] is True,
        "test authorization is not review-gated",
    )
    require(value["test_authorized"] is True, "test authorization disabled")
    require(
        value["authorization_sha256"] == sha256_file(AUTHORIZATION),
        "prevalue authorization changed before test",
    )
    require(
        value["candidate_freeze_sha256"] == sha256_file(CANDIDATE_FREEZE),
        "candidate freeze changed after test authorization",
    )
    require(
        value["development_result_sha256"] == sha256_file(DEV_DIR / "result.json"),
        "development result changed after test authorization",
    )
    frozen = validate_candidate()
    require(value["final_arm"] == frozen["final_arm"], "authorized arm changed")
    return value


def test_stage(confirm: str | None, limit: float) -> None:
    require(confirm == TEST_CONFIRM, "test confirmation missing")
    validate_authorization()
    frozen = validate_candidate()
    validate_test_authorization()
    require(frozen["development_passes"] and frozen["test_authorized"], "test denied")
    require(frozen["final_arm"] != "identity", "identity selected")
    require(not TEST_LOCK.exists() and not TEST_DIR.exists(), "test already attempted")
    data, (_, plan) = source_metadata(), validate_plan()
    atomic_json(
        TEST_LOCK,
        {
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "candidate_freeze_sha256": sha256_file(CANDIDATE_FREEZE),
        },
    )
    TEST_DIR.mkdir(parents=True)
    truth, source_stacks = load_stage(data["test"], limit)
    evidence = load_landmarks(
        data["test"], plan["test_source_indices"], plan["test_supported"], plan
    )
    arm = ARM_BY_NAME[frozen["final_arm"]]
    identity_similarity = {
        dataset: cosine_similarity(
            np.mean(source_stacks[PARENT.DATASET_GROUP[dataset]], axis=0),
            truth[dataset],
        )
        for dataset in PARENT.SCORED_DATASETS
    }
    score, invariants, _ = score_arm(
        arm,
        source_stacks,
        truth,
        evidence,
        plan["landmark_parent_indices"],
        identity_similarity,
    )
    fixed = aggregate(score, np.arange(TEST_COUNT))
    boot, boot_arrays = bootstrap(score, plan["test_blocks"], plan["test_bootstrap"])
    conc = concentration(score, plan["test_blocks"])
    target_values = run_nulls(
        "target",
        plan["test_target_nulls"],
        evidence,
        source_stacks,
        truth,
        plan["landmark_parent_indices"],
        None,
        arm,
        limit,
    )
    landmark_values = run_nulls(
        "landmark",
        plan["landmark_nulls"],
        evidence,
        source_stacks,
        truth,
        plan["landmark_parent_indices"],
        None,
        arm,
        limit,
    )
    nulls = null_report(fixed["macro"], target_values, landmark_values)
    invariant = all(
        v
        for group in invariants.values()
        for v in group.values()
        if isinstance(v, bool)
    )
    gate = checks(None, fixed, boot, nulls, conc, invariant)
    result = {
        "stage": "untouched_test",
        "arm": arm.name,
        "observed": {"datasets": fixed["datasets"], "macro": fixed["macro"]},
        "bootstrap": boot,
        "nulls": nulls,
        "concentration": conc,
        "invariants": invariants,
        "checks": gate,
        "passes": all(gate.values()),
        "decision": "PASS_FOR_REVIEW_ONLY" if all(gate.values()) else "NO-GO",
    }
    atomic_json(TEST_DIR / "result.json", result)
    atomic_npz(
        TEST_DIR / "statistics.npz",
        target_null_pds=target_values["pds"],
        target_null_credit=target_values["signed_credit"],
        landmark_null_pds=landmark_values["pds"],
        landmark_null_credit=landmark_values["signed_credit"],
        bootstrap_pds=boot_arrays["pds"],
        bootstrap_credit=boot_arrays["signed_credit"],
    )
    atomic_json(TEST_DIR / "artifact-manifest.json", manifest(TEST_DIR))
    print(json.dumps(result, indent=2, sort_keys=True))


def audit() -> None:
    payload = {
        "metadata": metadata_report(),
        "plan_exists": PLAN_DIR.exists(),
        "authorization_exists": AUTHORIZATION.exists(),
        "development_accessed": DEV_LOCK.exists(),
        "candidate_freeze_exists": CANDIDATE_FREEZE.exists(),
        "test_authorization_exists": TEST_AUTHORIZATION.exists(),
        "test_accessed": TEST_LOCK.exists(),
    }
    if PLAN_DIR.exists():
        payload["plan"] = validate_plan()[0]
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
        develop(args.confirm, args.max_rss_gib)
    elif args.command == "authorize-test":
        authorize_test(args.confirm)
    elif args.command == "test":
        test_stage(args.confirm, args.max_rss_gib)
    else:
        audit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

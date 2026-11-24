#!/usr/bin/env python3
"""Strict source-only broad-Pan neighborhood transfer for Experiment 129.

Stages are intentionally separated:

* ``plan`` reads labels, metadata, and STRING embeddings, but no effect arrays.
* ``freeze`` hashes every input and seals the implementation and deterministic plan.
* ``select`` may read only the 100 source-selection target rows and source-only
  predictor rows used to fit the incumbent geometry.
* ``evaluate`` creates a permanent access lock before reading the 55 untouched
  source-test target rows.

No command imports cell-eval2, emits counts, builds a challenge artifact, or
touches held HCT116/HEK293T/KOLF/Pan/DLD1/HCASMC/challenge outcomes.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
import resource
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")

import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
INVENTORY_PATH = HERE / "SOURCE_INVENTORY.json"
AUTHORIZATION_V1_PATH = HERE / "SOURCE_ONLY_AUTHORIZATION.json"
AUTHORIZATION_PATH = HERE / "SOURCE_ONLY_AUTHORIZATION_V2.json"
RECOVERY_001_PATH = HERE / "IMPLEMENTATION_RECOVERY_001.json"
RECOVERY_002_PATH = HERE / "IMPLEMENTATION_RECOVERY_002.json"
PLAN_DIR = HERE / "work" / "preregistered-freeze-v1"
SELECTION_DIR = HERE / "work" / "source-candidate-freeze-v1"
RESULT_DIR = HERE / "work" / "source-only-evaluation-v1"
TEST_ACCESS_LOCK = HERE / "work" / "source-test-access.lock"

PANEL_PATH = ROOT / "data/challenge/pert_counts.csv"
GENES_PATH = ROOT / "data/challenge/gene_names.csv"
STRING_PATH = ROOT / "data/external/slim/gene_string_embeddings.v0.3.h5"
PAN_ROOT = ROOT / "experiments/115-gse295214-h1/work/compact-donor"
EXP107_PLAN = ROOT / "experiments/107-jurkat-high-moi/work/strict-validation-plan.json"
EXP123_TARGETS = (
    ROOT
    / "experiments/123-slim-string-reproduction/work/preregistered-freeze-v1/all_targets.txt"
)

BASE_SOURCES = ("k562_gwps", "orion_hct116", "orion_hek293t")
OUTCOMES = ("rpe1", "nadig_jurkat")
PM_ROOT = ROOT / "data/external/perturbmatch"
TARGET_NAMESPACE = "experiment-129-broad-pan-knn-v1"
FOLD_NAMESPACE = "experiment-129-fold-v1"
EXPECTED_TARGETS = 155
SELECTION_COUNT = 100
TEST_COUNT = 55
FOLDS = 5
K_GRID = (5, 15, 30)
ANGLE_RADIANS = math.radians(5.0)
WHITEN_GAMMA = 0.65
PM_BLEND = 0.50
BASIS_CAP = 12_000
BASIS_NAMESPACE = "experiment-129-whitener-basis-v1"
NULL_DRAWS = 199
NULL_SEED = 12_920_261
BOOTSTRAP_DRAWS = 20_000
BOOTSTRAP_SEED = 12_920_262
CLUSTER_K = 5
CLUSTER_THRESHOLD = 0.80
EXPECTED_SCORE_GENES = 5730
EXPECTED_SCORE_SHA = "2efa164f9910158a591367fc237cf0dce0f2d05d03b5ce38ab76779ea21fd674"
EPS = 1e-12
GIB = 1024**3


class ContractError(RuntimeError):
    """Raised when a frozen source-only boundary is violated."""


@dataclass(frozen=True)
class Candidate:
    kind: str
    k: int

    @property
    def key(self) -> str:
        return "identity" if self.kind == "identity" else f"pan-knn:k={self.k}:angle=5deg"


@dataclass(frozen=True)
class Bank:
    name: str
    root: Path
    perts: tuple[str, ...]
    genes: tuple[str, ...]

    @classmethod
    def load_axis(cls, name: str, root: Path) -> "Bank":
        return cls(
            name,
            root,
            read_lines(root / "perts.txt"),
            read_lines(root / "genes.txt", require_unique=False),
        )

    @property
    def delta_path(self) -> Path:
        return self.root / "delta.npy"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if os.uname().sysname == "Darwin" else value / 1024**2


def memory_checkpoint(stage: str, limit_gib: float) -> None:
    peak = peak_rss_gib()
    print(f"[rss] {stage}: {peak:.3f} GiB", flush=True)
    if limit_gib > 0 and peak >= limit_gib:
        raise MemoryError(f"{stage}: peak RSS {peak:.3f} GiB reached {limit_gib:.3f} GiB")


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode("utf-8")).hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)


def read_lines(path: Path, *, require_unique: bool = True) -> tuple[str, ...]:
    values = tuple(line.strip() for line in path.read_text().splitlines() if line.strip())
    require(bool(values), f"{path}: empty axis")
    if require_unique:
        require(len(values) == len(set(values)), f"{path}: duplicate axis labels")
    return values


def write_lines(path: Path, values: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(values) + "\n", encoding="utf-8")


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames is not None and column in reader.fieldnames, f"{path}: missing {column}")
        values = tuple(
            str(row[column]).strip()
            for row in reader
            if row.get(column) and str(row[column]).strip()
        )
    require(bool(values), f"{path}: empty {column}")
    require(len(values) == len(set(values)), f"{path}: duplicate {column}")
    return values


def candidates() -> tuple[Candidate, ...]:
    return (Candidate("identity", 0), *(Candidate("knn", k) for k in K_GRID))


def stable_order(values: Iterable[str], namespace: str) -> tuple[str, ...]:
    prefix = namespace.encode("utf-8") + b"\0"
    return tuple(
        sorted(
            values,
            key=lambda value: (
                hashlib.sha256(prefix + value.encode("utf-8")).digest(),
                value.encode("utf-8"),
            ),
        )
    )


def deterministic_folds(targets: Sequence[str]) -> np.ndarray:
    prefix = FOLD_NAMESPACE.encode("utf-8") + b"\0"
    folds = np.asarray(
        [
            int.from_bytes(
                hashlib.sha256(prefix + target.encode("utf-8")).digest()[:8], "little"
            )
            % FOLDS
            for target in targets
        ],
        dtype=np.int8,
    )
    require(set(map(int, folds)) == set(range(FOLDS)), "not every selection fold is populated")
    return folds


def derangement(rng: np.random.Generator, size: int) -> np.ndarray:
    base = np.arange(size, dtype=np.int16)
    for _ in range(10_000):
        trial = rng.permutation(base)
        if np.all(trial != base):
            return trial.astype(np.int16)
    raise ContractError(f"failed to build derangement of size {size}")


def row_unit(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    norms = np.linalg.norm(array, axis=1, keepdims=True)
    return np.divide(array, norms, out=np.zeros_like(array), where=norms > EPS)


def reciprocal_knn_clusters(embeddings: np.ndarray) -> np.ndarray:
    unit = row_unit(embeddings)
    similarity = unit @ unit.T
    np.fill_diagonal(similarity, -np.inf)
    neighbors = np.argsort(-similarity, axis=1, kind="stable")[:, :CLUSTER_K]
    parent = np.arange(len(unit), dtype=np.int64)

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = int(parent[value])
        return value

    def union(left: int, right: int) -> None:
        a, b = find(left), find(right)
        if a != b:
            parent[max(a, b)] = min(a, b)

    neighbor_sets = [set(map(int, row)) for row in neighbors]
    for left in range(len(unit)):
        for right_value in neighbors[left]:
            right = int(right_value)
            if left in neighbor_sets[right] and similarity[left, right] >= CLUSTER_THRESHOLD:
                union(left, right)
    roots = [find(index) for index in range(len(unit))]
    mapping = {root: label for label, root in enumerate(sorted(set(roots)))}
    return np.asarray([mapping[root] for root in roots], dtype=np.int16)


def load_string_rows(targets: Sequence[str]) -> np.ndarray:
    values = np.empty((len(targets), 64), dtype=np.float64)
    with h5py.File(STRING_PATH, "r") as handle:
        for index, target in enumerate(targets):
            require(target in handle, f"STRING missing target {target}")
            row = np.asarray(handle[target][...], dtype=np.float64)
            require(row.shape == (64,) and np.all(np.isfinite(row)), f"{target}: invalid STRING row")
            values[index] = row
    return values


def pan_covered_genes() -> set[str]:
    result: set[str] = set()
    with (PAN_ROOT / "source_gene_map.tsv").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(reader.fieldnames is not None, "Pan source-gene map has no header")
        for row in reader:
            if row["mapped"] == "1":
                result.add(row["canonical_gene"])
    require(len(result) == 16_290, "Pan covered-gene count changed")
    return result


def consumed_targets() -> tuple[set[str], set[str]]:
    plan = read_json(EXP107_PLAN)
    used107 = set().union(*(set(gallery["targets"]) for gallery in plan["galleries"]))
    used123 = set(read_lines(EXP123_TARGETS))
    require(len(used107) == 1800, "Experiment 107 consumed-target count changed")
    require(len(used123) == 600, "Experiment 123 consumed-target count changed")
    return used107, used123


def source_roots() -> dict[str, Path]:
    return {
        "k562_gwps": ROOT / "data/donors/k562_gwps",
        "orion_hct116": ROOT / "data/donors/orion_hct116",
        "orion_hek293t": ROOT / "data/donors/orion_hek293t",
        "rpe1": ROOT / "data/donors/rpe1",
        "nadig_jurkat": ROOT / "data/donors/nadig_jurkat",
    }


def create_plan() -> None:
    require(not PLAN_DIR.exists(), f"{PLAN_DIR} exists; refusing overwrite")
    require(not AUTHORIZATION_PATH.exists(), "authorization already exists")
    require(not SELECTION_DIR.exists() and not RESULT_DIR.exists(), "later-stage output already exists")

    panel = set(read_csv_column(PANEL_PATH, "target_gene"))
    official = read_csv_column(GENES_PATH, "gene_name")
    roots = source_roots()
    banks = {name: Bank.load_axis(name, path) for name, path in roots.items()}
    pan_targets = set(read_lines(PAN_ROOT / "perts.txt"))
    with h5py.File(STRING_PATH, "r") as handle:
        string_targets = set(map(str, handle.keys()))
    used107, used123 = consumed_targets()

    eligible = set.intersection(
        pan_targets,
        *(set(banks[name].perts) for name in (*BASE_SOURCES, *OUTCOMES)),
        string_targets,
    )
    eligible -= panel | used107 | used123
    require(len(eligible) == EXPECTED_TARGETS, f"fresh target count {len(eligible)} != {EXPECTED_TARGETS}")
    ordered = stable_order(eligible, TARGET_NAMESPACE)
    selection_targets = ordered[:SELECTION_COUNT]
    test_targets = ordered[SELECTION_COUNT:]
    require(len(test_targets) == TEST_COUNT, "test target count changed")

    common_genes = pan_covered_genes()
    for name in (*BASE_SOURCES, *OUTCOMES):
        common_genes &= set(banks[name].genes)
    score_genes = tuple(
        gene
        for gene in official
        if gene in common_genes and gene not in panel and gene not in eligible
    )
    require(len(score_genes) == EXPECTED_SCORE_GENES, "score-gene count changed")
    require(sha256_lines(score_genes) == EXPECTED_SCORE_SHA, "score-gene axis changed")

    folds = deterministic_folds(selection_targets)
    rng = np.random.default_rng(NULL_SEED)
    selection_permutations = np.stack(
        [derangement(rng, SELECTION_COUNT) for _ in range(NULL_DRAWS)]
    )
    test_permutations = np.stack([derangement(rng, TEST_COUNT) for _ in range(NULL_DRAWS)])

    test_embeddings = load_string_rows(test_targets)
    cluster_labels = reciprocal_knn_clusters(test_embeddings)
    cluster_count = int(cluster_labels.max()) + 1
    bootstrap_rng = np.random.default_rng(BOOTSTRAP_SEED)
    bootstrap_draws = bootstrap_rng.integers(
        0,
        cluster_count,
        size=(BOOTSTRAP_DRAWS, cluster_count),
        dtype=np.uint16,
    )

    PLAN_DIR.mkdir(parents=True)
    write_lines(PLAN_DIR / "all_targets.txt", ordered)
    write_lines(PLAN_DIR / "selection_targets.txt", selection_targets)
    write_lines(PLAN_DIR / "test_targets.txt", test_targets)
    write_lines(PLAN_DIR / "score_genes.txt", score_genes)
    np.save(PLAN_DIR / "selection_folds.npy", folds, allow_pickle=False)
    np.save(PLAN_DIR / "selection_target_permutations.npy", selection_permutations, allow_pickle=False)
    np.save(PLAN_DIR / "test_target_permutations.npy", test_permutations, allow_pickle=False)
    np.save(PLAN_DIR / "test_cluster_labels.npy", cluster_labels, allow_pickle=False)
    np.save(PLAN_DIR / "bootstrap_cluster_draws.npy", bootstrap_draws, allow_pickle=False)

    files = {
        path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
        for path in sorted(PLAN_DIR.iterdir())
        if path.name != "plan-manifest.json"
    }
    manifest = {
        "schema_version": 1,
        "experiment": "129-gse295214-transfer",
        "created_before_source_effect_access": True,
        "eligible_target_count": len(eligible),
        "selection_target_count": len(selection_targets),
        "test_target_count": len(test_targets),
        "selection_pan_supported_count": sum(target in pan_targets for target in selection_targets),
        "test_pan_supported_count": sum(target in pan_targets for target in test_targets),
        "score_gene_count": len(score_genes),
        "test_cluster_count": cluster_count,
        "candidate_count": len(candidates()),
        "target_order_sha256": sha256_lines(ordered),
        "selection_target_order_sha256": sha256_lines(selection_targets),
        "test_target_order_sha256": sha256_lines(test_targets),
        "score_gene_order_sha256": sha256_lines(score_genes),
        "files": files,
    }
    write_json(PLAN_DIR / "plan-manifest.json", manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))


def inventory_paths() -> dict[str, str]:
    inventory = read_json(INVENTORY_PATH)
    result: dict[str, str] = {}
    for section in ("effect_inputs", "axis_and_metadata_inputs", "production_reference_inputs"):
        result.update({str(path): str(digest) for path, digest in inventory[section].items()})
    require(all(len(value) == 64 for value in result.values()), "inventory contains an unfrozen hash")
    return result


def validate_plan() -> dict[str, Any]:
    require(PLAN_DIR.is_dir(), "plan directory absent")
    manifest = read_json(PLAN_DIR / "plan-manifest.json")
    require(manifest["candidate_count"] == len(candidates()), "candidate count drifted")
    for name, expected in manifest["files"].items():
        path = PLAN_DIR / name
        require(path.is_file(), f"missing plan file {name}")
        require(path.stat().st_size == expected["bytes"], f"{name}: byte count changed")
        require(sha256_file(path) == expected["sha256"], f"{name}: hash changed")
    return manifest


def freeze_authorization() -> None:
    require(not AUTHORIZATION_PATH.exists(), "authorization already exists")
    require(not SELECTION_DIR.exists() and not RESULT_DIR.exists(), "effect-stage output exists")
    require(not TEST_ACCESS_LOCK.exists(), "source-test access lock already exists")
    require(AUTHORIZATION_V1_PATH.is_file(), "predecessor authorization absent")
    require(RECOVERY_001_PATH.is_file(), "recovery receipt 001 absent")
    require(RECOVERY_002_PATH.is_file(), "recovery receipt 002 absent")
    plan = validate_plan()
    for relative, expected in inventory_paths().items():
        path = ROOT / relative
        require(path.is_file(), f"missing frozen input {relative}")
        require(sha256_file(path) == expected, f"input hash changed: {relative}")
    predecessor = read_json(AUTHORIZATION_V1_PATH)
    recovery = read_json(RECOVERY_002_PATH)
    predecessor_sha256 = sha256_file(AUTHORIZATION_V1_PATH)
    require(
        predecessor_sha256
        == recovery["failed_implementation"]["SOURCE_ONLY_AUTHORIZATION.json_sha256"],
        "predecessor authorization does not match recovery receipt",
    )
    require(predecessor["source_only_authorized"] is True, "predecessor source access denied")
    require(predecessor["held_outcome_access_authorized"] is False, "predecessor held access enabled")
    require(recovery["access_boundary"]["untouched_source_test_effect_values_loaded"] is False, "recovery reports source-test access")
    require(recovery["result_boundary"]["accepted_result_produced"] is False, "recovery reports an accepted result")
    implementation = {
        str(path.relative_to(ROOT)): sha256_file(path)
        for path in (
            HERE / "README.md",
            PROTOCOL_PATH,
            INVENTORY_PATH,
            AUTHORIZATION_V1_PATH,
            RECOVERY_001_PATH,
            RECOVERY_002_PATH,
            HERE / "run.py",
            HERE / "test_run.py",
        )
    }
    authorization = {
        "schema_version": 1,
        "experiment": "129-gse295214-transfer",
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "authorization_revision": 2,
        "recovery_authorization": True,
        "created_before_source_effect_loading": False,
        "source_effect_loading_already_occurred": True,
        "created_before_source_test_effect_loading": True,
        "created_before_any_accepted_result": True,
        "source_only_authorized": True,
        "held_outcome_access_authorized": False,
        "cell_eval2_authorized": False,
        "count_generation_authorized": False,
        "challenge_build_or_submission_authorized": False,
        "commit_or_push_authorized": False,
        "predecessor_authorization_sha256": predecessor_sha256,
        "recovery_receipts": {
            "IMPLEMENTATION_RECOVERY_001.json": sha256_file(RECOVERY_001_PATH),
            "IMPLEMENTATION_RECOVERY_002.json": sha256_file(RECOVERY_002_PATH),
        },
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "inventory_sha256": sha256_file(INVENTORY_PATH),
        "plan_manifest_sha256": sha256_file(PLAN_DIR / "plan-manifest.json"),
        "implementation_files": implementation,
        "input_count": len(inventory_paths()),
        "plan_summary": {
            "eligible_target_count": plan["eligible_target_count"],
            "selection_target_count": plan["selection_target_count"],
            "test_target_count": plan["test_target_count"],
            "score_gene_count": plan["score_gene_count"],
            "candidate_count": plan["candidate_count"],
        },
    }
    write_json(AUTHORIZATION_PATH, authorization)
    print(json.dumps(authorization, indent=2, sort_keys=True))


def validate_authorization(*, hash_large_inputs: bool) -> dict[str, Any]:
    require(AUTHORIZATION_PATH.is_file(), "source-only authorization absent")
    authorization = read_json(AUTHORIZATION_PATH)
    require(authorization["source_only_authorized"] is True, "source-only access denied")
    require(authorization["authorization_revision"] == 2, "recovery authorization revision changed")
    require(authorization["recovery_authorization"] is True, "recovery authorization absent")
    require(authorization["created_before_source_test_effect_loading"] is True, "source-test boundary invalid")
    require(authorization["created_before_any_accepted_result"] is True, "accepted-result boundary invalid")
    require(authorization["held_outcome_access_authorized"] is False, "held access unexpectedly enabled")
    require(authorization["cell_eval2_authorized"] is False, "cell-eval2 unexpectedly enabled")
    require(
        authorization["challenge_build_or_submission_authorized"] is False,
        "challenge build unexpectedly enabled",
    )
    require(authorization["protocol_sha256"] == sha256_file(PROTOCOL_PATH), "protocol changed")
    require(authorization["inventory_sha256"] == sha256_file(INVENTORY_PATH), "inventory changed")
    require(
        authorization["plan_manifest_sha256"] == sha256_file(PLAN_DIR / "plan-manifest.json"),
        "plan manifest changed",
    )
    require(
        authorization["predecessor_authorization_sha256"]
        == sha256_file(AUTHORIZATION_V1_PATH),
        "predecessor authorization changed",
    )
    for name, expected in authorization["recovery_receipts"].items():
        require(sha256_file(HERE / name) == expected, f"recovery receipt changed: {name}")
    for relative, expected in authorization["implementation_files"].items():
        require(sha256_file(ROOT / relative) == expected, f"implementation changed: {relative}")
    if hash_large_inputs:
        for relative, expected in inventory_paths().items():
            require(sha256_file(ROOT / relative) == expected, f"input hash changed: {relative}")
    validate_plan()
    return authorization


def bank_indices(bank: Bank, targets: Sequence[str], genes: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    pert_index = {label: index for index, label in enumerate(bank.perts)}
    gene_index = {label: index for index, label in enumerate(bank.genes)}
    require(all(target in pert_index for target in targets), f"{bank.name}: target missing")
    require(all(gene in gene_index for gene in genes), f"{bank.name}: response gene missing")
    return (
        np.asarray([pert_index[target] for target in targets], dtype=np.int64),
        np.asarray([gene_index[gene] for gene in genes], dtype=np.int64),
    )


def load_rows(bank: Bank, targets: Sequence[str], genes: Sequence[str]) -> np.ndarray:
    rows, columns = bank_indices(bank, targets, genes)
    matrix = np.load(bank.delta_path, mmap_mode="r")
    values = np.asarray(matrix[rows][:, columns], dtype=np.float64)
    require(values.shape == (len(targets), len(genes)), f"{bank.name}: sliced shape changed")
    require(np.all(np.isfinite(values)), f"{bank.name}: non-finite values")
    return values


def apply_measured_mask(values: np.ndarray, support: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    measured = np.asarray(support, dtype=bool)
    require(array.shape == measured.shape, "Pan values and measured mask do not align")
    require(np.all(np.isfinite(array)), "Pan rows contain non-finite values")
    require(
        np.count_nonzero(array[~measured]) == 0,
        "Pan unsupported delta is not exact zero",
    )
    return np.where(measured, array, 0.0)


def load_pan_rows(targets: Sequence[str], score_genes: Sequence[str]) -> np.ndarray:
    perts = read_lines(PAN_ROOT / "perts.txt")
    genes = read_lines(PAN_ROOT / "genes.txt")
    pert_index = {label: index for index, label in enumerate(perts)}
    gene_index = {label: index for index, label in enumerate(genes)}
    rows = np.asarray([pert_index[target] for target in targets], dtype=np.int64)
    columns = np.asarray([gene_index[gene] for gene in score_genes], dtype=np.int64)
    delta = np.load(PAN_ROOT / "delta.npy", mmap_mode="r")
    measured = np.load(PAN_ROOT / "measured_mask.npy", mmap_mode="r")
    support = np.asarray(measured[rows][:, columns], dtype=bool)
    values = np.asarray(delta[rows][:, columns], dtype=np.float64)
    return apply_measured_mask(values, support)


def source_row_order(bank: Bank, excluded: set[str]) -> np.ndarray:
    eligible = [
        index
        for index, target in enumerate(bank.perts)
        if target not in excluded
    ]
    ordered = sorted(
        eligible,
        key=lambda index: (
            hashlib.sha256(
                BASIS_NAMESPACE.encode("utf-8")
                + b"\0"
                + bank.name.encode("utf-8")
                + b"\0"
                + bank.perts[index].encode("utf-8")
            ).digest(),
            bank.perts[index].encode("utf-8"),
        ),
    )
    return np.asarray(ordered[:BASIS_CAP], dtype=np.int64)


def source_mean_and_scale(
    bank: Bank,
    score_genes: Sequence[str],
    excluded: set[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[np.ndarray, float]:
    _, columns = bank_indices(bank, (), score_genes)
    matrix = np.load(bank.delta_path, mmap_mode="r")
    included = np.asarray([target not in excluded for target in bank.perts], dtype=bool)
    total = np.zeros(len(columns), dtype=np.float64)
    count = 0
    for lo in range(0, len(bank.perts), chunk_rows):
        hi = min(lo + chunk_rows, len(bank.perts))
        local = included[lo:hi]
        if not np.any(local):
            continue
        block = np.asarray(matrix[lo:hi], dtype=np.float32)[local][:, columns]
        total += block.sum(axis=0, dtype=np.float64)
        count += int(local.sum())
        memory_checkpoint(f"{bank.name} source mean", max_rss_gib)
    require(count > 0, f"{bank.name}: no rows for source mean")
    mean = total / count
    norms: list[np.ndarray] = []
    for lo in range(0, len(bank.perts), chunk_rows):
        hi = min(lo + chunk_rows, len(bank.perts))
        local = included[lo:hi]
        if not np.any(local):
            continue
        block = np.asarray(matrix[lo:hi], dtype=np.float32)[local][:, columns]
        norms.append(np.linalg.norm(block.astype(np.float64) - mean, axis=1))
    scale = float(np.median(np.concatenate(norms)))
    require(math.isfinite(scale) and scale > 0, f"{bank.name}: invalid source scale")
    return mean, scale


def pm_mean_and_scale(
    score_genes: Sequence[str],
    excluded: set[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[np.ndarray, np.ndarray, float, Bank]:
    bank = Bank.load_axis("perturbmatch", PM_ROOT)
    common = [gene for gene in score_genes if gene in set(bank.genes)]
    score_pos = {gene: index for index, gene in enumerate(score_genes)}
    destination = np.asarray([score_pos[gene] for gene in common], dtype=np.int64)
    mean, scale = source_mean_and_scale(
        bank,
        common,
        excluded,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    return destination, mean, scale, bank


def fit_source_geometry(
    score_genes: Sequence[str],
    excluded_targets: Sequence[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> dict[str, Any]:
    excluded = set(excluded_targets)
    banks = {
        name: Bank.load_axis(name, ROOT / "data/donors" / name)
        for name in BASE_SOURCES
    }
    means: dict[str, np.ndarray] = {}
    scales: dict[str, float] = {}
    for name, bank in banks.items():
        means[name], scales[name] = source_mean_and_scale(
            bank,
            score_genes,
            excluded,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
    reference_scale = float(np.median(list(scales.values())))
    require(reference_scale > 0, "invalid reference source scale")

    basis_rows = {name: source_row_order(bank, excluded) for name, bank in banks.items()}
    total_rows = sum(len(rows) for rows in basis_rows.values())
    basis_path = SELECTION_DIR.parent / ".experiment-129-source-basis.npy"
    require(not basis_path.exists(), f"stale scratch basis exists: {basis_path}")
    basis = np.lib.format.open_memmap(
        basis_path,
        mode="w+",
        dtype=np.float32,
        shape=(total_rows, len(score_genes)),
    )
    offset = 0
    for name, bank in banks.items():
        _, columns = bank_indices(bank, (), score_genes)
        matrix = np.load(bank.delta_path, mmap_mode="r")
        rows = basis_rows[name]
        multiplier = reference_scale / scales[name]
        for lo in range(0, len(rows), chunk_rows):
            hi = min(lo + chunk_rows, len(rows))
            selected = rows[lo:hi]
            block = np.asarray(matrix[selected], dtype=np.float32)[:, columns]
            basis[offset + lo : offset + hi] = (
                block.astype(np.float64) - means[name]
            ) * multiplier
            memory_checkpoint(f"{name} whitener basis", max_rss_gib)
        offset += len(rows)
    basis.flush()
    print(f"[whitener] basis={basis.shape}", flush=True)
    covariance = (basis.T @ basis).astype(np.float64) / total_rows
    del basis
    basis_path.unlink()
    evals, evecs = np.linalg.eigh(covariance)
    del covariance
    floor = max(float(evals[-1]) * 1e-12, np.finfo(np.float64).tiny)
    evals = np.maximum(evals, floor)
    weights = evals ** (-WHITEN_GAMMA / 2.0)
    memory_checkpoint("whitener eigendecomposition", max_rss_gib)

    pm_destination, pm_mean, pm_scale, _ = pm_mean_and_scale(
        score_genes,
        excluded,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    return {
        "means": means,
        "scales": scales,
        "reference_scale": reference_scale,
        "basis_rows": basis_rows,
        "evals": evals,
        "evecs": np.asarray(evecs, dtype=np.float32),
        "weights": weights,
        "pm_destination": pm_destination,
        "pm_mean": pm_mean,
        "pm_scale": pm_scale,
    }


def standardized_rows(
    bank: Bank,
    targets: Sequence[str],
    score_genes: Sequence[str],
    mean: np.ndarray,
    scale: float,
    reference_scale: float,
) -> np.ndarray:
    return (load_rows(bank, targets, score_genes) - mean) * (reference_scale / scale)


def build_incumbent_pre(
    targets: Sequence[str],
    score_genes: Sequence[str],
    geometry: Mapping[str, Any],
) -> np.ndarray:
    base_banks = {
        name: Bank.load_axis(name, ROOT / "data/donors" / name)
        for name in BASE_SOURCES
    }
    rows = {
        name: standardized_rows(
            bank,
            targets,
            score_genes,
            np.asarray(geometry["means"][name]),
            float(geometry["scales"][name]),
            float(geometry["reference_scale"]),
        )
        for name, bank in base_banks.items()
    }
    kslot = rows["k562_gwps"].copy()
    pm_bank = Bank.load_axis("perturbmatch", PM_ROOT)
    pm_pert_index = {target: index for index, target in enumerate(pm_bank.perts)}
    pm_targets = [target for target in targets if target in pm_pert_index]
    if pm_targets:
        pm_common = [score_genes[index] for index in np.asarray(geometry["pm_destination"], dtype=np.int64)]
        pm_values = standardized_rows(
            pm_bank,
            pm_targets,
            pm_common,
            np.asarray(geometry["pm_mean"]),
            float(geometry["pm_scale"]),
            float(geometry["reference_scale"]),
        )
        target_pos = {target: index for index, target in enumerate(targets)}
        destination = np.asarray(geometry["pm_destination"], dtype=np.int64)
        for source_row, target in enumerate(pm_targets):
            row = target_pos[target]
            kslot[row, destination] = (
                (1.0 - PM_BLEND) * kslot[row, destination]
                + PM_BLEND * pm_values[source_row]
            )
    return (kslot + rows["orion_hct116"] + rows["orion_hek293t"]) / 3.0


def apply_whitener(matrix: np.ndarray, geometry: Mapping[str, Any]) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64)
    before = float(np.median(np.linalg.norm(values, axis=1))) or 1.0
    evecs = np.asarray(geometry["evecs"], dtype=np.float64)
    transformed = ((values @ evecs) * np.asarray(geometry["weights"])) @ evecs.T
    after = float(np.median(np.linalg.norm(transformed, axis=1))) or 1.0
    return transformed * (before / after)


def truth_rows(targets: Sequence[str], score_genes: Sequence[str]) -> dict[str, np.ndarray]:
    return {
        name: load_rows(
            Bank.load_axis(name, ROOT / "data/donors" / name),
            targets,
            score_genes,
        )
        for name in OUTCOMES
    }


def tangent_residual(incumbent: np.ndarray, truth: np.ndarray) -> np.ndarray:
    base = row_unit(incumbent)
    observed = row_unit(truth)
    projection = np.sum(observed * base, axis=1, keepdims=True)
    return observed - projection * base


def knn_tangent_prediction(
    incumbent: np.ndarray,
    query_pan: np.ndarray,
    train_pan: np.ndarray,
    teacher_residual: np.ndarray,
    k: int,
) -> tuple[np.ndarray, np.ndarray]:
    base = np.asarray(incumbent, dtype=np.float64)
    query = row_unit(query_pan)
    train = row_unit(train_pan)
    residuals = np.asarray(teacher_residual, dtype=np.float64)
    require(train.shape[0] == residuals.shape[0], "training Pan/residual rows do not align")
    similarities = query @ train.T
    output = base.copy()
    used = np.zeros(len(base), dtype=bool)
    base_unit = row_unit(base)
    base_norm = np.linalg.norm(base, axis=1)
    for row in range(len(base)):
        order = np.argsort(-similarities[row], kind="stable")
        chosen = order[similarities[row, order] > 0][:k]
        if not len(chosen) or base_norm[row] <= EPS:
            continue
        tangent = residuals[chosen].mean(axis=0)
        tangent -= float(np.dot(tangent, base_unit[row])) * base_unit[row]
        norm = float(np.linalg.norm(tangent))
        if norm <= EPS or not math.isfinite(norm):
            continue
        direction = math.cos(ANGLE_RADIANS) * base_unit[row]
        direction += math.sin(ANGLE_RADIANS) * (tangent / norm)
        direction_norm = float(np.linalg.norm(direction))
        if direction_norm <= EPS or not math.isfinite(direction_norm):
            continue
        output[row] = direction * (base_norm[row] / direction_norm)
        used[row] = True
    return output, used


def pds_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    predicted = np.asarray(prediction, dtype=np.float64)
    observed = np.asarray(truth, dtype=np.float64)
    require(predicted.shape == observed.shape and predicted.ndim == 2, "PDS arrays do not align")
    require(len(predicted) >= 2, "PDS gallery too small")
    truth_unit = row_unit(observed)
    result = np.full(len(predicted), 0.5, dtype=np.float64)
    norms = np.linalg.norm(predicted, axis=1)
    rows = np.flatnonzero(norms > EPS)
    if len(rows):
        scores = predicted[rows] / norms[rows, None] @ truth_unit.T
        own = scores[np.arange(len(rows)), rows][:, None]
        rank = (scores > own).sum(axis=1).astype(np.float64)
        rank += 0.5 * ((scores == own).sum(axis=1) - 1)
        result[rows] = 1.0 - rank / (len(predicted) - 1)
    return result


def signed_credit_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    pred = np.asarray(prediction, dtype=np.float64)
    observed = np.asarray(truth, dtype=np.float64)
    denominator = np.abs(pred).sum(axis=1)
    numerator = (np.sign(observed) * pred).sum(axis=1)
    return np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > EPS)


def metric_bundle(prediction: np.ndarray, truth: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "pds": pds_per_target(prediction, truth),
        "signed_credit": signed_credit_per_target(prediction, truth),
    }


def selection_predictions(
    baseline: np.ndarray,
    pan: np.ndarray,
    truths: Mapping[str, np.ndarray],
    folds: np.ndarray,
    candidate: Candidate,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    predictions = {name: baseline.copy() for name in OUTCOMES}
    used = {name: np.zeros(len(baseline), dtype=bool) for name in OUTCOMES}
    if candidate.kind == "identity":
        return predictions, used
    for fold in range(FOLDS):
        train = folds != fold
        query = folds == fold
        for outcome, teacher in (("rpe1", "nadig_jurkat"), ("nadig_jurkat", "rpe1")):
            teacher_residual = tangent_residual(baseline[train], truths[teacher][train])
            predicted, local_used = knn_tangent_prediction(
                baseline[query],
                pan[query],
                pan[train],
                teacher_residual,
                candidate.k,
            )
            predictions[outcome][query] = predicted
            used[outcome][query] = local_used
    return predictions, used


def test_predictions(
    selection_baseline: np.ndarray,
    selection_pan: np.ndarray,
    selection_truths: Mapping[str, np.ndarray],
    test_baseline: np.ndarray,
    test_pan: np.ndarray,
    candidate: Candidate,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    predictions = {name: test_baseline.copy() for name in OUTCOMES}
    used = {name: np.zeros(len(test_baseline), dtype=bool) for name in OUTCOMES}
    if candidate.kind == "identity":
        return predictions, used
    for outcome, teacher in (("rpe1", "nadig_jurkat"), ("nadig_jurkat", "rpe1")):
        residual = tangent_residual(selection_baseline, selection_truths[teacher])
        predictions[outcome], used[outcome] = knn_tangent_prediction(
            test_baseline,
            test_pan,
            selection_pan,
            residual,
            candidate.k,
        )
    return predictions, used


def summarize_predictions(
    baseline: np.ndarray,
    predictions: Mapping[str, np.ndarray],
    truths: Mapping[str, np.ndarray],
    used: Mapping[str, np.ndarray],
) -> tuple[dict[str, Any], dict[str, dict[str, np.ndarray]]]:
    summaries: dict[str, Any] = {}
    per_context: dict[str, dict[str, np.ndarray]] = {}
    for outcome in OUTCOMES:
        base_metrics = metric_bundle(baseline, truths[outcome])
        candidate_metrics = metric_bundle(predictions[outcome], truths[outcome])
        arrays = {
            f"{metric}_baseline": base_metrics[metric]
            for metric in base_metrics
        }
        arrays.update(
            {f"{metric}_candidate": candidate_metrics[metric] for metric in candidate_metrics}
        )
        arrays.update(
            {
                f"{metric}_delta": candidate_metrics[metric] - base_metrics[metric]
                for metric in candidate_metrics
            }
        )
        arrays["used"] = np.asarray(used[outcome], dtype=bool)
        per_context[outcome] = arrays
        summaries[outcome] = {
            "baseline": {metric: float(values.mean()) for metric, values in base_metrics.items()},
            "candidate": {
                metric: float(values.mean()) for metric, values in candidate_metrics.items()
            },
            "delta": {
                metric: float((candidate_metrics[metric] - base_metrics[metric]).mean())
                for metric in candidate_metrics
            },
            "used_targets": int(np.sum(used[outcome])),
            "fallback_targets": int(np.sum(~used[outcome])),
            "fallback_bit_exact": bool(
                np.array_equal(predictions[outcome][~used[outcome]], baseline[~used[outcome]])
            ),
        }
    macro = {
        metric: float(np.mean([summaries[name]["delta"][metric] for name in OUTCOMES]))
        for metric in ("pds", "signed_credit")
    }
    return {"contexts": summaries, "macro_delta": macro}, per_context


def select_candidate(
    baseline: np.ndarray,
    pan: np.ndarray,
    truths: Mapping[str, np.ndarray],
    folds: np.ndarray,
) -> tuple[Candidate, dict[str, Any]]:
    table: list[dict[str, Any]] = []
    winner = candidates()[0]
    winner_key = (1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    for candidate in candidates():
        predictions, used = selection_predictions(baseline, pan, truths, folds, candidate)
        summary, _ = summarize_predictions(baseline, predictions, truths, used)
        pds = [summary["contexts"][name]["delta"]["pds"] for name in OUTCOMES]
        credit = [summary["contexts"][name]["delta"]["signed_credit"] for name in OUTCOMES]
        eligible = candidate.kind == "identity" or (min(pds) >= 0.0 and min(credit) >= 0.0)
        key = (
            float(eligible),
            float(min(pds)),
            float(np.mean(pds)),
            float(min(credit)),
            float(np.mean(credit)),
            -float(candidate.k),
        )
        if key > winner_key:
            winner, winner_key = candidate, key
        table.append(
            {
                "candidate": candidate.key,
                "eligible": bool(eligible),
                "selection_summary": summary,
                "selection_key": list(key),
            }
        )
    return winner, {
        "selected": winner.key,
        "selected_key": list(winner_key),
        "candidates": table,
    }


def serialize_geometry(path: Path, geometry: Mapping[str, Any]) -> None:
    arrays: dict[str, np.ndarray] = {
        "evals": np.asarray(geometry["evals"], dtype=np.float64),
        "evecs": np.asarray(geometry["evecs"], dtype=np.float32),
        "weights": np.asarray(geometry["weights"], dtype=np.float64),
        "pm_destination": np.asarray(geometry["pm_destination"], dtype=np.int64),
        "pm_mean": np.asarray(geometry["pm_mean"], dtype=np.float64),
        "pm_scale": np.asarray([geometry["pm_scale"]], dtype=np.float64),
        "reference_scale": np.asarray([geometry["reference_scale"]], dtype=np.float64),
    }
    for name in BASE_SOURCES:
        arrays[f"{name}__mean"] = np.asarray(geometry["means"][name], dtype=np.float64)
        arrays[f"{name}__scale"] = np.asarray([geometry["scales"][name]], dtype=np.float64)
        arrays[f"{name}__basis_rows"] = np.asarray(geometry["basis_rows"][name], dtype=np.int64)
    np.savez(path, **arrays)


def deserialize_geometry(path: Path) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as handle:
        return {
            "evals": np.asarray(handle["evals"], dtype=np.float64),
            "evecs": np.asarray(handle["evecs"], dtype=np.float32),
            "weights": np.asarray(handle["weights"], dtype=np.float64),
            "pm_destination": np.asarray(handle["pm_destination"], dtype=np.int64),
            "pm_mean": np.asarray(handle["pm_mean"], dtype=np.float64),
            "pm_scale": float(handle["pm_scale"][0]),
            "reference_scale": float(handle["reference_scale"][0]),
            "means": {
                name: np.asarray(handle[f"{name}__mean"], dtype=np.float64)
                for name in BASE_SOURCES
            },
            "scales": {
                name: float(handle[f"{name}__scale"][0])
                for name in BASE_SOURCES
            },
            "basis_rows": {
                name: np.asarray(handle[f"{name}__basis_rows"], dtype=np.int64)
                for name in BASE_SOURCES
            },
        }


def save_cache(path: Path, baseline: np.ndarray, pan: np.ndarray, truths: Mapping[str, np.ndarray]) -> None:
    np.savez_compressed(
        path,
        baseline=np.asarray(baseline, dtype=np.float64),
        pan=np.asarray(pan, dtype=np.float64),
        **{f"truth__{name}": np.asarray(truths[name], dtype=np.float64) for name in OUTCOMES},
    )


def load_cache(path: Path) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    with np.load(path, allow_pickle=False) as handle:
        return (
            np.asarray(handle["baseline"], dtype=np.float64),
            np.asarray(handle["pan"], dtype=np.float64),
            {
                name: np.asarray(handle[f"truth__{name}"], dtype=np.float64)
                for name in OUTCOMES
            },
        )


def selection_stage(args: argparse.Namespace) -> None:
    validate_authorization(hash_large_inputs=True)
    require(not SELECTION_DIR.exists(), f"{SELECTION_DIR} exists; refusing overwrite")
    require(not RESULT_DIR.exists() and not TEST_ACCESS_LOCK.exists(), "test stage already exists")
    all_targets = read_lines(PLAN_DIR / "all_targets.txt")
    selection_targets = read_lines(PLAN_DIR / "selection_targets.txt")
    score_genes = read_lines(PLAN_DIR / "score_genes.txt")
    folds = np.load(PLAN_DIR / "selection_folds.npy", allow_pickle=False)

    SELECTION_DIR.parent.mkdir(parents=True, exist_ok=True)
    geometry = fit_source_geometry(
        score_genes,
        all_targets,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    baseline = apply_whitener(
        build_incumbent_pre(selection_targets, score_genes, geometry),
        geometry,
    )
    pan = load_pan_rows(selection_targets, score_genes)
    truths = truth_rows(selection_targets, score_genes)
    selected, audit = select_candidate(baseline, pan, truths, folds)

    SELECTION_DIR.mkdir()
    geometry_path = SELECTION_DIR / "geometry.npz"
    cache_path = SELECTION_DIR / "selection-cache.npz"
    serialize_geometry(geometry_path, geometry)
    save_cache(cache_path, baseline, pan, truths)
    selection = {
        "schema_version": 1,
        "experiment": "129-gse295214-transfer",
        "sealed_before_source_test_effect_access": True,
        "selected": asdict(selected) | {"key": selected.key},
        "selection_audit": audit,
        "selection_target_order_sha256": sha256_lines(selection_targets),
        "score_gene_order_sha256": sha256_lines(score_genes),
        "geometry_sha256": sha256_file(geometry_path),
        "selection_cache_sha256": sha256_file(cache_path),
        "peak_rss_gib": peak_rss_gib(),
    }
    selection_path = SELECTION_DIR / "selection.json"
    write_json(selection_path, selection)
    manifest = {
        "schema_version": 1,
        "experiment": "129-gse295214-transfer",
        "selection_sha256": sha256_file(selection_path),
        "geometry_sha256": sha256_file(geometry_path),
        "selection_cache_sha256": sha256_file(cache_path),
        "authorization_sha256": sha256_file(AUTHORIZATION_PATH),
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "plan_manifest_sha256": sha256_file(PLAN_DIR / "plan-manifest.json"),
        "implementation_sha256": sha256_file(Path(__file__)),
    }
    write_json(SELECTION_DIR / "freeze-manifest.json", manifest)
    print(json.dumps(selection, indent=2, sort_keys=True))


def verify_selection() -> tuple[dict[str, Any], dict[str, Any], tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]]:
    validate_authorization(hash_large_inputs=False)
    manifest = read_json(SELECTION_DIR / "freeze-manifest.json")
    require(manifest["authorization_sha256"] == sha256_file(AUTHORIZATION_PATH), "authorization changed")
    require(manifest["protocol_sha256"] == sha256_file(PROTOCOL_PATH), "protocol changed")
    require(manifest["plan_manifest_sha256"] == sha256_file(PLAN_DIR / "plan-manifest.json"), "plan changed")
    require(manifest["implementation_sha256"] == sha256_file(Path(__file__)), "implementation changed")
    selection_path = SELECTION_DIR / "selection.json"
    geometry_path = SELECTION_DIR / "geometry.npz"
    cache_path = SELECTION_DIR / "selection-cache.npz"
    require(manifest["selection_sha256"] == sha256_file(selection_path), "selection changed")
    require(manifest["geometry_sha256"] == sha256_file(geometry_path), "geometry changed")
    require(manifest["selection_cache_sha256"] == sha256_file(cache_path), "selection cache changed")
    selection = read_json(selection_path)
    require(selection["sealed_before_source_test_effect_access"] is True, "selection seal invalid")
    return selection, deserialize_geometry(geometry_path), load_cache(cache_path)


def candidate_from_json(value: Mapping[str, Any]) -> Candidate:
    return Candidate(str(value["kind"]), int(value["k"]))


def null_distribution(
    selection_baseline: np.ndarray,
    selection_pan: np.ndarray,
    selection_truths: Mapping[str, np.ndarray],
    test_baseline: np.ndarray,
    test_pan: np.ndarray,
    test_truths: Mapping[str, np.ndarray],
    folds: np.ndarray,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    selection_permutations = np.load(
        PLAN_DIR / "selection_target_permutations.npy", allow_pickle=False
    )
    test_permutations = np.load(PLAN_DIR / "test_target_permutations.npy", allow_pickle=False)
    pds = np.empty(NULL_DRAWS, dtype=np.float64)
    credit = np.empty(NULL_DRAWS, dtype=np.float64)
    selected_keys: list[str] = []
    for draw in range(NULL_DRAWS):
        perm_selection_pan = selection_pan[selection_permutations[draw]]
        selected, _ = select_candidate(
            selection_baseline,
            perm_selection_pan,
            selection_truths,
            folds,
        )
        selected_keys.append(selected.key)
        predictions, used = test_predictions(
            selection_baseline,
            perm_selection_pan,
            selection_truths,
            test_baseline,
            test_pan[test_permutations[draw]],
            selected,
        )
        summary, _ = summarize_predictions(test_baseline, predictions, test_truths, used)
        pds[draw] = summary["macro_delta"]["pds"]
        credit[draw] = summary["macro_delta"]["signed_credit"]
    return (
        {
            "draws": NULL_DRAWS,
            "selected_keys": selected_keys,
            "pds": {
                "p95": float(np.quantile(pds, 0.95, method="linear")),
                "maximum": float(np.max(pds)),
                "mean": float(np.mean(pds)),
            },
            "signed_credit": {
                "p95": float(np.quantile(credit, 0.95, method="linear")),
                "maximum": float(np.max(credit)),
                "mean": float(np.mean(credit)),
            },
        },
        {"null__pds": pds, "null__signed_credit": credit},
    )


def bootstrap_intervals(
    per_context: Mapping[str, Mapping[str, np.ndarray]],
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    labels = np.load(PLAN_DIR / "test_cluster_labels.npy", allow_pickle=False)
    draws = np.load(PLAN_DIR / "bootstrap_cluster_draws.npy", allow_pickle=False)
    cluster_count = int(labels.max()) + 1
    arrays = {
        "pds": np.empty(BOOTSTRAP_DRAWS, dtype=np.float64),
        "signed_credit": np.empty(BOOTSTRAP_DRAWS, dtype=np.float64),
    }
    for draw_index, sampled in enumerate(draws):
        multiplicity = np.bincount(sampled, minlength=cluster_count)
        weights = multiplicity[labels].astype(np.float64)
        denominator = float(weights.sum())
        require(denominator > 0, "empty bootstrap draw")
        for metric in arrays:
            context_means = [
                float(np.dot(per_context[name][f"{metric}_delta"], weights) / denominator)
                for name in OUTCOMES
            ]
            arrays[metric][draw_index] = float(np.mean(context_means))
    summary = {
        metric: {
            "lower_95": float(np.quantile(values, 0.025, method="linear")),
            "median": float(np.quantile(values, 0.5, method="linear")),
            "upper_95": float(np.quantile(values, 0.975, method="linear")),
        }
        for metric, values in arrays.items()
    }
    return summary, {f"bootstrap__{name}": values for name, values in arrays.items()}


def evaluate_stage(args: argparse.Namespace) -> None:
    validate_authorization(hash_large_inputs=True)
    selection, geometry, selection_cache = verify_selection()
    require(not RESULT_DIR.exists(), f"{RESULT_DIR} exists; refusing overwrite")
    require(not TEST_ACCESS_LOCK.exists(), "source-test access already consumed")
    test_targets = read_lines(PLAN_DIR / "test_targets.txt")
    score_genes = read_lines(PLAN_DIR / "score_genes.txt")
    folds = np.load(PLAN_DIR / "selection_folds.npy", allow_pickle=False)
    write_json(
        TEST_ACCESS_LOCK,
        {
            "schema_version": 1,
            "experiment": "129-gse295214-transfer",
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "meaning": "the 55 untouched source-test target effects are permanently consumed",
            "selection_freeze_sha256": sha256_file(SELECTION_DIR / "freeze-manifest.json"),
            "authorization_sha256": sha256_file(AUTHORIZATION_PATH),
            "implementation_sha256": sha256_file(Path(__file__)),
        },
    )
    selection_baseline, selection_pan, selection_truths = selection_cache
    test_baseline = apply_whitener(
        build_incumbent_pre(test_targets, score_genes, geometry),
        geometry,
    )
    test_pan = load_pan_rows(test_targets, score_genes)
    test_truths = truth_rows(test_targets, score_genes)
    selected = candidate_from_json(selection["selected"])
    predictions, used = test_predictions(
        selection_baseline,
        selection_pan,
        selection_truths,
        test_baseline,
        test_pan,
        selected,
    )
    observed, per_context = summarize_predictions(
        test_baseline,
        predictions,
        test_truths,
        used,
    )
    bootstrap, bootstrap_arrays = bootstrap_intervals(per_context)
    nulls, null_arrays = null_distribution(
        selection_baseline,
        selection_pan,
        selection_truths,
        test_baseline,
        test_pan,
        test_truths,
        folds,
    )
    selection_contexts = {
        name: next(
            row["selection_summary"]["contexts"][name]
            for row in selection["selection_audit"]["candidates"]
            if row["candidate"] == selected.key
        )
        for name in OUTCOMES
    }
    checks = {
        "selected_candidate_nonidentity": selected.kind != "identity",
        "selection_pds_delta_each_context_nonnegative": all(
            selection_contexts[name]["delta"]["pds"] >= 0.0 for name in OUTCOMES
        ),
        "selection_signed_credit_delta_each_context_nonnegative": all(
            selection_contexts[name]["delta"]["signed_credit"] >= 0.0 for name in OUTCOMES
        ),
        "test_pds_delta_each_context_nonnegative": all(
            observed["contexts"][name]["delta"]["pds"] >= 0.0 for name in OUTCOMES
        ),
        "test_signed_credit_delta_each_context_nonnegative": all(
            observed["contexts"][name]["delta"]["signed_credit"] >= 0.0 for name in OUTCOMES
        ),
        "test_macro_pds_delta_at_least_0p003": observed["macro_delta"]["pds"] >= 0.003,
        "test_macro_signed_credit_delta_positive": observed["macro_delta"]["signed_credit"] > 0.0,
        "pds_bootstrap_lower_positive": bootstrap["pds"]["lower_95"] > 0.0,
        "signed_credit_bootstrap_lower_positive": (
            bootstrap["signed_credit"]["lower_95"] > 0.0
        ),
        "pds_above_target_label_null_p95": (
            observed["macro_delta"]["pds"] > nulls["pds"]["p95"]
        ),
        "signed_credit_above_target_label_null_p95": (
            observed["macro_delta"]["signed_credit"] > nulls["signed_credit"]["p95"]
        ),
        "all_fallback_rows_bit_exact": all(
            observed["contexts"][name]["fallback_bit_exact"] for name in OUTCOMES
        ),
    }
    passed = all(checks.values())

    RESULT_DIR.mkdir(parents=True)
    arrays: dict[str, np.ndarray] = {
        "test_targets": np.asarray(test_targets, dtype="U"),
        "test_pan": test_pan,
        "test_baseline": test_baseline,
    }
    for name, values in per_context.items():
        for metric, array in values.items():
            arrays[f"{name}__{metric}"] = np.asarray(array)
    arrays.update(bootstrap_arrays)
    arrays.update(null_arrays)
    arrays_path = RESULT_DIR / "result-arrays.npz"
    np.savez_compressed(arrays_path, **arrays)
    result = {
        "schema_version": 1,
        "experiment": "129-gse295214-transfer",
        "status": "GO" if passed else "NO-GO",
        "source_only": True,
        "materially_distinct_from_experiment_117": True,
        "direct_query_pan_effect_added": False,
        "held_outcomes_accessed": False,
        "selected": asdict(selected) | {"key": selected.key},
        "source_reproducibility": {
            "authorization_sha256": sha256_file(AUTHORIZATION_PATH),
            "selection_freeze_sha256": sha256_file(SELECTION_DIR / "freeze-manifest.json"),
            "source_test_access_lock_sha256": sha256_file(TEST_ACCESS_LOCK),
        },
        "source_only_selection": selection["selection_audit"],
        "untouched_source_test": observed,
        "target_block_bootstrap": bootstrap,
        "target_label_permutation_null": nulls,
        "promotion_checks": checks,
        "promotion_passed": passed,
        "independent_confirmation_required": (
            "A pass authorizes only a separately reviewed held evaluation on a surface not used "
            "by Experiments 107, 117, or 118; it does not authorize production."
        ),
        "peak_rss_gib": peak_rss_gib(),
        "result_arrays_sha256": sha256_file(arrays_path),
    }
    result_path = RESULT_DIR / "results.json"
    write_json(result_path, result)
    manifest = {
        "schema_version": 1,
        "experiment": "129-gse295214-transfer",
        "status": result["status"],
        "files": {
            "results.json": sha256_file(result_path),
            "result-arrays.npz": sha256_file(arrays_path),
        },
        "authorization_sha256": sha256_file(AUTHORIZATION_PATH),
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "selection_freeze_sha256": sha256_file(SELECTION_DIR / "freeze-manifest.json"),
        "source_test_access_lock_sha256": sha256_file(TEST_ACCESS_LOCK),
        "implementation_sha256": sha256_file(Path(__file__)),
    }
    write_json(RESULT_DIR / "result-manifest.json", manifest)
    print(json.dumps(result, indent=2, sort_keys=True))


def verify_result() -> dict[str, Any]:
    selection, _, _ = verify_selection()
    require(selection["sealed_before_source_test_effect_access"] is True, "selection seal invalid")
    manifest = read_json(RESULT_DIR / "result-manifest.json")
    require(manifest["authorization_sha256"] == sha256_file(AUTHORIZATION_PATH), "authorization changed")
    require(manifest["protocol_sha256"] == sha256_file(PROTOCOL_PATH), "protocol changed")
    require(manifest["selection_freeze_sha256"] == sha256_file(SELECTION_DIR / "freeze-manifest.json"), "selection changed")
    require(manifest["source_test_access_lock_sha256"] == sha256_file(TEST_ACCESS_LOCK), "test lock changed")
    require(manifest["implementation_sha256"] == sha256_file(Path(__file__)), "implementation changed")
    for name, digest in manifest["files"].items():
        require(sha256_file(RESULT_DIR / name) == digest, f"result file changed: {name}")
    result = read_json(RESULT_DIR / "results.json")
    require(result["held_outcomes_accessed"] is False, "held-access statement changed")
    require(result["direct_query_pan_effect_added"] is False, "direct Pan boundary changed")
    return result


def audit_stage() -> None:
    authorization = validate_authorization(hash_large_inputs=False)
    print(
        json.dumps(
            {
                "status": "GO-for-source-only-execution",
                "candidate_count": len(candidates()),
                "plan": validate_plan(),
                "authorization_sha256": sha256_file(AUTHORIZATION_PATH),
                "implementation_sha256": sha256_file(Path(__file__)),
                "held_outcome_access_authorized": authorization["held_outcome_access_authorized"],
                "cell_eval2_authorized": authorization["cell_eval2_authorized"],
                "challenge_build_or_submission_authorized": authorization[
                    "challenge_build_or_submission_authorized"
                ],
            },
            indent=2,
            sort_keys=True,
        )
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "command",
        choices=("plan", "freeze", "audit", "select", "evaluate", "verify"),
    )
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=32.0)
    args = parser.parse_args()
    if args.command == "plan":
        create_plan()
    elif args.command == "freeze":
        freeze_authorization()
    elif args.command == "audit":
        audit_stage()
    elif args.command == "select":
        selection_stage(args)
    elif args.command == "evaluate":
        evaluate_stage(args)
    else:
        print(json.dumps(verify_result(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

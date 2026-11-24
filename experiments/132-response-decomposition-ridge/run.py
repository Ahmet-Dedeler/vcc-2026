#!/usr/bin/env python3
"""Frozen source-only evaluation for Experiment 132.

Stages:

* plan: metadata and DepMap features only; no perturbation-effect values.
* freeze: hash inputs and implementation before effect-value access.
* select: read only the 613 source-selection target rows.
* evaluate: after a nonidentity selection freeze, read the 410 untouched rows.

Every output is written below this experiment directory.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import resource
import subprocess
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
from sklearn.decomposition import PCA


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL = HERE / "preregistered_protocol.json"
RECOVERY = HERE / "IMPLEMENTATION_RECOVERY_001.json"
PLAN_DIR = HERE / "work" / "preregistered-freeze-v2"
AUTHORIZATION = HERE / "SOURCE_ONLY_AUTHORIZATION.v2.json"
SELECTION_DIR = HERE / "work" / "selection-v2"
SELECTION_FREEZE = HERE / "SOURCE_CANDIDATE_FREEZE.v2.json"
TEST_ACCESS_LOCK = HERE / "work" / "source-test-access.v2.lock"
RESULT_DIR = HERE / "work" / "source-only-result-v2"

PANEL_PATH = ROOT / "data/challenge/pert_counts.csv"
CANONICAL_GENES_PATH = ROOT / "data/challenge/gene_names.csv"
DEPMAP_PATH = (
    ROOT
    / "experiments/114-depmap-coessentiality/work/source-features/depmap-features.npz"
)
DEPMAP_RECEIPT = (
    ROOT
    / "experiments/114-depmap-coessentiality/work/source-features/feature-receipt.json"
)
UPSTREAM = ROOT / "data/external/perturbation-decomposition"
UPSTREAM_COMMIT = "a15214780619736d393f40240e56ba992fd416a3"

SOURCE_ROOTS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "k562_essential": ROOT / "data/donors/k562_essential",
    "rpe1": ROOT / "data/donors/rpe1",
    "nadig_jurkat": ROOT / "data/donors/nadig_jurkat",
    "pan_h1": ROOT / "experiments/115-gse295214-h1/work/compact-donor",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
SCORED_DATASETS = ("k562_gwps", "k562_essential", "rpe1", "nadig_jurkat")
HELD_GROUPS = ("k562", "rpe1", "jurkat")
CONTEXTS = ("k562", "rpe1", "jurkat", "pan_h1", "orion_hct116", "orion_hek293t")
DATASET_GROUP = {
    "k562_gwps": "k562",
    "k562_essential": "k562",
    "rpe1": "rpe1",
    "nadig_jurkat": "jurkat",
}

EXPECTED_TARGETS = 1023
SELECTION_COUNT = 613
TEST_COUNT = 410
EXPECTED_RESPONSE_GENES = 4992
TARGET_NAMESPACE = "experiment-132-response-decomposition-ridge-target-v1"
NULL_DRAWS = 49
NULL_SEED = 13_220_261
BOOTSTRAP_DRAWS = 20_000
BOOTSTRAP_SEED = 13_220_262
PCA_COMPONENTS = 50
RIDGE_ALPHA = 10.0
PCA_SEED = 13_220_250
MIN_MACRO_PDS_GAIN = 0.003
CLUSTER_K = 5
CLUSTER_THRESHOLD = 0.40
EPS = 1e-12
GIB = 1024**3


class ContractError(RuntimeError):
    """Raised when a frozen experiment boundary is violated."""


@dataclass(frozen=True)
class Bank:
    name: str
    root: Path
    perts: tuple[str, ...]
    genes: tuple[str, ...]

    @classmethod
    def from_root(cls, name: str, root: Path) -> "Bank":
        return cls(
            name=name,
            root=root,
            perts=read_lines(root / "perts.txt"),
            genes=read_lines(root / "genes.txt", require_unique=False),
        )

    @property
    def delta_path(self) -> Path:
        return self.root / "delta.npy"


@dataclass(frozen=True)
class RidgeProjection:
    """PCA scores and the collapsed exact Ridge linear operator."""

    centered_scores: np.ndarray
    solve_left: np.ndarray

    def predict(self, response_mean: np.ndarray, permutation: np.ndarray | None = None) -> np.ndarray:
        values = np.asarray(response_mean, dtype=np.float64)
        require(values.ndim == 2, "response_mean must be a matrix")
        n = values.shape[0]
        if permutation is None:
            ordered = values
        else:
            perm = np.asarray(permutation, dtype=np.int64)
            require(perm.shape == (n,), "permutation shape changed")
            require(np.array_equal(np.sort(perm), np.arange(n)), "invalid permutation")
            inverse = np.argsort(perm)
            ordered = values[inverse]
        mean = ordered.mean(axis=0)
        weights = self.solve_left @ (ordered - mean)
        prediction_at_feature_position = self.centered_scores @ weights + mean
        if permutation is None:
            return prediction_at_feature_position
        return prediction_at_feature_position[np.asarray(permutation, dtype=np.int64)]


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


def read_lines(path: Path, *, require_unique: bool = True) -> tuple[str, ...]:
    values = tuple(line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip())
    require(bool(values), f"{path}: empty axis")
    if require_unique:
        require(len(values) == len(set(values)), f"{path}: duplicate axis")
    return values


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
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


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def atomic_write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".tmp")
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
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.", suffix=".npz")
    os.close(descriptor)
    try:
        np.savez_compressed(temporary, **arrays)
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def write_lines(path: Path, values: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text("\n".join(values) + "\n", encoding="utf-8")
    temporary.replace(path)


def stable_order(values: Iterable[str]) -> tuple[str, ...]:
    prefix = TARGET_NAMESPACE.encode("utf-8") + b"\0"
    return tuple(
        sorted(
            values,
            key=lambda value: (
                hashlib.sha256(prefix + value.encode("utf-8")).digest(),
                value.encode("utf-8"),
            ),
        )
    )


def derangement(rng: np.random.Generator, size: int) -> np.ndarray:
    base = np.arange(size, dtype=np.int16)
    for _ in range(10_000):
        trial = rng.permutation(base)
        if np.all(trial != base):
            return trial.astype(np.int16)
    raise ContractError(f"failed to create derangement of {size} rows")


def row_unit(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    norms = np.linalg.norm(array, axis=1, keepdims=True)
    return np.divide(array, norms, out=np.zeros_like(array), where=norms > EPS)


def reciprocal_feature_blocks(features: np.ndarray) -> np.ndarray:
    unit = row_unit(features)
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
    labels = {root: index for index, root in enumerate(sorted(set(roots)))}
    return np.asarray([labels[root] for root in roots], dtype=np.int16)


def load_depmap() -> tuple[tuple[str, ...], np.ndarray]:
    with np.load(DEPMAP_PATH, allow_pickle=False) as data:
        require(set(data.files) == {"genes", "matrix"}, "DepMap artifact schema changed")
        genes = tuple(map(str, data["genes"]))
        matrix = np.asarray(data["matrix"], dtype=np.float32)
    require(matrix.shape == (17_916, 1_178), "DepMap matrix shape changed")
    require(len(genes) == len(set(genes)), "DepMap target axis has duplicates")
    require(np.all(np.isfinite(matrix)), "DepMap matrix contains nonfinite values")
    return genes, matrix


def feature_rows(targets: Sequence[str]) -> np.ndarray:
    genes, matrix = load_depmap()
    index = {gene: row for row, gene in enumerate(genes)}
    require(all(target in index for target in targets), "target missing from DepMap")
    values = matrix[np.asarray([index[target] for target in targets], dtype=np.int64)]
    require(np.all(np.isfinite(values)), "selected DepMap features are nonfinite")
    return values


def source_banks() -> dict[str, Bank]:
    return {name: Bank.from_root(name, path) for name, path in SOURCE_ROOTS.items()}


def create_plan() -> None:
    require(not PLAN_DIR.exists(), f"{PLAN_DIR} exists; refusing overwrite")
    require(not AUTHORIZATION.exists(), "authorization already exists")
    require(not SELECTION_DIR.exists() and not RESULT_DIR.exists(), "later output already exists")

    panel = set(read_csv_column(PANEL_PATH, "target_gene"))
    canonical = read_csv_column(CANONICAL_GENES_PATH, "gene_name")
    banks = source_banks()
    depmap_genes, _ = load_depmap()
    target_set = set(depmap_genes) - panel
    for bank in banks.values():
        target_set &= set(bank.perts)
    require(len(target_set) == EXPECTED_TARGETS, f"target count {len(target_set)} != {EXPECTED_TARGETS}")
    targets = stable_order(target_set)
    selection_targets = targets[:SELECTION_COUNT]
    test_targets = targets[SELECTION_COUNT:]
    require(len(test_targets) == TEST_COUNT, "test target count changed")

    common_genes = set(canonical)
    for bank in banks.values():
        common_genes &= set(bank.genes)
    pan_bank = banks["pan_h1"]
    pan_measured = np.load(pan_bank.root / "measured_mask.npy", mmap_mode="r")
    require(
        pan_measured.shape == (len(pan_bank.perts), len(pan_bank.genes)),
        "Pan measured-mask shape changed",
    )
    pan_all_measured = np.all(pan_measured, axis=0)
    pan_measured_genes = {
        gene for gene, measured in zip(pan_bank.genes, pan_all_measured, strict=True) if measured
    }
    common_genes &= pan_measured_genes
    response_genes = tuple(
        gene for gene in canonical if gene in common_genes and gene not in panel and gene not in target_set
    )
    require(
        len(response_genes) == EXPECTED_RESPONSE_GENES,
        f"response gene count {len(response_genes)} != {EXPECTED_RESPONSE_GENES}",
    )

    selection_features = feature_rows(selection_targets)
    test_features = feature_rows(test_targets)
    rng = np.random.default_rng(NULL_SEED)
    selection_permutations = np.stack(
        [derangement(rng, SELECTION_COUNT) for _ in range(NULL_DRAWS)]
    )
    test_permutations = np.stack([derangement(rng, TEST_COUNT) for _ in range(NULL_DRAWS)])
    blocks = reciprocal_feature_blocks(test_features)
    block_count = int(blocks.max()) + 1
    block_sizes = np.bincount(blocks)
    require(block_count < TEST_COUNT, "target-block construction produced only singletons")
    require(int(block_sizes.max()) <= 20, "target-block construction produced an oversized block")
    bootstrap_rng = np.random.default_rng(BOOTSTRAP_SEED)
    bootstrap_draws = bootstrap_rng.integers(
        0,
        block_count,
        size=(BOOTSTRAP_DRAWS, block_count),
        dtype=np.uint16,
    )

    PLAN_DIR.mkdir(parents=True)
    write_lines(PLAN_DIR / "all_targets.txt", targets)
    write_lines(PLAN_DIR / "selection_targets.txt", selection_targets)
    write_lines(PLAN_DIR / "test_targets.txt", test_targets)
    write_lines(PLAN_DIR / "response_genes.txt", response_genes)
    atomic_savez(
        PLAN_DIR / "frozen-randomness.npz",
        selection_permutations=selection_permutations,
        test_permutations=test_permutations,
        test_block_labels=blocks,
        bootstrap_block_draws=bootstrap_draws,
    )
    files = {
        path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
        for path in sorted(PLAN_DIR.iterdir())
        if path.name != "plan-manifest.json"
    }
    manifest = {
        "schema_version": 1,
        "created_before_effect_value_access": True,
        "target_count": len(targets),
        "selection_count": len(selection_targets),
        "test_count": len(test_targets),
        "response_gene_count": len(response_genes),
        "target_sha256": sha256_lines(targets),
        "selection_target_sha256": sha256_lines(selection_targets),
        "test_target_sha256": sha256_lines(test_targets),
        "response_gene_sha256": sha256_lines(response_genes),
        "test_block_count": block_count,
        "test_multi_target_block_count": int(np.sum(block_sizes > 1)),
        "test_targets_in_multi_target_blocks": int(block_sizes[block_sizes > 1].sum()),
        "test_max_block_size": int(block_sizes.max()),
        "candidate_grid": ["identity_source_mean", "depmap_pca50_ridge10"],
        "null_draws_per_stage": NULL_DRAWS,
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "files": files,
    }
    atomic_write_json(PLAN_DIR / "plan-manifest.json", manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))


def validate_plan() -> dict[str, Any]:
    require(PLAN_DIR.is_dir(), "plan directory missing")
    manifest = read_json(PLAN_DIR / "plan-manifest.json")
    require(manifest["target_count"] == EXPECTED_TARGETS, "target count drifted")
    require(manifest["response_gene_count"] == EXPECTED_RESPONSE_GENES, "response axis drifted")
    for name, expected in manifest["files"].items():
        path = PLAN_DIR / name
        require(path.is_file(), f"plan file missing: {name}")
        require(path.stat().st_size == expected["bytes"], f"plan byte count changed: {name}")
        require(sha256_file(path) == expected["sha256"], f"plan hash changed: {name}")
    return manifest


def input_paths() -> tuple[Path, ...]:
    paths = [
        PANEL_PATH,
        CANONICAL_GENES_PATH,
        DEPMAP_PATH,
        DEPMAP_RECEIPT,
        UPSTREAM / "models/ridge.py",
        UPSTREAM / "benchmarks/cross_cl.py",
        UPSTREAM / "data/process.py",
        UPSTREAM / "evaluation/metrics.py",
        ROOT / "experiments/115-gse295214-h1/work/compact-donor/artifact-manifest.json",
        ROOT / "experiments/115-gse295214-h1/work/compact-donor/measured_mask.npy",
    ]
    for root in SOURCE_ROOTS.values():
        paths.extend((root / "perts.txt", root / "genes.txt", root / "delta.npy"))
    return tuple(paths)


def implementation_paths() -> tuple[Path, ...]:
    return (
        HERE / "DISTINCTNESS.md",
        HERE / "README.md",
        PROTOCOL,
        RECOVERY,
        HERE / "run.py",
        HERE / "test_run.py",
    )


def freeze() -> None:
    require(not AUTHORIZATION.exists(), "authorization already exists")
    require(not SELECTION_DIR.exists() and not RESULT_DIR.exists(), "effect-stage output exists")
    plan = validate_plan()
    commit = subprocess.check_output(
        ["git", "-C", str(UPSTREAM), "rev-parse", "HEAD"], text=True
    ).strip()
    require(commit == UPSTREAM_COMMIT, f"upstream commit changed: {commit}")
    bindings = {}
    for path in input_paths():
        require(path.is_file(), f"missing input: {path}")
        relative = str(path.relative_to(ROOT))
        bindings[relative] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
        print(f"[hash] {relative}", flush=True)
    implementation = {
        str(path.relative_to(ROOT)): sha256_file(path) for path in implementation_paths()
    }
    authorization = {
        "schema_version": 1,
        "experiment": "132-response-decomposition-ridge",
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "created_before_experiment_132_effect_value_loading": True,
        "source_only_authorized": True,
        "held_outcome_access_authorized": False,
        "cell_eval2_authorized": False,
        "count_generation_authorized": False,
        "challenge_build_or_submission_authorized": False,
        "commit_or_push_authorized": False,
        "upstream_commit": commit,
        "recovery_sha256": sha256_file(RECOVERY),
        "supersedes_authorization_sha256": (
            "f3f05fc88325b1d5643712c49ef0e1308a2de33f1b069a3551eb0e744e68038f"
        ),
        "protocol_sha256": sha256_file(PROTOCOL),
        "plan_manifest_sha256": sha256_file(PLAN_DIR / "plan-manifest.json"),
        "implementation": implementation,
        "inputs": bindings,
        "plan": {
            "target_count": plan["target_count"],
            "selection_count": plan["selection_count"],
            "test_count": plan["test_count"],
            "response_gene_count": plan["response_gene_count"],
            "candidate_grid": plan["candidate_grid"],
        },
    }
    atomic_write_json(AUTHORIZATION, authorization)
    print(json.dumps(authorization, indent=2, sort_keys=True))


def validate_authorization(*, check_large_hashes: bool = False) -> dict[str, Any]:
    require(AUTHORIZATION.is_file(), "source-only authorization missing")
    value = read_json(AUTHORIZATION)
    require(value["source_only_authorized"] is True, "source-only access disabled")
    require(value["held_outcome_access_authorized"] is False, "held access unexpectedly enabled")
    require(value["cell_eval2_authorized"] is False, "cell-eval2 unexpectedly enabled")
    require(value["challenge_build_or_submission_authorized"] is False, "challenge build enabled")
    require(value["protocol_sha256"] == sha256_file(PROTOCOL), "protocol changed after freeze")
    require(value["recovery_sha256"] == sha256_file(RECOVERY), "recovery changed after freeze")
    require(
        value["plan_manifest_sha256"] == sha256_file(PLAN_DIR / "plan-manifest.json"),
        "plan changed after freeze",
    )
    for relative, expected in value["implementation"].items():
        require(sha256_file(ROOT / relative) == expected, f"implementation changed: {relative}")
    for relative, expected in value["inputs"].items():
        path = ROOT / relative
        require(path.is_file() and path.stat().st_size == expected["bytes"], f"input changed: {relative}")
        if check_large_hashes:
            require(sha256_file(path) == expected["sha256"], f"input hash changed: {relative}")
    validate_plan()
    return value


def axis_indices(bank: Bank, targets: Sequence[str], genes: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    pert_index = {target: index for index, target in enumerate(bank.perts)}
    gene_index = {gene: index for index, gene in enumerate(bank.genes)}
    require(all(target in pert_index for target in targets), f"{bank.name}: missing target")
    require(all(gene in gene_index for gene in genes), f"{bank.name}: missing response gene")
    return (
        np.asarray([pert_index[target] for target in targets], dtype=np.int64),
        np.asarray([gene_index[gene] for gene in genes], dtype=np.int64),
    )


def load_effect_rows(bank: Bank, targets: Sequence[str], genes: Sequence[str]) -> np.ndarray:
    rows, columns = axis_indices(bank, targets, genes)
    if bank.name == "pan_h1":
        measured = np.load(bank.root / "measured_mask.npy", mmap_mode="r")
        support = np.asarray(measured[np.ix_(rows, columns)], dtype=bool)
        require(np.all(support), "Pan contains an unmeasured selected coordinate")
    matrix = np.load(bank.delta_path, mmap_mode="r")
    require(matrix.shape == (len(bank.perts), len(bank.genes)), f"{bank.name}: matrix shape changed")
    values = np.asarray(matrix[np.ix_(rows, columns)], dtype=np.float32)
    require(values.shape == (len(targets), len(genes)), f"{bank.name}: sliced shape changed")
    require(np.all(np.isfinite(values)), f"{bank.name}: nonfinite effect values")
    return values


def load_split_effects(targets: Sequence[str], response_genes: Sequence[str], limit_gib: float) -> dict[str, np.ndarray]:
    values: dict[str, np.ndarray] = {}
    for name, bank in source_banks().items():
        values[name] = load_effect_rows(bank, targets, response_genes)
        memory_checkpoint(f"loaded {name}", limit_gib)
    return values


def context_effects(source: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    return {
        "k562": 0.5 * (
            np.asarray(source["k562_gwps"], dtype=np.float64)
            + np.asarray(source["k562_essential"], dtype=np.float64)
        ),
        "rpe1": np.asarray(source["rpe1"], dtype=np.float64),
        "jurkat": np.asarray(source["nadig_jurkat"], dtype=np.float64),
        "pan_h1": np.asarray(source["pan_h1"], dtype=np.float64),
        "orion_hct116": np.asarray(source["orion_hct116"], dtype=np.float64),
        "orion_hek293t": np.asarray(source["orion_hek293t"], dtype=np.float64),
    }


def source_means_by_held(context: Mapping[str, np.ndarray]) -> dict[str, np.ndarray]:
    result = {}
    for held in HELD_GROUPS:
        allowed = [context[name] for name in CONTEXTS if name != held]
        require(len(allowed) == 5, "source-context count changed")
        result[held] = np.mean(np.stack(allowed, axis=0), axis=0)
    return result


def fit_projection(
    features: np.ndarray,
    *,
    n_sources: int = 5,
    n_components: int = PCA_COMPONENTS,
    alpha: float = RIDGE_ALPHA,
    seed: int = PCA_SEED,
) -> RidgeProjection:
    features = np.asarray(features, dtype=np.float64)
    n_comp = min(n_components, features.shape[0] - 1, features.shape[1])
    require(n_comp > 0, "PCA component count is zero")
    repeated = np.tile(features, (n_sources, 1))
    state = np.random.get_state()
    np.random.seed(seed)
    try:
        pca = PCA(n_comp).fit(repeated)
    finally:
        np.random.set_state(state)
    scores = np.asarray(pca.transform(features), dtype=np.float64)
    centered = scores - scores.mean(axis=0)
    gram = centered.T @ centered
    solve_left = np.linalg.solve(
        gram + (alpha / n_sources) * np.eye(n_comp, dtype=np.float64),
        centered.T,
    )
    return RidgeProjection(centered_scores=centered, solve_left=solve_left)


def ridge_predictions(
    projection: RidgeProjection,
    response_means: Mapping[str, np.ndarray],
    permutation: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    keys = tuple(HELD_GROUPS)
    width = next(iter(response_means.values())).shape[1]
    joined = np.concatenate([response_means[key] for key in keys], axis=1)
    predicted = projection.predict(joined, permutation=permutation)
    return {
        key: predicted[:, index * width : (index + 1) * width]
        for index, key in enumerate(keys)
    }


def cosine_matrix(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return row_unit(left) @ row_unit(right).T


def pds_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    require(prediction.shape == truth.shape, "PDS shape mismatch")
    count = prediction.shape[0]
    require(count > 1, "PDS gallery must have at least two targets")
    similarity = cosine_matrix(prediction, truth)
    own = np.diag(similarity)
    greater = np.sum(similarity > own[:, None], axis=1)
    equal = np.sum(similarity == own[:, None], axis=1)
    midrank = 1.0 + greater + 0.5 * (equal - 1)
    return 1.0 - (midrank - 1.0) / (count - 1.0)


def signed_credit_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    require(prediction.shape == truth.shape, "signed-credit shape mismatch")
    numerator = np.sum(np.sign(truth) * prediction, axis=1)
    denominator = np.sum(np.abs(prediction), axis=1)
    return np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator, dtype=np.float64),
        where=denominator > EPS,
    )


def summarize(
    candidate_by_group: Mapping[str, np.ndarray],
    identity_by_group: Mapping[str, np.ndarray],
    truth: Mapping[str, np.ndarray],
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    contexts: dict[str, Any] = {}
    arrays: dict[str, np.ndarray] = {}
    for dataset in SCORED_DATASETS:
        group = DATASET_GROUP[dataset]
        candidate_pds = pds_per_target(candidate_by_group[group], truth[dataset])
        identity_pds = pds_per_target(identity_by_group[group], truth[dataset])
        candidate_credit = signed_credit_per_target(candidate_by_group[group], truth[dataset])
        identity_credit = signed_credit_per_target(identity_by_group[group], truth[dataset])
        arrays[f"{dataset}__pds_delta"] = candidate_pds - identity_pds
        arrays[f"{dataset}__signed_credit_delta"] = candidate_credit - identity_credit
        contexts[dataset] = {
            "identity": {
                "pds": float(identity_pds.mean()),
                "signed_credit": float(identity_credit.mean()),
            },
            "candidate": {
                "pds": float(candidate_pds.mean()),
                "signed_credit": float(candidate_credit.mean()),
            },
            "delta": {
                "pds": float((candidate_pds - identity_pds).mean()),
                "signed_credit": float((candidate_credit - identity_credit).mean()),
            },
        }

    family_arrays = {}
    for metric in ("pds", "signed_credit"):
        k562 = 0.5 * (
            arrays[f"k562_gwps__{metric}_delta"]
            + arrays[f"k562_essential__{metric}_delta"]
        )
        family_arrays[metric] = np.mean(
            np.stack(
                [
                    k562,
                    arrays[f"rpe1__{metric}_delta"],
                    arrays[f"nadig_jurkat__{metric}_delta"],
                ],
                axis=0,
            ),
            axis=0,
        )
        arrays[f"macro_family__{metric}_delta"] = family_arrays[metric]

    summary = {
        "contexts": contexts,
        "macro_delta": {
            metric: float(family_arrays[metric].mean())
            for metric in ("pds", "signed_credit")
        },
    }
    return summary, arrays


def null_summary(
    projection: RidgeProjection,
    response_means: Mapping[str, np.ndarray],
    identity: Mapping[str, np.ndarray],
    truth: Mapping[str, np.ndarray],
    permutations: np.ndarray,
    limit_gib: float,
) -> dict[str, np.ndarray]:
    pds = np.empty(len(permutations), dtype=np.float64)
    credit = np.empty(len(permutations), dtype=np.float64)
    for draw, permutation in enumerate(permutations):
        candidate = ridge_predictions(projection, response_means, permutation)
        summary, _ = summarize(candidate, identity, truth)
        pds[draw] = summary["macro_delta"]["pds"]
        credit[draw] = summary["macro_delta"]["signed_credit"]
        if draw % 10 == 9:
            memory_checkpoint(f"null draw {draw + 1}/{len(permutations)}", limit_gib)
    return {"pds": pds, "signed_credit": credit}


def selection_passes(summary: Mapping[str, Any], nulls: Mapping[str, np.ndarray]) -> dict[str, bool]:
    contexts = summary["contexts"]
    checks = {
        "pds_nonnegative_every_dataset": all(
            contexts[name]["delta"]["pds"] >= 0.0 for name in SCORED_DATASETS
        ),
        "signed_credit_nonnegative_every_dataset": all(
            contexts[name]["delta"]["signed_credit"] >= 0.0 for name in SCORED_DATASETS
        ),
        "macro_pds_at_least_0.003": summary["macro_delta"]["pds"] >= MIN_MACRO_PDS_GAIN,
        "macro_signed_credit_positive": summary["macro_delta"]["signed_credit"] > 0.0,
        "macro_pds_above_null_p95": summary["macro_delta"]["pds"]
        > float(np.quantile(nulls["pds"], 0.95, method="linear")),
        "macro_signed_credit_above_null_p95": summary["macro_delta"]["signed_credit"]
        > float(np.quantile(nulls["signed_credit"], 0.95, method="linear")),
    }
    return checks


def load_plan_arrays() -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...], dict[str, np.ndarray]]:
    validate_plan()
    selection = read_lines(PLAN_DIR / "selection_targets.txt")
    test = read_lines(PLAN_DIR / "test_targets.txt")
    genes = read_lines(PLAN_DIR / "response_genes.txt")
    with np.load(PLAN_DIR / "frozen-randomness.npz", allow_pickle=False) as data:
        randomness = {key: np.asarray(data[key]) for key in data.files}
    return selection, test, genes, randomness


def output_manifest(directory: Path, *, exclude: set[str] | None = None) -> dict[str, Any]:
    excluded = exclude or set()
    files = {}
    for path in sorted(directory.iterdir()):
        if path.is_file() and path.name not in excluded:
            files[path.name] = {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
    return {"files": files}


def run_selection(limit_gib: float) -> None:
    validate_authorization()
    require(not SELECTION_DIR.exists(), "selection output exists; refusing overwrite")
    require(not SELECTION_FREEZE.exists(), "selection freeze exists")
    require(not TEST_ACCESS_LOCK.exists() and not RESULT_DIR.exists(), "test stage already exists")
    selection_targets, _, response_genes, randomness = load_plan_arrays()
    SELECTION_DIR.mkdir(parents=True)
    atomic_write_json(
        SELECTION_DIR / "effect-access-receipt.json",
        {
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "stage": "selection",
            "targets": len(selection_targets),
            "test_effect_values_accessed": False,
        },
    )
    effects = load_split_effects(selection_targets, response_genes, limit_gib)
    contexts = context_effects(effects)
    identity = source_means_by_held(contexts)
    features = feature_rows(selection_targets)
    projection = fit_projection(features)
    candidate = ridge_predictions(projection, identity)
    observed, per_target = summarize(candidate, identity, effects)
    nulls = null_summary(
        projection,
        identity,
        identity,
        effects,
        randomness["selection_permutations"],
        limit_gib,
    )
    checks = selection_passes(observed, nulls)
    selected = "depmap_pca50_ridge10" if all(checks.values()) else "identity_source_mean"
    result = {
        "schema_version": 1,
        "stage": "source_selection",
        "candidate": "depmap_pca50_ridge10",
        "selected": selected,
        "observed": observed,
        "null": {
            metric: {
                "draws": len(values),
                "p95": float(np.quantile(values, 0.95, method="linear")),
                "add_one_p": float((1 + np.sum(values >= observed["macro_delta"][metric])) / (len(values) + 1)),
            }
            for metric, values in nulls.items()
        },
        "checks": checks,
        "passes": selected == "depmap_pca50_ridge10",
        "peak_rss_gib": peak_rss_gib(),
    }
    atomic_write_json(SELECTION_DIR / "result.json", result)
    atomic_savez(
        SELECTION_DIR / "per-target-and-nulls.npz",
        **per_target,
        null__pds=nulls["pds"],
        null__signed_credit=nulls["signed_credit"],
    )
    manifest = output_manifest(SELECTION_DIR, exclude={"artifact-manifest.json"})
    atomic_write_json(SELECTION_DIR / "artifact-manifest.json", manifest)
    freeze_value = {
        "schema_version": 1,
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "selected_candidate": selected,
        "selection_passes": result["passes"],
        "selection_result_sha256": sha256_file(SELECTION_DIR / "result.json"),
        "selection_artifact_manifest_sha256": sha256_file(SELECTION_DIR / "artifact-manifest.json"),
        "authorization_sha256": sha256_file(AUTHORIZATION),
        "test_effect_access_authorized": result["passes"],
        "production_authorized": False,
    }
    atomic_write_json(SELECTION_FREEZE, freeze_value)
    print(json.dumps(result, indent=2, sort_keys=True))


def validate_selection_freeze() -> dict[str, Any]:
    require(SELECTION_FREEZE.is_file(), "selection freeze missing")
    value = read_json(SELECTION_FREEZE)
    require(
        value["authorization_sha256"] == sha256_file(AUTHORIZATION),
        "authorization changed after selection",
    )
    require(
        value["selection_result_sha256"] == sha256_file(SELECTION_DIR / "result.json"),
        "selection result changed",
    )
    require(
        value["selection_artifact_manifest_sha256"]
        == sha256_file(SELECTION_DIR / "artifact-manifest.json"),
        "selection manifest changed",
    )
    return value


def block_bootstrap(values: np.ndarray, labels: np.ndarray, draws: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    labels = np.asarray(labels, dtype=np.int64)
    count = int(labels.max()) + 1
    sums = np.bincount(labels, weights=values, minlength=count)
    sizes = np.bincount(labels, minlength=count).astype(np.float64)
    numerator = sums[draws].sum(axis=1)
    denominator = sizes[draws].sum(axis=1)
    return numerator / denominator


def run_evaluation(limit_gib: float) -> None:
    validate_authorization()
    selection = validate_selection_freeze()
    require(selection["selected_candidate"] == "depmap_pca50_ridge10", "selection chose identity")
    require(selection["test_effect_access_authorized"] is True, "test effect access denied")
    require(not TEST_ACCESS_LOCK.exists(), "test access lock exists")
    require(not RESULT_DIR.exists(), "result output exists; refusing overwrite")
    _, test_targets, response_genes, randomness = load_plan_arrays()
    TEST_ACCESS_LOCK.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_json(
        TEST_ACCESS_LOCK,
        {
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "selection_freeze_sha256": sha256_file(SELECTION_FREEZE),
            "targets": len(test_targets),
            "irreversible_test_access": True,
        },
    )
    RESULT_DIR.mkdir(parents=True)
    effects = load_split_effects(test_targets, response_genes, limit_gib)
    contexts = context_effects(effects)
    identity = source_means_by_held(contexts)
    features = feature_rows(test_targets)
    projection = fit_projection(features)
    candidate = ridge_predictions(projection, identity)
    observed, per_target = summarize(candidate, identity, effects)
    nulls = null_summary(
        projection,
        identity,
        identity,
        effects,
        randomness["test_permutations"],
        limit_gib,
    )
    blocks = randomness["test_block_labels"]
    draws = randomness["bootstrap_block_draws"].astype(np.int64)
    bootstrap = {
        metric: block_bootstrap(per_target[f"macro_family__{metric}_delta"], blocks, draws)
        for metric in ("pds", "signed_credit")
    }
    contexts_result = observed["contexts"]
    checks = {
        "selection_passed": selection["selection_passes"] is True,
        "test_pds_nonnegative_every_dataset": all(
            contexts_result[name]["delta"]["pds"] >= 0.0 for name in SCORED_DATASETS
        ),
        "test_signed_credit_nonnegative_every_dataset": all(
            contexts_result[name]["delta"]["signed_credit"] >= 0.0 for name in SCORED_DATASETS
        ),
        "test_macro_pds_at_least_0.003": observed["macro_delta"]["pds"] >= MIN_MACRO_PDS_GAIN,
        "test_macro_signed_credit_positive": observed["macro_delta"]["signed_credit"] > 0.0,
        "pds_bootstrap_lower_positive": float(
            np.quantile(bootstrap["pds"], 0.025, method="linear")
        ) > 0.0,
        "signed_credit_bootstrap_lower_positive": float(
            np.quantile(bootstrap["signed_credit"], 0.025, method="linear")
        ) > 0.0,
        "test_macro_pds_above_null_p95": observed["macro_delta"]["pds"]
        > float(np.quantile(nulls["pds"], 0.95, method="linear")),
        "test_macro_signed_credit_above_null_p95": observed["macro_delta"]["signed_credit"]
        > float(np.quantile(nulls["signed_credit"], 0.95, method="linear")),
    }
    result = {
        "schema_version": 1,
        "stage": "untouched_source_test",
        "candidate": selection["selected_candidate"],
        "observed": observed,
        "null": {
            metric: {
                "draws": len(values),
                "p95": float(np.quantile(values, 0.95, method="linear")),
                "add_one_p": float((1 + np.sum(values >= observed["macro_delta"][metric])) / (len(values) + 1)),
            }
            for metric, values in nulls.items()
        },
        "bootstrap": {
            metric: {
                "draws": len(values),
                "lower_95": float(np.quantile(values, 0.025, method="linear")),
                "upper_95": float(np.quantile(values, 0.975, method="linear")),
            }
            for metric, values in bootstrap.items()
        },
        "checks": checks,
        "passes": all(checks.values()),
        "decision": "PASS_FOR_SEPARATE_REVIEW_ONLY" if all(checks.values()) else "NO-GO",
        "production_authorized": False,
        "peak_rss_gib": peak_rss_gib(),
    }
    atomic_write_json(RESULT_DIR / "result.json", result)
    atomic_savez(
        RESULT_DIR / "per-target-nulls-bootstrap.npz",
        **per_target,
        null__pds=nulls["pds"],
        null__signed_credit=nulls["signed_credit"],
        bootstrap__pds=bootstrap["pds"],
        bootstrap__signed_credit=bootstrap["signed_credit"],
    )
    manifest = output_manifest(RESULT_DIR, exclude={"artifact-manifest.json"})
    atomic_write_json(RESULT_DIR / "artifact-manifest.json", manifest)
    print(json.dumps(result, indent=2, sort_keys=True))


def audit() -> None:
    authorization = validate_authorization()
    payload: dict[str, Any] = {
        "authorization_sha256": sha256_file(AUTHORIZATION),
        "upstream_commit": authorization["upstream_commit"],
        "selection_exists": SELECTION_DIR.exists(),
        "selection_freeze_exists": SELECTION_FREEZE.exists(),
        "test_accessed": TEST_ACCESS_LOCK.exists(),
        "result_exists": RESULT_DIR.exists(),
    }
    if SELECTION_FREEZE.exists():
        payload["selection"] = validate_selection_freeze()
    if RESULT_DIR.exists():
        manifest = read_json(RESULT_DIR / "artifact-manifest.json")
        for name, expected in manifest["files"].items():
            path = RESULT_DIR / name
            require(path.stat().st_size == expected["bytes"], f"result bytes changed: {name}")
            require(sha256_file(path) == expected["sha256"], f"result hash changed: {name}")
        payload["result"] = read_json(RESULT_DIR / "result.json")
    print(json.dumps(payload, indent=2, sort_keys=True))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("plan", "freeze", "select", "evaluate", "audit"))
    parser.add_argument("--max-rss-gib", type=float, default=18.0)
    args = parser.parse_args()
    if args.command == "plan":
        create_plan()
    elif args.command == "freeze":
        freeze()
    elif args.command == "select":
        run_selection(args.max_rss_gib)
    elif args.command == "evaluate":
        run_evaluation(args.max_rss_gib)
    else:
        audit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

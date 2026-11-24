#!/usr/bin/env python3
"""Preregistered exact-SLIM source validation for Experiment 123.

The staged boundary is deliberate:

* ``plan`` reads labels and the public STRING embedding only.
* ``select`` may read fit and tune source rows, then seals candidates.
* ``evaluate`` may read test source rows only after that seal exists.

No command imports or invokes cell-eval2, count generation, challenge builders,
or submission code.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
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
UPSTREAM_PATH = HERE / "UPSTREAM_SOURCE_AUDIT.json"
AUTHORIZATION_PATH = HERE / "SOURCE_ONLY_AUTHORIZATION.json"
PLAN_DIR = HERE / "work" / "preregistered-freeze-v1"
SELECTION_DIR = HERE / "work" / "source-candidate-freeze-v1"
RESULT_DIR = HERE / "work" / "source-only-evaluation-v1"
TEST_ACCESS_LOCK = HERE / "work" / "source-test-access.lock"

DONORS = ("k562_gwps", "rpe1", "nadig_jurkat")
HELD_FORBIDDEN = ("h1_2025", "orion_hct116", "orion_hek293t")
EMBEDDING_PATH = ROOT / "data/external/slim/gene_string_embeddings.v0.3.h5"
PANEL_PATH = ROOT / "data/challenge/pert_counts.csv"
TARGET_NAMESPACE = "experiment-123-split-v1"
TARGET_COUNT = 600
FIT_COUNT = 360
TUNE_COUNT = 120
TEST_COUNT = 120
MAX_K = 20
K_GRID = (5, 10, 20)
RIDGE_GRID = (0.01, 0.1, 1.0, 10.0)
BLEND_GRID = (0.05, 0.10, 0.20, 0.35, 0.50, 1.0)
PAIRING_NULL_DRAWS = 32
GAUSSIAN_NULL_DRAWS = 32
PAIRING_SEED = 12320261
GAUSSIAN_SEED = 12320262
BOOTSTRAP_SEED = 12320263
BOOTSTRAP_DRAWS = 20_000
NORM_BINS = 8
CLUSTER_K = 5
CLUSTER_THRESHOLD = 0.80
EPS = 1e-12


class ContractError(RuntimeError):
    """Raised when a frozen experiment boundary is violated."""


@dataclass(frozen=True, order=True)
class Candidate:
    kind: str
    K: int
    ridge: float
    blend: float

    @property
    def key(self) -> str:
        if self.kind == "identity":
            return "identity"
        return f"slim-blend:K={self.K}:ridge={self.ridge:g}:f={self.blend:g}"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode("utf-8")).hexdigest()


def sha256_array(values: np.ndarray) -> str:
    array = np.ascontiguousarray(values)
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(b"\0")
    digest.update(json.dumps(list(array.shape), separators=(",", ":")).encode("ascii"))
    digest.update(b"\0")
    digest.update(array.view(np.uint8))
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def write_lines(path: Path, values: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(values) + "\n")


def read_lines(path: Path) -> tuple[str, ...]:
    values = tuple(line.strip() for line in path.read_text().splitlines() if line.strip())
    require(bool(values), f"{path}: empty axis")
    require(len(values) == len(set(values)), f"{path}: duplicate labels")
    return values


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
    result = [Candidate("identity", 0, 0.0, 0.0)]
    result.extend(
        Candidate("blend", K, ridge, blend)
        for K in K_GRID
        for ridge in RIDGE_GRID
        for blend in BLEND_GRID
    )
    require(len(result) == 73, "candidate grid drift")
    return tuple(result)


def load_embedding_rows(targets: Sequence[str]) -> np.ndarray:
    values = np.empty((len(targets), 64), dtype=np.float64)
    with h5py.File(EMBEDDING_PATH, "r") as handle:
        require(len(handle.keys()) == 21_688, "STRING root-row count changed")
        for index, target in enumerate(targets):
            require(target in handle, f"missing frozen STRING target {target}")
            row = np.asarray(handle[target][...], dtype=np.float64)
            require(row.shape == (64,), f"{target}: STRING dimension changed")
            require(np.all(np.isfinite(row)), f"{target}: non-finite STRING row")
            values[index] = row
    return values


def stable_target_order(targets: Iterable[str]) -> tuple[str, ...]:
    return tuple(
        sorted(
            targets,
            key=lambda target: (
                hashlib.sha256(
                    TARGET_NAMESPACE.encode("utf-8") + b"\0" + target.encode("utf-8")
                ).digest(),
                target.encode("utf-8"),
            ),
        )
    )


def quantile_bins(values: np.ndarray, count: int) -> np.ndarray:
    order = np.argsort(values, kind="stable")
    labels = np.empty(len(values), dtype=np.int16)
    for label, rows in enumerate(np.array_split(order, count)):
        labels[rows] = label
    return labels


def derangement(rng: np.random.Generator, rows: np.ndarray) -> np.ndarray:
    rows = np.asarray(rows, dtype=np.int64)
    if len(rows) < 2:
        return rows.copy()
    for _ in range(10_000):
        trial = rng.permutation(rows)
        if np.all(trial != rows):
            return trial
    raise ContractError("failed to create deterministic derangement")


def reciprocal_knn_clusters(embeddings: np.ndarray) -> np.ndarray:
    unit = row_unit(embeddings)
    similarity = unit @ unit.T
    np.fill_diagonal(similarity, -np.inf)
    neighbors = np.argsort(-similarity, axis=1, kind="stable")[:, :CLUSTER_K]
    parent = np.arange(len(embeddings), dtype=np.int64)

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
    for left in range(len(embeddings)):
        for right in neighbors[left]:
            right_int = int(right)
            if (
                left in neighbor_sets[right_int]
                and similarity[left, right_int] >= CLUSTER_THRESHOLD
            ):
                union(left, right_int)
    roots = [find(index) for index in range(len(embeddings))]
    mapping = {root: label for label, root in enumerate(sorted(set(roots)))}
    return np.asarray([mapping[root] for root in roots], dtype=np.int16)


def create_plan() -> None:
    require(not PLAN_DIR.exists(), f"{PLAN_DIR} exists; refusing overwrite")
    require(not SELECTION_DIR.exists(), "selection already exists")
    require(not RESULT_DIR.exists(), "result already exists")
    panel = set(read_csv_column(PANEL_PATH, "target_gene"))
    donor_targets = []
    donor_genes = []
    for donor in DONORS:
        directory = ROOT / "data/donors" / donor
        donor_targets.append(set(read_lines(directory / "perts.txt")))
        donor_genes.append(set(read_lines(directory / "genes.txt")))
    with h5py.File(EMBEDDING_PATH, "r") as handle:
        embedding_targets = set(map(str, handle.keys()))
    eligible = set.intersection(*donor_targets) - panel
    eligible &= embedding_targets
    require(len(eligible) >= 2000, "source target pool unexpectedly small")
    all_targets = stable_target_order(eligible)[:TARGET_COUNT]
    fit_targets = all_targets[:FIT_COUNT]
    tune_targets = all_targets[FIT_COUNT : FIT_COUNT + TUNE_COUNT]
    test_targets = all_targets[FIT_COUNT + TUNE_COUNT :]
    require(len(test_targets) == TEST_COUNT, "target split size changed")

    score_genes = tuple(
        sorted(
            set.intersection(*donor_genes) - panel - set(all_targets),
            key=lambda value: value.encode("utf-8"),
        )
    )
    protocol = read_json(PROTOCOL_PATH)
    expected_axis = protocol["response_axis"]
    require(len(score_genes) == int(expected_axis["expected_gene_count"]), "response-axis size changed")
    require(
        sha256_lines(score_genes) == expected_axis["expected_ordered_sha256"],
        "response-axis order changed",
    )

    all_embeddings = load_embedding_rows(all_targets)
    fit_embeddings = all_embeddings[:FIT_COUNT]
    norm_bins = quantile_bins(np.linalg.norm(fit_embeddings, axis=1), NORM_BINS)
    rng = np.random.default_rng(PAIRING_SEED)
    permutations = np.tile(np.arange(FIT_COUNT, dtype=np.int16), (PAIRING_NULL_DRAWS, 1))
    for draw in range(PAIRING_NULL_DRAWS):
        for label in range(NORM_BINS):
            rows = np.flatnonzero(norm_bins == label)
            permutations[draw, rows] = derangement(rng, rows)
        require(np.all(permutations[draw] != np.arange(FIT_COUNT)), "pairing null has a fixed point")

    gaussian_rng = np.random.default_rng(GAUSSIAN_SEED)
    gaussian_seeds = gaussian_rng.integers(
        0, np.iinfo(np.uint64).max, size=GAUSSIAN_NULL_DRAWS, dtype=np.uint64
    )

    test_embeddings = all_embeddings[-TEST_COUNT:]
    cluster_labels = reciprocal_knn_clusters(test_embeddings)
    cluster_count = int(cluster_labels.max()) + 1
    bootstrap_rng = np.random.default_rng(BOOTSTRAP_SEED)
    bootstrap_clusters = bootstrap_rng.integers(
        0,
        cluster_count,
        size=(BOOTSTRAP_DRAWS, cluster_count),
        dtype=np.uint16,
    )

    PLAN_DIR.mkdir(parents=True)
    write_lines(PLAN_DIR / "all_targets.txt", all_targets)
    write_lines(PLAN_DIR / "fit_targets.txt", fit_targets)
    write_lines(PLAN_DIR / "tune_targets.txt", tune_targets)
    write_lines(PLAN_DIR / "test_targets.txt", test_targets)
    write_lines(PLAN_DIR / "score_genes.txt", score_genes)
    np.save(PLAN_DIR / "target_pairing_permutations.npy", permutations, allow_pickle=False)
    np.save(PLAN_DIR / "gaussian_seeds.npy", gaussian_seeds, allow_pickle=False)
    np.save(PLAN_DIR / "test_target_cluster_labels.npy", cluster_labels, allow_pickle=False)
    np.save(PLAN_DIR / "bootstrap_cluster_draws.npy", bootstrap_clusters, allow_pickle=False)
    np.save(PLAN_DIR / "all_target_embeddings.npy", all_embeddings, allow_pickle=False)
    files = {}
    for path in sorted(PLAN_DIR.iterdir()):
        if path.name == "plan-manifest.json":
            continue
        files[path.name] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    manifest = {
        "schema_version": 1,
        "experiment": "123-slim-string-reproduction",
        "created_before_source_effect_access": True,
        "target_pool_count": len(eligible),
        "all_target_count": len(all_targets),
        "fit_count": len(fit_targets),
        "tune_count": len(tune_targets),
        "test_count": len(test_targets),
        "score_gene_count": len(score_genes),
        "test_cluster_count": cluster_count,
        "candidate_count": len(candidates()),
        "target_order_sha256": sha256_lines(all_targets),
        "fit_target_order_sha256": sha256_lines(fit_targets),
        "tune_target_order_sha256": sha256_lines(tune_targets),
        "test_target_order_sha256": sha256_lines(test_targets),
        "score_gene_order_sha256": sha256_lines(score_genes),
        "files": files,
    }
    write_json(PLAN_DIR / "plan-manifest.json", manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))


def inventory_paths() -> dict[str, str]:
    inventory = read_json(INVENTORY_PATH)
    result: dict[str, str] = {}
    for donor, members in inventory["sources"].items():
        for filename, digest in members.items():
            result[f"data/donors/{donor}/{filename}"] = digest
    result.update(inventory["embedding"])
    result.update(inventory["canonical_axes"])
    return result


def validate_plan() -> dict[str, Any]:
    require(PLAN_DIR.is_dir(), "preregistered plan is absent")
    manifest_path = PLAN_DIR / "plan-manifest.json"
    manifest = read_json(manifest_path)
    require(manifest["candidate_count"] == len(candidates()), "candidate count changed")
    for name, expected in manifest["files"].items():
        path = PLAN_DIR / name
        require(path.is_file(), f"missing planned file {name}")
        require(path.stat().st_size == expected["bytes"], f"{name}: byte count changed")
        require(sha256_file(path) == expected["sha256"], f"{name}: hash changed")
    return manifest


def validate_authorization(*, hash_large_inputs: bool) -> dict[str, Any]:
    require(AUTHORIZATION_PATH.is_file(), "source-only authorization is absent")
    authorization = read_json(AUTHORIZATION_PATH)
    require(authorization.get("source_only_authorized") is True, "source-only run not authorized")
    require(authorization.get("held_evaluation_authorized") is False, "held evaluation unexpectedly authorized")
    require(authorization.get("cell_eval2_authorized") is False, "cell-eval2 unexpectedly authorized")
    require(authorization.get("challenge_build_authorized") is False, "challenge build unexpectedly authorized")
    require(authorization.get("submission_authorized") is False, "submission unexpectedly authorized")
    require(authorization["protocol_sha256"] == sha256_file(PROTOCOL_PATH), "protocol hash changed")
    require(authorization["source_inventory_sha256"] == sha256_file(INVENTORY_PATH), "inventory hash changed")
    require(authorization["upstream_audit_sha256"] == sha256_file(UPSTREAM_PATH), "upstream audit hash changed")
    require(
        authorization["plan_manifest_sha256"] == sha256_file(PLAN_DIR / "plan-manifest.json"),
        "plan manifest hash changed",
    )
    for relative, expected in authorization["implementation_files"].items():
        require(sha256_file(ROOT / relative) == expected, f"implementation hash changed: {relative}")
    if hash_large_inputs:
        for relative, expected in inventory_paths().items():
            require(sha256_file(ROOT / relative) == expected, f"input hash changed: {relative}")
    validate_plan()
    return authorization


def donor_axis(donor: str) -> tuple[tuple[str, ...], tuple[str, ...]]:
    directory = ROOT / "data/donors" / donor
    return read_lines(directory / "perts.txt"), read_lines(directory / "genes.txt")


def source_indices(donor: str, targets: Sequence[str], genes: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    perts, donor_genes = donor_axis(donor)
    pert_index = {label: index for index, label in enumerate(perts)}
    gene_index = {label: index for index, label in enumerate(donor_genes)}
    require(all(target in pert_index for target in targets), f"{donor}: missing target")
    require(all(gene in gene_index for gene in genes), f"{donor}: missing response gene")
    return (
        np.asarray([pert_index[target] for target in targets], dtype=np.int64),
        np.asarray([gene_index[gene] for gene in genes], dtype=np.int64),
    )


def load_delta_rows(donor: str, targets: Sequence[str], genes: Sequence[str]) -> np.ndarray:
    rows, columns = source_indices(donor, targets, genes)
    path = ROOT / "data/donors" / donor / "delta.npy"
    matrix = np.load(path, mmap_mode="r")
    values = np.asarray(matrix[rows][:, columns], dtype=np.float64)
    require(values.shape == (len(targets), len(genes)), f"{donor}: sliced delta shape changed")
    require(np.all(np.isfinite(values)), f"{donor}: non-finite source delta")
    return values


def load_control(donor: str, genes: Sequence[str]) -> np.ndarray:
    _, columns = source_indices(donor, (), genes)
    control = np.load(ROOT / "data/donors" / donor / "b_ctrl.npy", mmap_mode="r")
    values = np.asarray(control[columns], dtype=np.float64)
    require(values.shape == (len(genes),), f"{donor}: control shape changed")
    require(np.all(np.isfinite(values)), f"{donor}: non-finite control")
    return values


def stable_pca_scores(y_centered: np.ndarray, max_k: int) -> np.ndarray:
    values = np.asarray(y_centered, dtype=np.float64)
    require(values.ndim == 2, "PCA input is not two-dimensional")
    # sklearn PCA centers each perturbation-feature over the gene/sample axis.
    pca_input = values - values.mean(axis=0, keepdims=True)
    gram = pca_input.T @ pca_input
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    selected = order[:max_k]
    require(np.all(eigenvalues[selected] > max(float(eigenvalues[order[0]]), 1.0) * 1e-12), "PCA rank unstable")
    components = eigenvectors[:, selected]
    scores = pca_input @ components
    for column in range(scores.shape[1]):
        pivot = int(np.argmax(np.abs(scores[:, column])))
        if scores[pivot, column] < 0:
            scores[:, column] *= -1
    return scores


def build_factors(donor: str, fit_targets: Sequence[str], genes: Sequence[str]) -> dict[str, np.ndarray]:
    delta = load_delta_rows(donor, fit_targets, genes)
    control = load_control(donor, genes)
    absolute = np.maximum(control[None, :] + delta, 0.0)
    y = absolute.T
    bias = y.mean(axis=1)
    y_centered = y - bias[:, None]
    G = stable_pca_scores(y_centered, MAX_K)
    mean_delta = delta.mean(axis=0)
    centered = delta - mean_delta
    scale = float(np.median(np.linalg.norm(centered, axis=1)))
    require(math.isfinite(scale) and scale > 0, f"{donor}: invalid fit effect scale")
    return {
        "bias": bias,
        "control": control,
        "mean_delta": mean_delta,
        "G": G,
        "gt_y": G.T @ y_centered,
        "scale": np.asarray([scale], dtype=np.float64),
    }


def solve_w(factors: Mapping[str, np.ndarray], fit_embeddings: np.ndarray, K: int, ridge: float) -> np.ndarray:
    G = np.asarray(factors["G"][:, :K], dtype=np.float64)
    gt_y = np.asarray(factors["gt_y"][:K], dtype=np.float64)
    P = np.asarray(fit_embeddings, dtype=np.float64)
    z = np.linalg.solve(G.T @ G + ridge * np.eye(K), gt_y)
    return np.linalg.solve(P.T @ P + ridge * np.eye(P.shape[1]), (z @ P).T).T


def predict_effect(
    factors: Mapping[str, np.ndarray],
    fit_embeddings: np.ndarray,
    query_embeddings: np.ndarray,
    K: int,
    ridge: float,
    *,
    centered: bool,
) -> np.ndarray:
    W = solve_w(factors, fit_embeddings, K, ridge)
    shift = (
        np.asarray(query_embeddings, dtype=np.float64)
        @ W.T
        @ np.asarray(factors["G"][:, :K], dtype=np.float64).T
    )
    if centered:
        return shift
    absolute = np.asarray(factors["bias"], dtype=np.float64)[None, :] + shift
    absolute = np.maximum(absolute, 0.0)
    effect = absolute - np.asarray(factors["control"], dtype=np.float64)[None, :]
    return effect


def row_unit(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    norms = np.linalg.norm(array, axis=1, keepdims=True)
    return np.divide(array, norms, out=np.zeros_like(array), where=norms > EPS)


def direction_blend(incumbent: np.ndarray, slim: np.ndarray, fraction: float) -> np.ndarray:
    base = np.asarray(incumbent, dtype=np.float64)
    learned = np.asarray(slim, dtype=np.float64)
    require(base.shape == learned.shape and base.ndim == 2, "blend arrays do not align")
    if fraction == 0.0:
        return base.copy()
    base_norm = np.linalg.norm(base, axis=1)
    learned_norm = np.linalg.norm(learned, axis=1)
    usable = (base_norm > EPS) & (learned_norm > EPS)
    result = base.copy()
    mixed = (1.0 - fraction) * row_unit(base[usable]) + fraction * row_unit(learned[usable])
    mixed_norm = np.linalg.norm(mixed, axis=1)
    nonzero = mixed_norm > EPS
    rows = np.flatnonzero(usable)[nonzero]
    result[rows] = mixed[nonzero] / mixed_norm[nonzero, None] * base_norm[rows, None]
    return result


def pds_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    predicted = np.asarray(prediction, dtype=np.float64)
    observed = np.asarray(truth, dtype=np.float64)
    require(predicted.shape == observed.shape and predicted.ndim == 2, "PDS arrays do not align")
    require(predicted.shape[0] >= 2, "PDS gallery too small")
    truth_unit = row_unit(observed)
    result = np.full(predicted.shape[0], 0.5, dtype=np.float64)
    norms = np.linalg.norm(predicted, axis=1)
    rows = np.flatnonzero(norms > EPS)
    if len(rows):
        scores = predicted[rows] / norms[rows, None] @ truth_unit.T
        own = scores[np.arange(len(rows)), rows][:, None]
        rank = (scores > own).sum(axis=1).astype(np.float64)
        rank += 0.5 * ((scores == own).sum(axis=1) - 1)
        result[rows] = 1.0 - rank / (predicted.shape[0] - 1)
    return result


def pds_per_target_from_scores(scores: np.ndarray) -> np.ndarray:
    matrix = np.asarray(scores, dtype=np.float64)
    require(matrix.ndim == 2 and matrix.shape[0] == matrix.shape[1], "PDS score matrix must be square")
    require(matrix.shape[0] >= 2, "PDS gallery too small")
    own = np.diag(matrix)[:, None]
    rank = (matrix > own).sum(axis=1).astype(np.float64)
    rank += 0.5 * ((matrix == own).sum(axis=1) - 1)
    return 1.0 - rank / (matrix.shape[1] - 1)


def signed_credit_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    pred = np.asarray(prediction, dtype=np.float64)
    observed = np.asarray(truth, dtype=np.float64)
    denom = np.abs(pred).sum(axis=1)
    numer = (np.sign(observed) * pred).sum(axis=1)
    return np.divide(numer, denom, out=np.zeros_like(numer), where=denom > EPS)


def top20_precision_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    pred = np.asarray(prediction, dtype=np.float64)
    observed = np.asarray(truth, dtype=np.float64)
    count = max(1, int(math.ceil(0.20 * pred.shape[1])))
    output = np.empty(pred.shape[0], dtype=np.float64)
    for row in range(pred.shape[0]):
        order = np.argsort(-np.abs(pred[row]), kind="stable")[:count]
        weights = np.abs(observed[row, order])
        correct = (
            (np.sign(pred[row, order]) == np.sign(observed[row, order]))
            & (np.sign(pred[row, order]) != 0)
            & (np.sign(observed[row, order]) != 0)
        )
        total = float(weights.sum())
        output[row] = float(weights[correct].sum() / total) if total > EPS else 0.0
    return output


def metric_bundle(prediction: np.ndarray, truth: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "pds": pds_per_target(prediction, truth),
        "signed_credit": signed_credit_per_target(prediction, truth),
        "top20_precision": top20_precision_per_target(prediction, truth),
    }


def candidate_prediction(candidate: Candidate, incumbent: np.ndarray, slim: np.ndarray) -> np.ndarray:
    if candidate.kind == "identity":
        return np.asarray(incumbent, dtype=np.float64).copy()
    require(candidate.kind == "blend", f"unknown candidate kind {candidate.kind}")
    return direction_blend(incumbent, slim, candidate.blend)


def candidate_preference(candidate: Candidate) -> tuple[float, float, float]:
    if candidate.kind == "identity":
        return (0.0, 0.0, 0.0)
    return (
        -candidate.blend,
        -float(candidate.K),
        -abs(math.log10(candidate.ridge / 0.1)),
    )


def select_candidate_deltas(
    cell_deltas: Mapping[str, Mapping[str, float]],
) -> tuple[Candidate, dict[str, Any]]:
    rows = []
    winner = candidates()[0]
    winner_key = (1.0, 0.0, 0.0, *candidate_preference(winner))
    for candidate in candidates():
        deltas = [float(cell_deltas[cell][candidate.key]) for cell in sorted(cell_deltas)]
        eligible = min(deltas) >= 0.0
        key = (
            float(eligible),
            float(min(deltas)),
            float(np.mean(deltas)),
            *candidate_preference(candidate),
        )
        if key > winner_key:
            winner, winner_key = candidate, key
        rows.append(
            {
                "candidate": candidate.key,
                "eligible": bool(eligible),
                "worst_pds_delta": float(min(deltas)),
                "macro_pds_delta": float(np.mean(deltas)),
                "cell_pds_deltas": {
                    cell: float(value)
                    for cell, value in zip(sorted(cell_deltas), deltas, strict=True)
                },
            }
        )
    return winner, {
        "selected": winner.key,
        "selected_key": list(winner_key),
        "candidates": rows,
    }


def source_norm_factor(factors: Mapping[str, Mapping[str, np.ndarray]], sources: Sequence[str]) -> dict[str, float]:
    scales = {source: float(factors[source]["scale"][0]) for source in sources}
    reference = float(np.median(list(scales.values())))
    return {source: reference / scales[source] for source in sources}


def averaged_source_effect(
    values: Mapping[str, np.ndarray],
    factors: Mapping[str, Mapping[str, np.ndarray]],
    sources: Sequence[str],
) -> np.ndarray:
    weights = source_norm_factor(factors, sources)
    return np.mean(
        np.stack([weights[source] * values[source] for source in sources], axis=0),
        axis=0,
    )


def candidate_pds_deltas(
    factors: Mapping[str, np.ndarray],
    fit_embeddings: np.ndarray,
    query_embeddings: np.ndarray,
    incumbent: np.ndarray,
    truth: np.ndarray,
) -> dict[str, float]:
    base = np.asarray(incumbent, dtype=np.float64)
    observed = np.asarray(truth, dtype=np.float64)
    base_norm = np.linalg.norm(base, axis=1)
    truth_unit = row_unit(observed)
    base_unit = row_unit(base)
    base_scores = base_unit @ truth_unit.T
    base_pds = pds_per_target_from_scores(base_scores)
    result = {"identity": 0.0}
    for K in K_GRID:
        G = np.asarray(factors["G"][:, :K], dtype=np.float64)
        gtg = G.T @ G
        base_by_G = base_unit @ G
        for ridge in RIDGE_GRID:
            W = solve_w(factors, fit_embeddings, K, ridge)
            coefficient = np.asarray(query_embeddings, dtype=np.float64) @ W.T
            slim_norm_sq = np.einsum(
                "ik,kl,il->i", coefficient, gtg, coefficient, optimize=True
            )
            slim_norm = np.sqrt(np.maximum(slim_norm_sq, 0.0))
            slim_truth_numer = coefficient @ (G.T @ truth_unit.T)
            slim_scores = np.divide(
                slim_truth_numer,
                slim_norm[:, None],
                out=np.zeros_like(slim_truth_numer),
                where=slim_norm[:, None] > EPS,
            )
            base_slim_numer = np.sum(base_by_G * coefficient, axis=1)
            base_slim_cosine = np.divide(
                base_slim_numer,
                slim_norm,
                out=np.zeros_like(base_slim_numer),
                where=slim_norm > EPS,
            )
            usable = (base_norm > EPS) & (slim_norm > EPS)
            for fraction in BLEND_GRID:
                scores = base_scores.copy()
                mixed_norm = np.sqrt(
                    np.maximum(
                        (1.0 - fraction) ** 2
                        + fraction**2
                        + 2.0 * fraction * (1.0 - fraction) * base_slim_cosine,
                        0.0,
                    )
                )
                valid = usable & (mixed_norm > EPS)
                scores[valid] = (
                    (1.0 - fraction) * base_scores[valid]
                    + fraction * slim_scores[valid]
                ) / mixed_norm[valid, None]
                candidate = Candidate("blend", K, ridge, fraction)
                result[candidate.key] = float(
                    np.mean(pds_per_target_from_scores(scores) - base_pds)
                )
    require(set(result) == {candidate.key for candidate in candidates()}, "candidate PDS table incomplete")
    return result


def tune_cell_deltas(
    factors: Mapping[str, Mapping[str, np.ndarray]],
    fit_embeddings: np.ndarray,
    tune_embeddings: np.ndarray,
    tune_direct: Mapping[str, np.ndarray],
) -> dict[str, dict[str, float]]:
    output: dict[str, dict[str, float]] = {}
    for source in DONORS:
        incumbent = tune_direct[source] - factors[source]["mean_delta"][None, :]
        for held in DONORS:
            if source == held:
                continue
            output[f"{source}_to_{held}"] = candidate_pds_deltas(
                factors[source],
                fit_embeddings,
                tune_embeddings,
                incumbent,
                tune_direct[held],
            )
    return output


def outer_cell_deltas(
    outer: str,
    all_cells: Mapping[str, Mapping[str, float]],
) -> dict[str, Mapping[str, float]]:
    sources = [donor for donor in DONORS if donor != outer]
    return {
        f"{sources[0]}_to_{sources[1]}": all_cells[f"{sources[0]}_to_{sources[1]}"],
        f"{sources[1]}_to_{sources[0]}": all_cells[f"{sources[1]}_to_{sources[0]}"],
    }


def serialize_factors(path: Path, factors: Mapping[str, Mapping[str, np.ndarray]]) -> None:
    arrays: dict[str, np.ndarray] = {}
    for donor, values in factors.items():
        for name, value in values.items():
            arrays[f"{donor}__{name}"] = np.asarray(value)
    np.savez_compressed(path, **arrays)


def deserialize_factors(path: Path) -> dict[str, dict[str, np.ndarray]]:
    with np.load(path, allow_pickle=False) as handle:
        return {
            donor: {
                name: np.asarray(handle[f"{donor}__{name}"], dtype=np.float64)
                for name in ("bias", "control", "mean_delta", "G", "gt_y", "scale")
            }
            for donor in DONORS
        }


def selection_stage() -> None:
    validate_authorization(hash_large_inputs=True)
    require(not SELECTION_DIR.exists(), f"{SELECTION_DIR} exists; refusing overwrite")
    require(not RESULT_DIR.exists(), "source result already exists")
    fit_targets = read_lines(PLAN_DIR / "fit_targets.txt")
    tune_targets = read_lines(PLAN_DIR / "tune_targets.txt")
    score_genes = read_lines(PLAN_DIR / "score_genes.txt")
    embeddings = np.load(PLAN_DIR / "all_target_embeddings.npy", allow_pickle=False)
    fit_embeddings = np.asarray(embeddings[:FIT_COUNT], dtype=np.float64)
    tune_embeddings = np.asarray(embeddings[FIT_COUNT : FIT_COUNT + TUNE_COUNT], dtype=np.float64)

    factors = {donor: build_factors(donor, fit_targets, score_genes) for donor in DONORS}
    tune_direct = {donor: load_delta_rows(donor, tune_targets, score_genes) for donor in DONORS}
    all_cell_deltas = tune_cell_deltas(
        factors, fit_embeddings, tune_embeddings, tune_direct
    )

    outer_selection: dict[str, dict[str, Any]] = {}
    selected_candidates: dict[str, Candidate] = {}
    for outer in DONORS:
        selected, audit = select_candidate_deltas(
            outer_cell_deltas(outer, all_cell_deltas)
        )
        selected_candidates[outer] = selected
        outer_selection[outer] = audit
    deployment, deployment_audit = select_candidate_deltas(all_cell_deltas)

    SELECTION_DIR.mkdir(parents=True)
    factors_path = SELECTION_DIR / "source_factors.npz"
    serialize_factors(factors_path, factors)
    selection = {
        "schema_version": 1,
        "experiment": "123-slim-string-reproduction",
        "sealed_before_test_truth": True,
        "outer_selected": {
            donor: asdict(candidate) | {"key": candidate.key}
            for donor, candidate in selected_candidates.items()
        },
        "deployment_selected": asdict(deployment) | {"key": deployment.key},
        "outer_selection_audit": outer_selection,
        "deployment_selection_audit": deployment_audit,
        "fit_target_order_sha256": sha256_lines(fit_targets),
        "tune_target_order_sha256": sha256_lines(tune_targets),
        "score_gene_order_sha256": sha256_lines(score_genes),
        "source_factors_sha256": sha256_file(factors_path),
    }
    selection_path = SELECTION_DIR / "selection.json"
    write_json(selection_path, selection)
    manifest = {
        "schema_version": 1,
        "selection_sha256": sha256_file(selection_path),
        "source_factors_sha256": sha256_file(factors_path),
        "authorization_sha256": sha256_file(AUTHORIZATION_PATH),
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "plan_manifest_sha256": sha256_file(PLAN_DIR / "plan-manifest.json"),
        "implementation_sha256": sha256_file(Path(__file__)),
    }
    write_json(SELECTION_DIR / "freeze-manifest.json", manifest)
    print(json.dumps(selection, indent=2, sort_keys=True))


def verify_selection() -> tuple[dict[str, Any], dict[str, dict[str, np.ndarray]]]:
    validate_authorization(hash_large_inputs=False)
    manifest = read_json(SELECTION_DIR / "freeze-manifest.json")
    require(manifest["authorization_sha256"] == sha256_file(AUTHORIZATION_PATH), "selection authorization changed")
    require(manifest["protocol_sha256"] == sha256_file(PROTOCOL_PATH), "selection protocol changed")
    require(manifest["implementation_sha256"] == sha256_file(Path(__file__)), "selection implementation changed")
    require(
        manifest["plan_manifest_sha256"] == sha256_file(PLAN_DIR / "plan-manifest.json"),
        "selection plan changed",
    )
    selection_path = SELECTION_DIR / "selection.json"
    factors_path = SELECTION_DIR / "source_factors.npz"
    require(manifest["selection_sha256"] == sha256_file(selection_path), "selection hash changed")
    require(manifest["source_factors_sha256"] == sha256_file(factors_path), "source factors changed")
    selection = read_json(selection_path)
    require(selection.get("sealed_before_test_truth") is True, "selection was not sealed before test")
    return selection, deserialize_factors(factors_path)


def candidate_from_json(value: Mapping[str, Any]) -> Candidate:
    return Candidate(str(value["kind"]), int(value["K"]), float(value["ridge"]), float(value["blend"]))


def build_outer_prediction(
    outer: str,
    candidate: Candidate,
    factors: Mapping[str, Mapping[str, np.ndarray]],
    fit_embeddings: np.ndarray,
    query_embeddings: np.ndarray,
    query_direct: Mapping[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    sources = [donor for donor in DONORS if donor != outer]
    centered_direct = {
        source: query_direct[source] - factors[source]["mean_delta"][None, :]
        for source in sources
    }
    incumbent = averaged_source_effect(centered_direct, factors, sources)
    centered_slim = {
        source: predict_effect(
            factors[source],
            fit_embeddings,
            query_embeddings,
            candidate.K if candidate.kind != "identity" else 10,
            candidate.ridge if candidate.kind != "identity" else 0.1,
            centered=True,
        )
        for source in sources
    }
    slim_adapted = averaged_source_effect(centered_slim, factors, sources)
    upstream_full = np.mean(
        np.stack(
            [
                predict_effect(
                    factors[source],
                    fit_embeddings,
                    query_embeddings,
                    candidate.K if candidate.kind != "identity" else 10,
                    candidate.ridge if candidate.kind != "identity" else 0.1,
                    centered=False,
                )
                for source in sources
            ],
            axis=0,
        ),
        axis=0,
    )
    prediction = candidate_prediction(candidate, incumbent, slim_adapted)
    return incumbent, prediction, upstream_full, slim_adapted


def summarize_outer(
    incumbent: np.ndarray,
    prediction: np.ndarray,
    upstream: np.ndarray,
    truth: np.ndarray,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    base = metric_bundle(incumbent, truth)
    candidate = metric_bundle(prediction, truth)
    upstream_metrics = metric_bundle(upstream, truth)
    per_target = {
        f"{metric}_incumbent": base[metric]
        for metric in base
    }
    per_target.update({f"{metric}_candidate": candidate[metric] for metric in candidate})
    per_target.update({f"{metric}_upstream": upstream_metrics[metric] for metric in upstream_metrics})
    per_target.update(
        {
            f"{metric}_delta": candidate[metric] - base[metric]
            for metric in candidate
        }
    )
    summary = {
        "incumbent": {metric: float(np.mean(values)) for metric, values in base.items()},
        "candidate": {metric: float(np.mean(values)) for metric, values in candidate.items()},
        "exact_upstream": {
            metric: float(np.mean(values)) for metric, values in upstream_metrics.items()
        },
        "delta": {
            metric: float(np.mean(candidate[metric] - base[metric]))
            for metric in candidate
        },
    }
    return summary, per_target


def bootstrap_intervals(
    per_outer: Mapping[str, Mapping[str, np.ndarray]],
    cluster_labels: np.ndarray,
    cluster_draws: np.ndarray,
) -> tuple[dict[str, dict[str, float]], dict[str, np.ndarray]]:
    cluster_count = int(cluster_labels.max()) + 1
    metrics = ("pds", "signed_credit", "top20_precision")
    draws_by_metric = {
        metric: np.empty(len(cluster_draws), dtype=np.float64) for metric in metrics
    }
    for draw_index, sampled in enumerate(cluster_draws):
        multiplicity = np.bincount(sampled, minlength=cluster_count)
        weights = multiplicity[cluster_labels].astype(np.float64)
        denominator = float(weights.sum())
        require(denominator > 0, "empty bootstrap draw")
        for metric in metrics:
            donor_values = [
                float(np.sum(per_outer[donor][f"{metric}_delta"] * weights) / denominator)
                for donor in DONORS
            ]
            draws_by_metric[metric][draw_index] = float(np.mean(donor_values))
    intervals = {
        metric: {
            "lower_95": float(np.quantile(draws, 0.025, method="linear")),
            "median": float(np.quantile(draws, 0.5, method="linear")),
            "upper_95": float(np.quantile(draws, 0.975, method="linear")),
        }
        for metric, draws in draws_by_metric.items()
    }
    return intervals, draws_by_metric


def embeddings_for_null(
    kind: str,
    draw: int,
    real_embeddings: np.ndarray,
    permutations: np.ndarray,
    gaussian_seeds: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if kind == "target_pairing":
        fit = np.asarray(real_embeddings[:FIT_COUNT], dtype=np.float64)[permutations[draw]]
        tune = np.asarray(real_embeddings[FIT_COUNT : FIT_COUNT + TUNE_COUNT], dtype=np.float64)
        test = np.asarray(real_embeddings[-TEST_COUNT:], dtype=np.float64)
        return fit, tune, test
    require(kind == "gaussian_embedding", f"unknown null kind {kind}")
    rng = np.random.default_rng(int(gaussian_seeds[draw]))
    real_fit = np.asarray(real_embeddings[:FIT_COUNT], dtype=np.float64)
    mean = real_fit.mean(axis=0)
    scale = real_fit.std(axis=0)
    scale[scale < EPS] = 1.0
    generated = rng.normal(size=real_embeddings.shape) * scale + mean
    return (
        generated[:FIT_COUNT],
        generated[FIT_COUNT : FIT_COUNT + TUNE_COUNT],
        generated[-TEST_COUNT:],
    )


def evaluate_scenario(
    factors: Mapping[str, Mapping[str, np.ndarray]],
    fit_embeddings: np.ndarray,
    tune_embeddings: np.ndarray,
    test_embeddings: np.ndarray,
    tune_direct: Mapping[str, np.ndarray],
    test_direct: Mapping[str, np.ndarray],
    *,
    frozen_real_selection: Mapping[str, Any] | None,
) -> tuple[dict[str, Any], dict[str, dict[str, np.ndarray]]]:
    all_cell_deltas = tune_cell_deltas(
        factors, fit_embeddings, tune_embeddings, tune_direct
    )
    selected_by_outer: dict[str, Candidate] = {}
    selection_audit: dict[str, Any] = {}
    for outer in DONORS:
        if frozen_real_selection is None:
            selected, audit = select_candidate_deltas(
                outer_cell_deltas(outer, all_cell_deltas)
            )
        else:
            selected = candidate_from_json(frozen_real_selection["outer_selected"][outer])
            audit = {"selected": selected.key, "source": "frozen_real_selection"}
        selected_by_outer[outer] = selected
        selection_audit[outer] = audit

    summaries: dict[str, Any] = {}
    per_outer: dict[str, dict[str, np.ndarray]] = {}
    for outer in DONORS:
        incumbent, prediction, upstream, _ = build_outer_prediction(
            outer,
            selected_by_outer[outer],
            factors,
            fit_embeddings,
            test_embeddings,
            test_direct,
        )
        summary, per_target = summarize_outer(
            incumbent,
            prediction,
            upstream,
            test_direct[outer],
        )
        summaries[outer] = summary
        per_outer[outer] = per_target
    macro = {
        metric: float(np.mean([summaries[donor]["delta"][metric] for donor in DONORS]))
        for metric in ("pds", "signed_credit", "top20_precision")
    }
    return {
        "outer_selection": {
            donor: asdict(candidate) | {"key": candidate.key}
            for donor, candidate in selected_by_outer.items()
        },
        "selection_audit": selection_audit,
        "outer": summaries,
        "macro_delta": macro,
    }, per_outer


def null_distributions(
    factors: Mapping[str, Mapping[str, np.ndarray]],
    real_embeddings: np.ndarray,
    tune_direct: Mapping[str, np.ndarray],
    test_direct: Mapping[str, np.ndarray],
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    permutations = np.load(PLAN_DIR / "target_pairing_permutations.npy", allow_pickle=False)
    gaussian_seeds = np.load(PLAN_DIR / "gaussian_seeds.npy", allow_pickle=False)
    metrics = ("pds", "signed_credit", "top20_precision")
    arrays: dict[str, np.ndarray] = {}
    summary: dict[str, Any] = {}
    for kind, count in (
        ("target_pairing", PAIRING_NULL_DRAWS),
        ("gaussian_embedding", GAUSSIAN_NULL_DRAWS),
    ):
        draws = {metric: np.empty(count, dtype=np.float64) for metric in metrics}
        selected_keys: list[dict[str, str]] = []
        for draw in range(count):
            fit_emb, tune_emb, test_emb = embeddings_for_null(
                kind, draw, real_embeddings, permutations, gaussian_seeds
            )
            scenario, _ = evaluate_scenario(
                factors,
                fit_emb,
                tune_emb,
                test_emb,
                tune_direct,
                test_direct,
                frozen_real_selection=None,
            )
            selected_keys.append(
                {
                    donor: str(scenario["outer_selection"][donor]["key"])
                    for donor in DONORS
                }
            )
            for metric in metrics:
                draws[metric][draw] = scenario["macro_delta"][metric]
                arrays[f"{kind}__{metric}"] = draws[metric]
        summary[kind] = {
            "draws": count,
            "selected_keys": selected_keys,
            "statistics": {
                metric: {
                    "p95": float(np.quantile(values, 0.95, method="linear")),
                    "maximum": float(np.max(values)),
                    "mean": float(np.mean(values)),
                }
                for metric, values in draws.items()
            },
        }
    return summary, arrays


def evaluate_stage() -> None:
    validate_authorization(hash_large_inputs=True)
    selection, factors = verify_selection()
    require(not RESULT_DIR.exists(), f"{RESULT_DIR} exists; refusing overwrite")
    require(not TEST_ACCESS_LOCK.exists(), f"{TEST_ACCESS_LOCK} exists; source test access is consumed")
    tune_targets = read_lines(PLAN_DIR / "tune_targets.txt")
    test_targets = read_lines(PLAN_DIR / "test_targets.txt")
    score_genes = read_lines(PLAN_DIR / "score_genes.txt")
    real_embeddings = np.load(PLAN_DIR / "all_target_embeddings.npy", allow_pickle=False)
    fit_embeddings = np.asarray(real_embeddings[:FIT_COUNT], dtype=np.float64)
    tune_embeddings = np.asarray(real_embeddings[FIT_COUNT : FIT_COUNT + TUNE_COUNT], dtype=np.float64)
    test_embeddings = np.asarray(real_embeddings[-TEST_COUNT:], dtype=np.float64)
    write_json(
        TEST_ACCESS_LOCK,
        {
            "schema_version": 1,
            "experiment": "123-slim-string-reproduction",
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "meaning": "source test-row access is permanently consumed for the frozen candidate",
            "selection_freeze_sha256": sha256_file(SELECTION_DIR / "freeze-manifest.json"),
            "authorization_sha256": sha256_file(AUTHORIZATION_PATH),
            "implementation_sha256": sha256_file(Path(__file__)),
        },
    )
    # This is the first Experiment 123 access to source test rows.
    tune_direct = {donor: load_delta_rows(donor, tune_targets, score_genes) for donor in DONORS}
    test_direct = {donor: load_delta_rows(donor, test_targets, score_genes) for donor in DONORS}

    observed, per_outer = evaluate_scenario(
        factors,
        fit_embeddings,
        tune_embeddings,
        test_embeddings,
        tune_direct,
        test_direct,
        frozen_real_selection=selection,
    )
    cluster_labels = np.load(PLAN_DIR / "test_target_cluster_labels.npy", allow_pickle=False)
    cluster_draws = np.load(PLAN_DIR / "bootstrap_cluster_draws.npy", allow_pickle=False)
    intervals, bootstrap_arrays = bootstrap_intervals(per_outer, cluster_labels, cluster_draws)
    null_summary, null_arrays = null_distributions(
        factors, real_embeddings, tune_direct, test_direct
    )

    macro = observed["macro_delta"]
    checks: dict[str, bool] = {}
    checks["selected_nonidentity_every_outer"] = all(
        observed["outer_selection"][donor]["kind"] != "identity" for donor in DONORS
    )
    for donor in DONORS:
        checks[f"{donor}_pds_delta_nonnegative"] = observed["outer"][donor]["delta"]["pds"] >= 0.0
        checks[f"{donor}_signed_credit_delta_nonnegative"] = (
            observed["outer"][donor]["delta"]["signed_credit"] >= 0.0
        )
        checks[f"{donor}_top20_precision_delta_nonnegative"] = (
            observed["outer"][donor]["delta"]["top20_precision"] >= 0.0
        )
        checks[f"{donor}_upstream_pds_at_least_chance"] = (
            observed["outer"][donor]["exact_upstream"]["pds"] >= 0.5
        )
    checks["macro_pds_delta_at_least_0p003"] = macro["pds"] >= 0.003
    checks["macro_signed_credit_delta_positive"] = macro["signed_credit"] > 0.0
    checks["macro_top20_precision_delta_positive"] = macro["top20_precision"] > 0.0
    for metric in ("pds", "signed_credit", "top20_precision"):
        checks[f"{metric}_bootstrap_lower_positive"] = intervals[metric]["lower_95"] > 0.0
        for kind in ("target_pairing", "gaussian_embedding"):
            checks[f"{metric}_above_{kind}_p95"] = (
                macro[metric] > null_summary[kind]["statistics"][metric]["p95"]
            )
    checks["complete_gallery_and_coverage"] = all(
        len(per_outer[donor]["pds_candidate"]) == TEST_COUNT for donor in DONORS
    )
    passed = all(checks.values())

    RESULT_DIR.mkdir(parents=True)
    arrays: dict[str, np.ndarray] = {
        "test_targets": np.asarray(test_targets, dtype="U"),
        "score_genes": np.asarray(score_genes, dtype="U"),
        "cluster_labels": np.asarray(cluster_labels),
    }
    for donor, metrics in per_outer.items():
        for name, values in metrics.items():
            arrays[f"{donor}__{name}"] = values
    for metric, values in bootstrap_arrays.items():
        arrays[f"bootstrap__{metric}"] = values
    arrays.update(null_arrays)
    arrays_path = RESULT_DIR / "result-arrays.npz"
    np.savez_compressed(arrays_path, **arrays)
    result = {
        "schema_version": 1,
        "experiment": "123-slim-string-reproduction",
        "status": "GO" if passed else "NO-GO",
        "source_only": True,
        "held_h1_or_orion_accessed": False,
        "observed": observed,
        "cluster_bootstrap": intervals,
        "nulls": null_summary,
        "promotion_checks": checks,
        "promotion_passed": passed,
        "candidate_freeze_sha256": sha256_file(SELECTION_DIR / "freeze-manifest.json"),
        "source_test_access_lock_sha256": sha256_file(TEST_ACCESS_LOCK),
        "result_arrays_sha256": sha256_file(arrays_path),
    }
    result_path = RESULT_DIR / "results.json"
    write_json(result_path, result)
    manifest = {
        "schema_version": 1,
        "experiment": "123-slim-string-reproduction",
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
    selection, _ = verify_selection()
    require(selection.get("sealed_before_test_truth") is True, "candidate freeze invalid")
    manifest = read_json(RESULT_DIR / "result-manifest.json")
    require(manifest["authorization_sha256"] == sha256_file(AUTHORIZATION_PATH), "result authorization changed")
    require(manifest["protocol_sha256"] == sha256_file(PROTOCOL_PATH), "result protocol changed")
    require(manifest["implementation_sha256"] == sha256_file(Path(__file__)), "result implementation changed")
    require(
        manifest["selection_freeze_sha256"] == sha256_file(SELECTION_DIR / "freeze-manifest.json"),
        "result selection freeze changed",
    )
    require(
        manifest["source_test_access_lock_sha256"] == sha256_file(TEST_ACCESS_LOCK),
        "source test access lock changed",
    )
    for name, digest in manifest["files"].items():
        require(sha256_file(RESULT_DIR / name) == digest, f"result file changed: {name}")
    result = read_json(RESULT_DIR / "results.json")
    require(result["held_h1_or_orion_accessed"] is False, "held-access statement changed")
    return result


def audit_stage() -> None:
    authorization = validate_authorization(hash_large_inputs=True)
    report = {
        "status": "GO-for-source-only-execution",
        "candidate_count": len(candidates()),
        "plan": validate_plan(),
        "authorization_sha256": sha256_file(AUTHORIZATION_PATH),
        "implementation_sha256": sha256_file(Path(__file__)),
        "allowed_input_count": len(inventory_paths()),
        "held_evaluation_authorized": authorization["held_evaluation_authorized"],
        "cell_eval2_authorized": authorization["cell_eval2_authorized"],
        "challenge_build_authorized": authorization["challenge_build_authorized"],
        "submission_authorized": authorization["submission_authorized"],
    }
    print(json.dumps(report, indent=2, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("plan", "audit", "select", "evaluate", "verify"))
    args = parser.parse_args()
    if args.command == "plan":
        create_plan()
    elif args.command == "audit":
        audit_stage()
    elif args.command == "select":
        selection_stage()
    elif args.command == "evaluate":
        evaluate_stage()
    else:
        print(json.dumps(verify_result(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

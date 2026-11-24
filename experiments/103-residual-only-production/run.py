#!/usr/bin/env python3
"""Experiment 103: residual-only intervention in the exact production path.

The freeze command has no KOLF argument and cannot open KOLF. The score command
accepts only a sealed frozen bundle plus the held KOLF directory.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import fcntl
import hashlib
import importlib.util
import json
import os
import resource
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Callable, Iterable, Mapping, Sequence, TypeVar

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "4")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "4")

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
FROZEN = WORK / "frozen"
EXPERIMENT = "103-residual-only-production"

MAX_ALLOWED_RSS_GIB = 24.0
PANEL_TARGETS = 300
OFFICIAL_GENES = 18_533
PRODUCTION_OPERATOR_GENES = 17_942
KOLF_TARGETS = 273
KOLF_SCORE_GENES = 17_806
N_NULL = 20
BOOTSTRAPS = 20_000
BOOTSTRAP_SEED = 10_320_026
CLUSTER_K = 5
CLUSTER_MIN_COSINE = 0.80
MINIMUM_DELTA = 0.003
BASELINE_ALPHA = 0.982
BASELINE_CD4_TARGET_F = 0.10

EXPECTED_PANEL_ORDER_SHA256 = (
    "cd3d00f1cd5da660609ad41d85b53e2d5557f7f5edb3e1f64559f78f7fb3cbc4"
)
EXPECTED_KOLF_TARGET_ORDER_SHA256 = (
    "b366592072d1e0da04aaabafbcc78882d34c97ff26eb2fca5fd9a9af0084f22f"
)
EXPECTED_KOLF_SCORE_AXIS_SHA256 = (
    "280aab34a0e86f27db8ade213e000f9066c8d5f9b89edd3ea966aba061ba9ae9"
)
EXPECTED_PRODUCTION_OPERATOR_SHA256 = (
    "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
)
EXPECTED_OPERATOR_COLUMNS_SHA256 = (
    "b2cae0d381c3cd7065a3de33f1a6da5514577f2b3187b6dc84786ee789e27e4d"
)
EXPECTED_PM_SUPPORTED_TARGETS = 248
EXPECTED_PM_GENES = 9_373
EXPECTED_CD4_COMPLETE_TARGETS = 276
EXPECTED_BASELINE_PREWHITEN_SHA256 = (
    "0601336fdf1ec5c533dcad69f6cb6c09d6368182bc8c69d10d516a1ecb19bf65"
)
EXPECTED_REMOTE_BASELINE_EFFECT_SHA256 = (
    "2980d70655f2a9d4a81aa0d8cdb06f2389fca20dddf586d64c132cb05885e5f2"
)
EXPECTED_BASELINE_META_SHA256 = (
    "34223442e4e8816bb9be2b1fdfe6a563ed114c17a2893bee0ff603953f84287a"
)
EXPECTED_STAGED_BASELINE_FILE_SHA256 = (
    "6f2b0488d67a31b59aff2ce6c7868f1de443f5153f28386b6a07f7489e12eb6f"
)
EXPECTED_STAGED_BASELINE_METADATA_SHA256 = (
    "a1808f6af4fd803ed2084dbf02387d9f573b932a207683eab6110874b727ebb0"
)

PRODUCTION_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
K562 = "k562_gwps"
DONOR_PATHS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
BUNDLE = ROOT / "data/challenge"
PERTURBMATCH = ROOT / "data/external/perturbmatch"
CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
OPERATOR_PATH = ROOT / "data/donors/whitener_union_pruned.npz"
OPERATOR_METADATA_PATH = ROOT / "data/donors/whitener_union_pruned.json"
PRODUCTION_PATH = ROOT / "scripts/make_submission_pm_cd4.py"
BASELINE_META_PATH = ROOT / "submissions/a15-cd4f010.vcc.meta.json"
STAGED_BASELINE_PATH = WORK / "staged-a15-linux/a15-effects.npy"
STAGED_BASELINE_METADATA_PATH = Path(str(STAGED_BASELINE_PATH) + ".json")

EXPECTED_IMPLEMENTATION_HASHES = {
    "scripts/make_submission_pm_cd4.py": (
        "7562cf5a16fd751c55ddf248e05f2e4281f907e27be0a5bb9da959d61991b62b"
    ),
    "scripts/build_submission.py": (
        "3e0a7643f604c59b5ff8494af2fa36dc08c37303dae19bc7557ac95a8dde4e44"
    ),
    "scripts/make_submission.py": (
        "e85879effd33334cd6b1ddd189b02e1373ea6f918add30298e9cf9bebb6a0ab1"
    ),
    "src/vcc2026/donor_transfer.py": (
        "c91cb05e1204138900058be2c82a227f132e373493a3b3e15914cd18f807f3b9"
    ),
    "src/vcc2026/multi_donor.py": (
        "1daf35d97b50fb234a37255cf4136d34580f98f9d54faf175dca48461a33129a"
    ),
    "src/vcc2026/whiten.py": (
        "5dbddefde17192cda4bc7f8f7eceb202faf672de1833c42251366078e7e18bc1"
    ),
    "experiments/061-perturbmatch-hct-count/build_candidate.py": (
        "77e6fbf7edf1aeea6d39be4fe801591d0b0c2867f95042b697f0fc4042a2fae4"
    ),
    "experiments/062-cd4-post-whiten/run.py": (
        "16417a9b985d9aad8993bf2f24ce550ea1e6871d9401bb4672d11ec6d0046176"
    ),
}

EXPECTED_KOLF_HASHES = {
    "b_ctrl.npy": "dd4c30afee5b6bd7fc68c469da30a1eb771e367d186de41d13de16ca4aee7fc9",
    "delta.npy": "a8ed5c7eebb0b01a6a0a4effb0651fd0c6c706aee81cac87b6369d4d04c91d6b",
    "gene_map.csv": "d0cf4c941b26c35c769a8a8032564f81aeec73c1bf954fe6a69a3ea28109e4e6",
    "genes.txt": "2f077fda94cc07cf6cfddae28fc78440e80963e4e866f9904b91417f3189420b",
    "meta.json": "57737834c91663d154e84721ce73c0a33f4b670fea1d6bdab8ad35d0cf8b466a",
    "perts.txt": "b366592072d1e0da04aaabafbcc78882d34c97ff26eb2fca5fd9a9af0084f22f",
    "provenance.json": (
        "a35566f91d7ffedf6479c250f1cd7e00f56728b6131c758b79c1ce7c3c8bcc33"
    ),
    "source_rows_int8.npy": (
        "6c482697aae226d7dff4caf800c287893be0674a419c9d0d6537e3f3c16606d0"
    ),
    "target_row_map.csv": (
        "cab9ad14ce06d09c2dbef2ee3bedb7efcf1eaa25c4cad8488249e4ffe717aca6"
    ),
}

T = TypeVar("T")


def file_sha256(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_array(values: np.ndarray, dtype: str = "<f8") -> str:
    array = np.ascontiguousarray(values, dtype=dtype)
    return hashlib.sha256(array.tobytes()).hexdigest()


def sha256_labels(values: Iterable[str]) -> str:
    return hashlib.sha256(("\n".join(map(str, values)) + "\n").encode()).hexdigest()


def stable_seed(*parts: object) -> int:
    payload = "|".join(map(str, parts)).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def read_lines(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(line.strip() for line in path.read_text().splitlines() if line.strip())
    if not values:
        raise ValueError(f"{path}: empty axis")
    if unique and len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate axis")
    return values


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: missing {column}")
        values = tuple(str(row[column]) for row in reader)
    if not values or len(values) != len(set(values)):
        raise ValueError(f"{path}: empty or duplicate {column}")
    return values


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def validate_memory_limit(limit_gib: float) -> float:
    if not np.isfinite(limit_gib) or not 0 < limit_gib <= MAX_ALLOWED_RSS_GIB:
        raise ValueError(
            f"--max-rss-gib must be in (0, {MAX_ALLOWED_RSS_GIB:.0f}]"
        )
    return float(limit_gib)


def memory_checkpoint(stage: str, limit_gib: float) -> None:
    peak = peak_rss_gib()
    print(f"[rss] {stage}: {peak:.3f} GiB", flush=True)
    if peak >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak:.3f} GiB reached {limit_gib:.3f} GiB"
        )


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(text)
    temporary.replace(path)


def atomic_save_npy(path: Path, values: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("wb") as handle:
        np.save(handle, values, allow_pickle=False)
    temporary.replace(path)


def atomic_save_npz(path: Path, **values: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **values)
    temporary.replace(path)


@contextlib.contextmanager
def experiment_lock() -> Iterable[None]:
    WORK.mkdir(parents=True, exist_ok=True)
    with (WORK / "experiment.lock").open("a+") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("another experiment-103 process is active") from exc
        yield


@dataclass
class ProtocolState:
    phase: str = "pretruth"

    def require_pretruth_fit(self, operation: str) -> None:
        if self.phase != "pretruth":
            raise RuntimeError(f"{operation} is forbidden in protocol phase {self.phase}")

    def seal(self) -> None:
        if self.phase != "pretruth":
            raise RuntimeError(f"cannot seal from phase {self.phase}")
        self.phase = "frozen"

    def open_truth(self) -> None:
        if self.phase != "frozen":
            raise RuntimeError("truth may open only after a sealed freeze")
        self.phase = "truth-open"


def execute_after_freeze(
    freeze_builder: Callable[[], Mapping[str, object]],
    truth_loader: Callable[[Mapping[str, object]], T],
) -> tuple[Mapping[str, object], T]:
    frozen = freeze_builder()
    if not bool(frozen.get("sealed_before_truth")):
        raise ValueError("freeze is not sealed")
    return frozen, truth_loader(frozen)


def verify_hash_map(
    root: Path, expected: Mapping[str, str], label: str
) -> dict[str, str]:
    observed: dict[str, str] = {}
    for name, digest in expected.items():
        path = root / name
        if not path.is_file():
            raise FileNotFoundError(f"{label}: missing {path}")
        actual = file_sha256(path)
        if actual != digest:
            raise ValueError(f"{label}: {name} hash {actual} != {digest}")
        observed[name] = actual
    return observed


def verify_implementation_hashes() -> dict[str, str]:
    observed = verify_hash_map(ROOT, EXPECTED_IMPLEMENTATION_HASHES, "implementation")
    protocol_paths = (
        Path(__file__).resolve(),
        HERE / "README.md",
        HERE / "test_run.py",
        HERE / "export_a15_effect.py",
        ROOT / "experiments/102-compass-axis-transfer/README.md",
        ROOT / "experiments/102-compass-axis-transfer/RESULTS.md",
        ROOT / "experiments/102-compass-axis-transfer/work/summary.json",
    )
    for path in protocol_paths:
        if not path.is_file():
            raise FileNotFoundError(f"protocol provenance: missing {path}")
        observed[str(path.relative_to(ROOT))] = file_sha256(path)
    return observed


def load_production_module() -> ModuleType:
    verify_implementation_hashes()
    spec = importlib.util.spec_from_file_location("exp103_production", PRODUCTION_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import {PRODUCTION_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    previous = {
        name: os.environ.get(name)
        for name in ("VCC_ALPHA", "VCC_CD4_FRACTION", "VCC_CD4_TARGET_F")
    }
    os.environ["VCC_ALPHA"] = str(BASELINE_ALPHA)
    os.environ["VCC_CD4_FRACTION"] = "0.05"
    os.environ["VCC_CD4_TARGET_F"] = str(BASELINE_CD4_TARGET_F)
    try:
        spec.loader.exec_module(module)
    finally:
        for name, value in previous.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value
    module.assert_frozen_contract()
    expected = {
        "donors": PRODUCTION_DONORS,
        "alpha": BASELINE_ALPHA,
        "beta": 0.0,
        "gamma": 0.65,
        "lambda": 0.0,
        "pm_blend": 0.50,
        "cd4_mode": "gene_norm",
        "cd4_fraction": 0.05,
        "cd4_target_f": BASELINE_CD4_TARGET_F,
    }
    observed = {
        "donors": tuple(module.PRODUCTION_DONORS),
        "alpha": float(module.ALPHA),
        "beta": float(module.BETA),
        "gamma": float(module.WHITEN_GAMMA),
        "lambda": float(module.WHITEN_LAM),
        "pm_blend": float(module.PM_BLEND),
        "cd4_mode": str(module.CD4_MODE),
        "cd4_fraction": float(module.CD4_FRACTION),
        "cd4_target_f": module.CD4_TARGET_F,
    }
    if observed != expected:
        raise ValueError(
            "current production defaults drifted:\n"
            + json.dumps({"expected": expected, "observed": observed}, indent=2)
        )
    return module


def derive_score_columns(
    official: Sequence[str],
    panel: Sequence[str],
    hct_genes: Sequence[str],
    hek_genes: Sequence[str],
) -> np.ndarray:
    panel_set = set(map(str, panel))
    source_union = set(map(str, hct_genes)) | set(map(str, hek_genes))
    columns = np.asarray(
        [
            index
            for index, gene in enumerate(official)
            if str(gene) in source_union and str(gene) not in panel_set
        ],
        dtype=np.int64,
    )
    labels = tuple(str(official[index]) for index in columns)
    if len(columns) != KOLF_SCORE_GENES:
        raise ValueError(f"source-derived score axis {len(columns)} != {KOLF_SCORE_GENES}")
    if sha256_labels(labels) != EXPECTED_KOLF_SCORE_AXIS_SHA256:
        raise ValueError("source-derived score-axis hash drifted")
    return columns


def normalize(vector: np.ndarray, label: str) -> np.ndarray:
    value = np.asarray(vector, dtype=np.float64)
    norm = float(np.linalg.norm(value))
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError(f"{label} has zero or non-finite norm")
    return value / norm


def fit_axis_from_rows(
    raw_rows: np.ndarray,
    allowed_mask: np.ndarray,
    *,
    label: str,
    state: ProtocolState | None = None,
) -> np.ndarray:
    if state is not None:
        state.require_pretruth_fit(f"fit donor axis {label}")
    values = np.asarray(raw_rows, dtype=np.float64).copy()
    allowed = np.asarray(allowed_mask, dtype=bool)
    if values.ndim != 2 or values.shape[1] != len(allowed):
        raise ValueError("raw_rows and allowed_mask disagree")
    values[:, ~allowed] = 0.0
    return normalize(values.mean(axis=0), f"{label} donor-mean axis")


def fit_donor_axis_streaming(
    donor,
    official: Sequence[str],
    score_columns: np.ndarray,
    panel: Sequence[str],
    *,
    state: ProtocolState,
    max_rss_gib: float,
    block_rows: int,
) -> tuple[np.ndarray, dict[str, object]]:
    state.require_pretruth_fit(f"fit donor axis {getattr(donor, 'name', 'source')}")
    official_pos = {str(gene): index for index, gene in enumerate(official)}
    score_set = {str(official[index]) for index in score_columns}
    target_columns: list[int] = []
    donor_columns: list[int] = []
    for donor_column, gene in enumerate(donor.genes):
        target_column = official_pos.get(str(gene))
        if target_column is not None and str(gene) in score_set:
            target_columns.append(target_column)
            donor_columns.append(donor_column)
    tc = np.asarray(target_columns, dtype=np.int64)
    dc = np.asarray(donor_columns, dtype=np.int64)
    panel_set = set(map(str, panel))
    fit_rows = np.asarray(
        [index for index, target in enumerate(donor.perts) if str(target) not in panel_set],
        dtype=np.int64,
    )
    if len(fit_rows) < 2 or not len(tc):
        raise ValueError("donor axis has insufficient rows or coordinates")
    total = np.zeros(len(tc), dtype=np.float64)
    for start in range(0, len(fit_rows), block_rows):
        rows = fit_rows[start : start + block_rows]
        total += np.asarray(donor.delta[np.ix_(rows, dc)], dtype=np.float64).sum(0)
        memory_checkpoint("donor-axis streaming", max_rss_gib)
    mean = total / len(fit_rows)
    axis = np.zeros(len(official), dtype=np.float64)
    axis[tc] = mean
    axis = normalize(axis, "raw non-panel donor mean")
    return axis, {
        "fit_rows": int(len(fit_rows)),
        "fit_target_order_sha256": sha256_labels(
            str(donor.perts[index]) for index in fit_rows
        ),
        "mapped_score_coordinates": int(len(tc)),
        "axis_sha256": sha256_array(axis),
    }


def residualize(vector: np.ndarray, axis: np.ndarray) -> tuple[np.ndarray, float]:
    value = np.asarray(vector, dtype=np.float64)
    unit = normalize(axis, "residualization axis")
    loading = float(value @ unit)
    return value - loading * unit, loading


def source_panel_loadings(
    model,
    axes: Mapping[str, np.ndarray],
    targets: Sequence[str],
) -> tuple[np.ndarray, np.ndarray]:
    loadings = np.zeros((len(model.names), len(targets)), dtype=np.float64)
    carrier = np.zeros((len(model.names), len(targets)), dtype=bool)
    for source_index, name in enumerate(model.names):
        axis = np.asarray(axes[name], dtype=np.float64)
        for target_index, target in enumerate(targets):
            row = model.pert_index[name].get(str(target))
            if row is None:
                continue
            vector = np.zeros(len(model.target_genes), dtype=np.float64)
            vector[model.tcols[name]] = np.asarray(
                model.s[name][row][model.dcols[name]], dtype=np.float64
            )
            loadings[source_index, target_index] = float(vector @ axis)
            carrier[source_index, target_index] = True
    return loadings, carrier


def permute_source_loadings(
    loadings: np.ndarray,
    carrier: np.ndarray,
    draw: int,
) -> np.ndarray:
    values = np.asarray(loadings, dtype=np.float64)
    present = np.asarray(carrier, dtype=bool)
    if values.shape != present.shape:
        raise ValueError("loading and carrier shapes disagree")
    out = np.zeros_like(values)
    for source in range(len(values)):
        rows = np.flatnonzero(present[source])
        rng = np.random.default_rng(
            stable_seed("experiment-103-source-null", draw, source)
        )
        out[source, rows] = values[source, rng.permutation(rows)]
    return out


def projection_correction(
    model,
    targets: Sequence[str],
    axes: Mapping[str, np.ndarray],
    loadings: np.ndarray,
    carrier: np.ndarray,
    pm_supported: np.ndarray,
    pm_gene_mask: np.ndarray,
    *,
    alpha: float,
    pm_blend: float,
) -> np.ndarray:
    """Return the exact production-scaled real-donor projection to subtract."""

    correction = np.zeros(
        (len(targets), len(model.target_genes)), dtype=np.float64
    )
    pm_genes = np.asarray(pm_gene_mask, dtype=bool)
    if pm_genes.shape != (len(model.target_genes),):
        raise ValueError("PerturbMatch gene mask does not match the model axis")
    weights = {name: 1.0 for name in model.names}
    for target_index, target in enumerate(targets):
        have = [
            name for name in model.names if str(target) in model.pert_index[name]
        ]
        total_weight = float(sum(weights[name] for name in have))
        if total_weight <= 0:
            continue
        for source_index, name in enumerate(model.names):
            if not carrier[source_index, target_index]:
                continue
            coefficient = (
                alpha
                * float(model._global)
                * (weights[name] / total_weight)
                * float(model._renorm[name])
            )
            vector = coefficient * loadings[source_index, target_index] * axes[name]
            if name == K562 and bool(pm_supported[target_index]):
                retention = np.ones(len(vector), dtype=np.float64)
                retention[pm_genes] = 1.0 - pm_blend
                vector = vector * retention
            correction[target_index] += vector
    return correction


def apply_projection_correction(
    baseline_prewhiten: np.ndarray,
    correction: np.ndarray,
    *,
    enabled: bool,
) -> np.ndarray:
    baseline = np.asarray(baseline_prewhiten)
    if not enabled:
        return baseline.copy()
    result = np.asarray(baseline, dtype=np.float64) - np.asarray(
        correction, dtype=np.float64
    )
    if not np.isfinite(result).all():
        raise ValueError("projection-corrected pre-whitening matrix is non-finite")
    return result


def finish_prewhiten(
    production,
    prewhiten: np.ndarray,
    operator,
    panel_positions: np.ndarray,
    cd4,
) -> tuple[np.ndarray, dict[str, object]]:
    whitened, whitening_audit = production.apply_v9_whitening(
        prewhiten, operator.whitener, operator.columns
    )
    final, cd4_audit = production.apply_frozen_cd4_residual(
        whitened, cd4, operator.columns, panel_positions
    )
    return final, {
        "prewhiten_sha256": sha256_array(prewhiten),
        "post_whiten_sha256": sha256_array(whitened),
        "final_sha256": sha256_array(final),
        "whitening": whitening_audit,
        "cd4": cd4_audit,
    }


def verify_exact(label: str, observed: np.ndarray, expected: np.ndarray) -> None:
    if not np.array_equal(np.asarray(observed), np.asarray(expected)):
        error = float(
            np.max(
                np.abs(
                    np.asarray(observed, dtype=np.float64)
                    - np.asarray(expected, dtype=np.float64)
                )
            )
        )
        raise ValueError(f"{label} is not exact; maximum error {error:.17g}")


def row_normalize(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2 or not np.isfinite(array).all():
        raise ValueError("row normalization requires a finite matrix")
    norms = np.linalg.norm(array, axis=1, keepdims=True)
    if np.any(norms <= 1e-12):
        raise ValueError("cannot normalize a zero-norm row")
    return array / norms


def pds_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    predicted = np.asarray(prediction, dtype=np.float64)
    observed = np.asarray(truth, dtype=np.float64)
    if predicted.shape != observed.shape or predicted.ndim != 2:
        raise ValueError("PDS arrays do not align")
    similarity = row_normalize(predicted) @ row_normalize(observed).T
    own = np.diag(similarity)[:, None]
    rank = (similarity > own).sum(axis=1).astype(np.float64)
    rank += 0.5 * ((similarity == own).sum(axis=1) - 1)
    return 1.0 - rank / (len(predicted) - 1)


def source_effect_clusters(
    query: np.ndarray,
    *,
    reciprocal_k: int = CLUSTER_K,
    cosine_threshold: float = CLUSTER_MIN_COSINE,
) -> tuple[np.ndarray, dict[str, object]]:
    unit = row_normalize(query)
    similarity = unit @ unit.T
    np.fill_diagonal(similarity, -np.inf)
    k = min(max(1, reciprocal_k), len(unit) - 1)
    neighbors = np.argsort(-similarity, axis=1, kind="stable")[:, :k]
    directed = np.zeros((len(unit), len(unit)), dtype=bool)
    directed[np.arange(len(unit))[:, None], neighbors] = True
    reciprocal = directed & directed.T & (similarity >= cosine_threshold)
    parent = np.arange(len(unit), dtype=np.int64)

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = int(parent[value])
        return value

    def union(left: int, right: int) -> None:
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[max(root_left, root_right)] = min(root_left, root_right)

    left, right = np.where(np.triu(reciprocal, k=1))
    for a, b in zip(left.tolist(), right.tolist(), strict=True):
        union(a, b)
    roots = [find(index) for index in range(len(unit))]
    mapping = {root: index for index, root in enumerate(sorted(set(roots)))}
    labels = np.asarray([mapping[root] for root in roots], dtype=np.int64)
    sizes = np.bincount(labels)
    return labels, {
        "method": "baseline reciprocal-5NN connected components",
        "reciprocal_k": k,
        "minimum_cosine": cosine_threshold,
        "clusters": int(len(sizes)),
        "singletons": int(np.sum(sizes == 1)),
        "largest_cluster": int(sizes.max()),
    }


def cluster_bootstrap(
    deltas: np.ndarray,
    labels: np.ndarray,
    *,
    repetitions: int = BOOTSTRAPS,
    seed: int = BOOTSTRAP_SEED,
) -> dict[str, object]:
    values = np.asarray(deltas, dtype=np.float64)
    clusters = np.asarray(labels, dtype=np.int64)
    if values.shape != clusters.shape or values.ndim != 1:
        raise ValueError("bootstrap arrays do not align")
    unique = np.unique(clusters)
    if len(unique) < 2:
        return {
            "valid": False,
            "clusters": int(len(unique)),
            "point_estimate": float(values.mean()),
            "lower_95": None,
            "upper_95": None,
        }
    members = [np.flatnonzero(clusters == cluster) for cluster in unique]
    sums = np.asarray([values[index].sum() for index in members])
    sizes = np.asarray([len(index) for index in members])
    rng = np.random.default_rng(seed)
    draws = np.empty(repetitions, dtype=np.float64)
    for start in range(0, repetitions, 1_000):
        stop = min(start + 1_000, repetitions)
        selected = rng.integers(0, len(unique), size=(stop - start, len(unique)))
        draws[start:stop] = sums[selected].sum(1) / sizes[selected].sum(1)
    return {
        "valid": True,
        "clusters": int(len(unique)),
        "repetitions": repetitions,
        "seed": seed,
        "point_estimate": float(values.mean()),
        "lower_95": float(np.percentile(draws, 2.5)),
        "upper_95": float(np.percentile(draws, 97.5)),
    }


def write_freeze_bundle(
    staging: Path,
    *,
    official: Sequence[str],
    panel: Sequence[str],
    score_columns: np.ndarray,
    operator_columns: np.ndarray,
    baseline: np.ndarray,
    candidate: np.ndarray,
    nulls: Sequence[np.ndarray],
    axes: Mapping[str, np.ndarray],
    loadings: np.ndarray,
    carrier: np.ndarray,
    source_counts: np.ndarray,
    clusters: np.ndarray,
    pm_supported: np.ndarray,
    implementation_hashes: Mapping[str, str],
    input_hashes: Mapping[str, str],
    audit: Mapping[str, object],
    state: ProtocolState,
) -> dict[str, object]:
    state.require_pretruth_fit("write freeze bundle")
    staging.mkdir(parents=True, exist_ok=False)
    atomic_save_npy(staging / "baseline.npy", np.asarray(baseline, dtype=np.float64))
    atomic_save_npy(staging / "candidate.npy", np.asarray(candidate, dtype=np.float64))
    null_dir = staging / "null"
    null_dir.mkdir()
    for draw, values in enumerate(nulls):
        atomic_save_npy(
            null_dir / f"source-label-permutation-{draw:03d}.npy",
            np.asarray(values, dtype=np.float64),
        )
    atomic_save_npz(
        staging / "frozen-config.npz",
        official=np.asarray(official),
        panel=np.asarray(panel),
        score_columns=np.asarray(score_columns, dtype=np.int64),
        operator_columns=np.asarray(operator_columns, dtype=np.int64),
        source_names=np.asarray(PRODUCTION_DONORS),
        source_axes=np.stack([axes[name] for name in PRODUCTION_DONORS]),
        source_loadings=np.asarray(loadings, dtype=np.float64),
        source_carrier=np.asarray(carrier, dtype=bool),
        source_counts=np.asarray(source_counts, dtype=np.uint8),
        source_clusters=np.asarray(clusters, dtype=np.int64),
        pm_supported=np.asarray(pm_supported, dtype=bool),
    )
    config = {
        "experiment": EXPERIMENT,
        "targets": PANEL_TARGETS,
        "genes": OFFICIAL_GENES,
        "score_genes": KOLF_SCORE_GENES,
        "nulls": N_NULL,
        "bootstrap_repetitions": BOOTSTRAPS,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "minimum_delta": MINIMUM_DELTA,
        "cluster_k": CLUSTER_K,
        "cluster_minimum_cosine": CLUSTER_MIN_COSINE,
        "truth_representation": "stored raw KOLF effects; no target/gene centering",
        "official_compass_claim": False,
    }
    atomic_write_text(
        staging / "config.json", json.dumps(config, indent=2, sort_keys=True) + "\n"
    )
    files = {
        str(path.relative_to(staging)): file_sha256(path)
        for path in sorted(staging.rglob("*"))
        if path.is_file()
    }
    state.seal()
    manifest = {
        "schema_version": 1,
        "experiment": EXPERIMENT,
        "sealed_before_truth": True,
        "created_unix": time.time(),
        "kolf_files_opened_before_seal": 0,
        "post_truth_fitting_allowed": False,
        "implementation_hashes": dict(implementation_hashes),
        "non_kolf_input_hashes": dict(input_hashes),
        "files": files,
        "array_hashes": {
            "baseline": sha256_array(baseline),
            "candidate": sha256_array(candidate),
            "source_axes": sha256_array(
                np.stack([axes[name] for name in PRODUCTION_DONORS])
            ),
            "source_loadings": sha256_array(loadings),
            "source_carrier": sha256_array(carrier, "<u1"),
            "source_clusters": sha256_array(clusters, "<i8"),
        },
        "audit": dict(audit),
    }
    atomic_write_text(
        staging / "manifest.json",
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    return manifest


def verify_frozen_bundle(path: Path = FROZEN) -> dict[str, object]:
    manifest_path = path / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError("sealed experiment-103 freeze is absent")
    manifest = json.loads(manifest_path.read_text())
    if (
        manifest.get("experiment") != EXPERIMENT
        or manifest.get("sealed_before_truth") is not True
        or manifest.get("kolf_files_opened_before_seal") != 0
        or manifest.get("post_truth_fitting_allowed") is not False
    ):
        raise ValueError("freeze manifest protocol fields are invalid")
    expected = dict(manifest["files"])
    observed_names = {
        str(file.relative_to(path))
        for file in path.rglob("*")
        if file.is_file() and file.name != "manifest.json"
    }
    if observed_names != set(expected):
        raise ValueError("frozen file set changed")
    for name, digest in expected.items():
        if file_sha256(path / name) != digest:
            raise ValueError(f"frozen artifact hash changed: {name}")
    return manifest


def freeze(max_rss_gib: float, block_rows: int) -> None:
    if FROZEN.exists():
        verify_frozen_bundle()
        print(FROZEN / "manifest.json")
        return
    state = ProtocolState()
    production = load_production_module()
    implementation_hashes = verify_implementation_hashes()
    input_hashes = production.collect_input_hashes(
        BUNDLE,
        DONOR_PATHS,
        PERTURBMATCH,
        CD4_ROOT,
        OPERATOR_PATH,
        OPERATOR_METADATA_PATH,
        whitener_id="union_pruned",
    )
    if file_sha256(BASELINE_META_PATH) != EXPECTED_BASELINE_META_SHA256:
        raise ValueError("a15 baseline metadata hash changed")
    input_hashes[str(BASELINE_META_PATH.relative_to(ROOT))] = (
        EXPECTED_BASELINE_META_SHA256
    )
    if file_sha256(STAGED_BASELINE_PATH) != EXPECTED_STAGED_BASELINE_FILE_SHA256:
        raise ValueError("preserved staged-a15 effect file hash changed")
    if (
        file_sha256(STAGED_BASELINE_METADATA_PATH)
        != EXPECTED_STAGED_BASELINE_METADATA_SHA256
    ):
        raise ValueError("preserved staged-a15 effect metadata hash changed")
    input_hashes[str(STAGED_BASELINE_PATH.relative_to(ROOT))] = (
        EXPECTED_STAGED_BASELINE_FILE_SHA256
    )
    input_hashes[str(STAGED_BASELINE_METADATA_PATH.relative_to(ROOT))] = (
        EXPECTED_STAGED_BASELINE_METADATA_SHA256
    )
    memory_checkpoint("non-KOLF provenance", max_rss_gib)

    official = read_csv_column(BUNDLE / "gene_names.csv", "gene_name")
    panel = read_csv_column(BUNDLE / "pert_counts.csv", "target_gene")
    if len(official) != OFFICIAL_GENES or len(panel) != PANEL_TARGETS:
        raise ValueError("challenge axes changed")
    if sha256_labels(panel) != EXPECTED_PANEL_ORDER_SHA256:
        raise ValueError("panel order changed")
    score_columns = derive_score_columns(
        official,
        panel,
        read_lines(DONOR_PATHS["orion_hct116"] / "genes.txt", unique=False),
        read_lines(DONOR_PATHS["orion_hek293t"] / "genes.txt", unique=False),
    )
    score_mask = np.zeros(len(official), dtype=bool)
    score_mask[score_columns] = True

    donors = {
        name: production.Donor.load(path) for name, path in DONOR_PATHS.items()
    }
    model = production.MultiDonorModel(
        donors, np.asarray(official), gene_norm=False, eb=False
    )
    weights = production.donor_weights(model)
    axes: dict[str, np.ndarray] = {}
    axis_audit: dict[str, object] = {}
    for name in PRODUCTION_DONORS:
        axes[name], axis_audit[name] = fit_donor_axis_streaming(
            donors[name],
            official,
            score_columns,
            panel,
            state=state,
            max_rss_gib=max_rss_gib,
            block_rows=block_rows,
        )

    pm = production.EXP061.prepare_perturbmatch(
        PERTURBMATCH,
        tuple(official),
        chunk_rows=block_rows,
        max_rss_gib=max_rss_gib,
    )
    operator = production.load_production_operator(
        model,
        panel,
        OPERATOR_PATH,
        OPERATOR_METADATA_PATH,
        expected_sha256=EXPECTED_PRODUCTION_OPERATOR_SHA256,
    )
    if len(operator.columns) != PRODUCTION_OPERATOR_GENES:
        raise ValueError("production operator column count changed")
    if sha256_array(operator.columns, "<i8") != EXPECTED_OPERATOR_COLUMNS_SHA256:
        raise ValueError("production operator column mask changed")
    official_pos = {gene: index for index, gene in enumerate(official)}
    panel_positions = np.asarray(
        [official_pos[target] for target in panel], dtype=np.int64
    )
    if np.intersect1d(operator.columns, panel_positions).size:
        raise AssertionError("panel coordinate entered production operator")
    cd4_banks, cd4 = production.prepare_frozen_cd4(
        CD4_ROOT,
        panel,
        official,
        panel,
        chunk_rows=block_rows,
        max_rss_gib=max_rss_gib,
    )

    reference, reference_audit = production.build_production_effects(
        model,
        panel,
        operator.whitener,
        operator.columns,
        panel_positions,
        pm=pm,
        cd4=cd4,
        alpha=BASELINE_ALPHA,
        weights=weights,
    )
    baseline_pre = production.build_v9_pre_whiten_effects(
        model, panel, alpha=BASELINE_ALPHA, weights=weights
    )
    baseline_pm, pm_audit = production.apply_perturbmatch_k562_replacement(
        baseline_pre,
        model,
        pm,
        panel,
        alpha=BASELINE_ALPHA,
        weights=weights,
    )
    baseline, baseline_audit = finish_prewhiten(
        production, baseline_pm, operator, panel_positions, cd4
    )
    verify_exact("OFF baseline reconstruction", baseline, reference)
    if baseline_audit["prewhiten_sha256"] != reference_audit["candidate_pre_sha256"]:
        raise ValueError("baseline pre-whitening hash disagrees with production")
    if baseline_audit["post_whiten_sha256"] != reference_audit["post_whiten_sha256"]:
        raise ValueError("baseline whitening hash disagrees with production")
    if baseline_audit["final_sha256"] != reference_audit["final_sha256"]:
        raise ValueError("baseline final hash disagrees with production")
    if baseline_audit["prewhiten_sha256"] != EXPECTED_BASELINE_PREWHITEN_SHA256:
        raise ValueError(
            "reconstructed pre-whitening baseline differs from staged a15:\n"
            + json.dumps(
                {
                    "expected": EXPECTED_BASELINE_PREWHITEN_SHA256,
                    "observed": baseline_audit,
                    "production": reference_audit,
                },
                indent=2,
                sort_keys=True,
            )
        )
    expected_a15_scalars = {
        "post_whiten_norm": 4.2292956044543,
        "cd4_fraction": 0.09472706196346517,
        "cd4_preserve_scale": 0.9895210139518232,
    }
    observed_a15_scalars = {
        "post_whiten_norm": float(
            baseline_audit["whitening"]["after_all_median_norm"]
        ),
        "cd4_fraction": float(baseline_audit["cd4"]["fraction"]),
        "cd4_preserve_scale": float(
            baseline_audit["cd4"]["post_add_preserve"][
                "post_add_global_preserve_scale"
            ]
        ),
    }
    for name, expected in expected_a15_scalars.items():
        if not np.isclose(
            observed_a15_scalars[name], expected, rtol=0.0, atol=1e-12
        ):
            raise ValueError(
                f"staged a15 scalar {name} changed: "
                f"{observed_a15_scalars[name]:.17g} != {expected:.17g}"
            )
    staged_metadata = json.loads(STAGED_BASELINE_METADATA_PATH.read_text())
    if (
        staged_metadata.get("final_sha256")
        != EXPECTED_REMOTE_BASELINE_EFFECT_SHA256
        or staged_metadata.get("prewhiten_sha256")
        != EXPECTED_BASELINE_PREWHITEN_SHA256
        or staged_metadata.get("file_sha256")
        != EXPECTED_STAGED_BASELINE_FILE_SHA256
    ):
        raise ValueError("preserved staged-a15 metadata contract changed")
    staged_baseline = np.load(STAGED_BASELINE_PATH, mmap_mode="r")
    if staged_baseline.shape != (PANEL_TARGETS, OFFICIAL_GENES):
        raise ValueError("preserved staged-a15 effect shape changed")
    if sha256_array(staged_baseline) != EXPECTED_REMOTE_BASELINE_EFFECT_SHA256:
        raise ValueError("preserved staged-a15 effect array hash changed")
    baseline_difference = np.asarray(baseline) - np.asarray(staged_baseline)
    baseline_max_abs = float(np.max(np.abs(baseline_difference)))
    baseline_mean_abs = float(np.mean(np.abs(baseline_difference)))
    baseline_allclose = bool(
        np.allclose(baseline, staged_baseline, rtol=1e-12, atol=1e-12)
    )
    if not baseline_allclose:
        raise ValueError(
            "local production reconstruction is not elementwise-close to staged a15: "
            f"max_abs={baseline_max_abs:.17g}, mean_abs={baseline_mean_abs:.17g}"
        )

    pm_supported_targets = set(pm_audit["both_supported_modified"])
    pm_supported = np.asarray(
        [target in pm_supported_targets for target in panel], dtype=bool
    )
    if int(pm_supported.sum()) != EXPECTED_PM_SUPPORTED_TARGETS:
        raise ValueError("PerturbMatch target support changed")
    pm_gene_mask = np.zeros(len(official), dtype=bool)
    pm_gene_mask[np.asarray(pm["target_columns"], dtype=np.int64)] = True
    if int(pm_gene_mask.sum()) != EXPECTED_PM_GENES:
        raise ValueError("PerturbMatch gene support changed")
    if int(baseline_audit["cd4"]["fully_supported_targets"]) != (
        EXPECTED_CD4_COMPLETE_TARGETS
    ):
        raise ValueError("CD4 target support changed")

    loadings, carrier = source_panel_loadings(model, axes, panel)
    source_counts = carrier.sum(axis=0).astype(np.uint8)
    correction = projection_correction(
        model,
        panel,
        axes,
        loadings,
        carrier,
        pm_supported,
        pm_gene_mask,
        alpha=BASELINE_ALPHA,
        pm_blend=float(production.PM_BLEND),
    )
    candidate_pre = apply_projection_correction(
        baseline_pm, correction, enabled=True
    )
    candidate, candidate_audit = finish_prewhiten(
        production, candidate_pre, operator, panel_positions, cd4
    )

    nulls: list[np.ndarray] = []
    null_audits: list[dict[str, object]] = []
    for draw in range(N_NULL):
        permuted = permute_source_loadings(loadings, carrier, draw)
        null_correction = projection_correction(
            model,
            panel,
            axes,
            permuted,
            carrier,
            pm_supported,
            pm_gene_mask,
            alpha=BASELINE_ALPHA,
            pm_blend=float(production.PM_BLEND),
        )
        null_pre = apply_projection_correction(
            baseline_pm, null_correction, enabled=True
        )
        null, null_audit = finish_prewhiten(
            production, null_pre, operator, panel_positions, cd4
        )
        nulls.append(null)
        null_audits.append(
            {
                "draw": draw,
                "prewhiten_sha256": null_audit["prewhiten_sha256"],
                "final_sha256": null_audit["final_sha256"],
                "loading_sha256": sha256_array(permuted),
            }
        )
        memory_checkpoint(f"null {draw}", max_rss_gib)

    clusters, cluster_audit = source_effect_clusters(
        np.asarray(staged_baseline)[:, score_columns]
    )
    reconstruction_checks = {
        "reference_equals_explicit_baseline": np.array_equal(reference, baseline),
        "baseline_prewhiten_hash_matches": (
            baseline_audit["prewhiten_sha256"]
            == reference_audit["candidate_pre_sha256"]
        ),
        "baseline_postwhiten_hash_matches": (
            baseline_audit["post_whiten_sha256"]
            == reference_audit["post_whiten_sha256"]
        ),
        "baseline_final_hash_matches": (
            baseline_audit["final_sha256"] == reference_audit["final_sha256"]
        ),
        "baseline_prewhiten_matches_staged_a15": (
            baseline_audit["prewhiten_sha256"]
            == EXPECTED_BASELINE_PREWHITEN_SHA256
        ),
        "baseline_post_transform_scalars_match_staged_a15": all(
            np.isclose(
                observed_a15_scalars[name],
                expected,
                rtol=0.0,
                atol=1e-12,
            )
            for name, expected in expected_a15_scalars.items()
        ),
        "preserved_staged_a15_array_hash_matches": (
            sha256_array(staged_baseline)
            == EXPECTED_REMOTE_BASELINE_EFFECT_SHA256
        ),
        "local_baseline_elementwise_close_to_staged_a15": baseline_allclose,
        "source_coverage_unchanged": True,
        "panel_absent_from_score_axis": not bool(score_mask[panel_positions].any()),
        "panel_absent_from_operator": not bool(
            np.intersect1d(operator.columns, panel_positions).size
        ),
        "pm_support_unchanged": int(pm_supported.sum())
        == EXPECTED_PM_SUPPORTED_TARGETS,
        "pm_gene_support_unchanged": int(pm_gene_mask.sum()) == EXPECTED_PM_GENES,
        "cd4_support_unchanged": int(
            baseline_audit["cd4"]["fully_supported_targets"]
        )
        == EXPECTED_CD4_COMPLETE_TARGETS,
    }
    if not all(reconstruction_checks.values()):
        raise ValueError(f"reconstruction checks failed: {reconstruction_checks}")

    staging = WORK / f".freeze-{os.getpid()}-{time.time_ns()}"
    audit = {
        "production": {
            "reference": reference_audit,
            "explicit_baseline": baseline_audit,
            "candidate": candidate_audit,
            "pm": pm_audit,
            "staged_a15_remote_final_sha256": (
                EXPECTED_REMOTE_BASELINE_EFFECT_SHA256
            ),
            "local_final_sha256": baseline_audit["final_sha256"],
            "cross_architecture_bitwise_final_match": (
                baseline_audit["final_sha256"]
                == EXPECTED_REMOTE_BASELINE_EFFECT_SHA256
            ),
            "staged_a15_scalars_expected": expected_a15_scalars,
            "staged_a15_scalars_observed": observed_a15_scalars,
            "local_vs_staged_a15": {
                "rtol": 1e-12,
                "atol": 1e-12,
                "max_abs": baseline_max_abs,
                "mean_abs": baseline_mean_abs,
                "allclose": baseline_allclose,
            },
        },
        "source_axes": axis_audit,
        "nulls": null_audits,
        "clusters": cluster_audit,
        "source_count_histogram": {
            str(value): int(np.sum(source_counts == value))
            for value in sorted(set(source_counts.tolist()))
        },
        "score_axis_sha256": sha256_labels(
            official[index] for index in score_columns
        ),
        "operator_columns_sha256": sha256_array(operator.columns, "<i8"),
        "reconstruction_checks": reconstruction_checks,
    }
    manifest = write_freeze_bundle(
        staging,
        official=official,
        panel=panel,
        score_columns=score_columns,
        operator_columns=operator.columns,
        baseline=staged_baseline,
        candidate=candidate,
        nulls=nulls,
        axes=axes,
        loadings=loadings,
        carrier=carrier,
        source_counts=source_counts,
        clusters=clusters,
        pm_supported=pm_supported,
        implementation_hashes=implementation_hashes,
        input_hashes=input_hashes,
        audit=audit,
        state=state,
    )
    staging.replace(FROZEN)
    print(json.dumps(manifest, indent=2))


@dataclass
class HeldTruth:
    values: np.ndarray
    panel_rows: np.ndarray
    official_columns: np.ndarray
    targets: tuple[str, ...]
    score_genes: tuple[str, ...]
    hashes: dict[str, str]


def create_access_marker(manifest: Mapping[str, object]) -> Path:
    marker = WORK / "KOLF_ACCESS_CONSUMED.json"
    record = {
        "experiment": EXPERIMENT,
        "created_unix": time.time(),
        "freeze_manifest_sha256": file_sha256(FROZEN / "manifest.json"),
        "sealed_before_truth": bool(manifest["sealed_before_truth"]),
        "post_truth_fitting_allowed": False,
    }
    payload = (json.dumps(record, indent=2, sort_keys=True) + "\n").encode()
    descriptor = os.open(marker, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        os.write(descriptor, payload)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    return marker


def load_kolf_truth(
    kolf: Path,
    official: Sequence[str],
    panel: Sequence[str],
    score_columns: np.ndarray,
    *,
    state: ProtocolState,
    max_rss_gib: float,
) -> HeldTruth:
    state.open_truth()
    hashes = verify_hash_map(kolf, EXPECTED_KOLF_HASHES, "KOLF held truth")
    targets = read_lines(kolf / "perts.txt")
    genes = read_lines(kolf / "genes.txt")
    panel_set = set(panel)
    expected_targets = tuple(target for target in panel if target in set(targets))
    score_genes = tuple(
        gene for gene in official if gene in set(genes) and gene not in panel_set
    )
    if targets != expected_targets or len(targets) != KOLF_TARGETS:
        raise ValueError("KOLF target gallery changed")
    if sha256_labels(targets) != EXPECTED_KOLF_TARGET_ORDER_SHA256:
        raise ValueError("KOLF target order hash changed")
    if len(score_genes) != KOLF_SCORE_GENES:
        raise ValueError("KOLF score axis length changed")
    if sha256_labels(score_genes) != EXPECTED_KOLF_SCORE_AXIS_SHA256:
        raise ValueError("KOLF score axis hash changed")
    frozen_score_genes = tuple(official[index] for index in score_columns)
    if score_genes != frozen_score_genes:
        raise ValueError("post-freeze KOLF score axis differs from source-derived axis")

    panel_pos = {target: index for index, target in enumerate(panel)}
    official_pos = {gene: index for index, gene in enumerate(official)}
    kolf_pos = {gene: index for index, gene in enumerate(genes)}
    panel_rows = np.asarray([panel_pos[target] for target in targets], dtype=np.int64)
    official_columns = np.asarray(
        [official_pos[gene] for gene in score_genes], dtype=np.int64
    )
    kolf_columns = np.asarray([kolf_pos[gene] for gene in score_genes], dtype=np.int64)
    delta = np.load(kolf / "delta.npy", mmap_mode="r", allow_pickle=False)
    if delta.shape[0] != KOLF_TARGETS or delta.shape[1] != len(genes):
        raise ValueError("KOLF effect shape changed")
    values = np.asarray(delta[np.ix_(np.arange(KOLF_TARGETS), kolf_columns)], dtype=np.float64)
    if not np.isfinite(values).all():
        raise ValueError("KOLF truth contains non-finite values")
    memory_checkpoint("KOLF truth opened", max_rss_gib)
    return HeldTruth(
        values=values,
        panel_rows=panel_rows,
        official_columns=official_columns,
        targets=targets,
        score_genes=score_genes,
        hashes=hashes,
    )


def evaluate_frozen(
    baseline_full: np.ndarray,
    candidate_full: np.ndarray,
    null_full: Sequence[np.ndarray],
    clusters_full: np.ndarray,
    held: HeldTruth,
) -> dict[str, object]:
    rows, columns = held.panel_rows, held.official_columns
    baseline = np.asarray(baseline_full)[np.ix_(rows, columns)]
    candidate = np.asarray(candidate_full)[np.ix_(rows, columns)]
    pds_baseline = pds_per_target(baseline, held.values)
    pds_candidate = pds_per_target(candidate, held.values)
    per_target_delta = pds_candidate - pds_baseline
    null_means = []
    for values in null_full:
        null = np.asarray(values)[np.ix_(rows, columns)]
        null_means.append(
            float((pds_per_target(null, held.values) - pds_baseline).mean())
        )
    clusters = np.asarray(clusters_full, dtype=np.int64)[rows]
    bootstrap = cluster_bootstrap(per_target_delta, clusters)
    observed = float(per_target_delta.mean())
    null_values = np.asarray(null_means, dtype=np.float64)
    null_p95 = float(np.percentile(null_values, 95, method="higher"))
    null_one_sided_p = float(
        (1 + np.sum(null_values >= observed)) / (len(null_values) + 1)
    )
    return {
        "pds_baseline": float(pds_baseline.mean()),
        "pds_candidate": float(pds_candidate.mean()),
        "candidate_minus_baseline": observed,
        "per_target_delta": per_target_delta,
        "pds_baseline_per_target": pds_baseline,
        "pds_candidate_per_target": pds_candidate,
        "null_mean_deltas": null_means,
        "null_p95": null_p95,
        "null_percentile_method": "higher",
        "null_one_sided_p": null_one_sided_p,
        "bootstrap": bootstrap,
    }


def score(kolf: Path, max_rss_gib: float) -> None:
    if (WORK / "results.json").exists():
        raise FileExistsError("experiment-103 result already exists")
    manifest = verify_frozen_bundle()
    if not all(manifest["audit"]["reconstruction_checks"].values()):
        raise ValueError("frozen reconstruction/provenance checks did not all pass")
    if verify_implementation_hashes() != manifest["implementation_hashes"]:
        raise ValueError("implementation changed after freeze")
    create_access_marker(manifest)

    state = ProtocolState(phase="frozen")
    with np.load(FROZEN / "frozen-config.npz", allow_pickle=False) as config:
        official = tuple(config["official"].astype(str))
        panel = tuple(config["panel"].astype(str))
        score_columns = config["score_columns"].astype(np.int64)
        source_counts = config["source_counts"].astype(np.uint8)
        clusters = config["source_clusters"].astype(np.int64)
    held = load_kolf_truth(
        kolf,
        official,
        panel,
        score_columns,
        state=state,
        max_rss_gib=max_rss_gib,
    )
    baseline = np.load(FROZEN / "baseline.npy", mmap_mode="r")
    candidate = np.load(FROZEN / "candidate.npy", mmap_mode="r")
    nulls = [
        np.load(
            FROZEN / "null" / f"source-label-permutation-{draw:03d}.npy",
            mmap_mode="r",
        )
        for draw in range(N_NULL)
    ]
    evaluation = evaluate_frozen(baseline, candidate, nulls, clusters, held)
    coverage_unchanged = bool(
        manifest["audit"]["reconstruction_checks"]["source_coverage_unchanged"]
    )
    bootstrap = evaluation["bootstrap"]
    provenance_passed = (
        state.phase == "truth-open"
        and manifest["sealed_before_truth"] is True
        and manifest["post_truth_fitting_allowed"] is False
        and all(manifest["audit"]["reconstruction_checks"].values())
        and len(held.targets) == KOLF_TARGETS
        and len(held.score_genes) == KOLF_SCORE_GENES
    )
    checks = {
        "candidate_minus_baseline_at_least_0_003": (
            evaluation["candidate_minus_baseline"] >= MINIMUM_DELTA
        ),
        "target_cluster_bootstrap_lower_95_above_zero": (
            bool(bootstrap["valid"])
            and bootstrap["lower_95"] is not None
            and float(bootstrap["lower_95"]) > 0
        ),
        "observed_above_source_label_null_p95": (
            evaluation["candidate_minus_baseline"] > evaluation["null_p95"]
        ),
        "source_label_null_one_sided_p_at_most_0_05": (
            evaluation["null_one_sided_p"] <= 0.05
        ),
        "source_coverage_unchanged": coverage_unchanged,
        "all_provenance_and_reconstruction_checks_pass": provenance_passed,
    }
    result = {
        "experiment": EXPERIMENT,
        "status": "GO" if all(checks.values()) else "NO-GO",
        "interpretation": (
            "Evidence applies only to the donor-axis residual surrogate, not official "
            "COMPASS, whose target axis requires unavailable target perturbations. "
            "This is effect-level evidence only and does not authorize count generation "
            "or challenge submission."
        ),
        "metric": "complete-gallery raw uncentered PDS",
        "targets": len(held.targets),
        "score_genes": len(held.score_genes),
        "scores": {
            key: value
            for key, value in evaluation.items()
            if key
            not in {
                "per_target_delta",
                "pds_baseline_per_target",
                "pds_candidate_per_target",
            }
        },
        "checks": checks,
        "kolf_hashes": held.hashes,
        "freeze_manifest_sha256": file_sha256(FROZEN / "manifest.json"),
        "source_count_histogram_on_kolf": {
            str(value): int(np.sum(source_counts[held.panel_rows] == value))
            for value in sorted(set(source_counts[held.panel_rows].tolist()))
        },
        "post_truth_fitting_performed": False,
        "peak_rss_gib": peak_rss_gib(),
    }
    atomic_save_npz(
        WORK / "per-target-results.npz",
        targets=np.asarray(held.targets),
        pds_baseline=evaluation["pds_baseline_per_target"],
        pds_candidate=evaluation["pds_candidate_per_target"],
        delta=evaluation["per_target_delta"],
        clusters=clusters[held.panel_rows],
    )
    atomic_write_text(
        WORK / "results.json", json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(result, indent=2))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    freeze_parser = subparsers.add_parser("freeze")
    freeze_parser.add_argument("--max-rss-gib", type=float, default=24.0)
    freeze_parser.add_argument("--block-rows", type=int, default=128)
    score_parser = subparsers.add_parser("score")
    score_parser.add_argument("--kolf", type=Path, required=True)
    score_parser.add_argument("--max-rss-gib", type=float, default=24.0)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    limit = validate_memory_limit(args.max_rss_gib)
    with experiment_lock():
        if args.command == "freeze":
            if args.block_rows < 1:
                raise ValueError("--block-rows must be positive")
            freeze(limit, args.block_rows)
        else:
            score(args.kolf.expanduser().resolve(), limit)


if __name__ == "__main__":
    main()

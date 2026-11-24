#!/usr/bin/env python3
"""079 — untouched-KOLF confirmation of frozen reciprocal local centering.

This experiment independently tests experiment 075's already-selected
reciprocal-kNN local query-centering transform:

    k = 64
    beta = 0.5
    minimum reciprocal degree = 2, then deterministic top-k fallback

There is deliberately no candidate grid and no KOLF-facing hyperparameter.
The exact production-v9 effect query is built for all 300 VCC targets from
K562-GWPS, Orion HCT116, and Orion HEK293T, then transformed by the exact
three-donor union-axis gamma=0.65 whitener.  All 300 panel-gene coordinates
are absent from the whitening operator and from local centering.

Both identity and centered queries are frozen over the complete 300-query
panel before the compact KOLF2.1J artifact is opened.  KOLF is used only as
raw-centered held truth.  It contributes no source effect, covariance row,
query neighbor, normalization statistic, or hyperparameter decision.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import os
import platform
import resource
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Mapping, Sequence, TypeVar

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "4")

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from vcc2026.donor_transfer import Donor  # noqa: E402
from vcc2026.multi_donor import MultiDonorModel  # noqa: E402
from vcc2026.whiten import Whitener  # noqa: E402


EXPERIMENT = "079-kolf-local-centering"
PRODUCTION_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
ALPHA = 1.0796
MODEL_BETA = 0.0
WHITEN_GAMMA = 0.65
WHITEN_LAMBDA = 0.0
LOCAL_K = 64
LOCAL_BETA = 0.5
MINIMUM_MUTUAL = 2
BOOTSTRAP_REPLICATES = 20_000
BOOTSTRAP_SEED = 7902026
MINIMUM_RAW_PDS_DELTA = 0.0015

EXPECTED_PANEL_TARGETS = 300
EXPECTED_OFF_PANEL_GENES = 18_233
EXPECTED_OPERATOR_GENES = 17_942
EXPECTED_KOLF_TARGETS = 273
EXPECTED_KOLF_SCORE_GENES = 17_806

EXPECTED_PRODUCTION_WHITENER_SHA256 = (
    "5608e807230dea04f93f5e804ed26e8f7285d9487d16a2791283793ffb80fed3"
)
EXPECTED_PRODUCTION_WHITENER_JSON_SHA256 = (
    "a77f9425ac524cdcf3a53a1744520e1e63cb8b925ef77fae2864c4d8b72337b6"
)
EXPECTED_EXP075_RESULTS_SHA256 = (
    "72968b6ca3518e2c49c83365b661ce65942f50268c12587a04aa30fc41abbb4e"
)
EXPECTED_EXP075_SELECTION_FREEZE_SHA256 = (
    "d0aabafee320a8f9cd37aab2a1c2ecc88794b8be019ed3a865886173219b2988"
)
EXPECTED_PREREGISTRATION_SHA256 = (
    "4cc86fbdca015200b1d1f56aa526ccfb098fa8bcf25e4c878390f7a6fe2b464d"
)

EXPECTED_CHALLENGE_HASHES = {
    "gene_names.csv": "25bfa66715e186bebabce7ac788bbcea47e2bf59ca70be1f8f3a06f2f0e47201",
    "pert_counts.csv": "f57edd7b912ebd718efc7ee9d0f334772513e7cc418d133ce525470e373b3276",
}
EXPECTED_IMPLEMENTATION_HASHES = {
    "src/vcc2026/donor_transfer.py": (
        "cecb27035a6c301f96e7fe364032216b1211bb2a92fddd1f9439d9cf8730cb3e"
    ),
    "src/vcc2026/multi_donor.py": (
        "3bb75320ed54e438b0a1ea10acff413e9738fc3294056db06ac337ed4d51cd8f"
    ),
    "src/vcc2026/whiten.py": (
        "5dbddefde17192cda4bc7f8f7eceb202faf672de1833c42251366078e7e18bc1"
    ),
}
EXPECTED_DONOR_HASHES: dict[str, dict[str, str]] = {
    "k562_gwps": {
        "delta.npy": "0c9422280fed7afb81ca474e8b9024af788c9ffbde32e11f96cb52422e296aa9",
        "b_ctrl.npy": "d06cf51d51cde1ae999e86bbd3be0bdc15471392bae153a1aa1e21dc5c81b957",
        "genes.txt": "f499ffc3077e4f06c86cd4504bfeefd34741bc2fd515be93f7e436c5c08fa2ba",
        "perts.txt": "b784821e4a662aa8ab8d55da357fc712f109dda5a464c3e7e047b15b61408bec",
        "meta.json": "0bd3478aa6b63164bd4816cd8bbb439b805b6b93b2a267ddf8114b196c5e19aa",
    },
    "orion_hct116": {
        "delta.npy": "72c7a01dc2869bd3e94c37e89b158c3aaf4048a0dff63da6542afa2de4ec2b33",
        "b_ctrl.npy": "b3044d439cfac932024e0a48152196387714e4b460a339fb9a03e32fcb2fe069",
        "genes.txt": "db8ba814bd30dcc7fe48be659f4dee3f6955b626d7d7e385f4630232733f893e",
        "perts.txt": "fd2fbbbebdd43c360cd274a246c1eeeb9a2139711ee12a475a8e116021c5caa8",
        "meta.json": "7a1ec5a9323c6ebd4b586aa300b28bcd9ddf271e67f04bdaf353cb15b8b2d6c8",
    },
    "orion_hek293t": {
        "delta.npy": "d95202d984d9322dc3830f23b6ef55001ba8562dea1d72c43fb0c9d5064c7d1b",
        "b_ctrl.npy": "e04234cffbe102731cc562e8fb06bddb8ef0644b2af894b8a0ee6ce079ca4984",
        "genes.txt": "db8ba814bd30dcc7fe48be659f4dee3f6955b626d7d7e385f4630232733f893e",
        "perts.txt": "3846edf57cc7a0bd91a499ce9be4ed68c9c8ccc1f36ba062c3bb1301695ef547",
        "meta.json": "a7310e69e2e09ad42d88f9192a938b14bf3db44be92575123b86508f0818710f",
    },
}
EXPECTED_KOLF_HASHES = {
    "delta.npy": "a8ed5c7eebb0b01a6a0a4effb0651fd0c6c706aee81cac87b6369d4d04c91d6b",
    "b_ctrl.npy": "dd4c30afee5b6bd7fc68c469da30a1eb771e367d186de41d13de16ca4aee7fc9",
    "genes.txt": "2f077fda94cc07cf6cfddae28fc78440e80963e4e866f9904b91417f3189420b",
    "perts.txt": "b366592072d1e0da04aaabafbcc78882d34c97ff26eb2fca5fd9a9af0084f22f",
    "meta.json": "57737834c91663d154e84721ce73c0a33f4b670fea1d6bdab8ad35d0cf8b466a",
}

EXPECTED_PANEL_ORDER_SHA256 = (
    "cd3d00f1cd5da660609ad41d85b53e2d5557f7f5edb3e1f64559f78f7fb3cbc4"
)
EXPECTED_KOLF_TARGET_ORDER_SHA256 = (
    "b366592072d1e0da04aaabafbcc78882d34c97ff26eb2fca5fd9a9af0084f22f"
)
EXPECTED_KOLF_SCORE_AXIS_SHA256 = (
    "280aab34a0e86f27db8ade213e000f9066c8d5f9b89edd3ea966aba061ba9ae9"
)

T = TypeVar("T")


def peak_rss_gib() -> float:
    """Return process peak RSS in GiB on macOS and Linux."""
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


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


def sha256_labels(values: Iterable[str]) -> str:
    return hashlib.sha256(("\n".join(values) + "\n").encode()).hexdigest()


def sha256_array(values: np.ndarray) -> str:
    array = np.ascontiguousarray(values, dtype="<f8")
    return hashlib.sha256(array.tobytes()).hexdigest()


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: required column {column!r} is missing")
        values = tuple(
            str(row[column]).strip()
            for row in reader
            if row.get(column) and str(row[column]).strip()
        )
    if not values or len(values) != len(set(values)):
        raise ValueError(f"{path}: {column!r} is empty or contains duplicates")
    return values


def read_lines(path: Path) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if not values or len(values) != len(set(values)):
        raise ValueError(f"{path}: axis is empty or contains duplicates")
    return values


def verify_hash(path: Path, expected: str, label: str) -> str:
    actual = sha256_file(path)
    if actual != expected:
        raise ValueError(f"{label}: SHA-256 {actual} != frozen {expected}")
    return actual


def verify_hash_map(root: Path, expected: Mapping[str, str], label: str) -> None:
    for filename, digest in expected.items():
        verify_hash(root / filename, digest, f"{label}/{filename}")


def validate_exp075_freeze(path: Path) -> dict[str, object]:
    verify_hash(path, EXPECTED_EXP075_RESULTS_SHA256, "experiment 075 result")
    result = json.loads(path.read_text())
    selected = result["inner"]["selection_freeze"]["deployment_selected"]
    observed = {
        "experiment": result.get("experiment"),
        "selected": selected,
        "selection_freeze_sha256": result["inner"].get("selection_freeze_sha256"),
        "query_only": result["method"].get("query_only"),
        "hidden_gallery_used": result["method"].get("hidden_gallery_used"),
        "minimum_mutual": result["method"].get(
            "minimum_mutual_neighbors_before_knn_fallback"
        ),
    }
    expected = {
        "experiment": "075-postwhiten-local-centering",
        "selected": "reciprocal-local-center:k=64:beta=0.5",
        "selection_freeze_sha256": EXPECTED_EXP075_SELECTION_FREEZE_SHA256,
        "query_only": True,
        "hidden_gallery_used": False,
        "minimum_mutual": MINIMUM_MUTUAL,
    }
    if observed != expected:
        raise ValueError(
            "experiment 075 freeze contract drifted:\n"
            + json.dumps({"expected": expected, "observed": observed}, indent=2)
        )
    return observed


def validate_preregistration(path: Path) -> dict[str, object]:
    digest = verify_hash(
        path,
        EXPECTED_PREREGISTRATION_SHA256,
        "experiment 079 preregistration",
    )
    protocol = json.loads(path.read_text())
    observed = {
        "candidate": protocol.get("candidate"),
        "candidate_grid": protocol.get("candidate_grid"),
        "retuning_on_kolf": protocol.get("retuning_on_kolf"),
        "bootstrap": protocol.get("bootstrap"),
        "promotion_gate": protocol.get("promotion_gate"),
    }
    expected = {
        "candidate": {
            "k": LOCAL_K,
            "beta": LOCAL_BETA,
            "minimum_mutual": MINIMUM_MUTUAL,
        },
        "candidate_grid": None,
        "retuning_on_kolf": False,
        "bootstrap": {
            "unit": "target_gene",
            "replicates": BOOTSTRAP_REPLICATES,
            "seed": BOOTSTRAP_SEED,
            "interval": "percentile_95",
        },
        "promotion_gate": {
            "raw_pds_delta_positive": True,
            "lower_95_positive": True,
            "minimum_raw_pds_delta": MINIMUM_RAW_PDS_DELTA,
            "all_required": True,
        },
    }
    if observed != expected:
        raise ValueError(
            "preregistered protocol does not match executable constants:\n"
            + json.dumps({"expected": expected, "observed": observed}, indent=2)
        )
    return {"sha256": digest, **observed}


def row_normalize(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    norms = np.linalg.norm(array, axis=1, keepdims=True)
    return np.divide(
        array,
        np.maximum(norms, 1e-12),
        out=np.zeros_like(array),
    )


@dataclass(frozen=True)
class LocalCenterDiagnostics:
    rows: int
    genes: int
    k: int
    beta: float
    minimum_mutual: int
    fallback_rows: int
    mean_mutual_degree: float
    median_mutual_degree: float
    median_centroid_norm: float
    max_row_norm_error: float


def reciprocal_local_center(
    query: np.ndarray,
    *,
    k: int = LOCAL_K,
    beta: float = LOCAL_BETA,
    minimum_mutual: int = MINIMUM_MUTUAL,
) -> tuple[np.ndarray, LocalCenterDiagnostics]:
    """Independent implementation of experiment 075's frozen transform."""
    array = np.asarray(query, dtype=np.float64)
    if array.ndim != 2 or array.shape[0] < 3:
        raise ValueError("query must be a 2-D matrix with at least three rows")
    if k < 1 or k >= array.shape[0]:
        raise ValueError(f"k must be in [1, {array.shape[0] - 1}]")
    if not (0.0 < beta <= 1.0):
        raise ValueError("beta must be in (0, 1]")
    if minimum_mutual < 1:
        raise ValueError("minimum_mutual must be positive")
    if not np.isfinite(array).all():
        raise ValueError("query contains non-finite values")

    row_norms = np.linalg.norm(array, axis=1)
    unit = np.divide(
        array,
        np.maximum(row_norms[:, None], 1e-12),
        out=np.zeros_like(array),
    )
    similarities = unit @ unit.T
    np.fill_diagonal(similarities, -np.inf)
    neighbors = np.argsort(-similarities, axis=1, kind="stable")[:, :k]
    directed = np.zeros((len(array), len(array)), dtype=bool)
    directed[np.arange(len(array))[:, None], neighbors] = True
    reciprocal = directed & directed.T

    centroid = np.empty_like(unit)
    degrees = reciprocal.sum(axis=1).astype(np.int64)
    fallback_rows = 0
    for row in range(len(array)):
        local = np.flatnonzero(reciprocal[row])
        if len(local) < minimum_mutual:
            local = neighbors[row]
            fallback_rows += 1
        if row in local:
            raise AssertionError("self-neighbor entered reciprocal centroid")
        centroid[row] = unit[local].mean(axis=0)

    adjusted = unit - beta * centroid
    adjusted_norms = np.linalg.norm(adjusted, axis=1, keepdims=True)
    if np.any(adjusted_norms <= 1e-12):
        raise ValueError("local centering collapsed a nonzero query direction")
    transformed = adjusted / adjusted_norms * row_norms[:, None]
    transformed[row_norms <= 1e-12] = 0.0
    diagnostics = LocalCenterDiagnostics(
        rows=int(array.shape[0]),
        genes=int(array.shape[1]),
        k=int(k),
        beta=float(beta),
        minimum_mutual=int(minimum_mutual),
        fallback_rows=int(fallback_rows),
        mean_mutual_degree=float(np.mean(degrees)),
        median_mutual_degree=float(np.median(degrees)),
        median_centroid_norm=float(np.median(np.linalg.norm(centroid, axis=1))),
        max_row_norm_error=float(
            np.max(np.abs(np.linalg.norm(transformed, axis=1) - row_norms))
        ),
    )
    return transformed, diagnostics


def production_operator_columns(
    model: MultiDonorModel, panel: Sequence[str]
) -> np.ndarray:
    union: set[str] = set()
    for name in model.names:
        union |= set(map(str, model.target_genes[model.tcols[name]]))
    panel_set = set(map(str, panel))
    union -= panel_set
    target_index = {str(gene): index for index, gene in enumerate(model.target_genes)}
    columns = np.asarray(
        sorted(target_index[gene] for gene in union),
        dtype=np.int64,
    )
    if panel_set & set(map(str, model.target_genes[columns])):
        raise AssertionError("a VCC target coordinate entered the production operator")
    return columns


def apply_production_whitener(
    effects: np.ndarray,
    columns: np.ndarray,
    operator_path: Path,
    metadata_path: Path,
) -> tuple[np.ndarray, dict[str, object]]:
    operator_hash = verify_hash(
        operator_path,
        EXPECTED_PRODUCTION_WHITENER_SHA256,
        "production whitener",
    )
    metadata_hash = verify_hash(
        metadata_path,
        EXPECTED_PRODUCTION_WHITENER_JSON_SHA256,
        "production whitener metadata",
    )
    metadata = json.loads(metadata_path.read_text())
    expected_metadata = {
        "n_cols": EXPECTED_OPERATOR_GENES,
        "basis": [33_787, EXPECTED_OPERATOR_GENES],
        "donors": list(PRODUCTION_DONORS),
    }
    if metadata != expected_metadata:
        raise ValueError(
            "production whitener metadata changed:\n"
            + json.dumps(
                {"expected": expected_metadata, "observed": metadata},
                indent=2,
            )
        )
    if len(columns) != EXPECTED_OPERATOR_GENES:
        raise ValueError(
            f"operator axis has {len(columns)} genes, expected {EXPECTED_OPERATOR_GENES}"
        )

    output = np.asarray(effects, dtype=np.float64).copy()
    before_all = float(np.median(np.linalg.norm(output, axis=1)))
    before_operator = (
        float(np.median(np.linalg.norm(output[:, columns], axis=1))) or 1.0
    )
    whitener = Whitener.load(operator_path)
    if whitener.evecs.shape != (len(columns), len(columns)):
        raise ValueError(
            f"operator shape {whitener.evecs.shape} != {(len(columns), len(columns))}"
        )
    transformed = whitener.apply(
        output[:, columns],
        WHITEN_GAMMA,
        WHITEN_LAMBDA,
        renorm=False,
    )
    after_operator_unscaled = (
        float(np.median(np.linalg.norm(transformed, axis=1))) or 1.0
    )
    preserve_scale = before_operator / after_operator_unscaled
    output[:, columns] = transformed * preserve_scale
    del transformed, whitener
    gc.collect()
    if not np.isfinite(output).all():
        raise FloatingPointError("production-whitened query contains non-finite values")
    return output, {
        "sha256": operator_hash,
        "metadata_sha256": metadata_hash,
        "metadata": metadata,
        "columns": len(columns),
        "columns_index_sha256": hashlib.sha256(
            np.ascontiguousarray(columns, dtype="<i8").tobytes()
        ).hexdigest(),
        "gamma": WHITEN_GAMMA,
        "lambda": WHITEN_LAMBDA,
        "renormalization": "candidate-specific on-axis preserve",
        "before_operator_median_norm": before_operator,
        "after_operator_unscaled_median_norm": after_operator_unscaled,
        "operator_preserve_scale": preserve_scale,
        "before_all_median_norm": before_all,
        "after_all_median_norm": float(np.median(np.linalg.norm(output, axis=1))),
    }


def verify_production_inputs(
    donor_root: Path,
    genes_csv: Path,
    panel_csv: Path,
) -> None:
    verify_hash(
        genes_csv,
        EXPECTED_CHALLENGE_HASHES["gene_names.csv"],
        "challenge/gene_names.csv",
    )
    verify_hash(
        panel_csv,
        EXPECTED_CHALLENGE_HASHES["pert_counts.csv"],
        "challenge/pert_counts.csv",
    )
    for donor, hashes in EXPECTED_DONOR_HASHES.items():
        verify_hash_map(donor_root / donor, hashes, f"donor/{donor}")
    for relative, digest in EXPECTED_IMPLEMENTATION_HASHES.items():
        verify_hash(ROOT / relative, digest, f"implementation/{relative}")


def build_production_v9_query(
    donor_root: Path,
    official_genes: Sequence[str],
    panel: Sequence[str],
    operator_path: Path,
    operator_metadata_path: Path,
    *,
    max_rss_gib: float,
) -> tuple[np.ndarray, dict[str, object]]:
    """Build the exact v9 pre-gate effect query for all 300 panel targets."""
    donors = {name: Donor.load(donor_root / name) for name in PRODUCTION_DONORS}
    model = MultiDonorModel(
        donors,
        np.asarray(official_genes),
        gene_norm=False,
        eb=False,
    )
    if tuple(model.names) != PRODUCTION_DONORS:
        raise ValueError(
            f"production donor order {tuple(model.names)} != {PRODUCTION_DONORS}"
        )
    weights = {name: 1.0 for name in model.names}
    coverage = model.coverage(panel)
    pre_whiten = np.vstack(
        [
            model.effect(
                target,
                ALPHA,
                MODEL_BETA,
                weights=weights,
            )
            for target in panel
        ]
    ).astype(np.float64)
    if pre_whiten.shape != (EXPECTED_PANEL_TARGETS, len(official_genes)):
        raise ValueError(f"unexpected production query shape {pre_whiten.shape}")
    memory_checkpoint("three-donor pre-whitening query", max_rss_gib)

    operator_columns = production_operator_columns(model, panel)
    panel_positions = np.asarray(
        [
            {gene: index for index, gene in enumerate(official_genes)}[target]
            for target in panel
        ],
        dtype=np.int64,
    )
    if np.intersect1d(operator_columns, panel_positions).size:
        raise AssertionError(
            "one or more of the 300 panel coordinates entered whitening"
        )
    query, whitening_audit = apply_production_whitener(
        pre_whiten,
        operator_columns,
        operator_path,
        operator_metadata_path,
    )
    audit = {
        "donors": list(PRODUCTION_DONORS),
        "weights": weights,
        "gene_norm": False,
        "empirical_bayes": False,
        "alpha": ALPHA,
        "beta": MODEL_BETA,
        "targets": len(panel),
        "genes": len(official_genes),
        "coverage": coverage,
        "pre_whiten_sha256": sha256_array(pre_whiten),
        "post_whiten_sha256": sha256_array(query),
        "panel_coordinates_excluded_from_operator": int(
            len(panel_positions)
            - np.intersect1d(operator_columns, panel_positions).size
        ),
        "operator": whitening_audit,
    }
    del pre_whiten, model, donors
    gc.collect()
    memory_checkpoint("production query frozen", max_rss_gib)
    return query, audit


def freeze_identity_and_candidate(
    production_query: np.ndarray,
    official_genes: Sequence[str],
    panel: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, object]]:
    """Freeze both query arms before any held-truth artifact is opened."""
    panel_set = set(panel)
    score_positions = np.asarray(
        [index for index, gene in enumerate(official_genes) if gene not in panel_set],
        dtype=np.int64,
    )
    if len(score_positions) != EXPECTED_OFF_PANEL_GENES:
        raise ValueError(
            f"off-panel production axis has {len(score_positions)} genes, "
            f"expected {EXPECTED_OFF_PANEL_GENES}"
        )
    if panel_set & {official_genes[index] for index in score_positions}:
        raise AssertionError("a panel coordinate entered local centering")

    identity = np.asarray(production_query[:, score_positions], dtype=np.float64).copy()
    centered, diagnostics = reciprocal_local_center(identity)
    audit = {
        "formula": (
            "q_i' = ||q_i|| * normalize(normalize(q_i) - 0.5 * "
            "mean(reciprocal_64NN(normalized full-panel queries)))"
        ),
        "query_rows_seen_by_transform": int(len(identity)),
        "query_axis_genes": int(identity.shape[1]),
        "panel_coordinates_excluded": len(panel),
        "panel_coordinate_overlap": 0,
        "target_labels_used": False,
        "held_truth_used": False,
        "row_norm_preserved": diagnostics.max_row_norm_error < 1e-9,
        "diagnostics": diagnostics.__dict__,
        "identity_sha256": sha256_array(identity),
        "centered_sha256": sha256_array(centered),
        "score_positions_sha256": hashlib.sha256(
            np.ascontiguousarray(score_positions, dtype="<i8").tobytes()
        ).hexdigest(),
    }
    return identity, centered, score_positions, audit


@dataclass(frozen=True)
class HeldTruth:
    values: np.ndarray
    targets: tuple[str, ...]
    score_genes: tuple[str, ...]
    panel_rows: np.ndarray
    production_score_columns: np.ndarray
    audit: dict[str, object]


def load_kolf_held_truth(
    kolf_root: Path,
    official_genes: Sequence[str],
    panel: Sequence[str],
    production_score_positions: np.ndarray,
    *,
    max_rss_gib: float,
) -> HeldTruth:
    """Open KOLF only after query freeze; return raw-centered held truth."""
    verify_hash_map(kolf_root, EXPECTED_KOLF_HASHES, "held-truth/kolf2")
    metadata = json.loads((kolf_root / "meta.json").read_text())
    expected_marker = "delta = int8 source value * 0.04"
    if expected_marker not in str(metadata.get("effect_definition", "")):
        raise ValueError("KOLF effect definition is not the frozen quantized z-score")

    kolf_targets = read_lines(kolf_root / "perts.txt")
    kolf_genes = read_lines(kolf_root / "genes.txt")
    panel_set = set(panel)
    shared_targets = tuple(target for target in panel if target in set(kolf_targets))
    score_genes = tuple(
        gene
        for gene in official_genes
        if gene in set(kolf_genes) and gene not in panel_set
    )
    if len(shared_targets) != EXPECTED_KOLF_TARGETS:
        raise ValueError(
            f"KOLF/VCC target overlap {len(shared_targets)} != {EXPECTED_KOLF_TARGETS}"
        )
    if len(score_genes) != EXPECTED_KOLF_SCORE_GENES:
        raise ValueError(
            f"KOLF score axis {len(score_genes)} != {EXPECTED_KOLF_SCORE_GENES}"
        )
    if sha256_labels(panel) != EXPECTED_PANEL_ORDER_SHA256:
        raise ValueError("official panel order changed")
    if sha256_labels(shared_targets) != EXPECTED_KOLF_TARGET_ORDER_SHA256:
        raise ValueError("KOLF/VCC shared-target order changed")
    if sha256_labels(score_genes) != EXPECTED_KOLF_SCORE_AXIS_SHA256:
        raise ValueError("KOLF/VCC score-gene axis changed")
    if tuple(kolf_targets) != shared_targets:
        raise ValueError(
            "compact KOLF row order is not the official-panel subset order"
        )
    if panel_set & set(score_genes):
        raise AssertionError("a VCC panel coordinate entered held-truth scoring")

    panel_pos = {target: index for index, target in enumerate(panel)}
    score_pos = {
        official_genes[index]: local
        for local, index in enumerate(production_score_positions)
    }
    panel_rows = np.asarray(
        [panel_pos[target] for target in shared_targets], dtype=np.int64
    )
    production_score_columns = np.asarray(
        [score_pos[gene] for gene in score_genes],
        dtype=np.int64,
    )
    kolf_gene_pos = {gene: index for index, gene in enumerate(kolf_genes)}
    kolf_columns = np.asarray(
        [kolf_gene_pos[gene] for gene in score_genes],
        dtype=np.int64,
    )

    delta = np.load(kolf_root / "delta.npy", mmap_mode="r", allow_pickle=False)
    if delta.shape != (EXPECTED_KOLF_TARGETS, len(kolf_genes)):
        raise ValueError(f"unexpected KOLF delta shape {delta.shape}")
    raw = np.asarray(delta[:, kolf_columns], dtype=np.float64)
    values = raw - raw.mean(axis=0, keepdims=True)
    mmap = getattr(delta, "_mmap", None)
    if mmap is not None:
        mmap.close()
    del delta, raw
    gc.collect()
    if not np.isfinite(values).all():
        raise ValueError("KOLF held truth contains non-finite values")
    memory_checkpoint("KOLF held truth opened", max_rss_gib)
    return HeldTruth(
        values=values,
        targets=shared_targets,
        score_genes=score_genes,
        panel_rows=panel_rows,
        production_score_columns=production_score_columns,
        audit={
            "role": "untouched held truth only",
            "representation": (
                "decoded raw quantized z-score, per-gene centered across all "
                "273 compact KOLF VCC rows; no scaling or whitening"
            ),
            "source_repository": metadata.get("source_repository"),
            "source_commit": metadata.get("source_commit"),
            "targets": len(shared_targets),
            "genes": len(score_genes),
            "target_order_sha256": sha256_labels(shared_targets),
            "score_gene_order_sha256": sha256_labels(score_genes),
            "rows_used_as_source": 0,
            "rows_used_in_operator": 0,
            "rows_used_in_query_neighbors": 0,
            "used_for_hyperparameter_selection": False,
            "b_ctrl_used": False,
        },
    )


def freeze_queries_then_open_truth(
    query_builder: Callable[[], tuple[np.ndarray, dict[str, object]]],
    truth_loader: Callable[[np.ndarray], T],
    official_genes: Sequence[str],
    panel: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, object], T]:
    """Make the holdout-opening order executable and directly testable."""
    production_query, production_audit = query_builder()
    identity, centered, score_positions, transform_audit = (
        freeze_identity_and_candidate(production_query, official_genes, panel)
    )
    del production_query
    gc.collect()
    truth = truth_loader(score_positions)
    return (
        identity,
        centered,
        score_positions,
        {
            "production": production_audit,
            "transform": transform_audit,
        },
        truth,
    )


def pds_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    if prediction.shape != truth.shape:
        raise ValueError(f"PDS shape mismatch {prediction.shape} != {truth.shape}")
    if prediction.ndim != 2 or prediction.shape[0] < 2:
        raise ValueError(
            "PDS requires a two-dimensional gallery with at least two rows"
        )
    scores = row_normalize(prediction) @ row_normalize(truth).T
    if not np.isfinite(scores).all():
        raise FloatingPointError("non-finite cosine score")
    own = np.diag(scores)[:, None]
    ranks = (scores > own).sum(axis=1).astype(np.float64)
    ranks += 0.5 * ((scores == own).sum(axis=1) - 1)
    return 1.0 - ranks / (len(scores) - 1)


def target_cluster_bootstrap(
    deltas: np.ndarray,
    *,
    replicates: int = BOOTSTRAP_REPLICATES,
    seed: int = BOOTSTRAP_SEED,
    chunk: int = 1000,
) -> dict[str, float | int]:
    values = np.asarray(deltas, dtype=np.float64)
    if values.ndim != 1 or len(values) < 2:
        raise ValueError(
            "target-cluster bootstrap requires a one-dimensional target vector"
        )
    rng = np.random.default_rng(seed)
    draws = np.empty(replicates, dtype=np.float64)
    for lo in range(0, replicates, chunk):
        hi = min(lo + chunk, replicates)
        indices = rng.integers(0, len(values), size=(hi - lo, len(values)))
        draws[lo:hi] = values[indices].mean(axis=1)
    return {
        "unit": "target_gene",
        "replicates": int(replicates),
        "seed": int(seed),
        "mean": float(np.mean(draws)),
        "lower_95": float(np.quantile(draws, 0.025)),
        "upper_95": float(np.quantile(draws, 0.975)),
        "probability_positive": float(np.mean(draws > 0.0)),
    }


def promotion_gate(
    delta: float, bootstrap: Mapping[str, float | int]
) -> dict[str, object]:
    gate: dict[str, object] = {
        "predeclared": True,
        "raw_pds_delta_positive": delta > 0.0,
        "target_cluster_bootstrap_lower_95_positive": (
            float(bootstrap["lower_95"]) > 0.0
        ),
        "raw_pds_delta_at_least_0_0015": delta >= MINIMUM_RAW_PDS_DELTA,
        "minimum_raw_pds_delta": MINIMUM_RAW_PDS_DELTA,
    }
    gate["passed"] = all(
        bool(value)
        for key, value in gate.items()
        if key not in {"predeclared", "minimum_raw_pds_delta"}
    )
    return gate


def ensure_output_directory(path: Path) -> Path:
    output = path.expanduser().resolve()
    if HERE.resolve() not in (output, *output.parents):
        raise ValueError("output directory must stay inside experiment 079")
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(
            f"refusing to overwrite nonempty output directory {output}"
        )
    output.mkdir(parents=True, exist_ok=True)
    return output


def write_per_target(
    path: Path,
    targets: Sequence[str],
    identity: np.ndarray,
    centered: np.ndarray,
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "target_gene",
                "identity_pds",
                "local_center_pds",
                "delta",
            ),
        )
        writer.writeheader()
        for target, baseline, candidate in zip(
            targets, identity, centered, strict=True
        ):
            writer.writerow(
                {
                    "target_gene": target,
                    "identity_pds": float(baseline),
                    "local_center_pds": float(candidate),
                    "delta": float(candidate - baseline),
                }
            )


def run_experiment(args: argparse.Namespace) -> dict[str, object]:
    output = ensure_output_directory(args.output_dir)
    genes_csv = args.genes_csv.expanduser().resolve()
    panel_csv = args.panel_csv.expanduser().resolve()
    donor_root = args.donor_root.expanduser().resolve()
    kolf_root = args.kolf.expanduser().resolve()
    operator_path = args.whitener.expanduser().resolve()
    operator_metadata_path = args.whitener_metadata.expanduser().resolve()
    exp075_results = args.exp075_results.expanduser().resolve()

    preregistration = validate_preregistration(HERE / "preregistered_protocol.json")
    frozen075 = validate_exp075_freeze(exp075_results)
    verify_production_inputs(donor_root, genes_csv, panel_csv)
    official_genes = read_csv_column(genes_csv, "gene_name")
    panel = read_csv_column(panel_csv, "target_gene")
    if len(panel) != EXPECTED_PANEL_TARGETS:
        raise ValueError(f"expected 300 panel targets, got {len(panel)}")
    if not set(panel) <= set(official_genes):
        raise ValueError(
            "one or more VCC targets are absent from the official gene axis"
        )
    memory_checkpoint("production input verification", args.max_rss_gib)

    def query_builder() -> tuple[np.ndarray, dict[str, object]]:
        return build_production_v9_query(
            donor_root,
            official_genes,
            panel,
            operator_path,
            operator_metadata_path,
            max_rss_gib=args.max_rss_gib,
        )

    def truth_loader(score_positions: np.ndarray) -> HeldTruth:
        return load_kolf_held_truth(
            kolf_root,
            official_genes,
            panel,
            score_positions,
            max_rss_gib=args.max_rss_gib,
        )

    identity, centered, _, query_audit, held = freeze_queries_then_open_truth(
        query_builder,
        truth_loader,
        official_genes,
        panel,
    )
    identity_eval = identity[np.ix_(held.panel_rows, held.production_score_columns)]
    centered_eval = centered[np.ix_(held.panel_rows, held.production_score_columns)]
    if (
        identity_eval.shape != held.values.shape
        or centered_eval.shape != held.values.shape
    ):
        raise ValueError(
            "frozen-query/KOLF shape mismatch: "
            f"{identity_eval.shape}, {centered_eval.shape}, {held.values.shape}"
        )

    identity_per = pds_per_target(identity_eval, held.values)
    centered_per = pds_per_target(centered_eval, held.values)
    per_target_delta = centered_per - identity_per
    identity_pds = float(np.mean(identity_per))
    centered_pds = float(np.mean(centered_per))
    raw_delta = float(np.mean(per_target_delta))
    bootstrap = target_cluster_bootstrap(per_target_delta)
    gate = promotion_gate(raw_delta, bootstrap)
    write_per_target(
        output / "per_target_pds.csv",
        held.targets,
        identity_per,
        centered_per,
    )

    result: dict[str, object] = {
        "schema_version": 1,
        "experiment": EXPERIMENT,
        "status": "completed",
        "decision": "promote" if gate["passed"] else "reject",
        "frozen_method": {
            "preregistration": preregistration,
            "source_experiment": "075-postwhiten-local-centering",
            "source_result_sha256": EXPECTED_EXP075_RESULTS_SHA256,
            "source_selection_freeze_sha256": (EXPECTED_EXP075_SELECTION_FREEZE_SHA256),
            "source_contract": frozen075,
            "k": LOCAL_K,
            "beta": LOCAL_BETA,
            "minimum_mutual": MINIMUM_MUTUAL,
            "retuned_on_kolf": False,
            "candidate_count_evaluated_on_kolf": 1,
        },
        "protocol": {
            "metric": "raw effect-level rank PDS",
            "production_query_targets": len(panel),
            "production_query_axis": (
                "all official VCC genes before scorer target-coordinate exclusion"
            ),
            "local_center_query_rows": int(identity.shape[0]),
            "local_center_axis_genes": int(identity.shape[1]),
            "held_truth": "KOLF2.1J",
            "held_gallery_targets": len(held.targets),
            "held_score_genes": len(held.score_genes),
            "full_shared_gallery": True,
            "query_freeze_before_kolf_open": True,
            "count_emission": False,
            "cell_eval2": False,
            "submission": False,
        },
        "scores": {
            "identity_raw_pds": identity_pds,
            "local_center_raw_pds": centered_pds,
            "raw_pds_delta": raw_delta,
            "target_cluster_bootstrap": bootstrap,
        },
        "promotion_gate": gate,
        "query": query_audit,
        "held_truth": held.audit,
        "leakage_audit": {
            "kolf_opened_after_identity_and_candidate_frozen": True,
            "kolf_rows_in_source_model": 0,
            "kolf_rows_in_whitening_operator": 0,
            "kolf_rows_in_query_neighbor_graph": 0,
            "kolf_values_used_for_tuning": False,
            "kolf_labels_used_by_transform": False,
            "production_query_built_for_all_300_before_kolf_subset": True,
            "all_300_vcc_target_coordinates_excluded_from_whitener": (
                query_audit["production"]["panel_coordinates_excluded_from_operator"]
                == EXPECTED_PANEL_TARGETS
            ),
            "all_300_vcc_target_coordinates_excluded_from_local_centering": (
                query_audit["transform"]["panel_coordinates_excluded"]
                == EXPECTED_PANEL_TARGETS
                and query_audit["transform"]["panel_coordinate_overlap"] == 0
            ),
        },
        "array_hashes": {
            "identity_eval": sha256_array(identity_eval),
            "local_center_eval": sha256_array(centered_eval),
            "kolf_truth": sha256_array(held.values),
            "per_target_delta": sha256_array(per_target_delta),
        },
        "input_hashes": {
            "challenge": EXPECTED_CHALLENGE_HASHES,
            "implementation": EXPECTED_IMPLEMENTATION_HASHES,
            "production_donors": EXPECTED_DONOR_HASHES,
            "kolf_held_truth": EXPECTED_KOLF_HASHES,
            "production_whitener": EXPECTED_PRODUCTION_WHITENER_SHA256,
            "production_whitener_metadata": (EXPECTED_PRODUCTION_WHITENER_JSON_SHA256),
        },
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
        "peak_rss_gib": peak_rss_gib(),
    }
    (output / "results.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(result["scores"], indent=2), flush=True)
    print(json.dumps(result["promotion_gate"], indent=2), flush=True)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donor-root", type=Path, required=True)
    parser.add_argument("--kolf", type=Path, required=True)
    parser.add_argument("--whitener", type=Path, required=True)
    parser.add_argument("--whitener-metadata", type=Path, required=True)
    parser.add_argument("--exp075-results", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-rss-gib", type=float, default=32.0)
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing real-data execution without --guarded")
    if not (1.0 <= args.max_rss_gib <= 40.0):
        raise SystemExit("--max-rss-gib must be in [1, 40]")
    run_experiment(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

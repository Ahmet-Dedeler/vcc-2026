#!/usr/bin/env python3
"""Experiment 140: target-UMI-weighted batch-matched donor controls.

Stages:

* plan: freeze target folds, nulls, bootstrap draws, axes, and source schemas.
* freeze: hash all implementation and source inputs before outcome access.
* prepare: read controls plus target batch/depth metadata only; build corrections.
* identity: reproduce the eta-zero incumbent and freeze the outer whitened geometry.
* authorize-development: creates a lock only after explicit post-review confirmation.
* develop: opens the 613-target source-development outcomes and runs the frozen design.
* audit: validate all existing artifacts without opening an outcome.

There is intentionally no command that loads Experiment 132's 410-target test effects.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gc
import hashlib
import importlib.util
import json
import os
import re
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

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
README = HERE / "README.md"
PROTOCOL = HERE / "preregistered_protocol.json"
TESTS = HERE / "test_run.py"
PREFREEZE_AMENDMENT = HERE / "PREFREEZE_AMENDMENT_001.json"
PREFREEZE_AMENDMENT_002 = HERE / "PREFREEZE_AMENDMENT_002.json"
PREFREEZE_AMENDMENT_003 = HERE / "PREFREEZE_AMENDMENT_003.json"
PREFREEZE_AMENDMENT_004 = HERE / "PREFREEZE_AMENDMENT_004.json"
PREFREEZE_AMENDMENT_005 = HERE / "PREFREEZE_AMENDMENT_005.json"
PLAN_DIR = HERE / "work" / "preregistered-plan-v2"
SOURCE_HASH_MANIFEST = HERE / "SOURCE_HASH_MANIFEST.v5.json"
PREOUTCOME_AUTHORIZATION = HERE / "PREOUTCOME_AUTHORIZATION.v5.json"
PREPARE_DIR = HERE / "work" / "source-preparation-v5"
IDENTITY_DIR = HERE / "work" / "identity-reproduction-v5"
DEVELOPMENT_AUTHORIZATION = HERE / "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json"
DEVELOPMENT_LOCK = HERE / "work" / "development-outcome-access.v1.lock"
DEVELOPMENT_DIR = HERE / "work" / "development-result-v1"

EXP132 = ROOT / "experiments/132-response-decomposition-ridge"
PARENT_RUN = EXP132 / "run.py"
PARENT_AUTH = EXP132 / "SOURCE_ONLY_AUTHORIZATION.v2.json"
PARENT_PLAN = EXP132 / "work/preregistered-freeze-v2"
ALL_TARGETS = PARENT_PLAN / "all_targets.txt"
DEV_TARGETS = PARENT_PLAN / "selection_targets.txt"
TEST_TARGETS = PARENT_PLAN / "test_targets.txt"
RESPONSE_GENES = PARENT_PLAN / "response_genes.txt"
PANEL_PATH = ROOT / "data/challenge/pert_counts.csv"
CANONICAL_GENES_PATH = ROOT / "data/challenge/gene_names.csv"

K562_RAW = ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
ORION_DATA = ROOT / "data/external/orion/data"
ORION_GENE_METADATA = ROOT / "data/external/orion/metadata/gene_metadata.parquet"
ESSENTIAL_LIST = ROOT / "data/donors/k562_essential/perts.txt"
A15_META = ROOT / "submissions/a15-cd4f010.vcc.meta.json"
PRODUCTION_OPERATOR_META = ROOT / "data/donors/whitener_union_pruned.json"
PRODUCTION_OPERATOR = ROOT / "data/donors/whitener_union_pruned.npz"
PRODUCTION_OPERATOR_BUILDER = ROOT / "scripts/build_union_whitener.py"

DONOR_ROOTS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
PRODUCTION_DONORS = tuple(DONOR_ROOTS)
OUTER_DONORS = {
    "k562": ("orion_hct116", "orion_hek293t"),
    "rpe1": PRODUCTION_DONORS,
    "jurkat": PRODUCTION_DONORS,
}
SCORED_DATASETS = ("k562_gwps", "k562_essential", "rpe1", "nadig_jurkat")
DATASET_GROUP = {
    "k562_gwps": "k562",
    "k562_essential": "k562",
    "rpe1": "rpe1",
    "nadig_jurkat": "jurkat",
}

EXPECTED = {
    PARENT_RUN: "e78e079fc3b4a715e8b412b0a48bfe44a2f3f6cfcadb8a65a789351c53877778",
    PARENT_AUTH: "9bca6e620359e44b5440d0b8588412f40f274304ae67d9f55fd7e7d03341e197",
    PARENT_PLAN
    / "plan-manifest.json": "170d55ed5f812eb85f5503af4bb9d7dec8d17dce6e58fe966d85bf0eec5185d1",
    DEV_TARGETS: "74a070a538455ffc2a528dae3094d21e047905463779b2a5206d681f5ca6ebb5",
    TEST_TARGETS: "bd37bcacb4eea90da11902bbb8f4736068a7ff67c001739e3407de2b4bf73fc6",
    RESPONSE_GENES: "24743754eead5b8ddd81634bd12feac1df431eb2b9581197755b75c4c90949c8",
    PANEL_PATH: "f57edd7b912ebd718efc7ee9d0f334772513e7cc418d133ce525470e373b3276",
    CANONICAL_GENES_PATH: "25bfa66715e186bebabce7ac788bbcea47e2bf59ca70be1f8f3a06f2f0e47201",
    ESSENTIAL_LIST: "b1e5279534eb240f0f3487f0814d379e05b12799b21f6cb3818f5c7f55557e2d",
    A15_META: "34223442e4e8816bb9be2b1fdfe6a563ed114c17a2893bee0ff603953f84287a",
    PRODUCTION_OPERATOR_META: "838f263c6a4e1bed069d5faeb3bc3bb2b9234c5a5a933219a94241037a801f4d",
    PRODUCTION_OPERATOR: "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123",
    PRODUCTION_OPERATOR_BUILDER: "2337092687eaa973e2d141a869dfdb45969b4d4c48e9c0b166bb738c1eb1116e",
    DONOR_ROOTS["k562_gwps"]
    / "delta.npy": "0c9422280fed7afb81ca474e8b9024af788c9ffbde32e11f96cb52422e296aa9",
    DONOR_ROOTS["k562_gwps"]
    / "b_ctrl.npy": "d06cf51d51cde1ae999e86bbd3be0bdc15471392bae153a1aa1e21dc5c81b957",
    DONOR_ROOTS["orion_hct116"]
    / "delta.npy": "72c7a01dc2869bd3e94c37e89b158c3aaf4048a0dff63da6542afa2de4ec2b33",
    DONOR_ROOTS["orion_hct116"]
    / "b_ctrl.npy": "b3044d439cfac932024e0a48152196387714e4b460a339fb9a03e32fcb2fe069",
    DONOR_ROOTS["orion_hek293t"]
    / "delta.npy": "d95202d984d9322dc3830f23b6ef55001ba8562dea1d72c43fb0c9d5064c7d1b",
    DONOR_ROOTS["orion_hek293t"]
    / "b_ctrl.npy": "e04234cffbe102731cc562e8fb06bddb8ef0644b2af894b8a0ee6ce079ca4984",
}

TS = 50_000.0
DEV_COUNT = 613
TEST_COUNT = 410
ALL_TARGET_COUNT = 1023
RESPONSE_COUNT = 4992
FOLDS = 5
ETAS = (0.0, 0.25, 0.5, 1.0)
COMPONENT_MODES = ("all", "orion_only", "k562_only")
NULL_DRAWS = 49
BOOTSTRAP_DRAWS = 20_000
TARGET_NULL_SEED = 140_202_610
BATCH_NULL_SEED = 140_202_611
BOOTSTRAP_SEED = 140_202_612
WHITEN_GAMMA = 0.65
WHITEN_FLOOR = 1e-12
BASIS_PER_DONOR_CAP = 12_000
BASIS_SEED = 0
MIN_GAIN = 0.003
MAX_CONCENTRATION = 0.5
CONTROL_RECON_ATOL = 2e-6
EPS = 1e-12
GIB = 1024**3
DEVELOPMENT_CONFIRM = "I_CONFIRM_EXPERIMENT_140_DEVELOPMENT_OUTCOME"
ORION_CONTROL_LABEL = "Non-Targeting"
ORION_MISCLASSIFIED_TARGET = "NEGR1"
ORION_CONTROL_SCAN_LABELS = (
    ORION_CONTROL_LABEL,
    ORION_MISCLASSIFIED_TARGET,
)
ORION_LEGACY_CONTROL_PATTERNS = (
    r"^non.?targeting",
    r"^ntc",
    r"^control",
    r"^safe.?harbou?r",
    r"^neg",
)
ORION_EXPECTED_CONTROL_COUNTS = {
    "HCT116": {
        "exact_non_targeting": 165_777,
        "legacy_negr1": 154,
        "legacy_total": 165_931,
    },
    "HEK293T": {
        "exact_non_targeting": 218_838,
        "legacy_negr1": 233,
        "legacy_total": 219_071,
    },
}


class ContractError(RuntimeError):
    """Raised when a frozen experimental boundary is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def import_parent() -> Any:
    require(PARENT_RUN.is_file(), "Experiment 132 implementation is missing")
    digest = hashlib.sha256(PARENT_RUN.read_bytes()).hexdigest()
    require(
        digest == EXPECTED[PARENT_RUN],
        "Experiment 132 implementation changed before import",
    )
    spec = importlib.util.spec_from_file_location("experiment132_for_140", PARENT_RUN)
    require(spec is not None and spec.loader is not None, "cannot import Experiment 132")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


PARENT = import_parent()


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if sys.platform == "darwin" else value / 1024**2


def memory_check(label: str, limit: float) -> None:
    value = peak_rss_gib()
    print(f"[rss] {label}: {value:.3f} GiB", flush=True)
    if limit > 0 and value >= limit:
        raise MemoryError(f"{label}: peak RSS {value:.3f} GiB >= {limit:.3f} GiB")


def sha256_file(path: Path, block: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while piece := handle.read(block):
            digest.update(piece)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def read_lines(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(x.strip() for x in path.read_text().splitlines() if x.strip())
    require(bool(values), f"empty axis: {path}")
    if unique:
        require(len(values) == len(set(values)), f"duplicate axis: {path}")
    return values


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames is not None and column in reader.fieldnames, f"missing {column}")
        values = tuple(
            str(row[column]).strip()
            for row in reader
            if row.get(column) and str(row[column]).strip()
        )
    require(values and len(values) == len(set(values)), f"invalid {column}: {path}")
    return values


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


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
        np.savez_compressed(temporary, **arrays)
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def save_npy(path: Path, array: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as handle:
        np.save(handle, array)
    temporary.replace(path)


def file_record(path: Path, *, hash_file: bool = True) -> dict[str, Any]:
    result: dict[str, Any] = {"bytes": path.stat().st_size}
    if hash_file:
        result["sha256"] = sha256_file(path)
    return result


def stable_derangements(
    rng: np.random.Generator, draws: int, size: int
) -> np.ndarray:
    base = np.arange(size, dtype=np.int16)
    result = np.empty((draws, size), dtype=np.int16)
    for draw in range(draws):
        for _ in range(100_000):
            candidate = rng.permutation(base)
            if np.all(candidate != base):
                result[draw] = candidate
                break
        else:
            raise ContractError("failed to generate derangement")
    return result


def contiguous_folds(size: int) -> np.ndarray:
    labels = np.empty(size, dtype=np.int8)
    for fold, indices in enumerate(np.array_split(np.arange(size), FOLDS)):
        labels[indices] = fold
    return labels


def is_orion_control(value: str) -> bool:
    return str(value) == ORION_CONTROL_LABEL


def is_legacy_orion_control(value: str) -> bool:
    normalized = str(value).lower()
    return any(re.match(pattern, normalized) for pattern in ORION_LEGACY_CONTROL_PATTERNS)


def numbered_files(line: str) -> tuple[Path, ...]:
    files = sorted(
        ORION_DATA.glob(f"{line}_Batch*.parquet"),
        key=lambda path: int(re.search(r"Batch(\d+)", path.name)[1]),
    )
    require(bool(files), f"no Orion files for {line}")
    return tuple(files)


def basis_selection(donors: Sequence[str]) -> dict[str, np.ndarray]:
    """Reproduce a15 union_pruned row selection inside an outer source set."""
    essential = set(read_lines(ESSENTIAL_LIST))
    rng = np.random.default_rng(BASIS_SEED)
    result = {}
    for name in donors:
        perts = read_lines(DONOR_ROOTS[name] / "perts.txt")
        if len(perts) > BASIS_PER_DONOR_CAP:
            rows = np.sort(
                rng.choice(len(perts), BASIS_PER_DONOR_CAP, replace=False)
            ).astype(np.int64)
        else:
            rows = np.arange(len(perts), dtype=np.int64)
        pert_index = {pert: index for index, pert in enumerate(perts)}
        row_name: list[str | None] = [None] * len(perts)
        for pert, index in pert_index.items():
            row_name[index] = pert
        drop = np.asarray(
            [row_name[int(row)] in essential for row in rows], dtype=bool
        )
        result[name] = rows[~drop]
    return result


def basis_selection_report(donors: Sequence[str]) -> dict[str, Any]:
    selected = basis_selection(donors)
    source_rows = {
        name: len(read_lines(DONOR_ROOTS[name] / "perts.txt")) for name in donors
    }
    return {
        "donors": list(donors),
        "basis_per_donor_cap": BASIS_PER_DONOR_CAP,
        "seed": BASIS_SEED,
        "post_cap_prune": str(ESSENTIAL_LIST.relative_to(ROOT)),
        "source_rows": source_rows,
        "rows_kept": {name: int(len(rows)) for name, rows in selected.items()},
        "rows_dropped_after_cap": {
            name: int(min(source_rows[name], BASIS_PER_DONOR_CAP) - len(rows))
            for name, rows in selected.items()
        },
        "rows_kept_total": int(sum(len(rows) for rows in selected.values())),
        "row_index_sha256": {
            name: hashlib.sha256(rows.astype(np.int64).tobytes()).hexdigest()
            for name, rows in selected.items()
        },
    }


def donor_paths(name: str) -> tuple[Path, ...]:
    root = DONOR_ROOTS[name]
    return tuple(root / item for item in ("delta.npy", "b_ctrl.npy", "perts.txt", "genes.txt", "meta.json"))


def schema_report() -> dict[str, Any]:
    for path, expected in EXPECTED.items():
        require(path.is_file(), f"missing pinned input: {path}")
        require(sha256_file(path) == expected, f"pinned input changed: {path}")
    with h5py.File(K562_RAW, "r") as handle:
        perturbation = handle["obs"]["perturbation"]
        categories = np.asarray(perturbation["categories"]).astype(str)
        codes = np.asarray(perturbation["codes"], dtype=np.int64)
        batches = np.asarray(handle["obs"]["batch"], dtype=np.int64)
        batch_values = np.unique(batches)
        control_code = int(np.flatnonzero(categories == "control")[0])
        control = codes == control_code
        control_per_batch = np.asarray(
            [np.sum(control & (batches == batch)) for batch in batch_values],
            dtype=np.int64,
        )
        target_set = set(read_lines(DONOR_ROOTS["k562_gwps"] / "perts.txt"))
        target_spans = np.asarray(
            [
                len(np.unique(batches[codes == code]))
                for code, target in enumerate(categories)
                if target in target_set
            ],
            dtype=np.int64,
        )
        control_rows = np.flatnonzero(control)[:128]
        control_x_sum = np.asarray(
            handle["X"][control_rows, :], dtype=np.float64
        ).sum(axis=1)
        control_ncounts = np.asarray(
            handle["obs"]["ncounts"][control_rows], dtype=np.float64
        )
        control_umi_count = np.asarray(
            handle["obs"]["UMI_count"][control_rows], dtype=np.float64
        )
        report_k562 = {
            "bytes": K562_RAW.stat().st_size,
            "x_shape": list(handle["X"].shape),
            "x_dtype": str(handle["X"].dtype),
            "x_chunks": list(handle["X"].chunks or ()),
            "obs_keys": sorted(handle["obs"].keys()),
            "perturbation_categories": int(len(perturbation["categories"])),
            "batches": int(len(batch_values)),
            "control_cells": int(control.sum()),
            "control_batches": int(np.unique(batches[control]).size),
            "control_cells_per_batch_median": float(np.median(control_per_batch)),
            "target_batches_median": float(np.median(target_spans)),
            "depth_field": "ncounts",
            "excluded_depth_field": "UMI_count",
            "control_sample_ncounts_max_abs_vs_x_sum": float(
                np.max(np.abs(control_ncounts - control_x_sum))
            ),
            "control_sample_umi_count_max_abs_vs_x_sum": float(
                np.max(np.abs(control_umi_count - control_x_sum))
            ),
        }
    orion: dict[str, Any] = {}
    for line in ("HCT116", "HEK293T"):
        files = numbered_files(line)
        first = pq.ParquetFile(files[0])
        exact_count = 0
        legacy_count = 0
        negr1_count = 0
        sample_parts = []
        sample_remaining = 256
        for path in files:
            metadata = pq.read_table(
                path,
                columns=[
                    "gene_target",
                    "pass_guide_filter",
                ],
            ).to_pandas()
            metadata = metadata[metadata["pass_guide_filter"] == 1]
            labels = metadata["gene_target"].astype(str)
            exact = labels.map(is_orion_control).to_numpy()
            legacy = labels.map(is_legacy_orion_control).to_numpy()
            exact_count += int(exact.sum())
            legacy_count += int(legacy.sum())
            negr1_count += int((labels == ORION_MISCLASSIFIED_TARGET).sum())
            if sample_remaining > 0 and np.any(exact):
                frame = pq.read_table(
                    path,
                    columns=[
                        "gene_expression",
                        "gene_target",
                        "pass_guide_filter",
                        "total_counts",
                    ],
                ).to_pandas()
                frame = frame[frame["pass_guide_filter"] == 1]
                sample_exact = (
                    frame["gene_target"].astype(str).map(is_orion_control).to_numpy()
                )
                part = frame.loc[sample_exact].iloc[:sample_remaining]
                sample_parts.append(part)
                sample_remaining -= len(part)
        sample = pd.concat(sample_parts, ignore_index=True)
        expression_sum = np.fromiter(
            (float(np.sum(values)) for values in sample["gene_expression"]),
            dtype=np.float64,
            count=len(sample),
        )
        observed_counts = {
            "exact_non_targeting": exact_count,
            "legacy_negr1": negr1_count,
            "legacy_total": legacy_count,
        }
        require(
            observed_counts == ORION_EXPECTED_CONTROL_COUNTS[line],
            f"{line}: Orion control-count audit changed",
        )
        orion[line] = {
            "files": len(files),
            "first": files[0].name,
            "last": files[-1].name,
            "columns": first.schema_arrow.names,
            "depth_field": "total_counts",
            "batch_definition": "one Batch*.parquet file",
            "candidate_control_label": ORION_CONTROL_LABEL,
            "legacy_control_patterns": list(ORION_LEGACY_CONTROL_PATTERNS),
            "control_counts": observed_counts,
            "negr1_candidate_control": False,
            "negr1_targeting_row": True,
            "control_sample_total_counts_max_abs_vs_expression_sum": float(
                np.max(
                    np.abs(
                        sample["total_counts"].to_numpy(dtype=np.float64)
                        - expression_sum
                    )
                )
            ),
        }
    artifacts = {}
    for name, root in DONOR_ROOTS.items():
        perts = read_lines(root / "perts.txt")
        genes = read_lines(root / "genes.txt", unique=False)
        delta = np.load(root / "delta.npy", mmap_mode="r")
        ctrl = np.load(root / "b_ctrl.npy", mmap_mode="r")
        require(delta.shape == (len(perts), len(genes)), f"{name}: artifact shape")
        require(ctrl.shape == (len(genes),), f"{name}: control shape")
        artifacts[name] = {
            "targets": len(perts),
            "genes": len(genes),
            "delta_dtype": str(delta.dtype),
            "control_dtype": str(ctrl.dtype),
        }
    a15 = read_json(A15_META)
    operator = read_json(PRODUCTION_OPERATOR_META)
    require(a15["config"]["whitener_id"] == "union_pruned", "a15 whitener id changed")
    require(
        a15["operator"]["sha256"] == EXPECTED[PRODUCTION_OPERATOR],
        "a15 operator hash changed",
    )
    require(
        a15["operator"]["metadata_sha256"] == EXPECTED[PRODUCTION_OPERATOR_META],
        "a15 operator metadata hash changed",
    )
    require(operator["basis_per_donor_cap"] == BASIS_PER_DONOR_CAP, "basis cap changed")
    require(operator["seed"] == BASIS_SEED, "basis seed changed")
    require(operator["basis_prune"] == "essential", "basis prune changed")
    require(
        operator["rows_kept_per_donor"]
        == basis_selection_report(PRODUCTION_DONORS)["rows_kept"],
        "production basis selection changed",
    )
    return {
        "k562": report_k562,
        "orion": orion,
        "artifacts": artifacts,
        "incumbent_provenance": {
            "a15_meta": file_record(A15_META),
            "operator": file_record(PRODUCTION_OPERATOR),
            "operator_meta": file_record(PRODUCTION_OPERATOR_META),
            "operator_builder": file_record(PRODUCTION_OPERATOR_BUILDER),
            "essential_list": file_record(ESSENTIAL_LIST),
            "a15_whitener_id": a15["config"]["whitener_id"],
            "production_basis": basis_selection_report(PRODUCTION_DONORS),
            "outer_basis": {
                group: basis_selection_report(donors)
                for group, donors in {
                    "k562": OUTER_DONORS["k562"],
                    "non_k562": OUTER_DONORS["rpe1"],
                }.items()
            },
        },
    }


def create_plan() -> None:
    require(not PLAN_DIR.exists(), "plan already exists")
    require(not PREOUTCOME_AUTHORIZATION.exists(), "authorization already exists")
    require(not PREPARE_DIR.exists() and not IDENTITY_DIR.exists(), "later artifact exists")
    development = read_lines(DEV_TARGETS)
    untouched = read_lines(TEST_TARGETS)
    all_targets = read_lines(ALL_TARGETS)
    responses = read_lines(RESPONSE_GENES)
    panel = set(read_csv_column(PANEL_PATH, "target_gene"))
    require(len(development) == DEV_COUNT, "development target count changed")
    require(len(untouched) == TEST_COUNT, "test target count changed")
    require(len(all_targets) == ALL_TARGET_COUNT, "gallery target count changed")
    require(len(responses) == RESPONSE_COUNT, "response count changed")
    require(set(development).isdisjoint(untouched), "development/test overlap")
    require(set(development) | set(untouched) == set(all_targets), "gallery union changed")
    require(set(responses).isdisjoint(panel), "response axis includes VCC coordinate")
    require(set(responses).isdisjoint(all_targets), "response axis includes gallery coordinate")
    for name, root in DONOR_ROOTS.items():
        donor_targets = set(read_lines(root / "perts.txt"))
        require(set(development) <= donor_targets, f"{name}: development target missing")

    folds = contiguous_folds(DEV_COUNT)
    features = PARENT.feature_rows(development)
    blocks = PARENT.reciprocal_feature_blocks(features)
    block_count = int(blocks.max()) + 1
    target_rng = np.random.default_rng(TARGET_NULL_SEED)
    batch_rng = np.random.default_rng(BATCH_NULL_SEED)
    bootstrap_rng = np.random.default_rng(BOOTSTRAP_SEED)
    target_nulls = stable_derangements(target_rng, NULL_DRAWS, DEV_COUNT)
    batch_nulls = np.stack(
        [
            stable_derangements(batch_rng, NULL_DRAWS, DEV_COUNT)
            for _ in PRODUCTION_DONORS
        ],
        axis=1,
    )
    bootstrap = bootstrap_rng.integers(
        0,
        block_count,
        size=(BOOTSTRAP_DRAWS, block_count),
        dtype=np.uint16,
    )

    PLAN_DIR.mkdir(parents=True)
    atomic_npz(
        PLAN_DIR / "plan-arrays.npz",
        development_folds=folds,
        development_blocks=blocks,
        target_label_nulls=target_nulls,
        batch_mixture_nulls=batch_nulls,
        bootstrap_block_draws=bootstrap,
    )
    schema = schema_report()
    atomic_json(PLAN_DIR / "source-schema.json", schema)
    manifest = {
        "schema_version": 1,
        "created_before_development_outcome_access": True,
        "development_targets": DEV_COUNT,
        "untouched_test_targets": TEST_COUNT,
        "response_genes": RESPONSE_COUNT,
        "fold_sizes": [int(np.sum(folds == fold)) for fold in range(FOLDS)],
        "target_blocks": block_count,
        "targets_in_multi_target_blocks": int(
            np.bincount(blocks)[np.bincount(blocks) > 1].sum()
        ),
        "maximum_target_block_size": int(np.bincount(blocks).max()),
        "etas": list(ETAS),
        "null_draws_each": NULL_DRAWS,
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "seeds": {
            "target_label": TARGET_NULL_SEED,
            "batch_mixture": BATCH_NULL_SEED,
            "bootstrap": BOOTSTRAP_SEED,
        },
        "outer_operator_basis": {
            group: basis_selection_report(donors)
            for group, donors in {
                "k562": OUTER_DONORS["k562"],
                "non_k562": OUTER_DONORS["rpe1"],
            }.items()
        },
        "files": {
            path.name: file_record(path)
            for path in sorted(PLAN_DIR.iterdir())
            if path.name != "plan-manifest.json"
        },
    }
    atomic_json(PLAN_DIR / "plan-manifest.json", manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))


def validate_plan() -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    require(PLAN_DIR.is_dir(), "plan missing")
    manifest = read_json(PLAN_DIR / "plan-manifest.json")
    require(manifest["development_targets"] == DEV_COUNT, "development count drift")
    require(manifest["untouched_test_targets"] == TEST_COUNT, "test count drift")
    require(manifest["response_genes"] == RESPONSE_COUNT, "response count drift")
    require(manifest["etas"] == list(ETAS), "eta grid drift")
    for name, expected in manifest["files"].items():
        path = PLAN_DIR / name
        require(path.stat().st_size == expected["bytes"], f"plan bytes changed: {name}")
        require(sha256_file(path) == expected["sha256"], f"plan hash changed: {name}")
    with np.load(PLAN_DIR / "plan-arrays.npz", allow_pickle=False) as archive:
        arrays = {name: np.asarray(archive[name]) for name in archive.files}
    require(arrays["development_folds"].shape == (DEV_COUNT,), "fold shape")
    require(arrays["target_label_nulls"].shape == (NULL_DRAWS, DEV_COUNT), "target null shape")
    require(
        arrays["batch_mixture_nulls"].shape
        == (NULL_DRAWS, len(PRODUCTION_DONORS), DEV_COUNT),
        "batch null shape",
    )
    require(
        manifest["outer_operator_basis"]["non_k562"]["rows_kept"]
        == {
            "k562_gwps": 7784,
            "orion_hct116": 10807,
            "orion_hek293t": 10803,
        },
        "three-donor union_pruned basis changed",
    )
    require(
        manifest["outer_operator_basis"]["k562"]["rows_kept_total"] == 21610,
        "K562-held union_pruned basis changed",
    )
    return manifest, arrays


def source_input_paths() -> tuple[Path, ...]:
    paths = [
        K562_RAW,
        ORION_GENE_METADATA,
        PARENT_RUN,
        PARENT_AUTH,
        PARENT_PLAN / "plan-manifest.json",
        ALL_TARGETS,
        DEV_TARGETS,
        TEST_TARGETS,
        RESPONSE_GENES,
        PANEL_PATH,
        CANONICAL_GENES_PATH,
        ESSENTIAL_LIST,
        A15_META,
        PRODUCTION_OPERATOR_META,
        PRODUCTION_OPERATOR,
        PRODUCTION_OPERATOR_BUILDER,
    ]
    paths.extend(numbered_files("HCT116"))
    paths.extend(numbered_files("HEK293T"))
    for name in PRODUCTION_DONORS:
        paths.extend(donor_paths(name))
    for name in ("k562_essential", "rpe1", "nadig_jurkat"):
        root = PARENT.SOURCE_ROOTS[name]
        paths.extend((root / "delta.npy", root / "perts.txt", root / "genes.txt"))
    paths.extend((PARENT.DEPMAP_PATH, PARENT.DEPMAP_RECEIPT))
    return tuple(dict.fromkeys(paths))


def preparation_source_paths() -> tuple[Path, ...]:
    paths = [
        K562_RAW,
        ORION_GENE_METADATA,
        DEV_TARGETS,
        RESPONSE_GENES,
    ]
    paths.extend(numbered_files("HCT116"))
    paths.extend(numbered_files("HEK293T"))
    for name in PRODUCTION_DONORS:
        root = DONOR_ROOTS[name]
        paths.extend((root / "b_ctrl.npy", root / "perts.txt", root / "genes.txt"))
    return tuple(dict.fromkeys(paths))


def identity_source_paths() -> tuple[Path, ...]:
    paths = [
        DEV_TARGETS,
        RESPONSE_GENES,
        CANONICAL_GENES_PATH,
        ESSENTIAL_LIST,
        A15_META,
        PRODUCTION_OPERATOR_META,
        PRODUCTION_OPERATOR,
        PRODUCTION_OPERATOR_BUILDER,
    ]
    for name in PRODUCTION_DONORS:
        paths.extend(donor_paths(name))
    return tuple(dict.fromkeys(paths))


def development_source_paths() -> tuple[Path, ...]:
    paths = [
        PARENT_RUN,
        PARENT_AUTH,
        PARENT_PLAN / "plan-manifest.json",
        DEV_TARGETS,
        RESPONSE_GENES,
    ]
    for name in ("k562_gwps", "k562_essential", "rpe1", "nadig_jurkat"):
        root = PARENT.SOURCE_ROOTS[name]
        paths.extend((root / "delta.npy", root / "perts.txt", root / "genes.txt"))
    return tuple(dict.fromkeys(paths))


def validate_source_hashes(paths: Iterable[Path]) -> None:
    manifest = read_json(SOURCE_HASH_MANIFEST)["files"]
    for path in tuple(dict.fromkeys(paths)):
        relative = str(path.relative_to(ROOT))
        require(relative in manifest, f"source missing from manifest: {relative}")
        expected = manifest[relative]
        require(path.is_file(), f"source missing: {relative}")
        require(path.stat().st_size == expected["bytes"], f"source size changed: {relative}")
        require(sha256_file(path) == expected["sha256"], f"source hash changed: {relative}")


def freeze() -> None:
    require(not SOURCE_HASH_MANIFEST.exists(), "source hash manifest exists")
    require(not PREOUTCOME_AUTHORIZATION.exists(), "authorization exists")
    require(not PREPARE_DIR.exists() and not IDENTITY_DIR.exists(), "later artifact exists")
    schema_report()
    plan, _ = validate_plan()
    sources = {}
    for path in source_input_paths():
        require(path.is_file(), f"missing source input: {path}")
        relative = str(path.relative_to(ROOT))
        print(f"[hash] {relative}", flush=True)
        sources[relative] = file_record(path)
    atomic_json(
        SOURCE_HASH_MANIFEST,
        {
            "schema_version": 1,
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "created_before_development_outcome_access": True,
            "files": sources,
        },
    )
    implementation = {
        str(path.relative_to(ROOT)): sha256_file(path)
        for path in (
            README,
            PROTOCOL,
            HERE / "run.py",
            TESTS,
            PREFREEZE_AMENDMENT,
            PREFREEZE_AMENDMENT_002,
            PREFREEZE_AMENDMENT_003,
            PREFREEZE_AMENDMENT_004,
            PREFREEZE_AMENDMENT_005,
        )
    }
    authorization = {
        "schema_version": 1,
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "preoutcome_work_authorized": True,
        "development_outcome_authorized": False,
        "untouched_test_outcome_authorized": False,
        "cell_eval2_counts_build_submission_authorized": False,
        "protocol_sha256": sha256_file(PROTOCOL),
        "plan_manifest_sha256": sha256_file(PLAN_DIR / "plan-manifest.json"),
        "source_hash_manifest_sha256": sha256_file(SOURCE_HASH_MANIFEST),
        "implementation": implementation,
        "plan": plan,
        "peak_rss_gib": peak_rss_gib(),
    }
    atomic_json(PREOUTCOME_AUTHORIZATION, authorization)
    print(json.dumps(authorization, indent=2, sort_keys=True))


def validate_preoutcome(*, rehash_sources: bool = False) -> dict[str, Any]:
    require(PREOUTCOME_AUTHORIZATION.is_file(), "preoutcome authorization missing")
    auth = read_json(PREOUTCOME_AUTHORIZATION)
    require(auth["preoutcome_work_authorized"] is True, "preoutcome work disabled")
    require(auth["development_outcome_authorized"] is False, "development unexpectedly authorized")
    require(auth["untouched_test_outcome_authorized"] is False, "test unexpectedly authorized")
    require(auth["protocol_sha256"] == sha256_file(PROTOCOL), "protocol changed after freeze")
    require(
        auth["plan_manifest_sha256"] == sha256_file(PLAN_DIR / "plan-manifest.json"),
        "plan changed after freeze",
    )
    require(
        auth["source_hash_manifest_sha256"] == sha256_file(SOURCE_HASH_MANIFEST),
        "source manifest changed after freeze",
    )
    for relative, expected in auth["implementation"].items():
        require(sha256_file(ROOT / relative) == expected, f"implementation changed: {relative}")
    manifest = read_json(SOURCE_HASH_MANIFEST)
    for relative, expected in manifest["files"].items():
        path = ROOT / relative
        require(path.is_file(), f"source missing: {relative}")
        require(path.stat().st_size == expected["bytes"], f"source size changed: {relative}")
        if rehash_sources:
            require(sha256_file(path) == expected["sha256"], f"source hash changed: {relative}")
    validate_plan()
    return auth


def last_columns(source_genes: Sequence[str], wanted: Sequence[str]) -> np.ndarray:
    index = {gene: position for position, gene in enumerate(source_genes)}
    require(all(gene in index for gene in wanted), "response gene missing from donor")
    return np.asarray([index[gene] for gene in wanted], dtype=np.int64)


def h5_categorical(group: h5py.Group) -> tuple[np.ndarray, np.ndarray]:
    categories = np.asarray(group["categories"]).astype(str)
    codes = np.asarray(group["codes"], dtype=np.int64)
    require(np.all(codes >= 0), "categorical contains missing codes")
    return categories, codes


@dataclass
class BatchSummary:
    batch_labels: tuple[str, ...]
    matched_control_sums: np.ndarray
    matched_control_cells: np.ndarray
    reconstruction_control_sums: np.ndarray
    reconstruction_control_cells: np.ndarray
    target_umi: np.ndarray
    target_cells: np.ndarray
    control_policy: dict[str, Any]
    depth_validation: dict[str, Any]


def scan_k562(limit: float) -> BatchSummary:
    root = DONOR_ROOTS["k562_gwps"]
    perts = read_lines(root / "perts.txt")
    target_index = {target: row for row, target in enumerate(perts)}
    with h5py.File(K562_RAW, "r") as handle:
        categories, codes = h5_categorical(handle["obs"]["perturbation"])
        batches = np.asarray(handle["obs"]["batch"], dtype=np.int64)
        ncounts = np.asarray(handle["obs"]["ncounts"], dtype=np.float64)
        batch_values, batch_codes = np.unique(batches, return_inverse=True)
        n_batches = len(batch_values)
        category_to_target = np.asarray(
            [target_index.get(value, -1) for value in categories], dtype=np.int32
        )
        row_targets = category_to_target[codes]
        control_categories = np.asarray(
            [value == "control" for value in categories], dtype=bool
        )
        control_rows = control_categories[codes]
        target_umi = np.zeros((len(perts), n_batches), dtype=np.float64)
        target_cells = np.zeros((len(perts), n_batches), dtype=np.int32)
        active = row_targets >= 0
        flat = row_targets[active].astype(np.int64) * n_batches + batch_codes[active]
        np.add.at(target_umi.ravel(), flat, ncounts[active])
        np.add.at(target_cells.ravel(), flat, 1)

        x = handle["X"]
        control_sums = np.zeros((n_batches, x.shape[1]), dtype=np.float64)
        control_cells = np.zeros(n_batches, dtype=np.int64)
        max_depth_diff = 0.0
        block = 3886
        for lo in range(0, x.shape[0], block):
            hi = min(lo + block, x.shape[0])
            values = np.asarray(x[lo:hi], dtype=np.float64)
            row_sum = values.sum(axis=1)
            max_depth_diff = max(
                max_depth_diff,
                float(np.max(np.abs(row_sum - ncounts[lo:hi]))),
            )
            use = control_rows[lo:hi]
            if np.any(use):
                local_batch = batch_codes[lo:hi][use]
                indicator = sp.csr_matrix(
                    (
                        np.ones(int(use.sum()), dtype=np.float64),
                        (local_batch, np.arange(int(use.sum()))),
                    ),
                    shape=(n_batches, int(use.sum())),
                )
                control_sums += np.asarray(indicator @ values[use])
                control_cells += np.bincount(local_batch, minlength=n_batches)
            if (lo // block) % 40 == 39:
                print(f"[k562] {hi:,}/{x.shape[0]:,}", flush=True)
                memory_check("K562 control scan", limit)
    require(max_depth_diff == 0.0, "K562 ncounts does not exactly equal X row sum")
    require(np.all(control_cells > 0), "K562 batch without controls")
    return BatchSummary(
        batch_labels=tuple(map(str, batch_values)),
        matched_control_sums=control_sums,
        matched_control_cells=control_cells,
        reconstruction_control_sums=control_sums,
        reconstruction_control_cells=control_cells,
        target_umi=target_umi,
        target_cells=target_cells,
        control_policy={
            "candidate": "perturbation == control",
            "reconstruction": "same_as_candidate",
            "pools_identical": True,
        },
        depth_validation={
            "field": "ncounts",
            "max_abs_difference_vs_x_row_sum": max_depth_diff,
            "exact": True,
        },
    )


def orion_gene_axis() -> tuple[str, ...]:
    frame = pd.read_parquet(ORION_GENE_METADATA).sort_values("gene_token_id")
    size = int(frame["gene_token_id"].max()) + 1
    genes = np.asarray([""] * size, dtype=object)
    genes[frame["gene_token_id"].to_numpy(dtype=np.int64)] = frame["gene_name"].astype(str)
    return tuple(map(str, genes))


def release_arrow_memory() -> None:
    gc.collect()
    pa.default_memory_pool().release_unused()


def read_orion_control_sum(
    path: Path,
    label: str,
    n_genes: int,
) -> tuple[np.ndarray, int, float]:
    token_table = pq.read_table(
        path,
        columns=["gene_token_id"],
        filters=[
            ("pass_guide_filter", "=", 1),
            ("gene_target", "=", label),
        ],
    )
    rows = int(token_table.num_rows)
    if rows == 0:
        del token_table
        release_arrow_memory()
        return np.zeros(n_genes, dtype=np.float64), 0, 0.0
    token_lists = token_table["gene_token_id"].combine_chunks()
    tokens = np.array(token_lists.values, dtype=np.int32, copy=True)
    del token_table, token_lists
    release_arrow_memory()

    value_table = pq.read_table(
        path,
        columns=["gene_expression", "total_counts"],
        filters=[
            ("pass_guide_filter", "=", 1),
            ("gene_target", "=", label),
        ],
    )
    require(int(value_table.num_rows) == rows, f"{path.name}: filtered row-count drift")
    value_lists = value_table["gene_expression"].combine_chunks()
    values = np.asarray(value_lists.values, dtype=np.float64)
    require(len(tokens) == len(values), f"{path.name}: sparse token/value length drift")
    result = np.bincount(tokens, weights=values, minlength=n_genes)
    expression_total = float(values.sum())
    metadata_total = float(
        np.asarray(
            value_table["total_counts"].combine_chunks(),
            dtype=np.float64,
        ).sum()
    )
    difference = abs(expression_total - metadata_total)
    del value_table, value_lists, tokens, values
    release_arrow_memory()
    return result, rows, difference


def scan_orion(name: str, line: str, limit: float) -> BatchSummary:
    root = DONOR_ROOTS[name]
    perts = read_lines(root / "perts.txt")
    target_index = {target: row for row, target in enumerate(perts)}
    genes = orion_gene_axis()
    require(genes == read_lines(root / "genes.txt", unique=False), f"{name}: gene axis changed")
    files = numbered_files(line)
    target_umi = np.zeros((len(perts), len(files)), dtype=np.float64)
    target_cells = np.zeros((len(perts), len(files)), dtype=np.int32)
    matched_control_sums = np.zeros((len(files), len(genes)), dtype=np.float64)
    matched_control_cells = np.zeros(len(files), dtype=np.int64)
    reconstruction_control_sums = np.zeros((len(files), len(genes)), dtype=np.float64)
    reconstruction_control_cells = np.zeros(len(files), dtype=np.int64)
    max_matched_control_depth_diff = 0.0
    max_reconstruction_control_depth_diff = 0.0
    negr1_cells = 0
    for batch, path in enumerate(files):
        metadata_table = pq.read_table(
            path,
            columns=[
                "gene_target",
                "pass_guide_filter",
                "total_counts",
            ],
        )
        eligible = (
            np.asarray(
                metadata_table["pass_guide_filter"].combine_chunks(),
                dtype=np.int64,
            )
            == 1
        )
        labels = np.asarray(
            metadata_table["gene_target"].combine_chunks(),
            dtype=str,
        )[eligible]
        totals = np.asarray(
            metadata_table["total_counts"].combine_chunks(),
            dtype=np.float64,
        )[eligible]
        matched_control = np.fromiter(
            (is_orion_control(value) for value in labels),
            dtype=bool,
            count=len(labels),
        )
        reconstruction_control = np.fromiter(
            (is_legacy_orion_control(value) for value in labels),
            dtype=bool,
            count=len(labels),
        )
        require(
            np.all(~matched_control | reconstruction_control),
            f"{name}: exact control missing from legacy pool",
        )
        batch_negr1_cells = int(np.sum(labels == ORION_MISCLASSIFIED_TARGET))
        negr1_cells += batch_negr1_cells
        exact_sum, exact_cells, exact_depth_diff = read_orion_control_sum(
            path,
            ORION_CONTROL_LABEL,
            len(genes),
        )
        negr1_sum, expression_negr1_cells, negr1_depth_diff = read_orion_control_sum(
            path,
            ORION_MISCLASSIFIED_TARGET,
            len(genes),
        )
        require(
            exact_cells == int(matched_control.sum()),
            f"{name}: exact-control metadata/expression count mismatch",
        )
        require(
            expression_negr1_cells == batch_negr1_cells,
            f"{name}: NEGR1 metadata/expression count mismatch",
        )
        require(
            exact_cells + expression_negr1_cells == int(reconstruction_control.sum()),
            f"{name}: legacy-control metadata/expression count mismatch",
        )
        matched_control_sums[batch] = exact_sum
        matched_control_cells[batch] = exact_cells
        reconstruction_control_sums[batch] = exact_sum + negr1_sum
        reconstruction_control_cells[batch] = exact_cells + expression_negr1_cells
        max_matched_control_depth_diff = max(
            max_matched_control_depth_diff,
            exact_depth_diff,
        )
        max_reconstruction_control_depth_diff = max(
            max_reconstruction_control_depth_diff,
            exact_depth_diff + negr1_depth_diff,
        )
        targeting = ~matched_control
        if np.any(targeting):
            mapped = np.asarray(
                [target_index.get(value, -1) for value in labels[targeting]],
                dtype=np.int32,
            )
            keep = mapped >= 0
            target_umi[:, batch] = np.bincount(
                mapped[keep],
                weights=totals[targeting][keep],
                minlength=len(perts),
            )
            target_cells[:, batch] = np.bincount(
                mapped[keep], minlength=len(perts)
            )
        if (batch + 1) % 10 == 0 or batch + 1 == len(files):
            print(f"[{name}] {batch + 1}/{len(files)}", flush=True)
            memory_check(f"{name} control scan", limit)
        del (
            metadata_table,
            eligible,
            labels,
            totals,
            matched_control,
            reconstruction_control,
            exact_sum,
            negr1_sum,
        )
        release_arrow_memory()
    require(
        np.all(matched_control_cells > 0),
        f"{name}: batch without exact Non-Targeting controls",
    )
    require(
        np.all(reconstruction_control_cells > 0),
        f"{name}: batch without legacy reconstruction controls",
    )
    require(
        max_matched_control_depth_diff <= 1e-6,
        f"{name}: exact-control total_counts differs from expression sum",
    )
    require(
        max_reconstruction_control_depth_diff <= 1e-6,
        f"{name}: legacy-control total_counts differs from expression sum",
    )
    expected = ORION_EXPECTED_CONTROL_COUNTS[line]
    require(
        int(matched_control_cells.sum()) == expected["exact_non_targeting"],
        f"{name}: exact Non-Targeting count changed",
    )
    require(negr1_cells == expected["legacy_negr1"], f"{name}: NEGR1 count changed")
    require(
        int(reconstruction_control_cells.sum()) == expected["legacy_total"],
        f"{name}: legacy control count changed",
    )
    return BatchSummary(
        batch_labels=tuple(path.stem for path in files),
        matched_control_sums=matched_control_sums,
        matched_control_cells=matched_control_cells,
        reconstruction_control_sums=reconstruction_control_sums,
        reconstruction_control_cells=reconstruction_control_cells,
        target_umi=target_umi,
        target_cells=target_cells,
        control_policy={
            "candidate": f'gene_target == "{ORION_CONTROL_LABEL}"',
            "reconstruction": "legacy build_orion_donor regex",
            "legacy_patterns": list(ORION_LEGACY_CONTROL_PATTERNS),
            "exact_non_targeting_cells": int(matched_control_cells.sum()),
            "legacy_negr1_cells": negr1_cells,
            "legacy_total_cells": int(reconstruction_control_cells.sum()),
            "negr1_candidate_control": False,
            "negr1_targeting_row": True,
            "pools_identical": False,
        },
        depth_validation={
            "field": "total_counts",
            "matched_max_abs_batch_difference_vs_expression_sum": (
                max_matched_control_depth_diff
            ),
            "reconstruction_max_abs_batch_difference_vs_expression_sum": (
                max_reconstruction_control_depth_diff
            ),
            "exact_within_1e-6": True,
        },
    )


def build_correction(
    name: str,
    summary: BatchSummary,
    development: Sequence[str],
    responses: Sequence[str],
    output: Path,
    limit: float,
) -> dict[str, Any]:
    root = DONOR_ROOTS[name]
    perts = read_lines(root / "perts.txt")
    genes = read_lines(root / "genes.txt", unique=False)
    response_columns = last_columns(genes, responses)
    b_ctrl = np.asarray(np.load(root / "b_ctrl.npy", mmap_mode="r"), dtype=np.float64)
    matched_control_totals = summary.matched_control_sums.sum(axis=1)
    reconstruction_pooled = summary.reconstruction_control_sums.sum(axis=0)
    reconstructed = np.log1p(
        TS * reconstruction_pooled / reconstruction_pooled.sum()
    ).astype(np.float32)
    existing = np.asarray(b_ctrl, dtype=np.float32)
    control_error = np.abs(reconstructed - existing)
    require(
        float(control_error.max()) <= CONTROL_RECON_ATOL,
        f"{name}: global control reconstruction exceeds tolerance",
    )
    matched_pooled = summary.matched_control_sums.sum(axis=0)
    matched_global = np.log1p(TS * matched_pooled / matched_pooled.sum()).astype(
        np.float32
    )
    matched_vs_existing = np.abs(matched_global - existing)
    if summary.control_policy["pools_identical"]:
        require(
            float(matched_vs_existing.max()) == 0.0,
            f"{name}: identical control pools did not reproduce the incumbent exactly",
        )
    valid_control = matched_control_totals > 0
    target_totals = summary.target_umi.sum(axis=1)
    support = target_totals > 0
    support &= ~np.any((summary.target_umi > 0) & ~valid_control[None, :], axis=1)
    q = np.divide(
        summary.matched_control_sums[:, response_columns],
        matched_control_totals[:, None],
        out=np.zeros(
            (len(matched_control_totals), len(responses)), dtype=np.float64
        ),
        where=matched_control_totals[:, None] > 0,
    )
    correction_sum = np.zeros(len(responses), dtype=np.float64)
    development_index = {target: index for index, target in enumerate(development)}
    development_correction = np.zeros((len(development), len(responses)), dtype=np.float32)
    development_support = np.zeros(len(development), dtype=bool)
    development_umi = np.zeros(
        (len(development), len(matched_control_totals)), dtype=np.float64
    )
    development_cells = np.zeros(
        (len(development), len(matched_control_totals)), dtype=np.int32
    )
    row_map = {target: row for row, target in enumerate(perts)}
    for target, index in development_index.items():
        row = row_map[target]
        development_umi[index] = summary.target_umi[row]
        development_cells[index] = summary.target_cells[row]

    chunk = 256
    for lo in range(0, len(perts), chunk):
        hi = min(lo + chunk, len(perts))
        weights = np.divide(
            summary.target_umi[lo:hi],
            target_totals[lo:hi, None],
            out=np.zeros_like(summary.target_umi[lo:hi]),
            where=target_totals[lo:hi, None] > 0,
        )
        matched = weights @ q
        correction = (
            b_ctrl[response_columns][None, :] - np.log1p(TS * matched)
        )
        local_support = support[lo:hi]
        correction[~local_support] = 0.0
        correction_sum += correction.sum(axis=0)
        for local, target in enumerate(perts[lo:hi]):
            dev = development_index.get(target)
            if dev is not None:
                development_correction[dev] = correction[local].astype(np.float32)
                development_support[dev] = bool(local_support[local])
        if (lo // chunk) % 20 == 19:
            memory_check(f"{name} correction", limit)
    require(int(support.sum()) > 0, f"{name}: no supported source rows")
    mean_correction = correction_sum / int(support.sum())

    umi_weights = np.divide(
        development_umi,
        development_umi.sum(axis=1, keepdims=True),
        out=np.zeros_like(development_umi),
        where=development_umi.sum(axis=1, keepdims=True) > 0,
    )
    cell_weights = np.divide(
        development_cells,
        development_cells.sum(axis=1, keepdims=True),
        out=np.zeros_like(development_cells, dtype=np.float64),
        where=development_cells.sum(axis=1, keepdims=True) > 0,
    )
    weight_l1 = np.abs(umi_weights - cell_weights).sum(axis=1)

    output.mkdir(parents=True)
    save_npy(output / "development_correction.npy", development_correction)
    save_npy(output / "mean_correction.npy", mean_correction.astype(np.float64))
    save_npy(output / "development_supported.npy", development_support)
    save_npy(output / "development_batch_umi.npy", development_umi)
    save_npy(output / "development_batch_cells.npy", development_cells)
    save_npy(
        output / "batch_control_profiles.npy",
        q.astype(np.float32),
    )
    atomic_json(
        output / "receipt.json",
        {
            "source": name,
            "candidate_effect_values_used": False,
            "expression_access": (
                "K562 streams every X row only to verify ncounts against row sums, "
                "then accumulates expression from control rows only."
                if name == "k562_gwps"
                else "Orion reads exact Non-Targeting expression plus NEGR1 expression "
                "only to reconstruct the frozen legacy control; no other targeting "
                "expression vector is read or used."
            ),
            "candidate_construction_inputs": (
                "control expression plus target batch labels and total-UMI metadata"
            ),
            "batches": len(summary.batch_labels),
            "batch_labels": list(summary.batch_labels),
            "control_policy": summary.control_policy,
            "matched_control_cells": int(summary.matched_control_cells.sum()),
            "matched_control_cells_per_batch": {
                "min": int(summary.matched_control_cells.min()),
                "median": float(np.median(summary.matched_control_cells)),
                "max": int(summary.matched_control_cells.max()),
            },
            "reconstruction_control_cells": int(
                summary.reconstruction_control_cells.sum()
            ),
            "reconstruction_control_cells_per_batch": {
                "min": int(summary.reconstruction_control_cells.min()),
                "median": float(np.median(summary.reconstruction_control_cells)),
                "max": int(summary.reconstruction_control_cells.max()),
            },
            "retained_targets": len(perts),
            "supported_targets": int(support.sum()),
            "development_supported": int(development_support.sum()),
            "global_control_reconstruction": {
                "tolerance": CONTROL_RECON_ATOL,
                "max_abs": float(control_error.max()),
                "mean_abs": float(control_error.mean()),
                "within_tolerance": True,
            },
            "candidate_exact_control_vs_frozen_incumbent": {
                "equality_required": bool(
                    summary.control_policy["pools_identical"]
                ),
                "reason": (
                    "Candidate and reconstruction control pools are identical."
                    if summary.control_policy["pools_identical"]
                    else "Orion frozen b_ctrl includes NEGR1 through the legacy regex."
                ),
                "max_abs": float(matched_vs_existing.max()),
                "mean_abs": float(matched_vs_existing.mean()),
            },
            "target_umi_vs_cell_count_weight_l1": {
                "median": float(np.median(weight_l1)),
                "p95": float(np.quantile(weight_l1, 0.95)),
                "max": float(weight_l1.max()),
            },
            "depth_validation": summary.depth_validation,
            "exact_fallback_rows": int((~development_support).sum()),
            "peak_rss_gib": peak_rss_gib(),
        },
    )
    return read_json(output / "receipt.json")


def directory_manifest(directory: Path, output_name: str = "artifact-manifest.json") -> dict[str, Any]:
    files = {}
    for path in sorted(directory.rglob("*")):
        if path.is_file() and path.name != output_name:
            files[str(path.relative_to(directory))] = file_record(path)
    return {"files": files}


def prepare(limit: float) -> None:
    validate_preoutcome()
    validate_source_hashes(preparation_source_paths())
    require(not PREPARE_DIR.exists(), "source preparation exists")
    require(not IDENTITY_DIR.exists(), "identity artifact exists")
    development = read_lines(DEV_TARGETS)
    responses = read_lines(RESPONSE_GENES)
    PREPARE_DIR.mkdir(parents=True)
    reports = {}
    reports["k562_gwps"] = build_correction(
        "k562_gwps",
        scan_k562(limit),
        development,
        responses,
        PREPARE_DIR / "k562_gwps",
        limit,
    )
    reports["orion_hct116"] = build_correction(
        "orion_hct116",
        scan_orion("orion_hct116", "HCT116", limit),
        development,
        responses,
        PREPARE_DIR / "orion_hct116",
        limit,
    )
    reports["orion_hek293t"] = build_correction(
        "orion_hek293t",
        scan_orion("orion_hek293t", "HEK293T", limit),
        development,
        responses,
        PREPARE_DIR / "orion_hek293t",
        limit,
    )
    atomic_json(
        PREPARE_DIR / "summary.json",
        {
            "schema_version": 1,
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "development_outcomes_opened": False,
            "weighting": "target_total_umi",
            "reports": reports,
            "peak_rss_gib": peak_rss_gib(),
        },
    )
    atomic_json(PREPARE_DIR / "artifact-manifest.json", directory_manifest(PREPARE_DIR))
    print(json.dumps(read_json(PREPARE_DIR / "summary.json"), indent=2, sort_keys=True))


def validate_prepare() -> dict[str, Any]:
    require(PREPARE_DIR.is_dir(), "source preparation missing")
    manifest = read_json(PREPARE_DIR / "artifact-manifest.json")
    for relative, expected in manifest["files"].items():
        path = PREPARE_DIR / relative
        require(path.stat().st_size == expected["bytes"], f"prepared bytes changed: {relative}")
        require(sha256_file(path) == expected["sha256"], f"prepared hash changed: {relative}")
    summary = read_json(PREPARE_DIR / "summary.json")
    require(summary["development_outcomes_opened"] is False, "outcome boundary changed")
    require(summary["weighting"] == "target_total_umi", "weighting changed")
    return summary


@dataclass
class DonorGeometry:
    name: str
    perts: tuple[str, ...]
    genes: tuple[str, ...]
    response_columns: np.ndarray
    canonical_columns: np.ndarray
    mean_response: np.ndarray
    raw_median_norm: float
    development_rows: np.ndarray
    development_global_delta: np.ndarray
    development_global_centered: np.ndarray
    development_correction_centered: np.ndarray
    development_supported: np.ndarray


def blend_rows(
    global_rows: np.ndarray,
    matched_correction: np.ndarray,
    supported: np.ndarray,
    eta: float,
) -> np.ndarray:
    """Blend while making eta zero and unsupported rows literally incumbent."""
    base = np.asarray(global_rows)
    correction = np.asarray(matched_correction)
    support = np.asarray(supported, dtype=bool)
    require(base.shape == correction.shape, "blend shape mismatch")
    require(support.shape == (len(base),), "blend support shape mismatch")
    if eta == 0.0:
        return base.copy()
    output = base.copy()
    output[support] = base[support] + eta * correction[support]
    require(np.array_equal(output[~support], base[~support]), "fallback changed")
    return output


def load_donor_geometry(
    name: str,
    development: Sequence[str],
    responses: Sequence[str],
    canonical: set[str],
    limit: float,
) -> DonorGeometry:
    root = DONOR_ROOTS[name]
    perts = read_lines(root / "perts.txt")
    genes = read_lines(root / "genes.txt", unique=False)
    response_columns = last_columns(genes, responses)
    canonical_columns = np.asarray(
        [index for index, gene in enumerate(genes) if gene in canonical],
        dtype=np.int64,
    )
    target_index = {target: row for row, target in enumerate(perts)}
    development_rows = np.asarray(
        [target_index[target] for target in development], dtype=np.int64
    )
    delta = np.load(root / "delta.npy", mmap_mode="r")
    require(delta.shape == (len(perts), len(genes)), f"{name}: delta shape changed")
    mean_full = np.asarray(delta.mean(axis=0), dtype=np.float32)
    mean_response = mean_full[response_columns]
    norms = np.empty(len(perts), dtype=np.float32)
    chunk = 256
    for lo in range(0, len(perts), chunk):
        hi = min(lo + chunk, len(perts))
        values = np.asarray(delta[lo:hi, canonical_columns], dtype=np.float32)
        values -= mean_full[canonical_columns]
        norms[lo:hi] = np.linalg.norm(values, axis=1)
        if (lo // chunk) % 40 == 39:
            memory_check(f"{name} incumbent norm", limit)
    raw_median_norm = float(np.median(norms))
    require(raw_median_norm > 0.0, f"{name}: zero median norm")
    development_global_delta = np.asarray(
        delta[np.ix_(development_rows, response_columns)], dtype=np.float32
    )
    development_global_centered = development_global_delta - mean_response

    prepared = PREPARE_DIR / name
    correction = np.asarray(
        np.load(prepared / "development_correction.npy"), dtype=np.float64
    )
    mean_correction = np.asarray(
        np.load(prepared / "mean_correction.npy"), dtype=np.float64
    )
    supported = np.asarray(
        np.load(prepared / "development_supported.npy"), dtype=bool
    )
    require(correction.shape == development_global_delta.shape, f"{name}: correction shape")
    require(mean_correction.shape == (len(responses),), f"{name}: correction mean shape")
    require(supported.shape == (len(development),), f"{name}: support shape")
    centered_correction = np.zeros_like(correction)
    centered_correction[supported] = correction[supported] - mean_correction
    require(
        np.array_equal(centered_correction[~supported], np.zeros_like(centered_correction[~supported])),
        f"{name}: unsupported centered correction changed",
    )
    return DonorGeometry(
        name=name,
        perts=perts,
        genes=genes,
        response_columns=response_columns,
        canonical_columns=canonical_columns,
        mean_response=mean_response,
        raw_median_norm=raw_median_norm,
        development_rows=development_rows,
        development_global_delta=development_global_delta,
        development_global_centered=development_global_centered,
        development_correction_centered=centered_correction,
        development_supported=supported,
    )


@dataclass
class Operator:
    evals: np.ndarray
    evecs: np.ndarray

    def apply(self, matrix: np.ndarray) -> np.ndarray:
        values = np.asarray(matrix, dtype=np.float64)
        weights = self.evals ** (-WHITEN_GAMMA / 2.0)
        return ((values @ self.evecs) * weights) @ self.evecs.T


def fit_operator(
    group: str,
    donors: Sequence[str],
    geometry: Mapping[str, DonorGeometry],
    renorm: Mapping[str, float],
    output: Path,
    limit: float,
) -> tuple[Operator, dict[str, Any]]:
    require(not output.exists(), f"operator exists: {output}")
    scratch = IDENTITY_DIR / f".{group}-basis.npy"
    selection = basis_selection(donors)
    selection_report = basis_selection_report(donors)
    total_rows = int(sum(len(selection[name]) for name in donors))
    basis = np.lib.format.open_memmap(
        scratch,
        mode="w+",
        dtype=np.float64,
        shape=(total_rows, RESPONSE_COUNT),
    )
    offset = 0
    chunk = 256
    for name in donors:
        item = geometry[name]
        delta = np.load(DONOR_ROOTS[name] / "delta.npy", mmap_mode="r")
        selected_rows = selection[name]
        for lo in range(0, len(selected_rows), chunk):
            hi = min(lo + chunk, len(selected_rows))
            rows = selected_rows[lo:hi]
            centered = np.asarray(
                delta[rows][:, item.response_columns], dtype=np.float32
            )
            centered -= item.mean_response
            basis[offset + lo : offset + hi] = (
                renorm[name] * centered.astype(np.float64)
            )
            if (lo // chunk) % 40 == 39:
                memory_check(f"{group} whitener basis {name}", limit)
        offset += len(selected_rows)
    require(offset == total_rows, f"{group}: whitener basis row count mismatch")
    basis.flush()
    print(f"[whitener] {group}: covariance from {basis.shape}", flush=True)
    covariance = (basis.T @ basis) / total_rows
    del basis
    evals, evecs = np.linalg.eigh(covariance)
    del covariance
    evals = np.maximum(evals, float(evals.max()) * WHITEN_FLOOR)
    operator = Operator(evals=evals, evecs=evecs)
    np.savez(
        output,
        evals=evals.astype(np.float64),
        evecs=evecs.astype(np.float32),
    )
    atomic_json(output.with_suffix(".basis.json"), selection_report)
    scratch.unlink(missing_ok=True)
    memory_check(f"{group} whitener fit", limit)
    return operator, selection_report


def identity(limit: float) -> None:
    validate_preoutcome()
    validate_prepare()
    validate_source_hashes(identity_source_paths())
    require(not IDENTITY_DIR.exists(), "identity output exists")
    require(not DEVELOPMENT_AUTHORIZATION.exists(), "development already authorized")
    development = read_lines(DEV_TARGETS)
    responses = read_lines(RESPONSE_GENES)
    canonical = set(read_csv_column(CANONICAL_GENES_PATH, "gene_name"))
    geometry = {
        name: load_donor_geometry(name, development, responses, canonical, limit)
        for name in PRODUCTION_DONORS
    }
    require(
        all(np.all(item.development_supported) for item in geometry.values()),
        "all 613 development rows were preregistered as common-source supported",
    )

    IDENTITY_DIR.mkdir(parents=True)
    group_reports = {}
    for group, donors in OUTER_DONORS.items():
        group_key = "k562" if group == "k562" else "non_k562"
        if group_key in group_reports:
            continue
        scales = {name: geometry[name].raw_median_norm for name in donors}
        reference = float(np.median(list(scales.values())))
        renorm = {name: reference / scales[name] for name in donors}
        baseline_rows = [
            renorm[name] * geometry[name].development_global_centered
            for name in donors
        ]
        incumbent_pre = np.mean(np.stack(baseline_rows), axis=0)

        raw_identity_exact = True
        centered_identity_exact = True
        candidate_rows = []
        correction_rows = {}
        for name in donors:
            item = geometry[name]
            eta_zero_raw = blend_rows(
                item.development_global_delta,
                item.development_correction_centered,
                item.development_supported,
                0.0,
            )
            raw_identity_exact &= np.array_equal(
                eta_zero_raw, item.development_global_delta
            )
            eta_zero_centered = eta_zero_raw - item.mean_response
            centered_identity_exact &= np.array_equal(
                eta_zero_centered, item.development_global_centered
            )
            candidate_rows.append(renorm[name] * eta_zero_centered)
            correction_rows[name] = (
                renorm[name]
                * item.development_correction_centered
                / len(donors)
            )
        candidate_pre = np.mean(np.stack(candidate_rows), axis=0)
        pre_exact = np.array_equal(candidate_pre, incumbent_pre)
        require(raw_identity_exact, f"{group_key}: eta-zero raw identity failed")
        require(centered_identity_exact, f"{group_key}: eta-zero centered identity failed")
        require(pre_exact, f"{group_key}: eta-zero pre-whitening identity failed")

        operator_path = IDENTITY_DIR / f"whitener-{group_key}.npz"
        operator, basis_report = fit_operator(
            group_key, donors, geometry, renorm, operator_path, limit
        )
        incumbent_post = operator.apply(incumbent_pre)
        candidate_post = operator.apply(candidate_pre)
        post_exact = np.array_equal(candidate_post, incumbent_post)
        require(post_exact, f"{group_key}: eta-zero post-whitening identity failed")
        save_npy(IDENTITY_DIR / f"baseline-{group_key}.npy", incumbent_post)
        for name, values in correction_rows.items():
            save_npy(
                IDENTITY_DIR / f"correction-{group_key}-{name}.npy",
                operator.apply(values),
            )
        before = float(np.median(np.linalg.norm(incumbent_pre, axis=1)))
        after = float(np.median(np.linalg.norm(incumbent_post, axis=1)))
        group_reports[group_key] = {
            "donors": list(donors),
            "donor_raw_median_norms": scales,
            "donor_renormalization": renorm,
            "whitener_rows": basis_report["rows_kept_total"],
            "whitener_basis_selection": basis_report,
            "whitener_columns": RESPONSE_COUNT,
            "whitener_gamma": WHITEN_GAMMA,
            "preserve_scale": before / (after or 1.0),
            "raw_eta_zero_exact": raw_identity_exact,
            "centered_eta_zero_exact": centered_identity_exact,
            "pre_whitening_eta_zero_exact": pre_exact,
            "post_whitening_eta_zero_exact": post_exact,
            "development_supported_all_donors": True,
        }
        del operator, incumbent_post, candidate_post
        memory_check(f"{group_key} identity", limit)

    result = {
        "schema_version": 1,
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "development_outcomes_opened": False,
        "identity_reproduction_passes": all(
            all(
                report[key]
                for key in (
                    "raw_eta_zero_exact",
                    "centered_eta_zero_exact",
                    "pre_whitening_eta_zero_exact",
                    "post_whitening_eta_zero_exact",
                    "development_supported_all_donors",
                )
            )
            for report in group_reports.values()
        ),
        "groups": group_reports,
        "k562_interpretation": "processing-replicate diagnostic due to a15 PerturbMatch overlap",
        "orion_interpretation": "genuinely new production lever",
        "production_integration_authorized": False,
        "peak_rss_gib": peak_rss_gib(),
    }
    require(result["identity_reproduction_passes"], "identity gate failed")
    atomic_json(IDENTITY_DIR / "result.json", result)
    atomic_json(IDENTITY_DIR / "artifact-manifest.json", directory_manifest(IDENTITY_DIR))
    print(json.dumps(result, indent=2, sort_keys=True))


def validate_identity() -> dict[str, Any]:
    require(IDENTITY_DIR.is_dir(), "identity output missing")
    manifest = read_json(IDENTITY_DIR / "artifact-manifest.json")
    for relative, expected in manifest["files"].items():
        path = IDENTITY_DIR / relative
        require(path.stat().st_size == expected["bytes"], f"identity bytes changed: {relative}")
        require(sha256_file(path) == expected["sha256"], f"identity hash changed: {relative}")
    result = read_json(IDENTITY_DIR / "result.json")
    require(result["identity_reproduction_passes"] is True, "identity reproduction failed")
    require(result["development_outcomes_opened"] is False, "identity opened outcomes")
    return result


def authorize_development(confirm: str | None) -> None:
    require(confirm == DEVELOPMENT_CONFIRM, "explicit development confirmation missing")
    validate_preoutcome()
    validate_prepare()
    validate_identity()
    require(not DEVELOPMENT_AUTHORIZATION.exists(), "development authorization exists")
    require(not DEVELOPMENT_LOCK.exists() and not DEVELOPMENT_DIR.exists(), "development attempted")
    atomic_json(
        DEVELOPMENT_AUTHORIZATION,
        {
            "schema_version": 1,
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "explicit_post_review_confirmation": confirm,
            "development_613_authorized": True,
            "untouched_410_test_authorized": False,
            "preoutcome_authorization_sha256": sha256_file(PREOUTCOME_AUTHORIZATION),
            "source_preparation_manifest_sha256": sha256_file(
                PREPARE_DIR / "artifact-manifest.json"
            ),
            "identity_manifest_sha256": sha256_file(
                IDENTITY_DIR / "artifact-manifest.json"
            ),
            "production_integration_authorized": False,
            "counts_build_submission_authorized": False,
        },
    )


def validate_development_authorization() -> None:
    require(DEVELOPMENT_AUTHORIZATION.is_file(), "development authorization missing")
    value = read_json(DEVELOPMENT_AUTHORIZATION)
    require(value["development_613_authorized"] is True, "development access denied")
    require(value["untouched_410_test_authorized"] is False, "test access unexpectedly enabled")
    require(
        value["preoutcome_authorization_sha256"] == sha256_file(PREOUTCOME_AUTHORIZATION),
        "preoutcome authorization changed",
    )
    require(
        value["source_preparation_manifest_sha256"]
        == sha256_file(PREPARE_DIR / "artifact-manifest.json"),
        "source preparation changed",
    )
    require(
        value["identity_manifest_sha256"]
        == sha256_file(IDENTITY_DIR / "artifact-manifest.json"),
        "identity artifacts changed",
    )


def group_key(group: str) -> str:
    return "k562" if group == "k562" else "non_k562"


def load_geometry_arrays() -> tuple[
    dict[str, np.ndarray], dict[str, dict[str, np.ndarray]]
]:
    baseline = {}
    corrections = {}
    for group, donors in OUTER_DONORS.items():
        key = group_key(group)
        if key not in baseline:
            baseline[key] = np.asarray(
                np.load(IDENTITY_DIR / f"baseline-{key}.npy"), dtype=np.float64
            )
            corrections[key] = {
                name: np.asarray(
                    np.load(IDENTITY_DIR / f"correction-{key}-{name}.npy"),
                    dtype=np.float64,
                )
                for name in donors
            }
    return baseline, corrections


def candidate_matrix(
    base: np.ndarray,
    corrections: Mapping[str, np.ndarray],
    eta: float,
    permutations: Mapping[str, np.ndarray] | None = None,
) -> np.ndarray:
    if eta == 0.0:
        return base.copy()
    output = base.copy()
    for name, values in corrections.items():
        if permutations is None:
            output += eta * values
        else:
            output += eta * values[np.asarray(permutations[name], dtype=np.int64)]
    return output


def pds_subset(
    prediction: np.ndarray,
    truth: np.ndarray,
    query: np.ndarray,
    gallery: np.ndarray,
) -> np.ndarray:
    pred = np.asarray(prediction[query], dtype=np.float64)
    observed = np.asarray(truth[gallery], dtype=np.float64)
    pred /= np.maximum(np.linalg.norm(pred, axis=1, keepdims=True), EPS)
    observed /= np.maximum(np.linalg.norm(observed, axis=1, keepdims=True), EPS)
    similarity = pred @ observed.T
    positions = {int(value): index for index, value in enumerate(gallery)}
    own_columns = np.asarray([positions[int(value)] for value in query], dtype=np.int64)
    own = similarity[np.arange(len(query)), own_columns]
    greater = np.sum(similarity > own[:, None], axis=1)
    equal = np.sum(similarity == own[:, None], axis=1)
    midrank = 1.0 + greater + 0.5 * (equal - 1)
    return 1.0 - (midrank - 1.0) / (len(gallery) - 1.0)


def signed_credit(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    numerator = np.sum(np.sign(truth) * prediction, axis=1)
    denominator = np.sum(np.abs(prediction), axis=1)
    return np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator),
        where=denominator > EPS,
    )


def score_grid(
    baseline: Mapping[str, np.ndarray],
    corrections: Mapping[str, Mapping[str, np.ndarray]],
    truth: Mapping[str, np.ndarray],
    permutations: Mapping[str, np.ndarray] | None = None,
) -> tuple[dict[float, dict[str, dict[str, np.ndarray]]], dict[float, dict[str, np.ndarray]]]:
    all_indices = np.arange(DEV_COUNT, dtype=np.int64)
    identity_prediction = {
        group: baseline[group_key(group)] for group in OUTER_DONORS
    }
    identity_pds = {
        dataset: pds_subset(
            identity_prediction[DATASET_GROUP[dataset]],
            truth[dataset],
            all_indices,
            all_indices,
        )
        for dataset in SCORED_DATASETS
    }
    identity_credit = {
        dataset: signed_credit(
            identity_prediction[DATASET_GROUP[dataset]], truth[dataset]
        )
        for dataset in SCORED_DATASETS
    }
    scores: dict[float, dict[str, dict[str, np.ndarray]]] = {}
    predictions: dict[float, dict[str, np.ndarray]] = {}
    for eta in ETAS:
        predictions[eta] = {}
        for group in OUTER_DONORS:
            key = group_key(group)
            group_permutations = None
            if permutations is not None:
                group_permutations = {
                    name: permutations[name] for name in corrections[key]
                }
            predictions[eta][group] = candidate_matrix(
                baseline[key], corrections[key], eta, group_permutations
            )
        scores[eta] = {}
        for dataset in SCORED_DATASETS:
            group = DATASET_GROUP[dataset]
            scores[eta][dataset] = {
                "pds": pds_subset(
                    predictions[eta][group],
                    truth[dataset],
                    all_indices,
                    all_indices,
                )
                - identity_pds[dataset],
                "signed_credit": signed_credit(
                    predictions[eta][group], truth[dataset]
                )
                - identity_credit[dataset],
            }
    return scores, predictions


def aggregate(
    score: Mapping[str, Mapping[str, np.ndarray]], indices: np.ndarray
) -> dict[str, Any]:
    idx = np.asarray(indices, dtype=np.int64)
    datasets = {
        dataset: {
            metric: float(np.mean(score[dataset][metric][idx]))
            for metric in ("pds", "signed_credit")
        }
        for dataset in SCORED_DATASETS
    }
    families = {}
    per_target = {}
    macro = {}
    for metric in ("pds", "signed_credit"):
        family_arrays = {
            "k562": 0.5
            * (
                score["k562_gwps"][metric][idx]
                + score["k562_essential"][metric][idx]
            ),
            "rpe1": score["rpe1"][metric][idx],
            "jurkat": score["nadig_jurkat"][metric][idx],
        }
        families[metric] = {
            name: float(values.mean()) for name, values in family_arrays.items()
        }
        per_target[metric] = np.mean(np.stack(list(family_arrays.values())), axis=0)
        macro[metric] = float(per_target[metric].mean())
    return {
        "datasets": datasets,
        "families": families,
        "per_target": per_target,
        "macro": macro,
    }


def select_eta(
    scores: Mapping[float, Mapping[str, Mapping[str, np.ndarray]]],
    indices: np.ndarray,
    *,
    predictions: Mapping[float, Mapping[str, np.ndarray]] | None = None,
    truth: Mapping[str, np.ndarray] | None = None,
) -> tuple[float, dict[float, dict[str, Any]]]:
    idx = np.asarray(indices, dtype=np.int64)
    summaries: dict[float, dict[str, Any]] = {}
    for eta in ETAS:
        if predictions is None:
            summaries[eta] = aggregate(scores[eta], idx)
            continue
        require(truth is not None, "truth required for inner-gallery selection")
        local_score = {}
        for dataset in SCORED_DATASETS:
            group = DATASET_GROUP[dataset]
            candidate = pds_subset(
                predictions[eta][group], truth[dataset], idx, idx
            )
            identity = pds_subset(
                predictions[0.0][group], truth[dataset], idx, idx
            )
            local_score[dataset] = {
                "pds": candidate - identity,
                "signed_credit": scores[eta][dataset]["signed_credit"][idx],
            }
        summaries[eta] = aggregate(local_score, np.arange(len(idx), dtype=np.int64))
    eligible = [
        eta
        for eta in ETAS
        if all(
            value >= -EPS
            for metric in ("pds", "signed_credit")
            for value in summaries[eta]["families"][metric].values()
        )
    ]
    selected = 0.0
    for eta in eligible:
        candidate = summaries[eta]["macro"]
        incumbent = summaries[selected]["macro"]
        key = (candidate["pds"], candidate["signed_credit"], -eta)
        old = (incumbent["pds"], incumbent["signed_credit"], -selected)
        if key > old:
            selected = eta
    return selected, summaries


def nested_selection(
    scores: Mapping[float, Mapping[str, Mapping[str, np.ndarray]]],
    predictions: Mapping[float, Mapping[str, np.ndarray]],
    truth: Mapping[str, np.ndarray],
    folds: np.ndarray,
) -> dict[str, Any]:
    all_indices = np.arange(DEV_COUNT, dtype=np.int64)
    oof = {
        dataset: {
            metric: np.empty(DEV_COUNT, dtype=np.float64)
            for metric in ("pds", "signed_credit")
        }
        for dataset in SCORED_DATASETS
    }
    records = []
    for fold in range(FOLDS):
        held = all_indices[folds == fold]
        inner = all_indices[folds != fold]
        selected, summaries = select_eta(
            scores, inner, predictions=predictions, truth=truth
        )
        for dataset in SCORED_DATASETS:
            for metric in ("pds", "signed_credit"):
                oof[dataset][metric][held] = scores[selected][dataset][metric][held]
        records.append(
            {
                "fold": fold,
                "selected_eta": selected,
                "inner_macro": summaries[selected]["macro"],
                "held_complete_gallery": {
                    key: value
                    for key, value in aggregate(scores[selected], held).items()
                    if key != "per_target"
                },
            }
        )
    final_eta, summaries = select_eta(scores, all_indices)
    return {
        "folds": records,
        "oof": oof,
        "oof_summary": aggregate(oof, all_indices),
        "final_eta": final_eta,
        "final_summary": summaries[final_eta],
    }


def block_bootstrap(
    values: np.ndarray, labels: np.ndarray, draws: np.ndarray
) -> np.ndarray:
    labels = np.asarray(labels, dtype=np.int64)
    block_count = int(labels.max()) + 1
    sums = np.bincount(labels, weights=np.asarray(values, dtype=np.float64), minlength=block_count)
    sizes = np.bincount(labels, minlength=block_count).astype(np.float64)
    return sums[draws].sum(axis=1) / sizes[draws].sum(axis=1)


def positive_share(values: Iterable[float]) -> float:
    positive = np.maximum(np.asarray(tuple(values), dtype=np.float64), 0.0)
    return 1.0 if positive.sum() <= EPS else float(positive.max() / positive.sum())


def concentration(
    score: Mapping[str, Mapping[str, np.ndarray]], blocks: np.ndarray
) -> dict[str, Any]:
    summary = aggregate(score, np.arange(DEV_COUNT, dtype=np.int64))
    result = {}
    for metric in ("pds", "signed_credit"):
        block_values = np.bincount(
            blocks.astype(np.int64),
            weights=summary["per_target"][metric],
        )
        result[metric] = {
            "family": positive_share(summary["families"][metric].values()),
            "target_block": positive_share(block_values),
        }
    return result


def make_permutations(
    kind: str, draw: int, arrays: Mapping[str, np.ndarray]
) -> dict[str, np.ndarray]:
    if kind == "target_label":
        common = arrays["target_label_nulls"][draw]
        return {name: common for name in PRODUCTION_DONORS}
    require(kind == "batch_mixture", f"unknown null kind {kind}")
    values = arrays["batch_mixture_nulls"][draw]
    return {
        name: values[index] for index, name in enumerate(PRODUCTION_DONORS)
    }


def run_nulls(
    kind: str,
    baseline: Mapping[str, np.ndarray],
    corrections: Mapping[str, Mapping[str, np.ndarray]],
    truth: Mapping[str, np.ndarray],
    folds: np.ndarray,
    arrays: Mapping[str, np.ndarray],
    limit: float,
) -> dict[str, dict[str, np.ndarray]]:
    output = {
        scope: {
            metric: np.empty(NULL_DRAWS, dtype=np.float64)
            for metric in ("pds", "signed_credit")
        }
        for scope in ("nested", "fixed")
    }
    for draw in range(NULL_DRAWS):
        scores, predictions = score_grid(
            baseline,
            corrections,
            truth,
            make_permutations(kind, draw, arrays),
        )
        result = nested_selection(scores, predictions, truth, folds)
        for metric in ("pds", "signed_credit"):
            output["nested"][metric][draw] = result["oof_summary"]["macro"][metric]
            output["fixed"][metric][draw] = result["final_summary"]["macro"][metric]
        if draw % 5 == 4:
            print(f"[{kind}] {draw + 1}/{NULL_DRAWS}", flush=True)
            memory_check(f"{kind} null", limit)
    return output


def load_development_truth(limit: float) -> dict[str, np.ndarray]:
    targets = read_lines(DEV_TARGETS)
    responses = read_lines(RESPONSE_GENES)
    banks = PARENT.source_banks()
    truth = {}
    for dataset in SCORED_DATASETS:
        truth[dataset] = np.asarray(
            PARENT.load_effect_rows(banks[dataset], targets, responses),
            dtype=np.float64,
        )
        memory_check(f"loaded development truth {dataset}", limit)
    return truth


def bootstrap_report(
    score: Mapping[str, Mapping[str, np.ndarray]],
    blocks: np.ndarray,
    draws: np.ndarray,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    summary = aggregate(score, np.arange(DEV_COUNT, dtype=np.int64))
    report, values = {}, {}
    for metric in ("pds", "signed_credit"):
        samples = block_bootstrap(summary["per_target"][metric], blocks, draws)
        values[metric] = samples
        report[metric] = {
            "draws": len(samples),
            "lower_95": float(np.quantile(samples, 0.025, method="linear")),
            "upper_95": float(np.quantile(samples, 0.975, method="linear")),
        }
    return report, values


def result_checks(
    nested: Mapping[str, Any],
    fixed: Mapping[str, Any],
    boot: Mapping[str, Mapping[str, Any]],
    nulls: Mapping[str, Mapping[str, Any]],
    concentration_report: Mapping[str, Mapping[str, Any]],
    identity_result: Mapping[str, Any],
) -> dict[str, bool]:
    checks = {
        "nonidentity": float(nested["final_eta"]) > 0.0,
        "identity_exact": identity_result["identity_reproduction_passes"] is True,
    }
    for scope, summary in (("nested", nested["oof_summary"]), ("fixed", fixed)):
        for metric in ("pds", "signed_credit"):
            checks[f"{scope}_{metric}_families_nonnegative"] = all(
                value >= -EPS for value in summary["families"][metric].values()
            )
            checks[f"{scope}_{metric}_macro_at_least_0.003"] = (
                summary["macro"][metric] >= MIN_GAIN
            )
            checks[f"{scope}_bootstrap_{metric}_lower_positive"] = (
                boot[scope][metric]["lower_95"] > 0.0
            )
            checks[f"{scope}_target_label_null_{metric}"] = (
                summary["macro"][metric]
                > nulls["target_label"][scope][metric]["p95"]
            )
            checks[f"{scope}_batch_mixture_null_{metric}"] = (
                summary["macro"][metric]
                > nulls["batch_mixture"][scope][metric]["p95"]
            )
            checks[f"{scope}_family_concentration_{metric}"] = (
                concentration_report[scope][metric]["family"] <= MAX_CONCENTRATION
            )
            checks[f"{scope}_target_block_concentration_{metric}"] = (
                concentration_report[scope][metric]["target_block"] <= MAX_CONCENTRATION
            )
    return checks


def component_corrections(
    corrections: Mapping[str, Mapping[str, np.ndarray]],
    mode: str,
) -> dict[str, dict[str, np.ndarray]]:
    require(mode in COMPONENT_MODES, f"unknown component mode {mode}")
    if mode == "all":
        return {
            group: {name: values for name, values in donor_values.items()}
            for group, donor_values in corrections.items()
        }
    keep = (
        {"orion_hct116", "orion_hek293t"}
        if mode == "orion_only"
        else {"k562_gwps"}
    )
    return {
        group: {
            name: values
            for name, values in donor_values.items()
            if name in keep
        }
        for group, donor_values in corrections.items()
    }


def apply_primary_selection(
    scores: Mapping[float, Mapping[str, Mapping[str, np.ndarray]]],
    folds: np.ndarray,
    nested: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    all_indices = np.arange(DEV_COUNT, dtype=np.int64)
    oof = {
        dataset: {
            metric: np.empty(DEV_COUNT, dtype=np.float64)
            for metric in ("pds", "signed_credit")
        }
        for dataset in SCORED_DATASETS
    }
    for record in nested["folds"]:
        held = all_indices[folds == int(record["fold"])]
        eta = float(record["selected_eta"])
        for dataset in SCORED_DATASETS:
            for metric in ("pds", "signed_credit"):
                oof[dataset][metric][held] = scores[eta][dataset][metric][held]
    fixed_score = scores[float(nested["final_eta"])]
    return aggregate(oof, all_indices), aggregate(fixed_score, all_indices)


def develop(limit: float) -> None:
    validate_preoutcome()
    validate_prepare()
    identity_result = validate_identity()
    validate_development_authorization()
    validate_source_hashes(development_source_paths())
    require(not DEVELOPMENT_LOCK.exists() and not DEVELOPMENT_DIR.exists(), "development attempted")
    _, arrays = validate_plan()
    atomic_json(
        DEVELOPMENT_LOCK,
        {
            "created_utc": dt.datetime.now(dt.UTC).isoformat(),
            "development_targets": DEV_COUNT,
            "untouched_test_effects_accessed": False,
            "authorization_sha256": sha256_file(DEVELOPMENT_AUTHORIZATION),
        },
    )
    DEVELOPMENT_DIR.mkdir(parents=True)
    truth = load_development_truth(limit)
    baseline, corrections = load_geometry_arrays()
    scores, predictions = score_grid(baseline, corrections, truth)
    nested = nested_selection(
        scores, predictions, truth, arrays["development_folds"]
    )
    fixed = aggregate(
        scores[float(nested["final_eta"])],
        np.arange(DEV_COUNT, dtype=np.int64),
    )
    selected_score = nested["oof"]
    fixed_score = scores[float(nested["final_eta"])]
    boot = {}
    boot_arrays = {}
    concentration_report = {}
    for scope, selected in (("nested", selected_score), ("fixed", fixed_score)):
        boot[scope], boot_arrays[scope] = bootstrap_report(
            selected,
            arrays["development_blocks"],
            arrays["bootstrap_block_draws"].astype(np.int64),
        )
        concentration_report[scope] = concentration(
            selected, arrays["development_blocks"]
        )
    null_values = {
        kind: run_nulls(
            kind,
            baseline,
            corrections,
            truth,
            arrays["development_folds"],
            arrays,
            limit,
        )
        for kind in ("target_label", "batch_mixture")
    }
    null_report = {
        kind: {
            scope: {
                metric: {
                    "draws": len(values[scope][metric]),
                    "p95": float(
                        np.quantile(values[scope][metric], 0.95, method="linear")
                    ),
                    "add_one_p": float(
                        (
                            1
                            + np.sum(
                                values[scope][metric]
                                >= (
                                    nested["oof_summary"]["macro"][metric]
                                    if scope == "nested"
                                    else fixed["macro"][metric]
                                )
                            )
                        )
                        / (len(values[scope][metric]) + 1)
                    ),
                }
                for metric in ("pds", "signed_credit")
            }
            for scope in ("nested", "fixed")
        }
        for kind, values in null_values.items()
    }
    component_diagnostics = {}
    for mode in COMPONENT_MODES:
        mode_scores, _ = score_grid(
            baseline,
            component_corrections(corrections, mode),
            truth,
        )
        mode_oof, mode_fixed = apply_primary_selection(
            mode_scores,
            arrays["development_folds"],
            nested,
        )
        component_diagnostics[mode] = {
            "uses_primary_selected_etas_without_retuning": True,
            "nested_oof": {
                key: value for key, value in mode_oof.items() if key != "per_target"
            },
            "fixed_final_candidate": {
                key: value for key, value in mode_fixed.items() if key != "per_target"
            },
        }
    checks = result_checks(
        nested, fixed, boot, null_report, concentration_report, identity_result
    )
    result = {
        "schema_version": 1,
        "created_utc": dt.datetime.now(dt.UTC).isoformat(),
        "stage": "source_only_613_target_development",
        "selected_eta": nested["final_eta"],
        "nested_folds": nested["folds"],
        "nested_oof": {
            key: value
            for key, value in nested["oof_summary"].items()
            if key != "per_target"
        },
        "fixed_final_candidate": {
            key: value for key, value in fixed.items() if key != "per_target"
        },
        "bootstrap": boot,
        "nulls": null_report,
        "positive_gain_concentration": concentration_report,
        "component_diagnostics": component_diagnostics,
        "checks": checks,
        "passes": all(checks.values()),
        "decision": "PASS_FOR_SEPARATE_REVIEW_ONLY" if all(checks.values()) else "NO-GO",
        "k562_interpretation": "processing-replicate diagnostic with PerturbMatch overlap",
        "orion_interpretation": "genuinely new production lever",
        "a15_integration_authorized": False,
        "untouched_410_test_accessed": False,
        "peak_rss_gib": peak_rss_gib(),
    }
    atomic_json(DEVELOPMENT_DIR / "result.json", result)
    arrays_to_save = {}
    for dataset in SCORED_DATASETS:
        for metric in ("pds", "signed_credit"):
            arrays_to_save[f"oof__{dataset}__{metric}"] = selected_score[dataset][metric]
    for kind, values in null_values.items():
        for scope, metrics in values.items():
            for metric, array in metrics.items():
                arrays_to_save[f"null__{kind}__{scope}__{metric}"] = array
    for scope, metrics in boot_arrays.items():
        for metric, array in metrics.items():
            arrays_to_save[f"bootstrap__{scope}__{metric}"] = array
    atomic_npz(DEVELOPMENT_DIR / "per-target-nulls-bootstrap.npz", **arrays_to_save)
    atomic_json(
        DEVELOPMENT_DIR / "artifact-manifest.json",
        directory_manifest(DEVELOPMENT_DIR),
    )
    print(json.dumps(result, indent=2, sort_keys=True))


def audit() -> None:
    payload: dict[str, Any] = {
        "plan_exists": PLAN_DIR.exists(),
        "source_hash_manifest_exists": SOURCE_HASH_MANIFEST.exists(),
        "preoutcome_authorization_exists": PREOUTCOME_AUTHORIZATION.exists(),
        "source_preparation_exists": PREPARE_DIR.exists(),
        "identity_reproduction_exists": IDENTITY_DIR.exists(),
        "development_authorized": DEVELOPMENT_AUTHORIZATION.exists(),
        "development_accessed": DEVELOPMENT_LOCK.exists(),
        "development_result_exists": DEVELOPMENT_DIR.exists(),
        "untouched_test_code_path_exists": False,
    }
    if PREOUTCOME_AUTHORIZATION.exists():
        payload["preoutcome"] = validate_preoutcome()
    if PREPARE_DIR.exists():
        payload["source_preparation"] = validate_prepare()
    if IDENTITY_DIR.exists():
        payload["identity"] = validate_identity()
    if DEVELOPMENT_AUTHORIZATION.exists():
        validate_development_authorization()
    if DEVELOPMENT_DIR.exists():
        manifest = read_json(DEVELOPMENT_DIR / "artifact-manifest.json")
        for relative, expected in manifest["files"].items():
            path = DEVELOPMENT_DIR / relative
            require(path.stat().st_size == expected["bytes"], f"development bytes changed: {relative}")
            require(sha256_file(path) == expected["sha256"], f"development hash changed: {relative}")
        payload["development"] = read_json(DEVELOPMENT_DIR / "result.json")
    print(json.dumps(payload, indent=2, sort_keys=True))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "command",
        choices=(
            "plan",
            "freeze",
            "prepare",
            "identity",
            "authorize-development",
            "develop",
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
    elif args.command == "prepare":
        prepare(args.max_rss_gib)
    elif args.command == "identity":
        identity(args.max_rss_gib)
    elif args.command == "authorize-development":
        authorize_development(args.confirm)
    elif args.command == "develop":
        develop(args.max_rss_gib)
    else:
        audit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

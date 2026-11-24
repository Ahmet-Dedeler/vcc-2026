#!/usr/bin/env python3
"""Experiment 142: batch-disjoint donor reliability weighting.

Only ``plan`` is lightweight. ``freeze`` binds the pre-outcome implementation
and source provenance. ``source-build`` and ``source-evaluate`` are explicit
per-donor heavy stages. ``source-gate`` combines source-only evidence.

The ``develop`` command is a fail-closed interface boundary: it cannot open
outcomes without a separately added authorization artifact, and this
pre-outcome implementation contains no held-outcome loader.
"""

from __future__ import annotations

import argparse
import datetime as dt
import gc
import hashlib
import json
import os
import resource
import shutil
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator, Mapping, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")

import h5py
import numpy as np
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse as sp
from scipy.stats import rankdata


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
README = HERE / "README.md"
PROTOCOL = HERE / "preregistered_protocol.json"
TESTS = HERE / "test_run.py"
DEVELOPMENT_INTERFACE = HERE / "DEVELOPMENT_INTERFACE.v1.json"
PREFREEZE_AMENDMENT = HERE / "PREFREEZE_AMENDMENT_001.json"
PREFREEZE_AMENDMENT_2 = HERE / "PREFREEZE_AMENDMENT_002.json"
PREFREEZE_AMENDMENT_3 = HERE / "PREFREEZE_AMENDMENT_003.json"
PREFREEZE_AMENDMENT_4 = HERE / "PREFREEZE_AMENDMENT_004.json"
PREFREEZE_AMENDMENT_5 = HERE / "PREFREEZE_AMENDMENT_005.json"
PREFREEZE_AMENDMENT_6 = HERE / "PREFREEZE_AMENDMENT_006.json"
PROVENANCE_CLOCK_CORRECTION = HERE / "PROVENANCE_CLOCK_CORRECTION.v2.json"
SUPERSEDED_CLOCK_CORRECTION = HERE / "PROVENANCE_CLOCK_CORRECTION.v1.json"
SOURCE_ONLY_FREEZE_AUTHORIZATION = HERE / "SOURCE_ONLY_FREEZE_AUTHORIZATION.v6.json"
SUPERSEDED_SOURCE_AUTHORIZATION_V5 = HERE / "SOURCE_ONLY_FREEZE_AUTHORIZATION.v5.json"
SUPERSEDED_SOURCE_AUTHORIZATION_V4 = HERE / "SOURCE_ONLY_FREEZE_AUTHORIZATION.v4.json"
SUPERSEDED_SOURCE_AUTHORIZATION_V3 = HERE / "SOURCE_ONLY_FREEZE_AUTHORIZATION.v3.json"
SUPERSEDED_SOURCE_AUTHORIZATION_V2 = HERE / "SOURCE_ONLY_FREEZE_AUTHORIZATION.v2.json"

EXP132 = ROOT / "experiments/132-response-decomposition-ridge"
EXP132_PLAN = EXP132 / "work/preregistered-freeze-v2"
DEV_TARGETS = EXP132_PLAN / "selection_targets.txt"
RESPONSE_GENES = EXP132_PLAN / "response_genes.txt"
EXP132_SOURCE_AUTHORIZATION = EXP132 / "SOURCE_ONLY_AUTHORIZATION.v2.json"
EXP132_CANDIDATE_FREEZE = EXP132 / "SOURCE_CANDIDATE_FREEZE.v2.json"
EXP132_PLAN_MANIFEST = EXP132_PLAN / "plan-manifest.json"
EXP132_SELECTION_MANIFEST = EXP132 / "work/selection-v2/artifact-manifest.json"

EXP140 = ROOT / "experiments/140-batch-matched-controls"
EXP140_RUN = EXP140 / "run.py"
EXP140_PLAN = EXP140 / "work/preregistered-plan-v2"
EXP140_PREPARE = EXP140 / "work/source-preparation-v5"
EXP140_IDENTITY = EXP140 / "work/identity-reproduction-v5"
EXP140_SOURCE_MANIFEST = EXP140 / "SOURCE_HASH_MANIFEST.v5.json"
EXP140_PREOUTCOME = EXP140 / "PREOUTCOME_AUTHORIZATION.v5.json"

PLAN_DIR = HERE / "work/preregistered-plan-v2"
FREEZE_DIR = HERE / "work/preoutcome-freeze-v4"
SOURCE_HASH_MANIFEST = FREEZE_DIR / "SOURCE_HASH_MANIFEST.v1.json"
PREOUTCOME_AUTHORIZATION = FREEZE_DIR / "PREOUTCOME_AUTHORIZATION.v1.json"
SOURCE_BUILD_DIR = HERE / "work/source-build-v4"
SOURCE_EVAL_DIR = HERE / "work/source-evaluation-v4"
SOURCE_GATE_DIR = HERE / "work/source-gate-v4"
DEVELOPMENT_AUTHORIZATION = HERE / "DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json"
DEVELOPMENT_LOCK = HERE / "work/development-outcome-access.v1.lock"
HEAVY_STAGE_LOCK = HERE / "work/heavy-stage.active.lock"
SUPERSEDED_PREOUTCOME_AUTHORIZATION_V1 = (
    HERE / "work/preoutcome-freeze-v1/PREOUTCOME_AUTHORIZATION.v1.json"
)
SUPERSEDED_K562_BUILD_MANIFEST_V1 = (
    HERE / "work/source-build-v1/k562_gwps/artifact-manifest.json"
)
SUPERSEDED_PREOUTCOME_AUTHORIZATION_V2 = (
    HERE / "work/preoutcome-freeze-v2/PREOUTCOME_AUTHORIZATION.v1.json"
)
SUPERSEDED_K562_BUILD_MANIFEST_V2 = (
    HERE / "work/source-build-v2/k562_gwps/artifact-manifest.json"
)
SUPERSEDED_K562_EVALUATION_MANIFEST_V2 = (
    HERE / "work/source-evaluation-v2/k562_gwps/artifact-manifest.json"
)
SUPERSEDED_PREOUTCOME_AUTHORIZATION_V3 = (
    HERE / "work/preoutcome-freeze-v3/PREOUTCOME_AUTHORIZATION.v1.json"
)
SUPERSEDED_K562_BUILD_MANIFEST_V3 = (
    HERE / "work/source-build-v3/k562_gwps/artifact-manifest.json"
)
SUPERSEDED_K562_EVALUATION_MANIFEST_V3 = (
    HERE / "work/source-evaluation-v3/k562_gwps/artifact-manifest.json"
)

K562_RAW = ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
ORION_DATA = ROOT / "data/external/orion/data"
ORION_GENE_METADATA = ROOT / "data/external/orion/metadata/gene_metadata.parquet"
DONOR_ROOTS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
DONORS = tuple(DONOR_ROOTS)
DONOR_LINES = {
    "orion_hct116": "HCT116",
    "orion_hek293t": "HEK293T",
}
OPERATOR_ROLES = {
    "k562_gwps": ("non_k562",),
    "orion_hct116": ("non_k562", "k562"),
    "orion_hek293t": ("non_k562", "k562"),
}
OUTER_DONORS = {
    "non_k562": DONORS,
    "k562": ("orion_hct116", "orion_hek293t"),
}

EXPECTED = {
    DEV_TARGETS: "74a070a538455ffc2a528dae3094d21e047905463779b2a5206d681f5ca6ebb5",
    RESPONSE_GENES: "24743754eead5b8ddd81634bd12feac1df431eb2b9581197755b75c4c90949c8",
    EXP132_SOURCE_AUTHORIZATION: "9bca6e620359e44b5440d0b8588412f40f274304ae67d9f55fd7e7d03341e197",
    EXP132_CANDIDATE_FREEZE: "d1b14368f0aab9ee76798b3f298441af8c50cb6d0b5a8f5e8865e70953acc373",
    EXP132_PLAN_MANIFEST: "170d55ed5f812eb85f5503af4bb9d7dec8d17dce6e58fe966d85bf0eec5185d1",
    EXP132_SELECTION_MANIFEST: "822f624f4de6c6e5b14a0fefe16fbe865d343c17fcf14f52fe854acf19e0b7ed",
    EXP140_RUN: "3bea7cb91e7b5d86cb41d2921b156de164e79d9c3dc0837c20f416c4d3fc1a93",
    EXP140_PLAN
    / "plan-manifest.json": "2f91872b1318ab7f0ee608600f5acfaaad0397f508e833707e28c2f379e6007c",
    EXP140_PLAN
    / "plan-arrays.npz": "bdd17742d8691ff5174631c11c00c7d3dfb23f51c9e225ca672549608d466d5b",
    EXP140_PREPARE
    / "artifact-manifest.json": "90439ff15cbb18d2e283fdee247ba5cb75d471dd07930a8ce49d9448367df4b2",
    EXP140_IDENTITY
    / "artifact-manifest.json": "1967a66f463ba78b02c674c09274c9c399fb68b405ff9b0deb045ec4772a5649",
    EXP140_SOURCE_MANIFEST: "70c860e308bdf433fbd8e6812e95c5e24b090173f2cc3ac4e6d689ab9195a149",
    EXP140_PREOUTCOME: "395899b2a395fbbc9fec35f7a81fa20704eff95e83c84ba09b8ff8ec73314703",
}

DEV_COUNT = 613
RESPONSE_COUNT = 4992
SPLITS = 8
HALVES = 2
SPLIT_SEEDS = tuple(range(142_202_610, 142_202_618))
ETAS = (0.0, 0.25, 0.50)
MIN_CELLS_HALF = 5
MIN_BATCHES_HALF = 4
MIN_VALID_PAIRS = 6
ALL_EIGHT_FRACTION = 0.95
GLOBAL_SPLIT_SEARCH = 512
BOOTSTRAP_DRAWS = 20_000
BOOTSTRAP_SEED = 142_202_620
NULL_DRAWS = 49
NULL_SEED = 142_202_621
DONOR_NULL_SEED = 142_202_622
MATCH_BLOCK_SIZE = 8
TS = 50_000.0
WHITEN_GAMMA = 0.65
EPS = 1e-12
GIB = 1024**3
POOLED_RECON_ATOL = 2e-5
ORION_SCAN_BATCH_SIZE = 256
AUTHORIZED_ARTIFACT_DATE_UTC = dt.date(2026, 8, 31)


class ContractError(RuntimeError):
    """Raised when a frozen experimental boundary is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if sys.platform == "darwin" else value / 1024**2


def provenance_utc() -> str:
    """Return an explicit August 31 UTC timestamp; never trust the host clock."""
    override = os.environ.get("VCC_PROVENANCE_UTC")
    require(
        override is not None,
        "artifact stage requires an explicit VCC_PROVENANCE_UTC override",
    )
    parsed = dt.datetime.fromisoformat(override.replace("Z", "+00:00"))
    require(parsed.tzinfo is not None, "VCC_PROVENANCE_UTC must include a timezone")
    parsed_utc = parsed.astimezone(dt.UTC)
    require(
        parsed_utc.date() == AUTHORIZED_ARTIFACT_DATE_UTC,
        "VCC_PROVENANCE_UTC must resolve to the authorized August 31 UTC date",
    )
    return parsed_utc.isoformat()


def authorized_freeze_utc(source_authorization: Mapping[str, Any]) -> str:
    """Require the exact corrected timestamp bound by the source authorization."""
    override = os.environ.get("VCC_PROVENANCE_UTC")
    require(
        override is not None,
        "freeze requires the authorized VCC_PROVENANCE_UTC override",
    )
    actual = dt.datetime.fromisoformat(override.replace("Z", "+00:00"))
    require(actual.tzinfo is not None, "VCC_PROVENANCE_UTC must include a timezone")
    expected_text = source_authorization["bindings"]["freeze_provenance_utc"]
    expected = dt.datetime.fromisoformat(expected_text.replace("Z", "+00:00"))
    require(expected.tzinfo is not None, "authorized freeze timestamp lacks timezone")
    actual_utc = actual.astimezone(dt.UTC).isoformat()
    expected_utc = expected.astimezone(dt.UTC).isoformat()
    require(
        actual_utc == expected_utc,
        "VCC_PROVENANCE_UTC does not match the authorized freeze timestamp",
    )
    return actual_utc


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


def sha256_array(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def read_lines(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    require(bool(values), f"empty axis: {path}")
    if unique:
        require(len(values) == len(set(values)), f"duplicate axis: {path}")
    return values


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def file_record(path: Path, *, hash_file: bool = True) -> dict[str, Any]:
    record: dict[str, Any] = {"bytes": path.stat().st_size}
    if hash_file:
        record["sha256"] = sha256_file(path)
    return record


def atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    require(not path.exists(), f"refuse overwrite: {path}")
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
    require(not path.exists(), f"refuse overwrite: {path}")
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
    require(not path.exists(), f"refuse overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as handle:
        np.save(handle, array)
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


@contextmanager
def atomic_stage_directory(path: Path) -> Iterator[Path]:
    """Publish a complete stage directory with one same-filesystem rename."""
    require(not path.exists(), f"refuse overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{path.name}.", dir=path.parent))
    try:
        yield temporary
        os.replace(temporary, path)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise


@contextmanager
def heavy_stage_lock(stage: str, donor: str | None) -> Iterator[None]:
    """Serialize every source build, evaluation, and gate across the experiment."""
    require(
        stage in ("source-build", "source-evaluate", "source-gate"),
        f"invalid heavy stage: {stage}",
    )
    HEAVY_STAGE_LOCK.parent.mkdir(parents=True, exist_ok=True)
    try:
        descriptor = os.open(
            HEAVY_STAGE_LOCK,
            os.O_CREAT | os.O_EXCL | os.O_WRONLY,
            0o600,
        )
    except FileExistsError as error:
        raise ContractError(
            f"another Experiment 142 heavy stage is active: {HEAVY_STAGE_LOCK}"
        ) from error
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(f"{stage}\n{donor or '-'}\n{os.getpid()}\n")
            handle.flush()
            os.fsync(handle.fileno())
        yield
    finally:
        HEAVY_STAGE_LOCK.unlink(missing_ok=True)


def directory_manifest(
    directory: Path, output_name: str = "artifact-manifest.json"
) -> dict[str, Any]:
    files = {}
    for path in sorted(directory.rglob("*")):
        if path.is_file() and path.name != output_name:
            files[str(path.relative_to(directory))] = file_record(path)
    return {"schema_version": 1, "files": files}


def validate_directory_manifest(
    directory: Path, manifest_name: str = "artifact-manifest.json"
) -> dict[str, Any]:
    manifest_path = directory / manifest_name
    require(manifest_path.is_file(), f"missing artifact manifest: {manifest_path}")
    manifest = read_json(manifest_path)
    for relative, expected in manifest["files"].items():
        path = directory / relative
        require(path.is_file(), f"missing artifact: {path}")
        require(path.stat().st_size == expected["bytes"], f"bytes changed: {path}")
        require(sha256_file(path) == expected["sha256"], f"hash changed: {path}")
    return manifest


def validate_pinned_inputs() -> None:
    for path, expected in EXPECTED.items():
        require(path.is_file(), f"missing pinned input: {path}")
        require(sha256_file(path) == expected, f"pinned input changed: {path}")


def stable_uint64(*parts: object) -> int:
    text = "\x1f".join(map(str, parts)).encode()
    return int.from_bytes(hashlib.sha256(text).digest()[:8], "little")


def split_one_target(
    umi: np.ndarray,
    batch_labels: Sequence[str],
    donor: str,
    target: str,
    seed: int,
) -> np.ndarray:
    """Assign occupied batches to two UMI-balanced halves.

    The random-looking order comes only from SHA-256. Each entire batch remains
    in one half. Inactive batches are -1.
    """
    values = np.asarray(umi, dtype=np.float64)
    require(values.ndim == 1, "UMI vector must be one-dimensional")
    require(len(values) == len(batch_labels), "batch-label length mismatch")
    active = np.flatnonzero(values > 0)
    require(len(active) >= 2, f"{donor}/{target}: fewer than two occupied batches")
    order = sorted(
        map(int, active),
        key=lambda batch: stable_uint64(
            "split-order", seed, donor, target, batch_labels[batch]
        ),
    )
    assignment = np.full(len(values), -1, dtype=np.int8)
    totals = [0.0, 0.0]
    for position, batch in enumerate(order):
        if position < 2:
            half = position
        elif totals[0] < totals[1]:
            half = 0
        elif totals[1] < totals[0]:
            half = 1
        else:
            half = (
                stable_uint64("split-tie", seed, donor, target, batch_labels[batch]) & 1
            )
        assignment[batch] = half
        totals[half] += float(values[batch])
    for _ in range(32):
        difference = totals[0] - totals[1]
        current = abs(difference)
        best: tuple[float, str, int, int] | None = None
        heavy = 0 if difference > 0 else 1
        light = 1 - heavy
        heavy_batches = np.flatnonzero(assignment == heavy)
        light_batches = np.flatnonzero(assignment == light)
        if len(heavy_batches) > 1:
            for batch in map(int, heavy_batches):
                candidate = abs(current - 2.0 * values[batch])
                choice = (candidate, "move", batch, -1)
                if candidate + EPS < current and (best is None or choice < best):
                    best = choice
        target_difference = current / 2.0
        light_order = light_batches[np.argsort(values[light_batches], kind="stable")]
        light_weights = values[light_order]
        for heavy_batch in map(int, heavy_batches):
            wanted = values[heavy_batch] - target_difference
            position = int(np.searchsorted(light_weights, wanted))
            for local in (position - 1, position):
                if 0 <= local < len(light_order):
                    light_batch = int(light_order[local])
                    candidate = abs(
                        current - 2.0 * (values[heavy_batch] - values[light_batch])
                    )
                    choice = (candidate, "swap", heavy_batch, light_batch)
                    if candidate + EPS < current and (best is None or choice < best):
                        best = choice
        if best is None:
            break
        _, action, first, second = best
        if action == "move":
            assignment[first] = light
            totals[heavy] -= float(values[first])
            totals[light] += float(values[first])
        else:
            assignment[first] = light
            assignment[second] = heavy
            totals[heavy] += float(values[second] - values[first])
            totals[light] += float(values[first] - values[second])
    require(np.all(assignment[active] >= 0), "occupied batch was not assigned")
    require(np.all(assignment[values <= 0] == -1), "inactive batch was assigned")
    return assignment


def canonical_split_bytes(assignment: np.ndarray) -> bytes:
    """Canonicalize a split modulo swapping the A/B half labels."""
    values = np.asarray(assignment, dtype=np.int8).copy()
    active = np.flatnonzero(values >= 0)
    require(len(active) >= 2, "cannot canonicalize an empty split")
    if values[active[0]] == 1:
        values[active] = 1 - values[active]
    return values.tobytes()


def split_support(
    assignment: np.ndarray,
    target_umi: np.ndarray,
    target_cells: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Summarize one donor-wide batch partition on every target."""
    half_umi = np.stack(
        [target_umi[:, assignment == half].sum(axis=1) for half in range(HALVES)]
    )
    half_cells = np.stack(
        [target_cells[:, assignment == half].sum(axis=1) for half in range(HALVES)]
    ).astype(np.int32)
    half_batches = np.stack(
        [
            (target_cells[:, assignment == half] > 0).sum(axis=1)
            for half in range(HALVES)
        ]
    ).astype(np.int16)
    valid = np.all(half_cells >= MIN_CELLS_HALF, axis=0) & np.all(
        half_batches >= MIN_BATCHES_HALF, axis=0
    )
    return half_umi, half_cells, half_batches, valid


def global_split_score(
    assignment: np.ndarray,
    target_umi: np.ndarray,
    target_cells: np.ndarray,
) -> tuple[int, float, float, float, float]:
    """Lexicographic metadata-only objective for a coherent gallery split."""
    half_umi, _, _, valid = split_support(assignment, target_umi, target_cells)
    smaller = np.minimum(half_umi[0], half_umi[1])
    larger = np.maximum(half_umi[0], half_umi[1])
    balance = np.divide(
        smaller,
        larger,
        out=np.zeros_like(smaller),
        where=larger > 0,
    )
    return (
        int(valid.sum()),
        float(balance.min()),
        float(np.quantile(balance, 0.01)),
        float(np.quantile(balance, 0.05)),
        float(np.median(balance)),
    )


def make_split_plan(
    donor: str,
    targets: Sequence[str],
    batch_labels: Sequence[str],
    target_umi: np.ndarray,
    target_cells: np.ndarray,
) -> dict[str, np.ndarray]:
    require(target_umi.shape == target_cells.shape, f"{donor}: metadata shape mismatch")
    require(
        target_umi.shape == (len(targets), len(batch_labels)), f"{donor}: axis mismatch"
    )
    assignments = np.empty((SPLITS, len(batch_labels)), dtype=np.int8)
    half_umi = np.zeros((SPLITS, HALVES, len(targets)), dtype=np.float64)
    half_cells = np.zeros((SPLITS, HALVES, len(targets)), dtype=np.int32)
    half_batches = np.zeros((SPLITS, HALVES, len(targets)), dtype=np.int16)
    effective_seeds = np.empty(SPLITS, dtype=np.uint64)
    global_umi = target_umi.sum(axis=0)
    seen: set[bytes] = set()
    for split, seed in enumerate(SPLIT_SEEDS):
        best_assignment: np.ndarray | None = None
        best_score: tuple[int, float, float, float, float] | None = None
        best_seed = 0
        for attempt in range(GLOBAL_SPLIT_SEARCH):
            effective_seed = (
                seed
                if attempt == 0
                else stable_uint64("global-split-search", seed, attempt)
            )
            candidate = split_one_target(
                global_umi,
                batch_labels,
                donor,
                "__donor_wide__",
                effective_seed,
            )
            key = canonical_split_bytes(candidate)
            if key in seen:
                continue
            score = global_split_score(candidate, target_umi, target_cells)
            if best_score is None or score > best_score:
                best_assignment = candidate
                best_score = score
                best_seed = effective_seed
        require(
            best_assignment is not None,
            f"{donor}: unable to construct eight distinct donor-wide splits",
        )
        seen.add(canonical_split_bytes(best_assignment))
        assignments[split] = best_assignment
        effective_seeds[split] = best_seed
        local_umi, local_cells, local_batches, _ = split_support(
            best_assignment, target_umi, target_cells
        )
        half_umi[split] = local_umi
        half_cells[split] = local_cells
        half_batches[split] = local_batches
    valid = np.all(half_cells >= MIN_CELLS_HALF, axis=1) & np.all(
        half_batches >= MIN_BATCHES_HALF, axis=1
    )
    return {
        "assignments": assignments,
        "half_umi": half_umi,
        "half_cells": half_cells,
        "half_batches": half_batches,
        "valid": valid,
        "valid_pairs": valid.sum(axis=0).astype(np.int8),
        "occupied_batches": np.sum(target_cells > 0, axis=1).astype(np.int16),
        "batch_labels": np.asarray(batch_labels, dtype="U64"),
        "effective_seeds": effective_seeds,
    }


def donor_metadata_paths(donor: str) -> tuple[Path, Path, Path]:
    base = EXP140_PREPARE / donor
    return (
        base / "development_batch_umi.npy",
        base / "development_batch_cells.npy",
        base / "receipt.json",
    )


def summarize_split_plan(plan: Mapping[str, np.ndarray]) -> dict[str, Any]:
    valid_pairs = np.asarray(plan["valid_pairs"])
    all_eight = float(np.mean(valid_pairs == SPLITS))
    at_least_six = float(np.mean(valid_pairs >= MIN_VALID_PAIRS))
    half_cells = np.asarray(plan["half_cells"])
    half_umi = np.asarray(plan["half_umi"])
    smaller = np.minimum(half_umi[:, 0], half_umi[:, 1])
    larger = np.maximum(half_umi[:, 0], half_umi[:, 1])
    balance = np.divide(
        smaller,
        larger,
        out=np.zeros_like(smaller),
        where=larger > 0,
    )
    half_batches = np.asarray(plan["half_batches"])
    assignments = np.asarray(plan["assignments"])
    unique_splits = len(
        {
            canonical_split_bytes(assignments[split])
            for split in range(assignments.shape[0])
        }
    )
    return {
        "valid_pairs_min": int(valid_pairs.min()),
        "valid_pairs_median": float(np.median(valid_pairs)),
        "valid_pairs_max": int(valid_pairs.max()),
        "fraction_all_eight": all_eight,
        "fraction_at_least_six": at_least_six,
        "minimum_cells_per_half": int(half_cells.min()),
        "median_cells_per_half": float(np.median(half_cells)),
        "minimum_batches_per_half": int(half_batches.min()),
        "median_batches_per_half": float(np.median(half_batches)),
        "minimum_umi_balance": float(balance.min()),
        "median_umi_balance": float(np.median(balance)),
        "unique_donor_wide_splits": unique_splits,
        "support_gate_metadata_only": bool(
            all_eight >= ALL_EIGHT_FRACTION
            and at_least_six == 1.0
            and unique_splits == SPLITS
        ),
    }


def donor_assignment_nulls(donors: Sequence[str], draws: int, seed: int) -> np.ndarray:
    """Permute donor identities independently per target for each null draw."""
    rng = np.random.default_rng(seed)
    base = np.arange(len(donors), dtype=np.int8)
    result = np.empty((draws, DEV_COUNT, len(donors)), dtype=np.int8)
    for draw in range(draws):
        for target in range(DEV_COUNT):
            result[draw, target] = rng.permutation(base)
        require(
            np.any(result[draw] != base[None, :]),
            "donor-assignment null draw is identity",
        )
    return result


def create_plan() -> None:
    validate_pinned_inputs()
    require(not PLAN_DIR.exists(), f"plan already exists: {PLAN_DIR}")
    require(not SOURCE_HASH_MANIFEST.exists(), "source freeze already exists")
    targets = read_lines(DEV_TARGETS)
    responses = read_lines(RESPONSE_GENES)
    require(len(targets) == DEV_COUNT, "development target count changed")
    require(len(responses) == RESPONSE_COUNT, "response gene count changed")
    parent_arrays = np.load(EXP140_PLAN / "plan-arrays.npz", allow_pickle=False)
    require(
        parent_arrays["development_folds"].shape == (DEV_COUNT,), "parent folds changed"
    )
    require(
        parent_arrays["development_blocks"].shape == (DEV_COUNT,),
        "parent blocks changed",
    )

    arrays: dict[str, np.ndarray] = {
        "development_folds": np.asarray(
            parent_arrays["development_folds"], dtype=np.int8
        ),
        "development_blocks": np.asarray(
            parent_arrays["development_blocks"], dtype=np.int16
        ),
    }
    for outer, donors in OUTER_DONORS.items():
        arrays[f"donor_assignment_null_{outer}"] = donor_assignment_nulls(
            donors,
            NULL_DRAWS,
            stable_uint64(DONOR_NULL_SEED, outer),
        )
    donor_reports = {}
    source_schema = {
        "schema_version": 1,
        "created_before_source_expression_scan": True,
        "development_outcomes_opened": False,
        "sealed_410_effect_values_opened": False,
        "parent_hashes": {
            str(path.relative_to(ROOT)): digest for path, digest in EXPECTED.items()
        },
        "donors": {},
    }
    for donor in DONORS:
        umi_path, cells_path, receipt_path = donor_metadata_paths(donor)
        target_umi = np.asarray(np.load(umi_path, mmap_mode="r"), dtype=np.float64)
        target_cells = np.asarray(np.load(cells_path, mmap_mode="r"), dtype=np.int32)
        receipt = read_json(receipt_path)
        batch_labels = tuple(map(str, receipt["batch_labels"]))
        split_plan = make_split_plan(
            donor, targets, batch_labels, target_umi, target_cells
        )
        for key, value in split_plan.items():
            arrays[f"{donor}_{key}"] = value
        report = summarize_split_plan(split_plan)
        donor_reports[donor] = report
        source_schema["donors"][donor] = {
            "batches": len(batch_labels),
            "batch_labels_sha256": hashlib.sha256(
                "\n".join(batch_labels).encode()
            ).hexdigest(),
            "metadata": {
                str(umi_path.relative_to(ROOT)): file_record(umi_path),
                str(cells_path.relative_to(ROOT)): file_record(cells_path),
                str(receipt_path.relative_to(ROOT)): file_record(receipt_path),
            },
            "target_batch_umi_shape": list(target_umi.shape),
            "target_batch_cells_shape": list(target_cells.shape),
            "occupied_batches": {
                "min": int(np.sum(target_cells > 0, axis=1).min()),
                "median": float(np.median(np.sum(target_cells > 0, axis=1))),
                "max": int(np.sum(target_cells > 0, axis=1).max()),
            },
            "split_report": report,
        }
        require(
            report["support_gate_metadata_only"], f"{donor}: split support gate failed"
        )

    PLAN_DIR.mkdir(parents=True)
    atomic_npz(PLAN_DIR / "plan-arrays.npz", **arrays)
    atomic_json(PLAN_DIR / "source-schema.json", source_schema)
    plan_manifest = {
        "schema_version": 1,
        "created_utc": provenance_utc(),
        "created_before_source_expression_scan": True,
        "development_outcomes_opened": False,
        "sealed_410_effect_values_opened": False,
        "development_targets": len(targets),
        "response_genes": len(responses),
        "split_pairs": SPLITS,
        "split_seeds": list(SPLIT_SEEDS),
        "global_split_search": GLOBAL_SPLIT_SEARCH,
        "minimum_cells_per_half": MIN_CELLS_HALF,
        "minimum_batches_per_half": MIN_BATCHES_HALF,
        "minimum_valid_pairs": MIN_VALID_PAIRS,
        "etas": list(ETAS),
        "bootstrap_draws": BOOTSTRAP_DRAWS,
        "null_draws": NULL_DRAWS,
        "donor_assignment_null_seed": DONOR_NULL_SEED,
        "donors": donor_reports,
        "files": {
            "plan-arrays.npz": file_record(PLAN_DIR / "plan-arrays.npz"),
            "source-schema.json": file_record(PLAN_DIR / "source-schema.json"),
        },
        "peak_rss_gib": peak_rss_gib(),
    }
    atomic_json(PLAN_DIR / "plan-manifest.json", plan_manifest)
    print(json.dumps(plan_manifest, indent=2, sort_keys=True))


def validate_plan() -> dict[str, Any]:
    validate_pinned_inputs()
    manifest_path = PLAN_DIR / "plan-manifest.json"
    require(manifest_path.is_file(), "Experiment 142 plan missing")
    manifest = read_json(manifest_path)
    require(
        manifest["created_before_source_expression_scan"] is True, "plan timing changed"
    )
    require(manifest["development_outcomes_opened"] is False, "plan opened outcomes")
    require(manifest["sealed_410_effect_values_opened"] is False, "sealed split opened")
    for relative, expected in manifest["files"].items():
        path = PLAN_DIR / relative
        require(
            path.stat().st_size == expected["bytes"], f"plan bytes changed: {relative}"
        )
        require(
            sha256_file(path) == expected["sha256"], f"plan hash changed: {relative}"
        )
    arrays = np.load(PLAN_DIR / "plan-arrays.npz", allow_pickle=False)
    for donor in DONORS:
        assignments = np.asarray(arrays[f"{donor}_assignments"], dtype=np.int8)
        require(
            assignments.shape == (SPLITS, len(arrays[f"{donor}_batch_labels"])),
            f"{donor}: donor-wide split shape changed",
        )
        require(
            len({canonical_split_bytes(assignments[split]) for split in range(SPLITS)})
            == SPLITS,
            f"{donor}: donor-wide splits are not unique",
        )
        valid_pairs = arrays[f"{donor}_valid_pairs"]
        require(np.all(valid_pairs >= MIN_VALID_PAIRS), f"{donor}: support below six")
        require(
            np.mean(valid_pairs == SPLITS) >= ALL_EIGHT_FRACTION,
            f"{donor}: all-eight support below 95%",
        )
        require(
            np.all(arrays[f"{donor}_half_batches"] >= MIN_BATCHES_HALF),
            f"{donor}: split half has too few occupied batches",
        )
    for outer, donors in OUTER_DONORS.items():
        nulls = np.asarray(arrays[f"donor_assignment_null_{outer}"], dtype=np.int8)
        require(
            nulls.shape == (NULL_DRAWS, DEV_COUNT, len(donors)),
            f"{outer}: donor-assignment null shape changed",
        )
        expected = np.arange(len(donors), dtype=np.int8)
        require(
            np.all(np.sort(nulls, axis=2) == expected[None, None, :]),
            f"{outer}: donor-assignment null is not a permutation",
        )
    return manifest


def source_input_relatives() -> dict[str, tuple[str, ...]]:
    parent = read_json(EXP140_SOURCE_MANIFEST)["files"]
    result: dict[str, tuple[str, ...]] = {}
    k562 = (
        str(K562_RAW.relative_to(ROOT)),
        str((DONOR_ROOTS["k562_gwps"] / "genes.txt").relative_to(ROOT)),
        str((DONOR_ROOTS["k562_gwps"] / "perts.txt").relative_to(ROOT)),
        str((DONOR_ROOTS["k562_gwps"] / "delta.npy").relative_to(ROOT)),
        str((DONOR_ROOTS["k562_gwps"] / "b_ctrl.npy").relative_to(ROOT)),
    )
    require(
        all(item in parent for item in k562), "K562 source absent from parent manifest"
    )
    result["k562_gwps"] = k562
    gene_meta = str(ORION_GENE_METADATA.relative_to(ROOT))
    for donor, line in DONOR_LINES.items():
        prefix = f"data/external/orion/data/{line}_Batch"
        files = tuple(
            sorted(
                (
                    key
                    for key in parent
                    if key.startswith(prefix) and key.endswith(".parquet")
                ),
                key=lambda value: int(value.rsplit("Batch", 1)[1].split(".")[0]),
            )
        )
        extras = (
            gene_meta,
            str((DONOR_ROOTS[donor] / "genes.txt").relative_to(ROOT)),
            str((DONOR_ROOTS[donor] / "perts.txt").relative_to(ROOT)),
            str((DONOR_ROOTS[donor] / "delta.npy").relative_to(ROOT)),
            str((DONOR_ROOTS[donor] / "b_ctrl.npy").relative_to(ROOT)),
        )
        require(
            files and all(item in parent for item in extras),
            f"{donor}: source manifest gap",
        )
        result[donor] = files + extras
    return result


def freeze() -> None:
    plan = validate_plan()
    require(
        SOURCE_ONLY_FREEZE_AUTHORIZATION.is_file(),
        "source-only freeze authorization missing",
    )
    source_authorization = read_json(SOURCE_ONLY_FREEZE_AUTHORIZATION)
    require(
        source_authorization["source_only_freeze_authorized"] is True,
        "source-only freeze disabled",
    )
    require(
        source_authorization["source_builds_authorized"] is True,
        "source builds disabled",
    )
    require(
        source_authorization["source_evaluations_authorized"] is True,
        "source evaluations disabled",
    )
    require(
        source_authorization["source_gate_authorized"] is True,
        "source gate disabled",
    )
    require(
        source_authorization["development_613_authorized"] is False,
        "development unexpectedly enabled",
    )
    require(
        source_authorization["untouched_410_test_authorized"] is False,
        "sealed test unexpectedly enabled",
    )
    require(
        source_authorization["reserved_confirmation_expression_authorized"] is False,
        "reserved confirmation unexpectedly enabled",
    )
    require(
        source_authorization["count_generation_authorized"] is False,
        "count generation unexpectedly enabled",
    )
    require(
        source_authorization["cell_eval2_authorized"] is False,
        "cell-eval2 unexpectedly enabled",
    )
    require(
        source_authorization["challenge_build_authorized"] is False,
        "challenge build unexpectedly enabled",
    )
    require(
        source_authorization["submission_authorized"] is False,
        "submission unexpectedly enabled",
    )
    freeze_created_utc = authorized_freeze_utc(source_authorization)
    require(
        source_authorization["bindings"]["plan_manifest_sha256"]
        == sha256_file(PLAN_DIR / "plan-manifest.json"),
        "source authorization plan binding changed",
    )
    require(
        source_authorization["bindings"]["prefreeze_amendment_sha256"]
        == sha256_file(PREFREEZE_AMENDMENT),
        "source authorization amendment binding changed",
    )
    require(
        source_authorization["bindings"]["prefreeze_amendment_2_sha256"]
        == sha256_file(PREFREEZE_AMENDMENT_2),
        "source authorization second amendment binding changed",
    )
    require(
        source_authorization["bindings"]["prefreeze_amendment_3_sha256"]
        == sha256_file(PREFREEZE_AMENDMENT_3),
        "source authorization third amendment binding changed",
    )
    require(
        source_authorization["bindings"]["prefreeze_amendment_4_sha256"]
        == sha256_file(PREFREEZE_AMENDMENT_4),
        "source authorization fourth amendment binding changed",
    )
    require(
        source_authorization["bindings"]["prefreeze_amendment_5_sha256"]
        == sha256_file(PREFREEZE_AMENDMENT_5),
        "source authorization fifth amendment binding changed",
    )
    require(
        source_authorization["bindings"]["prefreeze_amendment_6_sha256"]
        == sha256_file(PREFREEZE_AMENDMENT_6),
        "source authorization sixth amendment binding changed",
    )
    require(
        source_authorization["bindings"]["provenance_clock_correction_sha256"]
        == sha256_file(PROVENANCE_CLOCK_CORRECTION),
        "source authorization clock-correction binding changed",
    )
    require(
        source_authorization["bindings"]["superseded_clock_correction_sha256"]
        == sha256_file(SUPERSEDED_CLOCK_CORRECTION),
        "superseded clock-correction artifact changed",
    )
    require(
        source_authorization["bindings"]["superseded_v5_authorization_sha256"]
        == sha256_file(SUPERSEDED_SOURCE_AUTHORIZATION_V5),
        "superseded v5 authorization changed",
    )
    require(
        source_authorization["bindings"]["superseded_v4_authorization_sha256"]
        == sha256_file(SUPERSEDED_SOURCE_AUTHORIZATION_V4),
        "superseded v4 authorization changed",
    )
    require(
        source_authorization["bindings"]["superseded_v3_authorization_sha256"]
        == sha256_file(SUPERSEDED_SOURCE_AUTHORIZATION_V3),
        "superseded v3 authorization changed",
    )
    require(
        source_authorization["bindings"]["superseded_v2_authorization_sha256"]
        == sha256_file(SUPERSEDED_SOURCE_AUTHORIZATION_V2),
        "superseded v2 authorization changed",
    )
    require(
        source_authorization["bindings"][
            "superseded_v1_preoutcome_authorization_sha256"
        ]
        == sha256_file(SUPERSEDED_PREOUTCOME_AUTHORIZATION_V1),
        "superseded v1 preoutcome authorization changed",
    )
    require(
        source_authorization["bindings"]["superseded_v1_k562_build_manifest_sha256"]
        == sha256_file(SUPERSEDED_K562_BUILD_MANIFEST_V1),
        "superseded v1 K562 build manifest changed",
    )
    require(
        source_authorization["bindings"][
            "superseded_v2_preoutcome_authorization_sha256"
        ]
        == sha256_file(SUPERSEDED_PREOUTCOME_AUTHORIZATION_V2),
        "superseded v2 preoutcome authorization changed",
    )
    require(
        source_authorization["bindings"]["superseded_v2_k562_build_manifest_sha256"]
        == sha256_file(SUPERSEDED_K562_BUILD_MANIFEST_V2),
        "superseded v2 K562 build manifest changed",
    )
    require(
        source_authorization["bindings"][
            "superseded_v2_k562_evaluation_manifest_sha256"
        ]
        == sha256_file(SUPERSEDED_K562_EVALUATION_MANIFEST_V2),
        "superseded v2 K562 evaluation manifest changed",
    )
    require(
        source_authorization["bindings"][
            "superseded_v3_preoutcome_authorization_sha256"
        ]
        == sha256_file(SUPERSEDED_PREOUTCOME_AUTHORIZATION_V3),
        "superseded v3 preoutcome authorization changed",
    )
    require(
        source_authorization["bindings"]["superseded_v3_k562_build_manifest_sha256"]
        == sha256_file(SUPERSEDED_K562_BUILD_MANIFEST_V3),
        "superseded v3 K562 build manifest changed",
    )
    require(
        source_authorization["bindings"][
            "superseded_v3_k562_evaluation_manifest_sha256"
        ]
        == sha256_file(SUPERSEDED_K562_EVALUATION_MANIFEST_V3),
        "superseded v3 K562 evaluation manifest changed",
    )
    for relative, expected in source_authorization["implementation"].items():
        path = ROOT / relative
        require(path.is_file(), f"authorized implementation missing: {relative}")
        require(
            sha256_file(path) == expected,
            f"authorized implementation changed: {relative}",
        )
    require(not FREEZE_DIR.exists(), "preoutcome freeze already exists")
    require(not SOURCE_BUILD_DIR.exists(), "source build already exists")
    validate_directory_manifest(EXP140_PREPARE)
    validate_directory_manifest(EXP140_IDENTITY)
    parent_manifest = read_json(EXP140_SOURCE_MANIFEST)
    raw_by_donor = source_input_relatives()
    records = {}
    for relatives in raw_by_donor.values():
        for relative in relatives:
            records[relative] = parent_manifest["files"][relative]
    supporting = (
        EXP140_PLAN / "plan-arrays.npz",
        EXP140_PLAN / "plan-manifest.json",
        EXP140_PREPARE / "artifact-manifest.json",
        EXP140_IDENTITY / "artifact-manifest.json",
        EXP140_IDENTITY / "whitener-k562.npz",
        EXP140_IDENTITY / "whitener-non_k562.npz",
        EXP132_SOURCE_AUTHORIZATION,
        EXP132_CANDIDATE_FREEZE,
        EXP132_PLAN_MANIFEST,
        EXP132_SELECTION_MANIFEST,
        DEV_TARGETS,
        RESPONSE_GENES,
    )
    for path in supporting:
        records[str(path.relative_to(ROOT))] = file_record(path)
    implementation = {
        str(path.relative_to(ROOT)): sha256_file(path)
        for path in (
            README,
            PROTOCOL,
            TESTS,
            HERE / "run.py",
            DEVELOPMENT_INTERFACE,
            PREFREEZE_AMENDMENT,
            PREFREEZE_AMENDMENT_2,
            PREFREEZE_AMENDMENT_3,
            PREFREEZE_AMENDMENT_4,
            PREFREEZE_AMENDMENT_5,
            PREFREEZE_AMENDMENT_6,
            PROVENANCE_CLOCK_CORRECTION,
            SUPERSEDED_CLOCK_CORRECTION,
            SOURCE_ONLY_FREEZE_AUTHORIZATION,
        )
    }
    FREEZE_DIR.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(
        tempfile.mkdtemp(prefix=f".{FREEZE_DIR.name}.", dir=FREEZE_DIR.parent)
    )
    try:
        source_manifest_path = temporary / SOURCE_HASH_MANIFEST.name
        preoutcome_path = temporary / PREOUTCOME_AUTHORIZATION.name
        atomic_json(
            source_manifest_path,
            {
                "schema_version": 1,
                "created_utc": freeze_created_utc,
                "created_before_source_expression_scan": True,
                "origin_manifest": {
                    "path": str(EXP140_SOURCE_MANIFEST.relative_to(ROOT)),
                    "sha256": EXPECTED[EXP140_SOURCE_MANIFEST],
                },
                "donor_files": {
                    key: list(value) for key, value in raw_by_donor.items()
                },
                "files": records,
            },
        )
        atomic_json(
            preoutcome_path,
            {
                "schema_version": 1,
                "created_utc": freeze_created_utc,
                "source_only_work_authorized": True,
                "development_613_authorized": False,
                "untouched_410_test_authorized": False,
                "reserved_confirmation_expression_authorized": False,
                "count_generation_authorized": False,
                "cell_eval2_authorized": False,
                "challenge_build_authorized": False,
                "submission_authorized": False,
                "source_only_freeze_authorization_sha256": sha256_file(
                    SOURCE_ONLY_FREEZE_AUTHORIZATION
                ),
                "protocol_sha256": sha256_file(PROTOCOL),
                "plan_manifest_sha256": sha256_file(PLAN_DIR / "plan-manifest.json"),
                "source_hash_manifest_sha256": sha256_file(source_manifest_path),
                "implementation": implementation,
                "plan": plan,
            },
        )
        os.replace(temporary, FREEZE_DIR)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    print(json.dumps(read_json(PREOUTCOME_AUTHORIZATION), indent=2, sort_keys=True))


def validate_preoutcome() -> dict[str, Any]:
    validate_plan()
    require(PREOUTCOME_AUTHORIZATION.is_file(), "preoutcome authorization missing")
    require(SOURCE_HASH_MANIFEST.is_file(), "source hash manifest missing")
    auth = read_json(PREOUTCOME_AUTHORIZATION)
    require(auth["source_only_work_authorized"] is True, "source-only work disabled")
    source_authorization = read_json(SOURCE_ONLY_FREEZE_AUTHORIZATION)
    require(
        source_authorization["source_builds_authorized"] is True,
        "source builds disabled",
    )
    require(
        source_authorization["source_evaluations_authorized"] is True,
        "source evaluations disabled",
    )
    require(
        source_authorization["source_gate_authorized"] is True,
        "source gate disabled",
    )
    require(
        source_authorization["development_613_authorized"] is False,
        "development unexpectedly enabled",
    )
    require(
        source_authorization["untouched_410_test_authorized"] is False,
        "sealed test unexpectedly enabled",
    )
    require(
        source_authorization["reserved_confirmation_expression_authorized"] is False,
        "reserved confirmation unexpectedly enabled",
    )
    require(
        source_authorization["count_generation_authorized"] is False,
        "count generation unexpectedly enabled",
    )
    require(
        source_authorization["cell_eval2_authorized"] is False,
        "cell-eval2 unexpectedly enabled",
    )
    require(
        source_authorization["challenge_build_authorized"] is False,
        "challenge build unexpectedly enabled",
    )
    require(
        source_authorization["submission_authorized"] is False,
        "submission unexpectedly enabled",
    )
    require(
        auth["development_613_authorized"] is False, "development unexpectedly enabled"
    )
    require(
        auth["untouched_410_test_authorized"] is False,
        "sealed test unexpectedly enabled",
    )
    require(
        auth["reserved_confirmation_expression_authorized"] is False,
        "reserved confirmation unexpectedly enabled",
    )
    require(auth["count_generation_authorized"] is False, "count generation enabled")
    require(auth["cell_eval2_authorized"] is False, "cell-eval2 enabled")
    require(auth["challenge_build_authorized"] is False, "challenge build enabled")
    require(auth["submission_authorized"] is False, "submission enabled")
    require(
        auth["source_only_freeze_authorization_sha256"]
        == sha256_file(SOURCE_ONLY_FREEZE_AUTHORIZATION),
        "source-only freeze authorization changed",
    )
    require(auth["protocol_sha256"] == sha256_file(PROTOCOL), "protocol changed")
    require(
        auth["plan_manifest_sha256"] == sha256_file(PLAN_DIR / "plan-manifest.json"),
        "plan changed",
    )
    require(
        auth["source_hash_manifest_sha256"] == sha256_file(SOURCE_HASH_MANIFEST),
        "source manifest changed",
    )
    for relative, expected in auth["implementation"].items():
        path = ROOT / relative
        require(path.is_file(), f"implementation missing: {relative}")
        require(sha256_file(path) == expected, f"implementation changed: {relative}")
    return auth


def validate_donor_source_hashes(donor: str) -> None:
    validate_preoutcome()
    manifest = read_json(SOURCE_HASH_MANIFEST)
    for relative in manifest["donor_files"][donor]:
        path = ROOT / relative
        expected = manifest["files"][relative]
        require(path.is_file(), f"source missing: {relative}")
        require(
            path.stat().st_size == expected["bytes"],
            f"source bytes changed: {relative}",
        )
        print(f"[hash] {relative}", flush=True)
        require(
            sha256_file(path) == expected["sha256"], f"source hash changed: {relative}"
        )


def validate_source_manifest_paths(paths: Sequence[Path]) -> None:
    manifest = read_json(SOURCE_HASH_MANIFEST)["files"]
    for path in paths:
        relative = str(path.relative_to(ROOT))
        require(relative in manifest, f"path absent from source manifest: {relative}")
        expected = manifest[relative]
        require(path.is_file(), f"source missing: {relative}")
        require(
            path.stat().st_size == expected["bytes"],
            f"source bytes changed: {relative}",
        )
        require(
            sha256_file(path) == expected["sha256"], f"source hash changed: {relative}"
        )


@contextmanager
def donor_lock(donor: str) -> Iterator[None]:
    SOURCE_BUILD_DIR.mkdir(parents=True, exist_ok=True)
    lock = SOURCE_BUILD_DIR / ".active-donor.lock"
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as error:
        raise ContractError("another Experiment 142 donor scan is active") from error
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(f"{donor}\n{os.getpid()}\n")
        yield
    finally:
        lock.unlink(missing_ok=True)


def h5_categorical(group: h5py.Group) -> tuple[np.ndarray, np.ndarray]:
    categories = np.asarray(group["categories"]).astype(str)
    codes = np.asarray(group["codes"], dtype=np.int64)
    require(np.all(codes >= 0), "categorical contains missing codes")
    return categories, codes


def numbered_orion_files(line: str) -> tuple[Path, ...]:
    paths = tuple(
        sorted(
            ORION_DATA.glob(f"{line}_Batch*.parquet"),
            key=lambda path: int(path.stem.rsplit("Batch", 1)[1]),
        )
    )
    require(bool(paths), f"no Orion files for {line}")
    return paths


def orion_gene_axis() -> tuple[str, ...]:
    table = pq.read_table(ORION_GENE_METADATA, columns=["gene_token_id", "gene_name"])
    tokens = np.asarray(table["gene_token_id"].combine_chunks(), dtype=np.int64)
    names = np.asarray(table["gene_name"].combine_chunks(), dtype=str)
    order = np.argsort(tokens, kind="stable")
    require(
        np.array_equal(tokens[order], np.arange(tokens.max() + 1)), "Orion token gap"
    )
    return tuple(map(str, names[order]))


def release_arrow_memory() -> None:
    gc.collect()
    pa.default_memory_pool().release_unused()


def validate_batch_axis(
    observed: Sequence[str], planned: Sequence[str], donor: str
) -> None:
    require(
        tuple(map(str, observed)) == tuple(map(str, planned)),
        f"{donor}: ordered batch labels changed",
    )


def response_column_plan(
    source_genes: Sequence[str], responses: Sequence[str]
) -> tuple[np.ndarray, np.ndarray]:
    index = {gene: position for position, gene in enumerate(source_genes)}
    require(
        all(gene in index for gene in responses), "response gene absent from source"
    )
    columns = np.asarray([index[gene] for gene in responses], dtype=np.int64)
    sorted_columns = np.sort(columns)
    inverse = np.searchsorted(sorted_columns, columns)
    require(np.array_equal(sorted_columns[inverse], columns), "column reorder failed")
    return sorted_columns, inverse


def accumulate_grouped_rows(
    destination: np.ndarray,
    values: np.ndarray,
    targets: np.ndarray,
    halves: np.ndarray,
) -> None:
    groups = targets.astype(np.int64) * HALVES + halves.astype(np.int64)
    indicator = sp.csr_matrix(
        (
            np.ones(len(groups), dtype=np.float64),
            (groups, np.arange(len(groups), dtype=np.int64)),
        ),
        shape=(DEV_COUNT * HALVES, len(groups)),
    )
    aggregated = np.asarray(indicator @ values, dtype=np.float64).reshape(
        DEV_COUNT, HALVES, RESPONSE_COUNT
    )
    destination[0] += aggregated[:, 0]
    destination[1] += aggregated[:, 1]


def scan_k562_source(
    assignments: np.ndarray,
    planned_batch_labels: Sequence[str],
    responses: Sequence[str],
    limit: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    targets = read_lines(DEV_TARGETS)
    target_index = {target: index for index, target in enumerate(targets)}
    source_genes = read_lines(DONOR_ROOTS["k562_gwps"] / "genes.txt")
    sorted_columns, inverse = response_column_plan(source_genes, responses)
    sums = np.zeros((SPLITS, HALVES, DEV_COUNT, RESPONSE_COUNT), dtype=np.float64)
    cells = np.zeros((SPLITS, HALVES, DEV_COUNT), dtype=np.int32)
    umi = np.zeros((SPLITS, HALVES, DEV_COUNT), dtype=np.float64)
    with h5py.File(K562_RAW, "r") as handle:
        categories, codes = h5_categorical(handle["obs"]["perturbation"])
        batches = np.asarray(handle["obs"]["batch"], dtype=np.int64)
        ncounts = np.asarray(handle["obs"]["ncounts"], dtype=np.float64)
        batch_values, batch_codes = np.unique(batches, return_inverse=True)
        validate_batch_axis(batch_values, planned_batch_labels, "k562_gwps")
        category_target = np.asarray(
            [target_index.get(value, -1) for value in categories], dtype=np.int32
        )
        row_target = category_target[codes]
        x = handle["X"]
        block = 2048
        for lo in range(0, x.shape[0], block):
            hi = min(lo + block, x.shape[0])
            local_targets = row_target[lo:hi]
            use = local_targets >= 0
            if np.any(use):
                values = np.asarray(x[lo:hi, sorted_columns], dtype=np.float64)
                values = values[use][:, inverse]
                selected_targets = local_targets[use]
                selected_batches = batch_codes[lo:hi][use]
                selected_umi = ncounts[lo:hi][use]
                for split in range(SPLITS):
                    selected_halves = assignments[split, selected_batches]
                    require(
                        np.all(selected_halves >= 0), "K562 active batch unassigned"
                    )
                    accumulate_grouped_rows(
                        sums[split], values, selected_targets, selected_halves
                    )
                    flat = selected_targets.astype(np.int64) * HALVES + selected_halves
                    cells[split] += (
                        np.bincount(flat, minlength=DEV_COUNT * HALVES)
                        .reshape(DEV_COUNT, HALVES)
                        .T.astype(np.int32)
                    )
                    umi[split] += (
                        np.bincount(
                            flat,
                            weights=selected_umi,
                            minlength=DEV_COUNT * HALVES,
                        )
                        .reshape(DEV_COUNT, HALVES)
                        .T
                    )
                del values
            if (lo // block) % 100 == 99 or hi == x.shape[0]:
                print(f"[k562_gwps] {hi:,}/{x.shape[0]:,}", flush=True)
                memory_check("K562 source scan", limit)
        require(len(batch_values) == assignments.shape[1], "K562 batch count changed")
    return sums, cells, umi


def sparse_orion_batch_sum(
    table: pa.Table,
    target_index: Mapping[str, int],
    token_to_response: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    labels = np.asarray(table["gene_target"].combine_chunks(), dtype=str)
    totals = np.asarray(table["total_counts"].combine_chunks(), dtype=np.float64)
    targets = np.asarray(
        [target_index.get(label, -1) for label in labels], dtype=np.int32
    )
    keep_rows = targets >= 0
    targets = targets[keep_rows]
    totals = totals[keep_rows]
    result = np.zeros((DEV_COUNT, RESPONSE_COUNT), dtype=np.float64)
    cells = np.bincount(targets, minlength=DEV_COUNT).astype(np.int32)
    umi = np.bincount(targets, weights=totals, minlength=DEV_COUNT)
    if not np.any(keep_rows):
        return result, cells, umi

    token_lists = table["gene_token_id"].combine_chunks()
    value_lists = table["gene_expression"].combine_chunks()
    offsets = np.asarray(token_lists.offsets, dtype=np.int64)
    lengths = np.diff(offsets)
    flat_row = np.repeat(np.arange(len(labels), dtype=np.int64), lengths)
    all_values = np.asarray(value_lists.values, dtype=np.float64)
    expression_totals = np.bincount(flat_row, weights=all_values, minlength=len(labels))
    require(
        float(np.max(np.abs(expression_totals - totals))) <= 1e-6,
        "Orion total_counts differs from sparse expression sum",
    )
    flat_keep = keep_rows[flat_row]
    tokens = np.asarray(token_lists.values, dtype=np.int64)[flat_keep]
    values = all_values[flat_keep]
    row_target = np.asarray(
        [target_index.get(label, -1) for label in labels], dtype=np.int32
    )
    flat_targets = row_target[flat_row[flat_keep]]
    response_columns = token_to_response[tokens]
    keep_values = response_columns >= 0
    matrix = sp.coo_matrix(
        (
            values[keep_values],
            (flat_targets[keep_values], response_columns[keep_values]),
        ),
        shape=(DEV_COUNT, RESPONSE_COUNT),
        dtype=np.float64,
    ).tocsr()
    result[:] = matrix.toarray()
    return result, cells, umi


def stream_orion_batch_sum(
    path: Path,
    target_index: Mapping[str, int],
    token_to_response: np.ndarray,
    target_values: Sequence[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Aggregate one Orion file without materializing the full filtered table."""
    total_sum = np.zeros((DEV_COUNT, RESPONSE_COUNT), dtype=np.float64)
    total_cells = np.zeros(DEV_COUNT, dtype=np.int32)
    total_umi = np.zeros(DEV_COUNT, dtype=np.float64)
    lazy = (
        pl.scan_parquet(
            path,
            parallel="none",
            low_memory=True,
            rechunk=False,
            cache=False,
        )
        .filter(
            (pl.col("pass_guide_filter") == 1)
            & pl.col("gene_target").is_in(list(target_values))
        )
        .select(
            "gene_target",
            "total_counts",
            "gene_token_id",
            "gene_expression",
        )
    )
    for frame in lazy.collect_batches(
        chunk_size=ORION_SCAN_BATCH_SIZE,
        maintain_order=True,
        engine="streaming",
    ):
        table = frame.to_arrow()
        batch_sum, batch_cells, batch_umi = sparse_orion_batch_sum(
            table,
            target_index,
            token_to_response,
        )
        total_sum += batch_sum
        total_cells += batch_cells
        total_umi += batch_umi
        del table, frame, batch_sum, batch_cells, batch_umi
        release_arrow_memory()
    del lazy
    release_arrow_memory()
    return total_sum, total_cells, total_umi


def scan_orion_source(
    donor: str,
    assignments: np.ndarray,
    planned_batch_labels: Sequence[str],
    responses: Sequence[str],
    limit: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    line = DONOR_LINES[donor]
    targets = read_lines(DEV_TARGETS)
    target_index = {target: index for index, target in enumerate(targets)}
    genes = orion_gene_axis()
    require(
        genes == read_lines(DONOR_ROOTS[donor] / "genes.txt", unique=False),
        f"{donor}: gene axis changed",
    )
    last_token = {gene: token for token, gene in enumerate(genes)}
    token_to_response = np.full(len(genes), -1, dtype=np.int32)
    for response_index, gene in enumerate(responses):
        require(gene in last_token, f"{donor}: response gene absent: {gene}")
        token_to_response[last_token[gene]] = response_index
    files = numbered_orion_files(line)
    validate_batch_axis(tuple(path.stem for path in files), planned_batch_labels, donor)
    require(len(files) == assignments.shape[1], f"{donor}: batch count changed")
    sums = np.zeros((SPLITS, HALVES, DEV_COUNT, RESPONSE_COUNT), dtype=np.float64)
    cells = np.zeros((SPLITS, HALVES, DEV_COUNT), dtype=np.int32)
    umi = np.zeros((SPLITS, HALVES, DEV_COUNT), dtype=np.float64)
    target_values = list(targets)
    for batch, path in enumerate(files):
        batch_sum, batch_cells, batch_umi = stream_orion_batch_sum(
            path,
            target_index,
            token_to_response,
            target_values,
        )
        active = batch_cells > 0
        active_targets = np.flatnonzero(active)
        for split in range(SPLITS):
            half = int(assignments[split, batch])
            require(half >= 0, f"{donor}: active batch unassigned")
            sums[split, half, active_targets] += batch_sum[active_targets]
            cells[split, half, active_targets] += batch_cells[active_targets]
            umi[split, half, active_targets] += batch_umi[active_targets]
        del batch_sum, batch_cells, batch_umi
        release_arrow_memory()
        if (batch + 1) % 10 == 0 or batch + 1 == len(files):
            print(f"[{donor}] {batch + 1}/{len(files)}", flush=True)
            memory_check(f"{donor} source scan", limit)
    return sums, cells, umi


def finalize_source_profiles(
    sums: np.ndarray,
    cells: np.ndarray,
    umi: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    require(sums.shape == (SPLITS, HALVES, DEV_COUNT, RESPONSE_COUNT), "sum shape")
    profiles = np.log1p(
        TS
        * np.divide(
            sums,
            umi[..., None],
            out=np.zeros_like(sums),
            where=umi[..., None] > 0,
        )
    )
    centered = np.zeros_like(profiles, dtype=np.float32)
    support = cells > 0
    for split in range(SPLITS):
        for half in range(HALVES):
            use = support[split, half]
            require(np.any(use), "empty split half")
            mean = profiles[split, half, use].mean(axis=0)
            centered[split, half, use] = (profiles[split, half, use] - mean).astype(
                np.float32
            )
    pooled_sums = sums[0, 0] + sums[0, 1]
    pooled_umi = umi[0, 0] + umi[0, 1]
    pooled_cells = cells[0, 0] + cells[0, 1]
    for split in range(1, SPLITS):
        require(
            np.allclose(sums[split, 0] + sums[split, 1], pooled_sums, rtol=0, atol=0),
            "split pooled expression changed",
        )
        require(
            np.array_equal(cells[split, 0] + cells[split, 1], pooled_cells),
            "split pooled cells changed",
        )
        require(
            np.allclose(umi[split, 0] + umi[split, 1], pooled_umi, rtol=0, atol=0),
            "split pooled UMI changed",
        )
    pooled = np.log1p(
        TS
        * np.divide(
            pooled_sums,
            pooled_umi[:, None],
            out=np.zeros_like(pooled_sums),
            where=pooled_umi[:, None] > 0,
        )
    )
    pooled_support = pooled_cells > 0
    pooled_centered = np.zeros_like(pooled, dtype=np.float32)
    pooled_centered[pooled_support] = (
        pooled[pooled_support] - pooled[pooled_support].mean(axis=0)
    ).astype(np.float32)
    valid = np.all(cells >= MIN_CELLS_HALF, axis=1)
    return centered, valid, pooled.astype(np.float32), pooled_centered, pooled_cells


def compact_pooled_profile(
    donor: str,
    targets: Sequence[str],
    responses: Sequence[str],
) -> np.ndarray:
    root = DONOR_ROOTS[donor]
    perts = read_lines(root / "perts.txt")
    genes = read_lines(root / "genes.txt", unique=False)
    pert_index = {target: index for index, target in enumerate(perts)}
    gene_index = {gene: index for index, gene in enumerate(genes)}
    require(all(target in pert_index for target in targets), f"{donor}: target absent")
    require(all(gene in gene_index for gene in responses), f"{donor}: response absent")
    rows = np.asarray([pert_index[target] for target in targets], dtype=np.int64)
    columns = np.asarray([gene_index[gene] for gene in responses], dtype=np.int64)
    delta = np.load(root / "delta.npy", mmap_mode="r")
    control = np.load(root / "b_ctrl.npy", mmap_mode="r")
    return (
        np.asarray(delta[np.ix_(rows, columns)], dtype=np.float32)
        + np.asarray(control[columns], dtype=np.float32)[None, :]
    )


def source_build(donor: str, limit: float) -> None:
    require(donor in DONORS, f"unknown donor: {donor}")
    validate_donor_source_hashes(donor)
    output = SOURCE_BUILD_DIR / donor
    require(not output.exists(), f"source build exists: {output}")
    plan = np.load(PLAN_DIR / "plan-arrays.npz", allow_pickle=False)
    assignments = np.asarray(plan[f"{donor}_assignments"], dtype=np.int8)
    planned_batch_labels = tuple(map(str, plan[f"{donor}_batch_labels"]))
    expected_cells = np.asarray(plan[f"{donor}_half_cells"], dtype=np.int32)
    expected_umi = np.asarray(plan[f"{donor}_half_umi"], dtype=np.float64)
    require(
        np.all(plan[f"{donor}_half_batches"] >= MIN_BATCHES_HALF),
        f"{donor}: plan has a split half below the batch minimum",
    )
    responses = read_lines(RESPONSE_GENES)
    with donor_lock(donor):
        if donor == "k562_gwps":
            sums, cells, umi = scan_k562_source(
                assignments, planned_batch_labels, responses, limit
            )
        else:
            sums, cells, umi = scan_orion_source(
                donor, assignments, planned_batch_labels, responses, limit
            )
        require(
            np.array_equal(cells, expected_cells), f"{donor}: cell metadata mismatch"
        )
        require(
            np.allclose(umi, expected_umi, rtol=0, atol=1e-6),
            f"{donor}: UMI metadata mismatch",
        )
        centered, valid, pooled_log, pooled, pooled_cells = finalize_source_profiles(
            sums, cells, umi
        )
        compact = compact_pooled_profile(donor, read_lines(DEV_TARGETS), responses)
        pooled_error = np.abs(
            pooled_log.astype(np.float64) - compact.astype(np.float64)
        )
        require(
            float(pooled_error.max()) <= POOLED_RECON_ATOL,
            f"{donor}: pooled raw reconstruction exceeds tolerance",
        )
        receipt = {
            "schema_version": 1,
            "created_utc": provenance_utc(),
            "donor": donor,
            "source_only": True,
            "development_outcomes_opened": False,
            "sealed_410_effect_values_opened": False,
            "split_rows_fit_whitener": False,
            "targets": DEV_COUNT,
            "response_genes": RESPONSE_COUNT,
            "split_pairs": SPLITS,
            "valid_pairs": {
                "min": int(valid.sum(axis=0).min()),
                "median": float(np.median(valid.sum(axis=0))),
                "max": int(valid.sum(axis=0).max()),
                "fraction_all_eight": float(np.mean(valid.sum(axis=0) == SPLITS)),
                "fraction_at_least_six": float(
                    np.mean(valid.sum(axis=0) >= MIN_VALID_PAIRS)
                ),
            },
            "pooled_cells": {
                "min": int(pooled_cells.min()),
                "median": float(np.median(pooled_cells)),
                "max": int(pooled_cells.max()),
            },
            "pooled_compact_reconstruction": {
                "tolerance": POOLED_RECON_ATOL,
                "max_abs": float(pooled_error.max()),
                "mean_abs": float(pooled_error.mean()),
                "within_tolerance": True,
            },
            "peak_rss_gib": peak_rss_gib(),
        }
        with atomic_stage_directory(output) as staging:
            save_npy(staging / "half-centered.npy", centered)
            save_npy(staging / "valid-pairs.npy", valid)
            save_npy(staging / "half-cells.npy", cells)
            save_npy(staging / "half-umi.npy", umi)
            save_npy(staging / "pooled-centered.npy", pooled)
            save_npy(staging / "pooled-cells.npy", pooled_cells)
            save_npy(staging / "pooled-umi.npy", umi[0].sum(axis=0))
            save_npy(
                staging / "pooled-zero-fraction.npy",
                np.mean(pooled_log == 0.0, axis=1).astype(np.float32),
            )
            atomic_json(staging / "receipt.json", receipt)
            atomic_json(staging / "artifact-manifest.json", directory_manifest(staging))
        print(json.dumps(receipt, indent=2, sort_keys=True))


def average_tie_ranks(values: np.ndarray) -> np.ndarray:
    data = np.asarray(values, dtype=np.float64)
    require(data.ndim == 1, "rank input must be one-dimensional")
    order = np.argsort(data, kind="mergesort")
    ranks = np.empty(len(data), dtype=np.float64)
    start = 0
    while start < len(data):
        end = start + 1
        while end < len(data) and data[order[end]] == data[order[start]]:
            end += 1
        ranks[order[start:end]] = 0.5 * (start + end - 1) + 1.0
        start = end
    return ranks


def midrank_z(values: np.ndarray, valid: np.ndarray | None = None) -> np.ndarray:
    data = np.asarray(values, dtype=np.float64)
    mask = np.isfinite(data) if valid is None else np.asarray(valid, dtype=bool)
    mask &= np.isfinite(data)
    result = np.full(len(data), np.nan, dtype=np.float64)
    count = int(mask.sum())
    require(count >= 2, "fewer than two valid values for midrank z")
    ranks = average_tie_ranks(data[mask])
    result[mask] = 2.0 * (ranks - 1.0) / (count - 1.0) - 1.0
    return result


def gated_midrank_z(values: np.ndarray, gate: np.ndarray) -> np.ndarray:
    """Rank only absolute-null survivors; null-level targets remain identity."""
    data = np.asarray(values, dtype=np.float64)
    mask = np.asarray(gate, dtype=bool) & np.isfinite(data)
    result = np.zeros(len(data), dtype=np.float64)
    count = int(mask.sum())
    if count == 1:
        result[mask] = 1.0
    elif count >= 2:
        result[mask] = midrank_z(data[mask])
    return result


def load_operator(role: str) -> tuple[np.ndarray, np.ndarray]:
    require(role in ("k562", "non_k562"), f"unknown operator role: {role}")
    path = EXP140_IDENTITY / f"whitener-{role}.npz"
    data = np.load(path, allow_pickle=False)
    evals = np.asarray(data["evals"], dtype=np.float64)
    evecs = np.asarray(data["evecs"], dtype=np.float32)
    require(evals.shape == (RESPONSE_COUNT,), "operator eigenvalue shape changed")
    require(
        evecs.shape == (RESPONSE_COUNT, RESPONSE_COUNT), "operator basis shape changed"
    )
    return evals, evecs


def spectral_coordinates(
    matrix: np.ndarray,
    evals: np.ndarray,
    evecs: np.ndarray,
    *,
    chunk: int = 64,
) -> np.ndarray:
    """Return coordinates equivalent for cosine to the full whitened vectors."""
    values = np.asarray(matrix, dtype=np.float32)
    output = np.empty(values.shape, dtype=np.float32)
    weights = evals ** (-WHITEN_GAMMA / 2.0)
    for lo in range(0, len(values), chunk):
        hi = min(lo + chunk, len(values))
        output[lo:hi] = (
            np.asarray(values[lo:hi], dtype=np.float64) @ evecs * weights[None, :]
        ).astype(np.float32)
    return output


def cosine_matrix(query: np.ndarray, gallery: np.ndarray) -> np.ndarray:
    q = np.asarray(query, dtype=np.float64)
    g = np.asarray(gallery, dtype=np.float64)
    q /= np.maximum(np.linalg.norm(q, axis=1, keepdims=True), EPS)
    g /= np.maximum(np.linalg.norm(g, axis=1, keepdims=True), EPS)
    return q @ g.T


def pds_from_similarity(similarity: np.ndarray) -> np.ndarray:
    result = np.empty_like(similarity)
    denominator = max(similarity.shape[1] - 1, 1)
    for row in range(similarity.shape[0]):
        values = similarity[row]
        order = np.argsort(-values, kind="mergesort")
        start = 0
        while start < len(values):
            end = start + 1
            while end < len(values) and values[order[end]] == values[order[start]]:
                end += 1
            midrank = 1.0 + start + 0.5 * (end - start - 1)
            result[row, order[start:end]] = 1.0 - (midrank - 1.0) / denominator
            start = end
    return result


def pds_matrix(query: np.ndarray, gallery: np.ndarray) -> np.ndarray:
    return pds_from_similarity(cosine_matrix(query, gallery))


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    x_rank = average_tie_ranks(np.asarray(x, dtype=np.float64))
    y_rank = average_tie_ranks(np.asarray(y, dtype=np.float64))
    x_rank -= x_rank.mean()
    y_rank -= y_rank.mean()
    denominator = float(np.linalg.norm(x_rank) * np.linalg.norm(y_rank))
    return 0.0 if denominator <= EPS else float(x_rank @ y_rank / denominator)


def bootstrap_spearman(
    x: np.ndarray, y: np.ndarray, draws: int, seed: int
) -> np.ndarray:
    x_values = np.asarray(x, dtype=np.float64)
    y_values = np.asarray(y, dtype=np.float64)
    rng = np.random.default_rng(seed)
    result = np.empty(draws, dtype=np.float64)
    chunk = 256
    for lo in range(0, draws, chunk):
        hi = min(lo + chunk, draws)
        indices = rng.integers(0, len(x_values), size=(hi - lo, len(x_values)))
        xb = rankdata(x_values[indices], axis=1, method="average")
        yb = rankdata(y_values[indices], axis=1, method="average")
        xb -= xb.mean(axis=1, keepdims=True)
        yb -= yb.mean(axis=1, keepdims=True)
        numerator = np.sum(xb * yb, axis=1)
        denominator = np.linalg.norm(xb, axis=1) * np.linalg.norm(yb, axis=1)
        result[lo:hi] = np.divide(
            numerator,
            denominator,
            out=np.zeros_like(numerator),
            where=denominator > EPS,
        )
    return result


def rank_scaled(values: np.ndarray) -> np.ndarray:
    ranks = average_tie_ranks(np.asarray(values, dtype=np.float64))
    return (ranks - 1.0) / max(len(ranks) - 1.0, 1.0)


def nearest_feature_blocks(
    features: np.ndarray, nominal_size: int = MATCH_BLOCK_SIZE
) -> tuple[np.ndarray, ...]:
    values = np.asarray(features, dtype=np.float64)
    require(values.ndim == 2 and len(values) >= 2, "invalid matching features")
    remaining = list(range(len(values)))
    blocks: list[np.ndarray] = []
    while remaining:
        if len(remaining) < 2:
            blocks[-1] = np.concatenate(
                [blocks[-1], np.asarray(remaining, dtype=np.int64)]
            )
            break
        anchor = remaining[0]
        candidates = np.asarray(remaining, dtype=np.int64)
        distance = np.sum((values[candidates] - values[anchor]) ** 2, axis=1)
        order = np.lexsort((candidates, distance))
        take = candidates[order[: min(nominal_size, len(candidates))]]
        if len(remaining) - len(take) == 1:
            take = candidates[order[: len(take) + 1]]
        chosen = set(map(int, take))
        remaining = [value for value in remaining if value not in chosen]
        blocks.append(np.asarray(sorted(chosen), dtype=np.int64))
    require(all(len(block) >= 2 for block in blocks), "singleton matching block")
    flattened = np.concatenate(blocks)
    require(
        np.array_equal(np.sort(flattened), np.arange(len(values))), "block coverage"
    )
    return tuple(blocks)


def stable_derangement(rng: np.random.Generator, size: int) -> np.ndarray:
    base = np.arange(size, dtype=np.int64)
    for _ in range(100_000):
        candidate = rng.permutation(base)
        if np.all(candidate != base):
            return candidate
    raise ContractError("failed to generate derangement")


def matched_null_permutations(
    features: np.ndarray, draws: int, seed: int
) -> tuple[np.ndarray, tuple[np.ndarray, ...]]:
    blocks = nearest_feature_blocks(features)
    rng = np.random.default_rng(seed)
    result = np.empty((draws, len(features)), dtype=np.int64)
    base = np.arange(len(features), dtype=np.int64)
    for draw in range(draws):
        permutation = base.copy()
        for block in blocks:
            permutation[block] = block[stable_derangement(rng, len(block))]
        require(np.all(permutation != base), "matched null contains a fixed point")
        result[draw] = permutation
    return result, blocks


def evaluate_role(
    donor: str,
    role: str,
    centered: np.ndarray,
    valid: np.ndarray,
    pooled: np.ndarray,
    pooled_cells: np.ndarray,
    pooled_umi: np.ndarray,
    pooled_zero_fraction: np.ndarray,
    occupied_batches: np.ndarray,
    limit: float,
) -> dict[str, Any]:
    evals, evecs = load_operator(role)
    pooled_spectral = spectral_coordinates(pooled, evals, evecs)
    pooled_norm = np.linalg.norm(pooled_spectral, axis=1)
    features = np.column_stack(
        [
            rank_scaled(np.log1p(pooled_cells)),
            rank_scaled(occupied_batches),
            rank_scaled(np.log1p(pooled_norm)),
        ]
    )
    null_permutations, blocks = matched_null_permutations(
        features,
        NULL_DRAWS,
        stable_uint64(NULL_SEED, donor, role),
    )
    directional = np.full((SPLITS, HALVES, DEV_COUNT), np.nan, dtype=np.float64)
    nearest_decoy_margin = np.full(
        (SPLITS, HALVES, DEV_COUNT), np.nan, dtype=np.float64
    )
    null_directional = np.full(
        (NULL_DRAWS, SPLITS, HALVES, DEV_COUNT), np.nan, dtype=np.float64
    )
    null_sum = np.zeros(NULL_DRAWS, dtype=np.float64)
    null_count = np.zeros(NULL_DRAWS, dtype=np.int64)
    all_indices = np.arange(DEV_COUNT, dtype=np.int64)
    for split in range(SPLITS):
        pair_valid = np.asarray(valid[split], dtype=bool)
        indices = all_indices[pair_valid]
        require(len(indices) >= 2, f"{donor}/{role}/split{split}: insufficient gallery")
        transformed = [
            spectral_coordinates(centered[split, half, indices], evals, evecs)
            for half in range(HALVES)
        ]
        for direction, (query, gallery) in enumerate(
            ((transformed[0], transformed[1]), (transformed[1], transformed[0]))
        ):
            similarity = cosine_matrix(query, gallery)
            matrix = pds_from_similarity(similarity)
            directional[split, direction, indices] = np.diag(matrix)
            decoys = similarity.copy()
            np.fill_diagonal(decoys, -np.inf)
            nearest_decoy_margin[split, direction, indices] = np.diag(
                similarity
            ) - np.max(decoys, axis=1)
            local_position = np.full(DEV_COUNT, -1, dtype=np.int64)
            local_position[indices] = np.arange(len(indices), dtype=np.int64)
            for draw in range(NULL_DRAWS):
                designated_global = null_permutations[draw, indices]
                designated_local = local_position[designated_global]
                usable = designated_local >= 0
                if np.any(usable):
                    rows = np.flatnonzero(usable)
                    values = matrix[rows, designated_local[usable]]
                    global_rows = indices[rows]
                    null_directional[draw, split, direction, global_rows] = values
                    null_sum[draw] += float(values.sum())
                    null_count[draw] += int(usable.sum())
        del transformed, similarity, decoys
        memory_check(f"{donor}/{role} split {split + 1}", limit)

    directional_flat = directional.reshape(SPLITS * HALVES, DEV_COUNT)
    directional_count = np.sum(np.isfinite(directional_flat), axis=0)
    reliability = np.nanquantile(directional_flat, 0.25, axis=0)
    reliability[directional_count < 2 * MIN_VALID_PAIRS] = np.nan
    null_reliability = np.nanquantile(
        null_directional.reshape(NULL_DRAWS, SPLITS * HALVES, DEV_COUNT),
        0.25,
        axis=1,
    )
    null_reliability[:, directional_count < 2 * MIN_VALID_PAIRS] = np.nan
    target_null_p95 = np.nanquantile(null_reliability, 0.95, axis=0)
    target_reliability_pass = (
        np.isfinite(reliability)
        & np.isfinite(target_null_p95)
        & (reliability > target_null_p95)
    )
    reliability_z = gated_midrank_z(reliability, target_reliability_pass)
    target_margin = np.nanquantile(
        nearest_decoy_margin.reshape(SPLITS * HALVES, DEV_COUNT), 0.25, axis=0
    )
    first = np.nanquantile(directional[:4].reshape(8, DEV_COUNT), 0.25, axis=0)
    second = np.nanquantile(directional[4:].reshape(8, DEV_COUNT), 0.25, axis=0)
    rank_valid = np.isfinite(first) & np.isfinite(second)
    require(int(rank_valid.sum()) >= 2, f"{donor}/{role}: no rank-stability targets")
    observed_rank_correlation = spearman(first[rank_valid], second[rank_valid])
    bootstrap = bootstrap_spearman(
        first[rank_valid],
        second[rank_valid],
        BOOTSTRAP_DRAWS,
        stable_uint64(BOOTSTRAP_SEED, donor, role),
    )
    null_values = np.divide(
        null_sum,
        null_count,
        out=np.full(NULL_DRAWS, np.nan, dtype=np.float64),
        where=null_count > 0,
    )
    observed_retrieval = float(np.nanmean(directional))
    null_p95 = float(np.nanquantile(null_values, 0.95))
    valid_pairs = np.sum(valid, axis=0)
    diagnostic_values = {
        "log1p_cell_count": np.log1p(pooled_cells),
        "log1p_total_umi": np.log1p(pooled_umi),
        "occupied_batch_count": occupied_batches,
        "log1p_pooled_effect_norm": np.log1p(pooled_norm),
        "pooled_zero_fraction": pooled_zero_fraction,
        "nearest_decoy_margin": target_margin,
    }
    finite_reliability = np.isfinite(reliability)
    correlations = {
        name: spearman(
            reliability[finite_reliability], np.asarray(values)[finite_reliability]
        )
        for name, values in diagnostic_values.items()
    }
    report = {
        "donor": donor,
        "operator_role": role,
        "source_only": True,
        "split_rows_fit_whitener": False,
        "valid_pairs": {
            "min": int(valid_pairs.min()),
            "median": float(np.median(valid_pairs)),
            "max": int(valid_pairs.max()),
            "fraction_all_eight": float(np.mean(valid_pairs == SPLITS)),
            "fraction_at_least_six": float(np.mean(valid_pairs >= MIN_VALID_PAIRS)),
        },
        "rank_stability": {
            "targets": int(rank_valid.sum()),
            "spearman": observed_rank_correlation,
            "bootstrap_draws": BOOTSTRAP_DRAWS,
            "lower_95": float(np.quantile(bootstrap, 0.025)),
            "upper_95": float(np.quantile(bootstrap, 0.975)),
        },
        "same_target_retrieval": {
            "observed_mean_pds": observed_retrieval,
            "matched_null_draws": NULL_DRAWS,
            "matched_null_p95": null_p95,
            "margin_over_null_p95": observed_retrieval - null_p95,
            "matching_blocks": len(blocks),
            "maximum_matching_block_size": max(map(len, blocks)),
        },
        "target_level_absolute_gate": {
            "passed_targets": int(target_reliability_pass.sum()),
            "fraction": float(np.mean(target_reliability_pass)),
            "null_draws": NULL_DRAWS,
            "strictly_above_target_matched_null_p95": True,
        },
        "technical_correlations": correlations,
    }
    report["passes"] = bool(
        report["valid_pairs"]["fraction_all_eight"] >= ALL_EIGHT_FRACTION
        and report["valid_pairs"]["fraction_at_least_six"] == 1.0
        and report["rank_stability"]["lower_95"] > 0.0
        and observed_retrieval > null_p95
    )
    return {
        "report": report,
        "directional_pds": directional,
        "reliability": reliability,
        "reliability_z": reliability_z,
        "target_null_p95": target_null_p95,
        "target_reliability_pass": target_reliability_pass,
        "cell_count_z": midrank_z(pooled_cells.astype(np.float64)),
        "total_umi_z": midrank_z(pooled_umi.astype(np.float64)),
        "occupied_batch_count_z": midrank_z(occupied_batches.astype(np.float64)),
        "pooled_norm": pooled_norm,
        "pooled_norm_z": midrank_z(pooled_norm),
        "pooled_zero_fraction": pooled_zero_fraction,
        "nearest_decoy_margin": target_margin,
        "rank_bootstrap": bootstrap,
        "matched_null": null_values,
        "target_matched_null_p95": target_null_p95,
        "matched_null_permutations": null_permutations.astype(np.int16),
    }


def source_evaluate(donor: str, limit: float) -> None:
    require(donor in DONORS, f"unknown donor: {donor}")
    validate_preoutcome()
    validate_source_manifest_paths(
        tuple(
            EXP140_IDENTITY / f"whitener-{role}.npz" for role in OPERATOR_ROLES[donor]
        )
    )
    source_dir = SOURCE_BUILD_DIR / donor
    validate_directory_manifest(source_dir)
    output = SOURCE_EVAL_DIR / donor
    require(not output.exists(), f"source evaluation exists: {output}")
    centered = np.load(source_dir / "half-centered.npy", mmap_mode="r")
    valid = np.asarray(np.load(source_dir / "valid-pairs.npy"), dtype=bool)
    pooled = np.asarray(np.load(source_dir / "pooled-centered.npy"), dtype=np.float32)
    pooled_cells = np.asarray(np.load(source_dir / "pooled-cells.npy"), dtype=np.int64)
    pooled_umi = np.asarray(np.load(source_dir / "pooled-umi.npy"), dtype=np.float64)
    pooled_zero_fraction = np.asarray(
        np.load(source_dir / "pooled-zero-fraction.npy"), dtype=np.float32
    )
    plan = np.load(PLAN_DIR / "plan-arrays.npz", allow_pickle=False)
    occupied_batches = np.asarray(plan[f"{donor}_occupied_batches"], dtype=np.int32)
    reports = {}
    with atomic_stage_directory(output) as staging:
        for role in OPERATOR_ROLES[donor]:
            result = evaluate_role(
                donor,
                role,
                centered,
                valid,
                pooled,
                pooled_cells,
                pooled_umi,
                pooled_zero_fraction,
                occupied_batches,
                limit,
            )
            reports[role] = result.pop("report")
            atomic_npz(staging / f"{role}.npz", **result)
        summary = {
            "schema_version": 1,
            "created_utc": provenance_utc(),
            "donor": donor,
            "source_only": True,
            "development_outcomes_opened": False,
            "sealed_410_effect_values_opened": False,
            "roles": reports,
            "passes": all(report["passes"] for report in reports.values()),
            "peak_rss_gib": peak_rss_gib(),
        }
        atomic_json(staging / "summary.json", summary)
        atomic_json(staging / "artifact-manifest.json", directory_manifest(staging))
    print(json.dumps(summary, indent=2, sort_keys=True))


def normalized_weights(
    z_by_donor: Mapping[str, np.ndarray],
    donors: Sequence[str],
    eta: float,
) -> np.ndarray:
    z = np.column_stack([np.asarray(z_by_donor[donor]) for donor in donors])
    complete = np.all(np.isfinite(z), axis=1)
    raw = np.ones_like(z, dtype=np.float64)
    raw[complete] = 1.0 + eta * z[complete]
    require(np.all(raw >= 0.5 - EPS), "weight fell below frozen lower bound")
    total = raw.sum(axis=1, keepdims=True)
    require(np.all(total > 0), "target has no available donor")
    result = raw / total
    common = complete & np.all(z == z[:, :1], axis=1)
    result[~complete | common] = 1.0 / len(donors)
    return result


def source_gate() -> None:
    validate_preoutcome()
    require(not SOURCE_GATE_DIR.exists(), f"source gate exists: {SOURCE_GATE_DIR}")
    for donor in DONORS:
        validate_directory_manifest(SOURCE_BUILD_DIR / donor)
        validate_directory_manifest(SOURCE_EVAL_DIR / donor)
    donor_summaries = {
        donor: read_json(SOURCE_EVAL_DIR / donor / "summary.json") for donor in DONORS
    }
    role_arrays: dict[tuple[str, str], dict[str, np.ndarray]] = {}
    for donor in DONORS:
        for role in OPERATOR_ROLES[donor]:
            data = np.load(SOURCE_EVAL_DIR / donor / f"{role}.npz", allow_pickle=False)
            role_arrays[(donor, role)] = {
                key: np.asarray(data[key]) for key in data.files
            }
    passes = all(summary["passes"] for summary in donor_summaries.values())
    result = {
        "schema_version": 1,
        "created_utc": provenance_utc(),
        "source_only": True,
        "development_outcomes_opened": False,
        "sealed_410_effect_values_opened": False,
        "donors": donor_summaries,
        "source_gate_passes": passes,
        "common_score_identity_control": True,
        "development_authorized": False,
        "production_integration_authorized": False,
        "count_generation_authorized": False,
        "submission_authorized": False,
        "peak_rss_gib": peak_rss_gib(),
    }
    with atomic_stage_directory(SOURCE_GATE_DIR) as staging:
        for outer, donors in OUTER_DONORS.items():
            payload: dict[str, np.ndarray] = {
                "etas": np.asarray(ETAS, dtype=np.float64)
            }
            for signal in (
                "reliability_z",
                "cell_count_z",
                "total_umi_z",
                "occupied_batch_count_z",
                "pooled_norm_z",
            ):
                z_by_donor = {
                    donor: role_arrays[(donor, outer)][signal] for donor in donors
                }
                payload[signal] = np.stack(
                    [normalized_weights(z_by_donor, donors, eta) for eta in ETAS],
                    axis=0,
                )
            common = np.linspace(-1.0, 1.0, DEV_COUNT)
            common_weights = normalized_weights(
                {donor: common for donor in donors}, donors, ETAS[-1]
            )
            identity = np.full_like(common_weights, 1.0 / len(donors))
            require(
                np.array_equal(common_weights, identity),
                f"{outer}: common-score identity control failed",
            )
            payload["donor_assignment_nulls"] = np.asarray(
                np.load(PLAN_DIR / "plan-arrays.npz", allow_pickle=False)[
                    f"donor_assignment_null_{outer}"
                ],
                dtype=np.int8,
            )
            payload["donors"] = np.asarray(donors, dtype="U32")
            atomic_npz(staging / f"candidate-weights-{outer}.npz", **payload)
        atomic_json(staging / "result.json", result)
        atomic_json(staging / "artifact-manifest.json", directory_manifest(staging))
    print(json.dumps(result, indent=2, sort_keys=True))


def validate_development_authorization() -> dict[str, Any]:
    """Validate a future additive authorization before any outcome import."""
    validate_preoutcome()
    validate_directory_manifest(SOURCE_GATE_DIR)
    gate = read_json(SOURCE_GATE_DIR / "result.json")
    require(gate["source_gate_passes"] is True, "source-only gate did not pass")
    require(
        DEVELOPMENT_AUTHORIZATION.is_file(),
        "separate DEVELOPMENT_ACCESS_AUTHORIZATION.v1.json is absent",
    )
    authorization = read_json(DEVELOPMENT_AUTHORIZATION)
    require(authorization["development_613_authorized"] is True, "development disabled")
    require(
        authorization["untouched_410_test_authorized"] is False, "sealed test enabled"
    )
    require(authorization["count_generation_authorized"] is False, "counts enabled")
    require(
        authorization["challenge_build_authorized"] is False, "challenge build enabled"
    )
    require(authorization["submission_authorized"] is False, "submission enabled")
    require(
        authorization["protocol_sha256"] == sha256_file(PROTOCOL),
        "authorization protocol binding changed",
    )
    require(
        authorization["plan_manifest_sha256"]
        == sha256_file(PLAN_DIR / "plan-manifest.json"),
        "authorization plan binding changed",
    )
    require(
        authorization["preoutcome_authorization_sha256"]
        == sha256_file(PREOUTCOME_AUTHORIZATION),
        "authorization preoutcome binding changed",
    )
    require(
        authorization["development_interface_sha256"]
        == sha256_file(DEVELOPMENT_INTERFACE),
        "authorization interface binding changed",
    )
    require(
        authorization["source_gate_manifest_sha256"]
        == sha256_file(SOURCE_GATE_DIR / "artifact-manifest.json"),
        "authorization source-gate binding changed",
    )
    return authorization


def develop_interface_only() -> None:
    validate_development_authorization()
    require(not DEVELOPMENT_LOCK.exists(), "development access lock already exists")
    raise ContractError(
        "authorization validated, but this pre-outcome implementation intentionally "
        "contains no held-development outcome loader; add and review a separate "
        "development-stage implementation before access"
    )


def audit() -> None:
    payload: dict[str, Any] = {"plan": validate_plan()}
    if PREOUTCOME_AUTHORIZATION.exists():
        payload["preoutcome"] = validate_preoutcome()
    for donor in DONORS:
        build = SOURCE_BUILD_DIR / donor
        evaluation = SOURCE_EVAL_DIR / donor
        if build.exists():
            validate_directory_manifest(build)
            payload.setdefault("source_builds", []).append(donor)
        if evaluation.exists():
            validate_directory_manifest(evaluation)
            payload.setdefault("source_evaluations", []).append(donor)
    if SOURCE_GATE_DIR.exists():
        validate_directory_manifest(SOURCE_GATE_DIR)
        payload["source_gate"] = read_json(SOURCE_GATE_DIR / "result.json")
    print(json.dumps(payload, indent=2, sort_keys=True))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "command",
        choices=(
            "plan",
            "freeze",
            "source-build",
            "source-evaluate",
            "source-gate",
            "develop",
            "audit",
        ),
    )
    parser.add_argument("--donor", choices=DONORS)
    parser.add_argument("--max-rss-gib", type=float, default=17.5)
    args = parser.parse_args()
    if args.command == "plan":
        require(args.donor is None, "plan does not accept --donor")
        create_plan()
    elif args.command == "freeze":
        require(args.donor is None, "freeze does not accept --donor")
        freeze()
    elif args.command == "source-build":
        require(args.donor is not None, "source-build requires --donor")
        with heavy_stage_lock(args.command, args.donor):
            source_build(args.donor, args.max_rss_gib)
    elif args.command == "source-evaluate":
        require(args.donor is not None, "source-evaluate requires --donor")
        with heavy_stage_lock(args.command, args.donor):
            source_evaluate(args.donor, args.max_rss_gib)
    elif args.command == "source-gate":
        require(args.donor is None, "source-gate does not accept --donor")
        with heavy_stage_lock(args.command, None):
            source_gate()
    elif args.command == "develop":
        require(args.donor is None, "develop does not accept --donor")
        develop_interface_only()
    else:
        audit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

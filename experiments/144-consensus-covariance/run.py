#!/usr/bin/env python3
"""Source-only v2 scaffold for Experiment 144.

The only mutable stage is ``prepare``. It fits cross-donor, matched-target
consensus covariance operators from explicitly allowed source rows. ``audit``
checks the frozen protocol, command surface, v1 quarantine, and any completed
v2 source-only artifact without opening biological outcome arrays.

There is deliberately no development, held-out, count, scoring, build, or
submission command in this file.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
import runpy
import select
import shutil
import tempfile
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, Mapping, Protocol, Sequence

import numpy as np
import scipy.linalg


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN_PATH = Path(__file__).resolve()
README_PATH = HERE / "README.md"
TEST_PATH = HERE / "test_run.py"
PROTOCOL_PATH = HERE / "preregistered_protocol.json"

PREPARE_AUTHORIZATION = HERE / "SOURCE_ONLY_PREPARE_AUTHORIZATION.v2.json"
OUTPUT_DIR = HERE / "work/source-only-v2"
PREPARE_LOCK = HERE / "work/source-only-v2.lock"
MEMORY_GUARD_SCRIPT = HERE / "memory_guard.py"

V1_AUTHORIZATION = HERE / "SOURCE_ONLY_PREPARE_AUTHORIZATION.v1.json"
V1_OUTPUT_DIR = HERE / "work/preoutcome-v1"

EXP132_PLAN = (
    ROOT / "experiments/132-response-decomposition-ridge/work/preregistered-freeze-v2"
)
ALL_TARGETS = EXP132_PLAN / "all_targets.txt"
RESPONSE_GENES = EXP132_PLAN / "response_genes.txt"
PANEL_TARGETS = ROOT / "data/challenge/pert_counts.csv"
ESSENTIAL_TARGETS = ROOT / "data/donors/k562_essential/perts.txt"

DONOR_ROOTS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
GROUP_DONORS = {
    "k562": ("orion_hct116", "orion_hek293t"),
    "non_k562": ("k562_gwps", "orion_hct116", "orion_hek293t"),
}

COMMANDS = ("prepare", "audit")
CANDIDATES = (
    {"name": "identity", "gamma": 0.0, "eta": 0.0},
    {"name": "consensus:g=0.65:eta=0.25", "gamma": 0.65, "eta": 0.25},
)

RESPONSE_COUNT = 4_992
WHITEN_GAMMA = 0.65
BLEND_ETA = 0.25
WHITEN_FLOOR_RELATIVE = 1e-12
NUMERICAL_RANK_EPS_MULTIPLIER = 10.0
MAX_RSS_CEILING_GIB = 18.0
MIN_FREE_FLOOR_PERCENT = 25.0
EPS = 1e-12
GIB = 1024**3
MANIFEST_NAME = "artifact-manifest.json"
RESULT_NAME = "result.json"
CANDIDATE_SPEC_NAME = "candidate-spec.json"
AUTHORIZATION_CONFIRMATION = "AUTHORIZE_SOURCE_ONLY_PREPARE_V2"
EXTERNAL_GUARD_ACTIVE_ENV = "VCC144_MEMORY_GUARD_ACTIVE"
EXTERNAL_GUARD_PATH_ENV = "VCC144_MEMORY_GUARD_PATH"
EXTERNAL_GUARD_MAX_RSS_ENV = "VCC144_MEMORY_GUARD_MAX_RSS_GIB"
EXTERNAL_GUARD_MIN_FREE_ENV = "VCC144_MEMORY_GUARD_MIN_FREE_PERCENT"
EXTERNAL_GUARD_REQUEST_FD_ENV = "VCC144_MEMORY_GUARD_REQUEST_FD"
EXTERNAL_GUARD_RESPONSE_FD_ENV = "VCC144_MEMORY_GUARD_RESPONSE_FD"
EXTERNAL_GUARD_TOKEN_ENV = "VCC144_MEMORY_GUARD_TOKEN"
EXTERNAL_GUARD_HANDSHAKE_TIMEOUT_SECONDS = 10.0

FORBIDDEN_FUNCTION_NAMES = {
    "build_submission",
    "develop",
    "evaluate",
    "load_confirmation_truth",
    "load_development_truth",
    "load_test_truth",
    "score",
    "submit",
}
FORBIDDEN_IMPORT_ROOTS = {
    "anndata",
    "h5py",
    "scanpy",
    "subprocess",
}
FORBIDDEN_SOURCE_PATH_PARTS = (
    "/context_A.h5ad",
    "/context_B.h5ad",
    "/context_C.h5ad",
    "/development-preoutcome",
    "/development-result",
    "/kolf",
    "/reserved",
    "/submissions/",
)


class ContractError(RuntimeError):
    """Raised when a frozen experiment boundary is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


class Guard(Protocol):
    """Minimal memory-guard interface used by every source data path."""

    def check(self, label: str) -> None: ...


@dataclass(frozen=True)
class MemorySnapshot:
    process_tree_rss_gib: float
    free_memory_percent: float


class MemoryGuard:
    """Fail-closed process-tree RSS and free-memory guard."""

    def __init__(self, max_rss_gib: float, min_free_percent: float) -> None:
        require(
            0.0 < max_rss_gib <= MAX_RSS_CEILING_GIB,
            f"max RSS must be in (0, {MAX_RSS_CEILING_GIB}] GiB",
        )
        require(
            MIN_FREE_FLOOR_PERCENT <= min_free_percent < 100.0,
            f"minimum free memory must be in [{MIN_FREE_FLOOR_PERCENT}, 100) percent",
        )
        self.max_rss_gib = float(max_rss_gib)
        self.min_free_percent = float(min_free_percent)
        self.peak_process_tree_rss_gib = 0.0
        self.minimum_free_memory_percent = 100.0

    @staticmethod
    def snapshot() -> MemorySnapshot:
        helpers = memory_guard_helpers()
        rss, process_count = helpers["process_rss"](os.getpid())
        require(process_count > 0, "memory guard lost the current process")
        free_percent = helpers["free_percent"]()
        require(free_percent is not None, "cannot measure system free memory")
        return MemorySnapshot(
            process_tree_rss_gib=float(rss) / GIB,
            free_memory_percent=float(free_percent),
        )

    def check(self, label: str) -> None:
        snapshot = self.snapshot()
        self.peak_process_tree_rss_gib = max(
            self.peak_process_tree_rss_gib,
            snapshot.process_tree_rss_gib,
        )
        self.minimum_free_memory_percent = min(
            self.minimum_free_memory_percent,
            snapshot.free_memory_percent,
        )
        print(
            f"[memory] {label}: tree_rss={snapshot.process_tree_rss_gib:.3f} GiB "
            f"free={snapshot.free_memory_percent:.1f}%",
            flush=True,
        )
        if snapshot.process_tree_rss_gib >= self.max_rss_gib:
            raise MemoryError(
                f"{label}: process-tree RSS "
                f"{snapshot.process_tree_rss_gib:.3f} GiB reached "
                f"{self.max_rss_gib:.3f} GiB"
            )
        if snapshot.free_memory_percent < self.min_free_percent:
            raise MemoryError(
                f"{label}: free memory {snapshot.free_memory_percent:.1f}% "
                f"fell below {self.min_free_percent:.1f}%"
            )

    def report(self) -> dict[str, float]:
        return {
            "max_rss_gib": self.max_rss_gib,
            "min_free_percent": self.min_free_percent,
            "peak_process_tree_rss_gib": self.peak_process_tree_rss_gib,
            "minimum_free_memory_percent": self.minimum_free_memory_percent,
        }


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"expected JSON object: {path}")
    return value


_MEMORY_GUARD_HELPERS: dict[str, Any] | None = None


def memory_guard_helpers() -> dict[str, Any]:
    global _MEMORY_GUARD_HELPERS
    if _MEMORY_GUARD_HELPERS is None:
        namespace = runpy.run_path(
            str(MEMORY_GUARD_SCRIPT),
            run_name="experiment144_memory_guard",
        )
        require(
            callable(namespace.get("process_rss")), "missing process-tree RSS helper"
        )
        require(callable(namespace.get("free_percent")), "missing free-memory helper")
        _MEMORY_GUARD_HELPERS = namespace
    return _MEMORY_GUARD_HELPERS


def read_lines(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    require(bool(values), f"empty axis: {path}")
    if unique:
        require(len(values) == len(set(values)), f"duplicate axis: {path}")
    return values


def read_csv_first_column(path: Path) -> tuple[str, ...]:
    lines = read_lines(path, unique=False)
    require(len(lines) > 1, f"empty CSV: {path}")
    values = tuple(line.split(",", 1)[0].strip() for line in lines[1:])
    require(len(values) == len(set(values)), f"duplicate CSV axis: {path}")
    return values


def sha256_file(path: Path, guard: Guard, block: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        block_index = 0
        while piece := handle.read(block):
            digest.update(piece)
            block_index += 1
            if block_index % 16 == 0:
                guard.check(f"hash {path.name}")
    guard.check(f"hash complete {path.name}")
    return digest.hexdigest()


def array_sha256(array: np.ndarray) -> str:
    value = np.ascontiguousarray(array)
    return hashlib.sha256(value.tobytes()).hexdigest()


def file_record(path: Path, guard: Guard) -> dict[str, Any]:
    return {
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path, guard),
    }


def directory_manifest(directory: Path, guard: Guard) -> dict[str, Any]:
    return {
        "schema_version": 2,
        "no_overwrite": True,
        "files": {
            path.name: file_record(path, guard)
            for path in sorted(directory.iterdir())
            if path.is_file() and path.name != MANIFEST_NAME
        },
    }


def validate_directory_manifest(directory: Path, guard: Guard) -> dict[str, Any]:
    manifest_path = directory / MANIFEST_NAME
    require(manifest_path.is_file(), f"missing artifact manifest: {directory}")
    manifest = read_json(manifest_path)
    require(manifest.get("schema_version") == 2, "v2 manifest schema changed")
    require(manifest.get("no_overwrite") is True, "manifest overwrite policy changed")
    expected_names = set(manifest["files"])
    actual_names = {
        path.name
        for path in directory.iterdir()
        if path.is_file() and path.name != MANIFEST_NAME
    }
    require(actual_names == expected_names, f"artifact set changed: {directory}")
    for name, expected in manifest["files"].items():
        path = directory / name
        require(path.stat().st_size == expected["bytes"], f"bytes changed: {path}")
        require(
            sha256_file(path, guard) == expected["sha256"],
            f"hash changed: {path}",
        )
    return manifest


def atomic_json(path: Path, payload: Mapping[str, Any], guard: Guard) -> None:
    require(not path.exists(), f"refuse overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    guard.check(f"before JSON serialization {path.name}")
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        guard.check(f"after JSON serialization {path.name}")
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def atomic_npz(path: Path, guard: Guard, **arrays: np.ndarray) -> None:
    require(not path.exists(), f"refuse overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    guard.check(f"before NPZ serialization {path.name}")
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".npz",
    )
    os.close(descriptor)
    try:
        np.savez_compressed(temporary, **arrays)
        os.replace(temporary, path)
        guard.check(f"after NPZ serialization {path.name}")
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


@contextmanager
def atomic_stage_directory(path: Path, guard: Guard) -> Iterator[Path]:
    require(not path.exists(), f"refuse overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{path.name}.", dir=path.parent))
    try:
        yield staging
        require(not path.exists(), f"refuse overwrite: {path}")
        guard.check(f"before finalizing {path.name}")
        os.replace(staging, path)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


@contextmanager
def exclusive_lock(path: Path) -> Iterator[None]:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor: int | None = None
    owner: tuple[int, int] | None = None
    try:
        descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        stat = os.fstat(descriptor)
        owner = (stat.st_dev, stat.st_ino)
        os.write(descriptor, f"{os.getpid()}\n".encode())
        os.fsync(descriptor)
        os.close(descriptor)
        descriptor = None
        yield
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if owner is not None:
            try:
                current = path.stat()
            except FileNotFoundError:
                pass
            else:
                if (current.st_dev, current.st_ino) == owner:
                    path.unlink()


def protocol() -> dict[str, Any]:
    value = read_json(PROTOCOL_PATH)
    require(value.get("schema_version") == 2, "protocol schema changed")
    require(value.get("experiment") == HERE.name, "protocol experiment changed")
    require(value.get("status") == "SOURCE_ONLY_V2_UNAUTHORIZED", "status changed")
    require(
        value["candidate_family"]["candidates"] == list(CANDIDATES),
        "candidate family drifted",
    )
    rank = value["numerical_rank"]
    require(
        rank["epsilon_multiplier"] == NUMERICAL_RANK_EPS_MULTIPLIER,
        "rank tolerance multiplier changed",
    )
    memory = value["memory"]
    require(
        memory["max_process_tree_rss_gib"] == MAX_RSS_CEILING_GIB,
        "RSS ceiling changed",
    )
    require(
        memory["min_free_memory_percent"] == MIN_FREE_FLOOR_PERCENT,
        "free-memory floor changed",
    )
    return value


def resolve_repo_path(value: str) -> Path:
    path = (ROOT / value).resolve()
    require(ROOT == path or ROOT in path.parents, f"path escapes repository: {value}")
    return path


def candidate_configs() -> tuple[dict[str, float | str], ...]:
    return tuple(dict(item) for item in CANDIDATES)


def source_paths() -> tuple[Path, ...]:
    paths = [
        RUN_PATH,
        README_PATH,
        TEST_PATH,
        PROTOCOL_PATH,
        MEMORY_GUARD_SCRIPT,
        ALL_TARGETS,
        RESPONSE_GENES,
        PANEL_TARGETS,
        ESSENTIAL_TARGETS,
    ]
    for root in DONOR_ROOTS.values():
        paths.extend(root / name for name in ("delta.npy", "perts.txt", "genes.txt"))
    return tuple(paths)


def ensure_source_allowlist() -> None:
    for path in source_paths():
        resolved = path.resolve()
        require(
            ROOT == resolved or ROOT in resolved.parents,
            f"source path escapes repository: {path}",
        )
        normalized = str(resolved).lower()
        for forbidden in FORBIDDEN_SOURCE_PATH_PARTS:
            require(
                forbidden.lower() not in normalized,
                f"forbidden source path: {path}",
            )


def authorization_bindings(authorization: Mapping[str, Any]) -> dict[Path, str]:
    raw = authorization.get("bindings")
    require(isinstance(raw, Mapping), "authorization lacks bindings")
    result = {}
    for name, value in raw.items():
        path = Path(str(name))
        resolved = (
            path.resolve() if path.is_absolute() else resolve_repo_path(str(name))
        )
        digest = value if isinstance(value, str) else value.get("sha256")
        require(isinstance(digest, str), f"invalid authorization binding: {name}")
        require(resolved not in result, f"duplicate authorization binding: {name}")
        result[resolved] = digest
    return result


def validate_bound_paths(
    authorization: Mapping[str, Any],
    required: Sequence[Path],
    guard: Guard,
) -> None:
    bindings = authorization_bindings(authorization)
    require(
        set(bindings) == {path.resolve() for path in required}, "binding set changed"
    )
    for path in required:
        resolved = path.resolve()
        require(resolved.is_file(), f"authorized source missing: {path}")
        require(
            sha256_file(resolved, guard) == bindings[resolved],
            f"authorization hash mismatch: {path}",
        )


def require_safety_flags(authorization: Mapping[str, Any]) -> None:
    required_false = (
        "development_outcomes_authorized",
        "held_outcomes_authorized",
        "sealed_613_410_values_authorized",
        "kolf_values_authorized",
        "reserved_expression_authorized",
        "challenge_controls_authorized",
        "cell_eval2_authorized",
        "counts_authorized",
        "build_authorized",
        "submission_authorized",
        "commit_authorized",
        "push_authorized",
    )
    for key in required_false:
        require(authorization.get(key) is False, f"unsafe authorization flag: {key}")


def validate_prepare_authorization_document(
    authorization: Mapping[str, Any],
) -> None:
    require(authorization.get("schema_version") == 2, "authorization schema changed")
    require(
        authorization.get("experiment") == HERE.name,
        "authorization experiment changed",
    )
    require(
        authorization.get("stage") == "source_only_prepare_v2",
        "authorization stage changed",
    )
    require(
        authorization.get("confirmation") == AUTHORIZATION_CONFIRMATION,
        "authorization confirmation changed",
    )
    require(
        authorization.get("source_only_prepare_v2_authorized") is True,
        "source-only v2 preparation is not authorized",
    )
    require(
        authorization.get("independent_no_truth_audit_passed") is True,
        "independent no-truth audit has not passed",
    )
    require_safety_flags(authorization)


def validate_prepare_authorization(
    requested_max_rss_gib: float,
    requested_min_free_percent: float,
) -> tuple[dict[str, Any], float, float]:
    require(
        PREPARE_AUTHORIZATION.is_file(),
        "SOURCE_ONLY_PREPARE_AUTHORIZATION.v2.json is absent; "
        "independent no-truth audit is required before preparation",
    )
    authorization = read_json(PREPARE_AUTHORIZATION)
    validate_prepare_authorization_document(authorization)
    authorized_max = float(authorization["max_process_tree_rss_gib"])
    authorized_floor = float(authorization["min_free_memory_percent"])
    effective_max = min(
        float(requested_max_rss_gib),
        authorized_max,
        MAX_RSS_CEILING_GIB,
    )
    effective_floor = max(
        float(requested_min_free_percent),
        authorized_floor,
        MIN_FREE_FLOOR_PERCENT,
    )
    require(effective_max > 0.0, "RSS limit cannot be bypassed")
    require(effective_floor < 100.0, "free-memory floor is invalid")
    return authorization, effective_max, effective_floor


def validate_external_guard(
    effective_max_rss_gib: float,
    effective_min_free_percent: float,
) -> dict[str, float | str]:
    require(
        os.environ.get(EXTERNAL_GUARD_ACTIVE_ENV) == "1",
        "prepare must run under the Experiment 144 external memory guard",
    )
    guard_path_raw = os.environ.get(EXTERNAL_GUARD_PATH_ENV)
    require(guard_path_raw is not None, "external guard path is absent")
    guard_path = Path(guard_path_raw).resolve()
    require(
        guard_path == MEMORY_GUARD_SCRIPT.resolve(),
        "unexpected external memory guard",
    )
    try:
        external_max = float(os.environ[EXTERNAL_GUARD_MAX_RSS_ENV])
        external_floor = float(os.environ[EXTERNAL_GUARD_MIN_FREE_ENV])
    except (KeyError, ValueError) as error:
        raise ContractError("external memory-guard limits are invalid") from error
    require(
        0.0 < external_max <= effective_max_rss_gib,
        "external RSS ceiling is looser than the authorized effective ceiling",
    )
    require(
        effective_min_free_percent <= external_floor < 100.0,
        "external free-memory floor is looser than the authorized effective floor",
    )
    try:
        request_fd = int(os.environ[EXTERNAL_GUARD_REQUEST_FD_ENV])
        response_fd = int(os.environ[EXTERNAL_GUARD_RESPONSE_FD_ENV])
        token = os.environ[EXTERNAL_GUARD_TOKEN_ENV]
    except (KeyError, ValueError) as error:
        raise ContractError("external memory-guard handshake is absent") from error
    require(
        len(token) == 64
        and all(character in "0123456789abcdef" for character in token),
        "external memory-guard handshake token is invalid",
    )
    try:
        os.fstat(request_fd)
        os.fstat(response_fd)
        request = f"{token}\n".encode()
        require(
            os.write(request_fd, request) == len(request),
            "external memory-guard handshake request was incomplete",
        )
        readable, _, _ = select.select(
            [response_fd],
            [],
            [],
            EXTERNAL_GUARD_HANDSHAKE_TIMEOUT_SECONDS,
        )
        require(bool(readable), "external memory-guard handshake timed out")
        response = os.read(response_fd, len(request))
        require(
            response == request,
            "external memory-guard handshake response was invalid",
        )
    except OSError as error:
        raise ContractError("external memory-guard handshake failed") from error
    finally:
        for descriptor in (request_fd, response_fd):
            try:
                os.close(descriptor)
            except OSError:
                pass
    return {
        "path": str(guard_path),
        "max_rss_gib": external_max,
        "min_free_percent": external_floor,
        "live_handshake_verified": "request-response-pipe",
    }


def verify_v1_quarantine(guard: Guard) -> dict[str, Any]:
    frozen = protocol()["v1_quarantine"]["sha256"]
    expected_v1_members = {
        relative for relative in frozen if relative.startswith("work/preoutcome-v1/")
    }
    actual_v1_members = {
        str(path.relative_to(HERE)) for path in V1_OUTPUT_DIR.rglob("*")
    }
    require(
        actual_v1_members == expected_v1_members,
        "quarantined v1 directory membership changed",
    )
    verified = {}
    for relative, expected in frozen.items():
        path = HERE / relative
        require(path.is_file(), f"quarantined v1 file missing: {path}")
        actual = sha256_file(path, guard)
        require(actual == expected, f"quarantined v1 file changed: {path}")
        verified[relative] = actual
    return {
        "files_verified": len(verified),
        "directory_membership_verified": True,
        "v1_authorization": str(V1_AUTHORIZATION.relative_to(HERE)),
        "v1_output": str(V1_OUTPUT_DIR.relative_to(HERE)),
    }


def forbidden_surface_audit() -> dict[str, Any]:
    tree = ast.parse(RUN_PATH.read_text(), filename=str(RUN_PATH))
    imported_roots = set()
    function_names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_roots.update(alias.name.split(".", 1)[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported_roots.add(node.module.split(".", 1)[0])
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            function_names.add(node.name)
    bad_imports = sorted(imported_roots & FORBIDDEN_IMPORT_ROOTS)
    bad_functions = sorted(function_names & FORBIDDEN_FUNCTION_NAMES)
    require(not bad_imports, f"forbidden imports: {bad_imports}")
    require(not bad_functions, f"forbidden functions: {bad_functions}")
    require(COMMANDS == ("prepare", "audit"), "command surface widened")
    ensure_source_allowlist()
    return {
        "commands": list(COMMANDS),
        "forbidden_imports_found": bad_imports,
        "forbidden_functions_found": bad_functions,
        "source_path_count": len(source_paths()),
    }


def last_columns(source: Sequence[str], wanted: Sequence[str]) -> np.ndarray:
    index = {gene: position for position, gene in enumerate(source)}
    require(all(gene in index for gene in wanted), "source misses a response gene")
    return np.asarray([index[gene] for gene in wanted], dtype=np.int64)


def allowed_target_names(
    donors: Sequence[str],
    donor_roots: Mapping[str, Path],
    forbidden: set[str],
) -> tuple[str, ...]:
    target_sets = [set(read_lines(donor_roots[name] / "perts.txt")) for name in donors]
    common = set.intersection(*target_sets)
    allowed = tuple(sorted(common - forbidden))
    require(bool(allowed), "no allowed matched targets")
    require(not (set(allowed) & forbidden), "forbidden target entered source basis")
    return allowed


@dataclass
class DonorStatistics:
    delta: np.ndarray
    rows: np.ndarray
    columns: np.ndarray
    mean: np.ndarray
    raw_median_norm: float


def fit_allowed_donor_statistics(
    delta: np.ndarray,
    rows: np.ndarray,
    columns: np.ndarray,
    guard: Guard,
    *,
    label: str,
    chunk: int = 256,
) -> tuple[np.ndarray, float]:
    """Fit center and norm scale using only the explicitly selected rows."""
    selected = np.asarray(rows, dtype=np.int64)
    response_columns = np.asarray(columns, dtype=np.int64)
    require(selected.ndim == 1 and len(selected) > 0, f"{label}: empty rows")
    require(
        response_columns.ndim == 1 and len(response_columns) > 0,
        f"{label}: empty columns",
    )
    require(
        int(selected.min()) >= 0 and int(selected.max()) < delta.shape[0],
        f"{label}: row index out of bounds",
    )
    require(
        int(response_columns.min()) >= 0
        and int(response_columns.max()) < delta.shape[1],
        f"{label}: column index out of bounds",
    )
    total = np.zeros(len(response_columns), dtype=np.float64)
    for lo in range(0, len(selected), chunk):
        chosen = selected[lo : lo + chunk]
        values = np.asarray(
            delta[np.ix_(chosen, response_columns)],
            dtype=np.float64,
        )
        total += values.sum(axis=0)
        guard.check(f"{label} allowed-row center")
    mean = total / len(selected)

    norms = np.empty(len(selected), dtype=np.float64)
    for lo in range(0, len(selected), chunk):
        hi = min(lo + chunk, len(selected))
        chosen = selected[lo:hi]
        values = np.asarray(
            delta[np.ix_(chosen, response_columns)],
            dtype=np.float64,
        )
        values -= mean
        norms[lo:hi] = np.linalg.norm(values, axis=1)
        guard.check(f"{label} allowed-row normalization")
    raw_median_norm = float(np.median(norms))
    require(
        np.isfinite(raw_median_norm) and raw_median_norm > 0.0,
        f"{label}: invalid allowed-row median norm",
    )
    return mean, raw_median_norm


def build_consensus_basis(
    donors: Sequence[str],
    donor_roots: Mapping[str, Path],
    responses: Sequence[str],
    forbidden: set[str],
    guard: Guard,
    *,
    chunk: int = 128,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Build a matched-target consensus basis with strictly local statistics."""
    targets = allowed_target_names(donors, donor_roots, forbidden)
    prepared: dict[str, DonorStatistics] = {}
    donor_reports = {}

    for donor in donors:
        root = donor_roots[donor]
        perts = read_lines(root / "perts.txt")
        genes = read_lines(root / "genes.txt", unique=False)
        row_index = {target: index for index, target in enumerate(perts)}
        rows = np.asarray([row_index[target] for target in targets], dtype=np.int64)
        columns = last_columns(genes, responses)
        delta = np.load(root / "delta.npy", mmap_mode="r", allow_pickle=False)
        require(delta.ndim == 2, f"{donor}: delta must be a matrix")
        require(
            delta.shape == (len(perts), len(genes)),
            f"{donor}: delta shape changed",
        )
        guard.check(f"{donor} source mmap")
        mean, raw_median_norm = fit_allowed_donor_statistics(
            delta,
            rows,
            columns,
            guard,
            label=donor,
            chunk=chunk,
        )
        prepared[donor] = DonorStatistics(
            delta=delta,
            rows=rows,
            columns=columns,
            mean=mean,
            raw_median_norm=raw_median_norm,
        )
        donor_reports[donor] = {
            "selected_rows": len(rows),
            "selected_row_index_sha256": array_sha256(rows),
            "allowed_row_center_sha256": array_sha256(mean),
            "allowed_row_raw_median_norm": raw_median_norm,
        }

    scales = {donor: prepared[donor].raw_median_norm for donor in donors}
    reference_scale = float(np.median(tuple(scales.values())))
    require(
        np.isfinite(reference_scale) and reference_scale > 0.0,
        "invalid role-local reference norm",
    )
    renormalization = {donor: reference_scale / scales[donor] for donor in donors}

    basis = np.empty((len(targets), len(responses)), dtype=np.float64)
    guard.check("consensus basis allocation")
    for lo in range(0, len(targets), chunk):
        hi = min(lo + chunk, len(targets))
        consensus = np.zeros((hi - lo, len(responses)), dtype=np.float64)
        for donor in donors:
            item = prepared[donor]
            chosen = item.rows[lo:hi]
            values = np.asarray(
                item.delta[np.ix_(chosen, item.columns)],
                dtype=np.float64,
            )
            values -= item.mean
            consensus += renormalization[donor] * values
        basis[lo:hi] = consensus / len(donors)
        guard.check("matched-target consensus rows")

    consensus_mean = np.asarray(basis.mean(axis=0), dtype=np.float64)
    basis -= consensus_mean
    residual_mean = np.asarray(basis.mean(axis=0), dtype=np.float64)
    guard.check("allowed-row consensus centering")

    for donor in donors:
        donor_reports[donor]["allowed_row_renormalization"] = float(
            renormalization[donor]
        )
    report = {
        "targets": targets,
        "target_axis_sha256": hashlib.sha256("\n".join(targets).encode()).hexdigest(),
        "basis_rows": len(targets),
        "basis_columns": len(responses),
        "donors": list(donors),
        "donor_statistics": donor_reports,
        "role_reference_median_norm": reference_scale,
        "consensus_mean_sha256": array_sha256(consensus_mean),
        "residual_mean_max_abs": float(np.max(np.abs(residual_mean))),
        "statistics_scope": "explicitly_allowed_matched_rows_only",
    }
    return basis, report


def regularize_spectrum(
    eigenvalues: np.ndarray,
    *,
    n_rows: int,
    n_columns: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    raw = np.asarray(eigenvalues, dtype=np.float64)
    require(raw.shape == (n_columns,), "eigenvalue count changed")
    require(np.all(np.isfinite(raw)), "non-finite covariance spectrum")
    maximum = float(raw.max())
    require(maximum > 0.0, "non-positive covariance spectrum")
    tolerance = (
        maximum
        * max(n_rows, n_columns)
        * np.finfo(np.float64).eps
        * NUMERICAL_RANK_EPS_MULTIPLIER
    )
    require(
        float(raw.min()) >= -tolerance,
        "covariance has a materially negative eigenvalue",
    )
    numerical_rank = int(np.count_nonzero(raw > tolerance))
    require(
        numerical_rank == n_columns,
        f"covariance numerical rank {numerical_rank}/{n_columns} "
        f"at tolerance {tolerance:.6e}",
    )
    floor = max(maximum * WHITEN_FLOOR_RELATIVE, tolerance)
    regularized = np.maximum(raw, floor)
    require(
        np.all(np.isfinite(regularized)) and np.all(regularized > 0.0),
        "regularized spectrum is not strictly positive",
    )
    return regularized, {
        "numerical_rank": numerical_rank,
        "expected_rank": n_columns,
        "rank_tolerance": tolerance,
        "rank_tolerance_formula": (
            "lambda_max * max(n_rows,n_columns) * eps(float64) * 10"
        ),
        "raw_eigenvalue_min": float(raw.min()),
        "raw_eigenvalue_max": maximum,
        "regularization_floor": floor,
        "regularized_eigenvalue_min": float(regularized.min()),
        "regularized_positive_spectrum": bool(np.all(regularized > 0.0)),
    }


def fit_consensus_operator(
    group: str,
    donors: Sequence[str],
    donor_roots: Mapping[str, Path],
    responses: Sequence[str],
    forbidden: set[str],
    output: Path,
    guard: Guard,
) -> dict[str, Any]:
    require(not output.exists(), f"operator exists: {output}")
    basis, basis_report = build_consensus_basis(
        donors,
        donor_roots,
        responses,
        forbidden,
        guard,
    )
    guard.check(f"before {group} covariance")
    covariance = np.asarray((basis.T @ basis) / len(basis), dtype=np.float64)
    guard.check(f"after {group} covariance")
    guard.check(f"before {group} eigensolver")
    raw_evals, evecs = scipy.linalg.eigh(
        covariance,
        overwrite_a=True,
        check_finite=False,
        driver="evr",
    )
    del covariance
    guard.check(f"after {group} eigensolver")
    evals, spectrum_report = regularize_spectrum(
        raw_evals,
        n_rows=len(basis),
        n_columns=basis.shape[1],
    )
    evecs = np.asarray(evecs, dtype=np.float32)
    atomic_npz(
        output,
        guard,
        evals=evals,
        raw_evals=np.asarray(raw_evals, dtype=np.float64),
        evecs=evecs,
    )
    report = {
        "schema_version": 2,
        "group": group,
        "candidate_gamma": WHITEN_GAMMA,
        "candidate_eta": BLEND_ETA,
        "basis": {
            key: value for key, value in basis_report.items() if key != "targets"
        },
        "spectrum": spectrum_report,
        "operator": file_record(output, guard),
        "excluded_rows_never_loaded_into_fitted_statistics": True,
    }
    atomic_json(output.with_suffix(".basis.json"), report, guard)
    del basis, raw_evals, evals, evecs
    guard.check(f"{group} operator complete")
    return report


def apply_operator(
    matrix: np.ndarray,
    evals: np.ndarray,
    evecs: np.ndarray,
    gamma: float,
    *,
    chunk: int = 64,
) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64)
    weights = np.asarray(evals, dtype=np.float64) ** (-gamma / 2.0)
    output = np.empty(values.shape, dtype=np.float32)
    for lo in range(0, len(values), chunk):
        hi = min(lo + chunk, len(values))
        spectral = (values[lo:hi] @ evecs) * weights[None, :]
        output[lo:hi] = (spectral @ evecs.T).astype(np.float32)
    return output


def row_norm_preserving_blend(
    baseline: np.ndarray,
    candidate: np.ndarray,
    eta: float,
) -> np.ndarray:
    base = np.asarray(baseline, dtype=np.float64)
    other = np.asarray(candidate, dtype=np.float64)
    require(base.shape == other.shape, "blend shape mismatch")
    if eta == 0.0:
        return np.asarray(baseline).copy()
    base_norm = np.linalg.norm(base, axis=1, keepdims=True)
    other_norm = np.linalg.norm(other, axis=1, keepdims=True)
    base_unit = np.divide(
        base,
        np.maximum(base_norm, EPS),
        out=np.zeros_like(base),
    )
    other_unit = np.divide(
        other,
        np.maximum(other_norm, EPS),
        out=np.zeros_like(other),
    )
    mixed = (1.0 - eta) * base_unit + eta * other_unit
    mixed_norm = np.linalg.norm(mixed, axis=1, keepdims=True)
    output = base_norm * np.divide(
        mixed,
        np.maximum(mixed_norm, EPS),
        out=np.zeros_like(mixed),
    )
    return output.astype(np.float32)


def candidate_spec() -> dict[str, Any]:
    return {
        "schema_version": 2,
        "candidates": list(CANDIDATES),
        "candidate_count": len(CANDIDATES),
        "operator_application": "q_consensus = W_consensus(gamma=0.65) @ q_pre",
        "blend": (
            "normalize((1-eta)*normalize(q_identity)"
            "+eta*normalize(q_consensus))*norm(q_identity)"
        ),
        "identity_is_exact_fallback": True,
        "development_evaluation_authorized": False,
    }


def prepare(max_rss_gib: float, min_free_percent: float) -> None:
    protocol()
    forbidden_surface_audit()
    authorization, effective_max, effective_floor = validate_prepare_authorization(
        max_rss_gib,
        min_free_percent,
    )
    external_guard = validate_external_guard(effective_max, effective_floor)
    guard = MemoryGuard(effective_max, effective_floor)
    guard.check("source-only v2 start")
    validate_bound_paths(authorization, source_paths(), guard)
    verify_v1_quarantine(guard)
    require(not OUTPUT_DIR.exists(), f"v2 output exists: {OUTPUT_DIR}")

    responses = read_lines(RESPONSE_GENES)
    require(len(responses) == RESPONSE_COUNT, "response axis changed")
    frozen_targets = set(read_lines(ALL_TARGETS))
    panel_targets = set(read_csv_first_column(PANEL_TARGETS))
    essential_targets = set(read_lines(ESSENTIAL_TARGETS))
    forbidden = frozen_targets | panel_targets | essential_targets
    guard.check("exclusion axes loaded")

    with (
        exclusive_lock(PREPARE_LOCK),
        atomic_stage_directory(OUTPUT_DIR, guard) as staging,
    ):
        atomic_json(staging / CANDIDATE_SPEC_NAME, candidate_spec(), guard)
        group_reports = {}
        for group, donors in GROUP_DONORS.items():
            operator_path = staging / f"consensus-operator-{group}.npz"
            group_reports[group] = fit_consensus_operator(
                group,
                donors,
                DONOR_ROOTS,
                responses,
                forbidden,
                operator_path,
                guard,
            )
        result = {
            "schema_version": 2,
            "stage": "source_only_consensus_covariance_v2",
            "status": "PREPARED_PENDING_INDEPENDENT_NO_TRUTH_AUDIT",
            "candidate_names": [item["name"] for item in CANDIDATES],
            "candidate_count": len(CANDIDATES),
            "groups": group_reports,
            "exclusions": {
                "frozen_613_410_target_names": len(frozen_targets),
                "vcc_panel_target_names": len(panel_targets),
                "essential_target_names": len(essential_targets),
                "union_target_names": len(forbidden),
            },
            "access": {
                "development_or_held_outcomes_opened": False,
                "sealed_613_410_values_opened": False,
                "kolf_values_opened": False,
                "reserved_expression_opened": False,
                "challenge_controls_opened": False,
                "cell_eval2_run": False,
                "counts_generated": False,
                "build_run": False,
                "submission_run": False,
            },
            "independent_no_truth_audit_required_before_evaluation": True,
            "prepare_authorization_sha256": sha256_file(
                PREPARE_AUTHORIZATION,
                guard,
            ),
            "external_memory_guard": external_guard,
            "run_sha256": sha256_file(RUN_PATH, guard),
            "protocol_sha256": sha256_file(PROTOCOL_PATH, guard),
            "memory": guard.report(),
        }
        atomic_json(staging / RESULT_NAME, result, guard)
        atomic_json(
            staging / MANIFEST_NAME,
            directory_manifest(staging, guard),
            guard,
        )
    print(json.dumps(result, indent=2, sort_keys=True))


def audit() -> None:
    protocol_value = protocol()
    guard = MemoryGuard(MAX_RSS_CEILING_GIB, MIN_FREE_FLOOR_PERCENT)
    guard.check("source-only audit start")
    surface = forbidden_surface_audit()
    quarantine = verify_v1_quarantine(guard)
    payload: dict[str, Any] = {
        "schema_version": 2,
        "experiment": HERE.name,
        "status": protocol_value["status"],
        "candidate_names": [item["name"] for item in CANDIDATES],
        "command_surface": surface,
        "v1_quarantine": quarantine,
        "source_only_v2_authorization_exists": PREPARE_AUTHORIZATION.exists(),
        "source_only_v2_output_exists": OUTPUT_DIR.exists(),
        "development_evaluation_command_exists": False,
        "outcome_arrays_opened_by_audit": False,
        "independent_no_truth_audit_still_required": True,
    }
    if OUTPUT_DIR.exists():
        validate_directory_manifest(OUTPUT_DIR, guard)
        result = read_json(OUTPUT_DIR / RESULT_NAME)
        require(
            result["status"] == "PREPARED_PENDING_INDEPENDENT_NO_TRUTH_AUDIT",
            "v2 prepared status changed",
        )
        require(
            all(value is False for value in result["access"].values()),
            "v2 result reports forbidden access",
        )
        payload["source_only_v2_manifest_valid"] = True
    payload["memory"] = guard.report()
    print(json.dumps(payload, indent=2, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=COMMANDS)
    parser.add_argument(
        "--max-rss-gib",
        type=float,
        default=MAX_RSS_CEILING_GIB,
    )
    parser.add_argument(
        "--min-free-percent",
        type=float,
        default=MIN_FREE_FLOOR_PERCENT,
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "prepare":
        prepare(args.max_rss_gib, args.min_free_percent)
    else:
        audit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

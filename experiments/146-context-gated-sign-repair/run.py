#!/usr/bin/env python3
"""Fail-closed source-only runner for Experiment 146.

``prepare`` authenticates and aligns public/source artifacts, excludes the
nine targets overlapping sealed 410 rows, copies only selected donor slices,
and freezes an 81-target source package. ``evaluate`` consumes only that
package and performs strict batch-by-target cross-fitting.

Neither command exists without a separately created stage authorization. This
module has no development, sealed-outcome, count-generation, cell-eval2,
challenge-build, or submission code path.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN_PATH = Path(__file__).resolve()
README_PATH = HERE / "README.md"
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
TEST_PATH = HERE / "test_run.py"

PREPARE_AUTHORIZATION = HERE / "SOURCE_ONLY_PREPARE_AUTHORIZATION.v4.json"
EVALUATE_AUTHORIZATION = HERE / "SOURCE_ONLY_EVALUATE_AUTHORIZATION.v1.json"
PREPARED_DIR = HERE / "work/source-prepared-v1"
EVALUATION_DIR = HERE / "work/source-evaluation-v1"
PREPARE_LOCK = HERE / "work/source-prepared-v1.lock"
EVALUATE_LOCK = HERE / "work/source-evaluation-v1.lock"

EXP111 = ROOT / "experiments/111-depmap-context-program-gain"
DEPMAP_MANIFEST = EXP111 / "SOURCE_MANIFEST.json"
DEPMAP_RECEIPT = EXP111 / "audit/acquisition_receipt.json"
DEPMAP_FREEZE = EXP111 / "SOURCE_MODEL_EXECUTION_FREEZE.v1.json"
DEPMAP_PSEUDOBULK = (
    ROOT / "data/external/depmap26q1/extracted/processed/pseudobulk_sum_by_ko.csv"
)
DEPMAP_FDR = ROOT / "data/external/depmap26q1/extracted/processed/fdr_matrix.csv"
SOURCE_LEDGER = (
    ROOT / "experiments/140-batch-matched-controls/SOURCE_HASH_MANIFEST.v5.json"
)

RESPONSE_AXIS = (
    ROOT / "experiments/132-response-decomposition-ridge/work/"
    "preregistered-freeze-v2/response_genes.txt"
)
GENE_AXIS = ROOT / "data/challenge/gene_names.csv"
PANEL_AXIS = ROOT / "data/challenge/pert_counts.csv"
CHALLENGE_MANIFEST = ROOT / "data/challenge/manifest.json"

DONOR_ROOTS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
MEMORY_GUARD = ROOT / "scripts/run_with_memory_guard.py"

PREPARE_CONFIRMATION = "I_CONFIRM_EXPERIMENT_146_SOURCE_PREPARE_V4"
EVALUATE_CONFIRMATION = "I_CONFIRM_EXPERIMENT_146_SOURCE_EVALUATE_V1"

CONTEXT_COUNT = 16
BATCH_COUNT = 6
ORIGINAL_TARGET_COUNT = 90
TARGET_COUNT = 81
CONTROL_TARGET_COUNT = 10
TARGET_FOLD_COUNT = 5
TARGET_FOLD_SIZES = (17, 16, 16, 16, 16)
CONTROL_PC_RANK = 8
TRAINING_ROW_CAP = 500_000
HGB_MAX_ITER = 100
MODEL_SEED = 146_202_601
BOOTSTRAP_DRAWS = 20_000
BOOTSTRAP_SEED = 146_202_602
TARGET_BOOTSTRAP_SEED = 146_202_604
PERMUTATION_NULL_DRAWS = 9
PERMUTATION_NULL_SEED = 146_202_603
MAX_PROCESS_TREE_RSS_GIB = 18.0
MIN_SYSTEM_FREE_PERCENT = 25
MEMORY_POLL_SECONDS = 0.25
SOURCE_TARGET_SUM = 50_000.0
A15_ALPHA = 0.982
EXPECTED_RESPONSE_AXIS_SHA256 = (
    "24743754eead5b8ddd81634bd12feac1df431eb2b9581197755b75c4c90949c8"
)
EXPECTED_SOURCE_RESPONSE_GENE_COUNT = 4_974
EXPECTED_SOURCE_RESPONSE_AXIS_SHA256 = (
    "fefbe82822b0cdbcc56f49231520e5f52b8bcbd7f400faa5aafdb0602fd8aa79"
)
FDR_AXIS_EXCLUSIONS = (
    "FOXN3",
    "PIGH",
    "CLN6",
    "PIF1",
    "FBXO33",
    "SPTSSA",
)
EXPECTED_AVAILABLE_SOURCE_PAIR_COUNT = 1_289
EXPECTED_FDR_NONFINITE_COUNT = 558_367
EXPECTED_FDR_ALL_MISSING_BLOCK_COUNT = 38
EXPECTED_FDR_ALL_MISSING_BLOCK_SHA256 = (
    "f0d86b37ac94bd29629b49f6d7e8854664df3e485b36c6006ea1928d3377fbeb"
)
SEALED_410_OVERLAP_TARGETS = (
    "COX5B",
    "METAP1",
    "MTPAP",
    "NRBP1",
    "POLR2D",
    "RAB18",
    "SLC39A9",
    "TRAIP",
    "UBE2D3",
)
SINGLE_DONOR_TARGETS = (
    "AURKB",
    "KIF11",
    "MIPEP",
    "MYC",
    "TYMS",
)
EXPECTED_EXPERIMENTAL_90_SHA256 = (
    "6e952a48ca1e8f7f7c5326f41b789d394fe23f311374e83e2b3edf023b4ce1b7"
)
EXPECTED_EXCLUDED_9_SHA256 = (
    "e29b0baa7a0254ce3f9a605b850ecd505fed968f45a517bae052b6a5794b04c3"
)
EXPECTED_EVALUATION_81_SHA256 = (
    "d3824212a27901a526136a42f403082fe5a5a4beca4df0d3410c9dafb78fe979"
)
SCORE_QUANTILE = 0.995
TOTAL_ENERGY_CAP = 0.0003
PDS_AXIS_ENERGY_CAP = 0.0003
SELECTED_COUNT_CAP = 64
SIGNED_L1_MASS_CAP = 0.005
TRAINING_FDR = 0.05
FDR_SENSITIVITY = (0.01, 0.05, 0.10)
WHITENER_GAMMA = 0.65
WHITENER_RIDGE_REL = 1e-3
EPS = 1e-12
GIB = 1024**3
MODEL_ARMS = ("full_hgb", "donor_hgb", "full_logistic", "donor_logistic")
ALL_ARMS = (
    "full_hgb",
    "zero_tail",
    "donor_agreement",
    "full_logistic",
    "donor_logistic",
    "donor_hgb",
)

FEATURE_NAMES = (
    "base_sign",
    "base_log_abs",
    "base_abs_rank",
    "donor_mean",
    "donor_median",
    "donor_mad",
    "donor_median_abs",
    "donor_positive_fraction",
    "donor_support_fraction",
    "donor_agreement_with_base",
    "similarity_weighted_donor_effect",
    "similarity_weighted_donor_sign_vote",
    "nearest_control_donor_effect",
    "destination_control_expression",
    "donor_weighted_control_expression",
    "destination_minus_donor_control_expression",
    "mean_donor_control_similarity",
    "maximum_donor_control_similarity",
    "target_base_norm",
    "target_donor_pairwise_cosine",
    *(f"destination_control_pc_{index}" for index in range(CONTROL_PC_RANK)),
    *(
        f"destination_minus_donor_control_pc_{index}"
        for index in range(CONTROL_PC_RANK)
    ),
)
DONOR_ONLY_FEATURE_INDICES = np.asarray(list(range(13)) + [18, 19], dtype=np.int64)


class ContractError(RuntimeError):
    """Raised when a frozen experiment boundary is violated."""


@dataclass(frozen=True)
class Candidate:
    score_quantile: float = SCORE_QUANTILE
    energy_cap: float = TOTAL_ENERGY_CAP
    pds_axis_energy_cap: float = PDS_AXIS_ENERGY_CAP
    selected_count_cap: int = SELECTED_COUNT_CAP
    signed_l1_mass_cap: float = SIGNED_L1_MASS_CAP

    @property
    def key(self) -> str:
        return "fixed_stage1"


FIXED_CANDIDATE = Candidate()
CANDIDATES = (FIXED_CANDIDATE,)


@dataclass(frozen=True)
class RepairAudit:
    threshold: float
    selected_count: int
    energy_fraction: float
    pds_axis_energy_fraction: float
    signed_l1_mass_fraction: float
    angle_degrees: float


_FIT_ACTIVE = False
_SOURCE_ACCESS_ACTIVE = False


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if sys.platform == "darwin" else value / 1024**2


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    require(
        path.name != "delta.npy",
        f"refuse to hash complete donor delta container: {path}",
    )
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256(
        "".join(f"{value}\n" for value in values).encode("utf-8")
    ).hexdigest()


def array_sha256(value: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def stable_seed(*parts: object, base: int = MODEL_SEED) -> int:
    payload = "|".join([str(base), *(str(part) for part in parts)]).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ContractError(f"cannot read JSON {path}: {error}") from error
    require(isinstance(value, dict), f"JSON root is not an object: {path}")
    return value


def read_lines(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    )
    require(bool(values), f"empty axis: {path}")
    if unique:
        require(len(values) == len(set(values)), f"duplicate axis: {path}")
    return values


def last_occurrence_index(values: Sequence[str]) -> dict[str, int]:
    """Map symbols to their final source-axis coordinate."""
    return {str(value): index for index, value in enumerate(values)}


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames is not None, f"missing CSV header: {path}")
        require(column in reader.fieldnames, f"{path}: missing column {column}")
        values = tuple(row[column] for row in reader)
    require(bool(values), f"empty CSV axis: {path}")
    require(all(values), f"blank value in {path}:{column}")
    require(len(values) == len(set(values)), f"duplicate value in {path}:{column}")
    return values


def read_multilevel_csv_index(
    path: Path,
    *,
    header_rows: int = 2,
) -> tuple[str, ...]:
    """Read only the first-column axis without materializing a wide table."""
    require(header_rows >= 1, "CSV header row count must be positive")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        for _ in range(header_rows):
            require(next(reader, None) is not None, f"truncated CSV header: {path}")
        values = tuple(row[0].strip() for row in reader if row and row[0].strip())
    require(bool(values), f"empty CSV index: {path}")
    require(len(values) == len(set(values)), f"duplicate CSV index: {path}")
    return values


def summarize_fdr_missingness(
    fdr: np.ndarray,
    available: np.ndarray,
    contexts: Sequence[str],
    targets: Sequence[str],
) -> dict[str, Any]:
    values = np.asarray(fdr)
    availability = np.asarray(available, dtype=bool)
    require(values.ndim == 3, "source FDR array must be three-dimensional")
    require(
        values.shape[:2] == availability.shape,
        "source FDR/availability shape mismatch",
    )
    require(
        availability.shape == (len(contexts), len(targets)),
        "source FDR labels do not match context-target axes",
    )
    require(
        np.all(np.isnan(values[~availability])),
        "unavailable source pairs contain FDR values",
    )
    finite = np.isfinite(values)
    finite_counts = np.sum(finite, axis=2)
    missing_count = int(np.sum(~finite[availability]))
    all_missing_indices = np.argwhere(availability & (finite_counts == 0))
    all_missing_pairs = tuple(
        f"{contexts[int(context)]}\t{targets[int(target)]}"
        for context, target in all_missing_indices
    )
    return {
        "available_context_target_pairs": int(np.sum(availability)),
        "nonfinite_entries_on_available_pairs": missing_count,
        "all_missing_available_blocks": len(all_missing_pairs),
        "all_missing_block_sha256": sha256_lines(all_missing_pairs),
        "all_missing_blocks": [
            {
                "context": contexts[int(context)],
                "target": targets[int(target)],
            }
            for context, target in all_missing_indices
        ],
        "policy": "nonfinite_is_unlabelled_never_imputed",
        "unavailable_pairs_all_nan": True,
    }


def validate_fdr_missingness(audit: Mapping[str, Any]) -> None:
    require(
        audit["available_context_target_pairs"] == EXPECTED_AVAILABLE_SOURCE_PAIR_COUNT,
        "available source-pair count changed",
    )
    require(
        audit["nonfinite_entries_on_available_pairs"] == EXPECTED_FDR_NONFINITE_COUNT,
        "source FDR missing-entry count changed",
    )
    require(
        audit["all_missing_available_blocks"] == EXPECTED_FDR_ALL_MISSING_BLOCK_COUNT,
        "all-missing FDR block count changed",
    )
    require(
        audit["all_missing_block_sha256"] == EXPECTED_FDR_ALL_MISSING_BLOCK_SHA256,
        "all-missing FDR block set changed",
    )
    require(
        audit["policy"] == "nonfinite_is_unlabelled_never_imputed",
        "source FDR missingness policy changed",
    )
    require(
        audit["unavailable_pairs_all_nan"] is True,
        "unavailable source-pair FDR values were retained",
    )


def protocol() -> dict[str, Any]:
    value = read_json(PROTOCOL_PATH)
    require(value.get("schema_version") == 1, "protocol schema changed")
    require(value.get("experiment") == HERE.name, "protocol experiment changed")
    require(
        value["source_schema"]["evaluation_targets"] == TARGET_COUNT,
        "evaluation target count changed",
    )
    require(
        tuple(value["source_schema"]["sealed_410_overlap_exclusions"])
        == SEALED_410_OVERLAP_TARGETS,
        "sealed-410 exclusion set changed",
    )
    require(
        value["intervention"]["score_quantile"] == SCORE_QUANTILE,
        "fixed score quantile changed",
    )
    require(
        value["intervention"]["squared_energy_cap"] == TOTAL_ENERGY_CAP,
        "energy cap changed",
    )
    require(
        value["intervention"]["pds_axis_squared_energy_cap"] == PDS_AXIS_ENERGY_CAP,
        "PDS-axis energy cap changed",
    )
    require(
        value["intervention"]["selected_count_cap"] == SELECTED_COUNT_CAP,
        "selected-count cap changed",
    )
    require(
        value["intervention"]["signed_l1_mass_cap"] == SIGNED_L1_MASS_CAP,
        "signed-L1 cap changed",
    )
    require(
        value["target_folds"]["count"] == TARGET_FOLD_COUNT,
        "target fold count changed",
    )
    require(
        tuple(value["target_folds"]["expected_sizes"]) == TARGET_FOLD_SIZES,
        "target fold sizes changed",
    )
    require(
        value["estimator"]["primary"]["max_iter"] == HGB_MAX_ITER,
        "HGB iteration cap changed",
    )
    require(
        value["estimator"]["training_row_cap"] == TRAINING_ROW_CAP,
        "training row cap changed",
    )
    require(
        value["source_schema"]["target_sum"] == SOURCE_TARGET_SUM,
        "source target_sum changed",
    )
    require(
        value["response_axis"]["source_count"] == EXPECTED_SOURCE_RESPONSE_GENE_COUNT,
        "source response-gene count changed",
    )
    require(
        value["response_axis"]["source_sha256"] == EXPECTED_SOURCE_RESPONSE_AXIS_SHA256,
        "source response-gene axis changed",
    )
    require(
        tuple(value["response_axis"]["fdr_axis_exclusions"]) == FDR_AXIS_EXCLUSIONS,
        "FDR-axis exclusion set changed",
    )
    fdr_missingness = value["source_schema"]["fdr_missingness"]
    require(
        fdr_missingness["available_context_target_pairs"]
        == EXPECTED_AVAILABLE_SOURCE_PAIR_COUNT,
        "available source-pair count changed",
    )
    require(
        fdr_missingness["nonfinite_entries_on_available_pairs"]
        == EXPECTED_FDR_NONFINITE_COUNT,
        "source FDR missing-entry count changed",
    )
    require(
        fdr_missingness["all_missing_available_blocks"]
        == EXPECTED_FDR_ALL_MISSING_BLOCK_COUNT,
        "all-missing FDR block count changed",
    )
    require(
        fdr_missingness["all_missing_block_sha256"]
        == EXPECTED_FDR_ALL_MISSING_BLOCK_SHA256,
        "all-missing FDR block set changed",
    )
    require(
        fdr_missingness["policy"] == "nonfinite_is_unlabelled_never_imputed",
        "source FDR missingness policy changed",
    )
    require(
        tuple(value["baseline"]["single_donor_targets"]) == SINGLE_DONOR_TARGETS,
        "single-donor target set changed",
    )
    require(
        value["baseline"]["one_donor_baseline_and_features_allowed"] is True,
        "one-donor baseline/feature policy changed",
    )
    permutation_policy = value["source_metrics"]["whole_target_refit_permutation"]
    require(
        permutation_policy["availability_pattern_stratified"] is True,
        "permutation availability stratification changed",
    )
    require(
        permutation_policy["singleton_class_policy"] == "fixed_and_recorded",
        "permutation singleton policy changed",
    )
    require(
        value["resource_policy"]["max_process_tree_rss_gib"]
        == MAX_PROCESS_TREE_RSS_GIB,
        "RSS ceiling changed",
    )
    require(
        value["resource_policy"]["min_system_free_percent"] == MIN_SYSTEM_FREE_PERCENT,
        "free-memory floor changed",
    )
    require(
        value["metric_masks"]["separate_mask_assertions_required"] is True,
        "separate metric-mask contract disabled",
    )
    return value


def repo_relative(path: Path) -> str:
    resolved = path.resolve()
    require(
        resolved == ROOT or ROOT in resolved.parents,
        f"path escapes repository: {resolved}",
    )
    return str(resolved.relative_to(ROOT))


def assert_write_path(path: Path) -> None:
    resolved = path.resolve()
    require(
        resolved == HERE or HERE in resolved.parents,
        f"write path escapes Experiment 146: {resolved}",
    )


def file_record(path: Path) -> dict[str, Any]:
    return {"bytes": path.stat().st_size, "sha256": sha256_file(path)}


def directory_manifest(directory: Path) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "files": {
            str(path.relative_to(directory)): file_record(path)
            for path in sorted(directory.rglob("*"))
            if path.is_file() and path.name != "artifact-manifest.json"
        },
    }


def validate_directory_manifest(directory: Path) -> dict[str, Any]:
    manifest_path = directory / "artifact-manifest.json"
    require(manifest_path.is_file(), f"missing artifact manifest: {directory}")
    manifest = read_json(manifest_path)
    expected = set(manifest.get("files", {}))
    observed = {
        str(path.relative_to(directory))
        for path in directory.rglob("*")
        if path.is_file() and path.name != "artifact-manifest.json"
    }
    require(observed == expected, f"artifact set changed: {directory}")
    for relative, record in manifest["files"].items():
        path = directory / relative
        require(path.stat().st_size == record["bytes"], f"bytes changed: {path}")
        require(sha256_file(path) == record["sha256"], f"hash changed: {path}")
    return manifest


def atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    assert_write_path(path)
    require(not path.exists(), f"refuse overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
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


def atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    assert_write_path(path)
    require(not path.exists(), f"refuse overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".npz",
    )
    os.close(descriptor)
    try:
        np.savez_compressed(temporary, **arrays)
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


@contextmanager
def atomic_stage_directory(path: Path) -> Iterator[Path]:
    assert_write_path(path)
    require(not path.exists(), f"refuse overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{path.name}.", dir=path.parent))
    try:
        yield staging
        os.replace(staging, path)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


@contextmanager
def exclusive_lock(path: Path) -> Iterator[None]:
    assert_write_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor: int | None = None
    identity: tuple[int, int] | None = None
    try:
        descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        stat = os.fstat(descriptor)
        identity = (stat.st_dev, stat.st_ino)
        os.write(descriptor, f"{os.getpid()}\n".encode())
        os.fsync(descriptor)
        os.close(descriptor)
        descriptor = None
        yield
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if identity is not None:
            try:
                current = path.stat()
            except FileNotFoundError:
                pass
            else:
                if (current.st_dev, current.st_ino) == identity:
                    path.unlink()


def validate_memory_report(report: Mapping[str, Any]) -> None:
    require(report.get("status") == "complete", "memory guard did not complete")
    require(
        float(report.get("max_rss_gib", -1)) == MAX_PROCESS_TREE_RSS_GIB,
        "memory guard RSS ceiling changed",
    )
    require(
        int(report.get("min_free_percent", -1)) == MIN_SYSTEM_FREE_PERCENT,
        "memory guard free-memory floor changed",
    )
    require(
        float(report.get("peak_rss_gib", math.inf)) <= MAX_PROCESS_TREE_RSS_GIB,
        "process-tree RSS exceeded the registered ceiling",
    )
    minimum_free = report.get("minimum_free_percent")
    require(
        minimum_free is not None and int(minimum_free) >= MIN_SYSTEM_FREE_PERCENT,
        "system free-memory floor was not measured or was violated",
    )
    require(int(report.get("samples", 0)) >= 1, "memory guard recorded no samples")


@contextmanager
def process_tree_memory_guard(stage: str) -> Iterator[dict[str, Any]]:
    """Monitor this process and every descendant during one heavy stage."""
    require(MEMORY_GUARD.is_file(), f"missing memory guard: {MEMORY_GUARD}")
    work = HERE / "work"
    assert_write_path(work)
    work.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f".{stage}-memory.", dir=work) as temporary:
        temp = Path(temporary)
        ready = temp / "ready.json"
        stop = temp / "stop"
        report_path = temp / "report.json"
        command = [
            sys.executable,
            str(MEMORY_GUARD),
            "--max-rss-gib",
            str(MAX_PROCESS_TREE_RSS_GIB),
            "--min-free-percent",
            str(MIN_SYSTEM_FREE_PERCENT),
            "--poll-seconds",
            str(MEMORY_POLL_SECONDS),
            "--monitor-pid",
            str(os.getpid()),
            "--ready-file",
            str(ready),
            "--stop-file",
            str(stop),
            "--report-file",
            str(report_path),
            "--quiet",
        ]
        monitor = subprocess.Popen(command)
        deadline = time.monotonic() + 30.0
        while not ready.exists():
            require(monitor.poll() is None, "memory monitor exited before ready")
            require(
                time.monotonic() < deadline,
                "memory monitor readiness timed out",
            )
            time.sleep(0.05)
        report: dict[str, Any] = {}
        try:
            yield report
        finally:
            if monitor.poll() is None:
                stop.touch(exist_ok=False)
            try:
                return_code = monitor.wait(timeout=30.0)
            except subprocess.TimeoutExpired as error:
                monitor.kill()
                monitor.wait(timeout=10.0)
                raise ContractError("memory monitor did not stop") from error
            require(return_code == 0, f"memory monitor failed: {return_code}")
            require(report_path.is_file(), "memory monitor report is missing")
            observed = read_json(report_path)
            validate_memory_report(observed)
            report.update(observed)


@contextmanager
def single_fit_guard(label: str) -> Iterator[None]:
    """Assert that model fits are sequential inside the evaluation process."""
    global _FIT_ACTIVE
    require(not _FIT_ACTIVE, f"concurrent model fit attempted: {label}")
    _FIT_ACTIVE = True
    try:
        yield
    finally:
        _FIT_ACTIVE = False


@contextmanager
def source_array_access() -> Iterator[None]:
    """Permit selected memory-mapped source slices only while the lock is held."""
    global _SOURCE_ACCESS_ACTIVE
    require(not _SOURCE_ACCESS_ACTIVE, "nested source-array access")
    _SOURCE_ACCESS_ACTIVE = True
    try:
        yield
    finally:
        _SOURCE_ACCESS_ACTIVE = False


def prepare_binding_paths() -> tuple[Path, ...]:
    paths = [
        RUN_PATH,
        README_PATH,
        PROTOCOL_PATH,
        TEST_PATH,
        DEPMAP_MANIFEST,
        DEPMAP_RECEIPT,
        DEPMAP_FREEZE,
        SOURCE_LEDGER,
        RESPONSE_AXIS,
        GENE_AXIS,
        PANEL_AXIS,
        CHALLENGE_MANIFEST,
        MEMORY_GUARD,
    ]
    for root in DONOR_ROOTS.values():
        paths.extend(
            root / name
            for name in (
                "genes.txt",
                "perts.txt",
                "meta.json",
            )
        )
    return tuple(paths)


def expected_bindings(paths: Sequence[Path]) -> dict[str, str]:
    for path in paths:
        require(path.is_file() and not path.is_symlink(), f"missing input: {path}")
        require(
            path.name != "delta.npy",
            f"complete donor delta must be ledger-bound, not hashed: {path}",
        )
    return {repo_relative(path): sha256_file(path) for path in paths}


def frozen_array_record(path: Path, ledger: Mapping[str, Any]) -> Mapping[str, Any]:
    files = ledger.get("files")
    require(isinstance(files, Mapping), "source ledger lacks file records")
    key = repo_relative(path)
    record = files.get(key)
    require(isinstance(record, Mapping), f"source ledger lacks {key}")
    require(path.is_file() and not path.is_symlink(), f"missing frozen array: {path}")
    require(
        int(record.get("bytes", -1)) == path.stat().st_size,
        f"frozen array byte count changed: {path}",
    )
    require(
        isinstance(record.get("sha256"), str) and len(record["sha256"]) == 64,
        f"frozen array ledger digest invalid: {path}",
    )
    return record


def selected_matrix_slice(
    path: Path,
    row_indices: Sequence[int],
    column_indices: Sequence[int],
    *,
    ledger: Mapping[str, Any],
) -> np.ndarray:
    """Copy only an explicitly selected two-dimensional memory-mapped slice."""
    require(_SOURCE_ACCESS_ACTIVE, "source array opened outside source-access lock")
    frozen_array_record(path, ledger)
    rows = np.asarray(row_indices, dtype=np.int64)
    columns = np.asarray(column_indices, dtype=np.int64)
    require(rows.ndim == columns.ndim == 1, "selected slice indices must be 1D")
    require(len(rows) <= TARGET_COUNT, "selected slice exceeds 81 target rows")
    matrix = np.load(path, mmap_mode="r", allow_pickle=False)
    require(isinstance(matrix, np.memmap), f"source array is not memory mapped: {path}")
    require(matrix.ndim == 2, f"source delta is not 2D: {path}")
    require(
        np.all((rows >= 0) & (rows < matrix.shape[0])),
        f"selected source rows out of range: {path}",
    )
    require(
        np.all((columns >= 0) & (columns < matrix.shape[1])),
        f"selected source columns out of range: {path}",
    )
    return np.asarray(matrix[np.ix_(rows, columns)], dtype=np.float32).copy()


def selected_vector_slice(
    path: Path,
    indices: Sequence[int],
    *,
    ledger: Mapping[str, Any],
) -> np.ndarray:
    """Copy only explicitly selected coordinates from a frozen vector."""
    require(_SOURCE_ACCESS_ACTIVE, "source vector opened outside source-access lock")
    frozen_array_record(path, ledger)
    columns = np.asarray(indices, dtype=np.int64)
    vector = np.load(path, mmap_mode="r", allow_pickle=False)
    require(
        isinstance(vector, np.memmap), f"source vector is not memory mapped: {path}"
    )
    require(vector.ndim == 1, f"source vector is not 1D: {path}")
    require(
        np.all((columns >= 0) & (columns < vector.shape[0])),
        f"selected source vector indices out of range: {path}",
    )
    return np.asarray(vector[columns], dtype=np.float32).copy()


def validate_authorization_document(
    document: Mapping[str, Any],
    *,
    stage: str,
    confirmation: str,
    bindings: Mapping[str, str],
) -> None:
    require(document.get("schema_version") == 1, "authorization schema changed")
    require(document.get("experiment") == HERE.name, "authorization experiment changed")
    require(document.get("stage") == stage, "authorization stage changed")
    require(document.get("status") == "authorized", "stage is not authorized")
    require(
        document.get("confirmation") == confirmation,
        "authorization confirmation changed",
    )
    require(
        document.get("bindings") == dict(bindings), "authorization bindings drifted"
    )
    permissions = document.get("permissions")
    require(isinstance(permissions, Mapping), "authorization permissions missing")
    allowed = (
        {"source_prepare": True, "source_evaluate": False}
        if stage == "source-only-prepare"
        else {"source_prepare": False, "source_evaluate": True}
    )
    for key, expected in allowed.items():
        require(permissions.get(key) is expected, f"authorization {key} changed")
    forbidden = (
        "development_613",
        "sealed_410",
        "challenge_control_expression",
        "cell_eval2",
        "count_generation",
        "challenge_build",
        "submission",
        "authorization_creation",
        "writes_outside_experiment_146",
    )
    for key in forbidden:
        require(permissions.get(key) is False, f"forbidden permission enabled: {key}")


def validate_stage_authorization(
    stage: str,
    confirmation: str,
    prepared_manifest: Path | None = None,
) -> dict[str, str]:
    if stage == "source-only-prepare":
        path = PREPARE_AUTHORIZATION
    else:
        path = EVALUATE_AUTHORIZATION
        require(prepared_manifest is not None, "evaluation requires prepared manifest")
    require(path.is_file() and not path.is_symlink(), f"missing authorization: {path}")
    if stage == "source-only-prepare":
        bindings = expected_bindings(prepare_binding_paths())
    else:
        bindings = expected_bindings(
            (
                RUN_PATH,
                README_PATH,
                PROTOCOL_PATH,
                TEST_PATH,
                MEMORY_GUARD,
                prepared_manifest,
            )
        )
    document = read_json(path)
    validate_authorization_document(
        document,
        stage=stage,
        confirmation=confirmation,
        bindings=bindings,
    )
    return bindings


def deterministic_target_folds(targets: Sequence[str]) -> np.ndarray:
    require(
        len(targets) == TARGET_COUNT, "target fold construction requires 81 targets"
    )
    require(len(set(targets)) == len(targets), "target axis contains duplicates")
    ordered = sorted(
        range(len(targets)),
        key=lambda index: (
            hashlib.sha256(targets[index].lower().encode()).hexdigest(),
            targets[index],
        ),
    )
    folds = np.empty(len(targets), dtype=np.int8)
    for position, target_index in enumerate(ordered):
        folds[target_index] = position % TARGET_FOLD_COUNT
    sizes = tuple(int(np.sum(folds == fold)) for fold in range(TARGET_FOLD_COUNT))
    require(
        sizes == TARGET_FOLD_SIZES,
        f"unbalanced target folds: {sizes}",
    )
    return folds


def permitted_evaluation_targets(
    experimental_targets: Sequence[str],
) -> tuple[str, ...]:
    """Freeze the 90 -> 81 isolation boundary before any outcome is opened."""
    experimental = tuple(map(str, experimental_targets))
    require(
        len(experimental) == ORIGINAL_TARGET_COUNT,
        "DepMap experimental target count changed",
    )
    require(
        sha256_lines(experimental) == EXPECTED_EXPERIMENTAL_90_SHA256,
        "DepMap experimental target axis changed",
    )
    require(
        sha256_lines(SEALED_410_OVERLAP_TARGETS) == EXPECTED_EXCLUDED_9_SHA256,
        "sealed-410 exclusion tuple changed",
    )
    require(
        set(SEALED_410_OVERLAP_TARGETS).issubset(experimental),
        "a sealed-410 overlap target is absent from the 90-target axis",
    )
    evaluation = tuple(
        target
        for target in experimental
        if target not in set(SEALED_410_OVERLAP_TARGETS)
    )
    require(len(evaluation) == TARGET_COUNT, "evaluation target count is not 81")
    require(
        sha256_lines(evaluation) == EXPECTED_EVALUATION_81_SHA256,
        "frozen 81-target evaluation axis changed",
    )
    require(
        not (set(evaluation) & set(SEALED_410_OVERLAP_TARGETS)),
        "sealed-410 row entered evaluation targets",
    )
    return evaluation


def experimental_batch_folds(batches: Sequence[str]) -> tuple[np.ndarray, ...]:
    require(len(batches) == CONTEXT_COUNT, "batch axis must contain 16 contexts")
    names = tuple(sorted(set(map(str, batches))))
    require(len(names) == BATCH_COUNT, "batch axis must contain six batches")
    folds = tuple(
        np.asarray(
            [index for index, value in enumerate(batches) if value == name],
            dtype=np.int64,
        )
        for name in names
    )
    require(all(len(fold) >= 1 for fold in folds), "empty experimental batch fold")
    require(
        sorted(np.concatenate(folds).tolist()) == list(range(CONTEXT_COUNT)),
        "experimental batch folds do not partition contexts",
    )
    return folds


def complete_gallery_contexts(
    available: np.ndarray,
    context_indices: Sequence[int],
    target_indices: Sequence[int],
) -> np.ndarray:
    """Return contexts with every requested gallery target observed."""
    matrix = np.asarray(available, dtype=bool)
    contexts = np.asarray(context_indices, dtype=np.int64)
    targets = np.asarray(target_indices, dtype=np.int64)
    require(matrix.ndim == 2, "availability matrix must be two-dimensional")
    require(contexts.ndim == targets.ndim == 1, "gallery indices must be 1D")
    require(len(targets) >= 2, "complete gallery requires at least two targets")
    selected = np.asarray(
        [context for context in contexts if np.all(matrix[context, targets])],
        dtype=np.int64,
    )
    require(len(selected) >= 2, "fewer than two complete-gallery contexts")
    return selected


def validate_target_truth_permutation(
    mapping: Mapping[int, int],
    target_indices: Sequence[int],
    available: np.ndarray,
    context_indices: Sequence[int],
) -> None:
    """Require a bijection that preserves every selected availability block."""
    targets = np.asarray(target_indices, dtype=np.int64)
    contexts = np.asarray(context_indices, dtype=np.int64)
    matrix = np.asarray(available, dtype=bool)
    require(matrix.ndim == 2, "availability matrix must be two-dimensional")
    require(targets.ndim == contexts.ndim == 1, "permutation indices must be 1D")
    target_set = set(map(int, targets))
    require(set(mapping) == target_set, "permutation source targets changed")
    require(set(mapping.values()) == target_set, "permutation is not bijective")
    for source, target in mapping.items():
        require(
            np.array_equal(matrix[contexts, source], matrix[contexts, target]),
            "permutation changed a target availability pattern",
        )


def availability_stratified_target_permutations(
    target_indices: Sequence[int],
    available: np.ndarray,
    context_indices: Sequence[int],
    *,
    draws: int = PERMUTATION_NULL_DRAWS,
    seed: int = PERMUTATION_NULL_SEED,
) -> tuple[tuple[dict[int, int], ...], dict[str, Any]]:
    """Permute whole targets only inside identical context-availability classes."""
    targets = np.asarray(target_indices, dtype=np.int64)
    contexts = np.asarray(context_indices, dtype=np.int64)
    matrix = np.asarray(available, dtype=bool)
    require(matrix.ndim == 2, "availability matrix must be two-dimensional")
    require(targets.ndim == contexts.ndim == 1, "permutation indices must be 1D")
    require(len(targets) >= 2, "target permutation requires at least two targets")
    require(len(set(map(int, targets))) == len(targets), "permutation targets repeat")
    require(draws >= 1, "permutation draw count must be positive")
    require(
        np.all((targets >= 0) & (targets < matrix.shape[1])),
        "permutation targets out of range",
    )
    require(
        np.all((contexts >= 0) & (contexts < matrix.shape[0])),
        "permutation contexts out of range",
    )
    classes: dict[tuple[bool, ...], list[int]] = {}
    for target in targets:
        pattern = tuple(bool(value) for value in matrix[contexts, target])
        classes.setdefault(pattern, []).append(int(target))
    ordered_classes = sorted(
        classes.items(),
        key=lambda item: (item[0], tuple(item[1])),
    )
    rng = np.random.default_rng(seed)
    mappings: list[dict[int, int]] = []
    for _ in range(draws):
        mapping: dict[int, int] = {}
        for _pattern, members in ordered_classes:
            local = np.asarray(members, dtype=np.int64)
            if len(local) == 1:
                mapping[int(local[0])] = int(local[0])
                continue
            for _attempt in range(10_000):
                permuted = rng.permutation(local)
                if np.all(permuted != local):
                    break
            else:
                raise ContractError(
                    "could not construct availability-stratified derangement"
                )
            mapping.update(
                {
                    int(target): int(permuted[index])
                    for index, target in enumerate(local)
                }
            )
        validate_target_truth_permutation(
            mapping,
            targets,
            matrix,
            contexts,
        )
        mappings.append(mapping)
    class_records = [
        {
            "availability_pattern": "".join("1" if value else "0" for value in pattern),
            "target_indices": members,
            "size": len(members),
            "permuted": len(members) >= 2,
        }
        for pattern, members in ordered_classes
    ]
    return tuple(mappings), {
        "availability_pattern_stratified": True,
        "context_indices": contexts.astype(int).tolist(),
        "class_count": len(class_records),
        "classes": class_records,
        "singleton_target_indices": [
            record["target_indices"][0]
            for record in class_records
            if record["size"] == 1
        ],
        "singleton_class_policy": "fixed_and_recorded",
        "all_training_blocks_preserved": True,
    }


def metric_masks(
    genes: Sequence[str],
    source_targets: Sequence[str],
    vcc_targets: Sequence[str],
) -> tuple[np.ndarray, np.ndarray]:
    """Return the PDS mask and per-target signed-credit masks."""
    require(len(set(genes)) == len(genes), "gene axis contains duplicates")
    source_set = set(source_targets)
    vcc_set = set(vcc_targets)
    pds = np.asarray(
        [gene not in source_set and gene not in vcc_set for gene in genes],
        dtype=bool,
    )
    signed = np.asarray(
        [
            [gene not in vcc_set and gene != target for gene in genes]
            for target in source_targets
        ],
        dtype=bool,
    )
    require(pds.shape == (len(genes),), "PDS mask shape")
    require(signed.shape == (len(source_targets), len(genes)), "signed mask shape")
    for index, gene in enumerate(genes):
        if gene in vcc_set:
            require(not pds[index], f"VCC target entered PDS: {gene}")
            require(
                not np.any(signed[:, index]), f"VCC target entered signed mask: {gene}"
            )
        if gene in source_set:
            require(not pds[index], f"source-panel target entered PDS: {gene}")
    positions = {gene: index for index, gene in enumerate(genes)}
    for row, target in enumerate(source_targets):
        position = positions.get(target)
        if position is not None:
            require(
                not signed[row, position], f"own target entered signed mask: {target}"
            )
            other_rows = np.arange(len(source_targets)) != row
            require(
                np.all(signed[other_rows, position]),
                f"non-own source target was removed from signed mask: {target}",
            )
    return pds, signed


def authenticated_source_member(
    receipt: Mapping[str, Any],
    member_name: str,
    path: Path,
) -> None:
    member = receipt.get("members", {}).get(member_name)
    require(isinstance(member, Mapping), f"DepMap receipt lacks {member_name}")
    require(path.is_file() and not path.is_symlink(), f"missing source member: {path}")
    require(
        path.stat().st_size == int(member["size_bytes"]),
        f"DepMap source member bytes changed: {member_name}",
    )
    require(
        sha256_file(path) == member["sha256"],
        f"DepMap source member hash changed: {member_name}",
    )


def reconstruct_log1p_profile(
    counts: np.ndarray,
    *,
    target_sum: float = SOURCE_TARGET_SUM,
) -> np.ndarray:
    values = np.asarray(counts, dtype=np.float64)
    require(values.ndim == 1, "pseudobulk profile must be one-dimensional")
    require(np.all(np.isfinite(values)) and np.all(values >= 0), "invalid counts")
    total = float(values.sum())
    require(total > 0, "empty pseudobulk profile")
    return np.log1p(values * (target_sum / total)).astype(np.float32)


def load_depmap_source() -> dict[str, Any]:
    """Open authenticated source outcomes only inside the authorized source lock."""
    require(_SOURCE_ACCESS_ACTIVE, "DepMap outcomes opened outside source-access lock")
    import pandas as pd

    manifest = read_json(DEPMAP_MANIFEST)
    receipt = read_json(DEPMAP_RECEIPT)
    freeze = read_json(DEPMAP_FREEZE)
    freeze_bindings = freeze.get("bindings", {})
    require(
        freeze_bindings.get("source_manifest", {}).get("sha256")
        == sha256_file(DEPMAP_MANIFEST),
        "DepMap source manifest no longer matches its frozen source contract",
    )
    require(
        freeze_bindings.get("acquisition_receipt", {}).get("sha256")
        == sha256_file(DEPMAP_RECEIPT),
        "DepMap acquisition receipt no longer matches its frozen source contract",
    )
    authenticated_source_member(
        receipt,
        "processed/pseudobulk_sum_by_ko.csv",
        DEPMAP_PSEUDOBULK,
    )
    authenticated_source_member(
        receipt,
        "processed/fdr_matrix.csv",
        DEPMAP_FDR,
    )

    contract = manifest["dataset_contract"]
    contexts = tuple(record["cell_line"] for record in contract["cell_lines"])
    batches = tuple(record["batch"] for record in contract["cell_lines"])
    controls = tuple(contract["control_targets"])
    experimental = tuple(
        target for target in contract["targets"] if target not in set(controls)
    )
    targets = permitted_evaluation_targets(experimental)
    require(len(contexts) == CONTEXT_COUNT, "DepMap context count changed")
    require(len(set(batches)) == BATCH_COUNT, "DepMap batch count changed")
    require(len(controls) == CONTROL_TARGET_COUNT, "DepMap control count changed")

    response_parent = read_lines(RESPONSE_AXIS)
    vcc_targets = read_csv_column(PANEL_AXIS, "target_gene")
    require(
        sha256_file(RESPONSE_AXIS) == EXPECTED_RESPONSE_AXIS_SHA256,
        "Experiment 142 response axis hash changed",
    )
    require(len(vcc_targets) == 300, "VCC panel target count changed")
    require(
        sha256_file(PANEL_AXIS) == manifest["challenge_contract"]["panel_sha256"],
        "VCC panel axis no longer matches the authenticated DepMap contract",
    )
    require(
        not (set(experimental) & set(vcc_targets)),
        "DepMap and VCC perturbation targets unexpectedly overlap",
    )
    require(
        not (set(response_parent) & set(vcc_targets)),
        "parent response axis contains VCC targets",
    )

    frame = pd.read_csv(
        DEPMAP_PSEUDOBULK,
        header=[0, 1],
        index_col=0,
        low_memory=False,
    )
    require(not frame.index.has_duplicates, "DepMap response axis contains duplicates")
    require(not frame.columns.has_duplicates, "DepMap context-target columns duplicate")
    all_columns = [(str(a), str(b)) for a, b in frame.columns]
    all_lookup = {pair: index for index, pair in enumerate(all_columns)}
    require(len(all_lookup) == len(all_columns), "DepMap column lookup is nonunique")
    permitted_pairs = [
        (context, target)
        for context in contexts
        for target in (*controls, *targets)
        if (context, target) in all_lookup
    ]
    require(
        all(target not in SEALED_410_OVERLAP_TARGETS for _, target in permitted_pairs),
        "sealed target entered permitted source columns",
    )
    frame = frame.iloc[
        :,
        [all_lookup[pair] for pair in permitted_pairs],
    ]
    complete = ~frame.isna().any(axis=1)
    complete_genes = set(frame.index[complete].astype(str))
    pseudobulk_genes = tuple(gene for gene in response_parent if gene in complete_genes)
    fdr_gene_axis = read_multilevel_csv_index(DEPMAP_FDR)
    fdr_gene_set = set(fdr_gene_axis)
    fdr_exclusions = tuple(
        gene for gene in pseudobulk_genes if gene not in fdr_gene_set
    )
    require(
        fdr_exclusions == FDR_AXIS_EXCLUSIONS,
        f"FDR-axis exclusions changed: {fdr_exclusions}",
    )
    genes = tuple(gene for gene in pseudobulk_genes if gene in fdr_gene_set)
    require(
        len(genes) == EXPECTED_SOURCE_RESPONSE_GENE_COUNT,
        f"source response-gene count changed: {len(genes)}",
    )
    require(
        sha256_lines(genes) == EXPECTED_SOURCE_RESPONSE_AXIS_SHA256,
        "source response-gene axis changed",
    )
    require(not (set(genes) & set(vcc_targets)), "VCC target entered response axis")

    columns = [(str(a), str(b)) for a, b in frame.columns]
    lookup = {pair: index for index, pair in enumerate(columns)}
    require(len(lookup) == len(columns), "DepMap column lookup is nonunique")
    values = frame.loc[list(genes)].to_numpy(dtype=np.float32, copy=True)
    del frame

    control_profiles = np.empty((CONTEXT_COUNT, len(genes)), dtype=np.float32)
    for context_index, context in enumerate(contexts):
        control_columns = [lookup[(context, target)] for target in controls]
        combined = values[:, control_columns].sum(axis=1, dtype=np.float64)
        control_profiles[context_index] = reconstruct_log1p_profile(combined)

    effects = np.zeros((CONTEXT_COUNT, TARGET_COUNT, len(genes)), dtype=np.float32)
    available = np.zeros((CONTEXT_COUNT, TARGET_COUNT), dtype=bool)
    for context_index, context in enumerate(contexts):
        for target_index, target in enumerate(targets):
            column = lookup.get((context, target))
            if column is None:
                continue
            effects[context_index, target_index] = (
                reconstruct_log1p_profile(values[:, column])
                - control_profiles[context_index]
            )
            available[context_index, target_index] = True
    del values

    fdr_frame = pd.read_csv(
        DEPMAP_FDR,
        header=[0, 1],
        index_col=0,
        low_memory=False,
    )
    require(not fdr_frame.index.has_duplicates, "DepMap FDR genes duplicate")
    require(not fdr_frame.columns.has_duplicates, "DepMap FDR columns duplicate")
    require(set(genes).issubset(set(map(str, fdr_frame.index))), "FDR genes incomplete")
    all_fdr_columns = [(str(a), str(b)) for a, b in fdr_frame.columns]
    all_fdr_lookup = {pair: index for index, pair in enumerate(all_fdr_columns)}
    permitted_fdr_pairs = [
        (context, target)
        for context_index, context in enumerate(contexts)
        for target_index, target in enumerate(targets)
        if available[context_index, target_index]
        and (context, target) in all_fdr_lookup
    ]
    require(
        len(permitted_fdr_pairs) == int(np.sum(available)),
        "an available source pair is absent from the FDR table",
    )
    fdr_frame = fdr_frame.iloc[
        :,
        [all_fdr_lookup[pair] for pair in permitted_fdr_pairs],
    ]
    fdr_columns = [(str(a), str(b)) for a, b in fdr_frame.columns]
    fdr_lookup = {pair: index for index, pair in enumerate(fdr_columns)}
    fdr_values = fdr_frame.loc[list(genes)].to_numpy(dtype=np.float32, copy=True)
    del fdr_frame
    fdr = np.full_like(effects, np.nan, dtype=np.float32)
    for context_index, context in enumerate(contexts):
        for target_index, target in enumerate(targets):
            if not available[context_index, target_index]:
                continue
            column = fdr_lookup.get((context, target))
            require(column is not None, "available source pair lacks aligned FDR")
            fdr[context_index, target_index] = fdr_values[:, column]
    del fdr_values

    missing_expected = {
        tuple(pair)
        for pair in contract["pseudobulk"]["missing_context_target_pairs"]
        if pair[1] in set(targets)
    }
    missing_observed = {
        (contexts[context], targets[target])
        for context, target in zip(*np.where(~available), strict=True)
    }
    require(missing_observed == missing_expected, "DepMap missing-pair set changed")
    fdr_missingness = summarize_fdr_missingness(
        fdr,
        available,
        contexts,
        targets,
    )
    validate_fdr_missingness(fdr_missingness)
    folds = deterministic_target_folds(targets)
    pds_mask, signed_masks = metric_masks(genes, targets, vcc_targets)
    return {
        "contexts": contexts,
        "batches": batches,
        "targets": targets,
        "controls": controls,
        "genes": genes,
        "fdr_axis_exclusions": FDR_AXIS_EXCLUSIONS,
        "fdr_missingness": fdr_missingness,
        "vcc_targets": vcc_targets,
        "effects": effects,
        "fdr": fdr,
        "available": available,
        "control_profiles": control_profiles,
        "target_folds": folds,
        "pds_mask": pds_mask,
        "signed_masks": signed_masks,
    }


@dataclass(frozen=True)
class HistoricalWhitener:
    """Low-rank signal-covariance power on a frozen safe coordinate mask."""

    mask: np.ndarray
    components: np.ndarray
    eigenvalues: np.ndarray
    ridge: float

    def apply(self, values: np.ndarray) -> np.ndarray:
        matrix = np.asarray(values, dtype=np.float64)
        output = matrix.copy()
        active = matrix[:, self.mask]
        before = float(np.median(np.linalg.norm(active, axis=1))) or 1.0
        base_weight = self.ridge ** (-WHITENER_GAMMA / 2.0)
        weights = (self.eigenvalues + self.ridge) ** (-WHITENER_GAMMA / 2.0)
        scores = active @ self.components.T
        transformed = active * base_weight
        transformed += (scores * (weights - base_weight)) @ self.components
        after = float(np.median(np.linalg.norm(transformed, axis=1))) or 1.0
        output[:, self.mask] = transformed * (before / after)
        return output.astype(np.float32)


def fit_historical_panel_whitener(
    donor_effects: np.ndarray,
    donor_available: np.ndarray,
    genes: Sequence[str],
    targets: Sequence[str],
    pds_mask: np.ndarray,
) -> tuple[HistoricalWhitener, dict[str, Any]]:
    mask = np.asarray(pds_mask, dtype=bool).copy()
    positions = {gene: index for index, gene in enumerate(genes)}
    for target in targets:
        position = positions.get(target)
        if position is not None:
            require(not mask[position], f"evaluation target entered whitener: {target}")
    require(np.any(mask), "historical whitener has no safe coordinates")
    rows = [
        np.asarray(donor_effects[donor, target, mask], dtype=np.float64)
        for donor in range(donor_effects.shape[0])
        for target in range(donor_effects.shape[1])
        if donor_available[donor, target]
    ]
    require(bool(rows), "historical whitener has no donor rows")
    basis = np.vstack(rows)
    require(np.all(np.isfinite(basis)), "historical whitener basis is nonfinite")
    _, singular, vt = np.linalg.svd(basis, full_matrices=False)
    eigenvalues = singular**2 / len(basis)
    positive = eigenvalues[eigenvalues > EPS]
    require(len(positive) > 0, "historical whitener basis has zero rank")
    ridge = max(float(np.median(positive)) * WHITENER_RIDGE_REL, EPS)
    state = HistoricalWhitener(
        mask=mask,
        components=np.asarray(vt, dtype=np.float64),
        eigenvalues=np.asarray(eigenvalues, dtype=np.float64),
        ridge=ridge,
    )
    return state, {
        "basis_rows": int(len(basis)),
        "safe_coordinates": int(np.sum(mask)),
        "evaluation_target_coordinates_excluded": True,
        "sealed_target_rows_used": False,
        "gamma": WHITENER_GAMMA,
        "ridge_relative": WHITENER_RIDGE_REL,
        "ridge": ridge,
        "rank": int(len(positive)),
    }


def load_donor_slice(
    name: str,
    root: Path,
    targets: Sequence[str],
    response_genes: Sequence[str],
    ledger: Mapping[str, Any],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    require(
        not (set(targets) & set(SEALED_410_OVERLAP_TARGETS)),
        f"{name}: sealed target requested",
    )
    perts = read_lines(root / "perts.txt")
    genes = read_lines(root / "genes.txt", unique=False)
    meta = read_json(root / "meta.json")
    pert_index = last_occurrence_index(perts)
    gene_index = last_occurrence_index(genes)
    require(
        len(genes) == int(meta["n_genes"]),
        f"{name}: gene axis length does not match donor metadata",
    )
    target_positions = [
        index for index, target in enumerate(targets) if target in pert_index
    ]
    source_rows = [pert_index[targets[index]] for index in target_positions]
    response_positions = [
        index for index, gene in enumerate(response_genes) if gene in gene_index
    ]
    source_columns = [gene_index[response_genes[index]] for index in response_positions]
    require(bool(source_rows), f"{name}: no permitted targets available")
    require(bool(source_columns), f"{name}: no response genes available")

    selected = selected_matrix_slice(
        root / "delta.npy",
        source_rows,
        source_columns,
        ledger=ledger,
    )
    effects = np.zeros((TARGET_COUNT, len(response_genes)), dtype=np.float32)
    effects[np.ix_(target_positions, response_positions)] = selected
    available = np.zeros(TARGET_COUNT, dtype=bool)
    available[target_positions] = True
    center = np.mean(effects[available], axis=0, dtype=np.float64)
    effects[available] -= center.astype(np.float32)

    selected_control = selected_vector_slice(
        root / "b_ctrl.npy",
        source_columns,
        ledger=ledger,
    )
    controls = np.zeros(len(response_genes), dtype=np.float32)
    controls[response_positions] = selected_control
    controls = reconstruct_log1p_profile(np.maximum(np.expm1(controls), 0.0))
    audit = {
        "name": name,
        "targets_requested": TARGET_COUNT,
        "targets_available": int(np.sum(available)),
        "response_genes_requested": len(response_genes),
        "response_genes_available": len(response_positions),
        "gene_axis_count": len(genes),
        "gene_axis_unique_count": len(gene_index),
        "duplicate_gene_symbol_count": len(genes) - len(gene_index),
        "duplicate_gene_symbol_policy": "last_occurrence",
        "delta_shape_from_meta": [
            int(meta["n_perturbations"]),
            int(meta["n_genes"]),
        ],
        "delta_ledger": dict(frozen_array_record(root / "delta.npy", ledger)),
        "control_ledger": dict(frozen_array_record(root / "b_ctrl.npy", ledger)),
        "complete_delta_hashed": False,
        "complete_delta_loaded": False,
        "sealed_target_rows_used": False,
    }
    return effects, available, controls, audit


def validate_donor_support(
    donor_available: np.ndarray,
    targets: Sequence[str],
) -> dict[str, Any]:
    availability = np.asarray(donor_available, dtype=bool)
    require(availability.ndim == 2, "donor availability must be two-dimensional")
    require(
        availability.shape[1] == len(targets),
        "donor availability target axis changed",
    )
    support = availability.sum(axis=0)
    require(np.all(support >= 1), "a permitted target has no production donor")
    single = tuple(str(targets[index]) for index in np.flatnonzero(support == 1))
    require(
        set(single) == set(SINGLE_DONOR_TARGETS)
        and len(single) == len(SINGLE_DONOR_TARGETS),
        f"single-donor target set changed: {sorted(single)}",
    )
    return {
        "minimum_support": int(np.min(support)),
        "support_histogram": {
            str(count): int(np.sum(support == count))
            for count in range(1, availability.shape[0] + 1)
        },
        "single_donor_targets": list(SINGLE_DONOR_TARGETS),
        "one_donor_baseline_and_features_allowed": True,
        "intervention_minimum_support": 2,
    }


def build_a15_source_baseline(
    source: Mapping[str, Any],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Build the selected-slice a15 core without opening complete donor arrays."""
    ledger = read_json(SOURCE_LEDGER)
    donor_rows = []
    availability = []
    controls = []
    donor_audits = []
    for name, root in DONOR_ROOTS.items():
        effects, available, control, audit = load_donor_slice(
            name,
            root,
            source["targets"],
            source["genes"],
            ledger,
        )
        donor_rows.append(effects)
        availability.append(available)
        controls.append(control)
        donor_audits.append(audit)
    donor_effects = np.asarray(donor_rows, dtype=np.float32)
    donor_available = np.asarray(availability, dtype=bool)
    donor_controls = np.asarray(controls, dtype=np.float32)
    support_audit = validate_donor_support(donor_available, source["targets"])

    donor_scales = []
    for donor in range(len(donor_effects)):
        norms = np.linalg.norm(
            donor_effects[donor, donor_available[donor]][:, source["pds_mask"]],
            axis=1,
        )
        donor_scales.append(float(np.median(norms)) or 1.0)
    reference_scale = float(np.median(donor_scales))
    for donor, scale in enumerate(donor_scales):
        donor_effects[donor, donor_available[donor]] *= reference_scale / scale

    baseline_pre = np.zeros((TARGET_COUNT, len(source["genes"])), dtype=np.float32)
    for target in range(TARGET_COUNT):
        have = np.flatnonzero(donor_available[:, target])
        baseline_pre[target] = np.mean(donor_effects[have, target], axis=0)
    whitener, whitening_audit = fit_historical_panel_whitener(
        donor_effects,
        donor_available,
        source["genes"],
        source["targets"],
        source["pds_mask"],
    )
    baseline = A15_ALPHA * whitener.apply(baseline_pre)
    require(np.all(np.isfinite(baseline)), "source baseline is nonfinite")
    return (
        baseline,
        donor_effects,
        donor_available,
        donor_controls,
        {
            "baseline": "selected_slice_a15_core",
            "donors": donor_audits,
            "donor_scales": donor_scales,
            "reference_scale": reference_scale,
            "donor_support": support_audit,
            "whitener": whitening_audit,
            "alpha": A15_ALPHA,
            "perturbmatch_used": False,
            "cd4_used": False,
            "complete_donor_delta_hashed": False,
            "complete_donor_delta_loaded": False,
            "sealed_target_rows_used": False,
        },
    )


def prepare_source(confirmation: str) -> None:
    protocol()
    require(confirmation == PREPARE_CONFIRMATION, "prepare confirmation mismatch")
    bindings = validate_stage_authorization("source-only-prepare", confirmation)
    require(not PREPARED_DIR.exists(), f"prepared output exists: {PREPARED_DIR}")
    require(not EVALUATION_DIR.exists(), "evaluation already exists")
    with exclusive_lock(PREPARE_LOCK), atomic_stage_directory(PREPARED_DIR) as stage:
        with process_tree_memory_guard("prepare") as memory_report:
            with source_array_access():
                source = load_depmap_source()
                (
                    baseline,
                    donor_effects,
                    donor_available,
                    donor_controls,
                    baseline_audit,
                ) = build_a15_source_baseline(source)
            require(
                baseline.shape == source["effects"].shape[1:],
                "baseline/source effect shape mismatch",
            )
            arrays_path = stage / "source-data.npz"
            atomic_npz(
                arrays_path,
                effects=np.asarray(source["effects"], dtype=np.float32),
                fdr=np.asarray(source["fdr"], dtype=np.float32),
                available=np.asarray(source["available"], dtype=bool),
                control_profiles=np.asarray(
                    source["control_profiles"], dtype=np.float32
                ),
                baseline=baseline,
                donor_effects=donor_effects,
                donor_available=donor_available,
                donor_controls=donor_controls,
                target_folds=np.asarray(source["target_folds"], dtype=np.int8),
                pds_mask=np.asarray(source["pds_mask"], dtype=bool),
                signed_masks=np.asarray(source["signed_masks"], dtype=bool),
            )
        plan = {
            "schema_version": 1,
            "experiment": HERE.name,
            "stage": "source-prepared",
            "contexts": list(source["contexts"]),
            "batches": list(source["batches"]),
            "targets": list(source["targets"]),
            "control_targets": list(source["controls"]),
            "genes": list(source["genes"]),
            "response_gene_count": len(source["genes"]),
            "response_axis_sha256": sha256_lines(source["genes"]),
            "fdr_axis_exclusions": list(source["fdr_axis_exclusions"]),
            "fdr_missingness": source["fdr_missingness"],
            "vcc_targets": list(source["vcc_targets"]),
            "donors": list(DONOR_ROOTS),
            "target_fold_sizes": [
                int(np.sum(source["target_folds"] == fold))
                for fold in range(TARGET_FOLD_COUNT)
            ],
            "source_target_sum": SOURCE_TARGET_SUM,
            "sealed_410_overlap_exclusions": list(SEALED_410_OVERLAP_TARGETS),
            "metric_masks": {
                "pds_true_count": int(np.sum(source["pds_mask"])),
                "signed_true_counts": np.sum(source["signed_masks"], axis=1)
                .astype(int)
                .tolist(),
                "separate": True,
            },
            "baseline_audit": baseline_audit,
            "input_bindings": bindings,
            "array_hashes": {
                "effects": array_sha256(source["effects"]),
                "fdr": array_sha256(source["fdr"]),
                "baseline": array_sha256(baseline),
                "donor_effects": array_sha256(donor_effects),
                "control_profiles": array_sha256(source["control_profiles"]),
                "pds_mask": array_sha256(source["pds_mask"]),
                "signed_masks": array_sha256(source["signed_masks"]),
            },
            "development_outcomes_opened": False,
            "challenge_control_expression_opened": False,
            "cell_eval2_used": False,
            "complete_donor_delta_hashed": False,
            "complete_donor_delta_loaded": False,
            "memory_report": "memory-report.json",
        }
        atomic_json(stage / "memory-report.json", memory_report)
        atomic_json(stage / "plan.json", plan)
        atomic_json(stage / "artifact-manifest.json", directory_manifest(stage))
    print(f"prepared immutable source package: {PREPARED_DIR}")


def validate_prepared() -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    validate_directory_manifest(PREPARED_DIR)
    plan = read_json(PREPARED_DIR / "plan.json")
    require(plan.get("experiment") == HERE.name, "prepared experiment changed")
    require(plan.get("stage") == "source-prepared", "prepared stage changed")
    require(
        plan.get("development_outcomes_opened") is False, "development boundary changed"
    )
    require(
        plan.get("challenge_control_expression_opened") is False,
        "challenge-control boundary changed",
    )
    with np.load(PREPARED_DIR / "source-data.npz", allow_pickle=False) as data:
        arrays = {name: np.asarray(data[name]) for name in data.files}
    contexts = plan["contexts"]
    targets = plan["targets"]
    genes = plan["genes"]
    require(len(contexts) == CONTEXT_COUNT, "prepared context count changed")
    require(len(targets) == TARGET_COUNT, "prepared target count changed")
    require(
        len(genes) == EXPECTED_SOURCE_RESPONSE_GENE_COUNT,
        "prepared response-gene count changed",
    )
    require(
        plan.get("response_gene_count") == EXPECTED_SOURCE_RESPONSE_GENE_COUNT,
        "prepared response-gene count audit changed",
    )
    require(
        sha256_lines(genes) == EXPECTED_SOURCE_RESPONSE_AXIS_SHA256,
        "prepared response-gene axis changed",
    )
    require(
        plan.get("response_axis_sha256") == EXPECTED_SOURCE_RESPONSE_AXIS_SHA256,
        "prepared response-axis digest changed",
    )
    require(
        tuple(plan.get("fdr_axis_exclusions", ())) == FDR_AXIS_EXCLUSIONS,
        "prepared FDR-axis exclusions changed",
    )
    require(
        arrays["effects"].shape == (CONTEXT_COUNT, TARGET_COUNT, len(genes)),
        "effect shape changed",
    )
    require(
        arrays["available"].shape == (CONTEXT_COUNT, TARGET_COUNT),
        "availability shape changed",
    )
    require(
        arrays["fdr"].shape == arrays["effects"].shape,
        "source FDR shape changed",
    )
    observed_fdr_missingness = summarize_fdr_missingness(
        arrays["fdr"],
        arrays["available"],
        contexts,
        targets,
    )
    validate_fdr_missingness(observed_fdr_missingness)
    require(
        plan.get("fdr_missingness") == observed_fdr_missingness,
        "prepared FDR missingness audit changed",
    )
    require(
        arrays["baseline"].shape == (TARGET_COUNT, len(genes)), "baseline shape changed"
    )
    require(
        arrays["donor_effects"].shape == (3, TARGET_COUNT, len(genes)),
        "donor effect shape changed",
    )
    require(
        arrays["donor_available"].shape == (3, TARGET_COUNT),
        "donor availability shape changed",
    )
    require(
        arrays["control_profiles"].shape == (CONTEXT_COUNT, len(genes)),
        "control shape changed",
    )
    require(
        arrays["donor_controls"].shape == (3, len(genes)), "donor control shape changed"
    )
    expected_pds, expected_signed = metric_masks(genes, targets, plan["vcc_targets"])
    require(np.array_equal(arrays["pds_mask"], expected_pds), "PDS mask changed")
    require(
        np.array_equal(arrays["signed_masks"], expected_signed),
        "signed-credit masks changed",
    )
    expected_folds = deterministic_target_folds(targets)
    require(
        np.array_equal(arrays["target_folds"], expected_folds), "target folds changed"
    )
    require(
        tuple(plan["sealed_410_overlap_exclusions"]) == SEALED_410_OVERLAP_TARGETS,
        "prepared sealed-target exclusion changed",
    )
    require(
        plan.get("complete_donor_delta_hashed") is False
        and plan.get("complete_donor_delta_loaded") is False,
        "prepared donor-array access boundary changed",
    )
    validate_memory_report(read_json(PREPARED_DIR / "memory-report.json"))
    return plan, arrays


def context_weights(indices: Sequence[int], batches: Sequence[str]) -> np.ndarray:
    indices = list(indices)
    represented = sorted({batches[index] for index in indices})
    require(bool(represented), "no represented batches")
    weights = np.zeros(len(indices), dtype=np.float64)
    for batch in represented:
        positions = [
            position
            for position, index in enumerate(indices)
            if batches[index] == batch
        ]
        for position in positions:
            weights[position] = 1.0 / (len(represented) * len(positions))
    require(np.isclose(weights.sum(), 1.0), "context weights do not sum to one")
    return weights


def fix_component_signs(components: np.ndarray) -> np.ndarray:
    fixed = np.asarray(components, dtype=np.float64).copy()
    for column in range(fixed.shape[1]):
        pivot = int(np.argmax(np.abs(fixed[:, column])))
        if fixed[pivot, column] < 0:
            fixed[:, column] *= -1.0
    return fixed


def fit_control_pca(
    controls: np.ndarray,
    batches: Sequence[str],
    train_contexts: Sequence[int],
    rank: int = CONTROL_PC_RANK,
) -> dict[str, np.ndarray]:
    indices = list(train_contexts)
    require(len(indices) > rank, "too few contexts for control PCA")
    weights = context_weights(indices, batches)
    matrix = np.asarray(controls[indices], dtype=np.float64)
    mean = np.sum(matrix * weights[:, None], axis=0)
    centered = matrix - mean
    weighted = centered * np.sqrt(weights[:, None])
    _, _, vt = np.linalg.svd(weighted, full_matrices=False)
    components = fix_component_signs(vt[:rank].T)
    scores = centered @ components
    scale = np.sqrt(np.sum(weights[:, None] * scores**2, axis=0))
    scale = np.where(scale > EPS, scale, 1.0)
    return {"mean": mean, "components": components, "scale": scale}


def transform_controls(
    controls: np.ndarray, state: Mapping[str, np.ndarray]
) -> np.ndarray:
    centered = np.asarray(controls, dtype=np.float64) - state["mean"]
    return (centered @ state["components"]) / state["scale"]


def normalized_control(value: np.ndarray) -> np.ndarray:
    centered = np.asarray(value, dtype=np.float64) - float(np.mean(value))
    norm = float(np.linalg.norm(centered))
    return centered / norm if norm > EPS else np.zeros_like(centered)


def average_pairwise_cosine(matrix: np.ndarray) -> float:
    values = np.asarray(matrix, dtype=np.float64)
    if len(values) < 2:
        return 0.0
    norms = np.linalg.norm(values, axis=1)
    total = 0.0
    count = 0
    for left in range(len(values)):
        for right in range(left + 1, len(values)):
            denominator = float(norms[left] * norms[right])
            total += (
                float(values[left] @ values[right]) / denominator
                if denominator > EPS
                else 0.0
            )
            count += 1
    return total / count


def rank_fraction(values: np.ndarray) -> np.ndarray:
    from scipy.stats import rankdata

    ranks = rankdata(np.asarray(values), method="average")
    denominator = max(len(ranks) - 1, 1)
    return ((ranks - 1.0) / denominator).astype(np.float64)


def build_feature_block(
    *,
    context_control: np.ndarray,
    context_score: np.ndarray,
    donor_controls: np.ndarray,
    donor_scores: np.ndarray,
    donor_effects: np.ndarray,
    donor_available: np.ndarray,
    baseline: np.ndarray,
    target_index: int,
) -> np.ndarray:
    available = np.flatnonzero(donor_available[:, target_index])
    require(len(available) >= 1, "target has no production donor")
    stack = np.asarray(donor_effects[available, target_index], dtype=np.float64)
    base = np.asarray(baseline[target_index], dtype=np.float64)
    context_unit = normalized_control(context_control)
    donor_units = np.vstack(
        [normalized_control(donor_controls[index]) for index in available]
    )
    similarities = donor_units @ context_unit
    logits = 10.0 * (similarities - float(np.max(similarities)))
    weights = np.exp(logits)
    weights /= float(weights.sum())
    weighted_effect = weights @ stack
    weighted_sign = weights @ np.sign(stack)
    nearest = stack[int(np.argmax(similarities))]
    weighted_control = weights @ donor_controls[available]
    weighted_score = weights @ donor_scores[available]
    median = np.median(stack, axis=0)
    mad = np.median(np.abs(stack - median), axis=0)
    base_sign = np.sign(base)
    donor_sign = np.sign(stack)
    agreement = np.mean(donor_sign == base_sign[None, :], axis=0)
    target_norm = math.log1p(float(np.linalg.norm(base)))
    pairwise = average_pairwise_cosine(stack)
    repeated_context_score = np.broadcast_to(
        np.asarray(context_score, dtype=np.float64),
        (len(base), CONTROL_PC_RANK),
    )
    repeated_score_difference = np.broadcast_to(
        np.asarray(context_score - weighted_score, dtype=np.float64),
        (len(base), CONTROL_PC_RANK),
    )
    result = np.column_stack(
        [
            base_sign,
            np.log1p(np.abs(base)),
            rank_fraction(np.abs(base)),
            np.mean(stack, axis=0),
            median,
            mad,
            np.median(np.abs(stack), axis=0),
            np.mean(stack > 0, axis=0),
            np.full(len(base), len(available) / donor_effects.shape[0]),
            agreement,
            weighted_effect,
            weighted_sign,
            nearest,
            context_control,
            weighted_control,
            context_control - weighted_control,
            np.full(len(base), float(np.mean(similarities))),
            np.full(len(base), float(np.max(similarities))),
            np.full(len(base), target_norm),
            np.full(len(base), pairwise),
            repeated_context_score,
            repeated_score_difference,
        ]
    ).astype(np.float32)
    require(result.shape[1] == len(FEATURE_NAMES), "feature width changed")
    require(np.all(np.isfinite(result)), "feature block contains nonfinite values")
    return result


def training_matrix(
    *,
    arrays: Mapping[str, np.ndarray],
    plan: Mapping[str, Any],
    contexts: Sequence[int],
    targets: Sequence[int],
    control_state: Mapping[str, np.ndarray],
    control_mapping: np.ndarray,
    row_cap: int = TRAINING_ROW_CAP,
    seed_label: object,
    target_truth_permutation: Mapping[int, int] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    require(0 < row_cap <= TRAINING_ROW_CAP, "training row cap exceeds 500,000")
    controls = arrays["control_profiles"]
    donor_controls = arrays["donor_controls"]
    mapped_controls = controls[control_mapping]
    context_scores = transform_controls(mapped_controls, control_state)
    donor_scores = transform_controls(donor_controls, control_state)
    blocks = [
        (context, target)
        for context in contexts
        for target in targets
        if arrays["available"][context, target]
    ]
    require(bool(blocks), "training matrix has no blocks")
    if target_truth_permutation is not None:
        validate_target_truth_permutation(
            target_truth_permutation,
            targets,
            arrays["available"],
            contexts,
        )
    quota = max(1, row_cap // len(blocks))
    x_parts: list[np.ndarray] = []
    y_parts: list[np.ndarray] = []
    w_parts: list[np.ndarray] = []
    for context, target in blocks:
        truth_target = (
            target_truth_permutation.get(target, target)
            if target_truth_permutation is not None
            else target
        )
        require(
            arrays["available"][context, truth_target],
            "permutation dropped an available training block",
        )
        features = build_feature_block(
            context_control=mapped_controls[context],
            context_score=context_scores[context],
            donor_controls=donor_controls,
            donor_scores=donor_scores,
            donor_effects=arrays["donor_effects"],
            donor_available=arrays["donor_available"],
            baseline=arrays["baseline"],
            target_index=target,
        )
        truth = arrays["effects"][context, truth_target]
        base = arrays["baseline"][target]
        eligible = (
            arrays["signed_masks"][target]
            & arrays["signed_masks"][truth_target]
            & np.isfinite(arrays["fdr"][context, truth_target])
            & (arrays["fdr"][context, truth_target] <= TRAINING_FDR)
            & np.isfinite(truth)
            & np.isfinite(base)
            & (truth != 0)
            & (base != 0)
        )
        indices = np.flatnonzero(eligible)
        if len(indices) > quota:
            rng = np.random.default_rng(stable_seed(seed_label, context, target))
            indices = np.sort(rng.choice(indices, size=quota, replace=False))
        if not len(indices):
            continue
        labels = (np.sign(truth[indices]) != np.sign(base[indices])).astype(np.uint8)
        weights = np.abs(base[indices]).astype(np.float64)
        total = float(weights.sum())
        if total <= EPS:
            continue
        weights /= total
        x_parts.append(features[indices])
        y_parts.append(labels)
        w_parts.append(weights)
    require(bool(x_parts), "no eligible training rows")
    x = np.vstack(x_parts)
    y = np.concatenate(y_parts)
    weight = np.concatenate(w_parts)
    if len(x) > row_cap:
        rng = np.random.default_rng(stable_seed(seed_label, "global-cap"))
        keep = np.sort(rng.choice(len(x), size=row_cap, replace=False))
        x, y, weight = x[keep], y[keep], weight[keep]
    require(set(np.unique(y).tolist()) == {0, 1}, "training labels lack a class")
    class_mass = np.asarray(
        [weight[y == value].sum() for value in (0, 1)], dtype=np.float64
    )
    require(np.all(class_mass > EPS), "a class has zero training weight")
    weight = weight / (2.0 * class_mass[y])
    weight *= len(weight) / float(weight.sum())
    return x, y, weight


def fit_estimators(
    x: np.ndarray,
    y: np.ndarray,
    weight: np.ndarray,
    seed: int,
) -> dict[str, Any]:
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    require(len(x) <= TRAINING_ROW_CAP, "estimator received more than 500,000 rows")
    sklearn_seed = int(seed % (2**32 - 1))
    models: dict[str, Any] = {}
    for name, columns in (
        ("full_hgb", np.arange(x.shape[1], dtype=np.int64)),
        ("donor_hgb", DONOR_ONLY_FEATURE_INDICES),
    ):
        estimator = HistGradientBoostingClassifier(
            max_depth=3,
            learning_rate=0.05,
            max_iter=HGB_MAX_ITER,
            l2_regularization=1.0,
            early_stopping=False,
            random_state=sklearn_seed,
        )
        with single_fit_guard(name):
            estimator.fit(x[:, columns], y, sample_weight=weight)
        models[name] = {
            "kind": "hgb",
            "estimator": estimator,
            "columns": columns,
        }
    for name, columns in (
        ("full_logistic", np.arange(x.shape[1], dtype=np.int64)),
        ("donor_logistic", DONOR_ONLY_FEATURE_INDICES),
    ):
        scaler = StandardScaler()
        scaled = scaler.fit_transform(x[:, columns], sample_weight=weight)
        estimator = LogisticRegression(
            C=1.0,
            solver="lbfgs",
            max_iter=200,
            tol=1e-6,
            random_state=sklearn_seed,
        )
        with single_fit_guard(name):
            estimator.fit(scaled, y, sample_weight=weight)
        models[name] = {
            "kind": "logistic",
            "estimator": estimator,
            "scaler": scaler,
            "columns": columns,
        }
    return models


def fit_full_hgb(
    x: np.ndarray,
    y: np.ndarray,
    weight: np.ndarray,
    seed: int,
) -> dict[str, Any]:
    from sklearn.ensemble import HistGradientBoostingClassifier

    require(len(x) <= TRAINING_ROW_CAP, "estimator received more than 500,000 rows")
    columns = np.arange(x.shape[1], dtype=np.int64)
    estimator = HistGradientBoostingClassifier(
        max_depth=3,
        learning_rate=0.05,
        max_iter=HGB_MAX_ITER,
        l2_regularization=1.0,
        early_stopping=False,
        random_state=int(seed % (2**32 - 1)),
    )
    with single_fit_guard("full_hgb_permutation"):
        estimator.fit(x, y, sample_weight=weight)
    return {
        "kind": "hgb",
        "estimator": estimator,
        "columns": columns,
    }


def predict_wrong_probability(
    model: Any, name: str, features: np.ndarray
) -> np.ndarray:
    columns = model["columns"]
    if model["kind"] == "hgb":
        probabilities = model["estimator"].predict_proba(features[:, columns])[:, 1]
    else:
        probabilities = model["estimator"].predict_proba(
            model["scaler"].transform(features[:, columns])
        )[:, 1]
    require(np.all(np.isfinite(probabilities)), f"{name}: nonfinite probabilities")
    return np.asarray(probabilities, dtype=np.float64)


def predict_probability_matrix(
    *,
    model: Any,
    model_name: str,
    arrays: Mapping[str, np.ndarray],
    context_index: int,
    target_indices: Sequence[int],
    control_state: Mapping[str, np.ndarray],
    control_mapping: np.ndarray,
) -> np.ndarray:
    mapped_controls = arrays["control_profiles"][control_mapping]
    context_scores = transform_controls(mapped_controls, control_state)
    donor_scores = transform_controls(arrays["donor_controls"], control_state)
    output = np.zeros(
        (len(target_indices), arrays["baseline"].shape[1]), dtype=np.float32
    )
    for row, target in enumerate(target_indices):
        features = build_feature_block(
            context_control=mapped_controls[context_index],
            context_score=context_scores[context_index],
            donor_controls=arrays["donor_controls"],
            donor_scores=donor_scores,
            donor_effects=arrays["donor_effects"],
            donor_available=arrays["donor_available"],
            baseline=arrays["baseline"],
            target_index=target,
        )
        output[row] = predict_wrong_probability(model, model_name, features).astype(
            np.float32
        )
    return output


def apply_sign_repair(
    base: np.ndarray,
    probabilities: np.ndarray,
    eligible: np.ndarray,
    candidate: Candidate | None,
    *,
    pds_mask: np.ndarray | None = None,
    threshold: float | None = None,
) -> tuple[np.ndarray, np.ndarray, RepairAudit]:
    base = np.asarray(base, dtype=np.float64)
    probabilities = np.asarray(probabilities, dtype=np.float64)
    eligible = np.asarray(eligible, dtype=bool)
    require(
        base.shape == probabilities.shape == eligible.shape, "repair shape mismatch"
    )
    pds = (
        np.ones_like(eligible) if pds_mask is None else np.asarray(pds_mask, dtype=bool)
    )
    require(pds.shape == base.shape, "PDS-axis mask shape mismatch")
    if candidate is None:
        return (
            base.copy(),
            np.zeros(len(base), dtype=bool),
            RepairAudit(
                threshold=float("inf"),
                selected_count=0,
                energy_fraction=0.0,
                pds_axis_energy_fraction=0.0,
                signed_l1_mass_fraction=0.0,
                angle_degrees=0.0,
            ),
        )
    norm2 = float(base @ base)
    if norm2 <= EPS:
        return (
            base.copy(),
            np.zeros(len(base), dtype=bool),
            RepairAudit(
                threshold=float("inf"),
                selected_count=0,
                energy_fraction=0.0,
                pds_axis_energy_fraction=0.0,
                signed_l1_mass_fraction=0.0,
                angle_degrees=0.0,
            ),
        )
    active = eligible & (base != 0) & np.isfinite(probabilities)
    if threshold is None:
        values = probabilities[active]
        threshold = (
            float(np.quantile(values, candidate.score_quantile, method="higher"))
            if len(values)
            else float("inf")
        )
    possible = np.flatnonzero(active & (probabilities >= threshold))
    order = sorted(possible.tolist(), key=lambda index: (-probabilities[index], index))
    limit = candidate.energy_cap * norm2
    pds_norm2 = float(base[pds] @ base[pds])
    pds_limit = candidate.pds_axis_energy_cap * pds_norm2
    signed_l1 = float(np.sum(np.abs(base[eligible])))
    signed_l1_limit = candidate.signed_l1_mass_cap * signed_l1
    selected: list[int] = []
    used = 0.0
    pds_used = 0.0
    l1_used = 0.0
    for index in order:
        cost = float(base[index] ** 2)
        pds_cost = cost if pds[index] else 0.0
        l1_cost = float(abs(base[index]))
        if len(selected) >= candidate.selected_count_cap:
            break
        if used + cost > limit + EPS:
            continue
        if pds_used + pds_cost > pds_limit + EPS:
            continue
        if l1_used + l1_cost > signed_l1_limit + EPS:
            continue
        selected.append(index)
        used += cost
        pds_used += pds_cost
        l1_used += l1_cost
    output = base.copy()
    mask = np.zeros(len(base), dtype=bool)
    if selected:
        indices = np.asarray(selected, dtype=np.int64)
        output[indices] *= -1.0
        mask[indices] = True
    require(np.array_equal(np.abs(output), np.abs(base)), "repair changed magnitude")
    require(
        np.isclose(output @ output, norm2, rtol=0.0, atol=1e-10), "repair changed norm"
    )
    energy_fraction = used / norm2
    require(
        energy_fraction <= candidate.energy_cap + 1e-12,
        "repair exceeded energy cap",
    )
    pds_energy_fraction = pds_used / pds_norm2 if pds_norm2 > EPS else 0.0
    require(
        pds_energy_fraction <= candidate.pds_axis_energy_cap + 1e-12,
        "repair exceeded PDS-axis energy cap",
    )
    signed_l1_fraction = l1_used / signed_l1 if signed_l1 > EPS else 0.0
    require(
        signed_l1_fraction <= candidate.signed_l1_mass_cap + 1e-12,
        "repair exceeded signed-L1-mass cap",
    )
    require(
        len(selected) <= candidate.selected_count_cap,
        "repair exceeded selected-count cap",
    )
    cosine = float(base @ output / norm2)
    angle = math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))
    return (
        output,
        mask,
        RepairAudit(
            threshold=float(threshold),
            selected_count=len(selected),
            energy_fraction=energy_fraction,
            pds_axis_energy_fraction=pds_energy_fraction,
            signed_l1_mass_fraction=signed_l1_fraction,
            angle_degrees=angle,
        ),
    )


def pds_per_target(
    prediction: np.ndarray,
    truth: np.ndarray,
    mask: np.ndarray,
) -> np.ndarray:
    pred = np.asarray(prediction[:, mask], dtype=np.float64)
    observed = np.asarray(truth[:, mask], dtype=np.float64)
    require(pred.shape == observed.shape, "PDS prediction/truth shape mismatch")
    require(len(pred) >= 2, "PDS gallery requires at least two rows")
    pred_norm = np.linalg.norm(pred, axis=1, keepdims=True)
    truth_norm = np.linalg.norm(observed, axis=1, keepdims=True)
    pred = np.divide(pred, pred_norm, out=np.zeros_like(pred), where=pred_norm > EPS)
    observed = np.divide(
        observed,
        truth_norm,
        out=np.zeros_like(observed),
        where=truth_norm > EPS,
    )
    similarity = pred @ observed.T
    own = np.diag(similarity)[:, None]
    greater = np.sum(similarity > own, axis=1)
    equal = np.sum(similarity == own, axis=1)
    midrank = 1.0 + greater + 0.5 * (equal - 1)
    return 1.0 - (midrank - 1.0) / (len(pred) - 1.0)


def signed_credit_per_target(
    prediction: np.ndarray,
    truth: np.ndarray,
    masks: np.ndarray,
) -> np.ndarray:
    require(
        prediction.shape == truth.shape == masks.shape, "signed-credit shape mismatch"
    )
    numerator = np.sum(np.where(masks, np.sign(truth) * prediction, 0.0), axis=1)
    denominator = np.sum(np.where(masks, np.abs(prediction), 0.0), axis=1)
    return np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator, dtype=np.float64),
        where=denominator > EPS,
    )


def weighted_flip_precision(
    base: np.ndarray,
    truth: np.ndarray,
    selected: np.ndarray,
    masks: np.ndarray,
) -> float:
    active = selected & masks & (base != 0) & (truth != 0)
    weights = np.where(active, np.abs(base), 0.0)
    denominator = float(weights.sum())
    if denominator <= EPS:
        return float("nan")
    wrong = np.sign(base) != np.sign(truth)
    return float(np.sum(weights * wrong) / denominator)


def repair_probability_matrix(
    *,
    arrays: Mapping[str, np.ndarray],
    target_indices: Sequence[int],
    probabilities: np.ndarray,
    candidate: Candidate | None,
    action: str = "flip",
) -> tuple[np.ndarray, np.ndarray, tuple[RepairAudit, ...]]:
    require(action in {"flip", "zero"}, f"unknown repair action: {action}")
    targets = np.asarray(target_indices, dtype=np.int64)
    require(
        probabilities.shape == (len(targets), arrays["baseline"].shape[1]),
        "probability matrix shape mismatch",
    )
    supports = arrays["donor_available"].sum(axis=0)
    eligible_rows = np.asarray(
        [
            arrays["signed_masks"][target] & (supports[target] >= 2)
            for target in targets
        ],
        dtype=bool,
    )
    score_parts = [
        probabilities[row][eligible]
        for row, eligible in enumerate(eligible_rows)
        if np.any(eligible)
    ]
    eligible_scores = (
        np.concatenate(score_parts) if score_parts else np.asarray([], dtype=np.float64)
    )
    threshold = (
        float(
            np.quantile(
                eligible_scores,
                candidate.score_quantile,
                method="higher",
            )
        )
        if candidate is not None and len(eligible_scores)
        else float("inf")
    )
    repaired_rows = []
    selected_rows = []
    audits = []
    for row, target in enumerate(targets):
        repaired, selected, audit = apply_sign_repair(
            arrays["baseline"][target],
            probabilities[row],
            eligible_rows[row],
            candidate,
            pds_mask=arrays["pds_mask"],
            threshold=threshold,
        )
        if action == "zero":
            repaired = np.asarray(arrays["baseline"][target], dtype=np.float64).copy()
            repaired[selected] = 0.0
        repaired_rows.append(repaired)
        selected_rows.append(selected)
        audits.append(audit)
    return (
        np.asarray(repaired_rows, dtype=np.float32),
        np.asarray(selected_rows, dtype=bool),
        tuple(audits),
    )


def score_candidate_context(
    *,
    arrays: Mapping[str, np.ndarray],
    context_index: int,
    target_indices: Sequence[int],
    gallery_indices: Sequence[int],
    probabilities: np.ndarray,
    candidate: Candidate | None,
    fdr_threshold: float = TRAINING_FDR,
    action: str = "flip",
) -> dict[str, Any]:
    target_indices = np.asarray(target_indices, dtype=np.int64)
    gallery_indices = np.asarray(gallery_indices, dtype=np.int64)
    available_gallery = gallery_indices[
        arrays["available"][context_index, gallery_indices]
    ]
    available_queries = target_indices[
        arrays["available"][context_index, target_indices]
    ]
    require(len(available_queries) >= 2, "too few available queries")
    require(
        set(available_queries).issubset(set(available_gallery)),
        "query absent from gallery",
    )
    probability_lookup = {
        int(target): probabilities[row] for row, target in enumerate(target_indices)
    }
    available_probabilities = np.asarray(
        [probability_lookup[int(target)] for target in available_queries],
        dtype=np.float32,
    )
    candidate_matrix, selected_matrix, audits = repair_probability_matrix(
        arrays=arrays,
        target_indices=available_queries,
        probabilities=available_probabilities,
        candidate=candidate,
        action=action,
    )
    base_matrix = arrays["baseline"][available_queries]
    truth_matrix = arrays["effects"][context_index, available_queries]
    signed_masks = np.asarray(
        [
            arrays["signed_masks"][target]
            & np.isfinite(arrays["fdr"][context_index, target])
            & (arrays["fdr"][context_index, target] <= fdr_threshold)
            for target in available_queries
        ],
        dtype=bool,
    )
    gallery_truth = arrays["effects"][context_index, available_gallery]
    gallery_position = {
        int(target): position for position, target in enumerate(available_gallery)
    }

    def pds_queries(matrix: np.ndarray) -> np.ndarray:
        pred = np.asarray(matrix[:, arrays["pds_mask"]], dtype=np.float64)
        truth = np.asarray(gallery_truth[:, arrays["pds_mask"]], dtype=np.float64)
        pred_norm = np.linalg.norm(pred, axis=1, keepdims=True)
        truth_norm = np.linalg.norm(truth, axis=1, keepdims=True)
        pred = np.divide(
            pred, pred_norm, out=np.zeros_like(pred), where=pred_norm > EPS
        )
        truth = np.divide(
            truth, truth_norm, out=np.zeros_like(truth), where=truth_norm > EPS
        )
        similarity = pred @ truth.T
        own_columns = np.asarray(
            [gallery_position[int(target)] for target in available_queries],
            dtype=np.int64,
        )
        own = similarity[np.arange(len(available_queries)), own_columns]
        greater = np.sum(similarity > own[:, None], axis=1)
        equal = np.sum(similarity == own[:, None], axis=1)
        midrank = 1.0 + greater + 0.5 * (equal - 1)
        return 1.0 - (midrank - 1.0) / (len(available_gallery) - 1.0)

    candidate_pds = pds_queries(candidate_matrix)
    base_pds = pds_queries(base_matrix)
    candidate_credit = signed_credit_per_target(
        candidate_matrix, truth_matrix, signed_masks
    )
    base_credit = signed_credit_per_target(base_matrix, truth_matrix, signed_masks)
    return {
        "targets": available_queries,
        "gallery_size": int(len(available_gallery)),
        "candidate": candidate_matrix,
        "selected": selected_matrix,
        "pds_delta": candidate_pds - base_pds,
        "signed_credit_delta": candidate_credit - base_credit,
        "flip_precision": weighted_flip_precision(
            base_matrix,
            truth_matrix,
            selected_matrix,
            signed_masks,
        ),
        "selected_count": int(selected_matrix.sum()),
        "max_selected_per_row": max(
            (audit.selected_count for audit in audits), default=0
        ),
        "max_energy_fraction": max(
            (audit.energy_fraction for audit in audits), default=0.0
        ),
        "max_pds_axis_energy_fraction": max(
            (audit.pds_axis_energy_fraction for audit in audits), default=0.0
        ),
        "max_signed_l1_mass_fraction": max(
            (audit.signed_l1_mass_fraction for audit in audits), default=0.0
        ),
        "max_angle_degrees": max(
            (audit.angle_degrees for audit in audits), default=0.0
        ),
        "score_threshold": audits[0].threshold if audits else float("inf"),
    }


def donor_agreement_probability(
    arrays: Mapping[str, np.ndarray],
    target_indices: Sequence[int],
) -> np.ndarray:
    output = np.zeros(
        (len(target_indices), arrays["baseline"].shape[1]), dtype=np.float32
    )
    for row, target in enumerate(target_indices):
        have = np.flatnonzero(arrays["donor_available"][:, target])
        require(len(have) >= 1, "donor-agreement arm has no donor")
        if len(have) == 1:
            continue
        stack = arrays["donor_effects"][have, target]
        agreement = np.mean(
            np.sign(stack) == np.sign(arrays["baseline"][target])[None, :],
            axis=0,
        )
        output[row] = (1.0 - agreement).astype(np.float32)
    return output


def select_fixed_primary(
    records: Mapping[str, Sequence[Mapping[str, Any]]],
) -> tuple[Candidate | None, dict[str, Any]]:
    summaries: dict[str, Any] = {}
    for arm in ALL_ARMS:
        local_records = list(records[arm])
        require(bool(local_records), f"no inner records for {arm}")
        batch_names = sorted({str(record["batch"]) for record in local_records})
        batch_credit = {
            batch: float(
                np.mean(
                    [
                        record["signed_credit_delta"]
                        for record in local_records
                        if record["batch"] == batch
                    ]
                )
            )
            for batch in batch_names
        }
        batch_pds = {
            batch: float(
                np.mean(
                    [
                        record["pds_delta"]
                        for record in local_records
                        if record["batch"] == batch
                    ]
                )
            )
            for batch in batch_names
        }
        summaries[arm] = {
            "batch_signed_credit_delta": batch_credit,
            "batch_pds_delta": batch_pds,
            "macro_signed_credit_delta": float(np.mean(list(batch_credit.values()))),
            "macro_pds_delta": float(np.mean(list(batch_pds.values()))),
        }
    primary = summaries["full_hgb"]
    enabled = (
        primary["macro_signed_credit_delta"] > 0
        and min(primary["batch_signed_credit_delta"].values()) >= -EPS
        and primary["macro_pds_delta"] >= -0.002
        and primary["macro_signed_credit_delta"]
        >= summaries["donor_hgb"]["macro_signed_credit_delta"] - EPS
    )
    return (
        FIXED_CANDIDATE if enabled else None,
        {
            "selected": FIXED_CANDIDATE.key if enabled else "identity",
            "fixed_candidate": asdict(FIXED_CANDIDATE),
            "summaries": summaries,
            "full_hgb_not_below_donor_hgb": (
                primary["macro_signed_credit_delta"]
                >= summaries["donor_hgb"]["macro_signed_credit_delta"] - EPS
            ),
        },
    )


def inner_complete_gallery_evaluation(
    *,
    outer_batch: int,
    outer_fold: int,
    arrays: Mapping[str, np.ndarray],
    plan: Mapping[str, Any],
    control_state: Mapping[str, np.ndarray],
) -> tuple[Candidate | None, dict[str, Any]]:
    batch_folds = experimental_batch_folds(plan["batches"])
    held_contexts = set(batch_folds[outer_batch].tolist())
    train_contexts = np.asarray(
        [index for index in range(CONTEXT_COUNT) if index not in held_contexts],
        dtype=np.int64,
    )
    target_folds = arrays["target_folds"]
    inner_targets = np.flatnonzero(target_folds != outer_fold)
    require(
        len(inner_targets) in (64, 65),
        "inner gallery must contain exactly 64 or 65 targets",
    )
    selection_contexts = complete_gallery_contexts(
        arrays["available"],
        train_contexts,
        inner_targets,
    )
    target_position = {
        int(target): position for position, target in enumerate(inner_targets)
    }
    identity_mapping = np.arange(CONTEXT_COUNT, dtype=np.int64)
    probabilities = {
        arm: {
            int(context): np.zeros(
                (len(inner_targets), arrays["baseline"].shape[1]),
                dtype=np.float32,
            )
            for context in selection_contexts
        }
        for arm in MODEL_ARMS
    }

    for validation_fold in range(TARGET_FOLD_COUNT):
        if validation_fold == outer_fold:
            continue
        fit_targets = np.flatnonzero(
            (target_folds != outer_fold) & (target_folds != validation_fold)
        )
        validation_targets = np.flatnonzero(target_folds == validation_fold)
        require(
            not (set(fit_targets) & set(validation_targets)),
            "inner fit and validation targets overlap",
        )
        x, y, weight = training_matrix(
            arrays=arrays,
            plan=plan,
            contexts=train_contexts,
            targets=fit_targets,
            control_state=control_state,
            control_mapping=identity_mapping,
            seed_label=(
                "inner",
                outer_batch,
                outer_fold,
                validation_fold,
            ),
        )
        models = fit_estimators(
            x,
            y,
            weight,
            stable_seed(
                "inner-model",
                outer_batch,
                outer_fold,
                validation_fold,
            ),
        )
        positions = np.asarray(
            [target_position[int(target)] for target in validation_targets],
            dtype=np.int64,
        )
        for context in selection_contexts:
            for arm in MODEL_ARMS:
                probabilities[arm][int(context)][positions] = (
                    predict_probability_matrix(
                        model=models[arm],
                        model_name=arm,
                        arrays=arrays,
                        context_index=int(context),
                        target_indices=validation_targets,
                        control_state=control_state,
                        control_mapping=identity_mapping,
                    )
                )
        del x, y, weight, models

    donor_scores = donor_agreement_probability(arrays, inner_targets)
    records: dict[str, list[dict[str, Any]]] = {arm: [] for arm in ALL_ARMS}
    for context in selection_contexts:
        for arm in ALL_ARMS:
            if arm == "zero_tail":
                arm_probabilities = probabilities["full_hgb"][int(context)]
                candidate = FIXED_CANDIDATE
                action = "zero"
            elif arm == "donor_agreement":
                arm_probabilities = donor_scores
                candidate = FIXED_CANDIDATE
                action = "flip"
            else:
                arm_probabilities = probabilities[arm][int(context)]
                candidate = FIXED_CANDIDATE
                action = "flip"
            scored = score_candidate_context(
                arrays=arrays,
                context_index=int(context),
                target_indices=inner_targets,
                gallery_indices=inner_targets,
                probabilities=arm_probabilities,
                candidate=candidate,
                action=action,
            )
            require(
                scored["gallery_size"] == len(inner_targets),
                "inner selection did not use the complete 64/65-target gallery",
            )
            records[arm].append(
                {
                    "context": plan["contexts"][context],
                    "batch": plan["batches"][context],
                    "gallery_size": len(
                        np.flatnonzero(arrays["available"][context, inner_targets])
                    ),
                    "pds_delta": float(np.mean(scored["pds_delta"])),
                    "signed_credit_delta": float(
                        np.mean(scored["signed_credit_delta"])
                    ),
                }
            )
    selected, audit = select_fixed_primary(records)
    audit["outer_batch"] = sorted(set(plan["batches"]))[outer_batch]
    audit["outer_fold"] = outer_fold
    audit["inner_gallery_size"] = len(inner_targets)
    audit["selection_contexts"] = [
        plan["contexts"][int(context)] for context in selection_contexts
    ]
    audit["excluded_incomplete_gallery_contexts"] = [
        plan["contexts"][int(context)]
        for context in train_contexts
        if context not in set(selection_contexts.tolist())
    ]
    audit["represented_selection_batches"] = sorted(
        {plan["batches"][int(context)] for context in selection_contexts}
    )
    audit["mini_gallery_scoring_used"] = False
    return selected, audit


def outer_batch_run(
    *,
    outer_batch: int,
    arrays: Mapping[str, np.ndarray],
    plan: Mapping[str, Any],
) -> dict[str, Any]:
    batch_folds = experimental_batch_folds(plan["batches"])
    held_contexts = batch_folds[outer_batch]
    held_set = set(held_contexts.tolist())
    train_contexts = np.asarray(
        [index for index in range(CONTEXT_COUNT) if index not in held_set],
        dtype=np.int64,
    )
    require(
        not (set(train_contexts) & held_set),
        "held experimental batch entered model fitting",
    )
    target_folds = arrays["target_folds"]
    identity_mapping = np.arange(CONTEXT_COUNT, dtype=np.int64)
    control_state = fit_control_pca(
        arrays["control_profiles"],
        plan["batches"],
        train_contexts,
    )
    base_cube = np.broadcast_to(
        arrays["baseline"],
        (len(held_contexts), *arrays["baseline"].shape),
    )
    crossfit = {arm: np.asarray(base_cube, dtype=np.float32).copy() for arm in ALL_ARMS}
    selected_masks = {arm: np.zeros_like(base_cube, dtype=bool) for arm in ALL_ARMS}
    permutation_crossfit = np.broadcast_to(
        arrays["baseline"],
        (
            PERMUTATION_NULL_DRAWS,
            len(held_contexts),
            *arrays["baseline"].shape,
        ),
    ).astype(np.float32, copy=True)
    repair_audits: dict[str, list[dict[str, Any]]] = {arm: [] for arm in ALL_ARMS}
    fold_records = []

    for outer_fold in range(TARGET_FOLD_COUNT):
        held_targets = np.flatnonzero(target_folds == outer_fold)
        refit_targets = np.flatnonzero(target_folds != outer_fold)
        require(
            not (set(held_targets) & set(refit_targets)),
            "outer target entered final model fit",
        )
        selected_primary, inner_audit = inner_complete_gallery_evaluation(
            outer_batch=outer_batch,
            outer_fold=outer_fold,
            arrays=arrays,
            plan=plan,
            control_state=control_state,
        )
        x, y, weight = training_matrix(
            arrays=arrays,
            plan=plan,
            contexts=train_contexts,
            targets=refit_targets,
            control_state=control_state,
            control_mapping=identity_mapping,
            seed_label=("refit", outer_batch, outer_fold),
        )
        models = fit_estimators(
            x,
            y,
            weight,
            stable_seed("refit-model", outer_batch, outer_fold),
        )
        donor_scores = donor_agreement_probability(arrays, held_targets)
        for context_position, context in enumerate(held_contexts):
            arm_probabilities = {
                arm: predict_probability_matrix(
                    model=models[arm],
                    model_name=arm,
                    arrays=arrays,
                    context_index=int(context),
                    target_indices=held_targets,
                    control_state=control_state,
                    control_mapping=identity_mapping,
                )
                for arm in MODEL_ARMS
            }
            arm_probabilities["donor_agreement"] = donor_scores
            arm_probabilities["zero_tail"] = arm_probabilities["full_hgb"]
            for arm in ALL_ARMS:
                candidate = (
                    selected_primary
                    if arm == "full_hgb"
                    else selected_primary
                    if arm == "zero_tail"
                    else FIXED_CANDIDATE
                )
                repaired, selected, audits = repair_probability_matrix(
                    arrays=arrays,
                    target_indices=held_targets,
                    probabilities=arm_probabilities[arm],
                    candidate=candidate,
                    action="zero" if arm == "zero_tail" else "flip",
                )
                crossfit[arm][context_position, held_targets] = repaired
                selected_masks[arm][context_position, held_targets] = selected
                repair_audits[arm].extend(asdict(audit) for audit in audits)

        permutations, permutation_audit = availability_stratified_target_permutations(
            refit_targets,
            arrays["available"],
            train_contexts,
            seed=stable_seed(
                "whole-target-permutation",
                outer_batch,
                outer_fold,
                base=PERMUTATION_NULL_SEED,
            ),
        )
        for draw, permutation in enumerate(permutations):
            null_x, null_y, null_weight = training_matrix(
                arrays=arrays,
                plan=plan,
                contexts=train_contexts,
                targets=refit_targets,
                control_state=control_state,
                control_mapping=identity_mapping,
                seed_label=("permutation-refit", outer_batch, outer_fold, draw),
                target_truth_permutation=permutation,
            )
            null_model = fit_full_hgb(
                null_x,
                null_y,
                null_weight,
                stable_seed("permutation-model", outer_batch, outer_fold, draw),
            )
            for context_position, context in enumerate(held_contexts):
                probabilities = predict_probability_matrix(
                    model=null_model,
                    model_name="full_hgb",
                    arrays=arrays,
                    context_index=int(context),
                    target_indices=held_targets,
                    control_state=control_state,
                    control_mapping=identity_mapping,
                )
                repaired, _, _ = repair_probability_matrix(
                    arrays=arrays,
                    target_indices=held_targets,
                    probabilities=probabilities,
                    candidate=selected_primary,
                )
                permutation_crossfit[draw, context_position, held_targets] = repaired
            del null_x, null_y, null_weight, null_model

        fold_records.append(
            {
                "outer_fold": outer_fold,
                "held_targets": [plan["targets"][index] for index in held_targets],
                "held_target_count": int(len(held_targets)),
                "fit_target_count": int(len(refit_targets)),
                "inner_gallery_size": int(len(refit_targets)),
                "inner_selection": inner_audit,
                "whole_target_permutation_refits": PERMUTATION_NULL_DRAWS,
                "permutation_availability_audit": {
                    **permutation_audit,
                    "classes": [
                        {
                            **record,
                            "targets": [
                                plan["targets"][index]
                                for index in record["target_indices"]
                            ],
                        }
                        for record in permutation_audit["classes"]
                    ],
                    "singleton_targets": [
                        plan["targets"][index]
                        for index in permutation_audit["singleton_target_indices"]
                    ],
                },
            }
        )
        del x, y, weight, models

    return {
        "outer_batch_index": outer_batch,
        "outer_batch": sorted(set(plan["batches"]))[outer_batch],
        "held_contexts": held_contexts,
        "predictions": crossfit,
        "selected_masks": selected_masks,
        "permutation_predictions": permutation_crossfit,
        "repair_audits": repair_audits,
        "fold_records": fold_records,
    }


def equal_batch_macro(values: np.ndarray, batches: Sequence[str]) -> float:
    batch_names = sorted(set(batches))
    return float(
        np.mean(
            [
                np.mean(
                    [
                        values[index]
                        for index, batch_value in enumerate(batches)
                        if batch_value == batch
                    ]
                )
                for batch in batch_names
            ]
        )
    )


def equal_batch_target_metric(
    records: Sequence[Mapping[str, Any]],
    *,
    arm: str,
    metric: str,
    target_count: int = TARGET_COUNT,
) -> np.ndarray:
    """Aggregate each target with the same equal-batch estimand as the macro score."""
    target_batch_values: list[dict[str, list[float]]] = [
        {} for _ in range(target_count)
    ]
    for record in records:
        available = np.asarray(record["available_targets"], dtype=np.int64)
        values = np.asarray(record["metrics"][arm][metric], dtype=np.float64)
        require(len(available) == len(values), "target metric alignment changed")
        for target, value in zip(available, values, strict=True):
            target_batch_values[int(target)].setdefault(
                str(record["batch"]), []
            ).append(float(value))
    result = np.asarray(
        [
            np.mean(
                [
                    np.mean(batch_values)
                    for batch_values in target_batch_values[target].values()
                ]
            )
            if target_batch_values[target]
            else np.nan
            for target in range(target_count)
        ],
        dtype=np.float64,
    )
    require(np.all(np.isfinite(result)), "a target has no held batch score")
    return result


def batch_bootstrap(
    values: np.ndarray,
    batches: Sequence[str],
    *,
    draws: int = BOOTSTRAP_DRAWS,
    seed: int = BOOTSTRAP_SEED,
) -> np.ndarray:
    batch_names = sorted(set(batches))
    batch_values = np.asarray(
        [
            np.mean(
                [
                    values[index]
                    for index, batch_value in enumerate(batches)
                    if batch_value == batch
                ]
            )
            for batch in batch_names
        ],
        dtype=np.float64,
    )
    rng = np.random.default_rng(seed)
    sampled = rng.integers(0, len(batch_values), size=(draws, len(batch_values)))
    return np.mean(batch_values[sampled], axis=1)


def target_block_bootstrap(
    target_values: np.ndarray,
    target_folds: np.ndarray,
    *,
    draws: int = BOOTSTRAP_DRAWS,
    seed: int = TARGET_BOOTSTRAP_SEED,
) -> np.ndarray:
    fold_values = np.asarray(
        [
            np.mean(target_values[target_folds == fold])
            for fold in range(TARGET_FOLD_COUNT)
        ],
        dtype=np.float64,
    )
    require(np.all(np.isfinite(fold_values)), "target-block bootstrap is nonfinite")
    rng = np.random.default_rng(seed)
    sampled = rng.integers(
        0,
        TARGET_FOLD_COUNT,
        size=(draws, TARGET_FOLD_COUNT),
    )
    return np.mean(fold_values[sampled], axis=1)


def summarize_evaluation(
    outer_results: Sequence[Mapping[str, Any]],
    arrays: Mapping[str, np.ndarray],
    plan: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    require(len(outer_results) == BATCH_COUNT, "expected six outer batch folds")
    context_records: list[dict[str, Any] | None] = [None] * CONTEXT_COUNT
    for batch_record in outer_results:
        held_contexts = np.asarray(batch_record["held_contexts"], dtype=np.int64)
        for position, context in enumerate(held_contexts):
            require(context_records[int(context)] is None, "context scored twice")
            available = np.flatnonzero(arrays["available"][context])
            truth = arrays["effects"][context, available]
            base = arrays["baseline"][available]
            record: dict[str, Any] = {
                "context_index": int(context),
                "context": plan["contexts"][context],
                "batch": plan["batches"][context],
                "available_targets": available,
                "metrics": {},
                "fdr_sensitivity": {},
                "null_signed_credit_delta": np.empty(
                    PERMUTATION_NULL_DRAWS, dtype=np.float64
                ),
            }
            fixed_masks = (
                arrays["signed_masks"][available]
                & np.isfinite(arrays["fdr"][context, available])
                & (arrays["fdr"][context, available] <= TRAINING_FDR)
            )
            base_pds = pds_per_target(base, truth, arrays["pds_mask"])
            base_credit = signed_credit_per_target(base, truth, fixed_masks)
            for arm in ALL_ARMS:
                prediction = batch_record["predictions"][arm][position, available]
                selected = batch_record["selected_masks"][arm][position, available]
                record["metrics"][arm] = {
                    "pds_delta": (
                        pds_per_target(prediction, truth, arrays["pds_mask"]) - base_pds
                    ),
                    "signed_credit_delta": (
                        signed_credit_per_target(prediction, truth, fixed_masks)
                        - base_credit
                    ),
                    "flip_precision": weighted_flip_precision(
                        base,
                        truth,
                        selected,
                        fixed_masks,
                    ),
                    "selected_count": int(selected.sum()),
                }
            full_prediction = batch_record["predictions"]["full_hgb"][
                position, available
            ]
            for threshold in FDR_SENSITIVITY:
                masks = (
                    arrays["signed_masks"][available]
                    & np.isfinite(arrays["fdr"][context, available])
                    & (arrays["fdr"][context, available] <= threshold)
                )
                record["fdr_sensitivity"][str(threshold)] = signed_credit_per_target(
                    full_prediction, truth, masks
                ) - signed_credit_per_target(base, truth, masks)
            for draw in range(PERMUTATION_NULL_DRAWS):
                null_prediction = batch_record["permutation_predictions"][
                    draw, position, available
                ]
                record["null_signed_credit_delta"][draw] = float(
                    np.mean(
                        signed_credit_per_target(
                            null_prediction,
                            truth,
                            fixed_masks,
                        )
                        - base_credit
                    )
                )
            context_records[int(context)] = record
    require(all(record is not None for record in context_records), "context omitted")
    records = [record for record in context_records if record is not None]
    batches = [record["batch"] for record in records]
    contexts = [record["context"] for record in records]
    summaries: dict[str, Any] = {}
    numeric: dict[str, np.ndarray] = {}
    for name in ALL_ARMS:
        context_pds = np.asarray(
            [np.mean(record["metrics"][name]["pds_delta"]) for record in records]
        )
        context_credit = np.asarray(
            [
                np.mean(record["metrics"][name]["signed_credit_delta"])
                for record in records
            ]
        )
        context_precision = np.asarray(
            [record["metrics"][name]["flip_precision"] for record in records]
        )
        bootstrap = batch_bootstrap(
            context_credit,
            batches,
            seed=BOOTSTRAP_SEED + ALL_ARMS.index(name),
        )
        summaries[name] = {
            "macro_pds_delta": equal_batch_macro(context_pds, batches),
            "macro_signed_credit_delta": equal_batch_macro(context_credit, batches),
            "signed_credit_batch_bootstrap_95": [
                float(np.quantile(bootstrap, 0.025)),
                float(np.quantile(bootstrap, 0.975)),
            ],
            "positive_contexts": int(np.sum(context_credit > 0)),
            "contexts_with_flip_precision_at_least_0_65": int(
                np.sum(context_precision >= 0.65)
            ),
            "minimum_context_pds_delta": float(np.min(context_pds)),
            "context_records": [
                {
                    "context": context,
                    "batch": batch,
                    "pds_delta": float(pds),
                    "signed_credit_delta": float(credit),
                    "flip_precision": (
                        None if not np.isfinite(precision) else float(precision)
                    ),
                }
                for context, batch, pds, credit, precision in zip(
                    contexts,
                    batches,
                    context_pds,
                    context_credit,
                    context_precision,
                    strict=True,
                )
            ],
        }
        numeric[f"{name}_context_pds_delta"] = context_pds
        numeric[f"{name}_context_signed_credit_delta"] = context_credit
        numeric[f"{name}_context_flip_precision"] = context_precision
        numeric[f"{name}_signed_credit_bootstrap"] = bootstrap

    primary_credit = numeric["full_hgb_context_signed_credit_delta"]
    target_credit = equal_batch_target_metric(
        records,
        arm="full_hgb",
        metric="signed_credit_delta",
    )
    target_fold_credit = np.zeros(TARGET_FOLD_COUNT, dtype=np.float64)
    for fold in range(TARGET_FOLD_COUNT):
        target_fold_credit[fold] = float(
            np.mean(target_credit[arrays["target_folds"] == fold])
        )
    target_bootstrap = target_block_bootstrap(
        target_credit,
        arrays["target_folds"],
    )
    numeric["full_hgb_target_signed_credit_delta"] = target_credit
    numeric["full_hgb_target_fold_signed_credit_delta"] = target_fold_credit
    numeric["full_hgb_target_block_bootstrap"] = target_bootstrap

    null_by_context = np.vstack(
        [record["null_signed_credit_delta"] for record in records]
    ).T
    null_macro = np.asarray(
        [equal_batch_macro(row, batches) for row in null_by_context],
        dtype=np.float64,
    )
    numeric["whole_target_permutation_null_macro_signed_credit_delta"] = null_macro

    control_differences = {}
    for control in ALL_ARMS:
        if control == "full_hgb":
            continue
        difference = primary_credit - numeric[f"{control}_context_signed_credit_delta"]
        bootstrap = batch_bootstrap(
            difference,
            batches,
            seed=BOOTSTRAP_SEED + 100 + len(control_differences),
        )
        control_differences[control] = {
            "macro_signed_credit_difference": equal_batch_macro(difference, batches),
            "batch_bootstrap_95": [
                float(np.quantile(bootstrap, 0.025)),
                float(np.quantile(bootstrap, 0.975)),
            ],
        }
        numeric[f"primary_minus_{control}_bootstrap"] = bootstrap

    positive = np.maximum(primary_credit, 0.0)
    concentration = (
        float(np.max(positive) / np.sum(positive))
        if float(np.sum(positive)) > EPS
        else 1.0
    )
    observed = summaries["full_hgb"]["macro_signed_credit_delta"]
    null_p95 = float(np.quantile(null_macro, 0.95))
    fdr_sensitivity = {}
    for threshold in FDR_SENSITIVITY:
        values = np.asarray(
            [np.mean(record["fdr_sensitivity"][str(threshold)]) for record in records],
            dtype=np.float64,
        )
        fdr_sensitivity[str(threshold)] = equal_batch_macro(values, batches)
        numeric[f"full_hgb_fdr_{threshold}_context_signed_credit_delta"] = values

    full_hgb_audits = [
        audit
        for batch_record in outer_results
        for audit in batch_record["repair_audits"]["full_hgb"]
    ]
    max_energy = max(
        (audit["energy_fraction"] for audit in full_hgb_audits), default=0.0
    )
    max_pds_energy = max(
        (audit["pds_axis_energy_fraction"] for audit in full_hgb_audits),
        default=0.0,
    )
    max_selected = max(
        (audit["selected_count"] for audit in full_hgb_audits), default=0
    )
    max_signed_l1 = max(
        (audit["signed_l1_mass_fraction"] for audit in full_hgb_audits),
        default=0.0,
    )
    checks = {
        "signed_credit_bootstrap_lower_positive": (
            summaries["full_hgb"]["signed_credit_batch_bootstrap_95"][0] > 0
        ),
        "target_block_bootstrap_lower_positive": (
            float(np.quantile(target_bootstrap, 0.025)) > 0
        ),
        "at_least_12_positive_contexts": summaries["full_hgb"]["positive_contexts"]
        >= 12,
        "every_target_fold_signed_credit_nonnegative": bool(
            np.all(target_fold_credit >= -EPS)
        ),
        "at_least_12_contexts_flip_precision_0_65": (
            summaries["full_hgb"]["contexts_with_flip_precision_at_least_0_65"] >= 12
        ),
        "macro_pds_noninferior": summaries["full_hgb"]["macro_pds_delta"] >= -0.002,
        "minimum_context_pds_noninferior": (
            summaries["full_hgb"]["minimum_context_pds_delta"] >= -0.010
        ),
        "beats_zero_tail_lower_95": (
            control_differences["zero_tail"]["batch_bootstrap_95"][0] > 0
        ),
        "beats_donor_hgb_lower_95": (
            control_differences["donor_hgb"]["batch_bootstrap_95"][0] > 0
        ),
        "beats_other_matched_controls_lower_95": all(
            control_differences[control]["batch_bootstrap_95"][0] > 0
            for control in ("donor_agreement", "full_logistic", "donor_logistic")
        ),
        "permutation_null_p95_below_observed": observed > null_p95,
        "all_fdr_sensitivity_positive": min(fdr_sensitivity.values()) > 0,
        "maximum_positive_context_concentration": concentration <= 0.30,
        "metric_masks_separate": not np.array_equal(
            arrays["signed_masks"][0], arrays["pds_mask"]
        ),
        "target_fold_sizes_exact": tuple(
            int(np.sum(arrays["target_folds"] == fold))
            for fold in range(TARGET_FOLD_COUNT)
        )
        == TARGET_FOLD_SIZES,
        "total_energy_cap": max_energy <= TOTAL_ENERGY_CAP + 1e-12,
        "pds_axis_energy_cap": (max_pds_energy <= PDS_AXIS_ENERGY_CAP + 1e-12),
        "selected_count_cap": max_selected <= SELECTED_COUNT_CAP,
        "signed_l1_mass_cap": max_signed_l1 <= SIGNED_L1_MASS_CAP + 1e-12,
    }
    result = {
        "schema_version": 1,
        "experiment": HERE.name,
        "stage": "source-evaluation",
        "status": (
            "PASS_FOR_SEPARATE_PRODUCTION_INTEGRATION_AUDIT_ONLY"
            if all(checks.values())
            else "NO_GO"
        ),
        "models": summaries,
        "primary_minus_controls": control_differences,
        "target_fold_signed_credit_delta": target_fold_credit.tolist(),
        "target_block_bootstrap_95": [
            float(np.quantile(target_bootstrap, 0.025)),
            float(np.quantile(target_bootstrap, 0.975)),
        ],
        "whole_target_refit_permutation_null": {
            "draws": PERMUTATION_NULL_DRAWS,
            "p95": null_p95,
            "observed": observed,
        },
        "fdr_sensitivity_macro_signed_credit_delta": fdr_sensitivity,
        "repair_caps_observed": {
            "max_energy_fraction": max_energy,
            "max_pds_axis_energy_fraction": max_pds_energy,
            "max_selected_count": max_selected,
            "max_signed_l1_mass_fraction": max_signed_l1,
        },
        "positive_context_gain_concentration": concentration,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
        "development_outcomes_opened": False,
        "challenge_control_expression_opened": False,
        "cell_eval2_used": False,
        "counts_generated": False,
        "build_or_submission_created": False,
    }
    return result, numeric


def evaluate_source(confirmation: str) -> None:
    protocol()
    require(confirmation == EVALUATE_CONFIRMATION, "evaluate confirmation mismatch")
    require(
        (PREPARED_DIR / "artifact-manifest.json").is_file(),
        "prepared manifest is missing",
    )
    bindings = validate_stage_authorization(
        "source-only-evaluate",
        confirmation,
        PREPARED_DIR / "artifact-manifest.json",
    )
    require(not EVALUATION_DIR.exists(), f"evaluation output exists: {EVALUATION_DIR}")
    with exclusive_lock(EVALUATE_LOCK), atomic_stage_directory(EVALUATION_DIR) as stage:
        with process_tree_memory_guard("evaluate") as memory_report:
            plan, arrays = validate_prepared()
            outer_results = []
            batch_names = sorted(set(plan["batches"]))
            for outer_batch in range(BATCH_COUNT):
                print(
                    f"[outer batch {outer_batch + 1}/{BATCH_COUNT}] "
                    f"{batch_names[outer_batch]}",
                    flush=True,
                )
                outer_results.append(
                    outer_batch_run(
                        outer_batch=outer_batch,
                        arrays=arrays,
                        plan=plan,
                    )
                )
            result, numeric = summarize_evaluation(outer_results, arrays, plan)
        result["input_bindings"] = bindings
        result["prepared_manifest_sha256"] = sha256_file(
            PREPARED_DIR / "artifact-manifest.json"
        )
        result["outer_folds"] = [
            {
                "batch": record["outer_batch"],
                "held_contexts": [
                    plan["contexts"][index] for index in record["held_contexts"]
                ],
                "fold_records": record["fold_records"],
            }
            for record in outer_results
        ]
        result["memory_report"] = "memory-report.json"
        atomic_json(stage / "memory-report.json", memory_report)
        validate_memory_report(memory_report)
        atomic_json(stage / "result.json", result)
        atomic_npz(stage / "metrics.npz", **numeric)
        atomic_json(stage / "artifact-manifest.json", directory_manifest(stage))
    print(f"wrote immutable source evaluation: {EVALUATION_DIR}")


def audit() -> None:
    value = protocol()
    manifest = read_json(DEPMAP_MANIFEST)
    controls = set(manifest["dataset_contract"]["control_targets"])
    experimental = tuple(
        target
        for target in manifest["dataset_contract"]["targets"]
        if target not in controls
    )
    evaluation_targets = permitted_evaluation_targets(experimental)
    folds = deterministic_target_folds(evaluation_targets)
    binding_paths = prepare_binding_paths()
    require(
        all(path.name != "delta.npy" for path in binding_paths),
        "prepare bindings include a complete donor delta",
    )
    source_text = RUN_PATH.read_text(encoding="utf-8")
    forbidden_command_present = any(
        ("add_parser(" + quote + name + quote) in source_text
        for name in ("develop", "submit")
        for quote in ('"', "'")
    )
    require(not forbidden_command_present, "forbidden CLI command exists")
    report = {
        "experiment": HERE.name,
        "protocol_status": value["status"],
        "prepare_authorization_exists": PREPARE_AUTHORIZATION.exists(),
        "evaluate_authorization_exists": EVALUATE_AUTHORIZATION.exists(),
        "prepared_exists": PREPARED_DIR.exists(),
        "evaluation_exists": EVALUATION_DIR.exists(),
        "write_root": str(HERE),
        "commands": ["audit", "prepare", "evaluate"],
        "forbidden_command_present": forbidden_command_present,
        "evaluation_target_count": len(evaluation_targets),
        "sealed_410_overlap_exclusions": list(SEALED_410_OVERLAP_TARGETS),
        "target_fold_sizes": [
            int(np.sum(folds == fold)) for fold in range(TARGET_FOLD_COUNT)
        ],
        "outer_fold_count": BATCH_COUNT * TARGET_FOLD_COUNT,
        "training_row_cap": TRAINING_ROW_CAP,
        "hgb_max_iter": HGB_MAX_ITER,
        "complete_delta_in_prepare_bindings": False,
        "source_target_sum": SOURCE_TARGET_SUM,
        "source_response_gene_count": EXPECTED_SOURCE_RESPONSE_GENE_COUNT,
        "source_response_axis_sha256": EXPECTED_SOURCE_RESPONSE_AXIS_SHA256,
        "fdr_axis_exclusions": list(FDR_AXIS_EXCLUSIONS),
        "fdr_missingness": {
            "available_context_target_pairs": EXPECTED_AVAILABLE_SOURCE_PAIR_COUNT,
            "nonfinite_entries_on_available_pairs": EXPECTED_FDR_NONFINITE_COUNT,
            "all_missing_available_blocks": EXPECTED_FDR_ALL_MISSING_BLOCK_COUNT,
            "all_missing_block_sha256": EXPECTED_FDR_ALL_MISSING_BLOCK_SHA256,
            "policy": "nonfinite_is_unlabelled_never_imputed",
            "unavailable_pairs_all_nan": True,
        },
        "single_donor_targets": list(SINGLE_DONOR_TARGETS),
        "one_donor_baseline_and_features_allowed": True,
        "intervention_minimum_donor_support": 2,
        "permutation_availability_pattern_stratified": True,
        "permutation_singleton_class_policy": "fixed_and_recorded",
        "memory_policy": {
            "max_process_tree_rss_gib": MAX_PROCESS_TREE_RSS_GIB,
            "min_system_free_percent": MIN_SYSTEM_FREE_PERCENT,
        },
    }
    if PREPARED_DIR.exists():
        validate_directory_manifest(PREPARED_DIR)
        report["prepared_manifest_valid"] = True
    if EVALUATION_DIR.exists():
        validate_directory_manifest(EVALUATION_DIR)
        report["evaluation_manifest_valid"] = True
    print(json.dumps(report, indent=2, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("audit", help="static protocol and immutable-artifact audit")
    prepare_parser = subparsers.add_parser(
        "prepare", help="authorized source-only preparation"
    )
    prepare_parser.add_argument("--confirm", required=True)
    evaluate_parser = subparsers.add_parser(
        "evaluate", help="authorized source-only cross-fitted evaluation"
    )
    evaluate_parser.add_argument("--confirm", required=True)
    args = parser.parse_args()
    if args.command == "audit":
        audit()
    elif args.command == "prepare":
        prepare_source(args.confirm)
    elif args.command == "evaluate":
        evaluate_source(args.confirm)
    else:
        raise AssertionError(args.command)


if __name__ == "__main__":
    main()

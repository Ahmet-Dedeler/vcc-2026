#!/usr/bin/env python3
"""Build Experiment 152's authorized source-only gene-consensus candidate.

The pure adapter functions in this file are intentionally testable without
opening AnnData. The executable build path refuses to proceed without a future
one-time source-build authorization. Once authorized, it may read only the
preflighted HCT non-targeting rows from the held harness.
"""

from __future__ import annotations

import gc
import hashlib
import importlib.util
import json
import os
import resource
import shutil
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PREREGISTRATION_PATH = HERE / "PREREGISTRATION.v1.json"
PREREGISTRATION_DIGEST_PATH = HERE / "PREREGISTRATION.v1.sha256"
BUILD_AUTHORIZATION_PATH = HERE / "SOURCE_BUILD_AUTHORIZATION.v1.json"
OUTPUT_DIR = HERE / "work/source-build-v1"
CONTROL = "non-targeting"
PERT_COL = "target_gene"
EXPECTED_CONTROL_START = 19_800
EXPECTED_CONTROL_STOP = 38_200
CONTROL_CHUNK_ROWS = 256
CPM_SUM = 1_000_000.0
GIB = 1024**3
RUNTIME_SOURCE_PATHS = (
    "experiments/149-dual-moment-hek293t/build_predictions.py",
    "experiments/147-dual-moment-emitter/build_candidates.py",
    "experiments/147-dual-moment-emitter/dual_moment.py",
    "scripts/make_submission_pm_cd4.py",
)
REQUIRED_BUILD_AUTHORIZATION_BINDINGS = frozenset(
    {
        "experiments/152-gene-consensus-ramp/PREREGISTRATION.v1.json",
        "experiments/152-gene-consensus-ramp/PREREGISTRATION.v1.sha256",
        "experiments/152-gene-consensus-ramp/test_experiment152.py",
        "experiments/149-dual-moment-hek293t/work/source-build-v1/effects.npz",
        "experiments/150-dual-moment-score-level-hct116/work/source-build-v1/effects.npz",
        "experiments/150-dual-moment-score-level-hct116/work/source-build-v1/build-report.json",
        "experiments/147-dual-moment-emitter/dual_moment.py",
        "experiments/147-dual-moment-emitter/build_candidates.py",
        "experiments/149-dual-moment-hek293t/build_predictions.py",
        "scripts/make_submission_pm_cd4.py",
        "data/harness/orion-hct116-thinned.h5ad",
        "data/harness/orion-hct116-thinned.h5ad.meta.json",
    }
)


@dataclass(frozen=True)
class GeneLambdaDiagnostics:
    eligible_gene_count: int
    agreeing_gene_count: int
    ramped_gene_count: int
    lambda_min: float
    lambda_max: float
    lambda_mean: float

    def to_dict(self) -> dict[str, int | float]:
        return asdict(self)


@dataclass(frozen=True)
class GeneRampDiagnostics:
    grand_total: int
    low_total: int
    high_total: int
    low_cells: int
    high_cells: int
    depth_contrast: float
    quasi_spread: float
    lambda_min: float
    lambda_max: float
    lambda_mean: float
    lambda_above_half_count: int
    lambda_one_count: int
    projected_low_gene_mass_l1: float
    clipped_gene_fraction: float
    desired_cpm_relative_l1: float
    feasible_cpm_relative_l1: float
    max_row_total: int
    output_nnz: int

    def to_dict(self) -> dict[str, int | float]:
        return asdict(self)


@dataclass(frozen=True)
class GeneRampPlan:
    primary_targets_hct_order: tuple[str, ...]
    hct_targets_order: tuple[str, ...]
    hct_only_targets: tuple[str, ...]
    lambda_by_target: dict[str, np.ndarray]
    per_target: dict[str, dict[str, Any]]
    hashes: dict[str, str]

    def summary(self) -> dict[str, Any]:
        return {
            "primary_targets_hct_order": self.primary_targets_hct_order,
            "hct_targets_order": self.hct_targets_order,
            "hct_only_targets": self.hct_only_targets,
            "per_target": self.per_target,
            "hashes": self.hashes,
        }


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_strings(values: tuple[str, ...] | list[str] | np.ndarray) -> str:
    payload = "\n".join(map(str, values)) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def array_sha256(value: np.ndarray) -> str:
    array = np.ascontiguousarray(value, dtype="<f8")
    return hashlib.sha256(array.tobytes()).hexdigest()


def stable_gene_ramp_seed(base: int, target: str) -> int:
    payload = f"{base}:{target}:gene-ramp-v1".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "little")


def peak_rss_gib() -> float:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if sys.platform == "darwin" else value * 1024 / GIB


def load_preregistration() -> dict[str, Any]:
    digest_line = PREREGISTRATION_DIGEST_PATH.read_text().strip().split()
    if len(digest_line) != 2 or digest_line[1] != PREREGISTRATION_PATH.name:
        raise ValueError("malformed immutable preregistration digest file")
    observed = sha256(PREREGISTRATION_PATH)
    if observed != digest_line[0]:
        raise ValueError(
            f"immutable preregistration hash changed: {observed} != {digest_line[0]}"
        )
    protocol = json.loads(PREREGISTRATION_PATH.read_text())
    if protocol.get("experiment") != "152-gene-consensus-ramp":
        raise ValueError("wrong preregistration experiment")
    if protocol.get("status") != "preregistered_before_hct_outcome_access":
        raise ValueError("Experiment 152 is not in its frozen preregistered state")
    return protocol


def _bound_path(relative: str) -> Path:
    path = (ROOT / relative).resolve()
    try:
        path.relative_to(ROOT.resolve())
    except ValueError as error:
        raise ValueError(f"authorization binding escapes repository: {relative}") from error
    return path


def verify_build_authorization(
    protocol: dict[str, Any],
    path: Path = BUILD_AUTHORIZATION_PATH,
) -> dict[str, Any]:
    if not path.is_file():
        raise PermissionError(
            "Experiment 152 source build requires a future one-time "
            "SOURCE_BUILD_AUTHORIZATION.v1.json"
        )
    authorization = json.loads(path.read_text())
    if authorization.get("experiment") != "152-gene-consensus-ramp":
        raise ValueError("wrong source-build authorization experiment")
    if authorization.get("status") != "authorized_once_source_only":
        raise PermissionError("Experiment 152 source build is not authorized once")
    builder = authorization.get("builder")
    if not isinstance(builder, dict):
        raise ValueError("source-build authorization has no builder binding")
    if builder.get("sha256") != sha256(Path(__file__).resolve()):
        raise ValueError("builder source does not match source-build authorization")
    bindings = authorization.get("bindings_sha256")
    if not isinstance(bindings, dict) or not bindings:
        raise ValueError("source-build authorization has no bindings")
    missing = sorted(REQUIRED_BUILD_AUTHORIZATION_BINDINGS - set(bindings))
    if missing:
        raise ValueError(f"source-build authorization misses bindings: {missing}")

    held_relative = protocol["source_artifacts"]["held_hct_harness"]["path"]
    held_expected = protocol["source_artifacts"]["held_hct_harness"]["sha256"]
    if bindings.get(held_relative) != held_expected:
        raise ValueError("source-build authorization has wrong held-HCT binding")
    if not _bound_path(held_relative).is_file():
        raise FileNotFoundError(_bound_path(held_relative))

    for relative, expected in bindings.items():
        if relative == held_relative:
            continue
        observed = sha256(_bound_path(relative))
        if observed != expected:
            raise ValueError(f"source-build authorization binding changed: {relative}")
    return authorization


def verify_source_artifacts(protocol: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Verify source files without scanning sealed held-perturbation bytes."""
    report: dict[str, dict[str, Any]] = {}
    for name, record in protocol["source_artifacts"].items():
        path = ROOT / record["path"]
        if not path.is_file():
            raise FileNotFoundError(path)
        if name == "held_hct_harness":
            report[name] = {
                "path": record["path"],
                "expected_sha256": record["sha256"],
                "verification": (
                    "authorization_binding_and_existence_only; expression bytes "
                    "not scanned by source verification"
                ),
            }
            continue
        observed = sha256(path)
        if observed != record["sha256"]:
            raise ValueError(
                f"{name} hash changed: {observed} != {record['sha256']}"
            )
        report[name] = {
            "path": record["path"],
            "sha256": observed,
            "verification": "sha256",
        }
    return report


def _load_effect_archive(
    path: Path,
    required_keys: tuple[str, ...],
) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        missing = sorted(set(required_keys) - set(archive.files))
        if missing:
            raise ValueError(f"{path}: missing required arrays {missing}")
        return {key: np.asarray(archive[key]) for key in archive.files}


def exact_midrank_percentiles(
    magnitudes: np.ndarray,
    eligible: np.ndarray,
) -> np.ndarray:
    """Return the frozen exact-tie midrank percentile on eligible coordinates."""
    values = np.asarray(magnitudes, dtype=np.float64)
    mask = np.asarray(eligible, dtype=bool)
    if values.ndim != 1 or mask.shape != values.shape:
        raise ValueError("magnitudes and eligibility must be same-shape vectors")
    if np.any(mask & ~np.isfinite(values)):
        raise ValueError("eligible magnitudes must be finite")
    positions = np.flatnonzero(mask)
    if not positions.size:
        raise ValueError("midrank percentile requires an eligible coordinate")
    ordered = positions[np.argsort(values[positions], kind="mergesort")]
    percentiles = np.zeros(values.size, dtype=np.float64)
    start = 0
    count = int(ordered.size)
    while start < count:
        stop = start + 1
        value = values[ordered[start]]
        while stop < count and values[ordered[stop]] == value:
            stop += 1
        group_size = stop - start
        percentile = (start + (group_size + 1) / 2.0) / count
        percentiles[ordered[start:stop]] = percentile
        start = stop
    return percentiles


def consensus_gene_lambdas(
    hek_relaxation: np.ndarray,
    hct_relaxation: np.ndarray,
    *,
    ineligible_lambda: float = 0.5,
) -> tuple[np.ndarray, GeneLambdaDiagnostics]:
    """Implement the frozen Experiment 152 gene-wise consensus formula exactly."""
    hek = np.asarray(hek_relaxation, dtype=np.float64)
    hct = np.asarray(hct_relaxation, dtype=np.float64)
    if hek.ndim != 1 or hct.shape != hek.shape:
        raise ValueError("HEK and HCT relaxation vectors must have equal shape")
    if ineligible_lambda != 0.5:
        raise ValueError("Experiment 152 freezes ineligible lambda at 0.5")

    magnitude = np.minimum(np.abs(hek), np.abs(hct))
    eligible = (
        np.isfinite(hek)
        & np.isfinite(hct)
        & np.isfinite(magnitude)
        & (magnitude > 0.0)
    )
    percentiles = exact_midrank_percentiles(magnitude, eligible)
    agreement = eligible & (hek * hct > 0.0)
    confidence = agreement.astype(np.float64) * np.clip(
        2.0 * percentiles - 1.0,
        0.0,
        1.0,
    )
    lambdas = np.full(hek.shape, ineligible_lambda, dtype=np.float64)
    lambdas[eligible] = 0.5 + 0.5 * confidence[eligible]
    if (
        not np.all(np.isfinite(lambdas))
        or np.any(lambdas < 0.5)
        or np.any(lambdas > 1.0)
    ):
        raise RuntimeError("gene-wise lambda construction left the frozen range")
    diagnostics = GeneLambdaDiagnostics(
        eligible_gene_count=int(eligible.sum()),
        agreeing_gene_count=int(agreement.sum()),
        ramped_gene_count=int(np.count_nonzero(lambdas > 0.5)),
        lambda_min=float(lambdas.min(initial=1.0)),
        lambda_max=float(lambdas.max(initial=0.0)),
        lambda_mean=float(lambdas.mean()),
    )
    return lambdas, diagnostics


def normalize_exact_cpm(values: np.ndarray) -> np.ndarray:
    """Normalize a finite nonnegative vector to an exact float64 CPM sum."""
    vector = np.asarray(values, dtype=np.float64)
    if (
        vector.ndim != 1
        or not np.all(np.isfinite(vector))
        or np.any(vector < 0.0)
    ):
        raise ValueError("desired arithmetic CPM must be finite and non-negative")
    total = float(vector.sum(dtype=np.float64))
    if total <= 0.0:
        raise ValueError("desired arithmetic CPM must have positive mass")
    normalized = CPM_SUM * vector / total
    pivot = int(np.argmax(normalized))
    for _ in range(4):
        residual = CPM_SUM - float(normalized.sum(dtype=np.float64))
        if residual == 0.0:
            break
        normalized[pivot] += residual
    if float(normalized.sum(dtype=np.float64)) != CPM_SUM:
        raise RuntimeError("failed exact float64 CPM normalization")
    if normalized[pivot] < 0.0 or not np.all(np.isfinite(normalized)):
        raise RuntimeError("exact CPM normalization produced an invalid coordinate")
    return normalized


def desired_gene_ramp_cpm(
    incumbent_cpm: np.ndarray,
    direction_cpm: np.ndarray,
    lambda_gene: np.ndarray,
) -> np.ndarray:
    """Apply the frozen coordinate interpolation followed by exact renormalization."""
    incumbent = np.asarray(incumbent_cpm, dtype=np.float64)
    direction = np.asarray(direction_cpm, dtype=np.float64)
    lambdas = np.asarray(lambda_gene, dtype=np.float64)
    if (
        incumbent.ndim != 1
        or direction.shape != incumbent.shape
        or lambdas.shape != incumbent.shape
    ):
        raise ValueError("CPM profiles and lambda_gene must be same-shape vectors")
    if (
        not np.all(np.isfinite(incumbent))
        or not np.all(np.isfinite(direction))
        or not np.all(np.isfinite(lambdas))
        or np.any(incumbent < 0.0)
        or np.any(direction < 0.0)
        or np.any(lambdas < 0.5)
        or np.any(lambdas > 1.0)
    ):
        raise ValueError("invalid Experiment 152 interpolation input")
    interpolated = incumbent + lambdas * (direction - incumbent)
    return normalize_exact_cpm(interpolated)


def emit_gene_ramp(
    incumbent_counts,
    direction_counts,
    lambda_gene: np.ndarray,
    *,
    dual_moment: ModuleType,
    depth_contrast: float = 0.5,
    quasi_spread: float = 0.25,
    seed: int,
) -> tuple[sp.csr_matrix, GeneRampDiagnostics]:
    """Use frozen exact-margin primitives with a gene-wise desired CPM vector."""
    incumbent = dual_moment.as_integer_csr(incumbent_counts)
    direction = dual_moment.as_integer_csr(direction_counts)
    if incumbent.shape != direction.shape:
        raise ValueError("incumbent and direction blocks must have identical shape")
    n_cells, n_genes = incumbent.shape
    lambdas = np.asarray(lambda_gene, dtype=np.float64)
    if lambdas.shape != (n_genes,):
        raise ValueError("lambda_gene does not match the emitted gene axis")
    if n_cells % 2:
        raise ValueError("gene-ramp emitter requires an even cell count")

    columns = np.asarray(incumbent.sum(axis=0), dtype=np.int64).ravel()
    grand_total = int(columns.sum())
    rows, low_mask = dual_moment.two_depth_row_margins(
        grand_total,
        n_cells,
        depth_contrast,
        seed=seed + 1,
    )
    incumbent_cpm = dual_moment.arithmetic_mean_cpm(incumbent)
    direction_cpm = dual_moment.arithmetic_mean_cpm(direction)
    desired_cpm = desired_gene_ramp_cpm(
        incumbent_cpm,
        direction_cpm,
        lambdas,
    )
    raw_low, low_columns, feasible_cpm = dual_moment.solve_low_group_totals(
        columns,
        desired_cpm,
        rows,
        low_mask,
        seed=seed + 2,
    )
    allocate_group = getattr(dual_moment, "_allocate_group", None)
    if not callable(allocate_group):
        raise RuntimeError("frozen dual_moment._allocate_group is unavailable")
    high_columns = columns - low_columns
    low_block = allocate_group(
        low_columns,
        rows[low_mask],
        quasi_spread=quasi_spread,
        seed=seed + 3,
    )
    high_block = allocate_group(
        high_columns,
        rows[~low_mask],
        quasi_spread=quasi_spread,
        seed=seed + 4,
    )
    emitted = sp.csr_matrix(np.vstack((low_block, high_block)))
    emitted.eliminate_zeros()
    emitted.sort_indices()

    got_rows = np.asarray(emitted.sum(axis=1), dtype=np.int64).ravel()
    got_columns = np.asarray(emitted.sum(axis=0), dtype=np.int64).ravel()
    if not np.array_equal(got_rows, rows):
        raise RuntimeError("final row margins are not exact")
    if not np.array_equal(got_columns, columns):
        raise RuntimeError("final column margins are not exact")
    if not np.array_equal(
        dual_moment.scorer_pseudobulk(emitted),
        dual_moment.scorer_pseudobulk(incumbent),
    ):
        raise RuntimeError("final scorer pseudobulk changed")
    if emitted.data.size and (
        np.any(emitted.data < 0)
        or not np.array_equal(emitted.data, np.rint(emitted.data))
    ):
        raise RuntimeError("final count matrix is not non-negative integer")
    if int(got_rows.max(initial=0)) > dual_moment.MAX_COUNTS_PER_CELL:
        raise RuntimeError("final row exceeds the VCC count cap")

    achieved_cpm = dual_moment.arithmetic_mean_cpm(emitted)
    diagnostics = GeneRampDiagnostics(
        grand_total=grand_total,
        low_total=int(rows[low_mask].sum()),
        high_total=int(rows[~low_mask].sum()),
        low_cells=int(low_mask.sum()),
        high_cells=int((~low_mask).sum()),
        depth_contrast=float(depth_contrast),
        quasi_spread=float(quasi_spread),
        lambda_min=float(lambdas.min(initial=1.0)),
        lambda_max=float(lambdas.max(initial=0.0)),
        lambda_mean=float(lambdas.mean()),
        lambda_above_half_count=int(np.count_nonzero(lambdas > 0.5)),
        lambda_one_count=int(np.count_nonzero(lambdas == 1.0)),
        projected_low_gene_mass_l1=float(
            np.abs(low_columns - raw_low).sum()
            / max(1.0, float(np.abs(raw_low).sum()))
        ),
        clipped_gene_fraction=float(
            np.mean((raw_low < 0.0) | (raw_low > columns))
        ),
        desired_cpm_relative_l1=float(
            np.abs(achieved_cpm - desired_cpm).sum()
            / float(np.abs(desired_cpm).sum())
        ),
        feasible_cpm_relative_l1=float(
            np.abs(achieved_cpm - feasible_cpm).sum()
            / float(np.abs(feasible_cpm).sum())
        ),
        max_row_total=int(got_rows.max(initial=0)),
        output_nnz=int(emitted.nnz),
    )
    return emitted, diagnostics


def compute_gene_ramp_plan(
    protocol: dict[str, Any],
    hek: dict[str, np.ndarray],
    hct: dict[str, np.ndarray],
) -> GeneRampPlan:
    genes = np.asarray(hct["genes"]).astype(str)
    if not np.array_equal(np.asarray(hek["genes"]).astype(str), genes):
        raise ValueError("HEK and HCT frozen effect gene axes differ")
    hek_targets = np.asarray(hek["targets"]).astype(str)
    hct_targets = np.asarray(hct["targets"]).astype(str)
    hek_index = {target: index for index, target in enumerate(hek_targets)}
    hct_index = {target: index for index, target in enumerate(hct_targets)}
    common = set(hek_index) & set(hct_index)
    surface = protocol["common_surface"]
    if len(common) != int(surface["common_target_count"]):
        raise ValueError("common target count changed")
    primary = tuple(target for target in hct_targets if target in common)
    primary_hash = sha256_strings(primary)
    if (
        primary_hash
        != surface["primary_targets_in_hct_order_sha256_newline_terminated"]
    ):
        raise ValueError("primary HCT target order changed")
    hct_only = tuple(sorted(set(hct_index) - set(hek_index)))
    expected_hct_only = tuple(
        sorted(surface["hct_only_targets_excluded_from_primary_analysis"])
    )
    if hct_only != expected_hct_only:
        raise ValueError("HCT-only target set changed")

    hek_relaxation = (
        np.asarray(hek["direction"], dtype=np.float64)
        - np.asarray(hek["incumbent"], dtype=np.float64)
    )
    hct_relaxation = (
        np.asarray(hct["direction"], dtype=np.float64)
        - np.asarray(hct["incumbent"], dtype=np.float64)
    )
    expected_eligible = int(surface["eligible_gene_count_per_target"])
    lambda_by_target: dict[str, np.ndarray] = {}
    per_target: dict[str, dict[str, Any]] = {}
    for target_value in hct_targets:
        target = str(target_value)
        if target in common:
            lambdas, diagnostics = consensus_gene_lambdas(
                hek_relaxation[hek_index[target]],
                hct_relaxation[hct_index[target]],
                ineligible_lambda=float(
                    protocol["gene_wise_selector"]["ineligible_lambda"]
                ),
            )
            if diagnostics.eligible_gene_count != expected_eligible:
                raise ValueError(
                    f"{target}: eligible gene count "
                    f"{diagnostics.eligible_gene_count} != {expected_eligible}"
                )
            group = "common_gene_ramp"
        else:
            lambdas = np.full(
                len(genes),
                float(protocol["gene_wise_selector"]["hct_only_lambda"]),
                dtype=np.float64,
            )
            diagnostics = GeneLambdaDiagnostics(
                eligible_gene_count=0,
                agreeing_gene_count=0,
                ramped_gene_count=len(genes),
                lambda_min=1.0,
                lambda_max=1.0,
                lambda_mean=1.0,
            )
            group = "hct_only_uniform_one"
        lambda_by_target[target] = lambdas
        per_target[target] = {
            "selector_group": group,
            "lambda_sha256_float64_le": array_sha256(lambdas),
            **diagnostics.to_dict(),
        }
    return GeneRampPlan(
        primary_targets_hct_order=primary,
        hct_targets_order=tuple(map(str, hct_targets)),
        hct_only_targets=hct_only,
        lambda_by_target=lambda_by_target,
        per_target=per_target,
        hashes={
            "primary_targets_in_hct_order_sha256_newline_terminated": primary_hash,
            "hct_targets_in_order_sha256_newline_terminated": sha256_strings(
                hct_targets
            ),
            "lambda_matrix_sha256_float64_le": array_sha256(
                np.vstack([lambda_by_target[str(target)] for target in hct_targets])
            ),
        },
    )


def load_plan_inputs(
    protocol: dict[str, Any],
) -> tuple[GeneRampPlan, dict[str, np.ndarray]]:
    required = tuple(protocol["required_effect_keys"])
    hek = _load_effect_archive(
        ROOT / protocol["source_artifacts"]["hek_effects"]["path"],
        required,
    )
    hct = _load_effect_archive(
        ROOT / protocol["source_artifacts"]["hct_effects"]["path"],
        required + ("bctrl",),
    )
    return compute_gene_ramp_plan(protocol, hek, hct), hct


def _load_build_runtime() -> ModuleType:
    return load_module(
        "exp152_build_base",
        ROOT / "experiments/149-dual-moment-hek293t/build_predictions.py",
    )


def runtime_source_hashes() -> dict[str, str]:
    return {relative: sha256(ROOT / relative) for relative in RUNTIME_SOURCE_PATHS}


def build_candidate(
    protocol: dict[str, Any],
    plan: GeneRampPlan,
    hct: dict[str, np.ndarray],
    staging: Path,
) -> dict[str, Any]:
    import anndata as ad
    import h5py

    base = _load_build_runtime()
    dual_moment = base.DM
    stream = base.BUILD
    real_path = ROOT / protocol["source_artifacts"]["held_hct_harness"]["path"]
    output_path = staging / "candidate.h5ad"
    metadata_path = staging / f".metadata-{os.getpid()}.h5ad"

    genes = np.asarray(hct["genes"]).astype(str)
    targets = np.asarray(hct["targets"]).astype(str)
    incumbent_effects = np.asarray(hct["incumbent"], dtype=np.float64)
    direction_effects = np.asarray(hct["direction"], dtype=np.float64)
    bctrl = np.asarray(hct["bctrl"], dtype=np.float64)
    fixed = protocol["fixed_emission"]
    if incumbent_effects.shape != direction_effects.shape:
        raise ValueError("HCT incumbent and direction effect shapes differ")
    if incumbent_effects.shape != (len(targets), len(genes)):
        raise ValueError("HCT effect arrays do not match their frozen axes")
    if int(fixed["cells_per_target"]) != base.CELLS_PER_TARGET:
        raise ValueError("Experiment 152 cell count differs from Experiment 150")

    real = ad.read_h5ad(real_path, backed="r")
    temporary_output = output_path.with_name(
        f".{output_path.name}.tmp-{os.getpid()}"
    )
    try:
        labels = real.obs[PERT_COL].astype(str).to_numpy()
        control_indices = np.flatnonzero(labels == CONTROL)
        expected_indices = np.arange(
            EXPECTED_CONTROL_START,
            EXPECTED_CONTROL_STOP,
            dtype=np.int64,
        )
        if not np.array_equal(control_indices, expected_indices):
            raise ValueError("HCT control-row preflight contract changed")
        if not np.array_equal(np.asarray(real.var_names).astype(str), genes):
            raise ValueError("HCT harness and frozen effect gene axes differ")
        if set(labels) - {CONTROL} != set(targets):
            raise ValueError("HCT metadata target set differs from frozen effects")

        libraries = base.control_libraries(
            real,
            EXPECTED_CONTROL_START,
            EXPECTED_CONTROL_STOP,
            chunk_rows=CONTROL_CHUNK_ROWS,
        )
        control_sums = np.zeros(len(genes), dtype=np.float64)
        for start in range(
            EXPECTED_CONTROL_START,
            EXPECTED_CONTROL_STOP,
            CONTROL_CHUNK_ROWS,
        ):
            stop = min(start + CONTROL_CHUNK_ROWS, EXPECTED_CONTROL_STOP)
            control_sums += np.asarray(
                stream.read_rows(real, slice(start, stop)).sum(axis=0),
                dtype=np.float64,
            ).ravel()
        observed_bctrl = np.log1p(50_000.0 * control_sums / control_sums.sum())
        if not np.array_equal(
            observed_bctrl.astype(np.float32),
            bctrl.astype(np.float32),
        ):
            raise ValueError("HCT controls differ from the frozen control profile")

        stream.write_metadata_template(
            metadata_path,
            target_order=list(map(str, targets)),
            cells_per_target=base.CELLS_PER_TARGET,
            control_cells=len(control_indices),
            var=real.var,
        )
        inc_pick_rng, inc_emission_rng = base.PROD.make_rngs(0)
        dir_pick_rng, dir_emission_rng = base.PROD.make_rngs(0)
        per_target: dict[str, dict[str, Any]] = {}
        with (
            h5py.File(metadata_path, "r") as metadata_h5,
            h5py.File(temporary_output, "w") as candidate_h5,
        ):
            stream.clone_non_x_metadata(metadata_h5, candidate_h5)
            n_rows = len(targets) * base.CELLS_PER_TARGET + len(control_indices)
            writer = stream.CSRStreamWriter(
                candidate_h5,
                n_rows=n_rows,
                n_cols=len(genes),
                data_dtype=np.dtype("float32"),
            )
            for index, target_value in enumerate(targets):
                target = str(target_value)
                incumbent = base.emitted_block(
                    incumbent_effects[index],
                    bctrl,
                    libraries,
                    inc_pick_rng,
                    inc_emission_rng,
                )
                direction = base.emitted_block(
                    direction_effects[index],
                    bctrl,
                    libraries,
                    dir_pick_rng,
                    dir_emission_rng,
                )
                lambdas = plan.lambda_by_target[target]
                candidate, diagnostics = emit_gene_ramp(
                    incumbent,
                    direction,
                    lambdas,
                    dual_moment=dual_moment,
                    depth_contrast=float(fixed["depth_contrast"]),
                    quasi_spread=float(fixed["within_group_quasi_spread"]),
                    seed=stable_gene_ramp_seed(
                        int(fixed["dual_moment_seed"]),
                        target,
                    ),
                )
                incumbent_columns = np.asarray(
                    incumbent.sum(axis=0), dtype=np.int64
                ).ravel()
                candidate_columns = np.asarray(
                    candidate.sum(axis=0), dtype=np.int64
                ).ravel()
                if not np.array_equal(incumbent_columns, candidate_columns):
                    raise RuntimeError(f"{target}: candidate column totals changed")
                if not np.array_equal(
                    dual_moment.scorer_pseudobulk(incumbent),
                    dual_moment.scorer_pseudobulk(candidate),
                ):
                    raise RuntimeError(f"{target}: candidate pseudobulk changed")
                row_totals = np.asarray(
                    candidate.sum(axis=1), dtype=np.int64
                ).ravel()
                if (
                    np.any(row_totals <= 0)
                    or row_totals.max(initial=0)
                    > dual_moment.MAX_COUNTS_PER_CELL
                ):
                    raise RuntimeError(f"{target}: invalid candidate row total")
                writer.append(candidate.astype(np.float32))
                per_target[target] = {
                    **plan.per_target[target],
                    "columns_exact": True,
                    "pseudobulk_exact": True,
                    "incumbent_columns_sha256": hashlib.sha256(
                        np.ascontiguousarray(
                            incumbent_columns,
                            dtype="<i8",
                        ).tobytes()
                    ).hexdigest(),
                    "candidate_nnz": int(candidate.nnz),
                    "incumbent_nnz": int(incumbent.nnz),
                    "max_row_total": int(row_totals.max(initial=0)),
                    "gene_ramp": diagnostics.to_dict(),
                }
                if (index + 1) % 10 == 0:
                    print(f"[emit] {index + 1}/{len(targets)} targets", flush=True)
                del incumbent, direction, candidate
                gc.collect()

            for start in range(
                EXPECTED_CONTROL_START,
                EXPECTED_CONTROL_STOP,
                CONTROL_CHUNK_ROWS,
            ):
                stop = min(start + CONTROL_CHUNK_ROWS, EXPECTED_CONTROL_STOP)
                control_block = stream.read_rows(real, slice(start, stop))
                writer.append(control_block.astype(np.float32))
                del control_block
            writer.finish()
            candidate_h5.flush()
            output_nnz = writer.nnz

        os.replace(temporary_output, output_path)
        metadata_path.unlink(missing_ok=True)
        built = ad.read_h5ad(output_path, backed="r")
        try:
            built_labels = built.obs[PERT_COL].astype(str).to_numpy()
            if built.n_vars != len(genes):
                raise RuntimeError("completed candidate gene count differs")
            if not np.array_equal(np.asarray(built.var_names).astype(str), genes):
                raise RuntimeError("completed candidate gene axis differs")
            if np.count_nonzero(built_labels == CONTROL) != len(control_indices):
                raise RuntimeError("completed candidate control count differs")
            for target in targets:
                if (
                    np.count_nonzero(built_labels == str(target))
                    != base.CELLS_PER_TARGET
                ):
                    raise RuntimeError(f"{target}: completed candidate row count differs")
        finally:
            built.file.close()

        return {
            "candidate": str((OUTPUT_DIR / "candidate.h5ad").relative_to(ROOT)),
            "candidate_sha256": sha256(output_path),
            "candidate_nnz": int(output_nnz),
            "target_order_sha256_newline_terminated": sha256_strings(targets),
            "targets": len(targets),
            "controls": len(control_indices),
            "real_control_rows_read": len(control_indices),
            "real_perturbation_rows_read": 0,
            "all_column_totals_exact": all(
                record["columns_exact"] for record in per_target.values()
            ),
            "all_pseudobulks_exact": all(
                record["pseudobulk_exact"] for record in per_target.values()
            ),
            "per_target": per_target,
        }
    finally:
        real.file.close()
        metadata_path.unlink(missing_ok=True)
        temporary_output.unlink(missing_ok=True)


def main() -> int:
    started = time.monotonic()
    protocol = load_preregistration()
    authorization = verify_build_authorization(protocol)
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"immutable build output already exists: {OUTPUT_DIR}")
    source_verification = verify_source_artifacts(protocol)
    plan, hct = load_plan_inputs(protocol)
    staging = OUTPUT_DIR.with_name(f".{OUTPUT_DIR.name}.tmp-{os.getpid()}")
    if staging.exists():
        raise FileExistsError(f"staging directory already exists: {staging}")
    staging.mkdir(parents=True)
    try:
        (staging / "gene-ramp.json").write_text(
            json.dumps(plan.summary(), indent=2, sort_keys=True) + "\n"
        )
        prediction = build_candidate(protocol, plan, hct, staging)
        report = {
            "experiment": "152-gene-consensus-ramp",
            "stage": "authorized_source_only_build",
            "source_build_authorization_sha256": sha256(BUILD_AUTHORIZATION_PATH),
            "preregistration_sha256": sha256(PREREGISTRATION_PATH),
            "preregistration_digest_file_sha256": sha256(
                PREREGISTRATION_DIGEST_PATH
            ),
            "builder_sha256": sha256(Path(__file__).resolve()),
            "runtime_source_hashes": runtime_source_hashes(),
            "source_verification": source_verification,
            "gene_ramp": plan.summary(),
            "prediction": prediction,
            "peak_rss_gib": peak_rss_gib(),
            "elapsed_seconds": time.monotonic() - started,
            "held_hct_expression_access": {
                "control_rows_read": prediction["real_control_rows_read"],
                "perturbation_rows_read": 0,
            },
            "cell_eval2_invoked": False,
            "challenge_build_created": False,
            "challenge_submission_created": False,
            "authorization_status": authorization["status"],
        }
        (staging / "build-report.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n"
        )
        os.replace(staging, OUTPUT_DIR)
        print(
            json.dumps(
                {
                    "candidate_sha256": prediction["candidate_sha256"],
                    "lambda_matrix_sha256": plan.hashes[
                        "lambda_matrix_sha256_float64_le"
                    ],
                    "all_column_totals_exact": prediction[
                        "all_column_totals_exact"
                    ],
                    "all_pseudobulks_exact": prediction[
                        "all_pseudobulks_exact"
                    ],
                    "real_perturbation_rows_read": 0,
                },
                indent=2,
            )
        )
        return 0
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Build Experiment 107's source-only strict validation artifacts.

This builder may read the Nadig target and gene axes, but it never opens or
hashes ``data/donors/nadig_jurkat/delta.npy``.  The expensive build command is
explicitly gated because each real whitening eigendecomposition is large.
"""

from __future__ import annotations

import argparse
import csv
import gc
import importlib.util
import json
import os
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any, Mapping, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.dont_write_bytecode = True
STRICT_PATH = HERE / "strict_nadig_validation.py"
PRODUCTION_PATH = ROOT / "scripts/make_submission_pm_cd4.py"
WHITENER_BUILDER_PATH = ROOT / "scripts/build_union_whitener.py"
OUTPUT_ROOT = HERE / "work/strict-validation"
TOP_LEVEL_MANIFEST = HERE / "work/strict-validation-inputs.json"

DONOR_PATHS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
PERTURBMATCH = ROOT / "data/external/perturbmatch"
CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
ESSENTIAL_PATH = ROOT / "data/donors/k562_essential/perts.txt"
CHALLENGE = ROOT / "data/challenge"
SOURCE_DIRECTORY = HERE / "work/frozen-source-model"
FORBIDDEN_NADIG_OUTCOME = ROOT / "data/donors/nadig_jurkat/delta.npy"
REPAIR_PROVENANCE_PATHS = (
    HERE / "BTBD11_OPERATOR_SUPPORT_REPAIR.json",
    HERE / "BTBD11_VALIDATION_PROVENANCE_AMENDMENT.json",
    HERE / "reseal_btbd11_repair.py",
    HERE / "finalize_btbd11_provenance.py",
    HERE / "work/btbd11-reseal-receipt.json",
    HERE / "work/btbd11-parent-binding-reconciliation-receipt.json",
)

BASE_PROVENANCE_PATHS = (
    HERE / "strict_validation_contract.json",
    HERE / "preregistered_protocol.json",
    HERE / "strict_artifact_builder.py",
    STRICT_PATH,
    *REPAIR_PROVENANCE_PATHS,
)
SOURCE_PROVENANCE_PATHS = (
    SOURCE_DIRECTORY / "source-model.npz",
    SOURCE_DIRECTORY / "manifest.json",
    HERE / "model_input.py",
    HERE / "modeling.py",
    HERE / "fit_source.py",
    HERE / "source_artifacts.py",
)
PRODUCTION_IMPLEMENTATION_PATHS = (
    ROOT / "scripts/make_submission_pm_cd4.py",
    ROOT / "scripts/build_submission.py",
    ROOT / "scripts/make_submission.py",
    ROOT / "src/vcc2026/donor_transfer.py",
    ROOT / "src/vcc2026/multi_donor.py",
    ROOT / "src/vcc2026/whiten.py",
    ROOT / "experiments/061-perturbmatch-hct-count/build_candidate.py",
    ROOT / "experiments/062-cd4-post-whiten/run.py",
)
DONOR_FILENAMES = ("delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt", "meta.json")
PERTURBMATCH_FILENAMES = ("delta.npy", "genes.txt", "perts.txt", "meta.json")
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
CD4_FILENAMES = ("delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt", "meta.json")


def _load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


STRICT = _load_module("exp107_strict_artifacts_validation", STRICT_PATH)
_PRODUCTION_MODULE: ModuleType | None = None


@dataclass(frozen=True)
class SurfacePlan:
    key: str
    targets: tuple[str, ...]
    score_genes: tuple[str, ...]
    score_indices_in_source: np.ndarray
    score_indices_in_official: np.ndarray
    operator_columns: np.ndarray
    basis_rows: Mapping[str, np.ndarray]
    donor_order: tuple[str, ...]
    source_target_indices: np.ndarray
    source_transferable: np.ndarray
    held_gene_axis: tuple[str, ...] | None = None
    harness_path: Path | None = None

    @property
    def role_suffix(self) -> str:
        return self.key


def _read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: missing column {column}")
        values = tuple(str(row[column]) for row in reader)
    if not values or any(not value for value in values):
        raise ValueError(f"{path}: invalid {column} axis")
    return values


def _read_held_metadata(path: Path) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Read only H5AD obs/var metadata; never materialize the expression matrix."""

    import anndata as ad

    handle = ad.read_h5ad(path, backed="r")
    try:
        genes = tuple(np.asarray(handle.var_names).astype(str))
        labels = np.asarray(handle.obs["target_gene"].astype(str))
        targets = tuple(sorted(set(labels) - {"non-targeting"}))
    finally:
        handle.file.close()
    return genes, targets


def _verify_surface_plan(plan: SurfacePlan, frozen: Mapping[str, Any]) -> None:
    if len(plan.targets) != int(frozen["target_count"]):
        raise ValueError(f"{plan.key}: target count changed")
    if STRICT.sha256_strings(plan.targets) != frozen["target_axis_sha256"]:
        raise ValueError(f"{plan.key}: target axis changed")
    if len(plan.score_genes) != int(frozen["score_gene_count"]):
        raise ValueError(f"{plan.key}: score-gene count changed")
    if STRICT.sha256_strings(plan.score_genes) != frozen["score_axis_sha256"]:
        raise ValueError(f"{plan.key}: score axis changed")
    if len(plan.operator_columns) != int(frozen["whitener_column_count"]):
        raise ValueError(f"{plan.key}: operator-column count changed")
    if (
        STRICT.sha256_int64(plan.operator_columns)
        != frozen["whitener_columns_int64_sha256"]
    ):
        raise ValueError(f"{plan.key}: operator columns changed")
    if tuple(plan.donor_order) != tuple(frozen["operator_donors"]):
        raise ValueError(f"{plan.key}: operator donor order changed")
    expected_rows = frozen["whitener_basis_rows_by_donor"]
    if set(plan.basis_rows) != set(plan.donor_order):
        raise ValueError(f"{plan.key}: basis donor set changed")
    for donor in plan.donor_order:
        rows = np.asarray(plan.basis_rows[donor], dtype=np.int64)
        if len(rows) != int(expected_rows[donor]["count"]):
            raise ValueError(f"{plan.key}: {donor} basis-row count changed")
        if STRICT.sha256_int64(rows) != expected_rows[donor]["int64_sha256"]:
            raise ValueError(f"{plan.key}: {donor} basis-row axis changed")
    if int(plan.source_transferable.sum()) != int(frozen["transferable_target_count"]):
        raise ValueError(f"{plan.key}: transferable target count changed")


def build_surface_plans() -> tuple[dict[str, SurfacePlan], Mapping[str, Any], Any]:
    """Reconstruct all frozen axes without touching the Nadig outcome matrix."""

    contract = STRICT.validate_contract()
    source = STRICT.load_frozen_source_model()
    official = _read_csv_column(CHALLENGE / "gene_names.csv", "gene_name")
    vcc_targets = _read_csv_column(CHALLENGE / "pert_counts.csv", "target_gene")
    official_index = {symbol: index for index, symbol in enumerate(official)}
    source_target_index = {
        symbol: index for index, symbol in enumerate(source.target_symbols)
    }
    source_response_index = {
        symbol: index for index, symbol in enumerate(source.response_symbols)
    }
    donor_genes = {
        name: STRICT._read_axis(  # noqa: SLF001 - shared strict implementation
            path / "genes.txt", f"{name} gene axis", require_unique=False
        )
        for name, path in DONOR_PATHS.items()
    }
    donor_perts = {
        name: STRICT._read_axis(  # noqa: SLF001 - shared strict implementation
            path / "perts.txt", f"{name} perturbation axis", require_unique=False
        )
        for name, path in DONOR_PATHS.items()
    }
    essential = set(STRICT._read_axis(ESSENTIAL_PATH, "essential target axis"))  # noqa: SLF001

    nadig_targets = STRICT._read_axis(  # noqa: SLF001
        ROOT / "data/donors/nadig_jurkat/perts.txt", "Nadig target axis"
    )
    nadig_genes = STRICT._read_axis(  # noqa: SLF001
        ROOT / "data/donors/nadig_jurkat/genes.txt", "Nadig gene axis"
    )
    gallery_axes = STRICT.reconstruct_gallery_axes(
        source,
        nadig_targets=nadig_targets,
        nadig_genes=nadig_genes,
        official_genes=official,
        vcc_targets=vcc_targets,
        donor_genes=donor_genes,
        donor_perts=donor_perts,
        essential_targets=essential,
        contract=contract,
    )
    plans: dict[str, SurfacePlan] = {}
    for axes in gallery_axes:
        key = f"nadig.gallery_{axes.number}"
        plan = SurfacePlan(
            key=key,
            targets=axes.targets,
            score_genes=axes.score_genes,
            score_indices_in_source=axes.score_indices_in_source,
            score_indices_in_official=axes.score_indices_in_official,
            operator_columns=axes.whitener_columns,
            basis_rows=axes.basis_rows,
            donor_order=tuple(STRICT.DONOR_ORDER),
            source_target_indices=axes.target_indices_in_source,
            source_transferable=axes.source_transferable,
        )
        _verify_surface_plan(
            plan,
            STRICT._surface_contract(  # noqa: SLF001
                f"operator.{key}", contract
            ),
        )
        plans[key] = plan

    combined_targets = tuple(
        target
        for number in range(1, 7)
        for target in plans[f"nadig.gallery_{number}"].targets
    )
    combined_set = set(combined_targets)
    perturbmatch_targets = set(
        STRICT._read_axis(  # noqa: SLF001
            PERTURBMATCH / "perts.txt", "PerturbMatch target axis"
        )
    )
    k562_targets = set(donor_perts["k562_gwps"])
    cd4_target_sets = [
        set(
            STRICT._read_axis(  # noqa: SLF001
                CD4_ROOT / condition / "perts.txt", f"{condition} target axis"
            )
        )
        for condition in CD4_NAMES
    ]
    known = contract["known_coverage"]
    if len(combined_set & perturbmatch_targets) != int(known["perturbmatch_targets"]):
        raise ValueError("PerturbMatch gallery coverage changed")
    if len(combined_set & perturbmatch_targets & k562_targets) != int(
        known["k562_perturbmatch_replacement_targets"]
    ):
        raise ValueError("K562/PerturbMatch replacement coverage changed")
    if len(combined_set & set.intersection(*cd4_target_sets)) != int(
        known["cd4_targets"]
    ):
        raise ValueError("complete CD4 gallery coverage changed")
    fallback_targets = tuple(
        target
        for target in combined_targets
        if not any(target in set(donor_perts[name]) for name in STRICT.DONOR_ORDER)
    )
    if fallback_targets != tuple(known["incumbent_mean_fallback_targets"]):
        raise ValueError(f"incumbent mean-fallback targets changed: {fallback_targets}")

    for surface in STRICT.AUTHORIZING_HELD_SURFACES:
        frozen = contract["held_context_gate"]["surfaces"][surface]
        harness = ROOT / frozen["harness_path"]
        held_genes, targets = _read_held_metadata(harness)
        if len(held_genes) != int(frozen["held_gene_count"]):
            raise ValueError(f"{surface}: held gene count changed")
        if STRICT.sha256_strings(held_genes) != frozen["held_gene_axis_sha256"]:
            raise ValueError(f"{surface}: held gene axis changed")
        donor_order = tuple(frozen["operator_donors"])
        score_genes = STRICT.score_axis(
            source.response_symbols,
            source.response_exportable,
            held_genes,
            vcc_targets,
            targets,
        )
        operator_columns = STRICT.whitener_columns(
            official,
            donor_genes,
            vcc_targets,
            targets,
            donor_order=donor_order,
        )
        basis_rows = STRICT.select_whitener_basis_rows(
            donor_perts,
            essential,
            set(targets),
            donor_order=donor_order,
        )
        target_indices = np.asarray(
            [source_target_index.get(symbol, -1) for symbol in targets], dtype=np.int64
        )
        transferable = np.asarray(
            [
                index >= 0 and bool(source.transferable[index])
                for index in target_indices
            ],
            dtype=bool,
        )
        plan = SurfacePlan(
            key=surface,
            targets=targets,
            score_genes=score_genes,
            score_indices_in_source=np.asarray(
                [source_response_index[symbol] for symbol in score_genes],
                dtype=np.int64,
            ),
            score_indices_in_official=np.asarray(
                [official_index[symbol] for symbol in score_genes], dtype=np.int64
            ),
            operator_columns=operator_columns,
            basis_rows=basis_rows,
            donor_order=donor_order,
            source_target_indices=target_indices,
            source_transferable=transferable,
            held_gene_axis=held_genes,
            harness_path=harness,
        )
        _verify_surface_plan(plan, frozen)
        plans[surface] = plan
    return plans, contract, source


def role_paths(role: str, output_root: Path = OUTPUT_ROOT) -> tuple[Path, Path]:
    return STRICT.expected_role_bundle_paths(role, bundle_root=output_root)


def existing_role_bundle(
    role: str,
    contract: Mapping[str, Any],
    *,
    output_root: Path = OUTPUT_ROOT,
) -> dict[str, Any] | None:
    paths = role_paths(role, output_root)
    present = tuple(path.is_file() and not path.is_symlink() for path in paths)
    if not any(present):
        return None
    if not all(present):
        raise RuntimeError(f"incomplete role bundle requires manual audit: {role}")
    return STRICT.validate_role_bundle(
        role,
        [path.relative_to(ROOT).as_posix() for path in paths],
        contract,
        root=ROOT,
    )


def _role_relpaths(role: str, output_root: Path = OUTPUT_ROOT) -> list[str]:
    return [path.relative_to(ROOT).as_posix() for path in role_paths(role, output_root)]


def _file_bindings(paths: Sequence[Path]) -> list[dict[str, Any]]:
    unique = sorted({Path(path).resolve() for path in paths}, key=str)
    if FORBIDDEN_NADIG_OUTCOME.resolve() in unique:
        raise ValueError("source-only builder refuses to open or hash Nadig delta.npy")
    return [STRICT.file_binding(path, root=ROOT) for path in unique]


def write_role_bundle(
    role: str,
    arrays: Mapping[str, np.ndarray],
    *,
    transitive_inputs: Sequence[Path],
    parent_roles: Sequence[str],
    parameters: Mapping[str, Any],
    contract: Mapping[str, Any],
    output_root: Path = OUTPUT_ROOT,
) -> dict[str, Any]:
    """Write and immediately revalidate one immutable role-owned bundle."""

    npz_path, manifest_path = role_paths(role, output_root)
    for path in (npz_path, manifest_path):
        resolved = path.resolve()
        if HERE.resolve() not in resolved.parents:
            raise ValueError(f"refusing write outside Experiment 107: {path}")
        if path.exists() or path.is_symlink():
            raise FileExistsError(path)
    npz_path.parent.mkdir(parents=True, exist_ok=True)
    contiguous = {
        name: np.ascontiguousarray(np.asarray(value)) for name, value in arrays.items()
    }
    if any(value.dtype.hasobject for value in contiguous.values()):
        raise ValueError(f"{role}: object arrays are forbidden")
    try:
        with npz_path.open("xb") as handle:
            np.savez(handle, **contiguous)
            handle.flush()
            os.fsync(handle.fileno())
        os.chmod(npz_path, 0o444)
        manifest = {
            "schema_version": STRICT.ROLE_BUNDLE_SCHEMA_VERSION,
            "experiment": "107-jurkat-high-moi",
            "status": STRICT.ROLE_BUNDLE_STATUS,
            "role": role,
            "kind": role.split(".", 1)[0],
            "npz": STRICT.file_binding(npz_path, root=ROOT),
            "arrays": {
                name: STRICT.array_binding(value)
                for name, value in sorted(contiguous.items())
            },
            "parent_roles": list(parent_roles),
            "transitive_inputs": _file_bindings(transitive_inputs),
            "parameters": dict(parameters),
        }
        STRICT.write_immutable_json(manifest_path, manifest)
        return STRICT.validate_role_bundle(
            role, _role_relpaths(role, output_root), contract, root=ROOT
        )
    except BaseException:
        try:
            npz_path.chmod(0o644)
        except OSError:
            pass
        npz_path.unlink(missing_ok=True)
        try:
            manifest_path.chmod(0o644)
        except OSError:
            pass
        manifest_path.unlink(missing_ok=True)
        raise


def _donor_inputs(name: str) -> tuple[Path, ...]:
    return tuple(DONOR_PATHS[name] / filename for filename in DONOR_FILENAMES)


def _common_incumbent_inputs(donor_order: Sequence[str]) -> list[Path]:
    paths = [
        *BASE_PROVENANCE_PATHS,
        *SOURCE_PROVENANCE_PATHS,
        CHALLENGE / "manifest.json",
        CHALLENGE / "gene_names.csv",
        CHALLENGE / "pert_counts.csv",
        *PRODUCTION_IMPLEMENTATION_PATHS,
    ]
    for donor in donor_order:
        paths.extend(_donor_inputs(donor))
    paths.extend(PERTURBMATCH / filename for filename in PERTURBMATCH_FILENAMES)
    for condition in CD4_NAMES:
        paths.extend(CD4_ROOT / condition / filename for filename in CD4_FILENAMES)
    return paths


def build_donor_bundle(
    name: str,
    contract: Mapping[str, Any],
    *,
    output_root: Path = OUTPUT_ROOT,
) -> dict[str, Any]:
    existing = existing_role_bundle(
        role=f"donor.{name}", contract=contract, output_root=output_root
    )
    if existing is not None:
        return existing
    directory = DONOR_PATHS[name]
    genes = np.asarray(
        STRICT._read_axis(  # noqa: SLF001
            directory / "genes.txt", f"{name} gene axis", require_unique=False
        ),
        dtype="U",
    )
    perts = np.asarray(
        STRICT._read_axis(  # noqa: SLF001
            directory / "perts.txt", f"{name} perturbation axis", require_unique=False
        ),
        dtype="U",
    )
    delta = np.load(directory / "delta.npy", mmap_mode="r", allow_pickle=False)
    b_ctrl = np.load(directory / "b_ctrl.npy", mmap_mode="r", allow_pickle=False)
    arrays = {
        "genes": genes,
        "perts": perts,
        "delta_shape": np.asarray(delta.shape, dtype=np.int64),
        "b_ctrl_shape": np.asarray(b_ctrl.shape, dtype=np.int64),
    }
    del delta, b_ctrl
    return write_role_bundle(
        f"donor.{name}",
        arrays,
        transitive_inputs=(*BASE_PROVENANCE_PATHS, *_donor_inputs(name)),
        parent_roles=(),
        parameters={"duplicate_gene_policy": "last_occurrence_on_official_axis"},
        contract=contract,
        output_root=output_root,
    )


def _load_production_module() -> ModuleType:
    global _PRODUCTION_MODULE
    if _PRODUCTION_MODULE is not None:
        return _PRODUCTION_MODULE
    previous = {
        key: os.environ.get(key)
        for key in ("VCC_ALPHA", "VCC_CD4_FRACTION", "VCC_CD4_TARGET_F")
    }
    os.environ["VCC_ALPHA"] = "0.982"
    os.environ["VCC_CD4_FRACTION"] = "0.05"
    os.environ.pop("VCC_CD4_TARGET_F", None)
    try:
        module = _load_module("exp107_strict_a14_production", PRODUCTION_PATH)
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
    module.assert_frozen_contract()
    observed = {
        "alpha": float(module.ALPHA),
        "cd4_fraction": float(module.CD4_FRACTION),
        "cd4_target_f": module.CD4_TARGET_F,
        "gamma": float(module.WHITEN_GAMMA),
        "lambda": float(module.WHITEN_LAM),
        "pm_blend": float(module.PM_BLEND),
    }
    expected = {
        "alpha": 0.982,
        "cd4_fraction": 0.05,
        "cd4_target_f": None,
        "gamma": 0.65,
        "lambda": 0.0,
        "pm_blend": 0.5,
    }
    if observed != expected:
        raise ValueError(f"a14 production settings changed: {observed} != {expected}")
    _PRODUCTION_MODULE = module
    return module


def _build_operator_arrays(
    plan: SurfacePlan,
    *,
    max_rss_gib: float,
    chunk_rows: int,
    scratch_root: Path,
) -> dict[str, np.ndarray]:
    """Fit one production-faithful operator. This is the intentionally heavy step."""

    builder = _load_module("exp107_union_whitener_builder", WHITENER_BUILDER_PATH)
    official = np.asarray(_read_csv_column(CHALLENGE / "gene_names.csv", "gene_name"))
    target_index = {str(gene): index for index, gene in enumerate(official)}
    views = {
        name: builder.DonorView(
            name,
            DONOR_PATHS[name],
            official,
            target_index,
            chunk_rows,
            max_rss_gib,
        )
        for name in plan.donor_order
    }
    reference = float(np.median([views[name].scale for name in plan.donor_order]))
    renorm = {name: reference / (views[name].scale or 1.0) for name in plan.donor_order}
    total = sum(len(plan.basis_rows[name]) for name in plan.donor_order)
    scratch_root.mkdir(parents=True, exist_ok=True)
    basis_path = scratch_root / f"{plan.key.replace('.', '-')}-basis.npy"
    basis = np.lib.format.open_memmap(
        basis_path,
        mode="w+",
        dtype=np.float64,
        shape=(total, len(plan.operator_columns)),
    )
    cursor = 0
    try:
        for name in plan.donor_order:
            rows = np.asarray(plan.basis_rows[name], dtype=np.int64)
            views[name].block(
                rows,
                len(official),
                plan.operator_columns,
                renorm[name],
                basis,
                cursor,
                chunk_rows,
                max_rss_gib,
            )
            cursor += len(rows)
        basis.flush()
        covariance = (basis.T @ basis) / basis.shape[0]
        builder.guard(max_rss_gib, f"{plan.key} covariance")
    finally:
        builder.drop_pages(basis)
        del basis
        gc.collect()
        basis_path.unlink(missing_ok=True)
    evals, evecs = fit_whitener_eigensystem_from_covariance(covariance)
    del covariance
    builder.guard(max_rss_gib, f"{plan.key} eigendecomposition")
    arrays: dict[str, np.ndarray] = {
        "evals": np.asarray(evals, dtype=np.float64),
        "evecs": np.asarray(evecs, dtype=np.float32),
        "columns": np.asarray(plan.operator_columns, dtype=np.int64),
    }
    for name in plan.donor_order:
        arrays[f"basis_{name}"] = np.asarray(plan.basis_rows[name], dtype=np.int64)
    return arrays


def fit_whitener_eigensystem_from_covariance(
    covariance: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply the production eigensystem and floor policy to a covariance matrix."""

    matrix = np.asarray(covariance, dtype=np.float64)
    if (
        matrix.ndim != 2
        or matrix.shape[0] != matrix.shape[1]
        or not np.all(np.isfinite(matrix))
    ):
        raise ValueError("whitener covariance must be one finite square matrix")
    evals, evecs = np.linalg.eigh(matrix)
    if evals.size == 0 or float(evals.max()) <= 0:
        raise ValueError("whitener covariance has no positive eigenvalue")
    evals = np.maximum(evals, float(evals.max()) * 1e-12)
    return np.asarray(evals, dtype=np.float64), np.asarray(evecs, dtype=np.float32)


def build_operator_bundle(
    plan: SurfacePlan,
    contract: Mapping[str, Any],
    *,
    output_root: Path = OUTPUT_ROOT,
    scratch_root: Path,
    max_rss_gib: float,
    chunk_rows: int,
) -> dict[str, Any]:
    existing = existing_role_bundle(
        role=f"operator.{plan.role_suffix}",
        contract=contract,
        output_root=output_root,
    )
    if existing is not None:
        return existing
    arrays = _build_operator_arrays(
        plan,
        max_rss_gib=max_rss_gib,
        chunk_rows=chunk_rows,
        scratch_root=scratch_root,
    )
    inputs = [
        *BASE_PROVENANCE_PATHS,
        *SOURCE_PROVENANCE_PATHS,
        CHALLENGE / "gene_names.csv",
        CHALLENGE / "pert_counts.csv",
        ESSENTIAL_PATH,
        WHITENER_BUILDER_PATH,
        ROOT / "src/vcc2026/whiten.py",
    ]
    for donor in plan.donor_order:
        inputs.extend(
            DONOR_PATHS[donor] / filename
            for filename in ("delta.npy", "genes.txt", "perts.txt")
        )
        inputs.extend(role_paths(f"donor.{donor}", output_root))
    if plan.key.startswith("nadig."):
        inputs.extend(
            (
                ROOT / "data/donors/nadig_jurkat/perts.txt",
                ROOT / "data/donors/nadig_jurkat/genes.txt",
            )
        )
    elif plan.harness_path is not None:
        inputs.append(plan.harness_path)
    return write_role_bundle(
        f"operator.{plan.role_suffix}",
        arrays,
        transitive_inputs=inputs,
        parent_roles=tuple(f"donor.{name}" for name in plan.donor_order),
        parameters={
            "donor_order": list(plan.donor_order),
            "basis_row_cap": 12_000,
            "basis_rng": "numpy_default_rng_seed_0_one_stream_in_donor_order",
            "basis_prune": "common_essential_after_sampling",
            "target_row_exclusion": list(plan.targets),
            "coordinate_exclusion": "VCC targets union surface targets",
            "eigenvalue_floor_relative": 1e-12,
            "serialization": {"evals": "<f8", "evecs": "<f4"},
        },
        contract=contract,
        output_root=output_root,
    )


def _load_operator_for_inference(
    plan: SurfacePlan,
    *,
    output_root: Path,
    production: ModuleType,
):
    npz_path, _ = role_paths(f"operator.{plan.role_suffix}", output_root)
    with np.load(npz_path, allow_pickle=False) as archive:
        evals = np.asarray(archive["evals"], dtype=np.float64)
        evecs = np.asarray(archive["evecs"], dtype=np.float32)
        columns = np.asarray(archive["columns"], dtype=np.int64)
    whitener = production.Whitener(evals, evecs.astype(np.float64))
    return whitener, columns


def partition_score_axis_by_operator(
    score_indices: np.ndarray,
    operator_columns: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map score coordinates into an operator without dropping scored genes.

    A held surface may score a source-response gene that is absent from every
    incumbent donor and therefore absent from the fitted whitening operator.
    That coordinate must remain in the evaluation axis, but the source channel
    cannot transform it.  Return the score-axis support mask, the matching
    positions inside the operator, and the unsupported official coordinates.
    """

    scores = np.asarray(score_indices, dtype=np.int64)
    columns = np.asarray(operator_columns, dtype=np.int64)
    if (
        scores.ndim != 1
        or columns.ndim != 1
        or len(np.unique(scores)) != len(scores)
        or len(np.unique(columns)) != len(columns)
    ):
        raise ValueError("score indices and operator columns must be unique vectors")
    operator_position = {int(column): index for index, column in enumerate(columns)}
    supported = np.asarray(
        [int(index) in operator_position for index in scores], dtype=bool
    )
    positions = np.asarray(
        [operator_position[int(index)] for index in scores[supported]], dtype=np.int64
    )
    return supported, positions, scores[~supported]


def build_incumbent_bundle(
    plan: SurfacePlan,
    contract: Mapping[str, Any],
    *,
    output_root: Path = OUTPUT_ROOT,
    chunk_rows: int,
    max_rss_gib: float,
) -> dict[str, Any]:
    existing = existing_role_bundle(
        role=f"incumbent.{plan.role_suffix}",
        contract=contract,
        output_root=output_root,
    )
    if existing is not None:
        return existing
    production = _load_production_module()
    spec = production.ChallengeSpec.load(CHALLENGE)
    donors = {
        name: production.Donor.load(DONOR_PATHS[name]) for name in plan.donor_order
    }
    model = production.MultiDonorModel(
        donors, np.asarray(spec.genes), gene_norm=False, eb=False
    )
    weights = {name: 1.0 for name in plan.donor_order}
    pm = production.EXP061.prepare_perturbmatch(
        PERTURBMATCH,
        tuple(map(str, spec.genes)),
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    whitener, columns = _load_operator_for_inference(
        plan, output_root=output_root, production=production
    )
    official = tuple(map(str, spec.genes))
    official_index = {symbol: index for index, symbol in enumerate(official)}
    excluded_symbols = set(map(str, spec.perts)) | set(plan.targets)
    excluded_positions = np.asarray(
        sorted(
            official_index[symbol]
            for symbol in excluded_symbols
            if symbol in official_index
        ),
        dtype=np.int64,
    )
    pre = production.build_v9_pre_whiten_effects(
        model, plan.targets, alpha=0.982, weights=weights
    )
    pre, pm_audit = production.apply_perturbmatch_k562_replacement(
        pre,
        model,
        pm,
        plan.targets,
        alpha=0.982,
        weights=weights,
        blend=0.5,
    )
    pre[:, excluded_positions] = 0.0
    whitened, whitening_audit = production.apply_v9_whitening(
        pre, whitener, columns, gamma=0.65, lam=0.0
    )
    _, cd4 = production.prepare_frozen_cd4(
        CD4_ROOT,
        plan.targets,
        official,
        tuple(sorted(excluded_symbols)),
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    final, cd4_audit = production.apply_frozen_cd4_residual(
        whitened, cd4, columns, excluded_positions
    )
    final[:, excluded_positions] = 0.0
    if not np.all(np.isfinite(final)):
        raise ValueError(f"{plan.key}: incumbent contains non-finite values")
    fallback = tuple(
        target
        for target in plan.targets
        if not any(target in model.pert_index[name] for name in model.names)
    )
    if plan.key.startswith("nadig."):
        allowed = set(contract["known_coverage"]["incumbent_mean_fallback_targets"])
        if not set(fallback).issubset(allowed):
            raise ValueError(
                f"{plan.key}: unexpected incumbent fallback targets {fallback}"
            )
    arrays = {
        "targets": np.asarray(plan.targets, dtype="U"),
        "genes": np.asarray(official, dtype="U"),
        "score_genes": np.asarray(plan.score_genes, dtype="U"),
        "score_indices": np.asarray(plan.score_indices_in_official, dtype=np.int64),
        "effects": np.asarray(final, dtype=np.float64),
    }
    parent_roles = (
        *(f"donor.{name}" for name in plan.donor_order),
        f"operator.{plan.role_suffix}",
    )
    inputs = _common_incumbent_inputs(plan.donor_order)
    for donor in plan.donor_order:
        inputs.extend(role_paths(f"donor.{donor}", output_root))
    inputs.extend(role_paths(f"operator.{plan.role_suffix}", output_root))
    if plan.key.startswith("nadig."):
        inputs.extend(
            (
                ROOT / "data/donors/nadig_jurkat/perts.txt",
                ROOT / "data/donors/nadig_jurkat/genes.txt",
            )
        )
    elif plan.harness_path is not None:
        inputs.append(plan.harness_path)
    return write_role_bundle(
        f"incumbent.{plan.role_suffix}",
        arrays,
        transitive_inputs=inputs,
        parent_roles=parent_roles,
        parameters={
            "alpha": 0.982,
            "donor_weights": weights,
            "gene_norm": False,
            "empirical_bayes": False,
            "perturbmatch_fraction": 0.5,
            "cd4_fraction": 0.05,
            "cd4_target_f": None,
            "whiten_gamma": 0.65,
            "whiten_lambda": 0.0,
            "excluded_coordinate_count": int(len(excluded_positions)),
            "mean_fallback_targets": list(fallback),
            "perturbmatch_audit": pm_audit,
            "whitening_audit": whitening_audit,
            "cd4_audit": cd4_audit,
        },
        contract=contract,
        output_root=output_root,
    )


def build_integration_bundle(
    plan: SurfacePlan,
    contract: Mapping[str, Any],
    source: Any,
    *,
    output_root: Path = OUTPUT_ROOT,
) -> dict[str, Any]:
    existing = existing_role_bundle(
        role=f"integration.{plan.role_suffix}",
        contract=contract,
        output_root=output_root,
    )
    if existing is not None:
        return existing
    incumbent_path, _ = role_paths(f"incumbent.{plan.role_suffix}", output_root)
    operator_path, _ = role_paths(f"operator.{plan.role_suffix}", output_root)
    with np.load(incumbent_path, allow_pickle=False) as archive:
        incumbent = np.asarray(archive["effects"], dtype=np.float64)
        genes = np.asarray(archive["genes"])
        targets = np.asarray(archive["targets"])
        score_genes = np.asarray(archive["score_genes"])
        score_indices = np.asarray(archive["score_indices"], dtype=np.int64)
    with np.load(operator_path, allow_pickle=False) as archive:
        evals = np.asarray(archive["evals"], dtype=np.float64)
        evecs = np.asarray(archive["evecs"], dtype=np.float32)
        columns = np.asarray(archive["columns"], dtype=np.int64)
    (
        operator_supported,
        supported_positions,
        unsupported_score_indices,
    ) = partition_score_axis_by_operator(
        score_indices,
        columns,
    )
    if unsupported_score_indices.size and np.any(
        incumbent[:, unsupported_score_indices]
    ):
        unsupported_genes = genes[unsupported_score_indices].astype(str).tolist()
        raise ValueError(
            f"{plan.key}: operator-missing score coordinates have nonzero incumbent "
            f"effects: {unsupported_genes}"
        )
    raw_source = np.zeros(
        (len(plan.targets), int(operator_supported.sum())), dtype=np.float64
    )
    transferable = np.asarray(plan.source_transferable, dtype=bool).copy()
    if "MEMO1" in plan.targets:
        transferable[plan.targets.index("MEMO1")] = False
    for row, source_row in enumerate(plan.source_target_indices):
        if source_row < 0 or not transferable[row]:
            continue
        raw_source[row] = source.exported_effects[
            int(source_row), plan.score_indices_in_source[operator_supported]
        ]
    transformed = STRICT.apply_supported_principal_whitener(
        raw_source,
        evals,
        evecs,
        supported_positions,
        gamma=0.65,
        lam=0.0,
    )
    source_full = np.zeros_like(incumbent)
    source_full[:, score_indices[operator_supported]] = transformed
    coordinate_supported = np.zeros(incumbent.shape[1], dtype=bool)
    coordinate_supported[score_indices] = True
    arrays: dict[str, np.ndarray] = {
        "targets": targets,
        "genes": genes,
        "score_genes": score_genes,
        "score_indices": score_indices,
        "incumbent_effects": incumbent,
        "source_effects": source_full,
        "transferable": transferable,
        "coordinate_supported": coordinate_supported,
    }
    for fraction in STRICT.FRACTIONS:
        suffix = f"f{int(round(fraction * 1000)):03d}"
        candidate, effective = STRICT.blend_with_source(
            incumbent,
            source_full,
            fraction,
            transferable,
            coordinate_supported,
        )
        arrays[f"candidate_{suffix}"] = candidate
        arrays[f"effective_{suffix}"] = effective
    parent_roles = (
        f"incumbent.{plan.role_suffix}",
        f"operator.{plan.role_suffix}",
    )
    inputs = [
        *BASE_PROVENANCE_PATHS,
        *SOURCE_PROVENANCE_PATHS,
        *role_paths(f"incumbent.{plan.role_suffix}", output_root),
        *role_paths(f"operator.{plan.role_suffix}", output_root),
    ]
    if plan.key.startswith("nadig."):
        inputs.extend(
            (
                ROOT / "data/donors/nadig_jurkat/perts.txt",
                ROOT / "data/donors/nadig_jurkat/genes.txt",
            )
        )
    elif plan.harness_path is not None:
        inputs.append(plan.harness_path)
    return write_role_bundle(
        f"integration.{plan.role_suffix}",
        arrays,
        transitive_inputs=inputs,
        parent_roles=parent_roles,
        parameters={
            "blend_grid": list(STRICT.FRACTIONS),
            "source_operator": "supported_principal_submatrix",
            "source_row_scaling": "match_incumbent_supported_row_norm",
            "unsupported_policy": "bit_exact_incumbent_fallback",
            "operator_missing_score_genes": genes[
                unsupported_score_indices
            ].astype(str).tolist(),
            "operator_missing_score_policy": (
                "retain_in_score_axis_require_zero_incumbent_and_source"
            ),
            "ambiguous_MEMO1_forced_unsupported": True,
        },
        contract=contract,
        output_root=output_root,
    )


def write_top_level_manifest(
    contract: Mapping[str, Any],
    *,
    output_root: Path = OUTPUT_ROOT,
    manifest_path: Path = TOP_LEVEL_MANIFEST,
) -> dict[str, Any]:
    if manifest_path.is_file() and not manifest_path.is_symlink():
        return STRICT.bind_artifact_manifest(
            manifest_path,
            required_roles=STRICT.REQUIRED_ARTIFACT_ROLES,
            root=ROOT,
            contract=contract,
        )
    if manifest_path.exists() or manifest_path.is_symlink():
        raise RuntimeError(f"unsafe top-level manifest path: {manifest_path}")
    artifacts = {
        role: _role_relpaths(role, output_root)
        for role in STRICT.REQUIRED_ARTIFACT_ROLES
    }
    payload = {
        "schema_version": 1,
        "experiment": "107-jurkat-high-moi",
        "status": STRICT.ROLE_BUNDLE_STATUS,
        "artifacts": artifacts,
    }
    try:
        STRICT.write_immutable_json(manifest_path, payload)
        return STRICT.bind_artifact_manifest(
            manifest_path,
            required_roles=STRICT.REQUIRED_ARTIFACT_ROLES,
            root=ROOT,
            contract=contract,
        )
    except BaseException:
        try:
            manifest_path.chmod(0o644)
        except OSError:
            pass
        manifest_path.unlink(missing_ok=True)
        raise


def metadata_report() -> dict[str, Any]:
    plans, contract, source = build_surface_plans()
    return {
        "status": "metadata_only_complete",
        "nadig_outcome_opened_or_hashed": False,
        "contract_sha256": STRICT.sha256_file(STRICT.CONTRACT_PATH),
        "source_model_sha256": source.model_sha256,
        "roles": list(STRICT.REQUIRED_ARTIFACT_ROLES),
        "surfaces": {
            key: {
                "targets": len(plan.targets),
                "score_genes": len(plan.score_genes),
                "operator_columns": len(plan.operator_columns),
                "basis_rows": {
                    name: len(plan.basis_rows[name]) for name in plan.donor_order
                },
                "donor_order": list(plan.donor_order),
                "transferable_targets": int(plan.source_transferable.sum()),
            }
            for key, plan in plans.items()
        },
        "held_surfaces": list(
            contract["held_context_gate"]["authorizing_directional_surfaces"]
        ),
    }


def build_all(
    *,
    output_root: Path,
    manifest_path: Path,
    scratch_root: Path,
    max_rss_gib: float,
    chunk_rows: int,
) -> dict[str, Any]:
    plans, contract, source = build_surface_plans()
    for donor in STRICT.DONOR_ORDER:
        build_donor_bundle(donor, contract, output_root=output_root)
    for key in (
        *(f"nadig.gallery_{number}" for number in range(1, 7)),
        *STRICT.AUTHORIZING_HELD_SURFACES,
    ):
        plan = plans[key]
        build_operator_bundle(
            plan,
            contract,
            output_root=output_root,
            scratch_root=scratch_root,
            max_rss_gib=max_rss_gib,
            chunk_rows=chunk_rows,
        )
        build_incumbent_bundle(
            plan,
            contract,
            output_root=output_root,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        build_integration_bundle(
            plan,
            contract,
            source,
            output_root=output_root,
        )
        gc.collect()
    return write_top_level_manifest(
        contract, output_root=output_root, manifest_path=manifest_path
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    commands.add_parser("check", help="metadata-only contract and geometry validation")
    build = commands.add_parser(
        "build", help="build all frozen source-only role bundles"
    )
    build.add_argument("--confirm-heavy-eigendecomposition", action="store_true")
    build.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    build.add_argument("--manifest", type=Path, default=TOP_LEVEL_MANIFEST)
    build.add_argument("--scratch-root", type=Path, default=OUTPUT_ROOT / ".scratch")
    build.add_argument("--max-rss-gib", type=float, default=24.0)
    build.add_argument("--chunk-rows", type=int, default=1024)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.command == "check":
        print(json.dumps(metadata_report(), indent=2, sort_keys=True))
        return 0
    if not args.confirm_heavy_eigendecomposition:
        raise SystemExit(
            "build requires --confirm-heavy-eigendecomposition; use `check` for metadata-only validation"
        )
    output_root = args.output_root.resolve()
    manifest = args.manifest.resolve()
    scratch = args.scratch_root.resolve()
    if output_root != OUTPUT_ROOT.resolve():
        raise SystemExit("output root is frozen by the strict contract")
    if manifest != TOP_LEVEL_MANIFEST.resolve():
        raise SystemExit("top-level manifest path is frozen by the strict contract")
    if HERE.resolve() not in output_root.parents:
        raise SystemExit("output root must stay inside Experiment 107")
    if HERE.resolve() not in manifest.parents:
        raise SystemExit("manifest must stay inside Experiment 107")
    if HERE.resolve() not in scratch.parents:
        raise SystemExit("scratch root must stay inside Experiment 107")
    try:
        result = build_all(
            output_root=output_root,
            manifest_path=manifest,
            scratch_root=scratch,
            max_rss_gib=args.max_rss_gib,
            chunk_rows=args.chunk_rows,
        )
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    print(json.dumps({"status": "built", "binding_sha256": result["binding_sha256"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

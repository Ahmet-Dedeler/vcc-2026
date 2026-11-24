#!/usr/bin/env python3
"""Build immutable, role-owned source inputs for Experiment 114.

The real build is deliberately separate from candidate freezing.  It consumes
only public/source artifacts plus Experiment 107's already-frozen role-specific
whitening operator.  It never opens a held expression matrix.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUN_PATH = HERE / "run.py"
ROLE_SOURCE_MANIFEST_PATH = HERE / "ROLE_SOURCE_MANIFEST.json"
ROLE_OUTPUT_ROOT = HERE / "work" / "roles"
EXPECTED_ROLE_SOURCE_MANIFEST_SHA256 = (
    "e05d877b5ea0a7aeae06f4c53ad363e984379678b5b61fd852e3e877d17d0248"
)
EXPECTED_ROLES = ("h1", "orion_hek293t")
SOURCE_BANK_ARMS = (
    "string_neighbor",
    "depmap_neighbor",
    "string_depmap_blend_w1",
)
FIXED_ZIP_TIMESTAMP = (1980, 1, 1, 0, 0, 0)
EFFECT_SCALE = 6.817


def _load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


RUN = _load_module("experiment_114_role_input_run", RUN_PATH)
ContractError = RUN.ContractError


def canonical_json_sha256(value: Mapping[str, Any]) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_int64(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values, dtype=np.int64))
    return hashlib.sha256(array.tobytes(order="C")).hexdigest()


def exp107_sha256_array(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values))
    digest = hashlib.sha256()
    digest.update(b"exp107-array-v1\0")
    digest.update(array.dtype.str.encode("ascii"))
    digest.update(b"\0")
    digest.update(json.dumps(list(array.shape), separators=(",", ":")).encode("ascii"))
    digest.update(b"\0")
    digest.update(memoryview(array).cast("B"))
    return digest.hexdigest()


def sha256_array_stream(values: np.ndarray, chunk_bytes: int = 64 << 20) -> str:
    """Match run.py's logical-array hash without materializing a giant bytes object."""

    array = np.asarray(values)
    if not array.flags.c_contiguous:
        raise ContractError("logical arrays must be C-contiguous before sealing")
    header = json.dumps(
        {"dtype": array.dtype.str, "shape": list(array.shape)},
        separators=(",", ":"),
        sort_keys=True,
    ).encode()
    digest = hashlib.sha256()
    digest.update(header)
    digest.update(b"\0")
    payload = memoryview(array).cast("B")
    for offset in range(0, len(payload), chunk_bytes):
        digest.update(payload[offset : offset + chunk_bytes])
    return digest.hexdigest()


def _record_path(record: Mapping[str, Any], label: str) -> Path:
    if set(record) != {"path", "sha256"}:
        raise ContractError(f"{label} must contain exactly path and sha256")
    return RUN.resolve_repository_path(str(record["path"]), label)


def verify_record(
    record: Mapping[str, Any],
    label: str,
    *,
    hash_cache: dict[Path, str] | None = None,
) -> Path:
    path = _record_path(record, label)
    if not path.is_file():
        raise ContractError(f"{label} is missing: {path}")
    resolved = path.resolve()
    observed = None if hash_cache is None else hash_cache.get(resolved)
    if observed is None:
        observed = RUN.sha256_file(path)
        if hash_cache is not None:
            hash_cache[resolved] = observed
    if observed != record["sha256"]:
        raise ContractError(f"{label} hash drifted")
    return path


def verify_axis_record(
    record: Mapping[str, Any],
    label: str,
    *,
    hash_cache: dict[Path, str] | None = None,
) -> Path:
    expected_keys = {"path", "sha256", "column", "count", "ordered_sha256"}
    if set(record) != expected_keys:
        raise ContractError(f"{label} schema drifted")
    return verify_record(
        {"path": record["path"], "sha256": record["sha256"]},
        label,
        hash_cache=hash_cache,
    )


def receipt_record(
    path: Path,
    *,
    hash_cache: dict[Path, str] | None = None,
) -> dict[str, Any]:
    resolved = path.resolve()
    if not resolved.is_relative_to(ROOT.resolve()):
        raise ContractError(f"receipt input escapes repository: {resolved}")
    RUN.reject_symlink_components(resolved, "receipt input")
    if not resolved.is_file():
        raise ContractError(f"receipt input is missing: {resolved}")
    digest = None if hash_cache is None else hash_cache.get(resolved)
    if digest is None:
        digest = RUN.sha256_file(resolved)
        if hash_cache is not None:
            hash_cache[resolved] = digest
    return {
        "path": resolved.relative_to(ROOT.resolve()).as_posix(),
        "sha256": digest,
        "source_only": True,
        "contains_held_truth": False,
    }


def binding_record(path: Path) -> dict[str, str]:
    record = receipt_record(path)
    return {"path": record["path"], "sha256": record["sha256"]}


def _iter_records(value: Any) -> Iterable[Mapping[str, Any]]:
    if isinstance(value, dict):
        if isinstance(value.get("path"), str) and isinstance(value.get("sha256"), str):
            yield {"path": value["path"], "sha256": value["sha256"]}
        else:
            for child in value.values():
                yield from _iter_records(child)
    elif isinstance(value, list):
        for child in value:
            yield from _iter_records(child)


def load_role_source_manifest() -> dict[str, Any]:
    if RUN.sha256_file(ROLE_SOURCE_MANIFEST_PATH) != (
        EXPECTED_ROLE_SOURCE_MANIFEST_SHA256
    ):
        raise ContractError("role source manifest changed")
    value = RUN.load_json(ROLE_SOURCE_MANIFEST_PATH)
    validate_role_source_manifest(value)
    return value


def validate_role_source_manifest(value: Mapping[str, Any]) -> None:
    expected_top = {
        "schema_version",
        "experiment",
        "date_frozen",
        "status",
        "canonical_axes",
        "eligibility",
        "a14",
        "source_bank",
        "upstream_contracts",
        "implementation_inputs",
        "source_inputs",
        "roles",
        "authorization",
    }
    if set(value) != expected_top:
        raise ContractError("role source manifest top-level schema drifted")
    if value.get("schema_version") != 1:
        raise ContractError("unsupported role source manifest schema")
    if value.get("experiment") != "114-depmap-coessentiality":
        raise ContractError("role source manifest belongs to another experiment")
    if (
        value.get("status")
        != "builders_defined_real_role_bundles_not_built_or_approved"
    ):
        raise ContractError("role source manifest status drifted")
    expected_axes = {
        "targets": {
            "path": "data/challenge/pert_counts.csv",
            "sha256": (
                "f57edd7b912ebd718efc7ee9d0f334772513e7cc418d133ce525470e373b3276"
            ),
            "column": "target_gene",
            "count": 300,
            "ordered_sha256": (
                "cd3d00f1cd5da660609ad41d85b53e2d5557f7f5edb3e1f64559f78f7fb3cbc4"
            ),
        },
        "genes": {
            "path": "data/challenge/gene_names.csv",
            "sha256": (
                "25bfa66715e186bebabce7ac788bbcea47e2bf59ca70be1f8f3a06f2f0e47201"
            ),
            "column": "gene_name",
            "count": 18_533,
            "ordered_sha256": (
                "c2b8dec58c02810d64018164155b655c4e7c7be451a3ae579965db193af33723"
            ),
        },
    }
    if value.get("canonical_axes") != expected_axes:
        raise ContractError("role source canonical axes drifted")
    a14 = value.get("a14")
    expected_a14 = {
        "alpha": 0.982,
        "effect_scale": EFFECT_SCALE,
        "donor_weights": "equal_over_role_donors_that_carry_the_target",
        "gene_norm": False,
        "empirical_bayes": False,
        "beta": 0.0,
        "perturbmatch_blend": 0.5,
        "perturbmatch_policy": (
            "replace_only_the_k562_slot_when_both_k562_and_perturbmatch_carry_the_target"
        ),
        "cd4_fraction": 0.05,
        "cd4_target_solve": None,
        "cd4_policy": "complete_three_condition_mean_residual_after_whitening",
        "whitening_gamma": 0.65,
        "whitening_lambda": 0.0,
        "basis_prune": "common_essential_after_sampling",
        "basis_row_cap_per_donor": 12_000,
        "basis_rng": "numpy_default_rng_seed_0_one_stream_in_role_donor_order",
        "basis_coordinate_exclusion": (
            "union_of_canonical_300_vcc_targets_and_role_held_surface_targets"
        ),
        "basis_target_row_exclusion": "all_role_held_surface_targets",
        "eigenvalue_floor_relative": 1e-12,
        "eigenvalue_dtype": "<f8",
        "eigenvector_dtype": "<f4",
        "duplicate_gene_policy": "last_occurrence_on_canonical_axis",
        "incumbent_effect_dtype": "<f8",
        "source_effect_dtype": "<f8",
    }
    if a14 != expected_a14:
        raise ContractError("exact a14 settings drifted")
    source_bank = value.get("source_bank")
    if not isinstance(source_bank, dict):
        raise ContractError("source-bank construction contract is missing")
    expected_source_fields = {
        "neighbours": 25,
        "similarity_floor": 0.0,
        "feature_blend_weight": 1.0,
        "perturbmatch_in_source_bank": False,
        "cd4_in_source_bank": False,
        "effect_space": ("a14_role_donor_deviation_role_whitened_no_global_renorm_v1"),
    }
    for key, expected in expected_source_fields.items():
        if source_bank.get(key) != expected:
            raise ContractError(f"source-bank contract drifted: {key}")
    authorization = value.get("authorization")
    if authorization != {
        "real_role_build": True,
        "source_candidate_freeze": False,
        "held_truth_access": False,
        "evaluation": False,
        "count_generation": False,
        "submission": False,
        "approved_roles": {},
    }:
        raise ContractError("role-source authorization map drifted")
    roles = value.get("roles")
    if not isinstance(roles, dict) or tuple(roles) != EXPECTED_ROLES:
        raise ContractError("role order or names drifted")
    expected_roles = {
        "h1": {
            "surface": "h1",
            "held": "h1_2025",
            "allowed": ["k562_gwps", "orion_hct116", "orion_hek293t"],
            "forbidden": ["h1_2025"],
            "operator_role": "operator.h1",
            "columns": 17_862,
            "column_hash": (
                "00242e6a8fdd1f5c6dd5b8ee6a511a949c730537b9352865d653a1361d80a622"
            ),
            "held_targets": 85,
            "held_hash": (
                "f49f8343fe877c378fc258884def4a66241adac43d4c2798e9257cecd6d6bd1f"
            ),
        },
        "orion_hek293t": {
            "surface": "orion_hek293t",
            "held": "orion_hek293t",
            "allowed": ["k562_gwps", "orion_hct116"],
            "forbidden": ["orion_hek293t"],
            "operator_role": "operator.orion_hek293t",
            "columns": 17_942,
            "column_hash": (
                "b2cae0d381c3cd7065a3de33f1a6da5514577f2b3187b6dc84786ee789e27e4d"
            ),
            "held_targets": 132,
            "held_hash": (
                "d4fbc4d5e296f9eb73a87a4df84e01ae91aead80aa52e110bd166782f87fccca"
            ),
        },
    }
    for role, expected in expected_roles.items():
        record = roles[role]
        if record.get("surface") != expected["surface"]:
            raise ContractError(f"{role}: surface drifted")
        if record.get("held_surface_source_id") != expected["held"]:
            raise ContractError(f"{role}: held source ID drifted")
        if record.get("allowed_source_ids") != expected["allowed"]:
            raise ContractError(f"{role}: source allowlist drifted")
        if record.get("forbidden_source_ids") != expected["forbidden"]:
            raise ContractError(f"{role}: source denylist drifted")
        if set(expected["allowed"]) & set(expected["forbidden"]):
            raise ContractError(f"{role}: source allowlist intersects denylist")
        operator = record.get("operator", {})
        if operator.get("role") != expected["operator_role"]:
            raise ContractError(f"{role}: operator role drifted")
        if operator.get("column_count") != expected["columns"]:
            raise ContractError(f"{role}: operator column count drifted")
        if operator.get("columns_int64_sha256") != expected["column_hash"]:
            raise ContractError(f"{role}: operator column hash drifted")
        if operator.get("target_exclusion_count") != expected["held_targets"]:
            raise ContractError(f"{role}: held target count drifted")
        if operator.get("target_exclusion_ordered_sha256") != expected["held_hash"]:
            raise ContractError(f"{role}: held target hash drifted")
        expected_output = f"experiments/114-depmap-coessentiality/work/roles/{role}"
        if record.get("output_directory") != expected_output:
            raise ContractError(f"{role}: output directory drifted")
        if record.get("incumbent_endpoint") != ("a14-role-heldout-canonical-300-v1"):
            raise ContractError(f"{role}: incumbent endpoint drifted")
        if record.get("source_effect_space") != source_bank["effect_space"]:
            raise ContractError(f"{role}: source effect space drifted")
    for index, record in enumerate(_iter_records(value["source_inputs"])):
        _record_path(record, f"source input record {index}")
    for index, record in enumerate(_iter_records(value["upstream_contracts"])):
        _record_path(record, f"upstream contract record {index}")
    implementation_inputs = value.get("implementation_inputs")
    if not isinstance(implementation_inputs, dict) or not implementation_inputs:
        raise ContractError("implementation input bindings are missing")
    for path, digest in implementation_inputs.items():
        _record_path({"path": path, "sha256": digest}, f"implementation {path}")


def validate_checked_in_locks() -> None:
    protocol = RUN.load_json(RUN.PROTOCOL_PATH)
    source_manifest = RUN.load_json(RUN.SOURCE_MANIFEST_PATH)
    boundary = RUN.load_json(RUN.SOURCE_CANDIDATE_FREEZE_BOUNDARY_PATH)
    RUN.validate_protocol(protocol)
    RUN.validate_source_manifest(source_manifest)
    RUN.validate_source_candidate_freeze_boundary(boundary)
    if protocol["current_authorization"]["source_candidate_freeze"] is not False:
        raise ContractError("source_candidate_freeze unexpectedly enabled in protocol")
    if boundary["source_candidate_freeze_authorized"] is not False:
        raise ContractError("source_candidate_freeze unexpectedly enabled in boundary")
    if boundary["approved_roles"] != {}:
        raise ContractError(
            "boundary contains role approvals before independent review"
        )


def active_exp107_build_pids() -> tuple[int, ...]:
    completed = subprocess.run(
        ["ps", "-axo", "pid=,command="],
        check=True,
        text=True,
        capture_output=True,
    )
    found: list[int] = []
    for line in completed.stdout.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        pid_text, _, command = stripped.partition(" ")
        if (
            "experiments/107-jurkat-high-moi/strict_artifact_builder.py build"
            in command
        ):
            found.append(int(pid_text))
    return tuple(found)


def read_axis(path: Path, *, unique: bool, label: str) -> tuple[str, ...]:
    values = tuple(line.strip() for line in path.read_text().splitlines())
    if not values or any(not value for value in values):
        raise ContractError(f"{label} contains blank values")
    if unique and len(values) != len(set(values)):
        raise ContractError(f"{label} contains duplicate values")
    return values


def last_occurrence_map(
    source_genes: Sequence[str],
    canonical_genes: Sequence[str],
) -> tuple[np.ndarray, np.ndarray]:
    canonical_index = {gene: index for index, gene in enumerate(canonical_genes)}
    last: dict[int, int] = {}
    for source_index, gene in enumerate(map(str, source_genes)):
        target_index = canonical_index.get(gene)
        if target_index is not None:
            last[target_index] = source_index
    target_columns = np.asarray(sorted(last), dtype=np.int64)
    source_columns = np.asarray(
        [last[int(target)] for target in target_columns], dtype=np.int64
    )
    return target_columns, source_columns


@dataclass
class RoleDonor:
    name: str
    delta_path: Path
    genes: np.ndarray
    perts: np.ndarray
    mean: np.ndarray
    tcols_all: np.ndarray
    dcols_all: np.ndarray
    tcols_unique: np.ndarray
    dcols_unique: np.ndarray
    pert_index: dict[str, int]
    scale: float
    delta: np.ndarray

    @classmethod
    def load(
        cls,
        name: str,
        records: Mapping[str, Mapping[str, Any]],
        canonical_genes: np.ndarray,
        *,
        chunk_rows: int,
        hash_cache: dict[Path, str] | None = None,
    ) -> "RoleDonor":
        delta_path = verify_record(
            records["delta"], f"{name} delta", hash_cache=hash_cache
        )
        genes_path = verify_record(
            records["genes"], f"{name} genes", hash_cache=hash_cache
        )
        perts_path = verify_record(
            records["perts"], f"{name} perts", hash_cache=hash_cache
        )
        verify_record(records["b_ctrl"], f"{name} control", hash_cache=hash_cache)
        verify_record(records["meta"], f"{name} metadata", hash_cache=hash_cache)
        genes = np.asarray(
            read_axis(genes_path, unique=False, label=f"{name} gene axis"),
            dtype=str,
        )
        perts = np.asarray(
            read_axis(perts_path, unique=True, label=f"{name} target axis"),
            dtype=str,
        )
        delta = np.load(delta_path, mmap_mode="r", allow_pickle=False)
        if delta.shape != (len(perts), len(genes)) or delta.dtype != np.float32:
            raise ContractError(f"{name}: donor delta shape or dtype drifted")
        mean = np.asarray(delta.mean(axis=0))
        target_index = {gene: index for index, gene in enumerate(canonical_genes)}
        tcols_all: list[int] = []
        dcols_all: list[int] = []
        for source_column, gene in enumerate(genes):
            target_column = target_index.get(str(gene))
            if target_column is not None:
                tcols_all.append(target_column)
                dcols_all.append(source_column)
        tcols = np.asarray(tcols_all, dtype=np.int64)
        dcols = np.asarray(dcols_all, dtype=np.int64)
        unique_tcols, unique_dcols = last_occurrence_map(genes, canonical_genes)
        norms = np.empty(len(perts), dtype=np.float32)
        centered_mean = mean[dcols]
        for start in range(0, len(perts), chunk_rows):
            stop = min(start + chunk_rows, len(perts))
            block = np.asarray(delta[start:stop][:, dcols], dtype=np.float32)
            block -= centered_mean
            norms[start:stop] = np.linalg.norm(block, axis=1)
        scale = float(np.median(norms))
        if not np.isfinite(scale) or scale <= 0:
            raise ContractError(f"{name}: invalid donor effect scale")
        return cls(
            name=name,
            delta_path=delta_path,
            genes=genes,
            perts=perts,
            mean=mean,
            tcols_all=tcols,
            dcols_all=dcols,
            tcols_unique=unique_tcols,
            dcols_unique=unique_dcols,
            pert_index={target: row for row, target in enumerate(perts.tolist())},
            scale=scale,
            delta=delta,
        )

    def centered_row(self, row: int) -> np.ndarray:
        return np.asarray(self.delta[row], dtype=np.float32) - self.mean

    def canonical_unscaled_row(self, row: int, n_genes: int) -> np.ndarray:
        output = np.zeros(n_genes, dtype=np.float64)
        centered = self.centered_row(row)
        output[self.tcols_all] = centered[self.dcols_all]
        return output

    def close(self) -> None:
        mmap = getattr(self.delta, "_mmap", None)
        if mmap is not None:
            mmap.close()


class _CenteredRows:
    def __init__(self, donor: RoleDonor):
        self.donor = donor

    def __getitem__(self, row: int) -> np.ndarray:
        return self.donor.centered_row(int(row))


class StreamingA14Model:
    """The exact gene_norm=False, EB=False MultiDonorModel arithmetic for 300 rows."""

    def __init__(self, donors: Mapping[str, RoleDonor], canonical_genes: np.ndarray):
        self.names = list(donors)
        self.target_genes = np.asarray(canonical_genes)
        self.effect_scale = EFFECT_SCALE
        self.tcols = {name: donor.tcols_all for name, donor in donors.items()}
        self.dcols = {name: donor.dcols_all for name, donor in donors.items()}
        self.s = {name: _CenteredRows(donor) for name, donor in donors.items()}
        self.pert_index = {name: donor.pert_index for name, donor in donors.items()}
        self._scale = {name: donor.scale for name, donor in donors.items()}
        reference = float(np.median(list(self._scale.values())))
        self._renorm = {
            name: reference / (self._scale[name] or 1.0) for name in self.names
        }
        self._global = EFFECT_SCALE / reference
        self._median_s_norm = EFFECT_SCALE
        self.m_full = np.zeros(len(canonical_genes), dtype=np.float32)
        for donor in donors.values():
            mapped = np.zeros(len(canonical_genes), dtype=np.float32)
            mapped[donor.tcols_all] = donor.mean[donor.dcols_all]
            self.m_full += mapped / len(donors)
        self._m_norm = float(np.linalg.norm(self.m_full)) or 1.0
        self._donors = dict(donors)

    def effect(
        self,
        pert: str,
        alpha: float,
        beta: float,
        *,
        weights: Mapping[str, float] | None = None,
    ) -> np.ndarray:
        effect = beta * self.m_full
        have = [name for name in self.names if pert in self.pert_index[name]]
        if not have:
            if alpha != 0.0:
                effect = (
                    effect + (alpha * self._median_s_norm / self._m_norm) * self.m_full
                )
            return effect
        active = (
            {name: float(weights.get(name, 0.0)) for name in have}
            if weights is not None
            else {name: 1.0 for name in have}
        )
        total = sum(active.values()) or 1.0
        accumulator = np.zeros(len(self.target_genes), dtype=np.float64)
        for name in have:
            donor = self._donors[name]
            vector = donor.canonical_unscaled_row(
                self.pert_index[name][pert], len(self.target_genes)
            )
            accumulator += (active[name] / total) * self._renorm[name] * vector
        return effect + (alpha * self._global) * accumulator.astype(np.float32)


def role_scaled_source_row(
    donor: RoleDonor,
    row: int,
    canonical_gene_count: int,
    role_reference_scale: float,
) -> np.ndarray:
    """Return one donor slot with MultiDonorModel's exact float32 rounding point."""

    raw = donor.canonical_unscaled_row(row, canonical_gene_count)
    renormalized = ((role_reference_scale / donor.scale) * raw).astype(np.float32)
    return ((0.982 * EFFECT_SCALE / role_reference_scale) * renormalized).astype(
        np.float64
    )


def _load_production_module() -> ModuleType:
    previous = {
        key: os.environ.get(key)
        for key in ("VCC_ALPHA", "VCC_CD4_FRACTION", "VCC_CD4_TARGET_F")
    }
    os.environ["VCC_ALPHA"] = "0.982"
    os.environ["VCC_CD4_FRACTION"] = "0.05"
    os.environ.pop("VCC_CD4_TARGET_F", None)
    try:
        module = _load_module(
            "experiment_114_role_a14_production",
            ROOT / "scripts" / "make_submission_pm_cd4.py",
        )
    finally:
        for key, value in previous.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
    observed = {
        "alpha": float(module.ALPHA),
        "beta": float(module.BETA),
        "pm": float(module.PM_BLEND),
        "cd4": float(module.CD4_FRACTION),
        "cd4_target": module.CD4_TARGET_F,
        "gamma": float(module.WHITEN_GAMMA),
        "lambda": float(module.WHITEN_LAM),
    }
    expected = {
        "alpha": 0.982,
        "beta": 0.0,
        "pm": 0.5,
        "cd4": 0.05,
        "cd4_target": None,
        "gamma": 0.65,
        "lambda": 0.0,
    }
    if observed != expected:
        raise ContractError(f"production a14 settings drifted: {observed}")
    return module


def _exp107_array_binding(value: np.ndarray) -> dict[str, Any]:
    array = np.ascontiguousarray(np.asarray(value))
    return {
        "dtype": array.dtype.str,
        "shape": list(array.shape),
        "c_order": True,
        "sha256_array": exp107_sha256_array(array),
    }


@dataclass(frozen=True)
class RoleOperator:
    evals: np.ndarray
    evecs: np.ndarray
    columns: np.ndarray
    held_targets: tuple[str, ...]
    artifact_path: Path
    manifest_path: Path
    artifact_sha256: str
    manifest_sha256: str


def load_role_operator(role: str, role_spec: Mapping[str, Any]) -> RoleOperator:
    spec = role_spec["operator"]
    artifact = RUN.resolve_repository_path(spec["artifact_path"], f"{role} operator")
    manifest_path = RUN.resolve_repository_path(
        spec["manifest_path"], f"{role} operator manifest"
    )
    for path, label in (
        (artifact, f"{role} operator"),
        (manifest_path, f"{role} operator manifest"),
    ):
        RUN.reject_symlink_components(path, label)
        if not path.is_file():
            raise ContractError(
                f"{label} is not ready; finish Experiment 107 before this role build"
            )
    manifest = RUN.load_json(manifest_path)
    expected_keys = {
        "schema_version",
        "experiment",
        "status",
        "role",
        "kind",
        "npz",
        "arrays",
        "parent_roles",
        "transitive_inputs",
        "parameters",
    }
    if set(manifest) != expected_keys:
        raise ContractError(f"{role}: upstream operator manifest schema drifted")
    if (
        manifest.get("schema_version") != 1
        or manifest.get("experiment") != "107-jurkat-high-moi"
        or manifest.get("status") != "frozen_before_experiment_107_nadig_outcome_access"
        or manifest.get("role") != spec["role"]
        or manifest.get("kind") != "operator"
    ):
        raise ContractError(f"{role}: upstream operator identity drifted")
    expected_parents = [
        f"donor.{source_id}" for source_id in role_spec["allowed_source_ids"]
    ]
    if manifest.get("parent_roles") != expected_parents:
        raise ContractError(f"{role}: upstream operator parent roles drifted")
    npz_binding = manifest.get("npz")
    if not isinstance(npz_binding, dict) or set(npz_binding) != {
        "path",
        "bytes",
        "sha256",
    }:
        raise ContractError(f"{role}: upstream operator NPZ binding is invalid")
    if npz_binding["path"] != spec["artifact_path"]:
        raise ContractError(f"{role}: upstream operator path drifted")
    if artifact.stat().st_size != int(npz_binding["bytes"]):
        raise ContractError(f"{role}: upstream operator byte count drifted")
    artifact_sha256 = RUN.sha256_file(artifact)
    if artifact_sha256 != npz_binding["sha256"]:
        raise ContractError(f"{role}: upstream operator hash drifted")
    parameters = manifest.get("parameters")
    expected_parameters = {
        "donor_order": role_spec["allowed_source_ids"],
        "basis_row_cap": 12_000,
        "basis_rng": "numpy_default_rng_seed_0_one_stream_in_donor_order",
        "basis_prune": "common_essential_after_sampling",
        "coordinate_exclusion": "VCC targets union surface targets",
        "eigenvalue_floor_relative": 1e-12,
        "serialization": {"evals": "<f8", "evecs": "<f4"},
    }
    if not isinstance(parameters, dict):
        raise ContractError(f"{role}: upstream operator parameters are missing")
    for key, expected in expected_parameters.items():
        if parameters.get(key) != expected:
            raise ContractError(f"{role}: upstream operator parameter drifted: {key}")
    held_targets = tuple(map(str, parameters.get("target_row_exclusion", [])))
    if len(held_targets) != int(spec["target_exclusion_count"]):
        raise ContractError(f"{role}: held target exclusion count drifted")
    if RUN.sha256_lines(held_targets) != spec["target_exclusion_ordered_sha256"]:
        raise ContractError(f"{role}: held target exclusion axis drifted")
    expected_arrays = {
        "evals",
        "evecs",
        "columns",
        *(f"basis_{name}" for name in role_spec["allowed_source_ids"]),
    }
    bindings = manifest.get("arrays")
    if not isinstance(bindings, dict) or set(bindings) != expected_arrays:
        raise ContractError(f"{role}: upstream operator array schema drifted")
    with np.load(artifact, allow_pickle=False) as archive:
        if set(archive.files) != expected_arrays:
            raise ContractError(f"{role}: upstream operator NPZ keys drifted")
        arrays = {name: np.asarray(archive[name]) for name in archive.files}
    for name, array in arrays.items():
        if bindings[name] != _exp107_array_binding(array):
            raise ContractError(f"{role}: upstream operator array drifted: {name}")
    evals = np.asarray(arrays["evals"], dtype=np.float64)
    evecs = np.asarray(arrays["evecs"], dtype=np.float64)
    columns = np.asarray(arrays["columns"], dtype=np.int64)
    count = int(spec["column_count"])
    if (
        evals.shape != (count,)
        or evecs.shape != (count, count)
        or columns.shape != (count,)
    ):
        raise ContractError(f"{role}: upstream operator dimensions drifted")
    if (
        not np.isfinite(evals).all()
        or np.any(evals <= 0)
        or not np.isfinite(evecs).all()
        or not np.array_equal(columns, np.unique(columns))
        or np.any(columns < 0)
        or np.any(columns >= 18_533)
    ):
        raise ContractError(f"{role}: upstream operator numerical contract failed")
    if sha256_int64(columns) != spec["columns_int64_sha256"]:
        raise ContractError(f"{role}: upstream operator column axis drifted")
    for donor, expected in spec["basis_rows_by_donor"].items():
        basis = np.asarray(arrays[f"basis_{donor}"], dtype=np.int64)
        if len(basis) != int(expected["count"]):
            raise ContractError(f"{role}: {donor} basis count drifted")
        if sha256_int64(basis) != expected["int64_sha256"]:
            raise ContractError(f"{role}: {donor} basis rows drifted")
    return RoleOperator(
        evals=evals,
        evecs=evecs,
        columns=columns,
        held_targets=held_targets,
        artifact_path=artifact,
        manifest_path=manifest_path,
        artifact_sha256=artifact_sha256,
        manifest_sha256=RUN.sha256_file(manifest_path),
    )


def source_target_closure(
    source_targets: Sequence[str],
    eligible_targets: Sequence[str],
    features_by_arm: Mapping[str, Any],
) -> tuple[tuple[str, ...], dict[str, Any]]:
    ordered_sources = np.asarray(sorted(set(map(str, source_targets))), dtype=str)
    selected_union: set[str] = set()
    arm_audit: dict[str, Any] = {}
    for arm in SOURCE_BANK_ARMS:
        features = features_by_arm[arm]
        source_matrix, source_available = RUN.feature_rows(features, ordered_sources)
        feature_index = {
            gene: row for row, gene in enumerate(features.genes.astype(str).tolist())
        }
        pairs: list[str] = []
        unsupported: list[str] = []
        for query in map(str, eligible_targets):
            query_row = feature_index.get(query)
            if (
                query_row is None
                or float(np.linalg.norm(features.matrix[query_row])) <= 0
            ):
                unsupported.append(query)
                continue
            similarity = source_matrix @ features.matrix[query_row]
            similarity[~source_available] = -np.inf
            similarity[ordered_sources == query] = -np.inf
            weights = RUN.neighbour_weights(similarity)
            chosen = np.flatnonzero(weights > 0)
            chosen = chosen[np.argsort(-weights[chosen], kind="stable")]
            for index in chosen:
                target = str(ordered_sources[index])
                selected_union.add(target)
                pairs.append(f"{query}\t{target}")
        arm_audit[arm] = {
            "query_neighbor_pair_count": len(pairs),
            "query_neighbor_pairs_ordered_sha256": RUN.sha256_lines(pairs),
            "unsupported_queries": unsupported,
            "unsupported_queries_ordered_sha256": RUN.sha256_lines(unsupported),
        }
    selected = tuple(sorted(selected_union))
    if not selected:
        raise ContractError("source neighbor closure is empty")
    return selected, {
        "full_source_target_count": len(ordered_sources),
        "full_source_target_ordered_sha256": RUN.sha256_lines(ordered_sources),
        "selected_target_count": len(selected),
        "selected_target_ordered_sha256": RUN.sha256_lines(selected),
        "arms": arm_audit,
    }


def _operator_apply(
    matrix: np.ndarray,
    operator: RoleOperator,
    *,
    gamma: float = 0.65,
    lam: float = 0.0,
) -> np.ndarray:
    values = np.asarray(matrix, dtype=np.float64)
    vectors = operator.evecs
    weights = (operator.evals + lam) ** (-gamma / 2.0)
    return ((values @ vectors) * weights) @ vectors.T


def build_source_bank_arrays(
    *,
    donors: Mapping[str, RoleDonor],
    allowed_source_ids: Sequence[str],
    canonical_genes: np.ndarray,
    canonical_targets: np.ndarray,
    held_targets: Sequence[str],
    operator: RoleOperator,
    selected_targets: Sequence[str],
    staging: Path,
    chunk_rows: int,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    allowed = tuple(map(str, allowed_source_ids))
    if tuple(donors) != allowed:
        raise ContractError("source donor order differs from the role allowlist")
    canonical_index = {gene: index for index, gene in enumerate(canonical_genes)}
    panel_positions = np.asarray(
        [canonical_index[target] for target in canonical_targets], dtype=np.int64
    )
    held_positions = np.asarray(
        sorted(
            canonical_index[target]
            for target in held_targets
            if target in canonical_index
        ),
        dtype=np.int64,
    )
    if np.intersect1d(operator.columns, panel_positions).size:
        raise ContractError("role operator contains a canonical panel coordinate")
    if np.intersect1d(operator.columns, held_positions).size:
        raise ContractError("role operator contains a held-surface coordinate")
    selected = tuple(sorted(set(map(str, selected_targets))))
    row_pairs = [
        (source_id, target)
        for source_id in allowed
        for target in selected
        if target in donors[source_id].pert_index
    ]
    if not row_pairs:
        raise ContractError("source-bank neighbor closure has no source rows")
    observed_ids = tuple(dict.fromkeys(source_id for source_id, _ in row_pairs))
    if observed_ids != allowed:
        raise ContractError(
            "not every allowlisted source contributes a source-bank row"
        )
    max_target = max(map(len, (target for _, target in row_pairs)))
    max_source = max(map(len, allowed))
    targets = np.asarray([target for _, target in row_pairs], dtype=f"U{max_target}")
    source_ids = np.asarray(
        [source_id for source_id, _ in row_pairs], dtype=f"U{max_source}"
    )
    genes = np.asarray(canonical_genes).astype(str)
    shape = (len(row_pairs), len(genes))
    effects = np.lib.format.open_memmap(
        staging / ".source-effects.npy",
        mode="w+",
        dtype=np.float64,
        shape=shape,
    )
    covered = np.lib.format.open_memmap(
        staging / ".source-covered.npy",
        mode="w+",
        dtype=bool,
        shape=shape,
    )
    effects[:] = 0.0
    covered[:] = False
    covered[:, operator.columns] = True
    operator_position = {
        int(column): index for index, column in enumerate(operator.columns)
    }
    row_cursor = 0
    per_source_rows: dict[str, int] = {}
    reference_scale = float(np.median([donor.scale for donor in donors.values()]))
    for source_id in allowed:
        donor = donors[source_id]
        donor_targets = [
            target
            for candidate_source, target in row_pairs
            if candidate_source == source_id
        ]
        per_source_rows[source_id] = len(donor_targets)
        for start in range(0, len(donor_targets), chunk_rows):
            stop = min(start + chunk_rows, len(donor_targets))
            batch_targets = donor_targets[start:stop]
            pre = np.zeros(
                (len(batch_targets), len(operator.columns)), dtype=np.float64
            )
            shared_target: list[int] = []
            for target_column in donor.tcols_unique:
                position = operator_position.get(int(target_column))
                if position is not None:
                    shared_target.append(position)
            target_columns = np.asarray(shared_target, dtype=np.int64)
            for row, target in enumerate(batch_targets):
                source_row = donor.pert_index[target]
                scaled = role_scaled_source_row(
                    donor,
                    source_row,
                    len(genes),
                    reference_scale,
                )
                pre[row, target_columns] = scaled[operator.columns[target_columns]]
            transformed = _operator_apply(pre, operator)
            if not np.isfinite(transformed).all():
                raise ContractError(f"{source_id}: nonfinite whitened source row")
            destination = slice(row_cursor + start, row_cursor + stop)
            full = np.zeros((len(batch_targets), len(genes)), dtype=np.float64)
            full[:, operator.columns] = transformed
            effects[destination] = full
        row_cursor += len(donor_targets)
    effects.flush()
    covered.flush()
    if np.any(effects[:, panel_positions] != 0):
        raise ContractError("source bank retained canonical panel coordinates")
    if np.any(covered[:, panel_positions]):
        raise ContractError("source bank covered canonical panel coordinates")
    if held_positions.size and (
        np.any(effects[:, held_positions] != 0) or np.any(covered[:, held_positions])
    ):
        raise ContractError("source bank retained held-surface coordinates")
    arrays: dict[str, np.ndarray] = {
        "targets": targets,
        "genes": genes,
        "effects": effects,
        "covered": covered,
        "source_ids": source_ids,
    }
    coverage_mask = np.zeros(len(genes), dtype=bool)
    coverage_mask[operator.columns] = True
    return arrays, {
        "row_count": len(row_pairs),
        "rows_by_source_id": per_source_rows,
        "target_axis_ordered_sha256": RUN.sha256_lines(targets),
        "unique_target_ordered_sha256": RUN.sha256_lines(sorted(set(targets.tolist()))),
        "coverage_mask_sha256": sha256_array_stream(coverage_mask),
        "covered_coordinates_per_row": int(coverage_mask.sum()),
        "panel_coordinates_zero_and_uncovered": True,
        "held_surface_coordinates_zero_and_uncovered": True,
    }


def derive_exact_eligibility_arrays(
    manifest: Mapping[str, Any],
    canonical_targets: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    donor_records = manifest["source_inputs"]["donors"]
    donor_sets = {
        name: set(
            read_axis(
                _record_path(record["perts"], f"{name} eligibility targets"),
                unique=True,
                label=f"{name} eligibility targets",
            )
        )
        for name, record in donor_records.items()
    }
    perturbmatch = set(
        read_axis(
            _record_path(
                manifest["source_inputs"]["perturbmatch"]["perts"],
                "PerturbMatch eligibility targets",
            ),
            unique=True,
            label="PerturbMatch eligibility targets",
        )
    )
    carriers = np.asarray(
        [
            sum(target in donor_sets[name] for name in donor_sets)
            for target in canonical_targets
        ],
        dtype=np.int16,
    )
    fallback = np.asarray(
        [
            target not in donor_sets["k562_gwps"] and target not in perturbmatch
            for target in canonical_targets
        ],
        dtype=bool,
    )
    expected = np.full(len(canonical_targets), 3, dtype=np.int16)
    low = (carriers < expected) & ~fallback
    eligible = fallback | low
    frozen = manifest["eligibility"]
    if int(fallback.sum()) != frozen["incumbent_uncovered_count"]:
        raise ContractError("incumbent fallback count differs from frozen eligibility")
    if (
        RUN.sha256_lines(canonical_targets[fallback])
        != (frozen["incumbent_uncovered_ordered_sha256"])
    ):
        raise ContractError("incumbent fallback axis differs from frozen eligibility")
    if int(low.sum()) != frozen["low_confidence_count"]:
        raise ContractError("low-confidence count differs from frozen eligibility")
    if (
        RUN.sha256_lines(canonical_targets[low])
        != (frozen["low_confidence_ordered_sha256"])
    ):
        raise ContractError("low-confidence axis differs from frozen eligibility")
    if int(eligible.sum()) != frozen["eligible_count"]:
        raise ContractError("eligible count differs from frozen eligibility")
    if (
        RUN.sha256_lines(canonical_targets[eligible])
        != (frozen["eligible_ordered_sha256"])
    ):
        raise ContractError("eligible axis differs from frozen eligibility")
    return fallback, carriers, expected


def build_incumbent_arrays(
    *,
    manifest: Mapping[str, Any],
    role_spec: Mapping[str, Any],
    donors: Mapping[str, RoleDonor],
    canonical_targets: np.ndarray,
    canonical_genes: np.ndarray,
    operator: RoleOperator,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    production = _load_production_module()
    model = StreamingA14Model(donors, canonical_genes)
    weights = {name: 1.0 for name in donors}
    baseline = np.vstack(
        [
            model.effect(target, 0.982, 0.0, weights=weights)
            for target in canonical_targets
        ]
    ).astype(np.float64)
    perturbmatch_root = ROOT / "data" / "external" / "perturbmatch"
    pm = production.EXP061.prepare_perturbmatch(
        perturbmatch_root,
        tuple(map(str, canonical_genes)),
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    candidate, pm_audit = production.apply_perturbmatch_k562_replacement(
        baseline,
        model,
        pm,
        canonical_targets,
        alpha=0.982,
        weights=weights,
        blend=0.5,
    )
    canonical_index = {gene: index for index, gene in enumerate(canonical_genes)}
    excluded_symbols = set(map(str, canonical_targets)) | set(operator.held_targets)
    excluded_positions = np.asarray(
        sorted(
            canonical_index[symbol]
            for symbol in excluded_symbols
            if symbol in canonical_index
        ),
        dtype=np.int64,
    )
    candidate[:, excluded_positions] = 0.0
    whitener = production.Whitener(
        operator.evals,
        operator.evecs,
    )
    whitened, whitening_audit = production.apply_v9_whitening(
        candidate,
        whitener,
        operator.columns,
        gamma=0.65,
        lam=0.0,
    )
    _, cd4 = production.prepare_frozen_cd4(
        ROOT / "data" / "external" / "cd4t" / "cloud-compact" / "cd4_compact",
        canonical_targets,
        canonical_genes,
        tuple(sorted(excluded_symbols)),
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    final, cd4_audit = production.apply_frozen_cd4_residual(
        whitened,
        cd4,
        operator.columns,
        excluded_positions,
    )
    final[:, excluded_positions] = 0.0
    if final.shape != (300, 18_533) or not np.isfinite(final).all():
        raise ContractError("canonical incumbent matrix is invalid")
    fallback, carriers, expected = derive_exact_eligibility_arrays(
        manifest, canonical_targets
    )
    arrays = {
        "targets": canonical_targets.astype(str),
        "genes": canonical_genes.astype(str),
        "effects": np.asarray(final, dtype=np.float64),
        "incumbent_fallback": fallback,
        "source_carrier_count": carriers,
        "expected_source_count": expected,
    }
    return arrays, {
        "a14_parameters": manifest["a14"],
        "role_donor_order": list(donors),
        "role_donor_scales": {name: donors[name].scale for name in donors},
        "role_reference_scale": float(np.median([d.scale for d in donors.values()])),
        "excluded_coordinate_count": len(excluded_positions),
        "excluded_coordinate_indices_sha256": sha256_int64(excluded_positions),
        "held_target_count": len(operator.held_targets),
        "held_target_ordered_sha256": RUN.sha256_lines(operator.held_targets),
        "perturbmatch_audit": pm_audit,
        "whitening_audit": whitening_audit,
        "cd4_audit": cd4_audit,
    }


def deterministic_npz(path: Path, arrays: Mapping[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.is_symlink():
        raise ContractError(f"refusing to overwrite NPZ artifact: {path}")
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("xb") as raw:
            with zipfile.ZipFile(raw, mode="w", allowZip64=True) as archive:
                for name, value in arrays.items():
                    array = np.asarray(value)
                    if array.dtype.hasobject or not array.flags.c_contiguous:
                        raise ContractError(
                            f"{path.name}:{name} must be non-object C-order"
                        )
                    info = zipfile.ZipInfo(
                        filename=f"{name}.npy",
                        date_time=FIXED_ZIP_TIMESTAMP,
                    )
                    info.compress_type = zipfile.ZIP_STORED
                    info.create_system = 3
                    info.external_attr = 0o600 << 16
                    with archive.open(info, mode="w", force_zip64=True) as member:
                        np.lib.format.write_array(
                            member,
                            array,
                            allow_pickle=False,
                        )
            raw.flush()
            os.fsync(raw.fileno())
        os.link(temporary, path)
        temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def deterministic_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.is_symlink():
        raise ContractError(f"refusing to overwrite JSON artifact: {path}")
    payload = (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False)
        + "\n"
    ).encode("utf-8")
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("xb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, path)
        temporary.unlink()
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _source_records_for_role(
    manifest: Mapping[str, Any],
    role_spec: Mapping[str, Any],
    *,
    include_incumbent_sources: bool,
    operator: RoleOperator,
    hash_cache: dict[Path, str] | None = None,
) -> list[dict[str, Any]]:
    paths: set[Path] = {
        ROLE_SOURCE_MANIFEST_PATH,
        Path(__file__).resolve(),
        RUN_PATH,
        RUN.PROTOCOL_PATH,
        RUN.SOURCE_MANIFEST_PATH,
        RUN.SOURCE_CANDIDATE_FREEZE_BOUNDARY_PATH,
        operator.artifact_path,
        operator.manifest_path,
    }
    for record in _iter_records(manifest["canonical_axes"]):
        paths.add(
            verify_record(
                record,
                "canonical axis input",
                hash_cache=hash_cache,
            )
        )
    for record in _iter_records(manifest["upstream_contracts"]):
        paths.add(
            verify_record(
                record,
                "upstream contract input",
                hash_cache=hash_cache,
            )
        )
    for relative, digest in manifest["implementation_inputs"].items():
        paths.add(
            verify_record(
                {"path": relative, "sha256": digest},
                f"implementation input {relative}",
                hash_cache=hash_cache,
            )
        )
    source_inputs = manifest["source_inputs"]
    for key in ("challenge_manifest", "eligibility_cross_check"):
        paths.add(verify_record(source_inputs[key], key, hash_cache=hash_cache))
    for key in ("depmap_features", "depmap_feature_receipt", "string_features"):
        paths.add(verify_record(source_inputs[key], key, hash_cache=hash_cache))
    for donor_name in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        donor = source_inputs["donors"][donor_name]
        if donor_name in role_spec["allowed_source_ids"]:
            for record in donor.values():
                paths.add(
                    verify_record(
                        record,
                        f"{donor_name} input",
                        hash_cache=hash_cache,
                    )
                )
        else:
            paths.add(
                verify_record(
                    donor["perts"],
                    f"{donor_name} eligibility-only target axis",
                    hash_cache=hash_cache,
                )
            )
    paths.add(
        verify_record(
            source_inputs["perturbmatch"]["perts"],
            "PM target axis",
            hash_cache=hash_cache,
        )
    )
    if include_incumbent_sources:
        for record in source_inputs["perturbmatch"].values():
            paths.add(
                verify_record(
                    record,
                    "PerturbMatch input",
                    hash_cache=hash_cache,
                )
            )
        for condition in source_inputs["cd4"].values():
            for record in condition.values():
                paths.add(verify_record(record, "CD4 input", hash_cache=hash_cache))
    forbidden_effect = (
        source_inputs["donors"]
        .get(role_spec["held_surface_source_id"], {})
        .get("delta")
    )
    if forbidden_effect is not None:
        forbidden_path = _record_path(forbidden_effect, "held donor effect")
        if forbidden_path in paths:
            raise ContractError("role source inputs include held-surface donor effects")
    return [
        receipt_record(path, hash_cache=hash_cache) for path in sorted(paths, key=str)
    ]


def seal_role_bundle(
    *,
    role: str,
    role_spec: Mapping[str, Any],
    incumbent_arrays: Mapping[str, np.ndarray],
    source_arrays: Mapping[str, np.ndarray],
    incumbent_audit: Mapping[str, Any],
    source_audit: Mapping[str, Any],
    neighbor_audit: Mapping[str, Any],
    operator: RoleOperator,
    manifest: Mapping[str, Any],
    output_directory: Path | None = None,
    incumbent_source_inputs: Sequence[Mapping[str, Any]] | None = None,
    source_bank_source_inputs: Sequence[Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    expected_output = RUN.resolve_repository_path(
        role_spec["output_directory"], f"{role} output directory"
    )
    output = expected_output if output_directory is None else output_directory.resolve()
    if output != expected_output:
        raise ContractError(f"{role}: noncanonical output directory")
    RUN.require_under_experiment(output)
    RUN.reject_symlink_components(output, f"{role} output")
    if output.exists() or output.is_symlink():
        raise ContractError(f"refusing to overwrite role directory: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{role}.", suffix=".tmp", dir=output.parent)
    )
    try:
        incumbent_path = staging / "incumbent.npz"
        source_path = staging / "source-bank.npz"
        deterministic_npz(incumbent_path, incumbent_arrays)
        deterministic_npz(source_path, source_arrays)
        builder = receipt_record(Path(__file__).resolve())
        incumbent_receipt = {
            "schema_version": 1,
            "experiment": "114-depmap-coessentiality",
            "stage": "role_incumbent_source_only",
            "role": role,
            "surface": role_spec["surface"],
            "source_only": True,
            "contains_held_truth": False,
            "sealed_before_held_truth": True,
            "artifact_path": (output / "incumbent.npz")
            .relative_to(ROOT.resolve())
            .as_posix(),
            "artifact_sha256": RUN.sha256_file(incumbent_path),
            "builder": builder,
            "role_source_manifest": binding_record(ROLE_SOURCE_MANIFEST_PATH),
            "source_inputs": list(incumbent_source_inputs)
            if incumbent_source_inputs is not None
            else _source_records_for_role(
                manifest,
                role_spec,
                include_incumbent_sources=True,
                operator=operator,
            ),
            "endpoint": role_spec["incumbent_endpoint"],
            "parameters": incumbent_audit,
            "axes": {
                "target_count": len(incumbent_arrays["targets"]),
                "target_ordered_sha256": RUN.sha256_lines(incumbent_arrays["targets"]),
                "gene_count": len(incumbent_arrays["genes"]),
                "gene_ordered_sha256": RUN.sha256_lines(incumbent_arrays["genes"]),
            },
            "array_sha256": {
                name: sha256_array_stream(np.asarray(incumbent_arrays[name]))
                for name in (
                    "targets",
                    "genes",
                    "effects",
                    "incumbent_fallback",
                    "source_carrier_count",
                    "expected_source_count",
                )
            },
        }
        source_receipt = {
            "schema_version": 1,
            "experiment": "114-depmap-coessentiality",
            "stage": "role_source_bank_source_only",
            "role": role,
            "surface": role_spec["surface"],
            "source_only": True,
            "contains_held_truth": False,
            "sealed_before_held_truth": True,
            "artifact_path": (output / "source-bank.npz")
            .relative_to(ROOT.resolve())
            .as_posix(),
            "artifact_sha256": RUN.sha256_file(source_path),
            "builder": builder,
            "role_source_manifest": binding_record(ROLE_SOURCE_MANIFEST_PATH),
            "source_inputs": list(source_bank_source_inputs)
            if source_bank_source_inputs is not None
            else _source_records_for_role(
                manifest,
                role_spec,
                include_incumbent_sources=False,
                operator=operator,
            ),
            "effect_space": role_spec["source_effect_space"],
            "source_ids": role_spec["allowed_source_ids"],
            "forbidden_source_ids": role_spec["forbidden_source_ids"],
            "held_surface_source_id": role_spec["held_surface_source_id"],
            "parameters": {
                "a14": manifest["a14"],
                "source_bank": manifest["source_bank"],
                "operator": {
                    "artifact_path": operator.artifact_path.relative_to(
                        ROOT.resolve()
                    ).as_posix(),
                    "artifact_sha256": operator.artifact_sha256,
                    "manifest_path": operator.manifest_path.relative_to(
                        ROOT.resolve()
                    ).as_posix(),
                    "manifest_sha256": operator.manifest_sha256,
                    "columns_int64_sha256": sha256_int64(operator.columns),
                    "held_target_ordered_sha256": RUN.sha256_lines(
                        operator.held_targets
                    ),
                },
                "neighbor_closure": neighbor_audit,
                "artifact": source_audit,
            },
            "axes": {
                "row_count": len(source_arrays["targets"]),
                "target_ordered_sha256": RUN.sha256_lines(source_arrays["targets"]),
                "gene_count": len(source_arrays["genes"]),
                "gene_ordered_sha256": RUN.sha256_lines(source_arrays["genes"]),
                "coverage_sha256": sha256_array_stream(source_arrays["covered"]),
            },
            "array_sha256": {
                name: sha256_array_stream(np.asarray(source_arrays[name]))
                for name in ("targets", "genes", "effects", "covered", "source_ids")
            },
        }
        incumbent_receipt_path = staging / "incumbent-receipt.json"
        source_receipt_path = staging / "source-bank-receipt.json"
        deterministic_json(incumbent_receipt_path, incumbent_receipt)
        deterministic_json(source_receipt_path, source_receipt)

        def output_record(staged: Path, final_name: str) -> dict[str, str]:
            return {
                "path": (output / final_name).relative_to(ROOT.resolve()).as_posix(),
                "sha256": RUN.sha256_file(staged),
            }

        role_manifest = {
            "schema_version": 1,
            "experiment": "114-depmap-coessentiality",
            "stage": "source_only_role_inputs",
            "role": role,
            "surface": role_spec["surface"],
            "role_directory": output.relative_to(ROOT.resolve()).as_posix(),
            "contains_held_truth": False,
            "source_only": True,
            "sealed_before_held_truth": True,
            "role_source_manifest": binding_record(ROLE_SOURCE_MANIFEST_PATH),
            "incumbent": {
                "artifact": output_record(incumbent_path, "incumbent.npz"),
                "receipt": output_record(
                    incumbent_receipt_path, "incumbent-receipt.json"
                ),
            },
            "source_bank": {
                "artifact": output_record(source_path, "source-bank.npz"),
                "receipt": output_record(
                    source_receipt_path, "source-bank-receipt.json"
                ),
            },
        }
        deterministic_json(staging / "manifest.json", role_manifest)
        for path in staging.iterdir():
            if path.is_file() and not path.name.startswith(".source-"):
                path.chmod(0o444)
        for scratch in (
            staging / ".source-effects.npy",
            staging / ".source-covered.npy",
        ):
            if scratch.exists():
                scratch.unlink()
        os.rename(staging, output)
        return role_manifest
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def _verify_manifest_inputs(manifest: Mapping[str, Any]) -> None:
    for index, record in enumerate(_iter_records(manifest["canonical_axes"])):
        verify_record(record, f"canonical input {index}")
    for index, record in enumerate(_iter_records(manifest["upstream_contracts"])):
        verify_record(record, f"upstream input {index}")
    for relative, digest in manifest["implementation_inputs"].items():
        verify_record(
            {"path": relative, "sha256": digest},
            f"implementation input {relative}",
        )


def build_role(
    role: str,
    *,
    confirmed: bool,
    chunk_rows: int,
    max_rss_gib: float,
) -> dict[str, Any]:
    if not confirmed:
        raise ContractError("real role build requires --confirm-source-only-real-build")
    active = active_exp107_build_pids()
    if active:
        raise ContractError(
            "Experiment 107 is still building; run role builds serially after PIDs "
            + ",".join(map(str, active))
            + " exit"
        )
    validate_checked_in_locks()
    manifest = load_role_source_manifest()
    if role not in manifest["roles"]:
        raise ContractError(f"unknown role: {role}")
    _verify_manifest_inputs(manifest)
    role_spec = manifest["roles"][role]
    operator = load_role_operator(role, role_spec)
    hash_cache: dict[Path, str] = {
        operator.artifact_path.resolve(): operator.artifact_sha256,
        operator.manifest_path.resolve(): operator.manifest_sha256,
    }
    incumbent_source_inputs = _source_records_for_role(
        manifest,
        role_spec,
        include_incumbent_sources=True,
        operator=operator,
        hash_cache=hash_cache,
    )
    source_bank_source_inputs = _source_records_for_role(
        manifest,
        role_spec,
        include_incumbent_sources=False,
        operator=operator,
        hash_cache=hash_cache,
    )
    target_record = manifest["canonical_axes"]["targets"]
    gene_record = manifest["canonical_axes"]["genes"]
    target_path = verify_axis_record(target_record, "canonical target axis")
    gene_path = verify_axis_record(gene_record, "canonical gene axis")
    canonical_targets = np.asarray(
        RUN.read_csv_column(target_path, target_record["column"]), dtype=str
    )
    canonical_genes = np.asarray(
        RUN.read_csv_column(gene_path, gene_record["column"]), dtype=str
    )
    if (
        len(canonical_targets) != 300
        or RUN.sha256_lines(canonical_targets) != target_record["ordered_sha256"]
        or len(canonical_genes) != 18_533
        or RUN.sha256_lines(canonical_genes) != gene_record["ordered_sha256"]
    ):
        raise ContractError("canonical 300x18533 axes drifted")
    donor_records = manifest["source_inputs"]["donors"]
    source_target_axes = {
        source_id: read_axis(
            _record_path(
                donor_records[source_id]["perts"],
                f"{source_id} source target axis",
            ),
            unique=True,
            label=f"{source_id} source target axis",
        )
        for source_id in role_spec["allowed_source_ids"]
    }
    boundary = RUN.load_json(RUN.SOURCE_CANDIDATE_FREEZE_BOUNDARY_PATH)
    source_manifest = RUN.load_json(RUN.SOURCE_MANIFEST_PATH)
    protocol = RUN.load_json(RUN.PROTOCOL_PATH)
    depmap = RUN.validate_depmap_feature_artifact(boundary, protocol, source_manifest)
    string = RUN.validate_string_feature_artifact(boundary, source_manifest)
    blend = RUN.blend_features(string, depmap, 1.0)
    eligible = np.asarray(protocol["eligibility"]["eligible_targets"], dtype=str)
    selected_targets, neighbor_audit = source_target_closure(
        [
            target
            for source_id in role_spec["allowed_source_ids"]
            for target in source_target_axes[source_id]
        ],
        eligible,
        {
            "string_neighbor": string,
            "depmap_neighbor": depmap,
            "string_depmap_blend_w1": blend,
        },
    )
    del depmap, string, blend
    gc.collect()
    donors: dict[str, RoleDonor] = {}
    try:
        for source_id in role_spec["allowed_source_ids"]:
            donors[source_id] = RoleDonor.load(
                source_id,
                donor_records[source_id],
                canonical_genes,
                chunk_rows=chunk_rows,
                hash_cache=hash_cache,
            )
        incumbent_arrays, incumbent_audit = build_incumbent_arrays(
            manifest=manifest,
            role_spec=role_spec,
            donors=donors,
            canonical_targets=canonical_targets,
            canonical_genes=canonical_genes,
            operator=operator,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        expected_output = RUN.resolve_repository_path(
            role_spec["output_directory"], f"{role} output"
        )
        scratch_parent = expected_output.parent
        scratch_parent.mkdir(parents=True, exist_ok=True)
        scratch = Path(tempfile.mkdtemp(prefix=f".{role}.arrays.", dir=scratch_parent))
        try:
            source_arrays, source_audit = build_source_bank_arrays(
                donors=donors,
                allowed_source_ids=role_spec["allowed_source_ids"],
                canonical_genes=canonical_genes,
                canonical_targets=canonical_targets,
                held_targets=operator.held_targets,
                operator=operator,
                selected_targets=selected_targets,
                staging=scratch,
                chunk_rows=chunk_rows,
            )
            result = seal_role_bundle(
                role=role,
                role_spec=role_spec,
                incumbent_arrays=incumbent_arrays,
                source_arrays=source_arrays,
                incumbent_audit=incumbent_audit,
                source_audit=source_audit,
                neighbor_audit=neighbor_audit,
                operator=operator,
                manifest=manifest,
                incumbent_source_inputs=incumbent_source_inputs,
                source_bank_source_inputs=source_bank_source_inputs,
            )
        finally:
            shutil.rmtree(scratch, ignore_errors=True)
        return result
    finally:
        for donor in donors.values():
            donor.close()
        gc.collect()


def metadata_report() -> dict[str, Any]:
    validate_checked_in_locks()
    manifest = load_role_source_manifest()
    roles: dict[str, Any] = {}
    for role, record in manifest["roles"].items():
        operator = record["operator"]
        artifact = RUN.resolve_repository_path(
            operator["artifact_path"], f"{role} operator"
        )
        upstream_manifest = RUN.resolve_repository_path(
            operator["manifest_path"], f"{role} operator manifest"
        )
        output = RUN.resolve_repository_path(
            record["output_directory"], f"{role} output"
        )
        roles[role] = {
            "allowed_source_ids": record["allowed_source_ids"],
            "forbidden_source_ids": record["forbidden_source_ids"],
            "held_surface_source_id": record["held_surface_source_id"],
            "operator_artifact_ready": artifact.is_file() and not artifact.is_symlink(),
            "operator_manifest_ready": (
                upstream_manifest.is_file() and not upstream_manifest.is_symlink()
            ),
            "role_bundle_exists": output.exists(),
            "future_command": (
                "uv run python experiments/114-depmap-coessentiality/"
                f"build_role_inputs.py build-role --role {role} "
                "--confirm-source-only-real-build"
            ),
        }
    return {
        "experiment": "114-depmap-coessentiality",
        "role_source_manifest_sha256": RUN.sha256_file(ROLE_SOURCE_MANIFEST_PATH),
        "canonical_shape": [300, 18_533],
        "eligible_count": 41,
        "source_candidate_freeze": False,
        "approved_roles": {},
        "active_experiment_107_build_pids": list(active_exp107_build_pids()),
        "roles": roles,
    }


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    commands.add_parser(
        "check",
        help="lightweight metadata and lock check; does not open source arrays",
    )
    build = commands.add_parser(
        "build-role",
        help="build one immutable real role bundle after Experiment 107 finishes",
    )
    build.add_argument("--role", choices=EXPECTED_ROLES, required=True)
    build.add_argument("--confirm-source-only-real-build", action="store_true")
    build.add_argument("--chunk-rows", type=int, default=64)
    build.add_argument("--max-rss-gib", type=float, default=24.0)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    if args.command == "check":
        print(json.dumps(metadata_report(), indent=2, sort_keys=True))
        return 0
    if args.chunk_rows <= 0:
        raise SystemExit("--chunk-rows must be positive")
    result = build_role(
        args.role,
        confirmed=args.confirm_source_only_real_build,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

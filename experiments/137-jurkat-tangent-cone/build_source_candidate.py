#!/usr/bin/env python3
"""Build the frozen Experiment 137 source candidate without Nadig truth access."""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.util
import json
import os
import sys
import tempfile
from pathlib import Path
from types import ModuleType
from typing import Any, Mapping

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXP107 = ROOT / "experiments/107-jurkat-high-moi"
METADATA_NPZ = HERE / "work/metadata-freeze-v1.npz"
METADATA_JSON = HERE / "work/metadata-freeze-v1.json"
PROTOCOL = HERE / "preregistered_protocol.json"
AUTHORIZATION = HERE / "PREVALUE_AUTHORIZATION.v1.json"
TANGENT_PATH = HERE / "tangent.py"
BUILDER_PATH = EXP107 / "strict_artifact_builder.py"
STRICT_PATH = EXP107 / "strict_nadig_validation.py"
SOURCE_DIRECTORY = EXP107 / "work/frozen-source-model"
OUTPUT_ROOT = HERE / "work/source-candidate-v1"
OPERATOR_NPZ = OUTPUT_ROOT / "operator.npz"
OPERATOR_JSON = OUTPUT_ROOT / "operator.json"
CANDIDATE_NPZ = OUTPUT_ROOT / "candidate.npz"
CANDIDATE_JSON = OUTPUT_ROOT / "candidate.json"
CHALLENGE = ROOT / "data/challenge"
DONOR_ROOT = ROOT / "data/donors"
PERTURBMATCH = ROOT / "data/external/perturbmatch"
CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"

EXPECTED_METADATA_NPZ_SHA256 = (
    "682afb9a32a2f2732486da04cf27359abf2d4b91c1ba92627b8d05f45dce9d1c"
)
EXPECTED_METADATA_JSON_SHA256 = (
    "45de70d1e9c687d8ab8967946f9c7f46a4982d41b6c107f9753905a0fe0ba1a0"
)
EXPECTED_PROTOCOL_SHA256 = (
    "f1f085296213cfcbfab93212cc963982f34b683e5ce6e7945edf5bee65ea65f7"
)
EXPECTED_SOURCE_MODEL_SHA256 = (
    "03759a7b2cafa86d2643a58669235aca7fc38a33599aa403c67439f8979a7e26"
)
DONOR_ORDER = ("k562_gwps", "orion_hct116", "orion_hek293t")


class BuildError(RuntimeError):
    """Raised when the frozen source-candidate contract is violated."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise BuildError(message)


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_bytes), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_array(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values))
    require(not array.dtype.hasobject, "object arrays cannot be hashed")
    digest = hashlib.sha256()
    digest.update(b"exp137-array-v1\0")
    digest.update(array.dtype.str.encode("ascii"))
    digest.update(b"\0")
    digest.update(json.dumps(array.shape).encode("ascii"))
    digest.update(b"\0")
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def atomic_savez(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".npz"
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        np.savez(temporary, **arrays)
        if path.exists():
            require(path.read_bytes() == temporary.read_bytes(), f"refusing overwrite: {path}")
            temporary.unlink()
        else:
            os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    encoded = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        require(path.read_bytes() == encoded, f"refusing overwrite: {path}")
        return
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp"
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(encoded)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def cap_key(cap: float) -> str:
    return f"cap_{int(round(float(cap) * 100)):04d}"


def file_binding(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(ROOT)),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def validate_authorization() -> dict[str, Any]:
    require(AUTHORIZATION.is_file(), "independent pre-value authorization is missing")
    payload = read_json(AUTHORIZATION)
    require(payload.get("schema_version") == 1, "authorization schema changed")
    require(payload.get("experiment") == "137-jurkat-tangent-cone", "wrong authorization")
    require(payload.get("status") == "approved_before_effect_access", "not authorized")
    bindings = payload.get("bindings")
    require(isinstance(bindings, dict), "authorization bindings are missing")
    expected_paths = {
        "metadata_npz": METADATA_NPZ,
        "metadata_json": METADATA_JSON,
        "protocol": PROTOCOL,
        "builder": Path(__file__).resolve(),
        "evaluator": HERE / "evaluate.py",
        "tangent": TANGENT_PATH,
        "builder_test": HERE / "test_build_source_candidate.py",
        "evaluator_test": HERE / "test_evaluate.py",
        "tangent_test": HERE / "test_tangent.py",
        "metadata_test": HERE / "test_freeze_metadata.py",
    }
    require(set(bindings) == set(expected_paths), "authorization binding keys changed")
    for role, path in expected_paths.items():
        record = bindings[role]
        require(record.get("path") == str(path.relative_to(ROOT)), f"{role} path changed")
        require(record.get("sha256") == sha256_file(path), f"{role} hash changed")
    return payload


def operator_report(
    plan: Any, arrays: Mapping[str, np.ndarray], builder: ModuleType
) -> dict[str, Any]:
    inputs = [
        METADATA_NPZ,
        METADATA_JSON,
        PROTOCOL,
        AUTHORIZATION,
        Path(__file__).resolve(),
        BUILDER_PATH,
        STRICT_PATH,
        ROOT / "scripts/build_union_whitener.py",
        ROOT / "src/vcc2026/whiten.py",
        CHALLENGE / "gene_names.csv",
        CHALLENGE / "pert_counts.csv",
        ROOT / "data/donors/k562_essential/perts.txt",
    ]
    for donor in DONOR_ORDER:
        inputs.extend(
            DONOR_ROOT / donor / filename
            for filename in ("delta.npy", "genes.txt", "perts.txt")
        )
    return {
        "schema_version": 1,
        "experiment": "137-jurkat-tangent-cone",
        "status": "source_only_operator_built_without_nadig_truth",
        "artifact": {
            "path": str(OPERATOR_NPZ.relative_to(ROOT)),
            "bytes": OPERATOR_NPZ.stat().st_size,
            "sha256": sha256_file(OPERATOR_NPZ),
        },
        "arrays": {name: sha256_array(value) for name, value in arrays.items()},
        "parameters": {
            "donor_order": list(DONOR_ORDER),
            "target_row_exclusion_count": len(plan.targets),
            "coordinate_exclusion": "all 300 VCC plus all 394 fresh targets",
            "gamma": 0.65,
            "lambda": 0.0,
        },
        "transitive_inputs": [file_binding(path) for path in inputs],
        "experiment_107_builder_version": sha256_file(Path(builder.__file__)),
    }


def candidate_transitive_inputs(builder: ModuleType) -> list[dict[str, Any]]:
    paths = [
        METADATA_NPZ,
        METADATA_JSON,
        PROTOCOL,
        AUTHORIZATION,
        Path(__file__).resolve(),
        TANGENT_PATH,
        BUILDER_PATH,
        STRICT_PATH,
        SOURCE_DIRECTORY / "source-model.npz",
        SOURCE_DIRECTORY / "manifest.json",
        CHALLENGE / "gene_names.csv",
        CHALLENGE / "pert_counts.csv",
        *builder.PRODUCTION_IMPLEMENTATION_PATHS,
    ]
    for donor in DONOR_ORDER:
        paths.extend(
            DONOR_ROOT / donor / filename for filename in builder.DONOR_FILENAMES
        )
    paths.extend(PERTURBMATCH / filename for filename in builder.PERTURBMATCH_FILENAMES)
    for condition in builder.CD4_NAMES:
        paths.extend(
            CD4_ROOT / condition / filename for filename in builder.CD4_FILENAMES
        )
    unique = tuple(dict.fromkeys(Path(path).resolve() for path in paths))
    return [file_binding(path) for path in unique]


def validate_frozen_inputs() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    require(
        sha256_file(METADATA_NPZ) == EXPECTED_METADATA_NPZ_SHA256,
        "metadata NPZ changed",
    )
    require(
        sha256_file(METADATA_JSON) == EXPECTED_METADATA_JSON_SHA256,
        "metadata JSON changed",
    )
    require(
        sha256_file(PROTOCOL) == EXPECTED_PROTOCOL_SHA256,
        "preregistered protocol changed",
    )
    require(
        sha256_file(SOURCE_DIRECTORY / "source-model.npz")
        == EXPECTED_SOURCE_MODEL_SHA256,
        "source model changed",
    )
    with np.load(METADATA_NPZ, allow_pickle=False) as archive:
        metadata = {name: archive[name] for name in archive.files}
    report = read_json(METADATA_JSON)
    require(
        report.get("status") == "metadata_frozen_before_fresh_target_outcome_access",
        "metadata status changed",
    )
    require(
        np.asarray(metadata["fresh_targets"]).shape == (394,),
        "fresh target axis changed",
    )
    require(
        np.asarray(metadata["test_targets"]).shape == (79,),
        "test target axis changed",
    )
    require(
        np.asarray(metadata["development_targets"]).shape == (315,),
        "development target axis changed",
    )
    return metadata, report


def build_plan(metadata: Mapping[str, np.ndarray], builder: ModuleType) -> Any:
    return builder.SurfacePlan(
        key="nadig.fresh394",
        targets=tuple(np.asarray(metadata["fresh_targets"]).astype(str)),
        score_genes=tuple(np.asarray(metadata["score_genes"]).astype(str)),
        score_indices_in_source=np.asarray(
            metadata["score_indices_in_source"], dtype=np.int64
        ),
        score_indices_in_official=np.asarray(
            metadata["score_indices_in_official"], dtype=np.int64
        ),
        operator_columns=np.asarray(metadata["operator_columns"], dtype=np.int64),
        basis_rows={
            donor: np.asarray(metadata[f"basis_rows_{donor}"], dtype=np.int64)
            for donor in DONOR_ORDER
        },
        donor_order=DONOR_ORDER,
        source_target_indices=np.asarray(
            metadata["fresh_source_indices"], dtype=np.int64
        ),
        source_transferable=np.asarray(metadata["fresh_transferable"], dtype=bool),
    )


def build_operator(
    plan: Any,
    builder: ModuleType,
    *,
    max_rss_gib: float,
    chunk_rows: int,
) -> dict[str, np.ndarray]:
    if OPERATOR_NPZ.exists():
        require(OPERATOR_JSON.is_file(), "operator manifest is missing")
        with np.load(OPERATOR_NPZ, allow_pickle=False) as archive:
            arrays = {name: archive[name] for name in archive.files}
        expected_keys = {
            "evals",
            "evecs",
            "columns",
            *(f"basis_{donor}" for donor in DONOR_ORDER),
        }
        require(set(arrays) == expected_keys, "stored operator arrays changed")
        require(
            np.array_equal(arrays["columns"], plan.operator_columns),
            "stored operator columns changed",
        )
        for donor in DONOR_ORDER:
            require(
                np.array_equal(arrays[f"basis_{donor}"], plan.basis_rows[donor]),
                f"stored {donor} basis rows changed",
            )
        require(
            read_json(OPERATOR_JSON) == operator_report(plan, arrays, builder),
            "stored operator manifest or provenance changed",
        )
        return arrays
    arrays = builder._build_operator_arrays(  # noqa: SLF001
        plan,
        max_rss_gib=max_rss_gib,
        chunk_rows=chunk_rows,
        scratch_root=OUTPUT_ROOT / "scratch",
    )
    atomic_savez(OPERATOR_NPZ, **arrays)
    atomic_write_json(OPERATOR_JSON, operator_report(plan, arrays, builder))
    return arrays


def build_incumbent_and_source(
    plan: Any,
    metadata: Mapping[str, np.ndarray],
    builder: ModuleType,
    strict: ModuleType,
    operator: Mapping[str, np.ndarray],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> dict[str, np.ndarray]:
    production = builder._load_production_module()  # noqa: SLF001
    spec = production.ChallengeSpec.load(CHALLENGE)
    official = tuple(map(str, spec.genes))
    official_index = {symbol: index for index, symbol in enumerate(official)}
    donors = {
        name: production.Donor.load(DONOR_ROOT / name) for name in DONOR_ORDER
    }
    model = production.MultiDonorModel(
        donors, np.asarray(spec.genes), gene_norm=False, eb=False
    )
    weights = {name: 1.0 for name in DONOR_ORDER}
    perturbmatch = production.EXP061.prepare_perturbmatch(
        PERTURBMATCH,
        official,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    excluded = set(map(str, spec.perts)) | set(plan.targets)
    excluded_positions = np.asarray(
        sorted(official_index[symbol] for symbol in excluded if symbol in official_index),
        dtype=np.int64,
    )
    pre = production.build_v9_pre_whiten_effects(
        model, plan.targets, alpha=0.982, weights=weights
    )
    pre, _ = production.apply_perturbmatch_k562_replacement(
        pre,
        model,
        perturbmatch,
        plan.targets,
        alpha=0.982,
        weights=weights,
        blend=0.5,
    )
    pre[:, excluded_positions] = 0.0
    whitener = production.Whitener(
        np.asarray(operator["evals"], dtype=np.float64),
        np.asarray(operator["evecs"], dtype=np.float64),
    )
    columns = np.asarray(operator["columns"], dtype=np.int64)
    whitened, _ = production.apply_v9_whitening(
        pre, whitener, columns, gamma=0.65, lam=0.0
    )
    _, cd4 = production.prepare_frozen_cd4(
        CD4_ROOT,
        plan.targets,
        official,
        tuple(sorted(excluded)),
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    incumbent, _ = production.apply_frozen_cd4_residual(
        whitened, cd4, columns, excluded_positions
    )
    incumbent[:, excluded_positions] = 0.0
    require(np.all(np.isfinite(incumbent)), "incumbent is non-finite")

    source_model = strict.load_frozen_source_model(SOURCE_DIRECTORY)
    score_indices = np.asarray(metadata["score_indices_in_official"], dtype=np.int64)
    operator_supported, supported_positions, unsupported_official = (
        builder.partition_score_axis_by_operator(score_indices, columns)
    )
    raw_source = np.zeros(
        (len(plan.targets), int(operator_supported.sum())), dtype=np.float64
    )
    transferable = np.asarray(plan.source_transferable, dtype=bool).copy()
    for row, source_row in enumerate(plan.source_target_indices):
        if not transferable[row]:
            continue
        raw_source[row] = source_model.exported_effects[
            int(source_row),
            plan.score_indices_in_source[operator_supported],
        ]
    transformed = strict.apply_supported_principal_whitener(
        raw_source,
        np.asarray(operator["evals"], dtype=np.float64),
        np.asarray(operator["evecs"], dtype=np.float32),
        supported_positions,
        gamma=0.65,
        lam=0.0,
    )
    incumbent_score = np.asarray(incumbent[:, score_indices], dtype=np.float64)
    source_score = np.zeros_like(incumbent_score)
    source_score[:, operator_supported] = transformed
    coordinate_supported = np.asarray(operator_supported, dtype=bool)
    if unsupported_official.size:
        require(
            np.array_equal(
                incumbent[:, unsupported_official],
                np.zeros_like(incumbent[:, unsupported_official]),
            ),
            "operator-missing score coordinates have nonzero incumbent effects",
        )
    del pre, whitened, incumbent, raw_source, transformed, source_model
    gc.collect()
    return {
        "targets": np.asarray(plan.targets, dtype="U"),
        "test_targets": np.asarray(metadata["test_targets"]).astype("U"),
        "development_targets": np.asarray(metadata["development_targets"]).astype("U"),
        "development_fold_id": np.asarray(
            metadata["development_fold_id"], dtype=np.int8
        ),
        "score_genes": np.asarray(plan.score_genes, dtype="U"),
        "incumbent_score": incumbent_score,
        "source_score": source_score,
        "transferable": transferable,
        "coordinate_supported": coordinate_supported,
    }


def build_candidates(
    arrays: dict[str, np.ndarray], tangent: ModuleType
) -> dict[str, np.ndarray]:
    result = dict(arrays)
    for cap in tangent.CAP_DEGREES:
        candidate, effective, realized = tangent.spherical_tangent_candidate(
            arrays["incumbent_score"],
            arrays["source_score"],
            arrays["transferable"],
            arrays["coordinate_supported"],
            cap,
        )
        key = cap_key(cap)
        result[f"candidate_{key}"] = candidate
        result[f"effective_{key}"] = effective
        result[f"realized_degrees_{key}"] = realized
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-rss-gib", type=float, default=52.0)
    parser.add_argument("--chunk-rows", type=int, default=256)
    args = parser.parse_args(argv)

    validate_authorization()
    metadata, _ = validate_frozen_inputs()
    builder = load_module("exp137_artifact_builder", BUILDER_PATH)
    strict = load_module("exp137_strict", STRICT_PATH)
    tangent = load_module("exp137_tangent_runtime", TANGENT_PATH)
    plan = build_plan(metadata, builder)
    operator = build_operator(
        plan,
        builder,
        max_rss_gib=args.max_rss_gib,
        chunk_rows=args.chunk_rows,
    )
    core = build_incumbent_and_source(
        plan,
        metadata,
        builder,
        strict,
        operator,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    arrays = build_candidates(core, tangent)
    atomic_savez(CANDIDATE_NPZ, **arrays)
    report = {
        "schema_version": 1,
        "experiment": "137-jurkat-tangent-cone",
        "status": "source_candidate_built_without_nadig_truth",
        "artifact": {
            "path": str(CANDIDATE_NPZ.relative_to(ROOT)),
            "bytes": CANDIDATE_NPZ.stat().st_size,
            "sha256": sha256_file(CANDIDATE_NPZ),
        },
        "operator": {
            "path": str(OPERATOR_NPZ.relative_to(ROOT)),
            "bytes": OPERATOR_NPZ.stat().st_size,
            "sha256": sha256_file(OPERATOR_NPZ),
        },
        "arrays": {name: sha256_array(value) for name, value in arrays.items()},
        "caps_degrees": list(tangent.CAP_DEGREES),
        "effective_targets": {
            cap_key(cap): int(np.asarray(arrays[f"effective_{cap_key(cap)}"]).sum())
            for cap in tangent.CAP_DEGREES
        },
        "maximum_realized_degrees": {
            cap_key(cap): float(
                np.asarray(arrays[f"realized_degrees_{cap_key(cap)}"]).max(initial=0.0)
            )
            for cap in tangent.CAP_DEGREES
        },
        "bindings": {
            "authorization_sha256": sha256_file(AUTHORIZATION),
            "metadata_npz_sha256": sha256_file(METADATA_NPZ),
            "metadata_json_sha256": sha256_file(METADATA_JSON),
            "protocol_sha256": sha256_file(PROTOCOL),
            "builder_sha256": sha256_file(Path(__file__)),
            "tangent_sha256": sha256_file(TANGENT_PATH),
            "experiment_107_builder_sha256": sha256_file(BUILDER_PATH),
            "experiment_107_strict_sha256": sha256_file(STRICT_PATH),
            "source_model_sha256": sha256_file(
                SOURCE_DIRECTORY / "source-model.npz"
            ),
            "source_manifest_sha256": sha256_file(
                SOURCE_DIRECTORY / "manifest.json"
            ),
            "operator_json_sha256": sha256_file(OPERATOR_JSON),
            "production_implementation_sha256": {
                str(path.relative_to(ROOT)): sha256_file(path)
                for path in builder.PRODUCTION_IMPLEMENTATION_PATHS
            },
        },
        "transitive_inputs": candidate_transitive_inputs(builder),
        "nadig_truth_opened": False,
    }
    atomic_write_json(CANDIDATE_JSON, report)
    print(json.dumps(report["artifact"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

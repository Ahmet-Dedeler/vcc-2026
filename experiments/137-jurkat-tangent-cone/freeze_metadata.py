#!/usr/bin/env python3
"""Freeze Experiment 137's target and geometry metadata without effect access."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXP107 = ROOT / "experiments/107-jurkat-high-moi"
SOURCE_MODEL = EXP107 / "work/frozen-source-model/source-model.npz"
SOURCE_MANIFEST = EXP107 / "work/frozen-source-model/manifest.json"
STRICT_PATH = EXP107 / "strict_nadig_validation.py"
CHALLENGE = ROOT / "data/challenge"
NADIG = ROOT / "data/donors/nadig_jurkat"
DONOR_ROOT = ROOT / "data/donors"
ESSENTIAL_TARGETS = DONOR_ROOT / "k562_essential/perts.txt"
OUTPUT_NPZ = HERE / "work/metadata-freeze-v1.npz"
OUTPUT_JSON = HERE / "work/metadata-freeze-v1.json"

DONOR_ORDER = ("k562_gwps", "orion_hct116", "orion_hek293t")
EXPECTED_SOURCE_MODEL_SHA256 = (
    "03759a7b2cafa86d2643a58669235aca7fc38a33599aa403c67439f8979a7e26"
)
EXPECTED_COMMON_NON_VCC = 2_194
EXPECTED_CONSUMED = 1_800
EXPECTED_FRESH = 394
EXPECTED_TEST = 79
EXPECTED_DEVELOPMENT = 315
EXPECTED_DEVELOPMENT_FOLDS = 5
EXPECTED_FOLD_SIZE = 63
EXPECTED_FRESH_SHA256 = (
    "666834b07d902ae0fad4b703eda3b16cd5fb58966cdaeedacc0be5f7d41c3442"
)
EXPECTED_TEST_SHA256 = (
    "2ca0d292b34938cdfa415af70bcd546b535c84a99d568668ba584c2bc092c7ab"
)
EXPECTED_DEVELOPMENT_SHA256 = (
    "445d92f6213a6806dd42e3216ea94149f303d401f9d718755dd7dee54c812c57"
)
EXPECTED_FOLD_SHA256 = (
    "5952059adf9081b73c662efd349b76a7bc046ec9cced0524a3e4e1e036bf038d",
    "0d40c57b1f964f089deba38f38fa5cad3714f913954321e81318110faa7f37b6",
    "77ce43d3a6b7d78d443ca43ecf43551e127d08a646bbd60e75b1ec11f42aa495",
    "ac3a3ec5142d815770293c235991166b94767b30ead5a66cf9487fab4874c7ab",
    "0c4aaba96690df95acbca4ea193c57c1f44cebd51f87db3b034ccc6f08845aea",
)
EXPECTED_SCORE_GENE_COUNT = 223
EXPECTED_SCORE_GENE_SHA256 = (
    "44d284f4aa0bc134b0a366d2b972b2db6d8cc3841fd4a45f89d8bafb882cccba"
)
EXPECTED_OPERATOR_COLUMN_COUNT = 17_585
EXPECTED_OPERATOR_COLUMN_SHA256 = (
    "6f8783a19bcdfc7e21e56214e0eb8b15725d35a0e0f832d92a9783d5ab97522d"
)
EXPECTED_BASIS_ROWS = {
    "k562_gwps": (
        7_721,
        "df506057644ce0feb7f4bc4c9c4fe4dd10e205ba980f18f601989166b4219df8",
    ),
    "orion_hct116": (
        10_771,
        "cb030eb8574ce8312e9ce6408b9a8686d7b3bfc9bd7b2879c4c95aa8cd34f76d",
    ),
    "orion_hek293t": (
        10_764,
        "b7e7fc76d9860cf2c48e6bb542fc69906bf2098d5c3462b5097ae3f7d181f5f3",
    ),
}


class FreezeError(RuntimeError):
    """Raised when a metadata-only experiment boundary changes."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise FreezeError(message)


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_bytes), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_strings(values: Sequence[str]) -> str:
    return hashlib.sha256(
        "".join(f"{value}\n" for value in values).encode("utf-8")
    ).hexdigest()


def sha256_int64(values: np.ndarray) -> str:
    array = np.ascontiguousarray(np.asarray(values, dtype=np.int64))
    return hashlib.sha256(array.tobytes(order="C")).hexdigest()


def read_lines(path: Path, *, require_unique: bool = True) -> tuple[str, ...]:
    values = tuple(path.read_text(encoding="utf-8").splitlines())
    require(bool(values) and all(values), f"{path}: empty axis entry")
    if require_unique:
        require(len(values) == len(set(values)), f"{path}: duplicate axis entry")
    return values


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == [column], f"{path}: header changed")
        values = tuple(str(row[column]) for row in reader)
    require(bool(values) and all(values), f"{path}: empty {column}")
    require(len(values) == len(set(values)), f"{path}: duplicate {column}")
    return values


def load_exp107_metadata() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    require(
        sha256_file(SOURCE_MODEL) == EXPECTED_SOURCE_MODEL_SHA256,
        "Experiment 107 frozen source model changed",
    )
    # Deliberately access only axes and support masks. exported_effects and all
    # coefficient arrays remain unopened in this metadata stage.
    with np.load(SOURCE_MODEL, allow_pickle=False) as archive:
        targets = np.asarray(archive["target_symbols"]).astype(str)
        transferable = np.asarray(archive["transferable"], dtype=bool)
        responses = np.asarray(archive["response_symbols"]).astype(str)
        exportable = np.asarray(archive["response_exportable"], dtype=bool)
    require(targets.shape == transferable.shape, "source target metadata misaligns")
    require(responses.shape == exportable.shape, "source response metadata misaligns")
    return targets, transferable, responses, exportable


def load_strict_module() -> Any:
    spec = importlib.util.spec_from_file_location("exp137_strict_metadata", STRICT_PATH)
    require(spec is not None and spec.loader is not None, "cannot load strict metadata code")
    module = importlib.util.module_from_spec(spec)
    sys.modules[str(spec.name)] = module
    spec.loader.exec_module(module)
    return module


def inherited_order(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(
        sorted(
            values,
            key=lambda symbol: (
                hashlib.sha256(f"exp107|{symbol}".encode("utf-8")).digest(),
                symbol,
            ),
        )
    )


def atomic_savez(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".npz"
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        np.savez_compressed(temporary, **arrays)
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


def build_freeze() -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    strict = load_strict_module()
    source_targets, transferable, responses, exportable = load_exp107_metadata()
    nadig_targets = read_lines(NADIG / "perts.txt")
    nadig_genes = read_lines(NADIG / "genes.txt")
    vcc_targets = read_csv_column(CHALLENGE / "pert_counts.csv", "target_gene")
    official_genes = read_csv_column(CHALLENGE / "gene_names.csv", "gene_name")
    donor_genes = {
        name: read_lines(DONOR_ROOT / name / "genes.txt", require_unique=False)
        for name in DONOR_ORDER
    }
    donor_targets = {
        name: read_lines(DONOR_ROOT / name / "perts.txt", require_unique=False)
        for name in DONOR_ORDER
    }
    essential = set(read_lines(ESSENTIAL_TARGETS))

    eligible = set(source_targets) & set(nadig_targets)
    eligible.difference_update(vcc_targets)
    ordered = inherited_order(eligible)
    require(len(ordered) == EXPECTED_COMMON_NON_VCC, "eligible target count changed")
    consumed = ordered[:EXPECTED_CONSUMED]
    fresh = ordered[EXPECTED_CONSUMED:]
    test = fresh[:EXPECTED_TEST]
    development = fresh[EXPECTED_TEST:]
    folds = tuple(
        development[index : index + EXPECTED_FOLD_SIZE]
        for index in range(0, len(development), EXPECTED_FOLD_SIZE)
    )
    require(len(fresh) == EXPECTED_FRESH, "fresh target count changed")
    require(len(test) == EXPECTED_TEST, "test target count changed")
    require(len(development) == EXPECTED_DEVELOPMENT, "development count changed")
    require(
        len(folds) == EXPECTED_DEVELOPMENT_FOLDS
        and all(len(fold) == EXPECTED_FOLD_SIZE for fold in folds),
        "development folds changed",
    )
    require(sha256_strings(fresh) == EXPECTED_FRESH_SHA256, "fresh order changed")
    require(sha256_strings(test) == EXPECTED_TEST_SHA256, "test order changed")
    require(
        sha256_strings(development) == EXPECTED_DEVELOPMENT_SHA256,
        "development order changed",
    )
    require(
        tuple(sha256_strings(fold) for fold in folds) == EXPECTED_FOLD_SHA256,
        "development fold order changed",
    )

    source_index = {symbol: index for index, symbol in enumerate(source_targets)}
    nadig_index = {symbol: index for index, symbol in enumerate(nadig_targets)}
    fresh_source_indices = np.asarray(
        [source_index[symbol] for symbol in fresh], dtype=np.int64
    )
    fresh_nadig_indices = np.asarray(
        [nadig_index[symbol] for symbol in fresh], dtype=np.int64
    )
    fresh_transferable = transferable[fresh_source_indices]
    require(
        int(fresh_transferable[:EXPECTED_TEST].sum()) == EXPECTED_TEST,
        "one or more frozen test targets lost support",
    )
    unsupported = tuple(
        symbol
        for symbol, supported in zip(fresh, fresh_transferable, strict=True)
        if not supported
    )
    require(
        unsupported == ("ZNF658", "LPA", "NOMO3"),
        f"unsupported fresh targets changed: {unsupported}",
    )

    score_genes = strict.score_axis(
        responses.tolist(),
        exportable,
        nadig_genes,
        vcc_targets,
        fresh,
    )
    require(len(score_genes) == EXPECTED_SCORE_GENE_COUNT, "score axis count changed")
    require(
        sha256_strings(score_genes) == EXPECTED_SCORE_GENE_SHA256,
        "score axis changed",
    )
    response_index = {symbol: index for index, symbol in enumerate(responses)}
    nadig_gene_index = {symbol: index for index, symbol in enumerate(nadig_genes)}
    official_index = {symbol: index for index, symbol in enumerate(official_genes)}
    operator_columns = strict.whitener_columns(
        official_genes,
        donor_genes,
        vcc_targets,
        fresh,
    )
    require(
        len(operator_columns) == EXPECTED_OPERATOR_COLUMN_COUNT,
        "operator column count changed",
    )
    require(
        sha256_int64(operator_columns) == EXPECTED_OPERATOR_COLUMN_SHA256,
        "operator columns changed",
    )
    basis_rows = strict.select_whitener_basis_rows(
        donor_targets,
        essential,
        set(fresh),
    )
    for donor in DONOR_ORDER:
        expected_count, expected_hash = EXPECTED_BASIS_ROWS[donor]
        require(len(basis_rows[donor]) == expected_count, f"{donor} basis count changed")
        require(
            sha256_int64(basis_rows[donor]) == expected_hash,
            f"{donor} basis rows changed",
        )

    fold_id = np.repeat(
        np.arange(1, EXPECTED_DEVELOPMENT_FOLDS + 1, dtype=np.int8),
        EXPECTED_FOLD_SIZE,
    )
    arrays = {
        "eligible_targets": np.asarray(ordered, dtype="U"),
        "consumed_targets": np.asarray(consumed, dtype="U"),
        "fresh_targets": np.asarray(fresh, dtype="U"),
        "test_targets": np.asarray(test, dtype="U"),
        "development_targets": np.asarray(development, dtype="U"),
        "development_fold_id": fold_id,
        "fresh_source_indices": fresh_source_indices,
        "fresh_nadig_indices": fresh_nadig_indices,
        "fresh_transferable": fresh_transferable,
        "score_genes": np.asarray(score_genes, dtype="U"),
        "score_indices_in_source": np.asarray(
            [response_index[symbol] for symbol in score_genes], dtype=np.int64
        ),
        "score_indices_in_nadig": np.asarray(
            [nadig_gene_index[symbol] for symbol in score_genes], dtype=np.int64
        ),
        "score_indices_in_official": np.asarray(
            [official_index[symbol] for symbol in score_genes], dtype=np.int64
        ),
        "operator_columns": np.asarray(operator_columns, dtype=np.int64),
        **{
            f"basis_rows_{donor}": np.asarray(basis_rows[donor], dtype=np.int64)
            for donor in DONOR_ORDER
        },
    }
    report: dict[str, Any] = {
        "schema_version": 1,
        "experiment": "137-jurkat-tangent-cone",
        "status": "metadata_frozen_before_fresh_target_outcome_access",
        "outcome_boundary": {
            "file": "data/donors/nadig_jurkat/delta.npy",
            "file_was_loaded_by_experiment_107": True,
            "fresh_rows_used_by_experiment_107_selection_or_reporting": False,
            "fresh_rows_parsed_by_experiment_137": False,
        },
        "split": {
            "eligible_common_non_vcc": len(ordered),
            "experiment_107_consumed": len(consumed),
            "fresh": len(fresh),
            "test": len(test),
            "development": len(development),
            "development_folds": len(folds),
            "development_fold_size": EXPECTED_FOLD_SIZE,
            "order": "SHA256 UTF-8 bytes of 'exp107|' plus target symbol, then symbol",
            "test_rule": "first 79 targets after the consumed first 1800",
            "development_rule": "remaining 315 targets in five contiguous 63-target folds",
            "fresh_sha256_newline_terminated": sha256_strings(fresh),
            "test_sha256_newline_terminated": sha256_strings(test),
            "development_sha256_newline_terminated": sha256_strings(development),
            "development_fold_sha256_newline_terminated": [
                sha256_strings(fold) for fold in folds
            ],
            "test_targets": list(test),
            "unsupported_fresh_targets": list(unsupported),
            "test_transferable_targets": int(fresh_transferable[:EXPECTED_TEST].sum()),
            "development_transferable_targets": int(
                fresh_transferable[EXPECTED_TEST:].sum()
            ),
        },
        "geometry": {
            "score_gene_count": len(score_genes),
            "score_gene_sha256_newline_terminated": sha256_strings(score_genes),
            "operator_column_count": len(operator_columns),
            "operator_columns_int64_sha256": sha256_int64(operator_columns),
            "basis_rows": {
                donor: {
                    "count": len(basis_rows[donor]),
                    "int64_sha256": sha256_int64(basis_rows[donor]),
                }
                for donor in DONOR_ORDER
            },
        },
        "bindings": {
            str(path.relative_to(ROOT)): {
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in (
                SOURCE_MODEL,
                SOURCE_MANIFEST,
                EXP107 / "work/metadata-crosswalk.npz",
                EXP107 / "work/target-support.npz",
                EXP107 / "strict_validation_contract.json",
                STRICT_PATH,
                NADIG / "perts.txt",
                NADIG / "genes.txt",
                CHALLENGE / "pert_counts.csv",
                CHALLENGE / "gene_names.csv",
                ESSENTIAL_TARGETS,
                *(DONOR_ROOT / donor / "perts.txt" for donor in DONOR_ORDER),
                *(DONOR_ROOT / donor / "genes.txt" for donor in DONOR_ORDER),
            )
        },
    }
    return arrays, report


def main() -> int:
    arrays, report = build_freeze()
    atomic_savez(OUTPUT_NPZ, **arrays)
    report["artifact"] = {
        "path": str(OUTPUT_NPZ.relative_to(ROOT)),
        "bytes": OUTPUT_NPZ.stat().st_size,
        "sha256": sha256_file(OUTPUT_NPZ),
    }
    report["implementation"] = {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "bytes": Path(__file__).stat().st_size,
        "sha256": sha256_file(Path(__file__)),
    }
    atomic_write_json(OUTPUT_JSON, report)
    print(
        json.dumps(
            {
                "status": report["status"],
                "fresh": report["split"]["fresh"],
                "test": report["split"]["test"],
                "development": report["split"]["development"],
                "artifact_sha256": report["artifact"]["sha256"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

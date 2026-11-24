#!/usr/bin/env python3
"""Streaming cell-call preparation for the GSE249595 Jurkat source.

The archive is intentionally processed one channel at a time.  This stage
does not fit perturbation effects or inspect any held perturbation truth.
"""

from __future__ import annotations

import argparse
import fcntl
import gzip
import hashlib
import json
import os
import sys
import tarfile
import tempfile
import time
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any, BinaryIO

import numpy as np
import scipy

from preprocess import (
    ambient_feature_counts,
    barcode_ranks_thresholds,
    emptydrops_fdr,
    emptydrops_monte_carlo,
    observed_data_log_probabilities,
    safe_good_turing_proportions,
    scan_transcriptome_summary,
)


HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "downloads" / "GSE249595_RAW.tar"
SCHEMA_REPORT = HERE / "audit" / "schema.json"
GZIP_REPORT = HERE / "audit" / "gzip-integrity.json"
SUPPORT_REPORT = HERE / "audit" / "support" / "COMPLETE.json"
ORACLE_RECEIPT = HERE / "R_ORACLE_RECEIPT.json"
PREFLIGHT_RECEIPT = HERE / "SOURCE_PREFLIGHT_RECEIPT.json"
PREFLIGHT_CODE = HERE / "source_preflight.py"
PROTOCOL = HERE / "preregistered_protocol.json"
PREFLIGHT_TEST_MODULES = (
    "test_preprocess.py",
    "test_source_stage.py",
    "test_guide_calling.py",
    "test_demultiplex.py",
)
PREFLIGHT_IMPLEMENTATIONS = (
    "source_stage.py",
    "preprocess.py",
    "guide_calling.py",
    "demultiplex.py",
)
WORK = HERE / "work"
LOWER = 100
FDR_THRESHOLD = 0.01
EXPECTED_RAW_SHA256 = "3c0c5fa8cab2983dfe0fa3b38b31662bb68e191e6dbaf3d03af02a67587a67a9"


class SourceStageError(RuntimeError):
    """A frozen source-stage invariant failed."""


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_bytes), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
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
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def atomic_savez(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            np.savez_compressed(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


@contextmanager
def exclusive_lock() -> Iterator[None]:
    WORK.mkdir(parents=True, exist_ok=True)
    path = WORK / ".source-stage.lock"
    with path.open("a+b") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise SourceStageError(
                "another experiment-107 source stage is active"
            ) from error
        yield


def load_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise SourceStageError(f"cannot read {path}: {error}") from error
    if not isinstance(payload, dict):
        raise SourceStageError(f"{path} is not a JSON object")
    return payload


def validate_preflight_receipt() -> dict[str, Any]:
    """Require a current passing synthetic receipt before source execution."""

    receipt = load_json(PREFLIGHT_RECEIPT)
    protocol = load_json(PROTOCOL)
    try:
        frozen_numpy = protocol["cell_and_assignment_freeze"][
            "transcriptome_cell_call"
        ]["rng_contract"]["numpy_version"]
        frozen_scipy = protocol["cell_and_assignment_freeze"]["guide_assignment"][
            "scipy_version"
        ]
    except (KeyError, TypeError) as error:
        raise SourceStageError(
            "protocol is missing frozen NumPy/SciPy versions"
        ) from error
    expected_scalars = {
        "status": "pass",
        "numpy_version": frozen_numpy,
        "scipy_version": frozen_scipy,
        "protocol_sha256": sha256_file(PROTOCOL),
        "preflight_implementation_sha256": sha256_file(PREFLIGHT_CODE),
        "held_outcomes_opened": False,
    }
    for field, expected in expected_scalars.items():
        if receipt.get(field) != expected:
            raise SourceStageError(f"source preflight receipt is stale at {field}")
    if np.__version__ != frozen_numpy or scipy.__version__ != frozen_scipy:
        raise SourceStageError("current NumPy/SciPy runtime differs from the freeze")
    implementations = receipt.get("implementation_sha256")
    expected_implementations = {
        filename: sha256_file(HERE / filename)
        for filename in PREFLIGHT_IMPLEMENTATIONS
    }
    if implementations != expected_implementations:
        raise SourceStageError("source preflight implementation hashes changed")
    tests = receipt.get("test_module_sha256")
    expected_tests = {
        filename: sha256_file(HERE / filename)
        for filename in PREFLIGHT_TEST_MODULES
    }
    if tests != expected_tests:
        raise SourceStageError("source preflight test hashes changed")
    if (
        not isinstance(receipt.get("tests_run"), int)
        or receipt["tests_run"] <= 0
        or receipt.get("failures") != 0
        or receipt.get("errors") != 0
        or receipt.get("skipped") != 0
    ):
        raise SourceStageError("source preflight synthetic tests did not fully pass")
    transcript_sha256 = receipt.get("transcript_sha256")
    if (
        not isinstance(transcript_sha256, str)
        or len(transcript_sha256) != 64
        or any(character not in "0123456789abcdef" for character in transcript_sha256)
    ):
        raise SourceStageError("source preflight transcript hash is invalid")
    return receipt


def validate_prerequisites(*, verify_raw_hash: bool) -> dict[str, Any]:
    for path in (
        ARCHIVE,
        SCHEMA_REPORT,
        GZIP_REPORT,
        SUPPORT_REPORT,
        ORACLE_RECEIPT,
        PREFLIGHT_RECEIPT,
    ):
        if not path.is_file() or path.is_symlink():
            raise SourceStageError(f"missing or unsafe prerequisite: {path}")
    validate_preflight_receipt()

    schema = load_json(SCHEMA_REPORT)
    gzip_report = load_json(GZIP_REPORT)
    support = load_json(SUPPORT_REPORT)
    oracle = load_json(ORACLE_RECEIPT)
    if schema.get("status") != "pass":
        raise SourceStageError("schema audit did not pass")
    if not schema.get("barcode_axes_verified"):
        raise SourceStageError("barcode axes were not verified")
    if not schema.get("assignment_coordinates_scanned"):
        raise SourceStageError("assignment coordinates were not scanned")
    if schema.get("transcriptome_coordinates_scanned"):
        raise SourceStageError("schema audit unexpectedly scanned transcriptome values")
    if gzip_report.get("status") != "pass" or not gzip_report.get(
        "all_members_read_to_gzip_eof"
    ):
        raise SourceStageError("full gzip integrity audit did not pass")
    if support.get("status") != "complete":
        raise SourceStageError("assignment support audit is incomplete")
    if oracle.get("status") != "pass":
        raise SourceStageError("the pinned R oracle did not pass")

    bindings = [
        schema.get("input_binding", {}),
        gzip_report.get("input_binding", {}),
        support.get("input_binding", {}),
    ]
    raw_hashes = {binding.get("raw_archive_sha256") for binding in bindings}
    if raw_hashes != {EXPECTED_RAW_SHA256}:
        raise SourceStageError(
            f"audit reports disagree on raw archive hash: {raw_hashes}"
        )
    if ARCHIVE.stat().st_size != 3_753_820_160:
        raise SourceStageError("raw archive byte count changed")
    if verify_raw_hash and sha256_file(ARCHIVE) != EXPECTED_RAW_SHA256:
        raise SourceStageError("raw archive SHA-256 changed")
    return schema


def channel_record(schema: dict[str, Any], channel: int) -> dict[str, Any]:
    matches = [
        record
        for record in schema.get("channels", [])
        if int(record.get("channel", -1)) == channel
    ]
    if len(matches) != 1:
        raise SourceStageError(f"expected one schema record for channel {channel}")
    return matches[0]


def member_name(record: dict[str, Any], modality: str, role: str) -> str:
    return (
        f"{record['accession']}_channel{record['channel']}_"
        f"{modality}_{role}"
        f"{'.mtx.gz' if role == 'matrix' else '.tsv.gz'}"
    )


@contextmanager
def gzip_member(
    archive: tarfile.TarFile,
    name: str,
) -> Iterator[BinaryIO]:
    extracted = archive.extractfile(name)
    if extracted is None:
        raise SourceStageError(f"archive member is missing: {name}")
    try:
        with gzip.GzipFile(fileobj=extracted, mode="rb") as decompressed:
            yield decompressed
    finally:
        extracted.close()


def read_feature_rows(
    archive: tarfile.TarFile,
    name: str,
) -> list[tuple[str, ...]]:
    rows: list[tuple[str, ...]] = []
    with gzip_member(archive, name) as stream:
        for raw in stream:
            fields = tuple(raw.decode("utf-8").rstrip("\r\n").split("\t"))
            if len(fields) < 2:
                raise SourceStageError(f"{name}: malformed feature row")
            rows.append(fields)
    return rows


def channel_directory(channel: int) -> Path:
    return WORK / f"channel{channel:02d}"


def prepare_channel(
    schema: dict[str, Any],
    channel: int,
) -> dict[str, Any]:
    record = channel_record(schema, channel)
    output_directory = channel_directory(channel)
    arrays_path = output_directory / "cell-call-inputs.npz"
    report_path = output_directory / "prepare-report.json"
    matrix_name = member_name(record, "transcriptome", "matrix")
    features_name = member_name(record, "transcriptome", "features")
    preflight = validate_preflight_receipt()

    started = time.time()
    with tarfile.open(ARCHIVE, mode="r:") as archive:
        feature_rows = read_feature_rows(archive, features_name)
        expected_rows = int(record["modalities"]["transcriptome"]["matrix"]["rows"])
        if len(feature_rows) != expected_rows:
            raise SourceStageError("transcriptome feature count changed")
        symbols = [row[1] for row in feature_rows]
        mitochondrial = np.flatnonzero(
            np.fromiter(
                (symbol.upper().startswith("MT-") for symbol in symbols),
                dtype=bool,
                count=len(symbols),
            )
        )

        with gzip_member(archive, matrix_name) as matrix:
            summary = scan_transcriptome_summary(
                matrix,
                mitochondrial_features=mitochondrial,
            )
        if (
            summary.rows != expected_rows
            or summary.columns
            != int(record["modalities"]["transcriptome"]["matrix"]["columns"])
            or summary.observed_nnz
            != int(record["modalities"]["transcriptome"]["matrix"]["nnz"])
        ):
            raise SourceStageError("transcriptome scan disagrees with frozen schema")

        rank = barcode_ranks_thresholds(summary.totals, lower=LOWER)
        with gzip_member(archive, matrix_name) as matrix:
            ambient_counts = ambient_feature_counts(
                matrix,
                totals=summary.totals,
                lower=LOWER,
            )
        feature_nonzero = summary.feature_totals > 0
        if not np.any(feature_nonzero):
            raise SourceStageError("transcriptome has no expressed features")
        ambient_probability = np.zeros(summary.rows, dtype=np.float64)
        ambient_probability[feature_nonzero] = safe_good_turing_proportions(
            ambient_counts[feature_nonzero]
        )

        tested = summary.totals > LOWER
        with gzip_member(archive, matrix_name) as matrix:
            observed = observed_data_log_probabilities(
                matrix,
                totals=summary.totals,
                ambient_probability=ambient_probability,
                tested=tested,
            )

    atomic_savez(
        arrays_path,
        totals=summary.totals,
        detected_features=summary.detected_features,
        mitochondrial_counts=summary.mitochondrial_counts,
        feature_totals=summary.feature_totals,
        ambient_counts=ambient_counts,
        ambient_probability=ambient_probability,
        tested_indices=np.flatnonzero(tested).astype(np.int32),
        observed_data_log_probability=observed,
    )
    report = {
        "status": "complete",
        "channel": channel,
        "accession": record["accession"],
        "matrix_member": matrix_name,
        "raw_archive_sha256": EXPECTED_RAW_SHA256,
        "source_stage_sha256": sha256_file(Path(__file__)),
        "preprocess_sha256": sha256_file(HERE / "preprocess.py"),
        "r_oracle_receipt_sha256": sha256_file(ORACLE_RECEIPT),
        "preflight_receipt_sha256": sha256_file(PREFLIGHT_RECEIPT),
        "preflight_implementation_sha256": preflight[
            "preflight_implementation_sha256"
        ],
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "arrays_path": str(arrays_path.relative_to(HERE)),
        "arrays_sha256": sha256_file(arrays_path),
        "shape": [summary.rows, summary.columns],
        "nnz": summary.observed_nnz,
        "coordinate_order": summary.coordinate_order,
        "nonzero_barcodes": int(np.count_nonzero(summary.totals)),
        "tested_barcodes": int(np.count_nonzero(tested)),
        "maximum_total": int(summary.totals.max()),
        "knee": rank.knee,
        "inflection": rank.inflection,
        "expressed_features": int(np.count_nonzero(feature_nonzero)),
        "mitochondrial_features": int(mitochondrial.size),
        "ambient_barcodes": int(np.count_nonzero(summary.totals <= LOWER)),
        "elapsed_seconds": time.time() - started,
    }
    atomic_write_json(report_path, report)
    return report


def call_channel(
    schema: dict[str, Any],
    channel: int,
    *,
    iterations: int,
    block_size: int,
) -> dict[str, Any]:
    channel_record(schema, channel)
    directory = channel_directory(channel)
    inputs_path = directory / "cell-call-inputs.npz"
    prepare_report_path = directory / "prepare-report.json"
    if not inputs_path.is_file() or not prepare_report_path.is_file():
        raise SourceStageError(f"channel {channel} has not completed prepare")
    prepare_report = load_json(prepare_report_path)
    preflight = validate_preflight_receipt()
    if prepare_report.get("arrays_sha256") != sha256_file(inputs_path):
        raise SourceStageError(f"channel {channel} prepared arrays changed")
    if prepare_report.get("source_stage_sha256") != sha256_file(Path(__file__)):
        raise SourceStageError(f"channel {channel} prepare used different stage code")
    if prepare_report.get("preprocess_sha256") != sha256_file(HERE / "preprocess.py"):
        raise SourceStageError(f"channel {channel} prepare used a different port")
    if prepare_report.get("r_oracle_receipt_sha256") != sha256_file(ORACLE_RECEIPT):
        raise SourceStageError(f"channel {channel} prepare used a different oracle")
    expected_preflight = {
        "preflight_receipt_sha256": sha256_file(PREFLIGHT_RECEIPT),
        "preflight_implementation_sha256": preflight[
            "preflight_implementation_sha256"
        ],
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
    }
    for field, expected in expected_preflight.items():
        if prepare_report.get(field) != expected:
            raise SourceStageError(
                f"channel {channel} prepare used a different {field}"
            )

    with np.load(inputs_path, allow_pickle=False) as inputs:
        tested_indices = inputs["tested_indices"].astype(np.int64, copy=False)
        tested_totals = inputs["totals"][tested_indices]
        observed = inputs["observed_data_log_probability"]
        ambient_probability = inputs["ambient_probability"]
        all_totals = inputs["totals"]
        if observed.shape != tested_indices.shape:
            raise SourceStageError("tested index and probability arrays disagree")

        started = time.time()
        monte_carlo = emptydrops_monte_carlo(
            tested_totals,
            observed,
            ambient_probability,
            channel_number=channel,
            iterations=iterations,
            iteration_block_size=block_size,
        )
        all_p = np.full(all_totals.shape, np.nan, dtype=np.float64)
        all_p[tested_indices] = monte_carlo.p_value
        rank = barcode_ranks_thresholds(all_totals, lower=LOWER)
        all_fdr = emptydrops_fdr(
            all_totals,
            all_p,
            lower=LOWER,
            retain=rank.knee,
        )
        accepted = np.flatnonzero(all_fdr <= FDR_THRESHOLD).astype(np.int32)
        tested_fdr = all_fdr[tested_indices]

    label = "production" if iterations == 10_000 else f"pilot-{iterations}"
    output_path = directory / f"cell-calls-{label}.npz"
    report_path = directory / f"cell-calls-{label}.json"
    atomic_savez(
        output_path,
        tested_indices=tested_indices.astype(np.int32),
        p_value=monte_carlo.p_value,
        fdr=tested_fdr,
        limited=monte_carlo.limited,
        accepted_indices=accepted,
    )
    report = {
        "status": "complete",
        "stage": label,
        "channel": channel,
        "iterations": iterations,
        "iteration_block_size": block_size,
        "prepare_report_sha256": sha256_file(prepare_report_path),
        "prepared_arrays_sha256": sha256_file(inputs_path),
        "source_stage_sha256": sha256_file(Path(__file__)),
        "preprocess_sha256": sha256_file(HERE / "preprocess.py"),
        "r_oracle_receipt_sha256": sha256_file(ORACLE_RECEIPT),
        **expected_preflight,
        "output_path": str(output_path.relative_to(HERE)),
        "output_sha256": sha256_file(output_path),
        "tested_barcodes": int(tested_indices.size),
        "accepted_barcodes": int(accepted.size),
        "forced_retained_barcodes": int(np.count_nonzero(tested_totals >= rank.knee)),
        "limited_barcodes": int(np.count_nonzero(monte_carlo.limited)),
        "knee": rank.knee,
        "elapsed_seconds": time.time() - started,
    }
    atomic_write_json(report_path, report)
    return report


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument("command", choices=("prepare", "call"))
    result.add_argument("--channel", type=int, required=True, choices=range(1, 17))
    result.add_argument("--iterations", type=int, default=10_000)
    result.add_argument("--block-size", type=int, default=64)
    result.add_argument("--verify-raw-hash", action="store_true")
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        with exclusive_lock():
            schema = validate_prerequisites(verify_raw_hash=args.verify_raw_hash)
            if args.command == "prepare":
                report = prepare_channel(schema, args.channel)
            else:
                report = call_channel(
                    schema,
                    args.channel,
                    iterations=args.iterations,
                    block_size=args.block_size,
                )
    except SourceStageError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

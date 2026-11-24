#!/usr/bin/env python3
"""Validate the experiment-061 prerequisite before experiment 069 dispatch."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()


def read_checksum(path: Path) -> tuple[str, str]:
    fields = path.read_text().strip().split()
    if len(fields) != 2 or not SHA256_RE.fullmatch(fields[0]):
        raise ValueError(f"{path}: invalid sha256sum file")
    return fields[0], fields[1].lstrip("*")


def validate_gate(
    *,
    results_path: Path,
    baseline_agg_path: Path,
    candidate_agg_path: Path,
    build_report_path: Path,
    candidate_checksum_path: Path,
    exit_code_path: Path,
    candidate_object_bytes: int | None = None,
) -> dict[str, object]:
    exit_code = exit_code_path.read_text().strip()
    if exit_code != "0":
        raise ValueError(f"experiment 061 exit code is {exit_code!r}, not '0'")

    results = json.loads(results_path.read_text())
    if results.get("experiment") != "061-perturbmatch-hct-count":
        raise ValueError("unexpected experiment-061 results identity")
    if results.get("status") != "scored":
        raise ValueError("experiment 061 is not marked scored")
    avg6_delta = float(results["avg6"]["candidate_minus_baseline"])
    if not math.isfinite(avg6_delta) or avg6_delta < 0:
        raise ValueError(
            f"experiment 061 count-level avg6 delta is {avg6_delta}, must be >= 0"
        )

    score_files = results.get("score_files", {})
    expected_baseline = str(score_files["baseline"]["sha256"])
    expected_candidate = str(score_files["candidate"]["sha256"])
    if not SHA256_RE.fullmatch(expected_baseline):
        raise ValueError("experiment 061 baseline aggregate hash is invalid")
    if not SHA256_RE.fullmatch(expected_candidate):
        raise ValueError("experiment 061 candidate aggregate hash is invalid")
    actual_baseline = file_sha256(baseline_agg_path)
    actual_candidate = file_sha256(candidate_agg_path)
    if actual_baseline != expected_baseline:
        raise ValueError("experiment 061 baseline aggregate hash mismatch")
    if actual_candidate != expected_candidate:
        raise ValueError("experiment 061 candidate aggregate hash mismatch")

    build = json.loads(build_report_path.read_text())
    if build.get("experiment") != "061-perturbmatch-hct-count":
        raise ValueError("unexpected experiment-061 build-report identity")
    if build.get("status") != "candidate_built":
        raise ValueError("experiment 061 candidate was not marked built")
    report_hash = str(build["candidate_h5ad"]["sha256"])
    report_bytes = int(build["candidate_h5ad"]["bytes"])
    checksum_hash, checksum_name = read_checksum(candidate_checksum_path)
    if report_hash != checksum_hash:
        raise ValueError("candidate H5AD checksum disagrees with build report")
    if Path(checksum_name).name != "candidate.h5ad":
        raise ValueError("candidate checksum does not name candidate.h5ad")
    if report_bytes <= 0:
        raise ValueError("candidate H5AD byte count is not positive")
    if candidate_object_bytes is not None and candidate_object_bytes != report_bytes:
        raise ValueError(
            "candidate H5AD S3 size disagrees with experiment-061 build report"
        )

    return {
        "status": "pass",
        "experiment": results["experiment"],
        "exit_code": 0,
        "count_level_avg6_delta": avg6_delta,
        "baseline_aggregate_sha256": actual_baseline,
        "candidate_aggregate_sha256": actual_candidate,
        "candidate_h5ad_sha256": report_hash,
        "candidate_h5ad_bytes": report_bytes,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--baseline-agg", required=True, type=Path)
    parser.add_argument("--candidate-agg", required=True, type=Path)
    parser.add_argument("--build-report", required=True, type=Path)
    parser.add_argument("--candidate-checksum", required=True, type=Path)
    parser.add_argument("--exit-code", required=True, type=Path)
    parser.add_argument("--candidate-object-bytes", type=int)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    result = validate_gate(
        results_path=args.results,
        baseline_agg_path=args.baseline_agg,
        candidate_agg_path=args.candidate_agg,
        build_report_path=args.build_report,
        candidate_checksum_path=args.candidate_checksum,
        exit_code_path=args.exit_code,
        candidate_object_bytes=args.candidate_object_bytes,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

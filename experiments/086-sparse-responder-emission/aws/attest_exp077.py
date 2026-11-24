#!/usr/bin/env python3
"""Emit a canonical, read-only terminal attestation for experiment 077.

This helper is intended to run through a dedicated AWS SSM Run Command after
the exp077 systemd unit has completed. Its one-line JSON stdout is retained by
SSM and therefore provides an AWS-authoritative binding between the launched
worker, the exact source/runner, and the terminal result hashes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path


RUN_ID_RE = re.compile(r"[0-9]{8}T[0-9]{6}Z(?:-[a-z0-9][a-z0-9-]{0,31})?")
PREFIX_RE = re.compile(r"[a-z0-9][a-z0-9._/-]{1,255}")
REQUIRED_FILES = {
    "source_archive": "source.tar.gz",
    "runner": "/root/exp077-run.sh",
    "grid": "repo/experiments/077-low-alpha-union/grid.json",
    "summarizer": "repo/experiments/077-low-alpha-union/summarize.py",
    "exit_code": "exit-code.txt",
    "summary": "output/summary.json",
    "results": "output/results.json",
    "artifact_manifest": "output/artifact-sha256.txt",
}


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def require_regular_file(path: Path) -> Path:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"required regular non-symlink file is absent: {path}")
    return path


def canonical_existing_directory(raw: str) -> Path:
    path = Path(raw)
    if not path.is_absolute():
        raise ValueError("work root must be absolute")
    resolved = path.resolve(strict=True)
    if str(resolved) != raw or path.is_symlink() or not path.is_dir():
        raise ValueError(f"work root is not a canonical real directory: {raw!r}")
    return path


def systemd_state(unit: str) -> dict[str, str]:
    if re.fullmatch(r"vcc-exp077-[0-9A-Za-z_-]+\.service", unit) is None:
        raise ValueError(f"unsafe systemd unit {unit!r}")
    process = subprocess.run(
        [
            "systemctl",
            "show",
            unit,
            "--property=ActiveState,SubState,Result,ExecMainStatus",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if process.returncode != 0:
        raise RuntimeError(f"systemctl show failed: {process.stderr.strip()}")
    record: dict[str, str] = {}
    for line in process.stdout.splitlines():
        key, separator, value = line.partition("=")
        if not separator or key in record:
            raise ValueError(f"malformed systemd property line {line!r}")
        record[key] = value
    expected_keys = {"ActiveState", "SubState", "Result", "ExecMainStatus"}
    if set(record) != expected_keys:
        raise ValueError(f"unexpected systemd properties: {sorted(record)}")
    if (
        record["ActiveState"] != "inactive"
        or record["Result"] != "success"
        or record["ExecMainStatus"] != "0"
    ):
        raise ValueError(
            f"experiment 077 unit is not successful and terminal: {record}"
        )
    return record


def build_attestation(args: argparse.Namespace) -> dict[str, object]:
    if RUN_ID_RE.fullmatch(args.run_id) is None:
        raise ValueError(f"unsafe run id {args.run_id!r}")
    if PREFIX_RE.fullmatch(args.source_prefix) is None:
        raise ValueError(f"unsafe source prefix {args.source_prefix!r}")
    if PREFIX_RE.fullmatch(args.output_prefix) is None:
        raise ValueError(f"unsafe output prefix {args.output_prefix!r}")
    expected_source = f"input/exp077-low-alpha-union/{args.run_id}"
    expected_output = f"output/exp077-low-alpha-union/{args.run_id}"
    if args.source_prefix != expected_source or args.output_prefix != expected_output:
        raise ValueError("experiment 077 prefixes are not canonical for the run id")

    work = canonical_existing_directory(args.work_root)
    expected_work = f"/work/exp077-{args.run_id}"
    if str(work) != expected_work:
        raise ValueError(f"work root {work} != {expected_work}")

    paths: dict[str, Path] = {}
    for label, relative in REQUIRED_FILES.items():
        path = Path(relative) if relative.startswith("/") else work / relative
        paths[label] = require_regular_file(path)
    if paths["exit_code"].read_text().strip() != "0":
        raise ValueError("experiment 077 exit code is not zero")

    unit = f"vcc-exp077-{args.run_id}.service"
    return {
        "schema_version": 1,
        "kind": "exp077-terminal-attestation",
        "run_id": args.run_id,
        "instance_id": args.instance_id,
        "region": args.region,
        "source_prefix": args.source_prefix,
        "output_prefix": args.output_prefix,
        "systemd_unit": unit,
        "systemd_state": systemd_state(unit),
        "sha256": {label: sha256_file(path) for label, path in sorted(paths.items())},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--instance-id", required=True)
    parser.add_argument("--region", required=True)
    parser.add_argument("--source-prefix", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--work-root", required=True)
    args = parser.parse_args()
    if re.fullmatch(r"i-[0-9a-f]{8,17}", args.instance_id) is None:
        raise SystemExit("unsafe instance id")
    if re.fullmatch(r"[a-z]{2}(?:-gov)?-[a-z]+-[0-9]", args.region) is None:
        raise SystemExit("unsafe AWS region")
    print(
        json.dumps(
            build_attestation(args),
            sort_keys=True,
            separators=(",", ":"),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

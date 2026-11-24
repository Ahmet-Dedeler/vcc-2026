#!/usr/bin/env python3
"""Verify every real-data dependency against experiment 083's manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Mapping


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def verify(path: Path, expected: str, label: str) -> None:
    actual = sha256_file(path)
    if actual != expected:
        raise ValueError(f"{label}: SHA-256 {actual} != {expected}")


def verify_group(root: Path, files: Mapping[str, str], label: str) -> None:
    for filename, expected in files.items():
        verify(root / filename, expected, f"{label}/{filename}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--seed-table", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--donor-root", type=Path, required=True)
    parser.add_argument("--perturbmatch", type=Path, required=True)
    parser.add_argument("--kolf", type=Path, required=True)
    parser.add_argument("--cd4-root", type=Path, required=True)
    parser.add_argument("--operator-dir", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    if manifest.get("schema_version") != 1:
        raise ValueError("unsupported manifest schema")
    if args.h5ad.stat().st_size != manifest["guide_source"]["bytes"]:
        raise ValueError("guide H5AD byte count mismatch")
    verify(
        args.h5ad,
        manifest["guide_source"]["sha256"],
        "guide source H5AD",
    )
    if args.seed_table.stat().st_size != manifest["seed_metadata"]["bytes"]:
        raise ValueError("seed table byte count mismatch")
    verify(
        args.seed_table,
        manifest["seed_metadata"]["sha256"],
        "seed metadata",
    )
    verify(
        args.genes_csv,
        manifest["challenge"]["gene_names.csv"],
        "gene_names.csv",
    )
    verify(
        args.panel_csv,
        manifest["challenge"]["pert_counts.csv"],
        "pert_counts.csv",
    )
    for donor in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        verify_group(
            args.donor_root / donor,
            manifest["donors"][donor],
            donor,
        )
    verify_group(args.kolf, manifest["donors"]["kolf"], "kolf")
    verify_group(args.perturbmatch, manifest["perturbmatch"], "perturbmatch")
    for donor, files in manifest["cd4"].items():
        verify_group(args.cd4_root / donor, files, donor)
    verify_group(args.operator_dir, manifest["operators"], "operators")
    print("all experiment 083 inputs match the pinned manifest")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

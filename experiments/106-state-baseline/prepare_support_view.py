#!/usr/bin/env python3
"""Create a copy-on-write STATE support view with consistent log1p metadata."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import h5py


TRAINING_FILES = (
    "competition_train.h5",
    "k562_gwps.h5",
    "rpe1.h5",
    "jurkat.h5",
    "k562.h5",
    "hepg2.h5",
)
PASSTHROUGH_FILES = (
    "competition_val_template.h5ad",
    "competition_test_template.h5ad",
    "starter.toml",
    "gene_names.csv",
    "ESM2_pert_features.pt",
)


def clone_file(source: Path, destination: Path) -> None:
    """Make an APFS clone and fail instead of silently duplicating tens of GiB."""
    subprocess.run(["cp", "-c", "-p", source, destination], check=True)


def ensure_log1p_marker(path: Path) -> bool:
    """Add the standard empty AnnData uns/log1p mapping if it is absent."""
    with h5py.File(path, "r+") as handle:
        uns = handle.require_group("uns")
        if "log1p" in uns:
            return False
        marker = uns.create_group("log1p")
        marker.attrs["encoding-type"] = "dict"
        marker.attrs["encoding-version"] = "0.1.0"
        return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()

    source = args.source.resolve()
    destination = args.destination.resolve()
    manifest_path = args.manifest.resolve()
    expected = set(TRAINING_FILES + PASSTHROUGH_FILES)
    observed = {path.name for path in source.iterdir() if path.is_file()}
    if observed != expected:
        raise SystemExit(
            f"support file set mismatch: missing={sorted(expected - observed)}, "
            f"extra={sorted(observed - expected)}"
        )
    if destination.exists():
        raise SystemExit(f"destination already exists: {destination}")

    staging = Path(
        tempfile.mkdtemp(prefix=f".{destination.name}.preparing-", dir=destination.parent)
    )
    records: list[dict[str, object]] = []
    try:
        for name in TRAINING_FILES + PASSTHROUGH_FILES:
            src = source / name
            dst = staging / name
            clone_file(src, dst)
            before = None
            added = False
            if name in TRAINING_FILES:
                with h5py.File(src, "r") as handle:
                    before = "uns" in handle and "log1p" in handle["uns"]
                added = ensure_log1p_marker(dst)
                with h5py.File(dst, "r") as handle:
                    if "uns" not in handle or "log1p" not in handle["uns"]:
                        raise RuntimeError(f"log1p marker missing after repair: {dst}")
            records.append(
                {
                    "name": name,
                    "source_bytes": src.stat().st_size,
                    "destination_bytes": dst.stat().st_size,
                    "log1p_before": before,
                    "log1p_added": added,
                }
            )
        staging.rename(destination)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise

    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source": str(source),
        "destination": str(destination),
        "copy_method": "macOS cp -c clonefile",
        "semantic_change": "add missing empty AnnData uns/log1p metadata to training files only",
        "files": records,
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

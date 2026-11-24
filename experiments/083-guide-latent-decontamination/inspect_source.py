#!/usr/bin/env python3
"""Read-only structural audit for the local Replogle K562 GWPS H5AD."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import BinaryIO

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_H5AD = ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
DEFAULT_PANEL = ROOT / "data/challenge/pert_counts.csv"
DEFAULT_PERTURBMATCH_PERTS = ROOT / "data/external/perturbmatch/perts.txt"
SOURCE_COMMIT = "198fe3d84263002e9cfa17e1a2178d4e8dd0a669"
SEED_TABLE_NAME = "replogle_k562_2kb_filt_no_lfc_cutoff_off_target_candidates.csv.gz"
SEED_TABLE_URL = (
    "https://raw.githubusercontent.com/AustinHartman/"
    f"perturb_seed_reproduce/{SOURCE_COMMIT}/data/{SEED_TABLE_NAME}"
)
SEED_TABLE_BYTES = 14_460_423
SEED_TABLE_SHA256 = "820f9e5ac81bbdbc1700e4520d7f95b0da6034f7640da1f9c0577e535454db8d"
MIN_PANEL_MULTIGUIDE = 100


def decode(values: np.ndarray) -> tuple[str, ...]:
    return tuple(
        value.decode() if isinstance(value, bytes) else str(value) for value in values
    )


def read_panel(path: Path) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "target_gene" not in reader.fieldnames:
            raise ValueError(f"{path}: missing target_gene")
        values = tuple(str(row["target_gene"]).strip() for row in reader)
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate targets")
    return values


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


class HashingReader:
    def __init__(self, raw: BinaryIO, maximum: int):
        self.raw = raw
        self.maximum = maximum
        self.count = 0
        self.digest = hashlib.sha256()

    def read(self, size: int = -1) -> bytes:
        value = self.raw.read(size)
        self.count += len(value)
        if self.count > self.maximum:
            raise ValueError("seed table exceeded byte cap")
        self.digest.update(value)
        return value

    def readable(self) -> bool:
        return True


def seed_guide_ids(
    *,
    local_path: Path | None,
    allow_network: bool,
) -> tuple[set[str], dict[str, object]]:
    if local_path is not None:
        raw: BinaryIO = local_path.open("rb")
        source = str(local_path.resolve())
    elif allow_network:
        request = urllib.request.Request(
            SEED_TABLE_URL, headers={"User-Agent": "vcc-exp083/1.0"}
        )
        raw = urllib.request.urlopen(request, timeout=120)
        source = SEED_TABLE_URL
    else:
        return set(), {
            "inspected": False,
            "reason": "no local seed table and --allow-network not set",
            "expected_sha256": SEED_TABLE_SHA256,
            "repository_commit": SOURCE_COMMIT,
        }
    hashing = HashingReader(raw, SEED_TABLE_BYTES + 1)
    guides: set[str] = set()
    rows = 0
    try:
        with gzip.GzipFile(fileobj=hashing) as compressed:
            text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
            reader = csv.DictReader(text)
            required = {"guide_id", "protospacer_sequence", "seed"}
            if reader.fieldnames is None or not required <= set(reader.fieldnames):
                raise ValueError("seed table schema mismatch")
            for row in reader:
                rows += 1
                if row["protospacer_sequence"] and row["seed"]:
                    guides.add(row["guide_id"])
    finally:
        raw.close()
    digest = hashing.digest.hexdigest()
    if hashing.count != SEED_TABLE_BYTES or digest != SEED_TABLE_SHA256:
        raise ValueError(
            f"seed table provenance mismatch: bytes={hashing.count}, sha256={digest}"
        )
    return guides, {
        "inspected": True,
        "source": source,
        "repository_commit": SOURCE_COMMIT,
        "compressed_bytes": hashing.count,
        "compressed_sha256": digest,
        "rows": rows,
        "guide_ids_with_sequence": len(guides),
    }


def audit(
    h5ad: Path,
    panel_csv: Path,
    *,
    seed_table: Path | None,
    allow_network: bool,
    hash_h5ad: bool,
    min_cells: int,
    perturbmatch_perts: Path | None,
) -> dict[str, object]:
    panel = read_panel(panel_csv)
    with h5py.File(h5ad, "r") as handle:
        if not isinstance(handle["X"], h5py.Dataset):
            raise ValueError("expected scPerturb dense HDF5 X dataset")
        x = handle["X"]
        targets = decode(handle["obs/perturbation/categories"][:])
        guides = decode(handle["obs/guide_id/categories"][:])
        target_codes_ds = handle["obs/perturbation/codes"]
        guide_codes_ds = handle["obs/guide_id/codes"]
        if (
            len(target_codes_ds) != len(guide_codes_ds)
            or len(target_codes_ds) != x.shape[0]
        ):
            raise ValueError("obs/X row mismatch")
        guide_target = np.full(len(guides), -1, dtype=np.int32)
        guide_cells = np.zeros(len(guides), dtype=np.int64)
        ambiguous: set[int] = set()
        block = 250_000
        for lo in range(0, len(target_codes_ds), block):
            hi = min(lo + block, len(target_codes_ds))
            target_codes = np.asarray(target_codes_ds[lo:hi], dtype=np.int32)
            qcodes = np.asarray(guide_codes_ds[lo:hi], dtype=np.int32)
            valid = (target_codes >= 0) & (qcodes >= 0)
            if not np.any(valid):
                continue
            guide_cells += np.bincount(qcodes[valid], minlength=len(guides))
            pairs = np.unique(
                np.stack((qcodes[valid], target_codes[valid]), axis=1), axis=0
            )
            for guide_code, target_code in pairs:
                current = int(guide_target[guide_code])
                if current >= 0 and current != int(target_code):
                    ambiguous.add(int(guide_code))
                else:
                    guide_target[guide_code] = int(target_code)
        if ambiguous:
            raise ValueError(f"{len(ambiguous)} guides map to multiple targets")
        target_to_guides: dict[str, list[int]] = defaultdict(list)
        for guide_code, target_code in enumerate(guide_target):
            if target_code >= 0 and guide_cells[guide_code] >= min_cells:
                target_to_guides[targets[int(target_code)]].append(guide_code)
        seed_ids, seed_meta = seed_guide_ids(
            local_path=seed_table, allow_network=allow_network
        )
        identifiable_by_target = {
            target: sum(guides[index] in seed_ids for index in rows)
            for target, rows in target_to_guides.items()
        }
        guide_histogram: dict[str, int] = defaultdict(int)
        perturbation_histogram: dict[str, int] = defaultdict(int)
        for rows in target_to_guides.values():
            guide_histogram[str(len(rows))] += 1
        for target, rows in target_to_guides.items():
            if target != "control":
                perturbation_histogram[str(len(rows))] += 1
        panel_covered = [target for target in panel if target in target_to_guides]
        panel_multiguide = [
            target for target in panel_covered if len(target_to_guides[target]) >= 2
        ]
        panel_identifiable = [
            target
            for target in panel_multiguide
            if identifiable_by_target.get(target, 0) >= 2
        ]
        if perturbmatch_perts is not None and perturbmatch_perts.exists():
            pm_targets = set(
                line.strip()
                for line in perturbmatch_perts.read_text().splitlines()
                if line.strip()
            )
            arm_d_evaluable = [
                target for target in panel_identifiable if target in pm_targets
            ]
            pm_metadata: dict[str, object] = {
                "path": str(perturbmatch_perts.resolve()),
                "sha256": sha256_file(perturbmatch_perts),
                "targets": len(pm_targets),
            }
        else:
            arm_d_evaluable = []
            pm_metadata = {
                "path": None,
                "reason": "PerturbMatch target list unavailable",
            }
        nonpanel_identifiable = [
            target
            for target in target_to_guides
            if target not in set(panel)
            and len(target_to_guides[target]) >= 2
            and identifiable_by_target.get(target, 0) >= 2
        ]
        result: dict[str, object] = {
            "experiment": "083-guide-latent-decontamination",
            "audit_kind": "read_only_h5ad_structure",
            "h5ad": {
                "path": str(h5ad.resolve()),
                "bytes": h5ad.stat().st_size,
                "sha256": sha256_file(h5ad) if hash_h5ad else None,
                "cells": int(x.shape[0]),
                "genes": int(x.shape[1]),
                "dtype": str(x.dtype),
                "chunks": list(x.chunks or ()),
                "compression": x.compression,
                "guide_categories": len(guides),
                "target_categories": len(targets),
                "ambiguous_guide_targets": 0,
                "minimum_cells_per_guide": min_cells,
            },
            "guide_structure": {
                "targets_with_cells": len(target_to_guides),
                "guide_count_histogram": dict(
                    sorted(guide_histogram.items(), key=lambda item: int(item[0]))
                ),
                "perturbation_guide_count_histogram": dict(
                    sorted(
                        perturbation_histogram.items(),
                        key=lambda item: int(item[0]),
                    )
                ),
                "multi_guide_targets": sum(
                    target != "control" and len(rows) >= 2
                    for target, rows in target_to_guides.items()
                ),
                "panel_targets": len(panel),
                "panel_covered_targets": len(panel_covered),
                "panel_multi_guide_targets": len(panel_multiguide),
                "panel_identifiable_multi_guide_targets": len(panel_identifiable),
                "panel_arm_d_evaluable_identifiable_multi_guide_targets": len(
                    arm_d_evaluable
                ),
                "nonpanel_identifiable_multi_guide_targets": len(nonpanel_identifiable),
                "panel_multi_guide_target_names": panel_multiguide,
                "panel_identifiable_multi_guide_target_names": panel_identifiable,
                "panel_arm_d_evaluable_identifiable_multi_guide_target_names": (
                    arm_d_evaluable
                ),
            },
            "seed_metadata": seed_meta,
            "perturbmatch_metadata": pm_metadata,
        }
    eligible = int(
        result["guide_structure"][  # type: ignore[index]
            "panel_arm_d_evaluable_identifiable_multi_guide_targets"
        ]
    )
    result["promotion_preflight"] = {
        "minimum_panel_multi_guide_targets": MIN_PANEL_MULTIGUIDE,
        "observed_panel_arm_d_evaluable_identifiable_multi_guide_targets": eligible,
        "passes": eligible >= MIN_PANEL_MULTIGUIDE,
        "heavy_compute_authorized_by_structure": eligible >= MIN_PANEL_MULTIGUIDE,
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--panel-csv", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--seed-table", type=Path)
    parser.add_argument("--allow-network", action="store_true")
    parser.add_argument("--hash-h5ad", action="store_true")
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument(
        "--perturbmatch-perts",
        type=Path,
        default=DEFAULT_PERTURBMATCH_PERTS,
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-promotion-coverage", action="store_true")
    args = parser.parse_args()
    result = audit(
        args.h5ad,
        args.panel_csv,
        seed_table=args.seed_table,
        allow_network=args.allow_network,
        hash_h5ad=args.hash_h5ad,
        min_cells=args.min_cells,
        perturbmatch_perts=args.perturbmatch_perts,
    )
    payload = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    print(payload, end="")
    if args.require_promotion_coverage and not result["promotion_preflight"]["passes"]:
        return 78
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

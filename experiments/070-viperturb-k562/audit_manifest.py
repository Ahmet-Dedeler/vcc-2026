#!/usr/bin/env python3
"""Validate the official VIPerturb genome-wide manifest against the VCC panel."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

EXPECTED_FILTERED_COVERAGE = 124
EXPECTED_SOURCE_PERTURBATIONS = 6_724
LOW_CELL_THRESHOLD = 30
CONTROL_PATTERN = re.compile(
    r"(^|[^a-z0-9])"
    r"(non?.?target|non?.?targeting|ntc|negative|control|safe.?harbou?r|nt)"
    r"([^a-z0-9]|$)"
)


def read_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: missing {column!r} column")
        values = tuple(row[column].strip() for row in reader if row[column].strip())
    if not values or len(values) != len(set(values)):
        raise ValueError(f"{path}: {column!r} is empty or contains duplicates")
    return values


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "t", "1"}:
        return True
    if normalized in {"false", "f", "0"}:
        return False
    raise ValueError(f"invalid boolean value {value!r}")


def read_manifest(path: Path) -> list[dict[str, object]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
        required = {"gene", "cell_count", "file", "is_in_filtered_obj"}
        if reader.fieldnames is None or set(reader.fieldnames) != required:
            raise ValueError(
                f"{path}: fields {reader.fieldnames!r} do not match "
                f"{sorted(required)!r}"
            )
        rows: list[dict[str, object]] = []
        for line_number, row in enumerate(reader, start=2):
            gene = row["gene"].strip()
            if not gene:
                raise ValueError(f"{path}:{line_number}: empty gene")
            cell_count = int(row["cell_count"])
            if cell_count < 0:
                raise ValueError(f"{path}:{line_number}: negative cell count")
            rows.append(
                {
                    "gene": gene,
                    "cell_count": cell_count,
                    "file": row["file"].strip(),
                    "is_in_filtered_obj": parse_bool(row["is_in_filtered_obj"]),
                }
            )
    if not rows:
        raise ValueError(f"{path}: empty manifest")
    return rows


def filtered_source_catalogue(
    manifest: list[dict[str, object]],
    *,
    expected_source_perturbations: int = EXPECTED_SOURCE_PERTURBATIONS,
) -> dict[str, object]:
    filtered = [row for row in manifest if bool(row["is_in_filtered_obj"])]
    labels = [str(row["gene"]).strip() for row in filtered]
    if not labels or any(not label for label in labels):
        raise ValueError("filtered source catalogue has empty perturbation labels")
    control_flags = [
        CONTROL_PATTERN.search(label.lower()) is not None for label in labels
    ]
    control_rows = [
        row
        for row, is_control in zip(filtered, control_flags, strict=True)
        if is_control
    ]
    if not control_rows:
        raise ValueError("filtered source catalogue has no control rows")
    control_labels = tuple(
        dict.fromkeys(
            label
            for label, is_control in zip(labels, control_flags, strict=True)
            if is_control
        )
    )
    if len(control_labels) != 1:
        raise ValueError(
            "filtered source catalogue has ambiguous control labels: "
            + ",".join(control_labels)
        )
    targets = tuple(
        label
        for label, is_control in zip(labels, control_flags, strict=True)
        if not is_control
    )
    seen: set[str] = set()
    duplicate_targets: list[str] = []
    for target in targets:
        if target in seen and target not in duplicate_targets:
            duplicate_targets.append(target)
        seen.add(target)
    if not targets or duplicate_targets:
        raise ValueError(
            "filtered source catalogue has empty or duplicated perturbations: "
            + ",".join(duplicate_targets[:10])
        )
    if len(targets) != expected_source_perturbations:
        raise ValueError(
            f"filtered source perturbation count {len(targets)} != "
            f"{expected_source_perturbations}"
        )
    return {
        "targets": targets,
        "control_labels": control_labels,
        "control_rows": len(control_rows),
    }


def audit_manifest(
    manifest_path: Path,
    panel_path: Path,
    *,
    expected_filtered_coverage: int = EXPECTED_FILTERED_COVERAGE,
    expected_source_perturbations: int = EXPECTED_SOURCE_PERTURBATIONS,
    low_cell_threshold: int = LOW_CELL_THRESHOLD,
) -> dict[str, object]:
    if low_cell_threshold < 1:
        raise ValueError("low-cell threshold must be positive")
    panel = read_column(panel_path, "target_gene")
    manifest = read_manifest(manifest_path)
    source_catalogue = filtered_source_catalogue(
        manifest,
        expected_source_perturbations=expected_source_perturbations,
    )
    panel_set = set(panel)
    panel_rows: dict[str, dict[str, object]] = {}
    for row in manifest:
        gene = str(row["gene"])
        if gene not in panel_set:
            continue
        if gene in panel_rows:
            raise ValueError(f"{manifest_path}: duplicate VCC target row {gene!r}")
        panel_rows[gene] = row
    filtered = [
        panel_rows[target]
        for target in panel
        if target in panel_rows and bool(panel_rows[target]["is_in_filtered_obj"])
    ]
    if len(filtered) != expected_filtered_coverage:
        raise ValueError(
            f"manifest filtered VCC coverage {len(filtered)} != "
            f"{expected_filtered_coverage}"
        )
    low = [row for row in filtered if int(row["cell_count"]) < low_cell_threshold]
    bins = {
        "lt_10": sum(int(row["cell_count"]) < 10 for row in filtered),
        "10_to_19": sum(10 <= int(row["cell_count"]) < 20 for row in filtered),
        "20_to_29": sum(20 <= int(row["cell_count"]) < 30 for row in filtered),
        "30_to_39": sum(30 <= int(row["cell_count"]) < 40 for row in filtered),
        "ge_40": sum(int(row["cell_count"]) >= 40 for row in filtered),
    }
    return {
        "manifest_rows": len(manifest),
        "source_perturbations": len(source_catalogue["targets"]),
        "manifest_control_labels": list(source_catalogue["control_labels"]),
        "manifest_control_rows": source_catalogue["control_rows"],
        "panel_targets": len(panel),
        "filtered_vcc_coverage": len(filtered),
        "filtered_targets": [str(row["gene"]) for row in filtered],
        "cell_count": {
            "minimum": min(int(row["cell_count"]) for row in filtered),
            "maximum": max(int(row["cell_count"]) for row in filtered),
            "bins": bins,
        },
        "low_cell_threshold": low_cell_threshold,
        "low_cell_count": len(low),
        "low_cell_targets": [
            {
                "target_gene": str(row["gene"]),
                "manifest_cell_count": int(row["cell_count"]),
            }
            for row in low
        ],
    }


def read_lines(path: Path) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if not values or len(values) != len(set(values)):
        raise ValueError(f"{path}: empty or duplicate labels")
    return values


def add_overlap_preflight(
    result: dict[str, object],
    panel_path: Path,
    *,
    k562_perts_path: Path,
    perturbmatch_perts_path: Path,
    hct116_perts_path: Path,
    hek293t_perts_path: Path,
    minimum_modified_targets: int,
) -> dict[str, object]:
    panel = read_column(panel_path, "target_gene")
    vip = set(map(str, result["filtered_targets"]))
    k562 = set(read_lines(k562_perts_path))
    perturbmatch = set(read_lines(perturbmatch_perts_path))
    hct116 = set(read_lines(hct116_perts_path))
    hek293t = set(read_lines(hek293t_perts_path))
    gallery = tuple(
        target for target in panel if target in hct116 and target in hek293t
    )
    modified = tuple(
        target
        for target in gallery
        if target in vip and target in k562 and target in perturbmatch
    )
    if len(modified) < minimum_modified_targets:
        raise ValueError(
            f"expected modified target count {len(modified)} < "
            f"{minimum_modified_targets}"
        )
    result["overlap_preflight"] = {
        "full_held_gallery_targets": len(gallery),
        "expected_modified_targets": len(modified),
        "minimum_modified_targets": minimum_modified_targets,
        "modified_targets": list(modified),
    }
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--panel-csv", required=True, type=Path)
    parser.add_argument(
        "--expected-filtered-coverage",
        type=int,
        default=EXPECTED_FILTERED_COVERAGE,
    )
    parser.add_argument(
        "--expected-source-perturbations",
        type=int,
        default=EXPECTED_SOURCE_PERTURBATIONS,
    )
    parser.add_argument("--low-cell-threshold", type=int, default=LOW_CELL_THRESHOLD)
    parser.add_argument("--k562-perts", type=Path)
    parser.add_argument("--perturbmatch-perts", type=Path)
    parser.add_argument("--hct116-perts", type=Path)
    parser.add_argument("--hek293t-perts", type=Path)
    parser.add_argument("--minimum-modified-targets", type=int, default=70)
    parser.add_argument("--output", type=Path)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    result = audit_manifest(
        args.manifest,
        args.panel_csv,
        expected_filtered_coverage=args.expected_filtered_coverage,
        expected_source_perturbations=args.expected_source_perturbations,
        low_cell_threshold=args.low_cell_threshold,
    )
    overlap_paths = (
        args.k562_perts,
        args.perturbmatch_perts,
        args.hct116_perts,
        args.hek293t_perts,
    )
    if any(path is not None for path in overlap_paths):
        if any(path is None for path in overlap_paths):
            raise ValueError("all four donor perturbation files are required")
        result = add_overlap_preflight(
            result,
            args.panel_csv,
            k562_perts_path=args.k562_perts,
            perturbmatch_perts_path=args.perturbmatch_perts,
            hct116_perts_path=args.hct116_perts,
            hek293t_perts_path=args.hek293t_perts,
            minimum_modified_targets=args.minimum_modified_targets,
        )
    payload = json.dumps(result, indent=2) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

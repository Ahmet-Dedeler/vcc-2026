#!/usr/bin/env python3
"""Audit Feng cell metadata without loading expression or guide-count matrices."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import resource
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
FILES = {
    "genomewide_fitness": {
        "name": "GenomeWideScreen_FitnessGenes_Cell-Metadata.tsv.gz",
        "file_id": 51_051_803,
        "bytes": 3_607_490,
        "md5": "15423c32d9122b2999bb7d657b8f8d18",
    },
    "genomewide_nonfitness": {
        "name": "GenomeWideScreen_NonFitnessGenes_Cell-Metadata.tsv.gz",
        "file_id": 51_051_806,
        "bytes": 4_238_782,
        "md5": "9f27a55d488d8ffae93d112d307c263b",
    },
    "targeted": {
        "name": "TargetedScreen_Cell-Metadata.tsv.gz",
        "file_id": 51_051_809,
        "bytes": 13_491_635,
        "md5": "980a440a8132a20567c4c3610e6a5b20",
    },
}


def peak_rss_gib() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def digest(path: Path) -> str:
    value = hashlib.md5()
    with path.open("rb") as handle:
        while block := handle.read(1 << 20):
            value.update(block)
    return value.hexdigest()


def read_panel(path: Path) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "target_gene" not in reader.fieldnames:
            raise ValueError(f"{path}: missing target_gene")
        panel = tuple(row["target_gene"].strip() for row in reader)
    if len(panel) != len(set(panel)):
        raise ValueError("VCC panel contains duplicate targets")
    return panel


def target_from_guide(call: str) -> str | None:
    if call == "unassigned" or call.startswith("NonTarget"):
        return None
    if "_" not in call:
        raise ValueError(f"malformed assigned guide call {call!r}")
    target, sequence = call.split("_", 1)
    if not target or not sequence:
        raise ValueError(f"malformed assigned guide call {call!r}")
    return target


def distribution(values: list[int]) -> dict[str, Any]:
    if not values:
        return {
            "n": 0,
            "min": None,
            "median": None,
            "max": None,
            "at_least_10": 0,
            "at_least_25": 0,
            "at_least_50": 0,
            "at_least_100": 0,
            "at_least_150": 0,
        }
    return {
        "n": len(values),
        "min": min(values),
        "median": statistics.median(values),
        "max": max(values),
        "at_least_10": sum(value >= 10 for value in values),
        "at_least_25": sum(value >= 25 for value in values),
        "at_least_50": sum(value >= 50 for value in values),
        "at_least_100": sum(value >= 100 for value in values),
        "at_least_150": sum(value >= 150 for value in values),
    }


def summarize(
    name: str,
    path: Path,
    panel: tuple[str, ...],
) -> tuple[dict[str, Any], Counter[str], dict[str, set[str]]]:
    specification = FILES[name]
    actual_size = path.stat().st_size
    actual_md5 = digest(path)
    if actual_size != specification["bytes"]:
        raise ValueError(f"{path}: size {actual_size} != {specification['bytes']}")
    if actual_md5 != specification["md5"]:
        raise ValueError(f"{path}: MD5 {actual_md5} != {specification['md5']}")

    statuses: Counter[str] = Counter()
    targets: Counter[str] = Counter()
    target_lines: dict[str, set[str]] = defaultdict(set)
    line_target_cells: Counter[tuple[str, str]] = Counter()
    lines: set[str] = set()

    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_header = ["Cell_ID", "Batch", "Guide_Call", "Cell_Line"]
        if reader.fieldnames != expected_header:
            raise ValueError(f"{path}: unexpected header {reader.fieldnames!r}")
        for row in reader:
            statuses["cells_total"] += 1
            line = row["Cell_Line"]
            call = row["Guide_Call"]
            lines.add(line)
            if call == "unassigned":
                statuses["unassigned"] += 1
                continue
            if call.startswith("NonTarget"):
                statuses["non_targeting"] += 1
                continue
            target = target_from_guide(call)
            assert target is not None
            statuses["assigned_targeting"] += 1
            targets[target] += 1
            target_lines[target].add(line)
            line_target_cells[(line, target)] += 1

    panel_set = set(panel)
    overlap = [target for target in panel if target in targets]
    panel_counts = [targets[target] for target in overlap]
    line_target_counts = [
        count
        for (line, target), count in line_target_cells.items()
        if target in panel_set
    ]
    result = {
        "screen": name,
        "source": {
            "figshare_file_id": specification["file_id"],
            "file_name": specification["name"],
            "bytes": actual_size,
            "md5": actual_md5,
        },
        "cells": dict(statuses),
        "cell_lines": {
            "count": len(lines),
            "symbols": sorted(lines),
        },
        "targets": {
            "all_assigned": len(targets),
            "all_assigned_cell_depth": distribution(list(targets.values())),
            "vcc_overlap": len(overlap),
            "vcc_overlap_symbols": overlap,
            "vcc_overlap_cell_depth": distribution(panel_counts),
            "vcc_line_target_cell_depth": distribution(line_target_counts),
            "vcc_cells": {target: targets[target] for target in overlap},
            "vcc_lines_per_target": {
                target: len(target_lines[target]) for target in overlap
            },
        },
    }
    return result, targets, target_lines


def merged_genomewide(
    panel: tuple[str, ...],
    summaries: dict[str, dict[str, Any]],
    counters: dict[str, Counter[str]],
    line_sets: dict[str, dict[str, set[str]]],
    effect_targets: set[str],
) -> dict[str, Any]:
    merged_counts = counters["genomewide_fitness"] + counters["genomewide_nonfitness"]
    merged_lines: dict[str, set[str]] = defaultdict(set)
    for screen in ("genomewide_fitness", "genomewide_nonfitness"):
        for target, lines in line_sets[screen].items():
            merged_lines[target].update(lines)

    panel_metadata_overlap = [target for target in panel if target in merged_counts]
    panel_effect_overlap = [target for target in panel if target in effect_targets]
    fitness_set = set(counters["genomewide_fitness"])
    nonfitness_set = set(counters["genomewide_nonfitness"])
    return {
        "cells_total": sum(
            summaries[name]["cells"]["cells_total"]
            for name in ("genomewide_fitness", "genomewide_nonfitness")
        ),
        "assigned_targeting": sum(
            summaries[name]["cells"]["assigned_targeting"]
            for name in ("genomewide_fitness", "genomewide_nonfitness")
        ),
        "cell_lines_union": sorted(
            set(summaries["genomewide_fitness"]["cell_lines"]["symbols"])
            | set(summaries["genomewide_nonfitness"]["cell_lines"]["symbols"])
        ),
        "metadata_targets": len(merged_counts),
        "vcc_metadata_overlap": len(panel_metadata_overlap),
        "vcc_effect_overlap_min_10_cells": len(panel_effect_overlap),
        "vcc_effect_overlap_symbols": panel_effect_overlap,
        "vcc_effect_depth": distribution(
            [merged_counts[target] for target in panel_effect_overlap]
        ),
        "vcc_effect_lines_per_target": distribution(
            [len(merged_lines[target]) for target in panel_effect_overlap]
        ),
        "vcc_effect_membership": {
            "fitness_only": [
                target
                for target in panel_effect_overlap
                if target in fitness_set and target not in nonfitness_set
            ],
            "nonfitness_only": [
                target
                for target in panel_effect_overlap
                if target in nonfitness_set and target not in fitness_set
            ],
            "both": [
                target
                for target in panel_effect_overlap
                if target in fitness_set and target in nonfitness_set
            ],
        },
        "vcc_effect_cells": {
            target: merged_counts[target] for target in panel_effect_overlap
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metadata-dir",
        type=Path,
        default=ROOT / "data/external/feng-ipsc",
    )
    parser.add_argument(
        "--panel-csv",
        type=Path,
        default=ROOT / "data/challenge/pert_counts.csv",
    )
    parser.add_argument(
        "--effect-meta",
        type=Path,
        default=(ROOT / "data/external/feng-ipsc/genomewide_vcc_panel/meta.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "data/external/feng-ipsc/metadata_coverage.json",
    )
    parser.add_argument("--max-rss-gib", type=float, default=0.90)
    args = parser.parse_args()

    panel = read_panel(args.panel_csv)
    effect_meta = json.loads(args.effect_meta.read_text())
    effect_targets = set(effect_meta["selection"]["covered_target_symbols"])

    summaries: dict[str, dict[str, Any]] = {}
    counters: dict[str, Counter[str]] = {}
    line_sets: dict[str, dict[str, set[str]]] = {}
    for name, specification in FILES.items():
        summary, counter, per_target_lines = summarize(
            name, args.metadata_dir / specification["name"], panel
        )
        summaries[name] = summary
        counters[name] = counter
        line_sets[name] = per_target_lines

    result = {
        "format_version": 1,
        "source_article": "https://doi.org/10.6084/m9.figshare.27989294",
        "screens": summaries,
        "genomewide_merged": merged_genomewide(
            panel, summaries, counters, line_sets, effect_targets
        ),
        "peak_rss_gib": peak_rss_gib(),
    }
    if result["peak_rss_gib"] >= args.max_rss_gib:
        raise MemoryError(
            f"peak RSS {result['peak_rss_gib']:.3f} GiB exceeded "
            f"{args.max_rss_gib:.3f} GiB"
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Freeze canonical experiment-109 axes without reading expression values."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import sys
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_H5AD = ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
DEFAULT_PANEL = ROOT / "data/challenge/pert_counts.csv"
DEFAULT_PM_PERTS = ROOT / "data/external/perturbmatch/perts.txt"
DEFAULT_PM_GENES = ROOT / "data/external/perturbmatch/genes.txt"
DEFAULT_SOURCE_AUDIT = (
    ROOT / "experiments/083-guide-latent-decontamination/source-audit.json"
)
DEFAULT_OUTPUT = HERE / "source-axes.json"
MIN_CELLS = 20


def _load_exp083_feature_builder():
    path = ROOT / "experiments/083-guide-latent-decontamination/build_seed_features.py"
    name = "exp109_reused_exp083_feature_builder"
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


EXP083_FEATURES = _load_exp083_feature_builder()


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_axis(values: tuple[str, ...]) -> str:
    return hashlib.sha256(("\n".join(values) + "\n").encode()).hexdigest()


def decode(values: np.ndarray) -> tuple[str, ...]:
    return tuple(
        value.decode() if isinstance(value, bytes) else str(value) for value in values
    )


def read_lines(path: Path) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def read_panel(path: Path) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "target_gene" not in reader.fieldnames:
            raise ValueError(f"{path}: missing target_gene")
        values = tuple(sorted(str(row["target_gene"]).strip() for row in reader))
    if len(values) != 300 or len(values) != len(set(values)):
        raise ValueError("challenge panel must contain 300 unique targets")
    return values


def build_axes(
    *,
    h5ad: Path,
    panel_csv: Path,
    perturbmatch_perts: Path,
    perturbmatch_genes: Path,
    source_audit_path: Path,
    seed_table: Path | None,
    allow_network: bool,
) -> dict[str, object]:
    source_audit = json.loads(source_audit_path.read_text())
    expected_h5ad = source_audit["h5ad"]
    if h5ad.stat().st_size != int(expected_h5ad["bytes"]):
        raise ValueError("K562 source byte count changed")
    if int(expected_h5ad["cells"]) != 1_989_578:
        raise ValueError("K562 source row count changed")
    if int(expected_h5ad["genes"]) != 8_248:
        raise ValueError("K562 source gene count changed")

    panel = read_panel(panel_csv)
    pm_targets = set(read_lines(perturbmatch_perts))
    pm_genes = set(read_lines(perturbmatch_genes))
    with h5py.File(h5ad, "r") as handle:
        if not isinstance(handle["X"], h5py.Dataset):
            raise TypeError("expected dense source X dataset")
        genes = decode(handle["var/gene_name"][:])
        targets = decode(handle["obs/perturbation/categories"][:])
        guides = decode(handle["obs/guide_id/categories"][:])
        target_codes_ds = handle["obs/perturbation/codes"]
        guide_codes_ds = handle["obs/guide_id/codes"]
        if len(genes) != handle["X"].shape[1]:
            raise ValueError("source gene axis and X differ")
        if (
            len(target_codes_ds) != len(guide_codes_ds)
            or len(target_codes_ds) != handle["X"].shape[0]
        ):
            raise ValueError("source observation axes differ")

        guide_target = np.full(len(guides), -1, dtype=np.int32)
        guide_cells = np.zeros(len(guides), dtype=np.int64)
        for lo in range(0, len(target_codes_ds), 250_000):
            hi = min(lo + 250_000, len(target_codes_ds))
            target_codes = np.asarray(target_codes_ds[lo:hi], dtype=np.int32)
            guide_codes = np.asarray(guide_codes_ds[lo:hi], dtype=np.int32)
            valid = (target_codes >= 0) & (guide_codes >= 0)
            if not np.any(valid):
                continue
            guide_cells += np.bincount(guide_codes[valid], minlength=len(guides))
            pairs = np.unique(
                np.stack((guide_codes[valid], target_codes[valid]), axis=1),
                axis=0,
            )
            for guide_code, target_code in pairs:
                current = int(guide_target[guide_code])
                if current >= 0 and current != int(target_code):
                    raise ValueError("one guide maps to multiple targets")
                guide_target[guide_code] = int(target_code)

    if len(genes) != len(set(genes)):
        raise ValueError("source gene axis contains duplicates")
    if seed_table is None and not allow_network:
        raise ValueError("sequence metadata is required to freeze source axes")
    _, feature_identifiable, feature_names, feature_provenance = (
        EXP083_FEATURES.build_features(
            guides,
            local_source=seed_table,
        )
    )
    if len(feature_names) != 229:
        raise ValueError("sequence feature count changed")
    seed_ids = {
        guide for guide, keep in zip(guides, feature_identifiable, strict=True) if keep
    }
    seed_metadata = {
        **feature_provenance,
        "inspected": True,
        "eligibility_predicate": "normalized_valid_protospacer_and_seed",
    }
    target_to_guides: dict[str, list[int]] = defaultdict(list)
    for guide_code, target_code in enumerate(guide_target):
        if target_code >= 0 and guide_cells[guide_code] >= MIN_CELLS:
            target_to_guides[targets[int(target_code)]].append(guide_code)

    panel_set = set(panel)
    identifiable_count = {
        target: sum(guides[index] in seed_ids for index in rows)
        for target, rows in target_to_guides.items()
    }
    development = tuple(
        sorted(
            target
            for target, rows in target_to_guides.items()
            if target != "control"
            and target not in panel_set
            and len(rows) >= 2
            and identifiable_count[target] >= 2
        )
    )
    singleton = tuple(
        sorted(
            target
            for target, rows in target_to_guides.items()
            if target != "control" and len(rows) == 1 and guides[rows[0]] in seed_ids
        )
    )
    response = tuple(
        gene for gene in genes if gene not in panel_set and gene not in set(development)
    )
    panel_sequence = tuple(
        target
        for target in panel
        if target in target_to_guides and identifiable_count[target] >= 1
    )
    panel_sequence_pm = tuple(
        target for target in panel_sequence if target in pm_targets
    )
    panel_sequence_pm_singleton = tuple(
        target for target in panel_sequence_pm if len(target_to_guides[target]) == 1
    )
    panel_sequence_pm_multiguide = tuple(
        target for target in panel_sequence_pm if len(target_to_guides[target]) >= 2
    )

    expected = {
        "development": 506,
        "singleton": 6689,
        "panel_sequence": 223,
        "panel_sequence_pm": 212,
        "panel_sequence_pm_singleton": 206,
        "panel_sequence_pm_multiguide": 6,
    }
    observed = {
        "development": len(development),
        "singleton": len(singleton),
        "panel_sequence": len(panel_sequence),
        "panel_sequence_pm": len(panel_sequence_pm),
        "panel_sequence_pm_singleton": len(panel_sequence_pm_singleton),
        "panel_sequence_pm_multiguide": len(panel_sequence_pm_multiguide),
    }
    if observed != expected:
        raise ValueError(f"source-axis audit counts changed: {observed} != {expected}")

    axes = {
        "gene_axis": genes,
        "panel_targets": panel,
        "development_targets": development,
        "singleton_centering_targets": singleton,
        "response_axis": response,
        "perturbmatch_targets": tuple(sorted(pm_targets)),
        "perturbmatch_genes": tuple(sorted(pm_genes)),
    }
    return {
        "schema_version": 1,
        "experiment": "109-global-guide-nuisance",
        "status": "canonical_source_axes_frozen_before_expression",
        "minimum_cells_per_guide": MIN_CELLS,
        "source": {
            "h5ad_path": str(h5ad.relative_to(ROOT)),
            "h5ad_bytes": h5ad.stat().st_size,
            "h5ad_sha256_from_pinned_source_audit": expected_h5ad["sha256"],
            "source_audit_path": str(source_audit_path.relative_to(ROOT)),
            "source_audit_sha256": sha256_file(source_audit_path),
            "seed_metadata": seed_metadata,
            "panel_csv_sha256": sha256_file(panel_csv),
            "perturbmatch_perts_sha256": sha256_file(perturbmatch_perts),
            "perturbmatch_genes_sha256": sha256_file(perturbmatch_genes),
        },
        "counts": observed,
        "axes": {
            name: {
                "count": len(values),
                "sha256": sha256_axis(values),
                "values": list(values),
            }
            for name, values in axes.items()
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", type=Path, default=DEFAULT_H5AD)
    parser.add_argument("--panel-csv", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--perturbmatch-perts", type=Path, default=DEFAULT_PM_PERTS)
    parser.add_argument("--perturbmatch-genes", type=Path, default=DEFAULT_PM_GENES)
    parser.add_argument("--source-audit", type=Path, default=DEFAULT_SOURCE_AUDIT)
    parser.add_argument("--seed-table", type=Path)
    parser.add_argument("--allow-network", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    result = build_axes(
        h5ad=args.h5ad.resolve(),
        panel_csv=args.panel_csv.resolve(),
        perturbmatch_perts=args.perturbmatch_perts.resolve(),
        perturbmatch_genes=args.perturbmatch_genes.resolve(),
        source_audit_path=args.source_audit.resolve(),
        seed_table=args.seed_table.resolve() if args.seed_table else None,
        allow_network=args.allow_network,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "status": result["status"],
                "output": str(args.output),
                "counts": result["counts"],
                "sha256": sha256_file(args.output),
                "expression_values_opened": False,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

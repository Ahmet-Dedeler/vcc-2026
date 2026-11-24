#!/usr/bin/env python3
"""Freeze Experiment 135 before TeloHAEC panel-response aggregation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "IMPLEMENTATION_FREEZE.json"

FILES = (
    "experiments/135-telohaec-atlas-confirmation/README.md",
    "experiments/135-telohaec-atlas-confirmation/preregistered_protocol.json",
    "experiments/135-telohaec-atlas-confirmation/run.py",
    "experiments/135-telohaec-atlas-confirmation/test_run.py",
    "experiments/135-telohaec-atlas-confirmation/freeze.py",
    "experiments/122-panel-response-atlas/work/source-only-phase1-v1/selected-panel-residual.npy",
    "experiments/122-panel-response-atlas/work/source-only-phase1-v1/selected-panel-confidence.npy",
    "experiments/122-panel-response-atlas/work/source-only-phase1-v1/all-candidate-panel-residuals.npy",
    "experiments/122-panel-response-atlas/work/source-only-phase1-v1/all-candidate-panel-confidence.npy",
    "experiments/122-panel-response-atlas/work/source-only-phase1-v1/all-candidate-panel-support.npy",
    "experiments/122-panel-response-atlas/work/source-only-phase1-v1/candidate-names.txt",
    "experiments/122-panel-response-atlas/work/source-only-phase1-v1/panel-targets.txt",
    "experiments/122-panel-response-atlas/work/source-only-phase1-v1/results.json",
    "experiments/122-panel-response-atlas/work/source-only-phase1-v1/artifact-manifest.json",
    "experiments/122-panel-response-atlas/work/source-only-phase1-v1/result-receipt.json",
    "experiments/128-panel-integer-budget-transport/IMPLEMENTATION_FREEZE.json",
    "experiments/128-panel-integer-budget-transport/run.py",
    "experiments/103-residual-only-production/work/staged-a15-linux/a15-effects.npy",
    "experiments/103-residual-only-production/work/staged-a15-linux/a15-effects.npy.json",
    "experiments/105-telohaec-residual/PROVENANCE.json",
    "experiments/130-telohaec-panel-confirmation/SOURCE_METADATA_RECEIPT.v1.json",
    "data/external/igvf-telohaec/guide_qc.csv",
    "data/challenge/gene_names.csv",
    "data/challenge/pert_counts.csv",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite {OUTPUT}")
    files = {}
    for relative in FILES:
        path = ROOT / relative
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(relative)
        files[relative] = {
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    protocol = json.loads((HERE / "preregistered_protocol.json").read_text())
    source = ROOT / protocol["source"]["matrix_path"]
    if source.stat().st_size != protocol["source"]["expected_matrix_size_bytes"]:
        raise RuntimeError("TeloHAEC matrix size changed")
    provenance = json.loads(
        (ROOT / "experiments/105-telohaec-residual/PROVENANCE.json").read_text()
    )
    recorded = provenance["local_files"][source.name]
    if recorded["sha256_observed"] != protocol["source"]["expected_matrix_sha256"]:
        raise RuntimeError("TeloHAEC provenance hash changed")
    payload = {
        "schema_version": 1,
        "experiment": "135-telohaec-atlas-confirmation",
        "created_at": "2026-08-30",
        "status": "frozen_before_experiment_135_telohaec_panel_response_aggregation_or_analysis",
        "files": files,
        "large_source": {
            "path": protocol["source"]["matrix_path"],
            "size_bytes": source.stat().st_size,
            "sha256_from_verified_experiment_105_provenance": recorded[
                "sha256_observed"
            ],
        },
        "candidate": "equal_mean_agreement",
        "dose": 1.0,
        "prior_raw_chunk_streaming_disclosed": True,
        "experiment_105_excluded_all_panel_coordinates_before_aggregation": True,
        "telohaec_panel_response_aggregation_or_analysis": False,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    OUTPUT.chmod(0o444)
    print(json.dumps({"freeze": str(OUTPUT), "files": len(files)}, indent=2))


if __name__ == "__main__":
    main()

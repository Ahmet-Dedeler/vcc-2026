#!/usr/bin/env python3
"""Build and validate held-context essential-pruned whitening operators."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BUILDER = ROOT / "scripts/build_union_whitener.py"
ESSENTIAL_LIST = ROOT / "data/donors/k562_essential/perts.txt"
DEFAULT_OUTPUT = HERE / "operators"

SPECS = {
    "hct116": {
        "donors": ("k562_gwps", "orion_hek293t"),
        "reference_metadata": (
            ROOT / "experiments/062-cd4-post-whiten/cloud-output/whiteners/"
            "held_orion_hct116.json"
        ),
        "filename": "held_orion_hct116_pruned.npz",
    },
    "hek293t": {
        "donors": ("k562_gwps", "orion_hct116"),
        "reference_metadata": (
            ROOT / "experiments/062-cd4-post-whiten/cloud-output/whiteners/"
            "held_orion_hek293t.json"
        ),
        "filename": "held_orion_hek293t_pruned.npz",
    },
}


def read_axis(path: Path) -> tuple[str, ...]:
    rows = tuple(line.strip().split(",")[0] for line in path.read_text().splitlines())
    if rows and rows[0] in {"gene_name", "target_gene"}:
        rows = rows[1:]
    return tuple(value for value in rows if value)


def expected_operator_axis(donors: tuple[str, ...]) -> np.ndarray:
    official = read_axis(ROOT / "data/challenge/gene_names.csv")
    panel = set(read_axis(ROOT / "data/challenge/pert_counts.csv"))
    donor_union: set[str] = set()
    for donor in donors:
        donor_union.update(read_axis(ROOT / f"data/donors/{donor}/genes.txt"))
    return np.asarray(
        [
            index
            for index, gene in enumerate(official)
            if gene in donor_union and gene not in panel
        ],
        dtype=np.int64,
    )


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def validate_operator(
    path: Path, *, donors: tuple[str, ...], reference_metadata: Path
) -> dict[str, object]:
    metadata_path = path.with_suffix(".json")
    if not path.is_file() or not metadata_path.is_file():
        raise FileNotFoundError(f"incomplete held operator pair: {path}")
    metadata = json.loads(metadata_path.read_text())
    reference = json.loads(reference_metadata.read_text())
    wcols = expected_operator_axis(donors)
    wcols_sha256 = hashlib.sha256(wcols.tobytes()).hexdigest()
    expected = {
        "axis": "union",
        "donors": list(donors),
        "n_cols": int(len(wcols)),
        "gamma": float(reference["gamma"]),
        "lambda": reference.get("lambda"),
        "wcols_sha256": wcols_sha256,
        "basis_prune": "essential",
        "dropped_in_panel": [],
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(
                f"{path}: metadata {key}={metadata.get(key)!r}, expected {value!r}"
            )
    observed = sha256_file(path)
    if metadata.get("sha256") != observed:
        raise ValueError(f"{path}: operator hash differs from metadata")
    if int(reference["n_cols"]) != len(wcols):
        raise ValueError("reference operator column count differs from recomputed axis")
    if str(reference["wcols_sha256"]) != wcols_sha256:
        raise ValueError("reference operator column hash differs from recomputed axis")
    with np.load(path, allow_pickle=False) as archive:
        evals = np.asarray(archive["evals"])
        evecs = np.asarray(archive["evecs"])
    if evals.shape != (len(wcols),):
        raise ValueError(f"{path}: eigenvalue shape {evals.shape} is invalid")
    if evecs.shape != (len(wcols), len(wcols)):
        raise ValueError(f"{path}: eigenvector shape {evecs.shape} is invalid")
    if not np.isfinite(evals).all() or np.any(evals <= 0):
        raise ValueError(f"{path}: eigenvalues must be finite and strictly positive")
    if not np.isfinite(evecs).all():
        raise ValueError(f"{path}: eigenvectors contain non-finite values")
    probe = np.linspace(0, len(wcols) - 1, num=min(32, len(wcols)), dtype=np.int64)
    gram = evecs[:, probe].T @ evecs[:, probe]
    if not np.allclose(gram, np.eye(len(probe)), rtol=5e-5, atol=5e-5):
        raise ValueError(f"{path}: sampled eigenvectors are not orthonormal")
    return metadata


def build_one(
    surface: str,
    output_dir: Path,
    *,
    max_rss_gib: float,
    chunk_rows: int,
) -> dict[str, object]:
    spec = SPECS[surface]
    donors = tuple(spec["donors"])
    reference_metadata = Path(spec["reference_metadata"])
    output = output_dir / str(spec["filename"])
    metadata_path = output.with_suffix(".json")
    if output.is_file() and metadata_path.is_file():
        return validate_operator(
            output,
            donors=donors,
            reference_metadata=reference_metadata,
        )
    if output.exists() or metadata_path.exists():
        output.unlink(missing_ok=True)
        metadata_path.unlink(missing_ok=True)

    output_dir.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{surface}-pruned-", dir=output_dir))
    staged_output = staging / output.name
    try:
        subprocess.run(
            [
                sys.executable,
                str(BUILDER),
                "--bundle",
                str(ROOT / "data/challenge"),
                "--donor-root",
                str(ROOT / "data/donors"),
                "--donors",
                ",".join(donors),
                "--basis-prune",
                "essential",
                "--essential-list",
                str(ESSENTIAL_LIST),
                "--out",
                str(staged_output),
                "--scratch",
                str(staging / "scratch"),
                "--chunk-rows",
                str(chunk_rows),
                "--max-rss-gib",
                str(max_rss_gib),
            ],
            check=True,
        )
        metadata = json.loads(staged_output.with_suffix(".json").read_text())
        reference = json.loads(reference_metadata.read_text())
        wcols = expected_operator_axis(donors)
        wcols_sha256 = hashlib.sha256(wcols.tobytes()).hexdigest()
        if int(metadata.get("n_cols", -1)) != len(wcols):
            raise ValueError("new operator column count differs from recomputed axis")
        if int(reference.get("n_cols", -1)) != len(wcols):
            raise ValueError(
                "reference operator column count differs from recomputed axis"
            )
        if reference.get("wcols_sha256") != wcols_sha256:
            raise ValueError(
                "reference operator column hash differs from recomputed axis"
            )
        for key in ("axis", "gamma", "lambda", "wcols_sha256"):
            metadata[key] = reference.get(key)
        metadata["source_builder"] = str(BUILDER.relative_to(ROOT))
        metadata["source_builder_sha256"] = sha256_file(BUILDER)
        metadata["essential_list_sha256"] = sha256_file(ESSENTIAL_LIST)
        metadata["reference_metadata"] = str(reference_metadata.relative_to(ROOT))
        metadata["reference_metadata_sha256"] = sha256_file(reference_metadata)
        metadata["sha256"] = sha256_file(staged_output)
        atomic_json(staged_output.with_suffix(".json"), metadata)
        os.replace(staged_output, output)
        os.replace(staged_output.with_suffix(".json"), metadata_path)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return validate_operator(
        output,
        donors=donors,
        reference_metadata=reference_metadata,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--max-rss-gib", type=float, default=20.0)
    parser.add_argument("--chunk-rows", type=int, default=2048)
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    results = {}
    for surface in ("hct116", "hek293t"):
        metadata = build_one(
            surface,
            output_dir,
            max_rss_gib=args.max_rss_gib,
            chunk_rows=args.chunk_rows,
        )
        results[surface] = {
            "path": str(
                (output_dir / str(SPECS[surface]["filename"])).relative_to(ROOT)
            ),
            "sha256": metadata["sha256"],
            "donors": metadata["donors"],
            "basis": metadata["basis"],
            "rows_dropped_total": metadata["rows_dropped_total"],
        }
    manifest = {
        "experiment": "106-state-baseline",
        "stage": "held-essential-pruned-operators",
        "selection": "fixed experiment-091 essential-list policy",
        "builder": str(BUILDER.relative_to(ROOT)),
        "builder_sha256": sha256_file(BUILDER),
        "essential_list_sha256": sha256_file(ESSENTIAL_LIST),
        "operators": results,
    }
    atomic_json(output_dir / "manifest.json", manifest)
    print(output_dir / "manifest.json")


if __name__ == "__main__":
    main()

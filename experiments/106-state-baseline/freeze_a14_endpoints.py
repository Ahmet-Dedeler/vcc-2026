#!/usr/bin/env python3
"""Freeze exact held-context a14-shaped effects before STATE scoring."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path
from types import ModuleType

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXP104 = ROOT / "experiments/104-context-heldout-credit-adapter/run.py"
A14_REFERENCE_META = ROOT / "submissions/a15-cd4f010.vcc.meta.json"
A14_REFERENCE_META_SHA256 = (
    "34223442e4e8816bb9be2b1fdfe6a563ed114c17a2893bee0ff603953f84287a"
)
HELPER_SOURCE_SHA256 = {
    "experiments/104-context-heldout-credit-adapter/run.py": (
        "7c0300992db400e5ad6a6e5be2e8d64955f4d00a06ab7d9738f14a8c5ad096dd"
    ),
    "experiments/100-sign-confidence/run.py": (
        "e47c5cff1671ebe49bc5c65f3af7f599c4e0a54907652f636b14c7143aa1492b"
    ),
    "experiments/062-cd4-post-whiten/run.py": (
        "16417a9b985d9aad8993bf2f24ce550ea1e6871d9401bb4672d11ec6d0046176"
    ),
    "experiments/064-pm-cd4-composition/run.py": (
        "39d0f5a29731c5633f3d8fd4035a5a9c3a91cbd435780a201059eaa25a36df81"
    ),
}
A14_ALPHA = 0.982
A14_CD4_FRACTION = 0.05
EXPECTED_CONTROL_ROWS = 18_400
HELD_PRUNED_OPERATORS = {
    "hct116": HERE / "operators/held_orion_hct116_pruned.npz",
    "hek293t": HERE / "operators/held_orion_hek293t_pruned.npz",
}
HELD_HARNESSES = {
    "hct116": ROOT / "data/harness/orion-hct116-thinned.h5ad",
    "hek293t": ROOT / "data/harness/orion-hek293t-thinned.h5ad",
}


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def atomic_npz(path: Path, **values: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with tmp.open("wb") as handle:
        np.savez_compressed(handle, **values)
    tmp.replace(path)


def atomic_json(path: Path, value: object) -> None:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def sha256_strings(values: tuple[str, ...] | list[str] | np.ndarray) -> str:
    return hashlib.sha256(("\n".join(map(str, values)) + "\n").encode()).hexdigest()


def a14_biological_inputs(surface: str) -> dict[str, Path]:
    other = "orion_hek293t" if surface == "hct116" else "orion_hct116"
    keys = ["bundle/gene_names.csv", "bundle/pert_counts.csv"]
    for donor in ("k562_gwps", other):
        keys.extend(
            f"donors/{donor}/{name}"
            for name in (
                "b_ctrl.npy",
                "delta.npy",
                "genes.txt",
                "meta.json",
                "perts.txt",
            )
        )
    keys.extend(
        f"perturbmatch/{name}"
        for name in ("delta.npy", "genes.txt", "meta.json", "perts.txt")
    )
    for condition in ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr"):
        keys.extend(
            f"cd4/{condition}/{name}"
            for name in (
                "b_ctrl.npy",
                "delta.npy",
                "genes.txt",
                "meta.json",
                "perts.txt",
            )
        )

    def resolve(key: str) -> Path:
        prefix, rest = key.split("/", 1)
        roots = {
            "bundle": ROOT / "data/challenge",
            "donors": ROOT / "data/donors",
            "perturbmatch": ROOT / "data/external/perturbmatch",
            "cd4": ROOT / "data/external/cd4t/cloud-compact/cd4_compact",
        }
        return roots[prefix] / rest

    return {key: resolve(key) for key in keys}


def validate_a14_reference(
    surface: str, controls_h5ad: Path | None = None
) -> dict[str, object]:
    if sha256_file(A14_REFERENCE_META) != A14_REFERENCE_META_SHA256:
        raise ValueError("frozen production reference manifest changed")
    reference = json.loads(A14_REFERENCE_META.read_text())
    expected_inputs = reference["hashes"]["inputs"]
    observed: dict[str, str] = {}
    for relative, expected in HELPER_SOURCE_SHA256.items():
        path = ROOT / relative
        actual = sha256_file(path)
        if actual != expected:
            raise ValueError(f"a14 helper source changed: {relative}")
        observed[relative] = actual
    for key, path in a14_biological_inputs(surface).items():
        expected = expected_inputs.get(key)
        if expected is None:
            raise ValueError(f"production reference lacks biological input: {key}")
        actual = sha256_file(path)
        if actual != expected:
            raise ValueError(f"a14 biological input changed: {key}")
        observed[str(path.relative_to(ROOT))] = actual
    held_paths = [
        HELD_HARNESSES[surface],
        HELD_HARNESSES[surface].with_suffix(".h5ad.meta.json"),
        HELD_PRUNED_OPERATORS[surface],
        HELD_PRUNED_OPERATORS[surface].with_suffix(".json"),
        ROOT / "scripts/build_union_whitener.py",
        ROOT / "data/donors/k562_essential/perts.txt",
    ]
    if controls_h5ad is not None:
        held_paths.extend(
            [
                controls_h5ad,
                controls_h5ad.with_suffix(controls_h5ad.suffix + ".manifest.json"),
            ]
        )
    for path in held_paths:
        if not path.is_file():
            raise FileNotFoundError(path)
        observed[str(path.relative_to(ROOT))] = sha256_file(path)
    return {
        "reference_manifest": str(A14_REFERENCE_META.relative_to(ROOT)),
        "reference_manifest_sha256": A14_REFERENCE_META_SHA256,
        "files": observed,
    }


def sum_matrix_rows(
    matrix: object, n_rows: int, n_genes: int, chunk_rows: int = 256
) -> np.ndarray:
    import scipy.sparse as sp

    sums = np.zeros(n_genes, dtype=np.float64)
    for start in range(0, n_rows, chunk_rows):
        block = matrix[start : start + chunk_rows, :]
        if sp.issparse(block):
            sums += np.asarray(block.sum(axis=0)).ravel()
        else:
            sums += np.asarray(block, dtype=np.float64).sum(axis=0)
    return sums


def load_control_profile(
    path: Path, expected_genes: np.ndarray, expected_surface: str
) -> tuple[np.ndarray, np.ndarray, int, dict[str, object]]:
    import anndata as ad

    manifest_path = path.with_suffix(path.suffix + ".manifest.json")
    if not path.is_file() or not manifest_path.is_file():
        raise FileNotFoundError(f"incomplete control-only artifact: {path}")
    manifest = json.loads(manifest_path.read_text())
    expected_manifest = {
        "experiment": "106-state-baseline",
        "stage": "held-control-only-input",
        "contains_perturbation_rows": False,
        "expression_read_policy": (
            "only predicate-filtered control rows are materialized by experiment code"
        ),
    }
    for key, value in expected_manifest.items():
        if manifest.get(key) != value:
            raise ValueError(
                f"{path}: control manifest {key}={manifest.get(key)!r}, "
                f"expected {value!r}"
            )
    if manifest.get("surface") != expected_surface:
        raise ValueError("control-only manifest has the wrong held surface")
    expected_label = {
        "hct116": {"Non-Targeting": 165_777},
        "hek293t": {"Non-Targeting": 218_838},
    }[expected_surface]
    if manifest.get("accepted_control_label_counts") != expected_label:
        raise ValueError("control-only manifest has the wrong accepted label contract")
    if manifest.get("output_sha256") != sha256_file(path):
        raise ValueError("control-only H5AD hash differs from its manifest")

    handle = ad.read_h5ad(path, backed="r")
    try:
        genes = np.asarray(handle.var_names).astype(str)
        if not np.array_equal(genes, np.asarray(expected_genes).astype(str)):
            raise ValueError("control-only and held endpoint gene axes differ")
        labels = np.asarray(handle.obs["target_gene"].astype(str))
        if len(labels) != EXPECTED_CONTROL_ROWS:
            raise ValueError(
                f"control-only input has {len(labels)} rows; "
                f"expected {EXPECTED_CONTROL_ROWS}"
            )
        if np.any(labels != "non-targeting"):
            raise ValueError("control-only input contains a non-control row")
        matrix = handle.X
        sums = sum_matrix_rows(matrix, len(labels), len(genes))
    finally:
        handle.file.close()
    if not np.isfinite(sums).all() or np.any(sums < 0) or sums.sum() <= 0:
        raise ValueError("invalid control-only aggregate")
    bctrl = np.log1p(50_000.0 * sums / sums.sum())
    ctrl_cpm = 1e6 * sums / sums.sum()
    return bctrl, ctrl_cpm, len(labels), manifest


def validate_held_operator(surface: str, path: Path) -> dict[str, object]:
    metadata_path = path.with_suffix(".json")
    if not path.is_file() or not metadata_path.is_file():
        raise FileNotFoundError(f"held pruned operator is missing: {path}")
    metadata = json.loads(metadata_path.read_text())
    expected_donors = {
        "hct116": ["k562_gwps", "orion_hek293t"],
        "hek293t": ["k562_gwps", "orion_hct116"],
    }[surface]
    expected = {
        "axis": "union",
        "donors": expected_donors,
        "gamma": 0.65,
        "basis_prune": "essential",
        "dropped_in_panel": [],
    }
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise ValueError(
                f"{path}: metadata {key}={metadata.get(key)!r}, expected {value!r}"
            )
    if metadata.get("sha256") != sha256_file(path):
        raise ValueError("held pruned operator hash differs from metadata")
    with np.load(path, allow_pickle=False) as archive:
        evals = np.asarray(archive["evals"])
        evecs = np.asarray(archive["evecs"])
    expected_cols = int(metadata["n_cols"])
    if evals.shape != (expected_cols,):
        raise ValueError("held operator eigenvalue shape is malformed")
    if evecs.shape != (expected_cols, expected_cols):
        raise ValueError("held operator eigenvector shape is malformed")
    if not np.isfinite(evals).all() or np.any(evals <= 0):
        raise ValueError("held operator eigenvalues must be finite and positive")
    if not np.isfinite(evecs).all():
        raise ValueError("held operator eigenvectors contain non-finite values")
    probe = np.linspace(
        0, expected_cols - 1, num=min(32, expected_cols), dtype=np.int64
    )
    gram = evecs[:, probe].T @ evecs[:, probe]
    if not np.allclose(gram, np.eye(len(probe)), rtol=5e-5, atol=5e-5):
        raise ValueError("held operator sampled eigenvectors are not orthonormal")
    return metadata


def freeze(
    surface: str,
    controls_h5ad: Path,
    output_dir: Path,
    limit_gib: float,
) -> None:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    provenance_before = validate_a14_reference(surface, controls_h5ad)
    exp104 = load_module("exp106_a14_builder", EXP104)
    if exp104.ALPHA != A14_ALPHA:
        raise ValueError(f"experiment 104 alpha changed: {exp104.ALPHA}")
    exp104.CD4_FRACTION = A14_CD4_FRACTION
    operator_path = HELD_PRUNED_OPERATORS[surface]
    operator_metadata = validate_held_operator(surface, operator_path)
    exp104.SURFACES[surface]["operator"] = operator_path

    genes, targets, _, _, baseline, audit = exp104.build_a15_shaped_baseline(
        surface, limit_gib
    )
    bctrl, ctrl_cpm, control_rows, control_manifest = load_control_profile(
        controls_h5ad, genes, surface
    )
    panel_mask = np.asarray(
        [str(gene) in exp104.PANEL_SET for gene in genes], dtype=bool
    )
    if np.any(baseline[:, panel_mask] != 0):
        raise AssertionError("a14 held endpoint retained a VCC-panel coordinate")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.", dir=output_dir.parent)
    )
    try:
        frozen_path = staging / "frozen.npz"
        atomic_npz(
            frozen_path,
            genes=np.asarray(genes, dtype="U"),
            targets=np.asarray(targets, dtype="U"),
            bctrl=np.asarray(bctrl, dtype=np.float32),
            ctrl_cpm=np.asarray(ctrl_cpm, dtype=np.float32),
            baseline=np.asarray(baseline, dtype=np.float32),
        )
        manifest = {
            "experiment": "106-state-baseline",
            "stage": "held-incumbent-freeze",
            "surface": surface,
            "sealed_before_perturbation_truth": True,
            "endpoint": "a14-shaped",
            "alpha": A14_ALPHA,
            "cd4_fraction": A14_CD4_FRACTION,
            "perturbmatch_fraction": 0.5,
            "whiten_gamma": 0.65,
            "basis_prune": "essential",
            "held_operator": str(operator_path.relative_to(ROOT)),
            "held_operator_sha256": operator_metadata["sha256"],
            "held_operator_metadata_sha256": sha256_file(
                operator_path.with_suffix(".json")
            ),
            "controls_h5ad": str(controls_h5ad.relative_to(ROOT)),
            "controls_h5ad_sha256": sha256_file(controls_h5ad),
            "controls_manifest_sha256": sha256_file(
                controls_h5ad.with_suffix(controls_h5ad.suffix + ".manifest.json")
            ),
            "controls_selection": control_manifest["selection"],
            "control_rows": control_rows,
            "panel_coordinates_zero": True,
            "targets": int(len(targets)),
            "target_order_sha256": sha256_strings(targets),
            "genes": int(len(genes)),
            "gene_axis_sha256": sha256_strings(genes),
            "builder_audit": audit,
            "builder_sha256": sha256_file(EXP104),
            "a14_reference_contract": provenance_before,
            "implementation_sha256": sha256_file(Path(__file__)),
            "frozen_sha256": sha256_file(frozen_path),
        }
        if validate_a14_reference(surface, controls_h5ad) != provenance_before:
            raise RuntimeError("a14 source or biological input changed during freeze")
        atomic_json(staging / "manifest.json", manifest)
        os.replace(staging, output_dir)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    print(output_dir / "manifest.json")


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--surface", choices=("hct116", "hek293t"), required=True)
    p.add_argument("--controls-h5ad", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--limit-gib", type=float, default=24.0)
    return p


def main() -> None:
    args = parser().parse_args()
    freeze(
        args.surface,
        args.controls_h5ad.resolve(),
        args.output_dir.resolve(),
        args.limit_gib,
    )


if __name__ == "__main__":
    main()

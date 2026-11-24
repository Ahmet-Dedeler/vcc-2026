#!/usr/bin/env python3
"""Run Experiment 147's frozen source-only feasibility gate.

This script reads count matrices only from the two prediction artifacts. It
uses the development reference solely for obs labels and the gene axis; it
never accesses ``real.X``.
"""

from __future__ import annotations

import gc
import hashlib
import importlib.util
import json
import resource
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SPEC = importlib.util.spec_from_file_location("dual_moment", HERE / "dual_moment.py")
assert SPEC and SPEC.loader
DM = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = DM
SPEC.loader.exec_module(DM)

CONTROL = "non-targeting"
PERT_COL = "target_gene"
LAMBDAS = (0.25, 0.5, 1.0)
LIMITS = {0.25: 0.002, 0.5: 0.005, 1.0: 0.02}
GIB = 1024**3
INCUMBENT_PATH = ROOT / (
    "experiments/095-clean-harness/out-h1test/runs/"
    "g0.65_q1.0_a1.08_spread0.25/pred.h5ad"
)
DIRECTION_PATH = ROOT / (
    "experiments/040-union-gamma/out/runs/"
    "g0.45_q1.0_a1.08_spread0.25/pred.h5ad"
)
REAL_PATH = ROOT / "data/harness/h1-fast.h5ad"
OUTPUT_PATH = (
    ROOT
    / "experiments/147-dual-moment-emitter/work/source-feasibility-v2/results.json"
)
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
SOURCE_FILES = (
    HERE / "README.md",
    PROTOCOL_PATH,
    HERE / "dual_moment.py",
    HERE / "test_dual_moment.py",
    HERE / "run_feasibility.py",
    HERE / "run_source_only_guarded.sh",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_seed(base: int, target: str, blend_lambda: float) -> int:
    payload = f"{base}:{target}:{blend_lambda:.2f}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "little")


def read_rows(backed: ad.AnnData, indices: np.ndarray) -> sp.csr_matrix:
    matrix = backed.X[indices]
    result = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)
    return DM.as_integer_csr(result)


def target_indices(labels: np.ndarray) -> dict[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    for label in np.unique(labels):
        if label == CONTROL:
            continue
        result[str(label)] = np.flatnonzero(labels == label)
    return result


def peak_rss_gib() -> float:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if sys.platform == "darwin" else value * 1024 / GIB


def verify_frozen_source_bindings() -> dict[str, str]:
    protocol = json.loads(PROTOCOL_PATH.read_text())
    expected = protocol.get("source_bindings_sha256")
    if not isinstance(expected, dict):
        raise ValueError("protocol has no source_bindings_sha256")
    observed = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in SOURCE_FILES
        if path != PROTOCOL_PATH
    }
    if observed != expected:
        raise ValueError(
            "current source does not match preregistered bindings: "
            f"observed={observed}, expected={expected}"
        )
    return observed


def main() -> int:
    started = time.time()
    if OUTPUT_PATH.exists():
        raise FileExistsError(
            f"immutable canonical output already exists: {OUTPUT_PATH}"
        )
    source_bindings = verify_frozen_source_bindings()
    for path in (INCUMBENT_PATH, DIRECTION_PATH, REAL_PATH):
        if not path.is_file():
            raise FileNotFoundError(path)

    incumbent = ad.read_h5ad(INCUMBENT_PATH, backed="r")
    direction = ad.read_h5ad(DIRECTION_PATH, backed="r")
    real = ad.read_h5ad(REAL_PATH, backed="r")
    try:
        if not np.array_equal(incumbent.var_names, direction.var_names):
            raise ValueError("prediction gene axes differ")
        if not np.array_equal(incumbent.var_names, real.var_names):
            raise ValueError("prediction and development gene axes differ")

        incumbent_labels = incumbent.obs[PERT_COL].astype(str).to_numpy()
        direction_labels = direction.obs[PERT_COL].astype(str).to_numpy()
        real_labels = real.obs[PERT_COL].astype(str).to_numpy()
        incumbent_rows = target_indices(incumbent_labels)
        direction_rows = target_indices(direction_labels)
        targets = [
            str(label)
            for label in dict.fromkeys(real_labels)
            if str(label) != CONTROL
        ]
        if len(targets) != 40:
            raise ValueError(f"expected 40 H1-fast targets, found {len(targets)}")
        missing = sorted(
            set(targets) - set(incumbent_rows) | set(targets) - set(direction_rows)
        )
        if missing:
            raise ValueError(f"prediction artifacts miss targets: {missing}")

        records: list[dict[str, object]] = []
        for target_number, target in enumerate(targets, 1):
            incumbent_index = incumbent_rows[target]
            direction_index = direction_rows[target]
            if incumbent_index.size != 400 or direction_index.size != 400:
                raise ValueError(f"{target}: prediction block is not 400 cells")
            incumbent_block = read_rows(incumbent, incumbent_index)
            direction_block = read_rows(direction, direction_index)

            for blend_lambda in LAMBDAS:
                seed = stable_seed(14700, target, blend_lambda)
                emitted, diagnostics = DM.emit_dual_moment(
                    incumbent_block,
                    direction_block,
                    blend_lambda=blend_lambda,
                    depth_contrast=0.5,
                    quasi_spread=0.25,
                    seed=seed,
                )
                repeated, repeated_diagnostics = DM.emit_dual_moment(
                    incumbent_block,
                    direction_block,
                    blend_lambda=blend_lambda,
                    depth_contrast=0.5,
                    quasi_spread=0.25,
                    seed=seed,
                )
                deterministic = bool(
                    (emitted != repeated).nnz == 0
                    and diagnostics == repeated_diagnostics
                )
                columns_exact = bool(
                    np.array_equal(
                        np.asarray(emitted.sum(axis=0), dtype=np.int64).ravel(),
                        np.asarray(
                            incumbent_block.sum(axis=0), dtype=np.int64
                        ).ravel(),
                    )
                )
                pseudobulk_exact = bool(
                    np.array_equal(
                        DM.scorer_pseudobulk(emitted),
                        DM.scorer_pseudobulk(incumbent_block),
                    )
                )
                row_sums = np.asarray(emitted.sum(axis=1), dtype=np.int64).ravel()
                records.append(
                    {
                        "target": target,
                        "lambda": blend_lambda,
                        "seed": seed,
                        "deterministic": deterministic,
                        "columns_exact": columns_exact,
                        "pseudobulk_exact": pseudobulk_exact,
                        "integer_nonnegative": bool(
                            np.issubdtype(emitted.dtype, np.integer)
                            and emitted.data.min(initial=0) >= 0
                        ),
                        "positive_rows": bool(np.all(row_sums > 0)),
                        "within_count_cap": bool(
                            row_sums.max(initial=0) <= DM.MAX_COUNTS_PER_CELL
                        ),
                        "diagnostics": diagnostics.to_dict(),
                    }
                )
                del emitted, repeated
            print(
                f"[{target_number:02d}/{len(targets)}] {target}",
                flush=True,
            )
            del incumbent_block, direction_block
            gc.collect()

        summaries: dict[str, dict[str, object]] = {}
        for blend_lambda in LAMBDAS:
            subset = [row for row in records if row["lambda"] == blend_lambda]
            errors = np.array(
                [
                    row["diagnostics"]["desired_cpm_relative_l1"]
                    for row in subset
                ],
                dtype=np.float64,
            )
            invariants = all(
                row["deterministic"]
                and row["columns_exact"]
                and row["pseudobulk_exact"]
                and row["integer_nonnegative"]
                and row["positive_rows"]
                and row["within_count_cap"]
                for row in subset
            )
            summaries[str(blend_lambda)] = {
                "targets": len(subset),
                "all_invariants": invariants,
                "desired_cpm_relative_l1_median": float(np.median(errors)),
                "desired_cpm_relative_l1_max": float(np.max(errors)),
                "limit": LIMITS[blend_lambda],
                "passes": bool(invariants and np.median(errors) <= LIMITS[blend_lambda]),
            }

        passed = all(summary["passes"] for summary in summaries.values())
        result = {
            "experiment": "147-dual-moment-emitter",
            "stage": "source_only_feasibility",
            "provenance_utc": "2026-08-31T18:45:00+00:00",
            "development_truth_opened": False,
            "development_control_counts_opened": False,
            "cell_eval2_invoked": False,
            "inputs": {
                "incumbent": str(INCUMBENT_PATH.relative_to(ROOT)),
                "direction": str(DIRECTION_PATH.relative_to(ROOT)),
                "development_metadata": str(REAL_PATH.relative_to(ROOT)),
                "incumbent_sha256": sha256(INCUMBENT_PATH),
                "direction_sha256": sha256(DIRECTION_PATH),
                "development_sha256_not_computed_to_avoid_X_scan": True
            },
            "source_bindings_sha256": source_bindings,
            "protocol_sha256": sha256(PROTOCOL_PATH),
            "parameters": {
                "lambdas": list(LAMBDAS),
                "depth_contrast": 0.5,
                "quasi_spread": 0.25,
                "base_seed": 14700,
            },
            "summaries": summaries,
            "source_only_gate_passed": passed,
            "records": records,
            "peak_rss_gib": peak_rss_gib(),
            "elapsed_seconds": time.time() - started,
        }
        OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps({k: result[k] for k in ("summaries", "source_only_gate_passed", "peak_rss_gib", "elapsed_seconds")}, indent=2))
        return 0 if passed else 2
    finally:
        incumbent.file.close()
        direction.file.close()
        real.file.close()


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Build the fixed lambda-1 dual-moment candidate for all 85 H1-test targets."""

from __future__ import annotations

import gc
import hashlib
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
INCUMBENT_PATH = ROOT / (
    "experiments/095-clean-harness/out-h1test/runs/"
    "g0.65_q1.0_a1.08_spread0.25/pred.h5ad"
)
DIRECTION_PATH = ROOT / (
    "experiments/040-union-gamma/out/runs/"
    "g0.45_q1.0_a1.08_spread0.25/pred.h5ad"
)
REAL_PATH = ROOT / "data/harness/h1-test.h5ad"
FAST_PATH = ROOT / "data/harness/h1-fast.h5ad"
OUTPUT_PATH = HERE / "work/candidate.h5ad"
REPORT_PATH = HERE / "work/build-report.json"
PERT_COL = "target_gene"
CONTROL = "non-targeting"

DM_SPEC = importlib.util.spec_from_file_location(
    "dual_moment_147",
    ROOT / "experiments/147-dual-moment-emitter/dual_moment.py",
)
assert DM_SPEC and DM_SPEC.loader
DM = importlib.util.module_from_spec(DM_SPEC)
sys.modules[DM_SPEC.name] = DM
DM_SPEC.loader.exec_module(DM)

BUILD_SPEC = importlib.util.spec_from_file_location(
    "build_147",
    ROOT / "experiments/147-dual-moment-emitter/build_candidates.py",
)
assert BUILD_SPEC and BUILD_SPEC.loader
BUILD = importlib.util.module_from_spec(BUILD_SPEC)
sys.modules[BUILD_SPEC.name] = BUILD
BUILD_SPEC.loader.exec_module(BUILD)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_seed(base: int, target: str) -> int:
    payload = f"{base}:{target}:1.00".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "little")


def contiguous_runs(labels: np.ndarray) -> list[tuple[str, int, int]]:
    values = np.asarray(labels, dtype=str)
    starts = np.r_[0, np.flatnonzero(values[1:] != values[:-1]) + 1]
    stops = np.r_[starts[1:], values.size]
    runs = [
        (str(values[start]), int(start), int(stop))
        for start, stop in zip(starts, stops, strict=True)
    ]
    names = [label for label, _, _ in runs]
    if len(names) != len(set(names)):
        raise ValueError("prediction labels are not contiguous")
    return runs


def verify_protocol() -> dict[str, object]:
    protocol = json.loads(PROTOCOL_PATH.read_text())
    if protocol.get("status") != "preregistered_no_truth":
        raise ValueError("experiment is not preregistered")
    builder = protocol.get("source_bindings_sha256", {}).get(
        str(Path(__file__).resolve().relative_to(ROOT))
    )
    if builder != sha256(Path(__file__).resolve()):
        raise ValueError("builder does not match preregistered source binding")
    return protocol


def main() -> int:
    started = time.time()
    if OUTPUT_PATH.exists() or REPORT_PATH.exists():
        raise FileExistsError("immutable candidate or build report already exists")
    protocol = verify_protocol()

    incumbent = ad.read_h5ad(INCUMBENT_PATH, backed="r")
    direction = ad.read_h5ad(DIRECTION_PATH, backed="r")
    real = ad.read_h5ad(REAL_PATH, backed="r")
    fast = ad.read_h5ad(FAST_PATH, backed="r")
    temporary = OUTPUT_PATH.with_name(f".{OUTPUT_PATH.name}.tmp-{os.getpid()}")
    temporary.unlink(missing_ok=True)
    try:
        if incumbent.shape != direction.shape:
            raise ValueError("prediction shapes differ")
        if not incumbent.obs.equals(direction.obs):
            raise ValueError("prediction obs metadata differ")
        if not incumbent.var.equals(direction.var):
            raise ValueError("prediction var metadata differ")
        if not np.array_equal(incumbent.var_names, real.var_names):
            raise ValueError("prediction and H1-test gene axes differ")

        labels = incumbent.obs[PERT_COL].astype(str).to_numpy()
        real_labels = real.obs[PERT_COL].astype(str).to_numpy()
        fast_labels = fast.obs[PERT_COL].astype(str).to_numpy()
        full_targets = [label for label, _, _ in contiguous_runs(labels) if label != CONTROL]
        real_targets = set(real_labels) - {CONTROL}
        if set(full_targets) != real_targets or len(full_targets) != 85:
            raise ValueError("prediction and H1-test target sets differ")
        held_targets = [target for target in full_targets if target not in set(fast_labels)]
        if set(held_targets) != set(protocol["held_targets"]):
            raise ValueError("held-target cohort differs from preregistration")
        held_targets = list(protocol["held_targets"])

        OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
        per_target: dict[str, dict[str, object]] = {}
        with (
            h5py.File(INCUMBENT_PATH, "r") as incumbent_h5,
            h5py.File(temporary, "w") as output_h5,
        ):
            BUILD.clone_non_x_metadata(incumbent_h5, output_h5)
            writer = BUILD.CSRStreamWriter(
                output_h5,
                n_rows=incumbent.n_obs,
                n_cols=incumbent.n_vars,
                data_dtype=np.dtype("float32"),
            )
            for run_number, (label, start, stop) in enumerate(
                contiguous_runs(labels),
                1,
            ):
                incumbent_block = BUILD.read_rows(incumbent, slice(start, stop))
                if label == CONTROL:
                    writer.append(incumbent_block.astype(np.float32))
                    del incumbent_block
                    continue
                direction_block = BUILD.read_rows(direction, slice(start, stop))
                emitted, diagnostics = DM.emit_dual_moment(
                    incumbent_block,
                    direction_block,
                    blend_lambda=1.0,
                    depth_contrast=0.5,
                    quasi_spread=0.25,
                    seed=stable_seed(14700, label),
                )
                source_columns = np.asarray(
                    incumbent_block.sum(axis=0),
                    dtype=np.int64,
                ).ravel()
                output_columns = np.asarray(
                    emitted.sum(axis=0),
                    dtype=np.int64,
                ).ravel()
                if not np.array_equal(source_columns, output_columns):
                    raise RuntimeError(f"{label}: column totals changed")
                if not np.array_equal(
                    DM.scorer_pseudobulk(incumbent_block),
                    DM.scorer_pseudobulk(emitted),
                ):
                    raise RuntimeError(f"{label}: scorer pseudobulk changed")
                writer.append(emitted.astype(np.float32))
                per_target[label] = {
                    "rows": stop - start,
                    "pseudobulk_exact": True,
                    "diagnostics": diagnostics.to_dict(),
                }
                if run_number % 20 == 0:
                    print(
                        f"[build] {run_number}/{len(full_targets) + 1} blocks",
                        flush=True,
                    )
                del incumbent_block, direction_block, emitted
                gc.collect()
            writer.finish()
            output_h5.flush()
            output_nnz = writer.nnz
        os.replace(temporary, OUTPUT_PATH)

        built = ad.read_h5ad(OUTPUT_PATH, backed="r")
        try:
            if built.shape != incumbent.shape:
                raise RuntimeError("completed candidate shape changed")
            if not built.obs.equals(incumbent.obs) or not built.var.equals(incumbent.var):
                raise RuntimeError("completed candidate metadata changed")
        finally:
            built.file.close()

        report = {
            "experiment": "148-dual-moment-held-target",
            "stage": "candidate_build",
            "provenance_utc": "2026-08-31T18:45:00+00:00",
            "protocol_sha256": sha256(PROTOCOL_PATH),
            "builder_sha256": sha256(Path(__file__).resolve()),
            "incumbent_sha256": sha256(INCUMBENT_PATH),
            "direction_sha256": sha256(DIRECTION_PATH),
            "candidate": str(OUTPUT_PATH.relative_to(ROOT)),
            "candidate_sha256": sha256(OUTPUT_PATH),
            "targets": len(per_target),
            "held_targets": held_targets,
            "all_pseudobulk_exact": all(
                record["pseudobulk_exact"] for record in per_target.values()
            ),
            "real_expression_opened": False,
            "fast_expression_opened": False,
            "output_nnz": int(output_nnz),
            "elapsed_seconds": time.time() - started,
            "per_target": per_target,
        }
        REPORT_PATH.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(
            json.dumps(
                {
                    "candidate": report["candidate"],
                    "candidate_sha256": report["candidate_sha256"],
                    "targets": report["targets"],
                    "all_pseudobulk_exact": report["all_pseudobulk_exact"],
                    "real_expression_opened": False,
                    "elapsed_seconds": report["elapsed_seconds"],
                },
                indent=2,
            )
        )
        return 0
    finally:
        temporary.unlink(missing_ok=True)
        incumbent.file.close()
        direction.file.close()
        real.file.close()
        fast.file.close()


if __name__ == "__main__":
    raise SystemExit(main())

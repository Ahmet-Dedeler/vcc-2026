#!/usr/bin/env python3
"""Build the frozen H1-fast incumbent and dual-moment candidate artifacts.

The development reference's perturbation matrix is never read. Its obs/var
metadata define the output layout, and only rows labelled ``non-targeting`` are
copied from ``real.X``. Every perturbation block comes from the two prediction
artifacts.
"""

from __future__ import annotations

import gc
import hashlib
import importlib.util
import json
import os
import resource
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
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
GIB = 1024**3
MAX_STORED_ENTRIES = 4_750_000_000
INCUMBENT_PATH = ROOT / (
    "experiments/095-clean-harness/out-h1test/runs/"
    "g0.65_q1.0_a1.08_spread0.25/pred.h5ad"
)
DIRECTION_PATH = ROOT / (
    "experiments/040-union-gamma/out/runs/"
    "g0.45_q1.0_a1.08_spread0.25/pred.h5ad"
)
REAL_PATH = ROOT / "data/harness/h1-fast.h5ad"
SOURCE_REPORT_PATH = (
    ROOT
    / "experiments/147-dual-moment-emitter/work/source-feasibility-v2/results.json"
)
SOURCE_MEMORY_REPORT_PATH = (
    ROOT
    / "experiments/147-dual-moment-emitter/work/source-feasibility-v2/memory-guard.json"
)
OUTDIR = ROOT / "experiments/147-dual-moment-emitter/work/h1-fast-v2"
DEVELOPMENT_PROTOCOL_PATH = HERE / "development_protocol.json"


def target_indices(labels: np.ndarray) -> dict[str, np.ndarray]:
    result: dict[str, np.ndarray] = {}
    for label in np.unique(labels):
        if label != CONTROL:
            result[str(label)] = np.flatnonzero(labels == label)
    return result


def read_rows(backed: ad.AnnData, selector: slice | np.ndarray) -> sp.csr_matrix:
    matrix = backed.X[selector]
    result = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)
    return DM.as_integer_csr(result)


def stable_seed(base: int, target: str, blend_lambda: float) -> int:
    payload = f"{base}:{target}:{blend_lambda:.2f}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "little")


def peak_rss_gib() -> float:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if sys.platform == "darwin" else value * 1024 / GIB


def clone_non_x_metadata(source: h5py.File, destination: h5py.File) -> None:
    for key, value in source.attrs.items():
        destination.attrs[key] = value
    for key in source:
        if key != "X":
            source.copy(key, destination, name=key)


def write_metadata_template(
    path: Path,
    *,
    target_order: list[str],
    cells_per_target: int,
    control_cells: int,
    var: pd.DataFrame,
) -> None:
    labels: list[str] = []
    indices: list[str] = []
    for target in target_order:
        labels.extend([target] * cells_per_target)
        indices.extend(
            f"pred-{target}-{cell:03d}" for cell in range(cells_per_target)
        )
    labels.extend([CONTROL] * control_cells)
    indices.extend(f"control-{cell:04d}" for cell in range(control_cells))
    obs = pd.DataFrame({PERT_COL: labels}, index=indices)
    empty = sp.csr_matrix((len(labels), len(var)), dtype=np.float32)
    ad.AnnData(X=empty, obs=obs, var=var.copy()).write_h5ad(
        path,
        compression="lzf",
    )


class CSRStreamWriter:
    def __init__(
        self,
        h5: h5py.File,
        *,
        n_rows: int,
        n_cols: int,
        data_dtype: np.dtype,
        chunk_nnz: int = 1 << 20,
    ) -> None:
        group = h5.create_group("X")
        group.attrs["encoding-type"] = "csr_matrix"
        group.attrs["encoding-version"] = "0.1.0"
        group.attrs["shape"] = np.array([n_rows, n_cols], dtype=np.int64)
        options = {
            "shape": (0,),
            "maxshape": (None,),
            "chunks": (chunk_nnz,),
            "compression": "lzf",
        }
        self.data = group.create_dataset("data", dtype=data_dtype, **options)
        self.indices = group.create_dataset("indices", dtype=np.int32, **options)
        self.group = group
        self.indptr = np.zeros(n_rows + 1, dtype=np.int64)
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.row = 0
        self.nnz = 0

    def append(self, matrix: sp.csr_matrix) -> None:
        block = matrix.tocsr(copy=False)
        block.sum_duplicates()
        block.eliminate_zeros()
        block.sort_indices()
        if block.shape[1] != self.n_cols:
            raise ValueError("block gene dimension does not match output")
        if self.row + block.shape[0] > self.n_rows:
            raise ValueError("too many rows appended")
        stop_nnz = self.nnz + block.nnz
        if stop_nnz > MAX_STORED_ENTRIES:
            raise ValueError("VCC stored-entry cap exceeded")
        self.data.resize((stop_nnz,))
        self.indices.resize((stop_nnz,))
        if block.nnz:
            self.data[self.nnz : stop_nnz] = block.data
            self.indices[self.nnz : stop_nnz] = block.indices.astype(
                np.int32,
                copy=False,
            )
        self.indptr[self.row + 1 : self.row + 1 + block.shape[0]] = (
            block.indptr[1:].astype(np.int64, copy=False) + self.nnz
        )
        self.row += block.shape[0]
        self.nnz = stop_nnz

    def finish(self) -> None:
        if self.row != self.n_rows:
            raise ValueError(f"wrote {self.row} rows, expected {self.n_rows}")
        self.group.create_dataset(
            "indptr",
            data=self.indptr,
            dtype=np.int64,
            compression="lzf",
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_source_only_report() -> tuple[dict[str, object], dict[str, object]]:
    payload = json.loads(SOURCE_REPORT_PATH.read_text())
    if payload.get("experiment") != "147-dual-moment-emitter":
        raise ValueError("wrong source-only report")
    if payload.get("development_truth_opened") is not False:
        raise ValueError("source-only report does not certify sealed truth")
    if payload.get("source_only_gate_passed") is not True:
        raise ValueError("source-only gate did not pass")
    memory = json.loads(SOURCE_MEMORY_REPORT_PATH.read_text())
    if memory.get("status") != "complete":
        raise ValueError("source-only memory guard did not complete")
    if float(memory.get("peak_rss_gib", 999.0)) > 6.0:
        raise ValueError("source-only run exceeded its memory limit")
    return payload, memory


def verify_development_protocol() -> dict[str, object]:
    protocol = json.loads(DEVELOPMENT_PROTOCOL_PATH.read_text())
    builder = protocol.get("builder", {})
    expected_path = builder.get("path")
    expected_hash = builder.get("sha256")
    observed_path = str(Path(__file__).resolve().relative_to(ROOT))
    observed_hash = sha256(Path(__file__).resolve())
    if expected_path != observed_path or expected_hash != observed_hash:
        raise ValueError(
            "builder does not match frozen development protocol: "
            f"path={observed_path}, sha256={observed_hash}"
        )
    return protocol


def build_one(
    *,
    incumbent: ad.AnnData,
    direction: ad.AnnData,
    real: ad.AnnData,
    target_order: list[str],
    control_indices: np.ndarray,
    incumbent_rows: dict[str, np.ndarray],
    direction_rows: dict[str, np.ndarray],
    output_path: Path,
    blend_lambda: float | None,
    depth_contrast: float,
    quasi_spread: float,
    base_seed: int,
) -> dict[str, object]:
    if output_path.exists():
        raise FileExistsError(f"immutable output already exists: {output_path}")
    temporary = output_path.with_name(f".{output_path.name}.tmp-{os.getpid()}")
    metadata_template = output_path.with_name(
        f".{output_path.name}.metadata-{os.getpid()}"
    )
    temporary.unlink(missing_ok=True)
    metadata_template.unlink(missing_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_metadata_template(
        metadata_template,
        target_order=target_order,
        cells_per_target=400,
        control_cells=int(control_indices.size),
        var=real.var,
    )

    per_target: dict[str, dict[str, object]] = {}
    real_control_rows_read = 0
    real_perturbation_rows_read = 0
    with (
        h5py.File(metadata_template, "r") as metadata_h5,
        h5py.File(temporary, "w") as output_h5,
    ):
        clone_non_x_metadata(metadata_h5, output_h5)
        writer = CSRStreamWriter(
            output_h5,
            n_rows=real.n_obs,
            n_cols=real.n_vars,
            data_dtype=np.dtype("float32"),
        )
        for target_number, target in enumerate(target_order, 1):
            incumbent_block = read_rows(incumbent, incumbent_rows[target])
            if incumbent_block.shape[0] != 400:
                raise ValueError(f"{target}: incumbent row count mismatch")
            if blend_lambda is None:
                output_block = incumbent_block
                diagnostics: dict[str, object] = {"mode": "exact_incumbent_copy"}
            else:
                direction_block = read_rows(direction, direction_rows[target])
                if direction_block.shape[0] != 400:
                    raise ValueError(f"{target}: direction row count mismatch")
                seed = stable_seed(base_seed, target, blend_lambda)
                output_block, emitted_diagnostics = DM.emit_dual_moment(
                    incumbent_block,
                    direction_block,
                    blend_lambda=blend_lambda,
                    depth_contrast=depth_contrast,
                    quasi_spread=quasi_spread,
                    seed=seed,
                )
                diagnostics = emitted_diagnostics.to_dict()

            incumbent_columns = np.asarray(
                incumbent_block.sum(axis=0),
                dtype=np.int64,
            ).ravel()
            output_columns = np.asarray(
                output_block.sum(axis=0),
                dtype=np.int64,
            ).ravel()
            columns_exact = bool(np.array_equal(incumbent_columns, output_columns))
            pseudobulk_exact = bool(
                np.array_equal(
                    DM.scorer_pseudobulk(incumbent_block),
                    DM.scorer_pseudobulk(output_block),
                )
            )
            if not columns_exact or not pseudobulk_exact:
                raise RuntimeError(f"{target}: pooled profile changed")
            row_totals = np.asarray(
                output_block.sum(axis=1),
                dtype=np.int64,
            ).ravel()
            if (
                np.any(row_totals <= 0)
                or row_totals.max(initial=0) > DM.MAX_COUNTS_PER_CELL
            ):
                raise RuntimeError(f"{target}: invalid output row total")
            writer.append(output_block.astype(np.float32))
            per_target[target] = {
                "rows": 400,
                "columns_exact": columns_exact,
                "pseudobulk_exact": pseudobulk_exact,
                "max_row_total": int(row_totals.max(initial=0)),
                "output_nnz": int(output_block.nnz),
                "diagnostics": diagnostics,
            }
            if target_number % 10 == 0:
                print(
                    f"[{output_path.name}] {target_number}/{len(target_order)} targets",
                    flush=True,
                )
            del incumbent_block, output_block
            if blend_lambda is not None:
                del direction_block
            gc.collect()

        # This is the only access to development-reference X, and the selector
        # was constructed solely from rows whose metadata label is CONTROL.
        control_block = read_rows(real, control_indices)
        real_control_rows_read += int(control_indices.size)
        writer.append(control_block.astype(np.float32))
        del control_block
        writer.finish()
        output_h5.flush()
        output_nnz = writer.nnz
    os.replace(temporary, output_path)
    metadata_template.unlink(missing_ok=True)

    # Metadata and shape can be verified without opening any real perturbation X.
    built = ad.read_h5ad(output_path, backed="r")
    try:
        if built.shape != real.shape:
            raise RuntimeError("completed candidate shape differs from development layout")
        if not built.var.equals(real.var):
            raise RuntimeError("completed candidate gene metadata differs")
        built_labels = built.obs[PERT_COL].astype(str).to_numpy()
        if np.count_nonzero(built_labels == CONTROL) != control_indices.size:
            raise RuntimeError("completed candidate control count differs")
        if any(np.count_nonzero(built_labels == target) != 400 for target in target_order):
            raise RuntimeError("completed candidate perturbation count differs")
    finally:
        built.file.close()

    return {
        "output": str(output_path.relative_to(ROOT)),
        "blend_lambda": blend_lambda,
        "shape": [int(real.n_obs), int(real.n_vars)],
        "targets": len(per_target),
        "real_control_rows_read": real_control_rows_read,
        "real_perturbation_rows_read": real_perturbation_rows_read,
        "output_nnz": int(output_nnz),
        "all_pseudobulk_exact": bool(
            all(record["pseudobulk_exact"] for record in per_target.values())
        ),
        "per_target": per_target,
        "peak_rss_gib": peak_rss_gib(),
    }


def main() -> int:
    started = time.time()
    if OUTDIR.exists():
        raise FileExistsError(f"immutable build directory already exists: {OUTDIR}")
    source_report, source_memory_report = verify_source_only_report()
    development_protocol = verify_development_protocol()

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
        target_order = [
            str(label)
            for label in dict.fromkeys(real_labels)
            if str(label) != CONTROL
        ]
        if len(target_order) != 40:
            raise ValueError("expected exactly 40 development perturbations")
        target_counts = {
            target: int(np.count_nonzero(real_labels == target))
            for target in target_order
        }
        if any(count != 400 for count in target_counts.values()):
            raise ValueError(f"development target counts differ: {target_counts}")
        control_indices = np.flatnonzero(real_labels == CONTROL)
        if control_indices.size != 4_000:
            raise ValueError(
                f"expected 4,000 H1-fast controls, found {control_indices.size}"
            )

        incumbent_rows = target_indices(incumbent_labels)
        direction_rows = target_indices(direction_labels)
        targets = set(target_order)
        missing = sorted(
            targets - set(incumbent_rows) | targets - set(direction_rows)
        )
        if missing:
            raise ValueError(f"prediction artifacts miss targets: {missing}")

        builds = []
        builds.append(
            build_one(
                incumbent=incumbent,
                direction=direction,
                real=real,
                target_order=target_order,
                control_indices=control_indices,
                incumbent_rows=incumbent_rows,
                direction_rows=direction_rows,
                output_path=OUTDIR / "incumbent.h5ad",
                blend_lambda=None,
                depth_contrast=0.5,
                quasi_spread=0.25,
                base_seed=14700,
            )
        )
        for blend_lambda in LAMBDAS:
            builds.append(
                build_one(
                    incumbent=incumbent,
                    direction=direction,
                    real=real,
                    target_order=target_order,
                    control_indices=control_indices,
                    incumbent_rows=incumbent_rows,
                    direction_rows=direction_rows,
                    output_path=OUTDIR / f"lambda-{blend_lambda:.2f}.h5ad",
                    blend_lambda=blend_lambda,
                    depth_contrast=0.5,
                    quasi_spread=0.25,
                    base_seed=14700,
                )
            )
        report = {
            "experiment": "147-dual-moment-emitter",
            "stage": "h1_fast_candidate_build",
            "provenance_utc": "2026-08-31T18:45:00+00:00",
            "source_only_gate_verified": True,
            "development_perturbation_truth_opened": False,
            "development_control_counts_opened": True,
            "cell_eval2_invoked": False,
            "source_report_sha256": sha256(SOURCE_REPORT_PATH),
            "source_memory_report_sha256": sha256(SOURCE_MEMORY_REPORT_PATH),
            "source_report_protocol_sha256": source_report["protocol_sha256"],
            "source_memory_peak_rss_gib": source_memory_report["peak_rss_gib"],
            "development_protocol_sha256": sha256(DEVELOPMENT_PROTOCOL_PATH),
            "builder_sha256": development_protocol["builder"]["sha256"],
            "parameters": {
                "lambdas": list(LAMBDAS),
                "depth_contrast": 0.5,
                "quasi_spread": 0.25,
                "base_seed": 14700,
            },
            "builds": builds,
            "peak_rss_gib": peak_rss_gib(),
            "elapsed_seconds": time.time() - started,
        }
        OUTDIR.mkdir(parents=True, exist_ok=True)
        (OUTDIR / "build-report.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n"
        )
        print(
            json.dumps(
                {
                    "outputs": [build["output"] for build in builds],
                    "all_pseudobulk_exact": all(
                        build["all_pseudobulk_exact"] for build in builds
                    ),
                    "development_perturbation_truth_opened": False,
                    "peak_rss_gib": report["peak_rss_gib"],
                    "elapsed_seconds": report["elapsed_seconds"],
                },
                indent=2,
            )
        )
        return 0
    finally:
        incumbent.file.close()
        direction.file.close()
        real.file.close()


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Stream a full HCT116 prediction through the fixed-margin template emitter.

The input and output have identical observation and variable metadata. Every
perturbation remains in the same 400-row position, but its counts are
reallocated across cells using real HCT116 control-cell templates. Per-gene
column totals and scorer pseudobulk are asserted equal before writing and again
from the completed output file.

No cell-eval2 code is imported or invoked.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import resource
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from template_emitter import (
    as_integer_csr,
    emit_template_fixed_margins,
    integer_margin,
    scorer_pseudobulk,
)

CONTROL = "non-targeting"
PERT_COL = "target_gene"
GIB = 1024**3
MAX_COUNTS_PER_CELL = 1_000_000
MAX_STORED_ENTRIES = 4_750_000_000


@dataclass(frozen=True)
class RowRun:
    label: str
    start: int
    stop: int

    @property
    def size(self) -> int:
        return self.stop - self.start


def peak_rss_bytes() -> int:
    """Return peak RSS in bytes on macOS and Linux."""
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def memory_checkpoint(stage: str, limit_bytes: int) -> None:
    rss = peak_rss_bytes()
    print(f"[rss] {stage}: {rss / GIB:.3f} GiB", flush=True)
    if rss >= limit_bytes:
        raise MemoryError(
            f"peak RSS {rss / GIB:.3f} GiB exceeded internal limit "
            f"{limit_bytes / GIB:.3f} GiB at {stage}"
        )


def read_rows(
    backed: ad.AnnData,
    selector: slice | np.ndarray,
) -> sp.csr_matrix:
    """Read only selected backed rows and return canonical CSR."""
    matrix = backed.X[selector]
    result = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)
    result.sum_duplicates()
    result.eliminate_zeros()
    result.sort_indices()
    return result


def row_runs(labels: np.ndarray) -> list[RowRun]:
    """Return contiguous label runs and reject labels split across runs."""
    values = np.asarray(labels, dtype=str)
    if values.ndim != 1 or values.size == 0:
        raise ValueError("target labels must be a non-empty vector")
    starts = np.r_[0, np.flatnonzero(values[1:] != values[:-1]) + 1]
    stops = np.r_[starts[1:], values.size]
    runs = [
        RowRun(str(values[start]), int(start), int(stop))
        for start, stop in zip(starts, stops, strict=True)
    ]
    seen: set[str] = set()
    duplicates: set[str] = set()
    for run in runs:
        if run.label in seen:
            duplicates.add(run.label)
        seen.add(run.label)
    if duplicates:
        raise ValueError(
            "streaming requires one contiguous run per label; split labels: "
            f"{sorted(duplicates)}"
        )
    return runs


def stratified_positions(values: np.ndarray, n: int) -> np.ndarray:
    """Select approximately equal library-size quantiles without replacement."""
    x = np.asarray(values)
    if n <= 0 or n > x.size:
        raise ValueError(f"cannot select {n} templates from {x.size} controls")
    order = np.argsort(x, kind="stable")
    positions = np.floor((np.arange(n) + 0.5) * x.size / n).astype(np.int64)
    selected = order[positions]
    if np.unique(selected).size != n:
        raise RuntimeError("quantile selection produced duplicate rows")
    return selected


def stream_row_sums(
    backed: ad.AnnData,
    indices: np.ndarray,
    *,
    chunk_rows: int,
    rss_limit: int,
) -> np.ndarray:
    """Compute row sums from sorted backed indices in bounded chunks."""
    ordered = np.asarray(indices, dtype=np.int64)
    if ordered.size == 0 or np.any(np.diff(ordered) < 0):
        raise ValueError("indices must be a non-empty sorted vector")
    result = np.empty(ordered.size, dtype=np.int64)
    for start in range(0, ordered.size, chunk_rows):
        stop = min(ordered.size, start + chunk_rows)
        block = read_rows(backed, ordered[start:stop])
        sums = np.asarray(block.sum(axis=1), dtype=np.float64).ravel()
        if not np.array_equal(sums, np.rint(sums)):
            raise ValueError("real control counts are not integer-valued")
        result[start:stop] = np.rint(sums).astype(np.int64)
        if start % (8 * chunk_rows) == 0:
            memory_checkpoint(f"control libraries {stop}/{ordered.size}", rss_limit)
    if np.any(result <= 0):
        raise ValueError("real control cells must have positive library sizes")
    return result


def stable_target_seed(base_seed: int, target: str) -> int:
    """Derive an order-independent uint32 seed for one target."""
    digest = hashlib.sha256(f"{base_seed}:{target}".encode()).digest()
    return int.from_bytes(digest[:4], "little")


def integer_column_totals(block: sp.csr_matrix) -> np.ndarray:
    """Validate a raw count block and return exact int64 column totals."""
    integer = as_integer_csr(block)
    columns = np.asarray(integer.sum(axis=0), dtype=np.int64).ravel()
    if int(columns.sum()) <= 0:
        raise ValueError("prediction block has zero total counts")
    row_sums = np.asarray(integer.sum(axis=1), dtype=np.int64).ravel()
    if np.any(row_sums <= 0):
        raise ValueError("prediction block contains an empty cell")
    if np.any(row_sums > MAX_COUNTS_PER_CELL):
        raise ValueError("prediction block exceeds 1,000,000 counts per cell")
    return columns


def columns_sha256(columns: np.ndarray) -> str:
    """Stable digest of a target's exact pseudobulk sufficient statistic."""
    canonical = np.asarray(columns, dtype="<i8")
    return hashlib.sha256(canonical.tobytes()).hexdigest()


def clone_non_x_metadata(source: h5py.File, destination: h5py.File) -> None:
    """Copy the complete AnnData container except for X."""
    for key, value in source.attrs.items():
        destination.attrs[key] = value
    for key in source:
        if key != "X":
            source.copy(key, destination, name=key)


class CSRStreamWriter:
    """Append CSR blocks into a growable H5AD X group."""

    def __init__(
        self,
        h5: h5py.File,
        *,
        n_rows: int,
        n_cols: int,
        data_dtype: np.dtype,
        chunk_nnz: int,
        compression: str | None,
    ) -> None:
        group = h5.create_group("X")
        group.attrs["encoding-type"] = "csr_matrix"
        group.attrs["encoding-version"] = "0.1.0"
        group.attrs["shape"] = np.array([n_rows, n_cols], dtype=np.int64)
        dataset_args: dict[str, object] = {
            "shape": (0,),
            "maxshape": (None,),
            "chunks": (chunk_nnz,),
        }
        if compression is not None:
            dataset_args["compression"] = compression
        self.data = group.create_dataset("data", dtype=data_dtype, **dataset_args)
        self.indices = group.create_dataset("indices", dtype=np.int32, **dataset_args)
        self.group = group
        self.indptr = np.zeros(n_rows + 1, dtype=np.int64)
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.row = 0
        self.nnz = 0
        self.compression = compression

    def append(self, matrix: sp.csr_matrix) -> None:
        block = matrix.tocsr(copy=False)
        block.sum_duplicates()
        block.eliminate_zeros()
        block.sort_indices()
        if block.shape[1] != self.n_cols:
            raise ValueError(
                f"block has {block.shape[1]} columns, expected {self.n_cols}"
            )
        if self.row + block.shape[0] > self.n_rows:
            raise ValueError("appended rows exceed declared output shape")
        if block.nnz and (
            np.any(block.indices < 0) or np.any(block.indices >= self.n_cols)
        ):
            raise ValueError("CSR block contains an out-of-range column index")

        stop_nnz = self.nnz + block.nnz
        if stop_nnz > MAX_STORED_ENTRIES:
            raise ValueError(
                f"stored-entry cap exceeded: {stop_nnz:,} > {MAX_STORED_ENTRIES:,}"
            )
        self.data.resize((stop_nnz,))
        self.indices.resize((stop_nnz,))
        if block.nnz:
            self.data[self.nnz : stop_nnz] = block.data
            self.indices[self.nnz : stop_nnz] = block.indices.astype(
                np.int32, copy=False
            )
        local = block.indptr[1:].astype(np.int64, copy=False) + self.nnz
        self.indptr[self.row + 1 : self.row + 1 + block.shape[0]] = local
        self.row += block.shape[0]
        self.nnz = stop_nnz

    def finish(self) -> None:
        if self.row != self.n_rows:
            raise ValueError(f"wrote {self.row} rows, expected {self.n_rows}")
        kwargs = {"compression": self.compression} if self.compression else {}
        self.group.create_dataset("indptr", data=self.indptr, dtype=np.int64, **kwargs)


def validate_inputs(
    prediction: ad.AnnData,
    real: ad.AnnData,
    *,
    pert_col: str,
    control: str,
    cells_per_perturbation: int,
) -> tuple[np.ndarray, np.ndarray, list[RowRun], np.ndarray]:
    """Validate the exact full-HCT streaming contract."""
    if pert_col not in prediction.obs or pert_col not in real.obs:
        raise ValueError(f"both files must contain obs[{pert_col!r}]")
    if prediction.n_vars != real.n_vars or not np.array_equal(
        np.asarray(prediction.var_names, dtype=str),
        np.asarray(real.var_names, dtype=str),
    ):
        raise ValueError(
            "prediction and real HCT116 files must have identical gene axes"
        )

    pred_labels = prediction.obs[pert_col].astype(str).to_numpy()
    real_labels = real.obs[pert_col].astype(str).to_numpy()
    runs = row_runs(pred_labels)
    pred_targets = set(pred_labels) - {control}
    real_targets = set(real_labels) - {control}
    if pred_targets != real_targets:
        missing = sorted(real_targets - pred_targets)
        extra = sorted(pred_targets - real_targets)
        raise ValueError(
            f"prediction target set differs from HCT116; missing={missing}, "
            f"extra={extra}"
        )
    if control not in set(pred_labels):
        raise ValueError("full scorer prediction must include control rows")
    for run in runs:
        if run.label != control and run.size != cells_per_perturbation:
            raise ValueError(
                f"{run.label} has {run.size} rows, expected {cells_per_perturbation}"
            )
    control_indices = np.flatnonzero(real_labels == control)
    if control_indices.size == 0:
        raise ValueError("real HCT116 file has no control cells")
    return pred_labels, real_labels, runs, control_indices


def verify_completed_output(
    prediction_path: Path,
    output_path: Path,
    *,
    runs: list[RowRun],
    pert_col: str,
    control: str,
    chunk_rows: int,
    rss_limit: int,
) -> dict[str, object]:
    """Re-read the output and verify metadata, controls, and every pseudobulk."""
    source = ad.read_h5ad(prediction_path, backed="r")
    output = ad.read_h5ad(output_path, backed="r")
    try:
        if source.shape != output.shape:
            raise RuntimeError(
                f"output shape {output.shape} differs from input {source.shape}"
            )
        if not source.obs.equals(output.obs):
            raise RuntimeError("output obs does not exactly match input obs")
        if not source.var.equals(output.var):
            raise RuntimeError("output var does not exactly match input var")
        if not np.array_equal(
            source.obs[pert_col].astype(str).to_numpy(),
            output.obs[pert_col].astype(str).to_numpy(),
        ):
            raise RuntimeError("output perturbation labels changed")

        target_checks: dict[str, dict[str, object]] = {}
        controls_equal = True
        for run_index, run in enumerate(runs):
            if run.label == control:
                for start in range(run.start, run.stop, chunk_rows):
                    stop = min(run.stop, start + chunk_rows)
                    source_block = read_rows(source, slice(start, stop))
                    output_block = read_rows(output, slice(start, stop))
                    if (source_block != output_block).nnz:
                        controls_equal = False
                        raise RuntimeError(f"control rows changed in [{start}, {stop})")
                continue

            source_block = read_rows(source, slice(run.start, run.stop))
            output_block = read_rows(output, slice(run.start, run.stop))
            source_columns = integer_column_totals(source_block)
            output_columns = integer_column_totals(output_block)
            column_equal = bool(np.array_equal(source_columns, output_columns))
            source_bulk = scorer_pseudobulk(source_block)
            output_bulk = scorer_pseudobulk(output_block)
            bulk_equal = bool(np.array_equal(source_bulk, output_bulk))
            max_column_difference = int(
                np.max(
                    np.abs(source_columns - output_columns),
                    initial=0,
                )
            )
            max_bulk_difference = float(
                np.max(np.abs(source_bulk - output_bulk), initial=0.0)
            )
            if not column_equal or not bulk_equal:
                raise RuntimeError(f"{run.label}: completed output changed pseudobulk")
            target_checks[run.label] = {
                "column_totals_exact": column_equal,
                "bulk_lognorm_array_equal": bulk_equal,
                "max_abs_column_total_difference": max_column_difference,
                "max_abs_bulk_lognorm_difference": max_bulk_difference,
                "columns_sha256": columns_sha256(source_columns),
            }
            if run_index % 16 == 0:
                memory_checkpoint(
                    f"post-verified run {run_index + 1}/{len(runs)}",
                    rss_limit,
                )
            del source_block, output_block
            gc.collect()

        return {
            "obs_exact": True,
            "var_exact": True,
            "controls_exact": controls_equal,
            "targets_checked": len(target_checks),
            "all_pseudobulk_exact": bool(
                all(
                    check["column_totals_exact"] and check["bulk_lognorm_array_equal"]
                    for check in target_checks.values()
                )
            ),
            "per_target": target_checks,
        }
    finally:
        source.file.close()
        output.file.close()


def transform_file(
    prediction_path: str | Path,
    real_path: str | Path,
    output_path: str | Path,
    *,
    report_path: str | Path | None = None,
    pert_col: str = PERT_COL,
    control: str = CONTROL,
    cells_per_perturbation: int = 400,
    template_cells: int = 400,
    chunk_rows: int = 256,
    chunk_nnz: int = 1 << 20,
    compression: str | None = "lzf",
    seed: int = 5300,
    max_rss_gib: float = 3.8,
    dry_run: bool = False,
    post_verify: bool = True,
    overwrite: bool = False,
) -> dict[str, object]:
    """Transform one complete HCT116 prediction using bounded memory."""
    started = time.time()
    prediction_path = Path(prediction_path)
    real_path = Path(real_path)
    output_path = Path(output_path)
    if report_path is None:
        report_path = output_path.with_suffix(".template-report.json")
    report_path = Path(report_path)
    rss_limit = int(max_rss_gib * GIB)
    memory_checkpoint("start", rss_limit)

    if not prediction_path.is_file():
        raise FileNotFoundError(prediction_path)
    if not real_path.is_file():
        raise FileNotFoundError(real_path)
    if output_path.exists() and not overwrite and not dry_run:
        raise FileExistsError(
            f"{output_path} exists; pass --overwrite to replace it atomically"
        )

    prediction = ad.read_h5ad(prediction_path, backed="r")
    real = ad.read_h5ad(real_path, backed="r")
    temporary_path = output_path.with_name(f".{output_path.name}.tmp-{os.getpid()}")
    try:
        pred_labels, _, runs, real_control_indices = validate_inputs(
            prediction,
            real,
            pert_col=pert_col,
            control=control,
            cells_per_perturbation=cells_per_perturbation,
        )
        prediction_shape = (int(prediction.n_obs), int(prediction.n_vars))
        control_libraries = stream_row_sums(
            real,
            real_control_indices,
            chunk_rows=chunk_rows,
            rss_limit=rss_limit,
        )
        selected_local = stratified_positions(control_libraries, template_cells)
        selected_global = np.sort(real_control_indices[selected_local])
        templates = as_integer_csr(read_rows(real, selected_global))
        template_libraries = np.asarray(templates.sum(axis=1), dtype=np.int64).ravel()
        memory_checkpoint("templates loaded", rss_limit)

        targets = [run for run in runs if run.label != control]
        if dry_run:
            first = targets[0]
            source_block = as_integer_csr(
                read_rows(prediction, slice(first.start, first.stop))
            )
            columns = integer_column_totals(source_block)
            target_seed = stable_target_seed(seed, first.label)
            rows = integer_margin(
                template_libraries,
                int(columns.sum()),
                seed=target_seed + 1,
            )
            emitted, diagnostics = emit_template_fixed_margins(
                templates,
                column_totals=columns,
                row_margins=rows,
                seed=target_seed,
            )
            emitted_columns = integer_column_totals(emitted)
            if not np.array_equal(columns, emitted_columns):
                raise RuntimeError("dry-run emission changed column totals")
            if not np.array_equal(
                scorer_pseudobulk(source_block),
                scorer_pseudobulk(emitted),
            ):
                raise RuntimeError("dry-run emission changed scorer pseudobulk")
            result = {
                "mode": "dry-run",
                "prediction": str(prediction_path),
                "real": str(real_path),
                "shape": list(prediction_shape),
                "targets": len(targets),
                "control_rows": int(np.count_nonzero(pred_labels == control)),
                "template_cells": template_cells,
                "exercised_target": first.label,
                "exact_pseudobulk": True,
                "emitter_diagnostics": diagnostics.to_dict(),
                "peak_rss_gib": peak_rss_bytes() / GIB,
                "elapsed_seconds": time.time() - started,
            }
            print(json.dumps(result, indent=2, sort_keys=True), flush=True)
            return result

        output_path.parent.mkdir(parents=True, exist_ok=True)
        if temporary_path.exists():
            temporary_path.unlink()
        per_target: dict[str, dict[str, object]] = {}
        copied_control_rows = 0
        try:
            with (
                h5py.File(prediction_path, "r") as source_h5,
                h5py.File(temporary_path, "w") as output_h5,
            ):
                clone_non_x_metadata(source_h5, output_h5)
                source_data_dtype = source_h5["X"]["data"].dtype
                writer = CSRStreamWriter(
                    output_h5,
                    n_rows=prediction.n_obs,
                    n_cols=prediction.n_vars,
                    data_dtype=source_data_dtype,
                    chunk_nnz=chunk_nnz,
                    compression=compression,
                )
                for run_index, run in enumerate(runs):
                    if run.label == control:
                        for start in range(run.start, run.stop, chunk_rows):
                            stop = min(run.stop, start + chunk_rows)
                            block = read_rows(prediction, slice(start, stop))
                            integer_column_totals(block)
                            writer.append(block)
                            copied_control_rows += stop - start
                        memory_checkpoint("controls copied", rss_limit)
                        continue

                    source_block = as_integer_csr(
                        read_rows(prediction, slice(run.start, run.stop))
                    )
                    columns = integer_column_totals(source_block)
                    target_seed = stable_target_seed(seed, run.label)
                    rows = integer_margin(
                        template_libraries,
                        int(columns.sum()),
                        seed=target_seed + 1,
                    )
                    emitted, diagnostics = emit_template_fixed_margins(
                        templates,
                        column_totals=columns,
                        row_margins=rows,
                        seed=target_seed,
                    )
                    emitted_columns = integer_column_totals(emitted)
                    source_bulk = scorer_pseudobulk(source_block)
                    emitted_bulk = scorer_pseudobulk(emitted)
                    column_equal = bool(np.array_equal(columns, emitted_columns))
                    bulk_equal = bool(np.array_equal(source_bulk, emitted_bulk))
                    if not column_equal or not bulk_equal:
                        raise RuntimeError(
                            f"{run.label}: in-memory exact-pseudobulk assertion failed"
                        )
                    writer.append(emitted)
                    per_target[run.label] = {
                        "rows": run.size,
                        "input_nnz": int(source_block.nnz),
                        "output_nnz": int(emitted.nnz),
                        "column_totals_exact": column_equal,
                        "bulk_lognorm_array_equal": bulk_equal,
                        "max_abs_column_total_difference": int(
                            np.max(
                                np.abs(columns - emitted_columns),
                                initial=0,
                            )
                        ),
                        "max_abs_bulk_lognorm_difference": float(
                            np.max(
                                np.abs(source_bulk - emitted_bulk),
                                initial=0.0,
                            )
                        ),
                        "columns_sha256": columns_sha256(columns),
                        "emitter_diagnostics": diagnostics.to_dict(),
                    }
                    if run_index % 8 == 0:
                        memory_checkpoint(
                            f"transformed run {run_index + 1}/{len(runs)}",
                            rss_limit,
                        )
                    del source_block, emitted
                    gc.collect()

                writer.finish()
                output_h5.flush()
                output_nnz = writer.nnz
            os.replace(temporary_path, output_path)
        except BaseException:
            if temporary_path.exists():
                temporary_path.unlink()
            raise

        memory_checkpoint("streaming transform complete", rss_limit)
    finally:
        prediction.file.close()
        real.file.close()

    verification_started = time.time()
    verification = (
        verify_completed_output(
            prediction_path,
            output_path,
            runs=runs,
            pert_col=pert_col,
            control=control,
            chunk_rows=chunk_rows,
            rss_limit=rss_limit,
        )
        if post_verify
        else {"skipped": True}
    )
    report = {
        "date": "2026-08-22",
        "experiment": "053-template-emitter-full-hct",
        "prediction": str(prediction_path),
        "real": str(real_path),
        "output": str(output_path),
        "shape": list(prediction_shape),
        "perturbations_transformed": len(per_target),
        "control_rows_copied_exactly": copied_control_rows,
        "template_cells": template_cells,
        "seed": seed,
        "output_nnz": int(output_nnz),
        "all_in_memory_pseudobulk_exact": bool(
            all(
                result["column_totals_exact"] and result["bulk_lognorm_array_equal"]
                for result in per_target.values()
            )
        ),
        "per_target": per_target,
        "post_write_verification": verification,
        "post_verify_seconds": time.time() - verification_started,
        "peak_rss_gib": peak_rss_bytes() / GIB,
        "elapsed_seconds": time.time() - started,
        "notes": [
            "All non-X H5AD objects and root attributes are copied from the input.",
            "Control X rows are copied unchanged; only perturbation rows are re-emitted.",
            "The output is written to a temporary sibling and atomically renamed.",
            "No cell-eval2 command or challenge submission is invoked.",
        ],
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    memory_checkpoint("done", rss_limit)
    print(
        f"wrote {output_path} with {len(per_target)} transformed targets; "
        f"peak RSS={report['peak_rss_gib']:.3f} GiB",
        flush=True,
    )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction", required=True)
    parser.add_argument("--real", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--report")
    parser.add_argument("--pert-col", default=PERT_COL)
    parser.add_argument("--control", default=CONTROL)
    parser.add_argument("--cells-per-perturbation", type=int, default=400)
    parser.add_argument("--template-cells", type=int, default=400)
    parser.add_argument("--chunk-rows", type=int, default=256)
    parser.add_argument("--chunk-nnz", type=int, default=1 << 20)
    parser.add_argument(
        "--compression",
        choices=("lzf", "gzip", "none"),
        default="lzf",
    )
    parser.add_argument("--seed", type=int, default=5300)
    parser.add_argument("--max-rss-gib", type=float, default=3.8)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-post-verify", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    transform_file(
        args.prediction,
        args.real,
        args.output,
        report_path=args.report,
        pert_col=args.pert_col,
        control=args.control,
        cells_per_perturbation=args.cells_per_perturbation,
        template_cells=args.template_cells,
        chunk_rows=args.chunk_rows,
        chunk_nnz=args.chunk_nnz,
        compression=None if args.compression == "none" else args.compression,
        seed=args.seed,
        max_rss_gib=args.max_rss_gib,
        dry_run=args.dry_run,
        post_verify=not args.skip_post_verify,
        overwrite=args.overwrite,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

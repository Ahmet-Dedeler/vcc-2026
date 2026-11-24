"""Extract VCC-panel donor tables from the public primary-CD4 DE atlas.

This script is intended to run on an ephemeral cloud worker beside
``GWCD4i.DE_stats.h5ad``.  It deliberately does not use
``anndata.read_h5ad(..., backed="r")``: current AnnData versions still load
``.layers`` eagerly in backed mode, and the six dense DE layers in this file
are much larger than the compact slice we need.  Instead, observation and
variable dataframes are decoded with AnnData's public ``read_elem`` API while
matrix layers remain HDF5-backed and are copied in small row blocks.

For each of Rest, Stim8hr, and Stim48hr the output directory contains:

* ``delta.npy``: finite float32 natural-log fold changes, compatible with the
  effect scale used by this repository (source ``log_fc`` is log2);
* ``b_ctrl.npy``: a compatibility-only basal proxy derived from the median
  DESeq2 ``baseMean`` profile, explicitly labelled as a proxy in ``meta.json``;
* ``perts.txt`` and ``genes.txt``: row and column labels;
* source statistic matrices and a finite-value mask; and
* ``obs.csv`` plus ``meta.json`` with provenance and QC.

No single-cell data or challenge submission is read or written.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import resource
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any

# Avoid accidental thread storms on a cloud VM or laptop.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem


SOURCE_URL = (
    "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com/"
    "marson2025_data/GWCD4i.DE_stats.h5ad"
)
CONDITIONS = ("Rest", "Stim8hr", "Stim48hr")
TARGET_SUM = 5e4
LN2 = math.log(2.0)

OBS_COLUMNS = (
    "target_contrast_gene_name",
    "culture_condition",
    "target_contrast",
    "chunk",
    "n_cells_target",
    "n_up_genes",
    "n_down_genes",
    "n_total_de_genes",
    "ontarget_effect_size",
    "ontarget_significant",
    "target_baseMean",
    "neighboring_gene_KD",
    "distal_offtarget_flag",
    "low_target_gex",
    "n_guides",
    "single_guide_estimate",
    "n_downstream",
    "guide_correlation_signif",
    "guide_correlation_all",
    "guide_n_signif_ontarget",
    "donor_correlation_all_mean",
    "donor_correlation_all_min",
    "donor_correlation_hits_mean",
    "donor_correlation_hits_min",
)

LAYER_OUTPUTS = {
    "log_fc": "log_fc.npy",
    "zscore": "zscore.npy",
    "lfcSE": "lfc_se.npy",
    "adj_p_value": "adj_p_value.npy",
    "p_value": "p_value.npy",
    "baseMean": "base_mean.npy",
}


def peak_rss_gib() -> float:
    """Return peak process RSS in GiB on macOS or Linux."""
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # macOS reports bytes; Linux reports KiB.
    return rss / 1024**3 if os.uname().sysname == "Darwin" else rss / 1024**2


def enforce_memory_limit(limit_gib: float) -> None:
    if limit_gib > 0 and peak_rss_gib() > limit_gib:
        raise MemoryError(
            f"peak RSS {peak_rss_gib():.2f} GiB exceeded limit {limit_gib:.2f} GiB"
        )


def sha256_file(path: Path, block_bytes: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def clean_string(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    text = str(value).strip()
    return "" if text.casefold() in {"nan", "none", "<na>"} else text


def clean_ensembl(value: Any) -> str:
    return clean_string(value).split(".", 1)[0]


def bool_value(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return clean_string(value).casefold() in {"1", "true", "t", "yes", "y"}


def numeric_value(value: Any, default: float = float("-inf")) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return default
    return number if np.isfinite(number) else default


def matrix_shape(elem: h5py.Dataset | h5py.Group) -> tuple[int, int]:
    if isinstance(elem, h5py.Dataset):
        return tuple(map(int, elem.shape))
    shape = elem.attrs.get("shape")
    if shape is None:
        raise ValueError(f"{elem.name}: sparse matrix has no shape attribute")
    return tuple(map(int, shape))


def encoding_type(elem: h5py.Dataset | h5py.Group) -> str:
    value = elem.attrs.get("encoding-type", "")
    return clean_string(value)


def read_matrix_subset_into(
    elem: h5py.Dataset | h5py.Group,
    rows: np.ndarray,
    columns: np.ndarray,
    output: np.ndarray,
    *,
    row_block: int,
    max_rss_gib: float,
) -> None:
    """Copy a 2-D HDF5 matrix subset into ``output`` without loading the matrix."""
    rows = np.asarray(rows, dtype=np.int64)
    columns = np.asarray(columns, dtype=np.int64)
    if output.shape != (len(rows), len(columns)):
        raise ValueError(
            f"output shape {output.shape} != requested {(len(rows), len(columns))}"
        )
    if len(np.unique(rows)) != len(rows):
        raise ValueError("row selection contains duplicates")
    if len(np.unique(columns)) != len(columns):
        raise ValueError("column selection contains duplicates")

    n_rows, n_cols = matrix_shape(elem)
    if rows.size and (rows.min() < 0 or rows.max() >= n_rows):
        raise IndexError("row selection is outside matrix bounds")
    if columns.size and (columns.min() < 0 or columns.max() >= n_cols):
        raise IndexError("column selection is outside matrix bounds")

    kind = encoding_type(elem)
    if isinstance(elem, h5py.Dataset) or kind in {"", "array"}:
        # h5py requires monotonically increasing fancy indices.  Preserve the
        # requested panel order by writing each sorted read back to its position.
        order = np.argsort(rows)
        for lo in range(0, len(rows), row_block):
            positions = order[lo : lo + row_block]
            source_rows = rows[positions]
            block = np.asarray(elem[source_rows, :], dtype=np.float32)
            output[positions] = block[:, columns]
            enforce_memory_limit(max_rss_gib)
        return

    if kind == "csr_matrix":
        indptr = np.asarray(elem["indptr"], dtype=np.int64)
        source_to_dest = np.full(n_cols, -1, dtype=np.int64)
        source_to_dest[columns] = np.arange(len(columns), dtype=np.int64)
        for out_row, source_row in enumerate(rows):
            lo, hi = int(indptr[source_row]), int(indptr[source_row + 1])
            source_cols = np.asarray(elem["indices"][lo:hi], dtype=np.int64)
            destination = source_to_dest[source_cols]
            keep = destination >= 0
            if keep.any():
                output[out_row, destination[keep]] = np.asarray(
                    elem["data"][lo:hi], dtype=np.float32
                )[keep]
            if out_row % row_block == 0:
                enforce_memory_limit(max_rss_gib)
        return

    if kind == "csc_matrix":
        indptr = np.asarray(elem["indptr"], dtype=np.int64)
        source_to_dest = np.full(n_rows, -1, dtype=np.int64)
        source_to_dest[rows] = np.arange(len(rows), dtype=np.int64)
        for out_col, source_col in enumerate(columns):
            lo, hi = int(indptr[source_col]), int(indptr[source_col + 1])
            source_rows = np.asarray(elem["indices"][lo:hi], dtype=np.int64)
            destination = source_to_dest[source_rows]
            keep = destination >= 0
            if keep.any():
                output[destination[keep], out_col] = np.asarray(
                    elem["data"][lo:hi], dtype=np.float32
                )[keep]
            if out_col % max(row_block, 1) == 0:
                enforce_memory_limit(max_rss_gib)
        return

    raise TypeError(f"{elem.name}: unsupported H5AD matrix encoding {kind!r}")


def read_matrix_subset(
    elem: h5py.Dataset | h5py.Group,
    rows: np.ndarray,
    columns: np.ndarray,
    *,
    row_block: int,
    max_rss_gib: float,
) -> np.ndarray:
    out = np.zeros((len(rows), len(columns)), dtype=np.float32)
    read_matrix_subset_into(
        elem,
        rows,
        columns,
        out,
        row_block=row_block,
        max_rss_gib=max_rss_gib,
    )
    return out


def resolve_target_rows(
    obs: pd.DataFrame,
    var: pd.DataFrame,
    panel: list[str],
) -> tuple[dict[str, dict[str, int]], list[dict[str, Any]]]:
    """Map source perturbation rows to VCC symbols with exact/Ensembl-safe rules."""
    required = {"target_contrast_gene_name", "culture_condition"}
    missing = required - set(obs.columns)
    if missing:
        raise ValueError(f"source .obs is missing required columns: {sorted(missing)}")

    panel_set = set(panel)
    casefold_panel: dict[str, list[str]] = defaultdict(list)
    for symbol in panel:
        casefold_panel[symbol.casefold()].append(symbol)

    var_symbols = (
        [clean_string(x) for x in var["gene_name"]]
        if "gene_name" in var.columns
        else [clean_string(x) for x in var.index]
    )
    var_ids = (
        [clean_ensembl(x) for x in var["gene_ids"]]
        if "gene_ids" in var.columns
        else [clean_ensembl(x) for x in var.index]
    )
    id_to_panel: dict[str, set[str]] = defaultdict(set)
    for symbol, gene_id in zip(var_symbols, var_ids, strict=True):
        if symbol in panel_set and gene_id:
            id_to_panel[gene_id].add(symbol)

    mapped: dict[tuple[str, str], list[tuple[int, str]]] = defaultdict(list)
    audit: list[dict[str, Any]] = []
    for row_number, (_, row) in enumerate(obs.iterrows()):
        condition = clean_string(row.get("culture_condition"))
        if condition not in CONDITIONS:
            continue
        raw_symbol = clean_string(row.get("target_contrast_gene_name"))
        raw_id = clean_ensembl(row.get("target_contrast"))
        resolved = ""
        method = ""
        if raw_symbol in panel_set:
            resolved, method = raw_symbol, "exact_symbol"
        else:
            folded = casefold_panel.get(raw_symbol.casefold(), [])
            if len(folded) == 1:
                resolved, method = folded[0], "casefold_symbol"
            elif raw_id and len(id_to_panel.get(raw_id, ())) == 1:
                resolved = next(iter(id_to_panel[raw_id]))
                method = "ensembl_id"
        if resolved:
            mapped[(condition, resolved)].append((row_number, method))

    selected: dict[str, dict[str, int]] = {condition: {} for condition in CONDITIONS}
    for (condition, target), candidates in mapped.items():

        def rank(candidate: tuple[int, str]) -> tuple[Any, ...]:
            i, _ = candidate
            row = obs.iloc[i]
            return (
                bool_value(row.get("ontarget_significant")),
                not bool_value(row.get("low_target_gex")),
                not bool_value(row.get("distal_offtarget_flag")),
                not bool_value(row.get("neighboring_gene_KD")),
                numeric_value(row.get("n_guides")),
                numeric_value(row.get("n_cells_target")),
                -i,
            )

        chosen = max(candidates, key=rank)
        selected[condition][target] = chosen[0]
        for row_number, method in candidates:
            audit.append(
                {
                    "condition": condition,
                    "target_gene": target,
                    "source_row": row_number,
                    "source_obs_name": clean_string(obs.index[row_number]),
                    "mapping_method": method,
                    "selected": row_number == chosen[0],
                    "duplicate_candidates": len(candidates),
                }
            )
    return selected, audit


def choose_gene_columns(
    var: pd.DataFrame,
    official_genes: list[str],
    base_mean: h5py.Dataset | h5py.Group,
    selected_rows: np.ndarray,
    *,
    row_block: int,
    max_rss_gib: float,
) -> tuple[list[str], np.ndarray, list[dict[str, Any]]]:
    """Map measured genes to unique official symbols.

    If several source variables map to one symbol, retain the feature with the
    largest median ``baseMean`` over selected panel rows.  This is deterministic
    and avoids averaging fold changes from distinct Ensembl features.
    """
    symbols = (
        np.asarray([clean_string(x) for x in var["gene_name"]], dtype=object)
        if "gene_name" in var.columns
        else np.asarray([clean_string(x) for x in var.index], dtype=object)
    )
    gene_ids = (
        np.asarray([clean_ensembl(x) for x in var["gene_ids"]], dtype=object)
        if "gene_ids" in var.columns
        else np.asarray([clean_ensembl(x) for x in var.index], dtype=object)
    )
    groups: dict[str, list[int]] = defaultdict(list)
    official_set = set(official_genes)
    for source_col, symbol in enumerate(symbols):
        if symbol in official_set:
            groups[str(symbol)].append(source_col)

    duplicate_source_cols = sorted(
        source_col
        for candidates in groups.values()
        if len(candidates) > 1
        for source_col in candidates
    )
    duplicate_medians: dict[int, float] = {}
    if duplicate_source_cols:
        values = read_matrix_subset(
            base_mean,
            selected_rows,
            np.asarray(duplicate_source_cols, dtype=np.int64),
            row_block=row_block,
            max_rss_gib=max_rss_gib,
        )
        with np.errstate(all="ignore"):
            medians = np.nanmedian(values, axis=0)
        duplicate_medians = {
            source_col: numeric_value(median)
            for source_col, median in zip(duplicate_source_cols, medians, strict=True)
        }

    output_genes: list[str] = []
    output_columns: list[int] = []
    audit: list[dict[str, Any]] = []
    for symbol in official_genes:
        candidates = groups.get(symbol, [])
        if not candidates:
            continue
        if len(candidates) == 1:
            chosen = candidates[0]
            method = "unique_symbol"
        else:
            chosen = max(
                candidates,
                key=lambda col: (duplicate_medians.get(col, float("-inf")), -col),
            )
            method = "max_median_baseMean"
        output_genes.append(symbol)
        output_columns.append(chosen)
        audit.append(
            {
                "gene_name": symbol,
                "source_var_index": chosen,
                "source_var_name": clean_string(var.index[chosen]),
                "source_gene_id": clean_ensembl(gene_ids[chosen]),
                "duplicate_count": len(candidates),
                "selection_method": method,
                "candidate_source_indices": ";".join(map(str, candidates)),
            }
        )
    return output_genes, np.asarray(output_columns, dtype=np.int64), audit


def condition_slug(condition: str) -> str:
    return {
        "Rest": "cd4_rest",
        "Stim8hr": "cd4_stim8hr",
        "Stim48hr": "cd4_stim48hr",
    }[condition]


def write_condition(
    source: h5py.File,
    obs: pd.DataFrame,
    condition: str,
    target_rows: dict[str, int],
    panel: list[str],
    genes: list[str],
    gene_columns: np.ndarray,
    out_dir: Path,
    *,
    row_block: int,
    max_rss_gib: float,
    min_targets: int,
    min_genes: int,
    min_finite_fraction: float,
) -> dict[str, Any]:
    ordered_perts = [target for target in panel if target in target_rows]
    rows = np.asarray([target_rows[target] for target in ordered_perts], dtype=np.int64)
    if len(ordered_perts) < min_targets:
        raise ValueError(
            f"{condition}: only {len(ordered_perts)} panel targets; require {min_targets}"
        )
    if len(genes) < min_genes:
        raise ValueError(
            f"{condition}: only {len(genes)} mapped genes; require {min_genes}"
        )

    out_dir.mkdir(parents=True, exist_ok=False)
    # Match the repository's existing Donor.load contract, which uses
    # ``split("\\n")`` and therefore cannot tolerate a trailing blank label.
    (out_dir / "perts.txt").write_text("\n".join(ordered_perts))
    (out_dir / "genes.txt").write_text("\n".join(genes))

    obs_columns = [column for column in OBS_COLUMNS if column in obs.columns]
    selected_obs = obs.iloc[rows][obs_columns].copy()
    selected_obs.insert(
        0, "source_obs_name", [clean_string(obs.index[i]) for i in rows]
    )
    selected_obs.insert(1, "target_gene", ordered_perts)
    selected_obs.to_csv(out_dir / "obs.csv", index=False)

    for source_name, output_name in LAYER_OUTPUTS.items():
        if source_name not in source["layers"]:
            raise ValueError(f"source is missing required layer {source_name!r}")
        destination = np.lib.format.open_memmap(
            out_dir / output_name,
            mode="w+",
            dtype=np.float32,
            shape=(len(rows), len(gene_columns)),
        )
        destination[:] = 0
        read_matrix_subset_into(
            source["layers"][source_name],
            rows,
            gene_columns,
            destination,
            row_block=row_block,
            max_rss_gib=max_rss_gib,
        )
        destination.flush()
        del destination
        print(f"  {condition:8} wrote {output_name}", flush=True)

    log_fc = np.load(out_dir / "log_fc.npy", mmap_mode="r")
    delta = np.lib.format.open_memmap(
        out_dir / "delta.npy",
        mode="w+",
        dtype=np.float32,
        shape=log_fc.shape,
    )
    valid = np.lib.format.open_memmap(
        out_dir / "valid_log_fc.npy",
        mode="w+",
        dtype=np.bool_,
        shape=log_fc.shape,
    )
    finite_count = 0
    for lo in range(0, len(rows), row_block):
        hi = min(lo + row_block, len(rows))
        block = np.asarray(log_fc[lo:hi], dtype=np.float32)
        finite = np.isfinite(block)
        valid[lo:hi] = finite
        delta[lo:hi] = np.where(finite, block * LN2, 0.0).astype(np.float32)
        finite_count += int(finite.sum())
        enforce_memory_limit(max_rss_gib)
    delta.flush()
    valid.flush()
    finite_fraction = finite_count / delta.size
    if finite_fraction < min_finite_fraction:
        raise ValueError(
            f"{condition}: finite log_fc fraction {finite_fraction:.4f} "
            f"< required {min_finite_fraction:.4f}"
        )

    # Donor.load expects b_ctrl.npy.  DESeq2 baseMean is not a control-only
    # profile, so this is explicitly a compatibility proxy and must not be used
    # for basal-context matching or empirical-Bayes noise calibration.
    base_mean = np.load(out_dir / "base_mean.npy", mmap_mode="r")
    with np.errstate(all="ignore"):
        base_proxy = np.nanmedian(np.asarray(base_mean, dtype=np.float64), axis=0)
    base_proxy = np.where(np.isfinite(base_proxy) & (base_proxy > 0), base_proxy, 0.0)
    if base_proxy.sum() <= 0:
        raise ValueError(f"{condition}: baseMean-derived basal proxy is degenerate")
    b_ctrl = np.log1p(TARGET_SUM * base_proxy / base_proxy.sum()).astype(np.float32)
    np.save(out_dir / "b_ctrl.npy", b_ctrl)

    effect_norms = np.linalg.norm(np.asarray(delta, dtype=np.float64), axis=1)
    cell_counts = (
        pd.to_numeric(selected_obs["n_cells_target"], errors="coerce")
        .dropna()
        .to_numpy()
        if "n_cells_target" in selected_obs
        else np.asarray([], dtype=np.float64)
    )
    meta = {
        "source": SOURCE_URL,
        "source_local_path": str(source.filename),
        "condition": condition,
        "n_perturbations": len(ordered_perts),
        "n_genes": len(genes),
        "panel_targets_missing": [
            target for target in panel if target not in target_rows
        ],
        "target_sum": TARGET_SUM,
        "effect_definition": "delta = ln(2) * source DESeq2 log2 fold change",
        "nonfinite_log_fc_policy": "stored in valid_log_fc.npy and replaced by 0 in delta.npy",
        "finite_log_fc_fraction": finite_fraction,
        "b_ctrl_kind": "median_DESeq2_baseMean_proxy_not_control_only",
        "b_ctrl_allowed_uses": ["Donor.load compatibility", "rough expression gating"],
        "b_ctrl_forbidden_uses": [
            "basal context matching",
            "empirical-Bayes count-noise calibration",
        ],
        "effect_norm": {
            "min": float(effect_norms.min()),
            "median": float(np.median(effect_norms)),
            "max": float(effect_norms.max()),
        },
        "cells_per_pert": (
            {
                "min": int(np.min(cell_counts)),
                "median": int(np.median(cell_counts)),
                "max": int(np.max(cell_counts)),
            }
            if cell_counts.size
            else None
        ),
        "layers": LAYER_OUTPUTS,
        "peak_rss_gib": peak_rss_gib(),
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    del log_fc, delta, valid, base_mean
    return meta


def run_extraction(args: argparse.Namespace) -> dict[str, Any]:
    input_path = Path(args.input).expanduser().resolve()
    panel_path = Path(args.panel_csv).expanduser().resolve()
    genes_path = Path(args.genes_csv).expanduser().resolve()
    out_root = Path(args.output).expanduser().resolve()
    for path in (input_path, panel_path, genes_path):
        if not path.exists():
            raise FileNotFoundError(path)

    panel = (
        pd.read_csv(panel_path)["target_gene"].dropna().astype(str).str.strip().tolist()
    )
    official_genes = (
        pd.read_csv(genes_path)["gene_name"].dropna().astype(str).str.strip().tolist()
    )
    if len(panel) != len(set(panel)):
        raise ValueError("VCC panel contains duplicate target symbols")
    if len(official_genes) != len(set(official_genes)):
        raise ValueError("VCC gene axis contains duplicate symbols")

    if not args.dry_run and out_root.exists():
        if not args.overwrite:
            raise FileExistsError(f"{out_root} exists; pass --overwrite to replace it")

    with h5py.File(input_path, "r") as source:
        for key in ("obs", "var", "layers"):
            if key not in source:
                raise ValueError(f"source H5AD is missing /{key}")
        obs = read_elem(source["obs"])
        var = read_elem(source["var"])
        if not isinstance(obs, pd.DataFrame) or not isinstance(var, pd.DataFrame):
            raise TypeError("source /obs and /var did not decode as pandas DataFrames")
        expected_shape = (len(obs), len(var))
        for layer in LAYER_OUTPUTS:
            if layer not in source["layers"]:
                raise ValueError(f"source is missing /layers/{layer}")
            if matrix_shape(source["layers"][layer]) != expected_shape:
                raise ValueError(
                    f"/layers/{layer} shape {matrix_shape(source['layers'][layer])} "
                    f"!= obs/var shape {expected_shape}"
                )

        selected, row_audit = resolve_target_rows(obs, var, panel)
        all_rows = np.asarray(
            sorted(
                {row for condition in selected.values() for row in condition.values()}
            ),
            dtype=np.int64,
        )
        if not len(all_rows):
            raise ValueError("no VCC-panel perturbations could be mapped")
        genes, gene_columns, gene_audit = choose_gene_columns(
            var,
            official_genes,
            source["layers"]["baseMean"],
            all_rows,
            row_block=args.row_block,
            max_rss_gib=args.max_rss_gib,
        )

        coverage = {condition: len(selected[condition]) for condition in CONDITIONS}
        union_targets = sorted(
            set().union(*(set(selected[condition]) for condition in CONDITIONS))
        )
        if len(union_targets) < args.min_union_targets:
            raise ValueError(
                f"only {len(union_targets)} panel targets map across any condition; "
                f"require {args.min_union_targets}"
            )
        below_target_threshold = {
            condition: count
            for condition, count in coverage.items()
            if count < args.min_targets
        }
        if below_target_threshold:
            raise ValueError(
                f"condition coverage below {args.min_targets}: {below_target_threshold}"
            )
        if len(genes) < args.min_genes:
            raise ValueError(
                f"only {len(genes)} official genes map; require {args.min_genes}"
            )
        print(
            f"source {expected_shape[0]:,} perturbation-condition rows x "
            f"{expected_shape[1]:,} genes",
            flush=True,
        )
        print(
            f"panel coverage {coverage}; mapped official genes {len(genes):,}",
            flush=True,
        )
        if args.dry_run:
            return {
                "dry_run": True,
                "source_shape": expected_shape,
                "coverage": coverage,
                "union_coverage": len(union_targets),
                "n_genes": len(genes),
                "peak_rss_gib": peak_rss_gib(),
            }

        if out_root.exists():
            shutil.rmtree(out_root)
        out_root.mkdir(parents=True)
        pd.DataFrame(row_audit).to_csv(out_root / "target_row_map.csv", index=False)
        pd.DataFrame(gene_audit).to_csv(out_root / "gene_map.csv", index=False)
        condition_meta = {}
        for condition in CONDITIONS:
            condition_meta[condition] = write_condition(
                source,
                obs,
                condition,
                selected[condition],
                panel,
                genes,
                gene_columns,
                out_root / condition_slug(condition),
                row_block=args.row_block,
                max_rss_gib=args.max_rss_gib,
                min_targets=args.min_targets,
                min_genes=args.min_genes,
                min_finite_fraction=args.min_finite_fraction,
            )

    manifest = {
        "source": SOURCE_URL,
        "source_local_path": str(input_path),
        "source_size_bytes": input_path.stat().st_size,
        "panel_csv": str(panel_path),
        "panel_csv_sha256": sha256_file(panel_path),
        "genes_csv": str(genes_path),
        "genes_csv_sha256": sha256_file(genes_path),
        "conditions": condition_meta,
        "panel_union_coverage": len(union_targets),
        "panel_union_missing": [
            target for target in panel if target not in union_targets
        ],
        "mapped_official_genes": len(genes),
        "row_block": args.row_block,
        "max_rss_gib": args.max_rss_gib,
        "peak_rss_gib": peak_rss_gib(),
    }
    (out_root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        required=True,
        help="cloud-local path to GWCD4i.DE_stats.h5ad; s3:// URLs are not downloaded",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="output root for cd4_rest/cd4_stim8hr/cd4_stim48hr",
    )
    parser.add_argument(
        "--panel-csv",
        default="data/challenge/pert_counts.csv",
        help="CSV with target_gene column",
    )
    parser.add_argument(
        "--genes-csv",
        default="data/challenge/gene_names.csv",
        help="CSV with canonical gene_name column",
    )
    parser.add_argument("--row-block", type=int, default=32)
    parser.add_argument(
        "--max-rss-gib",
        type=float,
        default=4.0,
        help="abort after a block if process peak RSS exceeds this value; 0 disables",
    )
    parser.add_argument("--min-targets", type=int, default=275)
    parser.add_argument("--min-union-targets", type=int, default=290)
    parser.add_argument("--min-genes", type=int, default=9_500)
    parser.add_argument("--min-finite-fraction", type=float, default=0.98)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="read only metadata/mappings; do not create output files",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.row_block < 1:
        raise SystemExit("--row-block must be positive")
    result = run_extraction(args)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()

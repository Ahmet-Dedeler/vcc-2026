#!/usr/bin/env python3
"""Low-memory source-vector ablation test for reported seed off-targets."""

from __future__ import annotations

import os

# Freeze numerical libraries before importing NumPy.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import argparse
import csv
import gc
import gzip
import hashlib
import io
import json
import math
import resource
import struct
import sys
import urllib.request
import zipfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
from numpy.lib import format as npformat


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
MAX_RSS_GIB = 0.75
STOP_RSS_GIB = 0.70
MAX_DOWNLOAD_BYTES = 50_000_000
EFFECT_SCALE = 6.817
WHITEN_GAMMA = 0.65
PM_FRACTION = 0.50
CD4_FRACTION = 0.05
NULL_CONTROLS = 64
NULL_REPS = 20_000
SEED = 6_802_026

PANEL_CSV = ROOT / "data/challenge/pert_counts.csv"
GENES_CSV = ROOT / "data/challenge/gene_names.csv"
K562_ROOT = ROOT / "data/donors/k562_gwps"
HCT_ROOT = ROOT / "data/donors/orion_hct116"
HEK_ROOT = ROOT / "data/donors/orion_hek293t"
PM_ROOT = ROOT / "data/external/perturbmatch"
KOLF_ROOT = ROOT / "data/external/kolf2"
CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
EXP064_JSON = (
    ROOT / "experiments/064-pm-cd4-composition/cloud-output/results.json"
)
EXP064_PDS = (
    ROOT / "experiments/064-pm-cd4-composition/cloud-output/per_target_pds.csv"
)
WHITENER_ROOT = (
    ROOT / "experiments/062-cd4-post-whiten/cloud-output/whiteners"
)

SOURCE_REPO = "AustinHartman/perturb_seed_reproduce"
SOURCE_COMMIT = "198fe3d84263002e9cfa17e1a2178d4e8dd0a669"
SOURCE_BASE = (
    "https://raw.githubusercontent.com/"
    f"{SOURCE_REPO}/{SOURCE_COMMIT}/data/"
)
DOI = "10.64898/2026.03.27.714658"
EXPECTED_TITLE = (
    "Systematic identification of seed-driven off-target effects in "
    "Perturb-seq experiments"
)

SOURCE_TABLES = {
    "replogle_k562": {
        "file": "replogle_k562_2kb_filt_no_lfc_cutoff_off_target_candidates.csv.gz",
        "git_blob_sha1": "186989c1eec7238eb55ef64053f2eb76aa227628",
        "git_bytes": 14_460_423,
        "pairs": (("PQBP1", "TIMM17B"),),
    },
    "xaira_hct116": {
        "file": "xaira_hct116_2kb_filt_no_lfc_cutoff_off_target_candidates.csv.gz",
        "git_blob_sha1": "e4feb6a8d219089f6dc68b32aaa87c3699ec7fee",
        "git_bytes": 21_669_263,
        "pairs": (
            ("PQBP1", "TIMM17B"),
            ("SHPRH", "CENPP"),
            ("FZD2", "P2RX7"),
        ),
    },
    "zhu_cd4": {
        "file": "zhu_cd4t_d1_8hr_2kb_filt_no_lfc_cutoff_off_target_candidates.csv.gz",
        "git_blob_sha1": "00941de24b282a9d2d2ccabfd13aaf8a14127189",
        "git_bytes": 15_925_175,
        "pairs": (("GLYR1", "UBN1"),),
    },
}

EXPECTED_WHITENER_SHA256 = {
    "hct_geometry": "7cecbf7f36e40d22a6bea1ff3adc86a8788b329def397044fb5b9ba675247cbc",
    "hek_geometry": "97e486e0d1a805e1f24a622f1cefbfb0917e9282c5d0491c4ddb41c8475f2af1",
}


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def enforce_rss(stage: str) -> None:
    peak = peak_rss_gib()
    if peak >= STOP_RSS_GIB:
        raise MemoryError(
            f"{stage}: peak RSS {peak:.6f} GiB reached stop threshold "
            f"{STOP_RSS_GIB:.6f} GiB (hard cap {MAX_RSS_GIB:.6f} GiB)"
        )


def sha256_file(path: Path, block_bytes: int = 4 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def read_lines(path: Path, *, require_unique: bool = True) -> tuple[str, ...]:
    values = tuple(line.strip() for line in path.read_text().splitlines() if line.strip())
    if require_unique and len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: missing column {column}")
        values = tuple(
            str(row[column]).strip()
            for row in reader
            if row.get(column) and str(row[column]).strip()
        )
    if len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate {column}")
    return values


def write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        return
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                fields.append(key)
                seen.add(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


@dataclass(frozen=True)
class NpyHeader:
    shape: tuple[int, ...]
    fortran_order: bool
    dtype: np.dtype
    offset: int


def read_npy_header(path: Path) -> NpyHeader:
    with path.open("rb") as handle:
        version = npformat.read_magic(handle)
        readers = {
            (1, 0): npformat.read_array_header_1_0,
            (2, 0): npformat.read_array_header_2_0,
            (3, 0): npformat.read_array_header_2_0,
        }
        if version not in readers:
            raise ValueError(f"{path}: unsupported npy version {version}")
        shape, fortran_order, dtype = readers[version](handle)
        return NpyHeader(
            shape=tuple(map(int, shape)),
            fortran_order=bool(fortran_order),
            dtype=np.dtype(dtype),
            offset=handle.tell(),
        )


@dataclass(frozen=True)
class Bank:
    name: str
    path: Path
    genes: tuple[str, ...]
    perts: tuple[str, ...]
    gene_pos: dict[str, int]
    pert_pos: dict[str, int]
    header: NpyHeader

    @classmethod
    def load(cls, name: str, path: Path) -> "Bank":
        path = path.resolve()
        genes = read_lines(path / "genes.txt", require_unique=False)
        perts = read_lines(path / "perts.txt")
        header = read_npy_header(path / "delta.npy")
        if header.shape != (len(perts), len(genes)):
            raise ValueError(
                f"{name}: shape {header.shape} != {(len(perts), len(genes))}"
            )
        return cls(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            gene_pos={gene: index for index, gene in enumerate(genes)},
            pert_pos={pert: index for index, pert in enumerate(perts)},
            header=header,
        )

    @property
    def delta_path(self) -> Path:
        return self.path / "delta.npy"


@dataclass(frozen=True)
class AxisMap:
    all_source_columns: np.ndarray
    unique_source_columns: np.ndarray
    unique_official_positions: np.ndarray
    unique_positions_in_all: np.ndarray
    duplicate_symbols: tuple[str, ...]


def map_bank(bank: Bank, official: Sequence[str]) -> AxisMap:
    official_pos = {gene: index for index, gene in enumerate(official)}
    all_source: list[int] = []
    last_source: dict[int, int] = {}
    duplicates: set[str] = set()
    seen: set[int] = set()
    for source_column, gene in enumerate(bank.genes):
        destination = official_pos.get(gene)
        if destination is None:
            continue
        all_source.append(source_column)
        if destination in seen:
            duplicates.add(gene)
        seen.add(destination)
        last_source[destination] = source_column
    unique_official = np.asarray(sorted(last_source), dtype=np.int64)
    unique_source = np.asarray(
        [last_source[position] for position in unique_official], dtype=np.int64
    )
    position_in_all = {source: index for index, source in enumerate(all_source)}
    unique_in_all = np.asarray(
        [position_in_all[int(source)] for source in unique_source], dtype=np.int64
    )
    return AxisMap(
        all_source_columns=np.asarray(all_source, dtype=np.int64),
        unique_source_columns=unique_source,
        unique_official_positions=unique_official,
        unique_positions_in_all=unique_in_all,
        duplicate_symbols=tuple(sorted(duplicates)),
    )


def iter_c_rows(bank: Bank, chunk_rows: int) -> Iterable[tuple[int, np.ndarray]]:
    if bank.header.fortran_order:
        raise ValueError(f"{bank.name}: expected C-order")
    n_rows, n_cols = bank.header.shape
    row_bytes = n_cols * bank.header.dtype.itemsize
    with bank.delta_path.open("rb") as handle:
        handle.seek(bank.header.offset)
        for lo in range(0, n_rows, chunk_rows):
            count = min(chunk_rows, n_rows - lo)
            raw = handle.read(count * row_bytes)
            if len(raw) != count * row_bytes:
                raise EOFError(f"{bank.delta_path}: truncated at row {lo}")
            yield lo, np.frombuffer(raw, dtype=bank.header.dtype).reshape(
                count, n_cols
            )


def read_fortran_column(bank: Bank, source_column: int) -> np.ndarray:
    if not bank.header.fortran_order:
        raise ValueError(f"{bank.name}: expected Fortran-order")
    n_rows, _ = bank.header.shape
    column_bytes = n_rows * bank.header.dtype.itemsize
    with bank.delta_path.open("rb") as handle:
        handle.seek(bank.header.offset + source_column * column_bytes)
        raw = handle.read(column_bytes)
    if len(raw) != column_bytes:
        raise EOFError(f"{bank.delta_path}: truncated column {source_column}")
    return np.frombuffer(raw, dtype=bank.header.dtype)


@dataclass
class PreparedMain:
    name: str
    effects: np.ndarray
    available: np.ndarray
    mean_on_axis: np.ndarray
    raw_median_norm: float
    mapped_gene_count: int
    duplicate_symbols: tuple[str, ...]


def prepare_main(
    bank: Bank,
    panel: Sequence[str],
    official: Sequence[str],
    panel_positions: np.ndarray,
    *,
    chunk_rows: int,
) -> PreparedMain:
    """Reproduce experiment-064 main-source centering and scaling."""
    mapping = map_bank(bank, official)
    n_rows = bank.header.shape[0]
    target_rows = np.asarray(
        [bank.pert_pos.get(target, -1) for target in panel], dtype=np.int64
    )
    available = target_rows >= 0
    raw_targets = np.zeros((len(panel), len(official)), dtype=np.float64)
    mean_on_axis = np.zeros(len(official), dtype=np.float64)

    if not bank.header.fortran_order:
        mean_all = np.zeros(len(mapping.all_source_columns), dtype=np.float64)
        for _, whole in iter_c_rows(bank, chunk_rows):
            mean_all += np.asarray(
                whole[:, mapping.all_source_columns], dtype=np.float64
            ).sum(axis=0, dtype=np.float64)
            enforce_rss(f"{bank.name} mean")
        mean_all /= n_rows

        norms = np.empty(n_rows, dtype=np.float64)
        for lo, whole in iter_c_rows(bank, chunk_rows):
            selected = np.asarray(
                whole[:, mapping.all_source_columns], dtype=np.float64
            )
            centered = selected - mean_all
            norms[lo : lo + len(whole)] = np.linalg.norm(centered, axis=1)
            in_block = np.flatnonzero(
                (target_rows >= lo) & (target_rows < lo + len(whole))
            )
            for destination in in_block.tolist():
                local = int(target_rows[destination] - lo)
                raw_targets[
                    destination, mapping.unique_official_positions
                ] = centered[local, mapping.unique_positions_in_all]
            enforce_rss(f"{bank.name} centered norms")
        unique_mean = mean_all[mapping.unique_positions_in_all]
    else:
        norm_sums = np.zeros(n_rows, dtype=np.float64)
        unique_source_to_index = {
            int(source): index
            for index, source in enumerate(mapping.unique_source_columns.tolist())
        }
        unique_mean = np.empty(len(mapping.unique_source_columns), dtype=np.float64)
        for source_column in mapping.all_source_columns.tolist():
            column = np.asarray(
                read_fortran_column(bank, int(source_column)), dtype=np.float64
            )
            mean = float(column.mean())
            centered = column - mean
            norm_sums += centered * centered
            unique_index = unique_source_to_index.get(int(source_column))
            if unique_index is not None:
                unique_mean[unique_index] = mean
                present = np.flatnonzero(available)
                raw_targets[
                    present,
                    mapping.unique_official_positions[unique_index],
                ] = centered[target_rows[present]]
            enforce_rss(f"{bank.name} column {source_column}")
        norms = np.sqrt(norm_sums)

    raw_median = float(np.median(norms)) or 1.0
    mean_on_axis[mapping.unique_official_positions] = unique_mean
    raw_targets *= EFFECT_SCALE / raw_median
    raw_targets[:, panel_positions] = 0.0
    mean_on_axis[panel_positions] = 0.0
    enforce_rss(f"prepared {bank.name}")
    return PreparedMain(
        name=bank.name,
        effects=raw_targets,
        available=available,
        mean_on_axis=mean_on_axis,
        raw_median_norm=raw_median,
        mapped_gene_count=len(mapping.unique_official_positions),
        duplicate_symbols=mapping.duplicate_symbols,
    )


@dataclass
class PreparedCD4:
    mean_effect: np.ndarray
    gene_votes: np.ndarray
    condition_available: dict[str, np.ndarray]
    details: list[dict[str, Any]]


def prepare_cd4(
    banks: Sequence[Bank],
    panel: Sequence[str],
    official: Sequence[str],
    panel_set: set[str],
    *,
    chunk_rows: int,
) -> PreparedCD4:
    official_pos = {gene: index for index, gene in enumerate(official)}
    value_sum = np.zeros((len(panel), len(official)), dtype=np.float64)
    gene_votes = np.zeros(value_sum.shape, dtype=np.uint8)
    condition_available: dict[str, np.ndarray] = {}
    details: list[dict[str, Any]] = []

    for bank in banks:
        columns: list[int] = []
        positions: list[int] = []
        for source_column, gene in enumerate(bank.genes):
            destination = official_pos.get(gene)
            if destination is None or gene in panel_set:
                continue
            columns.append(source_column)
            positions.append(destination)
        source_columns = np.asarray(columns, dtype=np.int64)
        official_positions = np.asarray(positions, dtype=np.int64)
        if len(np.unique(official_positions)) != len(official_positions):
            raise ValueError(f"{bank.name}: duplicate CD4 official genes")

        n_rows = bank.header.shape[0]
        means = np.zeros(len(source_columns), dtype=np.float64)
        for _, whole in iter_c_rows(bank, chunk_rows):
            means += np.asarray(whole[:, source_columns], dtype=np.float64).sum(
                axis=0, dtype=np.float64
            )
        means /= n_rows

        absolute = np.zeros(len(source_columns), dtype=np.float64)
        raw_norms = np.empty(n_rows, dtype=np.float64)
        for lo, whole in iter_c_rows(bank, chunk_rows):
            centered = (
                np.asarray(whole[:, source_columns], dtype=np.float64) - means
            )
            raw_norms[lo : lo + len(whole)] = np.linalg.norm(centered, axis=1)
            absolute += np.abs(centered).sum(axis=0, dtype=np.float64)
            enforce_rss(f"{bank.name} gene scale")
        mean_abs = absolute / n_rows
        live = mean_abs[mean_abs > 0]
        floor = float(np.percentile(live, 25.0)) if live.size else 1.0
        gene_scale = np.maximum(mean_abs, floor)

        normalized_norms = np.empty(n_rows, dtype=np.float64)
        target_rows = np.asarray(
            [bank.pert_pos.get(target, -1) for target in panel], dtype=np.int64
        )
        available = target_rows >= 0
        target_values = np.zeros((len(panel), len(source_columns)), dtype=np.float64)
        for lo, whole in iter_c_rows(bank, chunk_rows):
            centered = (
                np.asarray(whole[:, source_columns], dtype=np.float64) - means
            ) / gene_scale
            normalized_norms[lo : lo + len(whole)] = np.linalg.norm(
                centered, axis=1
            )
            in_block = np.flatnonzero(
                (target_rows >= lo) & (target_rows < lo + len(whole))
            )
            for destination in in_block.tolist():
                target_values[destination] = centered[
                    int(target_rows[destination] - lo)
                ]
            enforce_rss(f"{bank.name} normalized rows")
        normalized_median = float(np.median(normalized_norms)) or 1.0
        target_values *= EFFECT_SCALE / normalized_median

        present = np.flatnonzero(available)
        value_sum[np.ix_(present, official_positions)] += target_values[present]
        gene_votes[np.ix_(present, official_positions)] += 1
        condition_available[bank.name] = available
        details.append(
            {
                "condition": bank.name,
                "background_perturbations": n_rows,
                "mapped_non_panel_genes": len(source_columns),
                "covered_panel_targets": int(available.sum()),
                "raw_median_centered_norm": float(np.median(raw_norms)),
                "gene_normalized_median_before_rescale": normalized_median,
            }
        )
        del target_values
        enforce_rss(f"prepared {bank.name}")

    mean_effect = np.zeros_like(value_sum)
    live = gene_votes > 0
    mean_effect[live] = value_sum[live] / gene_votes[live]
    return PreparedCD4(
        mean_effect=mean_effect,
        gene_votes=gene_votes,
        condition_available=condition_available,
        details=details,
    )


@dataclass(frozen=True)
class StoredNpy:
    path: Path
    member: str
    shape: tuple[int, ...]
    fortran_order: bool
    dtype: np.dtype
    data_offset: int
    nbytes: int


def stored_npy_member(path: Path, member: str) -> StoredNpy:
    with zipfile.ZipFile(path) as archive:
        info = archive.getinfo(member)
        if info.compress_type != zipfile.ZIP_STORED:
            raise ValueError(f"{path}:{member} is compressed; bounded streaming unavailable")
        with path.open("rb") as handle:
            handle.seek(info.header_offset)
            local = handle.read(30)
            fields = struct.unpack("<IHHHHHIIIHH", local)
            signature, name_len, extra_len = fields[0], fields[-2], fields[-1]
            if signature != 0x04034B50:
                raise ValueError(f"{path}:{member} bad local ZIP header")
            handle.seek(name_len + extra_len, io.SEEK_CUR)
            if handle.tell() != info.header_offset + 30 + name_len + extra_len:
                raise AssertionError("ZIP member offset mismatch")
            version = npformat.read_magic(handle)
            readers = {
                (1, 0): npformat.read_array_header_1_0,
                (2, 0): npformat.read_array_header_2_0,
                (3, 0): npformat.read_array_header_2_0,
            }
            shape, fortran_order, dtype = readers[version](handle)
            data_offset = handle.tell()
    dtype = np.dtype(dtype)
    nbytes = int(math.prod(shape) * dtype.itemsize)
    if data_offset + nbytes > path.stat().st_size:
        raise ValueError(f"{path}:{member} exceeds archive size")
    return StoredNpy(
        path=path,
        member=member,
        shape=tuple(map(int, shape)),
        fortran_order=bool(fortran_order),
        dtype=dtype,
        data_offset=data_offset,
        nbytes=nbytes,
    )


@dataclass(frozen=True)
class Operator:
    name: str
    path: Path
    evals: np.ndarray
    evecs: StoredNpy
    wcols: np.ndarray
    preserve_scale: float
    metadata: dict[str, Any]


def load_operator(
    name: str,
    path: Path,
    source_banks: Sequence[Bank],
    official: Sequence[str],
    panel_set: set[str],
    preserve_scale: float,
) -> Operator:
    expected_hash = EXPECTED_WHITENER_SHA256[name]
    actual_hash = sha256_file(path)
    if actual_hash != expected_hash:
        raise ValueError(f"{name}: whitener hash {actual_hash} != {expected_hash}")
    metadata = json.loads(path.with_suffix(".json").read_text())
    donor_union = set().union(*(set(bank.genes) for bank in source_banks))
    wcols = np.asarray(
        [
            index
            for index, gene in enumerate(official)
            if gene in donor_union and gene not in panel_set
        ],
        dtype=np.int64,
    )
    if metadata.get("wcols_sha256") != hashlib.sha256(wcols.tobytes()).hexdigest():
        raise ValueError(f"{name}: operator axis mismatch")
    eval_header = stored_npy_member(path, "evals.npy")
    with path.open("rb") as handle:
        handle.seek(eval_header.data_offset)
        raw = handle.read(eval_header.nbytes)
    evals = np.frombuffer(raw, dtype=eval_header.dtype).copy()
    evecs = stored_npy_member(path, "evecs.npy")
    if evecs.fortran_order or evecs.shape != (len(wcols), len(wcols)):
        raise ValueError(f"{name}: unsupported eigenvector layout {evecs}")
    return Operator(
        name=name,
        path=path,
        evals=evals,
        evecs=evecs,
        wcols=wcols,
        preserve_scale=preserve_scale,
        metadata=metadata,
    )


def iter_evec_row_blocks(
    stored: StoredNpy, block_rows: int
) -> Iterable[tuple[int, np.ndarray]]:
    n_rows, n_cols = stored.shape
    row_bytes = n_cols * stored.dtype.itemsize
    with stored.path.open("rb") as handle:
        for lo in range(0, n_rows, block_rows):
            count = min(block_rows, n_rows - lo)
            handle.seek(stored.data_offset + lo * row_bytes)
            raw = handle.read(count * row_bytes)
            if len(raw) != count * row_bytes:
                raise EOFError(f"{stored.path}:{stored.member} truncated at row {lo}")
            yield lo, np.frombuffer(raw, dtype=stored.dtype).reshape(count, n_cols)


def transform_rows(
    matrix: np.ndarray,
    operator: Operator,
    *,
    block_rows: int,
) -> np.ndarray:
    """Apply V diag(lambda^(-gamma/2)) V.T with two bounded sequential scans."""
    values = np.asarray(matrix, dtype=np.float64)
    n_queries, n_cols = values.shape
    if n_cols != len(operator.evals):
        raise ValueError(f"{operator.name}: matrix width {n_cols} != {len(operator.evals)}")
    projected = np.zeros((n_queries, n_cols), dtype=np.float64)
    for lo, block in iter_evec_row_blocks(operator.evecs, block_rows):
        hi = lo + len(block)
        projected += values[:, lo:hi] @ block
        enforce_rss(f"{operator.name} forward block {lo}")
    floor = float(operator.evals.max()) * 1e-12
    weights = np.maximum(operator.evals, floor) ** (-WHITEN_GAMMA / 2.0)
    projected *= weights
    output = np.empty_like(projected)
    for lo, block in iter_evec_row_blocks(operator.evecs, block_rows):
        hi = lo + len(block)
        output[:, lo:hi] = projected @ block.T
        enforce_rss(f"{operator.name} inverse block {lo}")
    return output


def prepare_pm_slot(
    k562: PreparedMain,
    pm: PreparedMain,
    k562_bank: Bank,
    pm_bank: Bank,
    official: Sequence[str],
    panel_set: set[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    shared = np.asarray(
        [
            index
            for index, gene in enumerate(official)
            if gene not in panel_set
            and gene in k562_bank.gene_pos
            and gene in pm_bank.gene_pos
        ],
        dtype=np.int64,
    )
    effects = k562.effects.copy()
    effects[:, shared] = (
        (1.0 - PM_FRACTION) * k562.effects[:, shared]
        + PM_FRACTION * pm.effects[:, shared]
    )
    available = k562.available & pm.available
    return effects, available, shared


def build_baseline_rows(
    target_indices: Sequence[int],
    pm_slot: np.ndarray,
    pm_available: np.ndarray,
    k562_mean: np.ndarray,
    companion: PreparedMain,
    *,
    removed: frozenset[str],
) -> np.ndarray:
    output = np.zeros((len(target_indices), pm_slot.shape[1]), dtype=np.float64)
    fallback_mean = (k562_mean + companion.mean_on_axis) / 2.0
    fallback_norm = float(np.linalg.norm(fallback_mean)) or 1.0
    fallback = EFFECT_SCALE * fallback_mean / fallback_norm
    for destination, target_index in enumerate(target_indices):
        pieces: list[np.ndarray] = []
        if "k562_slot" not in removed and pm_available[target_index]:
            pieces.append(pm_slot[target_index])
        if companion.name not in removed and companion.available[target_index]:
            pieces.append(companion.effects[target_index])
        if pieces:
            output[destination] = sum(pieces) / len(pieces)
        else:
            output[destination] = fallback
    return output


def source_signature(
    target_index: int,
    cases: Sequence[dict[str, Any]],
    prepared: dict[str, PreparedMain],
    pm_available: np.ndarray,
    cd4: PreparedCD4,
) -> tuple[Any, ...]:
    signature: list[Any] = [
        bool(prepared["k562_gwps"].available[target_index]),
        bool(pm_available[target_index]),
    ]
    for case in cases:
        companion = (
            "orion_hek293t"
            if case["config"] == "hct_geometry"
            else "orion_hct116"
        )
        signature.extend(
            [
                case["truth"],
                bool(prepared[case["truth"]].available[target_index]),
                companion,
                bool(prepared[companion].available[target_index]),
                tuple(sorted(case["removed"])),
            ]
        )
    signature.append(
        tuple(
            bool(cd4.condition_available[name][target_index])
            for name in CD4_NAMES
        )
    )
    return tuple(signature)


def choose_null_controls(
    target: str,
    cases: Sequence[dict[str, Any]],
    panel: Sequence[str],
    panel_pos: dict[str, int],
    prepared: dict[str, PreparedMain],
    pm_available: np.ndarray,
    cd4: PreparedCD4,
    *,
    seed: int,
) -> tuple[str, ...]:
    target_index = panel_pos[target]
    expected = source_signature(
        target_index, cases, prepared, pm_available, cd4
    )
    excluded = {"PQBP1", "SHPRH", "FZD2", "GLYR1"}
    eligible = [
        candidate
        for candidate in panel
        if candidate not in excluded
        and source_signature(
            panel_pos[candidate], cases, prepared, pm_available, cd4
        )
        == expected
    ]
    if not eligible:
        raise ValueError(f"{target}: no exact-availability matched controls")
    if len(eligible) <= NULL_CONTROLS:
        return tuple(eligible)
    rng = np.random.default_rng(seed)
    chosen = rng.choice(np.asarray(eligible), size=NULL_CONTROLS, replace=False)
    return tuple(sorted(map(str, chosen.tolist())))


def prediction_rows_for_config(
    config: str,
    targets_by_mask: dict[frozenset[str], tuple[str, ...]],
    panel_pos: dict[str, int],
    pm_slot: np.ndarray,
    pm_available: np.ndarray,
    k562: PreparedMain,
    companion: PreparedMain,
    cd4: PreparedCD4,
    operator: Operator,
    *,
    block_rows: int,
) -> dict[tuple[str, frozenset[str]], np.ndarray]:
    labels: list[tuple[str, frozenset[str]]] = []
    pre_rows: list[np.ndarray] = []
    for removed, targets in targets_by_mask.items():
        indices = [panel_pos[target] for target in targets]
        baselines = build_baseline_rows(
            indices,
            pm_slot,
            pm_available,
            k562.mean_on_axis,
            companion,
            removed=removed,
        )
        for target, row in zip(targets, baselines, strict=True):
            labels.append((target, removed))
            pre_rows.append(row[operator.wcols])
    matrix = np.stack(pre_rows, axis=0)
    del pre_rows
    transformed = transform_rows(matrix, operator, block_rows=block_rows)
    del matrix
    transformed *= operator.preserve_scale
    result: dict[tuple[str, frozenset[str]], np.ndarray] = {}
    for index, (target, removed) in enumerate(labels):
        panel_index = panel_pos[target]
        vector = transformed[index].copy()
        if "cd4" not in removed:
            vector += CD4_FRACTION * cd4.mean_effect[
                panel_index, operator.wcols
            ]
        result[(target, removed)] = vector
    enforce_rss(f"{config} predictions")
    return result


def rank_query(
    prediction: np.ndarray,
    held: PreparedMain,
    gallery_indices: np.ndarray,
    target_index_in_gallery: int,
    retrieval_official: np.ndarray,
    wcols: np.ndarray,
) -> dict[str, float]:
    lookup = {int(position): index for index, position in enumerate(wcols.tolist())}
    local_columns = np.asarray(
        [lookup[int(position)] for position in retrieval_official], dtype=np.int64
    )
    observed = held.effects[np.ix_(gallery_indices, retrieval_official)]
    observed /= np.maximum(
        np.linalg.norm(observed, axis=1, keepdims=True), 1e-12
    )
    query = np.asarray(prediction[local_columns], dtype=np.float64)
    query /= max(float(np.linalg.norm(query)), 1e-12)
    similarities = observed @ query
    own = float(similarities[target_index_in_gallery])
    ahead = int(np.sum(similarities > own))
    ties_other = int(np.sum(similarities == own) - 1)
    rank_zero = float(ahead + 0.5 * ties_other)
    pds = float(1.0 - rank_zero / (len(gallery_indices) - 1))
    return {
        "rank_1based": rank_zero + 1.0,
        "ahead": ahead,
        "ties_other": ties_other,
        "own_cosine": own,
        "pds": pds,
    }


class HashingReader:
    def __init__(self, raw: Any, expected_bytes: int):
        self.raw = raw
        self.expected_bytes = expected_bytes
        self.count = 0
        self.sha256 = hashlib.sha256()
        self.git_sha1 = hashlib.sha1(
            f"blob {expected_bytes}\0".encode()
        )

    def read(self, size: int = -1) -> bytes:
        data = self.raw.read(size)
        self.count += len(data)
        if self.count > MAX_DOWNLOAD_BYTES:
            raise ValueError(
                f"download exceeded {MAX_DOWNLOAD_BYTES} bytes"
            )
        self.sha256.update(data)
        self.git_sha1.update(data)
        return data


def verify_source_tables() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    evidence: list[dict[str, Any]] = []
    table_meta: dict[str, Any] = {}
    for source, spec in SOURCE_TABLES.items():
        expected_bytes = int(spec["git_bytes"])
        if expected_bytes > MAX_DOWNLOAD_BYTES:
            raise ValueError(f"{source}: source exceeds download cap")
        url = SOURCE_BASE + str(spec["file"])
        request = urllib.request.Request(
            url, headers={"User-Agent": "vcc-exp068/1.0"}
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            content_length = response.headers.get("Content-Length")
            if content_length and int(content_length) > MAX_DOWNLOAD_BYTES:
                raise ValueError(f"{source}: HTTP object exceeds download cap")
            hashing = HashingReader(response, expected_bytes)
            with gzip.GzipFile(fileobj=hashing) as compressed:
                text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
                reader = csv.DictReader(text)
                fields = tuple(reader.fieldnames or ())
                required = {
                    "gene",
                    "guide_id",
                    "target_log2fc",
                    "offtarget_gene",
                    "offtarget_log2fc",
                    "neighbor_rank",
                    "direction",
                    "seed",
                    "seed_match_len",
                    "guide_hamming_dist",
                }
                if not required.issubset(fields):
                    raise ValueError(f"{source}: unexpected fields")
                pairs = {tuple(pair) for pair in spec["pairs"]}
                found: dict[tuple[str, str], list[dict[str, str]]] = {
                    pair: [] for pair in pairs
                }
                rows = 0
                for row in reader:
                    rows += 1
                    pair = (
                        str(row["gene"]).upper(),
                        str(row["offtarget_gene"]).upper(),
                    )
                    if pair in found:
                        found[pair].append(row)
                text.detach()
        if hashing.count != expected_bytes:
            raise ValueError(
                f"{source}: downloaded {hashing.count} != {expected_bytes}"
            )
        if hashing.git_sha1.hexdigest() != spec["git_blob_sha1"]:
            raise ValueError(f"{source}: Git blob SHA-1 mismatch")
        for pair, matches in found.items():
            if len(matches) != 1:
                raise ValueError(f"{source}: {pair} matched {len(matches)} rows")
            row = matches[0]
            evidence.append(
                {
                    "source": source,
                    "target_gene": pair[0],
                    "offtarget_gene": pair[1],
                    **row,
                }
            )
        table_meta[source] = {
            "url": url,
            "repository_commit": SOURCE_COMMIT,
            "git_blob_sha1": spec["git_blob_sha1"],
            "compressed_bytes": hashing.count,
            "compressed_sha256": hashing.sha256.hexdigest(),
            "data_rows": rows,
            "columns": list(fields),
        }
        enforce_rss(f"source verification {source}")

    crossref_url = f"https://api.crossref.org/works/{DOI}"
    request = urllib.request.Request(
        crossref_url, headers={"User-Agent": "vcc-exp068/1.0"}
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        metadata = json.load(response)["message"]
    title = str(metadata["title"][0])
    if title != EXPECTED_TITLE:
        raise ValueError(f"DOI title mismatch: {title}")
    publication = metadata.get("published", {}).get("date-parts", [[None]])[0]
    provenance = {
        "paper": {
            "doi": DOI,
            "title": title,
            "authors": [
                " ".join(
                    part
                    for part in (author.get("given"), author.get("family"))
                    if part
                )
                for author in metadata.get("author", [])
            ],
            "published_date_parts": publication,
            "crossref_url": crossref_url,
        },
        "repository": {
            "name": SOURCE_REPO,
            "commit": SOURCE_COMMIT,
            "readme": (
                "https://raw.githubusercontent.com/"
                f"{SOURCE_REPO}/{SOURCE_COMMIT}/README.md"
            ),
            "figure_notebook": (
                "https://raw.githubusercontent.com/"
                f"{SOURCE_REPO}/{SOURCE_COMMIT}/scripts/figure2.ipynb"
            ),
        },
        "tables": table_meta,
    }
    return evidence, provenance


def build_cases() -> dict[str, list[dict[str, Any]]]:
    return {
        "PQBP1": [
            {
                "case_id": "hct_truth__hct_geometry",
                "truth": "orion_hct116",
                "config": "hct_geometry",
                "removed": frozenset({"k562_slot"}),
                "tier": "native_064",
            },
            {
                "case_id": "hek_truth__hek_geometry",
                "truth": "orion_hek293t",
                "config": "hek_geometry",
                "removed": frozenset({"k562_slot", "orion_hct116"}),
                "tier": "native_064",
            },
            {
                "case_id": "kolf_truth__hct_geometry",
                "truth": "kolf",
                "config": "hct_geometry",
                "removed": frozenset({"k562_slot"}),
                "tier": "secondary_kolf",
            },
            {
                "case_id": "kolf_truth__hek_geometry",
                "truth": "kolf",
                "config": "hek_geometry",
                "removed": frozenset({"k562_slot", "orion_hct116"}),
                "tier": "secondary_kolf",
            },
        ],
        "SHPRH": [
            {
                "case_id": "hek_truth__hek_geometry",
                "truth": "orion_hek293t",
                "config": "hek_geometry",
                "removed": frozenset({"orion_hct116"}),
                "tier": "native_064",
            },
            {
                "case_id": "kolf_truth__hek_geometry",
                "truth": "kolf",
                "config": "hek_geometry",
                "removed": frozenset({"orion_hct116"}),
                "tier": "secondary_kolf",
            },
        ],
        "FZD2": [
            {
                "case_id": "hek_truth__hek_geometry",
                "truth": "orion_hek293t",
                "config": "hek_geometry",
                "removed": frozenset({"orion_hct116"}),
                "tier": "native_064",
            },
            {
                "case_id": "kolf_truth__hek_geometry",
                "truth": "kolf",
                "config": "hek_geometry",
                "removed": frozenset({"orion_hct116"}),
                "tier": "secondary_kolf",
            },
        ],
        "GLYR1": [
            {
                "case_id": "kolf_truth__hct_geometry",
                "truth": "kolf",
                "config": "hct_geometry",
                "removed": frozenset({"cd4"}),
                "tier": "secondary_kolf",
            },
            {
                "case_id": "kolf_truth__hek_geometry",
                "truth": "kolf",
                "config": "hek_geometry",
                "removed": frozenset({"cd4"}),
                "tier": "secondary_kolf",
            },
        ],
    }


def summarize_null(actual: float, controls: np.ndarray) -> dict[str, Any]:
    controls = np.asarray(controls, dtype=np.float64)
    return {
        "controls": int(len(controls)),
        "mean": float(controls.mean()),
        "median": float(np.median(controls)),
        "p05": float(np.percentile(controls, 5)),
        "p95": float(np.percentile(controls, 95)),
        "min": float(controls.min()),
        "max": float(controls.max()),
        "one_sided_p_ge_actual": float(
            (1 + np.sum(controls >= actual)) / (len(controls) + 1)
        ),
        "actual_percentile_midrank": float(
            100.0
            * (
                np.sum(controls < actual)
                + 0.5 * np.sum(controls == actual)
            )
            / len(controls)
        ),
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    evidence, public_provenance = verify_source_tables()
    panel = read_csv_column(PANEL_CSV, "target_gene")
    official = read_csv_column(GENES_CSV, "gene_name")
    panel_pos = {target: index for index, target in enumerate(panel)}
    official_pos = {gene: index for index, gene in enumerate(official)}
    panel_positions = np.asarray(
        [official_pos[target] for target in panel], dtype=np.int64
    )
    panel_set = set(panel)
    for target in ("PQBP1", "SHPRH", "FZD2", "GLYR1"):
        if target not in panel_pos:
            raise ValueError(f"reported target {target} is not in VCC panel")

    banks = {
        "k562_gwps": Bank.load("k562_gwps", K562_ROOT),
        "orion_hct116": Bank.load("orion_hct116", HCT_ROOT),
        "orion_hek293t": Bank.load("orion_hek293t", HEK_ROOT),
        "perturbmatch": Bank.load("perturbmatch", PM_ROOT),
        "kolf": Bank.load("kolf", KOLF_ROOT),
    }
    cd4_banks = [
        Bank.load(name, CD4_ROOT / name) for name in CD4_NAMES
    ]

    prepared = {
        name: prepare_main(
            bank,
            panel,
            official,
            panel_positions,
            chunk_rows=args.chunk_rows,
        )
        for name, bank in banks.items()
    }
    print(
        f"prepared main banks peak_rss={peak_rss_gib():.6f} GiB",
        flush=True,
    )
    cd4 = prepare_cd4(
        cd4_banks,
        panel,
        official,
        panel_set,
        chunk_rows=args.chunk_rows,
    )
    print(f"prepared CD4 peak_rss={peak_rss_gib():.6f} GiB", flush=True)
    pm_slot, pm_available, pm_shared = prepare_pm_slot(
        prepared["k562_gwps"],
        prepared["perturbmatch"],
        banks["k562_gwps"],
        banks["perturbmatch"],
        official,
        panel_set,
    )
    # The K562 and PerturbMatch full target matrices are no longer needed once
    # their frozen PM50 slot has been formed. Keep only availability and the
    # incumbent mean vector required by the existing missing-target fallback.
    k562_mean = prepared["k562_gwps"].mean_on_axis.copy()
    k562_available = prepared["k562_gwps"].available.copy()
    pm_source_available = prepared["perturbmatch"].available.copy()
    prepared["k562_gwps"].effects = np.empty((0, 0), dtype=np.float64)
    prepared["perturbmatch"].effects = np.empty((0, 0), dtype=np.float64)
    gc.collect()
    enforce_rss("released K562 and PerturbMatch target matrices")

    exp064 = json.loads(EXP064_JSON.read_text())
    preserve = {
        "hct_geometry": float(
            next(
                fold
                for fold in exp064["folds"]
                if fold["held"] == "orion_hct116"
            )["production_preserve"]["B"]["applied_scale"]
        ),
        "hek_geometry": float(
            next(
                fold
                for fold in exp064["folds"]
                if fold["held"] == "orion_hek293t"
            )["production_preserve"]["B"]["applied_scale"]
        ),
    }
    operators = {
        "hct_geometry": load_operator(
            "hct_geometry",
            WHITENER_ROOT / "held_orion_hct116.npz",
            [banks["k562_gwps"], banks["orion_hek293t"]],
            official,
            panel_set,
            preserve["hct_geometry"],
        ),
        "hek_geometry": load_operator(
            "hek_geometry",
            WHITENER_ROOT / "held_orion_hek293t.npz",
            [banks["k562_gwps"], banks["orion_hct116"]],
            official,
            panel_set,
            preserve["hek_geometry"],
        ),
    }
    if not np.array_equal(
        operators["hct_geometry"].wcols, operators["hek_geometry"].wcols
    ):
        raise ValueError("frozen experiment-064 operators do not share an axis")

    cases = build_cases()
    null_targets: dict[str, tuple[str, ...]] = {}
    for offset, (target, target_cases) in enumerate(cases.items()):
        null_targets[target] = choose_null_controls(
            target,
            target_cases,
            panel,
            panel_pos,
            prepared,
            pm_available,
            cd4,
            seed=SEED + offset,
        )

    targets_by_config_mask: dict[
        str, dict[frozenset[str], set[str]]
    ] = {
        "hct_geometry": {frozenset(): set()},
        "hek_geometry": {frozenset(): set()},
    }
    for target, target_cases in cases.items():
        controls = null_targets[target]
        for case in target_cases:
            config = case["config"]
            removed = case["removed"]
            targets_by_config_mask[config].setdefault(removed, set()).update(
                (target, *controls)
            )
            targets_by_config_mask[config][frozenset()].update(
                (target, *controls)
            )

    predictions: dict[
        str, dict[tuple[str, frozenset[str]], np.ndarray]
    ] = {}
    for config, grouped in targets_by_config_mask.items():
        companion_name = (
            "orion_hek293t"
            if config == "hct_geometry"
            else "orion_hct116"
        )
        predictions[config] = prediction_rows_for_config(
            config,
            {
                mask: tuple(sorted(targets))
                for mask, targets in grouped.items()
            },
            panel_pos,
            pm_slot,
            pm_available,
            PreparedMain(
                name="k562_gwps",
                effects=np.empty((0, 0), dtype=np.float64),
                available=k562_available,
                mean_on_axis=k562_mean,
                raw_median_norm=prepared["k562_gwps"].raw_median_norm,
                mapped_gene_count=prepared["k562_gwps"].mapped_gene_count,
                duplicate_symbols=prepared["k562_gwps"].duplicate_symbols,
            ),
            prepared[companion_name],
            cd4,
            operators[config],
            block_rows=args.evec_block_rows,
        )
        print(
            f"transformed {config} peak_rss={peak_rss_gib():.6f} GiB",
            flush=True,
        )

    galleries = {
        truth: tuple(
            target
            for target in panel
            if prepared[truth].available[panel_pos[target]]
        )
        for truth in ("orion_hct116", "orion_hek293t", "kolf")
    }
    retrieval = {
        (truth, config): np.asarray(
            [
                position
                for position in operators[config].wcols
                if official[position] in banks[truth].gene_pos
            ],
            dtype=np.int64,
        )
        for truth in galleries
        for config in operators
    }

    case_rows: list[dict[str, Any]] = []
    for report_target, target_cases in cases.items():
        for evaluated_target in (report_target, *null_targets[report_target]):
            is_reported = evaluated_target == report_target
            target_index = panel_pos[evaluated_target]
            for case in target_cases:
                truth = case["truth"]
                config = case["config"]
                removed = case["removed"]
                gallery = galleries[truth]
                gallery_indices = np.asarray(
                    [panel_pos[target] for target in gallery], dtype=np.int64
                )
                gallery_target_pos = gallery.index(evaluated_target)
                fixed = rank_query(
                    predictions[config][(evaluated_target, frozenset())],
                    prepared[truth],
                    gallery_indices,
                    gallery_target_pos,
                    retrieval[(truth, config)],
                    operators[config].wcols,
                )
                ablated = rank_query(
                    predictions[config][(evaluated_target, removed)],
                    prepared[truth],
                    gallery_indices,
                    gallery_target_pos,
                    retrieval[(truth, config)],
                    operators[config].wcols,
                )
                case_rows.append(
                    {
                        "report_target": report_target,
                        "evaluated_target": evaluated_target,
                        "is_reported": is_reported,
                        "tier": case["tier"],
                        "case_id": case["case_id"],
                        "held_truth": truth,
                        "config": config,
                        "removed_sources": "+".join(sorted(removed)),
                        "gallery_size": len(gallery),
                        "retrieval_genes": len(retrieval[(truth, config)]),
                        "fixed_rank_1based": fixed["rank_1based"],
                        "ablated_rank_1based": ablated["rank_1based"],
                        "rank_delta_ablated_minus_fixed": (
                            ablated["rank_1based"] - fixed["rank_1based"]
                        ),
                        "fixed_pds": fixed["pds"],
                        "ablated_pds": ablated["pds"],
                        "pds_delta": ablated["pds"] - fixed["pds"],
                        "fixed_own_cosine": fixed["own_cosine"],
                        "ablated_own_cosine": ablated["own_cosine"],
                        "fixed_ties_other": fixed["ties_other"],
                        "ablated_ties_other": ablated["ties_other"],
                    }
                )
                enforce_rss(f"scored {report_target} {case['case_id']}")

    target_rows: list[dict[str, Any]] = []
    target_control_deltas: dict[str, dict[str, float]] = {}
    for report_target in cases:
        actual_rows = [
            row
            for row in case_rows
            if row["report_target"] == report_target and row["is_reported"]
        ]
        actual = float(np.mean([row["pds_delta"] for row in actual_rows]))
        controls: dict[str, list[float]] = {}
        for row in case_rows:
            if row["report_target"] != report_target or row["is_reported"]:
                continue
            controls.setdefault(str(row["evaluated_target"]), []).append(
                float(row["pds_delta"])
            )
        control_means = {
            target: float(np.mean(values)) for target, values in controls.items()
        }
        target_control_deltas[report_target] = control_means
        null = summarize_null(
            actual, np.asarray(list(control_means.values()), dtype=np.float64)
        )
        native_rows = [row for row in actual_rows if row["tier"] == "native_064"]
        secondary_rows = [
            row for row in actual_rows if row["tier"] == "secondary_kolf"
        ]
        target_rows.append(
            {
                "target_gene": report_target,
                "cases": len(actual_rows),
                "native_cases": len(native_rows),
                "secondary_cases": len(secondary_rows),
                "mean_pds_delta_all_cases": actual,
                "mean_pds_delta_native": (
                    float(np.mean([row["pds_delta"] for row in native_rows]))
                    if native_rows
                    else ""
                ),
                "mean_pds_delta_secondary": (
                    float(np.mean([row["pds_delta"] for row in secondary_rows]))
                    if secondary_rows
                    else ""
                ),
                "null_controls": null["controls"],
                "null_mean": null["mean"],
                "null_p95": null["p95"],
                "null_p_ge_actual": null["one_sided_p_ge_actual"],
                "actual_percentile": null["actual_percentile_midrank"],
            }
        )

    actual_by_target = {
        row["target_gene"]: float(row["mean_pds_delta_all_cases"])
        for row in target_rows
    }
    overall_actual = float(np.mean(list(actual_by_target.values())))
    rng = np.random.default_rng(SEED)
    null_draws = np.empty(NULL_REPS, dtype=np.float64)
    report_targets = tuple(cases)
    for replicate in range(NULL_REPS):
        values = []
        for target in report_targets:
            candidates = tuple(target_control_deltas[target].values())
            values.append(float(candidates[rng.integers(0, len(candidates))]))
        null_draws[replicate] = float(np.mean(values))
    overall_null = {
        "reps": NULL_REPS,
        "seed": SEED,
        **summarize_null(overall_actual, null_draws),
    }

    native_actual = [
        float(row["pds_delta"])
        for row in case_rows
        if row["is_reported"] and row["tier"] == "native_064"
    ]

    exp064_rows = list(csv.DictReader(EXP064_PDS.open(newline="")))
    expected_064 = {
        (row["held"], row["target_gene"]): float(row["pds_D"])
        for row in exp064_rows
    }
    reconstruction_rows: list[dict[str, Any]] = []
    gallery_064 = tuple(exp064["axis"]["target_order"])
    gallery_064_indices = np.asarray(
        [panel_pos[target] for target in gallery_064], dtype=np.int64
    )
    for truth, config in (
        ("orion_hct116", "hct_geometry"),
        ("orion_hek293t", "hek_geometry"),
    ):
        for target in sorted(
            target
            for target, mask in predictions[config]
            if not mask and target in set(gallery_064)
        ):
            observed = rank_query(
                predictions[config][(target, frozenset())],
                prepared[truth],
                gallery_064_indices,
                gallery_064.index(target),
                retrieval[(truth, config)],
                operators[config].wcols,
            )
            expected = expected_064[(truth, target)]
            reconstruction_rows.append(
                {
                    "held": truth,
                    "target_gene": target,
                    "expected_pds_D": expected,
                    "reconstructed_pds_D": observed["pds"],
                    "absolute_error": abs(observed["pds"] - expected),
                }
            )
    max_reconstruction_error = max(
        row["absolute_error"] for row in reconstruction_rows
    )

    target_positive = sum(
        float(row["mean_pds_delta_all_cases"]) > 0 for row in target_rows
    )
    target_negative = sum(
        float(row["mean_pds_delta_all_cases"]) < 0 for row in target_rows
    )
    decision_checks = {
        "native_mean_positive": float(np.mean(native_actual)) > 0,
        "at_least_three_targets_positive": target_positive >= 3,
        "no_target_negative": target_negative == 0,
        "overall_exceeds_null_p95": overall_actual > overall_null["p95"],
        "exact_064_rank_reconstruction": max_reconstruction_error == 0.0,
    }
    decision = "SHIP" if all(decision_checks.values()) else "REJECT"

    source_by_target = {
        "PQBP1": "TIMM17B",
        "SHPRH": "CENPP",
        "FZD2": "P2RX7",
        "GLYR1": "UBN1",
    }
    result = {
        "experiment": "068-seed-offtarget-ablation",
        "status": "completed",
        "decision": decision,
        "decision_checks": decision_checks,
        "protocol": {
            "comparator": "frozen experiment-064 arm D",
            "reported_pairs": source_by_target,
            "primary_truths": ["orion_hct116", "orion_hek293t"],
            "secondary_truth": "kolf",
            "secondary_truth_role": (
                "independent held truth needed for GLYR1; never a source or "
                "operator input"
            ),
            "null_matching": (
                "same cases, source-ablation masks, and exact source/truth "
                "availability pattern; up to 64 deterministic controls"
            ),
            "panel_coordinates_excluded": True,
            "operator_refit": False,
            "hyperparameter_tuning": False,
            "count_emission": False,
            "cell_eval2": False,
            "cloud": False,
            "submission": False,
        },
        "source_evidence": evidence,
        "public_provenance": public_provenance,
        "axes": {
            "official_genes": len(official),
            "panel_targets": len(panel),
            "pm_shared_non_panel_genes": len(pm_shared),
            "operator_genes": {
                config: len(operator.wcols)
                for config, operator in operators.items()
            },
            "held_galleries": {
                truth: len(gallery) for truth, gallery in galleries.items()
            },
            "retrieval_genes": {
                f"{truth}__{config}": len(columns)
                for (truth, config), columns in retrieval.items()
            },
        },
        "target_results": target_rows,
        "actual_case_results": [
            row for row in case_rows if row["is_reported"]
        ],
        "native_064_summary": {
            "cases": len(native_actual),
            "mean_pds_delta": float(np.mean(native_actual)),
        },
        "four_target_summary": {
            "mean_pds_delta": overall_actual,
            "positive_targets": target_positive,
            "negative_targets": target_negative,
            "matched_random_target_null": overall_null,
        },
        "null_controls": {
            target: list(controls) for target, controls in null_targets.items()
        },
        "reconstruction": {
            "rows": len(reconstruction_rows),
            "max_absolute_pds_error": max_reconstruction_error,
            "exact": max_reconstruction_error == 0.0,
        },
        "leakage_audit": {
            "challenge_truth_used": False,
            "held_effect_used_for_rule_or_control_selection": False,
            "held_effect_used_only_for_final_rank_evaluation": True,
            "reported_rules_fixed_from_public_source_tables": True,
            "target_identity_used_only_for_predeclared_external_qc_rules": True,
            "null_controls_selected_without_outcomes": True,
            "operator_candidate_independent_and_frozen": True,
            "all_300_panel_coordinates_excluded": True,
            "direct_offtarget_effect_injected": False,
            "direct_offtarget_coordinate_zeroed_ad_hoc": False,
            "kolf_used_as_source": False,
            "manual_rule_scope": (
                "four source-row exclusions only; potential challenge-rule "
                "eligibility should be reviewed before any production use"
            ),
        },
        "resource_audit": {
            "processes_spawned_by_runner": 0,
            "blas_threads": 1,
            "hard_cap_gib": MAX_RSS_GIB,
            "internal_stop_threshold_gib": STOP_RSS_GIB,
            "peak_rss_gib": peak_rss_gib(),
            "max_download_bytes_per_object": MAX_DOWNLOAD_BYTES,
            "largest_download_bytes": max(
                meta["compressed_bytes"]
                for meta in public_provenance["tables"].values()
            ),
        },
        "input_hashes": {
            "run_py": sha256_file(Path(__file__)),
            "panel_csv": sha256_file(PANEL_CSV),
            "genes_csv": sha256_file(GENES_CSV),
            "exp064_results_json": sha256_file(EXP064_JSON),
            "exp064_per_target_pds_csv": sha256_file(EXP064_PDS),
            "whiteners": {
                config: sha256_file(operator.path)
                for config, operator in operators.items()
            },
            "banks": {
                name: {
                    "delta_npy": sha256_file(bank.delta_path),
                    "genes_txt": sha256_file(bank.path / "genes.txt"),
                    "perts_txt": sha256_file(bank.path / "perts.txt"),
                    "meta_json": sha256_file(bank.path / "meta.json"),
                }
                for name, bank in {
                    **banks,
                    **{bank.name: bank for bank in cd4_banks},
                }.items()
            },
        },
    }

    write_csv(HERE / "source_evidence.csv", evidence)
    write_csv(HERE / "per_case_results.csv", case_rows)
    write_csv(HERE / "per_target_results.csv", target_rows)
    write_csv(
        HERE / "null_controls.csv",
        [
            {
                "report_target": report_target,
                "control_target": control_target,
                "mean_pds_delta": delta,
            }
            for report_target, controls in target_control_deltas.items()
            for control_target, delta in controls.items()
        ],
    )
    write_csv(HERE / "reconstruction.csv", reconstruction_rows)
    result["resource_audit"]["peak_rss_gib"] = peak_rss_gib()
    if result["resource_audit"]["peak_rss_gib"] >= MAX_RSS_GIB:
        raise MemoryError("hard RSS cap was exceeded")
    (HERE / "metrics.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {
                "decision": decision,
                "native_064_summary": result["native_064_summary"],
                "four_target_summary": result["four_target_summary"],
                "peak_rss_gib": result["resource_audit"]["peak_rss_gib"],
                "reconstruction": result["reconstruction"],
            },
            indent=2,
        ),
        flush=True,
    )
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chunk-rows", type=int, default=4)
    parser.add_argument("--evec-block-rows", type=int, default=128)
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing unguarded real-data execution")
    if args.chunk_rows < 1 or args.evec_block_rows < 1:
        raise SystemExit("block sizes must be positive")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

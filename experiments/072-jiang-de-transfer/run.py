#!/usr/bin/env python3
"""Audit and evaluate the Jiang 2025 pathway DE archive without extracting it.

The archive contains one whitespace-delimited DE table per
``(perturbation target, stimulation pathway)``.  Each table reports effects for
six cell lines.  Pathway is an experimental batch, so Jiang effects are
centered within pathway using non-VCC targets only.  Non-panel targets use a
leave-one-target-out pathway center.

The transfer screen is deliberately separated into two evidential layers:

1. A development-only screen on non-panel perturbations shared with the three
   established donor banks.  It may select one representation and blend.
2. A held VCC-panel evaluation, allowed only when the archive contains enough
   panel targets for a stable gallery.

The shipped Jiang archive has nine panel targets, below the predeclared minimum
of 30.  Therefore the real run records a negative metadata gate and never reads
held-donor effects for those nine targets.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import resource
import sys
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

EXPECTED_MD5 = "f077cba680a1affc599f5153d99b0e45"
EXPECTED_DATA_FILES = 271
EXPECTED_UNIQUE_TARGETS = 218
EXPECTED_PATHWAYS = ("IFNB", "IFNG", "INS", "TGFB1", "TNFA")
CELL_LINES = ("A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7")
MAIN_BANK_PATHS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
HELD_BANKS = tuple(MAIN_BANK_PATHS)
REPRESENTATIONS = ("raw", "gene_norm")
DEFAULT_FRACTIONS = (0.02, 0.05, 0.10, 0.20, 0.35)


def peak_rss_gib() -> float:
    """Return peak process RSS in GiB on macOS or Linux."""
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def enforce_limit(limit_gib: float, stage: str) -> None:
    if limit_gib > 0 and peak_rss_gib() >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak_rss_gib():.3f} GiB >= {limit_gib:.3f} GiB"
        )


def digest_file(path: Path, algorithm: str, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256(
        "".join(f"{value}\n" for value in values).encode()
    ).hexdigest()


def read_csv_axis(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: missing column {column!r}")
        values = tuple(value for row in reader if (value := str(row[column]).strip()))
    if not values or len(values) != len(set(values)):
        raise ValueError(f"{path}: empty or duplicate {column!r} values")
    return values


def read_lines(path: Path, *, require_unique: bool = True) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if not values:
        raise ValueError(f"{path}: empty labels")
    if require_unique and len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    fields = list(rows[0])
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


@dataclass(frozen=True)
class Member:
    name: str
    target: str
    pathway: str
    compressed_size: int
    uncompressed_size: int
    crc: int


def parse_member(info: zipfile.ZipInfo) -> Member | None:
    """Parse one scientific DE table name; ignore directories and Finder files."""
    if info.is_dir() or not info.filename.endswith("_pathway_DE_results.txt"):
        return None
    parts = Path(info.filename).parts
    if len(parts) != 3 or not parts[1].startswith("Parse_"):
        raise ValueError(f"unexpected DE member path {info.filename!r}")
    pathway = parts[1].removeprefix("Parse_")
    suffix = f"_{pathway}_pathway_DE_results.txt"
    basename = parts[2]
    if not basename.endswith(suffix):
        raise ValueError(f"member/file pathway mismatch: {info.filename!r}")
    target = basename[: -len(suffix)]
    if not target:
        raise ValueError(f"empty target in {info.filename!r}")
    return Member(
        name=info.filename,
        target=target,
        pathway=pathway,
        compressed_size=int(info.compress_size),
        uncompressed_size=int(info.file_size),
        crc=int(getattr(info, "CRC", 0)),
    )


def expected_header() -> tuple[str, ...]:
    return (
        "gene_ID",
        *(f"log2FC_{cell}" for cell in CELL_LINES),
        *(f"beta_cell_type{cell}" for cell in CELL_LINES),
        *(f"p_cell_type{cell}" for cell in CELL_LINES),
    )


def parse_effect(value: str, *, member: str, line_number: int) -> float:
    if value == "NA":
        return math.nan
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(f"{member}:{line_number}: invalid effect {value!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{member}:{line_number}: non-finite effect")
    return parsed


def iter_effect_rows(
    archive: zipfile.ZipFile,
    member: Member,
) -> Iterable[tuple[str, np.ndarray]]:
    """Yield gene and six log2FC values from one member."""
    with archive.open(member.name, "r") as binary:
        with io.TextIOWrapper(binary, encoding="utf-8", newline="") as handle:
            header_line = handle.readline()
            if not header_line:
                raise ValueError(f"{member.name}: empty file")
            header = tuple(header_line.strip().split())
            if header != expected_header():
                raise ValueError(f"{member.name}: unexpected header {header!r}")
            for line_number, line in enumerate(handle, start=2):
                fields = line.strip().split()
                if not fields:
                    continue
                if len(fields) != len(header):
                    raise ValueError(
                        f"{member.name}:{line_number}: expected "
                        f"{len(header)} fields, got {len(fields)}"
                    )
                values = np.fromiter(
                    (
                        parse_effect(
                            fields[index],
                            member=member.name,
                            line_number=line_number,
                        )
                        for index in range(1, 1 + len(CELL_LINES))
                    ),
                    dtype=np.float64,
                    count=len(CELL_LINES),
                )
                yield fields[0], values


@dataclass(frozen=True)
class AuditState:
    members: tuple[Member, ...]
    targets: tuple[str, ...]
    target_pathways: dict[str, tuple[str, ...]]
    file_rows: tuple[dict[str, object], ...]
    panel_targets: tuple[str, ...]
    proxy_targets: tuple[str, ...]
    shared_panel_targets: tuple[str, ...]
    gene_seen: np.ndarray
    pathway_sum: np.ndarray
    pathway_count: np.ndarray
    archive_sha256: str
    archive_md5: str


def audit_archive(
    *,
    archive_path: Path,
    official_genes: tuple[str, ...],
    panel: tuple[str, ...],
    bank_perts: dict[str, set[str]],
    max_rss_gib: float,
) -> AuditState:
    """Validate all DE members and accumulate non-panel pathway centers."""
    if not archive_path.is_file():
        raise FileNotFoundError(archive_path)
    archive_md5 = digest_file(archive_path, "md5")
    if archive_md5 != EXPECTED_MD5:
        raise ValueError(f"archive MD5 {archive_md5} != expected {EXPECTED_MD5}")
    archive_sha256 = digest_file(archive_path, "sha256")
    panel_set = set(panel)
    gene_pos = {gene: index for index, gene in enumerate(official_genes)}
    pathway_pos = {pathway: index for index, pathway in enumerate(EXPECTED_PATHWAYS)}
    pathway_sum = np.zeros(
        (len(EXPECTED_PATHWAYS), len(CELL_LINES), len(official_genes)),
        dtype=np.float64,
    )
    pathway_count = np.zeros(pathway_sum.shape, dtype=np.uint16)
    gene_seen = np.zeros(len(official_genes), dtype=np.bool_)
    file_rows: list[dict[str, object]] = []

    with zipfile.ZipFile(archive_path) as archive:
        parsed = [
            member
            for info in archive.infolist()
            if (member := parse_member(info)) is not None
        ]
        members = tuple(sorted(parsed, key=lambda item: item.name))
        if len(members) != EXPECTED_DATA_FILES:
            raise ValueError(
                f"found {len(members)} DE files; expected {EXPECTED_DATA_FILES}"
            )
        if any(member.pathway not in pathway_pos for member in members):
            raise ValueError("archive contains an unexpected pathway")

        target_pathways_build: dict[str, set[str]] = defaultdict(set)
        for member_index, member in enumerate(members, start=1):
            target_pathways_build[member.target].add(member.pathway)
            seen_in_file: set[str] = set()
            mapped = 0
            duplicate_rows = 0
            finite_by_line = np.zeros(len(CELL_LINES), dtype=np.int64)
            common_finite = 0
            pathway_index = pathway_pos[member.pathway]
            is_panel = member.target in panel_set

            for gene, values in iter_effect_rows(archive, member):
                if gene in seen_in_file:
                    duplicate_rows += 1
                seen_in_file.add(gene)
                column = gene_pos.get(gene)
                if column is None:
                    continue
                mapped += 1
                gene_seen[column] = True
                finite = np.isfinite(values)
                finite_by_line += finite
                common_finite += int(finite.all())
                if not is_panel:
                    pathway_sum[pathway_index, finite, column] += values[finite]
                    pathway_count[pathway_index, finite, column] += 1

            row: dict[str, object] = {
                "archive_member": member.name,
                "target_gene": member.target,
                "pathway": member.pathway,
                "is_vcc_panel": is_panel,
                "rows": len(seen_in_file) + duplicate_rows,
                "unique_genes": len(seen_in_file),
                "duplicate_gene_rows": duplicate_rows,
                "mapped_vcc_genes": mapped,
                "common_finite_vcc_genes": common_finite,
                "compressed_bytes": member.compressed_size,
                "uncompressed_bytes": member.uncompressed_size,
                "crc32": f"{member.crc:08x}",
            }
            row.update(
                {
                    f"finite_log2fc_{cell}": int(finite_by_line[index])
                    for index, cell in enumerate(CELL_LINES)
                }
            )
            file_rows.append(row)
            if member_index % 25 == 0 or member_index == len(members):
                print(
                    f"[audit] {member_index}/{len(members)} members; "
                    f"peak_rss={peak_rss_gib():.3f} GiB",
                    flush=True,
                )
                enforce_limit(max_rss_gib, f"audit member {member_index}")

    targets = tuple(sorted(target_pathways_build))
    if len(targets) != EXPECTED_UNIQUE_TARGETS:
        raise ValueError(
            f"found {len(targets)} unique targets; expected {EXPECTED_UNIQUE_TARGETS}"
        )
    target_pathways = {
        target: tuple(sorted(pathways))
        for target, pathways in sorted(target_pathways_build.items())
    }
    panel_targets = tuple(target for target in panel if target in target_pathways)
    non_panel = set(targets) - panel_set
    proxy_targets = tuple(
        sorted(
            target
            for target in non_panel
            if all(target in values for values in bank_perts.values())
        )
    )
    shared_panel_targets = tuple(
        target
        for target in panel_targets
        if all(target in values for values in bank_perts.values())
    )
    return AuditState(
        members=members,
        targets=targets,
        target_pathways=target_pathways,
        file_rows=tuple(file_rows),
        panel_targets=panel_targets,
        proxy_targets=proxy_targets,
        shared_panel_targets=shared_panel_targets,
        gene_seen=gene_seen,
        pathway_sum=pathway_sum,
        pathway_count=pathway_count,
        archive_sha256=archive_sha256,
        archive_md5=archive_md5,
    )


@dataclass(frozen=True)
class JiangEffects:
    targets: tuple[str, ...]
    target_pos: dict[str, int]
    genes: tuple[str, ...]
    effects: np.ndarray
    observed: np.ndarray
    pathway_observations: np.ndarray


def build_centered_effects(
    *,
    archive_path: Path,
    audit: AuditState,
    official_genes: tuple[str, ...],
    panel: tuple[str, ...],
    max_rss_gib: float,
) -> JiangEffects:
    """Build pathway-centered effects with non-panel leave-one-out centering."""
    panel_set = set(panel)
    gene_pos = {gene: index for index, gene in enumerate(official_genes)}
    pathway_pos = {pathway: index for index, pathway in enumerate(EXPECTED_PATHWAYS)}
    target_pos = {target: index for index, target in enumerate(audit.targets)}
    shape = (len(audit.targets), len(CELL_LINES), len(official_genes))
    effect_sum = np.zeros(shape, dtype=np.float32)
    observations = np.zeros(shape, dtype=np.uint8)

    with zipfile.ZipFile(archive_path) as archive:
        for member_index, member in enumerate(audit.members, start=1):
            target_index = target_pos[member.target]
            pathway_index = pathway_pos[member.pathway]
            is_panel = member.target in panel_set
            for gene, values in iter_effect_rows(archive, member):
                column = gene_pos.get(gene)
                if column is None:
                    continue
                finite = np.isfinite(values)
                for line_index in np.flatnonzero(finite):
                    count = int(audit.pathway_count[pathway_index, line_index, column])
                    value = float(values[line_index])
                    if is_panel:
                        if count < 1:
                            continue
                        center = (
                            audit.pathway_sum[pathway_index, line_index, column] / count
                        )
                    else:
                        if count < 2:
                            continue
                        center = (
                            audit.pathway_sum[pathway_index, line_index, column] - value
                        ) / (count - 1)
                    effect_sum[target_index, line_index, column] += value - center
                    observations[target_index, line_index, column] += 1
            if member_index % 25 == 0 or member_index == len(audit.members):
                print(
                    f"[center] {member_index}/{len(audit.members)} members; "
                    f"peak_rss={peak_rss_gib():.3f} GiB",
                    flush=True,
                )
                enforce_limit(max_rss_gib, f"center member {member_index}")

    observed = observations > 0
    np.divide(
        effect_sum,
        observations,
        out=effect_sum,
        where=observed,
        casting="unsafe",
    )
    if not np.isfinite(effect_sum).all():
        raise ValueError("centered Jiang effects contain non-finite values")
    return JiangEffects(
        targets=audit.targets,
        target_pos=target_pos,
        genes=official_genes,
        effects=effect_sum,
        observed=observed,
        pathway_observations=observations,
    )


@dataclass(frozen=True)
class NpyLayout:
    path: Path
    shape: tuple[int, int]
    dtype: np.dtype[Any]
    data_offset: int
    row_bytes: int

    @classmethod
    def inspect(cls, path: Path) -> "NpyLayout":
        with path.open("rb") as handle:
            version = np.lib.format.read_magic(handle)
            if version == (1, 0):
                shape, fortran, dtype = np.lib.format.read_array_header_1_0(handle)
            elif version in ((2, 0), (3, 0)):
                shape, fortran, dtype = np.lib.format.read_array_header_2_0(handle)
            else:
                raise ValueError(f"{path}: unsupported NPY version {version}")
            offset = handle.tell()
        if fortran or len(shape) != 2:
            raise ValueError(f"{path}: expected C-order matrix")
        dtype = np.dtype(dtype)
        expected = offset + int(np.prod(shape)) * dtype.itemsize
        if path.stat().st_size != expected:
            raise ValueError(f"{path}: malformed NPY byte size")
        return cls(
            path=path,
            shape=(int(shape[0]), int(shape[1])),
            dtype=dtype,
            data_offset=offset,
            row_bytes=int(shape[1]) * dtype.itemsize,
        )

    def read_rows(self, rows: np.ndarray) -> np.ndarray:
        output = np.empty((len(rows), self.shape[1]), dtype=self.dtype)
        descriptor = os.open(self.path, os.O_RDONLY)
        try:
            for output_index, row in enumerate(rows):
                offset = self.data_offset + int(row) * self.row_bytes
                payload = os.pread(descriptor, self.row_bytes, offset)
                if len(payload) != self.row_bytes:
                    raise EOFError(f"{self.path}: short row read at {row}")
                output[output_index] = np.frombuffer(payload, dtype=self.dtype)
        finally:
            os.close(descriptor)
        return output


@dataclass(frozen=True)
class Bank:
    name: str
    genes: tuple[str, ...]
    perts: tuple[str, ...]
    gene_pos: dict[str, int]
    pert_pos: dict[str, int]
    layout: NpyLayout

    @classmethod
    def load(cls, name: str, path: Path) -> "Bank":
        genes = read_lines(path / "genes.txt", require_unique=False)
        perts = read_lines(path / "perts.txt")
        layout = NpyLayout.inspect(path / "delta.npy")
        if layout.shape != (len(perts), len(genes)):
            raise ValueError(f"{name}: donor shape mismatch")
        b_ctrl = np.load(path / "b_ctrl.npy", mmap_mode="r")
        gene_pos: dict[str, int] = {}
        for index, gene in enumerate(genes):
            previous = gene_pos.get(gene)
            if previous is None or float(b_ctrl[index]) > float(b_ctrl[previous]):
                gene_pos[gene] = index
        return cls(
            name=name,
            genes=genes,
            perts=perts,
            gene_pos=gene_pos,
            pert_pos={target: index for index, target in enumerate(perts)},
            layout=layout,
        )


def build_feature_axis(
    *,
    official_genes: tuple[str, ...],
    panel: set[str],
    proxy_targets: set[str],
    banks: dict[str, Bank],
    jiang_gene_seen: np.ndarray,
) -> tuple[str, ...]:
    """Construct a common axis without panel or proxy-target coordinates."""
    forbidden = panel | proxy_targets
    return tuple(
        gene
        for index, gene in enumerate(official_genes)
        if jiang_gene_seen[index]
        and gene not in forbidden
        and all(gene in bank.gene_pos for bank in banks.values())
    )


def prepare_bank_effects(
    *,
    bank: Bank,
    targets: tuple[str, ...],
    features: tuple[str, ...],
    excluded_center_targets: set[str],
    chunk_rows: int,
    max_rss_gib: float,
) -> np.ndarray:
    """Read selected rows and center on disjoint background perturbations."""
    columns = np.asarray([bank.gene_pos[gene] for gene in features])
    background = np.asarray(
        [
            index
            for index, target in enumerate(bank.perts)
            if target not in excluded_center_targets
        ],
        dtype=np.int64,
    )
    if len(background) < 100:
        raise ValueError(f"{bank.name}: insufficient background rows")
    center = np.zeros(len(features), dtype=np.float64)
    for start in range(0, len(background), chunk_rows):
        rows = background[start : start + chunk_rows]
        # Avoid mmap here. Sequentially touching the multi-GiB Orion arrays can
        # leave file-backed pages resident on macOS even though each numerical
        # chunk is tiny, causing a misleading and eventually unsafe RSS climb.
        chunk = bank.layout.read_rows(rows)
        center += np.asarray(chunk[:, columns], dtype=np.float64).sum(axis=0)
        enforce_limit(max_rss_gib, f"{bank.name} center")
    center /= len(background)
    target_rows = np.asarray([bank.pert_pos[target] for target in targets])
    selected = np.asarray(
        bank.layout.read_rows(target_rows)[:, columns],
        dtype=np.float64,
    )
    selected -= center
    return selected


def row_normalize(matrix: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    return matrix / np.maximum(norms, 1e-12)


def pds_per_query(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    similarity = row_normalize(prediction) @ row_normalize(truth).T
    if similarity.shape[0] < 2 or similarity.shape[0] != similarity.shape[1]:
        raise ValueError("PDS requires a square gallery of at least two targets")
    own = np.diag(similarity)[:, None]
    rank = (similarity > own).sum(axis=1).astype(np.float64)
    rank += 0.5 * ((similarity == own).sum(axis=1) - 1)
    return 1.0 - rank / (similarity.shape[1] - 1)


def norm_matched_mean(parts: Sequence[np.ndarray]) -> np.ndarray:
    norms = [float(np.median(np.linalg.norm(part, axis=1))) or 1.0 for part in parts]
    reference = float(np.median(norms)) or 1.0
    return sum(part * (reference / norm) for part, norm in zip(parts, norms)) / len(
        parts
    )


def jiang_line_indices(held_bank: str) -> tuple[int, ...]:
    """Prevent same-cell-line leakage when K562 is the held truth."""
    return tuple(
        index
        for index, line in enumerate(CELL_LINES)
        if not (held_bank == "k562_gwps" and line == "K562")
    )


def jiang_consensus(
    *,
    jiang: JiangEffects,
    targets: tuple[str, ...],
    feature_columns: np.ndarray,
    line_indices: tuple[int, ...],
) -> tuple[np.ndarray, np.ndarray]:
    target_rows = np.asarray([jiang.target_pos[target] for target in targets])
    values = jiang.effects[
        np.ix_(target_rows, np.asarray(line_indices), feature_columns)
    ].astype(np.float64)
    observed = jiang.observed[
        np.ix_(target_rows, np.asarray(line_indices), feature_columns)
    ]
    counts = observed.sum(axis=1)
    output = np.zeros((len(targets), len(feature_columns)), dtype=np.float64)
    np.divide(
        (values * observed).sum(axis=1),
        counts,
        out=output,
        where=counts > 0,
    )
    return output, counts > 0


def gene_normalize(
    values: np.ndarray, observed: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    counts = observed.sum(axis=0)
    scale = np.zeros(values.shape[1], dtype=np.float64)
    np.divide(
        (np.abs(values) * observed).sum(axis=0),
        counts,
        out=scale,
        where=counts > 0,
    )
    live = scale[scale > 0]
    floor = float(np.percentile(live, 25)) if live.size else 1.0
    scale = np.maximum(scale, floor)
    output = np.zeros_like(values)
    np.divide(values, scale, out=output, where=observed)
    return output, scale


def blend_candidate(
    *,
    baseline: np.ndarray,
    jiang: np.ndarray,
    observed: np.ndarray,
    fraction: float,
    frozen_scale: float | None = None,
) -> tuple[np.ndarray, float]:
    active_rows = observed.any(axis=1)
    if not active_rows.any():
        raise ValueError("Jiang candidate has no active targets")
    baseline_masked = np.where(observed[active_rows], baseline[active_rows], 0.0)
    jiang_masked = np.where(observed[active_rows], jiang[active_rows], 0.0)
    if frozen_scale is None:
        baseline_norm = float(np.median(np.linalg.norm(baseline_masked, axis=1)))
        jiang_norm = float(np.median(np.linalg.norm(jiang_masked, axis=1)))
        if jiang_norm <= 0:
            raise ValueError("Jiang candidate median norm is zero")
        scale = baseline_norm / jiang_norm
    else:
        scale = frozen_scale
    candidate = baseline.copy()
    candidate[observed] = (1.0 - fraction) * baseline[
        observed
    ] + fraction * scale * jiang[observed]
    return candidate, scale


def stable_fold(target: str, folds: int = 5) -> int:
    digest = hashlib.sha256(target.encode()).digest()
    return int.from_bytes(digest[:4], "big") % folds


def bootstrap_interval(
    deltas: np.ndarray, *, repetitions: int, seed: int
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    target_means = deltas.mean(axis=0)
    estimates = np.empty(repetitions, dtype=np.float64)
    for index in range(repetitions):
        sample = rng.integers(0, len(target_means), len(target_means))
        estimates[index] = target_means[sample].mean()
    low, high = np.quantile(estimates, [0.025, 0.975])
    return float(low), float(high)


def run_proxy_development(
    *,
    jiang: JiangEffects,
    audit: AuditState,
    official_genes: tuple[str, ...],
    panel: tuple[str, ...],
    banks: dict[str, Bank],
    fractions: tuple[float, ...],
    chunk_rows: int,
    max_rss_gib: float,
) -> dict[str, Any]:
    targets = audit.proxy_targets
    feature_axis = build_feature_axis(
        official_genes=official_genes,
        panel=set(panel),
        proxy_targets=set(targets),
        banks=banks,
        jiang_gene_seen=audit.gene_seen,
    )
    if len(targets) < 100:
        return {
            "status": "proxy_gate_failed",
            "reason": "fewer than 100 non-panel targets shared by all banks",
            "targets": len(targets),
            "features": len(feature_axis),
        }
    if len(feature_axis) < 5_000:
        return {
            "status": "proxy_gate_failed",
            "reason": "fewer than 5000 leakage-safe shared features",
            "targets": len(targets),
            "features": len(feature_axis),
        }

    excluded = set(panel) | set(targets)
    effects = {
        name: prepare_bank_effects(
            bank=bank,
            targets=targets,
            features=feature_axis,
            excluded_center_targets=excluded,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        for name, bank in banks.items()
    }
    official_pos = {gene: index for index, gene in enumerate(official_genes)}
    feature_columns = np.asarray([official_pos[gene] for gene in feature_axis])
    rows: list[dict[str, Any]] = []
    per_target: list[dict[str, Any]] = []
    candidate_deltas: dict[tuple[str, float], list[np.ndarray]] = defaultdict(list)

    for held in HELD_BANKS:
        source_names = tuple(name for name in HELD_BANKS if name != held)
        baseline = norm_matched_mean([effects[name] for name in source_names])
        truth = effects[held]
        baseline_pds = pds_per_query(baseline, truth)
        raw_jiang, observed = jiang_consensus(
            jiang=jiang,
            targets=targets,
            feature_columns=feature_columns,
            line_indices=jiang_line_indices(held),
        )
        representations: dict[str, np.ndarray] = {"raw": raw_jiang}
        representations["gene_norm"], _ = gene_normalize(raw_jiang, observed)

        for representation, values in representations.items():
            for fraction in fractions:
                candidate, scale = blend_candidate(
                    baseline=baseline,
                    jiang=values,
                    observed=observed,
                    fraction=fraction,
                )
                candidate_pds = pds_per_query(candidate, truth)
                delta = candidate_pds - baseline_pds
                candidate_deltas[(representation, fraction)].append(delta)
                fold_deltas = [
                    float(
                        delta[
                            np.asarray(
                                [stable_fold(target) == fold for target in targets]
                            )
                        ].mean()
                    )
                    for fold in range(5)
                ]
                rows.append(
                    {
                        "representation": representation,
                        "fraction": fraction,
                        "held_bank": held,
                        "baseline_pds": float(baseline_pds.mean()),
                        "candidate_pds": float(candidate_pds.mean()),
                        "delta_pds": float(delta.mean()),
                        "active_targets": int(observed.any(axis=1).sum()),
                        "available_coordinates": int(observed.sum()),
                        "jiang_scale": scale,
                        **{
                            f"target_fold_{fold}_delta": value
                            for fold, value in enumerate(fold_deltas)
                        },
                    }
                )
                for index, target in enumerate(targets):
                    per_target.append(
                        {
                            "representation": representation,
                            "fraction": fraction,
                            "held_bank": held,
                            "target_gene": target,
                            "target_fold": stable_fold(target),
                            "baseline_pds": float(baseline_pds[index]),
                            "candidate_pds": float(candidate_pds[index]),
                            "delta_pds": float(delta[index]),
                            "jiang_available_genes": int(observed[index].sum()),
                        }
                    )
        enforce_limit(max_rss_gib, f"proxy held {held}")

    summaries: list[dict[str, Any]] = []
    for (representation, fraction), fold_values in candidate_deltas.items():
        matrix = np.stack(fold_values)
        held_means = matrix.mean(axis=1)
        low, high = bootstrap_interval(matrix, repetitions=10_000, seed=7_202_026)
        summaries.append(
            {
                "representation": representation,
                "fraction": fraction,
                "macro_delta": float(matrix.mean()),
                "worst_held_delta": float(held_means.min()),
                "positive_held_folds": int((held_means >= 0).sum()),
                "bootstrap_95_low": low,
                "bootstrap_95_high": high,
            }
        )
    selected = max(
        summaries,
        key=lambda row: (
            row["worst_held_delta"],
            row["macro_delta"],
            -row["fraction"],
            row["representation"] == "raw",
        ),
    )
    selected = {
        **selected,
        "selection_rule": (
            "maximize worst held-donor delta, then macro delta, then smaller "
            "fraction, then raw representation"
        ),
        "passes_proxy_gate": bool(
            selected["positive_held_folds"] == len(HELD_BANKS)
            and selected["macro_delta"] >= 0.003
            and selected["bootstrap_95_low"] > 0
        ),
    }
    return {
        "status": "development_only",
        "targets": len(targets),
        "features": len(feature_axis),
        "feature_axis_sha256": sha256_lines(feature_axis),
        "removed_coordinate_count": len(set(panel) | set(targets)),
        "held_k562_jiang_lines": [
            CELL_LINES[index] for index in jiang_line_indices("k562_gwps")
        ],
        "held_k562_omits_jiang_k562": (
            "K562"
            not in [CELL_LINES[index] for index in jiang_line_indices("k562_gwps")]
        ),
        "rows": rows,
        "summaries": summaries,
        "selected": selected,
        "per_target": per_target,
    }


def save_compact(
    output_dir: Path,
    jiang: JiangEffects,
    audit: AuditState,
) -> dict[str, str]:
    compact = output_dir / "compact"
    compact.mkdir(parents=True, exist_ok=True)
    np.save(compact / "effects.npy", jiang.effects, allow_pickle=False)
    np.save(compact / "observed.npy", jiang.observed, allow_pickle=False)
    np.save(
        compact / "pathway_observations.npy",
        jiang.pathway_observations,
        allow_pickle=False,
    )
    (compact / "targets.txt").write_text(
        "".join(f"{target}\n" for target in jiang.targets)
    )
    (compact / "genes.txt").write_text("".join(f"{gene}\n" for gene in jiang.genes))
    metadata = {
        "shape": list(jiang.effects.shape),
        "cell_lines": list(CELL_LINES),
        "target_axis_sha256": sha256_lines(jiang.targets),
        "gene_axis_sha256": sha256_lines(jiang.genes),
        "centering": (
            "within pathway over non-panel targets; non-panel rows use "
            "leave-one-target-out centers; repeated target pathways averaged"
        ),
        "archive_md5": audit.archive_md5,
        "artifacts": {},
    }
    for name in (
        "effects.npy",
        "observed.npy",
        "pathway_observations.npy",
        "targets.txt",
        "genes.txt",
    ):
        metadata["artifacts"][name] = digest_file(compact / name, "sha256")
    write_json(compact / "meta.json", metadata)
    return {
        name: digest_file(compact / name, "sha256")
        for name in (*metadata["artifacts"], "meta.json")
    }


def build_output_rows(
    audit: AuditState,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    target_rows = [
        {
            "target_gene": target,
            "is_vcc_panel": target in set(audit.panel_targets),
            "n_pathways": len(pathways),
            "pathways": ",".join(pathways),
            "shared_by_all_established_banks": all(
                target in set(read_lines(path / "perts.txt"))
                for path in MAIN_BANK_PATHS.values()
            ),
        }
        for target, pathways in audit.target_pathways.items()
    ]
    panel_files = [row for row in audit.file_rows if bool(row["is_vcc_panel"])]
    return target_rows, panel_files


def run(args: argparse.Namespace) -> dict[str, Any]:
    panel = read_csv_axis(args.panel_csv, "target_gene")
    genes = read_csv_axis(args.genes_csv, "gene_name")
    banks = {name: Bank.load(name, path) for name, path in MAIN_BANK_PATHS.items()}
    bank_perts = {name: set(bank.perts) for name, bank in banks.items()}
    audit = audit_archive(
        archive_path=args.archive,
        official_genes=genes,
        panel=panel,
        bank_perts=bank_perts,
        max_rss_gib=args.max_rss_gib,
    )
    target_rows, panel_rows = build_output_rows(audit)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "file_inventory.csv", list(audit.file_rows))
    write_csv(args.output_dir / "target_pathways.csv", target_rows)
    write_csv(args.output_dir / "panel_coverage.csv", panel_rows)

    archive_audit = {
        "archive": str(args.archive.resolve()),
        "archive_bytes": args.archive.stat().st_size,
        "archive_md5": audit.archive_md5,
        "expected_md5": EXPECTED_MD5,
        "archive_sha256": audit.archive_sha256,
        "zip_entries_total": len(zipfile.ZipFile(args.archive).infolist()),
        "de_table_files": len(audit.members),
        "unique_targets": len(audit.targets),
        "pathways": {
            pathway: sum(member.pathway == pathway for member in audit.members)
            for pathway in EXPECTED_PATHWAYS
        },
        "cell_lines": list(CELL_LINES),
        "mapped_vcc_genes_union": int(audit.gene_seen.sum()),
        "panel_targets": list(audit.panel_targets),
        "panel_target_count": len(audit.panel_targets),
        "panel_target_files": len(panel_rows),
        "shared_panel_targets_all_banks": list(audit.shared_panel_targets),
        "shared_panel_target_count_all_banks": len(audit.shared_panel_targets),
        "proxy_targets_shared_all_banks": len(audit.proxy_targets),
        "panel_axis_sha256": sha256_lines(panel),
        "gene_axis_sha256": sha256_lines(genes),
        "peak_rss_gib_after_audit": peak_rss_gib(),
    }
    write_json(args.output_dir / "audit.json", archive_audit)

    jiang = build_centered_effects(
        archive_path=args.archive,
        audit=audit,
        official_genes=genes,
        panel=panel,
        max_rss_gib=args.max_rss_gib,
    )
    compact_hashes = save_compact(args.output_dir, jiang, audit)
    proxy = run_proxy_development(
        jiang=jiang,
        audit=audit,
        official_genes=genes,
        panel=panel,
        banks=banks,
        fractions=tuple(args.fractions),
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    if proxy.get("rows"):
        write_csv(args.output_dir / "development_rows.csv", proxy["rows"])
        write_csv(
            args.output_dir / "development_summary.csv",
            proxy["summaries"],
        )
        write_csv(
            args.output_dir / "development_per_target.csv",
            proxy["per_target"],
        )

    panel_gate_checks = {
        "panel_targets_at_least_minimum": (
            len(audit.panel_targets) >= args.min_panel_targets
        ),
        "shared_panel_targets_at_least_minimum": (
            len(audit.shared_panel_targets) >= args.min_panel_targets
        ),
        "mapped_vcc_genes_at_least_minimum": (
            int(audit.gene_seen.sum()) >= args.min_mapped_genes
        ),
        "proxy_targets_at_least_minimum": (
            len(audit.proxy_targets) >= args.min_proxy_targets
        ),
    }
    panel_gate_passed = all(panel_gate_checks.values())
    result = {
        "decision": (
            "eligible_for_frozen_panel_evaluation"
            if panel_gate_passed
            else "negative_metadata_gate"
        ),
        "verdict": (
            "The archive passes structural and proxy-development checks, "
            "but its VCC-panel coverage is too small for a valid held-panel "
            "retrieval decision. No held-donor panel outcomes were opened."
            if not panel_gate_passed
            else "Panel metadata gate passed; held-panel evaluation is required."
        ),
        "protocol": {
            "archive_access": "Python zipfile streaming; no extraction",
            "effect_columns": [f"log2FC_{cell}" for cell in CELL_LINES],
            "pathway_handling": (
                "pathway is a batch; center each pathway/cell/gene on non-panel targets"
            ),
            "non_panel_centering": "leave one target out",
            "repeated_target_handling": (
                "average centered effects equally across pathway files"
            ),
            "selection_targets": "non-panel targets shared by all three banks",
            "held_panel_truth_opened": False,
            "panel_and_proxy_target_coordinates_excluded": True,
            "jiang_k562_omitted_when_k562_held": True,
            "no_h5ad_download": True,
            "no_cloud": True,
            "no_cell_eval2": True,
            "no_submission": True,
        },
        "gate": {
            "minimum_panel_targets": args.min_panel_targets,
            "minimum_proxy_targets": args.min_proxy_targets,
            "minimum_mapped_genes": args.min_mapped_genes,
            "checks": panel_gate_checks,
            "passed": panel_gate_passed,
        },
        "audit": archive_audit,
        "development_proxy": {
            key: value for key, value in proxy.items() if key != "per_target"
        },
        "compact_artifact_sha256": compact_hashes,
        "peak_rss_gib": peak_rss_gib(),
    }
    write_json(args.output_dir / "results.json", result)
    enforce_limit(args.max_rss_gib, "completed experiment")
    print(json.dumps(result["gate"], indent=2))
    print(json.dumps(result["development_proxy"].get("selected", {}), indent=2))
    print(f"decision={result['decision']}")
    print(f"peak_rss_gib={peak_rss_gib():.3f}")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--archive",
        type=Path,
        default=ROOT / "data/external/jiang2025/DE_results_all_pathway.zip",
    )
    parser.add_argument(
        "--panel-csv",
        type=Path,
        default=ROOT / "data/challenge/pert_counts.csv",
    )
    parser.add_argument(
        "--genes-csv",
        type=Path,
        default=ROOT / "data/challenge/gene_names.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=HERE / "out")
    parser.add_argument(
        "--fractions",
        nargs="+",
        type=float,
        default=list(DEFAULT_FRACTIONS),
    )
    parser.add_argument("--min-panel-targets", type=int, default=30)
    parser.add_argument("--min-proxy-targets", type=int, default=100)
    parser.add_argument("--min-mapped-genes", type=int, default=5_000)
    parser.add_argument("--chunk-rows", type=int, default=4)
    parser.add_argument("--max-rss-gib", type=float, default=0.90)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if not 0 < args.max_rss_gib <= 1.0:
        raise SystemExit("--max-rss-gib must be in (0, 1]")
    if args.chunk_rows < 1:
        raise SystemExit("--chunk-rows must be positive")
    if args.min_panel_targets < 2 or args.min_proxy_targets < 2:
        raise SystemExit("target gates must be at least 2")
    if any(not 0 < fraction < 1 for fraction in args.fractions):
        raise SystemExit("blend fractions must be in (0, 1)")
    run(args)


if __name__ == "__main__":
    main()

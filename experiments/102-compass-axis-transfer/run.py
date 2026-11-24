#!/usr/bin/env python3
"""Experiment 102: strict-LODO donor-axis transfer.

This is an independent implementation of the equations preregistered in
README.md. The staged commands intentionally separate prediction freezing from
opening held perturbation counts:

    uv run python experiments/102-compass-axis-transfer/run.py prepare --surface h1
    uv run python experiments/102-compass-axis-transfer/run.py freeze --surface h1
    uv run python experiments/102-compass-axis-transfer/run.py score --surface h1
    uv run python experiments/102-compass-axis-transfer/run.py summarize

All writes remain below this experiment directory. This module never invokes
cell-eval2, a cloud service, or challenge build/submission code.
"""

from __future__ import annotations

import argparse
import contextlib
import fcntl
import hashlib
import json
import os
import resource
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
TS = 50_000.0
MAX_ALLOWED_RSS_GIB = 24.0
N_NULL = 20
N_BOOT = 2_000
ARM_NAMES = ("A", "B", "C", "D", "E", "F")

DONOR_PATHS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}

SURFACES = {
    "h1": {
        "truth_donor": "h1_2025",
        "harness": ROOT / "data/harness/h1-test.h5ad",
        "sources": ("k562_gwps", "orion_hek293t"),
    },
    "hct116": {
        "truth_donor": "orion_hct116",
        "harness": ROOT / "data/harness/orion-hct116-thinned.h5ad",
        "sources": ("k562_gwps", "orion_hek293t"),
    },
    "hek293t": {
        "truth_donor": "orion_hek293t",
        "harness": ROOT / "data/harness/orion-hek293t-thinned.h5ad",
        "sources": ("k562_gwps", "orion_hct116"),
    },
}


def read_lines(path: Path) -> tuple[str, ...]:
    return tuple(line.strip() for line in path.read_text().splitlines() if line.strip())


def last_lookup(values: Sequence[str]) -> dict[str, int]:
    out: dict[str, int] = {}
    for index, value in enumerate(values):
        out[str(value)] = index
    return out


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(map(str, values)).encode()).hexdigest()


def sha256_array(values: np.ndarray) -> str:
    array = np.ascontiguousarray(values)
    return hashlib.sha256(array.view(np.uint8)).hexdigest()


def file_sha256(path: Path, block_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def stable_seed(*parts: object) -> int:
    payload = "|".join(map(str, parts)).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def validate_memory_limit(limit_gib: float) -> float:
    if not np.isfinite(limit_gib) or limit_gib <= 0:
        raise ValueError("--max-rss-gib must be finite and positive")
    if limit_gib > MAX_ALLOWED_RSS_GIB:
        raise ValueError(
            f"experiment 102 is capped at {MAX_ALLOWED_RSS_GIB:.0f} GiB; "
            f"received {limit_gib:g}"
        )
    return float(limit_gib)


def guard(limit_gib: float, stage: str) -> None:
    peak = peak_rss_gib()
    if peak >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak:.2f} GiB >= configured {limit_gib:.2f} GiB"
        )


def json_dump_atomic(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def npz_dump_atomic(path: Path, **arrays: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    temporary.replace(path)


def npy_dump_atomic(path: Path, array: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("wb") as handle:
        np.save(handle, array, allow_pickle=False)
    temporary.replace(path)


@contextlib.contextmanager
def experiment_lock() -> Iterable[None]:
    WORK.mkdir(parents=True, exist_ok=True)
    path = WORK / "experiment.lock"
    with path.open("a+") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                "another experiment-102 data process holds work/experiment.lock"
            ) from exc
        yield


def load_panel() -> tuple[str, ...]:
    path = ROOT / "data/challenge/pert_counts.csv"
    frame = pd.read_csv(path)
    if "target_gene" not in frame:
        raise ValueError(f"{path} has no target_gene column")
    panel = tuple(frame["target_gene"].astype(str))
    if len(panel) != 300 or len(set(panel)) != 300:
        raise ValueError(f"expected 300 unique VCC panel targets, found {len(panel)}")
    return panel


def build_allowed_mask(
    genes: Sequence[str],
    panel: Iterable[str],
    held_targets: Iterable[str],
) -> np.ndarray:
    forbidden = set(map(str, panel)) | set(map(str, held_targets))
    return np.asarray([str(gene) not in forbidden for gene in genes], dtype=bool)


def common_fit_targets(
    source_perts: Mapping[str, Sequence[str]],
    panel: Iterable[str],
    held_targets: Iterable[str],
) -> tuple[str, ...]:
    if not source_perts:
        raise ValueError("no source donors")
    shared: set[str] | None = None
    for values in source_perts.values():
        current = set(map(str, values))
        shared = current if shared is None else shared & current
    assert shared is not None
    shared -= set(map(str, panel)) | set(map(str, held_targets))
    result = tuple(sorted(shared))
    if len(result) < 2:
        raise ValueError("fewer than two shared non-panel source perturbations")
    return result


def normalize(vector: np.ndarray, label: str) -> np.ndarray:
    value = np.asarray(vector, dtype=np.float64)
    norm = float(np.linalg.norm(value))
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError(f"{label} has zero or non-finite norm")
    return value / norm


@dataclass
class SourceGeometry:
    name: str
    mean_scaled: np.ndarray
    axis: np.ndarray
    beta_mean: float
    scale: float
    control: np.ndarray
    gene_mask: np.ndarray
    n_fit: int


def fit_source_geometry_from_rows(
    name: str,
    fit_delta: np.ndarray,
    control: np.ndarray,
    allowed_mask: np.ndarray,
) -> SourceGeometry:
    """Fit one source geometry after applying the forbidden-coordinate mask.

    This pure implementation is also the reference used by synthetic tests.
    """

    values = np.asarray(fit_delta, dtype=np.float64).copy()
    allowed = np.asarray(allowed_mask, dtype=bool)
    if values.ndim != 2 or values.shape[1] != len(allowed):
        raise ValueError("fit_delta and allowed_mask shapes do not agree")
    if values.shape[0] < 2:
        raise ValueError("at least two common fit rows are required")
    values[:, ~allowed] = 0.0
    mean = values.mean(axis=0)
    centered = values - mean
    norms = np.linalg.norm(centered, axis=1)
    positive = norms[np.isfinite(norms) & (norms > 0)]
    if not len(positive):
        raise ValueError(f"{name}: all common-row centered effects are zero")
    scale = float(np.median(positive))
    mean_scaled = mean / scale
    axis = normalize(mean_scaled, f"{name} source mean")
    beta_mean = float(np.mean((values / scale) @ axis))
    source_control = np.asarray(control, dtype=np.float64).copy()
    if source_control.shape != (values.shape[1],):
        raise ValueError("control and fit_delta gene axes do not agree")
    source_control[~allowed] = 0.0
    return SourceGeometry(
        name=name,
        mean_scaled=mean_scaled,
        axis=axis,
        beta_mean=beta_mean,
        scale=scale,
        control=source_control,
        gene_mask=allowed.copy(),
        n_fit=values.shape[0],
    )


class DonorArray:
    """Memory-mapped compact donor arrays with explicit symbol projections."""

    def __init__(self, name: str, path: Path):
        self.name = name
        self.path = path
        self.perts = read_lines(path / "perts.txt")
        self.genes = read_lines(path / "genes.txt")
        self.pert_index = last_lookup(self.perts)
        self.gene_index = last_lookup(self.genes)
        self.delta = np.load(path / "delta.npy", mmap_mode="r")
        self.b_ctrl = np.load(path / "b_ctrl.npy", mmap_mode="r")
        if self.delta.shape != (len(self.perts), len(self.genes)):
            raise ValueError(
                f"{name}: delta shape {self.delta.shape} disagrees with axes"
            )
        if self.b_ctrl.shape != (len(self.genes),):
            raise ValueError(
                f"{name}: b_ctrl shape {self.b_ctrl.shape} disagrees with genes"
            )

    def projection(
        self, target_genes: Sequence[str], allowed_mask: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        tcols: list[int] = []
        dcols: list[int] = []
        allowed = np.asarray(allowed_mask, dtype=bool)
        for target_col, gene in enumerate(target_genes):
            donor_col = self.gene_index.get(str(gene))
            if allowed[target_col] and donor_col is not None:
                tcols.append(target_col)
                dcols.append(donor_col)
        if not tcols:
            raise ValueError(f"{self.name}: no eligible genes overlap the target axis")
        return np.asarray(tcols, dtype=np.int64), np.asarray(dcols, dtype=np.int64)

    def project_control(
        self, target_genes: Sequence[str], allowed_mask: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        tcols, dcols = self.projection(target_genes, allowed_mask)
        out = np.zeros(len(target_genes), dtype=np.float64)
        out[tcols] = np.asarray(self.b_ctrl[dcols], dtype=np.float64)
        mask = np.zeros(len(target_genes), dtype=bool)
        mask[tcols] = True
        return out, mask

    def project_target(
        self,
        target: str,
        target_genes: Sequence[str],
        allowed_mask: np.ndarray,
    ) -> np.ndarray | None:
        row = self.pert_index.get(str(target))
        if row is None:
            return None
        tcols, dcols = self.projection(target_genes, allowed_mask)
        out = np.zeros(len(target_genes), dtype=np.float64)
        out[tcols] = np.asarray(self.delta[row, dcols], dtype=np.float64)
        return out


def fit_source_geometry_streaming(
    donor: DonorArray,
    fit_targets: Sequence[str],
    target_genes: Sequence[str],
    allowed_mask: np.ndarray,
    max_rss_gib: float,
    block_rows: int = 128,
) -> SourceGeometry:
    """Two-pass source fit without materializing the donor-by-gene matrix."""

    tcols, dcols = donor.projection(target_genes, allowed_mask)
    rows = np.asarray([donor.pert_index[str(target)] for target in fit_targets])
    mean_shared = np.zeros(len(tcols), dtype=np.float64)
    for start in range(0, len(rows), block_rows):
        row_block = rows[start : start + block_rows]
        block = np.asarray(
            donor.delta[np.ix_(row_block, dcols)],
            dtype=np.float64,
        )
        mean_shared += block.sum(axis=0)
        guard(max_rss_gib, f"{donor.name} source-mean pass")
    mean_shared /= len(rows)

    norms = np.empty(len(rows), dtype=np.float64)
    for start in range(0, len(rows), block_rows):
        stop = min(start + block_rows, len(rows))
        row_block = rows[start:stop]
        block = np.asarray(
            donor.delta[np.ix_(row_block, dcols)],
            dtype=np.float64,
        )
        block -= mean_shared
        norms[start:stop] = np.sqrt(np.einsum("ij,ij->i", block, block))
        guard(max_rss_gib, f"{donor.name} source-scale pass")
    positive = norms[np.isfinite(norms) & (norms > 0)]
    if not len(positive):
        raise ValueError(f"{donor.name}: all common-row centered effects are zero")
    scale = float(np.median(positive))

    mean_scaled = np.zeros(len(target_genes), dtype=np.float64)
    mean_scaled[tcols] = mean_shared / scale
    axis = normalize(mean_scaled, f"{donor.name} source mean")
    beta_total = 0.0
    for start in range(0, len(rows), block_rows):
        row_block = rows[start : start + block_rows]
        block = np.asarray(
            donor.delta[np.ix_(row_block, dcols)],
            dtype=np.float64,
        )
        beta_total += float(((block / scale) @ axis[tcols]).sum())
        guard(max_rss_gib, f"{donor.name} beta-mean pass")
    control, gene_mask = donor.project_control(target_genes, allowed_mask)
    return SourceGeometry(
        name=donor.name,
        mean_scaled=mean_scaled,
        axis=axis,
        beta_mean=beta_total / len(rows),
        scale=scale,
        control=control,
        gene_mask=gene_mask,
        n_fit=len(rows),
    )


def control_axis_weights(
    target_control: np.ndarray,
    geometries: Sequence[SourceGeometry],
    allowed_mask: np.ndarray,
) -> np.ndarray:
    allowed = np.asarray(allowed_mask, dtype=bool)
    target = np.asarray(target_control, dtype=np.float64).copy()
    target[~allowed] = 0.0
    similarities = []
    for geometry in geometries:
        use = allowed & geometry.gene_mask
        left = target[use]
        right = geometry.control[use]
        denom = float(np.linalg.norm(left) * np.linalg.norm(right))
        cosine = float(left @ right) / denom if denom > 0 else 0.0
        similarities.append(max(0.0, cosine))
    weights = np.asarray(similarities, dtype=np.float64)
    if not np.isfinite(weights).all() or weights.sum() <= 0:
        weights = np.ones(len(geometries), dtype=np.float64)
    return weights / weights.sum()


def carrier_mean(values: np.ndarray, carrier: np.ndarray) -> np.ndarray:
    """Average source-first arrays over available source rows per target."""

    array = np.asarray(values, dtype=np.float64)
    present = np.asarray(carrier, dtype=bool)
    if array.shape[:2] != present.shape:
        raise ValueError("values and carrier source/target dimensions disagree")
    weights = present.astype(np.float64)
    denominator = weights.sum(axis=0)
    if array.ndim == 3:
        numerator = (array * weights[:, :, None]).sum(axis=0)
        return np.divide(
            numerator,
            denominator[:, None],
            out=np.zeros_like(numerator),
            where=denominator[:, None] > 0,
        )
    if array.ndim == 2:
        numerator = (array * weights).sum(axis=0)
        return np.divide(
            numerator,
            denominator,
            out=np.zeros_like(numerator),
            where=denominator > 0,
        )
    raise ValueError("carrier_mean supports 2-D or 3-D source-first arrays")


def compose_arms(
    delta_scaled: np.ndarray,
    carrier: np.ndarray,
    geometries: Sequence[SourceGeometry],
    target_control: np.ndarray,
    allowed_mask: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Compose the six preregistered arms from already scaled source rows."""

    delta = np.asarray(delta_scaled, dtype=np.float64).copy()
    present = np.asarray(carrier, dtype=bool)
    allowed = np.asarray(allowed_mask, dtype=bool)
    if delta.ndim != 3 or delta.shape[:2] != present.shape:
        raise ValueError("delta_scaled must be source x target x gene")
    if delta.shape[0] != len(geometries) or delta.shape[2] != len(allowed):
        raise ValueError("geometry or gene dimensions disagree")
    delta[:, :, ~allowed] = 0.0

    means = np.stack([geometry.mean_scaled for geometry in geometries])
    axes = np.stack([geometry.axis for geometry in geometries])
    beta_means = np.asarray([geometry.beta_mean for geometry in geometries])
    means[:, ~allowed] = 0.0
    axes[:, ~allowed] = 0.0

    beta = np.einsum("ctg,cg->ct", delta, axes)
    residual = delta - beta[:, :, None] * axes[:, None, :]
    residual *= present[:, :, None]
    centered = (delta - means[:, None, :]) * present[:, :, None]
    raw = delta * present[:, :, None]

    residual_mean = carrier_mean(residual, present)
    centered_beta = beta - beta_means[:, None]
    primary_axis = normalize(axes.mean(axis=0), "donor-mean target axis")
    basal_weights = control_axis_weights(target_control, geometries, allowed)
    basal_axis = normalize(
        np.einsum("c,cg->g", basal_weights, axes),
        "basal-weighted target axis",
    )

    a = carrier_mean(centered, present)
    b = carrier_mean(raw, present)
    c = residual_mean
    d = c + carrier_mean(centered_beta, present)[:, None] * primary_axis
    e = c + carrier_mean(beta, present)[:, None] * primary_axis
    f = c + carrier_mean(centered_beta, present)[:, None] * basal_axis
    arms = {"A": a, "B": b, "C": c, "D": d, "E": e, "F": f}
    for name, value in arms.items():
        value[:, ~allowed] = 0.0
        if not np.isfinite(value).all():
            raise ValueError(f"arm {name} contains non-finite values")

    return arms, {
        "source_beta": beta,
        "source_centered_beta": centered_beta,
        "source_residual": residual,
        "residual_mean": residual_mean,
        "primary_axis": primary_axis,
        "basal_axis": basal_axis,
        "basal_weights": basal_weights,
        "coverage": present.sum(axis=0),
    }


def permuted_centered_beta(
    centered_beta: np.ndarray,
    carrier: np.ndarray,
    *,
    surface: str,
    draw: int,
) -> np.ndarray:
    values = np.asarray(centered_beta, dtype=np.float64)
    present = np.asarray(carrier, dtype=bool)
    if values.shape != present.shape:
        raise ValueError("centered_beta and carrier shapes disagree")
    out = values.copy()
    for source in range(len(values)):
        rows = np.flatnonzero(present[source])
        rng = np.random.default_rng(
            stable_seed("experiment-102-beta-null", surface, draw, source)
        )
        out[source, rows] = values[source, rng.permutation(rows)]
        out[source, ~present[source]] = 0.0
    return out


def beta_null_prediction(
    residual_mean: np.ndarray,
    centered_beta: np.ndarray,
    carrier: np.ndarray,
    primary_axis: np.ndarray,
    *,
    surface: str,
    draw: int,
    allowed_mask: np.ndarray,
) -> np.ndarray:
    permuted = permuted_centered_beta(
        centered_beta, carrier, surface=surface, draw=draw
    )
    loading = carrier_mean(permuted, carrier)
    out = np.asarray(residual_mean, dtype=np.float64) + loading[:, None] * np.asarray(
        primary_axis, dtype=np.float64
    )
    out[:, ~np.asarray(allowed_mask, dtype=bool)] = 0.0
    return out


def load_surface_controls(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read gene/target labels and control counts, never perturbation values."""

    handle = ad.read_h5ad(path, backed="r")
    try:
        genes = np.asarray(handle.var_names).astype(str)
        labels = np.asarray(handle.obs["target_gene"].astype(str))
        control_rows = np.flatnonzero(labels == "non-targeting")
        if not len(control_rows):
            raise ValueError(f"{path}: no non-targeting controls")
        xctrl = handle.X[control_rows, :]
        if sp.issparse(xctrl):
            sums = np.asarray(xctrl.sum(0)).ravel().astype(np.float64)
        else:
            sums = np.asarray(xctrl, dtype=np.float64).sum(axis=0)
    finally:
        handle.file.close()
    if sums.sum() <= 0:
        raise ValueError(f"{path}: control total is zero")
    bctrl = np.log1p(TS * sums / sums.sum())
    targets = np.asarray(sorted(set(labels) - {"non-targeting"}))
    return genes, targets, bctrl


def surface_dir(surface: str) -> Path:
    return WORK / surface


def prepared_paths(surface: str) -> tuple[Path, Path]:
    out = surface_dir(surface)
    return out / "prepared.npz", out / "prepare-audit.json"


def source_geometry_path(surface: str, source: str) -> Path:
    return surface_dir(surface) / f"source-{source}.npz"


def save_geometry(path: Path, geometry: SourceGeometry) -> None:
    npz_dump_atomic(
        path,
        name=np.asarray(geometry.name),
        mean_scaled=geometry.mean_scaled.astype(np.float64),
        axis=geometry.axis.astype(np.float64),
        beta_mean=np.asarray(geometry.beta_mean, dtype=np.float64),
        scale=np.asarray(geometry.scale, dtype=np.float64),
        control=geometry.control.astype(np.float64),
        gene_mask=geometry.gene_mask,
        n_fit=np.asarray(geometry.n_fit, dtype=np.int64),
    )


def load_geometry(path: Path) -> SourceGeometry:
    with np.load(path, allow_pickle=False) as data:
        return SourceGeometry(
            name=str(data["name"]),
            mean_scaled=data["mean_scaled"].astype(np.float64),
            axis=data["axis"].astype(np.float64),
            beta_mean=float(data["beta_mean"]),
            scale=float(data["scale"]),
            control=data["control"].astype(np.float64),
            gene_mask=data["gene_mask"].astype(bool),
            n_fit=int(data["n_fit"]),
        )


def prepare_surface(surface: str, max_rss_gib: float) -> None:
    spec = SURFACES[surface]
    out = surface_dir(surface)
    out.mkdir(parents=True, exist_ok=True)
    prepared_path, audit_path = prepared_paths(surface)
    sources = tuple(spec["sources"])
    source_paths = [source_geometry_path(surface, source) for source in sources]
    if (
        prepared_path.exists()
        and audit_path.exists()
        and all(path.exists() for path in source_paths)
    ):
        print(f"{surface}: prepare already complete")
        return

    genes, targets, target_control = load_surface_controls(Path(spec["harness"]))
    panel = load_panel()
    allowed = build_allowed_mask(genes, panel, targets)
    if int(allowed.sum()) < 2:
        raise ValueError(f"{surface}: fewer than two eligible response genes")
    donors = {name: DonorArray(name, DONOR_PATHS[name]) for name in sources}
    fit_targets = common_fit_targets(
        {name: donor.perts for name, donor in donors.items()},
        panel,
        targets,
    )
    if not prepared_path.exists():
        npz_dump_atomic(
            prepared_path,
            genes=genes,
            targets=targets,
            target_control=target_control.astype(np.float64),
            allowed_mask=allowed,
            panel=np.asarray(panel),
            fit_targets=np.asarray(fit_targets),
            sources=np.asarray(sources),
        )

    geometries: list[SourceGeometry] = []
    for source in sources:
        path = source_geometry_path(surface, source)
        if path.exists():
            geometry = load_geometry(path)
        else:
            geometry = fit_source_geometry_streaming(
                donors[source],
                fit_targets,
                genes,
                allowed,
                max_rss_gib,
            )
            save_geometry(path, geometry)
        geometries.append(geometry)
        guard(max_rss_gib, f"{surface} prepared {source}")

    audit = {
        "surface": surface,
        "held_truth_donor": spec["truth_donor"],
        "sources": list(sources),
        "genes": int(len(genes)),
        "targets": int(len(targets)),
        "eligible_coordinates": int(allowed.sum()),
        "forbidden_coordinates": int((~allowed).sum()),
        "panel_targets": len(panel),
        "shared_non_panel_fit_targets": len(fit_targets),
        "genes_sha256": sha256_lines(genes),
        "targets_sha256": sha256_lines(targets),
        "allowed_mask_sha256": sha256_array(allowed.astype(np.uint8)),
        "fit_targets_sha256": sha256_lines(fit_targets),
        "source_geometry": {
            geometry.name: {
                "n_fit": geometry.n_fit,
                "scale": geometry.scale,
                "beta_mean": geometry.beta_mean,
                "axis_sha256": sha256_array(geometry.axis),
                "eligible_shared_genes": int(geometry.gene_mask.sum()),
            }
            for geometry in geometries
        },
        "held_perturbation_values_read": False,
        "peak_rss_gib": peak_rss_gib(),
    }
    json_dump_atomic(audit_path, audit)
    print(json.dumps(audit, indent=2))


def load_prepared(
    surface: str,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    tuple[str, ...],
    list[SourceGeometry],
]:
    prepared_path, audit_path = prepared_paths(surface)
    if not prepared_path.exists() or not audit_path.exists():
        raise FileNotFoundError(f"{surface}: run prepare first")
    with np.load(prepared_path, allow_pickle=False) as data:
        genes = data["genes"].astype(str)
        targets = data["targets"].astype(str)
        target_control = data["target_control"].astype(np.float64)
        allowed = data["allowed_mask"].astype(bool)
        sources = tuple(data["sources"].astype(str))
    expected_sources = tuple(SURFACES[surface]["sources"])
    if sources != expected_sources:
        raise ValueError(f"{surface}: prepared sources changed")
    geometries = [
        load_geometry(source_geometry_path(surface, name)) for name in sources
    ]
    return genes, targets, target_control, allowed, sources, geometries


def freeze_surface(surface: str, max_rss_gib: float, n_null: int) -> None:
    out = surface_dir(surface)
    audit_path = out / "freeze-audit.json"
    null_dir = out / "null"
    required_nulls = [
        null_dir / f"beta-permutation-{draw:03d}.npy" for draw in range(n_null)
    ]
    if (
        audit_path.exists()
        and (out / "frozen-arms.npz").exists()
        and (out / "components.npz").exists()
        and all(path.exists() for path in required_nulls)
    ):
        print(f"{surface}: freeze already complete")
        return

    genes, targets, target_control, allowed, sources, geometries = load_prepared(
        surface
    )
    donors = {name: DonorArray(name, DONOR_PATHS[name]) for name in sources}
    n_sources, n_targets, n_genes = len(sources), len(targets), len(genes)
    delta_scaled = np.zeros((n_sources, n_targets, n_genes), dtype=np.float32)
    carrier = np.zeros((n_sources, n_targets), dtype=bool)
    for source_index, source in enumerate(sources):
        donor = donors[source]
        scale = geometries[source_index].scale
        for target_index, target in enumerate(targets):
            value = donor.project_target(str(target), genes, allowed)
            if value is None:
                continue
            delta_scaled[source_index, target_index] = (value / scale).astype(
                np.float32
            )
            carrier[source_index, target_index] = True
        guard(max_rss_gib, f"{surface} projected {source}")

    arms, components = compose_arms(
        delta_scaled, carrier, geometries, target_control, allowed
    )
    components_path = out / "components.npz"
    if not components_path.exists():
        npz_dump_atomic(
            components_path,
            source_centered_beta=components["source_centered_beta"].astype(np.float32),
            residual_mean=components["residual_mean"].astype(np.float32),
            primary_axis=components["primary_axis"].astype(np.float32),
            carrier=carrier,
            coverage=components["coverage"].astype(np.int16),
        )
    arms_path = out / "frozen-arms.npz"
    if not arms_path.exists():
        npz_dump_atomic(
            arms_path,
            genes=genes,
            targets=targets,
            allowed_mask=allowed,
            coverage=components["coverage"].astype(np.int16),
            basal_weights=components["basal_weights"].astype(np.float64),
            **{name: arms[name].astype(np.float32) for name in ARM_NAMES},
        )
    del delta_scaled
    guard(max_rss_gib, f"{surface} frozen real arms")

    null_dir.mkdir(parents=True, exist_ok=True)
    for draw, path in enumerate(required_nulls):
        if path.exists():
            continue
        prediction = beta_null_prediction(
            components["residual_mean"],
            components["source_centered_beta"],
            carrier,
            components["primary_axis"],
            surface=surface,
            draw=draw,
            allowed_mask=allowed,
        )
        npy_dump_atomic(path, prediction.astype(np.float32))
        guard(max_rss_gib, f"{surface} null draw {draw}")

    coverage = components["coverage"]
    coverage_hash = sha256_array(coverage.astype(np.int16))
    audit = {
        "surface": surface,
        "sources": list(sources),
        "targets": int(n_targets),
        "genes": int(n_genes),
        "covered_targets": int((coverage > 0).sum()),
        "uncovered_targets": targets[coverage == 0].tolist(),
        "source_coverage_sha256": coverage_hash,
        "coverage_identical_for_all_arms": True,
        "forbidden_nonzero_by_arm": {
            name: int(np.count_nonzero(arms[name][:, ~allowed])) for name in ARM_NAMES
        },
        "basal_weights": components["basal_weights"].tolist(),
        "n_null": n_null,
        "held_perturbation_values_read": False,
        "frozen_arms_sha256": file_sha256(arms_path),
        "components_sha256": file_sha256(components_path),
        "peak_rss_gib": peak_rss_gib(),
    }
    if any(audit["forbidden_nonzero_by_arm"].values()):
        raise AssertionError("a forbidden coordinate survived into a frozen arm")
    json_dump_atomic(audit_path, audit)
    print(json.dumps(audit, indent=2))


def load_truth_effects(
    path: Path, targets: np.ndarray, genes: np.ndarray, max_rss_gib: float
) -> np.ndarray:
    """Open held perturbation values only in the score stage."""

    data = ad.read_h5ad(path)
    try:
        if not np.array_equal(np.asarray(data.var_names).astype(str), genes):
            raise ValueError("held truth gene order changed")
        labels = np.asarray(data.obs["target_gene"].astype(str))
        x = data.X.tocsr() if sp.issparse(data.X) else sp.csr_matrix(data.X)
        control_sum = (
            np.asarray(x[labels == "non-targeting"].sum(axis=0))
            .ravel()
            .astype(np.float64)
        )
        if control_sum.sum() <= 0:
            raise ValueError("held control total is zero")
        control = np.log1p(TS * control_sum / control_sum.sum())
        rows = []
        for target in targets:
            selected = labels == str(target)
            if not selected.any():
                raise ValueError(f"held truth is missing target {target}")
            total = np.asarray(x[selected].sum(axis=0)).ravel().astype(np.float64)
            if total.sum() <= 0:
                raise ValueError(f"held truth target {target} has zero total")
            rows.append(np.log1p(TS * total / total.sum()) - control)
            guard(max_rss_gib, f"held truth aggregation {target}")
        return np.vstack(rows)
    finally:
        del data


def pds_per_target(
    predicted: np.ndarray, real: np.ndarray, score_mask: np.ndarray
) -> np.ndarray:
    predicted = np.asarray(predicted, dtype=np.float64)
    real = np.asarray(real, dtype=np.float64)
    mask = np.asarray(score_mask, dtype=bool)
    if predicted.shape != real.shape or predicted.shape[1] != len(mask):
        raise ValueError("predicted, real, and score-mask shapes disagree")
    if len(predicted) < 2:
        raise ValueError("PDS requires at least two gallery targets")
    p = predicted[:, mask]
    r = real[:, mask]
    pn = np.linalg.norm(p, axis=1)
    rn = np.linalg.norm(r, axis=1)
    similarity = (p @ r.T) / np.maximum(pn[:, None] * rn[None, :], 1e-300)
    similarity[pn == 0] = -1.0
    similarity[:, rn == 0] = -1.0
    own = np.diag(similarity)[:, None]
    ranks = (similarity > own).sum(axis=1)
    ties = (similarity == own).sum(axis=1) - 1
    return 1.0 - (ranks + 0.5 * ties) / (len(predicted) - 1)


def bootstrap_ci(values: np.ndarray, n_boot: int, seed: int) -> dict[str, float]:
    values = np.asarray(values, dtype=np.float64)
    if not len(values):
        raise ValueError("cannot bootstrap an empty target vector")
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot, dtype=np.float64)
    for draw in range(n_boot):
        rows = rng.integers(0, len(values), len(values))
        draws[draw] = values[rows].mean()
    return {
        "estimate": float(values.mean()),
        "lo": float(np.percentile(draws, 2.5)),
        "hi": float(np.percentile(draws, 97.5)),
    }


def score_surface(surface: str, max_rss_gib: float, n_boot: int) -> None:
    out = surface_dir(surface)
    metrics_path = out / "metrics.json"
    scores_path = out / "per-target-scores.npz"
    if metrics_path.exists() and scores_path.exists():
        print(f"{surface}: score already complete")
        return
    freeze_audit_path = out / "freeze-audit.json"
    if not freeze_audit_path.exists():
        raise FileNotFoundError(f"{surface}: run freeze first")
    freeze_audit = json.loads(freeze_audit_path.read_text())
    n_null = int(freeze_audit["n_null"])
    null_paths = [
        out / "null" / f"beta-permutation-{draw:03d}.npy" for draw in range(n_null)
    ]
    if not all(path.exists() for path in null_paths):
        raise FileNotFoundError(f"{surface}: frozen beta-null set is incomplete")

    with np.load(out / "frozen-arms.npz", allow_pickle=False) as frozen:
        genes = frozen["genes"].astype(str)
        targets = frozen["targets"].astype(str)
        allowed = frozen["allowed_mask"].astype(bool)
        coverage = frozen["coverage"].astype(np.int64)
        arms = {name: frozen[name].astype(np.float64) for name in ARM_NAMES}
    real = load_truth_effects(
        Path(SURFACES[surface]["harness"]), targets, genes, max_rss_gib
    )
    pds = {name: pds_per_target(arms[name], real, allowed) for name in ARM_NAMES}
    delta = pds["D"] - pds["A"]
    null_delta = np.empty((n_null, len(targets)), dtype=np.float64)
    for draw, path in enumerate(null_paths):
        null_prediction = np.load(path, mmap_mode="r")
        null_delta[draw] = pds_per_target(null_prediction, real, allowed) - pds["A"]
        guard(max_rss_gib, f"{surface} score null {draw}")

    coverage_hash = sha256_array(coverage.astype(np.int16))
    coverage_equal = coverage_hash == str(
        freeze_audit["source_coverage_sha256"]
    ) and bool(freeze_audit["coverage_identical_for_all_arms"])
    result = {
        "surface": surface,
        "metric": "raw_uncentered_complete_gallery_pds",
        "targets": int(len(targets)),
        "eligible_coordinates": int(allowed.sum()),
        "pds": {name: float(values.mean()) for name, values in pds.items()},
        "D_minus_A": float(delta.mean()),
        "D_minus_A_bootstrap": bootstrap_ci(
            delta, n_boot, stable_seed("experiment-102", surface, "bootstrap")
        ),
        "null_mean_deltas": null_delta.mean(axis=1).tolist(),
        "null_surface_p95": float(np.percentile(null_delta.mean(axis=1), 95)),
        "coverage_identical_A_D": coverage_equal,
        "source_coverage_sha256": coverage_hash,
        "held_perturbation_values_read": True,
        "peak_rss_gib": peak_rss_gib(),
    }
    npz_dump_atomic(
        scores_path,
        targets=targets,
        D_minus_A=delta,
        null_D_minus_A=null_delta,
        **{f"pds_{name}": values for name, values in pds.items()},
    )
    json_dump_atomic(metrics_path, result)
    print(json.dumps(result, indent=2))


def macro_bootstrap_ci(
    surface_deltas: Mapping[str, np.ndarray], n_boot: int, seed: int
) -> dict[str, float]:
    if set(surface_deltas) != set(SURFACES):
        raise ValueError("macro bootstrap requires all three frozen surfaces")
    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot, dtype=np.float64)
    for draw in range(n_boot):
        surface_means = []
        for surface in SURFACES:
            values = np.asarray(surface_deltas[surface], dtype=np.float64)
            rows = rng.integers(0, len(values), len(values))
            surface_means.append(float(values[rows].mean()))
        draws[draw] = float(np.mean(surface_means))
    estimate = float(
        np.mean([np.asarray(values).mean() for values in surface_deltas.values()])
    )
    return {
        "estimate": estimate,
        "lo": float(np.percentile(draws, 2.5)),
        "hi": float(np.percentile(draws, 97.5)),
    }


def evaluate_primary_gate(
    surface_deltas: Mapping[str, np.ndarray],
    null_surface_means: Mapping[str, np.ndarray],
    coverage_equal: Mapping[str, bool],
    *,
    n_boot: int,
) -> dict[str, object]:
    surface_means = {
        name: float(np.asarray(surface_deltas[name]).mean()) for name in SURFACES
    }
    bootstrap = macro_bootstrap_ci(
        surface_deltas,
        n_boot,
        stable_seed("experiment-102", "macro-bootstrap"),
    )
    null_lengths = {len(np.asarray(null_surface_means[name])) for name in SURFACES}
    if len(null_lengths) != 1:
        raise ValueError("surface null arrays have different draw counts")
    null_macro = np.mean(
        np.vstack([np.asarray(null_surface_means[name]) for name in SURFACES]),
        axis=0,
    )
    null_p95 = float(np.percentile(null_macro, 95))
    checks = {
        "D_beats_A_every_surface": all(value > 0 for value in surface_means.values()),
        "macro_delta_at_least_0_003": bootstrap["estimate"] >= 0.003,
        "macro_bootstrap_lower_bound_above_zero": bootstrap["lo"] > 0,
        "macro_delta_above_beta_null_p95": bootstrap["estimate"] > null_p95,
        "source_coverage_unchanged_every_surface": all(
            bool(coverage_equal[name]) for name in SURFACES
        ),
    }
    return {
        "decision": "GO" if all(checks.values()) else "NO-GO",
        "surface_D_minus_A": surface_means,
        "macro_D_minus_A": bootstrap["estimate"],
        "macro_bootstrap": bootstrap,
        "beta_permutation_null_macro_values": null_macro.tolist(),
        "beta_permutation_null_macro_p95": null_p95,
        "checks": checks,
    }


def summarize(n_boot: int) -> None:
    output = WORK / "summary.json"
    if output.exists():
        print(output.read_text(), end="")
        return
    surface_deltas: dict[str, np.ndarray] = {}
    null_means: dict[str, np.ndarray] = {}
    coverage: dict[str, bool] = {}
    diagnostics: dict[str, object] = {}
    for surface in SURFACES:
        metrics_path = surface_dir(surface) / "metrics.json"
        scores_path = surface_dir(surface) / "per-target-scores.npz"
        if not metrics_path.exists() or not scores_path.exists():
            raise FileNotFoundError(f"{surface}: run score before summarize")
        metrics = json.loads(metrics_path.read_text())
        with np.load(scores_path, allow_pickle=False) as scores:
            surface_deltas[surface] = scores["D_minus_A"].astype(np.float64)
            null_means[surface] = scores["null_D_minus_A"].mean(axis=1)
        coverage[surface] = bool(metrics["coverage_identical_A_D"])
        diagnostics[surface] = {
            "pds": metrics["pds"],
            "D_minus_A": metrics["D_minus_A"],
        }
    result = evaluate_primary_gate(surface_deltas, null_means, coverage, n_boot=n_boot)
    result["diagnostic_arms"] = diagnostics
    json_dump_atomic(output, result)
    print(json.dumps(result, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("prepare", "freeze", "score", "summarize"))
    parser.add_argument("--surface", choices=tuple(SURFACES))
    parser.add_argument("--max-rss-gib", type=float, default=MAX_ALLOWED_RSS_GIB)
    parser.add_argument("--n-null", type=int, default=N_NULL)
    parser.add_argument("--n-boot", type=int, default=N_BOOT)
    args = parser.parse_args()
    limit = validate_memory_limit(args.max_rss_gib)
    if args.n_null < 1:
        parser.error("--n-null must be positive")
    if args.n_boot < 100:
        parser.error("--n-boot must be at least 100")
    if args.n_null != N_NULL:
        parser.error(f"the preregistered beta-permutation count is exactly {N_NULL}")
    if args.n_boot != N_BOOT:
        parser.error(f"the preregistered target-bootstrap count is exactly {N_BOOT}")
    if args.command == "summarize":
        if args.surface:
            parser.error("--surface is not used by summarize")
        summarize(args.n_boot)
        return
    if not args.surface:
        parser.error("--surface is required for prepare/freeze/score")
    with experiment_lock():
        if args.command == "prepare":
            prepare_surface(args.surface, limit)
        elif args.command == "freeze":
            freeze_surface(args.surface, limit, args.n_null)
        else:
            score_surface(args.surface, limit, args.n_boot)


if __name__ == "__main__":
    main()

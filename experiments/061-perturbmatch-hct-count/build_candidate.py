#!/usr/bin/env python3
"""Build the single experiment-061 HCT116 count-level candidate.

The existing experiment-040 K562-GWPS + Orion-HEK293T model is reconstructed
without changing its source processing. Only its K562 slot is changed: on
PerturbMatch-observed targets and genes, the incumbent normalized K562 effect
is averaged 50:50 with the independently normalized PerturbMatch estimate.
Missing PerturbMatch targets and genes retain the incumbent slot exactly.

The already fitted experiment-040 union whitener is loaded and audited. It is
never refit on the candidate. All 300 VCC target coordinates must be absent.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import resource
import sys
from pathlib import Path
from typing import Iterable

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from vcc2026.donor_transfer import Donor, TS, emit_counts  # noqa: E402
from vcc2026.multi_donor import MultiDonorModel  # noqa: E402
from vcc2026.whiten import Whitener  # noqa: E402

CONTROL = "non-targeting"
PERT_COL = "target_gene"
K562 = "k562_gwps"
HEK = "orion_hek293t"
GAMMA = 0.65
ALPHA = 1.08
MIN_CPM = 5.0
SPREAD = 0.25
BLEND = 0.50
SEED = 0
GIB = 1024**3


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if sys.platform == "darwin" else value / 1024**2


def memory_checkpoint(stage: str, limit_gib: float) -> None:
    peak = peak_rss_gib()
    print(f"[rss] {stage}: {peak:.3f} GiB", flush=True)
    if limit_gib > 0 and peak >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak:.3f} GiB reached limit {limit_gib:.3f} GiB"
        )


def read_lines(path: Path) -> tuple[str, ...]:
    return tuple(line.strip() for line in path.read_text().splitlines() if line.strip())


def read_csv_column(path: Path, column: str) -> tuple[str, ...]:
    values = tuple(
        pd.read_csv(path)[column].dropna().astype(str).str.strip().tolist()
    )
    if not values or len(values) != len(set(values)):
        raise ValueError(f"{path}: {column} must be non-empty and unique")
    return values


def file_sha256(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def array_sha256(array: np.ndarray) -> str:
    canonical = np.ascontiguousarray(array, dtype="<f8")
    return hashlib.sha256(canonical.tobytes()).hexdigest()


def axis_mapping(
    source_genes: Iterable[str], target_genes: Iterable[str]
) -> tuple[np.ndarray, np.ndarray]:
    """Match genes by symbol, retaining the last source duplicate like assignment."""
    target_index = {gene: index for index, gene in enumerate(target_genes)}
    source_for_target: dict[int, int] = {}
    for source_index, gene in enumerate(source_genes):
        target_index_value = target_index.get(gene)
        if target_index_value is not None:
            source_for_target[target_index_value] = source_index
    target_columns = np.asarray(sorted(source_for_target), dtype=np.int64)
    source_columns = np.asarray(
        [source_for_target[index] for index in target_columns], dtype=np.int64
    )
    return target_columns, source_columns


def prepare_perturbmatch(
    path: Path,
    target_genes: tuple[str, ...],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> dict[str, object]:
    """Compute baseline-compatible raw centering and median-norm scaling."""
    genes = read_lines(path / "genes.txt")
    perts = read_lines(path / "perts.txt")
    if len(genes) != len(set(genes)) or len(perts) != len(set(perts)):
        raise ValueError("PerturbMatch genes and perturbations must be unique")
    metadata = json.loads((path / "meta.json").read_text())
    if metadata.get("source_md5") != "0017076a5cc7b9979febf3bd19bd2157":
        raise ValueError("PerturbMatch source MD5 is not the official artifact")

    delta_path = path / "delta.npy"
    delta = np.load(delta_path, mmap_mode="r")
    if delta.shape != (len(perts), len(genes)):
        raise ValueError(
            f"PerturbMatch shape {delta.shape} != {(len(perts), len(genes))}"
        )
    target_columns, source_columns = axis_mapping(genes, target_genes)
    if not len(target_columns):
        raise ValueError("PerturbMatch has no genes on the target axis")

    mean = np.zeros(len(source_columns), dtype=np.float64)
    for lo in range(0, len(perts), chunk_rows):
        hi = min(lo + chunk_rows, len(perts))
        mean += np.asarray(
            delta[np.ix_(np.arange(lo, hi), source_columns)], dtype=np.float64
        ).sum(axis=0)
        memory_checkpoint("PerturbMatch center", max_rss_gib)
    mean /= len(perts)

    norms = np.empty(len(perts), dtype=np.float64)
    for lo in range(0, len(perts), chunk_rows):
        hi = min(lo + chunk_rows, len(perts))
        block = np.asarray(
            delta[np.ix_(np.arange(lo, hi), source_columns)], dtype=np.float64
        )
        block -= mean
        norms[lo:hi] = np.linalg.norm(block, axis=1)
        memory_checkpoint("PerturbMatch scale", max_rss_gib)
    scale = float(np.median(norms)) or 1.0
    return {
        "genes": genes,
        "perts": perts,
        "pert_index": {pert: index for index, pert in enumerate(perts)},
        "target_columns": target_columns,
        "source_columns": source_columns,
        "mean": mean,
        "scale": scale,
        "delta": delta,
        "metadata": metadata,
        "delta_path": delta_path,
    }


def blend_k562_slot(
    baseline_effect: np.ndarray,
    incumbent_standardized: np.ndarray,
    perturbmatch_standardized: np.ndarray,
    perturbmatch_target_columns: np.ndarray,
    *,
    effect_scale: float,
    blend: float = BLEND,
    source_weight: float = 0.5,
) -> np.ndarray:
    """Replace only observed K562-slot coordinates, preserving all others exactly."""
    output = np.asarray(baseline_effect, dtype=np.float64).copy()
    columns = np.asarray(perturbmatch_target_columns, dtype=np.int64)
    output[columns] += (
        source_weight
        * effect_scale
        * blend
        * (
            np.asarray(perturbmatch_standardized, dtype=np.float64)
            - np.asarray(incumbent_standardized, dtype=np.float64)[columns]
        )
    )
    return output


def compute_effects(
    model: MultiDonorModel,
    pm: dict[str, object],
    targets: tuple[str, ...],
    *,
    max_rss_gib: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    if model.names != [K562, HEK]:
        raise ValueError(f"source order must be {[K562, HEK]}, got {model.names}")
    missing_incumbent = [
        target for target in targets if target not in model.pert_index[K562]
    ]
    missing_hek = [target for target in targets if target not in model.pert_index[HEK]]
    if missing_incumbent or missing_hek:
        raise ValueError(
            "production-shaped HCT comparison requires both incumbent sources for "
            f"every target; missing K562={missing_incumbent}, HEK={missing_hek}"
        )

    weights = {name: 1.0 for name in model.names}
    baseline = np.vstack(
        [model.effect(target, 1.0, 0.0, weights=weights) for target in targets]
    ).astype(np.float64)
    candidate = baseline.copy()

    k_scale = float(model._scale[K562]) or 1.0
    pm_scale = float(pm["scale"])
    pm_target_columns = np.asarray(pm["target_columns"], dtype=np.int64)
    pm_source_columns = np.asarray(pm["source_columns"], dtype=np.int64)
    pm_mean = np.asarray(pm["mean"], dtype=np.float64)
    pm_delta = pm["delta"]
    pm_pert_index = pm["pert_index"]

    pm_missing_targets: list[str] = []
    modified_targets: list[str] = []
    for row, target in enumerate(targets):
        pm_row = pm_pert_index.get(target)
        if pm_row is None:
            pm_missing_targets.append(target)
            if not np.array_equal(candidate[row], baseline[row]):
                raise AssertionError("missing PerturbMatch target did not preserve baseline")
            continue

        k_row = model.pert_index[K562][target]
        incumbent = np.zeros(len(model.target_genes), dtype=np.float64)
        incumbent[model.tcols[K562]] = (
            model.s[K562][k_row][model.dcols[K562]] / k_scale
        )
        pm_standardized = (
            np.asarray(pm_delta[pm_row, pm_source_columns], dtype=np.float64)
            - pm_mean
        ) / pm_scale
        before = candidate[row].copy()
        candidate[row] = blend_k562_slot(
            baseline[row],
            incumbent,
            pm_standardized,
            pm_target_columns,
            effect_scale=float(model.effect_scale),
        )
        missing_gene_mask = np.ones(len(model.target_genes), dtype=bool)
        missing_gene_mask[pm_target_columns] = False
        if not np.array_equal(
            candidate[row, missing_gene_mask], before[missing_gene_mask]
        ):
            raise AssertionError(
                f"{target}: a PerturbMatch-missing gene changed before whitening"
            )
        modified_targets.append(target)

    if not np.array_equal(
        candidate[
            [index for index, target in enumerate(targets) if target in pm_missing_targets]
        ],
        baseline[
            [index for index, target in enumerate(targets) if target in pm_missing_targets]
        ],
    ):
        raise AssertionError("PerturbMatch-missing target rows changed")
    memory_checkpoint("effect ensemble", max_rss_gib)
    return baseline, candidate, {
        "targets_total": len(targets),
        "targets_modified": len(modified_targets),
        "targets_preserved_exactly": len(pm_missing_targets),
        "missing_perturbmatch_targets": pm_missing_targets,
        "genes_observed_by_perturbmatch": int(len(pm_target_columns)),
        "genes_preserved_from_incumbent_slot": int(
            len(model.target_genes) - len(pm_target_columns)
        ),
        "incumbent_scale": k_scale,
        "perturbmatch_scale": pm_scale,
    }


def validate_and_load_operator(
    model: MultiDonorModel,
    panel: tuple[str, ...],
    operator_path: Path,
    metadata_path: Path,
) -> tuple[Whitener, np.ndarray, dict[str, object]]:
    union = set()
    for name in model.names:
        union |= set(model.target_genes[model.tcols[name]])
    panel_set = set(panel)
    wcols = np.asarray(
        sorted(
            index
            for index, gene in enumerate(model.target_genes)
            if gene in union and gene not in panel_set
        ),
        dtype=np.int64,
    )
    if any(model.target_genes[index] in panel_set for index in wcols):
        raise AssertionError("a VCC target coordinate entered the dense operator")
    actual = json.loads(metadata_path.read_text())
    expected = {
        "axis": "union",
        "donors": [K562, HEK],
        "n_cols": int(len(wcols)),
        "wcols_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
        "excluded_targets": sorted(panel_set),
    }
    if actual != expected:
        raise ValueError(
            "frozen whitener metadata mismatch:\n"
            + json.dumps({"expected": expected, "actual": actual}, indent=2)
        )
    whitener = Whitener.load(operator_path)
    if whitener.evecs.shape != (len(wcols), len(wcols)):
        raise ValueError(
            f"operator shape {whitener.evecs.shape} != {(len(wcols), len(wcols))}"
        )
    return whitener, wcols, actual


def apply_frozen_pipeline(
    baseline: np.ndarray,
    candidate: np.ndarray,
    whitener: Whitener,
    wcols: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    """Use baseline-derived renormalizers so only the K562 slot can change."""
    baseline_pre = np.asarray(baseline, dtype=np.float64)
    candidate_pre = np.asarray(candidate, dtype=np.float64)
    baseline_before = (
        float(np.median(np.linalg.norm(baseline_pre[:, wcols], axis=1))) or 1.0
    )
    baseline_w = whitener.apply(
        baseline_pre[:, wcols], GAMMA, 0.0, renorm=False
    )
    candidate_w = whitener.apply(
        candidate_pre[:, wcols], GAMMA, 0.0, renorm=False
    )
    baseline_after = float(np.median(np.linalg.norm(baseline_w, axis=1))) or 1.0
    operator_scale = baseline_before / baseline_after
    baseline_post = baseline_pre.copy()
    candidate_post = candidate_pre.copy()
    baseline_post[:, wcols] = baseline_w * operator_scale
    candidate_post[:, wcols] = candidate_w * operator_scale

    alpha_scale = ALPHA * 6.817 / (
        float(np.median(np.linalg.norm(baseline_post, axis=1))) or 1.0
    )
    baseline_post *= alpha_scale
    candidate_post *= alpha_scale
    return baseline_post, candidate_post, {
        "baseline_before_operator_median_norm": baseline_before,
        "baseline_after_operator_median_norm": baseline_after,
        "frozen_operator_scale": operator_scale,
        "frozen_alpha_scale": alpha_scale,
    }


def profile_from_effect(effect: np.ndarray, b_ctrl: np.ndarray) -> np.ndarray:
    expression = np.expm1(b_ctrl)
    cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
    gated = np.where((cpm > MIN_CPM) & (b_ctrl > 0), effect, 0.0)
    profile = np.clip(np.expm1(b_ctrl + gated), 0.0, None) / TS
    total = float(profile.sum())
    if total <= 0:
        raise ValueError("degenerate profile")
    return profile / total


def emit_prediction(
    real_path: Path,
    targets: tuple[str, ...],
    effects: np.ndarray,
    output: Path,
    *,
    max_rss_gib: float,
) -> dict[str, object]:
    real = ad.read_h5ad(real_path)
    matrix = real.X.tocsr() if not sp.isspmatrix_csr(real.X) else real.X
    labels = real.obs[PERT_COL].astype(str).to_numpy()
    genes = np.asarray(real.var_names, dtype=str)
    control_mask = labels == CONTROL
    controls = matrix[control_mask].tocsr()
    total = np.asarray(controls.sum(axis=0)).ravel().astype(np.float64)
    b_ctrl = np.log1p(TS * total / total.sum())
    libraries = np.asarray(controls.sum(axis=1)).ravel()
    if len(targets) != effects.shape[0] or effects.shape[1] != len(genes):
        raise ValueError("effect matrix does not match target/gene axes")

    rng = np.random.default_rng(SEED)
    data: list[np.ndarray] = []
    indices: list[np.ndarray] = []
    indptr = [0]
    obs: list[str] = []
    per_target_totals: dict[str, str] = {}
    for index, target in enumerate(targets):
        profile = profile_from_effect(effects[index], b_ctrl)
        picks = rng.choice(
            len(libraries), 400, replace=len(libraries) < 400
        )
        values, columns, pointers = emit_counts(
            profile,
            libraries[picks],
            0.0,
            rng,
            mode="quasi",
            spread=np.full(len(genes), SPREAD),
        )
        data.append(values)
        indices.append(columns)
        indptr.extend((pointers[1:] + indptr[-1]).tolist())
        obs.extend([target] * 400)
        gene_totals = np.bincount(
            columns, weights=values, minlength=len(genes)
        ).astype("<i8", copy=False)
        per_target_totals[target] = hashlib.sha256(
            gene_totals.tobytes()
        ).hexdigest()
        if index % 10 == 0 or index + 1 == len(targets):
            memory_checkpoint(
                f"count emission {index + 1}/{len(targets)}", max_rss_gib
            )

    data.append(controls.data)
    indices.append(controls.indices)
    indptr.extend((controls.indptr[1:] + indptr[-1]).tolist())
    obs.extend([CONTROL] * controls.shape[0])
    prediction = sp.csr_matrix(
        (
            np.concatenate(data).astype(np.float32),
            np.concatenate(indices).astype(np.int32),
            np.asarray(indptr, dtype=np.int64),
        ),
        shape=(len(obs), len(genes)),
    )
    prediction.sum_duplicates()
    prediction.eliminate_zeros()
    row_totals = np.asarray(prediction.sum(axis=1)).ravel()
    if np.any(row_totals <= 0) or np.any(row_totals > 1_000_000):
        raise ValueError("invalid output row totals")
    if prediction.nnz > 4_750_000_000:
        raise ValueError("output exceeds VCC stored-entry cap")

    output.parent.mkdir(parents=True, exist_ok=True)
    ad.AnnData(
        X=prediction,
        obs=pd.DataFrame(
            {PERT_COL: obs}, index=[f"c{index}" for index in range(len(obs))]
        ),
        var=pd.DataFrame(index=genes),
    ).write_h5ad(output, compression="lzf")
    memory_checkpoint("prediction write", max_rss_gib)
    return {
        "shape": list(prediction.shape),
        "nnz": int(prediction.nnz),
        "control_rows": int(controls.shape[0]),
        "perturbation_rows": int(len(targets) * 400),
        "per_target_column_total_sha256": per_target_totals,
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    panel = read_csv_column(args.panel_csv, "target_gene")
    official = read_csv_column(args.genes_csv, "gene_name")
    if len(panel) != 300:
        raise ValueError(f"expected 300 VCC targets, got {len(panel)}")

    real_backed = ad.read_h5ad(args.real, backed="r")
    try:
        genes = tuple(map(str, real_backed.var_names))
        targets = tuple(
            sorted(set(real_backed.obs[PERT_COL].astype(str)) - {CONTROL})
        )
    finally:
        real_backed.file.close()
    if genes != official:
        raise ValueError("HCT116 gene axis differs from the official VCC axis")
    if not set(targets) <= set(panel):
        raise ValueError("HCT116 targets are not a subset of the VCC panel")

    input_paths = {
        "real": args.real,
        "k562_delta": args.k562 / "delta.npy",
        "hek_delta": args.hek293t / "delta.npy",
        "perturbmatch_delta": args.perturbmatch / "delta.npy",
        "whitener": args.whitener,
        "whitener_metadata": args.whitener_metadata,
        "panel_csv": args.panel_csv,
        "genes_csv": args.genes_csv,
    }
    input_hashes = {name: file_sha256(path) for name, path in input_paths.items()}
    memory_checkpoint("input hashes", args.max_rss_gib)

    donors = {
        K562: Donor.load(args.k562),
        HEK: Donor.load(args.hek293t),
    }
    model = MultiDonorModel(donors, np.asarray(genes), gene_norm=False, eb=False)
    memory_checkpoint("incumbent model", args.max_rss_gib)
    pm = prepare_perturbmatch(
        args.perturbmatch,
        genes,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    baseline_pre, candidate_pre, blend_audit = compute_effects(
        model, pm, targets, max_rss_gib=args.max_rss_gib
    )
    whitener, wcols, operator_metadata = validate_and_load_operator(
        model, panel, args.whitener, args.whitener_metadata
    )
    memory_checkpoint("operator load", args.max_rss_gib)
    baseline_final, candidate_final, scaling = apply_frozen_pipeline(
        baseline_pre, candidate_pre, whitener, wcols
    )

    missing_rows = [
        index
        for index, target in enumerate(targets)
        if target in set(blend_audit["missing_perturbmatch_targets"])
    ]
    if missing_rows and not np.array_equal(
        baseline_final[missing_rows], candidate_final[missing_rows]
    ):
        raise AssertionError(
            "frozen scaling failed to preserve PerturbMatch-missing targets"
        )
    effect_audit = {
        "baseline_pre_sha256": array_sha256(baseline_pre),
        "candidate_pre_sha256": array_sha256(candidate_pre),
        "baseline_final_sha256": array_sha256(baseline_final),
        "candidate_final_sha256": array_sha256(candidate_final),
        "missing_target_final_rows_bitwise_equal": True,
        "changed_pre_entries": int(np.count_nonzero(candidate_pre != baseline_pre)),
        "changed_final_entries": int(np.count_nonzero(candidate_final != baseline_final)),
    }

    del donors, model, pm, whitener, baseline_pre, candidate_pre, baseline_final
    gc.collect()
    memory_checkpoint("before count emission", args.max_rss_gib)
    emission = emit_prediction(
        args.real,
        targets,
        candidate_final,
        args.output,
        max_rss_gib=args.max_rss_gib,
    )
    output_hash = file_sha256(args.output)
    report = {
        "experiment": "061-perturbmatch-hct-count",
        "status": "candidate_built",
        "protocol": {
            "held_out_truth": "Orion HCT116",
            "incumbent_sources": [K562, HEK],
            "candidate_change": (
                "replace K562 slot by 50:50 independently normalized raw-effect "
                "ensemble of incumbent K562-GWPS and PerturbMatch"
            ),
            "candidate_independent_operator": True,
            "operator_fit_sources": [K562, HEK],
            "all_300_vcc_target_coordinates_excluded": True,
            "gamma": GAMMA,
            "alpha": ALPHA,
            "cpm_gate": MIN_CPM,
            "emission": "quasi exact-total",
            "spread": SPREAD,
            "seed": SEED,
            "submission": False,
        },
        "axis": {
            "hct_targets": len(targets),
            "official_genes": len(official),
            "operator_genes": len(wcols),
            "operator_wcols_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
        },
        "blend_audit": blend_audit,
        "operator_metadata": operator_metadata,
        "scaling": scaling,
        "effect_audit": effect_audit,
        "emission": emission,
        "inputs_sha256": input_hashes,
        "candidate_h5ad": {
            "path": str(args.output),
            "bytes": args.output.stat().st_size,
            "sha256": output_hash,
        },
        "peak_rss_gib": peak_rss_gib(),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True), flush=True)
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--real", type=Path, required=True)
    parser.add_argument("--k562", type=Path, required=True)
    parser.add_argument("--hek293t", type=Path, required=True)
    parser.add_argument("--perturbmatch", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--whitener", type=Path, required=True)
    parser.add_argument("--whitener-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=28.0)
    return parser


if __name__ == "__main__":
    run(build_parser().parse_args())

"""060 — frozen CD4 rule under production-style union whitening, held K562.

Fit the whitening operator on Orion HCT116 + HEK293T only, transform the
two-source baseline and the frozen experiment-047 CD4 candidate, and then
evaluate once against untouched centered K562-GWPS effects.

This is an effect-level PDS experiment. It emits no counts, runs no scorer, and
cannot submit anything. The full 300-gene VCC panel is excluded before the
operator fit, transformation, and retrieval.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import resource
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

_BASE_SPEC = importlib.util.spec_from_file_location(
    "exp054_run", ROOT / "experiments/054-cd4-union-hct/run.py"
)
if _BASE_SPEC is None or _BASE_SPEC.loader is None:
    raise RuntimeError("cannot load experiment 054 support module")
BASE = importlib.util.module_from_spec(_BASE_SPEC)
sys.modules[_BASE_SPEC.name] = BASE
_BASE_SPEC.loader.exec_module(BASE)

SOURCE_NAMES = ("orion_hct116", "orion_hek293t")
CD4_NAMES = BASE.CD4_NAMES
EFFECT_SCALE = BASE.EFFECT_SCALE
BLEND_FRACTION = 0.50
WHITEN_GAMMA = 0.65


def peak_rss_gib() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Sequence[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def enforce_memory_limit(limit_gib: float, stage: str) -> None:
    if limit_gib > 0 and peak_rss_gib() >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak_rss_gib():.3f} GiB reached "
            f"{limit_gib:.3f} GiB"
        )


@dataclass(frozen=True)
class Operator:
    evals: np.ndarray
    evecs: np.ndarray
    metadata: dict[str, object]
    path: Path
    metadata_path: Path


def _last_source_columns(
    bank: BASE.Bank,
    official_genes: Sequence[str],
    selected_official_positions: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Map selected official positions to the donor's last duplicate column."""
    mapping = BASE.map_bank_to_axis(bank, official_genes)
    selected_to_axis = {
        int(position): index
        for index, position in enumerate(selected_official_positions.tolist())
    }
    source_columns: list[int] = []
    operator_columns: list[int] = []
    for source, destination in zip(
        mapping.unique_source_columns,
        mapping.unique_target_positions,
        strict=True,
    ):
        operator_column = selected_to_axis.get(int(destination))
        if operator_column is not None:
            source_columns.append(int(source))
            operator_columns.append(operator_column)
    return (
        np.asarray(source_columns, dtype=np.int64),
        np.asarray(operator_columns, dtype=np.int64),
    )


def fit_union_whitener(
    source_banks: Sequence[BASE.Bank],
    prepared_sources: Sequence[BASE.PreparedSource],
    official_genes: Sequence[str],
    panel: Sequence[str],
    output_path: Path,
    scratch: Path,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[Operator, np.ndarray]:
    """Fit the production covariance on source rows only, using a disk-backed basis."""
    if tuple(bank.name for bank in source_banks) != SOURCE_NAMES:
        raise ValueError(
            f"source order must be exactly {SOURCE_NAMES}, got "
            f"{tuple(bank.name for bank in source_banks)}"
        )
    panel_set = set(panel)
    donor_union = set().union(*(set(bank.genes) for bank in source_banks))
    wcols = np.asarray(
        [
            index
            for index, gene in enumerate(official_genes)
            if gene in donor_union and gene not in panel_set
        ],
        dtype=np.int64,
    )
    if any(official_genes[index] in panel_set for index in wcols):
        raise AssertionError("VCC target coordinate entered operator axis")

    metadata: dict[str, object] = {
        "axis": "union",
        "donors": list(SOURCE_NAMES),
        "n_cols": int(len(wcols)),
        "wcols_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
        "excluded_targets": sorted(panel_set),
        "basis_rows_by_source": {
            bank.name: int(bank.shape[0]) for bank in source_banks
        },
        "basis_normalization": (
            "production equal-source renorm: median source norm / donor norm"
        ),
        "fit_uses_held_k562": False,
        "gamma_selected_before_k562": WHITEN_GAMMA,
    }
    metadata_path = output_path.with_suffix(".json")
    if output_path.exists() and metadata_path.exists():
        actual = json.loads(metadata_path.read_text())
        if actual != metadata:
            raise ValueError("existing whitener metadata does not match frozen protocol")
        archive = np.load(output_path)
        evals = np.asarray(archive["evals"], dtype=np.float64)
        evecs = np.asarray(archive["evecs"], dtype=np.float64)
        return Operator(evals, evecs, actual, output_path, metadata_path), wcols

    scratch.mkdir(parents=True, exist_ok=True)
    total_rows = sum(bank.shape[0] for bank in source_banks)
    basis_path = scratch / "source_basis.npy"
    basis = np.lib.format.open_memmap(
        basis_path,
        mode="w+",
        dtype=np.float64,
        shape=(total_rows, len(wcols)),
    )
    source_scale = {
        source.name: source.raw_median_norm for source in prepared_sources
    }
    reference_norm = float(np.median(list(source_scale.values()))) or 1.0

    row_offset = 0
    for bank in source_banks:
        source_columns, operator_columns = _last_source_columns(
            bank, official_genes, wcols
        )
        mean = BASE._stream_column_mean(
            bank,
            source_columns,
            chunk_rows=chunk_rows,
            max_rss_gib=max_rss_gib,
        )
        delta = np.load(bank.delta_path, mmap_mode="r")
        multiplier = reference_norm / source_scale[bank.name]
        for lo in range(0, bank.shape[0], chunk_rows):
            hi = min(lo + chunk_rows, bank.shape[0])
            destination = basis[row_offset + lo : row_offset + hi]
            destination[:] = 0.0
            destination[:, operator_columns] = (
                np.asarray(delta[lo:hi, source_columns], dtype=np.float64) - mean
            ) * multiplier
            enforce_memory_limit(max_rss_gib, f"{bank.name} basis construction")
        row_offset += bank.shape[0]
        basis.flush()

    print(
        f"fitting source-only union whitener basis={basis.shape} "
        f"peak_rss={peak_rss_gib():.2f} GiB",
        flush=True,
    )
    covariance = (basis.T @ basis) / total_rows
    del basis
    evals, evecs = np.linalg.eigh(covariance)
    del covariance
    evals = np.maximum(evals, float(evals.max()) * 1e-12)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        output_path,
        evals=np.asarray(evals, dtype=np.float64),
        evecs=np.asarray(evecs, dtype=np.float32),
    )
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    try:
        basis_path.unlink()
    except FileNotFoundError:
        pass
    enforce_memory_limit(max_rss_gib, "whitener fit")
    return (
        Operator(evals, evecs, metadata, output_path, metadata_path),
        wcols,
    )


def apply_operator(
    matrix: np.ndarray,
    operator: Operator,
    wcols: np.ndarray,
) -> tuple[np.ndarray, dict[str, float]]:
    """Apply gamma=0.65 and production preserve-axis median renormalization."""
    output = np.asarray(matrix, dtype=np.float64).copy()
    before = float(np.median(np.linalg.norm(output[:, wcols], axis=1))) or 1.0
    weights = operator.evals ** (-WHITEN_GAMMA / 2.0)
    transformed = ((output[:, wcols] @ operator.evecs) * weights) @ operator.evecs.T
    after = float(np.median(np.linalg.norm(transformed, axis=1))) or 1.0
    scale = before / after
    output[:, wcols] = transformed * scale
    return output, {
        "before_operator_median_norm": before,
        "after_operator_median_norm": after,
        "applied_scale": scale,
    }


def load_k562_truth(
    bank: BASE.Bank,
    targets: Sequence[str],
    official_genes: Sequence[str],
    retrieval_columns: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> np.ndarray:
    """Load untouched centered K562 effects only after model fit/transformation."""
    source_columns, retrieval_positions = _last_source_columns(
        bank, official_genes, retrieval_columns
    )
    if len(source_columns) != len(retrieval_columns):
        raise ValueError("retrieval axis contains genes absent from K562 truth")
    if not np.array_equal(retrieval_positions, np.arange(len(retrieval_columns))):
        raise ValueError("K562 retrieval mapping is not complete and ordered")
    mean = BASE._stream_column_mean(
        bank,
        source_columns,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    rows = np.asarray([bank.pert_pos[target] for target in targets], dtype=np.int64)
    delta = np.load(bank.delta_path, mmap_mode="r")
    truth = np.asarray(delta[np.ix_(rows, source_columns)], dtype=np.float64) - mean
    enforce_memory_limit(max_rss_gib, "K562 truth load")
    return truth


def pds_per_target(prediction: np.ndarray, truth: np.ndarray) -> np.ndarray:
    prediction = np.asarray(prediction, dtype=np.float64)
    truth = np.asarray(truth, dtype=np.float64)
    prediction /= np.maximum(np.linalg.norm(prediction, axis=1, keepdims=True), 1e-12)
    truth /= np.maximum(np.linalg.norm(truth, axis=1, keepdims=True), 1e-12)
    similarities = prediction @ truth.T
    own = np.diag(similarities)[:, None]
    n = similarities.shape[0]
    ranks = (similarities > own).sum(axis=1)
    ranks = ranks + 0.5 * ((similarities == own).sum(axis=1) - 1)
    return 1.0 - ranks / (n - 1)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> dict[str, object]:
    panel = BASE.read_csv_column(args.panel_csv, "target_gene")
    official = BASE.read_csv_column(args.genes_csv, "gene_name")
    panel_set = set(panel)
    if len(panel_set) != args.expected_panel_targets:
        raise ValueError(
            f"expected {args.expected_panel_targets} panel targets, "
            f"got {len(panel_set)}"
        )

    source_banks = (
        BASE.Bank.load("orion_hct116", args.hct116),
        BASE.Bank.load("orion_hek293t", args.hek293t),
    )
    held_bank = BASE.Bank.load("k562_gwps", args.k562)
    cd4_banks = tuple(
        BASE.Bank.load(name, args.cd4_root / name, require_cd4_scale=True)
        for name in CD4_NAMES
    )
    all_banks = (*source_banks, held_bank, *cd4_banks)
    targets = tuple(
        target
        for target in panel
        if all(target in bank.pert_pos for bank in all_banks)
    )
    if len(targets) < args.min_shared_targets:
        raise ValueError(
            f"only {len(targets)} targets are shared by all six banks; "
            f"require {args.min_shared_targets}"
        )

    prepared_sources = [
        BASE.prepare_main_source(
            bank,
            targets,
            official,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        for bank in source_banks
    ]
    baseline_pre, source_votes = BASE.build_current_baseline(prepared_sources)

    primary_cd4 = BASE.prepare_cd4_mean(
        cd4_banks,
        targets,
        official,
        panel_set,
        "gene_norm",
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    if not np.all(primary_cd4.target_condition_count == len(CD4_NAMES)):
        raise ValueError("full shared gallery unexpectedly lacks a CD4 condition")
    primary_pre = BASE.blend_cd4(
        baseline_pre, primary_cd4, "preserve_baseline"
    )
    candidate_pre: dict[str, np.ndarray] = {
        "cd4_mean3_f0.5_gene_norm": primary_pre
    }
    if args.include_raw_diagnostic:
        raw_cd4 = BASE.prepare_cd4_mean(
            cd4_banks,
            targets,
            official,
            panel_set,
            "raw",
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        candidate_pre["cd4_mean3_f0.5_raw_diagnostic"] = BASE.blend_cd4(
            baseline_pre, raw_cd4, "preserve_baseline"
        )

    operator, wcols = fit_union_whitener(
        source_banks,
        prepared_sources,
        official,
        panel,
        args.whitener,
        args.scratch,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    if len(wcols) < args.min_operator_genes:
        raise ValueError(
            f"only {len(wcols)} operator genes; require {args.min_operator_genes}"
        )
    baseline_post, baseline_whitening = apply_operator(
        baseline_pre, operator, wcols
    )
    candidate_post = {
        name: apply_operator(matrix, operator, wcols)
        for name, matrix in candidate_pre.items()
    }

    # K562 values are deliberately first touched after both source-only fitting
    # and all prediction transforms are complete.
    held_gene_set = set(held_bank.genes)
    retrieval_columns = np.asarray(
        [index for index in wcols if official[index] in held_gene_set],
        dtype=np.int64,
    )
    if any(official[index] in panel_set for index in retrieval_columns):
        raise AssertionError("VCC target coordinate entered retrieval axis")
    if len(retrieval_columns) < args.min_retrieval_genes:
        raise ValueError(
            f"only {len(retrieval_columns)} K562-observed retrieval genes; "
            f"require {args.min_retrieval_genes}"
        )
    truth = load_k562_truth(
        held_bank,
        targets,
        official,
        retrieval_columns,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )

    baseline_per_target = pds_per_target(
        baseline_post[:, retrieval_columns], truth
    )
    baseline_pds = float(np.mean(baseline_per_target))
    rows = [
        {
            "candidate": "baseline_hct116_hek293t",
            "is_primary": False,
            "pds": baseline_pds,
            "delta_vs_baseline": 0.0,
        }
    ]
    per_target_rows: list[dict[str, object]] = []
    for target, value in zip(targets, baseline_per_target, strict=True):
        per_target_rows.append(
            {
                "target_gene": target,
                "candidate": "baseline_hct116_hek293t",
                "pds": float(value),
            }
        )

    candidates: list[dict[str, object]] = []
    for name, (matrix, whitening) in candidate_post.items():
        per_target = pds_per_target(matrix[:, retrieval_columns], truth)
        score = float(np.mean(per_target))
        item = {
            "candidate": name,
            "is_primary": name == "cd4_mean3_f0.5_gene_norm",
            "pds": score,
            "delta_vs_baseline": score - baseline_pds,
            "whitening": whitening,
        }
        candidates.append(item)
        rows.append(
            {
                "candidate": name,
                "is_primary": item["is_primary"],
                "pds": score,
                "delta_vs_baseline": item["delta_vs_baseline"],
            }
        )
        for target, value in zip(targets, per_target, strict=True):
            per_target_rows.append(
                {
                    "target_gene": target,
                    "candidate": name,
                    "pds": float(value),
                }
            )

    pre_whiten_diagnostic = {
        "baseline_pds": float(
            np.mean(pds_per_target(baseline_pre[:, retrieval_columns], truth))
        ),
        "candidates": {
            name: float(
                np.mean(pds_per_target(matrix[:, retrieval_columns], truth))
            )
            for name, matrix in candidate_pre.items()
        },
    }
    result: dict[str, object] = {
        "experiment": "060-cd4-union-k562",
        "status": "completed",
        "protocol": {
            "metric": "effect-level rank PDS",
            "held_out_truth": "centered raw K562-GWPS effects",
            "source_donors": list(SOURCE_NAMES),
            "target_gallery": (
                "VCC panel intersection across HCT116, HEK293T, K562, "
                "CD4 rest, CD4 stim8hr, and CD4 stim48hr"
            ),
            "gallery_size": len(targets),
            "whitening_axis": "source-donor union on official VCC axis",
            "whitening_gamma": WHITEN_GAMMA,
            "whitening_renormalization": (
                "production preserve: candidate-specific pre/post operator "
                "median norm on the operator axis"
            ),
            "cd4_rule": (
                "mean3 centered mean-absolute gene normalization, "
                "fraction 0.50, preserve baseline outside CD4 gene coverage"
            ),
            "selection_on_k562": False,
            "count_emission": False,
            "cell_eval2": False,
            "submission": False,
        },
        "inputs": {
            "hct116": str(args.hct116.resolve()),
            "hek293t": str(args.hek293t.resolve()),
            "k562_truth": str(args.k562.resolve()),
            "cd4_root": str(args.cd4_root.resolve()),
            "panel_csv": str(args.panel_csv.resolve()),
            "genes_csv": str(args.genes_csv.resolve()),
            "k562_delta_sha256": sha256_file(held_bank.delta_path),
        },
        "leakage_audit": {
            "held_k562_in_fit_sources": False,
            "held_k562_values_loaded_after_fit_and_transform": True,
            "hyperparameters_selected_before_k562": {
                "gamma": WHITEN_GAMMA,
                "cd4_fraction": BLEND_FRACTION,
                "cd4_mode": "gene_norm",
            },
            "challenge_targets_excluded_count": len(panel_set),
            "challenge_targets_in_operator_axis": len(
                panel_set & {official[index] for index in wcols}
            ),
            "challenge_targets_in_retrieval_axis": len(
                panel_set & {official[index] for index in retrieval_columns}
            ),
            "truth_transformed": False,
            "truth_normalized_gene_wise": False,
            "truth_used_for_candidate_selection": False,
            "full_gallery_used_for_every_query": len(targets),
        },
        "axis_audit": {
            "official_genes": len(official),
            "challenge_targets": len(panel),
            "operator_genes": len(wcols),
            "retrieval_genes": len(retrieval_columns),
            "panel_order_sha256": sha256_lines(panel),
            "target_order_sha256": sha256_lines(targets),
            "operator_columns_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
            "retrieval_columns_sha256": hashlib.sha256(
                retrieval_columns.tobytes()
            ).hexdigest(),
        },
        "whitener": {
            "path": str(operator.path.resolve()),
            "metadata_path": str(operator.metadata_path.resolve()),
            "npz_sha256": sha256_file(operator.path),
            "metadata": operator.metadata,
        },
        "baseline": {
            "pds": baseline_pds,
            "source_vote_counts": {
                str(vote): int(np.sum(source_votes == vote))
                for vote in sorted(set(source_votes.tolist()))
            },
            "source_details": [
                {
                    "name": source.name,
                    "raw_median_centered_norm": source.raw_median_norm,
                    "mapped_genes": source.mapped_gene_count,
                    "duplicate_symbols": list(source.duplicate_symbols),
                }
                for source in prepared_sources
            ],
            "whitening": baseline_whitening,
        },
        "candidates": candidates,
        "pre_whiten_diagnostic": pre_whiten_diagnostic,
        "infrastructure": {
            "provider": "AWS EC2",
            "instance_id": args.instance_id,
            "instance_type": args.instance_type,
            "s3_prefix": args.s3_prefix,
            "cleanup_needed": (
                "terminate the experiment worker after verified result download; "
                "retain the unique S3 output prefix until project-wide cleanup"
            ),
        },
        "peak_rss_gib": peak_rss_gib(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    _write_csv(args.output.with_name("summary.csv"), rows)
    _write_csv(args.output.with_name("per_target_pds.csv"), per_target_rows)
    shutil.rmtree(args.scratch, ignore_errors=True)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hct116", type=Path, required=True)
    parser.add_argument("--hek293t", type=Path, required=True)
    parser.add_argument("--k562", type=Path, required=True)
    parser.add_argument("--cd4-root", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    parser.add_argument(
        "--whitener", type=Path, default=HERE / "whitener_union.npz"
    )
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=55.0)
    parser.add_argument("--expected-panel-targets", type=int, default=300)
    parser.add_argument("--min-shared-targets", type=int, default=240)
    parser.add_argument("--min-operator-genes", type=int, default=17_000)
    parser.add_argument("--min-retrieval-genes", type=int, default=7_000)
    parser.add_argument(
        "--include-raw-diagnostic",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--instance-id", default="unknown")
    parser.add_argument("--instance-type", default="unknown")
    parser.add_argument("--s3-prefix", default="unknown")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.chunk_rows < 1:
        raise SystemExit("--chunk-rows must be positive")
    result = run(args)
    primary = next(item for item in result["candidates"] if item["is_primary"])
    print(
        json.dumps(
            {
                "baseline_pds": result["baseline"]["pds"],
                "primary_pds": primary["pds"],
                "primary_delta": primary["delta_vs_baseline"],
                "peak_rss_gib": result["peak_rss_gib"],
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()

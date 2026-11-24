#!/usr/bin/env python3
"""Leak-aware effect-level LODO test of KOLF2.1J as an added donor."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import resource
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
MAIN_DEFAULTS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}
DEVELOPMENT_HELD = ("orion_hct116", "orion_hek293t")
FROZEN_HELD = "k562_gwps"


def peak_rss_gib() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / 1024**3 if os.uname().sysname == "Darwin" else rss / 1024**2


def enforce_memory_limit(limit_gib: float) -> None:
    if limit_gib > 0 and peak_rss_gib() > limit_gib:
        raise MemoryError(
            f"peak RSS {peak_rss_gib():.2f} GiB exceeded {limit_gib:.2f} GiB"
        )


def read_lines(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines() if line]


@dataclass(frozen=True)
class Bank:
    name: str
    path: Path
    genes: tuple[str, ...]
    perts: tuple[str, ...]
    gene_pos: dict[str, int]
    pert_pos: dict[str, int]
    delta_path: Path
    shape: tuple[int, int]
    duplicate_gene_symbols: int

    @classmethod
    def load(cls, name: str, path: Path, *, require_kolf: bool = False) -> "Bank":
        genes = tuple(read_lines(path / "genes.txt"))
        perts = tuple(read_lines(path / "perts.txt"))
        delta_path = path / "delta.npy"
        delta = np.load(delta_path, mmap_mode="r")
        if delta.shape != (len(perts), len(genes)):
            raise ValueError(
                f"{name}: delta shape {delta.shape} != {(len(perts), len(genes))}"
            )
        if not np.issubdtype(delta.dtype, np.floating):
            raise TypeError(f"{name}: delta dtype must be floating point")
        if len(perts) != len(set(perts)):
            raise ValueError(f"{name}: duplicate perturbations")

        b_ctrl_path = path / "b_ctrl.npy"
        b_ctrl = np.load(b_ctrl_path, mmap_mode="r")
        if b_ctrl.shape != (len(genes),):
            raise ValueError(f"{name}: malformed b_ctrl.npy")
        gene_pos: dict[str, int] = {}
        duplicate_symbols: set[str] = set()
        for index, gene in enumerate(genes):
            previous = gene_pos.get(gene)
            if previous is None:
                gene_pos[gene] = index
                continue
            duplicate_symbols.add(gene)
            if float(b_ctrl[index]) > float(b_ctrl[previous]):
                gene_pos[gene] = index

        if require_kolf:
            meta = json.loads((path / "meta.json").read_text())
            expected = "delta = int8 source value * 0.04"
            if expected not in meta.get("effect_definition", ""):
                raise ValueError(f"{name}: unrecognized KOLF effect scale")
            if meta.get("vcc_target_coverage") != len(perts):
                raise ValueError(f"{name}: target coverage does not match labels")
            if duplicate_symbols:
                raise ValueError(f"{name}: compact KOLF axis must be unique")
        return cls(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            gene_pos=gene_pos,
            pert_pos={pert: index for index, pert in enumerate(perts)},
            delta_path=delta_path,
            shape=tuple(map(int, delta.shape)),
            duplicate_gene_symbols=len(duplicate_symbols),
        )


def streaming_center_and_scale(
    bank: Bank,
    columns: np.ndarray,
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute gene-wise center and mean-absolute scale without loading a bank."""
    delta = np.load(bank.delta_path, mmap_mode="r")
    mean = np.zeros(len(columns), dtype=np.float64)
    for lo in range(0, bank.shape[0], chunk_rows):
        hi = min(lo + chunk_rows, bank.shape[0])
        mean += np.asarray(delta[lo:hi, columns], dtype=np.float64).sum(axis=0)
        enforce_memory_limit(max_rss_gib)
    mean /= bank.shape[0]

    mean_abs = np.zeros(len(columns), dtype=np.float64)
    for lo in range(0, bank.shape[0], chunk_rows):
        hi = min(lo + chunk_rows, bank.shape[0])
        block = np.asarray(delta[lo:hi, columns], dtype=np.float64)
        mean_abs += np.abs(block - mean).sum(axis=0)
        enforce_memory_limit(max_rss_gib)
    mean_abs /= bank.shape[0]
    live = mean_abs[mean_abs > 0]
    floor = float(np.percentile(live, 25.0)) if live.size else 1.0
    return mean, np.maximum(mean_abs, floor)


def selected_effects(
    bank: Bank,
    perturbations: list[str],
    features: list[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
) -> tuple[np.ndarray, np.ndarray]:
    rows = np.asarray([bank.pert_pos[pert] for pert in perturbations], dtype=np.int64)
    columns = np.asarray([bank.gene_pos[gene] for gene in features], dtype=np.int64)
    mean, scale = streaming_center_and_scale(
        bank,
        columns,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    delta = np.load(bank.delta_path, mmap_mode="r")
    raw = np.asarray(delta[np.ix_(rows, columns)], dtype=np.float64) - mean
    normalized = raw / scale
    raw_norm = float(np.median(np.linalg.norm(raw, axis=1))) or 1.0
    normalized_norm = float(np.median(np.linalg.norm(normalized, axis=1))) or 1.0
    normalized *= raw_norm / normalized_norm
    enforce_memory_limit(max_rss_gib)
    return raw, normalized


def row_normalize(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.float64)
    return matrix / (np.linalg.norm(matrix, axis=1, keepdims=True) + 1e-12)


def pds(prediction: np.ndarray, truth: np.ndarray) -> float:
    """PDS over the complete shared-target retrieval gallery."""
    if prediction.shape != truth.shape:
        raise ValueError("prediction and truth must have identical shapes")
    scores = row_normalize(prediction) @ row_normalize(truth).T
    own = np.diag(scores)[:, None]
    n = scores.shape[1]
    if n < 2:
        raise ValueError("PDS requires at least two gallery targets")
    rank = (scores > own).sum(axis=1) + 0.5 * ((scores == own).sum(axis=1) - 1)
    return float(np.mean(1.0 - rank / (n - 1)))


def scale_to_reference(matrix: np.ndarray, reference_norm: float) -> np.ndarray:
    norm = float(np.median(np.linalg.norm(matrix, axis=1))) or 1.0
    return matrix * (reference_norm / norm)


def make_baseline(parts: dict[str, np.ndarray]) -> tuple[np.ndarray, float]:
    norms = [np.median(np.linalg.norm(value, axis=1)) for value in parts.values()]
    reference_norm = float(np.median(norms)) or 1.0
    scaled = [scale_to_reference(matrix, reference_norm) for matrix in parts.values()]
    return sum(scaled) / len(scaled), reference_norm


def blend_candidate(
    baseline: np.ndarray,
    auxiliary: np.ndarray,
    fraction: float,
    reference_norm: float,
) -> np.ndarray:
    auxiliary = scale_to_reference(auxiliary, reference_norm)
    return (1.0 - fraction) * baseline + fraction * auxiliary


def choose_development_candidate(
    rows: list[dict[str, object]],
) -> dict[str, object]:
    grouped: dict[tuple[str, float], list[float]] = {}
    for row in rows:
        key = (str(row["mode"]), float(row["fraction"]))
        grouped.setdefault(key, []).append(float(row["delta_vs_baseline"]))
    if any(len(values) != len(DEVELOPMENT_HELD) for values in grouped.values()):
        raise ValueError("development candidates do not have both held folds")
    ranked = [
        {
            "mode": mode,
            "fraction": fraction,
            "development_mean_delta": float(np.mean(values)),
            "development_min_delta": float(np.min(values)),
        }
        for (mode, fraction), values in grouped.items()
    ]
    # No K562 score exists yet. Prefer mean, then worst fold, then less KOLF.
    return max(
        ranked,
        key=lambda row: (
            row["development_mean_delta"],
            row["development_min_delta"],
            -row["fraction"],
            row["mode"] == "raw",
        ),
    )


def evaluate_development(
    effects: dict[str, dict[str, np.ndarray]],
    fractions: list[float],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    rows: list[dict[str, object]] = []
    for mode in ("raw", "gene_norm"):
        for held in DEVELOPMENT_HELD:
            source_names = [
                name
                for name in ("k562_gwps", "orion_hct116", "orion_hek293t")
                if name != held
            ]
            baseline, reference_norm = make_baseline(
                {name: effects[name][mode] for name in source_names}
            )
            truth = effects[held]["raw"]
            baseline_score = pds(baseline, truth)
            for fraction in fractions:
                prediction = blend_candidate(
                    baseline,
                    effects["kolf2"][mode],
                    fraction,
                    reference_norm,
                )
                score = pds(prediction, truth)
                rows.append(
                    {
                        "mode": mode,
                        "held": held,
                        "fraction": fraction,
                        "pds": score,
                        "baseline_pds": baseline_score,
                        "delta_vs_baseline": score - baseline_score,
                    }
                )
    return rows, choose_development_candidate(rows)


def evaluate_frozen_k562(
    effects: dict[str, dict[str, np.ndarray]],
    selected: dict[str, object],
) -> dict[str, object]:
    """Evaluate exactly one frozen development-selected K562 candidate."""
    mode = str(selected["mode"])
    fraction = float(selected["fraction"])
    baseline, reference_norm = make_baseline(
        {name: effects[name][mode] for name in ("orion_hct116", "orion_hek293t")}
    )
    truth = effects[FROZEN_HELD]["raw"]
    baseline_score = pds(baseline, truth)
    prediction = blend_candidate(
        baseline,
        effects["kolf2"][mode],
        fraction,
        reference_norm,
    )
    score = pds(prediction, truth)
    return {
        "mode": mode,
        "held": FROZEN_HELD,
        "fraction": fraction,
        "pds": score,
        "baseline_pds": baseline_score,
        "delta_vs_baseline": score - baseline_score,
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    panel = (
        pd.read_csv(args.panel_csv)["target_gene"]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )
    official = (
        pd.read_csv(args.genes_csv)["gene_name"]
        .dropna()
        .astype(str)
        .str.strip()
        .tolist()
    )
    panel_set = set(panel)
    if len(panel_set) != len(panel):
        raise ValueError("VCC panel contains duplicate symbols")

    paths = {
        **MAIN_DEFAULTS,
        "kolf2": Path(args.kolf),
    }
    banks = {
        name: Bank.load(name, path, require_kolf=name == "kolf2")
        for name, path in paths.items()
    }
    perturbations = [
        pert for pert in panel if all(pert in bank.pert_pos for bank in banks.values())
    ]
    if len(perturbations) < args.min_shared_targets:
        raise ValueError(
            f"only {len(perturbations)} panel targets shared by every bank; "
            f"require {args.min_shared_targets}"
        )
    features = [
        gene
        for gene in official
        if gene not in panel_set
        and all(gene in bank.gene_pos for bank in banks.values())
    ]
    if panel_set & set(features):
        raise AssertionError("target-coordinate leakage into feature axis")
    if len(features) < args.min_shared_genes:
        raise ValueError(
            f"only {len(features)} non-panel genes shared by every bank; "
            f"require {args.min_shared_genes}"
        )

    effects: dict[str, dict[str, np.ndarray]] = {}
    for name, bank in banks.items():
        raw, gene_norm = selected_effects(
            bank,
            perturbations,
            features,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        effects[name] = {"raw": raw, "gene_norm": gene_norm}
        print(
            f"prepared {name:16} rows={len(perturbations):3d} "
            f"genes={len(features):5d} peak_rss={peak_rss_gib():.2f} GiB",
            flush=True,
        )

    fractions = sorted(set(args.kolf_fractions))
    development_rows, selected = evaluate_development(effects, fractions)
    # The K562 target-role holdout is not computed until the development rule
    # has been frozen. Only that one selected candidate is evaluated.
    frozen = evaluate_frozen_k562(effects, selected)
    three_fold_macro_delta = float(
        np.mean(
            [
                selected["development_mean_delta"],
                selected["development_mean_delta"],
                frozen["delta_vs_baseline"],
            ]
        )
    )
    selected = {
        **selected,
        "selection_rule": (
            "maximize mean delta on held HCT116 and HEK293T only; "
            "tie-break by worst fold, then smaller fraction, then raw"
        ),
        "k562_pds": frozen["pds"],
        "k562_baseline_pds": frozen["baseline_pds"],
        "k562_delta": frozen["delta_vs_baseline"],
        "three_fold_macro_delta": three_fold_macro_delta,
        "passes_primary_gate": (
            float(frozen["delta_vs_baseline"]) >= 0.003
            and three_fold_macro_delta > 0
            and float(selected["development_min_delta"]) >= -0.002
        ),
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(development_rows[0]))
        writer.writeheader()
        writer.writerows(development_rows)
    result: dict[str, object] = {
        "protocol": {
            "perturbations": len(perturbations),
            "retrieval_gallery_size": len(perturbations),
            "queries_per_fold": len(perturbations),
            "genes": len(features),
            "perturbation_policy": "VCC panel intersection across all four banks",
            "gene_policy": (
                "official VCC-axis intersection across all banks minus all "
                "300 VCC target coordinates"
            ),
            "truth_space": "centered raw held-donor effects",
            "prediction_modes": [
                "centered raw",
                "centered mean-absolute gene normalized",
            ],
            "fractions": fractions,
            "development_folds": list(DEVELOPMENT_HELD),
            "frozen_holdout": FROZEN_HELD,
            "k562_candidate_evaluations_after_selection": 1,
            "no_whitening": True,
            "no_cell_eval2": True,
            "no_count_emission": True,
            "duplicate_gene_symbols_resolved_by_max_b_ctrl": {
                name: bank.duplicate_gene_symbols for name, bank in banks.items()
            },
        },
        "perturbations": perturbations,
        "feature_axis_sha256": hashlib.sha256(
            ("\n".join(features) + "\n").encode()
        ).hexdigest(),
        "development_rows": development_rows,
        "selected_on_hct_only": selected,
        "frozen_k562": frozen,
        "peak_rss_gib": peak_rss_gib(),
    }
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {output} and {output.with_suffix('.csv')}", flush=True)
    print(json.dumps(selected, indent=2), flush=True)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--kolf",
        default=str(ROOT / "data/external/kolf2"),
    )
    parser.add_argument(
        "--panel-csv",
        default=str(ROOT / "data/challenge/pert_counts.csv"),
    )
    parser.add_argument(
        "--genes-csv",
        default=str(ROOT / "data/challenge/gene_names.csv"),
    )
    parser.add_argument(
        "--kolf-fractions",
        nargs="+",
        type=float,
        default=[0.1, 0.25, 0.5],
    )
    parser.add_argument("--chunk-rows", type=int, default=64)
    parser.add_argument("--max-rss-gib", type=float, default=3.5)
    parser.add_argument("--min-shared-targets", type=int, default=240)
    parser.add_argument("--min-shared-genes", type=int, default=7_000)
    parser.add_argument(
        "--output",
        default=str(HERE / "results.json"),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.chunk_rows < 1:
        raise SystemExit("--chunk-rows must be positive")
    if not 0 < args.max_rss_gib <= 4:
        raise SystemExit("--max-rss-gib must be in (0, 4]")
    if sorted(set(args.kolf_fractions)) != [0.1, 0.25, 0.5]:
        raise SystemExit("--kolf-fractions must be exactly 0.1 0.25 0.5")
    run(args)


if __name__ == "__main__":
    main()

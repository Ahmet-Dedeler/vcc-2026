"""Memory-safe donor-held-out evaluation of direct primary-CD4 transfer.

The experiment compares the existing two-donor consensus with adding one
primary-CD4 condition as an extra donor.  It evaluates on VCC-panel
perturbations shared by every source and on a fixed official gene axis that
excludes *all* 300 panel target coordinates before centering, scaling, or
combining effects.

This is an effect-level PDS screen.  It does not run cell-eval2, generate
single-cell counts, fit a covariance matrix, or make a submission.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import resource
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

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
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")


def peak_rss_gib() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / 1024**3 if os.uname().sysname == "Darwin" else rss / 1024**2


def enforce_memory_limit(limit_gib: float) -> None:
    if limit_gib > 0 and peak_rss_gib() > limit_gib:
        raise MemoryError(
            f"peak RSS {peak_rss_gib():.2f} GiB exceeded limit {limit_gib:.2f} GiB"
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
    def load(cls, name: str, path: Path, *, require_cd4_scale: bool = False) -> "Bank":
        genes = tuple(read_lines(path / "genes.txt"))
        perts = tuple(read_lines(path / "perts.txt"))
        delta_path = path / "delta.npy"
        delta = np.load(delta_path, mmap_mode="r")
        if delta.shape != (len(perts), len(genes)):
            raise ValueError(
                f"{name}: delta shape {delta.shape} != labels {(len(perts), len(genes))}"
            )
        if not np.issubdtype(delta.dtype, np.floating):
            raise TypeError(f"{name}: delta dtype {delta.dtype} is not floating point")
        if len(perts) != len(set(perts)):
            raise ValueError(f"{name}: duplicate perturbations")
        b_ctrl = np.load(path / "b_ctrl.npy", mmap_mode="r")
        if b_ctrl.shape != (len(genes),):
            raise ValueError(
                f"{name}: b_ctrl shape {b_ctrl.shape} != gene axis {(len(genes),)}"
            )
        # Some Orion tables retain multiple Ensembl features that map to the
        # same gene symbol. The evaluator works on the official unique-symbol
        # axis, so choose the duplicate column with the strongest basal signal.
        # This mirrors the extractor's deterministic duplicate resolution and
        # avoids arbitrary first/last-column behavior.
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
        if require_cd4_scale:
            meta = json.loads((path / "meta.json").read_text())
            definition = meta.get("effect_definition", "")
            if "ln(2)" not in definition or "log2 fold change" not in definition:
                raise ValueError(
                    f"{name}: unrecognized CD4 effect scale in meta.json: {definition!r}"
                )
        return cls(
            name=name,
            path=path,
            genes=genes,
            perts=perts,
            gene_pos=gene_pos,
            pert_pos={pert: i for i, pert in enumerate(perts)},
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
    """Compute column mean and mean-absolute-deviation scale in two passes."""
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
    """Return centered raw and gene-normalized effects on one fixed axis."""
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
    before = float(np.median(np.linalg.norm(raw, axis=1))) or 1.0
    after = float(np.median(np.linalg.norm(normalized, axis=1))) or 1.0
    normalized *= before / after
    enforce_memory_limit(max_rss_gib)
    return raw, normalized


def row_normalize(matrix: np.ndarray) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.float64)
    return matrix / (np.linalg.norm(matrix, axis=1, keepdims=True) + 1e-12)


def pds(prediction: np.ndarray, truth: np.ndarray) -> float:
    scores = row_normalize(prediction) @ row_normalize(truth).T
    own = np.diag(scores)[:, None]
    n = scores.shape[1]
    if n < 2:
        raise ValueError("PDS requires at least two perturbations")
    rank = (scores > own).sum(axis=1) + 0.5 * ((scores == own).sum(axis=1) - 1)
    return float(np.mean(1.0 - rank / (n - 1)))


def scale_to_reference(matrix: np.ndarray, reference_norm: float) -> np.ndarray:
    norm = float(np.median(np.linalg.norm(matrix, axis=1))) or 1.0
    return matrix * (reference_norm / norm)


def candidate_predictions(
    baseline_parts: dict[str, np.ndarray],
    cd4_parts: dict[str, np.ndarray],
    fractions: Iterable[float],
) -> dict[str, np.ndarray]:
    """Build fixed direct-donor candidates with transparent convex weights."""
    reference_norm = (
        float(
            np.median(
                [
                    np.median(np.linalg.norm(matrix, axis=1))
                    for matrix in baseline_parts.values()
                ]
            )
        )
        or 1.0
    )
    baseline_scaled = [
        scale_to_reference(matrix, reference_norm) for matrix in baseline_parts.values()
    ]
    baseline = sum(baseline_scaled) / len(baseline_scaled)
    cd4_scaled = {
        name: scale_to_reference(matrix, reference_norm)
        for name, matrix in cd4_parts.items()
    }
    condition_matrices = list(cd4_scaled.values())
    cd4_scaled["cd4_mean3"] = sum(condition_matrices) / len(condition_matrices)

    candidates = {"baseline": baseline}
    for cd4_name, cd4_matrix in cd4_scaled.items():
        candidates[f"{cd4_name}:only"] = cd4_matrix
        for fraction in fractions:
            if fraction <= 0 or fraction >= 1:
                continue
            candidates[f"{cd4_name}:fraction={fraction:g}"] = (
                1.0 - fraction
            ) * baseline + fraction * cd4_matrix
    return candidates


def choose_without_k562(rows: list[dict[str, object]]) -> dict[str, object] | None:
    """Select one candidate on the two non-K562 folds only."""
    grouped: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        if row["held"] == "k562_gwps" or row["candidate"] == "baseline":
            continue
        key = (str(row["mode"]), str(row["candidate"]))
        grouped.setdefault(key, []).append(float(row["delta_vs_baseline"]))
    eligible = [
        (float(np.mean(values)), key, min(values), len(values))
        for key, values in grouped.items()
        if len(values) == 2
    ]
    if not eligible:
        return None
    mean_gain, (mode, candidate), min_gain, n_folds = max(
        eligible, key=lambda item: (item[0], item[2], item[1])
    )
    selected_rows = [
        row for row in rows if row["mode"] == mode and row["candidate"] == candidate
    ]
    k562 = next(row for row in selected_rows if row["held"] == "k562_gwps")
    macro_delta = float(
        np.mean([float(row["delta_vs_baseline"]) for row in selected_rows])
    )
    non_k562_min = min(
        float(row["delta_vs_baseline"])
        for row in selected_rows
        if row["held"] != "k562_gwps"
    )
    return {
        "selection_rule": "maximum mean delta on held HCT116 and HEK293T only",
        "mode": mode,
        "candidate": candidate,
        "selection_mean_delta": mean_gain,
        "selection_min_delta": min_gain,
        "selection_folds": n_folds,
        "k562_pds": float(k562["pds"]),
        "k562_baseline_pds": float(k562["baseline_pds"]),
        "k562_delta": float(k562["delta_vs_baseline"]),
        "three_fold_macro_delta": macro_delta,
        "non_k562_min_delta": non_k562_min,
        "passes_primary_gate": (
            float(k562["delta_vs_baseline"]) >= 0.003
            and macro_delta > 0
            and non_k562_min >= -0.002
        ),
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

    main_paths = {
        "k562_gwps": Path(args.k562),
        "orion_hct116": Path(args.hct116),
        "orion_hek293t": Path(args.hek293t),
    }
    cd4_root = Path(args.cd4_root)
    main = {name: Bank.load(name, path) for name, path in main_paths.items()}
    cd4 = {
        name: Bank.load(name, cd4_root / name, require_cd4_scale=True)
        for name in CD4_NAMES
    }
    all_banks = {**main, **cd4}

    perturbations = [
        pert
        for pert in panel
        if all(pert in bank.pert_pos for bank in all_banks.values())
    ]
    if len(perturbations) < args.min_shared_targets:
        raise ValueError(
            f"only {len(perturbations)} panel targets shared by every bank; "
            f"require {args.min_shared_targets}"
        )

    # This fixed axis is selected before any centering or normalization.  It
    # mirrors pds_cosine's panel-wide exclusion and prevents a target's direct
    # knockdown coordinate from mixing into or dominating any scored vector.
    features = [
        gene
        for gene in official
        if gene not in panel_set
        and all(gene in bank.gene_pos for bank in all_banks.values())
    ]
    if panel_set & set(features):
        raise AssertionError("target-coordinate leakage into LODO feature axis")
    if len(features) < args.min_shared_genes:
        raise ValueError(
            f"only {len(features)} non-panel genes shared by every bank; "
            f"require {args.min_shared_genes}"
        )

    effects: dict[str, dict[str, np.ndarray]] = {}
    for name, bank in all_banks.items():
        raw, normalized = selected_effects(
            bank,
            perturbations,
            features,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        effects[name] = {"raw": raw, "gene_norm": normalized}
        print(
            f"prepared {name:16} rows={len(perturbations):3d} "
            f"genes={len(features):5d} peak_rss={peak_rss_gib():.2f} GiB",
            flush=True,
        )

    rows: list[dict[str, object]] = []
    fractions = sorted(set(args.cd4_fractions))
    for mode in ("raw", "gene_norm"):
        for held, truth_bank in main.items():
            baseline_names = [name for name in main if name != held]
            candidates = candidate_predictions(
                {name: effects[name][mode] for name in baseline_names},
                {name: effects[name][mode] for name in CD4_NAMES},
                fractions,
            )
            truth = effects[truth_bank.name]["raw"]
            baseline_score = pds(candidates["baseline"], truth)
            for candidate, prediction in candidates.items():
                score = pds(prediction, truth)
                rows.append(
                    {
                        "mode": mode,
                        "held": held,
                        "candidate": candidate,
                        "pds": score,
                        "baseline_pds": baseline_score,
                        "delta_vs_baseline": score - baseline_score,
                    }
                )

    selected = choose_without_k562(rows)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    result: dict[str, object] = {
        "protocol": {
            "perturbations": len(perturbations),
            "genes": len(features),
            "perturbation_policy": "VCC panel intersection across all six banks",
            "gene_policy": "official-axis intersection across all banks minus all VCC targets",
            "truth_space": "centered raw donor effects",
            "prediction_modes": [
                "centered raw",
                "centered mean-absolute gene normalized",
            ],
            "cd4_fractions": fractions,
            "no_whitening": True,
            "no_cell_eval2": True,
            "duplicate_gene_symbols_resolved_by_max_b_ctrl": {
                name: bank.duplicate_gene_symbols
                for name, bank in all_banks.items()
            },
        },
        "perturbations": perturbations,
        "feature_axis_sha256": hashlib.sha256("\n".join(features).encode()).hexdigest(),
        "rows": rows,
        "somatic_selected_then_k562": selected,
        "peak_rss_gib": peak_rss_gib(),
    }
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {output} and {output.with_suffix('.csv')}", flush=True)
    if selected:
        print(json.dumps(selected, indent=2), flush=True)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cd4-root",
        required=True,
        help="root produced by extract_cd4.py",
    )
    parser.add_argument("--k562", default=str(MAIN_DEFAULTS["k562_gwps"]))
    parser.add_argument("--hct116", default=str(MAIN_DEFAULTS["orion_hct116"]))
    parser.add_argument("--hek293t", default=str(MAIN_DEFAULTS["orion_hek293t"]))
    parser.add_argument(
        "--panel-csv",
        default=str(ROOT / "data/challenge/pert_counts.csv"),
    )
    parser.add_argument(
        "--genes-csv",
        default=str(ROOT / "data/challenge/gene_names.csv"),
    )
    parser.add_argument(
        "--cd4-fractions",
        nargs="+",
        type=float,
        default=[0.1, 0.2, 1 / 3, 0.5, 0.67],
    )
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=4.0)
    parser.add_argument("--min-shared-targets", type=int, default=240)
    parser.add_argument("--min-shared-genes", type=int, default=6_300)
    parser.add_argument(
        "--output",
        default=str(HERE / "results.json"),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.chunk_rows < 1:
        raise SystemExit("--chunk-rows must be positive")
    if any(not 0 < fraction < 1 for fraction in args.cd4_fractions):
        raise SystemExit("--cd4-fractions must all be strictly between 0 and 1")
    run(args)


if __name__ == "__main__":
    main()

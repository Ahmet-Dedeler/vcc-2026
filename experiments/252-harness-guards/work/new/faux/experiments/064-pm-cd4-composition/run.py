#!/usr/bin/env python3
"""Frozen 2x2 PerturbMatch and post-whitening CD4 composition test."""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import importlib.util
import json
import os
import resource
import sys
from collections.abc import Iterable, Sequence
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
BASE_PATH = ROOT / "experiments/062-cd4-post-whiten/run.py"
SPEC = importlib.util.spec_from_file_location("exp062_composition_base", BASE_PATH)
if SPEC is None or SPEC.loader is None:
    raise ImportError(f"cannot import frozen experiment-062 implementation: {BASE_PATH}")
BASE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = BASE
SPEC.loader.exec_module(BASE)

FOLDS = ("orion_hct116", "orion_hek293t")
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
PM_FRACTION = 0.50
CD4_MODE = "gene_norm"
CD4_FRACTION = 0.05
WHITEN_GAMMA = 0.65
BOOTSTRAP_REPS = 20_000
BOOTSTRAP_SEED = 6_402_026
ARMS = ("A", "B", "C", "D")
# Every whitening operator a caller of `load_frozen_operator` may legitimately load, pinned by
# NAME (experiment 252, closing 242 finding #4).  The three held operators are 062's cloud fits
# (their metadata json carries the fit-time `wcols_sha256`); the two production operators are
# the builder's `WHITENER_REGISTRY` entries.  `plolo.OPERATOR_REGISTRY` carries the same shas
# (plus metadata sha, axis digest and max eigenvalue) and asserts agreement with this table at
# import, so the two cannot drift apart silently.
KNOWN_OPERATOR_SHA256: dict[str, str] = {
    "held_k562_gwps": "13d988fcc3b129d25d254524db73e8f8253143fa493ee64d9654f0f61bc3ac55",
    "held_orion_hct116": "7cecbf7f36e40d22a6bea1ff3adc86a8788b329def397044fb5b9ba675247cbc",
    "held_orion_hek293t": "97e486e0d1a805e1f24a622f1cefbfb0917e9282c5d0491c4ddb41c8475f2af1",
    "union": "5608e807230dea04f93f5e804ed26e8f7285d9487d16a2791283793ffb80fed3",
    "union_pruned": "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123",
}
EXPECTED_WHITENER_SHA256 = {
    "orion_hct116": KNOWN_OPERATOR_SHA256["held_orion_hct116"],
    "orion_hek293t": KNOWN_OPERATOR_SHA256["held_orion_hek293t"],
}
HEX64 = frozenset("0123456789abcdef")


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def enforce_limit(limit_gib: float, stage: str) -> None:
    if limit_gib > 0 and peak_rss_gib() >= limit_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak_rss_gib():.3f} GiB >= {limit_gib:.3f} GiB"
        )


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(values).encode()).hexdigest()


def load_perturbmatch(path: Path) -> BASE.Bank:
    """Load and authenticate the frozen official PerturbMatch compact bank."""
    path = path.expanduser().resolve()
    genes = BASE.read_lines(path / "genes.txt", require_unique=False)
    perts = BASE.read_lines(path / "perts.txt")
    if len(genes) != len(set(genes)):
        raise ValueError("PerturbMatch compact axis must have unique genes")
    delta_path = path / "delta.npy"
    delta = np.load(delta_path, mmap_mode="r")
    if not isinstance(delta, np.memmap):
        raise TypeError("PerturbMatch delta.npy is not memory mapped")
    if delta.shape != (len(perts), len(genes)):
        raise ValueError(
            f"PerturbMatch shape {delta.shape} != {(len(perts), len(genes))}"
        )
    metadata = json.loads((path / "meta.json").read_text())
    if metadata.get("source_md5") != "0017076a5cc7b9979febf3bd19bd2157":
        raise ValueError("PerturbMatch source checksum is not the official artifact")
    if metadata.get("vcc_target_coverage") != 251:
        raise ValueError("PerturbMatch VCC target coverage must be 251")
    return BASE.Bank(
        name="perturbmatch",
        path=path,
        genes=genes,
        perts=perts,
        pert_pos={pert: index for index, pert in enumerate(perts)},
        delta_path=delta_path,
        b_ctrl_path=path / "b_ctrl.npy",
        shape=tuple(map(int, delta.shape)),
    )


def shared_k562_pm_columns(
    incumbent: BASE.Bank,
    perturbmatch: BASE.Bank,
    official: Sequence[str],
    panel: set[str],
) -> np.ndarray:
    incumbent_genes = set(incumbent.genes)
    perturbmatch_genes = set(perturbmatch.genes)
    return np.asarray(
        [
            index
            for index, gene in enumerate(official)
            if gene not in panel
            and gene in incumbent_genes
            and gene in perturbmatch_genes
        ],
        dtype=np.int64,
    )


def build_pm_slot(
    incumbent: BASE.PreparedSource,
    perturbmatch: BASE.PreparedSource,
    shared_columns: np.ndarray,
) -> BASE.PreparedSource:
    """Apply exp058's frozen 50:50 blend on the shared K562/PM support."""
    if not np.all(incumbent.target_available & perturbmatch.target_available):
        raise ValueError("composition gallery contains an unavailable K562/PM target")
    effects = np.asarray(incumbent.effects, dtype=np.float64).copy()
    effects[:, shared_columns] = (
        (1.0 - PM_FRACTION) * incumbent.effects[:, shared_columns]
        + PM_FRACTION * perturbmatch.effects[:, shared_columns]
    )
    return BASE.PreparedSource(
        name="k562_gwps_pm50",
        effects=effects,
        target_available=(
            incumbent.target_available & perturbmatch.target_available
        ),
        mean_on_axis=np.asarray(incumbent.mean_on_axis, dtype=np.float64).copy(),
        raw_median_norm=incumbent.raw_median_norm,
        mapped_gene_count=incumbent.mapped_gene_count,
        duplicate_symbols=incumbent.duplicate_symbols,
    )


def operator_columns(
    source_banks: Sequence[BASE.Bank],
    official: Sequence[str],
    panel: set[str],
) -> np.ndarray:
    donor_union = set().union(*(set(bank.genes) for bank in source_banks))
    return np.asarray(
        [
            index
            for index, gene in enumerate(official)
            if gene in donor_union and gene not in panel
        ],
        dtype=np.int64,
    )


def load_frozen_operator(
    path: Path,
    source_banks: Sequence[BASE.Bank],
    official: Sequence[str],
    panel: set[str],
    *,
    expected_sha256: str,
) -> tuple[BASE.Operator, np.ndarray]:
    """Load an exp062 production whitener and verify its full axis contract."""
    path = path.expanduser().resolve()
    metadata_path = path.with_suffix(".json")
    metadata = json.loads(metadata_path.read_text())
    wcols = operator_columns(source_banks, official, panel)
    actual_hash = sha256_file(path)
    # 252 (242 #4): the expectation must be a real pin.  An empty expectation used to skip the
    # check, and four callers passed `sha256_file(path)` -- the hash of the file being checked
    # -- so the comparison below was `x != x` and could never fire.  Now the digest must be
    # well-formed, must match, AND must be one this module pins by name or one the operator's
    # own metadata json declares (`sha256`, written by build_union_whitener.py at fit time).
    if not (
        isinstance(expected_sha256, str)
        and len(expected_sha256) == 64
        and set(expected_sha256) <= HEX64
    ):
        raise ValueError(
            f"{path}: expected_sha256 must be a pinned 64-hex digest, got "
            f"{expected_sha256!r}; an empty expectation is not a check"
        )
    if actual_hash != expected_sha256:
        raise ValueError(
            f"{path}: whitener SHA-256 {actual_hash} != {expected_sha256}"
        )
    pinned_name = next(
        (name for name, digest in KNOWN_OPERATOR_SHA256.items() if digest == actual_hash),
        None,
    )
    if pinned_name is None and metadata.get("sha256") != actual_hash:
        raise ValueError(
            f"{path}: sha256 {actual_hash} is neither pinned in KNOWN_OPERATOR_SHA256 "
            f"({sorted(KNOWN_OPERATOR_SHA256)}) nor declared by its own metadata json "
            f"('sha256' = {metadata.get('sha256')!r}); an expectation computed from the "
            "file under check is not a check (242 finding #4, closed by 252)"
        )
    expected_donors = [bank.name for bank in source_banks]
    if metadata.get("axis") != "union":
        raise ValueError(f"{path}: operator is not union-axis")
    if metadata.get("donors") != expected_donors:
        raise ValueError(
            f"{path}: donors {metadata.get('donors')} != {expected_donors}"
        )
    if float(metadata.get("gamma", -1)) != WHITEN_GAMMA:
        raise ValueError(f"{path}: gamma is not {WHITEN_GAMMA}")
    if int(metadata.get("n_cols", -1)) != len(wcols):
        raise ValueError(f"{path}: operator-column count mismatch")
    if metadata.get("wcols_sha256") != hashlib.sha256(wcols.tobytes()).hexdigest():
        raise ValueError(f"{path}: operator-column hash mismatch")
    with np.load(path) as archive:
        evals = np.asarray(archive["evals"], dtype=np.float64)
        evecs = np.asarray(archive["evecs"], dtype=np.float32)
    if evecs.shape != (len(wcols), len(wcols)):
        raise ValueError(f"{path}: eigenvector shape {evecs.shape} is invalid")
    return (
        BASE.Operator(evals, evecs, metadata, path, metadata_path),
        wcols,
    )


def arm_statistics(scores: dict[str, float]) -> dict[str, float]:
    return {
        "pm_main_without_cd4_B_minus_A": scores["B"] - scores["A"],
        "cd4_main_without_pm_C_minus_A": scores["C"] - scores["A"],
        "pm_main_with_cd4_D_minus_C": scores["D"] - scores["C"],
        "cd4_main_with_pm_D_minus_B": scores["D"] - scores["B"],
        "combined_D_minus_A": scores["D"] - scores["A"],
        "interaction_D_minus_B_minus_C_plus_A": (
            scores["D"] - scores["B"] - scores["C"] + scores["A"]
        ),
    }


def paired_target_cluster_bootstrap(
    per_fold_delta: dict[str, np.ndarray],
    *,
    reps: int,
    seed: int,
) -> dict[str, object]:
    arrays = [np.asarray(per_fold_delta[held], dtype=np.float64) for held in FOLDS]
    if len({len(values) for values in arrays}) != 1:
        raise ValueError("folds do not share a common target gallery")
    cluster_delta = np.mean(np.stack(arrays, axis=1), axis=1)
    rng = np.random.default_rng(seed)
    draws = np.empty(reps, dtype=np.float64)
    chunk = 1_000
    for lo in range(0, reps, chunk):
        hi = min(lo + chunk, reps)
        indices = rng.integers(
            0, len(cluster_delta), size=(hi - lo, len(cluster_delta))
        )
        draws[lo:hi] = cluster_delta[indices].mean(axis=1)
    return {
        "method": (
            "paired target-cluster bootstrap: resample target genes once per "
            "replicate, preserving each target's A/D pairing and both held folds"
        ),
        "reps": reps,
        "seed": seed,
        "estimate": float(cluster_delta.mean()),
        "ci95_percentile": [
            float(np.percentile(draws, 2.5)),
            float(np.percentile(draws, 97.5)),
        ],
        "positive_probability": float(np.mean(draws > 0)),
        "targets": len(cluster_delta),
    }


def evaluate_fold(
    held: str,
    banks: dict[str, BASE.Bank],
    prepared: dict[str, BASE.PreparedSource],
    pm_slot: BASE.PreparedSource,
    cd4: BASE.CD4Prepared,
    official: Sequence[str],
    panel: Sequence[str],
    targets: Sequence[str],
    whitener_path: Path,
    *,
    expected_whitener_sha256: str,
    max_rss_gib: float,
    min_operator_genes: int,
    min_retrieval_genes: int,
) -> tuple[dict[str, object], list[dict[str, object]], np.ndarray]:
    companion_name = next(name for name in FOLDS if name != held)
    source_banks = (banks["k562_gwps"], banks[companion_name])
    operator, wcols = load_frozen_operator(
        whitener_path,
        source_banks,
        official,
        set(panel),
        expected_sha256=expected_whitener_sha256,
    )
    if len(wcols) < min_operator_genes:
        raise ValueError(f"{held}: only {len(wcols)} operator genes")

    baseline_a_pre, votes_a = BASE.build_baseline(
        [prepared["k562_gwps"], prepared[companion_name]]
    )
    baseline_b_pre, votes_b = BASE.build_baseline(
        [pm_slot, prepared[companion_name]]
    )
    if not np.array_equal(votes_a, votes_b):
        raise AssertionError("PerturbMatch changed target source-vote counts")

    arm_a, preserve_a = BASE.apply_operator(baseline_a_pre, operator, wcols)
    arm_b, preserve_b = BASE.apply_operator(baseline_b_pre, operator, wcols)
    arm_c, cd4_preserve_c = BASE.add_post_whiten_residual(
        arm_a, cd4, CD4_FRACTION, wcols
    )
    arm_d, cd4_preserve_d = BASE.add_post_whiten_residual(
        arm_b, cd4, CD4_FRACTION, wcols
    )

    held_genes = set(banks[held].genes)
    retrieval = np.asarray(
        [index for index in wcols if official[index] in held_genes],
        dtype=np.int64,
    )
    if len(retrieval) < min_retrieval_genes:
        raise ValueError(f"{held}: only {len(retrieval)} retrieval genes")
    if set(panel) & {official[index] for index in retrieval}:
        raise AssertionError("panel coordinate entered retrieval")

    truth = BASE.raw_truth(prepared[held])
    matrices = {"A": arm_a, "B": arm_b, "C": arm_c, "D": arm_d}
    per_arm = {
        arm: BASE.pds_per_target(matrix, truth, retrieval)
        for arm, matrix in matrices.items()
    }
    scores = {arm: float(np.mean(values)) for arm, values in per_arm.items()}
    effects = arm_statistics(scores)
    per_target_rows: list[dict[str, object]] = []
    for index, target in enumerate(targets):
        values = {arm: float(per_arm[arm][index]) for arm in ARMS}
        per_target_rows.append(
            {
                "held": held,
                "target_gene": target,
                **{f"pds_{arm}": values[arm] for arm in ARMS},
                "B_minus_A": values["B"] - values["A"],
                "C_minus_A": values["C"] - values["A"],
                "D_minus_A": values["D"] - values["A"],
                "interaction": (
                    values["D"] - values["B"] - values["C"] + values["A"]
                ),
            }
        )
    fold = {
        "held": held,
        "source_donors_for_operator": [bank.name for bank in source_banks],
        "scores": scores,
        "effects": effects,
        "operator_genes": len(wcols),
        "retrieval_genes": len(retrieval),
        "operator_columns_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
        "retrieval_columns_sha256": hashlib.sha256(retrieval.tobytes()).hexdigest(),
        "whitener_path": str(whitener_path.resolve()),
        "whitener_sha256": sha256_file(whitener_path),
        "whitener_metadata_sha256": sha256_file(
            whitener_path.with_suffix(".json")
        ),
        "production_preserve": {
            "A": preserve_a,
            "B": preserve_b,
            "C_post_cd4": cd4_preserve_c,
            "D_post_cd4": cd4_preserve_d,
        },
        "source_vote_counts": {
            str(value): int(np.sum(votes_a == value))
            for value in sorted(set(votes_a.tolist()))
        },
    }
    d_minus_a = per_arm["D"] - per_arm["A"]
    del (
        operator,
        baseline_a_pre,
        baseline_b_pre,
        arm_a,
        arm_b,
        arm_c,
        arm_d,
        truth,
        matrices,
        per_arm,
    )
    gc.collect()
    enforce_limit(max_rss_gib, f"{held} fold completion")
    return fold, per_target_rows, d_minus_a


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    if not rows:
        return
    flattened = [
        {
            key: json.dumps(value, sort_keys=True)
            if isinstance(value, (dict, list))
            else value
            for key, value in row.items()
        }
        for row in rows
    ]
    fields = sorted({key for row in flattened for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(flattened)


def bank_hashes(bank: BASE.Bank, *, include_b_ctrl: bool = True) -> dict[str, str]:
    names = ["delta.npy", "genes.txt", "perts.txt", "meta.json"]
    if include_b_ctrl:
        names.insert(1, "b_ctrl.npy")
    return {
        name: sha256_file(bank.path / name)
        for name in names
        if (bank.path / name).exists()
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    panel = BASE.read_csv_column(args.panel_csv, "target_gene")
    official = BASE.read_csv_column(args.genes_csv, "gene_name")
    if len(panel) != args.expected_panel_targets:
        raise ValueError(f"expected {args.expected_panel_targets} panel targets")
    official_pos = {gene: index for index, gene in enumerate(official)}
    if any(gene not in official_pos for gene in panel):
        raise ValueError("a panel target is absent from the official gene axis")
    panel_positions = np.asarray([official_pos[gene] for gene in panel], dtype=np.int64)
    panel_set = set(panel)

    banks = {
        "k562_gwps": BASE.Bank.load("k562_gwps", args.k562),
        "orion_hct116": BASE.Bank.load("orion_hct116", args.hct116),
        "orion_hek293t": BASE.Bank.load("orion_hek293t", args.hek293t),
    }
    perturbmatch = load_perturbmatch(args.perturbmatch)
    cd4_banks = tuple(
        BASE.Bank.load(name, args.cd4_root / name, require_cd4_scale=True)
        for name in CD4_NAMES
    )
    all_banks = [*banks.values(), perturbmatch, *cd4_banks]
    targets = tuple(
        target for target in panel if all(target in bank.pert_pos for bank in all_banks)
    )
    if len(targets) != args.expected_shared_targets:
        raise ValueError(
            f"shared gallery {len(targets)} != expected {args.expected_shared_targets}"
        )

    prepared = {
        name: BASE.prepare_main_source(
            bank,
            targets,
            official,
            panel_positions,
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        for name, bank in banks.items()
    }
    prepared_pm = BASE.prepare_main_source(
        perturbmatch,
        targets,
        official,
        panel_positions,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    shared_columns = shared_k562_pm_columns(
        banks["k562_gwps"], perturbmatch, official, panel_set
    )
    if len(shared_columns) < args.min_pm_shared_genes:
        raise ValueError(f"only {len(shared_columns)} shared K562/PM genes")
    pm_slot = build_pm_slot(prepared["k562_gwps"], prepared_pm, shared_columns)
    cd4 = BASE.prepare_cd4_mean(
        cd4_banks,
        targets,
        official,
        panel_set,
        CD4_MODE,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    if not np.all(cd4.target_condition_count == len(CD4_NAMES)):
        raise AssertionError("shared gallery lacks a CD4 condition")

    fold_rows: list[dict[str, object]] = []
    per_target_rows: list[dict[str, object]] = []
    d_minus_a: dict[str, np.ndarray] = {}
    for held in FOLDS:
        expected_hash = (
            args.expected_hct_whitener_sha256
            if held == "orion_hct116"
            else args.expected_hek_whitener_sha256
        )
        fold, rows, paired = evaluate_fold(
            held,
            banks,
            prepared,
            pm_slot,
            cd4,
            official,
            panel,
            targets,
            args.whitener_dir / f"held_{held}.npz",
            expected_whitener_sha256=expected_hash,
            max_rss_gib=args.max_rss_gib,
            min_operator_genes=args.min_operator_genes,
            min_retrieval_genes=args.min_retrieval_genes,
        )
        fold_rows.append(fold)
        per_target_rows.extend(rows)
        d_minus_a[held] = paired
        print(
            f"held={held} A={fold['scores']['A']:.9f} "
            f"D={fold['scores']['D']:.9f} "
            f"D-A={fold['effects']['combined_D_minus_A']:+.9f} "
            f"peak_rss={peak_rss_gib():.3f} GiB",
            flush=True,
        )

    macro_scores = {
        arm: float(np.mean([fold["scores"][arm] for fold in fold_rows]))
        for arm in ARMS
    }
    macro_effects = arm_statistics(macro_scores)
    bootstrap = paired_target_cluster_bootstrap(
        d_minus_a, reps=args.bootstrap_reps, seed=args.bootstrap_seed
    )
    if not np.isclose(
        bootstrap["estimate"], macro_effects["combined_D_minus_A"], atol=1e-14
    ):
        raise AssertionError("bootstrap point estimate does not match macro D-A")

    output_dir = args.output.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    result: dict[str, object] = {
        "experiment": "064-pm-cd4-composition",
        "status": "completed",
        "protocol": {
            "design": "frozen 2x2 factorial; no tuning",
            "arms": {
                "A": "incumbent production union-whitened baseline",
                "B": "A with fixed 50:50 incumbent/PerturbMatch K562 slot",
                "C": "A plus fixed post-whitening CD4 mean3 gene_norm f0.05",
                "D": "both B and C changes",
            },
            "valid_held_roles": list(FOLDS),
            "invalid_holdout_omitted": (
                "K562: PerturbMatch and incumbent K562 use the same assay"
            ),
            "gallery": "complete target intersection across all seven banks",
            "whitening_gamma": WHITEN_GAMMA,
            "whitener": (
                "frozen candidate-independent exp062 source-only production "
                "union operator from incumbent K562 and non-held Orion donor"
            ),
            "production_preserve_renorm": True,
            "perturbmatch_fraction": PM_FRACTION,
            "perturbmatch_support": (
                "blend on incumbent-K562/PerturbMatch shared official non-panel "
                "genes; retain incumbent K562 outside shared support"
            ),
            "cd4_mode": CD4_MODE,
            "cd4_fraction": CD4_FRACTION,
            "truth": "raw centered, unnormalized, untransformed",
            "count_emission": False,
            "cell_eval2": False,
            "submission": False,
        },
        "axis": {
            "official_genes": len(official),
            "challenge_targets": len(panel),
            "gallery_targets": len(targets),
            "target_order": list(targets),
            "target_order_sha256": sha256_lines(targets),
            "panel_order_sha256": sha256_lines(panel),
            "k562_perturbmatch_shared_non_panel_genes": len(shared_columns),
            "k562_perturbmatch_shared_columns_sha256": hashlib.sha256(
                shared_columns.tobytes()
            ).hexdigest(),
        },
        "folds": fold_rows,
        "macro": {
            "scores": macro_scores,
            "effects": macro_effects,
            "D_vs_A_bootstrap": bootstrap,
        },
        "leakage_audit": {
            "all_300_panel_coordinates_zeroed_in_main_predictions": True,
            "all_300_panel_coordinates_excluded_from_cd4": True,
            "all_300_panel_coordinates_excluded_from_operators": True,
            "all_300_panel_coordinates_excluded_from_retrieval": True,
            "operator_candidate_independent": True,
            "perturbmatch_used_in_operator": False,
            "truth_transformed": False,
            "hyperparameters_retuned": False,
            "k562_truth_evaluated": False,
        },
        "input_hashes": {
            "run_py": sha256_file(Path(__file__)),
            "base_exp062_run_py": sha256_file(BASE_PATH),
            "panel_csv": sha256_file(args.panel_csv),
            "genes_csv": sha256_file(args.genes_csv),
            "banks": {
                **{name: bank_hashes(bank) for name, bank in banks.items()},
                "perturbmatch": bank_hashes(
                    perturbmatch, include_b_ctrl=False
                ),
                **{bank.name: bank_hashes(bank) for bank in cd4_banks},
            },
        },
        "infrastructure": {
            "provider": "AWS EC2",
            "instance_id": args.instance_id,
            "instance_type": args.instance_type,
            "market": args.market,
            "s3_prefix": args.s3_prefix,
        },
        "peak_rss_gib": peak_rss_gib(),
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    write_csv(
        output_dir / "fold_summary.csv",
        [
            {
                "held": fold["held"],
                **{f"pds_{arm}": fold["scores"][arm] for arm in ARMS},
                **fold["effects"],
            }
            for fold in fold_rows
        ],
    )
    write_csv(output_dir / "per_target_pds.csv", per_target_rows)
    print(json.dumps(result["macro"], indent=2), flush=True)
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--k562", type=Path, required=True)
    parser.add_argument("--hct116", type=Path, required=True)
    parser.add_argument("--hek293t", type=Path, required=True)
    parser.add_argument("--perturbmatch", type=Path, required=True)
    parser.add_argument("--cd4-root", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--whitener-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=55.0)
    parser.add_argument("--expected-panel-targets", type=int, default=300)
    parser.add_argument("--expected-shared-targets", type=int, default=221)
    parser.add_argument("--min-pm-shared-genes", type=int, default=6_900)
    parser.add_argument("--min-operator-genes", type=int, default=17_000)
    parser.add_argument("--min-retrieval-genes", type=int, default=17_000)
    parser.add_argument("--bootstrap-reps", type=int, default=BOOTSTRAP_REPS)
    parser.add_argument("--bootstrap-seed", type=int, default=BOOTSTRAP_SEED)
    parser.add_argument(
        "--expected-hct-whitener-sha256",
        default=EXPECTED_WHITENER_SHA256["orion_hct116"],
    )
    parser.add_argument(
        "--expected-hek-whitener-sha256",
        default=EXPECTED_WHITENER_SHA256["orion_hek293t"],
    )
    parser.add_argument("--instance-id", default="unknown")
    parser.add_argument("--instance-type", default="unknown")
    parser.add_argument("--market", default="unknown")
    parser.add_argument("--s3-prefix", default="unknown")
    parser.add_argument("--guarded", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.guarded:
        raise SystemExit("refusing unguarded real-data execution")
    if args.chunk_rows < 1 or args.max_rss_gib <= 0:
        raise SystemExit("chunk rows and memory limit must be positive")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

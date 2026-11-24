#!/usr/bin/env python3
"""Build the frozen PerturbMatch + post-whitening CD4 HCT116 candidate.

The pre-whitening candidate is delegated to experiment 061 verbatim, including
its missing-target and missing-gene PerturbMatch behavior. The resulting exact
experiment-061 effect matrix is checked against that run's build report. The
experiment-062 CD4 mean3 gene-normalized residual is then added after whitening
at fraction 0.05, followed by the same one-scalar global norm preservation.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path
from types import ModuleType

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import anndata as ad
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
EXP061_PATH = ROOT / "experiments/061-perturbmatch-hct-count/build_candidate.py"
EXP062_PATH = ROOT / "experiments/062-cd4-post-whiten/run.py"


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


EXP061 = load_module("exp065_frozen_exp061", EXP061_PATH)
EXP062 = load_module("exp065_frozen_exp062", EXP062_PATH)

CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
CD4_MODE = "gene_norm"
CD4_FRACTION = 0.05
EXPECTED_GAMMA = 0.65
EXPECTED_ALPHA = 1.08
EXPECTED_MIN_CPM = 5.0
EXPECTED_SPREAD = 0.25
EXPECTED_PM_BLEND = 0.50
EXPECTED_SEED = 0


def assert_frozen_contract() -> None:
    observed = {
        "gamma": float(EXP061.GAMMA),
        "alpha": float(EXP061.ALPHA),
        "cpm_gate": float(EXP061.MIN_CPM),
        "spread": float(EXP061.SPREAD),
        "pm_blend": float(EXP061.BLEND),
        "seed": int(EXP061.SEED),
        "cd4_effect_scale": float(EXP062.EFFECT_SCALE),
        "cd4_gamma": float(EXP062.WHITEN_GAMMA),
    }
    expected = {
        "gamma": EXPECTED_GAMMA,
        "alpha": EXPECTED_ALPHA,
        "cpm_gate": EXPECTED_MIN_CPM,
        "spread": EXPECTED_SPREAD,
        "pm_blend": EXPECTED_PM_BLEND,
        "seed": EXPECTED_SEED,
        "cd4_effect_scale": 6.817,
        "cd4_gamma": EXPECTED_GAMMA,
    }
    if observed != expected:
        raise ValueError(
            "frozen source contract drifted:\n"
            + json.dumps({"expected": expected, "observed": observed}, indent=2)
        )


def file_sha256(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def load_exp061_build_report(path: Path) -> dict[str, object]:
    report = json.loads(path.read_text())
    if report.get("experiment") != "061-perturbmatch-hct-count":
        raise ValueError(f"{path}: not an experiment-061 build report")
    if report.get("status") != "candidate_built":
        raise ValueError(f"{path}: experiment-061 candidate was not built")
    protocol = report.get("protocol", {})
    expected = {
        "gamma": EXPECTED_GAMMA,
        "alpha": EXPECTED_ALPHA,
        "cpm_gate": EXPECTED_MIN_CPM,
        "spread": EXPECTED_SPREAD,
        "seed": EXPECTED_SEED,
    }
    observed = {key: protocol.get(key) for key in expected}
    if observed != expected:
        raise ValueError(
            f"{path}: experiment-061 frozen parameters differ: "
            f"{observed} != {expected}"
        )
    if protocol.get("candidate_independent_operator") is not True:
        raise ValueError(f"{path}: experiment-061 operator is not candidate-independent")
    if protocol.get("all_300_vcc_target_coordinates_excluded") is not True:
        raise ValueError(f"{path}: experiment-061 target-coordinate audit failed")
    return report


def freeze_complete_mean3(prepared):
    """Keep CD4 residuals only for targets observed in all three conditions."""
    complete = prepared.target_condition_count == len(CD4_NAMES)
    mean_effect = np.asarray(prepared.mean_effect, dtype=np.float64).copy()
    gene_votes = np.asarray(prepared.gene_votes, dtype=np.uint8).copy()
    mean_effect[~complete] = 0.0
    gene_votes[~complete] = 0
    return EXP062.CD4Prepared(
        mode=prepared.mode,
        mean_effect=mean_effect,
        gene_votes=gene_votes,
        target_condition_count=prepared.target_condition_count,
        condition_details=prepared.condition_details,
    )


def prepare_cd4(
    cd4_root: Path,
    targets: tuple[str, ...],
    official: tuple[str, ...],
    panel: tuple[str, ...],
    *,
    chunk_rows: int,
    max_rss_gib: float,
):
    banks = tuple(
        EXP062.Bank.load(name, cd4_root / name, require_cd4_scale=True)
        for name in CD4_NAMES
    )
    prepared = EXP062.prepare_cd4_mean(
        banks,
        targets,
        official,
        set(panel),
        CD4_MODE,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    # Experiment 064 defined "mean3" only on targets present in every CD4
    # condition. Extend that frozen candidate conservatively: incomplete rows
    # get no CD4 residual instead of averaging a newly chosen subset.
    return banks, freeze_complete_mean3(prepared)


def add_frozen_cd4_residual(
    pm_final: np.ndarray,
    cd4,
    wcols: np.ndarray,
    panel_positions: np.ndarray,
) -> tuple[np.ndarray, dict[str, float]]:
    """Apply exp062's exact post-whiten addition and verify excluded support."""
    panel_positions = np.asarray(panel_positions, dtype=np.int64)
    if np.any(cd4.gene_votes[:, panel_positions] != 0):
        raise AssertionError("a VCC target coordinate entered CD4 gene support")
    if np.any(cd4.mean_effect[:, panel_positions] != 0):
        raise AssertionError("a VCC target coordinate entered the CD4 residual")

    joint, preserve = EXP062.add_post_whiten_residual(
        pm_final, cd4, CD4_FRACTION, wcols
    )
    scale = float(preserve["post_add_global_preserve_scale"])
    np.testing.assert_allclose(
        joint[:, panel_positions],
        np.asarray(pm_final)[:, panel_positions] * scale,
        rtol=0.0,
        atol=0.0,
        err_msg="panel coordinates changed by more than the global preserve scalar",
    )
    unavailable = cd4.gene_votes == 0
    np.testing.assert_allclose(
        joint[unavailable],
        np.asarray(pm_final)[unavailable] * scale,
        rtol=0.0,
        atol=0.0,
        err_msg="unsupported CD4 coordinates changed by more than global scaling",
    )
    return joint, preserve


def bank_hashes(bank) -> dict[str, str]:
    return {
        name: file_sha256(bank.path / name)
        for name in ("delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt", "meta.json")
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    assert_frozen_contract()
    panel = EXP061.read_csv_column(args.panel_csv, "target_gene")
    official = EXP061.read_csv_column(args.genes_csv, "gene_name")
    if len(panel) != 300:
        raise ValueError(f"expected 300 VCC targets, got {len(panel)}")
    official_pos = {gene: index for index, gene in enumerate(official)}
    if any(gene not in official_pos for gene in panel):
        raise ValueError("a VCC panel target is absent from the official gene axis")
    panel_positions = np.asarray(
        [official_pos[gene] for gene in panel], dtype=np.int64
    )

    real_backed = ad.read_h5ad(args.real, backed="r")
    try:
        genes = tuple(map(str, real_backed.var_names))
        targets = tuple(
            sorted(
                set(real_backed.obs[EXP061.PERT_COL].astype(str))
                - {EXP061.CONTROL}
            )
        )
    finally:
        real_backed.file.close()
    if genes != official:
        raise ValueError("HCT116 gene axis differs from the official VCC axis")
    if not set(targets) <= set(panel):
        raise ValueError("HCT116 targets are not a subset of the VCC panel")

    exp061_report = load_exp061_build_report(args.exp061_build_report)
    input_paths = {
        "real": args.real,
        "k562_delta": args.k562 / "delta.npy",
        "hek_delta": args.hek293t / "delta.npy",
        "perturbmatch_delta": args.perturbmatch / "delta.npy",
        "whitener": args.whitener,
        "whitener_metadata": args.whitener_metadata,
        "panel_csv": args.panel_csv,
        "genes_csv": args.genes_csv,
        "exp061_build_report": args.exp061_build_report,
        "exp061_source": EXP061_PATH,
        "exp062_source": EXP062_PATH,
    }
    input_hashes = {name: file_sha256(path) for name, path in input_paths.items()}
    EXP061.memory_checkpoint("input hashes", args.max_rss_gib)

    donors = {
        EXP061.K562: EXP061.Donor.load(args.k562),
        EXP061.HEK: EXP061.Donor.load(args.hek293t),
    }
    model = EXP061.MultiDonorModel(
        donors, np.asarray(genes), gene_norm=False, eb=False
    )
    pm = EXP061.prepare_perturbmatch(
        args.perturbmatch,
        genes,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    baseline_pre, pm_pre, blend_audit = EXP061.compute_effects(
        model, pm, targets, max_rss_gib=args.max_rss_gib
    )
    whitener, wcols, operator_metadata = EXP061.validate_and_load_operator(
        model, panel, args.whitener, args.whitener_metadata
    )
    if np.intersect1d(wcols, panel_positions).size:
        raise AssertionError("a VCC target coordinate entered the fixed whitener")
    baseline_final, pm_final, scaling = EXP061.apply_frozen_pipeline(
        baseline_pre, pm_pre, whitener, wcols
    )

    expected_effect_audit = exp061_report["effect_audit"]
    reconstructed = {
        "baseline_pre_sha256": EXP061.array_sha256(baseline_pre),
        "candidate_pre_sha256": EXP061.array_sha256(pm_pre),
        "baseline_final_sha256": EXP061.array_sha256(baseline_final),
        "candidate_final_sha256": EXP061.array_sha256(pm_final),
    }
    for key, value in reconstructed.items():
        if expected_effect_audit.get(key) != value:
            raise ValueError(
                f"experiment-061 reconstruction mismatch for {key}: "
                f"{value} != {expected_effect_audit.get(key)}"
            )

    missing_pm = set(blend_audit["missing_perturbmatch_targets"])
    missing_rows = [
        index for index, target in enumerate(targets) if target in missing_pm
    ]
    if missing_rows and not np.array_equal(
        pm_final[missing_rows], baseline_final[missing_rows]
    ):
        raise AssertionError("PerturbMatch-missing rows changed before CD4 addition")

    cd4_banks, cd4 = prepare_cd4(
        args.cd4_root,
        targets,
        official,
        panel,
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    joint_final, cd4_preserve = add_frozen_cd4_residual(
        pm_final, cd4, wcols, panel_positions
    )
    effect_audit = {
        **reconstructed,
        "joint_final_sha256": EXP061.array_sha256(joint_final),
        "missing_pm_rows_equal_to_baseline_before_cd4": True,
        "missing_pm_targets": sorted(missing_pm),
        "changed_pm_pre_entries": int(np.count_nonzero(pm_pre != baseline_pre)),
        "changed_pm_final_entries": int(
            np.count_nonzero(pm_final != baseline_final)
        ),
        "changed_joint_vs_pm_entries": int(
            np.count_nonzero(joint_final != pm_final)
        ),
        "panel_cd4_residual_nonzero_entries": int(
            np.count_nonzero(cd4.mean_effect[:, panel_positions])
        ),
    }

    del (
        donors,
        model,
        pm,
        whitener,
        baseline_pre,
        pm_pre,
        baseline_final,
        pm_final,
    )
    gc.collect()
    EXP061.memory_checkpoint("before count emission", args.max_rss_gib)
    emission = EXP061.emit_prediction(
        args.real,
        targets,
        joint_final,
        args.output,
        max_rss_gib=args.max_rss_gib,
    )
    output_hash = file_sha256(args.output)
    report = {
        "experiment": "065-pm-cd4-hct-count",
        "status": "candidate_built",
        "protocol": {
            "held_out_truth": "Orion HCT116",
            "base_candidate": (
                "exact experiment-061 50:50 incumbent/PerturbMatch K562-slot "
                "candidate, verified by all four pre/post effect hashes"
            ),
            "candidate_change": (
                "add experiment-062 CD4 mean3 gene_norm residual after whitening "
                "at fraction 0.05, then apply its global preserve renormalizer"
            ),
            "missing_cd4_target_support": (
                "apply no CD4 residual unless the target exists in all three "
                "conditions, matching experiment 064's frozen mean3 gallery"
            ),
            "candidate_independent_operator": True,
            "operator_fit_sources": [EXP061.K562, EXP061.HEK],
            "all_300_vcc_target_coordinates_excluded_from_operator_and_cd4": True,
            "missing_perturbmatch_support_behavior": (
                "delegated verbatim to experiment 061; missing targets and genes "
                "retain incumbent K562 values before whitening"
            ),
            "gamma": EXPECTED_GAMMA,
            "alpha": EXPECTED_ALPHA,
            "cpm_gate": EXPECTED_MIN_CPM,
            "emission": "quasi exact-total",
            "spread": EXPECTED_SPREAD,
            "seed": EXPECTED_SEED,
            "cd4_conditions": list(CD4_NAMES),
            "cd4_mode": CD4_MODE,
            "cd4_fraction": CD4_FRACTION,
            "submission": False,
        },
        "axis": {
            "hct_targets": len(targets),
            "official_genes": len(official),
            "operator_genes": len(wcols),
            "operator_wcols_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
            "panel_positions_sha256": hashlib.sha256(
                panel_positions.tobytes()
            ).hexdigest(),
        },
        "blend_audit": blend_audit,
        "cd4_audit": {
            "mode": cd4.mode,
            "conditions": cd4.condition_details,
            "target_condition_count": {
                str(value): int(np.sum(cd4.target_condition_count == value))
                for value in sorted(set(cd4.target_condition_count.tolist()))
            },
            "supported_entries": int(np.count_nonzero(cd4.gene_votes)),
            "fully_supported_targets": int(
                np.sum(cd4.target_condition_count == len(CD4_NAMES))
            ),
            "incomplete_targets_preserved_before_global_scale": int(
                np.sum(cd4.target_condition_count != len(CD4_NAMES))
            ),
            "post_add_preserve": cd4_preserve,
        },
        "operator_metadata": operator_metadata,
        "exp061_scaling": scaling,
        "effect_audit": effect_audit,
        "emission": emission,
        "inputs_sha256": input_hashes,
        "cd4_inputs_sha256": {
            bank.name: bank_hashes(bank) for bank in cd4_banks
        },
        "candidate_h5ad": {
            "path": str(args.output),
            "bytes": args.output.stat().st_size,
            "sha256": output_hash,
        },
        "peak_rss_gib": EXP061.peak_rss_gib(),
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
    parser.add_argument("--cd4-root", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--genes-csv", type=Path, required=True)
    parser.add_argument("--whitener", type=Path, required=True)
    parser.add_argument("--whitener-metadata", type=Path, required=True)
    parser.add_argument("--exp061-build-report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=28.0)
    return parser


if __name__ == "__main__":
    run(build_parser().parse_args())

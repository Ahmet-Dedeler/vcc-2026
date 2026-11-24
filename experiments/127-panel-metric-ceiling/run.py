#!/usr/bin/env python3
"""Read-only panel-coordinate ceilings on already-consumed development surfaces."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
FREEZE_PATH = HERE / "IMPLEMENTATION_FREEZE.json"
TS = 50_000.0
ALPHA = 0.05
PURITY = 0.9


class AnalysisError(RuntimeError):
    pass


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AnalysisError(message)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 << 20):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text())
    require(isinstance(value, dict), f"{path}: JSON root must be an object")
    return value


def read_panel(path: Path) -> tuple[str, ...]:
    with path.open(newline="") as handle:
        values = tuple(row["target_gene"] for row in csv.DictReader(handle))
    require(len(values) == 300 and len(set(values)) == 300, "panel axis changed")
    return values


def verify_freeze(protocol: dict[str, Any]) -> None:
    freeze = load_json(FREEZE_PATH)
    require(
        freeze["status"] == "frozen_before_biological_value_array_access",
        "invalid freeze status",
    )
    require(freeze["protocol"]["sha256"] == sha256_file(PROTOCOL_PATH), "protocol changed")
    for relative, binding in freeze["inputs_verified_as_opaque_bytes"].items():
        path = ROOT / relative
        require(path.stat().st_size == binding["bytes"], f"input size changed: {relative}")
        require(sha256_file(path) == binding["sha256"], f"input hash changed: {relative}")
    require(protocol["hard_exclusions"]["cell_eval2"], "cell-eval2 exclusion changed")


@dataclass
class Reference:
    genes: np.ndarray
    gate: np.ndarray
    gate_genes: np.ndarray
    ref_lin: np.ndarray
    targets: np.ndarray
    lfc: np.ndarray
    sig: list[np.ndarray]


@dataclass
class Prediction:
    genes: np.ndarray
    targets: np.ndarray
    bctrl: np.ndarray
    ctrl_cpm: np.ndarray
    effect: np.ndarray


def load_reference(path: Path) -> Reference:
    with np.load(path, allow_pickle=False) as z:
        genes = z["genes"].astype(str)
        gate = z["gate"].astype(bool)
        gate_genes = z["gk"].astype(str)
        ref_lin = z["ref_lin"].astype(np.float64)
        targets = z["perts"].astype(str)
        lfc = z["lfc"].astype(np.float64)
        ptr = z["sig_ptr"].astype(np.int64)
        idx = z["sig_idx"].astype(np.int64)
    require(np.array_equal(genes[gate], gate_genes), "reference gated axis mismatch")
    require(lfc.shape == (len(targets), len(gate_genes)), "reference LFC shape mismatch")
    sig = [idx[ptr[i] : ptr[i + 1]] for i in range(len(targets))]
    return Reference(genes, gate, gate_genes, ref_lin, targets, lfc, sig)


def load_prediction(path: Path, effect_key: str) -> Prediction:
    with np.load(path, allow_pickle=False) as z:
        return Prediction(
            z["genes"].astype(str),
            z["targets"].astype(str),
            z["bctrl"].astype(np.float64),
            z["ctrl_cpm"].astype(np.float64),
            z[effect_key].astype(np.float64),
        )


def profile_from_effect(bctrl: np.ndarray, effect: np.ndarray) -> np.ndarray:
    profile = np.expm1(bctrl + effect) / TS
    np.clip(profile, 0.0, None, out=profile)
    total = float(profile.sum())
    require(total > 0 and np.isfinite(total), "invalid incumbent profile")
    return profile / total


def lfc_from_profile(profile: np.ndarray, gate: np.ndarray, ref_lin: np.ndarray) -> np.ndarray:
    return np.log2((1e6 * profile[gate] + 1e-9) / (ref_lin + 1e-9))


def nmae(pred: np.ndarray, truth: np.ndarray, sig: np.ndarray) -> float | None:
    finite = np.isfinite(truth[sig])
    selected = sig[finite]
    if selected.size < 10:
        return None
    denominator = float(np.abs(truth[selected]).sum())
    if not (denominator > 0 and np.isfinite(denominator)):
        return None
    values = np.nan_to_num(pred[selected], nan=0.0, posinf=0.0, neginf=0.0)
    return float(np.abs(values - truth[selected]).sum() / denominator)


def deepest_reach(matches: np.ndarray) -> float:
    if matches.size == 0:
        return 0.0
    purity = np.cumsum(matches, dtype=np.float64) / np.arange(1, matches.size + 1)
    good = np.flatnonzero(purity >= PURITY)
    return float((good[-1] + 1) / matches.size) if good.size else 0.0


def prediction_order(
    pred_sig: np.ndarray,
    pred_padj: np.ndarray,
    pred_p: np.ndarray,
    pred_lfc: np.ndarray,
    gene_names: np.ndarray,
) -> np.ndarray:
    return np.lexsort(
        (
            gene_names,
            -np.abs(np.nan_to_num(pred_lfc, nan=0.0)),
            np.nan_to_num(pred_p, nan=np.inf),
            np.nan_to_num(pred_padj, nan=np.inf),
            ~pred_sig,
        )
    )


def fidelity_oracle(
    fixed_n_pred: int,
    fixed_correct: int,
    n_real: int,
    mutable_adjudicable: int,
) -> tuple[float, int]:
    best = -1.0
    best_calls = 0
    for calls in range(mutable_adjudicable + 1):
        value = (fixed_correct + calls) / max(fixed_n_pred + calls, n_real, 1)
        if value > best:
            best = value
            best_calls = calls
    return best, best_calls


def fit_h1_nmae_scale(protocol: dict[str, Any]) -> dict[str, float]:
    records = []
    for item in protocol["surfaces"]["h1"]["scaled_nmae_calibration"]["sources"]:
        value = json.loads((ROOT / item["path"]).read_text())
        records.extend(value if isinstance(value, list) else [value])
    x = np.asarray([row["raw"]["de_wilcoxon_lfc_nmae"] for row in records], dtype=float)
    y = np.asarray([row["official"]["de_wilcoxon_lfc_nmae"] for row in records], dtype=float)
    slope, intercept = np.polyfit(x, y, 1)
    residual = float(np.max(np.abs(y - (slope * x + intercept))))
    require(slope < 0 and residual < 1e-10, "H1 NMAE affine scale was not exact")
    return {
        "slope": float(slope),
        "intercept": float(intercept),
        "max_absolute_fit_residual": residual,
        "points": int(len(x)),
        "raw_span_absolute": float(1.0 / abs(slope)),
    }


def load_exact_de(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as z:
        return {key: z[key] for key in z.files}


def surface_analysis(
    name: str,
    config: dict[str, Any],
    panel: tuple[str, ...],
    h1_scale: dict[str, float] | None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    ref = load_reference(ROOT / config["reference_de"]["path"])
    pred = load_prediction(
        ROOT / config["incumbent_prediction"]["path"],
        config["incumbent_prediction"]["effect_key"],
    )
    require(np.array_equal(ref.genes, pred.genes), f"{name}: gene axis mismatch")
    require(np.array_equal(ref.targets, pred.targets), f"{name}: target axis mismatch")
    panel_position = {gene: i for i, gene in enumerate(panel)}
    gene_position = {gene: i for i, gene in enumerate(ref.genes)}
    gate_position = {gene: i for i, gene in enumerate(ref.gate_genes)}
    require(all(gene in gene_position for gene in panel), f"{name}: panel gene missing")
    target_rows = [i for i, target in enumerate(ref.targets) if target in panel_position]
    exact = (
        load_exact_de(ROOT / config["exact_prediction_de_table"]["path"])
        if config["exact_prediction_de_table"]
        else None
    )
    if exact is not None:
        require(np.array_equal(exact["genes"].astype(str), ref.gate_genes), f"{name}: exact DE gene axis mismatch")
        require(np.array_equal(exact["perts"].astype(str), ref.targets), f"{name}: exact DE target axis mismatch")

    rows: list[dict[str, Any]] = []
    for row in target_rows:
        target = str(ref.targets[row])
        target_gate = gate_position.get(target)
        sig = np.asarray(ref.sig[row], dtype=np.int64)
        sig = sig[np.isfinite(ref.lfc[row, sig]) & (ref.lfc[row, sig] != 0)]
        if target_gate is not None:
            sig = sig[sig != target_gate]
        panel_gate_mask = np.asarray(
            [gene in panel_position and gene != target for gene in ref.gate_genes],
            dtype=bool,
        )
        panel_sig = sig[panel_gate_mask[sig]]
        denominator = float(np.abs(ref.lfc[row, sig]).sum()) if sig.size else 0.0
        panel_weight = float(np.abs(ref.lfc[row, panel_sig]).sum()) if panel_sig.size else 0.0

        baseline_profile = profile_from_effect(pred.bctrl, pred.effect[row])
        baseline_lfc = lfc_from_profile(baseline_profile, ref.gate, ref.ref_lin)
        incumbent_nmae = nmae(baseline_lfc, ref.lfc[row], sig)
        oracle_lfc = baseline_lfc.copy()
        oracle_lfc[panel_gate_mask] = ref.lfc[row, panel_gate_mask]
        oracle_nmae = nmae(oracle_lfc, ref.lfc[row], sig)

        record: dict[str, Any] = {
            "surface": name,
            "target": target,
            "n_reference_de": int(sig.size),
            "n_reference_de_panel": int(panel_sig.size),
            "reference_de_panel_fraction": float(panel_sig.size / sig.size) if sig.size else None,
            "reference_abs_lfc_panel_fraction": float(panel_weight / denominator) if denominator > 0 else None,
            "zero_to_perfect_panel_raw_nmae_ceiling": (
                float(panel_weight / denominator) if denominator > 0 else None
            ),
            "incumbent_nmae": incumbent_nmae,
            "oracle_perfect_panel_nmae": oracle_nmae,
            "oracle_perfect_panel_raw_nmae_gain": (
                float(incumbent_nmae - oracle_nmae)
                if incumbent_nmae is not None and oracle_nmae is not None
                else None
            ),
            "oracle_scaled_nmae_gain": (
                float((incumbent_nmae - oracle_nmae) * abs(h1_scale["slope"]))
                if h1_scale and incumbent_nmae is not None and oracle_nmae is not None
                else None
            ),
            "exact_prediction_de_available": exact is not None,
            "oracle_fidelity": None,
            "oracle_fidelity_raw_gain": None,
            "incumbent_fidelity": None,
            "oracle_jaccard": None,
            "oracle_jaccard_raw_gain": None,
            "incumbent_jaccard": None,
            "oracle_reach": None,
            "oracle_reach_raw_gain": None,
            "incumbent_reach": None,
        }

        if exact is not None:
            real_lfc = exact["real_lfc"][row].astype(np.float64)
            real_sig = (exact["real_padj"][row] < ALPHA) & np.isfinite(real_lfc) & (real_lfc != 0)
            pred_lfc = exact["pred_lfc"][row].astype(np.float64)
            pred_p = exact["pred_p"][row].astype(np.float64)
            pred_padj = exact["pred_padj"][row].astype(np.float64)
            pred_sig = pred_padj < ALPHA
            if target_gate is not None:
                real_sig[target_gate] = False
                pred_sig[target_gate] = False
            adjudicable = np.isfinite(real_lfc) & (real_lfc != 0)
            mutable = panel_gate_mask & adjudicable
            nonpanel = ~panel_gate_mask
            fixed_calls = pred_sig & adjudicable & nonpanel
            fixed_correct = fixed_calls & np.isfinite(pred_lfc) & (pred_lfc != 0) & (np.sign(pred_lfc) == np.sign(real_lfc))
            nr = int(real_sig.sum())
            npred = int((pred_sig & adjudicable).sum())
            correct = int(
                (
                    pred_sig
                    & adjudicable
                    & np.isfinite(pred_lfc)
                    & (pred_lfc != 0)
                    & (np.sign(pred_lfc) == np.sign(real_lfc))
                ).sum()
            )
            incumbent_fid = correct / max(npred, nr, 1)
            oracle_fid, oracle_calls = fidelity_oracle(
                int(fixed_calls.sum()),
                int(fixed_correct.sum()),
                nr,
                int(mutable.sum()),
            )
            incumbent_union = int((real_sig | pred_sig).sum())
            incumbent_jac = (
                float((real_sig & pred_sig).sum() / incumbent_union)
                if incumbent_union
                else 1.0
            )
            fixed_nonpanel_tp = int((real_sig & pred_sig & nonpanel).sum())
            fixed_nonpanel_fp = int((~real_sig & pred_sig & nonpanel).sum())
            oracle_jac = (
                float((fixed_nonpanel_tp + int((real_sig & panel_gate_mask).sum())) / (nr + fixed_nonpanel_fp))
                if nr + fixed_nonpanel_fp
                else 1.0
            )

            pool = np.flatnonzero(real_sig)
            order = prediction_order(
                pred_sig[pool],
                pred_padj[pool],
                pred_p[pool],
                pred_lfc[pool],
                ref.gate_genes[pool],
            )
            incumbent_matches = (
                np.isfinite(pred_lfc[pool])
                & (pred_lfc[pool] != 0)
                & (np.sign(pred_lfc[pool]) == np.sign(real_lfc[pool]))
            )[order]
            incumbent_reach = deepest_reach(incumbent_matches)
            nonpanel_pool = pool[~panel_gate_mask[pool]]
            nonpanel_order = prediction_order(
                pred_sig[nonpanel_pool],
                pred_padj[nonpanel_pool],
                pred_p[nonpanel_pool],
                pred_lfc[nonpanel_pool],
                ref.gate_genes[nonpanel_pool],
            )
            nonpanel_matches = (
                np.isfinite(pred_lfc[nonpanel_pool])
                & (pred_lfc[nonpanel_pool] != 0)
                & (np.sign(pred_lfc[nonpanel_pool]) == np.sign(real_lfc[nonpanel_pool]))
            )[nonpanel_order]
            oracle_matches = np.concatenate(
                [np.ones(int((real_sig & panel_gate_mask).sum()), dtype=bool), nonpanel_matches]
            )
            oracle_reach = deepest_reach(oracle_matches)

            record.update(
                {
                    "incumbent_fidelity": float(incumbent_fid),
                    "oracle_fidelity": float(oracle_fid),
                    "oracle_fidelity_raw_gain": float(oracle_fid - incumbent_fid),
                    "oracle_fidelity_panel_calls": int(oracle_calls),
                    "incumbent_jaccard": incumbent_jac,
                    "oracle_jaccard": oracle_jac,
                    "oracle_jaccard_raw_gain": float(oracle_jac - incumbent_jac),
                    "incumbent_reach": incumbent_reach,
                    "oracle_reach": oracle_reach,
                    "oracle_reach_raw_gain": float(oracle_reach - incumbent_reach),
                }
            )
        rows.append(record)

    def mean(key: str) -> float | None:
        values = [float(row[key]) for row in rows if row.get(key) is not None]
        return float(np.mean(values)) if values else None

    summary = {
        "surface": name,
        "evidence_class": "reused_development_evidence_not_independent_confirmation",
        "targets": len(rows),
        "nmae_scored_targets": sum(row["incumbent_nmae"] is not None for row in rows),
        "reference_de_panel_fraction_macro": mean("reference_de_panel_fraction"),
        "reference_abs_lfc_panel_fraction_macro": mean("reference_abs_lfc_panel_fraction"),
        "zero_to_perfect_panel_raw_nmae_ceiling_macro": mean(
            "zero_to_perfect_panel_raw_nmae_ceiling"
        ),
        "incumbent_nmae_macro": mean("incumbent_nmae"),
        "oracle_perfect_panel_nmae_macro": mean("oracle_perfect_panel_nmae"),
        "oracle_perfect_panel_raw_nmae_gain_macro": mean("oracle_perfect_panel_raw_nmae_gain"),
        "oracle_scaled_nmae_gain_macro": mean("oracle_scaled_nmae_gain"),
        "exact_prediction_de_available": exact is not None,
        "incumbent_fidelity_macro": mean("incumbent_fidelity"),
        "oracle_fidelity_macro": mean("oracle_fidelity"),
        "oracle_fidelity_raw_gain_macro": mean("oracle_fidelity_raw_gain"),
        "incumbent_jaccard_macro": mean("incumbent_jaccard"),
        "oracle_jaccard_macro": mean("oracle_jaccard"),
        "oracle_jaccard_raw_gain_macro": mean("oracle_jaccard_raw_gain"),
        "incumbent_reach_macro": mean("incumbent_reach"),
        "oracle_reach_macro": mean("oracle_reach"),
        "oracle_reach_raw_gain_macro": mean("oracle_reach_raw_gain"),
    }
    return summary, rows


def safe_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): safe_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe_json(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis", action="store_true")
    parser.add_argument("--output", type=Path, default=HERE / "work" / "analysis-v1")
    args = parser.parse_args()
    require(args.analysis, "run requires --analysis")
    require(not args.output.exists(), f"refusing to overwrite {args.output}")
    require(args.output.resolve().is_relative_to(HERE.resolve()), "output must stay inside Experiment 127")

    protocol = load_json(PROTOCOL_PATH)
    verify_freeze(protocol)
    panel = read_panel(ROOT / protocol["shared_inputs"]["panel_targets"]["path"])
    h1_scale = fit_h1_nmae_scale(protocol)
    surfaces = []
    per_target: list[dict[str, Any]] = []
    for name, config in protocol["surfaces"].items():
        surface, rows = surface_analysis(
            name,
            config,
            panel,
            h1_scale if name == "h1" else None,
        )
        surfaces.append(surface)
        per_target.extend(rows)

    def surface_mean(key: str) -> float | None:
        values = [surface[key] for surface in surfaces if surface[key] is not None]
        return float(np.mean(values)) if values else None

    summary = {
        "schema_version": 1,
        "experiment": "127-panel-metric-ceiling",
        "status": "analysis_complete",
        "evidence_class": "all_surfaces_are_reused_development_evidence",
        "surfaces": surfaces,
        "three_surface_equal_weight_macro": {
            "reference_de_panel_fraction": surface_mean("reference_de_panel_fraction_macro"),
            "reference_abs_lfc_panel_fraction": surface_mean("reference_abs_lfc_panel_fraction_macro"),
            "oracle_perfect_panel_raw_nmae_gain": surface_mean("oracle_perfect_panel_raw_nmae_gain_macro"),
        },
        "h1_scaled_nmae_calibration": h1_scale,
        "pan_exclusion": {
            "metric_ceiling_computed": False,
            "reason": protocol["pan_handling"]["reason_not_in_metric_ceiling"],
        },
        "interpretation": {
            "oracle_is_not_buildable": True,
            "buildable_transport_deferred_to_experiment_128": True,
            "experiment_124_used": False,
            "no_new_held_outcome_accessed": True,
            "cell_eval2_run": False,
            "counts_generated": False,
            "submission_built_or_sent": False,
        },
    }

    stage = args.output.parent / f".{args.output.name}.tmp-{os.getpid()}"
    if stage.exists():
        shutil.rmtree(stage)
    stage.mkdir(parents=True)
    summary_path = stage / "summary.json"
    per_target_path = stage / "per-target.csv"
    summary_path.write_text(json.dumps(safe_json(summary), indent=2, sort_keys=True) + "\n")
    write_csv(per_target_path, per_target)
    manifest = {
        "schema_version": 1,
        "experiment": "127-panel-metric-ceiling",
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "implementation_freeze_sha256": sha256_file(FREEZE_PATH),
        "files": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for path in (summary_path, per_target_path)
        },
    }
    manifest_path = stage / "artifact-manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    os.replace(stage, args.output)
    print(json.dumps(summary["three_surface_equal_weight_macro"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

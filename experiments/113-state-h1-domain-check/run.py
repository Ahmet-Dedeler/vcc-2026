#!/usr/bin/env python3
"""In-domain H1 diagnostic for the frozen experiment-106 STATE checkpoint."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import sys
import tempfile
from pathlib import Path
from typing import Sequence

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CHECKPOINT_DIR = ROOT / "experiments/106-state-baseline/frozen/state_sm-step20000"
CHECKPOINT_MANIFEST = CHECKPOINT_DIR / "checkpoint-manifest.json"
ADAPTER_PATH = ROOT / "experiments/106-state-baseline/state_effect_adapter.py"
BASELINE_PATH = ROOT / "experiments/100-sign-confidence/work/h1/frozen_predictions.npz"
HARNESS_PATH = ROOT / "data/harness/h1-test.h5ad"
TRUTH_DIR = ROOT / "data/donors/h1_2025"
REFDE_PATH = ROOT / "experiments/098-nmae/work/refde_h1.npz"
PANEL_PATH = ROOT / "data/challenge/pert_counts.csv"
WORK = HERE / "work"
FROZEN_PATH = WORK / "frozen.npz"
MANIFEST_PATH = WORK / "manifest.json"
RESULTS_PATH = WORK / "results.json"
CONTROL_LABEL = "non-targeting"
EXPECTED_CONTROLS = 18_400
N_BANKS = 4
SET_LEN = 128
TARGET_SUM = 50_000.0
N_BOOT = 2_000


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_strings(values: Sequence[str]) -> str:
    return hashlib.sha256(("\n".join(map(str, values)) + "\n").encode()).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temp.replace(path)


def atomic_npz(path: Path, **values: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent, delete=False
    ) as handle:
        temp = Path(handle.name)
        np.savez_compressed(handle, **values)
    temp.replace(path)


def read_lines(path: Path) -> tuple[str, ...]:
    return tuple(line.strip() for line in path.read_text().splitlines() if line.strip())


def read_panel() -> set[str]:
    return set(read_lines(PANEL_PATH)[1:])


def load_baseline() -> dict[str, np.ndarray]:
    with np.load(BASELINE_PATH, allow_pickle=False) as archive:
        required = {"genes", "targets", "bctrl", "ctrl_cpm", "dh_identity"}
        if not required.issubset(archive.files):
            raise ValueError("experiment-100 H1 freeze lacks required arrays")
        return {
            "genes": archive["genes"].astype(str),
            "targets": archive["targets"].astype(str),
            "bctrl": archive["bctrl"].astype(np.float64),
            "ctrl_cpm": archive["ctrl_cpm"].astype(np.float64),
            "incumbent": archive["dh_identity"].astype(np.float64),
        }


def freeze() -> None:
    import anndata as ad

    if WORK.exists():
        raise FileExistsError(f"refusing to overwrite {WORK}")
    adapter = load_module("exp106_state_adapter", ADAPTER_PATH)
    checkpoint = adapter.CheckpointSpec.load(CHECKPOINT_MANIFEST)
    baseline = load_baseline()
    state_genes = np.asarray(adapter.read_axis(adapter.STATE_GENE_PATH), dtype="U")
    if not np.array_equal(state_genes, baseline["genes"]):
        raise ValueError("H1 baseline axis differs from the frozen STATE axis")
    adapter.validate_runtime_provenance(checkpoint, CHECKPOINT_DIR, state_genes)

    handle = ad.read_h5ad(HARNESS_PATH, backed="r")
    try:
        if not np.array_equal(np.asarray(handle.var_names).astype(str), state_genes):
            raise ValueError("H1 harness axis differs from the frozen STATE axis")
        labels = np.asarray(handle.obs["target_gene"].astype(str))
        control_rows = np.flatnonzero(labels == CONTROL_LABEL)
        if len(control_rows) != EXPECTED_CONTROLS:
            raise ValueError(
                f"H1 control count {len(control_rows)} != {EXPECTED_CONTROLS}"
            )
        selected = adapter.sample_control_sets(
            control_rows,
            "H1",
            n_sets=N_BANKS,
            set_len=SET_LEN,
        )
        raw = handle.X[selected.reshape(-1), :]
        controls = (
            adapter.GeneAlignment.build(state_genes, state_genes)
            .align_and_log_controls(raw)
            .reshape(N_BANKS, SET_LEN, len(state_genes))
        )
        control_sum = adapter.sum_control_counts(
            handle.X, control_rows, len(state_genes)
        )
        control_mean = control_sum / len(control_rows)
    finally:
        handle.file.close()

    computed_bctrl = adapter.bulk_profile(control_mean)
    baseline_error = float(np.max(np.abs(computed_bctrl - baseline["bctrl"])))
    if baseline_error > 2e-5:
        raise ValueError(f"H1 baseline mismatch: max_abs={baseline_error}")

    model, pert_map, torch = adapter.load_runtime(
        checkpoint, CHECKPOINT_DIR, state_genes, "mps"
    )
    embeddings = adapter.load_target_embeddings(baseline["targets"], pert_map)
    forbidden_symbols = read_panel() | set(map(str, baseline["targets"]))
    forbidden = np.asarray(
        [gene in forbidden_symbols for gene in state_genes], dtype=bool
    )
    bank_effects = np.zeros(
        (N_BANKS, len(baseline["targets"]), len(state_genes)), dtype=np.float32
    )
    pooled_means = np.zeros(
        (len(baseline["targets"]), len(state_genes)), dtype=np.float64
    )
    output_min = np.inf
    output_max = -np.inf
    clipped_high = 0
    clipped_low = 0
    device = next(model.parameters()).device

    with torch.inference_mode():
        for row, target in enumerate(baseline["targets"]):
            embedding = embeddings[row]
            means = []
            for bank in range(N_BANKS):
                control_tensor = torch.from_numpy(controls[bank]).to(device)
                perturbation = adapter.build_perturbation_batch(embedding, 1, SET_LEN)
                prediction = model(
                    {
                        "ctrl_cell_emb": control_tensor,
                        "pert_emb": torch.from_numpy(perturbation).to(device),
                    },
                    padded=True,
                )
                values = (
                    prediction.detach().cpu().numpy().astype(np.float32, copy=False)
                )
                output_min = min(output_min, float(values.min()))
                output_max = max(output_max, float(values.max()))
                clipped_low += int(np.count_nonzero(values < 0.0))
                clipped_high += int(
                    np.count_nonzero(values > adapter.MAX_LOG_EXPRESSION)
                )
                mean_counts = adapter.streamed_mean_counts([values])
                means.append(mean_counts)
                effect = adapter.scorer_effect(mean_counts, control_mean)
                effect[forbidden] = 0.0
                bank_effects[bank, row] = effect.astype(np.float32)
            pooled_mean = np.mean(np.stack(means), axis=0)
            pooled_means[row] = pooled_mean

    pooled_effect = adapter.scorer_effect(pooled_means, control_mean)
    pooled_effect[:, forbidden] = 0.0
    mean_effect = np.mean(bank_effects.astype(np.float64), axis=0)
    mean_effect[:, forbidden] = 0.0
    incumbent = baseline["incumbent"].copy()
    incumbent[:, forbidden] = 0.0

    atomic_npz(
        FROZEN_PATH,
        genes=state_genes,
        targets=baseline["targets"],
        bctrl=baseline["bctrl"].astype(np.float32),
        ctrl_cpm=baseline["ctrl_cpm"].astype(np.float32),
        incumbent=incumbent.astype(np.float32),
        state_banks=bank_effects,
        state_mean=mean_effect.astype(np.float32),
        state_pooled=pooled_effect.astype(np.float32),
        sampled_control_rows=selected.astype(np.int64),
        forbidden_gene_mask=forbidden,
    )
    manifest = {
        "experiment": "113-state-h1-domain-check",
        "stage": "truth-free-inference-freeze",
        "sealed_before_h1_truth": True,
        "checkpoint_manifest": str(CHECKPOINT_MANIFEST.relative_to(ROOT)),
        "checkpoint_manifest_sha256": sha256_file(CHECKPOINT_MANIFEST),
        "checkpoint_sha256": checkpoint.sha256,
        "checkpoint_global_step": checkpoint.global_step,
        "checkpoint_profile": checkpoint.profile,
        "adapter_source": str(ADAPTER_PATH.relative_to(ROOT)),
        "adapter_source_sha256": sha256_file(ADAPTER_PATH),
        "implementation_sha256": sha256_file(Path(__file__)),
        "readme_sha256": sha256_file(HERE / "README.md"),
        "baseline_source": str(BASELINE_PATH.relative_to(ROOT)),
        "baseline_source_sha256": sha256_file(BASELINE_PATH),
        "harness": str(HARNESS_PATH.relative_to(ROOT)),
        "harness_sha256": sha256_file(HARNESS_PATH),
        "expression_rows_read": "non-targeting controls only",
        "control_rows": len(control_rows),
        "control_banks": N_BANKS,
        "control_set_len": SET_LEN,
        "sampled_control_rows_sha256": hashlib.sha256(selected.tobytes()).hexdigest(),
        "genes": len(state_genes),
        "gene_axis_sha256": sha256_strings(state_genes),
        "targets": len(baseline["targets"]),
        "target_axis_sha256": sha256_strings(baseline["targets"]),
        "forbidden_symbols": int(forbidden.sum()),
        "forbidden_policy": "union of VCC-2026 panel and H1 test target genes",
        "normalization": "log1p raw counts input; expm1 and mean in count space; 50k log-pseudobulk",
        "incumbent_bctrl_max_abs_diff": baseline_error,
        "output_bounds": {
            "observed_min": output_min,
            "observed_max": output_max,
            "clipped_low_values": clipped_low,
            "clipped_high_values": clipped_high,
            "clip_min": 0.0,
            "clip_max": adapter.MAX_LOG_EXPRESSION,
        },
    }
    manifest["frozen_sha256"] = sha256_file(FROZEN_PATH)
    atomic_json(MANIFEST_PATH, manifest)
    print(MANIFEST_PATH)


def pairwise_pds(
    predicted: np.ndarray, truth: np.ndarray, mask: np.ndarray
) -> np.ndarray:
    from cell_eval2.distances import pairwise_full

    distances = pairwise_full(predicted[:, mask], truth[:, mask], "cosine")
    own = np.diag(distances)[:, None]
    ranks = (distances < own).sum(axis=1)
    ties = (distances == own).sum(axis=1) - 1
    return 1.0 - (ranks + 0.5 * ties) / (len(predicted) - 1)


def cosine(predicted: np.ndarray, truth: np.ndarray, mask: np.ndarray) -> np.ndarray:
    left = predicted[:, mask]
    right = truth[:, mask]
    return np.einsum("ij,ij->i", left, right) / np.maximum(
        np.linalg.norm(left, axis=1) * np.linalg.norm(right, axis=1), 1e-300
    )


def nmse(predicted: np.ndarray, truth: np.ndarray, mask: np.ndarray) -> np.ndarray:
    return np.square(predicted[:, mask] - truth[:, mask]).sum(axis=1) / np.maximum(
        np.square(truth[:, mask]).sum(axis=1), 1e-300
    )


def bootstrap(values: np.ndarray, label: str) -> dict[str, float]:
    values = np.asarray(values, dtype=np.float64)
    rng = np.random.default_rng(
        int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "little")
    )
    draws = values[rng.integers(0, len(values), (N_BOOT, len(values)))].mean(axis=1)
    return {
        "estimate": float(values.mean()),
        "lo": float(np.percentile(draws, 2.5)),
        "hi": float(np.percentile(draws, 97.5)),
    }


def load_truth(targets: np.ndarray, genes: np.ndarray) -> np.ndarray:
    truth_targets = np.asarray(read_lines(TRUTH_DIR / "perts.txt"))
    truth_genes = np.asarray(read_lines(TRUTH_DIR / "genes.txt"))
    target_index = {value: index for index, value in enumerate(truth_targets)}
    gene_index = {value: index for index, value in enumerate(truth_genes)}
    if any(target not in target_index for target in targets):
        raise ValueError("H1 donor truth lacks a frozen target")
    if any(gene not in gene_index for gene in genes):
        raise ValueError("H1 donor truth lacks a frozen gene")
    rows = np.asarray([target_index[target] for target in targets], dtype=np.int64)
    cols = np.asarray([gene_index[gene] for gene in genes], dtype=np.int64)
    source = np.load(TRUTH_DIR / "delta.npy", mmap_mode="r")
    return np.asarray(source[np.ix_(rows, cols)], dtype=np.float64)


def directional_credit(
    predicted_lfc: np.ndarray,
    reference_lfc: np.ndarray,
    masks: Sequence[np.ndarray],
) -> np.ndarray:
    values = np.full(len(predicted_lfc), np.nan, dtype=np.float64)
    for row, mask in enumerate(masks):
        truth = reference_lfc[row, mask]
        denominator = float(np.abs(truth).sum())
        if denominator > 0:
            values[row] = float(np.sign(truth) @ predicted_lfc[row, mask]) / denominator
    return values


def score() -> None:
    if not FROZEN_PATH.is_file() or not MANIFEST_PATH.is_file():
        raise FileNotFoundError("run freeze before score")
    if RESULTS_PATH.exists():
        raise FileExistsError(f"refusing to overwrite {RESULTS_PATH}")
    manifest = json.loads(MANIFEST_PATH.read_text())
    if manifest.get("sealed_before_h1_truth") is not True:
        raise ValueError("inference was not sealed before H1 truth")
    if sha256_file(FROZEN_PATH) != manifest.get("frozen_sha256"):
        raise ValueError("frozen inference artifact changed")
    if sha256_file(Path(__file__)) != manifest.get("implementation_sha256"):
        raise ValueError("implementation changed after inference freeze")
    if sha256_file(HERE / "README.md") != manifest.get("readme_sha256"):
        raise ValueError("preregistration changed after inference freeze")

    with np.load(FROZEN_PATH, allow_pickle=False) as archive:
        genes = archive["genes"].astype(str)
        targets = archive["targets"].astype(str)
        bctrl = archive["bctrl"].astype(np.float64)
        incumbent = archive["incumbent"].astype(np.float64)
        state_banks = archive["state_banks"].astype(np.float64)
        state_mean = archive["state_mean"].astype(np.float64)
        state_pooled = archive["state_pooled"].astype(np.float64)
        forbidden = archive["forbidden_gene_mask"].astype(bool)
    if np.any(incumbent[:, forbidden]) or np.any(state_banks[:, :, forbidden]):
        raise ValueError("forbidden coordinate is nonzero")

    truth_hashes = {
        name: sha256_file(TRUTH_DIR / name)
        for name in ("genes.txt", "perts.txt", "delta.npy")
    }
    refde_hash = sha256_file(REFDE_PATH)
    truth = load_truth(targets, genes)
    score_mask = ~forbidden

    ref = np.load(REFDE_PATH, allow_pickle=False)
    ref_genes = ref["genes"].astype(str)
    if not np.array_equal(ref_genes, genes):
        raise ValueError("H1 reference-DE axis differs")
    gate = ref["gate"].astype(bool)
    gate_idx = np.flatnonzero(gate)
    ref_targets = ref["perts"].astype(str)
    ref_rows = {target: index for index, target in enumerate(ref_targets)}
    reference_lfc_all = ref["lfc"].astype(np.float64)
    sig_ptr = ref["sig_ptr"].astype(np.int64)
    sig_idx = ref["sig_idx"].astype(np.int64)
    ref.close()

    def effect_to_lfc(effect: np.ndarray) -> np.ndarray:
        profile = np.expm1(bctrl[None, :] + effect) / TARGET_SUM
        np.clip(profile, 0.0, None, out=profile)
        ref_lin = np.load(REFDE_PATH, allow_pickle=False)["ref_lin"].astype(np.float64)
        return np.log2((1e6 * profile[:, gate] + 1e-9) / (ref_lin[None, :] + 1e-9))

    candidates = {
        "incumbent": incumbent,
        "state_mean": state_mean,
        "state_pooled": state_pooled,
    }
    for bank in range(N_BANKS):
        candidates[f"state_bank_{bank}"] = state_banks[bank]
    active = score_mask
    inc_norm = np.linalg.norm(incumbent[:, active], axis=1)
    state_norm = np.linalg.norm(state_pooled[:, active], axis=1)
    state_matched = state_pooled * (inc_norm / np.maximum(state_norm, 1e-300))[:, None]
    for weight in (0.01, 0.02, 0.05, 0.10):
        mixed = (1.0 - weight) * incumbent + weight * state_matched
        mixed *= (
            inc_norm / np.maximum(np.linalg.norm(mixed[:, active], axis=1), 1e-300)
        )[:, None]
        mixed[:, forbidden] = 0.0
        candidates[f"direction_w{weight:.2f}"] = mixed

    per_target: dict[str, dict[str, np.ndarray]] = {}
    direction_masks = []
    reference_rows = []
    for target in targets:
        row = ref_rows[target]
        sig = sig_idx[sig_ptr[row] : sig_ptr[row + 1]]
        full = gate_idx[sig]
        keep = (
            (genes[full] != target)
            & np.isfinite(reference_lfc_all[row, sig])
            & (reference_lfc_all[row, sig] != 0)
        )
        direction_masks.append(sig[keep])
        reference_rows.append(reference_lfc_all[row])
    reference_lfc = np.stack(reference_rows)

    for name, effect in candidates.items():
        lfc = effect_to_lfc(effect)
        per_target[name] = {
            "pds": pairwise_pds(effect, truth, score_mask),
            "cosine": cosine(effect, truth, score_mask),
            "nmse": nmse(effect, truth, score_mask),
            "credit": directional_credit(lfc, reference_lfc, direction_masks),
        }

    inc_lfc = effect_to_lfc(incumbent)
    state_lfc = effect_to_lfc(state_pooled)
    gate_active = ~forbidden[gate]
    inc_lfc_norm = np.linalg.norm(inc_lfc[:, gate_active], axis=1)
    state_lfc_norm = np.linalg.norm(state_lfc[:, gate_active], axis=1)
    state_lfc_matched = (
        state_lfc * (inc_lfc_norm / np.maximum(state_lfc_norm, 1e-300))[:, None]
    )
    oracle = inc_lfc.copy()
    useful_pairs = 0
    disagreement_pairs = 0
    for row, sig in enumerate(direction_masks):
        truth_lfc = reference_lfc[row, sig]
        a = inc_lfc[row, sig]
        s = state_lfc_matched[row, sig]
        disagreement = np.sign(a) != np.sign(s)
        useful = (
            disagreement
            & (np.sign(s) == np.sign(truth_lfc))
            & (np.sign(a) != np.sign(truth_lfc))
        )
        oracle[row, sig[useful]] = s[useful]
        disagreement_pairs += int(disagreement.sum())
        useful_pairs += int(useful.sum())
    oracle_credit = directional_credit(oracle, reference_lfc, direction_masks)
    incumbent_credit = per_target["incumbent"]["credit"]

    baseline = per_target["incumbent"]
    result_candidates = {}
    for name, values in per_target.items():
        result_candidates[name] = {
            "means": {
                metric: float(np.nanmean(vector)) for metric, vector in values.items()
            },
            "delta_vs_incumbent": {
                metric: bootstrap(
                    vector - baseline[metric],
                    f"experiment-113:{name}:{metric}",
                )
                for metric, vector in values.items()
            },
        }

    result = {
        "experiment": "113-state-h1-domain-check",
        "status": "diagnostic-only",
        "targets": len(targets),
        "score_genes": int(score_mask.sum()),
        "pds_mask": "VCC-2026 panel union H1-test targets",
        "candidates": result_candidates,
        "oracle_disagreement_routing": {
            "credit_delta": bootstrap(
                oracle_credit - incumbent_credit,
                "experiment-113:oracle-credit",
            ),
            "disagreement_pairs": disagreement_pairs,
            "useful_pairs": useful_pairs,
            "useful_fraction": useful_pairs / max(disagreement_pairs, 1),
            "state_lfc_row_norm_matched": True,
        },
        "interpretation_thresholds": {
            "state_pds_minimum_for_learned_signal": 0.55,
            "experiment_110_oracle_credit_each_primary_minimum": 0.015,
        },
        "provenance": {
            "frozen_sha256": sha256_file(FROZEN_PATH),
            "manifest_sha256": sha256_file(MANIFEST_PATH),
            "truth": truth_hashes,
            "reference_de_sha256": refde_hash,
            "implementation_sha256": sha256_file(Path(__file__)),
            "readme_sha256": sha256_file(HERE / "README.md"),
        },
    }
    atomic_json(RESULTS_PATH, result)
    print(RESULTS_PATH)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("freeze", "score"))
    args = parser.parse_args()
    if args.command == "freeze":
        freeze()
    else:
        score()


if __name__ == "__main__":
    main()

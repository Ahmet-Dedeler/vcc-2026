#!/usr/bin/env python3
"""Repair NaN bootstrap reporting in experiment 113 without changing inference."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
RUN_PATH = HERE / "run.py"
RESULTS_V1 = HERE / "work/results.json"
RESULTS_V2 = HERE / "work/results-v2.json"


def load_run_module():
    spec = importlib.util.spec_from_file_location("experiment113_run", RUN_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(RUN_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules["experiment113_run"] = module
    spec.loader.exec_module(module)
    return module


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def finite_bootstrap(values: np.ndarray, label: str) -> dict[str, float | int]:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        raise ValueError(f"{label}: no finite values")
    rng = np.random.default_rng(
        int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "little")
    )
    draws = values[rng.integers(0, len(values), (2_000, len(values)))].mean(axis=1)
    return {
        "estimate": float(values.mean()),
        "lo": float(np.percentile(draws, 2.5)),
        "hi": float(np.percentile(draws, 97.5)),
        "n": int(len(values)),
    }


def main() -> None:
    if RESULTS_V2.exists():
        raise FileExistsError(f"refusing to overwrite {RESULTS_V2}")
    run = load_run_module()
    manifest = json.loads(run.MANIFEST_PATH.read_text())
    if sha256_file(run.FROZEN_PATH) != manifest["frozen_sha256"]:
        raise ValueError("frozen inference artifact changed")
    if sha256_file(RUN_PATH) != manifest["implementation_sha256"]:
        raise ValueError("inference implementation changed after freeze")

    with np.load(run.FROZEN_PATH, allow_pickle=False) as archive:
        genes = archive["genes"].astype(str)
        targets = archive["targets"].astype(str)
        bctrl = archive["bctrl"].astype(np.float64)
        incumbent = archive["incumbent"].astype(np.float64)
        state_banks = archive["state_banks"].astype(np.float64)
        state_mean = archive["state_mean"].astype(np.float64)
        state_pooled = archive["state_pooled"].astype(np.float64)
        forbidden = archive["forbidden_gene_mask"].astype(bool)

    with np.load(run.REFDE_PATH, allow_pickle=False) as ref:
        gate = ref["gate"].astype(bool)
        gate_idx = np.flatnonzero(gate)
        ref_targets = ref["perts"].astype(str)
        reference_lfc_all = ref["lfc"].astype(np.float64)
        sig_ptr = ref["sig_ptr"].astype(np.int64)
        sig_idx = ref["sig_idx"].astype(np.int64)
        ref_lin = ref["ref_lin"].astype(np.float64)

    def effect_to_lfc(effect: np.ndarray) -> np.ndarray:
        profile = np.expm1(bctrl[None, :] + effect) / run.TARGET_SUM
        np.clip(profile, 0.0, None, out=profile)
        return np.log2((1e6 * profile[:, gate] + 1e-9) / (ref_lin[None, :] + 1e-9))

    ref_rows = {target: index for index, target in enumerate(ref_targets)}
    direction_masks: list[np.ndarray] = []
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

    candidates = {
        "incumbent": incumbent,
        "state_mean": state_mean,
        "state_pooled": state_pooled,
    }
    for bank in range(run.N_BANKS):
        candidates[f"state_bank_{bank}"] = state_banks[bank]
    active = ~forbidden
    incumbent_norm = np.linalg.norm(incumbent[:, active], axis=1)
    state_norm = np.linalg.norm(state_pooled[:, active], axis=1)
    state_matched = (
        state_pooled * (incumbent_norm / np.maximum(state_norm, 1e-300))[:, None]
    )
    for weight in (0.01, 0.02, 0.05, 0.10):
        mixed = (1.0 - weight) * incumbent + weight * state_matched
        mixed *= (
            incumbent_norm
            / np.maximum(np.linalg.norm(mixed[:, active], axis=1), 1e-300)
        )[:, None]
        mixed[:, forbidden] = 0.0
        candidates[f"direction_w{weight:.2f}"] = mixed

    credit = {
        name: run.directional_credit(
            effect_to_lfc(effect), reference_lfc, direction_masks
        )
        for name, effect in candidates.items()
    }
    incumbent_credit = credit["incumbent"]

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
        incumbent_values = inc_lfc[row, sig]
        state_values = state_lfc_matched[row, sig]
        disagreement = np.sign(incumbent_values) != np.sign(state_values)
        useful = (
            disagreement
            & (np.sign(state_values) == np.sign(truth_lfc))
            & (np.sign(incumbent_values) != np.sign(truth_lfc))
        )
        oracle[row, sig[useful]] = state_values[useful]
        disagreement_pairs += int(disagreement.sum())
        useful_pairs += int(useful.sum())
    oracle_credit = run.directional_credit(oracle, reference_lfc, direction_masks)

    result = json.loads(RESULTS_V1.read_text())
    for name, values in credit.items():
        result["candidates"][name]["delta_vs_incumbent"]["credit"] = finite_bootstrap(
            values - incumbent_credit,
            f"experiment-113-v2:{name}:credit",
        )
    result["oracle_disagreement_routing"] = {
        "credit_delta": finite_bootstrap(
            oracle_credit - incumbent_credit,
            "experiment-113-v2:oracle-credit",
        ),
        "disagreement_pairs": disagreement_pairs,
        "useful_pairs": useful_pairs,
        "useful_fraction": useful_pairs / max(disagreement_pairs, 1),
        "state_lfc_row_norm_matched": True,
    }
    result["reporting_correction"] = {
        "reason": (
            "results.json bootstrapped NaN placeholders for three targets with "
            "no eligible reference-DE coordinates; results-v2 drops non-finite "
            "targets exactly as the reported nanmean already did"
        ),
        "inference_unchanged": True,
        "results_v1_sha256": sha256_file(RESULTS_V1),
        "rescore_source_sha256": sha256_file(Path(__file__)),
    }
    run.atomic_json(RESULTS_V2, result)
    print(RESULTS_V2)


if __name__ == "__main__":
    main()

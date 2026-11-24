#!/usr/bin/env python3
"""Score frozen STATE/incumbent effects without reopening inference choices."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from cell_eval2.distances import pairwise_full


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PANEL = tuple(
    row.strip()
    for row in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:]
    if row.strip()
)
PANEL_SET = set(PANEL)
TS = 50_000.0
MIN_CPM = 10.0
N_BOOT = 2_000
WEIGHTS = (0.10, 0.20, 0.35, 0.50, 0.65, 0.80)
EXPECTED_HELD_TARGET_AXES = {
    "hct116": {
        "count": 132,
        "sha256": "476f8fa9f34fe0865fe135eda61e0b3bb979bd5b4725fbcf2436dbc4cc68b17b",
    },
    "hek293t": {
        "count": 132,
        "sha256": "d4fbc4d5e296f9eb73a87a4df84e01ae91aead80aa52e110bd166782f87fccca",
    },
}
EXPECTED_COMMON_TARGET_COUNT = 120
EXPECTED_ORION_SCORE_GENES = 17_806
EXPECTED_COMMON_TARGET_SHA256 = (
    "ecbfb06f75a6c275cd1f06f4ea10a13cf4512ca279ef507b7bc0aca0fbbfa90a"
)

SURFACES = {
    "hct116": {
        "truth": ROOT / "data/donors/orion_hct116",
        "refde": ROOT / "experiments/098-nmae/work/refde_orion.npz",
    },
    "hek293t": {
        "truth": ROOT / "data/donors/orion_hek293t",
        "refde": ROOT / "experiments/099-donors/work/refde_hek293t.npz",
    },
}


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_strings(values: Sequence[str]) -> str:
    return hashlib.sha256(("\n".join(map(str, values)) + "\n").encode()).hexdigest()


def require_finite(name: str, values: np.ndarray, shape: tuple[int, ...]) -> None:
    array = np.asarray(values)
    if array.shape != shape:
        raise ValueError(f"{name} shape {array.shape} != {shape}")
    if not np.isfinite(array).all():
        raise ValueError(f"{name} contains non-finite values")


def truth_bundle_hashes(path: Path) -> dict[str, str]:
    required = ("genes.txt", "perts.txt", "delta.npy")
    result = {}
    for name in required:
        item = path / name
        if not item.is_file():
            raise FileNotFoundError(item)
        result[name] = sha256_file(item)
    return result


def stable_seed(*parts: object) -> int:
    return int.from_bytes(
        hashlib.sha256("\0".join(map(str, parts)).encode()).digest()[:8], "little"
    )


def atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def read_lines(path: Path) -> tuple[str, ...]:
    return tuple(line.strip() for line in path.read_text().splitlines() if line.strip())


def last_occurrence_lookup(
    values: Sequence[str],
) -> tuple[dict[str, int], dict[str, object]]:
    lookup: dict[str, int] = {}
    duplicate_symbols: set[str] = set()
    for index, value in enumerate(map(str, values)):
        if value in lookup:
            duplicate_symbols.add(value)
        lookup[value] = index
    return lookup, {
        "rows": len(values),
        "unique_symbols": len(lookup),
        "duplicate_rows": len(values) - len(lookup),
        "duplicate_symbols": sorted(duplicate_symbols),
        "duplicate_policy": "last occurrence",
    }


def validate_common_target_axis(
    targets: Sequence[str],
    *,
    expected_count: int = EXPECTED_COMMON_TARGET_COUNT,
    expected_sha256: str = EXPECTED_COMMON_TARGET_SHA256,
) -> None:
    values = tuple(map(str, targets))
    if len(values) != expected_count:
        raise ValueError(
            f"held common gallery has {len(values)} targets; expected {expected_count}"
        )
    if sha256_strings(values) != expected_sha256:
        raise ValueError("held common target axis changed")


def verify_manifest_file(
    manifest: dict[str, object], path_key: str, hash_key: str
) -> Path:
    path = (ROOT / str(manifest[path_key])).resolve()
    if ROOT not in path.parents:
        raise ValueError(f"manifest path escapes repository: {path}")
    if not path.is_file():
        raise FileNotFoundError(path)
    if sha256_file(path) != manifest.get(hash_key):
        raise ValueError(f"{path_key} differs from the frozen inference manifest")
    return path


@dataclass
class ReferenceDE:
    genes: np.ndarray
    gate: np.ndarray
    gate_idx: np.ndarray
    ref_lin: np.ndarray
    perts: np.ndarray
    lfc: np.ndarray
    sig: list[np.ndarray]

    @classmethod
    def load(cls, path: Path) -> "ReferenceDE":
        with np.load(path, allow_pickle=False) as z:
            genes = z["genes"].astype(str)
            gate = z["gate"].astype(bool)
            gate_idx = np.flatnonzero(gate)
            ref_lin = z["ref_lin"].astype(np.float64)
            perts = z["perts"].astype(str)
            lfc = z["lfc"].astype(np.float64)
            ptr = z["sig_ptr"].astype(np.int64)
            idx = z["sig_idx"].astype(np.int64)
        if len(gate_idx) != len(ref_lin):
            raise ValueError("reference gate and linear baseline differ")
        if len(genes) != len(gate):
            raise ValueError("reference gene axis and gate differ")
        if lfc.shape != (len(perts), len(gate_idx)):
            raise ValueError("reference LFC shape is malformed")
        if len(ptr) != len(perts) + 1 or ptr[0] != 0 or ptr[-1] != len(idx):
            raise ValueError("reference significance pointers are malformed")
        if np.any(np.diff(ptr) < 0) or np.any(idx < 0) or np.any(idx >= len(gate_idx)):
            raise ValueError("reference significance indices are malformed")
        if not np.isfinite(ref_lin).all() or np.any(ref_lin < 0):
            raise ValueError("reference linear baseline is invalid")
        sig = [idx[ptr[i] : ptr[i + 1]] for i in range(len(perts))]
        for row, selected in enumerate(sig):
            if not np.isfinite(lfc[row, selected]).all():
                raise ValueError("reference significant LFC values are non-finite")
        return cls(genes, gate, gate_idx, ref_lin, perts, lfc, sig)


def effect_to_profile(effect: np.ndarray, bctrl: np.ndarray) -> np.ndarray:
    profile = np.expm1(
        np.asarray(bctrl, dtype=np.float64)[None, :]
        + np.asarray(effect, dtype=np.float64)
    )
    np.clip(profile, 0.0, None, out=profile)
    totals = profile.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise ValueError("effect produced a zero-mass profile")
    return profile / totals


def profile_to_effect(profile: np.ndarray, bctrl: np.ndarray) -> np.ndarray:
    values = np.asarray(profile, dtype=np.float64)
    if values.ndim != 2 or not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("profile must be a finite nonnegative matrix")
    totals = values.sum(axis=1, keepdims=True)
    if np.any(totals <= 0):
        raise ValueError("profile has a zero-mass row")
    normalized = values / totals
    return np.log1p(TS * normalized) - np.asarray(bctrl, dtype=np.float64)[None, :]


def realize_effect(
    effect: np.ndarray,
    incumbent: np.ndarray,
    bctrl: np.ndarray,
    ctrl_cpm: np.ndarray,
    supported: np.ndarray,
    panel_mask: np.ndarray,
) -> np.ndarray:
    gated = np.asarray(effect, dtype=np.float64).copy()
    incumbent = np.asarray(incumbent, dtype=np.float64)
    gate = (np.asarray(ctrl_cpm) > MIN_CPM) & (np.asarray(bctrl) > 0)
    adjustable = np.asarray(supported, dtype=bool) & ~panel_mask & gate
    gated[:, ~adjustable] = incumbent[:, ~adjustable]

    incumbent_gated = incumbent.copy()
    incumbent_gated[:, ~gate] = 0.0
    incumbent_profile = effect_to_profile(incumbent_gated, bctrl)
    candidate_profile = effect_to_profile(gated, bctrl)
    output_profile = incumbent_profile.copy()
    desired_mass = incumbent_profile[:, adjustable].sum(axis=1, keepdims=True)
    candidate_mass = candidate_profile[:, adjustable].sum(axis=1, keepdims=True)
    if np.any((desired_mass > 0) & (candidate_mass <= 0)):
        raise ValueError("candidate has no mass on an adjustable gene block")
    output_profile[:, adjustable] = (
        candidate_profile[:, adjustable]
        * desired_mass
        / np.maximum(candidate_mass, 1e-300)
    )
    out = profile_to_effect(output_profile, bctrl)
    implied = np.expm1(np.asarray(bctrl)[None, :] + out).sum(axis=1)
    if not np.allclose(implied, TS, rtol=1e-10, atol=1e-6):
        raise AssertionError("realized effect does not close to the 50k profile")
    return out


def effect_blend(
    incumbent: np.ndarray,
    state: np.ndarray,
    supported: np.ndarray,
    weight: float,
) -> np.ndarray:
    out = np.asarray(incumbent, dtype=np.float64).copy()
    out[:, supported] = (1.0 - weight) * incumbent[:, supported] + weight * state[
        :, supported
    ]
    return out


def direction_blend(
    incumbent: np.ndarray,
    state: np.ndarray,
    supported: np.ndarray,
    active: np.ndarray,
    weight: float,
) -> np.ndarray:
    incumbent = np.asarray(incumbent, dtype=np.float64)
    state = np.asarray(state, dtype=np.float64)
    out = incumbent.copy()
    use = supported & active
    for row in range(len(out)):
        inc_norm = float(np.linalg.norm(incumbent[row, use]))
        state_norm = float(np.linalg.norm(state[row, use]))
        if inc_norm <= 0 or state_norm <= 0:
            continue
        matched = state[row, use] * (inc_norm / state_norm)
        mixed = (1.0 - weight) * incumbent[row, use] + weight * matched
        mixed_norm = float(np.linalg.norm(mixed))
        if mixed_norm > 0:
            mixed *= inc_norm / mixed_norm
        out[row, use] = mixed
    return out


def composition_blend(
    incumbent: np.ndarray,
    state: np.ndarray,
    bctrl: np.ndarray,
    supported: np.ndarray,
    weight: float,
) -> np.ndarray:
    q_inc = effect_to_profile(incumbent, bctrl)
    q_state = effect_to_profile(state, bctrl)
    out = q_inc.copy()
    common_mass = q_inc[:, supported].sum(axis=1, keepdims=True)
    mixed = (1.0 - weight) * q_inc[:, supported] + weight * q_state[:, supported]
    mixed_mass = mixed.sum(axis=1, keepdims=True)
    mixed *= common_mass / np.maximum(mixed_mass, 1e-300)
    out[:, supported] = mixed
    return np.log1p(TS * out) - np.asarray(bctrl, dtype=np.float64)[None, :]


def candidate_effects(
    incumbent: np.ndarray,
    state: np.ndarray,
    bctrl: np.ndarray,
    ctrl_cpm: np.ndarray,
    supported: np.ndarray,
    panel_mask: np.ndarray,
) -> dict[str, np.ndarray]:
    active = ~panel_mask
    state_fallback = incumbent.copy()
    state_fallback[:, supported] = state[:, supported]
    raw = {
        "incumbent": incumbent,
        "state_fallback": state_fallback,
    }
    for weight in WEIGHTS:
        suffix = f"{weight:.2f}"
        raw[f"effect_w{suffix}"] = effect_blend(incumbent, state, supported, weight)
        raw[f"direction_w{suffix}"] = direction_blend(
            incumbent, state, supported, active, weight
        )
        raw[f"composition_w{suffix}"] = composition_blend(
            incumbent, state, bctrl, supported, weight
        )
    realized = {
        name: realize_effect(
            value,
            incumbent,
            bctrl,
            ctrl_cpm,
            supported,
            panel_mask,
        )
        for name, value in raw.items()
    }
    return realized


def load_truth_effects(
    path: Path, targets: Sequence[str], genes: Sequence[str]
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    truth_genes = read_lines(path / "genes.txt")
    truth_targets = read_lines(path / "perts.txt")
    if len(truth_targets) != len(set(truth_targets)):
        raise ValueError("truth target axis contains duplicates")
    target_lookup = {target: index for index, target in enumerate(truth_targets)}
    missing_targets = [
        str(target) for target in targets if str(target) not in target_lookup
    ]
    if missing_targets:
        raise ValueError(f"truth lacks requested targets: {missing_targets[:5]}")
    rows = np.asarray(
        [target_lookup[str(target)] for target in targets], dtype=np.int64
    )
    delta = np.load(path / "delta.npy", mmap_mode="r")
    if delta.shape != (len(truth_targets), len(truth_genes)):
        raise ValueError("truth matrix shape differs from its axes")
    gene_lookup, axis_audit = last_occurrence_lookup(truth_genes)
    official = tuple(map(str, genes))
    if len(official) != len(set(official)):
        raise ValueError("candidate gene axis contains duplicates")
    measured = np.asarray([gene in gene_lookup for gene in official], dtype=bool)
    source_cols = np.asarray(
        [gene_lookup[gene] for gene in official if gene in gene_lookup],
        dtype=np.int64,
    )
    output_cols = np.flatnonzero(measured)
    result = np.zeros((len(targets), len(official)), dtype=np.float64)
    result[:, output_cols] = np.asarray(
        delta[np.ix_(rows, source_cols)], dtype=np.float64
    )
    require_finite("truth effects", result, (len(targets), len(official)))
    axis_audit = {
        **axis_audit,
        "official_genes": len(official),
        "official_measured_genes": int(measured.sum()),
        "official_missing_genes": int((~measured).sum()),
        "official_measured_axis_sha256": sha256_strings(np.asarray(official)[measured]),
    }
    return result, measured, axis_audit


def pds_per_target(
    predicted: np.ndarray, real: np.ndarray, score_mask: np.ndarray
) -> np.ndarray:
    p = np.asarray(predicted[:, score_mask], dtype=np.float64)
    r = np.asarray(real[:, score_mask], dtype=np.float64)
    distances = pairwise_full(p, r, "cosine")
    own = np.diag(distances)[:, None]
    ranks = (distances < own).sum(axis=1)
    ties = (distances == own).sum(axis=1) - 1
    if len(predicted) <= 1:
        return np.ones(len(predicted), dtype=np.float64)
    return 1.0 - (ranks + 0.5 * ties) / (len(predicted) - 1)


def cosine_per_target(
    predicted: np.ndarray, real: np.ndarray, score_mask: np.ndarray
) -> np.ndarray:
    p = np.asarray(predicted[:, score_mask], dtype=np.float64)
    r = np.asarray(real[:, score_mask], dtype=np.float64)
    return np.einsum("ij,ij->i", p, r) / np.maximum(
        np.linalg.norm(p, axis=1) * np.linalg.norm(r, axis=1), 1e-300
    )


def plugin_nmse_proxy_per_target(
    predicted: np.ndarray, real: np.ndarray, score_mask: np.ndarray
) -> np.ndarray:
    p = np.asarray(predicted[:, score_mask], dtype=np.float64)
    r = np.asarray(real[:, score_mask], dtype=np.float64)
    return np.square(p - r).sum(axis=1) / np.maximum(np.square(r).sum(axis=1), 1e-300)


def lfc_from_effect(
    effect: np.ndarray,
    bctrl: np.ndarray,
    ref: ReferenceDE,
) -> np.ndarray:
    profile = (
        np.expm1(
            np.asarray(bctrl, dtype=np.float64)[None, :]
            + np.asarray(effect, dtype=np.float64)
        )
        / TS
    )
    np.clip(profile, 0.0, None, out=profile)
    return np.log2((1e6 * profile[:, ref.gate] + 1e-9) / (ref.ref_lin[None, :] + 1e-9))


def direction_credit_per_target(
    effect: np.ndarray,
    bctrl: np.ndarray,
    targets: Sequence[str],
    ref: ReferenceDE,
) -> np.ndarray:
    if not np.array_equal(ref.genes, np.asarray(ref.genes)):
        raise AssertionError("unreachable")
    pred_lfc = lfc_from_effect(effect, bctrl, ref)
    ref_rows = {target: index for index, target in enumerate(ref.perts)}
    values = np.full(len(targets), np.nan, dtype=np.float64)
    for row, target in enumerate(targets):
        ref_row = ref_rows[str(target)]
        sig = np.asarray(ref.sig[ref_row], dtype=np.int64)
        full_idx = ref.gate_idx[sig]
        truth = ref.lfc[ref_row, sig]
        keep = (ref.genes[full_idx] != str(target)) & np.isfinite(truth) & (truth != 0)
        sig = sig[keep]
        truth = truth[keep]
        if len(sig) == 0:
            continue
        denominator = float(np.abs(truth).sum())
        values[row] = float(np.sign(truth) @ pred_lfc[row, sig]) / max(
            denominator, 1e-300
        )
    return values


def bootstrap(values: np.ndarray, label: str) -> dict[str, float]:
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1 or len(values) == 0 or not np.isfinite(values).all():
        raise ValueError(f"{label}: bootstrap values are empty or non-finite")
    rng = np.random.default_rng(stable_seed("exp106", label))
    draws = np.empty(N_BOOT, dtype=np.float64)
    for index in range(N_BOOT):
        draws[index] = values[rng.integers(0, len(values), len(values))].mean()
    return {
        "estimate": float(values.mean()),
        "lo": float(np.percentile(draws, 2.5)),
        "hi": float(np.percentile(draws, 97.5)),
    }


def common_targets(
    frozen_dirs: Sequence[Path],
    output: Path,
    *,
    expected_axes: dict[str, dict[str, object]] | None = EXPECTED_HELD_TARGET_AXES,
    expected_count: int | None = EXPECTED_COMMON_TARGET_COUNT,
    expected_sha256: str | None = EXPECTED_COMMON_TARGET_SHA256,
) -> None:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    axes = []
    for directory in frozen_dirs:
        manifest = json.loads((directory / "manifest.json").read_text())
        if manifest.get("stage") != "effect-freeze":
            raise ValueError(f"{directory}: not a frozen inference bundle")
        surface = str(manifest.get("surface", ""))
        with np.load(directory / "frozen.npz", allow_pickle=False) as archive:
            targets = tuple(archive["targets"].astype(str))
        if sha256_strings(targets) != manifest.get("target_order_sha256"):
            raise ValueError(f"{directory}: target axis differs from manifest")
        if expected_axes is not None:
            expected = expected_axes.get(surface)
            if expected is None:
                raise ValueError(f"{directory}: unexpected held surface {surface!r}")
            if len(targets) != expected["count"]:
                raise ValueError(
                    f"{directory}: {surface} target count {len(targets)} "
                    f"!= {expected['count']}"
                )
            if sha256_strings(targets) != expected["sha256"]:
                raise ValueError(f"{directory}: {surface} target axis changed")
        axes.append(targets)
    if len(axes) < 2:
        raise ValueError("common-target construction requires at least two surfaces")
    common = set(axes[0])
    for axis in axes[1:]:
        common.intersection_update(axis)
    ordered = tuple(target for target in axes[0] if target in common)
    if not ordered:
        raise ValueError("held surfaces have no common targets")
    if expected_count is not None and expected_sha256 is not None:
        validate_common_target_axis(
            ordered,
            expected_count=expected_count,
            expected_sha256=expected_sha256,
        )
    elif expected_count is not None and len(ordered) != expected_count:
        raise ValueError(
            f"held common gallery has {len(ordered)} targets; expected {expected_count}"
        )
    elif expected_sha256 is not None and sha256_strings(ordered) != expected_sha256:
        raise ValueError("held common target axis changed")
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    tmp.write_text("\n".join(ordered) + "\n")
    tmp.replace(output)
    print(output)


def score(surface: str, frozen_dir: Path, target_allowlist: Path) -> None:
    if surface not in SURFACES:
        raise ValueError(f"unknown surface: {surface}")
    manifest_path = frozen_dir / "manifest.json"
    frozen_path = frozen_dir / "frozen.npz"
    manifest = json.loads(manifest_path.read_text())
    expected_manifest = {
        "experiment": "106-state-baseline",
        "stage": "effect-freeze",
        "surface": surface,
        "sealed_before_perturbation_truth": True,
        "inference_api_has_no_truth_argument": True,
        "expression_rows_read": "non-targeting controls only",
    }
    for key, value in expected_manifest.items():
        if manifest.get(key) != value:
            raise ValueError(
                f"inference manifest {key}={manifest.get(key)!r}, expected {value!r}"
            )
    if manifest.get("sealed_before_perturbation_truth") is not True:
        raise ValueError("inference bundle was not sealed before truth")
    if sha256_file(frozen_path) != manifest.get("frozen_sha256"):
        raise ValueError("frozen inference bundle hash changed")
    checkpoint_manifest_path = verify_manifest_file(
        manifest, "checkpoint_manifest", "checkpoint_manifest_sha256"
    )
    incumbent_npz_path = verify_manifest_file(
        manifest, "incumbent_npz", "incumbent_npz_sha256"
    )
    incumbent_manifest_path = verify_manifest_file(
        manifest, "incumbent_manifest", "incumbent_manifest_sha256"
    )
    held_operator_verified = verify_manifest_file(
        manifest, "held_operator", "held_operator_sha256"
    )
    held_operator_path = (ROOT / str(manifest["held_operator"])).resolve()
    held_operator_metadata_path = held_operator_path.with_suffix(".json")
    if not held_operator_metadata_path.is_file() or sha256_file(
        held_operator_metadata_path
    ) != manifest.get("held_operator_metadata_sha256"):
        raise ValueError("held operator metadata changed after inference")
    controls_path = verify_manifest_file(
        manifest, "controls_h5ad", "controls_h5ad_sha256"
    )
    controls_manifest_path = controls_path.with_suffix(
        controls_path.suffix + ".manifest.json"
    )
    if not controls_manifest_path.is_file() or sha256_file(
        controls_manifest_path
    ) != manifest.get("controls_manifest_sha256"):
        raise ValueError("control-only manifest changed after inference")
    output_path = frozen_dir / "blend-metrics.json"
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite {output_path}")

    with np.load(frozen_path, allow_pickle=False) as z:
        genes = z["genes"].astype(str)
        targets = z["targets"].astype(str)
        bctrl = z["bctrl"].astype(np.float64)
        ctrl_cpm = z["ctrl_cpm"].astype(np.float64)
        incumbent = z["incumbent"].astype(np.float64)
        state = z["state"].astype(np.float64)
        supported = z["state_supported"].astype(bool)
    expected_matrix_shape = (len(targets), len(genes))
    require_finite("bctrl", bctrl, (len(genes),))
    require_finite("ctrl_cpm", ctrl_cpm, (len(genes),))
    require_finite("incumbent", incumbent, expected_matrix_shape)
    require_finite("state", state, expected_matrix_shape)
    if supported.shape != (len(genes),):
        raise ValueError("STATE support mask shape is malformed")
    if sha256_strings(targets) != manifest.get("target_order_sha256"):
        raise ValueError("frozen target axis differs from inference manifest")
    if sha256_strings(genes) != manifest.get("source_axis_sha256"):
        raise ValueError("frozen gene axis differs from inference manifest")

    allowed = read_lines(target_allowlist)
    allowed_set = set(allowed)
    if len(allowed) != len(allowed_set):
        raise ValueError("target allowlist contains duplicates")
    validate_common_target_axis(allowed)
    target_lookup = {target: index for index, target in enumerate(targets)}
    missing = [target for target in allowed if target not in target_lookup]
    if missing:
        raise ValueError(f"target allowlist is absent from {surface}: {missing[:5]}")
    rows = np.asarray([target_lookup[target] for target in allowed], dtype=np.int64)
    targets = targets[rows]
    incumbent = incumbent[rows]
    state = state[rows]
    panel_mask = np.asarray([gene in PANEL_SET for gene in genes], dtype=bool)
    candidates = candidate_effects(
        incumbent, state, bctrl, ctrl_cpm, supported, panel_mask
    )

    truth_hashes_before = truth_bundle_hashes(SURFACES[surface]["truth"])
    refde_hash_before = sha256_file(SURFACES[surface]["refde"])

    # This is the first perturbation-truth read.
    truth, truth_measured, truth_axis_audit = load_truth_effects(
        SURFACES[surface]["truth"], targets, genes
    )
    ref = ReferenceDE.load(SURFACES[surface]["refde"])
    if not np.array_equal(ref.genes, genes):
        raise ValueError("reference-DE and frozen gene axes differ")

    score_mask = ~panel_mask & truth_measured
    if int(score_mask.sum()) != EXPECTED_ORION_SCORE_GENES:
        raise ValueError(
            f"{surface}: measured non-panel score axis has {int(score_mask.sum())} "
            f"genes; expected {EXPECTED_ORION_SCORE_GENES}"
        )
    per_target = {}
    for name, effect in candidates.items():
        require_finite(name, effect, (len(targets), len(genes)))
        per_target[name] = {
            "pds": pds_per_target(effect, truth, score_mask),
            "cosine": cosine_per_target(effect, truth, score_mask),
            "plugin_nmse_proxy": plugin_nmse_proxy_per_target(
                effect, truth, score_mask
            ),
            "credit": direction_credit_per_target(effect, bctrl, targets, ref),
        }
        for metric, values in per_target[name].items():
            require_finite(
                f"{name}:{metric}",
                values,
                (len(targets),),
            )

    baseline = per_target["incumbent"]
    rows = {}
    for name, metrics in per_target.items():
        rows[name] = {
            "means": {
                metric: float(np.mean(values)) for metric, values in metrics.items()
            },
            "delta_vs_incumbent": {
                metric: bootstrap(
                    values - baseline[metric], f"{surface}:{name}:{metric}"
                )
                for metric, values in metrics.items()
            },
        }
    truth_hashes_after = truth_bundle_hashes(SURFACES[surface]["truth"])
    refde_hash_after = sha256_file(SURFACES[surface]["refde"])
    if truth_hashes_after != truth_hashes_before:
        raise RuntimeError("held truth changed while it was being scored")
    if refde_hash_after != refde_hash_before:
        raise RuntimeError("reference DE changed while it was being scored")
    result = {
        "experiment": "106-state-baseline",
        "evaluation_run_id": manifest["evaluation_run_id"],
        "surface": surface,
        "truth_read_after_freeze": True,
        "targets": len(targets),
        "target_order_sha256": sha256_strings(targets),
        "target_allowlist": str(target_allowlist.relative_to(ROOT)),
        "target_allowlist_sha256": sha256_file(target_allowlist),
        "provenance": {
            "inference_manifest_sha256": sha256_file(manifest_path),
            "frozen_inference_sha256": sha256_file(frozen_path),
            "checkpoint_manifest_sha256": sha256_file(checkpoint_manifest_path),
            "checkpoint_sha256": manifest["checkpoint_sha256"],
            "incumbent_npz_sha256": sha256_file(incumbent_npz_path),
            "incumbent_manifest_sha256": sha256_file(incumbent_manifest_path),
            "controls_h5ad_sha256": sha256_file(controls_path),
            "controls_manifest_sha256": sha256_file(controls_manifest_path),
            "held_operator_sha256": sha256_file(held_operator_verified),
            "held_operator_metadata_sha256": sha256_file(held_operator_metadata_path),
            "truth": truth_hashes_before,
            "reference_de_sha256": refde_hash_before,
            "panel_sha256": sha256_file(ROOT / "data/challenge/pert_counts.csv"),
            "scorer_source_sha256": sha256_file(Path(__file__)),
            "cell_eval2_version": importlib.metadata.version("cell-eval2"),
            "pairwise_full_source_sha256": sha256_file(
                Path(inspect.getsourcefile(pairwise_full) or "")
            ),
        },
        "genes": len(genes),
        "truth_axis": truth_axis_audit,
        "score_genes": int(score_mask.sum()),
        "score_gene_axis_sha256": sha256_strings(np.asarray(genes)[score_mask]),
        "state_supported_genes": int(supported.sum()),
        "pds_mask": "all 300 VCC panel coordinates",
        "pds_gallery": (
            "common held-context panel subset; diagnostic only, not the "
            "official 300-target challenge gallery"
        ),
        "direction_mask": "own perturbed gene only",
        "plugin_nmse_proxy": (
            "effect-level squared-error ratio; not official count-level "
            "expr_mse_unbiased_capped_norm"
        ),
        "delta_convention": "candidate minus incumbent",
        "metric_directions": {
            "pds": "higher",
            "cosine": "higher",
            "plugin_nmse_proxy": "lower",
            "credit": "higher",
        },
        "bootstrap_limitations": (
            "target bootstrap is descriptive and does not correct the 20-candidate "
            "screen or recompute gallery-coupled PDS"
        ),
        "min_cpm": MIN_CPM,
        "weights": list(WEIGHTS),
        "candidates": rows,
    }
    atomic_json(output_path, result)
    print(output_path)


def summarize(surface_dirs: Sequence[Path], output: Path) -> None:
    reports = []
    for directory in surface_dirs:
        report = json.loads((directory / "blend-metrics.json").read_text())
        if report["surface"] not in {"hct116", "hek293t"}:
            raise ValueError("primary summary accepts only unseen HCT116/HEK293T")
        reports.append(report)
    if {report["surface"] for report in reports} != {"hct116", "hek293t"}:
        raise ValueError("summary requires exactly HCT116 and HEK293T")
    if len({report["target_order_sha256"] for report in reports}) != 1:
        raise ValueError("held surfaces were not scored on one common target axis")
    shared_fields = {
        "evaluation_run_id": [report.get("evaluation_run_id") for report in reports],
        "target_allowlist_sha256": [
            report.get("target_allowlist_sha256") for report in reports
        ],
        "checkpoint_manifest_sha256": [
            report.get("provenance", {}).get("checkpoint_manifest_sha256")
            for report in reports
        ],
        "checkpoint_sha256": [
            report.get("provenance", {}).get("checkpoint_sha256") for report in reports
        ],
        "panel_sha256": [
            report.get("provenance", {}).get("panel_sha256") for report in reports
        ],
        "scorer_source_sha256": [
            report.get("provenance", {}).get("scorer_source_sha256")
            for report in reports
        ],
        "cell_eval2_version": [
            report.get("provenance", {}).get("cell_eval2_version") for report in reports
        ],
        "pairwise_full_source_sha256": [
            report.get("provenance", {}).get("pairwise_full_source_sha256")
            for report in reports
        ],
    }
    for name, values in shared_fields.items():
        if not values[0] or len(set(values)) != 1:
            raise ValueError(f"held reports do not share one {name}")
    names = set(reports[0]["candidates"])
    if any(set(report["candidates"]) != names for report in reports[1:]):
        raise ValueError("candidate grids differ across surfaces")
    candidates = {}
    for name in sorted(names):
        per_surface = {
            report["surface"]: report["candidates"][name] for report in reports
        }
        macro = {}
        for metric in ("pds", "cosine", "plugin_nmse_proxy", "credit"):
            estimates = [
                row["delta_vs_incumbent"][metric]["estimate"]
                for row in per_surface.values()
            ]
            macro[metric] = float(np.mean(estimates))
        candidates[name] = {
            "surfaces": per_surface,
            "diagnostic_macro_delta": macro,
            "authorization_delta": {
                "pds": macro["pds"],
                "cosine": macro["cosine"],
                "plugin_nmse_proxy": macro["plugin_nmse_proxy"],
            },
            "directional_authorization_delta": {
                "metric": "credit",
                "surface": "hek293t",
                "estimate": float(
                    per_surface["hek293t"]["delta_vs_incumbent"]["credit"]["estimate"]
                ),
                "pending": "H1 confirmation",
            },
        }
    result = {
        "experiment": "106-state-baseline",
        "selection_status": (
            "held-context diagnostic only; candidate ordering is not a "
            "production blend-weight selection and authorizes no submission"
        ),
        "multiplicity_adjusted": False,
        "bootstrap_limitations": (
            "descriptive only: no candidate multiplicity correction and no "
            "gallery-recomputed PDS bootstrap"
        ),
        "official_count_level_metrics": False,
        "directional_authorization": (
            "HEK293T-only pending H1; HCT116 credit is diagnostic-only and excluded"
        ),
        "hct116_directional_credit": "diagnostic-only",
        "surfaces": ["hct116", "hek293t"],
        "shared_provenance": {key: values[0] for key, values in shared_fields.items()},
        "surface_provenance": {
            report["surface"]: report["provenance"] for report in reports
        },
        "candidates": candidates,
    }
    atomic_json(output, result)
    print(output)


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)
    score_p = sub.add_parser("score")
    score_p.add_argument("--surface", choices=tuple(SURFACES), required=True)
    score_p.add_argument("--frozen-dir", type=Path, required=True)
    score_p.add_argument("--target-allowlist", type=Path, required=True)
    common_p = sub.add_parser("common-targets")
    common_p.add_argument("--frozen-dir", type=Path, action="append", required=True)
    common_p.add_argument("--output", type=Path, required=True)
    summary_p = sub.add_parser("summarize")
    summary_p.add_argument("--surface-dir", type=Path, action="append", required=True)
    summary_p.add_argument("--output", type=Path, required=True)
    return p


def main() -> None:
    args = parser().parse_args()
    if args.command == "score":
        score(
            args.surface,
            args.frozen_dir.resolve(),
            args.target_allowlist.resolve(),
        )
    elif args.command == "common-targets":
        common_targets(
            [path.resolve() for path in args.frozen_dir],
            args.output.resolve(),
        )
    else:
        summarize(
            [path.resolve() for path in args.surface_dir],
            args.output.resolve(),
        )


if __name__ == "__main__":
    main()

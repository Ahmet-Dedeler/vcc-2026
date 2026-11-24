#!/usr/bin/env python3
"""Run Experiment 125's strictly nested source-only evaluation."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import multiprocessing
import os
import shutil
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from scipy.special import ndtr

import freeze_plans


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
AUTHORIZATION_PATH = HERE / "SOURCE_ONLY_AUTHORIZATION.json"
IMPLEMENTATION_FREEZE_PATH = HERE / "IMPLEMENTATION_FREEZE.json"
PLAN_DIR = HERE / "work" / "preregistered-freeze-v1"
DEFAULT_OUTPUT = HERE / "work" / "source-only-evaluation-v1"
TS = 50000.0
FAMILY_NAMES = ("replogle_k562", "orion", "cd4")
TRANSPORT_NAMES = ("production_log1p", "abundance_exponential")


class ExperimentError(RuntimeError):
    """A scientific, nesting, source, or immutable-output invariant failed."""


@dataclass(frozen=True)
class PosteriorSetting:
    name: str
    pseudo_count: float
    variance_floor_fraction: float
    prior_scale: float


@dataclass(frozen=True)
class Candidate:
    name: str
    posterior: PosteriorSetting
    transport: str
    order: int


@dataclass(frozen=True)
class SourceData:
    panel: tuple[str, ...]
    genes: tuple[str, ...]
    family_names: tuple[str, ...]
    effects: np.ndarray
    valid: np.ndarray
    baseline_q: np.ndarray
    baseline_valid: np.ndarray


@dataclass(frozen=True)
class PreparedFamilies:
    effects: np.ndarray
    valid: np.ndarray
    normalized: np.ndarray
    row_scales: np.ndarray
    row_origin: np.ndarray


@dataclass(frozen=True)
class RawMoments:
    pair_count: np.ndarray
    prior_mean_raw: np.ndarray
    latent_raw: np.ndarray
    noise0_raw: np.ndarray
    noise1_raw: np.ndarray
    global_latent: float
    global_noise0: float
    global_noise1: float
    global_total: float
    fit_destination_rows: int
    fit_biological_rows: int


@dataclass(frozen=True)
class PosteriorParameters:
    prior_mean: np.ndarray
    latent_variance: np.ndarray
    noise0_variance: np.ndarray
    noise1_variance: np.ndarray
    floor: float


@dataclass(frozen=True)
class MetricResult:
    credit: np.ndarray
    precision: np.ndarray
    evaluable: np.ndarray


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ExperimentError(message)


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_array(array: np.ndarray) -> str:
    value = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(str(value.dtype).encode())
    digest.update(b"\0")
    digest.update(json.dumps(list(value.shape), separators=(",", ":")).encode())
    digest.update(b"\0")
    digest.update(value.tobytes(order="C"))
    return digest.hexdigest()


def json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    require(isinstance(value, dict), f"JSON root must be an object: {path}")
    return value


def read_axis(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(path.read_text(encoding="utf-8").splitlines())
    require(bool(values), f"empty axis: {path}")
    require(all(value and value == value.strip() for value in values), f"bad axis: {path}")
    if unique:
        require(len(values) == len(set(values)), f"duplicate axis values: {path}")
    return values


def read_csv_axis(path: Path, column: str) -> tuple[str, ...]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require(reader.fieldnames == [column], f"header changed: {path}")
        values = tuple(str(row[column]) for row in reader)
    require(len(values) == len(set(values)), f"duplicate values: {path}")
    return values


def candidates_from_protocol(protocol: dict[str, Any]) -> tuple[Candidate, ...]:
    settings = tuple(
        PosteriorSetting(
            name=str(item["name"]),
            pseudo_count=float(item["pseudo_count"]),
            variance_floor_fraction=float(item["variance_floor_fraction"]),
            prior_scale=float(item["prior_scale"]),
        )
        for item in protocol["candidate_grid"]
    )
    transports = tuple(
        str(item["name"]) for item in protocol["transport_semantics"]["arms"]
    )
    require(transports == TRANSPORT_NAMES, "transport order changed")
    candidates: list[Candidate] = []
    for setting in settings:
        for transport in transports:
            candidates.append(
                Candidate(
                    name=f"{setting.name}__{transport}",
                    posterior=setting,
                    transport=transport,
                    order=len(candidates),
                )
            )
    require(len(candidates) == 8, "candidate cross-product changed")
    return tuple(candidates)


def validate_implementation_freeze() -> dict[str, Any]:
    require(IMPLEMENTATION_FREEZE_PATH.is_file(), "implementation freeze missing")
    freeze = read_json(IMPLEMENTATION_FREEZE_PATH)
    if freeze.get("created_before_source_effect_loading") is not True:
        require(
            freeze.get("scientific_spec_frozen_before_source_effect_loading") is True,
            "implementation freeze timing changed",
        )
        recovery = freeze.get("implementation_recovery")
        require(isinstance(recovery, dict), "implementation recovery missing")
        recovery_path = ROOT / recovery["path"]
        require(
            recovery_path.is_file()
            and sha256_file(recovery_path) == recovery["sha256"],
            "implementation recovery hash mismatch",
        )
    require(
        freeze["plan_freeze_manifest_sha256"]
        == sha256_file(PLAN_DIR / "freeze-manifest.json"),
        "plan freeze binding changed",
    )
    for relative, digest in freeze["files"].items():
        path = ROOT / relative
        require(path.is_file() and not path.is_symlink(), f"missing code file: {relative}")
        require(sha256_file(path) == digest, f"implementation changed: {relative}")
    return freeze


def validate_authorization(*, hash_effect_arrays: bool) -> dict[str, Any]:
    freeze_plans.verify(PLAN_DIR)
    validate_implementation_freeze()
    authorization = read_json(AUTHORIZATION_PATH)
    require(authorization["source_only_authorized"] is True, "source-only not authorized")
    require(
        authorization["protocol_sha256"] == sha256_file(PROTOCOL_PATH),
        "protocol authorization hash mismatch",
    )
    denied = (
        "read_pan",
        "read_dld1",
        "read_h1",
        "read_orion_held_harness_outcomes",
        "read_challenge_truth",
        "read_challenge_context_expression_values",
        "run_cell_eval2",
        "generate_counts",
        "build_submission",
        "submit",
        "commit",
        "push",
    )
    for key in denied:
        require(authorization["authorization"].get(key) is False, f"denial changed: {key}")
    effect_paths = {
        "experiments/122-panel-response-atlas/work/source-only-phase1-v1/family-effects.npy",
        "experiments/122-panel-response-atlas/work/source-only-phase1-v1/family-valid.npy",
    }
    for relative, digest in authorization["allowed_input_files"].items():
        path = (ROOT / relative).resolve()
        require(path.is_relative_to(ROOT), f"input escapes repository: {relative}")
        require(path.is_file() and not path.is_symlink(), f"missing input: {relative}")
        if hash_effect_arrays or relative not in effect_paths:
            require(sha256_file(path) == digest, f"input hash changed: {relative}")
    return authorization


def panel_positions(axis: Sequence[str], panel: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    positions: dict[str, list[int]] = {}
    for index, label in enumerate(axis):
        positions.setdefault(label, []).append(index)
    duplicated = [label for label in panel if len(positions.get(label, [])) > 1]
    require(not duplicated, f"duplicated panel labels: {duplicated[:5]}")
    destination: list[int] = []
    source: list[int] = []
    for panel_index, label in enumerate(panel):
        if label in positions:
            destination.append(panel_index)
            source.append(positions[label][0])
    return np.asarray(destination, dtype=np.int64), np.asarray(source, dtype=np.int64)


def member_baseline(root: Path, panel: Sequence[str]) -> tuple[np.ndarray, np.ndarray]:
    genes = read_axis(root / "genes.txt", unique=False)
    destination, source = panel_positions(genes, panel)
    b_ctrl = np.load(root / "b_ctrl.npy", mmap_mode="r", allow_pickle=False)
    require(b_ctrl.shape == (len(genes),), f"baseline shape changed: {root}")
    q = np.zeros(len(panel), dtype=np.float64)
    valid = np.zeros(len(panel), dtype=bool)
    values = np.expm1(np.asarray(b_ctrl[source], dtype=np.float64)) / TS
    good = np.isfinite(values) & (values > 0.0)
    q[destination[good]] = values[good]
    valid[destination[good]] = True
    return q, valid


def collapse_baselines(
    members: Sequence[tuple[np.ndarray, np.ndarray]],
) -> tuple[np.ndarray, np.ndarray]:
    values = np.stack([item[0] for item in members])
    valid = np.stack([item[1] for item in members])
    masked = np.where(valid, values, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        collapsed = np.nanmedian(masked, axis=0)
    family_valid = np.any(valid, axis=0) & np.isfinite(collapsed) & (collapsed > 0.0)
    collapsed = np.where(family_valid, collapsed, 0.0)
    return collapsed, family_valid


def load_source_data() -> SourceData:
    validate_authorization(hash_effect_arrays=True)
    protocol = read_json(PROTOCOL_PATH)
    panel = read_axis(ROOT / protocol["inputs"]["panel_targets"])
    canonical_panel = read_csv_axis(ROOT / "data/challenge/pert_counts.csv", "target_gene")
    require(panel == canonical_panel, "panel axis mismatch")
    genes = read_csv_axis(ROOT / protocol["inputs"]["canonical_genes"], "gene_name")
    require(len(genes) == 18533, "canonical gene count changed")
    family_names = read_axis(ROOT / protocol["inputs"]["family_names"])
    require(family_names == FAMILY_NAMES, "family order changed")
    effects = np.load(ROOT / protocol["inputs"]["family_effects"], allow_pickle=False)
    valid = np.load(ROOT / protocol["inputs"]["family_valid"], allow_pickle=False)
    require(effects.shape == (3, 300, 300), "family effects shape changed")
    require(valid.shape == effects.shape and valid.dtype == np.bool_, "family validity changed")
    effects = np.asarray(effects, dtype=np.float64)
    valid = np.asarray(valid, dtype=bool)
    require(np.all(np.isfinite(effects)), "nonfinite source effect")
    require(np.all(effects[~valid] == 0.0), "invalid source effect is nonzero")
    diagonal = np.arange(300)
    require(np.all(effects[:, diagonal, diagonal] == 0.0), "source self effect changed")
    require(np.all(~valid[:, diagonal, diagonal]), "source self validity changed")

    k562 = collapse_baselines(
        (member_baseline(ROOT / "data/donors/k562_gwps", panel),)
    )
    orion = collapse_baselines(
        (
            member_baseline(ROOT / "data/donors/orion_hct116", panel),
            member_baseline(ROOT / "data/donors/orion_hek293t", panel),
        )
    )
    cd4_root = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
    cd4 = collapse_baselines(
        (
            member_baseline(cd4_root / "cd4_rest", panel),
            member_baseline(cd4_root / "cd4_stim8hr", panel),
            member_baseline(cd4_root / "cd4_stim48hr", panel),
        )
    )
    baseline_q = np.stack((k562[0], orion[0], cd4[0]))
    baseline_valid = np.stack((k562[1], orion[1], cd4[1]))
    require(np.all(np.isfinite(baseline_q)), "nonfinite family baseline")
    require(np.all(baseline_q >= 0.0), "negative family baseline")
    return SourceData(
        panel=panel,
        genes=genes,
        family_names=family_names,
        effects=effects,
        valid=valid,
        baseline_q=baseline_q,
        baseline_valid=baseline_valid,
    )


def transform_training(
    effects: np.ndarray,
    valid: np.ndarray,
    *,
    target_permutation: np.ndarray | None,
    response_permutation: np.ndarray | None,
) -> PreparedFamilies:
    n = effects.shape[1]
    rows = np.arange(n) if target_permutation is None else target_permutation
    columns = np.arange(n) if response_permutation is None else response_permutation
    transformed_effects = np.asarray(effects[:, rows][:, :, columns], dtype=np.float64).copy()
    transformed_valid = np.asarray(valid[:, rows][:, :, columns], dtype=bool).copy()
    diagonal = np.arange(n)
    transformed_effects[:, diagonal, diagonal] = 0.0
    transformed_valid[:, diagonal, diagonal] = False
    scales = np.zeros((effects.shape[0], n), dtype=np.float64)
    normalized = np.zeros_like(transformed_effects)
    for family in range(effects.shape[0]):
        for target in range(n):
            values = np.abs(
                transformed_effects[family, target, transformed_valid[family, target]]
            )
            values = values[np.isfinite(values) & (values > 0.0)]
            if values.size:
                scales[family, target] = float(np.median(values))
        usable = transformed_valid[family] & (scales[family, :, None] > 0.0)
        np.divide(
            transformed_effects[family],
            scales[family, :, None],
            out=normalized[family],
            where=usable,
        )
        transformed_valid[family] = usable
    return PreparedFamilies(
        effects=transformed_effects,
        valid=transformed_valid,
        normalized=normalized,
        row_scales=scales,
        row_origin=np.asarray(rows, dtype=np.int64),
    )


def fit_raw_moments(
    prepared: PreparedFamilies,
    fit_rows: np.ndarray,
    forbidden_biological_targets: np.ndarray,
) -> RawMoments:
    require(prepared.effects.shape[0] == 2, "moments require exactly two families")
    fit_rows = np.asarray(fit_rows, dtype=bool)
    forbidden = np.asarray(forbidden_biological_targets, dtype=bool)
    require(fit_rows.shape == (300,) and forbidden.shape == (300,), "fit mask changed")
    biological_allowed = ~forbidden[prepared.row_origin]
    rows = fit_rows & biological_allowed
    pair = prepared.valid[0] & prepared.valid[1] & rows[:, None]
    n = pair.sum(axis=0).astype(np.float64)
    x = prepared.normalized[0]
    y = prepared.normalized[1]
    sum_x = np.where(pair, x, 0.0).sum(axis=0)
    sum_y = np.where(pair, y, 0.0).sum(axis=0)
    mean_x = np.divide(sum_x, n, out=np.zeros(300), where=n > 0)
    mean_y = np.divide(sum_y, n, out=np.zeros(300), where=n > 0)
    prior_mean = 0.5 * (mean_x + mean_y)
    df = np.maximum(n - 1.0, 0.0)
    ss_x = np.where(pair, x * x, 0.0).sum(axis=0) - np.divide(
        sum_x * sum_x,
        n,
        out=np.zeros(300),
        where=n > 0,
    )
    ss_y = np.where(pair, y * y, 0.0).sum(axis=0) - np.divide(
        sum_y * sum_y,
        n,
        out=np.zeros(300),
        where=n > 0,
    )
    ss_xy = np.where(pair, x * y, 0.0).sum(axis=0) - np.divide(
        sum_x * sum_y,
        n,
        out=np.zeros(300),
        where=n > 0,
    )
    var_x = np.divide(ss_x, df, out=np.zeros(300), where=df > 0)
    var_y = np.divide(ss_y, df, out=np.zeros(300), where=df > 0)
    covariance = np.divide(ss_xy, df, out=np.zeros(300), where=df > 0)
    stable = df >= 2.0
    latent_raw = np.where(stable, np.maximum(covariance, 0.0), 0.0)
    noise0_raw = np.where(stable, np.maximum(var_x - covariance, 0.0), 0.0)
    noise1_raw = np.where(stable, np.maximum(var_y - covariance, 0.0), 0.0)
    weights = np.where(stable, df, 0.0)
    weight_sum = float(weights.sum())
    require(weight_sum > 0.0, "no paired degrees of freedom for nested fit")

    def weighted_global(values: np.ndarray) -> float:
        result = float(np.sum(weights * values) / weight_sum)
        return max(result, 0.0)

    global_latent = weighted_global(latent_raw)
    global_noise0 = weighted_global(noise0_raw)
    global_noise1 = weighted_global(noise1_raw)
    global_total = max(
        global_latent + 0.5 * (global_noise0 + global_noise1),
        1e-8,
    )
    return RawMoments(
        pair_count=n,
        prior_mean_raw=prior_mean,
        latent_raw=latent_raw,
        noise0_raw=noise0_raw,
        noise1_raw=noise1_raw,
        global_latent=global_latent,
        global_noise0=global_noise0,
        global_noise1=global_noise1,
        global_total=global_total,
        fit_destination_rows=int(np.count_nonzero(fit_rows)),
        fit_biological_rows=int(np.count_nonzero(rows)),
    )


def posterior_parameters(
    moments: RawMoments,
    setting: PosteriorSetting,
) -> PosteriorParameters:
    count = moments.pair_count
    pseudo = setting.pseudo_count
    mean_weight = np.divide(count, count + pseudo)
    component_weight = np.divide(np.maximum(count - 1.0, 0.0), np.maximum(count - 1.0, 0.0) + pseudo)
    prior_mean = mean_weight * moments.prior_mean_raw
    latent = (
        component_weight * moments.latent_raw
        + (1.0 - component_weight) * moments.global_latent
    ) * setting.prior_scale
    noise0 = (
        component_weight * moments.noise0_raw
        + (1.0 - component_weight) * moments.global_noise0
    )
    noise1 = (
        component_weight * moments.noise1_raw
        + (1.0 - component_weight) * moments.global_noise1
    )
    floor = setting.variance_floor_fraction * moments.global_total
    floor = max(float(floor), 1e-10)
    latent = np.maximum(latent, floor)
    noise0 = np.maximum(noise0, floor)
    noise1 = np.maximum(noise1, floor)
    require(np.all(np.isfinite(prior_mean)), "nonfinite prior mean")
    require(np.all(latent > 0.0), "nonpositive latent variance")
    require(np.all(noise0 > 0.0) and np.all(noise1 > 0.0), "nonpositive noise")
    return PosteriorParameters(prior_mean, latent, noise0, noise1, floor)


def output_scale(
    prepared: PreparedFamilies,
    observed: Sequence[int],
) -> tuple[np.ndarray, np.ndarray]:
    scales = prepared.row_scales[np.asarray(observed), :, None]
    valid = prepared.valid[np.asarray(observed)]
    scale_values = np.broadcast_to(scales, valid.shape)
    masked = np.where(valid, scale_values, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        result = np.nanmedian(masked, axis=0)
    result = np.where(np.isfinite(result), result, 0.0)
    support = valid.sum(axis=0).astype(np.int16)
    return result, support


def posterior_prediction(
    prepared: PreparedFamilies,
    observed: Sequence[int],
    parameters: PosteriorParameters,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    observed = tuple(observed)
    require(bool(observed), "posterior has no observations")
    prior_precision = 1.0 / parameters.latent_variance
    precision = np.broadcast_to(prior_precision, (300, 300)).copy()
    numerator = np.broadcast_to(
        parameters.prior_mean / parameters.latent_variance,
        (300, 300),
    ).copy()
    for local_index in observed:
        variance = (
            parameters.noise0_variance
            if local_index == 0
            else parameters.noise1_variance
        )
        valid = prepared.valid[local_index]
        weight = 1.0 / variance
        precision += valid * weight[None, :]
        numerator += valid * prepared.normalized[local_index] * weight[None, :]
    posterior_variance = 1.0 / precision
    posterior_mean = posterior_variance * numerator
    scale, support = output_scale(prepared, observed)
    allowed = support > 0
    effect = np.where(allowed, posterior_mean * scale, 0.0)
    z = np.divide(
        posterior_mean,
        np.sqrt(posterior_variance),
        out=np.zeros_like(posterior_mean),
        where=posterior_variance > 0.0,
    )
    confidence = np.where(allowed, np.abs(2.0 * ndtr(z) - 1.0), 0.0)
    diagonal = np.arange(300)
    effect[diagonal, diagonal] = 0.0
    confidence[diagonal, diagonal] = 0.0
    support[diagonal, diagonal] = 0
    require(np.all(np.isfinite(effect)), "nonfinite posterior effect")
    require(np.all(np.isfinite(confidence)), "nonfinite posterior confidence")
    require(np.all((confidence >= 0.0) & (confidence <= 1.0)), "bad confidence")
    return effect, confidence, support


def equal_mean_agreement(
    prepared: PreparedFamilies,
    observed: Sequence[int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    observed = np.asarray(tuple(observed), dtype=np.int64)
    valid = prepared.valid[observed]
    values = prepared.normalized[observed]
    support = valid.sum(axis=0).astype(np.int16)
    location = np.divide(
        np.where(valid, values, 0.0).sum(axis=0),
        support,
        out=np.zeros((300, 300), dtype=np.float64),
        where=support > 0,
    )
    nonzero = valid & (values != 0.0)
    nonzero_count = nonzero.sum(axis=0)
    sign_sum = np.where(nonzero, np.sign(values), 0.0).sum(axis=0)
    agreement = np.divide(
        np.abs(sign_sum),
        nonzero_count,
        out=np.zeros((300, 300), dtype=np.float64),
        where=nonzero_count > 0,
    )
    scale, _ = output_scale(prepared, observed)
    effect = np.where(support > 0, location * scale, 0.0)
    confidence = (
        np.abs(location)
        * (support.astype(np.float64) / len(observed))
        * agreement
    )
    diagonal = np.arange(300)
    effect[diagonal, diagonal] = 0.0
    confidence[diagonal, diagonal] = 0.0
    support[diagonal, diagonal] = 0
    return effect, confidence, support


def repair_to_sum(values: np.ndarray, target_sum: float, mutable: np.ndarray) -> None:
    indices = np.flatnonzero(mutable)
    require(indices.size > 0, "no mutable coordinates for repair")
    order = indices[np.argsort(values[indices], kind="stable")]
    for index_value in order:
        index = int(index_value)
        original = float(values[index])
        other_sum = math.fsum(
            float(values[other])
            for other in range(values.size)
            if other != index
        )
        candidate = target_sum - other_sum
        if candidate < 0.0 or not math.isfinite(candidate):
            continue
        values[index] = candidate
        observed = math.fsum(map(float, values))
        if observed == target_sum:
            require(np.all(values >= 0.0), "sum repair made negative abundance")
            return
        direction = math.inf if observed < target_sum else -math.inf
        current = float(values[index])
        for _ in range(128):
            current = math.nextafter(current, direction)
            if current < 0.0:
                break
            values[index] = current
            observed = math.fsum(map(float, values))
            if observed == target_sum:
                require(np.all(values >= 0.0), "sum repair made negative abundance")
                return
            crossed = (
                (direction == math.inf and observed > target_sum)
                or (direction == -math.inf and observed < target_sum)
            )
            if crossed:
                break
        values[index] = original
    raise ExperimentError("could not enforce exact compensated sum")


def transport_row(
    baseline_q: np.ndarray,
    baseline_valid: np.ndarray,
    residual: np.ndarray,
    support: np.ndarray,
    target: int,
    arm: str,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    baseline = np.asarray(baseline_q, dtype=np.float64)
    valid_baseline = np.asarray(baseline_valid, dtype=bool)
    candidate = baseline.copy()
    mutable = (
        valid_baseline
        & (baseline > 0.0)
        & np.isfinite(residual)
        & (support > 0)
    )
    mutable[target] = False
    preclip_negative = 0
    if np.any(mutable):
        original = baseline[mutable]
        target_sum = math.fsum(map(float, original))
        r = np.asarray(residual[mutable], dtype=np.float64)
        if arm == "production_log1p":
            transported = np.expm1(np.log1p(TS * original) + r) / TS
            preclip_negative = int(np.count_nonzero(transported < 0.0))
            transported = np.maximum(transported, 0.0)
        elif arm == "abundance_exponential":
            logits = r - float(np.max(r))
            transported = original * np.exp(logits)
        else:
            raise ExperimentError(f"unknown transport: {arm}")
        transported_sum = math.fsum(map(float, transported))
        require(transported_sum > 0.0, "transport collapsed mutable mass")
        transported *= target_sum / transported_sum
        repair_to_sum(transported, target_sum, np.ones(transported.size, dtype=bool))
        candidate[mutable] = transported
    original_mutable_sum = math.fsum(map(float, baseline[mutable]))
    candidate_mutable_sum = math.fsum(map(float, candidate[mutable]))
    original_immutable_sum = math.fsum(map(float, baseline[~mutable]))
    candidate_immutable_sum = math.fsum(map(float, candidate[~mutable]))
    original_panel_sum = math.fsum((original_mutable_sum, original_immutable_sum))
    candidate_panel_sum = math.fsum((candidate_mutable_sum, candidate_immutable_sum))
    realized = np.zeros(300, dtype=np.float64)
    positive = valid_baseline & (baseline > 0.0)
    realized[positive] = (
        np.log1p(TS * candidate[positive]) - np.log1p(TS * baseline[positive])
    )
    realized[target] = 0.0
    realized_support = mutable.astype(np.int16)
    invariants = {
        "finite": bool(np.all(np.isfinite(candidate)) and np.all(np.isfinite(realized))),
        "nonnegative": bool(np.all(candidate >= 0.0)),
        "target_self_abundance_identity": bool(candidate[target] == baseline[target]),
        "immutable_panel_identity": bool(np.array_equal(candidate[~mutable], baseline[~mutable])),
        "mutable_mass_conservation": bool(
            not np.any(mutable)
            or candidate_mutable_sum == original_mutable_sum
        ),
        "panel_mass_conservation": bool(candidate_panel_sum == original_panel_sum),
        "target_self_effect_zero": bool(realized[target] == 0.0),
        "panel_only_write": True,
    }
    require(all(invariants.values()), f"{arm}: transport invariant failed")
    return realized, realized_support, {
        "preclip_negative": preclip_negative,
        "mutable": int(np.count_nonzero(mutable)),
        "invariants": invariants,
    }


def realize_matrix(
    latent_effect: np.ndarray,
    confidence: np.ndarray,
    support: np.ndarray,
    baseline_q: np.ndarray,
    baseline_valid: np.ndarray,
    arm: str,
    target_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    effect = np.zeros_like(latent_effect)
    output_confidence = np.zeros_like(confidence)
    output_support = np.zeros_like(support)
    preclip_negative = 0
    mutable_count = 0
    invariant_rows = 0
    for target in np.flatnonzero(target_mask):
        row, row_support, audit = transport_row(
            baseline_q,
            baseline_valid,
            latent_effect[target],
            support[target],
            int(target),
            arm,
        )
        effect[target] = row
        output_support[target] = row_support
        output_confidence[target] = np.where(row_support > 0, confidence[target], 0.0)
        preclip_negative += int(audit["preclip_negative"])
        mutable_count += int(audit["mutable"])
        invariant_rows += int(all(audit["invariants"].values()))
    return effect, output_confidence, output_support, {
        "rows": int(np.count_nonzero(target_mask)),
        "invariant_rows": invariant_rows,
        "preclip_negative_coordinates": preclip_negative,
        "mutable_coordinates": mutable_count,
    }


def score_predictions(
    prediction: np.ndarray,
    confidence: np.ndarray,
    support: np.ndarray,
    held_effect: np.ndarray,
    held_valid: np.ndarray,
    target_mask: np.ndarray,
) -> MetricResult:
    credit = np.full(300, np.nan, dtype=np.float64)
    precision = np.full(300, np.nan, dtype=np.float64)
    evaluable = np.zeros(300, dtype=bool)
    canonical = np.arange(300)
    for target in np.flatnonzero(target_mask):
        valid = held_valid[target].copy()
        valid[target] = False
        indices = np.flatnonzero(valid)
        if indices.size < 20 or not np.any(support[target] > 0):
            continue
        truth = held_effect[target, indices]
        denominator = float(np.sum(np.abs(truth)))
        if not np.isfinite(denominator) or denominator <= 0.0:
            continue
        credit[target] = float(
            np.sum(np.sign(truth) * prediction[target, indices]) / denominator
        )
        top_n = max(1, int(math.ceil(0.20 * indices.size)))
        order = np.lexsort((canonical[indices], -confidence[target, indices]))[:top_n]
        chosen = indices[order]
        weights = np.abs(held_effect[target, chosen])
        weight_sum = float(np.sum(weights))
        if not np.isfinite(weight_sum) or weight_sum <= 0.0:
            credit[target] = np.nan
            continue
        predicted_sign = np.sign(prediction[target, chosen])
        truth_sign = np.sign(held_effect[target, chosen])
        correct = (
            (predicted_sign != 0.0)
            & (truth_sign != 0.0)
            & (predicted_sign == truth_sign)
        )
        precision[target] = float(np.sum(weights * correct) / weight_sum)
        evaluable[target] = True
    return MetricResult(credit, precision, evaluable)


def metric_delta(
    candidate: MetricResult,
    baseline: MetricResult,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    evaluable = candidate.evaluable & baseline.evaluable
    credit = np.full(300, np.nan, dtype=np.float64)
    precision = np.full(300, np.nan, dtype=np.float64)
    credit[evaluable] = candidate.credit[evaluable] - baseline.credit[evaluable]
    precision[evaluable] = (
        candidate.precision[evaluable] - baseline.precision[evaluable]
    )
    return credit, precision, evaluable


def select_candidate(
    summaries: Sequence[dict[str, Any]],
    candidates: Sequence[Candidate],
) -> Candidate:
    by_name = {candidate.name: candidate for candidate in candidates}
    require(set(by_name) == {str(item["name"]) for item in summaries}, "selection grid drift")

    def key(item: dict[str, Any]) -> tuple[float, ...]:
        eligible = float(
            min(item["family_credit_delta"]) >= 0.0
            and min(item["family_precision_delta"]) >= 0.0
        )
        return (
            eligible,
            round(min(item["family_credit_delta"]), 12),
            round(min(item["family_precision_delta"]), 12),
            round(float(item["macro_credit_delta"]), 12),
            round(float(item["macro_precision_delta"]), 12),
            -float(by_name[str(item["name"])].order),
        )

    selected_summary = max(summaries, key=key)
    return by_name[str(selected_summary["name"])]


def inner_select(
    prepared: PreparedFamilies,
    training_global_indices: tuple[int, int],
    family_baselines: np.ndarray,
    family_baseline_valid: np.ndarray,
    folds: np.ndarray,
    outer_test_fold: int,
    candidates: Sequence[Candidate],
) -> tuple[Candidate, list[dict[str, Any]]]:
    records: dict[str, dict[int, list[np.ndarray]]] = {
        candidate.name: {
            0: [[], []],
            1: [[], []],
        }
        for candidate in candidates
    }
    for inner_fold in range(5):
        if inner_fold == outer_test_fold:
            continue
        validation = folds == inner_fold
        forbidden = (folds == outer_test_fold) | validation
        fit_rows = ~forbidden
        moments = fit_raw_moments(prepared, fit_rows, forbidden)
        baseline_cache: dict[tuple[str, int], MetricResult] = {}
        for candidate in candidates:
            parameters = posterior_parameters(moments, candidate.posterior)
            for pseudoheld in (0, 1):
                observed = (1 - pseudoheld,)
                baseline_key = (candidate.transport, pseudoheld)
                if baseline_key not in baseline_cache:
                    base_latent, base_conf, base_support = equal_mean_agreement(
                        prepared,
                        observed,
                    )
                    base_effect, base_conf, base_support, _ = realize_matrix(
                        base_latent,
                        base_conf,
                        base_support,
                        family_baselines[training_global_indices[pseudoheld]],
                        family_baseline_valid[training_global_indices[pseudoheld]],
                        candidate.transport,
                        validation,
                    )
                    baseline_cache[baseline_key] = score_predictions(
                        base_effect,
                        base_conf,
                        base_support,
                        prepared.effects[pseudoheld],
                        prepared.valid[pseudoheld],
                        validation,
                    )
                cand_latent, cand_conf, cand_support = posterior_prediction(
                    prepared,
                    observed,
                    parameters,
                )
                cand_effect, cand_conf, cand_support, _ = realize_matrix(
                    cand_latent,
                    cand_conf,
                    cand_support,
                    family_baselines[training_global_indices[pseudoheld]],
                    family_baseline_valid[training_global_indices[pseudoheld]],
                    candidate.transport,
                    validation,
                )
                candidate_metric = score_predictions(
                    cand_effect,
                    cand_conf,
                    cand_support,
                    prepared.effects[pseudoheld],
                    prepared.valid[pseudoheld],
                    validation,
                )
                delta_credit, delta_precision, evaluable = metric_delta(
                    candidate_metric,
                    baseline_cache[baseline_key],
                )
                records[candidate.name][pseudoheld][0].append(
                    delta_credit[evaluable]
                )
                records[candidate.name][pseudoheld][1].append(
                    delta_precision[evaluable]
                )
    summaries: list[dict[str, Any]] = []
    for candidate in candidates:
        family_credit: list[float] = []
        family_precision: list[float] = []
        family_counts: list[int] = []
        for pseudoheld in (0, 1):
            credits = np.concatenate(records[candidate.name][pseudoheld][0])
            precisions = np.concatenate(records[candidate.name][pseudoheld][1])
            require(credits.size > 0 and precisions.size == credits.size, "empty inner fold")
            family_credit.append(float(np.mean(credits)))
            family_precision.append(float(np.mean(precisions)))
            family_counts.append(int(credits.size))
        summaries.append(
            {
                "name": candidate.name,
                "family_credit_delta": family_credit,
                "family_precision_delta": family_precision,
                "family_target_counts": family_counts,
                "macro_credit_delta": float(np.mean(family_credit)),
                "macro_precision_delta": float(np.mean(family_precision)),
            }
        )
    return select_candidate(summaries, candidates), summaries


def nested_evaluate(
    source: SourceData,
    folds: np.ndarray,
    candidates: Sequence[Candidate],
    *,
    target_permutation: np.ndarray | None = None,
    response_permutation: np.ndarray | None = None,
    detailed: bool,
) -> dict[str, Any]:
    family_credit_delta = np.full(3, np.nan, dtype=np.float64)
    family_precision_delta = np.full(3, np.nan, dtype=np.float64)
    per_target_credit_delta = np.full((3, 300), np.nan, dtype=np.float64)
    per_target_precision_delta = np.full((3, 300), np.nan, dtype=np.float64)
    candidate_effect = np.zeros((3, 300, 300), dtype=np.float64)
    candidate_confidence = np.zeros_like(candidate_effect)
    baseline_effect = np.zeros_like(candidate_effect)
    baseline_confidence = np.zeros_like(candidate_effect)
    selected_indices = np.full((3, 5), -1, dtype=np.int16)
    fold_records: list[dict[str, Any]] = []
    total_candidate_support = 0
    total_baseline_support = 0
    transport_invariant_rows = 0
    transport_total_rows = 0
    preclip_negative = 0

    for held_index in range(3):
        training_indices = tuple(index for index in range(3) if index != held_index)
        transformed = transform_training(
            source.effects[np.asarray(training_indices)],
            source.valid[np.asarray(training_indices)],
            target_permutation=target_permutation,
            response_permutation=response_permutation,
        )
        for outer_fold in range(5):
            test_mask = folds == outer_fold
            selected, inner_summaries = inner_select(
                transformed,
                training_indices,
                source.baseline_q,
                source.baseline_valid,
                folds,
                outer_fold,
                candidates,
            )
            selected_indices[held_index, outer_fold] = selected.order
            forbidden = test_mask
            moments = fit_raw_moments(transformed, ~test_mask, forbidden)
            parameters = posterior_parameters(moments, selected.posterior)
            latent, confidence, support = posterior_prediction(
                transformed,
                (0, 1),
                parameters,
            )
            realized, realized_conf, realized_support, candidate_audit = realize_matrix(
                latent,
                confidence,
                support,
                source.baseline_q[held_index],
                source.baseline_valid[held_index],
                selected.transport,
                test_mask,
            )
            base_latent, base_conf, base_support = equal_mean_agreement(
                transformed,
                (0, 1),
            )
            base_realized, base_conf, base_support, baseline_audit = realize_matrix(
                base_latent,
                base_conf,
                base_support,
                source.baseline_q[held_index],
                source.baseline_valid[held_index],
                selected.transport,
                test_mask,
            )
            cand_metric = score_predictions(
                realized,
                realized_conf,
                realized_support,
                source.effects[held_index],
                source.valid[held_index],
                test_mask,
            )
            base_metric = score_predictions(
                base_realized,
                base_conf,
                base_support,
                source.effects[held_index],
                source.valid[held_index],
                test_mask,
            )
            credit_delta, precision_delta, evaluable = metric_delta(
                cand_metric,
                base_metric,
            )
            per_target_credit_delta[held_index, evaluable] = credit_delta[evaluable]
            per_target_precision_delta[held_index, evaluable] = precision_delta[evaluable]
            candidate_effect[held_index, test_mask] = realized[test_mask]
            candidate_confidence[held_index, test_mask] = realized_conf[test_mask]
            baseline_effect[held_index, test_mask] = base_realized[test_mask]
            baseline_confidence[held_index, test_mask] = base_conf[test_mask]
            total_candidate_support += int(np.count_nonzero(realized_support[test_mask]))
            total_baseline_support += int(np.count_nonzero(base_support[test_mask]))
            transport_invariant_rows += int(candidate_audit["invariant_rows"])
            transport_invariant_rows += int(baseline_audit["invariant_rows"])
            transport_total_rows += 2 * int(candidate_audit["rows"])
            preclip_negative += int(candidate_audit["preclip_negative_coordinates"])
            preclip_negative += int(baseline_audit["preclip_negative_coordinates"])
            if detailed:
                fold_records.append(
                    {
                        "held_family": source.family_names[held_index],
                        "outer_target_fold": outer_fold,
                        "selected_candidate": selected.name,
                        "selected_candidate_index": selected.order,
                        "evaluable_targets": int(np.count_nonzero(evaluable)),
                        "signed_credit_delta": float(np.nanmean(credit_delta)),
                        "precision_delta": float(np.nanmean(precision_delta)),
                        "fit_destination_rows": moments.fit_destination_rows,
                        "fit_biological_rows": moments.fit_biological_rows,
                        "variance_floor": parameters.floor,
                        "candidate_transport": candidate_audit,
                        "baseline_transport": baseline_audit,
                        "inner_candidate_summaries": inner_summaries,
                    }
                )
        family_credit_delta[held_index] = float(
            np.nanmean(per_target_credit_delta[held_index])
        )
        family_precision_delta[held_index] = float(
            np.nanmean(per_target_precision_delta[held_index])
        )
    require(np.all(np.isfinite(family_credit_delta)), "missing family credit")
    require(np.all(np.isfinite(family_precision_delta)), "missing family precision")
    require(total_baseline_support > 0, "zero comparator support")
    return {
        "macro_credit_delta": float(np.mean(family_credit_delta)),
        "macro_precision_delta": float(np.mean(family_precision_delta)),
        "family_credit_delta": family_credit_delta,
        "family_precision_delta": family_precision_delta,
        "per_target_credit_delta": per_target_credit_delta,
        "per_target_precision_delta": per_target_precision_delta,
        "candidate_effect": candidate_effect,
        "candidate_confidence": candidate_confidence,
        "baseline_effect": baseline_effect,
        "baseline_confidence": baseline_confidence,
        "selected_indices": selected_indices,
        "support_retention": float(total_candidate_support / total_baseline_support),
        "candidate_supported_pairs": total_candidate_support,
        "baseline_supported_pairs": total_baseline_support,
        "transport_invariant_rows": transport_invariant_rows,
        "transport_total_rows": transport_total_rows,
        "preclip_negative_coordinates": preclip_negative,
        "fold_records": fold_records if detailed else None,
    }


def bootstrap_distribution(
    per_target: np.ndarray,
    bank: np.ndarray,
) -> np.ndarray:
    values = np.empty(len(bank), dtype=np.float64)
    for draw, indices in enumerate(bank):
        family_means = []
        for family in range(3):
            sampled = per_target[family, indices]
            finite = sampled[np.isfinite(sampled)]
            require(finite.size > 0, "bootstrap family has no evaluable target")
            family_means.append(float(np.mean(finite)))
        values[draw] = float(np.mean(family_means))
    return values


def concentration_report(per_target: np.ndarray) -> dict[str, Any]:
    family_delta = np.nanmean(per_target, axis=1)
    positive_family = np.maximum(family_delta, 0.0)
    family_denominator = float(positive_family.sum())
    family_share = (
        float(positive_family.max() / family_denominator)
        if family_denominator > 0.0
        else 1.0
    )
    contributions = np.zeros(300, dtype=np.float64)
    for family in range(3):
        valid = np.isfinite(per_target[family])
        count = int(np.count_nonzero(valid))
        require(count > 0, "family has no target contributions")
        contributions[valid] += per_target[family, valid] / (3.0 * count)
    positive_target = np.maximum(contributions, 0.0)
    target_denominator = float(positive_target.sum())
    target_share = (
        float(positive_target.max() / target_denominator)
        if target_denominator > 0.0
        else 1.0
    )
    return {
        "family_positive_gain_share_max": family_share,
        "target_cluster_positive_gain_share_max": target_share,
        "family_delta": family_delta.tolist(),
        "target_contribution_sum": float(contributions.sum()),
    }


_WORKER_SOURCE: SourceData | None = None
_WORKER_FOLDS: np.ndarray | None = None
_WORKER_CANDIDATES: tuple[Candidate, ...] | None = None


def init_worker(
    source: SourceData,
    folds: np.ndarray,
    candidates: tuple[Candidate, ...],
) -> None:
    global _WORKER_SOURCE, _WORKER_FOLDS, _WORKER_CANDIDATES
    _WORKER_SOURCE = source
    _WORKER_FOLDS = folds
    _WORKER_CANDIDATES = candidates


def null_worker(task: tuple[str, int, np.ndarray]) -> tuple[str, int, float, float]:
    kind, index, permutation = task
    require(_WORKER_SOURCE is not None, "worker source missing")
    require(_WORKER_FOLDS is not None, "worker folds missing")
    require(_WORKER_CANDIDATES is not None, "worker candidates missing")
    kwargs: dict[str, np.ndarray] = {}
    if kind == "target_label":
        kwargs["target_permutation"] = permutation
    elif kind == "response_label":
        kwargs["response_permutation"] = permutation
    else:
        raise ExperimentError(f"unknown null kind: {kind}")
    result = nested_evaluate(
        _WORKER_SOURCE,
        _WORKER_FOLDS,
        _WORKER_CANDIDATES,
        detailed=False,
        **kwargs,
    )
    return (
        kind,
        index,
        float(result["macro_credit_delta"]),
        float(result["macro_precision_delta"]),
    )


def run_nulls(
    source: SourceData,
    folds: np.ndarray,
    candidates: tuple[Candidate, ...],
    target_bank: np.ndarray,
    response_bank: np.ndarray,
    workers: int,
) -> dict[str, np.ndarray]:
    outputs = {
        "target_credit": np.empty(len(target_bank), dtype=np.float64),
        "target_precision": np.empty(len(target_bank), dtype=np.float64),
        "response_credit": np.empty(len(response_bank), dtype=np.float64),
        "response_precision": np.empty(len(response_bank), dtype=np.float64),
    }
    tasks = [
        ("target_label", index, row)
        for index, row in enumerate(target_bank)
    ] + [
        ("response_label", index, row)
        for index, row in enumerate(response_bank)
    ]
    completed = 0
    if workers <= 1:
        init_worker(source, folds, candidates)
        iterator = map(null_worker, tasks)
        for kind, index, credit, precision in iterator:
            prefix = "target" if kind == "target_label" else "response"
            outputs[f"{prefix}_credit"][index] = credit
            outputs[f"{prefix}_precision"][index] = precision
            completed += 1
            if completed % 20 == 0:
                print(f"  complete-procedure nulls: {completed}/{len(tasks)}", flush=True)
    else:
        context = multiprocessing.get_context("fork")
        with ProcessPoolExecutor(
            max_workers=workers,
            mp_context=context,
            initializer=init_worker,
            initargs=(source, folds, candidates),
        ) as executor:
            futures = [executor.submit(null_worker, task) for task in tasks]
            for future in as_completed(futures):
                kind, index, credit, precision = future.result()
                prefix = "target" if kind == "target_label" else "response"
                outputs[f"{prefix}_credit"][index] = credit
                outputs[f"{prefix}_precision"][index] = precision
                completed += 1
                if completed % 20 == 0:
                    print(
                        f"  complete-procedure nulls: {completed}/{len(tasks)}",
                        flush=True,
                    )
    return outputs


def quantile(array: np.ndarray, probability: float) -> float:
    return float(np.quantile(array, probability, method="linear"))


def p_value(null: np.ndarray, observed: float) -> float:
    return float((1 + np.count_nonzero(null >= observed)) / (len(null) + 1))


def expanded_pds_invariants(
    effects: np.ndarray,
    panel: Sequence[str],
    genes: Sequence[str],
) -> dict[str, Any]:
    gene_position = {gene: index for index, gene in enumerate(genes)}
    columns = np.asarray([gene_position[gene] for gene in panel], dtype=np.int64)
    nonpanel = np.ones(len(genes), dtype=bool)
    nonpanel[columns] = False
    expanded = np.zeros((effects.shape[0] * effects.shape[1], len(genes)), dtype=np.float32)
    expanded[:, columns] = effects.reshape(-1, 300).astype(np.float32)
    diagonal = np.arange(300)
    checks = {
        "finite": bool(np.all(np.isfinite(effects))),
        "target_self_effect_exact_zero": bool(
            np.all(effects[:, diagonal, diagonal] == 0.0)
        ),
        "expanded_nonpanel_effect_exact_zero": bool(
            np.all(expanded[:, nonpanel] == 0.0)
        ),
        "operation_writes_only_panel_coordinates": True,
    }
    checks["all_pass"] = all(checks.values())
    return checks


def safe_output(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    require(resolved != HERE and resolved.is_relative_to(HERE), "unsafe output path")
    return resolved


def write_npy(path: Path, value: np.ndarray) -> None:
    with path.open("xb") as handle:
        np.save(handle, np.ascontiguousarray(value), allow_pickle=False)


def runtime_versions() -> dict[str, str]:
    import scipy

    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }


def run(output: Path = DEFAULT_OUTPUT, *, workers: int = 1) -> Path:
    output = safe_output(output)
    require(not output.exists(), f"refusing to overwrite {output}")
    protocol = read_json(PROTOCOL_PATH)
    candidates = candidates_from_protocol(protocol)
    print("Verifying complete preregistration and implementation freeze...", flush=True)
    source = load_source_data()
    freeze_plans.verify(PLAN_DIR)
    folds = np.load(PLAN_DIR / "target-folds.npy", allow_pickle=False)
    target_bank = np.load(
        PLAN_DIR / "target-label-derangements.npy",
        allow_pickle=False,
    )
    response_bank = np.load(
        PLAN_DIR / "response-label-derangements.npy",
        allow_pickle=False,
    )
    bootstrap_bank = np.load(
        PLAN_DIR / "target-cluster-bootstrap.npy",
        allow_pickle=False,
    )
    print("Running observed nested leave-one-family-out evaluation...", flush=True)
    observed = nested_evaluate(source, folds, candidates, detailed=True)
    print(
        "Observed macro deltas: "
        f"credit={observed['macro_credit_delta']:+.8f}, "
        f"precision={observed['macro_precision_delta']:+.8f}",
        flush=True,
    )
    print(f"Running 400 complete-procedure nulls with {workers} worker(s)...", flush=True)
    nulls = run_nulls(
        source,
        folds,
        candidates,
        target_bank,
        response_bank,
        workers,
    )
    credit_bootstrap = bootstrap_distribution(
        observed["per_target_credit_delta"],
        bootstrap_bank,
    )
    precision_bootstrap = bootstrap_distribution(
        observed["per_target_precision_delta"],
        bootstrap_bank,
    )
    credit_concentration = concentration_report(
        observed["per_target_credit_delta"]
    )
    precision_concentration = concentration_report(
        observed["per_target_precision_delta"]
    )
    pds_invariants = expanded_pds_invariants(
        observed["candidate_effect"],
        source.panel,
        source.genes,
    )
    require(pds_invariants["all_pass"], "expanded PDS invariant failed")
    require(
        observed["transport_invariant_rows"] == observed["transport_total_rows"],
        "transport row invariant failure",
    )

    family_credit = np.asarray(observed["family_credit_delta"])
    family_precision = np.asarray(observed["family_precision_delta"])
    gates = {
        "nonnegative_signed_credit_delta_every_family": bool(np.all(family_credit >= 0.0)),
        "nonnegative_precision_delta_every_family": bool(np.all(family_precision >= 0.0)),
        "macro_signed_credit_gain_minimum": bool(observed["macro_credit_delta"] >= 0.005),
        "macro_precision_gain_minimum": bool(observed["macro_precision_delta"] >= 0.01),
        "signed_credit_bootstrap_lower_95_strictly_positive": bool(
            quantile(credit_bootstrap, 0.025) > 0.0
        ),
        "precision_bootstrap_lower_95_strictly_positive": bool(
            quantile(precision_bootstrap, 0.025) > 0.0
        ),
        "signed_credit_gain_strictly_above_target_null_p95": bool(
            observed["macro_credit_delta"] > quantile(nulls["target_credit"], 0.95)
        ),
        "signed_credit_gain_strictly_above_response_null_p95": bool(
            observed["macro_credit_delta"] > quantile(nulls["response_credit"], 0.95)
        ),
        "precision_gain_strictly_above_target_null_p95": bool(
            observed["macro_precision_delta"]
            > quantile(nulls["target_precision"], 0.95)
        ),
        "precision_gain_strictly_above_response_null_p95": bool(
            observed["macro_precision_delta"]
            > quantile(nulls["response_precision"], 0.95)
        ),
        "supported_pair_retention_minimum": bool(observed["support_retention"] >= 0.8),
        "maximum_family_positive_gain_share_credit": bool(
            credit_concentration["family_positive_gain_share_max"] <= 0.5
        ),
        "maximum_family_positive_gain_share_precision": bool(
            precision_concentration["family_positive_gain_share_max"] <= 0.5
        ),
        "maximum_target_cluster_positive_gain_share_credit": bool(
            credit_concentration["target_cluster_positive_gain_share_max"] <= 0.5
        ),
        "maximum_target_cluster_positive_gain_share_precision": bool(
            precision_concentration["target_cluster_positive_gain_share_max"] <= 0.5
        ),
        "transport_and_pds_invariants": bool(
            pds_invariants["all_pass"]
            and observed["transport_invariant_rows"] == observed["transport_total_rows"]
        ),
    }
    promotion = all(gates.values())
    selected_counts = {
        candidate.name: int(np.count_nonzero(observed["selected_indices"] == candidate.order))
        for candidate in candidates
    }
    results = {
        "schema_version": 1,
        "experiment": "125-hierarchical-panel-posterior",
        "phase": "source_only_nested_evaluation",
        "decision": "GO" if promotion else "NO-GO",
        "promotion_passes": promotion,
        "candidate_order": [candidate.name for candidate in candidates],
        "candidate_definitions": [
            {
                "name": candidate.name,
                "posterior": asdict(candidate.posterior),
                "transport": candidate.transport,
                "order": candidate.order,
            }
            for candidate in candidates
        ],
        "selected_candidate_counts_across_15_outer_target_folds": selected_counts,
        "observed": {
            "macro_signed_credit_delta": observed["macro_credit_delta"],
            "macro_precision_delta": observed["macro_precision_delta"],
            "family_signed_credit_delta": dict(
                zip(source.family_names, family_credit.tolist(), strict=True)
            ),
            "family_precision_delta": dict(
                zip(source.family_names, family_precision.tolist(), strict=True)
            ),
            "support_retention": observed["support_retention"],
            "candidate_supported_pairs": observed["candidate_supported_pairs"],
            "baseline_supported_pairs": observed["baseline_supported_pairs"],
            "preclip_negative_coordinates_candidate_plus_comparator": observed[
                "preclip_negative_coordinates"
            ],
            "folds": observed["fold_records"],
        },
        "bootstrap": {
            "replicates": len(bootstrap_bank),
            "signed_credit_delta_lower_95": quantile(credit_bootstrap, 0.025),
            "signed_credit_delta_median": quantile(credit_bootstrap, 0.5),
            "signed_credit_delta_upper_95": quantile(credit_bootstrap, 0.975),
            "precision_delta_lower_95": quantile(precision_bootstrap, 0.025),
            "precision_delta_median": quantile(precision_bootstrap, 0.5),
            "precision_delta_upper_95": quantile(precision_bootstrap, 0.975),
        },
        "complete_procedure_nulls": {
            "target_label": {
                "draws": len(target_bank),
                "signed_credit_delta_p95": quantile(nulls["target_credit"], 0.95),
                "precision_delta_p95": quantile(nulls["target_precision"], 0.95),
                "signed_credit_p_value": p_value(
                    nulls["target_credit"],
                    observed["macro_credit_delta"],
                ),
                "precision_p_value": p_value(
                    nulls["target_precision"],
                    observed["macro_precision_delta"],
                ),
            },
            "response_label_within_expression_bin": {
                "draws": len(response_bank),
                "signed_credit_delta_p95": quantile(nulls["response_credit"], 0.95),
                "precision_delta_p95": quantile(nulls["response_precision"], 0.95),
                "signed_credit_p_value": p_value(
                    nulls["response_credit"],
                    observed["macro_credit_delta"],
                ),
                "precision_p_value": p_value(
                    nulls["response_precision"],
                    observed["macro_precision_delta"],
                ),
            },
        },
        "concentration": {
            "signed_credit": credit_concentration,
            "precision": precision_concentration,
        },
        "invariants": {
            "transport_rows_passed": observed["transport_invariant_rows"],
            "transport_rows_total": observed["transport_total_rows"],
            "expanded_pds": pds_invariants,
        },
        "gates": gates,
        "safety": {
            "pan_accessed": False,
            "dld1_accessed": False,
            "h1_accessed": False,
            "orion_held_harness_outcomes_accessed": False,
            "challenge_truth_accessed": False,
            "challenge_context_expression_values_accessed": False,
            "cell_eval2_run": False,
            "counts_or_submission_built": False,
            "submission_sent": False,
            "commit_or_push": False,
        },
        "runtime_versions": runtime_versions(),
        "freeze_manifest_sha256": sha256_file(PLAN_DIR / "freeze-manifest.json"),
        "implementation_freeze_sha256": sha256_file(IMPLEMENTATION_FREEZE_PATH),
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    stage = output.parent / f".{output.name}.staging"
    require(not stage.exists(), f"staging directory exists: {stage}")
    stage.mkdir()
    try:
        arrays = {
            "observed-candidate-effect.npy": observed["candidate_effect"].astype(np.float32),
            "observed-candidate-confidence.npy": observed["candidate_confidence"].astype(np.float32),
            "observed-baseline-effect.npy": observed["baseline_effect"].astype(np.float32),
            "observed-baseline-confidence.npy": observed["baseline_confidence"].astype(np.float32),
            "per-target-credit-delta.npy": observed["per_target_credit_delta"],
            "per-target-precision-delta.npy": observed["per_target_precision_delta"],
            "selected-candidate-index.npy": observed["selected_indices"],
            "family-baseline-q.npy": source.baseline_q,
            "family-baseline-valid.npy": source.baseline_valid,
            "target-null-credit-delta.npy": nulls["target_credit"],
            "target-null-precision-delta.npy": nulls["target_precision"],
            "response-null-credit-delta.npy": nulls["response_credit"],
            "response-null-precision-delta.npy": nulls["response_precision"],
            "credit-bootstrap-delta.npy": credit_bootstrap,
            "precision-bootstrap-delta.npy": precision_bootstrap,
        }
        for name, array in arrays.items():
            write_npy(stage / name, array)
        (stage / "results.json").write_bytes(json_bytes(results))
        files = {
            path.name: {
                "sha256": sha256_file(path),
                "logical_array_sha256": (
                    sha256_array(arrays[path.name]) if path.name in arrays else None
                ),
                "size_bytes": path.stat().st_size,
            }
            for path in sorted(stage.iterdir())
        }
        manifest = {
            "schema_version": 1,
            "experiment": "125-hierarchical-panel-posterior",
            "phase": "source_only_nested_evaluation",
            "decision": results["decision"],
            "protocol_sha256": sha256_file(PROTOCOL_PATH),
            "authorization_sha256": sha256_file(AUTHORIZATION_PATH),
            "implementation_freeze_sha256": sha256_file(IMPLEMENTATION_FREEZE_PATH),
            "plan_freeze_manifest_sha256": sha256_file(
                PLAN_DIR / "freeze-manifest.json"
            ),
            "files": files,
        }
        (stage / "artifact-manifest.json").write_bytes(json_bytes(manifest))
        receipt = {
            "schema_version": 1,
            "experiment": "125-hierarchical-panel-posterior",
            "status": "immutable_source_only_result",
            "decision": results["decision"],
            "artifact_manifest_sha256": sha256_file(stage / "artifact-manifest.json"),
            "results_sha256": sha256_file(stage / "results.json"),
            "held_outcomes_accessed": False,
            "cell_eval2_run": False,
            "counts_or_submission_built": False,
            "submission_sent": False,
            "commit_or_push": False,
        }
        (stage / "result-receipt.json").write_bytes(json_bytes(receipt))
        (stage / "result-receipt.sha256").write_text(
            f"{sha256_file(stage / 'result-receipt.json')}  result-receipt.json\n",
            encoding="utf-8",
        )
        os.sync()
        stage.rename(output)
        for path in output.iterdir():
            path.chmod(0o444)
        output.chmod(0o555)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return output


def audit() -> dict[str, Any]:
    authorization = validate_authorization(hash_effect_arrays=False)
    protocol = read_json(PROTOCOL_PATH)
    candidates = candidates_from_protocol(protocol)
    return {
        "authorized": True,
        "allowed_input_count": len(authorization["allowed_input_files"]),
        "candidate_count": len(candidates),
        "candidate_names": [candidate.name for candidate in candidates],
        "family_names": list(FAMILY_NAMES),
        "plan": freeze_plans.verify(PLAN_DIR),
        "output_exists": DEFAULT_OUTPUT.exists(),
        "runtime_versions": runtime_versions(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--workers", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    parser.add_argument("--audit", action="store_true")
    args = parser.parse_args()
    if args.audit:
        print(json.dumps(audit(), indent=2, sort_keys=True))
    else:
        print(
            json.dumps(
                {"output": str(run(args.output, workers=max(1, args.workers)))},
                indent=2,
                sort_keys=True,
            )
        )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Source-only Experiment 118 builder.

This module intentionally has no held-evaluation or count-generation entry
point. Its only biological inputs are enumerated in SOURCE_ONLY_AUTHORIZATION.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
from scipy.special import erf


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
AUTHORIZATION_PATH = HERE / "SOURCE_ONLY_AUTHORIZATION.json"
OUTPUT = HERE / "work" / "source-freeze-v1"
CONDITION_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
ESTIMATORS = ("unweighted", "inverse_variance", "r_mu", "r2_mu")
ETAS = (0.25, 0.50, 1.00)
COMPLEXITY = {name: i for i, name in enumerate(ESTIMATORS)}
NULL_DRAWS = 100
NULL_SEED = 1182026


class ContractError(RuntimeError):
    """Raised when a frozen source-only contract is violated."""


@dataclass(frozen=True)
class Candidate:
    name: str
    estimator: str
    eta: float


@dataclass
class Condition:
    name: str
    lfc: np.ndarray
    se: np.ndarray
    fdr: np.ndarray
    source_valid: np.ndarray
    truth_valid: np.ndarray
    target_qc: np.ndarray


def candidates() -> list[Candidate]:
    out = [Candidate("identity", "identity", 0.0)]
    for estimator in ESTIMATORS:
        for eta in ETAS:
            eta_label = f"{eta:.2f}".replace(".", "p")
            out.append(Candidate(f"{estimator}_eta{eta_label}", estimator, eta))
    return out


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def sha256_array(array: np.ndarray) -> str:
    a = np.ascontiguousarray(array)
    h = hashlib.sha256()
    h.update(str(a.dtype).encode())
    h.update(b"\0")
    h.update(json.dumps(list(a.shape), separators=(",", ":")).encode())
    h.update(b"\0")
    h.update(a.tobytes(order="C"))
    return h.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ContractError(f"{path}: expected a JSON object")
    return value


def atomic_json(path: Path, value: Any) -> None:
    temp = path.with_name(path.name + ".tmp")
    if path.exists() or temp.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    with temp.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    temp.rename(path)


def atomic_npy(path: Path, value: np.ndarray) -> None:
    temp = path.with_name(path.name + ".tmp")
    if path.exists() or temp.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    with temp.open("xb") as handle:
        np.save(handle, np.ascontiguousarray(value), allow_pickle=False)
        handle.flush()
        os.fsync(handle.fileno())
    temp.rename(path)


def read_axis(path: Path, column: str) -> list[str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    values = [str(row[column]) for row in rows]
    if not values or len(values) != len(set(values)):
        raise ContractError(f"{path}: axis is empty or contains duplicates")
    return values


def parse_bool(value: str) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    raise ContractError(f"unexpected boolean value {value!r}")


def validate_protocol() -> dict[str, Any]:
    protocol = read_json(PROTOCOL_PATH)
    expected = [(c.estimator, c.eta) for c in candidates()]
    observed = [("identity", float(protocol["candidate_grid"]["identity"]["eta"]))]
    observed.extend(
        (estimator, float(eta))
        for estimator in protocol["candidate_grid"]["estimators"]
        for eta in protocol["candidate_grid"]["eta"]
    )
    if observed != expected:
        raise ContractError("protocol candidate grid does not match implementation")
    if protocol["nulls"]["seed"] != NULL_SEED:
        raise ContractError("protocol null seed mismatch")
    for family in ("target_label", "response_label_within_expression_bin"):
        if protocol["nulls"][family]["repetitions"] != NULL_DRAWS:
            raise ContractError(f"protocol {family} null count mismatch")
    return protocol


def validate_authorization() -> dict[str, Any]:
    protocol = validate_protocol()
    authorization = read_json(AUTHORIZATION_PATH)
    if authorization.get("source_only_authorized") is not True:
        raise ContractError("source-only build is not authorized")
    if authorization.get("scope") != "experiment_118_source_only_first_stage":
        raise ContractError("authorization scope mismatch")
    forbidden = authorization.get("external_actions", {})
    if not forbidden or any(value is not False for value in forbidden.values()):
        raise ContractError("every external/held action must remain false")

    allowed = authorization.get("allowed_input_files")
    if not isinstance(allowed, dict) or not allowed:
        raise ContractError("authorization has no pinned input files")
    for relative, expected in sorted(allowed.items()):
        path = (REPO / relative).resolve()
        try:
            path.relative_to(REPO)
        except ValueError as exc:
            raise ContractError(f"input escapes repository: {relative}") from exc
        if not path.is_file():
            raise ContractError(f"missing authorized input: {relative}")
        observed = sha256_file(path)
        if observed != expected:
            raise ContractError(
                f"authorized input hash mismatch: {relative}: {observed} != {expected}"
            )

    implementation = authorization.get("implementation_files")
    if not isinstance(implementation, dict) or not implementation:
        raise ContractError("authorization has no implementation hashes")
    for relative, expected in sorted(implementation.items()):
        path = (REPO / relative).resolve()
        if not path.is_file() or sha256_file(path) != expected:
            raise ContractError(f"implementation authorization mismatch: {relative}")

    if authorization.get("protocol_sha256") != sha256_file(PROTOCOL_PATH):
        raise ContractError("protocol authorization hash mismatch")
    if authorization.get("candidate_names") != [c.name for c in candidates()]:
        raise ContractError("authorization candidate ordering mismatch")
    if protocol["freeze"]["held_access_authorized"] is not False:
        raise ContractError("protocol unexpectedly authorizes held access")
    return authorization


def condition_root(name: str) -> Path:
    return (
        REPO
        / "data"
        / "external"
        / "cd4t"
        / "cloud-compact"
        / "cd4_compact"
        / name
    )


def target_qc_by_name(path: Path) -> dict[str, bool]:
    out: dict[str, bool] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            target = row["target_gene"]
            if target in out:
                raise ContractError(f"{path}: duplicate target {target}")
            out[target] = (
                parse_bool(row["ontarget_significant"])
                and not parse_bool(row["neighboring_gene_KD"])
                and not parse_bool(row["distal_offtarget_flag"])
                and not parse_bool(row["low_target_gex"])
                and float(row["n_guides"]) >= 2.0
            )
    return out


def source_axes(name: str) -> tuple[list[str], list[str]]:
    root = condition_root(name)
    return (
        (root / "perts.txt").read_text().splitlines(),
        (root / "genes.txt").read_text().splitlines(),
    )


def panel_source_indices(
    source_values: Sequence[str], canonical_values: Sequence[str]
) -> tuple[np.ndarray, np.ndarray]:
    lookup = {value: index for index, value in enumerate(source_values)}
    canonical_indices: list[int] = []
    source_indices: list[int] = []
    for canonical_index, value in enumerate(canonical_values):
        if value in lookup:
            canonical_indices.append(canonical_index)
            source_indices.append(lookup[value])
    return (
        np.asarray(canonical_indices, dtype=np.int64),
        np.asarray(source_indices, dtype=np.int64),
    )


def build_expression_bins(panel: Sequence[str], n_bins: int = 10) -> np.ndarray:
    pooled: list[list[np.ndarray]] = [[] for _ in panel]
    for condition_name in CONDITION_NAMES:
        root = condition_root(condition_name)
        _, source_genes = source_axes(condition_name)
        canonical_cols, source_cols = panel_source_indices(source_genes, panel)
        base = np.load(root / "base_mean.npy", mmap_mode="r", allow_pickle=False)
        if base.ndim != 2 or base.shape[1] != len(source_genes):
            raise ContractError(f"{condition_name}: base_mean shape mismatch")
        for canonical_col, source_col in zip(canonical_cols, source_cols, strict=True):
            values = np.asarray(base[:, source_col], dtype=np.float64)
            values = values[np.isfinite(values) & (values >= 0)]
            if values.size:
                pooled[int(canonical_col)].append(np.log1p(values))

    expression = np.full(len(panel), np.nan, dtype=np.float64)
    for index, blocks in enumerate(pooled):
        if blocks:
            expression[index] = float(np.median(np.concatenate(blocks)))
    observed = np.flatnonzero(np.isfinite(expression))
    order = observed[np.lexsort((observed, expression[observed]))]
    bins = np.full(len(panel), -1, dtype=np.int16)
    for bin_index, members in enumerate(np.array_split(order, n_bins)):
        bins[members] = bin_index
    if np.any(np.bincount(bins[bins >= 0], minlength=n_bins) < 2):
        raise ContractError("expression bins must each contain at least two genes")
    return bins


def random_derangement(rng: np.random.Generator, size: int) -> np.ndarray:
    identity = np.arange(size, dtype=np.int16)
    for _ in range(10000):
        candidate = rng.permutation(size).astype(np.int16)
        if np.all(candidate != identity):
            return candidate
    raise ContractError(f"failed to draw a derangement of size {size}")


def freeze_permutations(expression_bins: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(NULL_SEED)
    n = len(expression_bins)
    target = np.empty((NULL_DRAWS, n), dtype=np.int16)
    response = np.tile(np.arange(n, dtype=np.int16), (NULL_DRAWS, 1))
    for draw in range(NULL_DRAWS):
        target[draw] = random_derangement(rng, n)
    for draw in range(NULL_DRAWS):
        for bin_index in sorted(set(map(int, expression_bins)) - {-1}):
            members = np.flatnonzero(expression_bins == bin_index)
            local = random_derangement(rng, len(members))
            response[draw, members] = members[local]
    validate_permutations(target, response, expression_bins)
    return target, response


def validate_permutations(
    target: np.ndarray, response: np.ndarray, expression_bins: np.ndarray
) -> None:
    n = len(expression_bins)
    identity = np.arange(n)
    if target.shape != (NULL_DRAWS, n) or response.shape != (NULL_DRAWS, n):
        raise ContractError("permutation bank shape mismatch")
    for name, bank in (("target", target), ("response", response)):
        if not np.issubdtype(bank.dtype, np.integer):
            raise ContractError(f"{name} permutations must be integer")
        if not all(np.array_equal(np.sort(row), identity) for row in bank):
            raise ContractError(f"{name} permutation is not bijective")
    if np.any(target == identity):
        raise ContractError("target permutation fixed point")
    for row in response:
        if np.any(row[expression_bins < 0] != identity[expression_bins < 0]):
            raise ContractError("source-absent response gene moved")
        for bin_index in sorted(set(map(int, expression_bins)) - {-1}):
            members = np.flatnonzero(expression_bins == bin_index)
            if np.any(row[members] == members):
                raise ContractError("response permutation fixed point")
            if set(map(int, row[members])) != set(map(int, members)):
                raise ContractError("response permutation crossed expression bins")


def load_condition(
    name: str, panel: Sequence[str], *, permit_effect_values: bool
) -> Condition:
    if not permit_effect_values:
        raise ContractError("effect arrays require explicit source-only build authorization")
    root = condition_root(name)
    source_targets, source_genes = source_axes(name)
    target_cols, source_rows = panel_source_indices(source_targets, panel)
    gene_cols, source_cols = panel_source_indices(source_genes, panel)
    target_qc_map = target_qc_by_name(root / "obs.csv")

    shape = (len(panel), len(panel))
    lfc = np.full(shape, np.nan, dtype=np.float64)
    se = np.full(shape, np.nan, dtype=np.float64)
    fdr = np.full(shape, np.nan, dtype=np.float64)
    valid_flag = np.zeros(shape, dtype=bool)
    arrays = {
        "lfc": np.load(root / "log_fc.npy", mmap_mode="r", allow_pickle=False),
        "se": np.load(root / "lfc_se.npy", mmap_mode="r", allow_pickle=False),
        "fdr": np.load(root / "adj_p_value.npy", mmap_mode="r", allow_pickle=False),
        "valid": np.load(
            root / "valid_log_fc.npy", mmap_mode="r", allow_pickle=False
        ),
    }
    expected_shape = (len(source_targets), len(source_genes))
    if any(array.shape != expected_shape for array in arrays.values()):
        raise ContractError(f"{name}: compact array shape mismatch")
    destination = np.ix_(target_cols, gene_cols)
    source = np.ix_(source_rows, source_cols)
    lfc[destination] = np.asarray(arrays["lfc"][source], dtype=np.float64)
    se[destination] = np.asarray(arrays["se"][source], dtype=np.float64)
    fdr[destination] = np.asarray(arrays["fdr"][source], dtype=np.float64)
    valid_flag[destination] = np.asarray(arrays["valid"][source], dtype=bool)
    target_qc = np.asarray(
        [target_qc_map.get(target, False) for target in panel], dtype=bool
    )
    source_valid = (
        target_qc[:, None]
        & valid_flag
        & np.isfinite(lfc)
        & np.isfinite(se)
        & (se > 0)
        & np.isfinite(fdr)
    )
    diagonal = np.arange(len(panel))
    source_valid[diagonal, diagonal] = False
    truth_valid = source_valid & (fdr < 0.05) & (lfc != 0)
    return Condition(name, lfc, se, fdr, source_valid, truth_valid, target_qc)


def repair_zero_sum(row: np.ndarray, mutable: np.ndarray) -> None:
    indices = np.flatnonzero(mutable)
    if indices.size == 0:
        return
    # Prefer the final canonical mutable coordinate. If the correctly rounded
    # residual is below its ULP, continue deterministically through coordinates
    # in ascending magnitude so the correction remains representable.
    order = [int(indices[-1])]
    order.extend(
        int(index)
        for index in indices[np.argsort(np.abs(row[indices]), kind="stable")]
        if int(index) != order[0]
    )
    for repair_index in order:
        for _ in range(4):
            residual = math.fsum(map(float, row))
            if residual == 0.0:
                return
            before = float(row[repair_index])
            row[repair_index] -= residual
            if float(row[repair_index]) == before:
                break
    if math.fsum(map(float, row)) != 0.0:
        raise ContractError("failed exact compensated-float64 zero-sum repair")


def rows_have_exact_zero_sum(matrix: np.ndarray) -> bool:
    return all(math.fsum(map(float, row)) == 0.0 for row in matrix)


def estimate_base(
    conditions: Sequence[Condition],
    estimator: str,
    *,
    target_permutation: np.ndarray | None = None,
    response_permutation: np.ndarray | None = None,
) -> np.ndarray:
    if estimator not in ESTIMATORS:
        raise ContractError(f"unknown estimator {estimator}")
    n = conditions[0].lfc.shape[0]
    target_map = (
        np.arange(n, dtype=np.int64)
        if target_permutation is None
        else np.asarray(target_permutation, dtype=np.int64)
    )
    response_map = (
        np.arange(n, dtype=np.int64)
        if response_permutation is None
        else np.asarray(response_permutation, dtype=np.int64)
    )
    lfc = np.stack([c.lfc[target_map][:, response_map] for c in conditions])
    se = np.stack([c.se[target_map][:, response_map] for c in conditions])
    valid = np.stack(
        [c.source_valid[target_map][:, response_map] for c in conditions]
    )
    diagonal = np.arange(n)
    valid[:, diagonal, diagonal] = False
    count = valid.sum(axis=0)
    supported = count >= 2
    safe_lfc = np.where(valid, lfc, 0.0)

    if estimator == "unweighted":
        raw = np.divide(
            safe_lfc.sum(axis=0),
            count,
            out=np.zeros((n, n), dtype=np.float64),
            where=count > 0,
        )
    else:
        weights = np.where(valid, 1.0 / np.square(se), 0.0)
        weight_sum = weights.sum(axis=0)
        mu = np.divide(
            (weights * safe_lfc).sum(axis=0),
            weight_sum,
            out=np.zeros((n, n), dtype=np.float64),
            where=weight_sum > 0,
        )
        if estimator == "inverse_variance":
            raw = mu
        else:
            z = np.abs(mu) * np.sqrt(weight_sum)
            agreement = np.divide(
                np.abs((weights * np.sign(safe_lfc)).sum(axis=0)),
                weight_sum,
                out=np.zeros((n, n), dtype=np.float64),
                where=weight_sum > 0,
            )
            reliability = erf(z / math.sqrt(2.0)) * agreement
            raw = reliability * mu
            if estimator == "r2_mu":
                raw *= reliability

    out = np.where(supported, raw, 0.0)
    out[diagonal, diagonal] = 0.0
    supported[diagonal, diagonal] = False
    for target_index in range(n):
        mutable = supported[target_index]
        if np.any(mutable):
            out[target_index, mutable] -= float(np.mean(out[target_index, mutable]))
            repair_zero_sum(out[target_index], mutable)
    out[~supported] = 0.0
    out[diagonal, diagonal] = 0.0
    if not np.all(np.isfinite(out)):
        raise ContractError("candidate contains nonfinite values")
    if np.any(out[~supported] != 0):
        raise ContractError("candidate changed unsupported coordinates")
    if np.any(np.diag(out) != 0):
        raise ContractError("candidate changed target self coordinate")
    if not rows_have_exact_zero_sum(out):
        raise ContractError("candidate panel rows are not exact zero-sum")
    return out


def per_fold_metrics(prediction: np.ndarray, held: Condition) -> dict[str, Any]:
    credits: list[float] = []
    precisions: list[float] = []
    pair_count = 0
    n = prediction.shape[0]
    for target_index in range(n):
        valid = held.truth_valid[target_index].copy()
        valid[target_index] = False
        indices = np.flatnonzero(valid)
        if indices.size == 0:
            continue
        truth = held.lfc[target_index, indices]
        pred = prediction[target_index, indices]
        denom = float(np.abs(truth).sum())
        if not np.isfinite(denom) or denom <= 0:
            continue
        credits.append(float(np.sign(truth) @ pred) / denom)
        top_n = max(1, int(math.ceil(0.2 * len(indices))))
        order = np.argsort(-np.abs(pred), kind="stable")[:top_n]
        weights = np.abs(truth[order])
        correct = (np.sign(pred[order]) == np.sign(truth[order])).astype(np.float64)
        precisions.append(float(correct @ weights) / float(weights.sum()))
        pair_count += len(indices)
    if not credits:
        raise ContractError(f"{held.name}: no evaluable held targets")
    return {
        "condition": held.name,
        "target_count": len(credits),
        "truth_pair_count": pair_count,
        "signed_linear_credit": float(np.mean(credits)),
        "top20_weighted_sign_precision": float(np.mean(precisions)),
    }


def evaluate_grid(
    conditions: Sequence[Condition],
    *,
    target_permutation: np.ndarray | None = None,
    response_permutation: np.ndarray | None = None,
    include_all_metrics: bool = True,
) -> tuple[Candidate, list[dict[str, Any]]]:
    grid = candidates()
    metrics: dict[str, list[dict[str, Any]]] = {candidate.name: [] for candidate in grid}
    for held_index, held in enumerate(conditions):
        training = [c for index, c in enumerate(conditions) if index != held_index]
        base = {
            estimator: estimate_base(
                training,
                estimator,
                target_permutation=target_permutation,
                response_permutation=response_permutation,
            )
            for estimator in ESTIMATORS
        }
        metrics["identity"].append(
            per_fold_metrics(np.zeros_like(base[ESTIMATORS[0]]), held)
        )
        for candidate in grid[1:]:
            metrics[candidate.name].append(
                per_fold_metrics(candidate.eta * base[candidate.estimator], held)
            )

    summaries: list[dict[str, Any]] = []
    for candidate in grid:
        folds = metrics[candidate.name]
        credits = [float(fold["signed_linear_credit"]) for fold in folds]
        precisions = [
            float(fold["top20_weighted_sign_precision"]) for fold in folds
        ]
        eligible = candidate.name == "identity" or (
            all(value >= 0.0 for value in credits)
            and all(value >= 0.5 for value in precisions)
        )
        summary = {
            "name": candidate.name,
            "estimator": candidate.estimator,
            "eta": candidate.eta,
            "eligible": eligible,
            "worst_fold_signed_linear_credit": min(credits),
            "worst_fold_top20_weighted_sign_precision": min(precisions),
        }
        if include_all_metrics:
            summary["folds"] = folds
        summaries.append(summary)

    def selection_key(summary: dict[str, Any]) -> tuple[float, float, float, int]:
        if not summary["eligible"]:
            return (-math.inf, -math.inf, -math.inf, -10**9)
        estimator = str(summary["estimator"])
        complexity = -1 if estimator == "identity" else COMPLEXITY[estimator]
        return (
            float(summary["worst_fold_signed_linear_credit"]),
            float(summary["worst_fold_top20_weighted_sign_precision"]),
            -float(summary["eta"]),
            -complexity,
        )

    selected_summary = max(summaries, key=selection_key)
    selected = next(c for c in grid if c.name == selected_summary["name"])
    return selected, summaries


def build_full_candidates(conditions: Sequence[Condition]) -> np.ndarray:
    grid = candidates()
    out = np.zeros((len(grid), 300, 300), dtype=np.float64)
    bases = {
        estimator: estimate_base(conditions, estimator) for estimator in ESTIMATORS
    }
    for index, candidate in enumerate(grid[1:], start=1):
        out[index] = candidate.eta * bases[candidate.estimator]
    if np.any(out[0] != 0):
        raise ContractError("identity candidate is not exact zero")
    if np.any(np.diagonal(out, axis1=1, axis2=2) != 0):
        raise ContractError("candidate target-self coordinate is not exact zero")
    if not all(rows_have_exact_zero_sum(arm) for arm in out):
        raise ContractError("candidate panel residual is not exact zero-sum")
    return out


def null_selection(
    conditions: Sequence[Condition],
    bank: np.ndarray,
    family: str,
) -> dict[str, Any]:
    selected_names: list[str] = []
    statistics = np.empty(len(bank), dtype=np.float64)
    for draw, permutation in enumerate(bank):
        kwargs: dict[str, np.ndarray] = {}
        if family == "target_label":
            kwargs["target_permutation"] = permutation
        elif family == "response_label_within_expression_bin":
            kwargs["response_permutation"] = permutation
        else:
            raise ContractError(f"unknown null family {family}")
        selected, summaries = evaluate_grid(
            conditions, include_all_metrics=False, **kwargs
        )
        selected_summary = next(row for row in summaries if row["name"] == selected.name)
        selected_names.append(selected.name)
        statistics[draw] = float(
            selected_summary["worst_fold_signed_linear_credit"]
        )
    return {
        "family": family,
        "draws": len(bank),
        "selected_candidate_names": selected_names,
        "selected_worst_fold_signed_linear_credit": statistics.tolist(),
        "p95_linear": float(np.percentile(statistics, 95)),
    }


def full_residual(
    selected_name: str,
    candidate_panel: np.ndarray,
    panel: Sequence[str],
    genes: Sequence[str],
) -> np.ndarray:
    names = [candidate.name for candidate in candidates()]
    if selected_name not in names:
        raise ContractError(f"unknown selected candidate {selected_name}")
    panel_index = {gene: index for index, gene in enumerate(panel)}
    gene_index = {gene: index for index, gene in enumerate(genes)}
    if any(gene not in gene_index for gene in panel):
        raise ContractError("canonical panel is not a subset of canonical gene axis")
    output = np.zeros((len(panel), len(genes)), dtype=np.float64)
    values = candidate_panel[names.index(selected_name)]
    for gene, panel_col in panel_index.items():
        output[:, gene_index[gene]] = values[:, panel_col]
    return output


def invariant_report(
    candidate_panel: np.ndarray,
    selected_full: np.ndarray,
    panel: Sequence[str],
    genes: Sequence[str],
) -> dict[str, Any]:
    gene_index = {gene: index for index, gene in enumerate(genes)}
    panel_cols = np.asarray([gene_index[gene] for gene in panel], dtype=np.int64)
    visible = np.ones(len(genes), dtype=bool)
    visible[panel_cols] = False
    diagonal = np.arange(len(panel))
    report = {
        "identity_exact_zero": bool(np.all(candidate_panel[0] == 0)),
        "all_candidate_self_coordinates_exact_zero": bool(
            np.all(np.diagonal(candidate_panel, axis1=1, axis2=2) == 0)
        ),
        "all_candidate_panel_row_sums_exact_zero": bool(
            all(rows_have_exact_zero_sum(arm) for arm in candidate_panel)
        ),
        "selected_pds_visible_coordinates_exact_zero": bool(
            np.all(selected_full[:, visible] == 0)
        ),
        "selected_target_self_coordinates_exact_zero": bool(
            np.all(selected_full[diagonal, panel_cols] == 0)
        ),
        "selected_panel_row_sums_exact_zero": bool(
            rows_have_exact_zero_sum(selected_full[:, panel_cols])
        ),
        "selected_full_finite": bool(np.all(np.isfinite(selected_full))),
    }
    report["all_pass"] = all(report.values())
    return report


def runtime_versions() -> dict[str, str]:
    import scipy

    return {
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }


def audit() -> None:
    authorization = validate_authorization()
    report = {
        "authorized": True,
        "scope": authorization["scope"],
        "allowed_input_count": len(authorization["allowed_input_files"]),
        "candidate_names": [candidate.name for candidate in candidates()],
        "null_draws_each": NULL_DRAWS,
        "output_exists": OUTPUT.exists(),
        "runtime_versions": runtime_versions(),
    }
    print(json.dumps(report, indent=2, sort_keys=True))


def build() -> None:
    authorization = validate_authorization()
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite {OUTPUT}")
    work = OUTPUT.parent
    work.mkdir(parents=True, exist_ok=True)
    stage = work / ".source-freeze-v1.staging"
    if stage.exists():
        raise FileExistsError(f"staging directory already exists: {stage}")
    stage.mkdir()
    try:
        panel = read_axis(REPO / "data/challenge/pert_counts.csv", "target_gene")
        genes = read_axis(REPO / "data/challenge/gene_names.csv", "gene_name")
        if len(panel) != 300 or len(genes) != 18533:
            raise ContractError("canonical challenge axis length mismatch")

        # Null assignments are frozen before any LFC, SE, FDR, or validity
        # effect array is loaded. Expression bins use only allowed base_mean.
        expression_bins = build_expression_bins(panel)
        target_permutations, response_permutations = freeze_permutations(
            expression_bins
        )
        atomic_npy(stage / "response_expression_bins.npy", expression_bins)
        atomic_npy(stage / "target_label_permutations.npy", target_permutations)
        atomic_npy(
            stage / "response_label_permutations.npy", response_permutations
        )

        conditions = [
            load_condition(name, panel, permit_effect_values=True)
            for name in CONDITION_NAMES
        ]
        observed, candidate_metrics = evaluate_grid(conditions)
        candidate_panel = build_full_candidates(conditions)
        target_null = null_selection(
            conditions, target_permutations, "target_label"
        )
        response_null = null_selection(
            conditions,
            response_permutations,
            "response_label_within_expression_bin",
        )
        observed_summary = next(
            row for row in candidate_metrics if row["name"] == observed.name
        )
        observed_worst = float(
            observed_summary["worst_fold_signed_linear_credit"]
        )
        beats_target = observed_worst > float(target_null["p95_linear"])
        beats_response = observed_worst > float(response_null["p95_linear"])
        survives = (
            observed.name != "identity"
            and bool(observed_summary["eligible"])
            and beats_target
            and beats_response
        )
        final_name = observed.name if survives else "identity"
        selected_full = full_residual(final_name, candidate_panel, panel, genes)
        invariants = invariant_report(candidate_panel, selected_full, panel, genes)
        if not invariants["all_pass"]:
            raise ContractError(f"frozen output invariant failure: {invariants}")

        atomic_npy(stage / "candidate_panel_effects.npy", candidate_panel)
        atomic_npy(stage / "selected_residual_full.npy", selected_full)
        selection = {
            "candidate_order": [candidate.name for candidate in candidates()],
            "observed_selected_before_null_gates": observed.name,
            "authorized_selected_after_all_source_gates": final_name,
            "nonidentity_survives_source_only_gates": survives,
            "selection_order": [
                "worst_fold_signed_linear_credit",
                "worst_fold_top20_weighted_sign_precision",
                "smaller_eta",
                "simpler_estimator",
            ],
            "candidate_metrics": candidate_metrics,
        }
        results = {
            "experiment": "118-pds-null-panel-residual",
            "stage": "source_only_first_stage",
            "observed": observed_summary,
            "nulls": {
                "target_label": target_null,
                "response_label_within_expression_bin": response_null,
            },
            "source_gates": {
                "observed_nonidentity": observed.name != "identity",
                "nonnegative_credit_and_precision_at_least_half_every_fold": bool(
                    observed_summary["eligible"] and observed.name != "identity"
                ),
                "observed_worst_credit_strictly_above_target_null_p95": beats_target,
                "observed_worst_credit_strictly_above_response_null_p95": beats_response,
                "all_invariants": invariants["all_pass"],
            },
            "nonidentity_survives_source_only_gates": survives,
            "final_selected_candidate": final_name,
            "invariants": invariants,
            "held_accessed": False,
            "counts_generated": False,
            "submission_built_or_sent": False,
        }
        atomic_json(stage / "selection.json", selection)
        atomic_json(stage / "results.json", results)

        artifact_names = [
            "candidate_panel_effects.npy",
            "response_expression_bins.npy",
            "response_label_permutations.npy",
            "results.json",
            "selected_residual_full.npy",
            "selection.json",
            "target_label_permutations.npy",
        ]
        freeze_manifest = {
            "schema_version": 1,
            "experiment": "118-pds-null-panel-residual",
            "stage": "source_only_first_stage",
            "authorization_sha256": sha256_file(AUTHORIZATION_PATH),
            "protocol_sha256": sha256_file(PROTOCOL_PATH),
            "source_files": authorization["allowed_input_files"],
            "implementation_files": authorization["implementation_files"],
            "runtime_versions": runtime_versions(),
            "artifacts": {
                name: {
                    "sha256": sha256_file(stage / name),
                    "bytes": (stage / name).stat().st_size,
                }
                for name in artifact_names
            },
            "logical_array_hashes": {
                "candidate_panel_effects.npy": sha256_array(candidate_panel),
                "selected_residual_full.npy": sha256_array(selected_full),
                "response_expression_bins.npy": sha256_array(expression_bins),
                "target_label_permutations.npy": sha256_array(target_permutations),
                "response_label_permutations.npy": sha256_array(
                    response_permutations
                ),
            },
            "invariants": invariants,
            "nulls_frozen_before_effect_array_load": True,
            "held_accessed": False,
            "counts_generated": False,
            "submission_built_or_sent": False,
        }
        atomic_json(stage / "freeze-manifest.json", freeze_manifest)
        stage.rename(OUTPUT)
    except BaseException:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    print(json.dumps(read_json(OUTPUT / "results.json"), indent=2, sort_keys=True))


def verify() -> None:
    authorization = validate_authorization()
    if not OUTPUT.is_dir():
        raise ContractError(f"missing frozen output {OUTPUT}")
    manifest = read_json(OUTPUT / "freeze-manifest.json")
    if manifest["authorization_sha256"] != sha256_file(AUTHORIZATION_PATH):
        raise ContractError("freeze authorization hash mismatch")
    if manifest["protocol_sha256"] != sha256_file(PROTOCOL_PATH):
        raise ContractError("freeze protocol hash mismatch")
    if manifest["source_files"] != authorization["allowed_input_files"]:
        raise ContractError("freeze source-file manifest mismatch")
    if manifest["implementation_files"] != authorization["implementation_files"]:
        raise ContractError("freeze implementation manifest mismatch")
    for name, record in manifest["artifacts"].items():
        path = OUTPUT / name
        if sha256_file(path) != record["sha256"] or path.stat().st_size != record["bytes"]:
            raise ContractError(f"frozen artifact mismatch: {name}")

    panel = read_axis(REPO / "data/challenge/pert_counts.csv", "target_gene")
    genes = read_axis(REPO / "data/challenge/gene_names.csv", "gene_name")
    candidate_panel = np.load(
        OUTPUT / "candidate_panel_effects.npy", allow_pickle=False
    )
    selected_full = np.load(OUTPUT / "selected_residual_full.npy", allow_pickle=False)
    expression_bins = np.load(
        OUTPUT / "response_expression_bins.npy", allow_pickle=False
    )
    target_permutations = np.load(
        OUTPUT / "target_label_permutations.npy", allow_pickle=False
    )
    response_permutations = np.load(
        OUTPUT / "response_label_permutations.npy", allow_pickle=False
    )
    if candidate_panel.shape != (13, 300, 300):
        raise ContractError("candidate panel array shape mismatch")
    if selected_full.shape != (300, 18533):
        raise ContractError("selected full residual shape mismatch")
    validate_permutations(target_permutations, response_permutations, expression_bins)
    logical = {
        "candidate_panel_effects.npy": sha256_array(candidate_panel),
        "selected_residual_full.npy": sha256_array(selected_full),
        "response_expression_bins.npy": sha256_array(expression_bins),
        "target_label_permutations.npy": sha256_array(target_permutations),
        "response_label_permutations.npy": sha256_array(response_permutations),
    }
    if logical != manifest["logical_array_hashes"]:
        raise ContractError("logical array hash mismatch")
    invariants = invariant_report(candidate_panel, selected_full, panel, genes)
    if not invariants["all_pass"] or invariants != manifest["invariants"]:
        raise ContractError("verified invariants do not match frozen invariants")
    selection = read_json(OUTPUT / "selection.json")
    expected_full = full_residual(
        selection["authorized_selected_after_all_source_gates"],
        candidate_panel,
        panel,
        genes,
    )
    if not np.array_equal(expected_full, selected_full):
        raise ContractError("selected full residual does not match selected panel arm")
    report = {
        "verified": True,
        "output": str(OUTPUT),
        "freeze_manifest_sha256": sha256_file(OUTPUT / "freeze-manifest.json"),
        "invariants": invariants,
        "nonidentity_survives_source_only_gates": read_json(
            OUTPUT / "results.json"
        )["nonidentity_survives_source_only_gates"],
    }
    print(json.dumps(report, indent=2, sort_keys=True))


def main(argv: Iterable[str] | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("audit", "build", "verify"))
    args = parser.parse_args(argv)
    if args.command == "audit":
        audit()
    elif args.command == "build":
        build()
    else:
        verify()


if __name__ == "__main__":
    main()

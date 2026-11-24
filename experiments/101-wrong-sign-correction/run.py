#!/usr/bin/env python3
"""Experiment 101: cross-fitted selective wrong-sign correction.

The outer truth is never read by ``freeze``.  ``score`` opens it only after the
real selection, real outer predictions, null selections, and null model
parameters have been frozen below this experiment directory.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.util
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
N_FOLDS = 4
N_NULL = 20
N_BOOT = 2_000
RIDGES = (0.1, 1.0, 10.0)
TAUS = (0.10, 0.25, 0.50)
LAMBDAS = (0.5, 1.0)
MIN_GATE = 10


def _load_exp100():
    path = ROOT / "experiments/100-sign-confidence/run.py"
    spec = importlib.util.spec_from_file_location("exp100_reuse", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import experiment 100 from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    # Experiment 100's reusable loaders cache lazily through these globals.
    # Redirect them before constructing any object so experiment 101 is the
    # only directory written by this run.
    module.WORK = WORK
    module.STATS = WORK / "donor-stats"
    module.STRING_CACHE = WORK / "string-pca.npz"
    return module


E100 = _load_exp100()
PANEL = E100.PANEL
PANEL_SET = E100.PANEL_SET
DONOR_POOL = E100.DONOR_POOL
SURFACES = E100.SURFACES


def stable_seed(*parts: object) -> int:
    return int.from_bytes(
        hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:4], "little"
    )


def target_fold(target: str) -> int:
    return stable_seed("exp101-fold", target) % N_FOLDS


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(map(str, values)).encode()).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def json_dump(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


@dataclass(frozen=True)
class Candidate:
    ridge: float
    tau: float
    lam: float

    @property
    def key(self) -> str:
        return f"flip:r={self.ridge:g}:tau={self.tau:.2f}:lambda={self.lam:g}"

    @property
    def zero_key(self) -> str:
        return f"zero:r={self.ridge:g}:tau={self.tau:.2f}"

    def as_dict(self) -> dict[str, float]:
        return {"ridge": self.ridge, "tau": self.tau, "lambda": self.lam}


CANDIDATES = tuple(
    Candidate(ridge, tau, lam) for ridge in RIDGES for tau in TAUS for lam in LAMBDAS
)


@dataclass
class SignedModel:
    mean: np.ndarray
    scale: np.ndarray
    coef: np.ndarray
    intercept: float
    ridge: float

    def predict(self, x: np.ndarray) -> np.ndarray:
        z = (np.asarray(x, dtype=np.float64) - self.mean) / self.scale
        return z @ self.coef + self.intercept


def fit_ridge(
    x: np.ndarray, y: np.ndarray, weight: np.ndarray, ridge: float
) -> SignedModel:
    y = np.asarray(y, dtype=np.float64)
    if not np.all(np.isin(np.unique(y), (-1.0, 1.0))):
        raise ValueError("ridge labels must be exactly {-1,+1}")
    w = np.asarray(weight, dtype=np.float64).copy()
    if not np.all(np.isfinite(w)) or w.sum() <= 0:
        raise ValueError("invalid training weights")
    w *= len(w) / w.sum()
    scaler = StandardScaler()
    scaler.fit(x, sample_weight=w)
    scale = np.asarray(scaler.scale_, dtype=np.float64)
    scale[scale == 0] = 1.0
    z = (np.asarray(x, dtype=np.float64) - scaler.mean_) / scale
    model = Ridge(alpha=float(ridge), fit_intercept=True, solver="cholesky")
    model.fit(z, y, sample_weight=w)
    return SignedModel(
        mean=np.asarray(scaler.mean_, dtype=np.float64),
        scale=scale,
        coef=np.asarray(model.coef_, dtype=np.float64),
        intercept=float(model.intercept_),
        ridge=float(ridge),
    )


@dataclass
class SignedBlock:
    target: str
    context_row: int
    gate_indices: np.ndarray
    features: np.ndarray
    y: np.ndarray
    weight: np.ndarray


@dataclass
class ContextData:
    name: str
    sources: tuple[str, ...]
    ref: object
    genes: np.ndarray
    targets: np.ndarray
    bctrl: np.ndarray
    ctrl_cpm: np.ndarray
    active: np.ndarray
    panel_mask: np.ndarray
    builder: object
    base_effect: np.ndarray
    base_dh: np.ndarray
    base_lfc: np.ndarray
    real_effect: np.ndarray
    score_mask: np.ndarray
    blocks: list[SignedBlock]
    dense_y: np.ndarray
    target_to_ref: np.ndarray
    used_output_rows: np.ndarray
    used_ref_rows: np.ndarray


def build_unwhitened_effect(
    genes: np.ndarray,
    targets: np.ndarray,
    sources: Sequence[str],
    donors: dict[str, object],
) -> np.ndarray:
    allowed = np.asarray([g not in PANEL_SET for g in genes], dtype=bool)
    rows = []
    for target in targets:
        vectors = [donors[name].vector(str(target), genes, allowed) for name in sources]
        values = [value for value in vectors if value is not None]
        rows.append(
            np.mean(values, axis=0)
            if values
            else np.zeros(len(genes), dtype=np.float64)
        )
    effect = np.vstack(rows).astype(np.float64)
    effect[:, ~allowed] = 0.0
    scale = (
        E100.ALPHA
        * E100.EFFECT_SCALE
        / (float(np.median(np.linalg.norm(effect[:, allowed], axis=1))) or 1.0)
    )
    effect *= scale
    effect[:, ~allowed] = 0.0
    return effect


def build_context(
    outer: str,
    context: str,
    string: object,
    donors: dict[str, object],
    max_rss_gib: float,
) -> ContextData:
    if context == outer:
        raise AssertionError("outer truth cannot be an inner context")
    ref = E100.RefData.load(context)
    sources = E100.source_names_for_training(outer, context)
    blocked = {
        SURFACES[outer]["truth_donor"],
        SURFACES[context]["truth_donor"],
    }
    if blocked & set(sources):
        raise AssertionError("outer or held-inner truth donor entered source set")
    genes, targets, bctrl, ctrl_cpm = E100.load_controls_and_targets(
        SURFACES[context]["harness"]
    )
    if not np.array_equal(genes, ref.genes):
        raise ValueError(f"{context}: harness/ref gene mismatch")
    rpos = {str(p): i for i, p in enumerate(ref.perts)}
    if any(str(t) not in rpos for t in targets):
        raise ValueError(f"{context}: harness target missing from reference cache")
    target_to_ref = np.asarray([rpos[str(t)] for t in targets], dtype=np.int64)
    panel_mask = np.asarray([g in PANEL_SET for g in genes], dtype=bool)
    active = (~panel_mask) & (ctrl_cpm > E100.MIN_CPM)
    builder = E100.FeatureBuilder(genes, ctrl_cpm, sources, string, donors)
    base_effect = build_unwhitened_effect(genes, targets, sources, donors)
    base_dh = E100.dh_from_effect(base_effect, bctrl, ctrl_cpm, panel_mask)
    base_lfc = E100.lfc_from_dh(base_dh, bctrl, ref.gate, ref.ref_lin)
    real_effect = E100.load_truth_effects(SURFACES[context]["harness"], targets, genes)
    gallery_excluded = PANEL_SET | set(map(str, targets))
    score_mask = np.asarray([g not in gallery_excluded for g in genes], dtype=bool)

    dense_y = np.zeros((len(ref.perts), len(ref.gate_idx)), dtype=np.int8)
    for out_row, ref_row in enumerate(target_to_ref):
        dense_y[ref_row] = (
            np.sign(ref.lfc[ref_row]) * np.sign(base_lfc[out_row])
        ).astype(np.int8)

    blocks: list[SignedBlock] = []
    used_output_rows: list[int] = []
    used_ref_rows: list[int] = []
    eligible_gate = np.asarray(
        [ref.genes[i] not in PANEL_SET for i in ref.gate_idx], dtype=bool
    )
    out_pos = {str(t): i for i, t in enumerate(targets)}
    for ref_row, target in enumerate(ref.perts):
        out_row = out_pos.get(str(target))
        if out_row is None:
            continue
        sig = ref.sig[ref_row]
        if len(sig) < MIN_GATE:
            continue
        keep = (
            eligible_gate[sig]
            & np.isfinite(ref.lfc[ref_row, sig])
            & (ref.lfc[ref_row, sig] != 0)
            & np.isfinite(base_lfc[out_row, sig])
            & (base_lfc[out_row, sig] != 0)
        )
        sig = sig[keep]
        if len(sig) < MIN_GATE:
            continue
        full_idx = ref.gate_idx[sig]
        features = builder.features(str(target), full_idx, 16)
        if features is None:
            continue
        _, x16 = features
        y = (np.sign(ref.lfc[ref_row, sig]) * np.sign(base_lfc[out_row, sig])).astype(
            np.int8
        )
        valid = np.isin(y, (-1, 1))
        sig, x16, y = sig[valid], x16[valid], y[valid]
        if len(sig) < MIN_GATE:
            continue
        denom = float(np.abs(ref.lfc[ref_row, sig]).sum()) or 1.0
        weight = np.abs(base_lfc[out_row, sig]) / denom
        if weight.sum() <= 0:
            continue
        blocks.append(
            SignedBlock(
                target=str(target),
                context_row=ref_row,
                gate_indices=sig,
                features=x16,
                y=y,
                weight=weight.astype(np.float64),
            )
        )
        used_output_rows.append(out_row)
        used_ref_rows.append(ref_row)
    E100.guard(max_rss_gib, f"build inner context {context}")
    return ContextData(
        name=context,
        sources=tuple(sources),
        ref=ref,
        genes=genes,
        targets=targets,
        bctrl=bctrl,
        ctrl_cpm=ctrl_cpm,
        active=active,
        panel_mask=panel_mask,
        builder=builder,
        base_effect=base_effect,
        base_dh=base_dh,
        base_lfc=base_lfc,
        real_effect=real_effect,
        score_mask=score_mask,
        blocks=blocks,
        dense_y=dense_y,
        target_to_ref=target_to_ref,
        used_output_rows=np.asarray(used_output_rows, dtype=np.int64),
        used_ref_rows=np.asarray(used_ref_rows, dtype=np.int64),
    )


def collect_training(
    contexts: Sequence[ContextData],
    held_targets: set[str],
    null_seed: int | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    ws: list[np.ndarray] = []
    for context in contexts:
        eligible = [b for b in context.blocks if b.target not in held_targets]
        null_map: dict[int, int] = {}
        if null_seed is not None:
            rows = np.asarray([b.context_row for b in eligible], dtype=np.int64)
            rng = np.random.default_rng(
                stable_seed("exp101-null", null_seed, context.name)
            )
            perm = rng.permutation(rows)
            null_map = {int(a): int(b) for a, b in zip(rows, perm)}
        for block in eligible:
            x = block.features
            y = block.y
            w = block.weight
            if null_seed is not None:
                source_row = null_map[block.context_row]
                y = context.dense_y[source_row, block.gate_indices]
                valid = np.isin(y, (-1, 1))
                x, y, w = x[valid], y[valid], w[valid]
            if len(y):
                xs.append(x)
                ys.append(y)
                ws.append(w)
    if not xs:
        raise ValueError("no training rows remain after whole-target blocking")
    x = np.vstack(xs)
    y = np.concatenate(ys)
    w = np.concatenate(ws)
    if not np.all(np.isin(np.unique(y), (-1, 1))):
        raise AssertionError("collected labels escaped {-1,+1}")
    return x, y, w


def fit_fold_models(
    contexts: Sequence[ContextData],
    eval_targets: np.ndarray,
    ridge: float,
    null_seed: int | None,
) -> list[SignedModel]:
    models = []
    for fold in range(N_FOLDS):
        held = {
            str(target) for target in eval_targets if target_fold(str(target)) == fold
        }
        x, y, w = collect_training(contexts, held, null_seed)
        models.append(fit_ridge(x, y, w, ridge))
    return models


def predict_scores(
    builder: object,
    targets: np.ndarray,
    active: np.ndarray,
    models: Sequence[SignedModel],
) -> np.ndarray:
    scores = np.ones((len(targets), len(active)), dtype=np.float32)
    idx = np.flatnonzero(active)
    for row, target in enumerate(targets):
        features = builder.features(str(target), idx, 16)
        if features is None:
            continue
        _, x16 = features
        scores[row, idx] = models[target_fold(str(target))].predict(x16)
    scores[:, ~active] = 1.0
    return scores


def transform_effect(
    effect: np.ndarray,
    score: np.ndarray,
    active: np.ndarray,
    tau: float,
    lam: float,
    mode: str,
) -> tuple[np.ndarray, np.ndarray]:
    selected = (score < -float(tau)) & active[None, :] & (effect != 0)
    out = np.asarray(effect, dtype=np.float64).copy()
    if mode == "flip":
        out[selected] *= -float(lam)
    elif mode == "zero":
        out[selected] = 0.0
    elif mode != "identity":
        raise ValueError(mode)
    out[:, ~active] = 0.0
    return out, selected


def per_target_credit(
    context: ContextData, lfc: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    credit = []
    numerator = []
    denominator = []
    for out_row, ref_row in zip(context.used_output_rows, context.used_ref_rows):
        sig = context.ref.sig[ref_row]
        sig = sig[
            np.asarray(
                [
                    context.ref.genes[context.ref.gate_idx[i]] not in PANEL_SET
                    for i in sig
                ],
                dtype=bool,
            )
            & np.isfinite(context.ref.lfc[ref_row, sig])
            & (context.ref.lfc[ref_row, sig] != 0)
            & np.isfinite(context.base_lfc[out_row, sig])
            & (context.base_lfc[out_row, sig] != 0)
        ]
        truth = context.ref.lfc[ref_row, sig]
        pred = lfc[out_row, sig]
        denom = float(np.abs(truth).sum()) or 1.0
        credit.append(float(np.sign(truth) @ pred) / denom)
        numerator.append(float(np.sign(truth) @ pred))
        denominator.append(denom)
    return (
        np.asarray(credit, dtype=np.float64),
        np.asarray(numerator, dtype=np.float64),
        np.asarray(denominator, dtype=np.float64),
    )


def tail_error_parts(
    context: ContextData, selected: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    wrong_weight = []
    selected_weight = []
    selected_count = []
    for out_row, ref_row in zip(context.used_output_rows, context.used_ref_rows):
        sig = context.ref.sig[ref_row]
        eligible = np.asarray(
            [context.ref.genes[context.ref.gate_idx[i]] not in PANEL_SET for i in sig],
            dtype=bool,
        )
        sig = sig[
            eligible
            & np.isfinite(context.ref.lfc[ref_row, sig])
            & (context.ref.lfc[ref_row, sig] != 0)
            & np.isfinite(context.base_lfc[out_row, sig])
            & (context.base_lfc[out_row, sig] != 0)
        ]
        full_idx = context.ref.gate_idx[sig]
        choose = selected[out_row, full_idx]
        truth = context.ref.lfc[ref_row, sig]
        base = context.base_lfc[out_row, sig]
        denom = float(np.abs(truth).sum()) or 1.0
        weight = np.abs(base) / denom
        wrong = np.sign(truth) != np.sign(base)
        wrong_weight.append(float(weight[choose & wrong].sum()))
        selected_weight.append(float(weight[choose].sum()))
        selected_count.append(int(choose.sum()))
    return (
        np.asarray(wrong_weight),
        np.asarray(selected_weight),
        np.asarray(selected_count),
    )


def evaluate_effect(
    context: ContextData,
    effect: np.ndarray,
    selected: np.ndarray | None = None,
) -> dict[str, object]:
    dh = E100.dh_from_effect(
        effect, context.bctrl, context.ctrl_cpm, context.panel_mask
    )
    lfc = E100.lfc_from_dh(dh, context.bctrl, context.ref.gate, context.ref.ref_lin)
    credit, _, _ = per_target_credit(context, lfc)
    pds = E100.pds_per_target(dh, context.real_effect, context.score_mask)
    result: dict[str, object] = {
        "credit": credit,
        "linear_nmae": 1.0 - credit,
        "pds": pds,
    }
    if selected is not None:
        wrong, total, count = tail_error_parts(context, selected)
        result.update(
            {
                "tail_wrong_weight": wrong,
                "tail_selected_weight": total,
                "tail_selected_count": count,
            }
        )
    return result


def choose_candidate(
    cells: dict[str, list[dict[str, object]]],
) -> tuple[Candidate | None, dict[str, object]]:
    summaries: dict[str, object] = {}
    eligible: list[Candidate] = []
    for candidate in CANDIDATES:
        values = cells[candidate.key]
        credit = [float(v["credit_delta"]) for v in values]
        pds = [float(v["pds_delta"]) for v in values]
        summaries[candidate.key] = {
            "inner_credit_deltas": credit,
            "inner_pds_deltas": pds,
            "macro_credit_delta": float(np.mean(credit)),
            "minimum_pds_delta": float(min(pds)),
        }
        if min(pds) >= 0.0:
            eligible.append(candidate)
    positive = [
        candidate
        for candidate in eligible
        if float(summaries[candidate.key]["macro_credit_delta"]) > 0.0
    ]
    if not positive:
        return None, {
            "selected": "identity",
            "eligible_count": len(eligible),
            "candidate_summaries": summaries,
            "rule": "identity because no PDS-safe candidate had positive macro inner credit",
        }
    best = max(
        float(summaries[candidate.key]["macro_credit_delta"]) for candidate in positive
    )
    near = [
        candidate
        for candidate in positive
        if abs(float(summaries[candidate.key]["macro_credit_delta"]) - best) <= 1e-12
    ]
    selected = min(
        near,
        key=lambda value: (-value.ridge, -value.tau, value.lam, value.key),
    )
    return selected, {
        "selected": selected.key,
        "selected_candidate": selected.as_dict(),
        "eligible_count": len(eligible),
        "best_macro_credit_delta": best,
        "candidate_summaries": summaries,
        "rule": (
            "maximize macro inner signed-credit gain among candidates with "
            "nonnegative complete-gallery PDS in both inner cells"
        ),
    }


def nested_select(
    outer: str,
    contexts: dict[str, ContextData],
    max_rss_gib: float,
    null_seed: int | None,
) -> tuple[Candidate | None, dict[str, object]]:
    cells: dict[str, list[dict[str, object]]] = {
        candidate.key: [] for candidate in CANDIDATES
    }
    cell_audit = []
    names = tuple(contexts)
    for held_name in names:
        fit_names = tuple(name for name in names if name != held_name)
        if len(fit_names) != 1:
            raise AssertionError("each inner cell must fit exactly one context")
        held = contexts[held_name]
        identity = evaluate_effect(held, held.base_effect)
        models_by_ridge = {
            ridge: fit_fold_models(
                [contexts[fit_names[0]]],
                held.targets,
                ridge,
                null_seed,
            )
            for ridge in RIDGES
        }
        scores_by_ridge = {
            ridge: predict_scores(held.builder, held.targets, held.active, models)
            for ridge, models in models_by_ridge.items()
        }
        for candidate in CANDIDATES:
            effect, _ = transform_effect(
                held.base_effect,
                scores_by_ridge[candidate.ridge],
                held.active,
                candidate.tau,
                candidate.lam,
                "flip",
            )
            result = evaluate_effect(held, effect)
            cells[candidate.key].append(
                {
                    "held_context": held_name,
                    "fit_context": fit_names[0],
                    "credit_delta": float(
                        np.mean(result["credit"] - identity["credit"])
                    ),
                    "pds_delta": float(np.mean(result["pds"] - identity["pds"])),
                }
            )
        cell_audit.append(
            {
                "held_context": held_name,
                "fit_context": fit_names[0],
                "held_truth_donor": SURFACES[held_name]["truth_donor"],
                "prediction_sources": list(held.sources),
                "outer_truth_donor_absent": (
                    SURFACES[outer]["truth_donor"] not in held.sources
                ),
                "whole_target_folds": N_FOLDS,
            }
        )
        del models_by_ridge, scores_by_ridge, identity
        gc.collect()
        E100.guard(max_rss_gib, f"inner selection {outer}/{held_name}")
    selected, audit = choose_candidate(cells)
    audit["outer"] = outer
    audit["null_seed"] = null_seed
    audit["cells"] = cell_audit
    return selected, audit


def model_arrays(models: Sequence[SignedModel]) -> dict[str, np.ndarray]:
    return {
        "means": np.vstack([m.mean for m in models]).astype(np.float32),
        "scales": np.vstack([m.scale for m in models]).astype(np.float32),
        "coefs": np.vstack([m.coef for m in models]).astype(np.float32),
        "intercepts": np.asarray([m.intercept for m in models]),
        "ridges": np.asarray([m.ridge for m in models]),
    }


def models_from_arrays(
    values: dict[str, np.ndarray], draw: int | None = None
) -> list[SignedModel]:
    prefix = () if draw is None else (draw,)
    means = values["means"][prefix]
    scales = values["scales"][prefix]
    coefs = values["coefs"][prefix]
    intercepts = values["intercepts"][prefix]
    ridges = values["ridges"][prefix]
    return [
        SignedModel(
            mean=means[fold].astype(np.float64),
            scale=scales[fold].astype(np.float64),
            coef=coefs[fold].astype(np.float64),
            intercept=float(intercepts[fold]),
            ridge=float(ridges[fold]),
        )
        for fold in range(N_FOLDS)
    ]


def freeze_surface(surface: str, max_rss_gib: float, n_null: int) -> None:
    outdir = WORK / surface
    outdir.mkdir(parents=True, exist_ok=True)
    frozen_path = outdir / "frozen_predictions.npz"
    if frozen_path.exists():
        raise FileExistsError(f"{frozen_path} exists; refusing to overwrite")
    string = E100.StringTable()
    donors = {name: E100.DonorLite(name) for name in DONOR_POOL}
    train_names = tuple(name for name in SURFACES if name != surface)
    contexts = {
        name: build_context(surface, name, string, donors, max_rss_gib)
        for name in train_names
    }

    selected, selection_audit = nested_select(
        surface, contexts, max_rss_gib, null_seed=None
    )
    genes, targets, bctrl, ctrl_cpm = E100.load_controls_and_targets(
        SURFACES[surface]["harness"]
    )
    panel_mask = np.asarray([g in PANEL_SET for g in genes], dtype=bool)
    active = (~panel_mask) & (ctrl_cpm > E100.MIN_CPM)
    final_builder = E100.FeatureBuilder(
        genes,
        ctrl_cpm,
        SURFACES[surface]["sources"],
        string,
        donors,
    )
    identity_effect, identity_audit = E100.build_identity_effect(
        genes,
        targets,
        SURFACES[surface]["sources"],
        donors,
        SURFACES[surface]["whitener"],
        max_rss_gib,
    )
    dh_identity = E100.dh_from_effect(identity_effect, bctrl, ctrl_cpm, panel_mask)

    if selected is None:
        score = np.ones_like(identity_effect, dtype=np.float32)
        flip_effect = identity_effect.copy()
        zero_effect = identity_effect.copy()
        selected_mask = np.zeros_like(identity_effect, dtype=bool)
        real_models: list[SignedModel] = []
    else:
        real_models = fit_fold_models(
            list(contexts.values()),
            targets,
            selected.ridge,
            null_seed=None,
        )
        score = predict_scores(final_builder, targets, active, real_models)
        flip_effect, selected_mask = transform_effect(
            identity_effect,
            score,
            active,
            selected.tau,
            selected.lam,
            "flip",
        )
        zero_effect, zero_mask = transform_effect(
            identity_effect,
            score,
            active,
            selected.tau,
            selected.lam,
            "zero",
        )
        if not np.array_equal(selected_mask, zero_mask):
            raise AssertionError("flip and zero controls selected different tails")
    dh_flip = E100.dh_from_effect(flip_effect, bctrl, ctrl_cpm, panel_mask)
    dh_zero = E100.dh_from_effect(zero_effect, bctrl, ctrl_cpm, panel_mask)

    np.savez_compressed(
        frozen_path,
        genes=genes,
        targets=targets,
        bctrl=bctrl.astype(np.float64),
        ctrl_cpm=ctrl_cpm.astype(np.float64),
        panel_mask=panel_mask,
        active=active,
        identity_effect=identity_effect.astype(np.float32),
        score=score.astype(np.float32),
        selected_mask=selected_mask,
        dh_identity=dh_identity.astype(np.float32),
        dh_flip=dh_flip.astype(np.float32),
        dh_zero=dh_zero.astype(np.float32),
    )
    if real_models:
        np.savez_compressed(outdir / "real_models.npz", **model_arrays(real_models))
    else:
        json_dump(outdir / "real_models.identity.json", {"selected": "identity"})

    null_selections: list[Candidate | None] = []
    null_audits = []
    null_model_sets: list[list[SignedModel] | None] = []
    for draw in range(n_null):
        null_selected, null_audit = nested_select(
            surface, contexts, max_rss_gib, null_seed=draw
        )
        null_selections.append(null_selected)
        null_audits.append(null_audit)
        if null_selected is None:
            null_model_sets.append(None)
        else:
            null_model_sets.append(
                fit_fold_models(
                    list(contexts.values()),
                    targets,
                    null_selected.ridge,
                    null_seed=draw,
                )
            )
        E100.guard(max_rss_gib, f"freeze null {draw}")

    feature_dim = (
        len(real_models[0].coef)
        if real_models
        else next(
            (len(models[0].coef) for models in null_model_sets if models is not None),
            int(contexts[train_names[0]].blocks[0].features.shape[1]),
        )
    )
    means = np.zeros((n_null, N_FOLDS, feature_dim), dtype=np.float32)
    scales = np.ones_like(means)
    coefs = np.zeros_like(means)
    intercepts = np.ones((n_null, N_FOLDS), dtype=np.float64)
    ridges = np.full((n_null, N_FOLDS), np.nan, dtype=np.float64)
    identity_draw = np.zeros(n_null, dtype=bool)
    candidate_values = np.full((n_null, 3), np.nan, dtype=np.float64)
    for draw, (candidate, models) in enumerate(zip(null_selections, null_model_sets)):
        if candidate is None or models is None:
            identity_draw[draw] = True
            continue
        arrays = model_arrays(models)
        means[draw] = arrays["means"]
        scales[draw] = arrays["scales"]
        coefs[draw] = arrays["coefs"]
        intercepts[draw] = arrays["intercepts"]
        ridges[draw] = arrays["ridges"]
        candidate_values[draw] = [
            candidate.ridge,
            candidate.tau,
            candidate.lam,
        ]
    np.savez_compressed(
        outdir / "null_models.npz",
        means=means,
        scales=scales,
        coefs=coefs,
        intercepts=intercepts,
        ridges=ridges,
        identity_draw=identity_draw,
        candidates=candidate_values,
    )
    json_dump(
        outdir / "selection.json",
        {
            "real": selection_audit,
            "null": null_audits,
        },
    )
    audit = {
        "surface": surface,
        "outer_truth_donor": SURFACES[surface]["truth_donor"],
        "final_sources": list(SURFACES[surface]["sources"]),
        "training_contexts": list(train_names),
        "training_source_sets": {
            name: list(contexts[name].sources) for name in train_names
        },
        "training_blocks": {name: len(contexts[name].blocks) for name in train_names},
        "training_examples": {
            name: int(sum(len(block.y) for block in contexts[name].blocks))
            for name in train_names
        },
        "selected": "identity" if selected is None else selected.as_dict(),
        "target_folds": {str(t): target_fold(str(t)) for t in targets},
        "fold_hashes": {
            str(fold): sha256_lines(
                sorted(str(t) for t in targets if target_fold(str(t)) == fold)
            )
            for fold in range(N_FOLDS)
        },
        "identity": identity_audit,
        "panel_coordinates_removed": int(panel_mask.sum()),
        "active_response_coordinates": int(active.sum()),
        "n_null": n_null,
        "inner_base": "unwhitened_median_norm_scaled",
        "outer_base": "experiment100_surface_whitened",
        "peak_rss_gib": E100.peak_rss_gib(),
        "frozen_sha256": sha256_file(frozen_path),
        "selection_sha256": sha256_file(outdir / "selection.json"),
        "null_models_sha256": sha256_file(outdir / "null_models.npz"),
    }
    json_dump(outdir / "freeze-audit.json", audit)
    print(json.dumps(audit, indent=2))


def ratio_ci(
    numerator: np.ndarray,
    denominator: np.ndarray,
    n_boot: int,
    seed: int,
) -> dict[str, float | None]:
    numerator = np.asarray(numerator, dtype=np.float64)
    denominator = np.asarray(denominator, dtype=np.float64)
    total = float(denominator.sum())
    estimate = float(numerator.sum() / total) if total > 0 else None
    if estimate is None:
        return {"estimate": None, "lo": None, "hi": None}
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(n_boot):
        idx = rng.integers(0, len(numerator), len(numerator))
        d = float(denominator[idx].sum())
        if d > 0:
            values.append(float(numerator[idx].sum() / d))
    return {
        "estimate": estimate,
        "lo": float(np.percentile(values, 2.5)),
        "hi": float(np.percentile(values, 97.5)),
    }


def score_surface(surface: str, max_rss_gib: float, n_boot: int) -> None:
    outdir = WORK / surface
    metrics_path = outdir / "metrics.json"
    if metrics_path.exists():
        raise FileExistsError(f"{metrics_path} exists; refusing to overwrite")
    frozen = np.load(outdir / "frozen_predictions.npz", allow_pickle=False)
    audit = json.loads((outdir / "freeze-audit.json").read_text())
    genes = frozen["genes"].astype(str)
    targets = frozen["targets"].astype(str)
    bctrl = frozen["bctrl"].astype(np.float64)
    ctrl_cpm = frozen["ctrl_cpm"].astype(np.float64)
    panel_mask = frozen["panel_mask"].astype(bool)
    active = frozen["active"].astype(bool)
    selected_mask = frozen["selected_mask"].astype(bool)
    identity_effect = frozen["identity_effect"].astype(np.float64)
    dh = {
        "identity": frozen["dh_identity"].astype(np.float64),
        "flip": frozen["dh_flip"].astype(np.float64),
        "zero": frozen["dh_zero"].astype(np.float64),
    }

    # First held-surface truth read in this command.
    ref = E100.RefData.load(surface)
    if not np.array_equal(ref.genes, genes):
        raise ValueError("held reference gene order changed")
    rpos = {str(p): i for i, p in enumerate(ref.perts)}
    target_to_ref = np.asarray([rpos[str(t)] for t in targets], dtype=np.int64)
    lfc = {
        name: E100.lfc_from_dh(value, bctrl, ref.gate, ref.ref_lin)
        for name, value in dh.items()
    }
    eligible_gate = np.asarray(
        [ref.genes[i] not in PANEL_SET for i in ref.gate_idx], dtype=bool
    )
    used_rows = []
    used_ref_rows = []
    per_target: dict[str, list[float]] = {
        "identity_credit": [],
        "flip_credit": [],
        "zero_credit": [],
        "tail_wrong_weight": [],
        "tail_selected_weight": [],
        "tail_selected_count": [],
    }
    for out_row, ref_row in enumerate(target_to_ref):
        sig = ref.sig[ref_row]
        sig = sig[
            eligible_gate[sig]
            & np.isfinite(ref.lfc[ref_row, sig])
            & (ref.lfc[ref_row, sig] != 0)
            & np.isfinite(lfc["identity"][out_row, sig])
            & (lfc["identity"][out_row, sig] != 0)
        ]
        if len(sig) < MIN_GATE:
            continue
        truth = ref.lfc[ref_row, sig]
        base = lfc["identity"][out_row, sig]
        denom = float(np.abs(truth).sum()) or 1.0
        for name in ("identity", "flip", "zero"):
            per_target[f"{name}_credit"].append(
                float(np.sign(truth) @ lfc[name][out_row, sig]) / denom
            )
        full_idx = ref.gate_idx[sig]
        choose = selected_mask[out_row, full_idx]
        weight = np.abs(base) / denom
        wrong = np.sign(truth) != np.sign(base)
        per_target["tail_wrong_weight"].append(float(weight[choose & wrong].sum()))
        per_target["tail_selected_weight"].append(float(weight[choose].sum()))
        per_target["tail_selected_count"].append(int(choose.sum()))
        used_rows.append(out_row)
        used_ref_rows.append(ref_row)
    per = {name: np.asarray(value) for name, value in per_target.items()}
    if not len(used_rows):
        raise ValueError("no outer targets passed the DE gate")

    real_effect = E100.load_truth_effects(SURFACES[surface]["harness"], targets, genes)
    score_mask = np.asarray(
        [gene not in (PANEL_SET | set(map(str, targets))) for gene in genes],
        dtype=bool,
    )
    pds_rows = {
        name: E100.pds_per_target(value, real_effect, score_mask)
        for name, value in dh.items()
    }
    E100.guard(max_rss_gib, "real outer scoring")

    null_file = np.load(outdir / "null_models.npz", allow_pickle=False)
    null_values = {name: null_file[name] for name in null_file.files}
    string = E100.StringTable()
    donors = {name: E100.DonorLite(name) for name in DONOR_POOL}
    builder = E100.FeatureBuilder(
        genes, ctrl_cpm, SURFACES[surface]["sources"], string, donors
    )
    identity_credit = per["identity_credit"]
    identity_pds = pds_rows["identity"]
    null_credit_delta = np.zeros(len(null_values["identity_draw"]))
    null_pds_delta = np.zeros_like(null_credit_delta)
    null_selected = []
    for draw in range(len(null_credit_delta)):
        if bool(null_values["identity_draw"][draw]):
            null_selected.append("identity")
            continue
        candidate_raw = null_values["candidates"][draw]
        candidate = Candidate(*map(float, candidate_raw))
        models = models_from_arrays(null_values, draw)
        score = predict_scores(builder, targets, active, models)
        effect, _ = transform_effect(
            identity_effect,
            score,
            active,
            candidate.tau,
            candidate.lam,
            "flip",
        )
        null_dh = E100.dh_from_effect(effect, bctrl, ctrl_cpm, panel_mask)
        null_lfc = E100.lfc_from_dh(null_dh, bctrl, ref.gate, ref.ref_lin)
        credits = []
        for out_row, ref_row in zip(used_rows, used_ref_rows):
            sig = ref.sig[ref_row]
            sig = sig[
                eligible_gate[sig]
                & np.isfinite(ref.lfc[ref_row, sig])
                & (ref.lfc[ref_row, sig] != 0)
                & np.isfinite(lfc["identity"][out_row, sig])
                & (lfc["identity"][out_row, sig] != 0)
            ]
            truth = ref.lfc[ref_row, sig]
            denom = float(np.abs(truth).sum()) or 1.0
            credits.append(float(np.sign(truth) @ null_lfc[out_row, sig]) / denom)
        null_credit_delta[draw] = float(np.mean(credits)) - float(
            identity_credit.mean()
        )
        null_pds = E100.pds_per_target(null_dh, real_effect, score_mask)
        null_pds_delta[draw] = float(np.mean(null_pds - identity_pds))
        null_selected.append(candidate.key)
        E100.guard(max_rss_gib, f"null outer score {draw}")

    comparisons = {
        "flip_credit_minus_identity": per["flip_credit"] - per["identity_credit"],
        "zero_credit_minus_identity": per["zero_credit"] - per["identity_credit"],
        "flip_credit_minus_zero": per["flip_credit"] - per["zero_credit"],
        "flip_pds_minus_identity": pds_rows["flip"] - pds_rows["identity"],
        "zero_pds_minus_identity": pds_rows["zero"] - pds_rows["identity"],
    }
    bootstrap = {
        name: E100.ci(values, n_boot, stable_seed(surface, name))
        for name, values in comparisons.items()
    }
    tail = ratio_ci(
        per["tail_wrong_weight"],
        per["tail_selected_weight"],
        n_boot,
        stable_seed(surface, "tail-error"),
    )
    metrics = {
        "surface": surface,
        "selected": audit["selected"],
        "targets_total": int(len(targets)),
        "targets_de_scored": int(len(used_rows)),
        "linear_nmae_proxy": {
            name: float(1.0 - per[f"{name}_credit"].mean())
            for name in ("identity", "flip", "zero")
        },
        "linear_credit_gain": {
            "flip_minus_identity": float(
                np.mean(comparisons["flip_credit_minus_identity"])
            ),
            "zero_minus_identity": float(
                np.mean(comparisons["zero_credit_minus_identity"])
            ),
            "flip_minus_zero": float(np.mean(comparisons["flip_credit_minus_zero"])),
        },
        "pds": {name: float(np.mean(value)) for name, value in pds_rows.items()},
        "pds_delta": {
            "flip_minus_identity": float(
                np.mean(comparisons["flip_pds_minus_identity"])
            ),
            "zero_minus_identity": float(
                np.mean(comparisons["zero_pds_minus_identity"])
            ),
        },
        "flipped_tail": {
            "weighted_error_rate": tail,
            "selected_pairs": int(per["tail_selected_count"].sum()),
            "targets_with_selected_pairs": int(np.sum(per["tail_selected_count"] > 0)),
            "total_selected_weight": float(per["tail_selected_weight"].sum()),
        },
        "bootstrap": bootstrap,
        "permutation_null": {
            "linear_credit_gain": null_credit_delta.tolist(),
            "pds_delta": null_pds_delta.tolist(),
            "selected": null_selected,
        },
        "peak_rss_gib": E100.peak_rss_gib(),
    }
    np.savez_compressed(
        outdir / "per_target_metrics.npz",
        targets=targets[np.asarray(used_rows, dtype=np.int64)],
        **per,
        **{f"pds_{name}": values for name, values in pds_rows.items()},
        null_credit_delta=null_credit_delta,
        null_pds_delta=null_pds_delta,
    )
    json_dump(metrics_path, metrics)
    print(json.dumps(metrics, indent=2))


def prepare(max_rss_gib: float) -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    string = E100.StringTable()
    donors = {}
    for name in DONOR_POOL:
        donors[name] = E100.DonorLite(name)
        E100.guard(max_rss_gib, f"prepare donor {name}")
    audit = {
        "string_rows": int(len(string.genes)),
        "string_available": int(string.available.sum()),
        "string_rank": int(string.pcs.shape[1]),
        "feature_rank": 16,
        "ridge_grid": list(RIDGES),
        "tau_grid": list(TAUS),
        "lambda_grid": list(LAMBDAS),
        "candidate_count": len(CANDIDATES),
        "n_folds": N_FOLDS,
        "n_null": N_NULL,
        "n_boot": N_BOOT,
        "panel_genes": len(PANEL),
        "panel_sha256": sha256_lines(PANEL),
        "donors": {
            name: {
                "perts": len(donor.perts),
                "genes": len(donor.genes),
                "mean_sha256": E100.sha256_array(donor.mean),
            }
            for name, donor in donors.items()
        },
        "peak_rss_gib": E100.peak_rss_gib(),
    }
    json_dump(WORK / "prepare-audit.json", audit)
    print(json.dumps(audit, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("prepare", "freeze", "score"))
    parser.add_argument("--surface", choices=tuple(SURFACES))
    parser.add_argument("--max-rss-gib", type=float, default=28.0)
    parser.add_argument("--n-null", type=int, default=N_NULL)
    parser.add_argument("--n-boot", type=int, default=N_BOOT)
    args = parser.parse_args()
    if tuple(sorted(RIDGES)) != (0.1, 1.0, 10.0):
        raise SystemExit("ridge grid changed from preregistration")
    if tuple(sorted(TAUS)) != (0.10, 0.25, 0.50):
        raise SystemExit("tau grid changed from preregistration")
    if tuple(sorted(LAMBDAS)) != (0.5, 1.0):
        raise SystemExit("lambda grid changed from preregistration")
    if args.command == "prepare":
        prepare(args.max_rss_gib)
        return
    if args.surface is None:
        parser.error("--surface is required for freeze/score")
    if args.command == "freeze":
        freeze_surface(args.surface, args.max_rss_gib, args.n_null)
    else:
        score_surface(args.surface, args.max_rss_gib, args.n_boot)


if __name__ == "__main__":
    main()

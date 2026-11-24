#!/usr/bin/env python3
"""Experiment 104: context-held-out, target-visible direct-credit adapter.

This module never imports, enumerates, or opens KOLF. All generated artifacts
are constrained to this experiment directory.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import resource
import sys
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "4")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "4")

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work-v2"
EXPERIMENT = "104-context-heldout-credit-adapter"

MAX_RSS_GIB = 24.0
ALPHA = 0.982
CD4_FRACTION = 0.09472706196346517
BETA = 0.10
RIDGE_ALPHA = 1.0
N_FOLDS = 4
N_NULL = 20
N_BOOT = 2_000
MACRO_BOOT = 20_000
MIN_GATE = 10
TS = 50_000.0

PANEL = tuple(
    pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str)
)
PANEL_SET = set(PANEL)
OFFICIAL = tuple(
    pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str)
)
OFFICIAL_SET = set(OFFICIAL)

REAL_DONOR_PATHS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "h1_2025": ROOT / "data/donors/h1_2025",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
}

SURFACES = {
    "h1": {
        "truth_donor": "h1_2025",
        "harness": ROOT / "data/harness/h1-test.h5ad",
        "refde": ROOT / "experiments/098-nmae/work/refde_h1.npz",
        "sources": ("k562_gwps", "orion_hct116", "orion_hek293t"),
        "operator": ROOT / "data/donors/whitener_union_pruned.npz",
    },
    "hct116": {
        "truth_donor": "orion_hct116",
        "harness": ROOT / "data/harness/orion-hct116-thinned.h5ad",
        "refde": ROOT / "experiments/098-nmae/work/refde_orion.npz",
        "sources": ("k562_gwps", "orion_hek293t"),
        "operator": (
            ROOT / "experiments/062-cd4-post-whiten/cloud-output/whiteners/"
            "held_orion_hct116.npz"
        ),
    },
    "hek293t": {
        "truth_donor": "orion_hek293t",
        "harness": ROOT / "data/harness/orion-hek293t-thinned.h5ad",
        "refde": ROOT / "experiments/099-donors/work/refde_hek293t.npz",
        "sources": ("k562_gwps", "orion_hct116"),
        "operator": (
            ROOT / "experiments/062-cd4-post-whiten/cloud-output/whiteners/"
            "held_orion_hek293t.npz"
        ),
    },
}

PM_PATH = ROOT / "data/external/perturbmatch"
CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
STRING_PATH = ROOT / "data/external/slim/gene_string_embeddings.v0.3.h5"


def _load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


E100 = _load_module(
    "exp104_reuse_exp100", ROOT / "experiments/100-sign-confidence/run.py"
)
E062 = _load_module(
    "exp104_reuse_exp062", ROOT / "experiments/062-cd4-post-whiten/run.py"
)
E064 = _load_module(
    "exp104_reuse_exp064", ROOT / "experiments/064-pm-cd4-composition/run.py"
)

# Redirect every lazy cache used by imported helpers into the sole write scope
# and remove the forbidden donor from the reusable module's registry.
E100.WORK = WORK
E100.STATS = WORK / "donor-stats"
E100.STRING_CACHE = WORK / "string-pca.npz"
E100.DONOR_PATHS = dict(REAL_DONOR_PATHS)
E100.DONOR_POOL = tuple(REAL_DONOR_PATHS)
E100.SURFACES = SURFACES
E100.ALPHA = ALPHA


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def guard(stage: str, limit_gib: float) -> None:
    peak = peak_rss_gib()
    print(f"[rss] {stage}: {peak:.3f} GiB", flush=True)
    if not 0 < limit_gib <= MAX_RSS_GIB:
        raise ValueError(f"memory limit must be in (0, {MAX_RSS_GIB}]")
    if peak >= limit_gib:
        raise MemoryError(f"{stage}: {peak:.3f} GiB >= {limit_gib:.3f} GiB")


def stable_seed(*parts: object) -> int:
    return int.from_bytes(
        hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:8], "little"
    )


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_array(values: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(values).tobytes()).hexdigest()


def atomic_json(path: Path, value: object) -> None:
    if HERE not in path.resolve().parents:
        raise ValueError(f"refusing write outside experiment: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def atomic_npz(path: Path, **values: object) -> None:
    if HERE not in path.resolve().parents:
        raise ValueError(f"refusing write outside experiment: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with tmp.open("wb") as handle:
        np.savez_compressed(handle, **values)
    tmp.replace(path)


@dataclass
class CreditModel:
    mean: np.ndarray
    scale: np.ndarray
    coef: np.ndarray
    intercept: float

    def predict(self, x: np.ndarray) -> np.ndarray:
        z = (np.asarray(x, dtype=np.float64) - self.mean) / self.scale
        return z @ self.coef + self.intercept


def fit_credit_model(
    x: np.ndarray, y: np.ndarray, weight: np.ndarray, alpha: float = RIDGE_ALPHA
) -> CreditModel:
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    weight = np.asarray(weight, dtype=np.float64)
    if x.ndim != 2 or y.shape != (len(x),) or weight.shape != (len(x),):
        raise ValueError("training arrays do not align")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("non-finite training values")
    if set(np.unique(y)) - {-1.0, 1.0}:
        raise ValueError("credit labels must be -1/+1")
    if np.any(weight < 0) or float(weight.sum()) <= 0:
        raise ValueError("invalid sample weights")
    weight = weight * (len(weight) / weight.sum())
    # Compute weighted moments directly. StandardScaler's incremental weighted
    # variance can round a nearly constant column slightly below zero.
    mean = np.average(x, axis=0, weights=weight)
    variance = np.average((x - mean) ** 2, axis=0, weights=weight)
    scale = np.sqrt(np.maximum(variance, 0.0))
    scale[~np.isfinite(scale) | (scale < 1e-12)] = 1.0
    z = (x - mean) / scale
    model = Ridge(alpha=float(alpha), fit_intercept=True)
    model.fit(z, y, sample_weight=weight)
    return CreditModel(
        mean=np.asarray(mean, dtype=np.float64),
        scale=scale,
        coef=np.asarray(model.coef_, dtype=np.float64),
        intercept=float(model.intercept_),
    )


def serialize_models(models: Sequence[CreditModel]) -> dict[str, np.ndarray]:
    return {
        "means": np.stack([model.mean for model in models]).astype(np.float32),
        "scales": np.stack([model.scale for model in models]).astype(np.float32),
        "coefs": np.stack([model.coef for model in models]).astype(np.float32),
        "intercepts": np.asarray(
            [model.intercept for model in models], dtype=np.float64
        ),
    }


def deserialize_models(z: np.lib.npyio.NpzFile, prefix: str) -> list[CreditModel]:
    means = z[f"{prefix}_means"]
    scales = z[f"{prefix}_scales"]
    coefs = z[f"{prefix}_coefs"]
    intercepts = z[f"{prefix}_intercepts"]
    return [
        CreditModel(
            mean=means[i].astype(np.float64),
            scale=scales[i].astype(np.float64),
            coef=coefs[i].astype(np.float64),
            intercept=float(intercepts[i]),
        )
        for i in range(len(intercepts))
    ]


def target_fold(target: str) -> int:
    return stable_seed(EXPERIMENT, "target-fold", target) % N_FOLDS


def directional_mask(genes: Sequence[str], target: str) -> np.ndarray:
    """Directional metrics exclude only the perturbed gene's own coordinate."""
    return np.asarray(
        [str(gene) in OFFICIAL_SET and str(gene) != str(target) for gene in genes],
        dtype=bool,
    )


def pds_mask(genes: Sequence[str]) -> np.ndarray:
    """PDS excludes the complete fixed 300-gene VCC panel."""
    return np.asarray(
        [str(gene) in OFFICIAL_SET and str(gene) not in PANEL_SET for gene in genes],
        dtype=bool,
    )


def robust_score_multiplier(
    score: np.ndarray,
    effect: np.ndarray,
    active: np.ndarray,
    beta: float = BETA,
) -> np.ndarray:
    score = np.asarray(score, dtype=np.float64)
    effect = np.asarray(effect, dtype=np.float64)
    active = np.asarray(active, dtype=bool)
    use = active & np.isfinite(score) & (effect != 0)
    out = np.ones_like(effect, dtype=np.float64)
    if int(use.sum()) < 2 or beta == 0:
        return out
    values = score[use]
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)))
    robust_scale = max(1.4826 * mad, 1e-6)
    z = np.clip((values - median) / robust_scale, -2.0, 2.0)
    out[use] = np.exp(float(beta) * z)
    return out


def transform_effect(
    effect: np.ndarray,
    score: np.ndarray,
    active: np.ndarray,
    supported: np.ndarray,
    beta: float = BETA,
) -> np.ndarray:
    base = np.asarray(effect, dtype=np.float64)
    score = np.asarray(score, dtype=np.float64)
    active = np.asarray(active, dtype=bool)
    supported = np.asarray(supported, dtype=bool)
    if base.shape != score.shape or base.shape != supported.shape:
        raise ValueError("effect, score, and support must align")
    out = base.copy()
    for row in range(len(out)):
        row_active = active & supported[row]
        before = float(np.linalg.norm(out[row, active]))
        multiplier = robust_score_multiplier(
            score[row], out[row], row_active, beta=beta
        )
        out[row] *= multiplier
        after = float(np.linalg.norm(out[row, active]))
        if before > 0 and after > 0:
            out[row] *= before / after
    out[:, ~active] = 0.0
    if np.any(np.sign(out[base != 0]) != np.sign(base[base != 0])):
        raise AssertionError("sign-preserving transform flipped a sign")
    if np.any(out[(base == 0)] != 0):
        raise AssertionError("transform introduced a nonzero coordinate")
    return out


def model_sources(outer: str, truth_context: str | None = None) -> tuple[str, ...]:
    blocked = {SURFACES[outer]["truth_donor"]}
    if truth_context is not None:
        blocked.add(SURFACES[truth_context]["truth_donor"])
    return tuple(name for name in REAL_DONOR_PATHS if name not in blocked)


@dataclass
class TrainingBlock:
    context: str
    target: str
    x_agree: np.ndarray
    x_full: np.ndarray
    y: np.ndarray
    weight: np.ndarray


def _control_log_profile(ctrl_cpm: np.ndarray) -> np.ndarray:
    return np.log1p(TS * np.asarray(ctrl_cpm, dtype=np.float64) / 1e6)


def build_training_blocks(
    outer: str,
    contexts: Sequence[str],
    string: object,
    donors: dict[str, object],
    limit_gib: float,
) -> list[TrainingBlock]:
    blocks: list[TrainingBlock] = []
    for context in contexts:
        ref = E100.RefData.load(context)
        sources = model_sources(outer, context)
        if not sources:
            raise ValueError(f"{outer}/{context}: empty source set")
        builder = E100.FeatureBuilder(ref.genes, ref.ctrl_cpm, sources, string, donors)
        bctrl = _control_log_profile(ref.ctrl_cpm)
        panel_mask = np.asarray([gene in PANEL_SET for gene in ref.genes], dtype=bool)
        rows: list[np.ndarray] = []
        for target in ref.perts:
            state = builder.target_state(str(target))
            rows.append(
                np.zeros(len(ref.genes), dtype=np.float64)
                if state is None
                else state.base
            )
        effect = np.vstack(rows)
        active = pds_mask(ref.genes)
        median_norm = float(np.median(np.linalg.norm(effect[:, active], axis=1))) or 1.0
        effect *= ALPHA * E100.EFFECT_SCALE / median_norm
        effect[:, panel_mask] = 0.0
        dh = E100.dh_from_effect(effect, bctrl, ref.ctrl_cpm, panel_mask)
        pred_lfc = E100.lfc_from_dh(dh, bctrl, ref.gate, ref.ref_lin)

        for row, target in enumerate(ref.perts):
            sig = np.asarray(ref.sig[row], dtype=np.int64)
            if len(sig) < MIN_GATE:
                continue
            full_idx = ref.gate_idx[sig]
            own = ~directional_mask(ref.genes[full_idx], str(target))
            truth = ref.lfc[row, sig]
            base = pred_lfc[row, sig]
            keep = (
                ~own
                & np.isfinite(truth)
                & np.isfinite(base)
                & (truth != 0)
                & (base != 0)
            )
            sig, full_idx, truth, base = (
                sig[keep],
                full_idx[keep],
                truth[keep],
                base[keep],
            )
            if len(sig) < MIN_GATE:
                continue
            features = builder.features(str(target), full_idx, 8)
            if features is None:
                continue
            x_agree, x_full = features
            y = np.sign(truth) * np.sign(base)
            denom = float(np.abs(truth).sum()) or 1.0
            weight = np.abs(base) / denom
            if float(weight.sum()) <= 0:
                continue
            blocks.append(
                TrainingBlock(
                    context=context,
                    target=str(target),
                    x_agree=x_agree,
                    x_full=x_full,
                    y=y.astype(np.float64),
                    weight=weight.astype(np.float64),
                )
            )
        guard(f"training blocks {outer}/{context}", limit_gib)
    if not blocks:
        raise ValueError(f"{outer}: no training blocks")
    return blocks


def collect_blocks(
    blocks: Sequence[TrainingBlock],
    kind: str,
    *,
    held_fold: int | None = None,
    null_draw: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    chosen = [
        block
        for block in blocks
        if held_fold is None or target_fold(block.target) != held_fold
    ]
    if not chosen:
        raise ValueError("no blocks after fold exclusion")
    label_blocks = list(chosen)
    if null_draw is not None:
        remapped: list[TrainingBlock] = []
        for block in chosen:
            rng = np.random.default_rng(
                stable_seed(EXPERIMENT, "null", null_draw, block.context, block.target)
            )
            y = block.y.copy()
            edges = np.quantile(block.weight, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
            strata = np.digitize(block.weight, edges[1:-1], right=True)
            for stratum in range(5):
                rows = np.flatnonzero(strata == stratum)
                if len(rows) > 1:
                    y[rows] = y[rng.permutation(rows)]
            remapped.append(
                TrainingBlock(
                    context=block.context,
                    target=block.target,
                    x_agree=block.x_agree,
                    x_full=block.x_full,
                    y=y,
                    weight=block.weight,
                )
            )
        label_blocks = remapped
    x = np.vstack(
        [
            block.x_agree if kind == "agreement" else block.x_full
            for block in label_blocks
        ]
    )
    y = np.concatenate([block.y for block in label_blocks])
    w = np.concatenate([block.weight for block in label_blocks])
    return x, y, w


def fit_models(
    blocks: Sequence[TrainingBlock], null_draw: int | None = None
) -> dict[str, CreditModel]:
    output = {}
    for kind in ("agreement", "full"):
        x, y, weight = collect_blocks(blocks, kind, null_draw=null_draw)
        output[kind] = fit_credit_model(x, y, weight)
    return output


def crossfit_diagnostics(blocks: Sequence[TrainingBlock]) -> dict[str, object]:
    rows = []
    for fold in range(N_FOLDS):
        models = {}
        for kind in ("agreement", "full"):
            x, y, w = collect_blocks(blocks, kind, held_fold=fold)
            models[kind] = fit_credit_model(x, y, w)
        held = [block for block in blocks if target_fold(block.target) == fold]
        for kind in ("agreement", "full"):
            numerator = 0.0
            denominator = 0.0
            for block in held:
                x = block.x_agree if kind == "agreement" else block.x_full
                score = models[kind].predict(x)
                multiplier = robust_score_multiplier(
                    score, np.sign(block.y), np.ones(len(score), dtype=bool)
                )
                numerator += float(np.sum(block.weight * block.y * multiplier))
                denominator += float(np.sum(block.weight))
            rows.append(
                {
                    "fold": fold,
                    "kind": kind,
                    "targets": len(held),
                    "weighted_credit": numerator / max(denominator, 1e-12),
                }
            )
    return {"folds": rows}


def _operator(surface: str, banks: Sequence[object], official: Sequence[str]):
    path = Path(SURFACES[surface]["operator"])
    if surface == "h1":
        metadata = json.loads(path.with_suffix(".json").read_text())
        expected_donors = [bank.name for bank in banks]
        if metadata.get("donors") != expected_donors:
            raise ValueError("production-pruned operator donor set changed")
        wcols = E064.operator_columns(banks, official, PANEL_SET)
        if int(metadata.get("n_cols", -1)) != len(wcols):
            raise ValueError("production-pruned operator axis changed")
        if metadata.get("basis_prune") != "essential":
            raise ValueError("production-pruned operator policy changed")
        if metadata.get("sha256") != sha256_file(path):
            raise ValueError("production-pruned operator hash changed")
        with np.load(path) as archive:
            evals = np.asarray(archive["evals"], dtype=np.float64)
            evecs = np.asarray(archive["evecs"], dtype=np.float32)
        if evecs.shape != (len(wcols), len(wcols)):
            raise ValueError("production-pruned operator matrix shape changed")
        return (
            E062.Operator(evals, evecs, metadata, path, path.with_suffix(".json")),
            wcols,
        )
    return E064.load_frozen_operator(
        path,
        banks,
        official,
        PANEL_SET,
        expected_sha256=sha256_file(path),
    )


def _partial_pm_slot(k562, pm, shared_columns: np.ndarray):
    effects = np.asarray(k562.effects, dtype=np.float64).copy()
    both = k562.target_available & pm.target_available
    rows = np.flatnonzero(both)
    if len(rows):
        effects[np.ix_(rows, shared_columns)] = (
            0.5 * k562.effects[np.ix_(rows, shared_columns)]
            + 0.5 * pm.effects[np.ix_(rows, shared_columns)]
        )
    return E062.PreparedSource(
        name="k562_gwps_pm50",
        effects=effects,
        target_available=np.asarray(k562.target_available, dtype=bool).copy(),
        mean_on_axis=np.asarray(k562.mean_on_axis, dtype=np.float64).copy(),
        raw_median_norm=float(k562.raw_median_norm),
        mapped_gene_count=int(k562.mapped_gene_count),
        duplicate_symbols=tuple(k562.duplicate_symbols),
    ), both


def build_a15_shaped_baseline(
    surface: str, limit_gib: float
) -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, object]
]:
    held_genes, targets, bctrl, ctrl_cpm = E100.load_controls_and_targets(
        SURFACES[surface]["harness"]
    )
    official = OFFICIAL
    official_pos = {gene: index for index, gene in enumerate(official)}
    panel_positions = np.asarray(
        [official_pos[g] for g in PANEL if g in official_pos], dtype=np.int64
    )
    banks = {
        name: E062.Bank.load(name, REAL_DONOR_PATHS[name])
        for name in SURFACES[surface]["sources"]
    }
    prepared = {
        name: E062.prepare_main_source(
            bank,
            targets,
            official,
            panel_positions,
            chunk_rows=256,
            max_rss_gib=limit_gib,
        )
        for name, bank in banks.items()
    }
    pm_bank = E064.load_perturbmatch(PM_PATH)
    pm_prepared = E062.prepare_main_source(
        pm_bank,
        targets,
        official,
        panel_positions,
        chunk_rows=256,
        max_rss_gib=limit_gib,
    )
    shared = E064.shared_k562_pm_columns(
        banks["k562_gwps"], pm_bank, official, PANEL_SET
    )
    pm_slot, pm_supported = _partial_pm_slot(prepared["k562_gwps"], pm_prepared, shared)
    source_prepared = [
        pm_slot if name == "k562_gwps" else prepared[name]
        for name in SURFACES[surface]["sources"]
    ]
    pre, votes = E062.build_baseline(source_prepared)
    pre *= ALPHA
    source_banks = [banks[name] for name in SURFACES[surface]["sources"]]
    operator, wcols = _operator(surface, source_banks, official)
    whitened, whitening_audit = E062.apply_operator(pre, operator, wcols)
    cd4_banks = [
        E062.Bank.load(name, CD4_ROOT / name, require_cd4_scale=True)
        for name in CD4_NAMES
    ]
    cd4 = E062.prepare_cd4_mean(
        cd4_banks,
        targets,
        official,
        PANEL_SET,
        "gene_norm",
        chunk_rows=256,
        max_rss_gib=limit_gib,
    )
    final, cd4_audit = E062.add_post_whiten_residual(whitened, cd4, CD4_FRACTION, wcols)
    final[:, panel_positions] = 0.0
    official_pos = {gene: index for index, gene in enumerate(official)}
    held_columns = np.asarray(
        [official_pos.get(gene, -1) for gene in held_genes], dtype=np.int64
    )
    final_held = np.zeros((len(targets), len(held_genes)), dtype=np.float64)
    mapped = held_columns >= 0
    final_held[:, mapped] = final[:, held_columns[mapped]]
    audit = {
        "surface": surface,
        "sources": list(SURFACES[surface]["sources"]),
        "alpha": ALPHA,
        "pm_blend": 0.5,
        "pm_supported_targets": int(pm_supported.sum()),
        "cd4_fraction": CD4_FRACTION,
        "cd4_supported_targets": int(np.sum(cd4.target_condition_count > 0)),
        "operator": str(Path(SURFACES[surface]["operator"]).relative_to(ROOT)),
        "operator_sha256": sha256_file(Path(SURFACES[surface]["operator"])),
        "operator_columns": int(len(wcols)),
        "held_projection_genes": int(mapped.sum()),
        "held_unavailable_genes": int((~mapped).sum()),
        "whitening": whitening_audit,
        "cd4": cd4_audit,
        "vote_histogram": {str(v): int(np.sum(votes == v)) for v in np.unique(votes)},
    }
    guard(f"a15-shaped baseline {surface}", limit_gib)
    return held_genes, targets, bctrl, ctrl_cpm, final_held, audit


def predict_scores(
    targets: np.ndarray,
    genes: np.ndarray,
    ctrl_cpm: np.ndarray,
    sources: Sequence[str],
    donors: dict[str, object],
    string: object,
    model: CreditModel,
    kind: str,
) -> tuple[np.ndarray, np.ndarray]:
    builder = E100.FeatureBuilder(genes, ctrl_cpm, sources, string, donors)
    score = np.zeros((len(targets), len(genes)), dtype=np.float32)
    supported = np.zeros_like(score, dtype=bool)
    panel_mask = np.asarray([gene in PANEL_SET for gene in genes], dtype=bool)
    idx = np.flatnonzero(~panel_mask)
    for row, target in enumerate(targets):
        features = builder.features(str(target), idx, 8)
        if features is None:
            continue
        x_agree, x_full = features
        x = x_agree if kind == "agreement" else x_full
        score[row, idx] = model.predict(x).astype(np.float32)
        supported[row, idx] = True
    return score, supported


def freeze(surface: str, limit_gib: float) -> None:
    output = WORK / surface
    if output.exists():
        raise FileExistsError(f"{output} exists; refusing to overwrite")
    train_contexts = tuple(name for name in SURFACES if name != surface)
    donors = {name: E100.DonorLite(name) for name in REAL_DONOR_PATHS}
    string = E100.StringTable()
    blocks = build_training_blocks(surface, train_contexts, string, donors, limit_gib)
    models = fit_models(blocks)
    null_full = [fit_models(blocks, null_draw=draw)["full"] for draw in range(N_NULL)]
    crossfit = crossfit_diagnostics(blocks)
    genes, targets, bctrl, ctrl_cpm, baseline, baseline_audit = (
        build_a15_shaped_baseline(surface, limit_gib)
    )
    active = pds_mask(genes)
    score_agree, support_agree = predict_scores(
        targets,
        genes,
        ctrl_cpm,
        SURFACES[surface]["sources"],
        donors,
        string,
        models["agreement"],
        "agreement",
    )
    score_full, support_full = predict_scores(
        targets,
        genes,
        ctrl_cpm,
        SURFACES[surface]["sources"],
        donors,
        string,
        models["full"],
        "full",
    )
    if not np.array_equal(support_agree, support_full):
        raise AssertionError("matched ablation support differs")
    agreement = transform_effect(
        baseline, score_agree, active, support_agree, beta=BETA
    )
    full = transform_effect(baseline, score_full, active, support_full, beta=BETA)
    null_scores = []
    for model in null_full:
        values, support = predict_scores(
            targets,
            genes,
            ctrl_cpm,
            SURFACES[surface]["sources"],
            donors,
            string,
            model,
            "full",
        )
        if not np.array_equal(support, support_full):
            raise AssertionError("null support differs")
        null_scores.append(values)
    null_serialized = serialize_models(null_full)
    real_serialized = {
        **{
            f"agreement_{key}": value
            for key, value in serialize_models([models["agreement"]]).items()
        },
        **{
            f"full_{key}": value
            for key, value in serialize_models([models["full"]]).items()
        },
    }
    atomic_npz(
        output / "frozen.npz",
        genes=genes.astype("U"),
        targets=targets.astype("U"),
        bctrl=bctrl.astype(np.float32),
        ctrl_cpm=ctrl_cpm.astype(np.float32),
        baseline=baseline.astype(np.float32),
        agreement=agreement.astype(np.float32),
        full=full.astype(np.float32),
        score_agreement=score_agree,
        score_full=score_full,
        support=support_full,
        **{f"null_{key}": value for key, value in null_serialized.items()},
        **real_serialized,
    )
    frozen_path = output / "frozen.npz"
    manifest = {
        "experiment": EXPERIMENT,
        "surface": surface,
        "sealed_before_outer_truth": True,
        "outer_truth_loaded_during_freeze": False,
        "target_visible": True,
        "whole_target_folds": N_FOLDS,
        "n_null": N_NULL,
        "beta": BETA,
        "ridge_alpha": RIDGE_ALPHA,
        "train_contexts": list(train_contexts),
        "training_blocks": len(blocks),
        "training_targets": len({block.target for block in blocks}),
        "crossfit": crossfit,
        "baseline": baseline_audit,
        "masks": {
            "directional": "own target only at score time",
            "pds": "all 300 VCC panel coordinates",
        },
        "frozen_sha256": sha256_file(frozen_path),
        "input_hashes": {
            str(Path(SURFACES[name]["refde"]).relative_to(ROOT)): sha256_file(
                Path(SURFACES[name]["refde"])
            )
            for name in train_contexts
        },
        "implementation_sha256": sha256_file(Path(__file__)),
        "readme_sha256": sha256_file(HERE / "README.md"),
    }
    atomic_json(output / "manifest.json", manifest)
    print(output / "manifest.json")


def target_bootstrap(values: np.ndarray, label: str) -> dict[str, float]:
    values = np.asarray(values, dtype=np.float64)
    rng = np.random.default_rng(stable_seed(EXPERIMENT, label))
    draws = np.empty(N_BOOT, dtype=np.float64)
    for i in range(N_BOOT):
        draws[i] = values[rng.integers(0, len(values), len(values))].mean()
    return {
        "estimate": float(values.mean()),
        "lo": float(np.percentile(draws, 2.5)),
        "hi": float(np.percentile(draws, 97.5)),
    }


def weighted_top_precision(
    correct: np.ndarray, confidence: np.ndarray, weight: np.ndarray
) -> float:
    count = max(1, int(np.ceil(0.2 * len(correct))))
    order = np.argsort(-confidence, kind="stable")[:count]
    denom = float(weight[order].sum())
    return float((correct[order] * weight[order]).sum() / max(denom, 1e-12))


def score(surface: str) -> None:
    output = WORK / surface
    manifest_path = output / "manifest.json"
    frozen_path = output / "frozen.npz"
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("sealed_before_outer_truth") is not True:
        raise ValueError("bundle was not sealed before truth")
    if sha256_file(frozen_path) != manifest.get("frozen_sha256"):
        raise ValueError("frozen bundle hash changed")
    result_path = output / "metrics.json"
    if result_path.exists():
        raise FileExistsError(f"{result_path} exists; refusing overwrite")
    z = np.load(frozen_path, allow_pickle=False)
    genes = z["genes"].astype(str)
    targets = z["targets"].astype(str)
    bctrl = z["bctrl"].astype(np.float64)
    baseline_effect = z["baseline"].astype(np.float64)
    candidate_effects = {
        "identity": baseline_effect,
        "agreement": z["agreement"].astype(np.float64),
        "full": z["full"].astype(np.float64),
    }
    panel_mask = ~pds_mask(genes)
    dh = {
        name: E100.dh_from_effect(
            value, bctrl, z["ctrl_cpm"].astype(np.float64), panel_mask
        )
        for name, value in candidate_effects.items()
    }

    # This is the first held perturbation-truth read.
    ref = E100.RefData.load(surface)
    if not np.array_equal(ref.genes, genes):
        raise ValueError("held gene axis changed")
    ref_rows = {target: i for i, target in enumerate(ref.perts)}
    rows = np.asarray([ref_rows[str(target)] for target in targets], dtype=np.int64)
    pred_lfc = {
        name: E100.lfc_from_dh(value, bctrl, ref.gate, ref.ref_lin)
        for name, value in dh.items()
    }
    per = {
        "identity_credit": [],
        "agreement_credit": [],
        "full_credit": [],
        "agreement_precision": [],
        "full_precision": [],
    }
    used = []
    for out_row, ref_row in enumerate(rows):
        sig = np.asarray(ref.sig[ref_row], dtype=np.int64)
        full_idx = ref.gate_idx[sig]
        own = ~directional_mask(ref.genes[full_idx], str(targets[out_row]))
        truth = ref.lfc[ref_row, sig]
        keep = ~own & np.isfinite(truth) & (truth != 0)
        sig, full_idx, truth = sig[keep], full_idx[keep], truth[keep]
        base = pred_lfc["identity"][out_row, sig]
        keep = np.isfinite(base) & (base != 0)
        sig, full_idx, truth, base = (
            sig[keep],
            full_idx[keep],
            truth[keep],
            base[keep],
        )
        if len(sig) < MIN_GATE:
            continue
        denom = float(np.abs(truth).sum()) or 1.0
        weight = np.abs(base) / denom
        correct = (np.sign(base) == np.sign(truth)).astype(np.float64)
        for name in ("identity", "agreement", "full"):
            pred = pred_lfc[name][out_row, sig]
            per[f"{name}_credit"].append(float(np.sign(truth) @ pred) / denom)
        per["agreement_precision"].append(
            weighted_top_precision(
                correct, z["score_agreement"][out_row, full_idx], weight
            )
        )
        per["full_precision"].append(
            weighted_top_precision(correct, z["score_full"][out_row, full_idx], weight)
        )
        used.append(str(targets[out_row]))
    per_np = {name: np.asarray(value, dtype=np.float64) for name, value in per.items()}
    if not used:
        raise ValueError("no held targets passed directional gate")

    real_effect = E100.load_truth_effects(SURFACES[surface]["harness"], targets, genes)
    pds_columns = pds_mask(genes)
    pds = {
        name: E100.pds_per_target(value, real_effect, pds_columns)
        for name, value in dh.items()
    }

    null_models = deserialize_models(z, "null")
    donors = {name: E100.DonorLite(name) for name in REAL_DONOR_PATHS}
    string = E100.StringTable()
    null_credit = []
    for draw, model in enumerate(null_models):
        null_score, support = predict_scores(
            targets,
            genes,
            z["ctrl_cpm"].astype(np.float64),
            SURFACES[surface]["sources"],
            donors,
            string,
            model,
            "full",
        )
        null_effect = transform_effect(
            baseline_effect, null_score, pds_columns, support, beta=BETA
        )
        null_dh = E100.dh_from_effect(
            null_effect,
            bctrl,
            z["ctrl_cpm"].astype(np.float64),
            panel_mask,
        )
        null_lfc = E100.lfc_from_dh(null_dh, bctrl, ref.gate, ref.ref_lin)
        values = []
        for out_row, ref_row in enumerate(rows):
            sig = np.asarray(ref.sig[ref_row], dtype=np.int64)
            full_idx = ref.gate_idx[sig]
            own = ~directional_mask(ref.genes[full_idx], str(targets[out_row]))
            truth = ref.lfc[ref_row, sig]
            keep = ~own & np.isfinite(truth) & (truth != 0)
            sig, truth = sig[keep], truth[keep]
            base = pred_lfc["identity"][out_row, sig]
            keep = np.isfinite(base) & (base != 0)
            sig, truth = sig[keep], truth[keep]
            if len(sig) < MIN_GATE:
                continue
            denom = float(np.abs(truth).sum()) or 1.0
            values.append(float(np.sign(truth) @ null_lfc[out_row, sig]) / denom)
        null_credit.append(float(np.mean(values)))
        guard(f"null score {surface}/{draw}", MAX_RSS_GIB)

    deltas = {
        "full_minus_identity_credit": (
            per_np["full_credit"] - per_np["identity_credit"]
        ),
        "full_minus_agreement_credit": (
            per_np["full_credit"] - per_np["agreement_credit"]
        ),
        "full_minus_agreement_precision": (
            per_np["full_precision"] - per_np["agreement_precision"]
        ),
        "full_minus_identity_pds": pds["full"] - pds["identity"],
        "agreement_minus_identity_pds": pds["agreement"] - pds["identity"],
    }
    observed = float(np.mean(per_np["full_credit"] - per_np["agreement_credit"]))
    null_gain = np.asarray(null_credit) - float(np.mean(per_np["agreement_credit"]))
    metrics = {
        "surface": surface,
        "targets_total": int(len(targets)),
        "targets_directional": int(len(used)),
        "directional_mask": "reference significant genes excluding own target only",
        "pds_mask": "all 300 VCC panel coordinates excluded",
        "means": {
            name: float(value.mean())
            for name, value in {
                **per_np,
                **{f"pds_{key}": value for key, value in pds.items()},
            }.items()
        },
        "bootstrap": {
            name: target_bootstrap(value, f"{surface}:{name}")
            for name, value in deltas.items()
        },
        "null_full_minus_agreement_credit": {
            "values": null_gain.tolist(),
            "max": float(np.max(null_gain)),
            "p95_conservative": float(np.max(null_gain)),
            "empirical_p": float(
                (1 + np.sum(null_gain >= observed)) / (1 + len(null_gain))
            ),
        },
        "invariants": {
            "panel_coordinates_zero": all(
                not np.any(value[:, panel_mask]) for value in candidate_effects.values()
            ),
            "identity_support_preserved": bool(
                np.array_equal(candidate_effects["full"] == 0, baseline_effect == 0)
            ),
            "signs_preserved": bool(
                np.array_equal(
                    np.sign(candidate_effects["full"]),
                    np.sign(baseline_effect),
                )
            ),
        },
    }
    atomic_json(result_path, metrics)
    print(result_path)


def summarize() -> None:
    rows = {}
    for surface in SURFACES:
        path = WORK / surface / "metrics.json"
        if not path.is_file():
            raise FileNotFoundError(path)
        rows[surface] = json.loads(path.read_text())
    fields = (
        "full_minus_identity_credit",
        "full_minus_agreement_credit",
        "full_minus_agreement_precision",
        "full_minus_identity_pds",
    )
    summary = {"experiment": EXPERIMENT, "surfaces": rows, "macro": {}}
    for field in fields:
        estimates = np.asarray(
            [rows[s]["bootstrap"][field]["estimate"] for s in SURFACES],
            dtype=np.float64,
        )
        summary["macro"][field] = {"estimate": float(estimates.mean())}
    gates = {
        "credit_nonnegative_each_surface": all(
            rows[s]["bootstrap"]["full_minus_identity_credit"]["estimate"] >= 0
            for s in SURFACES
        ),
        "credit_macro_at_least_0_003": (
            summary["macro"]["full_minus_identity_credit"]["estimate"] >= 0.003
        ),
        "full_beats_agreement_credit": (
            summary["macro"]["full_minus_agreement_credit"]["estimate"] > 0
        ),
        "pds_nonnegative_each_surface": all(
            rows[s]["bootstrap"]["full_minus_identity_pds"]["estimate"] >= 0
            for s in SURFACES
        ),
        "precision_nonnegative_each_surface": all(
            rows[s]["bootstrap"]["full_minus_agreement_precision"]["estimate"] >= 0
            for s in SURFACES
        ),
        "precision_macro_at_least_0_01": (
            summary["macro"]["full_minus_agreement_precision"]["estimate"] >= 0.01
        ),
        "beats_each_surface_null_max": all(
            rows[s]["bootstrap"]["full_minus_agreement_credit"]["estimate"]
            > rows[s]["null_full_minus_agreement_credit"]["max"]
            for s in SURFACES
        ),
        "all_invariants": all(all(rows[s]["invariants"].values()) for s in SURFACES),
    }
    summary["gates"] = gates
    summary["verdict"] = "GO" if all(gates.values()) else "NO-GO"
    atomic_json(WORK / "summary.json", summary)
    print(json.dumps(summary["macro"], indent=2, sort_keys=True))
    print(summary["verdict"])


def audit_inputs() -> None:
    required = [
        HERE / "README.md",
        ROOT / "data/challenge/gene_names.csv",
        ROOT / "data/challenge/pert_counts.csv",
        STRING_PATH,
        PM_PATH / "delta.npy",
        PM_PATH / "genes.txt",
        PM_PATH / "perts.txt",
        PM_PATH / "meta.json",
    ]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    for path in REAL_DONOR_PATHS.values():
        for name in ("delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt", "meta.json"):
            if not (path / name).is_file():
                raise FileNotFoundError(path / name)
    for name in CD4_NAMES:
        for filename in (
            "delta.npy",
            "b_ctrl.npy",
            "genes.txt",
            "perts.txt",
            "meta.json",
        ):
            if not (CD4_ROOT / name / filename).is_file():
                raise FileNotFoundError(CD4_ROOT / name / filename)
    for config in SURFACES.values():
        for key in ("harness", "refde", "operator"):
            if not Path(config[key]).is_file():
                raise FileNotFoundError(config[key])
    report = {
        "experiment": EXPERIMENT,
        "kolf_paths_referenced": 0,
        "panel_targets": len(PANEL),
        "inputs": {
            str(path.relative_to(ROOT)): sha256_file(path)
            for path in required
            if path.suffix != ".h5"
        },
        "surfaces": {
            name: {
                "sources": list(config["sources"]),
                "harness": str(Path(config["harness"]).relative_to(ROOT)),
                "refde": str(Path(config["refde"]).relative_to(ROOT)),
                "operator": str(Path(config["operator"]).relative_to(ROOT)),
            }
            for name, config in SURFACES.items()
        },
    }
    atomic_json(WORK / "input-audit.json", report)
    print(WORK / "input-audit.json")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    sub = result.add_subparsers(dest="command", required=True)
    sub.add_parser("audit-inputs")
    freeze_parser = sub.add_parser("freeze")
    freeze_parser.add_argument("--surface", required=True, choices=tuple(SURFACES))
    freeze_parser.add_argument("--max-rss-gib", type=float, default=24.0)
    score_parser = sub.add_parser("score")
    score_parser.add_argument("--surface", required=True, choices=tuple(SURFACES))
    sub.add_parser("summarize")
    return result


def main() -> None:
    args = parser().parse_args()
    if args.command == "audit-inputs":
        audit_inputs()
    elif args.command == "freeze":
        freeze(args.surface, args.max_rss_gib)
    elif args.command == "score":
        score(args.surface)
    elif args.command == "summarize":
        summarize()
    else:
        raise AssertionError(args.command)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Experiment 100: leakage-safe, target-blocked sign-confidence models.

Commands:
    uv run python experiments/100-sign-confidence/run.py prepare
    uv run python experiments/100-sign-confidence/run.py freeze --surface h1
    uv run python experiments/100-sign-confidence/run.py score --surface h1

All generated files stay below this experiment directory. The implementation
does not invoke the cell-eval2 CLI or build a challenge artifact.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import os
import resource
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "8")

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
STATS = WORK / "donor-stats"
STRING_CACHE = WORK / "string-pca.npz"
TS = 50_000.0
ALPHA = 1.08
EFFECT_SCALE = 6.817
WHITEN_GAMMA = 0.65
MIN_CPM = 5.0
N_FOLDS = 4
N_NULL = 20
N_BOOT = 2_000
BETA = 0.15
PANEL = tuple(
    pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str)
)
PANEL_SET = set(PANEL)

DONOR_PATHS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
    "h1_2025": ROOT / "data/donors/h1_2025",
    "kolf2": ROOT / "data/external/kolf2",
}
DONOR_POOL = tuple(DONOR_PATHS)

SURFACES = {
    "h1": {
        "truth_donor": "h1_2025",
        "harness": ROOT / "data/harness/h1-test.h5ad",
        "refde": ROOT / "experiments/098-nmae/work/refde_h1.npz",
        "sources": ("k562_gwps", "orion_hek293t"),
        "whitener": ROOT / "experiments/095-clean-harness/out-h1test/whitener_union.npz",
    },
    "hct116": {
        "truth_donor": "orion_hct116",
        "harness": ROOT / "data/harness/orion-hct116-thinned.h5ad",
        "refde": ROOT / "experiments/098-nmae/work/refde_orion.npz",
        "sources": ("k562_gwps", "orion_hek293t"),
        "whitener": ROOT / "experiments/095-clean-harness/out/whitener_union.npz",
    },
    "hek293t": {
        "truth_donor": "orion_hek293t",
        "harness": ROOT / "data/harness/orion-hek293t-thinned.h5ad",
        "refde": ROOT / "experiments/099-donors/work/refde_hek293t.npz",
        "sources": ("k562_gwps", "orion_hct116"),
        "whitener": ROOT / "experiments/099-donors/work/whit/whitener_hek293t.npz",
    },
}


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def guard(limit_gib: float, stage: str) -> None:
    peak = peak_rss_gib()
    if limit_gib > 0 and peak >= limit_gib:
        raise MemoryError(f"{stage}: peak RSS {peak:.2f} GiB >= {limit_gib:.2f} GiB")


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256("\n".join(map(str, values)).encode()).hexdigest()


def sha256_array(values: np.ndarray) -> str:
    a = np.ascontiguousarray(values)
    return hashlib.sha256(a.view(np.uint8)).hexdigest()


def stable_seed(*parts: object) -> int:
    return int.from_bytes(
        hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:4], "little"
    )


def target_fold(target: str) -> int:
    return stable_seed("exp100-fold", target) % N_FOLDS


def read_lines(path: Path) -> tuple[str, ...]:
    return tuple(x.strip() for x in path.read_text().splitlines() if x.strip())


def last_lookup(values: Sequence[str]) -> dict[str, int]:
    out: dict[str, int] = {}
    for i, value in enumerate(values):
        out[str(value)] = i
    return out


def json_dump(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


@dataclass
class RefData:
    name: str
    genes: np.ndarray
    gate: np.ndarray
    gate_idx: np.ndarray
    ref_lin: np.ndarray
    perts: np.ndarray
    lfc: np.ndarray
    sig: list[np.ndarray]
    ctrl_cpm: np.ndarray

    @classmethod
    def load(cls, name: str) -> "RefData":
        z = np.load(SURFACES[name]["refde"], allow_pickle=False)
        genes = z["genes"].astype(str)
        gate = z["gate"].astype(bool)
        gate_idx = np.flatnonzero(gate)
        ref_lin = z["ref_lin"].astype(np.float64)
        if len(gate_idx) != len(ref_lin):
            raise ValueError(f"{name}: gate/ref_lin mismatch")
        perts = z["perts"].astype(str)
        lfc = z["lfc"].astype(np.float64)
        ptr, idx = z["sig_ptr"], z["sig_idx"]
        sig = [idx[ptr[i] : ptr[i + 1]].astype(np.int64) for i in range(len(perts))]
        ctrl = np.zeros(len(genes), dtype=np.float64)
        ctrl[gate_idx] = ref_lin
        return cls(name, genes, gate, gate_idx, ref_lin, perts, lfc, sig, ctrl)


class DonorLite:
    """Memory-mapped donor with cached exact centering and axis-specific norm."""

    def __init__(self, name: str):
        self.name = name
        self.path = DONOR_PATHS[name]
        self.perts = read_lines(self.path / "perts.txt")
        self.genes = read_lines(self.path / "genes.txt")
        self.pert_index = last_lookup(self.perts)
        self.gene_index = last_lookup(self.genes)
        self.delta = np.load(self.path / "delta.npy", mmap_mode="r")
        if self.delta.shape != (len(self.perts), len(self.genes)):
            raise ValueError(f"{name}: compact axis mismatch {self.delta.shape}")
        STATS.mkdir(parents=True, exist_ok=True)
        self.mean = self._load_mean()
        self._projections: dict[str, tuple[np.ndarray, np.ndarray, float]] = {}

    def _load_mean(self) -> np.ndarray:
        path = STATS / f"{self.name}-mean.npy"
        if path.exists():
            out = np.load(path)
            if out.shape == (self.delta.shape[1],):
                return out.astype(np.float32, copy=False)
        total = np.zeros(self.delta.shape[1], dtype=np.float64)
        for start in range(0, self.delta.shape[0], 256):
            total += np.asarray(self.delta[start : start + 256], dtype=np.float64).sum(0)
        out = (total / self.delta.shape[0]).astype(np.float32)
        np.save(path, out)
        return out

    def projection(
        self, axis: Sequence[str], allowed_mask: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, float]:
        axis_hash = sha256_lines(axis)[:16]
        mask_hash = sha256_array(np.asarray(allowed_mask, dtype=np.uint8))[:16]
        key = f"{axis_hash}-{mask_hash}"
        if key in self._projections:
            return self._projections[key]
        tcols, dcols = [], []
        for target_col, gene in enumerate(axis):
            donor_col = self.gene_index.get(str(gene))
            if donor_col is not None and allowed_mask[target_col]:
                tcols.append(target_col)
                dcols.append(donor_col)
        tc = np.asarray(tcols, dtype=np.int64)
        dc = np.asarray(dcols, dtype=np.int64)
        scale_path = STATS / f"{self.name}-scale-{key}.json"
        if scale_path.exists():
            scale = float(json.loads(scale_path.read_text())["median_centered_norm"])
        else:
            norms = np.empty(self.delta.shape[0], dtype=np.float64)
            mu = self.mean[dc]
            for start in range(0, self.delta.shape[0], 256):
                stop = min(start + 256, self.delta.shape[0])
                block = np.asarray(
                    self.delta[start:stop, dc], dtype=np.float32
                ).copy()
                block -= mu
                norms[start:stop] = np.sqrt(
                    np.einsum("ij,ij->i", block, block, dtype=np.float64)
                )
            scale = float(np.median(norms)) or 1.0
            json_dump(
                scale_path,
                {
                    "donor": self.name,
                    "axis_sha256": sha256_lines(axis),
                    "allowed_mask_sha256": sha256_array(
                        np.asarray(allowed_mask, dtype=np.uint8)
                    ),
                    "shared_genes": int(len(tc)),
                    "median_centered_norm": scale,
                },
            )
        self._projections[key] = (tc, dc, scale)
        return tc, dc, scale

    def vector(
        self, target: str, axis: Sequence[str], allowed_mask: np.ndarray
    ) -> np.ndarray | None:
        row = self.pert_index.get(target)
        if row is None:
            return None
        tc, dc, scale = self.projection(axis, allowed_mask)
        out = np.zeros(len(axis), dtype=np.float32)
        out[tc] = (np.asarray(self.delta[row, dc], dtype=np.float32) - self.mean[dc]) / scale
        return out


class StringTable:
    def __init__(self):
        self.genes: np.ndarray
        self.pcs: np.ndarray
        self.unit: np.ndarray
        self.available: np.ndarray
        if STRING_CACHE.exists():
            z = np.load(STRING_CACHE, allow_pickle=False)
            self.genes = z["genes"].astype(str)
            self.pcs = z["pcs"].astype(np.float32)
            self.unit = z["unit"].astype(np.float32)
            self.available = z["available"].astype(bool)
        else:
            self._build()
        self.index = {g: i for i, g in enumerate(self.genes)}

    def _build(self) -> None:
        path = ROOT / "data/external/slim/gene_string_embeddings.v0.3.h5"
        with h5py.File(path, "r") as handle:
            genes = np.asarray(sorted(handle.keys()))
            values = np.vstack(
                [np.asarray(handle[str(g)], dtype=np.float64) for g in genes]
            )
        available = np.linalg.norm(values, axis=1) > 0
        fit = values[available]
        mean = fit.mean(0)
        scale = fit.std(0)
        scale[scale == 0] = 1.0
        z = (fit - mean) / scale
        cov = (z.T @ z) / z.shape[0]
        vals, vecs = np.linalg.eigh(cov)
        components = vecs[:, np.argsort(vals)[::-1][:16]]
        pcs = np.zeros((len(genes), 16), dtype=np.float32)
        pcs[available] = (
            ((values[available] - mean) / scale) @ components
        ).astype(np.float32)
        unit = np.zeros_like(values, dtype=np.float32)
        norms = np.linalg.norm(values[available], axis=1, keepdims=True)
        unit[available] = (values[available] / np.maximum(norms, 1e-12)).astype(
            np.float32
        )
        STRING_CACHE.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            STRING_CACHE,
            genes=genes,
            pcs=pcs,
            unit=unit,
            available=available,
            mean=mean,
            scale=scale,
            components=components,
        )
        self.genes, self.pcs, self.unit, self.available = (
            genes,
            pcs,
            unit,
            available,
        )

    def on_axis(
        self, axis: Sequence[str]
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        rows = np.array([self.index.get(str(g), -1) for g in axis], dtype=np.int64)
        ok = rows >= 0
        pcs = np.zeros((len(axis), 16), dtype=np.float32)
        unit = np.zeros((len(axis), 64), dtype=np.float32)
        avail = np.zeros(len(axis), dtype=bool)
        pcs[ok] = self.pcs[rows[ok]]
        unit[ok] = self.unit[rows[ok]]
        avail[ok] = self.available[rows[ok]]
        return pcs, unit, avail

    def one(self, gene: str) -> tuple[np.ndarray, np.ndarray, bool]:
        row = self.index.get(gene)
        if row is None:
            return (
                np.zeros(16, dtype=np.float32),
                np.zeros(64, dtype=np.float32),
                False,
            )
        return self.pcs[row], self.unit[row], bool(self.available[row])


@dataclass
class TargetState:
    base: np.ndarray
    mean_abs: np.ndarray
    agreement_ratio: np.ndarray
    unanimous: np.ndarray
    min_ratio: np.ndarray
    std_ratio: np.ndarray
    n_carriers: int
    pair_cos: float


class FeatureBuilder:
    def __init__(
        self,
        axis: Sequence[str],
        ctrl_cpm: np.ndarray,
        source_names: Sequence[str],
        string: StringTable,
        donors: dict[str, DonorLite],
    ):
        self.axis = np.asarray(axis).astype(str)
        self.ctrl_cpm = np.asarray(ctrl_cpm, dtype=np.float64)
        self.source_names = tuple(source_names)
        self.string = string
        self.donors = donors
        self.allowed = np.asarray([g not in PANEL_SET for g in self.axis], dtype=bool)
        self.axis_pcs, self.axis_unit, self.axis_string_ok = string.on_axis(self.axis)
        self._state: dict[str, TargetState | None] = {}

    def target_state(self, target: str) -> TargetState | None:
        if target in self._state:
            return self._state[target]
        vectors = []
        for name in self.source_names:
            value = self.donors[name].vector(target, self.axis, self.allowed)
            if value is not None:
                value[~self.allowed] = 0.0
                vectors.append(value)
        if not vectors:
            self._state[target] = None
            return None
        stack = np.vstack(vectors).astype(np.float64)
        base = stack.mean(0)
        abs_stack = np.abs(stack)
        mean_abs = abs_stack.mean(0)
        eps = 1e-12
        signs = np.sign(stack)
        unanimous = (
            (np.abs(signs.sum(0)) == len(vectors)) & np.all(stack != 0, axis=0)
        ).astype(np.float64)
        pair = []
        active = self.allowed
        for i in range(len(vectors)):
            for j in range(i + 1, len(vectors)):
                a, b = stack[i, active], stack[j, active]
                pair.append(float(a @ b) / max(np.linalg.norm(a) * np.linalg.norm(b), eps))
        state = TargetState(
            base=base,
            mean_abs=mean_abs,
            agreement_ratio=np.abs(base) / (mean_abs + eps),
            unanimous=unanimous,
            min_ratio=abs_stack.min(0) / (mean_abs + eps),
            std_ratio=stack.std(0) / (mean_abs + eps),
            n_carriers=len(vectors),
            pair_cos=float(np.mean(pair)) if pair else 0.0,
        )
        self._state[target] = state
        return state

    def features(
        self, target: str, indices: np.ndarray, rank: int
    ) -> tuple[np.ndarray, np.ndarray] | None:
        state = self.target_state(target)
        if state is None:
            return None
        ix = np.asarray(indices, dtype=np.int64)
        base = state.base[ix]
        agree = np.column_stack(
            [
                np.log1p(np.abs(base)),
                state.agreement_ratio[ix],
                state.unanimous[ix],
                state.min_ratio[ix],
                state.std_ratio[ix],
                np.full(len(ix), state.n_carriers / len(DONOR_POOL)),
                np.full(len(ix), state.pair_cos),
            ]
        ).astype(np.float32)
        tpc, tunit, tok = self.string.one(target)
        rpc = self.axis_pcs[ix, :rank]
        runit = self.axis_unit[ix]
        rok = self.axis_string_ok[ix]
        cosine = runit @ tunit if tok else np.zeros(len(ix), dtype=np.float32)
        full = np.column_stack(
            [
                agree,
                np.sign(base),
                np.log1p(np.maximum(self.ctrl_cpm[ix], 0.0)),
                cosine,
                np.repeat(tpc[None, :rank], len(ix), axis=0),
                rpc,
                rpc * tpc[None, :rank],
                np.full(len(ix), float(not tok)),
                (~rok).astype(np.float32),
            ]
        ).astype(np.float32)
        return agree, full


@dataclass
class TrainBlock:
    context: str
    target: str
    context_row: int
    gate_indices: np.ndarray
    x_agree: np.ndarray
    x8: np.ndarray
    x16: np.ndarray
    y: np.ndarray
    weight: np.ndarray


@dataclass
class TrainContext:
    name: str
    blocks: list[TrainBlock]
    dense_y: np.ndarray
    eligible_rows: np.ndarray


@dataclass
class LinearModel:
    mean: np.ndarray
    scale: np.ndarray
    coef: np.ndarray
    intercept: float

    def predict(self, x: np.ndarray) -> np.ndarray:
        z = (np.asarray(x, dtype=np.float64) - self.mean) / self.scale
        logits = z @ self.coef + self.intercept
        logits = np.clip(logits, -30.0, 30.0)
        return 1.0 / (1.0 + np.exp(-logits))


def fit_logistic(x: np.ndarray, y: np.ndarray, weight: np.ndarray) -> LinearModel:
    if len(np.unique(y)) != 2:
        raise ValueError("training labels contain fewer than two classes")
    w = np.asarray(weight, dtype=np.float64)
    w *= len(w) / max(w.sum(), 1e-12)
    scaler = StandardScaler()
    scaler.fit(x, sample_weight=w)
    scale = np.asarray(scaler.scale_, dtype=np.float64)
    scale[scale == 0] = 1.0
    z = (np.asarray(x, dtype=np.float64) - scaler.mean_) / scale
    model = LogisticRegression(
        C=1.0,
        penalty="l2",
        solver="lbfgs",
        max_iter=200,
        tol=1e-6,
        random_state=0,
    )
    model.fit(z, y, sample_weight=w)
    return LinearModel(
        mean=np.asarray(scaler.mean_, dtype=np.float64),
        scale=scale,
        coef=np.asarray(model.coef_[0], dtype=np.float64),
        intercept=float(model.intercept_[0]),
    )


def source_names_for_training(outer: str, truth_context: str) -> tuple[str, ...]:
    blocked = {
        SURFACES[outer]["truth_donor"],
        SURFACES[truth_context]["truth_donor"],
    }
    out = tuple(name for name in DONOR_POOL if name not in blocked)
    if blocked & set(out):
        raise AssertionError("outer/truth donor entered a training source set")
    return out


def build_training_context(
    outer: str,
    context: str,
    string: StringTable,
    donors: dict[str, DonorLite],
    max_rss_gib: float,
) -> TrainContext:
    ref = RefData.load(context)
    sources = source_names_for_training(outer, context)
    builder = FeatureBuilder(ref.genes, ref.ctrl_cpm, sources, string, donors)
    eligible_gate = np.asarray(
        [ref.genes[i] not in PANEL_SET for i in ref.gate_idx], dtype=bool
    )
    dense_y = np.zeros((len(ref.perts), len(ref.gate_idx)), dtype=np.uint8)
    base_by_row: list[np.ndarray | None] = []
    eligible_rows = []
    for row, target in enumerate(ref.perts):
        state = builder.target_state(str(target))
        if state is None or state.n_carriers < 2:
            base_by_row.append(None)
            continue
        base_gate = state.base[ref.gate_idx]
        dense_y[row] = (
            (np.sign(base_gate) == np.sign(ref.lfc[row]))
            & (base_gate != 0)
            & np.isfinite(ref.lfc[row])
            & (ref.lfc[row] != 0)
            & eligible_gate
        )
        base_by_row.append(state.base)
        eligible_rows.append(row)

    blocks: list[TrainBlock] = []
    for row in eligible_rows:
        target = str(ref.perts[row])
        sig = ref.sig[row]
        if len(sig) < 10:
            continue
        keep = (
            eligible_gate[sig]
            & np.isfinite(ref.lfc[row, sig])
            & (ref.lfc[row, sig] != 0)
        )
        sig = sig[keep]
        if len(sig) < 10:
            continue
        full_idx = ref.gate_idx[sig]
        feats8 = builder.features(target, full_idx, 8)
        feats16 = builder.features(target, full_idx, 16)
        if feats8 is None or feats16 is None:
            continue
        xa, x8 = feats8
        _, x16 = feats16
        base = base_by_row[row]
        assert base is not None
        nonzero = base[full_idx] != 0
        sig, full_idx = sig[nonzero], full_idx[nonzero]
        xa, x8, x16 = xa[nonzero], x8[nonzero], x16[nonzero]
        if len(sig) < 10:
            continue
        truth = ref.lfc[row, sig]
        y = (np.sign(base[full_idx]) == np.sign(truth)).astype(np.uint8)
        denom = float(np.abs(truth).sum()) or 1.0
        weight = np.abs(base[full_idx]) / denom
        if weight.sum() <= 0:
            continue
        blocks.append(
            TrainBlock(
                context=context,
                target=target,
                context_row=row,
                gate_indices=sig,
                x_agree=xa,
                x8=x8,
                x16=x16,
                y=y,
                weight=weight.astype(np.float64),
            )
        )
    guard(max_rss_gib, f"training context {context}")
    return TrainContext(
        name=context,
        blocks=blocks,
        dense_y=dense_y,
        eligible_rows=np.asarray(eligible_rows, dtype=np.int64),
    )


def collect_training(
    contexts: Sequence[TrainContext],
    held_targets: set[str],
    kind: str,
    null_seed: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xs, ys, ws = [], [], []
    for ctx in contexts:
        eligible = [b for b in ctx.blocks if b.target not in held_targets]
        null_map: dict[int, int] = {}
        if null_seed is not None:
            rows = np.asarray([b.context_row for b in eligible], dtype=np.int64)
            rng = np.random.default_rng(stable_seed("null", null_seed, ctx.name))
            perm = rng.permutation(rows)
            null_map = {int(a): int(b) for a, b in zip(rows, perm)}
        for block in eligible:
            if kind == "agree":
                x = block.x_agree
            elif kind == "rank8":
                x = block.x8
            elif kind == "rank16":
                x = block.x16
            else:
                raise ValueError(kind)
            y = block.y
            if null_seed is not None:
                source_row = null_map[block.context_row]
                y = ctx.dense_y[source_row, block.gate_indices]
            xs.append(x)
            ys.append(y)
            ws.append(block.weight)
    if not xs:
        raise ValueError("no training examples after target blocking")
    return np.vstack(xs), np.concatenate(ys), np.concatenate(ws)


def fit_fold_models(
    contexts: Sequence[TrainContext],
    held_targets: set[str],
) -> dict[str, LinearModel]:
    out = {}
    for kind in ("agree", "rank8", "rank16"):
        x, y, w = collect_training(contexts, held_targets, kind)
        out[kind] = fit_logistic(x, y, w)
    return out


def load_controls_and_targets(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    handle = ad.read_h5ad(path, backed="r")
    genes = np.asarray(handle.var_names).astype(str)
    labels = np.asarray(handle.obs["target_gene"].astype(str))
    ctrl_rows = np.flatnonzero(labels == "non-targeting")
    xctrl = handle.X[ctrl_rows, :]
    if sp.issparse(xctrl):
        sums = np.asarray(xctrl.sum(0)).ravel().astype(np.float64)
    else:
        sums = np.asarray(xctrl, dtype=np.float64).sum(0)
    handle.file.close()
    bctrl = np.log1p(TS * sums / sums.sum())
    ctrl_cpm = 1e6 * sums / sums.sum()
    targets = np.asarray(sorted(set(labels) - {"non-targeting"}))
    return genes, targets, bctrl, ctrl_cpm


def derive_wcols(
    genes: np.ndarray,
    targets: np.ndarray,
    sources: Sequence[str],
    donors: dict[str, DonorLite],
    whitener_path: Path,
) -> np.ndarray:
    union = set()
    for source in sources:
        union |= set(donors[source].genes)
    excluded = set(PANEL) if set(targets) <= PANEL_SET else set(map(str, targets))
    union -= excluded
    lookup = {g: i for i, g in enumerate(genes)}
    wcols = np.asarray(sorted(lookup[g] for g in union if g in lookup), dtype=np.int64)
    meta_path = whitener_path.with_suffix(".json")
    meta = json.loads(meta_path.read_text())
    expected = meta.get("wcols_sha256")
    if expected is not None and sha256_array(wcols) != expected:
        raise ValueError(
            f"whitener axis mismatch: {len(wcols)} {sha256_array(wcols)[:16]} "
            f"!= {meta.get('n_cols')} {str(expected)[:16]}"
        )
    return wcols


def build_identity_effect(
    genes: np.ndarray,
    targets: np.ndarray,
    sources: Sequence[str],
    donors: dict[str, DonorLite],
    whitener_path: Path,
    max_rss_gib: float,
) -> tuple[np.ndarray, dict[str, object]]:
    allowed = np.asarray([g not in PANEL_SET for g in genes], dtype=bool)
    rows = []
    carriers = []
    for target in targets:
        values = [
            donors[name].vector(str(target), genes, allowed) for name in sources
        ]
        values = [v for v in values if v is not None]
        if not values:
            rows.append(np.zeros(len(genes), dtype=np.float64))
            carriers.append(0)
        else:
            rows.append(np.mean(values, axis=0))
            carriers.append(len(values))
    effect = np.vstack(rows).astype(np.float64)
    forbidden = ~allowed
    effect[:, forbidden] = 0.0
    wcols = derive_wcols(genes, targets, sources, donors, whitener_path)
    z = np.load(whitener_path)
    evals = z["evals"].astype(np.float64)
    evecs = z["evecs"].astype(np.float64)
    if evecs.shape != (len(wcols), len(wcols)):
        raise ValueError("whitener/effect axis mismatch")
    before = float(np.median(np.linalg.norm(effect[:, wcols], axis=1))) or 1.0
    weights = evals ** (-WHITEN_GAMMA / 2.0)
    transformed = ((effect[:, wcols] @ evecs) * weights) @ evecs.T
    after = float(np.median(np.linalg.norm(transformed, axis=1))) or 1.0
    effect[:, wcols] = transformed * (before / after)
    effect[:, forbidden] = 0.0
    scale = ALPHA * EFFECT_SCALE / (
        float(np.median(np.linalg.norm(effect[:, allowed], axis=1))) or 1.0
    )
    effect *= scale
    guard(max_rss_gib, "whitening")
    return effect, {
        "sources": list(sources),
        "carriers_min": int(np.min(carriers)),
        "carriers_median": float(np.median(carriers)),
        "carriers_max": int(np.max(carriers)),
        "wcols": int(len(wcols)),
        "wcols_sha256": sha256_array(wcols),
        "whitener": str(whitener_path),
        "before_whitening_median_norm": before,
        "after_whitening_median_norm": after,
        "final_scale": scale,
    }


def probability_multiplier(
    probability: np.ndarray,
    effect: np.ndarray,
    active: np.ndarray,
) -> np.ndarray:
    p = np.clip(np.asarray(probability, dtype=np.float64), 1e-5, 1 - 1e-5)
    logit = np.log(p / (1 - p))
    use = active & (effect != 0)
    out = np.ones_like(effect, dtype=np.float64)
    if use.sum() < 2:
        return out
    z = logit[use]
    z = (z - z.mean()) / (z.std() or 1.0)
    z = np.clip(z, -2.0, 2.0)
    out[use] = np.exp(BETA * z)
    return out


def transform_effect(
    effect: np.ndarray,
    probability: np.ndarray,
    active: np.ndarray,
) -> np.ndarray:
    out = np.asarray(effect, dtype=np.float64).copy()
    for row in range(len(out)):
        before = float(np.linalg.norm(out[row, active]))
        out[row] *= probability_multiplier(probability[row], out[row], active)
        after = float(np.linalg.norm(out[row, active])) or 1.0
        out[row] *= before / after
    out[:, ~active] = 0.0
    return out


def dh_from_effect(
    effect: np.ndarray, bctrl: np.ndarray, ctrl_cpm: np.ndarray, panel_mask: np.ndarray
) -> np.ndarray:
    gated = np.asarray(effect, dtype=np.float64).copy()
    gated[:, ctrl_cpm <= MIN_CPM] = 0.0
    gated[:, bctrl <= 0] = 0.0
    profile = np.expm1(bctrl[None, :] + gated)
    np.clip(profile, 0.0, None, out=profile)
    profile /= profile.sum(1, keepdims=True)
    out = np.log1p(TS * profile) - bctrl[None, :]
    out[:, panel_mask] = 0.0
    return out


def serialize_models(
    path: Path, models: list[list[LinearModel]], feature_dim: int
) -> None:
    means = np.zeros((len(models), N_FOLDS, feature_dim), dtype=np.float32)
    scales = np.zeros_like(means)
    coefs = np.zeros_like(means)
    intercepts = np.zeros((len(models), N_FOLDS), dtype=np.float64)
    for draw, fold_models in enumerate(models):
        for fold, model in enumerate(fold_models):
            means[draw, fold] = model.mean
            scales[draw, fold] = model.scale
            coefs[draw, fold] = model.coef
            intercepts[draw, fold] = model.intercept
    np.savez_compressed(
        path, means=means, scales=scales, coefs=coefs, intercepts=intercepts
    )


def freeze_surface(surface: str, max_rss_gib: float, n_null: int) -> None:
    outdir = WORK / surface
    outdir.mkdir(parents=True, exist_ok=True)
    frozen_path = outdir / "frozen_predictions.npz"
    if frozen_path.exists():
        raise FileExistsError(f"{frozen_path} already exists; refusing to overwrite")
    string = StringTable()
    donors = {name: DonorLite(name) for name in DONOR_POOL}
    train_names = tuple(name for name in SURFACES if name != surface)
    train_contexts = [
        build_training_context(surface, name, string, donors, max_rss_gib)
        for name in train_names
    ]

    genes, targets, bctrl, ctrl_cpm = load_controls_and_targets(
        SURFACES[surface]["harness"]
    )
    panel_mask = np.asarray([g in PANEL_SET for g in genes], dtype=bool)
    active = (~panel_mask) & (ctrl_cpm > MIN_CPM)
    final_builder = FeatureBuilder(
        genes,
        ctrl_cpm,
        SURFACES[surface]["sources"],
        string,
        donors,
    )

    fold_models: list[dict[str, LinearModel]] = []
    fold_audit = []
    for fold in range(N_FOLDS):
        held = {str(t) for t in targets if target_fold(str(t)) == fold}
        models = fit_fold_models(train_contexts, held)
        fold_models.append(models)
        fold_audit.append(
            {
                "fold": fold,
                "held_targets": sorted(held),
                "held_target_sha256": sha256_lines(sorted(held)),
                "coef_norms": {
                    kind: float(np.linalg.norm(model.coef))
                    for kind, model in models.items()
                },
            }
        )
        guard(max_rss_gib, f"fit fold {fold}")

    p_agree = np.full((len(targets), len(genes)), 0.5, dtype=np.float32)
    p8 = np.full_like(p_agree, 0.5)
    p16 = np.full_like(p_agree, 0.5)
    for row, target in enumerate(targets):
        fold = target_fold(str(target))
        idx = np.flatnonzero(active)
        f8 = final_builder.features(str(target), idx, 8)
        f16 = final_builder.features(str(target), idx, 16)
        if f8 is None or f16 is None:
            continue
        xa, x8 = f8
        _, x16 = f16
        p_agree[row, idx] = fold_models[fold]["agree"].predict(xa)
        p8[row, idx] = fold_models[fold]["rank8"].predict(x8)
        p16[row, idx] = fold_models[fold]["rank16"].predict(x16)
    p_agree[:, panel_mask] = 0.5
    p8[:, panel_mask] = 0.5
    p16[:, panel_mask] = 0.5

    identity_effect, identity_audit = build_identity_effect(
        genes,
        targets,
        SURFACES[surface]["sources"],
        donors,
        SURFACES[surface]["whitener"],
        max_rss_gib,
    )
    agreement_effect = transform_effect(identity_effect, p_agree, active)
    rank8_effect = transform_effect(identity_effect, p8, active)
    rank16_effect = transform_effect(identity_effect, p16, active)
    dh_identity = dh_from_effect(identity_effect, bctrl, ctrl_cpm, panel_mask)
    dh_agreement = dh_from_effect(agreement_effect, bctrl, ctrl_cpm, panel_mask)
    dh_rank8 = dh_from_effect(rank8_effect, bctrl, ctrl_cpm, panel_mask)
    dh_rank16 = dh_from_effect(rank16_effect, bctrl, ctrl_cpm, panel_mask)

    # Freeze the actual predictions before any held-surface truth values are read.
    np.savez_compressed(
        frozen_path,
        genes=genes,
        targets=targets,
        bctrl=bctrl.astype(np.float64),
        ctrl_cpm=ctrl_cpm.astype(np.float64),
        panel_mask=panel_mask,
        p_agree=p_agree,
        p_rank8=p8,
        p_rank16=p16,
        dh_identity=dh_identity.astype(np.float32),
        dh_agreement=dh_agreement.astype(np.float32),
        dh_rank8=dh_rank8.astype(np.float32),
        dh_rank16=dh_rank16.astype(np.float32),
    )

    # Freeze complete null model parameters. Held truth is not needed to fit them.
    null_draw_models: list[list[LinearModel]] = []
    for draw in range(n_null):
        per_fold = []
        for fold in range(N_FOLDS):
            held = {str(t) for t in targets if target_fold(str(t)) == fold}
            x, y, w = collect_training(
                train_contexts, held, "rank8", null_seed=draw
            )
            per_fold.append(fit_logistic(x, y, w))
        null_draw_models.append(per_fold)
        guard(max_rss_gib, f"null draw {draw}")
    if null_draw_models:
        feature_dim = len(null_draw_models[0][0].coef)
        serialize_models(outdir / "null_models.npz", null_draw_models, feature_dim)
    else:
        raise ValueError("at least one permutation null is required")

    audit = {
        "surface": surface,
        "outer_truth_donor": SURFACES[surface]["truth_donor"],
        "final_sources": list(SURFACES[surface]["sources"]),
        "training_contexts": list(train_names),
        "training_source_sets": {
            name: list(source_names_for_training(surface, name)) for name in train_names
        },
        "training_blocks": {
            ctx.name: len(ctx.blocks) for ctx in train_contexts
        },
        "training_examples": {
            ctx.name: int(sum(len(b.y) for b in ctx.blocks)) for ctx in train_contexts
        },
        "target_folds": {str(t): target_fold(str(t)) for t in targets},
        "folds": fold_audit,
        "identity": identity_audit,
        "panel_coordinates_removed": int(panel_mask.sum()),
        "active_response_coordinates": int(active.sum()),
        "n_null": n_null,
        "peak_rss_gib": peak_rss_gib(),
        "frozen_sha256": hashlib.sha256(frozen_path.read_bytes()).hexdigest(),
    }
    json_dump(outdir / "freeze-audit.json", audit)
    print(json.dumps(audit, indent=2))


def load_truth_effects(
    path: Path, targets: np.ndarray, genes: np.ndarray
) -> np.ndarray:
    data = ad.read_h5ad(path)
    if not np.array_equal(np.asarray(data.var_names).astype(str), genes):
        raise ValueError("held truth gene order changed")
    labels = np.asarray(data.obs["target_gene"].astype(str))
    x = data.X.tocsr() if sp.issparse(data.X) else sp.csr_matrix(data.X)
    ctrl_sum = np.asarray(x[labels == "non-targeting"].sum(0)).ravel().astype(np.float64)
    bctrl = np.log1p(TS * ctrl_sum / ctrl_sum.sum())
    rows = []
    for target in targets:
        total = np.asarray(x[labels == target].sum(0)).ravel().astype(np.float64)
        rows.append(np.log1p(TS * total / total.sum()) - bctrl)
    return np.vstack(rows)


def pds_per_target(
    predicted: np.ndarray, real: np.ndarray, score_mask: np.ndarray
) -> np.ndarray:
    p = np.asarray(predicted[:, score_mask], dtype=np.float64)
    r = np.asarray(real[:, score_mask], dtype=np.float64)
    pn = np.linalg.norm(p, axis=1)
    rn = np.linalg.norm(r, axis=1)
    scores = (p @ r.T) / np.maximum(pn[:, None] * rn[None, :], 1e-300)
    scores[pn == 0] = -1.0
    scores[:, rn == 0] = -1.0
    own = np.diag(scores)[:, None]
    ranks = (scores > own).sum(1)
    ties = (scores == own).sum(1) - 1
    return 1.0 - (ranks + 0.5 * ties) / (len(predicted) - 1)


def lfc_from_dh(
    dh: np.ndarray, bctrl: np.ndarray, gate: np.ndarray, ref_lin: np.ndarray
) -> np.ndarray:
    profile = np.expm1(bctrl[None, :] + np.asarray(dh, dtype=np.float64)) / TS
    # Frozen arrays are float32. At genes whose true profile is exactly zero,
    # roundoff in bctrl + dh can land a few ulps below zero after reload.
    np.clip(profile, 0.0, None, out=profile)
    return np.log2(
        (1e6 * profile[:, gate] + 1e-9) / (ref_lin[None, :] + 1e-9)
    )


def weighted_top_precision(
    correct: np.ndarray, confidence: np.ndarray, weight: np.ndarray, frac: float = 0.2
) -> float:
    n = max(1, int(math.ceil(frac * len(correct))))
    order = np.argsort(-confidence, kind="stable")[:n]
    denom = float(weight[order].sum())
    return (
        float((correct[order] * weight[order]).sum() / denom)
        if denom > 0
        else float(np.mean(correct[order]))
    )


def ci(values: np.ndarray, n_boot: int, seed: int) -> dict[str, float]:
    values = np.asarray(values, dtype=np.float64)
    rng = np.random.default_rng(seed)
    boot = np.empty(n_boot, dtype=np.float64)
    for i in range(n_boot):
        boot[i] = values[rng.integers(0, len(values), len(values))].mean()
    return {
        "estimate": float(values.mean()),
        "lo": float(np.percentile(boot, 2.5)),
        "hi": float(np.percentile(boot, 97.5)),
    }


def deserialize_null_models(path: Path) -> list[list[LinearModel]]:
    z = np.load(path)
    out = []
    for draw in range(z["means"].shape[0]):
        folds = []
        for fold in range(N_FOLDS):
            folds.append(
                LinearModel(
                    mean=z["means"][draw, fold].astype(np.float64),
                    scale=z["scales"][draw, fold].astype(np.float64),
                    coef=z["coefs"][draw, fold].astype(np.float64),
                    intercept=float(z["intercepts"][draw, fold]),
                )
            )
        out.append(folds)
    return out


def score_surface(surface: str, max_rss_gib: float, n_boot: int) -> None:
    outdir = WORK / surface
    metrics_path = outdir / "metrics.json"
    if metrics_path.exists():
        raise FileExistsError(f"{metrics_path} already exists; refusing to overwrite")
    z = np.load(outdir / "frozen_predictions.npz", allow_pickle=False)
    genes = z["genes"].astype(str)
    targets = z["targets"].astype(str)
    bctrl = z["bctrl"].astype(np.float64)
    panel_mask = z["panel_mask"].astype(bool)
    dh = {
        "identity": z["dh_identity"].astype(np.float64),
        "agreement": z["dh_agreement"].astype(np.float64),
        "rank8": z["dh_rank8"].astype(np.float64),
        "rank16": z["dh_rank16"].astype(np.float64),
    }
    confidence = {
        "agreement": z["p_agree"].astype(np.float64),
        "rank8": z["p_rank8"].astype(np.float64),
        "rank16": z["p_rank16"].astype(np.float64),
    }

    # Held reference values are first read here, after predictions and null
    # parameters have been frozen on disk.
    ref = RefData.load(surface)
    if not np.array_equal(ref.genes, genes):
        raise ValueError("held reference cache gene order changed")
    rpos = {p: i for i, p in enumerate(ref.perts)}
    rows = np.asarray([rpos[str(t)] for t in targets], dtype=np.int64)
    lfc = {
        name: lfc_from_dh(value, bctrl, ref.gate, ref.ref_lin)
        for name, value in dh.items()
    }
    null_models = deserialize_null_models(outdir / "null_models.npz")

    donors = {name: DonorLite(name) for name in DONOR_POOL}
    string = StringTable()
    final_builder = FeatureBuilder(
        genes,
        np.expm1(bctrl) / np.expm1(bctrl).sum() * 1e6,
        SURFACES[surface]["sources"],
        string,
        donors,
    )
    null_precision = np.zeros((len(null_models), len(targets)), dtype=np.float64)
    per_target: dict[str, list[float]] = {
        "identity_precision": [],
        "agreement_precision": [],
        "rank8_precision": [],
        "rank16_precision": [],
        "identity_credit": [],
        "agreement_credit": [],
        "rank8_credit": [],
        "rank16_credit": [],
    }

    eligible_gate = np.asarray(
        [ref.genes[i] not in PANEL_SET for i in ref.gate_idx], dtype=bool
    )
    used_targets = []
    for out_row, ref_row in enumerate(rows):
        sig = ref.sig[ref_row]
        sig = sig[
            eligible_gate[sig]
            & np.isfinite(ref.lfc[ref_row, sig])
            & (ref.lfc[ref_row, sig] != 0)
        ]
        if len(sig) < 10:
            continue
        full_idx = ref.gate_idx[sig]
        truth = ref.lfc[ref_row, sig]
        base = lfc["identity"][out_row, sig]
        valid = np.isfinite(base) & (base != 0)
        sig, full_idx, truth, base = (
            sig[valid],
            full_idx[valid],
            truth[valid],
            base[valid],
        )
        if len(sig) < 10:
            continue
        correct = (np.sign(base) == np.sign(truth)).astype(np.float64)
        denom = float(np.abs(truth).sum()) or 1.0
        weight = np.abs(base) / denom
        per_target["identity_precision"].append(
            float((correct * weight).sum() / max(weight.sum(), 1e-12))
        )
        for name in ("agreement", "rank8", "rank16"):
            per_target[f"{name}_precision"].append(
                weighted_top_precision(
                    correct, confidence[name][out_row, full_idx], weight
                )
            )
        for name in ("identity", "agreement", "rank8", "rank16"):
            pred = lfc[name][out_row, sig]
            per_target[f"{name}_credit"].append(
                float(np.sign(truth) @ pred) / denom
            )

        feats = final_builder.features(str(targets[out_row]), full_idx, 8)
        if feats is None:
            raise ValueError(f"no final features for {targets[out_row]}")
        _, x8 = feats
        fold = target_fold(str(targets[out_row]))
        for draw, fold_models in enumerate(null_models):
            pnull = fold_models[fold].predict(x8)
            null_precision[draw, len(used_targets)] = weighted_top_precision(
                correct, pnull, weight
            )
        used_targets.append(str(targets[out_row]))

    n_used = len(used_targets)
    if n_used == 0:
        raise ValueError("no held targets passed the DE gate")
    null_precision = null_precision[:, :n_used]
    per_target_np = {k: np.asarray(v, dtype=np.float64) for k, v in per_target.items()}

    real_effect = load_truth_effects(SURFACES[surface]["harness"], targets, genes)
    gallery_excluded = PANEL_SET | set(map(str, targets))
    score_mask = np.asarray([g not in gallery_excluded for g in genes], dtype=bool)
    pds_rows = {
        name: pds_per_target(value, real_effect, score_mask)
        for name, value in dh.items()
    }
    guard(max_rss_gib, "held scoring")

    result: dict[str, object] = {
        "surface": surface,
        "targets_total": int(len(targets)),
        "targets_de_scored": n_used,
        "weighted_sign_precision": {
            "identity_all": float(per_target_np["identity_precision"].mean()),
            "agreement_top20": float(per_target_np["agreement_precision"].mean()),
            "rank8_top20": float(per_target_np["rank8_precision"].mean()),
            "rank16_top20": float(per_target_np["rank16_precision"].mean()),
        },
        "linear_nmae_proxy": {
            name: float(1.0 - per_target_np[f"{name}_credit"].mean())
            for name in ("identity", "agreement", "rank8", "rank16")
        },
        "pds": {name: float(value.mean()) for name, value in pds_rows.items()},
        "bootstrap": {},
        "null_rank8_precision_gain_over_agreement": {
            "values": (
                null_precision.mean(1)
                - per_target_np["agreement_precision"].mean()
            ).tolist()
        },
        "used_targets": used_targets,
        "peak_rss_gib": peak_rss_gib(),
    }

    comparisons = {
        "rank8_precision_minus_agreement": (
            per_target_np["rank8_precision"]
            - per_target_np["agreement_precision"]
        ),
        "rank16_precision_minus_agreement": (
            per_target_np["rank16_precision"]
            - per_target_np["agreement_precision"]
        ),
        "rank8_linear_credit_minus_identity": (
            per_target_np["rank8_credit"] - per_target_np["identity_credit"]
        ),
        "rank16_linear_credit_minus_identity": (
            per_target_np["rank16_credit"] - per_target_np["identity_credit"]
        ),
        "agreement_linear_credit_minus_identity": (
            per_target_np["agreement_credit"] - per_target_np["identity_credit"]
        ),
        "rank8_pds_minus_identity": pds_rows["rank8"] - pds_rows["identity"],
        "rank16_pds_minus_identity": pds_rows["rank16"] - pds_rows["identity"],
        "agreement_pds_minus_identity": (
            pds_rows["agreement"] - pds_rows["identity"]
        ),
    }
    result["bootstrap"] = {
        name: ci(value, n_boot, stable_seed(surface, name))
        for name, value in comparisons.items()
    }
    np.savez_compressed(
        outdir / "per_target_metrics.npz",
        targets=np.asarray(used_targets),
        **per_target_np,
        **{f"pds_{k}": v for k, v in pds_rows.items()},
        null_precision=null_precision,
    )
    json_dump(metrics_path, result)
    print(json.dumps(result, indent=2))


def prepare(max_rss_gib: float) -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    string = StringTable()
    donors = {}
    for name in DONOR_POOL:
        donors[name] = DonorLite(name)
        guard(max_rss_gib, f"prepare donor {name}")
    audit = {
        "string_rows": int(len(string.genes)),
        "string_available": int(string.available.sum()),
        "string_rank": int(string.pcs.shape[1]),
        "donors": {
            name: {
                "perts": len(donor.perts),
                "genes": len(donor.genes),
                "mean_sha256": sha256_array(donor.mean),
            }
            for name, donor in donors.items()
        },
        "panel_genes": len(PANEL),
        "panel_sha256": sha256_lines(PANEL),
        "peak_rss_gib": peak_rss_gib(),
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
    if args.command == "prepare":
        prepare(args.max_rss_gib)
        return
    if not args.surface:
        parser.error("--surface is required for freeze/score")
    if args.command == "freeze":
        freeze_surface(args.surface, args.max_rss_gib, args.n_null)
    else:
        score_surface(args.surface, args.max_rss_gib, args.n_boot)


if __name__ == "__main__":
    main()

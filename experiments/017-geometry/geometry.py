"""Experiment 017 — the geometry of the prediction set, and target-basal modulation.

Two families, both aimed at the one thing that is binding: **directional accuracy**, read
through PDS. Everything downstream (mse, jac, per-perturbation magnitude) has already been
shown to be second-order at our cosine, so the only work worth doing is on the direction.

## Family A — spectral flattening of the panel prediction matrix

This is an inversion of experiment 004's negative result. Truncating the donor SVD lost
*monotonically* (rank 5 -> 0.039, rank 800 -> 0.052, full -> 0.073): the low-variance
directions carry signal. Truncation zeroes them. The opposite operation is to **amplify**
them, by flattening the singular-value spectrum of the prediction matrix:

    P = U S V^T   ->   P_g = U S^(1-gamma) V^T

gamma = 0 is the incumbent, gamma = 1 is full whitening. This is standard practice in
retrieval and in word-embedding post-processing ("all-but-the-top"), and PDS *is* a
retrieval problem: for each prediction, is its own true effect the nearest of 300?

The mechanism it targets is hubness. If every prediction carries a large shared component,
every prediction is close to whichever true effects lie near the panel mean, and
discrimination collapses. Flattening removes the shared component's dominance without
discarding anything -- unlike truncation, it is invertible.

Legality: this is a deterministic transform of our own 300 predictions, computed from our
own predictions. No target information enters.

## Family B — per-gene modulation by the target's basal profile

The control profile is the only context-specific signal the Challenge gives us. Experiment
010 showed it does not pick *donors* (one number per line). This tests it one level down, as
a *per-gene* modulation, which is a much richer use of the same signal:

    predicted effect on gene g  ->  p_target(g)^beta  x  donor effect on gene g

A gene that the target line does not express cannot move much in that line, however loudly
the donor screams about it. Every such gene is pure noise in the cosine.

    uv run python experiments/017-geometry/geometry.py --family A
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, _norm_scale  # noqa: E402

DONORS = {
    "k562_gwps": "data/donors/k562_gwps",
    "orion_hct116": "data/donors/orion_hct116",
    "orion_hek293t": "data/donors/orion_hek293t",
}
TS = 5e4


def pds_cosine(pred: np.ndarray, true: np.ndarray) -> np.ndarray:
    P = pred / (np.linalg.norm(pred, axis=1, keepdims=True) + 1e-12)
    T = true / (np.linalg.norm(true, axis=1, keepdims=True) + 1e-12)
    D = 1.0 - P @ T.T
    n = D.shape[1]
    self_d = np.diag(D)[:, None]
    lt = (D < self_d).sum(1)
    eq = (D == self_d).sum(1)
    return 1.0 - (lt + 0.5 * (eq - 1)) / (n - 1)


def flatten_spectrum(P: np.ndarray, gamma: float) -> np.ndarray:
    """P = U S V^T -> U S^(1-gamma) V^T. gamma=0 identity, gamma=1 full whitening."""
    if gamma == 0.0:
        return P
    U, S, Vt = np.linalg.svd(P, full_matrices=False)
    keep = S > S.max() * 1e-10          # numerical rank; do not resurrect null directions
    return (U[:, keep] * (S[keep] ** (1.0 - gamma))) @ Vt[keep]


def load_banks(treat: str = "gene_norm") -> dict:
    out = {}
    for name, path in DONORS.items():
        d = Donor.load(path)
        m, s = d.decompose()
        if treat == "gene_norm":
            s = (s / _norm_scale(s)).astype(np.float32)
        gl = [str(g) for g in d.genes]
        out[name] = dict(s=s, gpos={g: i for i, g in enumerate(gl)}, gset=set(gl),
                         idx={p: i for i, p in enumerate(d.perts)}, b_ctrl=d.b_ctrl)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", nargs="+", default=["A", "B"])
    ap.add_argument("--weights", default="4:1", help="predictor weights, DONORS order")
    ap.add_argument("--gammas", type=float, nargs="+",
                    default=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0])
    ap.add_argument("--betas", type=float, nargs="+",
                    default=[0.0, 0.1, 0.25, 0.5, 0.75, 1.0])
    ap.add_argument("--cpm-gates", type=float, nargs="+", default=[0.0, 1.0, 5.0, 20.0, 50.0])
    args = ap.parse_args()

    banks = load_banks()
    vcc = list(pd.read_csv("data/challenge/pert_counts.csv")["target_gene"])
    common = [g for g in vcc if all(g in b["idx"] for b in banks.values())]
    print(f"panel: {len(common)} perturbations in all {len(banks)} lines")

    # Pre-build every fold's aligned matrices once; the sweeps are then pure linear algebra.
    folds = {}
    for held in banks:
        others = [n for n in banks if n != held]
        shared = sorted(set.intersection(*[banks[n]["gset"] for n in [held] + others]))
        feat = [g for g in shared if g not in set(common)]
        col = {n: np.fromiter((banks[n]["gpos"][g] for g in feat), np.int64, len(feat))
               for n in [held] + others}
        T = np.vstack([banks[held]["s"][banks[held]["idx"][p]][col[held]] for p in common])
        ws = [float(x) for x in args.weights.split(":")]
        w = dict(zip(others, ws))
        scale = {n: float(np.median(np.linalg.norm(banks[n]["s"][:, col[n]], axis=1))) or 1.0
                 for n in others}
        ref = float(np.median([scale[n] for n in others]))
        P = np.zeros_like(T)
        for n in others:
            P += (w[n] / sum(w.values())) * (ref / scale[n]) * np.vstack(
                [banks[n]["s"][banks[n]["idx"][p]][col[n]] for p in common])
        # Target basal profile on the scored feature axis, as CPM.
        e = np.expm1(banks[held]["b_ctrl"][col[held]])
        cpm_t = 1e6 * e / max(e.sum(), 1e-9)
        folds[held] = dict(P=P.astype(np.float64), T=T.astype(np.float64), cpm=cpm_t)
        print(f"  fold {held:14} {P.shape[0]} perts x {P.shape[1]} scored genes")

    results = {}

    def record(key: str, transform) -> float:
        per = {h: float(pds_cosine(transform(f["P"], f), f["T"]).mean())
               for h, f in folds.items()}
        macro = float(np.mean(list(per.values())))
        results[key] = {**per, "macro": macro}
        folds_s = "  ".join(f"{k.replace('orion_','')[:9]}={v:.4f}" for k, v in per.items())
        print(f"  {key:28} macro PDS {macro:.4f}   [{folds_s}]", flush=True)
        return macro

    if "A" in args.family:
        print("\n== Family A: spectral flattening of the panel prediction matrix ==")
        for g in args.gammas:
            record(f"flatten gamma={g:g}", lambda P, f, g=g: flatten_spectrum(P, g))
        print("\n-- rank-1 removal (all-but-the-top), as a special case --")
        for c in [0.0, 0.25, 0.5, 0.75, 1.0]:
            record(f"demean c={c:g}", lambda P, f, c=c: P - c * P.mean(0, keepdims=True))

    if "B" in args.family:
        print("\n== Family B: per-gene modulation by the TARGET's basal profile ==")
        for b in args.betas:
            record(f"basal beta={b:g}",
                   lambda P, f, b=b: P * (np.maximum(f["cpm"], 1e-3) ** b))
        print("\n-- hard CPM gate on the target's own control profile --")
        for t in args.cpm_gates:
            record(f"cpm gate>{t:g}", lambda P, f, t=t: P * (f["cpm"] > t))

    best = max(results, key=lambda k: results[k]["macro"])
    print(f"\nbest: {best} -> {results[best]['macro']:.4f}")
    Path("experiments/017-geometry/results.json").write_text(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()

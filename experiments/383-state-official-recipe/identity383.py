#!/usr/bin/env python3
"""383 — identity check: does the checkpoint retrieve its own target?

Inputs:  work/truth_val.npz (true per-pert log1p deltas, 50 H1 val perts) and a
         STATE inference output (.h5ad, predicted cells with obs[target_gene]).
Method:  per-pert predicted delta = mean(pred cells) - mean(input control cells)
         (both in model log1p space); cosine-rank retrieval of each TRUE delta
         against all 50 PREDICTED deltas (pds direction: truth ranked vs preds),
         score = mean((n - rank) / (n - 1)) over the 50 truths (midrank ties).
Gate:    >= 0.7 PASS (identity); <= 0.55 STOP per PLAN (no identity, 106's trap).
Also reports: mean cosine(truth, own pred), and the same retrieval with our cap15
H1 design as the prediction set (reference leg, expected high).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"


def midrank_pds(truth: np.ndarray, pred: np.ndarray) -> tuple[float, np.ndarray, float]:
    """truth (P,G), pred (P,G) propagated in the same gene order. Returns
    (mean_score, per-target scores, mean own-cosine)."""
    tn = truth / np.maximum(np.linalg.norm(truth, axis=1, keepdims=True), 1e-12)
    pn = pred / np.maximum(np.linalg.norm(pred, axis=1, keepdims=True), 1e-12)
    sim = tn @ pn.T  # [P, P]: sim[i,j] = cos(truth_i, pred_j)
    p = sim.shape[0]
    own = np.diag(sim)
    scores = np.zeros(p)
    for i in range(p):
        greater = (sim[i] > sim[i, i]).sum()
        equal = (sim[i] == sim[i, i]).sum()
        rank = greater + (equal + 1) / 2  # midrank, 1-based
        scores[i] = (p - rank) / (p - 1)
    return float(scores.mean()), scores, float(own.mean())


def pred_deltas(pred_h5ad: str, pert_col: str = "target_gene",
                ctrl: str = "non-targeting") -> tuple[np.ndarray, np.ndarray]:
    import anndata as ad

    a = ad.read_h5ad(pred_h5ad)
    tg = a.obs[pert_col].astype(str).to_numpy()
    X = a.X.toarray() if not isinstance(a.X, np.ndarray) else a.X
    X = np.asarray(X, dtype=np.float64)
    perts = sorted(set(tg) - {ctrl})
    cm = X[tg == ctrl].mean(0)
    D = np.vstack([X[tg == p].mean(0) - cm for p in perts])
    return np.array(perts), D.astype(np.float32)


def main(pred_h5ad: str, out_name: str) -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    z = np.load(WORK / "truth_val.npz", allow_pickle=False)
    tperts = [str(t) for t in z["perts"]]
    truth = np.asarray(z["delta"], dtype=np.float64)
    pperts, pred = pred_deltas(pred_h5ad)
    common = [t for t in tperts if t in set(pperts)]
    print("common perts:", len(common), "/", len(tperts))
    ti = [tperts.index(t) for t in common]
    pi = [list(pperts).index(t) for t in common]
    score, per_t, own_cos = midrank_pds(truth[ti], np.asarray(pred[pi], dtype=np.float64))
    verdict = "PASS" if score >= 0.7 else ("STOP-NO-IDENTITY" if score <= 0.55 else "GRAY")
    out = {"n_common": len(common), "pds_identity": score,
           "mean_own_cosine": own_cos,
           "per_target": {t: float(s) for t, s in zip(common, per_t)},
           "gate": {"pass_ge": 0.7, "stop_le": 0.55}, "verdict": verdict,
           "pred_file": pred_h5ad}
    (WORK / out_name).write_text(json.dumps(out, indent=2))
    print(json.dumps({k: v for k, v in out.items() if k != "per_target"}, indent=2))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "identity.json")

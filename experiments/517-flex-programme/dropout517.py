#!/usr/bin/env python3
"""517 §3 — per-gene detection (dropout) in Flex controls: is excess dropout a shared probe/chemistry property across the
four Flex contexts (H1, A, B, C), and does it relate to H1's DE signs or DE calls?

Per context, control cells only: detection rate det_g = share of cells with a count > 0; Poisson expectation given the
gene's pooled fraction p_g and each cell's library L_i: E[det_g] = 1 - mean_i exp(-p_g L_i) (10 library bins);
excess dropout x_g = E[det_g] - det_g (> 0: fewer detections than Poisson sampling explains, e.g. a weak probe or bursty
expression). Then: agreement of x_g across contexts (genes >= 5 CPM in both), and on h1-test's truth, whether x_g relates to
(a) the sign of gate-gene lfc, (b) the rate of being called (gate membership) at matched expression.
usage: uv run python experiments/517-flex-programme/dropout517.py
"""
from __future__ import annotations

import csv
import glob
import json
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from stream517 import CSRStream  # noqa: E402


def cat(n):
    return n["categories"][:].astype(str)[n["codes"][:]] if isinstance(n, h5py.Group) else n[:].astype(str)


def excess_dropout(path, ctrl_only=True, nb=10):
    it = CSRStream(path)
    tv = cat(it.f["obs/target_gene"])
    rows = np.flatnonzero(tv == "non-targeting") if ctrl_only else np.arange(it.n)
    lib = it.libraries()
    lab = np.full(it.n, -1, np.int64); lab[rows] = 0
    counts = it.group_sums(lab, 1)[0]
    det = it.group_sums(lab, 1, binary=True)[0] / rows.size
    it.close()
    p = counts / counts.sum()
    L = lib[rows]
    edges = np.quantile(L, np.linspace(0, 1, nb + 1)); b = np.clip(np.searchsorted(edges, L, side="right") - 1, 0, nb - 1)
    Lb = np.array([L[b == k].mean() for k in range(nb)]); nbk = np.bincount(b, minlength=nb) / rows.size
    exp_det = 1.0 - (nbk[None, :] * np.exp(-p[:, None] * Lb[None, :])).sum(1)
    mean_cpm = 1e6 * p
    return dict(x=exp_det - det, det=det, exp_det=exp_det, mean_cpm=mean_cpm, n=int(rows.size))


def main():
    out = {}
    z = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=False)
    genes, keep, targets = z["genes"].astype(str), np.asarray(z["keep"], bool), list(z["targets"].astype(str))
    q_ref, lfc_ref = np.asarray(z["q"], float), np.asarray(z["lfc"], float)
    kg = genes[keep]
    H = excess_dropout(ROOT / "data/harness/h1-test.h5ad")
    official = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
    C = {c: excess_dropout(ROOT / f"data/challenge/context_{c}.h5ad", ctrl_only=True) for c in "ABC"}
    # agreement across contexts on the official axis
    opos = {g: i for i, g in enumerate(official)}
    hcol = np.array([opos.get(g, -1) for g in genes])
    X = {"H1": np.full(len(official), np.nan)}
    okh = (hcol >= 0) & (H["mean_cpm"] >= 5)
    X["H1"][hcol[okh]] = H["x"][okh]
    for c in "ABC":
        X[c] = np.where(C[c]["mean_cpm"] >= 5, C[c]["x"], np.nan)
    agree = {}
    for i, a in enumerate(X):
        for b_ in list(X)[i + 1:]:
            m = np.isfinite(X[a]) & np.isfinite(X[b_])
            agree[f"{a}-{b_}"] = float(np.corrcoef(X[a][m], X[b_][m])[0, 1])
    out["excess_dropout_cross_context_r"] = agree
    for k, d in {"H1": H, **C}.items():
        e = d["mean_cpm"] >= 5
        out[f"{k}_summary"] = dict(n_ctrl=d["n"], median_excess=float(np.median(d["x"][e])),
                                   q90_excess=float(np.quantile(d["x"][e], .9)), frac_excess_gt_0_05=float((d["x"][e] > 0.05).mean()))
    print("excess dropout agreement:", json.dumps({k: round(v, 3) for k, v in agree.items()}), flush=True)
    print("summaries:", json.dumps({k: v for k, v in out.items() if k.endswith("summary")}), flush=True)
    # relation to H1's truth, on the gated axis: sign of gate lfc, and gate membership at matched expression
    xk = H["x"][keep]; mk = H["mean_cpm"][keep]
    own = np.array([[g == t for g in kg] for t in targets])
    gate = (q_ref < 0.05) & np.isfinite(lfc_ref) & (lfc_ref != 0) & ~own
    down_rate = []; call_rate = []
    ebin = np.digitize(np.log10(np.maximum(mk, 1e-3)), np.quantile(np.log10(np.maximum(mk, 1e-3)), [.2, .4, .6, .8]))
    xq = np.digitize(xk, np.quantile(xk, [1 / 3, 2 / 3]))
    tab = {}
    for eb in range(5):
        for xb in range(3):
            m = (ebin == eb) & (xq == xb)
            g = gate[:, m]
            l_ = lfc_ref[:, m]
            n_g = int(g.sum())
            tab[f"expr{eb}_x{xb}"] = dict(genes=int(m.sum()), call_rate=float(g.mean()),
                                         down_frac=float((l_[g] < 0).mean()) if n_g else None, n_gate=n_g)
    out["h1_truth_by_expression_and_excess_dropout"] = tab
    for eb in range(5):
        print(f"expression quintile {eb}: " + "; ".join(
            f"x-tercile {xb}: call rate {tab[f'expr{eb}_x{xb}']['call_rate']:.4f}, down {tab[f'expr{eb}_x{xb}']['down_frac']:.3f}"
            for xb in range(3)), flush=True)
    (HERE / "work/dropout517.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

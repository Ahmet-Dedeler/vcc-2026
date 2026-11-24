#!/usr/bin/env python3
"""528 — transfer check on public cancer-line truths (480's folds; RELATIVE DIRECTION ONLY: other chemistries than Flex,
150 cells per target, ~12k UMI; the Orion truths carry a target-independent technical blob, AGENTS 6.5).

Does excess dropout in a line's OWN controls predict which genes that line's truth calls, at matched expression, as on H1?
Per fold: control features (feat528's rules on the fold's 18,400 controls: log10 mean per-cell CPM `lc`, excess dropout
`x`), the fold's gate genes (480's gate-<fold>.npz: q < .05, finite non-zero lfc, not the target), and 2,200-gene blocks
among non-panel genes ranked by lc (the CPM block), by lc + beta x (beta 22 = H1's trainval-selected value, and a grid), and
by lc + 22 x_shuf (x shuffled within 20 expression bins: the falsifier). Capture = share of gate genes in the block.

Folds: orion_hct116 / orion_hek293t (data/harness/*-thinned.h5ad; controls = its 18,400 non-targeting cells) and k562_gwps
(ReplogleWeissman2022 K562-GWPS; the same 18,400 controls k562harness480.py drew: default_rng(0).choice, streamed from the
dense source in its HDF5 row chunks). KA: each fold's control mean CPM on the gate axis equals the gate file's ctrl_cpm.

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/fold528.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
import h5py  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from feat528 import features  # noqa: E402

GATES = ROOT / "experiments/480-nmae-instrument/work"
K562 = ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
BETA_H1 = 22.0
GRID = [0.0, 5.0, 10.0, 15.0, 22.0, 30.0, 40.0, 60.0]


def k562_features(nb=10):
    h = h5py.File(K562, "r")
    g = h["obs"]["gene"]
    lab = g["categories"][:].astype(str)[g["codes"][:]]
    genes = h["var"]["gene_name"][:].astype(str)
    X = h["X"]
    n, G = X.shape
    rng = np.random.default_rng(0)                         # k562harness480's draw (its first rng call)
    ctrl = np.sort(rng.choice(np.flatnonzero(lab == "non-targeting"), 18400, replace=False))
    counts = np.zeros(G)
    dets = np.zeros(G)
    cpms = np.zeros(G)
    libs = []
    blk = X.chunks[0]
    for lo in range(0, n, blk):
        hi = min(lo + blk, n)
        sel = ctrl[(ctrl >= lo) & (ctrl < hi)] - lo
        if sel.size == 0:
            continue
        D = np.asarray(X[lo:hi, :][sel], dtype=np.float64)
        L = D.sum(1)
        libs.append(L)
        counts += D.sum(0)
        dets += (D > 0).sum(0)
        cpms += (D / np.maximum(L, 1.0)[:, None] * 1e6).sum(0)
    h.close()
    L = np.concatenate(libs)
    m = L.size
    p = counts / counts.sum()
    edges = np.quantile(L, np.linspace(0, 1, nb + 1))
    b = np.clip(np.searchsorted(edges, L, side="right") - 1, 0, nb - 1)
    Lb = np.array([L[b == k].mean() for k in range(nb)])
    nbk = np.bincount(b, minlength=nb) / m
    exp_det = 1.0 - (nbk[None, :] * np.exp(-p[:, None] * Lb[None, :])).sum(1)
    det = dets / m
    mean_cpm = cpms / m
    return dict(genes=genes, pooled_cpm=1e6 * p, mean_cpm=mean_cpm, det=det, exp_det=exp_det, x=exp_det - det,
                keep=mean_cpm > 5.0, n_ctrl=np.int64(m), lib_median=np.float64(np.median(L)))


def block(score, excl, n=2200):
    b = np.zeros(score.size, bool)
    b[np.argsort(-np.where(excl, -np.inf, score), kind="stable")[:n]] = True
    return b


def main():
    out = {}
    for fold in ("orion_hct116", "orion_hek293t", "k562_gwps"):
        t0 = time.time()
        if fold == "k562_gwps":
            F = k562_features()
        else:
            F = features(ROOT / f"data/harness/{fold.replace('_', '-')}-thinned.h5ad")
        np.savez(HERE / f"work/feat528_{fold}.npz", **F)
        z = np.load(GATES / f"gate-{fold}.npz", allow_pickle=True)
        gg = z["genes"].astype(str)
        pos = {g: i for i, g in enumerate(F["genes"])}
        col = np.array([pos[g] for g in gg])
        ka = float(np.max(np.abs(F["mean_cpm"][col] - np.asarray(z["ctrl_cpm"], float)) / np.asarray(z["ctrl_cpm"], float)))
        x = F["x"][col]
        lc = np.log10(np.maximum(F["mean_cpm"][col], 1e-3))
        targets = z["targets"].astype(str)
        q, lfc = np.asarray(z["q"], float), np.asarray(z["lfc"], float)
        own = np.equal.outer(targets, gg)
        gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0) & ~own
        excl = np.isin(gg, targets)
        rs = np.random.default_rng(52802)
        bins = np.digitize(lc, np.quantile(lc, np.linspace(0, 1, 21)[1:-1]))
        xs = x.copy()
        for b_ in np.unique(bins):
            idx = np.flatnonzero(bins == b_)
            xs[idx] = x[rs.permutation(idx)]
        rows = {}

        def cap(name, key):
            bl = block(key, excl)
            per = [gate[i, bl].sum() / gate[i].sum() for i in range(gate.shape[0]) if gate[i].sum() > 0]
            rows[name] = dict(capture=float(gate[:, bl].sum() / gate.sum()), capture_per_target=float(np.mean(per)),
                              swapped=int((bl & ~block(lc, excl)).sum()))
        for b_ in GRID:
            cap(f"beta{b_:g}", lc + b_ * x)
        cap("shuf22", lc + BETA_H1 * xs)
        e = F["keep"]
        # dropout at matched expression: call rate by expression quintile x dropout tercile (517's table)
        eb = np.digitize(lc, np.quantile(lc, [.2, .4, .6, .8]))
        xb = np.digitize(x, np.quantile(x, [1 / 3, 2 / 3]))
        tab = {f"e{a}": [float(gate[:, (eb == a) & (xb == c)].mean()) for c in range(3)] for a in range(5)}
        out[fold] = dict(ka_ctrl_cpm_max_rel=ka, n_ctrl=int(F["n_ctrl"]), lib_median=float(F["lib_median"]),
                         n_targets=int(gate.shape[0]), gate_per_target_median=float(np.median(gate.sum(1))),
                         x_median_expressed=float(np.median(F["x"][e])), frac_x_gt_0_05=float((F["x"][e] > 0.05).mean()),
                         rows=rows, call_rate_by_expr_quintile_x_tercile=tab, seconds=time.time() - t0)
        b0 = rows["beta0"]["capture"]
        print(f"{fold}: KA ctrl_cpm max rel {ka:.2e}; x median {out[fold]['x_median_expressed']:.4f}; gate/target median "
              f"{out[fold]['gate_per_target_median']:.0f}", flush=True)
        print("   " + "  ".join(f"{k} {v['capture']:.4f} ({v['capture'] - b0:+.4f}, sw {v['swapped']})" for k, v in rows.items()),
              flush=True)
        print("   call rate by expression quintile [x tercile lo/mid/hi]: " +
              "; ".join(f"{k} {v[0]:.4f}/{v[1]:.4f}/{v[2]:.4f}" for k, v in tab.items()), flush=True)
    (HERE / "work/fold528.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

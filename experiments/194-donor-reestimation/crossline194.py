"""194 — estimator-level geometry: cross-line agreement, the target-independent blob,
and the long-gene technical axis.

Three statistics per estimator, all computed on the SAME shared target set and the SAME
intersected gene axis so that only the estimator differs:

`same` / `diff`   mean cosine between two lines' rows for the SAME target, and for
                  mismatched targets, in the residual (per-table mean removed) space.
                  159 §3.1: attenuation multiplies both identically, so `same` rising
                  faster than `diff` is a real recovery of shared signal and not a
                  rescaling.  This is the quantity a better donor estimator has to move.

`blob`            share of total row energy carried by the single across-target mean
                  vector.  A target-independent component is the metric's own zero
                  (003) and it is what the harness warning in AGENTS.md 6.4 is about;
                  the `debias` variant predicts, in advance, that it should shrink,
                  because a concave `log1p` on a 10^2-cell target and a 10^5-cell
                  control leaves exactly such a common offset.

`longgene`        cosine between the table's top principal direction and the incumbent
                  table's top principal direction, i.e. does the arm still load the
                  X-Atlas long-gene / nuclear-RNA technical axis.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
SEED = 20260903

PANEL = {l.split(",")[0] for l in
         (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()}


def load(path: Path) -> tuple[np.ndarray, list[str], list[str]]:
    d = np.load(path / "delta.npy", mmap_mode="r")
    perts = (path / "perts.txt").read_text().split("\n")
    genes = (path / "genes.txt").read_text().split("\n")
    return d, perts, genes


def pair_geometry(pa: Path, pb: Path, max_targets: int = 3000) -> dict:
    da, ta, ga = load(pa)
    db, tb, gb = load(pb)
    shared_t = sorted(set(ta) & set(tb) - PANEL)
    rng = np.random.default_rng(SEED)
    if len(shared_t) > max_targets:
        shared_t = sorted(rng.choice(np.array(shared_t), max_targets, replace=False))
    pos_a = {g: i for i, g in enumerate(ga)}
    pos_b = {g: i for i, g in enumerate(gb)}
    shared_g = [g for g in ga if g and g in pos_b and g not in PANEL]
    ca = np.array([pos_a[g] for g in shared_g])
    cb = np.array([pos_b[g] for g in shared_g])
    ra = np.array([ta.index(t) for t in shared_t])
    rb = np.array([tb.index(t) for t in shared_t])
    A = np.asarray(da[ra][:, ca], dtype=np.float64)
    B = np.asarray(db[rb][:, cb], dtype=np.float64)
    A -= A.mean(0)
    B -= B.mean(0)
    An = A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)
    Bn = B / np.maximum(np.linalg.norm(B, axis=1, keepdims=True), 1e-300)
    S = An @ Bn.T
    n = len(shared_t)
    same = float(np.mean(np.diag(S)))
    off = ~np.eye(n, dtype=bool)
    diff = float(S[off].mean())
    return dict(n_targets=n, n_genes=len(shared_g), same=same, diff=diff,
                excess=same - diff, ratio=same / diff if diff else float("nan"),
                same_sd=float(np.std(np.diag(S), ddof=1)))


def table_geometry(path: Path, ref: Path | None = None,
                   max_targets: int = 3000) -> dict:
    d, perts, genes = load(path)
    keep = [i for i, t in enumerate(perts) if t not in PANEL]
    rng = np.random.default_rng(SEED)
    if len(keep) > max_targets:
        keep = sorted(rng.choice(np.array(keep), max_targets, replace=False))
    cols = np.array([i for i, g in enumerate(genes) if g and g not in PANEL])
    M = np.asarray(d[np.array(keep)][:, cols], dtype=np.float64)
    total = float((M * M).sum())
    mu = M.mean(0)
    blob = float(len(M) * (mu @ mu) / total)
    C = M - mu
    # top principal direction of the centred table
    u, s, vt = np.linalg.svd(C, full_matrices=False)
    top = vt[0]
    out = dict(n_targets=len(keep), n_genes=len(cols), blob_share=blob,
               pc1_share=float(s[0] ** 2 / (s ** 2).sum()),
               median_norm=float(np.median(np.linalg.norm(M, axis=1))))
    if ref is not None:
        dr, pr, gr = load(ref)
        posr = {g: i for i, g in enumerate(gr)}
        colr = np.array([posr[genes[c]] for c in cols])
        rowr = np.array([pr.index(perts[i]) for i in keep])
        R = np.asarray(dr[rowr][:, colr], dtype=np.float64)
        R -= R.mean(0)
        _, _, vtr = np.linalg.svd(R, full_matrices=False)
        out["pc1_cos_to_incumbent"] = float(abs(top @ vtr[0]))
        out["row_cos_to_incumbent"] = float(np.median(
            np.sum(C * R, 1) / np.maximum(np.linalg.norm(C, axis=1) *
                                          np.linalg.norm(R, axis=1), 1e-300)))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", nargs="*", required=True)
    ap.add_argument("--out", default=str(WORK / "crossline.json"))
    args = ap.parse_args()
    root = WORK / "tables"
    rows = []
    for v in args.variants:
        r = dict(variant=v)
        for a, b in (("orion_hct116", "orion_hek293t"),
                     ("k562_gwps", "orion_hct116"),
                     ("k562_gwps", "orion_hek293t")):
            r[f"{a[:4]}_{b[:8]}"] = pair_geometry(root / a / v, root / b / v)
        for d in ("k562_gwps", "orion_hct116", "orion_hek293t"):
            r[f"geom_{d}"] = table_geometry(root / d / v, ROOT / "data/donors" / d)
        rows.append(r)
        print(json.dumps(r, indent=2), flush=True)
    Path(args.out).write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()

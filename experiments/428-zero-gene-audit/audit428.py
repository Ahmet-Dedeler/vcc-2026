"""428: do our A/B/C packages carry genes that a perturbation's 400 cells never express (delivered arithmetic-mean
lfc = log2(1e-9 / ctrl) ~ -35)? Such a gene in the reference's gate costs ~33 in that target's nmae numerator.
Streams a prediction h5ad block by block; per (context, target): genes with control mean CPM > 5 (the scorer's gate),
target gene excluded, whose block total is 0, plus the delivered lfc distribution tail.
usage: uv run python audit428.py <pred.h5ad> [--bundle data/challenge]"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def str_col(h, path):
    obj = h[path]
    if isinstance(obj, h5py.Group):
        if "categories" in obj:
            return np.asarray(obj["categories"][:]).astype(str)[np.asarray(obj["codes"][:])]
        return np.asarray(obj["values"][:]).astype(str)
    return np.asarray(obj[:]).astype(str)


def ctrl_cpm(bundle, ctx, n_genes):
    with h5py.File(bundle / f"context_{ctx}.h5ad", "r") as h:
        ind = np.asarray(h["X"]["indices"][:], dtype=np.int64); dat = np.asarray(h["X"]["data"][:], dtype=np.float64)
        ip = np.asarray(h["X"]["indptr"][:]); lib = np.add.reduceat(dat, ip[:-1]); rid = np.repeat(np.arange(ip.size - 1), np.diff(ip))
        return np.bincount(ind, weights=dat / lib[rid] * 1e6, minlength=n_genes) / (ip.size - 1)


ap = argparse.ArgumentParser(); ap.add_argument("pred"); ap.add_argument("--bundle", default=str(ROOT / "data/challenge"))
a = ap.parse_args()
with h5py.File(a.pred, "r") as h:
    genes = str_col(h, "var/_index"); ctxs = str_col(h, "obs/context"); tg = str_col(h, "obs/target_gene")
    ip = np.asarray(h["X"]["indptr"][:])
    gpos = {g: i for i, g in enumerate(genes)}
    ctrl = {c: ctrl_cpm(Path(a.bundle), c, genes.size) for c in sorted(set(ctxs))}
    out, start = [], 0
    n = len(tg)
    while start < n:
        end = start
        while end < n and tg[end] == tg[start] and ctxs[end] == ctxs[start]:
            end += 1
        lo, hi = int(ip[start]), int(ip[end])
        ind = np.asarray(h["X"]["indices"][lo:hi], dtype=np.int64); dat = np.asarray(h["X"]["data"][lo:hi], dtype=np.float64)
        rows = np.repeat(np.arange(end - start), np.diff(ip[start:end + 1]))
        lib = np.bincount(rows, weights=dat, minlength=end - start)
        tot = np.bincount(ind, weights=dat, minlength=genes.size)
        mean_cpm = np.bincount(ind, weights=dat / lib[rows] * 1e6, minlength=genes.size) / (end - start)
        c = ctrl[ctxs[start]]; gate = c > 5.0
        t = gpos.get(tg[start], -1)
        if t >= 0:
            gate = gate.copy(); gate[t] = False
        lfc = np.log2((mean_cpm + 1e-9) / (c + 1e-9))
        zero = gate & (tot == 0)
        out.append(dict(context=str(ctxs[start]), target=str(tg[start]), n_zero=int(zero.sum()),
                        zero_genes=[str(genes[j]) for j in np.flatnonzero(zero)][:10],
                        n_lfc_below_m5=int((gate & (lfc < -5)).sum()), min_lfc=float(lfc[gate].min()),
                        n_lfc_above_p3=int((gate & (lfc > 3)).sum()), max_lfc=float(lfc[gate].max()),
                        q999_abs_lfc=float(np.quantile(np.abs(lfc[gate]), 0.999))))
        start = end
tot_zero = sum(r["n_zero"] for r in out)
print(f"blocks {len(out)}; gated genes with ZERO block total (target excluded): {tot_zero} in "
      f"{sum(r['n_zero'] > 0 for r in out)} blocks; lfc < -5: {sum(r['n_lfc_below_m5'] for r in out)}; "
      f"lfc > +3: {sum(r['n_lfc_above_p3'] for r in out)} (max {max(r['max_lfc'] for r in out):+.2f})")
for r in sorted(out, key=lambda r: -r["n_zero"])[:15]:
    if r["n_zero"] or r["n_lfc_below_m5"]:
        print(r)
Path(a.pred).with_suffix(".audit428.json").write_text(json.dumps(out, indent=1))

"""449: list gated (context, target, gene) cells with delivered arithmetic-mean lfc < -4.5 or > +3 in a restored package,
with the target's 410 key and weight. usage: uv run python audit_extremes.py <pred.h5ad>"""
import sys
import h5py
import numpy as np


def str_col(h, path):
    obj = h[path]
    if isinstance(obj, h5py.Group):
        if "categories" in obj:
            return np.asarray(obj["categories"][:]).astype(str)[np.asarray(obj["codes"][:])]
        return np.asarray(obj["values"][:]).astype(str)
    return np.asarray(obj[:]).astype(str)


def ctrl_cpm(ctx, n):
    with h5py.File(f"data/challenge/context_{ctx}.h5ad", "r") as h:
        ind = np.asarray(h["X"]["indices"][:], np.int64); dat = np.asarray(h["X"]["data"][:], np.float64)
        ip = np.asarray(h["X"]["indptr"][:]); lib = np.add.reduceat(dat, ip[:-1]); rid = np.repeat(np.arange(ip.size - 1), np.diff(ip))
        return np.bincount(ind, weights=dat / lib[rid] * 1e6, minlength=n) / (ip.size - 1)


k = np.load("experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi1e+09.npz", allow_pickle=True)
key = dict(zip(map(str, k["targets"]), k["key"])); w = dict(zip(map(str, k["targets"]), k["w"]))
rows = []
with h5py.File(sys.argv[1], "r") as h:
    genes = str_col(h, "var/_index"); ctxs = str_col(h, "obs/context"); tg = str_col(h, "obs/target_gene")
    ip = np.asarray(h["X"]["indptr"][:]); gpos = {g: i for i, g in enumerate(genes)}
    ctrl = {c: ctrl_cpm(c, genes.size) for c in sorted(set(ctxs))}
    start, n = 0, len(tg)
    while start < n:
        end = start
        while end < n and tg[end] == tg[start] and ctxs[end] == ctxs[start]:
            end += 1
        lo, hi = int(ip[start]), int(ip[end])
        ind = np.asarray(h["X"]["indices"][lo:hi], np.int64); dat = np.asarray(h["X"]["data"][lo:hi], np.float64)
        r = np.repeat(np.arange(end - start), np.diff(ip[start:end + 1])); lib = np.bincount(r, weights=dat, minlength=end - start)
        mean_cpm = np.bincount(ind, weights=dat / lib[r] * 1e6, minlength=genes.size) / (end - start)
        c = ctrl[ctxs[start]]; gate = c > 5.0; t = gpos.get(tg[start], -1)
        if t >= 0:
            gate[t] = False
        lfc = np.log2((mean_cpm + 1e-9) / (c + 1e-9))
        for j in np.flatnonzero(gate & ((lfc < -4.5) | (lfc > 3))):
            rows.append((ctxs[start], tg[start], genes[j], float(lfc[j]), float(c[j]), key.get(tg[start], np.nan), w.get(tg[start], np.nan)))
        start = end
rows.sort(key=lambda r: r[3])
print(len(rows), "extreme gated cells (lfc < -4.5 or > +3)")
for r in rows:
    print(f"  {r[0]} {r[1]:9} {r[2]:10} lfc {r[3]:+6.2f} ctrlCPM {r[4]:8.1f}  key {r[5]:+.4f} w {r[6]:.2f}")

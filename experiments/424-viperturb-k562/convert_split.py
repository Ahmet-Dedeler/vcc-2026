"""424 diagnostic — split-half pseudobulks from the VIPerturb bins, to read VIPerturb's own replicate reliability.

For a fixed target set (the 284 panel targets, the 300 strongest K562-GWPS targets shared with VIPerturb, and 700 random
shared targets; seed 424) and every bin, sums raw counts over:
  parity halves   cells split by crc32(barcode) & 1 (a random half of each target's cells, and of NO-TARGET)
  guide halves    the target's most frequent sgRNA vs its other sgRNAs (targets with >= 2 guides)
NO-TARGET parity halves come from binA only (the three bins carry the identical control cells).
Writes work/split/<bin>.npz (sums float64 [groups x features], group labels) with the same streaming reader and audits as
convert_bin.py (per-cell colSums == nCount_RNA).

    uv run python convert_split.py binA binB binC
"""
from __future__ import annotations

import json
import resource
import shutil
import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from convert_bin import column, dimnames, named, names_of
from rds_stream import RDSStream

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"


def target_set():
    cache = WORK / "split/targets.json"
    if cache.exists():
        return json.loads(cache.read_text())
    vp = [l for l in (ROOT / "data/donors/viperturb_k562/perts.txt").read_text().splitlines() if l]
    gp = [l for l in (ROOT / "data/donors/k562_gwps/perts.txt").read_text().splitlines() if l]
    panel = [l.split(",")[0] for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()]
    sh = sorted(set(vp) & set(gp))
    Dg = np.load(ROOT / "data/donors/k562_gwps/delta.npy", mmap_mode="r")
    ig = {p: i for i, p in enumerate(gp)}
    nrm = np.asarray([float(np.linalg.norm(np.asarray(Dg[ig[p]], dtype=np.float64))) for p in sh])
    strong = [sh[i] for i in np.argsort(-nrm)[:300]]
    rest = sorted(set(sh) - set(strong) - set(panel))
    rnd = sorted(np.random.default_rng(424).choice(rest, 700, replace=False).tolist())
    out = dict(panel=sorted(set(panel) & set(vp)), strong=strong, random=rnd)
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(out))
    return out


def main(b):
    t0 = time.time()
    rss = lambda: resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30
    log = lambda *a: print(f"[split {b} {time.time() - t0:5.0f}s rss {rss():.2f}G]", *a, flush=True)
    ts = target_set()
    sel = sorted(set(ts["panel"]) | set(ts["strong"]) | set(ts["random"]))
    spill = WORK / "spill" / f"split_{b}"
    rd = RDSStream(WORK / "dl" / f"genome_wide_{b}.RDS", spill, log=log)
    hdr, obj, eof = rd.read()
    assert eof
    rna = named(obj.attrs["assays"], "RNA")
    cnt = named(rna.attrs["layers"], "counts")
    features = dimnames(rna.attrs["features"])[0]
    meta = obj.attrs["meta.data"]
    cells = list(meta.attrs["row.names"])
    gene = column(meta, "gene").astype(str)
    guide = column(meta, "guide").astype(str)
    ncount = np.asarray(column(meta, "nCount_RNA"), dtype=np.float64)
    p = np.asarray(cnt.attrs["p"], dtype=np.int64)
    ispill, xspill = cnt.attrs["i"], cnt.attrs["x"]
    labels = [f"{t}|{h}" for t in sel for h in ("p0", "p1", "gA", "gB")] + ["NO-TARGET|p0", "NO-TARGET|p1"]
    li = {l: k for k, l in enumerate(labels)}
    par = np.asarray([zlib.crc32(c.encode()) & 1 for c in cells])
    cpar = np.full(len(cells), -1)
    cgui = np.full(len(cells), -1)
    df = pd.DataFrame(dict(gene=gene, guide=guide))
    top = df[df.gene.isin(sel)].groupby("gene").guide.agg(lambda s: s.value_counts().index[0]).to_dict()
    nguide = df[df.gene.isin(sel)].groupby("gene").guide.nunique().to_dict()
    for k, (g, gd) in enumerate(zip(gene, guide)):
        if g in top:
            cpar[k] = li[f"{g}|p{par[k]}"]
            if nguide[g] >= 2:
                cgui[k] = li[f"{g}|{'gA' if gd == top[g] else 'gB'}"]
        elif g == "NO-TARGET" and b == "binA":
            cpar[k] = li[f"NO-TARGET|p{par[k]}"]
    G, F = len(labels), len(features)
    acc = np.zeros((G, F))
    lib = np.zeros(len(cells))
    with open(ispill.path, "rb") as fi, open(xspill.path, "rb") as fx:
        for c0 in range(0, len(cells), 4000):
            c1 = min(c0 + 4000, len(cells))
            a, z = int(p[c0]), int(p[c1])
            fi.seek(4 * a); fx.seek(4 * a)
            ii = np.fromfile(fi, dtype=np.int32, count=z - a)
            xx = np.fromfile(fx, dtype=np.float32, count=z - a).astype(np.float64)
            local = np.repeat(np.arange(c1 - c0), np.diff(p[c0:c1 + 1]))
            lib[c0:c1] = np.bincount(local, weights=xx, minlength=c1 - c0)
            for codes in (cpar, cgui):
                cc = codes[c0:c1][local]
                keep = cc >= 0
                if keep.any():
                    m = sp.csr_matrix((xx[keep], (cc[keep], ii[keep])), shape=(G, F)).tocoo()
                    acc[m.row, m.col] += m.data
    assert float(np.abs(lib - ncount).max()) == 0.0, "colSums != nCount_RNA"
    ncell = np.bincount(cpar[cpar >= 0], minlength=G) + np.bincount(cgui[cgui >= 0], minlength=G)
    np.savez(WORK / "split" / f"{b}.npz", sums=acc, labels=np.asarray(labels), ncell=ncell, features=np.asarray(features))
    shutil.rmtree(spill)
    log(f"wrote split/{b}.npz: {G} groups ({len(top)} targets in this bin); peak RSS {rss():.2f} GiB")


if __name__ == "__main__":
    for b in sys.argv[1:]:
        main(b)

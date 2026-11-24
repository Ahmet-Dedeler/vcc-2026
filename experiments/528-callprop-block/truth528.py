#!/usr/bin/env python3
"""528 — a leak-free training truth for the call-propensity key: the scorer's DE calls for h1-trainval's 126 perturbations
(disjoint from h1-test's 85 targets; same H1 Flex screen), each of 400 cells tested against h1-test's own 18,400 control
cells (the reference's control set, so the calls have h1-test's power and threshold).

Scorer rules (exact_de323): per-cell CPM = counts / (library / 1e6); Wilcoxon rank-sum z of log1p CPM against the sorted
controls (production's kernel, rc_emit.wilcoxon_z, scanpy tie correction); p = 2 sf(|z|); BH per perturbation over the
5-CPM universe; lfc = log2 of the arithmetic mean CPM ratio (+1e-9).

Memory-bounded: the sorted control matrix is built one gene chunk at a time (6 chunks of ~1,800 universe genes, streamed
from the CSR file), and the perturbed cells' raw counts for the chunk are held as float32 (exact for integer counts).
Known answer first: the same code on 6 of h1-test's own targets must reproduce 323's reference q / lfc.

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/truth528.py
Writes work/truth528_trainval.npz: perts, genes (the h1-test 5-CPM universe, reference order), p, q, lfc (float64).
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
os.environ["NUMBA_CACHE_DIR"] = str(HERE / "work" / "numba-cache")
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT / "src"))

import h5py  # noqa: E402
import numpy as np  # noqa: E402
from scipy.stats import false_discovery_control, norm  # noqa: E402

from vcc2026 import rc_emit  # noqa: E402

TEST = ROOT / "data/harness/h1-test.h5ad"
TRAIN = ROOT / "data/harness/h1-trainval.h5ad"
REF = ROOT / "experiments/323-current-family-validation/work/reference.npz"
N_CHUNK = 6
ROWS = 2000


def labels(f):
    o = f["obs/target_gene"]
    if isinstance(o, h5py.Group):
        return o["categories"][:].astype(str)[o["codes"][:]] if "categories" in o else o["values"][:].astype(str)
    return o[:].astype(str)


class Rows:
    """row access to a CSR h5ad: full-row libraries, and a dense float32 (rows x cols) block of selected columns."""

    def __init__(self, path):
        self.f = h5py.File(path, "r")
        self.ipd = self.f["X/indptr"][:].astype(np.int64)
        self.n = self.ipd.size - 1
        self.G = int(self.f["X"].attrs["shape"][1])
        self.lab = labels(self.f)

    def libraries(self):
        out = np.zeros(self.n)
        for s in range(0, self.n, ROWS):
            e = min(self.n, s + ROWS)
            a, b = int(self.ipd[s]), int(self.ipd[e])
            d = self.f["X/data"][a:b].astype(np.float64)
            ip = self.ipd[s:e + 1] - a
            out[s:e] = np.add.reduceat(d, ip[:-1]) if d.size else 0.0
            empty = ip[1:] == ip[:-1]
            out[s:e][empty] = 0.0
        return out

    def collect(self, dest, n_out, colmap, ncol):
        """raw counts of every row r with dest[r] >= 0, scattered to out[dest[r]], in the columns colmap[g] >= 0, as
        float32 (exact for integer counts); one sequential pass over the file in 2,000-row chunks."""
        out = np.zeros((n_out, ncol), np.float32)
        for s in range(0, self.n, ROWS):
            e = min(self.n, s + ROWS)
            if not (dest[s:e] >= 0).any():
                continue
            a, b = int(self.ipd[s]), int(self.ipd[e])
            d = self.f["X/data"][a:b]
            ix = self.f["X/indices"][a:b]
            rr = np.repeat(np.arange(s, e), np.diff(self.ipd[s:e + 1]))
            dr = dest[rr]
            cc = colmap[ix]
            ok = (dr >= 0) & (cc >= 0)
            out[dr[ok], cc[ok]] = d[ok]
        return out


def main():
    t0 = time.time()
    ref = np.load(REF, allow_pickle=False)
    genes = ref["genes"].astype(str)
    keep = np.asarray(ref["keep"], bool)
    ref_lin_ref = np.asarray(ref["ref_lin"], float)
    gidx = np.flatnonzero(keep)
    n_keep = gidx.size
    Te, Tr = Rows(TEST), Rows(TRAIN)
    assert Te.G == Tr.G == genes.size
    ctrl_rows = np.flatnonzero(Te.lab == "non-targeting")
    assert ctrl_rows.size == 18400
    lib_te, lib_tr = Te.libraries(), Tr.libraries()
    lib_c = lib_te[ctrl_rows]
    # KA targets: 6 of h1-test's own (reference order), then the 126 trainval perturbations
    ka_targets = list(ref["targets"].astype(str))[::15][:6]
    perts = sorted(set(Tr.lab) - {"non-targeting"})
    jobs = [("KA", t, np.flatnonzero(Te.lab == t)) for t in ka_targets] + [("TR", p, np.flatnonzero(Tr.lab == p)) for p in perts]
    assert all(r.size == 400 for _, _, r in jobs), "not 400 cells per perturbation"
    # destination slots: controls in file order; each job's 400 cells in file order, jobs back to back per source file
    dest_ctrl = np.full(Te.n, -1, np.int64)
    dest_ctrl[ctrl_rows] = np.arange(ctrl_rows.size)
    dest_ka = np.full(Te.n, -1, np.int64)
    dest_tr = np.full(Tr.n, -1, np.int64)
    slot, n_ka_cells, n_tr_cells = [], 0, 0
    for src, name, rows in jobs:
        if src == "KA":
            dest_ka[rows] = n_ka_cells + np.arange(rows.size)
            slot.append(n_ka_cells)
            n_ka_cells += rows.size
        else:
            dest_tr[rows] = n_tr_cells + np.arange(rows.size)
            slot.append(n_tr_cells)
            n_tr_cells += rows.size
    Z = np.zeros((len(jobs), n_keep))
    LIN = np.zeros((len(jobs), n_keep))
    ref_lin = np.zeros(n_keep)
    chunks = np.array_split(np.arange(n_keep), N_CHUNK)
    for ci, ch in enumerate(chunks):
        colmap = np.full(genes.size, -1, np.int64)
        colmap[gidx[ch]] = np.arange(ch.size)
        C = Te.collect(dest_ctrl, ctrl_rows.size, colmap, ch.size).astype(np.float64)   # (18400, chunk) raw counts
        scale = lib_c / 1e6
        scale[scale == 0] = 1.0
        # the scorer's ref_lin: mean of counts / (lib / 1e6), accumulated in 1,000-row chunks like exact_de323
        acc = np.zeros(ch.size)
        for lo in range(0, C.shape[0], 1000):
            acc += (C[lo:lo + 1000] / scale[lo:lo + 1000, None]).sum(0)
        ref_lin[ch] = acc / C.shape[0]
        cs = np.ascontiguousarray(np.log1p(C / scale[:, None]).T)
        del C
        cs.sort(axis=1)
        tie = rc_emit.ctrl_tiesums(cs)
        XK = Te.collect(dest_ka, n_ka_cells, colmap, ch.size)
        XT = Tr.collect(dest_tr, n_tr_cells, colmap, ch.size)
        for j, (src, name, rows) in enumerate(jobs):
            X = (XK if src == "KA" else XT)[slot[j]:slot[j] + 400].astype(np.float64)
            lib = (lib_te if src == "KA" else lib_tr)[rows]
            sc = lib / 1e6
            sc[sc == 0] = 1.0
            lin = X / sc[:, None]
            Z[j, ch] = rc_emit.wilcoxon_z(cs, tie, np.ascontiguousarray(np.log1p(lin).T))
            LIN[j, ch] = lin.mean(0)
        del XK, XT
        del cs, tie
        print(f"chunk {ci + 1}/{N_CHUNK} ({ch.size} genes) {time.time() - t0:.0f}s", flush=True)
    P = 2 * norm.sf(np.abs(Z))
    Q = np.array([false_discovery_control(p, method="bh") for p in P])
    LFC = np.log2((LIN + 1e-9) / (ref_lin + 1e-9))
    # KA
    tpos = {t: i for i, t in enumerate(ref["targets"].astype(str))}
    rq, rl = np.asarray(ref["q"], float), np.asarray(ref["lfc"], float)
    ka = dict(ref_lin_max_rel=float(np.max(np.abs(ref_lin - ref_lin_ref) / ref_lin_ref)))
    for j, t in enumerate(ka_targets):
        i = tpos[t]
        ka[t] = dict(q_maxabs=float(np.nanmax(np.abs(Q[j] - rq[i]))), lfc_maxabs=float(np.nanmax(np.abs(LFC[j] - rl[i]))),
                     calls_equal=bool(np.array_equal(Q[j] < 0.05, rq[i] < 0.05)), n_calls=int((Q[j] < 0.05).sum()))
    print("KA (h1-test targets recomputed):", json.dumps(ka), flush=True)
    ok = all(v["calls_equal"] and v["q_maxabs"] < 1e-9 and v["lfc_maxabs"] < 1e-9 for k, v in ka.items() if k != "ref_lin_max_rel")
    n0 = len(ka_targets)
    np.savez(HERE / "work/truth528_trainval.npz", perts=np.array(perts), genes=genes[keep], p=P[n0:], q=Q[n0:], lfc=LFC[n0:])
    nr = (Q[n0:] < 0.05).sum(1)
    summ = dict(ka=ka, ka_ok=ok, n_perts=len(perts), n_real_quantiles=np.percentile(nr, [0, 10, 25, 50, 75, 90, 100]).tolist(),
                seconds=time.time() - t0)
    (HERE / "work/truth528.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps({k: v for k, v in summ.items() if k != "ka"}), flush=True)
    if not ok:
        raise SystemExit("KA FAILED")


if __name__ == "__main__":
    main()

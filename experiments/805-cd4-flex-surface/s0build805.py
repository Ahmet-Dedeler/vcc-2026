"""805 stage 0 (EC2 only): select cells by PREREGISTRATION-stage0.md and write the P and S truths at Arc-like depth.

    uv run python experiments/805-cd4-flex-surface/s0build805.py select  <condition> <srcdir> <outdir>
    uv run python experiments/805-cd4-flex-surface/s0build805.py extract <condition> <srcdir> <outdir> <donor>   (x4, parallel)
    uv run python experiments/805-cd4-flex-surface/s0build805.py merge   <condition> <srcdir> <outdir>

<srcdir> holds D{1..4}_<condition>.assigned_guide.h5ad (downloaded from the public bucket by the job). Writes
<outdir>/truth_P.h5ad, truth_S.h5ad (anndata: obs.target_gene with "non-targeting", var index = the official 18,533
genes, X CSR float32 integer counts), targets_P.txt, targets_S.txt and select.json (every number the prereg asks to be
recorded). Nothing about expression is computed here beyond the per-cell depth needed for thinning.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
SUP = ROOT / "data/external/cd4_cell2026/suppl_tables_2026-08-24"
DONORS = ("D1", "D2", "D3", "D4")
N_CTRL, N_CELLS, N_P, SEED, DEPTH = 18400, 400, 100, 805, 20000.0
if os.environ.get("T805_MOCK"):          # the mock test only (test805_mock.py): tiny sizes
    N_CTRL, N_CELLS, N_P, DEPTH = 40, 10, 5, 50.0
CTRL = "non-targeting"


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def col(h, path):
    """An obs/var column as a numpy array of str (categorical, nullable-string or plain), or bool/num as is."""
    node = h[path]
    if isinstance(node, h5py.Group):
        if "categories" in node:
            cats = node["categories"][:]
            cats = np.array([c.decode() if isinstance(c, bytes) else str(c) for c in cats], dtype=object)
            codes = node["codes"][:]
            out = np.where(codes >= 0, cats[np.clip(codes, 0, None)], "")
            return out
        if "values" in node:
            return node["values"].asstr()[:]
        raise ValueError(path)
    x = node[:]
    if x.dtype.kind in "SO":
        return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in x], dtype=object)
    return x


def largest_remainder(weights: np.ndarray, total: int) -> np.ndarray:
    q = total * weights / weights.sum()
    n = np.floor(q).astype(np.int64)
    n[np.argsort(-(q - n), kind="stable")[: total - int(n.sum())]] += 1
    return n


def select(cond: str, src: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    genes = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).to_numpy()
    gpos = {g: i for i, g in enumerate(genes)}
    ess = {ln.strip() for ln in open(ROOT / "data/donors/k562_essential/perts.txt") if ln.strip()}
    abc = pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str).tolist()
    de = pd.read_csv(SUP / "DE_stats.suppl_table.csv")
    de = de[de["culture_condition"] == cond]
    assert de["target_contrast_gene_name"].is_unique and de["target_contrast"].is_unique
    sg = pd.read_csv(SUP / "sgrna_library_metadata.suppl_table.csv")
    sg_tid = dict(zip(sg["sgRNA"].astype(str), sg["target_gene_id"].astype(str)))
    kd = pd.read_csv(SUP / "guide_kd_efficiency.suppl_table.csv")
    kd = kd[kd["culture_condition"] == cond]
    kd_ok = set(kd.loc[kd["signif_knockdown"].astype(bool), kd.columns[0]].astype(str))

    # ---- obs of the four files
    meta = {}
    for d in DONORS:
        f = src / f"{d}_{cond}.assigned_guide.h5ad"
        with h5py.File(f, "r") as h:
            gid = col(h, "obs/guide_id").astype(str)
            gtype = col(h, "obs/guide_type").astype(str)
            lq = col(h, "obs/low_quality")
            lq = lq.astype(str) == "True" if lq.dtype.kind in "OSU" else lq.astype(bool)
            pgid = col(h, "obs/perturbed_gene_id").astype(str)
            vname = col(h, "var/gene_name").astype(str)
            ip = h["X/indptr"][:]
        meta[d] = dict(gid=gid, gtype=gtype, lq=lq, pgid=pgid, vname=vname, n=len(gid), nnz=int(ip[-1]))
        log(d, "cells", len(gid), "guide types", dict(zip(*np.unique(gtype, return_counts=True))))
    # one var axis for all four files (asserted), mapped to the official genes by symbol
    # each file's own var axis, mapped to the official genes by symbol (the four files do not share one var order)
    axis, cms = {}, {}
    for d in DONORS:
        vn = meta[d]["vname"]
        first, dup = {}, 0
        for j, g in enumerate(vn):
            if g in first:
                dup += 1
                continue
            first[g] = j
        src_cols = np.array([first[g] for g in genes if g in first], dtype=np.int64)
        dst_cols = np.array([gpos[g] for g in genes if g in first], dtype=np.int32)
        cm = np.full(len(vn), -1, dtype=np.int32)
        cm[src_cols] = dst_cols
        cms[d] = cm
        axis[d] = dict(file_genes=int(len(vn)), dup_symbols=int(dup), overlap=int(len(src_cols)),
                       missing=int(len(genes) - len(src_cols)), present=set(g for g in genes if g in first))
    common = set.intersection(*(axis[d]["present"] for d in DONORS))
    union = set.union(*(axis[d]["present"] for d in DONORS))
    common_pos = np.array(sorted(gpos[g] for g in common), dtype=np.int32)
    for d in DONORS:     # a gene some donor lacks is zero for every cell (never a donor-composition artefact)
        cm = cms[d].copy()
        cm[(cm >= 0) & ~np.isin(cm, common_pos)] = -1
        np.save(out / f"colmap_{d}.npy", cm)
    log("axis", {d: {k: v for k, v in axis[d].items() if k != "present"} for d in DONORS},
        "common", len(common), "union", len(union))

    # ---- usable cells and labels
    ntc_types = sorted({t for d in DONORS for t in np.unique(meta[d]["gtype"]) if "non" in t.lower()})
    assert len(ntc_types) == 1, f"non-targeting guide_type is ambiguous: {ntc_types}"
    eid2sym = dict(zip(de["target_contrast"].astype(str), de["target_contrast_gene_name"].astype(str)))
    kd_all = set(kd[kd.columns[0]].astype(str))
    for d in DONORS:
        m = meta[d]
        tg = np.unique(m["gid"][(m["gtype"] != ntc_types[0]) & (m["gid"] != "multi-guide") & (m["gid"] != "")])
        in_sg, in_kd = np.mean([g in sg_tid for g in tg]), np.mean([g in kd_all for g in tg])
        log(d, "targeting guide ids", len(tg), "in sgRNA table", round(in_sg, 4), "in kd table", round(in_kd, 4), tg[:3])
        assert in_sg >= 0.9 and in_kd >= 0.9, "guide ids of obs do not match the suppl tables"
    lab = {}
    for d in DONORS:
        m = meta[d]
        usable = (~m["lq"]) & (m["gid"] != "multi-guide") & (m["gid"] != "") & (m["gid"] != "nan")
        tid = np.array([sg_tid.get(g, p) for g, p in zip(m["gid"], m["pgid"])], dtype=object)
        sym = np.array([eid2sym.get(t, "") for t in tid], dtype=object)
        okg = np.isin(m["gid"], list(kd_ok))
        L = np.full(m["n"], "", dtype=object)
        L[usable & (m["gtype"] == ntc_types[0])] = CTRL
        tsel = usable & (m["gtype"] != ntc_types[0]) & okg & (sym != "")
        L[tsel] = sym[tsel]
        lab[d] = L
    ctrl_n = np.array([(lab[d] == CTRL).sum() for d in DONORS], dtype=np.float64)
    pi = ctrl_n / ctrl_n.sum()
    q_ctrl, q_tgt = largest_remainder(pi, N_CTRL), largest_remainder(pi, N_CELLS)
    avail = {}
    for d in DONORS:
        u, c = np.unique(lab[d][lab[d] != ""], return_counts=True)
        for a, b in zip(u, c):
            avail.setdefault(a, np.zeros(len(DONORS), dtype=np.int64))[DONORS.index(d)] = b

    # ---- target sets
    pool = de[de["ontarget_significant"].astype(bool) & (de["n_cells_target"] >= N_CELLS)]
    pool = sorted(t for t in pool["target_contrast_gene_name"].astype(str) if t in gpos and t not in ess)
    order = np.random.default_rng(SEED).permutation(np.array(pool, dtype=object))
    P, skipped = [], []
    def qualifies(t) -> bool:      # Amendment 4: every donor holds its quota (no shortfall transfers)
        return t in avail and bool(np.all(avail[t] >= q_tgt))

    for t in order:
        if qualifies(t):
            P.append(str(t))
            if len(P) == N_P:
                break
        else:
            skipped.append(str(t))
    assert len(P) == N_P, "pool exhausted"
    pool_set = set(pool)
    S = sorted(t for t in abc if t in pool_set and qualifies(t))
    S_any400 = sum(1 for t in abc if t in pool_set and t in avail and avail[t].sum() >= N_CELLS)
    P = sorted(P)

    # ---- per-group rows (donor-stratified quotas, shortfall moved to donors with spare cells)
    rows = {d: [] for d in DONORS}   # (row, group_label)
    short = {}

    def take(label: str, total: int, quota: np.ndarray, k: int):
        av = avail[label].copy()
        n = np.minimum(quota, av)
        deficit = total - int(n.sum())
        if deficit:
            short[label] = int(deficit)
            assert label == CTRL, f'{label}: a qualified target cannot have a shortfall'
            spare = av - n
            add = largest_remainder(spare.astype(np.float64), deficit) if spare.sum() >= deficit else None
            assert add is not None and np.all(add <= spare), label
            n = n + add
        rng = np.random.default_rng([SEED, k])
        for i, d in enumerate(DONORS):
            if n[i]:
                cand = np.flatnonzero(lab[d] == label)
                rows[d].extend((int(r), label) for r in rng.choice(cand, size=int(n[i]), replace=False))
        return n

    take(CTRL, N_CTRL, q_ctrl, 0)
    groups = sorted(set(P) | set(S))
    for k, t in enumerate(groups, start=1):
        take(t, N_CELLS, q_tgt, k)

    for d in DONORS:
        rr = sorted(rows[d])
        np.savez(out / f"rows_{d}.npz", rows=np.array([r for r, _ in rr], dtype=np.int64),
                 labels=np.array([l for _, l in rr], dtype=str))
    info = dict(condition=cond, seed=SEED, donors=list(DONORS), ntc_guide_type=ntc_types[0],
                cells_per_donor={d: meta[d]["n"] for d in DONORS}, nnz_per_donor={d: meta[d]["nnz"] for d in DONORS},
                usable_controls_per_donor=dict(zip(DONORS, ctrl_n.astype(int).tolist())), pi=pi.tolist(),
                quota_controls=q_ctrl.tolist(), quota_target=q_tgt.tolist(), shortfalls=short,
                axis_per_donor={d: {k: v for k, v in axis[d].items() if k != "present"} for d in DONORS},
                axis_overlap_common=len(common), axis_overlap_union=len(union),
                axis_overlap=len(common), file_dup_symbols=int(sum(axis[d]["dup_symbols"] for d in DONORS)),
                axis_not_in_every_donor=sorted(union - common)[:200],
                axis_missing_genes_first50=[g for g in genes if g not in union][:50],
                pool_size=len(pool), P=P, P_sha=hashlib.sha256("\n".join(P).encode()).hexdigest(),
                P_skipped_for_cells=skipped, S=S, S_n=len(S), S_with_400_any_donor_mix=S_any400, rule='amendment4: every donor fills its quota', S_sha=hashlib.sha256("\n".join(S).encode()).hexdigest(),
                abc_in_pool=int(sum(t in pool_set for t in abc)),
                selected_rows_per_donor={d: len(rows[d]) for d in DONORS}, seconds_select=time.time() - t0)
    (out / "select.json").write_text(json.dumps(info, indent=1))
    log("select done", json.dumps({k: v for k, v in info.items() if k not in ("P", "S", "P_skipped_for_cells")}))


def extract(cond: str, src: Path, out: Path, d: str) -> None:
    """Read donor d's selected rows (sorted), map columns to the official axis, save a CSR npz."""
    t0 = time.time()
    z = np.load(out / f"rows_{d}.npz", allow_pickle=False)
    r_idx = z["rows"]
    colmap = np.load(out / f"colmap_{d}.npy")
    n_genes = len(pd.read_csv(ROOT / "data/challenge/gene_names.csv"))
    f = src / f"{d}_{cond}.assigned_guide.h5ad"
    with h5py.File(f, "r") as h:
        ip = h["X/indptr"][:]
        dd, di = h["X/data"], h["X/indices"]
        data, ind, ptr = [], [], [0]
        i = 0
        while i < len(r_idx):          # contiguous spans of sorted rows, read in windows of <= 64M entries
            a = int(ip[r_idx[i]])
            j = i
            while j + 1 < len(r_idx) and int(ip[r_idx[j + 1] + 1]) - a <= (64 << 20):
                j += 1
            b = int(ip[r_idx[j] + 1])
            bd, bi = dd[a:b], di[a:b]
            for r in r_idx[i:j + 1]:
                s, e = int(ip[r]) - a, int(ip[r + 1]) - a
                cm = colmap[bi[s:e]]
                ok = cm >= 0
                data.append(np.asarray(bd[s:e])[ok].astype(np.float32))
                ind.append(cm[ok])
                ptr.append(ptr[-1] + int(ok.sum()))
            i = j + 1
        M = sp.csr_matrix((np.concatenate(data), np.concatenate(ind).astype(np.int32), np.array(ptr, dtype=np.int64)),
                          shape=(len(r_idx), n_genes))
    M.sort_indices()
    assert np.all(M.data == np.round(M.data)) and (M.data.size == 0 or M.data.min() >= 0), "not integer counts"
    sp.save_npz(out / f"X_{d}.npz", M, compressed=False)
    log(d, "extracted", M.shape, M.nnz, f"{time.time() - t0:.0f}s")


def merge(cond: str, src: Path, out: Path) -> None:
    t0 = time.time()
    genes = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).to_numpy()
    info = json.loads((out / "select.json").read_text())
    P, S = info["P"], info["S"]
    blocks, labels = [], []
    for d in DONORS:
        blocks.append(sp.load_npz(out / f"X_{d}.npz"))
        labels.extend(np.load(out / f"rows_{d}.npz")["labels"].astype(str).tolist())
    X = sp.vstack(blocks).tocsr()
    labels = np.array(labels, dtype=object)
    assert np.all(X.data == np.round(X.data)) and X.data.min() >= 0, "counts are not non-negative integers"
    tot = np.asarray(X.sum(1)).ravel()
    med0 = float(np.median(tot))
    f_thin = min(1.0, DEPTH / med0)
    if f_thin < 1.0:
        rng = np.random.default_rng(SEED)
        X.data = rng.binomial(X.data.astype(np.int64), f_thin).astype(np.float32)
        X.eliminate_zeros()
    tot1 = np.asarray(X.sum(1)).ravel()

    def write(name: str, tset: list[str]):
        keep = np.isin(labels, tset + [CTRL])
        lab_k = labels[keep]
        o = np.argsort(np.where(lab_k == CTRL, "~", lab_k), kind="stable")   # targets sorted, controls last
        A = ad.AnnData(X=X[np.flatnonzero(keep)[o]].astype(np.float32),
                       obs=pd.DataFrame({"target_gene": pd.Categorical(lab_k[o].astype(str))},
                                        index=[f"c{i}" for i in range(int(keep.sum()))]),
                       var=pd.DataFrame(index=pd.Index(genes, name=None)))
        p = out / f"truth_{name}.h5ad"
        A.write_h5ad(p)
        (out / f"targets_{name}.txt").write_text("\n".join(tset) + "\n")
        h = hashlib.sha256(open(p, "rb").read()).hexdigest()
        log("wrote", p, A.shape, h)
        return dict(path=str(p), sha256=h, shape=list(A.shape),
                    cells_per_target=int(np.unique(lab_k[lab_k != CTRL], return_counts=True)[1].min()))

    info.update(depth_median_before=med0, thinning_f=f_thin, depth_median_after=float(np.median(tot1)),
                depth_median_after_ctrl=float(np.median(tot1[labels == CTRL])),
                depth_median_after_pert=float(np.median(tot1[labels != CTRL])),
                truth_P=write("P", P), truth_S=write("S", S), seconds_merge=time.time() - t0)
    (out / "select.json").write_text(json.dumps(info, indent=1))
    log("done", json.dumps({k: v for k, v in info.items() if k not in ("P", "S", "P_skipped_for_cells")}))


if __name__ == "__main__":
    cmd, cond, src, out = sys.argv[1], sys.argv[2], Path(sys.argv[3]), Path(sys.argv[4])
    out.mkdir(parents=True, exist_ok=True)
    if cmd == "select":
        select(cond, src, out)
    elif cmd == "extract":
        extract(cond, src, out, sys.argv[5])
    else:
        merge(cond, src, out)

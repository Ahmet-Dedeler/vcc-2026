"""466 — where does a perturbation programme's lab signature live? Honest (disjoint-control) agreement with single cells,
by target class, within a lab across lines and across labs.

Lines: Xaira's Orion HCT116 and HEK293T (one lab, genome-wide CRISPRi, parquet single cells) and Weissman's K562
genome-wide screen (scPerturb single cells). Targets: NON-ESSENTIAL (outside the Replogle K562/RPE1 essential sets) and
ESSENTIAL, each present in all three lines with >= MIN_CELLS cells, sampled with a fixed seed. Genes: the symbols all
three measure, then those >= 5 CPM in every line's controls.

Per line, each cell of a target or of the controls is assigned to half A or B at random (seeded). Programme = per-gene
P_down over the target set (sign of log1p(5e4 p/sum p) of the target's summed counts minus the controls'). Readings,
per target class:
  within-line reliability   P(half A vs controls A) against P(half B vs controls B)     (disjoint: honest)
  cross-line, same lab      HCT116 half A vs HEK293T half B (and the mirror)             (Xaira x Xaira)
  cross-lab                 each Xaira line's half vs K562's other half                  (Xaira x Weissman)
462h/i found on essential targets that cross-lab agreement is nearly as high as same-lab (biology), and 450h found
large lab effects on broad target sets. This measures both classes, with cells, in one design.

usage: uv run python experiments/466-lab-signature-cells/lab466.py [--smoke]   -> work/lab466.json
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse as sp

pa.set_memory_pool(pa.system_memory_pool())   # Arrow's default pool keeps freed batch memory: RSS grew 6 GB per batch

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
ORION = ROOT / "data/external/orion"
K562 = ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
TS, SEED, MIN_CELLS, N_NONESS, N_ESS = 5e4, 466, 50, 1500, 500
CPM_FLOORS = (5, 50)   # 50 CPM keeps ~15-cell halves out of the zero-count regime (the smoke test's artefact)
CONTROL = re.compile(r"^(non.?targeting|ntc|control|safe.?harbou?r|neg)", re.I)
t0 = time.time()
import resource
log = lambda *x: print("[466]", *x, f"({time.time() - t0:.0f}s, maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30:.1f} GB)", flush=True)


def lp(c):
    return np.log1p(TS * c / c.sum())


def pdn(M):
    nz = (M != 0).sum(0); dn = (M < 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)


def r(a, b, m):
    m = m & np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[m] - .5, b[m] - .5)[0, 1])


def orion_targets(line, files):
    """Pass 1 (cheap): cells per gene_target for a line (pass_guide_filter == 1)."""
    n = {}
    for f in files:
        t = pq.read_table(f, columns=["gene_target", "pass_guide_filter"]).to_pandas()
        t = t[t["pass_guide_filter"] == 1]
        for k, v in t["gene_target"].value_counts().items():
            n[k] = n.get(k, 0) + int(v)
    return n


def orion_sums(files, keys, gcols, n_tok):
    """Pass 2: summed counts per key (target or '__ctrl__') on the chosen gene columns, for halves A/B (full = A + B)."""
    kpos = {k: i for i, k in enumerate(keys)}
    S = np.zeros((len(keys), 2, len(gcols)))
    tok2col32 = np.full(n_tok, -1, dtype=np.int32); tok2col32[gcols] = np.arange(len(gcols), dtype=np.int32)
    rng = np.random.default_rng(SEED)
    import pyarrow.compute as pc
    for j, f in enumerate(files):
        # Arrow offsets/values directly (to_pandas() on list columns made one numpy object per cell); no filter copy;
        # int32/float32 arrays; the batch is freed before the next read (peak ~2x one batch otherwise)
        t = pq.read_table(f, columns=["gene_token_id", "gene_expression", "gene_target", "pass_guide_filter"])
        passing = t["pass_guide_filter"].to_numpy() == 1
        lab = np.asarray(t["gene_target"].to_pylist(), dtype=object).astype(str)
        lab = np.where([bool(CONTROL.match(x)) for x in lab], "__ctrl__", lab)
        sel = passing & np.array([x in kpos for x in lab])
        half_all = rng.integers(0, 2, len(lab))              # drawn for EVERY cell of the batch: stable per batch
        if not sel.any():
            del t; continue
        tok = t["gene_token_id"].combine_chunks(); exv = t["gene_expression"].combine_chunks()
        off = tok.offsets.to_numpy().astype(np.int64)
        col_all = tok2col32[tok.values.to_numpy()]; mask = col_all >= 0
        kept_before = np.concatenate([[0], np.cumsum(mask, dtype=np.int64)])
        indptr = kept_before[off]
        vals = exv.values.to_numpy()[mask].astype(np.float32)
        cols = col_all[mask]
        del t, tok, exv, col_all, mask, kept_before
        M = sp.csr_matrix((vals, cols, indptr), shape=(len(off) - 1, len(gcols)))
        del vals, cols
        cells = np.flatnonzero(sel)
        M = M[cells]
        lab = lab[sel]; half = half_all[sel]
        ki = np.array([kpos[x] for x in lab])
        for hv in (0, 1):
            m = half == hv
            G = sp.csr_matrix((np.ones(int(m.sum())), (ki[m], np.flatnonzero(m))), shape=(len(keys), len(cells)))
            S[:, hv] += np.asarray((G @ M).todense(), dtype=np.float64)
        if (j + 1) % 20 == 0:
            log(f"  {j + 1}/{len(files)} batches")
    return S


def k562_sums(keys, gcols, smoke):
    kpos = {k: i for i, k in enumerate(keys)}
    with h5py.File(K562, "r") as f:
        n = f["obs"]["perturbation"]
        lab = n["categories"][:].astype(str)[n["codes"][:]] if isinstance(n, h5py.Group) else n[:].astype(str)
        lab = np.where(lab == "control", "__ctrl__", lab)
        X = f["X"]; nrow = X.shape[0] if not smoke else 60000
        S = np.zeros((len(keys), 2, len(gcols)))
        rng = np.random.default_rng(SEED + 1); half = rng.integers(0, 2, X.shape[0])
        ki = np.array([kpos.get(x, -1) for x in lab])
        step = 3886 * 4
        for s in range(0, nrow, step):
            e = min(nrow, s + step)
            blk = np.asarray(X[s:e], dtype=np.float64)[:, gcols]
            for hv in (0, 1):
                m = (half[s:e] == hv) & (ki[s:e] >= 0)
                if m.any():
                    np.add.at(S[:, hv], ki[s:e][m], blk[m])
            if (s // step) % 25 == 0:
                log(f"  K562 rows {e}/{nrow}")
    return S


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--part", default="all", choices=("all", "HCT116", "HEK293T", "K562", "summary"),
                    help="run one line's sums in its own process (bounded memory; the allocator keeps ~18 GB) or the summary")
    a = ap.parse_args()
    gm = pd.read_parquet(ORION / "metadata/gene_metadata.parquet").sort_values("gene_token_id")
    n_tok = int(gm["gene_token_id"].max()) + 1
    tokname = np.array([""] * n_tok, dtype=object); tokname[gm["gene_token_id"].to_numpy()] = gm["gene_name"].to_numpy()
    with h5py.File(K562, "r") as f:
        g = f["var"]["gene_name"]; kg = (g["categories"][:].astype(str)[g["codes"][:]] if isinstance(g, h5py.Group) else g[:].astype(str))
        n = f["obs"]["perturbation"]
        klab = n["categories"][:].astype(str)[n["codes"][:]] if isinstance(n, h5py.Group) else n[:].astype(str)
    common = sorted(set(kg) & set(x for x in tokname if x))
    tokpos = {}
    for i, x in enumerate(tokname):
        if x and x not in tokpos:
            tokpos[x] = i
    kgpos = {x: i for i, x in enumerate(kg)}
    o_cols = np.array([tokpos[x] for x in common]); k_cols = np.array([kgpos[x] for x in common])
    log(f"common genes {len(common)}")
    files = {ln: sorted((ORION / "data").glob(f"{ln}_Batch*.parquet"), key=lambda p: int(re.search(r'Batch(\d+)', p.name)[1]))
             for ln in ("HCT116", "HEK293T")}
    if a.smoke:
        files = {ln: fs[:2] for ln, fs in files.items()}
    ncell = {ln: orion_targets(ln, fs) for ln, fs in files.items()}
    ku, kc = np.unique(klab, return_counts=True); ncell["K562"] = dict(zip(ku, kc))
    ess = set(l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt")) | \
          set(l.strip() for l in open(ROOT / "data/donors/rpe1/perts.txt"))
    ok = lambda t: all(ncell[ln].get(t, 0) >= (MIN_CELLS if not a.smoke else 1) for ln in ncell) and not CONTROL.match(t) and t in set(common) | set(ncell["K562"])
    cand = sorted(t for t in ncell["HCT116"] if ok(t))
    rng = np.random.default_rng(SEED)
    noness = [t for t in cand if t not in ess]; esst = [t for t in cand if t in ess]
    noness = sorted(rng.choice(noness, min(N_NONESS, len(noness)), replace=False).tolist())
    esst = sorted(rng.choice(esst, min(N_ESS, len(esst)), replace=False).tolist())
    keys = noness + esst + ["__ctrl__"]
    log(f"targets: non-essential {len(noness)}, essential {len(esst)} (candidates {len(cand)})")
    W = HERE / "work"; W.mkdir(exist_ok=True); tag = "_smoke" if a.smoke else ""
    (W / f"keys466{tag}.json").write_text(json.dumps({"noness": noness, "ess": esst, "genes": common}))
    def part(ln):
        S_ = orion_sums(files[ln], keys, o_cols, n_tok) if ln != "K562" else k562_sums(keys, k_cols, a.smoke)
        np.save(W / f"sums466_{ln}{tag}.npy", S_); log(f"{ln} sums saved")
    if a.part in ("HCT116", "HEK293T", "K562"):
        part(a.part); return
    if a.part == "all":
        for ln in ("HCT116", "HEK293T", "K562"):
            part(ln)
    S = {ln: np.load(W / f"sums466_{ln}{tag}.npy") for ln in ("HCT116", "HEK293T", "K562")}
    c = len(keys) - 1
    res = {"genes": len(common), "min_cells": MIN_CELLS, "classes": {}}
    for floor in CPM_FLOORS:
      expressed = np.ones(len(common), bool)
      for ln in S:
          C = S[ln][c, 0] + S[ln][c, 1]; expressed &= (1e6 * C / C.sum()) >= floor
      res[f"expressed_all_cpm{floor}"] = int(expressed.sum())
      for cls, T in (("non_essential", noness), ("essential", esst)):
        rows = [keys.index(t) for t in T]
        P = {}
        for ln in S:
            Sa, Sb = S[ln][:, 0], S[ln][:, 1]; Ca, Cb = Sa[c], Sb[c]
            P[ln] = {"A": pdn(np.stack([lp(Sa[i]) - lp(Ca) for i in rows])),
                     "B": pdn(np.stack([lp(Sb[i]) - lp(Cb) for i in rows]))}
        m = expressed
        within = {ln: r(P[ln]["A"], P[ln]["B"], m) for ln in P}
        x = lambda u, v: float(np.mean([r(P[u]["A"], P[v]["B"], m), r(P[u]["B"], P[v]["A"], m)]))
        rec = {"targets": len(T), "genes": int(m.sum()), "within_line_disjoint": within,
               "cross_line_same_lab (HCT116 x HEK293T)": x("HCT116", "HEK293T"),
               "cross_lab (HCT116 x K562)": x("HCT116", "K562"), "cross_lab (HEK293T x K562)": x("HEK293T", "K562")}
        wx = np.mean([within["HCT116"], within["HEK293T"]])
        rec["same_lab_over_within"] = rec["cross_line_same_lab (HCT116 x HEK293T)"] / wx
        rec["cross_lab_mean"] = float(np.mean([rec["cross_lab (HCT116 x K562)"], rec["cross_lab (HEK293T x K562)"]]))
        res["classes"][f"{cls}_cpm{floor}"] = rec
        log(f"{cls} cpm>={floor}: {json.dumps(rec)}")
    (HERE / "work").mkdir(exist_ok=True)
    (HERE / ("work/lab466_smoke.json" if a.smoke else "work/lab466.json")).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""377 — same-line L1000 response-subspace prior on cross-line CRISPRi direction.

Stages: base  (donor consensus predictions + truth per test line -> work/base377.npz)
         l1000 (row-chunked L1000 shRNA reads for HEPG2/MCF7/JURKAT -> work/l1000_*.npz)
         score (subspaces, corrections, endpoints, bootstrap -> work/scores377.json)

Reads only outside this directory. Writes only work/.
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
import sys, json, time, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
ROOT = os.path.dirname(os.path.dirname(HERE))
BANKS = {
    "k562_gwps": "data/donors/k562_gwps",
    "hepg2": "experiments/166-donor-acquisition/donors/nadig_hepg2",
    "hek293t": "data/donors/orion_hek293t",
    "hct116": "data/donors/orion_hct116",
    "rpe1": "data/donors/rpe1",
    "jurkat": "data/donors/nadig_jurkat",
}
TEST_LINES = ["hepg2", "jurkat"]
SEED = 3770909

def rms_direction(a):
    a = np.asarray(a, dtype=np.float64)
    rms = np.sqrt(np.mean(a * a, axis=-1))
    return a / np.maximum(rms[..., None], 1e-12), rms

def unit(a):
    a = np.asarray(a, dtype=np.float64)
    return a / np.maximum(np.linalg.norm(a, axis=-1, keepdims=True), 1e-12)

def load_bank(name):
    p = os.path.join(ROOT, BANKS[name])
    genes = [l.strip() for l in open(os.path.join(p, "genes.txt"))]
    perts = [l.strip() for l in open(os.path.join(p, "perts.txt"))]
    delta = np.load(os.path.join(p, "delta.npy"), mmap_mode="r")
    return genes, perts, delta

def stage_base():
    import pandas as pd
    gi = pd.read_csv(os.path.join(ROOT, "data/external/cmap/geneinfo_beta.txt"), sep="\t")
    lm = [s for s in gi[gi.feature_space == "landmark"].gene_symbol.tolist()]
    print(f"landmarks: {len(lm)}", flush=True)
    banks = {k: load_bank(k) for k in BANKS}
    out = {"landmarks": np.array(lm)}
    for L in TEST_LINES:
        Lg, Lp, Ld = banks[L]
        donors = [k for k in BANKS if k != L]
        # gene axis: landmarks present in L and every donor
        have = [set(banks[k][0]) for k in BANKS]
        axis = [g for g in lm if all(g in h for h in have)]
        print(f"{L}: axis={len(axis)}", flush=True)
        ax_idx = {k: np.array([banks[k][0].index(g) for g in axis]) for k in BANKS}
        # targets: L targets carried by >=2 donors
        Tp = {k: banks[k][1] for k in BANKS}
        Tpos = {k: {t: i for i, t in enumerate(Tp[k])} for k in BANKS}
        targets = [t for t in Lp if sum(t in Tpos[k] for k in donors) >= 2]
        print(f"{L}: targets={len(targets)}", flush=True)
        T, G = len(targets), len(axis)
        P = np.zeros((T, G))   # source_mean unit-RMS
        Y = np.zeros((T, G))   # truth unit-RMS
        cov = np.zeros(T, dtype=int)
        Lpos = Tpos[L]
        for i, t in enumerate(targets):
            Y[i] = rms_direction(Ld[Lpos[t]][ax_idx[L]][None, :])[0][0]
            vecs = []
            for k in donors:
                if t in Tpos[k]:
                    v = np.asarray(banks[k][2][Tpos[k][t]])[ax_idx[k]]
                    vecs.append(rms_direction(v[None, :])[0][0])
            vecs = np.array(vecs)
            cov[i] = len(vecs)
            m = vecs.mean(0)
            P[i] = m / max(np.sqrt(np.mean(m * m)), 1e-12)
        # 334-style transfer gain: pool (mean_others, y_d) over donor pseudo-recipients
        num = np.zeros(G); den = np.zeros(G)
        for d in donors:
            others = [k for k in donors if k != d]
            shared = [t for t in Tpos[d] if sum(t in Tpos[k] for k in others) >= 1]
            # cap pseudo-recipient targets for speed (deterministic first 3000)
            shared = shared[:3000]
            if not shared: continue
            M = np.zeros((len(shared), G)); Yd = np.zeros((len(shared), G))
            for i, t in enumerate(shared):
                vecs = []
                for k in others:
                    if t in Tpos[k]:
                        v = np.asarray(banks[k][2][Tpos[k][t]])[ax_idx[k]]
                        vecs.append(rms_direction(v[None, :])[0][0])
                M[i] = np.array(vecs).mean(0)
                v = np.asarray(banks[d][2][Tpos[d][t]])[ax_idx[d]]
                Yd[i] = rms_direction(v[None, :])[0][0]
            num += (M * Yd).sum(0); den += (M * M).sum(0)
        prior = float(np.median(den))
        gain = (num + prior) / np.maximum(den + prior, 1e-12)
        print(f"{L}: gain median={np.median(gain):.4f} iqr=[{np.quantile(gain,.25):.4f},{np.quantile(gain,.75):.4f}]", flush=True)
        out[f"{L}_axis"] = np.array(axis)
        out[f"{L}_targets"] = np.array(targets)
        out[f"{L}_Pmean"] = P
        out[f"{L}_Pdiag"] = P * gain[None, :]
        out[f"{L}_Y"] = Y
        out[f"{L}_cov"] = cov
    np.savez_compressed(os.path.join(WORK, "base377.npz"), **out)
    print("wrote base377.npz", flush=True)

def stage_l1000(which="HEPG2,MCF7,JURKAT"):
    import h5py
    cells_wanted = set(which.split(","))
    f = h5py.File(os.path.join(ROOT, "data/external/cmap/level5_beta_trt_sh.gctx"), "r")
    M = f["0/DATA/0/matrix"]
    N, Gtot = M.shape
    print(f"gctx: {N} x {Gtot}", flush=True)
    import pandas as pd
    gi = pd.read_csv(os.path.join(ROOT, "data/external/cmap/geneinfo_beta.txt"), sep="\t")
    row_ids = [x.decode() for x in f["0/META/ROW/id"][()]]
    lm_ids = set(str(r) for r in gi[gi.feature_space == "landmark"].gene_id)
    lm_pos = np.array([i for i, r in enumerate(row_ids) if r in lm_ids])
    print(f"landmark cols: {len(lm_pos)}", flush=True)
    # column cell labels (one full decode, cheap vs matrix)
    col_ids = [x.decode() for x in f["0/META/COL/id"][()]]
    col_cell = np.array([c.split("_")[1] for c in col_ids])
    col_gene = np.array([c.split(":")[1] if ":" in c else "" for c in col_ids])
    for c in sorted(cells_wanted):
        print(f"  {c}: {(col_cell == c).sum()}", flush=True)
    CH = 2000
    acc = {c: [] for c in cells_wanted}
    acc_gene = {c: [] for c in cells_wanted}
    t0 = time.time()
    for s in range(0, N, CH):
        e = min(s + CH, N)
        mask = np.array([c in cells_wanted for c in col_cell[s:e]])
        if not mask.any():
            continue
        block = M[s:e][mask][:, lm_pos].astype(np.float32)
        cc = col_cell[s:e][mask]
        gg = col_gene[s:e][mask]
        for c in cells_wanted:
            m = cc == c
            if m.any():
                acc[c].append(block[m])
                acc_gene[c].append(gg[m])
        if (s // CH) % 20 == 0:
            print(f"  rows {s}/{N} {time.time()-t0:.0f}s", flush=True)
    for c in cells_wanted:
        X = np.concatenate(acc[c], 0) if acc[c] else np.zeros((0, len(lm_pos)), np.float32)
        gg = np.concatenate(acc_gene[c], 0) if acc_gene[c] else np.array([])
        # landmark symbols in lm_pos order
        id2s = dict(zip(gi.gene_id.astype(str), gi.gene_symbol))
        syms = np.array([id2s[row_ids[i]] for i in lm_pos])
        np.savez_compressed(os.path.join(WORK, f"l1000_{c}.npz"), X=X, genes=syms, shrna_target=gg)
        print(f"wrote l1000_{c}.npz {X.shape}", flush=True)
    f.close()

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["base", "l1000", "score"])
    ap.add_argument("--cells", default="HEPG2,MCF7,JURKAT")
    a = ap.parse_args()
    t0 = time.time()
    if a.stage == "base":
        stage_base()
    elif a.stage == "l1000":
        stage_l1000(a.cells)
    print(f"done {time.time()-t0:.1f}s", flush=True)

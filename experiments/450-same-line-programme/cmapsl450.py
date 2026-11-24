"""450l: does LINCS L1000 (shRNA consensus signatures, CGS) carry SAME-LINE information? Truth = Nadig 2025 HepG2 CRISPRi
Perturb-seq (essential targets). Predictors: CMap CGS signatures of HEPG2 (same line, other technology) vs the other core
lines (A549, MCF7, HT29, PC3, A375, VCAP), on the targets every line covers, landmark genes only (176: inferred genes carry
no CRISPRi-shared signal). Per-target centred cosine (each source centred over the evaluated targets), and sign agreement
on each truth row's top-50 |delta| landmark genes. If HEPG2 >> the others, L1000 lines can join the same-line registry."""
from pathlib import Path

import h5py
import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
gi = pl.read_csv(ROOT / "data/external/cmap/geneinfo_beta.txt", separator="\t", infer_schema_length=0)
id2sym = dict(zip(gi["gene_id"].to_list(), gi["gene_symbol"].to_list())); id2fs = dict(zip(gi["gene_id"].to_list(), gi["feature_space"].to_list()))
nd = ROOT / "experiments/166-donor-acquisition/donors/nadig_hepg2"
TD = np.load(nd / "delta.npy", mmap_mode="r"); TP = [l.strip() for l in open(nd / "perts.txt")]; TG = [l.strip() for l in open(nd / "genes.txt")]
LINES = ["HEPG2", "A549", "MCF7", "HT29", "PC3", "A375", "VCAP"]
with h5py.File(ROOT / "data/external/cmap/level5_beta_trt_sh.gctx", "r") as h:
    cid = np.asarray([x.decode() for x in h["0/META/COL/id"][:]]); rid = np.asarray([x.decode() for x in h["0/META/ROW/id"][:]])
    is_cgs = np.char.startswith(cid, "CGS")
    parts = [c.split(":") for c in cid[is_cgs]]; line = np.asarray([p[0].split("_")[1] for p in parts]); tgt = np.asarray([p[1] for p in parts])
    sym = np.asarray([id2sym.get(r, "") for r in rid]); fs = np.asarray([id2fs.get(r, "") for r in rid])
    tgi = {g: i for i, g in enumerate(TG)}
    lm = np.array([(f == "landmark") and (s in tgi) for s, f in zip(sym, fs)])
    genes = sym[lm]; tcols = np.array([tgi[s] for s in genes])
    idx_all = np.flatnonzero(is_cgs)
    common = set(TP)
    per = {}
    for L in LINES:
        common &= set(tgt[line == L])
    common = sorted(common)
    print(f"landmark genes on Nadig's axis {lm.sum()}; targets in Nadig HepG2 and every core line: {len(common)}")
    for L in LINES:
        sel = np.flatnonzero(line == L); tpos = {t: k for k, t in enumerate(tgt[sel])}
        rows = idx_all[sel[[tpos[t] for t in common]]]; order = np.argsort(rows)
        M = np.empty((len(common), rid.size), dtype=np.float32); M[order] = h["0/DATA/0/matrix"][rows[order], :]
        per[L] = M[:, lm].astype(np.float64)
tpi = {p: i for i, p in enumerate(TP)}
T = np.asarray(TD[[tpi[t] for t in common]], dtype=np.float64)[:, tcols]
gpos = {g: j for j, g in enumerate(genes)}
for r, t in enumerate(common):          # the target's own gene out, on both sides
    if t in gpos:
        T[r, gpos[t]] = 0.0
        for L in LINES:
            per[L][r, gpos[t]] = 0.0
u = lambda A: A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-12)
Tc = T - T.mean(0); top = np.argsort(-np.abs(T), axis=1)[:, :50]
res = {}
for L in LINES:
    P = per[L]; Pc = P - P.mean(0)
    cc = np.sum(u(Pc) * u(Tc), 1); sg = np.mean(np.sign(np.take_along_axis(P, top, 1)) == np.sign(np.take_along_axis(T, top, 1)), 1)
    res[L] = (cc, sg)
    print(f"  CMap {L:6s}{' (same line)' if L == 'HEPG2' else '            '}  centred cos {cc.mean():+.4f}   sign@top50 {sg.mean():.4f}")
others = np.mean([res[L][0] for L in LINES if L != "HEPG2"], axis=0)
d = res["HEPG2"][0] - others
rng = np.random.default_rng(0)
bs = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(2000)]
print(f"HEPG2 - mean(other lines) centred cos per target: {d.mean():+.4f} [95% {np.percentile(bs, 2.5):+.4f}, {np.percentile(bs, 97.5):+.4f}]; "
      f"HEPG2 best of 7 on {int(np.mean([res['HEPG2'][0] >= max(res[L][0][i] for L in LINES) for i in range(len(common))]) * len(common))} / {len(common)} targets")

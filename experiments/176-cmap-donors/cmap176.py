"""176 -- CMap shRNA consensus signatures (trt_sh.cgs, 9 core lines) as sign predictors on H1.
Gates, in order: (1) platform agreement with our production donors on shared (target, landmark)
pairs; (2) sign accuracy vs the h1-test reference DE pool, landmark genes then inferred genes;
(3) added to the production trio."""
import h5py, numpy as np, polars as pl, sys, time, json
sys.path.insert(0, "experiments/153-de-bench"); sys.path.insert(0, "experiments/172-sign-confidence-key")
from debench import Surface
from donors172 import DonorHandle
t0 = time.time(); OUT = "experiments/176-cmap-donors/work"; import os; os.makedirs(OUT, exist_ok=True)
gi = pl.read_csv("data/external/cmap/geneinfo_beta.txt", separator="\t", infer_schema_length=0)
id2sym = dict(zip(gi["gene_id"].to_list(), gi["gene_symbol"].to_list())); id2fs = dict(zip(gi["gene_id"].to_list(), gi["feature_space"].to_list()))
S = Surface("experiments/164-callset-artefact/de164.a15.npz")
genes = np.asarray([str(g) for g in S.genes]); gpos = {g: i for i, g in enumerate(genes)}; targets = [str(p) for p in S.perts]
pool = S.real_sig & S.adj; sT = np.sign(S.real_lfc)
with h5py.File("data/external/cmap/level5_beta_trt_sh.gctx", "r") as h:
    cid = np.asarray([x.decode() for x in h["0/META/COL/id"][:]]); rid = np.asarray([x.decode() for x in h["0/META/ROW/id"][:]])
    is_cgs = np.char.startswith(cid, "CGS")
    print(f"{is_cgs.sum()} CGS signatures of {cid.size} [{time.time()-t0:.0f}s]")
    # parse "CGS001_A375_96H:A2M:1" -> line A375, target A2M
    parts = [c.split(":") for c in cid[is_cgs]]; line = np.asarray([p[0].split("_")[1] for p in parts]); tgt = np.asarray([p[1] for p in parts])
    # gene column mapping onto the DE-universe axis
    sym = np.asarray([id2sym.get(r, "") for r in rid]); fs = np.asarray([id2fs.get(r, "") for r in rid])
    on = np.asarray([s in gpos for s in sym]); col_axis = np.asarray([gpos[s] if s in gpos else -1 for s in sym])
    lm = (fs == "landmark") & on; inf = (fs != "landmark") & on
    print(f"L1000 genes on the DE-universe axis: landmark {lm.sum()}, inferred {inf.sum()}")
    idx_all = np.flatnonzero(is_cgs)
    lines = sorted(set(line)); Z = {}
    for L in lines:
        sel = np.flatnonzero(line == L); rows = idx_all[sel]; order = np.argsort(rows)
        M = np.empty((rows.size, rid.size), dtype=np.float32)
        for k in range(0, rows.size, 2000):
            rr = rows[order[k:k+2000]]; M[order[k:k+2000]] = h["0/DATA/0/matrix"][rr, :]
        tg = tgt[sel]
        # keep h1-test targets only, in S.perts order; NaN where absent
        Zl = np.full((85, len(genes)), np.nan, dtype=np.float32); have = np.zeros(85, bool)
        tpos = {t: i for i, t in enumerate(tg)}
        for i, t in enumerate(targets):
            if t in tpos: have[i] = True; Zl[i, col_axis[on]] = M[tpos[t]][on]
        Z[L] = (Zl, have); print(f"  {L:<8} {rows.size:5d} CGS sigs, {have.sum()}/85 h1 targets  [{time.time()-t0:.0f}s]", flush=True)
        np.savez_compressed(f"{OUT}/cgs_{L}.npz", Z=Zl, have=have, genes=genes, perts=np.asarray(targets))
np.save(f"{OUT}/landmark_mask.npy", np.isin(np.arange(len(genes)), col_axis[lm])); LM = np.isin(np.arange(len(genes)), col_axis[lm]); INF = np.isin(np.arange(len(genes)), col_axis[inf])
# production donors
D = {}
for n in ("k562_gwps", "orion_hct116", "orion_hek293t"):
    hd = DonorHandle(n); Sd, c = hd.rows(targets, genes); D[n] = (Sd.astype(np.float64), c)
def agree(a, b, m): m = m & (a != 0) & (b != 0) & np.isfinite(a) & np.isfinite(b); return (np.sign(a[m]) == np.sign(b[m])).mean(), int(m.sum())
print("\n=== GATE 1: platform agreement -- sign agreement of each CMap line with each production donor, shared (target, LANDMARK gene) pairs ===")
print(f"{'line':<8}" + "".join(f"{n:>16}" for n in D) + f"{'K562 x HCT116 (ref)':>22}")
kh = agree(D['k562_gwps'][0], D['orion_hct116'][0], np.ones_like(pool) & LM[None, :])
for L, (Zl, have) in Z.items():
    row = f"{L:<8}"
    for n, (Sd, c) in D.items():
        a, k = agree(Zl, Sd, have[:, None] & c[:, None] & LM[None, :]); row += f"{a:>10.4f}({k:>4d})"
    print(row + f"{kh[0]:>22.4f}")
print("\n=== GATE 2: sign accuracy vs the h1-test reference DE pool ===")
def acc(key, m): m = m & np.isfinite(key) & (key != 0); return (np.sign(key[m]) == sT[m]).mean(), int(m.sum())
print(f"{'source':<22}{'landmark acc':>14}{'n':>7}{'inferred acc':>14}{'n':>8}")
for n, (Sd, c) in D.items():
    a1, k1 = acc(Sd, pool & c[:, None] & LM[None, :]); a2, k2 = acc(Sd, pool & c[:, None] & INF[None, :]); print(f"{n:<22}{a1:>14.4f}{k1:>7d}{a2:>14.4f}{k2:>8d}")
cons3 = sum(Sd for Sd, _ in D.values()) / 3
a1, k1 = acc(cons3, pool & LM[None, :]); a2, k2 = acc(cons3, pool & INF[None, :]); print(f"{'trio equal mean':<22}{a1:>14.4f}{k1:>7d}{a2:>14.4f}{k2:>8d}")
for L, (Zl, have) in Z.items():
    # CMap z-scores: knockdown signature. Sign convention: z<0 means gene DOWN after knockdown -- same as our delta.
    a1, k1 = acc(Zl.astype(np.float64), pool & have[:, None] & LM[None, :]); a2, k2 = acc(Zl.astype(np.float64), pool & have[:, None] & INF[None, :])
    print(f"{'CMap ' + L:<22}{a1:>14.4f}{k1:>7d}{a2:>14.4f}{k2:>8d}")
# all CMap lines pooled (mean z where available)
Zs = np.stack([Zl for Zl, _ in Z.values()]); cm = np.nanmean(Zs, axis=0); ncov = np.isfinite(Zs).sum(0)
a1, k1 = acc(cm, pool & LM[None, :]); a2, k2 = acc(cm, pool & INF[None, :]); print(f"{'CMap 9-line mean':<22}{a1:>14.4f}{k1:>7d}{a2:>14.4f}{k2:>8d}")
print("\n=== GATE 3: added to the production trio (standardised per source, on pairs where CMap has coverage) ===")
def std(A): v = np.nanstd(A[np.isfinite(A) & (A != 0)]); return A / (v if v > 0 else 1)
cov = np.isfinite(cm)
for w in (0.0, 0.25, 0.5, 1.0, 2.0):
    comb = std(cons3) + w * np.nan_to_num(std(cm))
    a1, k1 = acc(comb, pool & cov & LM[None, :]); a2, k2 = acc(comb, pool & cov & INF[None, :]); a3, k3 = acc(comb, pool & cov)
    print(f"  trio + {w:.2f} x CMap-mean:   landmark {a1:.4f} (n={k1})   inferred {a2:.4f} (n={k2})   all covered {a3:.4f} (n={k3})")
print(f"\n[{time.time()-t0:.0f}s] done")

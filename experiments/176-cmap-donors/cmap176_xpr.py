"""176b -- CMap CRISPR-KO (trt_xpr) signatures, collapsed to one mean z per (line, target), on H1.
Same three gates as cmap176.py plus the own-target check.  Separates knockdown DEPTH (shRNA was
weak) from PLATFORM (bulk L1000, 96 h, inferred axis)."""
import h5py, numpy as np, polars as pl, sys, time, os
sys.path.insert(0, "experiments/153-de-bench"); sys.path.insert(0, "experiments/172-sign-confidence-key")
from debench import Surface
from donors172 import DonorHandle
t0 = time.time(); OUT = "experiments/176-cmap-donors/work"
gi = pl.read_csv("data/external/cmap/geneinfo_beta.txt", separator="\t", infer_schema_length=0)
id2sym = dict(zip(gi["gene_id"].to_list(), gi["gene_symbol"].to_list())); id2fs = dict(zip(gi["gene_id"].to_list(), gi["feature_space"].to_list()))
S = Surface("experiments/164-callset-artefact/de164.a15.npz")
genes = np.asarray([str(g) for g in S.genes]); gpos = {g: i for i, g in enumerate(genes)}; targets = [str(p) for p in S.perts]; tset = set(targets)
pool = S.real_sig & S.adj; sT = np.sign(S.real_lfc)
sig = pl.read_parquet("data/external/cmap/genetic_sigs.parquet").filter(pl.col("pert_type") == "trt_xpr")
sig = sig.filter(pl.col("cmap_name").is_in(list(tset)))
print(f"xpr sigs hitting h1-test targets: {sig.height} across {sig['cell_iname'].n_unique()} lines [{time.time()-t0:.0f}s]")
want = dict(zip(sig["sig_id"].to_list(), zip(sig["cell_iname"].to_list(), sig["cmap_name"].to_list())))
with h5py.File("data/external/cmap/level5_beta_trt_xpr.gctx", "r") as h:
    cid = np.asarray([x.decode() for x in h["0/META/COL/id"][:]]); rid = np.asarray([x.decode() for x in h["0/META/ROW/id"][:]])
    sym = np.asarray([id2sym.get(r, "") for r in rid]); fs = np.asarray([id2fs.get(r, "") for r in rid])
    on = np.asarray([s in gpos for s in sym]); col_axis = np.asarray([gpos[s] if s in gpos else -1 for s in sym])
    LM = np.zeros(len(genes), bool); LM[col_axis[(fs == "landmark") & on]] = True
    INF = np.zeros(len(genes), bool); INF[col_axis[(fs != "landmark") & on]] = True
    rows = np.asarray([i for i, c in enumerate(cid) if c in want]); print(f"matching columns in gctx: {rows.size}")
    M = np.empty((rows.size, rid.size), dtype=np.float32)
    for k in range(0, rows.size, 1000): M[k:k+1000] = h["0/DATA/0/matrix"][rows[k:k+1000], :]
    print(f"read {rows.size} x {rid.size} [{time.time()-t0:.0f}s]")
lines_of = np.asarray([want[cid[r]][0] for r in rows]); tg_of = np.asarray([want[cid[r]][1] for r in rows])
Z = {}
for L in sorted(set(lines_of)):
    Zl = np.full((85, len(genes)), np.nan, dtype=np.float32); have = np.zeros(85, bool); nsig = []
    for i, t in enumerate(targets):
        m = (lines_of == L) & (tg_of == t)
        if m.any(): have[i] = True; Zl[i, col_axis[on]] = M[m].mean(0)[on]; nsig.append(int(m.sum()))
    Z[L] = (Zl, have); np.savez_compressed(f"{OUT}/xpr_{L}.npz", Z=Zl, have=have, genes=genes, perts=np.asarray(targets))
    print(f"  {L:<8} {have.sum()}/85 h1 targets, median {np.median(nsig) if nsig else 0:.0f} sigs/target")
D = {n: DonorHandle(n).rows(targets, genes) for n in ("k562_gwps", "orion_hct116", "orion_hek293t")}
D = {n: (Sd.astype(np.float64), c) for n, (Sd, c) in D.items()}
def agree(a, b, m): m = m & (a != 0) & (b != 0) & np.isfinite(a) & np.isfinite(b); return ((np.sign(a[m]) == np.sign(b[m])).mean() if m.any() else np.nan), int(m.sum())
def acc(key, m): m = m & np.isfinite(key) & (key != 0); return ((np.sign(key[m]) == sT[m]).mean() if m.any() else np.nan), int(m.sum())
K = D["k562_gwps"][0]
print("\n=== own-target z (sign convention + knockout depth) ===")
for L, (Zl, have) in Z.items():
    own = np.array([Zl[i, gpos[t]] for i, t in enumerate(targets) if have[i] and t in gpos and np.isfinite(Zl[i, gpos[t]])])
    if own.size >= 5: print(f"  {L:<8} median own-z {np.median(own):+.2f}, frac<0 {np.mean(own<0):.2f}, n={own.size}   | |z| median {np.nanmedian(np.abs(Zl[have])):.2f}")
print("\n=== GATE 1: sign agreement with K562 on the H1 DE POOL (K562 x HCT116 reference: landmark 0.540) ===")
for L, (Zl, have) in Z.items():
    a1, k1 = agree(Zl.astype(float), K, pool & have[:, None] & LM[None, :]); a2, k2 = agree(Zl.astype(float), K, pool & have[:, None] & INF[None, :])
    if k1 > 300: print(f"  {L:<8} landmark {a1:.4f} (n={k1})   inferred {a2:.4f} (n={k2})")
print("\n=== GATE 2: sign accuracy vs the H1 reference DE pool (K562 reference: landmark 0.562 / inferred 0.568) ===")
for L, (Zl, have) in Z.items():
    a1, k1 = acc(Zl.astype(float), pool & have[:, None] & LM[None, :]); a2, k2 = acc(Zl.astype(float), pool & have[:, None] & INF[None, :])
    if k1 > 300: print(f"  CMap-KO {L:<8} landmark {a1:.4f} (n={k1})   inferred {a2:.4f} (n={k2})")
Zs = np.stack([Zl for Zl, _ in Z.values()]); cm = np.nanmean(Zs, axis=0)
a1, k1 = acc(cm, pool & LM[None, :]); a2, k2 = acc(cm, pool & INF[None, :]); print(f"  {'KO all-line mean':<16} landmark {a1:.4f} (n={k1})   inferred {a2:.4f} (n={k2})")
cons3 = sum(Sd for Sd, _ in D.values()) / 3
def std(A): v = np.nanstd(A[np.isfinite(A) & (A != 0)]); return A / (v if v > 0 else 1)
cov = np.isfinite(cm)
print("\n=== GATE 3: trio + w x KO-mean, on KO-covered pairs ===")
for w in (0.0, 0.5, 1.0, 2.0):
    comb = std(cons3) + w * np.nan_to_num(std(cm))
    a1, k1 = acc(comb, pool & cov & LM[None, :]); a2, k2 = acc(comb, pool & cov & INF[None, :]); a3, k3 = acc(comb, pool & cov)
    print(f"  w={w:.1f}:  landmark {a1:.4f} (n={k1})   inferred {a2:.4f} (n={k2})   all covered {a3:.4f} (n={k3})")
print(f"\n[{time.time()-t0:.0f}s] done")

"""516: does the H1 gain of the Stouffer-6 sign over production survive at the cancer folds' power? Same predictions
(435 production signs; Stouffer over K562-GWPS, Orion HCT116/HEK293T, CD4 x3), only the truth's power changes (480's
thinned H1 gates: cells per target and depth). Also the production analog used on the folds (donor mean + 0.34 CD4)
against production and against the Stouffer, to read the folds' numbers on H1's own scale."""
import numpy as np
import importlib.util as u
s = u.spec_from_file_location("v", "vote516.py"); V = u.module_from_spec(s); s.loader.exec_module(V)
from common516 import ROOT, sign_table

d = V.d; G = d["genes"]; T = d["targets"]; prod = d["prod"]; calls = d["calls"]
gi = {g: i for i, g in enumerate(G)}
blk = np.zeros(len(G), bool); blk[np.argsort(-d["ctrl"])[:2200]] = True
Z = V.stouffer(["K", "H", "E", "R", "S8", "S48"], False, 3.0)
x = {k: V.z[f"{V.SRC[k]}__x"].astype(float) if not V.SRC[k].startswith("cd4_") else V.z[f"{V.SRC[k]}__log_fc"].astype(float)
     for k in V.SRC}
cd4m = np.nanmean(np.stack([x["R"], x["S8"], x["S48"]]), 0)
an = np.sign(np.nan_to_num(np.nanmean(np.stack([x["K"], x["H"], x["E"]]), 0)) + 0.34 * np.nan_to_num(cd4m))
cands = {
    "stz6 |Z|>=0.5 uncalled vs prod": (prod, np.where(~calls & (np.abs(np.nan_to_num(Z)) >= 0.5), np.sign(np.nan_to_num(Z)), 0)),
    "stz6 |Z|>=0 rest vs prod": (prod, np.where(~blk[None, :], np.sign(np.nan_to_num(Z)), 0)),
    "analog(mean+.34cd4) rest vs prod": (prod, np.where(~blk[None, :], an, 0)),
    "stz6 |Z|>=0 rest vs analog": (np.where(an != 0, an, prod), np.where(~blk[None, :], np.sign(np.nan_to_num(Z)), 0)),
    "stz6 |Z|>=0.5 rest vs analog": (np.where(an != 0, an, prod),
                                     np.where(~blk[None, :] & (np.abs(np.nan_to_num(Z)) >= 0.5), np.sign(np.nan_to_num(Z)), 0)),
}
for lvl in ["h1_full_ka", "h1_c250", "h1_c400_d56", "h1_c150", "h1_c250_d56", "h1_thin", "h1_c150_d56", "h1_c100"]:
    g = np.load(ROOT / f"experiments/480-nmae-instrument/work/gate-{lvl}.npz", allow_pickle=True)
    gg = g["genes"].astype(str); assert (g["targets"].astype(str) == T).all()
    q = np.full((len(T), len(G)), 1.0); l = np.full((len(T), len(G)), np.nan)
    idx = np.array([gi.get(x_, -1) for x_ in gg]); ok = idx >= 0
    q[:, idx[ok]] = g["q"][:, ok]; l[:, idx[ok]] = g["lfc"][:, ok]
    own = G[None, :] == T[:, None]
    rs = (q < 0.05) & ~own; pool = rs & np.isfinite(l) & (l != 0)
    nr = rs.sum(1); band = np.select([nr < 20, nr < 100], [0, 1], 2); truth = np.sign(np.nan_to_num(l))
    print(f"\n## {lvl}: pool cells {pool.sum()}, targets per band {[int((band==b).sum()) for b in range(3)]}")
    for nm, (base, C) in cands.items():
        H = np.where(C != 0, C, base)
        for sn, sub in V.SUBS.items():
            tab = sign_table(truth, pool & sub[:, None], band, base, H)
            a = tab["all"]
            print(f"  {nm:34s} {sn:6s} {a[0]:.4f}->{a[1]:.4f} {a[2]:+.4f} [{a[3]:+.4f},{a[4]:+.4f}]  "
                  f"b<20 {tab['<20'][2]:+.3f} b20 {tab['20-100'][2]:+.3f} b100 {tab['>=100'][2]:+.4f}")

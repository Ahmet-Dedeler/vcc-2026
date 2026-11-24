"""516 exploration b: majority vote of six same-target sources (K562-GWPS, Orion HCT116, Orion HEK293T, CD4 x3) by
agreement level, vs production on the same pool cells, split block (top-2,200 H1 ctrl CPM) / rest."""
import numpy as np
from common516 import PANEL_DONORS, cd4_rows, donor_rows, known_answers, load_h1, price

d = load_h1(); got, ok, _ = known_answers(d); assert ok; print("KA PASS", got)
T, G = d["targets"], d["genes"]; pool, band, truth, prod = d["pool"], d["band"], d["truth"], d["prod"]
blk = np.zeros(len(G), bool); blk[np.argsort(-d["ctrl"])[:2200]] = True
rows = {dn: donor_rows(dn, T, G)[0] for dn in PANEL_DONORS}
cd4 = cd4_rows(T, G, layers=("log_fc", "zscore"))
for c, r in cd4.items():
    rows["cd4_" + c] = r["log_fc"]
S = np.stack([np.sign(np.nan_to_num(R)) for R in rows.values()])       # sources x T x G
npos = (S > 0).sum(0); nneg = (S < 0).sum(0); ncov = npos + nneg
maj = np.sign(npos - nneg)
M0 = pool & (prod != 0)
print(f"{'region':6s} {'ncov':>4s} {'agree':>5s} {'cells':>7s} {'prod':>6s} {'maj':>6s} {'diff':>7s}")
for rn, rm in (("block", blk), ("rest", ~blk)):
    for need in (2, 3, 4, 5, 6):
        m = M0 & rm[None, :] & (np.maximum(npos, nneg) >= need) & (np.minimum(npos, nneg) == 0)
        if m.sum() < 50: continue
        print(f"{rn:6s}  all>= {need}: {m.sum():7d} {(prod[m]==truth[m]).mean():.4f} {(maj[m]==truth[m]).mean():.4f} "
              f"{(maj[m]==truth[m]).mean()-(prod[m]==truth[m]).mean():+.4f}   prod agrees w/ maj {(prod[m]==maj[m]).mean():.3f}")
    for lo, hi in ((1, 1), (2, 2), (3, 3), (4, 4), (5, 6)):
        mg = np.abs(npos - nneg)
        m = M0 & rm[None, :] & (mg >= lo) & (mg <= hi)
        if m.sum() < 50: continue
        print(f"{rn:6s}  margin {lo}-{hi}: {m.sum():7d} {(prod[m]==truth[m]).mean():.4f} {(maj[m]==truth[m]).mean():.4f} "
              f"{(maj[m]==truth[m]).mean()-(prod[m]==truth[m]).mean():+.4f}")
# candidates priced end to end: majority vote with margin >= m on rest genes only / everywhere
for mm in (1, 2, 3):
    for region, rm in (("rest", ~blk), ("all", np.ones(len(G), bool))):
        C = np.where((np.abs(npos - nneg) >= mm) & rm[None, :], maj, 0)
        price(d, C, f"majority6 margin>={mm} on {region}")

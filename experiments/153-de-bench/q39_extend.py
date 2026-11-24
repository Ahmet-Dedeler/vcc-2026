"""Q39 — extend the shared-direction sweep, and price the profile error a19 died on.

q37/q38 stopped at beta=1.0 and I extrapolated the net linearly from beta=0.35, which was wrong:
`reach` is a THRESHOLD statistic (k* jumps once head purity crosses 0.90), so it grows
superlinearly while the `pds` cost is roughly linear.  beta=0.60 already nets +0.0067.

Also prices ||P - truth||^2 / ||truth||^2 in the reference's own log1p-pseudobulk space
(real_delta_hct116.npz), scale-matched -- the term `mse` reads, and the member that killed a19
because no screen measured it.
"""
import numpy as np, csv
from debench import Surface, overall_delta
from exact_pds import pds

hct = Surface("../095-clean-harness/de_thinned.npz")
hek = Surface("de_hek293t_ref.npz")
gh = {g: i for i, g in enumerate(hek.genes)}
col = np.array([gh.get(g, -1) for g in hct.genes]); ok = col >= 0
sh = np.zeros(hct.genes.size)
sh[ok] = np.nanmean(np.where(np.isfinite(hek.real_lfc), hek.real_lfc, np.nan), 0)[col[ok]]
u = np.nan_to_num(sh); u /= np.linalg.norm(u)

canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
panel = {r["target_gene"] for r in csv.DictReader(open("../../data/challenge/pert_counts.csv"))}
R = np.load("real_delta_hct116.npz", allow_pickle=True)
rg = {str(g): i for i, g in enumerate(R["genes"])}; rt = {str(t): i for i, t in enumerate(R["perts"])}
pg = np.array([i for i, g in enumerate(hct.genes) if g in rg and g not in panel])
pcol = np.array([rg[hct.genes[i]] for i in pg])
prow = np.array([rt.get(t, -1) for t in hct.perts]); keep = prow >= 0
Dreal = R["delta"][np.ix_(prow[keep], pcol)].astype(np.float64)
P0 = hct.pred_lfc.copy(); sc = float(np.median(np.linalg.norm(P0, axis=1)))

def prof_err(P):
    """scale-matched relative profile error, the core of mse's numerator"""
    A = P[np.ix_(keep, pg)].astype(np.float64)
    a = float((A * Dreal).sum() / max((A * A).sum(), 1e-12))     # optimal global scale
    return float(np.linalg.norm(a * A - Dreal) ** 2 / np.linalg.norm(Dreal) ** 2)

def row(beta):
    P = P0 + beta * sc * u[None, :]
    hct.pred_lfc = P
    hn = hd = 0
    for p in range(hct.perts.size):
        i = np.flatnonzero(hct.real_sig[p] & hct.adj[p])
        if i.size < 10: continue
        m = i[np.abs(P[p, i]) >= np.percentile(np.abs(P[p, i]), 90)]
        hn += int((np.sign(P[p, m]) == np.sign(hct.real_lfc[p, m])).sum()); hd += m.size
    return (hn/max(hd,1), hct.reach(key=[np.abs(P[p]) for p in range(hct.perts.size)]),
            hct.nmae(), pds(P[np.ix_(keep, pg)].astype(np.float64), Dreal), prof_err(P))

b = row(0.0)
print(f"{'beta':>6} {'head pur':>9} {'reach':>9} {'nmae':>8} {'pds':>8} {'prof err':>9} "
      f"{'Δreach':>8} {'Δnmae':>8} {'Δpds':>8} {'NET ovr':>9}")
print(f"{0.0:6.2f} {b[0]:9.4f} {b[1]:9.5f} {b[2]:8.5f} {b[3]:8.5f} {b[4]:9.4f} "
      f"{0:8.4f} {0:8.4f} {0:8.4f} {0:+9.5f}")
best = (None, -9)
for beta in (0.6, 1.0, 1.5, 2.0, 3.0, 5.0):
    s = row(beta)
    dr = overall_delta(d_reach=s[1]-b[1]); dn = -(s[2]-b[2])/0.60/6; dp = (s[3]-b[3])/0.452/6
    net = dr + dn + dp
    if net > best[1]: best = (beta, net)
    print(f"{beta:6.2f} {s[0]:9.4f} {s[1]:9.5f} {s[2]:8.5f} {s[3]:8.5f} {s[4]:9.4f} "
          f"{dr:+8.4f} {dn:+8.4f} {dp:+8.4f} {net:+9.5f}")
print(f"\nbest beta {best[0]} at {best[1]:+.5f} overall (reach + nmae + pds only)")
print("profile err rises with beta -> `mse` is the unpriced risk, exactly a19's failure mode.")

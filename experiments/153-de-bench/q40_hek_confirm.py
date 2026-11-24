"""Q40 — the shared-direction result on the second surface, donors swapped.

HCT116 was evaluated with HEK293T's mean as the shared vector.  Here HEK293T is evaluated with
HCT116's mean -- so neither run uses the evaluated line's own mean, and the two are independent.
"""
import numpy as np, csv
from debench import Surface, overall_delta
from exact_pds import pds

hek = Surface("de_hek293t_ref.npz")
hct = Surface("../095-clean-harness/de_thinned.npz")
gh = {g: i for i, g in enumerate(hct.genes)}
col = np.array([gh.get(g, -1) for g in hek.genes]); ok = col >= 0
sh = np.zeros(hek.genes.size)
sh[ok] = np.nanmean(np.where(np.isfinite(hct.real_lfc), hct.real_lfc, np.nan), 0)[col[ok]]
u = np.nan_to_num(sh); u /= np.linalg.norm(u)

panel = {r["target_gene"] for r in csv.DictReader(open("../../data/challenge/pert_counts.csv"))}
R = np.load("real_delta_hek293t.npz", allow_pickle=True)
rg = {str(g): i for i, g in enumerate(R["genes"])}; rt = {str(t): i for i, t in enumerate(R["perts"])}
pg = np.array([i for i, g in enumerate(hek.genes) if g in rg and g not in panel])
pcol = np.array([rg[hek.genes[i]] for i in pg])
prow = np.array([rt.get(t, -1) for t in hek.perts]); keep = prow >= 0
Dreal = R["delta"][np.ix_(prow[keep], pcol)].astype(np.float64)

# HEK293T has no shipped prediction in its table; use the 104 HEK baseline as P0
z = np.load("../104-context-heldout-credit-adapter/work/hek293t/frozen.npz", allow_pickle=True)
g104 = {str(g): i for i, g in enumerate(z["genes"])}; t104 = {str(t): i for i, t in enumerate(z["targets"])}
P0 = np.zeros_like(hek.real_lfc)
gc = np.array([g104.get(g, -1) for g in hek.genes]); pr = np.array([t104.get(t, -1) for t in hek.perts])
m1, m2 = pr >= 0, gc >= 0
P0[np.ix_(m1, m2)] = z["baseline"][np.ix_(pr[m1], gc[m2])]
sc = float(np.median(np.linalg.norm(P0, axis=1)))

def row(beta):
    P = P0 + beta * sc * u[None, :]
    hek.pred_lfc = P
    hn = hd = 0
    for p in range(hek.perts.size):
        i = np.flatnonzero(hek.real_sig[p] & hek.adj[p])
        if i.size < 10: continue
        m = i[np.abs(P[p, i]) >= np.percentile(np.abs(P[p, i]), 90)]
        hn += int((np.sign(P[p, m]) == np.sign(hek.real_lfc[p, m])).sum()); hd += m.size
    return (hn/max(hd,1), hek.reach(key=[np.abs(P[p]) for p in range(hek.perts.size)]),
            hek.nmae(), pds(P[np.ix_(keep, pg)].astype(np.float64), Dreal))

b = row(0.0)
print(f"HEK293T surface, shared vector = HCT116's mean (neither line uses its own)")
print(f"{'beta':>6} {'head pur':>9} {'reach':>9} {'nmae':>8} {'pds':>8} {'NET ovr':>9}")
print(f"{0.0:6.2f} {b[0]:9.4f} {b[1]:9.5f} {b[2]:8.5f} {b[3]:8.5f} {0:+9.5f}")
for beta in (0.35, 0.6, 1.0, 1.5, 2.0):
    s = row(beta)
    net = overall_delta(d_reach=s[1]-b[1]) - (s[2]-b[2])/0.60/6 + (s[3]-b[3])/0.452/6
    print(f"{beta:6.2f} {s[0]:9.4f} {s[1]:9.5f} {s[2]:8.5f} {s[3]:8.5f} {net:+9.5f}")

"""Q38 — the shared programme's pds cost comes from adding the SAME vector to all 300.

If instead each target carries its OWN loading on that programme -- beta_p * shared, with beta_p
varying by perturbation -- the predictions no longer shift identically, so the component becomes a
DIMENSION OF VARIATION rather than a constant, and `pds` should be far less damaged while the
per-gene sign benefit is retained.

That is biologically the right shape: the shared programme is growth arrest / stress, and
knockdowns differ enormously in how strongly they induce it.

Oracle first -- use the TRUE per-target loading (each target's own projection onto the shared
direction, from the donor's measured effects).  If the oracle does not recover the `pds`, the
idea is dead and no amount of work on predicting beta_p matters.  If it does, beta_p becomes a
supervised, per-target regression problem -- one scalar per target, not 18,533.
"""
import numpy as np, csv
from debench import Surface, overall_delta
from exact_pds import pds

hct = Surface("../095-clean-harness/de_thinned.npz")
hek = Surface("de_hek293t_ref.npz")
gh = {g: i for i, g in enumerate(hek.genes)}; th = {t: i for i, t in enumerate(hek.perts)}
col = np.array([gh.get(g, -1) for g in hct.genes]); ok = col >= 0
sh = np.zeros(hct.genes.size)
sh[ok] = np.nanmean(np.where(np.isfinite(hek.real_lfc), hek.real_lfc, np.nan), 0)[col[ok]]
sh = np.nan_to_num(sh); u = sh / np.linalg.norm(sh)

# true per-target loading: project the DONOR's own effect for that target onto the shared axis
load = np.ones(hct.perts.size)
for p, t in enumerate(hct.perts):
    if t in th:
        v = np.nan_to_num(hek.real_lfc[th[t]])[col.clip(0)] * ok
        load[p] = float(v @ u)
load = load / np.median(np.abs(load[load != 0]))          # scale-free, median 1

canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
panel = {r["target_gene"] for r in csv.DictReader(open("../../data/challenge/pert_counts.csv"))}
R = np.load("real_delta_hct116.npz", allow_pickle=True)
rg = {str(g): i for i, g in enumerate(R["genes"])}; rt = {str(t): i for i, t in enumerate(R["perts"])}
pg = np.array([i for i, g in enumerate(hct.genes) if g in rg and g not in panel])
pcol = np.array([rg[hct.genes[i]] for i in pg])
prow = np.array([rt.get(t, -1) for t in hct.perts]); keep = prow >= 0
Dreal = R["delta"][np.ix_(prow[keep], pcol)].astype(np.float64)
P0 = hct.pred_lfc.copy(); sc = float(np.median(np.linalg.norm(P0, axis=1)))

def score(P):
    hct.pred_lfc = P
    num = den = hn = hd = 0
    for p in range(hct.perts.size):
        i = np.flatnonzero(hct.real_sig[p] & hct.adj[p])
        if i.size < 10: continue
        num += int((np.sign(P[p, i]) == np.sign(hct.real_lfc[p, i])).sum()); den += i.size
        m = i[np.abs(P[p, i]) >= np.percentile(np.abs(P[p, i]), 90)]
        hn += int((np.sign(P[p, m]) == np.sign(hct.real_lfc[p, m])).sum()); hd += m.size
    return (num/den, hn/max(hd,1), hct.reach(key=[np.abs(P[p]) for p in range(hct.perts.size)]),
            hct.nmae(), pds(P[np.ix_(keep, pg)].astype(np.float64), Dreal))

b0 = score(P0.copy())
print(f"{'arm':34} {'sign':>7} {'head':>7} {'reach':>8} {'nmae':>8} {'pds':>8} {'net ovr':>9}")
def line(lbl, s):
    net = (overall_delta(d_reach=s[2]-b0[2]) - (s[3]-b0[3])/0.60/6 + (s[4]-b0[4])/0.452/6)
    print(f"{lbl:34} {s[0]:7.4f} {s[1]:7.4f} {s[2]:8.5f} {s[3]:8.5f} {s[4]:8.5f} {net:+9.5f}")
line("production", b0)
for beta in (0.35, 0.60):
    line(f"constant  beta={beta}", score(P0 + beta*sc*u[None, :]))
    line(f"ORACLE per-target  beta={beta}", score(P0 + beta*sc*load[:, None]*u[None, :]))

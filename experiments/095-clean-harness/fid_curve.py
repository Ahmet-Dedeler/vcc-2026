"""`fid` as a function of call-set size, on the repaired harness.

`fid = precision x min(1, n_pred/n_real)`. exp092 measured this curve on the
depth-confounded harness and it pointed at ~850; the board then returned the wrong sign
twice (v14, v15). This re-measures it, and answers the open question about arm a3:
--keep-spread 0.0 can grow n_pred by up to 31%, and nobody has measured which way that
moves fid.

The call set at size m is our top-m genes by predicted adjusted p-value, which is the
ordering the emission actually produces.
"""
import numpy as np

d = np.load("experiments/095-clean-harness/de_thinned.npz", allow_pickle=True)
genes, perts = d["genes"].astype(str), d["perts"].astype(str)
rq, rl, pq, pl = d["real_padj"], d["real_lfc"], d["pred_padj"], d["pred_lfc"]
G, ALPHA = genes.size, 0.05
gpos = {g: i for i, g in enumerate(genes)}
B_FID, R_FID = 0.51239, 0.80964
B_JAC, R_JAC = 0.030376, 0.397270

MS = [250, 500, 1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000]
fid = {m: [] for m in MS}
jac = {m: [] for m in MS}
prod_fid, prod_jac, npred = [], [], []
for i, t in enumerate(perts):
    ok = np.ones(G, bool)
    j = gpos.get(t)
    if j is not None: ok[j] = False
    rs = ok & (rq[i] < ALPHA)
    ind = ok & (rl[i] != 0) & ~np.isnan(rl[i])
    n_real = int(rs.sum())
    agree = (np.sign(pl[i]) == np.sign(rl[i])) & (pl[i] != 0) & ~np.isnan(pl[i])
    order = np.argsort(np.where(ok, pq[i], np.inf))     # our own confidence order
    ps = ok & (pq[i] < ALPHA)
    sc = ps & ind
    npred.append(int(sc.sum()))
    if max(int(sc.sum()), n_real) > 0:
        prod_fid.append(int((sc & agree).sum()) / max(int(sc.sum()), n_real))
    u = int((rs | ps).sum()); prod_jac.append(1.0 if u == 0 else int((rs & ps).sum()) / u)
    for m in MS:
        sel = np.zeros(G, bool); sel[order[:m]] = True
        s2 = sel & ind
        if max(int(s2.sum()), n_real) > 0:
            fid[m].append(int((s2 & agree).sum()) / max(int(s2.sum()), n_real))
        u2 = int((rs | (sel & ok)).sum())
        jac[m].append(1.0 if u2 == 0 else int((rs & sel & ok).sum()) / u2)

pf, pj = np.nanmean(prod_fid), np.mean(prod_jac)
print(f"production: n_pred median {np.median(npred):.0f}   raw fid {pf:.4f}   raw jac {pj:.4f}\n")
print(f"{'call-set size m':<18}{'raw fid':>10}{'raw jac':>10}{'Δfid+jac (overall)':>22}")
for m in MS:
    f, jj = np.nanmean(fid[m]), np.mean(jac[m])
    dd = ((f - pf) / (R_FID - B_FID) + (jj - pj) / (R_JAC - B_JAC)) / 6
    star = "  <-- production" if abs(m - np.median(npred)) < 600 else ""
    print(f"{m:<18}{f:>10.4f}{jj:>10.4f}{dd:>+22.5f}{star}")

import os, numpy as np
S = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(S, "de_tables.npz"), allow_pickle=True)
genes = d["genes"].astype(str); perts = d["perts"].astype(str)
rq, rl = d["real_padj"], d["real_lfc"]
pp, pq, pl = d["pred_p"], d["pred_padj"], d["pred_lfc"]
A = 0.05; G = len(genes)
gpos = {g: i for i, g in enumerate(genes)}
fr = np.empty(G, dtype=np.int64); fr[np.argsort(genes)] = np.arange(G)

pre = []
for i, t in enumerate(perts):
    ok = np.ones(G, bool); j = gpos.get(t)
    if j is not None: ok[j] = False
    rs = ok & (rq[i] < A); ps = ok & (pq[i] < A)
    ind = ok & (rl[i] != 0) & ~np.isnan(rl[i])
    idx = np.flatnonzero(ps)
    key = np.lexsort((fr[idx], -np.abs(pl[i][idx]), pp[i][idx], pq[i][idx]))
    r = idx[key]
    adj = ind[r]
    hit = adj & (pl[i][r] != 0) & (np.sign(pl[i][r]) == np.sign(rl[i][r]))
    pre.append((int(rs.sum()), np.cumsum(adj), np.cumsum(hit), np.cumsum(rs[r]), r.size))

# official anchors, verified exactly against 040 RESULTS.md and 031 results.json
F0, FS = 0.5135, 0.3
J0, JS = 0.029, 0.37

def eval_ratio(rat):
    F, J, Ac, M = [], [], [], []
    for n_real, ca, ch, ci, npa in pre:
        m = npa if rat is None else int(min(round(rat * n_real), npa))
        if m <= 0:
            F.append(0.0); J.append(1.0 if n_real == 0 else 0.0); Ac.append(np.nan); M.append(0); continue
        npd, k, inter = int(ca[m-1]), int(ch[m-1]), int(ci[m-1])
        F.append(k / max(npd, n_real) if max(npd, n_real) else np.nan)
        u = n_real + m - inter
        J.append(1.0 if u == 0 else inter / u)
        Ac.append(k / npd if npd else np.nan); M.append(m)
    return np.nanmean(F), np.nanmean(J), np.nanmean(Ac), np.mean(M)

grid = [0.25,0.5,0.6,0.7,0.8,0.9,0.95,1.0,1.05,1.1,1.2,1.3,1.4,1.5,1.75,2.0,2.5,3.0,4.0,5.0,
        6.0,7.0,8.0,10.0,12.0,15.0,20.0,25.0,30.0,50.0,None]
print(f"{'m/n_real':>9} {'mean m':>8} {'fid':>8} {'jac':>8} {'A(m)':>8} | {'off_fid':>8} {'off_jac':>8} {'d avg6':>8}")
base = None; rowsout=[]
for rat in grid:
    f, j, a, mm = eval_ratio(rat)
    of, oj = (f - F0) / FS, (j - J0) / JS
    if rat is None: base = (of, oj)
    rowsout.append((rat, mm, f, j, a, of, oj))
for rat, mm, f, j, a, of, oj in rowsout:
    d6 = ((of - base[0]) + (oj - base[1])) / 6.0
    lab = "all" if rat is None else f"{rat:g}"
    print(f"{lab:>9} {mm:8.1f} {f:8.4f} {j:8.4f} {a:8.4f} | {of:8.4f} {oj:8.4f} {d6:+8.4f}")

fs = np.array([r[2] for r in rowsout[:-1]]); js = np.array([r[3] for r in rowsout[:-1]])
d6s = np.array([((r[5]-base[0])+(r[6]-base[1]))/6 for r in rowsout[:-1]])
labs = [r[0] for r in rowsout[:-1]]
print(f"\nargmax fid: m/n_real={labs[int(np.argmax(fs))]:g}  fid={fs.max():.4f}")
print(f"argmax jac: m/n_real={labs[int(np.argmax(js))]:g}  jac={js.max():.4f}")
print(f"argmax delta-avg6: m/n_real={labs[int(np.argmax(d6s))]:g}  d_avg6={d6s.max():+.4f}")

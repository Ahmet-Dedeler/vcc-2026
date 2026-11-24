import os, numpy as np, csv
S = "/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/76a42756-c5fb-4398-aea0-8513243571e9/scratchpad"
d = np.load(os.path.join(S, "de_tables.npz"), allow_pickle=True)
genes = d["genes"].astype(str); perts = d["perts"].astype(str)
rq, rl = d["real_padj"], d["real_lfc"]
pp, pq, pl = d["pred_p"], d["pred_padj"], d["pred_lfc"]
A = 0.05; G = len(genes); P = len(perts)
gpos = {g: i for i, g in enumerate(genes)}
gene_order = np.argsort(genes)          # feature-asc tiebreak rank
feat_rank = np.empty(G, dtype=np.int64); feat_rank[gene_order] = np.arange(G)

RATIOS = [0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0]
per = []          # per-pert base stats
curves = {r: [] for r in RATIOS}; curves["all"] = []
best_fid_m, best_jac_m = [], []

for i, t in enumerate(perts):
    m_ok = np.ones(G, dtype=bool); j = gpos.get(t)
    if j is not None: m_ok[j] = False
    rs = m_ok & (rq[i] < A)
    ps = m_ok & (pq[i] < A)
    ind = m_ok & (rl[i] != 0) & ~np.isnan(rl[i])
    n_real = int(rs.sum())
    idx = np.flatnonzero(ps)
    # ranking key: p_adj asc, p_value asc, |lfc| desc, feature asc
    key = np.lexsort((feat_rank[idx], -np.abs(pl[i][idx]), pp[i][idx], pq[i][idx]))
    ranked = idx[key]
    ok_adj = ind[ranked]
    hit = ok_adj & (pl[i][ranked] != 0) & ~np.isnan(pl[i][ranked]) & (np.sign(pl[i][ranked]) == np.sign(rl[i][ranked]))
    inR = rs[ranked]
    cum_adj = np.cumsum(ok_adj); cum_hit = np.cumsum(hit); cum_in = np.cumsum(inR)
    n_pred_all = ranked.size
    def at(m):
        m = int(min(max(m, 0), n_pred_all))
        if m == 0: return 0, 0, 0
        return int(cum_adj[m-1]), int(cum_hit[m-1]), int(cum_in[m-1])
    per.append((t, n_real, n_pred_all, *at(n_pred_all)[:2]))
    for r in RATIOS + ["all"]:
        m = n_pred_all if r == "all" else int(round(r * n_real))
        npd, k, inter = at(m)
        mm = min(m, n_pred_all) if r != "all" else n_pred_all
        fid = k / max(npd, n_real) if max(npd, n_real) > 0 else np.nan
        union = n_real + mm - inter
        jac = 1.0 if union == 0 else inter / union
        acc = k / npd if npd > 0 else np.nan
        curves[r].append((fid, jac, acc, mm, npd, k))
    # per-pert argmax over every possible m
    fid_all = np.zeros(n_pred_all + 1); jac_all = np.zeros(n_pred_all + 1)
    mm = np.arange(1, n_pred_all + 1)
    fid_all[1:] = cum_hit / np.maximum(cum_adj, n_real)
    jac_all[1:] = cum_in / (n_real + mm - cum_in)
    jac_all[0] = 1.0 if n_real == 0 else 0.0
    best_fid_m.append((np.argmax(fid_all), fid_all.max(), n_real))
    best_jac_m.append((np.argmax(jac_all), jac_all.max(), n_real))

per = np.array(per, dtype=object)
n_real = per[:,1].astype(float); n_pred = per[:,2].astype(float); k = per[:,4].astype(float)
c = n_pred / n_real
def q(x): return np.nanpercentile(x, [0,25,50,75,100])

print("="*78)
print("1) COUNTS PER PERTURBATION (n=132)")
print("="*78)
for name, x in [("n_real", n_real), ("n_pred", n_pred), ("c = n_pred/n_real", c)]:
    a = q(x); print(f"{name:20s} min={a[0]:10.3f} q1={a[1]:10.3f} med={a[2]:10.3f} q3={a[3]:10.3f} max={a[4]:10.3f} mean={np.nanmean(x):10.3f}")
print(f"{'k':20s} min={q(k)[0]:10.0f} q1={q(k)[1]:10.0f} med={q(k)[2]:10.0f} q3={q(k)[3]:10.0f} max={q(k)[4]:10.0f} mean={k.mean():10.3f}")
prec = k / np.maximum(n_pred, 1)
a = q(prec); print(f"{'k/n_pred':20s} min={a[0]:10.4f} q1={a[1]:10.4f} med={a[2]:10.4f} q3={a[3]:10.4f} max={a[4]:10.4f} mean={prec.mean():10.4f}")
print(f"\ntotals: sum n_real={n_real.sum():.0f}  sum n_pred={n_pred.sum():.0f}  sum k={k.sum():.0f}  pooled c={n_pred.sum()/n_real.sum():.2f}")
print(f"perts with n_pred > n_real: {(n_pred>n_real).sum()}/132   n_pred < n_real: {(n_pred<n_real).sum()}/132")
print(f"perts with n_real == 0: {(n_real==0).sum()}")

print("\n" + "="*78)
print("4) TRUNCATION COUNTERFACTUAL  (top-m by our own p_adj; m = ratio * n_real)")
print("="*78)
print(f"{'m/n_real':>10} {'mean m':>9} {'mean fid':>10} {'mean jac':>10} {'mean A(m)':>10} {'d fid':>9} {'d jac':>9}")
base_f = np.nanmean([x[0] for x in curves['all']]); base_j = np.nanmean([x[1] for x in curves['all']])
for r in RATIOS + ["all"]:
    arr = curves[r]
    f = np.nanmean([x[0] for x in arr]); jj = np.nanmean([x[1] for x in arr])
    ac = np.nanmean([x[2] for x in arr]); mn = np.mean([x[3] for x in arr])
    lab = "all" if r == "all" else f"{r:g}"
    print(f"{lab:>10} {mn:9.1f} {f:10.4f} {jj:10.4f} {ac:10.4f} {f-base_f:+9.4f} {jj-base_j:+9.4f}")

bf = np.array([b[1] for b in best_fid_m]); bfm = np.array([b[0] for b in best_fid_m], float)
bj = np.array([b[1] for b in best_jac_m]); bjm = np.array([b[0] for b in best_jac_m], float)
rr = np.where(n_real > 0, 1.0, np.nan)
print(f"\nper-pert ORACLE truncation (best m for each pert individually):")
print(f"  fid: mean={bf.mean():.4f}  median m/n_real={np.nanmedian(bfm/np.maximum(n_real,1)):.3f}")
print(f"  jac: mean={bj.mean():.4f}  median m/n_real={np.nanmedian(bjm/np.maximum(n_real,1)):.3f}")

with open(os.path.join(S, "counterfactual.csv"), "w", newline="") as fo:
    w = csv.writer(fo); w.writerow(["ratio","mean_m","mean_fid","mean_jac","mean_A"])
    for r in RATIOS + ["all"]:
        arr = curves[r]
        w.writerow([r, np.mean([x[3] for x in arr]), np.nanmean([x[0] for x in arr]),
                    np.nanmean([x[1] for x in arr]), np.nanmean([x[2] for x in arr])])
np.save(os.path.join(S, "per_pert_counts.npy"), np.stack([n_real, n_pred, k, c]))

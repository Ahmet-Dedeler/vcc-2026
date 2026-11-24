import os, numpy as np
S = "/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/76a42756-c5fb-4398-aea0-8513243571e9/scratchpad"
d = np.load(os.path.join(S, "de_tables.npz"), allow_pickle=True)
genes = d["genes"].astype(str); perts = d["perts"].astype(str)
rp, rq, rl = d["real_p"], d["real_padj"], d["real_lfc"]
pp, pq, pl = d["pred_p"], d["pred_padj"], d["pred_lfc"]
A = 0.05
gpos = {g: i for i, g in enumerate(genes)}
P, G = len(perts), len(genes)
print(f"perts={P} gated genes={G}")
print("targets resolved in gated universe:", sum(1 for t in perts if t in gpos))

rows = []
for i, t in enumerate(perts):
    m = np.ones(G, dtype=bool)
    j = gpos.get(t)
    if j is not None: m[j] = False            # on-target exclusion (both sides)
    rs = m & (rq[i] < A)                      # real significant set R
    ps = m & (pq[i] < A)                      # pred significant set P
    ind = m & (rl[i] != 0) & ~np.isnan(rl[i]) # adjudicable (in_denom)
    scored = ps & ind
    match = scored & (pl[i] != 0) & ~np.isnan(pl[i]) & (np.sign(pl[i]) == np.sign(rl[i]))
    n_real, n_pred, k = int(rs.sum()), int(scored.sum()), int(match.sum())
    inter = int((rs & ps).sum()); union = int((rs | ps).sum())
    fid = k / max(n_pred, n_real) if max(n_pred, n_real) > 0 else float("nan")
    jac = 1.0 if union == 0 else inter / union
    rows.append((t, n_real, n_pred, k, inter, union, fid, jac,
                 int(ps.sum()), int((ps & ~ind).sum())))

import csv
with open(os.path.join(S, "per_pert.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["target","n_real","n_pred","k","inter","union","fid","jac","n_pred_sig_all","n_pred_sig_unadjudicable"])
    w.writerows(rows)

fid = np.array([r[6] for r in rows]); jac = np.array([r[7] for r in rows])
print(f"MEAN fid = {np.nanmean(fid):.10f}   (scorer: 0.5312814007682398)")
print(f"MEAN jac = {np.nanmean(jac):.10f}   (scorer: 0.0340283619023056)")
print(f"rel err fid = {abs(np.nanmean(fid)-0.5312814007682398)/0.5312814007682398*100:.3f}%")
print(f"rel err jac = {abs(np.nanmean(jac)-0.0340283619023056)/0.0340283619023056*100:.3f}%")
print("unadjudicable pred-sig rows total:", sum(r[9] for r in rows))

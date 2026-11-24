"""094 step 3 — is the anti-correlation real, and what causes it?

r2 found Spearman(our p_adj, sign correct) = +0.50: the scorer walks our ranking
from the most significant down, and that is precisely the WRONG end. Before
building anything on this, rule out the boring explanations.
"""
import os
import numpy as np
from scipy import stats

E092 = "~/Code/virtual-cell-challenge-2026/experiments/092-call-set-truncation"
d = np.load(os.path.join(E092, "de_tables.npz"), allow_pickle=True)
g = np.load(os.path.join(E092, "gate.npz"), allow_pickle=True)
genes = d["genes"].astype(str); perts = d["perts"].astype(str)
rq, rl = d["real_padj"], d["real_lfc"]
pq, pl = d["pred_padj"], d["pred_lfc"]
G, P = genes.size, perts.size
ALPHA = 0.05
gpos = {gg: i for i, gg in enumerate(genes)}

# control mean CPM, restricted to the 10,258 genes in the DE table
cpm_all = g["ctrl_mean_cpm"]; keep = g["keep"]
cpm = cpm_all[keep] if keep.sum() == G else None
print(f"gate keeps {keep.sum()} of {keep.size} genes; DE table has {G} -> cpm aligned: {cpm is not None}")

rows = []
for i, t in enumerate(perts):
    ok = np.ones(G, bool)
    j = gpos.get(t)
    if j is not None: ok[j] = False
    rs = ok & (rq[i] < ALPHA)
    gs = np.flatnonzero(rs)
    if gs.size == 0: continue
    rlv, plv = rl[i][gs], pl[i][gs]
    ind = (~np.isnan(rlv)) & (rlv != 0) & (~np.isnan(plv)) & (plv != 0)
    gs = gs[ind]
    if gs.size < 20: continue
    rlv, plv = rl[i][gs], pl[i][gs]
    rows.append(np.column_stack([
        np.full(gs.size, i), gs, (np.sign(plv) == np.sign(rlv)).astype(float),
        pq[i][gs], np.abs(plv), np.sign(plv), np.sign(rlv), rlv, plv,
        (cpm[gs] if cpm is not None else np.zeros(gs.size)),
    ]))
M = np.vstack(rows)
pert_i, gi, mt, qv, absl, sp, sr, rlv, plv, cpmv = [M[:, k] for k in range(10)]
print(f"pooled pairs: {M.shape[0]}\n")

print("=== 1. Is it a global sign flip? ===")
print(f"  overall sign accuracy               {mt.mean():.4f}")
print(f"  our predicted lfc: frac positive    {(sp>0).mean():.4f}")
print(f"  reference    lfc: frac positive     {(sr>0).mean():.4f}")
print("  -> a global flip would put accuracy far below 0.5 everywhere. It does not.\n")

print("=== 2. Split by the SIGN we predict (top decile of our confidence) ===")
top = qv <= np.quantile(qv, 0.10)
for name, m in (("our top-10% confident", top), ("all pairs", np.ones_like(top, bool))):
    up = m & (sp > 0); dn = m & (sp < 0)
    print(f"  {name:<24} n={m.sum():>6}  we predict UP: n={up.sum():>6} acc={mt[up].mean():.4f}"
          f"   we predict DOWN: n={dn.sum():>6} acc={mt[dn].mean():.4f}")
print("  -> if one predicted direction is reliably wrong, the confident set is one-sided.\n")

print("=== 3. Does control expression level explain it? ===")
if cpm is not None:
    lc = np.log10(cpmv + 1e-3)
    print(f"  Spearman(log ctrl CPM, our p_adj)   = {stats.spearmanr(lc, qv).statistic:+.4f}")
    print(f"  Spearman(log ctrl CPM, sign correct)= {stats.spearmanr(lc, mt).statistic:+.4f}")
    qs = np.quantile(lc, np.linspace(0, 1, 6))
    print(f"  {'ctrl CPM quintile':<20}{'n':>8}{'mean p_adj rank':>18}{'sign acc':>11}{'frac we call sig':>18}")
    qr = np.argsort(np.argsort(qv)) / qv.size
    for k in range(5):
        m = (lc >= qs[k]) & (lc <= qs[k+1])
        print(f"  Q{k+1} [{qs[k]:+.2f},{qs[k+1]:+.2f}]{m.sum():>8}{qr[m].mean():>18.4f}"
              f"{mt[m].mean():>11.4f}{(qv[m]<ALPHA).mean():>18.4f}")
    # partial: does p_adj still predict correctness after conditioning on CPM?
    keepm = np.isfinite(lc)
    resid_q = qr[keepm] - np.poly1d(np.polyfit(lc[keepm], qr[keepm], 3))(lc[keepm])
    resid_m = mt[keepm] - np.poly1d(np.polyfit(lc[keepm], mt[keepm], 3))(lc[keepm])
    print(f"\n  PARTIAL Spearman(p_adj rank, correct | ctrl CPM) = "
          f"{stats.spearmanr(resid_q, resid_m).statistic:+.4f}")
    print("  -> if this collapses to ~0, expression level is the whole story.")

print("\n=== 4. Per-perturbation consistency (is it universal or a few outliers?) ===")
per = []
for i in np.unique(pert_i):
    m = pert_i == i
    if m.sum() < 50: continue
    per.append(stats.spearmanr(qv[m], mt[m]).statistic)
per = np.array(per)
print(f"  perturbations with >=50 pairs: {per.size}")
print(f"  Spearman(p_adj, correct) per pert: mean {np.nanmean(per):+.4f}  median {np.nanmedian(per):+.4f}")
print(f"  fraction POSITIVE (i.e. our ranking inverted): {(per>0).mean():.3f}")
print(f"  quartiles: {np.nanpercentile(per,[10,25,50,75,90]).round(3)}")
np.savez(os.path.join(os.path.dirname(os.path.abspath(__file__)), "r3.npz"), M=M)

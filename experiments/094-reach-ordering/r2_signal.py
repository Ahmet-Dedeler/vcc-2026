"""094 step 2 — where does sign-correctness signal live?

r1 showed the ordering channel is worth up to +0.110 overall (oracle) and that
putting our significant calls LAST beats putting them first. Both say our own
significance is anti-informative about direction. This quantifies that and
searches for a confidence key that actually predicts sign correctness.
"""
import os
import numpy as np

E092 = "~/Code/virtual-cell-challenge-2026/experiments/092-call-set-truncation"
d = np.load(os.path.join(E092, "de_tables.npz"), allow_pickle=True)
genes = d["genes"].astype(str); perts = d["perts"].astype(str)
rq, rl = d["real_padj"], d["real_lfc"]
pp, pq, pl = d["pred_p"], d["pred_padj"], d["pred_lfc"]
G, P = genes.size, perts.size
ALPHA = 0.05
gpos = {g: i for i, g in enumerate(genes)}

# pooled over all (pert, gene) pairs inside the reference's DE set
sig_acc, nsig_acc = [], []
dec_lfc = [[] for _ in range(10)]
dec_padj = [[] for _ in range(10)]
pooled_mt, pooled_absl, pooled_pq, pooled_sig = [], [], [], []

for i, t in enumerate(perts):
    ok = np.ones(G, bool)
    j = gpos.get(t)
    if j is not None: ok[j] = False
    rs = ok & (rq[i] < ALPHA)
    gsel = np.flatnonzero(rs)
    if gsel.size == 0: continue
    rlv, plv = rl[i][gsel], pl[i][gsel]
    ind = (~np.isnan(rlv)) & (rlv != 0) & (~np.isnan(plv)) & (plv != 0)
    if ind.sum() < 20: continue
    gsel = gsel[ind]; rlv, plv = rl[i][gsel], pl[i][gsel]
    mt = (np.sign(plv) == np.sign(rlv)).astype(float)
    sg = (pq[i][gsel] < ALPHA)
    absl = np.abs(plv); qv = pq[i][gsel]
    if sg.sum(): sig_acc.append(mt[sg].mean())
    if (~sg).sum(): nsig_acc.append(mt[~sg].mean())
    # deciles WITHIN each perturbation so pert-level difficulty cannot confound
    for arr, box, rev in ((absl, dec_lfc, True), (qv, dec_padj, False)):
        r = np.argsort(np.argsort(-arr if rev else arr))
        b = (r * 10 // r.size).clip(0, 9)
        for k in range(10):
            if (b == k).any(): box[k].append(mt[b == k].mean())
    pooled_mt.append(mt); pooled_absl.append(absl); pooled_pq.append(qv); pooled_sig.append(sg)

print("=== Sign accuracy on the reference's DE genes, split by OUR significance call ===")
print(f"  genes we DO  call significant : {np.mean(sig_acc):.4f}   (n perts {len(sig_acc)})")
print(f"  genes we DON'T call significant: {np.mean(nsig_acc):.4f}   (n perts {len(nsig_acc)})")
print(f"  difference: {np.mean(sig_acc)-np.mean(nsig_acc):+.4f}")
from scipy import stats
common = min(len(sig_acc), len(nsig_acc))
t, pv = stats.ttest_rel(sig_acc[:common], nsig_acc[:common])
print(f"  paired t = {t:+.2f}, p = {pv:.2e}\n")

print("=== Sign accuracy by WITHIN-perturbation decile ===")
print(f"{'decile':<9}{'by |pred lfc| desc':>22}{'by pred p_adj asc':>21}")
for k in range(10):
    print(f"{k+1:<9}{np.mean(dec_lfc[k]):>22.4f}{np.mean(dec_padj[k]):>21.4f}")
print("\n(decile 1 = most confident under that key. A useful key is MONOTONE DECREASING.)")

mt = np.concatenate(pooled_mt); absl = np.concatenate(pooled_absl)
qv = np.concatenate(pooled_pq); sg = np.concatenate(pooled_sig)
print(f"\npooled pairs {mt.size},  overall sign accuracy {mt.mean():.4f}")
print(f"  Spearman(|pred lfc|, correct)  = {stats.spearmanr(absl, mt).statistic:+.4f}")
print(f"  Spearman(pred p_adj, correct)  = {stats.spearmanr(qv, mt).statistic:+.4f}")
print(f"  point-biserial(our sig, correct) = {stats.pointbiserialr(sg, mt).statistic:+.4f}")

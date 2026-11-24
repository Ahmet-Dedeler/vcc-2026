"""181 — independent reproduction of 177's ordering counterfactual and the 22:40 debench
correction, from the artifacts already on disk.

Known-answer rows (rule 10, and each one would read differently if the thing were absent):
  KA1  reach(production order) on a15's h1-test DE table          = 0.18283   (157/171/172/177)
  KA2  reach(key=|H|)                                              = 0.23533   (INBOX-4, 171, 172, 177)
  KA3  reach(sig-first, then |H| inside)                           = 0.21914
  KA4  reach(ORACLE)                                               = 0.64553
  KA5  median n_pred for a15                                       = 6034/6035
"""
import numpy as np, sys
sys.path.insert(0, 'experiments/153-de-bench')
from debench import Surface

DE = 'experiments/177-resample-emitter/de177.%s.npz'
S = Surface(DE % 'KA_a15_raw')

eff = np.load('experiments/169-sharpening-six-member/work/h1test-effect.npz', allow_pickle=True)
gi = {g: i for i, g in enumerate(eff['genes'])}
ti = {t: i for i, t in enumerate(eff['targets'])}
gcol = np.array([gi[g] for g in S.genes]); trow = np.array([ti[t] for t in S.perts])
H = eff['H'][np.ix_(trow, gcol)]
print(f'axis alignment: {len(gcol)}/{len(S.genes)} genes, {len(trow)}/{len(S.perts)} targets')

absH = np.abs(H)
print()
print('KA1 reach(production)      =', round(S.reach(), 5), ' want 0.18283')
print('KA2 reach(key=|H|)         =', round(S.reach(key=absH), 5), ' want 0.23533')
# sig-first then |H| inside: key = sig*BIG + rank(|H|)
r = np.zeros_like(absH)
for p in range(absH.shape[0]):
    r[p] = absH[p] / (absH[p].max() + 1e-30)
key_sigH = r + S.pred_sig.astype(float) * 10.0
print('KA3 reach(sig-first,|H|)   =', round(S.reach(key=key_sigH), 5), ' want 0.21914')
orc = (np.sign(S.pred_lfc) == np.sign(S.real_lfc)).astype(float)
print('KA4 reach(ORACLE)          =', round(S.reach(key=orc), 5), ' want 0.64553')
nr, npd = S.counts()
print('KA5 median n_pred          =', int(np.median(npd)), ' want 6034/6035')
print('    sign accuracy on pool  =', round(S.sign_accuracy(), 5), ' 169 measured 0.53435')
print()

# ---- the six-member debench table for the arms 177 priced
ARMS = ['KA_a15_raw','base_mf','base_mf_b','hy_kt8000','hy_kt4000','hy_kt2000','rs_all','quasi100']
print(f"{'arm':14} {'n_pred':>7} {'reach':>8} {'fid':>8} {'prec':>8} {'jac':>8} {'nmae':>8} {'recall':>7}")
res = {}
for a in ARMS:
    try:
        s = Surface(DE % a)
    except FileNotFoundError:
        continue
    nr_, np_ = s.counts()
    # precision = mean_p k/n_pred  (no n_real; rule 9)
    prec = []
    rec = []
    for p in range(len(s.perts)):
        sel = np.flatnonzero(s.pred_sig[p] & s.adj[p])
        if sel.size:
            k = int((np.sign(s.pred_lfc[p, sel]) == np.sign(s.real_lfc[p, sel])).sum())
            prec.append(k/sel.size)
        a_ = s.real_sig[p]; b_ = s.pred_sig[p]
        if a_.sum(): rec.append(int((a_ & b_).sum())/int(a_.sum()))
    row = dict(n=int(np.median(np_)), reach=s.reach(), fid=s.fid(), prec=float(np.mean(prec)),
               jac=s.jac(), nmae=s.nmae(), rec=float(np.mean(rec)))
    res[a] = row
    print(f"{a:14} {row['n']:7d} {row['reach']:8.5f} {row['fid']:8.5f} {row['prec']:8.5f} "
          f"{row['jac']:8.5f} {row['nmae']:8.5f} {row['rec']:7.4f}")

print()
print('--- the 22:40 correction: hy_kt8000 vs base_mf on h1-test REAL n_real ---')
b, h = res['base_mf'], res['hy_kt8000']
# 177/179 conversion: overall = d_raw / span / 6, spans from the board fit (181 board181.py)
SPAN = dict(reach=1/1.12927, fid=1/3.32482, jac=1/2.72081, nmae=1/1.64787)
for m, lab in [('reach','reach'), ('fid','fid'), ('jac','jac'), ('nmae','nmae')]:
    d = h[m]-b[m]
    sgn = -1 if m == 'nmae' else 1
    print(f'  {lab:6} {b[m]:.5f} -> {h[m]:.5f}  d_raw {d:+.5f}  d_overall {sgn*d/SPAN[m]/6:+.6f}')
print(f"  recall {b['rec']:.4f} -> {h['rec']:.4f}")
print()
# jac by n_real tercile
sB = Surface(DE % 'base_mf'); sH = Surface(DE % 'hy_kt8000')
nrv = np.array([int(sB.real_sig[p].sum()) for p in range(len(sB.perts))])
q = np.quantile(nrv, [1/3, 2/3])
print('  jac loss by n_real tercile (177 addendum: top tercile -0.042, others +0.0001/+0.0005):')
for lo, hi, nm in [(-1, q[0], 'low'), (q[0], q[1], 'mid'), (q[1], 1e9, 'high')]:
    idx = np.flatnonzero((nrv > lo) & (nrv <= hi))
    def jsub(s, idx):
        o=[]
        for p in idx:
            a_=s.real_sig[p]; b_=s.pred_sig[p]; u=int((a_|b_).sum())
            o.append(1.0 if u==0 else int((a_&b_).sum())/u)
        return float(np.mean(o))
    print(f'    {nm:5} n={len(idx):3d} n_real {nrv[idx].min()}-{nrv[idx].max():5d}  '
          f'jac {jsub(sB,idx):.5f} -> {jsub(sH,idx):.5f}  d {jsub(sH,idx)-jsub(sB,idx):+.5f}')

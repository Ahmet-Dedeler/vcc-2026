"""386 audit, row 3: INDEPENDENT re-implementation of the per-target amplitude ceiling.

Own code (not derived from mse_alloc_pi.py beyond the shared task spec):
own CSR load, own pseudobulk, own name-alignment, own oracle algebra.
Compares against PI-quoted: dS=-216 (gene-sum), G_eff~18,077, -0.012 raw,
top-10 share 58%, Spearman(design,truth norms)=0.19, D-weighted cosine 0.103.
Also judges the units conversion exactly (per-target G_p, no averaging).
"""
import json, time
import numpy as np
import h5py
import scipy.sparse as sp

t0 = time.time()
W = '~/Code/virtual-cell-challenge-2026/experiments/323-current-family-validation/work/'
z = np.load(W + 'effects.npz', allow_pickle=True)
genes = np.array([str(g) for g in z['genes']])
targets = np.array([str(t) for t in z['targets']])
D = np.asarray(z['cap15'], dtype=np.float64)          # 85 x 18533 design rows
m = json.load(open(W + 'metrics-cap15.json'))
num = m['mse_legs']['numerator']; den = m['mse_legs']['denominator']

# known-answer: legs ratio == recorded raw mse
N_all = np.array([num[t] for t in targets if t in num])
Dn_all = np.array([den[t] for t in targets if t in den])
print('KNOWN-ANSWER legs: sumN/sumD = %.7f  recorded raw mse = %.7f' % (
    N_all.sum() / Dn_all.sum(), m['raw']['mse']))

TS = 5e4
with h5py.File('~/Code/virtual-cell-challenge-2026/data/harness/h1-test.h5ad') as f:
    X = f['X']
    shape = tuple(X.attrs['shape'])
    csr = sp.csr_matrix((X['data'][:], X['indices'][:], X['indptr'][:]), shape=shape)
    def dec(x):
        return x.decode() if isinstance(x, (bytes, bytearray)) else str(x)
    hgenes = np.array([dec(g) for g in f['var']['_index']['values'][:]])
    obs = f['obs']['target_gene']
    cats = np.array([dec(c) for c in obs['categories'][:]])
    codes = obs['codes'][:]
print('h1 cells x genes:', csr.shape, 'ngenes_h1 =', len(hgenes), 'load %.0fs' % (time.time() - t0))

ctrl = [c for c in cats if 'non-targeting' in c or c == 'control' or c == 'NTC'][0]
n_ctrl = int((codes == np.where(cats == ctrl)[0][0]).sum())
P_ctrl = np.asarray(csr[codes == np.where(cats == ctrl)[0][0]].sum(axis=0)).ravel()
b_ctrl = np.log1p(TS * P_ctrl / P_ctrl.sum())
print('ctrl construct:', ctrl, 'n_ctrl =', n_ctrl)

# name alignment design -> h1
pos = {g: i for i, g in enumerate(hgenes)}
keep_mask = np.array([g in pos for g in genes])
hcol = np.array([pos[g] for g in genes[keep_mask]])
d_idx = np.where(keep_mask)[0]
print('aligned %d of %d design genes; h1 universe %d' % (len(d_idx), len(genes), len(hgenes)))

dd, tt, dt, NN, DD, Gp, TT = [], [], [], [], [], [], []
for k, t in enumerate(targets):
    if t not in set(cats.tolist()):
        continue
    row = D[k][d_idx]                                   # design, aligned coords
    P = np.asarray(csr[codes == np.where(cats == t)[0][0]].sum(axis=0)).ravel()
    delta = np.log1p(TS * P / P.sum()) - b_ctrl         # full h1 coords
    tru = delta[hcol]
    own = (genes[d_idx] == t)                           # own-target exclusion
    row = row.copy(); tru = tru.copy()
    row[own] = 0.0; tru[own] = 0.0
    dd.append(float(row @ row)); tt.append(float(tru @ tru)); dt.append(float(row @ tru))
    NN.append(num.get(t, np.nan)); DD.append(den.get(t, np.nan))
    Gp.append(float(len(hgenes) - 1))                   # harness sums all genes but own target
    TT.append(t)
dd = np.array(dd); tt = np.array(tt); dt = np.array(dt)
NN = np.array(NN); DD = np.array(DD); Gp = np.array(Gp)
n = len(TT)
cos = dt / np.sqrt(dd * tt)
nd = np.sqrt(dd); nt = np.sqrt(tt)
sumD = float(np.nansum(DD))
print('n targets: %d  harness raw mse sumN/sumD = %.5f' % (n, np.nansum(NN) / sumD))

# concentration
order = np.argsort(-DD)
cum = np.cumsum(DD[order]) / DD.sum()
for j in (5, 10, 20, 40):
    print('top-%d share of sumD: %.3f' % (j, cum[j - 1]))
print('D-weighted mean cos: %.4f  median cos: %.4f' % ((cos * DD).sum() / DD.sum(), np.median(cos)))

def spr(x, y):
    from scipy.stats import spearmanr
    return float(spearmanr(x, y)[0])
print('spearman(design norm, truth norm) = %.4f' % spr(nd, nt))

# oracle algebra
def Stot(a):
    return float(np.sum(a * a * dd - 2 * a * dt + tt))
a_glob = float(dt.sum() / dd.sum())
a_or = dt / dd
S_glob = Stot(np.full(n, a_glob)); S_or = Stot(a_or); S_1 = Stot(np.ones(n))
dS = S_or - S_glob
print('global alpha* = %.4f' % a_glob)
print('S(alpha=1)/... vs S(glob): %.1f vs %.1f; S(oracle) = %.1f; dS = %+.1f (PI -216)' % (S_1, S_glob, S_or, dS))

# units: exact conversion uses per-target G_p (per-gene means), PI used one G_eff
dS_per_gene_exact = float(np.sum(( (a_or**2) * dd - 2 * a_or * dt + tt - ((a_glob**2) * dd - 2 * a_glob * dt + tt) ) / Gp))
exact = dS_per_gene_exact / sumD
G_eff_pi = 18077.0
approx_pi = dS / (G_eff_pi * sumD)
G_eff_implied = dS / (exact * sumD)   # effective G my exact calc implies
print('sumD = %.6f' % sumD)
print('EXACT  per-target-G conversion: delta raw mse = %+.6f' % exact)
print('PI-STYLE G=18077 conversion   : delta raw mse = %+.6f (PI ~ -0.012)' % approx_pi)
print('my implied effective G = %.0f (PI 18077); n_aligned-1 = %d; harness G_p = %.0f' % (
    G_eff_implied, len(d_idx) - 1, Gp[0]))
print('scaled x(1/0.955)=x%.3f -> %+.6f ; overall /6 -> %+.6f' % (
    1 / 0.955, exact / 0.955, -exact / 0.955 / 6))
print('oracle a_p quantiles [.1,.25,.5,.75,.9]:', np.round(np.quantile(a_or, [.1, .25, .5, .75, .9]), 3).tolist())
print('done %.0fs' % (time.time() - t0))

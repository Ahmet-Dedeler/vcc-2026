"""Independent re-derivation of 344's context-A claims. h5py + numpy + scipy.sparse
only. Wilcoxon z, BH and Spearman hand-rolled (math.erf for the normal tail)."""
import h5py
import json
import math
import numpy as np
from scipy import sparse

F = '~/Code/virtual-cell-challenge-2026/data/challenge/context_A.h5ad'
OUT = {}
f = h5py.File(F, 'r')

n = f['X'].attrs['shape'][0]
indptr = f['X']['indptr'][:]
indices = f['X']['indices'][:]
data = f['X']['data'][:].astype(np.float64)
X = sparse.csr_matrix((data, indices, indptr), shape=tuple(f['X'].attrs['shape']))
libs = np.asarray(X.sum(axis=1)).ravel()          # library sizes
cats = [c.decode() if isinstance(c, bytes) else c
        for c in f['obs']['ntc_id']['categories'][:]]
codes = f['obs']['ntc_id']['codes'][:]
varnames = [c.decode() if isinstance(c, bytes) else c
            for c in f['var']['_index']['values'][:]]
OUT['n_cells'] = int(n)
OUT['n_constructs'] = len(set(codes.tolist()))
constructs = sorted(set(codes.tolist()))
OUT['cells_per_construct'] = sorted(int((codes == c).sum()) for c in constructs)
invlib = 1e6 / libs

# control mean CPM per gene (CSC column streaming, one gene at a time)
Xc = X.tocsc()
G = X.shape[1]
meancpm = np.empty(G)
for j in range(G):
    s, e = Xc.indptr[j], Xc.indptr[j + 1]
    meancpm[j] = (Xc.data[s:e] * invlib[Xc.indices[s:e]]).sum() / n
uni = np.where(meancpm > 5)[0]
OUT['universe_5cpm'] = int(len(uni))


def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def ranks_avg(v):
    """Average ranks (1-based) with tie handling, hand-rolled via argsort."""
    o = np.argsort(v, kind='mergesort')
    sv = v[o]
    r = np.empty_like(v, dtype=np.float64)
    i, m = 0, len(v)
    while i < m:
        j = i + 1
        while j < m and sv[j] == sv[i]:
            j += 1
        r[o[i:j]] = (i + 1 + j) / 2.0
        i = j
    return r


def spearman(x, y):
    rx, ry = ranks_avg(x), ranks_avg(y)
    rx -= rx.mean(); ry -= ry.mean()
    return float(rx.dot(ry) / math.sqrt(rx.dot(rx) * ry.dot(ry)))


N2 = n                      # pool size 18,400 (construct cells included, literal reading)
N1 = 400
N = N1 + N2
n1n2 = N1 * N2
mu = n1n2 / 2.0

Z = np.empty((len(constructs), len(uni)))   # z per construct x universe gene
P = np.empty_like(Z)
for k, j in enumerate(uni):
    s, e = Xc.indptr[j], Xc.indptr[j + 1]
    col = np.zeros(n)
    col[Xc.indices[s:e]] = Xc.data[s:e] * invlib[Xc.indices[s:e]]
    v = np.log1p(col)
    psorted = np.sort(v)
    pu, pc = np.unique(psorted, return_counts=True)
    pod = dict(zip(pu.tolist(), pc.tolist()))
    for ci, c in enumerate(constructs):
        vc = v[codes == c]
        lo = np.searchsorted(psorted, vc, side='left')
        hi = np.searchsorted(psorted, vc, side='right')
        U = float(np.sum(lo + 0.5 * (hi - lo)))
        # tie correction over the combined 18,800 multiset
        cu, cc = np.unique(vc, return_counts=True)
        tcorr = 0.0
        for val, cnt in zip(cu.tolist(), cc.tolist()):
            t = pod.get(val, 0) + cnt
            tcorr += t ** 3 - t
        var = n1n2 / 12.0 * ((N + 1) - tcorr / (N * (N - 1)))
        z = (U - mu) / math.sqrt(var)
        Z[ci, k] = z
        P[ci, k] = 2.0 * norm_sf(abs(z))

# BH step-up per construct, alpha 0.05 (hand-rolled)
ALPHA = 0.05
m = len(uni)
sig_counts = []
for ci in range(len(constructs)):
    o = np.argsort(P[ci])
    sp = P[ci][o]
    thresh = ALPHA * np.arange(1, m + 1) / m
    ok = np.where(sp <= thresh)[0]
    kmax = int(ok.max()) + 1 if len(ok) else 0
    sig = set(o[:kmax].tolist()) if kmax else set()
    sig_counts.append(len(sig))
sig_counts = np.array(sig_counts)
OUT['per_construct_hits'] = [
    {'construct': cats[c], 'n_sig': int(sig_counts[ci])}
    for ci, c in enumerate(constructs)]
OUT['constructs_with_any_hit'] = int((sig_counts > 0).sum())
OUT['median_sig_genes'] = float(np.median(sig_counts))
OUT['sig_counts_sorted'] = sorted(sig_counts.tolist())

key = np.abs(Z).mean(axis=0)
OUT['spearman_meanAbsZ_vs_ctrlCPM'] = spearman(key, meancpm[uni])
OUT['spearman_meanAbsZ_vs_logCPM'] = spearman(key, np.log(meancpm[uni]))

json.dump(OUT, open('~/Code/virtual-cell-challenge-2026/experiments/354-audit-343-344/work/audit344.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in OUT.items() if k != 'sig_counts_sorted'}, indent=1))
print('sig_counts_sorted =', OUT['sig_counts_sorted'])

"""603b — dose curve of the UP-only rank-neutral magnitude rule on H1: |lfc| <- max(|ol|, min(K * MDE_g, cap)) for genes
delivered UP; K in 0.5..1.5, cap in {inf, 1.0, 0.5} log2. Same state and scoring as ceiling603.py (3 seeds, f 1 / 0.35)."""
import sys, json
from pathlib import Path
import numpy as np
ROOT = Path('~/Code/virtual-cell-challenge-2026')
sys.path.insert(0, str(ROOT / 'experiments/521-strength-nreal')); sys.path.insert(0, str(ROOT / 'experiments/603-rank-neutral-ceiling'))
import lib521 as L
ZS = ROOT / 'experiments/435-zs-on-coarse/work/h1/stats-base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz'
import anndata as ad, scipy.sparse as sp
S0 = L.load_state(0)
a = ad.read_h5ad(ROOT / 'data/harness/h1-test.h5ad', backed='r')
ci = np.where(a.obs['target_gene'].values == 'non-targeting')[0]
X = sp.csr_matrix(a.X[ci], dtype=np.float64); mc_all = np.asarray(X.mean(0)).ravel()
pos = {g: i for i, g in enumerate(np.asarray(a.var_names))}
mc = np.array([mc_all[pos[g]] if g in pos else np.nan for g in S0['genes']])
mde = 0.24 / np.sqrt(np.maximum(mc, 1e-3))
grid = [(K, cap) for cap in (np.inf, 1.0, 0.5) for K in (0.5, 0.75, 1.0, 1.25, 1.5)]
acc = {(f, K, cap): [] for f in (1.0, 0.35) for K, cap in grid}
for seed in (0, 1, 2):
    S = L.load_state(seed); st = np.load(str(ZS).format(seed=seed)); ol = np.asarray(st['lfc'], float); rl = S['rl']
    up = ol > 0
    for f in (1.0, 0.35):
        pool = L.pools(S, f); base = np.nanmean(L.nmae_rows(ol, pool, rl))
        for K, cap in grid:
            mag = np.maximum(np.abs(ol), np.minimum(K * mde[None, :], cap))
            new = np.where(up, mag, ol)
            acc[(f, K, cap)].append(float(np.nanmean(L.nmae_rows(new, pool, rl)) - base))
out = {}
print('UP-only rank-neutral dose, mean dNMAE over 3 seeds [overall-equivalent]:')
for f in (1.0, 0.35):
    for cap in (np.inf, 1.0, 0.5):
        line = []
        for K in (0.5, 0.75, 1.0, 1.25, 1.5):
            d = float(np.mean(acc[(f, K, cap)])); sd = float(np.std(acc[(f, K, cap)]))
            out[f'f{f}_K{K}_cap{cap}'] = dict(dnmae=d, sd=sd, overall=-d * L.INV['nmae'] / 6)
            line.append(f'K{K}: {d:+.4f} [{-d * L.INV["nmae"] / 6:+.4f}]')
        print(f'  f {f:4} cap {cap:>4}: ' + '  '.join(line))
json.dump(out, open(ROOT / 'experiments/603-rank-neutral-ceiling/work/dose603.json', 'w'), indent=1)

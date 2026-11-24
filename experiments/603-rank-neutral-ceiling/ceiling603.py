"""603 — analytic ceiling of the rank-neutral nmae lever on H1 (council sweep #1, owner of the arm: d4cd5f).

nmae reads the arithmetic-mean lfc of the reference's gate genes; calls / reach read Wilcoxon ranks. If the emitter can
move a gene's mean without moving its ranks (UP: an outlier cell above every control value; DOWN: bounded), then the
lfc of every delivered gene can be re-sized with calls, order, signs, pds and mse pinned. This script holds 446bj2z's
DELIVERED H1 state (435's production-emitter stats WITH the zs pass, 3 seeds; calls = 521's reconstructed C = head|block)
and changes ONLY lfc magnitudes, then scores nmae exactly (521's nmae_rows) at f 1 and f 0.35 (PowerSurface).

Rules applied to genes delivered UP (ol > 0) — the free direction — and, separately as an unattainable bound, DOWN:
  mdeK  : |lfc| <- max(|ol|, K * MDE_g), MDE_g = 0.24 / sqrt(mean control count per cell), K in 1, 2, 4
  key   : |lfc| <- max(|ol|, |cap15 key| / ln 2)  (the bulk's own donor-predicted size, log-pseudobulk -> log2)
  oracle: |lfc| <- |reference lfc|                 (upper bound: truth magnitudes at OUR signs)
Truth is used only to score. No production file is read for writing; outputs go to this directory.
"""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path('~/Code/virtual-cell-challenge-2026')
sys.path.insert(0, str(ROOT / 'experiments/521-strength-nreal'))
sys.path.insert(0, str(ROOT / 'src'))
import lib521 as L

ZS = ROOT / 'experiments/435-zs-on-coarse/work/h1/stats-base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz'
OUT = ROOT / 'experiments/603-rank-neutral-ceiling/work'; OUT.mkdir(parents=True, exist_ok=True)

def mean_count_per_cell(genes):
    import anndata as ad, scipy.sparse as sp
    a = ad.read_h5ad(ROOT / 'data/harness/h1-test.h5ad', backed='r')
    ci = np.where(a.obs['target_gene'].values == 'non-targeting')[0]
    X = sp.csr_matrix(a.X[ci], dtype=np.float64)
    mc = np.asarray(X.mean(0)).ravel()
    pos = {g: i for i, g in enumerate(np.asarray(a.var_names))}
    return np.array([mc[pos[g]] if g in pos else np.nan for g in genes])

res = {}
mc = None
for seed in (0, 1, 2):
    S = L.load_state(seed)
    st = np.load(str(ZS).format(seed=seed))
    assert list(st['targets'].astype(str)) == list(S['targets'])
    ol = np.asarray(st['lfc'], float)                 # delivered lfc, zs form (446bj2z's DE construction)
    if mc is None:
        mc = mean_count_per_cell(S['genes'])
    mde = 0.24 / np.sqrt(np.maximum(mc, 1e-3))        # Poisson MDE law (425)
    keylfc = S['absk'] / np.log(2.0)
    rl = S['rl']
    up = (ol > 0); dn = (ol < 0)
    rules = {}
    for K in (1, 2, 4):
        rules[f'mde{K}'] = np.maximum(np.abs(ol), K * mde[None, :])
    rules['key'] = np.maximum(np.abs(ol), keylfc)
    rules['oracle'] = np.where(np.isfinite(rl), np.abs(rl), np.abs(ol))
    for f in (1.0, 0.35):
        pool = L.pools(S, f)
        base = np.nanmean(L.nmae_rows(ol, pool, rl))
        gate = pool & np.isfinite(rl)
        w = np.abs(np.where(gate, rl, 0))
        row = dict(base_nmae=float(base),
                   up_share_genes=float((gate & up).sum() / gate.sum()),
                   up_share_mass=float(w[gate & up].sum() / w[gate].sum()),
                   acc_up=float((np.sign(rl[gate & up]) > 0).mean()),
                   acc_dn=float((np.sign(rl[gate & dn]) < 0).mean()),
                   acc_all=float((np.sign(rl[gate & (up | dn)]) == np.sign(ol[gate & (up | dn)])).mean()))
        for name, mag in rules.items():
            for side, m in (('UP', up), ('DOWN', dn), ('BOTH', up | dn)):
                new = np.where(m, np.sign(ol) * mag, ol)
                v = np.nanmean(L.nmae_rows(new, pool, rl))
                row[f'{name}_{side}'] = float(v - base)
        res[f's{seed}_f{f}'] = row
        print(f'seed {seed} f {f}: base nmae {base:.5f} | UP share genes {row["up_share_genes"]:.3f} mass {row["up_share_mass"]:.3f} '
              f'acc_up {row["acc_up"]:.3f} acc_dn {row["acc_dn"]:.3f}', flush=True)
        for name in rules:
            print(f'    {name:7s} dNMAE  UP {row[name+"_UP"]:+.5f}   DOWN(bound) {row[name+"_DOWN"]:+.5f}   BOTH {row[name+"_BOTH"]:+.5f}')
json.dump(res, open(OUT / 'ceiling603.json', 'w'), indent=1)
# summary over seeds, in overall units (192's nmae slope via lib521's INV: overall = -dNMAE * 1.64787 / 6)
print('\nMEAN over 3 seeds (overall-equivalent = -dNMAE * 1.64787 / 6):')
for f in (1.0, 0.35):
    rows = [res[f's{s}_f{f}'] for s in (0, 1, 2)]
    for key in [k for k in rows[0] if k.endswith(('_UP', '_DOWN', '_BOTH'))]:
        d = np.mean([r[key] for r in rows])
        print(f'  f {f:4}  {key:14s} dNMAE {d:+.5f}  overall {-d * L.INV["nmae"] / 6:+.5f}')

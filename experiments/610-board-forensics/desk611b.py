# 611 desk check, part b: H1 prices of the nmae board pairs at 521's four regimes, each run scored with its OWN realized
# order (scorer's rule: sig first, q, p, |lfc|), 521's exact member functions. Pairs: zs (435 coarse vs coarse+zs, seeds
# 0-2) and 433a's coarse MDE (426 production-emitter no-MDE vs 435 coarse, seed 0; bulk sums must be equal).
import sys, json
import numpy as np
sys.path.insert(0, 'experiments/521-strength-nreal')
import lib521 as L
import table521 as TB
S = L.load_state(0)
own, rl = S['own'], S['rl']
adj = L.adjud(S)
H435 = 'experiments/435-zs-on-coarse/work/h1/stats-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse{z}.npz'
H426 = 'experiments/426-mde-magnitude/work/h1/stats-base_s0_bps_prodem.npz'
def load(p):
    z = np.load(p, allow_pickle=True)
    assert list(z['targets'].astype(str)) == list(S['targets'])
    return dict(q=np.asarray(z['q'], float), p=np.asarray(z['p'], float), lfc=np.asarray(z['lfc'], float), sums=z['sums'])
def score(run, reg):
    pool, _ = TB.regime_pool(S, reg)
    call = (run['q'] < 0.05) & ~own
    fr = L.fid_rows(call, run['lfc'], pool, adj, rl); jr = L.jac_rows(call, pool); nr = L.nmae_rows(run['lfc'], pool, rl)
    St = dict(S); St['oq'], St['op'], St['ol'] = run['q'], run['p'], run['lfc']
    rr = np.array([L.reach_row(L.realized_order_pos(St, i), pool[i], adj[i], run['lfc'][i], rl[i]) for i in range(len(S['targets']))])
    return dict(fid=np.nanmean(fr), jac=np.mean(jr), reach=np.nanmean(rr), nmae=np.nanmean(nr))
out = {}
pairs = {'zs (446bj2z-446bj2)': [(H435.format(s=s, z='_dnk250zs'), H435.format(s=s, z='')) for s in (0, 1, 2)],
         'coarse MDE (433a-417e)': [(H435.format(s=0, z=''), H426)]}
board = {'zs (446bj2z-446bj2)': dict(fid=0.0, jac=0.0, reach=0.0, nmae=-0.002179),
         'coarse MDE (433a-417e)': dict(fid=0.514327-0.514323, jac=0.040808-0.040829, reach=0.286842-0.284445, nmae=0.919618-0.936405)}
for name, prs in pairs.items():
    runs = [(load(a), load(b)) for a, b in prs]
    for a, b in runs:
        print('%s: bulk sums equal: %s' % (name, bool(np.array_equal(a['sums'], b['sums']))))
    out[name] = {}
    print(name, '| board raw delta:', ' '.join('%s %+.5f' % (m, v) for m, v in board[name].items()))
    for reg in TB.REGIMES:
        ds = []
        for a, b in runs:
            sa, sb = score(a, reg), score(b, reg)
            ds.append({m: sa[m] - sb[m] for m in sa})
        d = {m: float(np.mean([x[m] for x in ds])) for m in ds[0]}
        out[name][reg] = d
        rat = {m: (board[name][m] / d[m]) if abs(d[m]) > 1e-6 else float('nan') for m in d}
        print('  H1 %-4s delta: %s | board/H1: %s' % (reg, ' '.join('%s %+.5f' % (m, d[m]) for m in d), ' '.join('%s %+.2f' % (m, rat[m]) for m in d)))
json.dump(out, open('experiments/610-board-forensics/work/desk611b.json', 'w'), indent=1)

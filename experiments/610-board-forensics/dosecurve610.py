# 610: the PRIZE of rank-neutral mean delivery on uncalled genes (H1, 435 coarse+zs = 446bj2z's construction, seeds 0-2).
# Calls, p-values and order are held fixed (that is what rank-neutral delivery means), so fid/jac/reach do not move; only
# nmae does. Scale the delivered lfc of UNCALLED genes by c (sign kept), with three sign sources for them: delivered (ol),
# key (521 cap15), and a key/CLEAN4 hybrid (CLEAN4 where |key| is in the bottom half of the target's nonzero keys).
# Board expectation = H1 x 1.2-1.6 (611). Analysis only.
import sys, json
import numpy as np
sys.path.insert(0, 'experiments/521-strength-nreal')
import lib521 as L
import table521 as TB
c4 = np.load('experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz', allow_pickle=True)
res = {}
for seed in (0, 1, 2):
    S = L.load_state(seed)
    z = np.load('experiments/435-zs-on-coarse/work/h1/stats-base_s%d_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz' % seed)
    q, ol = np.asarray(z['q'], float), np.asarray(z['lfc'], float)
    own, rl, key = S['own'], S['rl'], S['key']
    unc = (q >= 0.05) & ~own
    g = list(S['genes']); c4g = {x: i for i, x in enumerate(c4['genes'])}
    pd4 = np.array([c4['P_down'][c4g[x]] if x in c4g else np.nan for x in g])
    s4 = np.where(np.isfinite(pd4), np.where(pd4 > 0.5, -1.0, 1.0), 0.0)
    T, G = ol.shape
    q50 = np.zeros((T, G), bool)
    for i in range(T):
        a = np.abs(key[i]); nz = a > 0
        if nz.sum(): q50[i, nz] = a[nz] < np.median(a[nz])
    mag = np.abs(ol)
    srcs = {'delivered': np.sign(ol), 'key': np.sign(key), 'hybrid': np.where(q50 & (s4 != 0), s4, np.sign(key))}
    for reg in ('1.0', '0.35'):
        pool, _ = TB.regime_pool(S, reg)
        base = np.nanmean(L.nmae_rows(ol, pool, rl))
        for sname, sg in srcs.items():
            for c in (1.0, 1.5, 2.0, 3.0, 4.0, 6.0):
                l2 = ol.copy()
                l2[unc] = sg[unc] * mag[unc] * c
                # genes with zero delivered magnitude stay zero (the MDE skipped them); keep that
                v = np.nanmean(L.nmae_rows(l2, pool, rl)) - base
                res.setdefault((reg, sname, c), []).append(v)
print('H1 nmae delta vs 446bj2z construction (mean of 3 seeds); board expectation = x1.2-1.6; overall = -delta x 1.64787/6')
for reg in ('1.0', '0.35'):
    print('regime f=%s' % reg)
    for sname in ('delivered', 'key', 'hybrid'):
        row = ['c%.1f %+.4f' % (c, np.mean(res[(reg, sname, c)])) for c in (1.0, 1.5, 2.0, 3.0, 4.0, 6.0)]
        best = min((np.mean(res[(reg, sname, c)]), c) for c in (1.0, 1.5, 2.0, 3.0, 4.0, 6.0))
        print('  %-9s %s | best c=%.1f -> board %+.4f..%+.4f raw, overall %+.4f..%+.4f' % (sname, ' '.join(row), best[1], best[0]*1.2, best[0]*1.6, -best[0]*1.2*1.64787/6, -best[0]*1.6*1.64787/6))
json.dump({'%s|%s|%s' % k: v for k, v in res.items()}, open('experiments/610-board-forensics/work/dosecurve610.json', 'w'))

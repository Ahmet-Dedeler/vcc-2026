# 612 diag b: where does idkg differ from the plain hook (c4h) on H1's real emitter (seed 0)? Per-target nmae contribution
# of (idkg - c4h), split by gene class: called / uncalled x (strong-key | weak-key) x (CLEAN4-covered at tau .05 | not).
import sys
import numpy as np
sys.path.insert(0, 'experiments/521-strength-nreal')
import lib521 as L
import table521 as TB
S = L.load_state(0)
own, rl, key = S['own'], S['rl'], S['key']
T = 'stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz'
runs = {n: np.load(p + T) for n, p in (('c4h', 'experiments/537-arm-pipeline/work/h1_cpm+c4h/'),
                                        ('idkg', 'experiments/612-key-gated-hook/work/h1_idkg/'),
                                        ('base', 'experiments/435-zs-on-coarse/work/h1/'))}
lf = {n: np.asarray(r['lfc'], float) for n, r in runs.items()}
call = (np.asarray(runs['base']['q'], float) < 0.05) & ~own
unc = ~call & ~own
a = np.abs(key); strong = np.zeros_like(unc)
for i in range(a.shape[0]):
    nz = a[i] > 0
    if nz.any(): strong[i] = nz & (a[i] >= np.quantile(a[i][nz], 0.5))
cp = np.load('experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz', allow_pickle=True)
names = list(cp['bank_names'].astype(str)); Pc = cp['P_banks'].astype(float)[[names.index(b) for b in ('k562_gwps','k562_essential','nadig_jurkat','nadig_hepg2')]].mean(0)
pos = {g: j for j, g in enumerate(cp['genes'].astype(str))}
Pg = np.array([Pc[pos[g]] if g in pos else np.nan for g in S['genes']])
c4cov = np.isfinite(Pg) & (np.abs(Pg - 0.5) >= 0.05)
for f in ('1.0', '0.35'):
    pool, _ = TB.regime_pool(S, f)
    gate = pool & np.isfinite(rl)
    cls = {'called': gate & call,
           'unc weak C4cov': gate & unc & ~strong & c4cov[None, :], 'unc weak notC4': gate & unc & ~strong & ~c4cov[None, :],
           'unc strong C4cov': gate & unc & strong & c4cov[None, :], 'unc strong notC4': gate & unc & strong & ~c4cov[None, :]}
    out = {c: [] for c in cls}; tot = []
    for i in range(rl.shape[0]):
        g = gate[i]
        if g.sum() < 10: continue
        den = np.abs(rl[i, g]).sum()
        for c, m in cls.items():
            mm = m[i]
            out[c].append((np.abs(lf['idkg'][i, mm] - rl[i, mm]).sum() - np.abs(lf['c4h'][i, mm] - rl[i, mm]).sum()) / den)
        tot.append(sum(out[c][-1] for c in cls))
    print(f'f={f}: d nmae (idkg - c4h) = {np.mean(tot):+.5f} = ' + ' + '.join(f'{c} {np.mean(v):+.5f}' for c, v in out.items()))
    for c, m in cls.items():
        if 'unc' in c:
            mm = m
            acc = {n: (np.sign(lf[n][mm]) == np.sign(rl[mm])).mean() for n in ('base', 'c4h', 'idkg')}
            print(f'    {c:17s} n {int(mm.sum()):6d}  sign acc base {acc["base"]:.3f} c4h {acc["c4h"]:.3f} idkg {acc["idkg"]:.3f}')

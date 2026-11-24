# 612 diagnosis: what does the REAL emitter deliver on uncalled gate genes under base / plain hook (c4h) / key-gated (c4kg)?
# seed 0, zs form; h1-test truth at f 1 and f 0.35. By |key| half: delivered sign accuracy, mean |lfc|, flip rate vs base.
import sys
import numpy as np
sys.path.insert(0, 'experiments/521-strength-nreal')
import lib521 as L
import table521 as TB
S = L.load_state(0)
own, rl, key = S['own'], S['rl'], S['key']
T = 'stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz'
runs = {n: np.load(p + T) for n, p in (('base', 'experiments/435-zs-on-coarse/work/h1/'),
                                        ('c4h', 'experiments/537-arm-pipeline/work/h1_cpm+c4h/'),
                                        ('c4kg', 'experiments/612-key-gated-hook/work/h1_c4kg/'))}
lf = {n: np.asarray(r['lfc'], float) for n, r in runs.items()}
unc = (np.asarray(runs['base']['q'], float) >= 0.05) & ~own
a = np.abs(key); strong = np.zeros_like(unc)
for i in range(a.shape[0]):
    nz = a[i] > 0
    if nz.any(): strong[i] = nz & (a[i] >= np.quantile(a[i][nz], 0.5))
for f in ('1.0', '0.35'):
    pool, _ = TB.regime_pool(S, f)
    gate = pool & np.isfinite(rl) & (rl != 0) & unc
    print(f'regime f={f}: uncalled gate genes {gate.sum()} (strong-key half {int((gate & strong).sum())})')
    for half, m in (('weak-key', gate & ~strong), ('strong-key', gate & strong)):
        line = f'  {half:10s}'
        for n in ('base', 'c4h', 'c4kg'):
            x = lf[n][m]; t = rl[m]
            acc = (np.sign(x) == np.sign(t)).mean(); mag = np.abs(x).mean()
            l1 = np.abs(x - t).sum() / np.abs(t).sum()
            flip = (np.sign(x) != np.sign(lf['base'][m])).mean()
            line += f' | {n}: acc {acc:.3f} |lfc| {mag:.4f} L1 {l1:.4f} flip {flip:.3f}'
        print(line)

# per-target nmae decomposition of (c4kg - c4h) into gene classes (nmae = mean over targets of sum|x-t| / sum|t| on the gate)
call = (np.asarray(runs['base']['q'], float) < 0.05) & ~own
for f in ('1.0', '0.35'):
    pool, _ = TB.regime_pool(S, f)
    gate = pool & np.isfinite(rl)
    classes = {'called': gate & call, 'unc weak': gate & unc & ~strong, 'unc strong': gate & unc & strong}
    contrib = {c: [] for c in classes}; tot = []
    for i in range(rl.shape[0]):
        g = gate[i]
        if g.sum() < 10: continue
        den = np.abs(rl[i, g]).sum()
        if not np.isfinite(den) or den == 0: continue
        for c, m in classes.items():
            mm = m[i]
            d = (np.abs(np.nan_to_num(lf['c4kg'][i, mm]) - rl[i, mm]).sum() - np.abs(np.nan_to_num(lf['c4h'][i, mm]) - rl[i, mm]).sum()) / den
            contrib[c].append(d)
        tot.append(sum(contrib[c][-1] for c in classes))
    print(f'f={f}: d nmae (c4kg - c4h) = {np.mean(tot):+.5f} = ' + ' + '.join(f'{c} {np.mean(v):+.5f}' for c, v in contrib.items()))
    same_called = np.array_equal(lf['c4kg'][call], lf['c4h'][call])
    print(f'   called genes delivered lfc identical between c4kg and c4h: {same_called}; '
          f'max |diff| on called {np.nanmax(np.abs(lf["c4kg"][call] - lf["c4h"][call])):.2e}')

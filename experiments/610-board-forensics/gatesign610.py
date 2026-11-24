# 610 (2026-10-02): on h1-test's GATE genes (the reference's significant set, what nmae and reach score), compare sign
# sources on the SAME genes: the delivered lfc (ol, 435's production emitter), our target-specific key (521's cap15
# effect key = the whitened donor composition), and CLEAN4 (391's P_down). Split by our delivery class: head (top-400
# |key|), CPM block, complement (uncalled). Weighted by |true lfc| too (nmae's weight). Analysis only.
import sys
import numpy as np
sys.path.insert(0, 'experiments/521-strength-nreal')
import lib521 as L
S = L.load_state(0)
own, rl, ol, key = S['own'], S['rl'], S['ol'], S['key']
genes = list(S['genes'])
c4 = np.load('experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz', allow_pickle=True)
c4g = {g: i for i, g in enumerate(c4['genes'])}
pd4 = np.array([c4['P_down'][c4g[g]] if g in c4g else np.nan for g in genes])
s4 = np.where(np.isfinite(pd4), np.where(pd4 > 0.5, -1.0, 1.0), np.nan)
T, G = rl.shape
S4 = np.tile(s4, (T, 1))
fin = np.isfinite(rl) & (rl != 0) & ~own
call = (S['oq'] < 0.05) & ~own
classes = {'head': S['head'], 'block': S['block'] & ~S['head'], 'complement': ~S['C']}
def acc(src, m, w=None):
    ok = m & np.isfinite(src) & (src != 0)
    agree = (np.sign(src[ok]) == np.sign(rl[ok]))
    if w is None: return agree.mean(), ok.sum()
    ww = np.abs(rl[ok]); return (agree * ww).sum() / ww.sum(), ok.sum()
for f in (1.0, 0.35):
    pool = L.pools(S, f) & fin
    print('regime f=%s: gate pairs %d' % (f, pool.sum()))
    for cname, cm in [('ALL gate', np.ones_like(pool))] + list(classes.items()):
        m = pool & cm
        both = m & np.isfinite(S4) & (key != 0)   # genes where all three sources carry a sign
        r = [acc(src, both) for src in (ol, key, S4)]
        rw = [acc(src, both, 1) for src in (ol, key, S4)]
        print('  %-11s n=%6d (called %4.1f%%) | sign acc delivered %.3f  key %.3f  CLEAN4 %.3f | |l|-weighted: %.3f %.3f %.3f' % (
            cname, both.sum(), 100 * (both & call).sum() / max(both.sum(), 1), r[0][0], r[1][0], r[2][0], rw[0][0], rw[1][0], rw[2][0]))
    # where key and CLEAN4 disagree on gate genes
    m = pool & np.isfinite(S4) & (key != 0)
    dis = m & (np.sign(key) != S4)
    print('  key vs CLEAN4 disagree on %.1f%% of gate genes; there key right %.3f, CLEAN4 right %.3f' % (
        100 * dis.sum() / m.sum(), acc(key, dis)[0], acc(S4, dis)[0]))
    # key confidence bins (|key| quantiles within target) on complement gate genes
    comp = pool & ~S['C'] & (key != 0)
    q = np.zeros_like(key)
    for i in range(T):
        a = np.abs(key[i]); nz = a > 0
        if nz.sum(): q[i, nz] = (np.argsort(np.argsort(a[nz])) + 0.5) / nz.sum()
    for lo, hi in [(0, .5), (.5, .8), (.8, .95), (.95, 1.01)]:
        mm = comp & (q >= lo) & (q < hi)
        print('    complement gate genes, |key| quantile [%.2f,%.2f): n=%5d key acc %.3f CLEAN4 acc %.3f delivered %.3f' % (lo, hi, mm.sum(), acc(key, mm)[0], acc(S4, mm & np.isfinite(S4))[0], acc(ol, mm)[0]))

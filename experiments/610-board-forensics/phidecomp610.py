# 610 (2026-10-02, for 7772dc's critique): decompose our H1 fid as fid ~ (1-phi) p_null + phi p_DE on the shipped DE
# construction (435's production-emitter stats via 521's load_state; zs leaves calls unchanged), at 521's power regimes,
# and compare composition alternatives at matched size: our head+CPM block vs a key-ranked (target-specific) call set.
# Signs for hypothetical call sets = the emitted lfc sign (ol). Analysis only, < 2 GB.
import sys
import numpy as np
sys.path.insert(0, 'experiments/521-strength-nreal')
import lib521 as L
S = L.load_state(0)
own, rl, ol, oq = S['own'], S['rl'], S['ol'], S['oq']
adj = L.adjud(S)
call = (oq < 0.05) & ~own
T, G = call.shape
def decomp(call, pool, label, sign_src=ol):
    rows = []
    for i in range(T):
        sel = call[i] & adj[i]; npred = int(sel.sum()); nreal = int(pool[i].sum())
        if npred == 0 and nreal == 0: continue
        agree = np.sign(sign_src[i, sel]) == np.sign(rl[i, sel])
        de = pool[i, sel]
        k = int(agree.sum())
        rows.append((npred, nreal, de.mean() if npred else np.nan, agree[de].mean() if de.sum() else np.nan,
                     agree[~de].mean() if (~de).sum() else np.nan, k / max(npred, nreal),
                     (call[i] & pool[i]).sum() / max((call[i] | pool[i]).sum(), 1)))
    a = np.array(rows, float)
    print('  %-34s n_pred %6.0f | phi %.4f  p_DE %.4f  p_null %.4f | fid %.4f  jac %.4f' % (
        label, np.nanmedian(a[:, 0]), np.nanmean(a[:, 2]), np.nanmean(a[:, 3]), np.nanmean(a[:, 4]), np.nanmean(a[:, 5]), np.nanmean(a[:, 6])))
for f in (1.0, 0.5, 0.35):
    pool = L.pools(S, f)
    print('regime f=%s: n_real median %d, mean %.0f' % (f, np.median(pool.sum(1)), pool.sum(1).mean()))
    decomp(call, pool, 'SHIPPED (emitted q<0.05)')
    decomp(call & S['head'], pool, '  of which head (top-400 |key|)')
    decomp(call & S['block'] & ~S['head'], pool, '  of which CPM block')
    absk = S['absk']
    for K in (400, 1000, 2600):
        M = np.zeros_like(call)
        for i in range(T):
            M[i, np.argsort(-(absk[i] + S['jit'][i]), kind='stable')[:K]] = True
        M &= ~own
        decomp(M, pool, 'key-ranked top-%d (target-specific)' % K)

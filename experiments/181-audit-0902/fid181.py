"""181 — three cheap checks on the call-set arithmetic.

(a) 164's tail-block decomposition (calls 41.8%, tail precision 0.4999, keep 0.5122).
(b) 179 S8c prices hy_kt8000's board `fid` at +0.0029 by carrying h1-test's precision delta
    (+0.0052) straight onto the board.  164 S4 and 171 S7a(iii) instead use the board's OWN
    five-point fid-vs-n_pred series, which has never once risen when the call set shrank.
    What does that series say for n_pred 4,474 -> 3,528?
(c) the h1-test -> board reach shrink 0.428 was calibrated on a8 vs a15.  How well determined
    is it, given h1-test's own reach twin spread?
"""
import numpy as np, sys
sys.path.insert(0,'experiments/153-de-bench')
from debench import Surface
P='experiments/177-resample-emitter/de177.%s.npz'
S=Surface(P%'KA_a15_raw'); z=np.load(P%'KA_a15_raw',allow_pickle=True)
spread=z['spread']
call = S.pred_sig & S.adj
keep = call & (spread==0.25); tail = call & (spread!=0.25)
def prec(m):
    ok=(np.sign(S.pred_lfc)==np.sign(S.real_lfc))
    per=[]
    for p in range(85):
        i=np.flatnonzero(m[p] & S.real_sig[p])
        if i.size: per.append(ok[p,i].mean())
    pooled=ok[m & S.real_sig].mean()
    return pooled, float(np.mean(per)), int((m & S.real_sig).sum())
print('(a) 164 S3.3 decomposition of the call set on the reference-adjudicable pool')
for nm,m in (('keep (spread .25)',keep),('tail (spread 1.0)',tail),('all',call)):
    p,q,n = prec(m); print(f'    {nm:20} calls {n:7d}  pooled {p:.4f}  per-pert {q:.4f}   (164: keep .5122 tail .4999 all .5070)')
print(f'    tail share of the scored call set: {tail[S.real_sig].sum()/call[S.real_sig].sum():.1%}  (164: 41.8%)')

print()
print('(b) the board fid-vs-n_pred series (164 S4), interpolated at hy_kt8000 n_pred 3,528')
series=[(1563,0.49126),(1911,0.49654),(3836,0.50210),(4424,0.50220),(5961,0.50178)]
xs=np.array([s[0] for s in series],float); ys=np.array([s[1] for s in series])
f=float(np.interp(3528,xs,ys)); f0=float(np.interp(4474,xs,ys))
print(f'    board fid at n_pred 4,474 (a23, measured 179) -> {f0:.5f} (interp); at 3,528 -> {f:.5f}')
print(f'    d raw fid {f-f0:+.5f}   d overall {(f-f0)*3.32482/6:+.6f}   (179 S8c claims +0.0029 raw +0.0052)')
print(f'    179 instead carries h1-test precision +0.00516 -> +0.00286 overall.')
print(f'    the two disagree by {abs((0.00516-(f-f0))*3.32482/6):.4f} overall and in SIGN.')
print('    a8 is the direct board test of removing the tail law: n_pred 3,836 fid 0.50210 vs')
print('    a15 n_pred 4,424 fid 0.50220 -> removing the tail block on the BOARD bought -0.0001.')

print()
print('(c) how well determined is the 0.428 h1-test->board reach shrink?')
print('    board a8 - a15 raw reach = +0.00790 (one pair, board has no run-to-run variance)')
print('    h1-test a8flat - a15, twin-averaged = +0.01847; twin spreads 0.00746 (a15) / 0.01068 (a8)')
se = np.hypot(0.00746/2, 0.01068/2)
print(f'    SE of the h1-test twin-mean difference ~ {se:.5f}  -> shrink 0.00790/(0.01847 +- {se:.4f})')
print(f'    = {0.00790/(0.01847+se):.3f} .. {0.00790/(0.01847-se):.3f}   (quoted as 0.428/0.393)')
print()
print('    n_pred cut, board vs h1-test, for the two arms:')
for nm,bh,ba,hh,ha in (('a8 vs a15',3836,4424,3581,6034),('hy_kt8000 vs a23',3528,4474,3304,6066)):
    cb=1-bh/ba; ch=1-hh/ha
    print(f'      {nm:18} board cut {cb:.1%}  h1-test cut {ch:.1%}  ratio {cb/ch:.3f}')
print('    -> the 0.428 shrink was calibrated on a pair whose own board/h1-test cut ratio was 0.33.')
print('       179 multiplies the shrink by ANOTHER 0.46 for the same phenomenon: double counting.')

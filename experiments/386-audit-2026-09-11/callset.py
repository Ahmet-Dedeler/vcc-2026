"""386 audit, row 4: call-set arithmetic from metric definitions (own code)."""
import numpy as np

J = 0.0409          # task-given raw jac (319p2 board raw 0.0409032)
n_pred = 2599       # task-given mean call-set size
fid = 0.508         # task-given raw fid (319p2 board raw 0.5079333)

print('== implied intersection I = J*(n_pred+n_real)/(1+J), recall = I/n_real ==')
for n_real in (100, 288, 765):
    I = J * (n_pred + n_real) / (1 + J)
    feasible = I <= min(n_pred, n_real)
    print('n_real=%4d  I=%7.2f  recall=%.3f  feasible=%s' % (n_real, I, I / n_real, feasible))
print('minimum feasible n_real: J*n_pred = %.1f' % (J * n_pred))

print('\n== fid: k = fid*max(n_pred,n_real); k vs I ==')
for n_real in (100, 288, 765):
    I = J * (n_pred + n_real) / (1 + J)
    k = fid * max(n_pred, n_real)
    print('n_real=%4d  k=%8.1f  I=%7.2f  k/I=%.2f  max in-intersection share=%.3f' % (
        n_real, k, I, k / I, I / k))
    # under uniform accuracy a over call set: a = k/n_pred = fid
    # under coin-flip outside intersection: a = (k - 0.5*(n_pred-I))/I
    a_out = (k - 0.5 * (n_pred - I)) / I
    print('           uniform-a => %.4f (=fid); coinflip-outside => a_in=%.3f' % (k / n_pred, a_out))

print('\n== yield term: call exactly K=400, ceiling = min(1, K/n_real) (k<=K, denom=max(K,n_real)) ==')
for n_real in (288, 1000):
    print('n_real=%4d  ceiling=%.4f' % (n_real, min(1.0, 400 / n_real)))

print('\n== jac break-even at K=400, n_real=288 (need J>=0.0409) ==')
# J = I/(400+288-I) >= 0.0409 -> I >= 0.0409*688/1.0409
I_be = J * 688 / (1 + J)
print('need I>=%.1f i.e. top-400 precision>=%.2f%% (current overall precision at 2599: %.2f%%)' % (
    I_be, 100 * I_be / 400, 100 * (J * (2599 + 288) / (1 + J)) / 2599))

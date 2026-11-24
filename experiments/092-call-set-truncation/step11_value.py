"""Step 11: value-of-an-estimator curve.
Simulate an estimator with lognormal error: nhat_i = n_real_i * exp(sigma * z_i).
sigma=0 is the oracle; sigma ~ 1.017 = sd(log n_real) is a useless estimator (= a constant).
For each sigma, find the best safety multiplier and report the achievable d_avg6."""
import numpy as np, trunc_lib as T
pre = T.precompute(); P = len(pre)
n_real = np.array([s["n_real"] for s in pre], float)
bf, bj, _ = T.fid_jac(pre, [None]*P)
bof, boj = T.scaled(float(np.nanmean(bf)), float(np.nanmean(bj)))
def ev(mv):
    f, j, m = T.fid_jac(pre, mv)
    of, oj = T.scaled(float(np.nanmean(f)), float(np.nanmean(j)))
    return ((of-bof)+(oj-boj))/6.0, m.mean()
rng = np.random.default_rng(7)
CONST = 0.0299
print(f"{'sigma':>7}{'implied R2(log)':>17}{'best safety':>13}{'mean m':>9}{'best d_avg6':>13}{'vs const':>10}{'% of oracle':>13}")
for sg in [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.017]:
    best=(-9,)
    for sf in np.arange(0.6,6.01,0.1):
        ds=[]; ms=[]
        for rep in range(8 if sg>0 else 1):
            z = rng.standard_normal(P) if sg>0 else np.zeros(P)
            d,m = ev(list(sf*n_real*np.exp(sg*z))); ds.append(d); ms.append(m)
        d=float(np.mean(ds))
        if d>best[0]: best=(d,sf,float(np.mean(ms)))
    r2 = max(0.0, 1 - sg**2/1.017**2)
    print(f"{sg:>7.3f}{r2:>17.3f}{best[1]:>13.1f}{best[2]:>9.0f}{best[0]:>+13.4f}{best[0]-CONST:>+10.4f}{100*best[0]/0.0388:>12.0f}%")

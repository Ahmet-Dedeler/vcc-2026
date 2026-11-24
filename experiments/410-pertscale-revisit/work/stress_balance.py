import sys, numpy as np
sys.path.insert(0, 'scripts')
import make_submission_pm_cd4 as M
rng = np.random.default_rng(7)
for t in range(3000):
    G = int(rng.integers(50, 3000))
    a = rng.poisson(rng.gamma(0.2, 5.0, G)).astype(np.int64); b = rng.poisson(rng.gamma(0.2, 5.0, G)).astype(np.int64)
    rho = float(rng.uniform(0, 0.95)); z = np.where(rng.random(G) < 0.5, -1, 1)
    eps = np.floor(rho * np.minimum(a, b)).astype(np.int64) * z
    eps = M._balance_pair_noise(eps)
    assert eps.sum() == 0 and (a + eps).min() >= 0 and (b - eps).min() >= 0
print("stress 3000 random sparse pairs PASS")

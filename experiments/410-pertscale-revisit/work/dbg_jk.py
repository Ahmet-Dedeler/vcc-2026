import numpy as np, sys
sys.path.insert(0, '.')
import run410b as B
rng = np.random.default_rng(0)
# synthetic ballast block resembling H1: K=2, ~995k each, carriers 398 x 20k
G = 18000
x = rng.gamma(0.3, 1.0, G); x /= x.sum()
car = np.vstack([rng.multinomial(20000, x) for _ in range(398)])
Q = rng.multinomial(2 * 995_000, x)
ball = np.vstack([Q // 2 + (Q % 2 > 0) * 0, Q - Q // 2])
block = np.vstack([car, ball]).astype(np.int64)
tgt = B.jk_sum(car) * 3
out, rho, v = B.jk_boost(block, 398, 2, tgt, np.random.default_rng(1))
print("rho", rho, "min", out.min(), "maxrow", out.sum(1).max(), "rows", out[398:].sum(1), "orig", block[398:].sum(1))

"""516: calibrated P(up) from the H1 out-of-fold 'reg' stack -> nmae at confidence-proportional magnitudes (analytic,
H1). Uncalled genes only (calls, called-gene lfc, order untouched). pred_uncalled = sign(P) * (|ours| + kappa |2P-1| m_g)
with m_g the shipped MDE law min(0.1/sqrt(cnt), 0.1) (cnt >= 0.3); kappa scanned in-sample (an upper bound) and by
5-fold target CV (honest)."""
import numpy as np
from common516 import WORK, boot_mean_delta, known_answers, load_h1
from common480 import nmae_rows

d = load_h1(); got, ok, _ = known_answers(d); assert ok; print("KA PASS", got)
P = np.load(WORK / "stackP-h1-reg.npy").astype(float)
unc = ~d["calls"]
cnt = d["ctrl"] * 20190.5 / 1e6
m = np.where(cnt >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(cnt, 1e-3)), 0.1), 0.0)[None, :]
ours = np.nan_to_num(d["ol"], nan=0.0, posinf=0.0, neginf=0.0)
base = nmae_rows(ours, d["gateT"])
s = np.sign(P - 0.5)
def pred(kappa, signswap=True):
    sg = np.where(unc, s, np.sign(ours)) if signswap else np.sign(ours)
    mag = np.abs(ours) + np.where(unc, kappa * np.abs(2 * P - 1) * m, 0.0)
    return np.where(unc, sg * mag, ours)
print(f"production nmae {np.nanmean(base):.5f}")
for kap in (0.0, 0.5, 1.0, 2.0, 3.0, 5.0):
    r = nmae_rows(pred(kap), d["gateT"]); dl = boot_mean_delta(base, r)
    r2 = nmae_rows(pred(kap, False), d["gateT"]); dl2 = boot_mean_delta(base, r2)
    print(f"kappa {kap:3.1f}: stack sign + conf magnitude {np.nanmean(r):.5f} {dl[0]:+.5f} [{dl[1]:+.5f},{dl[2]:+.5f}] | "
          f"production sign + same magnitude {np.nanmean(r2):.5f} {dl2[0]:+.5f}")
# target-CV kappa
T = len(d["targets"]); rng = np.random.default_rng(516); grid = np.array([0, 0.5, 1, 2, 3, 5, 8])
R = np.stack([nmae_rows(pred(k), d["gateT"]) for k in grid])
cv = np.full(T, np.nan)
for rep in range(20):
    perm = rng.permutation(T)
    for te in np.array_split(perm, 5):
        tr = np.setdiff1d(np.arange(T), te)
        kbest = int(np.nanargmin(np.nanmean(R[:, tr], 1)))
        cv[te] = np.where(np.isnan(cv[te]), 0, cv[te]) + R[kbest, te]
cv /= 20
dl = boot_mean_delta(base, cv)
print(f"target-CV kappa: nmae {np.nanmean(cv):.5f} {dl[0]:+.5f} [{dl[1]:+.5f},{dl[2]:+.5f}]  (in-sample best kappa "
      f"{grid[int(np.nanargmin(np.nanmean(R, 1)))]})")

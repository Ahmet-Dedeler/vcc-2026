import numpy as np, scipy.stats as ss, sys
sys.path.insert(0, "experiments/521-strength-nreal")
import lib521 as L
S = L.load_state(0)
rl, rq = S["rl"], S["rq"]; own = S["own"]; genes = np.asarray(S["genes"]).astype(str); targets = np.asarray(S["targets"]).astype(str)
ref_lin = np.asarray(S["ref_lin"], float); w = ref_lin / ref_lin.sum()
T = rl.shape[0]
contrib = np.zeros_like(rl)
m = np.zeros(T); nreal = np.zeros(T); fu = np.zeros(T)
for i in range(T):
    f = np.isfinite(rl[i]) & ~own[i]; de = f & (rq[i] < 0.05); null = f & ~de
    contrib[i, de] = w[de] * (2.0 ** rl[i, de] - 1)
    m[i] = contrib[i].sum(); nreal[i] = de.sum(); fu[i] = np.mean(rl[i, null] > 0)
print("spearman(m_real, n_real) %+.3f; spearman(n_real, frac_up_null) %+.3f" % (ss.spearmanr(m, nreal).correlation, ss.spearmanr(nreal, fu).correlation))
# which gene families carry the mass change
fam = {"RPL/RPS": np.char.startswith(genes, "RPL") | np.char.startswith(genes, "RPS"), "MT-": np.char.startswith(genes, "MT-")}
for k, msk in fam.items():
    part = contrib[:, msk].sum(1)
    print(f"{k:8s} n {msk.sum():4d} share of ctrl {w[msk].sum():.3f} | median share of |m| {np.median(np.abs(part) / np.maximum(np.abs(contrib).sum(1), 1e-12)):.3f} | spearman(part, frac_up) {ss.spearmanr(part, fu).correlation:+.3f}")
top = np.argsort(-np.abs(contrib).sum(0))[:12]
print("top genes by total |mass contribution|:", [(genes[j], round(float(np.abs(contrib[:, j]).sum()), 4)) for j in top])
# per target: frac_up null vs m; show extremes
o = np.argsort(m)
print("most mass LOST:", [(targets[i], round(m[i], 4), int(nreal[i]), round(fu[i], 3)) for i in o[:6]])
print("most mass GAINED:", [(targets[i], round(m[i], 4), int(nreal[i]), round(fu[i], 3)) for i in o[-6:]])

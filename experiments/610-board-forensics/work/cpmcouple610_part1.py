import numpy as np
r = np.load("experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
keep = r["keep"].astype(bool); genes = r["genes"].astype(str)[keep]; targets = r["targets"].astype(str)
lfc = np.asarray(r["lfc"], float); q = np.asarray(r["q"], float)
ref_lin = np.asarray(r["ref_lin"], float); ref_lin = ref_lin[keep] if ref_lin.shape[0] == keep.shape[0] else ref_lin
print("shapes", lfc.shape, q.shape, ref_lin.shape, "ref_lin sum", ref_lin.sum(), "min/med/max", ref_lin.min(), np.median(ref_lin), ref_lin.max())
print("keys", r.files)
w = ref_lin / ref_lin.sum()
own = np.equal.outer(targets, genes)
rows = []
for i, t in enumerate(targets):
    f = np.isfinite(lfc[i]) & ~own[i]
    de = f & (q[i] < 0.05)
    null = f & ~(q[i] < 0.05)
    # net library-share change from real-DE genes (and the target itself if present on the axis)
    dm = float(np.sum(w[de] * (np.power(2.0, lfc[i, de]) - 1.0)))
    tj = np.where(own[i])[0]
    dt = float(w[tj[0]] * (np.power(2.0, lfc[i, tj[0]]) - 1.0)) if len(tj) and np.isfinite(lfc[i, tj[0]]) else np.nan
    hi = null & (ref_lin >= np.quantile(ref_lin, 0.8))
    rows.append((t, int(de.sum()), dm, dt, float(np.mean(lfc[i, null] > 0)), float(np.mean(lfc[i, hi] > 0)), int(null.sum())))
import scipy.stats as ss
R = np.array([(x[2], (x[3] if np.isfinite(x[3]) else 0.0), x[4], x[5], x[1]) for x in rows])
mass = R[:, 0] + R[:, 1]
print("targets", len(rows), "| frac_up(null) mean %.4f sd %.4f range %.3f..%.3f" % (R[:, 2].mean(), R[:, 2].std(), R[:, 2].min(), R[:, 2].max()))
print("frac_up(null, top-20%% expressed) mean %.4f sd %.4f" % (R[:, 3].mean(), R[:, 3].std()))
for lab, x in (("DE mass change", R[:, 0]), ("target share change", R[:, 1]), ("DE+target", mass)):
    s1 = ss.spearmanr(x, R[:, 2]); s2 = ss.spearmanr(x, R[:, 3])
    print(f"spearman({lab}, frac_up null) {s1.correlation:+.3f} (p {s1.pvalue:.2g});  vs top-expressed null {s2.correlation:+.3f} (p {s2.pvalue:.2g})")
# sign rule: predict null genes' sign = -sign(mass change); accuracy over all null genes, per target then mean
acc = []
for (t, nde, dm, dt, fu, fh, nn) in rows:
    m = dm + (dt if np.isfinite(dt) else 0.0)
    acc.append(fu if m < 0 else 1 - fu)
print("oracle-mass sign rule p_null (mean over targets) %.4f; top-expressed %.4f" % (np.mean(acc), np.mean([ (x[5] if (x[2] + (x[3] if np.isfinite(x[3]) else 0)) < 0 else 1 - x[5]) for x in rows])))

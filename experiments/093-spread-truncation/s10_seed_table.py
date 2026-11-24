"""Step 10: seed-repeat table. Mean +/- sd of every metric over the independent seeds,
so the gap between the leading arms can be read against the emitter's own noise."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib093 as L  # noqa: E402
import s5_metrics as S  # noqa: E402

SEEDS_DIR = os.path.join(L.SCRATCH, "seeds093")
base = S.base
bs = {k: L.scaled(k, v) for k, v in zip(("fid", "jac", "reach", "nmae"), base[:4])}

arms = {}
for fn in sorted(os.listdir(SEEDS_DIR)):
    name, _, tag = fn[:-4].partition("__")
    arms.setdefault(name, []).append(os.path.join(SEEDS_DIR, fn))

# fold in the original sweep seed so each arm has four independent emissions
for name in arms:
    p = os.path.join(L.SCRATCH, "sweep093", f"{name}.npz")
    if os.path.exists(p):
        arms[name].append(p)

hdr = (f"{'policy':<16}{'seeds':>6}{'n_pred med':>12}{'fid':>19}{'jac':>19}"
       f"{'reach':>19}{'d_avg6':>19}")
print(hdr)
print("=" * len(hdr))
res = {}
for name, files in sorted(arms.items()):
    f, j, r, n, npm, d6 = [], [], [], [], [], []
    for path in files:
        z = np.load(path)
        a, b, c, e, npred = S.score(z["pred_p"], z["pred_padj"], z["pred_lfc"])
        sc = {k: L.scaled(k, v) for k, v in zip(("fid", "jac", "reach", "nmae"),
                                                (a, b, c, e))}
        f.append(a)
        j.append(b)
        r.append(c)
        n.append(e)
        npm.append(np.median(npred))
        d6.append(sum(sc[k] - bs[k] for k in sc) / 6)
    res[name] = (np.mean(d6), np.std(d6, ddof=1), len(files))
    print(f"{name:<16}{len(files):>6}{np.mean(npm):>12.0f}"
          f"{np.mean(f):>13.6f}+-{np.std(f, ddof=1):<4.4f}"
          f"{np.mean(j):>13.6f}+-{np.std(j, ddof=1):<4.4f}"
          f"{np.mean(r):>13.6f}+-{np.std(r, ddof=1):<4.4f}"
          f"{np.mean(d6):>+13.4f}+-{np.std(d6, ddof=1):<4.4f}")

print("\npairwise gaps against global_s1.0 (mean difference / pooled sd):")
g = res["global_s1.0"]
for name, (m, s, k) in sorted(res.items(), key=lambda x: -x[1][0]):
    if name == "global_s1.0":
        continue
    pooled = np.sqrt((s ** 2 + g[1] ** 2) / 2)
    print(f"  {name:<16} {m - g[0]:+.4f}   ({(m - g[0]) / pooled:+.2f} sd)")

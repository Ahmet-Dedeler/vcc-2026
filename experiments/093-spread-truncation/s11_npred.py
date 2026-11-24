"""Step 11: the achieved n_pred distribution -- the number that decides everything.

`fid = k / max(n_pred, n_real)`. While `n_pred > n_real` the denominator is `n_pred`, so
every gene we call over the reference's count is pure dilution. This prints, per policy,
the full distribution of `n_pred` and the ratio `n_pred / n_real` that the metric actually
divides by.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib093 as L  # noqa: E402
import s5_metrics as S  # noqa: E402

n_real = np.array([int((S.OKS[i] & (S.rq[i] < S.A)).sum()) for i in range(S.P)])
print(f"n_real: median {np.median(n_real):.0f}  mean {n_real.mean():.0f}  "
      f"[{n_real.min()}, {n_real.max()}]\n")

want = ["PRODUCTION", "global_s0.4", "global_s0.55", "global_s0.7", "global_s0.85",
        "global_s1.0", "effK850_s1.0", "nb_phi0.15", "nb_phi0.4"]
hdr = (f"{'policy':<15}{'min':>7}{'p10':>7}{'p25':>7}{'med':>7}{'p75':>7}{'p90':>7}"
       f"{'max':>7}{'mean':>8}{'med ratio':>11}{'over-call':>11}")
print(hdr)
print("=" * len(hdr))
for name in want:
    if name == "PRODUCTION":
        npred = S.base[4]
    else:
        for sub in ("sweep093", "sweep093nb"):
            p = os.path.join(L.SCRATCH, sub, f"{name}.npz")
            if os.path.exists(p):
                break
        z = np.load(p)
        npred = S.score(z["pred_p"], z["pred_padj"], z["pred_lfc"])[4]
    q = np.percentile(npred, [0, 10, 25, 50, 75, 90, 100])
    print(f"{name:<15}" + "".join(f"{v:>7.0f}" for v in q)
          + f"{npred.mean():>8.0f}{np.median(npred / n_real):>11.1f}x"
          + f"{int((npred > n_real).sum()):>8}/{S.P}")

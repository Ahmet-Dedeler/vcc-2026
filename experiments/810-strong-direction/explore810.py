"""810 -- per-target correlates of the strong tier's misdirection (local, light): Spearman of (cos_pre - cos_base), cos_fin,
and the realised cos with donor agreement, source strength, knockdown, noise, key weight; pooled over folds within tier."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import analyse810 as A  # noqa: E402

def main(prefix="L"):
    res, PT = A.main()
    rows = {t: [] for t in (0, 1, 2)}
    for st, (c, pt) in PT.items():
        if not st.startswith(prefix.rstrip("c")) or (prefix == "ABCc" and "kolf2" in st):
            continue                       # prefix "ABCc": the A/B/C cancer folds only (kolf2 out)
        srcs = c["sources"]
        held = st.split("__")[1]
        kd = np.nanmean(np.vstack([c[f"kd_{b}"] for b in srcs]), 0)
        kdmin = np.nanmax(np.vstack([c[f"kd_{b}"] for b in srcs]), 0)       # the weaker knockdown (less negative)
        noise = np.nanmean(np.vstack([c[f"noise_{b}"] for b in srcs]), 0)
        cnt = np.nanmean(np.vstack([np.log1p(c[f"cnt_{b}"]) for b in srcs]), 0)
        held_cnt = np.log1p(c.get(f"cnt_{held}", np.full(len(kd), np.nan)))
        held_noise = c.get(f"noise_{held}", np.full(len(kd), np.nan))
        for t in (0, 1, 2):
            m = (c["tier"] == t) & (pt["rd"] > 0)
            for i in np.flatnonzero(m):
                rows[t].append(dict(dwhite=pt["cos_pre"][i] - pt["cos_base"][i], cos_fin=pt["cos_fin"][i], cos_real=pt["cos_real"][i],
                                    agree=pt["cos_src01"][i], agree_raw=pt["cos_raw01"][i], kd=kd[i], kdweak=kdmin[i], noise=noise[i],
                                    cnt=cnt[i], held_cnt=held_cnt[i], held_noise=held_noise[i], kv=c["kv"][i], z=c["z"][i],
                                    e_h1=pt["energy_h1"][i], e_cd4=pt["energy_cd4"][i], dU=pt["cos_U"][i] - pt["cos_fin"][i],
                                    finU=pt["cos_fin_U"][i]))
    feats = ["agree", "agree_raw", "kd", "kdweak", "noise", "cnt", "held_cnt", "held_noise", "kv", "z", "e_h1", "e_cd4", "finU"]
    for t in (0, 1, 2):
        R = rows[t]
        print(f"tier {t}: n {len(R)}; mean dwhite {np.mean([r['dwhite'] for r in R]):+.4f}, share dwhite>0 {np.mean([r['dwhite'] > 0 for r in R]):.2f}; "
              f"mean dU {np.mean([r['dU'] for r in R]):+.4f}")
        for y in ("dwhite", "dU", "cos_fin", "cos_real"):
            line = []
            for f in feats:
                a = np.asarray([r[f] for r in R], float); b = np.asarray([r[y] for r in R], float)
                ok = np.isfinite(a) & np.isfinite(b)
                rho = spearmanr(a[ok], b[ok]).correlation if ok.sum() > 10 else np.nan
                line.append(f"{f} {rho:+.2f}")
            print(f"   {y:9s}: " + ", ".join(line))

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "L")

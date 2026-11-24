"""Diagnostics on `de_wilcoxon_lfc_nmae`'s gate.  Measurement only, no arms yet.

NMAE_p = sum_{g in S_p} |lfc_hat - lfc| / sum_{g in S_p} |lfc|, S_p chosen by the reference.
Predicting zero everywhere gives exactly 1.  So per gene our prediction is a CREDIT of
|lfc| - |lfc_hat - lfc|, which is negative wherever we do worse than silence.
"""
from __future__ import annotations
import sys
import numpy as np
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/098-nmae")
import nmaekit as K

for W_ in ("h1", "orion"):
    H = K.Harness(W_); eff, s = K.a12_eff(H)
    L = H.lfc_hat(eff[H.rows])
    hs, rs_, ws = [], [], []
    for k, j in enumerate(H.rrows):
        gs = H.sig[j]; hs.append(L[k][gs]); rs_.append(H.rlfc[j][gs])
        ws.append(np.full(gs.size, 1.0 / np.abs(H.rlfc[j][gs]).sum()))
    h = np.concatenate(hs); r = np.concatenate(rs_); w = np.concatenate(ws)
    w = w / w.sum() * h.size                      # per-gene weight in the row-mean NMAE
    nm = float(np.mean([np.abs(a - b).sum() / np.abs(b).sum() for a, b in zip(hs, rs_)]))
    print(f"\n===== {W_}: {H.rows.size} scored rows, {h.size} gate genes, NMAE {nm:.4f} =====")
    print(f"  gate size  median {np.median([x.size for x in rs_]):.0f}  "
          f"mean {np.mean([x.size for x in rs_]):.0f}")
    print(f"  |lfc| ref   mean {np.abs(r).mean():.4f}  median {np.median(np.abs(r)):.4f}")
    print(f"  |lfc| pred  mean {np.abs(h).mean():.4f}  median {np.median(np.abs(h)):.4f}"
          f"   ratio(mean) {np.abs(h).mean()/np.abs(r).mean():.4f}")
    print(f"  sign agreement on gate  {np.mean(np.sign(h) == np.sign(r)):.4f}"
          f"   (weighted {np.average(np.sign(h)==np.sign(r), weights=w):.4f})")
    print(f"  corr(h, r) pearson {np.corrcoef(h, r)[0,1]:.4f}")
    hurt = np.abs(h - r) > np.abs(r)
    print(f"  genes where we are WORSE than silence: {hurt.mean():.4f} of gate "
          f"(weighted {np.average(hurt, weights=w):.4f})")
    cred = (np.abs(r) - np.abs(h - r))
    print(f"  credit from helped genes {np.average(np.where(cred>0,cred,0),weights=w)/np.average(np.abs(r),weights=w):+.4f}"
          f"   damage from hurt genes {np.average(np.where(cred<0,cred,0),weights=w)/np.average(np.abs(r),weights=w):+.4f}")

    # ---- exact L1-optimal global scale on lfc, per row and pooled --------------------
    def nmae_scale(c):
        return float(np.mean([np.abs(c * a - b).sum() / np.abs(b).sum()
                              for a, b in zip(hs, rs_)]))
    grid = np.arange(0.0, 3.01, 0.05)
    vals = np.array([nmae_scale(c) for c in grid])
    b = int(vals.argmin())
    fine = np.arange(max(0, grid[b] - 0.05), grid[b] + 0.051, 0.005)
    fv = np.array([nmae_scale(c) for c in fine]); fb = int(fv.argmin())
    print(f"  --- global scale on lfc_hat ---")
    for c in (0.0, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0):
        print(f"      c={c:<5} NMAE {nmae_scale(c):.4f}")
    print(f"      OPTIMAL c* = {fine[fb]:.3f}  NMAE {fv[fb]:.4f}   (gain {nm-fv[fb]:+.4f})")
    # per-row optimum, as a ceiling for any per-row scale
    per = []
    for a, bb in zip(hs, rs_):
        g2 = np.arange(0.0, 3.01, 0.02)
        per.append(min(np.abs(c * a - bb).sum() / np.abs(bb).sum() for c in g2))
    print(f"      per-ROW optimal scale ceiling  NMAE {np.mean(per):.4f} (gain {nm-np.mean(per):+.4f})")
    # sign-only ceiling: predict k*sign(h)
    print(f"  --- what the direction alone is worth ---")
    for kk in (0.25, 0.5, 1.0):
        v = float(np.mean([np.abs(kk * np.sign(a) * np.abs(b).mean() - b).sum() / np.abs(b).sum()
                           for a, b in zip(hs, rs_)]))
        print(f"      predict {kk}*mean|lfc|*sign(ours): NMAE {v:.4f}")
    v = float(np.mean([np.abs(np.sign(b) * np.abs(b).mean() - b).sum() / np.abs(b).sum()
                       for a, b in zip(hs, rs_)]))
    print(f"      ORACLE sign at mean magnitude:      NMAE {v:.4f}")
    v = float(np.mean([np.abs(np.abs(b) * np.sign(a) - b).sum() / np.abs(b).sum()
                       for a, b in zip(hs, rs_)]))
    print(f"      our sign at ORACLE magnitude:       NMAE {v:.4f}")

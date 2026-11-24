"""The per-perturbation channel -- the only one that is EXACTLY free on `pds`.

`pds_cosine` normalises each row before the cosine, so scaling row p's effect by t_p leaves it
bit-identical (097 s9e, verified there to four decimals).  Under the linearisation

    NMAE = mean_p [ 1 - <sign(r_p), h_p> / sum_{S_p}|r_p| ] = mean_p [1 - c_p],

a per-row scale t_p changes it by -mean_p (t_p - 1) c_p.  So the question is only whether c_p
is predictable from something buildable.  Reported: the distribution of c_p, the ORACLE ceiling
(t_p from the sign of c_p), and every buildable predictor's rank correlation with c_p.
"""
from __future__ import annotations
import sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/098-nmae"); sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
import nmaekit as K
from vcc2026.donor_transfer import Donor
from scipy.stats import rankdata, spearmanr

DONS = {"h1": ["k562_gwps", "orion_hek293t", "orion_hct116"],
        "orion": ["k562_gwps", "orion_hek293t"]}
DEPTHS = (0.05, 0.35, 0.60)

for W_ in ("h1", "orion"):
    H = K.Harness(W_); eff, _ = K.a12_eff(H)
    z = np.load(f"{ROOT}/experiments/098-nmae/work/refde_{W_}.npz", allow_pickle=True)
    QQ, RL = z["q"], z["lfc"].astype(np.float64)
    L = H.lfc_hat(eff[H.rows])
    # per-perturbation donor cosine, exactly production's key
    M, HV = {}, {}
    for nd in DONS[W_]:
        dn = Donor.load(f"{ROOT}/data/donors/{nd}"); m, s_ = dn.decompose()
        gi = {g: i for i, g in enumerate(dn.genes)}
        cc = np.array([gi.get(str(g), -1) for g in H.genes]); ok = cc >= 0
        pi = {p: i for i, p in enumerate(dn.perts)}
        Mx = np.zeros((H.P, H.G)); hv = np.zeros(H.P, bool)
        for r, x in enumerate(H.tg):
            if x in pi: Mx[r, ok] = s_[pi[x]][cc[ok]]; hv[r] = True
        M[nd] = Mx; HV[nd] = hv; del dn, s_
    dc = np.zeros(H.P); nd_ = np.zeros(H.P, int)
    for i in range(H.P):
        have = [n for n in DONS[W_] if HV[n][i]]; nd_[i] = len(have)
        tot, cnt = 0.0, 0
        for a in range(len(have)):
            for b in range(a + 1, len(have)):
                x, y = M[have[a]][i], M[have[b]][i]
                tot += float(x @ y) / max(np.linalg.norm(x) * np.linalg.norm(y), 1e-12); cnt += 1
        dc[i] = tot / cnt if cnt else 0.0
    print(f"\n===== {W_} =====")
    for al in DEPTHS:
        rows, cp, den, hn, dn2, dcv = [], [], [], [], [], []
        for k, (i, j) in enumerate(zip(H.rows, H.rrows)):
            s = (QQ[j] < al) & np.isfinite(RL[j]) & (RL[j] != 0)
            gs = np.flatnonzero(s)
            if gs.size < 10: continue
            r = RL[j][gs]; h = L[k][gs]
            d0 = np.abs(r).sum()
            nm = np.abs(h - r).sum() / d0
            rows.append(k); cp.append(1.0 - nm); den.append(d0 / gs.size)
            hn.append(np.abs(h).mean()); dn2.append(float(np.linalg.norm(eff[i])))
            dcv.append(dc[i])
        cp = np.array(cp); den = np.array(den); hn = np.array(hn)
        dn2 = np.array(dn2); dcv = np.array(dcv)
        nm0 = 1.0 - cp.mean()
        # oracle ceiling for a per-row scale in [0, 2]
        orc = 1.0 - np.where(cp > 0, 2 * cp, 0.0).mean()
        print(f"  BH {al:<5} rows {cp.size:>4}  NMAE {nm0:.4f}   rows with c_p<0 "
              f"{np.mean(cp<0):.3f}   c_p sd {cp.std():.4f}")
        print(f"      oracle per-row scale in [0,2]: NMAE {orc:.4f} ({orc-nm0:+.4f})")
        for nmk, v in (("donor cosine", dcv), ("our ||dh_p||", dn2),
                       ("our mean|lfc| on gate", hn), ("ref mean|lfc| (ORACLE)", den)):
            rho = spearmanr(v, cp).statistic
            # buildable arm: t_p = 1 + a*(rank01(v) - 0.5), a = 1
            t = 1.0 + 1.0 * ((rankdata(v) - 0.5) / v.size - 0.5)
            print(f"      {nmk:<24} rho(.,c_p) {rho:+.4f}   linear-rank tilt a=1 -> "
                  f"NMAE {1.0-(t*cp).mean():.4f} ({1.0-(t*cp).mean()-nm0:+.4f})")

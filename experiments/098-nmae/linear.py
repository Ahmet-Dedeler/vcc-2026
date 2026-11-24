"""In OUR regime the metric is linear, and that changes what it wants.

|lfc_hat| is 0.14-0.17x |lfc| on the gate, so for almost every gene |h| <= |r| and

    |h - r| = |r| - sign(r) * h     ==>     NMAE_p = 1 - <sign(r), h> / sum|r|.

So `de_wilcoxon_lfc_nmae` is, here, a DIRECTION member with |h| as its own weight and the
reference's |lfc| appearing only in the normaliser.  This script (1) verifies the linearisation,
(2) reports the direction balance of the reference's gate, (3) measures how well the reference's
gate MEMBERSHIP is predictable, and (4) ranks per-gene confidence keys by the weighted sign
accuracy they buy -- the only buildable handle on a linear objective.
"""
from __future__ import annotations
import sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/098-nmae"); sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
import nmaekit as K
from vcc2026.donor_transfer import Donor
from scipy.stats import rankdata

DONS = {"h1": ["k562_gwps", "orion_hek293t", "orion_hct116"],
        "orion": ["k562_gwps", "orion_hek293t"]}


def donor_keys(H, which):
    """production's all-donor unanimity agreement, per (row, gene) -- exp097's `ag`."""
    P, G = H.P, H.G
    M, HV = {}, {}
    for nd in DONS[which]:
        dn = Donor.load(f"{ROOT}/data/donors/{nd}"); m, s_ = dn.decompose()
        gi = {g: i for i, g in enumerate(dn.genes)}
        cc = np.array([gi.get(str(g), -1) for g in H.genes]); ok = cc >= 0
        pi = {p: i for i, p in enumerate(dn.perts)}
        Mx = np.zeros((P, G)); hv = np.zeros(P, bool)
        for r, x in enumerate(H.tg):
            if x in pi: Mx[r, ok] = s_[pi[x]][cc[ok]]; hv[r] = True
        M[nd] = Mx / (float(np.median(np.linalg.norm(s_[:, cc[ok]], axis=1))) or 1.0)
        HV[nd] = hv; del dn, s_
    ag = np.zeros((P, G)); un = np.zeros((P, G))
    for i in range(P):
        have = [n for n in DONS[which] if HV[n][i]]
        if len(have) < 2: continue
        st = np.vstack([M[n][i] for n in have]); sg = np.sign(st)
        u = (np.abs(sg.sum(0)) == len(have)) & np.all(st != 0, axis=0)
        ag[i] = np.where(u, np.exp(np.log(np.abs(st) + 1e-300).mean(0)), 0.0)
        un[i] = u.astype(float)
    return ag, un


def _wa(v, w):
    s = w.sum()
    return float(np.average(v, weights=w)) if s > 0 else float(np.mean(v))


for W_ in ("h1", "orion"):
    H = K.Harness(W_); eff, _ = K.a12_eff(H)
    L = H.lfc_hat(eff[H.rows])
    ag, un = donor_keys(H, W_)
    g = H.gate
    hs, rs_, ws, ags, uns, exs, ridx = [], [], [], [], [], [], []
    for k, (i, j) in enumerate(zip(H.rows, H.rrows)):
        gs = H.sig[j]
        hs.append(L[k][gs]); rs_.append(H.rlfc[j][gs])
        ws.append(np.full(gs.size, 1.0 / np.abs(H.rlfc[j][gs]).sum() / H.rows.size))
        ags.append(ag[i][g][gs]); uns.append(un[i][g][gs]); exs.append(H.ref_lin[gs])
        ridx.append(np.full(gs.size, k))
    h = np.concatenate(hs); r = np.concatenate(rs_); w = np.concatenate(ws)
    agv = np.concatenate(ags); unv = np.concatenate(uns); ex = np.concatenate(exs)
    nm = float(np.mean([np.abs(a - b).sum() / np.abs(b).sum() for a, b in zip(hs, rs_)]))
    lin = 1.0 - float((np.sign(r) * h * w).sum())
    print(f"\n===== {W_}: NMAE {nm:.5f}   linearisation 1 - <sign(r),h>_w = {lin:.5f}"
          f"   residual {nm-lin:+.5f} =====")
    print(f"  genes with |h| <= |r|: {np.average(np.abs(h)<=np.abs(r), weights=w):.4f}")
    # (2) direction balance
    print(f"  reference gate direction: {np.average(r<0, weights=w*np.abs(r)):.4f} DOWN "
          f"(unweighted {np.mean(r<0):.4f});  our prediction on the gate: "
          f"{np.average(h<0, weights=w*np.abs(h)):.4f} DOWN (unweighted {np.mean(h<0):.4f})")
    # the value of a pure directional bias: h -> h + delta on the gate
    print("  --- a constant lfc shift applied to EVERY gated gene (the free directional bias) ---")
    for dl in (-0.05, -0.02, -0.01, 0.0, 0.01, 0.02, 0.05):
        v = float(np.mean([np.abs(a + dl - b).sum() / np.abs(b).sum() for a, b in zip(hs, rs_)]))
        print(f"      delta={dl:+.3f}  NMAE {v:.5f}  ({v-nm:+.5f})")
    # (3) is gate MEMBERSHIP predictable?
    allg = np.arange(g.sum())
    inn, out = [], []
    for k, (i, j) in enumerate(zip(H.rows, H.rrows)):
        m = np.zeros(g.sum(), bool); m[H.sig[j]] = True
        inn.append(H.ref_lin[m]); out.append(H.ref_lin[~m])
    inn = np.concatenate(inn); out = np.concatenate(out)
    print(f"  --- gate membership vs control expression --- "
          f"median CPM in-gate {np.median(inn):.1f}  out-of-gate {np.median(out):.1f}")
    # (4) confidence keys, ranked by the weighted sign accuracy of the top quintile
    print(f"  --- per-gene confidence keys on the gate (weight = |h|/sum|r|, the metric's own) ---")
    print(f"      {'key':<22}{'rho(key, sign ok)':>18}{'top-20% acc':>13}{'bot-20% acc':>13}"
          f"{'credit share top20':>20}")
    ok = (np.sign(h) == np.sign(r)).astype(float)
    wh = w * np.abs(h)
    cred = w * np.sign(r) * h
    keys = {"|lfc_hat| (self)": np.abs(h),
            "donor agreement": agv,
            "donor unanimity": unv,
            "control CPM": ex,
            "|lfc_hat| x agreement": np.abs(h) * agv,
            "-control CPM": -ex}
    for nmk, kv in keys.items():
        rr = rankdata(kv) / kv.size
        hi = rr > 0.8; lo = rr < 0.2
        rho = float(np.corrcoef(rankdata(kv), rankdata(ok))[0, 1])
        print(f"      {nmk:<22}{rho:>18.4f}"
              f"{_wa(ok[hi], wh[hi]):>13.4f}"
              f"{_wa(ok[lo], wh[lo]):>13.4f}"
              f"{cred[hi].sum()/max(cred.sum(),1e-30):>20.4f}")

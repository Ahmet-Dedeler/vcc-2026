"""Where does the fold-change error live?  Stratify the gate by our own |lfc_hat| and by
control expression, and ask (i) is our sign better in some strata, (ii) what global scale each
stratum wants, (iii) does a stratum-wise scale fitted on one harness transfer to the other.

Rule 10: a per-gene SLOPE is a cell-line property (transfers at 1.02x); a per-gene CONFIDENCE
is a model property (transfers at 1.17x).  A stratification on our OWN predicted magnitude or
on the control's expression is a model/experiment property, not a cell-line one -- which is the
shape that has transferred before.  The transfer test below is what decides it.
"""
from __future__ import annotations
import sys
import numpy as np
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/098-nmae")
import nmaekit as K

ND = 10
D = {}
for W_ in ("h1", "orion"):
    H = K.Harness(W_); eff, s = K.a12_eff(H)
    L = H.lfc_hat(eff[H.rows])
    cpm = H.ref_lin
    hs, rs_, ds_, es_ = [], [], [], []
    for k, j in enumerate(H.rrows):
        gs = H.sig[j]
        hh = L[k][gs]; rr = H.rlfc[j][gs]
        hs.append(hh); rs_.append(rr); es_.append(cpm[gs])
        # within-row decile of |lfc_hat|
        o = np.argsort(np.abs(hh)); dq = np.empty(gs.size, int)
        dq[o] = np.minimum((np.arange(gs.size) * ND) // max(gs.size, 1), ND - 1)
        ds_.append(dq)
    h = np.concatenate(hs); r = np.concatenate(rs_); dq = np.concatenate(ds_)
    ex = np.concatenate(es_)
    # weight each gene by 1/sum|r| of its row, so sums reproduce the row-mean NMAE
    wt = np.concatenate([np.full(x.size, 1.0 / np.abs(x).sum()) for x in rs_]) / len(rs_)
    D[W_] = dict(h=h, r=r, dq=dq, ex=ex, wt=wt, hs=hs, rs=rs_, ds=ds_)
    base = float((np.abs(h - r) * wt).sum())
    print(f"\n===== {W_}  NMAE {base:.4f}  ({len(rs_)} rows, {h.size} gate genes) =====")
    print(f"  {'|lfc_hat| decile':<18}{'n':>7}{'wt share':>10}{'sign acc':>10}"
          f"{'mean|h|':>9}{'mean|r|':>9}{'c* L1':>8}{'ΔNMAE at c*':>12}")
    for d in range(ND):
        m = dq == d
        cs = np.arange(0.0, 6.01, 0.02)
        f = np.array([((np.abs(c * h[m] - r[m]) - np.abs(h[m] - r[m])) * wt[m]).sum() for c in cs])
        b = int(f.argmin())
        print(f"  {d:<18}{m.sum():>7}{(wt[m]*np.abs(r[m])).sum():>10.4f}"
              f"{np.average(np.sign(h[m])==np.sign(r[m]), weights=wt[m]*np.abs(r[m])):>10.4f}"
              f"{np.abs(h[m]).mean():>9.4f}{np.abs(r[m]).mean():>9.4f}"
              f"{cs[b]:>8.2f}{f[b]:>12.5f}")
    # expression strata
    q = np.quantile(ex, np.linspace(0, 1, ND + 1)); q[0] -= 1
    eq = np.clip(np.searchsorted(q, ex, side="left") - 1, 0, ND - 1)
    D[W_]["eq"] = eq
    print(f"  {'ctrl-CPM decile':<18}{'n':>7}{'wt share':>10}{'sign acc':>10}"
          f"{'mean|h|':>9}{'mean|r|':>9}{'c* L1':>8}{'ΔNMAE at c*':>12}  cpm range")
    for d in range(ND):
        m = eq == d
        cs = np.arange(0.0, 6.01, 0.02)
        f = np.array([((np.abs(c * h[m] - r[m]) - np.abs(h[m] - r[m])) * wt[m]).sum() for c in cs])
        b = int(f.argmin())
        print(f"  {d:<18}{m.sum():>7}{(wt[m]*np.abs(r[m])).sum():>10.4f}"
              f"{np.average(np.sign(h[m])==np.sign(r[m]), weights=wt[m]*np.abs(r[m])):>10.4f}"
              f"{np.abs(h[m]).mean():>9.4f}{np.abs(r[m]).mean():>9.4f}"
              f"{cs[b]:>8.2f}{f[b]:>12.5f}   {ex[m].min():.0f}-{ex[m].max():.0f}")

# ---- transfer: fit per-stratum scale on one harness, apply on the other -----------------
def fit(dd, key):
    out = np.ones(ND)
    for d in range(ND):
        m = dd[key] == d
        cs = np.arange(0.0, 6.01, 0.01)
        f = np.array([(np.abs(c * dd["h"][m] - dd["r"][m]) * dd["wt"][m]).sum() for c in cs])
        out[d] = cs[int(f.argmin())]
    return out

def apply(dd, key, v):
    c = v[dd[key]]
    return float((np.abs(c * dd["h"] - dd["r"]) * dd["wt"]).sum())

print("\n===== TRANSFER: per-stratum L1 scale fitted on one harness, applied to the other =====")
for key, nm in (("dq", "|lfc_hat| decile"), ("eq", "ctrl-CPM decile")):
    vh, vo = fit(D["h1"], key), fit(D["orion"], key)
    for tgt, src, vsrc, vown in (("orion", "h1", vh, vo), ("h1", "orion", vo, vh)):
        dd = D[tgt]
        b = float((np.abs(dd["h"] - dd["r"]) * dd["wt"]).sum())
        g1 = float((np.abs(np.median(vsrc) * dd["h"] - dd["r"]) * dd["wt"]).sum())
        print(f"  {nm:<18} on {tgt:<6} base {b:.4f}  own-fit {apply(dd,key,vown):.4f}"
              f"  FROM-{src} {apply(dd,key,vsrc):.4f}   (flat scale from {src}: {g1:.4f})")
    print(f"     fitted scales  h1 {np.round(vh,2)}\n                    orion {np.round(vo,2)}")

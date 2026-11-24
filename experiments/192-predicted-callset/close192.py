"""192 — closing measurements: the v14 COMPOSITION test, the leakage re-fit, and the priced
optimum at the retrodicted board level.

Why the composition test exists.  `retro192.py` retrodicts the board's five-point `fid` curve
to 3% with one parameter but gets v14's `jac` SIGN WRONG.  v14 moved two variables -- the call
set's SIZE (5,961 -> 1,563) and, through `--tail-spread 1.0 --keep-top 850`, its COMPOSITION
relative to our own `|predicted effect|` ordering.  My size-only model cannot see the second.
Phase 1 says our `|lfc_pred|` head is DEPLETED of the reference's DE genes (fold 1.14 at n=25
falling to 0.72 at n=300), so a call set pushed OFF that head should gain `jac`.  This tests
that directly: three compositions at v14's own realised n_pred, one board number to match.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
E185 = os.path.join(REPO, "experiments/185-rank-controlled-emission")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, HERE)
import pandas as pd                                                       # noqa: E402
from debench import Surface                                               # noqa: E402
from frontier192 import A15, INV, PowerSurface, order_key, topk_call      # noqa: E402
from keys192 import gene_stats, h_to_lfc                                  # noqa: E402

BEST_F = 0.5
SHRINK = 0.4277
V14_BOARD = dict(fid=-0.035470, jac=+0.004271, reach=+0.013731)


def masked_topk(key, K, banned):
    k = np.where(banned, -np.inf, key)
    return topk_call(k, np.full(key.shape[0], K))


def main() -> None:
    E = np.load(os.path.join(E169, "work/h1test-effect.npz"), allow_pickle=True)
    pos = {str(g): i for i, g in enumerate(E["genes"])}
    tpos = {str(t): i for i, t in enumerate(E["targets"])}
    s = Surface(A15)
    z15 = np.load(A15, allow_pickle=True)
    gi = np.array([pos[g] for g in s.genes])
    ti = np.array([tpos[t] for t in s.perts])
    cpm = E["cpm"][gi].astype(np.float64)
    Hm = E["H"][np.ix_(ti, gi)].astype(np.float64)
    lfcP = np.abs(h_to_lfc(Hm, cpm))
    absH = np.abs(Hm)
    T = absH.shape[0]
    gs = gene_stats(os.path.join(E185, "work/ctrl_cs.npy"), os.path.join(WORK, "genestats.npz"))
    sig = s.real_sig & np.isfinite(s.real_lfc)
    flat = -np.log(np.maximum(
        np.load(os.path.join(E185, "work/de185.armC_a8_flat.npz"),
                allow_pickle=True)["pred_p"], 1e-300))       # the flat-0.25 law's own p-order
    inc = -np.log(np.maximum(z15["pred_p"], 1e-300))
    sf = PowerSurface(A15, BEST_F)
    res = {}

    # ---------------------------------------------------------------- v14 composition test
    print("=" * 100)
    print(f"v14 COMPOSITION TEST at the retrodicted level f = {BEST_F}.  v11 n_pred 5,961 -> "
          f"v14 n_pred 1,563.")
    base_call = topk_call(flat, np.full(T, 5961))
    b = dict(fid=sf.fid(call=base_call), jac=sf.jac(call=base_call),
             reach=sf.reach(key=order_key(base_call, absH)))
    print(f"  v11 model (top-5,961 by the FLAT-law p-order): fid {b['fid']:.5f} "
          f"jac {b['jac']:.5f} reach {b['reach']:.5f}")
    head850 = topk_call(lfcP, np.full(T, 850))
    comps = {
        "C size-only  (top-1,563 by flat p-order)":
            topk_call(flat, np.full(T, 1563)),
        "A head-first (top-850 by |lfc_pred| + 713 by flat p-order)":
            head850 | masked_topk(flat, 713, head850),
        "B tail-only  (top-1,563 by flat p-order OFF the |lfc_pred| top-850)":
            masked_topk(flat, 1563, head850),
    }
    print(f"  {'composition':<62}{'fid':>9}{'jac':>9}{'reach':>8}"
          f"{'d fid sc':>10}{'d jac sc':>10}{'d reach sc':>11}")
    for lab, call in comps.items():
        v = dict(fid=sf.fid(call=call), jac=sf.jac(call=call),
                 reach=sf.reach(key=order_key(call, absH)))
        d = {m: (v[m] - b[m]) * INV[m] * (SHRINK if m == "reach" else 1.0) for m in v}
        res[lab] = dict(raw=v, scaled=d)
        print(f"  {lab:<62}{v['fid']:>9.5f}{v['jac']:>9.5f}{v['reach']:>8.4f}"
              f"{d['fid']:>+10.5f}{d['jac']:>+10.5f}{d['reach']:>+11.5f}")
    print(f"  {'BOARD v14 - v11':<62}{'':>9}{'':>9}{'':>8}"
          f"{V14_BOARD['fid']:>+10.5f}{V14_BOARD['jac']:>+10.5f}{V14_BOARD['reach']:>+11.5f}")

    # ---------------------------------------------------------------- leakage audit / re-fit
    print(f"\n{'=' * 100}")
    panel = set(pd.read_csv(os.path.join(REPO, "data/challenge/pert_counts.csv"))
                ["target_gene"].astype(str))
    hit = sorted(set(s.perts.tolist()) & panel)
    print(f"FALSIFIER 4 — LEAKAGE AUDIT.  {len(hit)}/85 h1-test TARGETS are in the challenge's "
          f"300-gene panel: {hit}")
    print(f"                              {len(set(s.genes.tolist()) & panel)}/{len(s.genes)} "
          f"h1-test GENES are panel targets (they are features, never labels).")
    from sklearn.linear_model import LogisticRegression
    X = np.stack([np.broadcast_to(np.log(cpm + 1)[None, :], absH.shape),
                  np.broadcast_to(np.log(gs["nz"][:] + 1e-3)[None, :], absH.shape),
                  np.broadcast_to(np.log(gs["cv2"][:] + 1e-3)[None, :], absH.shape)], -1)
    Xz = (X - X.reshape(-1, 3).mean(0)) / X.reshape(-1, 3).std(0)
    banned_rows = np.array([p in panel for p in s.perts])
    fits = {}
    for lab, drop in (("LOTO (all 85 targets in the pool)", np.zeros(T, bool)),
                      ("LOTO with the 5 panel targets DELETED from every fit", banned_rows)):
        rngf = np.random.default_rng(7)
        out = np.zeros(absH.shape)
        for p in range(T):
            tr = np.array([q for q in range(T) if q != p and not drop[q]])
            Xt = Xz[tr].reshape(-1, 3)
            yt = sig[tr].reshape(-1)
            idx = rngf.choice(Xt.shape[0], size=min(300_000, Xt.shape[0]), replace=False)
            m = LogisticRegression(max_iter=400).fit(Xt[idx], yt[idx])
            out[p] = m.decision_function(Xz[p])
        fits[lab] = out

    # ---------------------------------------------------------------- the priced optimum
    print(f"\n{'=' * 100}")
    print(f"THE PRICED OPTIMUM at f = {BEST_F}, against a15's CHALLENGE operating point "
          f"(n_pred 4,474).")
    a15_call = topk_call(inc, np.full(T, 4474))
    a0 = dict(fid=sf.fid(call=a15_call), jac=sf.jac(call=a15_call),
              reach=sf.reach(key=order_key(a15_call, absH)))
    print(f"  a15 challenge-matched baseline: fid {a0['fid']:.5f} jac {a0['jac']:.5f} "
          f"reach {a0['reach']:.5f}")
    Ks = [800, 1200, 1600, 2200, 2800, 3600, 4474, 6000]
    for lab, key in list(fits.items()) + [("detectability (cpm) only",
                                           np.broadcast_to(cpm[None, :], absH.shape).copy())]:
        print(f"\n  key = {lab}")
        print(f"    {'K':>7}{'fid':>9}{'jac':>9}{'reach':>8}{'d fid':>10}{'d jac':>10}"
              f"{'d reach':>10}{'SUM(3) overall':>16}")
        for K in Ks:
            call = topk_call(key, np.full(T, K))
            v = dict(fid=sf.fid(call=call), jac=sf.jac(call=call),
                     reach=sf.reach(key=order_key(call, absH)))
            d = {m: (v[m] - a0[m]) * INV[m] * (SHRINK if m == "reach" else 1.0) / 6 for m in v}
            res[f"{lab}|{K}"] = dict(raw=v, overall=d)
            print(f"    {K:>7}{v['fid']:>9.5f}{v['jac']:>9.5f}{v['reach']:>8.4f}"
                  f"{d['fid']:>+10.5f}{d['jac']:>+10.5f}{d['reach']:>+10.5f}"
                  f"{sum(d.values()):>+16.5f}")

    with open(os.path.join(WORK, "close192.json"), "w") as fh:
        json.dump(dict(baseline=a0, v11=b, rows=res, panel_targets=hit), fh, indent=1,
                  default=float)
    print(f"\nwrote {WORK}/close192.json")


if __name__ == "__main__":
    main()

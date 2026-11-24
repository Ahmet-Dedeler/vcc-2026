"""171 — score every realized emission arm on all four DE members + pds, and show WHERE in
the walk each arm's reach comes from.

`mse` is not computable here (it needs the challenge panel's jackknife) and is priced by
`claim171.py`; every other member is exact given the cached DE table.

Reads `work/de171.*.npz` written by `realize171.py`.

KNOWN-ANSWER GATES: `a15` must read reach 0.18283 / median n_pred 6035 and `armC_a8_flat`
0.20291 / 3581 (157 `power_h1test.log`, 164 section 1).  Both are fresh emissions here, so
agreement is a real reproduction and not a cache hit.
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
sys.path.insert(0, f"{ROOT}/experiments/171-reach-prefix")
from debench import Surface, overall_delta  # noqa: E402
from diag171 import Walk, kstar, spec_order  # noqa: E402

W_DIR = f"{ROOT}/experiments/171-reach-prefix/work"
ORDER = ["a15", "armC_a8_flat", "armC_a8_flat_b", "ks000", "ks000_b", "ks010",
         "kt2000_ks000", "kt2000_ks025", "flat100_phi4", "kt10000", "kt20000", "kt20000_b"]
FACTORS = [1.0, 2.0, 3.5]
ALPHA = 0.05


def pool_at_factor(path, real_sig, self_idx, f):
    if f == 1.0:
        return real_sig.copy()
    z = np.load(path, allow_pickle=True)
    rp = np.where(np.isfinite(z["real_p"]), z["real_p"], np.inf)
    sig = np.zeros(rp.shape, bool)
    for i in range(rp.shape[0]):
        k = int(round(f * int(real_sig[i].sum())))
        k = max(0, min(k, int(np.isfinite(z["real_p"][i]).sum())))
        if k:
            sig[i, np.argsort(rp[i], kind="stable")[:k]] = True
    sig &= np.isfinite(z["real_p"])
    for p, j in enumerate(self_idx):
        if j >= 0:
            sig[p, j] = False
    return sig


def reach_pool(W, real_sig):
    out = []
    for p in range(W.perts.size):
        pool = np.flatnonzero(real_sig[p])
        if pool.size == 0:
            continue
        o = spec_order(pool, W.genes, W.pred_sig[p], W.pred_padj[p], W.pred_p[p],
                       W.pred_lfc[p])
        walked = pool[o]
        walked = walked[W.adj[p, walked]]
        if walked.size == 0:
            out.append(0.0)
            continue
        ok = np.sign(W.pred_lfc[p, walked]) == np.sign(W.real_lfc[p, walked])
        out.append(kstar(ok)[0] / pool.size)
    return float(np.mean(out))


def main():
    paths = {os.path.basename(p).split(".")[1]: p
             for p in glob.glob(f"{W_DIR}/de171.*.npz")}
    names = [n for n in ORDER if n in paths] + [n for n in paths if n not in ORDER]
    base = None

    print("=== 171 — realized emission arms on h1-test (85 perturbations x 10,779 genes)")
    print()
    print("KNOWN-ANSWER GATES")
    for nm, want_r, want_n in (("a15", 0.18283, 6035), ("armC_a8_flat", 0.20291, 3581)):
        if nm not in paths:
            continue
        S = Surface(paths[nm])
        npd = np.median([int((S.pred_sig[p]).sum()) for p in range(S.perts.size)])
        r = S.reach()
        print(f"  {nm:<14} reach {r:.5f} (expect {want_r:.5f}, "
              f"{'PASS' if abs(r-want_r) < 0.002 else 'FAIL'})   "
              f"median n_pred {npd:.0f} (expect {want_n}, "
              f"{'PASS' if abs(npd-want_n) < 40 else 'FAIL'})")
    print()

    print("SIX MEMBERS (mse from claim171.py; pds/nmae read pred_lfc, which the exact-total")
    print("repair pins to within re-seed noise -- 164 section 3.2 measured the re-seed twin as")
    print("MOVING MORE than any arm, so 'invariant' here means 'inside re-seed noise')")
    hdr = (f"  {'arm':<15}{'n_pred':>8}{'reach f=1':>11}{'f=2':>9}{'f=3.5':>9}"
           f"{'fid':>9}{'precis':>9}{'jac':>9}{'nmae':>9}{'pds':>9}"
           f"{'d(fid+reach+jac)':>18}")
    print(hdr)
    rows = {}
    for nm in names:
        S = Surface(paths[nm])
        W = Walk(paths[nm])
        rs = [reach_pool(W, pool_at_factor(paths[nm], W.real_sig, W.self_idx, f))
              for f in FACTORS]
        npd = np.median([int(S.pred_sig[p].sum()) for p in range(S.perts.size)])
        # 164: the challenge-readable form of fid is precision, since n_pred >> n_real there
        prec = []
        for p in range(S.perts.size):
            sel = np.flatnonzero(S.pred_sig[p] & S.adj[p])
            if sel.size:
                prec.append((np.sign(S.pred_lfc[p, sel])
                             == np.sign(S.real_lfc[p, sel])).mean())
        d = dict(n=npd, reach=rs, fid=S.fid(), prec=float(np.mean(prec)), jac=S.jac(),
                 nmae=S.nmae(), pds=S.pds())
        rows[nm] = d
        if base is None:
            base = d
        dov = overall_delta(d_fid=d["fid"] - base["fid"], d_reach=rs[0] - base["reach"][0],
                            d_jac=d["jac"] - base["jac"])
        print(f"  {nm:<15}{npd:>8.0f}{rs[0]:>11.5f}{rs[1]:>9.5f}{rs[2]:>9.5f}"
              f"{d['fid']:>9.5f}{d['prec']:>9.5f}{d['jac']:>9.5f}{d['nmae']:>9.5f}"
              f"{d['pds']:>9.5f}{dov:>+18.5f}")
    print()

    print("WALK ANATOMY PER ARM (f = 1)")
    print(f"  {'arm':<15}{'k* med':>8}{'k* mean':>9}{'k*=0 %':>8}{'P(1)':>8}{'P(3)':>8}"
          f"{'P(10)':>8}{'% head in TAIL block':>22}")
    for nm in names:
        W = Walk(paths[nm])
        rws = W.rows()
        ks = np.array([r["kstar"] for r in rws])
        p1 = np.mean([r["purity"][0] for r in rws if r["n_adj"] >= 1])
        p3 = np.mean([r["purity"][2] for r in rws if r["n_adj"] >= 3])
        p10 = np.mean([r["purity"][9] for r in rws if r["n_adj"] >= 10])
        t = n = 0
        for r in rws:
            m = min(10, r["n_adj"])
            n += m
            t += int(r["tail_walked"][:m].sum())
        print(f"  {nm:<15}{np.median(ks):>8.0f}{ks.mean():>9.2f}"
              f"{100*np.mean(ks == 0):>8.1f}{p1:>8.4f}{p3:>8.4f}{p10:>8.4f}"
              f"{100*t/max(n,1):>21.1f}%")
    print()

    print("REACH BY n_real TERCILE (f = 1)")
    W0 = Walk(paths[names[0]])
    nr = np.array([int(W0.real_sig[p].sum()) for p in range(W0.perts.size)])
    nr = nr[nr > 0]
    q = np.quantile(nr, [1 / 3, 2 / 3])
    print(f"  {'arm':<15}{'small':>10}{'mid':>10}{'large':>10}   "
          f"(n_real <= {q[0]:.0f} / <= {q[1]:.0f} / >)")
    for nm in names:
        W = Walk(paths[nm])
        per = []
        for p in range(W.perts.size):
            pool = np.flatnonzero(W.real_sig[p])
            if pool.size == 0:
                continue
            o = spec_order(pool, W.genes, W.pred_sig[p], W.pred_padj[p], W.pred_p[p],
                           W.pred_lfc[p])
            walked = pool[o]
            walked = walked[W.adj[p, walked]]
            ok = np.sign(W.pred_lfc[p, walked]) == np.sign(W.real_lfc[p, walked])
            per.append(kstar(ok)[0] / pool.size)
        per = np.array(per)
        print(f"  {nm:<15}{per[nr <= q[0]].mean():>10.5f}"
              f"{per[(nr > q[0]) & (nr <= q[1])].mean():>10.5f}"
              f"{per[nr > q[1]].mean():>10.5f}")
    print()

    if "ks000" in rows and "ks000_b" in rows:
        a, b = rows["ks000"], rows["ks000_b"]
        print("SEED SPREAD (ks000 vs its re-seeded twin) — every effect must beat this")
        print(f"  reach {abs(a['reach'][0]-b['reach'][0]):.5f}"
              f"   fid {abs(a['fid']-b['fid']):.5f}   jac {abs(a['jac']-b['jac']):.5f}"
              f"   nmae {abs(a['nmae']-b['nmae']):.5f}   pds {abs(a['pds']-b['pds']):.5f}")


if __name__ == "__main__":
    main()

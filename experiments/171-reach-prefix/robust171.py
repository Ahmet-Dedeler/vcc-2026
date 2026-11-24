"""171 — is the knee arm's reach gain carried by a handful of perturbations?

`reach` is a mean over perturbations of `k*/n_real`, and h1-test's `n_real` runs 1-5,484, so a
single perturbation with `n_real = 1` that flips from k* = 0 to k* = 1 moves the panel mean by
1/85 = 0.0118 -- more than the whole measured effect.  A mean that survives its seed twin can
still be one perturbation.

Three tests: the paired per-perturbation delta with a t statistic; the delta with the single
most influential perturbation deleted; and the delta computed on the perturbations with
`n_real >= 20` only (the ones the challenge's regime resembles).

Also prices `fid` two ways, because they DISAGREE and the disagreement is the arm's main risk:
164's conserved-headroom projection (optimistic) and the board's own five-point `fid`-vs-n_pred
series (direct, and pointing the other way).
"""
from __future__ import annotations

import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
sys.path.insert(0, f"{ROOT}/experiments/171-reach-prefix")
from debench import Surface  # noqa: E402
from diag171 import Walk, kstar, spec_order  # noqa: E402

W = f"{ROOT}/experiments/171-reach-prefix/work"
D164 = f"{ROOT}/experiments/164-callset-artefact"
REACH_SPAN = 0.96414 - 0.07902
FID_SPAN = 0.80964 - 0.51239


def per_pert(path, f=1.0):
    Wk = Walk(path)
    z = np.load(path, allow_pickle=True)
    if f == 1.0:
        sig = Wk.real_sig
    else:
        rp = np.where(np.isfinite(z["real_p"]), z["real_p"], np.inf)
        sig = np.zeros(rp.shape, bool)
        for i in range(rp.shape[0]):
            k = int(round(f * int(Wk.real_sig[i].sum())))
            k = max(0, min(k, int(np.isfinite(z["real_p"][i]).sum())))
            if k:
                sig[i, np.argsort(rp[i], kind="stable")[:k]] = True
        sig &= np.isfinite(z["real_p"])
        for p, j in enumerate(Wk.self_idx):
            if j >= 0:
                sig[p, j] = False
    R, NR, IDX = [], [], []
    for p in range(Wk.perts.size):
        pool = np.flatnonzero(sig[p])
        if pool.size == 0:
            continue
        o = spec_order(pool, Wk.genes, Wk.pred_sig[p], Wk.pred_padj[p], Wk.pred_p[p],
                       Wk.pred_lfc[p])
        w = pool[o]
        w = w[Wk.adj[p, w]]
        ok = (np.sign(Wk.pred_lfc[p, w]) == np.sign(Wk.real_lfc[p, w])) if w.size else \
            np.zeros(0, bool)
        R.append(kstar(ok)[0] / pool.size)
        NR.append(pool.size)
        IDX.append(p)
    return np.array(R), np.array(NR), np.array(IDX)


def main():
    for f in (1.0, 2.0):
        a1, nr, _ = per_pert(f"{W}/de171.a15.npz", f)
        a2, _, _ = per_pert(f"{D164}/de164.a15_b.npz", f)
        k1, _, _ = per_pert(f"{W}/de171.kt10000.npz", f)
        k2, _, _ = per_pert(f"{W}/de171.kt10000_b.npz", f)
        a, k = (a1 + a2) / 2, (k1 + k2) / 2
        d = k - a
        t = d.mean() / (d.std(ddof=1) / np.sqrt(d.size))
        j = int(np.argmax(np.abs(d)))
        keep = np.ones(d.size, bool)
        keep[j] = False
        big = nr >= 20
        print(f"=== kt10000 - a15, twin-averaged, proportional power f = {f}")
        print(f"  panel mean delta                {d.mean():+.5f}   "
              f"(n = {d.size} perturbations)")
        print(f"  paired t                        {t:+.2f}   "
              f"({int((d > 0).sum())} up / {int((d < 0).sum())} down / "
              f"{int((d == 0).sum())} unchanged)")
        print(f"  most influential perturbation   n_real {nr[j]}, delta {d[j]:+.5f} "
              f"({100*abs(d[j])/abs(d.sum()):.0f}% of the panel sum)")
        print(f"  delta with it deleted           {d[keep].mean():+.5f}")
        print(f"  delta on n_real >= 20 only      {d[big].mean():+.5f}  "
              f"(n = {int(big.sum())})")
        print(f"  top-5 |delta| share of the sum  "
              f"{100*np.abs(d)[np.argsort(-np.abs(d))[:5]].sum()/abs(d.sum()):.0f}%")
        print()

    print("=== `fid`, two models that DISAGREE")
    Sa, Sk = Surface(f"{W}/de171.a15.npz"), Surface(f"{W}/de171.kt10000.npz")
    na = np.median([int(Sa.pred_sig[p].sum()) for p in range(Sa.perts.size)])
    nk = np.median([int(Sk.pred_sig[p].sum()) for p in range(Sk.perts.size)])
    ratio = nk / na
    chal = 4424 * ratio
    print(f"  h1-test median n_pred  a15 {na:.0f} -> kt10000 {nk:.0f}  (x{ratio:.3f})")
    print(f"  implied challenge n_pred  4424 -> {chal:.0f}")
    p_after = 0.5 + 0.00220 * (4424 / chal)
    d1 = (p_after - 0.50220) / FID_SPAN / 6
    print(f"  (a) 164's conserved-headroom projection: fid 0.50220 -> {p_after:.5f}"
          f"   d overall {d1:+.6f}")
    series = np.array([[1563, 0.49126], [1911, 0.49654], [3836, 0.5021],
                       [4424, 0.5022], [5961, 0.50178]])
    f_at = np.interp(chal, series[:, 0], series[:, 1])
    d2 = (f_at - 0.5022) / FID_SPAN / 6
    print(f"  (b) the board's OWN five-point fid-vs-n_pred series, interpolated at "
          f"{chal:.0f}: {f_at:.5f}   d overall {d2:+.6f}")
    print(f"      series: " + "  ".join(f"{int(n)}->{v:.5f}" for n, v in series))
    print(f"  => the two models differ by {abs(d1-d2):.6f} overall and DISAGREE ON SIGN.")
    print(f"     164 section 4: 'shrinking the call set has never once raised board fid',")
    print(f"     measured five times.  (b) is a direct board measurement of this exact")
    print(f"     family (emission-only, fixed profile); (a) is a projection.  Take (b).")


if __name__ == "__main__":
    main()

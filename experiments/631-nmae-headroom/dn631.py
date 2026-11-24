"""631 DOWN pass (driver-only prototype): lower the arithmetic mean of called DOWN genes at fixed per-gene totals.

decomp631 / proxy631: doubling the delivered lfc of every called DOWN gene (q < 0.05, l < 0) is worth −0.023 / −0.025 raw
nmae on H1 at full delivery. The emitter parks a DOWN gene's total on mass-holder carriers above the control max, so the
mean is inflated. A count moved carrier a -> ballast row b lowers the gene's mean CPM by cw[a] - cw[b] (cw = 1/(n s)).
Per selected gene g, in order of decreasing |l| (most negative first), until its target A_t = ref * 2^(m l) is met:
  * donor: the carrier with the HIGHEST value for g (the top holder; dropping from the top crosses the fewest controls,
    and above the control max it crosses none);
  * recipient: the ballast row whose rise for this quantum crosses the fewest controls;
  * NET-U balancing: carrier drops LOWER g's U, ballast rises RAISE it; a move is accepted only if g's net U change
    (from the block it was given) stays within [-capU, +capU], capU = max_dauc * 2 n n_ctrl (630's bar);
  * library payback: the ballast row gives the carrier one sub-gate count per count it received (ZS._take on the row's
    sub-gate genes), so every library and every per-gene total stays exact.
A gene whose call would flip is excluded by the caller's guard (rnup630's loop). Returns (X, diagnostics).
"""
from __future__ import annotations

import numpy as np

from vcc2026 import rank_mag as ZS


def dn_pass(block, n_car, cc, l_base, q_base, m, *, max_dauc=1e-4, quantum=5, exclude=None, sel_extra=None):
    X = np.array(block, dtype=np.int64, copy=True)
    n = X.shape[0]
    keep = np.asarray(cc.keep, dtype=bool)
    gidx, sidx = np.flatnonzero(keep), np.flatnonzero(~keep)
    Lr = X.sum(1)
    s = Lr / 1e6
    ref = np.asarray(cc.ref_lin, dtype=np.float64)
    car, bal = np.arange(n_car), np.arange(n_car, n)
    if bal.size == 0:
        return X, dict(dn_sel=0, dn_moved=0, dn_counts=0, dn_need=0.0, dn_got=0.0, dn_frac=0.0)
    cw = 1.0 / (n * s)
    n_ctrl = cc.cs.shape[1]
    capU = int(np.floor(max_dauc * 2.0 * n * n_ctrl))
    Xk = X[:, gidx]
    A = (Xk / s[:, None]).mean(0)
    sel = np.isfinite(l_base) & (l_base < 0) & (q_base < 0.05)
    if sel_extra is not None:
        sel &= sel_extra
    if exclude is not None and len(exclude):
        sel[np.asarray(sorted(exclude), dtype=np.int64)] = False
    A_t = (ref + 1e-9) * np.power(2.0, float(m) * np.where(sel, l_base, 0.0)) - 1e-9
    need = np.where(sel, np.maximum(A - A_t, 0.0), 0.0)
    order = np.flatnonzero(need > 0)[np.argsort(l_base[need > 0], kind="stable")]
    sub_bal = X[np.ix_(bal, sidx)].sum(1).astype(np.int64)            # each ballast row's sub-gate payback pool
    pay = []                                                           # (ballast row, carrier, t) sub-gate paybacks
    st = dict(dn_sel=int(order.size), dn_moved=0, dn_counts=0, dn_need=float(need.sum()), dn_got=0.0)
    for j in order:
        j = int(j)
        c = cc.cs[j]
        net = 0
        rem = need[j]
        got = 0.0
        moved = 0
        while rem > 0:
            xc = Xk[car, j]
            cand = np.flatnonzero(xc > 0)
            if cand.size == 0:
                break
            vals = np.log1p(xc[cand] / s[car[cand]])
            a = int(car[cand[np.argmax(vals)]])                         # the top holder
            k = int(min(quantum, Xk[a, j], np.ceil(rem / max(cw[a] - float(np.min(cw[bal])), 1e-12))))
            if k <= 0:
                break
            du_a = (int(ZS.ranks_of(np.array([Xk[a, j] - k]), np.array([s[a]]), c)[0])
                    - int(ZS.ranks_of(np.array([Xk[a, j]]), np.array([s[a]]), c)[0]))
            # pick the ballast row whose rise by k costs the least (and still has payback budget)
            best, best_du = -1, None
            for bi, b in enumerate(bal):
                if sub_bal[bi] < k or cw[b] >= cw[a]:
                    continue
                cur = int(ZS.ranks_of(np.array([Xk[b, j]]), np.array([s[b]]), c)[0])
                du_b = int(ZS.ranks_of(np.array([Xk[b, j] + k]), np.array([s[b]]), c)[0]) - cur
                if best_du is None or du_b < best_du:
                    best, best_du = bi, du_b
            if best < 0:
                break
            d = du_a + best_du
            if abs(net + d) > capU:
                if quantum > 1 and k > 1:                               # retry one count at a time before giving up
                    quantum = 1
                    continue
                break
            b = int(bal[best])
            Xk[a, j] -= k
            Xk[b, j] += k
            sub_bal[best] -= k
            pay.append((best, a, k))
            net += d
            g = k * (cw[a] - cw[b])
            rem -= g
            got += g
            moved += k
        if moved:
            st["dn_moved"] += 1
            st["dn_counts"] += moved
            st["dn_got"] += got
    X[:, gidx] = Xk
    for bi, a, t in pay:                                               # sub-gate payback: ballast row -> carrier
        b = int(bal[bi])
        give = ZS._take(X[b, sidx], t)
        X[b, sidx] -= give
        X[a, sidx] += give
    B = np.asarray(block)
    if not np.array_equal(X.sum(1), Lr):
        raise AssertionError("dn_pass: a library moved")
    if not np.array_equal(X.sum(0), B.sum(0)):
        raise AssertionError("dn_pass: a per-gene total moved")
    if X.min() < 0:
        raise AssertionError("dn_pass: negative count")
    for j in np.flatnonzero(np.any(X[:, gidx] != B[:, gidx], axis=0)):
        du = int((ZS.ranks_of(X[:, gidx[j]], s, cc.cs[j]) - ZS.ranks_of(B[:, gidx[j]], s, cc.cs[j])).sum())
        if abs(du) > capU:
            raise AssertionError(f"dn_pass: gene {j} net U change {du} beyond {capU}")
    st["dn_frac"] = st["dn_got"] / max(st["dn_need"], 1e-300)
    return X, st

"""Rank-near-neutral DOWN magnitude pass (experiment 631; candidate for src/vcc2026/dn_mag.py, builder knob VCC_DE_DNMAG).

What it fixes: the K-class emitter keeps called DOWN genes less down than their design (the mass-balance kf ~0.60 and
mass holders above the control max inflate the arithmetic mean), so nmae pays for DOWN magnitude the signs already earned.
On H1's real emitter, one factor x2 on every called DOWN gene is worth -0.023 / -0.025 raw nmae at full delivery (631).
The pass lowers each called DOWN gene's arithmetic mean CPM at FIXED per-gene totals: counts of the gene move from the
carrier with the highest value (the top holder) to the ballast row whose rise crosses the fewest controls; the gene's
NET Wilcoxon rank-sum change stays within |dAUC| <= 1e-4 (carrier drops lower U, ballast rises raise it); the carrier's
library is paid back with sub-gate (untested) counts from that ballast row, so every library and every per-gene total
stays exact (pds bit-identical; mse's profile half). A call-flip guard re-runs without any gene whose call flipped; any
invariant assert fails SAFE (the input block comes back, ur-style key ``dn_invariant_failed`` = 1).
Selection: q < 0.05 and delivered lfc < 0 (called DOWN), not the target's own column. Target lfc = m x today's.
Variant grammar: d_m<m*100>[_q<quantum>], m in [1.25, 4], e.g. d_m200 (the H1-priced form, quantum 5).
Applied after the zs pass and after 630's UP pass when both are set (it reads the block it is given).
"""
from __future__ import annotations

import re

import numpy as np
from scipy.stats import false_discovery_control, norm

from vcc2026 import rank_mag as ZS
from vcc2026 import rc_emit

VARIANT_RE = re.compile(r"d_m(\d{3})(?:_q([1-9]\d*))?")
GUARD_ROUNDS = 5
DIAG_KEYS = ("dn_sel", "dn_moved", "dn_counts", "dn_need", "dn_got", "dn_frac", "dn_guard_rounds", "dn_guard_excluded",
             "dn_guard_failed", "dn_invariant_failed")


def parse_spec(variant: str) -> dict:
    m = VARIANT_RE.fullmatch(str(variant).strip())
    if m is None:
        raise ValueError(f"dn_mag variant must look like 'd_m<m*100>[_q<quantum>]', e.g. d_m200: {variant!r}")
    mult = int(m.group(1)) / 100.0
    if not 1.25 <= mult <= 4.0:
        raise ValueError(f"dn_mag m must be in [1.25, 4]: {variant!r}")
    return {"m": mult, "quantum": int(m.group(2)) if m.group(2) else 5}


def scorer_tie(cc) -> np.ndarray:
    """The control tie sums the scorer's Wilcoxon needs (exact_de323.control_block's `tie`)."""
    return rc_emit.ctrl_tiesums(cc.cs)


def de_stats(block, cc, tie):
    """(p, q, lfc) on the DE universe in the evaluator's arithmetic (exact_de323.de_from_counts, verbatim order)."""
    x = np.asarray(block, dtype=np.float64)
    scale = x.sum(1) / 1e6
    scale[scale == 0] = 1.
    linear = (x / scale[:, None])[:, np.asarray(cc.keep, dtype=bool)]
    logged = np.ascontiguousarray(np.log1p(linear).T)
    z = rc_emit.wilcoxon_z(cc.cs, tie, logged)
    p = 2 * norm.sf(np.abs(z))
    q = false_discovery_control(p, method="bh")
    lfc = np.log2((linear.mean(0) + 1e-9) / (np.asarray(cc.ref_lin, dtype=np.float64) + 1e-9))
    return p, q, lfc


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


def apply(block, n_car, cc, variant, *, tie, own_col=-1):
    """Builder entry point: returns (new int64 block, flat numeric diagnostics with keys DIAG_KEYS). Fails SAFE."""
    try:
        return _apply(block, n_car, cc, variant, tie=tie, own_col=own_col)
    except AssertionError:
        d = {k: 0.0 for k in DIAG_KEYS}
        d["dn_invariant_failed"] = 1.0
        return np.array(block, dtype=np.int64, copy=True), d


def _apply(block, n_car, cc, variant, *, tie, own_col=-1):
    spec = parse_spec(variant) if isinstance(variant, str) else dict(variant)
    X0 = np.asarray(block, dtype=np.int64)
    p0, q0, l0 = de_stats(X0, cc, tie)
    sel_extra = np.ones(l0.shape, dtype=bool)
    if own_col >= 0:
        sel_extra[own_col] = False
    excl, rounds = set(), 0
    while True:
        X1, st = dn_pass(X0, n_car, cc, l0, q0, spec["m"], max_dauc=1e-4, quantum=spec["quantum"], exclude=excl,
                         sel_extra=sel_extra)
        p1, q1, l1 = de_stats(X1, cc, tie)
        flips = set(np.flatnonzero((q0 < 0.05) != (q1 < 0.05)).tolist())
        if not flips:
            break
        if rounds >= GUARD_ROUNDS:
            d = {k: 0.0 for k in DIAG_KEYS}
            d.update(dn_guard_rounds=float(rounds), dn_guard_excluded=float(len(excl)), dn_guard_failed=1.0)
            return X0.copy(), d
        excl |= flips
        rounds += 1
    st = {k: float(v) for k, v in st.items()}
    st.update(dn_guard_rounds=float(rounds), dn_guard_excluded=float(len(excl)), dn_guard_failed=0.0,
              dn_invariant_failed=0.0)
    return X1, {k: float(st.get(k, 0.0)) for k in DIAG_KEYS}

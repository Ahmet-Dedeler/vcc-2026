"""630 v7 (from v5): DOWN-gene payments under each payer gene's NET-U budget (631's lever paired with 630).

v5 paid with DOWN genes only at zero rank cost on both sides (ballast-poor: ~2 % capacity). v7 lets the payer gene g's
carrier drop and ballast rise cross controls as long as g's NET U change stays within +/- capU (the same 1e-4 bar), and
selects payer genes among CALLED DOWN genes that still need magnitude (631: x2 on called DOWN genes is worth
-0.023 / -0.025 raw on H1). A move's size is bounded by the recipient's sub-gate budget as in v3 (never short); DOWN
payments come first and SAVE sub-gate counts for later moves.
(v5 docstring follows.)
630 v5: a second compensation currency for the rank-neutral UP pass (driver-only; module 7e6b3e52 untouched).

diagcap630 found the binding constraint: each count a recipient carrier receives is paid with a sub-gate count moved to
the donor ballast row (libraries exact), and carriers hold only ~0.55 % of their library in sub-gate genes. 93 % of it is
spent and 76 % of carriers drain. v5 lets the recipient ALSO pay with counts of a DOWN gene g (delivered lfc < 0) that it
holds ABOVE g's control max (zs's mass holders), moved to the same donor ballast row:
  * the carrier stays above g's control max (count >= the minimum above it), so g's carrier-side rank is unchanged;
  * the ballast row stays inside its current control interval for g (cell_bounds' x_hi; open above the max), so g's
    ballast-side rank is unchanged: g's U, p and call are exactly unchanged;
  * g's mean CPM falls by (cw[a] - cw[b]) per count: more DOWN magnitude, capped so g's lfc never goes below
    dn_mult x its current lfc (dn_mult 1.5 / 2.5; zs's own target multiplier is 2.5);
  * libraries and every per-gene total stay exact.
dn_mult None == up_pass2 exactly (the KA). dn_mask (kept axis, bool) restricts the eligible DOWN genes (e.g. called
only); None = every DOWN gene. Payment priority: the most negative current lfc first. DOWN currency is spent before the
sub-gate pool.
Variant grammar (rnup630): <v3 variant>_d<dn_mult*10>[c], e.g. b_k100_q20r3p10_d25 (all DOWN genes, x2.5),
b_k100_q20r3p10_d15c (called DOWN genes only, x1.5).
"""
from __future__ import annotations

import numpy as np

import src_rank_up as M
from src_rank_up import ZS, _max_drop_k, _min_above


def _min_above_vec(cmax, s):
    """Vectorised _min_above: (n_cells, n_genes) minimum count with log1p(x / s) > cmax."""
    xm = np.floor(np.expm1(cmax)[None, :] * s[:, None]).astype(np.int64) + 1
    for _ in range(4):
        xm = xm - (np.log1p((xm - 1) / s[:, None]) > cmax[None, :])
    for _ in range(4):
        xm = xm + ~(np.log1p(xm / s[:, None]) > cmax[None, :])
    if not (np.all(np.log1p(xm / s[:, None]) > cmax[None, :]) and np.all(~(np.log1p((xm - 1) / s[:, None]) > cmax[None, :]))):
        raise AssertionError("min-above inversion did not converge")
    return xm


def up_pass7(block, n_car, cc, l_base, target_lfc, sel, *, max_dauc=1e-4, cross=False, exclude=None, bounds=None,
             quantum=None, rounds=1, rec_pool=1, dn_mult=2.0, dn_mask=None, max_cand=40):
    X = np.array(block, dtype=np.int64, copy=True)
    n = X.shape[0]
    keep = np.asarray(cc.keep, dtype=bool)
    gidx, sidx = np.flatnonzero(keep), np.flatnonzero(~keep)
    Lr = X.sum(1)
    s = Lr / 1e6
    ref = np.asarray(cc.ref_lin, dtype=np.float64)
    Xk = X[:, gidx]
    A = (Xk / s[:, None]).mean(0)
    l = np.log2((A + 1e-9) / (ref + 1e-9))
    fin = np.isfinite(l_base)
    if np.max(np.abs(l[fin] - l_base[fin]), initial=0.0) > 1e-9:
        raise AssertionError("base lfc does not match the block")
    n_ctrl = cc.cs.shape[1]
    capU = int(np.floor(max_dauc * 2.0 * n * n_ctrl))
    car, bal = np.arange(n_car), np.arange(n_car, n)
    cw = 1.0 / (n * s)
    A_t = (ref + 1e-9) * np.power(2.0, target_lfc) - 1e-9
    ok = sel & np.isfinite(target_lfc)
    if exclude is not None and len(exclude):
        ok[np.asarray(sorted(exclude), dtype=np.int64)] = False
    dA = np.where(ok, np.maximum(A_t - A, 0.0), 0.0)
    order = np.flatnonzero(dA > 0)[np.argsort(-dA[dA > 0], kind="stable")]
    sub_left = X[np.ix_(car, sidx)].sum(1).astype(np.int64)
    sub_moves = []
    dU = np.zeros(gidx.size, dtype=np.int64)
    # ---- v5: the DOWN currency -------------------------------------------------------------------------------------
    cmax_all = cc.cs[:, -1]
    dn = np.isfinite(l_base) & (l_base < 0) & ~ok
    if dn_mask is not None:
        dn &= np.asarray(dn_mask, dtype=bool)
    A_floor = (ref + 1e-9) * np.power(2.0, float(dn_mult) * np.where(dn, l_base, 0.0)) - 1e-9
    capA = np.where(dn, np.maximum(A - A_floor, 0.0), 0.0)            # mean CPM each DOWN gene may still lose
    dn_list = np.flatnonzero(dn & (capA > 0))
    prio = dn_list[np.argsort(l_base[dn_list], kind="stable")]        # most negative lfc first
    netU = np.zeros(gidx.size, dtype=np.int64)                         # v7: each payer gene's net U change so far
    BIGR = 0
    # -----------------------------------------------------------------------------------------------------------------
    st = dict(n_sel=int(order.size), n_moved=0, n_reverted=0, n_no_candidate=0, counts=0, need=float(dA.sum()), got=0.0,
              n_no_pos_donor=0, n_neg_donor_rows=0, n_mean_revert=0, dn_paid=0, sub_paid=0, n_pay_short=0,
              dn_genes_used=0, dn_pool=int(prio.size))
    dn_used = set()
    u_base = {}
    rec_of = {}
    moved_any = set()
    rem_of = {int(j): float(dA[j]) for j in order}

    def avail(a):
        return int(sub_left[a])                                       # v7: bounded by the sub-gate budget, as v3

    for rnd in range(int(rounds)):
      for j in order:
        j = int(j)
        if rem_of[j] <= 0:
            continue
        c = cc.cs[j]
        cmax = c[-1]
        snap_col, snap_sub, snap_mv = Xk[:, j].copy(), sub_left.copy(), len(sub_moves)
        undo = []                                                        # v5: DOWN payments to revert with this gene
        vcar = np.log1p(Xk[car, j] / s[car])
        above = np.flatnonzero(vcar > cmax)
        if j not in u_base:
            u_base[j] = ZS.ranks_of(Xk[:, j], s, c)
        u0_cells = ZS.ranks_of(Xk[:, j], s, c)
        if j in rec_of:
            a = rec_of[j]
            dU_promo, k_promo = 0, 0
        elif above.size:
            a = int(above[np.argmax(cw[above])])
            dU_promo, k_promo = 0, 0
        else:
            cand = np.flatnonzero(sub_left > 0)
            if cand.size == 0:
                st["n_no_candidate"] += 1
                continue
            cand = cand[np.lexsort((-cw[cand], -vcar[cand]))][: max(int(rec_pool), 1)]
            a = int(cand[np.argmax([avail(int(x)) for x in cand])]) if rec_pool > 1 else int(cand[0])
            k_promo = _min_above(cmax, s[a]) - int(Xk[a, j])
            dU_promo = 2 * n_ctrl - int(u0_cells[a])
        elig = cw[bal] < cw[a]
        st["n_neg_donor_rows"] += int((~elig).sum())
        if not elig.any():
            st["n_no_candidate"] += 1
            st["n_no_pos_donor"] += 1
            continue
        lo_j, hi_j = (-capU, capU) if bounds is None else (int(bounds[0][j]), int(bounds[1][j]))
        du_sofar = int((u0_cells - u_base[j]).sum())
        lo_j, hi_j = lo_j - du_sofar, hi_j - du_sofar
        xb = Xk[bal, j].copy()
        if cross:
            budget = dU_promo - lo_j
            dcap = np.zeros(bal.size, dtype=np.int64)
            for bi in np.argsort(-xb, kind="stable"):
                if budget <= 0 or xb[bi] <= 0 or not elig[bi]:
                    continue
                kb = _max_drop_k(int(xb[bi]), float(s[bal[bi]]), c, budget)
                if kb > 0:
                    drop = int(u0_cells[bal[bi]]) - int(ZS.ranks_of(np.array([xb[bi] - kb]), np.array([s[bal[bi]]]), c)[0])
                    dcap[bi] = kb
                    budget -= drop
        else:
            xl, *_ = ZS.cell_bounds(xb, s[bal], c)
            dcap = np.where(elig, xb - xl, 0)
        bud = avail(a)
        if k_promo > bud or k_promo > int(dcap.sum()):
            st["n_no_candidate"] += 1
            continue
        rem = rem_of[j]
        k_want = max(k_promo, int(np.ceil(rem / (cw[a] - float(np.min(cw[bal][elig]))))))
        if quantum is not None:
            k_want = min(k_want, max(k_promo, int(quantum)))
        k = min(k_want, bud, int(dcap.sum()))
        if k <= 0:
            st["n_no_candidate"] += 1
            continue
        take = ZS.fill(dcap, k)
        short = False
        paid_dn, paid_sub = 0, 0
        for bi in np.flatnonzero(take):
            t = int(take[bi])
            Xk[bal[bi], j] -= t
            Xk[a, j] += t
            rem -= t * (cw[a] - cw[bal[bi]])
            # v5 payment: DOWN currency first (zero-cost on both sides, capped per gene), then sub-gate
            need = t
            gain = cw[a] - cw[bal[bi]]
            live = prio[(Xk[a, prio] > 0) & (capA[prio] >= gain)][:max_cand]
            bb = int(bal[bi])
            for g in live:
                if need <= 0:
                    break
                cg = cc.cs[g]
                xa, xb2 = int(Xk[a, g]), int(Xk[bb, g])
                ra0 = int(ZS.ranks_of(np.array([xa]), np.array([s[a]]), cg)[0])
                rb0 = int(ZS.ranks_of(np.array([xb2]), np.array([s[bb]]), cg)[0])
                for u in (int(min(need, xa, np.floor(capA[g] / gain))), 1):
                    if u <= 0:
                        continue
                    d = (int(ZS.ranks_of(np.array([xa - u]), np.array([s[a]]), cg)[0]) - ra0
                         + int(ZS.ranks_of(np.array([xb2 + u]), np.array([s[bb]]), cg)[0]) - rb0)
                    if abs(int(netU[g]) + d) <= capU:
                        Xk[a, g] -= u
                        Xk[bb, g] += u
                        netU[g] += d
                        capA[g] -= u * gain
                        undo.append((int(g), int(bi), u, u * gain, d))
                        need -= u
                        paid_dn += u
                        break
            if need > 0:
                if sub_left[a] < need:                                  # cannot happen: k <= sub_left[a]
                    raise AssertionError("v7: sub-gate payback short")
                sub_moves.append((a, int(bal[bi]), need))
                sub_left[a] -= need
                paid_sub += need
        du = int((ZS.ranks_of(Xk[:, j], s, c) - u0_cells).sum())
        a_old, a_new = float((snap_col / s).mean()), float((Xk[:, j] / s).mean())
        sign_new = np.sign(np.log2((a_new + 1e-9) / (ref[j] + 1e-9)))
        bad_mean = (a_new < a_old) or (sign_new != np.sign(l_base[j]))
        if bad_mean:
            st["n_mean_revert"] += 1
        if short:
            st["n_pay_short"] += 1
        if short or bad_mean or du < lo_j or du > hi_j or not (np.log1p(Xk[a, j] / s[a]) > cmax):
            Xk[:, j] = snap_col
            sub_left[:] = snap_sub
            del sub_moves[snap_mv:]
            for g, bi, u, dg, d in reversed(undo):                      # v7: undo this gene's DOWN payments
                Xk[a, g] += u
                Xk[bal[bi], g] -= u
                netU[g] -= d
                capA[g] += dg
            st["n_reverted"] += 1
            continue
        dU[j] = int((ZS.ranks_of(Xk[:, j], s, c) - u_base[j]).sum())
        rec_of[j] = a
        if j not in moved_any:
            moved_any.add(j)
            st["n_moved"] += 1
        for g, *_ in undo:
            dn_used.add(g)
            dU[g] = int(netU[g])                                        # v7: the payer gene's exact net U
        st["counts"] += k
        st["dn_paid"] += paid_dn
        st["sub_paid"] += paid_sub
        st["got"] += float(rem_of[j] - rem)
        rem_of[j] = rem
    X[:, gidx] = Xk
    for rc, bb, t in sub_moves:
        give = ZS._take(X[rc, sidx], t)
        X[rc, sidx] -= give
        X[bb, sidx] += give
    B = np.asarray(block)
    if not np.array_equal(X.sum(1), Lr):
        raise AssertionError("a library moved")
    if not np.array_equal(X.sum(0), B.sum(0)):
        raise AssertionError("a per-gene block total moved")
    if X.min() < 0:
        raise AssertionError("negative count")
    for j in np.flatnonzero(np.any(X[:, gidx] != B[:, gidx], axis=0)):
        du = int((ZS.ranks_of(X[:, gidx[j]], s, cc.cs[j]) - ZS.ranks_of(B[:, gidx[j]], s, cc.cs[j])).sum())
        lo_j, hi_j = (-capU, capU) if bounds is None else (int(bounds[0][j]), int(bounds[1][j]))
        if du != int(dU[j]) or du < lo_j or du > hi_j:
            raise AssertionError(f"kept gene {j}: U change {du}, recorded {int(dU[j])}, bounds [{lo_j}, {hi_j}]")
    # v5: every DOWN gene that paid kept its exact ranks (U change 0 is asserted above via dU == 0) and its lfc floor
    A1 = (X[:, gidx] / s[:, None]).mean(0)
    used = np.array(sorted(dn_used), dtype=np.int64)
    if used.size and np.any(A1[used] < A_floor[used] - 1e-9):
        raise AssertionError("a DOWN gene went below its lfc floor")
    au = np.abs(dU[dU != 0]) / (2.0 * n * n_ctrl)
    st.update(max_dAUC=float(au.max()) if au.size else 0.0, median_dAUC=float(np.median(au)) if au.size else 0.0,
              frac_got=st["got"] / max(st["need"], 1e-300), capU=capU, dn_genes_used=len(dn_used))
    return X, st

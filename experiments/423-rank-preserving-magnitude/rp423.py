"""423: rank-preserving magnitude post-pass on one emitted target block (carriers + ballast).

What the scorer reads from a 400-cell block, per DE-universe gene g (the reference-control CPM > 5 gate):
  * Wilcoxon z against the 18,400 reference controls, on log1p(count / (library / 1e6)).  The rank sum is
    U = sum_i (#controls < v_i + 0.5 #controls == v_i) (+ a constant), and the tie term depends on the tie groups.
  * the arithmetic-mean CPM A_g = mean_i count_ig / (library_i / 1e6), which sets lfc = log2(A / ref_lin).
  * pds/mse read only the per-gene block totals (and mse's jackknife).

The pass scales each selected gene's lfc toward k * lfc (A* = ref_lin * 2^(k lfc)) using only moves that keep
every carrier cell strictly inside its current interval between adjacent control values ("per-cell rank-free";
a cell above the control max has an open interval, so it can rise without limit; a cell tied with a control
value - e.g. a zero count against control zeros - is frozen):
  * carrier libraries are held fixed: every count added to (removed from) gene g in carrier cell i is taken
    from (given to) that cell's sub-gate genes (not tested, not in nmae) - so no other tested gene's CPM moves;
  * the ballast gives back / absorbs exactly the counts moved, gene by gene (design mode for touched genes),
    with its own library held fixed through the mirrored sub-gate counts - so every per-gene block total, and
    hence pds and the profile half of mse, is bit-identical.
The ballast is one (or K) very deep cell(s) and its value for gene g moves by the absorbed counts. In mode "pc"
that move is allowed to cross control values (its U change is measured); in mode "bal" carrier cells are
additionally moved ACROSS control values in the same direction as the magnitude change until the gene's total
U is restored exactly (the ballast's crossing is offset by the carriers'); a gene whose U cannot be restored
exactly is rolled back.

Variant grammar: <sel>_k<k*100>[_<mode>], sel in all|q50|q80|up|dn (quantiles of |lfc| over the DE universe),
mode in pc|bal (default pc).  Example: q80_k140_bal.
"""
from __future__ import annotations

import numpy as np

SEL = ("all", "q50", "q80", "up", "dn", "q80up", "q80dn")
BIG = np.int64(1) << 40


def parse_variant(v: str) -> dict:
    """<sel>_k<K*100>[u<Kup*100>][_<mode>]; the optional u part gives up-genes their own factor (sel 'all' only)."""
    parts = v.split("_")
    sel, kk = parts[0], parts[1]
    mode = parts[2] if len(parts) > 2 else "pc"
    if sel not in SEL or not kk.startswith("k") or mode not in ("pc", "zs", "bal", "ex", "zspc"):
        raise ValueError(f"bad variant {v}")
    k_up = None
    if "u" in kk:
        kk, ku = kk.split("u")
        k_up = int(ku) / 100.0
        if sel != "all":
            raise ValueError("a separate up factor needs sel 'all'")
    return {"sel": sel, "k": int(kk[1:]) / 100.0, "k_up": k_up, "mode": mode}


# ------------------------------------------------------------------ exact per-cell interval arithmetic
def cell_bounds(x, s, c):
    """For counts x (n,) at scales s (n,) (library/1e6, the scorer's float) against sorted control values c (m,):
    returns (x_lo, x_hi, lt, rt, tied). Any count in [x_lo, x_hi] keeps log1p(count/s) strictly inside the same
    open interval between adjacent control values (so #controls below, and no tie, are unchanged). x_hi = BIG
    above the control max. Tied cells get x_lo = x_hi = x (frozen)."""
    x = np.asarray(x, dtype=np.int64)
    v = np.log1p(x / s)
    m = c.size
    lt = np.searchsorted(c, v, "left")
    rt = np.searchsorted(c, v, "right")
    tied = rt > lt
    hi = np.where(rt < m, c[np.minimum(rt, m - 1)], np.inf)
    lo = np.where(lt > 0, c[np.maximum(lt - 1, 0)], -np.inf)
    xh = np.full(x.shape, BIG, dtype=np.int64)
    f = np.isfinite(hi)
    if f.any():
        sf, hf = s[f], hi[f]
        cand = np.floor(np.expm1(hf) * sf).astype(np.int64)
        for _ in range(4):
            cand = cand - (np.log1p(cand / sf) >= hf)
        for _ in range(4):
            cand = cand + (np.log1p((cand + 1) / sf) < hf)
        if not (np.all(np.log1p(cand / sf) < hf) and np.all(np.log1p((cand + 1) / sf) >= hf)):
            raise AssertionError("x_hi inversion did not converge")
        xh[f] = cand
    xl = np.zeros(x.shape, dtype=np.int64)
    f = np.isfinite(lo)
    if f.any():
        sf, lf = s[f], lo[f]
        cand = np.maximum(np.ceil(np.expm1(lf) * sf).astype(np.int64), 0)
        for _ in range(4):
            cand = cand + (np.log1p(cand / sf) <= lf)
        for _ in range(4):
            cand = cand - ((cand >= 1) & (np.log1p(np.maximum(cand - 1, 0) / sf) > lf))
        if not (np.all(np.log1p(cand / sf) > lf) and np.all((cand == 0) | (np.log1p((cand - 1) / sf) <= lf))):
            raise AssertionError("x_lo inversion did not converge")
        xl[f] = cand
    xl = np.where(tied, x, xl)
    xh = np.where(tied, x, xh)
    if np.any(xl > x) or np.any(xh < x):
        raise AssertionError("a cell lies outside its own interval")
    return xl, xh, lt, rt, tied


def ranks_of(x, s, c):
    """U contribution (x2, to stay integer) of each cell: 2*#controls below + #controls equal."""
    v = np.log1p(np.asarray(x, dtype=np.int64) / s)
    lt = np.searchsorted(c, v, "left")
    rt = np.searchsorted(c, v, "right")
    return 2 * lt + (rt - lt)


def _pos(x, s, c):
    v = np.log1p(np.asarray(x, dtype=np.int64) / s)
    lt = np.searchsorted(c, v, "left")
    return lt, np.searchsorted(c, v, "right") - lt


def _offset(x, s, c, R, room, cw, sgn, T=64, frozen=None):
    """Zero-sum carrier pairs on this gene (cell a +t, cell b -t, t <= T) whose control crossings sum to exactly R
    (2x rank units: + = carriers must gain). One side of a pair may be rank-free (a cell above the control max, or
    moving inside its interval), so a single cell can walk down into the sparse top tail of the controls and cancel
    the ballast's crossing count by count. Block totals and the ballast are untouched; each cell is used once;
    no cell lands on a control tie. Returns (x_new, n_pairs) or None when R cannot be hit exactly."""
    x = np.asarray(x, dtype=np.int64).copy()
    n = x.size
    lt0, eq0 = _pos(x, s, c)
    usable = eq0 == 0
    tt = np.arange(1, T + 1)[:, None]
    ltu, equ = _pos((x[None, :] + tt).ravel(), np.broadcast_to(s, (T, n)).ravel(), c)
    ltu, equ = ltu.reshape(T, n), equ.reshape(T, n)
    G = np.where(usable[None, :] & (equ == 0) & (room[None, :] >= tt), 2 * (ltu - lt0[None, :]), -1)
    xd = np.maximum(x[None, :] - tt, 0)
    ltd, eqd = _pos(xd.ravel(), np.broadcast_to(s, (T, n)).ravel(), c)
    ltd, eqd = ltd.reshape(T, n), eqd.reshape(T, n)
    Lo = np.where(usable[None, :] & (x[None, :] - tt >= 0) & (eqd == 0), 2 * (lt0[None, :] - ltd), -1)
    used = np.zeros(n, dtype=bool)
    if frozen is not None:
        used[frozen] = True
    npairs = 0
    for _ in range(n // 2):
        if R == 0:
            break
        best = None                                # (score tuple, t, a, b, net)
        for t in range(1, T + 1):
            ga = np.flatnonzero((G[t - 1] >= 0) & ~used)
            lb = np.flatnonzero((Lo[t - 1] >= 0) & ~used)
            if ga.size == 0 or lb.size == 0:
                continue
            lv = Lo[t - 1][lb]
            o = np.argsort(lv, kind="stable")
            lb, lv = lb[o], lv[o]
            gv = G[t - 1][ga]
            want = gv - R                          # exact: G_a - Lo_b == R
            i0 = np.searchsorted(lv, want, "left")
            i1 = np.searchsorted(lv, want, "right")
            for h in np.flatnonzero(i1 > i0):
                a = ga[h]
                for bi in range(i0[h], min(i1[h], i0[h] + 4)):
                    b = lb[bi]
                    if b == a:
                        continue
                    sc = (2, sgn * (cw[a] - cw[b]) * t, -t)
                    if best is None or sc > best[0]:
                        best = (sc, t, a, b, R)
                    break
            if best is not None and best[0][0] == 2:
                continue
            if R > 0:                              # progress: 0 < net <= R, as large as possible
                idx = np.searchsorted(lv, gv - R, "left")
                okp = idx < lv.size
                idx = np.minimum(idx, lv.size - 1)
                net = gv - lv[idx]
                okp &= (net > 0) & (net <= R) & (lb[idx] != ga)
            else:                                  # R <= net < 0
                idx = np.searchsorted(lv, gv - R, "right") - 1
                okp = idx >= 0
                idx = np.maximum(idx, 0)
                net = gv - lv[idx]
                okp &= (net < 0) & (net >= R) & (lb[idx] != ga)
            if okp.any():
                h = np.flatnonzero(okp)[np.argmax(np.abs(net[okp]))]
                sc = (1, abs(int(net[h])), -t)
                if best is None or sc > best[0]:
                    best = (sc, t, ga[h], lb[idx[h]], int(net[h]))
        if best is None:
            return None
        _, t, a, b, net = best
        x[a] += t
        x[b] -= t
        used[a] = used[b] = True
        R -= net
        npairs += 1
    if R != 0:
        return None
    return x, npairs


def _even_split(D, Kb):
    base = D // Kb
    return np.array([base + (1 if r < D % Kb else 0) for r in range(Kb)], dtype=np.int64)


def _drop_dump(x, xl, xh, room, xb, s_car, s_bal, c, cw, bw, rem, T=64):
    """Exact-U lowering for genes whose carriers cannot offset the ballast by pairs (e.g. a few mass-holder cells
    above the control max, every other carrier a zero tied with control zeros). One carrier h drops BELOW its
    interval to count y (losing lo = 2 x controls crossed); other carriers drop inside their intervals (m counts,
    lowest library first); of everything removed, D counts go to the ballast and the rest to one rank-free recipient
    (an above-max carrier in the highest-library class), so D is free and is chosen with ballast crossing
    d_b(D) == lo. When no D hits lo exactly, the closest configurations are finished with zero-sum pairs (_offset).
    Returns (x_new, xb_new) or None."""
    n = x.size
    Kb = xb.size
    lt0, eq0 = _pos(x, s_car, c)
    dn = x - xl
    effs = cw - float(np.mean(bw))
    Dcap = int(min(np.ceil(abs(rem) / max(float(effs.min()), 1e-12)) + 4 * T, dn.sum() + 4 * T))
    Ds = np.arange(Dcap + 1)
    rb0 = int(ranks_of(xb, s_bal, c).sum())
    db = np.zeros(Ds.size, dtype=np.int64)
    for r in range(Kb):
        share = Ds // Kb + (r < Ds % Kb)
        v = np.log1p((xb[r] + share) / s_bal[r])
        lt = np.searchsorted(c, v, "left"); rt = np.searchsorted(c, v, "right")
        db += 2 * lt + (rt - lt)
    db -= rb0
    if np.any(np.diff(db) < 0):
        return None
    order = np.argsort(-cw, kind="stable")                       # drops: lowest library first
    rec_ok = (xh >= BIG) & (room > 0)
    cands = []                                                    # (score, h, y, D, rho, m, residual)
    for h in np.flatnonzero((eq0 == 0) & (x >= 1)):
        ys = np.arange(max(1, int(xl[h]) - 256), int(xl[h]))  # every count strictly below the interval
        if ys.size == 0:
            continue
        lty, eqy = _pos(ys, np.full(ys.size, s_car[h]), c)
        lo = 2 * (lt0[h] - lty)
        rho_c = np.flatnonzero(rec_ok & (np.arange(n) != h))
        rho = int(rho_c[np.argmin(cw[rho_c])]) if rho_c.size else -1
        caps = dn.copy(); caps[h] = 0
        if rho >= 0:
            caps[rho] = 0
        capo = caps[order]
        cwo = cw[order]
        cumcap = np.cumsum(capo)
        cumA = np.cumsum(capo * cwo)
        F = int(cumcap[-1]) if cumcap.size else 0
        bw_mean = float(np.mean(bw))

        def dropA(m_):
            """CPM of the 400-cell mean removed by dropping m_ counts inside intervals, lowest library first."""
            if m_ <= 0:
                return 0.0
            i_ = int(np.searchsorted(cumcap, m_, "left"))
            prev_c = int(cumcap[i_ - 1]) if i_ > 0 else 0
            prev_a = float(cumA[i_ - 1]) if i_ > 0 else 0.0
            return prev_a + (m_ - prev_c) * float(cwo[i_])
        for k_ in np.flatnonzero(eqy == 0):
            y = int(ys[k_]); own = int(x[h]) - y
            i0 = int(np.searchsorted(db, lo[k_], "left")); i1 = int(np.searchsorted(db, lo[k_], "right"))
            if i1 > i0:
                Dopts, resid = [i0, i1 - 1], 0
            else:                                               # nearest ballast crossings either side
                Dopts = [d_ for d_ in range(i0 - 2, i0 + 2) if 0 <= d_ <= Dcap]
                resid = None
            for D in Dopts:
                res_ = int(db[D] - lo[k_]) if resid is None else 0
                if D > own + F and rho < 0:
                    continue
                # removed = own + m; ballast takes D; recipient takes own + m - D (needs rho if > 0)
                m_lo = max(0, D - own)
                m_hi = F if rho >= 0 else D - own
                if m_lo > m_hi:
                    continue
                # pick m by bisection on dA(m) against rem (dA decreases in m)
                bal_gain = (_even_split(D, Kb) * bw).sum() if Kb > 1 else D * float(bw[0])

                def dA_of(m_):
                    r_ = own + m_ - D
                    return (-(own * cw[h] + dropA(m_)) + bal_gain + (r_ * cw[rho] if rho >= 0 else 0.0))
                L_, H_ = m_lo, m_hi
                if dA_of(L_) > rem > dA_of(H_):
                    while H_ - L_ > 1:
                        M_ = (L_ + H_) // 2
                        if dA_of(M_) > rem:
                            L_ = M_
                        else:
                            H_ = M_
                for m_ in {L_, H_}:
                    if rho >= 0 and own + m_ - D > room[rho]:
                        continue
                    dA = dA_of(m_)
                    sc = abs(dA - rem) + (0.0 if dA >= rem else 0.25 * abs(dA - rem))
                    cands.append((abs(res_), sc, h, y, D, rho, m_, res_))
    if not cands:
        return None
    # exact first; then by closeness to the A target, preferring small residuals (easier to finish with pairs)
    cands.sort(key=lambda t: (t[0] > 0, t[1] + 0.02 * np.log1p(t[0])))
    for res_abs, sc, h, y, D, rho, m_, res_ in cands[:32]:
        caps = dn.copy(); caps[h] = 0
        if rho >= 0:
            caps[rho] = 0
        a_ = np.zeros(n, dtype=np.int64); left = m_
        for i in order:
            if left <= 0:
                break
            t_ = min(int(caps[i]), left); a_[i] = t_; left -= t_
        xn = x - a_
        xn[h] = y
        if rho >= 0:
            xn[rho] += int(x[h]) - y + m_ - D
        xbn = xb + _even_split(D, Kb)
        if res_ != 0:
            # carriers must change U by -(d_b - lo) = -res_ more: zero-sum pairs, h and rho frozen
            room2 = np.maximum(room - np.maximum(xn - x, 0), 0)
            room2[h] = 0
            fix = _offset(xn, s_car, c, -res_, room2, cw, -1, frozen=h)
            if fix is None:
                continue
            xn = fix[0]
        return xn, xbn
    return None


def fill(caps, need):
    """Max-min fair integer allocation of `need` units over cells with capacities caps (int64, BIG = open).
    Returns alloc (same shape) with alloc <= caps, sum(alloc) = min(need, sum(caps))."""
    caps = np.asarray(caps, dtype=np.int64)
    alloc = np.zeros_like(caps)
    need = int(need)
    live = np.flatnonzero(caps > 0)
    while need > 0 and live.size:
        share = need // live.size
        if share == 0:
            # hand out the last units one each to the largest-capacity cells
            order = live[np.argsort(-(caps[live] - alloc[live]), kind="stable")][:need]
            alloc[order] += 1
            need = 0
            break
        room = caps[live] - alloc[live]
        give = np.minimum(room, share)
        alloc[live] += give
        need -= int(give.sum())
        live = live[(caps[live] - alloc[live]) > 0]
    return alloc


# ------------------------------------------------------------------ the pass
def select(sel, l, ok):
    a = np.abs(l)
    if sel in ("all", "up", "dn"):
        base = ok.copy()
    else:
        q = 0.5 if sel == "q50" else 0.8
        thr = float(np.quantile(a[ok], q)) if ok.any() else np.inf
        base = ok & (a >= thr)
    if sel in ("up", "q80up"):
        base &= l > 0
    if sel in ("dn", "q80dn"):
        base &= l < 0
    if sel in ("q80up", "q80dn"):
        pass
    return base


def _levels(cw):
    """Carrier cells grouped by library class (the +1 remainder cells share their class's level)."""
    key = np.round(cw / cw.max(), 4)
    lv = np.unique(key)[::-1]                        # highest CPM-per-count (lowest library) first
    return [np.flatnonzero(key == v) for v in lv], np.array([cw[key == v].mean() for v in lv])


def _zs_step(x, xl, xh, room, cw, need, groups, gw):
    """Zero-sum carrier transfers of this gene between library classes, every cell staying inside its own
    control interval: counts move from high-library donors to low-library recipients (need > 0, A rises) or the
    reverse (need < 0). Block total, ballast and U are untouched. Returns the new x."""
    x = x.copy()
    rem = abs(need)
    up_cap = np.minimum(xh - x, room)
    dn_cap = x - xl
    pairs = []
    for a_ in range(len(groups)):
        for b_ in range(len(groups)):
            if a_ != b_:
                eff = gw[a_] - gw[b_]                # recipient a_, donor b_: A change per count
                if (need > 0 and eff > 0) or (need < 0 and eff < 0):
                    pairs.append((abs(eff), a_, b_))
    pairs.sort(reverse=True)
    for eff, ra, db in pairs:
        if rem <= 0:
            break
        # recipient group ra gains, donor group db loses: A moves by (gw[ra] - gw[db]) per count, the sign of need
        rc, dc = groups[ra], groups[db]
        want = int(np.floor(rem / eff + 0.5))
        tot = min(want, int(up_cap[rc].sum()), int(dn_cap[dc].sum()))
        if tot <= 0:
            continue
        ga = fill(up_cap[rc], tot)
        gd = fill(dn_cap[dc], tot)
        x[rc] += ga
        x[dc] -= gd
        up_cap[rc] -= ga
        dn_cap[dc] -= gd
        rem -= eff * tot
    return x


def _ballast_step(x, xl, xh, room, xb, cw, bw, need, groups, gw, rot=0):
    """Carriers move inside their intervals (lowest library first) and the ballast gives back / absorbs the same
    counts. Returns (x_new, xb_new) (xb_new may cross control values)."""
    up = need > 0
    cap = np.minimum(xh - x, room) if up else x - xl
    alloc = np.zeros(x.size, dtype=np.int64)
    rem = abs(need)
    ball_room = int(xb.sum()) if up else int(BIG)
    for cells, w in zip(groups, gw):
        eff = w - float(np.mean(bw))
        if eff <= 0 or rem <= 0:
            continue
        want = int(np.floor(rem / eff + 0.5))
        want = min(want, ball_room - int(alloc.sum()))
        if want <= 0:
            continue
        a_ = fill(cap[cells], want)
        alloc[cells] += a_
        rem -= eff * int(a_.sum())
    tot = int(alloc.sum())
    Kb = xb.size
    if Kb == 1:
        bsplit = np.array([tot], dtype=np.int64)
    else:
        wts = xb.astype(np.float64) if (up and xb.sum() > 0) else np.ones(Kb)
        raw = tot * wts / wts.sum()
        bsplit = np.floor(raw).astype(np.int64)
        r_ = tot - int(bsplit.sum())
        if r_:
            # remainders rotate across rows gene by gene (as ballast_rows does), so no row collects them all
            ordr = np.roll(np.arange(Kb), -(rot % Kb))
            bsplit[ordr[np.argsort(-(raw - bsplit)[ordr], kind="stable")][:r_]] += 1
        if up and np.any(bsplit > xb):
            spill = int((bsplit - xb)[bsplit > xb].sum()); bsplit = np.minimum(bsplit, xb)
            for q_ in np.argsort(-(xb - bsplit), kind="stable"):
                t_ = min(spill, int(xb[q_] - bsplit[q_])); bsplit[q_] += t_; spill -= t_
    sgn = 1 if up else -1
    return x + sgn * alloc, xb - sgn * bsplit, alloc, bsplit


MODES = ("pc", "zs", "bal", "ex", "zspc")


def postpass(block, n_car, cc, lfc_base, spec, *, sub_margin=0):
    """Apply the rank-preserving magnitude pass. Returns (new block, diagnostics).

    modes: pc   carriers per-cell rank-free, ballast absorbs (its control crossings are measured, not undone);
           zs   zero-sum carrier transfers between library classes only (ballast untouched; U exact per cell);
           bal  pc, then zero-sum carrier pairs that cross controls to undo the ballast's crossing exactly
                (gene rolled back if the offset cannot hit exactly);
           ex   zs, then bal on the remainder (a failed offset keeps the zs part);
           zspc zs, then pc on the remainder."""
    mode = spec["mode"]
    X = np.array(block, dtype=np.int64, copy=True)
    n = X.shape[0]
    keep = np.asarray(cc.keep, dtype=bool)
    gidx = np.flatnonzero(keep)
    sidx = np.flatnonzero(~keep)
    Lr = X.sum(1)
    s = Lr / 1e6                                     # the scorer's per-cell scale (float64 of an exact int)
    car = np.arange(n_car)
    bal = np.arange(n_car, n)
    ref = np.asarray(cc.ref_lin, dtype=np.float64)
    Xk = X[:, gidx]
    A = (Xk / s[:, None]).mean(0)
    l = np.asarray(lfc_base, dtype=np.float64)
    l_chk = np.log2((A + 1e-9) / (ref + 1e-9))
    fin = np.isfinite(l) & np.isfinite(l_chk)
    if np.max(np.abs(l[fin] - l_chk[fin])) > 1e-9:
        raise AssertionError("base lfc does not match the block")
    ok = fin & (A > 0) & (l != 0)
    k = float(spec["k"])
    touch = select(spec["sel"], l, ok)
    kvec = np.full(l.shape, k)
    if spec.get("k_up") is not None:
        kvec = np.where(l > 0, float(spec["k_up"]), k)
    A_t = (ref + 1e-9) * np.power(2.0, kvec * l) - 1e-9
    dA = np.where(touch, A_t - A, 0.0)
    # priority: largest |lfc| first (they carry the nmae and get the sub-gate budget first)
    order = np.flatnonzero(touch)[np.argsort(-np.abs(l[touch]), kind="stable")]
    sub_avail = X[np.ix_(car, sidx)].sum(1).astype(np.int64) - int(sub_margin)   # carrier sub-gate budget
    sub_delta = np.zeros(n_car, dtype=np.int64)       # net sub-gate change per carrier cell (+ = gained)
    # the ballast's sub-gate reservoir: carriers that lose kept counts are refilled with sub-gate counts from the
    # counts other carriers released, and the rest from the deepest ballast row (it holds ~0.5% of its library there)
    b_dep = bal[int(np.argmax(X[bal].sum(1)))]
    ball_sub = int(X[np.ix_(bal, sidx)].sum()) - int(sub_margin)     # every ballast row can refill carriers
    if bal.size > 1:
        # rows absorb evenly, so the poorest row bounds the total; keep headroom for settle_ballast's row balancing
        ball_sub = int(bal.size * (X[np.ix_(bal, sidx)].sum(1).min() - 64)) - int(sub_margin)
    lim_bsub = np.zeros(gidx.size, dtype=bool)
    cw = 1.0 / (n * s[car])                           # CPM of the 400-cell mean per count in each carrier cell
    bw = 1.0 / (n * s[bal])
    groups, gw = _levels(cw)
    G_ = gidx.size
    ach = np.zeros(G_); ach_zs = np.zeros(G_)
    rolled = np.zeros(G_, dtype=bool)
    dU_ball = np.zeros(G_, dtype=np.int64)
    dU_car = np.zeros(G_, dtype=np.int64)
    moved = np.zeros(G_, dtype=np.int64)
    pairs = np.zeros(G_, dtype=np.int64)
    dropdump = np.zeros(G_, dtype=bool)
    sign_guard = np.zeros(G_, dtype=bool)
    for j in order:
        need = dA[j]
        if need == 0:
            continue
        c = cc.cs[j]
        x0 = Xk[car, j].copy()
        xb0 = Xk[bal, j].copy()
        x, xb = x0, xb0
        room0 = np.maximum(sub_avail + sub_delta, 0)
        if mode in ("zs", "ex", "zspc"):
            xl, xh, *_ = cell_bounds(x, s[car], c)
            x = _zs_step(x, xl, xh, room0, cw, need, groups, gw)
            ach_zs[j] = ((x - x0) * cw).sum()
        x_zs = x.copy()
        rem = need - ((x - x0) * cw).sum()
        if mode in ("pc", "bal", "ex", "zspc") and abs(rem) > 0.5 * float(cw.min()):
            xl, xh, *_ = cell_bounds(x, s[car], c)
            room = np.maximum(room0 - np.maximum(x - x0, 0), 0)
            x1, xb1, alloc, bsplit = _ballast_step(x, xl, xh, room, xb, cw, bw, rem, groups, gw, rot=int(j))
            if mode in ("bal", "ex"):
                for _ in range(8):                   # the ballast must not land on a control tie (odd half-rank)
                    if not (_pos(xb1, s[bal], c)[1] > 0).any() or alloc.sum() == 0:
                        break
                    sg = 1 if rem > 0 else -1
                    i_ = int(np.argmax(alloc)); alloc[i_] -= 1
                    r_ = int(np.argmax(bsplit)); bsplit[r_] -= 1
                    x1 = x + sg * alloc; xb1 = xb - sg * bsplit
                d_b = int((ranks_of(xb1, s[bal], c) - ranks_of(xb, s[bal], c)).sum())
                if d_b != 0:
                    room1 = np.maximum(room0 - np.maximum(x1 - x0, 0), 0)
                    res = _offset(x1, s[car], c, -d_b, room1, cw, 1 if rem > 0 else -1)
                    if res is not None:
                        x1, pairs[j] = res
                    else:
                        dd = (_drop_dump(x, xl, xh, np.maximum(room0 - np.maximum(x - x0, 0), 0), xb, s[car], s[bal], c,
                                         cw, bw, rem) if rem < 0 else None)
                        if dd is not None:
                            x1, xb1 = dd
                            dropdump[j] = True
                        else:
                            rolled[j] = True
                            x1, xb1 = x, xb          # keep only the zs part (ex) / nothing (bal)
            x, xb = x1, xb1
        if np.array_equal(x, x0) and np.array_equal(xb, xb0):
            continue
        d_car = int((ranks_of(x, s[car], c) - ranks_of(x0, s[car], c)).sum())
        d_b = int((ranks_of(xb, s[bal], c) - ranks_of(xb0, s[bal], c)).sum())
        if mode in ("zs", "bal", "ex") and d_car + d_b != 0:
            raise AssertionError(f"mode {mode} moved U by {d_car + d_b}")
        if mode in ("pc", "zspc") and d_car != 0:
            raise AssertionError("a per-cell rank-free move crossed a control value")
        if int(x.sum() - x0.sum()) != -int(xb.sum() - xb0.sum()):
            raise AssertionError("gene total moved")
        dx = x - x0
        A_new = A[j] + (dx * cw).sum() + ((xb - xb0) * bw).sum()
        flip = np.sign(np.log2((A_new + 1e-9) / (ref[j] + 1e-9))) != np.sign(l[j])   # a magnitude pass never flips a sign
        over = int((sub_delta - dx).sum()) > ball_sub   # would draw more sub-gate counts than the ballast holds
        if flip or over:
            lim_bsub[j] |= over
            sign_guard[j] |= flip
            # fall back to the zero-sum part alone (it neither draws on the ballast nor moves U), if it is safe
            x, xb = x_zs, xb0
            dx = x - x0
            A_new = A[j] + (dx * cw).sum()
            if (not dx.any()) or np.sign(np.log2((A_new + 1e-9) / (ref[j] + 1e-9))) != np.sign(l[j]):
                continue
            d_car = int((ranks_of(x, s[car], c) - ranks_of(x0, s[car], c)).sum()); d_b = 0
            if d_car != 0:
                raise AssertionError("zero-sum fallback moved U")
        Xk[car, j] = x
        Xk[bal, j] = xb
        sub_delta -= dx
        if np.any(sub_avail + sub_delta < 0):
            raise AssertionError("sub-gate budget overdrawn")
        dU_car[j] = d_car; dU_ball[j] = d_b
        moved[j] = int(np.abs(dx).sum())
        ach[j] = (dx * cw).sum() + ((xb - xb0) * bw).sum()
    # ---- write the touched genes back, then realise the sub-gate compensation (libraries fixed)
    X[:, gidx] = Xk
    _subgate(X, car, bal, sidx, sub_delta, b_dep, Lr[bal])
    settle_ballast(X, bal, sidx, Lr[bal])
    if not np.array_equal(X.sum(1), Lr):
        raise AssertionError("a library moved")
    if not np.array_equal(X.sum(0), np.asarray(block).sum(0)):
        raise AssertionError("a per-gene block total moved")
    if X.min() < 0:
        raise AssertionError("negative count")
    t = np.flatnonzero(touch & (dA != 0))
    fr = np.clip(ach[t] / dA[t], -1, 2) if t.size else np.array([])
    w = np.abs(dA[t]) / max(np.abs(dA[t]).sum(), 1e-300) if t.size else np.array([])
    upm = dA[t] > 0
    diag = dict(n_touch=int(t.size), n_moved=int((moved > 0).sum()), counts_moved=int(moved.sum()),
                frac_mean=float(fr.mean()) if t.size else 0.0,
                frac_massw=float((fr * w).sum()) if t.size else 0.0,
                frac_massw_up=float((fr[upm] * w[upm]).sum() / max(w[upm].sum(), 1e-300)) if upm.any() else None,
                frac_massw_dn=float((fr[~upm] * w[~upm]).sum() / max(w[~upm].sum(), 1e-300)) if (~upm).any() else None,
                frac_massw_zs=float((np.clip(ach_zs[t] / dA[t], -1, 2) * w).sum()) if t.size else 0.0,
                frac_full=float(np.mean(fr >= 0.95)) if t.size else 0.0,
                frac_zero=float(np.mean(np.abs(fr) <= 0.05)) if t.size else 0.0,
                rolled=int(rolled.sum()), drop_dump=int(dropdump.sum()), sign_guard=int(sign_guard.sum()), lim_ballast_subgate=int(lim_bsub.sum()), ballast_sub_net=int(sub_delta.sum()),
                dU_ball_abs_sum=int(np.abs(dU_ball).sum()), dU_ball_abs_max=int(np.abs(dU_ball).max(initial=0)),
                n_dU_ball=int((dU_ball != 0).sum()), offset_pairs=int(pairs.sum()),
                dU_car_abs_sum=int(np.abs(dU_car).sum()), dU_net_abs_max=int(np.abs(dU_car + dU_ball).max(initial=0)),
                sub_used=int(-sub_delta[sub_delta < 0].sum()),
                sub_budget_exhausted_cells=int(((sub_avail + sub_delta) <= 0).sum()))
    diag["_fr_by_gene"] = (t, fr, l[t], dA[t])
    diag["_touched"] = moved > 0
    return X, diag


def _take(row, need):
    """Remove `need` counts from a non-negative integer row, largest holdings first, one count per gene per sweep."""
    take = np.zeros(row.size, dtype=np.int64)
    while need > 0:
        live = np.flatnonzero(row - take > 0)
        if live.size == 0:
            raise AssertionError("no counts left to take")
        o = live[np.argsort(-(row[live] - take[live]), kind="stable")][:need]
        take[o] += 1
        need -= o.size
    return take


def _subgate(X, car, bal, sidx, sub_delta, b_dep, L_bal):
    """Realise the carriers' sub-gate library compensation. Cells that gained kept counts give up sub-gate counts
    (largest holdings first) into a pool; cells that lost kept counts are refilled from that pool first, then from
    the ballast rows (the row with the most sub-gate counts left first). The deepest ballast row b_dep takes whatever
    is left in the pool, so every per-gene total is unchanged; settle_ballast then restores each row's library.
    With one ballast row this is exactly the single-row bookkeeping."""
    pool = np.zeros(sidx.size, dtype=np.int64)
    for i in np.flatnonzero(sub_delta < 0):
        t = _take(X[i, sidx], int(-sub_delta[i]))
        X[i, sidx] -= t
        pool += t
    Kb = bal.size
    drawn = np.zeros((Kb, sidx.size), dtype=np.int64)
    dk = X[bal].sum(1) - L_bal                           # each ballast row's kept-gene gain so far
    for i in np.flatnonzero(sub_delta > 0):
        d = int(sub_delta[i])
        from_pool = min(d, int(pool.sum()))
        g = _take(pool, from_pool) if from_pool else np.zeros(sidx.size, dtype=np.int64)
        pool -= g
        left = d - from_pool
        while left > 0:
            avail = X[np.ix_(bal, sidx)] - drawn
            owe = dk - drawn.sum(1)                      # rows that gained kept counts give sub-gate counts first
            r = int(np.argmax(owe)) if (owe > 0).any() and avail[int(np.argmax(owe))].sum() > 0 else int(np.argmax(avail.sum(1)))
            t_ = min(left, int(avail[r].sum()), int(owe[r]) if owe[r] > 0 else left)
            if t_ <= 0:
                raise AssertionError("ballast rows hold no sub-gate counts to give")
            gb = _take(avail[r], t_)
            drawn[r] += gb
            g = g + gb
            left -= t_
        X[i, sidx] += g
    for r in range(Kb):
        X[bal[r], sidx] -= drawn[r]
    X[b_dep, sidx] += pool
    if X[np.ix_(bal, sidx)].min() < 0:
        raise AssertionError("ballast sub-gate mirror went negative")
    return b_dep


def settle_ballast(X, bal, sidx, L_target):
    """Move sub-gate counts between ballast rows until every row's library equals L_target (same total)."""
    cur = X[bal].sum(1)
    if cur.sum() != L_target.sum():
        raise AssertionError("ballast mass drift")
    diff = cur - L_target                                # + = row holds too much
    for a in np.flatnonzero(diff > 0):
        for b in np.flatnonzero(diff < 0):
            t = int(min(diff[a], -diff[b]))
            if t <= 0:
                continue
            ra, rb = bal[a], bal[b]
            row = X[ra, sidx]
            take = np.zeros(sidx.size, dtype=np.int64)
            need = t
            while need > 0:
                live = np.flatnonzero(row - take > 0)
                if live.size == 0:
                    raise AssertionError("ballast row has no sub-gate counts to move")
                o = live[np.argsort(-(row[live] - take[live]), kind="stable")][:need]
                take[o] += 1
                need -= o.size
            X[ra, sidx] -= take
            X[rb, sidx] += take
            diff[a] -= t
            diff[b] += t
    if not np.array_equal(X[bal].sum(1), L_target):
        raise AssertionError("ballast rows not settled")


def summarise(diags):
    keys = [k for k in diags[0] if not k.startswith("_") and k != "target"]
    out = {}
    for k in keys:
        vals = [d[k] for d in diags if d.get(k) is not None]
        if vals and isinstance(vals[0], (int, float)):
            out[k + "_mean"] = float(np.mean(vals))
    return out


# ------------------------------------------------------------------ verification (scorer arithmetic)
def _r1_tie_kernel():
    from numba import njit, prange  # noqa: PLC0415

    @njit(parallel=True, cache=True)
    def r1_tie(ctrl_sorted, ctrl_tie, pert):
        G, m = ctrl_sorted.shape
        n = pert.shape[1]
        R1 = np.zeros(G, dtype=np.float64)
        TT = np.zeros(G, dtype=np.float64)
        for g in prange(G):
            ps = np.sort(pert[g])
            c = ctrl_sorted[g]
            r1 = 0.0
            tie = ctrl_tie[g]
            i = 0
            while i < n:
                j = i
                v = ps[i]
                while j + 1 < n and ps[j + 1] == v:
                    j += 1
                cnt_p = j - i + 1
                lo, hi = 0, m
                while lo < hi:
                    mid = (lo + hi) >> 1
                    if c[mid] < v:
                        lo = mid + 1
                    else:
                        hi = mid
                lt_c = lo
                lo, hi = 0, m
                while lo < hi:
                    mid = (lo + hi) >> 1
                    if c[mid] <= v:
                        lo = mid + 1
                    else:
                        hi = mid
                cnt_c = lo - lt_c
                tot = float(cnt_c + cnt_p)
                midrank = (lt_c + i) + (tot + 1.0) / 2.0
                r1 += cnt_p * midrank
                cc = float(cnt_c)
                tie += (tot * tot * tot - tot) - (cc * cc * cc - cc)
                i = j + 1
            R1[g] = r1
            TT[g] = tie
        return R1, TT
    return r1_tie


_R1T = None


def r1_tie(block, cc):
    """Rank sum R1 and tie sum per DE-universe gene, from exactly the scorer's per-cell values."""
    global _R1T
    if _R1T is None:
        _R1T = _r1_tie_kernel()
    x = np.asarray(block, dtype=np.float64)
    scale = x.sum(1) / 1e6
    scale[scale == 0] = 1.
    logged = np.ascontiguousarray(np.log1p((x / scale[:, None])[:, cc.keep]).T)
    return _R1T(cc.cs, np.ascontiguousarray(cc.tie), logged)


def verify(block, blk2, n_car, cc, s, p, q, l, s2, p2, q2, l2, touched=None, R0T0=None):
    R0, T0 = r1_tie(block, cc) if R0T0 is None else R0T0
    R1, T1 = r1_tie(blk2, cc)
    dR = R1 - R0
    ch = dR != 0
    fin = np.isfinite(l) & np.isfinite(l2)
    lp0 = np.log(np.maximum(p, 1e-300)); lp1 = np.log(np.maximum(p2, 1e-300))
    out = dict(sums_equal=bool(np.array_equal(s, s2)),
               libs_equal=bool(np.array_equal(np.asarray(block).sum(1), np.asarray(blk2).sum(1))),
               n_genes_U_changed=int(ch.sum()), U_abs_max=float(np.abs(dR).max(initial=0)),
               U_abs_sum=float(np.abs(dR).sum()),
               n_genes_tie_changed=int((T1 != T0).sum()),
               p_logmax=float(np.max(np.abs(lp1 - lp0))), p_equal=bool(np.array_equal(p, p2)),
               q_equal=bool(np.array_equal(q, q2)),
               call_flips=int(((q < .05) != (q2 < .05)).sum()),
               sign_flips=int((np.sign(l[fin]) != np.sign(l2[fin])).sum()),
               lfc_changed=int((l2 != l).sum()))
    if touched is not None:
        out["n_touched_U_changed"] = int((ch & touched).sum())
        out["n_untouched_U_changed"] = int((ch & ~touched).sum())
    return out

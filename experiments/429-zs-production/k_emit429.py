"""426 copy of 416's k_emit416.py + one optional `mde_add` key in the prior spec (MDE-structured magnitude, A-channel
only; see the 426 block). 416's doc: 416 copy of 400's k_emit400.py (389 K-class emitter) + one optional `gate` key in
the prior spec (see E1). Original doc:

Production `rc_hybrid_block` (src/vcc2026/rc_emit.py) welds every emitted cell to one
library total T = S // n.  With a single per-cell total the per-gene ARITHMETIC mean
of per-cell CPM — the quantity the scorer's DE table (LFC magnitude and sign) reads —
is a deterministic function of the per-gene panel total, i.e. the DE channel is welded
to the bulk (pds/mse) channel.  This module unwelds it: cells are assigned to K
library classes with levels L_k (sum still S), so at fixed per-gene total M the
arithmetic mean A = (1e6/n) * sum_cells c_i / L_{k(i)} is steerable by placement
(counts in low-library cells raise A, high-library cells lower it) while per-gene
totals — and hence pds and the profile half of mse — stay bit-identical.

Where the equal-T assumption binds, and how it is generalised per class
------------------------------------------------------------------------
`phi_tables` maps a COUNT level c to (#controls below, #controls equal) through the
value v = log1p(1e6*c/T).  With one T, count order == CPM order, so one table per
gene suffices.  With K libraries the same count means a different CPM value in each
class, so each class gets its OWN Phi table built at its own level
(v = log1p(1e6*c/L_k)); within a class count<->CPM monotonicity still holds and the
existing `design_gene` / `histogram_of` family scan runs unchanged against that
class's tables.  The pooled Wilcoxon z over classes needs care at exactly one point:
two (count, class) pairs can share one CPM value (e.g. c=2 at L and c=4 at 2L), and
the rank-sum groups them.  `hist_z_k` groups pairs by exact CPM value (float `==` on
the identically-computed `log1p` values — the same floats the scorer's kernel
compares, so the grouping is scorer-consistent by construction) and otherwise mirrors
`hist_z`'s scanpy-style tie-corrected computation.  An equivalence check
(`check_kernel_equivalence`, run once per panel) asserts `hist_z_k == wilcoxon_z` on
materialised per-cell CPM rows to < 1e-9 before any changed arm is interpreted.

The arithmetic mean depends ONLY on the per-class totals M_k
(A = (1e6/n) * sum_k M_k / L_k), never on the within-class arrangement, so steering
is exact up to 1-count granularity and is asserted as such.  `repair_rows` only ever
moves sub-gate (non-universe) genes, so it cannot disturb a steered gene's M_k.

K = 1 is degenerate (one class => A pinned by M, no steering freedom) and takes the
welded path verbatim: incumbent totals, `EME.design` / `EME.place`, identical RNG
draw order — bit-identical to `rc_hybrid_block` by construction.
"""
from __future__ import annotations

import numpy as np
from numba import njit, prange

from vcc2026 import rc_emit
from vcc2026.rc_emit import (
    CALL_SET_OFFSET,
    JITTER_SCALE,
    EME,
    assign_cells,
    de_from_block,
    design_gene,
    histogram_of,
    phi_tables,
    repair_rows,
    wilcoxon_z,
)

LAST_PRIOR_REC = []   # 391: per-call prior diagnostics, popped by run391
PMAX = 16          # max pooled (level, class) pairs: <= 4 classes x 4 levels
KMAX = 4
NLEV = 4


@njit(cache=True, inline="always")
def _pgroup(Lv, Mu, Cv, Ktot, Lk, lt_all, eq_all, og, gv_lt, gv_cc, gv_mult):
    """Group pooled (level, mult, class) pairs by exact CPM value.

    Sorts pair indices by v = log1p(1e6*c/L_k), then groups float-equal values.
    Fills gv_lt (controls strictly below), gv_cc (controls equal), gv_mult (pooled
    pert multiplicity) and returns the number of groups.  lt/eq for a group are read
    from any member's class table: equal CPM value => equal control counts.
    """
    order = np.empty(PMAX, np.int64)
    vv = np.empty(PMAX, np.float64)
    for t in range(Ktot):
        vv[t] = np.log1p(1e6 * float(Lv[t]) / float(Lk[Cv[t]]))
        order[t] = t
    for a in range(1, Ktot):          # insertion sort; Ktot <= 16
        t = order[a]
        vt = vv[t]
        b = a - 1
        while b >= 0 and vv[order[b]] > vt:
            order[b + 1] = order[b]
            b -= 1
        order[b + 1] = t
    ng = 0
    t = 0
    while t < Ktot:
        u = t + 1
        while u < Ktot and vv[order[u]] == vv[order[t]]:
            u += 1
        mult = 0
        for s in range(t, u):
            mult += Mu[order[s]]
        c0 = Lv[order[t]]
        k0 = Cv[order[t]]
        gv_lt[ng] = float(lt_all[og[k0] + c0])
        gv_cc[ng] = float(eq_all[og[k0] + c0])
        gv_mult[ng] = mult
        ng += 1
        t = u
    return ng


@njit(cache=True)
def rebalance_within_class(B, resid, cell_class, gated_cols, max_moves):
    """Move single GATED counts between cells of the SAME class until every
    residual is zero.  Per-class gene totals (hence steering targets) are
    preserved exactly; only the within-class arrangement moves.  Backstop for
    `repair_rows`, which stalls when an over-full cell holds no sub-gate counts
    (extreme steering fills low-library cells with gated counts only).  Each
    move strictly decreases sum|resid|, so the loop always terminates."""
    n = B.shape[0]
    moves = 0
    while moves < max_moves:
        # most over-full cell among classes that also hold an under-full cell
        src = -1
        smin = 0
        for i in range(n):
            if resid[i] < smin:
                k = cell_class[i]
                has = False
                for j in range(n):
                    if cell_class[j] == k and resid[j] > 0:
                        has = True
                        break
                if has:
                    smin = resid[i]
                    src = i
        if src < 0:
            break
        k = cell_class[src]
        dst = -1
        dmax = 0
        for i in range(n):
            if cell_class[i] == k and resid[i] > dmax:
                dmax = resid[i]
                dst = i
        if dst < 0:
            break
        Gg = gated_cols.shape[0]
        found = -1
        for jj in range(Gg):
            g = gated_cols[(moves + jj) % Gg]
            if B[src, g] > 0:
                found = g
                break
        if found < 0:
            break
        B[src, found] -= 1
        B[dst, found] += 1
        resid[src] += 1
        resid[dst] -= 1
        moves += 1
    unresolved = 0
    for i in range(n):
        unresolved += abs(resid[i])
    return moves, unresolved


@njit(cache=True)
def hist_z_k(Lv, Mu, Cv, Ktot, Lk, lt_all, eq_all, og, ctrl_tie, n, m):
    """Exact scanpy-style tie-corrected z for a pooled multi-class histogram."""
    gv_lt = np.empty(PMAX, np.float64)
    gv_cc = np.empty(PMAX, np.float64)
    gv_mult = np.empty(PMAX, np.int64)
    ng = _pgroup(Lv, Mu, Cv, Ktot, Lk, lt_all, eq_all, og, gv_lt, gv_cc, gv_mult)
    N = n + m
    R1 = 0.0
    tie = ctrl_tie
    below = 0
    for g in range(ng):
        cc = gv_cc[g]
        tot = cc + gv_mult[g]
        midrank = (gv_lt[g] + below) + (tot + 1.0) / 2.0
        R1 += gv_mult[g] * midrank
        tie += (tot * tot * tot - tot) - (cc * cc * cc - cc)
        below += gv_mult[g]
    tc = 1.0 - tie / (float(N) ** 3 - float(N))
    sd = np.sqrt(tc * n * m * (N + 1.0) / 12.0)
    if sd <= 0.0:
        return 0.0
    return (R1 - n * (N + 1.0) / 2.0) / sd


@njit(parallel=True, cache=True)
def design_all_k(Ms, nks, Lk, cmaxs, og, lt_all, eq_all, Phi_all, ctrl_ties, m,
                 ztargets, sd_pool, Kcls):
    """Per-gene multi-class histogram choice for signed z targets (mode "abs").

    Mirrors `design_all`'s mode 1: per class the two one-parameter families are
    scanned (reused `design_gene`) for the Fbar targets 0.5 +- |z*|*sd_pool; over the
    2^K sign patterns the pooled realisation (exact `hist_z_k`) whose |z| is closest
    to |z*| wins, ties broken toward sign(z*).  Returns per gene per class the
    winning (K, levels, mults) plus the pooled realised/even z.
    """
    G = Ms.shape[0]
    n = 0
    for k in range(Kcls):
        n += nks[k]
    KK = np.zeros((G, KMAX), np.int64)
    LL = np.zeros((G, KMAX, NLEV), np.int64)
    MM = np.zeros((G, KMAX, NLEV), np.int64)
    zr = np.zeros(G, np.float64)
    ze = np.zeros(G, np.float64)
    for g in prange(G):
        zt = ztargets[g]
        sd = sd_pool[g]
        fam = np.zeros((KMAX, 2), np.int64)
        par = np.zeros((KMAX, 2), np.int64)
        bsk = np.zeros(KMAX, np.int64)
        targets = np.empty(2, np.float64)
        targets[0] = 0.5 + abs(zt) * sd
        targets[1] = 0.5 - abs(zt) * sd
        for k in range(Kcls):
            fo = np.empty(2, np.int64)
            po = np.empty(2, np.int64)
            _fe, _fmn, _fmx, bs = design_gene(
                int(Ms[g, k]), int(nks[k]), int(cmaxs[g, k]), Phi_all,
                int(og[g, k]), lt_all, eq_all, ctrl_ties[g], m,
                targets, 2, fo, po)
            bsk[k] = bs
            for t in range(2):
                fam[k, t] = fo[t]
                par[k, t] = po[t]
        # even pooled z (the "natural" allocation)
        Lv = np.empty(PMAX, np.int64)
        Mu = np.empty(PMAX, np.int64)
        Cv = np.empty(PMAX, np.int64)
        Ktot = 0
        lv = np.zeros(NLEV, np.int64)
        mu = np.zeros(NLEV, np.int64)
        for k in range(Kcls):
            Ke = histogram_of(int(Ms[g, k]), int(nks[k]), int(cmaxs[g, k]),
                              0, 0, -1, lv, mu)
            for q in range(Ke):
                Lv[Ktot] = lv[q]
                Mu[Ktot] = mu[q]
                Cv[Ktot] = k
                Ktot += 1
        ze[g] = hist_z_k(Lv, Mu, Cv, Ktot, Lk, lt_all, eq_all, og[g],
                         ctrl_ties[g], n, m)
        bestd = 1e300
        bestz = 0.0
        bestpat = 0
        npat = 1 << Kcls
        for pat in range(npat):
            Ktot = 0
            for k in range(Kcls):
                t = (pat >> k) & 1
                Kh = histogram_of(int(Ms[g, k]), int(nks[k]), int(cmaxs[g, k]),
                                  int(fam[k, t]), int(par[k, t]), int(bsk[k]),
                                  lv, mu)
                for q in range(Kh):
                    Lv[Ktot] = lv[q]
                    Mu[Ktot] = mu[q]
                    Cv[Ktot] = k
                    Ktot += 1
            z = hist_z_k(Lv, Mu, Cv, Ktot, Lk, lt_all, eq_all, og[g],
                         ctrl_ties[g], n, m)
            d = abs(abs(z) - abs(zt))
            if (z > 0) != (zt > 0):
                d += 1e-9                  # tie-break toward sign(z*), as design_all
            if d < bestd:
                bestd = d
                bestz = z
                bestpat = pat
        zr[g] = bestz
        for k in range(Kcls):
            t = (bestpat >> k) & 1
            Kh = histogram_of(int(Ms[g, k]), int(nks[k]), int(cmaxs[g, k]),
                              int(fam[k, t]), int(par[k, t]), int(bsk[k]),
                              lv, mu)
            KK[g, k] = Kh
            for q in range(NLEV):
                LL[g, k, q] = lv[q]
                MM[g, k, q] = mu[q]
    return KK, LL, MM, zr, ze


# --------------------------------------------------------------------------
# Python drivers (all numpy; njit kernels above do the per-gene work)

def class_cells(n, classes, S):
    """Assign n cells to K library classes with levels scaled to sum exactly S.

    `classes`: sequence of (fraction, raw level).  Cell counts by largest
    remainder; levels scaled by S / sum(n_k q_k) and rounded; the integer residue
    is fixed on the first cells (the incumbent T+1 remainder convention).
    Returns (cell_class (n,), base levels Lk (K,), actual totals T_cells (n,),
    counts nk (K,)).
    """
    fracs = np.array([f for f, _ in classes], dtype=np.float64)
    raw = np.array([q for _, q in classes], dtype=np.float64)
    K = len(classes)
    if K > KMAX:
        raise ValueError(f"at most {KMAX} classes")
    if np.any(raw <= 0) or abs(fracs.sum() - 1.0) > 1e-9:
        raise ValueError("classes need positive levels and fractions summing to 1")
    base = np.floor(fracs * n).astype(np.int64)
    rem = int(n - base.sum())
    order = np.argsort(-(fracs * n - base), kind="stable")
    nk = base.copy()
    nk[order[:rem]] += 1
    scale = S / float((nk * raw).sum())
    Lk = np.maximum(1, np.round(raw * scale)).astype(np.int64)
    cell_class = np.repeat(np.arange(K, dtype=np.int64), nk)
    T_cells = Lk[cell_class].copy()
    diff = int(S - T_cells.sum())
    if diff > 0:
        T_cells[:diff] += 1
    elif diff < 0:
        if (-diff) > n or np.any(T_cells[:-diff if diff else n] < 2):
            raise AssertionError("remainder fix would empty a cell")
        T_cells[: -diff] -= 1
    if T_cells.sum() != S or np.any(T_cells < 1):
        raise AssertionError("class levels do not sum to S")
    return cell_class, Lk, T_cells, nk


def split_massprop(Ms, w):
    """Unsteered split: M_k proportional to class library mass (=> A == bulk)."""
    Ms = np.asarray(Ms, dtype=np.int64)
    w = np.asarray(w, dtype=np.float64)
    base = np.floor(Ms[:, None] * w[None, :]).astype(np.int64)
    rem = Ms - base.sum(1)
    frac = Ms[:, None] * w[None, :] - base
    rank = np.argsort(np.argsort(-frac, axis=1, kind="stable"), axis=1,
                      kind="stable")
    return base + (rank < rem[:, None]).astype(np.int64)


def feasible_dose(Ms, ref_lin, S, s_nom, mask, Lk, nk, caps, mass_k,
                  margin=6000, iters=8, perm=None):
    """Per-target feasible steering dose s* <= s_nom (class-level pigeonhole).

    Per-gene arithmetic targets are individually feasible, but jointly the steered
    demand can exceed a class's mass (sub-gate mass occupies the rest; even fully
    drained, gated_k <= mass_k is hard).  Bisection on s in [1, s_nom] keeps the
    largest feasible dose with `margin` counts of fragmentation slack per class
    (global repair moves sub-gate mass across classes exactly as in the unsteered
    arm).  Returns (s_star, Mk, A_t, A_ach, c_hi, c_lo) with A_t the full
    (possibly deranged) target array (NaN = unsteered) and A_ach the achieved
    means over the steered subset.  s = 1 (design dose, split for the bulk
    target) is feasible by the same margin whenever the unsteered split is; if
    not, falls back to the plain mass-proportional split (the S1 state, A_t all
    NaN, no steering asserts).
    """
    Ms = np.asarray(Ms, dtype=np.int64)
    Kcls = Lk.shape[0]
    w = (nk * Lk).astype(np.float64)
    w /= w.sum()
    base = split_massprop(Ms, w)
    bulk = 1e6 * Ms.astype(np.float64) / S
    with np.errstate(divide="ignore", invalid="ignore"):
        lfc_d = np.log2(np.maximum(bulk, 1e-12) / np.maximum(ref_lin, 1e-12))
    pm = None
    if perm is not None:
        pm = np.random.default_rng(int(perm)).permutation(Ms.shape[0])

    def demand(s):
        with np.errstate(divide="ignore", invalid="ignore"):
            A_t = ref_lin * np.power(2.0, s * lfc_d)
        A_t[Ms == 0] = 0.0
        A_t = np.where(mask, A_t, np.nan)
        if pm is not None:
            # shuffled-class control: derange resolved targets across genes;
            # feasibility is judged on the deranged (actually placed) demand.
            A_t = np.ascontiguousarray(A_t[pm])
        Mk = base.copy()
        st = np.isfinite(A_t)
        if st.any():
            Mk_s, A_ach, c_hi, c_lo = split_for_A(
                Ms[st], A_t[st], Lk, nk, caps[st])
            Mk[st] = Mk_s
        else:
            A_ach, c_hi, c_lo = np.array([]), np.array([]), np.array([])
        over = Mk.sum(0) - (mass_k - margin)
        return bool(np.all(over <= 0)), Mk, A_t, A_ach, c_hi, c_lo

    ok1, Mk1, A_t1, A_ach1, chi1, clo1 = demand(1.0)
    if not ok1:
        # even the design dose strains the margin: plain massprop (S1 state).
        return 1.0, base, np.full(Ms.shape, np.nan), np.array([]), np.array([])
    lo, best = 1.0, (1.0, Mk1, A_t1, A_ach1, chi1, clo1)
    hi = float(s_nom)
    okh, Mkh, A_th, A_achh, chih, cloh = demand(hi)
    if okh:
        return hi, Mkh, A_th, A_achh, chih, cloh
    for _ in range(iters):
        mid = (lo + hi) / 2
        okm, Mkm, A_tm, A_achm, chim, clom = demand(mid)
        if okm:
            lo, best = mid, (mid, Mkm, A_tm, A_achm, chim, clom)
        else:
            hi = mid
    s_star = best[0]
    # verify (non-monotone demand safety): step down until feasible.
    for _ in range(60):
        okv, Mkv, A_tv, A_achv, chiv, clov = demand(s_star)
        if okv:
            return s_star, Mkv, A_tv, A_achv, chiv, clov
        s_star = max(1.0, s_star - 0.005)
    return 1.0, base, np.full(Ms.shape, np.nan), np.array([]), np.array([])


def arith_of_split(Mk, Lk, n):
    """Delivered arithmetic-mean CPM for a per-class total split (exact)."""
    C = 1e6 / n
    return C * (Mk.astype(np.float64) / Lk[None, :]).sum(1)


def split_for_A(Ms, A_t, Lk, nk, caps, A_tol=1e-12):
    """Greedy per-class split hitting arithmetic-mean targets (vectorised bulk
    moves + single-count polish).  Returns (Mk, A_ach, clipped_hi, clipped_lo).

    `caps`: per-gene per-class capacity (n_k * L_k minus sub-gate mass).
    Unreachable targets clip to the feasible extreme and are flagged.
    """
    Ms = np.asarray(Ms, dtype=np.int64)
    A_t = np.asarray(A_t, dtype=np.float64)
    Lk = np.asarray(Lk, dtype=np.int64)
    G, K = Ms.shape[0], Lk.shape[0]
    n = int(nk.sum())
    C = 1e6 / n
    inv = 1.0 / Lk.astype(np.float64)
    w = (nk * Lk).astype(np.float64)
    w /= w.sum()
    Mk = split_massprop(Ms, w)
    Mk = np.minimum(Mk, caps)
    # capacity clipping can drop mass; restore from headroom, lowest-L first
    short = Ms - Mk.sum(1)
    lo_first = np.argsort(Lk, kind="stable")
    for k in lo_first:
        if np.all(short <= 0):
            break
        room = np.minimum(short, caps[:, k] - Mk[:, k])
        room = np.maximum(room, 0)
        Mk[:, k] += room
        short -= room
    if np.any(short != 0):
        raise AssertionError("class capacities cannot hold the totals")
    A = arith_of_split(Mk, Lk, n)
    clipped_hi = np.zeros(G, dtype=bool)
    clipped_lo = np.zeros(G, dtype=bool)
    asc = np.argsort(Lk, kind="stable")       # asc[0] = sink for raising
    desc = asc[::-1]                          # desc[0] = sink for lowering
    for gi in range(G):
        M = int(Ms[gi])
        if M == 0:
            continue
        t = float(A_t[gi])
        # feasible extremes (cap-aware fill cascades)
        cap = caps[gi].astype(np.int64)
        hi = np.zeros(K, dtype=np.int64)
        r = M
        for k in asc:
            x = min(r, int(cap[k]))
            hi[k] = x
            r -= x
        Amax = C * float((hi * inv).sum())
        lob = np.zeros(K, dtype=np.int64)
        r = M
        for k in desc:
            x = min(r, int(cap[k]))
            lob[k] = x
            r -= x
        Amin = C * float((lob * inv).sum())
        if t >= Amax:
            Mk[gi] = hi
            A[gi] = Amax
            if t > Amax + A_tol:
                clipped_hi[gi] = True
            continue
        if t <= Amin:
            Mk[gi] = lob
            A[gi] = Amin
            if t < Amin - A_tol:
                clipped_lo[gi] = True
            continue
        row = Mk[gi].astype(np.int64).copy()
        cur = float(arith_of_split(row[None, :], Lk, n)[0])
        if cur < t:                            # raise: sink = lowest-L class
            sink = int(asc[0])
            for s in asc[1:][::-1]:            # highest-L source first
                s = int(s)
                d = C * (inv[sink] - inv[s])
                if d <= 0 or row[s] == 0:
                    continue
                need = (t - cur) / d
                x = min(int(row[s]), int(cap[sink] - row[sink]),
                        int(np.floor(need + 1e-9)))
                if x > 0:
                    row[s] -= x
                    row[sink] += x
                    cur += x * d
                if cur >= t - A_tol:
                    break
            while cur < t - A_tol:             # single-count polish (<= K steps)
                moved = False
                for s in asc[1:][::-1]:
                    s = int(s)
                    if row[s] > 0 and row[sink] < cap[sink]:
                        d = C * (inv[sink] - inv[s])
                        if cur + d <= t + d / 2:
                            row[s] -= 1
                            row[sink] += 1
                            cur += d
                            moved = True
                            break
                if not moved:
                    break
        else:                                  # lower: sink = highest-L class
            sink = int(desc[0])
            for s in desc[1:][::-1]:
                s = int(s)
                d = C * (inv[s] - inv[sink])
                if d <= 0 or row[s] == 0:
                    continue
                need = (cur - t) / d
                x = min(int(row[s]), int(cap[sink] - row[sink]),
                        int(np.floor(need + 1e-9)))
                if x > 0:
                    row[s] -= x
                    row[sink] += x
                    cur -= x * d
                if cur <= t + A_tol:
                    break
            while cur > t + A_tol:
                moved = False
                for s in desc[1:][::-1]:
                    s = int(s)
                    if row[s] > 0 and row[sink] < cap[sink]:
                        d = C * (inv[s] - inv[sink])
                        if cur - d >= t - d / 2:
                            row[s] -= 1
                            row[sink] += 1
                            cur -= d
                            moved = True
                            break
                if not moved:
                    break
        Mk[gi] = row
        A[gi] = cur
    return Mk, A, clipped_hi, clipped_lo


def arith_range(Ms, Lk, nk, caps):
    """Feasible (A_min, A_max) arithmetic-mean CPM per gene (cap-aware fill)."""
    Ms = np.asarray(Ms, dtype=np.int64)
    G, K = Ms.shape[0], Lk.shape[0]
    n = int(nk.sum())
    C = 1e6 / n
    inv = 1.0 / Lk.astype(np.float64)
    asc = np.argsort(Lk, kind="stable")
    desc = asc[::-1]
    Amin = np.empty(G)
    Amax = np.empty(G)
    for gi in range(Ms.shape[0]):
        M = int(Ms[gi])
        cap = caps[gi].astype(np.int64)
        hi = np.zeros(K, dtype=np.int64)
        r = M
        for k in asc:
            x = min(r, int(cap[k]))
            hi[k] = x
            r -= x
        Amax[gi] = C * float((hi * inv).sum())
        lob = np.zeros(K, dtype=np.int64)
        r = M
        for k in desc:
            x = min(r, int(cap[k]))
            lob[k] = x
            r -= x
        Amin[gi] = C * float((lob * inv).sum())
    return Amin, Amax


def scale_targets(Ms, ref_lin, S, s, mask):
    """Arithmetic-CPM targets at s x the design LFC (log2 space).

    A_design = bulk CPM = 1e6*M/S; A_t = ref_lin * 2^(s * log2(A_design/ref_lin))
    on masked genes.  M = 0 -> A_t = 0.
    """
    Ms = np.asarray(Ms, dtype=np.float64)
    ref_lin = np.asarray(ref_lin, dtype=np.float64)
    bulk = 1e6 * Ms / S
    with np.errstate(divide="ignore", invalid="ignore"):
        lfc_d = np.log2(np.maximum(bulk, 1e-12) / np.maximum(ref_lin, 1e-12))
    A_t = ref_lin * np.power(2.0, s * lfc_d)
    A_t[Ms == 0] = 0.0
    return np.where(np.asarray(mask, dtype=bool), A_t, np.nan)


def prescribed_mask(key, control, call_set):
    """The CallSet membership recipe of `rc_hybrid_block` as a gated bool mask."""
    key_gated = np.asarray(key, dtype=np.float64)[control.keep]
    key_abs = np.abs(key_gated)
    n_keep = key_abs.size
    head = np.zeros(n_keep, dtype=bool)
    head[np.argsort(-key_abs, kind="stable")[: int(call_set.n_head)]] = True
    cpm_rank = np.where(head, -np.inf,
                        np.asarray(control.ref_lin, dtype=np.float64))
    called = head.copy()
    called[np.argsort(-cpm_rank, kind="stable")[: int(call_set.n_cpm)]] = True
    if int(called.sum()) != call_set.size:
        raise AssertionError("prescribed call set did not come out at its size")
    return called


def build_tables(cs, m, Mk, Lk):
    """Per-class Phi tables concatenated to one flat layout.

    Returns (cmaxs (G,K), og (G,K), lt_all, eq_all, Phi_all).  Class k's rows live
    at og[g,k] = coff[k] + offs_k[g].
    """
    G, K = Mk.shape
    offs_list, lt_list, eq_list, phi_list = [], [], [], []
    cmaxs = np.empty((G, K), dtype=np.int64)
    og = np.empty((G, K), dtype=np.int64)
    base = 0
    for k in range(K):
        cmax = np.minimum(Mk[:, k], Lk[k]).astype(np.int64)
        cmaxs[:, k] = cmax
        offs, lt, eq = phi_tables(np.ascontiguousarray(cs),
                                  float(Lk[k]),
                                  np.ascontiguousarray(cmax))
        og[:, k] = base + offs[:-1]
        base += int(offs[-1])
        offs_list.append(offs)
        lt_list.append(lt)
        eq_list.append(eq)
        phi_list.append((lt.astype(np.float64) + 0.5 * eq.astype(np.float64)) / m)
    lt_all = np.ascontiguousarray(np.concatenate(lt_list)).astype(np.int32)
    eq_all = np.ascontiguousarray(np.concatenate(eq_list)).astype(np.int32)
    Phi_all = np.ascontiguousarray(np.concatenate(phi_list))
    return cmaxs, og, lt_all, eq_all, Phi_all


def check_kernel_equivalence(cs, tie, seed=389, n_rep=3):
    """Assert hist_z_k == wilcoxon_z (< 1e-9) on synthetic multi-class histograms.

    Uses real control rows (so ties are realistic); per-cell CPM rows are
    materialised explicitly for the reference kernel.  Returns the max abs diff.
    """
    rng = np.random.default_rng(seed)
    G = 24
    m = cs.shape[1]
    css = np.ascontiguousarray(cs[:G])
    ties = np.ascontiguousarray(tie[:G])
    Lk = np.array([13986, 20191, 27234], dtype=np.int64)
    nks = np.array([5, 6, 5], dtype=np.int64)
    n = int(nks.sum())
    worst = 0.0
    for _ in range(n_rep):
        Ms = rng.integers(0, 60, size=(G, 3)).astype(np.int64)
        cmaxs, og, lt_all, eq_all, Phi_all = build_tables(css, m, Ms, Lk)
        zt = rng.normal(0, 3, size=G)
        sd = sd_pool_all(Ms, nks, Lk, cmaxs, og, lt_all, eq_all, ties, m, 3)
        KK, LL, MM, zr, _ze = design_all_k(Ms, nks, Lk, cmaxs, og, lt_all,
                                          eq_all, Phi_all, ties, m, zt, sd, 3)
        pert = np.empty((G, n))
        for g in range(G):
            vals = []
            for k in range(3):
                for q in range(int(KK[g, k])):
                    c = int(LL[g, k, q])
                    vals += [np.log1p(1e6 * c / float(Lk[k]))] * int(MM[g, k, q])
            assert len(vals) == n, (len(vals), n)
            pert[g] = vals
        zref = wilcoxon_z(css, ties, np.ascontiguousarray(pert))
        d = float(np.abs(zr - zref).max())
        worst = max(worst, d)
        if not np.all(np.asarray(Ms).sum(1) == np.array(
                [sum(int(MM[g, k, q]) * int(LL[g, k, q])
                     for k in range(3) for q in range(int(KK[g, k])))
                 for g in range(G)])):
            raise AssertionError("histogram totals do not match Ms")
    if not worst < 1e-9:
        raise AssertionError(f"hist_z_k != wilcoxon_z: max diff {worst}")
    return worst


def place_k(dsg, gated_idx, sub_block, T_cells, cell_class, Lk, nk, pass2_rng):
    """Assign per-class designed histograms to class cells, then global repair.

    Per class the call sequence mirrors `EME.place` (capacity, order by max level,
    `assign_cells` with fresh RNG keys); the residual repair is the shared
    `repair_rows` on the sub-gate block.  Returns (B, diagnostics).
    """
    K = int(len(nk))
    n = int(T_cells.shape[0])
    G_full = sub_block.shape[1]
    KK, LL, MM = dsg["K"], dsg["levels"], dsg["mults"]
    B = np.ascontiguousarray(sub_block).astype(np.int64).copy()
    forced = 0
    for k in range(K):
        rows = np.flatnonzero(cell_class == k)
        cap = T_cells[rows].astype(np.int64) - sub_block[rows].sum(1)
        Kk = np.ascontiguousarray(KK[:, k])
        Lk_ = np.ascontiguousarray(LL[:, k, :])
        Muk = np.ascontiguousarray(MM[:, k, :])
        maxlvl = np.array([Lk_[g, Kk[g] - 1] if Kk[g] > 0 else 0
                           for g in range(Kk.size)])
        order = np.argsort(-maxlvl, kind="stable").astype(np.int64)
        Bg = np.zeros((rows.size, Kk.size), np.int32)
        keys = pass2_rng.random((order.size, rows.size))
        forced += int(assign_cells(Kk, Lk_, Muk, order, cap, Bg, keys))
        B[np.ix_(rows, gated_idx)] += Bg.astype(np.int64)
    resid = T_cells.astype(np.int64) - B.sum(1)
    max_abs_resid = int(np.abs(resid).max())
    sub_idx = np.setdiff1d(np.arange(G_full), gated_idx)
    Bs = np.ascontiguousarray(B[:, sub_idx])
    perm = pass2_rng.permutation(sub_idx.size).astype(np.int64)
    moves, unresolved = repair_rows(Bs, resid, perm)
    B[:, sub_idx] = Bs
    rebalanced = 0
    if unresolved != 0:
        # sub-gate repair stalls when an over-full cell holds no sub-gate counts;
        # move gated counts within the class instead (M_k preserved exactly).
        max_moves = int(2 * np.abs(resid).sum() + n)
        rebalanced, unresolved = rebalance_within_class(
            B, np.ascontiguousarray(resid),
            np.ascontiguousarray(cell_class),
            np.ascontiguousarray(gated_idx), max_moves)
    return B, dict(forced=int(forced), max_abs_resid=max_abs_resid,
                   moves=int(moves), unresolved=int(unresolved),
                   rebalanced=int(rebalanced))


def rc_hybrid_block_k(
    profile,
    lib_sizes,
    dispersion,
    rng,
    *,
    spread,
    overdispersion,
    control,
    key,
    jitter_seed,
    call_set_first: bool = True,
    call_set=None,
    classes=None,
    lfc_target=None,
    pass2_rng=None,
    max_counts: int = 1_000_000,
    emit_counts=None,
    order_key=None,
):
    """K-class-library variant of `rc_emit.rc_hybrid_block`.

    Same pass 1 (incumbent block under the shipped quasi law, same rng stream) and
    same sub-gate block, so per-gene panel totals are bit-identical to
    `rc_hybrid_block`'s for the same inputs (asserted).  Only the per-cell totals
    and the pass-2 design/placement change: cells take class levels L_k (summing to
    the same grand total S) and gated histograms are designed per class against
    per-class Phi tables, then placed per class.

    `classes`: sequence of (fraction, raw library level); default None = welded
    (a single class, the verbatim incumbent path — bit-identical).  K == 1 with an
    explicit level is degenerate (no steering freedom) and also takes the welded
    path.  `lfc_target`: None = unsteered (mass-proportional split, delivered
    arithmetic mean == bulk CPM); a finite array over the gated (DE-universe)
    genes = per-gene arithmetic-mean CPM targets (NaN = unsteered); or
    {"scale": s, "mask": gated bool} resolved post-pass-1 to
    ref_lin * 2^(s * design-LFC); optional "permute": seed deranges the resolved
    targets across genes (shuffled-class control).  Steering with K == 1 raises
    ValueError.

    Invariants, asserted per block: (1) per-gene totals bit-identical to
    `rc_hybrid_block`; (2) per-cell totals equal the class levels (up to the S mod
    n remainder convention); (3) pooled Wilcoxon |z| realised in per-cell CPM
    space (measured exactly with `wilcoxon_z`); (4) steered genes' delivered
    arithmetic mean within half a one-count move of target (clips flagged).
    """
    if emit_counts is None:                       # pragma: no cover - always injected
        from vcc2026.donor_transfer import emit_counts as emit_counts  # noqa: PLC0415

    if pass2_rng is None:
        pass2_rng = rng
    if call_set is not None and not call_set_first:
        raise ValueError("a prescribed call set needs call_set_first=True")
    Kcls = 1 if classes is None else len(classes)
    if Kcls == 1 and lfc_target is not None:
        if isinstance(lfc_target, dict) or np.any(
                np.isfinite(np.asarray(lfc_target, dtype=np.float64))):
            raise ValueError("steering needs K > 1 library classes")
    n = int(len(lib_sizes))
    n_genes = int(np.asarray(profile).size)
    keep = control.keep
    if keep.size != n_genes:
        raise ValueError("control block gate and profile disagree on the gene axis")
    gidx = np.flatnonzero(keep)
    n_keep = int(gidx.size)
    spread_vector = np.asarray(spread, dtype=np.float64)
    if spread_vector.ndim == 0:
        spread_vector = np.full(n_genes, float(spread_vector))
    if spread_vector.size != n_genes:
        raise ValueError("spread must be a scalar or one value per gene")

    # ---- pass 1: verbatim incumbent (same rng stream => identical totals) --------
    values, columns, pointers = emit_counts(
        profile, lib_sizes, dispersion, rng, max_counts=max_counts,
        mode="quasi", spread=spread_vector, overdispersion=overdispersion,
    )
    incumbent = np.zeros((n, n_genes), dtype=np.int64)
    for row in range(n):
        lo, hi = int(pointers[row]), int(pointers[row + 1])
        incumbent[row, columns[lo:hi]] = values[lo:hi].astype(np.int64)
    totals = incumbent.sum(0)
    z1, _p1, q1, _lfc1 = de_from_block(incumbent, control)
    called = q1 < 0.05
    del incumbent

    # ---- per-cell totals: welded (verbatim) or class levels ----------------------
    grand = int(totals.sum())
    if Kcls == 1:
        per_cell = grand // n
        cell_totals = np.full(n, per_cell, dtype=np.int64)
        cell_totals[: grand % n] += 1
        cell_class = np.zeros(n, dtype=np.int64)
        Lk = np.array([per_cell], dtype=np.int64)
        nk = np.array([n], dtype=np.int64)
    else:
        cell_class, Lk, cell_totals, nk = class_cells(n, classes, grand)
        per_cell = grand // n

    # ---- assignment: verbatim incumbent re-key of the |z| multiset ----------------
    key_gated = np.asarray(key, dtype=np.float64)[keep]
    key_abs = np.abs(key_gated)
    if not np.all(np.isfinite(key_abs)):
        raise ValueError("the ordering key carries non-finite entries")
    if float(key_abs.max(initial=0.0)) >= CALL_SET_OFFSET:
        raise ValueError("the ordering key is too large for the call-set block offset")
    jitter = np.random.default_rng(jitter_seed).random(n_keep) * JITTER_SCALE
    called_pass1 = called
    if call_set is not None:
        if call_set.size > n_keep:
            raise ValueError("the prescribed call set is larger than the DE universe")
        head = np.zeros(n_keep, dtype=bool)
        head[np.argsort(-key_abs, kind="stable")[: int(call_set.n_head)]] = True
        cpm_rank = np.where(head, -np.inf, np.asarray(control.ref_lin, dtype=np.float64))
        called = head.copy()
        called[np.argsort(-cpm_rank, kind="stable")[: int(call_set.n_cpm)]] = True
        if int(called.sum()) != call_set.size:
            raise AssertionError("the prescribed call set did not come out at its size")
    block_bonus = (
        np.where(called, CALL_SET_OFFSET, 0.0) if call_set_first
        else np.zeros(n_keep, dtype=np.float64)
    )
    if order_key is None:
        order_abs = key_abs
    else:
        order_abs = np.asarray(order_key, dtype=np.float64)[keep]
        if not np.all(np.isfinite(order_abs)) or float(order_abs.min(initial=0.0)) < 0.0:
            raise ValueError("the order key must be finite and non-negative")
        if float(order_abs.max(initial=0.0)) >= CALL_SET_OFFSET:
            raise ValueError("the order key is too large for the call-set block offset")
    ranked = np.argsort(-(order_abs + jitter + block_bonus), kind="stable")
    position = np.empty(n_keep, dtype=np.int64)
    position[ranked] = np.arange(n_keep)
    z_target = np.sort(np.abs(z1))[::-1][position]
    if call_set is not None:
        complement = position >= call_set.size
        z_worst = float(z_target[complement].max(initial=0.0))
        if z_worst > call_set.z_cap:
            z_target = z_target.copy()
            z_target[complement] *= call_set.z_cap / z_worst
        if float(z_target[complement].max(initial=0.0)) > call_set.z_cap + 1e-9:
            raise AssertionError("the complement was not targeted below the cap")
    signs = np.sign(key_gated)
    signs[signs == 0] = 1.0

    def _emit_sub():
        # verbatim incumbent sub-gate block (same rng stream)
        sub_profile = np.where(keep, 0.0, np.asarray(profile, dtype=np.float64))
        v, col, ptr = emit_counts(
            sub_profile, lib_sizes, dispersion, pass2_rng, max_counts=max_counts,
            mode="quasi", spread=spread_vector, overdispersion=overdispersion,
        )
        sub = np.zeros((n, n_genes), dtype=np.int64)
        for row in range(n):
            lo, hi = int(ptr[row]), int(ptr[row + 1])
            sub[row, col[lo:hi]] = v[lo:hi].astype(np.int64)
        if sub[:, keep].any():
            raise AssertionError("the sub-gate block leaked into the DE universe")
        return sub

    # ---- pass-2 design ------------------------------------------------------------
    # K == 1 keeps the verbatim incumbent order (design, sub, place); K > 1 emits
    # the sub-block BEFORE design so splits use true per-class caps.  The order
    # swap is stream-safe: design takes no rng draws.
    steering = None
    dose_rec = {}
    if Kcls == 1:
        designed = EME(control.cs, control.tie, per_cell, n).design(
            totals[keep], signs * z_target, mode="abs", T=per_cell
        )
        Mk = None
        sub_block = _emit_sub()
    else:
        sub_block = _emit_sub()
        Ms = totals[keep].astype(np.int64)
        ref_lin = np.asarray(control.ref_lin, dtype=np.float64)
        if isinstance(lfc_target, dict) and "prior" in lfc_target:
            # 391 (PI, 2026-09-13): cross-line sign prior on the called CPM block at the delivered magnitude; down-signed
            # magnitudes shrunk by one factor so the block's arithmetic mass is unchanged (390 C5b). A feasible dose
            # lam in [0, 1] interpolates from the welded delivery (lam 0) to the full target, mirroring feasible_dose.
            submass = np.array([sub_block[cell_class == k].sum() for k in range(Kcls)], dtype=np.int64)
            cap1d = nk * Lk - submass
            if np.any(cap1d < 0):
                raise AssertionError("sub-gate mass overflows a class")
            caps = np.broadcast_to(cap1d.astype(np.int64), (n_keep, Kcls)).copy()
            A0 = 1e6 * Ms.astype(np.float64) / float(grand)
            with np.errstate(divide="ignore", invalid="ignore"):
                E0v = np.log2(np.maximum(A0, 1e-300) / np.maximum(ref_lin, 1e-300))
            kg = np.abs(np.asarray(key_gated, dtype=np.float64))
            head = np.zeros(n_keep, dtype=bool)
            head[np.argsort(-kg, kind="stable")[: int(call_set.n_head)]] = True
            blk = np.zeros(n_keep, dtype=bool)
            blk[np.argsort(-np.where(head, -np.inf, ref_lin), kind="stable")[: int(call_set.n_cpm)]] = True
            blk &= (Ms > 0) & np.isfinite(E0v)
            if "extra_mask" in lfc_target:
                # 400: uncalled, non-head genes whose neighbourhood sign is confident join the steered set (and the mass balance)
                xm = np.asarray(lfc_target["extra_mask"], dtype=bool)
                if xm.shape != (n_keep,):
                    raise ValueError("extra_mask must match the gated gene axis")
                blk |= xm & ~head & (Ms > 0) & np.isfinite(E0v)
            psign = np.asarray(lfc_target["prior"], dtype=np.float64)
            if psign.shape != (n_keep,):
                raise ValueError("prior must match the gated gene axis")
            if "permute" in lfc_target:
                psign = psign[np.random.default_rng(int(lfc_target["permute"])).permutation(n_keep)]
            floor = float(lfc_target.get("shrink_floor", 0.002))
            tgt = psign * np.abs(E0v)
            pq = lfc_target.get("protect_q")
            if pq is not None:
                # 394: flip only block genes with |E0| <= the block's q-quantile; the strongest keep their delivered LFC
                thr = float(np.quantile(np.abs(E0v[blk]), float(pq))) if blk.any() else 0.0
                tgt = np.where(blk & (np.abs(E0v) <= thr), tgt, E0v)
            if "gate" in lfc_target:
                # 416 (E1): the prior may flip a block gene only where the gate is True (CLEAN4 banks unanimous on the
                # side); every other block gene keeps its delivered LFC, exactly like a protected gene. Absent = 400.
                gm = np.asarray(lfc_target["gate"], dtype=bool)
                if gm.shape != (n_keep,):
                    raise ValueError("gate must match the gated gene axis")
                tgt = np.where(gm, tgt, E0v)
            side = blk & (tgt < 0)

            def _t2(kf):
                return np.where(side, np.minimum(tgt * kf, -floor), tgt)

            def _R(kf):
                t2 = _t2(kf)
                return float((A0[blk] * (np.power(2.0, t2[blk] - E0v[blk]) - 1.0)).sum())
            r_lo, r_hi = _R(0.0), _R(1.0)
            if r_lo * r_hi < 0:
                a_, b_ = 0.0, 1.0
                for _ in range(50):
                    m_ = 0.5 * (a_ + b_)
                    if _R(a_) * _R(m_) <= 0:
                        b_ = m_
                    else:
                        a_ = m_
                kf = 0.5 * (a_ + b_)
            else:
                kf = 0.0 if abs(r_lo) < abs(r_hi) else 1.0
            move = np.where(blk, _t2(kf) - E0v, 0.0)
            steer_set = blk
            move_m = None
            mde_rec = {}
            if "mde_add" in lfc_target:
                # 426: MDE-structured magnitude. nmae reads only the reference's significant genes, whose |lfc| is bounded
                # below by the minimum detectable effect (~ 1/sqrt(mean count)); adding sign x add to the delivered LFC gains
                # (2a - 1) x add per gene while add < |l| (a = sign accuracy). A-channel only: per-gene totals (bulk), the
                # pass-2 |z| design (calls, reach order) and every sign are unchanged. The gated genes' arithmetic mass is
                # fixed by the class levels, so the add is mass-balanced by one factor on the heavier side.
                madd = np.asarray(lfc_target["mde_add"], dtype=np.float64)
                if madd.shape != (n_keep,):
                    raise ValueError("mde_add must match the gated gene axis")
                scope = str(lfc_target.get("mde_scope", "all"))
                base_t = np.where(blk, _t2(kf), E0v)
                base_t0 = base_t.copy()
                sgn = np.sign(base_t)
                flipm = np.zeros(n_keep, dtype=bool)
                if "mde_sign" in lfc_target:
                    # 426b: a sign source for the NON-block genes (e.g. the CLEAN4 per-gene prior); 0 = keep the delivered
                    # sign. A gene whose delivered sign disagrees is flipped at its own magnitude before the add.
                    msg = np.asarray(lfc_target["mde_sign"], dtype=np.float64)
                    if msg.shape != (n_keep,):
                        raise ValueError("mde_sign must match the gated gene axis")
                    use = ~blk & (msg != 0) & (Ms > 0) & np.isfinite(E0v)
                    flipm = use & (np.sign(E0v) != msg)
                    base_t = np.where(use, msg * np.abs(np.where(np.isfinite(E0v), E0v, 0.0)), base_t)
                    sgn = np.where(use, msg, sgn)
                elig = (Ms > 0) & np.isfinite(E0v) & np.isfinite(madd) & (madd > 0) & (sgn != 0)
                if scope == "block":
                    elig &= blk
                elif scope == "nonblock":
                    elig &= ~blk
                elif scope != "all":
                    raise ValueError(f"mde_scope {scope!r}")
                if "mde_agree" in lfc_target:
                    # 426e: the add only where the delivered sign agrees with an independent per-gene sign (the CLEAN4
                    # prior where covered; 0 = silent, add kept). Where they disagree our sign is right < 1/2 of the
                    # time on H1 (0.45), so magnitude there costs nmae. No flips: fewer steered genes, no extra mass.
                    agr = np.asarray(lfc_target["mde_agree"], dtype=np.float64)
                    if agr.shape != (n_keep,):
                        raise ValueError("mde_agree must match the gated gene axis")
                    elig &= (agr == 0) | (agr == sgn)
                if bool(lfc_target.get("mde_concordant", False)):
                    # the add only where the delivered LFC sign agrees with the |z| design's sign (sign(key)): on a
                    # prior-flipped gene the add would push A against the z design and clip it (the reach order)
                    elig &= sgn == signs
                up = elig & (sgn > 0)
                dn = elig & (sgn < 0)
                Bm = A0 * np.power(2.0, np.where(np.isfinite(base_t - E0v), base_t - E0v, 0.0))
                # the flips' own mass change (0 without mde_sign) enters the balance with the adds
                Fm = float((A0[flipm] * (np.power(2.0, base_t[flipm] - E0v[flipm]) - 1.0)).sum()) if flipm.any() else 0.0

                def _U(f):
                    return float((Bm[up] * (np.power(2.0, f * madd[up]) - 1.0)).sum())

                def _D(f):
                    return float((Bm[dn] * (1.0 - np.power(2.0, -f * madd[dn]))).sum())
                U1, D1 = _U(1.0), _D(1.0)
                fu = fd = 1.0
                if U1 + Fm > D1 and U1 > 0:
                    a_, b_ = 0.0, 1.0
                    for _ in range(60):
                        m_ = 0.5 * (a_ + b_)
                        if _U(m_) + Fm > D1:
                            b_ = m_
                        else:
                            a_ = m_
                    fu = 0.5 * (a_ + b_)
                elif D1 > U1 + Fm and D1 > 0:
                    a_, b_ = 0.0, 1.0
                    for _ in range(60):
                        m_ = 0.5 * (a_ + b_)
                        if _D(m_) > U1 + Fm:
                            b_ = m_
                        else:
                            a_ = m_
                    fd = 0.5 * (a_ + b_)
                move_m = (base_t - base_t0) + np.where(up, fu * madd, 0.0) - np.where(dn, fd * madd, 0.0)
                move_m = np.where(np.isfinite(move_m), move_m, 0.0)
                steer_set = blk | elig | flipm
                mde_rec = {"mde_n_up": int(up.sum()), "mde_n_dn": int(dn.sum()), "mde_fu": float(fu), "mde_fd": float(fd),
                           "mde_mass_residual": float(_U(fu) - _D(fd) + Fm), "mde_n_flip": int(flipm.sum()),
                           "mde_scope_code": float({"all": 0, "block": 1, "nonblock": 2}[scope])}
            w_ = (nk * Lk).astype(np.float64)
            w_ /= w_.sum()
            base_ = split_massprop(Ms, w_)
            margin_ = 6000
            allow_m = None

            def _demand(lam, lam_m=0.0):
                if lam_m > 0.0:
                    A_t_ = np.where(steer_set, A0 * np.power(2.0, lam * move + lam_m * move_m), np.nan)
                else:
                    A_t_ = np.where(blk, A0 * np.power(2.0, lam * move), np.nan)
                Mk_ = base_.copy()
                st_ = np.isfinite(A_t_)
                if st_.any():
                    Mk_s_, A_ach_, chi_, clo_ = split_for_A(Ms[st_], A_t_[st_], Lk, nk, caps[st_])
                    Mk_[st_] = Mk_s_
                else:
                    A_ach_, chi_, clo_ = np.array([]), np.array([]), np.array([])
                ok_ = bool(np.all(Mk_.sum(0) - (nk * Lk - margin_) <= 0))
                if lam_m > 0.0 and allow_m is not None:
                    # 426: a class's over-fill (gated + sub-gate - capacity) must be moved out cell by cell through the
                    # sub-gate counts; past ~0.58 of the class's sub-gate mass the repair stalled on H1 ("did not meet
                    # the per-cell totals"). The MDE dose keeps every class at <= allow_m x its sub-gate mass.
                    ok_ = ok_ and bool(np.all(Mk_.sum(0) + submass - nk * Lk <= allow_m * submass))
                return ok_, Mk_, A_t_, A_ach_, chi_, clo_
            ok, Mk, A_t, A_ach, c_hi, c_lo = _demand(1.0)
            lam = 1.0
            if not ok:
                a_, b_, best = 0.0, 1.0, None
                for _ in range(10):
                    m_ = 0.5 * (a_ + b_)
                    res_ = _demand(m_)
                    if res_[0]:
                        a_, best = m_, (m_,) + res_[1:]
                    else:
                        b_ = m_
                if best is None:
                    lam, Mk, A_t = 0.0, base_, np.full(n_keep, np.nan)
                    A_ach, c_hi, c_lo = np.array([]), np.array([]), np.array([])
                else:
                    lam, Mk, A_t, A_ach, c_hi, c_lo = best
            if move_m is not None and lam > 0.0:
                # 426: the prior keeps its own dose; the MDE increment gets the largest feasible dose on top of it
                allow_m = lfc_target.get("mde_overfill", None)
                res_ = _demand(lam, 1.0)
                lam_m = 1.0
                if not res_[0]:
                    a_, b_, bm = 0.0, 1.0, None
                    for _ in range(10):
                        m_ = 0.5 * (a_ + b_)
                        r2_ = _demand(lam, m_)
                        if r2_[0]:
                            a_, bm = m_, r2_
                        else:
                            b_ = m_
                    res_, lam_m = (bm, a_) if bm is not None else (None, 0.0)
                if res_ is not None and lam_m > 0.0:
                    _, Mk, A_t, A_ach, c_hi, c_lo = res_
                mde_rec["mde_lambda"] = float(lam_m)
                mde_rec["mde_overfill_frac"] = float(np.max((np.asarray(Mk).sum(0) + submass - nk * Lk) / np.maximum(submass, 1)))
            dose_rec = {"dose_nominal": 1.0, "dose_delivered": float(lam)}
            LAST_PRIOR_REC.append({"prior_shrink_k": float(kf), "prior_lambda": float(lam),
                                   "prior_n_block": int(blk.sum()),
                                   "prior_n_flip_target": int((blk & (np.sign(_t2(kf)) != np.sign(E0v))).sum()),
                                   "prior_mass_residual": float(_R(kf)), **mde_rec})
            steered = np.isfinite(A_t)
            if A_ach.size:
                steering = dict(
                    n_steered=int(steered.sum()),
                    clip_hi=int(c_hi.sum()), clip_lo=int(c_lo.sum()),
                    target=A_ach, idx=np.flatnonzero(steered),
                    clipped=np.ascontiguousarray(c_hi | c_lo),
                    **dose_rec)
            else:
                steering = None
        elif isinstance(lfc_target, dict):
            # steered block: per-target feasible dose (class-level pigeonhole:
            # joint 1.48x demand can exceed a class's mass; bisect to s*).
            submass = np.array([sub_block[cell_class == k].sum()
                                for k in range(Kcls)], dtype=np.int64)
            cap1d = nk * Lk - submass
            if np.any(cap1d < 0):
                raise AssertionError("sub-gate mass overflows a class")
            caps = np.broadcast_to(cap1d.astype(np.int64),
                                   (n_keep, Kcls)).copy()
            s_nom = float(lfc_target["scale"])
            msk = np.asarray(lfc_target["mask"], dtype=bool)
            s_star, Mk, A_t, A_ach, c_hi, c_lo = feasible_dose(
                Ms, ref_lin, grand, s_nom, msk, Lk, nk, caps, nk * Lk,
                perm=lfc_target.get("permute"))
            dose_rec = {"dose_nominal": s_nom, "dose_delivered": float(s_star)}
            steered = np.isfinite(A_t)
            if A_ach.size:
                steering = dict(
                    n_steered=int(steered.sum()),
                    clip_hi=int(c_hi.sum()), clip_lo=int(c_lo.sum()),
                    target=A_ach, idx=np.flatnonzero(steered),
                    clipped=np.ascontiguousarray(c_hi | c_lo),
                    **dose_rec)
            else:
                steering = None
        elif lfc_target is None:
            # unsteered block (frozen S1 path: plain massprop, global repair).
            A_t = np.full(n_keep, np.nan)
            steered = np.isfinite(A_t)
            w = (nk * Lk).astype(np.float64)
            w /= w.sum()
            Mk = split_massprop(Ms, w)
        else:
            A_t = np.asarray(lfc_target, dtype=np.float64)
            if A_t.shape != (n_keep,):
                raise ValueError("lfc_target must match the gated gene axis")
            steered = np.isfinite(A_t)
            w = (nk * Lk).astype(np.float64)
            w /= w.sum()
            Mk = split_massprop(Ms, w)
            if steered.any():
                caps = np.broadcast_to((nk * Lk).astype(np.int64),
                                       (n_keep, Kcls)).copy()
                Mk_s, A_ach, c_hi, c_lo = split_for_A(
                    Ms[steered], A_t[steered], Lk, nk, caps[steered])
                Mk[steered] = Mk_s
                steering = dict(
                    n_steered=int(steered.sum()),
                    clip_hi=int(c_hi.sum()), clip_lo=int(c_lo.sum()),
                    target=A_ach, idx=np.flatnonzero(steered),
                    clipped=np.ascontiguousarray(c_hi | c_lo))
        cmaxs, og, lt_all, eq_all, Phi_all = build_tables(
            control.cs, control.n_controls, Mk, Lk)
        sd_pool = sd_pool_all(Mk, nk, Lk, cmaxs, og, lt_all, eq_all,
                              np.ascontiguousarray(control.tie),
                              control.n_controls, Kcls)
        KK, LL, MM, zr, ze = design_all_k(
            Mk, nk, Lk, cmaxs, og, lt_all, eq_all, Phi_all,
            np.ascontiguousarray(control.tie), control.n_controls,
            np.ascontiguousarray(signs * z_target), sd_pool, Kcls)
        designed = dict(K=KK, levels=LL, mults=MM, z=zr, z_even=ze)

    # ---- place ---------------------------------------------------------------------
    if Kcls == 1:
        block, diagnostics = EME(control.cs, control.tie, per_cell, n).place(
            designed, gidx, sub_block, cell_totals, pass2_rng
        )
    else:
        dsg = {"K": designed["K"], "levels": designed["levels"],
               "mults": designed["mults"]}
        block, diagnostics = place_k(dsg, gidx, sub_block, cell_totals,
                                     cell_class, Lk, nk, pass2_rng)

    delivered = block.sum(0)
    if not np.array_equal(delivered, totals):
        raise AssertionError("rc_hybrid_k moved a per-gene panel total")
    if not np.array_equal(block.sum(1), cell_totals):
        if Kcls > 1 and Mk is not None:
            _gm = np.asarray(Mk).sum(0)
            _sm = np.array([sub_block[cell_class == k].sum() for k in range(Kcls)])
            print(f"[426 dbg] class gated {_gm.tolist()} sub {_sm.tolist()} cap {(nk * Lk).tolist()} "
                  f"slack {((nk * Lk) - _gm - _sm).tolist()} diag {diagnostics}", flush=True)
        raise AssertionError("rc_hybrid_k did not meet the per-cell totals")
    if block.min() < 0:
        raise AssertionError("rc_hybrid_k emitted a negative count")
    realised, _p2, q2, lfc2 = de_from_block(block, control)
    if Kcls > 1 and steering is not None:
        # exact delivered class totals per steered gene from the placed design
        # (repair_rows never touches gated genes, so design M_k == delivered M_k)
        KK = designed["K"]
        LL = designed["levels"]
        MM = designed["mults"]
        idx = steering["idx"]
        Mk_del = np.array(
            [[int((LL[i, k, :int(KK[i, k])] * MM[i, k, :int(KK[i, k])]).sum())
              for k in range(Kcls)] for i in idx])
        if not np.array_equal(Mk_del, Mk[idx]):
            raise AssertionError("placed class totals differ from the split")
        A_del = arith_of_split(Mk_del, Lk, n)
        A_t = steering["target"]
        with np.errstate(divide="ignore", invalid="ignore"):
            err = np.abs(np.log2(np.maximum(A_del, 1e-12) / np.maximum(A_t, 1e-12)))
        err[(A_del == 0) & (A_t == 0)] = 0.0
        diagnostics.update(
            steer_max_log2_err=float(err.max()) if err.size else 0.0,
            steer_mean_log2_err=float(err.mean()) if err.size else 0.0,
            steer_clip_hi=int(steering["clip_hi"]),
            steer_clip_lo=int(steering["clip_lo"]),
        )
    if Kcls > 1 and dose_rec:
        diagnostics.update(dose_nominal=float(dose_rec["dose_nominal"]),
                           dose_delivered=float(dose_rec["dose_delivered"]))
        # tolerance: half a one-count move in log2 at the target's scale, plus the
        # cell-remainder lattice (levels are exact up to the S mod n fix-up)
        granule = (1e6 / n) * float(abs(1.0 / Lk.min() - 1.0 / Lk.max()))
        tol = granule / np.maximum(A_t, 1e-9) / np.log(2) / 2 + 0.02
        ok = np.isfinite(err) & ~steering["clipped"]
        bad = ok & (err > tol + 1e-9)
        if bad.any():
            raise AssertionError(
                f"steering tolerance exceeded on {int(bad.sum())} genes")
    diagnostics.update(
        call_set_first=float(bool(call_set_first)),
        n_pred_pass1=int(called_pass1.sum()),
        n_genes_de_universe=n_keep,
        per_cell_total=int(per_cell),
        n_classes=int(Kcls),
        z_realised_vs_designed_max=float(np.abs(realised - designed["z"]).max()),
        clip_frac=float(np.mean(np.abs(np.abs(designed["z"]) - z_target) > 0.5)),
        sign_flip_frac=float(
            np.mean((np.sign(designed["z"]) != signs) & (z_target > 0.5))
        ),
        sub_gate_spread_max=float(spread_vector[~keep].max(initial=0.0)),
        sub_gate_counts_per_cell=float(sub_block.sum(1).mean()),
    )
    if call_set is not None:
        delivered_call = q2 < 0.05
        n_delivered = int(delivered_call.sum())
        diagnostics.update(
            call_set_prescribed=int(call_set.size),
            call_set_delivered=n_delivered,
            call_set_jaccard=float(
                (delivered_call & called).sum() / max((delivered_call | called).sum(), 1)
            ),
            call_set_unsuppressible=int((delivered_call & ~called).sum()),
            call_set_missed=int((called & ~delivered_call).sum()),
            call_set_z_target_max=float(np.abs(z_target[complement]).max(initial=0.0)),
            call_set_z_scale=float(
                min(1.0, call_set.z_cap / z_worst) if z_worst > 0 else 1.0
            ),
        )
    return block, diagnostics


@njit(parallel=True, cache=True)
def sd_pool_all(Ms, nks, Lk, cmaxs, og, lt_all, eq_all, ctrl_ties, m, Kcls):
    """Per-gene Fbar sd of the pooled test under the even per-class allocation.

    sd_F = sqrt(tc (N+1) / (12 n m)) with tc from the pooled even histogram's exact
    tie accumulation (grouped across classes, same grouping as hist_z_k).
    """
    G = Ms.shape[0]
    n = 0
    for k in range(Kcls):
        n += nks[k]
    N = n + m
    out = np.empty(G, np.float64)
    for g in prange(G):
        Lv = np.empty(PMAX, np.int64)
        Mu = np.empty(PMAX, np.int64)
        Cv = np.empty(PMAX, np.int64)
        Ktot = 0
        lv = np.zeros(NLEV, np.int64)
        mu = np.zeros(NLEV, np.int64)
        for k in range(Kcls):
            Mk = Ms[g, k]
            Ke = histogram_of(Mk, int(nks[k]), int(cmaxs[g, k]), 0, 0, -1, lv, mu)
            for q in range(Ke):
                Lv[Ktot] = lv[q]
                Mu[Ktot] = mu[q]
                Cv[Ktot] = k
                Ktot += 1
        gv_lt = np.empty(PMAX, np.float64)
        gv_cc = np.empty(PMAX, np.float64)
        gv_mult = np.empty(PMAX, np.int64)
        ng = _pgroup(Lv, Mu, Cv, Ktot, Lk, lt_all, eq_all,
                     og[g], gv_lt, gv_cc, gv_mult)
        tie = ctrl_ties[g]
        for q in range(ng):
            tot = gv_cc[q] + gv_mult[q]
            tie += (tot ** 3 - tot) - (gv_cc[q] ** 3 - gv_cc[q])
        tc = 1.0 - tie / (float(N) ** 3 - float(N))
        out[g] = np.sqrt(max(tc, 1e-12) * (N + 1.0) / (12.0 * n * m))
    return out

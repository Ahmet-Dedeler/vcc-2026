"""185 — the exact-marginal emitter (EME): a count matrix with PRESCRIBED per-gene
histograms, bit-identical per-gene totals, and a CONSTANT per-cell total.

THE OBJECT
----------
The scorer reads three different functionals of our 400 x G count block:

  * pseudobulk (per-gene total)                  -> pds, mse's profile half
  * arithmetic mean of per-cell CPM             -> every lfc, hence nmae and every SIGN
  * per-gene Wilcoxon rank sum vs real controls -> fid, jac, reach (through the p-order)

With every emitted cell at the same library total T, the second is 1e6 * M_g / (n T): a
function of the per-gene total alone.  So at fixed {M_g} and fixed T the first two are
pinned and the third is free, because a rank-sum statistic reads the SHAPE of the per-cell
histogram and the mean does not.  This module realises that freedom.

WHAT IS EXACT
-------------
For gene g let Phi_g(c) = P_ctrl(cpm < 1e6 c/T) + 0.5 P_ctrl(cpm == 1e6 c/T), read off the
real control cells.  The Mann-Whitney U of our n cells against the m controls is
U = sum_i Phi_g(c_i) * m  EXACTLY (midranks make U independent of ties among our own cells),
so Fbar = U/(nm) = mean_i Phi_g(c_i): a separable function of the histogram.  The reachable
set of Fbar at total M is therefore the segment between the lower-convex and upper-concave
envelopes of Phi_g evaluated at the mean M/n, and each endpoint is attained by a two-level
histogram (Caratheodory on one moment constraint).  The tie correction enters only the
denominator and is computed exactly for the chosen histogram.  `hist_z` below replicates
`de_kernel.wilcoxon_z` (scanpy-identical) for a histogram and is checked against it.

THE COUPLING, AND THE REPAIR
----------------------------
Per-gene histograms are prescribed for the DE-visible genes (control mean CPM > 5).  Their
levels are assigned to cells greedily against remaining per-cell capacity (largest level to
the emptiest cell, random tie-break), so row sums land within a few counts of T.  The genes
BELOW the 5-CPM gate -- which no DE member can see -- are emitted with a15's own law and then
used as the slack: single counts of sub-gate genes are moved between cells until every row
sums to exactly T.  Per-gene totals are never touched by any step (a move keeps the column
total); the caller hashes the delivered column totals and requires a15's digest.
"""
from __future__ import annotations

import numpy as np
from numba import njit, prange


# ----------------------------------------------------------------------------- Phi tables
@njit(cache=True, inline="always")
def _bs_left(a, v):
    lo, hi = 0, a.shape[0]
    while lo < hi:
        mid = (lo + hi) >> 1
        if a[mid] < v:
            lo = mid + 1
        else:
            hi = mid
    return lo


@njit(cache=True, inline="always")
def _bs_right(a, v):
    lo, hi = 0, a.shape[0]
    while lo < hi:
        mid = (lo + hi) >> 1
        if a[mid] <= v:
            lo = mid + 1
        else:
            hi = mid
    return lo


@njit(parallel=True, cache=True)
def phi_tables(cs, T, cmax):
    """For gene g and count level c in 0..cmax[g]: lt = #ctrl strictly below our value
    log1p(1e6 c/T), eq = #ctrl exactly equal.  Flat arrays with offsets."""
    G = cs.shape[0]
    offs = np.zeros(G + 1, np.int64)
    for g in range(G):
        offs[g + 1] = offs[g] + cmax[g] + 1
    lt = np.empty(offs[G], np.int32)
    eq = np.empty(offs[G], np.int32)
    for g in prange(G):
        c = cs[g]
        o = offs[g]
        for k in range(cmax[g] + 1):
            v = np.log1p(1e6 * k / T)
            lo = _bs_left(c, v)
            hi = _bs_right(c, v)
            lt[o + k] = lo
            eq[o + k] = hi - lo
    return offs, lt, eq


# ----------------------------------------------------------------------------- exact z
@njit(cache=True)
def hist_z(levels, mults, K, lt, eq, o, ctrl_tie, n, m):
    """Exact scanpy-style tie-corrected z for a histogram: K levels (ascending, distinct)
    with multiplicities.  lt/eq are the gene's Phi tables at offset o."""
    N = n + m
    R1 = 0.0
    tie = ctrl_tie
    below = 0
    for k in range(K):
        c = levels[k]
        nk = mults[k]
        if nk == 0:
            continue
        cc = float(eq[o + c])
        tot = cc + nk
        midrank = (lt[o + c] + below) + (tot + 1.0) / 2.0
        R1 += nk * midrank
        tie += (tot * tot * tot - tot) - (cc * cc * cc - cc)
        below += nk
    tc = 1.0 - tie / (float(N) ** 3 - float(N))
    sd = np.sqrt(tc * n * m * (N + 1.0) / 12.0)
    if sd <= 0.0:
        return 0.0
    return (R1 - n * (N + 1.0) / 2.0) / sd


@njit(cache=True)
def _fbar3(n, n0, base, rem, n1, Phi, o, nb, b):
    # histogram: n0 at 0, nb at b, (n1-rem) at base, rem at base+1
    f = n0 * Phi[o] + nb * Phi[o + b] + (n1 - rem) * Phi[o + base]
    if rem > 0:
        f += rem * Phi[o + base + 1]
    return f / n


@njit(cache=True)
def _upper_hull_right(Phi, o, cmax, mu):
    """Right endpoint b* of the upper concave hull segment of (c, Phi(c)), c = 0..cmax,
    that contains mu.  Monotone chain.  Returns -1 if no hull vertex lies above mu."""
    # Andrew's monotone chain, upper hull, x ascending
    hx = np.empty(cmax + 1, np.int64)
    hy = np.empty(cmax + 1, np.float64)
    h = 0
    for c in range(cmax + 1):
        y = Phi[o + c]
        while h >= 2:
            # cross product of (h-2 -> h-1) x (h-2 -> c); pop while not a right turn
            x1 = hx[h - 1] - hx[h - 2]
            y1 = hy[h - 1] - hy[h - 2]
            x2 = c - hx[h - 2]
            y2 = y - hy[h - 2]
            if x1 * y2 - y1 * x2 >= 0.0:
                h -= 1
            else:
                break
        hx[h] = c
        hy[h] = y
        h += 1
    for k in range(h):
        if hx[k] > mu:
            return hx[k]
    return -1


@njit(cache=True)
def design_gene(M, n, cmax, Phi, o, lt, eq, ctrl_tie, m, targets, n_targets, family_out,
                param_out):
    """Scan the two one-parameter histogram families for gene g and return, for each
    target Fbar, the candidate closest to it; also the family-wide min/max/even Fbar.

    Family Z(n0): n0 cells at 0, the rest share M as evenly as integers allow.
    Family H(nb): nb cells at the hull vertex b*, the rest share the remainder evenly.
    Returns (F_even, F_min, F_max, b_star) and fills family_out/param_out per target.
    """
    if M == 0 or cmax == 0:
        f0 = Phi[o]
        for t in range(n_targets):
            family_out[t] = 0
            param_out[t] = 0
        return f0, f0, f0, -1
    base_e = M // n
    rem_e = M % n
    F_even = _fbar3(n, 0, base_e, rem_e, n, Phi, o, 0, 0)
    F_min = F_even
    F_max = F_even
    best = np.empty(n_targets, np.float64)
    for t in range(n_targets):
        best[t] = abs(F_even - targets[t])
        family_out[t] = 0
        param_out[t] = 0
    # ---- family Z
    n1_min = (M + cmax - 1) // cmax           # fewest cells that can hold M
    if n1_min < 1:
        n1_min = 1
    n0max = n - n1_min
    for n0 in range(1, n0max + 1):
        n1 = n - n0
        base = M // n1
        rem = M % n1
        if base + (1 if rem > 0 else 0) > cmax:
            continue
        f = _fbar3(n, n0, base, rem, n1, Phi, o, 0, 0)
        if f < F_min:
            F_min = f
        if f > F_max:
            F_max = f
        for t in range(n_targets):
            d = abs(f - targets[t])
            if d < best[t]:
                best[t] = d
                family_out[t] = 0
                param_out[t] = n0
    # ---- family H
    mu = M / n
    bstar = _upper_hull_right(Phi, o, cmax, mu)
    if bstar > 0:
        nbmax = M // bstar
        if nbmax > n - 1:
            nbmax = n - 1
        for nb in range(1, nbmax + 1):
            Mr = M - nb * bstar
            n1 = n - nb
            base = Mr // n1
            rem = Mr % n1
            if base + (1 if rem > 0 else 0) > cmax:
                continue
            f = _fbar3(n, 0, base, rem, n1, Phi, o, nb, bstar)
            if f < F_min:
                F_min = f
            if f > F_max:
                F_max = f
            for t in range(n_targets):
                d = abs(f - targets[t])
                if d < best[t]:
                    best[t] = d
                    family_out[t] = 1
                    param_out[t] = nb
    return F_even, F_min, F_max, bstar


@njit(cache=True)
def histogram_of(M, n, cmax, family, param, bstar, levels, mults):
    """Materialise the (levels ascending, mults) histogram of a family member.  Up to 4
    levels: 0, base, base+1, b*.  Returns the number of levels K."""
    K = 0
    if M == 0 or cmax == 0:
        levels[0] = 0
        mults[0] = n
        return 1
    if family == 0:
        n0 = param
        n1 = n - n0
        base = M // n1
        rem = M % n1
        if n0 > 0:
            levels[K] = 0
            mults[K] = n0
            K += 1
        # base may be 0 (then merge with the zero level)
        if base == 0:
            if K > 0:
                mults[K - 1] += n1 - rem
            else:
                levels[K] = 0
                mults[K] = n1 - rem
                K += 1
        else:
            levels[K] = base
            mults[K] = n1 - rem
            K += 1
        if rem > 0:
            levels[K] = base + 1
            mults[K] = rem
            K += 1
    else:
        nb = param
        Mr = M - nb * bstar
        n1 = n - nb
        base = Mr // n1
        rem = Mr % n1
        if base == 0:
            levels[K] = 0
            mults[K] = n1 - rem
            K += 1
        else:
            levels[K] = base
            mults[K] = n1 - rem
            K += 1
        if rem > 0:
            levels[K] = base + 1
            mults[K] = rem
            K += 1
        # b* is above base+1 by construction of the hull (b* > mu >= base); it can equal
        # base+1 only when rem == 0 and mu is an integer: merge then.
        if K > 0 and levels[K - 1] == bstar:
            mults[K - 1] += nb
        else:
            levels[K] = bstar
            mults[K] = nb
            K += 1
    # drop zero-multiplicity entries
    j = 0
    for k in range(K):
        if mults[k] > 0:
            levels[j] = levels[k]
            mults[j] = mults[k]
            j += 1
    return j


@njit(parallel=True, cache=True)
def design_all(Ms, n, cmaxs, Phi, offs, lt, eq, ctrl_ties, m, ztargets, mode, sd_even):
    """Choose a histogram per gene for signed z targets.

    mode 0 ("signed"): target Fbar = 0.5 + z* sd; clip to the reachable side.
    mode 1 ("abs")   : try both 0.5 +- |z*| sd and keep the one whose realised |z| is closer
                       to |z*|; ties broken toward sign(z*).
    Returns per gene: K, levels(4), mults(4), z_realised, F_even, F_min, F_max, z_even.
    """
    G = Ms.shape[0]
    Ks = np.zeros(G, np.int64)
    L = np.zeros((G, 4), np.int64)
    Mu = np.zeros((G, 4), np.int64)
    zr = np.zeros(G, np.float64)
    Fe = np.zeros(G, np.float64)
    Fmn = np.zeros(G, np.float64)
    Fmx = np.zeros(G, np.float64)
    ze = np.zeros(G, np.float64)
    for g in prange(G):
        o = offs[g]
        M = Ms[g]
        cmax = cmaxs[g]
        sd = sd_even[g]
        zt = ztargets[g]
        targets = np.empty(2, np.float64)
        fam = np.empty(2, np.int64)
        par = np.empty(2, np.int64)
        if mode == 0:
            targets[0] = 0.5 + zt * sd
            nt = 1
        else:
            targets[0] = 0.5 + abs(zt) * sd
            targets[1] = 0.5 - abs(zt) * sd
            nt = 2
        F_even, F_min, F_max, bstar = design_gene(M, n, cmax, Phi, o, lt, eq, ctrl_ties[g], m,
                                                  targets, nt, fam, par)
        lv = np.zeros(4, np.int64)
        mu_ = np.zeros(4, np.int64)
        # even histogram z (the "natural" allocation)
        Ke = histogram_of(M, n, cmax, 0, 0, bstar, lv, mu_)
        ze[g] = hist_z(lv, mu_, Ke, lt, eq, o, ctrl_ties[g], n, m)
        bestk = 0
        bestz = 0.0
        bestd = 1e300
        for t in range(nt):
            K = histogram_of(M, n, cmax, fam[t], par[t], bstar, lv, mu_)
            z = hist_z(lv, mu_, K, lt, eq, o, ctrl_ties[g], n, m)
            if mode == 0:
                d = abs(z - zt)
            else:
                d = abs(abs(z) - abs(zt))
                if (z > 0) != (zt > 0):
                    d += 1e-9          # tie-break toward sign(z*)
            if d < bestd:
                bestd = d
                bestk = t
                bestz = z
        K = histogram_of(M, n, cmax, fam[bestk], par[bestk], bstar, lv, mu_)
        Ks[g] = K
        for k in range(4):
            L[g, k] = lv[k]
            Mu[g, k] = mu_[k]
        zr[g] = bestz
        Fe[g] = F_even
        Fmn[g] = F_min
        Fmx[g] = F_max
    return Ks, L, Mu, zr, Fe, Fmn, Fmx, ze


@njit(parallel=True, cache=True)
def envelope_all(Ms, n, cmaxs, Phi, offs, lt, eq, ctrl_ties, m):
    """Part A: per gene, the exact z of the even, most-concentrated, min-Fbar and max-Fbar
    histograms, plus the z closest to 0 in the two families (the reachable |z| floor)."""
    G = Ms.shape[0]
    out = np.zeros((G, 6), np.float64)   # z_even, z_conc, z_min, z_max, z_floor, bstar
    for g in prange(G):
        o = offs[g]
        M = Ms[g]
        cmax = cmaxs[g]
        lv = np.zeros(4, np.int64)
        mu_ = np.zeros(4, np.int64)
        targets = np.empty(3, np.float64)
        fam = np.empty(3, np.int64)
        par = np.empty(3, np.int64)
        targets[0] = -10.0          # far below any reachable Fbar -> min-Fbar member
        targets[1] = 10.0           # far above -> max-Fbar member
        targets[2] = 0.5            # closest to the null
        F_even, F_min, F_max, bstar = design_gene(M, n, cmax, Phi, o, lt, eq, ctrl_ties[g], m,
                                                  targets, 3, fam, par)
        K = histogram_of(M, n, cmax, 0, 0, bstar, lv, mu_)
        out[g, 0] = hist_z(lv, mu_, K, lt, eq, o, ctrl_ties[g], n, m)
        # most concentrated: fill cells at cmax
        if M > 0 and cmax > 0:
            q = M // cmax
            r = M % cmax
            K = 0
            n0 = n - q - (1 if r > 0 else 0)
            if n0 > 0:
                lv[K] = 0
                mu_[K] = n0
                K += 1
            if r > 0:
                lv[K] = r
                mu_[K] = 1
                K += 1
            if q > 0:
                lv[K] = cmax
                mu_[K] = q
                K += 1
            out[g, 1] = hist_z(lv, mu_, K, lt, eq, o, ctrl_ties[g], n, m)
        else:
            out[g, 1] = out[g, 0]
        for t in range(3):
            K = histogram_of(M, n, cmax, fam[t], par[t], bstar, lv, mu_)
            out[g, 2 + t] = hist_z(lv, mu_, K, lt, eq, o, ctrl_ties[g], n, m)
        out[g, 5] = bstar
    return out


# ----------------------------------------------------------------------------- assignment
@njit(cache=True)
def assign_cells(Ks, L, Mu, order, cap, B, rng_keys):
    """Place each gene's histogram into cells against remaining capacity.

    `order`: gene processing order (descending max level).  `cap`: per-cell remaining
    capacity (modified in place).  `B`: (n, G) int32 output (gated columns only).
    `rng_keys`: (len(order), n) float64 random tie-breakers, one row per gene.
    Returns the number of placements that exceeded a cell's remaining capacity (forced)."""
    n = cap.shape[0]
    forced = 0
    key = np.empty(n, np.float64)
    for gi in range(order.shape[0]):
        g = order[gi]
        K = Ks[g]
        if K == 0:
            continue
        # sort cells by remaining capacity descending, random tie-break
        for i in range(n):
            key[i] = -cap[i] + rng_keys[gi, i] * 0.5
        idx = np.argsort(key)
        pos = 0
        # levels are ascending in L; place the largest first
        for k in range(K - 1, -1, -1):
            lvl = L[g, k]
            if lvl == 0:
                break
            for _ in range(Mu[g, k]):
                c = idx[pos]
                pos += 1
                B[c, g] = lvl
                if cap[c] < lvl:
                    forced += 1
                cap[c] -= lvl
    return forced


@njit(cache=True)
def repair_rows(B_sub, resid, rng_perm):
    """Move single counts of sub-gate genes from over-full cells (resid < 0) to under-full
    cells (resid > 0) until every residual is zero.  Column totals are preserved by
    construction.  Returns (moves, unresolved)."""
    n, Gs = B_sub.shape
    moves = 0
    # cells with positive residual, served in order of largest residual
    while True:
        src = -1
        smin = 0
        dst = -1
        dmax = 0
        for i in range(n):
            if resid[i] < smin:
                smin = resid[i]
                src = i
            if resid[i] > dmax:
                dmax = resid[i]
                dst = i
        if src < 0 or dst < 0:
            break
        # find a sub-gate gene with a count in src, scanning in a random order
        found = -1
        for jj in range(Gs):
            j = rng_perm[(moves + jj) % Gs]
            if B_sub[src, j] > 0:
                found = j
                break
        if found < 0:
            break
        B_sub[src, found] -= 1
        B_sub[dst, found] += 1
        resid[src] += 1
        resid[dst] -= 1
        moves += 1
    unresolved = 0
    for i in range(n):
        unresolved += abs(resid[i])
    return moves, unresolved


# ----------------------------------------------------------------------------- driver
class EME:
    """Holds the control-side tables for one harness and emits designed blocks."""

    def __init__(self, cs, ctrl_tie, T, n=400):
        self.cs = cs                       # (G, m) sorted log1p CPM of the real controls
        self.ctrl_tie = np.ascontiguousarray(ctrl_tie, dtype=np.float64)
        self.T = int(T)
        self.n = int(n)
        self.m = int(cs.shape[1])
        self._cache = {}

    def tables(self, Ms, T=None):
        """Phi tables for the given per-gene totals (cmax = min(T, M)) at per-cell total T
        (this TARGET's own S // n; the +1 remainder cells are ignored here and the realised
        z is re-measured exactly by the kernel afterwards)."""
        T = self.T if T is None else int(T)
        Ms = np.ascontiguousarray(Ms, dtype=np.int64)
        cmax = np.minimum(Ms, T).astype(np.int64)
        offs, lt, eq = phi_tables(self.cs, float(T), cmax)
        Phi = (lt.astype(np.float64) + 0.5 * eq.astype(np.float64)) / self.m
        return cmax, offs, lt, eq, Phi

    def sd_even(self, Ms, cmax, offs, lt, eq):
        """Fbar standard deviation under the even histogram's tie correction: converts a
        z target into an Fbar target.  sd_F = sqrt(tc (N+1) / (12 n m))."""
        n, m = self.n, self.m
        N = n + m
        G = Ms.size
        out = np.empty(G)
        lv = np.zeros(4, np.int64)
        mu = np.zeros(4, np.int64)
        for g in range(G):
            K = histogram_of(int(Ms[g]), n, int(cmax[g]), 0, 0, -1, lv, mu)
            o = offs[g]
            tie = self.ctrl_tie[g]
            for k in range(K):
                cc = float(eq[o + lv[k]])
                tot = cc + mu[k]
                tie += (tot ** 3 - tot) - (cc ** 3 - cc)
            tc = 1.0 - tie / (float(N) ** 3 - N)
            out[g] = np.sqrt(max(tc, 1e-12) * (N + 1.0) / (12.0 * n * m))
        return out

    def envelope(self, Ms, T=None):
        cmax, offs, lt, eq, Phi = self.tables(Ms, T)
        Ms = np.ascontiguousarray(Ms, dtype=np.int64)
        return envelope_all(Ms, self.n, cmax, Phi, offs, lt, eq, self.ctrl_tie, self.m)

    def design(self, Ms, ztargets, mode="abs", T=None):
        Ms = np.ascontiguousarray(Ms, dtype=np.int64)
        cmax, offs, lt, eq, Phi = self.tables(Ms, T)
        sd = self.sd_even(Ms, cmax, offs, lt, eq)
        zt = np.ascontiguousarray(ztargets, dtype=np.float64)
        Ks, L, Mu, zr, Fe, Fmn, Fmx, ze = design_all(
            Ms, self.n, cmax, Phi, offs, lt, eq, self.ctrl_tie, self.m, zt,
            0 if mode == "signed" else 1, sd)
        return dict(K=Ks, levels=L, mults=Mu, z=zr, F_even=Fe, F_min=Fmn, F_max=Fmx,
                    z_even=ze, sd=sd)

    def place(self, dsg, gated_idx, B_sub, T_cells, rng):
        """Assign the designed histograms of the gated genes to cells and repair with the
        sub-gate block.  Returns the full (n, G_full) int64 block and diagnostics."""
        n = self.n
        G_full = B_sub.shape[1]
        cap = T_cells.astype(np.int64) - B_sub.sum(1)
        Ks, L, Mu = dsg["K"], dsg["levels"], dsg["mults"]
        maxlvl = np.array([L[g, Ks[g] - 1] if Ks[g] > 0 else 0 for g in range(Ks.size)])
        order = np.argsort(-maxlvl, kind="stable").astype(np.int64)
        Bg = np.zeros((n, Ks.size), np.int32)
        keys = rng.random((order.size, n))
        forced = assign_cells(Ks, L, Mu, order, cap, Bg, keys)
        B = B_sub.astype(np.int64).copy()
        B[:, gated_idx] += Bg.astype(np.int64)
        resid = T_cells.astype(np.int64) - B.sum(1)
        max_abs_resid = int(np.abs(resid).max())
        sub_idx = np.setdiff1d(np.arange(G_full), gated_idx)
        Bs = np.ascontiguousarray(B[:, sub_idx])
        perm = rng.permutation(sub_idx.size).astype(np.int64)
        moves, unresolved = repair_rows(Bs, resid, perm)
        B[:, sub_idx] = Bs
        return B, dict(forced=int(forced), max_abs_resid=max_abs_resid, moves=int(moves),
                       unresolved=int(unresolved))

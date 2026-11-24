"""Exact replication of scanpy's wilcoxon(reference=..., tie_correct=True) rank-sum,
computed against a PRE-SORTED, fixed reference block (the real control), so the
O(N log N) joint rank is replaced by O(n1 log n2) per gene."""
import numpy as np
import numba
from numba import njit, prange


@njit(cache=True)
def _ctrl_tiesum_row(a):
    """sum over distinct values of (t^3 - t) for a SORTED row."""
    s = 0.0
    i = 0
    m = a.shape[0]
    while i < m:
        j = i
        while j + 1 < m and a[j + 1] == a[i]:
            j += 1
        t = float(j - i + 1)
        s += t * t * t - t
        i = j + 1
    return s


@njit(parallel=True, cache=True)
def ctrl_tiesums(ctrl_sorted):
    G = ctrl_sorted.shape[0]
    out = np.empty(G, dtype=np.float64)
    for g in prange(G):
        out[g] = _ctrl_tiesum_row(ctrl_sorted[g])
    return out


@njit(cache=True, inline="always")
def _bs_left(a, v):
    lo = 0
    hi = a.shape[0]
    while lo < hi:
        mid = (lo + hi) >> 1
        if a[mid] < v:
            lo = mid + 1
        else:
            hi = mid
    return lo


@njit(cache=True, inline="always")
def _bs_right(a, v):
    lo = 0
    hi = a.shape[0]
    while lo < hi:
        mid = (lo + hi) >> 1
        if a[mid] <= v:
            lo = mid + 1
        else:
            hi = mid
    return lo


@njit(parallel=True, cache=True)
def wilcoxon_z(ctrl_sorted, ctrl_tie, pert):
    """z-scores, scanpy-identical. ctrl_sorted (G,m) sorted rows; pert (G,n)."""
    G, m = ctrl_sorted.shape
    n = pert.shape[1]
    N = n + m
    z = np.zeros(G, dtype=np.float64)
    half = n * (N + 1.0) / 2.0
    base = n * m * (N + 1.0) / 12.0
    denom_tie = float(N) ** 3 - float(N)
    for g in prange(G):
        ps = np.sort(pert[g])
        c = ctrl_sorted[g]
        R1 = 0.0
        tie = ctrl_tie[g]
        i = 0
        while i < n:
            j = i
            v = ps[i]
            while j + 1 < n and ps[j + 1] == v:
                j += 1
            cnt_p = j - i + 1
            lt_c = _bs_left(c, v)
            cnt_c = _bs_right(c, v) - lt_c
            tot = float(cnt_c + cnt_p)
            midrank = (lt_c + i) + (tot + 1.0) / 2.0
            R1 += cnt_p * midrank
            cc = float(cnt_c)
            tie += (tot * tot * tot - tot) - (cc * cc * cc - cc)
            i = j + 1
        tc = 1.0 - tie / denom_tie
        sd = np.sqrt(tc * base)
        if sd > 0.0:
            val = (R1 - half) / sd
            if np.isnan(val):
                val = 0.0
            z[g] = val
        else:
            z[g] = 0.0
    return z

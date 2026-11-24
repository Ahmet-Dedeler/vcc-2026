"""Rank-preserving DE magnitude: the exact-U down-gene pass `dn_k250_zs` (experiments 423 and 429).

Vendored from experiments/423-rank-preserving-magnitude/rp423.py, zs mode only, by experiment 429. The production copy is
src/vcc2026/rank_mag.py (a verbatim copy of experiments/429-zs-production/src_rank_mag.py); the builder pins its sha256.

What the scorer reads from one emitted 400-cell block (carriers + ballast), per DE-universe gene g (reference-control mean
CPM > 5):
  * the Wilcoxon rank sum of v_i = log1p(count_ig / (library_i / 1e6)) against the context's control cells;
  * the arithmetic-mean CPM A_g = mean_i count_ig / (library_i / 1e6), which sets the delivered lfc = log2(A_g / ref_g);
  * pds and mse read only the per-gene block totals (and mse's jackknife).

The K-class carrier emitter delivers a strongly down-designed, low-expressed gene by parking its per-gene total on a few
"mass-holder" carrier cells above the control maximum; every other carrier sits at zero, tied with the control zeros. The
Wilcoxon then sees a strong down gene, but the arithmetic mean reports the bulk total, so the delivered |lfc| is small
(median |l|/|r| 0.46 on H1). On H1 these genes' signs are right 0.65-0.74 of the time (0.56 overall), so more magnitude
on them lowers nmae.

The pass (mode `zs`: zero-sum between library classes). For every selected gene, in order of decreasing |lfc|, counts
move between carrier cells of different library classes: a count in a low-library cell is worth more CPM than a count in
a high-library cell, so for a down gene counts move from low-library carriers to high-library carriers and A_g falls at a
fixed gene total. Every moved cell stays strictly inside its own open interval between adjacent sorted control values
(`cell_bounds`; a cell tied with a control value, e.g. a zero, is frozen; above the control max the interval is open).
Each carrier's library is held fixed: counts added to or removed from a tested gene are paid from, or returned to, that
cell's own sub-gate genes (control CPM <= 5, never tested), pooled across carriers so the pass is zero-sum there too.

Invariants (asserted inside `postpass`, in the scorer's float64 order):
  * every per-gene block total (so pds and mse's profile half are bit-identical), every library, non-negative counts;
  * the rank sum U of every gene (each carrier keeps its exact position among the controls; the ballast is untouched);
    hence the calls and every delivered sign are unchanged. Only the tie term among the block's own cells can move,
    so p moves at the 1e-3 level in log (423: 0 call flips);
  * the ballast rows are bit-identical (zs never touches them), so K = 0 blocks work too;
  * no delivered lfc sign flips (a move that would flip one is skipped).

H1 (423 section 7, production bulk, 3 seeds): nmae -0.00232 +- 0.00005 raw, the other five members bit-identical; target
bootstrap CI [-0.0043, -0.0007], 0 targets worse by more than 0.005. k does not matter above 2 (the rule is capacity-bound:
"all available class-shift capacity on down genes"); the shipped variant is dn_k250_zs. 429 measures it stacked with
426's MDE (VCC_DE_MDE) and 428's zero repair (VCC_DE_ZERO_REPAIR).

Builder use (`VCC_DE_MAG=dn_k250_zs`): per context, `scorer_control_block(control_matrix)` builds the scorer-order control
block (count / (lib / 1e6), as the evaluator ranks; the emitter's rc_control uses count * (1e6 / lib), which differs by
ulps); per block, after ballast_rows and the zero repair, `apply(block, n_carriers, scorer_control, "dn_k250_zs")`.
The pass consumes no RNG, so every other block of the build is bit-identical to a pass-off build.

Variant grammar (zs only): <sel>_k<k*100>_zs, sel in all|q50|q80|up|dn|q80up|q80dn. Example: dn_k250_zs.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

SEL = ("all", "q50", "q80", "up", "dn", "q80up", "q80dn")
BIG = np.int64(1) << 40
VARIANT_RE = re.compile(r"(all|q50|q80|up|dn|q80up|q80dn)_k(\d{3,4})_zs")


# ------------------------------------------------------------------ the scorer-order control block
@dataclass(frozen=True)
class ScorerControl:
    """One context's control cells in the evaluator's own arithmetic (experiments/323 exact_de323.control_block)."""

    keep: np.ndarray       # bool, full gene axis: control mean per-cell CPM > min_cpm (the DE universe)
    cs: np.ndarray         # (G_keep, n_controls) float64, each row sorted: log1p(count / (library / 1e6))
    ref_lin: np.ndarray    # (G_keep,) float64, mean per-cell CPM over the controls (the lfc denominator)
    n_controls: int
    min_cpm: float


def _normalized(block) -> np.ndarray:
    """count / (library / 1e6) in float64, a zero library scaled by 1 (scanpy's order; exact_de323.normalized_counts)."""
    x = np.asarray(block.toarray() if hasattr(block, "toarray") else block, dtype=np.float64)
    scale = x.sum(1) / 1e6
    scale[scale == 0] = 1.
    return x / scale[:, None]


def scorer_control_block(matrix, *, min_cpm: float = 5.0, chunk: int = 1000) -> ScorerControl:
    """The scorer-order control block from a context's control count matrix (CSR or dense), streamed in row chunks.

    Bit-identical to exact_de323.control_block (same chunking, same summation order): on the h1-test controls it must
    reproduce experiments/323-current-family-validation/work/reference.npz's `keep`, `cs` and `ref_lin` exactly
    (test_zs429.py). Peak extra memory: the (G_keep, n) float64 block, ~1.6 GB for 18,400 controls."""
    n, g = matrix.shape
    mean = np.zeros(g, dtype=np.float64)
    for lo in range(0, n, chunk):
        mean += _normalized(matrix[lo:lo + chunk]).sum(0)
    mean /= n
    keep = mean > min_cpm
    cs = np.empty((int(keep.sum()), n), dtype=np.float64)
    for lo in range(0, n, chunk):
        hi = min(lo + chunk, n)
        cs[:, lo:hi] = np.log1p(_normalized(matrix[lo:hi])[:, keep]).T
    cs.sort(axis=1)
    return ScorerControl(keep=keep, cs=cs, ref_lin=mean[keep], n_controls=int(n), min_cpm=float(min_cpm))


def block_lfc(block, cc) -> np.ndarray:
    """The delivered arithmetic-mean lfc on the DE universe, in the evaluator's arithmetic (exact_de323.de_from_counts)."""
    linear = _normalized(block)[:, np.asarray(cc.keep, dtype=bool)]
    return np.log2((linear.mean(0) + 1e-9) / (np.asarray(cc.ref_lin, dtype=np.float64) + 1e-9))


# ------------------------------------------------------------------ variant grammar
def parse_spec(variant: str) -> dict:
    """'dn_k250_zs' -> {'sel': 'dn', 'k': 2.5, 'k_up': None, 'mode': 'zs'}. Only the exact-U zs mode is vendored."""
    m = VARIANT_RE.fullmatch(str(variant).strip())
    if m is None:
        raise ValueError(f"rank_mag variant must look like '<sel>_k<k*100>_zs' (sel in {SEL}; zs mode only): {variant!r}")
    k = int(m.group(2)) / 100.0
    if not 1.0 <= k <= 10.0:
        raise ValueError(f"rank_mag k must be in [1, 10] (k100..k1000): {variant!r}")
    return {"sel": m.group(1), "k": k, "k_up": None, "mode": "zs"}


# ------------------------------------------------------------------ exact per-cell interval arithmetic (rp423, verbatim)
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


def select(sel, l, ok):
    """Genes the pass may move: sel over the finite nonzero-lfc DE universe (q50/q80 = |lfc| quantiles)."""
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


def _subgate_zs(X, sidx, sub_delta):
    """Realise the carriers' sub-gate library compensation for a zero-sum pass (rp423._subgate without the ballast
    legs, which zs never reaches). Carriers that gained kept counts give up sub-gate counts (largest holdings first)
    into a pool; carriers that lost kept counts are refilled from that pool in cell order. The pass is zero-sum per
    gene, so the pool is exactly used up."""
    pool = np.zeros(sidx.size, dtype=np.int64)
    for i in np.flatnonzero(sub_delta < 0):
        t = _take(X[i, sidx], int(-sub_delta[i]))
        X[i, sidx] -= t
        pool += t
    for i in np.flatnonzero(sub_delta > 0):
        d = int(sub_delta[i])
        if d > int(pool.sum()):
            raise AssertionError("zero-sum pass: the sub-gate pool ran short")
        g = _take(pool, d)
        pool -= g
        X[i, sidx] += g
    if pool.any():
        raise AssertionError("zero-sum pass: sub-gate counts left in the pool")


# ------------------------------------------------------------------ the pass
def postpass(block, n_car, cc, lfc_base, spec, *, sub_margin=0):
    """Apply the exact-U zero-sum magnitude pass to one block (rows [0, n_car) carriers, the rest ballast).

    `cc` needs .keep (full-axis bool), .cs (sorted scorer-order control values per kept gene) and .ref_lin; `lfc_base` is
    the block's delivered lfc on the kept axis (block_lfc). Returns (new block, diagnostics); the diagnostics carry
    `_fr_by_gene` (touched genes, achieved fraction, lfc, target dA) and `_touched` (moved-gene mask) for analysis.
    Bit-identical to rp423.postpass(..., mode 'zs') on every block (test_zs429.py)."""
    if spec.get("mode") != "zs":
        raise ValueError(f"only the zs mode is vendored: {spec!r}")
    if spec.get("k_up") is not None:
        raise ValueError("a separate up factor is not vendored")
    X = np.array(block, dtype=np.int64, copy=True)
    n = X.shape[0]
    n_car = int(n_car)
    if not 2 <= n_car <= n:
        raise ValueError(f"n_car {n_car} outside [2, {n}]")
    keep = np.asarray(cc.keep, dtype=bool)
    if keep.size != X.shape[1]:
        raise ValueError(f"control gate has {keep.size} genes, the block {X.shape[1]}")
    gidx = np.flatnonzero(keep)
    sidx = np.flatnonzero(~keep)
    Lr = X.sum(1)
    if np.any(Lr[:n_car] <= 0):
        raise AssertionError("a carrier cell has an empty library")
    s = Lr / 1e6                                     # the scorer's per-cell scale (float64 of an exact int)
    car = np.arange(n_car)
    ref = np.asarray(cc.ref_lin, dtype=np.float64)
    Xk = X[:, gidx]
    A = (Xk / s[:, None]).mean(0)
    l = np.asarray(lfc_base, dtype=np.float64)
    l_chk = np.log2((A + 1e-9) / (ref + 1e-9))
    fin = np.isfinite(l) & np.isfinite(l_chk)
    if np.max(np.abs(l[fin] - l_chk[fin]), initial=0.0) > 1e-9:
        raise AssertionError("base lfc does not match the block")
    ok = fin & (A > 0) & (l != 0)
    k = float(spec["k"])
    touch = select(spec["sel"], l, ok)
    kvec = np.full(l.shape, k)
    A_t = (ref + 1e-9) * np.power(2.0, kvec * l) - 1e-9
    dA = np.where(touch, A_t - A, 0.0)
    # priority: largest |lfc| first (they carry the nmae and get the sub-gate budget first)
    order = np.flatnonzero(touch)[np.argsort(-np.abs(l[touch]), kind="stable")]
    sub_avail = X[np.ix_(car, sidx)].sum(1).astype(np.int64) - int(sub_margin)   # carrier sub-gate budget
    sub_delta = np.zeros(n_car, dtype=np.int64)       # net sub-gate change per carrier cell (+ = gained)
    cw = 1.0 / (n * s[car])                           # CPM of the 400-cell mean per count in each carrier cell
    groups, gw = _levels(cw)
    G_ = gidx.size
    ach = np.zeros(G_)
    dU_car = np.zeros(G_, dtype=np.int64)
    moved = np.zeros(G_, dtype=np.int64)
    sign_guard = np.zeros(G_, dtype=bool)
    for j in order:
        need = dA[j]
        if need == 0:
            continue
        c = cc.cs[j]
        x0 = Xk[car, j].copy()
        room0 = np.maximum(sub_avail + sub_delta, 0)
        xl, xh, *_ = cell_bounds(x0, s[car], c)
        x = _zs_step(x0, xl, xh, room0, cw, need, groups, gw)
        if np.array_equal(x, x0):
            continue
        d_car = int((ranks_of(x, s[car], c) - ranks_of(x0, s[car], c)).sum())
        if d_car != 0:
            raise AssertionError(f"mode zs moved U by {d_car}")
        if int(x.sum() - x0.sum()) != 0:
            raise AssertionError("gene total moved")
        dx = x - x0
        A_new = A[j] + (dx * cw).sum()
        if np.sign(np.log2((A_new + 1e-9) / (ref[j] + 1e-9))) != np.sign(l[j]):   # a magnitude pass never flips a sign
            sign_guard[j] = True
            continue
        Xk[car, j] = x
        sub_delta -= dx
        if np.any(sub_avail + sub_delta < 0):
            raise AssertionError("sub-gate budget overdrawn")
        dU_car[j] = d_car
        moved[j] = int(np.abs(dx).sum())
        ach[j] = (dx * cw).sum()
    # ---- write the touched genes back, then realise the sub-gate compensation (libraries fixed)
    X[:, gidx] = Xk
    _subgate_zs(X, sidx, sub_delta)
    B = np.asarray(block)
    if not np.array_equal(X.sum(1), Lr):
        raise AssertionError("a library moved")
    if not np.array_equal(X.sum(0), B.sum(0)):
        raise AssertionError("a per-gene block total moved")
    if not np.array_equal(X[n_car:], B[n_car:]):
        raise AssertionError("a ballast row moved")
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
                frac_full=float(np.mean(fr >= 0.95)) if t.size else 0.0,
                frac_zero=float(np.mean(np.abs(fr) <= 0.05)) if t.size else 0.0,
                sign_guard=int(sign_guard.sum()), dU_car_abs_sum=int(np.abs(dU_car).sum()),
                sub_used=int(-sub_delta[sub_delta < 0].sum()),
                sub_budget_exhausted_cells=int(((sub_avail + sub_delta) <= 0).sum()))
    diag["_fr_by_gene"] = (t, fr, l[t], dA[t])
    diag["_touched"] = moved > 0
    return X, diag


def apply(block, n_car, cc, variant="dn_k250_zs"):
    """Builder entry point: the pass on one final block with its lfc computed in the evaluator's arithmetic.

    Returns (new int64 block, diagnostics). The diagnostics are flat and numeric (prefix `zs_`), as the builder's
    rc_diagnostics_summary needs; None values (e.g. the up-gene fraction of a down-only rule) are dropped."""
    spec = parse_spec(variant) if isinstance(variant, str) else dict(variant)
    X = np.asarray(block)
    out, diag = postpass(X, n_car, cc, block_lfc(X, cc), spec)
    diag.pop("_fr_by_gene")
    diag.pop("_touched")
    return out, {f"zs_{k}": float(v) for k, v in diag.items() if isinstance(v, (int, float, np.integer, np.floating))
                 and not isinstance(v, bool)}

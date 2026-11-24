#!/usr/bin/env python3
"""630 — the rank-neutral UP magnitude pass (council #1), applied after the production zs pass, scored on H1.

Why: nmae reads per-gene ARITHMETIC-mean CPM; calls, fid, jac and reach read Wilcoxon RANKS; pds and mse read per-gene
block totals (+ mse's jackknife). For an uncalled gene we sign UP, the L1-optimal delivered |lfc| is larger than today's
MDE-sized add whenever signs are good enough (632, 603, 610), but raising it by the emitter's own route leaks calls
(426/435). This pass raises such a gene's mean WITHOUT moving any rank or any per-gene total:
  * recipients: carrier cells whose value for the gene is already strictly above every control value (open interval:
    any larger count keeps the rank sum U unchanged);
  * donors: the ballast rows (deep cells), each staying inside its own interval between adjacent control values
    (rank_mag.cell_bounds), so U is unchanged there too;
  * libraries are held fixed and every per-gene total is preserved: for each count of gene j moved ballast -> carrier,
    the carrier gives one count of a sub-gate gene h (control CPM <= 5, never tested) back to the SAME ballast row's h.
  So per-gene block totals (pds bit-identical; mse's profile half), every library, and every kept gene's U are unchanged;
  only the tie term among the block's own cells can move p (as for zs, 423). Calls are verified unchanged after the pass.
Target per selected gene j: t_j = max(l_j, K * 0.24 / sqrt(ref_count_j)) (603's rule: MDE-scaled), with ref_count_j the
control mean count per cell. Selection: uncalled (base q >= 0.05), finite, not the target gene, UP by the sign source.

usage: uv run python experiments/630-rank-neutral-up/rnup630.py <seed> <variant> [<variant> ...]
  variant = up_k<K*100>[_<signsource>]   signsource: dlv (delivered lfc > 0; default)
Reads work/blocks/base_s<seed>_..._{i:02d}.npz (run630.py's dump), writes work/h1/stats-rnup_s<seed>_<variant>.npz,
metrics json, diag json. Known answer first: the zs pass re-applied to each dumped block must reproduce run630's own
zs stats (sums, jk, p, q, lfc) bit for bit, else abort.
"""
from __future__ import annotations

import glob
import json
import re
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/400-uncalled-knn"))
sys.path.insert(0, str(ROOT / "experiments/410-pertscale-revisit"))
sys.path.insert(0, str(ROOT / "src"))
import run400 as H  # noqa: E402  (read-only: R, E = the 323 reference + evaluator)
from vcc2026 import rank_mag as ZS  # noqa: E402  (production zs pass)

R, E = H.R, H.E
WORK = HERE / "work" / "h1"
TAG = "base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse"
BIG = ZS.BIG


def _min_above(cmax, s_a):
    """Smallest integer count x with log1p(x / s_a) > cmax (strictly above every control value)."""
    x = int(np.floor(np.expm1(cmax) * s_a)) + 1
    while x > 0 and np.log1p((x - 1) / s_a) > cmax:
        x -= 1
    while not np.log1p(x / s_a) > cmax:
        x += 1
    return x


def up_pass(block, n_car, cc, l_base, q_base, target_lfc, sel, promote=False, relax=False, max_rec=1, order_by="dA"):
    """One block. Returns (new block, diag). sel: kept-axis bool of genes to raise; target_lfc: kept-axis desired lfc.
    promote=False: strict (U unchanged for every kept gene; recipients must already sit above the control max).
    promote=True: if no carrier is above the control max, the highest carrier is first raised to just above it (one cell
    per gene; that gene's U rises by the controls it passes, recorded and asserted exactly; dAUC = dU / (2 n1 n2)).
    relax=True: ballast donors may give ALL their counts of the gene (crossing control values downward); each ballast
    cell is 1 of n cells, so its U effect is bounded by n_ctrl per cell (|dAUC| <= 1/n per ballast row); every U change is
    recorded and asserted exactly, and calls are verified unchanged by the caller.
    max_rec: with promote, up to this many carriers per gene are raised above the control max (the highest-valued carriers
    that still have sub-gate budget), each promotion's U change recorded. order_by: 'dA' (largest CPM ask first, the
    first version) or 'lowA' (lowest current mean first: most lfc per delivered count, since d log2 A / dA = 1 / (A ln 2))."""
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
    car, bal = np.arange(n_car), np.arange(n_car, n)
    cw = 1.0 / (n * s)                                   # CPM of the n-cell mean per count, per cell
    A_t = (ref + 1e-9) * np.power(2.0, target_lfc) - 1e-9
    dA = np.where(sel & np.isfinite(target_lfc), np.maximum(A_t - A, 0.0), 0.0)
    pos = np.flatnonzero(dA > 0)
    order = pos[np.argsort(-dA[pos], kind="stable")] if order_by == "dA" else pos[np.argsort(A[pos], kind="stable")]
    sub_left = X[np.ix_(car, sidx)].sum(1).astype(np.int64)   # each carrier's sub-gate budget
    sub_moves = []                                        # (carrier, ballast, n counts)
    stats = dict(n_sel=int(order.size), n_no_recipient=0, n_no_donor=0, n_moved=0, counts=0, need=float(dA.sum()),
                 got=0.0, n_budget_cut=0, n_promoted=0)
    dU = np.zeros(gidx.size, dtype=np.int64)            # recorded U change per kept gene (promotions only)
    for j in order:
        c = cc.cs[j]
        cmax = c[-1]
        vcar = np.log1p(Xk[car, j] / s[car])
        rec = np.flatnonzero(vcar > cmax)                 # strictly above every control value: open interval
        xb = Xk[bal, j]
        xl, xh, *_ = ZS.cell_bounds(xb, s[bal], c)
        dcap = xb.copy() if relax else xb - xl
        u_bal0 = int(ZS.ranks_of(xb, s[bal], c).sum()) if relax else 0
        if dcap.sum() <= 0:
            stats["n_no_donor"] += 1
            continue
        rem = dA[j]
        moved_j = 0
        if rec.size < max_rec and promote:
            # promote the highest-valued carriers that still have budget (ties: lowest library) to just above the max
            cand = np.setdiff1d(np.flatnonzero(sub_left > 0), rec)
            cand = cand[np.lexsort((-cw[car][cand], -vcar[cand]))]
            for a in cand[: max_rec - rec.size]:
                a = int(a)
                k = _min_above(cmax, s[car][a]) - int(Xk[car[a], j])
                if k <= 0 or k > int(sub_left[a]) or k > int(dcap.sum()):
                    continue
                u_before = int(ZS.ranks_of(Xk[car[a]:car[a] + 1, j], s[car[a]:car[a] + 1], c).sum())
                take = ZS.fill(dcap, k)
                for bi in np.flatnonzero(take):
                    t = int(take[bi])
                    Xk[bal[bi], j] -= t
                    Xk[car[a], j] += t
                    sub_moves.append((int(car[a]), int(bal[bi]), t))
                    rem -= t * (cw[car][a] - cw[bal][bi])
                dcap -= take
                sub_left[a] -= k
                moved_j += k
                u_after = int(ZS.ranks_of(Xk[car[a]:car[a] + 1, j], s[car[a]:car[a] + 1], c).sum())
                dU[j] += u_after - u_before
                stats["n_promoted"] += 1
                rec = np.append(rec, a)
        if rec.size == 0:
            stats["n_no_recipient"] += 1
            continue
        rec = rec[np.argsort(-cw[car][rec], kind="stable")]   # lowest library first (most CPM per count)
        for rc in rec:
            if rem <= 0:
                break
            budget = int(sub_left[rc])
            if budget <= 0:
                continue
            gain = cw[car][rc] - float(np.mean(cw[bal]))      # A per count, approximately (ballast weight is tiny)
            want = int(np.ceil(rem / max(gain, 1e-300)))
            k = min(want, budget, int(dcap.sum()))
            if k <= 0:
                stats["n_budget_cut"] += 1
                continue
            take = ZS.fill(dcap, k)                       # spread the donation over ballast rows by capacity
            for bi in np.flatnonzero(take):
                t = int(take[bi])
                Xk[bal[bi], j] -= t
                Xk[car[rc], j] += t
                sub_moves.append((int(car[rc]), int(bal[bi]), t))
                rem -= t * (cw[car][rc] - cw[bal][bi])
            dcap -= take
            sub_left[rc] -= k
            moved_j += k
        if relax:
            dU[j] += int(ZS.ranks_of(Xk[bal, j], s[bal], c).sum()) - u_bal0
        if moved_j:
            stats["n_moved"] += 1
            stats["counts"] += moved_j
            stats["got"] += float(dA[j] - rem)           # rem < 0 means a small overshoot (whole counts)
    X[:, gidx] = Xk
    # library compensation with matched sub-gate genes: carrier -> the same gene on the donor ballast row
    for rc, bb, t in sub_moves:
        row = X[rc, sidx]
        give = ZS._take(row, t)
        X[rc, sidx] -= give
        X[bb, sidx] += give
    # ---- invariants
    B = np.asarray(block)
    if not np.array_equal(X.sum(1), Lr):
        raise AssertionError("a library moved")
    if not np.array_equal(X.sum(0), B.sum(0)):
        raise AssertionError("a per-gene block total moved")
    if X.min() < 0:
        raise AssertionError("negative count")
    touched = np.flatnonzero(np.any(X[:, gidx] != B[:, gidx], axis=0))
    for j in touched:                                     # exact U check on every touched kept gene
        c = cc.cs[j]
        u0 = int(ZS.ranks_of(B[:, gidx[j]], s, c).sum())
        u1 = int(ZS.ranks_of(X[:, gidx[j]], s, c).sum())
        if u1 - u0 != int(dU[j]):
            raise AssertionError(f"U moved on kept gene {j}: {u0} -> {u1} (recorded promotion {int(dU[j])})")
    stats["n_touched_kept"] = int(touched.size)
    n_ctrl = cc.cs.shape[1]
    stats["max_dAUC"] = float(np.abs(dU).max() / (2.0 * n * n_ctrl)) if dU.size else 0.0   # ranks_of is in x2 units
    stats["net_dU_sum"] = int(dU.sum())
    stats["mean_dAUC_promoted"] = float(dU[dU > 0].mean() / (2.0 * n * n_ctrl)) if (dU > 0).any() else 0.0
    stats["frac_got"] = stats["got"] / max(stats["need"], 1e-300)
    return X, stats


def _max_drop_k(x, s_b, c, budget):
    """Largest k in [0, x] such that removing k counts from one cell (scale s_b) lowers its ranks_of value by <= budget."""
    r0 = int(ZS.ranks_of(np.array([x]), np.array([s_b]), c)[0])
    lo, hi = 0, int(x)
    while lo < hi:
        mid = (lo + hi + 1) // 2
        drop = r0 - int(ZS.ranks_of(np.array([x - mid]), np.array([s_b]), c)[0])
        if drop <= budget:
            lo = mid
        else:
            hi = mid - 1
    return lo


def up_pass2(block, n_car, cc, l_base, target_lfc, sel, *, max_dauc=1e-4, cross=False, exclude=None, bounds=None,
             quantum=None, rounds=1, rec_pool=1):
    """v2 (the coordinator's bar, 2026-10-02): per gene, one carrier is raised to just above the control max (if none is
    there yet) and counts of the gene move ballast -> that carrier, with matched sub-gate library compensation (per-gene
    totals and libraries exact). Each gene's NET U change (promotion up, ballast crossings down) is recorded exactly and
    kept within |dAUC| <= max_dauc; a gene that cannot meet it is reverted. cross=False: ballast donors stay inside their
    control intervals; cross=True: they may cross controls within the gene's remaining U budget (the promotion's upward
    move is offset by the donors' downward crossings). exclude: kept-axis indices never touched (the call-flip guard).
    bounds: optional (lo, hi) per kept gene in ranks_of x2 units (informational margin variant: each gene's net U change
    must stay in [lo_j, hi_j] instead of [-capU, +capU]).
    v3 load-spreading (defaults reproduce v2 exactly): quantum = max counts a gene receives per round (water-filling:
    lfc is concave in counts, log2(1 + k cw / A)); rounds = passes over the genes (a gene keeps its recipient across
    rounds); rec_pool = choose the recipient among the rec_pool highest-valued carriers by the most remaining sub-gate
    budget (spreads the load off the few shared top carriers)."""
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
    capU = int(np.floor(max_dauc * 2.0 * n * n_ctrl))     # ranks_of works in x2 units
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
    st = dict(n_sel=int(order.size), n_moved=0, n_reverted=0, n_no_candidate=0, counts=0, need=float(dA.sum()), got=0.0,
              n_no_pos_donor=0, n_neg_donor_rows=0, n_mean_revert=0)
    u_base = {}                                          # per gene: ranks_of of the ORIGINAL column (net U across rounds)
    rec_of = {}
    moved_any = set()
    rem_of = {int(j): float(dA[j]) for j in order}
    for rnd in range(int(rounds)):
      for j in order:
        j = int(j)
        if rem_of[j] <= 0:
            continue
        c = cc.cs[j]
        cmax = c[-1]
        snap_col, snap_sub, snap_mv = Xk[:, j].copy(), sub_left.copy(), len(sub_moves)
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
            a = int(cand[np.argmax(sub_left[cand])]) if rec_pool > 1 else int(cand[0])
            k_promo = _min_above(cmax, s[a]) - int(Xk[a, j])
            dU_promo = 2 * n_ctrl - int(u0_cells[a])
        # 605 review: only donors whose library is LARGER than the carrier's raise the mean (gain per count
        # cw[a] - cw[b] > 0); a move from a shallower donor lowers it. Production ballast rows are deep (~995k vs ~20k
        # carriers), so this never binds there; it keeps the pass correct for any block.
        elig = cw[bal] < cw[a]
        st["n_neg_donor_rows"] += int((~elig).sum())
        if not elig.any():
            st["n_no_candidate"] += 1
            st["n_no_pos_donor"] += 1
            continue
        # donor capacity per ballast row
        lo_j, hi_j = (-capU, capU) if bounds is None else (int(bounds[0][j]), int(bounds[1][j]))
        du_sofar = int((u0_cells - u_base[j]).sum())      # this gene's net U change from earlier rounds
        lo_j, hi_j = lo_j - du_sofar, hi_j - du_sofar
        xb = Xk[bal, j].copy()
        if cross:
            budget = dU_promo - lo_j                   # total downward ballast movement allowed (net stays >= lo_j)
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
        if k_promo > int(sub_left[a]) or k_promo > int(dcap.sum()):
            st["n_no_candidate"] += 1
            continue
        rem = rem_of[j]
        k_want = max(k_promo, int(np.ceil(rem / (cw[a] - float(np.min(cw[bal][elig]))))))
        if quantum is not None:
            k_want = min(k_want, max(k_promo, int(quantum)))
        k = min(k_want, int(sub_left[a]), int(dcap.sum()))
        if k <= 0:
            st["n_no_candidate"] += 1
            continue
        take = ZS.fill(dcap, k)
        for bi in np.flatnonzero(take):
            t = int(take[bi])
            Xk[bal[bi], j] -= t
            Xk[a, j] += t
            sub_moves.append((a, int(bal[bi]), t))
            rem -= t * (cw[a] - cw[bal[bi]])
        sub_left[a] -= k
        du = int((ZS.ranks_of(Xk[:, j], s, c) - u0_cells).sum())
        a_old, a_new = float((snap_col / s).mean()), float((Xk[:, j] / s).mean())
        sign_new = np.sign(np.log2((a_new + 1e-9) / (ref[j] + 1e-9)))
        bad_mean = (a_new < a_old) or (sign_new != np.sign(l_base[j]))   # 605: the mean may never fall or flip sign
        if bad_mean:
            st["n_mean_revert"] += 1
        if bad_mean or du < lo_j or du > hi_j or not (np.log1p(Xk[a, j] / s[a]) > cmax):
            Xk[:, j] = snap_col                         # revert this gene
            sub_left[:] = snap_sub
            del sub_moves[snap_mv:]
            st["n_reverted"] += 1
            continue
        dU[j] = int((ZS.ranks_of(Xk[:, j], s, c) - u_base[j]).sum())
        rec_of[j] = a
        if j not in moved_any:
            moved_any.add(j)
            st["n_moved"] += 1
        st["counts"] += k
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
    au = np.abs(dU[dU != 0]) / (2.0 * n * n_ctrl)
    st.update(max_dAUC=float(au.max()) if au.size else 0.0, median_dAUC=float(np.median(au)) if au.size else 0.0,
              frac_got=st["got"] / max(st["need"], 1e-300), capU=capU)
    return X, st


def load_block(path):
    import scipy.sparse as sp
    z = np.load(path, allow_pickle=False)
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
    return X, int(z["n_car"]), str(z["target"])


def main():
    args = sys.argv[1:]
    src_dir, name = HERE / "work", ""
    if "--src" in args:                      # stack on another construction's dumped blocks (e.g. 612's idkgd)
        i = args.index("--src"); src_dir = Path(args[i + 1]); del args[i:i + 2]
    if "--name" in args:
        i = args.index("--name"); name = args[i + 1] + "_"; del args[i:i + 2]
    global TAG
    if "--tag" in args:                      # another construction's block / stats tag (e.g. 654's x125 zr dumps)
        i = args.index("--tag"); TAG = args[i + 1]; del args[i:i + 2]
    seed = int(args[0])
    variants = args[1:]
    src_h1 = src_dir if (src_dir / f"stats-{TAG.format(seed=seed)}_dnk250zs.npz").exists() else src_dir / "h1"
    r, c, emit_c = E.load_reference()
    keep = np.asarray(c.keep, dtype=bool)
    G = r["genes"][keep]
    lib_med = float(np.median(r["libraries"]))
    ref_count = np.asarray(c.ref_lin, float) * lib_med / 1e6
    mde = 0.24 / np.sqrt(np.maximum(ref_count, 1e-12))
    tag = TAG.format(seed=seed)
    zs_ref = np.load(src_h1 / f"stats-{tag}_dnk250zs.npz", allow_pickle=False)
    # numeric order of the dump index ({tag}_{i:02d}.npz): lexicographic order breaks at >= 100 targets ([21dd4c], 618)
    files = sorted(glob.glob(str(src_dir / f"blocks/{tag}_*.npz")),
                   key=lambda f: int(re.search(r"_(\d+)\.npz$", f).group(1)))
    assert len(files) == len(r["targets"]), (len(files), len(r["targets"]))   # one block per reference target (H1: 85)
    spec = ZS.parse_spec("dn_k250_zs")
    acc = {v: dict(sums=[], jks=[], ps=[], qs=[], lfcs=[], diags=[]) for v in ["zsKA"] + variants}
    t0 = time.monotonic()
    for i, f in enumerate(files):
        X0, n_car, target = load_block(f)
        if str(target) != str(r["targets"][i]):
            raise AssertionError(f"block {f} holds target {target}, reference row {i} is {r['targets'][i]}")
        s0, j0, p0, q0, l0 = E.block_stats(X0, c)
        Xz, _ = ZS.postpass(X0, n_car, c, l0, spec)
        sz, jz, pz, qz, lz = E.block_stats(Xz, c)
        # KA: the re-applied zs pass reproduces run630's own zs stats for this target
        for nm, a, b in (("sums", sz, zs_ref["sums"][i]), ("p", pz, zs_ref["p"][i]), ("q", qz, zs_ref["q"][i]),
                         ("lfc", lz, zs_ref["lfc"][i])):
            if not np.array_equal(np.asarray(a), np.asarray(b), equal_nan=True):
                raise AssertionError(f"KA FAIL target {i} {target}: re-applied zs {nm} differs from run630's zs stats")
        if jz != float(zs_ref["jk"][i]):
            raise AssertionError(f"KA FAIL target {i} jk")
        for nm, val in (("sums", sz), ("jks", jz), ("ps", pz), ("qs", qz), ("lfcs", lz)):
            acc["zsKA"][nm].append(val)
        own = G == target
        uncalled = (qz >= 0.05) & np.isfinite(lz) & ~own
        for v in variants:
            m2 = re.fullmatch(r"(s|b|m)_k(\d+)(?:_q(\d+)r(\d+)p(\d+))?(?:_f([azr]*))?(?:_d(\d+)(c?))?(_L)?(?:_D(\d+))?(?:_P(\d+))?", v)   # v2: s/b/m; v3 suffix _q<quantum>r<rounds>p<pool>; v4 flags _f<azr> (up4_630.py); v5 DOWN currency _d<mult*10>[c] (up5_630.py)
            if m2:
                bnd = None
                if m2.group(1) == "m":                   # informational: per-gene z-margin bounds (not the shipped bar)
                    from scipy.stats import norm as _norm
                    sys.path.insert(0, str(ROOT / "experiments/323-current-family-validation"))
                    import exact_de323 as _X
                    zz = _X.de_from_counts(np.asarray(Xz, dtype=np.float64), c)[0]
                    n1, n2 = 400, c.cs.shape[1]
                    sdU = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)          # U sd (no tie term)
                    m_t = int(np.isfinite(pz).sum())
                    p_star = 0.05 * (int((qz < 0.05).sum()) + 1) / max(m_t, 1)
                    zc = float(_norm.isf(p_star / 2))
                    zmax = max(zc - 0.5, 0.0)
                    capm = int(np.floor(1e-3 * 2.0 * n1 * n2))
                    lo = np.maximum(np.floor((-zmax - zz) * sdU * 2.0), -capm)
                    hi = np.minimum(np.floor((zmax - zz) * sdU * 2.0), capm)
                    bnd = (np.minimum(lo, 0).astype(np.int64), np.maximum(hi, 0).astype(np.int64))
                K = int(m2.group(2)) / 100.0
                # K = 0 means NO UP pass (631's DN-alone arm): without this, float-noise asks (A_t == A up to rounding)
                # promoted ~425 carriers per target (found by test_dn_mag631's dn-mode KA, 2026-10-02)
                sel = (uncalled & (lz > 0)) if K > 0 else np.zeros_like(uncalled)
                t_lfc = np.maximum(lz, K * mde)
                excl, guard_rounds = set(), 0
                while True:                              # call-flip guard: re-run without any gene whose call flipped
                    v3 = dict(quantum=int(m2.group(3)), rounds=int(m2.group(4)), rec_pool=int(m2.group(5))) if m2.group(3) else {}
                    if m2.group(11):                     # v7: UP receipts paid with budgeted called-DOWN counts (up7_630.py)
                        import up7_630 as U7
                        Xu, dg = U7.up_pass7(Xz, n_car, c, lz, t_lfc, sel, max_dauc=1e-4, cross=m2.group(1) in ("b", "m"),
                                             exclude=excl, bounds=bnd, **v3, dn_mult=int(m2.group(11)) / 100.0,
                                             dn_mask=(qz < 0.05))
                    elif m2.group(9):                    # v6 PROBE: library drift as the compensation (up6_630.py)
                        import up6_630 as U6
                        Xu, dg = U6.up_pass6(Xz, n_car, c, lz, t_lfc, sel, max_dauc=1e-4, cross=m2.group(1) in ("b", "m"),
                                             exclude=excl, bounds=bnd, **v3)
                        _da = U6.all_gene_dauc(Xz, Xu, c)
                        dg.update(L_max_dAUC_all=float(_da.max()), L_n_over_1e4=int((_da > 1e-4).sum()),
                                  L_median_dAUC_touched=float(np.median(_da[_da > 0])) if (_da > 0).any() else 0.0,
                                  L_lib_change_max=float(np.max(np.abs(Xu.sum(1) / Xz.sum(1) - 1.0))))
                    elif m2.group(7) is not None:        # v5 DOWN-gene compensation currency (driver-only harness)
                        assert m2.group(6) is None, "v4 flags and v5 currency are not combined"
                        import up5_630 as U5
                        Xu, dg = U5.up_pass5(Xz, n_car, c, lz, t_lfc, sel, max_dauc=1e-4, cross=m2.group(1) in ("b", "m"),
                                             exclude=excl, bounds=bnd, **v3, dn_mult=int(m2.group(7)) / 10.0,
                                             dn_mask=(qz < 0.05) if m2.group(8) else None)
                    elif m2.group(6) is not None:        # v4 capacity flags (driver-only harness; the module is untouched)
                        import up4_630 as U4
                        Xu, dg = U4.up_pass4(Xz, n_car, c, lz, t_lfc, sel, max_dauc=1e-4, cross=m2.group(1) in ("b", "m"),
                                             exclude=excl, bounds=bnd, **v3, **U4.flag_kwargs(m2.group(6)))
                    else:
                        Xu, dg = up_pass2(Xz, n_car, c, lz, t_lfc, sel, max_dauc=1e-4, cross=m2.group(1) in ("b", "m"), exclude=excl,
                                          bounds=bnd, **v3)
                    if m2.group(10):                     # 631 DOWN pass after the UP pass (driver-only prototype)
                        sys.path.insert(0, str(ROOT / "experiments/631-nmae-headroom"))
                        import dn631 as DN
                        # shipping semantics (src_dn_mag631): selection from the stats of the block it is given, own gene out
                        _, _, _pD, _qD, _lD = E.block_stats(Xu, c)
                        Xu, dgd = DN.dn_pass(Xu, n_car, c, _lD, _qD, int(m2.group(10)) / 100.0, max_dauc=1e-4, exclude=excl,
                                             sel_extra=~own)
                        dg.update(dgd)
                    su, ju, pu, qu, lu = E.block_stats(Xu, c)
                    flips = set(np.flatnonzero((qz < 0.05) != (qu < 0.05)).tolist())
                    if not flips or guard_rounds >= 5:
                        break
                    excl |= flips
                    guard_rounds += 1
                dg.update(guard_excluded=len(excl), guard_rounds=guard_rounds)
            else:
                m = re.fullmatch(r"up_k(\d+)(?:_p(\d*))?(_r)?(?:_(lowA))?(?:_(dlv))?", v)
                K = int(m.group(1)) / 100.0
                sel = uncalled & (lz > 0)
                t_lfc = np.maximum(lz, K * mde)
                promote = m.group(2) is not None
                Xu, dg = up_pass(Xz, n_car, c, lz, qz, t_lfc, sel, promote=promote, relax=bool(m.group(3)),
                                 max_rec=int(m.group(2) or 1) if promote else 1, order_by="lowA" if m.group(4) else "dA")
                su, ju, pu, qu, lu = E.block_stats(Xu, c)
            calls0, calls1 = qz < 0.05, qu < 0.05
            dg.update(target=target, sums_equal=bool(np.array_equal(su, sz)), jk_delta=float(ju - jz),
                      calls_changed=int((calls0 != calls1).sum()), max_abs_dlogp=float(np.nanmax(np.abs(np.log10(np.maximum(pu, 1e-300)) - np.log10(np.maximum(pz, 1e-300))))))
            if dg["calls_changed"]:
                print(f"WARNING {target} {v}: {dg['calls_changed']} calls changed", flush=True)
            for nm, val in (("sums", su), ("jks", ju), ("ps", pu), ("qs", qu), ("lfcs", lu)):
                acc[v][nm].append(val)
            acc[v]["diags"].append(dg)
        if (i + 1) % 10 == 0:
            print(f"s{seed} {i + 1}/85 ({time.monotonic() - t0:.0f}s)", flush=True)
    out = {}
    for v, a in acc.items():
        raw, rows, legs = E.score(r, c, a["sums"], a["jks"], a["ps"], a["qs"], a["lfcs"])
        out[v] = raw
        if v != "zsKA":
            np.savez(WORK / f"stats-rnup_{name}s{seed}_{v}.npz", targets=r["targets"], sums=np.asarray(a["sums"]),
                     jk=np.asarray(a["jks"]), p=np.asarray(a["ps"]), q=np.asarray(a["qs"]), lfc=np.asarray(a["lfcs"]))
            d = a["diags"]
            keys = [k for k in ("n_sel", "n_no_recipient", "n_no_donor", "n_moved", "n_promoted", "n_reverted", "n_no_candidate",
                                "counts", "frac_got", "calls_changed", "jk_delta", "max_dAUC", "median_dAUC", "max_abs_dlogp",
                                "guard_excluded", "dn_paid", "sub_paid", "dn_genes_used", "n_pay_short", "dn_pool", "L_max_dAUC_all",
                                "L_n_over_1e4", "L_median_dAUC_touched", "L_lib_change_max", "dn_sel", "dn_moved", "dn_counts",
                                "dn_frac") if k in d[0]]
            summ = {k: float(np.mean([x[k] for x in d])) for k in keys}
            summ["max_dAUC_all"] = float(max(x["max_dAUC"] for x in d))
            summ["max_abs_dlogp_all"] = float(max(x["max_abs_dlogp"] for x in d))
            summ["all_sums_equal"] = all(x["sums_equal"] for x in d)
            summ["total_calls_changed"] = int(sum(x["calls_changed"] for x in d))
            out[v + "_summary"] = summ
            (WORK / f"diag-rnup_{name}s{seed}_{v}.json").write_text(json.dumps(d, indent=1, default=float))
    ka = json.loads((src_h1 / f"metrics-{tag}_dnk250zs.json").read_text())["raw"] if (src_h1 / f"metrics-{tag}_dnk250zs.json").exists() else None
    out["zsKA_matches_run630_metrics"] = (None if ka is None else all(abs(out["zsKA"][k] - ka[k]) < 1e-12 for k in ka))
    (WORK / f"rnup630_{name}s{seed}_{'+'.join(variants)}.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1, default=float), flush=True)


if __name__ == "__main__":
    sys.exit(main())

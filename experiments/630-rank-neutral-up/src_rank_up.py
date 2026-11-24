"""Rank-neutral UP magnitude pass (experiment 630; candidate for src/vcc2026/rank_up.py, builder knob VCC_DE_UPRN).

What the scorer reads from one emitted block (carriers + ballast), per DE-universe gene g (control mean CPM > 5):
  * Wilcoxon rank sums of log1p(count / (library / 1e6)) against the context's control cells -> calls, fid, jac, reach;
  * arithmetic-mean CPM -> the delivered lfc -> nmae;
  * per-gene block totals (+ mse's jackknife) -> pds, mse.
For uncalled genes the design signs UP, the L1-optimal delivered |lfc| is above today's MDE-sized add (experiments 603,
610, 632), but raising it through the emitter leaks calls (426/435). This pass raises such a gene's arithmetic mean after
the zs pass while keeping, for every gene, the block total and every library exact, and every gene's call unchanged:
  * one carrier per gene is raised to just above the control maximum (if none is there yet) and receives counts of the
    gene from the ballast rows (deep cells), so the mean CPM rises by ~count / (n * carrier scale);
  * each count the carrier receives is paid by one count of a sub-gate gene (control CPM <= 5, never tested) that moves to
    the SAME donor ballast row: libraries and every per-gene total are exact (pds bit-identical; mse's profile half);
  * the gene's NET rank-sum change (the promotion up; ballast crossings down, mode 'b') is kept within |dAUC| <= 1e-4
    or the gene is reverted; untouched genes keep identical counts and libraries in every cell, so identical ranks;
  * a call-flip guard recomputes the scorer's BH calls and re-runs without any gene whose call flipped (0 flips).
Target per selected gene: t = max(l, K * 0.24 / sqrt(control mean count per cell)); selection: uncalled (q >= 0.05),
finite lfc, delivered lfc > 0, not the target gene. Signs are today's delivered signs (the pass never flips one).

Variant grammar: <mode>_k<K*100>[_q<quantum>r<rounds>p<pool>], mode s (ballast donors stay inside their control
intervals) or b (budgeted crossings); the optional v3 suffix spreads the recipients' sub-gate budget over more genes
(per-gene count quantum per round, rounds, recipient chosen among the pool cheapest-to-promote carriers by remaining
budget). Examples: b_k100, b_k100_q20r3p10. H1 (experiment 630): see its RESULTS.md.
"""
from __future__ import annotations

import re

import numpy as np
from scipy.stats import false_discovery_control, norm

from vcc2026 import rank_mag as ZS
from vcc2026 import rc_emit

VARIANT_RE = re.compile(r"(s|b)_k(\d{2,3})(?:_q(\d+)r(\d+)p(\d+))?")
MAX_DAUC = 1e-4
GUARD_ROUNDS = 5


def parse_spec(variant: str) -> dict:
    m = VARIANT_RE.fullmatch(str(variant).strip())
    if m is None:
        raise ValueError(f"rank_up variant must look like '<s|b>_k<K*100>', e.g. b_k100: {variant!r}")
    k = int(m.group(2)) / 100.0
    if not 0.5 <= k <= 3.0:
        raise ValueError(f"rank_up K must be in [0.5, 3]: {variant!r}")
    spec = {"cross": m.group(1) == "b", "K": k}
    if m.group(3):                                      # v3 load-spreading (quantum, rounds, recipient pool)
        spec.update(quantum=int(m.group(3)), rounds=int(m.group(4)), rec_pool=int(m.group(5)))
    return spec


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


def _min_above(cmax, s_a):
    x = int(np.floor(np.expm1(cmax) * s_a)) + 1
    while x > 0 and np.log1p((x - 1) / s_a) > cmax:
        x -= 1
    while not np.log1p(x / s_a) > cmax:
        x += 1
    return x


def _max_drop_k(x, s_b, c, budget):
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


DIAG_KEYS = ("n_sel", "n_moved", "n_reverted", "n_no_candidate", "counts", "need", "got", "n_no_pos_donor",
             "n_neg_donor_rows", "n_mean_revert", "max_dAUC", "median_dAUC", "frac_got", "capU", "guard_rounds",
             "guard_excluded", "guard_failed", "max_abs_dlog10p", "invariant_failed")


def apply(block, n_car, cc, variant, *, tie, lib_med, own_col=-1):
    """Builder entry point (after rank_mag.apply): returns (new int64 block, flat numeric diagnostics, prefix 'ur_').
    tie = scorer_tie(cc) once per context; lib_med = median control library of the context; own_col = the target
    gene's index on the kept axis (-1 if not on it).
    Fails SAFE ([7772dc]'s review of 630): if any invariant assert inside the pass fires (libraries, per-gene totals,
    negative counts, the U record, rank_mag's helpers), this target keeps its input block and the meta records
    ur_invariant_failed = 1, with the same key set as a normal call. Bit-identical to _apply when nothing fires."""
    try:
        return _apply(block, n_car, cc, variant, tie=tie, lib_med=lib_med, own_col=own_col)
    except AssertionError:
        d = {f"ur_{k}": 0.0 for k in DIAG_KEYS}
        d["ur_invariant_failed"] = 1.0
        return np.array(block, dtype=np.int64, copy=True), d


def _apply(block, n_car, cc, variant, *, tie, lib_med, own_col=-1):
    spec = parse_spec(variant) if isinstance(variant, str) else dict(variant)
    X0 = np.asarray(block, dtype=np.int64)
    p0, q0, l0 = de_stats(X0, cc, tie)
    sel = (q0 >= 0.05) & np.isfinite(l0) & (l0 > 0)
    if own_col >= 0:
        sel[own_col] = False
    mde = 0.24 / np.sqrt(np.maximum(np.asarray(cc.ref_lin, dtype=np.float64) * float(lib_med) / 1e6, 1e-12))
    t = np.maximum(l0, spec["K"] * mde)
    excl, rounds = set(), 0
    while True:
        v3 = {k: spec[k] for k in ("quantum", "rounds", "rec_pool") if k in spec}
        X1, st = up_pass2(X0, n_car, cc, l0, t, sel, cross=spec["cross"], exclude=excl, **v3)
        p1, q1, l1 = de_stats(X1, cc, tie)
        flips = set(np.flatnonzero((q0 < 0.05) != (q1 < 0.05)).tolist())
        if not flips:
            break
        if rounds >= GUARD_ROUNDS:
            # 605 review: fail SAFE, not closed. A raise here would abort a D/E/F build; this target keeps its zs block
            # (the pass is skipped for it) and the meta records it (ur_guard_failed = 1).
            for k in ("n_moved", "n_reverted", "counts", "got", "max_dAUC", "median_dAUC", "frac_got"):
                st[k] = 0
            st.update(guard_rounds=rounds, guard_excluded=len(excl), guard_failed=1, max_abs_dlog10p=0.0, invariant_failed=0)
            return X0.copy(), {f"ur_{k}": float(v) for k, v in st.items()
                               if isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool)}
        excl |= flips
        rounds += 1
    st.update(guard_rounds=rounds, guard_excluded=len(excl), guard_failed=0, invariant_failed=0,
              max_abs_dlog10p=float(np.nanmax(np.abs(np.log10(np.maximum(p1, 1e-300)) - np.log10(np.maximum(p0, 1e-300))))))
    return X1, {f"ur_{k}": float(v) for k, v in st.items() if isinstance(v, (int, float, np.integer, np.floating))
                and not isinstance(v, bool)}

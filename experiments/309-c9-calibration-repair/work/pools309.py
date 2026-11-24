"""309 -- the matched-pool builder, in two configurations.

    hard_windows=False, global_excl=False   ->  C8's pools, byte-for-byte behaviour
                                                (50 hard cells: 10 basal deciles x 5 pi_off
                                                quintiles), used for the KA-4 / KA-6 regression rows.
    hard_windows=True,  global_excl=True    ->  C9's REPAIRED pools (PREREGISTRATION sec 6.2):
                                                150 hard cells (x 3 promoter-window classes) and
                                                the global L*-candidate exclusion that makes a pool
                                                recipient-independent.

Also the cross-recipient scoring machinery the KA-4 repair needs (sec 6.1).
"""
from __future__ import annotations

import numpy as np

import lib309 as X

LS = X.LS
M_POOL = 200
MIN_CELL = 50
MIN_WIDEN = 20
COV_NAMES = ["gc", "log1p_tss_density", "map_mult", "lowcomp", "mean_offset_100",
             "log1p_n_windows", "log1p_n_sites"]

_S: dict = {}


def setup(SP, cons_rows, gene_meta, lookup, pi_off,
          hard_windows: bool = False, global_excl: bool = False) -> None:
    GIX = SP["train"].gix
    NG = SP["train"].genes.size
    COV = np.zeros((NG, len(COV_NAMES)))
    HAS_WIN = np.zeros(NG, dtype=bool)
    NWIN = np.zeros(NG)
    for g, r in gene_meta.items():
        i = GIX[g]
        HAS_WIN[i] = True
        NWIN[i] = float(r["n_windows"])
        COV[i] = [float(r["gc"]), np.log1p(float(r["tss_density_10kb"])), float(r["map_mult"]),
                  float(r["lowcomp"]), float(r["mean_offset"]) / 100.0,
                  np.log1p(float(r["n_windows"])), np.log1p(float(r["n_sites"]))]
    PI = {L: np.zeros(NG) for L in LS}
    for L in LS:
        for g, v in pi_off[L].items():
            PI[L][GIX[g]] = v
    IS_LIB = np.zeros(NG, dtype=bool)
    for g in lookup["lib_targets"]:
        IS_LIB[GIX[g]] = True

    # every L* candidate of every one of the 65 constructs (134 genes) -- the global exclusion
    IS_CAND = np.zeros(NG, dtype=bool)
    for r in cons_rows:
        if r["candidates"]:
            for g in r["candidates"].split(","):
                IS_CAND[GIX[g]] = True

    BASE, PIQ, ZCOV, WINB, REF = {}, {}, {}, {}, {}
    for L in LS:
        # quantiles / z-scoring are taken over the SAME universe C8 used, so the strata do not move
        # because of the global exclusion; only pool MEMBERSHIP does.
        ref = np.flatnonzero(HAS_WIN & (PI[L] > 0) & ~IS_LIB)
        REF[L] = ref
        ok = HAS_WIN & (PI[L] > 0) & ~IS_LIB
        if global_excl:
            ok = ok & ~IS_CAND
        BASE[L] = np.flatnonzero(ok)
        q = np.quantile(PI[L][ref], [0.2, 0.4, 0.6, 0.8])
        PIQ[L] = np.searchsorted(q, PI[L], side="right")
        mu, sd = COV[ref].mean(0), COV[ref].std(0)
        sd[sd == 0] = 1.0
        ZCOV[L] = (COV - mu) / sd
        WINB[L] = np.clip(NWIN.astype(int), 1, 3) - 1          # 3 promoter-window classes
    BASDEC = {}
    for s in X.SPLITS:
        for L in LS:
            q = np.quantile(SP[s].basal[REF[L]], np.arange(1, 10) / 10.0)
            BASDEC[(s, L)] = np.searchsorted(q, SP[s].basal, side="right")

    _S.clear()
    _S.update(dict(SP=SP, GIX=GIX, COV=COV, PI=PI, BASE=BASE, PIQ=PIQ, ZCOV=ZCOV, BASDEC=BASDEC,
                   WINB=WINB, IS_CAND=IS_CAND, hard_windows=hard_windows, global_excl=global_excl,
                   CAND_L=lookup["cand_L"], REF=REF, NWIN=NWIN,
                   SPLIT_OF={r["construct"]: (r["split"] if r["unit"] == "targeting" else "train")
                             for r in cons_rows}))


def exclusions(c: str, row: dict) -> set:
    """C8's per-construct exclusions: own intended target + every gene its spacers hit at any L."""
    bad = set()
    if row["intended_target"]:
        bad.add(row["intended_target"])
    for L in LS:
        bad |= set(_S["CAND_L"][c][L])
    return {_S["GIX"][g] for g in bad if g in _S["GIX"]}


def build_pools(c: str, row: dict) -> list:
    L = int(row["L_star"]) if row["L_star"] else 0
    if L == 0:
        return []
    s = _S["SPLIT_OF"][c]
    ix = _S["BASE"][L]
    nb = 3 if _S["hard_windows"] else 1
    cell = (_S["BASDEC"][(s, L)] * 5 + _S["PIQ"][L]) * nb + (_S["WINB"][L] if nb == 3 else 0)
    bad = exclusions(c, row)
    out = []
    for g in row["candidates"].split(","):
        gi = _S["GIX"][g]
        widened = 0
        members = ix[cell[ix] == cell[gi]]
        members = np.asarray([m for m in members if m not in bad], dtype=np.int64)
        cell_n = int(members.size)
        if members.size < MIN_WIDEN:                     # widen across adjacent basal deciles
            dec, rest = cell[gi] // (5 * nb), cell[gi] % (5 * nb)
            base_dec = _S["BASDEC"][(s, L)][ix]
            keep = ((cell[ix] % (5 * nb)) == rest) & (np.abs(base_dec - dec) <= 1)
            members = np.asarray([m for m in ix[keep] if m not in bad], dtype=np.int64)
            widened = 1
        if members.size > MIN_CELL:
            d = np.linalg.norm(_S["ZCOV"][L][members] - _S["ZCOV"][L][gi], axis=1)
            members = members[np.argsort(d, kind="stable")[:M_POOL]]
        out.append(dict(gene=g, gi=int(gi), L=L, members=members, widened=widened,
                        cell_n=cell_n))
    return out


def percentiles(e, members):
    """Mid-rank percentile of every pool member within the pool (lower = more repressed)."""
    v = e[members]
    order = np.argsort(v, kind="stable")
    sv = v[order]
    lo = np.searchsorted(sv, v, side="left")
    hi = np.searchsorted(sv, v, side="right")
    return (lo + 0.5 * (hi - lo)) / v.size


def q_of(e, p) -> float:
    """Percentile of the candidate gene itself within its pool, in expression vector e."""
    v, pv = e[p["gi"]], e[p["members"]]
    return float((np.sum(pv < v) + 0.5 * np.sum(pv == v)) / pv.size)


def z_of(e, pools) -> float:
    if not pools:
        return 0.0
    qs = [q_of(e, p) for p in pools]
    return float(np.sqrt(12 * len(qs)) * (0.5 - np.mean(qs)))


def delta_of(e, pools) -> float:
    if not pools:
        return 0.0
    return float(np.mean([e[p["gi"]] - np.median(e[p["members"]]) for p in pools]))


def balance(constructs, POOLS, SP, SPLIT_OF, PI) -> dict:
    """Nine-axis standardized mean difference, real candidates vs their matched pools."""
    COV = _S["COV"]
    bal = {}
    for name, vec in [("pi_off_Lstar", None), ("basal_log1p_cp50k", None)] + \
                     [(n, COV[:, i]) for i, n in enumerate(COV_NAMES)]:
        a, bm = [], []
        for c in constructs:
            for p in POOLS[c]:
                if vec is None and name == "pi_off_Lstar":
                    v, w = np.log10(PI[p["L"]][p["gi"]]), np.log10(PI[p["L"]][p["members"]])
                elif vec is None:
                    v, w = SP[SPLIT_OF[c]].basal[p["gi"]], SP[SPLIT_OF[c]].basal[p["members"]]
                else:
                    v, w = vec[p["gi"]], vec[p["members"]]
                a.append(float(v))
                bm.append(np.asarray(w, dtype=float))
        a = np.asarray(a)
        means = np.asarray([x.mean() for x in bm])
        sd = np.sqrt((a.var() + np.mean([x.var() for x in bm])) / 2)
        bal[name] = dict(cand_mean=float(a.mean()), pool_mean=float(means.mean()),
                         smd=float((a.mean() - means.mean()) / sd) if sd > 0 else 0.0)
    return bal


# ------------------------------------------------------------------ C8's uncentred joint null
def joint_null_uncentred(constructs, POOLS, EFF, B, rng):
    Z = np.zeros((B, len(constructs)))
    for j, c in enumerate(constructs):
        pl = POOLS[c]
        k = len(pl)
        if k == 0:
            continue
        qp = [percentiles(EFF[c], p["members"]) for p in pl]
        mem = [p["members"] for p in pl]
        idx = np.stack([rng.integers(0, qp[i].size, B) for i in range(k)], axis=1)
        if k > 1:
            for _ in range(50):
                gid = np.stack([mem[i][idx[:, i]] for i in range(k)], axis=1)
                gs = np.sort(gid, axis=1)
                bad = np.flatnonzero((gs[:, 1:] == gs[:, :-1]).any(axis=1))
                if bad.size == 0:
                    break
                for i in range(k):
                    idx[bad, i] = rng.integers(0, qp[i].size, bad.size)
        q = np.stack([qp[i][idx[:, i]] for i in range(k)], axis=1)
        Z[:, j] = np.sqrt(12 * k) * (0.5 - q.mean(axis=1))
    return Z

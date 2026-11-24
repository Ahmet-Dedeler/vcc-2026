"""307 -- the matched-pool builder of `enrich307.py`, factored out so `sens307.py` can rebuild the
identical pools (`hard_windows=False`) and a POST-HOC repaired variant (`hard_windows=True`, which
adds promoter-window count to the hard strata -- the axis KA-6 failed on).

`sens307.py` asserts that `hard_windows=False` reproduces `enrich307.py`'s S to 1e-12.
"""
from __future__ import annotations

import numpy as np

import lib307 as X

LS = (12, 11, 10, 9)
M_POOL = 200
MIN_CELL = 50
MIN_WIDEN = 20
COV_NAMES = ["gc", "log1p_tss_density", "map_mult", "lowcomp", "mean_offset_100",
             "log1p_n_windows", "log1p_n_sites"]

_S: dict = {}


def setup(SP, cons_rows, gene_meta, lookup, pi_off, hard_windows: bool = False) -> None:
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

    BASE, PIQ, ZCOV, WINB = {}, {}, {}, {}
    for L in LS:
        ix = np.flatnonzero(HAS_WIN & (PI[L] > 0) & ~IS_LIB)
        BASE[L] = ix
        q = np.quantile(PI[L][ix], [0.2, 0.4, 0.6, 0.8])
        PIQ[L] = np.searchsorted(q, PI[L], side="right")
        mu, sd = COV[ix].mean(0), COV[ix].std(0)
        sd[sd == 0] = 1.0
        ZCOV[L] = (COV - mu) / sd
        # 3 promoter-window classes: 1, 2, >=3 (only used when hard_windows)
        WINB[L] = np.clip(NWIN.astype(int), 1, 3) - 1
    BASDEC = {}
    for s in X.SPLITS:
        for L in LS:
            q = np.quantile(SP[s].basal[BASE[L]], np.arange(1, 10) / 10.0)
            BASDEC[(s, L)] = np.searchsorted(q, SP[s].basal, side="right")

    _S.clear()
    _S.update(dict(SP=SP, GIX=GIX, COV=COV, PI=PI, BASE=BASE, PIQ=PIQ, ZCOV=ZCOV,
                   BASDEC=BASDEC, WINB=WINB, hard_windows=hard_windows,
                   CAND_L=lookup["cand_L"],
                   SPLIT_OF={r["construct"]: (r["split"] if r["unit"] == "targeting" else "train")
                             for r in cons_rows}))


def _exclusions(c: str, row: dict) -> set:
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
    cell = _S["BASDEC"][(s, L)] * 5 + _S["PIQ"][L]
    if _S["hard_windows"]:
        cell = cell * 3 + _S["WINB"][L]
    bad = _exclusions(c, row)
    out = []
    for g in row["candidates"].split(","):
        gi = _S["GIX"][g]
        widened = 0
        members = ix[cell[ix] == cell[gi]]
        members = np.asarray([m for m in members if m not in bad])
        if members.size < MIN_WIDEN:
            dec = cell[gi] // (15 if _S["hard_windows"] else 5)
            rest = cell[gi] % (15 if _S["hard_windows"] else 5)
            base_dec = _S["BASDEC"][(s, L)][ix]
            keep = ((cell[ix] % (15 if _S["hard_windows"] else 5)) == rest) & \
                   (np.abs(base_dec - dec) <= 1)
            members = np.asarray([m for m in ix[keep] if m not in bad])
            widened = 1
        if members.size > MIN_CELL:
            d = np.linalg.norm(_S["ZCOV"][L][members] - _S["ZCOV"][L][gi], axis=1)
            members = members[np.argsort(d, kind="stable")[:M_POOL]]
        out.append(dict(gene=g, gi=gi, L=L, members=members, widened=widened))
    return out


def _percentiles(e, members):
    v = e[members]
    order = np.argsort(v, kind="stable")
    sv = v[order]
    lo = np.searchsorted(sv, v, side="left")
    hi = np.searchsorted(sv, v, side="right")
    return (lo + 0.5 * (hi - lo)) / v.size


def z_of(e, pools) -> float:
    if not pools:
        return 0.0
    qs = []
    for p in pools:
        v, pv = e[p["gi"]], e[p["members"]]
        qs.append((np.sum(pv < v) + 0.5 * np.sum(pv == v)) / pv.size)
    return float(np.sqrt(12 * len(qs)) * (0.5 - np.mean(qs)))


def joint_null(constructs, POOLS, EFF, B, rng):
    Z = np.zeros((B, len(constructs)))
    for j, c in enumerate(constructs):
        pl = POOLS[c]
        k = len(pl)
        if k == 0:
            continue
        qp = [_percentiles(EFF[c], p["members"]) for p in pl]
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


def balance(constructs, POOLS, SP, SPLIT_OF, PI) -> dict:
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

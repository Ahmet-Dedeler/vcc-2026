"""480: the small analyses quoted in RESULTS.md that are not in the main scripts (saved so every number has a path).

    uv run python extras480.py aplus45      # §4.2: ours + k*cons on full-power H1 restricted to h1_thin's 45 targets
    uv run python extras480.py crossstudy   # §4.5: Orion folds with the cross-study donor vs the same-study sibling
    uv run python extras480.py riderclass   # §4.6: where the k 0.2 rider's gain sits (class contributions)
    uv run python extras480.py riderrobust  # §4.6: rider vs conversion (pc) and donor set, fixed k
    uv run python extras480.py updown       # §5: H1 truth P(up) by control-expression decile vs our and donor signs
Each writes work/extras480-<name>.json.
"""
from __future__ import annotations

import json
import sys

import numpy as np

import fold480 as F
from common480 import HERE, ROOT, boot_delta, gate_from, load_donor_rows, nmae, nmae_rows
from decomp480 import contributions, tss_classes


def h1_level(gate):
    if gate is None:
        S = F.load_surface("h1")
        return S["T"], S["genes"], S["targets"], S["cnt"]
    z = np.load(HERE / f"work/gate-{gate}.npz", allow_pickle=True)
    targets, genes = z["targets"].astype(str), z["genes"].astype(str)
    T = gate_from(z["q"].astype(float), z["lfc"].astype(float), genes, targets)
    return T, genes, targets, z["ctrl_cpm"].astype(float) * float(z["lib_med"]) / 1e6


def cons_of(donors, targets, genes, pc=1.0):
    return np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, targets, genes, pc_cpm=pc)[0] for d in donors]), 0))


def aplus45():
    keep = np.isfinite(F.load_surface("h1_thin")["T"]).any(1)
    T, genes, targets, _ = h1_level(None)
    T = T.copy()
    T[~keep] = np.nan
    L0, Lc = F.ours_on(targets, genes)
    cons = cons_of(F.DON, targets, genes)
    grid = np.linspace(0, 1.5, 61)
    v = [nmae(Lc + k * cons, T) for k in grid]
    k = float(grid[int(np.argmin(v))])
    out = dict(ours=nmae(Lc, T), k=k, ours_plus=min(v), delta=min(v) - nmae(Lc, T), cons1=nmae(cons, T))
    print(out)
    return out


def crossstudy():
    out = {}
    for name in ["orion_hct116", "orion_hek293t"]:
        S = F.load_surface(name)
        T, tg, gn, cnt = S["T"], S["targets"], S["genes"], S["cnt"]
        add = F.mde_add(cnt)[None, :]
        for donors in [["k562_gwps"], [d for d in S["donors"] if d != "k562_gwps"]]:
            c = cons_of(donors, tg, gn)
            an = np.sign(c) * add + F.K_H1_MDECONS * c
            d = boot_delta(nmae_rows(c, T), nmae_rows(an, T))
            row = dict(acc=F.sign_acc(c, T)[1], k=F.best_k(c, T, 0, 3), cons1=nmae(c, T), a_analog=nmae(an, T), an_c1=d)
            out[f"{name}|{'+'.join(donors)}"] = row
            print(name, donors, {k: (np.round(v, 4) if not isinstance(v, tuple) else np.round(v, 4)) for k, v in
                                 row.items()})
    return out


def riderclass():
    st = np.load(ROOT / "experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz")
    h1 = F.load_surface("h1")
    qc_full, genes_full = np.asarray(st["q"], float), h1["genes"]
    out = {}
    for label, gate in [("400c 20k (h1-test)", None), ("250c 20k", "h1_c250")]:
        T, genes, targets, cnt = h1_level(gate)
        L0, Lc = F.ours_on(targets, genes)
        col = {g: i for i, g in enumerate(genes_full)}
        idx = np.array([col.get(g, -1) for g in genes])
        qc = np.ones((len(targets), len(genes)))
        qc[:, idx >= 0] = qc_full[:, idx[idx >= 0]]
        cons = cons_of(F.DON, targets, genes)
        s = np.sign(Lc)
        s = np.where(s == 0, np.sign(cons), s)
        s = np.where(s == 0, 1.0, s)
        G = np.isfinite(T)
        Tz = np.where(G, T, 0.0)
        D = np.abs(Tz).sum(1)
        D[D == 0] = np.nan
        E0, E1 = np.abs(Lc - Tz), np.abs(Lc + 0.2 * cons - Tz)
        q1, q2 = np.quantile(np.broadcast_to(cnt, T.shape)[G], [1 / 3, 2 / 3])
        expr = np.broadcast_to(np.where(cnt < q1, 0, np.where(cnt < q2, 1, 2)), T.shape)
        classes = {"called by us": (np.where(qc < 0.05, 0, 1), ["called", "uncalled"]),
                   "truth direction": (np.where(Tz > 0, 0, 1), ["up", "down"]),
                   "our direction": (np.where(s > 0, 0, 1), ["we say up", "we say down"]),
                   "donor agrees with our sign": (np.where(cons == 0, 2, np.where(np.sign(cons) == s, 0, 1)),
                                                  ["agree", "disagree", "none"]),
                   "expression tertile": (expr, ["low", "mid", "high"]),
                   "cis <100 kb": (np.where(tss_classes(targets, genes) <= 1, 0, 1), ["<100 kb", "rest"]),
                   "rider flips our sign": (np.where(np.sign(Lc + 0.2 * cons) != np.sign(Lc), 0, 1),
                                            ["flipped", "kept"])}
        res = {"total": float(nmae(Lc + 0.2 * cons, T) - nmae(Lc, T))}
        for var, (cls, names) in classes.items():
            d = contributions(E1, D, G, cls, len(names)) - contributions(E0, D, G, cls, len(names))
            res[var] = {n: float(x) for n, x in zip(names, d)}
        out[label] = res
        print(label, json.dumps(res))
    return out


def riderrobust():
    out = {}
    for label, gate in [("400c 20k (h1-test)", None), ("250c 20k", "h1_c250")]:
        T, genes, targets, _ = h1_level(gate)
        L0, Lc = F.ours_on(targets, genes)
        r0 = nmae_rows(Lc, T)
        for pc in [0.5, 1.0, 5.0, 20.0]:
            R = {d: load_donor_rows(d, targets, genes, pc_cpm=pc)[0] for d in F.DON}
            for dset in [F.DON, ["k562_gwps"], ["orion_hct116", "orion_hek293t"]]:
                c = np.nan_to_num(np.nanmean(np.stack([R[d] for d in dset]), 0))
                for k in [0.1, 0.2]:
                    d = boot_delta(r0, nmae_rows(Lc + k * c, T))
                    out[f"{label}|pc{pc}|{'+'.join(dset)}|k{k}"] = d
                    print(label, pc, dset, k, np.round(d, 4))
    return out


def updown():
    T, genes, targets, cnt = h1_level(None)
    L0, Lc = F.ours_on(targets, genes)
    cons = cons_of(F.DON, targets, genes)
    G = np.isfinite(T)
    C = np.broadcast_to(cnt, T.shape)
    qs = np.quantile(C[G], np.linspace(0, 1, 11))
    b = np.digitize(C, qs[1:-1])
    Tz = np.nan_to_num(T)
    out = []
    for k in range(10):
        m = G & (b == k)
        out.append(dict(decile=k, cnt_lo=float(qs[k]), cnt_hi=float(qs[k + 1]), p_up=float((Tz[m] > 0).mean()),
                        p_we_up=float((Lc[m] > 0).mean()), p_cons_up=float((cons[m] > 0).mean())))
        print(out[-1])
    return out


if __name__ == "__main__":
    name = sys.argv[1]
    res = globals()[name]()
    json.dump(res, open(HERE / f"work/extras480-{name}.json", "w"), indent=1, default=float)

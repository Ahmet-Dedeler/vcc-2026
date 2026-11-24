"""480 step 3c: does any class structure in the MDE add transfer ACROSS surfaces (leave one surface out)?

Construction on every surface: the a_analog (sign(cons) x MDE add + 0.575 cons; cons = the surface's available panel
donors). The MDE add is power-normalised, add = c * sqrt(400 / n_cells) * min(1/sqrt(cnt), 1) * [cnt >= 0.3], so c
means the same thing at 400 cells (h1) and 150 (folds); shipped = c 0.1 everywhere.
For each class variable available on every surface, per-class c_C is fitted on the other three surfaces (sum of their
nmae) and applied to the held-out one. Reported on the held-out surface: Δ vs the shipped c 0.1, and Δ vs a single
global c fitted on the same training surfaces (the class structure's own value). Target bootstrap, 4,000 draws, unit =
scored target of the held-out surface.

    uv run python loso480.py     # -> work/loso480.json
"""
from __future__ import annotations

import json

import numpy as np

import fold480 as F
from common480 import HERE, boot_delta, load_donor_rows, nmae_rows
from decomp480 import tss_classes

SURF = ["h1", "orion_hct116", "orion_hek293t", "k562_gwps"]
C_GRID = np.linspace(0, 0.8, 41)


def prep(name):
    S = F.load_surface(name)
    T, cnt, genes, targets = S["T"], S["cnt"], S["genes"], S["targets"]
    G = np.isfinite(T)
    cons = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, targets, genes, pc_cpm=1.0)[0]
                                              for d in S["donors"]]), 0))
    ncell = 400.0 if name == "h1" else float(np.median(S["ncell"]))
    add = np.sqrt(400.0 / ncell) * (np.minimum(1 / np.sqrt(np.maximum(cnt, 1e-3)), 1.0) * (cnt >= 0.3))[None, :]
    s = np.sign(cons)
    q1, q2 = np.quantile(np.broadcast_to(cnt, T.shape)[G], [1 / 3, 2 / 3])
    expr = np.broadcast_to(np.where(cnt < q1, 0, np.where(cnt < q2, 1, 2)).astype(np.int8), T.shape)
    mag = np.digitize(np.abs(cons), np.quantile(np.abs(cons)[G], [1 / 3, 2 / 3])).astype(np.int8)
    strength = (np.abs(cons) > 0.2).sum(1)
    st = np.broadcast_to(np.digitize(strength, np.quantile(strength, [1 / 3, 2 / 3]))[:, None], T.shape)
    cis = np.where(tss_classes(targets, genes) <= 1, 0, 1).astype(np.int8)
    direc = np.where(s > 0, 0, 1).astype(np.int8)
    classes = {"expression tertile": (expr, 3), "predicted direction": (direc, 2), "donor |cons| tertile": (mag, 3),
               "target strength tertile": (st.astype(np.int8), 3), "cis <100 kb": (cis, 2),
               "direction x expression": ((direc * 3 + expr).astype(np.int8), 6)}
    base = F.K_H1_MDECONS * cons
    Tz = np.where(G, T, 0.0)
    D = np.abs(Tz).sum(1)
    D[D == 0] = np.nan
    return dict(T=T, G=G, Tz=Tz, D=D, s=s, add=add, base=base, classes=classes)


def table(P, cls, n):
    """tab[i_c, target, class] contributions for c in C_GRID."""
    tab = np.zeros((len(C_GRID), P["T"].shape[0], n))
    oh = [P["G"] & (cls == k) for k in range(n)]
    for i, c in enumerate(C_GRID):
        E = np.abs(P["base"] + P["s"] * c * P["add"] - P["Tz"])
        for k in range(n):
            tab[i, :, k] = np.where(oh[k], E, 0.0).sum(1) / P["D"]
    return tab


def main():
    Ps = {n: prep(n) for n in SURF}
    j0 = int(np.argmin(np.abs(C_GRID - 0.1)))
    out = {}
    for var in Ps["h1"]["classes"]:
        n = Ps["h1"]["classes"][var][1]
        tabs = {s: table(Ps[s], *Ps[s]["classes"][var]) for s in SURF}
        glob = {s: table(Ps[s], np.zeros(Ps[s]["T"].shape, np.int8), 1) for s in SURF}
        res = {}
        for held in SURF:
            train = [s for s in SURF if s != held]
            # per-class c: minimise the sum over training surfaces of their mean-over-targets contribution
            obj = sum(np.nanmean(tabs[s], axis=1) for s in train)          # c x class
            cC = obj.argmin(0)
            objg = sum(np.nanmean(glob[s], axis=1) for s in train)
            cg = int(objg.argmin(0)[0])
            rows_cls = np.nansum(np.stack([tabs[held][cC[k], :, k] for k in range(n)], 1), 1)
            sc = np.isfinite(Ps[held]["T"]).any(1)
            rows_cls[~sc] = np.nan
            rows_glob = glob[held][cg, :, 0].copy()
            rows_ship = glob[held][j0, :, 0].copy()
            res[held] = dict(c_class=C_GRID[cC].tolist(), c_global=float(C_GRID[cg]),
                             vs_shipped=boot_delta(rows_ship, rows_cls), vs_global=boot_delta(rows_glob, rows_cls),
                             global_vs_shipped=boot_delta(rows_ship, rows_glob))
            a, b = res[held]["vs_shipped"], res[held]["vs_global"]
            print(f"{var:24s} held {held:14s} c_C {np.round(C_GRID[cC], 2).tolist()} c_glob {C_GRID[cg]:.2f} | "
                  f"vs shipped {a[0]:+.4f} [{a[1]:+.4f},{a[2]:+.4f}] | vs global {b[0]:+.4f} [{b[1]:+.4f},{b[2]:+.4f}]",
                  flush=True)
        out[var] = res
    json.dump(out, open(HERE / "work/loso480.json", "w"), indent=1)


if __name__ == "__main__":
    main()

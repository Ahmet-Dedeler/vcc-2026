"""480 §4b: is the folds' preference for donor magnitude a power effect or a line effect? A power dose-response.

For each surface (a truth gate from gates480.py) report the donor-magnitude statistics of fold480:
  k_cons      in-sample L1-optimal global scale of the donor consensus
  an-c1       a_analog - cons@1 (the H1-calibrated small-magnitude construction vs pack-style scale 1; the board's one
              disclosed pair says this is negative on A/B/C)
  ours-c1     (H1 surfaces) our shipped construction - cons@1
  +cons       (H1 surfaces) gain of adding k*cons to our shipped construction, k in-sample
Surfaces: h1-test at 400/250/150/100 cells per target and at 1.0 / 0.56-0.6 depth; the K562-GWPS fold at 150 cells
(267 targets) and at 250 cells (the 55 targets with >= 250 cells), plus the 150-cell fold restricted to those 55; the
Orion folds. Target bootstrap (4,000 draws, unit = scored target).

    uv run python power480.py     # -> work/power480.json
"""
from __future__ import annotations

import json

import numpy as np

import fold480 as F
from common480 import HERE, boot_delta, gate_from, load_donor_rows, nmae, nmae_rows

SURF = [  # (surface name, gate file, line, target restriction from another gate or None)
    ("h1 400c 20k", None, "h1", None),
    ("h1 250c 20k", "h1_c250", "h1", None),
    ("h1 150c 20k", "h1_c150", "h1", None),
    ("h1 100c 20k", "h1_c100", "h1", None),
    ("h1 400c 11k", "h1_c400_d56", "h1", None),
    ("h1 250c 11k", "h1_c250_d56", "h1", None),
    ("h1 150c 11k", "h1_c150_d56", "h1", None),
    ("h1 150c 12k (h1_thin)", "h1_thin", "h1", None),
    ("K562 150c 11k (267 targets)", "k562_gwps", "k562_gwps", None),
    ("K562 150c 11k (55 targets)", "k562_gwps", "k562_gwps", "k562_c250"),
    ("K562 250c 11k (55 targets)", "k562_c250", "k562_gwps", None),
    ("HCT116 150c 12k", "orion_hct116", "orion_hct116", None),
    ("HEK293T 150c 14k", "orion_hek293t", "orion_hek293t", None),
]


def load(gate, line, restrict):
    if gate is None:
        S = F.load_surface("h1")
        T, cnt, genes, targets = S["T"], S["cnt"], S["genes"], S["targets"]
    else:
        z = np.load(HERE / f"work/gate-{gate}.npz", allow_pickle=True)
        targets, genes = z["targets"].astype(str), z["genes"].astype(str)
        T = gate_from(z["q"].astype(float), z["lfc"].astype(float), genes, targets)
        cnt = z["ctrl_cpm"].astype(float) * float(z["lib_med"]) / 1e6
    if restrict is not None:
        keep = set(np.load(HERE / f"work/gate-{restrict}.npz", allow_pickle=True)["targets"].astype(str))
        T[[t not in keep for t in targets]] = np.nan
    donors = F.DON if line == "h1" else [d for d in F.DON if d != line]
    return T, cnt, genes, targets, donors


def main():
    out = {}
    print(f"{'surface':30s} {'n':>3s} {'gate':>5s} {'acc':>5s} {'k':>5s} | {'an-c1':>24s} | {'ours-c1':>24s} | "
          f"{'+cons (k)':>28s}")
    for label, gate, line, restrict in SURF:
        T, cnt, genes, targets, donors = load(gate, line, restrict)
        G = np.isfinite(T)
        scored = G.any(1)
        cons = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, targets, genes, pc_cpm=1.0)[0]
                                                  for d in donors]), 0))
        s = np.sign(cons)
        add = F.mde_add(cnt)[None, :]
        k = F.best_k(cons, T, 0.0, 3.0)
        an = s * add + F.K_H1_MDECONS * cons
        r_c1 = nmae_rows(cons, T)
        d_an = boot_delta(r_c1, nmae_rows(an, T))
        acc = F.sign_acc(cons, T)[1]
        row = dict(n=int(scored.sum()), gate_med=float(np.median(G.sum(1)[scored])), acc=acc, k=k, an_c1=d_an)
        txt_o = txt_p = ""
        if line == "h1":
            L0, Lc = F.ours_on(targets, genes)
            d_o = boot_delta(r_c1, nmae_rows(Lc, T))
            grid = np.linspace(0, 1.5, 61)
            v = [nmae(Lc + kk * cons, T) for kk in grid]
            kp = float(grid[int(np.argmin(v))])
            d_p = boot_delta(nmae_rows(Lc, T), nmae_rows(Lc + kp * cons, T))
            row.update(ours_c1=d_o, plus_cons=d_p, k_plus=kp, ours=float(nmae(Lc, T)))
            txt_o = f"{d_o[0]:+.4f} [{d_o[1]:+.4f},{d_o[2]:+.4f}]"
            txt_p = f"{d_p[0]:+.4f} [{d_p[1]:+.4f},{d_p[2]:+.4f}] ({kp:.2f})"
        out[label] = row
        print(f"{label:30s} {row['n']:3d} {row['gate_med']:5.0f} {acc:5.3f} {k:5.2f} | "
              f"{d_an[0]:+.4f} [{d_an[1]:+.4f},{d_an[2]:+.4f}] | {txt_o:24s} | {txt_p:28s}", flush=True)
    json.dump(out, open(HERE / "work/power480.json", "w"), indent=1)


if __name__ == "__main__":
    main()

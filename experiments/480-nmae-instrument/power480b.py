"""480 §4c: H1 at graded power as the instrument for magnitude levers. Which power matches A/B/C?

For each H1 power level (cells per target x depth), on the SAME predictions:
  gate_all    median gate size over all 85 targets (the quantity the board-derived n_real estimates refer to:
              audit 181 from board reach: median 90-110; audit 192 from board fid coverage: median 31-45)
  mde         our real MDE add (435 coarse - 426 base); the board read the same lever at -0.0168 raw (1.76x full H1)
  b_k562      pack-style (b): K562-GWPS same-target response at scale 1 (pc1 and naive20) minus ours (coarse); the
              board's one disclosed pair read +0.044 (SQ W 09-07 vs our 09-07 construction)
  c1          donor consensus at scale 1 minus ours
  plus_cons   gain of ours + k*cons over ours (k in-sample; the most favourable reading of a donor-magnitude rider)
Target bootstrap, 4,000 draws, unit = scored target.

    uv run python power480b.py    # -> work/power480b.json
"""
from __future__ import annotations

import json

import numpy as np

import fold480 as F
from common480 import HERE, ROOT, boot_delta, gate_from, load_donor_rows, nmae, nmae_rows

LEVELS = [("400c 20k (h1-test)", None), ("250c 20k", "h1_c250"), ("400c 11k", "h1_c400_d56"),
          ("150c 20k", "h1_c150"), ("250c 11k", "h1_c250_d56"), ("100c 20k", "h1_c100"),
          ("150c 11k", "h1_c150_d56"), ("150c 12k", "h1_thin")]


def main():
    out = {}
    print(f"{'level':20s} {'gate_all':>8s} {'scored':>6s} {'ours':>6s} | {'mde':>24s} | {'b_k562 pc1':>24s} | "
          f"{'b_k562 naive20':>24s} | {'c1':>24s} | {'plus_cons (k)':>30s}")
    for label, gate in LEVELS:
        if gate is None:
            z = np.load(HERE / "work/gate-h1_full_ka.npz", allow_pickle=True)
        else:
            z = np.load(HERE / f"work/gate-{gate}.npz", allow_pickle=True)
        targets, genes = z["targets"].astype(str), z["genes"].astype(str)
        q, lfc = z["q"].astype(float), z["lfc"].astype(float)
        gate_all = [int(((q[k] < 0.05) & np.isfinite(lfc[k]) & (lfc[k] != 0) & (genes != t)).sum())
                    for k, t in enumerate(targets)]
        if gate is None:                                    # the validated 323 reference for the full surface
            S = F.load_surface("h1")
            T, genes, targets = S["T"], S["genes"], S["targets"]
        else:
            T = gate_from(q, lfc, genes, targets)
        L0, Lc = F.ours_on(targets, genes)
        r_o = nmae_rows(Lc, T)
        mde = boot_delta(nmae_rows(L0, T), r_o)
        kp = load_donor_rows("k562_gwps", targets, genes, pc_cpm=1.0)[0]
        kn = load_donor_rows("k562_gwps", targets, genes, pc_cpm=20.0)[0]
        cons = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, targets, genes, pc_cpm=1.0)[0]
                                                  for d in F.DON]), 0))
        b1 = boot_delta(r_o, nmae_rows(np.nan_to_num(kp), T))
        b2 = boot_delta(r_o, nmae_rows(np.nan_to_num(kn), T))
        c1 = boot_delta(r_o, nmae_rows(cons, T))
        grid = np.linspace(0, 1.5, 61)
        v = [nmae(Lc + kk * cons, T) for kk in grid]
        k = float(grid[int(np.argmin(v))])
        pc = boot_delta(r_o, nmae_rows(Lc + k * cons, T))
        row = dict(gate_all_median=float(np.median(gate_all)), scored=int(np.isfinite(T).any(1).sum()),
                   ours=float(np.nanmean(r_o)), mde=mde, b_k562_pc1=b1, b_k562_naive20=b2, c1=c1, plus_cons=pc,
                   k_plus=k)
        out[label] = row
        f = lambda d: f"{d[0]:+.4f} [{d[1]:+.4f},{d[2]:+.4f}]"
        print(f"{label:20s} {row['gate_all_median']:8.0f} {row['scored']:6d} {row['ours']:.4f} | {f(mde)} | {f(b1)} | "
              f"{f(b2)} | {f(c1)} | {f(pc)} ({k:.2f})", flush=True)
    json.dump(out, open(HERE / "work/power480b.json", "w"), indent=1)


if __name__ == "__main__":
    main()

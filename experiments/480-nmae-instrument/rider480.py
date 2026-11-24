"""480 §4e: the donor-consensus rider on our DE lfc, priced at the H1 power levels compatible with A/B/C.

Rider: lhat = ours (shipped coarse) + k * X, with X = cons (mixed: can flip our sign where |k cons| > |ours|) or
X = s * |cons| (magnitude only: our sign kept everywhere). cons = mean of the three panel donors' same-target log2
ratios (pc 1 CPM). For each power level: a fixed-k scan (k 0.1 ... 0.5, no fitting) and a target-CV k (fitted on 4/5
of targets, applied to the held-out fifth; 5-fold x 50 partitions; per-target held-out Δ averaged; target bootstrap,
4,000 draws). Power levels are those whose nmae replicate anchor sits in A/B/C's band (replicate480.py).

    uv run python rider480.py    # -> work/rider480.json
"""
from __future__ import annotations

import json

import numpy as np

import fold480 as F
from common480 import HERE, ROOT, gate_from, load_donor_rows, nmae_rows

LEVELS = [("400c 20k (h1-test)", None), ("250c 20k", "h1_c250"), ("400c 11k", "h1_c400_d56"),
          ("150c 20k", "h1_c150")]
K_FIX = [0.1, 0.2, 0.3, 0.4, 0.5]
K_GRID = np.linspace(0, 1.0, 41)
R, NB = 50, 4000


def boot(x, rng):
    x = x[np.isfinite(x)]
    bs = x[rng.integers(0, x.size, size=(NB, x.size))].mean(1)
    return float(x.mean()), float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))


def main():
    rng = np.random.default_rng(0)
    ess = set()
    for d in ["k562_essential", "rpe1", "nadig_jurkat"]:
        ess |= {l.strip() for l in open(ROOT / f"data/donors/{d}/perts.txt")}
    out = {}
    for label, gate in LEVELS:
        if gate is None:
            S = F.load_surface("h1")
            T, genes, targets = S["T"], S["genes"], S["targets"]
        else:
            z = np.load(HERE / f"work/gate-{gate}.npz", allow_pickle=True)
            targets, genes = z["targets"].astype(str), z["genes"].astype(str)
            T = gate_from(z["q"].astype(float), z["lfc"].astype(float), genes, targets)
        scored = np.isfinite(T).any(1)
        pl = np.array([t not in ess for t in targets])
        L0, Lc = F.ours_on(targets, genes)
        cons = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, targets, genes, pc_cpm=1.0)[0]
                                                  for d in F.DON]), 0))
        s = np.sign(Lc)
        s = np.where(s == 0, np.sign(cons), s)
        base_rows = nmae_rows(Lc, T)
        res = {}
        for kind, X in [("mixed", cons), ("magnitude-only", s * np.abs(cons))]:
            tab = np.stack([nmae_rows(Lc + k * X, T) - base_rows for k in K_GRID])      # k x target
            fixed = {f"{k:.1f}": boot(tab[int(np.argmin(np.abs(K_GRID - k)))].copy(), rng) for k in K_FIX}
            idx = np.flatnonzero(scored)
            acc = np.zeros(T.shape[0])
            cnt = np.zeros(T.shape[0])
            ks = []
            for _ in range(R):
                for f in np.array_split(rng.permutation(idx), 5):
                    tr = np.setdiff1d(idx, f)
                    j = int(np.argmin(np.nanmean(tab[:, tr], 1)))
                    ks.append(K_GRID[j])
                    acc[f] += tab[j, f]
                    cnt[f] += 1
            d = np.full(T.shape[0], np.nan)
            d[idx] = acc[idx] / cnt[idx]
            res[kind] = dict(fixed=fixed, cv=boot(d.copy(), rng), cv_panel_like=boot(d[pl].copy(), rng),
                             k_cv_mean=float(np.mean(ks)))
            fx = " ".join(f"k{k}:{v[0]:+.4f}" for k, v in fixed.items())
            c, cp = res[kind]["cv"], res[kind]["cv_panel_like"]
            print(f"{label:20s} {kind:15s} fixed {fx} | CV {c[0]:+.4f} [{c[1]:+.4f},{c[2]:+.4f}] (k {np.mean(ks):.2f}) "
                  f"| panel-like {cp[0]:+.4f} [{cp[1]:+.4f},{cp[2]:+.4f}]", flush=True)
        out[label] = res
    json.dump(out, open(HERE / "work/rider480.json", "w"), indent=1)


if __name__ == "__main__":
    main()

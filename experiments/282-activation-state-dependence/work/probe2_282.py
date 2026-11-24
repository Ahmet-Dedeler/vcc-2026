"""282 POST-HOC probe 2: the DECONFOUNDED across-state statistic at the donor-PAIR scale.

P-A showed the preregistered S1 (across state, donors AND construct held) is a (donor x guide) property:
swapping to the target's other guide in the same person destroys it.  So rebuild S1 / S2 / S12 with the
GUIDE SWAPPED on the two sides of every comparison -- the construct-specific component cannot then be
shared -- on the 74 targets that carry the same two guides in all 4 donors and all 3 states.

    PYTHONDONTWRITEBYTECODE=1 uv run python probe2_282.py
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, ".")
from lib282 import (CONDS, DORDER, PAIRINGS, PAIRS, SEED, STATE_PAIRS, W280, Bank, boot_median,  # noqa: E402
                    boot_median_diff, build_index, cosines, sd_index, thin)

P = lambda *a: print(*a, flush=True)  # noqa: E731
t0 = time.time()
bank = Bank(W280)
SC, CC = bank.scored, bank.cols_cd4
T, IDX = build_index(bank, "common")
TG = [t for t in T if all(len(IDX[(c, d, t)]) == 2 for c in CONDS for d in DORDER)]
P(f"{len(TG)} targets, same two guides in all 4 donors x 3 states")
rng = np.random.default_rng(SEED)
# per-target common budget over the 18 (donor pair, state) x 2 guide cells
bud = {t: min(float(IDX[(c, d, t)].total.values[k]) + float(IDX[(c, d2, t)].total.values[k])
              for c in CONDS for d, d2 in PAIRS for k in (0, 1)) for t in TG}
P(f"per-target budget for the matched arm: median {np.median(list(bud.values())):.3g} UMI")


def pair_guide_rows(matched):
    """(donor pair, guide rank, state) -> residual matrix over TG (both donors of the pair, ONE guide)."""
    out = {}
    for c in CONDS:
        for p in PAIRS:
            b = bank.ctrl_profile(c, p)[CC]
            for k in (0, 1):
                M = np.zeros((len(TG), len(CC)))
                for j, t in enumerate(TG):
                    S = bank.counts([int(IDX[(c, d, t)].i.values[k]) for d in p])
                    if matched:
                        S = thin(S, bud[t], rng)
                    M[j] = bank.profile(S)[CC] - b
                M -= M.mean(0, keepdims=True)
                M[:, ~SC] = 0.0
                out[("".join(p), k, c)] = M[:, SC]
    return out


rows = []
res = {}
for matched in (False, True):
    Rg = pair_guide_rows(matched)
    cells = {k: [] for k in ("S1x", "S2x", "S12x", "S0x", "S1", "S2", "S12")}
    for p in PAIRS:
        pk = "".join(p)
        for si, sj in STATE_PAIRS:
            cells["S1x"].append(cosines(Rg[(pk, 0, si)], Rg[(pk, 1, sj)]))   # same donors, GUIDE SWAPPED, diff state
            cells["S1x"].append(cosines(Rg[(pk, 1, si)], Rg[(pk, 0, sj)]))
            for k in (0, 1):
                cells["S1"].append(cosines(Rg[(pk, k, si)], Rg[(pk, k, sj)]))  # same donors, same guide (confounded)
        for c in CONDS:
            cells["S0x"].append(cosines(Rg[(pk, 0, c)], Rg[(pk, 1, c)]))     # same donors, GUIDE SWAPPED, same state
    for a, b_ in PAIRINGS:
        ak, bk = "".join(a), "".join(b_)
        for c in CONDS:
            cells["S2x"].append(cosines(Rg[(ak, 0, c)], Rg[(bk, 1, c)]))     # diff donors, GUIDE SWAPPED, same state
            cells["S2x"].append(cosines(Rg[(ak, 1, c)], Rg[(bk, 0, c)]))
            for k in (0, 1):
                cells["S2"].append(cosines(Rg[(ak, k, c)], Rg[(bk, k, c)]))  # diff donors, same guide
        for si, sj in STATE_PAIRS:
            cells["S12x"].append(cosines(Rg[(ak, 0, si)], Rg[(bk, 1, sj)]))
            cells["S12x"].append(cosines(Rg[(ak, 1, si)], Rg[(bk, 0, sj)]))
            for k in (0, 1):
                cells["S12"].append(cosines(Rg[(ak, k, si)], Rg[(bk, k, sj)]))
    agg = {k: np.nanmean(np.vstack(v), axis=0) for k, v in cells.items()}
    lab = {"S1": "S1  same donors, SAME guide,  diff state   [CONFOUNDED]",
           "S1x": "S1x same donors, guide SWAPPED, diff state  [deconfounded]",
           "S0x": "S0x same donors, guide SWAPPED, SAME state  [within-state comparator]",
           "S2x": "S2x diff donors, guide SWAPPED, same state  [deconfounded ceiling]",
           "S2": "S2  diff donors, SAME guide,  same state   [280's statistic]",
           "S12x": "S12x diff donors, guide SWAPPED, diff state [floor]",
           "S12": "S12 diff donors, SAME guide,  diff state"}
    P(f"\n--- {'depth-matched' if matched else 'full depth'} (~{107 if not matched else 0} cells/side) ---")
    for k in ("S1", "S1x", "S0x", "S2", "S2x", "S12", "S12x"):
        m, lo, hi, n = boot_median(agg[k])
        P(f"  {lab[k]:<52} n={n:3d} median {m:+.4f} [{lo:+.4f},{hi:+.4f}]  <=0 {np.nanmean(agg[k]<=0):.3f}")
        rows.append(dict(matched=matched, cls=k, label=lab[k], n=n, median=m, lo=lo, hi=hi,
                         frac_le0=float(np.nanmean(agg[k] <= 0))))
    for nm, x, y in (("S2x - S1x  (state costs, construct broken)", agg["S2x"], agg["S1x"]),
                     ("S0x - S1x  (state costs, same donors)", agg["S0x"], agg["S1x"]),
                     ("S2  - S1   (the preregistered contrast)", agg["S2"], agg["S1"]),
                     ("S2x - S12x (value of sharing the state)", agg["S2x"], agg["S12x"]),
                     ("S1x - S12x (value of sharing the donors)", agg["S1x"], agg["S12x"]),
                     ("S1  - S1x  (the construct component)", agg["S1"], agg["S1x"])):
        d, lo, hi, n = boot_median_diff(x, y)
        P(f"    {nm:<50} {d:+.4f} [{lo:+.4f},{hi:+.4f}]")
        rows.append(dict(matched=matched, cls=nm, label="contrast", n=n, median=d, lo=lo, hi=hi, frac_le0=np.nan))
    for nm, num, den in (("SD  = 1 - S1 /S2   [preregistered, confounded]", agg["S1"], agg["S2"]),
                         ("SDx = 1 - S1x/S2x  [deconfounded]", agg["S1x"], agg["S2x"]),
                         ("SDx0= 1 - S1x/S0x  [vs the within-state, same-donor comparator]",
                          agg["S1x"], agg["S0x"])):
        sd = sd_index(num, den)
        m, lo, hi, n = boot_median(sd)
        pooled = 1 - np.nanmedian(num) / np.nanmedian(den)
        P(f"    {nm:<50} n={n:3d} median {m:+.4f} [{lo:+.4f},{hi:+.4f}]  pooled {pooled:+.4f}")
        rows.append(dict(matched=matched, cls=nm, label="SD", n=n, median=m, lo=lo, hi=hi, frac_le0=float(pooled)))
    res[matched] = agg
    # depth diagnostics on this scale
    cellcount = np.array([sum(float(IDX[(c, d, t)].n_cells.sum()) for c in CONDS for d in DORDER) for t in TG])
    for k in ("S1", "S1x", "S2x"):
        v = agg[k]
        m = np.isfinite(v)
        P(f"    Spearman({k:<4}, cells per target) {spearmanr(cellcount[m], v[m]).statistic:+.3f}")
    del Rg
pd.DataFrame(rows).to_csv("pair_guide282.csv", index=False)
np.savez("pair_guide282.npz", targets=np.asarray(TG),
         **{f"{k}_{'m' if mm else 'f'}": v for mm, agg in res.items() for k, v in agg.items()})
json.dump(dict(n_targets=len(TG), elapsed_s=time.time() - t0), open("probe2_282.json", "w"), indent=1)
P(f"[{time.time()-t0:.0f}s]")

"""282 POST-HOC probes (labelled as such in RESULTS.md), added because the preregistered primary read
S1 (across state, donors held) FAR ABOVE S2 (across donors, state held) — the opposite of what the
preregistration anticipated.  Two things then decide the interpretation:

  P-A  the crossed 2x2x2 (donor x guide x state) at the (donor, guide, state) level: is the across-state
       agreement a property of the TARGET's response in that person, or of the (donor, guide) construct?
  P-B  a GLOBAL equal-depth arm: every (pair, state, target) thinned to one common UMI budget, so a
       target-specific depth/detection direction cannot be shared across states.
  P-C  depth diagnostics and the rank-1 share of the shared structure.

    PYTHONDONTWRITEBYTECODE=1 uv run python probe282.py
"""
from __future__ import annotations

import json
import sys
import time
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, ".")
from lib282 import (CONDS, DORDER, PAIRINGS, PAIRS, SEED, STATE_PAIRS, W280, Bank, boot_median,  # noqa: E402
                    boot_median_diff, build_index, build_rows, cosines, stat_S1, stat_S12, stat_S2, thin)

P = lambda *a: print(*a, flush=True)  # noqa: E731
t0 = time.time()
OUT: dict = {}
bank = Bank(W280)
SC, CC = bank.scored, bank.cols_cd4
T, IDX = build_index(bank, "common")

# ==================================================================== P-A the crossed 2x2x2
P("=== P-A  crossed donor x guide x state (identical guide ids in all 4 donors) ===")
TG = [t for t in T if all(len(IDX[(c, d, t)]) == 2 for c in CONDS for d in DORDER)]
gid = {t: list(IDX[("Rest", "D1", t)].guide.values) for t in TG}
assert all(len({tuple(sorted(IDX[(c, d, t)].guide.values)) for c in CONDS for d in DORDER}) == 1 for t in TG)
P(f"  {len(TG)} targets with the SAME two guides in all 4 donors x 3 states")

rng = np.random.default_rng(SEED)
budget = {t: min(float(IDX[(c, d, t)].total.values[k]) for c in CONDS for d in DORDER for k in (0, 1))
          for t in TG}
P(f"  per-target common budget = min over the 24 (donor, guide, state) cells: "
  f"median {np.median(list(budget.values())):.3g} UMI")


def cell_rows(matched: bool):
    """(donor, guide-rank, state) -> residual matrix over TG, centred on its own across-target mean."""
    out = {}
    for c in CONDS:
        for d in DORDER:
            b = bank.ctrl_profile(c, [d])[CC]
            for k in (0, 1):
                M = np.zeros((len(TG), len(CC)))
                for j, t in enumerate(TG):
                    S = bank.counts([int(IDX[(c, d, t)].i.values[k])])
                    if matched:
                        S = thin(S, budget[t], rng)
                    M[j] = bank.profile(S)[CC] - b
                M -= M.mean(0, keepdims=True)
                M[:, ~SC] = 0.0
                out[(d, k, c)] = M[:, SC]
    return out


rows_pa = []
for matched in (False, True):
    Rg = cell_rows(matched)
    cls: dict[str, list] = {k: [] for k in ("A_s", "A_g", "A_gs", "B_g", "B_s", "B_gg", "B_all")}
    for d in DORDER:
        for si, sj in STATE_PAIRS:
            for k in (0, 1):
                cls["A_s"].append(cosines(Rg[(d, k, si)], Rg[(d, k, sj)]))          # same D, same G, diff S
            cls["A_gs"].append(cosines(Rg[(d, 0, si)], Rg[(d, 1, sj)]))             # same D, diff G, diff S
            cls["A_gs"].append(cosines(Rg[(d, 1, si)], Rg[(d, 0, sj)]))
        for c in CONDS:
            cls["A_g"].append(cosines(Rg[(d, 0, c)], Rg[(d, 1, c)]))                # same D, diff G, same S
    for d1, d2 in combinations(DORDER, 2):
        for c in CONDS:
            for k in (0, 1):
                cls["B_g"].append(cosines(Rg[(d1, k, c)], Rg[(d2, k, c)]))          # diff D, same G, same S
            cls["B_gg"].append(cosines(Rg[(d1, 0, c)], Rg[(d2, 1, c)]))             # diff D, diff G, same S
            cls["B_gg"].append(cosines(Rg[(d1, 1, c)], Rg[(d2, 0, c)]))
        for si, sj in STATE_PAIRS:
            for k in (0, 1):
                cls["B_s"].append(cosines(Rg[(d1, k, si)], Rg[(d2, k, sj)]))        # diff D, same G, diff S
                cls["B_s"].append(cosines(Rg[(d1, k, sj)], Rg[(d2, k, si)]))
            cls["B_all"].append(cosines(Rg[(d1, 0, si)], Rg[(d2, 1, sj)]))          # all three differ
            cls["B_all"].append(cosines(Rg[(d1, 1, si)], Rg[(d2, 0, sj)]))
    agg = {k: np.nanmean(np.vstack(v), axis=0) for k, v in cls.items()}
    lab = {"A_s": "same donor, SAME guide, diff STATE", "A_g": "same donor, DIFF guide, same state",
           "A_gs": "same donor, DIFF guide, diff STATE", "B_g": "DIFF donor, SAME guide, same state",
           "B_s": "DIFF donor, SAME guide, diff STATE", "B_gg": "DIFF donor, DIFF guide, same state",
           "B_all": "DIFF donor, DIFF guide, diff STATE"}
    P(f"  --- {'depth-matched to the common per-target budget' if matched else 'full depth'} ---")
    for k in ("A_s", "A_g", "A_gs", "B_g", "B_gg", "B_s", "B_all"):
        m, lo, hi, n = boot_median(agg[k])
        P(f"    {k:<5} {lab[k]:<38} n={n:3d} median {m:+.4f} [{lo:+.4f},{hi:+.4f}]  <=0 {np.nanmean(agg[k]<=0):.3f}")
        rows_pa.append(dict(matched=matched, cls=k, label=lab[k], n=n, median=m, lo=lo, hi=hi,
                            frac_le0=float(np.nanmean(agg[k] <= 0))))
    for nm, x, y in (("GUIDE cost, state already differs   (A_s - A_gs)", agg["A_s"], agg["A_gs"]),
                     ("STATE cost, guide already differs   (A_g - A_gs)", agg["A_g"], agg["A_gs"]),
                     ("DONOR cost, guide held, state differs (A_s - B_s)", agg["A_s"], agg["B_s"]),
                     ("DONOR cost, guide held, same state  (A_g?? n/a) (B_g - B_gg)", agg["B_g"], agg["B_gg"]),
                     ("STATE cost, donor+guide held vs guide swap (A_s - A_g)", agg["A_s"], agg["A_g"])):
        d, lo, hi, n = boot_median_diff(x, y)
        P(f"      {nm:<58} {d:+.4f} [{lo:+.4f},{hi:+.4f}]")
        rows_pa.append(dict(matched=matched, cls=nm, label="contrast", n=n, median=d, lo=lo, hi=hi,
                            frac_le0=np.nan))
    if matched:
        np.savez("pa282_matched.npz", targets=np.asarray(TG), **agg)
    del Rg
pd.DataFrame(rows_pa).to_csv("pa282.csv", index=False)

# ==================================================================== P-B global equal depth
P(f"\n=== P-B  GLOBAL equal-depth arm  [{time.time()-t0:.0f}s] ===")
mins = {t: min(sum(float(IDX[(c, d, t)].total.sum()) for d in p) for c in CONDS for p in PAIRS) for t in T}
B = float(np.quantile(list(mins.values()), 0.10))
keep = [t for t in T if mins[t] >= B]
P(f"  global budget B = {B:.4g} UMI (10th pct of the per-target minima over the 18 cells); "
  f"{len(keep)}/{len(T)} targets kept, every cell thinned to exactly B")
Rb = build_rows(bank, IDX, keep, budget={t: B for t in keep}, seed=SEED)
s1, _ = stat_S1(Rb, keep)
s2, _ = stat_S2(Rb, keep)
s12, _ = stat_S12(Rb, keep)
rows_pb = []
for nm, x in (("S1  across STATE, same donors", s1), ("S2  across DONORS, same state", s2),
              ("S12 across state AND donors", s12)):
    m, lo, hi, n = boot_median(x)
    P(f"    {nm:<34} n={n:3d} median {m:+.4f} [{lo:+.4f},{hi:+.4f}]  <=0 {np.nanmean(x<=0):.3f}")
    rows_pb.append(dict(stat=nm, n=n, median=m, lo=lo, hi=hi, frac_le0=float(np.nanmean(x <= 0))))
d, lo, hi, n = boot_median_diff(s2, s1)
P(f"    S2 - S1                            {d:+.4f} [{lo:+.4f},{hi:+.4f}]")
rows_pb.append(dict(stat="S2 - S1", n=n, median=d, lo=lo, hi=hi, frac_le0=np.nan))
pd.DataFrame(rows_pb).to_csv("pb282.csv", index=False)
OUT["PB"] = dict(budget=B, n_targets=len(keep), S1=float(np.nanmedian(s1)), S2=float(np.nanmedian(s2)),
                 S12=float(np.nanmedian(s12)))
del Rb

# ==================================================================== P-C depth diagnostics
P(f"\n=== P-C  depth diagnostics on the primary arm  [{time.time()-t0:.0f}s] ===")
tab = pd.read_csv("state282.csv")
cells = np.array([sum(float(IDX[(c, d, t)].n_cells.sum()) for c in CONDS for d in DORDER) for t in tab.target])
tot = np.array([sum(float(IDX[(c, d, t)].total.sum()) for c in CONDS for d in DORDER) for t in tab.target])
for nm in ("S1", "S2", "S12", "frac_shared", "frac_state", "n_de_genes_Stim48hr"):
    v = tab[nm].values.astype(float)
    m = np.isfinite(v)
    r1 = spearmanr(cells[m], v[m]).statistic
    r2 = spearmanr(tot[m], v[m]).statistic
    P(f"    Spearman({nm:<20}, cells) {r1:+.3f}   (, total UMI) {r2:+.3f}   n={m.sum()}")
    OUT.setdefault("PC", {})[nm] = dict(rho_cells=float(r1), rho_umi=float(r2), n=int(m.sum()))
v = tab.n_de_genes_Stim48hr.values.astype(float)
for nm in ("S1", "S2"):
    x = tab[nm].values.astype(float)
    m = np.isfinite(v) & np.isfinite(x)
    P(f"    Spearman({nm}, n_DE_genes) {spearmanr(v[m], x[m]).statistic:+.3f}  n={m.sum()}")
OUT["elapsed_s"] = time.time() - t0
json.dump(OUT, open("probe282.json", "w"), indent=1, default=float)
P(f"[{time.time()-t0:.0f}s]")

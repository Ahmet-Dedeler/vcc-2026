"""280 Q1 -- across-donor reliability of a deep Flex response in primary human CD4 T cells.

Sections (all preregistered in ../PREREGISTRATION.md, sha256 7273ef3c…):
  3.2  donor halves {D1,D2}|{D3,D4}, {D1,D3}|{D2,D4}, {D1,D4}|{D2,D3}      -> headline reliability, P1, P4
  3.1  the 2x2 D / A / B / C at a common counting budget                    -> Δ_donor, P2
  P5   label-shuffled null (KA-4)
  P3   Spearman with 276's conservation
  KA-2 the pooled four-donor Rest row vs the on-disk cd4_compact bank
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, ".")
from lib280 import (CONDS, DORDER, PAIRINGS, SEED, Bank, boot_median, boot_median_diff, cosines,  # noqa: E402
                    residualise, thin)

ROOT = "~/Code/virtual-cell-challenge-2026"
KEEP = True          # primary: the authors' keep_for_DE guides (sensitivity: all guides, --all)
if "--all" in sys.argv:
    KEEP = False
TAG = "keepDE" if KEEP else "allguides"
P = lambda *a: print(*a, flush=True)  # noqa: E731
t0 = time.time()
bank = Bank(".")
SC = bank.scored
CC, CO = bank.cols_cd4, bank.cols_off
OUT: dict = dict(tag=TAG, seed=SEED)
P(f"bank: X {bank.X.shape}, mapped cols {len(CC)}, scored {SC.sum()}  [{time.time()-t0:.0f}s]")


def side_delta(cond, donors, rows_by_target, targets, thin_to=None, rng=None):
    """Rows for each target -> mapped-column delta matrix (n_targets, n_mapped)."""
    b = bank.ctrl_profile(cond, donors)[CC]
    out = np.zeros((len(targets), len(CC)))
    for k, t in enumerate(targets):
        S = bank.counts(rows_by_target[t])
        if thin_to is not None:
            S = thin(S, thin_to[t] if isinstance(thin_to, dict) else thin_to, rng)
        out[k] = bank.profile(S)[CC] - b
    return out


def resid(deltas):
    r = residualise(deltas)
    r[:, ~SC] = 0.0
    return r[:, SC]


# ------------------------------------------------------------------ target/guide index
o = bank.obs
tg = o[(o.gtype == "targeting")]
if KEEP:
    tg = tg[tg.keep_for_DE]
IDX: dict = {}                              # (cond, D, target) -> DataFrame of guide rows, guide-sorted
for (c, d, t), g in tg.groupby(["cond", "D", "pert"]):
    IDX[(c, d, t)] = g.sort_values("guide")
ALLT = sorted(set(tg.pert))

# ==================================================================== 3.2 donor halves
P("\n=== 3.2 donor halves (all guides of the side's donors, all cells) ===")
half_rows = []
percond = {}
for cond in CONDS:
    # targets present in all four donors
    tt = [t for t in ALLT if all((cond, d, t) in IDX for d in DORDER)]
    cells = {t: {d: float(IDX[(cond, d, t)].n_cells.sum()) for d in DORDER} for t in tt}
    store = {}
    for (a, b_) in PAIRINGS:
        rowsA = {t: np.concatenate([IDX[(cond, d, t)].i.values for d in a]) for t in tt}
        rowsB = {t: np.concatenate([IDX[(cond, d, t)].i.values for d in b_]) for t in tt}
        A = resid(side_delta(cond, a, rowsA, tt))
        B = resid(side_delta(cond, b_, rowsB, tt))
        cs = cosines(A, B)
        store["|".join(["".join(a), "".join(b_)])] = cs
        cA = np.array([cells[t][a[0]] + cells[t][a[1]] for t in tt])
        cB = np.array([cells[t][b_[0]] + cells[t][b_[1]] for t in tt])
        run_conf = (set(a) == {"D1", "D2"}) and cond in ("Rest", "Stim8hr")
        m, lo, hi, n = boot_median(cs)
        half_rows.append(dict(cond=cond, pairing=f"{'+'.join(a)} | {'+'.join(b_)}", n=n, median=m, lo=lo, hi=hi,
                              mean=float(np.nanmean(cs)), frac_le0=float(np.nanmean(cs <= 0)),
                              cells_A_med=float(np.median(cA)), cells_B_med=float(np.median(cB)),
                              run_confounded=bool(run_conf)))
        P(f"  {cond:<9} {'+'.join(a)}|{'+'.join(b_)}  n={n:3d}  median {m:+.4f} [{lo:+.4f},{hi:+.4f}]"
          f"  mean {np.nanmean(cs):+.4f}  <=0 {np.nanmean(cs<=0):.3f}"
          f"  cells/side {np.median(cA):.0f}/{np.median(cB):.0f}{'   RUN-CONFOUNDED' if run_conf else ''}")
    percond[cond] = dict(targets=tt, cos=store)
    pooled = np.concatenate([store[k] for k in store])
    m, lo, hi, n = boot_median(pooled)
    P(f"  {cond:<9} POOLED over the 3 pairings: n={n} median {m:+.4f} [{lo:+.4f},{hi:+.4f}]  <=0 {np.nanmean(pooled<=0):.3f}")
    half_rows.append(dict(cond=cond, pairing="POOLED(3 pairings)", n=n, median=m, lo=lo, hi=hi,
                          mean=float(np.nanmean(pooled)), frac_le0=float(np.nanmean(pooled <= 0)),
                          cells_A_med=np.nan, cells_B_med=np.nan, run_confounded=False))
pd.DataFrame(half_rows).to_csv(f"q1_halves_{TAG}.csv", index=False)

# P4: run-confounded vs run-balanced pairings
P("\n--- P4: run-confounded vs run-balanced pairings (paired over targets) ---")
p4 = []
for cond in CONDS:
    st = percond[cond]["cos"]
    conf = st["D1D2|D3D4"]
    bal = np.nanmean(np.vstack([st["D1D3|D2D4"], st["D1D4|D2D3"]]), axis=0)
    d, lo, hi, n = boot_median_diff(conf, bal)
    P(f"  {cond:<9} median(conf) - median(balanced) = {d:+.4f} [{lo:+.4f},{hi:+.4f}]  n={n}"
      f"   ({'donor==run' if cond != 'Stim48hr' else 'all four donors in one run'})")
    p4.append(dict(cond=cond, diff=d, lo=lo, hi=hi, n=n))
OUT["P4"] = p4

# P5 / KA-4: the label-shuffled null, on the primary condition
P("\n--- P5 / KA-4: label-shuffled null (Stim48hr, pairing D1D2|D3D4) ---")
cond = "Stim48hr"
tt = percond[cond]["targets"]
rowsA = {t: np.concatenate([IDX[(cond, d, t)].i.values for d in ("D1", "D2")]) for t in tt}
rowsB = {t: np.concatenate([IDX[(cond, d, t)].i.values for d in ("D3", "D4")]) for t in tt}
A = resid(side_delta(cond, ("D1", "D2"), rowsA, tt))
B = resid(side_delta(cond, ("D3", "D4"), rowsB, tt))
rng = np.random.default_rng(SEED)
perm = rng.permutation(len(tt))
while np.any(perm == np.arange(len(tt))):
    perm = rng.permutation(len(tt))
null = cosines(A, B[perm])
real = cosines(A, B)
P(f"  real   median {np.nanmedian(real):+.4f}   n={np.isfinite(real).sum()}")
P(f"  null   median {np.nanmedian(null):+.4f}  mean {np.nanmean(null):+.4f}  sd {np.nanstd(null):.4f}")
OUT["P5"] = dict(real_median=float(np.nanmedian(real)), null_median=float(np.nanmedian(null)),
                 null_mean=float(np.nanmean(null)), null_sd=float(np.nanstd(null)),
                 pass_=bool(abs(np.nanmean(null)) < 0.02))

# ==================================================================== 3.1 the 2x2
P("\n=== 3.1 the 2x2 at a common counting budget (D floor / A construct / B donor / C both) ===")
rows_2x2 = []
percell = {}
for cond in CONDS:
    # targets with 2 guides in >= 2 donors
    ok = [t for t in ALLT if sum(1 for d in DORDER if (cond, d, t) in IDX and len(IDX[(cond, d, t)]) >= 2) >= 2]
    P(f"  {cond}: {len(ok)} targets with 2 guides in >= 2 donors")
    # per target the available (donor, guide-rank) cells, and the common budget
    cellmap = {}
    budget = {}
    for t in ok:
        cm = {}
        for d in DORDER:
            g = IDX.get((cond, d, t))
            if g is not None and len(g) >= 2:
                cm[d] = [(int(g.i.values[0]), float(g.total.values[0]), float(g.n_cells.values[0])),
                         (int(g.i.values[1]), float(g.total.values[1]), float(g.n_cells.values[1]))]
        cellmap[t] = cm
        budget[t] = 0.5 * min(v[1] for lst in cm.values() for v in lst)
    percell[cond] = (ok, cellmap, budget)

    rng = np.random.default_rng(SEED)
    # the common centring subset: targets with 2 guides in ALL FOUR donors.  Every donor's matrix is
    # centred on the SAME rows, so a different-donor pair is not penalised by a different centring
    # (which would bias Delta_donor = A - C in the direction of the hypothesis).
    common = np.asarray([k for k, t in enumerate(ok) if len(cellmap[t]) == len(DORDER)])
    P(f"    common centring subset (2 guides in all 4 donors): {len(common)} targets")

    full, h1, h2 = {}, {}, {}
    for d in DORDER:
        b = bank.ctrl_profile(cond, [d])[CC]
        for r in (0, 1):
            Mf = np.full((len(ok), len(CC)), np.nan)
            M1 = np.full((len(ok), len(CC)), np.nan)
            M2 = np.full((len(ok), len(CC)), np.nan)
            for k, t in enumerate(ok):
                cm = cellmap[t].get(d)
                if cm is None:
                    continue
                S = bank.counts([cm[r][0]])
                Si = np.rint(S).astype(np.int64)
                s1 = rng.binomial(Si, 0.5).astype(np.float64)      # disjoint halves of the SAME counts
                s2 = Si.astype(np.float64) - s1
                Mf[k] = bank.profile(thin(S, budget[t], rng))[CC] - b
                M1[k] = bank.profile(thin(s1, budget[t], rng))[CC] - b
                M2[k] = bank.profile(thin(s2, budget[t], rng))[CC] - b
            full[(d, r)], h1[(d, r)], h2[(d, r)] = Mf, M1, M2

    def R(M):
        r = np.full((len(ok), int(SC.sum())), np.nan)
        mu = np.nanmean(M[common], axis=0)                          # the COMMON centring
        m = np.isfinite(M[:, 0])
        rr = M[m] - mu
        rr[:, ~SC] = 0.0
        r[m] = rr[:, SC]
        return r

    fullR = {k: R(v) for k, v in full.items()}
    h1R = {k: R(v) for k, v in h1.items()}
    h2R = {k: R(v) for k, v in h2.items()}

    cls = {"D": [], "A": [], "B": [], "C": []}
    for d in DORDER:
        for r in (0, 1):
            cls["D"].append(cosines(h1R[(d, r)], h2R[(d, r)]))
        cls["A"].append(cosines(fullR[(d, 0)], fullR[(d, 1)]))
    for i, d in enumerate(DORDER):
        for d2 in DORDER[i + 1:]:
            for r in (0, 1):
                cls["B"].append(cosines(fullR[(d, r)], fullR[(d2, r)]))
            cls["C"].append(cosines(fullR[(d, 0)], fullR[(d2, 1)]))
            cls["C"].append(cosines(fullR[(d, 1)], fullR[(d2, 0)]))
    agg = {k: np.nanmean(np.vstack(v), axis=0) for k, v in cls.items()}   # per target, averaged within class
    for k in "DABC":
        m, lo, hi, n = boot_median(agg[k])
        P(f"    class {k}: n={n:3d} median {m:+.4f} [{lo:+.4f},{hi:+.4f}]  mean {np.nanmean(agg[k]):+.4f}")
        rows_2x2.append(dict(cond=cond, cls=k, n=n, median=m, lo=lo, hi=hi, mean=float(np.nanmean(agg[k]))))
    for lab, x, y in (("Δ_donor = A - C", agg["A"], agg["C"]), ("D - B (donor beyond counting)", agg["D"], agg["B"]),
                      ("D - A (construct beyond counting)", agg["D"], agg["A"]),
                      ("B - C (construct, given donors differ)", agg["B"], agg["C"])):
        d, lo, hi, n = boot_median_diff(x, y)
        star = "  <-- PRIMARY" if lab.startswith("Δ_donor") and cond == "Stim48hr" else ""
        P(f"    {lab:<38} {d:+.4f} [{lo:+.4f},{hi:+.4f}]  n={n}{star}")
        rows_2x2.append(dict(cond=cond, cls=lab, n=n, median=d, lo=lo, hi=hi, mean=np.nan))
    np.savez(f"q1_2x2_{cond}_{TAG}.npz", targets=np.asarray(ok), **{k: agg[k] for k in "DABC"},
             budget=np.asarray([budget[t] for t in ok]))
pd.DataFrame(rows_2x2).to_csv(f"q1_2x2_{TAG}.csv", index=False)

# ==================================================================== P3 conservation
P("\n=== P3: across-donor reliability vs 276's conservation ===")
c276 = pd.read_csv(f"{ROOT}/experiments/276-conservation-map/work/conservation276.csv").set_index("target")
p3 = []
for cond in CONDS:
    tt = percond[cond]["targets"]
    pooled = np.nanmean(np.vstack([percond[cond]["cos"][k] for k in percond[cond]["cos"]]), axis=0)
    for col in ("C_raw", "C_rel"):
        v = c276.reindex(tt)[col].values.astype(float)
        m = np.isfinite(v) & np.isfinite(pooled)
        rho, p = spearmanr(v[m], pooled[m])
        P(f"  {cond:<9} Spearman(across-donor cos, {col}) = {rho:+.3f} (p={p:.2g}, n={m.sum()})")
        p3.append(dict(cond=cond, col=col, rho=float(rho), p=float(p), n=int(m.sum())))
    v = c276.reindex(tt)["C_raw"].values.astype(float)
    m = np.isfinite(v)
    q = np.quantile(v[m], [1 / 3, 2 / 3])
    ter = np.digitize(v, q)
    for k in range(3):
        s = (ter == k) & m & np.isfinite(pooled)
        P(f"      C_raw tercile {k}: n={s.sum():3d} median cos {np.nanmedian(pooled[s]):+.4f}")
        p3.append(dict(cond=cond, col=f"C_raw tercile {k}", rho=float(np.nanmedian(pooled[s])), p=np.nan, n=int(s.sum())))
    np.savez(f"q1_pertarget_{cond}_{TAG}.npz", targets=np.asarray(tt), cos=pooled,
             **{k: v for k, v in percond[cond]["cos"].items()})
pd.DataFrame(p3).to_csv(f"q1_conservation_{TAG}.csv", index=False)

OUT["elapsed_s"] = time.time() - t0
json.dump(OUT, open(f"q1_280_{TAG}.json", "w"), indent=1)
P(f"\n[{time.time()-t0:.0f}s]")

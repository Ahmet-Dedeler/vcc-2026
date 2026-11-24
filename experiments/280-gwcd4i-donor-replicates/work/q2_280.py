"""280 Q2 -- replicates vs cells: leave-one-donor-out, k in {1,2,3} training donors, at full cells and at a
matched TOTAL counting budget.  Preregistered in ../PREREGISTRATION.md §4 (P6, P7)."""
from __future__ import annotations

import itertools
import json
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from lib280 import CONDS, DORDER, SEED, Bank, boot_median, boot_median_diff, cosines, residualise, thin  # noqa: E402

ROOT = "~/Code/virtual-cell-challenge-2026"
KEEP = "--all" not in sys.argv
TAG = "keepDE" if KEEP else "allguides"
P = lambda *a: print(*a, flush=True)  # noqa: E731
t0 = time.time()
bank = Bank(".")
SC, CC = bank.scored, bank.cols_cd4

o = bank.obs
tg = o[o.gtype == "targeting"]
if KEEP:
    tg = tg[tg.keep_for_DE]
IDX = {(c, d, t): g.sort_values("guide") for (c, d, t), g in tg.groupby(["cond", "D", "pert"])}
ALLT = sorted(set(tg.pert))
rows = []
per = {}

for cond in CONDS:
    tt = [t for t in ALLT if all((cond, d, t) in IDX for d in DORDER)]
    tot = {t: {d: float(IDX[(cond, d, t)].total.sum()) for d in DORDER} for t in tt}
    cells = {t: {d: float(IDX[(cond, d, t)].n_cells.sum()) for d in DORDER} for t in tt}
    budget = {t: min(tot[t].values()) for t in tt}         # the k=1 counting budget, per target
    rng = np.random.default_rng(SEED)

    def donor_resid(d, share):
        """Residual rows for donor d, each target thinned to `share` * budget[t] counts (share=None: full)."""
        b = bank.ctrl_profile(cond, [d])[CC]
        M = np.zeros((len(tt), len(CC)))
        for k, t in enumerate(tt):
            S = bank.counts(IDX[(cond, d, t)].i.values)
            if share is not None:
                S = thin(S, share * budget[t], rng)
            M[k] = bank.profile(S)[CC] - b
        R = residualise(M)
        R[:, ~SC] = 0.0
        return R[:, SC]

    full = {d: donor_resid(d, None) for d in DORDER}
    matched = {k: {d: donor_resid(d, 1.0 / k) for d in DORDER} for k in (1, 2, 3)}

    for mode, src in (("full", lambda k, d: full[d]), ("matched", lambda k, d: matched[k][d])):
        curve = {}
        for k in (1, 2, 3):
            acc = []
            for h in DORDER:                                   # held-out donor = the truth
                tr = [d for d in DORDER if d != h]
                truth = full[h]
                for sub in itertools.combinations(tr, k):
                    pred = np.mean([src(k, d) for d in sub], axis=0)
                    acc.append(cosines(pred, truth))
            curve[k] = np.nanmean(np.vstack(acc), axis=0)
            m, lo, hi, n = boot_median(curve[k])
            P(f"  {cond:<9} {mode:<8} k={k}: n={n:3d} median cos-to-held-donor {m:+.4f} [{lo:+.4f},{hi:+.4f}]"
              f"   mean {np.nanmean(curve[k]):+.4f}")
            rows.append(dict(cond=cond, mode=mode, k=k, n=n, median=m, lo=lo, hi=hi,
                             mean=float(np.nanmean(curve[k]))))
        for a, b_ in ((3, 1), (2, 1), (3, 2)):
            d, lo, hi, n = boot_median_diff(curve[a], curve[b_])
            star = "  <-- P6" if (mode == "matched" and (a, b_) == (3, 1)) else ""
            P(f"  {cond:<9} {mode:<8} Δcos k={a} - k={b_}: {d:+.4f} [{lo:+.4f},{hi:+.4f}]  n={n}{star}")
            rows.append(dict(cond=cond, mode=mode, k=f"{a}-{b_}", n=n, median=d, lo=lo, hi=hi, mean=np.nan))
        per[(cond, mode)] = curve
    med_cells = np.median([sum(cells[t].values()) / 4 for t in tt])
    P(f"  {cond}: {len(tt)} targets in all 4 donors; median cells per donor {med_cells:.0f}; "
      f"median k=1 count budget {np.median([budget[t] for t in tt]):.3g}")
    np.savez(f"q2_{cond}_{TAG}.npz", targets=np.asarray(tt),
             **{f"{m}_k{k}": per[(cond, m)][k] for m in ("full", "matched") for k in (1, 2, 3)})

pd.DataFrame(rows).to_csv(f"q2_curve_{TAG}.csv", index=False)

# --- by 276 conservation tercile, matched mode, primary condition
c276 = pd.read_csv(f"{ROOT}/experiments/276-conservation-map/work/conservation276.csv").set_index("target")
out = []
for cond in CONDS:
    z = np.load(f"q2_{cond}_{TAG}.npz", allow_pickle=True)
    tt = [str(x) for x in z["targets"]]
    v = c276.reindex(tt)["C_raw"].values.astype(float)
    q = np.nanquantile(v, [1 / 3, 2 / 3])
    ter = np.digitize(v, q)
    for k in range(3):
        s = (ter == k) & np.isfinite(v)
        d, lo, hi, n = boot_median_diff(z["matched_k3"][s], z["matched_k1"][s])
        P(f"  {cond:<9} C_raw tercile {k}: n={n:3d} Δcos(k=3 - k=1, matched) {d:+.4f} [{lo:+.4f},{hi:+.4f}]"
          f"   median cos k=1 {np.nanmedian(z['matched_k1'][s]):+.4f}")
        out.append(dict(cond=cond, tercile=k, n=n, dcos=d, lo=lo, hi=hi,
                        cos_k1=float(np.nanmedian(z["matched_k1"][s]))))
pd.DataFrame(out).to_csv(f"q2_terciles_{TAG}.csv", index=False)
P(f"\n[{time.time()-t0:.0f}s]")

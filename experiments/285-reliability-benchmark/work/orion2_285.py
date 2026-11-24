"""285 -- the Orion / K562 BATCH-vs-CELLS contrast at a MATCHED UMI BUDGET, from 267's cached integer
pseudobulk sums (`work/half/<line>/acc_{batchA,batchB,randA}.npy`, mmapped; randB = batchA + batchB - randA).

267 found that on Orion HEK293T a batch split reads half a random split at the same cells; every other line
in this benchmark shows no batch term at all.  This row checks whether that exception survives a ladder in
the budget, and adds the guide row (`guide_accA`, the two-construct perturbations) on the same axis.
No new pass over the 118 GB of parquet is made and no raw cell is read.

    PYTHONDONTWRITEBYTECODE=1 uv run python orion2_285.py
"""
from __future__ import annotations
import csv, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib285 as L

P = L.Logger(L.WORK / "orion2_285.log")
E267 = L.ROOT / "experiments/267-replicate-validated-oracle/work/half"
NREP = 4


def main():
    t0 = time.time()
    P(f"285 orion2_285.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    with open(L.ROOT / "data/challenge/pert_counts.csv") as fh:
        PAN = {r["target_gene"] for r in csv.DictReader(fh)}
    out = {}
    for line in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        d = E267 / line
        z = np.load(d / "sums.npz", allow_pickle=True)
        names = np.asarray([str(x) for x in z["pert_names"]])
        keep_names = np.asarray([str(x) for x in z["keep_names"]])
        ncell = np.asarray(z["ncell"])          # (n_pert, 3) = full / batchA / batchB? -- halves key names it
        halves = [str(x) for x in z["halves"]]
        rows = np.flatnonzero(np.isin(names, sorted(PAN)))
        A = np.load(d / "acc_batchA.npy", mmap_mode="r")[rows].astype(np.float64)
        B = np.load(d / "acc_batchB.npy", mmap_mode="r")[rows].astype(np.float64)
        RA = np.load(d / "acc_randA.npy", mmap_mode="r")[rows].astype(np.float64)
        RB = A + B - RA
        tg = names[rows]
        gset = {g: i for i, g in enumerate(keep_names)}
        scored = np.asarray(sorted(set(range(len(keep_names))) - {gset[t] for t in tg if t in gset}), dtype=np.int64)
        tot = np.stack([X.sum(1) for X in (A, B, RA, RB)])
        P(f"  {line}: halves keys {halves}; {rows.size} panel targets, {len(keep_names)} mapped genes "
          f"({scored.size} scored); side totals median {np.median(tot, 1).round(0)}")
        out[line] = dict(n_targets=int(rows.size), budgets={})
        for T in (2e5, 5e5, 1e6):
            ok = np.flatnonzero((tot >= T).all(0))
            if ok.size < 20:
                continue
            cb = np.full((NREP, ok.size), np.nan); cr = np.full((NREP, ok.size), np.nan)
            for rep in range(NREP):
                r = np.random.default_rng(L.SEED + 31 * rep + int(T))
                M = [np.stack([L.thin(X[i], T / tot[k, i], r) for i in ok]) for k, X in enumerate((A, B, RA, RB))]
                Pr = [L.prof(m)[:, scored] for m in M]
                Pr = [q - q.mean(0, keepdims=True) for q in Pr]
                cb[rep] = L.row_cos(Pr[0], Pr[1])
                cr[rep] = L.row_cos(Pr[2], Pr[3])
            xb, xc = np.nanmean(cb, 0), np.nanmean(cr, 0)
            rng = np.random.default_rng(L.SEED)
            m = np.isfinite(xb) & np.isfinite(xc)
            ii = rng.integers(0, int(m.sum()), size=(L.NBOOT, int(m.sum())))
            rb = np.median(xb[m][ii], 1) / np.median(xc[m][ii], 1)
            out[line]["budgets"][f"{T:.0e}"] = dict(
                n=int(ok.size), batch=dict(**L.stats(xb), boot=L.boot_median(xb)),
                rand=dict(**L.stats(xc), boot=L.boot_median(xc)),
                rand_minus_batch=L.boot_median_diff(xc, xb),
                rho=dict(point=float(np.nanmedian(xb) / np.nanmedian(xc)),
                         lo=float(np.percentile(rb, 2.5)), hi=float(np.percentile(rb, 97.5))))
            v = out[line]["budgets"][f"{T:.0e}"]
            P(f"    T={T:.0e} UMIs/side, {ok.size:3d} targets: RANDOM {L.fmt(v['rand']['boot'])}  BATCH "
              f"{L.fmt(v['batch']['boot'])}  RANDOM − BATCH {L.fmt(v['rand_minus_batch'])}  rho "
              f"{v['rho']['point']:+.3f} [{v['rho']['lo']:+.3f}, {v['rho']['hi']:+.3f}]")
        del A, B, RA, RB
    out["seconds"] = round(time.time() - t0, 1); out["peak_rss_gib"] = L.peak_gib()
    L.Logger.dump(P, out, L.WORK / "orion2_285.json")
    P(f"\ndone {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB"); P.save()


if __name__ == "__main__":
    main()

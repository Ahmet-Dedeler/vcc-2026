"""285 -- the CD4 (GWCD4i) rows of the benchmark, from 280's (guide, donor, condition) pseudobulk bank.

Units the bank can carry (it has no cells):
  (b) different GUIDE construct, same target, same donor, same condition
  (c) different DONOR, same guide, same target, same condition -- FIXED donor pairs
  (s) different STATE, same donor, same guide -- FIXED condition pairs
  (a') the UMI budget: every side is thinned to the same total-UMI budget T, laddered over T
  (a*) the molecular COUNT-SPLIT of one pseudobulk into two independent binomial halves -- the Poisson part
       of the cell-sampling term, a LOWER BOUND on a true cell split (`csplit_285.py` calibrates the gap on
       h1-test and K562, where both exist).

CENTRING.  Every side is centred on ONE common reference: the mean, over the same targets, of the pooled
(all donors, all guides) profile of that condition thinned to the same budget T.  A first run centred each
side on its OWN mean and the donor row came out at -0.05 with 70 % of targets below zero: when side 0 is
donor A for some targets and donor B for others, a side-own mean is a donor MIXTURE, and subtracting it puts
+/- half the donor difference into the two sides with opposite signs -- a pure artefact of the assignment.
Fixed pairs + a common reference remove it.

    PYTHONDONTWRITEBYTECODE=1 uv run python cd4_285.py
"""
from __future__ import annotations

import sys
import time
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/280-gwcd4i-donor-replicates/work")
import lib285 as L  # noqa: E402
import lib280 as B280  # noqa: E402

P = L.Logger(L.WORK / "cd4_285.log")
FAIL: list[str] = []
E280 = "~/Code/virtual-cell-challenge-2026/experiments/280-gwcd4i-donor-replicates/work"
BUDGETS = (1e5, 2e5, 4e5)
NREP = 4
CONDS = ("Rest", "Stim8hr", "Stim48hr")


def main() -> None:
    t0 = time.time()
    P(f"285 cd4_285.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    bank = B280.Bank(E280)
    obs = bank.obs
    cc = bank.cols_cd4[bank.scored]
    panel_set = set(bank.panel)
    IDX = {(c, d, t): g for (c, d, t), g in obs[obs.gtype == "targeting"].groupby(["cond", "D", "pert"])}
    P(f"  bank {bank.X.shape}, {int(bank.scored.sum())} scored columns, {len(bank.panel)} panel targets")

    # ---------------------------------------------------------------- KA-3: 280's Stim48hr donor-half median
    COND = "Stim48hr"
    panel = [t for t in bank.panel if all((COND, d, t) in IDX for d in B280.DORDER)]
    acc = []
    for pr in B280.PAIRINGS:
        A = np.stack([bank.profile(bank.counts(np.concatenate([IDX[(COND, d, t)].i.values for d in pr[0]])))[cc]
                      for t in panel])
        Bm = np.stack([bank.profile(bank.counts(np.concatenate([IDX[(COND, d, t)].i.values for d in pr[1]])))[cc]
                       for t in panel])
        acc.append(L.row_cos(L.resid(A), L.resid(Bm)))
    dd = np.nanmean(np.vstack(acc), 0)
    got = float(np.nanmedian(dd))
    ok = abs(got - 0.0539) <= 1.5e-3
    ka = [dict(row="KA-3 280's Stim48hr donor-half median (3 disjoint 2v2 pairings, panel-centred, full depth)",
               required="|d| <= 1.5e-3 vs +0.0539", want=0.0539, got=got, ok=bool(ok), n=int(np.isfinite(dd).sum()))]
    P(f"KA-3: 280's Stim48hr donor-half median recomputed from the bank = {got:+.4f} (280: +0.0539) "
      f"{'PASS' if ok else 'FAIL'}, n {np.isfinite(dd).sum()}, {len(panel)} targets in all four donors")
    if not ok:
        FAIL.append("KA-3")

    tot = obs.total.values.astype(np.float64)
    ncl = obs.n_cells.values.astype(np.float64)
    gid = obs.guide.values.astype(str)
    X = np.asarray(bank.X[:, bank.cols_cd4[bank.scored]], dtype=np.float64)
    P(f"  scored bank materialised {X.shape}, peak {L.peak_gib():.2f} GiB, {time.time()-t0:.0f}s; "
      f"targeting row totals median {np.median(tot[obs.gtype.values == 'targeting']):.3g}, cells median "
      f"{np.median(ncl[obs.gtype.values == 'targeting']):.0f}")

    rows_ = {}          # (cond, donor, target) -> list of bank rows, sorted by guide id
    for (c, d, t), g in IDX.items():
        if t in panel_set:
            rows_[(c, d, t)] = sorted(g.i.values.tolist(), key=lambda i: gid[i])

    out: dict = {"ka": ka, "budgets": {}}
    arrays: dict[str, np.ndarray] = {}

    # Every contrast is a FIXED (donor, guide-slot, condition) triple on BOTH sides, and every side is
    # centred on its OWN mean over the targets of that cell.  A side-specific constant vector (the donor's
    # global deviation, a guide slot's global deviation) is then removed rather than left in the residual:
    # leaving it in put a common +delta_d in every target of a side, which is what made the first two runs
    # read donor = -0.05 (the four donors' deviations sum to zero, so two of them anti-correlate) and
    # guide ABOVE the count split (both guide sides share one donor's deviation).
    def slot(cond, d, t, k):
        rs = rows_.get((cond, d, t), [])
        return [rs[k]] if len(rs) > k else None

    for T in BUDGETS:
        P(f"\n  === budget T = {T:.0e} UMIs per side ===")
        key_T = f"{T:.0e}"
        out["budgets"][key_T] = {}
        cells_of: dict[str, list] = defaultdict(list)
        store: dict[str, list] = defaultdict(list)
        for rep in range(NREP):
            r = np.random.default_rng(L.SEED + 7919 * rep + int(T))

            def at_budget(rowset):
                if rowset is None:
                    return None, None
                s = tot[rowset].sum()
                if s < T:
                    return None, None
                S = X[rowset].sum(0) if len(rowset) > 1 else X[rowset[0]]
                return L.thin(S, T / s, r), float(ncl[rowset].sum() * T / s)

            def cell(tag, pairs):
                """pairs: (target, rowsetA, rowsetB) with a fixed side convention; per-side centring."""
                keep, A, Bm, nc = [], [], [], []
                for t, ra, rb in pairs:
                    va, na = at_budget(ra); vb, nb = at_budget(rb)
                    if va is None or vb is None:
                        continue
                    keep.append(t); A.append(va); Bm.append(vb); nc.append((na + nb) / 2)
                if len(keep) < 15:
                    return
                PA, PB = L.prof(np.stack(A)), L.prof(np.stack(Bm))
                c = L.row_cos(PA - PA.mean(0, keepdims=True), PB - PB.mean(0, keepdims=True))
                store[tag].append((keep, c))
                cells_of[tag].append(float(np.median(nc)))

            for cond in CONDS:
                tset = sorted({t for (c, d, t) in rows_ if c == cond})
                for d in B280.DORDER:
                    # (a*) the count split of ONE (donor, guide-slot-0) pseudobulk: the same construction as
                    # every other side, split into two independent binomial halves at budget T
                    pairs = []
                    for t in tset:
                        rs = slot(cond, d, t, 0)
                        if rs is None or tot[rs].sum() < 2 * T:
                            continue
                        h = L.thin(X[rs[0]], 2 * T / tot[rs].sum(), r)
                        a = r.binomial(np.rint(h).astype(np.int64), 0.5).astype(np.float64)
                        pairs.append((t, a, h - a))
                    if len(pairs) >= 15:
                        keep = [t for t, _, _ in pairs]
                        PA = L.prof(np.stack([a for _, a, _ in pairs]))
                        PB = L.prof(np.stack([b for _, _, b in pairs]))
                        c = L.row_cos(PA - PA.mean(0, keepdims=True), PB - PB.mean(0, keepdims=True))
                        store[f"{cond}|countsplit:{d}"].append((keep, c))
                        cells_of[f"{cond}|countsplit:{d}"].append(
                            float(np.median([ncl[slot(cond, d, t, 0)].sum() * T / tot[slot(cond, d, t, 0)].sum()
                                             for t in keep])))
                    # (b) guide slots 0 vs 1 inside ONE donor
                    cell(f"{cond}|guide:{d}", [(t, slot(cond, d, t, 0), slot(cond, d, t, 1)) for t in tset
                                               if slot(cond, d, t, 1) is not None])
                # (c) donor: fixed pairs, guide slot 0 on both sides, same guide id required
                for da, db in combinations(B280.DORDER, 2):
                    pr = []
                    for t in tset:
                        ra, rb = slot(cond, da, t, 0), slot(cond, db, t, 0)
                        if ra is not None and rb is not None and gid[ra[0]] == gid[rb[0]]:
                            pr.append((t, ra, rb))
                    cell(f"{cond}|donor:{da}|{db}", pr)
                # (s) state: fixed condition pairs, one donor, guide slot 0, same guide id
                for other in CONDS:
                    if other <= cond:
                        continue
                    for d in B280.DORDER:
                        pr = []
                        for t in tset:
                            ra, rb = slot(cond, d, t, 0), slot(other, d, t, 0)
                            if ra is not None and rb is not None and gid[ra[0]] == gid[rb[0]]:
                                pr.append((t, ra, rb))
                        cell(f"{cond}|state:{other}:{d}", pr)

        # aggregate reps, then collapse the donor / guide / state families
        agg: dict[str, tuple] = {}
        for tag, lst in store.items():
            common = sorted(set.intersection(*[set(k) for k, _ in lst]))
            if len(common) < 15:
                continue
            M = np.vstack([[c[k.index(t)] for t in common] for k, c in lst])
            agg[tag] = (common, np.nanmean(M, 0))
        for cond in CONDS:
            for fam in ("countsplit", "guide", "donor", "state"):
                keys = [k for k in agg if k.startswith(f"{cond}|{fam}:")]
                if not keys:
                    continue
                common = sorted(set.intersection(*[set(agg[k][0]) for k in keys]))
                if len(common) < 15:
                    continue
                M = np.vstack([[agg[k][1][agg[k][0].index(t)] for t in common] for k in keys])
                x = np.nanmean(M, 0)
                key = f"{cond}|{fam}"
                arrays[f"{key}|{key_T}"] = x
                out["budgets"][key_T][key] = dict(**L.stats(x), boot=L.boot_median(x), targets=common,
                                                  n_cells_families=len(keys),
                                                  cells_per_side=float(np.mean([np.mean(cells_of[k]) for k in keys
                                                                                if cells_of[k]])) if any(cells_of[k] for k in keys) else None)
        for key, v in sorted(out["budgets"][key_T].items()):
            P(f"    {key:26s} n_t {v['n']:3d}  ~{(v['cells_per_side'] or 0):5.1f} cells/side  {v['n_cells_families']} cells  "
              f"median {L.fmt(v['boot'])}  <=0 {v['share_le0']:.3f}")
        for cond in CONDS:
            base = f"{cond}|countsplit"
            if base not in out["budgets"][key_T]:
                continue
            for tag in ("guide", "donor", "state"):
                key = f"{cond}|{tag}"
                if key not in out["budgets"][key_T]:
                    continue
                ta, tb = out["budgets"][key_T][base]["targets"], out["budgets"][key_T][key]["targets"]
                com = sorted(set(ta) & set(tb))
                if len(com) < 15:
                    continue
                a = arrays[f"{base}|{key_T}"][[ta.index(t) for t in com]]
                b = arrays[f"{key}|{key_T}"][[tb.index(t) for t in com]]
                d = L.boot_median_diff(a, b)
                out["budgets"][key_T][f"{key}|vs_countsplit"] = dict(
                    n=len(com), diff=d, ratio=float(np.nanmedian(b) / max(np.nanmedian(a), 1e-9)))
                P(f"    COUNT-SPLIT − {key:22s} on {len(com):3d} shared targets: {L.fmt(d)}  "
                  f"(ratio {np.nanmedian(b)/max(np.nanmedian(a),1e-9):+.3f})")

    out["failures"] = FAIL
    out["seconds"] = round(time.time() - t0, 1)
    out["peak_rss_gib"] = L.peak_gib()
    np.savez_compressed(L.WORK / "cd4_285.npz", **arrays)
    L.Logger.dump(P, out, L.WORK / "cd4_285.json")
    P(f"\ndone: {len(FAIL)} KA failures, {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB")
    P.save()
    if FAIL:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

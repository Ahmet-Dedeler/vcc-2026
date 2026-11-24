"""285 -- the Broad-16 rows: (a) random cell halves vs (b) the two guide constructs, on 12 working lines,
98 targets each, from 268's cached `work/halves/<line>.npz` (per-target count sums for full / A / B / v1 / v2).

The halves are sums, so cells cannot be re-sampled; every side is instead thinned to a common TOTAL-UMI
budget T, which is the currency in which (a) and (b) are matched everywhere else in this benchmark.
(f) different LINES is read off the same object: the full-depth response of one line against another's.

    PYTHONDONTWRITEBYTECODE=1 uv run python broad_285.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib285 as L  # noqa: E402

P = L.Logger(L.WORK / "broad_285.log")
HV = L.ROOT / "experiments/268-broad16-context-residuals/work/halves"
BUDGETS = (2e5, 5e5)
NREP = 4


def main() -> None:
    t0 = time.time()
    P(f"285 broad_285.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    lines = sorted(p.stem for p in HV.glob("*.npz"))
    out: dict = {"lines": lines, "ka": [], "budgets": {}}
    data = {}
    for ln in lines:
        z = np.load(HV / f"{ln}.npz", allow_pickle=True)
        d = {k: np.asarray(z[k]) for k in ("targets", "genes", "ncell", "libmean", "full", "A", "B", "v1", "v2")}
        d["targets"] = d["targets"].astype(str); d["genes"] = d["genes"].astype(str)
        data[ln] = d
    tg = sorted(set.intersection(*[set(data[ln]["targets"].tolist()) for ln in lines]))
    genes = data[lines[0]]["genes"]
    assert all(list(data[ln]["genes"]) == list(genes) for ln in lines)
    for ln in lines:                                     # re-index every line onto the common target list
        pos = {t: i for i, t in enumerate(data[ln]["targets"])}
        sel = np.asarray([pos[t] for t in tg], dtype=np.int64)
        for k in ("targets", "ncell", "libmean", "full", "A", "B", "v1", "v2"):
            data[ln][k] = data[ln][k][sel]
    gset = {g: i for i, g in enumerate(genes)}
    scored = np.asarray(sorted(set(range(len(genes))) - {gset[t] for t in tg if t in gset}), dtype=np.int64)
    P(f"  {len(lines)} lines, {len(tg)} targets, {len(genes)} genes ({scored.size} scored); "
      f"cells/target mean {np.mean([d['ncell'][:, 0].mean() for d in data.values()]):.1f}")

    # KA-6: A + B == full and v1 + v2 <= full on every line (integer identity on the mapped genes)
    worst_ab, worst_v = 0.0, 0.0
    for ln in lines:
        d = data[ln]
        m = np.isfinite(d["full"]) & np.isfinite(d["A"]) & np.isfinite(d["B"])
        worst_ab = max(worst_ab, float(np.abs(d["A"][m] + d["B"][m] - d["full"][m]).max()))
        mv = np.isfinite(d["v1"]) & np.isfinite(d["v2"]) & np.isfinite(d["full"])
        worst_v = max(worst_v, float((d["v1"][mv] + d["v2"][mv] - d["full"][mv]).max()))
    ok = worst_ab == 0.0 and worst_v <= 0.0
    out["ka"].append(dict(row="KA-6 268's halves satisfy A + B == full and v1 + v2 <= full on all 12 lines",
                          required="0 / <= 0", got=[worst_ab, worst_v], ok=bool(ok)))
    P(f"KA-6: max|A + B - full| = {worst_ab:.1f}; max(v1 + v2 - full) = {worst_v:.1f} -> {'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit(1)

    for T in BUDGETS:
        key_T = f"{T:.0e}"
        out["budgets"][key_T] = {}
        cell_rows, guide_rows, per_line = [], [], {}
        for ln in lines:
            d = data[ln]
            nan = ~np.isfinite(d["full"])
            S = {k: np.nan_to_num(d[k], nan=0.0) for k in ("A", "B", "v1", "v2")}
            tot = {k: S[k].sum(1) for k in S}
            keep = np.flatnonzero((tot["A"] >= T) & (tot["B"] >= T) & (tot["v1"] >= T) & (tot["v2"] >= T))
            if keep.size < 15:
                continue
            ca = np.full((NREP, keep.size), np.nan); cg = np.full((NREP, keep.size), np.nan)
            for rep in range(NREP):
                r = np.random.default_rng(L.SEED + 101 * rep + int(T))
                M = {k: np.stack([L.thin(S[k][i], T / tot[k][i], r) for i in keep]) for k in S}
                PA, PB = L.prof(M["A"])[:, scored], L.prof(M["B"])[:, scored]
                QA, QB = L.prof(M["v1"])[:, scored], L.prof(M["v2"])[:, scored]
                ca[rep] = L.row_cos(PA - PA.mean(0, keepdims=True), PB - PB.mean(0, keepdims=True))
                cg[rep] = L.row_cos(QA - QA.mean(0, keepdims=True), QB - QB.mean(0, keepdims=True))
            xa, xg = np.nanmean(ca, 0), np.nanmean(cg, 0)
            per_line[ln] = dict(n=int(keep.size), cells=float(np.median(d["ncell"][keep, 0])),
                                cell=float(np.nanmedian(xa)), guide=float(np.nanmedian(xg)),
                                diff=L.boot_median_diff(xa, xg))
            cell_rows.append(xa); guide_rows.append(xg)
            P(f"  [T={key_T}] {ln:9s} n_t {keep.size:3d}  cells/target {np.median(d['ncell'][keep, 0]):5.0f}  "
              f"CELL halves {np.nanmedian(xa):+.4f}   GUIDE halves {np.nanmedian(xg):+.4f}   CELL − GUIDE "
              f"{L.fmt(per_line[ln]['diff'])}")
        if cell_rows:
            allc = np.concatenate(cell_rows); allg = np.concatenate(guide_rows)
            out["budgets"][key_T] = dict(per_line=per_line, pooled_cell=dict(**L.stats(allc), boot=L.boot_median(allc)),
                                         pooled_guide=dict(**L.stats(allg), boot=L.boot_median(allg)),
                                         pooled_diff=L.boot_median_diff(allc, allg),
                                         lines_with_positive_diff=int(sum(v["diff"]["point"] > 0 for v in per_line.values())),
                                         n_lines=len(per_line))
            # (f) different LINES at the SAME budget: line i's full row vs line j's, both thinned to T
            LP = {}
            for ln in lines:
                d = data[ln]
                S = np.nan_to_num(d["full"], nan=0.0)
                tt = S.sum(1)
                kp = np.flatnonzero(tt >= T)
                r = np.random.default_rng(L.SEED + 555 + int(T))
                M = np.stack([L.thin(S[i], T / tt[i], r) for i in kp])
                Pm = L.prof(M)[:, scored]
                LP[ln] = (kp, Pm - Pm.mean(0, keepdims=True))
            xl = []
            for i, a in enumerate(lines):
                for b in lines[i + 1:]:
                    com = np.intersect1d(LP[a][0], LP[b][0])
                    if com.size < 15:
                        continue
                    ia = np.searchsorted(LP[a][0], com); ib = np.searchsorted(LP[b][0], com)
                    xl.append(L.row_cos(LP[a][1][ia], LP[b][1][ib]))
            if xl:
                xl = np.concatenate(xl)
                out["budgets"][key_T]["lines_matched"] = dict(**L.stats(xl), boot=L.boot_median(xl))
                P(f"  [T={key_T}] DIFFERENT LINES at the same budget: {L.fmt(L.boot_median(xl))} "
                  f"(n {np.isfinite(xl).sum()} line-pair rows, <=0 {float(np.nanmean(xl <= 0)):.3f})")
            P(f"  [T={key_T}] POOLED over {len(per_line)} lines: CELL {L.fmt(out['budgets'][key_T]['pooled_cell']['boot'])}"
              f"   GUIDE {L.fmt(out['budgets'][key_T]['pooled_guide']['boot'])}   CELL − GUIDE "
              f"{L.fmt(out['budgets'][key_T]['pooled_diff'])}   (positive on "
              f"{out['budgets'][key_T]['lines_with_positive_diff']}/{len(per_line)} lines)")

    # (f) different LINES at full depth: line i's response vs line j's, same targets, panel-centred
    P("\n(f) different LINES, full depth (268's 12 working lines, 98 targets):")
    R = {}
    for ln in lines:
        d = data[ln]
        S = np.nan_to_num(d["full"], nan=0.0)
        ok = S.sum(1) > 0
        Pm = L.prof(S)[:, scored]
        R[ln] = (Pm - Pm[ok].mean(0, keepdims=True), ok)
    cr = []
    for i, a in enumerate(lines):
        for b in lines[i + 1:]:
            m = R[a][1] & R[b][1]
            cr.append(L.row_cos(R[a][0][m], R[b][0][m]))
    x = np.concatenate(cr)
    out["lines_full_depth"] = dict(**L.stats(x), boot=L.boot_median(x), n_pairs=len(cr))
    P(f"  {len(cr)} line pairs, pooled median {L.fmt(out['lines_full_depth']['boot'])}, "
      f"<=0 {out['lines_full_depth']['share_le0']:.3f}")

    out["seconds"] = round(time.time() - t0, 1)
    out["peak_rss_gib"] = L.peak_gib()
    L.Logger.dump(P, out, L.WORK / "broad_285.json")
    P(f"\ndone {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB")
    P.save()


if __name__ == "__main__":
    main()

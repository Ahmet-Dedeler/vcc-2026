"""296 -- the known-answer rows (PREREGISTRATION sec 1).  Nothing new is measured here: every number is
recomputed from 285's and 293's SAVED arrays / csv with this experiment's own code.

    PYTHONDONTWRITEBYTECODE=1 uv run python ka296.py
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib296 as X  # noqa: E402
from lib296 import L  # noqa: E402

P = L.Logger(X.WORK / "ka296.log")
E285 = X.ROOT / "experiments/285-reliability-benchmark/work"
E293 = X.ROOT / "experiments/293-challenge-construct-ceiling/work"
E268 = X.ROOT / "experiments/268-broad16-context-residuals/work/halves"
E280 = X.ROOT / "experiments/280-gwcd4i-donor-replicates/work"
FAIL: list[str] = []
ka: list[dict] = []


def check(row, required, want, got, ok, **kw):
    ka.append(dict(row=row, required=required, want=want, got=got, ok=bool(ok), **kw))
    P(f"{'PASS' if ok else 'FAIL'}  {row}\n        required {required}; got {got}")
    if not ok:
        FAIL.append(row)


def main() -> None:
    t0 = time.time()
    P(f"296 ka296.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")

    # ---------------------------------------------------------------- KA-1: 285's K562 guide row
    z = np.load(E285 / "k562_guide.npz", allow_pickle=True)
    j = json.load(open(E285 / "k562_guide.json"))
    got = {}
    for n in (25, 50):
        xg, xc = z[f"guide|{n}"], z[f"cells|{n}"]
        got[n] = dict(guide=float(np.nanmedian(xg)), cells=float(np.nanmedian(xc)),
                      n=int(np.isfinite(xg).sum()), rho=X.rho_of(xg, xc),
                      pds_guide=float(j["rows"][str(n)]["pds_guide"]),
                      pds_cells=float(j["rows"][str(n)]["pds_cells"]))
    want25 = dict(guide=j["rows"]["25"]["guide"]["median"], cells=j["rows"]["25"]["cells"]["median"],
                  rho=j["rows"]["25"]["rho"]["point"], n=j["rows"]["25"]["guide"]["n"])
    err = max(abs(got[25]["guide"] - want25["guide"]), abs(got[25]["cells"] - want25["cells"]),
              abs(got[25]["rho"] - want25["rho"]))
    check("KA-1  285's K562 guide row at n=25 (748 pairs) recomputed from k562_guide.npz",
          "max |delta| == 0.0e+00 vs the saved medians",
          {k: round(v, 6) for k, v in want25.items()},
          dict(guide=round(got[25]["guide"], 6), cells=round(got[25]["cells"], 6),
               rho=round(got[25]["rho"], 6), n=got[25]["n"], pds_guide=round(got[25]["pds_guide"], 4),
               pds_cells=round(got[25]["pds_cells"], 4), max_abs_err=err),
          err == 0.0, published="(b) +0.0006, (a) +0.0063, pds 0.633 -> 0.509")
    ka[-1]["n50"] = got[50]

    # ---------------------------------------------------------------- KA-2: 285 sec 6b's terciles
    zk = np.load(E285 / "k562_kd.npz", allow_pickle=True)
    jk = json.load(open(E285 / "k562_kd.json"))
    kd, xg, xc = zk["kd"], zk["guide"], zk["cells"]
    weak = np.nanmax(kd, axis=1)                     # the LESS effective construct (285's convention)
    q = np.nanpercentile(weak, [100 / 3, 200 / 3])
    ter = {"T1 weakest KD (least negative)": weak > q[1], "T2": (weak > q[0]) & (weak <= q[1]),
           "T3 strongest KD (most negative)": weak <= q[0]}
    worst, rows = 0.0, {}
    for k, m in ter.items():
        mm = m & np.isfinite(xg) & np.isfinite(xc)
        r = dict(n=int(mm.sum()), weak_kd_median=float(np.nanmedian(weak[mm])),
                 guide=float(np.nanmedian(xg[mm])), cells=float(np.nanmedian(xc[mm])),
                 rho=float(np.nanmedian(xg[mm]) / np.nanmedian(xc[mm])))
        w = jk["strata"][k]
        worst = max(worst, abs(r["n"] - w["n"]), abs(r["weak_kd_median"] - w["weak_kd_median"]),
                    abs(r["guide"] - w["guide"]["median"]), abs(r["cells"] - w["cells"]["median"]),
                    abs(r["rho"] - w["rho"]["point"]))
        rows[k] = r
        P(f"    {k:34s} n {r['n']:3d}  weakKD {r['weak_kd_median']:+.3f}  guide {r['guide']:+.4f}  "
          f"cells {r['cells']:+.4f}  rho {r['rho']:+.3f}")
    check("KA-2  285 sec 6b's three knockdown terciles recomputed from k562_kd.npz",
          "max |delta| == 0.0e+00 across n / weakKD / guide / cells / rho", "285's k562_kd.json",
          dict(max_abs_err=worst, rows=rows), worst == 0.0)
    ka[-1]["kd_medians"] = dict(per_construct=float(np.nanmedian(kd)), weaker=float(np.nanmedian(weak)),
                                both_negative=float(np.mean((kd < 0).all(1))), n=int(kd.shape[0]))

    # ---------------------------------------------------------------- KA-3: 293's three strata rows
    allr = list(csv.DictReader(open(E293 / "construct293.csv")))
    regimes = sorted({r["regime"] for r in allr})
    rows = [r for r in allr if int(r["n"]) == 200 and r["regime"] == "challenge_20k"]
    P(f"    construct293.csv: {len(allr)} rows, regimes {regimes}; using regime=challenge_20k, n=200 -> {len(rows)} rows")
    a = np.asarray([float(r["a_cos"]) for r in rows])
    b = np.asarray([float(r["b_cos"]) for r in rows])
    kdw = np.asarray([min(abs(float(r["kd_big"])), abs(float(r["kd_small"]))) for r in rows])
    sub = {"ALL": np.ones(len(rows), bool), "both >= 1.0": kdw >= 1.0, "weaker < 1.0": kdw < 1.0}
    s293 = json.load(open(E293 / "strat293.json"))["challenge"]["200"]
    want = {"ALL": s293["all"], "both >= 1.0": s293["both_work"], "weaker < 1.0": s293["one_fails"]}
    printed = {"ALL": (0.075, 0.312, 0.023, 16), "both >= 1.0": (0.687, 0.150, 0.103, 5),
               "weaker < 1.0": (-0.003, 0.317, -0.001, 11)}
    worst, got3 = 0.0, {}
    for k, m in sub.items():
        r = dict(k=int(m.sum()), a=float(np.median(a[m])), b=float(np.median(b[m])),
                 rho=X.rho_of(b[m], a[m]))
        w = want[k]
        worst = max(worst, abs(r["rho"] - w["rho"]["point"]), abs(r["a"] - w["a"]["point"]),
                    abs(r["b"] - w["b"]["point"]), abs(r["k"] - w["n_targets"]))
        got3[k] = {kk: (round(v, 9) if isinstance(v, float) else v) for kk, v in r.items()}
        P(f"    {k:14s} k {r['k']:2d}  (a) {r['a']:+.6f}  (b) {r['b']:+.6f}  rho {r['rho']:+.6f}   "
          f"[293 json: (a) {w['a']['point']:+.6f} (b) {w['b']['point']:+.6f} rho {w['rho']['point']:+.6f}]")
    check("KA-3  293's ALL / both->=1.0 / weaker-<1.0 rows re-derived from construct293.csv by THIS code",
          "every cell of 293's strat293.json reproduced to <= 1e-12",
          {k: dict(rho=v["rho"]["point"], a=v["a"]["point"], b=v["b"]["point"], k=v["n_targets"])
           for k, v in want.items()}, got3, worst <= 1e-12,
          published_printed=printed,
          note=("293's printed table shows (a) +0.150 for both->=1.0; the stored value is 0.1494559, "
                "a display rounding only"))

    # ---------------------------------------------------------------- KA-4: Broad-16 half identities
    lines = sorted(p.stem for p in E268.glob("*.npz"))
    wab, wv = 0.0, -np.inf
    for ln in lines:
        d = np.load(E268 / f"{ln}.npz", allow_pickle=True)
        m = np.isfinite(d["full"]) & np.isfinite(d["A"]) & np.isfinite(d["B"])
        wab = max(wab, float(np.abs(d["A"][m] + d["B"][m] - d["full"][m]).max()))
        mv = np.isfinite(d["v1"]) & np.isfinite(d["v2"]) & np.isfinite(d["full"])
        wv = max(wv, float((d["v1"][mv] + d["v2"][mv] - d["full"][mv]).max()))
    check("KA-4  268's halves satisfy A + B == full and v1 + v2 <= full on all 12 lines",
          "max|A+B-full| == 0 and max(v1+v2-full) <= 0", "0 / <= 0",
          dict(max_abs_ab=wab, max_v_minus_full=wv, n_lines=len(lines)), wab == 0.0 and wv <= 0.0)

    # ---------------------------------------------------------------- KA-5: the CD4 guide join
    o = np.load(E280 / "obs280.npz", allow_pickle=True)
    sel = np.load(E280 / "bank280_rows.npy")
    ob = pd.DataFrame({k: o[k][sel].astype(str) for k in ("guide", "pert", "cond", "donor", "gtype")})
    bank_guides = set(ob.loc[ob.gtype == "targeting", "guide"].unique())
    kdt = pd.read_csv(X.ROOT / "data/external/cd4_cell2026/suppl_tables_2026-08-24/"
                              "guide_kd_efficiency.suppl_table.csv", index_col=0)
    miss = sorted(bank_guides - set(kdt.index.astype(str)))
    check("KA-5a the CD4 bank's targeting guides all appear in the authors' guide_kd_efficiency table",
          "0 missing", 0, dict(n_bank_guides=len(bank_guides), missing=len(miss), examples=miss[:5]),
          len(miss) == 0)

    out = dict(ka=ka, failures=FAIL, seconds=round(time.time() - t0, 1))
    L.Logger.dump(P, out, X.WORK / "ka296.json")
    P(f"\ndone {out['seconds']}s -- {len(FAIL)} failures")
    P.save()
    if any(f.startswith(("KA-1", "KA-2", "KA-3", "KA-4")) for f in FAIL):
        raise SystemExit(1)


if __name__ == "__main__":
    main()

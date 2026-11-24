"""170 — turn the sweep into the six-member table, the slope, and the board translation.

Scaling anchors are the board-fitted ones this repo already uses (094 §1 / 153's `debench`,
and 096 §(d)'s recovered `mse` map). `nmae` is lower-is-better.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

ANCHOR = {"pds": (0.50020, 0.95228), "fid": (0.51239, 0.80964),
          "reach": (0.07902, 0.96414), "jac": (0.03038, 0.39727),
          "nmae": (1.0013, 0.400)}
MSE_B, MSE_SPAN = 0.987862879, 0.950162566
BASE = 0.65


def scaled(member, raw):
    b, r = ANCHOR[member]
    return (raw - b) / (r - b)


def mse_scaled(raw):
    return min(max((MSE_B - raw) / MSE_SPAN, 0.0), 1.0)


def fit_slope(x, y):
    """OLS slope with its standard error, plus the endpoint slope."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = x.size
    sx, sy = x.mean(), y.mean()
    sxx = ((x - sx) ** 2).sum()
    b1 = ((x - sx) * (y - sy)).sum() / sxx
    b0 = sy - b1 * sx
    resid = y - (b0 + b1 * x)
    se = np.sqrt((resid ** 2).sum() / max(n - 2, 1) / sxx)
    r2 = 1.0 - (resid ** 2).sum() / max(((y - sy) ** 2).sum(), 1e-300)
    return b1, se, r2


def main() -> None:
    rows = json.loads((WORK / "h1sweep.json").read_text())
    rows.sort(key=lambda r: r["gamma"])
    g = np.array([r["gamma"] for r in rows])
    base = [r for r in rows if abs(r["gamma"] - BASE) < 1e-9][0]

    print("\n===== h1-test, six members, raw =====")
    hdr = (f"{'gamma':>6}{'mse':>10}{'pds':>10}{'nmae':>10}{'reach':>10}{'fid':>10}"
           f"{'jac':>10}{'A':>10}{'X':>11}{'B':>10}{'claim':>9}{'rho':>7}{'cos':>9}")
    print(hdr)
    for r in rows:
        print(f"{r['gamma']:>6.2f}{r['mse']:>10.5f}{r['pds']:>10.5f}{r['nmae']:>10.5f}"
              f"{r['reach']:>10.5f}{r['fid']:>10.5f}{r['jac']:>10.5f}{r['A']:>10.5f}"
              f"{r['X']:>11.6f}{r['B']:>10.5f}{r['claim']:>9.5f}{r['rho']:>7.4f}"
              f"{r['cos_agg']:>9.5f}")

    print(f"\n===== h1-test, delta vs gamma {BASE}, in OVERALL units (member/6) =====")
    print(f"{'gamma':>6}{'mse':>10}{'pds':>10}{'nmae':>10}{'reach':>10}{'fid':>10}"
          f"{'jac':>10}{'sum':>11}{'sum ex-mse':>12}")
    for r in rows:
        d = {}
        d["mse"] = (mse_scaled(r["mse"]) - mse_scaled(base["mse"])) / 6
        d["pds"] = (scaled("pds", r["pds"]) - scaled("pds", base["pds"])) / 6
        d["nmae"] = (scaled("nmae", r["nmae"]) - scaled("nmae", base["nmae"])) / 6
        for m in ("reach", "fid", "jac"):
            d[m] = (scaled(m, r[m]) - scaled(m, base[m])) / 6
        tot = sum(d.values())
        print(f"{r['gamma']:>6.2f}{d['mse']:>+10.5f}{d['pds']:>+10.5f}{d['nmae']:>+10.5f}"
              f"{d['reach']:>+10.5f}{d['fid']:>+10.5f}{d['jac']:>+10.5f}{tot:>+11.5f}"
              f"{tot - d['mse']:>+12.5f}")
        r["_d"] = d

    print("\n===== the slope a19 implies, against the slope measured here =====")
    for key, label in (("mse", "raw mse"), ("X", "X (alignment)"), ("A", "A (amplitude)"),
                       ("pds", "raw pds"), ("nmae", "raw nmae")):
        y = np.array([r[key] for r in rows])
        b1, se, r2 = fit_slope(g, y)
        print(f"  d({label:<16})/d(gamma) = {b1:>+10.5f} +- {se:.5f}   R2 {r2:.4f}   "
              f"[{y.min():.5f} .. {y.max():.5f}]")
    # relative slope of X, the quantity that transfers
    lx = np.log(np.array([r["X"] for r in rows]))
    b1, se, r2 = fit_slope(g, lx)
    print(f"  d(log X)/d(gamma)          = {b1:>+10.5f} +- {se:.5f}   R2 {r2:.4f}")
    print(f"  a19's board-implied d(raw mse)/d(gamma) on the CHALLENGE panel: +0.6075 "
          f"(one point, lam-confounded, extrapolated in the opposite direction)")

    # -------- re-seeded twin: the arm-to-arm noise floor on this surface
    twin_path = WORK / "h1sweep_s1.json"
    if twin_path.exists():
        twins = {r["gamma"]: r for r in json.loads(twin_path.read_text())}
        print("\n===== RE-SEEDED TWIN (seed 1 vs seed 0): the noise floor =====")
        print(f"{'gamma':>6}{'d mse':>10}{'d pds':>10}{'d nmae':>10}{'d reach':>10}"
              f"{'d fid':>10}{'d jac':>10}{'d overall':>11}")
        for gam, tw in sorted(twins.items()):
            b = [r for r in rows if abs(r["gamma"] - gam) < 1e-9][0]
            d = {"mse": (mse_scaled(tw["mse"]) - mse_scaled(b["mse"])) / 6,
                 "pds": (scaled("pds", tw["pds"]) - scaled("pds", b["pds"])) / 6,
                 "nmae": (scaled("nmae", tw["nmae"]) - scaled("nmae", b["nmae"])) / 6}
            for m in ("reach", "fid", "jac"):
                d[m] = (scaled(m, tw[m]) - scaled(m, b[m])) / 6
            print(f"{gam:>6.2f}{d['mse']:>+10.5f}{d['pds']:>+10.5f}{d['nmae']:>+10.5f}"
                  f"{d['reach']:>+10.5f}{d['fid']:>+10.5f}{d['jac']:>+10.5f}"
                  f"{sum(d.values()):>+11.5f}")
            print(f"       raw: mse {tw['mse'] - b['mse']:+.5f}  pds "
                  f"{tw['pds'] - b['pds']:+.6f}  nmae {tw['nmae'] - b['nmae']:+.5f}  "
                  f"reach {tw['reach'] - b['reach']:+.5f}  fid {tw['fid'] - b['fid']:+.5f}"
                  f"  jac {tw['jac'] - b['jac']:+.5f}")

    # -------- challenge translation
    for tag in ("minnorm", "none"):
        f = WORK / f"chal_{tag}.json"
        s = WORK / f"chal_{tag}_solve.json"
        if not (f.exists() and s.exists()):
            continue
        crows = json.loads(f.read_text())
        crows.sort(key=lambda r: r["gamma"])
        sol = json.loads(s.read_text())
        X065, sumD = sol["X_at_065"], sol["sumD"]
        xr = {r["gamma"]: r["X"] / base["X"] for r in rows}
        print(f"\n===== CHALLENGE translation, MASS_FIX={tag} "
              f"(anchor {sol['anchor']}, raw {sol['target_raw']:.6f}) =====")
        print(f"{'gamma':>6}{'A':>10}{'B':>10}{'claim':>9}{'X_h1/X0':>9}{'X_c':>10}"
              f"{'raw mse':>10}{'score':>9}{'d overall':>11}")
        b0 = None
        for r in crows:
            gam = r["gamma"]
            if gam not in xr:
                continue
            Xc = X065 * xr[gam]
            per = []
            for c in ("A", "B", "C"):
                p = r["per"][c]
                per.append((p["A"] + p["Kc"] - min(p["claim"], p["B"]) - 2 * Xc) / sumD + 1.0)
            raw = float(np.mean(per))
            sc = mse_scaled(raw)
            if abs(gam - BASE) < 1e-9:
                b0 = sc
            print(f"{gam:>6.2f}{r['A']:>10.6f}{r['B']:>10.6f}{r['claim']:>9.5f}"
                  f"{xr[gam]:>9.4f}{Xc:>10.6f}{raw:>10.6f}{sc:>9.5f}"
                  f"{(sc - (b0 if b0 is not None else sc)) / 6:>+11.5f}")
        ys = []
        for r in crows:
            gam = r["gamma"]
            if gam not in xr:
                continue
            Xc = X065 * xr[gam]
            ys.append(float(np.mean([
                (r["per"][c]["A"] + r["per"][c]["Kc"]
                 - min(r["per"][c]["claim"], r["per"][c]["B"]) - 2 * Xc) / sumD + 1.0
                for c in ("A", "B", "C")])))
        gg = np.array([r["gamma"] for r in crows if r["gamma"] in xr])
        b1, se, r2 = fit_slope(gg, np.array(ys))
        print(f"  d(raw mse)/d(gamma) on the CHALLENGE panel = {b1:+.5f} +- {se:.5f} "
              f"(R2 {r2:.4f})  vs a19's implied +0.6075  -> ratio {0.6075 / b1:.2f}x")


if __name__ == "__main__":
    main()

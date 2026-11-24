"""233 — the retention table, the attribution, and the composite, from `score233.json`.

RETENTION.  For every member, `new panel / today's panel`, taken as a ratio of CREDIT above a
measured null rather than of the raw value, because several members have a large panel-
independent floor (`pds` is 0.5 at chance, `fid` 0.5, `nmae` 1.0) and a raw ratio would
report retention of the floor.  The null is the `null` arm: the shipped stack with its
predicted rows DERANGED, measured separately on each panel and fold, so the floor is measured
and not assumed.

  credit(pds)   = pds  - null_pds          credit(cos_agg) = cos_agg - null_cos_agg
  credit(mse*)  = null_mse* - mse*         (mse* = 1 - cos^2, lower is better)
  credit(nmae)  = null_nmae - nmae         (lower is better)
  credit(fid)   = fid  - null_fid
  credit(reach) = reach - null_reach
  credit(jac)   = jac  - null_jac

`pds` is additionally reported at MATCHED DECOY COUNT (`pds_at150`, `pds_at250`), because the
number of evaluable targets differs between panels (269/290/294 today, 300 everywhere on the
realistic draw, 150/291/293 on the pessimistic one) and `pds` is a retrieval statistic whose
value depends on how many decoys there are.  Every `pds` retention ratio quoted in RESULTS.md
is the matched one.

COMPOSITE.  Each member's credit is converted to board units with AGENTS.md §4's spans
(`/span/6`), which gives an interpretable "how much of our overall score is panel-fitted"
number.  It is a PROJECTION, not a board prediction: PLOLO levels are not board levels, and
the local-to-board factor is not identified (§8.4).  Only the RATIO is claimed.

    uv run python retain233.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panels233 as P  # noqa: E402

SPANS = dict(pds=0.452, mse=0.950162566, nmae=0.5963, fid=0.30, jac=0.365, reach=0.895)
DIRECTION = ("pds", "mse", "nmae")
ORDERING = ("fid", "reach", "jac")
PANELS = ("abc", "real", "pess")


def credits(row: dict, null: dict) -> dict:
    return dict(
        cos_agg=row["cos_agg"] - null["cos_agg"],
        cos=row["cos"] - null["cos"],
        pds=row["pds"] - null["pds"],
        pds_at140=row.get("pds_at140", float("nan")) - null.get("pds_at140", float("nan")),
        pds_at150=row.get("pds_at150", float("nan")) - null.get("pds_at150", float("nan")),
        pds_at250=row.get("pds_at250", float("nan")) - null.get("pds_at250", float("nan")),
        mse=null["mse_star"] - row["mse_star"],
        nmae=null["nmae"] - row["nmae"],
        fid=row["fid"] - null["fid"],
        reach=row["reach"] - null["reach"],
        jac=row["jac"]["89"] - null["jac"]["89"],
        jac200=row["jac"]["200"] - null["jac"]["200"],
    )


def main() -> None:
    rows = json.loads((P.WORK / "score233.json").read_text())
    R = {(r["panel"], r["fold"], r["arm"]): r for r in rows}
    # `pds@140` is computed separately (`pdsmatch233.py`) because 140 is the only decoy count
    # below every fold's evaluable-target count, including the pessimistic k562 fold's 150.
    pm_path = P.WORK / "pdsmatch233.json"
    if pm_path.exists():
        for k, v in json.loads(pm_path.read_text()).items():
            panel, fold, arm = k.split("|")
            if (panel, fold, arm) in R:
                R[(panel, fold, arm)].update({kk: vv for kk, vv in v.items()
                                              if kk.startswith("pds_at")})
    folds = sorted({r["fold"] for r in rows})
    panels = [p for p in PANELS if any(r["panel"] == p for r in rows)]
    out: dict[str, object] = {}

    # ---------------------------------------------------------------- levels
    print("LEVELS — the shipped stack (`BASE`) and its deranged null, per panel per fold\n")
    hdr = (f"{'panel':6}{'fold':15}{'n':>5}{'cos_agg':>10}{'pds':>9}{'pds@140':>9}"
           f"{'mse*':>9}{'nmae':>9}{'fid':>8}{'reach':>8}{'jac89':>8}")
    print(hdr)
    for panel in panels:
        for fold in folds:
            for arm in ("BASE", "null"):
                r = R.get((panel, fold, arm))
                if not r:
                    continue
                print(f"{panel if arm == 'BASE' else '':6}{(fold + ('' if arm == 'BASE' else ' (null)')):15}"
                      f"{r['n_targets']:5d}{r['cos_agg']:10.5f}{r['pds']:9.5f}"
                      f"{r.get('pds_at140', float('nan')):9.5f}{r['mse_star']:9.5f}"
                      f"{r['nmae']:9.5f}{r['fid']:8.4f}{r['reach']:8.4f}{r['jac']['89']:8.4f}")
        print()

    # ---------------------------------------------------------------- retention
    print("\nRETENTION — credit above the deranged null, new panel / today's panel\n")
    cred = {}
    for panel in panels:
        for fold in folds:
            b, n = R.get((panel, fold, "BASE")), R.get((panel, fold, "null"))
            if b and n:
                cred[(panel, fold)] = credits(b, n)
    keys = ("cos_agg", "pds_at140", "pds_at250", "pds", "mse", "nmae", "fid", "reach",
            "jac", "jac200")
    for panel in [p for p in panels if p != "abc"]:
        print(f"--- {panel} vs abc")
        print(f"{'member':12}" + "".join(f"{f[:13]:>15}" for f in folds) + f"{'MEAN':>12}")
        tab = {}
        for k in keys:
            vals = []
            for f in folds:
                a, b = cred.get(("abc", f)), cred.get((panel, f))
                vals.append(b[k] / a[k] if (a and b and a[k] != 0) else float("nan"))
            tab[k] = vals + [float(np.nanmean(vals))]
            print(f"{k:12}" + "".join(f"{100 * v:14.1f}%" for v in tab[k][:-1])
                  + f"{100 * tab[k][-1]:11.1f}%")
        out[f"retention|{panel}"] = tab
        # channel summary
        dirn = float(np.nanmean([tab["pds_at140"][-1], tab["mse"][-1], tab["nmae"][-1]]))
        ordn = float(np.nanmean([tab["fid"][-1], tab["reach"][-1], tab["jac"][-1]]))
        out[f"channel|{panel}"] = dict(direction=dirn, ordering=ordn,
                                       cos_agg=tab["cos_agg"][-1])
        print(f"\n  DIRECTION (pds@matched, mse*, nmae): {100 * dirn:.1f}%"
              f"     ORDERING (fid, reach, jac): {100 * ordn:.1f}%")
        print(f"  aggregate cosine:                    {100 * tab['cos_agg'][-1]:.1f}%\n")

    # ---------------------------------------------------------------- attribution
    print("\nATTRIBUTION — each component's contribution, per panel (BASE - ablation)\n")
    abls = sorted({r["arm"] for r in rows} - {"BASE", "null"})
    for panel in panels:
        print(f"--- {panel}")
        print(f"{'ablation':22}{'d cos_agg':>12}{'rel %':>9}{'d pds':>11}"
              f"{'d nmae':>10}{'d fid':>9}{'d reach':>9}{'d jac89':>9}   folds")
        for arm in abls:
            d = {k: [] for k in ("cos_agg", "rel", "pds", "nmae", "fid", "reach", "jac")}
            nf = 0
            for fold in folds:
                b, a = R.get((panel, fold, "BASE")), R.get((panel, fold, arm))
                if not (b and a):
                    continue
                nf += 1
                d["cos_agg"].append(b["cos_agg"] - a["cos_agg"])
                d["rel"].append(100 * (b["cos_agg"] / a["cos_agg"] - 1)
                                if a["cos_agg"] else float("nan"))
                d["pds"].append(b["pds"] - a["pds"])
                d["nmae"].append(a["nmae"] - b["nmae"])
                d["fid"].append(b["fid"] - a["fid"])
                d["reach"].append(b["reach"] - a["reach"])
                d["jac"].append(b["jac"]["89"] - a["jac"]["89"])
            if not nf:
                continue
            m = {k: float(np.nanmean(v)) for k, v in d.items()}
            out[f"contribution|{panel}|{arm}"] = m
            print(f"{arm:22}{m['cos_agg']:12.6f}{m['rel']:8.2f}%{m['pds']:11.6f}"
                  f"{m['nmae']:10.5f}{m['fid']:9.5f}{m['reach']:9.5f}{m['jac']:9.5f}   {nf}")
        print()

    # ---------------------------------------------------------------- composite
    print("\nCOMPOSITE PROJECTION — credit in board units (credit/span/6), and its ratio\n")
    print(f"{'panel':8}{'o_pds':>10}{'o_mse':>10}{'o_nmae':>10}{'o_fid':>10}"
          f"{'o_reach':>10}{'o_jac':>10}{'DIR':>10}{'ORD':>10}{'TOTAL':>10}")
    comp = {}
    for panel in panels:
        v = {}
        for f_ in folds:
            c = cred.get((panel, f_))
            if not c:
                continue
            for k, span in SPANS.items():
                kk = "pds_at140" if k == "pds" else k
                v.setdefault(k, []).append(c[kk] / span / 6.0)
        if not v:
            continue
        m = {k: float(np.nanmean(x)) for k, x in v.items()}
        m["DIR"] = m["pds"] + m["mse"] + m["nmae"]
        m["ORD"] = m["fid"] + m["reach"] + m["jac"]
        m["TOTAL"] = m["DIR"] + m["ORD"]
        comp[panel] = m
        print(f"{panel:8}" + "".join(f"{m[k]:10.5f}" for k in
                                     ("pds", "mse", "nmae", "fid", "reach", "jac",
                                      "DIR", "ORD", "TOTAL")))
    for panel in [p for p in panels if p != "abc"]:
        if panel in comp and "abc" in comp:
            print(f"  {panel} / abc:  DIR {100 * comp[panel]['DIR'] / comp['abc']['DIR']:.1f}%"
                  f"   ORD {100 * comp[panel]['ORD'] / comp['abc']['ORD']:.1f}%"
                  f"   TOTAL {100 * comp[panel]['TOTAL'] / comp['abc']['TOTAL']:.1f}%")
    out["composite"] = comp
    (P.WORK / "retain233.json").write_text(json.dumps(out, indent=2))
    print("\n-> work/retain233.json")


if __name__ == "__main__":
    main()

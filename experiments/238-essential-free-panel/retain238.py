"""238 — the RETENTION table on the essential-free panels, from `six238.json`.

`experiments/233-october-retention/retain233.py`'s conventions, unchanged, so the essential-free
rows are directly comparable to 233's published 65% / 0.137 projection:

  credit(pds)   = pds  - null_pds          credit(cos_agg) = cos_agg - null_cos_agg
  credit(mse*)  = null_mse* - mse*         (mse* = 1 - cos^2, lower is better)
  credit(nmae)  = null_nmae - nmae         (lower is better)
  credit(fid)   = fid  - null_fid          credit(reach)   = reach - null_reach
  credit(jac)   = jac  - null_jac

The null is the shipped stack with its predicted rows DERANGED, measured separately on every
panel and fold, so each member's floor (`pds` 0.5, `fid` 0.5, `nmae` 1.0) is measured and not
assumed.  `pds` is quoted at MATCHED DECOY COUNT (`pds_at140`; 140 is below every fold's
evaluable-target count on all five panels) because `pds` is a retrieval statistic whose value
depends on how many decoys there are.

§3.3's board projection is reproduced too -- our own live `h1union` row scaled by each member's
retention ratio.  **It is a projection, not a board prediction** (233 §3.3): it assumes each
board member's scaled credit scales like the local credit above the deranged null, which is a
different normalisation.  `fid` is held fixed because our board `fid` is NEGATIVE and a ratio
applied to a negative credit is meaningless.

    uv run python retain238.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panels238 as P  # noqa: E402

SPANS = dict(pds=0.452, mse=0.950162566, nmae=0.5963, fid=0.30, jac=0.365, reach=0.895)
PANELS = ("abc", "real", "pess", "efr", "efp")
# our live board row: `h1union`, entry iUziqE1gewfcPKzYUDaB, overall 0.21183096
BOARD = dict(pds=0.745374, mse=0.179569, nmae=0.110742, fid=-0.016476,
             reach=0.222939, jac=0.028837)
BOARD_OVERALL = 0.21183096


def credits(row: dict, null: dict) -> dict:
    return dict(
        cos_agg=row["cos_agg"] - null["cos_agg"],
        cos=row["cos"] - null["cos"],
        pds=row["pds"] - null["pds"],
        pds_at140=row.get("pds_at140", float("nan")) - null.get("pds_at140", float("nan")),
        pds_at250=row.get("pds_at250", float("nan")) - null.get("pds_at250", float("nan")),
        mse=null["mse_star"] - row["mse_star"],
        nmae=null["nmae"] - row["nmae"],
        fid=row["fid"] - null["fid"],
        reach=row["reach"] - null["reach"],
        jac=row["jac"]["89"] - null["jac"]["89"],
        jac200=row["jac"]["200"] - null["jac"]["200"],
    )


def main() -> None:
    rows = json.loads((P.WORK / "six238.json").read_text())
    R = {(r["panel"], r["fold"], r["arm"]): r for r in rows if not r.get("vacuous")}
    folds = [f for f in P.CANCER_FOLDS if any(r["fold"] == f for r in rows)]
    panels = [p for p in PANELS if (p, folds[0], "BASE") in R or
              any((p, f, "BASE") in R for f in folds)]
    out: dict[str, object] = {}

    print("LEVELS — the shipped stack (`BASE`) and its deranged null, per panel per fold\n")
    print(f"{'panel':6}{'fold':22}{'n':>5}{'cos_agg':>10}{'pds':>9}{'pds@140':>9}"
          f"{'mse*':>9}{'nmae':>9}{'fid':>8}{'reach':>8}{'jac89':>8}")
    for panel in panels:
        for fold in folds:
            for arm in ("BASE", "null"):
                r = R.get((panel, fold, arm))
                if not r:
                    continue
                lab = fold + ("" if arm == "BASE" else " (null)")
                print(f"{panel if arm == 'BASE' else '':6}{lab:22}{r['n_targets']:5d}"
                      f"{r['cos_agg']:10.5f}{r['pds']:9.5f}"
                      f"{r.get('pds_at140', float('nan')):9.5f}{r['mse_star']:9.5f}"
                      f"{r['nmae']:9.5f}{r['fid']:8.4f}{r['reach']:8.4f}"
                      f"{r['jac']['89']:8.4f}")
        print()

    cred = {}
    for panel in panels:
        for fold in folds:
            b, n = R.get((panel, fold, "BASE")), R.get((panel, fold, "null"))
            if b and n:
                cred[(panel, fold)] = credits(b, n)
    out["credits"] = {f"{k[0]}|{k[1]}": v for k, v in cred.items()}

    print("\nRETENTION — credit above the deranged null, new panel / today's panel\n")
    keys = ("cos_agg", "cos", "pds_at140", "pds_at250", "pds", "mse", "nmae", "fid",
            "reach", "jac", "jac200")
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
        dirn = float(np.nanmean([tab["pds_at140"][-1], tab["mse"][-1], tab["nmae"][-1]]))
        ordn = float(np.nanmean([tab["fid"][-1], tab["reach"][-1], tab["jac"][-1]]))
        out[f"channel|{panel}"] = dict(direction=dirn, ordering=ordn,
                                       cos_agg=tab["cos_agg"][-1])
        print(f"\n  DIRECTION (pds@140, mse*, nmae): {100 * dirn:.1f}%"
              f"     ORDERING (fid, reach, jac): {100 * ordn:.1f}%")
        print(f"  aggregate cosine:                {100 * tab['cos_agg'][-1]:.1f}%\n")

    print("\nCOMPOSITE — credit in board units (credit/span/6), and its ratio to abc\n")
    print(f"{'panel':8}{'o_pds':>10}{'o_mse':>10}{'o_nmae':>10}{'o_fid':>10}"
          f"{'o_reach':>10}{'o_jac':>10}{'DIR':>10}{'ORD':>10}{'TOTAL':>10}")
    comp = {}
    for panel in panels:
        v: dict[str, list] = {}
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

    # ------------------------------------------------------ 233 §3.3's board-unit projection
    print("\nBOARD-UNIT PROJECTION (233 §3.3's arithmetic; a PROJECTION, not a prediction)\n")
    print(f"{'':28}" + "".join(f"{k:>10}" for k in
                               ("pds", "mse", "nmae", "fid", "reach", "jac"))
          + f"{'OVERALL':>10}{'of today':>10}")
    print(f"{'today (h1union board)':28}"
          + "".join(f"{BOARD[k]:10.4f}" for k in
                    ("pds", "mse", "nmae", "fid", "reach", "jac"))
          + f"{BOARD_OVERALL:10.5f}{100.0:9.1f}%")
    proj = {}
    for panel in [p for p in panels if p != "abc"]:
        tab = out.get(f"retention|{panel}")
        if not tab:
            continue
        ratio = dict(pds=tab["pds_at140"][-1], mse=tab["mse"][-1], nmae=tab["nmae"][-1],
                     fid=1.0, reach=tab["reach"][-1], jac=tab["jac"][-1])
        v = {k: BOARD[k] * ratio[k] for k in BOARD}
        overall = float(np.mean(list(v.values())))
        proj[panel] = dict(members=v, overall=overall, of_today=overall / BOARD_OVERALL)
        print(f"{panel:28}" + "".join(f"{v[k]:10.4f}" for k in
                                      ("pds", "mse", "nmae", "fid", "reach", "jac"))
              + f"{overall:10.5f}{100 * overall / BOARD_OVERALL:9.1f}%")
    out["board_projection"] = proj
    print("\n(`fid` is held at its board value: our board `fid` is negative and a ratio "
          "applied to a negative credit is meaningless — 233 §3.3's own rule.)")

    (P.WORK / "retain238.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> {P.WORK / 'retain238.json'}")


if __name__ == "__main__":
    main()

"""534 Part B -- what each priced panel looks like to the instrument (light; reads the comp534/full534 jsons):
fold coverage, H1-union / CD4 rows, carriers, the percentile cap per fold and on the full panel, PLOLO retention ret_f,
level factors, operator identity and the axis check each operator passed.

    uv run python census534.py [panels...]        -> work/partB/census534.json (+ printed table)
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import json  # noqa: E402

import price534 as P5  # noqa: E402

HERE = P5.HERE


def main(panels):
    out = {}
    for p in panels:
        full = json.loads((HERE / p / "full534.json").read_text())
        lv = json.loads((HERE / p / "level534.json").read_text())
        P = P5.Panel(p)
        rec = dict(full=dict(op=full["op"], op_sha256=full["op_sha256"], axis_checked=full["op_axis_checked"],
                             coverage=full["coverage"], carriers=full["carriers"], pct_hi=full["keys"]["pct_hi"],
                             claims=full["claims"], rho_target=lv["rho_target"],
                             factors={n: lv["factors"][n] for n in ("h1c1|K3", "h2c1|K3+sh1.75", "h2c1|V1pct+sh1.75")},
                             peak_rss_gib=full["peak_rss_gib"]), folds={})
        for f in P5.CANCER:
            c = json.loads((HERE / p / f"comp534-{f}.json").read_text())
            rec["folds"][f] = dict(n_targets=c["n_targets"], op=c["op"], op_sha256=c["op_sha256"], axis_checked=c["op_axis_checked"],
                                   coverage=c["coverage"], carried_lodo=c["keys"]["carried_lodo"], pct_hi=c["keys"]["pct_hi"],
                                   above3=c["keys"]["n_above_cap3_under_pct"], ret=P.ret[f], peak_rss_gib=c["peak_rss_gib"])
        out[p] = rec
        print(f"{p}: full op {full['op_sha256'][:12]} ({full['op_axis_checked']}); donors {full['coverage']['donors']}; "
              f"H1 rows {full['coverage']['h1_union_rows']}; CD4 complete {full['coverage']['cd4_complete_rows']}; carriers "
              f"{full['carriers']}; pct hi {full['keys']['pct_hi']:.2f}")
        for f, r in rec["folds"].items():
            print(f"   {f:14s} n {r['n_targets']:3d}  op {r['op_sha256'][:12]} ({r['axis_checked']})  H1 {r['coverage']['h1_union_rows']:3d}  "
                  f"CD4 {r['coverage']['cd4_complete_rows']:3d}  carried {r['carried_lodo']:3d}  pct hi {r['pct_hi']:6.2f} ({r['above3']} > 3)  "
                  f"ret {r['ret']:.3f}  peak {r['peak_rss_gib']:.2f} GiB")
    (HERE / "census534.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(sys.argv[1:] or ["efp", "S1", "S2", "S3", "S4"])

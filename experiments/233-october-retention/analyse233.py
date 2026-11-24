"""233 — the two remaining tables: the coverage-thinning decomposition and the October sources."""
import json, sys
import numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import panels233 as P

S = {(r["panel"], r["fold"], r["arm"]): r for r in json.loads((P.WORK/"score233.json").read_text())}
for k, v in json.loads((P.WORK/"pdsmatch233.json").read_text()).items():
    a, b, c = k.split("|")
    if (a, b, c) in S: S[(a, b, c)].update({x: y for x, y in v.items() if x.startswith("pds_at")})
T = json.loads((P.WORK/"thin233.json").read_text())

print("COVERAGE THINNING — today's panel, each component's coverage cut to the pessimistic rate")
print(f"{'fold':15}{'arm':12}{'kept':>18}{'cos_agg':>10}{'ret cos':>9}{'pds':>9}{'ret pds':>9}")
out = {}
for fold in sorted({r['fold'] for r in T}):
    b, n = S[("abc", fold, "BASE")], S[("abc", fold, "null")]
    pb, pn = S[("pess", fold, "BASE")], S[("pess", fold, "null")]
    cb, cp = b["cos_agg"]-n["cos_agg"], pb["cos_agg"]-pn["cos_agg"]
    db, dp = b["pds"]-n["pds"], pb["pds"]-pn["pds"]
    print(f"{fold:15}{'abc BASE':12}{'':>18}{b['cos_agg']:10.5f}{100.0:8.1f}%{b['pds']:9.5f}{100.0:8.1f}%")
    for arm in ("thin_k562", "thin_cd4", "thin_h1", "thin_all"):
        rs = [r for r in T if r["fold"] == fold and r["arm"] == arm]
        if not rs: continue
        c = float(np.mean([r["cos_agg"] for r in rs])); p = float(np.mean([r["pds"] for r in rs]))
        sd = float(np.std([r["cos_agg"] for r in rs]))
        kept = f"{rs[0]['source_kept']}/{rs[0]['cd4_kept']}/{rs[0]['h1_kept']}"
        out[f"{fold}|{arm}"] = dict(cos_agg=c, cos_sd=sd, pds=p,
                                    ret_cos=100*(c-n["cos_agg"])/cb, ret_pds=100*(p-n["pds"])/db)
        print(f"{'':15}{arm:12}{kept:>18}{c:10.5f}{100*(c-n['cos_agg'])/cb:8.1f}%"
              f"{p:9.5f}{100*(p-n['pds'])/db:8.1f}%")
    print(f"{'':15}{'pess BASE':12}{'':>18}{pb['cos_agg']:10.5f}{100*cp/cb:8.1f}%"
          f"{pb['pds']:9.5f}{100*dp/db:8.1f}%")
    out[f"{fold}|pess"] = dict(ret_cos=100*cp/cb, ret_pds=100*dp/db)
    print()

print("\nOCTOBER-ONLY SOURCES — un-whitened post-whitening residual on the shipped stack")
rows = json.loads((P.WORK/"oct233.json").read_text())
base = {(r["panel"], r["fold"]): r for r in rows if r["arm"] == "BASE"}
agg = {}
for r in rows:
    if r["arm"] == "BASE": continue
    name, _, dose = r["arm"].partition("@")
    b = base[(r["panel"], r["fold"])]
    agg.setdefault((r["panel"], name, dose), []).append(
        (100*(r["cos_agg"]/b["cos_agg"]-1), r["pds"]-b["pds"], r["carrying_rows"],
         r["per_unit_carrying"], r["per_unit_all_rows"]))
print(f"{'panel':6}{'source':22}{'dose':>6}{'folds':>6}{'carry':>7}"
      f"{'rel cos %':>11}{'d pds':>10}{'per_unit':>10}{'all-rows':>10}")
for (panel, name, dose), v in sorted(agg.items()):
    if panel == "abc": continue
    a = np.mean([x[0] for x in v]); d = np.mean([x[1] for x in v])
    print(f"{panel:6}{name:22}{dose:>6}{len(v):6d}{int(np.mean([x[2] for x in v])):7d}"
          f"{a:10.2f}%{d:10.5f}{np.mean([x[3] for x in v]):10.4f}"
          f"{np.mean([x[4] for x in v]):10.4f}")
(P.WORK/"analyse233.json").write_text(json.dumps(
    {"thinning": out,
     "october": {f"{p}|{n}|{d}": dict(rel_cos_pct=float(np.mean([x[0] for x in v])),
                                      d_pds=float(np.mean([x[1] for x in v])),
                                      folds=len(v),
                                      carrying=float(np.mean([x[2] for x in v])),
                                      per_unit_carrying=float(np.mean([x[3] for x in v])),
                                      per_unit_all_rows=float(np.mean([x[4] for x in v])))
                 for (p, n, d), v in agg.items()}}, indent=2))
print("-> work/analyse233.json")

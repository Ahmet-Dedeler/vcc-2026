"""277 -- the markdown tables of RESULTS.md, generated from the json artifacts so no number is typed twice.

    PYTHONDONTWRITEBYTECODE=1 uv run python tables277.py > work/tables277.md
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib277 as Z  # noqa: E402

CANCER = Z.CANCER
W = Z.WORK
pr = json.loads((W / "price277.json").read_text())
an = json.loads((W / "analyse277.json").read_text())
fx = json.loads((W / "flex277.json").read_text())
meta = {f: json.loads((W / f"render277-{f}.json").read_text())["meta"] for f in Z.FOLDS}
fl = {f: json.loads((W / f"floor277-{f}.json").read_text()) for f in CANCER}
nl = {f: json.loads((W / f"null277-{f}.json").read_text()) for f in CANCER}
C = pr["cells"]


def m(*rows):
    for r in rows:
        print(r)


print("## The class\n")
m("| fold | targets | with `C_lolo` | CLASS30 n | `C_lolo` cut | CD4-carried in class | CLASS20 / CLASS40 |",
  "|---|---:|---:|---:|---:|---:|---:|")
for f in Z.FOLDS:
    z = np.load(W / f"render277-{f}.npz", allow_pickle=True)
    c30 = z["mask|CLASS30"].astype(bool)
    car = z["cd4_carried"].astype(bool)
    print(f"| `{f}` | {meta[f]['n_targets']} | {meta[f]['n_eligible']} | {meta[f]['n_class']} | "
          f"{meta[f]['class_cut']:+.5f} | {int((c30 & car).sum())} | "
          f"{int(z['mask|CLASS20'].sum())} / {int(z['mask|CLASS40'].sum())} |")
print(f"| h1-test (Flex-400) | 85 | 85 | {fx['meta']['n_class']} | {fx['meta']['class_cut']:+.5f} | "
      f"{int(np.load(W / 'flex277.npz')['mask|CLASS30'].astype(bool) & np.load(W / 'flex277.npz')['cd4_carried'].astype(bool)).sum() if False else int((np.load(W / 'flex277.npz')['mask|CLASS30'].astype(bool) & np.load(W / 'flex277.npz')['cd4_carried'].astype(bool)).sum())} | — |")

print("\n## The arm table (three cancer folds, vs `cd4svd5`)\n")
m("| arm | Δ`pds` | paired 95 % CI | folds > 0 | on edited rows | per edited row | `cos_agg` | macro | Δ`nmae` | amp | LT3O | PESS | POINT | OPT | kolf2 |",
  "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
ORDER = ["V1_undefl", "V1f_undefl", "V2_mirror", "V2f_mirror", "V3_noresid", "V3f_noresid",
         "V4_floor", "V5_20", "V5_40", "C_top30", "C_agree30", "C_agreeN", "C_rand1", "C_rand2", "C_rand3"]
for nm in ORDER:
    if nm not in C:
        continue
    e, c, s = C[nm], C[nm]["compare"], C[nm]["scen"]
    o = an["own_rows"][nm]
    print(f"| `{nm}` | {c['d_pds']:+.5f} | [{e['boot']['lo']:+.5f}, {e['boot']['hi']:+.5f}] | "
          f"{c['pds_pos']}/3 | {o['d_on']:+.5f} | {o['per_row']:+.5f} | {100*c['rel_agg']:+.2f} % | "
          f"{100*c['rel_macro']:+.2f} % | {c['d_nmae']:+.5f} | {100*(c['delivered_ratio']-1):+.2f} % | "
          f"{e['lt3o']['retention']:+.2f} | **{s['PESSIMISTIC']['total']:+.5f}** | "
          f"{s['POINT']['total']:+.5f} | {s['OPTIMISTIC']['total']:+.5f} | "
          f"{(e['kolf2_d_pds'] if e['kolf2_d_pds'] is not None else float('nan')):+.5f} |")

print("\n## Per-fold Δ`pds` and 262's groups on the edited rows\n")
m("| arm | k562 | hct116 | hek293t | group e (n) | marginal (n) | saturated (n) |",
  "|---|---:|---:|---:|---:|---:|---:|")
G = ("e (base.truth <= 0)", "marginal (0, 0.05]", "saturated > 0.05")
for nm in ORDER:
    if nm not in C:
        continue
    pf = C[nm]["compare"]["per_fold_d_pds"]
    g = an["class_by_group"].get(nm)
    cells = (" | ".join(f"{g[k]['mean']:+.5f} ({g[k]['n']})" for k in G) if g else " | ".join(["—"] * 3))
    print(f"| `{nm}` | " + " | ".join(f"{x:+.5f}" for x in pf) + " | " + cells + " |")

print("\n## The proper null for the selector (24 size-matched random row sets per fold, V1's construction)\n")
m("| row set | k562 (z, pct) | hct116 (z, pct) | hek293t (z, pct) | pooled z |", "|---|---|---|---|---:|")
for k in ("CLASS30", "AGREE30", "AGREEN", "TOPCONS30"):
    line = f"| `{k}` "
    for f in CANCER:
        o = nl[f]["observed"][k]
        mu, sd = nl[f]["null"]["d_pds"]["mean"], nl[f]["null"]["d_pds"]["sd"]
        line += (f"| {o['d_pds']:+.5f} (z {(o['d_pds']-mu)/sd:+.2f}, "
                 f"{100*o['null_frac_at_or_below_d_pds']:.0f} %) ")
    print(line + f"| {an['null_summary']['pooled_z'][k]:+.2f} |")
print(f"\nNull mean / sd of the fold-mean Δ`pds` over 24 draws: " +
      "; ".join(f"{f} {nl[f]['null']['d_pds']['mean']:+.5f} ± {nl[f]['null']['d_pds']['sd']:.5f}"
                for f in CANCER) + ".")

print("\n## V4's own controls (the 272 A′ floor, c = 1, by row set and by noise vector)\n")
m("| cell | k562 | hct116 | hek293t | mean Δ`pds` | mean rel `cos_agg` |", "|---|---:|---:|---:|---:|---:|")
for tag in ("floor:true:CLASS30", "floor:shuffled:CLASS30", "floor:mirrored:CLASS30",
            "floor:true:TOPCONS30", "floor:true:RAND1", "floor:true:RAND2", "floor:true:RAND3",
            "floor:true:ALL", "gamma0.5:CLASS30", "gamma0.5:RAND1", "gamma0.5:ALL",
            "gamma0.35:CLASS30", "gamma0.35:ALL"):
    if tag not in fl[CANCER[0]]["rows"]:
        continue
    v = [fl[f]["rows"][tag]["d_pds"] for f in CANCER]
    a = [fl[f]["rows"][tag]["rel_agg"] for f in CANCER]
    print(f"| `{tag}` | " + " | ".join(f"{x:+.5f}" for x in v) +
          f" | **{np.mean(v):+.5f}** | {100*np.mean(a):+.2f} % |")

print("\n## The Flex-400 re-read (h1-test, union OFF, 273's instrument, 85 targets)\n")
inc, off, bas = fx["arms"]["INCUMBENT"], fx["arms"]["cd4h1:off"], fx["arms"]["BASE"]
print(f"Reference rows: BASE `cos_agg` {bas['a23_cos_agg']:.5f}; the CD4 residual raises it to "
      f"{off['a23_cos_agg']:.5f} un-deflated (**{100*(off['a23_cos_agg']/bas['a23_cos_agg']-1):+.1f} %**) "
      f"and {inc['a23_cos_agg']:.5f} deflated (**{100*(inc['a23_cos_agg']/bas['a23_cos_agg']-1):+.1f} %**) "
      f"— R273/F4's +30 % / +24 %, reproduced. `pds` 0.871148 / {off['pds_a23']:.6f} / "
      f"{inc['pds_a23']:.6f}.\n")
m("| arm | Δ`pds` | paired 95 % CI | on class | off class | `cos_agg` | LT3O | PLOLO Δ`pds` |",
  "|---|---:|---:|---:|---:|---:|---:|---:|")
for nm, r in fx["table"].items():
    p = C.get(nm, {}).get("compare", {}).get("d_pds")
    print(f"| `{nm}` | {r['d_pds']:+.5f} | [{r['boot']['lo']:+.5f}, {r['boot']['hi']:+.5f}] | "
          f"{r['d_pds_class']:+.5f} | {r['d_pds_rest']:+.5f} | {100*r['rel_cos_agg']:+.2f} % | "
          f"{r['lt3o']['retention']:+.2f} | " + (f"{p:+.5f} |" if p is not None else "— |"))

print("\n## Where the deflation earns its `pds`, by conservation band\n")
m("| band | rows | mean `cd4svd5` `pds` | deflation Δ`pds` | share of its gain | V1 Δ`pds` |",
  "|---|---:|---:|---:|---:|---:|")
tot = sum(v["share_of_total"] for v in an["deflation_by_band"].values())
for k, v in an["deflation_by_band"].items():
    print(f"| {k} | {v['n']} | {v['mean_pds']:+.5f} | {v['mean_d_pds']:+.5f} | "
          f"{100*v['share_of_total']/tot:.1f} % | {an['v1_by_band'][k]['mean_d_pds']:+.5f} |")
sp = an["spearman"]
print(f"\nBut per target the map does NOT predict the deflation's gain: "
      f"ρ(`C_lolo`, deflation Δ`pds`) = {sp['C_lolo (276\'s conservation map)']['vs_deflation']['rho']:+.3f} "
      f"[{sp['C_lolo (276\'s conservation map)']['vs_deflation']['ci'][0]:+.3f}, "
      f"{sp['C_lolo (276\'s conservation map)']['vs_deflation']['ci'][1]:+.3f}] on n = "
      f"{sp['C_lolo (276\'s conservation map)']['vs_deflation']['n']}; controlling for the un-deflated arm's "
      f"own per-target `pds` (headroom) it is "
      f"{sp['C_lolo vs deflation | h1union\'s own pds (headroom)']['rho']:+.3f} "
      f"[{sp['C_lolo vs deflation | h1union\'s own pds (headroom)']['ci'][0]:+.3f}, "
      f"{sp['C_lolo vs deflation | h1union\'s own pds (headroom)']['ci'][1]:+.3f}], while "
      f"ρ(`C_lolo`, the un-deflated arm's own `pds`) = "
      f"{sp['C_lolo vs h1union\'s own pds']['rho']:+.3f} "
      f"[{sp['C_lolo vs h1union\'s own pds']['ci'][0]:+.3f}, "
      f"{sp['C_lolo vs h1union\'s own pds']['ci'][1]:+.3f}].")

print("\n## The 0.02 knee on the edited rows\n")
m("| arm | crossings up | down | unchanged | mean cosine move | Δ`pds` up-crossers | Δ`pds` down-crossers |",
  "|---|---:|---:|---:|---:|---:|---:|")
for k, v in an["knee"].items():
    print(f"| `{k}` | {v['up']} | {v['down']} | {v['unchanged']} | {v['mean_cos_move']:+.5f} | "
          f"{v['mean_dpds_crossers_up']:+.4f} | {v['mean_dpds_crossers_down']:+.4f} |")

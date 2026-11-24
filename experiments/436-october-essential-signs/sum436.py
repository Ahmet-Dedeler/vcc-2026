"""436 summary: the three-seed tables per variant (from ka436.py's work/tab436-<variant>.json), the per-block capacity
story, the per-covered-target gain split by K562-essential-list membership, and the October scaling
(work/size436.json coverage x the per-covered-target gain, in overall units via 435's board slopes).
Run ka436.py <variant> first for every variant. Writes work/sum436.json and prints the RESULTS.md tables."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
W = HERE / "work/h1"; BASE = "mde0.1cap0.1min0.3coarse"
MEM = ["pds", "mse", "nmae", "fid", "reach", "jac"]
PT = ["pds", "nmae", "fid", "reach", "jac"]
SLOPE = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532, "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}
VARS = [v for v in ("dsall", "dsallcp", "dsesscp") if (HERE / f"work/tab436-{v}.json").exists()]
ana = json.loads((HERE / "work/ana436.json").read_text())
offlist = set(ana["covered_gated_not_on_k562ess_list"])
out = {"variants": {}}


def ms(a):
    a = np.asarray(a, dtype=float)
    return float(a.mean()), (float(a.std(ddof=1)) if a.size > 1 else 0.0)


for v in VARS:
    t = json.loads((HERE / f"work/tab436-{v}.json").read_text())
    seeds = sorted(t["seeds"], key=int)
    rec = {"seeds": [int(s) for s in seeds]}
    for part, ks in (("all", MEM), ("covered", PT), ("uncovered", PT)):
        rec[part] = {k: ms([t["seeds"][s][part][k] for s in seeds]) for k in ks}
    rec["overall_H1"] = ms([t["seeds"][s]["overall_all"] for s in seeds])
    # per covered target, overall units on a 300-target panel: sum_k slope_k * d_k / 6 / 300 (covered-target means)
    per_cov = {k: rec["covered"][k][0] for k in PT}
    rec["per_covered_target_overall_all_members"] = sum(SLOPE[k] * per_cov[k] for k in PT) / 6 / 300
    rec["per_covered_target_overall_nmae_only"] = SLOPE["nmae"] * per_cov["nmae"] / 6 / 300
    rec["leak"] = {k: ms([t["seeds"][s]["leak"][k] for s in seeds]) for k in t["seeds"][seeds[0]]["leak"]}
    rec["uncovered_bit_identical"] = [f"{t['seeds'][s]['uncovered_bit_identical']}/{t['seeds'][s]['uncovered_n']}" for s in seeds]
    # per covered target nmae change (mean over seeds), split on / off the K562-essential list
    pt = {}
    for s in seeds:
        po = pl.read_csv(W / f"per-target-base_s{s}_bps_em436_{BASE}.csv").filter(pl.col("metric") == "de_wilcoxon_lfc_nmae")
        pn = pl.read_csv(W / f"per-target-base_s{s}_bps_em436_{BASE}_{v}.csv").filter(pl.col("metric") == "de_wilcoxon_lfc_nmae")
        j = po.join(pn, on="perturbation", suffix="_on")
        for tg_, a_, b_ in zip(j["perturbation"], j["value"], j["value_on"]):
            if b_ != a_:
                pt.setdefault(tg_, []).append(b_ - a_)
    pt = {k: float(np.mean(x)) for k, x in pt.items()}
    rec["per_target_nmae"] = pt
    on_l = [x for k, x in pt.items() if k not in offlist]; off_l = [x for k, x in pt.items() if k in offlist]
    rec["covered_nmae_on_list"] = (float(np.mean(on_l)), len(on_l)) if on_l else None
    rec["covered_nmae_off_list"] = (float(np.mean(off_l)), len(off_l)) if off_l else None
    rec["covered_nmae_worse_targets"] = sorted([k for k, x in pt.items() if x > 0])
    # per-block record, seed 0
    rec["blocks_s0"] = [b for b in t["seeds"][seeds[0]]["blocks"] if b["covered"]]
    out["variants"][v] = rec

print("## three-seed tables (raw deltas vs the same-seed option-off run; mean ± sd)")
for v, r in out["variants"].items():
    print(f"\n### {v} (seeds {r['seeds']}): H1 overall {r['overall_H1'][0]:+.5f} ± {r['overall_H1'][1]:.5f}")
    print("| part | " + " | ".join(MEM) + " |"); print("|---|" + "---:|" * len(MEM))
    for part in ("all", "covered", "uncovered"):
        print(f"| {part} | " + " | ".join(f"{r[part][k][0]:+.5f} ± {r[part][k][1]:.5f}" if k in r[part] else "·" for k in MEM) + " |")
    lk = r["leak"]
    print(f"leak: unsuppressible/covered block off {lk['unsupp_mean_off'][0]:.1f} on {lk['unsupp_mean_on'][0]:.1f}; blocks > 50 on "
          f"{lk['blocks_unsupp_gt50_on'][0]:.1f} (off {lk['blocks_unsupp_gt50_off'][0]:.0f}); extra calls sum {lk['extra_calls_sum'][0]:.0f}, "
          f"max/block {lk['extra_calls_max'][0]:.0f}; mde_lambda<1 on {lk['mde_lambda_lt1_on'][0]:.0f} (off {lk['mde_lambda_lt1_off'][0]:.0f})")
    print(f"uncovered bit-identical {r['uncovered_bit_identical']}; covered nmae on-list {r['covered_nmae_on_list']}, off-list "
          f"{r['covered_nmae_off_list']}; worse: {r['covered_nmae_worse_targets']}")
    print(f"per covered target (300-target panel): all members {r['per_covered_target_overall_all_members']:+.3e}, "
          f"nmae only {r['per_covered_target_overall_nmae_only']:+.3e} overall")

# October scaling
size = json.loads((HERE / "work/size436.json").read_text())
scen = {"Arc criterion as observed on A/B/C (point estimate)": 0, "rule-of-three 95 % upper bound": 3}
for nm, row in size["panels"].items():
    if nm.startswith("abc"):
        continue
    scen[f"{nm}: union {row['union']} ({row['panel_on_k562ess_list']} K562-essential genes on the panel)"] = row["union"]
scen["mandate's planning number"] = 45
out["october"] = {}
ref = "dsallcp" if "dsallcp" in out["variants"] else VARS[0]
r = out["variants"][ref]
print(f"\n## October scaling ({ref}; H1-priced, board DE transfer ~0.5-1.0x)")
# October's covered targets can only be OFF the K562-essential list (Arc's criterion); H1 has 3 such targets
off_rate = SLOPE["nmae"] * r["covered_nmae_off_list"][0] / 6 / 300
out["per_covered_target_overall_nmae_only_offlist"] = off_rate
print("| scenario | covered targets | overall, all members | overall, nmae only | overall, nmae only at H1's off-list rate |")
print("|---|---:|---:|---:|---:|")
for nm, n in scen.items():
    a = n * r["per_covered_target_overall_all_members"]; b = n * r["per_covered_target_overall_nmae_only"]; c = n * off_rate
    out["october"][nm] = {"n_covered": n, "overall_all_members": a, "overall_nmae_only": b, "overall_nmae_only_offlist": c}
    print(f"| {nm} | {n} | {a:+.5f} | {b:+.5f} | {c:+.5f} |")
(HERE / "work/sum436.json").write_text(json.dumps(out, indent=1, default=float))

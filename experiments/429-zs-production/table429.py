"""429 table: six members per arm, Δ vs the same-seed base, mean ± sd over seeds, overall at board slopes; the zs
increment on top of the MDE; a paired target bootstrap for the nmae increments; the zs verify census.
usage: uv run python table429.py [SEEDS]   (default 0,1,2) -> work/table429.md, work/table429.json"""
import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
W = HERE / "work/h1"
SEEDS = [int(s) for s in (sys.argv[1] if len(sys.argv) > 1 else "0,1,2").split(",")]
MEM = ("pds", "mse", "nmae", "fid", "reach", "jac")
SLOPE = dict(pds=2.2138, mse=-1.051, nmae=-1.6532, fid=3.3645, reach=1.1292, jac=2.7263)
MDE = "mde0.1cap0.1min0.3"
ARMS = {  # name -> tag template
    "base": "base_s{s}_bps",
    "base + zs": "base_s{s}_bps_dnk250zs",
    "MDE": f"base_s{{s}}_bps_{MDE}",
    "MDE + zs": f"base_s{{s}}_bps_{MDE}_dnk250zs",
    "MDE + repair": f"base_s{{s}}_bps_zr-5_{MDE}",
    "MDE + repair + zs": f"base_s{{s}}_bps_zr-5_{MDE}_dnk250zs",
}
CONTRASTS = [  # (label, arm, reference arm)
    ("base + zs − base (zs alone)", "base + zs", "base"),
    ("MDE − base", "MDE", "base"),
    ("MDE + zs − base (the stack)", "MDE + zs", "base"),
    ("MDE + repair + zs − base (the stacked arm)", "MDE + repair + zs", "base"),
    ("MDE + zs − MDE (zs on top of the MDE)", "MDE + zs", "MDE"),
    ("MDE + repair − MDE (repair on H1)", "MDE + repair", "MDE"),
    ("MDE + repair + zs − MDE + zs (repair on H1, with zs)", "MDE + repair + zs", "MDE + zs"),
]


def raw(tag):
    return json.loads((W / f"metrics-{tag}.json").read_text())["raw"]


def per(tag):
    d = pl.read_csv(W / f"per-target-{tag}.csv").filter(pl.col("metric") == "de_wilcoxon_lfc_nmae")
    return dict(zip(d["perturbation"].to_list(), d["value"].to_list()))


def have(arm):
    return all((W / f"metrics-{ARMS[arm].format(s=s)}.json").exists() for s in SEEDS)


def delta(arm, ref):
    D = np.array([[raw(ARMS[arm].format(s=s))[m] - raw(ARMS[ref].format(s=s))[m] for m in MEM] for s in SEEDS])
    ov = (D * np.array([SLOPE[m] for m in MEM])).sum(1) / 6
    return D, ov


def boot(arm, ref, rng):
    ds = []
    for s in SEEDS:
        a, b = per(ARMS[arm].format(s=s)), per(ARMS[ref].format(s=s))
        ks = sorted(set(a) & set(b))
        ds.append(np.array([a[k] - b[k] for k in ks]))
    D = np.mean(ds, axis=0)
    bs = np.array([rng.choice(D, D.size).mean() for _ in range(20000)])
    gain = D.sum()
    return dict(targets=int(D.size), mean=float(D.mean()), ci=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                p_ge0=float(np.mean(bs >= 0)), better=int((D < -1e-6).sum()), worse=int((D > 1e-6).sum()),
                worse_gt_0005=int((D > 0.005).sum()),
                top3_share=float(np.sort(D)[:3].sum() / gain) if gain < 0 else None)


def census(arm):
    out = {}
    rows = []
    for s in SEEDS:
        rows += json.loads((W / f"diag-{ARMS[arm].format(s=s)}.json").read_text())
    for k in ("n_genes_U_changed", "n_touched_U_changed", "n_untouched_U_changed", "call_flips", "sign_flips"):
        out[k + "_sum"] = int(sum(r[k] for r in rows))
    for k in ("sums_equal", "libs_equal", "ballast_equal", "rp423_equal"):
        out[k + "_all"] = bool(all(r.get(k, False) for r in rows))
    out["blocks"] = len(rows)
    out["n_moved_mean"] = float(np.mean([r["n_moved"] for r in rows]))
    out["counts_moved_mean"] = float(np.mean([r["counts_moved"] for r in rows]))
    out["frac_massw_mean"] = float(np.mean([r["frac_massw"] for r in rows]))
    out["n_genes_tie_changed_mean"] = float(np.mean([r["n_genes_tie_changed"] for r in rows]))
    out["p_logmax_max"] = float(np.max([r["p_logmax"] for r in rows]))
    out["sign_guard_sum"] = int(sum(r["sign_guard"] for r in rows))
    out["sub_budget_exhausted_cells_mean"] = float(np.mean([r["sub_budget_exhausted_cells"] for r in rows]))
    return out


def main():
    rng = np.random.default_rng(429)
    res = {"seeds": SEEDS, "absolute": {}, "contrasts": {}, "bootstrap_nmae": {}, "census": {}}
    lines = [f"# 429 table (H1, production bulk `--bulk ps`, seeds {SEEDS})", ""]
    lines += ["## Absolute (seed mean)", "", "| arm | " + " | ".join(MEM) + " |", "|---|" + "---:|" * len(MEM)]
    for arm in ARMS:
        if not have(arm):
            continue
        v = np.mean([[raw(ARMS[arm].format(s=s))[m] for m in MEM] for s in SEEDS], axis=0)
        res["absolute"][arm] = dict(zip(MEM, map(float, v)))
        lines.append(f"| {arm} | " + " | ".join(f"{x:.5f}" for x in v) + " |")
    lines += ["", "## Δ raw (mean ± sd over seeds), overall = Σ slope·Δ / 6", "",
              "| contrast | " + " | ".join(f"Δ{m}" for m in MEM) + " | overall |", "|---|" + "---:|" * (len(MEM) + 1)]
    for lab, arm, ref in CONTRASTS:
        if not (have(arm) and have(ref)):
            continue
        D, ov = delta(arm, ref)
        sd = D.std(0, ddof=1) if len(SEEDS) > 1 else np.zeros(len(MEM))
        res["contrasts"][lab] = {"mean": dict(zip(MEM, map(float, D.mean(0)))), "sd": dict(zip(MEM, map(float, sd))),
                                 "per_seed": D.tolist(), "overall_mean": float(ov.mean()),
                                 "overall_sd": float(ov.std(ddof=1)) if len(SEEDS) > 1 else 0.0}
        cells = [f"{mu:+.5f} ± {sg:.5f}" if abs(mu) > 0 or sg > 0 else "0" for mu, sg in zip(D.mean(0), sd)]
        lines.append(f"| {lab} | " + " | ".join(cells) + f" | **{ov.mean():+.5f} ± {ov.std(ddof=1) if len(SEEDS) > 1 else 0:.5f}** |")
    if have("base + zs") and have("MDE + zs") and have("MDE"):
        dz, _ = delta("base + zs", "base")
        dm, _ = delta("MDE", "base")
        dmz, _ = delta("MDE + zs", "base")
        inter = dmz[:, 2] - dm[:, 2] - dz[:, 2]
        res["nmae_interaction"] = {"per_seed": inter.tolist(), "mean": float(inter.mean()),
                                   "zs_retained_on_mde": float((dmz[:, 2] - dm[:, 2]).mean() / dz[:, 2].mean())}
        lines += ["", f"nmae interaction (MDE+zs) − MDE − zs: {inter.mean():+.5f} (per seed {', '.join(f'{x:+.5f}' for x in inter)}); "
                      f"zs keeps {res['nmae_interaction']['zs_retained_on_mde']:.2f} of its stand-alone nmae gain on top of the MDE"]
    lines += ["", "## nmae paired target bootstrap (seed-averaged per-target Δ, 20,000 resamples)", "",
              "| contrast | targets | mean | 95% CI | p(≥0) | better / worse | worse by > 0.005 | top-3 share |",
              "|---|---:|---:|---|---:|---|---:|---:|"]
    for lab, arm, ref in CONTRASTS:
        if not (have(arm) and have(ref)) or "repair" in lab and "−" in lab and arm.endswith("repair"):
            continue
        b = boot(arm, ref, rng)
        res["bootstrap_nmae"][lab] = b
        lines.append(f"| {lab} | {b['targets']} | {b['mean']:+.5f} | [{b['ci'][0]:+.5f}, {b['ci'][1]:+.5f}] | {b['p_ge0']:.4f} | "
                     f"{b['better']} / {b['worse']} | {b['worse_gt_0005']} | "
                     f"{'' if b['top3_share'] is None else format(b['top3_share'], '.2f')} |")
    lines += ["", "## zs verify census (423's verify, per target block, all seeds)", ""]
    for arm in ("base + zs", "MDE + zs", "MDE + repair + zs"):
        if have(arm):
            c = census(arm)
            res["census"][arm] = c
            lines.append(f"- **{arm}**: " + ", ".join(f"{k} {v:.4g}" if isinstance(v, float) else f"{k} {v}" for k, v in c.items()))
    (HERE / "work/table429.md").write_text("\n".join(lines) + "\n")
    (HERE / "work/table429.json").write_text(json.dumps(res, indent=1))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

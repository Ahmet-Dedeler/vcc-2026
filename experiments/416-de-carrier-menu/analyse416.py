#!/usr/bin/env python3
"""416 analysis: per design, mean +- sd over emission seeds of the six raw members; deltas vs base (unpaired means and
seed-paired); predicted overall at board slopes; noise floor from base's seeds and from the pooled within-design seed
variance; paired target bootstrap (seed-averaged per-target values) for the DE-member sum.
usage: uv run python analyse416.py [--bulk base|ps]   -> work/analysis416[-ps].json + markdown table on stdout
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
WORK = HERE / "work" / "h1"
SLOPE = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532, "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}
MEM = ("pds", "mse", "nmae", "fid", "reach", "jac")
DE = ("nmae", "fid", "reach", "jac")
FULL = {"pds": "pds_cosine", "nmae": "de_wilcoxon_lfc_nmae", "fid": "de_wilcoxon_direction_fidelity_yield_raw",
        "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard"}


def ov(d, members=MEM):
    return sum(SLOPE[m] * d[m] for m in members) / 6.0


def load(bulk, seeds=None):
    """Menu runs only (tags DESIGN_sN[_bps]); K-matched controls (_K...) are analysed separately (analyse416k.py)."""
    suf = "_bps" if bulk == "ps" else ""
    runs = {}
    for f in sorted(WORK.glob("metrics-*_s*.json")):
        m = re.fullmatch(r"metrics-([a-z0-9]+(?:_off)?)_s(\d+)" + re.escape(suf) + r"\.json", f.name)
        if not m or (seeds is not None and int(m.group(2)) not in seeds):
            continue
        runs.setdefault(m.group(1), {})[int(m.group(2))] = json.loads(f.read_text())
    return runs


def per_target(tag):
    df = pl.read_csv(WORK / f"per-target-{tag}.csv")
    out = {}
    names = sorted(set(df["perturbation"].to_list()))       # 85 targets; nmae omits 27 (NaN here, as the scorer's nanmean)
    for short in DE:
        sub = df.filter(pl.col("metric") == FULL[short])
        m = dict(zip(sub["perturbation"].to_list(), sub["value"].to_numpy().astype(float)))
        out[short] = (names, np.array([m.get(t, np.nan) for t in names], dtype=float))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bulk", default="base")
    ap.add_argument("--seeds", default=None, help="comma list, e.g. 0,1,2 (balanced menu); default all")
    a = ap.parse_args()
    runs = load(a.bulk, None if a.seeds is None else {int(x) for x in a.seeds.split(",")})
    base = runs.get("base", {})
    seeds_b = sorted(base)
    res = {"bulk": a.bulk, "slopes": SLOPE, "designs": {}}
    # pooled within-design seed variance (all designs with >= 2 seeds)
    pooled = {m: [] for m in MEM}; pooled_de = []; pooled_all = []
    for d, rs in runs.items():
        if len(rs) < 2:
            continue
        X = {m: np.array([rs[s]["raw"][m] for s in sorted(rs)]) for m in MEM}
        for m in MEM:
            pooled[m].append((X[m] - X[m].mean()) ** 2)
        de = np.array([ov(rs[s]["raw"], DE) for s in sorted(rs)])
        al = np.array([ov(rs[s]["raw"]) for s in sorted(rs)])
        pooled_de.append((de - de.mean()) ** 2); pooled_all.append((al - al.mean()) ** 2)
    def psd(lst):
        if not lst:
            return float("nan")
        ss = sum(float(x.sum()) for x in lst); dof = sum(len(x) - 1 for x in lst)
        return float(np.sqrt(ss / dof)) if dof else float("nan")
    res["noise_pooled_sd_single_seed"] = {**{m: psd(pooled[m]) for m in MEM}, "DE_overall": psd(pooled_de),
                                          "overall": psd(pooled_all),
                                          "dof": sum(len(x) - 1 for x in pooled_de)}
    if len(seeds_b) >= 2:
        bde = np.array([ov(base[s]["raw"], DE) for s in seeds_b])
        res["noise_base_sd_single_seed"] = {**{m: float(np.std([base[s]["raw"][m] for s in seeds_b], ddof=1)) for m in MEM},
                                            "DE_overall": float(bde.std(ddof=1)), "n": len(seeds_b)}
    bmean = {m: float(np.mean([base[s]["raw"][m] for s in seeds_b])) for m in MEM} if seeds_b else None
    bpt = None
    if seeds_b:
        bpt = {}
        for s in seeds_b:
            pt = per_target(base[s]["tag"])
            for k, (names, v) in pt.items():
                bpt.setdefault(k, []).append(v)
        bpt = {k: np.mean(v, 0) for k, v in bpt.items()}
    rng = np.random.default_rng(416)
    for d, rs in sorted(runs.items()):
        seeds = sorted(rs)
        X = {m: np.array([rs[s]["raw"][m] for s in seeds]) for m in MEM}
        rec = {"seeds": seeds, "mean": {m: float(X[m].mean()) for m in MEM},
               "sd": {m: (float(X[m].std(ddof=1)) if len(seeds) > 1 else float("nan")) for m in MEM},
               "K_hist": [rs[s].get("K_hist") for s in seeds], "clip_mass_total": [rs[s].get("clip_mass_total") for s in seeds],
               "prior_lambda_median": [rs[s].get("prior_lambda_median") for s in seeds],
               "prior_flips_median": [rs[s].get("prior_flips_median") for s in seeds],
               "seconds": [round(rs[s].get("seconds", 0)) for s in seeds]}
        if bmean is not None:
            dm = {m: rec["mean"][m] - bmean[m] for m in MEM}
            rec["delta_mean"] = dm
            rec["pred_overall"] = {"DE": ov(dm, DE), "bulk": ov(dm, ("pds", "mse")), "total": ov(dm),
                                   **{m: SLOPE[m] * dm[m] / 6 for m in MEM}}
            common = [s for s in seeds if s in base]
            if common:
                pde = np.array([ov({m: rs[s]["raw"][m] - base[s]["raw"][m] for m in MEM}, DE) for s in common])
                rec["paired_DE_overall"] = {"per_seed": pde.tolist(), "mean": float(pde.mean()),
                                            "sd": float(pde.std(ddof=1)) if len(common) > 1 else float("nan")}
            # paired target bootstrap on seed-averaged per-target DE members (target sampling, not realisation)
            try:
                pt = {}
                for s in seeds:
                    for k, (names, v) in per_target(rs[s]["tag"]).items():
                        pt.setdefault(k, []).append(v)
                pt = {k: np.mean(v, 0) for k, v in pt.items()}
                n = len(pt["fid"])
                diffs = {k: pt[k] - bpt[k] for k in DE}
                B = 4000
                idx = rng.integers(0, n, size=(B, n))
                boot = np.zeros(B)
                for k in DE:
                    dk = diffs[k]
                    ok = np.isfinite(dk)
                    # nanmean per bootstrap sample, as the scorer's nanmean over targets
                    vals = np.where(ok[idx], dk[idx], 0.0).sum(1) / np.maximum(ok[idx].sum(1), 1)
                    boot += SLOPE[k] * vals / 6.0
                rec["target_boot_DE_overall"] = {"ci95": np.quantile(boot, [0.025, 0.975]).tolist(), "p_le_0": float((boot <= 0).mean())}
            except Exception as exc:  # noqa: BLE001
                rec["target_boot_DE_overall"] = {"error": str(exc)}
        res["designs"][d] = rec
    nz = res["noise_pooled_sd_single_seed"]["DE_overall"]
    for d, rec in res["designs"].items():
        if d == "base" or "pred_overall" not in rec:
            continue
        k = len(rec["seeds"]); kb = len(seeds_b)
        se = nz * np.sqrt(1.0 / k + 1.0 / kb)
        rec["flag"] = {"gain_gt_2sd_single_seed": bool(rec["pred_overall"]["DE"] > 2 * nz),
                       "gain_gt_2se_of_mean_diff": bool(rec["pred_overall"]["DE"] > 2 * se), "se_mean_diff": float(se)}
    out = WORK.parent / (("analysis416" if a.bulk == "base" else "analysis416-ps")
                         + ("" if a.seeds is None else "-seeds" + a.seeds.replace(",", "")) + ".json")
    out.write_text(json.dumps(res, indent=1))
    # markdown
    print(f"noise (pooled within-design single-seed sd, dof {res['noise_pooled_sd_single_seed']['dof']}): " +
          ", ".join(f"{m} {res['noise_pooled_sd_single_seed'][m]:.5f}" for m in MEM) +
          f"; DE-sum overall {nz:.5f}; total overall {res['noise_pooled_sd_single_seed']['overall']:.5f}")
    if "noise_base_sd_single_seed" in res:
        nb = res["noise_base_sd_single_seed"]
        print("noise (base seeds only, n=%d): " % nb["n"] + ", ".join(f"{m} {nb[m]:.5f}" for m in MEM) + f"; DE-sum overall {nb['DE_overall']:.5f}")
    print("\n| design | seeds | " + " | ".join(MEM) + " | Δ DE-sum overall | Δ bulk overall | paired DE (sd) | target-boot 95% | flag 2sd / 2se |")
    print("|---|---|" + "---|" * (len(MEM) + 5))
    order = sorted(res["designs"], key=lambda d: -res["designs"][d].get("pred_overall", {}).get("DE", 0.0) if d != "base" else -1e9)
    for d in order:
        rec = res["designs"][d]
        cells = [f"{rec['mean'][m]:.5f} ± {rec['sd'][m]:.5f}" if d == "base" else
                 f"{rec['delta_mean'][m]:+.5f} (±{rec['sd'][m]:.5f})" for m in MEM]
        po = rec.get("pred_overall", {})
        pdd = rec.get("paired_DE_overall", {})
        tb = rec.get("target_boot_DE_overall", {})
        fl = rec.get("flag", {})
        print(f"| {d} | {len(rec['seeds'])} | " + " | ".join(cells) +
              f" | {po.get('DE', float('nan')):+.5f} | {po.get('bulk', float('nan')):+.5f} | "
              f"{pdd.get('mean', float('nan')):+.5f} ({pdd.get('sd', float('nan')):.5f}) | "
              f"{tb.get('ci95', [float('nan')] * 2)[0]:+.5f}…{tb.get('ci95', [float('nan')] * 2)[1]:+.5f} | "
              f"{'Y' if fl.get('gain_gt_2sd_single_seed') else '-'} / {'Y' if fl.get('gain_gt_2se_of_mean_diff') else '-'} |")


if __name__ == "__main__":
    main()

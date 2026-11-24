#!/usr/bin/env python3
"""416 secondary analyses (seed-paired, seeds 0-2): (1) K-matched control -- base carriers forced to g000's per-target
ballast count K; (2) the 2x2 {whitened 0.65, unwhitened} x {CD4, no CD4} on the carriers; (3) g000's gain by n_real
stratum (seed-averaged per-target values); (4) sign accuracy of the delivered LFC on the reference's DE genes.
Writes work/analysis416k.json and prints a summary.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
WORK = HERE / "work" / "h1"
REF = HERE.parents[1] / "experiments/323-current-family-validation/work"
S = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532, "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}
M = tuple(S); DE = ("nmae", "fid", "reach", "jac")
FULL = {"nmae": "de_wilcoxon_lfc_nmae", "fid": "de_wilcoxon_direction_fidelity_yield_raw",
        "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard"}
SEEDS = (0, 1, 2)


def raw(tag):
    return json.loads((WORK / f"metrics-{tag}.json").read_text())["raw"]


def ov(d, mem=M):
    return sum(S[m] * d[m] for m in mem) / 6


def paired(a_tags, b_tags):
    ds = [{m: raw(a)[m] - raw(b)[m] for m in M} for a, b in zip(a_tags, b_tags)]
    mean = {m: float(np.mean([d[m] for d in ds])) for m in M}
    return {"delta_raw": mean, "DE_overall": ov(mean, DE), "DE_overall_noreach": ov(mean, ("nmae", "fid", "jac")),
            "bulk_overall": ov(mean, ("pds", "mse")), "DE_overall_per_seed": [ov(d, DE) for d in ds]}


def per_target(tag):
    df = pl.read_csv(WORK / f"per-target-{tag}.csv")
    return {k: dict(zip(df.filter(pl.col("metric") == f)["perturbation"].to_list(),
                        df.filter(pl.col("metric") == f)["value"].to_list())) for k, f in FULL.items()}


def main():
    out = {}
    out["K_matched_minus_base"] = paired([f"base_s{s}_Kg000_s{s}" for s in SEEDS], [f"base_s{s}" for s in SEEDS])
    out["g000_minus_K_matched"] = paired([f"g000_s{s}" for s in SEEDS], [f"base_s{s}_Kg000_s{s}" for s in SEEDS])
    cells = {"whitened+CD4 (base)": "base", "whitened, no CD4": "w065nocd4", "unwhitened, no CD4 (pre)": "pre",
             "unwhitened+CD4 (g000)": "g000"}
    out["factorial_vs_base"] = {k: paired([f"{v}_s{s}" for s in SEEDS], [f"base_s{s}" for s in SEEDS]) for k, v in cells.items()
                                if v != "base"}
    out["effect_CD4_at_w065"] = paired([f"base_s{s}" for s in SEEDS], [f"w065nocd4_s{s}" for s in SEEDS])
    out["effect_CD4_at_g000"] = paired([f"g000_s{s}" for s in SEEDS], [f"pre_s{s}" for s in SEEDS])
    out["effect_unwhiten_with_CD4"] = paired([f"g000_s{s}" for s in SEEDS], [f"base_s{s}" for s in SEEDS])
    out["effect_unwhiten_no_CD4"] = paired([f"pre_s{s}" for s in SEEDS], [f"w065nocd4_s{s}" for s in SEEDS])
    # n_real strata
    nreal = json.loads((REF / "reference.json").read_text())["n_real_by_target"]
    T = sorted(nreal)
    def avg(design):
        runs = [per_target(f"{design}_s{s}") for s in SEEDS]
        return {k: np.array([np.mean([r[k].get(t, np.nan) for r in runs]) for t in T]) for k in FULL}
    b, g = avg("base"), avg("g000")
    nr = np.array([nreal[t] for t in T])
    strata = {}
    for lo, hi in ((0, 50), (50, 500), (500, 2000), (2000, 10 ** 9)):
        msk = (nr >= lo) & (nr < hi)
        strata[f"[{lo},{hi})"] = {"n": int(msk.sum()), **{k: float(np.nanmean((g[k] - b[k])[msk])) for k in FULL}}
    out["g000_strata_n_real"] = strata
    d = g["nmae"] - b["nmae"]; ok = np.isfinite(d)
    out["g000_targets"] = {"fid_up": int((g["fid"] > b["fid"]).sum()), "fid_down": int((g["fid"] < b["fid"]).sum()),
                           "nmae_better": int((d[ok] < 0).sum()), "nmae_worse": int((d[ok] > 0).sum())}
    # sign accuracy of delivered LFC on the reference's DE genes (own target excluded)
    ref = np.load(REF / "reference.npz")
    rl, rq = np.asarray(ref["lfc"], float), np.asarray(ref["q"], float)
    own = np.equal.outer(ref["targets"].astype(str), ref["genes"].astype(str)[ref["keep"]])
    de = (rq < 0.05) & ~own & np.isfinite(rl)
    acc = {}
    for design in ("base", "g000", "g030", "g090", "gm030", "pre", "w065nocd4", "signrep", "gain"):
        v = []
        for s in SEEDS:
            z = np.load(WORK / f"stats-{design}_s{s}.npz")
            l = z["lfc"]; okk = np.isfinite(l) & de
            v.append(float(((np.sign(l) == np.sign(rl)) & okk).sum() / okk.sum()))
        acc[design] = {"mean": float(np.mean(v)), "per_seed": v}
    out["sign_accuracy_on_reference_DE"] = acc
    (HERE / "work/analysis416k.json").write_text(json.dumps(out, indent=1))
    for k in ("K_matched_minus_base", "g000_minus_K_matched", "effect_CD4_at_w065", "effect_CD4_at_g000",
              "effect_unwhiten_with_CD4", "effect_unwhiten_no_CD4"):
        r = out[k]
        print(f"{k:28s} " + " ".join(f"{m} {r['delta_raw'][m]:+.5f}" for m in M) +
              f" | DE {r['DE_overall']:+.5f} (no reach {r['DE_overall_noreach']:+.5f}) bulk {r['bulk_overall']:+.5f}")
    for k, v in strata.items():
        print("n_real", k, json.dumps(v))
    print(json.dumps(out["g000_targets"]))
    print("sign accuracy on reference DE genes:", {k: round(v["mean"], 4) for k, v in acc.items()})


if __name__ == "__main__":
    main()

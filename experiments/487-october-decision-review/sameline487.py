"""487 decision 1: how robust is the same-line rows' H1 price (+0.0934 at w 2 + pt, 35 % coverage; 450 §9/§10)?

Reads only saved per-target outputs of 450's H1 emitter runs and 435's base (no builds, no scorer). Peak RSS < 1 GB.
- KA: per-arm overall in 450's board units (435's SLOPE, paired by seed) must reproduce read450.py's table.
- Per-target decomposition: covered (Pan's 30) vs uncovered, and covered essential vs covered non-essential
  (the 2026 panel is essential-free: 0 of 300 on the K562-essential list).
- Target bootstrap (4,000 draws, all 85 targets resampled; mse recomputed as a ratio of sums from the legs, nmae over
  its eligible targets), and a stratified version that keeps the covered share fixed.
- Coverage-share curve: synthetic 85-target panels with a share c of covered targets (drawn with replacement from
  the covered set, the rest from the uncovered set), optionally essential-free.
- Dose (w 1 / 2 / 4, 3 seeds) and bank quality (flip 0.2 / 0.3, seed 0): paired per-target contrasts, bootstrapped.
usage: uv run python experiments/487-october-decision-review/sameline487.py  -> work/sameline487.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
W450 = ROOT / "experiments/450-same-line-programme/work"
BASE = ROOT / "experiments/435-zs-on-coarse/work/h1"
TAG = "base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse"
SLOPE = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532, "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}  # 435 / read450
COL = {"pds": "pds_cosine", "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "nmae": "de_wilcoxon_lfc_nmae", "jac": "de_wilcoxon_sig_jaccard"}
NB = 4000


def load(d, s):
    """One run: per-target table (index = targets) with the five per-target members plus mse numerator/denominator."""
    p = d / f"per-target-{TAG.format(s=s)}.csv"
    m = d / f"metrics-{TAG.format(s=s)}.json"
    if not (p.exists() and m.exists()):
        return None
    df = pd.read_csv(p).pivot_table(index="perturbation", columns="metric", values="value")
    meta = json.loads(m.read_text())
    legs = meta["mse_legs"]
    df["mse_num"] = pd.Series(legs["numerator"])
    df["mse_den"] = pd.Series(legs["denominator"])
    df.attrs["raw"] = meta["raw"]
    return df


def members(df, idx):
    """Raw members over the target multiset idx (positions into df)."""
    out = {k: float(np.nanmean(df[COL[k]].values[idx])) for k in COL}
    out["mse"] = float(df["mse_num"].values[idx].sum() / df["mse_den"].values[idx].sum())
    return out


def overall(arm, base, idx):
    a, b = members(arm, idx), members(base, idx)
    return sum(SLOPE[k] * (a[k] - b[k]) for k in SLOPE) / 6, {k: a[k] - b[k] for k in SLOPE}


def main():
    ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
    targets = ref["targets"].astype(str)
    n_real = (ref["q"] < 0.05).sum(axis=1)
    pan = set(Path(ROOT / "data/external/gse295214/Pan2025_hESC_post_qc_targets.txt").read_text().split())
    ess = set(Path(ROOT / "data/donors/k562_essential/perts.txt").read_text().split())
    bases = {s: load(BASE, s) for s in (0, 1, 2)}
    order = bases[0].index.values
    assert set(order) == set(targets), "target sets differ"
    nr = pd.Series(n_real, index=targets).reindex(order).values
    cov = np.isin(order, list(pan))
    es = np.isin(order, list(ess))
    out = {"n_targets": int(len(order)), "covered": int(cov.sum()), "essential": int(es.sum()),
           "covered_essential": int((cov & es).sum()), "covered_nonessential": int((cov & ~es).sum()),
           "uncovered_essential": int((~cov & es).sum()),
           "median_n_real": {"covered_ess": float(np.median(nr[cov & es])), "covered_noness": float(np.median(nr[cov & ~es])),
                             "uncovered_noness": float(np.median(nr[~cov & ~es]))}}
    print(json.dumps(out))
    # mse denominator share (truth energy) by stratum: the ratio of sums is dominated by strong targets
    den = bases[0]["mse_den"].values
    out["mse_den_share"] = {"covered_ess": float(den[cov & es].sum() / den.sum()),
                            "covered_noness": float(den[cov & ~es].sum() / den.sum()),
                            "uncovered": float(den[~cov].sum() / den.sum())}
    print("mse denominator share", out["mse_den_share"])

    arms = {"rows w0.5": "h1b_w0.5_clean4", "rows w1": "h1b_w1_clean4", "rows w2": "h1b_w2_clean4",
            "rows w4": "h1b_w4_clean4", "rows w32": "h1b_w32_clean4",
            "full w1 (pt)": "h1c_w1_clean4", "full w2 (pt)": "h1c_w2_clean4", "full w4 (pt)": "h1c_w4_clean4",
            "flip0.2 w0.5": "h1c_w0.5_clean4_f0.2s0", "flip0.2 w1": "h1c_w1_clean4_f0.2s0", "flip0.2 w2": "h1c_w2_clean4_f0.2s0",
            "flip0.3 w0.5": "h1c_w0.5_clean4_f0.3s0", "flip0.3 w1": "h1c_w1_clean4_f0.3s0", "flip0.3 w2": "h1c_w2_clean4_f0.3s0"}
    runs = {k: {s: load(W450 / d, s) for s in (0, 1, 2)} for k, d in arms.items()}
    runs = {k: {s: v for s, v in r.items() if v is not None} for k, r in runs.items()}
    for k in runs:
        for s in runs[k]:
            runs[k][s] = runs[k][s].reindex(order)
    for s in bases:
        bases[s] = bases[s].reindex(order)
    allidx = np.arange(len(order))

    # --- KA: full-panel overall per arm (mean over its seeds), as read450.py prints it
    ka = {}
    for k, r in runs.items():
        ov = [overall(r[s], bases[s], allidx)[0] for s in sorted(r)]
        ka[k] = {"seeds": sorted(r), "overall": round(float(np.mean(ov)), 4),
                 "sd": round(float(np.std(ov, ddof=1)), 4) if len(ov) > 1 else None}
    out["ka_overall"] = ka
    print("KA (read450 reproduces: full w2 +0.0934 ± 0.0002, full w1 +0.0837, full w4 +0.0933 ± 0.0005, rows w2 +0.0790):")
    for k, v in ka.items():
        print(f"  {k:16s} seeds {v['seeds']}  overall {v['overall']:+.4f}  sd {v['sd']}")

    # --- decomposition of the w2 full-mode gain by stratum (per-target deltas; mse via its legs over the full denominator)
    def strata_contrib(arm_key, s):
        a, b = runs[arm_key][s], bases[s]
        res = {}
        for name, mask in {"covered_ess": cov & es, "covered_noness": cov & ~es, "uncovered": ~cov}.items():
            c = {}
            for k in ("pds", "fid", "reach", "jac"):
                c[k] = float((a[COL[k]].values[mask] - b[COL[k]].values[mask]).sum() / len(order))
            elig = np.isfinite(b[COL["nmae"]].values) & np.isfinite(a[COL["nmae"]].values)
            c["nmae"] = float((a[COL["nmae"]].values[mask & elig] - b[COL["nmae"]].values[mask & elig]).sum() / elig.sum())
            c["mse"] = float((a["mse_num"].values[mask] - b["mse_num"].values[mask]).sum() / b["mse_den"].values.sum())
            c["overall"] = sum(SLOPE[k] * c[k] for k in SLOPE) / 6
            res[name] = c
        return res
    dec = {s: strata_contrib("full w2 (pt)", s) for s in runs["full w2 (pt)"]}
    out["decomp_full_w2"] = {name: {k: float(np.mean([dec[s][name][k] for s in dec])) for k in list(SLOPE) + ["overall"]}
                             for name in ("covered_ess", "covered_noness", "uncovered")}
    print("full w2 gain split by stratum (contribution to the 85-target overall, mean of seeds):")
    for name, c in out["decomp_full_w2"].items():
        print(f"  {name:15s} " + " ".join(f"{k} {c[k]:+.4f}" for k in list(SLOPE) + ["overall"]))

    # per-covered-target gain (overall-equivalent per target, mse per-target normalised by its own share)
    rng = np.random.default_rng(487)

    def boot(arm_key, pool_idx, nb=NB, seeds=None):
        seeds = seeds or sorted(runs[arm_key])
        vals = []
        for s in seeds:
            a, b = runs[arm_key][s], bases[s]
            pt = overall(a, b, pool_idx)[0]
            bs = np.array([overall(a, b, pool_idx[rng.integers(0, len(pool_idx), len(pool_idx))])[0] for _ in range(nb // len(seeds))])
            vals.append((pt, bs))
        pts = float(np.mean([v[0] for v in vals]))
        bs = np.concatenate([v[1] for v in vals])
        return {"point": round(pts, 4), "ci95": [round(float(np.quantile(bs, .025)), 4), round(float(np.quantile(bs, .975)), 4)],
                "p_le_0": round(float((bs <= 0).mean()), 4)}

    out["bootstrap_all85"] = {k: boot(k, allidx) for k in ("full w2 (pt)", "full w1 (pt)", "full w4 (pt)", "rows w2")}
    print("target bootstrap, all 85 resampled:", json.dumps(out["bootstrap_all85"]))

    # stratified: keep 30 covered / 55 uncovered shares fixed
    def strat_boot(arm_key, groups, nb=NB):
        seeds = sorted(runs[arm_key])
        bs = []
        for s in seeds:
            a, b = runs[arm_key][s], bases[s]
            for _ in range(nb // len(seeds)):
                idx = np.concatenate([g[rng.integers(0, len(g), len(g))] for g in groups])
                bs.append(overall(a, b, idx)[0])
        bs = np.array(bs)
        return {"ci95": [round(float(np.quantile(bs, .025)), 4), round(float(np.quantile(bs, .975)), 4)], "p_le_0": round(float((bs <= 0).mean()), 4)}
    out["bootstrap_stratified"] = {k: strat_boot(k, [np.flatnonzero(cov), np.flatnonzero(~cov)]) for k in ("full w2 (pt)",)}
    print("stratified (covered/uncovered fixed):", out["bootstrap_stratified"])

    # --- essential-free H1 stand-in: only the 67 non-essential targets (covered non-essential + uncovered non-essential)
    ne = np.flatnonzero(~es)
    out["nonessential_only"] = {k: boot(k, ne) for k in ("full w2 (pt)", "full w1 (pt)", "full w4 (pt)", "rows w2")}
    for k in out["nonessential_only"]:
        s0 = sorted(runs[k])
        mem = {m: float(np.mean([overall(runs[k][s], bases[s], ne)[1][m] for s in s0])) for m in SLOPE}
        out["nonessential_only"][k]["members"] = {m: round(v, 4) for m, v in mem.items()}
    print("non-essential targets only (67; covered", int((cov & ~es).sum()), "):", json.dumps(out["nonessential_only"]))

    # --- coverage-share curve on synthetic 85-target panels
    def coverage_curve(arm_key, cpool, upool, shares=(0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0), nb=1000, n=85):
        seeds = sorted(runs[arm_key])
        res = {}
        for c in shares:
            k = int(round(c * n))
            vals, mems = [], []
            for s in seeds:
                a, b = runs[arm_key][s], bases[s]
                for _ in range(nb // len(seeds)):
                    idx = np.concatenate([cpool[rng.integers(0, len(cpool), k)], upool[rng.integers(0, len(upool), n - k)]])
                    o, m = overall(a, b, idx)
                    vals.append(o); mems.append(m)
            vals = np.array(vals)
            res[f"{c:.2f}"] = {"mean": round(float(vals.mean()), 4), "ci95": [round(float(np.quantile(vals, .025)), 4), round(float(np.quantile(vals, .975)), 4)],
                               "per_unit_share": round(float(vals.mean() / c), 4),
                               "mse": round(float(np.mean([m["mse"] for m in mems])), 4)}
        return res
    out["coverage_curve_all"] = coverage_curve("full w2 (pt)", np.flatnonzero(cov), np.flatnonzero(~cov))
    out["coverage_curve_essfree"] = coverage_curve("full w2 (pt)", np.flatnonzero(cov & ~es), np.flatnonzero(~cov & ~es))
    print("coverage curve (all targets):", json.dumps(out["coverage_curve_all"]))
    print("coverage curve (essential-free):", json.dumps(out["coverage_curve_essfree"]))

    # --- dose contrasts, paired per seed and target, bootstrapped over targets
    def contrast(k1, k2, pool, nb=NB):
        seeds = sorted(set(runs[k1]) & set(runs[k2]))
        pts, bs = [], []
        for s in seeds:
            a, b = runs[k1][s], runs[k2][s]
            pts.append(overall(a, b, pool)[0])
            bs += [overall(a, b, pool[rng.integers(0, len(pool), len(pool))])[0] for _ in range(nb // len(seeds))]
        bs = np.array(bs)
        return {"seeds": seeds, "point": round(float(np.mean(pts)), 4), "ci95": [round(float(np.quantile(bs, .025)), 4), round(float(np.quantile(bs, .975)), 4)],
                "p_le_0": round(float((bs <= 0).mean()), 4)}
    out["dose_contrasts_all85"] = {"w2 - w1": contrast("full w2 (pt)", "full w1 (pt)", allidx),
                                   "w4 - w2": contrast("full w4 (pt)", "full w2 (pt)", allidx),
                                   "flip0.2: w2 - w0.5": contrast("flip0.2 w2", "flip0.2 w0.5", allidx),
                                   "flip0.3: w2 - w0.5": contrast("flip0.3 w2", "flip0.3 w0.5", allidx),
                                   "flip0.2: w1 - w0.5": contrast("flip0.2 w1", "flip0.2 w0.5", allidx),
                                   "clean: w2(pt) - w0.5(rows only)": contrast("full w2 (pt)", "rows w0.5", allidx)}
    out["dose_contrasts_essfree"] = {k: contrast(*v, ne) for k, v in {
        "w2 - w1": ("full w2 (pt)", "full w1 (pt)"), "w4 - w2": ("full w4 (pt)", "full w2 (pt)"),
        "flip0.2: w2 - w0.5": ("flip0.2 w2", "flip0.2 w0.5"), "flip0.3: w2 - w0.5": ("flip0.3 w2", "flip0.3 w0.5")}.items()}
    print("dose contrasts (85):", json.dumps(out["dose_contrasts_all85"]))
    print("dose contrasts (essential-free 67):", json.dumps(out["dose_contrasts_essfree"]))
    # flip curve overall on both pools (seed 0), vs base
    fl = {}
    for k in [x for x in runs if x.startswith("flip") or x in ("rows w0.5", "full w1 (pt)", "full w2 (pt)")]:
        fl[k] = {"all85": round(overall(runs[k][0], bases[0], allidx)[0], 4), "essfree67": round(overall(runs[k][0], bases[0], ne)[0], 4)}
    out["flip_curve_seed0"] = fl
    print("flip curve seed 0:", json.dumps(fl))
    (HERE / "work").mkdir(exist_ok=True)
    (HERE / "work/sameline487.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

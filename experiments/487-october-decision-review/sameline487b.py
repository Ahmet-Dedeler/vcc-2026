"""487 decision 1b: per covered target, does the same-line gain depend on the target's response strength (n_real)?
Pan's 30 covered h1-test targets are strong responders (median n_real 537 non-essential / 1,862 essential, against 13
for the uncovered non-essential targets). A genome-wide October bank (Orion HCT116, K562-GWPS) covers every target,
weak ones included, so the H1 per-target gain must be read by strength before it is extrapolated to 100 % coverage.
Per-target overall-equivalent: the DE members' per-target deltas / 85 (nmae / its eligible count), pds / 85, and the
mse numerator delta over the base's full denominator (the ratio-of-sums contribution). usage: uv run python ... -> work/sameline487b.json
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
import importlib.util
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
spec = importlib.util.spec_from_file_location("s487", HERE / "sameline487.py"); S = importlib.util.module_from_spec(spec); spec.loader.exec_module(S)
ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
targets = ref["targets"].astype(str); n_real = pd.Series((ref["q"] < 0.05).sum(axis=1), index=targets)
pan = set((ROOT / "data/external/gse295214/Pan2025_hESC_post_qc_targets.txt").read_text().split())
ess = set((ROOT / "data/donors/k562_essential/perts.txt").read_text().split())
out = {}
for arm in ("h1c_w2_clean4", "h1c_w1_clean4", "h1b_w2_clean4", "h1c_w2_clean4_f0.2s0", "h1c_w2_clean4_f0.3s0", "h1c_w0.5_clean4_f0.3s0"):
    rows = []
    for s in (0, 1, 2):
        a = S.load(S.W450 / arm, s)
        if a is None:
            continue
        b = S.load(S.BASE, s).reindex(a.index)
        n = len(a); elig = np.isfinite(b[S.COL["nmae"]]).sum()
        per = pd.DataFrame(index=a.index)
        for k in ("pds", "fid", "reach", "jac"):
            per[k] = (a[S.COL[k]] - b[S.COL[k]]).fillna(0) / n
        per["nmae"] = (a[S.COL["nmae"]] - b[S.COL["nmae"]]).fillna(0) / elig
        per["mse"] = (a["mse_num"] - b["mse_num"]) / b["mse_den"].sum()
        per["ov"] = sum(S.SLOPE[k] * per[k] for k in S.SLOPE) / 6
        per["ov_de"] = sum(S.SLOPE[k] * per[k] for k in ("fid", "reach", "nmae", "jac")) / 6
        per["ov_bulk"] = sum(S.SLOPE[k] * per[k] for k in ("pds", "mse")) / 6
        per["seed"] = s
        rows.append(per)
    P = pd.concat(rows).groupby(level=0).mean(numeric_only=True)
    P["n_real"] = n_real.reindex(P.index); P["cov"] = P.index.isin(pan); P["ess"] = P.index.isin(ess)
    C = P[P["cov"]].copy()
    C["tercile"] = pd.qcut(C["n_real"].rank(method="first"), 3, labels=["weak", "mid", "strong"])
    tab = C.groupby("tercile", observed=True).agg(n=("ov", "size"), n_real_med=("n_real", "median"), ov_sum=("ov", "sum"),
                                                 ov_de_sum=("ov_de", "sum"), ov_bulk_sum=("ov_bulk", "sum"),
                                                 n_neg=("ov", lambda x: int((x < 0).sum())))
    from scipy.stats import spearmanr
    rho = spearmanr(C["n_real"], C["ov"] * 85).statistic
    rho_de = spearmanr(C["n_real"], C["ov_de"] * 85).statistic
    out[arm] = {"by_tercile": json.loads(tab.round(5).to_json(orient="index")), "spearman_nreal_ov": round(float(rho), 3),
                "spearman_nreal_ov_de": round(float(rho_de), 3),
                "covered_sum": round(float(C["ov"].sum()), 5), "uncovered_sum": round(float(P.loc[~P["cov"], "ov"].sum()), 5),
                "per_target": {t: {"n_real": int(C.loc[t, "n_real"]), "ess": bool(C.loc[t, "ess"]), "ov_x85": round(float(C.loc[t, "ov"] * 85), 4),
                                   "de_x85": round(float(C.loc[t, "ov_de"] * 85), 4), "bulk_x85": round(float(C.loc[t, "ov_bulk"] * 85), 4)}
                               for t in C.sort_values("n_real").index}}
    print(arm, "Spearman(n_real, per-target gain) =", out[arm]["spearman_nreal_ov"], "(DE part", out[arm]["spearman_nreal_ov_de"], ")",
          "covered sum", out[arm]["covered_sum"], "uncovered sum", out[arm]["uncovered_sum"])
    print(tab.round(4).to_string())
print("\nper covered target, w2 full (x85 = its gain if every target behaved like it):")
for t, v in out["h1c_w2_clean4"]["per_target"].items():
    print(f"  {t:10s} n_real {v['n_real']:5d} ess {int(v['ess'])}  total {v['ov_x85']:+.3f}  DE {v['de_x85']:+.3f}  bulk {v['bulk_x85']:+.3f}")
(HERE / "work/sameline487b.json").write_text(json.dumps(out, indent=1))

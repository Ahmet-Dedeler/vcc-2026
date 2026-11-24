"""288 -- coordinator check prompted by 287 (T3): does F8's 'conservation predicts per-target pds under LOLO' survive
matching on response MAGNITUDE?  Inputs: 276's conservation276.csv (C_lolo_<fold>, pds_<fold>, C_raw, half_A/half_B),
280's desc280_per_target.csv (GWCD4i n_total_de_genes per condition = an independent magnitude proxy), 276's rows276.npz
(row norms if present).  Spearman overall, within magnitude tertiles, and the partial Spearman given magnitude."""
import json, numpy as np, pandas as pd
from scipy.stats import spearmanr
rng = np.random.default_rng(20260906)
c = pd.read_csv("experiments/276-conservation-map/work/conservation276.csv")
d = pd.read_csv("experiments/280-gwcd4i-donor-replicates/work/desc280_per_target.csv")
print("cons cols:", list(c.columns)); print("desc cols:", list(d.columns)[:40])
magcol = next((k for k in d.columns if "n_total_de" in k or k.startswith("n_de")), None)
tcol = next((k for k in d.columns if k in ("target","gene","pert")), d.columns[0])
print("magnitude column:", magcol, "target col:", tcol)
dd = d.groupby(tcol)[magcol].median() if magcol else None
c["mag_de"] = c["target"].map(dd) if dd is not None else np.nan
# a second magnitude proxy: reliability-independent row norm from rows276.npz if present
try:
    z = np.load("experiments/276-conservation-map/work/rows276.npz", allow_pickle=True); print("rows276 keys:", z.files[:20])
except Exception as e: print("rows276:", e)
out = {}
def boot_rho(x, y, B=4000):
    x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]
    r = spearmanr(x, y).correlation; rs = []
    for _ in range(B):
        i = rng.integers(0, len(x), len(x)); rs.append(spearmanr(x[i], y[i]).correlation)
    return r, float(np.nanpercentile(rs, 2.5)), float(np.nanpercentile(rs, 97.5)), int(len(x))
def partial(x, y, zc):
    df = pd.DataFrame(dict(x=x, y=y, z=zc)).dropna(); R = df.rank()
    rx = np.polyfit(R.z, R.x, 1); ry = np.polyfit(R.z, R.y, 1)
    ex = R.x - np.polyval(rx, R.z); ey = R.y - np.polyval(ry, R.z)
    return float(np.corrcoef(ex, ey)[0, 1]), len(df)
for fold in ["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"]:
    cl, pc = f"C_lolo_{fold}", f"pds_{fold}"
    if cl not in c or pc not in c: continue
    r, lo, hi, n = boot_rho(c[cl], c[pc]); res = {"overall": dict(rho=r, lo=lo, hi=hi, n=n)}
    for mname, mcol in [("gwcd4i_de_count", "mag_de"), ("split_half_reliability", "half_A")]:
        if mcol not in c or c[mcol].isna().all(): continue
        q = pd.qcut(c[mcol].rank(method="first"), 3, labels=["low", "mid", "high"])
        res[mname] = {}
        for lab in ["low", "mid", "high"]:
            m = (q == lab).to_numpy(); rr, l2, h2, nn = boot_rho(c.loc[m, cl], c.loc[m, pc]); res[mname][lab] = dict(rho=rr, lo=l2, hi=h2, n=nn)
        pr, pn = partial(c[cl], c[pc], c[mcol]); res[mname]["partial_rho"] = dict(rho=pr, n=pn)
        res[mname]["rho_mag_vs_conservation"] = float(spearmanr(c[mcol], c[cl], nan_policy="omit").correlation)
        res[mname]["rho_mag_vs_pds"] = float(spearmanr(c[mcol], c[pc], nan_policy="omit").correlation)
    out[fold] = res
    print(f"\n== {fold}: overall rho(C_lolo, pds) = {r:+.3f} [{lo:+.3f}, {hi:+.3f}] n={n}")
    for mname in [k for k in res if k != "overall"]:
        rr = res[mname]; print(f"   by {mname}: " + " | ".join(f"{lab} {rr[lab]['rho']:+.3f} [{rr[lab]['lo']:+.3f},{rr[lab]['hi']:+.3f}] n={rr[lab]['n']}" for lab in ['low','mid','high']) + f" | partial {rr['partial_rho']['rho']:+.3f} (n={rr['partial_rho']['n']}); rho(mag,C)={rr['rho_mag_vs_conservation']:+.3f}, rho(mag,pds)={rr['rho_mag_vs_pds']:+.3f}")
json.dump(out, open("experiments/288-f8-magnitude-check/work/check288.json", "w"), indent=1)

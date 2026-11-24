"""262 -- three follow-up reads on the saved arrays (prints; appended to RESULTS.md by hand).
1. the own-cosine -> pds conversion curve (cancer folds pooled);
2. build-time-OBSERVABLE selectors vs the below-chance group (mean pds of the selected bin);
3. does CD4 / the union help or hurt on rows where it disagrees with the base?  (BASE vs CD4ONLY vs cd4svd5 on
   those rows); union carriers inside the marginal band."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
ARM = "cd4svd5"
out = {}

import importlib.util  # noqa: E402
spec = importlib.util.spec_from_file_location("an", HERE / "analyse262.py")
# re-use analyse262's TAB (per-fold covariate tables incl. cross-line agreement) without re-printing
import io, contextlib  # noqa: E402
an = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(an)
TAB = an.TAB


def pooled(key):
    return np.concatenate([TAB[f][1][key] for f in CANCER])


p = pooled(f"{ARM}|pds"); own = pooled(f"{ARM}|own"); bc = pooled("cov|base_cos_truth")
dd = pooled("agree_donor_donor"); hv = pooled("agree_held_vs_donors")
cb = pooled("cov|cd4_cos_base"); rb = pooled("cov|resid_cos_base"); ub = pooled("cov|union_count")
pb = pooled("BASE|pds"); pc = pooled("CD4ONLY|pds"); pu = pooled("UNIONONLY|pds"); ph = pooled("h1union|pds")
n = len(p)

print("1. CONVERSION CURVE: own cosine (delivered prediction vs its truth) -> mean pds, cancer folds pooled (n=%d)" % n)
edges = [-1, -0.02, -0.01, 0.0, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 1]
curve = []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (own > lo) & (own <= hi)
    if m.sum():
        curve.append(dict(bin=f"({lo:g}, {hi:g}]", n=int(m.sum()), share=float(m.mean()), mean_pds=float(p[m].mean()),
                          deficit_share=float((1 - p[m]).sum() / (1 - p).sum())))
        print(f"   own cos in ({lo:g}, {hi:g}]: n {m.sum():4d} ({100*m.mean():4.1f} %)  mean pds {p[m].mean():.3f}  deficit share {100*(1-p[m]).sum()/(1-p).sum():5.1f} %")
out["conversion_curve"] = curve
# slope: d pds / d cos in the band (0, 0.05]
band = (own > 0) & (own <= 0.05)
A = np.vstack([own[band], np.ones(band.sum())]).T
slope = float(np.linalg.lstsq(A, p[band], rcond=None)[0][0])
print(f"   least-squares slope of pds on own cosine inside (0, 0.05]: {slope:.1f} pds per unit cosine (= {slope/100:.3f} per +0.01)")
out["slope_band"] = slope

print("\n2. OBSERVABLE SELECTORS: mean pds of the selected bin (a bin below 0.5 would gain from a neutral row; none does)")
sel = {
    "donor-donor agree <= 0": dd <= 0,
    "donor-donor agree <= -0.02": dd <= -0.02,
    "CD4·base cos <= 0 (CD4 complete rows)": (cb <= 0),
    "residual·base cos <= 0": rb <= 0,
    "donor-donor <= 0 AND CD4·base <= 0": (dd <= 0) & (cb <= 0),
    "donor-donor <= 0 AND CD4·base <= 0 AND union absent": (dd <= 0) & (cb <= 0) & (ub == 0),
    "bottom 5 % by donor-donor agreement": dd <= np.nanquantile(dd, 0.05),
    "bottom 5 % by CD4·base cos": cb <= np.nanquantile(cb, 0.05),
    "bottom 5 % by min(donor-donor, CD4·base)": np.fmin(dd, cb) <= np.nanquantile(np.fmin(dd, cb), 0.05),
    "ORACLE: own cosine <= 0 (unobservable)": own <= 0,
    "ORACLE: base·truth <= 0 (unobservable)": bc <= 0,
}
out["selectors"] = {}
for k, m in sel.items():
    m = m & np.isfinite(p)
    out["selectors"][k] = dict(n=int(m.sum()), share=float(m.mean()), mean_pds=float(p[m].mean()), deficit_share=float((1 - p[m]).sum() / (1 - p).sum()))
    print(f"   {k:58s} n {m.sum():4d} ({100*m.mean():4.1f} %)  mean pds {p[m].mean():.3f}  deficit share {100*(1-p[m]).sum()/(1-p).sum():5.1f} %")

print("\n3. DOES THE RESIDUAL HELP OR HURT WHERE IT DISAGREES WITH THE BASE?  mean pds by object on the selected rows")
groups = {
    "all rows": np.ones(n, bool),
    "CD4·base cos <= 0": cb <= 0,
    "CD4·base cos in (0, 0.02]": (cb > 0) & (cb <= 0.02),
    "CD4·base cos > 0.02": cb > 0.02,
    "donor-donor agree <= 0": dd <= 0,
    "donor-donor agree > 0.05": dd > 0.05,
    "base·truth <= 0 (oracle group e)": bc <= 0,
    "base·truth in (0, 0.05] (the marginal band)": (bc > 0) & (bc <= 0.05),
    "base·truth > 0.05": bc > 0.05,
    "marginal band AND union carried": (bc > 0) & (bc <= 0.05) & (ub > 0),
    "marginal band AND union absent": (bc > 0) & (bc <= 0.05) & (ub == 0),
    "group e AND union carried": (bc <= 0) & (ub > 0),
    "group e AND union absent": (bc <= 0) & (ub == 0),
}
out["objects_by_group"] = {}
print(f"   {'rows':48s} {'n':>4s}  {'BASE':>6s} {'CD4ONLY':>8s} {'UNIONONLY':>9s} {'h1union':>8s} {'cd4svd5':>8s}   cd4svd5-BASE")
for k, m in groups.items():
    m = m & np.isfinite(p)
    r = dict(n=int(m.sum()), BASE=float(pb[m].mean()), CD4ONLY=float(pc[m].mean()), UNIONONLY=float(pu[m].mean()),
             h1union=float(ph[m].mean()), cd4svd5=float(p[m].mean()))
    out["objects_by_group"][k] = r
    print(f"   {k:48s} {m.sum():4d}  {r['BASE']:.3f}   {r['CD4ONLY']:.3f}    {r['UNIONONLY']:.3f}    {r['h1union']:.3f}    {r['cd4svd5']:.3f}   {r['cd4svd5']-r['BASE']:+.3f}")

# cross-tab of the deficit: held-vs-donor x donor-donor
print("\n4. CROSS-TAB of the pooled deficit: held-vs-donor agreement (rows) x donor-donor agreement (cols); cell = deficit share / n / mean pds")
cuts = [(-1, 0.0), (0.0, 0.05), (0.05, 1)]
xt = {}
for lo, hi in cuts:
    row = []
    for lo2, hi2 in cuts:
        m = (hv > lo) & (hv <= hi) & (dd > lo2) & (dd <= hi2)
        row.append(f"{100*(1-p[m]).sum()/(1-p).sum():5.1f} % / n {m.sum():3d} / {p[m].mean() if m.sum() else float('nan'):.2f}")
        xt[f"held({lo},{hi}]|donor({lo2},{hi2}]"] = dict(n=int(m.sum()), deficit_share=float((1 - p[m]).sum() / (1 - p).sum()),
                                                      mean_pds=float(p[m].mean()) if m.sum() else None)
    print(f"   held-vs-donor ({lo:g}, {hi:g}]:  " + "  |  ".join(row))
out["crosstab"] = xt

# the top decile's covariates vs the bottom's (pooled medians)
print("\n5. TOP (pds = 1) vs BOTTOM decile, pooled medians")
for k in ("cov|base_cos_truth", "agree_held_vs_donors", "agree_donor_donor", "cov|resid_cos_truth", "cov|cd4_cos_base", f"{ARM}|truth_norm",
          f"{ARM}|truth_nn", f"{ARM}|pred_nn", "cov|cd4_norm", "cov|union_count"):
    v = pooled(k)
    print(f"   {k:28s} top-block (pds=1) median {np.nanmedian(v[p == 1]):+.3f}   bottom-decile median {np.nanmedian(v[p <= np.quantile(p, 0.1)]):+.3f}   rest {np.nanmedian(v[(p < 1) & (p > np.quantile(p, 0.1))]):+.3f}")

(WORK / "extra262.json").write_text(json.dumps(out, indent=1, default=float))

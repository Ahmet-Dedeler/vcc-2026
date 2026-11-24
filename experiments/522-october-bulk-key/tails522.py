"""522 (4a) — where the percentile cap's October tail comes from, and whether the held-out truth backs it.

Tail T (per fold) = the targets the percentile cap (P_PROD) weights above 6: carried, key/median > 6 (only when hi > 6).
Per tail target:
  terms   the two LODO terms v_j / median (each fold carries exactly its two sources' terms); "both" = min term >= 3,
          "one donor" = min term < 1 (the other donor sits below the median target), "mixed" otherwise
  donors  raw donor-row cosine cos(d_1, d_2) (non-panel), each donor's row norm / its median over carried targets
  cover   each source's robust strength count and per-target noise (median |z| / 0.6745; few cells -> high noise),
          noise ratio = noise / the source's median noise over the fold's targets ("low coverage" = max ratio > 1.5)
  ess     DepMap mean Chronos over 17 lines (< -0.5 = essential-like), K562 essential list membership
  truth   held-out truth energy d / median d over carried; the local ideal-scale ratio u = (x/a) / (X/A) on the cap-6 arm
          (x = <P, T>, a = |P|^2 summed over frames; u > 1: the truth wants MORE weight than 6) and on the pct arm
          (u < 1: over-weighted)
Exact decomposition (fold522, |d| 0 route): arms 'cap 6 with class C raised to its pct weight' for each class, vs cap 6.
A/B/C folds: the same where the dump allows (no per-donor terms in 518's dumps).

    uv run python experiments/522-october-bulk-key/tails522.py        (~1-2 GB; after arms522 on every surface)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import fold522 as FD  # noqa: E402

WORK = HERE / "work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")


def spearman(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 4:
        return np.nan
    ra, rb = np.argsort(np.argsort(a[ok])), np.argsort(np.argsort(b[ok]))
    return float(np.corrcoef(ra, rb)[0, 1])


def ideal_u(arr):
    x, a = arr["x"].sum(0), arr["a"].sum(0)
    X, A = x.sum(), a.sum()
    return (x / np.maximum(a, 1e-300)) / (X / A), x / X - a / A


def run(panel, held, st):
    t0 = time.time()
    S = FD.load(panel, held, check=True)
    az = np.load(WORK / f"arms522-{panel}-{held}.npz")
    car = S.carried
    med = float(np.median(S.kv[car]))
    kn = np.where(car, S.kv / med, np.nan)
    hp = FD.pct_hi(S.kv, car, FD.P_PROD)
    w6, wp = az["w__cap6"], az["w__pctPROD"]
    T = car & (kn > 6.0) & (wp > 6.0 + 1e-12)
    panel_t = [str(t) for t in st[f"{panel}__targets"]]
    idx = np.asarray([panel_t.index(t) for t in S.targets])
    g = lambda k: np.asarray(st[f"{panel}__{k}"], float)[idx] if f"{panel}__{k}" in st else np.full(S.n, np.nan)
    dep = g("dep_mean")
    ess = np.asarray(st[f"{panel}__k562_essential"])[idx]
    cnt = np.vstack([g(f"{s}__count") for s in S.sources])
    noise = np.vstack([g(f"{s}__noise") for s in S.sources])
    with np.errstate(all="ignore"):
        nrat = noise / np.nanmedian(noise, axis=1, keepdims=True)
        lowcov = np.nanmax(nrat, 0) > 1.5
    arr6 = {k: az[f"cap6__{k}"] for k in ("x", "a", "d")}
    arrp = {k: az[f"pctPROD__{k}"] for k in ("x", "a", "d")}
    u6, g6 = ideal_u(arr6)
    up, gp = ideal_u(arrp)
    d = arr6["d"]
    drel = d / np.median(d[car])
    rec = dict(panel=panel, fold=held, n=S.n, n_carried=int(car.sum()), hi=hp, n_tail=int(T.sum()),
               frac_floor=float(np.mean(kn[car] <= 0.5)), frac_key_le0=float(np.mean(S.kv[car] <= 0)),
               key_median=med, sp_kn_d=spearman(kn[car], drel[car]), sp_kn_u6=spearman(kn[car], (u6 * w6)[car]),
               tail={}, classes={}, targets=[])
    if S.kind == "oct":
        Vn = S.V / med
        vmin, vmax = np.nanmin(Vn, 0), np.nanmax(Vn, 0)
        cls = np.where(vmin >= 3.0, "both", np.where(vmin < 1.0, "one donor", "mixed"))
        dn = S.dnorm / np.nanmedian(np.where(car[None, :], S.dnorm, np.nan), axis=1, keepdims=True)
    else:
        Vn = vmin = vmax = None
        cls = np.full(S.n, "n/a")
        dn = None
    essl = (dep < -0.5) | ess
    for i in np.flatnonzero(T):
        rec["targets"].append(dict(
            target=S.targets[i], kn=float(kn[i]), w_pct=float(wp[i]),
            terms=None if Vn is None else [float(v) for v in Vn[:, i]], cls=str(cls[i]),
            cos_dd=None if S.kind != "oct" else float(S.cos_dd[i]),
            dnorm_rel=None if dn is None else [float(v) for v in dn[:, i]],
            count=[float(v) for v in cnt[:, i]], noise_rel=[float(v) for v in nrat[:, i]], lowcov=bool(lowcov[i]),
            dep_mean=float(dep[i]), essential_like=bool(essl[i]), d_rel=float(drel[i]), u_cap6=float(u6[i]),
            u_pct=float(up[i]), grad_cap6=float(g6[i])))
    nt = car & ~T
    rec["tail"] = dict(
        median_kn=float(np.median(kn[T])) if T.any() else np.nan,
        median_d_rel_tail=float(np.median(drel[T])) if T.any() else np.nan,
        median_d_rel_nontail=float(np.median(drel[nt])),
        median_u6_tail=float(np.median(u6[T])) if T.any() else np.nan,
        frac_u6_gt1_tail=float(np.mean(u6[T] > 1)) if T.any() else np.nan,
        median_u6_nontail_top=float(np.median(u6[car & (kn > 3) & ~T])) if (car & (kn > 3) & ~T).any() else np.nan,
        median_upct_tail=float(np.median(up[T])) if T.any() else np.nan,
        frac_upct_lt1_tail=float(np.mean(up[T] < 1)) if T.any() else np.nan,
        first_order_gain_tail_cap6=float(np.sum(g6[T])),
        sp_kn_u6w6_tail=spearman(kn[T], (u6 * w6)[T]) if T.sum() >= 4 else np.nan,
        ideal_w_tail_median=float(np.median((u6 * w6)[T])) if T.any() else np.nan,
        frac_essential_like=float(np.mean(essl[T])) if T.any() else np.nan,
        frac_lowcov=float(np.mean(lowcov[T])) if T.any() else np.nan,
        cls_counts={c: int(np.sum(T & (cls == c))) for c in np.unique(cls[T])} if T.any() else {})
    # ---- exact decomposition: cap 6, then raise one class to its pct weight
    s6, _ = S.score(w6)
    sp, _ = S.score(wp)
    rec["r_cap6_vs_cap3"] = None
    rec["dcos_pct_vs_cap6"] = sp["cos_agg"] / s6["cos_agg"] - 1
    masks = {"all tail": T}
    if S.kind == "oct":
        for c in ("both", "mixed", "one donor"):
            masks[f"terms: {c}"] = T & (cls == c)
        masks["donor rows agree (cos_dd >= median of tail)"] = T & (S.cos_dd >= np.nanmedian(S.cos_dd[T])) if T.any() else T
    masks["essential-like"] = T & essl
    masks["not essential-like"] = T & ~essl
    masks["low coverage"] = T & lowcov
    masks["not low coverage"] = T & ~lowcov
    masks["truth-backed (u6 > 1)"] = T & (u6 > 1)          # an ORACLE split, for the anatomy only
    masks["truth-against (u6 <= 1)"] = T & (u6 <= 1)
    for nm, m in masks.items():
        if not m.any():
            rec["classes"][nm] = dict(n=0, r=0.0)
            continue
        w = np.where(m, wp, w6)
        s, _ = S.score(w)
        rec["classes"][nm] = dict(n=int(m.sum()), r=s["cos_agg"] / s6["cos_agg"] - 1, dpds=s["pds"] - s6["pds"])
    assert abs(rec["classes"]["all tail"]["r"] - rec["dcos_pct_vs_cap6"]) < 1e-12 if T.any() else True
    rec["seconds"] = time.time() - t0
    return rec


def main():
    st = np.load(WORK / "strength522.npz")
    out = {}
    for p, f in FD.surfaces(("abc", "efp", "efr", "real")):
        if not (WORK / f"arms522-{p}-{f}.npz").exists():
            continue
        r = run(p, f, st)
        out[f"{p}/{f}"] = r
        t = r["tail"]
        print(f"[{p}/{f}] hi {r['hi']:.2f} tail {r['n_tail']} ({t['cls_counts']}); floor {r['frac_floor']:.2f}; "
              f"d_rel tail {t['median_d_rel_tail']:.2f} vs non-tail {t['median_d_rel_nontail']:.2f}; u6 tail {t['median_u6_tail']:.2f} "
              f"(>1: {t['frac_u6_gt1_tail']:.2f}); ideal w tail {t['ideal_w_tail_median']:.2f}; ess-like {t['frac_essential_like']:.2f}; "
              f"lowcov {t['frac_lowcov']:.2f}; pct vs cap6 {100 * r['dcos_pct_vs_cap6']:+.2f} %", flush=True)
        for nm, c in r["classes"].items():
            print(f"      {nm:44s} n {c['n']:3d}  r vs cap6 {100 * c['r']:+6.2f} %", flush=True)
    (WORK / "tails522.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()

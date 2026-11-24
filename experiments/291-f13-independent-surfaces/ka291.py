"""291 step 0 — known-answer rows KA-1 (290's primary cell rebuilt from the per-source store),
KA-2 (279's c_held, read) and KA-5 (the half-1 aggregation reproduces feat290.csv).

    PYTHONDONTWRITEBYTECODE=1 uv run python ka291.py     -> work/ka291.json
"""
from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd

import lib291 as L

T0 = time.time()
LOG: list[str] = []


def P(m=""):
    line = f"[{time.time() - T0:6.1f}s] {m}" if m else ""
    print(line, flush=True)
    LOG.append(line)


def main():
    h = L.check_prereg()
    P(f"prereg {h} MATCH")
    st = L.Store()
    P(f"store: {len(st.genes)} genes, {len(st.keys)} sources, {len(set(st.lines))} lines")

    # ---- KA-5: the half-1 / half-2 aggregation reproduces feat290.csv -------------------------
    ref = pd.read_csv(L.FEAT290)
    assert list(ref["gene"]) == st.genes, "gene order drifted"
    ka5 = {}
    for tag, rows in (("h1", st.half1), ("h2", st.half2)):
        agg = st.aggregate(rows)
        d, npat = 0.0, True
        for c in L.FEATS290 + L.EXTRA_COLS:
            a, b = agg[c].to_numpy(float), ref[f"{c}_{tag}"].to_numpy(float)
            npat &= bool(np.array_equal(np.isfinite(a), np.isfinite(b)))
            m = np.isfinite(a) & np.isfinite(b)
            d = max(d, float(np.max(np.abs(a[m] - b[m]))) if m.any() else 0.0)
        ka5[tag] = dict(max_abs_diff=d, nan_pattern_identical=npat, pass_=bool(d < 1e-6 and npat))
        P(f"KA-5 {tag}: max|diff| vs feat290.csv {d:.2e}, NaN pattern identical {npat} "
          f"-> {'PASS' if ka5[tag]['pass_'] else 'FAIL'}")
    if not all(v["pass_"] for v in ka5.values()):
        raise SystemExit("KA-5 FAILED — abort (AGENTS.md:1375-1382)")

    # ---- KA-2: 279's c_held, read (not recomputed) --------------------------------------------
    z = np.load(L.HELD279, allow_pickle=True)
    c_held = np.asarray(z["c_held"], float)
    znf = pd.read_csv(L.ZNF279, low_memory=False)
    assert list(znf["gene"]) == st.genes, "znf279 gene order drifted"
    d_znf = float(np.nanmax(np.abs(c_held - znf["C_held"].to_numpy(float))))
    ka2 = dict(n=int(np.isfinite(c_held).sum()), median=float(np.nanmedian(c_held)),
               max_abs_diff_vs_znf279=d_znf,
               pass_=bool(int(np.isfinite(c_held).sum()) == 17106
                          and abs(float(np.nanmedian(c_held)) - 0.018605) < 5e-7
                          and d_znf < 1e-6))
    P(f"KA-2 c_held read: n {ka2['n']}, median {ka2['median']:+.6f}, "
      f"max|diff| vs znf279.C_held {d_znf:.1e} -> {'PASS' if ka2['pass_'] else 'FAIL'}")
    if not ka2["pass_"]:
        raise SystemExit("KA-2 FAILED — abort")

    # ---- KA-1: 290's primary cell rebuilt from the store ---------------------------------------
    rem = np.load(L.REM284, allow_pickle=True)
    ok = np.asarray(rem["ok"], bool)
    fold = np.asarray(rem["fold"], int)
    r2 = np.asarray(rem["r2"], float)
    assert [str(g) for g in rem["gene"]] == st.genes
    agg1 = st.aggregate(st.half1)
    X = agg1[L.FEATS290].to_numpy(float)[ok]
    y = r2[ok]
    f = fold[ok]
    P(f"KA-1 unit n {int(ok.sum())}, folds {np.bincount(f).tolist()}")
    pred, info = L.oof_ridge(X, y, f)
    rho = L.spearman(y, pred)
    ka1 = dict(n=int(ok.sum()), folds=np.bincount(f).tolist(), rho=rho,
               alphas={str(k): v["alpha"] for k, v in info.items()},
               required=0.317385, diff=abs(rho - 0.317385),
               pass_=bool(round(rho, 6) == 0.317385))
    P(f"KA-1 290's primary cell: rho {rho:+.6f} (required +0.317385, diff {ka1['diff']:.1e}) "
      f"-> {'PASS' if ka1['pass_'] else 'FAIL'}")
    if not ka1["pass_"]:
        raise SystemExit("KA-1 FAILED — abort (AGENTS.md:1375-1382)")

    out = dict(prereg_sha256=h, ka1=ka1, ka2=ka2, ka5=ka5, seconds=time.time() - T0)
    (L.WORK / "ka291.json").write_text(json.dumps(out, indent=1, default=float))
    (L.WORK / "ka291.log").write_text("\n".join(LOG))
    P("wrote work/ka291.json")


if __name__ == "__main__":
    main()

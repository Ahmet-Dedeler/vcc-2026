"""276 — two robustness reads the frozen definition does not cover, labelled POST-HOC.

 1. `C_out`: conservation over ONLY the sources that are NOT part of the incumbent's prediction
    (the incumbent uses k562_gwps + perturbmatch, orion_hct116, orion_hek293t, the three CD4
    conditions and the two H1 banks).  If `C_lolo` predicts transfer only because the
    prediction's own components agree with each other, `C_out` must fall to nothing.
 2. the mechanism: does conservation predict the base·truth cosine (262's `cov|base_cos_truth`),
    the quantity `pds` is a threshold in (F4)?

    PYTHONDONTWRITEBYTECODE=1 uv run python extra276.py > work/extra276.md
"""
from __future__ import annotations

import json

import numpy as np
from scipy.stats import rankdata, spearmanr

import lib276 as L

SEED = 20260906
B = 2000
E262 = L.ROOT / "experiments/262-where-pds-is-lost/work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
HELD = {"k562_gwps": "K562", "orion_hct116": "HCT116", "orion_hek293t": "HEK293T"}
IN_PREDICTION = ("k562_gwps", "perturbmatch", "orion_hct116", "orion_hek293t",
                 "cd4_rest", "cd4_stim8hr", "cd4_stim48hr", "gse295214_h1", "h1_2025")


def boot(x, y, b=B, seed=SEED):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[ok], np.asarray(y)[ok]
    rng = np.random.default_rng(seed)
    d = np.array([spearmanr(x[s := rng.integers(0, x.size, x.size)], y[s]).statistic
                  for _ in range(b)])
    return float(spearmanr(x, y).statistic), [float(v) for v in np.percentile(d, [2.5, 97.5])], int(x.size)


def main() -> None:
    zc = np.load(L.WORK / "cons276.npz", allow_pickle=True)
    keys = [str(k) for k in zc["keys"]]
    C = np.asarray(zc["Ccross"], dtype=np.float64)
    pairs = [tuple(p) for p in np.asarray(zc["pairs_cross"])]
    idx = {t: i for i, t in enumerate(L.PANEL)}

    def score(drop):
        d = {keys.index(k) for k in drop if k in keys}
        keep = [p for p, (i, j) in enumerate(pairs) if i not in d and j not in d]
        with np.errstate(invalid="ignore"):
            return np.nanmean(C[keep], axis=0), np.isfinite(C[keep]).sum(axis=0)

    variants = {
        "C_lolo (held line out)": lambda line: L.HELD_LINE_SOURCES.get(line, ()),
        "C_out (held line out AND every source the incumbent uses)":
            lambda line: tuple(set(L.HELD_LINE_SOURCES.get(line, ())) | set(IN_PREDICTION)),
        "C_out, no CD4 pairs only":
            lambda line: tuple(set(L.HELD_LINE_SOURCES.get(line, ()))
                               | {"cd4_rest", "cd4_stim8hr", "cd4_stim48hr"}),
    }
    print("## POST-HOC 1 — is it the prediction's own self-consistency?\n")
    print("| conservation score | sources in the map | median pairs / target | "
          "rho with `pds` (pooled, 3 cancer folds) | rho with base·truth cosine |")
    print("|---|---:|---:|---:|---:|")
    for name, fn in variants.items():
        px, py, pb, npair = [], [], [], []
        nsrc = None
        for fold in CANCER:
            z = np.load(E262 / f"diag262-{fold}.npz", allow_pickle=True)
            tg = [str(t) for t in z["targets"]]
            sel = np.asarray([idx[t] for t in tg])
            drop = fn(HELD[fold])
            nsrc = len(keys) - len({k for k in drop if k in keys})
            s, n = score(drop)
            s, n = s[sel], n[sel]
            pds = np.asarray(z["cd4svd5|pds"], dtype=np.float64)
            base = np.asarray(z["cov|base_cos_truth"], dtype=np.float64)
            ok = np.isfinite(s)
            r = lambda v: rankdata(v[ok]) / ok.sum()
            px.append(r(s)); py.append(r(pds)); pb.append(r(base)); npair.append(n[ok])
        px, py, pb = np.concatenate(px), np.concatenate(py), np.concatenate(pb)
        rp, cip, _ = boot(px, py)
        rb, cib, _ = boot(px, pb)
        print(f"| {name} | {nsrc} | {np.median(np.concatenate(npair)):.0f} | "
              f"{rp:+.3f} [{cip[0]:+.3f}, {cip[1]:+.3f}] | {rb:+.3f} [{cib[0]:+.3f}, {cib[1]:+.3f}] |")

    print("\n## POST-HOC 2 — the mechanism, by band of `C_out` (no incumbent source, held line out)\n")
    px, py, pa = [], [], []
    for fold in CANCER:
        z = np.load(E262 / f"diag262-{fold}.npz", allow_pickle=True)
        tg = [str(t) for t in z["targets"]]
        sel = np.asarray([idx[t] for t in tg])
        s, _ = score(tuple(set(L.HELD_LINE_SOURCES.get(HELD[fold], ())) | set(IN_PREDICTION)))
        s = s[sel]
        ok = np.isfinite(s)
        px.append(s[ok])
        py.append(np.asarray(z["cd4svd5|pds"], dtype=np.float64)[ok])
        pa.append((np.asarray(z["cov|base_cos_truth"], dtype=np.float64) <= 0)[ok])
    c, p, a = np.concatenate(px), np.concatenate(py), np.concatenate(pa)
    qs = np.percentile(c, [10, 30, 50, 70, 90])
    print("| `C_out` band | n | mean `pds` | share anti-aligned | median `C_out` |")
    print("|---|---:|---:|---:|---:|")
    for lo, hi, nm in ((-np.inf, qs[0], "bottom decile"), (qs[0], qs[1], "10–30 %"),
                       (qs[1], qs[2], "30–50 %"), (qs[2], qs[3], "50–70 %"),
                       (qs[3], qs[4], "70–90 %"), (qs[4], np.inf, "top decile")):
        m = (c > lo) & (c <= hi)
        print(f"| {nm} | {int(m.sum())} | {p[m].mean():.3f} | {a[m].mean()*100:.1f} % "
              f"| {np.median(c[m]):+.4f} |")
    print(f"| all | {c.size} | {p.mean():.3f} | {a.mean()*100:.1f} % | {np.median(c):+.4f} |")


if __name__ == "__main__":
    main()

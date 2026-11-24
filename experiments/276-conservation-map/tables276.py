"""276 — the reading tables: KA-3's sensitivity, the per-fold numbers, and the interpretable
tercile / decile view of what a conservation score buys.

    PYTHONDONTWRITEBYTECODE=1 uv run python tables276.py > work/tables276.md
"""
from __future__ import annotations

import json

import numpy as np
from scipy.stats import rankdata, spearmanr

import lib276 as L

E262 = L.ROOT / "experiments/262-where-pds-is-lost/work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
HELD = {"k562_gwps": "K562", "orion_hct116": "HCT116", "orion_hek293t": "HEK293T"}
DD = {"k562_gwps": ("orion_hct116", "orion_hek293t"),
      "orion_hct116": ("k562_gwps", "orion_hek293t"),
      "orion_hek293t": ("k562_gwps", "orion_hct116")}


def main() -> None:
    zc = np.load(L.WORK / "cons276.npz", allow_pickle=True)
    keys = [str(k) for k in zc["keys"]]
    Ccross = np.asarray(zc["Ccross"], dtype=np.float64)
    pairs = [tuple(p) for p in np.asarray(zc["pairs_cross"])]
    c_raw = np.asarray(zc["c_raw"])
    idx = {t: i for i, t in enumerate(L.PANEL)}

    # ---------- KA-3 sensitivity: the donor-donor baseline on three gene axes ---------------
    print("## KA-3 sensitivity — 262's `agree_donor_donor` recomputed on three gene axes\n")
    print("| axis | genes | rho(donor-donor, pds), pooled | folds pooled |")
    print("|---|---:|---:|---|")
    axes = {
        "276's common axis (26 sources)": [str(g) for g in np.load(L.WORK / "rows276.npz",
                                                                   allow_pickle=True)["axis"]],
        "the three donors' own shared axis": L.common_axis(["k562_gwps", "orion_hct116",
                                                            "orion_hek293t"]),
        "the two Orion lines' shared axis": L.common_axis(["orion_hct116", "orion_hek293t"]),
    }
    for name, ax in axes.items():
        rows, folds = {}, list(CANCER)
        axset = set(ax)
        for k in ("k562_gwps", "orion_hct116", "orion_hek293t"):
            _, g, _ = L.read_bank(L.SOURCES[k][0])
            if axset <= set(g):
                rows[k] = L.load_source(k, ax)
        folds = [f for f in CANCER if all(s in rows for s in DD[f])]
        px, py = [], []
        for fold in folds:
            z = np.load(E262 / f"diag262-{fold}.npz", allow_pickle=True)
            tg = [str(t) for t in z["targets"]]
            sel = np.asarray([idx[t] for t in tg])
            a, b = DD[fold]
            m = rows[a][1] & rows[b][1]
            v = np.full(len(L.PANEL), np.nan)
            v[m] = L.cosine_rows(rows[a][0][m].astype(np.float64), rows[b][0][m].astype(np.float64))
            v = v[sel]
            pds = np.asarray(z["cd4svd5|pds"], dtype=np.float64)
            ok = np.isfinite(v)
            rx = np.full(v.size, np.nan); ry = np.full(v.size, np.nan)
            rx[ok] = rankdata(v[ok]) / ok.sum(); ry[ok] = rankdata(pds[ok]) / ok.sum()
            px.append(rx[ok]); py.append(ry[ok])
        r = spearmanr(np.concatenate(px), np.concatenate(py)).statistic
        print(f"| {name} | {len(ax)} | {r:+.3f} | {", ".join(folds)} |")
    print("\n262 reports +0.381 pooled (R262 §4.1, `analysis262.json`).\n")

    # ---------- per-fold table --------------------------------------------------------------
    d = json.loads((L.WORK / "predict276.json").read_text())
    print("## Per-fold Spearman with `cd4svd5`'s per-target `pds`\n")
    print("| predictor | K562 fold | HCT116 fold | HEK293T fold | kolf2 fold | pooled (3 cancer) |")
    print("|---|---:|---:|---:|---:|---:|")
    for row in d["table"]:
        cells = []
        for f in ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"):
            v = row[f"rho|{f}"]
            cells.append(f"{v['rho']:+.3f}" if v["rho"] is not None else "—")
        p = row["rho|pooled3"]
        print(f"| {row['predictor']} | " + " | ".join(cells)
              + f" | **{p['rho']:+.3f}** [{p['ci'][0]:+.3f}, {p['ci'][1]:+.3f}] |")

    # ---------- tercile / decile view -------------------------------------------------------
    print("\n## What a conservation score buys, read as groups (LOLO, pooled over 3 cancer folds)\n")
    zp = np.load(L.WORK / "predict276.npz", allow_pickle=True)
    allc, allp, alla = [], [], []
    for fold in CANCER:
        allc.append(np.asarray(zp[f"{fold}|c_lolo"]))
        allp.append(np.asarray(zp[f"{fold}|pds"]))
        alla.append(np.asarray(zp[f"{fold}|anti"]))
    c = np.concatenate(allc); p = np.concatenate(allp); a = np.concatenate(alla).astype(bool)
    ok = np.isfinite(c)
    c, p, a = c[ok], p[ok], a[ok]
    qs = np.percentile(c, [10, 30, 50, 70, 90])
    bins = [(-np.inf, qs[0], "bottom decile"), (qs[0], qs[1], "10–30 %"),
            (qs[1], qs[2], "30–50 %"), (qs[2], qs[3], "50–70 %"),
            (qs[3], qs[4], "70–90 %"), (qs[4], np.inf, "top decile")]
    print("| `C_lolo` band | n (target, fold) | mean `pds` | share anti-aligned | median `C_lolo` |")
    print("|---|---:|---:|---:|---:|")
    for lo, hi, name in bins:
        m = (c > lo) & (c <= hi)
        print(f"| {name} | {int(m.sum())} | {p[m].mean():.3f} | {a[m].mean()*100:.1f} % "
              f"| {np.median(c[m]):+.4f} |")
    print(f"| all | {int(c.size)} | {p.mean():.3f} | {a.mean()*100:.1f} % | {np.median(c):+.4f} |")

    # the same for the h1-test instrument
    print("\n## The same on h1-test (85 Flex-400 targets, H1 sources removed)\n")
    hc = np.asarray(zp["h1_c_lolo"]); hp = np.asarray(zp["h1_pds"])
    ok = np.isfinite(hc)
    hc, hp = hc[ok], hp[ok]
    qs = np.percentile(hc, [33.3, 66.7])
    print("| `C_lolo` tercile | n | mean `pds` | median `C_lolo` |")
    print("|---|---:|---:|---:|")
    for lo, hi, name in ((-np.inf, qs[0], "bottom"), (qs[0], qs[1], "middle"), (qs[1], np.inf, "top")):
        m = (hc > lo) & (hc <= hi)
        print(f"| {name} | {int(m.sum())} | {hp[m].mean():.3f} | {np.median(hc[m]):+.4f} |")

    # ---------- distribution of the statistic ----------------------------------------------
    print("\n## The distribution of `C_raw` over the 300 panel targets\n")
    b = json.loads((L.WORK / "build276.json").read_text())
    r = json.loads((L.WORK / "rel276.json").read_text())
    print(f"- cross-line, raw: mean {b['c_raw']['mean']:+.4f}, median {b['c_raw']['median']:+.4f}, "
          f"p10 {b['c_raw']['p10']:+.4f}, p90 {b['c_raw']['p90']:+.4f}, "
          f"share <= 0 {b['c_raw']['share_le0']*100:.1f} %")
    print(f"- same-line pairs (reproducibility, not conservation): mean "
          f"{b['c_within']['mean']:+.4f} on {b['c_within']['n']} targets, {b['n_within_pairs']} pairs")
    print(f"- target-shuffled null: mean {b['null']['mean']:+.5f}, sd over targets "
          f"{b['null']['sd_over_targets']:.5f}, mean 95 % band "
          f"[{b['null']['mean_band'][0]:+.4f}, {b['null']['mean_band'][1]:+.4f}]")
    print(f"- reliability-corrected (4 sources with measured halves): median "
          f"{r['c_rel']['median']:+.3f}, p10 {r['c_rel']['p10']:+.3f}, p90 {r['c_rel']['p90']:+.3f}, "
          f"share <= 0 {r['c_rel']['share_le0']*100:.1f} %, n {r['c_rel']['n']}")

    # per-pair medians, the cross-line agreement matrix
    print("\n## Median cross-line pair cosine, by pair (the 20 pairs with the most targets)\n")
    n_ok = np.isfinite(Ccross).sum(axis=1)
    order = np.argsort(-n_ok)[:20]
    print("| pair | targets | median cosine |")
    print("|---|---:|---:|")
    for p_ in order:
        i, j = pairs[p_]
        v = Ccross[p_][np.isfinite(Ccross[p_])]
        print(f"| {keys[i]} ({L.SOURCES[keys[i]][1]}) — {keys[j]} ({L.SOURCES[keys[j]][1]}) "
              f"| {v.size} | {np.median(v):+.4f} |")


if __name__ == "__main__":
    main()

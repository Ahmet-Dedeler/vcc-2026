"""430 analysis: per fold and averaged over the three cancer folds, sign accuracy (pooled; per-target mean over targets
with >= 10 gate genes) and the L1 nmae proxy, each as a difference vs PP (the production composition), with a paired
bootstrap over targets (2,000 draws; the cancer mean resamples each fold independently and averages the fold deltas).

Target support: targets covered by the fold's PP sources AND by the cross-study donor set (so M3S/PPxs/BS* are defined
on every counted target; Orion folds lose the targets K562-GWPS does not cover). The nmae proxy's k was fitted by
fold430 on all covered targets and is held fixed here.

    uv run python analyse430.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FOLDS = (*CANCER, "kolf2")
CANDS = ("PP", "PPprior", "PPg0", "Wonly", "pre", "M3L", "M3Lc", "BLw1", "BLw1.5", "BLw2", "M3S", "M3Sc", "PPxs",
         "BSw1", "BSw1.5", "BSw2", "BScw1.5", "C4L", "C4S", "PPxsg0", "BXw1.5")
GATES = ("z1000", "z200", "top50", "top200", "top1000", "tau0.1", "tau0.2")
NB = 2000


def stats(z, key, m):
    h, n, nm = z[f"{key}__hits"][m], z[f"{key}__n"][m].astype(float), z[f"{key}__nmae"][m]
    return h, n, nm


def fold_boot(z, truth, gate, cand, m, idx, ref="PP"):
    """per-draw (Δpooled, Δpertarget, Δnmae) vs ref; idx (NB, n) resampled target positions within the mask."""
    kb, kc = f"{truth}|{ref}|{gate}", f"{truth}|{cand}|{gate}"
    hb, nb, mb = stats(z, kb, m)
    hc, nc, mc = stats(z, kc, m)
    ok = nb >= 10
    pool = hc[idx].sum(1) / nc[idx].sum(1) - hb[idx].sum(1) / nb[idx].sum(1)
    w = ok[idx].astype(float)
    pt = ((hc / np.maximum(nc, 1) - hb / np.maximum(nb, 1))[idx] * w).sum(1) / np.maximum(w.sum(1), 1)
    fin = np.isfinite(mc) & np.isfinite(mb)
    wf = fin[idx].astype(float)
    dn = (np.nan_to_num(mc - mb)[idx] * wf).sum(1) / np.maximum(wf.sum(1), 1)
    return pool, pt, dn


def point(z, truth, gate, cand, m):
    h, n, nm = stats(z, f"{truth}|{cand}|{gate}", m)
    ok = n >= 10
    return dict(pooled=float(h.sum() / n.sum()), pertarget=float(np.mean(h[ok] / n[ok])) if ok.any() else None,
                nmae=float(np.nanmean(nm)), k=float(z[f"{truth}|{cand}|{gate}__k"]), n_targets=int(m.sum()))


def main():
    rng = np.random.default_rng(430)
    Z = {f: np.load(WORK / f"fold430-{f}.npz") for f in FOLDS}
    masks = {f: np.asarray(Z[f]["cov_xs"], bool) for f in FOLDS}
    idx = {f: rng.integers(0, int(masks[f].sum()), size=(NB, int(masks[f].sum()))) for f in FOLDS}
    out = dict(support={f: int(masks[f].sum()) for f in FOLDS}, folds={}, cancer={})
    for truth in ("raw", "cen"):
        for gate in GATES:
            for cand in CANDS:
                draws = []
                for f in FOLDS:
                    key = f"{truth}|{cand}|{gate}"
                    if f"{key}__hits" not in Z[f].files:
                        continue
                    p = point(Z[f], truth, gate, cand, masks[f])
                    pb = point(Z[f], truth, gate, "PP", masks[f])
                    bo = fold_boot(Z[f], truth, gate, cand, masks[f], idx[f])
                    rec = dict(**p, d_pooled=p["pooled"] - pb["pooled"],
                               d_pertarget=(p["pertarget"] - pb["pertarget"]) if p["pertarget"] is not None else None,
                               d_nmae=p["nmae"] - pb["nmae"],
                               ci_pooled=np.quantile(bo[0], [.025, .975]).tolist(),
                               ci_pertarget=np.quantile(bo[1], [.025, .975]).tolist(),
                               ci_nmae=np.quantile(bo[2], [.025, .975]).tolist())
                    out["folds"][f"{f}|{key}"] = rec
                    if f in CANCER:
                        draws.append((rec, bo))
                if len(draws) == 3:
                    B = [np.mean([b[i] for _, b in draws], 0) for i in range(3)]
                    R = [r for r, _ in draws]
                    out["cancer"][f"{truth}|{cand}|{gate}"] = dict(
                        pooled=float(np.mean([r["pooled"] for r in R])),
                        pertarget=float(np.mean([r["pertarget"] for r in R])) if all(r["pertarget"] is not None for r in R) else None,
                        nmae=float(np.mean([r["nmae"] for r in R])),
                        d_pooled=float(np.mean([r["d_pooled"] for r in R])),
                        d_pertarget=float(np.mean([r["d_pertarget"] for r in R])) if all(r["d_pertarget"] is not None for r in R) else None,
                        d_nmae=float(np.mean([r["d_nmae"] for r in R])),
                        ci_pooled=np.quantile(B[0], [.025, .975]).tolist(),
                        ci_pertarget=np.quantile(B[1], [.025, .975]).tolist(),
                        ci_nmae=np.quantile(B[2], [.025, .975]).tolist(),
                        p_pooled_gt0=float((B[0] > 0).mean()),
                        folds_pos=int(sum(r["d_pooled"] > 0 for r in R)),
                        per_fold_d=[r["d_pooled"] for r in R], k=[r["k"] for r in R])
    # ------------------------------------------------------------------ cross-study view: every candidate built only from
    # donors of a different study than L (K562 fold: its Orion sources already are; Orion folds: the K562-only rows)
    XS = {"whitened": {"k562_gwps": "PP", "orion_hct116": "PPxs", "orion_hek293t": "PPxs"},
          "gamma0": {"k562_gwps": "PPg0", "orion_hct116": "PPxsg0", "orion_hek293t": "PPxsg0"},
          "rawmean": {"k562_gwps": "M3L", "orion_hct116": "M3S", "orion_hek293t": "M3S"},
          **{f"blend_w{w:g}": {"k562_gwps": f"BLw{w:g}", "orion_hct116": f"BXw{w:g}", "orion_hek293t": f"BXw{w:g}"}
             for w in (1.0, 1.5, 2.0)}}
    out["xstudy"] = {}
    for truth in ("raw", "cen"):
        for gate in GATES:
            for arm, mp in XS.items():
                R, B = [], []
                for f in CANCER:
                    ref, cand = XS["whitened"][f], mp[f]
                    p, pb = point(Z[f], truth, gate, cand, masks[f]), point(Z[f], truth, gate, ref, masks[f])
                    R.append((p, pb))
                    B.append(fold_boot(Z[f], truth, gate, cand, masks[f], idx[f], ref=ref))
                Bm = [np.mean([b[i] for b in B], 0) for i in range(3)]
                out["xstudy"][f"{truth}|{arm}|{gate}"] = dict(
                    per_fold_d=[p["pooled"] - pb["pooled"] for p, pb in R],
                    d_pooled=float(np.mean([p["pooled"] - pb["pooled"] for p, pb in R])),
                    ci_pooled=np.quantile(Bm[0], [.025, .975]).tolist(),
                    d_nmae=float(np.mean([p["nmae"] - pb["nmae"] for p, pb in R])),
                    ci_nmae=np.quantile(Bm[2], [.025, .975]).tolist())
    (WORK / "analyse430.json").write_text(json.dumps(out, indent=1))

    # ------------------------------------------------------------------ markdown tables
    L = []
    for truth in ("raw", "cen"):
        for gate in GATES:
            key0 = f"{truth}|PP|{gate}"
            if key0 not in out["cancer"]:
                continue
            L.append(f"\n### truth {truth}, gate {gate} (cancer-fold mean; Δ vs PP; 95% bootstrap CI over targets)\n")
            L.append("| candidate | acc pooled | Δ pooled [CI] | folds + | Δ per-fold k562 / hct116 / hek293t | acc per-target | "
                     "Δ per-target [CI] | nmae proxy | Δ nmae [CI] | k (per fold) |")
            L.append("|---|---:|---|---:|---|---:|---|---:|---|---|")
            for cand in CANDS:
                c = out["cancer"].get(f"{truth}|{cand}|{gate}")
                if c is None:
                    continue
                f3 = " / ".join(f"{x:+.4f}" for x in c["per_fold_d"])
                pt = f"{c['pertarget']:.4f}" if c["pertarget"] is not None else "—"
                dpt = (f"{c['d_pertarget']:+.4f} [{c['ci_pertarget'][0]:+.4f}, {c['ci_pertarget'][1]:+.4f}]"
                       if c["d_pertarget"] is not None else "—")
                L.append(f"| {cand} | {c['pooled']:.4f} | {c['d_pooled']:+.4f} [{c['ci_pooled'][0]:+.4f}, "
                         f"{c['ci_pooled'][1]:+.4f}] | {c['folds_pos']}/3 | {f3} | {pt} | {dpt} | {c['nmae']:.4f} | "
                         f"{c['d_nmae']:+.4f} [{c['ci_nmae'][0]:+.4f}, {c['ci_nmae'][1]:+.4f}] | "
                         f"{' / '.join(f'{x:.3g}' for x in c['k'])} |")
    # kolf2 (centred truth only; its raw table is a z-score matrix without a control-relative baseline)
    L.append("\n### kolf2 (non-cancer reference; centred truth; Δ vs PP [CI])\n")
    L.append("| candidate | " + " | ".join(f"{g} acc / Δ" for g in ("top50", "top200", "top1000")) + " |")
    L.append("|---|" + "---|" * 3)
    for cand in CANDS:
        cells = []
        for g in ("top50", "top200", "top1000"):
            r = out["folds"].get(f"kolf2|cen|{cand}|{g}")
            cells.append("—" if r is None else f"{r['pooled']:.4f} / {r['d_pooled']:+.4f} [{r['ci_pooled'][0]:+.4f}, {r['ci_pooled'][1]:+.4f}]")
        L.append(f"| {cand} | " + " | ".join(cells) + " |")
    L.append("\n### cross-study only (K562 fold as is; Orion folds with the K562-only base); Δ vs the whitened composition\n")
    L.append("| truth | gate | arm | Δ pooled [CI] | per fold k562 / hct116 / hek293t | Δ nmae proxy [CI] |")
    L.append("|---|---|---|---|---|---|")
    for truth in ("raw", "cen"):
        for gate in ("z1000", "z200", "top200", "tau0.1"):
            for arm in XS:
                if arm == "whitened":
                    continue
                x = out["xstudy"][f"{truth}|{arm}|{gate}"]
                L.append(f"| {truth} | {gate} | {arm} | {x['d_pooled']:+.4f} [{x['ci_pooled'][0]:+.4f}, {x['ci_pooled'][1]:+.4f}] | "
                         f"{' / '.join(f'{v:+.4f}' for v in x['per_fold_d'])} | {x['d_nmae']:+.4f} [{x['ci_nmae'][0]:+.4f}, {x['ci_nmae'][1]:+.4f}] |")
    (WORK / "table430.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    print("support:", out["support"])


if __name__ == "__main__":
    main()

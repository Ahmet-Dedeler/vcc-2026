"""286 -- Stage B analysis on saved arrays only: paired target bootstraps vs `cd4svd5` (per fold and pooled), by 262
group and by the held reference's split-half reliability tercile, the marginal-band crossing at 0.02, leave-top-3-out
retention, and each arm against its frozen controls.

    PYTHONDONTWRITEBYTECODE=1 uv run python analyseB286.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
WORK = ROOT / "experiments/286-crossproduct-estimator/work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t")
B = 20_000
SEED = 20260906
rng = np.random.default_rng(SEED)
BASE = "cd4svd5"
LINES: list[str] = []
MD: list[str] = []
OUT: dict = dict(seed=SEED, B=B, folds=list(FOLDS))


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def boot_pool(diffs: list[np.ndarray]) -> dict:
    acc = np.zeros(B)
    for diff in diffs:
        d = diff[np.isfinite(diff)]
        idx = rng.integers(0, d.size, size=(B, d.size))
        acc += d[idx].mean(1)
    means = acc / len(diffs)
    return dict(point=float(np.mean([np.nanmean(d) for d in diffs])), lo=float(np.percentile(means, 2.5)),
                hi=float(np.percentile(means, 97.5)), sd=float(means.std(ddof=1)),
                p_le0=float((means <= 0).mean()), per_fold=[float(np.nanmean(d)) for d in diffs],
                n=[int(np.isfinite(d).sum()) for d in diffs])


def fb(b: dict) -> str:
    return f"{b['point']:+.5f} [{b['lo']:+.5f}, {b['hi']:+.5f}]"


D = {}
for f in FOLDS:
    sw = json.loads((WORK / f"sweepB286-{f}.json").read_text())
    per = {k: np.asarray(v) for k, v in np.load(WORK / f"pertargetB286-{f}.npz", allow_pickle=True).items()}
    idn = {k: np.asarray(v, dtype=np.float64) for k, v in np.load(WORK / f"identityB286-{f}.npz").items()}
    base_cos = np.asarray(per["base262_cos_truth"], dtype=np.float64)
    group = np.where(base_cos <= 0, "e_anti", np.where(base_cos <= 0.05, "marginal", "saturated"))
    ref = np.asarray(per["ref_rel"], dtype=np.float64)
    cuts = np.nanpercentile(ref, [100 / 3, 200 / 3])
    terc = np.where(ref <= cuts[0], "T1 low", np.where(ref <= cuts[1], "T2 mid", "T3 high"))
    D[f] = dict(meta=sw["meta"], rows={r["arm"]: r for r in sw["rows"]}, per=per, idn=idn, group=group, terc=terc)
    m = sw["meta"]
    P(f"[{f}] {m['n_targets']} targets, donors {m['donors']}, KA {m['ka_pass']}/{m['ka_total']}, failures "
      f"{m['failures']}, {m['n_rows']} rows, {m['seconds']}s, peak {m['peak_rss_gib']:.2f} GiB; est286 "
      f"{m['est286_sha256'][:12]}")
    OUT.setdefault("ka", {})[f] = m["ka"]
    OUT.setdefault("est286_sha256", m["est286_sha256"])
    assert m["est286_sha256"] == OUT["est286_sha256"], "est286.py changed between folds"

ARMS = [a for a in D[FOLDS[0]]["rows"] if all(a in D[f]["rows"] for f in FOLDS)]
MAIN = [a for a in ARMS if a.startswith(("e1:", "e2:", "E3:", "E4:", "E5:")) and ":mirror" not in a
        and ":shufpair" not in a and ":wrong" not in a and ":rand" not in a]
CONTROLS = [a for a in ARMS if ":mirror" in a or ":shufpair" in a or ":wrong" in a or ":rand" in a]
P(f"\narms {len(ARMS)}: main {len(MAIN)}, controls {len(CONTROLS)}")


def nmae_concat(arm):
    U = np.concatenate([D[f]["idn"][f"{arm}|U"] for f in FOLDS])
    M = np.concatenate([D[f]["idn"][f"{arm}|M"] for f in FOLDS])
    Dd = np.concatenate([D[f]["idn"][f"{arm}|D"] for f in FOLDS])
    return float(np.mean(1 + (U - 2 * M) / Dd))


nm_base = nmae_concat(BASE)
MD.append("## M1. Every arm vs `cd4svd5`, three cancer folds (Δ`pds` per fold, pooled paired target bootstrap)\n")
MD.append("| arm | Δ`pds` K562 / HCT116 / HEK293T | **pooled [95 % CI]** | P(≤ 0) | Δ`cos_agg` % | Δ macro % | delivered | Δ`nmae` id | Δ prec89 | rows changed |")
MD.append("|---|---|---|---:|---:|---:|---:|---:|---:|---|")
OUT["arms"] = {}
for arm in ARMS:
    diffs = [D[f]["per"][f"{arm}|pds"] - D[f]["per"][f"{BASE}|pds"] for f in FOLDS]
    bp = boot_pool(diffs)
    ra = [D[f]["rows"][arm]["cos_agg"] / D[f]["rows"][BASE]["cos_agg"] - 1 for f in FOLDS]
    rm = [D[f]["rows"][arm]["a23_cos"] / D[f]["rows"][BASE]["a23_cos"] - 1 for f in FOLDS]
    dn = np.mean([D[f]["rows"][arm]["delivered_norm"] for f in FOLDS]) / np.mean([D[f]["rows"][BASE]["delivered_norm"] for f in FOLDS])
    dnm = nmae_concat(arm) - nm_base
    dprec = float(np.mean([D[f]["rows"][arm]["top89_precision"] - D[f]["rows"][BASE]["top89_precision"] for f in FOLDS]))
    chg = [D[f]["rows"][arm]["n_rows_changed"] for f in FOLDS]
    OUT["arms"][arm] = dict(pooled=bp, rel_agg=[float(x) for x in ra], rel_macro=[float(x) for x in rm],
                            delivered_ratio=float(dn), d_nmae=dnm, d_prec89=dprec, rows_changed=chg)
    MD.append(f"| `{arm}` | {' / '.join(f'{x:+.4f}' for x in bp['per_fold'])} | **{fb(bp)}** | {bp['p_le0']:.4f} | "
              f"{100*np.mean(ra):+.3f} | {100*np.mean(rm):+.3f} | {dn:.4f} | {dnm:+.5f} | {dprec:+.4f} | {'/'.join(map(str, chg))} |")
    P(f"  {arm:26s} d pds {' / '.join(f'{x:+.4f}' for x in bp['per_fold'])}  pooled {fb(bp)} P<=0 {bp['p_le0']:.3f}  "
      f"agg {100*np.mean(ra):+.3f}%  dn {dn:.4f}  d nmae {dnm:+.5f}")

# arm vs its controls
MD.append("\n## M2. Each estimator against its frozen controls (paired target bootstrap of the per-target difference)\n")
MD.append("| arm | vs `:mirror` (Cm) | vs `:shufpair` (Cs) | vs `:wrong` (Cw) | vs `:rand` (Cr) | vs `E4` | vs `E5` |")
MD.append("|---|---|---|---|---|---|---|")
OUT["vs_controls"] = {}
E4 = next((a for a in ARMS if a.startswith("E4")), None)
E5 = next((a for a in ARMS if a.startswith("E5")), None)
for arm in MAIN:
    rec, cells = {}, []
    for suf in (":mirror", ":shufpair", ":wrong", ":rand"):
        other = arm + suf if (arm + suf) in ARMS else None
        if other is None:
            rec[suf] = None
            cells.append("—")
            continue
        rec[suf] = boot_pool([D[f]["per"][f"{arm}|pds"] - D[f]["per"][f"{other}|pds"] for f in FOLDS])
        cells.append(fb(rec[suf]))
    for other in (E4, E5):
        if other is None or other == arm:
            rec[other] = None
            cells.append("—")
            continue
        rec[other] = boot_pool([D[f]["per"][f"{arm}|pds"] - D[f]["per"][f"{other}|pds"] for f in FOLDS])
        cells.append(fb(rec[other]))
    OUT["vs_controls"][arm] = rec
    MD.append(f"| `{arm}` | " + " | ".join(cells) + " |")

# strata
strata = [("group (e) base·truth ≤ 0", "group", "e_anti"), ("marginal band (0, 0.05]", "group", "marginal"),
          ("saturated > 0.05", "group", "saturated"), ("ref-rel T1 (low)", "terc", "T1 low"),
          ("ref-rel T2 (mid)", "terc", "T2 mid"), ("ref-rel T3 (high)", "terc", "T3 high")]
MD.append("\n## M3. Δ`pds` vs `cd4svd5` by 262 group and by the held reference's reliability tercile (pooled)\n")
MD.append("| arm | " + " | ".join(f"{s[0]} (n {sum(int((D[f][s[1]] == s[2]).sum()) for f in FOLDS)})" for s in strata) + " |")
MD.append("|---|" + "---|" * len(strata))
OUT["strata"] = {}
for arm in MAIN + CONTROLS:
    rec, cells = {}, []
    for label, key, val in strata:
        b = boot_pool([np.where(D[f][key] == val, D[f]["per"][f"{arm}|pds"] - D[f]["per"][f"{BASE}|pds"], np.nan) for f in FOLDS])
        rec[label] = b
        cells.append(fb(b))
    OUT["strata"][arm] = rec
    MD.append(f"| `{arm}` | " + " | ".join(cells) + " |")

# the marginal band / knee crossings
MD.append("\n## M4. The marginal band (`cd4svd5`'s own cosine in (0, 0.05]): knee crossings at 0.02\n")
MD.append("| arm | band n | mean Δ own cosine on the band | up | down | net share | Δ`pds` on the band | Δ own cosine, all |")
MD.append("|---|---:|---:|---:|---:|---:|---|---:|")
OUT["band"] = {}
for arm in MAIN + CONTROLS:
    nb = up = down = 0
    dc, dpds, dc_all = [], [], []
    for f in FOLDS:
        c0, c1 = D[f]["per"][f"{BASE}|cos"], D[f]["per"][f"{arm}|cos"]
        band = (c0 > 0) & (c0 <= 0.05)
        nb += int(band.sum())
        up += int((band & (c0 < 0.02) & (c1 >= 0.02)).sum())
        down += int((band & (c0 >= 0.02) & (c1 < 0.02)).sum())
        dc.append(np.where(band, c1 - c0, np.nan))
        dc_all.append(c1 - c0)
        dpds.append(np.where(band, D[f]["per"][f"{arm}|pds"] - D[f]["per"][f"{BASE}|pds"], np.nan))
    b = boot_pool(dpds)
    rec = dict(n=nb, up=up, down=down, net_share=(up - down) / nb, d_cos_band=float(np.nanmean(np.concatenate(dc))),
               d_pds_band=b, d_cos_all=float(np.nanmean(np.concatenate(dc_all))))
    OUT["band"][arm] = rec
    MD.append(f"| `{arm}` | {nb} | {rec['d_cos_band']:+.5f} | {up} | {down} | {100*rec['net_share']:+.1f} % | {fb(b)} | {rec['d_cos_all']:+.5f} |")

# leave-top-3-out retention
MD.append("\n## M5. Leave-top-3-out retention of the per-target Δ`pds` (261's statistic)\n")
MD.append("| arm | full | leave-top-3-out | retained | per fold |")
MD.append("|---|---:|---:|---:|---|")
OUT["retention"] = {}
for arm in MAIN + CONTROLS:
    full, loo = [], []
    for f in FOLDS:
        d = D[f]["per"][f"{arm}|pds"] - D[f]["per"][f"{BASE}|pds"]
        order = np.argsort(-d)
        keep = np.ones(d.size, bool)
        keep[order[:3]] = False
        full.append(d.mean())
        loo.append(d[keep].mean())
    fm, lm = float(np.mean(full)), float(np.mean(loo))
    rec = dict(full=fm, leave_top_out=lm, retained=(lm / fm if fm else float("nan")),
               per_fold=[float(l / g) if g else float("nan") for l, g in zip(loo, full)])
    OUT["retention"][arm] = rec
    MD.append(f"| `{arm}` | {fm:+.5f} | {lm:+.5f} | {100*rec['retained']:+.0f} % | "
              f"{' / '.join(f'{100*x:+.0f} %' for x in rec['per_fold'])} |")

(WORK / "analysisB286.json").write_text(json.dumps(OUT, indent=1, default=float))
(WORK / "tablesB286.md").write_text("\n".join(MD) + "\n")
(WORK / "analyseB286.log").write_text("\n".join(LINES) + "\n")
P("\nwritten work/analysisB286.json, work/tablesB286.md, work/analyseB286.log")

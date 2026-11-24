"""522 — tables, leave-one-panel-out and October pricing from arms522's per-surface jsons (no truth recomputation).

Per panel (A/B/C: the three cancer folds, kolf2 reported beside; efp / efr / real: their three cancer folds):
  r      relative aggregate cosine vs the cap-3 key (the October default), mean over the cancer folds; n+ folds positive
  dpds   literal PLOLO pds change (board: x 0.34, 540's key-step calibration)
Leave-one-panel-out (LOPO): each FAMILY's free parameter is chosen to maximise the mean of the training panels' mean r,
then scored on the held-out panel. Two versions: over the three October-like panels, and over all four surfaces.
Pricing (board units, overall), per panel with its retention ret (438: cos_agg(PP|K3) on the panel's fold / A/B/C's):
  overall ~= 0.0710 x k x r x ret + 0.1245 x dpds_lit   (518's identity: +1 % X = -0.0041 raw mse at k 1 on A/B/C;
                                                          pds x 0.34 / 0.455 / 6), k = 1.31 (key class) ... 1.

    uv run python experiments/522-october-bulk-key/analyse522.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
PANELS = ("abc", "efp", "efr", "real")
OCT = ("efp", "efr", "real")
# 438's retention (cos_agg(PP|K3) on the panel's fold / on A/B/C's fold), per cancer fold
RET = {"abc": {f: 1.0 for f in CANCER},
       "efp": {"k562_gwps": 0.775, "orion_hct116": 0.614, "orion_hek293t": 0.528},
       "efr": {"k562_gwps": 0.347, "orion_hct116": 0.327, "orion_hek293t": 0.315},
       "real": {"k562_gwps": 1.011, "orion_hct116": 0.852, "orion_hek293t": 0.903}}
PPROD = "pctPROD"                         # 497's unrounded P 96.64697522769474 (the candidate October form)


def load():
    R = {}
    for p in PANELS:
        for f in (CANCER + (("kolf2",) if p == "abc" else ())):
            path = WORK / f"arms522-{p}-{f}.json"
            if path.exists():
                R[(p, f)] = json.loads(path.read_text())
    return R


def families(names):
    fam = {
        "fixed cap": [("cap3", 3.0)] + [(n, float(n[3:])) for n in names if n.startswith("cap") and "x" not in n and n != "cap3"],
        "percentile cap": [(n, 96.64697522769474 if n == PPROD else float(n[3:])) for n in names
                           if n.startswith("pct") and "x" not in n],
        "strength on cap 3 (src)": [("cap3", 0.0)] + [(n, float(n[1:])) for n in names if n.startswith("S") and n[1:2].isdigit()],
        "strength+CD4 on cap 3": [("cap3", 0.0)] + [(n, float(n.split("_")[1])) for n in names if n.startswith("Scd4_")],
        "cap 6 x strength": [("cap6", 0.0)] + [(n, float(n.split("xS")[1])) for n in names if n.startswith("cap6xS")],
        "pct x strength": [(PPROD, 0.0)] + [(n, float(n.split("xS")[1])) for n in names if n.startswith("pctPRODxS")],
    }
    return fam


def panel_stat(R, p, nm, key="r"):
    vals = [R[(p, f)]["arms"][nm][key] for f in CANCER if (p, f) in R]
    return float(np.mean(vals)) if vals else np.nan, vals


def price(R, p, nm, k):
    out = []
    for f in CANCER:
        if (p, f) not in R:
            continue
        a = R[(p, f)]["arms"][nm]
        out.append(0.0710 * k * a["r"] * RET[p][f] + 0.1245 * a["dpds"])
    return float(np.mean(out)) if out else np.nan


def main():
    R = load()
    have = sorted({p for p, _ in R})
    names = list(next(iter(R.values()))["arms"])
    lines = []
    P = lambda *a: lines.append(" ".join(str(x) for x in a))
    # ---------------- per-panel tables
    for p in have:
        folds = [f for f in CANCER if (p, f) in R]
        P(f"\n## {p}: relative cos vs cap 3 (cancer mean, folds +), per fold {' / '.join(f[:6] for f in folds)}"
          + ("  | kolf2" if (p, "kolf2") in R else "") + "  | literal dpds | w_max | overall k1.31 (k1)")
        for nm in names:
            m, v = panel_stat(R, p, nm)
            dp, _ = panel_stat(R, p, nm, "dpds")
            wm = [R[(p, f)]["arms"][nm]["w_max"] for f in folds]
            kol = f"  | {100 * R[(p, 'kolf2')]['arms'][nm]['r']:+6.2f}" if (p, "kolf2") in R else ""
            P(f"{nm:22s} {100 * m:+7.2f} % ({sum(x > 0 for x in v)}/{len(v)})  " + " / ".join(f"{100 * x:+6.2f}" for x in v)
              + kol + f"  | {dp:+.5f} | {min(wm):5.2f}-{max(wm):5.2f} | {price(R, p, nm, 1.31):+.4f} ({price(R, p, nm, 1.0):+.4f})")
    # ---------------- hi of the percentile cap per surface
    P("\n## the percentile cap's hi (P 96.647) per surface")
    for (p, f), res in sorted(R.items()):
        P(f"{p:5s} {f:14s} hi {res['his']['hi_pct_546']:.3f}  n {res['n']} carried {res['n_carried']}")
    # ---------------- LOPO
    fam = families(names)
    out = dict(lopo={}, tables={})
    for label, pool in (("October-like (efp, efr, real)", [p for p in OCT if p in have]),
                        ("all four surfaces", [p for p in PANELS if p in have])):
        if len(pool) < 2:
            continue
        P(f"\n## LOPO over {label}: parameter chosen on the other panels (mean of panel means), scored held out")
        P(f"{'family':26s} " + " ".join(f"{h:>18s}" for h in pool) + "   mean held-out   worst")
        out["lopo"][label] = {}
        for crit, agg in (("mean", np.mean), ("minimax", np.min)):
          P(f"-- criterion: {crit} of the training panels' mean r")
          for fname, members in fam.items():
            row, held_vals = [], []
            for h in pool:
                train = [p for p in pool if p != h]
                best = max(members, key=lambda m: agg([panel_stat(R, p, m[0])[0] for p in train]))
                val = panel_stat(R, h, best[0])[0]
                held_vals.append(val)
                row.append(f"{100 * val:+6.2f}% @{best[1]:g}".rjust(18))
                out["lopo"][label].setdefault(crit, {}).setdefault(fname, {})[h] = dict(choice=best[0], param=best[1], heldout_r=val)
            P(f"{fname:26s} " + " ".join(row) + f"   {100 * np.mean(held_vals):+7.2f} %    {100 * min(held_vals):+7.2f} %")
            out["lopo"][label][crit][fname]["mean"] = float(np.mean(held_vals))
            out["lopo"][label][crit][fname]["worst"] = float(min(held_vals))
        # fixed choices for reference
        for nm in ("cap6", PPROD, "S1", "cap6xS0.5", "pctPRODxS0.5", "pctPRODxS0.25", "cap8", "cap1e+09"):
            if nm in names:
                v = [panel_stat(R, h, nm)[0] for h in pool]
                P(f"{'(fixed) ' + nm:26s} " + " ".join(f"{100 * x:+6.2f}%".rjust(18) for x in v)
                  + f"   {100 * np.mean(v):+7.2f} %    {100 * min(v):+7.2f} %")
    for p in have:
        out["tables"][p] = {nm: dict(r=panel_stat(R, p, nm)[0], folds=panel_stat(R, p, nm)[1],
                                     n_pos=int(sum(x > 0 for x in panel_stat(R, p, nm)[1])),
                                     dpds=panel_stat(R, p, nm, "dpds")[0],
                                     overall_k131=price(R, p, nm, 1.31), overall_k1=price(R, p, nm, 1.0),
                                     kolf2=R[(p, "kolf2")]["arms"][nm]["r"] if (p, "kolf2") in R else None)
                            for nm in names}
    txt = "\n".join(lines)
    print(txt)
    (WORK / "analyse522.txt").write_text(txt + "\n")
    (WORK / "analyse522.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()

"""228 — assemble the PLOLO arms: dose curves, the three cancer folds, `kolf2` apart, and
the leave-one-FOLD-out dose selection (218/207/223's rule: the dose is chosen on the other
folds and the price is read on the held one, so selection and reporting never overlap).

    uv run python analyse228.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
ALL = CANCER + ("kolf2",)
SPAN_PDS, CALIB = 0.452, 1.011
SIGMA_PDS = 0.001494


def load():
    rows = []
    for f in ALL:
        p = WORK / f"arms228-{f}.json"
        if p.exists():
            rows += json.loads(p.read_text())
    return rows


def main() -> None:
    rows = load()
    by = {(r["fold"], r["arm"], r["a"]): r for r in rows}
    base = {f: by[(f, "BASE", 1.0)] for f in ALL if (f, "BASE", 1.0) in by}
    folds = sorted(base)
    variants, doses = {}, set()
    for r in rows:
        if r["arm"] == "BASE" or r["a"] != 1.0 or r["variant"] == "none":
            continue
        variants.setdefault(r["variant"], set()).add(r["t"])
        doses.add(r["t"])

    def rel(f, arm, key="cos_agg"):
        b = base[f][key]
        return 100.0 * (by[(f, arm, 1.0)][key] - b) / abs(b)

    def d(f, arm, key):
        return by[(f, arm, 1.0)][key] - base[f][key]

    out = {"per_variant": {}}
    print(f"{'variant':<16}{'t':>6}" + "".join(f"{f[:11]:>13}" for f in folds)
          + f"{'CANCER rel%':>13}{'CANCER dpds':>13}{'kolf2 rel%':>12}")
    print("-" * (22 + 13 * len(folds) + 38))
    for v in sorted(variants):
        for t in sorted(variants[v]):
            arm = f"{v}@{t:g}"
            have = [f for f in folds if (f, arm, 1.0) in by]
            if not have:
                continue
            cell = "".join(f"{rel(f, arm):+13.2f}" if f in have else f"{'-':>13}"
                           for f in folds)
            cc = [rel(f, arm) for f in CANCER if f in have]
            cp = [d(f, arm, "pds") for f in CANCER if f in have]
            kk = rel("kolf2", arm) if "kolf2" in have else float("nan")
            print(f"{v:<16}{t:>6.2f}{cell}{np.mean(cc):+13.3f}"
                  f"{np.mean(cp):+13.5f}{kk:+12.2f}")
            out["per_variant"][arm] = dict(
                cancer_rel_cos_agg=float(np.mean(cc)),
                cancer_d_pds=float(np.mean(cp)),
                cancer_d_qw=float(np.mean([d(f, arm, "a15_sgn_w") for f in CANCER
                                           if f in have])),
                cancer_rel_macro_cos=float(np.mean([rel(f, arm, "cos") for f in CANCER
                                                    if f in have])),
                per_fold_rel_cos_agg={f: rel(f, arm) for f in have},
                per_fold_d_pds={f: d(f, arm, "pds") for f in have},
                kolf2_rel_cos_agg=None if "kolf2" not in have else kk,
                kolf2_d_pds=None if "kolf2" not in have else d("kolf2", arm, "pds"),
                pds_positive_cancer=int(sum(1 for x in cp if x > 0)),
                n_cancer=len(cc))

    # ---------------- leave-one-FOLD-out dose selection, on the three cancer folds
    print(f"\n{'LEAVE-ONE-FOLD-OUT dose selection (3 cancer folds)':^80}")
    print(f"{'variant':<16}{'held':<16}{'t* (other 2)':>14}{'rel% held':>12}"
          f"{'dpds held':>12}")
    lofo = {}
    for v in sorted(variants):
        ts = sorted(variants[v])
        picks = []
        for held in CANCER:
            others = [f for f in CANCER if f != held]
            avail = [t for t in ts if all((f, f"{v}@{t:g}", 1.0) in by
                                          for f in others + [held])]
            if not avail:
                continue
            tstar = max(avail, key=lambda t: np.mean(
                [rel(f, f"{v}@{t:g}") for f in others]))
            arm = f"{v}@{tstar:g}"
            picks.append(dict(held=held, t=tstar, rel=rel(held, arm),
                              d_pds=d(held, arm, "pds"),
                              d_qw=d(held, arm, "a15_sgn_w"),
                              rel_macro=rel(held, arm, "cos")))
            print(f"{v:<16}{held:<16}{tstar:>14.2f}{picks[-1]['rel']:>+12.3f}"
                  f"{picks[-1]['d_pds']:>+12.5f}")
        if picks:
            lofo[v] = dict(
                rel_cos_agg=float(np.mean([p["rel"] for p in picks])),
                rel_macro_cos=float(np.mean([p["rel_macro"] for p in picks])),
                d_pds=float(np.mean([p["d_pds"] for p in picks])),
                d_qw=float(np.mean([p["d_qw"] for p in picks])),
                pds_positive=int(sum(1 for p in picks if p["d_pds"] > 0)),
                picks=picks)
            print(f"{'':<16}{'LOFO MEAN':<16}{'':>14}"
                  f"{lofo[v]['rel_cos_agg']:>+12.3f}{lofo[v]['d_pds']:>+12.5f}"
                  f"   pds+ {lofo[v]['pds_positive']}/3")
    out["lofo"] = lofo

    # ---------------- the VCC_ALPHA interaction
    print(f"\n{'VCC_ALPHA interaction (a = multiplier on ALPHA = 0.982)':^80}")
    alphas = sorted({r["a"] for r in rows})
    print(f"{'arm':<16}" + "".join(f"{a:>13.3f}" for a in alphas))
    for arm in ("none@0", "gain@1", "gain_reg@1"):
        line, vals = "", {}
        for a in alphas:
            v = [by[(f, arm, a)]["cos_agg"] for f in CANCER if (f, arm, a) in by]
            r = [by[(f, "none@0", a)]["cos_agg"] if (f, "none@0", a) in by else None
                 for f in CANCER]
            if not v or any(x is None for x in r):
                line += f"{'-':>13}"
                continue
            g = 100.0 * (np.mean(v) - np.mean(r)) / abs(np.mean(r))
            vals[a] = dict(cos_agg=float(np.mean(v)), rel_vs_same_alpha=float(g))
            line += f"{np.mean(v):>13.6f}"
        print(f"{arm:<16}{line}")
        out.setdefault("alpha", {})[arm] = vals
    print(f"{'':<16}" + "  relative to `none@0` at the SAME alpha:")
    for arm in ("gain@1", "gain_reg@1"):
        line = "".join(f"{out['alpha'][arm][a]['rel_vs_same_alpha']:>+13.3f}"
                       if a in out.get("alpha", {}).get(arm, {}) else f"{'-':>13}"
                       for a in alphas)
        print(f"{arm:<16}{line}")

    (WORK / "analyse228.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> {WORK / 'analyse228.json'}")


if __name__ == "__main__":
    main()

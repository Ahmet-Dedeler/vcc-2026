"""223 — reduce `focus223.json`: dose curves, all falsifiers, LOFO pricing for the survivors.

    uv run python analyse_focus.py
"""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
RATE_FIVE = 0.003634 + 0.002446 + 0.000980 + 0.000847 + 0.000100   # 218 section 7.4
PDS_SPAN = 0.452


def price(d_cos_abs: float, d_pds_raw: float, k: float = 1.0) -> float:
    return RATE_FIVE * (d_cos_abs * k / 0.01) + (d_pds_raw / PDS_SPAN) / 6.0 * 1.011 * 2.0


def main() -> None:
    rows = json.loads((WORK / "focus223.json").read_text())
    base = {(r["fold"], r["base"]): r for r in rows if r["kind"] == "baseline"}
    arms = [r for r in rows if r["kind"] == "arm"]
    idx = defaultdict(dict)
    for r in arms:
        idx[(r["base"], r["source"], r["g"])][r["fold"]] = r
    folds = sorted({r["fold"] for r in arms})
    doses = sorted({r["g"] for r in arms})
    sources = sorted({r["source"] for r in arms})
    real = [s for s in sources if "|" not in s]
    out: dict[str, object] = {}

    def m(bn, s, g, fs, key):
        v = [idx[(bn, s, g)][f][key] for f in fs if f in idx[(bn, s, g)]]
        return float(np.mean(v)) if v else float("nan")

    for bn in ("base222", "cd4"):
        title = ("MARGINAL over cd4@0.34 + h1pair@0.28 (222's arm)" if bn == "base222"
                 else "ALONE on top of cd4@0.34 (222 not shipped)")
        print(f"\n### DOSE CURVE — {title}, three cancer folds")
        h = (f"{'source':16s}{'g':>6s}{'rel dcos%':>11s}{'dcos':>11s}{'dpds':>10s}"
             f"{'dq_w':>10s}{'SHUF%':>9s}{'SHUF2%':>9s}{'CD4rows%':>10s}"
             f"{'carry%':>9s}{'noncar pp':>10s}{'price':>10s}")
        print(h)
        print("-" * len(h))
        for s in real:
            for g in doses:
                if not idx[(bn, s, g)]:
                    continue
                dc = m(bn, s, g, CANCER, "d_cos")
                dp = m(bn, s, g, CANCER, "d_pds")
                rc = m(bn, s, g, CANCER, "rel_cos_pct")
                sh = m(bn, f"{s}|SHUF", g, CANCER, "rel_cos_pct")
                s2 = m(bn, f"{s}|SHUF2", g, CANCER, "rel_cos_pct")
                cr = m(bn, f"{s}|CD4ROWS", g, CANCER, "rel_cos_pct")
                car = float(np.mean([idx[(bn, s, g)][f]["rowsplit"]["carrying_rel_pct"]
                                     for f in CANCER if f in idx[(bn, s, g)]]))
                nc = float(np.max([abs(idx[(bn, s, g)][f]["rowsplit"]
                                       .get("non_carrying_rel_pct", 0.0))
                                   for f in CANCER if f in idx[(bn, s, g)]]))
                p = price(dc, dp)
                out.setdefault(f"{bn}_curve", []).append(
                    dict(base=bn, source=s, g=g, rel_cos=rc, d_cos=dc, d_pds=dp,
                         d_qw=m(bn, s, g, CANCER, "d_qw"), shuf=sh, shuf2=s2,
                         cd4rows=cr, carry_rel=car, noncarry_pp=nc, price_k1=p))
                print(f"{s:16s}{g:>6.2f}{rc:>11.3f}{dc:>11.2e}{dp:>10.5f}"
                      f"{m(bn, s, g, CANCER, 'd_qw'):>10.5f}{sh:>9.3f}{s2:>9.3f}"
                      f"{cr:>10.3f}{car:>9.3f}{nc:>10.4f}{p:>10.5f}")

    # ------------------------------------------------------------- per-fold detail
    print("\n### PER-FOLD, marginal over base222, at each source's cancer-fold argmax")
    best = {}
    for s in real:
        cand = [(price(m("base222", s, g, CANCER, "d_cos"),
                       m("base222", s, g, CANCER, "d_pds")), g) for g in doses
                if idx[("base222", s, g)]]
        best[s] = max(cand)[1]
    hh = f"{'source':16s}{'g*':>6s}" + "".join(f"{f[:13]:>15s}" for f in folds)
    print(hh)
    print("-" * len(hh))
    for s in real:
        g = best[s]
        line = f"{s:16s}{g:>6.2f}"
        for f in folds:
            r = idx[("base222", s, g)].get(f)
            line += (f"{r['rel_cos_pct']:>8.3f}/{r['d_pds']:>+6.4f}" if r
                     else "".rjust(15))
        print(line)
    out["argmax_dose"] = best

    # ------------------------------------------------------------------- LOFO
    print("\n### LEAVE-ONE-FOLD-OUT (dose chosen on the other three, priced on the held one)")
    print(f"{'source':16s}{'held':16s}{'g chosen':>10s}{'rel dcos%':>11s}"
          f"{'dpds':>10s}{'price':>10s}")
    lofo = []
    for s in real:
        tot = []
        for held in folds:
            others = [f for f in folds if f != held]
            cand = [(price(m("base222", s, g, others, "d_cos"),
                           m("base222", s, g, others, "d_pds")), g)
                    for g in doses if idx[("base222", s, g)]]
            if not cand:
                continue
            g = max(cand)[1]
            r = idx[("base222", s, g)].get(held)
            if r is None:
                continue
            p = price(r["d_cos"], r["d_pds"])
            tot.append(p)
            lofo.append(dict(source=s, held=held, g=g, rel_cos=r["rel_cos_pct"],
                             d_pds=r["d_pds"], price=p))
            print(f"{s:16s}{held:16s}{g:>10.2f}{r['rel_cos_pct']:>11.3f}"
                  f"{r['d_pds']:>10.5f}{p:>10.5f}")
        if tot:
            print(f"{s:16s}{'MEAN(4 folds)':16s}{'':>10s}{'':>11s}{'':>10s}"
                  f"{np.mean(tot):>10.5f}   ({sum(t > 0 for t in tot)}/{len(tot)} positive)")
    out["lofo"] = lofo
    (WORK / "analyse_focus.json").write_text(json.dumps(out, indent=2))
    print(f"\n-> {WORK / 'analyse_focus.json'}")


if __name__ == "__main__":
    main()

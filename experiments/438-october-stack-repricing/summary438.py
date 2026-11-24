"""438 summary: the cross-panel tables of RESULTS.md (every number is price438's; this only arranges them).

    uv run python summary438.py [--panels abc,abcdef,real,efp,efr] [--full efr=efr_s2]  -> work/summary438.{json,md}
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import argparse  # noqa: E402
import json  # noqa: E402

import numpy as np  # noqa: E402

from price438 import CANCER, REF, SLOTS, WORK, Panel, october_union  # noqa: E402

STACK = "h2c1.5|K3+sh1.25"
H1S, CD4S = (0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5), (0.75, 1, 1.25, 1.5, 1.75, 2)
# calibrated sharpening transfer from 417a's board row (FINDINGS 0: k_sh 0.42 from pds, 0.65 from mse)
KSH, KSH_PDS = 0.65, 0.42


def cell(P, arm, ref=REF, level=True, k=1.0, k_pds=1.0, boot=False):
    p = {f: P.price(f, arm, k, ref, level, k_pds=k_pds)["ov"] for f in P.folds}
    m = float(np.mean(list(p.values())))
    s = f"{m:+.5f} {sum(v > 0 for v in p.values())}/3"
    b = None
    if boot:
        b = P.boot(arm, ref, level, k=k)
        s += f" [{b['ci'][0]:+.4f}, {b['ci'][1]:+.4f}]"
    return s, dict(per_fold=p, mean=m, boot=b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", default="abc,abcdef,real,efp,efr")
    ap.add_argument("--full", default="", help="panel=full_tag overrides (efr's surrogate)")
    ap.add_argument("--fold-ops", default="held")
    ap.add_argument("--no-ret", action="store_true")
    ap.add_argument("--out", default="summary438")
    ap.add_argument("--exclude-oct4", action="store_true", help="drop the four essential screens' targets (Arc-like)")
    ap.add_argument("--fold-energy", default="", help="panels priced on the fold-route energies (no 3-donor operator)")
    a = ap.parse_args()
    fulls = dict(kv.split("=") for kv in a.full.split(",") if kv)
    Ps = {}
    for t in a.panels.split(","):
        try:
            Ps[t] = Panel(t, a.fold_ops, full_tag=fulls.get(t), ret=not a.no_ret,
                          exclude=october_union() if a.exclude_oct4 else None,
                          energy="fold" if t in a.fold_energy.split(",") else "full")
            assert len(Ps[t].folds) == 3
        except (FileNotFoundError, AssertionError) as e:
            print(f"skip {t}: {e}")
            Ps.pop(t, None)
    tags = list(Ps)
    J, md = {}, []
    hdr = "| | " + " | ".join(tags) + " |\n|---|" + "---|" * len(tags)

    def table(title, rows, level, boot=False):
        md.append(f"\n**{title}**\n\n" + hdr)
        for lab, arm, ref, kw in rows:
            cells = []
            for t in tags:
                P = Ps[t]
                if not all(arm + (":L" if level else "") in P.R[f] and ref + (":L" if level else "") in P.R[f] for f in P.folds):
                    cells.append("-")
                    continue
                s, d = cell(P, arm, ref, level, boot=boot, **kw)
                J.setdefault(title, {}).setdefault(lab, {})[t] = d
                cells.append(s)
            md.append(f"| {lab} | " + " | ".join(cells) + " |")

    for level in (True, False):
        lv = "LEVELLED (D/E/F default)" if level else "unlevelled (420's primary route)"
        table(f"Stack and halves vs each panel's 417a analogue, {lv}, k = 1: cancer mean n+ [95% CI of the mean]", [
            ("H1 x2 half", "h2c1|K3+sh1.25", REF, {}), ("CD4 x1.5 half", "h1c1.5|K3+sh1.25", REF, {}),
            ("stack (H1 x2 + CD4 x1.5)", STACK, REF, {}), ("stack, k = 1.407", STACK, REF, dict(k=1.407181520135006)),
            ("CD4 x1.5 on top of H1 x2", STACK, "h2c1|K3+sh1.25", {}), ("H1 x2 on top of CD4 x1.5", STACK, "h1c1.5|K3+sh1.25", {}),
            ("CONTROL: V1 key permuted (PP)", "CTRL:h1c1|V1shuf+sh1.25", REF, {})], level, boot=True)
        table(f"Stack neighbours vs the A/B/C stack (H1 x2, CD4 x1.5, q 1.25), {lv}, k = 1", [
            *[(f"H1 x{m:g}, CD4 x1.5", f"h{m:g}c1.5|K3+sh1.25", STACK, {}) for m in H1S if m != 2],
            *[(f"H1 x2, CD4 x{m:g}", f"h2c{m:g}|K3+sh1.25", STACK, {}) for m in CD4S if m != 1.5]], level)
        table(f"Sharpening on the stack, {lv}: q step vs the stack at q 1.25", [
            *[(f"q {q} (k = 1)", f"h2c1.5|K3+sh{q}", STACK, {}) for q in ("1.35", "1.5", "1.75", "2")],
            *[(f"q {q} (417a-calibrated k 0.65 / k_pds 0.42)", f"h2c1.5|K3+sh{q}", STACK, dict(k=KSH, k_pds=KSH_PDS))
              for q in ("1.35", "1.5", "1.75", "2")],
            ("PP: q 1.5 vs q 1.25 (k = 1)", "h1c1|K3+sh1.5", REF, {})], level, boot=False)
        table(f"Keys on the stack, {lv}: vs the stack with K3 (410 key clip 0.3-3), k = 1", [
            ("LODO V1 clip (0.5, 3), q 1.25", "h2c1.5|V1_5_3+sh1.25", STACK, {}),
            ("LODO V1 clip (0.3, 3), q 1.25", "h2c1.5|V1+sh1.25", STACK, {}),
            ("410 key clip (0.5, 3), q 1.25", "h2c1.5|K5_3+sh1.25", STACK, {}),
            ("LODO V1 clip (0.5, 3), q 1.5 (vs K3 q 1.5)", "h2c1.5|V1_5_3+sh1.5", "h2c1.5|K3+sh1.5", {}),
            ("LODO V1 clip (0.5, 3) on PP, q 1.25", "h1c1|V1_5_3+sh1.25", REF, {})], level, boot=level)
    # ---------------------------------------------------------------- surfaces: argmax and the cost of A/B/C's weights
    md.append("\n**H1 x CD4 surface (K3), LEVELLED k = 1: argmax per panel, and what keeping A/B/C's weights (x2, x1.5) costs**\n\n"
              "| panel | q | argmax (H1, CD4) | argmax vs 417a analogue | A/B/C weights vs 417a analogue | cost of keeping A/B/C weights "
              "[paired 95% CI] (n+) | best on the H1 axis at CD4 x1.5 | best on the CD4 axis at H1 x2 |\n|---|---|---|---:|---:|---|---|---|")
    for t in tags:
        P = Ps[t]
        for q in ("1.25", "1.5"):
            inc = f"h2c1.5|K3+sh{q}"
            mean = {(mh, mc): float(np.mean([P.price(f, f"h{mh:g}c{mc:g}|K3+sh{q}", 1.0, REF, True)["ov"] for f in P.folds]))
                    for mh in H1S for mc in CD4S}
            best = max(mean, key=mean.get)
            bn = f"h{best[0]:g}c{best[1]:g}|K3+sh{q}"
            if best != (2, 1.5):
                b = P.boot(bn, inc, True)
                pf = {f: P.price(f, bn, 1.0, inc, True)["ov"] for f in P.folds}
                cost = f"{np.mean(list(pf.values())):+.5f} [{b['ci'][0]:+.4f}, {b['ci'][1]:+.4f}] ({sum(v > 0 for v in pf.values())}/3)"
            else:
                b, cost = None, "0 (A/B/C weights are the argmax)"
            bh = max(H1S, key=lambda m: mean[(m, 1.5)])
            bc = max(CD4S, key=lambda m: mean[(2, m)])
            J.setdefault("surface", {})[f"{t}|q{q}"] = dict(mean={f"{k[0]:g},{k[1]:g}": v for k, v in mean.items()}, best=best,
                                                              boot=b, best_h1_axis=bh, best_cd4_axis=bc)
            md.append(f"| {t} | {q} | x{best[0]:g}, x{best[1]:g} | {mean[best]:+.5f} | {mean[(2, 1.5)]:+.5f} | {cost} | "
                      f"x{bh:g} ({mean[(bh, 1.5)]:+.5f}) | x{bc:g} ({mean[(2, bc)]:+.5f}) |")
    # ---------------------------------------------------------------- energies
    md.append("\n**Energies: ρ = B/claim per context slot (A/B/C, or D/E/F on the synthetic bundles), unlevelled, anchored route**\n\n"
              "| panel | 417a analogue | H1 x2 half | CD4 x1.5 half | stack | stack q 1.5 | stack levelled (one ntclevel step) |\n|---|---|---|---|---|---|---|")
    for t in tags:
        E = Ps[t].E
        r = lambda n, k="rho": "/".join(f"{E[n][k][c]:.3f}" for c in SLOTS)
        md.append(f"| {t} | {r(REF)} | {r('h2c1|K3+sh1.25')} | {r('h1c1.5|K3+sh1.25')} | {r(STACK)} | {r('h2c1.5|K3+sh1.5')} | "
                  f"{r(STACK, 'rho_level')} |")
    (WORK / f"{a.out}.json").write_text(json.dumps(J, indent=1, default=lambda o: list(o) if isinstance(o, tuple) else str(o)))
    (WORK / f"{a.out}.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()

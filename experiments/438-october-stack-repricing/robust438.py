"""438 robustness of a step's price: 236's leave-top-3-out falsifier, and the essential / four-screen exclusions.

For each (panel, arm vs ref), levelled route, k = 1: the cancer-mean price (i) on all targets, (ii) without the union
of each fold's top-3 contributing targets (by the step's summed dx), (iii) without the four essential screens' targets
(the Arc-like panel, 436), (iv) without K562-essential targets.

    uv run python robust438.py [--fold-energy efr] [--full efr=efr_s1]      -> work/robust438.{json,md}
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import argparse  # noqa: E402
import json  # noqa: E402

import numpy as np  # noqa: E402

from price438 import CANCER, REF, ROOT, WORK, Panel, october_union  # noqa: E402

STEPS = [("H1 x2 half", "h2c1|K3+sh1.25", REF), ("CD4 x1.5 half", "h1c1.5|K3+sh1.25", REF),
         ("stack", "h2c1.5|K3+sh1.25", REF), ("CD4 x1.5 on top of H1 x2", "h2c1.5|K3+sh1.25", "h2c1|K3+sh1.25"),
         ("H1 x2 -> x3 (CD4 x1)", "h3c1|K3+sh1.25", "h2c1|K3+sh1.25"), ("H1 x3 -> x5 (CD4 x1)", "h5c1|K3+sh1.25", "h3c1|K3+sh1.25"),
         ("CD4 x1.5 -> x1 at H1 x2", "h2c1|K3+sh1.25", "h2c1.5|K3+sh1.25"),
         ("q 1.25 -> 1.5 on the stack", "h2c1.5|K3+sh1.5", "h2c1.5|K3+sh1.25"),
         ("key K3 -> LODO V1 (0.5, 3) on the stack", "h2c1.5|V1_5_3+sh1.25", "h2c1.5|K3+sh1.25")]
ESS = {l.strip() for l in (ROOT / "data/donors/k562_essential/perts.txt").read_text().splitlines() if l.strip()}


def top3(P, arm, ref, level=True):
    L = ":L" if level else ""
    out = set()
    for f in P.folds:
        Z = np.load(P.Z[f])
        tg = [str(t) for t in np.load(WORK / P.tag / f"fold-{f}-{P.fold_ops}.npz", allow_pickle=True)["targets"]]
        dx = (Z[f"{arm}{L}@R16__x"] - Z[f"{ref}{L}@R16__x"]).sum(0)
        out |= {tg[i] for i in np.argsort(-dx)[:3]}
    return out


def mean_price(P, arm, ref, level=True):
    p = [P.price(f, arm, 1.0, ref, level)["ov"] for f in P.folds]
    return float(np.mean(p)), int(sum(v > 0 for v in p))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", default="abc,real,efp,efr")
    ap.add_argument("--full", default="efr=efr_s1")
    ap.add_argument("--fold-energy", default="efr")
    a = ap.parse_args()
    fulls = dict(kv.split("=") for kv in a.full.split(",") if kv)
    fe = a.fold_energy.split(",")
    U = october_union()
    out, md = {}, ["| step | panel | all targets | minus each fold's top-3 (targets) | minus the 4 screens' targets (Arc-like) "
                   "| minus K562-essential |", "|---|---|---|---|---|---|"]
    for lab, arm, ref in STEPS:
        for t in a.panels.split(","):
            kw = dict(full_tag=fulls.get(t), energy="fold" if t in fe else "full")
            try:
                P = Panel(t, **kw)
            except FileNotFoundError:
                continue
            if len(P.folds) < 3 or not all(f"{arm}:L" in P.R[f] for f in P.folds):
                continue
            t3 = top3(P, arm, ref)
            r = dict(all=mean_price(P, arm, ref), minus_top3=mean_price(Panel(t, exclude=t3, **kw), arm, ref),
                     top3=sorted(t3), arc_like=mean_price(Panel(t, exclude=U, **kw), arm, ref),
                     minus_essential=mean_price(Panel(t, exclude=ESS, **kw), arm, ref))
            out.setdefault(lab, {})[t] = r
            fm = lambda v: f"{v[0]:+.5f} {v[1]}/3"
            md.append(f"| {lab} | {t} | {fm(r['all'])} | {fm(r['minus_top3'])} ({', '.join(sorted(t3))}) | {fm(r['arc_like'])} | "
                      f"{fm(r['minus_essential'])} |")
            print(md[-1], flush=True)
    (WORK / "robust438.json").write_text(json.dumps(out, indent=1))
    (WORK / "robust438.md").write_text("\n".join(md) + "\n")


if __name__ == "__main__":
    main()

"""640c: does 438's CD4 / H1 / key verdict on the October-like draws depend on the draws' strength mix?

Per fold, targets are put into abc's terciles of the held-out truth's response energy d (weak / mid / strong). Each 438
step is priced (levelled, k = 1, cancer-fold mean) on each tercile, and on the whole draw RE-WEIGHTED to abc's tercile
mix (1/3 each): weights w_i = (1/3) / frac(tercile of i) enter cos_agg = |sum w x| / sqrt(sum w a * sum w d) and the
pds means, energies fixed at the full panel (438's subset convention). abc itself is the known answer: its re-weighted
price must equal its plain price (weights all 1).
    uv run python experiments/640-arc-regime-slice/strength_slice640.py
"""
from __future__ import annotations

import sys
sys.dont_write_bytecode = True
from pathlib import Path
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/438-october-stack-repricing"))
from price438 import CANCER, WORK, Panel  # noqa: E402
from robust438 import STEPS  # noqa: E402

OUT = Path(__file__).resolve().parent / "work"


def terc_edges(f):
    d = np.asarray(np.load(WORK / "abc" / f"arrays438-{f}-held.npz")["d"], float)
    return np.quantile(np.log10(d), [1 / 3, 2 / 3])


EDGES = {f: terc_edges(f) for f in CANCER}


def mk(tag):
    kw = dict(full_tag="efr_s1" if tag == "efr" else None, energy="fold" if tag == "efr" else "full")
    return Panel(tag, **kw)


def weighted_stats(P, f, names, w):
    Z = P._arrays(f)
    d = Z["d"]
    out = {}
    for nm in names:
        t = f"{nm}@R16" if f"{nm}@R16__x" in Z.files else f"{nm}@D"
        x, a, pl, pp = Z[f"{t}__x"], Z[f"{t}__a"], Z[f"{t}__pdsl"], Z[f"{t}__pds"]
        cos = float(np.mean([abs((w * x[c]).sum()) / np.sqrt((w * a[c]).sum() * (w * d).sum()) for c in range(3)]))
        out[nm] = dict(cos=cos, pdsl=float((pl * w).sum(1).mean() / w.sum()), pds=float((pp * w).sum(1).mean() / w.sum()))
    return out


def price_w(P, arm, ref, wfun):
    names = {arm + ":L", ref + ":L", "h1c1|K3", "h1c1|none"}
    per = {}
    for f in P.folds:
        Z = P._arrays(f)
        ld = np.log10(np.asarray(Z["d"], float))
        w = wfun(f, ld)
        S = weighted_stats(P, f, names, w)
        per[f] = P.price(f, arm, 1.0, ref, True, S=S)["ov"]
    return float(np.mean(list(per.values()))), int(sum(v > 0 for v in per.values())), per


def terc(f, ld):
    e = EDGES[f]
    return np.where(ld <= e[0], 0, np.where(ld <= e[1], 1, 2))


def main():
    out, lines = {}, ["| step | panel | plain | weak tercile | mid | strong | re-weighted to abc's mix |", "|---|---|---|---|---|---|---|"]
    for lab, arm, ref in STEPS:
        for tag in ("abc", "real", "efp", "efr"):
            P = mk(tag)
            if len(P.folds) < 3 or not all(f"{arm}:L" in P.R[f] for f in P.folds):
                continue
            plain = price_w(P, arm, ref, lambda f, ld: np.ones_like(ld))
            ka = np.mean([P.price(f, arm, 1.0, ref, True)["ov"] for f in P.folds])
            assert abs(plain[0] - ka) < 1e-12, (lab, tag, plain[0], ka)          # KA: unit weights == 438's price
            cells = {}
            for k, nm in enumerate(("weak", "mid", "strong")):
                cells[nm] = price_w(P, arm, ref, lambda f, ld, k=k: (terc(f, ld) == k).astype(float))

            def rew(f, ld):
                t = terc(f, ld)
                fr = np.array([np.mean(t == k) for k in range(3)])
                return np.array([(1 / 3) / fr[k] for k in t])
            rw = price_w(P, arm, ref, rew)
            out.setdefault(lab, {})[tag] = dict(plain=plain[:2], **{k: v[:2] for k, v in cells.items()}, reweighted=rw[:2],
                                                reweighted_per_fold=rw[2])
            fm = lambda v: f"{v[0]:+.5f} {v[1]}/3"
            lines.append(f"| {lab} | {tag} | {fm(plain)} | {fm(cells['weak'])} | {fm(cells['mid'])} | {fm(cells['strong'])} | {fm(rw)} |")
            print(lines[-1], flush=True)
    (OUT / "strength_slice640.json").write_text(json.dumps(out, indent=1))
    (OUT / "strength_slice640.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

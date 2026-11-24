"""522 — October value of the key forms, in board units, on every surface (full panels and the Arc-like emulations).

overall ~= 0.0710 x k x r x ret + 0.1245 x dpds_lit    (analyse522's route: 518's identity, +1 % X = -0.0041 raw mse at
k 1 on A/B/C, X scaled by 438's per-fold retention ret; literal pds x 0.34 / 0.455 / 6), k = 1.31 (the key class) and 1.
r and dpds are vs the October default (cap 3), per cancer fold, then the panel mean. Sub-panel modes use the same folds'
retention (an assumption: an Arc-like re-draw keeps the fold's transfer level).

    uv run python experiments/522-october-bulk-key/price522.py
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from analyse522 import CANCER, RET  # noqa: E402

WORK = HERE / "work"
FORMS = ("cap6", "pctPROD", "pctPRODxS0.25", "pctPRODxS0.5", "S1", "cap8", "uncapped")
ALIAS = {"uncapped": "cap1e+09"}


def main():
    rows, out = [], {}
    P = rows.append
    # ---- full panels (arms522)
    surf = {}
    for p in ("abc", "efp", "efr", "real"):
        for f in CANCER:
            a = json.loads((WORK / f"arms522-{p}-{f}.json").read_text())["arms"]
            surf.setdefault(("full", p), {})[f] = {nm: (a[ALIAS.get(nm, nm)]["r"], a[ALIAS.get(nm, nm)]["dpds"]) for nm in FORMS}
    # ---- sub-panel modes (subpanel522*.json; newest record of a surface wins)
    S = {}
    for fp in sorted(WORK.glob("subpanel522*.json"), key=lambda q: q.stat().st_mtime):
        S.update(json.loads(fp.read_text()))
    for key, rec in S.items():
        p, f = key.split("/")
        for mode, m in rec["modes"].items():
            A = m["arms"]
            surf.setdefault((mode, p), {})[f] = {nm: (A[nm]["r_mean"], A[nm]["dpds_mean"]) for nm in FORMS if nm in A}
    P(f"{'surface':44s} " + " ".join(f"{nm:>22s}" for nm in FORMS))
    for (mode, p), folds in surf.items():
        cells = []
        for nm in FORMS:
            if not all(nm in folds[f] for f in folds):
                cells.append(f"{'':>22s}")
                continue
            r = np.mean([folds[f][nm][0] for f in folds])
            v131 = np.mean([0.0710 * 1.31 * folds[f][nm][0] * RET[p][f] + 0.1245 * folds[f][nm][1] for f in folds])
            v1 = np.mean([0.0710 * 1.0 * folds[f][nm][0] * RET[p][f] + 0.1245 * folds[f][nm][1] for f in folds])
            worst = min(folds[f][nm][0] for f in folds)
            out.setdefault(f"{mode}/{p}", {})[nm] = dict(r=r, overall_k131=v131, overall_k1=v1, worst_fold_r=worst, n_folds=len(folds))
            cells.append(f"{100 * r:+6.2f}% {v131:+.4f}/{v1:+.4f}")
        P(f"{mode + ' / ' + p + ' (' + str(len(folds)) + ')':44s} " + " ".join(cells))
    P("\ncell = mean r vs cap 3 (relative cos) and overall at k 1.31 / k 1 (board units, retention-scaled)")
    # ---- worst cases across surfaces
    P("\nworst case across surfaces (panel mean r; worst single fold r):")
    for nm in FORMS:
        vals = {k: v[nm] for k, v in out.items() if nm in v}
        wk = min(vals, key=lambda k: vals[k]["r"])
        wf = min(vals, key=lambda k: vals[k]["worst_fold_r"])
        P(f"  {nm:16s} worst panel {wk:40s} {100 * vals[wk]['r']:+6.2f} %   worst fold ({wf}) {100 * vals[wf]['worst_fold_r']:+6.2f} %")
    txt = "\n".join(rows)
    print(txt)
    (WORK / "price522.txt").write_text(txt + "\n")
    (WORK / "price522.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

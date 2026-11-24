"""438 diagnostic: how well does each residual (CD4 deflated, H1 union) and the whitened donor base point at the PLOLO
truth, per carrying row, on each panel? Row cosine over the fold's scored genes, residual restricted to its support.

    uv run python diag438.py <panel> [<panel> ...]        (one panel per process: bind() is global)
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402

import numpy as np  # noqa: E402

import lib438 as L  # noqa: E402


def rowcos(R, T):
    n = np.linalg.norm(R, axis=1) * np.linalg.norm(T, axis=1)
    ok = n > 0
    return (np.einsum("ij,ij->i", R, T)[ok] / n[ok]), ok


def main(panel, folds):
    ctx = L.bind(panel)
    out = {}
    for held in folds:
        fold = L.build_fold438(ctx, held, "held")
        sc, targets = fold["scored"], list(fold["targets"])
        T = fold["truth"][:, sc]
        cd4, partial, h1, _ = L.residuals(ctx, targets)
        dn = L.Donors438(targets, list(fold["sources"]))
        pre = L.pre_normrel(dn, list(fold["sources"]), "cap", -1.5, fallback_rows=fold["pre"])
        whit = fold["op"].apply(pre * L.ALPHA, fold["wcols"], L.GAMMA, L.LAM)
        r = {}
        for tag, M in (("cd4", np.where(np.asarray(cd4.gene_votes) > 0, np.asarray(cd4.mean_effect), 0.0)),
                       ("h1", np.where(np.asarray(h1.gene_votes) > 0, np.asarray(h1.mean_effect), 0.0)),
                       ("base", whit)):
            c, ok = rowcos(M[:, sc], T)
            r[tag] = dict(n=int(ok.sum()), mean=float(c.mean()), median=float(np.median(c)), frac_pos=float((c > 0).mean()))
        out[held] = r
        print(f"{panel} {held}: " + " | ".join(f"{k} n {v['n']} mean cos {v['mean']:+.4f} median {v['median']:+.4f} "
                                             f"pos {v['frac_pos']:.2f}" for k, v in r.items()), flush=True)
    f = ctx.work / "diag438.json"
    prev = json.loads(f.read_text()) if f.exists() else {}
    prev.update(out)
    f.write_text(json.dumps(prev, indent=1))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:] or list(L.CANCER))     # one fold per process keeps two operators out of RAM

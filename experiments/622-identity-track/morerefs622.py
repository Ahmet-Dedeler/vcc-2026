#!/usr/bin/env python3
"""622d — do MORE Flex reference profiles strengthen the identity axis's 3-reference board form?

Board form: each A/B/C context vs the other two (+ H1) -> only 2-3 references. Extra Flex profiles on disk: Satija's
VIPerturb pilot NTC pseudobulks (597: K562, HAP1, HEK293; Flex v1, another lab) and viperturb_k562's controls
(Satija). Read on 636's gates (gate sign acc of z(C4) + lam z(I), paired vs CLEAN4, target bootstrap) with references
{A,B,C}[+h1] vs {A,B,C}[+h1] + Satija's three. Analysis only, < 2 GB.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402

RNG = np.random.default_rng(6223)


def main() -> None:
    I = Inst()
    prof = {}
    for s, f in SURF.items():
        z = np.load(G480 / f, allow_pickle=True)
        prof[s] = dict(zip(z["genes"].astype(str), z["ctrl_cpm"].astype(float)))
    z = np.load(HERE / "work/abc_cpm622.npz", allow_pickle=True)
    for k in ("ctxA", "ctxB", "ctxC"):
        prof[k] = dict(zip(z["genes"].astype(str), z[k].astype(float)))
    sat = np.load(ROOT / "experiments/597-flex-hek293-floor/work/pseudobulk597.npz", allow_pickle=True)
    for i, ln in enumerate(sat["lines"].astype(str)):
        c = sat["counts"][i].astype(float)
        prof[f"sat_{ln}"] = dict(zip(sat["genes"].astype(str), 1e6 * c / c.sum()))  # pooled CPM (pseudobulk proxy)
    out = {}
    for s in I.S:
        D = I.S[s]
        C4, _ = clean4_lolo(I, s)
        zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
        gate, dn = D["gate"], D["lfc"] < 0
        base = [k for k in ("ctxA", "ctxB", "ctxC", "h1") if k != s]
        sets = {"flex": base, "flex+satija": base + ["sat_K562", "sat_HAP1", "sat_HEK293"]}
        res = {}
        for name, refs in sets.items():
            common = sorted(set(prof[s]).intersection(*(set(prof[r]) for r in refs)))
            Ls = np.log1p(np.array([prof[s][g] for g in common]))
            Lr = np.mean([np.log1p(np.array([prof[r][g] for g in common])) for r in refs], 0)
            zi = I.gene_vec(s, common, Ls - Lr)
            zi = (zi - np.nanmean(zi)) / np.nanstd(zi)
            for lam in (0.5, 1.0):
                sc = zc + lam * zi
                both = np.isfinite(zc) & np.isfinite(sc)
                m = gate & both[None, :]
                okS = (((sc > 0)[None, :] == dn) & m).sum(1).astype(float)
                okC = (((zc > 0)[None, :] == dn) & m).sum(1).astype(float)
                n = m.sum(1).astype(float)
                idx = np.flatnonzero(n > 0)
                bs = [(okS[b].sum() - okC[b].sum()) / n[b].sum() for b in (RNG.choice(idx, idx.size) for _ in range(2000))]
                res[f"{name}|lam{lam}"] = dict(genes=len(common), d=float((okS.sum() - okC.sum()) / n.sum()),
                                               ci=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))])
                r = res[f"{name}|lam{lam}"]
                print(f"{s:8s} {name:12s} λ{lam}: genes {len(common):5d}  Δ {r['d']:+.4f} [{r['ci'][0]:+.4f}, {r['ci'][1]:+.4f}]", flush=True)
        out[s] = res
    (HERE / "work/morerefs622.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()

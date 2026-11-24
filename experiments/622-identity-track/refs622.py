#!/usr/bin/env python3
"""622b — which REFERENCE SET should the identity axis use? (d4cd5f's (b)/(c))

For each 636 surface s, I_s = log1p(ctrl CPM_s) - mean over references r of log1p(ctrl CPM_r), on genes all chosen
profiles measure. Reference sets:
  mixed      every other profile (639's construction: the other 636 surfaces + A/B/C)
  flex       Flex profiles only: A, B, C (+ h1 when s != h1)     <- D/E/F's day-of set is Flex-heavy
  threep     public 3' profiles only (the other 3' surfaces)
  mixed-minus-r   leave-one-reference-out of the mixed set
Read: gate sign accuracy of z(P_C4 - .5) + lam z(I) on 636's gate, paired vs CLEAN4 (LOLO) on genes both cover, 2,000-draw
target bootstrap, lam 0.5 and 1. A/B/C control CPM cached in work/abc_cpm622.npz. Analysis only, < 3 GB.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402

RNG = np.random.default_rng(6221)
LAMS = (0.5, 1.0)
CACHE = WORK / "abc_cpm622.npz"


def abc_profiles() -> dict[str, dict[str, float]]:
    if not CACHE.exists():
        import importlib.util
        spec = importlib.util.spec_from_file_location("s622", HERE / "stab622.py")
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)  # type: ignore[union-attr]
        genes, abc = m.abc_cpm(halves=False)
        np.savez_compressed(CACHE, genes=np.array(genes), **{k: v for k, v in abc.items()})
    z = np.load(CACHE, allow_pickle=True)
    genes = z["genes"].astype(str)
    return {k: dict(zip(genes, z[k].astype(float))) for k in ("ctxA", "ctxB", "ctxC")}


def main() -> None:
    I = Inst()
    prof = {}
    for s, f in SURF.items():
        z = np.load(G480 / f, allow_pickle=True)
        prof[s] = dict(zip(z["genes"].astype(str), z["ctrl_cpm"].astype(float)))
    prof.update(abc_profiles())
    flex_all = ["ctxA", "ctxB", "ctxC", "h1"]
    threep_all = ["hct116", "hek293t", "k562"]
    out = {}
    for s in I.S:
        D = I.S[s]
        C4, _ = clean4_lolo(I, s)
        zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
        gate, dn = D["gate"], D["lfc"] < 0
        others = [k for k in prof if k != s]
        sets = {"mixed": others,
                "flex": [k for k in flex_all if k != s],
                "threep": [k for k in threep_all if k != s]}
        for r in others:
            sets[f"mixed-minus-{r}"] = [k for k in others if k != r]

        def paired(score):
            both = np.isfinite(zc) & np.isfinite(score)
            m = gate & both[None, :]
            okS = (((score > 0)[None, :] == dn) & m).sum(1).astype(float)
            okC = (((zc > 0)[None, :] == dn) & m).sum(1).astype(float)
            n = m.sum(1).astype(float)
            idx = np.flatnonzero(n > 0)
            bs = [(okS[b].sum() - okC[b].sum()) / n[b].sum() for b in (RNG.choice(idx, idx.size) for _ in range(2000))]
            return dict(d=float((okS.sum() - okC.sum()) / n.sum()), ci=[float(np.quantile(bs, .025)), float(np.quantile(bs, .975))],
                        n=int(n.sum()))

        res = {}
        for name, refs in sets.items():
            common = sorted(set(prof[s]).intersection(*(set(prof[r]) for r in refs)))
            Ls = np.log1p(np.array([prof[s][g] for g in common]))
            Lr = np.mean([np.log1p(np.array([prof[r][g] for g in common])) for r in refs], 0)
            zi = I.gene_vec(s, common, Ls - Lr)
            zi = (zi - np.nanmean(zi)) / np.nanstd(zi)
            res[name] = {"refs": refs, "genes": len(common), **{f"lam{l}": paired(zc + l * zi) for l in LAMS}}
            r1 = res[name]
            print(f"{s:8s} {name:24s} genes {len(common):5d} | " + "  ".join(
                f"λ{l}: {r1[f'lam{l}']['d']:+.4f} [{r1[f'lam{l}']['ci'][0]:+.4f},{r1[f'lam{l}']['ci'][1]:+.4f}]" for l in LAMS), flush=True)
        out[s] = res
    (WORK / "refs622.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
